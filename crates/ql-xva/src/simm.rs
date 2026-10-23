use ql_ad::Real;

/// The six ISDA SIMM risk classes (in the ordering used by the v2.6
/// cross-risk-class correlation matrix).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
#[repr(usize)]
pub enum SimmRiskClass {
    InterestRate = 0,
    CreditQualifying = 1,
    CreditNonQualifying = 2,
    Equity = 3,
    Commodity = 4,
    Fx = 5,
}

/// Aggregate class margins using the ISDA SIMM v2.6 cross-risk-class matrix.
/// The input ordering is `[IR, CreditQ, CreditNQ, Equity, Commodity, FX]`.
pub fn compute_isda_simm_6_risk_classes<T: Real>(class_margins: &[T; 6]) -> T {
    // ISDA SIMM v2.6, Table 2.2.  The matrix is symmetric with unit diagonal.
    const PSI: [[f64; 6]; 6] = [
        [1.00, 0.28, 0.21, 0.18, 0.26, 0.18],
        [0.28, 1.00, 0.42, 0.20, 0.26, 0.22],
        [0.21, 0.42, 1.00, 0.20, 0.24, 0.24],
        [0.18, 0.20, 0.20, 1.00, 0.27, 0.32],
        [0.26, 0.26, 0.24, 0.27, 1.00, 0.28],
        [0.18, 0.22, 0.24, 0.32, 0.28, 1.00],
    ];
    let mut total = T::zero();
    for i in 0..6 {
        for j in 0..6 {
            total += T::from_f64(PSI[i][j]) * class_margins[i] * class_margins[j];
        }
    }
    // Numerical roundoff can make a nearly-zero portfolio marginally negative
    // for scalar f64.  Keep AD values in-domain without stripping tangents.
    total.max_val(T::zero()).sqrt()
}
#[derive(Clone, Debug)]
pub struct SimmBucket<T: Real> {
    pub sensitivities: Vec<T>,
    pub risk_weights: Vec<f64>,
    pub concentration_threshold: f64,
    pub intra_corr: Vec<Vec<f64>>,
}
pub fn compute_isda_simm<T: Real>(b: &[SimmBucket<T>], corr: &[Vec<f64>]) -> T {
    assert_eq!(
        corr.len(),
        b.len(),
        "inter-bucket correlation dimension mismatch"
    );
    assert!(corr.iter().all(|row| row.len() == b.len()));
    for x in b {
        assert_eq!(x.sensitivities.len(), x.risk_weights.len());
        assert_eq!(x.intra_corr.len(), x.sensitivities.len());
        assert!(
            x.intra_corr
                .iter()
                .all(|row| row.len() == x.sensitivities.len())
        );
        assert!(x.concentration_threshold.is_finite() && x.concentration_threshold > 0.0);
    }
    let mut ks = Vec::new();
    let mut ss = Vec::new();
    for x in b {
        let concentration_ratios: Vec<T> = x
            .sensitivities
            .iter()
            .map(|s| {
                let threshold = x.concentration_threshold.max(1e-12);
                let abs_ratio = s.abs() / threshold;
                // Keep the concentration multiplier in the AD domain.  Using
                // only `s.val()` here would silently drop the derivative of
                // the ISDA concentration adjustment for large sensitivities.
                if abs_ratio.val() > 1.0 {
                    abs_ratio.sqrt()
                } else if abs_ratio.val() == 1.0 {
                    // At the kink use the symmetric finite-difference
                    // subgradient (the mean of the two one-sided slopes).
                    T::one() + (abs_ratio - T::one()) * 0.25
                } else {
                    T::one()
                }
            })
            .collect();
        let ws: Vec<T> = x
            .sensitivities
            .iter()
            .zip(&x.risk_weights)
            .zip(&concentration_ratios)
            .map(|((s, w), concentration)| *s * *concentration * *w)
            .collect();
        let mut k2 = T::zero();
        for i in 0..ws.len() {
            for j in 0..ws.len() {
                // SIMM v2.6 scales off-diagonal correlations by the ratio of
                // concentration factors.  Compute this in T so AD sensitivities
                // remain exact through the concentration adjustment.
                let concentration_corr = if i == j {
                    T::one()
                } else {
                    let lo = concentration_ratios[i].min_val(concentration_ratios[j]);
                    let hi = concentration_ratios[i].max_val(concentration_ratios[j]);
                    if hi.val().abs() > 1e-14 {
                        lo / hi
                    } else {
                        T::one()
                    }
                };
                k2 += T::from_f64(
                    x.intra_corr
                        .get(i)
                        .and_then(|r| r.get(j))
                        .copied()
                        .unwrap_or(if i == j { 1. } else { 0. }),
                ) * concentration_corr
                    * ws[i]
                    * ws[j];
            }
        }
        let k = k2.sqrt();
        ks.push(k);
        ss.push(
            ws.iter()
                .copied()
                .fold(T::zero(), |a, z| a + z)
                .min_val(k)
                .max_val(-k),
        );
    }
    let mut total = T::zero();
    for i in 0..ks.len() {
        total += ks[i] * ks[i];
        for j in 0..ks.len() {
            if i != j {
                total += T::from_f64(corr.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.))
                    * ss[i]
                    * ss[j];
            }
        }
    }
    total.sqrt()
}
