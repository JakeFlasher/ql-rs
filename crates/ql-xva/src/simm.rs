use ql_ad::Real;
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
        let ws: Vec<T> = x
            .sensitivities
            .iter()
            .zip(&x.risk_weights)
            .map(|(s, w)| {
                let threshold = x.concentration_threshold.max(1e-12);
                let abs_ratio = s.abs() / threshold;
                // Keep the concentration multiplier in the AD domain.  Using
                // only `s.val()` here would silently drop the derivative of
                // the ISDA concentration adjustment for large sensitivities.
                let concentration = if abs_ratio.val() > 1.0 {
                    abs_ratio.sqrt()
                } else if abs_ratio.val() == 1.0 {
                    // At the kink use the symmetric finite-difference
                    // subgradient (the mean of the two one-sided slopes).
                    T::one() + (abs_ratio - T::one()) * 0.25
                } else {
                    T::one()
                };
                *s * concentration * *w
            })
            .collect();
        let mut k2 = T::zero();
        for i in 0..ws.len() {
            for j in 0..ws.len() {
                k2 += T::from_f64(
                    x.intra_corr
                        .get(i)
                        .and_then(|r| r.get(j))
                        .copied()
                        .unwrap_or(if i == j { 1. } else { 0. }),
                ) * ws[i]
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
