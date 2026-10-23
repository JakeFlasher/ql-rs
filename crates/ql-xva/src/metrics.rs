use crate::NpvCube;
use ql_ad::Real;
use ql_curves::YieldTermStructure;
#[derive(Copy, Clone, Debug)]
pub struct CsaSpec {
    pub threshold_cpty: f64,
    pub threshold_self: f64,
    pub mta: f64,
}
impl CsaSpec {
    /// Remaining signed uncollateralized exposure for a netting set.
    ///
    /// Variation margin is called only when the positive/negative exposure
    /// exceeds the threshold plus the minimum transfer amount.  Once a call
    /// is made the residual is the threshold (not the collateral call itself).
    pub fn collateralized_exposure<T: Real>(&self, x: T) -> T {
        if x.val() >= 0.0 {
            collateralized_exposure(x, self.threshold_cpty, self.mta)
        } else {
            -collateralized_exposure(-x, self.threshold_self, self.mta)
        }
    }
}

/// Positive remaining exposure after a one-way CSA margin call.
///
/// The branch is selected using the primal value while all arithmetic in the
/// selected branch remains in `T`, preserving Dual/Dual2/TapeVar tangents.
/// For a zero-threshold, zero-MTA CSA this is identically zero for positive
/// mark-to-market values, as expected for continuous variation margin.
pub fn collateralized_exposure<T: Real>(mtm: T, threshold_cpty: f64, mta: f64) -> T {
    let e = mtm.max_val(T::zero());
    let threshold = T::from_f64(threshold_cpty.max(0.0));
    let excess = e - threshold;
    if excess.val() > mta.max(0.0) {
        e.min_val(threshold)
    } else {
        e
    }
}
#[derive(Clone, Debug)]
pub struct XvaReport<T: Real> {
    pub ee: Vec<T>,
    pub ene: Vec<T>,
    pub pfe_95: Vec<f64>,
    pub epe: T,
    pub cva: T,
    pub dva: T,
    pub fva: T,
}

/// Expected positive exposure by simulation time, aggregated across all
/// netting sets.  The result is intentionally `f64` for charting and report
/// consumers; the AD-preserving XVA valuation remains available through
/// `compute_xva`.
pub fn expected_positive_exposure<T: Real>(cube: &NpvCube<T>) -> Vec<f64> {
    (0..cube.times.len())
        .map(|k| {
            let mut sum = 0.0;
            for path in 0..cube.n_paths {
                let mut value = 0.0;
                for ns in 0..cube.n_netting_sets {
                    // Legal netting sets may not offset one another.  Floor
                    // each set before summing to the counterparty profile.
                    value += cube.get(ns, k, path).val().max(0.0);
                }
                sum += value;
            }
            sum / cube.n_paths.max(1) as f64
        })
        .collect()
}

/// Empirical potential future exposure by simulation time.
pub fn potential_future_exposure<T: Real>(cube: &NpvCube<T>, quantile: f64) -> Vec<f64> {
    let q = quantile.clamp(0.0, 1.0);
    (0..cube.times.len())
        .map(|k| {
            let mut values = Vec::with_capacity(cube.n_paths);
            for path in 0..cube.n_paths {
                let mut value = 0.0;
                for ns in 0..cube.n_netting_sets {
                    value += cube.get(ns, k, path).val().max(0.0);
                }
                values.push(value);
            }
            values.sort_by(f64::total_cmp);
            let index = (q * values.len().saturating_sub(1) as f64).round() as usize;
            values
                .get(index.min(values.len().saturating_sub(1)))
                .copied()
                .unwrap_or(0.0)
        })
        .collect()
}
#[allow(clippy::too_many_arguments)]
pub fn compute_xva<T: Real>(
    cube: &NpvCube<T>,
    curve: &impl YieldTermStructure<T>,
    csa: Option<&CsaSpec>,
    hc: T,
    lc: T,
    hs: T,
    ls: T,
    fs: T,
) -> XvaReport<T> {
    let mut ee = vec![];
    let mut ene = vec![];
    let mut pfe = vec![];
    let mut cva = T::zero();
    let mut dva = T::zero();
    let mut fva = T::zero();
    for k in 0..cube.times.len() {
        // Keep positive and negative exposures separate all the way through
        // aggregation.  Summing signed MTMs first would illegally net two
        // distinct legal netting sets (e.g. +10 and -8 -> +2 instead of +10).
        let mut pos_vals = Vec::with_capacity(cube.n_paths);
        let mut neg_vals = Vec::with_capacity(cube.n_paths);
        for p in 0..cube.n_paths {
            let mut pos = T::zero();
            let mut neg = T::zero();
            for ns in 0..cube.n_netting_sets {
                let v = cube.get(ns, k, p);
                match csa {
                    Some(c) => {
                        pos += collateralized_exposure(v, c.threshold_cpty, c.mta);
                        neg -= collateralized_exposure(-v, c.threshold_self, c.mta);
                    }
                    None => {
                        pos += v.max_val(T::zero());
                        neg += v.min_val(T::zero());
                    }
                }
            }
            pos_vals.push(pos);
            neg_vals.push(neg);
        }
        let e = pos_vals.iter().copied().fold(T::zero(), |a, b| a + b)
            / T::from_f64(cube.n_paths.max(1) as f64);
        let n = neg_vals.iter().copied().fold(T::zero(), |a, b| a + b)
            / T::from_f64(cube.n_paths.max(1) as f64);
        let mut raw = pos_vals.iter().map(|x| x.val()).collect::<Vec<_>>();
        raw.sort_by(|a, b| a.total_cmp(b));
        let idx = ((0.95 * (raw.len().saturating_sub(1) as f64)).round() as usize)
            .min(raw.len().saturating_sub(1));
        ee.push(e);
        ene.push(n);
        pfe.push(raw.get(idx).copied().unwrap_or(0.));
        let prev = if k == 0 { 0. } else { cube.times[k - 1] };
        let t = cube.times[k];
        let df = curve.discount_factor(t);
        let sc0 = (-hc * prev).exp();
        let sc = (-hc * t).exp();
        let sb0 = (-hs * prev).exp();
        let sb = (-hs * t).exp();
        cva += lc * e * df * (sc0 - sc) * sb0;
        dva -= ls * n * df * (sb0 - sb) * sc0;
        fva += fs * (e + n) * df * sc * sb * T::from_f64(t - prev);
    }
    // EPE is a time integral, not an unweighted arithmetic mean.  Include
    // EE(0)=0 and use trapezoidal integration on the supplied irregular grid.
    let epe = if let Some(&t_last) = cube.times.last() {
        if t_last > 0.0 && ee.len() == cube.times.len() {
            let mut integral = T::zero();
            let mut prev_t = 0.0;
            let mut prev_ee = T::zero();
            for (&t, &e_t) in cube.times.iter().zip(&ee) {
                let dt = (t - prev_t).max(0.0);
                integral += (prev_ee + e_t) * (0.5 * dt);
                prev_t = t;
                prev_ee = e_t;
            }
            integral / T::from_f64(t_last)
        } else {
            T::zero()
        }
    } else {
        T::zero()
    };
    XvaReport {
        ee,
        ene,
        pfe_95: pfe,
        epe,
        cva,
        dva,
        fva,
    }
}
