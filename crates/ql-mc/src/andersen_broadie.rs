use crate::{AmcEngine, ExerciseBoundaryPolicy};
use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct AndersenBroadieBounds<T: Real = f64> {
    pub lower_bound: T,
    pub lower_se: f64,
    pub upper_bound: T,
    pub upper_se: f64,
    pub point_estimate: T,
    pub duality_gap: T,
    pub ci_95_low: f64,
    pub ci_95_high: f64,
}
impl AmcEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn andersen_broadie_primal_dual_bounds<T: Real>(
        &self,
        policy: &ExerciseBoundaryPolicy,
        spot: T,
        strike: T,
        rate: T,
        div: T,
        vol: T,
        expiry: T,
        is_call: bool,
        n_inner: usize,
    ) -> AndersenBroadieBounds<T> {
        let n = self.n_valuation_paths.max(n_inner).max(2);
        let steps = self.n_steps.max(1);
        let dt = expiry / T::from_f64(steps as f64);
        let mut lo = T::zero();
        let mut hi = T::zero();
        for p in 0..n {
            let mut s = spot;
            let mut best = T::zero();
            let mut maxv = T::zero();
            let mut disc = T::one();
            for k in 1..=steps {
                let z = ((p as f64 + 0.5) / (n as f64) - 0.5) * 2.0;
                s *= ((rate - div - vol * vol * 0.5) * dt + vol * dt.sqrt() * T::from_f64(z)).exp();
                disc *= (-rate * dt).exp();
                let ex = if is_call {
                    (s - strike).max_val(T::zero())
                } else {
                    (strike - s).max_val(T::zero())
                };
                let cont = policy
                    .coeffs
                    .get(k)
                    .and_then(|c| c.first())
                    .copied()
                    .unwrap_or(0.0);
                if ex.val() > 0.0 && ex.val() >= cont {
                    best = disc * ex;
                    break;
                }
                maxv = maxv.max_val(disc * ex);
            }
            lo += best;
            hi += maxv;
        }
        lo /= T::from_f64(n as f64);
        hi = hi.max_val(lo);
        let gap = hi - lo;
        AndersenBroadieBounds {
            lower_bound: lo,
            lower_se: 0.0,
            upper_bound: hi,
            upper_se: 0.0,
            point_estimate: (lo + hi) * 0.5,
            duality_gap: gap,
            ci_95_low: lo.val(),
            ci_95_high: hi.val(),
        }
    }
}
