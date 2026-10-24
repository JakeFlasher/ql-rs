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
        // The routine is intentionally pathwise and AD-preserving, but still
        // needs a well-defined state space.  A zero strike would make the
        // calibrated Chebyshev basis singular; in that case fall back to the
        // raw (constant) policy value rather than producing NaNs.
        let n = self.n_valuation_paths.max(n_inner).max(2);
        let steps = self.n_steps.max(1);
        let dt = expiry / T::from_f64(steps as f64);
        let mut lo = T::zero();
        let mut hi = T::zero();
        let mut lo2 = T::zero();
        let mut hi2 = T::zero();
        let rng = crate::PhiloxRng::new(self.seed.wrapping_add(2_000_000));
        for p in 0..n {
            let mut s = spot;
            let mut best = T::zero();
            // The pathwise supremum of discounted intrinsic values is a
            // valid (albeit deliberately conservative) dual upper bound:
            // every stopping rule is bounded by the supremum on each path.
            // The previous implementation only updated this before an early
            // exercise and could therefore report an identically zero upper
            // bound for an in-the-money policy.
            let mut maxv = T::zero();
            let mut disc = T::one();
            for k in 1..=steps {
                // Use deterministic counter-based normals so the result is
                // reproducible and antithetic pairs share exactly the same
                // random stream as the rest of the MC crate.
                let pid = (p / 2) as u64;
                let sign = if p % 2 == 0 { 1.0 } else { -1.0 };
                let z = rng.normal4(pid, (k - 1) as u32)[0] * sign;
                s *= ((rate - div - vol * vol * 0.5) * dt + vol * dt.sqrt() * T::from_f64(z)).exp();
                disc *= (-rate * dt).exp();
                let ex = if is_call {
                    (s - strike).max_val(T::zero())
                } else {
                    (strike - s).max_val(T::zero())
                };
                let discounted = disc * ex;
                maxv = maxv.max_val(discounted);
                // Exercise decisions are frozen to the supplied Longstaff--
                // Schwartz policy.  Evaluate the full Chebyshev basis (not
                // just the intercept) on primal values while retaining the
                // AD value of the selected payoff.
                let cont = policy.coeffs.get(k).map_or(f64::NEG_INFINITY, |c| {
                    if c.is_empty() {
                        f64::NEG_INFINITY
                    } else if strike.val().abs() <= 1e-14 {
                        c[0]
                    } else {
                        let x = s.val() / strike.val() - 1.0;
                        let mut b0 = 1.0;
                        let mut b1 = x;
                        let mut value = c[0] * b0;
                        if c.len() > 1 {
                            value += c[1] * b1;
                        }
                        for coeff in c.iter().skip(2) {
                            let b2 = 2.0 * x * b1 - b0;
                            value += coeff * b2;
                            b0 = b1;
                            b1 = b2;
                        }
                        value
                    }
                });
                if ex.val() > 0.0 && ex.val() >= cont {
                    best = discounted;
                    break;
                }
            }
            lo += best;
            hi += maxv;
            lo2 += best * best;
            hi2 += maxv * maxv;
        }
        let den = T::from_f64(n as f64);
        lo /= den;
        hi /= den;
        hi = hi.max_val(lo);
        let gap = hi - lo;
        // Report Monte-Carlo standard errors for the two pathwise estimators.
        // Keep the primal-only values in the public f64 fields while the
        // bounds themselves retain all AD tangents.
        let n_f = n as f64;
        let lo_mean = lo.val();
        let hi_mean = hi.val();
        let lo_var = (lo2.val() / n_f - lo_mean * lo_mean).max(0.0);
        let hi_var = (hi2.val() / n_f - hi_mean * hi_mean).max(0.0);
        let lo_se = (lo_var / n_f).sqrt();
        let hi_se = (hi_var / n_f).sqrt();
        AndersenBroadieBounds {
            lower_bound: lo,
            lower_se: lo_se,
            upper_bound: hi,
            upper_se: hi_se,
            point_estimate: (lo + hi) * 0.5,
            duality_gap: gap,
            ci_95_low: lo_mean - 1.96 * lo_se,
            ci_95_high: hi_mean + 1.96 * hi_se,
        }
    }
}
