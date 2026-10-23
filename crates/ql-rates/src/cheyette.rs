use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_mc::PhiloxRng;
#[derive(Copy, Clone, Debug)]
pub struct Cheyette1F<T: Real> {
    pub mu: T,
}
impl<T: Real> Cheyette1F<T> {
    pub fn g_fn(&self, t: f64, big_t: f64) -> T {
        let tau = (big_t - t).max(0.);
        Self::g(self.mu, tau)
    }
    fn g(mu: T, tau: f64) -> T {
        if (mu.val() * tau).abs() < 1e-6 {
            T::from_f64(tau) - mu * tau * tau * 0.5 + mu * mu * tau * tau * tau / 6.
        } else {
            (T::one() - (-mu * tau).exp()) / mu
        }
    }
    pub fn discount_bond<C: YieldTermStructure<T>>(
        &self,
        curve: &C,
        t: f64,
        big_t: f64,
        x: T,
        y: T,
    ) -> T {
        let g = self.g_fn(t, big_t);
        curve.discount_factor(big_t) / curve.discount_factor(t)
            * (-g * x - T::from_f64(0.5) * g * g * y).exp()
    }
    pub fn explicit_local_variance_first_order(
        mu: T,
        k: T,
        w: T,
        dw_dt: T,
        dw_dk: T,
        d2w_dk2: T,
    ) -> T {
        let den = (T::one() - k * dw_dk / (T::from_f64(2.) * w)).powf(2.)
            + T::from_f64(0.5) * (d2w_dk2 - dw_dk * dw_dk / (T::from_f64(2.) * w));
        (dw_dt + mu * (T::from_f64(2.) * w - k * dw_dk) + w * dw_dk) / den
    }
    pub fn explicit_local_variance_third_order(
        mu: T,
        k: T,
        w: T,
        dw_dt: T,
        dw_dk: T,
        d2w_dk2: T,
    ) -> T {
        Self::explicit_local_variance_first_order(mu, k, w, dw_dt, dw_dk, d2w_dk2)
            + dw_dk * dw_dk * dw_dk
    }
    pub fn simulate_short_rate_call(
        &self,
        expiry: f64,
        strike: f64,
        n_paths: usize,
        n_steps: usize,
        seed: u64,
        total_var_fn: impl Fn(f64, f64) -> (f64, f64, f64, f64),
    ) -> f64 {
        let rng = PhiloxRng::new(seed);
        let dt = expiry / n_steps.max(1) as f64;
        let mut sum = 0.;
        for p in 0..n_paths.max(1) {
            let sign = if p % 2 == 0 { 1. } else { -1. };
            let id = (p / 2) as u64;
            let mut x = 0.;
            let mut y = 0.;
            for j in 0..n_steps.max(1) {
                let t = (j + 1) as f64 * dt;
                let (wv, dwv, dkv, d2v) = total_var_fn(t, x);
                let w = T::from_f64(wv.max(1e-10));
                let var = Self::explicit_local_variance_third_order(
                    self.mu,
                    T::from_f64(x),
                    w,
                    T::from_f64(dwv),
                    T::from_f64(dkv),
                    T::from_f64(d2v),
                )
                .val()
                .max(1e-10);
                let z = rng.normal4(id, j as u32)[0] * sign;
                x += (y - self.mu.val() * x) * dt + var.sqrt() * dt.sqrt() * z;
                y += (var - 2. * self.mu.val() * y) * dt;
            }
            sum += (x - strike).max(0.);
        }
        sum / n_paths.max(1) as f64
    }
}
