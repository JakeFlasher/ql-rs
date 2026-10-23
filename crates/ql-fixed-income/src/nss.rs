use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct NelsonSiegelSvensson<T: Real> {
    pub beta0: T,
    pub beta1: T,
    pub beta2: T,
    pub beta3: T,
    pub tau1: T,
    pub tau2: T,
}
impl<T: Real> NelsonSiegelSvensson<T> {
    fn h1(x: T) -> T {
        if x.val().abs() < 1e-6 {
            T::one() - x * 0.5 + x * x / 6. - x * x * x / 24.
        } else {
            (T::one() - (-x).exp()) / x
        }
    }
    pub fn zero_rate(&self, t: T) -> T {
        let x1 = t / self.tau1;
        let x2 = t / self.tau2;
        let h_1 = Self::h1(x1);
        let h_2 = h_1 - (-x1).exp();
        self.beta0 + self.beta1 * h_1 + self.beta2 * h_2 + self.beta3 * (Self::h1(x2) - (-x2).exp())
    }
    pub fn forward_rate(&self, t: T) -> T {
        let e1 = (-t / self.tau1).exp();
        let e2 = (-t / self.tau2).exp();
        self.beta0
            + self.beta1 * e1
            + self.beta2 * (t / self.tau1) * e1
            + self.beta3 * (t / self.tau2) * e2
    }
    pub fn discount_factor(&self, t: T) -> T {
        (-self.zero_rate(t) * t).exp()
    }
    pub fn fit_ols_for_taus(times: &[f64], rates: &[f64], tau1: f64, tau2: f64) -> Self {
        use nalgebra::{DMatrix, DVector};
        let mut a = DMatrix::zeros(times.len(), 4);
        for (i, &t) in times.iter().enumerate() {
            let x1 = t / tau1;
            let x2 = t / tau2;
            // Keep the OLS design matrix finite at a zero/very short
            // maturity, matching the AD-safe `h1` used by `zero_rate`.
            let h1 = if x1.abs() < 1e-6 {
                1. - x1 * 0.5 + x1 * x1 / 6. - x1 * x1 * x1 / 24.
            } else {
                (1. - (-x1).exp()) / x1
            };
            let h2 = h1 - (-x1).exp();
            let h1_2 = if x2.abs() < 1e-6 {
                1. - x2 * 0.5 + x2 * x2 / 6. - x2 * x2 * x2 / 24.
            } else {
                (1. - (-x2).exp()) / x2
            };
            let h3 = h1_2 - (-x2).exp();
            a[(i, 0)] = 1.;
            a[(i, 1)] = h1;
            a[(i, 2)] = h2;
            a[(i, 3)] = h3;
        }
        let y = DVector::from_column_slice(rates);
        let b = a
            .svd(true, true)
            .solve(&y, 1e-12)
            .expect("NSS least-squares solve");
        Self {
            beta0: T::from_f64(b[0]),
            beta1: T::from_f64(b[1]),
            beta2: T::from_f64(b[2]),
            beta3: T::from_f64(b[3]),
            tau1: T::from_f64(tau1),
            tau2: T::from_f64(tau2),
        }
    }
}
