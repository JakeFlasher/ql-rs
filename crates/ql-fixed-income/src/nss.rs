use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_time::{Date, DayCounter};
#[derive(Clone, Copy, Debug)]
pub struct NelsonSiegelSvensson<T: Real> {
    pub reference_date: Date,
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
        let tau1 = self.tau1.max_val(T::from_f64(1e-12));
        let tau2 = self.tau2.max_val(T::from_f64(1e-12));
        let x1 = t / tau1;
        let x2 = t / tau2;
        let h_1 = Self::h1(x1);
        let h_2 = h_1 - (-x1).exp();
        self.beta0 + self.beta1 * h_1 + self.beta2 * h_2 + self.beta3 * (Self::h1(x2) - (-x2).exp())
    }
    pub fn forward_rate(&self, t: T) -> T {
        let tau1 = self.tau1.max_val(T::from_f64(1e-12));
        let tau2 = self.tau2.max_val(T::from_f64(1e-12));
        let e1 = (-t / tau1).exp();
        let e2 = (-t / tau2).exp();
        self.beta0 + self.beta1 * e1 + self.beta2 * (t / tau1) * e1 + self.beta3 * (t / tau2) * e2
    }
    pub fn discount_factor(&self, t: T) -> T {
        (-self.zero_rate(t) * t).exp()
    }
    pub fn fit_ols_for_taus(times: &[f64], rates: &[f64], tau1: f64, tau2: f64) -> Self {
        use nalgebra::{DMatrix, DVector};
        assert_eq!(
            times.len(),
            rates.len(),
            "NSS times/rates dimension mismatch"
        );
        assert!(!times.is_empty(), "NSS fit requires observations");
        assert!(
            times.iter().chain(rates.iter()).all(|x| x.is_finite()),
            "NSS fit inputs must be finite"
        );
        let tau1 = tau1.max(1e-12);
        let tau2 = tau2.max(1e-12);
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
            reference_date: Date(0),
            beta0: T::from_f64(b[0]),
            beta1: T::from_f64(b[1]),
            beta2: T::from_f64(b[2]),
            beta3: T::from_f64(b[3]),
            tau1: T::from_f64(tau1),
            tau2: T::from_f64(tau2),
        }
    }

    pub fn with_reference_date(mut self, reference_date: Date) -> Self {
        self.reference_date = reference_date;
        self
    }
}

impl<T: Real> YieldTermStructure<T> for NelsonSiegelSvensson<T> {
    fn reference_date(&self) -> Date {
        self.reference_date
    }
    fn discount_factor(&self, t: f64) -> T {
        NelsonSiegelSvensson::discount_factor(self, T::from_f64(t))
    }
    fn zero_rate(&self, t: f64) -> T {
        NelsonSiegelSvensson::zero_rate(self, T::from_f64(t))
    }
    fn instantaneous_forward(&self, t: f64) -> T {
        NelsonSiegelSvensson::forward_rate(self, T::from_f64(t))
    }
    fn forward_rate(&self, d1: Date, d2: Date, dc: DayCounter) -> T {
        let t1 = dc.year_fraction(self.reference_date, d1, None, None);
        let t2 = dc.year_fraction(self.reference_date, d2, None, None);
        let tau = dc.year_fraction(d1, d2, None, None);
        (<Self as YieldTermStructure<T>>::discount_factor(self, t1)
            / <Self as YieldTermStructure<T>>::discount_factor(self, t2)
            - T::one())
            / tau
    }
}
