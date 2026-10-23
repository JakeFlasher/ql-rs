use ql_ad::Real;
use ql_math::{Interpolator1D, MonotoneCubicInterpolator};
use ql_time::{Date, DayCounter};

pub trait YieldTermStructure<T: Real = f64>: Send + Sync {
    fn reference_date(&self) -> Date;
    fn discount_factor(&self, t: f64) -> T;
    fn zero_rate(&self, t: f64) -> T;
    fn instantaneous_forward(&self, t: f64) -> T;
    fn forward_rate(&self, d1: Date, d2: Date, dc: DayCounter) -> T {
        let t1 = (d1 - self.reference_date()) as f64 / 365.0;
        let t2 = (d2 - self.reference_date()) as f64 / 365.0;
        let tau = dc.year_fraction(d1, d2, None, None);
        (self.discount_factor(t1) / self.discount_factor(t2) - T::one()) / tau
    }
}
#[derive(Debug, Clone)]
pub struct FlatForwardCurve<T: Real = f64> {
    pub reference_date: Date,
    pub rate: T,
    pub day_counter: DayCounter,
}
impl<T: Real> FlatForwardCurve<T> {
    pub fn new(reference_date: Date, rate: T) -> Self {
        Self {
            reference_date,
            rate,
            day_counter: DayCounter::Act365Fixed,
        }
    }
    pub fn with_day_counter(reference_date: Date, rate: T, day_counter: DayCounter) -> Self {
        Self {
            reference_date,
            rate,
            day_counter,
        }
    }
    pub fn from_rate(reference_date: Date, rate: T) -> Self {
        Self::new(reference_date, rate)
    }
}
impl<T: Real> YieldTermStructure<T> for FlatForwardCurve<T> {
    fn reference_date(&self) -> Date {
        self.reference_date
    }
    fn discount_factor(&self, t: f64) -> T {
        (-self.rate * t).exp()
    }
    fn zero_rate(&self, _: f64) -> T {
        self.rate
    }
    fn instantaneous_forward(&self, _: f64) -> T {
        self.rate
    }
}
#[derive(Debug, Clone)]
pub struct DiscountCurve<T: Real = f64> {
    pub reference_date: Date,
    pub times: Vec<f64>,
    pub discounts: Vec<T>,
}
impl<T: Real> DiscountCurve<T> {
    pub fn new(reference_date: Date, times: Vec<f64>, discounts: Vec<T>) -> Self {
        assert_eq!(times.len(), discounts.len());
        assert!(!times.is_empty());
        assert!(times.iter().all(|t| *t > 0.0));
        assert!(times.windows(2).all(|w| w[1] > w[0]));
        assert!(discounts.iter().all(|d| d.val() > 0.0));
        Self {
            reference_date,
            times,
            discounts,
        }
    }
    pub fn from_nodes(reference_date: Date, times: Vec<f64>, discounts: Vec<T>) -> Self {
        Self::new(reference_date, times, discounts)
    }
    pub fn to_ad_curve<S: Real>(&self, node_vars: &[S]) -> DiscountCurve<S> {
        assert_eq!(node_vars.len(), self.times.len());
        DiscountCurve::new(self.reference_date, self.times.clone(), node_vars.to_vec())
    }
}
impl<T: Real> YieldTermStructure<T> for DiscountCurve<T> {
    fn reference_date(&self) -> Date {
        self.reference_date
    }
    fn discount_factor(&self, t: f64) -> T {
        if t <= 0.0 {
            return T::one();
        }
        if t <= self.times[0] {
            // log-linear interpolation from anchor D(0)=1 to first node
            let w = t / self.times[0];
            return (self.discounts[0].ln() * T::from_f64(w)).exp();
        }
        if t >= *self.times.last().unwrap() {
            let n = self.times.len();
            if n == 1 {
                return (self.discounts[0].ln() * T::from_f64(t / self.times[0])).exp();
            }
            let slope = (self.discounts[n - 1].ln() - self.discounts[n - 2].ln())
                / (self.times[n - 1] - self.times[n - 2]);
            return (self.discounts[n - 1].ln() + slope * T::from_f64(t - self.times[n - 1])).exp();
        }
        let i = self.times.partition_point(|x| *x <= t).saturating_sub(1);
        let l = self.discounts[i].ln();
        let r = self.discounts[i + 1].ln();
        (l + (r - l) * ((t - self.times[i]) / (self.times[i + 1] - self.times[i]))).exp()
    }
    fn zero_rate(&self, t: f64) -> T {
        if t.abs() < 1e-14 {
            return T::zero();
        }
        -self.discount_factor(t).ln() / t
    }
    fn instantaneous_forward(&self, t: f64) -> T {
        if self.times.len() == 1 {
            // With a single node the curve is log-linear from the anchor and
            // flat-forward extrapolated thereafter, so the instantaneous
            // forward is the node's continuously compounded zero rate at all
            // strictly positive maturities.
            return -self.discounts[0].ln() / self.times[0];
        }
        if t <= self.times[0] {
            return -self.discounts[0].ln() / self.times[0];
        }
        if t >= *self.times.last().unwrap() {
            let n = self.times.len();
            return -(self.discounts[n - 1].ln() - self.discounts[n - 2].ln())
                / (self.times[n - 1] - self.times[n - 2]);
        }
        let i = self.times.partition_point(|x| *x <= t).saturating_sub(1);
        -(self.discounts[i + 1].ln() - self.discounts[i].ln()) / (self.times[i + 1] - self.times[i])
    }
}
#[derive(Debug, Clone)]
pub struct ZeroCurve<T: Real = f64> {
    pub reference_date: Date,
    pub times: Vec<f64>,
    pub zero_rates: Vec<T>,
}
impl<T: Real> ZeroCurve<T> {
    pub fn new(reference_date: Date, times: Vec<f64>, zero_rates: Vec<T>) -> Self {
        assert_eq!(times.len(), zero_rates.len());
        assert!(!times.is_empty());
        assert!(times.iter().all(|t| *t > 0.0));
        assert!(times.windows(2).all(|w| w[1] > w[0]));
        Self {
            reference_date,
            times,
            zero_rates,
        }
    }
    pub fn to_ad_curve<S: Real>(&self, node_vars: &[S]) -> ZeroCurve<S> {
        ZeroCurve::new(self.reference_date, self.times.clone(), node_vars.to_vec())
    }
}
impl<T: Real> YieldTermStructure<T> for ZeroCurve<T> {
    fn reference_date(&self) -> Date {
        self.reference_date
    }
    fn zero_rate(&self, t: f64) -> T {
        linear_slice_eval(&self.times, &self.zero_rates, t)
    }
    fn discount_factor(&self, t: f64) -> T {
        (-self.zero_rate(t) * t).exp()
    }
    fn instantaneous_forward(&self, t: f64) -> T {
        let (z, dz) = linear_slice_value_derivative(&self.times, &self.zero_rates, t);
        z + T::from_f64(t) * dz
    }
}
pub type FlatYieldCurve<T = f64> = FlatForwardCurve<T>;
pub type PiecewiseLogLinearDiscountCurve<T = f64> = DiscountCurve<T>;
pub type PiecewiseCubicZeroCurve<T = f64> = ZeroCurve<T>;
#[derive(Debug, Clone)]
pub struct MonotoneCubicZeroCurve<T: Real = f64> {
    pub reference_date: Date,
    pub interpolator: MonotoneCubicInterpolator<T>,
}
impl<T: Real> MonotoneCubicZeroCurve<T> {
    pub fn new(reference_date: Date, times: Vec<f64>, zero_rates: Vec<T>) -> Self {
        Self {
            reference_date,
            interpolator: MonotoneCubicInterpolator::new(times, zero_rates),
        }
    }
    pub fn to_ad_curve<S: Real>(&self, node_vars: &[S]) -> MonotoneCubicZeroCurve<S> {
        MonotoneCubicZeroCurve {
            reference_date: self.reference_date,
            interpolator: MonotoneCubicInterpolator::new(
                self.interpolator.x.clone(),
                node_vars.to_vec(),
            ),
        }
    }
}
impl<T: Real> YieldTermStructure<T> for MonotoneCubicZeroCurve<T> {
    fn reference_date(&self) -> Date {
        self.reference_date
    }
    fn zero_rate(&self, t: f64) -> T {
        self.interpolator.eval(T::from_f64(t))
    }
    fn discount_factor(&self, t: f64) -> T {
        (-self.zero_rate(t) * t).exp()
    }
    fn instantaneous_forward(&self, t: f64) -> T {
        let z = self.zero_rate(t);
        z + T::from_f64(t) * self.interpolator.derivative(T::from_f64(t))
    }
}
fn linear_slice_eval<T: Real>(xs: &[f64], ys: &[T], t: f64) -> T {
    linear_slice_value_derivative(xs, ys, t).0
}
fn linear_slice_value_derivative<T: Real>(xs: &[f64], ys: &[T], t: f64) -> (T, T) {
    if xs.len() == 1 {
        return (ys[0], T::zero());
    }
    let i = if t <= xs[0] {
        0
    } else if t >= xs[xs.len() - 1] {
        xs.len() - 2
    } else {
        xs.partition_point(|x| *x <= t).saturating_sub(1)
    };
    let h = xs[i + 1] - xs[i];
    (
        ys[i] + (T::from_f64(t - xs[i]) / h) * (ys[i + 1] - ys[i]),
        (ys[i + 1] - ys[i]) / h,
    )
}
