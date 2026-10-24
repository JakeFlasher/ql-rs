use ql_ad::Real;
use ql_math::{Interpolator1D, MonotoneCubicInterpolator};
use ql_time::{Date, DayCounter};

pub trait YieldTermStructure<T: Real = f64>: Send + Sync {
    fn reference_date(&self) -> Date;
    fn discount_factor(&self, t: f64) -> T;
    fn zero_rate(&self, t: f64) -> T;
    fn instantaneous_forward(&self, t: f64) -> T;
    fn forward_rate(&self, d1: Date, d2: Date, dc: DayCounter) -> T {
        // Use the caller's day-count convention for *both* curve lookup
        // times and the accrual denominator.  The old hard-coded `/365`
        // silently mis-priced ACT/360 and 30/360 forwards.
        let t1 = dc.year_fraction(self.reference_date(), d1, None, None);
        let t2 = dc.year_fraction(self.reference_date(), d2, None, None);
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
pub type PiecewiseCubicZeroCurve<T = f64> = MonotoneCubicZeroCurve<T>;
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

/// Hagan--West monotone-convex interpolation of continuously compounded zero
/// rates.  The implementation stores the discrete forwards and nodal forwards
/// explicitly; integrating the local polynomial therefore recovers every input
/// discount factor exactly (up to the final exponential evaluation).
#[derive(Debug, Clone)]
pub struct MonotoneConvexZeroCurve<T: Real = f64> {
    pub reference_date: Date,
    pub times: Vec<f64>,
    pub zero_rates: Vec<T>,
    pub discrete_forwards: Vec<T>,
    pub nodal_forwards: Vec<T>,
    pub enforce_positivity: bool,
}
impl<T: Real> MonotoneConvexZeroCurve<T> {
    pub fn new(reference_date: Date, times: Vec<f64>, zero_rates: Vec<T>) -> Self {
        Self::with_positivity(reference_date, times, zero_rates, true)
    }
    pub fn with_positivity(
        reference_date: Date,
        times: Vec<f64>,
        zero_rates: Vec<T>,
        enforce_positivity: bool,
    ) -> Self {
        assert!(!times.is_empty() && times.len() == zero_rates.len());
        assert!(times[0] > 0.0 && times.windows(2).all(|w| w[1] > w[0]));
        let mut fd = Vec::with_capacity(times.len());
        let mut prev_t = 0.0;
        let mut prev_rt = T::zero();
        for (t, r) in times.iter().copied().zip(zero_rates.iter().copied()) {
            fd.push((r * t - prev_rt) / (t - prev_t));
            prev_t = t;
            prev_rt = r * t;
        }
        let mut f = vec![T::zero(); fd.len() + 1];
        if fd.len() == 1 {
            f[0] = fd[0];
            f[1] = fd[0];
        } else {
            for i in 0..fd.len() - 1 {
                let d0 = if i == 0 {
                    times[0]
                } else {
                    times[i] - times[i - 1]
                };
                let d1 = times[i + 1] - times[i];
                f[i + 1] = (fd[i + 1] * d0 + fd[i] * d1) / (d0 + d1);
            }
            f[0] = fd[0] - (f[1] - fd[0]) * 0.5;
            let n = fd.len();
            f[n] = fd[n - 1] - (f[n - 1] - fd[n - 1]) * 0.5;
        }
        if enforce_positivity {
            for i in 1..fd.len() {
                let cap = (fd[i - 1].min_val(fd[i]) * 2.0).max_val(T::zero());
                f[i] = f[i].max_val(T::zero()).min_val(cap);
            }
            f[0] = f[0]
                .max_val(T::zero())
                .min_val((fd[0] * 2.0).max_val(T::zero()));
            let n = fd.len();
            f[n] = f[n]
                .max_val(T::zero())
                .min_val((fd[n - 1] * 2.0).max_val(T::zero()));
        }
        Self {
            reference_date,
            times,
            zero_rates,
            discrete_forwards: fd,
            nodal_forwards: f,
            enforce_positivity,
        }
    }
    pub fn to_ad_curve<S: Real>(&self, zero_rates: &[S]) -> MonotoneConvexZeroCurve<S> {
        MonotoneConvexZeroCurve::with_positivity(
            self.reference_date,
            self.times.clone(),
            zero_rates.to_vec(),
            self.enforce_positivity,
        )
    }
    fn integrated(&self, t: f64) -> T {
        if t <= 0.0 {
            return T::zero();
        }
        let n = self.times.len();
        if t >= self.times[n - 1] {
            return self.zero_rates[n - 1] * self.times[n - 1]
                + self.nodal_forwards[n] * (t - self.times[n - 1]);
        }
        // Include a knot in the interval ending at that knot so that the
        // integrated curve reproduces every input discount exactly.
        // `i` denotes the right-hand knot of the active interval.  The
        // partition point counts knots strictly below `t`; advance once so a
        // point immediately after knot j uses the (j,j+1) interval rather
        // than clamping to the interval that ended at j.
        let i = self
            .times
            .partition_point(|x| *x < t)
            .saturating_add(1)
            .min(n)
            .max(1);
        let t0 = if i == 1 { 0.0 } else { self.times[i - 2] };
        let t1 = self.times[i - 1];
        let dt = t1 - t0;
        let x = ((t - t0) / dt).clamp(0.0, 1.0);
        let g0 = self.nodal_forwards[i - 1] - self.discrete_forwards[i - 1];
        let g1 = self.nodal_forwards[i] - self.discrete_forwards[i - 1];
        // Integrating the Hagan--West correction from 0 to x gives
        // g0*x*(1-x)^2 - g1*x^2*(1-x).  The correction integrates to zero
        // at x=1, so every input discount factor is reproduced exactly.
        let gx = g0 * x * (1.0 - x) * (1.0 - x) - g1 * x * x * (1.0 - x);
        let base = if i == 1 {
            T::zero()
        } else {
            self.zero_rates[i - 2] * t0
        };
        base + self.discrete_forwards[i - 1] * (t - t0) + T::from_f64(dt) * gx
    }
}
impl<T: Real> YieldTermStructure<T> for MonotoneConvexZeroCurve<T> {
    fn reference_date(&self) -> Date {
        self.reference_date
    }
    fn discount_factor(&self, t: f64) -> T {
        (-self.integrated(t)).exp()
    }
    fn zero_rate(&self, t: f64) -> T {
        if t.abs() < 1e-14 {
            T::zero()
        } else {
            self.integrated(t) / t
        }
    }
    fn instantaneous_forward(&self, t: f64) -> T {
        if t <= 0.0 {
            return self.nodal_forwards[0];
        }
        let i = self
            .times
            .partition_point(|x| *x < t)
            .saturating_add(1)
            .min(self.times.len())
            .max(1);
        let t0 = if i == 1 { 0.0 } else { self.times[i - 2] };
        let dt = self.times[i - 1] - t0;
        let x = ((t - t0) / dt).clamp(0.0, 1.0);
        let g0 = self.nodal_forwards[i - 1] - self.discrete_forwards[i - 1];
        let g1 = self.nodal_forwards[i] - self.discrete_forwards[i - 1];
        self.discrete_forwards[i - 1]
            + g0 * (1.0 - 4.0 * x + 3.0 * x * x)
            + g1 * (-2.0 * x + 3.0 * x * x)
    }
}

#[derive(Debug, Clone)]
pub struct MeetingDateJumpOisCurve<T: Real = f64> {
    pub reference_date: Date,
    pub meeting_times: Vec<f64>,
    pub initial_rate: T,
    pub meeting_jumps: Vec<T>,
}
impl<T: Real> MeetingDateJumpOisCurve<T> {
    pub fn overnight_rate_after_meeting(&self, idx: usize) -> T {
        self.initial_rate
            + self
                .meeting_jumps
                .iter()
                .take(idx + 1)
                .copied()
                .fold(T::zero(), |a, z| a + z)
    }
    pub fn to_ad_curve<S: Real>(&self, initial_rate: S, jumps: &[S]) -> MeetingDateJumpOisCurve<S> {
        MeetingDateJumpOisCurve {
            reference_date: self.reference_date,
            meeting_times: self.meeting_times.clone(),
            initial_rate,
            meeting_jumps: jumps.to_vec(),
        }
    }
    fn integrated_rate(&self, t: f64) -> T {
        let mut last = 0.0;
        let mut rate = self.initial_rate;
        let mut out = T::zero();
        for (i, &meeting) in self.meeting_times.iter().enumerate() {
            let end = t.min(meeting);
            if end > last {
                out += rate * (end - last);
            }
            if t <= meeting {
                return out;
            }
            last = meeting;
            rate += self.meeting_jumps.get(i).copied().unwrap_or(T::zero());
        }
        if t > last {
            out + rate * (t - last)
        } else {
            out
        }
    }
}
impl<T: Real> YieldTermStructure<T> for MeetingDateJumpOisCurve<T> {
    fn reference_date(&self) -> Date {
        self.reference_date
    }
    fn discount_factor(&self, t: f64) -> T {
        (-self.integrated_rate(t)).exp()
    }
    fn zero_rate(&self, t: f64) -> T {
        if t.abs() < 1e-14 {
            T::zero()
        } else {
            self.integrated_rate(t) / t
        }
    }
    fn instantaneous_forward(&self, t: f64) -> T {
        for (i, &m) in self.meeting_times.iter().enumerate() {
            if t < m {
                return if i == 0 {
                    self.initial_rate
                } else {
                    self.overnight_rate_after_meeting(i - 1)
                };
            }
        }
        self.initial_rate
            + self
                .meeting_jumps
                .iter()
                .copied()
                .fold(T::zero(), |a, z| a + z)
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
