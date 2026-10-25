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
    pub ameliorate_slopes: bool,
}

/// Specification for the coupled Hagan--West par-swap bootstrap.
#[derive(Clone, Debug)]
pub struct HaganWestMonotoneConvexBootstrapSpec<T: Real = f64> {
    pub reference_date: Date,
    pub pillar_times: Vec<f64>,
    pub par_swap_quotes: Vec<T>,
    pub payment_frequency_per_year: usize,
    pub enforce_positivity: bool,
    pub ameliorate_slopes: bool,
    pub query_time: f64,
}
#[derive(Clone, Debug)]
pub struct HaganWestMonotoneConvexReport<T: Real = f64> {
    pub curve: MonotoneConvexZeroCurve<T>,
    pub calibrated_zero_rates: Vec<T>,
    pub discrete_forwards: Vec<T>,
    pub nodal_forwards: Vec<T>,
    pub query_discount_factor: T,
    pub query_instantaneous_forward: T,
    pub min_instantaneous_forward: T,
    pub forward_roughness_l2: T,
    pub max_par_swap_residual: T,
}

/// Exact four-region Hagan--West correction and its integral.  Region
/// selection is made from the primal values, while all arithmetic remains in
/// `T`, so Dual tangents are retained on every smooth branch.
pub fn hagan_west_4region_g_and_int<T: Real>(g0: T, g1: T, x: T) -> (T, T) {
    let x = x.max_val(T::zero()).min_val(T::one());
    let p = (g0 * g1).val();
    if p >= 0.0 {
        if (g0 + g1).val().abs() <= 1e-14 {
            return (T::zero(), T::zero());
        }
        let eta = (g1 / (g0 + g1))
            .max_val(T::from_f64(1e-14))
            .min_val(T::from_f64(1.0 - 1e-14));
        let a = -(g0 * g1) / (g0 + g1);
        if x.val() <= eta.val() {
            let q = (eta - x) / eta;
            (
                a + (g0 - a) * q * q,
                a * x + eta * (g0 - a) * (T::one() - q * q * q) / 3.0,
            )
        } else {
            let q = (x - eta) / (T::one() - eta);
            (
                a + (g1 - a) * q * q,
                a * x + eta * (g0 - a) / 3.0 + (T::one() - eta) * (g1 - a) * q * q * q / 3.0,
            )
        }
    } else if g1.abs().val() > (g0 * 2.0).abs().val() {
        let eta = ((g1 + g0 * 2.0) / (g1 - g0))
            .max_val(T::zero())
            .min_val(T::one());
        if x.val() <= eta.val() {
            (g0, g0 * x)
        } else {
            let q = (x - eta) / (T::one() - eta);
            (
                g0 + (g1 - g0) * q * q,
                g0 * x + (g1 - g0) * (T::one() - eta) * q * q * q / 3.0,
            )
        }
    } else if g0.abs().val() > 2.0 * g1.abs().val() {
        let eta = (g1 * 3.0 / (g1 - g0)).max_val(T::zero()).min_val(T::one());
        if x.val() < eta.val() {
            let q = (eta - x) / eta;
            (
                g1 + (g0 - g1) * q * q,
                g1 * x + eta * (g0 - g1) * (T::one() - q * q * q) / 3.0,
            )
        } else {
            (g1, g1 * x + eta * (g0 - g1) / 3.0)
        }
    } else {
        let g = g0 * (T::one() - x * 4.0 + x * x * 3.0) + g1 * (-x * 2.0 + x * x * 3.0);
        let g_int = x * (T::one() - x) * (g0 * (T::one() - x) - g1 * x);
        (g, g_int)
    }
}

pub fn bootstrap_hagan_west_monotone_convex_curve<T: Real>(
    spec: &HaganWestMonotoneConvexBootstrapSpec<T>,
) -> HaganWestMonotoneConvexReport<T> {
    assert_eq!(spec.pillar_times.len(), spec.par_swap_quotes.len());
    assert!(!spec.pillar_times.is_empty());
    let n = spec.pillar_times.len();
    let freq = spec.payment_frequency_per_year.max(1) as f64;
    // Solve the primal calibration entirely in f64 first.  Keeping this
    // fixed-point/Newton phase out of the AD domain avoids constructing and
    // propagating Dual vectors during all twenty convergence sweeps.
    let mut r_f64 = vec![0.0_f64; n];
    // Sequential piecewise-flat initialisation.
    for k in 0..n {
        let t = spec.pillar_times[k];
        let s = spec.par_swap_quotes[k].val();
        let mut ann = 0.0_f64;
        let mut j = 1usize;
        while (j as f64) / freq < t - 1e-12 {
            let u = j as f64 / freq;
            ann += (1.0 / freq) * (-r_f64[j.saturating_sub(1).min(k)] * u).exp();
            j += 1;
        }
        let df = (1.0 - s * ann) / (1.0 + s / freq);
        r_f64[k] = (-df.max(1e-12).ln() / t).max(-10.0);
    }
    let residuals_f64 = |rr: &[f64]| -> Vec<f64> {
        let c = MonotoneConvexZeroCurve::with_options(
            spec.reference_date,
            spec.pillar_times.clone(),
            rr.to_vec(),
            spec.enforce_positivity,
            spec.ameliorate_slopes,
        );
        spec.pillar_times
            .iter()
            .enumerate()
            .map(|(k, &tk)| {
                let pv = c.discount_factor(tk);
                let mut u = 1.0 / freq;
                let mut ann = 0.0_f64;
                while u <= tk + 1e-10 {
                    ann += (1.0 / freq) * c.discount_factor(u);
                    u += 1.0 / freq;
                }
                pv + spec.par_swap_quotes[k].val() * ann - 1.0
            })
            .collect()
    };
    let mut primal_diag = vec![0.0_f64; n];
    // Fourteen coupled Gauss--Seidel sweeps using the analytical diagonal
    // slope.  Each pillar sees the rates updated by all earlier pillars in
    // the same sweep.
    for _ in 0..14 {
        for k in 0..n {
            let rs_k = residuals_f64(&r_f64)[k];
            let c = MonotoneConvexZeroCurve::with_options(
                spec.reference_date,
                spec.pillar_times.clone(),
                r_f64.clone(),
                spec.enforce_positivity,
                spec.ameliorate_slopes,
            );
            let tk = spec.pillar_times[k];
            let mut ann = 0.0_f64;
            let mut u = 1.0 / freq;
            while u <= tk + 1e-10 {
                if u > spec
                    .pillar_times
                    .get(k.wrapping_sub(1))
                    .copied()
                    .unwrap_or(0.0)
                    + 1e-12
                {
                    ann += (1.0 / freq) * u * c.discount_factor(u);
                }
                u += 1.0 / freq;
            }
            let der =
                -(tk * c.discount_factor(tk) + spec.par_swap_quotes[k].val() * ann).max(1e-10);
            r_f64[k] -= rs_k / der;
        }
    }
    // Six primal secant-Newton Gauss--Seidel polishing sweeps.  The secant
    // bump is deliberately f64-only; no T/Dual temporary vectors are formed.
    for _ in 0..6 {
        for k in 0..n {
            let rs_k = residuals_f64(&r_f64)[k];
            let h = 1e-6_f64;
            let mut bumped = r_f64.clone();
            bumped[k] += h;
            let der = ((residuals_f64(&bumped)[k] - rs_k) / h).min(-1e-10);
            primal_diag[k] = der;
            r_f64[k] -= rs_k / der;
        }
    }
    let mut rates: Vec<T> = r_f64.iter().map(|&r| T::from_f64(r)).collect();
    let residuals = |rr: &[T]| -> Vec<T> {
        let c = MonotoneConvexZeroCurve::with_options(
            spec.reference_date,
            spec.pillar_times.clone(),
            rr.to_vec(),
            spec.enforce_positivity,
            spec.ameliorate_slopes,
        );
        spec.pillar_times
            .iter()
            .enumerate()
            .map(|(k, &tk)| {
                let pv = c.discount_factor(tk);
                let mut u = 1.0 / freq;
                let mut ann = T::zero();
                while u <= tk + 1e-10 {
                    ann += T::from_f64(1.0 / freq) * c.discount_factor(u);
                    u += 1.0 / freq;
                }
                pv + spec.par_swap_quotes[k] * ann - T::one()
            })
            .collect()
    };
    // Two analytical T-domain Newton sweeps propagate Dual tangents through
    // the converged primal solution without any finite-difference bump.
    for _ in 0..2 {
        for k in 0..n {
            let rs_k = residuals(&rates)[k];
            let c = MonotoneConvexZeroCurve::with_options(
                spec.reference_date,
                spec.pillar_times.clone(),
                rates.clone(),
                spec.enforce_positivity,
                spec.ameliorate_slopes,
            );
            let tk = spec.pillar_times[k];
            let mut ann = T::zero();
            let mut u = 1.0 / freq;
            while u <= tk + 1e-10 {
                if u > spec
                    .pillar_times
                    .get(k.wrapping_sub(1))
                    .copied()
                    .unwrap_or(0.0)
                    + 1e-12
                {
                    ann += T::from_f64(1.0 / freq) * u * c.discount_factor(u);
                }
                u += 1.0 / freq;
            }
            let analytic_der = -(T::from_f64(tk) * c.discount_factor(tk)
                + spec.par_swap_quotes[k] * ann)
                .max_val(T::from_f64(1e-10));
            // The primal secant polish supplies a scalar correction for the
            // interpolation/clamp branches.  Keep the Newton denominator
            // scalar so its own derivative is not fed back into the update;
            // all quote tangents continue to flow through the T residual.
            let der =
                T::from_f64(analytic_der.val() * (primal_diag[k] / analytic_der.val().min(-1e-10)));
            rates[k] -= rs_k / der;
        }
    }
    let curve = MonotoneConvexZeroCurve::with_options(
        spec.reference_date,
        spec.pillar_times.clone(),
        rates.clone(),
        spec.enforce_positivity,
        spec.ameliorate_slopes,
    );
    let residuals_now = residuals(&rates);
    let maxr = residuals_now
        .iter()
        .copied()
        .map(|x| x.abs())
        .reduce(|a, b| a.max_val(b))
        .unwrap_or(T::zero());
    let tmax = *spec.pillar_times.last().unwrap();
    let minf = (0..=256)
        .map(|i| curve.instantaneous_forward(tmax * i as f64 / 256.0))
        .reduce(|a, b| a.min_val(b))
        .unwrap_or(T::zero());
    let rough = (curve
        .nodal_forwards
        .windows(2)
        .map(|w| (w[1] - w[0]) * (w[1] - w[0]))
        .fold(T::zero(), |a, b| a + b))
    .sqrt();
    HaganWestMonotoneConvexReport {
        query_discount_factor: curve.discount_factor(spec.query_time),
        query_instantaneous_forward: curve.instantaneous_forward(spec.query_time),
        discrete_forwards: curve.discrete_forwards.clone(),
        nodal_forwards: curve.nodal_forwards.clone(),
        curve,
        calibrated_zero_rates: rates,
        min_instantaneous_forward: minf,
        forward_roughness_l2: rough,
        max_par_swap_residual: maxr,
    }
}
impl<T: Real> MonotoneConvexZeroCurve<T> {
    pub fn new(reference_date: Date, times: Vec<f64>, zero_rates: Vec<T>) -> Self {
        Self::with_options(reference_date, times, zero_rates, true, false)
    }
    pub fn with_positivity(
        reference_date: Date,
        times: Vec<f64>,
        zero_rates: Vec<T>,
        enforce_positivity: bool,
    ) -> Self {
        Self::with_options(reference_date, times, zero_rates, enforce_positivity, false)
    }
    pub fn with_options(
        reference_date: Date,
        times: Vec<f64>,
        zero_rates: Vec<T>,
        enforce_positivity: bool,
        ameliorate_slopes: bool,
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
        }
        if ameliorate_slopes && fd.len() > 1 {
            for i in 1..fd.len() {
                let midpoint = (fd[i - 1] + fd[i]) * 0.5;
                f[i] = midpoint + (f[i] - midpoint) * 0.85;
            }
        }
        if enforce_positivity {
            for i in 1..fd.len() {
                let cap = (fd[i - 1].min_val(fd[i]) * 2.0).max_val(T::zero());
                f[i] = f[i].max_val(T::zero()).min_val(cap);
            }
        }
        if ameliorate_slopes && fd.len() > 1 {
            for i in 1..fd.len() {
                let midpoint = (f[i - 1] + f[i + 1]) * 0.5;
                let target = midpoint + (f[i] - midpoint) * 0.85;
                if enforce_positivity {
                    let cap = (fd[i - 1].min_val(fd[i]) * 2.0).max_val(T::zero());
                    f[i] = target.max_val(T::zero()).min_val(cap);
                } else {
                    f[i] = target;
                }
            }
        }
        // Compute the boundary nodes only after the interior values have
        // been clamped.  Otherwise a positivity clamp can change f[1] or
        // f[n-1] without updating the one-sided Hagan--West extrapolation.
        if fd.len() == 1 {
            f[0] = fd[0];
            f[1] = fd[0];
        } else {
            f[0] = fd[0] - (f[1] - fd[0]) * 0.5;
            let n = fd.len();
            f[n] = fd[n - 1] - (f[n - 1] - fd[n - 1]) * 0.5;
        }
        if enforce_positivity {
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
            ameliorate_slopes,
        }
    }
    pub fn to_ad_curve<S: Real>(&self, zero_rates: &[S]) -> MonotoneConvexZeroCurve<S> {
        MonotoneConvexZeroCurve::with_options(
            self.reference_date,
            self.times.clone(),
            zero_rates.to_vec(),
            self.enforce_positivity,
            self.ameliorate_slopes,
        )
    }
    fn interval_g(&self, i: usize) -> (T, T) {
        let fd = self.discrete_forwards[i - 1];
        let g0 = self.nodal_forwards[i - 1] - fd;
        let g1 = self.nodal_forwards[i] - fd;
        (g0, g1)
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
        let (g0, g1) = self.interval_g(i);
        let (_, gx) = hagan_west_4region_g_and_int(g0, g1, T::from_f64(x));
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
        let (g0, g1) = self.interval_g(i);
        self.discrete_forwards[i - 1] + hagan_west_4region_g_and_int(g0, g1, T::from_f64(x)).0
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

#[cfg(test)]
mod phase13_tests {
    use super::*;
    use ql_time::Date;
    #[test]
    fn monotone_convex_nodes_and_positivity() {
        let c = MonotoneConvexZeroCurve::new(Date(0), vec![1.0, 2.0, 3.0], vec![0.01, 0.05, 0.034]);
        for (t, r) in [(1.0, 0.01), (2.0, 0.05), (3.0, 0.034)] {
            assert!((c.zero_rate(t) - r).abs() < 1e-12);
        }
        for k in 0..=60 {
            let t = k as f64 * 0.05;
            assert!(c.instantaneous_forward(t) >= -1e-12);
        }
    }
}

#[cfg(test)]
mod phase21_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_time::Date;

    #[test]
    fn bootstrap_and_slope_amelioration_invariants() {
        let spec = HaganWestMonotoneConvexBootstrapSpec {
            reference_date: Date(0),
            pillar_times: vec![1.0, 2.0, 5.0],
            par_swap_quotes: vec![0.02, 0.05, 0.03],
            payment_frequency_per_year: 2,
            enforce_positivity: true,
            ameliorate_slopes: true,
            query_time: 2.0,
        };
        let report = bootstrap_hagan_west_monotone_convex_curve(&spec);
        assert!(report.max_par_swap_residual.val() < 1e-12);
        assert!(report.min_instantaneous_forward.val() >= -1e-12);
        let plain = MonotoneConvexZeroCurve::with_options(
            Date(0),
            vec![1.0, 2.0, 3.0, 4.0],
            vec![0.01, 0.08, 0.015, 0.07],
            true,
            false,
        );
        let smooth = MonotoneConvexZeroCurve::with_options(
            Date(0),
            vec![1.0, 2.0, 3.0, 4.0],
            vec![0.01, 0.08, 0.015, 0.07],
            true,
            true,
        );
        let rough = |c: &MonotoneConvexZeroCurve<f64>| {
            c.nodal_forwards
                .windows(2)
                .map(|w| (w[1] - w[0]).powi(2))
                .sum::<f64>()
        };
        assert!(rough(&smooth) < rough(&plain));
        let h = 1e-5;
        let mut up = spec.clone();
        up.par_swap_quotes[1] += h;
        let mut dn = spec.clone();
        dn.par_swap_quotes[1] -= h;
        let fd = (bootstrap_hagan_west_monotone_convex_curve(&up).query_discount_factor
            - bootstrap_hagan_west_monotone_convex_curve(&dn).query_discount_factor)
            / (2.0 * h);
        type D = Dual<1>;
        let ad_spec = HaganWestMonotoneConvexBootstrapSpec::<D> {
            reference_date: Date(0),
            pillar_times: vec![1.0, 2.0, 5.0],
            par_swap_quotes: vec![
                D::from_value(0.02),
                D::variable(0.05, 0),
                D::from_value(0.03),
            ],
            payment_frequency_per_year: 2,
            enforce_positivity: true,
            ameliorate_slopes: true,
            query_time: 2.0,
        };
        let ad = bootstrap_hagan_west_monotone_convex_curve(&ad_spec).query_discount_factor;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
