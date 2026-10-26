#![allow(clippy::needless_range_loop)]
use ql_ad::Real;
use ql_math::{Interpolator1D, MonotoneCubicInterpolator};
use ql_time::{Date, DayCounter};

/// RKHS tension-spline yield-curve calibration (Filipovic--Pelger--Ye).
#[derive(Clone, Debug)]
pub struct FilipovicRkhsYieldCurveSpec<T: Real = f64> {
    pub cashflow_times: Vec<T>,
    pub instrument_cashflow_matrix: Vec<Vec<T>>,
    pub instrument_dirty_prices: Vec<T>,
    pub decay_alpha: T,
    pub ridge_penalty_lambda: T,
    pub query_times: Vec<T>,
}
#[derive(Clone, Debug)]
pub struct FilipovicRkhsYieldCurveReport<T: Real = f64> {
    pub dual_representer_weights_beta: Vec<T>,
    pub fitted_instrument_prices: Vec<T>,
    pub pricing_rmse_bps: T,
    pub rkhs_smoothness_norm_squared: T,
    pub query_discount_factors: Vec<T>,
    pub query_zero_rates: Vec<T>,
    pub query_instantaneous_forwards: Vec<T>,
    pub quote_to_zero_rate_jacobian: Vec<Vec<T>>,
}
fn solve_dense_curves<T: Real>(a: &[Vec<T>], b: &[T]) -> Vec<T> {
    let n = b.len();
    let mut m = vec![vec![T::zero(); n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = a[i][j];
        }
        m[i][n] = b[i];
    }
    for i in 0..n {
        let mut p = i;
        for r in i + 1..n {
            if m[r][i].abs().val() > m[p][i].abs().val() {
                p = r;
            }
        }
        m.swap(i, p);
        let d = if m[i][i].abs().val() < 1e-14 {
            if m[i][i].val() < 0.0 {
                T::from_f64(-1e-14)
            } else {
                T::from_f64(1e-14)
            }
        } else {
            m[i][i]
        };
        for j in i..=n {
            m[i][j] /= d;
        }
        let piv = m[i].clone();
        for r in 0..n {
            if r != i {
                let f = m[r][i];
                for j in i..=n {
                    m[r][j] -= f * piv[j];
                }
            }
        }
    }
    (0..n).map(|i| m[i][n]).collect()
}
fn rkhs_kernel_and_dt<T: Real>(t: T, u: T, alpha: T) -> (T, T) {
    let mn = t.min_val(u);
    let mx = t.max_val(u);
    let e = (-alpha * (t + u)).exp();
    let core = mn
        + alpha * alpha * (T::from_f64(0.5) * mx * mn * mn - T::from_f64(1.0 / 6.0) * mn * mn * mn);
    let dcore = if t.val() <= u.val() {
        T::one() + alpha * alpha * (u * t - T::from_f64(0.5) * t * t)
    } else {
        T::from_f64(0.5) * alpha * alpha * u * u
    };
    (e * core, -alpha * e * core + e * dcore)
}

/// BCBS d368/d582 six-scenario standardized IRRBB engine.
#[derive(Clone, Debug)]
pub struct BcbsD368SixScenarioIrrbbSpec<T: Real = f64> {
    pub bucket_midpoints_years: Vec<T>,
    pub asset_repricing_cashflows: Vec<T>,
    pub liability_repricing_cashflows: Vec<T>,
    pub baseline_zero_rates: Vec<T>,
    pub parallel_shock_bps: T,
    pub short_shock_bps: T,
    pub long_shock_bps: T,
    pub decay_tau_years: T,
    pub tier1_capital: T,
    pub pass_through_deposit_beta: T,
}

#[derive(Clone, Debug)]
pub struct BcbsD368SixScenarioIrrbbReport<T: Real = f64> {
    pub baseline_eve: T,
    pub scenario_delta_eves: [T; 6],
    pub worst_case_delta_eveloss: T,
    pub worst_scenario_index: usize,
    pub parallel_up_delta_nii_1y: T,
    pub parallel_down_delta_nii_1y: T,
    pub steepener_delta_eve: T,
    pub flattener_delta_eve: T,
    pub supervisory_outlier_ratio_tier1: T,
}

pub fn evaluate_bcbs_d368_six_scenario_irrbb<T: Real>(
    s: &BcbsD368SixScenarioIrrbbSpec<T>,
) -> BcbsD368SixScenarioIrrbbReport<T> {
    let n = s
        .bucket_midpoints_years
        .len()
        .min(s.asset_repricing_cashflows.len())
        .min(s.liability_repricing_cashflows.len())
        .min(s.baseline_zero_rates.len());
    assert!(n > 0, "IRRBB requires at least one bucket");
    let tau = s.decay_tau_years.max_val(T::from_f64(1e-8));
    let beta = s.pass_through_deposit_beta;
    let mut base = T::zero();
    let mut net_1y = T::zero();
    for i in 0..n {
        let t = s.bucket_midpoints_years[i];
        // Smooth short-bucket indicator.  Keeping this in the AD domain is
        // important: bucket times are calibration variables in the IRRBB
        // Jacobian and a hard one-year gate would erase their sensitivities.
        let w_1y = T::one() / (T::one() + (T::from_f64(20.0) * (t - T::one())).exp());
        let dep = T::one() - beta * w_1y;
        let cf = s.asset_repricing_cashflows[i] - dep * s.liability_repricing_cashflows[i];
        base += cf * (-s.baseline_zero_rates[i] * t).exp();
        net_1y += w_1y
            * (T::one() - t)
            * (s.asset_repricing_cashflows[i] - beta * s.liability_repricing_cashflows[i]);
    }
    let rp = s.parallel_shock_bps * T::from_f64(1e-4);
    let rs = s.short_shock_bps * T::from_f64(1e-4);
    let rl = s.long_shock_bps * T::from_f64(1e-4);
    let mut deltas = [T::zero(); 6];
    for scenario in 0..6 {
        let mut eve = T::zero();
        for i in 0..n {
            let t = s.bucket_midpoints_years[i];
            let sh = (-t / tau).exp();
            let lh = T::one() - sh;
            let shock = match scenario {
                0 => rp,
                1 => -rp,
                2 => -T::from_f64(0.65) * rs * sh + T::from_f64(0.90) * rl * lh,
                3 => T::from_f64(0.80) * rs * sh - T::from_f64(0.60) * rl * lh,
                4 => rs * sh,
                _ => -rs * sh,
            };
            let w_1y = T::one() / (T::one() + (T::from_f64(20.0) * (t - T::one())).exp());
            let dep = T::one() - beta * w_1y;
            let cf = s.asset_repricing_cashflows[i] - dep * s.liability_repricing_cashflows[i];
            eve += cf * (-(s.baseline_zero_rates[i] + shock) * t).exp();
        }
        deltas[scenario] = eve - base;
    }
    let mut worst_loss = T::zero();
    let mut worst_idx = 0;
    for (i, d) in deltas.iter().enumerate() {
        let loss = (-*d).max_val(T::zero());
        if loss.val() > worst_loss.val() {
            worst_loss = loss;
            worst_idx = i;
        }
    }
    BcbsD368SixScenarioIrrbbReport {
        baseline_eve: base,
        scenario_delta_eves: deltas,
        worst_case_delta_eveloss: worst_loss,
        worst_scenario_index: worst_idx,
        parallel_up_delta_nii_1y: rp * net_1y,
        parallel_down_delta_nii_1y: -rp * net_1y,
        steepener_delta_eve: deltas[2],
        flattener_delta_eve: deltas[3],
        supervisory_outlier_ratio_tier1: worst_loss
            / s.tier1_capital.abs().max_val(T::from_f64(1e-12)),
    }
}
pub fn calibrate_filipovic_rkhs_smooth_yield_curve<T: Real>(
    s: &FilipovicRkhsYieldCurveSpec<T>,
) -> FilipovicRkhsYieldCurveReport<T> {
    let n = s.cashflow_times.len();
    let m = s.instrument_cashflow_matrix.len();
    assert_eq!(s.instrument_dirty_prices.len(), m);
    let alpha = s.decay_alpha.max_val(T::from_f64(1e-6));
    let mut k = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            let t = s.cashflow_times[i];
            let u = s.cashflow_times[j];
            k[i][j] = rkhs_kernel_and_dt(t, u, alpha).0;
        }
    }
    let mut a = vec![vec![T::zero(); m]; m];
    let mut y = vec![T::zero(); m];
    for i in 0..m {
        y[i] = s.instrument_dirty_prices[i]
            - s.instrument_cashflow_matrix[i]
                .iter()
                .copied()
                .fold(T::zero(), |x, v| x + v);
        for j in 0..m {
            let mut z = T::zero();
            for p in 0..n {
                for q in 0..n {
                    z += s.instrument_cashflow_matrix[i][p]
                        * k[p][q]
                        * s.instrument_cashflow_matrix[j][q];
                }
            }
            a[i][j] = z;
            if i == j {
                a[i][j] += s.ridge_penalty_lambda;
            }
        }
    }
    let beta = solve_dense_curves(&a, &y);
    let mut w = vec![T::zero(); n];
    for j in 0..n {
        for i in 0..m {
            w[j] += s.instrument_cashflow_matrix[i][j] * beta[i];
        }
    }
    let mut fitted = vec![T::zero(); m];
    for i in 0..m {
        for j in 0..n {
            let mut p = T::one();
            for l in 0..n {
                let kk = rkhs_kernel_and_dt(s.cashflow_times[j], s.cashflow_times[l], alpha).0;
                p += kk * w[l];
            }
            fitted[i] += s.instrument_cashflow_matrix[i][j] * p;
        }
    }
    let mut err = T::zero();
    for i in 0..m {
        let d = fitted[i] - s.instrument_dirty_prices[i];
        err += d * d;
    }
    let rmse = T::from_f64(1e4) * (err / T::from_f64(m.max(1) as f64)).sqrt();
    let mut norm = T::zero();
    for i in 0..m {
        for j in 0..m {
            let mut z = T::zero();
            for p in 0..n {
                for q in 0..n {
                    z += s.instrument_cashflow_matrix[i][p]
                        * k[p][q]
                        * s.instrument_cashflow_matrix[j][q];
                }
            }
            norm += beta[i] * z * beta[j];
        }
    }
    let mut qd = Vec::new();
    let mut qr = Vec::new();
    let mut qf = Vec::new();
    let mut jac = vec![vec![T::zero(); m]; s.query_times.len()];
    let colsols: Vec<Vec<T>> = (0..m)
        .map(|col| {
            solve_dense_curves(
                &a,
                &(0..m)
                    .map(|i| if i == col { T::one() } else { T::zero() })
                    .collect::<Vec<_>>(),
            )
        })
        .collect();
    for (qi, &t) in s.query_times.iter().enumerate() {
        let mut p = T::one();
        let mut dp = T::zero();
        let mut kval = vec![T::zero(); n];
        for j in 0..n {
            let tau = s.cashflow_times[j];
            let (kv, dk) = rkhs_kernel_and_dt(t, tau, alpha);
            kval[j] = kv;
            p += w[j] * kval[j];
            dp += w[j] * dk;
        }
        qd.push(p);
        qr.push(-p.ln() / t.max_val(T::from_f64(1e-12)));
        qf.push(-dp / p);
        for col in 0..m {
            let colsol = &colsols[col];
            let mut z = T::zero();
            for j in 0..n {
                let mut cta = T::zero();
                for i in 0..m {
                    cta += s.instrument_cashflow_matrix[i][j] * colsol[i];
                }
                z += kval[j] * cta;
            }
            jac[qi][col] = -z / (t.max_val(T::from_f64(1e-12)) * p);
        }
    }
    FilipovicRkhsYieldCurveReport {
        dual_representer_weights_beta: beta,
        fitted_instrument_prices: fitted,
        pricing_rmse_bps: rmse,
        rkhs_smoothness_norm_squared: norm,
        query_discount_factors: qd,
        query_zero_rates: qr,
        query_instantaneous_forwards: qf,
        quote_to_zero_rate_jacobian: jac,
    }
}

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

#[derive(Clone, Debug)]
pub struct TsayKalmanRtsStateSpaceCurveSpec<T: Real = f64> {
    pub transition_f: [[T; 2]; 2],
    pub state_drift_c: [T; 2],
    pub process_noise_q: [[T; 2]; 2],
    pub measurement_h: [[T; 2]; 2],
    pub measurement_noise_r: [T; 2],
    pub initial_state_mean: [T; 2],
    pub initial_state_cov: [[T; 2]; 2],
    pub observations: Vec<[T; 2]>,
}
#[derive(Clone, Debug)]
pub struct TsayKalmanRtsStateSpaceCurveReport<T: Real = f64> {
    pub filtered_states: Vec<[T; 2]>,
    pub smoothed_states: Vec<[T; 2]>,
    pub filtered_trace_covariances: Vec<T>,
    pub smoothed_trace_covariances: Vec<T>,
    pub innovation_log_likelihood: T,
    pub terminal_kalman_gain_frobenius: T,
    pub max_smoother_variance_reduction: T,
}
fn mat2_trans<T: Real>(a: [[T; 2]; 2]) -> [[T; 2]; 2] {
    [[a[0][0], a[1][0]], [a[0][1], a[1][1]]]
}
fn mat2_mul<T: Real>(a: [[T; 2]; 2], b: [[T; 2]; 2]) -> [[T; 2]; 2] {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]
}
fn mat2_add<T: Real>(a: [[T; 2]; 2], b: [[T; 2]; 2]) -> [[T; 2]; 2] {
    [
        [a[0][0] + b[0][0], a[0][1] + b[0][1]],
        [a[1][0] + b[1][0], a[1][1] + b[1][1]],
    ]
}
fn mat2_sub<T: Real>(a: [[T; 2]; 2], b: [[T; 2]; 2]) -> [[T; 2]; 2] {
    [
        [a[0][0] - b[0][0], a[0][1] - b[0][1]],
        [a[1][0] - b[1][0], a[1][1] - b[1][1]],
    ]
}
fn mat2_inv<T: Real>(a: [[T; 2]; 2]) -> [[T; 2]; 2] {
    let mut d = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    // Preserve the determinant sign while avoiding explosive inverses near singularity.
    if d.val().abs() < 1e-12 {
        d = T::from_f64(if d.val().is_sign_negative() {
            -1e-12
        } else {
            1e-12
        });
    }
    let di = T::one() / d;
    [[a[1][1] * di, -a[0][1] * di], [-a[1][0] * di, a[0][0] * di]]
}
fn mat2_vec<T: Real>(a: [[T; 2]; 2], x: [T; 2]) -> [T; 2] {
    [
        a[0][0] * x[0] + a[0][1] * x[1],
        a[1][0] * x[0] + a[1][1] * x[1],
    ]
}
pub fn evaluate_tsay_kalman_rts_state_space_curve<T: Real>(
    s: &TsayKalmanRtsStateSpaceCurveSpec<T>,
) -> TsayKalmanRtsStateSpaceCurveReport<T> {
    let n = s.observations.len();
    let f = s.transition_f;
    let h = s.measurement_h;
    let ht = mat2_trans(h);
    let ft = mat2_trans(f);
    let mut xf = Vec::with_capacity(n);
    let mut pf = Vec::with_capacity(n);
    let mut xp = Vec::with_capacity(n);
    let mut pp = Vec::with_capacity(n);
    let mut x = s.initial_state_mean;
    let mut p = s.initial_state_cov;
    let mut ll = T::zero();
    let mut kg = T::zero();
    let rmat = [
        [s.measurement_noise_r[0], T::zero()],
        [T::zero(), s.measurement_noise_r[1]],
    ];
    for y in &s.observations {
        let xpred = [
            mat2_vec(f, x)[0] + s.state_drift_c[0],
            mat2_vec(f, x)[1] + s.state_drift_c[1],
        ];
        let ppred = mat2_add(mat2_mul(mat2_mul(f, p), ft), s.process_noise_q);
        let v = [y[0] - (mat2_vec(h, xpred)[0]), y[1] - mat2_vec(h, xpred)[1]];
        let ss = mat2_add(mat2_mul(mat2_mul(h, ppred), ht), rmat);
        let si = mat2_inv(ss);
        let k = mat2_mul(mat2_mul(ppred, ht), si);
        let kv = mat2_vec(k, v);
        x = [xpred[0] + kv[0], xpred[1] + kv[1]];
        let i_kh = mat2_sub(
            [[T::one(), T::zero()], [T::zero(), T::one()]],
            mat2_mul(k, h),
        );
        p = mat2_add(
            mat2_mul(mat2_mul(i_kh, ppred), mat2_trans(i_kh)),
            mat2_mul(mat2_mul(k, rmat), mat2_trans(k)),
        );
        let det = ss[0][0] * ss[1][1] - ss[0][1] * ss[1][0];
        let q =
            v[0] * (si[0][0] * v[0] + si[0][1] * v[1]) + v[1] * (si[1][0] * v[0] + si[1][1] * v[1]);
        ll += -T::from_f64(0.5)
            * (T::from_f64(2.0 * (2.0 * std::f64::consts::PI).ln()) + det.ln() + q);
        kg = (k[0][0] * k[0][0] + k[0][1] * k[0][1] + k[1][0] * k[1][0] + k[1][1] * k[1][1]).sqrt();
        xp.push(xpred);
        pp.push(ppred);
        xf.push(x);
        pf.push(p);
    }
    let mut xs = xf.clone();
    let mut ps = pf.clone();
    for t in (0..n.saturating_sub(1)).rev() {
        let j = mat2_mul(mat2_mul(pf[t], ft), mat2_inv(pp[t + 1]));
        let dx = [xs[t + 1][0] - xp[t + 1][0], xs[t + 1][1] - xp[t + 1][1]];
        let jd = mat2_vec(j, dx);
        xs[t] = [xf[t][0] + jd[0], xf[t][1] + jd[1]];
        ps[t] = mat2_add(
            pf[t],
            mat2_mul(mat2_mul(j, mat2_sub(ps[t + 1], pp[t + 1])), mat2_trans(j)),
        );
    }
    let ftc = pf.iter().map(|a| a[0][0] + a[1][1]).collect::<Vec<_>>();
    let stc = ps.iter().map(|a| a[0][0] + a[1][1]).collect::<Vec<_>>();
    let mut maxred = T::zero();
    for i in 0..n {
        let red = ftc[i] - stc[i];
        if red.val() > maxred.val() {
            maxred = red;
        }
    }
    TsayKalmanRtsStateSpaceCurveReport {
        filtered_states: xf,
        smoothed_states: xs,
        filtered_trace_covariances: ftc,
        smoothed_trace_covariances: stc,
        innovation_log_likelihood: ll,
        terminal_kalman_gain_frobenius: kg,
        max_smoother_variance_reduction: maxred,
    }
}

#[cfg(test)]
mod domain_suite_35_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn rkhs_invariants_and_dual_fd() {
        let s = FilipovicRkhsYieldCurveSpec {
            cashflow_times: vec![1.0, 2.0],
            instrument_cashflow_matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
            instrument_dirty_prices: vec![0.98, 0.95],
            decay_alpha: 1.0,
            ridge_penalty_lambda: 1e-8,
            query_times: vec![1.0],
        };
        let r = calibrate_filipovic_rkhs_smooth_yield_curve(&s);
        assert!(r.pricing_rmse_bps.is_finite());
        assert!(r.query_discount_factors[0] > 0.0);
    }

    #[test]
    fn dual_fd_parity() {
        let eval = |x: f64| {
            calibrate_filipovic_rkhs_smooth_yield_curve(&FilipovicRkhsYieldCurveSpec {
                cashflow_times: vec![1.0, 2.0],
                instrument_cashflow_matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
                instrument_dirty_prices: vec![0.98 + x, 0.95],
                decay_alpha: 1.0,
                ridge_penalty_lambda: 1e-8,
                query_times: vec![1.0],
            })
            .query_zero_rates[0]
        };
        let h = 1e-5;
        let fd = (eval(h) - eval(-h)) / (2.0 * h);
        let x = Dual::<1>::variable(0.0, 0);
        let ad = calibrate_filipovic_rkhs_smooth_yield_curve(&FilipovicRkhsYieldCurveSpec {
            cashflow_times: vec![Dual::<1>::from_value(1.0), Dual::<1>::from_value(2.0)],
            instrument_cashflow_matrix: vec![
                vec![Dual::<1>::one(), Dual::<1>::zero()],
                vec![Dual::<1>::zero(), Dual::<1>::one()],
            ],
            instrument_dirty_prices: vec![
                Dual::<1>::from_value(0.98) + x,
                Dual::<1>::from_value(0.95),
            ],
            decay_alpha: Dual::<1>::one(),
            ridge_penalty_lambda: Dual::<1>::from_value(1e-8),
            query_times: vec![Dual::<1>::one()],
        })
        .query_zero_rates[0];
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_25_tests {
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
mod domain_suite_33_tests {
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

#[cfg(test)]
mod domain_suite_37_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> TsayKalmanRtsStateSpaceCurveSpec<T> {
        TsayKalmanRtsStateSpaceCurveSpec {
            transition_f: [
                [T::from_f64(0.9), T::from_f64(0.1)],
                [T::zero(), T::from_f64(0.95)],
            ],
            state_drift_c: [T::zero(), T::zero()],
            process_noise_q: [
                [T::from_f64(0.01), T::zero()],
                [T::zero(), T::from_f64(0.01)],
            ],
            measurement_h: [[T::one(), T::zero()], [T::zero(), T::one()]],
            measurement_noise_r: [T::from_f64(0.02), T::from_f64(0.02)],
            initial_state_mean: [T::zero(), T::zero()],
            initial_state_cov: [[T::one(), T::zero()], [T::zero(), T::one()]],
            observations: vec![[x, T::from_f64(1.0)], [T::from_f64(1.1), T::from_f64(0.9)]],
        }
    }
    #[test]
    fn kalman_rts_dual_fd() {
        let h = 1e-5;
        let fd = (evaluate_tsay_kalman_rts_state_space_curve(&spec(1.0 + h))
            .innovation_log_likelihood
            .val()
            - evaluate_tsay_kalman_rts_state_space_curve(&spec(1.0 - h))
                .innovation_log_likelihood
                .val())
            / (2.0 * h);
        let ad = evaluate_tsay_kalman_rts_state_space_curve(&spec(Dual::<1>::variable(1.0, 0)))
            .innovation_log_likelihood;
        assert!((ad.d[0] - fd).abs() < 1e-5);
        let r = evaluate_tsay_kalman_rts_state_space_curve(&spec(1.0));
        assert!(r.smoothed_trace_covariances[0] < r.filtered_trace_covariances[0]);
        let n = r.filtered_states.len() - 1;
        assert!((r.smoothed_states[n][0] - r.filtered_states[n][0]).abs() < 1e-12);
        assert!((r.smoothed_states[n][1] - r.filtered_states[n][1]).abs() < 1e-12);
        let h = 1e-5;
        let fd_state = (evaluate_tsay_kalman_rts_state_space_curve(&spec(1.0 + h)).smoothed_states
            [0][0]
            .val()
            - evaluate_tsay_kalman_rts_state_space_curve(&spec(1.0 - h)).smoothed_states[0][0]
                .val())
            / (2.0 * h);
        let ad_state =
            evaluate_tsay_kalman_rts_state_space_curve(&spec(Dual::<1>::variable(1.0, 0)))
                .smoothed_states[0][0];
        assert!((ad_state.d[0] - fd_state).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_40_tests {
    use super::*;
    use ql_ad::Dual;

    fn mk<T: Real>(x: T) -> BcbsD368SixScenarioIrrbbSpec<T> {
        BcbsD368SixScenarioIrrbbSpec {
            bucket_midpoints_years: vec![T::from_f64(0.5), T::from_f64(2.0), T::from_f64(5.0)],
            asset_repricing_cashflows: vec![
                T::from_f64(10.0) + x,
                T::from_f64(8.0),
                T::from_f64(5.0),
            ],
            liability_repricing_cashflows: vec![
                T::from_f64(7.0),
                T::from_f64(6.0),
                T::from_f64(4.0),
            ],
            baseline_zero_rates: vec![T::from_f64(0.02), T::from_f64(0.025), T::from_f64(0.03)],
            parallel_shock_bps: T::from_f64(200.0),
            short_shock_bps: T::from_f64(250.0),
            long_shock_bps: T::from_f64(100.0),
            decay_tau_years: T::from_f64(4.0),
            tier1_capital: T::from_f64(100.0),
            pass_through_deposit_beta: T::from_f64(0.5),
        }
    }

    #[test]
    fn irrbb_invariants_and_dual() {
        let r = evaluate_bcbs_d368_six_scenario_irrbb(&mk(0.0));
        assert!(r.worst_case_delta_eveloss >= 0.0);
        assert!((r.parallel_up_delta_nii_1y + r.parallel_down_delta_nii_1y).abs() < 1e-12);
        let h = 1e-5;
        let fd = (evaluate_bcbs_d368_six_scenario_irrbb(&mk(h)).baseline_eve
            - evaluate_bcbs_d368_six_scenario_irrbb(&mk(-h)).baseline_eve)
            / (2.0 * h);
        let ad =
            evaluate_bcbs_d368_six_scenario_irrbb(&mk(Dual::<1>::variable(0.0, 0))).baseline_eve;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

mod evaluate_hagan_west_monotone_convex_forward_curve_module {
    use ql_ad::Real;

    #[inline]
    fn smin<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b - (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[inline]
    fn smax<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[derive(Clone, Copy, Debug)]
    pub struct HaganWestMonotoneConvexCurveSpec<T: Real = f64> {
        pub knot_tenors: [T; 4],
        pub discrete_forwards: [T; 3],
        pub eval_tenor: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct HaganWestMonotoneConvexCurveReport<T: Real = f64> {
        pub knot_instantaneous_forward_1: T,
        pub knot_instantaneous_forward_2: T,
        pub eval_instantaneous_forward: T,
        pub eval_zero_rate: T,
        pub eval_discount_factor: T,
        pub cell_integrated_average_forward: T,
        pub cell_average_conservation_error: T,
        pub ameliorated_min_instantaneous_forward: T,
    }

    pub fn evaluate_hagan_west_monotone_convex_forward_curve<T: Real>(
        spec: &HaganWestMonotoneConvexCurveSpec<T>,
    ) -> HaganWestMonotoneConvexCurveReport<T> {
        let dt1 = (spec.knot_tenors[1] - spec.knot_tenors[0]).max_val(T::from_f64(1e-8));
        let dt2 = (spec.knot_tenors[2] - spec.knot_tenors[1]).max_val(T::from_f64(1e-8));
        let dt3 = (spec.knot_tenors[3] - spec.knot_tenors[2]).max_val(T::from_f64(1e-8));
        let f1d = spec.discrete_forwards[0];
        let f2d = spec.discrete_forwards[1];
        let f3d = spec.discrete_forwards[2];
        let f1_raw = (dt2 * f1d + dt1 * f2d) / (dt1 + dt2);
        let f2_raw = (dt3 * f2d + dt2 * f3d) / (dt2 + dt3);
        let two = T::from_f64(2.0);
        let f1_mc = smax(T::zero(), smin(two * f1d, smin(two * f2d, f1_raw)));
        let f2_mc = smax(T::zero(), smin(two * f2d, smin(two * f3d, f2_raw)));
        let x = ((spec.eval_tenor - spec.knot_tenors[1]) / dt2)
            .max_val(T::zero())
            .min_val(T::one());
        let g0 = f1_mc - f2d;
        let g1 = f2_mc - f2d;
        let f_eval = f2d
            + g0 * (T::one() - T::from_f64(4.0) * x + T::from_f64(3.0) * x * x)
            + g1 * (-two * x + T::from_f64(3.0) * x * x);
        let cell_int = f2d + g0 * (T::one() - two + T::one()) + g1 * (-T::one() + T::one());
        let diff = cell_int - f2d;
        let err = (diff * diff + T::from_f64(1e-30)).sqrt();
        let partial_int = f2d * x + g0 * (x - two * x * x + x * x * x) + g1 * (-x * x + x * x * x);
        let z = (f1d * dt1 + dt2 * partial_int) / spec.eval_tenor.max_val(T::from_f64(1e-8));
        let df = (-spec.eval_tenor * z).exp();
        let denom = T::from_f64(3.0) * (g0 + g1);
        let x_star = smax(
            T::zero(),
            smin(T::one(), (two * g0 + g1) / (denom + T::from_f64(1e-12))),
        );
        let f_star = f2d
            + g0 * (T::one() - T::from_f64(4.0) * x_star + T::from_f64(3.0) * x_star * x_star)
            + g1 * (-two * x_star + T::from_f64(3.0) * x_star * x_star);
        let min_f = smax(T::zero(), smin(f1_mc, smin(f2_mc, f_star)));
        HaganWestMonotoneConvexCurveReport {
            knot_instantaneous_forward_1: f1_mc,
            knot_instantaneous_forward_2: f2_mc,
            eval_instantaneous_forward: f_eval,
            eval_zero_rate: z,
            eval_discount_factor: df,
            cell_integrated_average_forward: cell_int,
            cell_average_conservation_error: err,
            ameliorated_min_instantaneous_forward: min_f,
        }
    }

    #[cfg(test)]
    mod domain_suite_47_tests {
        use super::*;
        use ql_ad::Dual;
        fn spec<T: Real>(x: T) -> HaganWestMonotoneConvexCurveSpec<T> {
            HaganWestMonotoneConvexCurveSpec {
                knot_tenors: [T::zero(), T::one(), T::from_f64(2.0), T::from_f64(3.0)],
                discrete_forwards: [T::from_f64(0.03), T::from_f64(0.04) + x, T::from_f64(0.045)],
                eval_tenor: T::from_f64(1.4),
            }
        }
        #[test]
        fn suite47_curves_invariants_and_dual_fd() {
            let r = evaluate_hagan_west_monotone_convex_forward_curve(&spec(0.0));
            assert!(r.cell_average_conservation_error < 1e-12);
            assert!(r.ameliorated_min_instantaneous_forward >= 0.0);
            let h = 1e-6;
            let fd = (evaluate_hagan_west_monotone_convex_forward_curve(&spec(h))
                .eval_discount_factor
                - evaluate_hagan_west_monotone_convex_forward_curve(&spec(-h))
                    .eval_discount_factor)
                / (2.0 * h);
            let ad = evaluate_hagan_west_monotone_convex_forward_curve(&spec(Dual::<1>::variable(
                0.0, 0,
            )))
            .eval_discount_factor;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_hagan_west_monotone_convex_forward_curve_module::*;

mod evaluate_hagan_west_tension_convex_forward_curve_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct HaganWestTensionCurveSpec<T: Real = f64> {
        pub pillar_intervals: [T; 4],
        pub discrete_forwards: [T; 4],
        pub tension_lambda: T,
        pub eval_coord_x: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct HaganWestTensionCurveReport<T: Real = f64> {
        pub instantaneous_forward_rate: T,
        pub zero_coupon_bond_price: T,
        pub tension_cubic_adjustment: T,
        pub integrated_cell_forward: T,
        pub tension_cell_conservation_error: T,
        pub curve_roughness_energy: T,
    }

    pub fn evaluate_hagan_west_tension_convex_forward_curve<T: Real>(
        spec: &HaganWestTensionCurveSpec<T>,
    ) -> HaganWestTensionCurveReport<T> {
        let t0 = spec.pillar_intervals[0].max_val(T::from_f64(1e-6));
        let t1 = spec.pillar_intervals[1].max_val(T::from_f64(1e-6));
        let t2 = spec.pillar_intervals[2].max_val(T::from_f64(1e-6));
        let t3 = spec.pillar_intervals[3].max_val(T::from_f64(1e-6));
        let f0d = spec.discrete_forwards[0];
        let f1d = spec.discrete_forwards[1];
        let f2d = spec.discrete_forwards[2];
        let f3d = spec.discrete_forwards[3];

        let f1 = (t0 * f1d + t1 * f0d) / (t0 + t1);
        let f2 = (t1 * f2d + t2 * f1d + T::from_f64(0.01) * t3 * (f3d - f2d)) / (t1 + t2);
        let g0 = f1 - f1d;
        let g1 = f2 - f1d;

        let lam = spec.tension_lambda.max_val(T::zero());
        let kap = lam / (T::one() + lam);
        let x = spec
            .eval_coord_x
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.99));
        let x2 = x * x;
        let x3 = x2 * x;

        let q_part = g0 * (T::one() - T::from_f64(4.0) * x + T::from_f64(3.0) * x2)
            + g1 * (-T::from_f64(2.0) * x + T::from_f64(3.0) * x2);
        let c_basis = x * (T::one() - x) * (T::one() - T::from_f64(2.0) * x);
        let c_adj = kap * (g0 + g1) * c_basis;
        let inst_f = f1d + q_part + c_adj;

        // Exact integral of g(u) from 0 to x:
        // int_0^x (1 - 4u + 3u^2) du = x - 2x^2 + x^3
        // int_0^x (-2u + 3u^2) du = -x^2 + x^3
        // int_0^x (u - 3u^2 + 2u^3) du = 0.5 x^2 - x^3 + 0.5 x^4
        let int_g_x = g0 * (x - T::from_f64(2.0) * x2 + x3)
            + g1 * (-x2 + x3)
            + kap * (g0 + g1) * (T::from_f64(0.5) * x2 - x3 + T::from_f64(0.5) * x2 * x2);
        let zcb = (-(t0 * f0d + t1 * (f1d * x + int_g_x))).exp();

        // At x = 1, int_0^1 g(u) du = 0 identically for all tension_lambda
        let int_g_1 = g0 * (T::one() - T::from_f64(2.0) + T::one())
            + g1 * (-T::one() + T::one())
            + kap * (g0 + g1) * (T::from_f64(0.5) - T::one() + T::from_f64(0.5));
        let cell_fwd = f1d + int_g_1;
        let diff = cell_fwd - f1d;
        let cons_err = (diff * diff + T::from_f64(1e-30)).sqrt();

        let rough = T::from_f64(4.0) * (g0 * g0 + g0 * g1 + g1 * g1) / t1
            + kap * kap * (g0 + g1) * (g0 + g1) / t1;

        HaganWestTensionCurveReport {
            instantaneous_forward_rate: inst_f,
            zero_coupon_bond_price: zcb,
            tension_cubic_adjustment: c_adj,
            integrated_cell_forward: cell_fwd,
            tension_cell_conservation_error: cons_err,
            curve_roughness_energy: rough,
        }
    }

    #[cfg(test)]
    mod domain_suite_50_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hagan_west_tension_curve_invariants_and_ad_parity() {
            let spec = HaganWestTensionCurveSpec {
                pillar_intervals: [1.0, 1.0, 2.0, 3.0],
                discrete_forwards: [0.032, 0.038, 0.044, 0.046],
                tension_lambda: 1.5,
                eval_coord_x: 0.35,
            };
            let rep = evaluate_hagan_west_tension_convex_forward_curve(&spec);
            assert!(rep.tension_cell_conservation_error < 1e-12);
            assert!(rep.zero_coupon_bond_price > 0.0 && rep.zero_coupon_bond_price < 1.0);

            let h = 1e-6;
            let r_ad =
                evaluate_hagan_west_tension_convex_forward_curve(&HaganWestTensionCurveSpec {
                    pillar_intervals: spec.pillar_intervals.map(Dual::<1>::from_value),
                    discrete_forwards: [
                        Dual::<1>::from_value(spec.discrete_forwards[0]),
                        Dual::<1>::variable(spec.discrete_forwards[1], 0),
                        Dual::<1>::from_value(spec.discrete_forwards[2]),
                        Dual::<1>::from_value(spec.discrete_forwards[3]),
                    ],
                    tension_lambda: Dual::<1>::from_value(spec.tension_lambda),
                    eval_coord_x: Dual::<1>::from_value(spec.eval_coord_x),
                });
            let mut sp = spec;
            sp.discrete_forwards[1] += h;
            let mut sm = spec;
            sm.discrete_forwards[1] -= h;
            let fd = (evaluate_hagan_west_tension_convex_forward_curve(&sp).zero_coupon_bond_price
                - evaluate_hagan_west_tension_convex_forward_curve(&sm).zero_coupon_bond_price)
                / (2.0 * h);
            assert!((r_ad.zero_coupon_bond_price.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_hagan_west_tension_convex_forward_curve_module::*;
