use crate::curves::{DiscountCurve, YieldTermStructure};
use nalgebra::DMatrix;
use ql_ad::{Dual, IftSystemSensitivity, Real, TapeScope, TapeVar};
use ql_time::{Date, DayCounter, Frequency, Period};
use thiserror::Error;

#[derive(Debug, Error, Clone)]
pub enum BootstrapError {
    #[error("invalid bootstrap input")]
    InvalidInput,
    #[error("solver error: {0}")]
    Solver(String),
}
pub trait CalibrationHelper<T: Real = f64>: Send + Sync {
    fn residual(&self, curve: &impl YieldTermStructure<T>) -> T;
    fn quote(&self) -> f64;
}
#[derive(Debug, Clone)]
pub struct DepositInstrument {
    pub start: Date,
    pub end: Date,
    pub quote: f64,
    pub day_counter: DayCounter,
}
impl DepositInstrument {
    pub fn new(start: Date, end: Date, quote: f64, day_counter: DayCounter) -> Self {
        Self {
            start,
            end,
            quote,
            day_counter,
        }
    }
    pub fn par_rate<T: Real>(&self, c: &impl YieldTermStructure<T>) -> T {
        (c.discount_factor((self.start - c.reference_date()) as f64 / 365.0)
            / c.discount_factor((self.end - c.reference_date()) as f64 / 365.0)
            - T::one())
            / self
                .day_counter
                .year_fraction(self.start, self.end, None, None)
    }
}
impl<T: Real> CalibrationHelper<T> for DepositInstrument {
    fn residual(&self, c: &impl YieldTermStructure<T>) -> T {
        self.par_rate(c) - self.quote
    }
    fn quote(&self) -> f64 {
        self.quote
    }
}
#[derive(Debug, Clone)]
pub struct FraInstrument {
    pub fixing_start: Date,
    pub fixing_end: Date,
    pub quote: f64,
    pub day_counter: DayCounter,
}
impl<T: Real> CalibrationHelper<T> for FraInstrument {
    fn residual(&self, c: &impl YieldTermStructure<T>) -> T {
        let t1 = (self.fixing_start - c.reference_date()) as f64 / 365.0;
        let t2 = (self.fixing_end - c.reference_date()) as f64 / 365.0;
        (c.discount_factor(t1) / c.discount_factor(t2) - T::one())
            / self
                .day_counter
                .year_fraction(self.fixing_start, self.fixing_end, None, None)
            - self.quote
    }
    fn quote(&self) -> f64 {
        self.quote
    }
}
#[derive(Debug, Clone)]
pub struct OisSwapInstrument {
    pub start: Date,
    pub end: Date,
    pub fixed_rate: f64,
    pub frequency: Frequency,
    pub day_counter: DayCounter,
}
impl<T: Real> CalibrationHelper<T> for OisSwapInstrument {
    fn residual(&self, c: &impl YieldTermStructure<T>) -> T {
        let mut ann = T::zero();
        let mut prev = self.start;
        let period = Period::from(self.frequency);
        while prev < self.end {
            let next = prev.add_period(period, false).min(self.end);
            let tau = self.day_counter.year_fraction(prev, next, None, None);
            let t = (next - c.reference_date()) as f64 / 365.0;
            ann += c.discount_factor(t) * tau;
            prev = next;
        }
        let p0 = c.discount_factor((self.start - c.reference_date()) as f64 / 365.0);
        let pe = c.discount_factor((self.end - c.reference_date()) as f64 / 365.0);
        (p0 - pe) / ann - self.fixed_rate
    }
    fn quote(&self) -> f64 {
        self.fixed_rate
    }
}
#[derive(Debug, Clone)]
pub struct IrsSwapInstrument {
    pub start: Date,
    pub end: Date,
    pub fixed_rate: f64,
    pub frequency: Frequency,
    pub day_counter: DayCounter,
}
impl<T: Real> CalibrationHelper<T> for IrsSwapInstrument {
    fn residual(&self, c: &impl YieldTermStructure<T>) -> T {
        OisSwapInstrument {
            start: self.start,
            end: self.end,
            fixed_rate: self.fixed_rate,
            frequency: self.frequency,
            day_counter: self.day_counter.clone(),
        }
        .residual(c)
    }
    fn quote(&self) -> f64 {
        self.fixed_rate
    }
}
#[derive(Debug, Clone)]
pub struct BasisSwapInstrument {
    pub start: Date,
    pub end: Date,
    pub spread: f64,
    pub day_counter: DayCounter,
}
impl<T: Real> CalibrationHelper<T> for BasisSwapInstrument {
    fn residual(&self, c: &impl YieldTermStructure<T>) -> T {
        let t1 = (self.start - c.reference_date()) as f64 / 365.0;
        let t2 = (self.end - c.reference_date()) as f64 / 365.0;
        (c.discount_factor(t1) / c.discount_factor(t2) - T::one())
            / self
                .day_counter
                .year_fraction(self.start, self.end, None, None)
            - self.spread
    }
    fn quote(&self) -> f64 {
        self.spread
    }
}

pub fn bootstrap_piecewise_curve(
    reference_date: Date,
    times: &[f64],
    quotes: &[f64],
) -> Result<DiscountCurve<f64>, BootstrapError> {
    if times.len() != quotes.len()
        || times.is_empty()
        || times.iter().any(|t| !t.is_finite() || *t <= 0.0)
        || times.windows(2).any(|w| w[1] <= w[0])
        || quotes.iter().any(|q| !q.is_finite())
    {
        return Err(BootstrapError::InvalidInput);
    }
    let mut ds = Vec::with_capacity(times.len());
    for (t, q) in times.iter().zip(quotes) {
        ds.push(1.0 / (1.0 + q * t));
    }
    Ok(DiscountCurve::new(reference_date, times.to_vec(), ds))
}
pub fn bootstrap_global_multicurve(
    reference_date: Date,
    times: &[f64],
    quotes: &[f64],
) -> Result<DiscountCurve<f64>, BootstrapError> {
    bootstrap_piecewise_curve(reference_date, times, quotes)
}
#[derive(Debug, Clone)]
pub struct CalibratedCurveSystem {
    pub parameters: Vec<f64>,
    pub quotes: Vec<f64>,
    pub curves: Vec<DiscountCurve<f64>>,
    pub sensitivity: IftSystemSensitivity,
}
impl CalibratedCurveSystem {
    pub fn from_single_curve_helpers(
        reference_date: Date,
        helpers: &[&dyn SingleCurveHelper],
    ) -> Result<Self, BootstrapError> {
        let curve = bootstrap_piecewise_curve_helpers(reference_date, helpers)?;
        let theta = curve.discounts.clone();
        let k = theta.len();
        let mut jt = DMatrix::zeros(k, k);
        for j in 0..k {
            let nodes = theta
                .iter()
                .enumerate()
                .map(|(i, &x)| {
                    if i == j {
                        Dual::<1>::variable(x, 0)
                    } else {
                        Dual::from_value(x)
                    }
                })
                .collect::<Vec<_>>();
            let cd = DiscountCurve::new(reference_date, curve.times.clone(), nodes);
            for i in 0..k {
                jt[(i, j)] = helpers[i].residual_curve_dual(&cd).d[0];
            }
        }
        let jq = -DMatrix::identity(k, k);
        Ok(Self::new(
            theta.clone(),
            helpers.iter().map(|h| h.quote_value()).collect(),
            jt,
            jq,
            vec![curve],
        ))
    }
    pub fn from_multicurve_helpers(
        reference_date: Date,
        ois_helpers: &[&dyn SingleCurveHelper],
        fwd_helpers: &[IrsSwapInstrument],
    ) -> Result<Self, BootstrapError> {
        let ois = bootstrap_piecewise_curve_helpers(reference_date, ois_helpers)?;
        let times: Vec<f64> = fwd_helpers
            .iter()
            .map(|h| (h.end - reference_date) as f64 / 365.0)
            .collect();
        let mut dfs = vec![1.0; times.len()];
        for k in 0..fwd_helpers.len() {
            let h = &fwd_helpers[k];
            let prev = dfs.clone();
            let f = |d: f64| {
                let mut t = prev.clone();
                t[k] = d;
                let c = DiscountCurve::new(reference_date, times.clone(), t);
                h.residual_with_discount(&c, &ois).val()
            };
            dfs[k] = ql_math::brent_solve(f, 1e-6, 2.0, 1e-14, 200)
                .map_err(|e| BootstrapError::Solver(e.to_string()))?;
        }
        let fwd = DiscountCurve::new(reference_date, times, dfs);
        let mut params = ois.discounts.clone();
        params.extend(fwd.discounts.clone());
        let ko = ois.discounts.len();
        let k = params.len();
        let mut jt = DMatrix::zeros(k, k);
        for j in 0..k {
            let nodes = params
                .iter()
                .enumerate()
                .map(|(i, &x)| {
                    if i == j {
                        Dual::<1>::variable(x, 0)
                    } else {
                        Dual::from_value(x)
                    }
                })
                .collect::<Vec<_>>();
            let oc = DiscountCurve::new(reference_date, ois.times.clone(), nodes[..ko].to_vec());
            let fc = DiscountCurve::new(reference_date, fwd.times.clone(), nodes[ko..].to_vec());
            for (i, helper) in ois_helpers.iter().enumerate() {
                jt[(i, j)] = helper.residual_curve_dual(&oc).d[0];
            }
            for (i, helper) in fwd_helpers.iter().enumerate() {
                jt[(ko + i, j)] = helper.residual_with_discount(&fc, &oc).d[0];
            }
        }
        let jq = -DMatrix::identity(k, k);
        let mut quotes = ois_helpers
            .iter()
            .map(|h| h.quote_value())
            .collect::<Vec<_>>();
        quotes.extend(fwd_helpers.iter().map(|h| h.fixed_rate));
        Ok(Self::new(params, quotes, jt, jq, vec![ois, fwd]))
    }
}
impl CalibratedCurveSystem {
    pub fn new(
        parameters: Vec<f64>,
        quotes: Vec<f64>,
        j_theta: DMatrix<f64>,
        j_q: DMatrix<f64>,
        curves: Vec<DiscountCurve<f64>>,
    ) -> Self {
        Self {
            parameters,
            quotes,
            curves,
            sensitivity: IftSystemSensitivity::new(j_theta, j_q),
        }
    }
    pub fn vjp(&self, p: &[f64]) -> Vec<f64> {
        self.sensitivity.vjp(p)
    }
    pub fn jvp(&self, q: &[f64]) -> Vec<f64> {
        self.sensitivity.jvp(q)
    }
    pub fn value_and_par_dv01<F>(&self, pricer: F) -> (f64, Vec<f64>)
    where
        F: Fn(&[TapeVar<f64>]) -> TapeVar<f64>,
    {
        TapeScope::start_recording();
        let vars: Vec<_> = self
            .parameters
            .iter()
            .map(|&x| TapeScope::with_active(|t| t.input(x)))
            .collect();
        let out = pricer(&vars);
        let adj = TapeScope::gradient_from(out);
        let pillar: Vec<f64> = vars.iter().map(|v| adj.gradient(*v)).collect();
        let value = out.val;
        TapeScope::stop_recording();
        (value, self.sensitivity.vjp(&pillar))
    }
    pub fn value_par_dv01_and_gamma<F>(&self, pricer: F) -> (f64, Vec<f64>, DMatrix<f64>)
    where
        F: Fn(&[TapeVar<f64>]) -> TapeVar<f64>,
    {
        let eval_grad = |params: &[f64]| {
            TapeScope::start_recording();
            let vars: Vec<_> = params
                .iter()
                .map(|&x| TapeScope::with_active(|t| t.input(x)))
                .collect::<Vec<_>>();
            let out = pricer(&vars);
            let adj = TapeScope::gradient_from(out);
            let g: Vec<f64> = vars.iter().map(|v| adj.gradient(*v)).collect();
            let value = out.val;
            TapeScope::stop_recording();
            (value, g)
        };
        let (value, g) = eval_grad(&self.parameters);
        let k = g.len();
        let mut h = DMatrix::zeros(k, k);
        let bump = 1e-5;
        for jcol in 0..k {
            let mut plus = self.parameters.clone();
            let mut minus = self.parameters.clone();
            plus[jcol] += bump;
            minus[jcol] -= bump;
            let (_, gp) = eval_grad(&plus);
            let (_, gm) = eval_grad(&minus);
            for i in 0..k {
                h[(i, jcol)] = (gp[i] - gm[i]) / (2.0 * bump);
            }
        }
        let j = self.sensitivity.jacobian_dtheta_dq();
        let gamma = j.transpose() * h * j;
        let dv = self.sensitivity.vjp(&g);
        (value, dv, gamma)
    }
}
pub trait SingleCurveHelper: Send + Sync {
    fn pillar_time(&self, reference_date: Date) -> f64;
    fn quote_value(&self) -> f64;
    fn residual_curve(&self, curve: &DiscountCurve<f64>) -> f64;
    fn residual_curve_dual(&self, curve: &DiscountCurve<Dual<1>>) -> Dual<1> {
        let values = curve.discounts.iter().map(|x| x.val).collect::<Vec<_>>();
        let base = DiscountCurve::new(curve.reference_date, curve.times.clone(), values.clone());
        let value = self.residual_curve(&base);
        let Some(j) = curve.discounts.iter().position(|x| x.d[0].abs() > 0.5) else {
            return Dual::from_value(value);
        };
        let h = 1e-6;
        let mut plus = values.clone();
        let mut minus = values;
        plus[j] += h;
        minus[j] -= h;
        let rp = self.residual_curve(&DiscountCurve::new(
            curve.reference_date,
            curve.times.clone(),
            plus,
        ));
        let rm = self.residual_curve(&DiscountCurve::new(
            curve.reference_date,
            curve.times.clone(),
            minus,
        ));
        Dual {
            val: value,
            d: [(rp - rm) / (2.0 * h)],
        }
    }
}
impl SingleCurveHelper for DepositInstrument {
    fn pillar_time(&self, reference_date: Date) -> f64 {
        (self.end - reference_date) as f64 / 365.0
    }
    fn quote_value(&self) -> f64 {
        self.quote
    }
    fn residual_curve(&self, c: &DiscountCurve<f64>) -> f64 {
        self.par_rate(c).val() - self.quote
    }
    fn residual_curve_dual(&self, c: &DiscountCurve<Dual<1>>) -> Dual<1> {
        self.par_rate(c) - self.quote
    }
}
impl SingleCurveHelper for FraInstrument {
    fn pillar_time(&self, reference_date: Date) -> f64 {
        (self.fixing_end - reference_date) as f64 / 365.0
    }
    fn quote_value(&self) -> f64 {
        self.quote
    }
    fn residual_curve(&self, c: &DiscountCurve<f64>) -> f64 {
        let t1 = (self.fixing_start - c.reference_date()) as f64 / 365.0;
        let t2 = (self.fixing_end - c.reference_date()) as f64 / 365.0;
        (c.discount_factor(t1) / c.discount_factor(t2) - 1.0)
            / self
                .day_counter
                .year_fraction(self.fixing_start, self.fixing_end, None, None)
            - self.quote
    }
    fn residual_curve_dual(&self, c: &DiscountCurve<Dual<1>>) -> Dual<1> {
        let t1 = (self.fixing_start - c.reference_date()) as f64 / 365.0;
        let t2 = (self.fixing_end - c.reference_date()) as f64 / 365.0;
        (c.discount_factor(t1) / c.discount_factor(t2) - Dual::one())
            / self
                .day_counter
                .year_fraction(self.fixing_start, self.fixing_end, None, None)
            - self.quote
    }
}
pub fn bootstrap_piecewise_curve_helpers(
    reference_date: Date,
    helpers: &[&dyn SingleCurveHelper],
) -> Result<DiscountCurve<f64>, BootstrapError> {
    if helpers.is_empty() {
        return Err(BootstrapError::InvalidInput);
    }
    let times: Vec<f64> = helpers
        .iter()
        .map(|h| h.pillar_time(reference_date))
        .collect();
    if times.iter().any(|t| !t.is_finite() || *t <= 0.0)
        || times.windows(2).any(|w| w[1] <= w[0])
        || helpers.iter().any(|h| !h.quote_value().is_finite())
    {
        return Err(BootstrapError::InvalidInput);
    }
    let mut dfs = vec![1.0; times.len()];
    for k in 0..helpers.len() {
        let h = helpers[k];
        let prev = dfs.clone();
        let f = |d: f64| {
            let mut trial = prev.clone();
            trial[k] = d;
            let c = DiscountCurve::new(reference_date, times.clone(), trial);
            h.residual_curve(&c)
        };
        dfs[k] = ql_math::brent_solve(f, 1e-6, 2.0, 1e-14, 200)
            .map_err(|e| BootstrapError::Solver(e.to_string()))?;
    }
    Ok(DiscountCurve::new(reference_date, times, dfs))
}
impl IrsSwapInstrument {
    pub fn residual_with_discount<T: Real>(
        &self,
        projection: &impl YieldTermStructure<T>,
        discount: &impl YieldTermStructure<T>,
    ) -> T {
        let mut fixed = T::zero();
        let mut floating = T::zero();
        let mut prev = self.start;
        let period = Period::from(self.frequency);
        while prev < self.end {
            let next = prev.add_period(period, false).min(self.end);
            let tau = self.day_counter.year_fraction(prev, next, None, None);
            let t = (next - projection.reference_date()) as f64 / 365.0;
            let p = discount.discount_factor(t);
            fixed += p * tau;
            let a = (prev - projection.reference_date()) as f64 / 365.0;
            floating +=
                (projection.discount_factor(a) / projection.discount_factor(t) - T::one()) * p;
            prev = next;
        }
        floating / fixed - self.fixed_rate
    }
}
impl BasisSwapInstrument {
    pub fn residual_two_curves<T: Real>(
        &self,
        leg1: &impl YieldTermStructure<T>,
        leg2: &impl YieldTermStructure<T>,
        discount: &impl YieldTermStructure<T>,
    ) -> T {
        let n = ((self.end - self.start) as f64 / 90.0).round() as usize;
        let mut x = T::zero();
        for i in 1..=n {
            let t = (self.end - discount.reference_date()) as f64 / 365.0 - (n - i) as f64 * 0.25;
            let p = discount.discount_factor(t);
            let a = (t - 0.25).max(0.0);
            x += (leg1.discount_factor(a) / leg1.discount_factor(t)
                - leg2.discount_factor(a) / leg2.discount_factor(t))
                * p;
        }
        x - T::from_f64(self.spread)
    }
}
impl SingleCurveHelper for OisSwapInstrument {
    fn pillar_time(&self, reference_date: Date) -> f64 {
        (self.end - reference_date) as f64 / 365.0
    }
    fn quote_value(&self) -> f64 {
        self.fixed_rate
    }
    fn residual_curve(&self, c: &DiscountCurve<f64>) -> f64 {
        self.residual(c).val()
    }
    fn residual_curve_dual(&self, c: &DiscountCurve<Dual<1>>) -> Dual<1> {
        self.residual(c)
    }
}
