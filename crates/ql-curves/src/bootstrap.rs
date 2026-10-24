use crate::curves::{DiscountCurve, YieldTermStructure};
use nalgebra::{DMatrix, DVector};
use ql_ad::{Dual, Dual2, IftSystemSensitivity, Real, TapeScope, TapeVar};
use ql_time::{Date, DayCounter, Frequency, Period};
use std::any::Any;
use thiserror::Error;

#[inline]
fn curve_time(dc: &DayCounter, reference_date: Date, date: Date) -> f64 {
    dc.year_fraction(reference_date, date, None, None)
}

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
        (c.discount_factor(curve_time(
            &self.day_counter,
            c.reference_date(),
            self.start,
        )) / c.discount_factor(curve_time(&self.day_counter, c.reference_date(), self.end))
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
        let t1 = curve_time(&self.day_counter, c.reference_date(), self.fixing_start);
        let t2 = curve_time(&self.day_counter, c.reference_date(), self.fixing_end);
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
            let t = curve_time(&self.day_counter, c.reference_date(), next);
            ann += c.discount_factor(t) * tau;
            prev = next;
        }
        let p0 = c.discount_factor(curve_time(
            &self.day_counter,
            c.reference_date(),
            self.start,
        ));
        let pe = c.discount_factor(curve_time(&self.day_counter, c.reference_date(), self.end));
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
        let t1 = curve_time(&self.day_counter, c.reference_date(), self.start);
        let t2 = curve_time(&self.day_counter, c.reference_date(), self.end);
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
    for (i, (&t, &q)) in times.iter().zip(quotes).enumerate() {
        let d = if t <= 1.0 {
            1.0 / (1.0 + q * t)
        } else {
            // Sequential par-swap bootstrap on the supplied payment grid:
            // q * sum(delta_j P_j) + P_i = 1.  All prior discount factors
            // are already calibrated, leaving a one-dimensional linear
            // solve for the current pillar.
            let mut ann_prev = 0.0;
            let mut prev_t = 0.0;
            for j in 0..i {
                ann_prev += ds[j] * (times[j] - prev_t);
                prev_t = times[j];
            }
            let tau = t - prev_t;
            (1.0 - q * ann_prev) / (1.0 + q * tau)
        };
        if !d.is_finite() || d <= 0.0 {
            return Err(BootstrapError::Solver(
                "non-positive discount factor".into(),
            ));
        }
        ds.push(d);
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
    /// Calibration instruments retained so second-order quote curvature can
    /// include the implicit-function constraint Hessians.
    pub instruments: Vec<RateInstrumentSpec>,
    pub ois_pillars: usize,
}

#[derive(Debug, Clone)]
pub enum RateInstrumentKind {
    Deposit(DepositInstrument),
    Fra(FraInstrument),
    OisSwap(OisSwapInstrument),
    IrsSwap(IrsSwapInstrument),
    BasisSwap(BasisSwapInstrument),
}

#[derive(Debug, Clone)]
pub struct RateInstrumentSpec {
    pub kind: RateInstrumentKind,
    pub reference_date: Date,
    pub ois_times: Vec<f64>,
    pub fwd_times: Vec<f64>,
}

impl RateInstrumentSpec {
    fn residual<T: Real>(&self, parameters: &[T], ois_pillars: usize) -> T {
        let ois = DiscountCurve::new(
            self.reference_date,
            self.ois_times.clone(),
            parameters[..ois_pillars].to_vec(),
        );
        match &self.kind {
            RateInstrumentKind::Deposit(x) => x.residual(&ois),
            RateInstrumentKind::Fra(x) => x.residual(&ois),
            RateInstrumentKind::OisSwap(x) => x.residual(&ois),
            RateInstrumentKind::IrsSwap(x) => {
                let fwd = DiscountCurve::new(
                    self.reference_date,
                    self.fwd_times.clone(),
                    parameters[ois_pillars..].to_vec(),
                );
                x.residual_with_discount(&fwd, &ois)
            }
            RateInstrumentKind::BasisSwap(x) => x.residual(&ois),
        }
    }
}

pub fn instrument_residual_with_pillars<T: Real>(
    instrument: &RateInstrumentSpec,
    parameters: &[T],
    ois_pillars: usize,
) -> T {
    instrument.residual(parameters, ois_pillars)
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
        let mut out = Self::new(
            theta.clone(),
            helpers.iter().map(|h| h.quote_value()).collect(),
            jt,
            jq,
            vec![curve],
        );
        let times = out.curves[0].times.clone();
        out.instruments = helpers
            .iter()
            .filter_map(|h| {
                let kind = if let Some(x) = h.as_any().downcast_ref::<DepositInstrument>() {
                    RateInstrumentKind::Deposit(x.clone())
                } else if let Some(x) = h.as_any().downcast_ref::<FraInstrument>() {
                    RateInstrumentKind::Fra(x.clone())
                } else {
                    if let Some(x) = h.as_any().downcast_ref::<OisSwapInstrument>() {
                        RateInstrumentKind::OisSwap(x.clone())
                    } else {
                        let x = h.as_any().downcast_ref::<BasisSwapInstrument>()?;
                        RateInstrumentKind::BasisSwap(x.clone())
                    }
                };
                Some(RateInstrumentSpec {
                    kind,
                    reference_date,
                    ois_times: times.clone(),
                    fwd_times: Vec::new(),
                })
            })
            .collect();
        out.ois_pillars = helpers.len();
        Ok(out)
    }
    pub fn from_multicurve_helpers(
        reference_date: Date,
        ois_helpers: &[&dyn SingleCurveHelper],
        fwd_helpers: &[IrsSwapInstrument],
    ) -> Result<Self, BootstrapError> {
        let ois = bootstrap_piecewise_curve_helpers(reference_date, ois_helpers)?;
        let times: Vec<f64> = fwd_helpers
            .iter()
            .map(|h| {
                h.day_counter
                    .year_fraction(reference_date, h.end, None, None)
            })
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
        let ois_times = ois.times.clone();
        let fwd_times = fwd.times.clone();
        let mut instruments = Vec::with_capacity(k);
        for h in ois_helpers {
            let kind = if let Some(x) = h.as_any().downcast_ref::<DepositInstrument>() {
                RateInstrumentKind::Deposit(x.clone())
            } else if let Some(x) = h.as_any().downcast_ref::<FraInstrument>() {
                RateInstrumentKind::Fra(x.clone())
            } else if let Some(x) = h.as_any().downcast_ref::<OisSwapInstrument>() {
                RateInstrumentKind::OisSwap(x.clone())
            } else if let Some(x) = h.as_any().downcast_ref::<BasisSwapInstrument>() {
                RateInstrumentKind::BasisSwap(x.clone())
            } else {
                continue;
            };
            instruments.push(RateInstrumentSpec {
                kind,
                reference_date,
                ois_times: ois_times.clone(),
                fwd_times: fwd_times.clone(),
            });
        }
        instruments.extend(fwd_helpers.iter().cloned().map(|x| RateInstrumentSpec {
            kind: RateInstrumentKind::IrsSwap(x),
            reference_date,
            ois_times: ois_times.clone(),
            fwd_times: fwd_times.clone(),
        }));
        let mut out = Self::new(params, quotes, jt, jq, vec![ois, fwd]);
        out.instruments = instruments;
        out.ois_pillars = ois_helpers.len();
        Ok(out)
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
        let ois_pillars = curves.first().map_or(0, |c| c.times.len());
        Self {
            parameters,
            quotes,
            curves,
            sensitivity: IftSystemSensitivity::new(j_theta, j_q),
            instruments: Vec::new(),
            ois_pillars,
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
        let guard = TapeScope::guard();
        let vars: Vec<_> = self
            .parameters
            .iter()
            .map(|&x| TapeScope::with_active(|t| t.input(x)))
            .collect();
        let out = pricer(&vars);
        let adj = TapeScope::gradient_from(out);
        let pillar: Vec<f64> = vars.iter().map(|v| adj.gradient(*v)).collect();
        let value = out.val;
        let _ = guard.stop();
        (value, self.sensitivity.vjp(&pillar))
    }
    pub fn value_par_dv01_and_gamma<F>(&self, pricer: F) -> (f64, Vec<f64>, DMatrix<f64>)
    where
        F: Fn(&[TapeVar<f64>]) -> TapeVar<f64>,
    {
        let eval_grad = |params: &[f64]| {
            let guard = TapeScope::guard();
            let vars: Vec<_> = params
                .iter()
                .map(|&x| TapeScope::with_active(|t| t.input(x)))
                .collect::<Vec<_>>();
            let out = pricer(&vars);
            let adj = TapeScope::gradient_from(out);
            let g: Vec<f64> = vars.iter().map(|v| adj.gradient(*v)).collect();
            let value = out.val;
            let _ = guard.stop();
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
        let mut h_eff = h;
        if self.instruments.len() == k && k <= 64 {
            let lambda = self
                .sensitivity
                .lu_t
                .solve(&DVector::from_vec(g.clone()))
                .expect("singular IFT Jacobian");
            // A compact second-order jet keeps the calibration Hessian off the
            // stack in ordinary strips; very wide systems should use the
            // pairwise helper API rather than allocating an N² jet per node.
            let seeded: Vec<Dual2<64>> = self
                .parameters
                .iter()
                .enumerate()
                .map(|(i, &x)| Dual2::<64>::variable(x, i))
                .collect();
            for (m, instrument) in self.instruments.iter().enumerate() {
                let residual =
                    instrument_residual_with_pillars(instrument, &seeded, self.ois_pillars);
                for i in 0..k {
                    for jcol in 0..k {
                        h_eff[(i, jcol)] -= lambda[m] * residual.d2[i][jcol];
                    }
                }
            }
        } else if self.instruments.len() == k {
            let lambda = self
                .sensitivity
                .lu_t
                .solve(&DVector::from_vec(g.clone()))
                .expect("singular IFT Jacobian");
            for (m, instrument) in self.instruments.iter().enumerate() {
                for i in 0..k {
                    for jcol in i..k {
                        let seeded: Vec<Dual2<2>> = self
                            .parameters
                            .iter()
                            .enumerate()
                            .map(|(idx, &x)| {
                                if idx == i {
                                    Dual2::<2>::variable(x, 0)
                                } else if idx == jcol {
                                    Dual2::<2>::variable(x, if i == jcol { 0 } else { 1 })
                                } else {
                                    Dual2::from_value(x)
                                }
                            })
                            .collect();
                        let residual =
                            instrument_residual_with_pillars(instrument, &seeded, self.ois_pillars);
                        let hij = if i == jcol {
                            residual.d2[0][0]
                        } else {
                            residual.d2[0][1]
                        };
                        h_eff[(i, jcol)] -= lambda[m] * hij;
                        if i != jcol {
                            h_eff[(jcol, i)] -= lambda[m] * hij;
                        }
                    }
                }
            }
        }
        let gamma = j.transpose() * h_eff * j;
        let dv = self.sensitivity.vjp(&g);
        (value, dv, gamma)
    }
}
pub trait SingleCurveHelper: Send + Sync + Any {
    fn as_any(&self) -> &dyn Any;
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
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn pillar_time(&self, reference_date: Date) -> f64 {
        self.day_counter
            .year_fraction(reference_date, self.end, None, None)
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
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn pillar_time(&self, reference_date: Date) -> f64 {
        self.day_counter
            .year_fraction(reference_date, self.fixing_end, None, None)
    }
    fn quote_value(&self) -> f64 {
        self.quote
    }
    fn residual_curve(&self, c: &DiscountCurve<f64>) -> f64 {
        let t1 = self
            .day_counter
            .year_fraction(c.reference_date(), self.fixing_start, None, None);
        let t2 = self
            .day_counter
            .year_fraction(c.reference_date(), self.fixing_end, None, None);
        (c.discount_factor(t1) / c.discount_factor(t2) - 1.0)
            / self
                .day_counter
                .year_fraction(self.fixing_start, self.fixing_end, None, None)
            - self.quote
    }
    fn residual_curve_dual(&self, c: &DiscountCurve<Dual<1>>) -> Dual<1> {
        let t1 = self
            .day_counter
            .year_fraction(c.reference_date(), self.fixing_start, None, None);
        let t2 = self
            .day_counter
            .year_fraction(c.reference_date(), self.fixing_end, None, None);
        (c.discount_factor(t1) / c.discount_factor(t2) - Dual::one())
            / self
                .day_counter
                .year_fraction(self.fixing_start, self.fixing_end, None, None)
            - self.quote
    }
}
impl SingleCurveHelper for BasisSwapInstrument {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn pillar_time(&self, reference_date: Date) -> f64 {
        self.day_counter
            .year_fraction(reference_date, self.end, None, None)
    }
    fn quote_value(&self) -> f64 {
        self.spread
    }
    fn residual_curve(&self, c: &DiscountCurve<f64>) -> f64 {
        self.residual(c).val()
    }
    fn residual_curve_dual(&self, c: &DiscountCurve<Dual<1>>) -> Dual<1> {
        self.residual(c)
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
            let t = self
                .day_counter
                .year_fraction(projection.reference_date(), next, None, None);
            let p = discount.discount_factor(t);
            fixed += p * tau;
            let a = self
                .day_counter
                .year_fraction(projection.reference_date(), prev, None, None);
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
        let mut annuity = T::zero();
        for i in 1..=n {
            let t = DayCounter::Act365Fixed.year_fraction(
                discount.reference_date(),
                self.end,
                None,
                None,
            ) - (n - i) as f64 * 0.25;
            let p = discount.discount_factor(t);
            let a = (t - 0.25).max(0.0);
            annuity += p * 0.25;
            x += (leg1.discount_factor(a) / leg1.discount_factor(t)
                - leg2.discount_factor(a) / leg2.discount_factor(t))
                * p;
        }
        x / annuity.max_val(T::from_f64(1e-14)) - T::from_f64(self.spread)
    }
}
impl SingleCurveHelper for OisSwapInstrument {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn pillar_time(&self, reference_date: Date) -> f64 {
        self.day_counter
            .year_fraction(reference_date, self.end, None, None)
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
