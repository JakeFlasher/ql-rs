use ql_ad::{Real, TapeScope, TapeVar};
use ql_curves::YieldTermStructure;
use ql_math::brent_solve;
use ql_time::Date;
use ql_vol::{bachelier_price, black76_price};
use thiserror::Error;

#[derive(Clone, Debug)]
pub struct CochraneFiscalTheoryPriceLevelSpec<T: Real = f64> {
    pub debt_maturities_years: Vec<T>,
    pub nominal_debt_face_values: Vec<T>,
    pub nominal_zero_yields: Vec<T>,
    pub real_zero_yields: Vec<T>,
    pub primary_real_surpluses: Vec<T>,
    pub terminal_surplus_growth_rate: T,
    pub fiscal_surplus_shock_fraction: T,
    pub expected_future_inflation_rise: T,
    pub inflation_risk_beta: T,
    pub inflation_volatility: T,
    pub tips_liquidity_spread: T,
}

#[derive(Clone, Debug)]
pub struct CochraneFiscalTheoryPriceLevelReport<T: Real = f64> {
    pub nominal_debt_market_value: T,
    pub real_surplus_present_value: T,
    pub equilibrium_price_level_p0: T,
    pub fiscal_valuation_identity_error: T,
    pub debt_portfolio_duration_years: T,
    pub short_debt_immediate_inflation_jump: T,
    pub long_debt_smoothed_inflation_jump: T,
    pub maturity_structure_shock_absorption_ratio: T,
    pub model_breakeven_inflation_10y: T,
    pub inflation_risk_premium_10y: T,
}

pub fn evaluate_cochrane_ftpl_and_maturity_inflation_smoothing<T: Real>(
    s: &CochraneFiscalTheoryPriceLevelSpec<T>,
) -> CochraneFiscalTheoryPriceLevelReport<T> {
    let mut v = T::zero();
    let mut dur_num = T::zero();
    for i in 0..s.debt_maturities_years.len() {
        let tau = s.debt_maturities_years[i];
        let b = *s.nominal_debt_face_values.get(i).unwrap_or(&T::zero());
        let y = *s.nominal_zero_yields.get(i).unwrap_or(&T::zero());
        let pv = b * (-y * tau).exp();
        v += pv;
        dur_num += tau * pv;
    }
    let h = s.primary_real_surpluses.len();
    let mut pv_s = T::zero();
    for (i, &sur) in s.primary_real_surpluses.iter().enumerate() {
        let r = *s.real_zero_yields.get(i).unwrap_or(&T::zero());
        pv_s += sur / (T::one() + r).powf((i + 1) as f64);
    }
    if h > 0 {
        let r = *s.real_zero_yields.get(h - 1).unwrap_or(&T::zero());
        let den = (r - s.terminal_surplus_growth_rate).max_val(T::from_f64(1e-4));
        let sh = s.primary_real_surpluses[h - 1];
        pv_s += sh * (T::one() + s.terminal_surplus_growth_rate)
            / (den * (T::one() + r).powf(h as f64));
    }
    let p0 = v / pv_s.max_val(T::from_f64(1e-16));
    let err = (p0 * pv_s - v).abs();
    let shock = s.fiscal_surplus_shock_fraction;
    let short = T::one() / (T::one() + shock) - T::one();
    let mut vshock = T::zero();
    for i in 0..s.debt_maturities_years.len() {
        let tau = s.debt_maturities_years[i];
        let b = *s.nominal_debt_face_values.get(i).unwrap_or(&T::zero());
        let y = *s.nominal_zero_yields.get(i).unwrap_or(&T::zero());
        vshock += b * (-(y + s.expected_future_inflation_rise) * tau).exp();
    }
    let long = (vshock / v) / (T::one() + shock) - T::one();
    let ratio = T::one() - long / short;
    let irp =
        T::from_f64(10.0) * s.inflation_risk_beta * s.inflation_volatility * s.inflation_volatility;
    let yn = *s.nominal_zero_yields.last().unwrap_or(&T::zero());
    let yr = *s.real_zero_yields.last().unwrap_or(&T::zero());
    CochraneFiscalTheoryPriceLevelReport {
        nominal_debt_market_value: v,
        real_surplus_present_value: pv_s,
        equilibrium_price_level_p0: p0,
        fiscal_valuation_identity_error: err,
        debt_portfolio_duration_years: dur_num / v.max_val(T::from_f64(1e-16)),
        short_debt_immediate_inflation_jump: short,
        long_debt_smoothed_inflation_jump: long,
        maturity_structure_shock_absorption_ratio: ratio,
        model_breakeven_inflation_10y: yn - yr + irp - s.tips_liquidity_spread,
        inflation_risk_premium_10y: irp,
    }
}

#[cfg(test)]
mod phase27_tests {
    use super::*;
    use ql_ad::Dual;
    fn make<T: Real>(x: T) -> CochraneFiscalTheoryPriceLevelSpec<T> {
        CochraneFiscalTheoryPriceLevelSpec {
            debt_maturities_years: vec![T::one(), T::from_f64(5.0)],
            nominal_debt_face_values: vec![x, T::from_f64(50.0)],
            nominal_zero_yields: vec![T::from_f64(0.03), T::from_f64(0.04)],
            real_zero_yields: vec![T::from_f64(0.01), T::from_f64(0.015)],
            primary_real_surpluses: vec![T::from_f64(5.0), T::from_f64(5.2)],
            terminal_surplus_growth_rate: T::from_f64(0.02),
            fiscal_surplus_shock_fraction: T::from_f64(-0.05),
            expected_future_inflation_rise: T::from_f64(0.01),
            inflation_risk_beta: T::from_f64(0.4),
            inflation_volatility: T::from_f64(0.1),
            tips_liquidity_spread: T::from_f64(0.002),
        }
    }
    #[test]
    fn dual_fd_parity() {
        let h = 1e-5;
        let fd = (evaluate_cochrane_ftpl_and_maturity_inflation_smoothing(&make(100.0 + h))
            .equilibrium_price_level_p0
            - evaluate_cochrane_ftpl_and_maturity_inflation_smoothing(&make(100.0 - h))
                .equilibrium_price_level_p0)
            / (2.0 * h);
        let d = evaluate_cochrane_ftpl_and_maturity_inflation_smoothing(&make(
            Dual::<1>::variable(100.0, 0),
        ))
        .equilibrium_price_level_p0;
        assert!((d.d[0] - fd).abs() < 1e-5);
    }
}

#[derive(Debug, Error, Clone)]
pub enum InflationError {
    #[error("invalid seasonality factors")]
    InvalidSeasonality,
    #[error("invalid inflation curve")]
    InvalidCurve,
}

#[derive(Clone, Debug)]
pub struct BrigoMercurioInflationMarketModelSpec<T: Real = f64> {
    pub notional: T,
    pub base_cpi: T,
    pub payment_times: Vec<T>,
    pub forward_cpis: Vec<T>,
    pub nominal_forward_rates: Vec<T>,
    pub nominal_discount_factors: Vec<T>,
    pub forward_cpi_vols: Vec<T>,
    pub nominal_forward_vols: Vec<T>,
    pub inter_maturity_cpi_corr: T,
    pub nominal_cpi_corr: T,
    pub lpi_floor: T,
    pub lpi_cap: T,
    pub fixed_swap_rate: T,
}
#[derive(Clone, Debug)]
pub struct BrigoMercurioInflationMarketModelReport<T: Real = f64> {
    pub unadjusted_yoy_forwards: Vec<T>,
    pub convexity_adjusted_yoy_forwards: Vec<T>,
    pub total_yoy_convexity_adjustment_bps: T,
    pub yoy_effective_total_vols: Vec<T>,
    pub fimm_lpi_expected_rates: Vec<T>,
    pub unadjusted_lpi_swap_npv: T,
    pub fimm_convexity_adjusted_lpi_swap_npv: T,
    pub compounded_lpi_pension_terminal_index: T,
}
pub fn price_brigo_mercurio_yoy_and_lpi_market_model<T: Real>(
    s: &BrigoMercurioInflationMarketModelSpec<T>,
) -> BrigoMercurioInflationMarketModelReport<T> {
    let m = s.payment_times.len();
    let mut un = Vec::with_capacity(m);
    let mut adj = Vec::with_capacity(m);
    let mut vols = Vec::with_capacity(m);
    let mut lpi = Vec::with_capacity(m);
    let mut prev = s.base_cpi;
    let mut tprev = T::zero();
    let mut n_un = T::zero();
    let mut n_ad = T::zero();
    let mut terminal = s.base_cpi;
    for i in 0..m {
        let t = s.payment_times[i];
        let dt = (t - tprev).max_val(T::from_f64(1e-8));
        let cpi = *s.forward_cpis.get(i).unwrap_or(&prev);
        let g = cpi / prev;
        let f = g - T::one();
        let si = *s.forward_cpi_vols.get(i).unwrap_or(&T::zero());
        let sip = if i == 0 {
            T::zero()
        } else {
            *s.forward_cpi_vols.get(i - 1).unwrap_or(&T::zero())
        };
        let sn = *s.nominal_forward_vols.get(i).unwrap_or(&T::zero());
        let fnr = *s.nominal_forward_rates.get(i).unwrap_or(&T::zero());
        let c = if i == 0 {
            T::zero()
        } else {
            tprev
                * (sip * sip - s.inter_maturity_cpi_corr * si * sip
                    + (dt * fnr / (T::one() + dt * fnr)) * sn * s.nominal_cpi_corr * (si - sip))
        };
        let ga = g * c.exp();
        let fa = ga - T::one();
        let v = ((if i == 0 {
            T::zero()
        } else {
            tprev * (si * si + sip * sip - T::from_f64(2.0) * s.inter_maturity_cpi_corr * si * sip)
        }) + dt * si * si)
            .max_val(T::from_f64(1e-16))
            .sqrt();
        let lpi_rate = |gg: T, ff: T| {
            let kf = T::one() + s.lpi_floor;
            let kc = T::one() + s.lpi_cap;
            let d1f = ((gg / kf).ln() + v * v * 0.5) / v;
            let d2f = d1f - v;
            let put = kf * (-d2f).norm_cdf() - gg * (-d1f).norm_cdf();
            let d1c = ((gg / kc).ln() + v * v * 0.5) / v;
            let d2c = d1c - v;
            let call = gg * d1c.norm_cdf() - kc * d2c.norm_cdf();
            ff + put - call
        };
        let lu = lpi_rate(g, f);
        let la = lpi_rate(ga, fa);
        un.push(f);
        adj.push(fa);
        vols.push(v);
        lpi.push(la);
        let df = *s.nominal_discount_factors.get(i).unwrap_or(&T::one());
        n_un += s.notional * df * dt * (lu - s.fixed_swap_rate);
        n_ad += s.notional * df * dt * (la - s.fixed_swap_rate);
        terminal *= T::one() + la;
        prev = cpi;
        tprev = t;
    }
    let total = adj
        .iter()
        .zip(un.iter())
        .map(|(a, b)| (*a - *b) * 10000.0)
        .fold(T::zero(), |x, y| x + y);
    BrigoMercurioInflationMarketModelReport {
        unadjusted_yoy_forwards: un,
        convexity_adjusted_yoy_forwards: adj,
        total_yoy_convexity_adjustment_bps: total,
        yoy_effective_total_vols: vols,
        fimm_lpi_expected_rates: lpi,
        unadjusted_lpi_swap_npv: n_un,
        fimm_convexity_adjusted_lpi_swap_npv: n_ad,
        compounded_lpi_pension_terminal_index: terminal,
    }
}
#[derive(Clone, Debug)]
pub struct MultiplicativeSeasonality {
    pub factors: [f64; 12],
}
impl MultiplicativeSeasonality {
    pub fn new(factors: [f64; 12]) -> Result<Self, InflationError> {
        if factors.iter().any(|x| !x.is_finite() || *x <= 0.)
            || (factors.iter().product::<f64>() - 1.).abs() > 1e-6
        {
            return Err(InflationError::InvalidSeasonality);
        }
        Ok(Self { factors })
    }
    pub fn factor_for_month(&self, m: u8) -> f64 {
        self.factors[(m.clamp(1, 12) - 1) as usize]
    }
}
#[derive(Clone, Debug)]
pub struct ZeroInflationCurve<T: Real = f64> {
    pub base_date: Date,
    pub base_cpi: f64,
    pub pillar_times: Vec<f64>,
    pub zero_rates: Vec<T>,
    pub seasonality: Option<MultiplicativeSeasonality>,
}
impl<T: Real> ZeroInflationCurve<T> {
    pub fn new(
        base_date: Date,
        base_cpi: f64,
        pillar_times: Vec<f64>,
        zero_rates: Vec<T>,
        seasonality: Option<MultiplicativeSeasonality>,
    ) -> Self {
        assert_eq!(pillar_times.len(), zero_rates.len());
        Self {
            base_date,
            base_cpi,
            pillar_times,
            zero_rates,
            seasonality,
        }
    }
    fn z(&self, t: f64) -> T {
        if self.pillar_times.is_empty() {
            return T::zero();
        }
        // A zero curve is flat-extrapolated before its first quoted pillar.
        // Using the first interval's slope here creates an artificial shock
        // for lagged CPI dates just before the base pillar.
        if t <= self.pillar_times[0] {
            return self.zero_rates[0];
        }
        let i = if t <= self.pillar_times[0] {
            0
        } else if t >= *self.pillar_times.last().unwrap() {
            self.pillar_times.len() - 1
        } else {
            self.pillar_times.partition_point(|x| *x <= t) - 1
        };
        if i + 1 >= self.pillar_times.len() {
            self.zero_rates[i]
        } else {
            let w = (t - self.pillar_times[i]) / (self.pillar_times[i + 1] - self.pillar_times[i]);
            self.zero_rates[i] + (self.zero_rates[i + 1] - self.zero_rates[i]) * w
        }
    }
    pub fn forward_cpi(&self, t: f64, target_month: Option<u8>) -> T {
        let mut out = T::from_f64(self.base_cpi)
            * (T::one() + self.z(t))
                .ln()
                .mul_add(T::from_f64(t), T::zero())
                .exp();
        if let Some(s) = &self.seasonality {
            let m = target_month.unwrap_or_else(|| {
                self.base_date
                    .add_months((t * 12.0).round() as i32, false)
                    .month()
            });
            out *= T::from_f64(s.factor_for_month(m) / s.factor_for_month(self.base_date.month()));
        }
        out
    }
    pub fn interpolated_cpi(cpi_lag3: T, cpi_lag2: T, day: u8, dim: u8) -> T {
        if dim == 0 {
            return cpi_lag3;
        }
        cpi_lag3 + (cpi_lag2 - cpi_lag3) * (day.saturating_sub(1) as f64 / dim as f64)
    }
    pub fn to_ad<S: Real>(&self, nodes: &[S]) -> ZeroInflationCurve<S> {
        ZeroInflationCurve::new(
            self.base_date,
            self.base_cpi,
            self.pillar_times.clone(),
            nodes.to_vec(),
            self.seasonality.clone(),
        )
    }
}
#[derive(Clone, Debug)]
pub struct ZeroCouponInflationSwap {
    pub notional: f64,
    pub fixed_rate: f64,
    pub maturity_years: f64,
    pub target_month: Option<u8>,
    pub pay_fixed: bool,
}
impl ZeroCouponInflationSwap {
    pub fn npv<T: Real>(&self, d: &impl YieldTermStructure<T>, i: &ZeroInflationCurve<T>) -> T {
        let ratio = i.forward_cpi(self.maturity_years, self.target_month) / T::from_f64(i.base_cpi);
        let fixed = (T::one() + self.fixed_rate).powf(self.maturity_years) - T::one();
        let n = T::from_f64(self.notional) * d.discount_factor(self.maturity_years);
        let x = n * (ratio - T::one() - fixed);
        if self.pay_fixed { x } else { -x }
    }
    pub fn fair_rate<T: Real>(&self, i: &ZeroInflationCurve<T>) -> T {
        if self.maturity_years <= 0.0 {
            return T::zero();
        }
        (i.forward_cpi(self.maturity_years, self.target_month) / T::from_f64(i.base_cpi))
            .powf(1. / self.maturity_years)
            - T::one()
    }
}
#[derive(Clone, Debug)]
pub struct YearOnYearInflationSwap {
    pub notional: f64,
    pub fixed_rate: f64,
    pub maturity_years: usize,
    pub pay_fixed: bool,
}
impl YearOnYearInflationSwap {
    pub fn npv<T: Real>(&self, d: &impl YieldTermStructure<T>, i: &ZeroInflationCurve<T>) -> T {
        let mut x = T::zero();
        let mut prev = T::from_f64(i.base_cpi);
        for j in 1..=self.maturity_years {
            let cur = i.forward_cpi(j as f64, None);
            x += d.discount_factor(j as f64) * (cur / prev - T::one() - self.fixed_rate);
            prev = cur;
        }
        let x = x * T::from_f64(self.notional);
        if self.pay_fixed { x } else { -x }
    }
}
#[derive(Clone, Debug)]
pub struct CalibratedInflationCurve {
    pub curve: ZeroInflationCurve<f64>,
    pub quotes: Vec<(f64, f64)>,
    pub sensitivities: nalgebra::DMatrix<f64>,
    pub discount: ql_curves::DiscountCurve<f64>,
}
impl CalibratedInflationCurve {
    pub fn bootstrap(
        base: Date,
        base_cpi: f64,
        d: ql_curves::DiscountCurve<f64>,
        quotes: &[(f64, f64)],
    ) -> Self {
        Self::bootstrap_with_seasonality(base, base_cpi, d, quotes, None)
            .expect("valid inflation quotes")
    }

    /// Bootstrap zero-inflation pillars and cache the exact IFT quote
    /// sensitivities.  The one-dimensional solves are sequential because a
    /// pillar's CPI is determined by its own zero rate once earlier pillars
    /// have been fixed.
    pub fn bootstrap_with_seasonality(
        base: Date,
        base_cpi: f64,
        d: ql_curves::DiscountCurve<f64>,
        quotes: &[(f64, f64)],
        seasonality: Option<MultiplicativeSeasonality>,
    ) -> Result<Self, InflationError> {
        if quotes.is_empty()
            || quotes.iter().any(|(t, _)| *t <= 0.0)
            || quotes.windows(2).any(|w| w[1].0 <= w[0].0)
        {
            return Err(InflationError::InvalidCurve);
        }
        let times: Vec<f64> = quotes.iter().map(|(t, _)| *t).collect();
        let mut rates = vec![0.02; quotes.len()];
        for k in 0..quotes.len() {
            let prior = rates.clone();
            let maturity = quotes[k].0;
            let fixed = quotes[k].1;
            let f = |z: f64| {
                let mut trial = prior.clone();
                trial[k] = z;
                let curve = ZeroInflationCurve::new(
                    base,
                    base_cpi,
                    times.clone(),
                    trial,
                    seasonality.clone(),
                );
                let swap = ZeroCouponInflationSwap {
                    notional: 1.0,
                    fixed_rate: fixed,
                    maturity_years: maturity,
                    target_month: None,
                    pay_fixed: true,
                };
                swap.npv(&d, &curve).val()
            };
            rates[k] =
                brent_solve(f, -0.99, 2.0, 1e-14, 200).map_err(|_| InflationError::InvalidCurve)?;
        }

        let curve = ZeroInflationCurve::new(
            base,
            base_cpi,
            times.clone(),
            rates.clone(),
            seasonality.clone(),
        );
        let m = quotes.len();
        let mut jz = nalgebra::DMatrix::zeros(m, m);
        for j in 0..m {
            let ad_rates: Vec<ql_ad::Dual<1>> = rates
                .iter()
                .enumerate()
                .map(|(i, &x)| {
                    if i == j {
                        ql_ad::Dual::<1>::variable(x, 0)
                    } else {
                        ql_ad::Dual::<1>::from_value(x)
                    }
                })
                .collect();
            let ad_curve = curve.to_ad(&ad_rates);
            let ad_discount = d.to_ad_curve(
                &d.discounts
                    .iter()
                    .map(|&x| ql_ad::Dual::<1>::from_value(x))
                    .collect::<Vec<_>>(),
            );
            for i in 0..m {
                let swap = ZeroCouponInflationSwap {
                    notional: 1.0,
                    fixed_rate: quotes[i].1,
                    maturity_years: quotes[i].0,
                    target_month: None,
                    pay_fixed: true,
                };
                jz[(i, j)] = swap.npv(&ad_discount, &ad_curve).d[0];
            }
        }
        let mut jk = nalgebra::DMatrix::zeros(m, m);
        for i in 0..m {
            let t = quotes[i].0;
            jk[(i, i)] = -d.discount_factor(t).val() * t * (1.0 + quotes[i].1).powf(t - 1.0) * 1.0;
        }
        let sensitivities = jz.lu().solve(&(-jk)).ok_or(InflationError::InvalidCurve)?;
        Ok(Self {
            curve,
            quotes: quotes.to_vec(),
            sensitivities,
            discount: d,
        })
    }
    pub fn value_and_par_inflation_delta<F>(&self, f: F) -> (f64, Vec<f64>)
    where
        F: Fn(&ZeroInflationCurve<TapeVar>) -> TapeVar,
    {
        let guard = TapeScope::guard();
        let vars: Vec<_> = self
            .curve
            .zero_rates
            .iter()
            .map(|&x| TapeScope::with_active(|t| t.input(x)))
            .collect();
        let curve = self.curve.to_ad(&vars);
        let out = f(&curve);
        let adj = TapeScope::gradient_from(out);
        let pillar: Vec<f64> = vars.iter().map(|v| adj.gradient(*v)).collect();
        let _ = guard.stop();
        let q = &self.sensitivities.transpose() * nalgebra::DVector::from_vec(pillar);
        (out.val, q.iter().copied().collect())
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum YoYInflationCapFloorType {
    Cap,
    Floor,
    Collar,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum YoYOptionletModel {
    Black,
    UnitDisplacedBlack,
    Bachelier,
}

#[derive(Copy, Clone, Debug)]
pub struct YoYInflationCapFloor<T: Real = f64> {
    pub notional: T,
    pub maturity_years: usize,
    pub cap_rate: Option<T>,
    pub floor_rate: Option<T>,
    pub gearing: T,
    pub spread: T,
    pub cap_floor_type: YoYInflationCapFloorType,
}

/// Three-factor Jarrow--Yildirim nominal/real/CPI Gaussian hybrid.
#[derive(Copy, Clone, Debug)]
pub struct JarrowYildirimParams<T: Real = f64> {
    pub a_nom: T,
    pub sigma_nom: T,
    pub a_real: T,
    pub sigma_real: T,
    pub sigma_cpi: T,
    pub rho_nom_real: T,
    pub rho_nom_cpi: T,
    pub rho_real_cpi: T,
}

fn jy_b<T: Real>(a: T, tau: f64) -> T {
    let t = T::from_f64(tau.max(0.0));
    let at = a * t;
    if at.val().abs() < 1e-6 {
        t * (T::one() - at * 0.5 + at * at / 6.0)
    } else {
        (T::one() - (-at).exp()) / a
    }
}
fn jy_psi<T: Real>(a1: T, a2: T, tau: f64) -> T {
    let t = tau.max(0.0);
    if t <= 0.0 {
        return T::zero();
    }
    let nodes = [
        -0.9602898565,
        -0.7966664774,
        -0.5255324099,
        -0.1834346425,
        0.1834346425,
        0.5255324099,
        0.7966664774,
        0.9602898565,
    ];
    let weights = [
        0.1012285363,
        0.2223810345,
        0.3137066459,
        0.3626837834,
        0.3626837834,
        0.3137066459,
        0.2223810345,
        0.1012285363,
    ];
    let mut sum = T::zero();
    for (&x, &w) in nodes.iter().zip(weights.iter()) {
        let s = 0.5 * t * (x + 1.0);
        sum += jy_b(a1, s) * jy_b(a2, s) * (0.5 * t * w);
    }
    sum
}
fn jy_phi<T: Real>(a: T, tau: f64) -> T {
    jy_psi(a, a, tau)
}
fn jy_int_b<T: Real>(a: T, tau: f64) -> T {
    let t = T::from_f64(tau.max(0.0));
    let at = a * t;
    if at.val().abs() < 1e-5 {
        t * t * (T::from_f64(0.5) - at / 6.0 + at * at / 24.0)
    } else {
        (t - jy_b(a, tau)) / a
    }
}
pub fn jarrow_yildirim_yoy_moments<T: Real>(
    p: &JarrowYildirimParams<T>,
    t_start: f64,
    t_end: f64,
) -> (T, T) {
    let t0 = t_start.max(0.0);
    let dt = (t_end - t0).max(0.0);
    if dt <= 0.0 {
        return (T::zero(), T::zero());
    }
    let bn = jy_b(p.a_nom, dt);
    let br = jy_b(p.a_real, dt);
    let pre = if t0 > 0.0 {
        p.sigma_nom * p.sigma_nom * bn * bn * jy_b(p.a_nom * 2.0, t0)
            + p.sigma_real * p.sigma_real * br * br * jy_b(p.a_real * 2.0, t0)
            - p.rho_nom_real
                * p.sigma_nom
                * p.sigma_real
                * bn
                * br
                * 2.0
                * jy_b(p.a_nom + p.a_real, t0)
    } else {
        T::zero()
    };
    let acc = p.sigma_cpi * p.sigma_cpi * dt
        + p.sigma_nom * p.sigma_nom * jy_phi(p.a_nom, dt)
        + p.sigma_real * p.sigma_real * jy_phi(p.a_real, dt)
        - p.rho_nom_real * p.sigma_nom * p.sigma_real * jy_psi(p.a_nom, p.a_real, dt) * 2.0
        + p.rho_nom_cpi * p.sigma_nom * p.sigma_cpi * jy_int_b(p.a_nom, dt) * 2.0
        - p.rho_real_cpi * p.sigma_real * p.sigma_cpi * jy_int_b(p.a_real, dt) * 2.0;
    let var = (pre + acc).max_val(T::zero());
    let conv = if t0 > 0.0 {
        let nodes = [
            -0.9602898565,
            -0.7966664774,
            -0.5255324099,
            -0.1834346425,
            0.1834346425,
            0.5255324099,
            0.7966664774,
            0.9602898565,
        ];
        let weights = [
            0.1012285363,
            0.2223810345,
            0.3137066459,
            0.3626837834,
            0.3626837834,
            0.3137066459,
            0.2223810345,
            0.1012285363,
        ];
        let mut integral = T::zero();
        for (&x, &w) in nodes.iter().zip(weights.iter()) {
            let s = 0.5 * t0 * (x + 1.0);
            let rem0 = t0 - s;
            let rem1 = t_end - s;
            let term = p.sigma_real * p.sigma_real * jy_b(p.a_real, rem0)
                - p.rho_real_cpi * p.sigma_real * p.sigma_cpi
                - p.rho_nom_real * p.sigma_nom * p.sigma_real * jy_b(p.a_nom, rem1);
            integral += term * (-p.a_real * rem0).exp() * (0.5 * t0 * w);
        }
        br * integral
    } else {
        T::zero()
    };
    (conv, var)
}
pub fn price_yoy_capfloor_jarrow_yildirim<T: Real>(
    instrument: &YoYInflationCapFloor<T>,
    discount: &impl YieldTermStructure<T>,
    inflation: &ZeroInflationCurve<T>,
    params: &JarrowYildirimParams<T>,
) -> (T, Vec<T>) {
    if instrument.maturity_years == 0
        || !instrument.gearing.val().is_finite()
        || instrument.gearing.val() < 0.0
        || !instrument.spread.val().is_finite()
    {
        return (T::zero(), Vec::new());
    }
    let mut cap_npv = T::zero();
    let mut floor_npv = T::zero();
    let mut fwds = Vec::with_capacity(instrument.maturity_years);
    let g = instrument.gearing;
    for y in 1..=instrument.maturity_years {
        let t0 = (y - 1) as f64;
        let t1 = y as f64;
        let ratio = inflation.forward_cpi(t1, None) / inflation.forward_cpi(t0, None);
        let (conv, var) = jarrow_yildirim_yoy_moments(params, t0, t1);
        let gross = ratio * conv.exp();
        fwds.push(gross - 1.0);
        let df = discount.discount_factor(t1);
        let price_let = |k_rate: T, omega: f64| {
            if g.val().abs() < 1e-14 {
                return (T::from_f64(omega) * (instrument.spread - k_rate)).max_val(T::zero());
            }
            let strike = ((k_rate - instrument.spread) / g + 1.0).max_val(T::from_f64(1e-8));
            let st = var.sqrt();
            let unit = if st.val() <= 1e-14 {
                (T::from_f64(omega) * (gross - strike)).max_val(T::zero())
            } else {
                let d1 = ((gross / strike).ln() + var * 0.5) / st;
                let d2 = d1 - st;
                T::from_f64(omega)
                    * (gross * (T::from_f64(omega) * d1).norm_cdf()
                        - strike * (T::from_f64(omega) * d2).norm_cdf())
            };
            g * unit
        };
        if let Some(k) = instrument.cap_rate {
            cap_npv += instrument.notional * df * price_let(k, 1.0);
        }
        if let Some(k) = instrument.floor_rate {
            floor_npv += instrument.notional * df * price_let(k, -1.0);
        }
    }
    let npv = match instrument.cap_floor_type {
        YoYInflationCapFloorType::Cap => cap_npv,
        YoYInflationCapFloorType::Floor => floor_npv,
        YoYInflationCapFloorType::Collar => cap_npv - floor_npv,
    };
    (npv, fwds)
}

fn yoy_optionlet<T: Real>(
    f: T,
    k: T,
    sigma: T,
    t: T,
    df: T,
    model: YoYOptionletModel,
    cap: bool,
) -> T {
    match model {
        YoYOptionletModel::Black => {
            // Black is undefined for non-positive inflation forwards/strikes;
            // use the continuous Bachelier extension in that regime.
            if f.val() <= 0.0 || k.val() <= 0.0 {
                bachelier_price(f, k, sigma, t, df, cap)
            } else {
                black76_price(f, k, sigma, t, df, cap)
            }
        }
        YoYOptionletModel::UnitDisplacedBlack => {
            let fd = f + T::one();
            let kd = k + T::one();
            if fd.val() <= 0.0 || kd.val() <= 0.0 {
                bachelier_price(f, k, sigma, t, df, cap)
            } else {
                black76_price(fd, kd, sigma, t, df, cap)
            }
        }
        YoYOptionletModel::Bachelier => bachelier_price(f, k, sigma, t, df, cap),
    }
}

pub fn price_yoy_inflation_capfloor<T: Real>(
    instrument: &YoYInflationCapFloor<T>,
    discount: &impl YieldTermStructure<T>,
    inflation: &ZeroInflationCurve<T>,
    vol: T,
    model: YoYOptionletModel,
) -> T {
    if instrument.maturity_years == 0
        || !instrument.gearing.val().is_finite()
        || instrument.gearing.val() < 0.0
        || !instrument.spread.val().is_finite()
        || !vol.val().is_finite()
        || vol.val() < 0.0
        || instrument
            .cap_rate
            .zip(instrument.floor_rate)
            .map(|(cap, floor)| cap.val() < floor.val())
            .unwrap_or(false)
    {
        return T::zero();
    }
    if instrument.gearing.val().abs() < 1e-14 {
        let mut cap = T::zero();
        let mut floor = T::zero();
        for j in 1..=instrument.maturity_years {
            let df = discount.discount_factor(j as f64);
            let spread = instrument.spread;
            if let Some(k) = instrument.cap_rate {
                cap += instrument.notional * df * (spread - k).max_val(T::zero());
            }
            if let Some(k) = instrument.floor_rate {
                floor += instrument.notional * df * (k - spread).max_val(T::zero());
            }
        }
        return match instrument.cap_floor_type {
            YoYInflationCapFloorType::Cap => cap,
            YoYInflationCapFloorType::Floor => floor,
            YoYInflationCapFloorType::Collar => cap - floor,
        };
    }
    let mut cap = T::zero();
    let mut floor = T::zero();
    let mut prev = T::from_f64(inflation.base_cpi);
    for j in 1..=instrument.maturity_years {
        let cur = inflation.forward_cpi(j as f64, None);
        let f = cur / prev - T::one();
        let df = discount.discount_factor(j as f64);
        let t = T::from_f64(j as f64);
        // The quoted option is on the un-geared YoY rate.  Rewriting
        // g(F-K)+s = g(F-(K-s)/g) keeps the option model and cap/floor
        // parity consistent for non-zero gearing and spread.
        let g = instrument.gearing;
        let k_eff = |k: T| {
            if g.val().abs() < 1e-14 {
                T::zero()
            } else {
                (k - instrument.spread) / g
            }
        };
        if let Some(k) = instrument.cap_rate {
            cap += instrument.notional * g * yoy_optionlet(f, k_eff(k), vol, t, df, model, true);
        }
        if let Some(k) = instrument.floor_rate {
            floor += instrument.notional * g * yoy_optionlet(f, k_eff(k), vol, t, df, model, false);
        }
        prev = cur;
    }
    match instrument.cap_floor_type {
        YoYInflationCapFloorType::Cap => cap,
        YoYInflationCapFloorType::Floor => floor,
        YoYInflationCapFloorType::Collar => cap - floor,
    }
}

#[derive(Copy, Clone, Debug)]
pub struct LpiInflationSwapSpec<T: Real = f64> {
    pub notional: T,
    pub fixed_rate: T,
    pub floor_rate: T,
    pub cap_rate: T,
    pub maturity_years: usize,
    pub is_compounded_zero_coupon: bool,
    pub pay_fixed: bool,
}
pub fn lpi_inflation_swap_npv<T: Real>(
    s: &LpiInflationSwapSpec<T>,
    discount: &impl YieldTermStructure<T>,
    inflation: &ZeroInflationCurve<T>,
    vol: T,
    model: YoYOptionletModel,
) -> (T, T) {
    if s.maturity_years == 0
        || s.cap_rate.val() < s.floor_rate.val()
        || !vol.val().is_finite()
        || vol.val() < 0.0
    {
        return (T::zero(), T::zero());
    }
    if s.is_compounded_zero_coupon {
        let mut index_ratio = T::one();
        for j in 1..=s.maturity_years {
            let f = inflation.forward_cpi(j as f64, None)
                / inflation.forward_cpi((j - 1) as f64, None)
                - T::one();
            let floorlet = yoy_optionlet(
                f,
                s.floor_rate,
                vol,
                T::from_f64(j as f64),
                T::one(),
                model,
                false,
            );
            let caplet = yoy_optionlet(
                f,
                s.cap_rate,
                vol,
                T::from_f64(j as f64),
                T::one(),
                model,
                true,
            );
            let y = (f + floorlet - caplet)
                .max_val(s.floor_rate)
                .min_val(s.cap_rate);
            index_ratio *= T::one() + y;
        }
        let t = s.maturity_years as f64;
        let fixed_ratio = (T::one() + s.fixed_rate).powf(t);
        let df = discount.discount_factor(t);
        let npv = s.notional * df * (index_ratio - fixed_ratio);
        let fair = index_ratio.powf(1.0 / t) - T::one();
        return (if s.pay_fixed { npv } else { -npv }, fair);
    }
    let mut pv = T::zero();
    let mut annuity = T::zero();
    let mut weighted = T::zero();
    for j in 1..=s.maturity_years {
        let f = inflation.forward_cpi(j as f64, None) / inflation.forward_cpi((j - 1) as f64, None)
            - T::one();
        let floorlet = yoy_optionlet(
            f,
            s.floor_rate,
            vol,
            T::from_f64(j as f64),
            T::one(),
            model,
            false,
        );
        let caplet = yoy_optionlet(
            f,
            s.cap_rate,
            vol,
            T::from_f64(j as f64),
            T::one(),
            model,
            true,
        );
        let y = (f + floorlet - caplet)
            .max_val(s.floor_rate)
            .min_val(s.cap_rate);
        let df = discount.discount_factor(j as f64);
        pv += df * (y - s.fixed_rate);
        annuity += df;
        weighted += df * y;
    }
    let fair = weighted / annuity.max_val(T::from_f64(1e-14));
    let npv = pv * s.notional;
    (if s.pay_fixed { npv } else { -npv }, fair)
}

#[derive(Copy, Clone, Debug)]
pub struct ZciisSwaptionSpec<T: Real = f64> {
    pub notional: T,
    pub strike_rate: T,
    pub option_expiry_years: f64,
    pub swap_maturity_years: f64,
    pub target_month: Option<u8>,
    pub is_payer: bool,
}

#[derive(Clone, Debug)]
pub struct YoYInflationSwaptionSpec<T: Real = f64> {
    pub notional: T,
    pub strike_rate: T,
    pub option_expiry_years: f64,
    pub start_year: usize,
    pub end_year: usize,
    pub is_payer: bool,
}

pub fn price_zciis_swaption<T: Real>(
    spec: &ZciisSwaptionSpec<T>,
    discount: &impl YieldTermStructure<T>,
    inflation: &ZeroInflationCurve<T>,
    vol: T,
    model: YoYOptionletModel,
) -> (T, T) {
    let te = spec.option_expiry_years;
    let tm = spec.swap_maturity_years;
    if te <= 0.0 || tm <= 0.0 || te > tm || !vol.val().is_finite() || vol.val() < 0.0 {
        return (T::zero(), T::zero());
    }
    let ratio = inflation.forward_cpi(tm, spec.target_month) / T::from_f64(inflation.base_cpi);
    let z = ratio.powf(1.0 / tm) - T::one();
    let dz = z - spec.strike_rate;
    let comp_factor = if dz.val().abs() > 1e-7 {
        ((T::one() + z).powf(tm) - (T::one() + spec.strike_rate).powf(tm)) / dz
    } else {
        T::from_f64(tm) * (T::one() + z).powf(tm - 1.0)
            - T::from_f64(0.5 * tm * (tm - 1.0)) * (T::one() + z).powf(tm - 2.0) * dz
    };
    let ann = spec.notional * comp_factor * discount.discount_factor(tm);
    (
        yoy_optionlet(
            z,
            spec.strike_rate,
            vol,
            T::from_f64(te),
            ann,
            model,
            spec.is_payer,
        ),
        z,
    )
}

pub fn price_yoy_inflation_swaption<T: Real>(
    spec: &YoYInflationSwaptionSpec<T>,
    discount: &impl YieldTermStructure<T>,
    inflation: &ZeroInflationCurve<T>,
    vol: T,
    model: YoYOptionletModel,
) -> (T, T) {
    if spec.end_year <= spec.start_year
        || spec.option_expiry_years <= 0.0
        || spec.option_expiry_years > spec.start_year as f64 + 1e-12
        || !vol.val().is_finite()
        || vol.val() < 0.0
    {
        return (T::zero(), T::zero());
    }
    let mut ann = T::zero();
    let mut weighted = T::zero();
    for y in (spec.start_year + 1)..=spec.end_year {
        let df = discount.discount_factor(y as f64);
        let f = inflation.forward_cpi(y as f64, None) / inflation.forward_cpi((y - 1) as f64, None)
            - T::one();
        ann += df;
        weighted += df * f;
    }
    let fwd = weighted / ann.max_val(T::from_f64(1e-14));
    (
        yoy_optionlet(
            fwd,
            spec.strike_rate,
            vol,
            T::from_f64(spec.option_expiry_years),
            spec.notional * ann,
            model,
            spec.is_payer,
        ),
        fwd,
    )
}

pub fn estimate_multiplicative_seasonality(
    observations: &[(u8, f64)],
) -> Result<MultiplicativeSeasonality, InflationError> {
    if observations.len() < 13
        || observations
            .iter()
            .any(|(m, x)| *m == 0 || *m > 12 || !x.is_finite() || *x <= 0.0)
    {
        return Err(InflationError::InvalidSeasonality);
    }
    let mut seen = [false; 13];
    let mut logs = Vec::with_capacity(observations.len() - 1);
    for w in observations.windows(2) {
        let g = (w[1].1 / w[0].1).ln();
        logs.push((w[1].0, g));
        seen[w[1].0 as usize] = true;
    }
    if (1..=12).any(|m| !seen[m]) {
        return Err(InflationError::InvalidSeasonality);
    }
    let mean = logs.iter().map(|(_, g)| *g).sum::<f64>() / logs.len() as f64;
    let mut gamma = [0.0; 12];
    let mut count = [0usize; 12];
    for (m, g) in logs {
        gamma[m as usize - 1] += g - mean;
        count[m as usize - 1] += 1;
    }
    for i in 0..12 {
        gamma[i] /= count[i].max(1) as f64;
    }
    let center = gamma.iter().sum::<f64>() / 12.0;
    let mut factors = [0.0; 12];
    for i in 0..12 {
        factors[i] = (gamma[i] - center).exp();
    }
    MultiplicativeSeasonality::new(factors)
}

#[derive(Copy, Clone, Debug)]
pub struct ZciisCapFloorSpec<T: Real = f64> {
    pub notional: T,
    pub maturity_years: f64,
    pub strike_annual_rate: T,
    pub target_month: Option<u8>,
    pub is_cap: bool,
}

#[derive(Copy, Clone, Debug)]
pub struct ZciisCapFloorValuation<T: Real> {
    pub npv: T,
    pub forward_cpi_ratio: T,
    pub strike_growth_factor: T,
    pub log_cpi_variance: T,
    pub implied_annual_inflation_forward: T,
}

pub fn price_zciis_capfloor_jarrow_yildirim<T: Real>(
    spec: &ZciisCapFloorSpec<T>,
    discount: &impl YieldTermStructure<T>,
    inflation: &ZeroInflationCurve<T>,
    params: &JarrowYildirimParams<T>,
) -> ZciisCapFloorValuation<T> {
    let t = spec.maturity_years.max(0.0);
    if t <= 0.0 {
        return ZciisCapFloorValuation {
            npv: T::zero(),
            forward_cpi_ratio: T::one(),
            strike_growth_factor: T::one(),
            log_cpi_variance: T::zero(),
            implied_annual_inflation_forward: T::zero(),
        };
    }
    let df = discount.discount_factor(t);
    let f_ratio = (inflation.forward_cpi(t, spec.target_month)
        / T::from_f64(inflation.base_cpi.max(1e-12)))
    .max_val(T::from_f64(1e-8));
    let one_plus_k = (T::one() + spec.strike_annual_rate).max_val(T::from_f64(1e-6));
    let k_tot = (one_plus_k.ln() * t).exp();
    let (_, var_cpi) = jarrow_yildirim_yoy_moments(params, 0.0, t);
    let std = var_cpi.max_val(T::from_f64(1e-16)).sqrt();
    let d1 = ((f_ratio / k_tot).ln() + var_cpi * 0.5) / std;
    let d2 = d1 - std;
    let undiscounted = if spec.is_cap {
        f_ratio * d1.norm_cdf() - k_tot * d2.norm_cdf()
    } else {
        k_tot * (-d2).norm_cdf() - f_ratio * (-d1).norm_cdf()
    };
    let imp_ann = (f_ratio.ln() / t).exp() - T::one();
    ZciisCapFloorValuation {
        npv: spec.notional * df * undiscounted,
        forward_cpi_ratio: f_ratio,
        strike_growth_factor: k_tot,
        log_cpi_variance: var_cpi,
        implied_annual_inflation_forward: imp_ann,
    }
}

#[cfg(test)]
mod phase14_tests {
    use super::*;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;
    #[test]
    fn phase14_inflation_swaption_and_seasonality() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let disc = FlatForwardCurve::new(d, 0.02);
        let inf = ZeroInflationCurve::new(
            d,
            100.0,
            vec![1.0, 2.0, 5.0],
            vec![0.02, 0.021, 0.023],
            None,
        );
        let z = ZciisSwaptionSpec {
            notional: 1.0,
            strike_rate: 0.02,
            option_expiry_years: 1.0,
            swap_maturity_years: 3.0,
            target_month: None,
            is_payer: true,
        };
        let (v, f) = price_zciis_swaption(&z, &disc, &inf, 0.01, YoYOptionletModel::Bachelier);
        assert!(v.is_finite() && f.is_finite());
        let obs: Vec<(u8, f64)> = (0..=12)
            .map(|i| ((i % 12 + 1) as u8, 100.0 + i as f64))
            .collect();
        let s = estimate_multiplicative_seasonality(&obs).unwrap();
        assert!((s.factors.iter().product::<f64>() - 1.0).abs() < 1e-10);
    }
}

#[cfg(test)]
mod phase15_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;
    #[test]
    fn jarrow_yildirim_moments_parity_and_ad() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let disc = FlatForwardCurve::new(d, 0.025);
        let inf = ZeroInflationCurve::new(d, 100.0, vec![1.0, 2.0, 5.0], vec![0.02; 3], None);
        let p = JarrowYildirimParams {
            a_nom: 0.1,
            sigma_nom: 0.01,
            a_real: 0.08,
            sigma_real: 0.008,
            sigma_cpi: 0.012,
            rho_nom_real: 0.2,
            rho_nom_cpi: 0.15,
            rho_real_cpi: -0.2,
        };
        let mut p_pos = p;
        p_pos.rho_nom_cpi = 0.5;
        let mut p_neg = p;
        p_neg.rho_nom_cpi = -0.5;
        assert!(
            jarrow_yildirim_yoy_moments(&p_pos, 1.0, 2.0).1
                > jarrow_yildirim_yoy_moments(&p_neg, 1.0, 2.0).1
        );
        let cap = YoYInflationCapFloor {
            notional: 100.0,
            maturity_years: 3,
            cap_rate: Some(0.022),
            floor_rate: Some(0.022),
            gearing: 1.5,
            spread: 0.002,
            cap_floor_type: YoYInflationCapFloorType::Cap,
        };
        let mut flr = cap;
        flr.cap_floor_type = YoYInflationCapFloorType::Floor;
        let mut col = cap;
        col.cap_floor_type = YoYInflationCapFloorType::Collar;
        let (vc, _) = price_yoy_capfloor_jarrow_yildirim(&cap, &disc, &inf, &p);
        let (vf, _) = price_yoy_capfloor_jarrow_yildirim(&flr, &disc, &inf, &p);
        let (vcol, _) = price_yoy_capfloor_jarrow_yildirim(&col, &disc, &inf, &p);
        assert!((vc - vf - vcol).abs() < 1e-12);
        let h = 1e-5;
        let mut p_up = p;
        let mut p_dn = p;
        p_up.sigma_real += h;
        p_dn.sigma_real -= h;
        let fd = (price_yoy_capfloor_jarrow_yildirim(&cap, &disc, &inf, &p_up).0
            - price_yoy_capfloor_jarrow_yildirim(&cap, &disc, &inf, &p_dn).0)
            / (2.0 * h);
        type D = Dual<1>;
        let disc_d = FlatForwardCurve::new(d, D::from_value(0.025));
        let inf_d = ZeroInflationCurve::new(
            d,
            100.0,
            vec![1.0, 2.0, 5.0],
            vec![D::from_value(0.02); 3],
            None,
        );
        let cap_d = YoYInflationCapFloor {
            notional: D::from_value(100.0),
            maturity_years: 3,
            cap_rate: Some(D::from_value(0.022)),
            floor_rate: Some(D::from_value(0.022)),
            gearing: D::from_value(1.5),
            spread: D::from_value(0.002),
            cap_floor_type: YoYInflationCapFloorType::Cap,
        };
        let pd = JarrowYildirimParams {
            a_nom: D::from_value(0.1),
            sigma_nom: D::from_value(0.01),
            a_real: D::from_value(0.08),
            sigma_real: D::variable(0.008, 0),
            sigma_cpi: D::from_value(0.012),
            rho_nom_real: D::from_value(0.2),
            rho_nom_cpi: D::from_value(0.15),
            rho_real_cpi: D::from_value(-0.2),
        };
        let ad = price_yoy_capfloor_jarrow_yildirim(&cap_d, &disc_d, &inf_d, &pd).0;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase19_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;

    #[test]
    fn zciis_capfloor_jarrow_yildirim_parity_and_ad() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let disc = FlatForwardCurve::new(d, 0.03);
        let infl = ZeroInflationCurve::new(
            d,
            100.0,
            vec![1.0, 3.0, 5.0],
            vec![0.022, 0.024, 0.025],
            None,
        );
        let jy = JarrowYildirimParams {
            a_nom: 0.05,
            sigma_nom: 0.008,
            a_real: 0.04,
            sigma_real: 0.006,
            sigma_cpi: 0.012,
            rho_nom_real: 0.3,
            rho_nom_cpi: 0.2,
            rho_real_cpi: -0.1,
        };
        let cap_spec = ZciisCapFloorSpec {
            notional: 100.0,
            maturity_years: 5.0,
            strike_annual_rate: 0.02,
            target_month: None,
            is_cap: true,
        };
        let mut flr_spec = cap_spec;
        flr_spec.is_cap = false;
        let cap = price_zciis_capfloor_jarrow_yildirim(&cap_spec, &disc, &infl, &jy);
        let flr = price_zciis_capfloor_jarrow_yildirim(&flr_spec, &disc, &infl, &jy);
        let h = 1e-5;
        assert!(cap.log_cpi_variance > 0.0);
        assert!(flr.npv > 0.0);
        let zciis = ZeroCouponInflationSwap {
            notional: 100.0,
            fixed_rate: 0.02,
            maturity_years: 5.0,
            target_month: None,
            pay_fixed: true,
        };
        assert!((cap.npv - flr.npv - zciis.npv(&disc, &infl)).abs() < 1e-10);

        let eval_sig = |sc: f64| {
            let mut p = jy;
            p.sigma_cpi = sc;
            price_zciis_capfloor_jarrow_yildirim(&cap_spec, &disc, &infl, &p).npv
        };
        let fd = (eval_sig(0.012 + h) - eval_sig(0.012 - h)) / (2.0 * h);
        type D = Dual<1>;
        let disc_d = FlatForwardCurve::new(d, D::from_value(0.03));
        let infl_d = infl.to_ad(&[
            D::from_value(0.022),
            D::from_value(0.024),
            D::from_value(0.025),
        ]);
        let jy_d = JarrowYildirimParams {
            a_nom: D::from_value(0.05),
            sigma_nom: D::from_value(0.008),
            a_real: D::from_value(0.04),
            sigma_real: D::from_value(0.006),
            sigma_cpi: D::variable(0.012, 0),
            rho_nom_real: D::from_value(0.3),
            rho_nom_cpi: D::from_value(0.2),
            rho_real_cpi: D::from_value(-0.1),
        };
        let cap_d = price_zciis_capfloor_jarrow_yildirim(
            &ZciisCapFloorSpec {
                notional: D::from_value(100.0),
                maturity_years: 5.0,
                strike_annual_rate: D::from_value(0.02),
                target_month: None,
                is_cap: true,
            },
            &disc_d,
            &infl_d,
            &jy_d,
        );
        assert!(cap_d.npv.d[0] > 0.0);
        assert!((cap_d.npv.d[0] - fd).abs() < 2e-6);
    }
}

#[cfg(test)]
mod phase22_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn brigo_mercurio_dual_fd_parity_and_lpi_bounds() {
        let make = |k| BrigoMercurioInflationMarketModelSpec {
            notional: 1_000_000.0,
            base_cpi: 100.0,
            payment_times: vec![1.0, 2.0, 3.0],
            forward_cpis: vec![102.0, 104.0, 106.0],
            nominal_forward_rates: vec![0.03; 3],
            nominal_discount_factors: vec![0.97, 0.94, 0.91],
            forward_cpi_vols: vec![0.01; 3],
            nominal_forward_vols: vec![0.01; 3],
            inter_maturity_cpi_corr: 0.5,
            nominal_cpi_corr: 0.2,
            lpi_floor: -0.02,
            lpi_cap: 0.03,
            fixed_swap_rate: k,
        };
        let r = price_brigo_mercurio_yoy_and_lpi_market_model(&make(0.02));
        assert!(
            r.fimm_lpi_expected_rates
                .iter()
                .all(|x| x.val() >= -0.02 - 1e-8 && x.val() <= 0.03 + 1e-8)
        );
        let h = 1e-5;
        let fd = (price_brigo_mercurio_yoy_and_lpi_market_model(&make(0.02 + h))
            .fimm_convexity_adjusted_lpi_swap_npv
            - price_brigo_mercurio_yoy_and_lpi_market_model(&make(0.02 - h))
                .fimm_convexity_adjusted_lpi_swap_npv)
            / (2.0 * h);
        let ad =
            price_brigo_mercurio_yoy_and_lpi_market_model(&BrigoMercurioInflationMarketModelSpec {
                notional: Dual::<1>::from_value(1_000_000.0),
                base_cpi: Dual::<1>::from_value(100.0),
                payment_times: vec![
                    Dual::<1>::from_value(1.0),
                    Dual::<1>::from_value(2.0),
                    Dual::<1>::from_value(3.0),
                ],
                forward_cpis: vec![
                    Dual::<1>::from_value(102.0),
                    Dual::<1>::from_value(104.0),
                    Dual::<1>::from_value(106.0),
                ],
                nominal_forward_rates: vec![Dual::<1>::from_value(0.03); 3],
                nominal_discount_factors: vec![
                    Dual::<1>::from_value(0.97),
                    Dual::<1>::from_value(0.94),
                    Dual::<1>::from_value(0.91),
                ],
                forward_cpi_vols: vec![Dual::<1>::from_value(0.01); 3],
                nominal_forward_vols: vec![Dual::<1>::from_value(0.01); 3],
                inter_maturity_cpi_corr: Dual::<1>::from_value(0.5),
                nominal_cpi_corr: Dual::<1>::from_value(0.2),
                lpi_floor: Dual::<1>::from_value(-0.02),
                lpi_cap: Dual::<1>::from_value(0.03),
                fixed_swap_rate: Dual::<1>::variable(0.02, 0),
            });
        assert!((ad.fimm_convexity_adjusted_lpi_swap_npv.d[0] - fd).abs() < 1e-5);
    }
}
