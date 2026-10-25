use ql_ad::{Real, TapeScope, TapeVar};
use ql_curves::YieldTermStructure;
use ql_math::brent_solve;
use ql_time::Date;
use ql_vol::{bachelier_price, black76_price};
use thiserror::Error;

#[derive(Debug, Error, Clone)]
pub enum InflationError {
    #[error("invalid seasonality factors")]
    InvalidSeasonality,
    #[error("invalid inflation curve")]
    InvalidCurve,
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
        let vc = price_yoy_capfloor_jarrow_yildirim(&cap, &disc, &inf, &p).0;
        let vf = price_yoy_capfloor_jarrow_yildirim(&flr, &disc, &inf, &p).0;
        let vcol = price_yoy_capfloor_jarrow_yildirim(&col, &disc, &inf, &p).0;
        assert!((vc - vf - vcol).abs() < 1e-12);
        let h = 1e-5;
        let mut pu = p;
        let mut pd = p;
        pu.sigma_real += h;
        pd.sigma_real -= h;
        let fd = (price_yoy_capfloor_jarrow_yildirim(&cap, &disc, &inf, &pu).0
            - price_yoy_capfloor_jarrow_yildirim(&cap, &disc, &inf, &pd).0)
            / (2.0 * h);
        type D = Dual<1>;
        let dd = FlatForwardCurve::new(d, D::from_value(0.025));
        let id = ZeroInflationCurve::new(
            d,
            100.0,
            vec![1.0, 2.0, 5.0],
            vec![D::from_value(0.02); 3],
            None,
        );
        let cd = YoYInflationCapFloor {
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
        assert!((price_yoy_capfloor_jarrow_yildirim(&cd, &dd, &id, &pd).0.d[0] - fd).abs() < 1e-5);
    }
}
