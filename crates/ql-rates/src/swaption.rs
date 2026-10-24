use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_math::solve_tridiagonal;

fn normal_sabr_vol<T: Real>(f: T, k: T, expiry: f64, alpha: T, beta: f64, rho: T, nu: T) -> T {
    let fk = f * k;
    let fb = fk.powf((1.0 - beta) * 0.5);
    let z = nu / alpha * fb * (f - k);
    let ratio = if z.val().abs() < 1e-7 {
        T::one() - rho * z * 0.5 + (T::from_f64(2.0) - rho * rho * 3.0) * z * z / 12.0
    } else {
        let disc = (T::one() - rho * z * 2.0 + z * z).sqrt();
        z / ((disc + z - rho) / (T::one() - rho)).ln()
    };
    let omb = T::from_f64(1.0 - beta);
    let corr = T::one()
        + T::from_f64(expiry)
            * (-T::from_f64(beta * (2.0 - beta)) * alpha * alpha / (fk.pow(omb) * 24.0)
                + rho * T::from_f64(beta) * nu * alpha / (fk.pow(omb * 0.5) * 4.0)
                + (T::from_f64(2.0) - rho * rho * 3.0) * nu * nu / 24.0);
    alpha * fk.powf(beta * 0.5) * ratio * corr
}
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SwaptionSettlement {
    Physical,
    CashParYield,
}
#[derive(Clone, Debug)]
pub struct EuropeanSwaption {
    pub notional: f64,
    pub strike: f64,
    pub expiry: f64,
    pub swap_times: Vec<f64>,
    pub is_payer: bool,
    pub settlement: SwaptionSettlement,
}
impl EuropeanSwaption {
    fn forward_annuity<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        ois: &D,
        fwd: &F,
    ) -> (T, T, T) {
        if self.swap_times.len() < 2 {
            return (T::zero(), T::zero(), T::zero());
        }
        let mut a = T::zero();
        let mut num = T::zero();
        for w in self.swap_times.windows(2) {
            let tau = w[1] - w[0];
            let ff = (fwd.discount_factor(w[0]) / fwd.discount_factor(w[1]) - T::one()) / tau;
            let df = ois.discount_factor(w[1]);
            a += T::from_f64(tau) * df;
            num += T::from_f64(tau) * ff * df;
        }
        let s = num / a;
        let ann = match self.settlement {
            SwaptionSettlement::Physical => a,
            SwaptionSettlement::CashParYield => {
                let mut x = T::zero();
                for (i, w) in self.swap_times.windows(2).enumerate() {
                    let tau = w[1] - w[0];
                    x += ois.discount_factor(self.swap_times[0]) * T::from_f64(tau)
                        / (T::one() + T::from_f64(tau) * s).powf((i + 1) as f64);
                }
                x
            }
        };
        (s, ann, num)
    }
    pub fn price_bachelier<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        ois: &D,
        fwd: &F,
        normal_vol: T,
    ) -> T {
        let (s, a, _) = self.forward_annuity(ois, fwd);
        let st = normal_vol * self.expiry.sqrt();
        let om = if self.is_payer { 1.0 } else { -1.0 };
        if st.val().abs() <= 1e-14 {
            return T::from_f64(self.notional)
                * a
                * (T::from_f64(om) * (s - self.strike)).max_val(T::zero());
        }
        let d = if st.val().abs() > 1e-14 {
            (s - self.strike) / st
        } else {
            T::zero()
        };
        T::from_f64(self.notional)
            * a
            * (T::from_f64(om) * (s - self.strike) * (T::from_f64(om) * d).norm_cdf()
                + st * d.norm_pdf())
    }
    pub fn price_shifted_black76<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        ois: &D,
        fwd: &F,
        vol: T,
        shift: f64,
    ) -> T {
        let (s, a, _) = self.forward_annuity(ois, fwd);
        let st = vol * self.expiry.sqrt();
        let f = s + shift;
        let k = T::from_f64(self.strike + shift);
        if st.val().abs() <= 1e-14 || f.val() <= 0.0 || k.val() <= 0.0 {
            let om = if self.is_payer { 1.0 } else { -1.0 };
            return T::from_f64(self.notional)
                * a
                * (T::from_f64(om) * (s - self.strike)).max_val(T::zero());
        }
        let d1 = ((f / k).ln() + T::from_f64(0.5) * st * st) / st;
        let d2 = d1 - st;
        let om = if self.is_payer { 1.0 } else { -1.0 };
        T::from_f64(self.notional)
            * a
            * (T::from_f64(om)
                * (f * (T::from_f64(om) * d1).norm_cdf() - k * (T::from_f64(om) * d2).norm_cdf()))
    }
    pub fn price_sabr_bachelier<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        ois: &D,
        fwd: &F,
        alpha: T,
        beta: f64,
        rho: T,
        nu: T,
    ) -> T {
        let (s, _, _) = self.forward_annuity(ois, fwd);
        let vol = normal_sabr_vol(
            s,
            T::from_f64(self.strike),
            self.expiry,
            alpha,
            beta,
            rho,
            nu,
        );
        self.price_bachelier(ois, fwd, vol)
    }
}

#[derive(Clone, Debug)]
pub struct BermudanSwaption {
    pub notional: f64,
    pub strike: f64,
    pub exercise_times: Vec<f64>,
    pub swap_times: Vec<f64>,
    pub is_payer: bool,
    pub settlement: SwaptionSettlement,
}
impl BermudanSwaption {
    pub fn from_european(s: &EuropeanSwaption) -> Self {
        Self {
            notional: s.notional,
            strike: s.strike,
            exercise_times: vec![s.expiry],
            swap_times: s.swap_times.clone(),
            is_payer: s.is_payer,
            settlement: s.settlement,
        }
    }
    pub fn from_resolved_swap(
        swap: &ql_instruments::ResolvedVanillaSwap,
        exercise_times: Vec<f64>,
        settlement: SwaptionSettlement,
    ) -> Self {
        let mut times = Vec::with_capacity(swap.fixed_leg.payment_times.len() + 1);
        if let Some(&t0) = swap.fixed_leg.accrual_start_times.first() {
            times.push(t0);
        }
        times.extend(swap.fixed_leg.payment_times.iter().copied());
        let notional = swap.fixed_leg.notionals.first().copied().unwrap_or(0.0);
        let strike = swap
            .fixed_leg
            .rates_or_spreads
            .first()
            .copied()
            .unwrap_or(0.0);
        Self {
            notional,
            strike,
            exercise_times,
            swap_times: times,
            is_payer: matches!(swap.swap_type, ql_instruments::SwapType::Payer),
            settlement,
        }
    }
}

#[derive(Clone, Debug)]
pub struct SwaptionSabrNode<T: Real> {
    pub expiry: f64,
    pub swap_tenor: f64,
    pub forward: T,
    pub alpha: T,
    pub beta: f64,
    pub rho: T,
    pub nu: T,
}
#[derive(Clone, Debug)]
pub struct SwaptionSabrCube<T: Real> {
    pub expiries: Vec<f64>,
    pub swap_tenors: Vec<f64>,
    pub beta: f64,
    pub nodes: Vec<Vec<SwaptionSabrNode<T>>>,
    pub use_normal_vol: bool,
}
impl<T: Real> SwaptionSabrCube<T> {
    fn bracket(xs: &[f64], x: f64) -> (usize, usize, f64) {
        if xs.len() <= 1 {
            return (0, 0, 0.0);
        }
        if x <= xs[0] {
            return (0, 0, 0.0);
        }
        if x >= xs[xs.len() - 1] {
            let n = xs.len() - 1;
            return (n, n, 0.0);
        }
        let i = xs.partition_point(|v| *v <= x).saturating_sub(1);
        (
            i,
            i + 1,
            ((x - xs[i]) / (xs[i + 1] - xs[i])).clamp(0.0, 1.0),
        )
    }
    pub fn interpolate_params(&self, expiry: f64, tenor: f64) -> (T, T, T) {
        let (i0, i1, we) = Self::bracket(&self.expiries, expiry);
        let (j0, j1, ws) = Self::bracket(&self.swap_tenors, tenor);
        let w = |a: usize, b: usize| {
            T::from_f64(if a == 0 { 1.0 - we } else { we } * if b == 0 { 1.0 - ws } else { ws })
        };
        let alpha = self.nodes[i0][j0].alpha * w(0, 0)
            + self.nodes[i0][j1].alpha * w(0, 1)
            + self.nodes[i1][j0].alpha * w(1, 0)
            + self.nodes[i1][j1].alpha * w(1, 1);
        let rho = (self.nodes[i0][j0].rho * w(0, 0)
            + self.nodes[i0][j1].rho * w(0, 1)
            + self.nodes[i1][j0].rho * w(1, 0)
            + self.nodes[i1][j1].rho * w(1, 1))
        .max_val(T::from_f64(-0.999))
        .min_val(T::from_f64(0.999));
        let nu = self.nodes[i0][j0].nu * w(0, 0)
            + self.nodes[i0][j1].nu * w(0, 1)
            + self.nodes[i1][j0].nu * w(1, 0)
            + self.nodes[i1][j1].nu * w(1, 1);
        (
            alpha.max_val(T::from_f64(1e-8)),
            rho,
            nu.max_val(T::from_f64(1e-8)),
        )
    }
    pub fn implied_vol(&self, expiry: f64, tenor: f64, forward: T, strike: T) -> T {
        let (a, r, n) = self.interpolate_params(expiry, tenor);
        ql_vol::AdSabrSmile::new(
            forward,
            T::from_f64(expiry.max(1e-8)),
            a,
            T::from_f64(self.beta),
            r,
            n,
        )
        .vol(strike)
    }
    pub fn alpha_from_atm_vol(
        forward: T,
        expiry: f64,
        atm: T,
        beta: f64,
        rho: T,
        nu: T,
        normal: bool,
    ) -> T {
        let mut a = if normal {
            atm
        } else {
            atm * forward.powf(1.0 - beta)
        };
        for _ in 0..6 {
            let sm = ql_vol::AdSabrSmile::new(
                forward,
                T::from_f64(expiry.max(1e-8)),
                a,
                T::from_f64(beta),
                rho,
                nu,
            );
            let h = 1e-5_f64.max(a.val().abs() * 1e-5);
            let up = ql_vol::AdSabrSmile::new(
                forward,
                T::from_f64(expiry.max(1e-8)),
                a + T::from_f64(h),
                T::from_f64(beta),
                rho,
                nu,
            )
            .vol(forward);
            let dn = ql_vol::AdSabrSmile::new(
                forward,
                T::from_f64(expiry.max(1e-8)),
                a - T::from_f64(h),
                T::from_f64(beta),
                rho,
                nu,
            )
            .vol(forward);
            a -= (sm.vol(forward) - atm) / ((up - dn) / (2.0 * h)).max_val(T::from_f64(1e-10));
            a = a.max_val(T::from_f64(1e-8));
        }
        a
    }
}
#[derive(Clone, Debug)]
pub struct CalibratedSwaptionSabrCube {
    pub cube: SwaptionSabrCube<f64>,
    pub node_calibrations: Vec<Vec<ql_vol::CalibratedSabrSmile>>,
}
impl CalibratedSwaptionSabrCube {
    pub fn calibrate(
        expiries: Vec<f64>,
        tenors: Vec<f64>,
        beta: f64,
        initial: (f64, f64, f64),
        quotes: &[Vec<(f64, Vec<ql_vol::SmileStrikeQuote>)>],
    ) -> Result<Self, ql_math::MathError> {
        let mut nodes = Vec::new();
        let mut cal = Vec::new();
        for (i, &e) in expiries.iter().enumerate() {
            let mut nr = Vec::new();
            let mut nc = Vec::new();
            for (j, &tenor) in tenors.iter().enumerate() {
                let (f, q) = quotes
                    .get(i)
                    .and_then(|r| r.get(j))
                    .cloned()
                    .ok_or(ql_math::MathError::NotBracketed)?;
                let c = ql_vol::CalibratedSabrSmile::calibrate(f, e, beta, initial, &q)?;
                nr.push(SwaptionSabrNode {
                    expiry: e,
                    swap_tenor: tenor,
                    forward: f,
                    alpha: c.smile.alpha,
                    beta,
                    rho: c.smile.rho,
                    nu: c.smile.nu,
                });
                nc.push(c);
            }
            nodes.push(nr);
            cal.push(nc);
        }
        Ok(Self {
            cube: SwaptionSabrCube {
                expiries,
                swap_tenors: tenors,
                beta,
                nodes,
                use_normal_vol: false,
            },
            node_calibrations: cal,
        })
    }
    pub fn vol_and_quote_vega(
        &self,
        expiry: f64,
        tenor: f64,
        forward: f64,
        strike: f64,
    ) -> (f64, Vec<f64>) {
        let (i0, i1, we) = SwaptionSabrCube::<f64>::bracket(&self.cube.expiries, expiry);
        let (j0, j1, ws) = SwaptionSabrCube::<f64>::bracket(&self.cube.swap_tenors, tenor);
        let n = self
            .node_calibrations
            .iter()
            .flat_map(|r| r.iter().map(|c| c.quotes.len()))
            .sum();
        let mut g = vec![0.0; n];
        let mut off = 0;
        for (i, row) in self.node_calibrations.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                let w = if i == i0 {
                    1.0 - we
                } else if i == i1 {
                    we
                } else {
                    0.0
                } * if j == j0 {
                    1.0 - ws
                } else if j == j1 {
                    ws
                } else {
                    0.0
                };
                if w != 0.0 {
                    let p = ql_vol::AdSabrSmile::new(
                        ql_ad::Dual::<3>::from_value(forward),
                        ql_ad::Dual::<3>::from_value(expiry.max(1e-8)),
                        ql_ad::Dual::<3>::variable(c.smile.alpha, 0),
                        ql_ad::Dual::<3>::from_value(c.smile.beta),
                        ql_ad::Dual::<3>::variable(c.smile.rho, 1),
                        ql_ad::Dual::<3>::variable(c.smile.nu, 2),
                    );
                    let d = p.vol(ql_ad::Dual::<3>::from_value(strike));
                    let v = c.sensitivity.vjp(&[w * d.d[0], w * d.d[1], w * d.d[2]]);
                    for k in 0..v.len() {
                        g[off + k] += v[k];
                    }
                }
                off += c.quotes.len();
            }
        }
        let (aa, rr, nn) = self.cube.interpolate_params(expiry, tenor);
        let pd = ql_vol::AdSabrSmile::new(
            ql_ad::Dual::<3>::from_value(forward),
            ql_ad::Dual::<3>::from_value(expiry.max(1e-8)),
            ql_ad::Dual::from_value(aa),
            ql_ad::Dual::from_value(self.cube.beta),
            ql_ad::Dual::from_value(rr),
            ql_ad::Dual::from_value(nn),
        );
        (pd.vol(ql_ad::Dual::<3>::from_value(strike)).val(), g)
    }
}

#[derive(Clone, Debug)]
pub struct MidCurveSwaptionSpec {
    pub notional: f64,
    pub strike: f64,
    pub option_expiry: f64,
    pub forward_swap_times: Vec<f64>,
    pub is_payer: bool,
}
pub fn price_midcurve_swaption_hw1f<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
    model: &ql_models::HullWhite1F<T>,
    ois: &D,
    fwd: &F,
    s: &MidCurveSwaptionSpec,
) -> T {
    if s.option_expiry <= 0.0
        || s.forward_swap_times.len() < 2
        || s.forward_swap_times[0] < s.option_expiry - 1e-12
    {
        return T::zero();
    }
    let vx = (model.sigma * model.sigma * hw_var(model, T::from_f64(s.option_expiry)))
        .max_val(T::zero())
        .sqrt();
    let mut out = T::zero();
    for (z, w) in [
        (-1.7320508075688772, 1.0 / 6.0),
        (0.0, 2.0 / 3.0),
        (1.7320508075688772, 1.0 / 6.0),
    ] {
        let x = vx * z;
        let mut ann = T::zero();
        let mut fl = T::zero();
        for win in s.forward_swap_times.windows(2) {
            let tau = win[1] - win[0];
            let p1 = conditional_bond(model, ois, s.option_expiry, win[1], x);
            let b0 = conditional_bond(model, fwd, s.option_expiry, win[0], x);
            let b1 = conditional_bond(model, fwd, s.option_expiry, win[1], x);
            ann += p1 * tau;
            fl += (b0 / b1 - T::one()) * p1;
        }
        let sign = T::from_f64(if s.is_payer { 1.0 } else { -1.0 });
        out += T::from_f64(w)
            * T::from_f64(s.notional)
            * (sign * (fl - T::from_f64(s.strike) * ann)).max_val(T::zero());
    }
    out
}

pub fn price_midcurve_swaption_g2pp<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
    model: &crate::G2PlusPlus<T>,
    ois: &D,
    _fwd: &F,
    s: &MidCurveSwaptionSpec,
    n_quad: usize,
) -> T {
    if s.option_expiry <= 0.0
        || s.forward_swap_times.len() < 2
        || s.forward_swap_times[0] < s.option_expiry - 1e-12
    {
        return T::zero();
    }
    let sw = EuropeanSwaption {
        notional: s.notional,
        strike: s.strike,
        expiry: s.option_expiry,
        swap_times: s.forward_swap_times.clone(),
        is_payer: s.is_payer,
        settlement: SwaptionSettlement::Physical,
    };
    model.swaption_price(ois, &sw, n_quad)
}

fn hw_var<T: Real>(model: &ql_models::HullWhite1F<T>, t: T) -> T {
    if (model.a * t).val().abs() < 1e-6 {
        t - model.a * t * t + model.a * model.a * t * t * t * (2.0 / 3.0)
    } else {
        (T::one() - (-model.a * t * 2.0).exp()) / (model.a * 2.0)
    }
}
fn conditional_bond<T: Real, C: YieldTermStructure<T>>(
    model: &ql_models::HullWhite1F<T>,
    curve: &C,
    t: f64,
    big_t: f64,
    x: T,
) -> T {
    if big_t <= t {
        T::one()
    } else {
        let b = model.b_factor(t, big_t);
        let p0t = curve.discount_factor(t);
        let p0t_big = curve.discount_factor(big_t);
        let bt0 = model.b_factor(0.0, t);
        let vr = hw_var(model, T::from_f64(t));
        p0t_big / p0t
            * (-T::from_f64(0.5) * model.sigma * model.sigma * b * (b * vr + bt0 * bt0) - b * x)
                .exp()
    }
}
fn bermudan_payoff<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
    model: &ql_models::HullWhite1F<T>,
    ois: &D,
    fwd: &F,
    s: &BermudanSwaption,
    t: f64,
    x: T,
) -> T {
    let mut start = 0usize;
    while start + 1 < s.swap_times.len() && s.swap_times[start] < t - 1e-10 {
        start += 1;
    }
    if start + 1 >= s.swap_times.len() {
        return T::zero();
    }
    let p_start = conditional_bond(model, ois, t, s.swap_times[start], x);
    let mut ann = T::zero();
    let mut float = T::zero();
    for w in s.swap_times[start..].windows(2) {
        let p_end = conditional_bond(model, ois, t, w[1], x);
        let af = conditional_bond(model, fwd, t, w[0], x);
        let bf = conditional_bond(model, fwd, t, w[1], x);
        let tau = T::from_f64(w[1] - w[0]);
        ann += tau * p_end;
        float += (af / bf - T::one()) * p_end;
    }
    let sign = T::from_f64(if s.is_payer { 1.0 } else { -1.0 });
    let intrinsic = match s.settlement {
        SwaptionSettlement::Physical => float - T::from_f64(s.strike) * ann,
        SwaptionSettlement::CashParYield => {
            let rate = float / ann.max_val(T::from_f64(1e-14));
            let mut cash_ann = T::zero();
            let mut df_par = T::one();
            for w in s.swap_times[start..].windows(2) {
                let tau = T::from_f64(w[1] - w[0]);
                df_par /= (T::one() + tau * rate).max_val(T::from_f64(1e-12));
                cash_ann += p_start * tau * df_par;
            }
            cash_ann * (rate - T::from_f64(s.strike))
        }
    };
    T::from_f64(s.notional) * (sign * intrinsic).max_val(T::zero())
}

#[derive(Copy, Clone, Debug)]
pub struct HullWhiteBermudanTreeEngine {
    pub steps: usize,
}
#[derive(Copy, Clone, Debug)]
pub struct Fd1dHullWhiteBermudanEngine {
    pub time_steps: usize,
    pub state_points: usize,
}

impl HullWhiteBermudanTreeEngine {
    pub fn price<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        model: &ql_models::HullWhite1F<T>,
        ois: &D,
        fwd: &F,
        s: &BermudanSwaption,
    ) -> T {
        let mut ex = s
            .exercise_times
            .iter()
            .copied()
            .filter(|x| *x > 0.0)
            .collect::<Vec<_>>();
        ex.sort_by(f64::total_cmp);
        ex.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        if ex.is_empty() {
            return T::zero();
        }
        let nbase = self.steps.max(1);
        let max_t = *ex.last().unwrap();
        let h = max_t / nbase as f64;
        let mut times = vec![0.0];
        for i in 1..=nbase {
            times.push(i as f64 * h);
        }
        times.extend(ex.iter().copied());
        times.sort_by(f64::total_cmp);
        times.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        let msteps = times.len() - 1;
        let max_dt = times
            .windows(2)
            .map(|w| w[1] - w[0])
            .fold(0.0, f64::max)
            .max(1e-12);
        let v_max = model.sigma * model.sigma * hw_var(model, T::from_f64(max_dt));
        let dx = (v_max * 3.0)
            .max_val(T::from_f64(1e-16))
            .sqrt()
            .max_val(T::from_f64(1e-8));
        let jmax = if model.a.val() * max_dt > 1e-4 {
            (0.184 / (model.a.val() * max_dt)).ceil() as isize
        } else {
            msteps as isize
        };
        let jb = |m: usize| (m as isize).min(jmax);
        let mut q = times
            .iter()
            .enumerate()
            .map(|(m, _)| vec![T::zero(); (2 * jb(m) + 1) as usize])
            .collect::<Vec<_>>();
        q[0][0] = T::one();
        let mut alpha = vec![T::zero(); msteps];
        let trans = |m: usize, j: isize, qj: T, target: &mut Vec<T>, alpha_m: T| {
            let dt = T::from_f64((times[m + 1] - times[m]).max(1e-12));
            let var = model.sigma * model.sigma * hw_var(model, dt);
            let mu = dx * T::from_f64(j as f64) * (-model.a * dt).exp();
            let next_b = jb(m + 1);
            let k = ((mu / dx).val().round() as isize).clamp(-next_b + 1, next_b - 1);
            let e = mu - dx * T::from_f64(k as f64);
            let dx2 = dx * dx;
            let pd = (var + e * e - e * dx) / (dx2 * 2.0);
            let pm = T::one() - (var + e * e) / dx2;
            let pu = (var + e * e + e * dx) / (dx2 * 2.0);
            let disc = (-(alpha_m + dx * T::from_f64(j as f64)) * dt).exp() * qj;
            for (node, p) in [(k - 1, pd), (k, pm), (k + 1, pu)] {
                let idx = (node + next_b).clamp(0, 2 * next_b) as usize;
                target[idx] += disc * p;
            }
        };
        for m in 0..msteps {
            let dt = T::from_f64((times[m + 1] - times[m]).max(1e-12));
            let sum = q[m].iter().enumerate().fold(T::zero(), |acc, (idx, qj)| {
                let j = idx as isize - jb(m);
                acc + *qj * (-dx * T::from_f64(j as f64) * dt).exp()
            });
            alpha[m] = (sum / ois.discount_factor(times[m + 1])).ln() / dt;
            for idx in 0..q[m].len() {
                let j = idx as isize - jb(m);
                trans(m, j, q[m][idx], &mut q[m + 1], alpha[m]);
            }
        }
        let mut values = times
            .iter()
            .enumerate()
            .map(|(m, _)| vec![T::zero(); (2 * jb(m) + 1) as usize])
            .collect::<Vec<_>>();
        let last = msteps;
        for idx in 0..values[last].len() {
            let j = idx as isize - jb(last);
            values[last][idx] =
                bermudan_payoff(model, ois, fwd, s, times[last], dx * T::from_f64(j as f64));
        }
        for m in (0..msteps).rev() {
            let dt = T::from_f64((times[m + 1] - times[m]).max(1e-12));
            let mut next = vec![T::zero(); values[m].len()];
            for (idx, value) in next.iter_mut().enumerate() {
                let j = idx as isize - jb(m);
                let var = model.sigma * model.sigma * hw_var(model, dt);
                let mu = dx * T::from_f64(j as f64) * (-model.a * dt).exp();
                let next_b = jb(m + 1);
                let k = ((mu / dx).val().round() as isize).clamp(-next_b + 1, next_b - 1);
                let e = mu - dx * T::from_f64(k as f64);
                let dx2 = dx * dx;
                let pd = (var + e * e - e * dx) / (dx2 * 2.0);
                let pm = T::one() - (var + e * e) / dx2;
                let pu = (var + e * e + e * dx) / (dx2 * 2.0);
                let get =
                    |node: isize| values[m + 1][(node + next_b).clamp(0, 2 * next_b) as usize];
                *value = (-(alpha[m] + dx * T::from_f64(j as f64)) * dt).exp()
                    * (pd * get(k - 1) + pm * get(k) + pu * get(k + 1));
            }
            if ex.iter().any(|e| (times[m] - *e).abs() < 1e-10) {
                for (idx, value) in next.iter_mut().enumerate() {
                    let j = idx as isize - jb(m);
                    let p =
                        bermudan_payoff(model, ois, fwd, s, times[m], dx * T::from_f64(j as f64));
                    *value = (*value).max_val(p);
                }
            }
            values[m] = next;
        }
        values[0][0]
    }
    pub fn price_single_curve<T: Real, C: YieldTermStructure<T>>(
        &self,
        m: &ql_models::HullWhite1F<T>,
        c: &C,
        s: &BermudanSwaption,
    ) -> T {
        self.price(m, c, c, s)
    }
    pub fn price_from_context<T: Real + 'static>(
        &self,
        m: &ql_models::HullWhite1F<T>,
        ctx: &ql_curves::MarketContext<T>,
        d: ql_curves::DiscountCurveId,
        f: ql_curves::ForwardIndexId,
        s: &BermudanSwaption,
    ) -> Option<T> {
        Some(self.price(m, ctx.discount_curve(d)?, ctx.forward_curve(f)?, s))
    }
}
impl Fd1dHullWhiteBermudanEngine {
    pub fn price<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        model: &ql_models::HullWhite1F<T>,
        ois: &D,
        fwd: &F,
        s: &BermudanSwaption,
    ) -> T {
        let mut ex = s
            .exercise_times
            .iter()
            .copied()
            .filter(|x| *x > 0.0)
            .collect::<Vec<_>>();
        ex.sort_by(f64::total_cmp);
        ex.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        if ex.is_empty() {
            return T::zero();
        }
        let nbase = self.time_steps.max(2);
        let last_t = *ex.last().unwrap();
        let h = last_t / nbase as f64;
        let mut times = vec![0.0];
        for i in 1..=nbase {
            times.push(i as f64 * h);
        }
        times.extend(ex.iter().copied());
        times.sort_by(f64::total_cmp);
        times.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        let nx = self.state_points.max(51) | 1;
        let xmax = (model.sigma * hw_var(model, T::from_f64(last_t)).sqrt() * 5.5)
            .max_val(T::from_f64(1e-4));
        let dx = xmax * 2.0 / T::from_f64((nx - 1) as f64);
        let xgrid = (0..nx)
            .map(|i| -xmax + dx * T::from_f64(i as f64))
            .collect::<Vec<_>>();
        let msteps = times.len() - 1;
        let mut values = vec![T::zero(); nx];
        for (i, x) in xgrid.iter().enumerate() {
            values[i] = bermudan_payoff(model, ois, fwd, s, last_t, *x);
        }
        let nint = nx - 2;
        let mut lower = vec![T::zero(); nint.saturating_sub(1)];
        let mut diag = vec![T::zero(); nint];
        let mut upper = vec![T::zero(); nint.saturating_sub(1)];
        let mut rhs = vec![T::zero(); nint];
        let mut out = vec![T::zero(); nint];
        let mut scratch = vec![T::zero(); nint];
        for m in (0..msteps).rev() {
            let dt_f = (times[m + 1] - times[m]).max(1e-12);
            let dt = T::from_f64(dt_f);
            let theta = if msteps - m <= 2 {
                T::one()
            } else {
                T::from_f64(0.5)
            };
            let phi = -((ois.discount_factor(times[m + 1]) / ois.discount_factor(times[m])).ln())
                / dt
                + model.sigma
                    * model.sigma
                    * (model.b_factor(0.0, times[m]).powf(2.0)
                        + model.b_factor(0.0, times[m + 1]).powf(2.0))
                    * 0.25;
            let dx2 = dx * dx;
            let sig2 = model.sigma * model.sigma;
            let mut new_values = vec![T::zero(); nx];
            let bc_low = bermudan_payoff(model, ois, fwd, s, times[m], xgrid[0]);
            let bc_high = bermudan_payoff(model, ois, fwd, s, times[m], xgrid[nx - 1]);
            let bc_low_next = bermudan_payoff(model, ois, fwd, s, times[m + 1], xgrid[0]);
            let bc_high_next = bermudan_payoff(model, ois, fwd, s, times[m + 1], xgrid[nx - 1]);
            new_values[0] = bc_low;
            new_values[nx - 1] = bc_high;
            for i in 1..nx - 1 {
                let ell = sig2 / (dx2 * 2.0) + model.a * xgrid[i] / (dx * 2.0);
                let upp = sig2 / (dx2 * 2.0) - model.a * xgrid[i] / (dx * 2.0);
                let c = -sig2 / dx2 - (xgrid[i] + phi);
                let ii = i - 1;
                let a_left = -theta * dt * ell;
                let a_right = -theta * dt * upp;
                diag[ii] = T::one() - theta * dt * c;
                if ii > 0 {
                    lower[ii - 1] = a_left;
                }
                if ii + 1 < nint {
                    upper[ii] = a_right;
                }
                rhs[ii] = values[i]
                    + (T::one() - theta)
                        * dt
                        * (ell * values[i - 1] + c * values[i] + upp * values[i + 1]);
                if ii == 0 {
                    rhs[ii] -= a_left * bc_low;
                }
                if ii + 1 == nint {
                    rhs[ii] -= a_right * bc_high;
                }
            }
            solve_tridiagonal(&lower, &diag, &upper, &rhs, &mut out, &mut scratch);
            new_values[1..nx - 1].copy_from_slice(&out);
            if ex.iter().any(|e| (times[m] - *e).abs() < 1e-10) {
                for i in 1..nx - 1 {
                    new_values[i] = new_values[i]
                        .max_val(bermudan_payoff(model, ois, fwd, s, times[m], xgrid[i]));
                }
            }
            // The unused next boundary values are deliberately evaluated in
            // T so their AD dependence is retained when the grid is refined.
            let _ = (bc_low_next, bc_high_next);
            values = new_values;
            if values.iter().any(|z| !z.val().is_finite()) {
                return T::zero();
            }
        }
        values[nx / 2]
    }
    pub fn price_single_curve<T: Real, C: YieldTermStructure<T>>(
        &self,
        m: &ql_models::HullWhite1F<T>,
        c: &C,
        s: &BermudanSwaption,
    ) -> T {
        self.price(m, c, c, s)
    }
    pub fn price_from_context<T: Real + 'static>(
        &self,
        m: &ql_models::HullWhite1F<T>,
        ctx: &ql_curves::MarketContext<T>,
        d: ql_curves::DiscountCurveId,
        f: ql_curves::ForwardIndexId,
        s: &BermudanSwaption,
    ) -> Option<T> {
        Some(self.price(m, ctx.discount_curve(d)?, ctx.forward_curve(f)?, s))
    }
}
