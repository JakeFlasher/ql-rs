use ql_ad::Real;
use ql_curves::YieldTermStructure;

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
        return T::one();
    }
    let b = model.b_factor(t, big_t);
    let p0t = curve.discount_factor(t);
    let p0t_big = curve.discount_factor(big_t);
    let bt0 = model.b_factor(0.0, t);
    let vr = hw_var(model, T::from_f64(t));
    p0t_big / p0t
        * (-T::from_f64(0.5) * model.sigma * model.sigma * b * (b * vr + bt0 * bt0) - b * x).exp()
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
    let mut ann = T::zero();
    let mut float = T::zero();
    for w in s.swap_times[start..].windows(2) {
        let a = conditional_bond(model, ois, t, w[0], x);
        let b = conditional_bond(model, ois, t, w[1], x);
        let af = conditional_bond(model, fwd, t, w[0], x);
        let bf = conditional_bond(model, fwd, t, w[1], x);
        let tau = w[1] - w[0];
        ann += T::from_f64(tau) * b;
        float += (af / bf - T::one()) * b;
        let _ = a;
    }
    let sign = if s.is_payer { 1.0 } else { -1.0 };
    let value = match s.settlement {
        SwaptionSettlement::Physical => float - T::from_f64(s.strike) * ann,
        SwaptionSettlement::CashParYield => {
            let rate = float / ann.max_val(T::from_f64(1e-14));
            T::from_f64(sign) * ann * (rate - s.strike)
        }
    };
    T::from_f64(s.notional) * T::from_f64(sign) * value.max_val(T::zero())
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
        let n = self.steps.max(1);
        let max_t = *ex.last().unwrap();
        let dt = max_t / n as f64;
        let mut times = vec![0.0];
        for i in 1..=n {
            times.push(i as f64 * dt);
        }
        for &e in &ex {
            if !times.iter().any(|x| (*x - e).abs() < 1e-10) {
                times.push(e);
            }
        }
        times.sort_by(f64::total_cmp);
        times.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        let max_dt = times.windows(2).map(|w| w[1] - w[0]).fold(0.0, f64::max);
        let dx = (3.0 * model.sigma.val() * model.sigma.val() * max_dt)
            .sqrt()
            .max(1e-8);
        let mut values = vec![T::zero(); 2 * n + 3];
        let last = *times.last().unwrap();
        for j in -(n as isize)..=(n as isize) {
            values[(j + n as isize) as usize] =
                bermudan_payoff(model, ois, fwd, s, last, T::from_f64(j as f64 * dx));
        }
        for m in (0..times.len() - 1).rev() {
            let dtime = times[m + 1] - times[m];
            let var = model.sigma * model.sigma * hw_var(model, T::from_f64(dtime));
            let ratio = (var.val() / dx.powi(2)).max(0.0);
            let pm = (1.0 - ratio).clamp(0.0, 1.0);
            let pu = (ratio * 0.5).clamp(0.0, 0.5);
            let pd = pu;
            let mut next = vec![T::zero(); values.len()];
            for j in -(n as isize)..=(n as isize) {
                let idx = (j + n as isize) as usize;
                let cont = (-T::from_f64(j as f64 * dx) * dtime).exp()
                    * (values[idx.saturating_sub(1)] * pd
                        + values[idx] * pm
                        + values[(idx + 1).min(values.len() - 1)] * pu);
                next[idx] = cont;
            }
            for &e in &ex {
                if (times[m] - e).abs() < 1e-10 {
                    for j in -(n as isize)..=(n as isize) {
                        let idx = (j + n as isize) as usize;
                        next[idx] = next[idx].max_val(bermudan_payoff(
                            model,
                            ois,
                            fwd,
                            s,
                            times[m],
                            T::from_f64(j as f64 * dx),
                        ));
                    }
                }
            }
            values = next;
        }
        values[n]
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
        m: &ql_models::HullWhite1F<T>,
        ois: &D,
        fwd: &F,
        s: &BermudanSwaption,
    ) -> T {
        HullWhiteBermudanTreeEngine {
            steps: self.time_steps.max(1),
        }
        .price(m, ois, fwd, s)
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
