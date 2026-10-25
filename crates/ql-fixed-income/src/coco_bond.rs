use ql_ad::Real;
use ql_vol::black_scholes_price;

#[derive(Copy, Clone, Debug)]
pub enum CocoTriggerPayoff<T: Real> {
    WriteDown { recovery_fraction: T },
    EquityConversion { conversion_price: T },
}
#[derive(Clone, Debug)]
pub struct CocoBondSpec<T: Real> {
    pub face_value: T,
    pub coupon_rate: T,
    pub coupon_times: Vec<f64>,
    pub maturity: f64,
    pub rate: T,
    pub credit_spread: T,
    pub spot: T,
    pub trigger_barrier: T,
    pub div_yield: T,
    pub eq_vol: T,
    pub trigger_payoff: CocoTriggerPayoff<T>,
}
#[derive(Copy, Clone, Debug)]
pub struct CocoBondValuation<T: Real> {
    pub dirty_price: T,
    pub straight_bond_floor: T,
    pub coupon_loss_pv: T,
    pub principal_loss_pv: T,
    pub conversion_or_recovery_pv: T,
    pub trigger_probability: T,
}
fn hit_probability<T: Real>(s: &CocoBondSpec<T>, t: f64, stock_weight: bool) -> T {
    if s.trigger_barrier.val() <= 0.0 {
        return T::zero();
    }
    if s.spot.val() <= s.trigger_barrier.val() {
        return T::one();
    }
    if t <= 0.0 || s.eq_vol.val() <= 1e-14 {
        return T::zero();
    }
    let sig = s.eq_vol * T::from_f64(t.sqrt());
    let nu = s.rate - s.div_yield - s.eq_vol * s.eq_vol * 0.5;
    let drift = if stock_weight {
        nu + s.eq_vol * s.eq_vol
    } else {
        nu
    };
    let ell = (s.trigger_barrier / s.spot).ln();
    let a = ((ell - drift * t) / sig).norm_cdf();
    let exp_arg = (T::from_f64(2.0) * drift * ell / (s.eq_vol * s.eq_vol))
        .max_val(T::from_f64(-80.0))
        .min_val(T::from_f64(80.0));
    let b = exp_arg.exp() * ((ell + drift * t) / sig).norm_cdf();
    (a + b).max_val(T::zero()).min_val(T::one())
}
pub fn price_coco_bond_ds<T: Real>(s: &CocoBondSpec<T>) -> CocoBondValuation<T> {
    if s.maturity <= 0.0 || s.face_value.val() <= 0.0 {
        return CocoBondValuation {
            dirty_price: T::zero(),
            straight_bond_floor: T::zero(),
            coupon_loss_pv: T::zero(),
            principal_loss_pv: T::zero(),
            conversion_or_recovery_pv: T::zero(),
            trigger_probability: T::zero(),
        };
    }
    let lambda = s.credit_spread.max_val(T::zero());
    let mut floor = s.face_value * (-(s.rate + lambda) * s.maturity).exp();
    let mut coupon_loss = T::zero();
    let mut dt_prev = 0.0;
    for &t in &s.coupon_times {
        let tau = (t - dt_prev).max(0.0);
        let c = s.face_value * s.coupon_rate * tau;
        let df = (-(s.rate + lambda) * t).exp();
        floor += c * df;
        coupon_loss += c * df * hit_probability(s, t, false);
        dt_prev = t;
    }
    let df_m = (-(s.rate + lambda) * s.maturity).exp();
    let ph = hit_probability(s, s.maturity, false);
    let principal_loss = s.face_value * df_m * ph;
    let conversion = match s.trigger_payoff {
        CocoTriggerPayoff::WriteDown { recovery_fraction } => {
            s.face_value * recovery_fraction * df_m * ph
        }
        CocoTriggerPayoff::EquityConversion { conversion_price } => {
            s.face_value / conversion_price
                * s.spot
                * (-s.div_yield * s.maturity).exp()
                * hit_probability(s, s.maturity, true)
        }
    };
    CocoBondValuation {
        dirty_price: floor - coupon_loss - principal_loss + conversion,
        straight_bond_floor: floor,
        coupon_loss_pv: coupon_loss,
        principal_loss_pv: principal_loss,
        conversion_or_recovery_pv: conversion,
        trigger_probability: ph,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct DeSpiegeleerCoCoMandatoryConvSpec<T: Real = f64> {
    pub notional: T,
    pub spot: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub credit_spread: T,
    pub maturity: T,
    pub coupon_rate: T,
    pub coupon_freq: usize,
    pub decs_lower_strike: T,
    pub decs_upper_strike: T,
    pub vol_lower_strike: T,
    pub vol_upper_strike: T,
    pub coco_cet1_trigger_spot: T,
    pub coco_ponv_trigger_spot: T,
    pub coco_conversion_floor_price: T,
    pub coco_equity_vol: T,
}
#[derive(Clone, Debug)]
pub struct DeSpiegeleerCoCoMandatoryConvReport<T: Real = f64> {
    pub decs_mandatory_convertible_value: T,
    pub percs_capped_convertible_value: T,
    pub decs_embedded_collar_value: T,
    pub decs_equity_delta: T,
    pub dual_trigger_coco_value: T,
    pub coco_cet1_trigger_probability: T,
    pub coco_ponv_writedown_probability: T,
    pub coco_implied_yield_spread_bps: T,
}
fn hit<T: Real>(s: T, b: T, mu: T, sig: T, t: T) -> T {
    if b.val() >= s.val() {
        return T::one();
    }
    if b.val() <= 0. || sig.val() <= 1e-12 || t.val() <= 0. {
        return T::zero();
    }
    let st = sig * t.sqrt();
    let l = (b / s).ln();
    let a = ((l - mu * t) / st).norm_cdf();
    let e = (T::from_f64(2.0) * mu * l / (sig * sig)).exp();
    (a + e * ((l + mu * t) / st).norm_cdf())
        .max_val(T::zero())
        .min_val(T::one())
}
pub fn price_despiegeleer_coco_and_mandatory_convertible<T: Real>(
    s: &DeSpiegeleerCoCoMandatoryConvSpec<T>,
) -> DeSpiegeleerCoCoMandatoryConvReport<T> {
    let n = (s.coupon_freq.max(1) as f64 * s.maturity.val().round().max(1.0)) as usize;
    let dt = s.maturity / T::from_f64(n as f64);
    let rr = s.risk_free_rate + s.credit_spread;
    let mut coupons = T::zero();
    for j in 1..=n {
        coupons += s.notional * s.coupon_rate * dt * (-rr * dt * j as f64).exp();
    }
    let cl = black_scholes_price(
        s.spot,
        s.decs_lower_strike,
        s.vol_lower_strike,
        s.maturity,
        s.risk_free_rate,
        s.dividend_yield,
        true,
    );
    let cu = black_scholes_price(
        s.spot,
        s.decs_upper_strike,
        s.vol_upper_strike,
        s.maturity,
        s.risk_free_rate,
        s.dividend_yield,
        true,
    );
    let ru = s.notional / s.decs_lower_strike;
    let rl = s.notional / s.decs_upper_strike;
    let decs = coupons + ru * s.spot * (-s.dividend_yield * s.maturity).exp() - ru * cl + rl * cu;
    let percs = coupons + ru * s.spot * (-s.dividend_yield * s.maturity).exp() - ru * cu;
    let d1 = |k: T, vol: T| {
        if vol.val() <= 1e-12 || s.maturity.val() <= 1e-12 {
            if s.spot.val() > k.val() {
                T::from_f64(1e6)
            } else {
                T::from_f64(-1e6)
            }
        } else {
            ((s.spot / k).ln()
                + (s.risk_free_rate - s.dividend_yield + vol * vol * 0.5) * s.maturity)
                / (vol * s.maturity.sqrt())
        }
    };
    let d1l = d1(s.decs_lower_strike, s.vol_lower_strike);
    let d1u = d1(s.decs_upper_strike, s.vol_upper_strike);
    let delta = (-s.dividend_yield * s.maturity).exp()
        * (ru * (T::one() - d1l.norm_cdf()) + rl * d1u.norm_cdf());
    let straight = coupons + s.notional * (-rr * s.maturity).exp();
    let mut loss = T::zero();
    for j in 1..=n {
        let t = dt * j as f64;
        loss += s.notional
            * s.coupon_rate
            * dt
            * (-rr * t).exp()
            * hit(
                s.spot,
                s.coco_cet1_trigger_spot,
                s.risk_free_rate - s.dividend_yield - s.coco_equity_vol * s.coco_equity_vol * 0.5,
                s.coco_equity_vol,
                t,
            );
    }
    let mu = s.risk_free_rate - s.dividend_yield - s.coco_equity_vol * s.coco_equity_vol * 0.5;
    let mus = mu + s.coco_equity_vol * s.coco_equity_vol;
    let pc = hit(
        s.spot,
        s.coco_cet1_trigger_spot,
        mu,
        s.coco_equity_vol,
        s.maturity,
    );
    let pp = hit(
        s.spot,
        s.coco_ponv_trigger_spot,
        mu,
        s.coco_equity_vol,
        s.maturity,
    );
    let pcs = hit(
        s.spot,
        s.coco_cet1_trigger_spot,
        mus,
        s.coco_equity_vol,
        s.maturity,
    );
    let pps = hit(
        s.spot,
        s.coco_ponv_trigger_spot,
        mus,
        s.coco_equity_vol,
        s.maturity,
    );
    let cp = s
        .coco_cet1_trigger_spot
        .max_val(s.coco_conversion_floor_price);
    let principal = s.notional * (-rr * s.maturity).exp() * pc
        - s.notional / cp
            * s.spot
            * (-(s.dividend_yield + s.credit_spread) * s.maturity).exp()
            * (pcs - pps);
    let coco = straight - loss - principal;
    let spread = if coco.val() > 0. && straight.val() > 0. {
        T::from_f64(10000.0) * (s.credit_spread - (coco / straight).ln() / s.maturity)
    } else {
        T::zero()
    };
    DeSpiegeleerCoCoMandatoryConvReport {
        decs_mandatory_convertible_value: decs,
        percs_capped_convertible_value: percs,
        decs_embedded_collar_value: -ru * cl + rl * cu,
        decs_equity_delta: delta,
        dual_trigger_coco_value: coco,
        coco_cet1_trigger_probability: pc,
        coco_ponv_writedown_probability: pp,
        coco_implied_yield_spread_bps: spread,
    }
}

#[cfg(test)]
mod phase14_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn phase14_coco_barrier_and_spot_delta() {
        let s = CocoBondSpec {
            face_value: 100.0,
            coupon_rate: 0.05,
            coupon_times: vec![1.0, 2.0],
            maturity: 2.0,
            rate: 0.03,
            credit_spread: 0.01,
            spot: 100.0,
            trigger_barrier: 1e-6,
            div_yield: 0.0,
            eq_vol: 0.25,
            trigger_payoff: CocoTriggerPayoff::WriteDown {
                recovery_fraction: 0.4,
            },
        };
        let v = price_coco_bond_ds(&s);
        assert!((v.dirty_price - v.straight_bond_floor).abs() < 1e-2);
        let sd = CocoBondSpec {
            spot: Dual::<1>::variable(100.0, 0),
            face_value: Dual::from_value(100.0),
            coupon_rate: Dual::from_value(0.05),
            coupon_times: vec![1.0, 2.0],
            maturity: 2.0,
            rate: Dual::from_value(0.03),
            credit_spread: Dual::from_value(0.01),
            trigger_barrier: Dual::from_value(60.0),
            div_yield: Dual::from_value(0.0),
            eq_vol: Dual::from_value(0.25),
            trigger_payoff: CocoTriggerPayoff::WriteDown {
                recovery_fraction: Dual::from_value(0.4),
            },
        };
        assert!(price_coco_bond_ds(&sd).dirty_price.d[0] > 0.0);
    }
}
