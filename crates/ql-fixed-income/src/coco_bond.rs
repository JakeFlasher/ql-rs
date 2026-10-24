use ql_ad::Real;

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
    if t <= 0.0 || s.spot.val() <= s.trigger_barrier.val() {
        return T::one();
    }
    if s.eq_vol.val() <= 1e-14 {
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
    let b = (T::from_f64(2.0) * drift * ell / (s.eq_vol * s.eq_vol)).exp()
        * ((ell + drift * t) / sig).norm_cdf();
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
