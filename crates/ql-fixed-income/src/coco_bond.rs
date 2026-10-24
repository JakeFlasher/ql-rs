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
