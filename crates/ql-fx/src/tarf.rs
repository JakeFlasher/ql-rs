use ql_ad::Real;
use ql_mc::PhiloxRng;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum FxTarfCapConvention {
    ExactCap,
    FullCouponAtKnockOut,
}
#[derive(Clone, Debug)]
pub struct FxTarfSpec<T: Real> {
    pub notional: T,
    pub spot: T,
    pub strike: T,
    pub r_dom: T,
    pub r_for: T,
    pub vol: T,
    pub fixing_times: Vec<f64>,
    pub target_profit: T,
    pub downside_leverage: T,
    pub is_exporter_call: bool,
    pub cap_convention: FxTarfCapConvention,
    pub smoothing_eps: f64,
}
#[derive(Clone, Debug)]
pub struct FxTarfValuation<T: Real> {
    pub npv: T,
    pub expected_knockout_time: T,
    pub expected_accumulated_profit: T,
}

pub fn price_fx_tarf_mc_ad<T: Real>(
    spec: &FxTarfSpec<T>,
    n_paths: usize,
    seed: u64,
) -> FxTarfValuation<T> {
    if spec.fixing_times.is_empty() || n_paths == 0 {
        return FxTarfValuation {
            npv: T::zero(),
            expected_knockout_time: T::zero(),
            expected_accumulated_profit: T::zero(),
        };
    }
    let rng = PhiloxRng::new(seed);
    let eps_s = spec.smoothing_eps.max(1e-6);
    let eps_t = spec.smoothing_eps.max(1e-6);
    let eta = T::from_f64(if spec.is_exporter_call { 1.0 } else { -1.0 });
    let mut pv = T::zero();
    let mut kt = T::zero();
    let mut ap = T::zero();
    for p in 0..n_paths {
        let sign = if p % 2 == 0 { 1.0 } else { -1.0 };
        let pair = (p / 2) as u64;
        let mut spot = spec.spot;
        let mut acc = T::zero();
        let mut path_pv = T::zero();
        let mut knocked = T::zero();
        for (m, &time) in spec.fixing_times.iter().enumerate() {
            let prev_t = if m == 0 {
                0.0
            } else {
                spec.fixing_times[m - 1]
            };
            let dt = (time - prev_t).max(0.0);
            let z = rng.normal4(pair, m as u32)[0] * sign;
            spot *= (((spec.r_dom - spec.r_for - spec.vol * spec.vol * 0.5) * T::from_f64(dt))
                + spec.vol * T::from_f64(dt.sqrt() * z))
            .exp();
            let x = eta * (spot - spec.strike);
            let root = (x * x + T::from_f64(eps_s * eps_s)).sqrt();
            let gain = (x + root - T::from_f64(eps_s)).max_val(T::zero()) * 0.5;
            let loss = (-x + root - T::from_f64(eps_s)).max_val(T::zero()) * 0.5;
            let room = spec.target_profit - acc;
            let room_root = (room * room + T::from_f64(eps_t * eps_t)).sqrt();
            let alive = (T::one() + room / room_root) * 0.5;
            let positive_room = (room + room_root - T::from_f64(eps_t)).max_val(T::zero()) * 0.5;
            let capped = match spec.cap_convention {
                FxTarfCapConvention::ExactCap => {
                    let d = positive_room - gain;
                    let excess = (d + (d * d + T::from_f64(eps_t * eps_t)).sqrt()) * 0.5;
                    alive * (positive_room - excess).max_val(T::zero())
                }
                FxTarfCapConvention::FullCouponAtKnockOut => alive * gain,
            };
            let cf = spec.notional
                * (-spec.r_dom * T::from_f64(time)).exp()
                * (capped - alive * spec.downside_leverage * loss);
            path_pv += cf;
            acc += capped;
            knocked += alive * T::from_f64(dt);
        }
        pv += path_pv;
        kt += knocked;
        ap += acc;
    }
    let inv = T::from_f64(1.0 / n_paths as f64);
    FxTarfValuation {
        npv: pv * inv,
        expected_knockout_time: kt * inv,
        expected_accumulated_profit: ap * inv,
    }
}

#[cfg(test)]
mod phase14_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn phase14_tarf_cap_and_delta_match_fd() {
        let eval = |spot: f64, target: f64| {
            let s = FxTarfSpec {
                notional: 1.0,
                spot,
                strike: 1.0,
                r_dom: 0.02,
                r_for: 0.01,
                vol: 0.12,
                fixing_times: vec![0.25, 0.5, 0.75, 1.0],
                target_profit: target,
                downside_leverage: 1.5,
                is_exporter_call: true,
                cap_convention: FxTarfCapConvention::ExactCap,
                smoothing_eps: 1e-3,
            };
            price_fx_tarf_mc_ad(&s, 64, 9)
        };
        let v_cap = eval(1.05, 0.10);
        assert!(v_cap.expected_accumulated_profit <= 0.10 + 1e-3);
        let h = 1e-5;
        let fd_spot = (eval(1.05 + h, 0.15).npv - eval(1.05 - h, 0.15).npv) / (2.0 * h);
        let fd_target = (eval(1.05, 0.15 + h).npv - eval(1.05, 0.15 - h).npv) / (2.0 * h);
        let sd = FxTarfSpec {
            notional: Dual::<2>::from_value(1.0),
            spot: Dual::<2>::variable(1.05, 0),
            strike: Dual::<2>::from_value(1.0),
            r_dom: Dual::<2>::from_value(0.02),
            r_for: Dual::<2>::from_value(0.01),
            vol: Dual::<2>::from_value(0.12),
            fixing_times: vec![0.25, 0.5, 0.75, 1.0],
            target_profit: Dual::<2>::variable(0.15, 1),
            downside_leverage: Dual::<2>::from_value(1.5),
            is_exporter_call: true,
            cap_convention: FxTarfCapConvention::ExactCap,
            smoothing_eps: 1e-3,
        };
        let ad = price_fx_tarf_mc_ad(&sd, 64, 9).npv;
        assert!((ad.d[0] - fd_spot).abs() < 1e-6);
        assert!((ad.d[1] - fd_target).abs() < 1e-6);
    }
}

#[derive(Clone, Debug)]
pub struct FxAccumulatorSpec<T: Real = f64> {
    pub base_notional_per_fixing: T,
    pub spot: T,
    pub strike: T,
    pub knockout_barrier: T,
    pub r_dom: T,
    pub r_for: T,
    pub vol: T,
    pub fixing_times: Vec<f64>,
    pub guaranteed_fixings: usize,
    pub downside_gearing: T,
    pub is_accumulator: bool,
    pub smoothing_eps: f64,
}

#[derive(Copy, Clone, Debug)]
pub struct FxAccumulatorValuation<T: Real> {
    pub npv: T,
    pub expected_active_fixings: T,
    pub no_knockout_strip_npv: T,
}

pub fn price_fx_accumulator_strip_no_ko<T: Real>(spec: &FxAccumulatorSpec<T>) -> T {
    let mut sum = T::zero();
    for &t in &spec.fixing_times {
        if t <= 0.0 {
            continue;
        }
        let t_t = T::from_f64(t);
        let c = crate::garman_kohlhagen_price(
            spec.spot,
            spec.strike,
            spec.vol,
            t_t,
            spec.r_dom,
            spec.r_for,
            true,
        );
        let p = crate::garman_kohlhagen_price(
            spec.spot,
            spec.strike,
            spec.vol,
            t_t,
            spec.r_dom,
            spec.r_for,
            false,
        );
        let unit = if spec.is_accumulator {
            c - spec.downside_gearing * p
        } else {
            p - spec.downside_gearing * c
        };
        sum += spec.base_notional_per_fixing * unit;
    }
    sum
}

pub fn price_fx_accumulator_mc_ad<T: Real>(
    spec: &FxAccumulatorSpec<T>,
    n_paths: usize,
    seed: u64,
) -> FxAccumulatorValuation<T> {
    let no_ko = price_fx_accumulator_strip_no_ko(spec);
    if spec.fixing_times.is_empty() || n_paths == 0 {
        return FxAccumulatorValuation {
            npv: T::zero(),
            expected_active_fixings: T::zero(),
            no_knockout_strip_npv: no_ko,
        };
    }
    let paths = n_paths.max(2);
    let rng = PhiloxRng::new(seed);
    let eps = T::from_f64(spec.smoothing_eps.max(1e-4));
    let omega = if spec.is_accumulator { 1.0 } else { -1.0 };
    let mut sum_pv = T::zero();
    let mut sum_active = T::zero();

    for p in 0..paths {
        let sign = if p % 2 == 0 { 1.0 } else { -1.0 };
        let pair = (p / 2) as u64;
        let mut s = spec.spot;
        let mut alive = T::one();
        let mut prev_t = 0.0;
        let mut path_pv = T::zero();
        let mut path_active = T::zero();

        for (m, &t) in spec.fixing_times.iter().enumerate() {
            let dt = (t - prev_t).max(1e-8);
            let z = rng.normal4(pair, m as u32)[0] * sign;
            s *= ((spec.r_dom - spec.r_for - spec.vol * spec.vol * 0.5) * dt
                + spec.vol * (dt.sqrt() * z))
                .exp();

            let bar_dist = (spec.knockout_barrier - s) * omega;
            let step_surv =
                T::from_f64(0.5) * (T::one() + bar_dist / (bar_dist * bar_dist + eps * eps).sqrt());
            alive *= step_surv;
            let w_active = if m < spec.guaranteed_fixings {
                T::one()
            } else {
                alive
            };

            let x = (s - spec.strike) * omega;
            let root = (x * x + T::from_f64(1e-12)).sqrt();
            let gain = (x + root) * 0.5;
            let loss = (-x + root) * 0.5;
            let df = (-spec.r_dom * t).exp();
            path_pv += spec.base_notional_per_fixing
                * df
                * w_active
                * (gain - spec.downside_gearing * loss);
            path_active += w_active;
            prev_t = t;
        }
        sum_pv += path_pv;
        sum_active += path_active;
    }

    let inv = 1.0 / (paths as f64);
    FxAccumulatorValuation {
        npv: sum_pv * inv,
        expected_active_fixings: sum_active * inv,
        no_knockout_strip_npv: no_ko,
    }
}

#[cfg(test)]
mod phase17_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn fx_accumulator_forward_strip_limit_and_ad() {
        let s = FxAccumulatorSpec {
            base_notional_per_fixing: 10.0,
            spot: 1.10,
            strike: 1.05,
            knockout_barrier: 1.20,
            r_dom: 0.03,
            r_for: 0.01,
            vol: 0.12,
            fixing_times: vec![0.25, 0.50, 0.75, 1.0],
            guaranteed_fixings: 1,
            downside_gearing: 1.0,
            is_accumulator: true,
            smoothing_eps: 1e-2,
        };
        let strip = price_fx_accumulator_strip_no_ko(&s);
        let fwd_strip: f64 = s
            .fixing_times
            .iter()
            .map(|&t| 10.0 * (1.10 * (-0.01 * t).exp() - 1.05 * (-0.03 * t).exp()))
            .sum();
        assert!((strip - fwd_strip).abs() < 1e-12);

        let mut s_geared = s.clone();
        s_geared.downside_gearing = 2.0;
        assert!(price_fx_accumulator_strip_no_ko(&s_geared) < strip);

        let h = 1e-5;
        let eval_spot = |sp_fx: f64| {
            let mut sp = s_geared.clone();
            sp.spot = sp_fx;
            price_fx_accumulator_mc_ad(&sp, 64, 17).npv
        };
        let fd = (eval_spot(1.10 + h) - eval_spot(1.10 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sd = FxAccumulatorSpec {
            base_notional_per_fixing: D::from_value(10.0),
            spot: D::variable(1.10, 0),
            strike: D::from_value(1.05),
            knockout_barrier: D::from_value(1.20),
            r_dom: D::from_value(0.03),
            r_for: D::from_value(0.01),
            vol: D::from_value(0.12),
            fixing_times: vec![0.25, 0.50, 0.75, 1.0],
            guaranteed_fixings: 1,
            downside_gearing: D::from_value(2.0),
            is_accumulator: true,
            smoothing_eps: 1e-2,
        };
        let ad = price_fx_accumulator_mc_ad(&sd, 64, 17).npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
