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
