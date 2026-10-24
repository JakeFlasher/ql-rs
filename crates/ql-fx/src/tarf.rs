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
    let eps_s = spec.smoothing_eps.max(1e-6) * spec.strike.val().abs().max(1.0);
    let eps_t = spec.smoothing_eps.max(1e-6) * spec.target_profit.val().abs().max(1e-4);
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
            let gain = (x + root - T::from_f64(eps_s)) * 0.5;
            let loss = (-x + root - T::from_f64(eps_s)) * 0.5;
            let room = spec.target_profit - acc;
            let room_root = (room * room + T::from_f64(eps_t * eps_t)).sqrt();
            let alive = (T::one() + room / room_root) * 0.5;
            let positive_room = (room + room_root) * 0.5;
            let capped = match spec.cap_convention {
                FxTarfCapConvention::ExactCap => {
                    let d = positive_room - gain;
                    positive_room - (d + (d * d + T::from_f64(eps_t * eps_t)).sqrt()) * 0.5
                }
                FxTarfCapConvention::FullCouponAtKnockOut => alive * gain,
            };
            let cf = spec.notional
                * (-spec.r_dom * T::from_f64(time)).exp()
                * (capped - alive * spec.downside_leverage * loss);
            path_pv += cf;
            acc += capped;
            knocked += (T::one() - alive) * T::from_f64(time);
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
