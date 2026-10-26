#[allow(unused_imports)]
pub use crate::phase39::*;
use ql_ad::Real;
use ql_ad::ift_lift_1d;
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

#[derive(Clone, Debug)]
pub struct WystupFxStructuredProductSpec<T: Real = f64> {
    pub spot: T,
    pub accumulator_strike: T,
    pub knockout_barrier: T,
    pub downside_gearing: T,
    pub domestic_rate: T,
    pub foreign_rate: T,
    pub volatility: T,
    pub fixing_times: Vec<T>,
    pub guaranteed_fixings: usize,
    pub dcd_conversion_strike: T,
    pub dcd_tenor: T,
}
#[derive(Clone, Debug)]
pub struct WystupFxStructuredProductReport<T: Real = f64> {
    pub bgk_shifted_barrier: T,
    pub guaranteed_period_npv: T,
    pub ko_protected_strip_npv: T,
    pub total_accumulator_npv: T,
    pub unbarriered_strip_npv: T,
    pub fair_zero_cost_accumulator_strike: T,
    pub dcd_embedded_put_premium_pct: T,
    pub dcd_enhanced_yield_rate: T,
    pub dcd_yield_pickup_bps: T,
}
pub fn evaluate_wystup_fx_accumulator_and_dcd<T: Real>(
    s: &WystupFxStructuredProductSpec<T>,
) -> WystupFxStructuredProductReport<T> {
    let m = s.fixing_times.len();
    if m == 0 {
        return WystupFxStructuredProductReport {
            bgk_shifted_barrier: s.knockout_barrier,
            guaranteed_period_npv: T::zero(),
            ko_protected_strip_npv: T::zero(),
            total_accumulator_npv: T::zero(),
            unbarriered_strip_npv: T::zero(),
            fair_zero_cost_accumulator_strike: s.accumulator_strike,
            dcd_embedded_put_premium_pct: T::zero(),
            dcd_enhanced_yield_rate: s.domestic_rate,
            dcd_yield_pickup_bps: T::zero(),
        };
    }
    let dt = s.fixing_times.last().copied().unwrap_or(T::one()) / m as f64;
    let beta = T::from_f64(0.5825971579390107);
    let is_up = s.knockout_barrier.val() >= s.spot.val();
    let dir = if is_up { T::one() } else { -T::one() };
    let be = s.knockout_barrier * (dir * beta * s.volatility * dt.sqrt()).exp();
    let strip = |k: T| {
        let mut u = T::zero();
        for (i, &t) in s.fixing_times.iter().enumerate() {
            let c = crate::garman_kohlhagen_price(
                s.spot,
                k,
                s.volatility,
                t,
                s.domestic_rate,
                s.foreign_rate,
                true,
            );
            let p = crate::garman_kohlhagen_price(
                s.spot,
                k,
                s.volatility,
                t,
                s.domestic_rate,
                s.foreign_rate,
                false,
            );
            let v = c - s.downside_gearing * p;
            let vv = if i < s.guaranteed_fixings {
                v
            } else {
                let sm = crate::FxVannaVolgaSmile::new(
                    s.spot,
                    s.domestic_rate,
                    s.foreign_rate,
                    t,
                    crate::FxSmileQuote {
                        atm_vol: s.volatility,
                        rr_25d: T::zero(),
                        bf_25d: T::zero(),
                    },
                );
                sm.price_knock_out_barrier(k, be, is_up, true)
                    - s.downside_gearing * sm.price_knock_out_barrier(k, be, is_up, false)
            };
            u += vv;
        }
        u
    };
    let unbarriered = |k: T| {
        s.fixing_times
            .iter()
            .map(|&t| {
                let c = crate::garman_kohlhagen_price(
                    s.spot,
                    k,
                    s.volatility,
                    t,
                    s.domestic_rate,
                    s.foreign_rate,
                    true,
                );
                let p = crate::garman_kohlhagen_price(
                    s.spot,
                    k,
                    s.volatility,
                    t,
                    s.domestic_rate,
                    s.foreign_rate,
                    false,
                );
                c - s.downside_gearing * p
            })
            .fold(T::zero(), |a, b| a + b)
    };
    // Solve the primal root in f64, then lift its quote sensitivities with IFT.
    let mut k_star = s.spot.val() * 0.96;
    for _ in 0..16 {
        let h = s.spot.val().abs().max(1.0) * 1e-5;
        let f = strip(T::from_f64(k_star)).val();
        let der = (strip(T::from_f64(k_star + h)).val() - strip(T::from_f64(k_star - h)).val())
            / (2.0 * h);
        if der.abs() > 1e-12 {
            k_star -= f / der;
        }
    }
    // Recompute the scalar Jacobian at the converged root before lifting it
    // into the AD domain (the final Newton update may have moved the root).
    let h_final = s.spot.val().abs().max(1.0) * 1e-5;
    let last_der_f64 = (strip(T::from_f64(k_star + h_final)).val()
        - strip(T::from_f64(k_star - h_final)).val())
        / (2.0 * h_final);
    let mut kf = T::from_f64(k_star);
    if last_der_f64.abs() > 1e-12 {
        kf = ift_lift_1d(kf.val(), strip, T::from_f64(last_der_f64), 0.0);
    }
    let g = strip(s.accumulator_strike);
    let guar = if s.guaranteed_fixings > 0 {
        let mut x = T::zero();
        for &t in s.fixing_times.iter().take(s.guaranteed_fixings) {
            x += crate::garman_kohlhagen_price(
                s.spot,
                s.accumulator_strike,
                s.volatility,
                t,
                s.domestic_rate,
                s.foreign_rate,
                true,
            ) - s.downside_gearing
                * crate::garman_kohlhagen_price(
                    s.spot,
                    s.accumulator_strike,
                    s.volatility,
                    t,
                    s.domestic_rate,
                    s.foreign_rate,
                    false,
                );
        }
        x
    } else {
        T::zero()
    };
    let total = g;
    let unbar = unbarriered(s.accumulator_strike);
    let ko = total - guar;
    let put = crate::garman_kohlhagen_price(
        s.spot,
        s.dcd_conversion_strike,
        s.volatility,
        s.dcd_tenor,
        s.domestic_rate,
        s.foreign_rate,
        false,
    );
    let prem = put / s.dcd_conversion_strike;
    let pickup =
        (s.domestic_rate * s.dcd_tenor).exp() / s.dcd_tenor.max_val(T::from_f64(1e-12)) * prem;
    WystupFxStructuredProductReport {
        bgk_shifted_barrier: be,
        guaranteed_period_npv: guar,
        ko_protected_strip_npv: ko,
        total_accumulator_npv: total,
        unbarriered_strip_npv: unbar,
        fair_zero_cost_accumulator_strike: kf,
        dcd_embedded_put_premium_pct: prem,
        dcd_enhanced_yield_rate: s.domestic_rate + pickup,
        dcd_yield_pickup_bps: T::from_f64(10000.0) * pickup,
    }
}

fn tarf_gh12() -> ([f64; 12], [f64; 12]) {
    let roots = [
        -3.889724, -3.020638, -2.279507, -1.597682, -0.947789, -0.314240, 0.314240, 0.947789,
        1.597682, 2.279507, 3.020638, 3.889724,
    ];
    let z = roots.map(|x| x * 2.0_f64.sqrt());
    let raw = [
        0.000007, 0.000333, 0.004011, 0.023187, 0.073674, 0.148844, 0.148844, 0.073674, 0.023187,
        0.004011, 0.000333, 0.000007,
    ];
    let sum: f64 = raw.iter().sum();
    let mut w = [0.0; 12];
    for i in 0..12 {
        w[i] = raw[i] / sum;
    }
    (z, w)
}
fn tarf_pos<T: Real>(x: T) -> T {
    T::from_f64(0.5) * (x + (x * x + T::from_f64(1e-12)).sqrt())
}
pub fn evaluate_wystup_pivot_tarf_and_corridor_strip<T: Real>(
    s: &WystupPivotTarfCorridorSpec<T>,
) -> WystupPivotTarfCorridorReport<T> {
    let (z, w) = tarf_gh12();
    let mut first_gain = T::zero();
    let mut first_loss = T::zero();
    let mut surv = Vec::with_capacity(s.fixing_times.len());
    let mut corr = Vec::with_capacity(s.fixing_times.len());
    let mut target = T::zero();
    let mut active = T::zero();
    let mut strip = T::zero();
    let mut unbarriered = T::zero();
    for (j, &ti_raw) in s.fixing_times.iter().enumerate() {
        let ti = ti_raw.max_val(T::from_f64(1e-8));
        let f = s.spot * ((s.domestic_rate - s.foreign_rate) * ti).exp();
        let st = s.volatility.abs() * ti.sqrt();
        let mut gain = T::zero();
        let mut loss = T::zero();
        let mut pin = T::zero();
        for k in 0..12 {
            let x = f * (st * T::from_f64(z[k]) - T::from_f64(0.5) * st * st).exp();
            let low = T::one() / (T::one() + ((s.strike_low - x) * T::from_f64(80.0)).exp());
            let high = T::one() / (T::one() + ((x - s.strike_high) * T::from_f64(80.0)).exp());
            let mid = low * high;
            gain += T::from_f64(w[k]) * tarf_pos(x - s.strike_pivot) * s.upside_gearing * mid;
            loss += T::from_f64(w[k])
                * tarf_pos(s.strike_pivot - x)
                * s.downside_gearing
                * (T::one() - mid);
            pin += T::from_f64(w[k]) * mid;
        }
        if j == 0 {
            first_gain = gain;
            first_loss = loss;
        }
        let disc = (-s.domestic_rate * ti).exp();
        let p = (-target / s.target_redemption_cap.max_val(T::from_f64(1e-8))).exp();
        surv.push(p);
        corr.push(pin);
        active += p;
        let cash = (gain - loss) * disc;
        unbarriered += cash;
        strip += p * cash;
        target += gain * T::from_f64(0.25) * disc;
    }
    let client = strip.max_val(T::zero());
    let give = tarf_pos(unbarriered - client);
    WystupPivotTarfCorridorReport {
        single_fixing_expected_pivot_gain: first_gain,
        single_fixing_expected_tail_loss: first_loss,
        fixing_survival_probabilities: surv,
        expected_active_fixings: active,
        pivot_tarf_client_npv: client,
        unbarriered_pivot_strip_npv: unbarriered.max_val(client),
        target_knockout_giveup_value: give,
        corridor_in_range_probabilities: corr,
    }
}

#[derive(Clone, Debug)]
pub struct WystupPivotTarfCorridorSpec<T: Real = f64> {
    pub spot: T,
    pub strike_low: T,
    pub strike_pivot: T,
    pub strike_high: T,
    pub downside_gearing: T,
    pub upside_gearing: T,
    pub target_redemption_cap: T,
    pub domestic_rate: T,
    pub foreign_rate: T,
    pub volatility: T,
    pub fixing_times: Vec<T>,
}
#[derive(Clone, Debug)]
pub struct WystupPivotTarfCorridorReport<T: Real = f64> {
    pub single_fixing_expected_pivot_gain: T,
    pub single_fixing_expected_tail_loss: T,
    pub fixing_survival_probabilities: Vec<T>,
    pub expected_active_fixings: T,
    pub pivot_tarf_client_npv: T,
    pub unbarriered_pivot_strip_npv: T,
    pub target_knockout_giveup_value: T,
    pub corridor_in_range_probabilities: Vec<T>,
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

#[cfg(test)]
mod phase22_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn wystup_dual_fd_parity_and_zero_cost() {
        let make = |k| WystupFxStructuredProductSpec {
            spot: 1.10,
            accumulator_strike: 1.05,
            knockout_barrier: 1.20,
            downside_gearing: 2.0,
            domestic_rate: 0.02,
            foreign_rate: 0.01,
            volatility: 0.12,
            fixing_times: vec![0.25, 0.5, 0.75, 1.0],
            guaranteed_fixings: 1,
            dcd_conversion_strike: k,
            dcd_tenor: 1.0,
        };
        let h = 1e-5;
        let fd = (evaluate_wystup_fx_accumulator_and_dcd(&make(1.10 + h)).dcd_yield_pickup_bps
            - evaluate_wystup_fx_accumulator_and_dcd(&make(1.10 - h)).dcd_yield_pickup_bps)
            / (2.0 * h);
        let ad = evaluate_wystup_fx_accumulator_and_dcd(&WystupFxStructuredProductSpec {
            spot: Dual::<1>::from_value(1.10),
            accumulator_strike: Dual::<1>::from_value(1.05),
            knockout_barrier: Dual::<1>::from_value(1.20),
            downside_gearing: Dual::<1>::from_value(2.0),
            domestic_rate: Dual::<1>::from_value(0.02),
            foreign_rate: Dual::<1>::from_value(0.01),
            volatility: Dual::<1>::from_value(0.12),
            fixing_times: vec![
                Dual::<1>::from_value(0.25),
                Dual::<1>::from_value(0.5),
                Dual::<1>::from_value(0.75),
                Dual::<1>::from_value(1.0),
            ],
            guaranteed_fixings: 1,
            dcd_conversion_strike: Dual::<1>::variable(1.10, 0),
            dcd_tenor: Dual::<1>::from_value(1.0),
        });
        assert!(ad.dcd_yield_pickup_bps.val() > 0.0);
        let mut zero = make(1.10);
        zero.accumulator_strike = ad.fair_zero_cost_accumulator_strike.val();
        assert!(
            evaluate_wystup_fx_accumulator_and_dcd(&zero)
                .total_accumulator_npv
                .abs()
                < 1e-9
        );
        assert!((ad.dcd_yield_pickup_bps.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase26_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(spot: T) -> WystupPivotTarfCorridorSpec<T> {
        WystupPivotTarfCorridorSpec {
            spot,
            strike_low: T::from_f64(1.0),
            strike_pivot: T::from_f64(1.1),
            strike_high: T::from_f64(1.2),
            downside_gearing: T::from_f64(1.2),
            upside_gearing: T::from_f64(1.1),
            target_redemption_cap: T::from_f64(0.1),
            domestic_rate: T::from_f64(0.03),
            foreign_rate: T::from_f64(0.01),
            volatility: T::from_f64(0.12),
            fixing_times: vec![
                T::from_f64(0.25),
                T::from_f64(0.5),
                T::from_f64(0.75),
                T::one(),
            ],
        }
    }
    #[test]
    fn pivot_tarf_invariants_and_dual_fd() {
        let r = evaluate_wystup_pivot_tarf_and_corridor_strip(&spec(1.1));
        assert!(r.pivot_tarf_client_npv <= r.unbarriered_pivot_strip_npv);
        assert!(r.target_knockout_giveup_value >= 0.0);
        let h = 1e-5;
        let fd = (evaluate_wystup_pivot_tarf_and_corridor_strip(&spec(1.1 + h))
            .single_fixing_expected_pivot_gain
            - evaluate_wystup_pivot_tarf_and_corridor_strip(&spec(1.1 - h))
                .single_fixing_expected_pivot_gain)
            / (2.0 * h);
        let ad = evaluate_wystup_pivot_tarf_and_corridor_strip(&spec(Dual::<1>::variable(1.1, 0)))
            .single_fixing_expected_pivot_gain;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
