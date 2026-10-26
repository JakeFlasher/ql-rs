#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Clone, Copy, Debug)]
pub struct BouzoubaaTimerOptionVarianceClockSpec<T: Real = f64> {
    pub spot_price: T,
    pub strike_price: T,
    pub risk_free_rate: T,
    pub target_budget_vol: T,
    pub nominal_horizon_years: T,
    pub market_implied_vol: T,
    pub expected_realized_vol: T,
    pub vol_of_vol_xi: T,
}

#[derive(Clone, Copy, Debug)]
pub struct BouzoubaaTimerOptionVarianceClockReport<T: Real = f64> {
    pub variance_budget_bv: T,
    pub expected_stopping_time_years: T,
    pub stopping_time_variance: T,
    pub zero_rate_dds_timer_call: T,
    pub stopping_time_rate_convexity_adj: T,
    pub timer_call_fair_value: T,
    pub vanilla_market_implied_call: T,
    pub timer_savings_over_vanilla: T,
}

pub fn evaluate_bouzoubaa_timer_option_variance_clock<T: Real>(
    spec: &BouzoubaaTimerOptionVarianceClockSpec<T>,
) -> BouzoubaaTimerOptionVarianceClockReport<T> {
    let s0 = spec.spot_price.max_val(T::from_f64(1e-6));
    let k = spec.strike_price.max_val(T::from_f64(1e-6));
    let t_nom = spec.nominal_horizon_years.max_val(T::from_f64(1e-6));
    let bv = (spec.target_budget_vol * spec.target_budget_vol * t_nom).max_val(T::from_f64(1e-8));
    let real_v2 =
        (spec.expected_realized_vol * spec.expected_realized_vol).max_val(T::from_f64(1e-6));
    let exp_tau = bv / real_v2;
    let var_tau = (spec.vol_of_vol_xi * spec.vol_of_vol_xi * bv) / (real_v2 * real_v2 * real_v2);

    let sqrt_bv = bv.sqrt();
    let d1_0 = ((s0 / k).ln() + T::from_f64(0.5) * bv) / sqrt_bv;
    let d2_0 = d1_0 - sqrt_bv;
    let dds_zero_call = s0 * norm_cdf(d1_0) - k * norm_cdf(d2_0);

    let r = spec.risk_free_rate;
    let d1_tau = ((s0 / k).ln() + r * exp_tau + T::from_f64(0.5) * bv) / sqrt_bv;
    let d2_tau = d1_tau - sqrt_bv;
    let base_timer = s0 * norm_cdf(d1_tau) - k * (-r * exp_tau).exp() * norm_cdf(d2_tau);
    let conv_adj = T::from_f64(0.5) * r * r * var_tau * k * (-r * exp_tau).exp() * norm_cdf(d2_tau);
    let timer_pv = base_timer + conv_adj;

    let st_van = (spec.market_implied_vol * t_nom.sqrt()).max_val(T::from_f64(1e-6));
    let d1_v = ((s0 / k).ln()
        + (r + T::from_f64(0.5) * spec.market_implied_vol * spec.market_implied_vol) * t_nom)
        / st_van;
    let d2_v = d1_v - st_van;
    let van_pv = s0 * norm_cdf(d1_v) - k * (-r * t_nom).exp() * norm_cdf(d2_v);

    BouzoubaaTimerOptionVarianceClockReport {
        variance_budget_bv: bv,
        expected_stopping_time_years: exp_tau,
        stopping_time_variance: var_tau,
        zero_rate_dds_timer_call: dds_zero_call,
        stopping_time_rate_convexity_adj: conv_adj,
        timer_call_fair_value: timer_pv,
        vanilla_market_implied_call: van_pv,
        timer_savings_over_vanilla: van_pv - timer_pv,
    }
}

#[cfg(test)]
mod phase37_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn bouzoubaa_timer_option_variance_clock_invariants_and_ad_parity() {
        let spec = BouzoubaaTimerOptionVarianceClockSpec {
            spot_price: 100.0,
            strike_price: 100.0,
            risk_free_rate: 0.02,
            target_budget_vol: 0.16,
            nominal_horizon_years: 1.0,
            market_implied_vol: 0.22,
            expected_realized_vol: 0.20,
            vol_of_vol_xi: 0.30,
        };
        let rep = evaluate_bouzoubaa_timer_option_variance_clock(&spec);
        assert!(rep.timer_savings_over_vanilla > 0.0);
        assert!(rep.expected_stopping_time_years < spec.nominal_horizon_years);

        let h = 1e-6;
        let r_ad = evaluate_bouzoubaa_timer_option_variance_clock(
            &BouzoubaaTimerOptionVarianceClockSpec {
                spot_price: Dual::<1>::variable(spec.spot_price, 0),
                strike_price: Dual::<1>::from_value(spec.strike_price),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                target_budget_vol: Dual::<1>::from_value(spec.target_budget_vol),
                nominal_horizon_years: Dual::<1>::from_value(spec.nominal_horizon_years),
                market_implied_vol: Dual::<1>::from_value(spec.market_implied_vol),
                expected_realized_vol: Dual::<1>::from_value(spec.expected_realized_vol),
                vol_of_vol_xi: Dual::<1>::from_value(spec.vol_of_vol_xi),
            },
        );
        let mut sp = spec;
        sp.spot_price += h;
        let mut sm = spec;
        sm.spot_price -= h;
        let fd = (evaluate_bouzoubaa_timer_option_variance_clock(&sp).timer_call_fair_value
            - evaluate_bouzoubaa_timer_option_variance_clock(&sm).timer_call_fair_value)
            / (2.0 * h);
        assert!((r_ad.timer_call_fair_value.d[0] - fd).abs() < 1e-5);
    }
}
