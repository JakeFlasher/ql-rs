#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;
use ql_math::{norm_cdf, norm_pdf};

#[derive(Clone, Copy, Debug)]
pub struct GatheralBergomiStickyLocalVolDeltaSpec<T: Real = f64> {
    pub spot_price: T,
    pub strike_price: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub forward_start_t1: T,
    pub maturity_t2: T,
    pub atm_vol: T,
    pub log_moneyness_skew_s0: T,
    pub sv_mean_reversion_kappa: T,
}

#[derive(Clone, Copy, Debug)]
pub struct GatheralBergomiStickyLocalVolDeltaReport<T: Real = f64> {
    pub bs_call_price: T,
    pub bs_vega: T,
    pub sticky_strike_delta: T,
    pub sticky_delta_moneyness_delta: T,
    pub dupire_local_vol_delta: T,
    pub shadow_vega_delta_adjustment: T,
    pub sticky_local_vol_symmetry_error: T,
    pub local_vol_forward_start_skew: T,
    pub stochastic_vol_forward_start_skew: T,
}

pub fn evaluate_gatheral_bergomi_sticky_vs_local_vol_delta<T: Real>(
    spec: &GatheralBergomiStickyLocalVolDeltaSpec<T>,
) -> GatheralBergomiStickyLocalVolDeltaReport<T> {
    let s0 = spec.spot_price.max_val(T::from_f64(1e-6));
    let k = spec.strike_price.max_val(T::from_f64(1e-6));
    let t2 = spec.maturity_t2.max_val(T::from_f64(1e-6));
    let sqrt_t = t2.sqrt();
    let df_q = (-spec.dividend_yield * t2).exp();
    let df_r = (-spec.risk_free_rate * t2).exp();
    let f = s0 * ((spec.risk_free_rate - spec.dividend_yield) * t2).exp();
    let k_log = (k / f).ln();
    let sig = (spec.atm_vol + spec.log_moneyness_skew_s0 * k_log).max_val(T::from_f64(0.02));
    let st = sig * sqrt_t;
    let d1 = (-k_log + T::from_f64(0.5) * st * st) / st;
    let d2 = d1 - st;

    let call = s0 * df_q * norm_cdf(d1) - k * df_r * norm_cdf(d2);
    let delta_ss = df_q * norm_cdf(d1);
    let vega = s0 * df_q * sqrt_t * norm_pdf(d1);
    let shadow_adj = (vega / s0) * spec.log_moneyness_skew_s0;

    let delta_sd = delta_ss - shadow_adj;
    let delta_lv = delta_ss + shadow_adj;
    let d_sym = (delta_sd + delta_lv) - T::from_f64(2.0) * delta_ss;
    let sym_err = (d_sym * d_sym + T::from_f64(1e-30)).sqrt();

    let dt_fwd = (t2 - spec.forward_start_t1).max_val(T::from_f64(1e-4));
    let lv_fwd_skew = spec.log_moneyness_skew_s0 * (dt_fwd / t2);
    let k_dt = (spec.sv_mean_reversion_kappa * dt_fwd).max_val(T::from_f64(1e-6));
    let sv_fwd_skew = spec.log_moneyness_skew_s0 * (T::one() - (-k_dt).exp()) / k_dt;

    GatheralBergomiStickyLocalVolDeltaReport {
        bs_call_price: call,
        bs_vega: vega,
        sticky_strike_delta: delta_ss,
        sticky_delta_moneyness_delta: delta_sd,
        dupire_local_vol_delta: delta_lv,
        shadow_vega_delta_adjustment: shadow_adj,
        sticky_local_vol_symmetry_error: sym_err,
        local_vol_forward_start_skew: lv_fwd_skew,
        stochastic_vol_forward_start_skew: sv_fwd_skew,
    }
}

#[cfg(test)]
mod phase37_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn gatheral_bergomi_sticky_vs_local_vol_delta_invariants_and_ad_parity() {
        let spec = GatheralBergomiStickyLocalVolDeltaSpec {
            spot_price: 100.0,
            strike_price: 100.0,
            risk_free_rate: 0.03,
            dividend_yield: 0.01,
            forward_start_t1: 0.5,
            maturity_t2: 1.5,
            atm_vol: 0.20,
            log_moneyness_skew_s0: -0.15,
            sv_mean_reversion_kappa: 1.2,
        };
        let rep = evaluate_gatheral_bergomi_sticky_vs_local_vol_delta(&spec);
        assert!(rep.sticky_local_vol_symmetry_error < 1e-12);
        assert!(rep.sticky_delta_moneyness_delta > rep.sticky_strike_delta);
        assert!(rep.sticky_strike_delta > rep.dupire_local_vol_delta);

        let h = 1e-6;
        let r_ad = evaluate_gatheral_bergomi_sticky_vs_local_vol_delta(
            &GatheralBergomiStickyLocalVolDeltaSpec {
                spot_price: Dual::<1>::variable(spec.spot_price, 0),
                strike_price: Dual::<1>::from_value(spec.strike_price),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                dividend_yield: Dual::<1>::from_value(spec.dividend_yield),
                forward_start_t1: Dual::<1>::from_value(spec.forward_start_t1),
                maturity_t2: Dual::<1>::from_value(spec.maturity_t2),
                atm_vol: Dual::<1>::from_value(spec.atm_vol),
                log_moneyness_skew_s0: Dual::<1>::from_value(spec.log_moneyness_skew_s0),
                sv_mean_reversion_kappa: Dual::<1>::from_value(spec.sv_mean_reversion_kappa),
            },
        );
        let mut sp = spec;
        sp.spot_price += h;
        let mut sm = spec;
        sm.spot_price -= h;
        let fd = (evaluate_gatheral_bergomi_sticky_vs_local_vol_delta(&sp).dupire_local_vol_delta
            - evaluate_gatheral_bergomi_sticky_vs_local_vol_delta(&sm).dupire_local_vol_delta)
            / (2.0 * h);
        assert!((r_ad.dupire_local_vol_delta.d[0] - fd).abs() < 1e-5);
    }
}
