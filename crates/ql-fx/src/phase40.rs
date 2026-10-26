//! Wystup (2017) & Clark (2011) FX Power-Reverse Dual-Currency (`PRDC`) Coupon &
//! 3-Factor FX-IR Correlation Convexity Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct WystupPrdcThreeFactorCouponSpec<T: Real = f64> {
    pub initial_spot_fx_s0: T,
    pub current_spot_fx_st: T,
    pub foreign_coupon_rate_cf: T,
    pub domestic_strike_coupon_cd: T,
    pub domestic_rate_rd: T,
    pub foreign_rate_rf: T,
    pub spot_fx_vol_sigma_s: T,
    pub dom_rate_vol_sigma_rd: T,
    pub for_rate_vol_sigma_rf: T,
    pub corr_fx_dom_rho_sd: T,
    pub corr_fx_for_rho_sf: T,
    pub maturity_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct WystupPrdcThreeFactorCouponReport<T: Real = f64> {
    pub prdc_fx_strike_k: T,
    pub prdc_fx_leverage_notional: T,
    pub forward_fx_rate: T,
    pub three_factor_hybrid_fx_vol: T,
    pub ir_fx_correlation_vol_uplift: T,
    pub prdc_coupon_present_value: T,
}

pub fn evaluate_wystup_prdc_three_factor_coupon<T: Real>(
    spec: &WystupPrdcThreeFactorCouponSpec<T>,
) -> WystupPrdcThreeFactorCouponReport<T> {
    let s0 = spec.initial_spot_fx_s0.max_val(T::from_f64(1e-4));
    let st = spec.current_spot_fx_st.max_val(T::from_f64(1e-4));
    let cf = spec.foreign_coupon_rate_cf.max_val(T::from_f64(1e-4));
    let cd = spec.domestic_strike_coupon_cd.max_val(T::from_f64(1e-4));
    let t = spec.maturity_years_t.max_val(T::from_f64(0.1));

    // PRDC coupon payoff max(0, c_f * (S_T / S_0) - c_d) = (c_f / S_0) * max(0, S_T - K_PRDC)
    let k_prdc = s0 * (cd / cf);
    let lev_notional = cf / s0;
    let fwd_fx = st * ((spec.domestic_rate_rd - spec.foreign_rate_rf) * t).exp();

    // 3-Factor FX-IR hybrid forward FX variance sigma_PRDC^2(T):
    let s_fx = spec.spot_fx_vol_sigma_s;
    let s_rd = spec.dom_rate_vol_sigma_rd;
    let s_rf = spec.for_rate_vol_sigma_rf;
    let ir_var_term = (t * t / T::from_f64(3.0)) * (s_rd * s_rd + s_rf * s_rf);
    let cross_fx_ir =
        t * (spec.corr_fx_dom_rho_sd * s_fx * s_rd - spec.corr_fx_for_rho_sf * s_fx * s_rf);
    let hybrid_var = (s_fx * s_fx + ir_var_term + cross_fx_ir).max_val(T::from_f64(1e-8));
    let hybrid_vol = hybrid_var.sqrt();
    let vol_uplift = hybrid_vol - s_fx;

    let std = hybrid_vol * t.sqrt();
    let d1 = ((fwd_fx / k_prdc).ln() + T::from_f64(0.5) * std * std) / std;
    let d2 = d1 - std;
    let df_d = (-spec.domestic_rate_rd * t).exp();
    let pv = lev_notional * df_d * (fwd_fx * norm_cdf(d1) - k_prdc * norm_cdf(d2));

    WystupPrdcThreeFactorCouponReport {
        prdc_fx_strike_k: k_prdc,
        prdc_fx_leverage_notional: lev_notional,
        forward_fx_rate: fwd_fx,
        three_factor_hybrid_fx_vol: hybrid_vol,
        ir_fx_correlation_vol_uplift: vol_uplift,
        prdc_coupon_present_value: pv,
    }
}

#[cfg(test)]
mod phase40_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn wystup_prdc_three_factor_coupon_invariants_and_ad_parity() {
        let spec = WystupPrdcThreeFactorCouponSpec {
            initial_spot_fx_s0: 105.0,
            current_spot_fx_st: 108.0,
            foreign_coupon_rate_cf: 0.055,
            domestic_strike_coupon_cd: 0.035,
            domestic_rate_rd: 0.010,
            foreign_rate_rf: 0.045,
            spot_fx_vol_sigma_s: 0.11,
            dom_rate_vol_sigma_rd: 0.007,
            for_rate_vol_sigma_rf: 0.012,
            corr_fx_dom_rho_sd: 0.25,
            corr_fx_for_rho_sf: -0.40,
            maturity_years_t: 5.0,
        };
        let rep = evaluate_wystup_prdc_three_factor_coupon(&spec);
        assert!(rep.ir_fx_correlation_vol_uplift > 0.0);
        assert!(rep.prdc_coupon_present_value > 0.0);

        let h = 1e-6;
        let r_ad = evaluate_wystup_prdc_three_factor_coupon(&WystupPrdcThreeFactorCouponSpec {
            initial_spot_fx_s0: Dual::<1>::from_value(spec.initial_spot_fx_s0),
            current_spot_fx_st: Dual::<1>::variable(spec.current_spot_fx_st, 0),
            foreign_coupon_rate_cf: Dual::<1>::from_value(spec.foreign_coupon_rate_cf),
            domestic_strike_coupon_cd: Dual::<1>::from_value(spec.domestic_strike_coupon_cd),
            domestic_rate_rd: Dual::<1>::from_value(spec.domestic_rate_rd),
            foreign_rate_rf: Dual::<1>::from_value(spec.foreign_rate_rf),
            spot_fx_vol_sigma_s: Dual::<1>::from_value(spec.spot_fx_vol_sigma_s),
            dom_rate_vol_sigma_rd: Dual::<1>::from_value(spec.dom_rate_vol_sigma_rd),
            for_rate_vol_sigma_rf: Dual::<1>::from_value(spec.for_rate_vol_sigma_rf),
            corr_fx_dom_rho_sd: Dual::<1>::from_value(spec.corr_fx_dom_rho_sd),
            corr_fx_for_rho_sf: Dual::<1>::from_value(spec.corr_fx_for_rho_sf),
            maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
        });
        let mut sp = spec;
        sp.current_spot_fx_st += h;
        let mut sm = spec;
        sm.current_spot_fx_st -= h;
        let fd = (evaluate_wystup_prdc_three_factor_coupon(&sp).prdc_coupon_present_value
            - evaluate_wystup_prdc_three_factor_coupon(&sm).prdc_coupon_present_value)
            / (2.0 * h);
        assert!((r_ad.prdc_coupon_present_value.d[0] - fd).abs() < 1e-6);
    }
}
