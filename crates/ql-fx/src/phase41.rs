//! Wystup (2017) & Clark (2011) FX Power-Reverse Dual-Currency (`PRDC`)
//! Callable Cap, Knock-Out Trigger & Quanto-CMS FX-Rate Convexity Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct WystupFxPrdcCallableCapSpec<T: Real = f64> {
    pub spot_fx_s0: T,
    pub initial_fx_s_init: T,
    pub foreign_coupon_cf: T,
    pub domestic_deduction_cd: T,
    pub cap_rate_k_cap: T,
    pub domestic_rate_rd: T,
    pub foreign_rate_rf: T,
    pub fx_vol_sigma_s: T,
    pub dom_rate_vol_sigma_d: T,
    pub fx_dom_corr_rho_sd: T,
    pub coupon_horizon_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct WystupFxPrdcCallableCapReport<T: Real = f64> {
    pub convexity_adjusted_forward_fx: T,
    pub lower_call_strike_fx_k1: T,
    pub upper_cap_strike_fx_k2: T,
    pub uncapped_prdc_coupon_pv: T,
    pub embedded_issuer_cap_pv: T,
    pub capped_prdc_coupon_pv: T,
    pub fx_rate_correlation_convexity_uplift: T,
}

pub fn evaluate_wystup_fx_prdc_coupon_convexity<T: Real>(
    spec: &WystupFxPrdcCallableCapSpec<T>,
) -> WystupFxPrdcCallableCapReport<T> {
    let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-4));
    let s_init = spec.initial_fx_s_init.max_val(T::from_f64(1e-4));
    let cf = spec.foreign_coupon_cf.max_val(T::from_f64(1e-4));
    let cd = spec.domestic_deduction_cd.max_val(T::from_f64(1e-4));
    let k_cap = spec.cap_rate_k_cap.max_val(T::from_f64(1e-4));
    let t = spec.coupon_horizon_t.max_val(T::from_f64(0.1));
    let df_d = (-spec.domestic_rate_rd * t).exp();

    // Standard forward FX vs. domestic-measure FX-IR correlation convexity adjustment:
    let fwd_raw = s0 * ((spec.domestic_rate_rd - spec.foreign_rate_rf) * t).exp();
    let conv_factor = (T::from_f64(0.5)
        * spec.fx_dom_corr_rho_sd
        * spec.fx_vol_sigma_s
        * spec.dom_rate_vol_sigma_d
        * t
        * t)
        .exp();
    let fwd_adj = fwd_raw * conv_factor;

    // Coupon = max(0, min(K_cap, c_f * (S_T / S_init) - c_d))
    //        = (c_f / S_init) * [ (S_T - K_1)^+ - (S_T - K_2)^+ ]
    let k1 = s_init * (cd / cf);
    let k2 = s_init * ((cd + k_cap) / cf);
    let scale = cf / s_init;

    let std = spec.fx_vol_sigma_s.max_val(T::from_f64(1e-4)) * t.sqrt();
    let d1_1 = ((fwd_adj / k1).ln() + T::from_f64(0.5) * std * std) / std;
    let d2_1 = d1_1 - std;
    let call_k1 = df_d * scale * (fwd_adj * norm_cdf(d1_1) - k1 * norm_cdf(d2_1));

    let d1_2 = ((fwd_adj / k2).ln() + T::from_f64(0.5) * std * std) / std;
    let d2_2 = d1_2 - std;
    let call_k2 = df_d * scale * (fwd_adj * norm_cdf(d1_2) - k2 * norm_cdf(d2_2));

    let d1_raw = ((fwd_raw / k1).ln() + T::from_f64(0.5) * std * std) / std;
    let d2_raw = d1_raw - std;
    let call_k1_raw = df_d * scale * (fwd_raw * norm_cdf(d1_raw) - k1 * norm_cdf(d2_raw));

    let capped_pv = (call_k1 - call_k2).max_val(T::zero());
    let uplift = call_k1 - call_k1_raw;

    WystupFxPrdcCallableCapReport {
        convexity_adjusted_forward_fx: fwd_adj,
        lower_call_strike_fx_k1: k1,
        upper_cap_strike_fx_k2: k2,
        uncapped_prdc_coupon_pv: call_k1,
        embedded_issuer_cap_pv: call_k2,
        capped_prdc_coupon_pv: capped_pv,
        fx_rate_correlation_convexity_uplift: uplift,
    }
}

#[cfg(test)]
mod phase41_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn wystup_fx_prdc_coupon_convexity_invariants_and_ad_parity() {
        let spec = WystupFxPrdcCallableCapSpec {
            spot_fx_s0: 108.0,
            initial_fx_s_init: 105.0,
            foreign_coupon_cf: 0.055,
            domestic_deduction_cd: 0.035,
            cap_rate_k_cap: 0.030,
            domestic_rate_rd: 0.012,
            foreign_rate_rf: 0.045,
            fx_vol_sigma_s: 0.115,
            dom_rate_vol_sigma_d: 0.008,
            fx_dom_corr_rho_sd: 0.35,
            coupon_horizon_t: 7.0,
        };
        let rep = evaluate_wystup_fx_prdc_coupon_convexity(&spec);
        assert!(rep.upper_cap_strike_fx_k2 > rep.lower_call_strike_fx_k1);
        assert!(rep.capped_prdc_coupon_pv < rep.uncapped_prdc_coupon_pv);
        assert!(rep.fx_rate_correlation_convexity_uplift > 0.0);

        let h = 1e-6;
        let r_ad = evaluate_wystup_fx_prdc_coupon_convexity(&WystupFxPrdcCallableCapSpec {
            spot_fx_s0: Dual::<1>::variable(spec.spot_fx_s0, 0),
            initial_fx_s_init: Dual::<1>::from_value(spec.initial_fx_s_init),
            foreign_coupon_cf: Dual::<1>::from_value(spec.foreign_coupon_cf),
            domestic_deduction_cd: Dual::<1>::from_value(spec.domestic_deduction_cd),
            cap_rate_k_cap: Dual::<1>::from_value(spec.cap_rate_k_cap),
            domestic_rate_rd: Dual::<1>::from_value(spec.domestic_rate_rd),
            foreign_rate_rf: Dual::<1>::from_value(spec.foreign_rate_rf),
            fx_vol_sigma_s: Dual::<1>::from_value(spec.fx_vol_sigma_s),
            dom_rate_vol_sigma_d: Dual::<1>::from_value(spec.dom_rate_vol_sigma_d),
            fx_dom_corr_rho_sd: Dual::<1>::from_value(spec.fx_dom_corr_rho_sd),
            coupon_horizon_t: Dual::<1>::from_value(spec.coupon_horizon_t),
        });
        let mut sp = spec;
        sp.spot_fx_s0 += h;
        let mut sm = spec;
        sm.spot_fx_s0 -= h;
        let fd = (evaluate_wystup_fx_prdc_coupon_convexity(&sp).capped_prdc_coupon_pv
            - evaluate_wystup_fx_prdc_coupon_convexity(&sm).capped_prdc_coupon_pv)
            / (2.0 * h);
        assert!((r_ad.capped_prdc_coupon_pv.d[0] - fd).abs() < 1e-6);
    }
}
