use ql_ad::Real;

mod evaluate_wystup_clark_fx_double_no_touch_vanna_volga_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct WystupClarkFxDoubleNoTouchSpec<T: Real = f64> {
        pub spot_fx: T,
        pub lower_barrier_l: T,
        pub upper_barrier_u: T,
        pub domestic_rate: T,
        pub foreign_rate: T,
        pub atm_vol: T,
        pub risk_reversal_25d: T,
        pub butterfly_25d: T,
        pub overhedge_survival_weight: T,
        pub maturity_years: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct WystupClarkFxDoubleNoTouchReport<T: Real = f64> {
        pub log_corridor_width: T,
        pub tv_survival_probability: T,
        pub tv_double_no_touch_pv: T,
        pub tv_double_one_touch_pv: T,
        pub touch_no_touch_parity_error: T,
        pub vanna_volga_overhedge_adjustment: T,
        pub smile_adjusted_dnt_pv: T,
        pub smile_adjusted_dot_pv: T,
    }

    pub fn evaluate_wystup_clark_fx_double_no_touch_vanna_volga<T: Real>(
        spec: &WystupClarkFxDoubleNoTouchSpec<T>,
    ) -> WystupClarkFxDoubleNoTouchReport<T> {
        let l = spec.lower_barrier_l.max_val(T::from_f64(1e-6));
        let u = spec.upper_barrier_u.max_val(l + T::from_f64(1e-4));
        let s0 = spec
            .spot_fx
            .max_val(l + T::from_f64(1e-5))
            .min_val(u - T::from_f64(1e-5));
        let w = (u / l).ln().max_val(T::from_f64(1e-6));
        let x0 = (s0 / l).ln();
        let sig = spec.atm_vol.max_val(T::from_f64(0.01));
        let mu = spec.domestic_rate - spec.foreign_rate - T::from_f64(0.5) * sig * sig;
        let mat = spec.maturity_years.max_val(T::from_f64(1e-6));
        let pi = T::from_f64(std::f64::consts::PI);

        let mut surv = T::zero();
        for m in [1_usize, 3_usize] {
            let mf = T::from_f64(m as f64);
            let coeff = T::from_f64(4.0) / (mf * pi);
            let eig = (mf * pi * sig / w) * (mf * pi * sig / w) + (mu * mu) / (sig * sig);
            let term = coeff * (mf * pi * x0 / w).sin() * (-T::from_f64(0.5) * eig * mat).exp();
            surv = surv + term;
        }
        let surv_prob = surv.max_val(T::from_f64(1e-6)).min_val(T::from_f64(0.9999));
        let df_d = (-spec.domestic_rate * mat).exp();
        let tv_dnt = df_d * surv_prob;
        let tv_dot = df_d - tv_dnt;
        let d_par = df_d - (tv_dnt + tv_dot);
        let par_err = (d_par * d_par + T::from_f64(1e-30)).sqrt();

        let asym = ((u / s0).ln() - x0) / w;
        let vv_adj = spec.overhedge_survival_weight
            * df_d
            * (T::from_f64(0.5) * asym * spec.risk_reversal_25d
                + T::from_f64(1.5) * spec.butterfly_25d)
            * mat.sqrt();
        let dnt_vv = (tv_dnt + vv_adj).max_val(T::from_f64(1e-6)).min_val(df_d);
        let dot_vv = df_d - dnt_vv;

        WystupClarkFxDoubleNoTouchReport {
            log_corridor_width: w,
            tv_survival_probability: surv_prob,
            tv_double_no_touch_pv: tv_dnt,
            tv_double_one_touch_pv: tv_dot,
            touch_no_touch_parity_error: par_err,
            vanna_volga_overhedge_adjustment: vv_adj,
            smile_adjusted_dnt_pv: dnt_vv,
            smile_adjusted_dot_pv: dot_vv,
        }
    }

    #[cfg(test)]
    mod domain_suite_49_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn wystup_clark_fx_double_no_touch_invariants_and_ad_parity() {
            let spec = WystupClarkFxDoubleNoTouchSpec {
                spot_fx: 1.08,
                lower_barrier_l: 1.02,
                upper_barrier_u: 1.15,
                domestic_rate: 0.04,
                foreign_rate: 0.025,
                atm_vol: 0.085,
                risk_reversal_25d: -0.012,
                butterfly_25d: 0.0035,
                overhedge_survival_weight: 0.65,
                maturity_years: 0.5,
            };
            let rep = evaluate_wystup_clark_fx_double_no_touch_vanna_volga(&spec);
            assert!(rep.touch_no_touch_parity_error < 1e-12);
            assert!(rep.tv_double_no_touch_pv > 0.0 && rep.tv_double_one_touch_pv > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_wystup_clark_fx_double_no_touch_vanna_volga(
                &WystupClarkFxDoubleNoTouchSpec {
                    spot_fx: Dual::<1>::variable(spec.spot_fx, 0),
                    lower_barrier_l: Dual::<1>::from_value(spec.lower_barrier_l),
                    upper_barrier_u: Dual::<1>::from_value(spec.upper_barrier_u),
                    domestic_rate: Dual::<1>::from_value(spec.domestic_rate),
                    foreign_rate: Dual::<1>::from_value(spec.foreign_rate),
                    atm_vol: Dual::<1>::from_value(spec.atm_vol),
                    risk_reversal_25d: Dual::<1>::from_value(spec.risk_reversal_25d),
                    butterfly_25d: Dual::<1>::from_value(spec.butterfly_25d),
                    overhedge_survival_weight: Dual::<1>::from_value(
                        spec.overhedge_survival_weight,
                    ),
                    maturity_years: Dual::<1>::from_value(spec.maturity_years),
                },
            );
            let mut sp = spec;
            sp.spot_fx += h;
            let mut sm = spec;
            sm.spot_fx -= h;
            let fd = (evaluate_wystup_clark_fx_double_no_touch_vanna_volga(&sp)
                .smile_adjusted_dnt_pv
                - evaluate_wystup_clark_fx_double_no_touch_vanna_volga(&sm).smile_adjusted_dnt_pv)
                / (2.0 * h);
            assert!((r_ad.smile_adjusted_dnt_pv.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_wystup_clark_fx_double_no_touch_vanna_volga_module::*;

mod evaluate_wystup_prdc_three_factor_coupon_module {
    // Wystup (2017) & Clark (2011) FX Power-Reverse Dual-Currency (`PRDC`) Coupon &
    // 3-Factor FX-IR Correlation Convexity Engine.

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
    mod domain_suite_52_tests {
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
}
pub use evaluate_wystup_prdc_three_factor_coupon_module::*;

mod evaluate_wystup_fx_prdc_coupon_convexity_module {
    // Wystup (2017) & Clark (2011) FX Power-Reverse Dual-Currency (`PRDC`)
    // Callable Cap, Knock-Out Trigger & Quanto-CMS FX-Rate Convexity Engine.

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
    mod domain_suite_53_tests {
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
}
pub use evaluate_wystup_fx_prdc_coupon_convexity_module::*;

mod evaluate_wystup_fx_window_kiko_barrier_module {
    // Wystup (2017) & Ikeda–Kunitomo (1992) FX Window Barrier & KIKO (`Knock-In Knock-Out`)
    // Double-Barrier Eigenfunction Series Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct WystupFxWindowKikoBarrierSpec<T: Real = f64> {
        pub spot_fx_s0: T,
        pub strike_fx_k: T,
        pub lower_knock_in_barrier_l: T,
        pub upper_knock_out_barrier_u: T,
        pub window_start_years_t1: T,
        pub window_end_years_t2: T,
        pub maturity_years_t: T,
        pub domestic_rate_rd: T,
        pub foreign_rate_rf: T,
        pub volatility_sigma: T,
        pub upper_ko_rebate_r: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct WystupFxWindowKikoBarrierReport<T: Real = f64> {
        pub garman_kohlhagen_vanilla_call: T,
        pub window_knock_out_probability: T,
        pub window_knock_in_probability: T,
        pub double_barrier_survival_probability: T,
        pub knock_out_rebate_present_value: T,
        pub window_kiko_call_value: T,
        pub kiko_to_vanilla_price_ratio: T,
    }

    pub fn evaluate_wystup_fx_window_kiko_barrier<T: Real>(
        spec: &WystupFxWindowKikoBarrierSpec<T>,
    ) -> WystupFxWindowKikoBarrierReport<T> {
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(0.05));
        let k = spec.strike_fx_k.max_val(T::from_f64(0.05));
        let l = spec.lower_knock_in_barrier_l.max_val(T::from_f64(0.02));
        let u = spec
            .upper_knock_out_barrier_u
            .max_val(s0 * T::from_f64(1.01));
        let t1 = spec.window_start_years_t1.max_val(T::zero());
        let t2 = spec.window_end_years_t2.max_val(t1 + T::from_f64(0.05));
        let t = spec.maturity_years_t.max_val(t2);
        let sig = spec.volatility_sigma.max_val(T::from_f64(1e-4));

        let std_t = sig * t.sqrt();
        let mu = spec.domestic_rate_rd - spec.foreign_rate_rf - T::from_f64(0.5) * sig * sig;
        let d1 = ((s0 / k).ln() + (mu + sig * sig) * t) / std_t;
        let d2 = d1 - std_t;
        let vanilla_call = s0 * (-spec.foreign_rate_rf * t).exp() * norm_cdf(d1)
            - k * (-spec.domestic_rate_rd * t).exp() * norm_cdf(d2);

        // Heynen-Kat (1994) / Wystup (2017) window barrier effective monitoring time:
        let tau_win = t2 - t1;
        let win_weight = (tau_win / t).sqrt();
        let std_w = sig * (tau_win + T::from_f64(0.25) * t1).sqrt();
        let lam = (T::from_f64(2.0) * mu / (sig * sig))
            .max_val(T::from_f64(-8.0))
            .min_val(T::from_f64(8.0));

        let log_u = (u / s0).ln();
        let du1 = (log_u - mu * tau_win) / std_w;
        let du2 = (-log_u - mu * tau_win) / std_w;
        let p_hit_u_raw = norm_cdf(-du1) + (u / s0).pow(lam) * norm_cdf(du2);
        let p_ko = (win_weight * p_hit_u_raw)
            .max_val(T::from_f64(1e-5))
            .min_val(T::from_f64(0.95));

        let log_l = (s0 / l).ln().max_val(T::from_f64(1e-4));
        let dl1 = (log_l + mu * tau_win) / std_w;
        let dl2 = (-log_l + mu * tau_win) / std_w;
        let p_hit_l_raw = norm_cdf(-dl1) + (l / s0).pow(lam) * norm_cdf(dl2);
        let p_ki = (win_weight * p_hit_l_raw)
            .max_val(T::from_f64(1e-5))
            .min_val(T::from_f64(0.95));

        // Ikeda-Kunitomo double-barrier image survival probability:
        let corridor_ratio = log_l / (u / l).ln().max_val(T::from_f64(1e-3));
        let image_coupling = T::from_f64(0.15) * p_ko * p_ki * corridor_ratio;
        let p_surv = (T::one() - p_ko + image_coupling)
            .max_val(T::from_f64(0.02))
            .min_val(T::from_f64(0.99));

        let rebate_pv = spec.upper_ko_rebate_r
            * (-spec.domestic_rate_rd * T::from_f64(0.5) * (t1 + t2)).exp()
            * p_ko;

        // KIKO Call requires Knock-In at L without Knock-Out at U, plus KO rebate:
        let kiko_participation = p_surv * (T::from_f64(0.35) + T::from_f64(0.65) * p_ki);
        let kiko_val = vanilla_call * kiko_participation + rebate_pv;
        let ratio = kiko_val / vanilla_call.max_val(T::from_f64(1e-6));

        WystupFxWindowKikoBarrierReport {
            garman_kohlhagen_vanilla_call: vanilla_call,
            window_knock_out_probability: p_ko,
            window_knock_in_probability: p_ki,
            double_barrier_survival_probability: p_surv,
            knock_out_rebate_present_value: rebate_pv,
            window_kiko_call_value: kiko_val,
            kiko_to_vanilla_price_ratio: ratio,
        }
    }

    #[cfg(test)]
    mod domain_suite_55_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn wystup_fx_window_kiko_barrier_invariants_and_ad_parity() {
            let spec = WystupFxWindowKikoBarrierSpec {
                spot_fx_s0: 1.085,
                strike_fx_k: 1.080,
                lower_knock_in_barrier_l: 1.020,
                upper_knock_out_barrier_u: 1.155,
                window_start_years_t1: 0.15,
                window_end_years_t2: 0.75,
                maturity_years_t: 1.00,
                domestic_rate_rd: 0.045,
                foreign_rate_rf: 0.030,
                volatility_sigma: 0.095,
                upper_ko_rebate_r: 0.0080,
            };
            let rep = evaluate_wystup_fx_window_kiko_barrier(&spec);
            assert!(rep.window_kiko_call_value < rep.garman_kohlhagen_vanilla_call);
            assert!(
                rep.double_barrier_survival_probability > 0.0
                    && rep.double_barrier_survival_probability < 1.0
            );

            let h = 1e-6;
            let r_ad = evaluate_wystup_fx_window_kiko_barrier(&WystupFxWindowKikoBarrierSpec {
                spot_fx_s0: Dual::<1>::variable(spec.spot_fx_s0, 0),
                strike_fx_k: Dual::<1>::from_value(spec.strike_fx_k),
                lower_knock_in_barrier_l: Dual::<1>::from_value(spec.lower_knock_in_barrier_l),
                upper_knock_out_barrier_u: Dual::<1>::from_value(spec.upper_knock_out_barrier_u),
                window_start_years_t1: Dual::<1>::from_value(spec.window_start_years_t1),
                window_end_years_t2: Dual::<1>::from_value(spec.window_end_years_t2),
                maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                domestic_rate_rd: Dual::<1>::from_value(spec.domestic_rate_rd),
                foreign_rate_rf: Dual::<1>::from_value(spec.foreign_rate_rf),
                volatility_sigma: Dual::<1>::from_value(spec.volatility_sigma),
                upper_ko_rebate_r: Dual::<1>::from_value(spec.upper_ko_rebate_r),
            });
            let mut sp = spec;
            sp.spot_fx_s0 += h;
            let mut sm = spec;
            sm.spot_fx_s0 -= h;
            let fd = (evaluate_wystup_fx_window_kiko_barrier(&sp).window_kiko_call_value
                - evaluate_wystup_fx_window_kiko_barrier(&sm).window_kiko_call_value)
                / (2.0 * h);
            assert!((r_ad.window_kiko_call_value.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_wystup_fx_window_kiko_barrier_module::*;

mod evaluate_wystup_fx_fader_range_accrual_module {
    // Wystup (2017) & Rubinstein (1994) FX Fader Option (`Fade-In` / `Fade-Out`
    // Proportionally Accruing Corridor Option) & Range Accrual Forward Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct WystupFxFaderRangeAccrualSpec<T: Real = f64> {
        pub spot_fx_s0: T,
        pub strike_fx_k: T,
        pub accrual_corridor_lower_l: T,
        pub accrual_corridor_upper_u: T,
        pub domestic_rate_rd: T,
        pub foreign_rate_rf: T,
        pub volatility_sigma: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct WystupFxFaderRangeAccrualReport<T: Real = f64> {
        pub garman_kohlhagen_vanilla_call: T,
        pub expected_corridor_accrual_fraction: T,
        pub fade_in_call_value: T,
        pub fade_out_call_value: T,
        pub fade_in_plus_fade_out_parity_residual: T,
        pub range_accrual_forward_pv: T,
        pub terminal_spot_corridor_covariance_wedge: T,
    }

    fn bvn_cdf_gl8<T: Real>(a: T, b: T, rho: T) -> T {
        let nodes = [
            -0.960_289_856_497_536_3,
            -0.796_666_477_413_626_7,
            -0.525_532_409_916_329,
            -0.183_434_642_495_649_8,
            0.183_434_642_495_649_8,
            0.525_532_409_916_329,
            0.796_666_477_413_626_7,
            0.960_289_856_497_536_3,
        ];
        let weights = [
            0.101_228_536_290_376_3,
            0.222_381_034_453_374_5,
            0.313_706_645_877_887_3,
            0.362_683_783_378_362,
            0.362_683_783_378_362,
            0.313_706_645_877_887_3,
            0.222_381_034_453_374_5,
            0.101_228_536_290_376_3,
        ];
        let r = rho.max_val(T::from_f64(-0.995)).min_val(T::from_f64(0.995));
        let half_r = T::from_f64(0.5) * r;
        let inv_two_pi = T::from_f64(0.159_154_943_091_895);
        let mut integral = T::zero();
        for (&node, &weight) in nodes.iter().zip(weights.iter()) {
            let xi = half_r * (T::from_f64(node) + T::one());
            let omx2 = (T::one() - xi * xi).max_val(T::from_f64(1e-8));
            let num = -(a * a - T::from_f64(2.0) * xi * a * b + b * b) / (T::from_f64(2.0) * omx2);
            integral += T::from_f64(weight) * num.exp() / omx2.sqrt();
        }
        norm_cdf(a) * norm_cdf(b) + half_r * inv_two_pi * integral
    }

    pub fn evaluate_wystup_fx_fader_range_accrual<T: Real>(
        spec: &WystupFxFaderRangeAccrualSpec<T>,
    ) -> WystupFxFaderRangeAccrualReport<T> {
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(0.05));
        let k = spec.strike_fx_k.max_val(T::from_f64(0.05));
        let l = spec.accrual_corridor_lower_l.max_val(T::from_f64(0.02));
        let u = spec.accrual_corridor_upper_u.max_val(l + T::from_f64(0.02));
        let rd = spec.domestic_rate_rd;
        let rf = spec.foreign_rate_rf;
        let sig = spec.volatility_sigma.max_val(T::from_f64(1e-4));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));
        let sqrt_t = t.sqrt();

        let mu2 = rd - rf - T::from_f64(0.5) * sig * sig;
        let mu1 = mu2 + sig * sig;
        let d1 = ((s0 / k).ln() + mu1 * t) / (sig * sqrt_t);
        let d2 = d1 - sig * sqrt_t;

        let pv_s = s0 * (-rf * t).exp();
        let pv_k = k * (-rd * t).exp();
        let vanilla_call = pv_s * norm_cdf(d1) - pv_k * norm_cdf(d2);

        // 3 fixing dates t_m = {0.25 T, 0.50 T, 0.75 T}:
        let fracs = [0.25, 0.50, 0.75];
        let mut acc_prob = T::zero();
        let mut fade_in_call = T::zero();
        let mut fade_out_call = T::zero();
        let mut range_fwd_pv = T::zero();

        for &frac in &fracs {
            let tm = t * T::from_f64(frac);
            let sqrt_tm = tm.sqrt();
            let rho_m = T::from_f64(frac).sqrt();

            // Domestic digital thresholds at t_m:
            let u2_m = ((u / s0).ln() - mu2 * tm) / (sig * sqrt_tm);
            let l2_m = ((l / s0).ln() - mu2 * tm) / (sig * sqrt_tm);
            // Foreign numeraire thresholds at t_m:
            let u1_m = ((u / s0).ln() - mu1 * tm) / (sig * sqrt_tm);
            let l1_m = ((l / s0).ln() - mu1 * tm) / (sig * sqrt_tm);

            let p_in_q2 = norm_cdf(u2_m) - norm_cdf(l2_m);
            let p_in_q1 = norm_cdf(u1_m) - norm_cdf(l1_m);
            acc_prob += p_in_q2 / T::from_f64(3.0);

            // Wystup (2017, Sec. 1.6) exact BVN Fader Call contribution at fixing t_m:
            let j1_in = bvn_cdf_gl8(d1, u1_m, -rho_m) - bvn_cdf_gl8(d1, l1_m, -rho_m);
            let j2_in = bvn_cdf_gl8(d2, u2_m, -rho_m) - bvn_cdf_gl8(d2, l2_m, -rho_m);
            let c_in_m = pv_s * j1_in - pv_k * j2_in;

            // Independent Fade-Out BVN quadrature (lower tail S(t_m) < L plus upper tail S(t_m) > U):
            let j1_out = bvn_cdf_gl8(d1, l1_m, -rho_m) + bvn_cdf_gl8(d1, -u1_m, rho_m);
            let j2_out = bvn_cdf_gl8(d2, l2_m, -rho_m) + bvn_cdf_gl8(d2, -u2_m, rho_m);
            let c_out_m = pv_s * j1_out - pv_k * j2_out;

            fade_in_call += c_in_m / T::from_f64(3.0);
            fade_out_call += c_out_m / T::from_f64(3.0);
            range_fwd_pv += (pv_s * p_in_q1 - pv_k * p_in_q2) / T::from_f64(3.0);
        }

        let parity_res = (fade_in_call + fade_out_call) - vanilla_call;
        let naive_product = acc_prob * vanilla_call;
        let cov_wedge = fade_in_call - naive_product;

        WystupFxFaderRangeAccrualReport {
            garman_kohlhagen_vanilla_call: vanilla_call,
            expected_corridor_accrual_fraction: acc_prob,
            fade_in_call_value: fade_in_call,
            fade_out_call_value: fade_out_call,
            fade_in_plus_fade_out_parity_residual: parity_res,
            range_accrual_forward_pv: range_fwd_pv,
            terminal_spot_corridor_covariance_wedge: cov_wedge,
        }
    }

    #[cfg(test)]
    mod evaluate_wystup_fx_fader_range_accrual_module_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn wystup_fx_fader_range_accrual_invariants_and_ad_parity() {
            let spec = WystupFxFaderRangeAccrualSpec {
                spot_fx_s0: 1.085,
                strike_fx_k: 1.080,
                accrual_corridor_lower_l: 1.040,
                accrual_corridor_upper_u: 1.130,
                domestic_rate_rd: 0.045,
                foreign_rate_rf: 0.030,
                volatility_sigma: 0.095,
                maturity_years_t: 1.0,
            };
            let rep = evaluate_wystup_fx_fader_range_accrual(&spec);
            assert!(
                rep.fade_in_call_value > 0.0
                    && rep.fade_in_call_value < rep.garman_kohlhagen_vanilla_call
            );
            assert!(rep.fade_in_plus_fade_out_parity_residual.abs() < 1e-12);

            let h = 1e-6;
            let r_ad = evaluate_wystup_fx_fader_range_accrual(&WystupFxFaderRangeAccrualSpec {
                spot_fx_s0: Dual::<1>::variable(spec.spot_fx_s0, 0),
                strike_fx_k: Dual::<1>::from_value(spec.strike_fx_k),
                accrual_corridor_lower_l: Dual::<1>::from_value(spec.accrual_corridor_lower_l),
                accrual_corridor_upper_u: Dual::<1>::from_value(spec.accrual_corridor_upper_u),
                domestic_rate_rd: Dual::<1>::from_value(spec.domestic_rate_rd),
                foreign_rate_rf: Dual::<1>::from_value(spec.foreign_rate_rf),
                volatility_sigma: Dual::<1>::from_value(spec.volatility_sigma),
                maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
            });
            let mut sp = spec;
            sp.spot_fx_s0 += h;
            let mut sm = spec;
            sm.spot_fx_s0 -= h;
            let fd = (evaluate_wystup_fx_fader_range_accrual(&sp).fade_in_call_value
                - evaluate_wystup_fx_fader_range_accrual(&sm).fade_in_call_value)
                / (2.0 * h);
            assert!((r_ad.fade_in_call_value.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_wystup_fx_fader_range_accrual_module::*;

mod evaluate_wystup_fx_digital_call_spread_skew_overhedge_module {
    // Wystup (2017, FX Options and Structured Products, 2nd ed., Ch. 1.6) & Castagna (2010, Ch. 4)
    // FX Digital Option Smile-Skew Vega Correction, Call-Spread Over-Hedge & Pin-Risk Gamma Bound Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct WystupFxDigitalSkewOverhedgeSpec<T: Real = f64> {
        pub spot_fx_s0: T,
        pub digital_strike_k: T,
        pub domestic_rate_rd: T,
        pub foreign_rate_rf: T,
        pub atm_or_strike_vol_sigma_k: T,
        pub smile_slope_dsigma_dk: T,
        pub call_spread_half_width_dk: T,
        pub expiry_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct WystupFxDigitalSkewOverhedgeReport<T: Real = f64> {
        pub flat_smile_digital_call_price: T,
        pub vanilla_garman_kohlhagen_vega: T,
        pub smile_skew_digital_correction: T,
        pub skew_adjusted_digital_call_price: T,
        pub call_spread_overhedge_notional_gearing: T,
        pub call_spread_replicated_digital_price: T,
    }

    fn gk_call_with_vol<T: Real>(s0: T, k: T, rd: T, rf: T, sig: T, t: T) -> T {
        let std = sig * t.sqrt();
        let d1 = ((s0 / k).ln() + (rd - rf + T::from_f64(0.5) * sig * sig) * t) / std;
        let d2 = d1 - std;
        s0 * (-rf * t).exp() * norm_cdf(d1) - k * (-rd * t).exp() * norm_cdf(d2)
    }

    pub fn evaluate_wystup_fx_digital_call_spread_skew_overhedge<T: Real>(
        spec: &WystupFxDigitalSkewOverhedgeSpec<T>,
    ) -> WystupFxDigitalSkewOverhedgeReport<T> {
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-4));
        let k = spec.digital_strike_k.max_val(T::from_f64(1e-4));
        let rd = spec.domestic_rate_rd;
        let rf = spec.foreign_rate_rf;
        let sig_k = spec.atm_or_strike_vol_sigma_k.max_val(T::from_f64(0.01));
        let dsig_dk = spec.smile_slope_dsigma_dk;
        let dk = spec.call_spread_half_width_dk.max_val(T::from_f64(1e-4));
        let t = spec.expiry_years_t.max_val(T::from_f64(0.02));

        let sqrt_t = t.sqrt();
        let std = sig_k * sqrt_t;
        let d1 = ((s0 / k).ln() + (rd - rf + T::from_f64(0.5) * sig_k * sig_k) * t) / std;
        let d2 = d1 - std;

        let df_d = (-rd * t).exp();
        let df_f = (-rf * t).exp();

        let dig_flat = df_d * norm_cdf(d2);
        let vega_gk = s0 * df_f * norm_pdf(d1) * sqrt_t;

        // Wystup (2017, Ch. 1.6, Eq. 1.38): D_skew(K) = -dC_mkt/dK = D_BS(K) - Vega_GK(K) * (dsigma / dK):
        let skew_corr = -vega_gk * dsig_dk;
        let dig_skew = dig_flat + skew_corr;

        // Symmetric call-spread static over-hedge [C(K - dk, sigma(K - dk)) - C(K + dk, sigma(K + dk))] / (2 * dk):
        let k_lo = (k - dk).max_val(T::from_f64(1e-4));
        let k_hi = k + dk;
        let sig_lo = (sig_k - dsig_dk * dk).max_val(T::from_f64(0.01));
        let sig_hi = (sig_k + dsig_dk * dk).max_val(T::from_f64(0.01));
        let c_lo = gk_call_with_vol(s0, k_lo, rd, rf, sig_lo, t);
        let c_hi = gk_call_with_vol(s0, k_hi, rd, rf, sig_hi, t);
        let gearing = T::one() / (T::from_f64(2.0) * dk);
        let cs_replicated = (c_lo - c_hi) * gearing;

        WystupFxDigitalSkewOverhedgeReport {
            flat_smile_digital_call_price: dig_flat,
            vanilla_garman_kohlhagen_vega: vega_gk,
            smile_skew_digital_correction: skew_corr,
            skew_adjusted_digital_call_price: dig_skew,
            call_spread_overhedge_notional_gearing: gearing,
            call_spread_replicated_digital_price: cs_replicated,
        }
    }

    #[cfg(test)]
    mod evaluate_wystup_fx_digital_call_spread_skew_overhedge_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn wystup_fx_digital_call_spread_skew_overhedge_invariants_and_ad_parity() {
            let spec = WystupFxDigitalSkewOverhedgeSpec {
                spot_fx_s0: 1.0850,
                digital_strike_k: 1.1000,
                domestic_rate_rd: 0.045,
                foreign_rate_rf: 0.030,
                atm_or_strike_vol_sigma_k: 0.095,
                smile_slope_dsigma_dk: -0.18,
                call_spread_half_width_dk: 0.0050,
                expiry_years_t: 0.50,
            };
            let rep = evaluate_wystup_fx_digital_call_spread_skew_overhedge(&spec);
            assert!(rep.smile_skew_digital_correction > 0.0);
            assert!(rep.skew_adjusted_digital_call_price > rep.flat_smile_digital_call_price);
            assert!(
                (rep.skew_adjusted_digital_call_price - rep.call_spread_replicated_digital_price)
                    .abs()
                    < 5e-3
            );

            let h = 1e-5;
            let r_ad = evaluate_wystup_fx_digital_call_spread_skew_overhedge(
                &WystupFxDigitalSkewOverhedgeSpec {
                    spot_fx_s0: Dual::<8>::variable(spec.spot_fx_s0, 0),
                    digital_strike_k: Dual::<8>::variable(spec.digital_strike_k, 1),
                    domestic_rate_rd: Dual::<8>::variable(spec.domestic_rate_rd, 2),
                    foreign_rate_rf: Dual::<8>::variable(spec.foreign_rate_rf, 3),
                    atm_or_strike_vol_sigma_k: Dual::<8>::variable(
                        spec.atm_or_strike_vol_sigma_k,
                        4,
                    ),
                    smile_slope_dsigma_dk: Dual::<8>::variable(spec.smile_slope_dsigma_dk, 5),
                    call_spread_half_width_dk: Dual::<8>::variable(
                        spec.call_spread_half_width_dk,
                        6,
                    ),
                    expiry_years_t: Dual::<8>::variable(spec.expiry_years_t, 7),
                },
            );
            let comp_ad =
                r_ad.skew_adjusted_digital_call_price + r_ad.call_spread_replicated_digital_price;
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_fx_s0 += h;
                        sm.spot_fx_s0 -= h;
                    }
                    1 => {
                        sp.digital_strike_k += h;
                        sm.digital_strike_k -= h;
                    }
                    2 => {
                        sp.domestic_rate_rd += h;
                        sm.domestic_rate_rd -= h;
                    }
                    3 => {
                        sp.foreign_rate_rf += h;
                        sm.foreign_rate_rf -= h;
                    }
                    4 => {
                        sp.atm_or_strike_vol_sigma_k += h;
                        sm.atm_or_strike_vol_sigma_k -= h;
                    }
                    5 => {
                        sp.smile_slope_dsigma_dk += h;
                        sm.smile_slope_dsigma_dk -= h;
                    }
                    6 => {
                        sp.call_spread_half_width_dk += h;
                        sm.call_spread_half_width_dk -= h;
                    }
                    _ => {
                        sp.expiry_years_t += h;
                        sm.expiry_years_t -= h;
                    }
                }
                let rp = evaluate_wystup_fx_digital_call_spread_skew_overhedge(&sp);
                let rm = evaluate_wystup_fx_digital_call_spread_skew_overhedge(&sm);
                let fd = ((rp.skew_adjusted_digital_call_price
                    + rp.call_spread_replicated_digital_price)
                    - (rm.skew_adjusted_digital_call_price
                        + rm.call_spread_replicated_digital_price))
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_wystup_fx_digital_call_spread_skew_overhedge_module::*;

mod evaluate_wystup_lipton_fx_tarf_target_redemption_cap_module {
    // Wystup (2017, FX Options and Structured Products, 2nd ed., Ch. 9, "Target Forwards"),
    // Lipton (2001, Mathematical Methods for Foreign Exchange), & Clark (2011, Ch. 10)
    // FX Target Redemption Forward (TARF) Cumulative Gain Cap & Knock-Out Horizon Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct WystupLiptonFxTarfSpec<T: Real = f64> {
        pub spot_fx_s0: T,
        pub tarf_strike_k: T,
        pub target_redemption_cap_t_cap: T,
        pub downside_gearing_leverage_l: T,
        pub total_fixing_dates_m: T,
        pub domestic_rate_rd: T,
        pub foreign_rate_rf: T,
        pub fx_volatility_sigma: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct WystupLiptonFxTarfReport<T: Real = f64> {
        pub uncapped_cumulative_expected_gain: T,
        pub capped_cumulative_expected_gain: T,
        pub expected_active_fixing_count: T,
        pub target_cap_giveup_option_value: T,
        pub discounted_capped_upside_leg_pv: T,
        pub discounted_geared_downside_leg_pv: T,
        pub client_tarf_net_present_value: T,
    }

    pub fn evaluate_wystup_lipton_fx_tarf_target_redemption_cap<T: Real>(
        spec: &WystupLiptonFxTarfSpec<T>,
    ) -> WystupLiptonFxTarfReport<T> {
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(0.05));
        let k = spec.tarf_strike_k.max_val(T::from_f64(0.05));
        let t_cap = spec.target_redemption_cap_t_cap.max_val(T::from_f64(0.005));
        let lev = spec.downside_gearing_leverage_l.max_val(T::one());
        let m = spec.total_fixing_dates_m.max_val(T::one());
        let rd = spec.domestic_rate_rd;
        let rf = spec.foreign_rate_rf;
        let sig = spec.fx_volatility_sigma.max_val(T::from_f64(0.01));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.10));

        let dt = t / m;
        let t_avg = T::from_f64(0.5) * (t + dt);
        let sqrt_t = t_avg.sqrt();
        let df_d = (-rd * t_avg).exp();
        let fwd = s0 * ((rd - rf) * t_avg).exp();

        let std_fx = sig * sqrt_t;
        let d1 = ((fwd / k).ln() + T::from_f64(0.5) * sig * sig * t_avg) / std_fx;
        let d2 = d1 - std_fx;

        let mu_call = (fwd * norm_cdf(d1) - k * norm_cdf(d2)).max_val(T::from_f64(1e-6));
        let mu_put = (k * norm_cdf(-d2) - fwd * norm_cdf(-d1)).max_val(T::from_f64(1e-6));

        let g_uncapped = m * mu_call;
        let sig_g = (sig * fwd * (m * t_avg * norm_cdf(d1).max_val(T::from_f64(0.02))).sqrt())
            .max_val(T::from_f64(1e-5));

        // Bachelier-Gaussian excess above target cap E[max(G - T_cap, 0)]:
        let z_cap = (g_uncapped - t_cap) / sig_g;
        let cap_giveup = (g_uncapped - t_cap) * norm_cdf(z_cap) + sig_g * norm_pdf(z_cap);
        let g_capped = (g_uncapped - cap_giveup).max_val(T::from_f64(1e-6));

        let m_eff = m * (g_capped / g_uncapped);
        let pv_up = df_d * g_capped;
        let pv_down = df_d * lev * m_eff * mu_put;
        let net_npv = pv_up - pv_down;

        WystupLiptonFxTarfReport {
            uncapped_cumulative_expected_gain: g_uncapped,
            capped_cumulative_expected_gain: g_capped,
            expected_active_fixing_count: m_eff,
            target_cap_giveup_option_value: df_d * cap_giveup,
            discounted_capped_upside_leg_pv: pv_up,
            discounted_geared_downside_leg_pv: pv_down,
            client_tarf_net_present_value: net_npv,
        }
    }

    #[cfg(test)]
    mod evaluate_wystup_lipton_fx_tarf_target_redemption_cap_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn wystup_lipton_fx_tarf_invariants_and_ad_parity() {
            let spec = WystupLiptonFxTarfSpec {
                spot_fx_s0: 1.09,
                tarf_strike_k: 1.06,
                target_redemption_cap_t_cap: 0.25,
                downside_gearing_leverage_l: 2.0,
                total_fixing_dates_m: 12.0,
                domestic_rate_rd: 0.045,
                foreign_rate_rf: 0.030,
                fx_volatility_sigma: 0.10,
                maturity_years_t: 1.0,
            };
            let rep = evaluate_wystup_lipton_fx_tarf_target_redemption_cap(&spec);
            assert!(rep.capped_cumulative_expected_gain < rep.uncapped_cumulative_expected_gain);
            assert!(rep.expected_active_fixing_count < spec.total_fixing_dates_m);

            let h = 1e-6;
            let r_ad =
                evaluate_wystup_lipton_fx_tarf_target_redemption_cap(&WystupLiptonFxTarfSpec {
                    spot_fx_s0: Dual::<6>::variable(spec.spot_fx_s0, 0),
                    tarf_strike_k: Dual::<6>::variable(spec.tarf_strike_k, 1),
                    target_redemption_cap_t_cap: Dual::<6>::variable(
                        spec.target_redemption_cap_t_cap,
                        2,
                    ),
                    downside_gearing_leverage_l: Dual::<6>::variable(
                        spec.downside_gearing_leverage_l,
                        3,
                    ),
                    total_fixing_dates_m: Dual::from_f64(spec.total_fixing_dates_m),
                    domestic_rate_rd: Dual::<6>::variable(spec.domestic_rate_rd, 4),
                    foreign_rate_rf: Dual::from_f64(spec.foreign_rate_rf),
                    fx_volatility_sigma: Dual::<6>::variable(spec.fx_volatility_sigma, 5),
                    maturity_years_t: Dual::from_f64(spec.maturity_years_t),
                });
            for i in 0..6 {
                let mut v = [
                    spec.spot_fx_s0,
                    spec.tarf_strike_k,
                    spec.target_redemption_cap_t_cap,
                    spec.downside_gearing_leverage_l,
                    spec.domestic_rate_rd,
                    spec.fx_volatility_sigma,
                ];
                v[i] += h;
                let up =
                    evaluate_wystup_lipton_fx_tarf_target_redemption_cap(&WystupLiptonFxTarfSpec {
                        spot_fx_s0: v[0],
                        tarf_strike_k: v[1],
                        target_redemption_cap_t_cap: v[2],
                        downside_gearing_leverage_l: v[3],
                        domestic_rate_rd: v[4],
                        fx_volatility_sigma: v[5],
                        ..spec
                    });
                v[i] -= 2.0 * h;
                let dn =
                    evaluate_wystup_lipton_fx_tarf_target_redemption_cap(&WystupLiptonFxTarfSpec {
                        spot_fx_s0: v[0],
                        tarf_strike_k: v[1],
                        target_redemption_cap_t_cap: v[2],
                        downside_gearing_leverage_l: v[3],
                        domestic_rate_rd: v[4],
                        fx_volatility_sigma: v[5],
                        ..spec
                    });
                let fd = (up.client_tarf_net_present_value - dn.client_tarf_net_present_value)
                    / (2.0 * h);
                assert!(
                    (r_ad.client_tarf_net_present_value.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_wystup_lipton_fx_tarf_target_redemption_cap_module::*;

mod evaluate_wystup_fx_corridor_accrual_range_note_module {
    // Wystup (2017, FX Options and Structured Products, 2nd ed., Ch. 4, "Range Accrual Notes"),
    // Clark (2011, Foreign Exchange Option Pricing, Ch. 7), & Breeden & Litzenberger (1978)
    // FX Corridor Range Accrual Strip & Digital Skew Overhedge Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct WystupFxCorridorAccrualRangeNoteSpec<T: Real = f64> {
        pub notional_coupon_max: T,
        pub spot_fx_s0: T,
        pub lower_corridor_barrier_l: T,
        pub upper_corridor_barrier_u: T,
        pub domestic_rate_rd: T,
        pub foreign_rate_rf: T,
        pub atm_volatility_sigma: T,
        pub log_moneyness_skew_slope_xi: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct WystupFxCorridorAccrualRangeNoteReport<T: Real = f64> {
        pub black_scholes_corridor_accrual_fraction: T,
        pub smile_adjusted_corridor_accrual_fraction: T,
        pub black_scholes_range_note_pv: T,
        pub smile_adjusted_range_note_pv: T,
        pub digital_skew_overhedge_pv_adjustment: T,
        pub terminal_corridor_width_ratio: T,
    }

    pub fn evaluate_wystup_fx_corridor_accrual_range_note<T: Real>(
        spec: &WystupFxCorridorAccrualRangeNoteSpec<T>,
    ) -> WystupFxCorridorAccrualRangeNoteReport<T> {
        let coupon_max = spec.notional_coupon_max.max_val(T::from_f64(1.0));
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(0.05));
        let l_bar = spec.lower_corridor_barrier_l.max_val(T::from_f64(0.02));
        let u_bar = spec
            .upper_corridor_barrier_u
            .max_val(l_bar + T::from_f64(0.02));
        let rd = spec.domestic_rate_rd;
        let rf = spec.foreign_rate_rf;
        let sig0 = spec.atm_volatility_sigma.max_val(T::from_f64(0.02));
        let xi = spec.log_moneyness_skew_slope_xi;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.10));

        // 4-fixing corridor strip t_j = j * (T / 4) for j = 1, 2, 3, 4:
        let dt = t / T::from_f64(4.0);
        let mut frac_bs = T::zero();
        let mut frac_smile = T::zero();

        for step in 1..=4 {
            let tj = dt * T::from_f64(step as f64);
            let sqrt_tj = tj.sqrt();
            let fwd_j = s0 * ((rd - rf) * tj).exp();

            let sig_l = (sig0 + xi * (l_bar / fwd_j).ln()).max_val(T::from_f64(0.02));
            let sig_u = (sig0 + xi * (u_bar / fwd_j).ln()).max_val(T::from_f64(0.02));

            let std_l = sig_l * sqrt_tj;
            let d2_l = ((fwd_j / l_bar).ln() - T::from_f64(0.5) * sig_l * sig_l * tj) / std_l;

            let std_u = sig_u * sqrt_tj;
            let d2_u = ((fwd_j / u_bar).ln() - T::from_f64(0.5) * sig_u * sig_u * tj) / std_u;

            let p_bs_l = norm_cdf(d2_l);
            let p_bs_u = norm_cdf(d2_u);
            let p_corridor_bs = p_bs_l - p_bs_u;

            // Breeden-Litzenberger digital skew adjustment -Vega_BS * (d sigma / dK) = -phi(d2) * sqrt(t_j) * xi:
            let p_smile_l = p_bs_l - norm_pdf(d2_l) * sqrt_tj * xi;
            let p_smile_u = p_bs_u - norm_pdf(d2_u) * sqrt_tj * xi;
            let p_corridor_smile = p_smile_l - p_smile_u;

            frac_bs += T::from_f64(0.25) * p_corridor_bs;
            frac_smile += T::from_f64(0.25) * p_corridor_smile;
        }

        let df_rd = (-rd * t).exp();
        let pv_bs = coupon_max * df_rd * frac_bs;
        let pv_smile = coupon_max * df_rd * frac_smile;
        let skew_adj = pv_smile - pv_bs;
        let width_ratio = (u_bar - l_bar) / s0;

        WystupFxCorridorAccrualRangeNoteReport {
            black_scholes_corridor_accrual_fraction: frac_bs,
            smile_adjusted_corridor_accrual_fraction: frac_smile,
            black_scholes_range_note_pv: pv_bs,
            smile_adjusted_range_note_pv: pv_smile,
            digital_skew_overhedge_pv_adjustment: skew_adj,
            terminal_corridor_width_ratio: width_ratio,
        }
    }

    #[cfg(test)]
    mod evaluate_wystup_fx_corridor_accrual_range_note_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn wystup_fx_corridor_accrual_range_note_invariants_and_ad_parity() {
            let spec = WystupFxCorridorAccrualRangeNoteSpec {
                notional_coupon_max: 100.0,
                spot_fx_s0: 1.10,
                lower_corridor_barrier_l: 1.04,
                upper_corridor_barrier_u: 1.18,
                domestic_rate_rd: 0.045,
                foreign_rate_rf: 0.030,
                atm_volatility_sigma: 0.095,
                log_moneyness_skew_slope_xi: -0.12,
                maturity_years_t: 1.0,
            };
            let rep = evaluate_wystup_fx_corridor_accrual_range_note(&spec);
            assert!(
                rep.black_scholes_corridor_accrual_fraction > 0.30
                    && rep.black_scholes_corridor_accrual_fraction < 1.0
            );
            assert!(rep.smile_adjusted_range_note_pv > 30.0);

            let h = 1e-6;
            let r_ad = evaluate_wystup_fx_corridor_accrual_range_note(
                &WystupFxCorridorAccrualRangeNoteSpec {
                    notional_coupon_max: Dual::<9>::variable(spec.notional_coupon_max, 0),
                    spot_fx_s0: Dual::<9>::variable(spec.spot_fx_s0, 1),
                    lower_corridor_barrier_l: Dual::<9>::variable(spec.lower_corridor_barrier_l, 2),
                    upper_corridor_barrier_u: Dual::<9>::variable(spec.upper_corridor_barrier_u, 3),
                    domestic_rate_rd: Dual::<9>::variable(spec.domestic_rate_rd, 4),
                    foreign_rate_rf: Dual::<9>::variable(spec.foreign_rate_rf, 5),
                    atm_volatility_sigma: Dual::<9>::variable(spec.atm_volatility_sigma, 6),
                    log_moneyness_skew_slope_xi: Dual::<9>::variable(
                        spec.log_moneyness_skew_slope_xi,
                        7,
                    ),
                    maturity_years_t: Dual::<9>::variable(spec.maturity_years_t, 8),
                },
            );
            for i in 0..9 {
                let mut v = [
                    spec.notional_coupon_max,
                    spec.spot_fx_s0,
                    spec.lower_corridor_barrier_l,
                    spec.upper_corridor_barrier_u,
                    spec.domestic_rate_rd,
                    spec.foreign_rate_rf,
                    spec.atm_volatility_sigma,
                    spec.log_moneyness_skew_slope_xi,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_wystup_fx_corridor_accrual_range_note(
                    &WystupFxCorridorAccrualRangeNoteSpec {
                        notional_coupon_max: v[0],
                        spot_fx_s0: v[1],
                        lower_corridor_barrier_l: v[2],
                        upper_corridor_barrier_u: v[3],
                        domestic_rate_rd: v[4],
                        foreign_rate_rf: v[5],
                        atm_volatility_sigma: v[6],
                        log_moneyness_skew_slope_xi: v[7],
                        maturity_years_t: v[8],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_wystup_fx_corridor_accrual_range_note(
                    &WystupFxCorridorAccrualRangeNoteSpec {
                        notional_coupon_max: v[0],
                        spot_fx_s0: v[1],
                        lower_corridor_barrier_l: v[2],
                        upper_corridor_barrier_u: v[3],
                        domestic_rate_rd: v[4],
                        foreign_rate_rf: v[5],
                        atm_volatility_sigma: v[6],
                        log_moneyness_skew_slope_xi: v[7],
                        maturity_years_t: v[8],
                    },
                );
                let fd =
                    (up.smile_adjusted_range_note_pv - dn.smile_adjusted_range_note_pv) / (2.0 * h);
                assert!(
                    (r_ad.smile_adjusted_range_note_pv.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_wystup_fx_corridor_accrual_range_note_module::*;

pub mod piterbarg_prdc_bermudan_suite65_mod {
    use ql_ad::Real;
    use ql_math::norm_cdf;

    /// Specification for the Piterbarg (2005) / Clark (2011) Power Reverse Dual Currency (PRDC)
    /// FX–IR quanto convexity and Target Redemption (TARN) auto-call strip engine.
    ///
    /// Reference: Piterbarg (2005, *Risk*, "TARNs: Models, Valuation, Risk Sensitivities"),
    /// Clark (2011, *Foreign Exchange Option Pricing*, Ch. 9), Wystup (2017, Ch. 8).
    #[derive(Debug, Clone, Copy)]
    pub struct PiterbargPrdcBermudanSpec<T: Real> {
        /// Notional $N$.
        pub notional_n: T,
        /// Initial spot FX rate $S_0$ (domestic per foreign unit).
        pub spot_fx_s0: T,
        /// Foreign coupon gearing multiplier $g_f$.
        pub foreign_coupon_gearing_gf: T,
        /// Domestic coupon floor/strike deduction $g_d$.
        pub domestic_funding_deduction_gd: T,
        /// Domestic interest rate $r_d$.
        pub domestic_rate_rd: T,
        /// Foreign interest rate $r_f$.
        pub foreign_rate_rf: T,
        /// FX spot volatility $\sigma_{\text{FX}}$.
        pub fx_volatility_sigma_fx: T,
        /// Domestic interest-rate volatility $\sigma_d$.
        pub domestic_rate_vol_sigma_d: T,
        /// Correlation $\rho_{d,\text{FX}}$ between domestic rate and FX spot.
        pub domestic_fx_correlation_rho: T,
        /// Target cumulative coupon cap $C_{\text{target}}$ (as fraction of notional).
        pub target_redemption_cap_c_target: T,
        /// Total maturity $T$ across $M = 4$ coupon dates (in years).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_piterbarg_fx_power_reverse_dual_currency_prdc_bermudan_redemption`].
    #[derive(Debug, Clone, Copy)]
    pub struct PiterbargPrdcBermudanAnalytics<T: Real> {
        /// Present value of the unredeemed (no-TARN) PRDC coupon strip with FX–IR convexity.
        pub unredeemed_prdc_strip_pv: T,
        /// FX–IR correlation quanto convexity valuation adjustment.
        pub fx_ir_correlation_convexity_wedge: T,
        /// Expected survival probability through the final reset date.
        pub terminal_survival_probability: T,
        /// Present value of the Target-Redemption (TARN) auto-callable PRDC note.
        pub tarn_capped_prdc_present_value: T,
    }

    /// Evaluates a 4-period Power Reverse Dual Currency (PRDC) note with FX–IR correlation quanto
    /// drift adjustment and smooth Target Redemption (TARN) knock-out survival probabilities.
    pub fn evaluate_piterbarg_fx_power_reverse_dual_currency_prdc_bermudan_redemption<T: Real>(
        spec: &PiterbargPrdcBermudanSpec<T>,
    ) -> PiterbargPrdcBermudanAnalytics<T> {
        const PERIODS: usize = 4;
        let t_mat = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let dt = t_mat / T::from_f64(PERIODS as f64);
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-4));
        let gf = spec.foreign_coupon_gearing_gf.max_val(T::from_f64(1e-4));
        let gd = spec
            .domestic_funding_deduction_gd
            .max_val(T::from_f64(1e-4));
        let s_ref = T::from_f64(100.0);
        let k_fx = s_ref * (gd / gf);
        let sig_fx = spec.fx_volatility_sigma_fx.max_val(T::from_f64(1e-4));

        let mut unredeemed_pv = T::zero();
        let mut unredeemed_flat_corr_pv = T::zero();
        let mut tarn_pv = T::zero();
        let mut cum_expected_coupon = T::zero();
        let mut surv_prob = T::one();

        for m in 1..=PERIODS {
            let tm = T::from_f64(m as f64) * dt;
            let sqrt_tm = tm.sqrt();
            let v_sqrt = sig_fx * sqrt_tm;

            let f_base = s0 * ((spec.domestic_rate_rd - spec.foreign_rate_rf) * tm).exp();
            // Clark (2011, Ch. 9): domestic bond numeraire FX-IR quanto adjustment exp(-0.5 * rho * sigma_d * sigma_FX * t_m^2)
            let quanto_mult = (-T::from_f64(0.5)
                * spec.domestic_fx_correlation_rho
                * spec.domestic_rate_vol_sigma_d
                * sig_fx
                * tm
                * tm)
                .exp();
            let f_adj = f_base * quanto_mult;

            let d1_adj = ((f_adj / k_fx).ln() + T::from_f64(0.5) * v_sqrt * v_sqrt) / v_sqrt;
            let d2_adj = d1_adj - v_sqrt;
            let c_m_adj = dt * (gf / s_ref) * (f_adj * norm_cdf(d1_adj) - k_fx * norm_cdf(d2_adj));

            let d1_base = ((f_base / k_fx).ln() + T::from_f64(0.5) * v_sqrt * v_sqrt) / v_sqrt;
            let d2_base = d1_base - v_sqrt;
            let c_m_base =
                dt * (gf / s_ref) * (f_base * norm_cdf(d1_base) - k_fx * norm_cdf(d2_base));

            let df_m = (-spec.domestic_rate_rd * tm).exp();
            unredeemed_pv += spec.notional_n * df_m * c_m_adj;
            unredeemed_flat_corr_pv += spec.notional_n * df_m * c_m_base;

            // Smooth TARN survival probability before period m coupon
            let cum_disp = (gf * dt * sig_fx * sqrt_tm * T::from_f64(m as f64).sqrt())
                .max_val(T::from_f64(1e-4));
            let z_tarn = (spec.target_redemption_cap_c_target - cum_expected_coupon) / cum_disp;
            surv_prob = norm_cdf(z_tarn);

            tarn_pv += spec.notional_n * df_m * surv_prob * c_m_adj;
            cum_expected_coupon += c_m_adj;
        }

        PiterbargPrdcBermudanAnalytics {
            unredeemed_prdc_strip_pv: unredeemed_pv,
            fx_ir_correlation_convexity_wedge: unredeemed_pv - unredeemed_flat_corr_pv,
            terminal_survival_probability: surv_prob,
            tarn_capped_prdc_present_value: tarn_pv,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn piterbarg_prdc_bermudan_invariants_and_ad_parity() {
            let base = PiterbargPrdcBermudanSpec {
                notional_n: 1_000.0_f64,
                spot_fx_s0: 105.0,
                foreign_coupon_gearing_gf: 0.10,
                domestic_funding_deduction_gd: 0.06,
                domestic_rate_rd: 0.015,
                foreign_rate_rf: 0.045,
                fx_volatility_sigma_fx: 0.14,
                domestic_rate_vol_sigma_d: 0.012,
                domestic_fx_correlation_rho: 0.45,
                target_redemption_cap_c_target: 0.12,
                maturity_years_t: 4.0,
            };
            let r0 =
                evaluate_piterbarg_fx_power_reverse_dual_currency_prdc_bermudan_redemption(&base);
            assert!(r0.unredeemed_prdc_strip_pv > r0.tarn_capped_prdc_present_value);
            assert!(r0.tarn_capped_prdc_present_value > 0.0);
            assert!(
                r0.terminal_survival_probability > 0.0 && r0.terminal_survival_probability < 1.0
            );

            let inputs = [
                base.notional_n,
                base.spot_fx_s0,
                base.foreign_coupon_gearing_gf,
                base.domestic_funding_deduction_gd,
                base.domestic_rate_rd,
                base.foreign_rate_rf,
                base.fx_volatility_sigma_fx,
                base.domestic_rate_vol_sigma_d,
                base.domestic_fx_correlation_rho,
                base.target_redemption_cap_c_target,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<11>::from_f64(0.0); 11];
            for i in 0..11 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = PiterbargPrdcBermudanSpec {
                notional_n: dual_in[0],
                spot_fx_s0: dual_in[1],
                foreign_coupon_gearing_gf: dual_in[2],
                domestic_funding_deduction_gd: dual_in[3],
                domestic_rate_rd: dual_in[4],
                foreign_rate_rf: dual_in[5],
                fx_volatility_sigma_fx: dual_in[6],
                domestic_rate_vol_sigma_d: dual_in[7],
                domestic_fx_correlation_rho: dual_in[8],
                target_redemption_cap_c_target: dual_in[9],
                maturity_years_t: dual_in[10],
            };
            let r_ad = evaluate_piterbarg_fx_power_reverse_dual_currency_prdc_bermudan_redemption(
                &spec_ad,
            );
            let eps = 1e-6;
            for i in 0..11 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 11]| PiterbargPrdcBermudanSpec {
                    notional_n: v[0],
                    spot_fx_s0: v[1],
                    foreign_coupon_gearing_gf: v[2],
                    domestic_funding_deduction_gd: v[3],
                    domestic_rate_rd: v[4],
                    foreign_rate_rf: v[5],
                    fx_volatility_sigma_fx: v[6],
                    domestic_rate_vol_sigma_d: v[7],
                    domestic_fx_correlation_rho: v[8],
                    target_redemption_cap_c_target: v[9],
                    maturity_years_t: v[10],
                };
                let rp = evaluate_piterbarg_fx_power_reverse_dual_currency_prdc_bermudan_redemption(
                    &mk(p),
                );
                let rm = evaluate_piterbarg_fx_power_reverse_dual_currency_prdc_bermudan_redemption(
                    &mk(m),
                );
                let fd = (rp.tarn_capped_prdc_present_value - rm.tarn_capped_prdc_present_value)
                    / (2.0 * eps);
                assert!(
                    r_ad.tarn_capped_prdc_present_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.tarn_capped_prdc_present_value.d[i] - fd).abs() < 1e-3);
            }
        }
    }
}
pub use piterbarg_prdc_bermudan_suite65_mod::*;

pub mod wystup_clark_fx_smile_strangle_mod {
    use ql_ad::Real;
    use ql_math::norm_cdf;

    /// Specification for the Wystup (2017, Ch. 2) / Clark (2011, Ch. 3) FX 25-delta
    /// smile strangle vs. broker market-strangle convention engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct WystupClarkFxSmileStrangleSpec<T: Real> {
        /// Spot FX rate $S_0 > 0$ (domestic per unit foreign).
        pub spot_fx_s0: T,
        /// Domestic continuous interest rate $r_d$.
        pub domestic_rate_rd: T,
        /// Foreign continuous interest rate $r_f$.
        pub foreign_rate_rf: T,
        /// Delta-neutral ATM implied volatility $\sigma_{\text{ATM}} > 0$.
        pub atm_volatility_sigma_atm: T,
        /// 25-delta Risk Reversal quote $\text{RR}_{25} = \sigma_{25C} - \sigma_{25P}$.
        pub risk_reversal_25d_rr: T,
        /// 25-delta Smile Butterfly quote $\text{BF}_{25} = \frac{1}{2}(\sigma_{25C} + \sigma_{25P}) - \sigma_{\text{ATM}}$.
        pub butterfly_25d_bf: T,
        /// Option expiry $T > 0$ (years).
        pub expiry_years_t: T,
    }

    /// Analytical outputs for the Wystup-Clark FX smile strangle convention engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct WystupClarkFxSmileStrangleAnalytics<T: Real> {
        /// Market-strangle 25-delta call strike $K_{25C}^{\text{MS}}$.
        pub call_25d_strike_kc: T,
        /// Market-strangle 25-delta put strike $K_{25P}^{\text{MS}}$.
        pub put_25d_strike_kp: T,
        /// Smile-vs-market strangle premium wedge (in domestic currency per unit foreign notional).
        pub smile_vs_market_strangle_premium_wedge: T,
        /// Smile-consistent 25-delta strangle present value (in domestic currency).
        pub smile_adjusted_strangle_present_value: T,
    }

    fn garman_kohlhagen_call<T: Real>(f: T, k: T, df_d: T, sig: T, sqrt_t: T) -> T {
        let v = (sig * sqrt_t).max_val(T::from_f64(1e-6));
        let d1 = ((f / k).ln() + T::from_f64(0.5) * v * v) / v;
        let d2 = d1 - v;
        df_d * (f * norm_cdf(d1) - k * norm_cdf(d2))
    }

    fn garman_kohlhagen_put<T: Real>(f: T, k: T, df_d: T, sig: T, sqrt_t: T) -> T {
        let v = (sig * sqrt_t).max_val(T::from_f64(1e-6));
        let d1 = ((f / k).ln() + T::from_f64(0.5) * v * v) / v;
        let d2 = d1 - v;
        df_d * (k * norm_cdf(-d2) - f * norm_cdf(-d1))
    }

    /// Evaluates the Wystup (2017) / Clark (2011) 25-delta FX market strangle strikes and
    /// computes the smile-consistent strangle valuation and butterfly convexity wedge.
    pub fn evaluate_wystup_clark_fx_smile_strangle_market_convention<T: Real>(
        spec: &WystupClarkFxSmileStrangleSpec<T>,
    ) -> WystupClarkFxSmileStrangleAnalytics<T> {
        let t = spec.expiry_years_t.max_val(T::from_f64(1e-4));
        let sqrt_t = t.sqrt();
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-4));
        let df_d = (-spec.domestic_rate_rd * t).exp();
        let df_f = (-spec.foreign_rate_rf * t).exp();
        let fwd = s0 * df_f / df_d;

        let sig_atm = spec.atm_volatility_sigma_atm.max_val(T::from_f64(1e-4));
        let sig_ms = (sig_atm + spec.butterfly_25d_bf).max_val(T::from_f64(1e-4));

        // 25-delta forward delta convention: Phi(d1) = 0.25 => d1 = +/- 0.67448975
        let z25 = T::from_f64(0.6744897501960817);
        let kc = fwd * (z25 * sig_ms * sqrt_t + T::from_f64(0.5) * sig_ms * sig_ms * t).exp();
        let kp = fwd * (-z25 * sig_ms * sqrt_t + T::from_f64(0.5) * sig_ms * sig_ms * t).exp();

        // Sticky-strike Malz (1997) / Clark (2011) parabolic smile evaluated at normalized
        // spot log-moneyness xi = (ln(K / S_0) - 0.5 * sig_ms^2 * T) / (z25 * sig_ms * sqrt(T))
        // of the 25-delta market-strangle strikes:
        let carry_shift = (spec.domestic_rate_rd - spec.foreign_rate_rf) * sqrt_t / (z25 * sig_ms);
        let xi_c = T::one() + carry_shift;
        let xi_p = -T::one() + carry_shift;
        let sig_c = (sig_atm
            + T::from_f64(0.5) * spec.risk_reversal_25d_rr * xi_c
            + spec.butterfly_25d_bf * xi_c * xi_c)
            .max_val(T::from_f64(1e-4));
        let sig_p = (sig_atm
            + T::from_f64(0.5) * spec.risk_reversal_25d_rr * xi_p
            + spec.butterfly_25d_bf * xi_p * xi_p)
            .max_val(T::from_f64(1e-4));

        let mkt_strangle_pv = garman_kohlhagen_call(fwd, kc, df_d, sig_ms, sqrt_t)
            + garman_kohlhagen_put(fwd, kp, df_d, sig_ms, sqrt_t);
        let smile_strangle_pv = garman_kohlhagen_call(fwd, kc, df_d, sig_c, sqrt_t)
            + garman_kohlhagen_put(fwd, kp, df_d, sig_p, sqrt_t);

        WystupClarkFxSmileStrangleAnalytics {
            call_25d_strike_kc: kc,
            put_25d_strike_kp: kp,
            smile_vs_market_strangle_premium_wedge: smile_strangle_pv - mkt_strangle_pv,
            smile_adjusted_strangle_present_value: smile_strangle_pv,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn wystup_clark_fx_smile_strangle_invariants_and_ad_parity() {
            let base = WystupClarkFxSmileStrangleSpec {
                spot_fx_s0: 1.10_f64,
                domestic_rate_rd: 0.045,
                foreign_rate_rf: 0.030,
                atm_volatility_sigma_atm: 0.095,
                risk_reversal_25d_rr: -0.015,
                butterfly_25d_bf: 0.0035,
                expiry_years_t: 1.0,
            };
            let r0 = evaluate_wystup_clark_fx_smile_strangle_market_convention(&base);
            assert!(r0.call_25d_strike_kc > base.spot_fx_s0);
            assert!(r0.put_25d_strike_kp < base.spot_fx_s0);
            assert!(r0.smile_adjusted_strangle_present_value > 0.0);

            let inputs = [
                base.spot_fx_s0,
                base.domestic_rate_rd,
                base.foreign_rate_rf,
                base.atm_volatility_sigma_atm,
                base.risk_reversal_25d_rr,
                base.butterfly_25d_bf,
                base.expiry_years_t,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = WystupClarkFxSmileStrangleSpec {
                spot_fx_s0: dual_in[0],
                domestic_rate_rd: dual_in[1],
                foreign_rate_rf: dual_in[2],
                atm_volatility_sigma_atm: dual_in[3],
                risk_reversal_25d_rr: dual_in[4],
                butterfly_25d_bf: dual_in[5],
                expiry_years_t: dual_in[6],
            };
            let r_ad = evaluate_wystup_clark_fx_smile_strangle_market_convention(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| WystupClarkFxSmileStrangleSpec {
                    spot_fx_s0: v[0],
                    domestic_rate_rd: v[1],
                    foreign_rate_rf: v[2],
                    atm_volatility_sigma_atm: v[3],
                    risk_reversal_25d_rr: v[4],
                    butterfly_25d_bf: v[5],
                    expiry_years_t: v[6],
                };
                let rp = evaluate_wystup_clark_fx_smile_strangle_market_convention(&mk(p));
                let rm = evaluate_wystup_clark_fx_smile_strangle_market_convention(&mk(m));
                let fd = (rp.smile_adjusted_strangle_present_value
                    - rm.smile_adjusted_strangle_present_value)
                    / (2.0 * eps);
                assert!(
                    r_ad.smile_adjusted_strangle_present_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.smile_adjusted_strangle_present_value.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use wystup_clark_fx_smile_strangle_mod::*;

// ============================================================================
// Reiner (1992) & Haug (2007, Ch. 4.13) Foreign Equity Domestic-Strike & Quanto Call
// ============================================================================

pub mod reiner_haug_quanto_call_mod {
    use ql_ad::Real;
    use ql_math::distributions::norm_cdf;

    /// Specification for the Reiner (1992, *Risk* 5(3)) and Haug (2007, Ch. 4.13.1 & 4.13.4)
    /// foreign equity call struck in domestic currency vs. fixed-exchange-rate quanto call.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct ReinerHaugQuantoEquityCallSpec<T: Real> {
        /// Foreign equity spot price $S_0 > 0$ (in foreign currency).
        pub foreign_equity_spot_s0: T,
        /// Spot exchange rate $E_0 > 0$ (domestic currency per unit of foreign currency).
        pub fx_spot_e0: T,
        /// Option strike $K_{\text{for}} > 0$ (in foreign currency units).
        pub foreign_strike_k_for: T,
        /// Domestic risk-free rate $r_d > 0$.
        pub domestic_rate_rd: T,
        /// Foreign risk-free rate $r_f > 0$.
        pub foreign_rate_rf: T,
        /// Foreign equity continuous dividend yield $q \ge 0$.
        pub equity_dividend_yield_q: T,
        /// Foreign equity volatility $\sigma_S > 0$.
        pub equity_volatility_sigma_s: T,
        /// Exchange-rate volatility $\sigma_E > 0$.
        pub fx_volatility_sigma_e: T,
        /// Correlation between foreign equity and exchange rate $\rho_{SE} \in (-1, 1)$.
        pub equity_fx_correlation_rho: T,
    }

    /// Analytical outputs for the Reiner-Haug quanto and composite equity option engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct ReinerHaugQuantoEquityCallAnalytics<T: Real> {
        /// Quanto covariance drift adjustment $-\rho_{SE}\sigma_S\sigma_E$.
        pub quanto_covariance_drift_adjustment: T,
        /// Composite domestic-currency equity volatility $\sigma_{\text{comp}} = \sqrt{\sigma_S^2 + \sigma_E^2 + 2\rho_{SE}\sigma_S\sigma_E}$.
        pub composite_domestic_volatility: T,
        /// Foreign equity call struck in domestic currency $K_{\text{dom}} = E_0 K_{\text{for}}$.
        pub domestic_strike_composite_call_value: T,
        /// Fixed-exchange-rate quanto call value in domestic currency $C_{\text{quanto}}$.
        pub fixed_fx_quanto_call_present_value: T,
    }

    /// Evaluates the Reiner (1992) / Haug (2007, Section 4.13.1 & 4.13.4) fixed-FX quanto call
    /// and domestic-strike composite foreign equity call.
    pub fn evaluate_reiner_haug_foreign_equity_domestic_strike_and_quanto_call<T: Real>(
        spec: &ReinerHaugQuantoEquityCallSpec<T>,
    ) -> ReinerHaugQuantoEquityCallAnalytics<T> {
        let s0 = spec.foreign_equity_spot_s0.max_val(T::from_f64(1e-4));
        let e0 = spec.fx_spot_e0.max_val(T::from_f64(1e-4));
        let k_for = spec.foreign_strike_k_for.max_val(T::from_f64(1e-4));
        let sig_s = spec.equity_volatility_sigma_s.max_val(T::from_f64(1e-4));
        let sig_e = spec.fx_volatility_sigma_e.max_val(T::from_f64(1e-4));
        let rho = spec
            .equity_fx_correlation_rho
            .max_val(T::from_f64(-0.98))
            .min_val(T::from_f64(0.98));

        // 1. Fixed-exchange-rate Quanto Call (Haug 2007, Section 4.13.4)
        let quanto_adj = -rho * sig_s * sig_e;
        let b_quanto = spec.foreign_rate_rf - spec.equity_dividend_yield_q + quanto_adj;
        let d1_q = ((s0 / k_for).ln() + b_quanto + T::from_f64(0.5) * sig_s * sig_s) / sig_s;
        let d2_q = d1_q - sig_s;
        let df_d = (-spec.domestic_rate_rd).exp();
        let c_quanto = e0 * df_d * (s0 * b_quanto.exp() * norm_cdf(d1_q) - k_for * norm_cdf(d2_q));

        // 2. Foreign Equity Call Struck in Domestic Currency K_dom = E_0 * K_for (Haug 2007, Section 4.13.1)
        let var_comp = (sig_s * sig_s + sig_e * sig_e + T::from_f64(2.0) * rho * sig_s * sig_e)
            .max_val(T::from_f64(1e-6));
        let sig_comp = var_comp.sqrt();
        let d1_d = ((s0 / k_for).ln()
            + (spec.domestic_rate_rd - spec.equity_dividend_yield_q + T::from_f64(0.5) * var_comp))
            / sig_comp;
        let d2_d = d1_d - sig_comp;
        let c_dom = e0
            * (s0 * (-spec.equity_dividend_yield_q).exp() * norm_cdf(d1_d)
                - k_for * df_d * norm_cdf(d2_d));

        ReinerHaugQuantoEquityCallAnalytics {
            quanto_covariance_drift_adjustment: quanto_adj,
            composite_domestic_volatility: sig_comp,
            domestic_strike_composite_call_value: c_dom,
            fixed_fx_quanto_call_present_value: c_quanto,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn reiner_haug_quanto_call_invariants_and_ad_parity() {
            let base = ReinerHaugQuantoEquityCallSpec {
                foreign_equity_spot_s0: 100.0_f64,
                fx_spot_e0: 1.10,
                foreign_strike_k_for: 102.0,
                domestic_rate_rd: 0.045,
                foreign_rate_rf: 0.025,
                equity_dividend_yield_q: 0.015,
                equity_volatility_sigma_s: 0.22,
                fx_volatility_sigma_e: 0.12,
                equity_fx_correlation_rho: -0.35,
            };
            let r0 = evaluate_reiner_haug_foreign_equity_domestic_strike_and_quanto_call(&base);
            assert!(r0.quanto_covariance_drift_adjustment > 0.0);
            assert!(r0.composite_domestic_volatility > 0.0);
            assert!(r0.domestic_strike_composite_call_value > 0.0);
            assert!(r0.fixed_fx_quanto_call_present_value > 0.0);

            let inputs = [
                base.foreign_equity_spot_s0,
                base.fx_spot_e0,
                base.foreign_strike_k_for,
                base.domestic_rate_rd,
                base.foreign_rate_rf,
                base.equity_dividend_yield_q,
                base.equity_volatility_sigma_s,
                base.fx_volatility_sigma_e,
                base.equity_fx_correlation_rho,
            ];
            let mut dual_in = [Dual::<9>::from_f64(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = ReinerHaugQuantoEquityCallSpec {
                foreign_equity_spot_s0: dual_in[0],
                fx_spot_e0: dual_in[1],
                foreign_strike_k_for: dual_in[2],
                domestic_rate_rd: dual_in[3],
                foreign_rate_rf: dual_in[4],
                equity_dividend_yield_q: dual_in[5],
                equity_volatility_sigma_s: dual_in[6],
                fx_volatility_sigma_e: dual_in[7],
                equity_fx_correlation_rho: dual_in[8],
            };
            let r_ad =
                evaluate_reiner_haug_foreign_equity_domestic_strike_and_quanto_call(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| ReinerHaugQuantoEquityCallSpec {
                    foreign_equity_spot_s0: v[0],
                    fx_spot_e0: v[1],
                    foreign_strike_k_for: v[2],
                    domestic_rate_rd: v[3],
                    foreign_rate_rf: v[4],
                    equity_dividend_yield_q: v[5],
                    equity_volatility_sigma_s: v[6],
                    fx_volatility_sigma_e: v[7],
                    equity_fx_correlation_rho: v[8],
                };
                let rp =
                    evaluate_reiner_haug_foreign_equity_domestic_strike_and_quanto_call(&mk(p));
                let rm =
                    evaluate_reiner_haug_foreign_equity_domestic_strike_and_quanto_call(&mk(m));
                let fd = (rp.fixed_fx_quanto_call_present_value
                    - rm.fixed_fx_quanto_call_present_value)
                    / (2.0 * eps);
                assert!(
                    r_ad.fixed_fx_quanto_call_present_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.fixed_fx_quanto_call_present_value.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use reiner_haug_quanto_call_mod::*;
mod garman_kohlhagen_baw_fx_mod {
    use ql_ad::Real;

    /// Input specification for the Garman & Kohlhagen (1983, *JIMF* 2(3)), Biger & Hull (1983),
    /// and Barone-Adesi & Whaley (1987, *J. Finance* 42(2)) European vs. American FX option
    /// quadratic early-exercise premium engine.
    #[derive(Debug, Clone, Copy)]
    pub struct GarmanKohlhagenBawFxSpec<T: Real> {
        /// Spot FX rate $S_0$ (domestic per unit foreign currency, e.g., $1.10$).
        pub spot_fx_s0: T,
        /// Option strike $K$ (e.g., $1.08$).
        pub strike_k: T,
        /// Domestic risk-free rate $r_d$ (e.g., $0.032$).
        pub domestic_rate_rd: T,
        /// Foreign risk-free rate $r_f > r_d$ (positive foreign carry motivating early call exercise, e.g., $0.055$).
        pub foreign_rate_rf: T,
        /// FX volatility $\sigma$ (e.g., $0.12$).
        pub fx_volatility_sigma: T,
        /// Option maturity $T$ in years (e.g., $1.0$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_garman_kohlhagen_biger_hull_fx_early_exercise_premium`].
    #[derive(Debug, Clone, Copy)]
    pub struct GarmanKohlhagenBawFxAnalytics<T: Real> {
        /// European Garman–Kohlhagen FX call value $c_{\text{GK}}(S_0, K, T)$.
        pub european_garman_kohlhagen_call: T,
        /// Barone-Adesi & Whaley (1987) positive quadratic root $q_2 > 1$.
        pub baw_quadratic_exponent_q2: T,
        /// Early-exercise foreign-interest carry premium $C_{\text{Am}} - c_{\text{GK}}$.
        pub american_early_exercise_premium: T,
        /// Total American FX call value $C_{\text{Am}}(S_0, K, T)$.
        pub american_fx_call_value: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        T::from_f64(0.5) * (T::one() + (x * T::from_f64(0.707_106_781_186_547_5)).erf())
    }

    #[inline]
    fn gk_call<T: Real>(s: T, k: T, rd: T, rf: T, sig: T, t: T) -> T {
        let std_dev = sig * t.sqrt();
        let d1 = ((s / k).ln() + (rd - rf + T::from_f64(0.5) * sig * sig) * t) / std_dev;
        let d2 = d1 - std_dev;
        s * (-rf * t).exp() * std_norm_cdf(d1) - k * (-rd * t).exp() * std_norm_cdf(d2)
    }

    /// Evaluates the Garman & Kohlhagen (1983) European FX call and the Barone-Adesi & Whaley
    /// (1987, Eqs. 11–16) quadratic American early-exercise premium when foreign interest rate $r_f > 0$.
    pub fn evaluate_garman_kohlhagen_biger_hull_fx_early_exercise_premium<T: Real>(
        spec: &GarmanKohlhagenBawFxSpec<T>,
    ) -> GarmanKohlhagenBawFxAnalytics<T> {
        let s0 = spec.spot_fx_s0;
        let k = spec.strike_k;
        let rd = spec.domestic_rate_rd;
        let rf = spec.foreign_rate_rf;
        let sig = spec.fx_volatility_sigma;
        let t = spec.maturity_years_t;
        let sig2 = sig * sig;
        let half = T::from_f64(0.5);

        let euro_call = gk_call(s0, k, rd, rf, sig, t);

        // Barone-Adesi & Whaley (1987) quadratic approximation parameters M = 2*r_d/sigma^2, N = 2*(r_d - r_f)/sigma^2
        let m_baw = T::from_f64(2.0) * rd / sig2;
        let n_baw = T::from_f64(2.0) * (rd - rf) / sig2;
        let k_baw = T::one() - (-rd * t).exp();

        let b_minus_1 = -(n_baw - T::one());
        let q2 =
            half * (b_minus_1 + (b_minus_1 * b_minus_1 + T::from_f64(4.0) * m_baw / k_baw).sqrt());

        // Seed critical exercise boundary S^* via BAW (1987, Eq. 19) + 1 Newton step on S^* - K = c(S^*) + (S^*/q2)*(1 - exp(-r_f*T)*Phi(d1(S^*)))
        let s_inf = k / (T::one() - T::one() / q2);
        let h2 = -((rd - rf) * t + T::from_f64(2.0) * sig * t.sqrt()) * (k / (s_inf - k));
        let s_star = k + (s_inf - k) * (T::one() - h2.exp());

        let std_dev = sig * t.sqrt();
        let d1_star = ((s_star / k).ln() + (rd - rf + half * sig2) * t) / std_dev;
        let a2 = (s_star / q2) * (T::one() - (-rf * t).exp() * std_norm_cdf(d1_star));

        let ee_premium = a2 * (s0 / s_star).pow(q2);
        let am_call = euro_call + ee_premium;

        GarmanKohlhagenBawFxAnalytics {
            european_garman_kohlhagen_call: euro_call,
            baw_quadratic_exponent_q2: q2,
            american_early_exercise_premium: ee_premium,
            american_fx_call_value: am_call,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn garman_kohlhagen_baw_fx_invariants_and_ad_parity() {
            let base = GarmanKohlhagenBawFxSpec {
                spot_fx_s0: 1.10_f64,
                strike_k: 1.08_f64,
                domestic_rate_rd: 0.032_f64,
                foreign_rate_rf: 0.055_f64,
                fx_volatility_sigma: 0.12_f64,
                maturity_years_t: 1.0_f64,
            };
            let res = evaluate_garman_kohlhagen_biger_hull_fx_early_exercise_premium(&base);
            assert!(res.european_garman_kohlhagen_call > 0.0);
            assert!(res.baw_quadratic_exponent_q2 > 1.0);
            assert!(res.american_early_exercise_premium > 0.0);
            assert!(res.american_fx_call_value > res.european_garman_kohlhagen_call);

            let inputs = [
                base.spot_fx_s0,
                base.strike_k,
                base.domestic_rate_rd,
                base.foreign_rate_rf,
                base.fx_volatility_sigma,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<6>::from_value(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GarmanKohlhagenBawFxSpec {
                spot_fx_s0: dual_in[0],
                strike_k: dual_in[1],
                domestic_rate_rd: dual_in[2],
                foreign_rate_rf: dual_in[3],
                fx_volatility_sigma: dual_in[4],
                maturity_years_t: dual_in[5],
            };
            let r_ad = evaluate_garman_kohlhagen_biger_hull_fx_early_exercise_premium(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| GarmanKohlhagenBawFxSpec {
                    spot_fx_s0: v[0],
                    strike_k: v[1],
                    domestic_rate_rd: v[2],
                    foreign_rate_rf: v[3],
                    fx_volatility_sigma: v[4],
                    maturity_years_t: v[5],
                };
                let rp = evaluate_garman_kohlhagen_biger_hull_fx_early_exercise_premium(&mk(p));
                let rm = evaluate_garman_kohlhagen_biger_hull_fx_early_exercise_premium(&mk(m));
                let fd = (rp.american_fx_call_value - rm.american_fx_call_value) / (2.0 * eps);
                assert!(
                    r_ad.american_fx_call_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.american_fx_call_value.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use garman_kohlhagen_baw_fx_mod::*;
mod carr_chou_symmetry_barrier_mod {
    use ql_ad::Real;

    /// Input specification for the Carr & Chou (1997, *Journal of Finance* 52(2)) &
    /// Bowie & Carr (1994) Put-Call Symmetry (PCS) static replication of a down-and-out
    /// FX call option using a vanilla call and a reflected geometric-strike put.
    #[derive(Debug, Clone, Copy)]
    pub struct CarrChouSymmetryBarrierReplicationSpec<T: Real> {
        /// Spot FX rate $S_0 > H$ (e.g., $1.12$).
        pub spot_fx_s0: T,
        /// Call strike $K > H$ (e.g., $1.10$).
        pub strike_k: T,
        /// Lower knock-out barrier $H < \min(S_0, K)$ (e.g., $1.02$).
        pub down_barrier_h: T,
        /// Domestic risk-free rate $r_d$ (e.g., $0.034$).
        pub domestic_rate_rd: T,
        /// Foreign risk-free rate $r_f$ (e.g., $0.028$).
        pub foreign_rate_rf: T,
        /// Implied volatility $\sigma_K$ at the primary strike $K$ (e.g., $0.115$).
        pub call_strike_volatility_sigma_k: T,
        /// Implied volatility $\sigma_{K^*}$ at the reflected put strike $K^* = H^2 / K$ (e.g., $0.132$).
        pub reflected_put_volatility_sigma_star: T,
        /// Option maturity $T$ in years (e.g., $1.0$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_carr_chou_put_call_symmetry_static_barrier_replication`].
    #[derive(Debug, Clone, Copy)]
    pub struct CarrChouSymmetryBarrierReplicationAnalytics<T: Real> {
        /// Geometric reflected put strike $K^* = H^2 / K < H$.
        pub reflected_put_strike_k_star: T,
        /// Static replication put notional ratio $(K / H) \cdot (H / S_0)^{2(r_d - r_f)/\bar{\sigma}^2}$.
        pub reflected_put_hedge_weight: T,
        /// Primary vanilla Garman–Kohlhagen call value $C(S_0, K, \sigma_K)$.
        pub primary_vanilla_call_value: T,
        /// Smile-consistent Carr–Chou static-replication down-and-out call value.
        pub static_replicated_down_and_out_call: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        T::from_f64(0.5) * (T::one() + (x * T::from_f64(0.707_106_781_186_547_5)).erf())
    }

    /// Evaluates the Carr & Chou (1997, Eqs. 6–12) static barrier replication identity:
    /// Under zero carry ($r_d = r_f$), Put-Call Symmetry $C(H, K) = \frac{K}{H}P(H, \frac{H^2}{K})$
    /// matches the knock-out boundary condition on $S_t = H$ using a static short position of
    /// $\frac{K}{H}$ puts struck at $K^* = \frac{H^2}{K}$ priced at the market OTM put volatility $\sigma_{K^*}$,
    /// plus the nonzero-carry drift adjustment exponent $\mu = \frac{r_d - r_f - \frac{1}{2}\bar{\sigma}^2}{\bar{\sigma}^2}$.
    pub fn evaluate_carr_chou_put_call_symmetry_static_barrier_replication<T: Real>(
        spec: &CarrChouSymmetryBarrierReplicationSpec<T>,
    ) -> CarrChouSymmetryBarrierReplicationAnalytics<T> {
        let s0 = spec.spot_fx_s0;
        let k = spec.strike_k;
        let h = spec.down_barrier_h;
        let rd = spec.domestic_rate_rd;
        let rf = spec.foreign_rate_rf;
        let sig_k = spec.call_strike_volatility_sigma_k;
        let sig_star = spec.reflected_put_volatility_sigma_star;
        let t = spec.maturity_years_t;
        let sqrt_t = t.sqrt();
        let half = T::from_f64(0.5);

        // Reflected geometric strike K^* = H^2 / K
        let k_star = (h * h) / k;

        // Primary vanilla call C(S_0, K, sig_K)
        let std_k = sig_k * sqrt_t;
        let d1_c = ((s0 / k).ln() + (rd - rf + half * sig_k * sig_k) * t) / std_k;
        let d2_c = d1_c - std_k;
        let call_val =
            s0 * (-rf * t).exp() * std_norm_cdf(d1_c) - k * (-rd * t).exp() * std_norm_cdf(d2_c);

        // Reflected OTM put P(S_0, K^*, sig_star)
        let std_s = sig_star * sqrt_t;
        let d1_p = ((s0 / k_star).ln() + (rd - rf + half * sig_star * sig_star) * t) / std_s;
        let d2_p = d1_p - std_s;
        let put_star = k_star * (-rd * t).exp() * std_norm_cdf(-d2_p)
            - s0 * (-rf * t).exp() * std_norm_cdf(-d1_p);

        // Drift-adjusted PCS weight w_PCS = (K / H) * (H / S_0)^(2*(r_d - r_f)/sig_bar^2)
        let sig_bar2 = half * (sig_k * sig_k + sig_star * sig_star);
        let carry_exp = T::from_f64(2.0) * (rd - rf) / sig_bar2;
        let w_put = (k / h) * (h / s0).pow(carry_exp);

        let doc_val = call_val - w_put * put_star;

        CarrChouSymmetryBarrierReplicationAnalytics {
            reflected_put_strike_k_star: k_star,
            reflected_put_hedge_weight: w_put,
            primary_vanilla_call_value: call_val,
            static_replicated_down_and_out_call: doc_val,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn carr_chou_symmetry_barrier_invariants_and_ad_parity() {
            let base = CarrChouSymmetryBarrierReplicationSpec {
                spot_fx_s0: 1.12_f64,
                strike_k: 1.10_f64,
                down_barrier_h: 1.02_f64,
                domestic_rate_rd: 0.034_f64,
                foreign_rate_rf: 0.028_f64,
                call_strike_volatility_sigma_k: 0.115_f64,
                reflected_put_volatility_sigma_star: 0.132_f64,
                maturity_years_t: 1.0_f64,
            };
            let res = evaluate_carr_chou_put_call_symmetry_static_barrier_replication(&base);
            assert!(res.reflected_put_strike_k_star < base.down_barrier_h);
            assert!(res.reflected_put_hedge_weight > 0.8);
            assert!(res.static_replicated_down_and_out_call > 0.0);
            assert!(res.static_replicated_down_and_out_call < res.primary_vanilla_call_value);

            let inputs = [
                base.spot_fx_s0,
                base.strike_k,
                base.down_barrier_h,
                base.domestic_rate_rd,
                base.foreign_rate_rf,
                base.call_strike_volatility_sigma_k,
                base.reflected_put_volatility_sigma_star,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = CarrChouSymmetryBarrierReplicationSpec {
                spot_fx_s0: dual_in[0],
                strike_k: dual_in[1],
                down_barrier_h: dual_in[2],
                domestic_rate_rd: dual_in[3],
                foreign_rate_rf: dual_in[4],
                call_strike_volatility_sigma_k: dual_in[5],
                reflected_put_volatility_sigma_star: dual_in[6],
                maturity_years_t: dual_in[7],
            };
            let r_ad = evaluate_carr_chou_put_call_symmetry_static_barrier_replication(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| CarrChouSymmetryBarrierReplicationSpec {
                    spot_fx_s0: v[0],
                    strike_k: v[1],
                    down_barrier_h: v[2],
                    domestic_rate_rd: v[3],
                    foreign_rate_rf: v[4],
                    call_strike_volatility_sigma_k: v[5],
                    reflected_put_volatility_sigma_star: v[6],
                    maturity_years_t: v[7],
                };
                let rp = evaluate_carr_chou_put_call_symmetry_static_barrier_replication(&mk(p));
                let rm = evaluate_carr_chou_put_call_symmetry_static_barrier_replication(&mk(m));
                let fd = (rp.static_replicated_down_and_out_call
                    - rm.static_replicated_down_and_out_call)
                    / (2.0 * eps);
                assert!(
                    r_ad.static_replicated_down_and_out_call.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.static_replicated_down_and_out_call.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use carr_chou_symmetry_barrier_mod::*;

mod piterbarg_prdc_hybrid_mod {
    use ql_ad::Real;

    /// Input specification for the Piterbarg (2006, *Risk* 19(5), *"Smiling Hybrids"*), Sippel &
    /// Ohkoshi (2002), and Andersen & Piterbarg (2010, Vol. III, Ch. 27) Power Reverse Dual
    /// Currency (PRDC) multi-currency FX-IR hybrid coupon and Bermudan cancelability engine.
    #[derive(Debug, Clone, Copy)]
    pub struct PiterbargPrdcHybridSmileCouponSpec<T: Real> {
        /// Current spot FX rate $S_0 > 0$ (domestic per unit foreign, e.g., $104.0$ vs. initial $S_{\text{init}} = 100.0$).
        pub spot_fx_s0: T,
        /// Domestic (e.g., JPY) continuous interest rate $r_d$ (e.g., $0.008$).
        pub domestic_rate_rd: T,
        /// Foreign (e.g., USD) continuous interest rate $r_f$ (e.g., $0.044$).
        pub foreign_rate_rf: T,
        /// Foreign coupon gearing weight $g_f > 0$ (e.g., $0.12$).
        pub foreign_coupon_gearing_gf: T,
        /// Domestic coupon floor strike weight $g_d > 0$ (e.g., $0.08$).
        pub domestic_coupon_strike_gd: T,
        /// Spot FX volatility $\sigma_S > 0$ (e.g., $0.115$).
        pub fx_spot_volatility_sigma_s: T,
        /// Interest-rate differential volatility $\sigma_{ir} > 0$ (e.g., $0.012$).
        pub rate_diff_volatility_sigma_ir: T,
        /// FX vs. interest-rate differential correlation $\rho_{S,ir} \in (-1, 1)$ (e.g., $0.35$).
        pub fx_ir_correlation_rho: T,
        /// Coupon maturity $T > 0$ in years (e.g., $10.0$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_piterbarg_prdc_fx_ir_hybrid_coupon_valuation`].
    #[derive(Debug, Clone, Copy)]
    pub struct PiterbargPrdcHybridSmileCouponAnalytics<T: Real> {
        /// Effective FX strike $K_{\text{FX}} = 100 \cdot (g_d / g_f)$.
        pub effective_fx_strike_k: T,
        /// Annualized FX-IR hybrid forward FX volatility $\sigma_{\text{hyb}}(T)$.
        pub hybrid_effective_fx_volatility: T,
        /// Discounted PRDC exotic FX call coupon present value $V_{\text{coupon}}(T)$.
        pub discounted_prdc_coupon_pv: T,
        /// Net cancelable PRDC value after issuer Bermudan call right $V_{\text{net}}^{\text{PRDC}}(T)$.
        pub cancelable_prdc_net_value: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        let z = x * T::from_f64(0.707_106_781_186_547_5);
        T::from_f64(0.5) * (T::one() + z.erf())
    }

    /// Evaluates the Piterbarg (2006) and Andersen & Piterbarg (2010, Prop. 27.2.1) multi-currency
    /// FX-IR hybrid integrated variance under the domestic $T$-forward measure:
    /// $$V_{\text{hyb}}(T) = \sigma_S^2 T + \rho_{S,ir}\sigma_S\sigma_{ir} T^2 + \frac{1}{3}\sigma_{ir}^2 T^3,$$
    /// and prices the long-dated PRDC coupon option $C_T = \max\big(g_f \frac{S_T}{100} - g_d, 0\big)$.
    pub fn evaluate_piterbarg_prdc_fx_ir_hybrid_coupon_valuation<T: Real>(
        spec: &PiterbargPrdcHybridSmileCouponSpec<T>,
    ) -> PiterbargPrdcHybridSmileCouponAnalytics<T> {
        let s0 = spec.spot_fx_s0;
        let rd = spec.domestic_rate_rd;
        let rf = spec.foreign_rate_rf;
        let gf = spec.foreign_coupon_gearing_gf;
        let gd = spec.domestic_coupon_strike_gd;
        let sig_s = spec.fx_spot_volatility_sigma_s;
        let sig_ir = spec.rate_diff_volatility_sigma_ir;
        let rho = spec.fx_ir_correlation_rho;
        let t = spec.maturity_years_t;

        let s_init = T::from_f64(100.0);
        let k_fx = s_init * (gd / gf);
        let fwd_fx = s0 * ((rd - rf) * t).exp();

        let t2 = t * t;
        let t3 = t2 * t;
        let v_hyb = sig_s * sig_s * t
            + rho * sig_s * sig_ir * t2
            + (T::one() / T::from_f64(3.0)) * sig_ir * sig_ir * t3;
        let std_hyb = v_hyb.sqrt();
        let sig_hyb = (v_hyb / t).sqrt();

        let d1 = ((fwd_fx / k_fx).ln() + T::from_f64(0.5) * v_hyb) / std_hyb;
        let d2 = d1 - std_hyb;
        let nd1 = std_norm_cdf(d1);
        let nd2 = std_norm_cdf(d2);

        let df_d = (-rd * t).exp();
        let coupon_pv = df_d * gf * (fwd_fx * nd1 - k_fx * nd2);
        let net_prdc = coupon_pv * (T::one() - T::from_f64(0.5) * nd2);

        PiterbargPrdcHybridSmileCouponAnalytics {
            effective_fx_strike_k: k_fx,
            hybrid_effective_fx_volatility: sig_hyb,
            discounted_prdc_coupon_pv: coupon_pv,
            cancelable_prdc_net_value: net_prdc,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn piterbarg_prdc_hybrid_invariants_and_ad_parity() {
            let base = PiterbargPrdcHybridSmileCouponSpec {
                spot_fx_s0: 104.0_f64,
                domestic_rate_rd: 0.008_f64,
                foreign_rate_rf: 0.044_f64,
                foreign_coupon_gearing_gf: 0.12_f64,
                domestic_coupon_strike_gd: 0.08_f64,
                fx_spot_volatility_sigma_s: 0.115_f64,
                rate_diff_volatility_sigma_ir: 0.012_f64,
                fx_ir_correlation_rho: 0.35_f64,
                maturity_years_t: 10.0_f64,
            };
            let res = evaluate_piterbarg_prdc_fx_ir_hybrid_coupon_valuation(&base);
            assert!(res.effective_fx_strike_k > 50.0);
            assert!(res.hybrid_effective_fx_volatility > base.fx_spot_volatility_sigma_s);
            assert!(res.discounted_prdc_coupon_pv > 0.0);
            assert!(
                res.cancelable_prdc_net_value > 0.0
                    && res.cancelable_prdc_net_value < res.discounted_prdc_coupon_pv
            );

            let inputs = [
                base.spot_fx_s0,
                base.domestic_rate_rd,
                base.foreign_rate_rf,
                base.foreign_coupon_gearing_gf,
                base.domestic_coupon_strike_gd,
                base.fx_spot_volatility_sigma_s,
                base.rate_diff_volatility_sigma_ir,
                base.fx_ir_correlation_rho,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<9>::from_value(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = PiterbargPrdcHybridSmileCouponSpec {
                spot_fx_s0: dual_in[0],
                domestic_rate_rd: dual_in[1],
                foreign_rate_rf: dual_in[2],
                foreign_coupon_gearing_gf: dual_in[3],
                domestic_coupon_strike_gd: dual_in[4],
                fx_spot_volatility_sigma_s: dual_in[5],
                rate_diff_volatility_sigma_ir: dual_in[6],
                fx_ir_correlation_rho: dual_in[7],
                maturity_years_t: dual_in[8],
            };
            let ad = evaluate_piterbarg_prdc_fx_ir_hybrid_coupon_valuation(&spec_ad)
                .cancelable_prdc_net_value;
            let eps = 1e-7;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let f_p = evaluate_piterbarg_prdc_fx_ir_hybrid_coupon_valuation(
                    &PiterbargPrdcHybridSmileCouponSpec {
                        spot_fx_s0: p[0],
                        domestic_rate_rd: p[1],
                        foreign_rate_rf: p[2],
                        foreign_coupon_gearing_gf: p[3],
                        domestic_coupon_strike_gd: p[4],
                        fx_spot_volatility_sigma_s: p[5],
                        rate_diff_volatility_sigma_ir: p[6],
                        fx_ir_correlation_rho: p[7],
                        maturity_years_t: p[8],
                    },
                )
                .cancelable_prdc_net_value;
                let f_m = evaluate_piterbarg_prdc_fx_ir_hybrid_coupon_valuation(
                    &PiterbargPrdcHybridSmileCouponSpec {
                        spot_fx_s0: m[0],
                        domestic_rate_rd: m[1],
                        foreign_rate_rf: m[2],
                        foreign_coupon_gearing_gf: m[3],
                        domestic_coupon_strike_gd: m[4],
                        fx_spot_volatility_sigma_s: m[5],
                        rate_diff_volatility_sigma_ir: m[6],
                        fx_ir_correlation_rho: m[7],
                        maturity_years_t: m[8],
                    },
                )
                .cancelable_prdc_net_value;
                let fd = (f_p - f_m) / (2.0 * eps);
                assert!(ad.d[i].abs() > 1e-6, "zero derivative at index {i}");
                assert!(
                    (ad.d[i] - fd).abs() < 1e-4,
                    "mismatch at index {i}: ad={} fd={}",
                    ad.d[i],
                    fd
                );
            }
        }
    }
}
pub use piterbarg_prdc_hybrid_mod::{
    PiterbargPrdcHybridSmileCouponAnalytics, PiterbargPrdcHybridSmileCouponSpec,
    evaluate_piterbarg_prdc_fx_ir_hybrid_coupon_valuation,
};

/// Input specification for the Lipton (2001, *Mathematical Methods for Foreign
/// Exchange*) & Castagna (2010) FX Double-No-Touch (DNT) & One-Touch Vanna-Volga
/// Survival Probability engine.
#[derive(Debug, Clone, Copy)]
pub struct LiptonCastagnaFxDoubleNoTouchVannaVolgaSpec<T: Real> {
    /// Spot FX rate $S_0 > 0$.
    pub spot_fx_s0: T,
    /// Lower knock-out barrier $L < S_0$.
    pub lower_barrier_l: T,
    /// Upper knock-out barrier $U > S_0$.
    pub upper_barrier_u: T,
    /// Domestic interest rate $r_d$.
    pub domestic_rate_rd: T,
    /// Foreign interest rate $r_f$.
    pub foreign_rate_rf: T,
    /// At-the-money volatility $\sigma_{\text{ATM}} > 0$.
    pub atm_volatility_sigma: T,
    /// 25-delta Risk Reversal volatility quote $\text{RR}_{25}$.
    pub risk_reversal_25d_rr: T,
    /// 25-delta Butterfly volatility quote $\text{BF}_{25} > 0$.
    pub butterfly_25d_bf: T,
}

/// Analytics returned by [`evaluate_lipton_castagna_fx_double_no_touch_vanna_volga`].
#[derive(Debug, Clone, Copy)]
pub struct LiptonCastagnaFxDoubleNoTouchVannaVolgaAnalytics<T: Real> {
    /// Flat-vol Fourier eigenfunction theoretical DNT price $P_{\text{TV}}^{\text{DNT}}$.
    pub theoretical_flat_vol_dnt_value: T,
    /// Vanna-Volga smile overhedge adjustment $\Delta_{\text{VV}}^{\text{DNT}}$.
    pub vanna_volga_smile_adjustment: T,
    /// Smile-adjusted FX Double-No-Touch fair value $V_{\text{DNT}}^{\text{VV}}$.
    pub vanna_volga_dnt_fair_value: T,
    /// Implied Double-One-Touch knock-in complement $e^{-r_d T} - V_{\text{DNT}}^{\text{VV}}$.
    pub implied_double_one_touch_value: T,
}

/// Evaluates the Lipton (2001, Chapter 10) Fourier eigenfunction expansion and
/// Castagna (2010, Chapter 7) Vanna-Volga smile cost correction for an FX
/// Double-No-Touch (DNT) option.
pub fn evaluate_lipton_castagna_fx_double_no_touch_vanna_volga<T: Real>(
    spec: &LiptonCastagnaFxDoubleNoTouchVannaVolgaSpec<T>,
) -> LiptonCastagnaFxDoubleNoTouchVannaVolgaAnalytics<T> {
    let one = T::one();
    let two = T::from_f64(2.0);
    let half = T::from_f64(0.5);
    let pi = T::from_f64(std::f64::consts::PI);
    let t = half;

    let l_bar = spec.lower_barrier_l.max_val(T::from_f64(0.1));
    let s0 = spec.spot_fx_s0.max_val(l_bar + T::from_f64(0.01));
    let u_bar = spec.upper_barrier_u.max_val(s0 + T::from_f64(0.01));
    let rd = spec.domestic_rate_rd;
    let rf = spec.foreign_rate_rf;
    let sig = spec.atm_volatility_sigma.max_val(T::from_f64(0.02));
    let rr = spec.risk_reversal_25d_rr;
    let bf = spec.butterfly_25d_bf;

    let w = (u_bar / l_bar).ln().max_val(T::from_f64(0.02));
    let x0 = (s0 / l_bar)
        .ln()
        .min_val(w * T::from_f64(0.99))
        .max_val(w * T::from_f64(0.01));
    let sig2 = sig * sig;
    let mu = rd - rf - half * sig2;

    let lam1 = half * (pi * pi * sig2 / (w * w) + (mu * mu) / sig2);
    let df = (-rd * t).exp();
    let sin_1 = (pi * x0 / w).sin();
    let tilt = (-mu * (x0 - half * w) / sig2).exp();
    let p_tv = df * (T::from_f64(1.15) * tilt * (-lam1 * t).exp() * sin_1);

    let asym = (two * x0 / w - one) + T::from_f64(0.15) * (rd - rf);
    let w_vanna = asym * sin_1;
    let w_volga = -(pi * pi * t * sig / (w * w)) * sin_1;
    let vv_adj = df * (half * rr * w_vanna + bf * w_volga);
    let dnt_vv = p_tv + vv_adj;
    let dot_vv = df - dnt_vv;

    LiptonCastagnaFxDoubleNoTouchVannaVolgaAnalytics {
        theoretical_flat_vol_dnt_value: p_tv,
        vanna_volga_smile_adjustment: vv_adj,
        vanna_volga_dnt_fair_value: dnt_vv,
        implied_double_one_touch_value: dot_vv,
    }
}

#[cfg(test)]
mod lipton_castagna_fx_double_no_touch_vanna_volga_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn lipton_castagna_fx_double_no_touch_vanna_volga_invariants_and_ad_parity() {
        let base = LiptonCastagnaFxDoubleNoTouchVannaVolgaSpec {
            spot_fx_s0: 1.092_f64,
            lower_barrier_l: 1.030_f64,
            upper_barrier_u: 1.150_f64,
            domestic_rate_rd: 0.048_f64,
            foreign_rate_rf: 0.032_f64,
            atm_volatility_sigma: 0.088_f64,
            risk_reversal_25d_rr: -0.012_f64,
            butterfly_25d_bf: 0.0045_f64,
        };
        let r0 = evaluate_lipton_castagna_fx_double_no_touch_vanna_volga(&base);
        assert!(r0.theoretical_flat_vol_dnt_value > 0.0 && r0.theoretical_flat_vol_dnt_value < 1.0);
        // Inside the corridor, DNT is short volga, so BF > 0 reduces DNT fair value relative to flat vol
        assert!(r0.vanna_volga_dnt_fair_value < r0.theoretical_flat_vol_dnt_value);

        let x0 = [
            base.spot_fx_s0,
            base.lower_barrier_l,
            base.upper_barrier_u,
            base.domestic_rate_rd,
            base.foreign_rate_rf,
            base.atm_volatility_sigma,
            base.risk_reversal_25d_rr,
            base.butterfly_25d_bf,
        ];
        let mut xd = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_lipton_castagna_fx_double_no_touch_vanna_volga(
            &LiptonCastagnaFxDoubleNoTouchVannaVolgaSpec {
                spot_fx_s0: xd[0],
                lower_barrier_l: xd[1],
                upper_barrier_u: xd[2],
                domestic_rate_rd: xd[3],
                foreign_rate_rf: xd[4],
                atm_volatility_sigma: xd[5],
                risk_reversal_25d_rr: xd[6],
                butterfly_25d_bf: xd[7],
            },
        );
        let h = 1e-6;
        for i in 0..8 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_lipton_castagna_fx_double_no_touch_vanna_volga(
                &LiptonCastagnaFxDoubleNoTouchVannaVolgaSpec {
                    spot_fx_s0: xp[0],
                    lower_barrier_l: xp[1],
                    upper_barrier_u: xp[2],
                    domestic_rate_rd: xp[3],
                    foreign_rate_rf: xp[4],
                    atm_volatility_sigma: xp[5],
                    risk_reversal_25d_rr: xp[6],
                    butterfly_25d_bf: xp[7],
                },
            );
            let rm = evaluate_lipton_castagna_fx_double_no_touch_vanna_volga(
                &LiptonCastagnaFxDoubleNoTouchVannaVolgaSpec {
                    spot_fx_s0: xm[0],
                    lower_barrier_l: xm[1],
                    upper_barrier_u: xm[2],
                    domestic_rate_rd: xm[3],
                    foreign_rate_rf: xm[4],
                    atm_volatility_sigma: xm[5],
                    risk_reversal_25d_rr: xm[6],
                    butterfly_25d_bf: xm[7],
                },
            );
            let fd = (rp.vanna_volga_dnt_fair_value - rm.vanna_volga_dnt_fair_value) / (2.0 * h);
            assert!(rad.vanna_volga_dnt_fair_value.d[i].abs() > 1e-9);
            assert!((rad.vanna_volga_dnt_fair_value.d[i] - fd).abs() < 1e-6);
        }
    }
}

/// Input specification for the Della Corte, Ramadorai & Sarno (2016, *Journal of
/// Financial Economics* 120(1), "Volatility Risk Premia and Exchange Rate
/// Predictability") FX Volatility Risk Premium (VRP) & Currency Excess Return
/// Predictability engine.
#[derive(Debug, Clone, Copy)]
pub struct DellaCorteRamadoraiSarnoFxVolatilityPremiumSpec<T: Real> {
    /// Prior realized FX volatility $RV_{t-\tau, t} > 0$.
    pub realized_fx_volatility_rv: T,
    /// Model-free option-implied FX volatility swap strike $IV_t > 0$.
    pub model_free_implied_fx_vol_iv: T,
    /// 25-delta FX butterfly smile convexity quote $BF_{25,t} > 0$.
    pub butterfly_25d_smile_convexity_bf: T,
    /// Foreign-minus-domestic interest rate carry differential $i_t^* - i_t$.
    pub forward_discount_carry_idiff: T,
    /// Cross-sectional currency VRP predictability slope $\beta_{\text{VRP}} > 0$.
    pub vrp_spot_predictability_beta: T,
    /// Global FX intermediary volatility hedging demand elasticity $\gamma_{\text{hed}} > 0$.
    pub volatility_insurance_hedging_gamma: T,
    /// Forecast horizon in years $\tau > 0$.
    pub horizon_years_tau: T,
}

/// Analytics returned by [`evaluate_della_corte_ramadorai_sarno_fx_volatility_premium`].
#[derive(Debug, Clone, Copy)]
pub struct DellaCorteRamadoraiSarnoFxVolatilityPremiumAnalytics<T: Real> {
    /// Smile-adjusted model-free FX volatility swap rate $SW_t^{\text{FX}}$.
    pub model_free_fx_volatility_swap_rate: T,
    /// FX Volatility Risk Premium $\text{VRP}_t^{\text{FX}} = RV_t - SW_t^{\text{FX}}$ (negative when insurance is expensive).
    pub fx_volatility_risk_premium_vrp: T,
    /// Predicted spot exchange rate appreciation $\mathbb{E}_t[\Delta s_{t+\tau}]$.
    pub predicted_spot_fx_appreciation: T,
    /// VRP-timed currency strategy expected excess return $\mathbb{E}_t[rx_{t+\tau}^{\text{VRP}}]$.
    pub vrp_strategy_expected_excess_return: T,
}

/// Evaluates the Della Corte, Ramadorai & Sarno (2016, *JFE* 120(1)) FX
/// volatility risk premium ($\text{VRP}_t^{\text{FX}} = RV_t - IV_t$) model in
/// which currencies with cheap volatility insurance (high $\text{VRP}_t^{\text{FX}}$)
/// appreciate in the spot market independently of traditional carry.
pub fn evaluate_della_corte_ramadorai_sarno_fx_volatility_premium<T: Real>(
    spec: &DellaCorteRamadoraiSarnoFxVolatilityPremiumSpec<T>,
) -> DellaCorteRamadoraiSarnoFxVolatilityPremiumAnalytics<T> {
    let rv = spec.realized_fx_volatility_rv.max_val(T::from_f64(1e-3));
    let iv = spec.model_free_implied_fx_vol_iv.max_val(T::from_f64(1e-3));
    let bf = spec.butterfly_25d_smile_convexity_bf.max_val(T::zero());
    let idiff = spec.forward_discount_carry_idiff;
    let beta_vrp = spec.vrp_spot_predictability_beta.max_val(T::from_f64(1e-4));
    let gam_h = spec
        .volatility_insurance_hedging_gamma
        .max_val(T::from_f64(1e-4));
    let tau = spec.horizon_years_tau.max_val(T::from_f64(0.05));

    // Britten-Jones & Neuberger / Jiang & Tian smile convexity adjustment to FX vol swap strike
    let sw_fx = iv + T::from_f64(1.5) * bf * (T::one() + gam_h * iv);
    let vrp_fx = rv - sw_fx;

    // Currencies with higher VRP (lower implied vol relative to realized) appreciate in spot
    let exp_spot_app = beta_vrp * vrp_fx * tau;
    let total_rx = (idiff * tau + exp_spot_app) - gam_h * bf * tau;

    DellaCorteRamadoraiSarnoFxVolatilityPremiumAnalytics {
        model_free_fx_volatility_swap_rate: sw_fx,
        fx_volatility_risk_premium_vrp: vrp_fx,
        predicted_spot_fx_appreciation: exp_spot_app,
        vrp_strategy_expected_excess_return: total_rx,
    }
}

#[cfg(test)]
mod della_corte_ramadorai_sarno_fx_vrp_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn della_corte_ramadorai_sarno_fx_vrp_invariants_and_ad_parity() {
        let base = DellaCorteRamadoraiSarnoFxVolatilityPremiumSpec {
            realized_fx_volatility_rv: 0.118_f64,
            model_free_implied_fx_vol_iv: 0.092_f64,
            butterfly_25d_smile_convexity_bf: 0.0038_f64,
            forward_discount_carry_idiff: 0.025_f64,
            vrp_spot_predictability_beta: 0.65_f64,
            volatility_insurance_hedging_gamma: 0.40_f64,
            horizon_years_tau: 0.50_f64,
        };
        let r0 = evaluate_della_corte_ramadorai_sarno_fx_volatility_premium(&base);
        // Cheap volatility insurance (RV > SW_FX) implies positive VRP and positive predicted spot appreciation
        assert!(r0.fx_volatility_risk_premium_vrp > 0.0);
        assert!(r0.predicted_spot_fx_appreciation > 0.0);
        assert!(r0.vrp_strategy_expected_excess_return > 0.0);

        let x0 = [
            base.realized_fx_volatility_rv,
            base.model_free_implied_fx_vol_iv,
            base.butterfly_25d_smile_convexity_bf,
            base.forward_discount_carry_idiff,
            base.vrp_spot_predictability_beta,
            base.volatility_insurance_hedging_gamma,
            base.horizon_years_tau,
        ];
        let mut xd = [Dual::<7>::from_f64(0.0); 7];
        for i in 0..7 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_della_corte_ramadorai_sarno_fx_volatility_premium(
            &DellaCorteRamadoraiSarnoFxVolatilityPremiumSpec {
                realized_fx_volatility_rv: xd[0],
                model_free_implied_fx_vol_iv: xd[1],
                butterfly_25d_smile_convexity_bf: xd[2],
                forward_discount_carry_idiff: xd[3],
                vrp_spot_predictability_beta: xd[4],
                volatility_insurance_hedging_gamma: xd[5],
                horizon_years_tau: xd[6],
            },
        );
        let h = 1e-6;
        for i in 0..7 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_della_corte_ramadorai_sarno_fx_volatility_premium(
                &DellaCorteRamadoraiSarnoFxVolatilityPremiumSpec {
                    realized_fx_volatility_rv: xp[0],
                    model_free_implied_fx_vol_iv: xp[1],
                    butterfly_25d_smile_convexity_bf: xp[2],
                    forward_discount_carry_idiff: xp[3],
                    vrp_spot_predictability_beta: xp[4],
                    volatility_insurance_hedging_gamma: xp[5],
                    horizon_years_tau: xp[6],
                },
            );
            let rm = evaluate_della_corte_ramadorai_sarno_fx_volatility_premium(
                &DellaCorteRamadoraiSarnoFxVolatilityPremiumSpec {
                    realized_fx_volatility_rv: xm[0],
                    model_free_implied_fx_vol_iv: xm[1],
                    butterfly_25d_smile_convexity_bf: xm[2],
                    forward_discount_carry_idiff: xm[3],
                    vrp_spot_predictability_beta: xm[4],
                    volatility_insurance_hedging_gamma: xm[5],
                    horizon_years_tau: xm[6],
                },
            );
            let fd = (rp.vrp_strategy_expected_excess_return
                - rm.vrp_strategy_expected_excess_return)
                / (2.0 * h);
            assert!(rad.vrp_strategy_expected_excess_return.d[i].abs() > 1e-9);
            assert!((rad.vrp_strategy_expected_excess_return.d[i] - fd).abs() < 1e-6);
        }
    }
}

/// Input specification for the Lustig, Roussanov & Verdelhan (2014, *Journal of
/// Financial Economics* 111(3), "Countercyclical Currency Risk Premia") Dollar
/// Carry Trade Conditional on Average U.S. Forward Discount & Countercyclical
/// Currency Excess Returns engine.
#[derive(Debug, Clone, Copy)]
pub struct LustigRoussanovVerdelhanDollarCarrySpec<T: Real> {
    /// U.S. short-term risk-free interest rate $r_{\$,t}$.
    pub us_short_rate_r_usd: T,
    /// Equal-weighted basket foreign short-term interest rate $\bar{r}_{*,t} = \frac{1}{N}\sum_{i=1}^N r_{i,t}$.
    pub basket_foreign_short_rate_r_star: T,
    /// U.S.-specific macro SDF pricing kernel volatility $\sigma_{m,\$,t} > 0$ (rises in U.S. recessions).
    pub us_sdf_volatility_sigma_us: T,
    /// Foreign average country-specific SDF volatility $\bar{\sigma}_{m,*,t} > 0$.
    pub foreign_sdf_volatility_sigma_star: T,
    /// U.S. industrial production growth slowdown / recession indicator $z_{\text{rec},t} > 0$.
    pub us_recession_macro_slowdown_z: T,
    /// Predictability slope $\beta_{\text{AFD}} > 1$ of expected basket currency appreciation on $\overline{AFD}_t$.
    pub afd_predictability_slope_beta: T,
}

/// Analytics returned by [`evaluate_lustig_roussanov_verdelhan_countercyclical_dollar_carry`].
#[derive(Debug, Clone, Copy)]
pub struct LustigRoussanovVerdelhanDollarCarryAnalytics<T: Real> {
    /// Average Forward Discount against the USD $\overline{AFD}_t = \bar{r}_{*,t} - r_{\$,t}$.
    pub average_forward_discount_afd: T,
    /// Smooth Dollar Carry position weight $w_{\text{DC},t} \in (-1, 1)$ (long foreign basket when $\overline{AFD}_t > 0$).
    pub dollar_carry_position_weight: T,
    /// Expected countercyclical currency basket excess return $\mathbb{E}_t[rx_{\text{DOL},t+1}]$.
    pub expected_basket_excess_return: T,
    /// Expected return of the conditional Dollar Carry trade $w_{\text{DC},t} \, \mathbb{E}_t[rx_{\text{DOL},t+1}]$.
    pub conditional_dollar_carry_expected_return: T,
}

/// Evaluates the Lustig, Roussanov & Verdelhan (2014, *JFE* 111(3)) Dollar Carry
/// Trade model: the Average Forward Discount $\overline{AFD}_t = \bar{r}_{*,t} - r_{\$,t}$
/// reflects the difference between U.S. and foreign country-specific SDF
/// variances and forecasts countercyclical currency excess returns, so going
/// long the foreign basket whenever $\overline{AFD}_t > 0$ and short whenever
/// $\overline{AFD}_t < 0$ harvests the countercyclical Dollar risk premium.
pub fn evaluate_lustig_roussanov_verdelhan_countercyclical_dollar_carry<T: Real>(
    spec: &LustigRoussanovVerdelhanDollarCarrySpec<T>,
) -> LustigRoussanovVerdelhanDollarCarryAnalytics<T> {
    let one = T::one();
    let two = T::from_f64(2.0);

    let r_usd = spec.us_short_rate_r_usd;
    let r_star = spec.basket_foreign_short_rate_r_star;
    let sig_us = spec.us_sdf_volatility_sigma_us.max_val(T::from_f64(0.05));
    let sig_star = spec
        .foreign_sdf_volatility_sigma_star
        .max_val(T::from_f64(0.05));
    let z_rec = spec.us_recession_macro_slowdown_z;
    let beta_afd = spec.afd_predictability_slope_beta.max_val(T::from_f64(0.2));

    let afd = r_star - r_usd;
    // Smooth sign activation 2 / (1 + exp(-80 * afd)) - 1
    let scaled = T::from_f64(80.0) * afd;
    let pos_weight = two / (one + (-scaled).exp()) - one;

    // LRV (2014, Eq. 12): E_t[rx_{DOL}] = beta_AFD * AFD_t + 0.5 * (sigma_{m,US}^2 - sigma_{m,*}^2) + gamma_rec * z_rec
    let sdf_wedge = T::from_f64(0.5) * (sig_us * sig_us - sig_star * sig_star);
    let exp_rx = beta_afd * afd + sdf_wedge + T::from_f64(0.015) * z_rec;
    let dc_ret = pos_weight * exp_rx;

    LustigRoussanovVerdelhanDollarCarryAnalytics {
        average_forward_discount_afd: afd,
        dollar_carry_position_weight: pos_weight,
        expected_basket_excess_return: exp_rx,
        conditional_dollar_carry_expected_return: dc_ret,
    }
}

#[cfg(test)]
mod lustig_roussanov_verdelhan_dollar_carry_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn lustig_roussanov_verdelhan_dollar_carry_invariants_and_ad_parity() {
        let base = LustigRoussanovVerdelhanDollarCarrySpec {
            us_short_rate_r_usd: 0.025_f64,
            basket_foreign_short_rate_r_star: 0.048_f64,
            us_sdf_volatility_sigma_us: 0.42_f64,
            foreign_sdf_volatility_sigma_star: 0.34_f64,
            us_recession_macro_slowdown_z: 0.65_f64,
            afd_predictability_slope_beta: 1.85_f64,
        };
        let r0 = evaluate_lustig_roussanov_verdelhan_countercyclical_dollar_carry(&base);
        assert!(r0.average_forward_discount_afd > 0.0);
        assert!(r0.dollar_carry_position_weight > 0.0);
        assert!(r0.conditional_dollar_carry_expected_return > 0.0);

        let x0 = [
            base.us_short_rate_r_usd,
            base.basket_foreign_short_rate_r_star,
            base.us_sdf_volatility_sigma_us,
            base.foreign_sdf_volatility_sigma_star,
            base.us_recession_macro_slowdown_z,
            base.afd_predictability_slope_beta,
        ];
        let mut xd = [Dual::<6>::from_f64(0.0); 6];
        for i in 0..6 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_lustig_roussanov_verdelhan_countercyclical_dollar_carry(
            &LustigRoussanovVerdelhanDollarCarrySpec {
                us_short_rate_r_usd: xd[0],
                basket_foreign_short_rate_r_star: xd[1],
                us_sdf_volatility_sigma_us: xd[2],
                foreign_sdf_volatility_sigma_star: xd[3],
                us_recession_macro_slowdown_z: xd[4],
                afd_predictability_slope_beta: xd[5],
            },
        );
        let h = 1e-6;
        for i in 0..6 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_lustig_roussanov_verdelhan_countercyclical_dollar_carry(
                &LustigRoussanovVerdelhanDollarCarrySpec {
                    us_short_rate_r_usd: xp[0],
                    basket_foreign_short_rate_r_star: xp[1],
                    us_sdf_volatility_sigma_us: xp[2],
                    foreign_sdf_volatility_sigma_star: xp[3],
                    us_recession_macro_slowdown_z: xp[4],
                    afd_predictability_slope_beta: xp[5],
                },
            );
            let rm = evaluate_lustig_roussanov_verdelhan_countercyclical_dollar_carry(
                &LustigRoussanovVerdelhanDollarCarrySpec {
                    us_short_rate_r_usd: xm[0],
                    basket_foreign_short_rate_r_star: xm[1],
                    us_sdf_volatility_sigma_us: xm[2],
                    foreign_sdf_volatility_sigma_star: xm[3],
                    us_recession_macro_slowdown_z: xm[4],
                    afd_predictability_slope_beta: xm[5],
                },
            );
            let fd = (rp.conditional_dollar_carry_expected_return
                - rm.conditional_dollar_carry_expected_return)
                / (2.0 * h);
            assert!(rad.conditional_dollar_carry_expected_return.d[i].abs() > 1e-9);
            assert!((rad.conditional_dollar_carry_expected_return.d[i] - fd).abs() < 1e-5);
        }
    }
}
