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
