mod evaluate_ljungqvist_new_keynesian_three_equation_ree_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct LjungqvistNewKeynesianReeSpec<T: Real = f64> {
        pub risk_aversion_sigma: T,
        pub discount_factor_beta: T,
        pub calvo_phillips_kappa: T,
        pub taylor_inflation_phi_pi: T,
        pub taylor_output_phi_y: T,
        pub natural_rate_persistence_rho: T,
        pub natural_rate_shock_rn: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct LjungqvistNewKeynesianReeReport<T: Real = f64> {
        pub blanchard_kahn_determinacy_margin: T,
        pub output_gap_loading_psi_y: T,
        pub inflation_loading_psi_pi: T,
        pub equilibrium_output_gap: T,
        pub equilibrium_inflation: T,
        pub equilibrium_policy_rate: T,
        pub ex_ante_real_interest_rate: T,
        pub ree_is_nkpc_conservation_error: T,
    }

    pub fn evaluate_ljungqvist_new_keynesian_three_equation_ree<T: Real>(
        spec: &LjungqvistNewKeynesianReeSpec<T>,
    ) -> LjungqvistNewKeynesianReeReport<T> {
        let sig = spec.risk_aversion_sigma.max_val(T::from_f64(1e-6));
        let beta = spec.discount_factor_beta;
        let kappa = spec.calvo_phillips_kappa;
        let phi_pi = spec.taylor_inflation_phi_pi;
        let phi_y = spec.taylor_output_phi_y;
        let rho = spec.natural_rate_persistence_rho;
        let rn = spec.natural_rate_shock_rn;

        let bk_margin = kappa * (phi_pi - T::one()) + (T::one() - beta) * phi_y;
        let one_m_brho = T::one() - beta * rho;
        let lam_nk = ((sig * (T::one() - rho) + phi_y) * one_m_brho + kappa * (phi_pi - rho))
            .max_val(T::from_f64(1e-8));
        let psi_y = one_m_brho / lam_nk;
        let psi_pi = kappa / lam_nk;

        let y_t = psi_y * rn;
        let pi_t = psi_pi * rn;
        let i_t = phi_pi * pi_t + phi_y * y_t;
        let ey_tp1 = rho * y_t;
        let epi_tp1 = rho * pi_t;
        let real_r = i_t - epi_tp1;

        let res_is = y_t - (ey_tp1 - (real_r - rn) / sig);
        let res_nkpc = pi_t - (beta * epi_tp1 + kappa * y_t);
        let err = (res_is * res_is + res_nkpc * res_nkpc + T::from_f64(1e-30)).sqrt();

        LjungqvistNewKeynesianReeReport {
            blanchard_kahn_determinacy_margin: bk_margin,
            output_gap_loading_psi_y: psi_y,
            inflation_loading_psi_pi: psi_pi,
            equilibrium_output_gap: y_t,
            equilibrium_inflation: pi_t,
            equilibrium_policy_rate: i_t,
            ex_ante_real_interest_rate: real_r,
            ree_is_nkpc_conservation_error: err,
        }
    }

    #[cfg(test)]
    mod domain_suite_49_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn ljungqvist_new_keynesian_ree_invariants_and_ad_parity() {
            let spec = LjungqvistNewKeynesianReeSpec {
                risk_aversion_sigma: 1.5,
                discount_factor_beta: 0.99,
                calvo_phillips_kappa: 0.15,
                taylor_inflation_phi_pi: 1.50,
                taylor_output_phi_y: 0.125,
                natural_rate_persistence_rho: 0.80,
                natural_rate_shock_rn: 0.015,
            };
            let rep = evaluate_ljungqvist_new_keynesian_three_equation_ree(&spec);
            assert!(rep.blanchard_kahn_determinacy_margin > 0.0);
            assert!(rep.ree_is_nkpc_conservation_error < 1e-12);

            let h = 1e-6;
            let r_ad = evaluate_ljungqvist_new_keynesian_three_equation_ree(
                &LjungqvistNewKeynesianReeSpec {
                    risk_aversion_sigma: Dual::<1>::from_value(spec.risk_aversion_sigma),
                    discount_factor_beta: Dual::<1>::from_value(spec.discount_factor_beta),
                    calvo_phillips_kappa: Dual::<1>::from_value(spec.calvo_phillips_kappa),
                    taylor_inflation_phi_pi: Dual::<1>::variable(spec.taylor_inflation_phi_pi, 0),
                    taylor_output_phi_y: Dual::<1>::from_value(spec.taylor_output_phi_y),
                    natural_rate_persistence_rho: Dual::<1>::from_value(
                        spec.natural_rate_persistence_rho,
                    ),
                    natural_rate_shock_rn: Dual::<1>::from_value(spec.natural_rate_shock_rn),
                },
            );
            let mut sp = spec;
            sp.taylor_inflation_phi_pi += h;
            let mut sm = spec;
            sm.taylor_inflation_phi_pi -= h;
            let fd = (evaluate_ljungqvist_new_keynesian_three_equation_ree(&sp)
                .equilibrium_policy_rate
                - evaluate_ljungqvist_new_keynesian_three_equation_ree(&sm)
                    .equilibrium_policy_rate)
                / (2.0 * h);
            assert!((r_ad.equilibrium_policy_rate.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_ljungqvist_new_keynesian_three_equation_ree_module::*;

mod evaluate_broadie_glasserman_kou_discrete_barrier_lookback_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Clone, Copy, Debug)]
    pub struct BroadieGlassermanKouDiscreteBarrierSpec<T: Real = f64> {
        pub spot_price: T,
        pub strike_price: T,
        pub down_barrier_h: T,
        pub risk_free_rate: T,
        pub dividend_yield: T,
        pub volatility: T,
        pub maturity_years: T,
        pub monitoring_steps_m: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct BroadieGlassermanKouDiscreteBarrierReport<T: Real = f64> {
        pub shifted_effective_barrier: T,
        pub continuous_down_and_out_call: T,
        pub discrete_bgk_down_and_out_call: T,
        pub continuity_correction_price_uplift: T,
        pub continuous_lookback_call: T,
        pub discrete_bgk_lookback_call: T,
    }

    fn bs_call_generic<T: Real>(s: T, k: T, r: T, q: T, sig: T, mat: T) -> T {
        let st = (sig * mat.sqrt()).max_val(T::from_f64(1e-6));
        let d1 = ((s / k).ln() + (r - q + T::from_f64(0.5) * sig * sig) * mat) / st;
        let d2 = d1 - st;
        s * (-q * mat).exp() * norm_cdf(d1) - k * (-r * mat).exp() * norm_cdf(d2)
    }

    fn down_and_out_call<T: Real>(s: T, k: T, b: T, r: T, q: T, sig: T, mat: T) -> T {
        let c_bs = bs_call_generic(s, k, r, q, sig, mat);
        let sig2 = (sig * sig).max_val(T::from_f64(1e-6));
        let lam = (r - q - T::from_f64(0.5) * sig2) / sig2;
        let refl_s = (b * b) / s;
        let factor = ((b / s).ln() * T::from_f64(2.0) * lam).exp();
        let c_refl = bs_call_generic(refl_s, k, r, q, sig, mat);
        (c_bs - factor * c_refl).max_val(T::from_f64(1e-6))
    }

    pub fn evaluate_broadie_glasserman_kou_discrete_barrier_lookback<T: Real>(
        spec: &BroadieGlassermanKouDiscreteBarrierSpec<T>,
    ) -> BroadieGlassermanKouDiscreteBarrierReport<T> {
        let beta1 = T::from_f64(0.582_597_157_939_010_7);
        let s = spec.spot_price.max_val(T::from_f64(1e-4));
        let k = spec.strike_price.max_val(T::from_f64(1e-4));
        let h_bar = spec.down_barrier_h.max_val(T::from_f64(1e-4));
        let sig = spec.volatility.max_val(T::from_f64(0.01));
        let mat = spec.maturity_years.max_val(T::from_f64(1e-4));
        let m_steps = spec.monitoring_steps_m.max_val(T::from_f64(1.0));
        let dt = mat / m_steps;
        let shift = (-beta1 * sig * dt.sqrt()).exp();
        let h_eff = h_bar * shift;

        let c_cont = down_and_out_call(
            s,
            k,
            h_bar,
            spec.risk_free_rate,
            spec.dividend_yield,
            sig,
            mat,
        );
        let c_disc = down_and_out_call(
            s,
            k,
            h_eff,
            spec.risk_free_rate,
            spec.dividend_yield,
            sig,
            mat,
        );

        let atm_call = bs_call_generic(s, s, spec.risk_free_rate, spec.dividend_yield, sig, mat);
        let lb_cont = atm_call + T::from_f64(0.5) * s * sig * mat.sqrt();
        let lb_disc = shift * lb_cont + (T::one() - shift) * atm_call;

        BroadieGlassermanKouDiscreteBarrierReport {
            shifted_effective_barrier: h_eff,
            continuous_down_and_out_call: c_cont,
            discrete_bgk_down_and_out_call: c_disc,
            continuity_correction_price_uplift: c_disc - c_cont,
            continuous_lookback_call: lb_cont,
            discrete_bgk_lookback_call: lb_disc,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct OhlsonPenmanLinearInfoDynamicsSpec<T: Real = f64> {
        pub current_book_value_bt: T,
        pub current_earnings_xt: T,
        pub current_dividends_dt: T,
        pub other_information_shock_nu: T,
        pub cost_of_equity_re: T,
        pub abnormal_earnings_persistence_omega: T,
        pub other_info_persistence_gamma: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct OhlsonPenmanLinearInfoDynamicsReport<T: Real = f64> {
        pub prior_clean_surplus_book_value: T,
        pub current_abnormal_earnings: T,
        pub abnormal_earnings_multiplier_alpha1: T,
        pub other_information_multiplier_alpha2: T,
        pub earnings_valuation_weight_k: T,
        pub lid_intrinsic_equity_value: T,
        pub weighted_earnings_book_value: T,
        pub ohlson_weighted_book_earnings_error: T,
        pub implied_price_to_book_ratio: T,
    }

    pub fn evaluate_ohlson_penman_linear_information_dynamics<T: Real>(
        spec: &OhlsonPenmanLinearInfoDynamicsSpec<T>,
    ) -> OhlsonPenmanLinearInfoDynamicsReport<T> {
        let bt = spec.current_book_value_bt.max_val(T::from_f64(1e-4));
        let xt = spec.current_earnings_xt;
        let dt = spec.current_dividends_dt;
        let nu = spec.other_information_shock_nu;
        let re = spec.cost_of_equity_re.max_val(T::from_f64(1e-4));
        let big_r = T::one() + re;
        let omega = spec.abnormal_earnings_persistence_omega;
        let gamma = spec.other_info_persistence_gamma;

        let b_tm1 = bt + dt - xt;
        let x_a = xt - re * b_tm1;

        let r_m_om = (big_r - omega).max_val(T::from_f64(1e-6));
        let r_m_gm = (big_r - gamma).max_val(T::from_f64(1e-6));
        let alpha1 = omega / r_m_om;
        let alpha2 = big_r / (r_m_om * r_m_gm);
        let k_wt = re * alpha1;
        let phi = big_r / re;

        let v_lid = bt + alpha1 * x_a + alpha2 * nu;
        let v_weighted = k_wt * (phi * xt - dt) + (T::one() - k_wt) * bt + alpha2 * nu;
        let d_eq = v_lid - v_weighted;
        let eq_err = (d_eq * d_eq + T::from_f64(1e-30)).sqrt();

        OhlsonPenmanLinearInfoDynamicsReport {
            prior_clean_surplus_book_value: b_tm1,
            current_abnormal_earnings: x_a,
            abnormal_earnings_multiplier_alpha1: alpha1,
            other_information_multiplier_alpha2: alpha2,
            earnings_valuation_weight_k: k_wt,
            lid_intrinsic_equity_value: v_lid,
            weighted_earnings_book_value: v_weighted,
            ohlson_weighted_book_earnings_error: eq_err,
            implied_price_to_book_ratio: v_lid / bt,
        }
    }

    #[cfg(test)]
    mod domain_suite_50_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bgk_barrier_and_ohlson_lid_invariants_and_ad_parity() {
            let h = 1e-6;

            let bgk_spec = BroadieGlassermanKouDiscreteBarrierSpec {
                spot_price: 100.0,
                strike_price: 100.0,
                down_barrier_h: 90.0,
                risk_free_rate: 0.04,
                dividend_yield: 0.01,
                volatility: 0.22,
                maturity_years: 1.0,
                monitoring_steps_m: 50.0,
            };
            let bgk_rep = evaluate_broadie_glasserman_kou_discrete_barrier_lookback(&bgk_spec);
            assert!(bgk_rep.shifted_effective_barrier < bgk_spec.down_barrier_h);
            assert!(bgk_rep.discrete_bgk_down_and_out_call > bgk_rep.continuous_down_and_out_call);
            let bgk_ad = evaluate_broadie_glasserman_kou_discrete_barrier_lookback(
                &BroadieGlassermanKouDiscreteBarrierSpec {
                    spot_price: Dual::<1>::variable(bgk_spec.spot_price, 0),
                    strike_price: Dual::<1>::from_value(bgk_spec.strike_price),
                    down_barrier_h: Dual::<1>::from_value(bgk_spec.down_barrier_h),
                    risk_free_rate: Dual::<1>::from_value(bgk_spec.risk_free_rate),
                    dividend_yield: Dual::<1>::from_value(bgk_spec.dividend_yield),
                    volatility: Dual::<1>::from_value(bgk_spec.volatility),
                    maturity_years: Dual::<1>::from_value(bgk_spec.maturity_years),
                    monitoring_steps_m: Dual::<1>::from_value(bgk_spec.monitoring_steps_m),
                },
            );
            let mut bp = bgk_spec;
            bp.spot_price += h;
            let mut bm = bgk_spec;
            bm.spot_price -= h;
            let fd_bgk = (evaluate_broadie_glasserman_kou_discrete_barrier_lookback(&bp)
                .discrete_bgk_down_and_out_call
                - evaluate_broadie_glasserman_kou_discrete_barrier_lookback(&bm)
                    .discrete_bgk_down_and_out_call)
                / (2.0 * h);
            assert!((bgk_ad.discrete_bgk_down_and_out_call.d[0] - fd_bgk).abs() < 1e-5);

            let ohl_spec = OhlsonPenmanLinearInfoDynamicsSpec {
                current_book_value_bt: 100.0,
                current_earnings_xt: 14.5,
                current_dividends_dt: 4.5,
                other_information_shock_nu: 1.8,
                cost_of_equity_re: 0.09,
                abnormal_earnings_persistence_omega: 0.62,
                other_info_persistence_gamma: 0.45,
            };
            let ohl_rep = evaluate_ohlson_penman_linear_information_dynamics(&ohl_spec);
            assert!(ohl_rep.ohlson_weighted_book_earnings_error < 1e-10);
            assert!(ohl_rep.abnormal_earnings_multiplier_alpha1 > 0.0);
            let ohl_ad = evaluate_ohlson_penman_linear_information_dynamics(
                &OhlsonPenmanLinearInfoDynamicsSpec {
                    current_book_value_bt: Dual::<1>::from_value(ohl_spec.current_book_value_bt),
                    current_earnings_xt: Dual::<1>::variable(ohl_spec.current_earnings_xt, 0),
                    current_dividends_dt: Dual::<1>::from_value(ohl_spec.current_dividends_dt),
                    other_information_shock_nu: Dual::<1>::from_value(
                        ohl_spec.other_information_shock_nu,
                    ),
                    cost_of_equity_re: Dual::<1>::from_value(ohl_spec.cost_of_equity_re),
                    abnormal_earnings_persistence_omega: Dual::<1>::from_value(
                        ohl_spec.abnormal_earnings_persistence_omega,
                    ),
                    other_info_persistence_gamma: Dual::<1>::from_value(
                        ohl_spec.other_info_persistence_gamma,
                    ),
                },
            );
            let mut op = ohl_spec;
            op.current_earnings_xt += h;
            let mut om = ohl_spec;
            om.current_earnings_xt -= h;
            let fd_ohl = (evaluate_ohlson_penman_linear_information_dynamics(&op)
                .lid_intrinsic_equity_value
                - evaluate_ohlson_penman_linear_information_dynamics(&om)
                    .lid_intrinsic_equity_value)
                / (2.0 * h);
            assert!((ohl_ad.lid_intrinsic_equity_value.d[0] - fd_ohl).abs() < 1e-5);
        }
    }
}
pub use evaluate_broadie_glasserman_kou_discrete_barrier_lookback_module::*;

mod evaluate_campbell_cochrane_external_habit_sdf_module {
    // Cochrane (2005) & Campbell (2018) Campbell–Cochrane (1999) External Habit-Formation
    // Surplus-Consumption SDF, Countercyclical Sharpe Ratio & Constant Risk-Free Rate Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct CampbellCochraneExternalHabitSpec<T: Real = f64> {
        pub utility_curvature_gamma: T,
        pub consumption_growth_mean_gc: T,
        pub consumption_growth_vol_sigmac: T,
        pub log_surplus_persistence_phis: T,
        pub subjective_discount_rate_rho_time: T,
        pub current_log_surplus_deviation_ds: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct CampbellCochraneExternalHabitReport<T: Real = f64> {
        pub steady_state_surplus_ratio_sbar: T,
        pub max_log_surplus_smax: T,
        pub sensitivity_function_lambda_s: T,
        pub local_effective_risk_aversion: T,
        pub real_risk_free_rate_rf: T,
        pub conditional_max_sharpe_ratio: T,
        pub conditional_equity_risk_premium: T,
    }

    pub fn evaluate_campbell_cochrane_external_habit_sdf<T: Real>(
        spec: &CampbellCochraneExternalHabitSpec<T>,
    ) -> CampbellCochraneExternalHabitReport<T> {
        let gamma = spec.utility_curvature_gamma.max_val(T::from_f64(0.5));
        let sig_c = spec
            .consumption_growth_vol_sigmac
            .max_val(T::from_f64(1e-4));
        let one_minus_phi =
            (T::one() - spec.log_surplus_persistence_phis).max_val(T::from_f64(1e-4));

        // Steady-state surplus consumption ratio S_bar = sigma_c * sqrt(gamma / (1 - phi_s))
        let s_bar = sig_c * (gamma / one_minus_phi).sqrt();
        let s_bar_log = s_bar.ln();
        let s_max = s_bar_log + T::from_f64(0.5) * (T::one() - s_bar * s_bar);

        // Sensitivity function lambda(s_t) = (1 / S_bar) * sqrt(1 - 2 * (s_t - s_bar)) - 1
        let inside = (T::one() - T::from_f64(2.0) * spec.current_log_surplus_deviation_ds)
            .max_val(T::from_f64(1e-6));
        let lambda_s = inside.sqrt() / s_bar - T::one();

        // Local effective relative risk aversion gamma_t = gamma * (1 + lambda(s_t))
        let local_rra = gamma * (T::one() + lambda_s);

        // Exact constant real risk-free rate under Campbell-Cochrane baseline (B = 0):
        // r_f = rho_time + gamma * g_c - 0.5 * gamma * (1 - phi_s)
        let rf = spec.subjective_discount_rate_rho_time + gamma * spec.consumption_growth_mean_gc
            - T::from_f64(0.5) * gamma * one_minus_phi;

        // Conditional maximum Sharpe ratio SR_t = gamma * sigma_c * (1 + lambda(s_t))
        let sr_t = local_rra * sig_c;
        // Conditional equity risk premium for claim with return volatility ~ sigma_c * (1 + lambda(s_t)):
        let erp_t = sr_t * sig_c * (T::one() + lambda_s);

        CampbellCochraneExternalHabitReport {
            steady_state_surplus_ratio_sbar: s_bar,
            max_log_surplus_smax: s_max,
            sensitivity_function_lambda_s: lambda_s,
            local_effective_risk_aversion: local_rra,
            real_risk_free_rate_rf: rf,
            conditional_max_sharpe_ratio: sr_t,
            conditional_equity_risk_premium: erp_t,
        }
    }

    #[cfg(test)]
    mod domain_suite_51_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn campbell_cochrane_external_habit_invariants_and_ad_parity() {
            let boom = CampbellCochraneExternalHabitSpec {
                utility_curvature_gamma: 2.0,
                consumption_growth_mean_gc: 0.0189,
                consumption_growth_vol_sigmac: 0.015,
                log_surplus_persistence_phis: 0.87,
                subjective_discount_rate_rho_time: 0.0116,
                current_log_surplus_deviation_ds: 0.05,
            };
            let mut recession = boom;
            recession.current_log_surplus_deviation_ds = -0.15;

            let r_boom = evaluate_campbell_cochrane_external_habit_sdf(&boom);
            let r_rec = evaluate_campbell_cochrane_external_habit_sdf(&recession);
            // Countercyclical risk premium: Sharpe ratio and ERP rise in recessions (ds < 0)
            // while real risk-free rate remains invariant to surplus consumption ds!
            assert!(r_rec.conditional_max_sharpe_ratio > r_boom.conditional_max_sharpe_ratio);
            assert!(r_rec.conditional_equity_risk_premium > r_boom.conditional_equity_risk_premium);
            assert!((r_rec.real_risk_free_rate_rf - r_boom.real_risk_free_rate_rf).abs() < 1e-12);

            let h = 1e-6;
            let r_ad =
                evaluate_campbell_cochrane_external_habit_sdf(&CampbellCochraneExternalHabitSpec {
                    utility_curvature_gamma: Dual::<1>::from_value(
                        recession.utility_curvature_gamma,
                    ),
                    consumption_growth_mean_gc: Dual::<1>::from_value(
                        recession.consumption_growth_mean_gc,
                    ),
                    consumption_growth_vol_sigmac: Dual::<1>::from_value(
                        recession.consumption_growth_vol_sigmac,
                    ),
                    log_surplus_persistence_phis: Dual::<1>::from_value(
                        recession.log_surplus_persistence_phis,
                    ),
                    subjective_discount_rate_rho_time: Dual::<1>::from_value(
                        recession.subjective_discount_rate_rho_time,
                    ),
                    current_log_surplus_deviation_ds: Dual::<1>::variable(
                        recession.current_log_surplus_deviation_ds,
                        0,
                    ),
                });
            let mut sp = recession;
            sp.current_log_surplus_deviation_ds += h;
            let mut sm = recession;
            sm.current_log_surplus_deviation_ds -= h;
            let fd = (evaluate_campbell_cochrane_external_habit_sdf(&sp)
                .conditional_max_sharpe_ratio
                - evaluate_campbell_cochrane_external_habit_sdf(&sm).conditional_max_sharpe_ratio)
                / (2.0 * h);
            assert!((r_ad.conditional_max_sharpe_ratio.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_campbell_cochrane_external_habit_sdf_module::*;

mod evaluate_bansal_yaron_long_run_risks_sdf_module {
    // Cochrane (2005) & Bansal–Yaron (2004) Long-Run Risks (`LRR`) Epstein–Zin Recursive Utility SDF,
    // Persistent Growth A_1 & Stochastic Volatility A_2 Equity Premium Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BansalYaronLongRunRisksSpec<T: Real = f64> {
        pub relative_risk_aversion_gamma: T,
        pub elasticity_intertemporal_sub_psi: T,
        pub campbell_shiller_kappa1: T,
        pub growth_persistence_rho_x: T,
        pub volatility_persistence_nu_sigma: T,
        pub growth_shock_scale_phi_e: T,
        pub baseline_consumption_vol_sigma: T,
        pub vol_of_vol_sigma_w: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BansalYaronLongRunRisksReport<T: Real = f64> {
        pub epstein_zin_theta: T,
        pub wealth_growth_loading_a1: T,
        pub wealth_volatility_loading_a2: T,
        pub market_price_of_short_run_risk_lambda_m: T,
        pub market_price_of_long_run_growth_risk_lambda_e: T,
        pub market_price_of_volatility_risk_lambda_w: T,
        pub total_consumption_claim_risk_premium: T,
    }

    pub fn evaluate_bansal_yaron_long_run_risks_sdf<T: Real>(
        spec: &BansalYaronLongRunRisksSpec<T>,
    ) -> BansalYaronLongRunRisksReport<T> {
        let gamma = spec.relative_risk_aversion_gamma;
        let psi = spec
            .elasticity_intertemporal_sub_psi
            .max_val(T::from_f64(0.1));
        let inv_psi = T::one() / psi;

        // Epstein-Zin preference parameter theta = (1 - gamma) / (1 - 1/psi):
        let theta = (T::one() - gamma) / (T::one() - inv_psi);
        let k1 = spec.campbell_shiller_kappa1;

        // Bansal-Yaron (2004, Eq. 9-10) wealth-consumption loadings:
        // A_1 = (1 - 1/psi) / (1 - kappa_1 * rho_x)
        let a1 = (T::one() - inv_psi) / (T::one() - k1 * spec.growth_persistence_rho_x);

        // A_2 = 0.5 * ((theta - theta/psi)^2 + (theta * kappa_1 * A_1 * phi_e)^2) / (theta * (1 - kappa_1 * nu_sigma))
        let term_c = theta - theta * inv_psi;
        let term_x = theta * k1 * a1 * spec.growth_shock_scale_phi_e;
        let a2 = T::from_f64(0.5) * (term_c * term_c + term_x * term_x)
            / (theta * (T::one() - k1 * spec.volatility_persistence_nu_sigma));

        // Market prices of short-run consumption shock (eta), long-run growth shock (e), and vol shock (w):
        let lam_m = gamma;
        let lam_e = (T::one() - theta) * k1 * a1 * spec.growth_shock_scale_phi_e;
        let lam_w = (T::one() - theta) * k1 * a2 * spec.vol_of_vol_sigma_w;

        // Consumption-wealth claim betas and total risk premium:
        let beta_m_e = k1 * a1 * spec.growth_shock_scale_phi_e;
        let beta_m_w = k1 * a2 * spec.vol_of_vol_sigma_w;
        let sig2 = spec.baseline_consumption_vol_sigma * spec.baseline_consumption_vol_sigma;
        let rp = (lam_m + lam_e * beta_m_e) * sig2 + lam_w * beta_m_w;

        BansalYaronLongRunRisksReport {
            epstein_zin_theta: theta,
            wealth_growth_loading_a1: a1,
            wealth_volatility_loading_a2: a2,
            market_price_of_short_run_risk_lambda_m: lam_m,
            market_price_of_long_run_growth_risk_lambda_e: lam_e,
            market_price_of_volatility_risk_lambda_w: lam_w,
            total_consumption_claim_risk_premium: rp,
        }
    }

    #[cfg(test)]
    mod domain_suite_52_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bansal_yaron_long_run_risks_invariants_and_ad_parity() {
            let spec = BansalYaronLongRunRisksSpec {
                relative_risk_aversion_gamma: 7.5,
                elasticity_intertemporal_sub_psi: 1.5,
                campbell_shiller_kappa1: 0.997,
                growth_persistence_rho_x: 0.979,
                volatility_persistence_nu_sigma: 0.987,
                growth_shock_scale_phi_e: 0.044,
                baseline_consumption_vol_sigma: 0.027,
                vol_of_vol_sigma_w: 0.00023,
            };
            let rep = evaluate_bansal_yaron_long_run_risks_sdf(&spec);
            // For gamma > 1 and psi > 1 (preference for early resolution of uncertainty):
            // theta < 0, A_1 > 0, A_2 < 0, and long-run growth risk commands a positive risk price!
            assert!(rep.epstein_zin_theta < 0.0);
            assert!(rep.wealth_growth_loading_a1 > 0.0);
            assert!(rep.wealth_volatility_loading_a2 < 0.0);
            assert!(rep.market_price_of_long_run_growth_risk_lambda_e > 0.0);
            assert!(rep.total_consumption_claim_risk_premium > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_bansal_yaron_long_run_risks_sdf(&BansalYaronLongRunRisksSpec {
                relative_risk_aversion_gamma: Dual::<1>::variable(
                    spec.relative_risk_aversion_gamma,
                    0,
                ),
                elasticity_intertemporal_sub_psi: Dual::<1>::from_value(
                    spec.elasticity_intertemporal_sub_psi,
                ),
                campbell_shiller_kappa1: Dual::<1>::from_value(spec.campbell_shiller_kappa1),
                growth_persistence_rho_x: Dual::<1>::from_value(spec.growth_persistence_rho_x),
                volatility_persistence_nu_sigma: Dual::<1>::from_value(
                    spec.volatility_persistence_nu_sigma,
                ),
                growth_shock_scale_phi_e: Dual::<1>::from_value(spec.growth_shock_scale_phi_e),
                baseline_consumption_vol_sigma: Dual::<1>::from_value(
                    spec.baseline_consumption_vol_sigma,
                ),
                vol_of_vol_sigma_w: Dual::<1>::from_value(spec.vol_of_vol_sigma_w),
            });
            let mut sp = spec;
            sp.relative_risk_aversion_gamma += h;
            let mut sm = spec;
            sm.relative_risk_aversion_gamma -= h;
            let fd = (evaluate_bansal_yaron_long_run_risks_sdf(&sp)
                .total_consumption_claim_risk_premium
                - evaluate_bansal_yaron_long_run_risks_sdf(&sm)
                    .total_consumption_claim_risk_premium)
                / (2.0 * h);
            assert!((r_ad.total_consumption_claim_risk_premium.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_bansal_yaron_long_run_risks_sdf_module::*;

mod evaluate_gabaix_rare_disaster_resilience_module {
    // Gabaix (2012) & Barro (2006) Time-Varying Rare Disaster Resilience (`Linearity-Generating`
    // Process), Equity Risk Premium & Disaster Put Skew Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct GabaixRareDisasterResilienceSpec<T: Real = f64> {
        pub disaster_probability_p: T,
        pub consumption_recovery_in_disaster_b: T,
        pub dividend_recovery_in_disaster_f: T,
        pub crra_risk_aversion_gamma: T,
        pub time_preference_rate_rho: T,
        pub normal_growth_rate_g: T,
        pub resilience_deviation_h_hat: T,
        pub resilience_mean_reversion_phi_h: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GabaixRareDisasterResilienceReport<T: Real = f64> {
        pub baseline_equity_resilience_h_star: T,
        pub effective_equity_discount_rate_delta_d: T,
        pub linearity_generating_price_dividend_ratio: T,
        pub disaster_depressed_risk_free_rate: T,
        pub conditional_equity_risk_premium: T,
        pub risk_neutral_disaster_intensity: T,
        pub disaster_jump_variance_skew_premium: T,
    }

    pub fn evaluate_gabaix_rare_disaster_resilience<T: Real>(
        spec: &GabaixRareDisasterResilienceSpec<T>,
    ) -> GabaixRareDisasterResilienceReport<T> {
        let p = spec.disaster_probability_p.max_val(T::from_f64(1e-4));
        let bc = spec
            .consumption_recovery_in_disaster_b
            .max_val(T::from_f64(0.20))
            .min_val(T::from_f64(0.99));
        let fd = spec
            .dividend_recovery_in_disaster_f
            .max_val(T::from_f64(0.10))
            .min_val(T::from_f64(0.99));
        let gam = spec.crra_risk_aversion_gamma.max_val(T::from_f64(1.0));
        let rho_pref = spec.time_preference_rate_rho.max_val(T::from_f64(0.002));
        let g = spec.normal_growth_rate_g;
        let h_hat = spec.resilience_deviation_h_hat;
        let phi_h = spec
            .resilience_mean_reversion_phi_h
            .max_val(T::from_f64(0.02));

        let sdf_disaster_mult = bc.pow(-gam);
        // Gabaix (2012, QJE Eq. 11 & 14) baseline equity resilience & discount rate:
        let h_star = p * (sdf_disaster_mult * fd - T::one());
        let delta_d = (rho_pref + (gam - T::one()) * g - h_star).max_val(T::from_f64(0.01));

        // Gabaix Theorem 1 (Eq. 16) exact Linearity-Generating price-dividend ratio:
        let pd_ratio = (T::one() + h_hat / (delta_d + phi_h)) / delta_d;

        // Barro (2006) / Gabaix precautionary-savings risk-free rate & conditional ERP:
        let rf = rho_pref + gam * g - p * (sdf_disaster_mult - T::one());
        let erp = p * (sdf_disaster_mult - T::one()) * (T::one() - T::from_f64(0.5) * fd)
            - h_hat * (delta_d / (delta_d + phi_h));

        let p_rn = p * sdf_disaster_mult;
        let log_fd = fd.ln();
        let skew_var_prem = (p_rn - p) * log_fd * log_fd;

        GabaixRareDisasterResilienceReport {
            baseline_equity_resilience_h_star: h_star,
            effective_equity_discount_rate_delta_d: delta_d,
            linearity_generating_price_dividend_ratio: pd_ratio,
            disaster_depressed_risk_free_rate: rf,
            conditional_equity_risk_premium: erp,
            risk_neutral_disaster_intensity: p_rn,
            disaster_jump_variance_skew_premium: skew_var_prem,
        }
    }

    #[cfg(test)]
    mod domain_suite_55_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gabaix_rare_disaster_resilience_invariants_and_ad_parity() {
            let spec = GabaixRareDisasterResilienceSpec {
                disaster_probability_p: 0.018,
                consumption_recovery_in_disaster_b: 0.72,
                dividend_recovery_in_disaster_f: 0.55,
                crra_risk_aversion_gamma: 4.0,
                time_preference_rate_rho: 0.025,
                normal_growth_rate_g: 0.020,
                resilience_deviation_h_hat: -0.008,
                resilience_mean_reversion_phi_h: 0.20,
            };
            let rep = evaluate_gabaix_rare_disaster_resilience(&spec);
            assert!(rep.linearity_generating_price_dividend_ratio > 10.0);
            assert!(rep.conditional_equity_risk_premium > 0.03);
            assert!(rep.risk_neutral_disaster_intensity > spec.disaster_probability_p);

            let h = 1e-6;
            let r_ad =
                evaluate_gabaix_rare_disaster_resilience(&GabaixRareDisasterResilienceSpec {
                    disaster_probability_p: Dual::<1>::variable(spec.disaster_probability_p, 0),
                    consumption_recovery_in_disaster_b: Dual::<1>::from_value(
                        spec.consumption_recovery_in_disaster_b,
                    ),
                    dividend_recovery_in_disaster_f: Dual::<1>::from_value(
                        spec.dividend_recovery_in_disaster_f,
                    ),
                    crra_risk_aversion_gamma: Dual::<1>::from_value(spec.crra_risk_aversion_gamma),
                    time_preference_rate_rho: Dual::<1>::from_value(spec.time_preference_rate_rho),
                    normal_growth_rate_g: Dual::<1>::from_value(spec.normal_growth_rate_g),
                    resilience_deviation_h_hat: Dual::<1>::from_value(
                        spec.resilience_deviation_h_hat,
                    ),
                    resilience_mean_reversion_phi_h: Dual::<1>::from_value(
                        spec.resilience_mean_reversion_phi_h,
                    ),
                });
            let mut sp = spec;
            sp.disaster_probability_p += h;
            let mut sm = spec;
            sm.disaster_probability_p -= h;
            let fd = (evaluate_gabaix_rare_disaster_resilience(&sp)
                .conditional_equity_risk_premium
                - evaluate_gabaix_rare_disaster_resilience(&sm).conditional_equity_risk_premium)
                / (2.0 * h);
            assert!((r_ad.conditional_equity_risk_premium.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_gabaix_rare_disaster_resilience_module::*;

mod evaluate_laubach_williams_natural_rate_r_star_module {
    // Laubach & Williams (2003) & Holston–Laubach–Williams (`HLW`, 2017) Natural Rate of
    // Interest (`r^*`), Output-Gap Kalman Filter & Real Monetary Stance Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct LaubachWilliamsNaturalRateRStarSpec<T: Real = f64> {
        pub prior_trend_growth_g0: T,
        pub prior_other_factor_z0: T,
        pub nominal_policy_rate_i: T,
        pub expected_inflation_pi_e: T,
        pub realized_core_inflation_pi: T,
        pub is_curve_surprise_dy: T,
        pub is_rate_sensitivity_a_r: T,
        pub phillips_slope_b_y: T,
        pub signal_to_noise_ratio_lambda_z: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct LaubachWilliamsNaturalRateRStarReport<T: Real = f64> {
        pub prior_natural_real_rate_r_star: T,
        pub posterior_natural_real_rate_r_star: T,
        pub ex_ante_real_policy_rate: T,
        pub real_monetary_policy_stance_gap: T,
        pub filtered_output_gap_pct: T,
        pub kalman_gain_on_r_star: T,
        pub taylor_rule_implied_nominal_rate: T,
    }

    pub fn evaluate_laubach_williams_natural_rate_r_star<T: Real>(
        spec: &LaubachWilliamsNaturalRateRStarSpec<T>,
    ) -> LaubachWilliamsNaturalRateRStarReport<T> {
        let g0 = spec.prior_trend_growth_g0;
        let z0 = spec.prior_other_factor_z0;
        let r_star_prior = g0 + z0;

        let real_rate = spec.nominal_policy_rate_i - spec.expected_inflation_pi_e;
        let a_r = spec.is_rate_sensitivity_a_r.max_val(T::from_f64(0.02));
        let b_y = spec.phillips_slope_b_y.max_val(T::from_f64(0.02));
        let lam_z = spec
            .signal_to_noise_ratio_lambda_z
            .max_val(T::from_f64(0.005));

        // Steady-state Riccati Kalman gain K_z from signal-to-noise ratio lambda_z = a_r * sigma_z / sigma_y:
        let p_ss = T::from_f64(0.5)
            * (lam_z * lam_z
                + (lam_z * lam_z * lam_z * lam_z + T::from_f64(4.0) * lam_z * lam_z).sqrt());
        let k_is = (a_r * p_ss) / (T::one() + a_r * a_r * p_ss + b_y * b_y * p_ss);
        let k_pc = (b_y * p_ss) / (T::one() + a_r * a_r * p_ss + b_y * b_y * p_ss);

        // IS-curve and Phillips-curve innovations:
        let pi_surprise = spec.realized_core_inflation_pi - spec.expected_inflation_pi_e;
        let is_implied_gap = spec.is_curve_surprise_dy + a_r * (real_rate - r_star_prior);
        let filtered_y_gap =
            T::from_f64(0.65) * is_implied_gap + T::from_f64(0.35) * (pi_surprise / b_y);

        // Posterior natural rate r*_t|t updated positively when output/inflation exceed IS prior:
        let r_star_post = r_star_prior + k_is * is_implied_gap + k_pc * pi_surprise;
        let stance_gap = real_rate - r_star_post;

        // Laubach-Williams augmented Taylor rule: i_TR = r*_post + pi + 0.5 * (pi - 0.02) + 0.5 * y_gap
        let taylor_nom = r_star_post
            + spec.realized_core_inflation_pi
            + T::from_f64(0.5) * (spec.realized_core_inflation_pi - T::from_f64(0.02))
            + T::from_f64(0.5) * filtered_y_gap;

        LaubachWilliamsNaturalRateRStarReport {
            prior_natural_real_rate_r_star: r_star_prior,
            posterior_natural_real_rate_r_star: r_star_post,
            ex_ante_real_policy_rate: real_rate,
            real_monetary_policy_stance_gap: stance_gap,
            filtered_output_gap_pct: filtered_y_gap * T::from_f64(100.0),
            kalman_gain_on_r_star: k_is,
            taylor_rule_implied_nominal_rate: taylor_nom,
        }
    }

    #[cfg(test)]
    mod evaluate_laubach_williams_natural_rate_r_star_module_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn laubach_williams_natural_rate_r_star_invariants_and_ad_parity() {
            let spec = LaubachWilliamsNaturalRateRStarSpec {
                prior_trend_growth_g0: 0.018,
                prior_other_factor_z0: -0.004,
                nominal_policy_rate_i: 0.0525,
                expected_inflation_pi_e: 0.024,
                realized_core_inflation_pi: 0.028,
                is_curve_surprise_dy: 0.006,
                is_rate_sensitivity_a_r: 0.09,
                phillips_slope_b_y: 0.15,
                signal_to_noise_ratio_lambda_z: 0.042,
            };
            let rep = evaluate_laubach_williams_natural_rate_r_star(&spec);
            assert!(rep.posterior_natural_real_rate_r_star > rep.prior_natural_real_rate_r_star);
            assert!(rep.real_monetary_policy_stance_gap > 0.0);
            assert!(rep.kalman_gain_on_r_star > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_laubach_williams_natural_rate_r_star(
                &LaubachWilliamsNaturalRateRStarSpec {
                    prior_trend_growth_g0: Dual::<1>::variable(spec.prior_trend_growth_g0, 0),
                    prior_other_factor_z0: Dual::<1>::from_value(spec.prior_other_factor_z0),
                    nominal_policy_rate_i: Dual::<1>::from_value(spec.nominal_policy_rate_i),
                    expected_inflation_pi_e: Dual::<1>::from_value(spec.expected_inflation_pi_e),
                    realized_core_inflation_pi: Dual::<1>::from_value(
                        spec.realized_core_inflation_pi,
                    ),
                    is_curve_surprise_dy: Dual::<1>::from_value(spec.is_curve_surprise_dy),
                    is_rate_sensitivity_a_r: Dual::<1>::from_value(spec.is_rate_sensitivity_a_r),
                    phillips_slope_b_y: Dual::<1>::from_value(spec.phillips_slope_b_y),
                    signal_to_noise_ratio_lambda_z: Dual::<1>::from_value(
                        spec.signal_to_noise_ratio_lambda_z,
                    ),
                },
            );
            let mut sp = spec;
            sp.prior_trend_growth_g0 += h;
            let mut sm = spec;
            sm.prior_trend_growth_g0 -= h;
            let fd = (evaluate_laubach_williams_natural_rate_r_star(&sp)
                .taylor_rule_implied_nominal_rate
                - evaluate_laubach_williams_natural_rate_r_star(&sm)
                    .taylor_rule_implied_nominal_rate)
                / (2.0 * h);
            assert!((r_ad.taylor_rule_implied_nominal_rate.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_laubach_williams_natural_rate_r_star_module::*;
