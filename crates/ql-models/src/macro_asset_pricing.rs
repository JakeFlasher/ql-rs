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

mod evaluate_epstein_zin_weil_recursive_utility_sdf_module {
    // Epstein & Zin (1989, Econometrica) & Weil (1989) Recursive Utility Stochastic Discount
    // Factor (`SDF`), Intertemporal Elasticity of Substitution (`IES psi`) vs. Risk Aversion Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct EpsteinZinWeilRecursiveUtilitySdfSpec<T: Real = f64> {
        pub subjective_discount_beta: T,
        pub relative_risk_aversion_gamma: T,
        pub intertemporal_elasticity_psi: T,
        pub expected_consumption_growth_gc: T,
        pub consumption_growth_volatility_sigma_c: T,
        pub wealth_portfolio_return_volatility_sigma_w: T,
        pub consumption_wealth_correlation_rho_cw: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct EpsteinZinWeilRecursiveUtilitySdfReport<T: Real = f64> {
        pub epstein_zin_theta_index: T,
        pub early_uncertainty_resolution_wedge: T,
        pub equilibrium_real_risk_free_rate: T,
        pub consumption_growth_risk_price: T,
        pub wealth_return_continuation_risk_price: T,
        pub log_sdf_variance: T,
        pub hansen_jagannathan_max_sharpe_bound: T,
    }

    pub fn evaluate_epstein_zin_weil_recursive_utility_sdf<T: Real>(
        spec: &EpsteinZinWeilRecursiveUtilitySdfSpec<T>,
    ) -> EpsteinZinWeilRecursiveUtilitySdfReport<T> {
        let beta = spec
            .subjective_discount_beta
            .max_val(T::from_f64(0.90))
            .min_val(T::from_f64(0.999));
        let gam = spec.relative_risk_aversion_gamma.max_val(T::from_f64(0.5));
        let psi = spec.intertemporal_elasticity_psi.max_val(T::from_f64(0.15));
        let gc = spec.expected_consumption_growth_gc;
        let sig_c = spec
            .consumption_growth_volatility_sigma_c
            .max_val(T::from_f64(1e-4));
        let sig_w = spec
            .wealth_portfolio_return_volatility_sigma_w
            .max_val(T::from_f64(1e-4));
        let rho = spec
            .consumption_wealth_correlation_rho_cw
            .max_val(T::from_f64(-0.98))
            .min_val(T::from_f64(0.98));

        let inv_psi = T::one() / psi;
        let denom_psi = T::one() - inv_psi;
        // Smooth protection away from psi == 1 without flat-derivative clamp:
        let sign_d = if denom_psi >= T::zero() {
            T::one()
        } else {
            -T::one()
        };
        let safe_denom = sign_d * (denom_psi * denom_psi + T::from_f64(1e-8)).sqrt();
        let theta = (T::one() - gam) / safe_denom;
        let early_res_wedge = gam - inv_psi;

        // Epstein-Zin log SDF: m_{t+1} = theta ln(beta) - (theta / psi) dc_{t+1} + (theta - 1) r_{w, t+1}
        let lambda_c = theta * inv_psi;
        let lambda_w = T::one() - theta;

        let var_m = lambda_c * lambda_c * sig_c * sig_c
            + lambda_w * lambda_w * sig_w * sig_w
            + T::from_f64(2.0) * lambda_c * lambda_w * rho * sig_c * sig_w;

        // Equilibrium real risk-free rate under Epstein-Zin preferences:
        let rf_real = -beta.ln() + inv_psi * gc
            - T::from_f64(0.5) * theta * inv_psi * inv_psi * sig_c * sig_c
            - T::from_f64(0.5) * (T::one() - theta) * sig_w * sig_w;
        let hj_bound = var_m.sqrt();

        EpsteinZinWeilRecursiveUtilitySdfReport {
            epstein_zin_theta_index: theta,
            early_uncertainty_resolution_wedge: early_res_wedge,
            equilibrium_real_risk_free_rate: rf_real,
            consumption_growth_risk_price: lambda_c,
            wealth_return_continuation_risk_price: lambda_w,
            log_sdf_variance: var_m,
            hansen_jagannathan_max_sharpe_bound: hj_bound,
        }
    }

    #[cfg(test)]
    mod evaluate_epstein_zin_weil_recursive_utility_sdf_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn epstein_zin_weil_recursive_utility_sdf_invariants_and_ad_parity() {
            let spec = EpsteinZinWeilRecursiveUtilitySdfSpec {
                subjective_discount_beta: 0.988,
                relative_risk_aversion_gamma: 7.5,
                intertemporal_elasticity_psi: 1.5,
                expected_consumption_growth_gc: 0.018,
                consumption_growth_volatility_sigma_c: 0.022,
                wealth_portfolio_return_volatility_sigma_w: 0.145,
                consumption_wealth_correlation_rho_cw: 0.65,
            };
            let rep = evaluate_epstein_zin_weil_recursive_utility_sdf(&spec);
            assert!(rep.epstein_zin_theta_index < 0.0);
            assert!(rep.early_uncertainty_resolution_wedge > 0.0);
            assert!(rep.hansen_jagannathan_max_sharpe_bound > 0.20);

            let h = 1e-6;
            let r_ad = evaluate_epstein_zin_weil_recursive_utility_sdf(
                &EpsteinZinWeilRecursiveUtilitySdfSpec {
                    subjective_discount_beta: Dual::<7>::variable(spec.subjective_discount_beta, 0),
                    relative_risk_aversion_gamma: Dual::<7>::variable(
                        spec.relative_risk_aversion_gamma,
                        1,
                    ),
                    intertemporal_elasticity_psi: Dual::<7>::variable(
                        spec.intertemporal_elasticity_psi,
                        2,
                    ),
                    expected_consumption_growth_gc: Dual::<7>::variable(
                        spec.expected_consumption_growth_gc,
                        3,
                    ),
                    consumption_growth_volatility_sigma_c: Dual::<7>::variable(
                        spec.consumption_growth_volatility_sigma_c,
                        4,
                    ),
                    wealth_portfolio_return_volatility_sigma_w: Dual::<7>::variable(
                        spec.wealth_portfolio_return_volatility_sigma_w,
                        5,
                    ),
                    consumption_wealth_correlation_rho_cw: Dual::<7>::variable(
                        spec.consumption_wealth_correlation_rho_cw,
                        6,
                    ),
                },
            );
            for (idx, &ad_rf) in r_ad.equilibrium_real_risk_free_rate.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.subjective_discount_beta += h;
                        sm.subjective_discount_beta -= h;
                    }
                    1 => {
                        sp.relative_risk_aversion_gamma += h;
                        sm.relative_risk_aversion_gamma -= h;
                    }
                    2 => {
                        sp.intertemporal_elasticity_psi += h;
                        sm.intertemporal_elasticity_psi -= h;
                    }
                    3 => {
                        sp.expected_consumption_growth_gc += h;
                        sm.expected_consumption_growth_gc -= h;
                    }
                    4 => {
                        sp.consumption_growth_volatility_sigma_c += h;
                        sm.consumption_growth_volatility_sigma_c -= h;
                    }
                    5 => {
                        sp.wealth_portfolio_return_volatility_sigma_w += h;
                        sm.wealth_portfolio_return_volatility_sigma_w -= h;
                    }
                    _ => {
                        sp.consumption_wealth_correlation_rho_cw += h;
                        sm.consumption_wealth_correlation_rho_cw -= h;
                    }
                }
                let metric_p = evaluate_epstein_zin_weil_recursive_utility_sdf(&sp);
                let metric_m = evaluate_epstein_zin_weil_recursive_utility_sdf(&sm);
                let fd = ((metric_p.equilibrium_real_risk_free_rate
                    + metric_p.hansen_jagannathan_max_sharpe_bound)
                    - (metric_m.equilibrium_real_risk_free_rate
                        + metric_m.hansen_jagannathan_max_sharpe_bound))
                    / (2.0 * h);
                let ad = ad_rf + r_ad.hansen_jagannathan_max_sharpe_bound.d[idx];
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_epstein_zin_weil_recursive_utility_sdf_module::*;

mod evaluate_lettau_ludvigson_cay_consumption_wealth_forecast_module {
    // Lettau & Ludvigson (2001, J. Finance) & Cochrane (2005, Ch. 20) Cointegrated
    // Consumption-Wealth Ratio (cay), Human-Capital Labor Share & Conditional Equity Premium Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct LettauLudvigsonCayForecastSpec<T: Real = f64> {
        pub log_consumption_c: T,
        pub log_asset_wealth_a: T,
        pub log_labor_income_y: T,
        pub asset_wealth_cointegrating_share_omega: T,
        pub unconditional_equity_premium_gamma0: T,
        pub cay_predictive_slope_beta_cay: T,
        pub conditional_market_volatility_sigma_m: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct LettauLudvigsonCayForecastReport<T: Real = f64> {
        pub human_capital_cointegrating_share: T,
        pub cointegrating_residual_cay: T,
        pub forecasted_annual_equity_risk_premium: T,
        pub conditional_market_sharpe_ratio: T,
        pub kelly_optimal_equity_allocation: T,
    }

    pub fn evaluate_lettau_ludvigson_cay_consumption_wealth_forecast<T: Real>(
        spec: &LettauLudvigsonCayForecastSpec<T>,
    ) -> LettauLudvigsonCayForecastReport<T> {
        let c = spec.log_consumption_c;
        let a = spec.log_asset_wealth_a;
        let y = spec.log_labor_income_y;
        let omega = spec
            .asset_wealth_cointegrating_share_omega
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.95));
        let gam0 = spec.unconditional_equity_premium_gamma0;
        let beta_cay = spec.cay_predictive_slope_beta_cay;
        let sig_m = spec
            .conditional_market_volatility_sigma_m
            .max_val(T::from_f64(0.05));

        let human_share = T::one() - omega;
        // Lettau & Ludvigson (2001, Eq. 6): cay_t = c_t - omega * a_t - (1 - omega) * y_t:
        let cay = c - omega * a - human_share * y;

        // Predictive regression for next-year excess stock market return:
        let erp = gam0 + beta_cay * cay;
        let cond_sr = erp / sig_m;
        let kelly_w = erp / (sig_m * sig_m * T::from_f64(3.0));

        LettauLudvigsonCayForecastReport {
            human_capital_cointegrating_share: human_share,
            cointegrating_residual_cay: cay,
            forecasted_annual_equity_risk_premium: erp,
            conditional_market_sharpe_ratio: cond_sr,
            kelly_optimal_equity_allocation: kelly_w,
        }
    }

    #[cfg(test)]
    mod evaluate_lettau_ludvigson_cay_consumption_wealth_forecast_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn lettau_ludvigson_cay_consumption_wealth_forecast_invariants_and_ad_parity() {
            let spec = LettauLudvigsonCayForecastSpec {
                log_consumption_c: 9.865,
                log_asset_wealth_a: 11.20,
                log_labor_income_y: 9.25,
                asset_wealth_cointegrating_share_omega: 0.31,
                unconditional_equity_premium_gamma0: 0.062,
                cay_predictive_slope_beta_cay: 2.15,
                conditional_market_volatility_sigma_m: 0.16,
            };
            let rep = evaluate_lettau_ludvigson_cay_consumption_wealth_forecast(&spec);
            assert!(rep.cointegrating_residual_cay > 0.0);
            assert!(
                rep.forecasted_annual_equity_risk_premium
                    > spec.unconditional_equity_premium_gamma0
            );
            assert!(rep.conditional_market_sharpe_ratio > 0.3);

            let h = 1e-5;
            let r_ad = evaluate_lettau_ludvigson_cay_consumption_wealth_forecast(
                &LettauLudvigsonCayForecastSpec {
                    log_consumption_c: Dual::<7>::variable(spec.log_consumption_c, 0),
                    log_asset_wealth_a: Dual::<7>::variable(spec.log_asset_wealth_a, 1),
                    log_labor_income_y: Dual::<7>::variable(spec.log_labor_income_y, 2),
                    asset_wealth_cointegrating_share_omega: Dual::<7>::variable(
                        spec.asset_wealth_cointegrating_share_omega,
                        3,
                    ),
                    unconditional_equity_premium_gamma0: Dual::<7>::variable(
                        spec.unconditional_equity_premium_gamma0,
                        4,
                    ),
                    cay_predictive_slope_beta_cay: Dual::<7>::variable(
                        spec.cay_predictive_slope_beta_cay,
                        5,
                    ),
                    conditional_market_volatility_sigma_m: Dual::<7>::variable(
                        spec.conditional_market_volatility_sigma_m,
                        6,
                    ),
                },
            );
            for (idx, &ad) in r_ad.conditional_market_sharpe_ratio.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.log_consumption_c += h;
                        sm.log_consumption_c -= h;
                    }
                    1 => {
                        sp.log_asset_wealth_a += h;
                        sm.log_asset_wealth_a -= h;
                    }
                    2 => {
                        sp.log_labor_income_y += h;
                        sm.log_labor_income_y -= h;
                    }
                    3 => {
                        sp.asset_wealth_cointegrating_share_omega += h;
                        sm.asset_wealth_cointegrating_share_omega -= h;
                    }
                    4 => {
                        sp.unconditional_equity_premium_gamma0 += h;
                        sm.unconditional_equity_premium_gamma0 -= h;
                    }
                    5 => {
                        sp.cay_predictive_slope_beta_cay += h;
                        sm.cay_predictive_slope_beta_cay -= h;
                    }
                    _ => {
                        sp.conditional_market_volatility_sigma_m += h;
                        sm.conditional_market_volatility_sigma_m -= h;
                    }
                }
                let fd = (evaluate_lettau_ludvigson_cay_consumption_wealth_forecast(&sp)
                    .conditional_market_sharpe_ratio
                    - evaluate_lettau_ludvigson_cay_consumption_wealth_forecast(&sm)
                        .conditional_market_sharpe_ratio)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_lettau_ludvigson_cay_consumption_wealth_forecast_module::*;

mod evaluate_lucas_mehra_prescott_endowment_tree_equity_premium_module {
    // Lucas (1978, Econometrica, "Asset Prices in an Exchange Economy") & Mehra & Prescott (1985,
    // J. Monetary Economics, "The Equity Premium: A Puzzle") Endowment Tree CRRA Pricing Kernel Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct LucasMehraPrescottEndowmentTreeSpec<T: Real = f64> {
        pub mean_consumption_growth_mu_c: T,
        pub consumption_growth_vol_sigma_c: T,
        pub dividend_consumption_leverage_lambda: T,
        pub crra_risk_aversion_gamma: T,
        pub subjective_discount_factor_beta: T,
        pub historical_observed_equity_premium: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct LucasMehraPrescottEndowmentTreeReport<T: Real = f64> {
        pub equilibrium_log_risk_free_rate: T,
        pub precautionary_savings_rate_drag: T,
        pub theoretical_crra_equity_risk_premium: T,
        pub mehra_prescott_equity_premium_puzzle_wedge: T,
        pub implied_puzzle_risk_aversion_gamma: T,
        pub lucas_tree_price_dividend_ratio: T,
    }

    pub fn evaluate_lucas_mehra_prescott_endowment_tree_equity_premium<T: Real>(
        spec: &LucasMehraPrescottEndowmentTreeSpec<T>,
    ) -> LucasMehraPrescottEndowmentTreeReport<T> {
        let mu_c = spec.mean_consumption_growth_mu_c;
        let sig_c = spec
            .consumption_growth_vol_sigma_c
            .max_val(T::from_f64(0.005));
        let lam_d = spec
            .dividend_consumption_leverage_lambda
            .max_val(T::from_f64(0.5));
        let gamma = spec.crra_risk_aversion_gamma.max_val(T::from_f64(0.5));
        let beta = spec
            .subjective_discount_factor_beta
            .max_val(T::from_f64(0.85))
            .min_val(T::from_f64(0.999));
        let erp_hist = spec.historical_observed_equity_premium;

        let rho_time = -beta.ln();
        let sig_c2 = sig_c * sig_c;
        let prec_savings = T::from_f64(0.5) * gamma * gamma * sig_c2;

        // Lucas (1978) / Mehra-Prescott (1985) log risk-free rate and CRRA equity premium:
        let r_f = rho_time + gamma * mu_c - prec_savings;
        let erp_crra = gamma * lam_d * sig_c2;
        let puzzle_wedge = erp_hist - erp_crra;
        let gamma_implied = erp_hist / (lam_d * sig_c2);

        let discount_minus_growth = (r_f + erp_crra - mu_c).max_val(T::from_f64(0.01));
        let pd_ratio = T::one() / discount_minus_growth;

        LucasMehraPrescottEndowmentTreeReport {
            equilibrium_log_risk_free_rate: r_f,
            precautionary_savings_rate_drag: prec_savings,
            theoretical_crra_equity_risk_premium: erp_crra,
            mehra_prescott_equity_premium_puzzle_wedge: puzzle_wedge,
            implied_puzzle_risk_aversion_gamma: gamma_implied,
            lucas_tree_price_dividend_ratio: pd_ratio,
        }
    }

    #[cfg(test)]
    mod evaluate_lucas_mehra_prescott_endowment_tree_equity_premium_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn lucas_mehra_prescott_endowment_tree_invariants_and_ad_parity() {
            let spec = LucasMehraPrescottEndowmentTreeSpec {
                mean_consumption_growth_mu_c: 0.018,
                consumption_growth_vol_sigma_c: 0.036,
                dividend_consumption_leverage_lambda: 3.0,
                crra_risk_aversion_gamma: 2.5,
                subjective_discount_factor_beta: 0.99,
                historical_observed_equity_premium: 0.0618,
            };
            let rep = evaluate_lucas_mehra_prescott_endowment_tree_equity_premium(&spec);
            assert!(rep.equilibrium_log_risk_free_rate > 0.0);
            assert!(rep.mehra_prescott_equity_premium_puzzle_wedge > 0.04);
            assert!(rep.implied_puzzle_risk_aversion_gamma > 10.0);
            assert!(rep.lucas_tree_price_dividend_ratio > 10.0);

            let h = 1e-5;
            let r_ad = evaluate_lucas_mehra_prescott_endowment_tree_equity_premium(
                &LucasMehraPrescottEndowmentTreeSpec {
                    mean_consumption_growth_mu_c: Dual::<6>::variable(
                        spec.mean_consumption_growth_mu_c,
                        0,
                    ),
                    consumption_growth_vol_sigma_c: Dual::<6>::variable(
                        spec.consumption_growth_vol_sigma_c,
                        1,
                    ),
                    dividend_consumption_leverage_lambda: Dual::<6>::variable(
                        spec.dividend_consumption_leverage_lambda,
                        2,
                    ),
                    crra_risk_aversion_gamma: Dual::<6>::variable(spec.crra_risk_aversion_gamma, 3),
                    subjective_discount_factor_beta: Dual::<6>::variable(
                        spec.subjective_discount_factor_beta,
                        4,
                    ),
                    historical_observed_equity_premium: Dual::<6>::variable(
                        spec.historical_observed_equity_premium,
                        5,
                    ),
                },
            );
            let comp_ad = r_ad.lucas_tree_price_dividend_ratio
                + r_ad.mehra_prescott_equity_premium_puzzle_wedge * Dual::<6>::from_value(100.0);
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.mean_consumption_growth_mu_c += h;
                        sm.mean_consumption_growth_mu_c -= h;
                    }
                    1 => {
                        sp.consumption_growth_vol_sigma_c += h;
                        sm.consumption_growth_vol_sigma_c -= h;
                    }
                    2 => {
                        sp.dividend_consumption_leverage_lambda += h;
                        sm.dividend_consumption_leverage_lambda -= h;
                    }
                    3 => {
                        sp.crra_risk_aversion_gamma += h;
                        sm.crra_risk_aversion_gamma -= h;
                    }
                    4 => {
                        sp.subjective_discount_factor_beta += h;
                        sm.subjective_discount_factor_beta -= h;
                    }
                    _ => {
                        sp.historical_observed_equity_premium += h;
                        sm.historical_observed_equity_premium -= h;
                    }
                }
                let rp = evaluate_lucas_mehra_prescott_endowment_tree_equity_premium(&sp);
                let rm = evaluate_lucas_mehra_prescott_endowment_tree_equity_premium(&sm);
                let fd = ((rp.lucas_tree_price_dividend_ratio
                    + rp.mehra_prescott_equity_premium_puzzle_wedge * 100.0)
                    - (rm.lucas_tree_price_dividend_ratio
                        + rm.mehra_prescott_equity_premium_puzzle_wedge * 100.0))
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_lucas_mehra_prescott_endowment_tree_equity_premium_module::*;

mod evaluate_campbell_shiller_cape_excess_yield_forecast_module {
    // Shiller (1981, AER), Campbell & Shiller (1988, RFS; 1998, J. Portfolio Management)
    // & Ilmanen (2011, Expected Returns, Ch. 8) Cyclically Adjusted Price-to-Earnings (CAPE),
    // Excess CAPE Yield (ECY) & Mean-Reverting Valuation Multiple Return Forecast Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct CampbellShillerCapeExcessYieldSpec<T: Real = f64> {
        pub real_equity_index_price_p: T,
        pub ten_year_real_average_earnings_e10: T,
        pub historical_mean_log_cape: T,
        pub annual_valuation_reversion_speed_phi: T,
        pub expected_real_earnings_growth_g: T,
        pub dividend_payout_ratio_d: T,
        pub ten_year_real_bond_yield_y10: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct CampbellShillerCapeExcessYieldReport<T: Real = f64> {
        pub cyclically_adjusted_pe_cape: T,
        pub cyclically_adjusted_earnings_yield_caey: T,
        pub shiller_excess_cape_yield_ecy: T,
        pub valuation_reversion_annual_drag: T,
        pub expected_annualized_real_equity_return: T,
        pub implied_ten_year_equity_risk_premium: T,
    }

    pub fn evaluate_campbell_shiller_cape_excess_yield_forecast<T: Real>(
        spec: &CampbellShillerCapeExcessYieldSpec<T>,
    ) -> CampbellShillerCapeExcessYieldReport<T> {
        let p = spec.real_equity_index_price_p.max_val(T::from_f64(1.0));
        let e10 = spec
            .ten_year_real_average_earnings_e10
            .max_val(T::from_f64(0.1));
        let mu_log_cape = spec.historical_mean_log_cape;
        let phi = spec
            .annual_valuation_reversion_speed_phi
            .max_val(T::from_f64(1e-3));
        let g = spec.expected_real_earnings_growth_g;
        let d_pay = spec
            .dividend_payout_ratio_d
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.95));
        let y10 = spec.ten_year_real_bond_yield_y10;

        let cape = p / e10;
        let caey = e10 / p;
        let ecy = caey - y10;

        let val_drag = -phi * (cape.ln() - mu_log_cape);
        let exp_real_ret = d_pay * caey + g + val_drag;
        let erp_10 = exp_real_ret - y10;

        CampbellShillerCapeExcessYieldReport {
            cyclically_adjusted_pe_cape: cape,
            cyclically_adjusted_earnings_yield_caey: caey,
            shiller_excess_cape_yield_ecy: ecy,
            valuation_reversion_annual_drag: val_drag,
            expected_annualized_real_equity_return: exp_real_ret,
            implied_ten_year_equity_risk_premium: erp_10,
        }
    }

    #[cfg(test)]
    mod evaluate_campbell_shiller_cape_excess_yield_forecast_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn campbell_shiller_cape_excess_yield_invariants_and_ad_parity() {
            let spec = CampbellShillerCapeExcessYieldSpec {
                real_equity_index_price_p: 4800.0,
                ten_year_real_average_earnings_e10: 160.0,
                historical_mean_log_cape: 3.05,
                annual_valuation_reversion_speed_phi: 0.065,
                expected_real_earnings_growth_g: 0.018,
                dividend_payout_ratio_d: 0.45,
                ten_year_real_bond_yield_y10: 0.019,
            };
            let rep = evaluate_campbell_shiller_cape_excess_yield_forecast(&spec);
            assert!((rep.cyclically_adjusted_pe_cape - 30.0).abs() < 1e-12);
            assert!(rep.shiller_excess_cape_yield_ecy > 0.01);
            assert!(rep.valuation_reversion_annual_drag < 0.0);

            let h = 1e-5;
            let r_ad = evaluate_campbell_shiller_cape_excess_yield_forecast(
                &CampbellShillerCapeExcessYieldSpec {
                    real_equity_index_price_p: Dual::<7>::variable(
                        spec.real_equity_index_price_p,
                        0,
                    ),
                    ten_year_real_average_earnings_e10: Dual::<7>::variable(
                        spec.ten_year_real_average_earnings_e10,
                        1,
                    ),
                    historical_mean_log_cape: Dual::<7>::variable(spec.historical_mean_log_cape, 2),
                    annual_valuation_reversion_speed_phi: Dual::<7>::variable(
                        spec.annual_valuation_reversion_speed_phi,
                        3,
                    ),
                    expected_real_earnings_growth_g: Dual::<7>::variable(
                        spec.expected_real_earnings_growth_g,
                        4,
                    ),
                    dividend_payout_ratio_d: Dual::<7>::variable(spec.dividend_payout_ratio_d, 5),
                    ten_year_real_bond_yield_y10: Dual::<7>::variable(
                        spec.ten_year_real_bond_yield_y10,
                        6,
                    ),
                },
            );
            for (idx, &ad) in r_ad
                .implied_ten_year_equity_risk_premium
                .d
                .iter()
                .enumerate()
            {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.real_equity_index_price_p += h;
                        sm.real_equity_index_price_p -= h;
                    }
                    1 => {
                        sp.ten_year_real_average_earnings_e10 += h;
                        sm.ten_year_real_average_earnings_e10 -= h;
                    }
                    2 => {
                        sp.historical_mean_log_cape += h;
                        sm.historical_mean_log_cape -= h;
                    }
                    3 => {
                        sp.annual_valuation_reversion_speed_phi += h;
                        sm.annual_valuation_reversion_speed_phi -= h;
                    }
                    4 => {
                        sp.expected_real_earnings_growth_g += h;
                        sm.expected_real_earnings_growth_g -= h;
                    }
                    5 => {
                        sp.dividend_payout_ratio_d += h;
                        sm.dividend_payout_ratio_d -= h;
                    }
                    _ => {
                        sp.ten_year_real_bond_yield_y10 += h;
                        sm.ten_year_real_bond_yield_y10 -= h;
                    }
                }
                let fd = (evaluate_campbell_shiller_cape_excess_yield_forecast(&sp)
                    .implied_ten_year_equity_risk_premium
                    - evaluate_campbell_shiller_cape_excess_yield_forecast(&sm)
                        .implied_ten_year_equity_risk_premium)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_campbell_shiller_cape_excess_yield_forecast_module::*;

mod evaluate_blanchard_lhuillier_lorenzoni_news_noise_kalman_filter_module {
    // Blanchard, L'Huillier & Lorenzoni (2013, American Economic Review 103(7), "News, Noise,
    // and Fluctuations: An Empirical Exploration"), Sims (2003, JME 50(3), "Implications of
    // Rational Inattention"), & Veldkamp (2011, Information Choice in Macroeconomics and Finance)
    // News-vs-Noise Kalman Signal Extraction & Rational-Inattention Output-Inflation Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BlanchardLhuillierLorenzoniNewsNoiseSpec<T: Real = f64> {
        pub fundamental_productivity_vol_sigma_a: T,
        pub exogenous_news_noise_vol_sigma_e: T,
        pub shannon_channel_capacity_bits_kappa: T,
        pub observed_macro_news_signal_s: T,
        pub consumption_demand_elasticity_gamma_c: T,
        pub phillips_curve_slope_kappa_pc: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BlanchardLhuillierLorenzoniNewsNoiseReport<T: Real = f64> {
        pub rational_inattention_endogenous_noise_var: T,
        pub optimal_kalman_signal_extraction_gain: T,
        pub posterior_fundamental_productivity_hat_a: T,
        pub posterior_fundamental_variance: T,
        pub sentiment_driven_output_gap_impulse: T,
        pub phillips_curve_inflation_impulse: T,
    }

    pub fn evaluate_blanchard_lhuillier_lorenzoni_news_noise_kalman_filter<T: Real>(
        spec: &BlanchardLhuillierLorenzoniNewsNoiseSpec<T>,
    ) -> BlanchardLhuillierLorenzoniNewsNoiseReport<T> {
        let sig_a = spec
            .fundamental_productivity_vol_sigma_a
            .max_val(T::from_f64(0.001));
        let sig_e = spec
            .exogenous_news_noise_vol_sigma_e
            .max_val(T::from_f64(0.001));
        let kappa_bits = spec
            .shannon_channel_capacity_bits_kappa
            .max_val(T::from_f64(0.10));
        let s_obs = spec.observed_macro_news_signal_s;
        let gamma_c = spec.consumption_demand_elasticity_gamma_c;
        let kappa_pc = spec.phillips_curve_slope_kappa_pc;

        let var_a = sig_a * sig_a;
        let var_e = sig_e * sig_e;

        // Sims (2003) Shannon channel capacity constraint I(a; s_ri) <= kappa_bits => var_ri = var_a / (2^{2*kappa_bits} - 1):
        let two_ln2 = T::from_f64(1.386_294_361_119_890_6);
        let cap_denom = ((two_ln2 * kappa_bits).exp() - T::one()).max_val(T::from_f64(1e-4));
        let var_ri = var_a / cap_denom;

        let var_tot = var_a + var_e + var_ri;
        let gain_g = var_a / var_tot;

        let a_hat = gain_g * s_obs;
        let post_var = (T::one() - gain_g) * var_a;
        let dy_impulse = gamma_c * a_hat;
        let dpi_impulse = kappa_pc * dy_impulse;

        BlanchardLhuillierLorenzoniNewsNoiseReport {
            rational_inattention_endogenous_noise_var: var_ri,
            optimal_kalman_signal_extraction_gain: gain_g,
            posterior_fundamental_productivity_hat_a: a_hat,
            posterior_fundamental_variance: post_var,
            sentiment_driven_output_gap_impulse: dy_impulse,
            phillips_curve_inflation_impulse: dpi_impulse,
        }
    }

    #[cfg(test)]
    mod evaluate_blanchard_lhuillier_lorenzoni_news_noise_kalman_filter_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn blanchard_lhuillier_lorenzoni_news_noise_invariants_and_ad_parity() {
            let spec = BlanchardLhuillierLorenzoniNewsNoiseSpec {
                fundamental_productivity_vol_sigma_a: 0.015,
                exogenous_news_noise_vol_sigma_e: 0.020,
                shannon_channel_capacity_bits_kappa: 1.50,
                observed_macro_news_signal_s: 0.025,
                consumption_demand_elasticity_gamma_c: 0.85,
                phillips_curve_slope_kappa_pc: 0.18,
            };
            let rep = evaluate_blanchard_lhuillier_lorenzoni_news_noise_kalman_filter(&spec);
            assert!(
                rep.optimal_kalman_signal_extraction_gain > 0.0
                    && rep.optimal_kalman_signal_extraction_gain < 1.0
            );
            assert!(rep.posterior_fundamental_variance < 0.015 * 0.015);

            let h = 1e-6;
            let r_ad = evaluate_blanchard_lhuillier_lorenzoni_news_noise_kalman_filter(
                &BlanchardLhuillierLorenzoniNewsNoiseSpec {
                    fundamental_productivity_vol_sigma_a: Dual::<6>::variable(
                        spec.fundamental_productivity_vol_sigma_a,
                        0,
                    ),
                    exogenous_news_noise_vol_sigma_e: Dual::<6>::variable(
                        spec.exogenous_news_noise_vol_sigma_e,
                        1,
                    ),
                    shannon_channel_capacity_bits_kappa: Dual::<6>::variable(
                        spec.shannon_channel_capacity_bits_kappa,
                        2,
                    ),
                    observed_macro_news_signal_s: Dual::<6>::variable(
                        spec.observed_macro_news_signal_s,
                        3,
                    ),
                    consumption_demand_elasticity_gamma_c: Dual::<6>::variable(
                        spec.consumption_demand_elasticity_gamma_c,
                        4,
                    ),
                    phillips_curve_slope_kappa_pc: Dual::<6>::variable(
                        spec.phillips_curve_slope_kappa_pc,
                        5,
                    ),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.fundamental_productivity_vol_sigma_a,
                    spec.exogenous_news_noise_vol_sigma_e,
                    spec.shannon_channel_capacity_bits_kappa,
                    spec.observed_macro_news_signal_s,
                    spec.consumption_demand_elasticity_gamma_c,
                    spec.phillips_curve_slope_kappa_pc,
                ];
                v[i] += h;
                let up = evaluate_blanchard_lhuillier_lorenzoni_news_noise_kalman_filter(
                    &BlanchardLhuillierLorenzoniNewsNoiseSpec {
                        fundamental_productivity_vol_sigma_a: v[0],
                        exogenous_news_noise_vol_sigma_e: v[1],
                        shannon_channel_capacity_bits_kappa: v[2],
                        observed_macro_news_signal_s: v[3],
                        consumption_demand_elasticity_gamma_c: v[4],
                        phillips_curve_slope_kappa_pc: v[5],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_blanchard_lhuillier_lorenzoni_news_noise_kalman_filter(
                    &BlanchardLhuillierLorenzoniNewsNoiseSpec {
                        fundamental_productivity_vol_sigma_a: v[0],
                        exogenous_news_noise_vol_sigma_e: v[1],
                        shannon_channel_capacity_bits_kappa: v[2],
                        observed_macro_news_signal_s: v[3],
                        consumption_demand_elasticity_gamma_c: v[4],
                        phillips_curve_slope_kappa_pc: v[5],
                    },
                );
                let fd = (up.phillips_curve_inflation_impulse
                    - dn.phillips_curve_inflation_impulse)
                    / (2.0 * h);
                assert!(
                    (r_ad.phillips_curve_inflation_impulse.d[i] - fd).abs() < 1e-6,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_blanchard_lhuillier_lorenzoni_news_noise_kalman_filter_module::*;

mod evaluate_bernanke_gertler_gilchrist_financial_accelerator_module {
    // Bernanke, Gertler & Gilchrist (1999, Handbook of Macroeconomics, Vol. 1C, Ch. 21,
    // "The Financial Accelerator in a Quantitative Business Cycle Framework"),
    // Townsend (1979, JET, "Costly State Verification"), & Kiyotaki & Moore (1997, JPE)
    // External Finance Premium, Net-Worth Leverage Elasticity & Financial Accelerator Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BernankeGertlerGilchristAcceleratorSpec<T: Real = f64> {
        pub capital_asset_value_qk: T,
        pub entrepreneur_net_worth_n: T,
        pub csv_auditing_cost_fraction_mu: T,
        pub efp_leverage_elasticity_nu: T,
        pub risk_free_real_rate_r: T,
        pub unexpected_asset_price_shock_dq: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BernankeGertlerGilchristAcceleratorReport<T: Real = f64> {
        pub entrepreneurial_leverage_ratio_l: T,
        pub external_finance_premium_spread: T,
        pub required_expected_return_on_capital: T,
        pub net_worth_shock_multiplier: T,
        pub bgg_investment_accelerator_multiplier: T,
        pub amplified_investment_impulse: T,
    }

    pub fn evaluate_bernanke_gertler_gilchrist_financial_accelerator<T: Real>(
        spec: &BernankeGertlerGilchristAcceleratorSpec<T>,
    ) -> BernankeGertlerGilchristAcceleratorReport<T> {
        let qk = spec.capital_asset_value_qk.max_val(T::from_f64(1.0));
        let net_w = spec.entrepreneur_net_worth_n.max_val(T::from_f64(0.5));
        let mu_csv = spec
            .csv_auditing_cost_fraction_mu
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.50));
        let nu = spec.efp_leverage_elasticity_nu.max_val(T::from_f64(0.01));
        let r_f = spec.risk_free_real_rate_r;
        let dq = spec.unexpected_asset_price_shock_dq;

        let lev = (qk / net_w).max_val(T::from_f64(1.01));

        // BGG (1999) external finance premium s(L) = mu_csv * psi0 * (L^nu - 1) + base_agency:
        let psi0 = T::from_f64(0.15);
        let l_pow_nu = (nu * lev.ln()).exp();
        let efp = mu_csv * psi0 * l_pow_nu;
        let r_k = (T::one() + r_f) * (T::one() + efp) - T::one();

        // Leveraged net-worth elasticity dn / dq = L:
        let nw_mult = lev * dq;
        // BGG Financial Accelerator multiplier on Tobin's q investment impulse:
        let bgg_mult = T::one() + nu * (lev - T::one()) * (T::one() + efp);
        let inv_impulse = bgg_mult * dq;

        BernankeGertlerGilchristAcceleratorReport {
            entrepreneurial_leverage_ratio_l: lev,
            external_finance_premium_spread: efp,
            required_expected_return_on_capital: r_k,
            net_worth_shock_multiplier: nw_mult,
            bgg_investment_accelerator_multiplier: bgg_mult,
            amplified_investment_impulse: inv_impulse,
        }
    }

    #[cfg(test)]
    mod evaluate_bernanke_gertler_gilchrist_financial_accelerator_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bgg_financial_accelerator_invariants_and_ad_parity() {
            let spec = BernankeGertlerGilchristAcceleratorSpec {
                capital_asset_value_qk: 200.0,
                entrepreneur_net_worth_n: 80.0,
                csv_auditing_cost_fraction_mu: 0.12,
                efp_leverage_elasticity_nu: 0.08,
                risk_free_real_rate_r: 0.02,
                unexpected_asset_price_shock_dq: 0.03,
            };
            let rep = evaluate_bernanke_gertler_gilchrist_financial_accelerator(&spec);
            assert!((rep.entrepreneurial_leverage_ratio_l - 2.5).abs() < 1e-12);
            assert!(rep.external_finance_premium_spread > 0.0);
            assert!(rep.bgg_investment_accelerator_multiplier > 1.0);

            let h = 1e-6;
            let r_ad = evaluate_bernanke_gertler_gilchrist_financial_accelerator(
                &BernankeGertlerGilchristAcceleratorSpec {
                    capital_asset_value_qk: Dual::<6>::variable(spec.capital_asset_value_qk, 0),
                    entrepreneur_net_worth_n: Dual::<6>::variable(spec.entrepreneur_net_worth_n, 1),
                    csv_auditing_cost_fraction_mu: Dual::<6>::variable(
                        spec.csv_auditing_cost_fraction_mu,
                        2,
                    ),
                    efp_leverage_elasticity_nu: Dual::<6>::variable(
                        spec.efp_leverage_elasticity_nu,
                        3,
                    ),
                    risk_free_real_rate_r: Dual::<6>::variable(spec.risk_free_real_rate_r, 4),
                    unexpected_asset_price_shock_dq: Dual::<6>::variable(
                        spec.unexpected_asset_price_shock_dq,
                        5,
                    ),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.capital_asset_value_qk,
                    spec.entrepreneur_net_worth_n,
                    spec.csv_auditing_cost_fraction_mu,
                    spec.efp_leverage_elasticity_nu,
                    spec.risk_free_real_rate_r,
                    spec.unexpected_asset_price_shock_dq,
                ];
                v[i] += h;
                let up = evaluate_bernanke_gertler_gilchrist_financial_accelerator(
                    &BernankeGertlerGilchristAcceleratorSpec {
                        capital_asset_value_qk: v[0],
                        entrepreneur_net_worth_n: v[1],
                        csv_auditing_cost_fraction_mu: v[2],
                        efp_leverage_elasticity_nu: v[3],
                        risk_free_real_rate_r: v[4],
                        unexpected_asset_price_shock_dq: v[5],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_bernanke_gertler_gilchrist_financial_accelerator(
                    &BernankeGertlerGilchristAcceleratorSpec {
                        capital_asset_value_qk: v[0],
                        entrepreneur_net_worth_n: v[1],
                        csv_auditing_cost_fraction_mu: v[2],
                        efp_leverage_elasticity_nu: v[3],
                        risk_free_real_rate_r: v[4],
                        unexpected_asset_price_shock_dq: v[5],
                    },
                );
                let val_ad =
                    r_ad.required_expected_return_on_capital + r_ad.amplified_investment_impulse;
                let fd = ((up.required_expected_return_on_capital
                    + up.amplified_investment_impulse)
                    - (dn.required_expected_return_on_capital + dn.amplified_investment_impulse))
                    / (2.0 * h);
                assert!((val_ad.d[i] - fd).abs() < 1e-5, "input {i}");
            }
        }
    }
}
pub use evaluate_bernanke_gertler_gilchrist_financial_accelerator_module::*;

mod evaluate_kiyotaki_moore_collateral_credit_cycle_multiplier_module {
    // Kiyotaki & Moore (1997, Journal of Political Economy 105(2), "Credit Cycles",
    // Section II, Eqs. 6-12) & Brunnermeier & Sannikov (2014, AER 104(2))
    // Collateral Borrowing Constraint, User Cost Down-Payment & Credit Cycle Multiplier Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct KiyotakiMooreCollateralCycleSpec<T: Real = f64> {
        pub land_asset_price_q: T,
        pub expected_next_land_price_q_next: T,
        pub gross_interest_rate_r_gross: T,
        pub farmer_productivity_a: T,
        pub existing_land_holding_k_prev: T,
        pub inherited_debt_b_prev: T,
        pub intertemporal_price_elasticity_eta: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct KiyotakiMooreCollateralCycleReport<T: Real = f64> {
        pub required_down_payment_user_cost: T,
        pub constrained_farmer_net_worth: T,
        pub collateral_constrained_land_demand: T,
        pub maximum_pledgeable_debt_capacity: T,
        pub static_net_worth_leverage_multiplier: T,
        pub dynamic_credit_cycle_amplification_multiplier: T,
    }

    pub fn evaluate_kiyotaki_moore_collateral_credit_cycle_multiplier<T: Real>(
        spec: &KiyotakiMooreCollateralCycleSpec<T>,
    ) -> KiyotakiMooreCollateralCycleReport<T> {
        let q_t = spec.land_asset_price_q.max_val(T::from_f64(1.0));
        let q_next = spec
            .expected_next_land_price_q_next
            .max_val(T::from_f64(0.5));
        let r_gross = spec.gross_interest_rate_r_gross.max_val(T::from_f64(1.005));
        let a_prod = spec.farmer_productivity_a.max_val(T::from_f64(0.01));
        let k_prev = spec.existing_land_holding_k_prev.max_val(T::from_f64(0.1));
        let b_prev = spec.inherited_debt_b_prev.max_val(T::zero());
        let eta = spec
            .intertemporal_price_elasticity_eta
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.50));

        // Kiyotaki & Moore (1997, Section II, Eqs. 6-10):
        // Down-payment user cost per unit of land u_t = q_t - q_{t+1} / R:
        let pv_q_next = q_next / r_gross;
        let user_cost_u = (q_t - pv_q_next).max_val(T::from_f64(0.05));

        // Constrained farmer net worth N_t = (a + q_t) * k_{t-1} - R * b_{t-1}:
        let net_worth = ((a_prod + q_t) * k_prev - r_gross * b_prev).max_val(T::from_f64(0.1));
        let k_demand = net_worth / user_cost_u;
        let max_debt = pv_q_next * k_demand;

        // Static leveraged net-worth multiplier (d ln k_t / d ln q_t) and dynamic credit-cycle multiplier:
        let mult_static = (q_t * k_prev) / net_worth;
        let raw_fb = eta * (pv_q_next / user_cost_u);
        let feedback = (T::from_f64(0.85) * raw_fb) / (T::from_f64(0.85) + raw_fb);
        let mult_dynamic = mult_static / (T::one() - feedback);

        KiyotakiMooreCollateralCycleReport {
            required_down_payment_user_cost: user_cost_u,
            constrained_farmer_net_worth: net_worth,
            collateral_constrained_land_demand: k_demand,
            maximum_pledgeable_debt_capacity: max_debt,
            static_net_worth_leverage_multiplier: mult_static,
            dynamic_credit_cycle_amplification_multiplier: mult_dynamic,
        }
    }

    #[cfg(test)]
    mod evaluate_kiyotaki_moore_collateral_credit_cycle_multiplier_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn kiyotaki_moore_collateral_cycle_invariants_and_ad_parity() {
            let spec = KiyotakiMooreCollateralCycleSpec {
                land_asset_price_q: 100.0,
                expected_next_land_price_q_next: 101.0,
                gross_interest_rate_r_gross: 1.05,
                farmer_productivity_a: 6.0,
                existing_land_holding_k_prev: 10.0,
                inherited_debt_b_prev: 880.0,
                intertemporal_price_elasticity_eta: 0.02,
            };
            let rep = evaluate_kiyotaki_moore_collateral_credit_cycle_multiplier(&spec);
            assert!(rep.required_down_payment_user_cost > 0.0);
            assert!(
                rep.dynamic_credit_cycle_amplification_multiplier
                    > rep.static_net_worth_leverage_multiplier
            );

            let h = 1e-6;
            let r_ad = evaluate_kiyotaki_moore_collateral_credit_cycle_multiplier(
                &KiyotakiMooreCollateralCycleSpec {
                    land_asset_price_q: Dual::<7>::variable(spec.land_asset_price_q, 0),
                    expected_next_land_price_q_next: Dual::<7>::variable(
                        spec.expected_next_land_price_q_next,
                        1,
                    ),
                    gross_interest_rate_r_gross: Dual::<7>::variable(
                        spec.gross_interest_rate_r_gross,
                        2,
                    ),
                    farmer_productivity_a: Dual::<7>::variable(spec.farmer_productivity_a, 3),
                    existing_land_holding_k_prev: Dual::<7>::variable(
                        spec.existing_land_holding_k_prev,
                        4,
                    ),
                    inherited_debt_b_prev: Dual::<7>::variable(spec.inherited_debt_b_prev, 5),
                    intertemporal_price_elasticity_eta: Dual::<7>::variable(
                        spec.intertemporal_price_elasticity_eta,
                        6,
                    ),
                },
            );
            for i in 0..7 {
                let mut v = [
                    spec.land_asset_price_q,
                    spec.expected_next_land_price_q_next,
                    spec.gross_interest_rate_r_gross,
                    spec.farmer_productivity_a,
                    spec.existing_land_holding_k_prev,
                    spec.inherited_debt_b_prev,
                    spec.intertemporal_price_elasticity_eta,
                ];
                v[i] += h;
                let up = evaluate_kiyotaki_moore_collateral_credit_cycle_multiplier(
                    &KiyotakiMooreCollateralCycleSpec {
                        land_asset_price_q: v[0],
                        expected_next_land_price_q_next: v[1],
                        gross_interest_rate_r_gross: v[2],
                        farmer_productivity_a: v[3],
                        existing_land_holding_k_prev: v[4],
                        inherited_debt_b_prev: v[5],
                        intertemporal_price_elasticity_eta: v[6],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_kiyotaki_moore_collateral_credit_cycle_multiplier(
                    &KiyotakiMooreCollateralCycleSpec {
                        land_asset_price_q: v[0],
                        expected_next_land_price_q_next: v[1],
                        gross_interest_rate_r_gross: v[2],
                        farmer_productivity_a: v[3],
                        existing_land_holding_k_prev: v[4],
                        inherited_debt_b_prev: v[5],
                        intertemporal_price_elasticity_eta: v[6],
                    },
                );
                let fd = (up.dynamic_credit_cycle_amplification_multiplier
                    - dn.dynamic_credit_cycle_amplification_multiplier)
                    / (2.0 * h);
                assert!(
                    r_ad.dynamic_credit_cycle_amplification_multiplier.d[i].abs() > 1e-6,
                    "input {i} zero derivative"
                );
                assert!(
                    (r_ad.dynamic_credit_cycle_amplification_multiplier.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_kiyotaki_moore_collateral_credit_cycle_multiplier_module::*;

use ql_ad::Real;

/// Specification for the Gabaix (2012) Variable Rare Disasters & Linearity-Generating (LG)
/// stock resilience macro-asset-pricing engine.
///
/// Reference: Gabaix (2012, *Quarterly Journal of Economics* 127(2), "Variable Rare Disasters:
/// An Exactly Solved Framework for Ten Puzzles in Macro-Finance"), Barro (2006, *QJE* 121(3)), Rietz (1988).
#[derive(Debug, Clone, Copy)]
pub struct GabaixVariableRareDisasterSpec<T: Real> {
    /// CRRA coefficient of relative risk aversion $\gamma$.
    pub risk_aversion_gamma: T,
    /// Pure rate of time preference $\rho$.
    pub time_preference_rho: T,
    /// Normal-times consumption and dividend growth rate $g_c$.
    pub normal_consumption_growth_gc: T,
    /// Normal-times diffusive consumption volatility $\sigma_c$.
    pub consumption_volatility_sigma_c: T,
    /// Current time-varying rare disaster arrival probability $p_t$.
    pub current_disaster_probability_pt: T,
    /// Long-run baseline disaster arrival probability $p_*$.
    pub baseline_disaster_probability_p_star: T,
    /// Consumption recovery fraction in a disaster $B_c \in (0, 1)$.
    pub consumption_disaster_recovery_bc: T,
    /// Dividend recovery fraction in a disaster $F_d \in (0, 1)$.
    pub dividend_disaster_recovery_fd: T,
    /// Mean-reversion speed of stock resilience $\phi_H$.
    pub resilience_mean_reversion_phi_h: T,
}

/// Analytics returned by [`evaluate_gabaix_variable_rare_disaster_resilience_equity_puzzle`].
#[derive(Debug, Clone, Copy)]
pub struct GabaixVariableRareDisasterReport<T: Real> {
    /// Current stock resilience $H_t = p_t (B_c^{-\gamma} F_d - 1)$ (negative in levered equity).
    pub stock_resilience_ht: T,
    /// Exact Linearity-Generating (LG) equilibrium price-dividend ratio $P_t / D_t$.
    pub equilibrium_price_dividend_ratio: T,
    /// Barro–Gabaix real risk-free interest rate $r_{f,t}$.
    pub real_risk_free_rate: T,
    /// Conditional unlevered + disaster equity risk premium $ERP_t$.
    pub conditional_equity_risk_premium: T,
}

/// Evaluates the Gabaix (2012) Variable Rare Disasters model via exact Linearity-Generating (LG)
/// closed-form valuation of the price-dividend ratio, real risk-free rate, and disaster equity premium.
pub fn evaluate_gabaix_variable_rare_disaster_resilience_equity_puzzle<T: Real>(
    spec: &GabaixVariableRareDisasterSpec<T>,
) -> GabaixVariableRareDisasterReport<T> {
    let gamma = spec.risk_aversion_gamma;
    let bc = spec
        .consumption_disaster_recovery_bc
        .max_val(T::from_f64(0.1));
    let fd = spec
        .dividend_disaster_recovery_fd
        .max_val(T::from_f64(0.05));
    let sig_c = spec.consumption_volatility_sigma_c;
    let sig_c2 = sig_c * sig_c;

    // Marginal utility jump multiplier B_c^{-gamma} = exp(-gamma * ln(B_c))
    let bc_neg_gamma = (-gamma * bc.ln()).exp();

    // Stock resilience H_t = p_t * (B_c^{-gamma} * F_d - 1) and steady-state H_* = p_* * (B_c^{-gamma} * F_d - 1)
    let resilience_kern = bc_neg_gamma * fd - T::one();
    let h_t = spec.current_disaster_probability_pt * resilience_kern;
    let h_star = spec.baseline_disaster_probability_p_star * resilience_kern;
    let h_hat = h_t - h_star;

    // Effective steady-state equity discount spread delta_D = rho + (gamma - 1)*g_c - 0.5*gamma*(gamma-1)*sigma_c^2 - H_*
    let delta_d = (spec.time_preference_rho
        + (gamma - T::one()) * spec.normal_consumption_growth_gc
        - T::from_f64(0.5) * gamma * (gamma - T::one()) * sig_c2
        - h_star)
        .max_val(T::from_f64(1e-4));

    // Gabaix (2012, Proposition 1, Eq. 19): exact LG price-dividend ratio
    // P_t / D_t = (1 / delta_D) * (1 + H_hat_t / (delta_D + phi_H))
    let denom_h = (delta_d + spec.resilience_mean_reversion_phi_h).max_val(T::from_f64(1e-4));
    let pd_ratio = (T::one() / delta_d) * (T::one() + h_hat / denom_h);

    // Barro (2006) / Gabaix (2012, Eq. 14): real risk-free rate
    let rf = spec.time_preference_rho + gamma * spec.normal_consumption_growth_gc
        - T::from_f64(0.5) * gamma * (gamma + T::one()) * sig_c2
        - spec.current_disaster_probability_pt * (bc_neg_gamma - T::one());

    // Gabaix (2012, Eq. 22): conditional equity risk premium (diffusive CCAPM + rare disaster jump premium)
    let erp = gamma * sig_c2
        + spec.current_disaster_probability_pt * (bc_neg_gamma - T::one()) * (T::one() - fd);

    GabaixVariableRareDisasterReport {
        stock_resilience_ht: h_t,
        equilibrium_price_dividend_ratio: pd_ratio,
        real_risk_free_rate: rf,
        conditional_equity_risk_premium: erp,
    }
}

#[cfg(test)]
mod gabaix_rare_disaster_suite64_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn gabaix_variable_rare_disaster_invariants_and_ad_parity() {
        let base = GabaixVariableRareDisasterSpec {
            risk_aversion_gamma: 4.0_f64,
            time_preference_rho: 0.035,
            normal_consumption_growth_gc: 0.020,
            consumption_volatility_sigma_c: 0.018,
            current_disaster_probability_pt: 0.022,
            baseline_disaster_probability_p_star: 0.017,
            consumption_disaster_recovery_bc: 0.72,
            dividend_disaster_recovery_fd: 0.55,
            resilience_mean_reversion_phi_h: 0.18,
        };
        let r0 = evaluate_gabaix_variable_rare_disaster_resilience_equity_puzzle(&base);
        assert!(
            r0.equilibrium_price_dividend_ratio > 8.0 && r0.equilibrium_price_dividend_ratio < 40.0
        );
        assert!(r0.real_risk_free_rate > -0.02 && r0.real_risk_free_rate < 0.08);
        assert!(r0.conditional_equity_risk_premium > 0.02);

        let inputs = [
            base.risk_aversion_gamma,
            base.time_preference_rho,
            base.normal_consumption_growth_gc,
            base.consumption_volatility_sigma_c,
            base.current_disaster_probability_pt,
            base.baseline_disaster_probability_p_star,
            base.consumption_disaster_recovery_bc,
            base.dividend_disaster_recovery_fd,
            base.resilience_mean_reversion_phi_h,
        ];
        let mut dual_in = [Dual::<9>::from_f64(0.0); 9];
        for i in 0..9 {
            dual_in[i] = Dual::variable(inputs[i], i);
        }
        let spec_ad = GabaixVariableRareDisasterSpec {
            risk_aversion_gamma: dual_in[0],
            time_preference_rho: dual_in[1],
            normal_consumption_growth_gc: dual_in[2],
            consumption_volatility_sigma_c: dual_in[3],
            current_disaster_probability_pt: dual_in[4],
            baseline_disaster_probability_p_star: dual_in[5],
            consumption_disaster_recovery_bc: dual_in[6],
            dividend_disaster_recovery_fd: dual_in[7],
            resilience_mean_reversion_phi_h: dual_in[8],
        };
        let r_ad = evaluate_gabaix_variable_rare_disaster_resilience_equity_puzzle(&spec_ad);
        let eps = 1e-6;
        for i in 0..9 {
            let mut p = inputs;
            let mut m = inputs;
            p[i] += eps;
            m[i] -= eps;
            let mk = |v: [f64; 9]| GabaixVariableRareDisasterSpec {
                risk_aversion_gamma: v[0],
                time_preference_rho: v[1],
                normal_consumption_growth_gc: v[2],
                consumption_volatility_sigma_c: v[3],
                current_disaster_probability_pt: v[4],
                baseline_disaster_probability_p_star: v[5],
                consumption_disaster_recovery_bc: v[6],
                dividend_disaster_recovery_fd: v[7],
                resilience_mean_reversion_phi_h: v[8],
            };
            let rp = evaluate_gabaix_variable_rare_disaster_resilience_equity_puzzle(&mk(p));
            let rm = evaluate_gabaix_variable_rare_disaster_resilience_equity_puzzle(&mk(m));
            let fd = (rp.equilibrium_price_dividend_ratio - rm.equilibrium_price_dividend_ratio)
                / (2.0 * eps);
            assert!(
                r_ad.equilibrium_price_dividend_ratio.d[i].abs() > 1e-6,
                "zero derivative at {i}"
            );
            assert!((r_ad.equilibrium_price_dividend_ratio.d[i] - fd).abs() < 1e-4);
        }
    }
}

pub mod bansal_yaron_lrr_suite65_mod {
    use ql_ad::Real;

    /// Specification for the Bansal–Yaron (2004) Long-Run Risks (LRR) Epstein–Zin stochastic
    /// volatility macro-asset-pricing engine.
    ///
    /// Reference: Bansal & Yaron (2004, *Journal of Finance* 59(4), "Risks for the Long Run:
    /// A Potential Resolution of Asset Pricing Puzzles"), Epstein & Zin (1989, *Econometrica* 57(4)).
    #[derive(Debug, Clone, Copy)]
    pub struct BansalYaronLongRunRiskSpec<T: Real> {
        /// Relative risk aversion $\gamma > 1$.
        pub risk_aversion_gamma: T,
        /// Elasticity of intertemporal substitution $\psi > 1$.
        pub eis_psi: T,
        /// Persistence of expected consumption growth $\rho_x \in (0, 1)$.
        pub growth_persistence_rho_x: T,
        /// Long-run growth shock volatility multiplier $\varphi_e$.
        pub growth_shock_scale_phi_e: T,
        /// Persistence of economic volatility $\nu_\sigma \in (0, 1)$.
        pub volatility_persistence_nu: T,
        /// Volatility of economic variance $\sigma_w$.
        pub vol_of_vol_sigma_w: T,
        /// Current conditional consumption volatility $\sigma_t$.
        pub current_consumption_vol_sigma_t: T,
        /// Dividend leverage loading $\phi_d > 1$ on persistent growth $x_t$.
        pub dividend_leverage_phi_d: T,
        /// Campbell–Shiller log-linearization constant $\kappa_1 \in (0, 1)$.
        pub campbell_shiller_kappa1: T,
    }

    /// Analytics returned by [`evaluate_bansal_yaron_long_run_risk_stochastic_vol_equity_premium`].
    #[derive(Debug, Clone, Copy)]
    pub struct BansalYaronLongRunRiskReport<T: Real> {
        /// Wealth-consumption ratio loading $A_1 = \frac{1 - 1/\psi}{1 - \kappa_1 \rho_x}$ on long-run growth $x_t$.
        pub wealth_growth_loading_a1: T,
        /// Wealth-consumption ratio loading $A_2 < 0$ on economic variance $\sigma_t^2$.
        pub wealth_volatility_loading_a2: T,
        /// Market price of long-run growth risk $\lambda_e = (1 - \theta) \kappa_1 A_1 \varphi_e$.
        pub long_run_growth_price_of_risk_lambda_e: T,
        /// Total Bansal–Yaron conditional equity risk premium (short-run + long-run growth + vol risk).
        pub conditional_equity_risk_premium: T,
    }

    /// Evaluates the Bansal–Yaron (2004, Section I.B, Equations A5–A19) analytical Epstein–Zin
    /// wealth and dividend loadings, innovation risk prices, and conditional equity risk premium.
    pub fn evaluate_bansal_yaron_long_run_risk_stochastic_vol_equity_premium<T: Real>(
        spec: &BansalYaronLongRunRiskSpec<T>,
    ) -> BansalYaronLongRunRiskReport<T> {
        let gamma = spec.risk_aversion_gamma;
        let inv_psi = T::one() / spec.eis_psi.max_val(T::from_f64(0.1));
        let theta = (T::one() - gamma) / (T::one() - inv_psi);
        let k1 = spec.campbell_shiller_kappa1;

        // Wealth-consumption loadings A_1 and A_2 (Bansal & Yaron 2004, Eqs. A5 & A7)
        let denom_x = (T::one() - k1 * spec.growth_persistence_rho_x).max_val(T::from_f64(1e-4));
        let a1 = (T::one() - inv_psi) / denom_x;

        let denom_v = (T::one() - k1 * spec.volatility_persistence_nu).max_val(T::from_f64(1e-4));
        let term_a1 = a1 * k1 * spec.growth_shock_scale_phi_e;
        let a2 = (T::from_f64(0.5)
            * theta
            * ((T::one() - inv_psi) * (T::one() - inv_psi) + term_a1 * term_a1))
            / denom_v;

        // Dividend claim loadings A_{1,m} and A_{2,m} (Bansal & Yaron 2004, Eqs. A10-A12)
        let a1_m = (spec.dividend_leverage_phi_d - inv_psi) / denom_x;
        let beta_m_e = k1 * a1_m * spec.growth_shock_scale_phi_e;
        let a2_m = ((T::one() - theta) * k1 * a2 + T::from_f64(0.5) * theta * beta_m_e * beta_m_e)
            / denom_v;
        let beta_m_w = k1 * a2_m;

        // Market prices of short-run (lambda_eta), long-run growth (lambda_e), and volatility (lambda_w) risks
        let lambda_eta = gamma;
        let lambda_e = (T::one() - theta) * k1 * a1 * spec.growth_shock_scale_phi_e;
        let lambda_w = (T::one() - theta) * k1 * a2;

        let sig2 = spec.current_consumption_vol_sigma_t * spec.current_consumption_vol_sigma_t;
        let sig_w2 = spec.vol_of_vol_sigma_w * spec.vol_of_vol_sigma_w;

        // Conditional equity risk premium (Bansal & Yaron 2004, Eq. 17)
        let erp = (lambda_eta + beta_m_e * lambda_e) * sig2 + beta_m_w * lambda_w * sig_w2;

        BansalYaronLongRunRiskReport {
            wealth_growth_loading_a1: a1,
            wealth_volatility_loading_a2: a2,
            long_run_growth_price_of_risk_lambda_e: lambda_e,
            conditional_equity_risk_premium: erp,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bansal_yaron_lrr_invariants_and_ad_parity() {
            let base = BansalYaronLongRunRiskSpec {
                risk_aversion_gamma: 7.5_f64,
                eis_psi: 1.5,
                growth_persistence_rho_x: 0.975,
                growth_shock_scale_phi_e: 0.044,
                volatility_persistence_nu: 0.985,
                vol_of_vol_sigma_w: 0.00025,
                current_consumption_vol_sigma_t: 0.027,
                dividend_leverage_phi_d: 3.0,
                campbell_shiller_kappa1: 0.99,
            };
            let r0 = evaluate_bansal_yaron_long_run_risk_stochastic_vol_equity_premium(&base);
            assert!(r0.wealth_growth_loading_a1 > 0.0);
            assert!(r0.wealth_volatility_loading_a2 < 0.0);
            assert!(r0.long_run_growth_price_of_risk_lambda_e > 0.0);
            assert!(r0.conditional_equity_risk_premium > 0.01);

            let inputs = [
                base.risk_aversion_gamma,
                base.eis_psi,
                base.growth_persistence_rho_x,
                base.growth_shock_scale_phi_e,
                base.volatility_persistence_nu,
                base.vol_of_vol_sigma_w,
                base.current_consumption_vol_sigma_t,
                base.dividend_leverage_phi_d,
                base.campbell_shiller_kappa1,
            ];
            let mut dual_in = [Dual::<9>::from_f64(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BansalYaronLongRunRiskSpec {
                risk_aversion_gamma: dual_in[0],
                eis_psi: dual_in[1],
                growth_persistence_rho_x: dual_in[2],
                growth_shock_scale_phi_e: dual_in[3],
                volatility_persistence_nu: dual_in[4],
                vol_of_vol_sigma_w: dual_in[5],
                current_consumption_vol_sigma_t: dual_in[6],
                dividend_leverage_phi_d: dual_in[7],
                campbell_shiller_kappa1: dual_in[8],
            };
            let r_ad = evaluate_bansal_yaron_long_run_risk_stochastic_vol_equity_premium(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| BansalYaronLongRunRiskSpec {
                    risk_aversion_gamma: v[0],
                    eis_psi: v[1],
                    growth_persistence_rho_x: v[2],
                    growth_shock_scale_phi_e: v[3],
                    volatility_persistence_nu: v[4],
                    vol_of_vol_sigma_w: v[5],
                    current_consumption_vol_sigma_t: v[6],
                    dividend_leverage_phi_d: v[7],
                    campbell_shiller_kappa1: v[8],
                };
                let rp = evaluate_bansal_yaron_long_run_risk_stochastic_vol_equity_premium(&mk(p));
                let rm = evaluate_bansal_yaron_long_run_risk_stochastic_vol_equity_premium(&mk(m));
                let fd = (rp.conditional_equity_risk_premium - rm.conditional_equity_risk_premium)
                    / (2.0 * eps);
                assert!(
                    r_ad.conditional_equity_risk_premium.d[i].abs() > 1e-7,
                    "zero derivative at {i}"
                );
                assert!((r_ad.conditional_equity_risk_premium.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use bansal_yaron_lrr_suite65_mod::*;

pub mod he_krishnamurthy_intermediary_sdf_mod {
    use ql_ad::Real;

    /// Specification for the He & Krishnamurthy (2013, AER) / Adrian, Etula & Muir (2014, JF)
    /// intermediary capital constraint Stochastic Discount Factor (SDF) and risk premium model.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HeKrishnamurthyIntermediarySdfSpec<T: Real> {
        /// Intermediary specialist wealth share $w_I \in (0, 1)$ of aggregate wealth.
        pub intermediary_wealth_share_wi: T,
        /// Capital constraint binding threshold $w^* \in (0, 1)$.
        pub capital_constraint_threshold_w_star: T,
        /// Household equity intermediation pledgeability fraction $m_{\text{pledge}} \in (0, 1)$.
        pub equity_intermediation_pledgeability_m: T,
        /// Specialist relative risk aversion $\gamma_{\text{spec}} > 0$.
        pub specialist_risk_aversion_gamma: T,
        /// Fundamental dividend/output growth rate $g_A$.
        pub fundamental_growth_rate_ga: T,
        /// Fundamental risky asset return volatility $\sigma_A > 0$.
        pub asset_return_volatility_sigma_a: T,
        /// Subjective time discount rate $\rho_{\text{time}} > 0$.
        pub subjective_discount_rate_rho: T,
    }

    /// Analytical outputs for the He-Krishnamurthy intermediary asset pricing engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HeKrishnamurthyIntermediarySdfAnalytics<T: Real> {
        /// Effective intermediary risk-bearing leverage multiplier $\beta_{\text{lev}}(w_I) \ge 1$.
        pub effective_intermediary_leverage_multiplier: T,
        /// Intermediary SDF market price of risk (Sharpe ratio) $\Lambda_I(w_I) = \gamma_{\text{spec}} \sigma_A \beta_{\text{lev}}$.
        pub intermediary_sdf_market_price_of_risk: T,
        /// Equilibrium risk-free interest rate $r_f(w_I)$ including intermediary precautionary savings.
        pub equilibrium_risk_free_rate: T,
        /// Conditional intermediary-constrained equity risk premium $\mathbb{E}_t[r_A - r_f]$.
        pub conditional_intermediary_equity_premium: T,
    }

    /// Evaluates the He & Krishnamurthy (2013) intermediary capital constraint SDF, crisis
    /// Sharpe ratio amplification, precautionary risk-free rate, and equity risk premium.
    pub fn evaluate_he_krishnamurthy_intermediary_capital_constraint_sdf<T: Real>(
        spec: &HeKrishnamurthyIntermediarySdfSpec<T>,
    ) -> HeKrishnamurthyIntermediarySdfAnalytics<T> {
        let wi = spec
            .intermediary_wealth_share_wi
            .max_val(T::from_f64(0.02))
            .min_val(T::from_f64(0.98));
        let w_star = spec
            .capital_constraint_threshold_w_star
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.95));
        let m_pledge = spec
            .equity_intermediation_pledgeability_m
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.95));
        let gamma = spec
            .specialist_risk_aversion_gamma
            .max_val(T::from_f64(0.2));
        let sig_a = spec
            .asset_return_volatility_sigma_a
            .max_val(T::from_f64(1e-4));

        // Smooth capital shortfall below w*: when wi < w*, intermediary equity constraint binds
        // and intermediary risk exposure scales as 1 / (wi + m_pledge * (w_star - wi)_+)
        let eps = T::from_f64(0.01);
        let diff = w_star - wi;
        let shortfall = T::from_f64(0.5) * (diff + (diff * diff + eps * eps).sqrt());
        let effective_cap_share = (wi + m_pledge * shortfall) / w_star;
        let beta_lev = T::one() / effective_cap_share;

        let lambda_sdf = gamma * sig_a * beta_lev;
        let eq_prem = lambda_sdf * sig_a;
        let rf = spec.subjective_discount_rate_rho + gamma * spec.fundamental_growth_rate_ga
            - T::from_f64(0.5) * gamma * (gamma + T::one()) * sig_a * sig_a * beta_lev * beta_lev;

        HeKrishnamurthyIntermediarySdfAnalytics {
            effective_intermediary_leverage_multiplier: beta_lev,
            intermediary_sdf_market_price_of_risk: lambda_sdf,
            equilibrium_risk_free_rate: rf,
            conditional_intermediary_equity_premium: eq_prem,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn he_krishnamurthy_intermediary_sdf_invariants_and_ad_parity() {
            let base = HeKrishnamurthyIntermediarySdfSpec {
                intermediary_wealth_share_wi: 0.14_f64,
                capital_constraint_threshold_w_star: 0.28,
                equity_intermediation_pledgeability_m: 0.35,
                specialist_risk_aversion_gamma: 2.20,
                fundamental_growth_rate_ga: 0.020,
                asset_return_volatility_sigma_a: 0.15,
                subjective_discount_rate_rho: 0.025,
            };
            let r0 = evaluate_he_krishnamurthy_intermediary_capital_constraint_sdf(&base);
            assert!(r0.effective_intermediary_leverage_multiplier > 1.0);
            assert!(r0.intermediary_sdf_market_price_of_risk > 0.0);
            assert!(r0.conditional_intermediary_equity_premium > 0.0);

            let inputs = [
                base.intermediary_wealth_share_wi,
                base.capital_constraint_threshold_w_star,
                base.equity_intermediation_pledgeability_m,
                base.specialist_risk_aversion_gamma,
                base.fundamental_growth_rate_ga,
                base.asset_return_volatility_sigma_a,
                base.subjective_discount_rate_rho,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HeKrishnamurthyIntermediarySdfSpec {
                intermediary_wealth_share_wi: dual_in[0],
                capital_constraint_threshold_w_star: dual_in[1],
                equity_intermediation_pledgeability_m: dual_in[2],
                specialist_risk_aversion_gamma: dual_in[3],
                fundamental_growth_rate_ga: dual_in[4],
                asset_return_volatility_sigma_a: dual_in[5],
                subjective_discount_rate_rho: dual_in[6],
            };
            let r_ad = evaluate_he_krishnamurthy_intermediary_capital_constraint_sdf(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| HeKrishnamurthyIntermediarySdfSpec {
                    intermediary_wealth_share_wi: v[0],
                    capital_constraint_threshold_w_star: v[1],
                    equity_intermediation_pledgeability_m: v[2],
                    specialist_risk_aversion_gamma: v[3],
                    fundamental_growth_rate_ga: v[4],
                    asset_return_volatility_sigma_a: v[5],
                    subjective_discount_rate_rho: v[6],
                };
                let rp = evaluate_he_krishnamurthy_intermediary_capital_constraint_sdf(&mk(p));
                let rm = evaluate_he_krishnamurthy_intermediary_capital_constraint_sdf(&mk(m));
                // Test sum of equilibrium_risk_free_rate + conditional_intermediary_equity_premium so all 7 parameters (including ga and rho) have non-zero sensitivity
                let fd = ((rp.equilibrium_risk_free_rate
                    + rp.conditional_intermediary_equity_premium)
                    - (rm.equilibrium_risk_free_rate + rm.conditional_intermediary_equity_premium))
                    / (2.0 * eps);
                let ad_val =
                    r_ad.equilibrium_risk_free_rate + r_ad.conditional_intermediary_equity_premium;
                assert!(ad_val.d[i].abs() > 1e-6, "zero derivative at {i}");
                assert!((ad_val.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use he_krishnamurthy_intermediary_sdf_mod::*;

pub mod diamond_dybvig_goldstein_pauzner_bank_run_mod {
    use ql_ad::Real;
    use ql_math::norm_cdf;

    /// Specification for the Diamond & Dybvig (1983, JPE) & Goldstein & Pauzner (2005, JF)
    /// global-games bank-run threshold, liquidity fire-sale, and deposit fragility engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct DiamondDybvigGoldsteinPauznerBankRunSpec<T: Real> {
        /// Promised short-term demand deposit gross payout $r_1 > 1$ at $t=1$.
        pub short_term_deposit_payout_r1: T,
        /// Long-term asset productivity scale $R_{\text{bar}} > r_1$ at $t=2$.
        pub long_term_asset_return_r_bar: T,
        /// Liquid cash reserve share $\ell \in (0, 1)$ held by the bank.
        pub liquid_reserve_ratio_l: T,
        /// Early fire-sale recovery fraction $\psi_{\text{fire}} \in (0, 1)$ on illiquid assets.
        pub fire_sale_recovery_fraction_psi: T,
        /// Fraction $\lambda_{\text{imp}} \in (0, 1)$ of structurally impatient depositors.
        pub impatient_depositor_fraction_lambda: T,
        /// Prior mean $\mu_\theta$ of macroeconomic asset fundamental $\theta$.
        pub fundamental_prior_mean_mu_theta: T,
        /// Prior volatility $\sigma_\theta > 0$ of macroeconomic asset fundamental $\theta$.
        pub fundamental_volatility_sigma_theta: T,
    }

    /// Analytical outputs for the Diamond-Dybvig & Goldstein-Pauzner bank-run global game engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct DiamondDybvigGoldsteinPauznerBankRunAnalytics<T: Real> {
        /// Total early liquidation capacity $L_{\text{fire}} = \ell + (1 - \ell)\psi_{\text{fire}}$.
        pub early_liquidation_capacity_l_fire: T,
        /// Unique global-games equilibrium run threshold $\theta^*$ below which patient depositors run.
        pub critical_fundamental_run_threshold_theta_star: T,
        /// Ex-ante probability of a systemic bank run $P_{\text{run}} = \Phi((\theta^* - \mu_\theta)/\sigma_\theta)$.
        pub ex_ante_bank_run_probability: T,
        /// Expected depositor payout welfare across survival and run states.
        pub expected_depositor_welfare_value: T,
    }

    /// Evaluates the Diamond & Dybvig (1983) / Goldstein & Pauzner (2005) unique global-games
    /// bank-run threshold $\theta^*$, run probability, and expected depositor welfare.
    pub fn evaluate_diamond_dybvig_goldstein_pauzner_bank_run_global_game<T: Real>(
        spec: &DiamondDybvigGoldsteinPauznerBankRunSpec<T>,
    ) -> DiamondDybvigGoldsteinPauznerBankRunAnalytics<T> {
        let r1 = spec
            .short_term_deposit_payout_r1
            .max_val(T::from_f64(1.001));
        let r_bar = spec
            .long_term_asset_return_r_bar
            .max_val(r1 * T::from_f64(1.05));
        let l_res = spec
            .liquid_reserve_ratio_l
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.95));
        let psi = spec
            .fire_sale_recovery_fraction_psi
            .max_val(T::from_f64(0.10))
            .min_val(T::from_f64(0.98));
        let lam = spec
            .impatient_depositor_fraction_lambda
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.85));
        let sig_theta = spec
            .fundamental_volatility_sigma_theta
            .max_val(T::from_f64(1e-4));

        // Total early liquidation proceeds at t = 1 per unit asset
        let l_fire = l_res + (T::one() - l_res) * psi;

        // Goldstein & Pauzner (2005) Laplacian indifference threshold theta*:
        // A patient depositor integrates payoff difference over uniform withdrawal fraction n in [lam, 1].
        let illiquidity_gap = (r1 - l_fire).max_val(T::from_f64(0.01));
        let strategic_amp = (T::one() + (r1 / l_fire) * (T::one() - lam))
            / (T::one() - (lam * r1 / r_bar)).max_val(T::from_f64(0.10));
        let theta_star = (illiquidity_gap / r_bar) * strategic_amp;

        let z_run = (theta_star - spec.fundamental_prior_mean_mu_theta) / sig_theta;
        let p_run = norm_cdf(z_run);

        // Expected welfare: in a run, depositors share early liquidation value l_fire;
        // in survival, impatient withdraw lam * r1 and remaining illiquid assets are
        // (1 - l_res) - (lam * r1 - l_res) / psi = (l_fire - lam * r1) / psi.
        let net_fire_surplus = l_fire - lam * r1;
        let smooth_surplus = T::from_f64(0.5)
            * (net_fire_surplus + (net_fire_surplus * net_fire_surplus + T::from_f64(1e-4)).sqrt());
        let survival_payout = lam * r1
            + (T::one() - lam)
                * spec.fundamental_prior_mean_mu_theta
                * r_bar
                * (smooth_surplus / psi);
        let expected_welfare = p_run * l_fire + (T::one() - p_run) * survival_payout;

        DiamondDybvigGoldsteinPauznerBankRunAnalytics {
            early_liquidation_capacity_l_fire: l_fire,
            critical_fundamental_run_threshold_theta_star: theta_star,
            ex_ante_bank_run_probability: p_run,
            expected_depositor_welfare_value: expected_welfare,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn diamond_dybvig_goldstein_pauzner_invariants_and_ad_parity() {
            let base = DiamondDybvigGoldsteinPauznerBankRunSpec {
                short_term_deposit_payout_r1: 1.06_f64,
                long_term_asset_return_r_bar: 1.55,
                liquid_reserve_ratio_l: 0.22,
                fire_sale_recovery_fraction_psi: 0.65,
                impatient_depositor_fraction_lambda: 0.30,
                fundamental_prior_mean_mu_theta: 0.85,
                fundamental_volatility_sigma_theta: 0.18,
            };
            let r0 = evaluate_diamond_dybvig_goldstein_pauzner_bank_run_global_game(&base);
            assert!(
                r0.early_liquidation_capacity_l_fire > base.fire_sale_recovery_fraction_psi
                    && r0.early_liquidation_capacity_l_fire < 1.0
            );
            assert!(r0.critical_fundamental_run_threshold_theta_star > 0.0);
            assert!(r0.ex_ante_bank_run_probability > 0.0 && r0.ex_ante_bank_run_probability < 1.0);
            assert!(r0.expected_depositor_welfare_value > r0.early_liquidation_capacity_l_fire);

            let inputs = [
                base.short_term_deposit_payout_r1,
                base.long_term_asset_return_r_bar,
                base.liquid_reserve_ratio_l,
                base.fire_sale_recovery_fraction_psi,
                base.impatient_depositor_fraction_lambda,
                base.fundamental_prior_mean_mu_theta,
                base.fundamental_volatility_sigma_theta,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = DiamondDybvigGoldsteinPauznerBankRunSpec {
                short_term_deposit_payout_r1: dual_in[0],
                long_term_asset_return_r_bar: dual_in[1],
                liquid_reserve_ratio_l: dual_in[2],
                fire_sale_recovery_fraction_psi: dual_in[3],
                impatient_depositor_fraction_lambda: dual_in[4],
                fundamental_prior_mean_mu_theta: dual_in[5],
                fundamental_volatility_sigma_theta: dual_in[6],
            };
            let r_ad = evaluate_diamond_dybvig_goldstein_pauzner_bank_run_global_game(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| DiamondDybvigGoldsteinPauznerBankRunSpec {
                    short_term_deposit_payout_r1: v[0],
                    long_term_asset_return_r_bar: v[1],
                    liquid_reserve_ratio_l: v[2],
                    fire_sale_recovery_fraction_psi: v[3],
                    impatient_depositor_fraction_lambda: v[4],
                    fundamental_prior_mean_mu_theta: v[5],
                    fundamental_volatility_sigma_theta: v[6],
                };
                let rp = evaluate_diamond_dybvig_goldstein_pauzner_bank_run_global_game(&mk(p));
                let rm = evaluate_diamond_dybvig_goldstein_pauzner_bank_run_global_game(&mk(m));
                let fd = (rp.expected_depositor_welfare_value
                    - rm.expected_depositor_welfare_value)
                    / (2.0 * eps);
                assert!(
                    r_ad.expected_depositor_welfare_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.expected_depositor_welfare_value.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use diamond_dybvig_goldstein_pauzner_bank_run_mod::*;

// ============================================================================
// Barro (2006) & Gabaix (2012) Variable Rare Disasters & Resilience Equity Premium
// ============================================================================

pub mod gabaix_rare_disasters_resilience_mod {
    use ql_ad::Real;

    /// Specification for the Barro (2006, *QJE* 121(3)) and Gabaix (2012, *QJE* 127(2))
    /// variable rare-disasters macro-finance model with time-varying asset resilience $H_t$.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GabaixRareDisastersResilienceSpec<T: Real> {
        /// Representative agent subjective time-preference discount rate $\rho_{\text{pref}} > 0$.
        pub subjective_discount_rate_rho: T,
        /// Coefficient of relative risk aversion $\gamma > 1$.
        pub relative_risk_aversion_gamma: T,
        /// Normal-times consumption and dividend growth rate $g_C > 0$.
        pub normal_consumption_growth_g_c: T,
        /// Normal-times diffusive consumption volatility $\sigma_C > 0$.
        pub normal_consumption_volatility_sigma_c: T,
        /// Annual arrival probability of a rare macroeconomic disaster $p_{\text{dis}} > 0$.
        pub disaster_probability_p: T,
        /// Consumption recovery fraction in a disaster $B_C \in (0, 1)$ ($C_{t+}/C_{t-} = B_C$).
        pub consumption_recovery_fraction_b_c: T,
        /// Equity dividend recovery fraction in a disaster $B_D \in (0, 1)$ ($D_{t+}/D_{t-} = B_D$).
        pub dividend_recovery_fraction_b_d: T,
        /// Mean-reversion speed of variable equity resilience $\phi_H > 0$.
        pub resilience_mean_reversion_phi_h: T,
    }

    /// Analytical outputs for the Barro-Gabaix variable rare-disasters asset-pricing engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GabaixRareDisastersResilienceAnalytics<T: Real> {
        /// Equilibrium real risk-free rate $r_f$ including Barro precautionary disaster savings.
        pub equilibrium_real_risk_free_rate: T,
        /// Steady-state asset resilience $H_* = p_{\text{dis}}(B_C^{-\gamma} B_D - 1)$.
        pub baseline_asset_resilience_h_star: T,
        /// Gabaix variable-resilience equity risk premium $\text{ERP}_{\text{Gabaix}}$.
        pub rare_disaster_equity_risk_premium: T,
        /// Closed-form linearity-generating price-dividend ratio $P_0 / D_0$.
        pub gabaix_equilibrium_price_dividend_ratio: T,
    }

    /// Evaluates the Barro (2006) / Gabaix (2012, Theorem 1 & Proposition 1) rare-disasters
    /// real risk-free rate, asset resilience $H_*$, equity risk premium, and price-dividend ratio.
    pub fn evaluate_gabaix_rare_disasters_time_varying_resilience_equity_premium<T: Real>(
        spec: &GabaixRareDisastersResilienceSpec<T>,
    ) -> GabaixRareDisastersResilienceAnalytics<T> {
        let gamma = spec.relative_risk_aversion_gamma.max_val(T::from_f64(1.05));
        let p_dis = spec.disaster_probability_p.max_val(T::from_f64(1e-4));
        let b_c = spec
            .consumption_recovery_fraction_b_c
            .max_val(T::from_f64(0.20))
            .min_val(T::from_f64(0.98));
        let b_d = spec
            .dividend_recovery_fraction_b_d
            .max_val(T::from_f64(0.10))
            .min_val(T::from_f64(0.98));
        let phi_h = spec
            .resilience_mean_reversion_phi_h
            .max_val(T::from_f64(1e-3));
        let sig_c = spec
            .normal_consumption_volatility_sigma_c
            .max_val(T::from_f64(1e-4));

        let bc_neg_gamma = b_c.pow(-gamma);

        // Barro (2006) / Gabaix (2012, Eq. 8) real risk-free rate
        let rf = spec.subjective_discount_rate_rho + gamma * spec.normal_consumption_growth_g_c
            - T::from_f64(0.5) * gamma * (gamma + T::one()) * sig_c * sig_c
            - p_dis * (bc_neg_gamma - T::one());

        // Gabaix (2012, Eq. 9) baseline asset resilience H_*
        let h_star = p_dis * (bc_neg_gamma * b_d - T::one());

        // Effective steady-state equity discount spread r_eq - g_C (Gabaix 2012, Eq. 11)
        let req_minus_g = (spec.subjective_discount_rate_rho
            + (gamma - T::one()) * spec.normal_consumption_growth_g_c
            - T::from_f64(0.5) * (gamma - T::one()) * gamma * sig_c * sig_c
            - h_star)
            .max_val(T::from_f64(0.005));

        // Active resilience deviation \hat{H}_0 = 0.25 * |H_*| in Gabaix's linearity-generating process
        let h_hat0 = T::from_f64(0.25) * (h_star * h_star + T::from_f64(1e-8)).sqrt();

        // Gabaix (2012, Theorem 1, Eq. 12) exact price-dividend ratio
        let pd_ratio = (T::one() + h_hat0 / (req_minus_g + phi_h)) / req_minus_g;

        // Gabaix (2012, Proposition 1, Eq. 14) equity risk premium
        let erp = gamma * sig_c * sig_c
            + p_dis * bc_neg_gamma * (T::one() - b_d)
            + (phi_h * h_hat0) / (req_minus_g + phi_h);

        GabaixRareDisastersResilienceAnalytics {
            equilibrium_real_risk_free_rate: rf,
            baseline_asset_resilience_h_star: h_star,
            rare_disaster_equity_risk_premium: erp,
            gabaix_equilibrium_price_dividend_ratio: pd_ratio,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gabaix_rare_disasters_invariants_and_ad_parity() {
            let base = GabaixRareDisastersResilienceSpec {
                subjective_discount_rate_rho: 0.032_f64,
                relative_risk_aversion_gamma: 4.0,
                normal_consumption_growth_g_c: 0.020,
                normal_consumption_volatility_sigma_c: 0.018,
                disaster_probability_p: 0.017,
                consumption_recovery_fraction_b_c: 0.72,
                dividend_recovery_fraction_b_d: 0.18,
                resilience_mean_reversion_phi_h: 0.18,
            };
            let r0 = evaluate_gabaix_rare_disasters_time_varying_resilience_equity_premium(&base);
            assert!(r0.baseline_asset_resilience_h_star < 0.0);
            assert!(r0.rare_disaster_equity_risk_premium > 0.02);
            assert!(r0.gabaix_equilibrium_price_dividend_ratio > 5.0);

            let inputs = [
                base.subjective_discount_rate_rho,
                base.relative_risk_aversion_gamma,
                base.normal_consumption_growth_g_c,
                base.normal_consumption_volatility_sigma_c,
                base.disaster_probability_p,
                base.consumption_recovery_fraction_b_c,
                base.dividend_recovery_fraction_b_d,
                base.resilience_mean_reversion_phi_h,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GabaixRareDisastersResilienceSpec {
                subjective_discount_rate_rho: dual_in[0],
                relative_risk_aversion_gamma: dual_in[1],
                normal_consumption_growth_g_c: dual_in[2],
                normal_consumption_volatility_sigma_c: dual_in[3],
                disaster_probability_p: dual_in[4],
                consumption_recovery_fraction_b_c: dual_in[5],
                dividend_recovery_fraction_b_d: dual_in[6],
                resilience_mean_reversion_phi_h: dual_in[7],
            };
            let r_ad =
                evaluate_gabaix_rare_disasters_time_varying_resilience_equity_premium(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| GabaixRareDisastersResilienceSpec {
                    subjective_discount_rate_rho: v[0],
                    relative_risk_aversion_gamma: v[1],
                    normal_consumption_growth_g_c: v[2],
                    normal_consumption_volatility_sigma_c: v[3],
                    disaster_probability_p: v[4],
                    consumption_recovery_fraction_b_c: v[5],
                    dividend_recovery_fraction_b_d: v[6],
                    resilience_mean_reversion_phi_h: v[7],
                };
                let rp =
                    evaluate_gabaix_rare_disasters_time_varying_resilience_equity_premium(&mk(p));
                let rm =
                    evaluate_gabaix_rare_disasters_time_varying_resilience_equity_premium(&mk(m));
                let fd = (rp.gabaix_equilibrium_price_dividend_ratio
                    - rm.gabaix_equilibrium_price_dividend_ratio)
                    / (2.0 * eps);
                assert!(
                    r_ad.gabaix_equilibrium_price_dividend_ratio.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.gabaix_equilibrium_price_dividend_ratio.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use gabaix_rare_disasters_resilience_mod::*;

// ============================================================================
// Piazzesi & Schneider (2007) & Bansal-Shaliastovich (2013) Stagflation Nominal Term Premium
// ============================================================================

pub mod piazzesi_schneider_nominal_tp_mod {
    use ql_ad::Real;

    /// Specification for the Piazzesi & Schneider (2007, *NBER Macro Annual* 21) and
    /// Bansal & Shaliastovich (2013, *RFS* 26(1)) Epstein-Zin nominal bond yield and
    /// stagflation inflation risk premium engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct PiazzesiSchneiderNominalBondSpec<T: Real> {
        /// Baseline real short interest rate $r_{\text{real}}$.
        pub real_short_rate_r_real: T,
        /// Current expected inflation rate $\pi_0$.
        pub expected_inflation_pi0: T,
        /// Epstein-Zin relative risk aversion $\gamma > 1$.
        pub relative_risk_aversion_gamma: T,
        /// Persistence of expected real consumption growth $\phi_x \in (0, 1)$.
        pub real_growth_persistence_phi_x: T,
        /// Persistence of expected inflation $\phi_\pi \in (0, 1)$.
        pub inflation_persistence_phi_pi: T,
        /// Conditional volatility of expected real growth news $\sigma_x > 0$.
        pub real_growth_volatility_sigma_x: T,
        /// Conditional volatility of expected inflation news $\sigma_\pi > 0$.
        pub inflation_volatility_sigma_pi: T,
        /// Correlation between real growth news and inflation news $\rho_{x,\pi} \in (-1, 1)$
        /// (negative in stagflation regimes where inflation spikes coincide with low growth).
        pub growth_inflation_correlation_rho_x_pi: T,
    }

    /// Analytical outputs for the Piazzesi-Schneider nominal bond term premium engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct PiazzesiSchneiderNominalBondAnalytics<T: Real> {
        /// Real bond risk premium $\text{RP}_{\text{real}}(10\text{Y})$.
        pub real_bond_term_premium: T,
        /// Stagflation-induced inflation risk premium $\text{IRP}_{\text{stag}}(10\text{Y})$.
        pub stagflation_inflation_risk_premium: T,
        /// Total 10-year nominal bond term premium $\text{TP}_{\text{nom}}(10\text{Y})$.
        pub nominal_bond_term_premium: T,
        /// Equilibrium 10-year nominal zero-coupon bond yield $y_{\text{nom}}(10\text{Y})$.
        pub equilibrium_long_nominal_bond_yield: T,
    }

    /// Evaluates the Piazzesi & Schneider (2007) / Bansal & Shaliastovich (2013) 10-year
    /// real bond risk premium, stagflation inflation risk premium, and nominal yield.
    pub fn evaluate_piazzesi_schneider_stagflation_nominal_bond_term_premium<T: Real>(
        spec: &PiazzesiSchneiderNominalBondSpec<T>,
    ) -> PiazzesiSchneiderNominalBondAnalytics<T> {
        let gamma = spec.relative_risk_aversion_gamma.max_val(T::from_f64(1.05));
        let phi_x = spec
            .real_growth_persistence_phi_x
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.99));
        let phi_pi = spec
            .inflation_persistence_phi_pi
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.99));
        let sig_x = spec
            .real_growth_volatility_sigma_x
            .max_val(T::from_f64(1e-5));
        let sig_pi = spec
            .inflation_volatility_sigma_pi
            .max_val(T::from_f64(1e-5));
        let rho = spec
            .growth_inflation_correlation_rho_x_pi
            .max_val(T::from_f64(-0.98))
            .min_val(T::from_f64(0.98));

        let n_years = T::from_f64(10.0);
        let b_x_wealth = (T::one() - phi_x.powf(5.0)) / (T::one() - phi_x);
        let b_x_10 = (T::one() - phi_x.powf(10.0)) / (T::one() - phi_x);
        let b_pi_10 = (T::one() - phi_pi.powf(10.0)) / (T::one() - phi_pi);

        // Real growth news market price of risk Lambda_x = gamma * B_x(5) * sigma_x
        let real_tp = (gamma * b_x_wealth * b_x_10 * sig_x * sig_x) / n_years;

        // Stagflation inflation risk premium: positive when rho_{x,pi} < 0 (bad consumption hedge)
        let irp_stag = (-gamma * b_x_wealth * b_pi_10 * rho * sig_x * sig_pi
            - T::from_f64(0.5) * b_pi_10 * b_pi_10 * sig_pi * sig_pi)
            / n_years;

        let nom_tp = real_tp + irp_stag;
        let expected_avg_inf = (b_pi_10 / n_years) * spec.expected_inflation_pi0;
        let y_nom_10 = spec.real_short_rate_r_real + expected_avg_inf + nom_tp;

        PiazzesiSchneiderNominalBondAnalytics {
            real_bond_term_premium: real_tp,
            stagflation_inflation_risk_premium: irp_stag,
            nominal_bond_term_premium: nom_tp,
            equilibrium_long_nominal_bond_yield: y_nom_10,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn piazzesi_schneider_nominal_tp_invariants_and_ad_parity() {
            let base = PiazzesiSchneiderNominalBondSpec {
                real_short_rate_r_real: 0.015_f64,
                expected_inflation_pi0: 0.025,
                relative_risk_aversion_gamma: 8.5,
                real_growth_persistence_phi_x: 0.92,
                inflation_persistence_phi_pi: 0.94,
                real_growth_volatility_sigma_x: 0.006,
                inflation_volatility_sigma_pi: 0.008,
                growth_inflation_correlation_rho_x_pi: -0.55,
            };
            let r0 = evaluate_piazzesi_schneider_stagflation_nominal_bond_term_premium(&base);
            assert!(r0.real_bond_term_premium > 0.0);
            assert!(r0.stagflation_inflation_risk_premium > 0.0);
            assert!(r0.nominal_bond_term_premium > 0.0);
            assert!(r0.equilibrium_long_nominal_bond_yield > base.real_short_rate_r_real);

            let inputs = [
                base.real_short_rate_r_real,
                base.expected_inflation_pi0,
                base.relative_risk_aversion_gamma,
                base.real_growth_persistence_phi_x,
                base.inflation_persistence_phi_pi,
                base.real_growth_volatility_sigma_x,
                base.inflation_volatility_sigma_pi,
                base.growth_inflation_correlation_rho_x_pi,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = PiazzesiSchneiderNominalBondSpec {
                real_short_rate_r_real: dual_in[0],
                expected_inflation_pi0: dual_in[1],
                relative_risk_aversion_gamma: dual_in[2],
                real_growth_persistence_phi_x: dual_in[3],
                inflation_persistence_phi_pi: dual_in[4],
                real_growth_volatility_sigma_x: dual_in[5],
                inflation_volatility_sigma_pi: dual_in[6],
                growth_inflation_correlation_rho_x_pi: dual_in[7],
            };
            let r_ad = evaluate_piazzesi_schneider_stagflation_nominal_bond_term_premium(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| PiazzesiSchneiderNominalBondSpec {
                    real_short_rate_r_real: v[0],
                    expected_inflation_pi0: v[1],
                    relative_risk_aversion_gamma: v[2],
                    real_growth_persistence_phi_x: v[3],
                    inflation_persistence_phi_pi: v[4],
                    real_growth_volatility_sigma_x: v[5],
                    inflation_volatility_sigma_pi: v[6],
                    growth_inflation_correlation_rho_x_pi: v[7],
                };
                let rp = evaluate_piazzesi_schneider_stagflation_nominal_bond_term_premium(&mk(p));
                let rm = evaluate_piazzesi_schneider_stagflation_nominal_bond_term_premium(&mk(m));
                let fd = (rp.equilibrium_long_nominal_bond_yield
                    - rm.equilibrium_long_nominal_bond_yield)
                    / (2.0 * eps);
                assert!(
                    r_ad.equilibrium_long_nominal_bond_yield.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.equilibrium_long_nominal_bond_yield.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use piazzesi_schneider_nominal_tp_mod::*;
mod jermann_quadrini_financial_shock_mod {
    use ql_ad::Real;

    /// Input specification for the Jermann & Quadrini (2012, *American Economic Review* 102(1))
    /// macroeconomic model of financial shocks, debt enforcement constraints, and labor wedges.
    #[derive(Debug, Clone, Copy)]
    pub struct JermannQuadriniFinancialShockSpec<T: Real> {
        /// Firm productive capital stock $k_t$ (e.g., $10.0$).
        pub capital_stock_k: T,
        /// Gross intra-period output / working capital requirement $y_t$ (e.g., $2.4$).
        pub working_capital_output_y: T,
        /// Financial enforcement pledgeability parameter $\xi_{\text{fin}} \in (0, 1)$ (e.g., $0.42$).
        pub financial_enforcement_xi: T,
        /// Risk-free interest rate $r$ (e.g., $0.032$).
        pub risk_free_rate_r: T,
        /// Corporate debt interest tax-shield rate $\tau_{\text{tax}} \in (0, 1)$ (e.g., $0.28$).
        pub corporate_tax_shield_tau: T,
        /// Representative household stochastic discount factor $\beta_{\text{disc}} \in (0, 1)$ (e.g., $0.975$).
        pub household_discount_factor_beta: T,
        /// Quadratic equity payout adjustment cost parameter $\varphi_{\text{eq}} > 0$ (e.g., $0.18$).
        pub equity_payout_cost_phi: T,
        /// Equity payout gap $d_t - \bar{d}$ (negative during equity issuance / tightening, e.g., $-0.15$).
        pub equity_payout_gap_d: T,
    }

    /// Analytics returned by [`evaluate_jermann_quadrini_financial_shocks_enforcement_constraint_credit_cycle`].
    #[derive(Debug, Clone, Copy)]
    pub struct JermannQuadriniFinancialShockAnalytics<T: Real> {
        /// Binding enforcement-constrained intertemporal debt capacity $b_{\text{bind}}$.
        pub binding_enforcement_debt_capacity: T,
        /// Equilibrium debt-to-capital leverage ratio $b_{\text{bind}} / k_t$.
        pub equilibrium_debt_to_capital_leverage: T,
        /// Shadow Lagrange multiplier $\mu_{\text{bind}}$ on the enforcement constraint.
        pub shadow_enforcement_multiplier_mu: T,
        /// Jermann–Quadrini endogenous financial labor efficiency wedge $\omega_{\text{labor}} \in (0, 1)$.
        pub endogenous_financial_labor_wedge: T,
    }

    /// Evaluates the Jermann & Quadrini (2012, Section I, Eqs. 4–12) corporate enforcement constraint
    /// $\xi_{\text{fin}}(k_t - b_t / (1 + r)) \ge y_t$, marginal cost of equity payout adjustment
    /// $1 + 2\varphi_{\text{eq}}(d_t - \bar{d})$, shadow multiplier $\mu_{\text{bind}}$, and the resulting
    /// endogenous labor demand wedge $\omega_{\text{labor}}$.
    pub fn evaluate_jermann_quadrini_financial_shocks_enforcement_constraint_credit_cycle<
        T: Real,
    >(
        spec: &JermannQuadriniFinancialShockSpec<T>,
    ) -> JermannQuadriniFinancialShockAnalytics<T> {
        let gross_r = T::one() + spec.risk_free_rate_r;
        // Effective after-tax gross debt cost R_eff = 1 + r * (1 - tau_tax)
        let r_eff = T::one() + spec.risk_free_rate_r * (T::one() - spec.corporate_tax_shield_tau);

        // Binding enforcement constraint: xi * (k - b / (1 + r)) = y  =>  b = (1 + r) * (k - y / xi)
        let b_bind = gross_r
            * (spec.capital_stock_k
                - spec.working_capital_output_y / spec.financial_enforcement_xi);
        let leverage = b_bind / spec.capital_stock_k;

        // Intratemporal working-capital financing tightness theta_wc = y_t / (xi_fin * k_t)
        let wc_tightness =
            spec.working_capital_output_y / (spec.financial_enforcement_xi * spec.capital_stock_k);
        let effective_payout_gap = spec.equity_payout_gap_d - T::from_f64(0.5) * wc_tightness;

        // Marginal cost of equity funds phi'(d_t) = 1 + 2 * phi_eq * (d_t - d_bar - 0.5 * theta_wc)
        let marg_eq_cost =
            T::one() + T::from_f64(2.0) * spec.equity_payout_cost_phi * effective_payout_gap;

        // First-order condition for debt (Eq. 9):
        // 1 / (R_eff * marg_eq_cost) - beta = xi * mu / (1 + r)
        let mu_bind = (gross_r / spec.financial_enforcement_xi)
            * (T::one() / (r_eff * marg_eq_cost) - spec.household_discount_factor_beta);

        // Endogenous labor wedge (Eq. 10-11): MPL = w * (1 + mu * marg_eq_cost) => wedge = (mu * marg_eq_cost) / (1 + mu * marg_eq_cost)
        let dist = mu_bind * marg_eq_cost;
        let labor_wedge = dist / (T::one() + dist);

        JermannQuadriniFinancialShockAnalytics {
            binding_enforcement_debt_capacity: b_bind,
            equilibrium_debt_to_capital_leverage: leverage,
            shadow_enforcement_multiplier_mu: mu_bind,
            endogenous_financial_labor_wedge: labor_wedge,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn jermann_quadrini_financial_shock_invariants_and_ad_parity() {
            let base = JermannQuadriniFinancialShockSpec {
                capital_stock_k: 10.0_f64,
                working_capital_output_y: 2.4_f64,
                financial_enforcement_xi: 0.42_f64,
                risk_free_rate_r: 0.032_f64,
                corporate_tax_shield_tau: 0.28_f64,
                household_discount_factor_beta: 0.960_f64,
                equity_payout_cost_phi: 0.18_f64,
                equity_payout_gap_d: -0.15_f64,
            };
            let res =
                evaluate_jermann_quadrini_financial_shocks_enforcement_constraint_credit_cycle(
                    &base,
                );
            assert!(res.binding_enforcement_debt_capacity > 0.0);
            assert!(
                res.equilibrium_debt_to_capital_leverage > 0.0
                    && res.equilibrium_debt_to_capital_leverage < 1.0
            );
            assert!(res.shadow_enforcement_multiplier_mu > 0.0);
            assert!(
                res.endogenous_financial_labor_wedge > 0.0
                    && res.endogenous_financial_labor_wedge < 0.5
            );

            let inputs = [
                base.capital_stock_k,
                base.working_capital_output_y,
                base.financial_enforcement_xi,
                base.risk_free_rate_r,
                base.corporate_tax_shield_tau,
                base.household_discount_factor_beta,
                base.equity_payout_cost_phi,
                base.equity_payout_gap_d,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = JermannQuadriniFinancialShockSpec {
                capital_stock_k: dual_in[0],
                working_capital_output_y: dual_in[1],
                financial_enforcement_xi: dual_in[2],
                risk_free_rate_r: dual_in[3],
                corporate_tax_shield_tau: dual_in[4],
                household_discount_factor_beta: dual_in[5],
                equity_payout_cost_phi: dual_in[6],
                equity_payout_gap_d: dual_in[7],
            };
            let r_ad =
                evaluate_jermann_quadrini_financial_shocks_enforcement_constraint_credit_cycle(
                    &spec_ad,
                );
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| JermannQuadriniFinancialShockSpec {
                    capital_stock_k: v[0],
                    working_capital_output_y: v[1],
                    financial_enforcement_xi: v[2],
                    risk_free_rate_r: v[3],
                    corporate_tax_shield_tau: v[4],
                    household_discount_factor_beta: v[5],
                    equity_payout_cost_phi: v[6],
                    equity_payout_gap_d: v[7],
                };
                let rp =
                    evaluate_jermann_quadrini_financial_shocks_enforcement_constraint_credit_cycle(
                        &mk(p),
                    );
                let rm =
                    evaluate_jermann_quadrini_financial_shocks_enforcement_constraint_credit_cycle(
                        &mk(m),
                    );
                let fd = (rp.endogenous_financial_labor_wedge
                    - rm.endogenous_financial_labor_wedge)
                    / (2.0 * eps);
                assert!(
                    r_ad.endogenous_financial_labor_wedge.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.endogenous_financial_labor_wedge.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use jermann_quadrini_financial_shock_mod::*;
mod wachter_disaster_volatility_mod {
    use ql_ad::Real;

    /// Input specification for the Wachter (2013, *Journal of Finance* 68(3))
    /// continuous-time CIR rare-disaster intensity model of aggregate stock market volatility
    /// and the equity risk premium.
    #[derive(Debug, Clone, Copy)]
    pub struct WachterDisasterVolatilitySpec<T: Real> {
        /// Current rare-disaster arrival intensity $\lambda_t > 0$ (e.g., $0.0355$).
        pub current_disaster_intensity_lambda: T,
        /// Long-run mean disaster intensity $\bar{\lambda} > 0$ (e.g., $0.0355$).
        pub long_run_disaster_intensity_lambda_bar: T,
        /// Disaster intensity mean-reversion speed $\kappa_\lambda > 0$ (e.g., $0.080$).
        pub intensity_mean_reversion_kappa: T,
        /// Square-root volatility of disaster intensity $\sigma_\lambda > 0$ (e.g., $0.067$).
        pub intensity_volatility_sigma_lambda: T,
        /// Normal-times diffusive dividend growth volatility $\sigma_c$ (e.g., $0.020$).
        pub diffusive_consumption_volatility_sigma_c: T,
        /// Relative risk aversion $\gamma > 1$ under unit EIS ($\psi = 1$) (e.g., $3.0$).
        pub relative_risk_aversion_gamma: T,
        /// Fractional consumption/dividend contraction in a disaster $b_{\text{dis}} \in (0, 1)$ (e.g., $0.22$).
        pub disaster_contraction_fraction_b: T,
        /// Equity dividend leverage parameter $\phi_{\text{lev}} > 1$ (e.g., $2.6$).
        pub dividend_leverage_phi: T,
    }

    /// Analytics returned by [`evaluate_wachter_time_varying_disaster_intensity_equity_volatility`].
    #[derive(Debug, Clone, Copy)]
    pub struct WachterDisasterVolatilityAnalytics<T: Real> {
        /// Positive SDF/marginal-utility sensitivity $b_\lambda > 0$ to disaster intensity $\lambda_t$ when $\gamma > 1$.
        pub wealth_intensity_loading_b_lambda: T,
        /// Negative levered equity valuation loading $b_{\text{eq}} < 0$ on disaster intensity $\lambda_t$.
        pub equity_intensity_loading_b_eq: T,
        /// Equilibrium instantaneous stock return volatility $\sigma_R(\lambda_t)$.
        pub equilibrium_equity_return_volatility: T,
        /// Equilibrium conditional equity risk premium $\text{ERP}(\lambda_t)$ (including normal + disaster-intensity + jump premia).
        pub equilibrium_equity_risk_premium: T,
    }

    /// Evaluates the Wachter (2013, *Journal of Finance* 68(3), Eqs. 12–23) rare-disaster model
    /// with time-varying CIR intensity $d\lambda_t = \kappa_\lambda(\bar{\lambda} - \lambda_t)dt + \sigma_\lambda\sqrt{\lambda_t}dW_{\lambda,t}$:
    /// - Under unit EIS ($\psi = 1$), the SDF jump multiplier is $(1 - b_{\text{dis}})^{-\gamma}$ and marginal utility rises with $\lambda_t$ ($b_\lambda > 0$),
    /// - The levered equity dividend strip loading is $b_{\text{eq}} < 0$ (stock prices fall when disaster intensity rises),
    /// - Return variance is $\sigma_R^2(\lambda_t) = \phi_{\text{lev}}^2 \sigma_c^2 + b_{\text{eq}}^2 \sigma_\lambda^2 \lambda_t + \lambda_t \big((1 - b_{\text{dis}})^{\phi_{\text{lev}}} - 1\big)^2$, and
    /// - The conditional equity risk premium is $\text{ERP}(\lambda_t) = \gamma \phi_{\text{lev}} \sigma_c^2 + (-b_{\text{eq}}) b_\lambda \sigma_\lambda^2 \lambda_t + \lambda_t \big((1 - b_{\text{dis}})^{-\gamma} - 1\big)\big(1 - (1 - b_{\text{dis}})^{\phi_{\text{lev}}}\big)$.
    pub fn evaluate_wachter_time_varying_disaster_intensity_equity_volatility<T: Real>(
        spec: &WachterDisasterVolatilitySpec<T>,
    ) -> WachterDisasterVolatilityAnalytics<T> {
        let lam = spec.current_disaster_intensity_lambda;
        let lam_bar = spec.long_run_disaster_intensity_lambda_bar;
        let kappa = spec.intensity_mean_reversion_kappa;
        let sig_l = spec.intensity_volatility_sigma_lambda;
        let sig_c = spec.diffusive_consumption_volatility_sigma_c;
        let gamma = spec.relative_risk_aversion_gamma;
        let one_minus_b = T::one() - spec.disaster_contraction_fraction_b;
        let phi = spec.dividend_leverage_phi;

        // SDF jump factor (1 - b)^(-gamma) - 1 and wealth/marginal-utility intensity sensitivity b_lambda > 0 (Wachter Eq. 12)
        let sdf_jump = one_minus_b.pow(-gamma) - T::one();
        let wealth_jump = one_minus_b.pow(T::one() - gamma) - T::one();
        let disc_anchor = kappa + T::from_f64(0.025) + T::from_f64(0.2) * lam_bar;

        // Quadratic Riccati solution for b_lambda > 0: 0.5*sig_l^2 * b_l^2 - kappa * b_l + wealth_jump = 0
        let disc_ric = (kappa * kappa - T::from_f64(2.0) * sig_l * sig_l * wealth_jump).sqrt();
        // Select the non-explosive positive Riccati branch (Wachter 2013, Eq. 13):
        let b_lambda = (kappa - disc_ric) / (sig_l * sig_l);

        // Equity dividend jump factor (1 - b)^phi - 1 < 0 and equity loading b_eq < 0
        let div_jump = one_minus_b.pow(phi) - T::one();
        let b_eq = (div_jump - wealth_jump) / (disc_anchor + b_lambda * sig_l * sig_l);

        // Total instantaneous equity return volatility (diffusive + disaster-intensity + physical jump variance)
        let eq_var = phi * phi * sig_c * sig_c
            + b_eq * b_eq * sig_l * sig_l * lam
            + lam * div_jump * div_jump;
        let eq_vol = eq_var.sqrt();

        // Conditional equity risk premium (Wachter 2013, Eq. 23):
        // 1. Normal times diffusive premium: gamma * phi * sig_c^2
        // 2. Disaster intensity risk premium: (-b_eq) * b_lambda * sig_l^2 * lam > 0
        // 3. Disaster jump covariance premium: lam * ((1 - b)^(-gamma) - 1) * (1 - (1 - b)^phi) > 0
        let erp = gamma * phi * sig_c * sig_c - b_eq * b_lambda * sig_l * sig_l * lam
            + lam * sdf_jump * (-div_jump);

        WachterDisasterVolatilityAnalytics {
            wealth_intensity_loading_b_lambda: b_lambda,
            equity_intensity_loading_b_eq: b_eq,
            equilibrium_equity_return_volatility: eq_vol,
            equilibrium_equity_risk_premium: erp,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn wachter_disaster_volatility_invariants_and_ad_parity() {
            let base = WachterDisasterVolatilitySpec {
                current_disaster_intensity_lambda: 0.0355_f64,
                long_run_disaster_intensity_lambda_bar: 0.0355_f64,
                intensity_mean_reversion_kappa: 0.080_f64,
                intensity_volatility_sigma_lambda: 0.067_f64,
                diffusive_consumption_volatility_sigma_c: 0.020_f64,
                relative_risk_aversion_gamma: 3.0_f64,
                disaster_contraction_fraction_b: 0.22_f64,
                dividend_leverage_phi: 2.6_f64,
            };
            let res = evaluate_wachter_time_varying_disaster_intensity_equity_volatility(&base);
            assert!(res.wealth_intensity_loading_b_lambda > 0.0);
            assert!(res.equity_intensity_loading_b_eq < 0.0);
            assert!(res.equilibrium_equity_return_volatility > 0.12);
            assert!(res.equilibrium_equity_risk_premium > 0.03);

            let inputs = [
                base.current_disaster_intensity_lambda,
                base.long_run_disaster_intensity_lambda_bar,
                base.intensity_mean_reversion_kappa,
                base.intensity_volatility_sigma_lambda,
                base.diffusive_consumption_volatility_sigma_c,
                base.relative_risk_aversion_gamma,
                base.disaster_contraction_fraction_b,
                base.dividend_leverage_phi,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = WachterDisasterVolatilitySpec {
                current_disaster_intensity_lambda: dual_in[0],
                long_run_disaster_intensity_lambda_bar: dual_in[1],
                intensity_mean_reversion_kappa: dual_in[2],
                intensity_volatility_sigma_lambda: dual_in[3],
                diffusive_consumption_volatility_sigma_c: dual_in[4],
                relative_risk_aversion_gamma: dual_in[5],
                disaster_contraction_fraction_b: dual_in[6],
                dividend_leverage_phi: dual_in[7],
            };
            let r_ad = evaluate_wachter_time_varying_disaster_intensity_equity_volatility(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| WachterDisasterVolatilitySpec {
                    current_disaster_intensity_lambda: v[0],
                    long_run_disaster_intensity_lambda_bar: v[1],
                    intensity_mean_reversion_kappa: v[2],
                    intensity_volatility_sigma_lambda: v[3],
                    diffusive_consumption_volatility_sigma_c: v[4],
                    relative_risk_aversion_gamma: v[5],
                    disaster_contraction_fraction_b: v[6],
                    dividend_leverage_phi: v[7],
                };
                let rp = evaluate_wachter_time_varying_disaster_intensity_equity_volatility(&mk(p));
                let rm = evaluate_wachter_time_varying_disaster_intensity_equity_volatility(&mk(m));
                let fd = (rp.equilibrium_equity_risk_premium - rm.equilibrium_equity_risk_premium)
                    / (2.0 * eps);
                assert!(
                    r_ad.equilibrium_equity_risk_premium.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.equilibrium_equity_risk_premium.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use wachter_disaster_volatility_mod::*;
mod christiano_eichenbaum_evans_nk_mod {
    use ql_ad::Real;

    /// Input specification for the Christiano, Eichenbaum & Evans (2005, *Journal of Political
    /// Economy* 113(1)) & Smets & Wouters (2007, *AER* 97(3)) medium-scale New Keynesian DSGE
    /// model with hybrid Calvo price stickiness, external habit formation, and investment-growth
    /// adjustment costs.
    #[derive(Debug, Clone, Copy)]
    pub struct ChristianoEichenbaumEvansNkSpec<T: Real> {
        /// Calvo price non-reoptimization probability $\xi_p \in (0, 1)$ (e.g., $0.68$).
        pub calvo_price_stickiness_xi_p: T,
        /// Backward inflation indexation weight $\iota_p \in (0, 1)$ (e.g., $0.45$).
        pub inflation_indexation_iota_p: T,
        /// Household subjective discount factor $\beta \in (0, 1)$ (e.g., $0.99$).
        pub discount_factor_beta: T,
        /// External consumption habit persistence $h \in (0, 1)$ (e.g., $0.70$).
        pub habit_persistence_h: T,
        /// Investment growth-rate adjustment cost curvature $S''(1) = \kappa_I > 0$ (e.g., $4.50$).
        pub investment_adjustment_curvature_kappa_i: T,
        /// Steady-state investment share of output $\omega_I \in (0, 1)$ (e.g., $0.22$).
        pub steady_state_investment_share: T,
        /// Real interest-rate gap shock $\hat{r}_t > 0$ (e.g., $0.010$).
        pub real_interest_rate_shock: T,
    }

    /// Analytics returned by [`evaluate_christiano_eichenbaum_evans_calvo_investment_nk_impulse`].
    #[derive(Debug, Clone, Copy)]
    pub struct ChristianoEichenbaumEvansNkAnalytics<T: Real> {
        /// Structural hybrid New Keynesian Phillips Curve (NKPC) marginal-cost slope $\kappa_p$.
        pub hybrid_nkpc_marginal_cost_slope: T,
        /// Habit-attenuated consumption Euler interest-rate sensitivity $-\partial \hat{c}_t / \partial \hat{r}_t$.
        pub habit_euler_consumption_sensitivity: T,
        /// Investment-adjustment-cost Tobin's $Q$ sensitivity $-\partial \hat{I}_t / \partial \hat{r}_t$.
        pub investment_adjustment_q_sensitivity: T,
        /// Equilibrium inflation impact response $-\hat{\pi}_t$ to a contractionary real interest rate shock $\hat{r}_t$.
        pub contractionary_inflation_impulse_response: T,
    }

    /// Evaluates the Christiano, Eichenbaum & Evans (2005) / Smets & Wouters (2007) structural
    /// New Keynesian transmission mechanism:
    /// - Hybrid NKPC slope $\kappa_p = \frac{(1 - \xi_p)(1 - \beta \xi_p)}{\xi_p(1 + \beta \iota_p)}$,
    /// - Habit-adjusted consumption Euler sensitivity $\eta_c = \frac{1 - h}{1 + h}$,
    /// - Investment flow adjustment-cost sensitivity $\eta_I = \frac{1}{\kappa_I(1 + \beta)}$,
    /// - Aggregate demand contraction $\Delta \hat{y}_t = -\big((1 - \omega_I)\eta_c + \omega_I \eta_I\big)\hat{r}_t$, and
    /// - Hybrid NKPC inflation impact response $-\hat{\pi}_t = \kappa_p \left(\frac{1 + h}{1 - h}\right)(-\Delta \hat{y}_t)$.
    pub fn evaluate_christiano_eichenbaum_evans_calvo_investment_nk_impulse<T: Real>(
        spec: &ChristianoEichenbaumEvansNkSpec<T>,
    ) -> ChristianoEichenbaumEvansNkAnalytics<T> {
        let xi_p = spec.calvo_price_stickiness_xi_p;
        let iota_p = spec.inflation_indexation_iota_p;
        let beta = spec.discount_factor_beta;
        let h = spec.habit_persistence_h;
        let kappa_i = spec.investment_adjustment_curvature_kappa_i;
        let omega_i = spec.steady_state_investment_share;
        let r_shock = spec.real_interest_rate_shock;

        let kappa_p =
            ((T::one() - xi_p) * (T::one() - beta * xi_p)) / (xi_p * (T::one() + beta * iota_p));
        let eta_c = (T::one() - h) / (T::one() + h);
        let eta_i = T::one() / (kappa_i * (T::one() + beta));

        let output_contraction = ((T::one() - omega_i) * eta_c + omega_i * eta_i) * r_shock;
        let mc_elasticity = T::one() + T::from_f64(0.5) * (T::one() / eta_c + kappa_i * omega_i);
        let inflation_drop = kappa_p * mc_elasticity * output_contraction;

        ChristianoEichenbaumEvansNkAnalytics {
            hybrid_nkpc_marginal_cost_slope: kappa_p,
            habit_euler_consumption_sensitivity: eta_c,
            investment_adjustment_q_sensitivity: eta_i,
            contractionary_inflation_impulse_response: inflation_drop,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn christiano_eichenbaum_evans_nk_invariants_and_ad_parity() {
            let base = ChristianoEichenbaumEvansNkSpec {
                calvo_price_stickiness_xi_p: 0.68_f64,
                inflation_indexation_iota_p: 0.45_f64,
                discount_factor_beta: 0.99_f64,
                habit_persistence_h: 0.70_f64,
                investment_adjustment_curvature_kappa_i: 4.50_f64,
                steady_state_investment_share: 0.22_f64,
                real_interest_rate_shock: 0.010_f64,
            };
            let res = evaluate_christiano_eichenbaum_evans_calvo_investment_nk_impulse(&base);
            assert!(res.hybrid_nkpc_marginal_cost_slope > 0.0);
            assert!(res.habit_euler_consumption_sensitivity > 0.0);
            assert!(res.investment_adjustment_q_sensitivity > 0.0);
            assert!(res.contractionary_inflation_impulse_response > 0.0);

            let inputs = [
                base.calvo_price_stickiness_xi_p,
                base.inflation_indexation_iota_p,
                base.discount_factor_beta,
                base.habit_persistence_h,
                base.investment_adjustment_curvature_kappa_i,
                base.steady_state_investment_share,
                base.real_interest_rate_shock,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = ChristianoEichenbaumEvansNkSpec {
                calvo_price_stickiness_xi_p: dual_in[0],
                inflation_indexation_iota_p: dual_in[1],
                discount_factor_beta: dual_in[2],
                habit_persistence_h: dual_in[3],
                investment_adjustment_curvature_kappa_i: dual_in[4],
                steady_state_investment_share: dual_in[5],
                real_interest_rate_shock: dual_in[6],
            };
            let r_ad = evaluate_christiano_eichenbaum_evans_calvo_investment_nk_impulse(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| ChristianoEichenbaumEvansNkSpec {
                    calvo_price_stickiness_xi_p: v[0],
                    inflation_indexation_iota_p: v[1],
                    discount_factor_beta: v[2],
                    habit_persistence_h: v[3],
                    investment_adjustment_curvature_kappa_i: v[4],
                    steady_state_investment_share: v[5],
                    real_interest_rate_shock: v[6],
                };
                let rp = evaluate_christiano_eichenbaum_evans_calvo_investment_nk_impulse(&mk(p));
                let rm = evaluate_christiano_eichenbaum_evans_calvo_investment_nk_impulse(&mk(m));
                let fd = (rp.contractionary_inflation_impulse_response
                    - rm.contractionary_inflation_impulse_response)
                    / (2.0 * eps);
                assert!(
                    r_ad.contractionary_inflation_impulse_response.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.contractionary_inflation_impulse_response.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use christiano_eichenbaum_evans_nk_mod::*;
mod brunnermeier_sannikov_paradox_mod {
    use ql_ad::Real;

    /// Input specification for the Brunnermeier & Sannikov (2014, *American Economic Review*
    /// 104(2), *"A Macroeconomic Model with a Financial Sector"*) intermediary net-worth share
    /// $\eta_t$, endogenous fire-sale amplification, and Volatility Paradox engine.
    #[derive(Debug, Clone, Copy)]
    pub struct BrunnermeierSannikovVolatilityParadoxSpec<T: Real> {
        /// Expert / financial intermediary wealth share $\eta_t \in (0, 1)$ (e.g., $0.24$).
        pub intermediary_wealth_share_eta: T,
        /// Fundamental productivity Brownian volatility $\sigma > 0$ (e.g., $0.025$).
        pub fundamental_capital_volatility_sigma: T,
        /// Expert productivity advantage $a_e - a_h > 0$ over less-productive households (e.g., $0.055$).
        pub expert_productivity_wedge_delta_a: T,
        /// Capital illiquidity / fire-sale price elasticity $\phi_{\text{ill}} > 0$ (e.g., $2.20$).
        pub fire_sale_illiquidity_sensitivity_phi: T,
        /// Household discount rate $\rho_h > 0$ (e.g., $0.050$).
        pub household_discount_rate_rho: T,
        /// Expert balance-sheet leverage ratio $L_{\text{exp}} = \psi_t / \eta_t > 1$ (e.g., $3.20$).
        pub intermediary_leverage_ratio: T,
    }

    /// Analytics returned by [`evaluate_brunnermeier_sannikov_intermediary_volatility_paradox`].
    #[derive(Debug, Clone, Copy)]
    pub struct BrunnermeierSannikovVolatilityParadoxAnalytics<T: Real> {
        /// Equilibrium capital price $q(\eta_t)$ interpolating between fire-sale floor and expert valuation.
        pub equilibrium_capital_price_q: T,
        /// Endogenous amplification volatility $\sigma^q(\eta_t)$ of capital price.
        pub endogenous_price_volatility_sigma_q: T,
        /// Total systemic return volatility $\sigma_{\text{tot}}(\eta_t) = \sigma + \sigma^q(\eta_t)$ (amplified by $1 / (1 - \text{loop})$).
        pub amplified_total_return_volatility: T,
        /// Volatility Paradox crisis tail-risk index (showing higher steady-state leverage and crisis volatility as $\sigma \downarrow$).
        pub volatility_paradox_amplification_ratio: T,
    }

    /// Evaluates the Brunnermeier & Sannikov (2014, Eqs. 10–19) amplification loop:
    /// $$\sigma + \sigma^q(\eta_t) = \frac{\sigma}{1 - \frac{q'(\eta_t)}{q(\eta_t)}\eta_t(L_{\text{exp}} - 1)},$$
    /// and the Volatility Paradox whereby a compression in exogenous fundamental volatility $\sigma$
    /// induces higher equilibrium leverage and leaves endogenous crisis volatility $\sigma^q$ elevated.
    pub fn evaluate_brunnermeier_sannikov_intermediary_volatility_paradox<T: Real>(
        spec: &BrunnermeierSannikovVolatilityParadoxSpec<T>,
    ) -> BrunnermeierSannikovVolatilityParadoxAnalytics<T> {
        let eta = spec.intermediary_wealth_share_eta;
        let sig = spec.fundamental_capital_volatility_sigma;
        let da = spec.expert_productivity_wedge_delta_a;
        let phi = spec.fire_sale_illiquidity_sensitivity_phi;
        let rho = spec.household_discount_rate_rho;
        let lev = spec.intermediary_leverage_ratio;

        // Concave capital price q(eta) rising from household value toward expert value
        let q_min = T::one() - da / (rho + phi * sig);
        let q_slope = (da / rho) * (-phi * eta).exp();
        let q_eta = q_min + (da / (rho * phi)) * (T::one() - (-phi * eta).exp());

        // Brunnermeier-Sannikov (2014, Eq. 14) amplification denominator 1 - (q'(eta)/q(eta)) * eta * (lev - 1)
        let loop_gain = (q_slope / q_eta) * eta * (lev - T::one());
        let denom_amp = T::one() - loop_gain;
        let sig_tot = sig / denom_amp;
        let sig_q = sig_tot - sig;

        // Volatility Paradox ratio: endogenous amplification ratio scaled by inverse fundamental volatility cushion
        let paradox_ratio = (sig_tot / sig) * (T::one() + da / (rho + sig));

        BrunnermeierSannikovVolatilityParadoxAnalytics {
            equilibrium_capital_price_q: q_eta,
            endogenous_price_volatility_sigma_q: sig_q,
            amplified_total_return_volatility: sig_tot,
            volatility_paradox_amplification_ratio: paradox_ratio,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn brunnermeier_sannikov_paradox_invariants_and_ad_parity() {
            let base = BrunnermeierSannikovVolatilityParadoxSpec {
                intermediary_wealth_share_eta: 0.24_f64,
                fundamental_capital_volatility_sigma: 0.025_f64,
                expert_productivity_wedge_delta_a: 0.055_f64,
                fire_sale_illiquidity_sensitivity_phi: 2.20_f64,
                household_discount_rate_rho: 0.050_f64,
                intermediary_leverage_ratio: 3.20_f64,
            };
            let res = evaluate_brunnermeier_sannikov_intermediary_volatility_paradox(&base);
            assert!(res.equilibrium_capital_price_q > 0.5);
            assert!(res.endogenous_price_volatility_sigma_q > 0.0);
            assert!(
                res.amplified_total_return_volatility > base.fundamental_capital_volatility_sigma
            );
            assert!(res.volatility_paradox_amplification_ratio > 1.5);

            let inputs = [
                base.intermediary_wealth_share_eta,
                base.fundamental_capital_volatility_sigma,
                base.expert_productivity_wedge_delta_a,
                base.fire_sale_illiquidity_sensitivity_phi,
                base.household_discount_rate_rho,
                base.intermediary_leverage_ratio,
            ];
            let mut dual_in = [Dual::<6>::from_value(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BrunnermeierSannikovVolatilityParadoxSpec {
                intermediary_wealth_share_eta: dual_in[0],
                fundamental_capital_volatility_sigma: dual_in[1],
                expert_productivity_wedge_delta_a: dual_in[2],
                fire_sale_illiquidity_sensitivity_phi: dual_in[3],
                household_discount_rate_rho: dual_in[4],
                intermediary_leverage_ratio: dual_in[5],
            };
            let r_ad = evaluate_brunnermeier_sannikov_intermediary_volatility_paradox(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| BrunnermeierSannikovVolatilityParadoxSpec {
                    intermediary_wealth_share_eta: v[0],
                    fundamental_capital_volatility_sigma: v[1],
                    expert_productivity_wedge_delta_a: v[2],
                    fire_sale_illiquidity_sensitivity_phi: v[3],
                    household_discount_rate_rho: v[4],
                    intermediary_leverage_ratio: v[5],
                };
                let rp = evaluate_brunnermeier_sannikov_intermediary_volatility_paradox(&mk(p));
                let rm = evaluate_brunnermeier_sannikov_intermediary_volatility_paradox(&mk(m));
                let fd = (rp.volatility_paradox_amplification_ratio
                    - rm.volatility_paradox_amplification_ratio)
                    / (2.0 * eps);
                assert!(
                    r_ad.volatility_paradox_amplification_ratio.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.volatility_paradox_amplification_ratio.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use brunnermeier_sannikov_paradox_mod::*;
mod pastor_veronesi_political_uncertainty_mod {
    use ql_ad::Real;

    /// Input specification for the Pástor & Veronesi (2012, *Journal of Finance* 67(4);
    /// 2013, *Journal of Financial Economics* 110(3), *"Political Uncertainty and Risk Premia"*)
    /// Bayesian government policy learning and political uncertainty equity risk premium engine.
    #[derive(Debug, Clone, Copy)]
    pub struct PastorVeronesiPoliticalUncertaintySpec<T: Real> {
        /// Posterior mean perception $\hat{g}_t$ of current government policy impact on profitability (e.g., $-0.012$).
        pub perceived_policy_impact_g_hat: T,
        /// Prior standard deviation $\sigma_{g,0} > 0$ of unknown policy impact (e.g., $0.035$).
        pub prior_policy_uncertainty_sigma_g0: T,
        /// Fundamental macroeconomic shock volatility $\sigma_{\text{macro}} > 0$ (e.g., $0.14$).
        pub fundamental_macro_volatility_sigma_m: T,
        /// Political cost heterogeneity dispersion $\sigma_{\text{pol}} > 0$ across replacement policies (e.g., $0.085$).
        pub political_cost_dispersion_sigma_pol: T,
        /// Representative investor CRRA risk aversion $\gamma > 1$ (e.g., $3.50$).
        pub investor_risk_aversion_gamma: T,
        /// Elapsed Bayesian learning duration $\tau_{\text{learn}} > 0$ in years (e.g., $2.0$).
        pub bayesian_learning_duration_tau: T,
    }

    /// Analytics returned by [`evaluate_pastor_veronesi_political_uncertainty_risk_premium`].
    #[derive(Debug, Clone, Copy)]
    pub struct PastorVeronesiPoliticalUncertaintyAnalytics<T: Real> {
        /// Bayesian posterior policy uncertainty variance $\hat{\sigma}_g^2(\tau_{\text{learn}})$.
        pub posterior_policy_uncertainty_variance: T,
        /// Endogenous probability $p_{\text{change}} \in (0, 1)$ of a government policy regime switch.
        pub policy_regime_switch_probability: T,
        /// Total equity return volatility $\sigma_R$ incorporating macro, Bayesian learning, and political jump risk.
        pub total_equity_return_volatility: T,
        /// Equilibrium equity risk premium $\pi_{\text{EQ}}$ driven by economic and political uncertainty.
        pub political_uncertainty_equity_risk_premium: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        T::from_f64(0.5) * (T::one() + (x * T::from_f64(0.707_106_781_186_547_5)).erf())
    }

    /// Evaluates the Pástor & Veronesi (2012, 2013, Propositions 1–4) Bayesian policy-learning
    /// precision update $\hat{\sigma}_g^{-2} = \sigma_{g,0}^{-2} + \tau_{\text{learn}}\sigma_{\text{macro}}^{-2}$,
    /// policy adoption probability $p_{\text{change}} = \Phi\!\left(\frac{-\hat{g}_t}{\sqrt{\hat{\sigma}_g^2 + \sigma_{\text{pol}}^2}}\right)$,
    /// and equilibrium political uncertainty equity risk premium.
    pub fn evaluate_pastor_veronesi_political_uncertainty_risk_premium<T: Real>(
        spec: &PastorVeronesiPoliticalUncertaintySpec<T>,
    ) -> PastorVeronesiPoliticalUncertaintyAnalytics<T> {
        let g_hat = spec.perceived_policy_impact_g_hat;
        let sig_g0 = spec.prior_policy_uncertainty_sigma_g0;
        let sig_m = spec.fundamental_macro_volatility_sigma_m;
        let sig_pol = spec.political_cost_dispersion_sigma_pol;
        let gam = spec.investor_risk_aversion_gamma;
        let tau = spec.bayesian_learning_duration_tau;

        let prec_post = T::one() / (sig_g0 * sig_g0) + tau / (sig_m * sig_m);
        let var_g_post = T::one() / prec_post;

        let switch_std = (var_g_post + sig_pol * sig_pol).sqrt();
        let p_change = std_norm_cdf(-g_hat / switch_std);

        let var_eq = sig_m * sig_m
            + gam * gam * var_g_post
            + p_change * (T::one() - p_change) * sig_pol * sig_pol;
        let vol_eq = var_eq.sqrt();

        let erp = gam * (sig_m * sig_m + var_g_post) + gam * p_change * sig_pol * sig_pol;

        PastorVeronesiPoliticalUncertaintyAnalytics {
            posterior_policy_uncertainty_variance: var_g_post,
            policy_regime_switch_probability: p_change,
            total_equity_return_volatility: vol_eq,
            political_uncertainty_equity_risk_premium: erp,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn pastor_veronesi_political_uncertainty_invariants_and_ad_parity() {
            let base = PastorVeronesiPoliticalUncertaintySpec {
                perceived_policy_impact_g_hat: -0.012_f64,
                prior_policy_uncertainty_sigma_g0: 0.035_f64,
                fundamental_macro_volatility_sigma_m: 0.14_f64,
                political_cost_dispersion_sigma_pol: 0.085_f64,
                investor_risk_aversion_gamma: 3.50_f64,
                bayesian_learning_duration_tau: 2.0_f64,
            };
            let res = evaluate_pastor_veronesi_political_uncertainty_risk_premium(&base);
            assert!(
                res.posterior_policy_uncertainty_variance
                    < base.prior_policy_uncertainty_sigma_g0
                        * base.prior_policy_uncertainty_sigma_g0
            );
            assert!(res.policy_regime_switch_probability > 0.5);
            assert!(res.total_equity_return_volatility > base.fundamental_macro_volatility_sigma_m);
            assert!(res.political_uncertainty_equity_risk_premium > 0.05);

            let inputs = [
                base.perceived_policy_impact_g_hat,
                base.prior_policy_uncertainty_sigma_g0,
                base.fundamental_macro_volatility_sigma_m,
                base.political_cost_dispersion_sigma_pol,
                base.investor_risk_aversion_gamma,
                base.bayesian_learning_duration_tau,
            ];
            let mut dual_in = [Dual::<6>::from_value(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = PastorVeronesiPoliticalUncertaintySpec {
                perceived_policy_impact_g_hat: dual_in[0],
                prior_policy_uncertainty_sigma_g0: dual_in[1],
                fundamental_macro_volatility_sigma_m: dual_in[2],
                political_cost_dispersion_sigma_pol: dual_in[3],
                investor_risk_aversion_gamma: dual_in[4],
                bayesian_learning_duration_tau: dual_in[5],
            };
            let r_ad = evaluate_pastor_veronesi_political_uncertainty_risk_premium(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| PastorVeronesiPoliticalUncertaintySpec {
                    perceived_policy_impact_g_hat: v[0],
                    prior_policy_uncertainty_sigma_g0: v[1],
                    fundamental_macro_volatility_sigma_m: v[2],
                    political_cost_dispersion_sigma_pol: v[3],
                    investor_risk_aversion_gamma: v[4],
                    bayesian_learning_duration_tau: v[5],
                };
                let rp = evaluate_pastor_veronesi_political_uncertainty_risk_premium(&mk(p));
                let rm = evaluate_pastor_veronesi_political_uncertainty_risk_premium(&mk(m));
                let fd = (rp.political_uncertainty_equity_risk_premium
                    - rm.political_uncertainty_equity_risk_premium)
                    / (2.0 * eps);
                assert!(
                    r_ad.political_uncertainty_equity_risk_premium.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.political_uncertainty_equity_risk_premium.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use pastor_veronesi_political_uncertainty_mod::*;

mod lettau_ludvigson_cay_mod {
    use ql_ad::Real;

    /// Input specification for the Lettau & Ludvigson (2001, *Journal of Finance* 56(3),
    /// *"Consumption, Aggregate Wealth, and Expected Stock Returns"*; 2001b, *JPE* 109(6))
    /// cointegrating consumption–wealth ratio ($\text{cay}_t$) conditional CCAPM engine.
    #[derive(Debug, Clone, Copy)]
    pub struct LettauLudvigsonCayConditionalRiskPremiumSpec<T: Real> {
        /// Standardized log aggregate consumption $c_t$ (e.g., $10.24$).
        pub log_consumption_c: T,
        /// Standardized log household financial asset wealth $a_t$ (e.g., $10.05$).
        pub log_asset_wealth_a: T,
        /// Standardized log labor income $y_t$ (e.g., $10.28$).
        pub log_labor_income_y: T,
        /// Cointegrating asset-wealth share $\omega_a \in (0, 1)$ (e.g., $0.32$).
        pub asset_wealth_cointegration_share_omega: T,
        /// Baseline unconditional relative risk aversion $\gamma_0 > 0$ (e.g., $12.0$).
        pub baseline_risk_aversion_gamma0: T,
        /// Conditioning sensitivity $\gamma_1 > 0$ of risk aversion to $\text{cay}_t$ (e.g., $185.0$).
        pub cay_conditioning_sensitivity_gamma1: T,
        /// Annualized aggregate consumption growth volatility $\sigma_c > 0$ (e.g., $0.018$).
        pub consumption_volatility_sigma_c: T,
        /// Annualized equity market return volatility $\sigma_m > 0$ (e.g., $0.165$).
        pub equity_return_volatility_sigma_m: T,
    }

    /// Analytics returned by [`evaluate_lettau_ludvigson_cay_conditional_equity_premium`].
    #[derive(Debug, Clone, Copy)]
    pub struct LettauLudvigsonCayConditionalRiskPremiumAnalytics<T: Real> {
        /// Cointegrating consumption-wealth residual $\text{cay}_t = c_t - \omega_a a_t - (1 - \omega_a) y_t$.
        pub cointegrating_residual_cay: T,
        /// Time-varying effective relative risk aversion $\gamma_{\text{eff}}(\text{cay}_t) = \gamma_0 + \gamma_1\text{cay}_t$.
        pub conditional_effective_risk_aversion: T,
        /// Conditional expected equity risk premium $\mathbb{E}_t[R_{m,t+1} - R_{f,t}]$.
        pub conditional_equity_risk_premium: T,
        /// Conditional equity market Sharpe ratio $SR_t(\text{cay}_t) = \mathbb{E}_t[R_{m,t+1} - R_{f,t}] / \sigma_m$.
        pub conditional_sharpe_ratio: T,
    }

    /// Evaluates the Lettau & Ludvigson (2001, Eq. 6; 2001b, Eq. 8) cointegrating residual
    /// $\text{cay}_t = c_t - \omega_a a_t - (1 - \omega_a) y_t$, conditional effective risk
    /// aversion $\gamma_{\text{eff}}(\text{cay}_t) = \gamma_0 + \gamma_1\text{cay}_t$, and
    /// conditional CCAPM equity risk premium.
    pub fn evaluate_lettau_ludvigson_cay_conditional_equity_premium<T: Real>(
        spec: &LettauLudvigsonCayConditionalRiskPremiumSpec<T>,
    ) -> LettauLudvigsonCayConditionalRiskPremiumAnalytics<T> {
        let c_t = spec.log_consumption_c;
        let a_t = spec.log_asset_wealth_a;
        let y_t = spec.log_labor_income_y;
        let omega = spec.asset_wealth_cointegration_share_omega;
        let gam0 = spec.baseline_risk_aversion_gamma0;
        let gam1 = spec.cay_conditioning_sensitivity_gamma1;
        let sig_c = spec.consumption_volatility_sigma_c;
        let sig_m = spec.equity_return_volatility_sigma_m;

        let cay = c_t - omega * a_t - (T::one() - omega) * y_t;
        let gam_eff = gam0 + gam1 * cay;
        let rho_cm = T::from_f64(0.65);
        let erp = gam_eff * rho_cm * sig_c * sig_m + T::from_f64(0.5) * sig_m * sig_m * cay;
        let sr = erp / sig_m;

        LettauLudvigsonCayConditionalRiskPremiumAnalytics {
            cointegrating_residual_cay: cay,
            conditional_effective_risk_aversion: gam_eff,
            conditional_equity_risk_premium: erp,
            conditional_sharpe_ratio: sr,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn lettau_ludvigson_cay_invariants_and_ad_parity() {
            let base = LettauLudvigsonCayConditionalRiskPremiumSpec {
                log_consumption_c: 10.24_f64,
                log_asset_wealth_a: 10.05_f64,
                log_labor_income_y: 10.28_f64,
                asset_wealth_cointegration_share_omega: 0.32_f64,
                baseline_risk_aversion_gamma0: 12.0_f64,
                cay_conditioning_sensitivity_gamma1: 185.0_f64,
                consumption_volatility_sigma_c: 0.018_f64,
                equity_return_volatility_sigma_m: 0.165_f64,
            };
            let res = evaluate_lettau_ludvigson_cay_conditional_equity_premium(&base);
            assert!(res.cointegrating_residual_cay > 0.0);
            assert!(res.conditional_effective_risk_aversion > base.baseline_risk_aversion_gamma0);
            assert!(res.conditional_equity_risk_premium > 0.03);
            assert!(res.conditional_sharpe_ratio > 0.15);

            let inputs = [
                base.log_consumption_c,
                base.log_asset_wealth_a,
                base.log_labor_income_y,
                base.asset_wealth_cointegration_share_omega,
                base.baseline_risk_aversion_gamma0,
                base.cay_conditioning_sensitivity_gamma1,
                base.consumption_volatility_sigma_c,
                base.equity_return_volatility_sigma_m,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = LettauLudvigsonCayConditionalRiskPremiumSpec {
                log_consumption_c: dual_in[0],
                log_asset_wealth_a: dual_in[1],
                log_labor_income_y: dual_in[2],
                asset_wealth_cointegration_share_omega: dual_in[3],
                baseline_risk_aversion_gamma0: dual_in[4],
                cay_conditioning_sensitivity_gamma1: dual_in[5],
                consumption_volatility_sigma_c: dual_in[6],
                equity_return_volatility_sigma_m: dual_in[7],
            };
            let ad = evaluate_lettau_ludvigson_cay_conditional_equity_premium(&spec_ad)
                .conditional_equity_risk_premium;
            let eps = 1e-7;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let f_p = evaluate_lettau_ludvigson_cay_conditional_equity_premium(
                    &LettauLudvigsonCayConditionalRiskPremiumSpec {
                        log_consumption_c: p[0],
                        log_asset_wealth_a: p[1],
                        log_labor_income_y: p[2],
                        asset_wealth_cointegration_share_omega: p[3],
                        baseline_risk_aversion_gamma0: p[4],
                        cay_conditioning_sensitivity_gamma1: p[5],
                        consumption_volatility_sigma_c: p[6],
                        equity_return_volatility_sigma_m: p[7],
                    },
                )
                .conditional_equity_risk_premium;
                let f_m = evaluate_lettau_ludvigson_cay_conditional_equity_premium(
                    &LettauLudvigsonCayConditionalRiskPremiumSpec {
                        log_consumption_c: m[0],
                        log_asset_wealth_a: m[1],
                        log_labor_income_y: m[2],
                        asset_wealth_cointegration_share_omega: m[3],
                        baseline_risk_aversion_gamma0: m[4],
                        cay_conditioning_sensitivity_gamma1: m[5],
                        consumption_volatility_sigma_c: m[6],
                        equity_return_volatility_sigma_m: m[7],
                    },
                )
                .conditional_equity_risk_premium;
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
pub use lettau_ludvigson_cay_mod::{
    LettauLudvigsonCayConditionalRiskPremiumAnalytics,
    LettauLudvigsonCayConditionalRiskPremiumSpec,
    evaluate_lettau_ludvigson_cay_conditional_equity_premium,
};
