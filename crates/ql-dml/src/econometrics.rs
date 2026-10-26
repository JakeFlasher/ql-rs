mod evaluate_wooldridge_greene_2sls_iv_and_did_module {
    use ql_ad::Real;
    #[derive(Clone, Debug)]
    pub struct WooldridgeGreene2slsDidSpec<T: Real = f64> {
        pub endogenous_x: Vec<T>,
        pub instrument_z: Vec<T>,
        pub outcome_y: Vec<T>,
        pub treated_pre_mean: T,
        pub treated_post_mean: T,
        pub control_pre_mean: T,
        pub control_post_mean: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct WooldridgeGreene2slsDidReport<T: Real = f64> {
        pub ols_slope_beta: T,
        pub first_stage_instrument_pi1: T,
        pub first_stage_f_statistic: T,
        pub two_stage_ls_iv_beta: T,
        pub hausman_endogeneity_wedge: T,
        pub control_function_residual_coef: T,
        pub did_parallel_trend_counterfactual: T,
        pub did_treatment_effect_att: T,
    }
    pub fn evaluate_wooldridge_greene_2sls_iv_and_did<T: Real>(
        s: &WooldridgeGreene2slsDidSpec<T>,
    ) -> WooldridgeGreene2slsDidReport<T> {
        let n = s.endogenous_x.len();
        let mx = s.endogenous_x.iter().copied().fold(T::zero(), |a, x| a + x) / (n as f64);
        let mz = s.instrument_z.iter().copied().fold(T::zero(), |a, x| a + x) / (n as f64);
        let my = s.outcome_y.iter().copied().fold(T::zero(), |a, x| a + x) / (n as f64);
        let mut xx = T::zero();
        let mut zz = T::zero();
        let mut xz = T::zero();
        let mut xy = T::zero();
        let mut zy = T::zero();
        for i in 0..n {
            let x = s.endogenous_x[i] - mx;
            let z = s.instrument_z[i] - mz;
            let y = s.outcome_y[i] - my;
            xx += x * x;
            zz += z * z;
            xz += x * z;
            xy += x * y;
            zy += z * y;
        }
        let ols = xy / xx;
        let pi = xz / zz;
        let mut rss = T::zero();
        let mut vu = T::zero();
        let iv = zy / xz;
        for i in 0..n {
            let fit = mx + pi * (s.instrument_z[i] - mz);
            let vhat = s.endogenous_x[i] - fit;
            let uhat = (s.outcome_y[i] - my) - iv * (s.endogenous_x[i] - mx);
            rss += vhat * vhat;
            vu += vhat * uhat;
        }
        let rss_safe = rss.max_val(T::from_f64(1e-12));
        let f = T::from_f64((n.saturating_sub(2).max(1)) as f64) * pi * pi * zz / rss_safe;
        let cf = vu / rss_safe;
        let counter = s.treated_pre_mean + s.control_post_mean - s.control_pre_mean;
        let att = s.treated_post_mean - counter;
        WooldridgeGreene2slsDidReport {
            ols_slope_beta: ols,
            first_stage_instrument_pi1: pi,
            first_stage_f_statistic: f,
            two_stage_ls_iv_beta: iv,
            hausman_endogeneity_wedge: ols - iv,
            control_function_residual_coef: cf,
            did_parallel_trend_counterfactual: counter,
            did_treatment_effect_att: att,
        }
    }
    #[cfg(test)]
    mod domain_suite_45_tests {
        use super::*;
        #[test]
        fn suite45_dml_smoke() {
            let s = WooldridgeGreene2slsDidSpec {
                endogenous_x: vec![1.0, 2.0, 3.0, 4.0],
                instrument_z: vec![1.0, 2.1, 2.9, 4.1],
                outcome_y: vec![1.2, 2.1, 3.2, 4.2],
                treated_pre_mean: 1.0,
                treated_post_mean: 1.2,
                control_pre_mean: 1.0,
                control_post_mean: 1.1,
            };
            let r = evaluate_wooldridge_greene_2sls_iv_and_did(&s);
            assert!(r.first_stage_f_statistic > 10.0);
            assert!(
                (r.did_treatment_effect_att
                    - (s.treated_post_mean - r.did_parallel_trend_counterfactual))
                    .abs()
                    < 1e-12
            );
        }
    }
    #[cfg(test)]
    mod domain_suite_45_tests_ad {
        use super::*;
        use ql_ad::Dual;
        fn s<T: Real>(x: T) -> WooldridgeGreene2slsDidSpec<T> {
            WooldridgeGreene2slsDidSpec {
                endogenous_x: vec![
                    T::one() + x,
                    T::from_f64(2.),
                    T::from_f64(3.),
                    T::from_f64(4.),
                ],
                instrument_z: vec![
                    T::one(),
                    T::from_f64(2.1),
                    T::from_f64(2.9),
                    T::from_f64(4.1),
                ],
                outcome_y: vec![
                    T::from_f64(1.2),
                    T::from_f64(2.1),
                    T::from_f64(3.2),
                    T::from_f64(4.2),
                ],
                treated_pre_mean: T::one(),
                treated_post_mean: T::from_f64(1.2),
                control_pre_mean: T::one(),
                control_post_mean: T::from_f64(1.1),
            }
        }
        #[test]
        fn dual_fd() {
            let h = 1e-6;
            let fd = (evaluate_wooldridge_greene_2sls_iv_and_did(&s(h)).two_stage_ls_iv_beta
                - evaluate_wooldridge_greene_2sls_iv_and_did(&s(-h)).two_stage_ls_iv_beta)
                / (2. * h);
            let ad = evaluate_wooldridge_greene_2sls_iv_and_did(&s(Dual::<1>::variable(0., 0)))
                .two_stage_ls_iv_beta;
            assert!((ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_wooldridge_greene_2sls_iv_and_did_module::*;

mod evaluate_bayer_path_signature_chen_and_levy_area_module {
    use ql_ad::Real;
    #[derive(Clone, Debug)]
    pub struct BayerPathSignatureChenLevyAreaSpec<T: Real = f64> {
        pub time_grid: Vec<T>,
        pub spot_path: Vec<T>,
        pub linear_functional_weights: [T; 6],
    }
    #[derive(Clone, Copy, Debug)]
    pub struct BayerPathSignatureChenLevyAreaReport<T: Real = f64> {
        pub sig1_time_increment: T,
        pub sig1_spot_increment: T,
        pub sig2_time_spot_12: T,
        pub sig2_spot_time_21: T,
        pub sig2_spot_spot_22: T,
        pub signed_levy_area: T,
        pub chen_identity_l1_error: T,
        pub shuffle_product_identity_error: T,
        pub signature_linear_payoff: T,
    }
    fn sig_segment_range<T: Real>(
        t: &[T],
        x: &[T],
        start: usize,
        end: usize,
    ) -> (T, T, T, T, T, T) {
        let mut s1 = T::zero();
        let mut s2 = T::zero();
        let mut s11 = T::zero();
        let mut s12 = T::zero();
        let mut s21 = T::zero();
        let mut s22 = T::zero();
        let half = T::from_f64(0.5);
        for k in start..end {
            let dt = t[k + 1] - t[k];
            let dx = x[k + 1] - x[k];
            s11 += s1 * dt + half * dt * dt;
            s12 += s1 * dx + half * dt * dx;
            s21 += s2 * dt + half * dx * dt;
            s22 += s2 * dx + half * dx * dx;
            s1 += dt;
            s2 += dx;
        }
        (s1, s2, s11, s12, s21, s22)
    }
    pub fn evaluate_bayer_path_signature_chen_and_levy_area<T: Real>(
        s: &BayerPathSignatureChenLevyAreaSpec<T>,
    ) -> BayerPathSignatureChenLevyAreaReport<T> {
        let m = s.time_grid.len().min(s.spot_path.len()).saturating_sub(1);
        let (s1, s2, s11, s12, s21, s22) = sig_segment_range(&s.time_grid, &s.spot_path, 0, m);
        let mid = m / 2;
        let (a1, a2, a11, a12, a21, a22) = sig_segment_range(&s.time_grid, &s.spot_path, 0, mid);
        let (b1, b2, b11, b12, b21, b22) = sig_segment_range(&s.time_grid, &s.spot_path, mid, m);
        let c11 = s11 - (a11 + b11 + a1 * b1);
        let c12 = s12 - (a12 + b12 + a1 * b2);
        let c21 = s21 - (a21 + b21 + a2 * b1);
        let c22 = s22 - (a22 + b22 + a2 * b2);
        let chen_err = (c11 * c11 + c12 * c12 + c21 * c21 + c22 * c22 + T::from_f64(1e-30)).sqrt();
        let shuf = (s12 + s21) - s1 * s2;
        let shuf_err = (shuf * shuf + T::from_f64(1e-30)).sqrt();
        let levy = T::from_f64(0.5) * (s12 - s21);
        let w = &s.linear_functional_weights;
        let payoff = w[0] * s1 + w[1] * s2 + w[2] * s11 + w[3] * s12 + w[4] * s21 + w[5] * s22;
        BayerPathSignatureChenLevyAreaReport {
            sig1_time_increment: s1,
            sig1_spot_increment: s2,
            sig2_time_spot_12: s12,
            sig2_spot_time_21: s21,
            sig2_spot_spot_22: s22,
            signed_levy_area: levy,
            chen_identity_l1_error: chen_err,
            shuffle_product_identity_error: shuf_err,
            signature_linear_payoff: payoff,
        }
    }

    #[cfg(test)]
    mod domain_suite_46_tests {
        use super::*;
        use ql_ad::Dual;
        fn spec<T: Real>(x: T) -> BayerPathSignatureChenLevyAreaSpec<T> {
            BayerPathSignatureChenLevyAreaSpec {
                time_grid: vec![T::zero(), T::from_f64(0.5), T::one()],
                spot_path: vec![
                    T::from_f64(100.0) + x,
                    T::from_f64(104.0),
                    T::from_f64(102.0),
                ],
                linear_functional_weights: [
                    T::from_f64(0.1),
                    T::from_f64(0.2),
                    T::from_f64(0.3),
                    T::from_f64(0.4),
                    T::from_f64(-0.2),
                    T::from_f64(0.1),
                ],
            }
        }
        #[test]
        fn suite46_dml_invariants_and_dual_fd() {
            let r = evaluate_bayer_path_signature_chen_and_levy_area(&spec(0.0));
            assert!(
                r.signed_levy_area.abs() > 1e-6
                    && r.chen_identity_l1_error < 1e-12
                    && r.shuffle_product_identity_error < 1e-12
            );
            let h = 1e-6;
            let fd = (evaluate_bayer_path_signature_chen_and_levy_area(&spec(h))
                .signature_linear_payoff
                - evaluate_bayer_path_signature_chen_and_levy_area(&spec(-h))
                    .signature_linear_payoff)
                / (2.0 * h);
            let ad = evaluate_bayer_path_signature_chen_and_levy_area(&spec(Dual::<1>::variable(
                0.0, 0,
            )))
            .signature_linear_payoff;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_bayer_path_signature_chen_and_levy_area_module::*;

mod evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;
    use ql_math::{inv_norm_cdf, norm_cdf};

    #[derive(Clone, Copy, Debug)]
    pub struct LopezDePradoDsrPurgedCvSpec<T: Real = f64> {
        pub observed_sharpe_ratio: T,
        pub sample_length_t: T,
        pub return_skewness: T,
        pub return_kurtosis: T,
        pub num_independent_trials_k: T,
        pub trials_sharpe_std: T,
        pub total_cv_samples_n: T,
        pub cv_folds_k: T,
        pub label_overlap_horizon_h: T,
        pub embargo_fraction_pct: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct LopezDePradoDsrPurgedCvReport<T: Real = f64> {
        pub expected_max_sharpe_hurdle_sr0: T,
        pub non_gaussian_sharpe_std_error: T,
        pub probabilistic_sharpe_ratio_zero: T,
        pub deflated_sharpe_ratio_dsr: T,
        pub non_gaussian_variance_inflation: T,
        pub unpurged_train_sample_fraction: T,
        pub purged_embargoed_train_fraction: T,
        pub purged_leakage_reduction_ratio: T,
    }

    pub fn evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv<T: Real>(
        spec: &LopezDePradoDsrPurgedCvSpec<T>,
    ) -> LopezDePradoDsrPurgedCvReport<T> {
        let gamma_em = T::from_f64(0.577_215_664_901_532_9);
        let e_val = T::one().exp();
        let k_tr = spec.num_independent_trials_k.max_val(T::from_f64(2.0));
        let p1 = (T::one() - T::one() / k_tr).min_val(T::from_f64(0.9999));
        let p2 = (T::one() - T::one() / (k_tr * e_val)).min_val(T::from_f64(0.9999));
        let sr0 = spec.trials_sharpe_std
            * ((T::one() - gamma_em) * inv_norm_cdf(p1) + gamma_em * inv_norm_cdf(p2));

        let sr = spec.observed_sharpe_ratio;
        let var_infl = (T::one() - spec.return_skewness * sr
            + T::from_f64(0.25) * (spec.return_kurtosis - T::one()) * sr * sr)
            .max_val(T::from_f64(1e-6));
        let t_m1 = (spec.sample_length_t - T::one()).max_val(T::from_f64(1.0));
        let se_sr = (var_infl / t_m1).sqrt();

        let psr0 = norm_cdf(sr / se_sr);
        let dsr = norm_cdf((sr - sr0) / se_sr);

        let k_folds = spec.cv_folds_k.max_val(T::from_f64(2.0));
        let n_obs = spec.total_cv_samples_n.max_val(T::from_f64(10.0));
        let unpurged_frac = (k_folds - T::one()) / k_folds;
        let purge_drop = T::from_f64(2.0) * k_folds * spec.label_overlap_horizon_h / n_obs;
        let purged_frac =
            (unpurged_frac - purge_drop - spec.embargo_fraction_pct).max_val(T::from_f64(0.05));
        let leak_red = (unpurged_frac - purged_frac) / unpurged_frac;

        LopezDePradoDsrPurgedCvReport {
            expected_max_sharpe_hurdle_sr0: sr0,
            non_gaussian_sharpe_std_error: se_sr,
            probabilistic_sharpe_ratio_zero: psr0,
            deflated_sharpe_ratio_dsr: dsr,
            non_gaussian_variance_inflation: var_infl,
            unpurged_train_sample_fraction: unpurged_frac,
            purged_embargoed_train_fraction: purged_frac,
            purged_leakage_reduction_ratio: leak_red,
        }
    }

    #[cfg(test)]
    mod domain_suite_48_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn lopez_de_prado_dsr_purged_cv_invariants_and_ad_parity() {
            let spec = LopezDePradoDsrPurgedCvSpec {
                observed_sharpe_ratio: 1.45,
                sample_length_t: 60.0,
                return_skewness: -0.65,
                return_kurtosis: 4.5,
                num_independent_trials_k: 25.0,
                trials_sharpe_std: 0.45,
                total_cv_samples_n: 1000.0,
                cv_folds_k: 5.0,
                label_overlap_horizon_h: 5.0,
                embargo_fraction_pct: 0.01,
            };
            let rep = evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv(&spec);
            assert!(rep.expected_max_sharpe_hurdle_sr0 > 0.0);
            assert!(rep.deflated_sharpe_ratio_dsr < rep.probabilistic_sharpe_ratio_zero);
            assert!(rep.purged_embargoed_train_fraction < rep.unpurged_train_sample_fraction);

            let h = 1e-6;
            let r_ad = evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv(
                &LopezDePradoDsrPurgedCvSpec {
                    observed_sharpe_ratio: Dual::<1>::variable(spec.observed_sharpe_ratio, 0),
                    sample_length_t: Dual::<1>::from_value(spec.sample_length_t),
                    return_skewness: Dual::<1>::from_value(spec.return_skewness),
                    return_kurtosis: Dual::<1>::from_value(spec.return_kurtosis),
                    num_independent_trials_k: Dual::<1>::from_value(spec.num_independent_trials_k),
                    trials_sharpe_std: Dual::<1>::from_value(spec.trials_sharpe_std),
                    total_cv_samples_n: Dual::<1>::from_value(spec.total_cv_samples_n),
                    cv_folds_k: Dual::<1>::from_value(spec.cv_folds_k),
                    label_overlap_horizon_h: Dual::<1>::from_value(spec.label_overlap_horizon_h),
                    embargo_fraction_pct: Dual::<1>::from_value(spec.embargo_fraction_pct),
                },
            );
            let mut sp = spec;
            sp.observed_sharpe_ratio += h;
            let mut sm = spec;
            sm.observed_sharpe_ratio -= h;
            let _fd = (evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv(&sp)
                .deflated_sharpe_ratio_dsr
                - evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv(&sm)
                    .deflated_sharpe_ratio_dsr)
                / (2.0 * h);
            assert!((r_ad.deflated_sharpe_ratio_dsr.d[0] - _fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv_module::*;

mod evaluate_tsay_johansen_vecm_cointegration_module {
    // Tsay (2010) & Hamilton (1994) Johansen Cointegration Vector Error-Correction (`VECM`),
    // Ornstein–Uhlenbeck Half-Life & Trace Statistic Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct TsayJohansenVecmCointegrationSpec<T: Real = f64> {
        pub cointegrating_hedge_ratio_beta: T,
        pub error_correction_speed_alpha1: T,
        pub error_correction_speed_alpha2: T,
        pub current_spread_deviation_zt: T,
        pub spread_innovation_vol_sigma_eps: T,
        pub canonical_eigenvalues_lambda: [T; 2],
        pub sample_size_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct TsayJohansenVecmCointegrationReport<T: Real = f64> {
        pub spread_ar1_persistence_phi: T,
        pub ou_mean_reversion_speed_kappa: T,
        pub mean_reversion_half_life_periods: T,
        pub stationary_spread_volatility: T,
        pub expected_next_period_spread: T,
        pub johansen_trace_statistic: T,
        pub johansen_max_eigen_statistic: T,
    }

    pub fn evaluate_tsay_johansen_vecm_cointegration<T: Real>(
        spec: &TsayJohansenVecmCointegrationSpec<T>,
    ) -> TsayJohansenVecmCointegrationReport<T> {
        // Under VECM: dy_{1,t} = alpha_1 z_{t-1} + e_{1,t}, dy_{2,t} = alpha_2 z_{t-1} + e_{2,t}
        // Therefore dz_t = dy_{1,t} - beta * dy_{2,t} = (alpha_1 - beta * alpha_2) z_{t-1} + e_{z,t}
        // => z_t = phi_z * z_{t-1} + e_{z,t} with phi_z = 1 + alpha_1 - beta * alpha_2:
        let phi = (T::one() + spec.error_correction_speed_alpha1
            - spec.cointegrating_hedge_ratio_beta * spec.error_correction_speed_alpha2)
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.995));
        let kappa = -phi.ln();
        let half_life = T::from_f64(std::f64::consts::LN_2) / kappa;

        let stat_vol = spec.spread_innovation_vol_sigma_eps
            / (T::one() - phi * phi).max_val(T::from_f64(1e-6)).sqrt();
        let next_z = phi * spec.current_spread_deviation_zt;

        // Johansen (1988/1991) likelihood-ratio statistics:
        let n_obs = spec.sample_size_t.max_val(T::from_f64(10.0));
        let l1 = spec.canonical_eigenvalues_lambda[0]
            .max_val(T::from_f64(1e-6))
            .min_val(T::from_f64(0.99));
        let l2 = spec.canonical_eigenvalues_lambda[1]
            .max_val(T::from_f64(1e-6))
            .min_val(T::from_f64(0.99));
        let max_eig = -n_obs * (T::one() - l1).ln();
        let trace_stat = -n_obs * ((T::one() - l1).ln() + (T::one() - l2).ln());

        TsayJohansenVecmCointegrationReport {
            spread_ar1_persistence_phi: phi,
            ou_mean_reversion_speed_kappa: kappa,
            mean_reversion_half_life_periods: half_life,
            stationary_spread_volatility: stat_vol,
            expected_next_period_spread: next_z,
            johansen_trace_statistic: trace_stat,
            johansen_max_eigen_statistic: max_eig,
        }
    }

    #[cfg(test)]
    mod domain_suite_52_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn tsay_johansen_vecm_cointegration_invariants_and_ad_parity() {
            let spec = TsayJohansenVecmCointegrationSpec {
                cointegrating_hedge_ratio_beta: 1.15,
                error_correction_speed_alpha1: -0.18,
                error_correction_speed_alpha2: 0.08,
                current_spread_deviation_zt: 1.4,
                spread_innovation_vol_sigma_eps: 0.25,
                canonical_eigenvalues_lambda: [0.16, 0.03],
                sample_size_t: 250.0,
            };
            let rep = evaluate_tsay_johansen_vecm_cointegration(&spec);
            assert!(rep.spread_ar1_persistence_phi > 0.0 && rep.spread_ar1_persistence_phi < 1.0);
            assert!(rep.johansen_trace_statistic > rep.johansen_max_eigen_statistic);
            assert!(rep.stationary_spread_volatility > spec.spread_innovation_vol_sigma_eps);

            let h = 1e-6;
            let r_ad =
                evaluate_tsay_johansen_vecm_cointegration(&TsayJohansenVecmCointegrationSpec {
                    cointegrating_hedge_ratio_beta: Dual::<1>::variable(
                        spec.cointegrating_hedge_ratio_beta,
                        0,
                    ),
                    error_correction_speed_alpha1: Dual::<1>::from_value(
                        spec.error_correction_speed_alpha1,
                    ),
                    error_correction_speed_alpha2: Dual::<1>::from_value(
                        spec.error_correction_speed_alpha2,
                    ),
                    current_spread_deviation_zt: Dual::<1>::from_value(
                        spec.current_spread_deviation_zt,
                    ),
                    spread_innovation_vol_sigma_eps: Dual::<1>::from_value(
                        spec.spread_innovation_vol_sigma_eps,
                    ),
                    canonical_eigenvalues_lambda: spec
                        .canonical_eigenvalues_lambda
                        .map(Dual::<1>::from_value),
                    sample_size_t: Dual::<1>::from_value(spec.sample_size_t),
                });
            let mut sp = spec;
            sp.cointegrating_hedge_ratio_beta += h;
            let mut sm = spec;
            sm.cointegrating_hedge_ratio_beta -= h;
            let fd = (evaluate_tsay_johansen_vecm_cointegration(&sp)
                .mean_reversion_half_life_periods
                - evaluate_tsay_johansen_vecm_cointegration(&sm).mean_reversion_half_life_periods)
                / (2.0 * h);
            assert!((r_ad.mean_reversion_half_life_periods.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_tsay_johansen_vecm_cointegration_module::*;

mod evaluate_tsay_star_threshold_autoregression_band_module {
    // Tsay (2010) & Teräsvirta (1994) Exponential Smooth-Transition Threshold
    // Autoregression (`ESTAR`) No-Trade Band & Regime Half-Life Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct TsayEstarThresholdAutoregressionSpec<T: Real = f64> {
        pub current_mispricing_z: T,
        pub transaction_cost_threshold_c: T,
        pub transition_smoothness_gamma: T,
        pub inner_band_reversion_speed_kappa_in: T,
        pub outer_arbitrage_reversion_speed_kappa_out: T,
        pub innovation_volatility_sigma_eps: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct TsayEstarThresholdAutoregressionReport<T: Real = f64> {
        pub transition_activation_weight_g: T,
        pub effective_local_reversion_speed_kappa: T,
        pub expected_one_step_spread_change: T,
        pub expected_next_period_mispricing: T,
        pub local_regime_half_life_periods: T,
        pub inner_vs_outer_half_life_ratio: T,
        pub local_stationary_variance_proxy: T,
    }

    pub fn evaluate_tsay_star_threshold_autoregression_band<T: Real>(
        spec: &TsayEstarThresholdAutoregressionSpec<T>,
    ) -> TsayEstarThresholdAutoregressionReport<T> {
        let z = spec.current_mispricing_z;
        let c = spec.transaction_cost_threshold_c.max_val(T::from_f64(1e-4));
        let gam = spec.transition_smoothness_gamma.max_val(T::from_f64(0.05));
        let k_in = spec
            .inner_band_reversion_speed_kappa_in
            .max_val(T::from_f64(0.01));
        let k_out = spec.outer_arbitrage_reversion_speed_kappa_out.max_val(k_in);

        // Smooth excess squared deviation above transaction-cost band c^2:
        let d_sq = z * z - c * c;
        let excess_sq = T::from_f64(0.5) * (d_sq + (d_sq * d_sq + T::from_f64(0.01)).sqrt());

        // ESTAR transition function G(z; gamma, c) = 1 - exp(-gamma * excess_sq) in (0, 1):
        let g = T::one() - (-gam * excess_sq).exp();
        let k_eff = k_in * (T::one() - g) + k_out * g;

        let dz = -k_eff * z;
        let z_next = z + dz;
        let phi_eff = (T::one() - k_eff)
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.995));
        let half_life = -T::from_f64(std::f64::consts::LN_2) / phi_eff.ln();

        let phi_in = (T::one() - k_in)
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.995));
        let phi_out = (T::one() - k_out)
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.995));
        let hl_ratio = phi_out.ln() / phi_in.ln();

        let stat_var = spec.innovation_volatility_sigma_eps * spec.innovation_volatility_sigma_eps
            / (T::one() - phi_eff * phi_eff).max_val(T::from_f64(1e-4));

        TsayEstarThresholdAutoregressionReport {
            transition_activation_weight_g: g,
            effective_local_reversion_speed_kappa: k_eff,
            expected_one_step_spread_change: dz,
            expected_next_period_mispricing: z_next,
            local_regime_half_life_periods: half_life,
            inner_vs_outer_half_life_ratio: hl_ratio,
            local_stationary_variance_proxy: stat_var,
        }
    }

    #[cfg(test)]
    mod domain_suite_54_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn tsay_star_threshold_autoregression_invariants_and_ad_parity() {
            let spec = TsayEstarThresholdAutoregressionSpec {
                current_mispricing_z: 1.85,
                transaction_cost_threshold_c: 0.90,
                transition_smoothness_gamma: 0.85,
                inner_band_reversion_speed_kappa_in: 0.04,
                outer_arbitrage_reversion_speed_kappa_out: 0.38,
                innovation_volatility_sigma_eps: 0.22,
            };
            let rep = evaluate_tsay_star_threshold_autoregression_band(&spec);
            assert!(rep.transition_activation_weight_g > 0.5);
            assert!(
                rep.effective_local_reversion_speed_kappa
                    > spec.inner_band_reversion_speed_kappa_in
            );
            assert!(rep.inner_vs_outer_half_life_ratio > 2.0);

            let h = 1e-6;
            let r_ad = evaluate_tsay_star_threshold_autoregression_band(
                &TsayEstarThresholdAutoregressionSpec {
                    current_mispricing_z: Dual::<1>::variable(spec.current_mispricing_z, 0),
                    transaction_cost_threshold_c: Dual::<1>::from_value(
                        spec.transaction_cost_threshold_c,
                    ),
                    transition_smoothness_gamma: Dual::<1>::from_value(
                        spec.transition_smoothness_gamma,
                    ),
                    inner_band_reversion_speed_kappa_in: Dual::<1>::from_value(
                        spec.inner_band_reversion_speed_kappa_in,
                    ),
                    outer_arbitrage_reversion_speed_kappa_out: Dual::<1>::from_value(
                        spec.outer_arbitrage_reversion_speed_kappa_out,
                    ),
                    innovation_volatility_sigma_eps: Dual::<1>::from_value(
                        spec.innovation_volatility_sigma_eps,
                    ),
                },
            );
            let mut sp = spec;
            sp.current_mispricing_z += h;
            let mut sm = spec;
            sm.current_mispricing_z -= h;
            let fd = (evaluate_tsay_star_threshold_autoregression_band(&sp)
                .local_regime_half_life_periods
                - evaluate_tsay_star_threshold_autoregression_band(&sm)
                    .local_regime_half_life_periods)
                / (2.0 * h);
            assert!((r_ad.local_regime_half_life_periods.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_tsay_star_threshold_autoregression_band_module::*;

mod evaluate_newey_west_andrews_hac_variance_module {
    // Newey & West (1987) & Andrews (1991) Automatic-Bandwidth (`ARMA(1,1)` Plug-In) `HAC`
    // Long-Run Variance, Bartlett & Quadratic-Spectral (`QS`) Kernel Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct NeweyWestAndrewsHacVarianceSpec<T: Real = f64> {
        pub sample_variance_gamma0: T,
        pub autocovariance_lag1_gamma1: T,
        pub autocovariance_lag2_gamma2: T,
        pub autocovariance_lag3_gamma3: T,
        pub sample_size_n: T,
        pub ols_point_estimate_beta: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct NeweyWestAndrewsHacVarianceReport<T: Real = f64> {
        pub newey_west_optimal_bandwidth: T,
        pub andrews_qs_optimal_bandwidth: T,
        pub bartlett_long_run_variance: T,
        pub quadratic_spectral_long_run_variance: T,
        pub hac_variance_inflation_factor: T,
        pub hac_robust_standard_error: T,
        pub hac_robust_t_statistic: T,
    }

    fn qs_kernel<T: Real>(x: T) -> T {
        let x_pos = x.max_val(T::from_f64(1e-4));
        let x2 = x_pos * x_pos;
        (T::from_f64(3.0) / x2) * (x_pos.sin() / x_pos - x_pos.cos())
    }

    pub fn evaluate_newey_west_andrews_hac_variance<T: Real>(
        spec: &NeweyWestAndrewsHacVarianceSpec<T>,
    ) -> NeweyWestAndrewsHacVarianceReport<T> {
        let g0 = spec.sample_variance_gamma0.max_val(T::from_f64(1e-6));
        let g1 = spec.autocovariance_lag1_gamma1;
        let g2 = spec.autocovariance_lag2_gamma2;
        let g3 = spec.autocovariance_lag3_gamma3;
        let n = spec.sample_size_n.max_val(T::from_f64(20.0));

        let rho1 = (g1 / g0)
            .max_val(T::from_f64(-0.92))
            .min_val(T::from_f64(0.92));
        let om_r2 = (T::one() - rho1 * rho1).max_val(T::from_f64(1e-4));
        let om_r = (T::one() - rho1).max_val(T::from_f64(1e-4));

        // Andrews (1991, Econometrica Eq. 6.2 & 6.4) AR(1) plug-in smoothness parameters:
        let alpha1 = (T::from_f64(4.0) * rho1 * rho1) / (om_r2 * om_r2);
        let alpha2 = (T::from_f64(4.0) * rho1 * rho1) / (om_r * om_r * om_r * om_r);

        let s_nw = T::from_f64(1.1447) * (alpha1 * n).pow(T::from_f64(1.0 / 3.0));
        let s_qs = T::from_f64(1.3221) * (alpha2 * n).pow(T::from_f64(0.20));

        let eps = T::from_f64(0.02);
        let smooth_pos = |u: T| T::from_f64(0.5) * (u + (u * u + eps * eps).sqrt());
        let w1_nw = smooth_pos(T::one() - T::one() / s_nw);
        let w2_nw = smooth_pos(T::one() - T::from_f64(2.0) / s_nw);
        let w3_nw = smooth_pos(T::one() - T::from_f64(3.0) / s_nw);
        let omega_nw = g0 + T::from_f64(2.0) * (w1_nw * g1 + w2_nw * g2 + w3_nw * g3);

        let scale_qs = T::from_f64(6.0 * std::f64::consts::PI / 5.0) / s_qs;
        let w1_qs = qs_kernel(scale_qs);
        let w2_qs = qs_kernel(T::from_f64(2.0) * scale_qs);
        let w3_qs = qs_kernel(T::from_f64(3.0) * scale_qs);
        let omega_qs = (g0 + T::from_f64(2.0) * (w1_qs * g1 + w2_qs * g2 + w3_qs * g3))
            .max_val(T::from_f64(1e-6));

        let vif = omega_qs / g0;
        let se_hac = (omega_qs / n).sqrt();
        let t_hac = spec.ols_point_estimate_beta / se_hac;

        NeweyWestAndrewsHacVarianceReport {
            newey_west_optimal_bandwidth: s_nw,
            andrews_qs_optimal_bandwidth: s_qs,
            bartlett_long_run_variance: omega_nw,
            quadratic_spectral_long_run_variance: omega_qs,
            hac_variance_inflation_factor: vif,
            hac_robust_standard_error: se_hac,
            hac_robust_t_statistic: t_hac,
        }
    }

    #[cfg(test)]
    mod domain_suite_55_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn newey_west_andrews_hac_variance_invariants_and_ad_parity() {
            let spec = NeweyWestAndrewsHacVarianceSpec {
                sample_variance_gamma0: 0.040,
                autocovariance_lag1_gamma1: 0.022,
                autocovariance_lag2_gamma2: 0.012,
                autocovariance_lag3_gamma3: 0.006,
                sample_size_n: 250.0,
                ols_point_estimate_beta: 0.035,
            };
            let rep = evaluate_newey_west_andrews_hac_variance(&spec);
            assert!(rep.quadratic_spectral_long_run_variance > spec.sample_variance_gamma0);
            assert!(rep.hac_variance_inflation_factor > 1.5);
            assert!(rep.andrews_qs_optimal_bandwidth > 1.0);

            let h = 1e-6;
            let r_ad = evaluate_newey_west_andrews_hac_variance(&NeweyWestAndrewsHacVarianceSpec {
                sample_variance_gamma0: Dual::<1>::variable(spec.sample_variance_gamma0, 0),
                autocovariance_lag1_gamma1: Dual::<1>::from_value(spec.autocovariance_lag1_gamma1),
                autocovariance_lag2_gamma2: Dual::<1>::from_value(spec.autocovariance_lag2_gamma2),
                autocovariance_lag3_gamma3: Dual::<1>::from_value(spec.autocovariance_lag3_gamma3),
                sample_size_n: Dual::<1>::from_value(spec.sample_size_n),
                ols_point_estimate_beta: Dual::<1>::from_value(spec.ols_point_estimate_beta),
            });
            let mut sp = spec;
            sp.sample_variance_gamma0 += h;
            let mut sm = spec;
            sm.sample_variance_gamma0 -= h;
            let fd = (evaluate_newey_west_andrews_hac_variance(&sp).hac_robust_t_statistic
                - evaluate_newey_west_andrews_hac_variance(&sm).hac_robust_t_statistic)
                / (2.0 * h);
            assert!((r_ad.hac_robust_t_statistic.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_newey_west_andrews_hac_variance_module::*;

mod evaluate_diebold_mariano_hln_forecast_test_module {
    // Diebold & Mariano (1995) & Harvey–Leybourne–Newbold (`HLN`, 1997) Predictive Accuracy
    // Loss-Differential Test & Bates–Granger Optimal Forecast Combination Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct DieboldMarianoHlnForecastTestSpec<T: Real = f64> {
        pub model1_mse_loss_l1: T,
        pub model2_mse_loss_l2: T,
        pub forecast_error_correlation_rho12: T,
        pub loss_differential_autocov_lag1: T,
        pub forecast_horizon_steps_h: T,
        pub evaluation_sample_size_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct DieboldMarianoHlnForecastTestReport<T: Real = f64> {
        pub mean_loss_differential_dbar: T,
        pub long_run_loss_differential_variance: T,
        pub raw_diebold_mariano_statistic: T,
        pub hln_small_sample_correction_factor: T,
        pub hln_corrected_dm_statistic: T,
        pub bates_granger_optimal_weight_model1: T,
        pub combined_forecast_mse_reduction: T,
    }

    pub fn evaluate_diebold_mariano_hln_forecast_test<T: Real>(
        spec: &DieboldMarianoHlnForecastTestSpec<T>,
    ) -> DieboldMarianoHlnForecastTestReport<T> {
        let l1 = spec.model1_mse_loss_l1.max_val(T::from_f64(1e-6));
        let l2 = spec.model2_mse_loss_l2.max_val(T::from_f64(1e-6));
        let rho = spec
            .forecast_error_correlation_rho12
            .max_val(T::from_f64(-0.95))
            .min_val(T::from_f64(0.95));
        let gam1 = spec.loss_differential_autocov_lag1;
        let h_step = spec.forecast_horizon_steps_h.max_val(T::one());
        let t_obs = spec.evaluation_sample_size_t.max_val(T::from_f64(15.0));

        let d_bar = l1 - l2;
        let sig1 = l1.sqrt();
        let sig2 = l2.sqrt();
        let cov12 = rho * sig1 * sig2;

        // Variance of squared-error loss differential d_t = e_{1,t}^2 - e_{2,t}^2 under Gaussian errors:
        // Var(d_t) = 2 (sigma_1^4 + sigma_2^4 - 2 cov_{12}^2)
        let gamma0 = T::from_f64(2.0) * (l1 * l1 + l2 * l2 - T::from_f64(2.0) * cov12 * cov12);
        let m_trunc = (h_step - T::one()).max_val(T::zero());
        let w_h = m_trunc / (m_trunc + T::one());
        let lrv = (gamma0 + T::from_f64(2.0) * w_h * gam1).max_val(T::from_f64(1e-6));

        let se_dm = (lrv / t_obs).sqrt();
        let dm_raw = d_bar / se_dm;

        // Harvey, Leybourne & Newbold (1997, IJF Eq. 3) small-sample bias correction:
        let hln_num = (t_obs + T::one() - T::from_f64(2.0) * h_step
            + (h_step * (h_step - T::one())) / t_obs)
            .max_val(T::from_f64(1.0));
        let hln_factor = (hln_num / t_obs).sqrt();
        let dm_hln = dm_raw * hln_factor;

        // Bates & Granger (1969) optimal combination weight w1* minimizing combined MSE:
        let denom_bg = (l1 + l2 - T::from_f64(2.0) * cov12).max_val(T::from_f64(1e-6));
        let w1_opt = ((l2 - cov12) / denom_bg)
            .max_val(T::zero())
            .min_val(T::one());
        let w2_opt = T::one() - w1_opt;
        let mse_comb = w1_opt * w1_opt * l1
            + w2_opt * w2_opt * l2
            + T::from_f64(2.0) * w1_opt * w2_opt * cov12;
        let mse_min = l1.min_val(l2);
        let mse_red = (mse_min - mse_comb).max_val(T::zero());

        DieboldMarianoHlnForecastTestReport {
            mean_loss_differential_dbar: d_bar,
            long_run_loss_differential_variance: lrv,
            raw_diebold_mariano_statistic: dm_raw,
            hln_small_sample_correction_factor: hln_factor,
            hln_corrected_dm_statistic: dm_hln,
            bates_granger_optimal_weight_model1: w1_opt,
            combined_forecast_mse_reduction: mse_red,
        }
    }

    #[cfg(test)]
    mod evaluate_diebold_mariano_hln_forecast_test_module_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn diebold_mariano_hln_forecast_test_invariants_and_ad_parity() {
            let spec = DieboldMarianoHlnForecastTestSpec {
                model1_mse_loss_l1: 0.025,
                model2_mse_loss_l2: 0.036,
                forecast_error_correlation_rho12: 0.45,
                loss_differential_autocov_lag1: 0.0004,
                forecast_horizon_steps_h: 3.0,
                evaluation_sample_size_t: 80.0,
            };
            let rep = evaluate_diebold_mariano_hln_forecast_test(&spec);
            assert!(rep.raw_diebold_mariano_statistic < 0.0);
            assert!(rep.hln_small_sample_correction_factor < 1.0);
            assert!(rep.bates_granger_optimal_weight_model1 > 0.5);
            assert!(rep.combined_forecast_mse_reduction > 0.0);

            let h = 1e-6;
            let r_ad =
                evaluate_diebold_mariano_hln_forecast_test(&DieboldMarianoHlnForecastTestSpec {
                    model1_mse_loss_l1: Dual::<1>::variable(spec.model1_mse_loss_l1, 0),
                    model2_mse_loss_l2: Dual::<1>::from_value(spec.model2_mse_loss_l2),
                    forecast_error_correlation_rho12: Dual::<1>::from_value(
                        spec.forecast_error_correlation_rho12,
                    ),
                    loss_differential_autocov_lag1: Dual::<1>::from_value(
                        spec.loss_differential_autocov_lag1,
                    ),
                    forecast_horizon_steps_h: Dual::<1>::from_value(spec.forecast_horizon_steps_h),
                    evaluation_sample_size_t: Dual::<1>::from_value(spec.evaluation_sample_size_t),
                });
            let mut sp = spec;
            sp.model1_mse_loss_l1 += h;
            let mut sm = spec;
            sm.model1_mse_loss_l1 -= h;
            let fd = (evaluate_diebold_mariano_hln_forecast_test(&sp).hln_corrected_dm_statistic
                - evaluate_diebold_mariano_hln_forecast_test(&sm).hln_corrected_dm_statistic)
                / (2.0 * h);
            assert!((r_ad.hln_corrected_dm_statistic.d[0] - fd).abs() < 1e-4);
        }
    }
}
pub use evaluate_diebold_mariano_hln_forecast_test_module::*;

mod evaluate_engle_granger_cointegration_ecm_half_life_module {
    // Engle & Granger (1987, Econometrica) & Tsay (2010, Ch. 8) Residual Cointegration
    // ADF Test, Error-Correction Model (`ECM`) Speed & Ornstein–Uhlenbeck Half-Life Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct EngleGrangerCointegrationEcmHalfLifeSpec<T: Real = f64> {
        pub variance_asset_x: T,
        pub covariance_xy: T,
        pub spread_lag1_autocorrelation_phi: T,
        pub spread_innovation_sigma_e: T,
        pub current_disequilibrium_spread_z: T,
        pub sample_observations_t: T,
        pub time_step_years_dt: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct EngleGrangerCointegrationEcmHalfLifeReport<T: Real = f64> {
        pub cointegrating_hedge_ratio_beta: T,
        pub error_correction_speed_alpha: T,
        pub engle_granger_adf_t_statistic: T,
        pub ornstein_uhlenbeck_kappa_annual: T,
        pub mean_reversion_half_life_periods: T,
        pub stationary_spread_variance: T,
        pub expected_next_step_spread_change: T,
    }

    pub fn evaluate_engle_granger_cointegration_ecm_half_life<T: Real>(
        spec: &EngleGrangerCointegrationEcmHalfLifeSpec<T>,
    ) -> EngleGrangerCointegrationEcmHalfLifeReport<T> {
        let var_x = spec.variance_asset_x.max_val(T::from_f64(1e-6));
        let cov_xy = spec.covariance_xy;
        let phi = spec
            .spread_lag1_autocorrelation_phi
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.995));
        let sig_e = spec.spread_innovation_sigma_e.max_val(T::from_f64(1e-4));
        let z_curr = spec.current_disequilibrium_spread_z;
        let t_obs = spec.sample_observations_t.max_val(T::from_f64(20.0));
        let dt = spec.time_step_years_dt.max_val(T::from_f64(1e-4));

        let beta_coint = cov_xy / var_x;
        let alpha_ecm = phi - T::one();

        // Asymptotic SE(phi_hat) = sqrt((1 - phi^2) / T) and Engle-Granger ADF t-stat = (phi - 1) / SE(phi):
        let se_phi = ((T::one() - phi * phi) / t_obs).sqrt();
        let adf_t = alpha_ecm / se_phi;

        let ou_kappa = -phi.ln() / dt;
        let ln2 = T::from_f64(std::f64::consts::LN_2);
        let half_life_periods = -ln2 / phi.ln();
        let stat_var = (sig_e * sig_e) / (T::one() - phi * phi);
        let z_coint = z_curr - (beta_coint - T::one()) * var_x.sqrt();
        let exp_dz = alpha_ecm * z_coint;

        EngleGrangerCointegrationEcmHalfLifeReport {
            cointegrating_hedge_ratio_beta: beta_coint,
            error_correction_speed_alpha: alpha_ecm,
            engle_granger_adf_t_statistic: adf_t,
            ornstein_uhlenbeck_kappa_annual: ou_kappa,
            mean_reversion_half_life_periods: half_life_periods,
            stationary_spread_variance: stat_var,
            expected_next_step_spread_change: exp_dz,
        }
    }

    #[cfg(test)]
    mod evaluate_engle_granger_cointegration_ecm_half_life_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn engle_granger_cointegration_ecm_half_life_invariants_and_ad_parity() {
            let spec = EngleGrangerCointegrationEcmHalfLifeSpec {
                variance_asset_x: 0.042,
                covariance_xy: 0.038,
                spread_lag1_autocorrelation_phi: 0.86,
                spread_innovation_sigma_e: 0.015,
                current_disequilibrium_spread_z: 0.025,
                sample_observations_t: 250.0,
                time_step_years_dt: 1.0 / 252.0,
            };
            let rep = evaluate_engle_granger_cointegration_ecm_half_life(&spec);
            assert!(rep.engle_granger_adf_t_statistic < -3.0);
            assert!(rep.mean_reversion_half_life_periods > 1.0);
            assert!(rep.expected_next_step_spread_change < 0.0);

            let h = 1e-6;
            let r_ad = evaluate_engle_granger_cointegration_ecm_half_life(
                &EngleGrangerCointegrationEcmHalfLifeSpec {
                    variance_asset_x: Dual::<7>::variable(spec.variance_asset_x, 0),
                    covariance_xy: Dual::<7>::variable(spec.covariance_xy, 1),
                    spread_lag1_autocorrelation_phi: Dual::<7>::variable(
                        spec.spread_lag1_autocorrelation_phi,
                        2,
                    ),
                    spread_innovation_sigma_e: Dual::<7>::variable(
                        spec.spread_innovation_sigma_e,
                        3,
                    ),
                    current_disequilibrium_spread_z: Dual::<7>::variable(
                        spec.current_disequilibrium_spread_z,
                        4,
                    ),
                    sample_observations_t: Dual::<7>::variable(spec.sample_observations_t, 5),
                    time_step_years_dt: Dual::<7>::variable(spec.time_step_years_dt, 6),
                },
            );
            let composite_ad = r_ad.engle_granger_adf_t_statistic * Dual::<7>::from_value(0.01)
                + r_ad.ornstein_uhlenbeck_kappa_annual * Dual::<7>::from_value(1e-3)
                + r_ad.stationary_spread_variance
                + r_ad.expected_next_step_spread_change;
            for (idx, &ad) in composite_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.variance_asset_x += h;
                        sm.variance_asset_x -= h;
                    }
                    1 => {
                        sp.covariance_xy += h;
                        sm.covariance_xy -= h;
                    }
                    2 => {
                        sp.spread_lag1_autocorrelation_phi += h;
                        sm.spread_lag1_autocorrelation_phi -= h;
                    }
                    3 => {
                        sp.spread_innovation_sigma_e += h;
                        sm.spread_innovation_sigma_e -= h;
                    }
                    4 => {
                        sp.current_disequilibrium_spread_z += h;
                        sm.current_disequilibrium_spread_z -= h;
                    }
                    5 => {
                        sp.sample_observations_t += h;
                        sm.sample_observations_t -= h;
                    }
                    _ => {
                        sp.time_step_years_dt += h;
                        sm.time_step_years_dt -= h;
                    }
                }
                let rp = evaluate_engle_granger_cointegration_ecm_half_life(&sp);
                let rm = evaluate_engle_granger_cointegration_ecm_half_life(&sm);
                let fp = rp.engle_granger_adf_t_statistic * 0.01
                    + rp.ornstein_uhlenbeck_kappa_annual * 1e-3
                    + rp.stationary_spread_variance
                    + rp.expected_next_step_spread_change;
                let fm = rm.engle_granger_adf_t_statistic * 0.01
                    + rm.ornstein_uhlenbeck_kappa_annual * 1e-3
                    + rm.stationary_spread_variance
                    + rm.expected_next_step_spread_change;
                let fd = (fp - fm) / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_engle_granger_cointegration_ecm_half_life_module::*;

mod evaluate_barndorff_nielsen_shephard_bipower_jump_har_cj_module {
    // Barndorff-Nielsen & Shephard (BNS, 2004 / 2006, J. Financial Econometrics) &
    // Andersen, Bollerslev & Diebold (HAR-CJ, 2007, Econometrica) Realized Bipower Variation (BV),
    // Jump Detection Z-Test & Continuous-Jump HAR-CJ Volatility Forecast Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BarndorffNielsenShephardHarCjSpec<T: Real = f64> {
        pub daily_realized_variance_rv: T,
        pub daily_bipower_variation_bv: T,
        pub daily_tripower_quarticity_tq: T,
        pub weekly_continuous_bv_w: T,
        pub monthly_continuous_bv_m: T,
        pub intraday_sampling_frequency_m: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BarndorffNielsenShephardHarCjReport<T: Real = f64> {
        pub continuous_integrated_variance_est: T,
        pub realized_jump_quadratic_variation: T,
        pub relative_jump_contribution_rj: T,
        pub bns_asymptotic_jump_z_statistic: T,
        pub har_cj_next_day_variance_forecast: T,
        pub har_cj_annualized_volatility_forecast: T,
    }

    pub fn evaluate_barndorff_nielsen_shephard_bipower_jump_har_cj<T: Real>(
        spec: &BarndorffNielsenShephardHarCjSpec<T>,
    ) -> BarndorffNielsenShephardHarCjReport<T> {
        let rv = spec.daily_realized_variance_rv.max_val(T::from_f64(1e-6));
        let bv = spec
            .daily_bipower_variation_bv
            .max_val(T::from_f64(1e-6))
            .min_val(rv);
        let tq = spec.daily_tripower_quarticity_tq.max_val(bv * bv);
        let bv_w = spec.weekly_continuous_bv_w.max_val(T::from_f64(1e-6));
        let bv_m = spec.monthly_continuous_bv_m.max_val(T::from_f64(1e-6));
        let m_obs = spec
            .intraday_sampling_frequency_m
            .max_val(T::from_f64(12.0));

        let raw_jump = rv - bv;
        let eps = T::from_f64(1e-6);
        let jump_qv = T::from_f64(0.5) * (raw_jump + (raw_jump * raw_jump + eps * eps).sqrt());
        let rj = (rv - bv) / rv;

        // Barndorff-Nielsen & Shephard (2006) relative jump test variance constant:
        // theta_BNS = (pi / 2)^2 + pi - 5 = 0.60906898
        let theta_bns = T::from_f64(0.609_068_98);
        let quart_ratio = (tq / (bv * bv)).max_val(T::one());
        let z_bns = rj / ((theta_bns / m_obs) * quart_ratio).sqrt();

        // Andersen, Bollerslev & Diebold (2007) HAR-CJ forecast:
        // E_t[RV_{t+1}] = beta_0 + beta_d * C_d + beta_w * C_w + beta_m * C_m + beta_j * J_d
        let har_cj_var = T::from_f64(0.000015)
            + T::from_f64(0.38) * bv
            + T::from_f64(0.34) * bv_w
            + T::from_f64(0.20) * bv_m
            + T::from_f64(0.08) * jump_qv;
        let har_cj_ann_vol = (har_cj_var * T::from_f64(252.0)).sqrt();

        BarndorffNielsenShephardHarCjReport {
            continuous_integrated_variance_est: bv,
            realized_jump_quadratic_variation: jump_qv,
            relative_jump_contribution_rj: rj,
            bns_asymptotic_jump_z_statistic: z_bns,
            har_cj_next_day_variance_forecast: har_cj_var,
            har_cj_annualized_volatility_forecast: har_cj_ann_vol,
        }
    }

    #[cfg(test)]
    mod evaluate_barndorff_nielsen_shephard_bipower_jump_har_cj_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn barndorff_nielsen_shephard_bipower_jump_har_cj_invariants_and_ad_parity() {
            let spec = BarndorffNielsenShephardHarCjSpec {
                daily_realized_variance_rv: 0.00042,
                daily_bipower_variation_bv: 0.00028,
                daily_tripower_quarticity_tq: 1.15e-7,
                weekly_continuous_bv_w: 0.00024,
                monthly_continuous_bv_m: 0.00020,
                intraday_sampling_frequency_m: 78.0,
            };
            let rep = evaluate_barndorff_nielsen_shephard_bipower_jump_har_cj(&spec);
            assert!(rep.realized_jump_quadratic_variation > 0.0);
            assert!(rep.bns_asymptotic_jump_z_statistic > 2.0);
            assert!(rep.har_cj_annualized_volatility_forecast > 0.10);

            let h = 1e-8;
            let r_ad = evaluate_barndorff_nielsen_shephard_bipower_jump_har_cj(
                &BarndorffNielsenShephardHarCjSpec {
                    daily_realized_variance_rv: Dual::<6>::variable(
                        spec.daily_realized_variance_rv,
                        0,
                    ),
                    daily_bipower_variation_bv: Dual::<6>::variable(
                        spec.daily_bipower_variation_bv,
                        1,
                    ),
                    daily_tripower_quarticity_tq: Dual::<6>::variable(
                        spec.daily_tripower_quarticity_tq,
                        2,
                    ),
                    weekly_continuous_bv_w: Dual::<6>::variable(spec.weekly_continuous_bv_w, 3),
                    monthly_continuous_bv_m: Dual::<6>::variable(spec.monthly_continuous_bv_m, 4),
                    intraday_sampling_frequency_m: Dual::<6>::variable(
                        spec.intraday_sampling_frequency_m,
                        5,
                    ),
                },
            );
            let comp_ad = r_ad.bns_asymptotic_jump_z_statistic * Dual::<6>::from_value(1e-3)
                + r_ad.har_cj_annualized_volatility_forecast;
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                let step = if idx == 2 { 1e-11 } else { h };
                match idx {
                    0 => {
                        sp.daily_realized_variance_rv += step;
                        sm.daily_realized_variance_rv -= step;
                    }
                    1 => {
                        sp.daily_bipower_variation_bv += step;
                        sm.daily_bipower_variation_bv -= step;
                    }
                    2 => {
                        sp.daily_tripower_quarticity_tq += step;
                        sm.daily_tripower_quarticity_tq -= step;
                    }
                    3 => {
                        sp.weekly_continuous_bv_w += step;
                        sm.weekly_continuous_bv_w -= step;
                    }
                    4 => {
                        sp.monthly_continuous_bv_m += step;
                        sm.monthly_continuous_bv_m -= step;
                    }
                    _ => {
                        sp.intraday_sampling_frequency_m += step;
                        sm.intraday_sampling_frequency_m -= step;
                    }
                }
                let rp = evaluate_barndorff_nielsen_shephard_bipower_jump_har_cj(&sp);
                let rm = evaluate_barndorff_nielsen_shephard_bipower_jump_har_cj(&sm);
                let fp = rp.bns_asymptotic_jump_z_statistic * 1e-3
                    + rp.har_cj_annualized_volatility_forecast;
                let fm = rm.bns_asymptotic_jump_z_statistic * 1e-3
                    + rm.har_cj_annualized_volatility_forecast;
                let fd = (fp - fm) / (2.0 * step);
                assert!(ad.abs() > 1e-5);
                assert!((ad - fd).abs() / ad.abs().max(1.0) < 1e-3);
            }
        }
    }
}
pub use evaluate_barndorff_nielsen_shephard_bipower_jump_har_cj_module::*;

mod evaluate_pesaran_timmermann_directional_timing_hit_ratio_test_module {
    // Pesaran & Timmermann (1992, J. Business & Economic Statistics, "A Simple Nonparametric Test
    // of Predictive Performance") & Anatolyev-Gerko (2005) Directional Market-Timing Hit-Ratio Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct PesaranTimmermannDirectionalTimingSpec<T: Real = f64> {
        pub true_positive_up_fraction_p11: T,
        pub true_negative_down_fraction_p00: T,
        pub false_positive_up_fraction_p01: T,
        pub sample_observations_t: T,
        pub mean_absolute_return_move_bps: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct PesaranTimmermannDirectionalTimingReport<T: Real = f64> {
        pub empirical_directional_hit_ratio_p_hat: T,
        pub independent_benchmark_hit_ratio_p_star: T,
        pub kuipers_henriksson_merton_skill_score: T,
        pub pesaran_timmermann_z_statistic: T,
        pub directional_timing_expected_edge_bps: T,
    }

    pub fn evaluate_pesaran_timmermann_directional_timing_hit_ratio_test<T: Real>(
        spec: &PesaranTimmermannDirectionalTimingSpec<T>,
    ) -> PesaranTimmermannDirectionalTimingReport<T> {
        let p11 = spec
            .true_positive_up_fraction_p11
            .max_val(T::from_f64(0.05));
        let p00 = spec
            .true_negative_down_fraction_p00
            .max_val(T::from_f64(0.05));
        let p01 = spec
            .false_positive_up_fraction_p01
            .max_val(T::from_f64(0.02));
        let n_obs = spec.sample_observations_t.max_val(T::from_f64(20.0));
        let abs_move_bps = spec.mean_absolute_return_move_bps.max_val(T::from_f64(1.0));

        let p10 = (T::one() - p11 - p00 - p01).max_val(T::from_f64(0.02));
        let p_hat = p11 + p00;
        let p_y = p11 + p10;
        let p_yhat = p11 + p01;

        // Null benchmark hit rate under independence between y_t and \hat{y}_t:
        let p_star = p_y * p_yhat + (T::one() - p_y) * (T::one() - p_yhat);

        // Pesaran & Timmermann (1992, Eqs. 5-6) exact finite-sample variances:
        let var_p_hat = (p_star * (T::one() - p_star)) / n_obs;
        let two_py_m1 = T::from_f64(2.0) * p_y - T::one();
        let two_pyhat_m1 = T::from_f64(2.0) * p_yhat - T::one();
        let var_p_star = (two_py_m1 * two_py_m1 * p_yhat * (T::one() - p_yhat)
            + two_pyhat_m1 * two_pyhat_m1 * p_y * (T::one() - p_y)
            + (T::from_f64(4.0) / n_obs) * p_y * p_yhat * (T::one() - p_y) * (T::one() - p_yhat))
            / n_obs;

        let var_diff = (var_p_hat - var_p_star).max_val(T::from_f64(1e-8));
        let pt_stat = (p_hat - p_star) / var_diff.sqrt();

        let kss = (p11 / p_y) - (p01 / (T::one() - p_y));
        let edge_bps = (T::from_f64(2.0) * (p_hat - p_star)) * abs_move_bps;

        PesaranTimmermannDirectionalTimingReport {
            empirical_directional_hit_ratio_p_hat: p_hat,
            independent_benchmark_hit_ratio_p_star: p_star,
            kuipers_henriksson_merton_skill_score: kss,
            pesaran_timmermann_z_statistic: pt_stat,
            directional_timing_expected_edge_bps: edge_bps,
        }
    }

    #[cfg(test)]
    mod evaluate_pesaran_timmermann_directional_timing_hit_ratio_test_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn pesaran_timmermann_directional_timing_invariants_and_ad_parity() {
            let spec = PesaranTimmermannDirectionalTimingSpec {
                true_positive_up_fraction_p11: 0.34,
                true_negative_down_fraction_p00: 0.25,
                false_positive_up_fraction_p01: 0.20,
                sample_observations_t: 250.0,
                mean_absolute_return_move_bps: 85.0,
            };
            let rep = evaluate_pesaran_timmermann_directional_timing_hit_ratio_test(&spec);
            assert!(
                rep.empirical_directional_hit_ratio_p_hat
                    > rep.independent_benchmark_hit_ratio_p_star
            );
            assert!(rep.kuipers_henriksson_merton_skill_score > 0.10);
            assert!(rep.pesaran_timmermann_z_statistic > 2.0);

            let h = 1e-5;
            let r_ad = evaluate_pesaran_timmermann_directional_timing_hit_ratio_test(
                &PesaranTimmermannDirectionalTimingSpec {
                    true_positive_up_fraction_p11: Dual::<5>::variable(
                        spec.true_positive_up_fraction_p11,
                        0,
                    ),
                    true_negative_down_fraction_p00: Dual::<5>::variable(
                        spec.true_negative_down_fraction_p00,
                        1,
                    ),
                    false_positive_up_fraction_p01: Dual::<5>::variable(
                        spec.false_positive_up_fraction_p01,
                        2,
                    ),
                    sample_observations_t: Dual::<5>::variable(spec.sample_observations_t, 3),
                    mean_absolute_return_move_bps: Dual::<5>::variable(
                        spec.mean_absolute_return_move_bps,
                        4,
                    ),
                },
            );
            let comp_ad =
                r_ad.pesaran_timmermann_z_statistic + r_ad.directional_timing_expected_edge_bps;
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.true_positive_up_fraction_p11 += h;
                        sm.true_positive_up_fraction_p11 -= h;
                    }
                    1 => {
                        sp.true_negative_down_fraction_p00 += h;
                        sm.true_negative_down_fraction_p00 -= h;
                    }
                    2 => {
                        sp.false_positive_up_fraction_p01 += h;
                        sm.false_positive_up_fraction_p01 -= h;
                    }
                    3 => {
                        sp.sample_observations_t += h;
                        sm.sample_observations_t -= h;
                    }
                    _ => {
                        sp.mean_absolute_return_move_bps += h;
                        sm.mean_absolute_return_move_bps -= h;
                    }
                }
                let rp = evaluate_pesaran_timmermann_directional_timing_hit_ratio_test(&sp);
                let rm = evaluate_pesaran_timmermann_directional_timing_hit_ratio_test(&sm);
                let fd = ((rp.pesaran_timmermann_z_statistic
                    + rp.directional_timing_expected_edge_bps)
                    - (rm.pesaran_timmermann_z_statistic
                        + rm.directional_timing_expected_edge_bps))
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_pesaran_timmermann_directional_timing_hit_ratio_test_module::*;
