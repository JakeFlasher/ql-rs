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
        let gamma_em = T::from_f64(std::f64::consts::EULER_GAMMA);
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

mod evaluate_bollerslev_tsay_garch_term_structure_kurtosis_module {
    // Bollerslev (1986, J. Econometrics, "Generalized Autoregressive Conditional Heteroskedasticity")
    // & Tsay (2010, Analysis of Financial Time Series, 3rd ed., Ch. 3) GARCH(1,1) Multi-Step
    // Variance Term Structure, Shock Half-Life & Fourth-Moment Unconditional Excess Kurtosis Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BollerslevTsayGarchTermStructureSpec<T: Real = f64> {
        pub garch_intercept_omega: T,
        pub arch_shock_alpha: T,
        pub garch_persistence_beta: T,
        pub current_daily_variance_h1: T,
        pub innovation_kurtosis_kappa_z: T,
        pub forecast_horizon_days_k: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BollerslevTsayGarchTermStructureReport<T: Real = f64> {
        pub total_persistence_phi: T,
        pub unconditional_daily_variance: T,
        pub k_step_ahead_daily_variance: T,
        pub annualized_term_structure_volatility: T,
        pub volatility_shock_half_life_days: T,
        pub unconditional_return_excess_kurtosis: T,
    }

    pub fn evaluate_bollerslev_tsay_garch_term_structure_kurtosis<T: Real>(
        spec: &BollerslevTsayGarchTermStructureSpec<T>,
    ) -> BollerslevTsayGarchTermStructureReport<T> {
        let omega = spec.garch_intercept_omega.max_val(T::from_f64(1e-8));
        let alpha = spec
            .arch_shock_alpha
            .max_val(T::from_f64(1e-3))
            .min_val(T::from_f64(0.30));
        let beta = spec
            .garch_persistence_beta
            .max_val(T::from_f64(0.10))
            .min_val(T::from_f64(0.95));
        let h1 = spec.current_daily_variance_h1.max_val(T::from_f64(1e-6));
        let kappa_z = spec.innovation_kurtosis_kappa_z.max_val(T::from_f64(3.0));
        let k_days = spec.forecast_horizon_days_k.max_val(T::one());

        let phi = (alpha + beta).min_val(T::from_f64(0.995));
        let one_m_phi = T::one() - phi;
        let h_inf = omega / one_m_phi;

        let phi_km1 = phi.pow(k_days - T::one());
        let phi_k = phi.pow(k_days);
        let h_k = h_inf + phi_km1 * (h1 - h_inf);
        let h_avg = h_inf + (h1 - h_inf) * (T::one() - phi_k) / (k_days * one_m_phi);
        let sig_ann = (T::from_f64(252.0) * h_avg).sqrt();

        let half_life = T::from_f64(2.0_f64.ln()) / (-phi.ln()).max_val(T::from_f64(1e-6));

        // Bollerslev (1986, Theorem 1) unconditional fourth moment kurtosis:
        // kappa_r = kappa_z * (1 - phi^2) / (1 - phi^2 - (kappa_z - 1) * alpha^2)
        let one_m_phi2 = T::one() - phi * phi;
        let denom = (one_m_phi2 - (kappa_z - T::one()) * alpha * alpha).max_val(T::from_f64(1e-4));
        let kappa_r = kappa_z * one_m_phi2 / denom;
        let ex_kurt = kappa_r - T::from_f64(3.0);

        BollerslevTsayGarchTermStructureReport {
            total_persistence_phi: phi,
            unconditional_daily_variance: h_inf,
            k_step_ahead_daily_variance: h_k,
            annualized_term_structure_volatility: sig_ann,
            volatility_shock_half_life_days: half_life,
            unconditional_return_excess_kurtosis: ex_kurt,
        }
    }

    #[cfg(test)]
    mod evaluate_bollerslev_tsay_garch_term_structure_kurtosis_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bollerslev_tsay_garch_term_structure_kurtosis_invariants_and_ad_parity() {
            let spec = BollerslevTsayGarchTermStructureSpec {
                garch_intercept_omega: 0.000_004,
                arch_shock_alpha: 0.085,
                garch_persistence_beta: 0.895,
                current_daily_variance_h1: 0.000_32,
                innovation_kurtosis_kappa_z: 3.50,
                forecast_horizon_days_k: 21.0,
            };
            let rep = evaluate_bollerslev_tsay_garch_term_structure_kurtosis(&spec);
            assert!(rep.annualized_term_structure_volatility > 0.20);
            assert!(rep.volatility_shock_half_life_days > 30.0);
            assert!(rep.unconditional_return_excess_kurtosis > 0.50);

            let h = 1e-6;
            let r_ad = evaluate_bollerslev_tsay_garch_term_structure_kurtosis(
                &BollerslevTsayGarchTermStructureSpec {
                    garch_intercept_omega: Dual::<6>::variable(spec.garch_intercept_omega, 0),
                    arch_shock_alpha: Dual::<6>::variable(spec.arch_shock_alpha, 1),
                    garch_persistence_beta: Dual::<6>::variable(spec.garch_persistence_beta, 2),
                    current_daily_variance_h1: Dual::<6>::variable(
                        spec.current_daily_variance_h1,
                        3,
                    ),
                    innovation_kurtosis_kappa_z: Dual::<6>::variable(
                        spec.innovation_kurtosis_kappa_z,
                        4,
                    ),
                    forecast_horizon_days_k: Dual::<6>::variable(spec.forecast_horizon_days_k, 5),
                },
            );
            let comp_ad = r_ad.annualized_term_structure_volatility * Dual::from_f64(10.0)
                + r_ad.unconditional_return_excess_kurtosis;
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                let step = if idx == 0 { 1e-8 } else { h };
                match idx {
                    0 => {
                        sp.garch_intercept_omega += step;
                        sm.garch_intercept_omega -= step;
                    }
                    1 => {
                        sp.arch_shock_alpha += step;
                        sm.arch_shock_alpha -= step;
                    }
                    2 => {
                        sp.garch_persistence_beta += step;
                        sm.garch_persistence_beta -= step;
                    }
                    3 => {
                        sp.current_daily_variance_h1 += step;
                        sm.current_daily_variance_h1 -= step;
                    }
                    4 => {
                        sp.innovation_kurtosis_kappa_z += step;
                        sm.innovation_kurtosis_kappa_z -= step;
                    }
                    _ => {
                        sp.forecast_horizon_days_k += step;
                        sm.forecast_horizon_days_k -= step;
                    }
                }
                let rp = evaluate_bollerslev_tsay_garch_term_structure_kurtosis(&sp);
                let rm = evaluate_bollerslev_tsay_garch_term_structure_kurtosis(&sm);
                let fd = ((rp.annualized_term_structure_volatility * 10.0
                    + rp.unconditional_return_excess_kurtosis)
                    - (rm.annualized_term_structure_volatility * 10.0
                        + rm.unconditional_return_excess_kurtosis))
                    / (2.0 * step);
                assert!(ad.abs() > 1e-5);
                assert!((ad - fd).abs() / ad.abs().max(1.0) < 1e-3);
            }
        }
    }
}
pub use evaluate_bollerslev_tsay_garch_term_structure_kurtosis_module::*;

mod evaluate_corsi_har_realized_volatility_multi_horizon_forecast_module {
    // Corsi (2009, Journal of Financial Econometrics 7(2), "A Simple Approximate Long-Memory
    // Model of Realized Volatility") & Andersen, Bollerslev, Diebold & Labys (2003, Econometrica 71(2))
    // 3-Cascade HAR-RV Multi-Horizon Term Structure & Realized Quarticity Band Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct CorsiHarRealizedVolatilitySpec<T: Real = f64> {
        pub har_intercept_c: T,
        pub daily_rv_weight_beta_d: T,
        pub weekly_rv_weight_beta_w: T,
        pub monthly_rv_weight_beta_m: T,
        pub current_daily_rv: T,
        pub current_weekly_rv: T,
        pub current_monthly_rv: T,
        pub realized_quarticity_rq: T,
        pub forecast_horizon_days_h: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct CorsiHarRealizedVolatilityReport<T: Real = f64> {
        pub har_total_persistence_phi: T,
        pub unconditional_long_run_rv: T,
        pub one_day_ahead_har_rv_forecast: T,
        pub h_day_ahead_har_rv_forecast: T,
        pub cumulative_average_h_day_rv: T,
        pub bns_asymptotic_rv_standard_error: T,
    }

    pub fn evaluate_corsi_har_realized_volatility_multi_horizon_forecast<T: Real>(
        spec: &CorsiHarRealizedVolatilitySpec<T>,
    ) -> CorsiHarRealizedVolatilityReport<T> {
        let c = spec.har_intercept_c.max_val(T::from_f64(1e-6));
        let bd = spec.daily_rv_weight_beta_d.max_val(T::from_f64(0.01));
        let bw = spec.weekly_rv_weight_beta_w.max_val(T::from_f64(0.01));
        let bm = spec.monthly_rv_weight_beta_m.max_val(T::from_f64(0.01));
        let rv_d = spec.current_daily_rv.max_val(T::from_f64(1e-6));
        let rv_w = spec.current_weekly_rv.max_val(T::from_f64(1e-6));
        let rv_m = spec.current_monthly_rv.max_val(T::from_f64(1e-6));
        let rq = spec.realized_quarticity_rq.max_val(T::from_f64(1e-8));
        let h = spec.forecast_horizon_days_h.max_val(T::one());

        let phi = (bd + bw + bm).min_val(T::from_f64(0.995));
        let one_minus_phi = (T::one() - phi).max_val(T::from_f64(0.005));
        let rv_inf = c / one_minus_phi;

        let rv_1 = c + bd * rv_d + bw * rv_w + bm * rv_m;
        let ln_phi = phi.ln();
        let phi_hm1 = ((h - T::one()) * ln_phi).exp();
        let phi_h = (h * ln_phi).exp();

        let rv_h = rv_inf + phi_hm1 * (rv_1 - rv_inf);
        let rv_avg = rv_inf + (rv_1 - rv_inf) * ((T::one() - phi_h) / (h * one_minus_phi));

        // Barndorff-Nielsen & Shephard (2002) asymptotic standard error of RV with M = 78 (5-min) intraday returns:
        let se_rv = (T::from_f64(2.0 / 78.0) * rq).sqrt();

        CorsiHarRealizedVolatilityReport {
            har_total_persistence_phi: phi,
            unconditional_long_run_rv: rv_inf,
            one_day_ahead_har_rv_forecast: rv_1,
            h_day_ahead_har_rv_forecast: rv_h,
            cumulative_average_h_day_rv: rv_avg,
            bns_asymptotic_rv_standard_error: se_rv,
        }
    }

    #[cfg(test)]
    mod evaluate_corsi_har_realized_volatility_multi_horizon_forecast_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn corsi_har_realized_volatility_invariants_and_ad_parity() {
            let spec = CorsiHarRealizedVolatilitySpec {
                har_intercept_c: 0.000_012,
                daily_rv_weight_beta_d: 0.36,
                weekly_rv_weight_beta_w: 0.34,
                monthly_rv_weight_beta_m: 0.25,
                current_daily_rv: 0.000_42,
                current_weekly_rv: 0.000_31,
                current_monthly_rv: 0.000_26,
                realized_quarticity_rq: 0.000_000_45,
                forecast_horizon_days_h: 22.0,
            };
            let rep = evaluate_corsi_har_realized_volatility_multi_horizon_forecast(&spec);
            assert!((rep.har_total_persistence_phi - 0.95).abs() < 1e-12);
            assert!(
                rep.h_day_ahead_har_rv_forecast < rep.one_day_ahead_har_rv_forecast
                    && rep.h_day_ahead_har_rv_forecast > rep.unconditional_long_run_rv
            );

            let h = 1e-6;
            let r_ad = evaluate_corsi_har_realized_volatility_multi_horizon_forecast(
                &CorsiHarRealizedVolatilitySpec {
                    har_intercept_c: Dual::from_f64(spec.har_intercept_c),
                    daily_rv_weight_beta_d: Dual::<6>::variable(spec.daily_rv_weight_beta_d, 0),
                    weekly_rv_weight_beta_w: Dual::<6>::variable(spec.weekly_rv_weight_beta_w, 1),
                    monthly_rv_weight_beta_m: Dual::<6>::variable(spec.monthly_rv_weight_beta_m, 2),
                    current_daily_rv: Dual::<6>::variable(spec.current_daily_rv, 3),
                    current_weekly_rv: Dual::<6>::variable(spec.current_weekly_rv, 4),
                    current_monthly_rv: Dual::from_f64(spec.current_monthly_rv),
                    realized_quarticity_rq: Dual::<6>::variable(spec.realized_quarticity_rq, 5),
                    forecast_horizon_days_h: Dual::from_f64(spec.forecast_horizon_days_h),
                },
            );
            for i in 0..6 {
                // Keep the quarticity bump below its positive floor so the
                // central difference remains in the smooth branch.
                let step = if i == 5 { 1e-10 } else { h };
                let mut v = [
                    spec.daily_rv_weight_beta_d,
                    spec.weekly_rv_weight_beta_w,
                    spec.monthly_rv_weight_beta_m,
                    spec.current_daily_rv,
                    spec.current_weekly_rv,
                    spec.realized_quarticity_rq,
                ];
                v[i] += step;
                let up = evaluate_corsi_har_realized_volatility_multi_horizon_forecast(
                    &CorsiHarRealizedVolatilitySpec {
                        daily_rv_weight_beta_d: v[0],
                        weekly_rv_weight_beta_w: v[1],
                        monthly_rv_weight_beta_m: v[2],
                        current_daily_rv: v[3],
                        current_weekly_rv: v[4],
                        realized_quarticity_rq: v[5],
                        ..spec
                    },
                );
                v[i] -= 2.0 * step;
                let dn = evaluate_corsi_har_realized_volatility_multi_horizon_forecast(
                    &CorsiHarRealizedVolatilitySpec {
                        daily_rv_weight_beta_d: v[0],
                        weekly_rv_weight_beta_w: v[1],
                        monthly_rv_weight_beta_m: v[2],
                        current_daily_rv: v[3],
                        current_weekly_rv: v[4],
                        realized_quarticity_rq: v[5],
                        ..spec
                    },
                );
                let val_ad =
                    r_ad.cumulative_average_h_day_rv + r_ad.bns_asymptotic_rv_standard_error;
                let fd = ((up.cumulative_average_h_day_rv + up.bns_asymptotic_rv_standard_error)
                    - (dn.cumulative_average_h_day_rv + dn.bns_asymptotic_rv_standard_error))
                    / (2.0 * step);
                assert!((val_ad.d[i] - fd).abs() < 1e-5, "input {i}");
            }
        }
    }
}
pub use evaluate_corsi_har_realized_volatility_multi_horizon_forecast_module::*;

mod evaluate_bai_ng_approximate_factor_icp2_diffusion_index_module {
    // Bai & Ng (2002, Econometrica 70(1), "Determining the Number of Factors in Approximate
    // Factor Models", Eq. 9) & Stock & Watson (2002, JASA 97(460), "Forecasting Using Principal
    // Components from a Large Number of Predictors")
    // Approximate Factor Model IC_p2 Criterion & Diffusion-Index Attenuation Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BaiNgApproximateFactorIcp2Spec<T: Real = f64> {
        pub cross_section_dimension_n: T,
        pub time_series_dimension_t: T,
        pub candidate_factor_count_k: T,
        pub total_panel_variance_sigma_x2: T,
        pub explained_variance_fraction_r2_k: T,
        pub diffusion_forecast_beta_f: T,
        pub factor_shock_signal_s_f: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BaiNgApproximateFactorIcp2Report<T: Real = f64> {
        pub residual_idiosyncratic_mse_v_k: T,
        pub bai_ng_complexity_penalty: T,
        pub bai_ng_ic_p2_criterion: T,
        pub bai_ng_pc_p2_criterion: T,
        pub finite_n_factor_reliability_ratio: T,
        pub attenuated_diffusion_index_forecast: T,
    }

    pub fn evaluate_bai_ng_approximate_factor_icp2_diffusion_index<T: Real>(
        spec: &BaiNgApproximateFactorIcp2Spec<T>,
    ) -> BaiNgApproximateFactorIcp2Report<T> {
        let n = spec.cross_section_dimension_n.max_val(T::from_f64(10.0));
        let t_dim = spec.time_series_dimension_t.max_val(T::from_f64(10.0));
        let k = spec.candidate_factor_count_k.max_val(T::one());
        let sig_x2 = spec
            .total_panel_variance_sigma_x2
            .max_val(T::from_f64(1e-4));
        let r2_k = spec
            .explained_variance_fraction_r2_k
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.99));
        let beta_f = spec.diffusion_forecast_beta_f;
        let s_f = spec.factor_shock_signal_s_f;

        let v_k = (T::one() - r2_k) * sig_x2;

        // Smooth C_NT^2 = min(N, T) for Bai-Ng (2002, Eq. 9) penalty:
        let diff_nt = n - t_dim;
        let c_nt2 = T::from_f64(0.5) * (n + t_dim - (diff_nt * diff_nt + T::from_f64(1.0)).sqrt());
        let penalty_rate = ((n + t_dim) / (n * t_dim)) * c_nt2.ln();
        let penalty = k * penalty_rate;

        let ic_p2 = v_k.ln() + penalty;
        let pc_p2 = v_k + sig_x2 * penalty;

        // Stock & Watson (2002) finite-N signal-to-noise reliability of PCA factor estimate:
        let sig_f2 = (r2_k / k) * sig_x2;
        let reliability = (n * sig_f2) / (n * sig_f2 + v_k);
        let forecast = reliability * beta_f * s_f;

        BaiNgApproximateFactorIcp2Report {
            residual_idiosyncratic_mse_v_k: v_k,
            bai_ng_complexity_penalty: penalty,
            bai_ng_ic_p2_criterion: ic_p2,
            bai_ng_pc_p2_criterion: pc_p2,
            finite_n_factor_reliability_ratio: reliability,
            attenuated_diffusion_index_forecast: forecast,
        }
    }

    #[cfg(test)]
    mod evaluate_bai_ng_approximate_factor_icp2_diffusion_index_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bai_ng_approximate_factor_icp2_invariants_and_ad_parity() {
            let spec = BaiNgApproximateFactorIcp2Spec {
                cross_section_dimension_n: 120.0,
                time_series_dimension_t: 240.0,
                candidate_factor_count_k: 4.0,
                total_panel_variance_sigma_x2: 1.0,
                explained_variance_fraction_r2_k: 0.58,
                diffusion_forecast_beta_f: 0.35,
                factor_shock_signal_s_f: 1.40,
            };
            let rep = evaluate_bai_ng_approximate_factor_icp2_diffusion_index(&spec);
            assert!(rep.bai_ng_complexity_penalty > 0.0);
            assert!(
                rep.finite_n_factor_reliability_ratio > 0.90
                    && rep.finite_n_factor_reliability_ratio < 1.0
            );

            let h = 1e-5;
            let r_ad = evaluate_bai_ng_approximate_factor_icp2_diffusion_index(
                &BaiNgApproximateFactorIcp2Spec {
                    cross_section_dimension_n: Dual::<7>::variable(
                        spec.cross_section_dimension_n,
                        0,
                    ),
                    time_series_dimension_t: Dual::<7>::variable(spec.time_series_dimension_t, 1),
                    candidate_factor_count_k: Dual::<7>::variable(spec.candidate_factor_count_k, 2),
                    total_panel_variance_sigma_x2: Dual::<7>::variable(
                        spec.total_panel_variance_sigma_x2,
                        3,
                    ),
                    explained_variance_fraction_r2_k: Dual::<7>::variable(
                        spec.explained_variance_fraction_r2_k,
                        4,
                    ),
                    diffusion_forecast_beta_f: Dual::<7>::variable(
                        spec.diffusion_forecast_beta_f,
                        5,
                    ),
                    factor_shock_signal_s_f: Dual::<7>::variable(spec.factor_shock_signal_s_f, 6),
                },
            );
            for i in 0..7 {
                let mut v = [
                    spec.cross_section_dimension_n,
                    spec.time_series_dimension_t,
                    spec.candidate_factor_count_k,
                    spec.total_panel_variance_sigma_x2,
                    spec.explained_variance_fraction_r2_k,
                    spec.diffusion_forecast_beta_f,
                    spec.factor_shock_signal_s_f,
                ];
                v[i] += h;
                let up = evaluate_bai_ng_approximate_factor_icp2_diffusion_index(
                    &BaiNgApproximateFactorIcp2Spec {
                        cross_section_dimension_n: v[0],
                        time_series_dimension_t: v[1],
                        candidate_factor_count_k: v[2],
                        total_panel_variance_sigma_x2: v[3],
                        explained_variance_fraction_r2_k: v[4],
                        diffusion_forecast_beta_f: v[5],
                        factor_shock_signal_s_f: v[6],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_bai_ng_approximate_factor_icp2_diffusion_index(
                    &BaiNgApproximateFactorIcp2Spec {
                        cross_section_dimension_n: v[0],
                        time_series_dimension_t: v[1],
                        candidate_factor_count_k: v[2],
                        total_panel_variance_sigma_x2: v[3],
                        explained_variance_fraction_r2_k: v[4],
                        diffusion_forecast_beta_f: v[5],
                        factor_shock_signal_s_f: v[6],
                    },
                );
                let val_ad = r_ad.bai_ng_pc_p2_criterion + r_ad.attenuated_diffusion_index_forecast;
                let fd = ((up.bai_ng_pc_p2_criterion + up.attenuated_diffusion_index_forecast)
                    - (dn.bai_ng_pc_p2_criterion + dn.attenuated_diffusion_index_forecast))
                    / (2.0 * h);
                assert!((val_ad.d[i] - fd).abs() < 1e-5, "input {i}");
            }
        }
    }
}
pub use evaluate_bai_ng_approximate_factor_icp2_diffusion_index_module::*;

mod evaluate_barndorff_nielsen_shephard_realized_quarticity_iv_ci_module {
    // Barndorff-Nielsen & Shephard (2002, JRSS Series B 64(2), "Econometric Analysis of
    // Realized Volatility and Its Use in Estimating Stochastic Volatility Models", Section 2.4),
    // Andersen, Bollerslev, Diebold & Labys (2003, Econometrica 71(2)), & Hansen & Lunde (2006)
    // Realized Quarticity Asymptotic CLT & Log-Transformed Integrated Variance Confidence Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BarndorffNielsenShephardRealizedQuarticitySpec<T: Real = f64> {
        pub realized_variance_rv: T,
        pub realized_quarticity_rq: T,
        pub intraday_sampling_frequency_m: T,
        pub first_order_autocovariance_gamma1: T,
        pub normal_critical_z_alpha: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BarndorffNielsenShephardRealizedQuarticityReport<T: Real = f64> {
        pub asymptotic_variance_of_rv: T,
        pub log_rv_asymptotic_standard_error: T,
        pub log_transformed_iv_lower_bound: T,
        pub log_transformed_iv_upper_bound: T,
        pub newey_west_microstructure_corrected_rv: T,
        pub relative_quarticity_dispersion_ratio: T,
    }

    pub fn evaluate_barndorff_nielsen_shephard_realized_quarticity_iv_ci<T: Real>(
        spec: &BarndorffNielsenShephardRealizedQuarticitySpec<T>,
    ) -> BarndorffNielsenShephardRealizedQuarticityReport<T> {
        let rv = spec.realized_variance_rv.max_val(T::from_f64(1e-6));
        let rq = spec.realized_quarticity_rq.max_val(T::from_f64(1e-9));
        let m = spec.intraday_sampling_frequency_m.max_val(T::from_f64(2.0));
        let gamma1 = spec.first_order_autocovariance_gamma1;
        let z_crit = spec.normal_critical_z_alpha.max_val(T::from_f64(0.5));

        // Barndorff-Nielsen & Shephard (2002, Theorem 1 & Eq. 11):
        // AsymVar(RV_M) = 2 * RQ_M / M, and SE(ln RV_M) = sqrt(2 * RQ_M / (M * RV_M^2)):
        let asym_var = (T::from_f64(2.0) * rq) / m;
        let se_log_rv = (asym_var / (rv * rv)).sqrt();
        let half_width = z_crit * se_log_rv;

        let iv_lower = rv * (-half_width).exp();
        let iv_upper = rv * half_width.exp();

        // Hansen & Lunde (2006) lag-1 Newey-West microstructure-noise-adjusted realized variance:
        let rv_nw1 = (rv + T::from_f64(2.0) * gamma1).max_val(T::from_f64(1e-6));
        let rq_over_rv2 = rq / (rv * rv);

        BarndorffNielsenShephardRealizedQuarticityReport {
            asymptotic_variance_of_rv: asym_var,
            log_rv_asymptotic_standard_error: se_log_rv,
            log_transformed_iv_lower_bound: iv_lower,
            log_transformed_iv_upper_bound: iv_upper,
            newey_west_microstructure_corrected_rv: rv_nw1,
            relative_quarticity_dispersion_ratio: rq_over_rv2,
        }
    }

    #[cfg(test)]
    mod evaluate_barndorff_nielsen_shephard_realized_quarticity_iv_ci_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn barndorff_nielsen_shephard_realized_quarticity_invariants_and_ad_parity() {
            let spec = BarndorffNielsenShephardRealizedQuarticitySpec {
                realized_variance_rv: 0.032,
                realized_quarticity_rq: 0.0014,
                intraday_sampling_frequency_m: 78.0,
                first_order_autocovariance_gamma1: -0.0008,
                normal_critical_z_alpha: 1.96,
            };
            let rep = evaluate_barndorff_nielsen_shephard_realized_quarticity_iv_ci(&spec);
            assert!(
                rep.log_transformed_iv_lower_bound < spec.realized_variance_rv
                    && rep.log_transformed_iv_upper_bound > spec.realized_variance_rv
            );

            let h = 1e-6;
            let r_ad = evaluate_barndorff_nielsen_shephard_realized_quarticity_iv_ci(
                &BarndorffNielsenShephardRealizedQuarticitySpec {
                    realized_variance_rv: Dual::<5>::variable(spec.realized_variance_rv, 0),
                    realized_quarticity_rq: Dual::<5>::variable(spec.realized_quarticity_rq, 1),
                    intraday_sampling_frequency_m: Dual::<5>::variable(
                        spec.intraday_sampling_frequency_m,
                        2,
                    ),
                    first_order_autocovariance_gamma1: Dual::<5>::variable(
                        spec.first_order_autocovariance_gamma1,
                        3,
                    ),
                    normal_critical_z_alpha: Dual::<5>::variable(spec.normal_critical_z_alpha, 4),
                },
            );
            for i in 0..5 {
                let mut v = [
                    spec.realized_variance_rv,
                    spec.realized_quarticity_rq,
                    spec.intraday_sampling_frequency_m,
                    spec.first_order_autocovariance_gamma1,
                    spec.normal_critical_z_alpha,
                ];
                v[i] += h;
                let up = evaluate_barndorff_nielsen_shephard_realized_quarticity_iv_ci(
                    &BarndorffNielsenShephardRealizedQuarticitySpec {
                        realized_variance_rv: v[0],
                        realized_quarticity_rq: v[1],
                        intraday_sampling_frequency_m: v[2],
                        first_order_autocovariance_gamma1: v[3],
                        normal_critical_z_alpha: v[4],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_barndorff_nielsen_shephard_realized_quarticity_iv_ci(
                    &BarndorffNielsenShephardRealizedQuarticitySpec {
                        realized_variance_rv: v[0],
                        realized_quarticity_rq: v[1],
                        intraday_sampling_frequency_m: v[2],
                        first_order_autocovariance_gamma1: v[3],
                        normal_critical_z_alpha: v[4],
                    },
                );
                let val_ad = r_ad.log_transformed_iv_upper_bound
                    + r_ad.newey_west_microstructure_corrected_rv;
                let fd = ((up.log_transformed_iv_upper_bound
                    + up.newey_west_microstructure_corrected_rv)
                    - (dn.log_transformed_iv_upper_bound
                        + dn.newey_west_microstructure_corrected_rv))
                    / (2.0 * h);
                assert!((val_ad.d[i] - fd).abs() < 1e-4, "input {i}");
            }
        }
    }
}
pub use evaluate_barndorff_nielsen_shephard_realized_quarticity_iv_ci_module::*;

use ql_ad::Real;

/// Specification for the Engle–Kroner (1995) positive-definite bivariate BEKK(1,1,1) multivariate GARCH
/// covariance and volatility spillover engine.
///
/// Reference: Engle & Kroner (1995, *Econometric Theory* 11(1), "Multivariate Simultaneous
/// Generalized ARCH"), Bollerslev, Engle & Wooldridge (1988, *JPE*), Tsay (2010, Ch. 10).
#[derive(Debug, Clone, Copy)]
pub struct EngleKronerBekkGarchSpec<T: Real> {
    /// Upper-triangular intercept entries $(c_{11}, c_{12}, c_{22})$ of $C$.
    pub cholesky_intercept_c: [T; 3],
    /// ARCH shock matrix $A = [[a_{11}, a_{12}], [a_{21}, a_{22}]]$ including cross-market spillover $a_{12}, a_{21}$.
    pub arch_spillover_matrix_a: [T; 4],
    /// Diagonal GARCH persistence parameters $(b_{11}, b_{22})$ of $B$.
    pub garch_persistence_diag_b: [T; 2],
    /// Lagged bivariate return innovations $(\varepsilon_{1,t-1}, \varepsilon_{2,t-1})$.
    pub lagged_innovations_eps: [T; 2],
    /// Lagged conditional covariance entries $(h_{11,t-1}, h_{12,t-1}, h_{22,t-1})$.
    pub lagged_covariance_h: [T; 3],
}

/// Analytics returned by [`evaluate_engle_kroner_bekk_bivariate_garch_covariance`].
#[derive(Debug, Clone, Copy)]
pub struct EngleKronerBekkGarchAnalytics<T: Real> {
    /// Updated conditional variance of asset 1, $h_{11,t}$.
    pub conditional_variance_h11: T,
    /// Updated conditional covariance $h_{12,t}$.
    pub conditional_covariance_h12: T,
    /// Updated conditional correlation $\rho_{12,t} = h_{12,t} / \sqrt{h_{11,t} h_{22,t}}$.
    pub conditional_correlation_rho12: T,
    /// Conditional time-varying beta $\beta_{1|2,t} = h_{12,t} / h_{22,t}$ plus persistence trace.
    pub conditional_beta_and_persistence_metric: T,
}

/// Evaluates the Engle–Kroner (1995) symmetric positive-definite BEKK(1,1,1) recursion
/// $H_t = C^\top C + A^\top \varepsilon_{t-1} \varepsilon_{t-1}^\top A + B^\top H_{t-1} B$.
pub fn evaluate_engle_kroner_bekk_bivariate_garch_covariance<T: Real>(
    spec: &EngleKronerBekkGarchSpec<T>,
) -> EngleKronerBekkGarchAnalytics<T> {
    let c11 = spec.cholesky_intercept_c[0];
    let c12 = spec.cholesky_intercept_c[1];
    let c22 = spec.cholesky_intercept_c[2];

    // C^T C where C = [[c11, c12], [0, c22]]
    let cc11 = c11 * c11;
    let cc12 = c11 * c12;
    let cc22 = c12 * c12 + c22 * c22;

    // v = A^T eps_{t-1}, so A^T eps eps^T A = v v^T
    let a11 = spec.arch_spillover_matrix_a[0];
    let a12 = spec.arch_spillover_matrix_a[1];
    let a21 = spec.arch_spillover_matrix_a[2];
    let a22 = spec.arch_spillover_matrix_a[3];
    let e1 = spec.lagged_innovations_eps[0];
    let e2 = spec.lagged_innovations_eps[1];

    let v1 = a11 * e1 + a21 * e2;
    let v2 = a12 * e1 + a22 * e2;

    // B^T H_{t-1} B for diagonal B = diag(b11, b22)
    let b11 = spec.garch_persistence_diag_b[0];
    let b22 = spec.garch_persistence_diag_b[1];
    let h11_prev = spec.lagged_covariance_h[0];
    let h12_prev = spec.lagged_covariance_h[1];
    let h22_prev = spec.lagged_covariance_h[2];

    let h11 = cc11 + v1 * v1 + b11 * b11 * h11_prev;
    let h12 = cc12 + v1 * v2 + b11 * b22 * h12_prev;
    let h22 = cc22 + v2 * v2 + b22 * b22 * h22_prev;

    let denom = (h11 * h22).max_val(T::from_f64(1e-12)).sqrt();
    let rho12 = h12 / denom;
    let beta_1_on_2 = h12 / h22.max_val(T::from_f64(1e-8));
    let persistence_trace = a11 * a11 + a12 * a12 + a21 * a21 + a22 * a22 + b11 * b11 + b22 * b22;

    EngleKronerBekkGarchAnalytics {
        conditional_variance_h11: h11,
        conditional_covariance_h12: h12,
        conditional_correlation_rho12: rho12,
        conditional_beta_and_persistence_metric: beta_1_on_2 + T::from_f64(0.1) * persistence_trace,
    }
}

#[cfg(test)]
mod engle_kroner_bekk_suite64_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn engle_kroner_bekk_invariants_and_ad_parity() {
        let base = EngleKronerBekkGarchSpec {
            cholesky_intercept_c: [0.04_f64, 0.015, 0.035],
            arch_spillover_matrix_a: [0.28, 0.08, 0.06, 0.25],
            garch_persistence_diag_b: [0.92, 0.90],
            lagged_innovations_eps: [-0.025, -0.018],
            lagged_covariance_h: [0.032, 0.018, 0.028],
        };
        let r0 = evaluate_engle_kroner_bekk_bivariate_garch_covariance(&base);
        assert!(r0.conditional_variance_h11 > 0.0);
        assert!(r0.conditional_correlation_rho12 > -1.0 && r0.conditional_correlation_rho12 < 1.0);
        assert!(r0.conditional_beta_and_persistence_metric > 0.0);

        let inputs = [
            base.cholesky_intercept_c[0],
            base.cholesky_intercept_c[1],
            base.cholesky_intercept_c[2],
            base.arch_spillover_matrix_a[0],
            base.arch_spillover_matrix_a[1],
            base.garch_persistence_diag_b[0],
            base.lagged_innovations_eps[0],
            base.lagged_covariance_h[1],
        ];
        let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            dual_in[i] = Dual::variable(inputs[i], i);
        }
        let spec_ad = EngleKronerBekkGarchSpec {
            cholesky_intercept_c: [dual_in[0], dual_in[1], dual_in[2]],
            arch_spillover_matrix_a: [
                dual_in[3],
                dual_in[4],
                Dual::from_f64(base.arch_spillover_matrix_a[2]),
                Dual::from_f64(base.arch_spillover_matrix_a[3]),
            ],
            garch_persistence_diag_b: [
                dual_in[5],
                Dual::from_f64(base.garch_persistence_diag_b[1]),
            ],
            lagged_innovations_eps: [dual_in[6], Dual::from_f64(base.lagged_innovations_eps[1])],
            lagged_covariance_h: [
                Dual::from_f64(base.lagged_covariance_h[0]),
                dual_in[7],
                Dual::from_f64(base.lagged_covariance_h[2]),
            ],
        };
        let r_ad = evaluate_engle_kroner_bekk_bivariate_garch_covariance(&spec_ad);
        let eps = 1e-6;
        for i in 0..8 {
            let mut p = inputs;
            let mut m = inputs;
            p[i] += eps;
            m[i] -= eps;
            let mk = |v: [f64; 8]| EngleKronerBekkGarchSpec {
                cholesky_intercept_c: [v[0], v[1], v[2]],
                arch_spillover_matrix_a: [
                    v[3],
                    v[4],
                    base.arch_spillover_matrix_a[2],
                    base.arch_spillover_matrix_a[3],
                ],
                garch_persistence_diag_b: [v[5], base.garch_persistence_diag_b[1]],
                lagged_innovations_eps: [v[6], base.lagged_innovations_eps[1]],
                lagged_covariance_h: [
                    base.lagged_covariance_h[0],
                    v[7],
                    base.lagged_covariance_h[2],
                ],
            };
            let rp = evaluate_engle_kroner_bekk_bivariate_garch_covariance(&mk(p));
            let rm = evaluate_engle_kroner_bekk_bivariate_garch_covariance(&mk(m));
            let fd =
                (rp.conditional_correlation_rho12 - rm.conditional_correlation_rho12) / (2.0 * eps);
            assert!(
                r_ad.conditional_correlation_rho12.d[i].abs() > 1e-5,
                "zero derivative at {i}"
            );
            assert!((r_ad.conditional_correlation_rho12.d[i] - fd).abs() < 1e-4);
        }
    }
}

pub mod hamilton_markov_switching_suite65_mod {
    use ql_ad::Real;

    /// Specification for the Hamilton (1989, 1994) two-state Markov-switching regime filter engine.
    ///
    /// Reference: Hamilton (1989, *Econometrica* 57(2), "A New Approach to the Economic Analysis
    /// of Nonstationary Time Series and the Business Cycle"; 1994, *Time Series Analysis*, Ch. 22).
    #[derive(Debug, Clone, Copy)]
    pub struct HamiltonMarkovSwitchingSpec<T: Real> {
        /// Current observed macro/return realization $y_t$.
        pub observed_realization_yt: T,
        /// Expansion regime ($s_t = 1$) conditional mean $\mu_1$.
        pub expansion_mean_mu1: T,
        /// Recession/crisis regime ($s_t = 2$) conditional mean $\mu_2$.
        pub recession_mean_mu2: T,
        /// Expansion regime standard deviation $\sigma_1$.
        pub expansion_volatility_sigma1: T,
        /// Recession/crisis regime standard deviation $\sigma_2$.
        pub recession_volatility_sigma2: T,
        /// Transition probability $p_{11} = \mathbb{P}(s_t = 1 \mid s_{t-1} = 1)$ of staying in expansion.
        pub stay_expansion_prob_p11: T,
        /// Transition probability $p_{22} = \mathbb{P}(s_t = 2 \mid s_{t-1} = 2)$ of staying in recession.
        pub stay_recession_prob_p22: T,
        /// Prior filtered probability $\xi_{1, t-1 \mid t-1} = \mathbb{P}(s_{t-1} = 1 \mid \mathcal{F}_{t-1})$.
        pub prior_expansion_prob_xi1: T,
    }

    /// Analytics returned by [`evaluate_hamilton_markov_switching_regime_filter_step`].
    #[derive(Debug, Clone, Copy)]
    pub struct HamiltonMarkovSwitchingAnalytics<T: Real> {
        /// One-step-ahead predictive expansion probability $\xi_{1, t \mid t-1}$.
        pub predictive_expansion_probability: T,
        /// Updated posterior filtered recession probability $\xi_{2, t \mid t} = \mathbb{P}(s_t = 2 \mid \mathcal{F}_t)$.
        pub filtered_recession_probability: T,
        /// Marginal log-likelihood contribution $\ell_t = \ln f(y_t \mid \mathcal{F}_{t-1})$.
        pub marginal_log_likelihood_increment: T,
        /// One-step-ahead regime-filtered conditional expectation $\mathbb{E}[y_{t+1} \mid \mathcal{F}_t]$.
        pub next_step_conditional_forecast: T,
    }

    /// Evaluates one exact Bayes–Hamilton (1989, Section 4) nonlinear Markov-switching filter step.
    pub fn evaluate_hamilton_markov_switching_regime_filter_step<T: Real>(
        spec: &HamiltonMarkovSwitchingSpec<T>,
    ) -> HamiltonMarkovSwitchingAnalytics<T> {
        let p11 = spec.stay_expansion_prob_p11;
        let p22 = spec.stay_recession_prob_p22;
        let xi1_prev = spec.prior_expansion_prob_xi1;
        let xi2_prev = T::one() - xi1_prev;

        // 1. Predictive regime probabilities xi_{j, t|t-1}
        let xi1_pred = p11 * xi1_prev + (T::one() - p22) * xi2_prev;
        let xi2_pred = T::one() - xi1_pred;

        // 2. Regime-conditional Gaussian densities eta_{j, t}
        let inv_sqrt_2pi = T::from_f64(0.3989422804014327);
        let sig1 = spec.expansion_volatility_sigma1.max_val(T::from_f64(1e-4));
        let sig2 = spec.recession_volatility_sigma2.max_val(T::from_f64(1e-4));
        let z1 = (spec.observed_realization_yt - spec.expansion_mean_mu1) / sig1;
        let z2 = (spec.observed_realization_yt - spec.recession_mean_mu2) / sig2;

        let eta1 = (inv_sqrt_2pi / sig1) * (-T::from_f64(0.5) * z1 * z1).exp();
        let eta2 = (inv_sqrt_2pi / sig2) * (-T::from_f64(0.5) * z2 * z2).exp();

        // 3. Marginal density f(y_t | F_{t-1}) and posterior filtered probabilities xi_{j, t|t}
        let f_yt = (xi1_pred * eta1 + xi2_pred * eta2).max_val(T::from_f64(1e-15));
        let log_lik = f_yt.ln();
        let xi1_filt = (xi1_pred * eta1) / f_yt;
        let xi2_filt = (xi2_pred * eta2) / f_yt;

        // 4. Next-step predictive probabilities xi_{j, t+1|t} and conditional mean E[y_{t+1} | F_t]
        let xi1_next = p11 * xi1_filt + (T::one() - p22) * xi2_filt;
        let xi2_next = T::one() - xi1_next;
        let forecast_next = xi1_next * spec.expansion_mean_mu1 + xi2_next * spec.recession_mean_mu2;

        HamiltonMarkovSwitchingAnalytics {
            predictive_expansion_probability: xi1_pred,
            filtered_recession_probability: xi2_filt,
            marginal_log_likelihood_increment: log_lik,
            next_step_conditional_forecast: forecast_next,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hamilton_markov_switching_invariants_and_ad_parity() {
            let base = HamiltonMarkovSwitchingSpec {
                observed_realization_yt: -0.018_f64,
                expansion_mean_mu1: 0.008,
                recession_mean_mu2: -0.022,
                expansion_volatility_sigma1: 0.012,
                recession_volatility_sigma2: 0.028,
                stay_expansion_prob_p11: 0.95,
                stay_recession_prob_p22: 0.85,
                prior_expansion_prob_xi1: 0.80,
            };
            let r0 = evaluate_hamilton_markov_switching_regime_filter_step(&base);
            assert!(r0.predictive_expansion_probability > 0.70);
            assert!(
                r0.filtered_recession_probability > 0.30 && r0.filtered_recession_probability < 1.0
            );
            assert!(r0.marginal_log_likelihood_increment.is_finite());

            let inputs = [
                base.observed_realization_yt,
                base.expansion_mean_mu1,
                base.recession_mean_mu2,
                base.expansion_volatility_sigma1,
                base.recession_volatility_sigma2,
                base.stay_expansion_prob_p11,
                base.stay_recession_prob_p22,
                base.prior_expansion_prob_xi1,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HamiltonMarkovSwitchingSpec {
                observed_realization_yt: dual_in[0],
                expansion_mean_mu1: dual_in[1],
                recession_mean_mu2: dual_in[2],
                expansion_volatility_sigma1: dual_in[3],
                recession_volatility_sigma2: dual_in[4],
                stay_expansion_prob_p11: dual_in[5],
                stay_recession_prob_p22: dual_in[6],
                prior_expansion_prob_xi1: dual_in[7],
            };
            let r_ad = evaluate_hamilton_markov_switching_regime_filter_step(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| HamiltonMarkovSwitchingSpec {
                    observed_realization_yt: v[0],
                    expansion_mean_mu1: v[1],
                    recession_mean_mu2: v[2],
                    expansion_volatility_sigma1: v[3],
                    recession_volatility_sigma2: v[4],
                    stay_expansion_prob_p11: v[5],
                    stay_recession_prob_p22: v[6],
                    prior_expansion_prob_xi1: v[7],
                };
                let rp = evaluate_hamilton_markov_switching_regime_filter_step(&mk(p));
                let rm = evaluate_hamilton_markov_switching_regime_filter_step(&mk(m));
                let fd = (rp.next_step_conditional_forecast - rm.next_step_conditional_forecast)
                    / (2.0 * eps);
                assert!(
                    r_ad.next_step_conditional_forecast.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.next_step_conditional_forecast.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use hamilton_markov_switching_suite65_mod::*;

pub mod hasbrouck_gonzalo_granger_info_share_mod {
    use ql_ad::Real;

    /// Specification for the Hasbrouck (1995, J. Finance) & Gonzalo-Granger (1995, JBES)
    /// cointegrated bivariate VECM price discovery Information Share and Component Share engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HasbrouckGonzaloGrangerSpec<T: Real> {
        /// Error-correction speed $\alpha_1 > 0$ of Market 1 (magnitude of downward adjustment to $P_1 - P_2$).
        pub venue1_error_correction_speed_alpha1: T,
        /// Error-correction speed $\alpha_2 > 0$ of Market 2 (magnitude of upward adjustment to $P_1 - P_2$).
        pub venue2_error_correction_speed_alpha2: T,
        /// Reduced-form price innovation volatility $\sigma_1 > 0$ of Market 1.
        pub venue1_innovation_volatility_sigma1: T,
        /// Reduced-form price innovation volatility $\sigma_2 > 0$ of Market 2.
        pub venue2_innovation_volatility_sigma2: T,
        /// Contemporaneous innovation correlation $\rho \in (-1, 1)$ between Markets 1 and 2.
        pub cross_venue_innovation_correlation_rho: T,
    }

    /// Analytical outputs for the Hasbrouck & Gonzalo-Granger price discovery engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HasbrouckGonzaloGrangerAnalytics<T: Real> {
        /// Gonzalo-Granger (1995) permanent-transitory Component Share $\psi_1 = \alpha_2 / (\alpha_1 + \alpha_2)$.
        pub gonzalo_granger_component_share_market1: T,
        /// Long-run efficient random-walk fundamental variance $V_{\text{eff}} = \psi^\top \Omega \psi$.
        pub long_run_efficient_price_variance: T,
        /// Hasbrouck (1995) Cholesky upper-bound Information Share $IS_1^{\text{upper}}$ for Market 1.
        pub hasbrouck_information_share_upper_bound: T,
        /// Hasbrouck (1995) Cholesky midpoint Information Share $IS_1^{\text{mid}} = \frac{1}{2}(IS_1^{\text{upper}} + IS_1^{\text{lower}})$.
        pub hasbrouck_midpoint_information_share: T,
    }

    /// Evaluates the Hasbrouck (1995) Information Shares and Gonzalo-Granger (1995) Component Share
    /// for two cointegrated trading venues.
    pub fn evaluate_hasbrouck_gonzalo_granger_price_discovery_information_share<T: Real>(
        spec: &HasbrouckGonzaloGrangerSpec<T>,
    ) -> HasbrouckGonzaloGrangerAnalytics<T> {
        let a1 = spec
            .venue1_error_correction_speed_alpha1
            .max_val(T::from_f64(1e-5));
        let a2 = spec
            .venue2_error_correction_speed_alpha2
            .max_val(T::from_f64(1e-5));
        let sig1 = spec
            .venue1_innovation_volatility_sigma1
            .max_val(T::from_f64(1e-5));
        let sig2 = spec
            .venue2_innovation_volatility_sigma2
            .max_val(T::from_f64(1e-5));
        let rho = spec
            .cross_venue_innovation_correlation_rho
            .max_val(T::from_f64(-0.98))
            .min_val(T::from_f64(0.98));

        let a_sum = a1 + a2;
        let psi1 = a2 / a_sum;
        let psi2 = a1 / a_sum;

        let v_eff = (psi1 * psi1 * sig1 * sig1
            + psi2 * psi2 * sig2 * sig2
            + T::from_f64(2.0) * rho * psi1 * psi2 * sig1 * sig2)
            .max_val(T::from_f64(1e-12));

        let num_upper = psi1 * sig1 + psi2 * rho * sig2;
        let is1_upper = (num_upper * num_upper) / v_eff;
        let is1_lower = (psi1 * psi1 * sig1 * sig1 * (T::one() - rho * rho)) / v_eff;
        let is1_mid = T::from_f64(0.5) * (is1_upper + is1_lower);

        HasbrouckGonzaloGrangerAnalytics {
            gonzalo_granger_component_share_market1: psi1,
            long_run_efficient_price_variance: v_eff,
            hasbrouck_information_share_upper_bound: is1_upper,
            hasbrouck_midpoint_information_share: is1_mid,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hasbrouck_gonzalo_granger_invariants_and_ad_parity() {
            let base = HasbrouckGonzaloGrangerSpec {
                venue1_error_correction_speed_alpha1: 0.08_f64,
                venue2_error_correction_speed_alpha2: 0.24,
                venue1_innovation_volatility_sigma1: 0.015,
                venue2_innovation_volatility_sigma2: 0.018,
                cross_venue_innovation_correlation_rho: 0.55,
            };
            let r0 = evaluate_hasbrouck_gonzalo_granger_price_discovery_information_share(&base);
            assert!(
                r0.gonzalo_granger_component_share_market1 > 0.50
                    && r0.gonzalo_granger_component_share_market1 < 1.0
            );
            assert!(
                r0.hasbrouck_information_share_upper_bound
                    > r0.hasbrouck_midpoint_information_share
            );
            assert!(
                r0.hasbrouck_midpoint_information_share > 0.0
                    && r0.hasbrouck_midpoint_information_share < 1.0
            );

            let inputs = [
                base.venue1_error_correction_speed_alpha1,
                base.venue2_error_correction_speed_alpha2,
                base.venue1_innovation_volatility_sigma1,
                base.venue2_innovation_volatility_sigma2,
                base.cross_venue_innovation_correlation_rho,
            ];
            let mut dual_in = [Dual::<5>::from_f64(0.0); 5];
            for i in 0..5 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HasbrouckGonzaloGrangerSpec {
                venue1_error_correction_speed_alpha1: dual_in[0],
                venue2_error_correction_speed_alpha2: dual_in[1],
                venue1_innovation_volatility_sigma1: dual_in[2],
                venue2_innovation_volatility_sigma2: dual_in[3],
                cross_venue_innovation_correlation_rho: dual_in[4],
            };
            let r_ad =
                evaluate_hasbrouck_gonzalo_granger_price_discovery_information_share(&spec_ad);
            let eps = 1e-6;
            for i in 0..5 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 5]| HasbrouckGonzaloGrangerSpec {
                    venue1_error_correction_speed_alpha1: v[0],
                    venue2_error_correction_speed_alpha2: v[1],
                    venue1_innovation_volatility_sigma1: v[2],
                    venue2_innovation_volatility_sigma2: v[3],
                    cross_venue_innovation_correlation_rho: v[4],
                };
                let rp =
                    evaluate_hasbrouck_gonzalo_granger_price_discovery_information_share(&mk(p));
                let rm =
                    evaluate_hasbrouck_gonzalo_granger_price_discovery_information_share(&mk(m));
                let fd = (rp.hasbrouck_midpoint_information_share
                    - rm.hasbrouck_midpoint_information_share)
                    / (2.0 * eps);
                assert!(
                    r_ad.hasbrouck_midpoint_information_share.d[i].abs() > 1e-5,
                    "zero derivative at {i}"
                );
                assert!((r_ad.hasbrouck_midpoint_information_share.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use hasbrouck_gonzalo_granger_info_share_mod::*;

pub mod hansen_lunde_overnight_realized_variance_mod {
    use ql_ad::Real;

    /// Specification for the Hansen & Lunde (2005, JBES 23(4)) / Andersen, Bollerslev, Diebold & Labys (2003)
    /// optimal minimum-MSE combination of overnight squared return and intraday realized variance.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HansenLundeOvernightRvSpec<T: Real> {
        /// Realized squared overnight close-to-open return $r_{\text{on}}^2 > 0$.
        pub realized_overnight_sq_return_rv1: T,
        /// Realized intraday open-to-close high-frequency variance $RV_{\text{intra}} > 0$.
        pub realized_intraday_variance_rv2: T,
        /// Unconditional mean $\mu_1 = \mathbb{E}[r_{\text{on}}^2] > 0$.
        pub unconditional_overnight_mean_mu1: T,
        /// Unconditional mean $\mu_2 = \mathbb{E}[RV_{\text{intra}}] > 0$ (including microstructure bias).
        pub unconditional_intraday_mean_mu2: T,
        /// Unconditional target full-day integrated variance $\mu_0 = \mathbb{E}[IV] > 0$.
        pub unconditional_full_day_iv_mu0: T,
        /// Variance $\eta_1^2 = \text{Var}(r_{\text{on}}^2) > 0$ of overnight squared return.
        pub overnight_estimator_variance_eta1_sq: T,
        /// Variance $\eta_2^2 = \text{Var}(RV_{\text{intra}}) > 0$ of intraday realized variance.
        pub intraday_estimator_variance_eta2_sq: T,
        /// Covariance $\eta_{12} = \text{Cov}(r_{\text{on}}^2, RV_{\text{intra}})$.
        pub estimator_covariance_eta12: T,
    }

    /// Analytical outputs for the Hansen-Lunde overnight + intraday realized variance engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HansenLundeOvernightRvAnalytics<T: Real> {
        /// Optimal Hansen-Lunde (2005, Theorem 1) weight $\omega_1^*$ on overnight squared return.
        pub optimal_overnight_weight_w1: T,
        /// Optimal Hansen-Lunde (2005, Theorem 1) weight $\omega_2^*$ on intraday realized variance.
        pub optimal_intraday_weight_w2: T,
        /// Relative importance factor $\phi^*$ minimizing conditional MSE.
        pub optimal_relative_importance_phi: T,
        /// Combined minimum-MSE full-day realized volatility $\sqrt{\omega_1^* r_{\text{on}}^2 + \omega_2^* RV_{\text{intra}}}$.
        pub hansen_lunde_combined_realized_volatility: T,
    }

    /// Evaluates the Hansen & Lunde (2005, Theorem 1) minimum-MSE linear combination of
    /// overnight squared return and noisy intraday realized variance.
    pub fn evaluate_hansen_lunde_overnight_intraday_realized_variance_combination<T: Real>(
        spec: &HansenLundeOvernightRvSpec<T>,
    ) -> HansenLundeOvernightRvAnalytics<T> {
        let mu1 = spec
            .unconditional_overnight_mean_mu1
            .max_val(T::from_f64(1e-8));
        let mu2 = spec
            .unconditional_intraday_mean_mu2
            .max_val(T::from_f64(1e-8));
        let mu0 = spec
            .unconditional_full_day_iv_mu0
            .max_val(T::from_f64(1e-8));
        let v1 = spec
            .overnight_estimator_variance_eta1_sq
            .max_val(T::from_f64(1e-12));
        let v2 = spec
            .intraday_estimator_variance_eta2_sq
            .max_val(T::from_f64(1e-12));
        let c12 = spec.estimator_covariance_eta12;

        // Hansen & Lunde (2005, Eq. 6): relative importance parameter phi*
        // phi* = (mu2^2 * v1 - mu1 * mu2 * c12) / (mu2^2 * v1 + mu1^2 * v2 - 2 * mu1 * mu2 * c12)
        let num_phi = mu2 * mu2 * v1 - mu1 * mu2 * c12;
        let den_phi = (mu2 * mu2 * v1 + mu1 * mu1 * v2 - T::from_f64(2.0) * mu1 * mu2 * c12)
            // The variance inputs are naturally on the 1e-8 scale, so the
            // quadratic denominator is around 1e-15 for realistic quotes.
            // Keep the numerical guard below that scale to preserve the
            // estimator-variance sensitivities in the AD domain.
            .max_val(T::from_f64(1e-18));
        let phi_star = (num_phi / den_phi)
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.99));

        // Optimal weights w1* = (1 - phi*) * (mu0 / mu1), w2* = phi* * (mu0 / mu2)
        // satisfying the unbiasedness constraint w1* * mu1 + w2* * mu2 = mu0
        let w1 = (T::one() - phi_star) * (mu0 / mu1);
        let w2 = phi_star * (mu0 / mu2);

        let rv_comb = (w1 * spec.realized_overnight_sq_return_rv1
            + w2 * spec.realized_intraday_variance_rv2)
            .max_val(T::from_f64(1e-10));

        HansenLundeOvernightRvAnalytics {
            optimal_overnight_weight_w1: w1,
            optimal_intraday_weight_w2: w2,
            optimal_relative_importance_phi: phi_star,
            hansen_lunde_combined_realized_volatility: rv_comb.sqrt(),
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hansen_lunde_overnight_rv_invariants_and_ad_parity() {
            let base = HansenLundeOvernightRvSpec {
                realized_overnight_sq_return_rv1: 0.00012_f64,
                realized_intraday_variance_rv2: 0.00038,
                unconditional_overnight_mean_mu1: 0.00010,
                unconditional_intraday_mean_mu2: 0.00035,
                unconditional_full_day_iv_mu0: 0.00042,
                overnight_estimator_variance_eta1_sq: 2.5e-8,
                intraday_estimator_variance_eta2_sq: 1.8e-8,
                estimator_covariance_eta12: 0.4e-8,
            };
            let r0 = evaluate_hansen_lunde_overnight_intraday_realized_variance_combination(&base);
            let unbiased_sum = r0.optimal_overnight_weight_w1
                * base.unconditional_overnight_mean_mu1
                + r0.optimal_intraday_weight_w2 * base.unconditional_intraday_mean_mu2;
            assert!((unbiased_sum - base.unconditional_full_day_iv_mu0).abs() < 1e-12);
            assert!(
                r0.optimal_relative_importance_phi > 0.0
                    && r0.optimal_relative_importance_phi < 1.0
            );
            assert!(r0.hansen_lunde_combined_realized_volatility > 0.0);

            let inputs = [
                base.realized_overnight_sq_return_rv1,
                base.realized_intraday_variance_rv2,
                base.unconditional_overnight_mean_mu1,
                base.unconditional_intraday_mean_mu2,
                base.unconditional_full_day_iv_mu0,
                base.overnight_estimator_variance_eta1_sq,
                base.intraday_estimator_variance_eta2_sq,
                base.estimator_covariance_eta12,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HansenLundeOvernightRvSpec {
                realized_overnight_sq_return_rv1: dual_in[0],
                realized_intraday_variance_rv2: dual_in[1],
                unconditional_overnight_mean_mu1: dual_in[2],
                unconditional_intraday_mean_mu2: dual_in[3],
                unconditional_full_day_iv_mu0: dual_in[4],
                overnight_estimator_variance_eta1_sq: dual_in[5],
                intraday_estimator_variance_eta2_sq: dual_in[6],
                estimator_covariance_eta12: dual_in[7],
            };
            let r_ad =
                evaluate_hansen_lunde_overnight_intraday_realized_variance_combination(&spec_ad);
            let eps = 1e-10;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| HansenLundeOvernightRvSpec {
                    realized_overnight_sq_return_rv1: v[0],
                    realized_intraday_variance_rv2: v[1],
                    unconditional_overnight_mean_mu1: v[2],
                    unconditional_intraday_mean_mu2: v[3],
                    unconditional_full_day_iv_mu0: v[4],
                    overnight_estimator_variance_eta1_sq: v[5],
                    intraday_estimator_variance_eta2_sq: v[6],
                    estimator_covariance_eta12: v[7],
                };
                let rp =
                    evaluate_hansen_lunde_overnight_intraday_realized_variance_combination(&mk(p));
                let rm =
                    evaluate_hansen_lunde_overnight_intraday_realized_variance_combination(&mk(m));
                let fd = (rp.hansen_lunde_combined_realized_volatility
                    - rm.hansen_lunde_combined_realized_volatility)
                    / (2.0 * eps);
                assert!(
                    r_ad.hansen_lunde_combined_realized_volatility.d[i].abs() > 1e-4,
                    "zero derivative at {i}"
                );
                let rel_err = (r_ad.hansen_lunde_combined_realized_volatility.d[i] - fd).abs()
                    / fd.abs().max(1.0);
                assert!(rel_err < 1e-4);
            }
        }
    }
}
pub use hansen_lunde_overnight_realized_variance_mod::*;

// ============================================================================
// Engle & Granger (1987), Johansen (1991) & Vidyamurthy (2004) VECM Cointegration Spread
// ============================================================================

pub mod engle_granger_johansen_cointegration_mod {
    use ql_ad::Real;

    /// Specification for the Engle & Granger (1987, *Econometrica* 55(2)), Johansen (1991,
    /// *Econometrica* 59(6)), and Vidyamurthy (2004, Ch. 4–6) VECM cointegration pairs engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct EngleGrangerJohansenCointegrationSpec<T: Real> {
        /// Current log price of target asset $y_t = \ln P_{Y,t}$.
        pub current_log_price_y: T,
        /// Current log price of cointegrating peer asset $x_t = \ln P_{X,t}$.
        pub current_log_price_x: T,
        /// Johansen cointegrating vector slope $\beta_{\text{coint}} > 0$.
        pub cointegration_hedge_ratio_beta: T,
        /// Long-run equilibrium spread intercept $\mu_{\text{eq}}$.
        pub equilibrium_intercept_mu: T,
        /// VECM error-correction adjustment speed $\alpha_{\text{ecm}} \in (0, 1)$.
        pub error_correction_speed_alpha: T,
        /// Residual innovation volatility $\sigma_\varepsilon > 0$.
        pub innovation_volatility_sigma_eps: T,
        /// Johansen squared canonical correlation $R_{\text{canon}}^2 \in (0, 1)$.
        pub canonical_correlation_r_sq: T,
    }

    /// Analytical outputs for the Engle-Granger-Johansen VECM cointegration engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct EngleGrangerJohansenCointegrationAnalytics<T: Real> {
        /// Current disequilibrium error-correction spread $z_t = y_t - \mu_{\text{eq}} - \beta_{\text{coint}} x_t$.
        pub disequilibrium_error_correction_spread: T,
        /// Standardized stationary spread z-score $s_t = z_t / \sigma_z$.
        pub standardized_stationary_z_score: T,
        /// Exact discrete-time mean-reversion half-life $t_{1/2} = \frac{\ln 2}{-\ln(1 - \alpha_{\text{ecm}})}$.
        pub mean_reversion_half_life_periods: T,
        /// Johansen-eigenvalue-weighted error-correction Sharpe signal.
        pub cointegration_confidence_weighted_sharpe: T,
    }

    /// Evaluates the Engle & Granger (1987) / Johansen (1991) error-correction spread,
    /// stationary z-score, half-life, and canonical-correlation-weighted Sharpe ratio.
    pub fn evaluate_engle_granger_johansen_cointegration_error_correction_pairs_spread<T: Real>(
        spec: &EngleGrangerJohansenCointegrationSpec<T>,
    ) -> EngleGrangerJohansenCointegrationAnalytics<T> {
        let alpha = spec
            .error_correction_speed_alpha
            .max_val(T::from_f64(1e-4))
            .min_val(T::from_f64(0.95));
        let sig_eps = spec
            .innovation_volatility_sigma_eps
            .max_val(T::from_f64(1e-5));
        let r_sq = spec
            .canonical_correlation_r_sq
            .max_val(T::from_f64(1e-4))
            .min_val(T::from_f64(0.99));

        let z_t = spec.current_log_price_y
            - spec.equilibrium_intercept_mu
            - spec.cointegration_hedge_ratio_beta * spec.current_log_price_x;

        // Stationary AR(1) variance of z_t under phi = 1 - alpha: Var(z) = sigma_eps^2 / (alpha * (2 - alpha))
        let var_z = (sig_eps * sig_eps) / (alpha * (T::from_f64(2.0) - alpha));
        let sig_z = var_z.sqrt();
        let z_score = z_t / sig_z;

        let kappa_ecm = -(T::one() - alpha).ln();
        let half_life = T::from_f64(std::f64::consts::LN_2) / kappa_ecm;

        // Johansen canonical eigenvalue statistic lambda_J = -ln(1 - R_canon^2)
        let johansen_lambda = -(T::one() - r_sq).ln();
        let smooth_abs_z = (z_t * z_t + T::from_f64(1e-8)).sqrt();
        let ecm_sharpe = (alpha * smooth_abs_z / sig_eps) * johansen_lambda.sqrt();

        EngleGrangerJohansenCointegrationAnalytics {
            disequilibrium_error_correction_spread: z_t,
            standardized_stationary_z_score: z_score,
            mean_reversion_half_life_periods: half_life,
            cointegration_confidence_weighted_sharpe: ecm_sharpe,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn engle_granger_johansen_cointegration_invariants_and_ad_parity() {
            let base = EngleGrangerJohansenCointegrationSpec {
                current_log_price_y: 4.65_f64,
                current_log_price_x: 4.52,
                cointegration_hedge_ratio_beta: 1.04,
                equilibrium_intercept_mu: -0.08,
                error_correction_speed_alpha: 0.28,
                innovation_volatility_sigma_eps: 0.032,
                canonical_correlation_r_sq: 0.24,
            };
            let r0 =
                evaluate_engle_granger_johansen_cointegration_error_correction_pairs_spread(&base);
            assert!(r0.disequilibrium_error_correction_spread.abs() > 0.0);
            assert!(r0.standardized_stationary_z_score.abs() > 0.0);
            assert!(r0.mean_reversion_half_life_periods > 0.0);
            assert!(r0.cointegration_confidence_weighted_sharpe > 0.0);

            let inputs = [
                base.current_log_price_y,
                base.current_log_price_x,
                base.cointegration_hedge_ratio_beta,
                base.equilibrium_intercept_mu,
                base.error_correction_speed_alpha,
                base.innovation_volatility_sigma_eps,
                base.canonical_correlation_r_sq,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = EngleGrangerJohansenCointegrationSpec {
                current_log_price_y: dual_in[0],
                current_log_price_x: dual_in[1],
                cointegration_hedge_ratio_beta: dual_in[2],
                equilibrium_intercept_mu: dual_in[3],
                error_correction_speed_alpha: dual_in[4],
                innovation_volatility_sigma_eps: dual_in[5],
                canonical_correlation_r_sq: dual_in[6],
            };
            let r_ad = evaluate_engle_granger_johansen_cointegration_error_correction_pairs_spread(
                &spec_ad,
            );
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| EngleGrangerJohansenCointegrationSpec {
                    current_log_price_y: v[0],
                    current_log_price_x: v[1],
                    cointegration_hedge_ratio_beta: v[2],
                    equilibrium_intercept_mu: v[3],
                    error_correction_speed_alpha: v[4],
                    innovation_volatility_sigma_eps: v[5],
                    canonical_correlation_r_sq: v[6],
                };
                let rp =
                    evaluate_engle_granger_johansen_cointegration_error_correction_pairs_spread(
                        &mk(p),
                    );
                let rm =
                    evaluate_engle_granger_johansen_cointegration_error_correction_pairs_spread(
                        &mk(m),
                    );
                let fd = (rp.cointegration_confidence_weighted_sharpe
                    - rm.cointegration_confidence_weighted_sharpe)
                    / (2.0 * eps);
                assert!(
                    r_ad.cointegration_confidence_weighted_sharpe.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.cointegration_confidence_weighted_sharpe.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use engle_granger_johansen_cointegration_mod::*;

// ============================================================================
// Diebold & Mariano (1995), West (1996) & Harvey-Leybourne-Newbold (1997) Forecast Test
// ============================================================================

pub mod diebold_mariano_west_hln_mod {
    use ql_ad::Real;

    /// Specification for the Diebold & Mariano (1995, *JBES* 13(3)), West (1996, *Econometrica*
    /// 64(5)), and Harvey, Leybourne & Newbold (1997, *IJF* 13(2)) forecast evaluation test.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct DieboldMarianoWestForecastSpec<T: Real> {
        /// Out-of-sample mean loss differential $\bar{d} = \frac{1}{P}\sum (L_{1,t} - L_{2,t})$.
        pub mean_loss_differential_d_bar: T,
        /// Lag-0 sample variance of the loss differential $\hat{\gamma}_0 > 0$.
        pub lag0_loss_diff_variance_gamma0: T,
        /// Lag-1 sample autocovariance of the loss differential $\hat{\gamma}_1$.
        pub lag1_loss_diff_autocovariance_gamma1: T,
        /// Multi-step forecast horizon $h \ge 1$.
        pub forecast_horizon_steps_h: T,
        /// Out-of-sample evaluation sample size $P > h$.
        pub out_of_sample_observations_p: T,
        /// In-sample estimation sample size $R > 0$ for West's (1996) parameter uncertainty factor.
        pub in_sample_estimation_observations_r: T,
        /// Asymptotic parameter-estimation uncertainty contribution $V_\theta \ge 0$.
        pub parameter_estimation_variance_v_theta: T,
    }

    /// Analytical outputs for the Diebold-Mariano-West-HLN forecast accuracy test engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct DieboldMarianoWestForecastAnalytics<T: Real> {
        /// Newey-West / Bartlett + West (1996) long-run loss-differential variance $\hat{\Omega}_{\text{DMW}}$.
        pub dmw_long_run_loss_differential_variance: T,
        /// Standard asymptotic Diebold-Mariano-West statistic $\text{DMW}$.
        pub asymptotic_dmw_t_statistic: T,
        /// Harvey-Leybourne-Newbold (1997) finite-sample bias correction factor $\kappa_{\text{HLN}}$.
        pub harvey_leybourne_newbold_correction_factor: T,
        /// Small-sample corrected HLN-DMW test statistic $\text{DMW}_{\text{HLN}}$.
        pub hln_adjusted_dmw_t_statistic: T,
    }

    /// Evaluates the Diebold & Mariano (1995), West (1996), and Harvey, Leybourne & Newbold (1997)
    /// out-of-sample predictive accuracy test statistic.
    pub fn evaluate_diebold_mariano_west_hln_forecast_comparison_test<T: Real>(
        spec: &DieboldMarianoWestForecastSpec<T>,
    ) -> DieboldMarianoWestForecastAnalytics<T> {
        let g0 = spec
            .lag0_loss_diff_variance_gamma0
            .max_val(T::from_f64(1e-6));
        let h = spec.forecast_horizon_steps_h.max_val(T::one());
        let p_obs = spec
            .out_of_sample_observations_p
            .max_val(h + T::from_f64(5.0));
        let r_obs = spec.in_sample_estimation_observations_r.max_val(T::one());

        // Bartlett kernel weight at lag 1 and West (1996, Table I) recursive estimation multiplier
        let w_bartlett = T::one() - T::one() / (h + T::one());
        let pi_ratio = p_obs / r_obs;
        let west_mult = T::from_f64(2.0) * (T::one() - (T::one() + pi_ratio).ln() / pi_ratio);

        let omega_dmw = (g0
            + T::from_f64(2.0) * w_bartlett * spec.lag1_loss_diff_autocovariance_gamma1
            + west_mult * spec.parameter_estimation_variance_v_theta)
            .max_val(T::from_f64(1e-6));

        let dmw_stat = spec.mean_loss_differential_d_bar / (omega_dmw / p_obs).sqrt();

        // Harvey, Leybourne & Newbold (1997, Eq. 5) small-sample degrees-of-freedom factor
        let hln_num = (p_obs + T::one() - T::from_f64(2.0) * h + h * (h - T::one()) / p_obs)
            .max_val(T::one());
        let kappa_hln = (hln_num / p_obs).sqrt();
        let hln_stat = kappa_hln * dmw_stat;

        DieboldMarianoWestForecastAnalytics {
            dmw_long_run_loss_differential_variance: omega_dmw,
            asymptotic_dmw_t_statistic: dmw_stat,
            harvey_leybourne_newbold_correction_factor: kappa_hln,
            hln_adjusted_dmw_t_statistic: hln_stat,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn diebold_mariano_west_hln_invariants_and_ad_parity() {
            let base = DieboldMarianoWestForecastSpec {
                mean_loss_differential_d_bar: 0.028_f64,
                lag0_loss_diff_variance_gamma0: 0.045,
                lag1_loss_diff_autocovariance_gamma1: 0.014,
                forecast_horizon_steps_h: 3.0,
                out_of_sample_observations_p: 120.0,
                in_sample_estimation_observations_r: 240.0,
                parameter_estimation_variance_v_theta: 0.008,
            };
            let r0 = evaluate_diebold_mariano_west_hln_forecast_comparison_test(&base);
            assert!(
                r0.dmw_long_run_loss_differential_variance > base.lag0_loss_diff_variance_gamma0
            );
            assert!(r0.asymptotic_dmw_t_statistic > 0.0);
            assert!(
                r0.harvey_leybourne_newbold_correction_factor > 0.0
                    && r0.harvey_leybourne_newbold_correction_factor < 1.0
            );
            assert!(r0.hln_adjusted_dmw_t_statistic > 0.0);

            let inputs = [
                base.mean_loss_differential_d_bar,
                base.lag0_loss_diff_variance_gamma0,
                base.lag1_loss_diff_autocovariance_gamma1,
                base.forecast_horizon_steps_h,
                base.out_of_sample_observations_p,
                base.in_sample_estimation_observations_r,
                base.parameter_estimation_variance_v_theta,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = DieboldMarianoWestForecastSpec {
                mean_loss_differential_d_bar: dual_in[0],
                lag0_loss_diff_variance_gamma0: dual_in[1],
                lag1_loss_diff_autocovariance_gamma1: dual_in[2],
                forecast_horizon_steps_h: dual_in[3],
                out_of_sample_observations_p: dual_in[4],
                in_sample_estimation_observations_r: dual_in[5],
                parameter_estimation_variance_v_theta: dual_in[6],
            };
            let r_ad = evaluate_diebold_mariano_west_hln_forecast_comparison_test(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| DieboldMarianoWestForecastSpec {
                    mean_loss_differential_d_bar: v[0],
                    lag0_loss_diff_variance_gamma0: v[1],
                    lag1_loss_diff_autocovariance_gamma1: v[2],
                    forecast_horizon_steps_h: v[3],
                    out_of_sample_observations_p: v[4],
                    in_sample_estimation_observations_r: v[5],
                    parameter_estimation_variance_v_theta: v[6],
                };
                let rp = evaluate_diebold_mariano_west_hln_forecast_comparison_test(&mk(p));
                let rm = evaluate_diebold_mariano_west_hln_forecast_comparison_test(&mk(m));
                let fd = (rp.hln_adjusted_dmw_t_statistic - rm.hln_adjusted_dmw_t_statistic)
                    / (2.0 * eps);
                assert!(
                    r_ad.hln_adjusted_dmw_t_statistic.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.hln_adjusted_dmw_t_statistic.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use diebold_mariano_west_hln_mod::*;
mod andrews_ploberger_wald_break_mod {
    use ql_ad::Real;

    /// Input specification for the Andrews (1993, *Econometrica* 61(4)), Andrews & Ploberger
    /// (1994, *Econometrica* 62(6)), and Hansen (1997, *JBES* 15(1)) structural break test suite.
    #[derive(Debug, Clone, Copy)]
    pub struct AndrewsPlobergerStructuralBreakSpec<T: Real> {
        /// Total sample size $T_{\text{obs}}$ (e.g., $240.0$).
        pub sample_size_t_obs: T,
        /// Pre-break parameter estimate $\hat{\beta}_1$ at interior candidate break (e.g., $0.45$).
        pub pre_break_slope_beta1: T,
        /// Post-break parameter estimate $\hat{\beta}_2$ at interior candidate break (e.g., $0.92$).
        pub post_break_slope_beta2: T,
        /// Slope drift dispersion $\delta_{\text{profile}}$ across the trimmed candidate grid $[\pi_0, 1-\pi_0]$ (e.g., $0.14$).
        pub break_profile_dispersion_delta: T,
        /// Regression error variance $\hat{\sigma}_\varepsilon^2$ (e.g., $0.042$).
        pub residual_variance_sigma2: T,
        /// Regressor second moment $\hat{Q}_{XX} = \frac{1}{T}\sum X_t^2$ (e.g., $1.15$).
        pub regressor_second_moment_q_xx: T,
        /// Symmetric trimming fraction $\pi_0 \in (0.05, 0.30)$ (e.g., $0.15$).
        pub trimming_fraction_pi0: T,
    }

    /// Analytics returned by [`evaluate_andrews_ploberger_hansen_sup_exp_ave_wald_structural_break_test`].
    #[derive(Debug, Clone, Copy)]
    pub struct AndrewsPlobergerStructuralBreakAnalytics<T: Real> {
        /// Andrews & Ploberger (1994) Average Wald statistic $\text{AveW}_T$.
        pub average_wald_statistic_ave_w: T,
        /// Andrews & Ploberger (1994) optimal Exponential Wald statistic $\text{ExpW}_T$.
        pub exponential_wald_statistic_exp_w: T,
        /// Smooth Andrews (1993) Supremum Wald statistic $\text{SupW}_T$.
        pub supremum_wald_statistic_sup_w: T,
        /// Hansen (1997) asymptotic 5% critical-value excess $\text{ExpW}_T - 2.06$.
        pub hansen_exp_wald_5pct_excess: T,
    }

    /// Evaluates the Andrews (1993) and Andrews & Ploberger (1994) $\text{AveW}_T$, $\text{ExpW}_T$,
    /// and smooth $\text{SupW}_T$ structural change test statistics across a 4-point Gauss–Legendre
    /// quadrature grid on the trimmed break interval $[\pi_0, 1 - \pi_0]$.
    pub fn evaluate_andrews_ploberger_hansen_sup_exp_ave_wald_structural_break_test<T: Real>(
        spec: &AndrewsPlobergerStructuralBreakSpec<T>,
    ) -> AndrewsPlobergerStructuralBreakAnalytics<T> {
        let t_obs = spec.sample_size_t_obs;
        let diff_base = spec.post_break_slope_beta2 - spec.pre_break_slope_beta1;
        let scale = t_obs * spec.regressor_second_moment_q_xx / spec.residual_variance_sigma2;
        let pi0 = spec.trimming_fraction_pi0;
        let span = T::one() - T::from_f64(2.0) * pi0;

        let gl_nodes = [
            0.069_431_844_202_973_71,
            0.330_009_478_207_571_87,
            0.669_990_521_792_428_1,
            0.930_568_155_797_026_2,
        ];
        let gl_weights = [
            0.173_927_422_568_726_92,
            0.326_072_577_431_273_1,
            0.326_072_577_431_273_1,
            0.173_927_422_568_726_92,
        ];

        let half = T::from_f64(0.5);
        let eta_sup = T::from_f64(2.0);
        let mut ave_w = T::zero();
        let mut exp_half_sum = T::zero();
        let mut sup_exp_sum = T::zero();

        for idx in 0..4 {
            let u = T::from_f64(gl_nodes[idx]);
            let w = T::from_f64(gl_weights[idx]);
            let pi_k = pi0 + u * span;
            let centered = pi_k - half;
            let local_diff = diff_base
                + spec.break_profile_dispersion_delta * (T::from_f64(0.25) - centered * centered);
            let w_k = scale * pi_k * (T::one() - pi_k) * local_diff * local_diff;

            ave_w += w * w_k;
            exp_half_sum += w * (half * w_k).exp();
            sup_exp_sum += w * (eta_sup * w_k).exp();
        }

        let exp_w = exp_half_sum.ln();
        let sup_w = sup_exp_sum.ln() / eta_sup;
        // Andrews & Ploberger (1994, Table I, p=1, pi0=0.15) 5% critical value for ExpW is 2.06
        let excess_5pct = exp_w - T::from_f64(2.06);

        AndrewsPlobergerStructuralBreakAnalytics {
            average_wald_statistic_ave_w: ave_w,
            exponential_wald_statistic_exp_w: exp_w,
            supremum_wald_statistic_sup_w: sup_w,
            hansen_exp_wald_5pct_excess: excess_5pct,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn andrews_ploberger_wald_break_invariants_and_ad_parity() {
            let base = AndrewsPlobergerStructuralBreakSpec {
                sample_size_t_obs: 120.0_f64,
                pre_break_slope_beta1: 0.45_f64,
                post_break_slope_beta2: 0.72_f64,
                break_profile_dispersion_delta: 0.14_f64,
                residual_variance_sigma2: 0.55_f64,
                regressor_second_moment_q_xx: 1.15_f64,
                trimming_fraction_pi0: 0.15_f64,
            };
            let res =
                evaluate_andrews_ploberger_hansen_sup_exp_ave_wald_structural_break_test(&base);
            assert!(res.average_wald_statistic_ave_w > 0.0);
            assert!(res.exponential_wald_statistic_exp_w > 0.0);
            assert!(res.supremum_wald_statistic_sup_w > res.average_wald_statistic_ave_w);

            let inputs = [
                base.sample_size_t_obs,
                base.pre_break_slope_beta1,
                base.post_break_slope_beta2,
                base.break_profile_dispersion_delta,
                base.residual_variance_sigma2,
                base.regressor_second_moment_q_xx,
                base.trimming_fraction_pi0,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = AndrewsPlobergerStructuralBreakSpec {
                sample_size_t_obs: dual_in[0],
                pre_break_slope_beta1: dual_in[1],
                post_break_slope_beta2: dual_in[2],
                break_profile_dispersion_delta: dual_in[3],
                residual_variance_sigma2: dual_in[4],
                regressor_second_moment_q_xx: dual_in[5],
                trimming_fraction_pi0: dual_in[6],
            };
            let r_ad =
                evaluate_andrews_ploberger_hansen_sup_exp_ave_wald_structural_break_test(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| AndrewsPlobergerStructuralBreakSpec {
                    sample_size_t_obs: v[0],
                    pre_break_slope_beta1: v[1],
                    post_break_slope_beta2: v[2],
                    break_profile_dispersion_delta: v[3],
                    residual_variance_sigma2: v[4],
                    regressor_second_moment_q_xx: v[5],
                    trimming_fraction_pi0: v[6],
                };
                let rp = evaluate_andrews_ploberger_hansen_sup_exp_ave_wald_structural_break_test(
                    &mk(p),
                );
                let rm = evaluate_andrews_ploberger_hansen_sup_exp_ave_wald_structural_break_test(
                    &mk(m),
                );
                let fd = (rp.exponential_wald_statistic_exp_w
                    - rm.exponential_wald_statistic_exp_w)
                    / (2.0 * eps);
                assert!(
                    r_ad.exponential_wald_statistic_exp_w.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.exponential_wald_statistic_exp_w.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use andrews_ploberger_wald_break_mod::*;
mod ait_sahalia_jacod_jump_test_mod {
    use ql_ad::Real;

    /// Input specification for the Aït-Sahalia & Jacod (2009, *Annals of Statistics* 37(1);
    /// 2014, *High-Frequency Financial Econometrics*, Ch. 6 & Ch. 10) multi-scale power-variation
    /// jump test and Blumenthal–Getoor jump activity index estimator.
    #[derive(Debug, Clone, Copy)]
    pub struct AitSahaliaJacodJumpTestSpec<T: Real> {
        /// Fine-grid realized $p$-power variation $\hat{B}(p, \Delta_n)$ (e.g., $0.0085$).
        pub fine_grid_power_variation_b_fine: T,
        /// Coarse-grid ($k\Delta_n$) realized $p$-power variation $\hat{B}(p, k\Delta_n)$ (e.g., $0.0092$).
        pub coarse_grid_power_variation_b_coarse: T,
        /// Power exponent $p > 2$ (typically $p = 4.0$).
        pub power_exponent_p: T,
        /// Subsampling frequency multiplier $k \ge 2$ (typically $k = 2.0$).
        pub subsampling_multiplier_k: T,
        /// Integrated $2p$-power multipower variation $\hat{A}(2p, \Delta_n)$ for CLT standardization (e.g., $0.00018$).
        pub multipower_quarticity_scale_a2p: T,
        /// High-frequency observation count $n$ (e.g., $390.0$).
        pub high_frequency_sample_size_n: T,
    }

    /// Analytics returned by [`evaluate_ait_sahalia_jacod_power_variation_jump_activity_test`].
    #[derive(Debug, Clone, Copy)]
    pub struct AitSahaliaJacodJumpTestAnalytics<T: Real> {
        /// Observed coarse-to-fine power variation ratio $\hat{S}(p, k, \Delta_n) = \hat{B}(p, k\Delta_n) / \hat{B}(p, \Delta_n)$.
        pub power_variation_ratio_s_pk: T,
        /// Null no-jump diffusive limit $k^{p/2 - 1}$ (equals $2.0$ for $p=4, k=2$).
        pub diffusive_null_ratio_limit: T,
        /// Implied Blumenthal–Getoor jump activity index $\hat{\beta}_{\text{BG}} = 2 - \frac{2 \ln \hat{S}(p, k)}{\ln k}$ when jumps are present.
        pub implied_jump_activity_index_beta: T,
        /// Standardized Aït-Sahalia & Jacod (2009, Theorem 2) jump test $Z$-statistic.
        pub standardized_asj_jump_z_statistic: T,
    }

    /// Evaluates the Aït-Sahalia & Jacod (2009, Eqs. 8–15) coarse-to-fine power-variation ratio
    /// $\hat{S}(p, k, \Delta_n) = \frac{\hat{B}(p, k\Delta_n)}{\hat{B}(p, \Delta_n)}$, which converges
    /// to $1$ when jumps are present and to $k^{p/2 - 1}$ for continuous semimartingales, together
    /// with its standardized CLT test statistic.
    pub fn evaluate_ait_sahalia_jacod_power_variation_jump_activity_test<T: Real>(
        spec: &AitSahaliaJacodJumpTestSpec<T>,
    ) -> AitSahaliaJacodJumpTestAnalytics<T> {
        let b_fine = spec.fine_grid_power_variation_b_fine;
        let b_coarse = spec.coarse_grid_power_variation_b_coarse;
        let p = spec.power_exponent_p;
        let k = spec.subsampling_multiplier_k;
        let half = T::from_f64(0.5);

        let s_ratio = b_coarse / b_fine;
        let null_diffusive = k.pow(half * p - T::one());

        // Implied activity exponent from coarse/fine scaling
        let beta_bg = T::from_f64(2.0) * (T::one() - s_ratio.ln() / k.ln());

        // Asymptotic variance V_n(p, k) under the jump alternative (Ait-Sahalia & Jacod 2009, Theorem 2):
        // M(p, k) = (1 + k^(p-1)) / k^(p-1)
        let k_p_minus_1 = k.pow(p - T::one());
        let m_pk = (T::one() + k_p_minus_1) / k_p_minus_1;
        let asymp_var = m_pk * spec.multipower_quarticity_scale_a2p
            / (b_fine * b_fine * spec.high_frequency_sample_size_n.sqrt());
        let z_stat = (null_diffusive - s_ratio) / asymp_var.sqrt();

        AitSahaliaJacodJumpTestAnalytics {
            power_variation_ratio_s_pk: s_ratio,
            diffusive_null_ratio_limit: null_diffusive,
            implied_jump_activity_index_beta: beta_bg,
            standardized_asj_jump_z_statistic: z_stat,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn ait_sahalia_jacod_jump_test_invariants_and_ad_parity() {
            let base = AitSahaliaJacodJumpTestSpec {
                fine_grid_power_variation_b_fine: 0.0085_f64,
                coarse_grid_power_variation_b_coarse: 0.0092_f64,
                power_exponent_p: 4.0_f64,
                subsampling_multiplier_k: 2.0_f64,
                multipower_quarticity_scale_a2p: 0.00018_f64,
                high_frequency_sample_size_n: 390.0_f64,
            };
            let res = evaluate_ait_sahalia_jacod_power_variation_jump_test(&base);
            assert!(res.power_variation_ratio_s_pk > 1.0 && res.power_variation_ratio_s_pk < 1.5);
            assert!((res.diffusive_null_ratio_limit - 2.0).abs() < 1e-12);
            assert!(
                res.implied_jump_activity_index_beta > 1.0
                    && res.implied_jump_activity_index_beta < 2.0
            );
            assert!(res.standardized_asj_jump_z_statistic > 1.96);

            let inputs = [
                base.fine_grid_power_variation_b_fine,
                base.coarse_grid_power_variation_b_coarse,
                base.power_exponent_p,
                base.subsampling_multiplier_k,
                base.multipower_quarticity_scale_a2p,
                base.high_frequency_sample_size_n,
            ];
            let mut dual_in = [Dual::<6>::from_value(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = AitSahaliaJacodJumpTestSpec {
                fine_grid_power_variation_b_fine: dual_in[0],
                coarse_grid_power_variation_b_coarse: dual_in[1],
                power_exponent_p: dual_in[2],
                subsampling_multiplier_k: dual_in[3],
                multipower_quarticity_scale_a2p: dual_in[4],
                high_frequency_sample_size_n: dual_in[5],
            };
            let r_ad = evaluate_ait_sahalia_jacod_power_variation_jump_activity_test(&spec_ad);
            let eps = 1e-8;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| AitSahaliaJacodJumpTestSpec {
                    fine_grid_power_variation_b_fine: v[0],
                    coarse_grid_power_variation_b_coarse: v[1],
                    power_exponent_p: v[2],
                    subsampling_multiplier_k: v[3],
                    multipower_quarticity_scale_a2p: v[4],
                    high_frequency_sample_size_n: v[5],
                };
                let rp = evaluate_ait_sahalia_jacod_power_variation_jump_activity_test(&mk(p));
                let rm = evaluate_ait_sahalia_jacod_power_variation_jump_activity_test(&mk(m));
                let fd = (rp.standardized_asj_jump_z_statistic
                    - rm.standardized_asj_jump_z_statistic)
                    / (2.0 * eps);
                assert!(
                    r_ad.standardized_asj_jump_z_statistic.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.standardized_asj_jump_z_statistic.d[i] - fd).abs() < 1e-4);
            }
        }

        fn evaluate_ait_sahalia_jacod_power_variation_jump_test<T: Real>(
            spec: &AitSahaliaJacodJumpTestSpec<T>,
        ) -> AitSahaliaJacodJumpTestAnalytics<T> {
            evaluate_ait_sahalia_jacod_power_variation_jump_activity_test(spec)
        }
    }
}
pub use ait_sahalia_jacod_jump_test_mod::*;
mod pesaran_shin_smith_ardl_bounds_mod {
    use ql_ad::Real;

    /// Input specification for the Pesaran, Shin & Smith (2001, *Journal of Applied Econometrics*
    /// 16(3)) Autoregressive Distributed Lag (ARDL) bounds test for level relationships
    /// across mixed $I(0)/I(1)$ regressors.
    #[derive(Debug, Clone, Copy)]
    pub struct PesaranShinSmithArdlBoundsSpec<T: Real> {
        /// Error-correction speed of adjustment $\pi_{yy} < 0$ on lagged dependent level $y_{t-1}$ (e.g., $-0.24$).
        pub lagged_dependent_coef_pi_yy: T,
        /// Long-run forcing coefficient $\pi_{yx}$ on lagged regressor level $x_{t-1}$ (e.g., $0.36$).
        pub lagged_regressor_coef_pi_yx: T,
        /// Sample variance $s_{yy}$ of demeaned $y_{t-1}$ (e.g., $1.45$).
        pub lagged_y_variance_s_yy: T,
        /// Sample variance $s_{xx}$ of demeaned $x_{t-1}$ (e.g., $1.80$).
        pub lagged_x_variance_s_xx: T,
        /// Sample covariance $s_{yx}$ between demeaned $y_{t-1}$ and $x_{t-1}$ (e.g., $0.95$).
        pub lagged_yx_covariance_s_yx: T,
        /// Conditional ECM residual variance $\hat{\sigma}_\varepsilon^2$ (e.g., $0.085$).
        pub residual_variance_sigma2: T,
        /// Effective sample observations $T_{\text{obs}}$ (e.g., $120.0$).
        pub sample_size_t_obs: T,
    }

    /// Analytics returned by [`evaluate_pesaran_shin_smith_ardl_bounds_cointegration_test`].
    #[derive(Debug, Clone, Copy)]
    pub struct PesaranShinSmithArdlBoundsAnalytics<T: Real> {
        /// Implied long-run level equilibrium multiplier $\hat{\theta}_{\text{LR}} = -\pi_{yx} / \pi_{yy}$.
        pub long_run_level_multiplier_theta: T,
        /// Banerjee–Dolado–Mestre $t$-bounds statistic $t_{\text{BDM}}$ on $\pi_{yy} = 0$.
        pub bdm_t_bounds_statistic: T,
        /// Pesaran–Shin–Smith joint $F_{\text{PSS}}$-bounds statistic testing $H_0: \pi_{yy} = \pi_{yx} = 0$.
        pub pss_f_bounds_statistic: T,
        /// Long-run multiplier delta-method standard error $\text{SE}(\hat{\theta}_{\text{LR}})$.
        pub long_run_multiplier_standard_error: T,
    }

    /// Evaluates the Pesaran, Shin & Smith (2001, Eqs. 16–22) ARDL conditional error-correction
    /// bounds test statistics $F_{\text{PSS}}$ and $t_{\text{BDM}}$, together with the long-run
    /// equilibrium level multiplier $\hat{\theta}_{\text{LR}} = -\pi_{yx}/\pi_{yy}$ and its delta-method SE.
    pub fn evaluate_pesaran_shin_smith_ardl_bounds_cointegration_test<T: Real>(
        spec: &PesaranShinSmithArdlBoundsSpec<T>,
    ) -> PesaranShinSmithArdlBoundsAnalytics<T> {
        let pi_yy = spec.lagged_dependent_coef_pi_yy;
        let pi_yx = spec.lagged_regressor_coef_pi_yx;
        let s_yy = spec.lagged_y_variance_s_yy;
        let s_xx = spec.lagged_x_variance_s_xx;
        let s_yx = spec.lagged_yx_covariance_s_yx;
        let sig2 = spec.residual_variance_sigma2;
        let t_obs = spec.sample_size_t_obs;

        let theta_lr = -pi_yx / pi_yy;

        // Joint quadratic form pi^T (T * S_zz) pi / (2 * sig2) for H_0: pi_yy = pi_yx = 0
        let quad =
            pi_yy * pi_yy * s_yy + T::from_f64(2.0) * pi_yy * pi_yx * s_yx + pi_yx * pi_yx * s_xx;
        let f_pss = (t_obs * quad) / (T::from_f64(2.0) * sig2);

        // Inverse covariance of (pi_yy, pi_yx): Var(pi) = (sig2 / (T * det_S)) * [[s_xx, -s_yx], [-s_yx, s_yy]]
        let det_s = s_yy * s_xx - s_yx * s_yx;
        let var_pi_yy = (sig2 * s_xx) / (t_obs * det_s);
        let var_pi_yx = (sig2 * s_yy) / (t_obs * det_s);
        let cov_pi = -(sig2 * s_yx) / (t_obs * det_s);

        let t_bdm = pi_yy / var_pi_yy.sqrt();

        // Delta-method variance of theta_LR = -pi_yx / pi_yy:
        // d(theta)/d(pi_yy) = pi_yx / pi_yy^2 = -theta / pi_yy, d(theta)/d(pi_yx) = -1 / pi_yy
        let inv_pi_yy2 = T::one() / (pi_yy * pi_yy);
        let var_theta = inv_pi_yy2
            * (theta_lr * theta_lr * var_pi_yy + T::from_f64(2.0) * theta_lr * cov_pi + var_pi_yx);
        let se_theta = var_theta.sqrt();

        PesaranShinSmithArdlBoundsAnalytics {
            long_run_level_multiplier_theta: theta_lr,
            bdm_t_bounds_statistic: t_bdm,
            pss_f_bounds_statistic: f_pss,
            long_run_multiplier_standard_error: se_theta,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn pesaran_shin_smith_ardl_bounds_invariants_and_ad_parity() {
            let base = PesaranShinSmithArdlBoundsSpec {
                lagged_dependent_coef_pi_yy: -0.24_f64,
                lagged_regressor_coef_pi_yx: 0.36_f64,
                lagged_y_variance_s_yy: 1.45_f64,
                lagged_x_variance_s_xx: 1.80_f64,
                lagged_yx_covariance_s_yx: 0.95_f64,
                residual_variance_sigma2: 0.085_f64,
                sample_size_t_obs: 120.0_f64,
            };
            let res = evaluate_pesaran_shin_smith_ardl_bounds_cointegration_test(&base);
            assert!((res.long_run_level_multiplier_theta - 1.5).abs() < 1e-12);
            assert!(res.bdm_t_bounds_statistic < -3.5);
            // Exceeds Pesaran-Shin-Smith (2001) upper I(1) critical bound (~5.73 at 1%)
            assert!(res.pss_f_bounds_statistic > 6.0);
            assert!(res.long_run_multiplier_standard_error > 0.0);

            let inputs = [
                base.lagged_dependent_coef_pi_yy,
                base.lagged_regressor_coef_pi_yx,
                base.lagged_y_variance_s_yy,
                base.lagged_x_variance_s_xx,
                base.lagged_yx_covariance_s_yx,
                base.residual_variance_sigma2,
                base.sample_size_t_obs,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = PesaranShinSmithArdlBoundsSpec {
                lagged_dependent_coef_pi_yy: dual_in[0],
                lagged_regressor_coef_pi_yx: dual_in[1],
                lagged_y_variance_s_yy: dual_in[2],
                lagged_x_variance_s_xx: dual_in[3],
                lagged_yx_covariance_s_yx: dual_in[4],
                residual_variance_sigma2: dual_in[5],
                sample_size_t_obs: dual_in[6],
            };
            let r_ad = evaluate_pesaran_shin_smith_ardl_bounds_cointegration_test(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| PesaranShinSmithArdlBoundsSpec {
                    lagged_dependent_coef_pi_yy: v[0],
                    lagged_regressor_coef_pi_yx: v[1],
                    lagged_y_variance_s_yy: v[2],
                    lagged_x_variance_s_xx: v[3],
                    lagged_yx_covariance_s_yx: v[4],
                    residual_variance_sigma2: v[5],
                    sample_size_t_obs: v[6],
                };
                let rp = evaluate_pesaran_shin_smith_ardl_bounds_cointegration_test(&mk(p));
                let rm = evaluate_pesaran_shin_smith_ardl_bounds_cointegration_test(&mk(m));
                let fd = (rp.long_run_multiplier_standard_error
                    - rm.long_run_multiplier_standard_error)
                    / (2.0 * eps);
                assert!(
                    r_ad.long_run_multiplier_standard_error.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.long_run_multiplier_standard_error.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use pesaran_shin_smith_ardl_bounds_mod::*;
mod hansen_white_spa_test_mod {
    use ql_ad::Real;

    /// Input specification for the Hansen (2005, *Journal of Business & Economic Statistics*
    /// 23(4)) Test for Superior Predictive Ability ($\text{SPA}$) & White (2000, *Econometrica*
    /// 68(5)) Reality Check across competing forecasting models.
    #[derive(Debug, Clone, Copy)]
    pub struct HansenWhiteSuperiorPredictiveAbilitySpec<T: Real> {
        /// Sample mean loss differential $\bar{d}_1 > 0$ of candidate model 1 vs. benchmark (e.g., $0.018$).
        pub model_1_mean_loss_reduction_d1: T,
        /// Sample mean loss differential $\bar{d}_2 < 0$ of clearly inferior model 2 vs. benchmark (e.g., $-0.024$).
        pub model_2_mean_loss_reduction_d2: T,
        /// HAC long-run variance $\hat{\omega}_1^2$ of model 1 loss differentials (e.g., $0.065$).
        pub model_1_hac_variance_w1: T,
        /// HAC long-run variance $\hat{\omega}_2^2$ of model 2 loss differentials (e.g., $0.080$).
        pub model_2_hac_variance_w2: T,
        /// Cross-model loss differential correlation $\rho_{12} \in (-1, 1)$ (e.g., $0.35$).
        pub cross_model_correlation_rho12: T,
        /// Out-of-sample evaluation length $n$ (e.g., $150.0$).
        pub out_of_sample_observations_n: T,
    }

    /// Analytics returned by [`evaluate_hansen_white_superior_predictive_ability_test`].
    #[derive(Debug, Clone, Copy)]
    pub struct HansenWhiteSuperiorPredictiveAbilityAnalytics<T: Real> {
        /// Studentized $t$-statistic $t_1 = \sqrt{n}\bar{d}_1 / \hat{\omega}_1$ of candidate model 1.
        pub model_1_studentized_t_statistic: T,
        /// Hansen (2005) recentered mean $\hat{\mu}_2^c \le 0$ for the poor alternative model 2.
        pub model_2_recentered_null_mean: T,
        /// Smooth Hansen $\text{SPA}$ studentized max statistic $T_n^{\text{SPA}}$.
        pub hansen_spa_test_statistic: T,
        /// Multiplicity-adjusted Hansen $\text{SPA}$ asymptotic $p$-value.
        pub hansen_spa_adjusted_p_value: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        T::from_f64(0.5) * (T::one() + (x * T::from_f64(0.707_106_781_186_547_5)).erf())
    }

    /// Evaluates the Hansen (2005, Eqs. 5–9) studentized $\text{SPA}$ test statistic and
    /// recentered null distribution, which eliminates the conservative bias of White's (2000)
    /// Reality Check when poor/irrelevant models ($\sqrt{n}\bar{d}_k / \hat{\omega}_k < -\sqrt{2\ln\ln n}$)
    /// are included in the model universe.
    pub fn evaluate_hansen_white_superior_predictive_ability_test<T: Real>(
        spec: &HansenWhiteSuperiorPredictiveAbilitySpec<T>,
    ) -> HansenWhiteSuperiorPredictiveAbilityAnalytics<T> {
        let n = spec.out_of_sample_observations_n;
        let sqrt_n = n.sqrt();
        let w1 = spec.model_1_hac_variance_w1.sqrt();
        let w2 = spec.model_2_hac_variance_w2.sqrt();

        let t1 = sqrt_n * spec.model_1_mean_loss_reduction_d1 / w1;
        let t2 = sqrt_n * spec.model_2_mean_loss_reduction_d2 / w2;

        // Hansen (2005, Eq. 7) threshold -sqrt(2 * ln(ln(n))) for recentering negative mean alternatives
        let thresh = -(T::from_f64(2.0) * n.ln().ln()).sqrt();
        // Smooth indicator for t2 <= thresh
        let weight_recenter = std_norm_cdf((thresh - t2) * T::from_f64(2.0));
        let mu2_c = spec.model_2_mean_loss_reduction_d2 * weight_recenter;

        // Smooth studentized SPA statistic max(t1, t2, 0) via log-sum-exp with temperature beta = 8.0
        let beta = T::from_f64(8.0);
        let t_spa = t1 + ((T::one() + (beta * (t2 - t1)).exp() + (-beta * t1).exp()).ln()) / beta;

        // Effective number of independent tests under correlation rho12 and Hansen recentering
        let rho = spec.cross_model_correlation_rho12;
        let m_eff =
            T::one() + (T::one() - weight_recenter * T::from_f64(0.5)) * (T::one() - rho * rho);
        let p_single = T::one() - std_norm_cdf(t_spa);
        let p_spa = T::one() - (T::one() - p_single).pow(m_eff);

        HansenWhiteSuperiorPredictiveAbilityAnalytics {
            model_1_studentized_t_statistic: t1,
            model_2_recentered_null_mean: mu2_c,
            hansen_spa_test_statistic: t_spa,
            hansen_spa_adjusted_p_value: p_spa,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hansen_white_spa_test_invariants_and_ad_parity() {
            let base = HansenWhiteSuperiorPredictiveAbilitySpec {
                model_1_mean_loss_reduction_d1: 0.045_f64,
                model_2_mean_loss_reduction_d2: -0.024_f64,
                model_1_hac_variance_w1: 0.065_f64,
                model_2_hac_variance_w2: 0.080_f64,
                cross_model_correlation_rho12: 0.35_f64,
                out_of_sample_observations_n: 150.0_f64,
            };
            let res = evaluate_hansen_white_superior_predictive_ability_test(&base);
            assert!(res.model_1_studentized_t_statistic > 1.96);
            assert!(res.model_2_recentered_null_mean < 0.0);
            assert!(res.hansen_spa_test_statistic >= res.model_1_studentized_t_statistic);
            assert!(
                res.hansen_spa_adjusted_p_value > 0.0 && res.hansen_spa_adjusted_p_value < 0.05
            );

            let inputs = [
                base.model_1_mean_loss_reduction_d1,
                base.model_2_mean_loss_reduction_d2,
                base.model_1_hac_variance_w1,
                base.model_2_hac_variance_w2,
                base.cross_model_correlation_rho12,
                base.out_of_sample_observations_n,
            ];
            let mut dual_in = [Dual::<6>::from_value(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HansenWhiteSuperiorPredictiveAbilitySpec {
                model_1_mean_loss_reduction_d1: dual_in[0],
                model_2_mean_loss_reduction_d2: dual_in[1],
                model_1_hac_variance_w1: dual_in[2],
                model_2_hac_variance_w2: dual_in[3],
                cross_model_correlation_rho12: dual_in[4],
                out_of_sample_observations_n: dual_in[5],
            };
            let r_ad = evaluate_hansen_white_superior_predictive_ability_test(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| HansenWhiteSuperiorPredictiveAbilitySpec {
                    model_1_mean_loss_reduction_d1: v[0],
                    model_2_mean_loss_reduction_d2: v[1],
                    model_1_hac_variance_w1: v[2],
                    model_2_hac_variance_w2: v[3],
                    cross_model_correlation_rho12: v[4],
                    out_of_sample_observations_n: v[5],
                };
                let rp = evaluate_hansen_white_superior_predictive_ability_test(&mk(p));
                let rm = evaluate_hansen_white_superior_predictive_ability_test(&mk(m));
                let fd =
                    (rp.hansen_spa_adjusted_p_value - rm.hansen_spa_adjusted_p_value) / (2.0 * eps);
                assert!(
                    r_ad.hansen_spa_adjusted_p_value.d[i].abs() > 1e-7,
                    "zero derivative at {i}"
                );
                assert!((r_ad.hansen_spa_adjusted_p_value.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use hansen_white_spa_test_mod::*;
