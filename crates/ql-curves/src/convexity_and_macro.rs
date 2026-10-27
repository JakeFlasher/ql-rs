mod evaluate_romer_veronesi_taylor_rule_macro_yield_curve_module {
    use ql_ad::Real;

    #[derive(Clone, Debug)]
    pub struct RomerVeronesiTaylorRuleMacroCurveSpec<T: Real = f64> {
        pub natural_real_rate_r_star: T,
        pub inflation_target_pi_star: T,
        pub current_inflation: T,
        pub current_output_gap: T,
        pub taylor_inflation_weight_phi_pi: T,
        pub taylor_output_weight_phi_y: T,
        pub inflation_mean_reversion: T,
        pub output_gap_mean_reversion: T,
        pub inflation_vol: T,
        pub output_gap_vol: T,
        pub market_price_of_inflation_risk: T,
        pub zlb_floor_rate: T,
        pub maturities: Vec<T>,
    }
    #[derive(Clone, Debug)]
    pub struct RomerVeronesiTaylorRuleMacroCurveReport<T: Real = f64> {
        pub unconstrained_taylor_shadow_rate: T,
        pub zlb_constrained_policy_rate: T,
        pub macro_affine_zero_yields: Vec<T>,
        pub expectations_hypothesis_yields: Vec<T>,
        pub macro_term_premia_bps: Vec<T>,
        pub yield_inflation_beta_10y: T,
        pub yield_output_gap_beta_10y: T,
    }
    pub fn evaluate_romer_veronesi_taylor_rule_macro_yield_curve<T: Real>(
        s: &RomerVeronesiTaylorRuleMacroCurveSpec<T>,
    ) -> RomerVeronesiTaylorRuleMacroCurveReport<T> {
        let shadow = s.natural_real_rate_r_star
            + s.inflation_target_pi_star
            + s.taylor_inflation_weight_phi_pi * (s.current_inflation - s.inflation_target_pi_star)
            + s.taylor_output_weight_phi_y * s.current_output_gap;
        let policy = shadow.max_val(s.zlb_floor_rate);
        let mut affine = Vec::with_capacity(s.maturities.len());
        let mut eh = Vec::with_capacity(s.maturities.len());
        let mut premia = Vec::with_capacity(s.maturities.len());
        let b_load = |kappa: T, tau: T| -> T {
            let kt = (kappa * tau).max_val(T::from_f64(1e-8));
            (T::one() - (-kt).exp()) / kt
        };
        for &t in &s.maturities {
            let bi_t = b_load(s.inflation_mean_reversion, t);
            let by_t = b_load(s.output_gap_mean_reversion, t);
            let expected = s.natural_real_rate_r_star
                + s.inflation_target_pi_star
                + s.taylor_inflation_weight_phi_pi
                    * (s.current_inflation - s.inflation_target_pi_star)
                    * bi_t
                + s.taylor_output_weight_phi_y * s.current_output_gap * by_t;
            let risk = s.market_price_of_inflation_risk
                * s.inflation_vol
                * s.inflation_vol
                * (T::one() - bi_t)
                + s.output_gap_vol * s.output_gap_vol * (T::one() - by_t) * T::from_f64(0.5);
            let y = expected + risk;
            affine.push(y);
            eh.push(expected);
            premia.push(risk * T::from_f64(10_000.0));
        }
        let horizon = T::from_f64(10.0);
        let bi = s.taylor_inflation_weight_phi_pi * b_load(s.inflation_mean_reversion, horizon);
        let by = s.taylor_output_weight_phi_y * b_load(s.output_gap_mean_reversion, horizon);
        RomerVeronesiTaylorRuleMacroCurveReport {
            unconstrained_taylor_shadow_rate: shadow,
            zlb_constrained_policy_rate: policy,
            macro_affine_zero_yields: affine,
            expectations_hypothesis_yields: eh,
            macro_term_premia_bps: premia,
            yield_inflation_beta_10y: bi,
            yield_output_gap_beta_10y: by,
        }
    }

    #[cfg(test)]
    mod domain_suite_42_tests {
        use super::*;
        use ql_ad::Dual;
        fn make<T: Real>(x: T) -> RomerVeronesiTaylorRuleMacroCurveSpec<T> {
            RomerVeronesiTaylorRuleMacroCurveSpec {
                natural_real_rate_r_star: T::from_f64(0.02) + x,
                inflation_target_pi_star: T::from_f64(0.02),
                current_inflation: T::from_f64(0.025),
                current_output_gap: T::from_f64(0.01),
                taylor_inflation_weight_phi_pi: T::from_f64(1.5),
                taylor_output_weight_phi_y: T::from_f64(0.5),
                inflation_mean_reversion: T::from_f64(0.4),
                output_gap_mean_reversion: T::from_f64(0.3),
                inflation_vol: T::from_f64(0.2),
                output_gap_vol: T::from_f64(0.1),
                market_price_of_inflation_risk: T::from_f64(0.2),
                zlb_floor_rate: T::zero(),
                maturities: vec![T::one(), T::from_f64(10.0)],
            }
        }
        #[test]
        fn suite42_taylor_invariants_and_dual_fd() {
            let r = evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&make(0.0));
            assert!(r.zlb_constrained_policy_rate >= 0.0);
            let mut s0 = make(0.0);
            s0.current_inflation = s0.inflation_target_pi_star;
            s0.current_output_gap = 0.0;
            let r0 = evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&s0);
            assert!(
                (r0.unconstrained_taylor_shadow_rate
                    - (s0.natural_real_rate_r_star + s0.inflation_target_pi_star))
                    .abs()
                    < 1e-12
            );
            let h = 1e-6;
            let fd = (evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&make(h))
                .unconstrained_taylor_shadow_rate
                - evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&make(-h))
                    .unconstrained_taylor_shadow_rate)
                / (2.0 * h);
            let ad = evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&make(
                Dual::<1>::variable(0.0, 0),
            ))
            .unconstrained_taylor_shadow_rate;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
            let mut sp = make(0.0);
            let mut sm = make(0.0);
            sp.current_inflation += h;
            sm.current_inflation -= h;
            let fd_y = (evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&sp)
                .macro_affine_zero_yields[0]
                - evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&sm)
                    .macro_affine_zero_yields[0])
                / (2.0 * h);
            let mut sad = make(Dual::<1>::from_f64(0.0));
            sad.current_inflation = Dual::<1>::variable(0.025, 0);
            let ad_y = evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&sad)
                .macro_affine_zero_yields[0];
            assert!(ad_y.d[0].abs() > 1e-6 && (ad_y.d[0] - fd_y).abs() < 1e-5);
        }
    }
}
pub use evaluate_romer_veronesi_taylor_rule_macro_yield_curve_module::*;

mod evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization_module {
    use ql_ad::Real;

    #[derive(Clone, Debug)]
    pub struct FabozziKeyRateLscImmunizationSpec<T: Real = f64> {
        pub cashflow_tenors: Vec<T>,
        pub cashflows: Vec<T>,
        pub zero_rates: Vec<T>,
        pub key_rate_tenors: Vec<T>,
        pub ns_decay_lambda: T,
        pub liability_horizon_years: T,
        pub barbell_short_tenor: T,
        pub barbell_long_tenor: T,
        pub twist_slope_shock_bps: T,
        pub butterfly_curvature_shock_bps: T,
    }

    #[derive(Clone, Debug)]
    pub struct FabozziKeyRateLscImmunizationReport<T: Real = f64> {
        pub bond_present_value: T,
        pub fisher_weil_duration: T,
        pub key_rate_durations: Vec<T>,
        pub krd_partition_unity_error: T,
        pub lsc_level_duration: T,
        pub lsc_slope_duration: T,
        pub lsc_curvature_duration: T,
        pub barbell_short_weight: T,
        pub barbell_convexity_advantage: T,
        pub non_parallel_shock_pnl_pct: T,
    }

    pub fn evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization<T: Real>(
        s: &FabozziKeyRateLscImmunizationSpec<T>,
    ) -> FabozziKeyRateLscImmunizationReport<T> {
        let n = s
            .cashflow_tenors
            .len()
            .min(s.cashflows.len())
            .min(s.zero_rates.len());
        assert!(n > 0 && !s.key_rate_tenors.is_empty());
        let mut pv = T::zero();
        let mut dur_num = T::zero();
        let mut conv_num = T::zero();
        let mut lsc_slope = T::zero();
        let mut lsc_curv = T::zero();
        let lambda = s.ns_decay_lambda.max_val(T::from_f64(1e-12));
        let mut krd_num = vec![T::zero(); s.key_rate_tenors.len()];
        for i in 0..n {
            let tau = s.cashflow_tenors[i];
            let pv_cf = s.cashflows[i] * (-s.zero_rates[i] * tau).exp();
            pv += pv_cf;
            dur_num += tau * pv_cf;
            conv_num += tau * tau * pv_cf;
            let x = lambda * tau;
            let e = (-x).exp();
            let l2 = (T::one() - e) / x.max_val(T::from_f64(1e-12));
            lsc_slope += tau * l2 * pv_cf;
            lsc_curv += tau * (l2 - e) * pv_cf;
            for k in 0..s.key_rate_tenors.len() {
                let tk = s.key_rate_tenors[k];
                let w = if (k == 0 && tau <= tk) || (k + 1 == s.key_rate_tenors.len() && tau >= tk)
                {
                    T::one()
                } else if k > 0 && tau >= s.key_rate_tenors[k - 1] && tau <= tk {
                    (tau - s.key_rate_tenors[k - 1])
                        / (tk - s.key_rate_tenors[k - 1]).max_val(T::from_f64(1e-12))
                } else if k + 1 < s.key_rate_tenors.len()
                    && tau >= tk
                    && tau <= s.key_rate_tenors[k + 1]
                {
                    (s.key_rate_tenors[k + 1] - tau)
                        / (s.key_rate_tenors[k + 1] - tk).max_val(T::from_f64(1e-12))
                } else {
                    T::zero()
                };
                krd_num[k] += tau * w * pv_cf;
            }
        }
        let den = pv.abs().max_val(T::from_f64(1e-12));
        let duration = dur_num / den;
        let krd: Vec<T> = krd_num.into_iter().map(|v| v / den).collect();
        let sum_krd = krd.iter().fold(T::zero(), |a, &b| a + b);
        let short = s.barbell_short_tenor;
        let long = s.barbell_long_tenor;
        let horizon = s.liability_horizon_years;
        let ws = (long - horizon) / (long - short).max_val(T::from_f64(1e-12));
        let barbell_conv = ws * short * short + (T::one() - ws) * long * long;
        let bullet_conv = horizon * horizon;
        let non_parallel = -(lsc_slope / den * s.twist_slope_shock_bps * T::from_f64(1e-4)
            + lsc_curv / den * s.butterfly_curvature_shock_bps * T::from_f64(1e-4));
        FabozziKeyRateLscImmunizationReport {
            bond_present_value: pv,
            fisher_weil_duration: duration,
            key_rate_durations: krd,
            krd_partition_unity_error: (duration - sum_krd).abs(),
            lsc_level_duration: duration,
            lsc_slope_duration: lsc_slope / den,
            lsc_curvature_duration: lsc_curv / den,
            barbell_short_weight: ws,
            barbell_convexity_advantage: barbell_conv - bullet_conv,
            non_parallel_shock_pnl_pct: non_parallel,
        }
    }

    #[cfg(test)]
    mod domain_suite_43_tests {
        use super::*;
        use ql_ad::Dual;

        fn spec<T: Real>(x: T) -> FabozziKeyRateLscImmunizationSpec<T> {
            FabozziKeyRateLscImmunizationSpec {
                cashflow_tenors: vec![T::one(), T::from_f64(5.0), T::from_f64(10.0)],
                cashflows: vec![T::from_f64(2.0) + x, T::from_f64(2.0), T::from_f64(102.0)],
                zero_rates: vec![T::from_f64(0.03), T::from_f64(0.035), T::from_f64(0.04)],
                key_rate_tenors: vec![T::one(), T::from_f64(5.0), T::from_f64(10.0)],
                ns_decay_lambda: T::from_f64(0.5),
                liability_horizon_years: T::from_f64(5.0),
                barbell_short_tenor: T::one(),
                barbell_long_tenor: T::from_f64(10.0),
                twist_slope_shock_bps: T::from_f64(20.0),
                butterfly_curvature_shock_bps: T::from_f64(10.0),
            }
        }

        #[test]
        fn invariants_and_dual_fd() {
            let r = evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization(&spec(0.0));
            assert!(r.krd_partition_unity_error < 1e-12);
            assert!(r.barbell_convexity_advantage > 0.0);
            let h = 1e-6;
            let fd = (evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization(&spec(h))
                .bond_present_value
                - evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization(&spec(-h))
                    .bond_present_value)
                / (2.0 * h);
            let ad = evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization(&spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .bond_present_value;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization_module::*;

mod evaluate_brigo_mercurio_rebonato_swaption_and_cms_module {
    use ql_ad::Real;
    #[derive(Clone, Debug)]
    pub struct BrigoMercurioRebonatoCmsSpec<T: Real = f64> {
        pub fixing_tenors: Vec<T>,
        pub accrual_taus: Vec<T>,
        pub forward_rates: Vec<T>,
        pub discount_factors: Vec<T>,
        pub forward_vols: Vec<T>,
        pub corr_decay_beta: T,
        pub option_expiry: T,
        pub swaption_strike: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct BrigoMercurioRebonatoCmsReport<T: Real = f64> {
        pub swap_annuity_pv01: T,
        pub forward_swap_rate: T,
        pub rebonato_lmm_swaption_vol: T,
        pub hagan_annuity_convexity_ratio: T,
        pub cms_convexity_adjustment_bps: T,
        pub convexity_adjusted_cms_rate: T,
        pub black76_payer_swaption_pv: T,
        pub weights_sum_conservation_error: T,
    }
    pub fn evaluate_brigo_mercurio_rebonato_swaption_and_cms<T: Real>(
        s: &BrigoMercurioRebonatoCmsSpec<T>,
    ) -> BrigoMercurioRebonatoCmsReport<T> {
        let n = s.accrual_taus.len();
        let mut a = T::zero();
        for i in 0..n {
            a += s.accrual_taus[i] * s.discount_factors[i];
        }
        let mut sw = T::zero();
        for i in 0..n {
            sw += s.accrual_taus[i] * s.discount_factors[i] * s.forward_rates[i];
        }
        let sw = sw / a;
        let mut var = T::zero();
        let mut ws = T::zero();
        for i in 0..n {
            let wi = s.accrual_taus[i] * s.discount_factors[i] / a;
            ws += wi;
            for j in 0..n {
                let wj = s.accrual_taus[j] * s.discount_factors[j] / a;
                let d = if i == j {
                    T::zero()
                } else {
                    let z = s.fixing_tenors[i] - s.fixing_tenors[j];
                    (z * z + T::from_f64(1e-30)).sqrt()
                };
                let rho = (-(s.corr_decay_beta * d)).exp();
                var += wi
                    * wj
                    * s.forward_rates[i]
                    * s.forward_rates[j]
                    * rho
                    * s.forward_vols[i]
                    * s.forward_vols[j]
                    / (sw * sw);
            }
        }
        let vol = var.sqrt();
        let x = sw.max_val(T::from_f64(1e-8));
        let one = T::one();
        let nf = T::from_f64(n as f64);
        let u = one + x;
        let g = (one - u.powf(-(n as f64))) / x;
        let gp = (nf * u.powf(-(n as f64) - 1.0)) / x - g / x;
        let gpp = -(nf * (nf + one) * u.powf(-(n as f64) - 2.0)) / x - T::from_f64(2.0) * gp / x;
        let ratio = -(gpp / gp);
        let adj = T::from_f64(0.5) * x * x * ratio * ((vol * vol * s.option_expiry).exp() - one);
        let cms = x + adj;
        let d1 = ((sw / s.swaption_strike).ln() + T::from_f64(0.5) * vol * vol * s.option_expiry)
            / (vol * s.option_expiry.sqrt());
        let d2 = d1 - vol * s.option_expiry.sqrt();
        let pv = a * (sw * d1.norm_cdf() - s.swaption_strike * d2.norm_cdf());
        BrigoMercurioRebonatoCmsReport {
            swap_annuity_pv01: a,
            forward_swap_rate: sw,
            rebonato_lmm_swaption_vol: vol,
            hagan_annuity_convexity_ratio: ratio,
            cms_convexity_adjustment_bps: adj * T::from_f64(10000.0),
            convexity_adjusted_cms_rate: cms,
            black76_payer_swaption_pv: pv,
            weights_sum_conservation_error: (ws - one).abs(),
        }
    }
    #[cfg(test)]
    mod domain_suite_45_tests {
        use super::*;
        use ql_ad::Dual;
        fn spec<T: Real>(x: T) -> BrigoMercurioRebonatoCmsSpec<T> {
            BrigoMercurioRebonatoCmsSpec {
                fixing_tenors: vec![T::one(), T::from_f64(2.)],
                accrual_taus: vec![T::one(), T::one()],
                forward_rates: vec![T::from_f64(0.03) + x, T::from_f64(0.035)],
                discount_factors: vec![T::from_f64(0.98), T::from_f64(0.95)],
                forward_vols: vec![T::from_f64(0.2), T::from_f64(0.2)],
                corr_decay_beta: T::from_f64(0.1),
                option_expiry: T::one(),
                swaption_strike: T::from_f64(0.03),
            }
        }
        #[test]
        fn suite45_rebonato_smoke() {
            let r = evaluate_brigo_mercurio_rebonato_swaption_and_cms(&spec(0.));
            assert!(
                r.weights_sum_conservation_error < 1e-12 && r.cms_convexity_adjustment_bps > 0.
            );
            let h = 1e-6;
            let fd = (evaluate_brigo_mercurio_rebonato_swaption_and_cms(&spec(h))
                .black76_payer_swaption_pv
                - evaluate_brigo_mercurio_rebonato_swaption_and_cms(&spec(-h))
                    .black76_payer_swaption_pv)
                / (2. * h);
            let ad = evaluate_brigo_mercurio_rebonato_swaption_and_cms(&spec(Dual::<1>::variable(
                0., 0,
            )))
            .black76_payer_swaption_pv;
            assert!((ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_brigo_mercurio_rebonato_swaption_and_cms_module::*;

mod evaluate_brigo_tuckman_nss_and_key_rate_durations_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct BrigoTuckmanNssKeyRateDurationSpec<T: Real = f64> {
        pub beta0_level: T,
        pub beta1_slope: T,
        pub beta2_hump1: T,
        pub beta3_hump2: T,
        pub tau1_decay: T,
        pub tau2_decay: T,
        pub eval_maturity: T,
        pub cashflow_weights_2y_5y_10y: [T; 3],
    }

    #[derive(Clone, Copy, Debug)]
    pub struct BrigoTuckmanNssKeyRateDurationReport<T: Real = f64> {
        pub short_rate_limit_r0: T,
        pub asymptotic_long_rate: T,
        pub nss_zero_yield_at_eval: T,
        pub nss_instantaneous_forward_at_eval: T,
        pub forward_zero_differential_error: T,
        pub key_rate_durations_2y_5y_10y: [T; 3],
        pub total_fisher_weil_duration: T,
        pub krd_sum_conservation_error: T,
    }

    fn nss_zero_and_fwd<T: Real>(b0: T, b1: T, b2: T, b3: T, tau1: T, tau2: T, m: T) -> (T, T) {
        let u1 = m / tau1.max_val(T::from_f64(1e-8));
        let u2 = m / tau2.max_val(T::from_f64(1e-8));
        let e1 = (-u1).exp();
        let e2 = (-u2).exp();
        let l1_1 = (T::one() - e1) / u1.max_val(T::from_f64(1e-8));
        let l2_1 = l1_1 - e1;
        let l1_2 = (T::one() - e2) / u2.max_val(T::from_f64(1e-8));
        let l2_2 = l1_2 - e2;
        let y = b0 + b1 * l1_1 + b2 * l2_1 + b3 * l2_2;
        let f = b0 + b1 * e1 + b2 * u1 * e1 + b3 * u2 * e2;
        (y, f)
    }

    pub fn evaluate_brigo_tuckman_nss_and_key_rate_durations<T: Real>(
        spec: &BrigoTuckmanNssKeyRateDurationSpec<T>,
    ) -> BrigoTuckmanNssKeyRateDurationReport<T> {
        let m = spec.eval_maturity.max_val(T::from_f64(1e-6));
        let (y_eval, f_eval) = nss_zero_and_fwd(
            spec.beta0_level,
            spec.beta1_slope,
            spec.beta2_hump1,
            spec.beta3_hump2,
            spec.tau1_decay,
            spec.tau2_decay,
            m,
        );
        let dm = T::from_f64(1e-6);
        let (y_p, _) = nss_zero_and_fwd(
            spec.beta0_level,
            spec.beta1_slope,
            spec.beta2_hump1,
            spec.beta3_hump2,
            spec.tau1_decay,
            spec.tau2_decay,
            m + dm,
        );
        let (y_m, _) = nss_zero_and_fwd(
            spec.beta0_level,
            spec.beta1_slope,
            spec.beta2_hump1,
            spec.beta3_hump2,
            spec.tau1_decay,
            spec.tau2_decay,
            m - dm,
        );
        let d_my_dm = ((m + dm) * y_p - (m - dm) * y_m) / (T::from_f64(2.0) * dm);
        let diff_f = f_eval - d_my_dm;
        let fwd_err = (diff_f * diff_f + T::from_f64(1e-30)).sqrt();

        let tenors = [T::from_f64(2.0), T::from_f64(5.0), T::from_f64(10.0)];
        let mut pvs = [T::zero(); 3];
        let mut p_tot = T::zero();
        for k in 0..3 {
            let (yk, _) = nss_zero_and_fwd(
                spec.beta0_level,
                spec.beta1_slope,
                spec.beta2_hump1,
                spec.beta3_hump2,
                spec.tau1_decay,
                spec.tau2_decay,
                tenors[k],
            );
            pvs[k] = spec.cashflow_weights_2y_5y_10y[k] * (-tenors[k] * yk).exp();
            p_tot = p_tot + pvs[k];
        }
        let p_denom = p_tot.max_val(T::from_f64(1e-12));
        let krds = [
            tenors[0] * pvs[0] / p_denom,
            tenors[1] * pvs[1] / p_denom,
            tenors[2] * pvs[2] / p_denom,
        ];
        let fw_dur = (tenors[0] * pvs[0] + tenors[1] * pvs[1] + tenors[2] * pvs[2]) / p_denom;
        let d_krd = (krds[0] + krds[1] + krds[2]) - fw_dur;
        let krd_err = (d_krd * d_krd + T::from_f64(1e-30)).sqrt();

        BrigoTuckmanNssKeyRateDurationReport {
            short_rate_limit_r0: spec.beta0_level + spec.beta1_slope,
            asymptotic_long_rate: spec.beta0_level,
            nss_zero_yield_at_eval: y_eval,
            nss_instantaneous_forward_at_eval: f_eval,
            forward_zero_differential_error: fwd_err,
            key_rate_durations_2y_5y_10y: krds,
            total_fisher_weil_duration: fw_dur,
            krd_sum_conservation_error: krd_err,
        }
    }

    #[cfg(test)]
    mod domain_suite_49_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn brigo_tuckman_nss_and_krd_invariants_and_ad_parity() {
            let spec = BrigoTuckmanNssKeyRateDurationSpec {
                beta0_level: 0.045,
                beta1_slope: -0.015,
                beta2_hump1: 0.018,
                beta3_hump2: -0.008,
                tau1_decay: 1.8,
                tau2_decay: 5.5,
                eval_maturity: 4.0,
                cashflow_weights_2y_5y_10y: [0.25, 0.35, 0.40],
            };
            let rep = evaluate_brigo_tuckman_nss_and_key_rate_durations(&spec);
            assert!(rep.forward_zero_differential_error < 1e-9);
            assert!(rep.krd_sum_conservation_error < 1e-12);

            let h = 1e-6;
            let r_ad = evaluate_brigo_tuckman_nss_and_key_rate_durations(
                &BrigoTuckmanNssKeyRateDurationSpec {
                    beta0_level: Dual::<1>::variable(spec.beta0_level, 0),
                    beta1_slope: Dual::<1>::from_value(spec.beta1_slope),
                    beta2_hump1: Dual::<1>::from_value(spec.beta2_hump1),
                    beta3_hump2: Dual::<1>::from_value(spec.beta3_hump2),
                    tau1_decay: Dual::<1>::from_value(spec.tau1_decay),
                    tau2_decay: Dual::<1>::from_value(spec.tau2_decay),
                    eval_maturity: Dual::<1>::from_value(spec.eval_maturity),
                    cashflow_weights_2y_5y_10y: spec
                        .cashflow_weights_2y_5y_10y
                        .map(Dual::<1>::from_value),
                },
            );
            let mut sp = spec;
            sp.beta0_level += h;
            let mut sm = spec;
            sm.beta0_level -= h;
            let fd = (evaluate_brigo_tuckman_nss_and_key_rate_durations(&sp)
                .total_fisher_weil_duration
                - evaluate_brigo_tuckman_nss_and_key_rate_durations(&sm)
                    .total_fisher_weil_duration)
                / (2.0 * h);
            assert!((r_ad.total_fisher_weil_duration.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_brigo_tuckman_nss_and_key_rate_durations_module::*;

mod evaluate_andersen_piterbarg_cms_replication_convexity_module {
    // Andersen & Piterbarg (2010) & Brigo–Mercurio (2006) Constant Maturity Swap (`CMS`)
    // Cash-vs-Physical Annuity Ratio & Swaption-Replication Convexity Adjustment Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct AndersenPiterbargCmsReplicationSpec<T: Real = f64> {
        pub forward_swap_rate_s0: T,
        pub swap_tenor_years_n: T,
        pub payment_frequency_m: T,
        pub fixing_expiry_te: T,
        pub normal_swaption_vol_sigma_n: T,
        pub ois_spread_delta_r: T,
        pub mean_reversion_speed_a: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct AndersenPiterbargCmsReplicationReport<T: Real = f64> {
        pub cash_settled_annuity_gc: T,
        pub physical_ois_annuity_ap: T,
        pub annuity_derivative_ratio_g1: T,
        pub annuity_second_deriv_ratio_g2: T,
        pub replication_cms_convexity_adj_bps: T,
        pub convexity_adjusted_cms_rate: T,
        pub in_arrears_swap_rate_convexity_bps: T,
    }

    fn cash_annuity<T: Real>(s: T, n_years: T, m_freq: T) -> T {
        let s_pos = s.max_val(T::from_f64(1e-4));
        let m = m_freq.max_val(T::one());
        let total_periods = m * n_years.max_val(T::one());
        let base = T::one() + s_pos / m;
        (T::one() - base.pow(-total_periods)) / s_pos
    }

    pub fn evaluate_andersen_piterbarg_cms_replication_convexity<T: Real>(
        spec: &AndersenPiterbargCmsReplicationSpec<T>,
    ) -> AndersenPiterbargCmsReplicationReport<T> {
        let s0 = spec.forward_swap_rate_s0.max_val(T::from_f64(1e-4));
        let n = spec.swap_tenor_years_n.max_val(T::one());
        let m = spec.payment_frequency_m.max_val(T::one());
        let te = spec.fixing_expiry_te.max_val(T::from_f64(0.05));
        let sig_n = spec.normal_swaption_vol_sigma_n.max_val(T::from_f64(1e-5));

        let gc = cash_annuity(s0, n, m);
        let s_ois = (s0 - spec.ois_spread_delta_r).max_val(T::from_f64(1e-4));
        let ap = cash_annuity(s_ois, n, m);

        // Analytic finite-step derivatives of cash annuity G(S) in T: Real:
        let ds = T::from_f64(1e-4);
        let g_up = cash_annuity(s0 + ds, n, m);
        let g_dn = cash_annuity(s0 - ds, n, m);
        let g_prime = (g_up - g_dn) / (T::from_f64(2.0) * ds);
        let g_second = (g_up - T::from_f64(2.0) * gc + g_dn) / (ds * ds);

        // Andersen-Piterbarg (2010, Vol. I Ch. 16) annuity-slope ratio -G'(S_0)/G(S_0) > 0:
        let g1_ratio = (-g_prime / gc).max_val(T::from_f64(1e-4));
        let g2_ratio = (g_second / gc).max_val(T::from_f64(1e-4));

        // Mean-reversion damping of swap-rate variance over [0, T_e]:
        let a = spec.mean_reversion_speed_a.max_val(T::from_f64(1e-3));
        let mr_factor = (T::one() - (T::from_f64(-2.0) * a * te).exp()) / (T::from_f64(2.0) * a);
        let var_s = sig_n * sig_n * mr_factor * (ap / gc);

        // CMS convexity adjustment = (-G'(S_0) / G(S_0)) * Var(S_{T_e}):
        let cms_conv = g1_ratio * var_s;
        let cms_rate = s0 + cms_conv;
        let in_arrears_bps =
            (cms_conv + T::from_f64(0.5) * g2_ratio * s0 * var_s) * T::from_f64(1e4);

        AndersenPiterbargCmsReplicationReport {
            cash_settled_annuity_gc: gc,
            physical_ois_annuity_ap: ap,
            annuity_derivative_ratio_g1: g1_ratio,
            annuity_second_deriv_ratio_g2: g2_ratio,
            replication_cms_convexity_adj_bps: cms_conv * T::from_f64(1e4),
            convexity_adjusted_cms_rate: cms_rate,
            in_arrears_swap_rate_convexity_bps: in_arrears_bps,
        }
    }

    #[cfg(test)]
    mod domain_suite_53_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn andersen_piterbarg_cms_replication_invariants_and_ad_parity() {
            let spec = AndersenPiterbargCmsReplicationSpec {
                forward_swap_rate_s0: 0.038,
                swap_tenor_years_n: 10.0,
                payment_frequency_m: 2.0,
                fixing_expiry_te: 5.0,
                normal_swaption_vol_sigma_n: 0.0085,
                ois_spread_delta_r: 0.0015,
                mean_reversion_speed_a: 0.04,
            };
            let rep = evaluate_andersen_piterbarg_cms_replication_convexity(&spec);
            assert!(rep.cash_settled_annuity_gc > 0.0);
            assert!(rep.replication_cms_convexity_adj_bps > 0.0);
            assert!(rep.convexity_adjusted_cms_rate > spec.forward_swap_rate_s0);

            let h = 1e-6;
            let r_ad = evaluate_andersen_piterbarg_cms_replication_convexity(
                &AndersenPiterbargCmsReplicationSpec {
                    forward_swap_rate_s0: Dual::<1>::variable(spec.forward_swap_rate_s0, 0),
                    swap_tenor_years_n: Dual::<1>::from_value(spec.swap_tenor_years_n),
                    payment_frequency_m: Dual::<1>::from_value(spec.payment_frequency_m),
                    fixing_expiry_te: Dual::<1>::from_value(spec.fixing_expiry_te),
                    normal_swaption_vol_sigma_n: Dual::<1>::from_value(
                        spec.normal_swaption_vol_sigma_n,
                    ),
                    ois_spread_delta_r: Dual::<1>::from_value(spec.ois_spread_delta_r),
                    mean_reversion_speed_a: Dual::<1>::from_value(spec.mean_reversion_speed_a),
                },
            );
            let mut sp = spec;
            sp.forward_swap_rate_s0 += h;
            let mut sm = spec;
            sm.forward_swap_rate_s0 -= h;
            let fd = (evaluate_andersen_piterbarg_cms_replication_convexity(&sp)
                .convexity_adjusted_cms_rate
                - evaluate_andersen_piterbarg_cms_replication_convexity(&sm)
                    .convexity_adjusted_cms_rate)
                / (2.0 * h);
            assert!((r_ad.convexity_adjusted_cms_rate.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_andersen_piterbarg_cms_replication_convexity_module::*;

mod evaluate_litterman_scheinkman_butterfly_pca_carry_rolldown_module {
    // Litterman & Scheinkman (1991, J. Fixed Income) & Piazzesi (2010) Yield Curve
    // Principal-Component Butterfly (Duration & PCA-Slope Neutral Barbell) & Carry-Roll-Down Convexity Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct LittermanScheinkmanButterflyPcaSpec<T: Real = f64> {
        pub short_wing_duration_d1: T,
        pub belly_duration_d2: T,
        pub long_wing_duration_d3: T,
        pub short_wing_slope_loading_s1: T,
        pub belly_slope_loading_s2: T,
        pub long_wing_slope_loading_s3: T,
        pub short_wing_convexity_c1: T,
        pub belly_convexity_c2: T,
        pub long_wing_convexity_c3: T,
        pub parallel_yield_shift_bps: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct LittermanScheinkmanButterflyPcaReport<T: Real = f64> {
        pub optimal_short_wing_weight_w1: T,
        pub optimal_long_wing_weight_w3: T,
        pub duration_neutral_fifty_fifty_w1: T,
        pub duration_neutral_fifty_fifty_w3: T,
        pub net_butterfly_convexity_exposure: T,
        pub convexity_pnl_advantage_bps: T,
        pub residual_duration_exposure: T,
    }

    pub fn evaluate_litterman_scheinkman_butterfly_pca_carry_rolldown<T: Real>(
        spec: &LittermanScheinkmanButterflyPcaSpec<T>,
    ) -> LittermanScheinkmanButterflyPcaReport<T> {
        let d1 = spec.short_wing_duration_d1.max_val(T::from_f64(0.5));
        let d2 = spec.belly_duration_d2.max_val(T::from_f64(1.0));
        let d3 = spec.long_wing_duration_d3.max_val(T::from_f64(2.0));
        let s1 = spec.short_wing_slope_loading_s1;
        let s2 = spec.belly_slope_loading_s2;
        let s3 = spec.long_wing_slope_loading_s3;
        let dy = spec.parallel_yield_shift_bps * T::from_f64(1e-4);

        // Solve 2x2 system for unit short belly (-1 unit of leg 2):
        // w1 * d1 + w3 * d3 = d2  (duration neutrality)
        // w1 * d1 * s1 + w3 * d3 * s3 = d2 * s2  (PCA slope-factor dollar neutrality)
        let ds31 = (s3 - s1).max_val(T::from_f64(0.05));
        let w1 = (d2 * (s3 - s2)) / (d1 * ds31);
        let w3 = (d2 * (s2 - s1)) / (d3 * ds31);

        // Standard 50-50 duration-weighted butterfly weights:
        let w1_50 = T::from_f64(0.5) * d2 / d1;
        let w3_50 = T::from_f64(0.5) * d2 / d3;

        let net_conv = w1 * spec.short_wing_convexity_c1 - spec.belly_convexity_c2
            + w3 * spec.long_wing_convexity_c3;
        let conv_pnl_bps = T::from_f64(0.5) * net_conv * dy * dy * T::from_f64(10_000.0);
        let res_dur = w1 * d1 - d2 + w3 * d3;

        LittermanScheinkmanButterflyPcaReport {
            optimal_short_wing_weight_w1: w1,
            optimal_long_wing_weight_w3: w3,
            duration_neutral_fifty_fifty_w1: w1_50,
            duration_neutral_fifty_fifty_w3: w3_50,
            net_butterfly_convexity_exposure: net_conv,
            convexity_pnl_advantage_bps: conv_pnl_bps,
            residual_duration_exposure: res_dur,
        }
    }

    #[cfg(test)]
    mod evaluate_litterman_scheinkman_butterfly_pca_carry_rolldown_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn litterman_scheinkman_butterfly_pca_invariants_and_ad_parity() {
            let spec = LittermanScheinkmanButterflyPcaSpec {
                short_wing_duration_d1: 1.92,
                belly_duration_d2: 4.55,
                long_wing_duration_d3: 8.40,
                short_wing_slope_loading_s1: -0.65,
                belly_slope_loading_s2: -0.10,
                long_wing_slope_loading_s3: 0.48,
                short_wing_convexity_c1: 4.2,
                belly_convexity_c2: 23.5,
                long_wing_convexity_c3: 82.0,
                parallel_yield_shift_bps: 50.0,
            };
            let rep = evaluate_litterman_scheinkman_butterfly_pca_carry_rolldown(&spec);
            assert!(rep.optimal_short_wing_weight_w1 > 0.0);
            assert!(rep.optimal_long_wing_weight_w3 > 0.0);
            assert!(rep.net_butterfly_convexity_exposure > 0.0);
            assert!(rep.residual_duration_exposure.abs() < 1e-6);

            let h = 1e-5;
            let r_ad = evaluate_litterman_scheinkman_butterfly_pca_carry_rolldown(
                &LittermanScheinkmanButterflyPcaSpec {
                    short_wing_duration_d1: Dual::<10>::variable(spec.short_wing_duration_d1, 0),
                    belly_duration_d2: Dual::<10>::variable(spec.belly_duration_d2, 1),
                    long_wing_duration_d3: Dual::<10>::variable(spec.long_wing_duration_d3, 2),
                    short_wing_slope_loading_s1: Dual::<10>::variable(
                        spec.short_wing_slope_loading_s1,
                        3,
                    ),
                    belly_slope_loading_s2: Dual::<10>::variable(spec.belly_slope_loading_s2, 4),
                    long_wing_slope_loading_s3: Dual::<10>::variable(
                        spec.long_wing_slope_loading_s3,
                        5,
                    ),
                    short_wing_convexity_c1: Dual::<10>::variable(spec.short_wing_convexity_c1, 6),
                    belly_convexity_c2: Dual::<10>::variable(spec.belly_convexity_c2, 7),
                    long_wing_convexity_c3: Dual::<10>::variable(spec.long_wing_convexity_c3, 8),
                    parallel_yield_shift_bps: Dual::<10>::variable(
                        spec.parallel_yield_shift_bps,
                        9,
                    ),
                },
            );
            for (idx, &ad) in r_ad.convexity_pnl_advantage_bps.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.short_wing_duration_d1 += h;
                        sm.short_wing_duration_d1 -= h;
                    }
                    1 => {
                        sp.belly_duration_d2 += h;
                        sm.belly_duration_d2 -= h;
                    }
                    2 => {
                        sp.long_wing_duration_d3 += h;
                        sm.long_wing_duration_d3 -= h;
                    }
                    3 => {
                        sp.short_wing_slope_loading_s1 += h;
                        sm.short_wing_slope_loading_s1 -= h;
                    }
                    4 => {
                        sp.belly_slope_loading_s2 += h;
                        sm.belly_slope_loading_s2 -= h;
                    }
                    5 => {
                        sp.long_wing_slope_loading_s3 += h;
                        sm.long_wing_slope_loading_s3 -= h;
                    }
                    6 => {
                        sp.short_wing_convexity_c1 += h;
                        sm.short_wing_convexity_c1 -= h;
                    }
                    7 => {
                        sp.belly_convexity_c2 += h;
                        sm.belly_convexity_c2 -= h;
                    }
                    8 => {
                        sp.long_wing_convexity_c3 += h;
                        sm.long_wing_convexity_c3 -= h;
                    }
                    _ => {
                        sp.parallel_yield_shift_bps += h;
                        sm.parallel_yield_shift_bps -= h;
                    }
                }
                let fd = (evaluate_litterman_scheinkman_butterfly_pca_carry_rolldown(&sp)
                    .convexity_pnl_advantage_bps
                    - evaluate_litterman_scheinkman_butterfly_pca_carry_rolldown(&sm)
                        .convexity_pnl_advantage_bps)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_litterman_scheinkman_butterfly_pca_carry_rolldown_module::*;

mod evaluate_gurkaynak_sack_wright_svensson_treasury_forward_curve_module {
    // Gürkaynak, Sack & Wright (2007, J. Monetary Economics, "The U.S. Treasury Yield Curve:
    // 1961 to the Present") & Svensson (1994, NBER) Off-the-Run Svensson Instantaneous Forward Curve,
    // Two-Hump Curvature Decomposition & Off-the-Run Liquidity Spread Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct GswSvenssonTreasuryCurveSpec<T: Real = f64> {
        pub svensson_beta0: T,
        pub svensson_beta1: T,
        pub svensson_beta2: T,
        pub svensson_beta3: T,
        pub svensson_tau1: T,
        pub svensson_tau2: T,
        pub maturity_years_t: T,
        pub on_the_run_par_yield_y_otr: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GswSvenssonTreasuryCurveReport<T: Real = f64> {
        pub svensson_zero_coupon_yield: T,
        pub svensson_instantaneous_forward_rate: T,
        pub forward_minus_zero_slope_wedge_bps: T,
        pub second_hump_curvature_contribution: T,
        pub off_the_run_implied_par_yield: T,
        pub off_the_run_liquidity_spread_bps: T,
    }

    fn svensson_zero_yield_inner<T: Real>(
        beta0: T,
        beta1: T,
        beta2: T,
        beta3: T,
        tau1: T,
        tau2: T,
        t: T,
    ) -> (T, T, T) {
        let x1 = t / tau1;
        let x2 = t / tau2;
        let e1 = (-x1).exp();
        let e2 = (-x2).exp();
        let l1 = (T::one() - e1) / x1;
        let c1 = l1 - e1;
        let c2 = (T::one() - e2) / x2 - e2;
        let y_zero = beta0 + beta1 * l1 + beta2 * c1 + beta3 * c2;
        let f_inst = beta0 + beta1 * e1 + beta2 * x1 * e1 + beta3 * x2 * e2;
        let hump2 = beta3 * c2;
        (y_zero, f_inst, hump2)
    }

    pub fn evaluate_gurkaynak_sack_wright_svensson_treasury_forward_curve<T: Real>(
        spec: &GswSvenssonTreasuryCurveSpec<T>,
    ) -> GswSvenssonTreasuryCurveReport<T> {
        let b0 = spec.svensson_beta0;
        let b1 = spec.svensson_beta1;
        let b2 = spec.svensson_beta2;
        let b3 = spec.svensson_beta3;
        let tau1 = spec.svensson_tau1.max_val(T::from_f64(0.10));
        let tau2 = spec.svensson_tau2.max_val(T::from_f64(0.25));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.50));
        let y_otr = spec.on_the_run_par_yield_y_otr;

        let (y_zero, f_inst, hump2) = svensson_zero_yield_inner(b0, b1, b2, b3, tau1, tau2, t);
        const GL_NODES: [f64; 4] = [
            -0.861_136_311_594_052_6,
            -0.339_981_043_584_856_26,
            0.339_981_043_584_856_26,
            0.861_136_311_594_052_6,
        ];
        const GL_WEIGHTS: [f64; 4] = [
            0.347_854_845_137_453_85,
            0.652_145_154_862_546_1,
            0.652_145_154_862_546_1,
            0.347_854_845_137_453_85,
        ];
        let half_t = T::from_f64(0.5) * t;
        let mut annuity = T::zero();
        for (&xi, &wi) in GL_NODES.iter().zip(GL_WEIGHTS.iter()) {
            let u = half_t * (T::one() + T::from_f64(xi));
            let (y_u, _, _) = svensson_zero_yield_inner(b0, b1, b2, b3, tau1, tau2, u);
            annuity += half_t * T::from_f64(wi) * (-y_u * u).exp();
        }
        let p_full = (-y_zero * t).exp();
        let c_par = (T::one() - p_full) / annuity.max_val(T::from_f64(1e-6));

        let slope_wedge_bps = (f_inst - y_zero) * T::from_f64(10_000.0);
        let liq_spread_bps = (c_par - y_otr) * T::from_f64(10_000.0);

        GswSvenssonTreasuryCurveReport {
            svensson_zero_coupon_yield: y_zero,
            svensson_instantaneous_forward_rate: f_inst,
            forward_minus_zero_slope_wedge_bps: slope_wedge_bps,
            second_hump_curvature_contribution: hump2,
            off_the_run_implied_par_yield: c_par,
            off_the_run_liquidity_spread_bps: liq_spread_bps,
        }
    }

    #[cfg(test)]
    mod evaluate_gurkaynak_sack_wright_svensson_treasury_forward_curve_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gsw_svensson_treasury_forward_curve_invariants_and_ad_parity() {
            let spec = GswSvenssonTreasuryCurveSpec {
                svensson_beta0: 0.045,
                svensson_beta1: -0.012,
                svensson_beta2: 0.018,
                svensson_beta3: -0.010,
                svensson_tau1: 1.80,
                svensson_tau2: 6.50,
                maturity_years_t: 10.0,
                on_the_run_par_yield_y_otr: 0.0415,
            };
            let rep = evaluate_gurkaynak_sack_wright_svensson_treasury_forward_curve(&spec);
            assert!(rep.svensson_zero_coupon_yield > 0.02 && rep.svensson_zero_coupon_yield < 0.08);
            assert!(rep.svensson_instantaneous_forward_rate > 0.02);
            assert!(rep.second_hump_curvature_contribution < 0.0);
            assert!(
                (rep.off_the_run_implied_par_yield - rep.svensson_zero_coupon_yield).abs() < 0.005
            );

            let h = 1e-5;
            let r_ad = evaluate_gurkaynak_sack_wright_svensson_treasury_forward_curve(
                &GswSvenssonTreasuryCurveSpec {
                    svensson_beta0: Dual::<8>::variable(spec.svensson_beta0, 0),
                    svensson_beta1: Dual::<8>::variable(spec.svensson_beta1, 1),
                    svensson_beta2: Dual::<8>::variable(spec.svensson_beta2, 2),
                    svensson_beta3: Dual::<8>::variable(spec.svensson_beta3, 3),
                    svensson_tau1: Dual::<8>::variable(spec.svensson_tau1, 4),
                    svensson_tau2: Dual::<8>::variable(spec.svensson_tau2, 5),
                    maturity_years_t: Dual::<8>::variable(spec.maturity_years_t, 6),
                    on_the_run_par_yield_y_otr: Dual::<8>::variable(
                        spec.on_the_run_par_yield_y_otr,
                        7,
                    ),
                },
            );
            for (idx, &ad) in r_ad.off_the_run_liquidity_spread_bps.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.svensson_beta0 += h;
                        sm.svensson_beta0 -= h;
                    }
                    1 => {
                        sp.svensson_beta1 += h;
                        sm.svensson_beta1 -= h;
                    }
                    2 => {
                        sp.svensson_beta2 += h;
                        sm.svensson_beta2 -= h;
                    }
                    3 => {
                        sp.svensson_beta3 += h;
                        sm.svensson_beta3 -= h;
                    }
                    4 => {
                        sp.svensson_tau1 += h;
                        sm.svensson_tau1 -= h;
                    }
                    5 => {
                        sp.svensson_tau2 += h;
                        sm.svensson_tau2 -= h;
                    }
                    6 => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                    _ => {
                        sp.on_the_run_par_yield_y_otr += h;
                        sm.on_the_run_par_yield_y_otr -= h;
                    }
                }
                let fd = (evaluate_gurkaynak_sack_wright_svensson_treasury_forward_curve(&sp)
                    .off_the_run_liquidity_spread_bps
                    - evaluate_gurkaynak_sack_wright_svensson_treasury_forward_curve(&sm)
                        .off_the_run_liquidity_spread_bps)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_gurkaynak_sack_wright_svensson_treasury_forward_curve_module::*;

mod evaluate_joslin_priebsch_singleton_unspanned_macro_term_structure_module {
    // Joslin, Priebsch & Singleton (2014, Journal of Finance 69(3), "Risk Premiums in Dynamic
    // Term Structure Models with Unspanned Macro Risks"), Duffee (2011, RFS 24(9)),
    // & Piazzesi (2010) Unspanned Macro Risk Term Structure & Bond Excess Return Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct JoslinPriebschSingletonUnspannedMacroSpec<T: Real = f64> {
        pub yield_curve_level_l: T,
        pub yield_curve_slope_s: T,
        pub yield_curve_curvature_c: T,
        pub real_output_growth_g: T,
        pub core_inflation_gap_pi: T,
        pub macro_on_yields_projection_beta: T,
        pub spanned_slope_risk_price_lambda_p: T,
        pub unspanned_macro_risk_price_lambda_m: T,
        pub q_measure_decay_kappa_q: T,
        pub bond_maturity_years_n: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct JoslinPriebschSingletonUnspannedMacroReport<T: Real = f64> {
        pub spanned_model_zero_yield: T,
        pub orthogonal_unspanned_macro_impulse: T,
        pub spanned_yield_curve_bond_risk_premium: T,
        pub unspanned_macro_bond_risk_premium: T,
        pub total_expected_bond_excess_return: T,
        pub unspanned_risk_premium_variance_share: T,
    }

    pub fn evaluate_joslin_priebsch_singleton_unspanned_macro_term_structure<T: Real>(
        spec: &JoslinPriebschSingletonUnspannedMacroSpec<T>,
    ) -> JoslinPriebschSingletonUnspannedMacroReport<T> {
        let l_pc = spec.yield_curve_level_l;
        let s_pc = spec.yield_curve_slope_s;
        let c_pc = spec.yield_curve_curvature_c;
        let g_macro = spec.real_output_growth_g;
        let pi_macro = spec.core_inflation_gap_pi;
        let beta_proj = spec.macro_on_yields_projection_beta;
        let lam_p = spec.spanned_slope_risk_price_lambda_p;
        let lam_m = spec.unspanned_macro_risk_price_lambda_m;
        let kappa_q = spec.q_measure_decay_kappa_q.max_val(T::from_f64(0.05));
        let n = spec.bond_maturity_years_n.max_val(T::from_f64(1.25));

        let kn = kappa_q * n;
        let e_kn = (-kn).exp();
        let b_slope = (T::one() - e_kn) / kn;
        let b_curv = b_slope - e_kn;

        // Under JPS (2014), the Q-measure yield curve is spanned solely by (L, S, C):
        let y_spanned = l_pc + b_slope * s_pc + b_curv * c_pc;

        // Unspanned macro factor orthogonal to the spanned yield-curve principal components:
        let raw_macro = T::from_f64(0.60) * pi_macro + T::from_f64(0.40) * g_macro;
        let spanned_macro_fit =
            beta_proj * (l_pc + T::from_f64(0.5) * s_pc + T::from_f64(0.25) * c_pc);
        let m_orth = raw_macro - spanned_macro_fit;

        let nm1 = n - T::one();
        let knm1 = kappa_q * nm1;
        let e_knm1 = (-knm1).exp();
        let b_slope_nm1 = (T::one() - e_knm1) / knm1;
        let b_curv_nm1 = b_slope_nm1 - e_knm1;

        // 1-year expected bond excess return decomposes into spanned (slope/curvature) + unspanned macro risk premia:
        let rp_spanned = -nm1 * (b_slope_nm1 * s_pc - T::from_f64(0.5) * b_curv_nm1 * c_pc) * lam_p;
        let rp_unspanned = -nm1 * lam_m * m_orth;
        let rp_total = rp_spanned + rp_unspanned;

        let var_sp = rp_spanned * rp_spanned + T::from_f64(1e-8);
        let var_un = rp_unspanned * rp_unspanned;
        let share_un = var_un / (var_sp + var_un);

        JoslinPriebschSingletonUnspannedMacroReport {
            spanned_model_zero_yield: y_spanned,
            orthogonal_unspanned_macro_impulse: m_orth,
            spanned_yield_curve_bond_risk_premium: rp_spanned,
            unspanned_macro_bond_risk_premium: rp_unspanned,
            total_expected_bond_excess_return: rp_total,
            unspanned_risk_premium_variance_share: share_un,
        }
    }

    #[cfg(test)]
    mod evaluate_joslin_priebsch_singleton_unspanned_macro_term_structure_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn jps_unspanned_macro_term_structure_invariants_and_ad_parity() {
            let spec = JoslinPriebschSingletonUnspannedMacroSpec {
                yield_curve_level_l: 0.042,
                yield_curve_slope_s: -0.014,
                yield_curve_curvature_c: 0.008,
                real_output_growth_g: 0.024,
                core_inflation_gap_pi: 0.015,
                macro_on_yields_projection_beta: 0.25,
                spanned_slope_risk_price_lambda_p: 0.18,
                unspanned_macro_risk_price_lambda_m: -0.22,
                q_measure_decay_kappa_q: 0.45,
                bond_maturity_years_n: 10.0,
            };
            let rep = evaluate_joslin_priebsch_singleton_unspanned_macro_term_structure(&spec);
            assert!(rep.spanned_model_zero_yield > 0.02 && rep.spanned_model_zero_yield < 0.07);
            assert!(
                rep.unspanned_risk_premium_variance_share > 0.0
                    && rep.unspanned_risk_premium_variance_share < 1.0
            );

            let h = 1e-6;
            let r_ad = evaluate_joslin_priebsch_singleton_unspanned_macro_term_structure(
                &JoslinPriebschSingletonUnspannedMacroSpec {
                    yield_curve_level_l: Dual::<6>::variable(spec.yield_curve_level_l, 0),
                    yield_curve_slope_s: Dual::<6>::variable(spec.yield_curve_slope_s, 1),
                    yield_curve_curvature_c: Dual::<6>::variable(spec.yield_curve_curvature_c, 2),
                    real_output_growth_g: Dual::<6>::variable(spec.real_output_growth_g, 3),
                    core_inflation_gap_pi: Dual::<6>::variable(spec.core_inflation_gap_pi, 4),
                    macro_on_yields_projection_beta: Dual::from_f64(
                        spec.macro_on_yields_projection_beta,
                    ),
                    spanned_slope_risk_price_lambda_p: Dual::from_f64(
                        spec.spanned_slope_risk_price_lambda_p,
                    ),
                    unspanned_macro_risk_price_lambda_m: Dual::<6>::variable(
                        spec.unspanned_macro_risk_price_lambda_m,
                        5,
                    ),
                    q_measure_decay_kappa_q: Dual::from_f64(spec.q_measure_decay_kappa_q),
                    bond_maturity_years_n: Dual::from_f64(spec.bond_maturity_years_n),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.yield_curve_level_l,
                    spec.yield_curve_slope_s,
                    spec.yield_curve_curvature_c,
                    spec.real_output_growth_g,
                    spec.core_inflation_gap_pi,
                    spec.unspanned_macro_risk_price_lambda_m,
                ];
                v[i] += h;
                let up = evaluate_joslin_priebsch_singleton_unspanned_macro_term_structure(
                    &JoslinPriebschSingletonUnspannedMacroSpec {
                        yield_curve_level_l: v[0],
                        yield_curve_slope_s: v[1],
                        yield_curve_curvature_c: v[2],
                        real_output_growth_g: v[3],
                        core_inflation_gap_pi: v[4],
                        unspanned_macro_risk_price_lambda_m: v[5],
                        ..spec
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_joslin_priebsch_singleton_unspanned_macro_term_structure(
                    &JoslinPriebschSingletonUnspannedMacroSpec {
                        yield_curve_level_l: v[0],
                        yield_curve_slope_s: v[1],
                        yield_curve_curvature_c: v[2],
                        real_output_growth_g: v[3],
                        core_inflation_gap_pi: v[4],
                        unspanned_macro_risk_price_lambda_m: v[5],
                        ..spec
                    },
                );
                let val_ad = r_ad.spanned_model_zero_yield + r_ad.total_expected_bond_excess_return;
                let fd = ((up.spanned_model_zero_yield + up.total_expected_bond_excess_return)
                    - (dn.spanned_model_zero_yield + dn.total_expected_bond_excess_return))
                    / (2.0 * h);
                assert!((val_ad.d[i] - fd).abs() < 1e-5, "input {i}");
            }
        }
    }
}
pub use evaluate_joslin_priebsch_singleton_unspanned_macro_term_structure_module::*;

use ql_ad::Real;

/// Specification for the Joslin–Singleton–Zhu (2011) canonical Gaussian Dynamic Term Structure Model (GDTSM).
///
/// Reference: Joslin, Singleton & Zhu (2011, *Review of Financial Studies* 24(3), "A New Perspective on
/// Gaussian Dynamic Term Structure Models"), Piazzesi (2010), Duffie & Kan (1996).
#[derive(Debug, Clone, Copy)]
pub struct JoslinSingletonZhuGdtsmSpec<T: Real> {
    /// Long-run risk-neutral short-rate mean $r_\infty^{\mathbb{Q}}$.
    pub risk_neutral_long_run_mean_r_inf_q: T,
    /// Long-run physical short-rate mean $r_\infty^{\mathbb{P}}$.
    pub physical_long_run_mean_r_inf_p: T,
    /// Three risk-neutral mean-reversion decay speeds $(\kappa_1^{\mathbb{Q}}, \kappa_2^{\mathbb{Q}}, \kappa_3^{\mathbb{Q}})$.
    pub risk_neutral_decay_speeds_kappa_q: [T; 3],
    /// Three physical mean-reversion decay speeds $(\kappa_1^{\mathbb{P}}, \kappa_2^{\mathbb{P}}, \kappa_3^{\mathbb{P}})$.
    pub physical_decay_speeds_kappa_p: [T; 3],
    /// Three canonical pricing factors (level, slope, curvature states) $(x_1, x_2, x_3)$.
    pub canonical_factor_states_x: [T; 3],
    /// Three factor innovation volatilities $(\sigma_1, \sigma_2, \sigma_3)$.
    pub factor_volatilities_sigma: [T; 3],
    /// Target bond maturity $\tau$ (in years).
    pub bond_maturity_years_tau: T,
}

/// Analytics returned by [`evaluate_joslin_singleton_zhu_canonical_gdtsm_term_premium`].
#[derive(Debug, Clone, Copy)]
pub struct JoslinSingletonZhuGdtsmAnalytics<T: Real> {
    /// Model-implied affine zero-coupon bond yield $y^{\mathbb{Q}}(\tau)$ including convexity.
    pub canonical_zero_coupon_yield_q: T,
    /// Expectations-hypothesis average physical short rate $y^{\mathbb{P},\text{EH}}(\tau)$.
    pub expectations_hypothesis_yield_p: T,
    /// Canonical yield term premium $TP(\tau) = y^{\mathbb{Q}}(\tau) - y^{\mathbb{P},\text{EH}}(\tau)$.
    pub canonical_yield_term_premium: T,
    /// Instantaneous maximum Sharpe ratio (market price of risk norm) $\|\Lambda_t\|_2$.
    pub instantaneous_market_price_of_risk_norm: T,
}

/// Evaluates the Joslin–Singleton–Zhu (2011) canonical 3-factor GDTSM yield, expectations-hypothesis path,
/// and term premium decomposition.
pub fn evaluate_joslin_singleton_zhu_canonical_gdtsm_term_premium<T: Real>(
    spec: &JoslinSingletonZhuGdtsmSpec<T>,
) -> JoslinSingletonZhuGdtsmAnalytics<T> {
    let tau = spec.bond_maturity_years_tau.max_val(T::from_f64(1e-4));
    let r_inf_q_third = spec.risk_neutral_long_run_mean_r_inf_q / T::from_f64(3.0);
    let r_inf_p_third = spec.physical_long_run_mean_r_inf_p / T::from_f64(3.0);

    let mut log_bond_price_q = T::zero();
    let mut eh_yield_p = T::zero();
    let mut lambda_sq_sum = T::zero();

    for i in 0..3 {
        let kq = spec.risk_neutral_decay_speeds_kappa_q[i].max_val(T::from_f64(1e-4));
        let kp = spec.physical_decay_speeds_kappa_p[i].max_val(T::from_f64(1e-4));
        let xi = spec.canonical_factor_states_x[i];
        let sig = spec.factor_volatilities_sigma[i].max_val(T::from_f64(1e-6));
        let sig2 = sig * sig;

        // Exact Vasicek / JSZ risk-neutral factor loading B_i^Q(tau) = (1 - exp(-kappa_i^Q * tau)) / kappa_i^Q
        let b_q = (T::one() - (-kq * tau).exp()) / kq;
        // Exact affine intercept A_i^Q(tau) including Jensen convexity
        let a_q = (r_inf_q_third - sig2 / (T::from_f64(2.0) * kq * kq)) * (b_q - tau)
            - sig2 * b_q * b_q / (T::from_f64(4.0) * kq);
        log_bond_price_q = log_bond_price_q + a_q - b_q * xi;

        // Physical expectations-hypothesis average short rate over [0, tau]
        let b_p = (T::one() - (-kp * tau).exp()) / kp;
        eh_yield_p = eh_yield_p + r_inf_p_third + (b_p / tau) * (xi - r_inf_p_third);

        // Instantaneous market price of risk Lambda_i(x_i) = (drift_P - drift_Q) / sigma_i
        let drift_p = kp * (r_inf_p_third - xi);
        let drift_q = kq * (r_inf_q_third - xi);
        let lam_i = (drift_p - drift_q) / sig;
        lambda_sq_sum += lam_i * lam_i;
    }

    let y_q = -log_bond_price_q / tau;
    let tp = y_q - eh_yield_p;
    let mpr_norm = lambda_sq_sum.sqrt();

    JoslinSingletonZhuGdtsmAnalytics {
        canonical_zero_coupon_yield_q: y_q,
        expectations_hypothesis_yield_p: eh_yield_p,
        canonical_yield_term_premium: tp,
        instantaneous_market_price_of_risk_norm: mpr_norm,
    }
}

#[cfg(test)]
mod jsz_gdtsm_suite64_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn jsz_canonical_gdtsm_invariants_and_ad_parity() {
        let base = JoslinSingletonZhuGdtsmSpec {
            risk_neutral_long_run_mean_r_inf_q: 0.045_f64,
            physical_long_run_mean_r_inf_p: 0.032,
            risk_neutral_decay_speeds_kappa_q: [0.06, 0.45, 1.20],
            physical_decay_speeds_kappa_p: [0.18, 0.65, 1.45],
            canonical_factor_states_x: [0.015, 0.008, 0.012],
            factor_volatilities_sigma: [0.0065, 0.0080, 0.0095],
            bond_maturity_years_tau: 10.0,
        };
        let r0 = evaluate_joslin_singleton_zhu_canonical_gdtsm_term_premium(&base);
        assert!(r0.canonical_zero_coupon_yield_q > 0.02 && r0.canonical_zero_coupon_yield_q < 0.08);
        assert!(r0.expectations_hypothesis_yield_p > 0.01);
        assert!(r0.instantaneous_market_price_of_risk_norm > 0.0);

        let inputs = [
            base.risk_neutral_long_run_mean_r_inf_q,
            base.physical_long_run_mean_r_inf_p,
            base.risk_neutral_decay_speeds_kappa_q[0],
            base.physical_decay_speeds_kappa_p[0],
            base.canonical_factor_states_x[0],
            base.factor_volatilities_sigma[0],
            base.bond_maturity_years_tau,
        ];
        let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
        for i in 0..7 {
            dual_in[i] = Dual::variable(inputs[i], i);
        }
        let spec_ad = JoslinSingletonZhuGdtsmSpec {
            risk_neutral_long_run_mean_r_inf_q: dual_in[0],
            physical_long_run_mean_r_inf_p: dual_in[1],
            risk_neutral_decay_speeds_kappa_q: [
                dual_in[2],
                Dual::from_f64(base.risk_neutral_decay_speeds_kappa_q[1]),
                Dual::from_f64(base.risk_neutral_decay_speeds_kappa_q[2]),
            ],
            physical_decay_speeds_kappa_p: [
                dual_in[3],
                Dual::from_f64(base.physical_decay_speeds_kappa_p[1]),
                Dual::from_f64(base.physical_decay_speeds_kappa_p[2]),
            ],
            canonical_factor_states_x: [
                dual_in[4],
                Dual::from_f64(base.canonical_factor_states_x[1]),
                Dual::from_f64(base.canonical_factor_states_x[2]),
            ],
            factor_volatilities_sigma: [
                dual_in[5],
                Dual::from_f64(base.factor_volatilities_sigma[1]),
                Dual::from_f64(base.factor_volatilities_sigma[2]),
            ],
            bond_maturity_years_tau: dual_in[6],
        };
        let r_ad = evaluate_joslin_singleton_zhu_canonical_gdtsm_term_premium(&spec_ad);
        let eps = 1e-6;
        for i in 0..7 {
            let mut p = inputs;
            let mut m = inputs;
            p[i] += eps;
            m[i] -= eps;
            let mk = |v: [f64; 7]| JoslinSingletonZhuGdtsmSpec {
                risk_neutral_long_run_mean_r_inf_q: v[0],
                physical_long_run_mean_r_inf_p: v[1],
                risk_neutral_decay_speeds_kappa_q: [
                    v[2],
                    base.risk_neutral_decay_speeds_kappa_q[1],
                    base.risk_neutral_decay_speeds_kappa_q[2],
                ],
                physical_decay_speeds_kappa_p: [
                    v[3],
                    base.physical_decay_speeds_kappa_p[1],
                    base.physical_decay_speeds_kappa_p[2],
                ],
                canonical_factor_states_x: [
                    v[4],
                    base.canonical_factor_states_x[1],
                    base.canonical_factor_states_x[2],
                ],
                factor_volatilities_sigma: [
                    v[5],
                    base.factor_volatilities_sigma[1],
                    base.factor_volatilities_sigma[2],
                ],
                bond_maturity_years_tau: v[6],
            };
            let rp = evaluate_joslin_singleton_zhu_canonical_gdtsm_term_premium(&mk(p));
            let rm = evaluate_joslin_singleton_zhu_canonical_gdtsm_term_premium(&mk(m));
            let fd =
                (rp.canonical_yield_term_premium - rm.canonical_yield_term_premium) / (2.0 * eps);
            assert!(
                r_ad.canonical_yield_term_premium.d[i].abs() > 1e-7,
                "zero derivative at {i}"
            );
            assert!((r_ad.canonical_yield_term_premium.d[i] - fd).abs() < 1e-5);
        }
    }
}

pub mod andersen_piterbarg_lsrm_cms_spread_mod {
    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    /// Specification for the Andersen & Piterbarg (2010, Vol. III, Ch. 16-17) / Hagan (1999)
    /// Linear Swap Rate Model (LSRM) CMS spread call option.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct AndersenPiterbargCmsSpreadSpec<T: Real> {
        /// Forward par swap rate $S_{1,0}$ of the long-tenor CMS leg (e.g. 10Y swap rate).
        pub cms1_forward_swap_rate_s1: T,
        /// Swap tenor $\tau_1$ (years) of the first CMS rate.
        pub cms1_tenor_years_tau1: T,
        /// Normal volatility $\sigma_{n,1}$ of the first swap rate.
        pub cms1_normal_volatility_sigma1: T,
        /// Forward par swap rate $S_{2,0}$ of the short-tenor CMS leg (e.g. 2Y swap rate).
        pub cms2_forward_swap_rate_s2: T,
        /// Swap tenor $\tau_2$ (years) of the second CMS rate.
        pub cms2_tenor_years_tau2: T,
        /// Normal volatility $\sigma_{n,2}$ of the second swap rate.
        pub cms2_normal_volatility_sigma2: T,
        /// Instantaneous correlation $\rho_{12} \in (-1, 1)$ between the two swap rates.
        pub swap_rate_correlation_rho12: T,
        /// Gearing weight $w_1 > 0$ on the first CMS rate.
        pub leg1_gearing_w1: T,
        /// Gearing weight $w_2 > 0$ on the second CMS rate.
        pub leg2_gearing_w2: T,
        /// Spread call strike $K$.
        pub spread_strike_k: T,
        /// Short-rate mean-reversion speed $\kappa > 0$ governing LSRM annuity slope.
        pub mean_reversion_kappa: T,
        /// Continuous OIS discount rate $r_{\text{ois}}$ to fixing/payment date $T_0$.
        pub ois_discount_rate_r: T,
        /// Option fixing and payment horizon $T_0 > 0$ (years).
        pub fixing_horizon_years_t0: T,
    }

    /// Analytical outputs for the Andersen-Piterbarg LSRM CMS spread option engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct AndersenPiterbargCmsSpreadAnalytics<T: Real> {
        /// Convexity-adjusted forward CMS rate $\mu_1 = \mathbb{E}^{T_0}[S_1(T_0)]$.
        pub convexity_adjusted_cms1_forward: T,
        /// Convexity-adjusted forward CMS rate $\mu_2 = \mathbb{E}^{T_0}[S_2(T_0)]$.
        pub convexity_adjusted_cms2_forward: T,
        /// Net LSRM convexity drift wedge $w_1(\mu_1 - S_{1,0}) - w_2(\mu_2 - S_{2,0})$.
        pub cms_spread_convexity_wedge: T,
        /// Discounted Bachelier CMS spread call option present value.
        pub cms_spread_call_present_value: T,
    }

    /// Evaluates the Andersen & Piterbarg (2010) Linear Swap Rate Model (LSRM) CMS spread call option.
    pub fn evaluate_andersen_piterbarg_lsrm_cms_spread_option<T: Real>(
        spec: &AndersenPiterbargCmsSpreadSpec<T>,
    ) -> AndersenPiterbargCmsSpreadAnalytics<T> {
        let t0 = spec.fixing_horizon_years_t0.max_val(T::from_f64(1e-4));
        let kappa = spec.mean_reversion_kappa.max_val(T::from_f64(1e-4));
        let tau1 = spec.cms1_tenor_years_tau1.max_val(T::from_f64(0.25));
        let tau2 = spec.cms2_tenor_years_tau2.max_val(T::from_f64(0.25));

        // Hagan (1999) / Andersen-Piterbarg (2010) LSRM annuity mapping coefficient q_i:
        // d(A(0)/A(S_i))/dS_i|_{S_{i,0}} approximated smoothly via mean-reverting duration factor
        let b1 = (T::one() - (-kappa * tau1).exp()) / (kappa * tau1);
        let b2 = (T::one() - (-kappa * tau2).exp()) / (kappa * tau2);
        let q1 = (tau1 / (T::one() + tau1 * spec.cms1_forward_swap_rate_s1))
            * (T::from_f64(0.5) * (T::one() + b1));
        let q2 = (tau2 / (T::one() + tau2 * spec.cms2_forward_swap_rate_s2))
            * (T::from_f64(0.5) * (T::one() + b2));

        let sig1 = spec
            .cms1_normal_volatility_sigma1
            .max_val(T::from_f64(1e-5));
        let sig2 = spec
            .cms2_normal_volatility_sigma2
            .max_val(T::from_f64(1e-5));

        let conv1 = q1 * sig1 * sig1 * t0;
        let conv2 = q2 * sig2 * sig2 * t0;
        let mu1 = spec.cms1_forward_swap_rate_s1 + conv1;
        let mu2 = spec.cms2_forward_swap_rate_s2 + conv2;

        let spread_wedge = spec.leg1_gearing_w1 * conv1 - spec.leg2_gearing_w2 * conv2;
        let mu_spread = spec.leg1_gearing_w1 * mu1 - spec.leg2_gearing_w2 * mu2;

        let rho = spec
            .swap_rate_correlation_rho12
            .max_val(T::from_f64(-0.99))
            .min_val(T::from_f64(0.99));
        let w1_s1 = spec.leg1_gearing_w1 * sig1;
        let w2_s2 = spec.leg2_gearing_w2 * sig2;
        let var_spread = (w1_s1 * w1_s1 + w2_s2 * w2_s2 - T::from_f64(2.0) * rho * w1_s1 * w2_s2)
            .max_val(T::from_f64(1e-10))
            * t0;
        let std_spread = var_spread.sqrt();

        let d = (mu_spread - spec.spread_strike_k) / std_spread;
        let df = (-spec.ois_discount_rate_r * t0).exp();
        let pv = df * ((mu_spread - spec.spread_strike_k) * norm_cdf(d) + std_spread * norm_pdf(d));

        AndersenPiterbargCmsSpreadAnalytics {
            convexity_adjusted_cms1_forward: mu1,
            convexity_adjusted_cms2_forward: mu2,
            cms_spread_convexity_wedge: spread_wedge,
            cms_spread_call_present_value: pv,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn andersen_piterbarg_lsrm_cms_spread_invariants_and_ad_parity() {
            let base = AndersenPiterbargCmsSpreadSpec {
                cms1_forward_swap_rate_s1: 0.042_f64,
                cms1_tenor_years_tau1: 10.0,
                cms1_normal_volatility_sigma1: 0.0085,
                cms2_forward_swap_rate_s2: 0.034,
                cms2_tenor_years_tau2: 2.0,
                cms2_normal_volatility_sigma2: 0.0095,
                swap_rate_correlation_rho12: 0.78,
                leg1_gearing_w1: 1.0,
                leg2_gearing_w2: 1.0,
                spread_strike_k: 0.005,
                mean_reversion_kappa: 0.08,
                ois_discount_rate_r: 0.035,
                fixing_horizon_years_t0: 3.0,
            };
            let r0 = evaluate_andersen_piterbarg_lsrm_cms_spread_option(&base);
            assert!(r0.convexity_adjusted_cms1_forward > base.cms1_forward_swap_rate_s1);
            assert!(r0.convexity_adjusted_cms2_forward > base.cms2_forward_swap_rate_s2);
            assert!(r0.cms_spread_convexity_wedge > 0.0);
            assert!(r0.cms_spread_call_present_value > 0.0);

            let inputs = [
                base.cms1_forward_swap_rate_s1,
                base.cms1_tenor_years_tau1,
                base.cms1_normal_volatility_sigma1,
                base.cms2_forward_swap_rate_s2,
                base.cms2_tenor_years_tau2,
                base.cms2_normal_volatility_sigma2,
                base.swap_rate_correlation_rho12,
                base.leg1_gearing_w1,
                base.leg2_gearing_w2,
                base.spread_strike_k,
                base.mean_reversion_kappa,
                base.ois_discount_rate_r,
                base.fixing_horizon_years_t0,
            ];
            let mut dual_in = [Dual::<13>::from_f64(0.0); 13];
            for i in 0..13 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = AndersenPiterbargCmsSpreadSpec {
                cms1_forward_swap_rate_s1: dual_in[0],
                cms1_tenor_years_tau1: dual_in[1],
                cms1_normal_volatility_sigma1: dual_in[2],
                cms2_forward_swap_rate_s2: dual_in[3],
                cms2_tenor_years_tau2: dual_in[4],
                cms2_normal_volatility_sigma2: dual_in[5],
                swap_rate_correlation_rho12: dual_in[6],
                leg1_gearing_w1: dual_in[7],
                leg2_gearing_w2: dual_in[8],
                spread_strike_k: dual_in[9],
                mean_reversion_kappa: dual_in[10],
                ois_discount_rate_r: dual_in[11],
                fixing_horizon_years_t0: dual_in[12],
            };
            let r_ad = evaluate_andersen_piterbarg_lsrm_cms_spread_option(&spec_ad);
            let eps = 1e-6;
            for i in 0..13 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 13]| AndersenPiterbargCmsSpreadSpec {
                    cms1_forward_swap_rate_s1: v[0],
                    cms1_tenor_years_tau1: v[1],
                    cms1_normal_volatility_sigma1: v[2],
                    cms2_forward_swap_rate_s2: v[3],
                    cms2_tenor_years_tau2: v[4],
                    cms2_normal_volatility_sigma2: v[5],
                    swap_rate_correlation_rho12: v[6],
                    leg1_gearing_w1: v[7],
                    leg2_gearing_w2: v[8],
                    spread_strike_k: v[9],
                    mean_reversion_kappa: v[10],
                    ois_discount_rate_r: v[11],
                    fixing_horizon_years_t0: v[12],
                };
                let rp = evaluate_andersen_piterbarg_lsrm_cms_spread_option(&mk(p));
                let rm = evaluate_andersen_piterbarg_lsrm_cms_spread_option(&mk(m));
                let fd = (rp.cms_spread_call_present_value - rm.cms_spread_call_present_value)
                    / (2.0 * eps);
                assert!(
                    r_ad.cms_spread_call_present_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.cms_spread_call_present_value.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use andersen_piterbarg_lsrm_cms_spread_mod::*;

// ============================================================================
// Mercurio (2010, 2013) & Henrard (2014) Post-Crisis Multi-Curve FVA-Adjusted Swap Rate
// ============================================================================

pub mod mercurio_post_crisis_fva_multi_curve_swap_mod {
    use ql_ad::Real;

    /// Specification for the Mercurio (2010, 2013) / Henrard (2014) post-crisis multi-curve
    /// swap valuation with multiplicative IBOR-OIS basis and uncollateralized FVA spread tilt.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct MercurioFvaMultiCurveSwapSpec<T: Real> {
        /// Initial instantaneous OIS short rate $r_0^{\text{OIS}}$.
        pub ois_short_rate_r0: T,
        /// OIS term-structure slope parameter $\theta_{\text{OIS}}$.
        pub ois_curve_slope_theta: T,
        /// Initial multiplicative IBOR-OIS basis spread $b_0 > 0$.
        pub initial_multiplicative_ibor_basis_b0: T,
        /// Exponential decay speed of the IBOR-OIS term basis $\kappa_b > 0$.
        pub basis_term_decay_kappa_b: T,
        /// Unsecured bank funding borrowing spread $s_f > 0$ over OIS.
        pub funding_borrowing_spread_s_f: T,
        /// Bachelier swap-rate volatility $\sigma_s > 0$ driving uncollateralized EPE.
        pub swap_rate_volatility_sigma_s: T,
        /// Swap maturity $T > 0$ (years).
        pub swap_maturity_years_t: T,
    }

    /// Analytical outputs for the Mercurio-Henrard multi-curve and FVA-adjusted swap engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct MercurioFvaMultiCurveSwapAnalytics<T: Real> {
        /// OIS-discounted multi-curve swap annuity $A_{\text{OIS}}(0, T)$.
        pub ois_swap_annuity_pv01: T,
        /// Collateralized (CSA) multi-curve par swap rate $S_{\text{CSA}}(0, T)$.
        pub collateralized_multi_curve_swap_rate: T,
        /// Asymmetric Funding Valuation Adjustment charge $\text{FVA}$ per unit notional.
        pub uncollateralized_fva_charge: T,
        /// FVA-adjusted uncollateralized breakeven swap rate $S_{\text{FVA}}(0, T)$.
        pub fva_adjusted_uncollateralized_swap_rate: T,
    }

    /// Evaluates the Mercurio (2010, *Risk* 23(3); 2013) and Henrard (2014, Ch. 2 & 5)
    /// post-crisis multi-curve forward-IBOR vs. OIS swap rate and FVA spread adjustment.
    pub fn evaluate_mercurio_post_crisis_fva_adjusted_multi_curve_swap_rate<T: Real>(
        spec: &MercurioFvaMultiCurveSwapSpec<T>,
    ) -> MercurioFvaMultiCurveSwapAnalytics<T> {
        const PERIODS: usize = 4;
        let t_mat = spec.swap_maturity_years_t.max_val(T::from_f64(1e-4));
        let tau = t_mat / T::from_f64(PERIODS as f64);
        let kappa_b = spec.basis_term_decay_kappa_b.max_val(T::from_f64(1e-4));
        let sig_s = spec.swap_rate_volatility_sigma_s.max_val(T::from_f64(1e-5));
        let sf = spec.funding_borrowing_spread_s_f.max_val(T::from_f64(1e-5));

        let mut p_ois = [T::one(); PERIODS + 1];
        for (m, p) in p_ois.iter_mut().enumerate().skip(1) {
            let tm = T::from_f64(m as f64) * tau;
            let slope_factor = (T::one() - (-T::from_f64(0.5) * tm).exp()) / T::from_f64(0.5);
            let y_ois = spec.ois_short_rate_r0 + spec.ois_curve_slope_theta * slope_factor;
            *p = (-y_ois * tm).exp();
        }

        let mut annuity = T::zero();
        let mut float_leg_pv = T::zero();
        let mut l_fwd = [T::zero(); PERIODS];

        for m in 1..=PERIODS {
            let tm = T::from_f64(m as f64) * tau;
            let f_ois_m = (p_ois[m - 1] / p_ois[m] - T::one()) / tau;
            let b_m = spec.initial_multiplicative_ibor_basis_b0 * (-kappa_b * tm).exp();
            // Henrard (2014, Eq. 2.11) multiplicative basis: 1 + tau * L_m = (1 + tau * b_m) * (1 + tau * F_m^OIS)
            let l_m = f_ois_m + b_m * (T::one() + tau * f_ois_m);
            l_fwd[m - 1] = l_m;
            annuity += tau * p_ois[m];
            float_leg_pv += tau * p_ois[m] * l_m;
        }

        let s_csa = float_leg_pv / annuity;
        let inv_sqrt_2pi = T::from_f64(0.398_942_280_401_432_7);

        let mut fva = T::zero();
        for m in 1..=PERIODS {
            let t_exp = (T::from_f64(m as f64) - T::from_f64(0.5)) * tau;
            let mut rem_ann = T::zero();
            for p in p_ois.iter().skip(m) {
                rem_ann += tau * *p;
            }
            let epe_disc = rem_ann
                * (inv_sqrt_2pi * sig_s * t_exp.sqrt()
                    + T::from_f64(0.25) * tau * (l_fwd[m - 1] - s_csa));
            fva += tau * sf * epe_disc;
        }

        let s_fva = s_csa + fva / annuity;

        MercurioFvaMultiCurveSwapAnalytics {
            ois_swap_annuity_pv01: annuity,
            collateralized_multi_curve_swap_rate: s_csa,
            uncollateralized_fva_charge: fva,
            fva_adjusted_uncollateralized_swap_rate: s_fva,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn mercurio_fva_multi_curve_swap_invariants_and_ad_parity() {
            let base = MercurioFvaMultiCurveSwapSpec {
                ois_short_rate_r0: 0.038_f64,
                ois_curve_slope_theta: 0.004,
                initial_multiplicative_ibor_basis_b0: 0.0025,
                basis_term_decay_kappa_b: 0.25,
                funding_borrowing_spread_s_f: 0.0065,
                swap_rate_volatility_sigma_s: 0.0085,
                swap_maturity_years_t: 4.0,
            };
            let r0 = evaluate_mercurio_post_crisis_fva_adjusted_multi_curve_swap_rate(&base);
            assert!(r0.ois_swap_annuity_pv01 > 0.0);
            assert!(r0.collateralized_multi_curve_swap_rate > base.ois_short_rate_r0);
            assert!(r0.uncollateralized_fva_charge > 0.0);
            assert!(
                r0.fva_adjusted_uncollateralized_swap_rate
                    > r0.collateralized_multi_curve_swap_rate
            );

            let inputs = [
                base.ois_short_rate_r0,
                base.ois_curve_slope_theta,
                base.initial_multiplicative_ibor_basis_b0,
                base.basis_term_decay_kappa_b,
                base.funding_borrowing_spread_s_f,
                base.swap_rate_volatility_sigma_s,
                base.swap_maturity_years_t,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = MercurioFvaMultiCurveSwapSpec {
                ois_short_rate_r0: dual_in[0],
                ois_curve_slope_theta: dual_in[1],
                initial_multiplicative_ibor_basis_b0: dual_in[2],
                basis_term_decay_kappa_b: dual_in[3],
                funding_borrowing_spread_s_f: dual_in[4],
                swap_rate_volatility_sigma_s: dual_in[5],
                swap_maturity_years_t: dual_in[6],
            };
            let r_ad = evaluate_mercurio_post_crisis_fva_adjusted_multi_curve_swap_rate(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| MercurioFvaMultiCurveSwapSpec {
                    ois_short_rate_r0: v[0],
                    ois_curve_slope_theta: v[1],
                    initial_multiplicative_ibor_basis_b0: v[2],
                    basis_term_decay_kappa_b: v[3],
                    funding_borrowing_spread_s_f: v[4],
                    swap_rate_volatility_sigma_s: v[5],
                    swap_maturity_years_t: v[6],
                };
                let rp = evaluate_mercurio_post_crisis_fva_adjusted_multi_curve_swap_rate(&mk(p));
                let rm = evaluate_mercurio_post_crisis_fva_adjusted_multi_curve_swap_rate(&mk(m));
                let fd = (rp.fva_adjusted_uncollateralized_swap_rate
                    - rm.fva_adjusted_uncollateralized_swap_rate)
                    / (2.0 * eps);
                assert!(
                    r_ad.fva_adjusted_uncollateralized_swap_rate.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!(
                    (r_ad.fva_adjusted_uncollateralized_swap_rate.d[i] - fd).abs() < 1e-5,
                    "mismatch at {i}"
                );
            }
        }
    }
}
pub use mercurio_post_crisis_fva_multi_curve_swap_mod::*;
mod smith_wilson_ufr_extrapolation_mod {
    use ql_ad::Real;

    /// Input specification for the Smith & Wilson (2001) / EIOPA (2015) Ultimate Forward Rate (UFR)
    /// yield curve extrapolation model beyond the Last Liquid Point (LLP).
    #[derive(Debug, Clone, Copy)]
    pub struct SmithWilsonUfrExtrapolationSpec<T: Real> {
        /// Intermediate liquid pillar maturity $u_1$ in years (e.g., $10.0$).
        pub intermediate_pillar_years_u1: T,
        /// Last Liquid Point (LLP) maturity $u_2 = T_{\text{LLP}}$ in years (e.g., $20.0$).
        pub last_liquid_point_years_u2: T,
        /// Extrapolated ultra-long maturity $T_{\text{ext}} > u_2$ in years (e.g., $50.0$).
        pub extrapolated_maturity_years_t_ext: T,
        /// Continuously compounded zero rate at intermediate pillar $u_1$ (e.g., $0.028$).
        pub intermediate_zero_rate_r1: T,
        /// Continuously compounded zero rate at the Last Liquid Point $u_2$ (e.g., $0.031$).
        pub last_liquid_point_zero_rate_r2: T,
        /// Continuously compounded Ultimate Forward Rate $\omega_{\text{UFR}} = \ln(1 + \text{UFR})$ (e.g., $0.0345$).
        pub ultimate_forward_rate_omega: T,
        /// Smith–Wilson convergence speed parameter $\alpha > 0$ (e.g., $0.12$).
        pub convergence_speed_alpha: T,
    }

    /// Analytics returned by [`evaluate_smith_wilson_eiopa_ultimate_forward_rate_extrapolation`].
    #[derive(Debug, Clone, Copy)]
    pub struct SmithWilsonUfrExtrapolationAnalytics<T: Real> {
        /// Extrapolated zero-coupon discount factor $P_{\text{SW}}(0, T_{\text{ext}})$.
        pub extrapolated_discount_factor: T,
        /// Extrapolated continuously compounded zero yield $y_{\text{SW}}(T_{\text{ext}})$.
        pub extrapolated_zero_yield: T,
        /// Extrapolated instantaneous forward rate $f_{\text{SW}}(T_{\text{ext}})$.
        pub extrapolated_instantaneous_forward_rate: T,
        /// Signed UFR forward convergence gap $\omega_{\text{UFR}} - f_{\text{SW}}(T_{\text{ext}})$.
        pub ufr_forward_convergence_gap: T,
    }

    /// Evaluates the Smith & Wilson (2001) / EIOPA Solvency II Ultimate Forward Rate (UFR)
    /// yield curve extrapolation beyond the Last Liquid Point (LLP) using the exact symmetric
    /// Wilson kernel $W(t, u)$ calibrated to two liquid zero-coupon pillars $(u_1, u_2)$.
    pub fn evaluate_smith_wilson_eiopa_ultimate_forward_rate_extrapolation<T: Real>(
        spec: &SmithWilsonUfrExtrapolationSpec<T>,
    ) -> SmithWilsonUfrExtrapolationAnalytics<T> {
        let u1 = spec.intermediate_pillar_years_u1;
        let u2 = spec.last_liquid_point_years_u2;
        let t_ext = spec.extrapolated_maturity_years_t_ext;
        let r1 = spec.intermediate_zero_rate_r1;
        let r2 = spec.last_liquid_point_zero_rate_r2;
        let omega = spec.ultimate_forward_rate_omega;
        let alpha = spec.convergence_speed_alpha;
        let half = T::from_f64(0.5);

        // Wilson kernel H(min(t,u), max(t,u)) = alpha * t_min - 0.5 * exp(-alpha * t_max) * (exp(alpha * t_min) - exp(-alpha * t_min))
        // Ordered pillars: 0 < u1 < u2 < t_ext
        let sinh2_u1 = (alpha * u1).exp() - (-alpha * u1).exp();
        let sinh2_u2 = (alpha * u2).exp() - (-alpha * u2).exp();

        let h11 = alpha * u1 - half * (-alpha * u1).exp() * sinh2_u1;
        let h12 = alpha * u1 - half * (-alpha * u2).exp() * sinh2_u1;
        let h22 = alpha * u2 - half * (-alpha * u2).exp() * sinh2_u2;

        let w11 = (-omega * (u1 + u1)).exp() * h11;
        let w12 = (-omega * (u1 + u2)).exp() * h12;
        let w22 = (-omega * (u2 + u2)).exp() * h22;

        // Observed liquid discount factors minus UFR baseline discount factors
        let d1 = (-r1 * u1).exp() - (-omega * u1).exp();
        let d2 = (-r2 * u2).exp() - (-omega * u2).exp();

        // Exact 2x2 symmetric kernel linear system W * zeta = d
        let det = w11 * w22 - w12 * w12;
        let zeta1 = (w22 * d1 - w12 * d2) / det;
        let zeta2 = (w11 * d2 - w12 * d1) / det;

        // Evaluate kernel and time derivative at t_ext > u2 > u1
        let exp_neg_alpha_t = (-alpha * t_ext).exp();
        let h_t1 = alpha * u1 - half * exp_neg_alpha_t * sinh2_u1;
        let h_t2 = alpha * u2 - half * exp_neg_alpha_t * sinh2_u2;

        let dh_dt1 = half * alpha * exp_neg_alpha_t * sinh2_u1;
        let dh_dt2 = half * alpha * exp_neg_alpha_t * sinh2_u2;

        let e_t1 = (-omega * (t_ext + u1)).exp();
        let e_t2 = (-omega * (t_ext + u2)).exp();

        let w_t1 = e_t1 * h_t1;
        let w_t2 = e_t2 * h_t2;

        let dw_dt1 = e_t1 * (-omega * h_t1 + dh_dt1);
        let dw_dt2 = e_t2 * (-omega * h_t2 + dh_dt2);

        let p_ufr_t = (-omega * t_ext).exp();
        let p_sw = p_ufr_t + zeta1 * w_t1 + zeta2 * w_t2;
        let dp_sw_dt = -omega * p_ufr_t + zeta1 * dw_dt1 + zeta2 * dw_dt2;

        let y_sw = -p_sw.ln() / t_ext;
        let f_sw = -dp_sw_dt / p_sw;

        SmithWilsonUfrExtrapolationAnalytics {
            extrapolated_discount_factor: p_sw,
            extrapolated_zero_yield: y_sw,
            extrapolated_instantaneous_forward_rate: f_sw,
            ufr_forward_convergence_gap: omega - f_sw,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn smith_wilson_ufr_extrapolation_invariants_and_ad_parity() {
            let base = SmithWilsonUfrExtrapolationSpec {
                intermediate_pillar_years_u1: 10.0_f64,
                last_liquid_point_years_u2: 20.0_f64,
                extrapolated_maturity_years_t_ext: 50.0_f64,
                intermediate_zero_rate_r1: 0.026_f64,
                last_liquid_point_zero_rate_r2: 0.029_f64,
                ultimate_forward_rate_omega: 0.036_f64,
                convergence_speed_alpha: 0.12_f64,
            };
            let res = evaluate_smith_wilson_eiopa_ultimate_forward_rate_extrapolation(&base);
            assert!(
                res.extrapolated_discount_factor > 0.0 && res.extrapolated_discount_factor < 1.0
            );
            assert!(res.extrapolated_zero_yield > 0.029 && res.extrapolated_zero_yield < 0.036);
            assert!(res.ufr_forward_convergence_gap.abs() < 0.005);

            let inputs = [
                base.intermediate_pillar_years_u1,
                base.last_liquid_point_years_u2,
                base.extrapolated_maturity_years_t_ext,
                base.intermediate_zero_rate_r1,
                base.last_liquid_point_zero_rate_r2,
                base.ultimate_forward_rate_omega,
                base.convergence_speed_alpha,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = SmithWilsonUfrExtrapolationSpec {
                intermediate_pillar_years_u1: dual_in[0],
                last_liquid_point_years_u2: dual_in[1],
                extrapolated_maturity_years_t_ext: dual_in[2],
                intermediate_zero_rate_r1: dual_in[3],
                last_liquid_point_zero_rate_r2: dual_in[4],
                ultimate_forward_rate_omega: dual_in[5],
                convergence_speed_alpha: dual_in[6],
            };
            let r_ad = evaluate_smith_wilson_eiopa_ultimate_forward_rate_extrapolation(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| SmithWilsonUfrExtrapolationSpec {
                    intermediate_pillar_years_u1: v[0],
                    last_liquid_point_years_u2: v[1],
                    extrapolated_maturity_years_t_ext: v[2],
                    intermediate_zero_rate_r1: v[3],
                    last_liquid_point_zero_rate_r2: v[4],
                    ultimate_forward_rate_omega: v[5],
                    convergence_speed_alpha: v[6],
                };
                let rp = evaluate_smith_wilson_eiopa_ultimate_forward_rate_extrapolation(&mk(p));
                let rm = evaluate_smith_wilson_eiopa_ultimate_forward_rate_extrapolation(&mk(m));
                let fd = (rp.extrapolated_zero_yield - rm.extrapolated_zero_yield) / (2.0 * eps);
                assert!(
                    r_ad.extrapolated_zero_yield.d[i].abs() > 1e-7,
                    "zero derivative at {i}"
                );
                assert!((r_ad.extrapolated_zero_yield.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use smith_wilson_ufr_extrapolation_mod::*;
mod hagan_west_monotone_convex_mod {
    use ql_ad::Real;

    /// Input specification for the Hagan & West (2006, *Applied Mathematical Finance* 13(2);
    /// 2008, *Wilmott Magazine*) monotone convex forward-rate curve interpolation engine.
    #[derive(Debug, Clone, Copy)]
    pub struct HaganWestMonotoneConvexSpec<T: Real> {
        /// Left interval start maturity $t_{i-1}$ in years (e.g., $2.0$).
        pub interval_start_years_t0: T,
        /// Right interval end maturity $t_i$ in years (e.g., $5.0$).
        pub interval_end_years_t1: T,
        /// Target evaluation maturity $t \in (t_{i-1}, t_i)$ in years (e.g., $3.5$).
        pub target_maturity_years_t: T,
        /// Continuously compounded zero rate $r(t_{i-1})$ at $t_{i-1}$ (e.g., $0.031$).
        pub start_zero_rate_r0: T,
        /// Discrete forward rate $f_{i-1}^d$ on the preceding interval $[t_{i-2}, t_{i-1}]$ (e.g., $0.032$).
        pub previous_discrete_forward_f_prev: T,
        /// Discrete forward rate $f_i^d$ on the active interval $[t_{i-1}, t_i]$ (e.g., $0.038$).
        pub current_discrete_forward_f_curr: T,
        /// Discrete forward rate $f_{i+1}^d$ on the succeeding interval $[t_i, t_{i+1}]$ (e.g., $0.043$).
        pub next_discrete_forward_f_next: T,
        /// Adjacent-interval blending weight $w \in (0, 1)$ at interior pillar nodes (e.g., $0.50$).
        pub node_blending_weight_w: T,
    }

    /// Analytics returned by [`evaluate_hagan_west_monotone_convex_forward_interpolation`].
    #[derive(Debug, Clone, Copy)]
    pub struct HaganWestMonotoneConvexAnalytics<T: Real> {
        /// Left pillar instantaneous forward rate $f(t_{i-1})$.
        pub left_node_instantaneous_forward: T,
        /// Right pillar instantaneous forward rate $f(t_i)$.
        pub right_node_instantaneous_forward: T,
        /// Monotone-convex instantaneous forward rate $f(t)$ at $t \in (t_{i-1}, t_i)$.
        pub instantaneous_forward_rate: T,
        /// Monotone-convex interpolated zero-coupon yield $r(t)$ at $t$.
        pub interpolated_zero_rate: T,
    }

    /// Evaluates the Hagan & West (2006, Eqs. 16–25) monotone convex quadratic forward spline
    /// on normalized coordinate $x = \frac{t - t_{i-1}}{t_i - t_{i-1}} \in (0, 1)$:
    /// $$f(x) = f_i^d + g_0(1 - 4x + 3x^2) + g_1(-2x + 3x^2),$$
    /// where $g_0 = f(t_{i-1}) - f_i^d$ and $g_1 = f(t_i) - f_i^d$ satisfy $\int_0^1 (f(x) - f_i^d)\,dx = 0$,
    /// and recovers the exact zero rate $r(t) = \frac{t_{i-1} r(t_{i-1}) + (t_i - t_{i-1})\int_0^x f(u)\,du}{t}$.
    pub fn evaluate_hagan_west_monotone_convex_forward_interpolation<T: Real>(
        spec: &HaganWestMonotoneConvexSpec<T>,
    ) -> HaganWestMonotoneConvexAnalytics<T> {
        let t0 = spec.interval_start_years_t0;
        let t1 = spec.interval_end_years_t1;
        let t = spec.target_maturity_years_t;
        let dt = t1 - t0;
        let x = (t - t0) / dt;

        let w = spec.node_blending_weight_w;
        let f_prev = spec.previous_discrete_forward_f_prev;
        let f_curr = spec.current_discrete_forward_f_curr;
        let f_next = spec.next_discrete_forward_f_next;

        // Hagan & West (2006, Eq. 16) node instantaneous forwards
        let f_left = w * f_prev + (T::one() - w) * f_curr;
        let f_right = w * f_curr + (T::one() - w) * f_next;

        let g0 = f_left - f_curr;
        let g1 = f_right - f_curr;
        let x2 = x * x;
        let x3 = x2 * x;

        // Quadratic zero-integral basis g(x) = g0*(1 - 4x + 3x^2) + g1*(-2x + 3x^2)
        let g_x = g0 * (T::one() - T::from_f64(4.0) * x + T::from_f64(3.0) * x2)
            + g1 * (-T::from_f64(2.0) * x + T::from_f64(3.0) * x2);
        let f_inst = f_curr + g_x;

        // Exact integral int_0^x f(u) du = f_curr*x + g0*(x - 2x^2 + x^3) + g1*(-x^2 + x^3)
        let int_f = f_curr * x + g0 * (x - T::from_f64(2.0) * x2 + x3) + g1 * (-x2 + x3);

        let r_interp = (t0 * spec.start_zero_rate_r0 + dt * int_f) / t;

        HaganWestMonotoneConvexAnalytics {
            left_node_instantaneous_forward: f_left,
            right_node_instantaneous_forward: f_right,
            instantaneous_forward_rate: f_inst,
            interpolated_zero_rate: r_interp,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hagan_west_monotone_convex_invariants_and_ad_parity() {
            let base = HaganWestMonotoneConvexSpec {
                interval_start_years_t0: 2.0_f64,
                interval_end_years_t1: 5.0_f64,
                target_maturity_years_t: 3.4_f64,
                start_zero_rate_r0: 0.031_f64,
                previous_discrete_forward_f_prev: 0.032_f64,
                current_discrete_forward_f_curr: 0.038_f64,
                next_discrete_forward_f_next: 0.044_f64,
                node_blending_weight_w: 0.50_f64,
            };
            let res = evaluate_hagan_west_monotone_convex_forward_interpolation(&base);
            assert!(res.left_node_instantaneous_forward < res.right_node_instantaneous_forward);
            assert!(
                res.instantaneous_forward_rate >= res.left_node_instantaneous_forward
                    && res.instantaneous_forward_rate <= res.right_node_instantaneous_forward
            );
            assert!(res.interpolated_zero_rate > base.start_zero_rate_r0);

            let inputs = [
                base.interval_start_years_t0,
                base.interval_end_years_t1,
                base.target_maturity_years_t,
                base.start_zero_rate_r0,
                base.previous_discrete_forward_f_prev,
                base.current_discrete_forward_f_curr,
                base.next_discrete_forward_f_next,
                base.node_blending_weight_w,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HaganWestMonotoneConvexSpec {
                interval_start_years_t0: dual_in[0],
                interval_end_years_t1: dual_in[1],
                target_maturity_years_t: dual_in[2],
                start_zero_rate_r0: dual_in[3],
                previous_discrete_forward_f_prev: dual_in[4],
                current_discrete_forward_f_curr: dual_in[5],
                next_discrete_forward_f_next: dual_in[6],
                node_blending_weight_w: dual_in[7],
            };
            let r_ad = evaluate_hagan_west_monotone_convex_forward_interpolation(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| HaganWestMonotoneConvexSpec {
                    interval_start_years_t0: v[0],
                    interval_end_years_t1: v[1],
                    target_maturity_years_t: v[2],
                    start_zero_rate_r0: v[3],
                    previous_discrete_forward_f_prev: v[4],
                    current_discrete_forward_f_curr: v[5],
                    next_discrete_forward_f_next: v[6],
                    node_blending_weight_w: v[7],
                };
                let rp = evaluate_hagan_west_monotone_convex_forward_interpolation(&mk(p));
                let rm = evaluate_hagan_west_monotone_convex_forward_interpolation(&mk(m));
                let fd = (rp.interpolated_zero_rate - rm.interpolated_zero_rate) / (2.0 * eps);
                assert!(
                    r_ad.interpolated_zero_rate.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.interpolated_zero_rate.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use hagan_west_monotone_convex_mod::*;
mod filipovic_larsson_trolle_linear_rational_mod {
    use ql_ad::Real;

    /// Input specification for the Filipović, Larsson & Trolle (2017, *Journal of Finance*
    /// 72(2), *"Linear-Rational Term Structure Models"*) two-factor Zero-Lower-Bound (ZLB)
    /// discount curve and yield volatility compression engine.
    #[derive(Debug, Clone, Copy)]
    pub struct FilipovicLarssonTrolleLinearRationalCurveSpec<T: Real> {
        /// Asymptotic long-horizon discount rate $\alpha > 0$ (e.g., $0.042$).
        pub asymptotic_rate_alpha: T,
        /// Current level of positive state factor 1 $X_{1,0} > 0$ (e.g., $0.35$).
        pub factor_1_state_x1: T,
        /// Current level of positive state factor 2 $X_{2,0} > 0$ (e.g., $0.25$).
        pub factor_2_state_x2: T,
        /// Long-run mean $\theta_1 > 0$ of state factor 1 under $\mathbb{P}$ (e.g., $0.30$).
        pub factor_1_mean_theta1: T,
        /// Long-run mean $\theta_2 > 0$ of state factor 2 under $\mathbb{P}$ (e.g., $0.20$).
        pub factor_2_mean_theta2: T,
        /// Mean-reversion speed $\kappa_1 > 0$ of state factor 1 (e.g., $0.45$).
        pub factor_1_speed_kappa1: T,
        /// Mean-reversion speed $\kappa_2 > 0$ of state factor 2 (e.g., $0.18$).
        pub factor_2_speed_kappa2: T,
        /// Diffusion volatility scale $\sigma_1 > 0$ of state factor 1 (e.g., $0.14$).
        pub factor_1_volatility_sigma1: T,
        /// Maturity $T$ in years (e.g., $5.0$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_filipovic_larsson_trolle_linear_rational_zlb_curve`].
    #[derive(Debug, Clone, Copy)]
    pub struct FilipovicLarssonTrolleLinearRationalCurveAnalytics<T: Real> {
        /// Instantaneous ZLB-consistent short rate $r_0 \ge 0$.
        pub instantaneous_short_rate_r0: T,
        /// Linear-rational zero-coupon bond price $P(0, T) \in (0, 1)$.
        pub linear_rational_discount_factor: T,
        /// Continuously compounded zero-coupon yield $y(0, T) = -\frac{1}{T}\ln P(0, T)$.
        pub linear_rational_zero_yield: T,
        /// Endogenous state-dependent Black yield volatility $\sigma_{\text{Black}}(0, T) = \sigma_y(0, T) / y(0, T)$.
        pub state_dependent_yield_volatility: T,
    }

    /// Evaluates the Filipović, Larsson & Trolle (2017, Eqs. 4–12) linear-rational state-price
    /// deflator $\zeta_t = e^{-\alpha t}(1 + X_{1,t} + X_{2,t})$, yielding closed-form bond price:
    /// $$P(0, T) = e^{-\alpha T} \frac{1 + \mathbb{E}_0^\mathbb{P}[X_{1,T}] + \mathbb{E}_0^\mathbb{P}[X_{2,T}]}{1 + X_{1,0} + X_{2,0}}, \qquad r_0 = \alpha - \frac{\kappa_1(\theta_1 - X_{1,0}) + \kappa_2(\theta_2 - X_{2,0})}{1 + X_{1,0} + X_{2,0}},$$
    /// and state-dependent relative yield volatility $\sigma_{\text{Black}}(0, T) = \frac{\sqrt{(e^{-\kappa_1 T}\sigma_1)^2 X_{1,0} + (0.5 e^{-\kappa_2 T}\sigma_1)^2 X_{2,0}}}{T \, y(0, T)\big(1 + \mathbb{E}_0^\mathbb{P}[X_{1,T}] + \mathbb{E}_0^\mathbb{P}[X_{2,T}]\big)}$.
    pub fn evaluate_filipovic_larsson_trolle_linear_rational_zlb_curve<T: Real>(
        spec: &FilipovicLarssonTrolleLinearRationalCurveSpec<T>,
    ) -> FilipovicLarssonTrolleLinearRationalCurveAnalytics<T> {
        let alpha = spec.asymptotic_rate_alpha;
        let x1 = spec.factor_1_state_x1;
        let x2 = spec.factor_2_state_x2;
        let th1 = spec.factor_1_mean_theta1;
        let th2 = spec.factor_2_mean_theta2;
        let k1 = spec.factor_1_speed_kappa1;
        let k2 = spec.factor_2_speed_kappa2;
        let sig1 = spec.factor_1_volatility_sigma1;
        let t = spec.maturity_years_t;

        let denom_0 = T::one() + x1 + x2;
        let drift_0 = k1 * (th1 - x1) + k2 * (th2 - x2);
        let r0 = alpha - drift_0 / denom_0;

        let e1 = (-k1 * t).exp();
        let e2 = (-k2 * t).exp();
        let ex1_t = th1 + (x1 - th1) * e1;
        let ex2_t = th2 + (x2 - th2) * e2;
        let num_t = T::one() + ex1_t + ex2_t;

        let df = (-alpha * t).exp() * (num_t / denom_0);
        let y_t = -df.ln() / t;

        let sig2 = T::from_f64(0.5) * sig1;
        let var_num = (e1 * sig1) * (e1 * sig1) * x1 + (e2 * sig2) * (e2 * sig2) * x2;
        let sig_y = var_num.sqrt() / (t * num_t * y_t);

        FilipovicLarssonTrolleLinearRationalCurveAnalytics {
            instantaneous_short_rate_r0: r0,
            linear_rational_discount_factor: df,
            linear_rational_zero_yield: y_t,
            state_dependent_yield_volatility: sig_y,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn filipovic_larsson_trolle_linear_rational_invariants_and_ad_parity() {
            let base = FilipovicLarssonTrolleLinearRationalCurveSpec {
                asymptotic_rate_alpha: 0.042_f64,
                factor_1_state_x1: 0.35_f64,
                factor_2_state_x2: 0.25_f64,
                factor_1_mean_theta1: 0.30_f64,
                factor_2_mean_theta2: 0.20_f64,
                factor_1_speed_kappa1: 0.45_f64,
                factor_2_speed_kappa2: 0.18_f64,
                factor_1_volatility_sigma1: 0.14_f64,
                maturity_years_t: 5.0_f64,
            };
            let res = evaluate_filipovic_larsson_trolle_linear_rational_zlb_curve(&base);
            assert!(res.instantaneous_short_rate_r0 > 0.0);
            assert!(
                res.linear_rational_discount_factor > 0.0
                    && res.linear_rational_discount_factor < 1.0
            );
            assert!(res.linear_rational_zero_yield > 0.0);
            assert!(res.state_dependent_yield_volatility > 0.0);

            let inputs = [
                base.asymptotic_rate_alpha,
                base.factor_1_state_x1,
                base.factor_2_state_x2,
                base.factor_1_mean_theta1,
                base.factor_2_mean_theta2,
                base.factor_1_speed_kappa1,
                base.factor_2_speed_kappa2,
                base.factor_1_volatility_sigma1,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<9>::from_value(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = FilipovicLarssonTrolleLinearRationalCurveSpec {
                asymptotic_rate_alpha: dual_in[0],
                factor_1_state_x1: dual_in[1],
                factor_2_state_x2: dual_in[2],
                factor_1_mean_theta1: dual_in[3],
                factor_2_mean_theta2: dual_in[4],
                factor_1_speed_kappa1: dual_in[5],
                factor_2_speed_kappa2: dual_in[6],
                factor_1_volatility_sigma1: dual_in[7],
                maturity_years_t: dual_in[8],
            };
            let r_ad = evaluate_filipovic_larsson_trolle_linear_rational_zlb_curve(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| FilipovicLarssonTrolleLinearRationalCurveSpec {
                    asymptotic_rate_alpha: v[0],
                    factor_1_state_x1: v[1],
                    factor_2_state_x2: v[2],
                    factor_1_mean_theta1: v[3],
                    factor_2_mean_theta2: v[4],
                    factor_1_speed_kappa1: v[5],
                    factor_2_speed_kappa2: v[6],
                    factor_1_volatility_sigma1: v[7],
                    maturity_years_t: v[8],
                };
                let rp = evaluate_filipovic_larsson_trolle_linear_rational_zlb_curve(&mk(p));
                let rm = evaluate_filipovic_larsson_trolle_linear_rational_zlb_curve(&mk(m));
                let fd = (rp.state_dependent_yield_volatility
                    - rm.state_dependent_yield_volatility)
                    / (2.0 * eps);
                assert!(
                    r_ad.state_dependent_yield_volatility.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.state_dependent_yield_volatility.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use filipovic_larsson_trolle_linear_rational_mod::*;
