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
