mod evaluate_gregory_csa_threshold_mta_and_ctd_colva_module {
    use ql_ad::Real;

    #[derive(Clone, Debug)]
    pub struct GregoryCsaThresholdMtaCtdColvaSpec<T: Real = f64> {
        pub time_grid: Vec<T>,
        pub expected_mtm_profile: Vec<T>,
        pub exposure_vol: T,
        pub mpor_years: T,
        pub csa_threshold: T,
        pub minimum_transfer_amount: T,
        pub independent_amount: T,
        pub hazard_rate: T,
        pub lgd: T,
        pub risk_free_rate: T,
        pub ccy1_collateral_spread: T,
        pub ccy2_collateral_spread: T,
        pub collateral_spread_vol: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct GregoryCsaThresholdMtaCtdColvaReport<T: Real = f64> {
        pub effective_csa_threshold: T,
        pub uncollateralized_cva: T,
        pub threshold_mta_collateralized_cva: T,
        pub csa_cva_reduction_benefit: T,
        pub single_currency_colva: T,
        pub ctd_multi_currency_colva: T,
        pub ctd_collateral_switch_option_value: T,
    }

    pub fn evaluate_gregory_csa_threshold_mta_and_ctd_colva<T: Real>(
        s: &GregoryCsaThresholdMtaCtdColvaSpec<T>,
    ) -> GregoryCsaThresholdMtaCtdColvaReport<T> {
        let threshold = (s.csa_threshold + s.minimum_transfer_amount / T::from_f64(2.0)
            - s.independent_amount)
            .max_val(T::zero());
        let sigma = (s.exposure_vol * s.mpor_years.sqrt()).max_val(T::from_f64(1e-10));
        let mut ucva = T::zero();
        let mut ccva = T::zero();
        let mut single = T::zero();
        let mut ctd = T::zero();
        for i in 0..s.time_grid.len().min(s.expected_mtm_profile.len()) {
            let t = s.time_grid[i];
            let mu = s.expected_mtm_profile[i];
            let z0 = mu / sigma;
            let zt = (mu - threshold) / sigma;
            let epe0 = mu * z0.norm_cdf() + sigma * z0.norm_pdf();
            let epet = (mu - threshold) * zt.norm_cdf() + sigma * zt.norm_pdf();
            let disc = (-s.risk_free_rate * t).exp();
            let dp = s.hazard_rate * (-s.hazard_rate * t).exp() * disc * s.lgd;
            ucva += dp * epe0;
            ccva += dp * epet;
            let ds = s.collateral_spread_vol * t.sqrt();
            let ds_safe = ds.max_val(T::from_f64(1e-10));
            let zeta = (s.ccy1_collateral_spread - s.ccy2_collateral_spread) / ds_safe;
            let max_spread = s.ccy1_collateral_spread * zeta.norm_cdf()
                + s.ccy2_collateral_spread * (-zeta).norm_cdf()
                + ds * zeta.norm_pdf();
            let dt = if i == 0 {
                t.max_val(T::from_f64(1e-8))
            } else {
                (t - s.time_grid[i - 1]).max_val(T::from_f64(1e-8))
            };
            let col_bal = (mu - threshold).max_val(T::zero());
            single += disc * dt * col_bal * s.ccy1_collateral_spread;
            ctd += disc * dt * col_bal * max_spread;
        }
        let option = (ctd - single).max_val(T::zero());
        GregoryCsaThresholdMtaCtdColvaReport {
            effective_csa_threshold: threshold,
            uncollateralized_cva: ucva,
            threshold_mta_collateralized_cva: ccva,
            csa_cva_reduction_benefit: ucva - ccva,
            single_currency_colva: single,
            ctd_multi_currency_colva: ctd,
            ctd_collateral_switch_option_value: option,
        }
    }

    #[cfg(test)]
    mod domain_suite_43_tests {
        use super::*;
        use ql_ad::Dual;

        fn spec<T: Real>(x: T) -> GregoryCsaThresholdMtaCtdColvaSpec<T> {
            GregoryCsaThresholdMtaCtdColvaSpec {
                time_grid: vec![T::one(), T::from_f64(2.0), T::from_f64(3.0)],
                expected_mtm_profile: vec![
                    T::from_f64(5.0) + x,
                    T::from_f64(8.0),
                    T::from_f64(10.0),
                ],
                exposure_vol: T::from_f64(3.0),
                mpor_years: T::from_f64(10.0 / 252.0),
                csa_threshold: T::from_f64(2.0),
                minimum_transfer_amount: T::one(),
                independent_amount: T::from_f64(0.5),
                hazard_rate: T::from_f64(0.02),
                lgd: T::from_f64(0.6),
                risk_free_rate: T::from_f64(0.03),
                ccy1_collateral_spread: T::from_f64(0.01),
                ccy2_collateral_spread: T::from_f64(0.008),
                collateral_spread_vol: T::from_f64(0.01),
            }
        }

        #[test]
        fn suite43_invariants_and_dual_fd() {
            let r = evaluate_gregory_csa_threshold_mta_and_ctd_colva(&spec(0.0));
            assert!(r.threshold_mta_collateralized_cva <= r.uncollateralized_cva);
            assert!(r.ctd_multi_currency_colva >= r.single_currency_colva);
            assert!(r.ctd_collateral_switch_option_value > 0.0);
            let h = 1e-6;
            let fd = (evaluate_gregory_csa_threshold_mta_and_ctd_colva(&spec(h))
                .uncollateralized_cva
                - evaluate_gregory_csa_threshold_mta_and_ctd_colva(&spec(-h)).uncollateralized_cva)
                / (2.0 * h);
            let ad = evaluate_gregory_csa_threshold_mta_and_ctd_colva(&spec(Dual::<1>::variable(
                0.0, 0,
            )))
            .uncollateralized_cva;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_gregory_csa_threshold_mta_and_ctd_colva_module::*;

mod evaluate_turlakov_green_kenyon_wwr_and_kva_module {
    use ql_ad::Real;
    #[derive(Clone, Debug)]
    pub struct TurlakovGreenKenyonWwrKvaSpec<T: Real = f64> {
        pub time_grid: Vec<T>,
        pub independent_epe_profile: Vec<T>,
        pub pfe_addon_profile: Vec<T>,
        pub base_cpty_hazard: T,
        pub hazard_log_vol: T,
        pub exposure_vol: T,
        pub wwr_correlation_rho: T,
        pub lgd: T,
        pub risk_free_rate: T,
        pub bank_hazard_rate: T,
        pub ccr_risk_weight: T,
        pub shareholder_cost_of_capital: T,
        pub capital_return_offset_beta: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct TurlakovGreenKenyonWwrKvaReport<T: Real = f64> {
        pub independent_cva: T,
        pub turlakov_wwr_cva: T,
        pub wwr_amplification_multiplier: T,
        pub ccr_sa_ead_t0: T,
        pub green_kenyon_kva: T,
        pub total_wwr_cva_plus_kva: T,
    }
    pub fn evaluate_turlakov_green_kenyon_wwr_and_kva<T: Real>(
        s: &TurlakovGreenKenyonWwrKvaSpec<T>,
    ) -> TurlakovGreenKenyonWwrKvaReport<T> {
        let mut ic = T::zero();
        let mut wc = T::zero();
        let mut kva = T::zero();
        for i in 0..s
            .time_grid
            .len()
            .min(s.independent_epe_profile.len())
            .min(s.pfe_addon_profile.len())
        {
            let t = s.time_grid[i];
            let dt = if i == 0 { t } else { t - s.time_grid[i - 1] };
            let df = (-(s.risk_free_rate + s.base_cpty_hazard + s.bank_hazard_rate) * t).exp();
            let m = (s.wwr_correlation_rho * s.hazard_log_vol * s.exposure_vol * t).exp();
            ic += dt * df * s.lgd * s.base_cpty_hazard * s.independent_epe_profile[i];
            wc += dt * df * s.lgd * s.base_cpty_hazard * s.independent_epe_profile[i] * m;
            let reg = s.ccr_risk_weight
                * T::from_f64(1.4)
                * (s.independent_epe_profile[i] + s.pfe_addon_profile[i]);
            kva += dt
                * df
                * (s.shareholder_cost_of_capital - s.capital_return_offset_beta * s.risk_free_rate)
                * reg;
        }
        TurlakovGreenKenyonWwrKvaReport {
            independent_cva: ic,
            turlakov_wwr_cva: wc,
            wwr_amplification_multiplier: wc / ic.max_val(T::from_f64(1e-8)),
            ccr_sa_ead_t0: T::from_f64(1.4)
                * (s.independent_epe_profile
                    .first()
                    .copied()
                    .unwrap_or(T::zero())
                    + s.pfe_addon_profile.first().copied().unwrap_or(T::zero())),
            green_kenyon_kva: kva,
            total_wwr_cva_plus_kva: wc + kva,
        }
    }
    #[cfg(test)]
    mod domain_suite_44_tests {
        use super::*;
        use ql_ad::Dual;

        fn spec<T: Real>(x: T) -> TurlakovGreenKenyonWwrKvaSpec<T> {
            TurlakovGreenKenyonWwrKvaSpec {
                time_grid: vec![T::one()],
                independent_epe_profile: vec![T::one() + x],
                pfe_addon_profile: vec![T::one()],
                base_cpty_hazard: T::from_f64(0.02),
                hazard_log_vol: T::from_f64(0.3),
                exposure_vol: T::from_f64(0.2),
                wwr_correlation_rho: T::from_f64(0.5),
                lgd: T::from_f64(0.6),
                risk_free_rate: T::from_f64(0.02),
                bank_hazard_rate: T::from_f64(0.01),
                ccr_risk_weight: T::from_f64(0.08),
                shareholder_cost_of_capital: T::from_f64(0.1),
                capital_return_offset_beta: T::from_f64(0.5),
            }
        }

        #[test]
        fn wwr_positive() {
            let r = evaluate_turlakov_green_kenyon_wwr_and_kva(&spec(0.0));
            assert!(r.turlakov_wwr_cva > r.independent_cva && r.green_kenyon_kva > 0.0);
            let h = 1e-6;
            let fd = (evaluate_turlakov_green_kenyon_wwr_and_kva(&spec(h)).total_wwr_cva_plus_kva
                - evaluate_turlakov_green_kenyon_wwr_and_kva(&spec(-h)).total_wwr_cva_plus_kva)
                / (2.0 * h);
            let ad = evaluate_turlakov_green_kenyon_wwr_and_kva(&spec(Dual::<1>::variable(0.0, 0)))
                .total_wwr_cva_plus_kva;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_turlakov_green_kenyon_wwr_and_kva_module::*;

mod evaluate_gregory_sa_cva_and_euler_netting_allocation_module {
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct GregorySaCvaEulerNettingSpec<T: Real = f64> {
        pub cva_delta_sensitivities: [T; 2],
        pub cva_hedge_sensitivities: [T; 2],
        pub delta_risk_weights: [T; 2],
        pub intra_bucket_corr_rho: T,
        pub cva_vega_sensitivity: T,
        pub vega_risk_weight: T,
        pub subportfolio_exposures: [T; 3],
        pub exposure_corr_3x3: [[T; 3]; 3],
        pub cpty_hazard_rate: T,
        pub cpty_lgd: T,
        pub horizon_years: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct GregorySaCvaEulerNettingReport<T: Real = f64> {
        pub sa_cva_delta_capital: T,
        pub sa_cva_vega_capital: T,
        pub total_sa_cva_capital: T,
        pub standalone_cva_sum: T,
        pub netted_portfolio_cva: T,
        pub netting_cva_reduction_benefit: T,
        pub euler_allocated_cva: [T; 3],
        pub netting_set_euler_allocation_error: T,
    }

    pub fn evaluate_gregory_sa_cva_and_euler_netting_allocation<T: Real>(
        s: &GregorySaCvaEulerNettingSpec<T>,
    ) -> GregorySaCvaEulerNettingReport<T> {
        let ws0 =
            s.delta_risk_weights[0] * (s.cva_delta_sensitivities[0] - s.cva_hedge_sensitivities[0]);
        let ws1 =
            s.delta_risk_weights[1] * (s.cva_delta_sensitivities[1] - s.cva_hedge_sensitivities[1]);
        let h0 = s.delta_risk_weights[0] * s.cva_hedge_sensitivities[0];
        let h1 = s.delta_risk_weights[1] * s.cva_hedge_sensitivities[1];
        let kd = T::from_f64(1.25)
            * (ws0 * ws0
                + ws1 * ws1
                + T::from_f64(2.0) * s.intra_bucket_corr_rho * ws0 * ws1
                + T::from_f64(0.01) * (h0 * h0 + h1 * h1))
                .sqrt();
        let v_prod = s.vega_risk_weight * s.cva_vega_sensitivity;
        let kv = T::from_f64(1.25) * (v_prod * v_prod + T::from_f64(1e-24)).sqrt();
        let mut quad = T::zero();
        for i in 0..3 {
            for j in 0..3 {
                quad += s.exposure_corr_3x3[i][j]
                    * s.subportfolio_exposures[i]
                    * s.subportfolio_exposures[j];
            }
        }
        let e_net = quad.max_val(T::from_f64(1e-12)).sqrt();
        let alpha_cva = s.cpty_lgd * (T::one() - (-s.cpty_hazard_rate * s.horizon_years).exp());
        let mut euler = [T::zero(); 3];
        for i in 0..3 {
            let mut row = T::zero();
            for j in 0..3 {
                row += s.exposure_corr_3x3[i][j] * s.subportfolio_exposures[j];
            }
            euler[i] = alpha_cva * s.subportfolio_exposures[i] * row / e_net;
        }
        let cva_net = alpha_cva * e_net;
        let cva_stand = alpha_cva
            * (s.subportfolio_exposures[0]
                + s.subportfolio_exposures[1]
                + s.subportfolio_exposures[2]);
        let diff = cva_net - (euler[0] + euler[1] + euler[2]);
        GregorySaCvaEulerNettingReport {
            sa_cva_delta_capital: kd,
            sa_cva_vega_capital: kv,
            total_sa_cva_capital: kd + kv,
            standalone_cva_sum: cva_stand,
            netted_portfolio_cva: cva_net,
            netting_cva_reduction_benefit: cva_stand - cva_net,
            euler_allocated_cva: euler,
            netting_set_euler_allocation_error: (diff * diff + T::from_f64(1e-30)).sqrt(),
        }
    }

    #[cfg(test)]
    mod domain_suite_47_tests {
        use super::*;
        use ql_ad::Dual;
        fn spec<T: Real>(x: T) -> GregorySaCvaEulerNettingSpec<T> {
            GregorySaCvaEulerNettingSpec {
                cva_delta_sensitivities: [T::from_f64(10.0) + x, T::from_f64(8.0)],
                cva_hedge_sensitivities: [T::from_f64(6.0), T::from_f64(5.0)],
                delta_risk_weights: [T::from_f64(0.05), T::from_f64(0.06)],
                intra_bucket_corr_rho: T::from_f64(0.5),
                cva_vega_sensitivity: T::from_f64(4.0),
                vega_risk_weight: T::from_f64(0.1),
                subportfolio_exposures: [
                    T::from_f64(20.0) + x,
                    T::from_f64(15.0),
                    T::from_f64(10.0),
                ],
                exposure_corr_3x3: [
                    [T::one(), T::from_f64(0.3), T::from_f64(0.1)],
                    [T::from_f64(0.3), T::one(), T::from_f64(0.2)],
                    [T::from_f64(0.1), T::from_f64(0.2), T::one()],
                ],
                cpty_hazard_rate: T::from_f64(0.02),
                cpty_lgd: T::from_f64(0.6),
                horizon_years: T::from_f64(3.0),
            }
        }
        #[test]
        fn suite47_xva_invariants_and_dual_fd() {
            let r = evaluate_gregory_sa_cva_and_euler_netting_allocation(&spec(0.0));
            assert!(r.netting_set_euler_allocation_error < 1e-12);
            assert!(r.netting_cva_reduction_benefit > 0.0);
            let h = 1e-6;
            let fd = (evaluate_gregory_sa_cva_and_euler_netting_allocation(&spec(h))
                .total_sa_cva_capital
                - evaluate_gregory_sa_cva_and_euler_netting_allocation(&spec(-h))
                    .total_sa_cva_capital)
                / (2.0 * h);
            let ad = evaluate_gregory_sa_cva_and_euler_netting_allocation(&spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .total_sa_cva_capital;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_gregory_sa_cva_and_euler_netting_allocation_module::*;

mod evaluate_gregory_hull_white_wrong_way_risk_cva_module {
    // Gregory (2020) & Brigo–Morini–Pallavicini (2013) Hull–White Wrong-Way / Right-Way Risk
    // (`WWR`/`RWR`) Hazard-Exposure Co-Movement & CVA Multiplier Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryHullWhiteWwrCvaSpec<T: Real = f64> {
        pub baseline_exposure_mu: T,
        pub exposure_volatility_sigma_x: T,
        pub baseline_hazard_rate_lambda0: T,
        pub log_intensity_vol_sigma_lambda: T,
        pub exposure_credit_corr_rho: T,
        pub recovery_rate: T,
        pub risk_free_rate: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryHullWhiteWwrCvaReport<T: Real = f64> {
        pub independent_epe: T,
        pub wwr_conditional_epe: T,
        pub independent_cva: T,
        pub wwr_adjusted_cva: T,
        pub wwr_cva_multiplier: T,
        pub effective_wwr_hazard_rate: T,
    }

    fn bachelier_epe<T: Real>(mean: T, std: T) -> T {
        let s = std.max_val(T::from_f64(1e-8));
        let d = mean / s;
        mean * norm_cdf(d) + s * norm_pdf(d)
    }

    pub fn evaluate_gregory_hull_white_wrong_way_risk_cva<T: Real>(
        spec: &GregoryHullWhiteWwrCvaSpec<T>,
    ) -> GregoryHullWhiteWwrCvaReport<T> {
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));
        let half_t = T::from_f64(0.5) * t;
        let std_x = spec.exposure_volatility_sigma_x * half_t.sqrt();

        // Independent Expected Positive Exposure at representative midpoint t/2:
        let epe_ind = bachelier_epe(spec.baseline_exposure_mu, std_x);

        // Under Hull-White / Brigo lognormal hazard linkage lambda(t) = lambda_0 exp(sigma_lambda W^lambda_t - 0.5 sigma_lambda^2 t^2),
        // conditioning on default at time u shifts the Gaussian exposure drift by
        // Delta_mu(u) = rho_{X,lambda} * sigma_X * sigma_lambda * u:
        let drift_shift = spec.exposure_credit_corr_rho
            * spec.exposure_volatility_sigma_x
            * spec.log_intensity_vol_sigma_lambda
            * half_t;
        let epe_wwr = bachelier_epe(spec.baseline_exposure_mu + drift_shift, std_x);

        let lgd = T::one() - spec.recovery_rate;
        let df = (-spec.risk_free_rate * half_t).exp();
        let pd = T::one() - (-spec.baseline_hazard_rate_lambda0 * t).exp();

        let cva_ind = df * lgd * pd * epe_ind;
        let cva_wwr = df * lgd * pd * epe_wwr;
        let mult = cva_wwr / cva_ind.max_val(T::from_f64(1e-12));
        let eff_lambda = spec.baseline_hazard_rate_lambda0 * mult;

        GregoryHullWhiteWwrCvaReport {
            independent_epe: epe_ind,
            wwr_conditional_epe: epe_wwr,
            independent_cva: cva_ind,
            wwr_adjusted_cva: cva_wwr,
            wwr_cva_multiplier: mult,
            effective_wwr_hazard_rate: eff_lambda,
        }
    }

    #[cfg(test)]
    mod domain_suite_51_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gregory_hull_white_wwr_cva_invariants_and_ad_parity() {
            let spec = GregoryHullWhiteWwrCvaSpec {
                baseline_exposure_mu: 12.0,
                exposure_volatility_sigma_x: 18.0,
                baseline_hazard_rate_lambda0: 0.025,
                log_intensity_vol_sigma_lambda: 0.65,
                exposure_credit_corr_rho: 0.55,
                recovery_rate: 0.40,
                risk_free_rate: 0.035,
                maturity_years_t: 3.0,
            };
            let rep = evaluate_gregory_hull_white_wrong_way_risk_cva(&spec);
            // Positive exposure-credit correlation (Wrong-Way Risk) strictly amplifies EPE and CVA:
            assert!(rep.wwr_conditional_epe > rep.independent_epe);
            assert!(rep.wwr_cva_multiplier > 1.0);

            let h = 1e-6;
            let r_ad =
                evaluate_gregory_hull_white_wrong_way_risk_cva(&GregoryHullWhiteWwrCvaSpec {
                    baseline_exposure_mu: Dual::<1>::from_value(spec.baseline_exposure_mu),
                    exposure_volatility_sigma_x: Dual::<1>::from_value(
                        spec.exposure_volatility_sigma_x,
                    ),
                    baseline_hazard_rate_lambda0: Dual::<1>::from_value(
                        spec.baseline_hazard_rate_lambda0,
                    ),
                    log_intensity_vol_sigma_lambda: Dual::<1>::from_value(
                        spec.log_intensity_vol_sigma_lambda,
                    ),
                    exposure_credit_corr_rho: Dual::<1>::variable(spec.exposure_credit_corr_rho, 0),
                    recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
                    risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                    maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                });
            let mut sp = spec;
            sp.exposure_credit_corr_rho += h;
            let mut sm = spec;
            sm.exposure_credit_corr_rho -= h;
            let fd = (evaluate_gregory_hull_white_wrong_way_risk_cva(&sp).wwr_adjusted_cva
                - evaluate_gregory_hull_white_wrong_way_risk_cva(&sm).wwr_adjusted_cva)
                / (2.0 * h);
            assert!((r_ad.wwr_adjusted_cva.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_gregory_hull_white_wrong_way_risk_cva_module::*;

mod evaluate_gregory_albanese_fva_fca_fba_overlap_module {
    // Gregory (2020) & Albanese–Andersen (2014) Asymmetric Funding Valuation Adjustment
    // (`FCA`/`FBA`) & DVA Bond-CDS Funding-Spread Overlap Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryAlbaneseFvaOverlapSpec<T: Real = f64> {
        pub expected_positive_exposure_epe: T,
        pub expected_negative_exposure_ene: T,
        pub bank_cds_hazard_rate_lambda_b: T,
        pub bank_recovery_rate_rb: T,
        pub bond_cds_liquidity_basis_s_liq: T,
        pub treasury_deposit_spread_s_dep: T,
        pub counterparty_survival_weight: T,
        pub risk_free_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryAlbaneseFvaOverlapReport<T: Real = f64> {
        pub bank_default_credit_spread: T,
        pub all_in_unsecured_borrowing_spread: T,
        pub gross_funding_cost_adjustment_fca: T,
        pub gross_funding_benefit_adjustment_fba: T,
        pub unilateral_dva_benefit: T,
        pub dva_fba_double_count_overlap: T,
        pub non_overlapping_net_fva_plus_dva: T,
    }

    pub fn evaluate_gregory_albanese_fva_fca_fba_overlap<T: Real>(
        spec: &GregoryAlbaneseFvaOverlapSpec<T>,
    ) -> GregoryAlbaneseFvaOverlapReport<T> {
        let epe = spec.expected_positive_exposure_epe.max_val(T::zero());
        let ene = spec.expected_negative_exposure_ene.max_val(T::zero());
        let lgd_b = (T::one() - spec.bank_recovery_rate_rb).max_val(T::from_f64(0.05));
        let lam_b = spec.bank_cds_hazard_rate_lambda_b.max_val(T::zero());
        let s_def = lam_b * lgd_b;
        let s_borrow = s_def + spec.bond_cds_liquidity_basis_s_liq;
        let s_lend = spec.treasury_deposit_spread_s_dep;

        let tot_disc = (spec.risk_free_rate_r + lam_b).max_val(T::from_f64(1e-4));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));
        let risky_annuity =
            spec.counterparty_survival_weight * (T::one() - (-tot_disc * t).exp()) / tot_disc;

        let fca = s_borrow * epe * risky_annuity;
        let fba = s_lend * ene * risky_annuity;
        let dva = s_def * ene * risky_annuity;

        // Gregory (2020, Ch. 15) / Burgard-Kjaer: smooth harmonic overlap between bank default spread and lending spread:
        let s_def_pos = s_def.max_val(T::from_f64(1e-6));
        let s_lend_pos = s_lend.max_val(T::from_f64(1e-6));
        let overlap_spread = (s_def_pos * s_lend_pos) / (s_def_pos + s_lend_pos);
        let overlap = overlap_spread * ene * risky_annuity;
        let net_fva_dva = fca - fba - dva + overlap;

        GregoryAlbaneseFvaOverlapReport {
            bank_default_credit_spread: s_def,
            all_in_unsecured_borrowing_spread: s_borrow,
            gross_funding_cost_adjustment_fca: fca,
            gross_funding_benefit_adjustment_fba: fba,
            unilateral_dva_benefit: dva,
            dva_fba_double_count_overlap: overlap,
            non_overlapping_net_fva_plus_dva: net_fva_dva,
        }
    }

    #[cfg(test)]
    mod domain_suite_54_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gregory_albanese_fva_fca_fba_overlap_invariants_and_ad_parity() {
            let spec = GregoryAlbaneseFvaOverlapSpec {
                expected_positive_exposure_epe: 42.0,
                expected_negative_exposure_ene: 28.0,
                bank_cds_hazard_rate_lambda_b: 0.015,
                bank_recovery_rate_rb: 0.40,
                bond_cds_liquidity_basis_s_liq: 0.0045,
                treasury_deposit_spread_s_dep: 0.0035,
                counterparty_survival_weight: 0.94,
                risk_free_rate_r: 0.038,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_gregory_albanese_fva_fca_fba_overlap(&spec);
            assert!(
                rep.gross_funding_cost_adjustment_fca > rep.gross_funding_benefit_adjustment_fba
            );
            assert!(rep.dva_fba_double_count_overlap > 0.0);

            let h = 1e-6;
            let r_ad =
                evaluate_gregory_albanese_fva_fca_fba_overlap(&GregoryAlbaneseFvaOverlapSpec {
                    expected_positive_exposure_epe: Dual::<1>::variable(
                        spec.expected_positive_exposure_epe,
                        0,
                    ),
                    expected_negative_exposure_ene: Dual::<1>::from_value(
                        spec.expected_negative_exposure_ene,
                    ),
                    bank_cds_hazard_rate_lambda_b: Dual::<1>::from_value(
                        spec.bank_cds_hazard_rate_lambda_b,
                    ),
                    bank_recovery_rate_rb: Dual::<1>::from_value(spec.bank_recovery_rate_rb),
                    bond_cds_liquidity_basis_s_liq: Dual::<1>::from_value(
                        spec.bond_cds_liquidity_basis_s_liq,
                    ),
                    treasury_deposit_spread_s_dep: Dual::<1>::from_value(
                        spec.treasury_deposit_spread_s_dep,
                    ),
                    counterparty_survival_weight: Dual::<1>::from_value(
                        spec.counterparty_survival_weight,
                    ),
                    risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
                    maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                });
            let mut sp = spec;
            sp.expected_positive_exposure_epe += h;
            let mut sm = spec;
            sm.expected_positive_exposure_epe -= h;
            let fd = (evaluate_gregory_albanese_fva_fca_fba_overlap(&sp)
                .non_overlapping_net_fva_plus_dva
                - evaluate_gregory_albanese_fva_fca_fba_overlap(&sm)
                    .non_overlapping_net_fva_plus_dva)
                / (2.0 * h);
            assert!((r_ad.non_overlapping_net_fva_plus_dva.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_gregory_albanese_fva_fca_fba_overlap_module::*;

mod evaluate_gregory_pykhtin_csa_threshold_mta_cva_module {
    // Gregory (2020) & Pykhtin–Rosen (2010) Bilateral CSA Threshold, Minimum Transfer
    // Amount (`MTA`), Independent Amount (`IA`) & MPOR Collateralized `CVA` Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryPykhtinCsaThresholdMtaCvaSpec<T: Real = f64> {
        pub uncollateralized_exposure_vol_sigma_v: T,
        pub current_mtm_value_v0: T,
        pub csa_threshold_h: T,
        pub minimum_transfer_amount_mta: T,
        pub independent_amount_ia: T,
        pub margin_period_of_risk_days: T,
        pub counterparty_hazard_rate_lambda: T,
        pub recovery_rate_r: T,
        pub ois_discount_rate_r: T,
        pub horizon_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryPykhtinCsaThresholdMtaCvaReport<T: Real = f64> {
        pub uncollateralized_epe: T,
        pub effective_csa_threshold_buffer: T,
        pub mpor_residual_epe_floor: T,
        pub mta_leakage_epe_addon: T,
        pub collateralized_epe: T,
        pub uncollateralized_cva: T,
        pub collateralized_cva: T,
        pub csa_cva_reduction_ratio: T,
    }

    fn folded_normal_pos_mean<T: Real>(mu: T, s: T) -> T {
        let std = s.max_val(T::from_f64(1e-6));
        let z = mu / std;
        mu * norm_cdf(z) + std * norm_pdf(z)
    }

    pub fn evaluate_gregory_pykhtin_csa_threshold_mta_cva<T: Real>(
        spec: &GregoryPykhtinCsaThresholdMtaCvaSpec<T>,
    ) -> GregoryPykhtinCsaThresholdMtaCvaReport<T> {
        let sig_v = spec
            .uncollateralized_exposure_vol_sigma_v
            .max_val(T::from_f64(1e-4));
        let v0 = spec.current_mtm_value_v0;
        let t = spec.horizon_years_t.max_val(T::from_f64(0.25));
        let dt_mpor = spec.margin_period_of_risk_days.max_val(T::one()) / T::from_f64(252.0);
        let s_mpor = sig_v * dt_mpor.sqrt();

        // Smooth effective threshold H_eff = softplus(H + 0.5 * MTA - IA)
        let raw_h = spec.csa_threshold_h + T::from_f64(0.5) * spec.minimum_transfer_amount_mta
            - spec.independent_amount_ia;
        // Use a capital-scale smoothing width so the leverage constraint remains
        // differentiable even when the RWA leg is materially dominant.
        let eps = T::from_f64(1.0);
        let h_eff = T::from_f64(0.5) * (raw_h + (raw_h * raw_h + eps * eps).sqrt());
        let raw_h_no_mta = spec.csa_threshold_h - spec.independent_amount_ia;
        let h_no_mta =
            T::from_f64(0.5) * (raw_h_no_mta + (raw_h_no_mta * raw_h_no_mta + eps * eps).sqrt());

        let time_fracs = [0.25, 0.55, 0.85];
        let weights = [0.30, 0.40, 0.30];
        let inv_sqrt_2pi = T::from_f64(0.398_942_280_401_432_7);

        let mut epe_uncoll = T::zero();
        let mut epe_coll = T::zero();
        let mut epe_no_mta = T::zero();
        let mut epe_mpor_only = T::zero();

        for (&tf, &wf) in time_fracs.iter().zip(weights.iter()) {
            let tk = t * T::from_f64(tf);
            let wk = T::from_f64(wf);
            let sk = sig_v * tk.sqrt();

            let ee_u = folded_normal_pos_mean(v0, sk);
            let p_above_h = norm_cdf((v0 - h_eff) / sk);
            let p_above_h0 = norm_cdf((v0 - h_no_mta) / sk);
            let p_above_0 = norm_cdf(v0 / sk);

            let ee_c =
                ee_u - folded_normal_pos_mean(v0 - h_eff, sk) + p_above_h * s_mpor * inv_sqrt_2pi;
            let ee_c0 = ee_u - folded_normal_pos_mean(v0 - h_no_mta, sk)
                + p_above_h0 * s_mpor * inv_sqrt_2pi;
            let ee_mpor = p_above_0 * s_mpor * inv_sqrt_2pi;

            epe_uncoll += wk * ee_u;
            epe_coll += wk * ee_c;
            epe_no_mta += wk * ee_c0;
            epe_mpor_only += wk * ee_mpor;
        }

        let mta_addon = (epe_coll - epe_no_mta).max_val(T::zero());
        let lgd = (T::one() - spec.recovery_rate_r).max_val(T::from_f64(0.05));
        let lam = spec
            .counterparty_hazard_rate_lambda
            .max_val(T::from_f64(1e-4));
        let pd = T::one() - (-lam * t).exp();
        let df = (-spec.ois_discount_rate_r * t * T::from_f64(0.5)).exp();

        let cva_uncoll = lgd * pd * df * epe_uncoll;
        let cva_coll = lgd * pd * df * epe_coll;
        let red_ratio = (cva_uncoll - cva_coll) / cva_uncoll.max_val(T::from_f64(1e-6));

        GregoryPykhtinCsaThresholdMtaCvaReport {
            uncollateralized_epe: epe_uncoll,
            effective_csa_threshold_buffer: h_eff,
            mpor_residual_epe_floor: epe_mpor_only,
            mta_leakage_epe_addon: mta_addon,
            collateralized_epe: epe_coll,
            uncollateralized_cva: cva_uncoll,
            collateralized_cva: cva_coll,
            csa_cva_reduction_ratio: red_ratio,
        }
    }

    #[cfg(test)]
    mod domain_suite_55_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gregory_pykhtin_csa_threshold_mta_cva_invariants_and_ad_parity() {
            let spec = GregoryPykhtinCsaThresholdMtaCvaSpec {
                uncollateralized_exposure_vol_sigma_v: 35.0,
                current_mtm_value_v0: 12.0,
                csa_threshold_h: 15.0,
                minimum_transfer_amount_mta: 2.5,
                independent_amount_ia: 4.0,
                margin_period_of_risk_days: 10.0,
                counterparty_hazard_rate_lambda: 0.022,
                recovery_rate_r: 0.40,
                ois_discount_rate_r: 0.038,
                horizon_years_t: 3.0,
            };
            let rep = evaluate_gregory_pykhtin_csa_threshold_mta_cva(&spec);
            assert!(rep.collateralized_epe < rep.uncollateralized_epe);
            assert!(rep.collateralized_epe > rep.mpor_residual_epe_floor);
            assert!(rep.mta_leakage_epe_addon > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_gregory_pykhtin_csa_threshold_mta_cva(
                &GregoryPykhtinCsaThresholdMtaCvaSpec {
                    uncollateralized_exposure_vol_sigma_v: Dual::<1>::from_value(
                        spec.uncollateralized_exposure_vol_sigma_v,
                    ),
                    current_mtm_value_v0: Dual::<1>::variable(spec.current_mtm_value_v0, 0),
                    csa_threshold_h: Dual::<1>::from_value(spec.csa_threshold_h),
                    minimum_transfer_amount_mta: Dual::<1>::from_value(
                        spec.minimum_transfer_amount_mta,
                    ),
                    independent_amount_ia: Dual::<1>::from_value(spec.independent_amount_ia),
                    margin_period_of_risk_days: Dual::<1>::from_value(
                        spec.margin_period_of_risk_days,
                    ),
                    counterparty_hazard_rate_lambda: Dual::<1>::from_value(
                        spec.counterparty_hazard_rate_lambda,
                    ),
                    recovery_rate_r: Dual::<1>::from_value(spec.recovery_rate_r),
                    ois_discount_rate_r: Dual::<1>::from_value(spec.ois_discount_rate_r),
                    horizon_years_t: Dual::<1>::from_value(spec.horizon_years_t),
                },
            );
            let mut sp = spec;
            sp.current_mtm_value_v0 += h;
            let mut sm = spec;
            sm.current_mtm_value_v0 -= h;
            let fd = (evaluate_gregory_pykhtin_csa_threshold_mta_cva(&sp).collateralized_cva
                - evaluate_gregory_pykhtin_csa_threshold_mta_cva(&sm).collateralized_cva)
                / (2.0 * h);
            assert!((r_ad.collateralized_cva.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_gregory_pykhtin_csa_threshold_mta_cva_module::*;

mod evaluate_gregory_break_clause_extinguishable_cva_module {
    // Gregory (2020) & Brigo–Mercurio (2006) Mandatory/Optional Early-Termination
    // Break-Clause (`ATB`) & Extinguishable Swap `CVA` Reduction Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryBreakClauseExtinguishableCvaSpec<T: Real = f64> {
        pub unmitigated_swap_epe_profile: [T; 3],
        pub break_clause_date_years_tb: T,
        pub final_maturity_years_t: T,
        pub counterparty_hazard_rate_lambda: T,
        pub break_execution_probability_p_exec: T,
        pub break_settlement_slippage_fraction: T,
        pub recovery_rate_r: T,
        pub ois_discount_rate_r: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryBreakClauseExtinguishableCvaReport<T: Real = f64> {
        pub unmitigated_full_horizon_cva: T,
        pub pre_break_truncated_cva: T,
        pub post_break_unexercised_leakage_cva: T,
        pub break_clause_mitigated_cva: T,
        pub extinguishable_swap_survival_discount: T,
        pub cva_reduction_from_break_clause: T,
        pub break_clause_hedge_effectiveness: T,
    }

    pub fn evaluate_gregory_break_clause_extinguishable_cva<T: Real>(
        spec: &GregoryBreakClauseExtinguishableCvaSpec<T>,
    ) -> GregoryBreakClauseExtinguishableCvaReport<T> {
        let e0 = spec.unmitigated_swap_epe_profile[0].max_val(T::from_f64(1e-4));
        let e1 = spec.unmitigated_swap_epe_profile[1].max_val(T::from_f64(1e-4));
        let e2 = spec.unmitigated_swap_epe_profile[2].max_val(T::from_f64(1e-4));
        let tb = spec.break_clause_date_years_tb.max_val(T::from_f64(0.25));
        let t_mat = spec.final_maturity_years_t.max_val(tb + T::from_f64(0.5));
        let lam = spec
            .counterparty_hazard_rate_lambda
            .max_val(T::from_f64(1e-4));
        let p_exec = spec
            .break_execution_probability_p_exec
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.99));
        let slip = spec
            .break_settlement_slippage_fraction
            .max_val(T::zero())
            .min_val(T::from_f64(0.50));
        let lgd = (T::one() - spec.recovery_rate_r).max_val(T::from_f64(0.05));
        let r = spec.ois_discount_rate_r.max_val(T::from_f64(1e-4));

        let t1 = T::from_f64(0.5) * tb;
        let t2 = tb;
        let t3 = T::from_f64(0.5) * (tb + t_mat);

        let s0 = T::one();
        let s1 = (-lam * t1).exp();
        let s2 = (-lam * t2).exp();
        let s3 = (-lam * t_mat).exp();

        let d1 = (-r * t1).exp() * (s0 - s1);
        let d2 = (-r * t2).exp() * (s1 - s2);
        let d3 = (-r * t3).exp() * (s2 - s3);

        let cva_p1 = lgd * e0 * d1;
        let cva_p2 = lgd * e1 * d2;
        let cva_p3 = lgd * e2 * d3;

        let cva_unmit = cva_p1 + cva_p2 + cva_p3;
        let cva_pre_break = cva_p1 + cva_p2;

        // Post-break exposure survives with probability (1 - p_exec) plus close-out settlement slippage scaled by EPE(tb)/EPE_post:
        let closeout_ratio = (T::from_f64(2.0) * e1) / (e1 + e2);
        let post_weight =
            ((T::one() - p_exec) + p_exec * slip * closeout_ratio).min_val(T::from_f64(0.99));
        let cva_post_leakage = cva_p3 * post_weight;
        let cva_mitigated = cva_pre_break + cva_post_leakage;

        let cva_saved = cva_unmit - cva_mitigated;
        let eff = cva_saved / cva_unmit.max_val(T::from_f64(1e-6));
        // Exact risky OIS survival/extinguishment discount factor at T:
        let ext_surv = (-r * t_mat).exp() * (p_exec * s2 + (T::one() - p_exec) * s3);

        GregoryBreakClauseExtinguishableCvaReport {
            unmitigated_full_horizon_cva: cva_unmit,
            pre_break_truncated_cva: cva_pre_break,
            post_break_unexercised_leakage_cva: cva_post_leakage,
            break_clause_mitigated_cva: cva_mitigated,
            extinguishable_swap_survival_discount: ext_surv,
            cva_reduction_from_break_clause: cva_saved,
            break_clause_hedge_effectiveness: eff,
        }
    }

    #[cfg(test)]
    mod evaluate_gregory_break_clause_extinguishable_cva_module_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gregory_break_clause_extinguishable_cva_invariants_and_ad_parity() {
            let spec = GregoryBreakClauseExtinguishableCvaSpec {
                unmitigated_swap_epe_profile: [18.0, 32.0, 45.0],
                break_clause_date_years_tb: 2.0,
                final_maturity_years_t: 5.0,
                counterparty_hazard_rate_lambda: 0.025,
                break_execution_probability_p_exec: 0.85,
                break_settlement_slippage_fraction: 0.08,
                recovery_rate_r: 0.40,
                ois_discount_rate_r: 0.038,
            };
            let rep = evaluate_gregory_break_clause_extinguishable_cva(&spec);
            assert!(rep.break_clause_mitigated_cva < rep.unmitigated_full_horizon_cva);
            assert!(rep.cva_reduction_from_break_clause > 0.0);
            assert!(rep.break_clause_hedge_effectiveness > 0.30);

            let h = 1e-6;
            let r_ad = evaluate_gregory_break_clause_extinguishable_cva(
                &GregoryBreakClauseExtinguishableCvaSpec {
                    unmitigated_swap_epe_profile: [
                        Dual::<10>::variable(spec.unmitigated_swap_epe_profile[0], 0),
                        Dual::<10>::variable(spec.unmitigated_swap_epe_profile[1], 1),
                        Dual::<10>::variable(spec.unmitigated_swap_epe_profile[2], 2),
                    ],
                    break_clause_date_years_tb: Dual::<10>::variable(
                        spec.break_clause_date_years_tb,
                        3,
                    ),
                    final_maturity_years_t: Dual::<10>::variable(spec.final_maturity_years_t, 4),
                    counterparty_hazard_rate_lambda: Dual::<10>::variable(
                        spec.counterparty_hazard_rate_lambda,
                        5,
                    ),
                    break_execution_probability_p_exec: Dual::<10>::variable(
                        spec.break_execution_probability_p_exec,
                        6,
                    ),
                    break_settlement_slippage_fraction: Dual::<10>::variable(
                        spec.break_settlement_slippage_fraction,
                        7,
                    ),
                    recovery_rate_r: Dual::<10>::variable(spec.recovery_rate_r, 8),
                    ois_discount_rate_r: Dual::<10>::variable(spec.ois_discount_rate_r, 9),
                },
            );
            for (idx, &ad_deriv) in r_ad.break_clause_mitigated_cva.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.unmitigated_swap_epe_profile[0] += h;
                        sm.unmitigated_swap_epe_profile[0] -= h;
                    }
                    1 => {
                        sp.unmitigated_swap_epe_profile[1] += h;
                        sm.unmitigated_swap_epe_profile[1] -= h;
                    }
                    2 => {
                        sp.unmitigated_swap_epe_profile[2] += h;
                        sm.unmitigated_swap_epe_profile[2] -= h;
                    }
                    3 => {
                        sp.break_clause_date_years_tb += h;
                        sm.break_clause_date_years_tb -= h;
                    }
                    4 => {
                        sp.final_maturity_years_t += h;
                        sm.final_maturity_years_t -= h;
                    }
                    5 => {
                        sp.counterparty_hazard_rate_lambda += h;
                        sm.counterparty_hazard_rate_lambda -= h;
                    }
                    6 => {
                        sp.break_execution_probability_p_exec += h;
                        sm.break_execution_probability_p_exec -= h;
                    }
                    7 => {
                        sp.break_settlement_slippage_fraction += h;
                        sm.break_settlement_slippage_fraction -= h;
                    }
                    8 => {
                        sp.recovery_rate_r += h;
                        sm.recovery_rate_r -= h;
                    }
                    _ => {
                        sp.ois_discount_rate_r += h;
                        sm.ois_discount_rate_r -= h;
                    }
                }
                let fd = (evaluate_gregory_break_clause_extinguishable_cva(&sp)
                    .break_clause_mitigated_cva
                    - evaluate_gregory_break_clause_extinguishable_cva(&sm)
                        .break_clause_mitigated_cva)
                    / (2.0 * h);
                assert!((ad_deriv - fd).abs() < 1e-5);
            }
        }
    }
}
pub use evaluate_gregory_break_clause_extinguishable_cva_module::*;

mod evaluate_gregory_kenyon_simm_mva_funding_module {
    // Gregory (2020, Ch. 16) & Green–Kenyon (2015) Initial-Margin Valuation Adjustment
    // (`MVA`) under ISDA SIMM vs. CCP Clearing with Joint Survival & Custodian Haircut Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryKenyonSimmMvaFundingSpec<T: Real = f64> {
        pub initial_simm_margin_im0: T,
        pub ccp_cleared_margin_ratio: T,
        pub bank_funding_spread_s_f: T,
        pub custodian_remuneration_spread_s_c: T,
        pub bank_hazard_rate_lambda_b: T,
        pub counterparty_hazard_rate_lambda_c: T,
        pub ois_discount_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryKenyonSimmMvaFundingReport<T: Real = f64> {
        pub bilateral_simm_expected_im_pv: T,
        pub bilateral_simm_mva_charge: T,
        pub ccp_cleared_mva_charge: T,
        pub mva_savings_from_central_clearing: T,
        pub static_peak_im_mva_overstatement: T,
        pub effective_annualized_mva_drag_bps: T,
    }

    pub fn evaluate_gregory_kenyon_simm_mva_funding<T: Real>(
        spec: &GregoryKenyonSimmMvaFundingSpec<T>,
    ) -> GregoryKenyonSimmMvaFundingReport<T> {
        let im0 = spec.initial_simm_margin_im0.max_val(T::from_f64(1e-4));
        let ccp_ratio = spec
            .ccp_cleared_margin_ratio
            .max_val(T::from_f64(0.20))
            .min_val(T::from_f64(1.20));
        let s_f = spec.bank_funding_spread_s_f.max_val(T::from_f64(1e-4));
        let s_c = spec.custodian_remuneration_spread_s_c;
        let net_im_spread = (s_f - s_c).max_val(T::from_f64(1e-4));
        let lb = spec.bank_hazard_rate_lambda_b.max_val(T::from_f64(1e-4));
        let lc = spec
            .counterparty_hazard_rate_lambda_c
            .max_val(T::from_f64(1e-4));
        let r = spec.ois_discount_rate_r.max_val(T::from_f64(1e-4));
        let t_mat = spec.maturity_years_t.max_val(T::from_f64(0.5));

        // 4-step Gauss-midpoint integration of decaying SIMM profile EIM(t) = IM0 * sqrt((T - t)/T) * (1 + 0.15 t/T):
        let dt = t_mat / T::from_f64(4.0);
        let mut eim_pv = T::zero();
        let mut static_im_pv = T::zero();
        let eff_disc = r + lb + lc;

        for k in 0..4 {
            let tk = dt * (T::from_f64(k as f64) + T::from_f64(0.5));
            let rem_frac = ((t_mat - tk) / t_mat).max_val(T::from_f64(0.01));
            let profile_scale = rem_frac.sqrt() * (T::one() + T::from_f64(0.15) * (tk / t_mat));
            let surv_df = (-eff_disc * tk).exp();
            eim_pv += dt * surv_df * im0 * profile_scale;
            static_im_pv += dt * surv_df * im0;
        }

        let mva_bilat = eim_pv * net_im_spread;
        let mva_ccp = eim_pv * ccp_ratio * net_im_spread;
        let mva_saved = mva_bilat - mva_ccp;
        let static_over = (static_im_pv - eim_pv) * net_im_spread;
        let drag_bps = (mva_bilat / (im0 * t_mat)) * T::from_f64(10_000.0);

        GregoryKenyonSimmMvaFundingReport {
            bilateral_simm_expected_im_pv: eim_pv,
            bilateral_simm_mva_charge: mva_bilat,
            ccp_cleared_mva_charge: mva_ccp,
            mva_savings_from_central_clearing: mva_saved,
            static_peak_im_mva_overstatement: static_over,
            effective_annualized_mva_drag_bps: drag_bps,
        }
    }

    #[cfg(test)]
    mod evaluate_gregory_kenyon_simm_mva_funding_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gregory_kenyon_simm_mva_funding_invariants_and_ad_parity() {
            let spec = GregoryKenyonSimmMvaFundingSpec {
                initial_simm_margin_im0: 42.0,
                ccp_cleared_margin_ratio: 0.62,
                bank_funding_spread_s_f: 0.0095,
                custodian_remuneration_spread_s_c: 0.0015,
                bank_hazard_rate_lambda_b: 0.012,
                counterparty_hazard_rate_lambda_c: 0.018,
                ois_discount_rate_r: 0.038,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_gregory_kenyon_simm_mva_funding(&spec);
            assert!(rep.bilateral_simm_mva_charge > rep.ccp_cleared_mva_charge);
            assert!(rep.mva_savings_from_central_clearing > 0.0);
            assert!(rep.static_peak_im_mva_overstatement > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_gregory_kenyon_simm_mva_funding(&GregoryKenyonSimmMvaFundingSpec {
                initial_simm_margin_im0: Dual::<8>::variable(spec.initial_simm_margin_im0, 0),
                ccp_cleared_margin_ratio: Dual::<8>::variable(spec.ccp_cleared_margin_ratio, 1),
                bank_funding_spread_s_f: Dual::<8>::variable(spec.bank_funding_spread_s_f, 2),
                custodian_remuneration_spread_s_c: Dual::<8>::variable(
                    spec.custodian_remuneration_spread_s_c,
                    3,
                ),
                bank_hazard_rate_lambda_b: Dual::<8>::variable(spec.bank_hazard_rate_lambda_b, 4),
                counterparty_hazard_rate_lambda_c: Dual::<8>::variable(
                    spec.counterparty_hazard_rate_lambda_c,
                    5,
                ),
                ois_discount_rate_r: Dual::<8>::variable(spec.ois_discount_rate_r, 6),
                maturity_years_t: Dual::<8>::variable(spec.maturity_years_t, 7),
            });
            for (idx, &ad) in r_ad.mva_savings_from_central_clearing.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.initial_simm_margin_im0 += h;
                        sm.initial_simm_margin_im0 -= h;
                    }
                    1 => {
                        sp.ccp_cleared_margin_ratio += h;
                        sm.ccp_cleared_margin_ratio -= h;
                    }
                    2 => {
                        sp.bank_funding_spread_s_f += h;
                        sm.bank_funding_spread_s_f -= h;
                    }
                    3 => {
                        sp.custodian_remuneration_spread_s_c += h;
                        sm.custodian_remuneration_spread_s_c -= h;
                    }
                    4 => {
                        sp.bank_hazard_rate_lambda_b += h;
                        sm.bank_hazard_rate_lambda_b -= h;
                    }
                    5 => {
                        sp.counterparty_hazard_rate_lambda_c += h;
                        sm.counterparty_hazard_rate_lambda_c -= h;
                    }
                    6 => {
                        sp.ois_discount_rate_r += h;
                        sm.ois_discount_rate_r -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_gregory_kenyon_simm_mva_funding(&sp)
                    .mva_savings_from_central_clearing
                    - evaluate_gregory_kenyon_simm_mva_funding(&sm)
                        .mva_savings_from_central_clearing)
                    / (2.0 * h);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_gregory_kenyon_simm_mva_funding_module::*;

mod evaluate_brigo_capponi_cds_contagion_gap_cva_module {
    // Brigo, Capponi & Pallavicini (2014, Math. Finance) & Gregory (2020, Ch. 13)
    // Credit Default Swap (CDS) Counterparty Contagion Jump-to-Default & Gap-Risk CVA Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BrigoCapponiCdsContagionGapCvaSpec<T: Real = f64> {
        pub reference_hazard_rate_lambda_ref: T,
        pub protection_seller_hazard_rate_lambda_c: T,
        pub contractual_running_spread_s_cds: T,
        pub contagion_spread_jump_gamma: T,
        pub simultaneous_default_prob_p_sim: T,
        pub reference_recovery_rate_r_ref: T,
        pub counterparty_recovery_rate_r_c: T,
        pub ois_discount_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BrigoCapponiCdsContagionGapCvaReport<T: Real = f64> {
        pub standalone_par_cds_spread: T,
        pub post_contagion_reference_spread: T,
        pub independent_replacement_cva: T,
        pub contagion_jump_gap_cva_addon: T,
        pub simultaneous_double_default_loss: T,
        pub total_contagion_adjusted_cds_cva: T,
        pub effective_protection_haircut_pct: T,
    }

    pub fn evaluate_brigo_capponi_cds_contagion_gap_cva<T: Real>(
        spec: &BrigoCapponiCdsContagionGapCvaSpec<T>,
    ) -> BrigoCapponiCdsContagionGapCvaReport<T> {
        let lam_r = spec
            .reference_hazard_rate_lambda_ref
            .max_val(T::from_f64(1e-4));
        let lam_c = spec
            .protection_seller_hazard_rate_lambda_c
            .max_val(T::from_f64(1e-4));
        let s_contract = spec.contractual_running_spread_s_cds;
        let gam_cont = spec.contagion_spread_jump_gamma.max_val(T::zero());
        let p_sim = spec
            .simultaneous_default_prob_p_sim
            .max_val(T::zero())
            .min_val(T::from_f64(0.50));
        let lgd_r = (T::one() - spec.reference_recovery_rate_r_ref).max_val(T::from_f64(0.05));
        let lgd_c = (T::one() - spec.counterparty_recovery_rate_r_c).max_val(T::from_f64(0.05));
        let r = spec.ois_discount_rate_r;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.5));

        let s_par = lam_r * lgd_r;
        let lam_r_post = lam_r * (T::one() + gam_cont);
        let s_post = lam_r_post * lgd_r;

        // Risky annuity before and after protection seller default at mid-horizon T/2:
        let tau_rem = T::from_f64(0.5) * t;
        let rpv_pre = (T::one() - (-(r + lam_r) * tau_rem).exp()) / (r + lam_r);
        let rpv_post = (T::one() - (-(r + lam_r_post) * tau_rem).exp()) / (r + lam_r_post);

        let mtm_pre = ((s_par - s_contract) * rpv_pre).max_val(T::from_f64(1e-4));
        let mtm_post = ((s_post - s_contract) * rpv_post).max_val(mtm_pre);

        let lam_tot = lam_c + lam_r;
        let pd_seller_first = (lam_c / lam_tot) * (T::one() - (-lam_tot * t).exp());
        let df_mid = (-r * tau_rem).exp();

        let indep_cva = lgd_c * pd_seller_first * df_mid * mtm_pre;
        let contagion_addon =
            lgd_c * pd_seller_first * (T::one() - p_sim) * df_mid * (mtm_post - mtm_pre);
        let sim_loss = lgd_c * pd_seller_first * p_sim * df_mid * lgd_r;
        let total_cva = indep_cva + contagion_addon + sim_loss;

        let prot_pv = s_par * (T::one() - (-(r + lam_r) * t).exp()) / (r + lam_r);
        let haircut_pct = (total_cva / prot_pv.max_val(T::from_f64(1e-4))) * T::from_f64(100.0);

        BrigoCapponiCdsContagionGapCvaReport {
            standalone_par_cds_spread: s_par,
            post_contagion_reference_spread: s_post,
            independent_replacement_cva: indep_cva,
            contagion_jump_gap_cva_addon: contagion_addon,
            simultaneous_double_default_loss: sim_loss,
            total_contagion_adjusted_cds_cva: total_cva,
            effective_protection_haircut_pct: haircut_pct,
        }
    }

    #[cfg(test)]
    mod evaluate_brigo_capponi_cds_contagion_gap_cva_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn brigo_capponi_cds_contagion_gap_cva_invariants_and_ad_parity() {
            let spec = BrigoCapponiCdsContagionGapCvaSpec {
                reference_hazard_rate_lambda_ref: 0.025,
                protection_seller_hazard_rate_lambda_c: 0.020,
                contractual_running_spread_s_cds: 0.010,
                contagion_spread_jump_gamma: 0.85,
                simultaneous_default_prob_p_sim: 0.06,
                reference_recovery_rate_r_ref: 0.40,
                counterparty_recovery_rate_r_c: 0.40,
                ois_discount_rate_r: 0.038,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_brigo_capponi_cds_contagion_gap_cva(&spec);
            assert!(rep.contagion_jump_gap_cva_addon > 0.0);
            assert!(rep.total_contagion_adjusted_cds_cva > rep.independent_replacement_cva);
            assert!(rep.effective_protection_haircut_pct > 1.0);

            let h = 1e-5;
            let r_ad =
                evaluate_brigo_capponi_cds_contagion_gap_cva(&BrigoCapponiCdsContagionGapCvaSpec {
                    reference_hazard_rate_lambda_ref: Dual::<9>::variable(
                        spec.reference_hazard_rate_lambda_ref,
                        0,
                    ),
                    protection_seller_hazard_rate_lambda_c: Dual::<9>::variable(
                        spec.protection_seller_hazard_rate_lambda_c,
                        1,
                    ),
                    contractual_running_spread_s_cds: Dual::<9>::variable(
                        spec.contractual_running_spread_s_cds,
                        2,
                    ),
                    contagion_spread_jump_gamma: Dual::<9>::variable(
                        spec.contagion_spread_jump_gamma,
                        3,
                    ),
                    simultaneous_default_prob_p_sim: Dual::<9>::variable(
                        spec.simultaneous_default_prob_p_sim,
                        4,
                    ),
                    reference_recovery_rate_r_ref: Dual::<9>::variable(
                        spec.reference_recovery_rate_r_ref,
                        5,
                    ),
                    counterparty_recovery_rate_r_c: Dual::<9>::variable(
                        spec.counterparty_recovery_rate_r_c,
                        6,
                    ),
                    ois_discount_rate_r: Dual::<9>::variable(spec.ois_discount_rate_r, 7),
                    maturity_years_t: Dual::<9>::variable(spec.maturity_years_t, 8),
                });
            for (idx, &ad) in r_ad.total_contagion_adjusted_cds_cva.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.reference_hazard_rate_lambda_ref += h;
                        sm.reference_hazard_rate_lambda_ref -= h;
                    }
                    1 => {
                        sp.protection_seller_hazard_rate_lambda_c += h;
                        sm.protection_seller_hazard_rate_lambda_c -= h;
                    }
                    2 => {
                        sp.contractual_running_spread_s_cds += h;
                        sm.contractual_running_spread_s_cds -= h;
                    }
                    3 => {
                        sp.contagion_spread_jump_gamma += h;
                        sm.contagion_spread_jump_gamma -= h;
                    }
                    4 => {
                        sp.simultaneous_default_prob_p_sim += h;
                        sm.simultaneous_default_prob_p_sim -= h;
                    }
                    5 => {
                        sp.reference_recovery_rate_r_ref += h;
                        sm.reference_recovery_rate_r_ref -= h;
                    }
                    6 => {
                        sp.counterparty_recovery_rate_r_c += h;
                        sm.counterparty_recovery_rate_r_c -= h;
                    }
                    7 => {
                        sp.ois_discount_rate_r += h;
                        sm.ois_discount_rate_r -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_brigo_capponi_cds_contagion_gap_cva(&sp)
                    .total_contagion_adjusted_cds_cva
                    - evaluate_brigo_capponi_cds_contagion_gap_cva(&sm)
                        .total_contagion_adjusted_cds_cva)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_brigo_capponi_cds_contagion_gap_cva_module::*;

mod evaluate_pykhtin_brownian_bridge_mpor_cashflow_spike_epe_module {
    // Pykhtin (2009, J. Credit Risk) & Andersen, Pykhtin & Sokol (2017, Risk, "Rethinking Margin Period of Risk")
    // Margin Period of Risk (MPoR) Brownian-Bridge Unpaid Cashflow Spike & Collateralized EPE Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct PykhtinMporCashflowSpikeEpeSpec<T: Real = f64> {
        pub annualized_portfolio_mtm_volatility_sigma: T,
        pub margin_period_of_risk_years_tau_mpor: T,
        pub csa_threshold_plus_mta: T,
        pub independent_amount_held_ia: T,
        pub contractual_cashflow_paid_in_mpor: T,
        pub post_payment_lookback_weight_theta: T,
        pub counterparty_hazard_rate_lambda: T,
        pub counterparty_recovery_rate_r: T,
        pub ois_discount_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct PykhtinMporCashflowSpikeEpeReport<T: Real = f64> {
        pub mpor_exposure_volatility_sigma_mpor: T,
        pub classical_continuous_collateralized_epe: T,
        pub unpaid_cashflow_spike_collateralized_epe: T,
        pub epe_spike_amplification_ratio: T,
        pub incremental_cashflow_spike_cva: T,
        pub total_mpor_collateralized_cva: T,
    }

    fn bachelier_positive_part_expectation<T: Real>(mu: T, sigma: T) -> T {
        let d = mu / sigma;
        mu * norm_cdf(d) + sigma * norm_pdf(d)
    }

    pub fn evaluate_pykhtin_brownian_bridge_mpor_cashflow_spike_epe<T: Real>(
        spec: &PykhtinMporCashflowSpikeEpeSpec<T>,
    ) -> PykhtinMporCashflowSpikeEpeReport<T> {
        let sig = spec
            .annualized_portfolio_mtm_volatility_sigma
            .max_val(T::from_f64(0.1));
        let tau = spec
            .margin_period_of_risk_years_tau_mpor
            .max_val(T::from_f64(1.0 / 365.0));
        let h_mta = spec.csa_threshold_plus_mta;
        let ia = spec.independent_amount_held_ia;
        let cf_paid = spec.contractual_cashflow_paid_in_mpor.max_val(T::zero());
        let theta = spec
            .post_payment_lookback_weight_theta
            .max_val(T::from_f64(0.05))
            .min_val(T::one());
        let lam = spec
            .counterparty_hazard_rate_lambda
            .max_val(T::from_f64(1e-4));
        let lgd = (T::one() - spec.counterparty_recovery_rate_r).max_val(T::from_f64(0.05));
        let r = spec.ois_discount_rate_r;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

        let sig_mpor = sig * tau.sqrt();
        let mu_classical = h_mta - ia;
        let mu_spike = mu_classical + theta * cf_paid;

        let epe_classical = bachelier_positive_part_expectation(mu_classical, sig_mpor);
        let epe_spike = bachelier_positive_part_expectation(mu_spike, sig_mpor);
        let amp_ratio = epe_spike / epe_classical.max_val(T::from_f64(1e-4));

        let pd = T::one() - (-lam * t).exp();
        let df_mid = (-r * T::from_f64(0.5) * t).exp();
        let cva_classical = lgd * pd * df_mid * epe_classical;
        let cva_total = lgd * pd * df_mid * epe_spike;
        let cva_spike_addon = cva_total - cva_classical;

        PykhtinMporCashflowSpikeEpeReport {
            mpor_exposure_volatility_sigma_mpor: sig_mpor,
            classical_continuous_collateralized_epe: epe_classical,
            unpaid_cashflow_spike_collateralized_epe: epe_spike,
            epe_spike_amplification_ratio: amp_ratio,
            incremental_cashflow_spike_cva: cva_spike_addon,
            total_mpor_collateralized_cva: cva_total,
        }
    }

    #[cfg(test)]
    mod evaluate_pykhtin_brownian_bridge_mpor_cashflow_spike_epe_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn pykhtin_brownian_bridge_mpor_cashflow_spike_invariants_and_ad_parity() {
            let spec = PykhtinMporCashflowSpikeEpeSpec {
                annualized_portfolio_mtm_volatility_sigma: 28.0,
                margin_period_of_risk_years_tau_mpor: 10.0 / 252.0,
                csa_threshold_plus_mta: 2.5,
                independent_amount_held_ia: 1.0,
                contractual_cashflow_paid_in_mpor: 6.0,
                post_payment_lookback_weight_theta: 0.60,
                counterparty_hazard_rate_lambda: 0.022,
                counterparty_recovery_rate_r: 0.40,
                ois_discount_rate_r: 0.038,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_pykhtin_brownian_bridge_mpor_cashflow_spike_epe(&spec);
            assert!(
                rep.unpaid_cashflow_spike_collateralized_epe
                    > rep.classical_continuous_collateralized_epe
            );
            assert!(rep.epe_spike_amplification_ratio > 1.2);
            assert!(rep.incremental_cashflow_spike_cva > 0.0);

            let h = 1e-5;
            let r_ad = evaluate_pykhtin_brownian_bridge_mpor_cashflow_spike_epe(
                &PykhtinMporCashflowSpikeEpeSpec {
                    annualized_portfolio_mtm_volatility_sigma: Dual::<10>::variable(
                        spec.annualized_portfolio_mtm_volatility_sigma,
                        0,
                    ),
                    margin_period_of_risk_years_tau_mpor: Dual::<10>::variable(
                        spec.margin_period_of_risk_years_tau_mpor,
                        1,
                    ),
                    csa_threshold_plus_mta: Dual::<10>::variable(spec.csa_threshold_plus_mta, 2),
                    independent_amount_held_ia: Dual::<10>::variable(
                        spec.independent_amount_held_ia,
                        3,
                    ),
                    contractual_cashflow_paid_in_mpor: Dual::<10>::variable(
                        spec.contractual_cashflow_paid_in_mpor,
                        4,
                    ),
                    post_payment_lookback_weight_theta: Dual::<10>::variable(
                        spec.post_payment_lookback_weight_theta,
                        5,
                    ),
                    counterparty_hazard_rate_lambda: Dual::<10>::variable(
                        spec.counterparty_hazard_rate_lambda,
                        6,
                    ),
                    counterparty_recovery_rate_r: Dual::<10>::variable(
                        spec.counterparty_recovery_rate_r,
                        7,
                    ),
                    ois_discount_rate_r: Dual::<10>::variable(spec.ois_discount_rate_r, 8),
                    maturity_years_t: Dual::<10>::variable(spec.maturity_years_t, 9),
                },
            );
            for (idx, &ad) in r_ad.total_mpor_collateralized_cva.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.annualized_portfolio_mtm_volatility_sigma += h;
                        sm.annualized_portfolio_mtm_volatility_sigma -= h;
                    }
                    1 => {
                        sp.margin_period_of_risk_years_tau_mpor += h;
                        sm.margin_period_of_risk_years_tau_mpor -= h;
                    }
                    2 => {
                        sp.csa_threshold_plus_mta += h;
                        sm.csa_threshold_plus_mta -= h;
                    }
                    3 => {
                        sp.independent_amount_held_ia += h;
                        sm.independent_amount_held_ia -= h;
                    }
                    4 => {
                        sp.contractual_cashflow_paid_in_mpor += h;
                        sm.contractual_cashflow_paid_in_mpor -= h;
                    }
                    5 => {
                        sp.post_payment_lookback_weight_theta += h;
                        sm.post_payment_lookback_weight_theta -= h;
                    }
                    6 => {
                        sp.counterparty_hazard_rate_lambda += h;
                        sm.counterparty_hazard_rate_lambda -= h;
                    }
                    7 => {
                        sp.counterparty_recovery_rate_r += h;
                        sm.counterparty_recovery_rate_r -= h;
                    }
                    8 => {
                        sp.ois_discount_rate_r += h;
                        sm.ois_discount_rate_r -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_pykhtin_brownian_bridge_mpor_cashflow_spike_epe(&sp)
                    .total_mpor_collateralized_cva
                    - evaluate_pykhtin_brownian_bridge_mpor_cashflow_spike_epe(&sm)
                        .total_mpor_collateralized_cva)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_pykhtin_brownian_bridge_mpor_cashflow_spike_epe_module::*;

mod evaluate_gregory_incremental_marginal_netting_set_cva_module {
    // Gregory (2020, The xVA Challenge, 4th ed., Ch. 7 & 14) & Brigo-Morini-Pallavicini (2013)
    // Incremental vs. Euler Marginal Netting-Set CVA Decomposition & Bilateral Netting Benefit Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryIncrementalMarginalNettingCvaSpec<T: Real = f64> {
        pub existing_portfolio_mtm_m1: T,
        pub new_trade_mtm_m2: T,
        pub existing_portfolio_vol_sigma1: T,
        pub new_trade_vol_sigma2: T,
        pub exposure_correlation_rho12: T,
        pub counterparty_hazard_rate_lambda: T,
        pub counterparty_recovery_rate_r: T,
        pub ois_discount_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryIncrementalMarginalNettingCvaReport<T: Real = f64> {
        pub standalone_existing_portfolio_cva: T,
        pub standalone_new_trade_cva: T,
        pub netted_portfolio_total_cva: T,
        pub incremental_trade_cva: T,
        pub euler_marginal_existing_cva: T,
        pub euler_marginal_trade_cva: T,
        pub bilateral_netting_benefit_ratio: T,
    }

    fn bachelier_epe_inner<T: Real>(m: T, s: T) -> T {
        let z = m / s;
        m * norm_cdf(z) + s * norm_pdf(z)
    }

    pub fn evaluate_gregory_incremental_marginal_netting_set_cva<T: Real>(
        spec: &GregoryIncrementalMarginalNettingCvaSpec<T>,
    ) -> GregoryIncrementalMarginalNettingCvaReport<T> {
        let m1 = spec.existing_portfolio_mtm_m1;
        let m2 = spec.new_trade_mtm_m2;
        let sig1 = spec.existing_portfolio_vol_sigma1.max_val(T::from_f64(0.1));
        let sig2 = spec.new_trade_vol_sigma2.max_val(T::from_f64(0.1));
        let rho = spec
            .exposure_correlation_rho12
            .max_val(T::from_f64(-0.98))
            .min_val(T::from_f64(0.98));
        let lam = spec
            .counterparty_hazard_rate_lambda
            .max_val(T::from_f64(1e-4));
        let rec = spec
            .counterparty_recovery_rate_r
            .max_val(T::zero())
            .min_val(T::from_f64(0.95));
        let r = spec.ois_discount_rate_r.max_val(T::from_f64(1e-4));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

        let sqrt_tau = (T::from_f64(0.5) * t).sqrt();
        let s1 = sig1 * sqrt_tau;
        let s2 = sig2 * sqrt_tau;
        let s_net = (s1 * s1 + s2 * s2 + T::from_f64(2.0) * rho * s1 * s2)
            .max_val(T::from_f64(1e-6))
            .sqrt();

        let epe1 = bachelier_epe_inner(m1, s1);
        let epe2 = bachelier_epe_inner(m2, s2);
        let m_net = m1 + m2;
        let epe_net = bachelier_epe_inner(m_net, s_net);

        let rpv01 = (T::one() - (-(r + lam) * t).exp()) / (r + lam);
        let lgd_mult = (T::one() - rec) * lam * rpv01;

        let cva1 = lgd_mult * epe1;
        let cva2 = lgd_mult * epe2;
        let cva_net = lgd_mult * epe_net;
        let inc_cva2 = cva_net - cva1;

        // Euler homogeneous allocation: CVA_1^Euler + CVA_2^Euler == CVA_net
        let z_net = m_net / s_net;
        let cdf_net = norm_cdf(z_net);
        let pdf_net = norm_pdf(z_net);
        let euler1 = lgd_mult * (m1 * cdf_net + ((s1 * s1 + rho * s1 * s2) / s_net) * pdf_net);
        let euler2 = lgd_mult * (m2 * cdf_net + ((s2 * s2 + rho * s1 * s2) / s_net) * pdf_net);

        let nbr = T::one() - epe_net / (epe1 + epe2).max_val(T::from_f64(1e-6));

        GregoryIncrementalMarginalNettingCvaReport {
            standalone_existing_portfolio_cva: cva1,
            standalone_new_trade_cva: cva2,
            netted_portfolio_total_cva: cva_net,
            incremental_trade_cva: inc_cva2,
            euler_marginal_existing_cva: euler1,
            euler_marginal_trade_cva: euler2,
            bilateral_netting_benefit_ratio: nbr,
        }
    }

    #[cfg(test)]
    mod evaluate_gregory_incremental_marginal_netting_set_cva_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gregory_incremental_marginal_netting_cva_invariants_and_ad_parity() {
            let spec = GregoryIncrementalMarginalNettingCvaSpec {
                existing_portfolio_mtm_m1: 12.0,
                new_trade_mtm_m2: -4.0,
                existing_portfolio_vol_sigma1: 18.0,
                new_trade_vol_sigma2: 10.0,
                exposure_correlation_rho12: -0.35,
                counterparty_hazard_rate_lambda: 0.024,
                counterparty_recovery_rate_r: 0.40,
                ois_discount_rate_r: 0.038,
                maturity_years_t: 4.0,
            };
            let rep = evaluate_gregory_incremental_marginal_netting_set_cva(&spec);
            // Euler additivity identity:
            assert!(
                (rep.euler_marginal_existing_cva + rep.euler_marginal_trade_cva
                    - rep.netted_portfolio_total_cva)
                    .abs()
                    < 1e-12
            );
            // Risk-reducing offsetting trade produces negative incremental CVA:
            assert!(rep.incremental_trade_cva < 0.0);
            assert!(rep.bilateral_netting_benefit_ratio > 0.15);

            let h = 1e-5;
            let r_ad = evaluate_gregory_incremental_marginal_netting_set_cva(
                &GregoryIncrementalMarginalNettingCvaSpec {
                    existing_portfolio_mtm_m1: Dual::<9>::variable(
                        spec.existing_portfolio_mtm_m1,
                        0,
                    ),
                    new_trade_mtm_m2: Dual::<9>::variable(spec.new_trade_mtm_m2, 1),
                    existing_portfolio_vol_sigma1: Dual::<9>::variable(
                        spec.existing_portfolio_vol_sigma1,
                        2,
                    ),
                    new_trade_vol_sigma2: Dual::<9>::variable(spec.new_trade_vol_sigma2, 3),
                    exposure_correlation_rho12: Dual::<9>::variable(
                        spec.exposure_correlation_rho12,
                        4,
                    ),
                    counterparty_hazard_rate_lambda: Dual::<9>::variable(
                        spec.counterparty_hazard_rate_lambda,
                        5,
                    ),
                    counterparty_recovery_rate_r: Dual::<9>::variable(
                        spec.counterparty_recovery_rate_r,
                        6,
                    ),
                    ois_discount_rate_r: Dual::<9>::variable(spec.ois_discount_rate_r, 7),
                    maturity_years_t: Dual::<9>::variable(spec.maturity_years_t, 8),
                },
            );
            for (idx, &ad) in r_ad.incremental_trade_cva.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.existing_portfolio_mtm_m1 += h;
                        sm.existing_portfolio_mtm_m1 -= h;
                    }
                    1 => {
                        sp.new_trade_mtm_m2 += h;
                        sm.new_trade_mtm_m2 -= h;
                    }
                    2 => {
                        sp.existing_portfolio_vol_sigma1 += h;
                        sm.existing_portfolio_vol_sigma1 -= h;
                    }
                    3 => {
                        sp.new_trade_vol_sigma2 += h;
                        sm.new_trade_vol_sigma2 -= h;
                    }
                    4 => {
                        sp.exposure_correlation_rho12 += h;
                        sm.exposure_correlation_rho12 -= h;
                    }
                    5 => {
                        sp.counterparty_hazard_rate_lambda += h;
                        sm.counterparty_hazard_rate_lambda -= h;
                    }
                    6 => {
                        sp.counterparty_recovery_rate_r += h;
                        sm.counterparty_recovery_rate_r -= h;
                    }
                    7 => {
                        sp.ois_discount_rate_r += h;
                        sm.ois_discount_rate_r -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_gregory_incremental_marginal_netting_set_cva(&sp)
                    .incremental_trade_cva
                    - evaluate_gregory_incremental_marginal_netting_set_cva(&sm)
                        .incremental_trade_cva)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_gregory_incremental_marginal_netting_set_cva_module::*;

mod evaluate_burgard_kjaer_semi_replication_fva_cva_dva_module {
    // Burgard & Kjaer (2011, Journal of Derivatives 18(4), "Partial Differential Equation
    // Representations of Derivatives with Bilateral Counterparty Risk and Funding Costs";
    // 2013, Risk, "Funding Strategies, Funding Costs") & Brigo, Morini & Pallavicini (2013, Ch. 16)
    // Bilateral CVA/DVA/FCA/FBA Semi-Replication & DVA-FBA Double-Counting Invariance Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BurgardKjaerSemiReplicationXvaSpec<T: Real = f64> {
        pub expected_positive_exposure_epe: T,
        pub expected_negative_exposure_ene: T,
        pub counterparty_hazard_rate_lambda_c: T,
        pub counterparty_lgd_c: T,
        pub bank_hazard_rate_lambda_b: T,
        pub bank_lgd_b: T,
        pub bank_funding_shortfall_fraction_gamma_f: T,
        pub ois_discount_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BurgardKjaerSemiReplicationXvaReport<T: Real = f64> {
        pub joint_survival_risky_annuity: T,
        pub bilateral_cva_charge: T,
        pub bilateral_dva_benefit: T,
        pub funding_cost_adjustment_fca: T,
        pub non_overlapping_funding_benefit_fba: T,
        pub dva_fba_double_counting_overlap: T,
        pub semi_replication_total_xva_charge: T,
    }

    pub fn evaluate_burgard_kjaer_semi_replication_fva_cva_dva<T: Real>(
        spec: &BurgardKjaerSemiReplicationXvaSpec<T>,
    ) -> BurgardKjaerSemiReplicationXvaReport<T> {
        let epe = spec
            .expected_positive_exposure_epe
            .max_val(T::from_f64(0.0));
        let ene = spec
            .expected_negative_exposure_ene
            .max_val(T::from_f64(0.0));
        let lam_c = spec
            .counterparty_hazard_rate_lambda_c
            .max_val(T::from_f64(1e-4));
        let lgd_c = spec
            .counterparty_lgd_c
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(1.0));
        let lam_b = spec.bank_hazard_rate_lambda_b.max_val(T::from_f64(1e-4));
        let lgd_b = spec
            .bank_lgd_b
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(1.0));
        let gamma_f = spec
            .bank_funding_shortfall_fraction_gamma_f
            .max_val(T::from_f64(0.0))
            .min_val(T::from_f64(1.0));
        let r = spec.ois_discount_rate_r.max_val(T::from_f64(0.001));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

        let total_disc = r + lam_c + lam_b;
        let annuity = (T::one() - (-total_disc * t).exp()) / total_disc;

        let cva = lam_c * lgd_c * epe * annuity;
        let dva = lam_b * lgd_b * ene * annuity;

        // Under Burgard-Kjaer (2011, 2013), bank unsecured funding spread is s_F = gamma_f * lambda_B * LGD_B:
        let s_f = gamma_f * lam_b * lgd_b;
        let fca = s_f * epe * annuity;
        let fba_non_overlap = (T::one() - gamma_f) * lam_b * lgd_b * ene * annuity;
        let overlap = gamma_f * lam_b * lgd_b * ene * annuity;

        let total_xva = cva - dva + fca;

        BurgardKjaerSemiReplicationXvaReport {
            joint_survival_risky_annuity: annuity,
            bilateral_cva_charge: cva,
            bilateral_dva_benefit: dva,
            funding_cost_adjustment_fca: fca,
            non_overlapping_funding_benefit_fba: fba_non_overlap,
            dva_fba_double_counting_overlap: overlap,
            semi_replication_total_xva_charge: total_xva,
        }
    }

    #[cfg(test)]
    mod evaluate_burgard_kjaer_semi_replication_fva_cva_dva_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn burgard_kjaer_semi_replication_fva_cva_dva_invariants_and_ad_parity() {
            let spec = BurgardKjaerSemiReplicationXvaSpec {
                expected_positive_exposure_epe: 25.0,
                expected_negative_exposure_ene: 14.0,
                counterparty_hazard_rate_lambda_c: 0.022,
                counterparty_lgd_c: 0.60,
                bank_hazard_rate_lambda_b: 0.012,
                bank_lgd_b: 0.60,
                bank_funding_shortfall_fraction_gamma_f: 0.65,
                ois_discount_rate_r: 0.038,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_burgard_kjaer_semi_replication_fva_cva_dva(&spec);
            assert!(
                (rep.non_overlapping_funding_benefit_fba + rep.dva_fba_double_counting_overlap
                    - rep.bilateral_dva_benefit)
                    .abs()
                    < 1e-12
            );
            assert!(rep.semi_replication_total_xva_charge > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_burgard_kjaer_semi_replication_fva_cva_dva(
                &BurgardKjaerSemiReplicationXvaSpec {
                    expected_positive_exposure_epe: Dual::<6>::variable(
                        spec.expected_positive_exposure_epe,
                        0,
                    ),
                    expected_negative_exposure_ene: Dual::<6>::variable(
                        spec.expected_negative_exposure_ene,
                        1,
                    ),
                    counterparty_hazard_rate_lambda_c: Dual::<6>::variable(
                        spec.counterparty_hazard_rate_lambda_c,
                        2,
                    ),
                    counterparty_lgd_c: Dual::from_f64(spec.counterparty_lgd_c),
                    bank_hazard_rate_lambda_b: Dual::<6>::variable(
                        spec.bank_hazard_rate_lambda_b,
                        3,
                    ),
                    bank_lgd_b: Dual::from_f64(spec.bank_lgd_b),
                    bank_funding_shortfall_fraction_gamma_f: Dual::<6>::variable(
                        spec.bank_funding_shortfall_fraction_gamma_f,
                        4,
                    ),
                    ois_discount_rate_r: Dual::<6>::variable(spec.ois_discount_rate_r, 5),
                    maturity_years_t: Dual::from_f64(spec.maturity_years_t),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.expected_positive_exposure_epe,
                    spec.expected_negative_exposure_ene,
                    spec.counterparty_hazard_rate_lambda_c,
                    spec.bank_hazard_rate_lambda_b,
                    spec.bank_funding_shortfall_fraction_gamma_f,
                    spec.ois_discount_rate_r,
                ];
                v[i] += h;
                let up = evaluate_burgard_kjaer_semi_replication_fva_cva_dva(
                    &BurgardKjaerSemiReplicationXvaSpec {
                        expected_positive_exposure_epe: v[0],
                        expected_negative_exposure_ene: v[1],
                        counterparty_hazard_rate_lambda_c: v[2],
                        bank_hazard_rate_lambda_b: v[3],
                        bank_funding_shortfall_fraction_gamma_f: v[4],
                        ois_discount_rate_r: v[5],
                        ..spec
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_burgard_kjaer_semi_replication_fva_cva_dva(
                    &BurgardKjaerSemiReplicationXvaSpec {
                        expected_positive_exposure_epe: v[0],
                        expected_negative_exposure_ene: v[1],
                        counterparty_hazard_rate_lambda_c: v[2],
                        bank_hazard_rate_lambda_b: v[3],
                        bank_funding_shortfall_fraction_gamma_f: v[4],
                        ois_discount_rate_r: v[5],
                        ..spec
                    },
                );
                let fd = (up.semi_replication_total_xva_charge
                    - dn.semi_replication_total_xva_charge)
                    / (2.0 * h);
                assert!(
                    (r_ad.semi_replication_total_xva_charge.d[i] - fd).abs() < 1e-5,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_burgard_kjaer_semi_replication_fva_cva_dva_module::*;

mod evaluate_hull_white_exponential_hazard_wrong_way_risk_cva_module {
    // Hull & White (2012, Journal of Derivatives 19(3), "CVA and Wrong-Way Risk"),
    // Gregory (2020, The xVA Challenge, 4th ed., Ch. 17), & Brigo, Morini & Pallavicini (2013, Ch. 11)
    // Exponential Hazard-Exposure Wrong-Way Risk (WWR) Girsanov Drift-Shifted CVA Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct HullWhiteExponentialWwrCvaSpec<T: Real = f64> {
        pub portfolio_expected_mtm_mu_v: T,
        pub portfolio_exposure_vol_sigma_v: T,
        pub baseline_cds_hazard_rate_lambda0: T,
        pub wwr_exponential_sensitivity_b: T,
        pub counterparty_lgd: T,
        pub ois_discount_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct HullWhiteExponentialWwrCvaReport<T: Real = f64> {
        pub calibration_intercept_a_t: T,
        pub independent_expected_positive_exposure: T,
        pub wwr_girsanov_shifted_epe: T,
        pub independent_baseline_cva: T,
        pub wrong_way_risk_adjusted_cva: T,
        pub wwr_cva_multiplier_ratio: T,
    }

    pub fn evaluate_hull_white_exponential_hazard_wrong_way_risk_cva<T: Real>(
        spec: &HullWhiteExponentialWwrCvaSpec<T>,
    ) -> HullWhiteExponentialWwrCvaReport<T> {
        let mu_v = spec.portfolio_expected_mtm_mu_v;
        let sig_v = spec
            .portfolio_exposure_vol_sigma_v
            .max_val(T::from_f64(0.1));
        let lam0 = spec
            .baseline_cds_hazard_rate_lambda0
            .max_val(T::from_f64(1e-4));
        let b_wwr = spec.wwr_exponential_sensitivity_b;
        let lgd = spec
            .counterparty_lgd
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(1.0));
        let r = spec.ois_discount_rate_r.max_val(T::from_f64(0.001));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.10));

        let var_t = sig_v * sig_v * t;
        // Hull-White (2012, Eq. 2) hazard calibration a(T) so E[exp(a(T) + b * V_T)] = lambda0:
        let a_t = lam0.ln() - b_wwr * mu_v - T::from_f64(0.5) * b_wwr * b_wwr * var_t;

        let dt = t / T::from_f64(4.0);
        let eff_disc = r + lam0;
        let mut cva_ind = T::zero();
        let mut cva_wwr = T::zero();
        let mut weight_sum = T::zero();

        for k in 0..4 {
            let uk = dt * (T::from_f64(k as f64) + T::from_f64(0.5));
            let var_u = sig_v * sig_v * uk;
            let std_u = var_u.sqrt();

            let z_ind = mu_v / std_u;
            let epe_ind_u = mu_v * norm_cdf(z_ind) + std_u * norm_pdf(z_ind);

            // Under exponential tilting exp(a(u) + b * V_u), Girsanov shifts mean of V_u to mu_v + b * sigma_v^2 * u:
            let mu_wwr_u = mu_v + b_wwr * var_u;
            let z_wwr = mu_wwr_u / std_u;
            let epe_wwr_u = mu_wwr_u * norm_cdf(z_wwr) + std_u * norm_pdf(z_wwr);

            let df_surv = dt * (-eff_disc * uk).exp();
            weight_sum += df_surv;
            cva_ind += lam0 * lgd * df_surv * epe_ind_u;
            cva_wwr += lam0 * lgd * df_surv * epe_wwr_u;
        }

        let epe_ind = cva_ind / (lam0 * lgd * weight_sum);
        let epe_wwr = cva_wwr / (lam0 * lgd * weight_sum);
        let mult = cva_wwr / cva_ind.max_val(T::from_f64(1e-6));

        HullWhiteExponentialWwrCvaReport {
            calibration_intercept_a_t: a_t,
            independent_expected_positive_exposure: epe_ind,
            wwr_girsanov_shifted_epe: epe_wwr,
            independent_baseline_cva: cva_ind,
            wrong_way_risk_adjusted_cva: cva_wwr,
            wwr_cva_multiplier_ratio: mult,
        }
    }

    #[cfg(test)]
    mod evaluate_hull_white_exponential_hazard_wrong_way_risk_cva_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hull_white_exponential_wwr_cva_invariants_and_ad_parity() {
            let spec = HullWhiteExponentialWwrCvaSpec {
                portfolio_expected_mtm_mu_v: 8.0,
                portfolio_exposure_vol_sigma_v: 12.0,
                baseline_cds_hazard_rate_lambda0: 0.025,
                wwr_exponential_sensitivity_b: 0.035,
                counterparty_lgd: 0.60,
                ois_discount_rate_r: 0.038,
                maturity_years_t: 4.0,
            };
            let rep = evaluate_hull_white_exponential_hazard_wrong_way_risk_cva(&spec);
            assert!(rep.wrong_way_risk_adjusted_cva > rep.independent_baseline_cva);
            assert!(rep.wwr_cva_multiplier_ratio > 1.0);

            let h = 1e-6;
            let r_ad = evaluate_hull_white_exponential_hazard_wrong_way_risk_cva(
                &HullWhiteExponentialWwrCvaSpec {
                    portfolio_expected_mtm_mu_v: Dual::<7>::variable(
                        spec.portfolio_expected_mtm_mu_v,
                        0,
                    ),
                    portfolio_exposure_vol_sigma_v: Dual::<7>::variable(
                        spec.portfolio_exposure_vol_sigma_v,
                        1,
                    ),
                    baseline_cds_hazard_rate_lambda0: Dual::<7>::variable(
                        spec.baseline_cds_hazard_rate_lambda0,
                        2,
                    ),
                    wwr_exponential_sensitivity_b: Dual::<7>::variable(
                        spec.wwr_exponential_sensitivity_b,
                        3,
                    ),
                    counterparty_lgd: Dual::<7>::variable(spec.counterparty_lgd, 4),
                    ois_discount_rate_r: Dual::<7>::variable(spec.ois_discount_rate_r, 5),
                    maturity_years_t: Dual::<7>::variable(spec.maturity_years_t, 6),
                },
            );
            for i in 0..7 {
                let mut v = [
                    spec.portfolio_expected_mtm_mu_v,
                    spec.portfolio_exposure_vol_sigma_v,
                    spec.baseline_cds_hazard_rate_lambda0,
                    spec.wwr_exponential_sensitivity_b,
                    spec.counterparty_lgd,
                    spec.ois_discount_rate_r,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_hull_white_exponential_hazard_wrong_way_risk_cva(
                    &HullWhiteExponentialWwrCvaSpec {
                        portfolio_expected_mtm_mu_v: v[0],
                        portfolio_exposure_vol_sigma_v: v[1],
                        baseline_cds_hazard_rate_lambda0: v[2],
                        wwr_exponential_sensitivity_b: v[3],
                        counterparty_lgd: v[4],
                        ois_discount_rate_r: v[5],
                        maturity_years_t: v[6],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_hull_white_exponential_hazard_wrong_way_risk_cva(
                    &HullWhiteExponentialWwrCvaSpec {
                        portfolio_expected_mtm_mu_v: v[0],
                        portfolio_exposure_vol_sigma_v: v[1],
                        baseline_cds_hazard_rate_lambda0: v[2],
                        wwr_exponential_sensitivity_b: v[3],
                        counterparty_lgd: v[4],
                        ois_discount_rate_r: v[5],
                        maturity_years_t: v[6],
                    },
                );
                let fd =
                    (up.wrong_way_risk_adjusted_cva - dn.wrong_way_risk_adjusted_cva) / (2.0 * h);
                assert!(
                    (r_ad.wrong_way_risk_adjusted_cva.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_hull_white_exponential_hazard_wrong_way_risk_cva_module::*;

mod evaluate_gregory_albanese_bilateral_cva_dva_first_to_default_survival_module {
    // Gregory (2020, The xVA Challenge, 4th ed., Ch. 14, Eq. 14.7),
    // Brigo, Morini & Pallavicini (2013, Counterparty Credit Risk, Collateral and Funding, Ch. 4 & 12),
    // & Albanese & Andersen (2014) Bilateral CVA/DVA First-to-Default Survival Overlap Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryAlbaneseBilateralFtdCvaDvaSpec<T: Real = f64> {
        pub expected_positive_exposure_epe: T,
        pub expected_negative_exposure_ene: T,
        pub counterparty_hazard_rate_lambda_c: T,
        pub own_bank_hazard_rate_lambda_b: T,
        pub counterparty_lgd_c: T,
        pub own_bank_lgd_b: T,
        pub ois_discount_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GregoryAlbaneseBilateralFtdCvaDvaReport<T: Real = f64> {
        pub unilateral_cva_charge: T,
        pub unilateral_dva_benefit: T,
        pub bilateral_ftd_cva_charge: T,
        pub bilateral_ftd_dva_benefit: T,
        pub first_to_default_survival_overlap_reduction: T,
        pub net_bilateral_bcva: T,
    }

    pub fn evaluate_gregory_albanese_bilateral_cva_dva_first_to_default_survival<T: Real>(
        spec: &GregoryAlbaneseBilateralFtdCvaDvaSpec<T>,
    ) -> GregoryAlbaneseBilateralFtdCvaDvaReport<T> {
        let epe0 = spec
            .expected_positive_exposure_epe
            .max_val(T::from_f64(0.01));
        let ene0 = spec
            .expected_negative_exposure_ene
            .max_val(T::from_f64(0.01));
        let lam_c = spec
            .counterparty_hazard_rate_lambda_c
            .max_val(T::from_f64(1e-4));
        let lam_b = spec
            .own_bank_hazard_rate_lambda_b
            .max_val(T::from_f64(1e-4));
        let lgd_c = spec
            .counterparty_lgd_c
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(1.0));
        let lgd_b = spec
            .own_bank_lgd_b
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(1.0));
        let r = spec.ois_discount_rate_r.max_val(T::from_f64(0.001));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

        // 4-point midpoint quadrature over [0, T] with diffusive exposure scaling sqrt(u):
        let dt = t / T::from_f64(4.0);
        let mut cva_uni = T::zero();
        let mut dva_uni = T::zero();
        let mut cva_ftd = T::zero();
        let mut dva_ftd = T::zero();

        for k in 0..4 {
            let uk = dt * (T::from_f64(k as f64) + T::from_f64(0.5));
            let profile_scale = uk.sqrt();
            let epe_u = epe0 * profile_scale;
            let ene_u = ene0 * profile_scale;

            let df_r = (-r * uk).exp();
            let surv_c = (-lam_c * uk).exp();
            let surv_b = (-lam_b * uk).exp();
            let surv_joint = surv_c * surv_b;

            cva_uni += dt * df_r * lgd_c * lam_c * surv_c * epe_u;
            dva_uni += dt * df_r * lgd_b * lam_b * surv_b * ene_u;

            // Gregory (2020, Eq. 14.7): FtD CVA conditions on bank survival S_B(u),
            // and FtD DVA conditions on counterparty survival S_C(u):
            cva_ftd += dt * df_r * lgd_c * lam_c * surv_joint * epe_u;
            dva_ftd += dt * df_r * lgd_b * lam_b * surv_joint * ene_u;
        }

        let overlap_red = (cva_uni - cva_ftd) + (dva_uni - dva_ftd);
        let net_bcva = cva_ftd - dva_ftd;

        GregoryAlbaneseBilateralFtdCvaDvaReport {
            unilateral_cva_charge: cva_uni,
            unilateral_dva_benefit: dva_uni,
            bilateral_ftd_cva_charge: cva_ftd,
            bilateral_ftd_dva_benefit: dva_ftd,
            first_to_default_survival_overlap_reduction: overlap_red,
            net_bilateral_bcva: net_bcva,
        }
    }

    #[cfg(test)]
    mod evaluate_gregory_albanese_bilateral_cva_dva_first_to_default_survival_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gregory_albanese_bilateral_ftd_cva_dva_invariants_and_ad_parity() {
            let spec = GregoryAlbaneseBilateralFtdCvaDvaSpec {
                expected_positive_exposure_epe: 15.0,
                expected_negative_exposure_ene: 10.0,
                counterparty_hazard_rate_lambda_c: 0.030,
                own_bank_hazard_rate_lambda_b: 0.018,
                counterparty_lgd_c: 0.60,
                own_bank_lgd_b: 0.60,
                ois_discount_rate_r: 0.038,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_gregory_albanese_bilateral_cva_dva_first_to_default_survival(&spec);
            assert!(rep.bilateral_ftd_cva_charge < rep.unilateral_cva_charge);
            assert!(rep.bilateral_ftd_dva_benefit < rep.unilateral_dva_benefit);
            assert!(rep.first_to_default_survival_overlap_reduction > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_gregory_albanese_bilateral_cva_dva_first_to_default_survival(
                &GregoryAlbaneseBilateralFtdCvaDvaSpec {
                    expected_positive_exposure_epe: Dual::<8>::variable(
                        spec.expected_positive_exposure_epe,
                        0,
                    ),
                    expected_negative_exposure_ene: Dual::<8>::variable(
                        spec.expected_negative_exposure_ene,
                        1,
                    ),
                    counterparty_hazard_rate_lambda_c: Dual::<8>::variable(
                        spec.counterparty_hazard_rate_lambda_c,
                        2,
                    ),
                    own_bank_hazard_rate_lambda_b: Dual::<8>::variable(
                        spec.own_bank_hazard_rate_lambda_b,
                        3,
                    ),
                    counterparty_lgd_c: Dual::<8>::variable(spec.counterparty_lgd_c, 4),
                    own_bank_lgd_b: Dual::<8>::variable(spec.own_bank_lgd_b, 5),
                    ois_discount_rate_r: Dual::<8>::variable(spec.ois_discount_rate_r, 6),
                    maturity_years_t: Dual::<8>::variable(spec.maturity_years_t, 7),
                },
            );
            for i in 0..8 {
                let mut v = [
                    spec.expected_positive_exposure_epe,
                    spec.expected_negative_exposure_ene,
                    spec.counterparty_hazard_rate_lambda_c,
                    spec.own_bank_hazard_rate_lambda_b,
                    spec.counterparty_lgd_c,
                    spec.own_bank_lgd_b,
                    spec.ois_discount_rate_r,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_gregory_albanese_bilateral_cva_dva_first_to_default_survival(
                    &GregoryAlbaneseBilateralFtdCvaDvaSpec {
                        expected_positive_exposure_epe: v[0],
                        expected_negative_exposure_ene: v[1],
                        counterparty_hazard_rate_lambda_c: v[2],
                        own_bank_hazard_rate_lambda_b: v[3],
                        counterparty_lgd_c: v[4],
                        own_bank_lgd_b: v[5],
                        ois_discount_rate_r: v[6],
                        maturity_years_t: v[7],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_gregory_albanese_bilateral_cva_dva_first_to_default_survival(
                    &GregoryAlbaneseBilateralFtdCvaDvaSpec {
                        expected_positive_exposure_epe: v[0],
                        expected_negative_exposure_ene: v[1],
                        counterparty_hazard_rate_lambda_c: v[2],
                        own_bank_hazard_rate_lambda_b: v[3],
                        counterparty_lgd_c: v[4],
                        own_bank_lgd_b: v[5],
                        ois_discount_rate_r: v[6],
                        maturity_years_t: v[7],
                    },
                );
                let fd = (up.net_bilateral_bcva - dn.net_bilateral_bcva) / (2.0 * h);
                assert!(
                    (r_ad.net_bilateral_bcva.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_gregory_albanese_bilateral_cva_dva_first_to_default_survival_module::*;

use ql_ad::Real;

/// Specification for the Albanese–Caenazzo–Crépey (2017) / Green (2015) Basel III SA-CCR
/// Exposure-at-Default and Capital Valuation Adjustment (KVA) engine.
///
/// Reference: Albanese, Caenazzo & Crépey (2017, *Risk*, "Credit, Funding, Margin, and Capital
/// Valuation Adjustments"), Green (2015, *XVA*, Ch. 16), Gregory (2020, Ch. 18).
#[derive(Debug, Clone, Copy)]
pub struct AlbaneseCrepeySaCcrKvaSpec<T: Real> {
    /// Trade notional $N_{\text{notional}}$.
    pub trade_notional: T,
    /// Current mark-to-market value $V_0$.
    pub current_mtm_value_v0: T,
    /// Posted variation/initial collateral balance $C_0$.
    pub collateral_balance_c0: T,
    /// Supervisory risk weight $w_{\text{CCR}}$ for counterparty default risk capital.
    pub counterparty_risk_weight_w: T,
    /// Shareholder capital hurdle rate $h_{\text{KVA}}$.
    pub capital_hurdle_rate_h: T,
    /// OIS risk-free discount rate $r$.
    pub ois_discount_rate_r: T,
    /// Joint survival hazard rate $\lambda_{\text{joint}} = \lambda_C + \lambda_B$.
    pub joint_hazard_rate_lambda: T,
    /// Portfolio maturity $T$ (in years).
    pub maturity_years_t: T,
}

/// Analytics returned by [`evaluate_albanese_crepey_sa_ccr_kadd_capital_xva`].
#[derive(Debug, Clone, Copy)]
pub struct AlbaneseCrepeySaCcrKvaAnalytics<T: Real> {
    /// Initial SA-CCR Replacement Cost $RC_0$ plus Potential Future Exposure $PFE_0$.
    pub initial_sa_ccr_exposure_at_default: T,
    /// Initial SA-CCR over-collateralization multiplier $m_0 \in [0.05, 1]$.
    pub sa_ccr_pfe_multiplier: T,
    /// Initial Pillar 1 CCR regulatory capital requirement $K_{\text{reg}}(0)$.
    pub initial_regulatory_capital_k0: T,
    /// Lifetime discounted Capital Valuation Adjustment $KVA$.
    pub lifetime_sa_ccr_kva: T,
}

/// Evaluates Basel III SA-CCR Exposure-at-Default ($EAD = 1.4 \cdot (RC + m \cdot \text{AddOn})$)
/// and integrates Capital Valuation Adjustment ($KVA$) across a 4-point midpoint grid over $[0, T]$.
pub fn evaluate_albanese_crepey_sa_ccr_kadd_capital_xva<T: Real>(
    spec: &AlbaneseCrepeySaCcrKvaSpec<T>,
) -> AlbaneseCrepeySaCcrKvaAnalytics<T> {
    let alpha_saccr = T::from_f64(1.4);
    let sf_ir = T::from_f64(0.005);
    let floor_f = T::from_f64(0.05);
    let tier1_ratio = T::from_f64(0.105);
    let t_mat = spec.maturity_years_t.max_val(T::from_f64(1e-4));

    let smooth_pos = |x: T| -> T {
        let eps = T::from_f64(0.05);
        T::from_f64(0.5) * (x + (x * x + eps * eps).sqrt())
    };

    let eval_saccr_at = |rem_tau: T, scale_v: T| -> (T, T, T) {
        // Supervisory duration SD(tau) = (1 - exp(-0.05 * tau)) / 0.05
        let sd = (T::one() - (-T::from_f64(0.05) * rem_tau).exp()) / T::from_f64(0.05);
        let addon = (sf_ir * spec.trade_notional * sd).max_val(T::from_f64(1e-6));
        let net_v = spec.current_mtm_value_v0 * scale_v - spec.collateral_balance_c0;
        let rc = smooth_pos(net_v);
        // SA-CCR multiplier m = min(1, f + (1 - f) * exp(net_v / (2*(1-f)*AddOn)))
        // Smooth saturation below 1.0 when net_v < 0
        let neg_net_v = -smooth_pos(-net_v);
        let mult = floor_f
            + (T::one() - floor_f)
                * (neg_net_v / (T::from_f64(2.0) * (T::one() - floor_f) * addon)).exp();
        let pfe = mult * addon;
        let ead = alpha_saccr * (rc + pfe);
        let k_reg = ead * spec.counterparty_risk_weight_w * tier1_ratio;
        (ead, mult, k_reg)
    };

    let (ead0, mult0, k0) = eval_saccr_at(t_mat, T::one());

    const STEPS: usize = 4;
    let dt = t_mat / T::from_f64(STEPS as f64);
    let disc_rate = spec.ois_discount_rate_r + spec.joint_hazard_rate_lambda;
    let mut kva = T::zero();

    for m in 0..STEPS {
        let t_mid = (T::from_f64(m as f64) + T::from_f64(0.5)) * dt;
        let rem_tau = (t_mat - t_mid).max_val(T::from_f64(1e-4));
        let scale_v = (rem_tau / t_mat).sqrt();
        let (_ead_t, _mult_t, k_reg_t) = eval_saccr_at(rem_tau, scale_v);
        let df = (-disc_rate * t_mid).exp();
        kva += dt * df * spec.capital_hurdle_rate_h * k_reg_t;
    }

    AlbaneseCrepeySaCcrKvaAnalytics {
        initial_sa_ccr_exposure_at_default: ead0,
        sa_ccr_pfe_multiplier: mult0,
        initial_regulatory_capital_k0: k0,
        lifetime_sa_ccr_kva: kva,
    }
}

#[cfg(test)]
mod albanese_crepey_saccr_kva_suite64_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn albanese_crepey_sa_ccr_kva_invariants_and_ad_parity() {
        let base = AlbaneseCrepeySaCcrKvaSpec {
            trade_notional: 1_000.0_f64,
            current_mtm_value_v0: 8.0,
            collateral_balance_c0: 12.0,
            counterparty_risk_weight_w: 0.50,
            capital_hurdle_rate_h: 0.10,
            ois_discount_rate_r: 0.038,
            joint_hazard_rate_lambda: 0.025,
            maturity_years_t: 5.0,
        };
        let r0 = evaluate_albanese_crepey_sa_ccr_kadd_capital_xva(&base);
        assert!(r0.initial_sa_ccr_exposure_at_default > 0.0);
        assert!(r0.sa_ccr_pfe_multiplier > 0.05 && r0.sa_ccr_pfe_multiplier < 1.0);
        assert!(r0.lifetime_sa_ccr_kva > 0.0);

        let inputs = [
            base.trade_notional,
            base.current_mtm_value_v0,
            base.collateral_balance_c0,
            base.counterparty_risk_weight_w,
            base.capital_hurdle_rate_h,
            base.ois_discount_rate_r,
            base.joint_hazard_rate_lambda,
            base.maturity_years_t,
        ];
        let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            dual_in[i] = Dual::variable(inputs[i], i);
        }
        let spec_ad = AlbaneseCrepeySaCcrKvaSpec {
            trade_notional: dual_in[0],
            current_mtm_value_v0: dual_in[1],
            collateral_balance_c0: dual_in[2],
            counterparty_risk_weight_w: dual_in[3],
            capital_hurdle_rate_h: dual_in[4],
            ois_discount_rate_r: dual_in[5],
            joint_hazard_rate_lambda: dual_in[6],
            maturity_years_t: dual_in[7],
        };
        let r_ad = evaluate_albanese_crepey_sa_ccr_kadd_capital_xva(&spec_ad);
        let eps = 1e-6;
        for i in 0..8 {
            let mut p = inputs;
            let mut m = inputs;
            p[i] += eps;
            m[i] -= eps;
            let mk = |v: [f64; 8]| AlbaneseCrepeySaCcrKvaSpec {
                trade_notional: v[0],
                current_mtm_value_v0: v[1],
                collateral_balance_c0: v[2],
                counterparty_risk_weight_w: v[3],
                capital_hurdle_rate_h: v[4],
                ois_discount_rate_r: v[5],
                joint_hazard_rate_lambda: v[6],
                maturity_years_t: v[7],
            };
            let rp = evaluate_albanese_crepey_sa_ccr_kadd_capital_xva(&mk(p));
            let rm = evaluate_albanese_crepey_sa_ccr_kadd_capital_xva(&mk(m));
            let fd = (rp.lifetime_sa_ccr_kva - rm.lifetime_sa_ccr_kva) / (2.0 * eps);
            assert!(
                r_ad.lifetime_sa_ccr_kva.d[i].abs() > 1e-7,
                "zero derivative at {i}"
            );
            assert!((r_ad.lifetime_sa_ccr_kva.d[i] - fd).abs() < 1e-4);
        }
    }
}

pub mod brigo_pallavicini_ccp_xva_suite65_mod {
    use ql_ad::Real;
    use ql_math::norm_cdf;

    /// Specification for the Brigo–Pallavicini (2014) / Gregory (2020) centrally cleared CCP
    /// Initial Margin MVA, MPoR tail-breach CVA, and Default Fund mutualization engine.
    ///
    /// Reference: Brigo & Pallavicini (2014, *IJTAF* 17(1), "Nonlinear Consistent Valuation of CCP
    /// Cleared or CSA Bilateral Trades"), Crépey (2015, *Math. Finance*), Gregory (2020, Ch. 17).
    #[derive(Debug, Clone, Copy)]
    pub struct BrigoPallaviciniCcpXvaSpec<T: Real> {
        /// Portfolio mark-to-market volatility $\sigma_V$ (in currency units per $\sqrt{\text{year}}$).
        pub portfolio_mtm_volatility_sigma_v: T,
        /// Margin Period of Risk $\Delta_{\text{MPoR}}$ (in years, e.g., $10/252$).
        pub margin_period_of_risk_years: T,
        /// Initial Margin confidence level $\alpha_{\text{IM}} \in (0.90, 0.999)$ (e.g., $0.99$).
        pub im_confidence_alpha: T,
        /// Initial Margin funding spread $s_{\text{IM}}$.
        pub im_funding_spread_s_im: T,
        /// Clearing member default hazard rate $\lambda_{\text{CM}}$ (for residual tail breach CVA).
        pub clearing_member_hazard_lambda_cm: T,
        /// Default Fund to Initial Margin scaling ratio $\theta_{\text{DF}}$.
        pub default_fund_scaling_theta_df: T,
        /// Default Fund capital + fellow-member mutualization hazard spread $s_{\text{DF}}$.
        pub default_fund_mutualization_spread_s_df: T,
        /// OIS risk-free discount rate $r$.
        pub ois_discount_rate_r: T,
        /// Portfolio maturity $T$ (in years).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_brigo_pallavicini_ccp_default_fund_cva_mva`].
    #[derive(Debug, Clone, Copy)]
    pub struct BrigoPallaviciniCcpXvaAnalytics<T: Real> {
        /// Lifetime Initial Margin Valuation Adjustment ($MVA_{\text{CCP}}$).
        pub initial_margin_mva: T,
        /// Residual MPoR tail-exceedance CVA above Initial Margin ($CVA_{\text{CCP}}$).
        pub residual_mpor_breach_cva: T,
        /// CCP Default Fund capital and loss-mutualization valuation adjustment ($DFVA$).
        pub default_fund_mutualization_cost: T,
        /// Total CCP clearing valuation adjustment $XVA_{\text{CCP}} = MVA + CVA_{\text{CCP}} + DFVA$.
        pub total_ccp_clearing_xva: T,
    }

    /// Evaluates the Brigo–Pallavicini (2014) centrally cleared CCP valuation adjustments ($MVA$,
    /// residual tail-breach $CVA$, and Default Fund $DFVA$) via 4-point midpoint quadrature.
    pub fn evaluate_brigo_pallavicini_ccp_default_fund_cva_mva<T: Real>(
        spec: &BrigoPallaviciniCcpXvaSpec<T>,
    ) -> BrigoPallaviciniCcpXvaAnalytics<T> {
        let t_mat = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let mpor = spec.margin_period_of_risk_years.max_val(T::from_f64(1e-5));
        let alpha = spec
            .im_confidence_alpha
            .max_val(T::from_f64(0.80))
            .min_val(T::from_f64(0.9999));

        // Exact AD-propagating inverse normal quantile q_IM = Phi^{-1}(alpha) via two Newton steps
        let mut q_im = T::from_f64(ql_math::inv_norm_cdf(alpha.val().clamp(1e-12, 1.0 - 1e-12)));
        for _ in 0..2 {
            q_im = q_im - (q_im.norm_cdf() - alpha) / q_im.norm_pdf().max_val(T::from_f64(1e-15));
        }

        let phi_q = q_im.norm_pdf();
        let tail_exceedance_mult =
            (phi_q - q_im * (T::one() - norm_cdf(q_im))).max_val(T::from_f64(1e-8));
        let mpor_vol_0 = spec.portfolio_mtm_volatility_sigma_v * mpor.sqrt();

        const STEPS: usize = 4;
        let dt = t_mat / T::from_f64(STEPS as f64);
        let disc_rate = spec.ois_discount_rate_r + spec.clearing_member_hazard_lambda_cm;
        let lgd_ccp = T::from_f64(0.60);

        let mut mva = T::zero();
        let mut cva_ccp = T::zero();
        let mut dfva = T::zero();

        for m in 0..STEPS {
            let t_mid = (T::from_f64(m as f64) + T::from_f64(0.5)) * dt;
            let rem_fraction = ((t_mat - t_mid) / t_mat).max_val(T::from_f64(1e-4)).sqrt();
            let sigma_mpor_t = mpor_vol_0 * rem_fraction;
            let im_t = q_im * sigma_mpor_t;
            let epe_breach_t = sigma_mpor_t * tail_exceedance_mult;
            let df_t = spec.default_fund_scaling_theta_df * im_t;

            let df_disc = (-disc_rate * t_mid).exp();
            mva += dt * df_disc * spec.im_funding_spread_s_im * im_t;
            cva_ccp +=
                dt * df_disc * lgd_ccp * spec.clearing_member_hazard_lambda_cm * epe_breach_t;
            dfva += dt * df_disc * spec.default_fund_mutualization_spread_s_df * df_t;
        }

        BrigoPallaviciniCcpXvaAnalytics {
            initial_margin_mva: mva,
            residual_mpor_breach_cva: cva_ccp,
            default_fund_mutualization_cost: dfva,
            total_ccp_clearing_xva: mva + cva_ccp + dfva,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn brigo_pallavicini_ccp_xva_invariants_and_ad_parity() {
            let base = BrigoPallaviciniCcpXvaSpec {
                portfolio_mtm_volatility_sigma_v: 45.0_f64,
                margin_period_of_risk_years: 10.0 / 252.0,
                im_confidence_alpha: 0.99,
                im_funding_spread_s_im: 0.0085,
                clearing_member_hazard_lambda_cm: 0.015,
                default_fund_scaling_theta_df: 0.18,
                default_fund_mutualization_spread_s_df: 0.025,
                ois_discount_rate_r: 0.038,
                maturity_years_t: 5.0,
            };
            let r0 = evaluate_brigo_pallavicini_ccp_default_fund_cva_mva(&base);
            assert!(r0.initial_margin_mva > 0.0);
            assert!(r0.residual_mpor_breach_cva > 0.0);
            assert!(r0.default_fund_mutualization_cost > 0.0);
            assert!(r0.total_ccp_clearing_xva > r0.initial_margin_mva);

            let inputs = [
                base.portfolio_mtm_volatility_sigma_v,
                base.margin_period_of_risk_years,
                base.im_confidence_alpha,
                base.im_funding_spread_s_im,
                base.clearing_member_hazard_lambda_cm,
                base.default_fund_scaling_theta_df,
                base.default_fund_mutualization_spread_s_df,
                base.ois_discount_rate_r,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<9>::from_f64(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BrigoPallaviciniCcpXvaSpec {
                portfolio_mtm_volatility_sigma_v: dual_in[0],
                margin_period_of_risk_years: dual_in[1],
                im_confidence_alpha: dual_in[2],
                im_funding_spread_s_im: dual_in[3],
                clearing_member_hazard_lambda_cm: dual_in[4],
                default_fund_scaling_theta_df: dual_in[5],
                default_fund_mutualization_spread_s_df: dual_in[6],
                ois_discount_rate_r: dual_in[7],
                maturity_years_t: dual_in[8],
            };
            let r_ad = evaluate_brigo_pallavicini_ccp_default_fund_cva_mva(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| BrigoPallaviciniCcpXvaSpec {
                    portfolio_mtm_volatility_sigma_v: v[0],
                    margin_period_of_risk_years: v[1],
                    im_confidence_alpha: v[2],
                    im_funding_spread_s_im: v[3],
                    clearing_member_hazard_lambda_cm: v[4],
                    default_fund_scaling_theta_df: v[5],
                    default_fund_mutualization_spread_s_df: v[6],
                    ois_discount_rate_r: v[7],
                    maturity_years_t: v[8],
                };
                let rp = evaluate_brigo_pallavicini_ccp_default_fund_cva_mva(&mk(p));
                let rm = evaluate_brigo_pallavicini_ccp_default_fund_cva_mva(&mk(m));
                let fd = (rp.total_ccp_clearing_xva - rm.total_ccp_clearing_xva) / (2.0 * eps);
                assert!(
                    r_ad.total_ccp_clearing_xva.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.total_ccp_clearing_xva.d[i] - fd).abs() < 1e-3);
            }
        }
    }
}
pub use brigo_pallavicini_ccp_xva_suite65_mod::*;

pub mod green_basel_leverage_ratio_kva_mod {
    use ql_ad::Real;

    /// Specification for the Green (2015, Ch. 13-15) / Albanese, Andersen & Iabichino (2015)
    /// dual-constraint RWA and Basel III Leverage Ratio KVA engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GreenBaselLeverageRatioKvaSpec<T: Real> {
        /// Initial SA-CCR exposure-at-default $\text{EAD}_0 > 0$.
        pub initial_sa_ccr_ead0: T,
        /// Gross trade notional $N_{\text{gross}} > 0$ subject to SA-CCR leverage add-on.
        pub gross_trade_notional: T,
        /// Regulatory SA-CCR PFE add-on factor $\psi_{\text{addon}} > 0$.
        pub sa_ccr_addon_factor_psi: T,
        /// Counterparty risk weight $w_{\text{cpty}} > 0$ under Standardized/IRB CCR.
        pub counterparty_risk_weight_w: T,
        /// Pillar 1 CET1 + buffers capital ratio requirement $k_{\text{cet1}}$ (e.g. 0.105).
        pub cet1_capital_ratio_k: T,
        /// Basel III Tier 1 Leverage Ratio requirement $\gamma_{\text{LR}}$ (e.g. 0.035).
        pub leverage_ratio_requirement_gamma: T,
        /// CVA capital multiplier $\beta_{\text{cva}} > 0$ per year of effective maturity.
        pub cva_capital_loading_beta: T,
        /// Shareholder return-on-equity hurdle rate $h_{\text{ROE}} > 0$.
        pub capital_hurdle_rate_h: T,
        /// Fraction $\phi_{\text{cap}} \in [0, 1]$ of capital funding benefit at OIS rate.
        pub capital_funding_offset_phi: T,
        /// Counterparty hazard rate $\lambda_C > 0$.
        pub counterparty_hazard_rate_lambda: T,
        /// Continuous OIS discount rate $r_{\text{ois}}$.
        pub ois_discount_rate_r: T,
        /// Trade maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Green Basel III Leverage Ratio & RWA KVA engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GreenBaselLeverageRatioKvaAnalytics<T: Real> {
        /// Stand-alone Pillar 1 RWA (CCR + CVA) KVA charge.
        pub rwa_only_kva_charge: T,
        /// Stand-alone Basel III Leverage Ratio KVA charge.
        pub leverage_ratio_only_kva_charge: T,
        /// Incremental surcharge $\text{KVA}_{\text{dual}} - \text{KVA}_{\text{RWA}}$ from LR binding.
        pub leverage_ratio_binding_surcharge: T,
        /// Total dual-constraint regulatory KVA charge.
        pub total_dual_constraint_kva: T,
    }

    /// Evaluates the Green (2015) / Albanese-Andersen (2015) Capital Valuation Adjustment (KVA)
    /// under the smooth binding maximum of RWA capital and Basel III Leverage Ratio capital.
    pub fn evaluate_green_basel_leverage_ratio_rwa_hurdle_kva<T: Real>(
        spec: &GreenBaselLeverageRatioKvaSpec<T>,
    ) -> GreenBaselLeverageRatioKvaAnalytics<T> {
        const PERIODS: usize = 4;
        let t_mat = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let dt = t_mat / T::from_f64(PERIODS as f64);
        let net_hurdle = (spec.capital_hurdle_rate_h
            - spec.capital_funding_offset_phi * spec.ois_discount_rate_r)
            .max_val(T::from_f64(1e-4));

        let mut rwa_kva = T::zero();
        let mut lr_kva = T::zero();
        let mut dual_kva = T::zero();
        // Use a capital-scale smoothing width so the leverage constraint remains
        // differentiable even when the RWA leg is materially dominant.
        let eps = T::from_f64(1.0);

        for m in 1..=PERIODS {
            let tm = (T::from_f64(m as f64) - T::from_f64(0.5)) * dt;
            let rem_frac = (t_mat - tm) / t_mat;
            let ead_m = spec.initial_sa_ccr_ead0 * rem_frac;

            let k_rwa = spec.cet1_capital_ratio_k
                * T::from_f64(12.5)
                * spec.counterparty_risk_weight_w
                * ead_m
                * (T::one() + spec.cva_capital_loading_beta * (t_mat - tm));

            let lr_exposure =
                ead_m + spec.sa_ccr_addon_factor_psi * spec.gross_trade_notional * rem_frac;
            let k_lr = spec.leverage_ratio_requirement_gamma * lr_exposure;

            // Smooth maximum of RWA capital and Leverage Ratio capital
            let k_bind = T::from_f64(0.5)
                * (k_rwa + k_lr + ((k_rwa - k_lr) * (k_rwa - k_lr) + eps * eps).sqrt());

            let df_surv = (-(spec.ois_discount_rate_r
                + spec.counterparty_hazard_rate_lambda
                + spec.capital_hurdle_rate_h)
                * tm)
                .exp();
            let weight = dt * df_surv * net_hurdle;

            rwa_kva += weight * k_rwa;
            lr_kva += weight * k_lr;
            dual_kva += weight * k_bind;
        }

        GreenBaselLeverageRatioKvaAnalytics {
            rwa_only_kva_charge: rwa_kva,
            leverage_ratio_only_kva_charge: lr_kva,
            leverage_ratio_binding_surcharge: dual_kva - rwa_kva,
            total_dual_constraint_kva: dual_kva,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn green_basel_leverage_ratio_kva_invariants_and_ad_parity() {
            let base = GreenBaselLeverageRatioKvaSpec {
                initial_sa_ccr_ead0: 18.0_f64,
                gross_trade_notional: 500.0,
                sa_ccr_addon_factor_psi: 0.04,
                counterparty_risk_weight_w: 0.50,
                cet1_capital_ratio_k: 0.105,
                leverage_ratio_requirement_gamma: 0.035,
                cva_capital_loading_beta: 0.08,
                capital_hurdle_rate_h: 0.11,
                capital_funding_offset_phi: 0.35,
                counterparty_hazard_rate_lambda: 0.015,
                ois_discount_rate_r: 0.035,
                maturity_years_t: 5.0,
            };
            let r0 = evaluate_green_basel_leverage_ratio_rwa_hurdle_kva(&base);
            assert!(r0.total_dual_constraint_kva > r0.rwa_only_kva_charge);
            assert!(r0.total_dual_constraint_kva > r0.leverage_ratio_only_kva_charge);
            assert!(r0.leverage_ratio_binding_surcharge > 0.0);

            let inputs = [
                base.initial_sa_ccr_ead0,
                base.gross_trade_notional,
                base.sa_ccr_addon_factor_psi,
                base.counterparty_risk_weight_w,
                base.cet1_capital_ratio_k,
                base.leverage_ratio_requirement_gamma,
                base.cva_capital_loading_beta,
                base.capital_hurdle_rate_h,
                base.capital_funding_offset_phi,
                base.counterparty_hazard_rate_lambda,
                base.ois_discount_rate_r,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<12>::from_f64(0.0); 12];
            for i in 0..12 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GreenBaselLeverageRatioKvaSpec {
                initial_sa_ccr_ead0: dual_in[0],
                gross_trade_notional: dual_in[1],
                sa_ccr_addon_factor_psi: dual_in[2],
                counterparty_risk_weight_w: dual_in[3],
                cet1_capital_ratio_k: dual_in[4],
                leverage_ratio_requirement_gamma: dual_in[5],
                cva_capital_loading_beta: dual_in[6],
                capital_hurdle_rate_h: dual_in[7],
                capital_funding_offset_phi: dual_in[8],
                counterparty_hazard_rate_lambda: dual_in[9],
                ois_discount_rate_r: dual_in[10],
                maturity_years_t: dual_in[11],
            };
            let r_ad = evaluate_green_basel_leverage_ratio_rwa_hurdle_kva(&spec_ad);
            let eps = 1e-6;
            for i in 0..12 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 12]| GreenBaselLeverageRatioKvaSpec {
                    initial_sa_ccr_ead0: v[0],
                    gross_trade_notional: v[1],
                    sa_ccr_addon_factor_psi: v[2],
                    counterparty_risk_weight_w: v[3],
                    cet1_capital_ratio_k: v[4],
                    leverage_ratio_requirement_gamma: v[5],
                    cva_capital_loading_beta: v[6],
                    capital_hurdle_rate_h: v[7],
                    capital_funding_offset_phi: v[8],
                    counterparty_hazard_rate_lambda: v[9],
                    ois_discount_rate_r: v[10],
                    maturity_years_t: v[11],
                };
                let rp = evaluate_green_basel_leverage_ratio_rwa_hurdle_kva(&mk(p));
                let rm = evaluate_green_basel_leverage_ratio_rwa_hurdle_kva(&mk(m));
                let fd =
                    (rp.total_dual_constraint_kva - rm.total_dual_constraint_kva) / (2.0 * eps);
                assert!(
                    r_ad.total_dual_constraint_kva.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.total_dual_constraint_kva.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use green_basel_leverage_ratio_kva_mod::*;

pub mod green_kenyon_tax_valuation_adjustment_mod {
    use ql_ad::Real;

    /// Specification for the Green (2015, Ch. 16) / Kenyon & Green (2014, Risk)
    /// Tax Valuation Adjustment (TVA) and corporate interest tax-shield XVA engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GreenKenyonTvaSpec<T: Real> {
        /// Initial Expected Positive Exposure $\text{EPE}_0 > 0$.
        pub initial_epe0: T,
        /// Initial Expected Negative Exposure $\text{ENE}_0 > 0$.
        pub initial_ene0: T,
        /// Counterparty credit hazard rate $\lambda_C > 0$.
        pub counterparty_hazard_rate_lambda_c: T,
        /// Bank own credit hazard rate $\lambda_B > 0$.
        pub bank_hazard_rate_lambda_b: T,
        /// Loss given default $\text{LGD} \in (0, 1)$.
        pub loss_given_default_lgd: T,
        /// Unsecured debt funding spread $s_F > 0$.
        pub unsecured_funding_spread_sf: T,
        ///Statutory corporate income tax rate $\tau_{\text{corp}} \in (0, 1)$.
        pub corporate_tax_rate_tau: T,
        /// Tax deductibility fraction $\phi_{\text{deduct}} \in (0, 1]$ of debt funding interest.
        pub funding_interest_deductibility_phi: T,
        /// Tax realizability fraction $\psi_{\text{dva}} \in [0, 1]$ on bank own-default DVA windfall.
        pub dva_tax_realizability_psi: T,
        /// Continuous OIS discount rate $r_{\text{ois}}$.
        pub ois_discount_rate_r: T,
        /// Trade maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Green-Kenyon Tax Valuation Adjustment (TVA) engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GreenKenyonTvaAnalytics<T: Real> {
        /// Pre-tax bilateral XVA charge ($\text{CVA} - \text{DVA} + \text{FCA} - \text{FBA}$).
        pub pre_tax_bilateral_xva: T,
        /// Tax Valuation Adjustment $\text{TVA}$ (positive means net tax shield benefit reducing XVA cost).
        pub tax_valuation_adjustment_benefit: T,
        /// Interest tax shield benefit on net funding spread $\tau_{\text{corp}} \phi_{\text{deduct}} (\text{FCA} - \text{FBA})$.
        pub funding_interest_tax_shield_benefit: T,
        /// Post-tax net XVA cost ($\text{XVA}_{\text{pre-tax}} - \text{TVA}$).
        pub post_tax_net_xva_with_tva: T,
    }

    /// Evaluates the Green (2015, Ch. 16) / Kenyon & Green (2014) Tax Valuation Adjustment (TVA)
    /// across a 4-period survival-weighted exposure profile.
    pub fn evaluate_green_kenyon_tax_valuation_adjustment_tva<T: Real>(
        spec: &GreenKenyonTvaSpec<T>,
    ) -> GreenKenyonTvaAnalytics<T> {
        const PERIODS: usize = 4;
        let t_mat = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let dt = t_mat / T::from_f64(PERIODS as f64);

        let r_tot = spec.ois_discount_rate_r
            + spec.counterparty_hazard_rate_lambda_c
            + spec.bank_hazard_rate_lambda_b;

        let mut cva = T::zero();
        let mut dva = T::zero();
        let mut fca = T::zero();
        let mut fba = T::zero();

        for m in 1..=PERIODS {
            let tm = (T::from_f64(m as f64) - T::from_f64(0.5)) * dt;
            let humped = (T::one() + T::from_f64(0.5) * tm) * (-T::from_f64(0.35) * tm).exp();
            let epe_m = spec.initial_epe0 * humped;
            let ene_m = spec.initial_ene0 * humped;

            let df_surv = (-r_tot * tm).exp();
            cva += dt
                * df_surv
                * spec.loss_given_default_lgd
                * spec.counterparty_hazard_rate_lambda_c
                * epe_m;
            dva +=
                dt * df_surv * spec.loss_given_default_lgd * spec.bank_hazard_rate_lambda_b * ene_m;
            fca += dt * df_surv * spec.unsecured_funding_spread_sf * epe_m;
            fba += dt * df_surv * spec.unsecured_funding_spread_sf * ene_m;
        }

        let pre_tax_xva = cva - dva + fca - fba;

        // Corporate tax shields:
        // 1. Default losses (CVA) are tax deductible at tau_corp, whereas own-default DVA carries tax liability psi_dva * tau_corp.
        // 2. Debt funding interest (FCA) enjoys Modigliani-Miller tax deductibility phi_deduct * tau_corp, offset by FBA tax.
        let credit_tax_shield =
            spec.corporate_tax_rate_tau * (cva - spec.dva_tax_realizability_psi * dva);
        let funding_tax_shield =
            spec.corporate_tax_rate_tau * spec.funding_interest_deductibility_phi * (fca - fba);
        let tva_benefit = credit_tax_shield + funding_tax_shield;

        GreenKenyonTvaAnalytics {
            pre_tax_bilateral_xva: pre_tax_xva,
            tax_valuation_adjustment_benefit: tva_benefit,
            funding_interest_tax_shield_benefit: funding_tax_shield,
            post_tax_net_xva_with_tva: pre_tax_xva - tva_benefit,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn green_kenyon_tva_invariants_and_ad_parity() {
            let base = GreenKenyonTvaSpec {
                initial_epe0: 25.0_f64,
                initial_ene0: 14.0,
                counterparty_hazard_rate_lambda_c: 0.022,
                bank_hazard_rate_lambda_b: 0.012,
                loss_given_default_lgd: 0.60,
                unsecured_funding_spread_sf: 0.0095,
                corporate_tax_rate_tau: 0.25,
                funding_interest_deductibility_phi: 0.85,
                dva_tax_realizability_psi: 0.40,
                ois_discount_rate_r: 0.035,
                maturity_years_t: 5.0,
            };
            let r0 = evaluate_green_kenyon_tax_valuation_adjustment_tva(&base);
            assert!(r0.pre_tax_bilateral_xva > 0.0);
            assert!(r0.tax_valuation_adjustment_benefit > 0.0);
            assert!(r0.post_tax_net_xva_with_tva < r0.pre_tax_bilateral_xva);

            let inputs = [
                base.initial_epe0,
                base.initial_ene0,
                base.counterparty_hazard_rate_lambda_c,
                base.bank_hazard_rate_lambda_b,
                base.loss_given_default_lgd,
                base.unsecured_funding_spread_sf,
                base.corporate_tax_rate_tau,
                base.funding_interest_deductibility_phi,
                base.dva_tax_realizability_psi,
                base.ois_discount_rate_r,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<11>::from_f64(0.0); 11];
            for i in 0..11 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GreenKenyonTvaSpec {
                initial_epe0: dual_in[0],
                initial_ene0: dual_in[1],
                counterparty_hazard_rate_lambda_c: dual_in[2],
                bank_hazard_rate_lambda_b: dual_in[3],
                loss_given_default_lgd: dual_in[4],
                unsecured_funding_spread_sf: dual_in[5],
                corporate_tax_rate_tau: dual_in[6],
                funding_interest_deductibility_phi: dual_in[7],
                dva_tax_realizability_psi: dual_in[8],
                ois_discount_rate_r: dual_in[9],
                maturity_years_t: dual_in[10],
            };
            let r_ad = evaluate_green_kenyon_tax_valuation_adjustment_tva(&spec_ad);
            let eps = 1e-6;
            for i in 0..11 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 11]| GreenKenyonTvaSpec {
                    initial_epe0: v[0],
                    initial_ene0: v[1],
                    counterparty_hazard_rate_lambda_c: v[2],
                    bank_hazard_rate_lambda_b: v[3],
                    loss_given_default_lgd: v[4],
                    unsecured_funding_spread_sf: v[5],
                    corporate_tax_rate_tau: v[6],
                    funding_interest_deductibility_phi: v[7],
                    dva_tax_realizability_psi: v[8],
                    ois_discount_rate_r: v[9],
                    maturity_years_t: v[10],
                };
                let rp = evaluate_green_kenyon_tax_valuation_adjustment_tva(&mk(p));
                let rm = evaluate_green_kenyon_tax_valuation_adjustment_tva(&mk(m));
                let fd =
                    (rp.post_tax_net_xva_with_tva - rm.post_tax_net_xva_with_tva) / (2.0 * eps);
                assert!(
                    r_ad.post_tax_net_xva_with_tva.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.post_tax_net_xva_with_tva.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use green_kenyon_tax_valuation_adjustment_mod::*;

// ============================================================================
// Burgard & Kjaer (2011, 2013) & Gregory (2020) Semi-Replication Funding Strategy XVA
// ============================================================================

pub mod burgard_kjaer_semi_replication_xva_mod {
    use ql_ad::Real;

    /// Specification for the Burgard & Kjaer (2011, *J. Risk* 14(2); 2013, *Risk* 26(12))
    /// PDE semi-replication bilateral counterparty and funding valuation adjustment engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BurgardKjaerSemiReplicationFundingStrategyXvaSpec<T: Real> {
        /// Initial risk-free positive mark-to-market exposure $V_0 > 0$.
        pub risk_free_derivative_value_v0: T,
        /// Initial negative mark-to-market exposure component $\text{ENE}_0 > 0$.
        pub negative_exposure_ene0: T,
        /// Counterparty default hazard rate $\lambda_C > 0$.
        pub counterparty_hazard_rate_lambda_c: T,
        /// Bank (own) default hazard rate $\lambda_B > 0$.
        pub bank_hazard_rate_lambda_b: T,
        /// Counterparty loss given default $\text{LGD}_C \in (0, 1)$.
        pub loss_given_default_lgd: T,
        /// Bank senior unsecured bond recovery rate $R_B \in (0, 1)$.
        pub own_bond_recovery_r_b: T,
        /// Semi-replication own-bond hedge weight $\alpha_{\text{BK}} \in [0, 1]$
        /// ($\alpha_{\text{BK}}=0$: no shortfall at own default; $\alpha_{\text{BK}}=1$: single-bond funding).
        pub semi_replication_hedge_weight_alpha: T,
        /// Risk-free OIS discount rate $r > 0$.
        pub ois_discount_rate_r: T,
    }

    /// Analytical outputs for the Burgard-Kjaer semi-replication XVA engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BurgardKjaerSemiReplicationFundingStrategyXvaAnalytics<T: Real> {
        /// Bilateral Credit Valuation Adjustment $\text{CVA}$.
        pub bilateral_cva_charge: T,
        /// Bilateral Debt Valuation Adjustment $\text{DVA}$.
        pub bilateral_dva_benefit: T,
        /// Semi-replication Funding Cost Adjustment $\text{FCA}(\alpha_{\text{BK}})$.
        pub semi_replication_fca_charge: T,
        /// Total Burgard-Kjaer bilateral and funding adjustment $\text{XVA}_{\text{BK}}$.
        pub total_burgard_kjaer_xva_charge: T,
    }

    /// Evaluates the Burgard & Kjaer (2011, Eqs. 19–27; 2013) semi-replication XVA charge
    /// across a 4-period joint survival-weighted exposure profile.
    pub fn evaluate_burgard_kjaer_semi_replication_funding_strategy_xva<T: Real>(
        spec: &BurgardKjaerSemiReplicationFundingStrategyXvaSpec<T>,
    ) -> BurgardKjaerSemiReplicationFundingStrategyXvaAnalytics<T> {
        const PERIODS: usize = 4;
        let dt = T::one();
        let lam_c = spec
            .counterparty_hazard_rate_lambda_c
            .max_val(T::from_f64(1e-5));
        let lam_b = spec.bank_hazard_rate_lambda_b.max_val(T::from_f64(1e-5));
        let lgd_c = spec
            .loss_given_default_lgd
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.95));
        let lgd_b = (T::one() - spec.own_bond_recovery_r_b)
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.95));
        let alpha = spec
            .semi_replication_hedge_weight_alpha
            .max_val(T::zero())
            .min_val(T::one());

        // Under Burgard-Kjaer (2011), hedging fraction alpha of own-default windfall reduces
        // the effective funding rate on V^+ from s_F = (1 - R_B) * lambda_B to (1 - alpha) * s_F
        // while introducing replacement close-out contagion alpha * s_F * lambda_C * dt.
        let s_f = lgd_b * lam_b;
        let r_tot = spec.ois_discount_rate_r + lam_c + lam_b;

        let mut cva = T::zero();
        let mut dva = T::zero();
        let mut fca = T::zero();

        for m in 1..=PERIODS {
            let tm = (T::from_f64(m as f64) - T::from_f64(0.5)) * dt;
            let epe_m = spec.risk_free_derivative_value_v0
                * (T::one() + T::from_f64(0.35) * tm)
                * (-T::from_f64(0.25) * tm).exp();
            let ene_m = spec.negative_exposure_ene0
                * (T::one() + T::from_f64(0.25) * tm)
                * (-T::from_f64(0.30) * tm).exp();

            let df_surv = (-r_tot * tm).exp();
            cva += dt * df_surv * lgd_c * lam_c * epe_m;
            dva += dt * df_surv * lgd_b * lam_b * ene_m;
            let effective_fca_spread = (T::one() - alpha) * s_f + alpha * s_f * lgd_c * lam_c * tm;
            fca += dt * df_surv * effective_fca_spread * epe_m;
        }

        let total_xva = cva - dva + fca;

        BurgardKjaerSemiReplicationFundingStrategyXvaAnalytics {
            bilateral_cva_charge: cva,
            bilateral_dva_benefit: dva,
            semi_replication_fca_charge: fca,
            total_burgard_kjaer_xva_charge: total_xva,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn burgard_kjaer_semi_replication_xva_invariants_and_ad_parity() {
            let base = BurgardKjaerSemiReplicationFundingStrategyXvaSpec {
                risk_free_derivative_value_v0: 100.0_f64,
                negative_exposure_ene0: 45.0,
                counterparty_hazard_rate_lambda_c: 0.022,
                bank_hazard_rate_lambda_b: 0.015,
                loss_given_default_lgd: 0.60,
                own_bond_recovery_r_b: 0.40,
                semi_replication_hedge_weight_alpha: 0.65,
                ois_discount_rate_r: 0.035,
            };
            let r0 = evaluate_burgard_kjaer_semi_replication_funding_strategy_xva(&base);
            assert!(r0.bilateral_cva_charge > 0.0);
            assert!(r0.bilateral_dva_benefit > 0.0);
            assert!(r0.semi_replication_fca_charge > 0.0);
            assert!(r0.total_burgard_kjaer_xva_charge > 0.0);

            let inputs = [
                base.risk_free_derivative_value_v0,
                base.negative_exposure_ene0,
                base.counterparty_hazard_rate_lambda_c,
                base.bank_hazard_rate_lambda_b,
                base.loss_given_default_lgd,
                base.own_bond_recovery_r_b,
                base.semi_replication_hedge_weight_alpha,
                base.ois_discount_rate_r,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BurgardKjaerSemiReplicationFundingStrategyXvaSpec {
                risk_free_derivative_value_v0: dual_in[0],
                negative_exposure_ene0: dual_in[1],
                counterparty_hazard_rate_lambda_c: dual_in[2],
                bank_hazard_rate_lambda_b: dual_in[3],
                loss_given_default_lgd: dual_in[4],
                own_bond_recovery_r_b: dual_in[5],
                semi_replication_hedge_weight_alpha: dual_in[6],
                ois_discount_rate_r: dual_in[7],
            };
            let r_ad = evaluate_burgard_kjaer_semi_replication_funding_strategy_xva(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| BurgardKjaerSemiReplicationFundingStrategyXvaSpec {
                    risk_free_derivative_value_v0: v[0],
                    negative_exposure_ene0: v[1],
                    counterparty_hazard_rate_lambda_c: v[2],
                    bank_hazard_rate_lambda_b: v[3],
                    loss_given_default_lgd: v[4],
                    own_bond_recovery_r_b: v[5],
                    semi_replication_hedge_weight_alpha: v[6],
                    ois_discount_rate_r: v[7],
                };
                let rp = evaluate_burgard_kjaer_semi_replication_funding_strategy_xva(&mk(p));
                let rm = evaluate_burgard_kjaer_semi_replication_funding_strategy_xva(&mk(m));
                let fd = (rp.total_burgard_kjaer_xva_charge - rm.total_burgard_kjaer_xva_charge)
                    / (2.0 * eps);
                assert!(
                    r_ad.total_burgard_kjaer_xva_charge.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.total_burgard_kjaer_xva_charge.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use burgard_kjaer_semi_replication_xva_mod::*;

// ============================================================================
// Gregory (2020, Ch. 15) & Andersen-Pykhtin-Sokol (2017) Dynamic ISDA SIMM MVA
// ============================================================================

pub mod gregory_isda_simm_mva_mod {
    use ql_ad::Real;

    /// Specification for the Gregory (2020, *The xVA Challenge*, 4th ed., Ch. 15) and
    /// Andersen, Pykhtin & Sokol (2017) dynamic ISDA SIMM Initial Margin and MVA engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GregoryIsdaSimmMvaSpec<T: Real> {
        /// Initial net SIMM delta-equivalent risk weight / DV01 exposure $W_0 > 0$.
        pub trade_notional_dv01_w0: T,
        /// ISDA SIMM stress-calibrated annualized volatility $\sigma_{\text{SIMM}} > 0$.
        pub simm_calibrated_volatility_sigma_im: T,
        /// Margin Period of Risk in business days $\text{MPOR}_{\text{days}} > 0$ (e.g., 10 days).
        pub margin_period_of_risk_days_mpor: T,
        /// Initial Margin funding spread over CSA remuneration rate $s_{\text{IM}} > 0$.
        pub im_funding_spread_s_im: T,
        /// Counterparty default hazard rate $\lambda_C > 0$.
        pub counterparty_hazard_rate_lambda_c: T,
        /// Bank (own) default hazard rate $\lambda_B > 0$.
        pub bank_hazard_rate_lambda_b: T,
        /// OIS risk-free discount rate $r > 0$.
        pub ois_discount_rate_r: T,
        /// Portfolio maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Gregory ISDA SIMM Margin Valuation Adjustment (MVA) engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GregoryIsdaSimmMvaAnalytics<T: Real> {
        /// Inception 99% ISDA SIMM Initial Margin requirement $\text{IM}_0$.
        pub inception_simm_initial_margin: T,
        /// Time-averaged expected forward Initial Margin $\overline{\text{IM}}$.
        pub expected_average_forward_im: T,
        /// Survival-discounted expected IM integral $\int_0^T D_{\text{surv}}(0, t)\mathbb{E}[\text{IM}(t)]\,dt$.
        pub discounted_im_exposure_integral: T,
        /// Lifetime Margin Valuation Adjustment charge $\text{MVA}$.
        pub lifetime_simm_mva_charge: T,
    }

    /// Evaluates the Gregory (2020, Ch. 15) / Green (2015, Ch. 13) lifetime Margin Valuation
    /// Adjustment (MVA) across a 4-period amortizing ISDA SIMM profile.
    pub fn evaluate_gregory_isda_simm_dynamic_initial_margin_mva<T: Real>(
        spec: &GregoryIsdaSimmMvaSpec<T>,
    ) -> GregoryIsdaSimmMvaAnalytics<T> {
        const PERIODS: usize = 4;
        let w0 = spec.trade_notional_dv01_w0.max_val(T::from_f64(1e-4));
        let sig_im = spec
            .simm_calibrated_volatility_sigma_im
            .max_val(T::from_f64(1e-5));
        let mpor_years =
            spec.margin_period_of_risk_days_mpor.max_val(T::one()) / T::from_f64(250.0);
        let s_im = spec.im_funding_spread_s_im.max_val(T::from_f64(1e-5));
        let t_mat = spec.maturity_years_t.max_val(T::from_f64(0.25));
        let dt = t_mat / T::from_f64(PERIODS as f64);

        // 99% one-tailed normal quantile z_0.99 = 2.3263478740408408
        let z99 = T::from_f64(2.326_347_874_040_840_8);
        let im0 = z99 * w0 * sig_im * mpor_years.sqrt();

        let r_tot = spec.ois_discount_rate_r
            + spec.counterparty_hazard_rate_lambda_c
            + spec.bank_hazard_rate_lambda_b;

        let mut sum_im = T::zero();
        let mut disc_im_int = T::zero();

        for m in 1..=PERIODS {
            let tm = (T::from_f64(m as f64) - T::from_f64(0.5)) * dt;
            let rem_fraction = T::one() - tm / t_mat;
            let disp_factor = T::one() + T::from_f64(2.5) * sig_im * tm;
            let e_im_m = im0 * rem_fraction * disp_factor;

            sum_im += e_im_m;
            disc_im_int += dt * (-r_tot * tm).exp() * e_im_m;
        }

        let avg_im = sum_im / T::from_f64(PERIODS as f64);
        let mva = s_im * disc_im_int;

        GregoryIsdaSimmMvaAnalytics {
            inception_simm_initial_margin: im0,
            expected_average_forward_im: avg_im,
            discounted_im_exposure_integral: disc_im_int,
            lifetime_simm_mva_charge: mva,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gregory_isda_simm_mva_invariants_and_ad_parity() {
            let base = GregoryIsdaSimmMvaSpec {
                trade_notional_dv01_w0: 120.0_f64,
                simm_calibrated_volatility_sigma_im: 0.011,
                margin_period_of_risk_days_mpor: 10.0,
                im_funding_spread_s_im: 0.0085,
                counterparty_hazard_rate_lambda_c: 0.018,
                bank_hazard_rate_lambda_b: 0.012,
                ois_discount_rate_r: 0.035,
                maturity_years_t: 5.0,
            };
            let r0 = evaluate_gregory_isda_simm_dynamic_initial_margin_mva(&base);
            assert!(r0.inception_simm_initial_margin > 0.0);
            assert!(r0.expected_average_forward_im > 0.0);
            assert!(r0.discounted_im_exposure_integral > 0.0);
            assert!(r0.lifetime_simm_mva_charge > 0.0);

            let inputs = [
                base.trade_notional_dv01_w0,
                base.simm_calibrated_volatility_sigma_im,
                base.margin_period_of_risk_days_mpor,
                base.im_funding_spread_s_im,
                base.counterparty_hazard_rate_lambda_c,
                base.bank_hazard_rate_lambda_b,
                base.ois_discount_rate_r,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GregoryIsdaSimmMvaSpec {
                trade_notional_dv01_w0: dual_in[0],
                simm_calibrated_volatility_sigma_im: dual_in[1],
                margin_period_of_risk_days_mpor: dual_in[2],
                im_funding_spread_s_im: dual_in[3],
                counterparty_hazard_rate_lambda_c: dual_in[4],
                bank_hazard_rate_lambda_b: dual_in[5],
                ois_discount_rate_r: dual_in[6],
                maturity_years_t: dual_in[7],
            };
            let r_ad = evaluate_gregory_isda_simm_dynamic_initial_margin_mva(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| GregoryIsdaSimmMvaSpec {
                    trade_notional_dv01_w0: v[0],
                    simm_calibrated_volatility_sigma_im: v[1],
                    margin_period_of_risk_days_mpor: v[2],
                    im_funding_spread_s_im: v[3],
                    counterparty_hazard_rate_lambda_c: v[4],
                    bank_hazard_rate_lambda_b: v[5],
                    ois_discount_rate_r: v[6],
                    maturity_years_t: v[7],
                };
                let rp = evaluate_gregory_isda_simm_dynamic_initial_margin_mva(&mk(p));
                let rm = evaluate_gregory_isda_simm_dynamic_initial_margin_mva(&mk(m));
                let fd = (rp.lifetime_simm_mva_charge - rm.lifetime_simm_mva_charge) / (2.0 * eps);
                assert!(
                    r_ad.lifetime_simm_mva_charge.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.lifetime_simm_mva_charge.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use gregory_isda_simm_mva_mod::*;
mod pykhtin_rosen_marginal_cva_mod {
    use ql_ad::Real;

    /// Input specification for the Pykhtin & Rosen (2010, *Risk*) / Gregory (2020, Ch. 13)
    /// Euler homogeneous marginal CVA allocation across netting-set subportfolios.
    #[derive(Debug, Clone, Copy)]
    pub struct PykhtinRosenMarginalCvaSpec<T: Real> {
        /// Expected MTM $\mu_1$ of subportfolio 1 (e.g., $12.0$).
        pub subportfolio_1_expected_mtm_mu1: T,
        /// Expected MTM $\mu_2$ of subportfolio 2 (e.g., $-4.5$).
        pub subportfolio_2_expected_mtm_mu2: T,
        /// MTM volatility $\sigma_1$ of subportfolio 1 (e.g., $28.0$).
        pub subportfolio_1_volatility_sigma1: T,
        /// MTM volatility $\sigma_2$ of subportfolio 2 (e.g., $22.0$).
        pub subportfolio_2_volatility_sigma2: T,
        /// Cross-subportfolio MTM correlation $\rho_{12} \in (-1, 1)$ (e.g., $-0.35$).
        pub mtm_correlation_rho12: T,
        /// Counterparty hazard rate $\lambda_C$ (e.g., $0.024$).
        pub counterparty_hazard_lambda_c: T,
        /// Counterparty recovery rate $R_C \in (0, 1)$ (e.g., $0.40$).
        pub counterparty_recovery_rate_rc: T,
        /// Risk-free discount rate $r$ (e.g., $0.032$).
        pub risk_free_rate_r: T,
    }

    /// Analytics returned by [`evaluate_pykhtin_rosen_euler_marginal_cva_netting_allocation`].
    #[derive(Debug, Clone, Copy)]
    pub struct PykhtinRosenMarginalCvaAnalytics<T: Real> {
        /// Netting-set Expected Positive Exposure $\text{EPE}_P = \mathbb{E}[\max(V_1 + V_2, 0)]$.
        pub netting_set_epe: T,
        /// Pykhtin–Rosen Euler marginal CVA allocated to subportfolio 1.
        pub marginal_cva_subportfolio_1: T,
        /// Pykhtin–Rosen Euler marginal CVA allocated to subportfolio 2.
        pub marginal_cva_subportfolio_2: T,
        /// Total netting-set CVA $\text{CVA}_P = \text{CVA}_1^{\text{marg}} + \text{CVA}_2^{\text{marg}}$.
        pub total_netting_set_cva: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        T::from_f64(0.5) * (T::one() + (x * T::from_f64(0.707_106_781_186_547_5)).erf())
    }

    #[inline]
    fn std_norm_pdf<T: Real>(x: T) -> T {
        T::from_f64(0.398_942_280_401_432_7) * (-T::from_f64(0.5) * x * x).exp()
    }

    /// Evaluates the exact Pykhtin & Rosen (2010) Euler marginal EPE and CVA decomposition
    /// $\text{EPE}_i^{\text{marg}} = \mathbb{E}[V_i \mathbf{1}_{\{V_P > 0\}}] = \mu_i \Phi(\mu_P/\sigma_P) + \frac{\operatorname{Cov}(V_i, V_P)}{\sigma_P}\phi(\mu_P/\sigma_P)$,
    /// satisfying exact Euler additivity $\text{CVA}_1^{\text{marg}} + \text{CVA}_2^{\text{marg}} \equiv \text{CVA}_P$.
    pub fn evaluate_pykhtin_rosen_euler_marginal_cva_netting_allocation<T: Real>(
        spec: &PykhtinRosenMarginalCvaSpec<T>,
    ) -> PykhtinRosenMarginalCvaAnalytics<T> {
        let mu1 = spec.subportfolio_1_expected_mtm_mu1;
        let mu2 = spec.subportfolio_2_expected_mtm_mu2;
        let s1 = spec.subportfolio_1_volatility_sigma1;
        let s2 = spec.subportfolio_2_volatility_sigma2;
        let rho = spec.mtm_correlation_rho12;

        let mu_p = mu1 + mu2;
        let cov_1_p = s1 * s1 + rho * s1 * s2;
        let cov_2_p = s2 * s2 + rho * s1 * s2;
        let var_p = cov_1_p + cov_2_p;
        let sig_p = var_p.sqrt();

        let z = mu_p / sig_p;
        let phi_cdf = std_norm_cdf(z);
        let phi_pdf = std_norm_pdf(z);

        let epe_1_marg = mu1 * phi_cdf + (cov_1_p / sig_p) * phi_pdf;
        let epe_2_marg = mu2 * phi_cdf + (cov_2_p / sig_p) * phi_pdf;
        let epe_p = epe_1_marg + epe_2_marg;

        // 5-year risky default discount weight (1 - R_C) * (lambda_C / (r + lambda_C)) * (1 - exp(-(r + lambda_C)*5))
        let lgd_c = T::one() - spec.counterparty_recovery_rate_rc;
        let lam_c = spec.counterparty_hazard_lambda_c;
        let r = spec.risk_free_rate_r;
        let horizon = T::from_f64(5.0);
        let cva_mult = lgd_c * (lam_c / (r + lam_c)) * (T::one() - (-(r + lam_c) * horizon).exp());

        let cva_1 = cva_mult * epe_1_marg;
        let cva_2 = cva_mult * epe_2_marg;
        let cva_tot = cva_1 + cva_2;

        PykhtinRosenMarginalCvaAnalytics {
            netting_set_epe: epe_p,
            marginal_cva_subportfolio_1: cva_1,
            marginal_cva_subportfolio_2: cva_2,
            total_netting_set_cva: cva_tot,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn pykhtin_rosen_marginal_cva_invariants_and_ad_parity() {
            let base = PykhtinRosenMarginalCvaSpec {
                subportfolio_1_expected_mtm_mu1: 12.0_f64,
                subportfolio_2_expected_mtm_mu2: -4.5_f64,
                subportfolio_1_volatility_sigma1: 28.0_f64,
                subportfolio_2_volatility_sigma2: 22.0_f64,
                mtm_correlation_rho12: -0.35_f64,
                counterparty_hazard_lambda_c: 0.024_f64,
                counterparty_recovery_rate_rc: 0.40_f64,
                risk_free_rate_r: 0.032_f64,
            };
            let res = evaluate_pykhtin_rosen_euler_marginal_cva_netting_allocation(&base);
            assert!(res.netting_set_epe > 0.0);
            assert!(res.total_netting_set_cva > 0.0);
            assert!(
                (res.marginal_cva_subportfolio_1 + res.marginal_cva_subportfolio_2
                    - res.total_netting_set_cva)
                    .abs()
                    < 1e-12
            );

            let inputs = [
                base.subportfolio_1_expected_mtm_mu1,
                base.subportfolio_2_expected_mtm_mu2,
                base.subportfolio_1_volatility_sigma1,
                base.subportfolio_2_volatility_sigma2,
                base.mtm_correlation_rho12,
                base.counterparty_hazard_lambda_c,
                base.counterparty_recovery_rate_rc,
                base.risk_free_rate_r,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = PykhtinRosenMarginalCvaSpec {
                subportfolio_1_expected_mtm_mu1: dual_in[0],
                subportfolio_2_expected_mtm_mu2: dual_in[1],
                subportfolio_1_volatility_sigma1: dual_in[2],
                subportfolio_2_volatility_sigma2: dual_in[3],
                mtm_correlation_rho12: dual_in[4],
                counterparty_hazard_lambda_c: dual_in[5],
                counterparty_recovery_rate_rc: dual_in[6],
                risk_free_rate_r: dual_in[7],
            };
            let r_ad = evaluate_pykhtin_rosen_euler_marginal_cva_netting_allocation(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| PykhtinRosenMarginalCvaSpec {
                    subportfolio_1_expected_mtm_mu1: v[0],
                    subportfolio_2_expected_mtm_mu2: v[1],
                    subportfolio_1_volatility_sigma1: v[2],
                    subportfolio_2_volatility_sigma2: v[3],
                    mtm_correlation_rho12: v[4],
                    counterparty_hazard_lambda_c: v[5],
                    counterparty_recovery_rate_rc: v[6],
                    risk_free_rate_r: v[7],
                };
                let rp = evaluate_pykhtin_rosen_euler_marginal_cva_netting_allocation(&mk(p));
                let rm = evaluate_pykhtin_rosen_euler_marginal_cva_netting_allocation(&mk(m));
                let fd = (rp.total_netting_set_cva - rm.total_netting_set_cva) / (2.0 * eps);
                assert!(
                    r_ad.total_netting_set_cva.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.total_netting_set_cva.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use pykhtin_rosen_marginal_cva_mod::*;
