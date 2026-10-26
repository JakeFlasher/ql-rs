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
        let eps = T::from_f64(0.05);
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
