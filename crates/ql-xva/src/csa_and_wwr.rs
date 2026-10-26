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
