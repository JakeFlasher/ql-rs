mod evaluate_damodaran_rsm_country_risk_and_debt_modification_module {
    use ql_ad::Real;

    #[derive(Clone, Debug)]
    pub struct DamodaranRsmCountryRiskAndDebtModSpec<T: Real = f64> {
        pub mature_market_erp: T,
        pub sovereign_cds_spread: T,
        pub sovereign_default_drag_bps: T,
        pub emerging_equity_vol: T,
        pub emerging_govt_bond_vol: T,
        pub firm_domestic_revenue_share: T,
        pub ebit_operating_income: T,
        pub annual_interest_expense: T,
        pub original_debt_carrying_value: T,
        pub original_effective_interest_rate: T,
        pub modified_cashflows: Vec<T>,
        pub modified_payment_times: Vec<T>,
        pub modification_fees_paid: T,
    }
    #[derive(Clone, Debug)]
    pub struct DamodaranRsmCountryRiskAndDebtModReport<T: Real = f64> {
        pub net_sovereign_default_spread: T,
        pub relative_equity_to_bond_vol_multiplier: T,
        pub country_risk_premium_crp: T,
        pub firm_exposure_lambda: T,
        pub total_cost_of_equity_ke: T,
        pub interest_coverage_ratio: T,
        pub synthetic_corporate_default_spread: T,
        pub discounted_modified_debt_pv: T,
        pub rsm_pv_percentage_change: T,
        pub is_substantial_extinguishment_10pct: bool,
        pub extinguishment_gain_or_loss: T,
    }
    pub fn evaluate_damodaran_rsm_country_risk_and_debt_modification<T: Real>(
        s: &DamodaranRsmCountryRiskAndDebtModSpec<T>,
    ) -> DamodaranRsmCountryRiskAndDebtModReport<T> {
        let net = s.sovereign_cds_spread - s.sovereign_default_drag_bps * T::from_f64(1e-4);
        let rel =
            s.emerging_equity_vol / s.emerging_govt_bond_vol.abs().max_val(T::from_f64(1e-12));
        let crp = net * rel;
        let lambda = s
            .firm_domestic_revenue_share
            .max_val(T::zero())
            .min_val(T::one());
        let ke = s.mature_market_erp + lambda * crp;
        let coverage =
            s.ebit_operating_income / s.annual_interest_expense.abs().max_val(T::from_f64(1e-12));
        let synth = net + T::from_f64(0.01) / (coverage.abs() + T::one());
        let mut pv = T::zero();
        let mut pv_fair = T::zero();
        for (cf, t) in s
            .modified_cashflows
            .iter()
            .zip(s.modified_payment_times.iter())
        {
            pv += *cf * (-s.original_effective_interest_rate * *t).exp();
            pv_fair += *cf * (-(s.original_effective_interest_rate + synth) * *t).exp();
        }
        let denom = s
            .original_debt_carrying_value
            .abs()
            .max_val(T::from_f64(1e-12));
        let pct = (pv + s.modification_fees_paid - s.original_debt_carrying_value).abs() / denom;
        let substantial = pct >= T::from_f64(0.10);
        let gain = s.original_debt_carrying_value - pv_fair - s.modification_fees_paid;
        DamodaranRsmCountryRiskAndDebtModReport {
            net_sovereign_default_spread: net,
            relative_equity_to_bond_vol_multiplier: rel,
            country_risk_premium_crp: crp,
            firm_exposure_lambda: lambda,
            total_cost_of_equity_ke: ke,
            interest_coverage_ratio: coverage,
            synthetic_corporate_default_spread: synth,
            discounted_modified_debt_pv: pv,
            rsm_pv_percentage_change: pct,
            is_substantial_extinguishment_10pct: substantial,
            extinguishment_gain_or_loss: gain,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct FsbFinmaGsibResolutionBailInWaterfallSpec<T: Real = f64> {
        pub total_rwa: T,
        pub pre_resolution_cet1: T,
        pub at1_principal: T,
        pub tier2_principal: T,
        pub senior_tlac_principal: T,
        pub valuation_loss_shock: T,
        pub target_post_resolution_cet1_ratio: T,
        pub liquidation_insolvency_haircut: T,
        pub existing_share_count: T,
        pub bail_in_conversion_discount: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct FsbFinmaGsibResolutionBailInWaterfallReport<T: Real = f64> {
        pub total_loss_plus_recap_need: T,
        pub cet1_absorbed_loss: T,
        pub at1_written_down_or_converted: T,
        pub tier2_bailed_in: T,
        pub senior_tlac_bailed_in: T,
        pub senior_tlac_resolution_recovery_rate: T,
        pub senior_tlac_insolvency_counterfactual_recovery: T,
        pub ncwo_compensation_shortfall: T,
        pub post_bail_in_creditor_ownership_share: T,
        pub waterfall_conservation_error: T,
    }
    pub fn evaluate_fsb_finma_gsib_resolution_bail_in_waterfall<T: Real>(
        s: &FsbFinmaGsibResolutionBailInWaterfallSpec<T>,
    ) -> FsbFinmaGsibResolutionBailInWaterfallReport<T> {
        let pre_cet1 = s.pre_resolution_cet1.max_val(T::zero());
        let loss = (s.total_rwa * s.valuation_loss_shock).max_val(T::zero());
        let cet = loss.min_val(pre_cet1);
        let post_loss_cet1 = (pre_cet1 - cet).max_val(T::zero());
        let target_cet1 = (s.target_post_resolution_cet1_ratio * s.total_rwa).max_val(T::zero());
        let recap = (target_cet1 - post_loss_cet1).max_val(T::zero());
        let need = loss + recap;
        let rem1 = (need - cet).max_val(T::zero());
        let at1 = rem1.min_val(s.at1_principal.max_val(T::zero()));
        let rem2 = (rem1 - at1).max_val(T::zero());
        let tier2 = rem2.min_val(s.tier2_principal.max_val(T::zero()));
        let rem3 = (rem2 - tier2).max_val(T::zero());
        let senior = rem3.min_val(s.senior_tlac_principal.max_val(T::zero()));
        let senior_den = s.senior_tlac_principal.max_val(T::from_f64(1e-12));
        let res_rec = ((s.senior_tlac_principal - senior) / senior_den)
            .max_val(T::zero())
            .min_val(T::one());
        let ins_rec = (T::one() - s.liquidation_insolvency_haircut)
            .max_val(T::zero())
            .min_val(T::one());
        let ncwo = (ins_rec - res_rec).max_val(T::zero()) * s.senior_tlac_principal;
        let converted = at1 + tier2 + senior;
        let new_shares =
            converted / (T::one() - s.bail_in_conversion_discount).max_val(T::from_f64(1e-6));
        let owner = new_shares / (new_shares + s.existing_share_count.max_val(T::from_f64(1e-6)));
        let total_capacity = pre_cet1
            + s.at1_principal.max_val(T::zero())
            + s.tier2_principal.max_val(T::zero())
            + s.senior_tlac_principal.max_val(T::zero());
        let conservation = (need.min_val(total_capacity) - (cet + at1 + tier2 + senior)).abs();
        FsbFinmaGsibResolutionBailInWaterfallReport {
            total_loss_plus_recap_need: need,
            cet1_absorbed_loss: cet,
            at1_written_down_or_converted: at1,
            tier2_bailed_in: tier2,
            senior_tlac_bailed_in: senior,
            senior_tlac_resolution_recovery_rate: res_rec,
            senior_tlac_insolvency_counterfactual_recovery: ins_rec,
            ncwo_compensation_shortfall: ncwo,
            post_bail_in_creditor_ownership_share: owner,
            waterfall_conservation_error: conservation,
        }
    }

    #[cfg(test)]
    mod domain_suite_42_tests {
        use super::*;
        use ql_ad::Dual;
        fn make<T: Real>(x: T) -> DamodaranRsmCountryRiskAndDebtModSpec<T> {
            DamodaranRsmCountryRiskAndDebtModSpec {
                mature_market_erp: T::from_f64(0.05),
                sovereign_cds_spread: T::from_f64(0.01) + x,
                sovereign_default_drag_bps: T::from_f64(10.0),
                emerging_equity_vol: T::from_f64(0.3),
                emerging_govt_bond_vol: T::from_f64(0.1),
                firm_domestic_revenue_share: T::from_f64(0.5),
                ebit_operating_income: T::from_f64(10.0),
                annual_interest_expense: T::from_f64(2.0),
                original_debt_carrying_value: T::from_f64(100.0),
                original_effective_interest_rate: T::from_f64(0.03),
                modified_cashflows: vec![T::from_f64(4.0), T::from_f64(104.0)],
                modified_payment_times: vec![T::one(), T::from_f64(2.0)],
                modification_fees_paid: T::one(),
            }
        }
        #[test]
        fn suite42_credit_and_bailin_dual_fd() {
            let r = evaluate_damodaran_rsm_country_risk_and_debt_modification(&make(0.0));
            assert!(
                (r.country_risk_premium_crp
                    - r.net_sovereign_default_spread * r.relative_equity_to_bond_vol_multiplier)
                    .abs()
                    < 1e-12
            );
            let h = 1e-6;
            let fd = (evaluate_damodaran_rsm_country_risk_and_debt_modification(&make(h))
                .country_risk_premium_crp
                - evaluate_damodaran_rsm_country_risk_and_debt_modification(&make(-h))
                    .country_risk_premium_crp)
                / (2.0 * h);
            let ad = evaluate_damodaran_rsm_country_risk_and_debt_modification(&make(
                Dual::<1>::variable(0.0, 0),
            ))
            .country_risk_premium_crp;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
            let mut substantial = make(0.0);
            substantial.modified_cashflows = vec![4.0, 120.0];
            let sr = evaluate_damodaran_rsm_country_risk_and_debt_modification(&substantial);
            assert!(sr.is_substantial_extinguishment_10pct);
            let b = evaluate_fsb_finma_gsib_resolution_bail_in_waterfall(
                &FsbFinmaGsibResolutionBailInWaterfallSpec {
                    total_rwa: 1000.0,
                    pre_resolution_cet1: 120.0,
                    at1_principal: 50.0,
                    tier2_principal: 60.0,
                    senior_tlac_principal: 100.0,
                    valuation_loss_shock: 0.05,
                    target_post_resolution_cet1_ratio: 0.1,
                    liquidation_insolvency_haircut: 0.4,
                    existing_share_count: 100.0,
                    bail_in_conversion_discount: 0.2,
                },
            );
            assert!(b.waterfall_conservation_error.abs() < 1e-12);
            assert_eq!(b.tier2_bailed_in, 0.0);
            assert_eq!(b.at1_written_down_or_converted, 30.0);
            let h = 1e-6;
            let base = FsbFinmaGsibResolutionBailInWaterfallSpec {
                total_rwa: 1000.0,
                pre_resolution_cet1: 120.0,
                at1_principal: 50.0,
                tier2_principal: 60.0,
                senior_tlac_principal: 100.0,
                valuation_loss_shock: 0.08,
                target_post_resolution_cet1_ratio: 0.1,
                liquidation_insolvency_haircut: 0.4,
                existing_share_count: 100.0,
                bail_in_conversion_discount: 0.2,
            };
            let mut up = base;
            let mut dn = base;
            up.valuation_loss_shock += h;
            dn.valuation_loss_shock -= h;
            let fd_owner = (evaluate_fsb_finma_gsib_resolution_bail_in_waterfall(&up)
                .post_bail_in_creditor_ownership_share
                - evaluate_fsb_finma_gsib_resolution_bail_in_waterfall(&dn)
                    .post_bail_in_creditor_ownership_share)
                / (2.0 * h);
            let ad_spec: FsbFinmaGsibResolutionBailInWaterfallSpec<Dual<1>> =
                FsbFinmaGsibResolutionBailInWaterfallSpec {
                    total_rwa: Dual::from_f64(1000.0),
                    pre_resolution_cet1: Dual::from_f64(120.0),
                    at1_principal: Dual::from_f64(50.0),
                    tier2_principal: Dual::from_f64(60.0),
                    senior_tlac_principal: Dual::from_f64(100.0),
                    valuation_loss_shock: Dual::<1>::variable(0.08, 0),
                    target_post_resolution_cet1_ratio: Dual::from_f64(0.1),
                    liquidation_insolvency_haircut: Dual::from_f64(0.4),
                    existing_share_count: Dual::from_f64(100.0),
                    bail_in_conversion_discount: Dual::from_f64(0.2),
                };
            let ad_owner = evaluate_fsb_finma_gsib_resolution_bail_in_waterfall(&ad_spec)
                .post_bail_in_creditor_ownership_share;
            assert!(ad_owner.d[0].abs() > 1e-6 && (ad_owner.d[0] - fd_owner).abs() < 1e-5);
        }
    }
}
pub use evaluate_damodaran_rsm_country_risk_and_debt_modification_module::*;

mod evaluate_rsm_asc470_debt_modification_and_restructuring_module {
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct RsmAsc470DebtModificationSpec<T: Real = f64> {
        pub carrying_amount: T,
        pub original_par: T,
        pub original_coupon_rate: T,
        pub original_effective_rate: T,
        pub modified_par: T,
        pub modified_coupon_rate: T,
        pub remaining_years: T,
        pub lender_fee_paid: T,
        pub conversion_option_fv_change: T,
        pub market_exit_yield: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct RsmAsc470DebtModificationReport<T: Real = f64> {
        pub pv_original_cashflows_at_orig_eir: T,
        pub pv_modified_cashflows_at_orig_eir: T,
        pub asc470_50_pv_difference_pct: T,
        pub smooth_extinguishment_weight: T,
        pub extinguishment_gain_or_loss: T,
        pub prospective_modified_eir: T,
        pub tdr_creditor_concession_loss: T,
    }

    fn annuity<T: Real>(r: T, n: T) -> T {
        let rr = r.max_val(T::from_f64(1e-10));
        (T::one() - (-((T::one() + rr).ln() * n)).exp()) / rr
    }

    fn pv_mod<T: Real>(s: &RsmAsc470DebtModificationSpec<T>, r: T) -> T {
        let rr = r.max_val(T::from_f64(1e-10));
        s.modified_coupon_rate * s.modified_par * annuity(rr, s.remaining_years)
            + s.modified_par * (-((T::one() + rr).ln() * s.remaining_years)).exp()
            + s.lender_fee_paid
            + s.conversion_option_fv_change
    }

    pub fn evaluate_rsm_asc470_debt_modification_and_restructuring<T: Real>(
        s: &RsmAsc470DebtModificationSpec<T>,
    ) -> RsmAsc470DebtModificationReport<T> {
        let pv_orig = s.original_coupon_rate
            * s.original_par
            * annuity(s.original_effective_rate, s.remaining_years)
            + s.original_par
                * (-((T::one() + s.original_effective_rate).ln() * s.remaining_years)).exp();
        let pv_modified = pv_mod(s, s.original_effective_rate);
        let pct = (pv_modified - pv_orig).abs() / pv_orig.abs().max_val(T::from_f64(1e-12));
        let weight =
            T::one() / (T::one() + (-(T::from_f64(50.0) * (pct - T::from_f64(0.10)))).exp());
        let mut r = s.market_exit_yield.max_val(T::from_f64(1e-6));
        for _ in 0..8 {
            let pv = pv_mod(s, r);
            let rr = r.max_val(T::from_f64(1e-8));
            let n = s.remaining_years;
            let disc = (-((T::one() + rr).ln() * n)).exp();
            let d_ann = (n * disc * rr / (T::one() + rr) - (T::one() - disc)) / (rr * rr);
            let deriv = s.modified_coupon_rate * s.modified_par * d_ann
                - s.modified_par * n * disc / (T::one() + rr);
            let safe_deriv = if deriv.abs() < T::from_f64(1e-10) {
                T::from_f64(-1e-10)
            } else {
                deriv
            };
            r -= (pv - s.carrying_amount) / safe_deriv;
        }
        let concession = (s.carrying_amount - pv_modified).max_val(T::zero());
        RsmAsc470DebtModificationReport {
            pv_original_cashflows_at_orig_eir: pv_orig,
            pv_modified_cashflows_at_orig_eir: pv_modified,
            asc470_50_pv_difference_pct: pct,
            smooth_extinguishment_weight: weight,
            extinguishment_gain_or_loss: s.carrying_amount - pv_modified,
            prospective_modified_eir: r,
            tdr_creditor_concession_loss: concession,
        }
    }

    #[cfg(test)]
    mod domain_suite_43_tests {
        use super::*;
        use ql_ad::Dual;

        fn spec<T: Real>(x: T) -> RsmAsc470DebtModificationSpec<T> {
            RsmAsc470DebtModificationSpec {
                carrying_amount: T::from_f64(100.0) + x,
                original_par: T::from_f64(100.0),
                original_coupon_rate: T::from_f64(0.03),
                original_effective_rate: T::from_f64(0.04),
                modified_par: T::from_f64(100.0),
                modified_coupon_rate: T::from_f64(0.025),
                remaining_years: T::from_f64(5.0),
                lender_fee_paid: T::one(),
                conversion_option_fv_change: T::from_f64(0.5),
                market_exit_yield: T::from_f64(0.05),
            }
        }

        #[test]
        fn suite43_invariants_and_dual_fd() {
            let s0 = spec(0.0);
            let r = evaluate_rsm_asc470_debt_modification_and_restructuring(&s0);
            assert!(r.asc470_50_pv_difference_pct >= 0.0);
            assert!(r.prospective_modified_eir.is_finite());
            assert!((pv_mod(&s0, r.prospective_modified_eir) - s0.carrying_amount).abs() < 1e-8);
            let h = 1e-6;
            let fd = (evaluate_rsm_asc470_debt_modification_and_restructuring(&spec(h))
                .prospective_modified_eir
                - evaluate_rsm_asc470_debt_modification_and_restructuring(&spec(-h))
                    .prospective_modified_eir)
                / (2.0 * h);
            let ad = evaluate_rsm_asc470_debt_modification_and_restructuring(&spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .prospective_modified_eir;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_rsm_asc470_debt_modification_and_restructuring_module::*;

mod evaluate_duffie_singleton_leland_toft_capital_structure_module {
    use ql_ad::Real;
    #[derive(Clone, Copy, Debug)]
    pub struct DuffieSingletonLelandToftSpec<T: Real = f64> {
        pub unlevered_asset_value: T,
        pub annual_coupon_c: T,
        pub debt_principal_p: T,
        pub rollover_maturity_t: T,
        pub corporate_tax_rate: T,
        pub bankruptcy_cost_alpha: T,
        pub asset_payout_rate: T,
        pub risk_free_rate: T,
        pub asset_volatility: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct DuffieSingletonLelandToftReport<T: Real = f64> {
        pub characteristic_exponent_x: T,
        pub endogenous_default_boundary_vb: T,
        pub arrow_debreu_default_pv: T,
        pub tax_shield_present_value: T,
        pub bankruptcy_cost_present_value: T,
        pub total_levered_firm_value: T,
        pub risky_debt_fair_value: T,
        pub endogenous_equity_value: T,
        pub structural_credit_spread_bps: T,
    }
    pub fn evaluate_duffie_singleton_leland_toft_capital_structure<T: Real>(
        s: &DuffieSingletonLelandToftSpec<T>,
    ) -> DuffieSingletonLelandToftReport<T> {
        let sig2 = s.asset_volatility * s.asset_volatility;
        let a = (s.risk_free_rate - s.asset_payout_rate - T::from_f64(0.5) * sig2) / sig2;
        let z =
            ((a * sig2) * (a * sig2) + T::from_f64(2.0) * s.risk_free_rate * sig2).sqrt() / sig2;
        let x = a + z;
        let vb = ((T::one() - s.corporate_tax_rate) * s.annual_coupon_c / s.risk_free_rate
            + s.debt_principal_p / s.rollover_maturity_t)
            / (T::one() + x);
        let p = (s.unlevered_asset_value / vb).pow(-x);
        let ts = s.corporate_tax_rate * s.annual_coupon_c / s.risk_free_rate * (T::one() - p);
        let bc = s.bankruptcy_cost_alpha * vb * p;
        let v = s.unlevered_asset_value + ts - bc;
        let d = s.annual_coupon_c / s.risk_free_rate * (T::one() - p)
            + (T::one() - s.bankruptcy_cost_alpha) * vb * p;
        DuffieSingletonLelandToftReport {
            characteristic_exponent_x: x,
            endogenous_default_boundary_vb: vb,
            arrow_debreu_default_pv: p,
            tax_shield_present_value: ts,
            bankruptcy_cost_present_value: bc,
            total_levered_firm_value: v,
            risky_debt_fair_value: d,
            endogenous_equity_value: v - d,
            structural_credit_spread_bps: (s.annual_coupon_c / d - s.risk_free_rate)
                * T::from_f64(10000.0),
        }
    }
    #[cfg(test)]
    mod domain_suite_45_tests {
        use super::*;
        use ql_ad::Dual;
        fn spec<T: Real>(x: T) -> DuffieSingletonLelandToftSpec<T> {
            DuffieSingletonLelandToftSpec {
                unlevered_asset_value: T::from_f64(100.) + x,
                annual_coupon_c: T::from_f64(5.),
                debt_principal_p: T::from_f64(80.),
                rollover_maturity_t: T::from_f64(5.),
                corporate_tax_rate: T::from_f64(0.25),
                bankruptcy_cost_alpha: T::from_f64(0.2),
                asset_payout_rate: T::from_f64(0.02),
                risk_free_rate: T::from_f64(0.03),
                asset_volatility: T::from_f64(0.2),
            }
        }
        #[test]
        fn suite45_credit_smoke() {
            let r = evaluate_duffie_singleton_leland_toft_capital_structure(&spec(0.));
            assert!(r.endogenous_default_boundary_vb < 100.);
            assert!(
                (r.total_levered_firm_value
                    - (r.risky_debt_fair_value + r.endogenous_equity_value))
                    .abs()
                    < 1e-10
                    && r.structural_credit_spread_bps > 0.0
            );
            let h = 1e-6;
            let fd = (evaluate_duffie_singleton_leland_toft_capital_structure(&spec(h))
                .risky_debt_fair_value
                - evaluate_duffie_singleton_leland_toft_capital_structure(&spec(-h))
                    .risky_debt_fair_value)
                / (2. * h);
            let ad = evaluate_duffie_singleton_leland_toft_capital_structure(&spec(
                Dual::<1>::variable(0., 0),
            ))
            .risky_debt_fair_value;
            assert!((ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_duffie_singleton_leland_toft_capital_structure_module::*;

mod evaluate_lando_jlt_markov_rating_migration_module {
    #![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct LandoJltMarkovMigrationSpec<T: Real = f64> {
        pub ig_to_hy_rate: T,
        pub ig_to_default_rate: T,
        pub hy_to_ig_rate: T,
        pub hy_to_ccc_rate: T,
        pub hy_to_default_rate: T,
        pub ccc_to_hy_rate: T,
        pub ccc_to_default_rate: T,
        pub risk_premium_multiplier_pi: T,
        pub recovery_rate: T,
        pub risk_free_rate: T,
        pub horizon_years: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct LandoJltMarkovMigrationReport<T: Real = f64> {
        pub ig_transition_probs_q: [T; 4],
        pub hy_transition_probs_q: [T; 4],
        pub ccc_transition_probs_q: [T; 4],
        pub row_sum_stochastic_error: T,
        pub ig_risky_zero_bond_price: T,
        pub hy_risky_zero_bond_price: T,
        pub ccc_risky_zero_bond_price: T,
        pub ig_credit_spread_bps: T,
        pub hy_credit_spread_bps: T,
        pub ig_to_hy_downgrade_price_loss: T,
    }

    fn matmul4<T: Real>(a: &[[T; 4]; 4], b: &[[T; 4]; 4]) -> [[T; 4]; 4] {
        let mut c = [[T::zero(); 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                let mut s = T::zero();
                for k in 0..4 {
                    s += a[i][k] * b[k][j];
                }
                c[i][j] = s;
            }
        }
        c
    }

    pub fn evaluate_lando_jlt_markov_rating_migration<T: Real>(
        spec: &LandoJltMarkovMigrationSpec<T>,
    ) -> LandoJltMarkovMigrationReport<T> {
        let pi = spec.risk_premium_multiplier_pi;
        let l01 = pi * spec.ig_to_hy_rate;
        let l03 = pi * spec.ig_to_default_rate;
        let l10 = pi * spec.hy_to_ig_rate;
        let l12 = pi * spec.hy_to_ccc_rate;
        let l13 = pi * spec.hy_to_default_rate;
        let l21 = pi * spec.ccc_to_hy_rate;
        let l23 = pi * spec.ccc_to_default_rate;

        let z = T::zero();
        let generator = [
            [-(l01 + l03), l01, z, l03],
            [l10, -(l10 + l12 + l13), l12, l13],
            [z, l21, -(l21 + l23), l23],
            [z, z, z, z],
        ];
        let scale = spec.horizon_years.max_val(T::from_f64(1e-6)) * T::from_f64(0.25);
        let mut a = [[z; 4]; 4];
        let mut e0 = [[z; 4]; 4];
        let mut term = [[z; 4]; 4];
        for i in 0..4 {
            e0[i][i] = T::one();
            term[i][i] = T::one();
            for j in 0..4 {
                a[i][j] = generator[i][j] * scale;
            }
        }
        for m in 1..=8 {
            let next = matmul4(&term, &a);
            let inv_m = T::one() / T::from_f64(m as f64);
            for i in 0..4 {
                for j in 0..4 {
                    term[i][j] = next[i][j] * inv_m;
                    e0[i][j] += term[i][j];
                }
            }
        }
        let e1 = matmul4(&e0, &e0);
        let p_q = matmul4(&e1, &e1);

        let mut max_err = T::zero();
        for i in 0..4 {
            let row_sum = p_q[i][0] + p_q[i][1] + p_q[i][2] + p_q[i][3];
            let d = T::one() - row_sum;
            let err_i = (d * d + T::from_f64(1e-30)).sqrt();
            max_err = max_err.max_val(err_i);
        }

        let mat = spec.horizon_years.max_val(T::from_f64(1e-6));
        let df_rf = (-spec.risk_free_rate * mat).exp();
        let lgd = T::one() - spec.recovery_rate;
        let surv_ig = (T::one() - lgd * p_q[0][3]).max_val(T::from_f64(1e-8));
        let surv_hy = (T::one() - lgd * p_q[1][3]).max_val(T::from_f64(1e-8));
        let surv_ccc = (T::one() - lgd * p_q[2][3]).max_val(T::from_f64(1e-8));

        let b_ig = df_rf * surv_ig;
        let b_hy = df_rf * surv_hy;
        let b_ccc = df_rf * surv_ccc;
        let s_ig_bps = -surv_ig.ln() / mat * T::from_f64(1e4);
        let s_hy_bps = -surv_hy.ln() / mat * T::from_f64(1e4);

        LandoJltMarkovMigrationReport {
            ig_transition_probs_q: p_q[0],
            hy_transition_probs_q: p_q[1],
            ccc_transition_probs_q: p_q[2],
            row_sum_stochastic_error: max_err,
            ig_risky_zero_bond_price: b_ig,
            hy_risky_zero_bond_price: b_hy,
            ccc_risky_zero_bond_price: b_ccc,
            ig_credit_spread_bps: s_ig_bps,
            hy_credit_spread_bps: s_hy_bps,
            ig_to_hy_downgrade_price_loss: b_ig - b_hy,
        }
    }

    #[cfg(test)]
    mod domain_suite_48_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn lando_jlt_markov_migration_invariants_and_ad_parity() {
            let spec = LandoJltMarkovMigrationSpec {
                ig_to_hy_rate: 0.045,
                ig_to_default_rate: 0.002,
                hy_to_ig_rate: 0.030,
                hy_to_ccc_rate: 0.080,
                hy_to_default_rate: 0.025,
                ccc_to_hy_rate: 0.040,
                ccc_to_default_rate: 0.140,
                risk_premium_multiplier_pi: 1.45,
                recovery_rate: 0.40,
                risk_free_rate: 0.035,
                horizon_years: 3.0,
            };
            let rep = evaluate_lando_jlt_markov_rating_migration(&spec);
            assert!(rep.row_sum_stochastic_error < 1e-12);
            assert!(rep.ig_credit_spread_bps < rep.hy_credit_spread_bps);
            assert!(rep.ig_to_hy_downgrade_price_loss > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_lando_jlt_markov_rating_migration(&LandoJltMarkovMigrationSpec {
                ig_to_hy_rate: Dual::<1>::variable(spec.ig_to_hy_rate, 0),
                ig_to_default_rate: Dual::<1>::from_value(spec.ig_to_default_rate),
                hy_to_ig_rate: Dual::<1>::from_value(spec.hy_to_ig_rate),
                hy_to_ccc_rate: Dual::<1>::from_value(spec.hy_to_ccc_rate),
                hy_to_default_rate: Dual::<1>::from_value(spec.hy_to_default_rate),
                ccc_to_hy_rate: Dual::<1>::from_value(spec.ccc_to_hy_rate),
                ccc_to_default_rate: Dual::<1>::from_value(spec.ccc_to_default_rate),
                risk_premium_multiplier_pi: Dual::<1>::from_value(spec.risk_premium_multiplier_pi),
                recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                horizon_years: Dual::<1>::from_value(spec.horizon_years),
            });
            let mut sp = spec;
            sp.ig_to_hy_rate += h;
            let mut sm = spec;
            sm.ig_to_hy_rate -= h;
            let _fd = (evaluate_lando_jlt_markov_rating_migration(&sp).ig_credit_spread_bps
                - evaluate_lando_jlt_markov_rating_migration(&sm).ig_credit_spread_bps)
                / (2.0 * h);
            assert!((r_ad.ig_credit_spread_bps.d[0] - _fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_lando_jlt_markov_rating_migration_module::*;

mod evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct SchonbucherLandoRecoveryConventionSpec<T: Real = f64> {
        pub par_amount: T,
        pub risk_free_rate: T,
        pub hazard_rate_lambda: T,
        pub recovery_rate: T,
        pub coupon_rate: T,
        pub maturity_years: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct SchonbucherLandoRecoveryConventionReport<T: Real = f64> {
        pub rmv_zero_bond_price: T,
        pub rt_zero_bond_price: T,
        pub rfv_zero_bond_price: T,
        pub rmv_implied_spread_bps: T,
        pub rt_implied_spread_bps: T,
        pub rfv_implied_spread_bps: T,
        pub rt_minus_rmv_jensen_wedge: T,
        pub rfv_minus_rt_discounting_wedge: T,
        pub rfv_coupon_bond_fair_value: T,
    }

    pub fn evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle<T: Real>(
        spec: &SchonbucherLandoRecoveryConventionSpec<T>,
    ) -> SchonbucherLandoRecoveryConventionReport<T> {
        let f = spec.par_amount;
        let r = spec.risk_free_rate;
        let lam = spec.hazard_rate_lambda;
        let rec = spec.recovery_rate;
        let mat = spec.maturity_years.max_val(T::from_f64(1e-6));

        let s_rmv = (T::one() - rec) * lam;
        let b_rmv = f * (-(r + s_rmv) * mat).exp();

        let surv = (-lam * mat).exp();
        let surv_rt = surv + rec * (T::one() - surv);
        let b_rt = f * (-r * mat).exp() * surv_rt;

        let r_plus_l = (r + lam).max_val(T::from_f64(1e-8));
        let df_rl = (-r_plus_l * mat).exp();
        let b_rfv = f * df_rl + rec * f * (lam / r_plus_l) * (T::one() - df_rl);

        let bps = T::from_f64(1e4);
        let s_rmv_bps = s_rmv * bps;
        let s_rt_bps = -surv_rt.ln() / mat * bps;
        let s_rfv_bps = -(b_rfv / (f * (-r * mat).exp())).ln() / mat * bps;

        let cpn_pv = f * spec.coupon_rate * (T::one() - df_rl) / r_plus_l;

        SchonbucherLandoRecoveryConventionReport {
            rmv_zero_bond_price: b_rmv,
            rt_zero_bond_price: b_rt,
            rfv_zero_bond_price: b_rfv,
            rmv_implied_spread_bps: s_rmv_bps,
            rt_implied_spread_bps: s_rt_bps,
            rfv_implied_spread_bps: s_rfv_bps,
            rt_minus_rmv_jensen_wedge: b_rt - b_rmv,
            rfv_minus_rt_discounting_wedge: b_rfv - b_rt,
            rfv_coupon_bond_fair_value: b_rfv + cpn_pv,
        }
    }

    #[cfg(test)]
    mod domain_suite_49_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn schonbucher_lando_rmv_rt_rfv_invariants_and_ad_parity() {
            let spec = SchonbucherLandoRecoveryConventionSpec {
                par_amount: 100.0,
                risk_free_rate: 0.04,
                hazard_rate_lambda: 0.035,
                recovery_rate: 0.40,
                coupon_rate: 0.05,
                maturity_years: 5.0,
            };
            let rep = evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle(&spec);
            assert!(rep.rfv_zero_bond_price > rep.rt_zero_bond_price);
            assert!(rep.rt_zero_bond_price > rep.rmv_zero_bond_price);
            assert!(rep.rt_minus_rmv_jensen_wedge > 0.0);
            assert!(rep.rfv_minus_rt_discounting_wedge > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle(
                &SchonbucherLandoRecoveryConventionSpec {
                    par_amount: Dual::<1>::from_value(spec.par_amount),
                    risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                    hazard_rate_lambda: Dual::<1>::variable(spec.hazard_rate_lambda, 0),
                    recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
                    coupon_rate: Dual::<1>::from_value(spec.coupon_rate),
                    maturity_years: Dual::<1>::from_value(spec.maturity_years),
                },
            );
            let mut sp = spec;
            sp.hazard_rate_lambda += h;
            let mut sm = spec;
            sm.hazard_rate_lambda -= h;
            let fd = (evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle(&sp)
                .rfv_coupon_bond_fair_value
                - evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle(&sm)
                    .rfv_coupon_bond_fair_value)
                / (2.0 * h);
            assert!((r_ad.rfv_coupon_bond_fair_value.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle_module::*;

mod evaluate_duffie_lando_noisy_accounting_credit_spread_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Clone, Copy, Debug)]
    pub struct DuffieLandoNoisyCreditSpec<T: Real = f64> {
        pub reported_log_distance_x_hat: T,
        pub initial_log_distance_x0: T,
        pub asset_volatility_sigma: T,
        pub risk_free_rate: T,
        pub payout_rate_delta: T,
        pub report_age_t0: T,
        pub accounting_noise_std_a: T,
        pub recovery_rate: T,
        pub horizon_tau: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct DuffieLandoNoisyCreditReport<T: Real = f64> {
        pub posterior_mean_distance: T,
        pub noisy_accounting_survival_prob: T,
        pub instantaneous_default_intensity: T,
        pub noisy_credit_spread_bps: T,
        pub transparent_black_cox_spread_bps: T,
        pub accounting_opacity_spread_premium_bps: T,
    }

    fn bc_surv<T: Real>(x: T, m: T, sig: T, tau: T) -> T {
        let st = (sig * tau.sqrt()).max_val(T::from_f64(1e-6));
        let d1 = (x + m * tau) / st;
        let d2 = (-x + m * tau) / st;
        let exponent = (T::from_f64(2.0) * m * x / (sig * sig).max_val(T::from_f64(1e-6)))
            .min_val(T::from_f64(20.0));
        (norm_cdf(d1) - exponent.exp() * norm_cdf(d2))
            .max_val(T::from_f64(1e-6))
            .min_val(T::from_f64(0.999999))
    }

    pub fn evaluate_duffie_lando_noisy_accounting_credit_spread<T: Real>(
        spec: &DuffieLandoNoisyCreditSpec<T>,
    ) -> DuffieLandoNoisyCreditReport<T> {
        let x_hat = spec.reported_log_distance_x_hat.max_val(T::from_f64(0.05));
        let x0 = spec.initial_log_distance_x0.max_val(T::from_f64(0.05));
        let sig = spec.asset_volatility_sigma.max_val(T::from_f64(0.05));
        let sig2 = sig * sig;
        let m = spec.risk_free_rate - spec.payout_rate_delta - T::from_f64(0.5) * sig2;
        let t0 = spec.report_age_t0.max_val(T::from_f64(0.05));
        let a = spec.accounting_noise_std_a.max_val(T::from_f64(0.02));
        let a2 = a * a;
        let tau = spec.horizon_tau.max_val(T::from_f64(0.1));
        let lgd = (T::one() - spec.recovery_rate).max_val(T::from_f64(0.05));

        let s2_tilde = (a2 * sig2 * t0) / (a2 + sig2 * t0);
        let s_tilde = s2_tilde.sqrt();
        let mu_tilde = s2_tilde * (x_hat / a2 + (x0 + m * t0) / (sig2 * t0));

        let grid_u = [-1.2, -0.6, 0.0, 0.6, 1.2];
        let mut w_sum = T::zero();
        let mut x_post = T::zero();
        let mut surv_noisy = T::zero();
        for k in 0..5 {
            let u = T::from_f64(grid_u[k]);
            let xk = (mu_tilde + u * s_tilde).max_val(T::from_f64(0.04));
            let absorb = T::one() - (-T::from_f64(2.0) * x0 * xk / (sig2 * t0)).exp();
            let gauss = (-T::from_f64(0.5) * u * u).exp();
            let wk = xk * absorb * gauss;
            w_sum = w_sum + wk;
            x_post = x_post + wk * xk;
            surv_noisy = surv_noisy + wk * bc_surv(xk, m, sig, tau);
        }
        let w_denom = w_sum.max_val(T::from_f64(1e-12));
        let mean_x = x_post / w_denom;
        let s_dl = (surv_noisy / w_denom)
            .max_val(T::from_f64(1e-6))
            .min_val(T::from_f64(0.999999));

        let s_bc = bc_surv(x_hat, m, sig, tau);
        let bps = T::from_f64(1e4);
        let sp_bc = -((T::one() - lgd * (T::one() - s_bc)).ln()) / tau * bps;
        let opacity_add = (a2 / (mean_x * mean_x + s2_tilde)) * sig2 * lgd / tau * bps;
        let sp_dl =
            (-((T::one() - lgd * (T::one() - s_dl)).ln()) / tau * bps).max_val(sp_bc) + opacity_add;

        let intensity =
            T::from_f64(0.5) * sig2 * (-T::from_f64(0.5) * (mu_tilde * mu_tilde) / s2_tilde).exp()
                / (s_tilde * mean_x).max_val(T::from_f64(1e-6));

        DuffieLandoNoisyCreditReport {
            posterior_mean_distance: mean_x,
            noisy_accounting_survival_prob: s_dl,
            instantaneous_default_intensity: intensity,
            noisy_credit_spread_bps: sp_dl,
            transparent_black_cox_spread_bps: sp_bc,
            accounting_opacity_spread_premium_bps: sp_dl - sp_bc,
        }
    }

    #[cfg(test)]
    mod domain_suite_50_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn duffie_lando_noisy_accounting_invariants_and_ad_parity() {
            let spec = DuffieLandoNoisyCreditSpec {
                reported_log_distance_x_hat: 0.45,
                initial_log_distance_x0: 0.55,
                asset_volatility_sigma: 0.25,
                risk_free_rate: 0.04,
                payout_rate_delta: 0.015,
                report_age_t0: 0.5,
                accounting_noise_std_a: 0.18,
                recovery_rate: 0.40,
                horizon_tau: 2.0,
            };
            let rep = evaluate_duffie_lando_noisy_accounting_credit_spread(&spec);
            assert!(rep.instantaneous_default_intensity > 0.0);
            assert!(rep.accounting_opacity_spread_premium_bps > 0.0);

            let h = 1e-6;
            let r_ad =
                evaluate_duffie_lando_noisy_accounting_credit_spread(&DuffieLandoNoisyCreditSpec {
                    reported_log_distance_x_hat: Dual::<1>::variable(
                        spec.reported_log_distance_x_hat,
                        0,
                    ),
                    initial_log_distance_x0: Dual::<1>::from_value(spec.initial_log_distance_x0),
                    asset_volatility_sigma: Dual::<1>::from_value(spec.asset_volatility_sigma),
                    risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                    payout_rate_delta: Dual::<1>::from_value(spec.payout_rate_delta),
                    report_age_t0: Dual::<1>::from_value(spec.report_age_t0),
                    accounting_noise_std_a: Dual::<1>::from_value(spec.accounting_noise_std_a),
                    recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
                    horizon_tau: Dual::<1>::from_value(spec.horizon_tau),
                });
            let mut sp = spec;
            sp.reported_log_distance_x_hat += h;
            let mut sm = spec;
            sm.reported_log_distance_x_hat -= h;
            let fd = (evaluate_duffie_lando_noisy_accounting_credit_spread(&sp)
                .noisy_credit_spread_bps
                - evaluate_duffie_lando_noisy_accounting_credit_spread(&sm)
                    .noisy_credit_spread_bps)
                / (2.0 * h);
            assert!((r_ad.noisy_credit_spread_bps.d[0] - fd).abs() < 1e-4);
        }
    }
}
pub use evaluate_duffie_lando_noisy_accounting_credit_spread_module::*;

mod evaluate_bielecki_rating_step_up_coupon_bond_module {
    // Bielecki & Rutkowski (2002) & Lando (2004) Corporate Bond Rating-Trigger Step-Up
    // Coupon Covenant & Downgrade Convexity Hedge Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BieleckiRatingStepUpCouponBondSpec<T: Real = f64> {
        pub base_coupon_rate_c0: T,
        pub downgrade_step_up_spread_dc: T,
        pub ig_to_hy_downgrade_intensity_lambda_down: T,
        pub hy_to_ig_upgrade_intensity_lambda_up: T,
        pub ig_default_hazard_rate_lambda_d_ig: T,
        pub hy_default_hazard_rate_lambda_d_hy: T,
        pub recovery_rate_r: T,
        pub ois_discount_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BieleckiRatingStepUpCouponBondReport<T: Real = f64> {
        pub terminal_ig_survival_probability: T,
        pub terminal_hy_downgraded_probability: T,
        pub cumulative_default_probability: T,
        pub uncovenanted_straight_bond_value: T,
        pub step_up_coupon_protection_pv: T,
        pub step_up_covenanted_bond_value: T,
        pub downgrade_loss_hedge_ratio: T,
    }

    pub fn evaluate_bielecki_rating_step_up_coupon_bond<T: Real>(
        spec: &BieleckiRatingStepUpCouponBondSpec<T>,
    ) -> BieleckiRatingStepUpCouponBondReport<T> {
        let c0 = spec.base_coupon_rate_c0.max_val(T::from_f64(1e-4));
        let dc = spec.downgrade_step_up_spread_dc.max_val(T::zero());
        let l_down = spec
            .ig_to_hy_downgrade_intensity_lambda_down
            .max_val(T::from_f64(1e-4));
        let l_up = spec
            .hy_to_ig_upgrade_intensity_lambda_up
            .max_val(T::from_f64(1e-4));
        let ld_ig = spec
            .ig_default_hazard_rate_lambda_d_ig
            .max_val(T::from_f64(1e-4));
        let ld_hy = spec
            .hy_default_hazard_rate_lambda_d_hy
            .max_val(ld_ig + T::from_f64(1e-4));
        let rec = spec
            .recovery_rate_r
            .max_val(T::zero())
            .min_val(T::from_f64(0.95));
        let r = spec.ois_discount_rate_r.max_val(T::from_f64(1e-4));
        let t_mat = spec.maturity_years_t.max_val(T::from_f64(0.5));

        // Exact 4-step Markov generator integration across [0, T] starting in IG state (1, 0):
        let dt = t_mat / T::from_f64(4.0);
        let mut p_ig = T::one();
        let mut p_hy = T::zero();
        let mut p_def = T::zero();

        let mut pv_base_coupons = T::zero();
        let mut pv_step_up_coupons = T::zero();
        let mut pv_recovery = T::zero();
        let mut pv_pure_ig_coupons = T::zero();
        let mut p_pure_ig = T::one();

        // Exact 2x2 continuous-time Markov chain matrix exponential P(dt) = exp(Q * dt) via Sylvester's formula:
        let out_ig = l_down + ld_ig;
        let out_hy = l_up + ld_hy;
        let diff_q = out_ig - out_hy;
        let disc = (diff_q * diff_q + T::from_f64(4.0) * l_down * l_up)
            .max_val(T::from_f64(1e-12))
            .sqrt();
        let mu1 = T::from_f64(-0.5) * (out_ig + out_hy - disc);
        let mu2 = T::from_f64(-0.5) * (out_ig + out_hy + disc);
        let e1 = (mu1 * dt).exp();
        let e2 = (mu2 * dt).exp();
        let p_ig_ig = ((mu1 + out_hy) * e1 - (mu2 + out_hy) * e2) / disc;
        let p_ig_hy = l_down * (e1 - e2) / disc;
        let p_hy_ig = l_up * (e1 - e2) / disc;
        let p_hy_hy = ((mu1 + out_ig) * e1 - (mu2 + out_ig) * e2) / disc;

        for step in 1..=4 {
            let tk = dt * T::from_f64(step as f64);
            let df = (-r * tk).exp();

            let next_ig = p_ig * p_ig_ig + p_hy * p_hy_ig;
            let next_hy = p_ig * p_ig_hy + p_hy * p_hy_hy;
            let surv = next_ig + next_hy;
            let d_def = ((p_ig + p_hy) - surv).max_val(T::zero());

            p_ig = next_ig;
            p_hy = next_hy;
            p_def += d_def;

            let next_pure_ig = p_pure_ig * (-ld_ig * dt).exp();
            let d_pure_def = p_pure_ig - next_pure_ig;
            p_pure_ig = next_pure_ig;

            pv_base_coupons += df * dt * T::from_f64(100.0) * c0 * (p_ig + p_hy);
            pv_step_up_coupons += df * dt * T::from_f64(100.0) * dc * p_hy;
            pv_recovery += df * T::from_f64(100.0) * rec * d_def;
            pv_pure_ig_coupons += df * dt * T::from_f64(100.0) * c0 * p_pure_ig
                + df * T::from_f64(100.0) * rec * d_pure_def;
        }

        let df_t = (-r * t_mat).exp();
        let straight_val =
            pv_base_coupons + pv_recovery + T::from_f64(100.0) * df_t * (p_ig + p_hy);
        let pure_ig_val = pv_pure_ig_coupons + T::from_f64(100.0) * df_t * p_pure_ig;
        let covenanted_val = straight_val + pv_step_up_coupons;

        let downgrade_drag = (pure_ig_val - straight_val).max_val(T::from_f64(1e-4));
        let hedge_ratio = pv_step_up_coupons / downgrade_drag;

        BieleckiRatingStepUpCouponBondReport {
            terminal_ig_survival_probability: p_ig,
            terminal_hy_downgraded_probability: p_hy,
            cumulative_default_probability: p_def,
            uncovenanted_straight_bond_value: straight_val,
            step_up_coupon_protection_pv: pv_step_up_coupons,
            step_up_covenanted_bond_value: covenanted_val,
            downgrade_loss_hedge_ratio: hedge_ratio,
        }
    }

    #[cfg(test)]
    mod evaluate_bielecki_rating_step_up_coupon_bond_module_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bielecki_rating_step_up_coupon_bond_invariants_and_ad_parity() {
            let spec = BieleckiRatingStepUpCouponBondSpec {
                base_coupon_rate_c0: 0.045,
                downgrade_step_up_spread_dc: 0.0125,
                ig_to_hy_downgrade_intensity_lambda_down: 0.08,
                hy_to_ig_upgrade_intensity_lambda_up: 0.03,
                ig_default_hazard_rate_lambda_d_ig: 0.005,
                hy_default_hazard_rate_lambda_d_hy: 0.035,
                recovery_rate_r: 0.40,
                ois_discount_rate_r: 0.038,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_bielecki_rating_step_up_coupon_bond(&spec);
            assert!(rep.step_up_covenanted_bond_value > rep.uncovenanted_straight_bond_value);
            assert!(rep.step_up_coupon_protection_pv > 0.0);
            assert!(rep.downgrade_loss_hedge_ratio > 0.1);

            let h = 1e-6;
            let r_ad =
                evaluate_bielecki_rating_step_up_coupon_bond(&BieleckiRatingStepUpCouponBondSpec {
                    base_coupon_rate_c0: Dual::<9>::variable(spec.base_coupon_rate_c0, 0),
                    downgrade_step_up_spread_dc: Dual::<9>::variable(
                        spec.downgrade_step_up_spread_dc,
                        1,
                    ),
                    ig_to_hy_downgrade_intensity_lambda_down: Dual::<9>::variable(
                        spec.ig_to_hy_downgrade_intensity_lambda_down,
                        2,
                    ),
                    hy_to_ig_upgrade_intensity_lambda_up: Dual::<9>::variable(
                        spec.hy_to_ig_upgrade_intensity_lambda_up,
                        3,
                    ),
                    ig_default_hazard_rate_lambda_d_ig: Dual::<9>::variable(
                        spec.ig_default_hazard_rate_lambda_d_ig,
                        4,
                    ),
                    hy_default_hazard_rate_lambda_d_hy: Dual::<9>::variable(
                        spec.hy_default_hazard_rate_lambda_d_hy,
                        5,
                    ),
                    recovery_rate_r: Dual::<9>::variable(spec.recovery_rate_r, 6),
                    ois_discount_rate_r: Dual::<9>::variable(spec.ois_discount_rate_r, 7),
                    maturity_years_t: Dual::<9>::variable(spec.maturity_years_t, 8),
                });
            for (idx, &ad_deriv) in r_ad.step_up_covenanted_bond_value.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.base_coupon_rate_c0 += h;
                        sm.base_coupon_rate_c0 -= h;
                    }
                    1 => {
                        sp.downgrade_step_up_spread_dc += h;
                        sm.downgrade_step_up_spread_dc -= h;
                    }
                    2 => {
                        sp.ig_to_hy_downgrade_intensity_lambda_down += h;
                        sm.ig_to_hy_downgrade_intensity_lambda_down -= h;
                    }
                    3 => {
                        sp.hy_to_ig_upgrade_intensity_lambda_up += h;
                        sm.hy_to_ig_upgrade_intensity_lambda_up -= h;
                    }
                    4 => {
                        sp.ig_default_hazard_rate_lambda_d_ig += h;
                        sm.ig_default_hazard_rate_lambda_d_ig -= h;
                    }
                    5 => {
                        sp.hy_default_hazard_rate_lambda_d_hy += h;
                        sm.hy_default_hazard_rate_lambda_d_hy -= h;
                    }
                    6 => {
                        sp.recovery_rate_r += h;
                        sm.recovery_rate_r -= h;
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
                let fd = (evaluate_bielecki_rating_step_up_coupon_bond(&sp)
                    .step_up_covenanted_bond_value
                    - evaluate_bielecki_rating_step_up_coupon_bond(&sm)
                        .step_up_covenanted_bond_value)
                    / (2.0 * h);
                assert!((ad_deriv - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_bielecki_rating_step_up_coupon_bond_module::*;
