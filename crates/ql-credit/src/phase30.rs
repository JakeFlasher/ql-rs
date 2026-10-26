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
    let rel = s.emerging_equity_vol / s.emerging_govt_bond_vol.abs().max_val(T::from_f64(1e-12));
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
    for (cf, t) in s
        .modified_cashflows
        .iter()
        .zip(s.modified_payment_times.iter())
    {
        pv += *cf * (-(s.original_effective_interest_rate + synth) * *t).exp();
    }
    let denom = s
        .original_debt_carrying_value
        .abs()
        .max_val(T::from_f64(1e-12));
    let pct = (pv - s.original_debt_carrying_value).abs() / denom;
    let gain = s.original_debt_carrying_value - pv - s.modification_fees_paid;
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
        is_substantial_extinguishment_10pct: pct >= T::from_f64(0.10),
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
    let loss = (s.total_rwa * s.valuation_loss_shock).max_val(T::zero());
    let recap = (s.target_post_resolution_cet1_ratio * s.total_rwa - s.pre_resolution_cet1)
        .max_val(T::zero());
    let need = loss + recap;
    let cet = loss.min_val(s.pre_resolution_cet1.max_val(T::zero()));
    let rem1 = (need - cet).max_val(T::zero());
    let at1 = rem1.min_val(s.at1_principal.max_val(T::zero()));
    let rem2 = (rem1 - at1).max_val(T::zero());
    let tier2 = rem2.min_val(s.tier2_principal.max_val(T::zero()));
    let rem3 = (rem2 - tier2).max_val(T::zero());
    let senior = rem3.min_val(s.senior_tlac_principal.max_val(T::zero()));
    let senior_principal = s.senior_tlac_principal.abs().max_val(T::from_f64(1e-12));
    let recovery = (T::one() - senior / senior_principal).max_val(T::zero());
    let insolvency_recovery = (T::one() - s.liquidation_insolvency_haircut).max_val(T::zero());
    let counterfactual = s.senior_tlac_principal * insolvency_recovery;
    let resolution_recovery = s.senior_tlac_principal - senior;
    let shortfall = (counterfactual - resolution_recovery).max_val(T::zero());
    let converted = at1 + tier2 + senior;
    let ownership = converted
        / (s.existing_share_count.abs().max_val(T::from_f64(1e-12))
            * (T::one() - s.bail_in_conversion_discount).max_val(T::from_f64(1e-6))
            + converted);
    let conservation = need - (cet + at1 + tier2 + senior);
    FsbFinmaGsibResolutionBailInWaterfallReport {
        total_loss_plus_recap_need: need,
        cet1_absorbed_loss: cet,
        at1_written_down_or_converted: at1,
        tier2_bailed_in: tier2,
        senior_tlac_bailed_in: senior,
        senior_tlac_resolution_recovery_rate: recovery,
        senior_tlac_insolvency_counterfactual_recovery: counterfactual,
        ncwo_compensation_shortfall: shortfall,
        post_bail_in_creditor_ownership_share: ownership,
        waterfall_conservation_error: conservation,
    }
}

#[cfg(test)]
mod phase30_tests {
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
    fn phase30_credit_and_bailin_dual_fd() {
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
    }
}
