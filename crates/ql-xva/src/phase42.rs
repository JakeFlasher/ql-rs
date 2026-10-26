//! Gregory (2020) & Albanese–Andersen (2014) Asymmetric Funding Valuation Adjustment
//! (`FCA`/`FBA`) & DVA Bond-CDS Funding-Spread Overlap Engine.

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
mod phase42_tests {
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
        assert!(rep.gross_funding_cost_adjustment_fca > rep.gross_funding_benefit_adjustment_fba);
        assert!(rep.dva_fba_double_count_overlap > 0.0);

        let h = 1e-6;
        let r_ad = evaluate_gregory_albanese_fva_fca_fba_overlap(&GregoryAlbaneseFvaOverlapSpec {
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
            counterparty_survival_weight: Dual::<1>::from_value(spec.counterparty_survival_weight),
            risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
            maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
        });
        let mut sp = spec;
        sp.expected_positive_exposure_epe += h;
        let mut sm = spec;
        sm.expected_positive_exposure_epe -= h;
        let fd = (evaluate_gregory_albanese_fva_fca_fba_overlap(&sp)
            .non_overlapping_net_fva_plus_dva
            - evaluate_gregory_albanese_fva_fca_fba_overlap(&sm).non_overlapping_net_fva_plus_dva)
            / (2.0 * h);
        assert!((r_ad.non_overlapping_net_fva_plus_dva.d[0] - fd).abs() < 1e-6);
    }
}
