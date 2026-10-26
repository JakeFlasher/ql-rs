use crate::Real;

#[derive(Clone, Copy, Debug)]
pub struct BeneishWahlenHkas32ForensicBifurcationSpec<T: Real = f64> {
    pub dsri: T,
    pub gmi: T,
    pub aqi: T,
    pub sgi: T,
    pub depi: T,
    pub sgai: T,
    pub lvgi: T,
    pub tata: T,
    pub net_income: T,
    pub operating_cash_flow: T,
    pub average_total_assets: T,
    pub convertible_issue_proceeds: T,
    pub convertible_face_value: T,
    pub contractual_coupon_rate: T,
    pub comparable_straight_debt_yield: T,
    pub direct_issuance_transaction_costs: T,
    pub maturity_years: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct BeneishWahlenHkas32ForensicBifurcationReport<T: Real = f64> {
    pub beneish_8var_m_score: T,
    pub manipulation_probit_probability: T,
    pub sloan_balance_sheet_accrual_ratio: T,
    pub hkas32_gross_liability_component: T,
    pub hkas32_gross_equity_reserve_component: T,
    pub hkas32_net_liability_carrying_amount: T,
    pub hkas32_net_equity_carrying_amount: T,
    pub hkas32_proceeds_conservation_error: T,
    pub year1_effective_interest_expense: T,
}
pub fn evaluate_beneish_wahlen_forensic_mscore_and_hkas32_bifurcation<T: Real>(
    s: &BeneishWahlenHkas32ForensicBifurcationSpec<T>,
) -> BeneishWahlenHkas32ForensicBifurcationReport<T> {
    let m = T::from_f64(-4.84)
        + T::from_f64(0.92) * s.dsri
        + T::from_f64(0.528) * s.gmi
        + T::from_f64(0.404) * s.aqi
        + T::from_f64(0.892) * s.sgi
        + T::from_f64(0.115) * s.depi
        - T::from_f64(0.172) * s.sgai
        + T::from_f64(4.679) * s.tata
        - T::from_f64(0.327) * s.lvgi;
    let prob = T::one() / (T::one() + (-m).exp());
    let accrual = (s.net_income - s.operating_cash_flow)
        / s.average_total_assets.abs().max_val(T::from_f64(1e-12));
    let n = s.maturity_years.max(1) as f64;
    let mut liability = T::zero();
    for i in 1..=s.maturity_years.max(1) {
        let t = T::from_f64(i as f64);
        liability += s.convertible_face_value
            * s.contractual_coupon_rate
            * (-(s.comparable_straight_debt_yield) * t).exp();
    }
    liability +=
        s.convertible_face_value * (-(s.comparable_straight_debt_yield) * T::from_f64(n)).exp();
    let equity = s.convertible_issue_proceeds - liability;
    let net_liability = liability
        - s.direct_issuance_transaction_costs * liability
            / (liability + equity).abs().max_val(T::from_f64(1e-12));
    let net_equity = equity
        - s.direct_issuance_transaction_costs * equity
            / (liability + equity).abs().max_val(T::from_f64(1e-12));
    let conservation = net_liability + net_equity
        - (s.convertible_issue_proceeds - s.direct_issuance_transaction_costs);
    BeneishWahlenHkas32ForensicBifurcationReport {
        beneish_8var_m_score: m,
        manipulation_probit_probability: prob,
        sloan_balance_sheet_accrual_ratio: accrual,
        hkas32_gross_liability_component: liability,
        hkas32_gross_equity_reserve_component: equity,
        hkas32_net_liability_carrying_amount: net_liability,
        hkas32_net_equity_carrying_amount: net_equity,
        hkas32_proceeds_conservation_error: conservation,
        year1_effective_interest_expense: net_liability * s.comparable_straight_debt_yield,
    }
}

#[cfg(test)]
mod phase30_tests {
    use super::*;
    use crate::Dual;
    fn make<T: Real>(x: T) -> BeneishWahlenHkas32ForensicBifurcationSpec<T> {
        BeneishWahlenHkas32ForensicBifurcationSpec {
            dsri: T::from_f64(1.1) + x,
            gmi: T::from_f64(1.05),
            aqi: T::from_f64(1.2),
            sgi: T::from_f64(1.15),
            depi: T::from_f64(1.03),
            sgai: T::from_f64(1.02),
            lvgi: T::from_f64(1.1),
            tata: T::from_f64(0.03),
            net_income: T::from_f64(10.0),
            operating_cash_flow: T::from_f64(8.0),
            average_total_assets: T::from_f64(200.0),
            convertible_issue_proceeds: T::from_f64(100.0),
            convertible_face_value: T::from_f64(100.0),
            contractual_coupon_rate: T::from_f64(0.03),
            comparable_straight_debt_yield: T::from_f64(0.04),
            direct_issuance_transaction_costs: T::one(),
            maturity_years: 5,
        }
    }
    #[test]
    fn phase30_beneish_hkas32_invariants_and_dual_fd() {
        let r = evaluate_beneish_wahlen_forensic_mscore_and_hkas32_bifurcation(&make(0.0));
        assert!(r.hkas32_proceeds_conservation_error.abs() < 1e-12);
        let h = 1e-6;
        let fd = (evaluate_beneish_wahlen_forensic_mscore_and_hkas32_bifurcation(&make(h))
            .beneish_8var_m_score
            - evaluate_beneish_wahlen_forensic_mscore_and_hkas32_bifurcation(&make(-h))
                .beneish_8var_m_score)
            / (2.0 * h);
        let ad = evaluate_beneish_wahlen_forensic_mscore_and_hkas32_bifurcation(&make(
            Dual::<1>::variable(0.0, 0),
        ))
        .beneish_8var_m_score;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
