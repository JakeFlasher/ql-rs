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
    let z = ((a * sig2) * (a * sig2) + T::from_f64(2.0) * s.risk_free_rate * sig2).sqrt() / sig2;
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
mod phase33_tests {
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
    fn phase33_credit_smoke() {
        let r = evaluate_duffie_singleton_leland_toft_capital_structure(&spec(0.));
        assert!(r.endogenous_default_boundary_vb < 100.);
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
