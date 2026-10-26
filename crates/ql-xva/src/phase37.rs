#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct GregoryBilateralSimmVsCcpMvaSpec<T: Real = f64> {
    pub portfolio_annual_volatility: T,
    pub bilateral_mpor_years: T,
    pub ccp_mpor_years: T,
    pub exposure_decay_kappa: T,
    pub im_funding_spread: T,
    pub cpty_hazard_rate: T,
    pub cpty_lgd: T,
    pub risk_free_rate: T,
    pub ccp_default_fund_contribution: T,
    pub capital_hurdle_rate: T,
    pub maturity_years: T,
}

#[derive(Clone, Copy, Debug)]
pub struct GregoryBilateralSimmVsCcpMvaReport<T: Real = f64> {
    pub bilateral_simm_initial_margin: T,
    pub ccp_initial_margin: T,
    pub mpor_sqrt_scaling_ratio: T,
    pub bilateral_mva_charge: T,
    pub bilateral_residual_mpor_cva: T,
    pub total_bilateral_uncleared_xva: T,
    pub ccp_mva_charge: T,
    pub ccp_default_fund_kva: T,
    pub total_ccp_cleared_xva: T,
    pub central_clearing_xva_savings: T,
}

pub fn evaluate_gregory_bilateral_simm_vs_ccp_mva<T: Real>(
    spec: &GregoryBilateralSimmVsCcpMvaSpec<T>,
) -> GregoryBilateralSimmVsCcpMvaReport<T> {
    let z99 = T::from_f64(2.326_347_874_040_840_8);
    let tau_b = spec.bilateral_mpor_years.max_val(T::from_f64(1e-6));
    let tau_c = spec.ccp_mpor_years.max_val(T::from_f64(1e-6));
    let im_b = z99 * spec.portfolio_annual_volatility * tau_b.sqrt();
    let im_c = z99 * spec.portfolio_annual_volatility * tau_c.sqrt();
    let mpor_ratio = (tau_b / tau_c).sqrt();

    let tot_disc = (spec.risk_free_rate + spec.cpty_hazard_rate + spec.exposure_decay_kappa)
        .max_val(T::from_f64(1e-6));
    let mat = spec.maturity_years.max_val(T::from_f64(1e-6));
    let ann = (T::one() - (-tot_disc * mat).exp()) / tot_disc;

    let mva_b = spec.im_funding_spread * im_b * ann;
    let inv_sqrt_2pi = T::from_f64(0.398_942_280_401_432_7);
    let res_cva_b = spec.cpty_lgd
        * spec.cpty_hazard_rate
        * (spec.portfolio_annual_volatility * tau_b.sqrt() * inv_sqrt_2pi)
        * ann;
    let tot_b = mva_b + res_cva_b;

    let mva_c = spec.im_funding_spread * im_c * ann;
    let df_kva_c = spec.capital_hurdle_rate * spec.ccp_default_fund_contribution * ann;
    let tot_c = mva_c + df_kva_c;

    GregoryBilateralSimmVsCcpMvaReport {
        bilateral_simm_initial_margin: im_b,
        ccp_initial_margin: im_c,
        mpor_sqrt_scaling_ratio: mpor_ratio,
        bilateral_mva_charge: mva_b,
        bilateral_residual_mpor_cva: res_cva_b,
        total_bilateral_uncleared_xva: tot_b,
        ccp_mva_charge: mva_c,
        ccp_default_fund_kva: df_kva_c,
        total_ccp_cleared_xva: tot_c,
        central_clearing_xva_savings: tot_b - tot_c,
    }
}

#[cfg(test)]
mod phase37_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn gregory_bilateral_simm_vs_ccp_mva_invariants_and_ad_parity() {
        let spec = GregoryBilateralSimmVsCcpMvaSpec {
            portfolio_annual_volatility: 15.0,
            bilateral_mpor_years: 10.0 / 250.0,
            ccp_mpor_years: 5.0 / 250.0,
            exposure_decay_kappa: 0.12,
            im_funding_spread: 0.0085,
            cpty_hazard_rate: 0.015,
            cpty_lgd: 0.60,
            risk_free_rate: 0.035,
            ccp_default_fund_contribution: 0.35,
            capital_hurdle_rate: 0.10,
            maturity_years: 5.0,
        };
        let rep = evaluate_gregory_bilateral_simm_vs_ccp_mva(&spec);
        assert!(rep.bilateral_simm_initial_margin > rep.ccp_initial_margin);
        assert!(rep.bilateral_mva_charge > rep.ccp_mva_charge);

        let h = 1e-6;
        let r_ad = evaluate_gregory_bilateral_simm_vs_ccp_mva(&GregoryBilateralSimmVsCcpMvaSpec {
            portfolio_annual_volatility: Dual::<1>::variable(spec.portfolio_annual_volatility, 0),
            bilateral_mpor_years: Dual::<1>::from_value(spec.bilateral_mpor_years),
            ccp_mpor_years: Dual::<1>::from_value(spec.ccp_mpor_years),
            exposure_decay_kappa: Dual::<1>::from_value(spec.exposure_decay_kappa),
            im_funding_spread: Dual::<1>::from_value(spec.im_funding_spread),
            cpty_hazard_rate: Dual::<1>::from_value(spec.cpty_hazard_rate),
            cpty_lgd: Dual::<1>::from_value(spec.cpty_lgd),
            risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
            ccp_default_fund_contribution: Dual::<1>::from_value(
                spec.ccp_default_fund_contribution,
            ),
            capital_hurdle_rate: Dual::<1>::from_value(spec.capital_hurdle_rate),
            maturity_years: Dual::<1>::from_value(spec.maturity_years),
        });
        let mut sp = spec;
        sp.portfolio_annual_volatility += h;
        let mut sm = spec;
        sm.portfolio_annual_volatility -= h;
        let fd = (evaluate_gregory_bilateral_simm_vs_ccp_mva(&sp).central_clearing_xva_savings
            - evaluate_gregory_bilateral_simm_vs_ccp_mva(&sm).central_clearing_xva_savings)
            / (2.0 * h);
        assert!((r_ad.central_clearing_xva_savings.d[0] - fd).abs() < 1e-5);
    }
}
