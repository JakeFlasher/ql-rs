//! Cochrane (2005) & Bansal–Yaron (2004) Long-Run Risks (`LRR`) Epstein–Zin Recursive Utility SDF,
//! Persistent Growth A_1 & Stochastic Volatility A_2 Equity Premium Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct BansalYaronLongRunRisksSpec<T: Real = f64> {
    pub relative_risk_aversion_gamma: T,
    pub elasticity_intertemporal_sub_psi: T,
    pub campbell_shiller_kappa1: T,
    pub growth_persistence_rho_x: T,
    pub volatility_persistence_nu_sigma: T,
    pub growth_shock_scale_phi_e: T,
    pub baseline_consumption_vol_sigma: T,
    pub vol_of_vol_sigma_w: T,
}

#[derive(Debug, Clone, Copy)]
pub struct BansalYaronLongRunRisksReport<T: Real = f64> {
    pub epstein_zin_theta: T,
    pub wealth_growth_loading_a1: T,
    pub wealth_volatility_loading_a2: T,
    pub market_price_of_short_run_risk_lambda_m: T,
    pub market_price_of_long_run_growth_risk_lambda_e: T,
    pub market_price_of_volatility_risk_lambda_w: T,
    pub total_consumption_claim_risk_premium: T,
}

pub fn evaluate_bansal_yaron_long_run_risks_sdf<T: Real>(
    spec: &BansalYaronLongRunRisksSpec<T>,
) -> BansalYaronLongRunRisksReport<T> {
    let gamma = spec.relative_risk_aversion_gamma;
    let psi = spec
        .elasticity_intertemporal_sub_psi
        .max_val(T::from_f64(0.1));
    let inv_psi = T::one() / psi;

    // Epstein-Zin preference parameter theta = (1 - gamma) / (1 - 1/psi):
    let theta = (T::one() - gamma) / (T::one() - inv_psi);
    let k1 = spec.campbell_shiller_kappa1;

    // Bansal-Yaron (2004, Eq. 9-10) wealth-consumption loadings:
    // A_1 = (1 - 1/psi) / (1 - kappa_1 * rho_x)
    let a1 = (T::one() - inv_psi) / (T::one() - k1 * spec.growth_persistence_rho_x);

    // A_2 = 0.5 * ((theta - theta/psi)^2 + (theta * kappa_1 * A_1 * phi_e)^2) / (theta * (1 - kappa_1 * nu_sigma))
    let term_c = theta - theta * inv_psi;
    let term_x = theta * k1 * a1 * spec.growth_shock_scale_phi_e;
    let a2 = T::from_f64(0.5) * (term_c * term_c + term_x * term_x)
        / (theta * (T::one() - k1 * spec.volatility_persistence_nu_sigma));

    // Market prices of short-run consumption shock (eta), long-run growth shock (e), and vol shock (w):
    let lam_m = gamma;
    let lam_e = (T::one() - theta) * k1 * a1 * spec.growth_shock_scale_phi_e;
    let lam_w = (T::one() - theta) * k1 * a2 * spec.vol_of_vol_sigma_w;

    // Consumption-wealth claim betas and total risk premium:
    let beta_m_e = k1 * a1 * spec.growth_shock_scale_phi_e;
    let beta_m_w = k1 * a2 * spec.vol_of_vol_sigma_w;
    let sig2 = spec.baseline_consumption_vol_sigma * spec.baseline_consumption_vol_sigma;
    let rp = (lam_m + lam_e * beta_m_e) * sig2 + lam_w * beta_m_w;

    BansalYaronLongRunRisksReport {
        epstein_zin_theta: theta,
        wealth_growth_loading_a1: a1,
        wealth_volatility_loading_a2: a2,
        market_price_of_short_run_risk_lambda_m: lam_m,
        market_price_of_long_run_growth_risk_lambda_e: lam_e,
        market_price_of_volatility_risk_lambda_w: lam_w,
        total_consumption_claim_risk_premium: rp,
    }
}

#[cfg(test)]
mod phase40_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn bansal_yaron_long_run_risks_invariants_and_ad_parity() {
        let spec = BansalYaronLongRunRisksSpec {
            relative_risk_aversion_gamma: 7.5,
            elasticity_intertemporal_sub_psi: 1.5,
            campbell_shiller_kappa1: 0.997,
            growth_persistence_rho_x: 0.979,
            volatility_persistence_nu_sigma: 0.987,
            growth_shock_scale_phi_e: 0.044,
            baseline_consumption_vol_sigma: 0.027,
            vol_of_vol_sigma_w: 0.00023,
        };
        let rep = evaluate_bansal_yaron_long_run_risks_sdf(&spec);
        // For gamma > 1 and psi > 1 (preference for early resolution of uncertainty):
        // theta < 0, A_1 > 0, A_2 < 0, and long-run growth risk commands a positive risk price!
        assert!(rep.epstein_zin_theta < 0.0);
        assert!(rep.wealth_growth_loading_a1 > 0.0);
        assert!(rep.wealth_volatility_loading_a2 < 0.0);
        assert!(rep.market_price_of_long_run_growth_risk_lambda_e > 0.0);
        assert!(rep.total_consumption_claim_risk_premium > 0.0);

        let h = 1e-6;
        let r_ad = evaluate_bansal_yaron_long_run_risks_sdf(&BansalYaronLongRunRisksSpec {
            relative_risk_aversion_gamma: Dual::<1>::variable(spec.relative_risk_aversion_gamma, 0),
            elasticity_intertemporal_sub_psi: Dual::<1>::from_value(
                spec.elasticity_intertemporal_sub_psi,
            ),
            campbell_shiller_kappa1: Dual::<1>::from_value(spec.campbell_shiller_kappa1),
            growth_persistence_rho_x: Dual::<1>::from_value(spec.growth_persistence_rho_x),
            volatility_persistence_nu_sigma: Dual::<1>::from_value(
                spec.volatility_persistence_nu_sigma,
            ),
            growth_shock_scale_phi_e: Dual::<1>::from_value(spec.growth_shock_scale_phi_e),
            baseline_consumption_vol_sigma: Dual::<1>::from_value(
                spec.baseline_consumption_vol_sigma,
            ),
            vol_of_vol_sigma_w: Dual::<1>::from_value(spec.vol_of_vol_sigma_w),
        });
        let mut sp = spec;
        sp.relative_risk_aversion_gamma += h;
        let mut sm = spec;
        sm.relative_risk_aversion_gamma -= h;
        let fd = (evaluate_bansal_yaron_long_run_risks_sdf(&sp)
            .total_consumption_claim_risk_premium
            - evaluate_bansal_yaron_long_run_risks_sdf(&sm).total_consumption_claim_risk_premium)
            / (2.0 * h);
        assert!((r_ad.total_consumption_claim_risk_premium.d[0] - fd).abs() < 1e-6);
    }
}
