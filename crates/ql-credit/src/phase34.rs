use ql_ad::Real;
#[derive(Clone, Debug)]
pub struct CrepeyBurgardKjaerBilateralFvaSpec<T: Real = f64> {
    pub time_grid: Vec<T>,
    pub epe_profile: Vec<T>,
    pub ene_profile: Vec<T>,
    pub cpty_hazard_rate: T,
    pub bank_hazard_rate: T,
    pub cpty_recovery_rate: T,
    pub bank_recovery_rate: T,
    pub risk_free_rate: T,
    pub liquidity_funding_spread: T,
    pub burgard_kjaer_alpha: T,
}
#[derive(Clone, Copy, Debug)]
pub struct CrepeyBurgardKjaerBilateralFvaReport<T: Real = f64> {
    pub bilateral_cva: T,
    pub bilateral_dva: T,
    pub burgard_kjaer_effective_funding_spread: T,
    pub funding_cost_adjustment_fca: T,
    pub funding_benefit_adjustment_fba: T,
    pub dva_fca_double_counting_overlap: T,
    pub total_burgard_kjaer_xva_charge: T,
}
pub fn evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva<T: Real>(
    s: &CrepeyBurgardKjaerBilateralFvaSpec<T>,
) -> CrepeyBurgardKjaerBilateralFvaReport<T> {
    let mut c = T::zero();
    let mut d = T::zero();
    let mut f = T::zero();
    let mut b = T::zero();
    for i in 0..s
        .time_grid
        .len()
        .min(s.epe_profile.len())
        .min(s.ene_profile.len())
    {
        let t = s.time_grid[i];
        let dt = if i == 0 { t } else { t - s.time_grid[i - 1] };
        let w = dt * (-(s.risk_free_rate + s.bank_hazard_rate + s.cpty_hazard_rate) * t).exp();
        c += w * (T::one() - s.cpty_recovery_rate) * s.cpty_hazard_rate * s.epe_profile[i];
        d += w * (T::one() - s.bank_recovery_rate) * s.bank_hazard_rate * s.ene_profile[i];
        f += w
            * (s.liquidity_funding_spread
                + s.burgard_kjaer_alpha * (T::one() - s.bank_recovery_rate) * s.bank_hazard_rate)
            * s.epe_profile[i];
        b += w * s.liquidity_funding_spread * s.ene_profile[i];
    }
    let eff = s.liquidity_funding_spread
        + s.burgard_kjaer_alpha * (T::one() - s.bank_recovery_rate) * s.bank_hazard_rate;
    CrepeyBurgardKjaerBilateralFvaReport {
        bilateral_cva: c,
        bilateral_dva: d,
        burgard_kjaer_effective_funding_spread: eff,
        funding_cost_adjustment_fca: f,
        funding_benefit_adjustment_fba: b,
        dva_fca_double_counting_overlap: s.burgard_kjaer_alpha * d,
        total_burgard_kjaer_xva_charge: c - d + f - b,
    }
}

#[cfg(test)]
mod phase34_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> CrepeyBurgardKjaerBilateralFvaSpec<T> {
        CrepeyBurgardKjaerBilateralFvaSpec {
            time_grid: vec![T::from_f64(0.5), T::one(), T::from_f64(2.0)],
            epe_profile: vec![T::from_f64(10.0) + x, T::from_f64(12.0), T::from_f64(8.0)],
            ene_profile: vec![T::from_f64(6.0), T::from_f64(7.0), T::from_f64(5.0)],
            cpty_hazard_rate: T::from_f64(0.02),
            bank_hazard_rate: T::from_f64(0.015),
            cpty_recovery_rate: T::from_f64(0.4),
            bank_recovery_rate: T::from_f64(0.4),
            risk_free_rate: T::from_f64(0.03),
            liquidity_funding_spread: T::from_f64(0.008),
            burgard_kjaer_alpha: T::from_f64(0.5),
        }
    }
    #[test]
    fn phase34_credit_invariants_and_dual_fd() {
        let r = evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva(&spec(0.0));
        assert!(r.bilateral_cva > 0.0 && r.bilateral_dva > 0.0);
        let h = 1e-6;
        let fd = (evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva(&spec(h))
            .total_burgard_kjaer_xva_charge
            - evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva(&spec(-h))
                .total_burgard_kjaer_xva_charge)
            / (2.0 * h);
        let ad =
            evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva(&spec(Dual::<1>::variable(0.0, 0)))
                .total_burgard_kjaer_xva_charge;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
