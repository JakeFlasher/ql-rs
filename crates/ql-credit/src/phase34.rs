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
