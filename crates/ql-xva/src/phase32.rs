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
mod phase32_tests {
    use super::*;
    #[test]
    fn wwr_positive() {
        let s = TurlakovGreenKenyonWwrKvaSpec {
            time_grid: vec![1.0],
            independent_epe_profile: vec![1.0],
            pfe_addon_profile: vec![1.0],
            base_cpty_hazard: 0.02,
            hazard_log_vol: 0.3,
            exposure_vol: 0.2,
            wwr_correlation_rho: 0.5,
            lgd: 0.6,
            risk_free_rate: 0.02,
            bank_hazard_rate: 0.01,
            ccr_risk_weight: 0.08,
            shareholder_cost_of_capital: 0.1,
            capital_return_offset_beta: 0.5,
        };
        let r = evaluate_turlakov_green_kenyon_wwr_and_kva(&s);
        assert!(r.turlakov_wwr_cva > r.independent_cva);
    }
}
