//! Cochrane (2005) & Campbell (2018) Campbell–Cochrane (1999) External Habit-Formation
//! Surplus-Consumption SDF, Countercyclical Sharpe Ratio & Constant Risk-Free Rate Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct CampbellCochraneExternalHabitSpec<T: Real = f64> {
    pub utility_curvature_gamma: T,
    pub consumption_growth_mean_gc: T,
    pub consumption_growth_vol_sigmac: T,
    pub log_surplus_persistence_phis: T,
    pub subjective_discount_rate_rho_time: T,
    pub current_log_surplus_deviation_ds: T,
}

#[derive(Debug, Clone, Copy)]
pub struct CampbellCochraneExternalHabitReport<T: Real = f64> {
    pub steady_state_surplus_ratio_sbar: T,
    pub max_log_surplus_smax: T,
    pub sensitivity_function_lambda_s: T,
    pub local_effective_risk_aversion: T,
    pub real_risk_free_rate_rf: T,
    pub conditional_max_sharpe_ratio: T,
    pub conditional_equity_risk_premium: T,
}

pub fn evaluate_campbell_cochrane_external_habit_sdf<T: Real>(
    spec: &CampbellCochraneExternalHabitSpec<T>,
) -> CampbellCochraneExternalHabitReport<T> {
    let gamma = spec.utility_curvature_gamma.max_val(T::from_f64(0.5));
    let sig_c = spec
        .consumption_growth_vol_sigmac
        .max_val(T::from_f64(1e-4));
    let one_minus_phi = (T::one() - spec.log_surplus_persistence_phis).max_val(T::from_f64(1e-4));

    // Steady-state surplus consumption ratio S_bar = sigma_c * sqrt(gamma / (1 - phi_s))
    let s_bar = sig_c * (gamma / one_minus_phi).sqrt();
    let s_bar_log = s_bar.ln();
    let s_max = s_bar_log + T::from_f64(0.5) * (T::one() - s_bar * s_bar);

    // Sensitivity function lambda(s_t) = (1 / S_bar) * sqrt(1 - 2 * (s_t - s_bar)) - 1
    let inside = (T::one() - T::from_f64(2.0) * spec.current_log_surplus_deviation_ds)
        .max_val(T::from_f64(1e-6));
    let lambda_s = inside.sqrt() / s_bar - T::one();

    // Local effective relative risk aversion gamma_t = gamma * (1 + lambda(s_t))
    let local_rra = gamma * (T::one() + lambda_s);

    // Exact constant real risk-free rate under Campbell-Cochrane baseline (B = 0):
    // r_f = rho_time + gamma * g_c - 0.5 * gamma * (1 - phi_s)
    let rf = spec.subjective_discount_rate_rho_time + gamma * spec.consumption_growth_mean_gc
        - T::from_f64(0.5) * gamma * one_minus_phi;

    // Conditional maximum Sharpe ratio SR_t = gamma * sigma_c * (1 + lambda(s_t))
    let sr_t = local_rra * sig_c;
    // Conditional equity risk premium for claim with return volatility ~ sigma_c * (1 + lambda(s_t)):
    let erp_t = sr_t * sig_c * (T::one() + lambda_s);

    CampbellCochraneExternalHabitReport {
        steady_state_surplus_ratio_sbar: s_bar,
        max_log_surplus_smax: s_max,
        sensitivity_function_lambda_s: lambda_s,
        local_effective_risk_aversion: local_rra,
        real_risk_free_rate_rf: rf,
        conditional_max_sharpe_ratio: sr_t,
        conditional_equity_risk_premium: erp_t,
    }
}

#[cfg(test)]
mod phase39_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn campbell_cochrane_external_habit_invariants_and_ad_parity() {
        let boom = CampbellCochraneExternalHabitSpec {
            utility_curvature_gamma: 2.0,
            consumption_growth_mean_gc: 0.0189,
            consumption_growth_vol_sigmac: 0.015,
            log_surplus_persistence_phis: 0.87,
            subjective_discount_rate_rho_time: 0.0116,
            current_log_surplus_deviation_ds: 0.05,
        };
        let mut recession = boom;
        recession.current_log_surplus_deviation_ds = -0.15;

        let r_boom = evaluate_campbell_cochrane_external_habit_sdf(&boom);
        let r_rec = evaluate_campbell_cochrane_external_habit_sdf(&recession);
        // Countercyclical risk premium: Sharpe ratio and ERP rise in recessions (ds < 0)
        // while real risk-free rate remains invariant to surplus consumption ds!
        assert!(r_rec.conditional_max_sharpe_ratio > r_boom.conditional_max_sharpe_ratio);
        assert!(r_rec.conditional_equity_risk_premium > r_boom.conditional_equity_risk_premium);
        assert!((r_rec.real_risk_free_rate_rf - r_boom.real_risk_free_rate_rf).abs() < 1e-12);

        let h = 1e-6;
        let r_ad =
            evaluate_campbell_cochrane_external_habit_sdf(&CampbellCochraneExternalHabitSpec {
                utility_curvature_gamma: Dual::<1>::from_value(recession.utility_curvature_gamma),
                consumption_growth_mean_gc: Dual::<1>::from_value(
                    recession.consumption_growth_mean_gc,
                ),
                consumption_growth_vol_sigmac: Dual::<1>::from_value(
                    recession.consumption_growth_vol_sigmac,
                ),
                log_surplus_persistence_phis: Dual::<1>::from_value(
                    recession.log_surplus_persistence_phis,
                ),
                subjective_discount_rate_rho_time: Dual::<1>::from_value(
                    recession.subjective_discount_rate_rho_time,
                ),
                current_log_surplus_deviation_ds: Dual::<1>::variable(
                    recession.current_log_surplus_deviation_ds,
                    0,
                ),
            });
        let mut sp = recession;
        sp.current_log_surplus_deviation_ds += h;
        let mut sm = recession;
        sm.current_log_surplus_deviation_ds -= h;
        let fd = (evaluate_campbell_cochrane_external_habit_sdf(&sp).conditional_max_sharpe_ratio
            - evaluate_campbell_cochrane_external_habit_sdf(&sm).conditional_max_sharpe_ratio)
            / (2.0 * h);
        assert!((r_ad.conditional_max_sharpe_ratio.d[0] - fd).abs() < 1e-6);
    }
}
