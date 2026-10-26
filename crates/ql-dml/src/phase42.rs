//! Tsay (2010) & Teräsvirta (1994) Exponential Smooth-Transition Threshold
//! Autoregression (`ESTAR`) No-Trade Band & Regime Half-Life Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct TsayEstarThresholdAutoregressionSpec<T: Real = f64> {
    pub current_mispricing_z: T,
    pub transaction_cost_threshold_c: T,
    pub transition_smoothness_gamma: T,
    pub inner_band_reversion_speed_kappa_in: T,
    pub outer_arbitrage_reversion_speed_kappa_out: T,
    pub innovation_volatility_sigma_eps: T,
}

#[derive(Debug, Clone, Copy)]
pub struct TsayEstarThresholdAutoregressionReport<T: Real = f64> {
    pub transition_activation_weight_g: T,
    pub effective_local_reversion_speed_kappa: T,
    pub expected_one_step_spread_change: T,
    pub expected_next_period_mispricing: T,
    pub local_regime_half_life_periods: T,
    pub inner_vs_outer_half_life_ratio: T,
    pub local_stationary_variance_proxy: T,
}

pub fn evaluate_tsay_star_threshold_autoregression_band<T: Real>(
    spec: &TsayEstarThresholdAutoregressionSpec<T>,
) -> TsayEstarThresholdAutoregressionReport<T> {
    let z = spec.current_mispricing_z;
    let c = spec.transaction_cost_threshold_c.max_val(T::from_f64(1e-4));
    let gam = spec.transition_smoothness_gamma.max_val(T::from_f64(0.05));
    let k_in = spec
        .inner_band_reversion_speed_kappa_in
        .max_val(T::from_f64(0.01));
    let k_out = spec.outer_arbitrage_reversion_speed_kappa_out.max_val(k_in);

    // Smooth excess squared deviation above transaction-cost band c^2:
    let d_sq = z * z - c * c;
    let excess_sq = T::from_f64(0.5) * (d_sq + (d_sq * d_sq + T::from_f64(0.01)).sqrt());

    // ESTAR transition function G(z; gamma, c) = 1 - exp(-gamma * excess_sq) in (0, 1):
    let g = T::one() - (-gam * excess_sq).exp();
    let k_eff = k_in * (T::one() - g) + k_out * g;

    let dz = -k_eff * z;
    let z_next = z + dz;
    let phi_eff = (T::one() - k_eff)
        .max_val(T::from_f64(0.05))
        .min_val(T::from_f64(0.995));
    let half_life = -T::from_f64(std::f64::consts::LN_2) / phi_eff.ln();

    let phi_in = (T::one() - k_in)
        .max_val(T::from_f64(0.05))
        .min_val(T::from_f64(0.995));
    let phi_out = (T::one() - k_out)
        .max_val(T::from_f64(0.05))
        .min_val(T::from_f64(0.995));
    let hl_ratio = phi_out.ln() / phi_in.ln();

    let stat_var = spec.innovation_volatility_sigma_eps * spec.innovation_volatility_sigma_eps
        / (T::one() - phi_eff * phi_eff).max_val(T::from_f64(1e-4));

    TsayEstarThresholdAutoregressionReport {
        transition_activation_weight_g: g,
        effective_local_reversion_speed_kappa: k_eff,
        expected_one_step_spread_change: dz,
        expected_next_period_mispricing: z_next,
        local_regime_half_life_periods: half_life,
        inner_vs_outer_half_life_ratio: hl_ratio,
        local_stationary_variance_proxy: stat_var,
    }
}

#[cfg(test)]
mod phase42_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn tsay_star_threshold_autoregression_invariants_and_ad_parity() {
        let spec = TsayEstarThresholdAutoregressionSpec {
            current_mispricing_z: 1.85,
            transaction_cost_threshold_c: 0.90,
            transition_smoothness_gamma: 0.85,
            inner_band_reversion_speed_kappa_in: 0.04,
            outer_arbitrage_reversion_speed_kappa_out: 0.38,
            innovation_volatility_sigma_eps: 0.22,
        };
        let rep = evaluate_tsay_star_threshold_autoregression_band(&spec);
        assert!(rep.transition_activation_weight_g > 0.5);
        assert!(
            rep.effective_local_reversion_speed_kappa > spec.inner_band_reversion_speed_kappa_in
        );
        assert!(rep.inner_vs_outer_half_life_ratio > 2.0);

        let h = 1e-6;
        let r_ad = evaluate_tsay_star_threshold_autoregression_band(
            &TsayEstarThresholdAutoregressionSpec {
                current_mispricing_z: Dual::<1>::variable(spec.current_mispricing_z, 0),
                transaction_cost_threshold_c: Dual::<1>::from_value(
                    spec.transaction_cost_threshold_c,
                ),
                transition_smoothness_gamma: Dual::<1>::from_value(
                    spec.transition_smoothness_gamma,
                ),
                inner_band_reversion_speed_kappa_in: Dual::<1>::from_value(
                    spec.inner_band_reversion_speed_kappa_in,
                ),
                outer_arbitrage_reversion_speed_kappa_out: Dual::<1>::from_value(
                    spec.outer_arbitrage_reversion_speed_kappa_out,
                ),
                innovation_volatility_sigma_eps: Dual::<1>::from_value(
                    spec.innovation_volatility_sigma_eps,
                ),
            },
        );
        let mut sp = spec;
        sp.current_mispricing_z += h;
        let mut sm = spec;
        sm.current_mispricing_z -= h;
        let fd = (evaluate_tsay_star_threshold_autoregression_band(&sp)
            .local_regime_half_life_periods
            - evaluate_tsay_star_threshold_autoregression_band(&sm).local_regime_half_life_periods)
            / (2.0 * h);
        assert!((r_ad.local_regime_half_life_periods.d[0] - fd).abs() < 1e-5);
    }
}
