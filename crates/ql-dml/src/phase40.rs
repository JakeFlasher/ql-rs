//! Tsay (2010) & Hamilton (1994) Johansen Cointegration Vector Error-Correction (`VECM`),
//! Ornstein–Uhlenbeck Half-Life & Trace Statistic Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct TsayJohansenVecmCointegrationSpec<T: Real = f64> {
    pub cointegrating_hedge_ratio_beta: T,
    pub error_correction_speed_alpha1: T,
    pub error_correction_speed_alpha2: T,
    pub current_spread_deviation_zt: T,
    pub spread_innovation_vol_sigma_eps: T,
    pub canonical_eigenvalues_lambda: [T; 2],
    pub sample_size_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct TsayJohansenVecmCointegrationReport<T: Real = f64> {
    pub spread_ar1_persistence_phi: T,
    pub ou_mean_reversion_speed_kappa: T,
    pub mean_reversion_half_life_periods: T,
    pub stationary_spread_volatility: T,
    pub expected_next_period_spread: T,
    pub johansen_trace_statistic: T,
    pub johansen_max_eigen_statistic: T,
}

pub fn evaluate_tsay_johansen_vecm_cointegration<T: Real>(
    spec: &TsayJohansenVecmCointegrationSpec<T>,
) -> TsayJohansenVecmCointegrationReport<T> {
    // Under VECM: dy_{1,t} = alpha_1 z_{t-1} + e_{1,t}, dy_{2,t} = alpha_2 z_{t-1} + e_{2,t}
    // Therefore dz_t = dy_{1,t} - beta * dy_{2,t} = (alpha_1 - beta * alpha_2) z_{t-1} + e_{z,t}
    // => z_t = phi_z * z_{t-1} + e_{z,t} with phi_z = 1 + alpha_1 - beta * alpha_2:
    let phi = (T::one() + spec.error_correction_speed_alpha1
        - spec.cointegrating_hedge_ratio_beta * spec.error_correction_speed_alpha2)
        .max_val(T::from_f64(0.05))
        .min_val(T::from_f64(0.995));
    let kappa = -phi.ln();
    let half_life = T::from_f64(std::f64::consts::LN_2) / kappa;

    let stat_vol = spec.spread_innovation_vol_sigma_eps
        / (T::one() - phi * phi).max_val(T::from_f64(1e-6)).sqrt();
    let next_z = phi * spec.current_spread_deviation_zt;

    // Johansen (1988/1991) likelihood-ratio statistics:
    let n_obs = spec.sample_size_t.max_val(T::from_f64(10.0));
    let l1 = spec.canonical_eigenvalues_lambda[0]
        .max_val(T::from_f64(1e-6))
        .min_val(T::from_f64(0.99));
    let l2 = spec.canonical_eigenvalues_lambda[1]
        .max_val(T::from_f64(1e-6))
        .min_val(T::from_f64(0.99));
    let max_eig = -n_obs * (T::one() - l1).ln();
    let trace_stat = -n_obs * ((T::one() - l1).ln() + (T::one() - l2).ln());

    TsayJohansenVecmCointegrationReport {
        spread_ar1_persistence_phi: phi,
        ou_mean_reversion_speed_kappa: kappa,
        mean_reversion_half_life_periods: half_life,
        stationary_spread_volatility: stat_vol,
        expected_next_period_spread: next_z,
        johansen_trace_statistic: trace_stat,
        johansen_max_eigen_statistic: max_eig,
    }
}

#[cfg(test)]
mod phase40_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn tsay_johansen_vecm_cointegration_invariants_and_ad_parity() {
        let spec = TsayJohansenVecmCointegrationSpec {
            cointegrating_hedge_ratio_beta: 1.15,
            error_correction_speed_alpha1: -0.18,
            error_correction_speed_alpha2: 0.08,
            current_spread_deviation_zt: 1.4,
            spread_innovation_vol_sigma_eps: 0.25,
            canonical_eigenvalues_lambda: [0.16, 0.03],
            sample_size_t: 250.0,
        };
        let rep = evaluate_tsay_johansen_vecm_cointegration(&spec);
        assert!(rep.spread_ar1_persistence_phi > 0.0 && rep.spread_ar1_persistence_phi < 1.0);
        assert!(rep.johansen_trace_statistic > rep.johansen_max_eigen_statistic);
        assert!(rep.stationary_spread_volatility > spec.spread_innovation_vol_sigma_eps);

        let h = 1e-6;
        let r_ad = evaluate_tsay_johansen_vecm_cointegration(&TsayJohansenVecmCointegrationSpec {
            cointegrating_hedge_ratio_beta: Dual::<1>::variable(
                spec.cointegrating_hedge_ratio_beta,
                0,
            ),
            error_correction_speed_alpha1: Dual::<1>::from_value(
                spec.error_correction_speed_alpha1,
            ),
            error_correction_speed_alpha2: Dual::<1>::from_value(
                spec.error_correction_speed_alpha2,
            ),
            current_spread_deviation_zt: Dual::<1>::from_value(spec.current_spread_deviation_zt),
            spread_innovation_vol_sigma_eps: Dual::<1>::from_value(
                spec.spread_innovation_vol_sigma_eps,
            ),
            canonical_eigenvalues_lambda: spec
                .canonical_eigenvalues_lambda
                .map(Dual::<1>::from_value),
            sample_size_t: Dual::<1>::from_value(spec.sample_size_t),
        });
        let mut sp = spec;
        sp.cointegrating_hedge_ratio_beta += h;
        let mut sm = spec;
        sm.cointegrating_hedge_ratio_beta -= h;
        let fd = (evaluate_tsay_johansen_vecm_cointegration(&sp).mean_reversion_half_life_periods
            - evaluate_tsay_johansen_vecm_cointegration(&sm).mean_reversion_half_life_periods)
            / (2.0 * h);
        assert!((r_ad.mean_reversion_half_life_periods.d[0] - fd).abs() < 1e-6);
    }
}
