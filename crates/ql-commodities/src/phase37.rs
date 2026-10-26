#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct SchwartzSmithTwoFactorKalmanSpec<T: Real = f64> {
    pub prior_short_term_chi0: T,
    pub prior_equilibrium_xi0: T,
    pub mean_reversion_kappa: T,
    pub equilibrium_drift_mu: T,
    pub short_term_vol_sigma_chi: T,
    pub equilibrium_vol_sigma_xi: T,
    pub factor_correlation_rho: T,
    pub maturity_years: T,
    pub prior_var_p11: T,
    pub prior_var_p22: T,
    pub observed_log_futures: T,
    pub measurement_noise_var_r: T,
}

#[derive(Clone, Copy, Debug)]
pub struct SchwartzSmithTwoFactorKalmanReport<T: Real = f64> {
    pub prior_model_log_futures: T,
    pub prior_model_futures_price: T,
    pub integrated_log_variance: T,
    pub kalman_innovation_residual: T,
    pub kalman_gain_chi: T,
    pub kalman_gain_xi: T,
    pub posterior_short_term_chi: T,
    pub posterior_equilibrium_xi: T,
    pub posterior_fitted_log_futures: T,
    pub state_variance_trace_reduction: T,
}

pub fn evaluate_schwartz_smith_two_factor_and_kalman_update<T: Real>(
    spec: &SchwartzSmithTwoFactorKalmanSpec<T>,
) -> SchwartzSmithTwoFactorKalmanReport<T> {
    let kappa = spec.mean_reversion_kappa.max_val(T::from_f64(1e-6));
    let mat = spec.maturity_years.max_val(T::from_f64(1e-6));
    let two = T::from_f64(2.0);
    let h0 = (-kappa * mat).exp();
    let var_chi = (T::one() - (-two * kappa * mat).exp()) / (two * kappa)
        * spec.short_term_vol_sigma_chi
        * spec.short_term_vol_sigma_chi;
    let var_xi = spec.equilibrium_vol_sigma_xi * spec.equilibrium_vol_sigma_xi * mat;
    let cov_term = two
        * spec.factor_correlation_rho
        * spec.short_term_vol_sigma_chi
        * spec.equilibrium_vol_sigma_xi
        * (T::one() - h0)
        / kappa;
    let int_var = var_chi + var_xi + cov_term;
    let a_ss = spec.equilibrium_drift_mu * mat + T::from_f64(0.5) * int_var;

    let prior_log_f = h0 * spec.prior_short_term_chi0 + spec.prior_equilibrium_xi0 + a_ss;
    let prior_f = prior_log_f.exp();

    let v_inn = spec.observed_log_futures - prior_log_f;
    let s_inn = (h0 * h0 * spec.prior_var_p11 + spec.prior_var_p22 + spec.measurement_noise_var_r)
        .max_val(T::from_f64(1e-12));
    let k1 = h0 * spec.prior_var_p11 / s_inn;
    let k2 = spec.prior_var_p22 / s_inn;
    let post_chi = spec.prior_short_term_chi0 + k1 * v_inn;
    let post_xi = spec.prior_equilibrium_xi0 + k2 * v_inn;
    let post_log_f = h0 * post_chi + post_xi + a_ss;
    let trace_red = (h0 * h0 * h0 * h0 * spec.prior_var_p11 * spec.prior_var_p11
        + spec.prior_var_p22 * spec.prior_var_p22)
        / s_inn;

    SchwartzSmithTwoFactorKalmanReport {
        prior_model_log_futures: prior_log_f,
        prior_model_futures_price: prior_f,
        integrated_log_variance: int_var,
        kalman_innovation_residual: v_inn,
        kalman_gain_chi: k1,
        kalman_gain_xi: k2,
        posterior_short_term_chi: post_chi,
        posterior_equilibrium_xi: post_xi,
        posterior_fitted_log_futures: post_log_f,
        state_variance_trace_reduction: trace_red,
    }
}

#[cfg(test)]
mod phase37_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn schwartz_smith_two_factor_kalman_invariants_and_ad_parity() {
        let spec = SchwartzSmithTwoFactorKalmanSpec {
            prior_short_term_chi0: 0.08,
            prior_equilibrium_xi0: 4.20,
            mean_reversion_kappa: 1.4,
            equilibrium_drift_mu: 0.01,
            short_term_vol_sigma_chi: 0.28,
            equilibrium_vol_sigma_xi: 0.14,
            factor_correlation_rho: 0.30,
            maturity_years: 1.0,
            prior_var_p11: 0.04,
            prior_var_p22: 0.02,
            observed_log_futures: 4.32,
            measurement_noise_var_r: 0.005,
        };
        let rep = evaluate_schwartz_smith_two_factor_and_kalman_update(&spec);
        assert!(rep.state_variance_trace_reduction > 0.0);
        assert!(
            (rep.posterior_fitted_log_futures - spec.observed_log_futures).abs()
                < rep.kalman_innovation_residual.abs()
        );

        let h = 1e-6;
        let r_ad = evaluate_schwartz_smith_two_factor_and_kalman_update(
            &SchwartzSmithTwoFactorKalmanSpec {
                prior_short_term_chi0: Dual::<1>::variable(spec.prior_short_term_chi0, 0),
                prior_equilibrium_xi0: Dual::<1>::from_value(spec.prior_equilibrium_xi0),
                mean_reversion_kappa: Dual::<1>::from_value(spec.mean_reversion_kappa),
                equilibrium_drift_mu: Dual::<1>::from_value(spec.equilibrium_drift_mu),
                short_term_vol_sigma_chi: Dual::<1>::from_value(spec.short_term_vol_sigma_chi),
                equilibrium_vol_sigma_xi: Dual::<1>::from_value(spec.equilibrium_vol_sigma_xi),
                factor_correlation_rho: Dual::<1>::from_value(spec.factor_correlation_rho),
                maturity_years: Dual::<1>::from_value(spec.maturity_years),
                prior_var_p11: Dual::<1>::from_value(spec.prior_var_p11),
                prior_var_p22: Dual::<1>::from_value(spec.prior_var_p22),
                observed_log_futures: Dual::<1>::from_value(spec.observed_log_futures),
                measurement_noise_var_r: Dual::<1>::from_value(spec.measurement_noise_var_r),
            },
        );
        let mut sp = spec;
        sp.prior_short_term_chi0 += h;
        let mut sm = spec;
        sm.prior_short_term_chi0 -= h;
        let fd = (evaluate_schwartz_smith_two_factor_and_kalman_update(&sp)
            .posterior_fitted_log_futures
            - evaluate_schwartz_smith_two_factor_and_kalman_update(&sm)
                .posterior_fitted_log_futures)
            / (2.0 * h);
        assert!((r_ad.posterior_fitted_log_futures.d[0] - fd).abs() < 1e-5);
    }
}
