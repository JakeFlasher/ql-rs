#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Clone, Copy, Debug)]
pub struct DuffieLandoNoisyCreditSpec<T: Real = f64> {
    pub reported_log_distance_x_hat: T,
    pub initial_log_distance_x0: T,
    pub asset_volatility_sigma: T,
    pub risk_free_rate: T,
    pub payout_rate_delta: T,
    pub report_age_t0: T,
    pub accounting_noise_std_a: T,
    pub recovery_rate: T,
    pub horizon_tau: T,
}

#[derive(Clone, Copy, Debug)]
pub struct DuffieLandoNoisyCreditReport<T: Real = f64> {
    pub posterior_mean_distance: T,
    pub noisy_accounting_survival_prob: T,
    pub instantaneous_default_intensity: T,
    pub noisy_credit_spread_bps: T,
    pub transparent_black_cox_spread_bps: T,
    pub accounting_opacity_spread_premium_bps: T,
}

fn bc_surv<T: Real>(x: T, m: T, sig: T, tau: T) -> T {
    let st = (sig * tau.sqrt()).max_val(T::from_f64(1e-6));
    let d1 = (x + m * tau) / st;
    let d2 = (-x + m * tau) / st;
    let exponent = (T::from_f64(2.0) * m * x / (sig * sig).max_val(T::from_f64(1e-6)))
        .min_val(T::from_f64(20.0));
    (norm_cdf(d1) - exponent.exp() * norm_cdf(d2))
        .max_val(T::from_f64(1e-6))
        .min_val(T::from_f64(0.999999))
}

pub fn evaluate_duffie_lando_noisy_accounting_credit_spread<T: Real>(
    spec: &DuffieLandoNoisyCreditSpec<T>,
) -> DuffieLandoNoisyCreditReport<T> {
    let x_hat = spec.reported_log_distance_x_hat.max_val(T::from_f64(0.05));
    let x0 = spec.initial_log_distance_x0.max_val(T::from_f64(0.05));
    let sig = spec.asset_volatility_sigma.max_val(T::from_f64(0.05));
    let sig2 = sig * sig;
    let m = spec.risk_free_rate - spec.payout_rate_delta - T::from_f64(0.5) * sig2;
    let t0 = spec.report_age_t0.max_val(T::from_f64(0.05));
    let a = spec.accounting_noise_std_a.max_val(T::from_f64(0.02));
    let a2 = a * a;
    let tau = spec.horizon_tau.max_val(T::from_f64(0.1));
    let lgd = (T::one() - spec.recovery_rate).max_val(T::from_f64(0.05));

    let s2_tilde = (a2 * sig2 * t0) / (a2 + sig2 * t0);
    let s_tilde = s2_tilde.sqrt();
    let mu_tilde = s2_tilde * (x_hat / a2 + (x0 + m * t0) / (sig2 * t0));

    let grid_u = [-1.2, -0.6, 0.0, 0.6, 1.2];
    let mut w_sum = T::zero();
    let mut x_post = T::zero();
    let mut surv_noisy = T::zero();
    for k in 0..5 {
        let u = T::from_f64(grid_u[k]);
        let xk = (mu_tilde + u * s_tilde).max_val(T::from_f64(0.04));
        let absorb = T::one() - (-T::from_f64(2.0) * x0 * xk / (sig2 * t0)).exp();
        let gauss = (-T::from_f64(0.5) * u * u).exp();
        let wk = xk * absorb * gauss;
        w_sum = w_sum + wk;
        x_post = x_post + wk * xk;
        surv_noisy = surv_noisy + wk * bc_surv(xk, m, sig, tau);
    }
    let w_denom = w_sum.max_val(T::from_f64(1e-12));
    let mean_x = x_post / w_denom;
    let s_dl = (surv_noisy / w_denom)
        .max_val(T::from_f64(1e-6))
        .min_val(T::from_f64(0.999999));

    let s_bc = bc_surv(x_hat, m, sig, tau);
    let bps = T::from_f64(1e4);
    let sp_bc = -((T::one() - lgd * (T::one() - s_bc)).ln()) / tau * bps;
    let opacity_add = (a2 / (mean_x * mean_x + s2_tilde)) * sig2 * lgd / tau * bps;
    let sp_dl =
        (-((T::one() - lgd * (T::one() - s_dl)).ln()) / tau * bps).max_val(sp_bc) + opacity_add;

    let intensity =
        T::from_f64(0.5) * sig2 * (-T::from_f64(0.5) * (mu_tilde * mu_tilde) / s2_tilde).exp()
            / (s_tilde * mean_x).max_val(T::from_f64(1e-6));

    DuffieLandoNoisyCreditReport {
        posterior_mean_distance: mean_x,
        noisy_accounting_survival_prob: s_dl,
        instantaneous_default_intensity: intensity,
        noisy_credit_spread_bps: sp_dl,
        transparent_black_cox_spread_bps: sp_bc,
        accounting_opacity_spread_premium_bps: sp_dl - sp_bc,
    }
}

#[cfg(test)]
mod phase38_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn duffie_lando_noisy_accounting_invariants_and_ad_parity() {
        let spec = DuffieLandoNoisyCreditSpec {
            reported_log_distance_x_hat: 0.45,
            initial_log_distance_x0: 0.55,
            asset_volatility_sigma: 0.25,
            risk_free_rate: 0.04,
            payout_rate_delta: 0.015,
            report_age_t0: 0.5,
            accounting_noise_std_a: 0.18,
            recovery_rate: 0.40,
            horizon_tau: 2.0,
        };
        let rep = evaluate_duffie_lando_noisy_accounting_credit_spread(&spec);
        assert!(rep.instantaneous_default_intensity > 0.0);
        assert!(rep.accounting_opacity_spread_premium_bps > 0.0);

        let h = 1e-6;
        let r_ad =
            evaluate_duffie_lando_noisy_accounting_credit_spread(&DuffieLandoNoisyCreditSpec {
                reported_log_distance_x_hat: Dual::<1>::variable(
                    spec.reported_log_distance_x_hat,
                    0,
                ),
                initial_log_distance_x0: Dual::<1>::from_value(spec.initial_log_distance_x0),
                asset_volatility_sigma: Dual::<1>::from_value(spec.asset_volatility_sigma),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                payout_rate_delta: Dual::<1>::from_value(spec.payout_rate_delta),
                report_age_t0: Dual::<1>::from_value(spec.report_age_t0),
                accounting_noise_std_a: Dual::<1>::from_value(spec.accounting_noise_std_a),
                recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
                horizon_tau: Dual::<1>::from_value(spec.horizon_tau),
            });
        let mut sp = spec;
        sp.reported_log_distance_x_hat += h;
        let mut sm = spec;
        sm.reported_log_distance_x_hat -= h;
        let fd = (evaluate_duffie_lando_noisy_accounting_credit_spread(&sp)
            .noisy_credit_spread_bps
            - evaluate_duffie_lando_noisy_accounting_credit_spread(&sm).noisy_credit_spread_bps)
            / (2.0 * h);
        assert!((r_ad.noisy_credit_spread_bps.d[0] - fd).abs() < 1e-4);
    }
}
