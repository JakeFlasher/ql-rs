//! Gatheral (2006) & Demeterfi–Derman–Kamal–Zou (1999) Variance-Swap Log-Strip
//! Replication, Wing-Truncation Bias & Jump-to-Default Variance Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct GatheralVarSwapTruncationJumpSpec<T: Real = f64> {
    pub forward_price_f0: T,
    pub diffusive_vol_sigma: T,
    pub lower_cutoff_moneyness_kl_over_f: T,
    pub upper_cutoff_moneyness_ku_over_f: T,
    pub jump_intensity_lambda_j: T,
    pub mean_jump_log_size_mu_j: T,
    pub jump_log_vol_sigma_j: T,
    pub maturity_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct GatheralVarSwapTruncationJumpReport<T: Real = f64> {
    pub exact_diffusive_variance: T,
    pub corridor_replicated_diffusive_variance: T,
    pub wing_truncation_variance_loss: T,
    pub jump_quadratic_variation_addon: T,
    pub third_order_jump_skew_error: T,
    pub total_fair_variance_swap_strike: T,
    pub fair_volatility_swap_convexity_adjusted: T,
}

pub fn evaluate_gatheral_varswap_log_strip_tail_truncation<T: Real>(
    spec: &GatheralVarSwapTruncationJumpSpec<T>,
) -> GatheralVarSwapTruncationJumpReport<T> {
    let f0 = spec.forward_price_f0.max_val(T::from_f64(1e-4));
    let f_scale = T::from_f64(100.0) / f0;
    let sig = spec.diffusive_vol_sigma.max_val(T::from_f64(1e-4));
    let t = spec.maturity_years_t.max_val(T::from_f64(0.05));
    let std = sig * t.sqrt();

    let kl_ratio = (spec.lower_cutoff_moneyness_kl_over_f * f_scale)
        .max_val(T::from_f64(0.1))
        .min_val(T::from_f64(0.95));
    let ku_ratio = (spec.upper_cutoff_moneyness_ku_over_f * f_scale).max_val(T::from_f64(1.05));

    // Exact log-normal probability mass captured within [K_L, K_U]:
    let dl = (kl_ratio.ln() + T::from_f64(0.5) * std * std) / std;
    let du = (ku_ratio.ln() + T::from_f64(0.5) * std * std) / std;
    let captured_fraction = (norm_cdf(du) - norm_cdf(dl)).max_val(T::from_f64(1e-4));

    let var_diff = sig * sig;
    let var_corridor = var_diff * captured_fraction;
    let trunc_loss = var_diff - var_corridor;

    // Demeterfi-Derman-Kamal-Zou (1999) Eq. 29 jump contribution to variance swap:
    // E[J^2] = mu_J^2 + sigma_J^2, plus third-order log-contract vs quadratic-variation gap:
    let lam_j = spec.jump_intensity_lambda_j.max_val(T::zero());
    let muj = spec.mean_jump_log_size_mu_j;
    let sigj = spec.jump_log_vol_sigma_j.max_val(T::zero());
    let jump_qv = lam_j * (muj * muj + sigj * sigj);
    let jump_skew_err = T::from_f64(2.0)
        * lam_j
        * ((muj + T::from_f64(0.5) * sigj * sigj).exp()
            - T::one()
            - muj
            - T::from_f64(0.5) * (muj * muj + sigj * sigj));

    let fair_var = var_corridor + jump_qv + jump_skew_err;
    let vol_of_var = T::from_f64(0.5) * (jump_qv / fair_var.max_val(T::from_f64(1e-4)));
    let fair_vol_swap = fair_var.max_val(T::from_f64(1e-6)).sqrt()
        * (T::one() - T::from_f64(0.125) * vol_of_var * t);

    GatheralVarSwapTruncationJumpReport {
        exact_diffusive_variance: var_diff,
        corridor_replicated_diffusive_variance: var_corridor,
        wing_truncation_variance_loss: trunc_loss,
        jump_quadratic_variation_addon: jump_qv,
        third_order_jump_skew_error: jump_skew_err,
        total_fair_variance_swap_strike: fair_var,
        fair_volatility_swap_convexity_adjusted: fair_vol_swap,
    }
}

#[cfg(test)]
mod phase42_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn gatheral_varswap_log_strip_tail_truncation_invariants_and_ad_parity() {
        let spec = GatheralVarSwapTruncationJumpSpec {
            forward_price_f0: 100.0,
            diffusive_vol_sigma: 0.20,
            lower_cutoff_moneyness_kl_over_f: 0.65,
            upper_cutoff_moneyness_ku_over_f: 1.40,
            jump_intensity_lambda_j: 0.35,
            mean_jump_log_size_mu_j: -0.12,
            jump_log_vol_sigma_j: 0.08,
            maturity_years_t: 0.5,
        };
        let rep = evaluate_gatheral_varswap_log_strip_tail_truncation(&spec);
        assert!(rep.corridor_replicated_diffusive_variance < rep.exact_diffusive_variance);
        assert!(rep.wing_truncation_variance_loss > 0.0);
        assert!(rep.jump_quadratic_variation_addon > 0.0);

        let h = 1e-6;
        let r_ad = evaluate_gatheral_varswap_log_strip_tail_truncation(
            &GatheralVarSwapTruncationJumpSpec {
                forward_price_f0: Dual::<1>::from_value(spec.forward_price_f0),
                diffusive_vol_sigma: Dual::<1>::variable(spec.diffusive_vol_sigma, 0),
                lower_cutoff_moneyness_kl_over_f: Dual::<1>::from_value(
                    spec.lower_cutoff_moneyness_kl_over_f,
                ),
                upper_cutoff_moneyness_ku_over_f: Dual::<1>::from_value(
                    spec.upper_cutoff_moneyness_ku_over_f,
                ),
                jump_intensity_lambda_j: Dual::<1>::from_value(spec.jump_intensity_lambda_j),
                mean_jump_log_size_mu_j: Dual::<1>::from_value(spec.mean_jump_log_size_mu_j),
                jump_log_vol_sigma_j: Dual::<1>::from_value(spec.jump_log_vol_sigma_j),
                maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
            },
        );
        let mut sp = spec;
        sp.diffusive_vol_sigma += h;
        let mut sm = spec;
        sm.diffusive_vol_sigma -= h;
        let fd = (evaluate_gatheral_varswap_log_strip_tail_truncation(&sp)
            .total_fair_variance_swap_strike
            - evaluate_gatheral_varswap_log_strip_tail_truncation(&sm)
                .total_fair_variance_swap_strike)
            / (2.0 * h);
        assert!((r_ad.total_fair_variance_swap_strike.d[0] - fd).abs() < 1e-6);
    }
}
