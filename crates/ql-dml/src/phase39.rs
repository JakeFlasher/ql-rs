//! Hastie, Tibshirani & Friedman (2009) Positive-Part James–Stein Empirical-Bayes Shrinkage
//! & Stein's Unbiased Risk Estimate (`SURE`) Dominance Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct HastieJamesSteinSureSpec<T: Real = f64> {
    pub sample_means_y: [T; 4],
    pub grand_prior_target_mu0: T,
    pub noise_variance_sigma2: T,
    pub singular_values_squared: [T; 4],
    pub ridge_penalty_lambda: T,
}

#[derive(Debug, Clone, Copy)]
pub struct HastieJamesSteinSureReport<T: Real = f64> {
    pub james_stein_shrinkage_factor_b: T,
    pub shrunk_estimates_theta_js: [T; 4],
    pub mle_expected_risk: T,
    pub stein_unbiased_risk_estimate_sure: T,
    pub sure_risk_reduction_pct: T,
    pub ridge_effective_degrees_of_freedom: T,
}

pub fn evaluate_hastie_james_stein_empirical_bayes_sure<T: Real>(
    spec: &HastieJamesSteinSureSpec<T>,
) -> HastieJamesSteinSureReport<T> {
    let p = T::from_f64(4.0);
    let p_minus_2 = T::from_f64(2.0);
    let sig2 = spec.noise_variance_sigma2.max_val(T::from_f64(1e-8));

    let mut sq_norm = T::zero();
    for i in 0..4 {
        let diff = spec.sample_means_y[i] - spec.grand_prior_target_mu0;
        sq_norm += diff * diff;
    }
    let denom = sq_norm.max_val(T::from_f64(1e-8));

    // James-Stein shrinkage factor B = (p - 2) * sigma^2 / ||y - mu_0||^2
    let b_raw = (p_minus_2 * sig2) / denom;
    // Positive-part James-Stein multiplier (1 - B)_+:
    let mult = (T::one() - b_raw).max_val(T::zero());
    let mut theta_js = [T::zero(); 4];
    for i in 0..4 {
        theta_js[i] = spec.grand_prior_target_mu0
            + mult * (spec.sample_means_y[i] - spec.grand_prior_target_mu0);
    }

    // Stein's Unbiased Risk Estimate (SURE) for p = 4 >= 3:
    // SURE = p * sigma^2 - (p - 2)^2 * sigma^4 / ||y - mu_0||^2
    let mle_risk = p * sig2;
    let sure_risk = (mle_risk - (p_minus_2 * p_minus_2 * sig2 * sig2) / denom).max_val(T::zero());
    let risk_red_pct = ((mle_risk - sure_risk) / mle_risk) * T::from_f64(100.0);

    // Effective degrees of freedom df(lambda) = sum_{j=1}^p d_j^2 / (d_j^2 + lambda)
    let lam = spec.ridge_penalty_lambda.max_val(T::zero());
    let mut df = T::zero();
    for &sv in &spec.singular_values_squared {
        let d2 = sv.max_val(T::from_f64(1e-8));
        df += d2 / (d2 + lam);
    }

    HastieJamesSteinSureReport {
        james_stein_shrinkage_factor_b: b_raw,
        shrunk_estimates_theta_js: theta_js,
        mle_expected_risk: mle_risk,
        stein_unbiased_risk_estimate_sure: sure_risk,
        sure_risk_reduction_pct: risk_red_pct,
        ridge_effective_degrees_of_freedom: df,
    }
}

#[cfg(test)]
mod phase39_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn hastie_james_stein_sure_invariants_and_ad_parity() {
        let spec = HastieJamesSteinSureSpec {
            sample_means_y: [0.08, 0.03, -0.02, 0.11],
            grand_prior_target_mu0: 0.04,
            noise_variance_sigma2: 0.0025,
            singular_values_squared: [4.2, 2.5, 1.1, 0.4],
            ridge_penalty_lambda: 0.8,
        };
        let rep = evaluate_hastie_james_stein_empirical_bayes_sure(&spec);
        // Stein's phenomenon (p = 4 >= 3): SURE risk strictly dominates MLE risk p * sigma^2
        assert!(rep.stein_unbiased_risk_estimate_sure < rep.mle_expected_risk);
        assert!(rep.sure_risk_reduction_pct > 0.0);
        assert!(
            rep.ridge_effective_degrees_of_freedom > 0.0
                && rep.ridge_effective_degrees_of_freedom < 4.0
        );

        let h = 1e-6;
        let r_ad = evaluate_hastie_james_stein_empirical_bayes_sure(&HastieJamesSteinSureSpec {
            sample_means_y: [
                Dual::<1>::variable(spec.sample_means_y[0], 0),
                Dual::<1>::from_value(spec.sample_means_y[1]),
                Dual::<1>::from_value(spec.sample_means_y[2]),
                Dual::<1>::from_value(spec.sample_means_y[3]),
            ],
            grand_prior_target_mu0: Dual::<1>::from_value(spec.grand_prior_target_mu0),
            noise_variance_sigma2: Dual::<1>::from_value(spec.noise_variance_sigma2),
            singular_values_squared: spec.singular_values_squared.map(Dual::<1>::from_value),
            ridge_penalty_lambda: Dual::<1>::from_value(spec.ridge_penalty_lambda),
        });
        let mut sp = spec;
        sp.sample_means_y[0] += h;
        let mut sm = spec;
        sm.sample_means_y[0] -= h;
        let fd = (evaluate_hastie_james_stein_empirical_bayes_sure(&sp)
            .stein_unbiased_risk_estimate_sure
            - evaluate_hastie_james_stein_empirical_bayes_sure(&sm)
                .stein_unbiased_risk_estimate_sure)
            / (2.0 * h);
        assert!((r_ad.stein_unbiased_risk_estimate_sure.d[0] - fd).abs() < 1e-6);
    }
}
