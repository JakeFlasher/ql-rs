#![allow(clippy::excessive_precision)]
use ql_ad::Real;
use ql_math::{inv_norm_cdf, norm_cdf};

#[derive(Clone, Copy, Debug)]
pub struct LandoSchonbucherBasketDefaultSpec<T: Real = f64> {
    pub obligor_hazard_rates: [T; 4],
    pub factor_loading_beta: T,
    pub recovery_rate: T,
    pub risk_free_rate: T,
    pub maturity_years: T,
}

#[derive(Clone, Copy, Debug)]
pub struct LandoSchonbucherBasketDefaultReport<T: Real = f64> {
    pub unconditional_default_count_probs: [T; 5],
    pub loss_distribution_sum_error: T,
    pub first_to_default_prob: T,
    pub second_to_default_prob: T,
    pub first_to_default_spread_bps: T,
    pub second_to_default_spread_bps: T,
    pub correlation_diversification_ratio: T,
}

pub fn evaluate_lando_schonbucher_kth_to_default_basket<T: Real>(
    spec: &LandoSchonbucherBasketDefaultSpec<T>,
) -> LandoSchonbucherBasketDefaultReport<T> {
    let nodes = [
        -2.8569700138728056,
        -1.3556261799742659,
        0.0,
        1.3556261799742659,
        2.8569700138728056,
    ];
    let weights = [
        0.011257411327720689,
        0.22207592200561264,
        0.5333333333333333,
        0.22207592200561264,
        0.011257411327720689,
    ];
    let t = spec.maturity_years.max_val(T::from_f64(1e-6));
    let beta = spec.factor_loading_beta;
    let denom = (T::one() - beta * beta).max_val(T::from_f64(1e-8)).sqrt();
    let mut c = [T::zero(); 4];
    let mut sum_lambda = T::zero();
    for i in 0..4 {
        let p_i = (T::one() - (-spec.obligor_hazard_rates[i] * t).exp())
            .max_val(T::from_f64(1e-8))
            .min_val(T::from_f64(1.0 - 1e-8));
        c[i] = inv_norm_cdf(p_i);
        sum_lambda += spec.obligor_hazard_rates[i];
    }
    let mut probs = [T::zero(); 5];
    for q_idx in 0..5 {
        let m = T::from_f64(nodes[q_idx]);
        let w = T::from_f64(weights[q_idx]);
        let mut cond = [T::zero(); 5];
        cond[0] = T::one();
        for i in 0..4 {
            let pi_m = norm_cdf((c[i] - beta * m) / denom);
            let mut next = [T::zero(); 5];
            for k in 0..=i {
                next[k] += cond[k] * (T::one() - pi_m);
                next[k + 1] += cond[k] * pi_m;
            }
            cond = next;
        }
        for k in 0..5 {
            probs[k] += w * cond[k];
        }
    }
    let total_p = probs[0] + probs[1] + probs[2] + probs[3] + probs[4];
    let diff = T::one() - total_p;
    let err = (diff * diff + T::from_f64(1e-30)).sqrt();
    let q1 = probs[1] + probs[2] + probs[3] + probs[4];
    let q2 = probs[2] + probs[3] + probs[4];
    let lam1 = -(T::one() - q1).max_val(T::from_f64(1e-12)).ln() / t;
    let lam2 = -(T::one() - q2).max_val(T::from_f64(1e-12)).ln() / t;
    let df_adj = T::one() + T::from_f64(0.05) * spec.risk_free_rate * t;
    let lgd = T::one() - spec.recovery_rate;
    let s1 = lgd * lam1 * df_adj * T::from_f64(10000.0);
    let s2 = lgd * lam2 * df_adj * T::from_f64(10000.0);
    LandoSchonbucherBasketDefaultReport {
        unconditional_default_count_probs: probs,
        loss_distribution_sum_error: err,
        first_to_default_prob: q1,
        second_to_default_prob: q2,
        first_to_default_spread_bps: s1,
        second_to_default_spread_bps: s2,
        correlation_diversification_ratio: lam1 / sum_lambda.max_val(T::from_f64(1e-8)),
    }
}

#[cfg(test)]
mod phase35_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> LandoSchonbucherBasketDefaultSpec<T> {
        LandoSchonbucherBasketDefaultSpec {
            obligor_hazard_rates: [
                T::from_f64(0.015) + x,
                T::from_f64(0.02),
                T::from_f64(0.025),
                T::from_f64(0.03),
            ],
            factor_loading_beta: T::from_f64(0.5),
            recovery_rate: T::from_f64(0.4),
            risk_free_rate: T::from_f64(0.03),
            maturity_years: T::from_f64(5.0),
        }
    }
    #[test]
    fn phase35_credit_invariants_and_dual_fd() {
        let r = evaluate_lando_schonbucher_kth_to_default_basket(&spec(0.0));
        assert!(r.loss_distribution_sum_error < 1e-12);
        assert!(r.first_to_default_spread_bps > r.second_to_default_spread_bps);
        let h = 1e-6;
        let fd = (evaluate_lando_schonbucher_kth_to_default_basket(&spec(h))
            .first_to_default_spread_bps
            - evaluate_lando_schonbucher_kth_to_default_basket(&spec(-h))
                .first_to_default_spread_bps)
            / (2.0 * h);
        let ad =
            evaluate_lando_schonbucher_kth_to_default_basket(&spec(Dual::<1>::variable(0.0, 0)))
                .first_to_default_spread_bps;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-4);
    }
}
