#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;
use ql_math::{inv_norm_cdf, norm_cdf};

#[derive(Clone, Copy, Debug)]
pub struct LopezDePradoDsrPurgedCvSpec<T: Real = f64> {
    pub observed_sharpe_ratio: T,
    pub sample_length_t: T,
    pub return_skewness: T,
    pub return_kurtosis: T,
    pub num_independent_trials_k: T,
    pub trials_sharpe_std: T,
    pub total_cv_samples_n: T,
    pub cv_folds_k: T,
    pub label_overlap_horizon_h: T,
    pub embargo_fraction_pct: T,
}

#[derive(Clone, Copy, Debug)]
pub struct LopezDePradoDsrPurgedCvReport<T: Real = f64> {
    pub expected_max_sharpe_hurdle_sr0: T,
    pub non_gaussian_sharpe_std_error: T,
    pub probabilistic_sharpe_ratio_zero: T,
    pub deflated_sharpe_ratio_dsr: T,
    pub non_gaussian_variance_inflation: T,
    pub unpurged_train_sample_fraction: T,
    pub purged_embargoed_train_fraction: T,
    pub purged_leakage_reduction_ratio: T,
}

pub fn evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv<T: Real>(
    spec: &LopezDePradoDsrPurgedCvSpec<T>,
) -> LopezDePradoDsrPurgedCvReport<T> {
    let gamma_em = T::from_f64(0.577_215_664_901_532_9);
    let e_val = T::one().exp();
    let k_tr = spec.num_independent_trials_k.max_val(T::from_f64(2.0));
    let p1 = (T::one() - T::one() / k_tr).min_val(T::from_f64(0.9999));
    let p2 = (T::one() - T::one() / (k_tr * e_val)).min_val(T::from_f64(0.9999));
    let sr0 = spec.trials_sharpe_std
        * ((T::one() - gamma_em) * inv_norm_cdf(p1) + gamma_em * inv_norm_cdf(p2));

    let sr = spec.observed_sharpe_ratio;
    let var_infl = (T::one() - spec.return_skewness * sr
        + T::from_f64(0.25) * (spec.return_kurtosis - T::one()) * sr * sr)
        .max_val(T::from_f64(1e-6));
    let t_m1 = (spec.sample_length_t - T::one()).max_val(T::from_f64(1.0));
    let se_sr = (var_infl / t_m1).sqrt();

    let psr0 = norm_cdf(sr / se_sr);
    let dsr = norm_cdf((sr - sr0) / se_sr);

    let k_folds = spec.cv_folds_k.max_val(T::from_f64(2.0));
    let n_obs = spec.total_cv_samples_n.max_val(T::from_f64(10.0));
    let unpurged_frac = (k_folds - T::one()) / k_folds;
    let purge_drop = T::from_f64(2.0) * k_folds * spec.label_overlap_horizon_h / n_obs;
    let purged_frac =
        (unpurged_frac - purge_drop - spec.embargo_fraction_pct).max_val(T::from_f64(0.05));
    let leak_red = (unpurged_frac - purged_frac) / unpurged_frac;

    LopezDePradoDsrPurgedCvReport {
        expected_max_sharpe_hurdle_sr0: sr0,
        non_gaussian_sharpe_std_error: se_sr,
        probabilistic_sharpe_ratio_zero: psr0,
        deflated_sharpe_ratio_dsr: dsr,
        non_gaussian_variance_inflation: var_infl,
        unpurged_train_sample_fraction: unpurged_frac,
        purged_embargoed_train_fraction: purged_frac,
        purged_leakage_reduction_ratio: leak_red,
    }
}

#[cfg(test)]
mod phase36_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn lopez_de_prado_dsr_purged_cv_invariants_and_ad_parity() {
        let spec = LopezDePradoDsrPurgedCvSpec {
            observed_sharpe_ratio: 1.45,
            sample_length_t: 60.0,
            return_skewness: -0.65,
            return_kurtosis: 4.5,
            num_independent_trials_k: 25.0,
            trials_sharpe_std: 0.45,
            total_cv_samples_n: 1000.0,
            cv_folds_k: 5.0,
            label_overlap_horizon_h: 5.0,
            embargo_fraction_pct: 0.01,
        };
        let rep = evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv(&spec);
        assert!(rep.expected_max_sharpe_hurdle_sr0 > 0.0);
        assert!(rep.deflated_sharpe_ratio_dsr < rep.probabilistic_sharpe_ratio_zero);
        assert!(rep.purged_embargoed_train_fraction < rep.unpurged_train_sample_fraction);

        let h = 1e-6;
        let r_ad =
            evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv(&LopezDePradoDsrPurgedCvSpec {
                observed_sharpe_ratio: Dual::<1>::variable(spec.observed_sharpe_ratio, 0),
                sample_length_t: Dual::<1>::from_value(spec.sample_length_t),
                return_skewness: Dual::<1>::from_value(spec.return_skewness),
                return_kurtosis: Dual::<1>::from_value(spec.return_kurtosis),
                num_independent_trials_k: Dual::<1>::from_value(spec.num_independent_trials_k),
                trials_sharpe_std: Dual::<1>::from_value(spec.trials_sharpe_std),
                total_cv_samples_n: Dual::<1>::from_value(spec.total_cv_samples_n),
                cv_folds_k: Dual::<1>::from_value(spec.cv_folds_k),
                label_overlap_horizon_h: Dual::<1>::from_value(spec.label_overlap_horizon_h),
                embargo_fraction_pct: Dual::<1>::from_value(spec.embargo_fraction_pct),
            });
        let mut sp = spec;
        sp.observed_sharpe_ratio += h;
        let mut sm = spec;
        sm.observed_sharpe_ratio -= h;
        let _fd = (evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv(&sp)
            .deflated_sharpe_ratio_dsr
            - evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv(&sm).deflated_sharpe_ratio_dsr)
            / (2.0 * h);
        assert!((r_ad.deflated_sharpe_ratio_dsr.d[0] - _fd).abs() < 1e-5);
    }
}
