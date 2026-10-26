#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[inline]
fn smax<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
}

#[inline]
fn sabs<T: Real>(z: T) -> T {
    (z * z + T::from_f64(1e-16)).sqrt()
}

#[inline]
fn soft_thresh<T: Real>(r: T, lam1: T) -> T {
    let mag = sabs(r);
    r * smax(mag - lam1, T::zero()) / mag
}

#[derive(Clone, Copy, Debug)]
pub struct HastieZouElasticNetGroupingSpec<T: Real = f64> {
    pub collinear_rho01: T,
    pub feature_target_corrs: [T; 3],
    pub lasso_penalty_lambda1: T,
    pub ridge_penalty_lambda2: T,
}

#[derive(Clone, Copy, Debug)]
pub struct HastieZouElasticNetGroupingReport<T: Real = f64> {
    pub elastic_net_betas: [T; 3],
    pub pure_ridge_betas: [T; 3],
    pub collinear_pair_beta_gap: T,
    pub zou_hastie_grouping_upper_bound: T,
    pub grouping_theorem_bound_slack: T,
    pub elastic_net_objective_value: T,
    pub coordinate_fixed_point_error: T,
}

pub fn evaluate_hastie_zou_elastic_net_and_grouping<T: Real>(
    spec: &HastieZouElasticNetGroupingSpec<T>,
) -> HastieZouElasticNetGroupingReport<T> {
    let rho = spec.collinear_rho01;
    let c = &spec.feature_target_corrs;
    let l1 = spec.lasso_penalty_lambda1.max_val(T::zero());
    let l2 = spec.ridge_penalty_lambda2.max_val(T::from_f64(1e-6));
    let denom = T::one() + l2;

    let mut b = [T::zero(); 3];
    for _ in 0..24 {
        let r0 = c[0] - rho * b[1];
        b[0] = soft_thresh(r0, l1) / denom;
        let r1 = c[1] - rho * b[0];
        b[1] = soft_thresh(r1, l1) / denom;
        let r2 = c[2];
        b[2] = soft_thresh(r2, l1) / denom;
    }

    let b0_next = soft_thresh(c[0] - rho * b[1], l1) / denom;
    let b1_next = soft_thresh(c[1] - rho * b[0], l1) / denom;
    let b2_next = soft_thresh(c[2], l1) / denom;
    let fp_err = sabs(b0_next - b[0]) + sabs(b1_next - b[1]) + sabs(b2_next - b[2]);

    let det01 = denom * denom - rho * rho;
    let ridge_b = [
        (denom * c[0] - rho * c[1]) / det01,
        (denom * c[1] - rho * c[0]) / det01,
        c[2] / denom,
    ];

    let gap = sabs(b[0] - b[1]);
    let upper = (sabs(c[0] - c[1]) + (T::one() - rho) * (sabs(b[0]) + sabs(b[1]))) / l2;
    let slack = upper - gap;

    let quad = b[0] * b[0] + b[1] * b[1] + b[2] * b[2] + T::from_f64(2.0) * rho * b[0] * b[1];
    let lin = c[0] * b[0] + c[1] * b[1] + c[2] * b[2];
    let l1_norm = sabs(b[0]) + sabs(b[1]) + sabs(b[2]);
    let l2_norm2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    let obj = T::from_f64(0.5) * (T::one() - T::from_f64(2.0) * lin + quad)
        + l1 * l1_norm
        + T::from_f64(0.5) * l2 * l2_norm2;

    HastieZouElasticNetGroupingReport {
        elastic_net_betas: b,
        pure_ridge_betas: ridge_b,
        collinear_pair_beta_gap: gap,
        zou_hastie_grouping_upper_bound: upper,
        grouping_theorem_bound_slack: slack,
        elastic_net_objective_value: obj,
        coordinate_fixed_point_error: fp_err,
    }
}

#[cfg(test)]
mod phase37_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn hastie_zou_elastic_net_grouping_invariants_and_ad_parity() {
        let spec = HastieZouElasticNetGroupingSpec {
            collinear_rho01: 0.92,
            feature_target_corrs: [0.68, 0.65, 0.42],
            lasso_penalty_lambda1: 0.08,
            ridge_penalty_lambda2: 0.25,
        };
        let rep = evaluate_hastie_zou_elastic_net_and_grouping(&spec);
        assert!(rep.coordinate_fixed_point_error < 1e-6);
        assert!(rep.grouping_theorem_bound_slack >= -1e-10);

        let h = 1e-6;
        let r_ad = evaluate_hastie_zou_elastic_net_and_grouping(&HastieZouElasticNetGroupingSpec {
            collinear_rho01: Dual::<1>::from_value(spec.collinear_rho01),
            feature_target_corrs: [
                Dual::<1>::variable(spec.feature_target_corrs[0], 0),
                Dual::<1>::from_value(spec.feature_target_corrs[1]),
                Dual::<1>::from_value(spec.feature_target_corrs[2]),
            ],
            lasso_penalty_lambda1: Dual::<1>::from_value(spec.lasso_penalty_lambda1),
            ridge_penalty_lambda2: Dual::<1>::from_value(spec.ridge_penalty_lambda2),
        });
        let mut sp = spec;
        sp.feature_target_corrs[0] += h;
        let mut sm = spec;
        sm.feature_target_corrs[0] -= h;
        let fd = (evaluate_hastie_zou_elastic_net_and_grouping(&sp).elastic_net_objective_value
            - evaluate_hastie_zou_elastic_net_and_grouping(&sm).elastic_net_objective_value)
            / (2.0 * h);
        assert!((r_ad.elastic_net_objective_value.d[0] - fd).abs() < 1e-5);
    }
}
