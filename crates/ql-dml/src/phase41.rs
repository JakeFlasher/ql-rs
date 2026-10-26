//! Hastie, Tibshirani & Friedman (2009) Partial Least Squares (`PLS` Krylov Subspace)
//! vs. Principal Component Regression (`PCR`) & Ridge Shrinkage Factors Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct HtfPlsKrylovShrinkageSpec<T: Real = f64> {
    pub pc_eigenvalues_d2: [T; 3],
    pub ols_pc_projections_u: [T; 3],
    pub ridge_penalty_lambda: T,
    pub noise_variance_sigma2: T,
}

#[derive(Debug, Clone, Copy)]
pub struct HtfPlsKrylovShrinkageReport<T: Real = f64> {
    pub pls1_direction_weights: [T; 3],
    pub pls1_shrinkage_factors: [T; 3],
    pub ridge_shrinkage_factors: [T; 3],
    pub pcr1_retained_projection: T,
    pub pls1_fitted_norm_sq: T,
    pub super_unitary_expansion_factor_max: T,
}

pub fn evaluate_htf_pls_vs_pcr_krylov_shrinkage<T: Real>(
    spec: &HtfPlsKrylovShrinkageSpec<T>,
) -> HtfPlsKrylovShrinkageReport<T> {
    let d0 = spec.pc_eigenvalues_d2[0].max_val(T::from_f64(1e-4));
    let d1 = spec.pc_eigenvalues_d2[1].max_val(T::from_f64(1e-4));
    let d2 = spec.pc_eigenvalues_d2[2].max_val(T::from_f64(1e-4));
    let d = [d0, d1, d2];
    let u = spec.ols_pc_projections_u;
    let lam = spec.ridge_penalty_lambda.max_val(T::from_f64(1e-4));

    // In the principal-component basis of X^T X = diag(d_j^2), X^T y has coordinates c_j = d_j^2 u_j.
    // First PLS direction (Helland 1988; HTF Eq. 3.64) has shrinkage factors:
    // f_j^{(1)} = (sum_k d_k^2 c_k^2 / sum_k d_k^4 c_k^2) * d_j^2
    let mut num = T::zero();
    let mut den = T::zero();
    let mut c_sq = [T::zero(); 3];
    for j in 0..3 {
        let cj = d[j] * u[j];
        c_sq[j] = cj * cj;
        num += c_sq[j];
        den += d[j] * c_sq[j];
    }
    let scale = num / den.max_val(T::from_f64(1e-8));
    let pls_f = [scale * d[0], scale * d[1], scale * d[2]];
    let ridge_f = [
        d[0] / (d[0] + lam),
        d[1] / (d[1] + lam),
        d[2] / (d[2] + lam),
    ];

    let w_sum = (c_sq[0] + c_sq[1] + c_sq[2]).max_val(T::from_f64(1e-8));
    let pls_w = [c_sq[0] / w_sum, c_sq[1] / w_sum, c_sq[2] / w_sum];

    // Noise-regularized fitted response energy:
    let mut pls_norm_sq = T::zero();
    for j in 0..3 {
        let sh = pls_f[j] * u[j];
        pls_norm_sq += d[j] * sh * sh + spec.noise_variance_sigma2 * pls_f[j] * pls_f[j];
    }
    let max_f = pls_f[0].max_val(pls_f[1]).max_val(pls_f[2]);

    HtfPlsKrylovShrinkageReport {
        pls1_direction_weights: pls_w,
        pls1_shrinkage_factors: pls_f,
        ridge_shrinkage_factors: ridge_f,
        pcr1_retained_projection: u[0],
        pls1_fitted_norm_sq: pls_norm_sq,
        super_unitary_expansion_factor_max: max_f,
    }
}

#[cfg(test)]
mod phase41_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn htf_pls_vs_pcr_krylov_shrinkage_invariants_and_ad_parity() {
        let spec = HtfPlsKrylovShrinkageSpec {
            pc_eigenvalues_d2: [4.5, 1.8, 0.4],
            ols_pc_projections_u: [0.8, 1.4, 1.1],
            ridge_penalty_lambda: 0.6,
            noise_variance_sigma2: 0.05,
        };
        let rep = evaluate_htf_pls_vs_pcr_krylov_shrinkage(&spec);
        // HTF Theorem 3.3: PLS shrinkage factor on the largest eigenvalue exceeds 1.0!
        assert!(rep.super_unitary_expansion_factor_max > 1.0);
        assert!(rep.ridge_shrinkage_factors[0] < 1.0);

        let h = 1e-6;
        let r_ad = evaluate_htf_pls_vs_pcr_krylov_shrinkage(&HtfPlsKrylovShrinkageSpec {
            pc_eigenvalues_d2: [
                Dual::<1>::variable(spec.pc_eigenvalues_d2[0], 0),
                Dual::<1>::from_value(spec.pc_eigenvalues_d2[1]),
                Dual::<1>::from_value(spec.pc_eigenvalues_d2[2]),
            ],
            ols_pc_projections_u: spec.ols_pc_projections_u.map(Dual::<1>::from_value),
            ridge_penalty_lambda: Dual::<1>::from_value(spec.ridge_penalty_lambda),
            noise_variance_sigma2: Dual::<1>::from_value(spec.noise_variance_sigma2),
        });
        let mut sp = spec;
        sp.pc_eigenvalues_d2[0] += h;
        let mut sm = spec;
        sm.pc_eigenvalues_d2[0] -= h;
        let fd = (evaluate_htf_pls_vs_pcr_krylov_shrinkage(&sp).pls1_fitted_norm_sq
            - evaluate_htf_pls_vs_pcr_krylov_shrinkage(&sm).pls1_fitted_norm_sq)
            / (2.0 * h);
        assert!((r_ad.pls1_fitted_norm_sq.d[0] - fd).abs() < 1e-6);
    }
}
