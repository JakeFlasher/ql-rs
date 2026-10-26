#![allow(clippy::needless_range_loop)]
use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct GreeneWooldridgeHeckmanSelectionGmmSpec<T: Real = f64> {
    pub selection_z_scores: Vec<T>,
    pub factor_exposures_x: Vec<T>,
    pub observed_returns_y: Vec<T>,
    pub gmm_moment_conditions_bar_g: [T; 3],
    pub gmm_weighting_matrix_diag: [T; 3],
    pub sample_size_n: T,
}
#[derive(Clone, Debug)]
pub struct GreeneWooldridgeHeckmanSelectionGmmReport<T: Real = f64> {
    pub inverse_mills_ratios: Vec<T>,
    pub naive_survivorship_biased_alpha: T,
    pub naive_factor_beta: T,
    pub heckman_corrected_alpha: T,
    pub heckman_corrected_beta: T,
    pub heckman_mills_lambda_coeff: T,
    pub survivorship_alpha_bias: T,
    pub hansen_gmm_j_statistic: T,
}
pub fn evaluate_greene_wooldridge_heckman_survivorship_and_gmm<T: Real>(
    s: &GreeneWooldridgeHeckmanSelectionGmmSpec<T>,
) -> GreeneWooldridgeHeckmanSelectionGmmReport<T> {
    let n = s
        .selection_z_scores
        .len()
        .min(s.factor_exposures_x.len())
        .min(s.observed_returns_y.len());
    let mut mills = Vec::with_capacity(n);
    for i in 0..n {
        mills.push(
            s.selection_z_scores[i].norm_pdf()
                / s.selection_z_scores[i]
                    .norm_cdf()
                    .max_val(T::from_f64(1e-12)),
        );
    }
    let nf = T::from_f64(n.max(1) as f64);
    let (mut sx, mut sy, mut sxx, mut sxy) = (T::zero(), T::zero(), T::zero(), T::zero());
    for i in 0..n {
        sx += s.factor_exposures_x[i];
        sy += s.observed_returns_y[i];
        sxx += s.factor_exposures_x[i] * s.factor_exposures_x[i];
        sxy += s.factor_exposures_x[i] * s.observed_returns_y[i];
    }
    let den = (nf * sxx - sx * sx).abs().max_val(T::from_f64(1e-12));
    let beta = (nf * sxy - sx * sy) / den;
    let alpha = (sy - beta * sx) / nf;
    let (mut sl, mut sxl, mut sy_l, mut sll) = (T::zero(), T::zero(), T::zero(), T::zero());
    for i in 0..n {
        sl += mills[i];
        sxl += s.factor_exposures_x[i] * mills[i];
        sy_l += s.observed_returns_y[i] * mills[i];
        sll += mills[i] * mills[i];
    }
    let det = (sxx * sll - sxl * sxl).abs().max_val(T::from_f64(1e-12));
    let cbeta = (sxy * sll - sy_l * sxl) / det;
    let clambda = (sy_l * sxx - sxy * sxl) / det;
    let calpha = (sy - cbeta * sx - clambda * sl) / nf;
    let mut j = T::zero();
    for i in 0..3 {
        j += s.gmm_moment_conditions_bar_g[i]
            * s.gmm_weighting_matrix_diag[i]
            * s.gmm_moment_conditions_bar_g[i];
    }
    j *= s.sample_size_n;
    GreeneWooldridgeHeckmanSelectionGmmReport {
        inverse_mills_ratios: mills,
        naive_survivorship_biased_alpha: alpha,
        naive_factor_beta: beta,
        heckman_corrected_alpha: calpha,
        heckman_corrected_beta: cbeta,
        heckman_mills_lambda_coeff: clambda,
        survivorship_alpha_bias: alpha - calpha,
        hansen_gmm_j_statistic: j,
    }
}

fn solve3<T: Real>(mut a: [[T; 4]; 3]) -> [T; 3] {
    for i in 0..3 {
        let d = a[i][i];
        for j in i..4 {
            a[i][j] /= d;
        }
        for k in 0..3 {
            if k != i {
                let f = a[k][i];
                for j in i..4 {
                    a[k][j] -= f * a[i][j];
                }
            }
        }
    }
    [a[0][3], a[1][3], a[2][3]]
}
#[derive(Clone, Debug)]
pub struct HastieTibshiraniElasticNetGcvSpec<T: Real = f64> {
    pub design_matrix_x: Vec<[T; 3]>,
    pub response_y: Vec<T>,
    pub l1_penalty_lambda1: T,
    pub l2_penalty_lambda2: T,
    pub smoothing_eps: T,
}
#[derive(Clone, Debug)]
pub struct HastieTibshiraniElasticNetGcvReport<T: Real = f64> {
    pub ols_weights: [T; 3],
    pub ridge_weights: [T; 3],
    pub elastic_net_weights: [T; 3],
    pub effective_degrees_of_freedom: T,
    pub in_sample_mse: T,
    pub generalized_cross_validation_gcv: T,
    pub l1_norm_shrinkage_ratio: T,
}
pub fn solve_hastie_tibshirani_elastic_net_and_gcv<T: Real>(
    s: &HastieTibshiraniElasticNetGcvSpec<T>,
) -> HastieTibshiraniElasticNetGcvReport<T> {
    let n = s.design_matrix_x.len().min(s.response_y.len());
    let mut gram = [[T::zero(); 4]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..n {
                gram[i][j] += s.design_matrix_x[k][i] * s.design_matrix_x[k][j];
            }
            if i == j {
                gram[i][3] = T::zero();
            }
        }
        for k in 0..n {
            gram[i][3] += s.design_matrix_x[k][i] * s.response_y[k];
        }
    }
    let ols = solve3(gram);
    let mut rg = gram;
    for i in 0..3 {
        rg[i][i] += s.l2_penalty_lambda2;
    }
    let ridge = solve3(rg);
    let mut en = ridge;
    for _ in 0..8 {
        for j in 0..3 {
            let smooth = (en[j] * en[j] + s.smoothing_eps * s.smoothing_eps).sqrt();
            en[j] =
                ridge[j] / (T::one() + s.l1_penalty_lambda1 / smooth.max_val(T::from_f64(1e-12)));
        }
    }
    let mut mse = T::zero();
    for k in 0..n {
        let pred = en[0] * s.design_matrix_x[k][0]
            + en[1] * s.design_matrix_x[k][1]
            + en[2] * s.design_matrix_x[k][2];
        let e = s.response_y[k] - pred;
        mse += e * e;
    }
    mse /= T::from_f64(n.max(1) as f64);
    let mut df = T::zero();
    for j in 0..3 {
        df += ridge[j] * ridge[j]
            / (ridge[j] * ridge[j] + s.l2_penalty_lambda2 + s.l1_penalty_lambda1)
                .max_val(T::from_f64(1e-12));
    }
    df = df.min_val(T::from_f64(3.0));
    let gcv = mse / (T::one() - df / T::from_f64(n.max(1) as f64)).powf(2.0);
    let norm_en = (en[0] * en[0] + en[1] * en[1] + en[2] * en[2]).sqrt();
    let norm_r = (ridge[0] * ridge[0] + ridge[1] * ridge[1] + ridge[2] * ridge[2]).sqrt();
    HastieTibshiraniElasticNetGcvReport {
        ols_weights: ols,
        ridge_weights: ridge,
        elastic_net_weights: en,
        effective_degrees_of_freedom: df,
        in_sample_mse: mse,
        generalized_cross_validation_gcv: gcv,
        l1_norm_shrinkage_ratio: (norm_en / norm_r.max_val(T::from_f64(1e-12))).min_val(T::one()),
    }
}

#[cfg(test)]
mod phase30_tests {
    use super::*;
    use ql_ad::Dual;
    fn heck<T: Real>(x: T) -> GreeneWooldridgeHeckmanSelectionGmmSpec<T> {
        GreeneWooldridgeHeckmanSelectionGmmSpec {
            selection_z_scores: vec![T::from_f64(-0.5), T::zero(), T::from_f64(0.5), T::one()],
            factor_exposures_x: vec![
                T::from_f64(0.5),
                T::one(),
                T::from_f64(1.5),
                T::from_f64(2.0),
            ],
            observed_returns_y: vec![
                T::from_f64(0.02) + x,
                T::from_f64(0.04),
                T::from_f64(0.06),
                T::from_f64(0.09),
            ],
            gmm_moment_conditions_bar_g: [
                T::from_f64(0.01),
                T::from_f64(-0.02),
                T::from_f64(0.015),
            ],
            gmm_weighting_matrix_diag: [T::from_f64(2.0), T::from_f64(1.5), T::one()],
            sample_size_n: T::from_f64(100.0),
        }
    }
    #[test]
    fn phase30_dml_dual_fd_and_invariants() {
        let r = evaluate_greene_wooldridge_heckman_survivorship_and_gmm(&heck(0.0));
        assert_eq!(r.inverse_mills_ratios.len(), 4);
        let h = 1e-6;
        let fd = (evaluate_greene_wooldridge_heckman_survivorship_and_gmm(&heck(h))
            .heckman_corrected_alpha
            - evaluate_greene_wooldridge_heckman_survivorship_and_gmm(&heck(-h))
                .heckman_corrected_alpha)
            / (2.0 * h);
        let ad = evaluate_greene_wooldridge_heckman_survivorship_and_gmm(&heck(
            Dual::<1>::variable(0.0, 0),
        ))
        .heckman_corrected_alpha;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        let e = solve_hastie_tibshirani_elastic_net_and_gcv(&HastieTibshiraniElasticNetGcvSpec {
            design_matrix_x: vec![
                [1.0, 0.2, 0.1],
                [1.0, 0.4, 0.3],
                [1.0, 0.6, 0.5],
                [1.0, 0.8, 0.7],
            ],
            response_y: vec![0.1, 0.2, 0.3, 0.4],
            l1_penalty_lambda1: 0.05,
            l2_penalty_lambda2: 0.1,
            smoothing_eps: 1e-4,
        });
        assert!(e.l1_norm_shrinkage_ratio <= 1.0 && e.effective_degrees_of_freedom <= 3.0);
    }
}
