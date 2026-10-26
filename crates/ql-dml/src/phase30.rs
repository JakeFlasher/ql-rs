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
    let mut sl = T::zero();
    for i in 0..n {
        sl += mills[i];
    }
    let mx = sx / nf;
    let my = sy / nf;
    let ml = sl / nf;
    let (mut cxx, mut cll, mut cxl, mut cxy, mut cyl) =
        (T::zero(), T::zero(), T::zero(), T::zero(), T::zero());
    for i in 0..n {
        let dx = s.factor_exposures_x[i] - mx;
        let dy = s.observed_returns_y[i] - my;
        let dl = mills[i] - ml;
        cxx += dx * dx;
        cll += dl * dl;
        cxl += dx * dl;
        cxy += dx * dy;
        cyl += dl * dy;
    }
    let det = (cxx * cll - cxl * cxl).abs().max_val(T::from_f64(1e-12));
    let cbeta = (cxy * cll - cyl * cxl) / det;
    let clambda = (cyl * cxx - cxy * cxl) / det;
    let calpha = my - cbeta * mx - clambda * ml;
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
    for col in 0..3 {
        let piv = if a[col][col] >= T::zero() {
            a[col][col].max_val(T::from_f64(1e-10))
        } else {
            a[col][col].min_val(T::from_f64(-1e-10))
        };
        for j in col..4 {
            a[col][j] /= piv;
        }
        for k in 0..3 {
            if k != col {
                let f = a[k][col];
                for j in col..4 {
                    a[k][j] -= f * a[col][j];
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
        let mut eng = rg;
        for j in 0..3 {
            let smooth = (en[j] * en[j] + s.smoothing_eps * s.smoothing_eps).sqrt();
            eng[j][j] += s.l1_penalty_lambda1 / smooth.max_val(T::from_f64(1e-12));
        }
        en = solve3(eng);
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
        let smooth = (en[j] * en[j] + s.smoothing_eps * s.smoothing_eps).sqrt();
        let pen = s.l2_penalty_lambda2 + s.l1_penalty_lambda1 / smooth.max_val(T::from_f64(1e-12));
        df += gram[j][j] / (gram[j][j] + pen).max_val(T::from_f64(1e-12));
    }
    df = df.min_val(T::from_f64(3.0));
    let gcv = mse / (T::one() - df / T::from_f64(n.max(1) as f64)).powf(2.0);
    let l1_en = (en[0] * en[0] + s.smoothing_eps * s.smoothing_eps).sqrt()
        + (en[1] * en[1] + s.smoothing_eps * s.smoothing_eps).sqrt()
        + (en[2] * en[2] + s.smoothing_eps * s.smoothing_eps).sqrt();
    let l1_ols = (ols[0] * ols[0] + s.smoothing_eps * s.smoothing_eps).sqrt()
        + (ols[1] * ols[1] + s.smoothing_eps * s.smoothing_eps).sqrt()
        + (ols[2] * ols[2] + s.smoothing_eps * s.smoothing_eps).sqrt();
    let shrink = (l1_en / l1_ols.max_val(T::from_f64(1e-12))).min_val(T::one());
    HastieTibshiraniElasticNetGcvReport {
        ols_weights: ols,
        ridge_weights: ridge,
        elastic_net_weights: en,
        effective_degrees_of_freedom: df,
        in_sample_mse: mse,
        generalized_cross_validation_gcv: gcv,
        l1_norm_shrinkage_ratio: shrink,
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
    fn phase30_heckman_exact_synthetic_recovery() {
        let z = [-1.0_f64, -0.25, 0.5, 1.25, 2.0];
        let x = [0.2_f64, 0.7, 1.1, 1.8, 2.4];
        let mills: Vec<f64> = z.iter().map(|v| v.norm_pdf() / v.norm_cdf()).collect();
        let y: Vec<f64> = x
            .iter()
            .zip(mills.iter())
            .map(|(x, m)| 0.01 + 0.02 * x + 0.03 * m)
            .collect();
        let r = evaluate_greene_wooldridge_heckman_survivorship_and_gmm(
            &GreeneWooldridgeHeckmanSelectionGmmSpec {
                selection_z_scores: z.to_vec(),
                factor_exposures_x: x.to_vec(),
                observed_returns_y: y,
                gmm_moment_conditions_bar_g: [0.0; 3],
                gmm_weighting_matrix_diag: [1.0; 3],
                sample_size_n: 5.0,
            },
        );
        assert!((r.heckman_corrected_alpha - 0.01).abs() < 1e-10);
        assert!((r.heckman_corrected_beta - 0.02).abs() < 1e-10);
        assert!((r.heckman_mills_lambda_coeff - 0.03).abs() < 1e-10);
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
        let mut ep = HastieTibshiraniElasticNetGcvSpec {
            design_matrix_x: vec![
                [1.0, 0.2, 0.1],
                [1.0, 0.4, 0.5],
                [1.0, 0.6, 0.2],
                [1.0, 0.8, 0.9],
            ],
            response_y: vec![0.1, 0.2, 0.3, 0.4],
            l1_penalty_lambda1: 0.05,
            l2_penalty_lambda2: 0.1,
            smoothing_eps: 1e-4,
        };
        let e = solve_hastie_tibshirani_elastic_net_and_gcv(&ep);
        assert!(e.l1_norm_shrinkage_ratio <= 1.0 && e.effective_degrees_of_freedom <= 3.0);
        let h = 1e-6;
        ep.l1_penalty_lambda1 += h;
        let up = solve_hastie_tibshirani_elastic_net_and_gcv(&ep).elastic_net_weights[0];
        ep.l1_penalty_lambda1 -= 2.0 * h;
        let dn = solve_hastie_tibshirani_elastic_net_and_gcv(&ep).elastic_net_weights[0];
        let fd_en = (up - dn) / (2.0 * h);
        let ad_spec: HastieTibshiraniElasticNetGcvSpec<Dual<1>> =
            HastieTibshiraniElasticNetGcvSpec {
                design_matrix_x: vec![
                    [
                        Dual::from_f64(1.0),
                        Dual::from_f64(0.2),
                        Dual::from_f64(0.1),
                    ],
                    [
                        Dual::from_f64(1.0),
                        Dual::from_f64(0.4),
                        Dual::from_f64(0.5),
                    ],
                    [
                        Dual::from_f64(1.0),
                        Dual::from_f64(0.6),
                        Dual::from_f64(0.2),
                    ],
                    [
                        Dual::from_f64(1.0),
                        Dual::from_f64(0.8),
                        Dual::from_f64(0.9),
                    ],
                ],
                response_y: vec![
                    Dual::from_f64(0.1),
                    Dual::from_f64(0.2),
                    Dual::from_f64(0.3),
                    Dual::from_f64(0.4),
                ],
                l1_penalty_lambda1: Dual::<1>::variable(0.05, 0),
                l2_penalty_lambda2: Dual::from_f64(0.1),
                smoothing_eps: Dual::from_f64(1e-4),
            };
        let ad_en = solve_hastie_tibshirani_elastic_net_and_gcv(&ad_spec).elastic_net_weights[0];
        assert!(ad_en.d[0].abs() > 1e-6 && (ad_en.d[0] - fd_en).abs() < 1e-5);
    }
}
