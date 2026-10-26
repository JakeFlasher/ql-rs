#![allow(clippy::needless_range_loop)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct HastieKernelRidgeLoocvSpec<T: Real = f64> {
    pub feature_x: [T; 4],
    pub target_y: [T; 4],
    pub rbf_bandwidth_ell: T,
    pub ridge_penalty_lambda: T,
}

#[derive(Clone, Copy, Debug)]
pub struct HastieKernelRidgeLoocvReport<T: Real = f64> {
    pub dual_weights_alpha: [T; 4],
    pub fitted_values: [T; 4],
    pub hat_matrix_leverages: [T; 4],
    pub effective_degrees_of_freedom: T,
    pub in_sample_mse: T,
    pub loocv_shortcut_mse: T,
    pub loocv_bruteforce_mse: T,
    pub loocv_shortcut_identity_error: T,
    pub gcv_score: T,
}

fn solve_4x4<T: Real>(mut a: [[T; 4]; 4], mut b: [T; 4]) -> [T; 4] {
    for k in 0..4 {
        let piv = a[k][k];
        for j in k..4 {
            a[k][j] /= piv;
        }
        b[k] /= piv;
        for i in 0..4 {
            if i != k {
                let f = a[i][k];
                for j in k..4 {
                    a[i][j] -= f * a[k][j];
                }
                b[i] -= f * b[k];
            }
        }
    }
    b
}

fn solve_3x3<T: Real>(mut a: [[T; 3]; 3], mut b: [T; 3]) -> [T; 3] {
    for k in 0..3 {
        let piv = a[k][k];
        for j in k..3 {
            a[k][j] /= piv;
        }
        b[k] /= piv;
        for i in 0..3 {
            if i != k {
                let f = a[i][k];
                for j in k..3 {
                    a[i][j] -= f * a[k][j];
                }
                b[i] -= f * b[k];
            }
        }
    }
    b
}

pub fn evaluate_hastie_kernel_ridge_and_exact_loocv<T: Real>(
    s: &HastieKernelRidgeLoocvSpec<T>,
) -> HastieKernelRidgeLoocvReport<T> {
    let mut k_mat = [[T::zero(); 4]; 4];
    let mut a_mat = [[T::zero(); 4]; 4];
    let two_l2 = T::from_f64(2.0) * s.rbf_bandwidth_ell * s.rbf_bandwidth_ell;
    for i in 0..4 {
        for j in 0..4 {
            let d = s.feature_x[i] - s.feature_x[j];
            k_mat[i][j] = (-d * d / two_l2).exp();
            a_mat[i][j] = k_mat[i][j]
                + if i == j {
                    s.ridge_penalty_lambda
                } else {
                    T::zero()
                };
        }
    }
    let alpha = solve_4x4(a_mat, s.target_y);
    let mut fit = [T::zero(); 4];
    let mut lev = [T::zero(); 4];
    let mut df = T::zero();
    for i in 0..4 {
        let mut yi = T::zero();
        for j in 0..4 {
            yi += k_mat[i][j] * alpha[j];
        }
        fit[i] = yi;
        let mut ei = [T::zero(); 4];
        ei[i] = T::one();
        let col = solve_4x4(a_mat, ei);
        let mut sii = T::zero();
        for j in 0..4 {
            sii += k_mat[i][j] * col[j];
        }
        lev[i] = sii;
        df += sii;
    }
    let mut in_mse = T::zero();
    let mut loocv_short = T::zero();
    let mut loocv_brute = T::zero();
    let mut id_err_sq = T::zero();
    for i in 0..4 {
        let res = s.target_y[i] - fit[i];
        in_mse += res * res;
        let e_short = res / (T::one() - lev[i]);
        loocv_short += e_short * e_short;
        let mut sub_a = [[T::zero(); 3]; 3];
        let mut sub_y = [T::zero(); 3];
        let mut r_idx = 0;
        for r in 0..4 {
            if r != i {
                sub_y[r_idx] = s.target_y[r];
                let mut c_idx = 0;
                for c in 0..4 {
                    if c != i {
                        sub_a[r_idx][c_idx] = a_mat[r][c];
                        c_idx += 1;
                    }
                }
                r_idx += 1;
            }
        }
        let sub_alpha = solve_3x3(sub_a, sub_y);
        let mut y_hat_minus_i = T::zero();
        let mut c_idx = 0;
        for c in 0..4 {
            if c != i {
                y_hat_minus_i += k_mat[i][c] * sub_alpha[c_idx];
                c_idx += 1;
            }
        }
        let e_brute = s.target_y[i] - y_hat_minus_i;
        loocv_brute += e_brute * e_brute;
        let diff = e_brute - e_short;
        id_err_sq += diff * diff;
    }
    let four = T::from_f64(4.0);
    in_mse /= four;
    loocv_short /= four;
    loocv_brute /= four;
    let gcv_denom = T::one() - df / four;
    let gcv = in_mse / (gcv_denom * gcv_denom);
    HastieKernelRidgeLoocvReport {
        dual_weights_alpha: alpha,
        fitted_values: fit,
        hat_matrix_leverages: lev,
        effective_degrees_of_freedom: df,
        in_sample_mse: in_mse,
        loocv_shortcut_mse: loocv_short,
        loocv_bruteforce_mse: loocv_brute,
        loocv_shortcut_identity_error: (id_err_sq + T::from_f64(1e-30)).sqrt(),
        gcv_score: gcv,
    }
}

#[cfg(test)]
mod phase35_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> HastieKernelRidgeLoocvSpec<T> {
        HastieKernelRidgeLoocvSpec {
            feature_x: [
                T::from_f64(-1.0),
                T::from_f64(-0.3),
                T::from_f64(0.4),
                T::from_f64(1.1),
            ],
            target_y: [
                T::from_f64(0.8) + x,
                T::from_f64(0.2),
                T::from_f64(0.5),
                T::from_f64(1.3),
            ],
            rbf_bandwidth_ell: T::from_f64(0.8),
            ridge_penalty_lambda: T::from_f64(0.15),
        }
    }
    #[test]
    fn phase35_dml_invariants_and_dual_fd() {
        let r = evaluate_hastie_kernel_ridge_and_exact_loocv(&spec(0.0));
        assert!(r.loocv_shortcut_identity_error < 1e-10);
        assert!(r.effective_degrees_of_freedom > 0.0 && r.effective_degrees_of_freedom < 4.0);
        assert!(r.loocv_shortcut_mse >= r.in_sample_mse);
        let h = 1e-6;
        let fd = (evaluate_hastie_kernel_ridge_and_exact_loocv(&spec(h)).loocv_shortcut_mse
            - evaluate_hastie_kernel_ridge_and_exact_loocv(&spec(-h)).loocv_shortcut_mse)
            / (2.0 * h);
        let ad = evaluate_hastie_kernel_ridge_and_exact_loocv(&spec(Dual::<1>::variable(0.0, 0)))
            .loocv_shortcut_mse;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
