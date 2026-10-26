#![allow(clippy::needless_range_loop)]
use nalgebra::{DMatrix, DVector};
use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct AgrawalBoydDifferentiableQpSpec<T: Real = f64> {
    pub q_matrix: Vec<Vec<T>>,
    pub q_linear: Vec<T>,
    pub a_eq: Vec<Vec<T>>,
    pub b_eq: Vec<T>,
    pub g_ineq: Vec<Vec<T>>,
    pub h_ineq: Vec<T>,
    pub barrier_tau: T,
}
#[derive(Clone, Debug)]
pub struct AgrawalBoydDifferentiableQpReport<T: Real = f64> {
    pub optimal_primal_z: Vec<T>,
    pub optimal_ineq_dual_lambda: Vec<T>,
    pub optimal_eq_dual_nu: Vec<T>,
    pub optimal_objective_value: T,
    pub kkt_stationarity_max_error: T,
    pub primal_eq_feasibility_error: T,
    pub reverse_kkt_grad_q_linear: Vec<T>,
    pub reverse_kkt_grad_b_eq: Vec<T>,
    pub reverse_kkt_grad_h_ineq: Vec<T>,
}
fn solve_kkt_linear<T: Real>(mat: &[Vec<T>], rhs: &[T]) -> Vec<T> {
    let n = rhs.len();
    let mut a = vec![vec![T::zero(); n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = mat[i][j];
        }
        a[i][n] = rhs[i];
    }
    for i in 0..n {
        let mut p = i;
        for r in (i + 1)..n {
            if a[r][i].abs().val() > a[p][i].abs().val() {
                p = r;
            }
        }
        a.swap(i, p);
        let piv = if a[i][i].abs().val() < 1e-14 {
            T::from_f64(1e-14)
        } else {
            a[i][i]
        };
        for j in i..=n {
            a[i][j] /= piv;
        }
        let pivot_row = a[i].clone();
        for r in 0..n {
            if r == i {
                continue;
            }
            let f = a[r][i];
            for j in i..=n {
                a[r][j] -= f * pivot_row[j];
            }
        }
    }
    (0..n).map(|i| a[i][n]).collect()
}

pub fn solve_agrawal_boyd_differentiable_qp_layer<T: Real>(
    s: &AgrawalBoydDifferentiableQpSpec<T>,
) -> AgrawalBoydDifferentiableQpReport<T> {
    let n = s.q_linear.len();
    let m = s.g_ineq.len();
    let p = s.a_eq.len();
    let mut z = vec![T::zero(); n];
    if p > 0 {
        let mut aa = vec![vec![T::zero(); p]; p];
        for i in 0..p {
            for j in 0..p {
                aa[i][j] = (0..n).fold(T::zero(), |a, k| {
                    a + s.a_eq[i].get(k).copied().unwrap_or(T::zero())
                        * s.a_eq[j].get(k).copied().unwrap_or(T::zero())
                });
            }
        }
        let y = solve_kkt_linear(&aa, &s.b_eq);
        for k in 0..n {
            z[k] = (0..p).fold(T::zero(), |a, i| {
                a + s.a_eq[i].get(k).copied().unwrap_or(T::zero()) * y[i]
            });
        }
    }
    let mut lambda = vec![T::zero(); m];
    for j in 0..m {
        let gz = s.g_ineq[j]
            .iter()
            .zip(z.iter())
            .map(|(&a, &b)| a * b)
            .fold(T::zero(), |a, b| a + b);
        lambda[j] = s.barrier_tau
            / (s.h_ineq.get(j).copied().unwrap_or(T::zero()) - gz).max_val(T::from_f64(0.5));
    }
    let mut nu = vec![T::zero(); p];
    for _ in 0..25 {
        let dim = n + m + p;
        let mut mat = vec![vec![T::zero(); dim]; dim];
        let mut rhs = vec![T::zero(); dim];
        for i in 0..n {
            for k in 0..n {
                mat[i][k] = s
                    .q_matrix
                    .get(i)
                    .and_then(|r| r.get(k))
                    .copied()
                    .unwrap_or(T::zero());
            }
            for j in 0..m {
                mat[i][n + j] = s.g_ineq[j].get(i).copied().unwrap_or(T::zero());
            }
            for l in 0..p {
                mat[i][n + m + l] = s.a_eq[l].get(i).copied().unwrap_or(T::zero());
            }
            let mut r = s.q_linear[i];
            for k in 0..n {
                r += mat[i][k] * z[k];
            }
            for j in 0..m {
                r += mat[i][n + j] * lambda[j];
            }
            for l in 0..p {
                r += mat[i][n + m + l] * nu[l];
            }
            rhs[i] = -r;
        }
        for j in 0..m {
            let row = n + j;
            let gz = s.g_ineq[j]
                .iter()
                .zip(z.iter())
                .map(|(&a, &b)| a * b)
                .fold(T::zero(), |a, b| a + b);
            let slack =
                (s.h_ineq.get(j).copied().unwrap_or(T::zero()) - gz).max_val(T::from_f64(1e-8));
            for k in 0..n {
                mat[row][k] = lambda[j] * s.g_ineq[j].get(k).copied().unwrap_or(T::zero());
            }
            mat[row][n + j] = -slack;
            rhs[row] =
                -(lambda[j] * (gz - s.h_ineq.get(j).copied().unwrap_or(T::zero())) + s.barrier_tau);
        }
        for l in 0..p {
            let row = n + m + l;
            for k in 0..n {
                mat[row][k] = s.a_eq[l].get(k).copied().unwrap_or(T::zero());
            }
            let az = s.a_eq[l]
                .iter()
                .zip(z.iter())
                .map(|(&a, &b)| a * b)
                .fold(T::zero(), |a, b| a + b);
            rhs[row] = -(az - s.b_eq.get(l).copied().unwrap_or(T::zero()));
        }
        let step = solve_kkt_linear(&mat, &rhs);
        let mut alpha = T::one();
        for j in 0..m {
            if step[n + j].val() < 0.0 {
                alpha = alpha.min_val(T::from_f64(0.99) * (-lambda[j] / step[n + j]));
            }
            let gdz = s.g_ineq[j]
                .iter()
                .enumerate()
                .map(|(k, &g)| g * step[k])
                .fold(T::zero(), |a, b| a + b);
            if gdz.val() > 0.0 {
                let slack = (s.h_ineq[j]
                    - s.g_ineq[j]
                        .iter()
                        .zip(z.iter())
                        .map(|(&g, &x)| g * x)
                        .fold(T::zero(), |a, b| a + b))
                .max_val(T::from_f64(1e-8));
                alpha = alpha.min_val(T::from_f64(0.99) * slack / gdz);
            }
        }
        alpha = alpha.max_val(T::from_f64(1e-4));
        for k in 0..n {
            z[k] += alpha * step[k];
        }
        for j in 0..m {
            lambda[j] += alpha * step[n + j];
        }
        for l in 0..p {
            nu[l] += alpha * step[n + m + l];
        }
    }
    let mut obj = T::zero();
    for i in 0..n {
        for j in 0..n {
            obj += T::from_f64(0.5)
                * z[i]
                * s.q_matrix
                    .get(i)
                    .and_then(|r| r.get(j))
                    .copied()
                    .unwrap_or(T::zero())
                * z[j];
        }
        obj += s.q_linear[i] * z[i];
    }
    let mut stationarity = T::zero();
    for i in 0..n {
        let mut r = s.q_linear[i];
        for j in 0..n {
            r += s
                .q_matrix
                .get(i)
                .and_then(|row| row.get(j))
                .copied()
                .unwrap_or(T::zero())
                * z[j];
        }
        for j in 0..m {
            r += s.g_ineq[j].get(i).copied().unwrap_or(T::zero()) * lambda[j];
        }
        for l in 0..p {
            r += s.a_eq[l].get(i).copied().unwrap_or(T::zero()) * nu[l];
        }
        stationarity = stationarity.max_val(r.abs());
    }
    let feas = (0..p).fold(T::zero(), |acc, l| {
        let az = s.a_eq[l]
            .iter()
            .zip(z.iter())
            .map(|(&a, &b)| a * b)
            .fold(T::zero(), |a, b| a + b);
        acc.max_val((az - s.b_eq.get(l).copied().unwrap_or(T::zero())).abs())
    });
    AgrawalBoydDifferentiableQpReport {
        optimal_primal_z: z.clone(),
        optimal_ineq_dual_lambda: lambda.clone(),
        optimal_eq_dual_nu: nu.clone(),
        optimal_objective_value: obj,
        kkt_stationarity_max_error: stationarity,
        primal_eq_feasibility_error: feas,
        reverse_kkt_grad_q_linear: z,
        reverse_kkt_grad_b_eq: nu.iter().map(|&x| -x).collect(),
        reverse_kkt_grad_h_ineq: lambda.iter().map(|&x| -x).collect(),
    }
}

#[derive(Clone, Debug)]
pub struct DiffPolynomialRegressor {
    pub mean: Vec<f64>,
    pub std: Vec<f64>,
    pub exponents: Vec<Vec<usize>>,
    pub coefficients: Vec<f64>,
}
fn enumerate(d: usize, degree: usize) -> Vec<Vec<usize>> {
    fn rec(i: usize, d: usize, left: usize, a: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if i == d {
            out.push(a.clone());
            return;
        }
        for k in 0..=left {
            a.push(k);
            rec(i + 1, d, left - k, a, out);
            a.pop();
        }
    }
    let mut out = Vec::new();
    rec(0, d, degree, &mut Vec::new(), &mut out);
    out
}
impl DiffPolynomialRegressor {
    fn basis(&self, x: &[f64]) -> (Vec<f64>, Vec<Vec<f64>>) {
        let d = self.mean.len();
        assert_eq!(x.len(), d, "regressor state dimension mismatch");
        let z: Vec<f64> = x
            .iter()
            .enumerate()
            .map(|(i, &v)| (v - self.mean[i]) / self.std[i])
            .collect();
        let mut p = Vec::with_capacity(self.exponents.len());
        let mut der = vec![vec![0.; self.exponents.len()]; d];
        for e in &self.exponents {
            let mut v = 1.;
            for i in 0..d {
                v *= z[i].powi(e[i] as i32);
            }
            let idx = p.len();
            p.push(v);
            for i in 0..d {
                if e[i] > 0 {
                    let mut q = 1.;
                    for j in 0..d {
                        if j != i {
                            q *= z[j].powi(e[j] as i32);
                        }
                    }
                    der[i][idx] = e[i] as f64 * z[i].powi((e[i] - 1) as i32) * q;
                }
            }
        }
        (p, der)
    }
    pub fn fit(
        states: &[Vec<f64>],
        payoffs: &[f64],
        diffs: &[Vec<f64>],
        degree: usize,
        diff_weight: f64,
        ridge: f64,
    ) -> Self {
        let m = states.len();
        assert!(
            m > 0
                && payoffs.len() == m
                && diffs.len() == m
                && states.iter().all(|x| x.len() == states[0].len())
                && diffs.iter().all(|x| x.len() == states[0].len())
                && states.iter().flatten().all(|x| x.is_finite())
                && payoffs.iter().all(|x| x.is_finite())
                && diffs.iter().flatten().all(|x| x.is_finite())
                && diff_weight.is_finite()
                && diff_weight >= 0.0
                && ridge.is_finite()
                && ridge >= 0.0,
            "invalid differential-regression dimensions or parameters"
        );
        let d = states[0].len();
        assert!(d > 0, "regression requires at least one state variable");
        let mean = (0..d)
            .map(|j| states.iter().map(|x| x[j]).sum::<f64>() / m as f64)
            .collect::<Vec<_>>();
        let std = (0..d)
            .map(|j| {
                (states.iter().map(|x| (x[j] - mean[j]).powi(2)).sum::<f64>() / m as f64)
                    .sqrt()
                    .max(1e-12)
            })
            .collect::<Vec<_>>();
        let exponents = enumerate(d, degree);
        let k = exponents.len();
        let tmp = Self {
            mean: mean.clone(),
            std: std.clone(),
            exponents: exponents.clone(),
            coefficients: vec![0.; k],
        };
        let mut lhs = DMatrix::zeros(k, k);
        let mut rhs = DVector::zeros(k);
        let y2 = payoffs.iter().map(|x| x * x).sum::<f64>() / m as f64;
        let lam: Vec<f64> = (0..d)
            .map(|q| {
                let z2 = diffs
                    .iter()
                    .map(|row| row[q] * std[q])
                    .map(|v| v * v)
                    .sum::<f64>()
                    / m as f64;
                diff_weight * y2 / z2.max(1e-12)
            })
            .collect();
        for s in 0..m {
            let (phi, der) = tmp.basis(&states[s]);
            for i in 0..k {
                rhs[i] += phi[i] * payoffs[s];
                for j in 0..k {
                    lhs[(i, j)] += phi[i] * phi[j];
                }
            }
            if diff_weight > 0. {
                for q in 0..d {
                    let z = diffs[s][q] * std[q];
                    let lam = lam[q];
                    for i in 0..k {
                        rhs[i] += lam * der[q][i] * z;
                        for j in 0..k {
                            lhs[(i, j)] += lam * der[q][i] * der[q][j];
                        }
                    }
                }
            }
        }
        for i in 0..k {
            lhs[(i, i)] += ridge * m as f64;
        }
        let sol = lhs
            .svd(true, true)
            .solve(&rhs, 1e-12)
            .expect("regression solve");
        Self {
            mean,
            std,
            exponents,
            coefficients: sol.iter().copied().collect(),
        }
    }
    pub fn predict_value(&self, x: &[f64]) -> f64 {
        let (p, _) = self.basis(x);
        p.iter().zip(&self.coefficients).map(|(a, b)| a * b).sum()
    }
    pub fn predict_gradient(&self, x: &[f64]) -> Vec<f64> {
        let (_, d) = self.basis(x);
        d.iter()
            .enumerate()
            .map(|(i, v)| {
                v.iter()
                    .zip(&self.coefficients)
                    .map(|(a, b)| a * b)
                    .sum::<f64>()
                    / self.std[i]
            })
            .collect()
    }
}

#[cfg(test)]
mod domain_suite_39_tests {
    use super::*;
    use ql_ad::Dual;
    fn make<T: Real>(x: T) -> AgrawalBoydDifferentiableQpSpec<T> {
        AgrawalBoydDifferentiableQpSpec {
            q_matrix: vec![vec![T::one()]],
            q_linear: vec![x],
            a_eq: vec![],
            b_eq: vec![],
            g_ineq: vec![],
            h_ineq: vec![],
            barrier_tau: T::from_f64(1e-4),
        }
    }
    #[test]
    fn dual_fd_parity() {
        let h = 1e-5;
        let fd = (solve_agrawal_boyd_differentiable_qp_layer(&make(1.0 + h))
            .optimal_objective_value
            - solve_agrawal_boyd_differentiable_qp_layer(&make(1.0 - h)).optimal_objective_value)
            / (2.0 * h);
        let d = solve_agrawal_boyd_differentiable_qp_layer(&make(Dual::<1>::variable(1.0, 0)))
            .optimal_objective_value;
        assert!(d.d[0].abs() > 1e-6 && (d.d[0] - fd).abs() < 1e-5);
    }

    #[test]
    fn suite39_invariants() {
        let r = solve_agrawal_boyd_differentiable_qp_layer(&AgrawalBoydDifferentiableQpSpec {
            q_matrix: vec![vec![2.0, 0.5], vec![0.5, 2.0]],
            q_linear: vec![0.0, 0.0],
            a_eq: vec![vec![1.0, 1.0]],
            b_eq: vec![1.0],
            g_ineq: vec![vec![-1.0, 0.0], vec![0.0, -1.0]],
            h_ineq: vec![0.0, 0.0],
            barrier_tau: 1e-4,
        });
        assert!(r.kkt_stationarity_max_error < 1e-9);
        assert!(r.primal_eq_feasibility_error < 1e-10);
        assert!((r.reverse_kkt_grad_b_eq[0] + r.optimal_eq_dual_nu[0]).abs() < 1e-6);
    }
}

mod evaluate_greene_wooldridge_heckman_survivorship_and_gmm_module {
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
            let pen =
                s.l2_penalty_lambda2 + s.l1_penalty_lambda1 / smooth.max_val(T::from_f64(1e-12));
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
    mod domain_suite_42_tests {
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
        fn suite42_heckman_exact_synthetic_recovery() {
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
        fn suite42_dml_dual_fd_and_invariants() {
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
            let ad_en =
                solve_hastie_tibshirani_elastic_net_and_gcv(&ad_spec).elastic_net_weights[0];
            assert!(ad_en.d[0].abs() > 1e-6 && (ad_en.d[0] - fd_en).abs() < 1e-5);
        }
    }
}
pub use evaluate_greene_wooldridge_heckman_survivorship_and_gmm_module::*;

mod evaluate_hastie_kernel_ridge_and_exact_loocv_module {
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
    mod domain_suite_47_tests {
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
        fn suite47_dml_invariants_and_dual_fd() {
            let r = evaluate_hastie_kernel_ridge_and_exact_loocv(&spec(0.0));
            assert!(r.loocv_shortcut_identity_error < 1e-10);
            assert!(r.effective_degrees_of_freedom > 0.0 && r.effective_degrees_of_freedom < 4.0);
            assert!(r.loocv_shortcut_mse >= r.in_sample_mse);
            let h = 1e-6;
            let fd = (evaluate_hastie_kernel_ridge_and_exact_loocv(&spec(h)).loocv_shortcut_mse
                - evaluate_hastie_kernel_ridge_and_exact_loocv(&spec(-h)).loocv_shortcut_mse)
                / (2.0 * h);
            let ad =
                evaluate_hastie_kernel_ridge_and_exact_loocv(&spec(Dual::<1>::variable(0.0, 0)))
                    .loocv_shortcut_mse;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_hastie_kernel_ridge_and_exact_loocv_module::*;

mod evaluate_hastie_zou_elastic_net_and_grouping_module {
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
    mod domain_suite_49_tests {
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
            let r_ad =
                evaluate_hastie_zou_elastic_net_and_grouping(&HastieZouElasticNetGroupingSpec {
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
            let fd = (evaluate_hastie_zou_elastic_net_and_grouping(&sp)
                .elastic_net_objective_value
                - evaluate_hastie_zou_elastic_net_and_grouping(&sm).elastic_net_objective_value)
                / (2.0 * h);
            assert!((r_ad.elastic_net_objective_value.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_hastie_zou_elastic_net_and_grouping_module::*;

mod evaluate_hastie_james_stein_empirical_bayes_sure_module {
    // Hastie, Tibshirani & Friedman (2009) Positive-Part James–Stein Empirical-Bayes Shrinkage
    // & Stein's Unbiased Risk Estimate (`SURE`) Dominance Engine.

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
        let sure_risk =
            (mle_risk - (p_minus_2 * p_minus_2 * sig2 * sig2) / denom).max_val(T::zero());
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
    mod domain_suite_51_tests {
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
            let r_ad =
                evaluate_hastie_james_stein_empirical_bayes_sure(&HastieJamesSteinSureSpec {
                    sample_means_y: [
                        Dual::<1>::variable(spec.sample_means_y[0], 0),
                        Dual::<1>::from_value(spec.sample_means_y[1]),
                        Dual::<1>::from_value(spec.sample_means_y[2]),
                        Dual::<1>::from_value(spec.sample_means_y[3]),
                    ],
                    grand_prior_target_mu0: Dual::<1>::from_value(spec.grand_prior_target_mu0),
                    noise_variance_sigma2: Dual::<1>::from_value(spec.noise_variance_sigma2),
                    singular_values_squared: spec
                        .singular_values_squared
                        .map(Dual::<1>::from_value),
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
}
pub use evaluate_hastie_james_stein_empirical_bayes_sure_module::*;

mod evaluate_htf_pls_vs_pcr_krylov_shrinkage_module {
    // Hastie, Tibshirani & Friedman (2009) Partial Least Squares (`PLS` Krylov Subspace)
    // vs. Principal Component Regression (`PCR`) & Ridge Shrinkage Factors Engine.

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
    mod domain_suite_53_tests {
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
}
pub use evaluate_htf_pls_vs_pcr_krylov_shrinkage_module::*;
