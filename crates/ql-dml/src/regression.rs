#![allow(clippy::needless_range_loop)]

pub use crate::phase30::*;
pub use crate::phase33::*;
pub use crate::phase34::*;
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
mod phase27_tests {
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
    fn phase27_invariants() {
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
