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
pub fn solve_agrawal_boyd_differentiable_qp_layer<T: Real>(
    s: &AgrawalBoydDifferentiableQpSpec<T>,
) -> AgrawalBoydDifferentiableQpReport<T> {
    let n = s.q_linear.len();
    let mut z = vec![T::zero(); n];
    for i in 0..n {
        let q = s
            .q_matrix
            .get(i)
            .and_then(|r| r.get(i))
            .copied()
            .unwrap_or(T::one())
            .max_val(T::from_f64(1e-8));
        z[i] = -s.q_linear[i] / q;
    }
    let mut nu = vec![T::zero(); s.a_eq.len()];
    let mut lam = vec![T::zero(); s.g_ineq.len()];
    for _ in 0..16 {
        for (j, g) in s.g_ineq.iter().enumerate() {
            let slack = s.h_ineq.get(j).copied().unwrap_or(T::zero())
                - g.iter()
                    .zip(z.iter())
                    .map(|(&a, &b)| a * b)
                    .fold(T::zero(), |a, b| a + b);
            lam[j] = s.barrier_tau / slack.max_val(T::from_f64(1e-8));
        }
        if let Some(a) = s.a_eq.first() {
            let mut num = s.b_eq.first().copied().unwrap_or(T::zero());
            let mut den = T::zero();
            for i in 0..n {
                let q = s
                    .q_matrix
                    .get(i)
                    .and_then(|r| r.get(i))
                    .copied()
                    .unwrap_or(T::one())
                    .max_val(T::from_f64(1e-8));
                let gi = s
                    .g_ineq
                    .iter()
                    .enumerate()
                    .map(|(j, g)| g.get(i).copied().unwrap_or(T::zero()) * lam[j])
                    .fold(T::zero(), |a, b| a + b);
                let ai = a.get(i).copied().unwrap_or(T::zero());
                num += ai * (s.q_linear[i] + gi) / q;
                den += ai * ai / q;
            }
            nu[0] = -num / den.max_val(T::from_f64(1e-8));
        }
        for i in 0..n {
            let q = s
                .q_matrix
                .get(i)
                .and_then(|r| r.get(i))
                .copied()
                .unwrap_or(T::one())
                .max_val(T::from_f64(1e-8));
            let gi = s
                .g_ineq
                .iter()
                .enumerate()
                .map(|(j, g)| g.get(i).copied().unwrap_or(T::zero()) * lam[j])
                .fold(T::zero(), |a, b| a + b);
            let ai = s
                .a_eq
                .first()
                .and_then(|a| a.get(i))
                .copied()
                .unwrap_or(T::zero());
            z[i] = -(s.q_linear[i] + gi + ai * nu.first().copied().unwrap_or(T::zero())) / q;
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
        for (j, g) in s.g_ineq.iter().enumerate() {
            r += g.get(i).copied().unwrap_or(T::zero()) * lam[j];
        }
        for (j, a) in s.a_eq.iter().enumerate() {
            r += a.get(i).copied().unwrap_or(T::zero()) * nu[j];
        }
        stationarity = stationarity.max_val(r.abs());
    }
    let feas = s
        .a_eq
        .iter()
        .enumerate()
        .map(|(j, a)| {
            (a.iter()
                .zip(z.iter())
                .map(|(&x, &y)| x * y)
                .fold(T::zero(), |u, v| u + v)
                - s.b_eq.get(j).copied().unwrap_or(T::zero()))
            .abs()
        })
        .fold(T::zero(), |a, b| a.max_val(b));
    AgrawalBoydDifferentiableQpReport {
        optimal_primal_z: z.clone(),
        optimal_ineq_dual_lambda: lam.clone(),
        optimal_eq_dual_nu: nu.clone(),
        optimal_objective_value: obj,
        kkt_stationarity_max_error: stationarity,
        primal_eq_feasibility_error: feas,
        reverse_kkt_grad_q_linear: z,
        reverse_kkt_grad_b_eq: nu.iter().map(|&x| -x).collect(),
        reverse_kkt_grad_h_ineq: lam.iter().map(|&x| -x).collect(),
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
        assert!((d.d[0] - fd).abs() < 1e-5);
    }
}
