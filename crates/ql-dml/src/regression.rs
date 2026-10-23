use nalgebra::{DMatrix, DVector};
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
        assert!(m > 0 && payoffs.len() == m && diffs.len() == m);
        let d = states[0].len();
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
