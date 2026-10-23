use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct BlackLittermanSpec<T: Real> {
    pub covariance: Vec<Vec<T>>,
    pub market_weights: Vec<T>,
    pub risk_aversion: T,
    pub tau: T,
    pub view_matrix: Vec<Vec<T>>,
    pub view_returns: Vec<T>,
    pub omega: Vec<Vec<T>>,
}
#[derive(Clone, Debug)]
pub struct BlackLittermanPosterior<T: Real> {
    pub expected_returns: Vec<T>,
    pub covariance: Vec<Vec<T>>,
    pub weights: Vec<T>,
}
impl<T: Real> BlackLittermanSpec<T> {
    pub fn new(
        covariance: Vec<Vec<T>>,
        market_weights: Vec<T>,
        risk_aversion: T,
        tau: T,
        view_matrix: Vec<Vec<T>>,
        view_returns: Vec<T>,
        omega: Vec<Vec<T>>,
    ) -> Self {
        Self {
            covariance,
            market_weights,
            risk_aversion,
            tau,
            view_matrix,
            view_returns,
            omega,
        }
    }
}

#[allow(clippy::needless_range_loop)]
fn solve_f64(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut m = vec![vec![0.; n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = a[i][j];
        }
        m[i][n] = b[i];
    }
    for k in 0..n {
        let p = (k..n)
            .max_by(|&i, &j| m[i][k].abs().partial_cmp(&m[j][k].abs()).unwrap())
            .unwrap();
        m.swap(k, p);
        let d = m[k][k].max(1e-14);
        for j in k..=n {
            m[k][j] /= d;
        }
        for i in 0..n {
            if i != k {
                let x = m[i][k];
                for j in k..=n {
                    m[i][j] -= x * m[k][j];
                }
            }
        }
    }
    m.into_iter().map(|r| r[n]).collect()
}
pub fn implied_equilibrium_returns<T: Real>(
    cov: &[Vec<T>],
    weights: &[T],
    risk_aversion: T,
) -> Vec<T> {
    (0..cov.len())
        .map(|i| {
            risk_aversion
                * cov[i]
                    .iter()
                    .zip(weights)
                    .map(|(a, b)| *a * *b)
                    .fold(T::zero(), |a, b| a + b)
        })
        .collect()
}
pub fn black_litterman_posterior<T: Real>(s: &BlackLittermanSpec<T>) -> BlackLittermanPosterior<T> {
    let n = s.covariance.len();
    let k = s.view_returns.len();
    let pi = implied_equilibrium_returns(&s.covariance, &s.market_weights, s.risk_aversion);
    if k == 0 {
        return BlackLittermanPosterior {
            expected_returns: pi.clone(),
            covariance: s.covariance.clone(),
            weights: s.market_weights.clone(),
        };
    }
    let cov = nalgebra::DMatrix::from_fn(n, n, |i, j| s.covariance[i][j].val());
    let p = nalgebra::DMatrix::from_fn(k, n, |i, j| s.view_matrix[i][j].val());
    let om = nalgebra::DMatrix::from_fn(k, k, |i, j| s.omega[i][j].val());
    let tau = s.tau.val();
    let middle = (p.clone() * (&cov * tau) * p.transpose() + om)
        .try_inverse()
        .unwrap();
    let piv = nalgebra::DVector::from_iterator(
        k,
        s.view_returns.iter().enumerate().map(|(i, q)| {
            q.val()
                - p.row(i).dot(&nalgebra::DVector::from_iterator(
                    n,
                    pi.iter().map(|x| x.val()),
                ))
        }),
    );
    let mu = nalgebra::DVector::from_iterator(n, pi.iter().map(|x| x.val()))
        + &cov * p.transpose() * middle * piv;
    let weighted = cov.clone() * s.risk_aversion.val();
    let weighted_vec = (0..n)
        .map(|i| (0..n).map(|j| weighted[(i, j)]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let w = solve_f64(&weighted_vec, mu.as_slice());
    BlackLittermanPosterior {
        expected_returns: mu.iter().map(|&x| T::from_f64(x)).collect(),
        covariance: s.covariance.clone(),
        weights: w.iter().map(|&x| T::from_f64(x)).collect(),
    }
}
pub fn markowitz_minimum_variance_weights<T: Real>(cov: &[Vec<T>]) -> Vec<T> {
    let n = cov.len();
    let ones = vec![1.; n];
    let x = solve_f64(
        &cov.iter()
            .map(|r| r.iter().map(|x| x.val()).collect())
            .collect::<Vec<Vec<_>>>(),
        &ones,
    );
    let z = x.iter().sum::<f64>().max(1e-14);
    x.into_iter().map(|x| T::from_f64(x / z)).collect()
}
pub fn markowitz_tangency_weights<T: Real>(cov: &[Vec<T>], excess: &[T]) -> Vec<T> {
    let a = cov
        .iter()
        .map(|r| r.iter().map(|x| x.val()).collect())
        .collect::<Vec<Vec<_>>>();
    let b = excess.iter().map(|x| x.val()).collect::<Vec<_>>();
    let x = solve_f64(&a, &b);
    let z = x.iter().sum::<f64>().max(1e-14);
    x.into_iter().map(|x| T::from_f64(x / z)).collect()
}

fn ccd(c: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut y = b.iter().map(|x| x.sqrt()).collect::<Vec<_>>();
    for _ in 0..5000 {
        let mut max: f64 = 0.0;
        for i in 0..n {
            let off = (0..n)
                .filter(|&j| j != i)
                .map(|j| c[i][j] * y[j])
                .sum::<f64>();
            let disc = off * off + 4. * c[i][i].max(1e-14) * b[i];
            let ny = (-off + disc.sqrt()) / (2. * c[i][i].max(1e-14));
            max = max.max((ny - y[i]).abs());
            y[i] = ny;
        }
        if max < 1e-14 {
            break;
        }
    }
    for _ in 0..5 {
        let f = (0..n)
            .map(|i| (0..n).map(|j| c[i][j] * y[j]).sum::<f64>() - b[i] / y[i])
            .collect::<Vec<_>>();
        let h = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            c[i][j] + b[i] / (y[i] * y[i])
                        } else {
                            c[i][j]
                        }
                    })
                    .collect()
            })
            .collect::<Vec<Vec<_>>>();
        let d = solve_f64(&h, &f);
        for i in 0..n {
            y[i] -= d[i];
        }
    }
    y
}
pub fn risk_budgeting_weights<T: Real>(cov: &[Vec<T>], budgets: &[T]) -> Vec<T> {
    let n = cov.len();
    let sum = budgets.iter().map(|x| x.val()).sum::<f64>().max(1e-14);
    let b = budgets.iter().map(|x| x.val() / sum).collect::<Vec<_>>();
    let cf = cov
        .iter()
        .map(|r| r.iter().map(|x| x.val()).collect())
        .collect::<Vec<Vec<_>>>();
    let y = ccd(&cf, &b);
    let h = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        cf[i][j] + b[i] / (y[i] * y[i])
                    } else {
                        cf[i][j]
                    }
                })
                .collect()
        })
        .collect::<Vec<Vec<_>>>();
    let inv = (0..n)
        .map(|i| {
            solve_f64(
                &h,
                &(0..n)
                    .map(|j| if i == j { 1. } else { 0. })
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<Vec<_>>();
    let budget_sum_t = budgets.iter().fold(T::zero(), |a, x| a + *x);
    let mut yt = vec![T::zero(); n];
    for i in 0..n {
        let d = (0..n)
            .map(|j| {
                T::from_f64(-inv[j][i])
                    * ((0..n)
                        .map(|k| cov[j][k] * y[k])
                        .fold(T::zero(), |a, x| a + x)
                        - budgets[j] / budget_sum_t / y[j])
            })
            .fold(T::zero(), |a, x| a + x);
        yt[i] = T::from_f64(y[i]) + d;
    }
    let z = yt.iter().fold(T::zero(), |a, x| a + *x);
    yt.into_iter().map(|x| x / z).collect()
}
pub fn equal_risk_contribution_weights<T: Real>(cov: &[Vec<T>]) -> Vec<T> {
    let n = cov.len();
    risk_budgeting_weights(cov, &vec![T::one() / n as f64; n])
}
pub fn euler_risk_contributions<T: Real>(cov: &[Vec<T>], w: &[T]) -> (T, Vec<T>) {
    let mut q = T::zero();
    for i in 0..w.len() {
        for j in 0..w.len() {
            q += w[i] * cov[i][j] * w[j];
        }
    }
    let sig = q.sqrt();
    let rc = (0..w.len())
        .map(|i| {
            w[i] * (0..w.len())
                .map(|j| cov[i][j] * w[j])
                .fold(T::zero(), |a, x| a + x)
                / sig
        })
        .collect();
    (sig, rc)
}
