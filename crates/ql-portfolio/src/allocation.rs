use ql_ad::Real;
use thiserror::Error;

#[derive(Debug, Error, Clone, Copy, PartialEq, Eq)]
pub enum AllocationError {
    #[error("singular covariance or view matrix")]
    SingularCovariance,
    #[error("invalid matrix dimensions")]
    InvalidDimensions,
    #[error("infeasible box bounds")]
    InfeasibleBounds,
}

#[derive(Clone, Debug)]
pub struct BoxConstrainedMarkowitzSpec<T: Real> {
    pub covariance: Vec<Vec<T>>,
    pub expected_returns: Vec<T>,
    pub risk_aversion: T,
    pub lower_bounds: Vec<T>,
    pub upper_bounds: Vec<T>,
}

pub fn try_box_constrained_markowitz_weights<T: Real>(
    s: &BoxConstrainedMarkowitzSpec<T>,
) -> Result<Vec<T>, AllocationError> {
    let n = s.expected_returns.len();
    if n == 0
        || s.covariance.len() != n
        || s.covariance.iter().any(|r| r.len() != n)
        || s.lower_bounds.len() != n
        || s.upper_bounds.len() != n
        || !s.risk_aversion.val().is_finite()
        || s.covariance
            .iter()
            .flat_map(|r| r.iter())
            .any(|x| !x.val().is_finite())
        || s.expected_returns.iter().any(|x| !x.val().is_finite())
        || s.lower_bounds
            .iter()
            .chain(s.upper_bounds.iter())
            .any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    if s.risk_aversion.val() <= 0.0 {
        return Err(AllocationError::InvalidDimensions);
    }
    if s.lower_bounds
        .iter()
        .zip(&s.upper_bounds)
        .any(|(l, u)| l.val() > u.val() + 1e-12)
        || s.lower_bounds.iter().map(|x| x.val()).sum::<f64>() > 1.0 + 1e-12
        || s.upper_bounds.iter().map(|x| x.val()).sum::<f64>() < 1.0 - 1e-12
    {
        return Err(AllocationError::InfeasibleBounds);
    }
    // Determine the active set on primal f64 data.  No iterative solve is
    // recorded on an AD tape; the reduced KKT system is lifted once below.
    let mut state = vec![0i8; n]; // -1 lower, +1 upper, 0 free
    let cf = s
        .covariance
        .iter()
        .map(|r| r.iter().map(|x| x.val()).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let muf = s
        .expected_returns
        .iter()
        .map(|x| x.val())
        .collect::<Vec<_>>();
    let lf = s.lower_bounds.iter().map(|x| x.val()).collect::<Vec<_>>();
    let uf = s.upper_bounds.iter().map(|x| x.val()).collect::<Vec<_>>();
    let gam = s.risk_aversion.val();
    let solve_primal = |state: &[i8]| -> Result<(Vec<f64>, f64), AllocationError> {
        let free = (0..n).filter(|&i| state[i] == 0).collect::<Vec<_>>();
        if free.is_empty() {
            let w = state
                .iter()
                .enumerate()
                .map(|(i, z)| if *z < 0 { lf[i] } else { uf[i] })
                .collect::<Vec<_>>();
            let g_fn = |i: usize| muf[i] - gam * (0..n).map(|j| cf[i][j] * w[j]).sum::<f64>();
            let lo_lam = (0..n)
                .filter(|&i| state[i] > 0)
                .map(g_fn)
                .fold(f64::NEG_INFINITY, f64::max);
            let hi_lam = (0..n)
                .filter(|&i| state[i] < 0)
                .map(g_fn)
                .fold(f64::INFINITY, f64::min);
            let lambda = match (lo_lam.is_finite(), hi_lam.is_finite()) {
                (true, true) => 0.5 * (lo_lam + hi_lam),
                (true, false) => lo_lam,
                (false, true) => hi_lam,
                (false, false) => 0.0,
            };
            return Ok((w, lambda));
        }
        let m = free.len() + 1;
        let mut a = vec![vec![0.0; m]; m];
        let mut b = vec![0.0; m];
        let fixed_sum = (0..n)
            .filter(|&i| state[i] != 0)
            .map(|i| if state[i] < 0 { lf[i] } else { uf[i] })
            .sum::<f64>();
        for (ii, &i) in free.iter().enumerate() {
            for (jj, &j) in free.iter().enumerate() {
                a[ii][jj] = gam * cf[i][j];
            }
            a[ii][m - 1] = 1.0;
            a[m - 1][ii] = 1.0;
            b[ii] = muf[i]
                - (0..n)
                    .filter(|&j| state[j] != 0)
                    .map(|j| gam * cf[i][j] * if state[j] < 0 { lf[j] } else { uf[j] })
                    .sum::<f64>();
        }
        b[m - 1] = 1.0 - fixed_sum;
        let sol = solve_linear_generic::<f64>(&a, &b)?;
        let mut w = vec![0.0; n];
        for i in 0..n {
            w[i] = if state[i] < 0 {
                lf[i]
            } else if state[i] > 0 {
                uf[i]
            } else {
                0.0
            };
        }
        for (ii, &i) in free.iter().enumerate() {
            w[i] = sol[ii];
        }
        Ok((w, sol[m - 1]))
    };
    for _ in 0..(4 * n + 16) {
        let (w, lambda) = solve_primal(&state)?;
        let fixed_sum = (0..n).filter(|&i| state[i] != 0).map(|i| w[i]).sum::<f64>();
        if state.iter().all(|z| *z != 0) && (fixed_sum - 1.0).abs() > 1e-10 {
            // An all-active partition can only be feasible when the bounds
            // themselves sum to one.  Release one bound and re-solve the
            // equality KKT system instead of returning an infeasible budget.
            let want_upper = fixed_sum > 1.0;
            if let Some(i) =
                (0..n).find(|&i| (want_upper && state[i] > 0) || (!want_upper && state[i] < 0))
            {
                state[i] = 0;
                continue;
            }
        }
        let mut worst: Option<(usize, f64)> = None;
        for i in 0..n {
            if state[i] == 0 {
                let v = if w[i] < lf[i] - 1e-11 {
                    lf[i] - w[i]
                } else if w[i] > uf[i] + 1e-11 {
                    w[i] - uf[i]
                } else {
                    0.0
                };
                if v > 0.0 && worst.as_ref().is_none_or(|(_, x)| v > *x) {
                    worst = Some((i, v));
                }
            }
        }
        if let Some((i, _)) = worst {
            state[i] = if w[i] < lf[i] { -1 } else { 1 };
            continue;
        }
        // Release an incorrectly active bound using the KKT multiplier.
        let mut release: Option<(usize, f64)> = None;
        for i in 0..n {
            if state[i] == 0 {
                continue;
            }
            let grad = gam * (0..n).map(|j| cf[i][j] * w[j]).sum::<f64>() - muf[i] + lambda;
            let violation = if state[i] < 0 {
                (-grad).max(0.0)
            } else {
                grad.max(0.0)
            };
            if violation > 1e-10 && release.is_none_or(|(_, x)| violation > x) {
                release = Some((i, violation));
            }
        }
        if let Some((i, _)) = release {
            state[i] = 0;
            continue;
        }
        break;
    }
    let _ = solve_primal(&state)?;
    let free = (0..n).filter(|&i| state[i] == 0).collect::<Vec<_>>();
    let mut w = vec![T::zero(); n];
    let fixed = (0..n).filter(|&i| state[i] != 0).fold(T::zero(), |a, i| {
        a + if state[i] < 0 {
            s.lower_bounds[i]
        } else {
            s.upper_bounds[i]
        }
    });
    for i in 0..n {
        if state[i] < 0 {
            w[i] = s.lower_bounds[i];
        } else if state[i] > 0 {
            w[i] = s.upper_bounds[i];
        }
    }
    if free.is_empty() {
        return Ok(w);
    }
    if free.len() == 1 {
        w[free[0]] = T::one() - fixed;
        return Ok(w);
    }
    let m = free.len() + 1;
    let mut a = vec![vec![T::zero(); m]; m];
    let mut b = vec![T::zero(); m];
    for (ii, &i) in free.iter().enumerate() {
        for (jj, &j) in free.iter().enumerate() {
            a[ii][jj] = s.risk_aversion * s.covariance[i][j];
        }
        a[ii][m - 1] = T::one();
        a[m - 1][ii] = T::one();
        b[ii] = s.expected_returns[i]
            - (0..n)
                .filter(|&j| state[j] != 0)
                .map(|j| s.risk_aversion * s.covariance[i][j] * w[j])
                .fold(T::zero(), |x, y| x + y);
    }
    b[m - 1] = T::one() - fixed;
    let sol = solve_linear_generic(&a, &b)?;
    for (ii, &i) in free.iter().enumerate() {
        w[i] = sol[ii];
    }
    Ok(w)
}
pub fn box_constrained_markowitz_weights<T: Real>(s: &BoxConstrainedMarkowitzSpec<T>) -> Vec<T> {
    try_box_constrained_markowitz_weights(s).expect("feasible Markowitz bounds")
}

pub fn try_long_only_markowitz_weights<T: Real>(
    covariance: &[Vec<T>],
    expected_returns: &[T],
    risk_aversion: T,
) -> Result<Vec<T>, AllocationError> {
    let n = expected_returns.len();
    try_box_constrained_markowitz_weights(&BoxConstrainedMarkowitzSpec {
        covariance: covariance.to_vec(),
        expected_returns: expected_returns.to_vec(),
        risk_aversion,
        lower_bounds: vec![T::zero(); n],
        upper_bounds: vec![T::one(); n],
    })
}

pub fn long_only_markowitz_weights<T: Real>(
    covariance: &[Vec<T>],
    expected_returns: &[T],
    risk_aversion: T,
) -> Vec<T> {
    try_long_only_markowitz_weights(covariance, expected_returns, risk_aversion)
        .expect("feasible long-only Markowitz inputs")
}

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

/// AD-preserving Gaussian elimination with partial pivoting.  Pivot choice is
/// primal-only, but every row operation is performed in `T`, so Dual/Dual2
/// and TapeVar derivatives survive the solve.
pub fn solve_linear_generic<T: Real>(a: &[Vec<T>], b: &[T]) -> Result<Vec<T>, AllocationError> {
    let n = b.len();
    if n == 0
        || a.len() != n
        || a.iter().any(|r| r.len() != n)
        || a.iter().flatten().any(|x| !x.val().is_finite())
        || b.iter().any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let mut m = vec![vec![T::zero(); n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = a[i][j];
        }
        m[i][n] = b[i];
    }
    for k in 0..n {
        let Some(p) = (k..n).max_by(|&i, &j| m[i][k].val().abs().total_cmp(&m[j][k].val().abs()))
        else {
            return Err(AllocationError::SingularCovariance);
        };
        if !m[p][k].val().is_finite() || m[p][k].val().abs() <= 1e-14 {
            return Err(AllocationError::SingularCovariance);
        }
        m.swap(k, p);
        let pivot = m[k][k];
        for value in m[k].iter_mut().take(n + 1).skip(k) {
            *value /= pivot;
        }
        for i in 0..n {
            if i == k {
                continue;
            }
            let f = m[i][k];
            let pivot_row = m[k].clone();
            for j in k..=n {
                m[i][j] -= f * pivot_row[j];
            }
        }
    }
    Ok(m.into_iter().map(|r| r[n]).collect())
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
pub fn try_black_litterman_posterior<T: Real>(
    s: &BlackLittermanSpec<T>,
) -> Result<BlackLittermanPosterior<T>, AllocationError> {
    let n = s.covariance.len();
    let k = s.view_returns.len();
    if n == 0
        || s.market_weights.len() != n
        || s.covariance.iter().any(|r| r.len() != n)
        || s.view_matrix.len() != k
        || s.view_matrix.iter().any(|r| r.len() != n)
        || s.omega.len() != k
        || s.omega.iter().any(|r| r.len() != k)
        || !s.risk_aversion.val().is_finite()
        || s.risk_aversion.val() <= 0.0
        || !s.tau.val().is_finite()
        || s.tau.val() < 0.0
        || s.market_weights.iter().any(|x| !x.val().is_finite())
        || s.covariance
            .iter()
            .flat_map(|r| r.iter())
            .any(|x| !x.val().is_finite())
        || s.view_matrix
            .iter()
            .flat_map(|r| r.iter())
            .chain(s.view_returns.iter())
            .chain(s.omega.iter().flat_map(|r| r.iter()))
            .any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let pi = implied_equilibrium_returns(&s.covariance, &s.market_weights, s.risk_aversion);
    if k == 0 {
        return Ok(BlackLittermanPosterior {
            expected_returns: pi.clone(),
            covariance: s.covariance.clone(),
            weights: s.market_weights.clone(),
        });
    }
    let tau_cov = s
        .covariance
        .iter()
        .map(|row| row.iter().map(|x| *x * s.tau).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut middle = vec![vec![T::zero(); k]; k];
    for i in 0..k {
        for j in 0..k {
            middle[i][j] = s.omega[i][j]
                + (0..n)
                    .flat_map(|a| (0..n).map(move |b| (a, b)))
                    .map(|(a, b)| s.view_matrix[i][a] * tau_cov[a][b] * s.view_matrix[j][b])
                    .fold(T::zero(), |x, y| x + y);
        }
    }
    let mut rhs_y = vec![vec![T::zero(); n]; k];
    for i in 0..k {
        for j in 0..n {
            rhs_y[i][j] = (0..n)
                .map(|a| s.view_matrix[i][a] * tau_cov[a][j])
                .fold(T::zero(), |x, y| x + y);
        }
    }
    let mut ymat = vec![vec![T::zero(); n]; k];
    for j in 0..n {
        let rhs = (0..k).map(|i| rhs_y[i][j]).collect::<Vec<_>>();
        let sol = solve_linear_generic(&middle, &rhs)?;
        for i in 0..k {
            ymat[i][j] = sol[i];
        }
    }
    let rhs_v = (0..k)
        .map(|i| {
            s.view_returns[i]
                - (0..n)
                    .map(|j| s.view_matrix[i][j] * pi[j])
                    .fold(T::zero(), |x, y| x + y)
        })
        .collect::<Vec<_>>();
    let v = solve_linear_generic(&middle, &rhs_v)?;
    let mu = (0..n)
        .map(|i| {
            pi[i]
                + (0..k)
                    .map(|a| {
                        (0..n)
                            .map(|j| tau_cov[i][j] * s.view_matrix[a][j])
                            .fold(T::zero(), |x, y| x + y)
                            * v[a]
                    })
                    .fold(T::zero(), |x, y| x + y)
        })
        .collect::<Vec<_>>();
    let mut post_cov = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            post_cov[i][j] = s.covariance[i][j] + tau_cov[i][j]
                - (0..k)
                    .map(|a| {
                        (0..n)
                            .map(|b| tau_cov[i][b] * s.view_matrix[a][b])
                            .fold(T::zero(), |x, y| x + y)
                            * ymat[a][j]
                    })
                    .fold(T::zero(), |x, y| x + y);
        }
    }
    let weighted = s
        .covariance
        .iter()
        .map(|row| row.iter().map(|x| *x * s.risk_aversion).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let raw = solve_linear_generic(&weighted, &mu)?;
    let z = raw.iter().fold(T::zero(), |a, x| a + *x);
    if !z.val().is_finite() || z.val().abs() <= 1e-14 {
        return Err(AllocationError::SingularCovariance);
    }
    let w = raw.iter().map(|x| *x / z).collect::<Vec<_>>();
    Ok(BlackLittermanPosterior {
        expected_returns: mu,
        covariance: post_cov,
        weights: w,
    })
}
pub fn black_litterman_posterior<T: Real>(s: &BlackLittermanSpec<T>) -> BlackLittermanPosterior<T> {
    try_black_litterman_posterior(s).expect("valid Black-Litterman covariance/view system")
}
pub fn try_markowitz_minimum_variance_weights<T: Real>(
    cov: &[Vec<T>],
) -> Result<Vec<T>, AllocationError> {
    let n = cov.len();
    if n == 0 || cov.iter().any(|r| r.len() != n) {
        return Err(AllocationError::InvalidDimensions);
    }
    let ones = vec![T::one(); n];
    let x = solve_linear_generic(cov, &ones)?;
    let z = x.iter().fold(T::zero(), |a, b| a + *b);
    if !z.val().is_finite() || z.val().abs() <= 1e-14 {
        return Err(AllocationError::SingularCovariance);
    }
    Ok(x.into_iter().map(|x| x / z).collect())
}
pub fn markowitz_minimum_variance_weights<T: Real>(cov: &[Vec<T>]) -> Vec<T> {
    try_markowitz_minimum_variance_weights(cov).expect("valid covariance matrix")
}
pub fn try_markowitz_tangency_weights<T: Real>(
    cov: &[Vec<T>],
    excess: &[T],
) -> Result<Vec<T>, AllocationError> {
    let n = cov.len();
    if n == 0 || excess.len() != n || cov.iter().any(|r| r.len() != n) {
        return Err(AllocationError::InvalidDimensions);
    }
    let x = solve_linear_generic(cov, excess)?;
    let z = x.iter().fold(T::zero(), |a, b| a + *b);
    if !z.val().is_finite() || z.val().abs() <= 1e-14 {
        return Err(AllocationError::SingularCovariance);
    }
    Ok(x.into_iter().map(|x| x / z).collect())
}
pub fn markowitz_tangency_weights<T: Real>(cov: &[Vec<T>], excess: &[T]) -> Vec<T> {
    try_markowitz_tangency_weights(cov, excess).expect("valid covariance matrix")
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
        let d = solve_linear_generic::<f64>(&h, &f).expect("risk-budgeting Hessian is SPD");
        for i in 0..n {
            y[i] -= d[i];
        }
    }
    y
}
pub fn try_risk_budgeting_weights<T: Real>(
    cov: &[Vec<T>],
    budgets: &[T],
) -> Result<Vec<T>, AllocationError> {
    let n = cov.len();
    if n == 0
        || budgets.len() != n
        || cov.iter().any(|r| r.len() != n)
        || budgets
            .iter()
            .any(|x| !x.val().is_finite() || x.val() <= 0.0)
        || cov
            .iter()
            .flat_map(|r| r.iter())
            .any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    // Risk budgeting is defined for a symmetric positive-definite covariance
    // matrix.  Without this guard the coordinate descent below can generate
    // NaNs (negative diagonal or an asymmetric input) and those NaNs used to
    // leak out as apparently valid weights.
    for i in 0..n {
        if cov[i][i].val() <= 0.0 {
            return Err(AllocationError::SingularCovariance);
        }
        for j in (i + 1)..n {
            if (cov[i][j].val() - cov[j][i].val()).abs() > 1e-10 {
                return Err(AllocationError::InvalidDimensions);
            }
        }
    }
    let sum = budgets.iter().map(|x| x.val()).sum::<f64>().max(1e-14);
    let b = budgets.iter().map(|x| x.val() / sum).collect::<Vec<_>>();
    let cf = cov
        .iter()
        .map(|r| r.iter().map(|x| x.val()).collect())
        .collect::<Vec<Vec<_>>>();
    let y = ccd(&cf, &b);
    if y.iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err(AllocationError::SingularCovariance);
    }
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
            solve_linear_generic::<f64>(
                &h,
                &(0..n)
                    .map(|j| if i == j { 1. } else { 0. })
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let budget_sum_t = budgets.iter().fold(T::zero(), |a, x| a + *x);
    let mut yt = vec![T::zero(); n];
    for i in 0..n {
        let d = (0..n)
            .map(|j| {
                T::from_f64(-inv[j][i])
                    * (((0..n)
                        .map(|k| cov[j][k] * y[k])
                        .fold(T::zero(), |a, x| a + x)
                        - budgets[j] / budget_sum_t / y[j])
                        - T::from_f64((0..n).map(|k| cf[j][k] * y[k]).sum::<f64>() - b[j] / y[j]))
            })
            .fold(T::zero(), |a, x| a + x);
        yt[i] = T::from_f64(y[i]) + d;
    }
    let z = yt.iter().fold(T::zero(), |a, x| a + *x);
    if !z.val().is_finite() || z.val().abs() <= 1e-14 {
        return Err(AllocationError::SingularCovariance);
    }
    Ok(yt.into_iter().map(|x| x / z).collect())
}
pub fn risk_budgeting_weights<T: Real>(cov: &[Vec<T>], budgets: &[T]) -> Vec<T> {
    try_risk_budgeting_weights(cov, budgets).expect("valid risk-budgeting covariance and budgets")
}
pub fn equal_risk_contribution_weights<T: Real>(cov: &[Vec<T>]) -> Vec<T> {
    let n = cov.len();
    if n == 0 {
        return Vec::new();
    }
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

#[derive(Clone, Debug)]
pub struct MeanCvarSpec<T: Real> {
    pub scenario_returns: Vec<Vec<T>>,
    pub expected_returns: Vec<T>,
    pub risk_aversion: T,
    pub confidence: f64,
    pub smoothing_eps: f64,
    pub lower_bounds: Vec<T>,
    pub upper_bounds: Vec<T>,
}
#[derive(Clone, Debug)]
pub struct MeanCvarResult<T: Real> {
    pub weights: Vec<T>,
    pub var_threshold: T,
    pub cvar: T,
    pub expected_return: T,
}
fn project_capped_simplex<T: Real>(
    v: &[T],
    lower: &[T],
    upper: &[T],
) -> Result<Vec<T>, AllocationError> {
    let n = v.len();
    if n == 0
        || lower.len() != n
        || upper.len() != n
        || lower
            .iter()
            .zip(upper)
            .any(|(l, u)| l.val() > u.val() + 1e-12)
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let lo_sum = lower.iter().copied().fold(T::zero(), |a, z| a + z);
    let hi_sum = upper.iter().copied().fold(T::zero(), |a, z| a + z);
    if lo_sum.val() > 1.0 + 1e-12 || hi_sum.val() < 1.0 - 1e-12 {
        return Err(AllocationError::InfeasibleBounds);
    }
    let mut lo = v
        .iter()
        .zip(upper)
        .map(|(x, u)| x.val() - u.val())
        .fold(f64::INFINITY, f64::min);
    let mut hi = v
        .iter()
        .zip(lower)
        .map(|(x, l)| x.val() - l.val())
        .fold(f64::NEG_INFINITY, f64::max);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let sum: f64 = v
            .iter()
            .zip(lower)
            .zip(upper)
            .map(|((x, l), u)| (x.val() - mid).max(l.val()).min(u.val()))
            .sum();
        if sum > 1.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mid = 0.5 * (lo + hi);
    let mut free_sum = T::zero();
    let mut clamped_sum = T::zero();
    let mut free_cnt = 0usize;
    for i in 0..v.len() {
        let cand = v[i].val() - mid;
        if cand <= lower[i].val() {
            clamped_sum += lower[i];
        } else if cand >= upper[i].val() {
            clamped_sum += upper[i];
        } else {
            free_sum += v[i];
            free_cnt += 1;
        }
    }
    let shift = if free_cnt > 0 {
        (free_sum + clamped_sum - T::one()) / (free_cnt as f64)
    } else {
        T::from_f64(mid)
    };
    Ok((0..v.len())
        .map(|i| {
            let cand = v[i].val() - mid;
            if cand <= lower[i].val() {
                lower[i]
            } else if cand >= upper[i].val() {
                upper[i]
            } else {
                (v[i] - shift).max_val(lower[i]).min_val(upper[i])
            }
        })
        .collect())
}

pub fn try_mean_cvar_weights<T: Real>(
    s: &MeanCvarSpec<T>,
) -> Result<MeanCvarResult<T>, AllocationError> {
    let n = s.expected_returns.len();
    let m_scen = s.scenario_returns.len();
    if n == 0
        || m_scen == 0
        || s.scenario_returns.iter().any(|r| r.len() != n)
        || s.lower_bounds.len() != n
        || s.upper_bounds.len() != n
        || !(0.5..1.0).contains(&s.confidence)
        || !s.risk_aversion.val().is_finite()
        || s.risk_aversion.val() < 0.0
        || !s.smoothing_eps.is_finite()
        || s.smoothing_eps < 0.0
        || s.scenario_returns
            .iter()
            .flat_map(|r| r.iter())
            .chain(s.expected_returns.iter())
            .chain(s.lower_bounds.iter())
            .chain(s.upper_bounds.iter())
            .any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let eps = T::from_f64(s.smoothing_eps.max(1e-5));
    let tail_mass = T::from_f64(((1.0 - s.confidence) * m_scen as f64).max(1e-12));
    let idx = ((s.confidence * m_scen as f64).floor() as usize).min(m_scen - 1);
    let init_guess: Vec<T> = (0..n).map(|_| T::from_f64(1.0 / n as f64)).collect();
    let mut w = project_capped_simplex(&init_guess, &s.lower_bounds, &s.upper_bounds)?;
    let solve_zeta = |weights: &[T]| -> (T, Vec<T>) {
        let losses: Vec<T> = s
            .scenario_returns
            .iter()
            .map(|r| {
                -r.iter()
                    .zip(weights)
                    .map(|(x, y)| *x * *y)
                    .fold(T::zero(), |a, z| a + z)
            })
            .collect();
        let mut sorted = losses.clone();
        sorted.sort_by(|a, b| a.val().total_cmp(&b.val()));
        let mut zeta = sorted[idx];
        for _ in 0..4 {
            let mut f = -tail_mass;
            let mut df = T::zero();
            for &l in &losses {
                let u = l - zeta;
                let rad = (u * u + eps * eps).sqrt();
                f += (T::one() + u / rad) * 0.5;
                df += eps * eps / (rad * rad * rad) * 0.5;
            }
            zeta += f / df.max_val(T::from_f64(1e-12));
        }
        (zeta, losses)
    };
    let eval_obj = |weights: &[T]| -> T {
        let (z, ls) = solve_zeta(weights);
        let mut h_sum = T::zero();
        for &l in &ls {
            let u = l - z;
            h_sum += (u + (u * u + eps * eps).sqrt()) * 0.5;
        }
        let cv = z + h_sum / tail_mass;
        let mu = s
            .expected_returns
            .iter()
            .zip(weights)
            .map(|(a, b)| *a * *b)
            .fold(T::zero(), |a, x| a + x);
        -mu + s.risk_aversion * cv
    };
    for _ in 0..32 {
        let (zeta, losses) = solve_zeta(&w);
        let mut grad = vec![T::zero(); n];
        let mut g_norm2 = 0.0;
        for i in 0..n {
            let mut tail_g = T::zero();
            for m in 0..m_scen {
                let u = losses[m] - zeta;
                let dphi = (T::one() + u / (u * u + eps * eps).sqrt()) * 0.5;
                tail_g -= (dphi / tail_mass) * s.scenario_returns[m][i];
            }
            grad[i] = -s.expected_returns[i] + s.risk_aversion * tail_g;
            g_norm2 += grad[i].val() * grad[i].val();
        }
        let g_norm = g_norm2.sqrt().max(1e-12);
        let base_obj = eval_obj(&w).val();
        let mut step = 0.5 / g_norm;
        for _ in 0..6 {
            let trial_in: Vec<T> = (0..n).map(|i| w[i] - grad[i] * step).collect();
            let cand = project_capped_simplex(&trial_in, &s.lower_bounds, &s.upper_bounds)?;
            if eval_obj(&cand).val() <= base_obj + 1e-12 {
                w = cand;
                break;
            }
            step *= 0.5;
        }
    }
    let (zeta, losses) = solve_zeta(&w);
    let mut smooth_hinge_sum = T::zero();
    for &l in &losses {
        let u = l - zeta;
        smooth_hinge_sum += (u + (u * u + eps * eps).sqrt()) * 0.5;
    }
    let cvar = zeta + smooth_hinge_sum / tail_mass;
    let er = s
        .expected_returns
        .iter()
        .zip(&w)
        .map(|(a, b)| *a * *b)
        .fold(T::zero(), |a, z| a + z);
    Ok(MeanCvarResult {
        weights: w,
        var_threshold: zeta,
        cvar,
        expected_return: er,
    })
}

pub fn mean_cvar_weights<T: Real>(s: &MeanCvarSpec<T>) -> MeanCvarResult<T> {
    try_mean_cvar_weights(s).expect("valid mean-CVaR inputs")
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ShrinkageTarget {
    ScaledIdentity,
    ConstantVarianceIdentity,
    DiagonalVariance,
    ConstantCorrelation,
}
#[derive(Clone, Debug)]
pub struct LedoitWolfResult<T: Real> {
    pub shrunk_covariance: Vec<Vec<T>>,
    pub sample_covariance: Vec<Vec<T>>,
    pub target_matrix: Vec<Vec<T>>,
    pub shrinkage_intensity: T,
}
pub fn ledoit_wolf_shrinkage<T: Real>(
    r: &[Vec<T>],
    target: ShrinkageTarget,
) -> Result<LedoitWolfResult<T>, AllocationError> {
    let n = r.len();
    if n < 2
        || r.iter().any(|x| {
            x.is_empty() || x.len() != r[0].len() || x.iter().any(|v| !v.val().is_finite())
        })
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let p = r[0].len();
    let mean = (0..p)
        .map(|j| r.iter().map(|x| x[j]).fold(T::zero(), |a, z| a + z) / n as f64)
        .collect::<Vec<_>>();
    let s = (0..p)
        .map(|i| {
            (0..p)
                .map(|j| {
                    r.iter()
                        .map(|x| (x[i] - mean[i]) * (x[j] - mean[j]))
                        .fold(T::zero(), |a, z| a + z)
                        / n as f64
                })
                .collect()
        })
        .collect::<Vec<Vec<_>>>();
    let mut f = vec![vec![T::zero(); p]; p];
    match target {
        ShrinkageTarget::ScaledIdentity | ShrinkageTarget::ConstantVarianceIdentity => {
            let m = (0..p).map(|i| s[i][i]).fold(T::zero(), |a, z| a + z) / p as f64;
            for (i, row) in f.iter_mut().enumerate().take(p) {
                row[i] = m;
            }
        }
        ShrinkageTarget::DiagonalVariance => {
            for (i, row) in f.iter_mut().enumerate().take(p) {
                row[i] = s[i][i];
            }
        }
        ShrinkageTarget::ConstantCorrelation => {
            let mut c = T::zero();
            let mut cnt: f64 = 0.;
            for i in 0..p {
                for j in 0..i {
                    c += s[i][j] / (s[i][i] * s[j][j]).sqrt().max_val(T::from_f64(1e-12));
                    cnt += 1.;
                }
            }
            let c = c / cnt.max(1.);
            for i in 0..p {
                for j in 0..p {
                    f[i][j] = if i == j {
                        s[i][i]
                    } else {
                        c * (s[i][i] * s[j][j]).sqrt()
                    };
                }
            }
        }
    }
    // Ledoit--Wolf intensity: estimate the variance of the sample
    // covariance elementwise (pi-hat) and divide by the squared distance to
    // the target (gamma).  The former fixed `0.1/gamma` heuristic was not a
    // shrinkage estimator and could materially over/under-shrink depending
    // only on units.
    let mut pi = T::zero();
    for row in r {
        for i in 0..p {
            for j in 0..p {
                let outer = (row[i] - mean[i]) * (row[j] - mean[j]);
                let err = outer - s[i][j];
                pi += err * err;
            }
        }
    }
    pi = pi / n as f64;
    let d = (0..p)
        .map(|i| {
            (0..p)
                .map(|j| (s[i][j] - f[i][j]) * (s[i][j] - f[i][j]))
                .fold(T::zero(), |a, z| a + z)
        })
        .fold(T::zero(), |a, z| a + z);
    let delta = (pi / (d + T::from_f64(1e-12)))
        .max_val(T::zero())
        .min_val(T::one());
    let shr = (0..p)
        .map(|i| {
            (0..p)
                .map(|j| delta * f[i][j] + (T::one() - delta) * s[i][j])
                .collect()
        })
        .collect();
    Ok(LedoitWolfResult {
        shrunk_covariance: shr,
        sample_covariance: s,
        target_matrix: f,
        shrinkage_intensity: delta,
    })
}
#[derive(Clone, Debug)]
pub struct RobustMarkowitzSpec<T: Real> {
    pub covariance: Vec<Vec<T>>,
    pub return_uncertainty_cov: Vec<Vec<T>>,
    pub expected_returns: Vec<T>,
    pub risk_aversion: T,
    pub uncertainty_radius: T,
    pub lower_bounds: Vec<T>,
    pub upper_bounds: Vec<T>,
}
pub fn try_robust_markowitz_weights<T: Real>(
    s: &RobustMarkowitzSpec<T>,
) -> Result<Vec<T>, AllocationError> {
    let n = s.expected_returns.len();
    // Validate every matrix/vector before the fixed-point loop.  The old
    // implementation only checked `n`, so malformed uncertainty matrices
    // could panic while indexing (and negative radii produced a nonsensical
    // adversarial covariance adjustment).
    if n == 0
        || s.covariance.len() != n
        || s.covariance.iter().any(|r| r.len() != n)
        || s.return_uncertainty_cov.len() != n
        || s.return_uncertainty_cov.iter().any(|r| r.len() != n)
        || s.lower_bounds.len() != n
        || s.upper_bounds.len() != n
        || !s.risk_aversion.val().is_finite()
        || s.risk_aversion.val() <= 0.0
        || !s.uncertainty_radius.val().is_finite()
        || s.uncertainty_radius.val() < 0.0
        || s.expected_returns
            .iter()
            .chain(s.lower_bounds.iter())
            .chain(s.upper_bounds.iter())
            .any(|x| !x.val().is_finite())
        || s.covariance
            .iter()
            .chain(s.return_uncertainty_cov.iter())
            .flat_map(|r| r.iter())
            .any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let mut cov = s.covariance.clone();
    let mut w = try_box_constrained_markowitz_weights(&BoxConstrainedMarkowitzSpec {
        covariance: cov.clone(),
        expected_returns: s.expected_returns.clone(),
        risk_aversion: s.risk_aversion,
        lower_bounds: s.lower_bounds.clone(),
        upper_bounds: s.upper_bounds.clone(),
    })?;
    for _ in 0..8 {
        let mut ow = T::zero();
        for i in 0..n {
            for j in 0..n {
                ow += w[i] * s.return_uncertainty_cov[i][j] * w[j];
            }
        }
        let scale = s.uncertainty_radius / (ow + T::from_f64(1e-12)).sqrt();
        for i in 0..n {
            for j in 0..n {
                cov[i][j] = s.covariance[i][j]
                    + scale * s.return_uncertainty_cov[i][j]
                        / s.risk_aversion.max_val(T::from_f64(1e-12));
            }
        }
        let nw = try_box_constrained_markowitz_weights(&BoxConstrainedMarkowitzSpec {
            covariance: cov.clone(),
            expected_returns: s.expected_returns.clone(),
            risk_aversion: s.risk_aversion,
            lower_bounds: s.lower_bounds.clone(),
            upper_bounds: s.upper_bounds.clone(),
        })?;
        if nw
            .iter()
            .zip(&w)
            .all(|(a, b)| (a.val() - b.val()).abs() < 1e-12)
        {
            w = nw;
            break;
        }
        w = nw;
    }
    Ok(w)
}
pub fn robust_markowitz_weights<T: Real>(s: &RobustMarkowitzSpec<T>) -> Vec<T> {
    try_robust_markowitz_weights(s).expect("valid robust Markowitz inputs")
}
