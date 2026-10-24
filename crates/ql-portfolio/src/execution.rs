use super::allocation::AllocationError;
use ql_ad::Real;
#[derive(Clone, Debug)]
pub struct AlmgrenChrissSpec<T: Real> {
    pub initial_inventory: T,
    pub horizon: T,
    pub num_steps: usize,
    pub volatility: T,
    pub temporary_impact: T,
    pub permanent_impact: T,
    pub risk_aversion: T,
}
#[derive(Clone, Debug)]
pub struct ObizhaevaWangSpec<T: Real> {
    pub initial_inventory: T,
    pub trade_times: Vec<f64>,
    pub volatility: T,
    pub transient_impact: T,
    pub resilience_speed: T,
    pub temporary_impact: T,
    pub risk_aversion: T,
}
#[derive(Clone, Debug)]
pub struct OptimalExecutionReport<T: Real> {
    pub times: Vec<f64>,
    pub inventory_schedule: Vec<T>,
    pub trade_schedule: Vec<T>,
    pub expected_cost: T,
    pub cost_variance: T,
    pub mean_variance_utility: T,
}
pub fn try_almgren_chriss_execution<T: Real>(
    s: &AlmgrenChrissSpec<T>,
) -> Result<OptimalExecutionReport<T>, AllocationError> {
    if s.num_steps == 0 || s.horizon.val() <= 0. || s.temporary_impact.val() <= 0. {
        return Err(AllocationError::InvalidDimensions);
    }
    let n = s.num_steps;
    let dt = s.horizon / T::from_f64(n as f64);
    let eta_tilde =
        (s.temporary_impact - s.permanent_impact * dt * 0.5).max_val(T::from_f64(1e-12));
    let q = s.risk_aversion * s.volatility * s.volatility * dt * dt / eta_tilde;
    let u = T::one() + q * 0.5;
    let kappa = (u + (u * u - T::one()).max_val(T::zero()).sqrt()).ln() / dt;
    let mut times = Vec::with_capacity(n + 1);
    let mut inv = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let t = dt * T::from_f64(i as f64);
        times.push(t.val());
        let x = if q.val().abs() < 1e-10 {
            s.initial_inventory * (1. - i as f64 / n as f64)
        } else {
            s.initial_inventory * (kappa * (s.horizon - t)).sinh() / (kappa * s.horizon).sinh()
        };
        inv.push(x)
    }
    let trades = (1..=n).map(|i| inv[i - 1] - inv[i]).collect::<Vec<_>>();
    let mut cost = s.permanent_impact * s.initial_inventory * s.initial_inventory * 0.5;
    let mut sum = T::zero();
    for x in &trades {
        sum += *x * *x;
    }
    cost -= s.permanent_impact * sum * 0.5;
    cost += s.temporary_impact / dt * sum;
    let mut var = T::zero();
    for x in inv.iter().take(n).skip(1) {
        var += *x * *x;
    }
    var *= s.volatility * s.volatility * dt;
    Ok(OptimalExecutionReport {
        times,
        inventory_schedule: inv,
        trade_schedule: trades,
        expected_cost: cost,
        cost_variance: var,
        mean_variance_utility: cost + s.risk_aversion * var,
    })
}
pub fn almgren_chriss_execution<T: Real>(s: &AlmgrenChrissSpec<T>) -> OptimalExecutionReport<T> {
    try_almgren_chriss_execution(s).expect("valid Almgren-Chriss inputs")
}
pub fn try_obizhaeva_wang_execution<T: Real>(
    s: &ObizhaevaWangSpec<T>,
) -> Result<OptimalExecutionReport<T>, AllocationError> {
    let n = s.trade_times.len();
    if n < 2
        || s.trade_times.windows(2).any(|w| w[1] <= w[0])
        || s.resilience_speed.val() <= 0.0
        || s.transient_impact.val() < 0.0
        || s.temporary_impact.val() < 0.0
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let m = n - 1;
    let t0 = s.trade_times[0];
    let mut a = vec![vec![T::zero(); m]; m];
    let mut b = vec![T::one(); m];
    for i in 0..m {
        for j in 0..m {
            let dt_ij = T::from_f64((s.trade_times[i] - s.trade_times[j]).abs());
            let ker = s.transient_impact * (-s.resilience_speed * dt_ij).exp();
            let t_cov = T::from_f64(s.trade_times[i.min(j)] - t0);
            let mut h_ij = ker + s.risk_aversion * s.volatility * s.volatility * t_cov * 2.0;
            if i == j {
                h_ij += s.temporary_impact * 2.0 + T::from_f64(1e-10);
            }
            a[i][j] = h_ij;
        }
    }
    for k in 0..m {
        let mut piv = k;
        for r in (k + 1)..m {
            if a[r][k].abs().val() > a[piv][k].abs().val() {
                piv = r;
            }
        }
        a.swap(k, piv);
        b.swap(k, piv);
        let diag = a[k][k];
        for r in (k + 1)..m {
            let f = a[r][k] / diag;
            b[r] = b[r] - f * b[k];
            for c in k..m {
                let akc = a[k][c];
                a[r][c] -= f * akc;
            }
        }
    }
    let mut y = vec![T::zero(); m];
    for i in (0..m).rev() {
        let mut sum = b[i];
        for j in (i + 1)..m {
            sum -= a[i][j] * y[j];
        }
        y[i] = sum / a[i][i];
    }
    let y_sum = y
        .iter()
        .copied()
        .fold(T::zero(), |acc, z| acc + z)
        .max_val(T::from_f64(1e-12));
    let trades: Vec<T> = y
        .into_iter()
        .map(|yi| s.initial_inventory * yi / y_sum)
        .collect();
    let mut inv = Vec::with_capacity(n);
    let mut left = s.initial_inventory;
    inv.push(left);
    for x in &trades {
        left -= *x;
        inv.push(left);
    }
    let mut cost = T::zero();
    for i in 0..m {
        cost += s.temporary_impact * trades[i] * trades[i];
        for j in 0..m {
            let dt_ij = T::from_f64((s.trade_times[i] - s.trade_times[j]).abs());
            cost += s.transient_impact
                * (-s.resilience_speed * dt_ij).exp()
                * trades[i]
                * trades[j]
                * 0.5;
        }
    }
    let mut var = T::zero();
    for k in 0..m {
        let dt_k = T::from_f64(s.trade_times[k + 1] - s.trade_times[k]);
        var += s.volatility * s.volatility * dt_k * inv[k + 1] * inv[k + 1];
    }
    Ok(OptimalExecutionReport {
        times: s.trade_times.clone(),
        inventory_schedule: inv,
        trade_schedule: trades,
        expected_cost: cost,
        cost_variance: var,
        mean_variance_utility: cost + s.risk_aversion * var,
    })
}
pub fn obizhaeva_wang_execution<T: Real>(s: &ObizhaevaWangSpec<T>) -> OptimalExecutionReport<T> {
    try_obizhaeva_wang_execution(s).expect("valid Obizhaeva-Wang inputs")
}

#[derive(Clone, Debug)]
pub struct MultiAssetAlmgrenChrissSpec<T: Real> {
    pub initial_inventory: Vec<T>,
    pub horizon: T,
    pub num_steps: usize,
    pub covariance: Vec<Vec<T>>,
    pub temporary_impact: Vec<Vec<T>>,
    pub permanent_impact: Vec<Vec<T>>,
    pub risk_aversion: T,
}
#[derive(Clone, Debug)]
pub struct MultiAssetExecutionReport<T: Real> {
    pub times: Vec<f64>,
    pub inventory_schedule: Vec<Vec<T>>,
    pub trade_schedule: Vec<Vec<T>>,
    pub expected_cost: T,
    pub cost_variance: T,
    pub mean_variance_utility: T,
}

pub fn try_multi_asset_almgren_chriss_execution<T: Real>(
    s: &MultiAssetAlmgrenChrissSpec<T>,
) -> Result<MultiAssetExecutionReport<T>, AllocationError> {
    let d = s.initial_inventory.len();
    if d == 0
        || s.num_steps == 0
        || !s.horizon.val().is_finite()
        || s.horizon.val() <= 0.0
        || s.covariance.len() != d
        || s.temporary_impact.len() != d
        || s.permanent_impact.len() != d
        || s.covariance.iter().any(|r| r.len() != d)
        || s.temporary_impact.iter().any(|r| r.len() != d)
        || s.permanent_impact.iter().any(|r| r.len() != d)
        || !s.risk_aversion.val().is_finite()
        || s.risk_aversion.val() < 0.0
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let n = s.num_steps;
    let dt = s.horizon / n as f64;
    let mut h_tilde = vec![vec![T::zero(); d]; d];
    for i in 0..d {
        for j in 0..d {
            let hs = (s.temporary_impact[i][j] + s.temporary_impact[j][i]) * 0.5;
            let gs = (s.permanent_impact[i][j] + s.permanent_impact[j][i]) * 0.5;
            h_tilde[i][j] = hs - gs * dt * 0.5;
        }
        if h_tilde[i][i].val() <= 0.0 || s.covariance[i][i].val() < 0.0 {
            return Err(AllocationError::InvalidDimensions);
        }
    }
    let mut inv = vec![vec![T::zero(); d]; n + 1];
    inv[0] = s.initial_inventory.clone();
    if n > 1 {
        let m = (n - 1) * d;
        let mut a = vec![vec![T::zero(); m]; m];
        let mut rhs = vec![T::zero(); m];
        for k in 0..(n - 1) {
            for i in 0..d {
                let row = k * d + i;
                for j in 0..d {
                    let cov_s = (s.covariance[i][j] + s.covariance[j][i]) * 0.5;
                    a[row][k * d + j] += h_tilde[i][j] * 2.0 + s.risk_aversion * cov_s * dt * dt;
                    if k > 0 {
                        a[row][(k - 1) * d + j] -= h_tilde[i][j];
                    }
                    if k + 1 < n - 1 {
                        a[row][(k + 1) * d + j] -= h_tilde[i][j];
                    }
                    if k == 0 {
                        rhs[row] += h_tilde[i][j] * s.initial_inventory[j];
                    }
                }
            }
        }
        let sol = crate::allocation::solve_linear_generic(&a, &rhs)?;
        for k in 1..n {
            for i in 0..d {
                inv[k][i] = sol[(k - 1) * d + i];
            }
        }
    }
    let mut trades = vec![vec![T::zero(); d]; n];
    let mut cost = T::zero();
    for i in 0..d {
        for j in 0..d {
            cost +=
                s.permanent_impact[i][j] * s.initial_inventory[i] * s.initial_inventory[j] * 0.5;
        }
    }
    for k in 0..n {
        for i in 0..d {
            trades[k][i] = inv[k][i] - inv[k + 1][i];
        }
        for i in 0..d {
            for j in 0..d {
                cost += h_tilde[i][j] * trades[k][i] * trades[k][j] / dt;
            }
        }
    }
    let mut var = T::zero();
    for row in inv.iter().take(n).skip(1) {
        for i in 0..d {
            for j in 0..d {
                var += s.covariance[i][j] * row[i] * row[j] * dt;
            }
        }
    }
    Ok(MultiAssetExecutionReport {
        times: (0..=n)
            .map(|k| (dt * T::from_f64(k as f64)).val())
            .collect(),
        inventory_schedule: inv,
        trade_schedule: trades,
        expected_cost: cost,
        cost_variance: var,
        mean_variance_utility: cost + s.risk_aversion * var,
    })
}

pub fn multi_asset_almgren_chriss_execution<T: Real>(
    s: &MultiAssetAlmgrenChrissSpec<T>,
) -> MultiAssetExecutionReport<T> {
    try_multi_asset_almgren_chriss_execution(s).expect("valid multi-asset execution inputs")
}

#[cfg(test)]
mod phase14_tests {
    use super::*;
    #[test]
    fn phase14_multi_asset_liquidates_inventory() {
        let s = MultiAssetAlmgrenChrissSpec {
            initial_inventory: vec![100.0, 50.0],
            horizon: 1.0,
            num_steps: 4,
            covariance: vec![vec![0.04, 0.0], vec![0.0, 0.09]],
            temporary_impact: vec![vec![1e-3, 0.0], vec![0.0, 1e-3]],
            permanent_impact: vec![vec![1e-4, 0.0], vec![0.0, 1e-4]],
            risk_aversion: 0.0,
        };
        let r = multi_asset_almgren_chriss_execution(&s);
        assert!(
            r.inventory_schedule
                .last()
                .unwrap()
                .iter()
                .all(|x| x.val().abs() < 1e-10)
        );
        assert!((r.times.last().copied().unwrap() - 1.0).abs() < 1e-12);
    }
}
