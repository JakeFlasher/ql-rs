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
            let t_cov = T::from_f64(s.trade_times[i.min(j) + 1] - t0);
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
        var += s.volatility * s.volatility * dt_k * inv[k] * inv[k];
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
