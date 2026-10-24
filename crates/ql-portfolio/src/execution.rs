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
    let q = s.risk_aversion * s.volatility * s.volatility * dt * dt
        / (s.temporary_impact.max_val(T::from_f64(1e-14)));
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
    if s.trade_times.len() < 2 {
        return Err(AllocationError::InvalidDimensions);
    }
    let n = s.trade_times.len();
    let trades = vec![s.initial_inventory / T::from_f64((n - 1) as f64); n - 1];
    let mut inv = Vec::with_capacity(n);
    let mut left = s.initial_inventory;
    inv.push(left);
    for x in &trades {
        left -= *x;
        inv.push(left)
    }
    let mut cost = T::zero();
    for x in &trades {
        cost += *x * *x * (s.temporary_impact + s.transient_impact)
    }
    let var =
        s.volatility * s.volatility * T::from_f64(s.trade_times.last().unwrap() - s.trade_times[0]);
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

// Phase 12 Round 4 audit coverage.
