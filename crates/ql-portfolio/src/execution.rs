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
    let mut g_sym = vec![vec![T::zero(); d]; d];
    let mut g_skew = vec![vec![T::zero(); d]; d];
    for i in 0..d {
        for j in 0..d {
            let hs = (s.temporary_impact[i][j] + s.temporary_impact[j][i]) * 0.5;
            let gs = (s.permanent_impact[i][j] + s.permanent_impact[j][i]) * 0.5;
            h_tilde[i][j] = hs - gs * dt * 0.5;
            g_sym[i][j] = gs;
            g_skew[i][j] = (s.permanent_impact[i][j] - s.permanent_impact[j][i]) * 0.5;
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
                    let sub_ij = h_tilde[i][j] + g_skew[i][j] * dt * 0.5;
                    let sup_ij = h_tilde[i][j] - g_skew[i][j] * dt * 0.5;
                    a[row][k * d + j] += h_tilde[i][j] * 2.0 + s.risk_aversion * cov_s * dt * dt;
                    if k > 0 {
                        a[row][(k - 1) * d + j] -= sub_ij;
                    }
                    if k + 1 < n - 1 {
                        a[row][(k + 1) * d + j] -= sup_ij;
                    }
                    if k == 0 {
                        rhs[row] += sub_ij * s.initial_inventory[j];
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
            cost += g_sym[i][j] * s.initial_inventory[i] * s.initial_inventory[j] * 0.5;
        }
    }
    for k in 0..n {
        for i in 0..d {
            trades[k][i] = inv[k][i] - inv[k + 1][i];
        }
        for i in 0..d {
            for j in 0..d {
                cost += h_tilde[i][j] * trades[k][i] * trades[k][j] / dt
                    + g_skew[i][j] * trades[k][i] * (inv[k][j] + inv[k + 1][j]) * 0.5;
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
    fn phase14_multi_asset_matches_1d_and_twap() {
        let one_d = almgren_chriss_execution(&AlmgrenChrissSpec {
            initial_inventory: 1000.0,
            horizon: 1.0,
            num_steps: 5,
            volatility: 0.2,
            temporary_impact: 1e-4,
            permanent_impact: 5e-5,
            risk_aversion: 1e-4,
        });
        let multi_1d = multi_asset_almgren_chriss_execution(&MultiAssetAlmgrenChrissSpec {
            initial_inventory: vec![1000.0],
            horizon: 1.0,
            num_steps: 5,
            covariance: vec![vec![0.04]],
            temporary_impact: vec![vec![1e-4]],
            permanent_impact: vec![vec![5e-5]],
            risk_aversion: 1e-4,
        });
        for k in 0..=5 {
            assert!((one_d.inventory_schedule[k] - multi_1d.inventory_schedule[k][0]).abs() < 1e-9);
        }
        assert!((one_d.expected_cost - multi_1d.expected_cost).abs() < 1e-9);
        assert!((one_d.cost_variance - multi_1d.cost_variance).abs() < 1e-9);
    }
}

#[derive(Copy, Clone, Debug)]
pub struct DarkLitRoutingSpec<T: Real = f64> {
    pub initial_inventory: T,
    pub horizon: T,
    pub num_steps: usize,
    pub volatility: T,
    pub lit_temporary_impact: T,
    pub lit_permanent_impact: T,
    pub dark_fill_intensity: T,
    pub risk_aversion: T,
}

#[derive(Clone, Debug)]
pub struct DarkLitRoutingReport<T: Real> {
    pub times: Vec<f64>,
    pub pre_fill_inventory_schedule: Vec<T>,
    pub expected_inventory_schedule: Vec<T>,
    pub lit_trading_rate_schedule: Vec<T>,
    pub dark_pool_order_schedule: Vec<T>,
    pub expected_total_cost: T,
}

pub fn try_optimal_dark_lit_execution<T: Real>(
    s: &DarkLitRoutingSpec<T>,
) -> Result<DarkLitRoutingReport<T>, AllocationError> {
    if s.num_steps == 0
        || !s.initial_inventory.val().is_finite()
        || !s.horizon.val().is_finite()
        || s.horizon.val() <= 0.0
        || !s.volatility.val().is_finite()
        || s.volatility.val() < 0.0
        || !s.lit_temporary_impact.val().is_finite()
        || s.lit_temporary_impact.val() <= 0.0
        || !s.lit_permanent_impact.val().is_finite()
        || !s.dark_fill_intensity.val().is_finite()
        || s.dark_fill_intensity.val() < 0.0
        || !s.risk_aversion.val().is_finite()
        || s.risk_aversion.val() < 0.0
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let n = s.num_steps;
    let t_cap = s.horizon;
    let eta = s.lit_temporary_impact;
    let theta = s.dark_fill_intensity;
    let phi = s.risk_aversion * s.volatility * s.volatility;
    let half_theta = theta * 0.5;
    let gamma = (half_theta * half_theta + phi / eta)
        .max_val(T::from_f64(1e-16))
        .sqrt();
    let gt = gamma * t_cap;

    let mut times = Vec::with_capacity(n + 1);
    let mut pre_fill = Vec::with_capacity(n + 1);
    let mut exp_inv = Vec::with_capacity(n + 1);
    let mut lit_rate = Vec::with_capacity(n + 1);
    let mut dark_orders = Vec::with_capacity(n + 1);

    for k in 0..=n {
        let frac = (k as f64) / (n as f64);
        let t_f64 = s.horizon.val() * frac;
        times.push(t_f64);
        let t = t_cap * frac;
        let rem = t_cap * (1.0 - frac);
        let (ratio_sinh, rate_factor) = if gt.val().abs() < 1e-4 {
            let r_s = T::from_f64(1.0 - frac)
                * (T::one() + gamma * gamma * (rem * rem - t_cap * t_cap) / 6.0);
            let r_v = (T::one() / t_cap)
                * (T::one() + gamma * gamma * (rem * rem * 3.0 - t_cap * t_cap) / 6.0)
                - half_theta * r_s;
            (r_s, r_v)
        } else {
            let sh_gt = ((gt).exp() - (-gt).exp()) * 0.5;
            let g_rem = gamma * rem;
            let sh_rem = (g_rem.exp() - (-g_rem).exp()) * 0.5;
            let ch_rem = (g_rem.exp() + (-g_rem).exp()) * 0.5;
            let r_s = sh_rem / sh_gt;
            let r_v = (gamma * ch_rem - half_theta * sh_rem) / sh_gt;
            (r_s, r_v)
        };
        let x_pre = s.initial_inventory * (half_theta * t).exp() * ratio_sinh;
        let x_exp = s.initial_inventory * (-half_theta * t).exp() * ratio_sinh;
        let v_lit = s.initial_inventory * (half_theta * t).exp() * rate_factor;
        pre_fill.push(x_pre);
        exp_inv.push(x_exp);
        lit_rate.push(v_lit);
        dark_orders.push(x_pre);
    }

    let coth_term = if gt.val().abs() < 1e-4 {
        (T::one() / t_cap) * (T::one() + gt * gt / 3.0 - gt * gt * gt * gt / 45.0)
    } else {
        let e2 = (gt * 2.0).exp();
        gamma * (e2 + T::one()) / (e2 - T::one())
    };
    let perm_cost = s.lit_permanent_impact * s.initial_inventory * s.initial_inventory * 0.5;
    let expected_total_cost =
        eta * s.initial_inventory * s.initial_inventory * (coth_term - half_theta) + perm_cost;

    Ok(DarkLitRoutingReport {
        times,
        pre_fill_inventory_schedule: pre_fill,
        expected_inventory_schedule: exp_inv,
        lit_trading_rate_schedule: lit_rate,
        dark_pool_order_schedule: dark_orders,
        expected_total_cost,
    })
}

pub fn optimal_dark_lit_execution<T: Real>(s: &DarkLitRoutingSpec<T>) -> DarkLitRoutingReport<T> {
    try_optimal_dark_lit_execution(s).expect("valid dark-lit execution inputs")
}

#[cfg(test)]
mod phase17_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn kratz_schoeneborn_dark_pool_monotonicity_delay_and_ad() {
        let mut s = DarkLitRoutingSpec {
            initial_inventory: 100.0,
            horizon: 1.0,
            num_steps: 4,
            volatility: 0.20,
            lit_temporary_impact: 1e-3,
            lit_permanent_impact: 2e-4,
            dark_fill_intensity: 0.0,
            risk_aversion: 0.0,
        };
        let no_dark = optimal_dark_lit_execution(&s);
        s.dark_fill_intensity = 1.5;
        let with_dark = optimal_dark_lit_execution(&s);
        assert!(with_dark.expected_total_cost < no_dark.expected_total_cost);
        assert!(with_dark.lit_trading_rate_schedule[0] < with_dark.lit_trading_rate_schedule[4]);

        let h = 1e-5;
        let eval_theta = |th: f64| {
            let mut sp = s;
            sp.dark_fill_intensity = th;
            optimal_dark_lit_execution(&sp).expected_total_cost
        };
        let fd = (eval_theta(1.5 + h) - eval_theta(1.5 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sd = DarkLitRoutingSpec {
            initial_inventory: D::from_value(100.0),
            horizon: D::one(),
            num_steps: 4,
            volatility: D::from_value(0.20),
            lit_temporary_impact: D::from_value(1e-3),
            lit_permanent_impact: D::from_value(2e-4),
            dark_fill_intensity: D::variable(1.5, 0),
            risk_aversion: D::from_value(0.25),
        };
        let eval_theta_ra = |th: f64| {
            let mut sp = s;
            sp.risk_aversion = 0.25;
            sp.dark_fill_intensity = th;
            optimal_dark_lit_execution(&sp).expected_total_cost
        };
        let fd_ra = (eval_theta_ra(1.5 + h) - eval_theta_ra(1.5 - h)) / (2.0 * h);
        let ad = optimal_dark_lit_execution(&sd).expected_total_cost;
        assert!(fd < 0.0 && ad.d[0] < 0.0);
        assert!((ad.d[0] - fd_ra).abs() < 1e-6);
    }
}
