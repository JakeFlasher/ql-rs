use super::allocation::AllocationError;
use ql_ad::Real;

#[derive(Copy, Clone, Debug)]
pub struct CarteaJaimungalAlphaExecutionSpec<T: Real = f64>{pub initial_inventory:T,pub horizon:T,pub num_steps:usize,pub temporary_impact_k:T,pub permanent_impact_b:T,pub inventory_penalty_phi:T,pub terminal_penalty_a:T,pub initial_alpha_signal:T,pub alpha_mean_reversion_kappa:T,pub alpha_volatility:T}
#[derive(Clone, Debug)]
pub struct CarteaJaimungalAlphaExecutionReport<T: Real = f64>{pub times:Vec<T>,pub inventory_trajectory:Vec<T>,pub baseline_ac_inventory_trajectory:Vec<T>,pub optimal_trading_rates:Vec<T>,pub alpha_speculative_tilt_rates:Vec<T>,pub expected_alpha_capture_pnl:T,pub expected_execution_cost_net_of_alpha:T,pub value_function_h0_alpha_variance_bonus:T}
pub fn solve_cartea_jaimungal_alpha_execution<T: Real>(s:&CarteaJaimungalAlphaExecutionSpec<T>)->CarteaJaimungalAlphaExecutionReport<T>{
 let n=s.num_steps.max(4);let dt=s.horizon/T::from_f64(n as f64);let k=s.temporary_impact_k.max_val(T::from_f64(1e-8));let phi=s.inventory_penalty_phi.max_val(T::from_f64(1e-8));let gam=(phi/k).sqrt();let kap=s.alpha_mean_reversion_kappa.max_val(T::from_f64(1e-6));let den=(gam*s.horizon).sinh().max_val(T::from_f64(1e-12));let mut times=Vec::with_capacity(n+1);let mut ac=Vec::with_capacity(n+1);let mut q=Vec::with_capacity(n+1);for j in 0..=n{let tt=dt*(j as f64);let qa=s.initial_inventory*(gam*(s.horizon-tt)).sinh()/den;let tilt=s.initial_alpha_signal/(k*den*2.0)*(gam*(s.horizon-tt)).sinh()*(T::one()-(-kap*tt).exp())/kap;times.push(tt);ac.push(qa);q.push(qa+tilt);}
 let mut rates=Vec::with_capacity(n);let mut tilt_rates=Vec::with_capacity(n);let mut cap=T::zero();let mut cost=T::zero();for j in 0..n{let nu=-(q[j+1]-q[j])/dt;let nua=-(ac[j+1]-ac[j])/dt;rates.push(nu);tilt_rates.push(nu-nua);let t=times[j];let alpha=s.initial_alpha_signal*(-kap*t).exp();cap+=alpha*q[j]*dt;cost+=(k*nu*nu+phi*q[j]*q[j])*dt;}cap+=s.alpha_volatility*s.alpha_volatility/(k*kap*kap*4.0)*s.horizon;cost+=s.permanent_impact_b*0.5*s.initial_inventory*s.initial_inventory-cap;CarteaJaimungalAlphaExecutionReport{times,inventory_trajectory:q,baseline_ac_inventory_trajectory:ac,optimal_trading_rates:rates,alpha_speculative_tilt_rates:tilt_rates,expected_alpha_capture_pnl:cap,expected_execution_cost_net_of_alpha:cost,value_function_h0_alpha_variance_bonus:s.alpha_volatility*s.alpha_volatility/(k*kap*kap*4.0)*s.horizon}
}
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

#[derive(Copy, Clone, Debug)]
pub struct AvellanedaStoikovMmSpec<T: Real = f64> {
    pub mid_price: T,
    pub initial_inventory: T,
    pub horizon: T,
    pub num_steps: usize,
    pub volatility: T,
    pub risk_aversion: T,
    pub order_book_kappa: T,
    pub base_arrival_intensity: T,
    pub alpha_drift: T,
}

#[derive(Clone, Debug)]
pub struct AvellanedaStoikovMmReport<T: Real> {
    pub reservation_price: T,
    pub optimal_bid_quote: T,
    pub optimal_ask_quote: T,
    pub optimal_bid_ask_spread: T,
    pub bid_fill_intensity: T,
    pub ask_fill_intensity: T,
    pub expected_inventory_trajectory: Vec<T>,
    pub expected_spread_capture_pnl: T,
}

pub fn solve_avellaneda_stoikov_market_making<T: Real>(
    spec: &AvellanedaStoikovMmSpec<T>,
) -> AvellanedaStoikovMmReport<T> {
    let n = spec.num_steps.max(1);
    let tau = spec.horizon.max_val(T::from_f64(1e-6));
    let dt = tau / n as f64;
    let gamma = spec.risk_aversion.max_val(T::from_f64(1e-8));
    let kappa = spec.order_book_kappa.max_val(T::from_f64(1e-8));
    let sig2 = spec.volatility * spec.volatility;

    // Reservation (indifference) price at t = 0: r(s, q, T) = s + alpha * T - q * gamma * sigma^2 * T
    let res_price = spec.mid_price + spec.alpha_drift * tau - spec.initial_inventory * gamma * sig2 * tau;
    let base_half_spread = (T::one() + gamma / kappa).ln() / kappa;
    let half_spread = base_half_spread + gamma * sig2 * tau * 0.5;
    let opt_bid = res_price - half_spread;
    let opt_ask = res_price + half_spread;
    let opt_spread = opt_ask - opt_bid;

    let delta_bid0 = spec.mid_price - opt_bid;
    let delta_ask0 = opt_ask - spec.mid_price;
    let lam_bid0 = spec.base_arrival_intensity * (-kappa * delta_bid0).exp();
    let lam_ask0 = spec.base_arrival_intensity * (-kappa * delta_ask0).exp();

    let mut q = spec.initial_inventory;
    let mut traj = Vec::with_capacity(n + 1);
    traj.push(q);
    let mut pnl = T::zero();

    for step in 0..n {
        let rem = tau * ((n - step) as f64 / n as f64);
        let r_t = spec.mid_price + spec.alpha_drift * rem - q * gamma * sig2 * rem;
        let hs_t = base_half_spread + gamma * sig2 * rem * 0.5;
        let db = (spec.mid_price - (r_t - hs_t))
            .max_val(T::from_f64(-20.0))
            .min_val(T::from_f64(20.0));
        let da = ((r_t + hs_t) - spec.mid_price)
            .max_val(T::from_f64(-20.0))
            .min_val(T::from_f64(20.0));
        let lb = spec.base_arrival_intensity * (-kappa * db).exp();
        let la = spec.base_arrival_intensity * (-kappa * da).exp();
        pnl += (lb * db + la * da + q * spec.alpha_drift - q * q * gamma * sig2 * 0.5) * dt;
        q += (lb - la) * dt;
        traj.push(q);
    }

    AvellanedaStoikovMmReport {
        reservation_price: res_price,
        optimal_bid_quote: opt_bid,
        optimal_ask_quote: opt_ask,
        optimal_bid_ask_spread: opt_spread,
        bid_fill_intensity: lam_bid0,
        ask_fill_intensity: lam_ask0,
        expected_inventory_trajectory: traj,
        expected_spread_capture_pnl: pnl,
    }
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

#[cfg(test)]
mod phase19_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn avellaneda_stoikov_inventory_skew_symmetry_and_ad() {
        let zero_q = AvellanedaStoikovMmSpec {
            mid_price: 100.0,
            initial_inventory: 0.0,
            horizon: 1.0,
            num_steps: 10,
            volatility: 0.20,
            risk_aversion: 0.5,
            order_book_kappa: 1.5,
            base_arrival_intensity: 10.0,
            alpha_drift: 0.0,
        };
        let rep0 = solve_avellaneda_stoikov_market_making(&zero_q);
        assert!((rep0.reservation_price - 100.0).abs() < 1e-12);
        assert!((rep0.bid_fill_intensity - rep0.ask_fill_intensity).abs() < 1e-12);

        let mut long_q = zero_q;
        long_q.initial_inventory = 5.0;
        let rep_long = solve_avellaneda_stoikov_market_making(&long_q);
        // Long inventory skews reservation price downward and increases ask fill intensity vs bid
        assert!(rep_long.reservation_price < 100.0);
        assert!(rep_long.ask_fill_intensity > rep_long.bid_fill_intensity);
        assert!(rep_long.expected_inventory_trajectory.last().copied().unwrap() < 5.0);

        let h = 1e-5;
        let eval_q = |q: f64| {
            let mut s = long_q;
            s.initial_inventory = q;
            solve_avellaneda_stoikov_market_making(&s).optimal_bid_quote
        };
        let fd = (eval_q(5.0 + h) - eval_q(5.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let spec_d = AvellanedaStoikovMmSpec {
            mid_price: D::from_value(100.0),
            initial_inventory: D::variable(5.0, 0),
            horizon: D::one(),
            num_steps: 10,
            volatility: D::from_value(0.20),
            risk_aversion: D::from_value(0.5),
            order_book_kappa: D::from_value(1.5),
            base_arrival_intensity: D::from_value(10.0),
            alpha_drift: D::zero(),
        };
        let ad = solve_avellaneda_stoikov_market_making(&spec_d).optimal_bid_quote;
        assert!(ad.d[0] < 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}


#[derive(Copy, Clone, Debug)]
pub struct AlfonsiFruthSchiedLobSpec<T: Real = f64> {
    pub total_order_size: T,
    pub horizon: T,
    pub num_intervals: usize,
    pub lob_density_q: T,
    /// Power-law shape exponent $\alpha \le 0.95$ in $f(x) = q / (|x| + 1)^\alpha$ ($\alpha = 0$ is flat Block LOB).
    pub shape_exponent_alpha: T,
    pub resilience_rho: T,
}

#[derive(Clone, Debug)]
pub struct AlfonsiFruthSchiedLobReport<T: Real = f64> {
    pub volume_resilience_trades: Vec<T>,
    pub spread_resilience_trades: Vec<T>,
    pub volume_resilience_cost: T,
    pub spread_resilience_cost: T,
    pub twap_cost_volume_model: T,
    pub cost_savings_vs_twap: T,
    pub initial_extra_spread_d0: T,
}

fn afs_f_density<T: Real>(q: T, alpha: T, x: T) -> T {
    let xp1 = (x + T::one()).max_val(T::from_f64(1e-8));
    q / xp1.pow(alpha)
}

fn afs_f_cum<T: Real>(q: T, alpha: T, x: T) -> T {
    let xp1 = (x + T::one()).max_val(T::from_f64(1e-8));
    let oma = T::one() - alpha;
    if oma.val().abs() < 1e-6 {
        q * xp1.ln()
    } else {
        q * (xp1.pow(oma) - T::one()) / oma
    }
}

fn afs_f_inv<T: Real>(q: T, alpha: T, y: T) -> T {
    let y_pos = y.max_val(T::zero());
    let oma = T::one() - alpha;
    if oma.val().abs() < 1e-6 {
        (y_pos / q).exp() - T::one()
    } else {
        let base = (T::one() + oma * y_pos / q).max_val(T::from_f64(1e-12));
        base.pow(T::one() / oma) - T::one()
    }
}

fn afs_f_tilde_cost<T: Real>(q: T, alpha: T, x: T) -> T {
    let xp1 = (x + T::one()).max_val(T::from_f64(1e-8));
    let oma = T::one() - alpha;
    let tma = T::from_f64(2.0) - alpha;
    if oma.val().abs() < 1e-6 {
        q * (x - xp1.ln())
    } else {
        q * ((xp1.pow(tma) - T::one()) / tma - (xp1.pow(oma) - T::one()) / oma)
    }
}

fn afs_g_cost<T: Real>(q: T, alpha: T, y: T) -> T {
    afs_f_tilde_cost(q, alpha, afs_f_inv(q, alpha, y))
}

/// Alfonsi, Fruth & Schied (2010) Optimal Execution in Limit Order Books with General
/// Power-Law Shape Functions $f(x) = q / (|x| + 1)^\alpha$ (`papers/55_Alfonsi_Fruth_Schied_2010_Optimal_Execution_Transient_Impact.pdf`, Theorems 4.1 & 5.1).
pub fn solve_alfonsi_fruth_schied_lob_execution<T: Real>(
    spec: &AlfonsiFruthSchiedLobSpec<T>,
) -> AlfonsiFruthSchiedLobReport<T> {
    let n_int = spec.num_intervals.max(1);
    let n_f = n_int as f64;
    let tau = spec.horizon.max_val(T::from_f64(1e-6)) / n_f;
    let a = (-spec.resilience_rho.max_val(T::from_f64(1e-6)) * tau).exp();
    let one_m_a = (T::one() - a).max_val(T::from_f64(1e-8));
    let q = spec.lob_density_q.max_val(T::from_f64(1e-6));
    let alpha = spec.shape_exponent_alpha.min_val(T::from_f64(0.95));
    let x0 = spec.total_order_size.max_val(T::from_f64(1e-8));

    // Block-LOB exact solution (Corollary 6.1) used as starting point for Newton
    let xi0_block = x0 / (one_m_a * (n_f - 1.0) + T::from_f64(2.0));

    // Model 1 (Volume Resilience): solve h_1(y) = F^{-1}(y) - a F^{-1}(a y) - (1-a) F^{-1}(X_0 - N(1-a)y) = 0
    let mut y1 = xi0_block;
    for _ in 0..8 {
        let rem_vol = (x0 - one_m_a * y1 * n_f).max_val(T::from_f64(1e-8));
        let inv_y = afs_f_inv(q, alpha, y1);
        let inv_ay = afs_f_inv(q, alpha, a * y1);
        let inv_rem = afs_f_inv(q, alpha, rem_vol);
        let h1 = inv_y - a * inv_ay - one_m_a * inv_rem;
        let dh1 = T::one() / afs_f_density(q, alpha, inv_y)
            - (a * a) / afs_f_density(q, alpha, inv_ay)
            + (one_m_a * one_m_a * n_f) / afs_f_density(q, alpha, inv_rem);
        y1 = (y1 - h1 / dh1.max_val(T::from_f64(1e-12)))
            .max_val(x0 * 0.01)
            .min_val(x0 / (one_m_a * n_f + T::from_f64(0.05)));
    }
    let mut vol_trades = Vec::with_capacity(n_int + 1);
    vol_trades.push(y1);
    for _ in 1..n_int {
        vol_trades.push(one_m_a * y1);
    }
    let xi_n1 = x0 - y1 - one_m_a * y1 * (n_f - 1.0);
    vol_trades.push(xi_n1);
    let vol_cost = (afs_g_cost(q, alpha, y1) - afs_g_cost(q, alpha, a * y1)) * n_f
        + afs_g_cost(q, alpha, x0 - one_m_a * y1 * n_f);

    // Model 2 (Spread Resilience): solve h_2(d_0) = d_0 (f(d_0) - a^2 f(a d_0))/(f(d_0) - a f(a d_0)) - F^{-1}(X_0 - N(F(d_0) - F(a d_0))) = 0
    // Model 2 (Spread Resilience): solve h_2(d_0) = d_0 (f(d_0) - a^2 f(a d_0))/(f(d_0) - a f(a d_0)) - F^{-1}(X_0 - N(F(d_0) - F(a d_0))) = 0
    let mut d0 = afs_f_inv(q, alpha, xi0_block);
    for _ in 0..8 {
        let fd = afs_f_density(q, alpha, d0);
        let fad = afs_f_density(q, alpha, a * d0);
        let dfd = -alpha * q * (d0 + T::one()).pow(-alpha - T::one());
        let dfad = -alpha * q * (a * d0 + T::one()).pow(-alpha - T::one());
        let u_num = d0 * (fd - a * a * fad);
        let v_den = (fd - a * fad).max_val(T::from_f64(1e-12));
        let du_dd = (fd - a * a * fad) + d0 * (dfd - a * a * a * dfad);
        let dv_dd = dfd - a * a * dfad;
        let d_ratio = (du_dd * v_den - u_num * dv_dd) / (v_den * v_den);

        let step_vol = afs_f_cum(q, alpha, d0) - afs_f_cum(q, alpha, a * d0);
        let rem = (x0 - step_vol * n_f).max_val(T::from_f64(1e-8));
        let inv_rem = afs_f_inv(q, alpha, rem);
        let r0 = u_num / v_den - inv_rem;
        let d_inv_rem = -T::from_f64(n_f) * v_den
            / afs_f_density(q, alpha, inv_rem).max_val(T::from_f64(1e-12));
        let dr = d_ratio - d_inv_rem;
        d0 = (d0 - r0 / dr.max_val(T::from_f64(1e-12))).max_val(T::from_f64(1e-6));
    }
    let f_d0 = afs_f_cum(q, alpha, d0);
    let f_ad0 = afs_f_cum(q, alpha, a * d0);
    let mid_trade2 = f_d0 - f_ad0;
    let mut spread_trades = Vec::with_capacity(n_int + 1);
    spread_trades.push(f_d0);
    for _ in 1..n_int {
        spread_trades.push(mid_trade2);
    }
    spread_trades.push(x0 - f_d0 - mid_trade2 * (n_f - 1.0));
    let spread_cost = (afs_f_tilde_cost(q, alpha, d0) - afs_f_tilde_cost(q, alpha, a * d0)) * n_f
        + afs_g_cost(q, alpha, x0 - mid_trade2 * n_f);

    // TWAP benchmark cost under Volume Resilience
    let twap_slice = x0 / (n_f + 1.0);
    let mut e_twap = T::zero();
    let mut twap_cost = T::zero();
    for k in 0..=n_int {
        let e_after = e_twap + twap_slice;
        twap_cost += afs_g_cost(q, alpha, e_after) - afs_g_cost(q, alpha, e_twap);
        if k < n_int {
            e_twap = a * e_after;
        }
    }

    AlfonsiFruthSchiedLobReport {
        volume_resilience_trades: vol_trades,
        spread_resilience_trades: spread_trades,
        volume_resilience_cost: vol_cost,
        spread_resilience_cost: spread_cost,
        twap_cost_volume_model: twap_cost,
        cost_savings_vs_twap: twap_cost - vol_cost,
        initial_extra_spread_d0: d0,
    }
}



#[cfg(test)]
mod phase20_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn alfonsi_fruth_schied_block_collapse_twap_outperformance_and_ad() {
        // When shape_exponent_alpha == 0.0 (Block LOB), Model 1 and Model 2 coincide and match Obizhaeva-Wang
        let block_spec = AlfonsiFruthSchiedLobSpec {
            total_order_size: 100.0,
            horizon: 1.0,
            num_intervals: 4,
            lob_density_q: 50.0,
            shape_exponent_alpha: 0.0,
            resilience_rho: 1.2,
        };
        let r_block = solve_alfonsi_fruth_schied_lob_execution(&block_spec);
        for k in 0..=4 {
            assert!(
                (r_block.volume_resilience_trades[k] - r_block.spread_resilience_trades[k]).abs()
                    < 1e-10
            );
        }
        assert!(
            (r_block.volume_resilience_trades[0] - r_block.volume_resilience_trades[4]).abs()
                < 1e-10
        );

        // Power-law decaying LOB (alpha = 0.5): strictly outperforms TWAP and satisfies volume conservation
        let power_spec = AlfonsiFruthSchiedLobSpec {
            shape_exponent_alpha: 0.5,
            ..block_spec
        };
        let r_pow = solve_alfonsi_fruth_schied_lob_execution(&power_spec);
        let sum_v: f64 = r_pow.volume_resilience_trades.iter().sum();
        let sum_s: f64 = r_pow.spread_resilience_trades.iter().sum();
        assert!((sum_v - 100.0).abs() < 1e-12);
        assert!((sum_s - 100.0).abs() < 1e-12);
        assert!(r_pow.cost_savings_vs_twap > 0.0);

        let h = 1e-5;
        let eval_rho = |rho: f64| {
            let mut s = power_spec;
            s.resilience_rho = rho;
            solve_alfonsi_fruth_schied_lob_execution(&s).volume_resilience_cost
        };
        let fd = (eval_rho(1.2 + h) - eval_rho(1.2 - h)) / (2.0 * h);
        type D = Dual<1>;
        let spec_d = AlfonsiFruthSchiedLobSpec {
            total_order_size: D::from_value(100.0),
            horizon: D::one(),
            num_intervals: 4,
            lob_density_q: D::from_value(50.0),
            shape_exponent_alpha: D::from_value(0.5),
            resilience_rho: D::variable(1.2, 0),
        };
        let ad = solve_alfonsi_fruth_schied_lob_execution(&spec_d).volume_resilience_cost;
        assert!(ad.d[0] < 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
