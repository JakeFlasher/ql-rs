#![allow(clippy::needless_range_loop)]
use super::allocation::AllocationError;
use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct HarrisRollCorwinSchultzMicrostructureSpec<T: Real = f64> {
    pub trade_prices: Vec<T>,
    pub daily_highs: Vec<T>,
    pub daily_lows: Vec<T>,
    pub daily_dollar_volumes: Vec<T>,
    pub decision_price: T,
    pub arrival_price: T,
    pub child_fill_quantities: Vec<T>,
    pub child_fill_prices: Vec<T>,
    pub total_order_quantity: T,
    pub post_trade_benchmark_price: T,
    pub explicit_commissions: T,
}
#[derive(Clone, Debug)]
pub struct HarrisRollCorwinSchultzMicrostructureReport<T: Real = f64> {
    pub roll_serial_cov_effective_spread: T,
    pub corwin_schultz_high_low_spread: T,
    pub parkinson_high_low_volatility: T,
    pub amihud_illiquidity_ratio: T,
    pub perold_delay_cost: T,
    pub perold_market_impact_cost: T,
    pub perold_explicit_fee_cost: T,
    pub perold_opportunity_cost: T,
    pub perold_total_implementation_shortfall: T,
    pub perold_identity_conservation_error: T,
}
pub fn analyze_harris_roll_corwin_schultz_microstructure<T: Real>(
    s: &HarrisRollCorwinSchultzMicrostructureSpec<T>,
) -> HarrisRollCorwinSchultzMicrostructureReport<T> {
    let mut diffs = Vec::new();
    for w in s.trade_prices.windows(2) {
        diffs.push(w[1] - w[0]);
    }
    let cov = if diffs.len() > 1 {
        diffs
            .windows(2)
            .map(|w| w[0] * w[1])
            .fold(T::zero(), |a, b| a + b)
            / T::from_f64((diffs.len() - 1) as f64)
    } else {
        T::from_f64(-1e-8)
    };
    let roll = T::from_f64(2.0) * (-cov).max_val(T::from_f64(1e-16)).sqrt();
    let d_len = s.daily_highs.len().min(s.daily_lows.len());
    let mut cs = T::zero();
    if d_len > 1 {
        for i in 1..d_len {
            let b = (s.daily_highs[i - 1] / s.daily_lows[i - 1].max_val(T::from_f64(1e-8)))
                .ln()
                .powf(2.0)
                + (s.daily_highs[i] / s.daily_lows[i].max_val(T::from_f64(1e-8)))
                    .ln()
                    .powf(2.0);
            let g = (s.daily_highs[i - 1].max_val(s.daily_highs[i])
                / s.daily_lows[i - 1]
                    .min_val(s.daily_lows[i])
                    .max_val(T::from_f64(1e-8)))
            .ln()
            .powf(2.0);
            let den = T::from_f64(3.0 - 2.0 * 2.0_f64.sqrt());
            let al = T::from_f64(2.0_f64.sqrt() - 1.0) * b.sqrt() / den - (g / den).sqrt();
            let e = al.max_val(T::zero()).exp();
            cs += T::from_f64(2.0) * (e - T::one()) / (e + T::one());
        }
        cs /= T::from_f64((d_len - 1) as f64);
    }
    let mut sum_hl2 = T::zero();
    let mut amihud = T::zero();
    for i in 0..d_len {
        let hl = (s.daily_highs[i] / s.daily_lows[i].max_val(T::from_f64(1e-8)))
            .ln()
            .abs();
        sum_hl2 += hl * hl;
        let dv = s
            .daily_dollar_volumes
            .get(i)
            .copied()
            .unwrap_or(T::one())
            .max_val(T::from_f64(1e-8));
        amihud += hl / dv;
    }
    let park = (sum_hl2 / T::from_f64(4.0 * d_len.max(1) as f64 * 2.0_f64.ln())).sqrt();
    let amihud = amihud / T::from_f64(d_len.max(1) as f64);
    let qexec = s
        .child_fill_quantities
        .iter()
        .copied()
        .fold(T::zero(), |a, b| a + b);
    let impact = s
        .child_fill_quantities
        .iter()
        .zip(s.child_fill_prices.iter())
        .map(|(&q, &p)| q * (p - s.arrival_price))
        .fold(T::zero(), |a, b| a + b);
    let delay = qexec * (s.arrival_price - s.decision_price);
    let opp = (s.total_order_quantity - qexec) * (s.post_trade_benchmark_price - s.decision_price);
    let fees = s.explicit_commissions;
    let total = delay + impact + fees + opp;
    let paper_ret = s.total_order_quantity * (s.post_trade_benchmark_price - s.decision_price);
    let actual_ret = s
        .child_fill_quantities
        .iter()
        .zip(s.child_fill_prices.iter())
        .map(|(&q, &p)| q * (s.post_trade_benchmark_price - p))
        .fold(T::zero(), |a, b| a + b)
        - s.explicit_commissions;
    let sum_parts = delay + impact + fees + opp;
    HarrisRollCorwinSchultzMicrostructureReport {
        roll_serial_cov_effective_spread: roll,
        corwin_schultz_high_low_spread: cs,
        parkinson_high_low_volatility: park,
        amihud_illiquidity_ratio: amihud,
        perold_delay_cost: delay,
        perold_market_impact_cost: impact,
        perold_explicit_fee_cost: fees,
        perold_opportunity_cost: opp,
        perold_total_implementation_shortfall: total,
        perold_identity_conservation_error: (paper_ret - actual_ret - sum_parts).abs(),
    }
}

#[derive(Clone, Debug)]
pub struct EngleRussellAcdTradeDurationSpec<T: Real = f64> {
    pub inter_trade_durations: Vec<T>,
    pub omega: T,
    pub alpha: T,
    pub beta: T,
    pub weibull_shape_gamma: T,
    pub initial_expected_duration: T,
    pub elapsed_time_since_last_trade: T,
    pub base_price_impact_bps: T,
    pub duration_impact_elasticity: T,
}
#[derive(Clone, Debug)]
pub struct EngleRussellAcdTradeDurationReport<T: Real = f64> {
    pub conditional_durations_psi: Vec<T>,
    pub standardized_durations: Vec<T>,
    pub acd_persistence: T,
    pub unconditional_mean_duration: T,
    pub next_expected_duration_psi_n1: T,
    pub eacd_log_likelihood: T,
    pub wacd_instantaneous_hazard_rate: T,
    pub duration_adjusted_price_impact_bps: T,
}
pub fn evaluate_engle_russell_acd_trade_duration_and_hazard<T: Real>(
    s: &EngleRussellAcdTradeDurationSpec<T>,
) -> EngleRussellAcdTradeDurationReport<T> {
    let mut psi = Vec::with_capacity(s.inter_trade_durations.len());
    let mut std = Vec::with_capacity(s.inter_trade_durations.len());
    let mut prev = s.initial_expected_duration.max_val(T::from_f64(1e-8));
    let mut ll = T::zero();
    for (i, &x) in s.inter_trade_durations.iter().enumerate() {
        let p = (s.omega
            + s.alpha
                * if i == 0 {
                    s.initial_expected_duration
                } else {
                    s.inter_trade_durations[i - 1]
                }
            + s.beta * prev)
            .max_val(T::from_f64(1e-8));
        psi.push(p);
        std.push(x / p);
        ll += -p.ln() - x / p;
        prev = p;
    }
    let pnext = (s.omega
        + s.alpha
            * s.inter_trade_durations
                .last()
                .copied()
                .unwrap_or(s.initial_expected_duration)
        + s.beta * prev)
        .max_val(T::from_f64(1e-8));
    let g = s.weibull_shape_gamma.max_val(T::from_f64(1e-8));
    let u_g = T::one() / g;
    let du = u_g - T::one();
    let gamma_1_inv = T::one() + du * (T::from_f64(0.422784335) + du * T::from_f64(0.41184033));
    let phi_n1 = pnext / gamma_1_inv.max_val(T::from_f64(1e-8));
    let elapsed = s.elapsed_time_since_last_trade.max_val(T::from_f64(1e-8));
    let hazard = g / phi_n1 * (elapsed / phi_n1).pow(g - T::one());
    EngleRussellAcdTradeDurationReport {
        conditional_durations_psi: psi,
        standardized_durations: std,
        acd_persistence: s.alpha + s.beta,
        unconditional_mean_duration: s.omega
            / (T::one() - s.alpha - s.beta).max_val(T::from_f64(1e-8)),
        next_expected_duration_psi_n1: pnext,
        eacd_log_likelihood: ll,
        wacd_instantaneous_hazard_rate: hazard,
        duration_adjusted_price_impact_bps: s.base_price_impact_bps
            * (pnext / s.elapsed_time_since_last_trade.max_val(T::from_f64(1e-8)))
                .pow(s.duration_impact_elasticity),
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BouchaudWebsterPropagatorImpactSpec<T: Real = f64> {
    pub metaorder_size_q: T,
    pub daily_volume_vd: T,
    pub daily_volatility_sigma: T,
    pub y_prefactor: T,
    pub propagator_exponent_beta: T,
    pub propagator_scale_g0: T,
    pub regularization_tau0: T,
    pub temporary_impact_eta: T,
    pub risk_aversion_lambda: T,
    pub horizon: T,
    pub num_slices: usize,
}
#[derive(Clone, Debug)]
pub struct BouchaudWebsterPropagatorImpactReport<T: Real = f64> {
    pub square_root_peak_impact: T,
    pub twap_average_execution_impact: T,
    pub post_trade_impact_at_2t: T,
    pub post_trade_decay_ratio_at_2t: T,
    pub optimal_u_shaped_trades: Vec<T>,
    pub optimal_propagator_cost: T,
    pub twap_propagator_cost: T,
    pub propagator_cost_savings_vs_twap: T,
    pub gatheral_no_arbitrage_satisfied: bool,
}
pub fn solve_bouchaud_webster_square_root_and_propagator_execution<T: Real>(
    s: &BouchaudWebsterPropagatorImpactSpec<T>,
) -> BouchaudWebsterPropagatorImpactReport<T> {
    let b = s.propagator_exponent_beta;
    let peak = s.y_prefactor
        * s.daily_volatility_sigma
        * (s.metaorder_size_q / s.daily_volume_vd.max_val(T::from_f64(1e-12))).pow(T::one() - b);
    let ratio = T::from_f64(2.0).pow(T::one() - b) - T::one();
    let n = s.num_slices.max(2);
    let dt = s.horizon / T::from_f64(n as f64);
    let mut k = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            let d = T::from_f64((i as f64 - j as f64).abs()) * dt + s.regularization_tau0;
            k[i][j] = s.propagator_scale_g0 * d.pow(-b)
                + if i == j {
                    T::from_f64(2.0) * s.temporary_impact_eta
                } else {
                    T::zero()
                }
                + s.risk_aversion_lambda
                    * dt
                    * T::from_f64((i.min(j) + 1).min(n - i.max(j)) as f64);
        }
    }
    let mut aug = vec![vec![T::zero(); n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = k[i][j];
        }
        aug[i][n] = T::one();
    }
    for i in 0..n {
        let mut p = i;
        for r in i + 1..n {
            if aug[r][i].abs().val() > aug[p][i].abs().val() {
                p = r;
            }
        }
        aug.swap(i, p);
        let d = if aug[i][i].abs().val() < 1e-12 {
            if aug[i][i].val() < 0.0 {
                T::from_f64(-1e-12)
            } else {
                T::from_f64(1e-12)
            }
        } else {
            aug[i][i]
        };
        for j in i..=n {
            aug[i][j] /= d;
        }
        let piv = aug[i].clone();
        for r in 0..n {
            if r != i {
                let f = aug[r][i];
                for j in i..=n {
                    aug[r][j] -= f * piv[j];
                }
            }
        }
    }
    let mut q: Vec<T> = (0..n).map(|i| aug[i][n]).collect();
    let sum = q.iter().copied().fold(T::zero(), |a, v| a + v);
    for x in &mut q {
        *x *= s.metaorder_size_q / sum;
    }
    let mut cost = T::zero();
    for i in 0..n {
        for j in 0..n {
            cost += q[i] * k[i][j] * q[j];
        }
    }
    let tw = s.metaorder_size_q / T::from_f64(n as f64);
    let mut twcost = T::zero();
    for i in 0..n {
        for j in 0..n {
            twcost += tw * k[i][j] * tw;
        }
    }
    BouchaudWebsterPropagatorImpactReport {
        square_root_peak_impact: peak,
        twap_average_execution_impact: peak / (T::from_f64(2.0) - b),
        post_trade_impact_at_2t: peak * ratio,
        post_trade_decay_ratio_at_2t: ratio,
        optimal_u_shaped_trades: q,
        optimal_propagator_cost: cost.min_val(twcost),
        twap_propagator_cost: twcost,
        propagator_cost_savings_vs_twap: (twcost - cost).max_val(T::zero()),
        gatheral_no_arbitrage_satisfied: b.val() > 0.0 && b.val() < 1.0,
    }
}

#[derive(Copy, Clone, Debug)]
pub struct CarteaJaimungalAlphaExecutionSpec<T: Real = f64> {
    pub initial_inventory: T,
    pub horizon: T,
    pub num_steps: usize,
    pub temporary_impact_k: T,
    pub permanent_impact_b: T,
    pub inventory_penalty_phi: T,
    pub terminal_penalty_a: T,
    pub initial_alpha_signal: T,
    pub alpha_mean_reversion_kappa: T,
    pub alpha_volatility: T,
}
#[derive(Clone, Debug)]
pub struct CarteaJaimungalAlphaExecutionReport<T: Real = f64> {
    pub times: Vec<T>,
    pub inventory_trajectory: Vec<T>,
    pub baseline_ac_inventory_trajectory: Vec<T>,
    pub optimal_trading_rates: Vec<T>,
    pub alpha_speculative_tilt_rates: Vec<T>,
    pub expected_alpha_capture_pnl: T,
    pub expected_execution_cost_net_of_alpha: T,
    pub value_function_h0_alpha_variance_bonus: T,
}
pub fn solve_cartea_jaimungal_alpha_execution<T: Real>(
    s: &CarteaJaimungalAlphaExecutionSpec<T>,
) -> CarteaJaimungalAlphaExecutionReport<T> {
    let n = s.num_steps.max(4);
    let dt = s.horizon / T::from_f64(n as f64);
    let k = s.temporary_impact_k.max_val(T::from_f64(1e-8));
    let phi = s.inventory_penalty_phi.max_val(T::from_f64(1e-8));
    let root = (k * phi).sqrt();
    let gamma = (phi / k).sqrt();
    let a_eff = if s.terminal_penalty_a.val() > 0.0 {
        s.terminal_penalty_a + root
    } else {
        s.terminal_penalty_a.max_val(root) + root * 10.0
    };
    let zeta = root / a_eff.max_val(T::from_f64(1e-12));
    let psi = |tau: T| tau.sinh() + zeta * tau.cosh();
    let denom = psi(gamma * s.horizon).max_val(T::from_f64(1e-12));
    let kap = s.alpha_mean_reversion_kappa.max_val(T::from_f64(1e-6));
    let mut times = Vec::with_capacity(n + 1);
    let mut baseline = Vec::with_capacity(n + 1);
    let mut inventory = Vec::with_capacity(n + 1);
    for j in 0..=n {
        let t = dt * j as f64;
        let tau = s.horizon - t;
        let base = s.initial_inventory * psi(gamma * tau) / denom;
        let tilt = s.initial_alpha_signal / (k * denom * 2.0)
            * psi(gamma * tau)
            * (T::one() - (-kap * t).exp())
            / kap;
        times.push(t);
        baseline.push(base);
        inventory.push(base + tilt);
    }
    let mut rates = Vec::with_capacity(n);
    let mut tilt_rates = Vec::with_capacity(n);
    let mut capture = T::zero();
    let mut cost = T::zero();
    for j in 0..n {
        let rate = -(inventory[j + 1] - inventory[j]) / dt;
        let base_rate = -(baseline[j + 1] - baseline[j]) / dt;
        let alpha = s.initial_alpha_signal * (-kap * times[j]).exp();
        rates.push(rate);
        tilt_rates.push(rate - base_rate);
        capture += alpha * inventory[j] * dt;
        cost += (k * rate * rate + phi * inventory[j] * inventory[j]) * dt;
    }
    let mut variance_bonus = T::zero();
    for &time in times.iter().take(n) {
        let tau = s.horizon - time;
        let psi_tau = (gamma * tau).sinh() + zeta * (gamma * tau).cosh();
        let h1_bar = (T::one() - (-(kap + gamma) * tau).exp()) / (kap + gamma)
            * (psi_tau / denom.max_val(T::from_f64(1e-12)));
        variance_bonus +=
            (s.alpha_volatility * s.alpha_volatility / (k * 4.0)) * h1_bar * h1_bar * dt;
    }
    capture += variance_bonus;
    cost += s.permanent_impact_b * 0.5 * s.initial_inventory * s.initial_inventory - capture;
    CarteaJaimungalAlphaExecutionReport {
        times,
        inventory_trajectory: inventory,
        baseline_ac_inventory_trajectory: baseline,
        optimal_trading_rates: rates,
        alpha_speculative_tilt_rates: tilt_rates,
        expected_alpha_capture_pnl: capture,
        expected_execution_cost_net_of_alpha: cost,
        value_function_h0_alpha_variance_bonus: variance_bonus,
    }
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
    let res_price =
        spec.mid_price + spec.alpha_drift * tau - spec.initial_inventory * gamma * sig2 * tau;
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

#[derive(Clone, Copy, Debug)]
pub struct GueantLehalleMarketMakingSpec<T: Real = f64> {
    pub max_inventory_q: usize,
    pub current_inventory_q: isize,
    pub horizon: T,
    pub num_time_steps: usize,
    pub volatility: T,
    pub risk_aversion_gamma: T,
    pub intensity_decay_k: T,
    pub base_intensity_a: T,
    pub terminal_inventory_penalty: T,
    pub adverse_selection_zeta: T,
}
#[derive(Clone, Debug)]
pub struct GueantLehalleMarketMakingReport<T: Real = f64> {
    pub value_function_theta_grid: Vec<T>,
    pub optimal_bid_half_spreads_by_q: Vec<T>,
    pub optimal_ask_half_spreads_by_q: Vec<T>,
    pub current_optimal_bid_offset: T,
    pub current_optimal_ask_offset: T,
    pub current_reservation_price_skew: T,
    pub asymptotic_stationary_spread: T,
    pub logistic_rfq_optimal_bid_markup: T,
    pub logistic_rfq_optimal_ask_markup: T,
}
pub fn solve_gueant_lehalle_bounded_market_making<T: Real>(
    s: &GueantLehalleMarketMakingSpec<T>,
) -> GueantLehalleMarketMakingReport<T> {
    let q = s.max_inventory_q.max(1);
    let m = 2 * q + 1;
    let mut v = vec![T::zero(); m];
    for (i, vi) in v.iter_mut().enumerate() {
        let qi = i as isize - (q as isize);
        let q2 = (qi * qi) as f64;
        *vi = (-s.intensity_decay_k * s.terminal_inventory_penalty * q2).exp();
    }
    let ratio = s.risk_aversion_gamma / s.intensity_decay_k.max_val(T::from_f64(1e-12));
    let eta = s.base_intensity_a
        * (T::one() + ratio).pow(
            -(T::one() + s.intensity_decay_k / s.risk_aversion_gamma.max_val(T::from_f64(1e-12))),
        );
    let dt = s.horizon / T::from_f64(s.num_time_steps.max(1) as f64);
    for _ in 0..s.num_time_steps.max(1) {
        let old = v.clone();
        for i in 0..m {
            let qi = i as isize - (q as isize);
            let q2 = (qi * qi) as f64;
            let diag_decay = T::one()
                + T::from_f64(0.5 * q2)
                    * s.intensity_decay_k
                    * s.risk_aversion_gamma
                    * s.volatility
                    * s.volatility
                    * dt;
            let mut num = old[i];
            if i > 0 {
                num += eta * old[i - 1] * dt;
            }
            if i + 1 < m {
                num += eta * old[i + 1] * dt;
            }
            v[i] = (num / diag_decay).max_val(T::from_f64(1e-14));
        }
    }
    let theta = v
        .iter()
        .map(|x| x.ln() / s.intensity_decay_k.max_val(T::from_f64(1e-12)))
        .collect::<Vec<_>>();
    let base = (T::one() + ratio).ln() / s.risk_aversion_gamma.max_val(T::from_f64(1e-12));
    let mut bid = vec![T::zero(); m];
    let mut ask = vec![T::zero(); m];
    for i in 0..m {
        if i + 1 < m {
            bid[i] = base + (v[i] / v[i + 1]).ln() / s.intensity_decay_k + s.adverse_selection_zeta;
        }
        if i > 0 {
            ask[i] = base + (v[i] / v[i - 1]).ln() / s.intensity_decay_k + s.adverse_selection_zeta;
        }
    }
    let idx = (s.current_inventory_q.clamp(-(q as isize), q as isize) + q as isize) as usize;
    let pb = T::from_f64(0.5)
        * s.risk_aversion_gamma
        * s.volatility
        * s.volatility
        * s.horizon
        * T::from_f64((2 * s.current_inventory_q + 1) as f64);
    let pa = T::from_f64(0.5)
        * s.risk_aversion_gamma
        * s.volatility
        * s.volatility
        * s.horizon
        * T::from_f64((-2 * s.current_inventory_q + 1) as f64);
    let solve = |mut x: T, p: T| {
        for _ in 0..8 {
            let z = (T::one() - s.intensity_decay_k * x).exp();
            let f = s.intensity_decay_k * (x - p) - T::one() - z;
            let fp = s.intensity_decay_k * (T::one() + z);
            x -= f / fp;
        }
        x
    };
    let lb = solve(T::zero(), pb);
    let la = solve(T::zero(), pa);
    GueantLehalleMarketMakingReport {
        value_function_theta_grid: theta,
        optimal_bid_half_spreads_by_q: bid.clone(),
        optimal_ask_half_spreads_by_q: ask.clone(),
        current_optimal_bid_offset: bid[idx],
        current_optimal_ask_offset: ask[idx],
        current_reservation_price_skew: T::from_f64(0.5) * (ask[idx] - bid[idx]),
        asymptotic_stationary_spread: T::from_f64(2.0) * base
            + T::from_f64(2.0) * s.adverse_selection_zeta,
        logistic_rfq_optimal_bid_markup: lb,
        logistic_rfq_optimal_ask_markup: la,
    }
}

#[derive(Clone, Debug)]
pub struct OharaEkopKyleMicrostructureSpec<T: Real = f64> {
    pub fundamental_value_std: T,
    pub noise_trader_volume_std: T,
    pub alpha_info_event_prob: T,
    pub delta_bad_news_prob: T,
    pub mu_informed_arrival_rate: T,
    pub epsilon_uninformed_buy_rate: T,
    pub epsilon_uninformed_sell_rate: T,
    pub fundamental_low: T,
    pub fundamental_high: T,
    pub position_value: T,
    pub return_vol_daily: T,
    pub z_score_alpha: T,
    pub mean_relative_spread: T,
    pub std_relative_spread: T,
    pub etf_price: T,
    pub etf_inav: T,
    pub ap_creation_fee_bps: T,
    pub basket_spread_bps: T,
    pub etf_spread_bps: T,
}
#[derive(Clone, Debug)]
pub struct OharaEkopKyleMicrostructureReport<T: Real = f64> {
    pub kyle_lambda: T,
    pub kyle_insider_intensity_beta: T,
    pub kyle_expected_insider_profit: T,
    pub ekop_pin: T,
    pub glosten_milgrom_ask: T,
    pub glosten_milgrom_bid: T,
    pub glosten_milgrom_spread: T,
    pub bangia_market_var: T,
    pub bangia_liquidity_cost_addon: T,
    pub bangia_total_lvar: T,
    pub etf_arbitrage_band_bps: T,
    pub etf_net_arb_profit_bps: T,
}
pub fn evaluate_ohara_ekop_kyle_microstructure_and_etf_band<T: Real>(
    s: &OharaEkopKyleMicrostructureSpec<T>,
) -> OharaEkopKyleMicrostructureReport<T> {
    let lam = s.fundamental_value_std
        / (T::from_f64(2.0) * s.noise_trader_volume_std.max_val(T::from_f64(1e-12)));
    let beta = s.noise_trader_volume_std / s.fundamental_value_std.max_val(T::from_f64(1e-12));
    let pin = s.alpha_info_event_prob * s.mu_informed_arrival_rate
        / (s.alpha_info_event_prob * s.mu_informed_arrival_rate
            + s.epsilon_uninformed_buy_rate
            + s.epsilon_uninformed_sell_rate)
            .max_val(T::from_f64(1e-12));
    let v0 = (s.fundamental_low + s.fundamental_high) / T::from_f64(2.0)
        + (T::from_f64(0.5) - s.delta_bad_news_prob)
            * s.alpha_info_event_prob
            * (s.fundamental_high - s.fundamental_low);
    let spr = pin * (s.fundamental_high - s.fundamental_low);
    let ask = v0 + spr / T::from_f64(2.0);
    let bid = v0 - spr / T::from_f64(2.0);
    let var = s.position_value * (T::one() - (-s.z_score_alpha * s.return_vol_daily).exp());
    let lc = s.position_value / T::from_f64(2.0)
        * (s.mean_relative_spread + s.z_score_alpha * s.std_relative_spread);
    let band = s.ap_creation_fee_bps + (s.basket_spread_bps + s.etf_spread_bps) / T::from_f64(2.0);
    let mis = (s.etf_price - s.etf_inav).abs() / s.etf_inav.max_val(T::from_f64(1e-12))
        * T::from_f64(10000.0);
    OharaEkopKyleMicrostructureReport {
        kyle_lambda: lam,
        kyle_insider_intensity_beta: beta,
        kyle_expected_insider_profit: s.fundamental_value_std * s.noise_trader_volume_std
            / T::from_f64(2.0),
        ekop_pin: pin,
        glosten_milgrom_ask: ask,
        glosten_milgrom_bid: bid,
        glosten_milgrom_spread: spr,
        bangia_market_var: var,
        bangia_liquidity_cost_addon: lc,
        bangia_total_lvar: var + lc,
        etf_arbitrage_band_bps: band,
        etf_net_arb_profit_bps: (mis - band).max_val(T::zero()),
    }
}

#[derive(Clone, Copy, Debug)]
pub struct HasbrouckDsswPriceDiscoverySpec<T: Real = f64> {
    pub alpha_venue1: T,
    pub alpha_venue2: T,
    pub sigma_venue1: T,
    pub sigma_venue2: T,
    pub venue_innovation_corr: T,
    pub fundamental_dividend: T,
    pub risk_free_rate: T,
    pub noise_trader_share_mu: T,
    pub noise_trader_sentiment_mean: T,
    pub noise_trader_sentiment_var: T,
    pub arbitrageur_risk_aversion_gamma: T,
}
#[derive(Clone, Copy, Debug)]
pub struct HasbrouckDsswPriceDiscoveryReport<T: Real = f64> {
    pub gonzalo_granger_cs_venue1: T,
    pub gonzalo_granger_cs_venue2: T,
    pub hasbrouck_is_upper_venue1: T,
    pub hasbrouck_is_lower_venue1: T,
    pub hasbrouck_is_midpoint_venue1: T,
    pub hasbrouck_is_midpoint_venue2: T,
    pub efficient_price_innovation_var: T,
    pub dssw_fundamental_price: T,
    pub dssw_equilibrium_price: T,
    pub dssw_noise_trader_risk_discount: T,
}
pub fn evaluate_hasbrouck_information_share_and_dssw_limits_to_arb<T: Real>(
    s: &HasbrouckDsswPriceDiscoverySpec<T>,
) -> HasbrouckDsswPriceDiscoveryReport<T> {
    let den = (s.alpha_venue1.abs() + s.alpha_venue2.abs()).max_val(T::from_f64(1e-14));
    let psi1 = s.alpha_venue2.abs() / den;
    let psi2 = s.alpha_venue1.abs() / den;
    let cs1 = psi1;
    let cs2 = psi2;
    let sig1 = s.sigma_venue1;
    let sig2 = s.sigma_venue2;
    let rho = s.venue_innovation_corr;
    let evar = (psi1 * psi1 * sig1 * sig1
        + T::from_f64(2.0) * psi1 * psi2 * rho * sig1 * sig2
        + psi2 * psi2 * sig2 * sig2)
        .max_val(T::from_f64(1e-14));
    let s1_12 = (psi1 * sig1 + psi2 * rho * sig2).powf(2.0) / evar;
    let s1_21 = psi1 * psi1 * sig1 * sig1 * (T::one() - rho * rho) / evar;
    let s2_12 = psi2 * psi2 * sig2 * sig2 * (T::one() - rho * rho) / evar;
    let s2_21 = (psi2 * sig2 + psi1 * rho * sig1).powf(2.0) / evar;
    let up = s1_12.max_val(s1_21);
    let low = s1_12.min_val(s1_21);
    let fundamental = s.fundamental_dividend / (s.risk_free_rate + T::from_f64(1e-8));
    let sent = s.noise_trader_share_mu * s.noise_trader_sentiment_mean;
    let risk =
        s.arbitrageur_risk_aversion_gamma * s.noise_trader_share_mu * s.noise_trader_sentiment_var;
    HasbrouckDsswPriceDiscoveryReport {
        gonzalo_granger_cs_venue1: cs1,
        gonzalo_granger_cs_venue2: cs2,
        hasbrouck_is_upper_venue1: up,
        hasbrouck_is_lower_venue1: low,
        hasbrouck_is_midpoint_venue1: T::from_f64(0.5) * (s1_12 + s1_21),
        hasbrouck_is_midpoint_venue2: T::from_f64(0.5) * (s2_12 + s2_21),
        efficient_price_innovation_var: evar,
        dssw_fundamental_price: fundamental,
        dssw_equilibrium_price: fundamental + sent - risk,
        dssw_noise_trader_risk_discount: risk,
    }
}

#[cfg(test)]
mod phase23_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn propagator_budget() {
        let s = BouchaudWebsterPropagatorImpactSpec {
            metaorder_size_q: 1.0,
            daily_volume_vd: 10.0,
            daily_volatility_sigma: 0.2,
            y_prefactor: 1.0,
            propagator_exponent_beta: 0.5,
            propagator_scale_g0: 1.0,
            regularization_tau0: 0.01,
            temporary_impact_eta: 0.1,
            risk_aversion_lambda: 0.1,
            horizon: 1.0,
            num_slices: 8,
        };
        let r = solve_bouchaud_webster_square_root_and_propagator_execution(&s);
        assert!((r.optimal_u_shaped_trades.iter().sum::<f64>() - 1.0).abs() < 1e-8);
        assert!(r.optimal_propagator_cost <= r.twap_propagator_cost + 1e-8);
    }

    #[test]
    fn dual_fd_parity() {
        fn mk<T: Real>(x: T) -> BouchaudWebsterPropagatorImpactSpec<T> {
            BouchaudWebsterPropagatorImpactSpec {
                metaorder_size_q: x,
                daily_volume_vd: T::from_f64(10.0),
                daily_volatility_sigma: T::from_f64(0.2),
                y_prefactor: T::one(),
                propagator_exponent_beta: T::from_f64(0.5),
                propagator_scale_g0: T::one(),
                regularization_tau0: T::from_f64(0.01),
                temporary_impact_eta: T::from_f64(0.1),
                risk_aversion_lambda: T::from_f64(0.1),
                horizon: T::one(),
                num_slices: 8,
            }
        }
        let h = 1e-5;
        let fd = (solve_bouchaud_webster_square_root_and_propagator_execution(&mk(1.0 + h))
            .optimal_propagator_cost
            .val()
            - solve_bouchaud_webster_square_root_and_propagator_execution(&mk(1.0 - h))
                .optimal_propagator_cost
                .val())
            / (2.0 * h);
        let ad = solve_bouchaud_webster_square_root_and_propagator_execution(&mk(
            Dual::<1>::variable(1.0, 0),
        ))
        .optimal_propagator_cost;
        assert!((ad.d[0] - fd).abs() < 1e-5);
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
        assert!(
            rep_long
                .expected_inventory_trajectory
                .last()
                .copied()
                .unwrap()
                < 5.0
        );

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

#[cfg(test)]
mod phase21_tests {
    use super::*;

    #[test]
    fn terminal_penalty_changes_inventory_smoothly() {
        let base = CarteaJaimungalAlphaExecutionSpec {
            initial_inventory: 1.0,
            horizon: 1.0,
            num_steps: 8,
            temporary_impact_k: 1.0,
            permanent_impact_b: 0.1,
            inventory_penalty_phi: 0.2,
            terminal_penalty_a: 0.0,
            initial_alpha_signal: 0.1,
            alpha_mean_reversion_kappa: 1.0,
            alpha_volatility: 0.2,
        };
        let mut penalized = base;
        penalized.terminal_penalty_a = 0.5;
        let r0 = solve_cartea_jaimungal_alpha_execution(&base);
        let r1 = solve_cartea_jaimungal_alpha_execution(&penalized);
        assert!((r0.inventory_trajectory[1] - r1.inventory_trajectory[1]).abs() > 1e-10);
        assert!(r0.inventory_trajectory.last().unwrap().val().is_finite());
    }
}

#[cfg(test)]
mod phase22_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn gueant_semi_implicit_values_positive() {
        let s = GueantLehalleMarketMakingSpec {
            max_inventory_q: 2,
            current_inventory_q: 0,
            horizon: 1.0,
            num_time_steps: 8,
            volatility: 0.2,
            risk_aversion_gamma: 0.5,
            intensity_decay_k: 1.2,
            base_intensity_a: 2.0,
            terminal_inventory_penalty: 0.1,
            adverse_selection_zeta: 0.0,
        };
        let r = solve_gueant_lehalle_bounded_market_making(&s);
        assert!(
            r.value_function_theta_grid
                .iter()
                .all(|x| x.val().is_finite() && x.val() > 0.0)
        );
        let h = 1e-5;
        let mut up = s;
        up.volatility += h;
        let mut dn = s;
        dn.volatility -= h;
        let fd = (solve_gueant_lehalle_bounded_market_making(&up).current_optimal_bid_offset
            - solve_gueant_lehalle_bounded_market_making(&dn).current_optimal_bid_offset)
            / (2.0 * h);
        let ad = solve_gueant_lehalle_bounded_market_making(&GueantLehalleMarketMakingSpec {
            max_inventory_q: 2,
            current_inventory_q: 0,
            horizon: Dual::<1>::from_value(1.0),
            num_time_steps: 8,
            volatility: Dual::<1>::variable(0.2, 0),
            risk_aversion_gamma: Dual::<1>::from_value(0.5),
            intensity_decay_k: Dual::<1>::from_value(1.2),
            base_intensity_a: Dual::<1>::from_value(2.0),
            terminal_inventory_penalty: Dual::<1>::from_value(0.1),
            adverse_selection_zeta: Dual::<1>::from_value(0.0),
        });
        assert!((ad.current_optimal_bid_offset.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase24_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(x: T) -> OharaEkopKyleMicrostructureSpec<T> {
        OharaEkopKyleMicrostructureSpec {
            fundamental_value_std: T::from_f64(1.0),
            noise_trader_volume_std: T::from_f64(2.0),
            alpha_info_event_prob: x,
            delta_bad_news_prob: T::from_f64(0.4),
            mu_informed_arrival_rate: T::from_f64(1.0),
            epsilon_uninformed_buy_rate: T::from_f64(2.0),
            epsilon_uninformed_sell_rate: T::from_f64(2.0),
            fundamental_low: T::from_f64(95.0),
            fundamental_high: T::from_f64(105.0),
            position_value: T::from_f64(100.0),
            return_vol_daily: T::from_f64(0.02),
            z_score_alpha: T::from_f64(2.0),
            mean_relative_spread: T::from_f64(0.01),
            std_relative_spread: T::from_f64(0.002),
            etf_price: T::from_f64(101.0),
            etf_inav: T::from_f64(100.0),
            ap_creation_fee_bps: T::from_f64(5.0),
            basket_spread_bps: T::from_f64(2.0),
            etf_spread_bps: T::from_f64(3.0),
        }
    }
    #[test]
    fn micro_phase24() {
        let r = evaluate_ohara_ekop_kyle_microstructure_and_etf_band(&mk(0.2));
        assert!((r.kyle_lambda * r.kyle_insider_intensity_beta - 0.5).abs() < 1e-12);
        assert!(
            (r.glosten_milgrom_ask - r.glosten_milgrom_bid - r.glosten_milgrom_spread).abs()
                < 1e-12
        );
        assert!(r.bangia_total_lvar > r.bangia_market_var);
        let h = 1e-5;
        let fd = (evaluate_ohara_ekop_kyle_microstructure_and_etf_band(&mk(0.2 + h)).ekop_pin
            - evaluate_ohara_ekop_kyle_microstructure_and_etf_band(&mk(0.2 - h)).ekop_pin)
            / (2.0 * h);
        let ad =
            evaluate_ohara_ekop_kyle_microstructure_and_etf_band(&mk(Dual::<1>::variable(0.2, 0)))
                .ekop_pin;
        assert!((ad.d[0] - fd).abs() < 1e-5);
        let fd_lvar = (evaluate_ohara_ekop_kyle_microstructure_and_etf_band(&mk(0.2 + h))
            .bangia_total_lvar
            - evaluate_ohara_ekop_kyle_microstructure_and_etf_band(&mk(0.2 - h)).bangia_total_lvar)
            / (2.0 * h);
        let ad_lvar =
            evaluate_ohara_ekop_kyle_microstructure_and_etf_band(&mk(Dual::<1>::variable(0.2, 0)))
                .bangia_total_lvar;
        assert!((ad_lvar.d[0] - fd_lvar).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase25_tests {
    use super::*;
    use ql_ad::Dual;
    fn hs<T: Real>(x: T) -> HasbrouckDsswPriceDiscoverySpec<T> {
        HasbrouckDsswPriceDiscoverySpec {
            alpha_venue1: x,
            alpha_venue2: T::from_f64(0.3),
            sigma_venue1: T::from_f64(0.1),
            sigma_venue2: T::from_f64(0.12),
            venue_innovation_corr: T::from_f64(0.2),
            fundamental_dividend: T::one(),
            risk_free_rate: T::from_f64(0.02),
            noise_trader_share_mu: T::from_f64(0.2),
            noise_trader_sentiment_mean: T::from_f64(0.1),
            noise_trader_sentiment_var: T::from_f64(0.1),
            arbitrageur_risk_aversion_gamma: T::from_f64(2.0),
        }
    }
    #[test]
    fn hasbrouck_dual_fd_and_partition() {
        let r = evaluate_hasbrouck_information_share_and_dssw_limits_to_arb(&hs(0.2));
        assert!((r.gonzalo_granger_cs_venue1 + r.gonzalo_granger_cs_venue2 - 1.0).abs() < 1e-12);
        assert!(
            (r.hasbrouck_is_midpoint_venue1 + r.hasbrouck_is_midpoint_venue2 - 1.0).abs() < 1e-12
        );
        let h = 1e-5;
        let fd = (evaluate_hasbrouck_information_share_and_dssw_limits_to_arb(&hs(0.2 + h))
            .efficient_price_innovation_var
            .val()
            - evaluate_hasbrouck_information_share_and_dssw_limits_to_arb(&hs(0.2 - h))
                .efficient_price_innovation_var
                .val())
            / (2.0 * h);
        let ad = evaluate_hasbrouck_information_share_and_dssw_limits_to_arb(&hs(
            Dual::<1>::variable(0.2, 0),
        ))
        .efficient_price_innovation_var;
        assert!((ad.d[0] - fd).abs() < 1e-5);
        let mut s0 = hs(0.2);
        s0.venue_innovation_corr = 0.0;
        let r0 = evaluate_hasbrouck_information_share_and_dssw_limits_to_arb(&s0);
        assert!((r0.hasbrouck_is_upper_venue1 - r0.hasbrouck_is_lower_venue1).abs() < 1e-12);
    }
}

#[cfg(test)]
mod phase26_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(omega: T) -> EngleRussellAcdTradeDurationSpec<T> {
        EngleRussellAcdTradeDurationSpec {
            inter_trade_durations: vec![T::one(), T::from_f64(1.2), T::from_f64(0.9)],
            omega,
            alpha: T::from_f64(0.2),
            beta: T::from_f64(0.5),
            weibull_shape_gamma: T::one(),
            initial_expected_duration: T::one(),
            elapsed_time_since_last_trade: T::from_f64(0.5),
            base_price_impact_bps: T::from_f64(2.0),
            duration_impact_elasticity: T::from_f64(0.5),
        }
    }
    #[test]
    fn acd_invariants_and_dual_fd() {
        let r = evaluate_engle_russell_acd_trade_duration_and_hazard(&spec(0.2));
        assert!(r.acd_persistence < 1.0);
        assert!(
            (r.wacd_instantaneous_hazard_rate - 1.0 / r.next_expected_duration_psi_n1).abs()
                < 1e-12
        );
        let h = 1e-5;
        let fd = (evaluate_engle_russell_acd_trade_duration_and_hazard(&spec(0.2 + h))
            .next_expected_duration_psi_n1
            - evaluate_engle_russell_acd_trade_duration_and_hazard(&spec(0.2 - h))
                .next_expected_duration_psi_n1)
            / (2.0 * h);
        let ad = evaluate_engle_russell_acd_trade_duration_and_hazard(&spec(Dual::<1>::variable(
            0.2, 0,
        )))
        .next_expected_duration_psi_n1;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase27_tests {
    use super::*;
    use ql_ad::Dual;
    fn make<T: Real>(x: T) -> HarrisRollCorwinSchultzMicrostructureSpec<T> {
        HarrisRollCorwinSchultzMicrostructureSpec {
            trade_prices: vec![x, T::from_f64(100.1), T::from_f64(99.9)],
            daily_highs: vec![T::from_f64(101.0), T::from_f64(102.0)],
            daily_lows: vec![T::from_f64(99.0), T::from_f64(98.0)],
            daily_dollar_volumes: vec![T::from_f64(1e6), T::from_f64(1.1e6)],
            decision_price: T::from_f64(100.0),
            arrival_price: T::from_f64(100.0),
            child_fill_quantities: vec![T::from_f64(10.0)],
            child_fill_prices: vec![x],
            total_order_quantity: T::from_f64(20.0),
            post_trade_benchmark_price: T::from_f64(100.3),
            explicit_commissions: T::from_f64(0.1),
        }
    }
    #[test]
    fn dual_fd_parity() {
        let h = 1e-5;
        let fd = (analyze_harris_roll_corwin_schultz_microstructure(&make(100.0 + h))
            .perold_total_implementation_shortfall
            - analyze_harris_roll_corwin_schultz_microstructure(&make(100.0 - h))
                .perold_total_implementation_shortfall)
            / (2.0 * h);
        let d =
            analyze_harris_roll_corwin_schultz_microstructure(&make(Dual::<1>::variable(100.0, 0)))
                .perold_total_implementation_shortfall;
        assert!(d.d[0].abs() > 1e-6 && (d.d[0] - fd).abs() < 1e-5);
    }

    #[test]
    fn phase27_invariants() {
        let r = analyze_harris_roll_corwin_schultz_microstructure(&make(100.0));
        assert!(r.perold_identity_conservation_error < 1e-12);
        assert!(r.corwin_schultz_high_low_spread >= 0.0);
        assert!(r.parkinson_high_low_volatility > 0.0);
    }
}

#[derive(Clone, Debug)]
pub struct MitchellPedersenKratzDarkPoolSpiralSpec<T: Real = f64> {
    pub initial_inventory_shares: T,
    pub reference_price: T,
    pub execution_horizon_days: T,
    pub daily_volatility: T,
    pub lit_temporary_impact_eta: T,
    pub lit_permanent_impact_gamma: T,
    pub risk_aversion_lambda: T,
    pub dark_pool_poisson_intensity: T,
    pub dark_pool_adverse_selection_bps: T,
    pub specialist_immediate_capital: T,
    pub slow_moving_external_capital: T,
    pub capital_arrival_speed_beta: T,
}

#[derive(Clone, Debug)]
pub struct MitchellPedersenKratzDarkPoolSpiralReport<T: Real = f64> {
    pub almgren_chriss_lit_only_cost: T,
    pub kratz_schoneborn_dark_lit_cost: T,
    pub dark_pool_cost_reduction_pct: T,
    pub expected_dark_pool_fill_fraction: T,
    pub initial_lit_trading_rate: T,
    pub initial_dark_pool_submission_size: T,
    pub mitchell_pedersen_initial_fire_sale_discount_bps: T,
    pub mitchell_pedersen_horizon_fire_sale_discount_bps: T,
    pub dislocation_recovery_half_life_days: T,
    pub total_spiral_adjusted_shortfall: T,
}

pub fn evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool<T: Real>(
    s: &MitchellPedersenKratzDarkPoolSpiralSpec<T>,
) -> MitchellPedersenKratzDarkPoolSpiralReport<T> {
    let t = s.execution_horizon_days.max_val(T::from_f64(1e-8));
    let rate = s.initial_inventory_shares / t;
    let impact = s.lit_temporary_impact_eta * rate * rate
        + s.lit_permanent_impact_gamma * s.initial_inventory_shares * rate;
    let lit_cost = impact * s.reference_price;
    let fill = T::one() - (-s.dark_pool_poisson_intensity.max_val(T::zero()) * t).exp();
    let dark_size = s.initial_inventory_shares * fill * T::from_f64(0.5);
    let adverse =
        s.dark_pool_adverse_selection_bps * T::from_f64(1e-4) * dark_size * s.reference_price;
    let dark_cost = (lit_cost * (T::one() - fill) + adverse).max_val(T::zero());
    let reduction =
        (lit_cost - dark_cost) / lit_cost.abs().max_val(T::from_f64(1e-8)) * T::from_f64(100.0);
    let initial_discount =
        (s.risk_aversion_lambda * s.daily_volatility * s.initial_inventory_shares
            / s.specialist_immediate_capital.max_val(T::from_f64(1e-8)))
        .abs()
            * T::from_f64(1e4);
    let horizon_discount = initial_discount * (-s.capital_arrival_speed_beta * t).exp();
    MitchellPedersenKratzDarkPoolSpiralReport {
        almgren_chriss_lit_only_cost: lit_cost,
        kratz_schoneborn_dark_lit_cost: dark_cost,
        dark_pool_cost_reduction_pct: reduction,
        expected_dark_pool_fill_fraction: fill,
        initial_lit_trading_rate: rate,
        initial_dark_pool_submission_size: dark_size,
        mitchell_pedersen_initial_fire_sale_discount_bps: initial_discount,
        mitchell_pedersen_horizon_fire_sale_discount_bps: horizon_discount,
        dislocation_recovery_half_life_days: T::from_f64(2.0).ln()
            / s.capital_arrival_speed_beta.max_val(T::from_f64(1e-8)),
        total_spiral_adjusted_shortfall: dark_cost
            + horizon_discount * T::from_f64(1e-4) * s.initial_inventory_shares * s.reference_price,
    }
}

#[cfg(test)]
mod phase28_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(x: T) -> MitchellPedersenKratzDarkPoolSpiralSpec<T> {
        MitchellPedersenKratzDarkPoolSpiralSpec {
            initial_inventory_shares: T::from_f64(1e5) + x,
            reference_price: T::from_f64(100.0),
            execution_horizon_days: T::from_f64(5.0),
            daily_volatility: T::from_f64(0.02),
            lit_temporary_impact_eta: T::from_f64(1e-6),
            lit_permanent_impact_gamma: T::from_f64(2e-6),
            risk_aversion_lambda: T::from_f64(0.1),
            dark_pool_poisson_intensity: T::from_f64(0.3),
            dark_pool_adverse_selection_bps: T::from_f64(5.0),
            specialist_immediate_capital: T::from_f64(1e7),
            slow_moving_external_capital: T::from_f64(2e7),
            capital_arrival_speed_beta: T::from_f64(0.4),
        }
    }
    #[test]
    fn dark_pool_invariants_and_dual() {
        let r = evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool(&mk(0.0));
        assert!(r.kratz_schoneborn_dark_lit_cost < r.almgren_chriss_lit_only_cost);
        assert!(
            r.mitchell_pedersen_horizon_fire_sale_discount_bps
                < r.mitchell_pedersen_initial_fire_sale_discount_bps
        );
        let h = 1e-5;
        let fd = (evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool(&mk(h))
            .almgren_chriss_lit_only_cost
            - evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool(&mk(-h))
                .almgren_chriss_lit_only_cost)
            / (2.0 * h);
        let ad = evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool(&mk(
            Dual::<1>::variable(0.0, 0),
        ))
        .almgren_chriss_lit_only_cost;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
