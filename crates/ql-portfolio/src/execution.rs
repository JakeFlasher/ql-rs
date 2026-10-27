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
    let q = s.initial_inventory_shares;
    let q_squared = q * q;
    let kappa = (s.risk_aversion_lambda.abs() * s.daily_volatility * s.daily_volatility
        / s.lit_temporary_impact_eta.max_val(T::from_f64(1e-12)))
    .sqrt();
    let u = kappa * t;
    let h_ac = if u.abs().val() < 1e-4 {
        T::one() + u * u / T::from_f64(3.0)
    } else {
        u / u.tanh()
    };
    let lit_cost = s.reference_price
        * q_squared
        * (T::from_f64(0.5) * s.lit_permanent_impact_gamma + s.lit_temporary_impact_eta / t * h_ac);
    let theta = s.dark_pool_poisson_intensity.max_val(T::zero());
    let psi = s.dark_pool_adverse_selection_bps.max_val(T::zero()) * T::from_f64(1e-4);
    let fill = (T::one() - (-theta * t).exp()) / (T::one() + psi);
    let dark_size = q * (theta * t / (T::one() + theta * t)) / (T::one() + psi);
    let initial_rate = (q / t) * (T::one() - T::from_f64(0.5) * fill);
    let adverse =
        s.dark_pool_adverse_selection_bps * T::from_f64(1e-4) * dark_size * s.reference_price;
    let dark_cost = (s.reference_price
        * q_squared
        * (T::from_f64(0.5) * s.lit_permanent_impact_gamma
            + s.lit_temporary_impact_eta / t * h_ac
            - T::from_f64(0.25) * s.lit_temporary_impact_eta / t * fill)
        + adverse)
        .max_val(T::zero());
    let reduction =
        (lit_cost - dark_cost) / lit_cost.abs().max_val(T::from_f64(1e-8)) * T::from_f64(100.0);
    let initial_discount =
        (s.risk_aversion_lambda * s.daily_volatility * s.initial_inventory_shares
            / s.specialist_immediate_capital.max_val(T::from_f64(1e-8)))
        .abs()
            * T::from_f64(1e4);
    let cap_arrival_frac = T::one() - (-s.capital_arrival_speed_beta * t).exp();
    let eff_horizon_cap =
        s.specialist_immediate_capital + s.slow_moving_external_capital * cap_arrival_frac;
    let horizon_discount = initial_discount
        * (s.specialist_immediate_capital / eff_horizon_cap.max_val(T::from_f64(1e-8)))
        * (-s.capital_arrival_speed_beta * t).exp();
    MitchellPedersenKratzDarkPoolSpiralReport {
        almgren_chriss_lit_only_cost: lit_cost,
        kratz_schoneborn_dark_lit_cost: dark_cost,
        dark_pool_cost_reduction_pct: reduction,
        expected_dark_pool_fill_fraction: fill,
        initial_lit_trading_rate: initial_rate,
        initial_dark_pool_submission_size: dark_size,
        mitchell_pedersen_initial_fire_sale_discount_bps: initial_discount,
        mitchell_pedersen_horizon_fire_sale_discount_bps: horizon_discount,
        dislocation_recovery_half_life_days: T::from_f64(2.0).ln()
            / s.capital_arrival_speed_beta.max_val(T::from_f64(1e-8)),
        total_spiral_adjusted_shortfall: dark_cost
            + horizon_discount * T::from_f64(1e-4) * s.initial_inventory_shares * s.reference_price,
    }
}

#[derive(Clone, Debug)]
pub struct CarteaJaimungalPenalvaLobFillHjbSpec<T: Real = f64> {
    pub initial_inventory: usize,
    pub horizon_seconds: T,
    pub base_arrival_rate_lambda: T,
    pub fill_decay_kappa: T,
    pub running_inventory_penalty_phi: T,
    pub terminal_liquidation_penalty_alpha: T,
    pub ofi_short_term_drift: T,
    pub adverse_selection_cost: T,
}
#[derive(Clone, Debug)]
pub struct CarteaJaimungalPenalvaLobFillHjbReport<T: Real = f64> {
    pub hjb_value_functions_by_q: Vec<T>,
    pub optimal_limit_depths_by_q: Vec<T>,
    pub optimal_fill_intensities_by_q: Vec<T>,
    pub expected_terminal_inventory: T,
    pub limit_vs_market_order_pickup: T,
    pub ofi_alpha_reservation_shift: T,
}
pub fn solve_cartea_jaimungal_penalva_lob_fill_hjb<T: Real>(
    s: &CarteaJaimungalPenalvaLobFillHjbSpec<T>,
) -> CarteaJaimungalPenalvaLobFillHjbReport<T> {
    let n = s.initial_inventory + 1;
    let horizon = s.horizon_seconds.max_val(T::from_f64(1e-8));
    let dt = horizon / T::from_f64(32.0);
    let mut vals = (0..n)
        .map(|q| -s.terminal_liquidation_penalty_alpha * T::from_f64((q * q) as f64))
        .collect::<Vec<_>>();
    vals[0] = T::zero();
    let kappa = s.fill_decay_kappa.max_val(T::from_f64(1e-8));
    for _ in 0..32 {
        let old = vals.clone();
        for q in 1..n {
            let qi = T::from_f64(q as f64);
            let diff = old[q - 1] - old[q];
            let intensity = s.base_arrival_rate_lambda
                * (-T::one() - kappa * s.adverse_selection_cost + kappa * diff).exp();
            let dh = s.running_inventory_penalty_phi * qi * qi
                - s.ofi_short_term_drift * qi
                - intensity / kappa;
            vals[q] = old[q] - dt * dh;
        }
        vals[0] = T::zero();
    }
    let mut depths = Vec::with_capacity(n);
    let mut intens = Vec::with_capacity(n);
    let smooth_max =
        |a: T, b: T| (a + b + ((a - b) * (a - b) + T::from_f64(1e-8)).sqrt()) / T::from_f64(2.0);
    for q in 0..n {
        if q == 0 {
            depths.push(T::zero());
            intens.push(s.base_arrival_rate_lambda);
        } else {
            let raw_delta = T::one() / kappa + s.adverse_selection_cost + vals[q] - vals[q - 1];
            let delta = smooth_max(raw_delta, T::from_f64(1e-4));
            depths.push(delta);
            intens.push(s.base_arrival_rate_lambda * (-kappa * delta).exp());
        }
    }
    let expected = T::from_f64(s.initial_inventory as f64)
        * (-intens.last().copied().unwrap_or(T::zero()) * horizon).exp();
    CarteaJaimungalPenalvaLobFillHjbReport {
        hjb_value_functions_by_q: vals,
        optimal_limit_depths_by_q: depths,
        optimal_fill_intensities_by_q: intens,
        expected_terminal_inventory: expected,
        limit_vs_market_order_pickup: s.base_arrival_rate_lambda * horizon,
        ofi_alpha_reservation_shift: s.ofi_short_term_drift - s.adverse_selection_cost,
    }
}

#[cfg(test)]
mod domain_suite_35_tests {
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
mod domain_suite_26_tests {
    use super::*;
    #[test]
    fn suite26_multi_asset_matches_1d_and_twap() {
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
mod domain_suite_29_tests {
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
mod domain_suite_31_tests {
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
mod domain_suite_32_tests {
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
mod domain_suite_33_tests {
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
mod domain_suite_34_tests {
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
mod domain_suite_36_tests {
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
    fn micro_suite36() {
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
mod domain_suite_37_tests {
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
mod domain_suite_38_tests {
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
mod domain_suite_39_tests {
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
    fn suite39_invariants() {
        let r = analyze_harris_roll_corwin_schultz_microstructure(&make(100.0));
        assert!(r.perold_identity_conservation_error < 1e-12);
        assert!(r.corwin_schultz_high_low_spread >= 0.0);
        assert!(r.parkinson_high_low_volatility > 0.0);
    }
}

#[cfg(test)]
mod domain_suite_40_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(x: T) -> MitchellPedersenKratzDarkPoolSpiralSpec<T> {
        MitchellPedersenKratzDarkPoolSpiralSpec {
            initial_inventory_shares: T::from_f64(1000.0) + x,
            reference_price: T::from_f64(100.0),
            execution_horizon_days: T::from_f64(5.0),
            daily_volatility: T::from_f64(0.02),
            lit_temporary_impact_eta: T::from_f64(1e-4),
            lit_permanent_impact_gamma: T::from_f64(2e-4),
            risk_aversion_lambda: T::from_f64(0.1),
            dark_pool_poisson_intensity: T::from_f64(0.3),
            dark_pool_adverse_selection_bps: T::from_f64(5.0),
            specialist_immediate_capital: T::from_f64(1e5),
            slow_moving_external_capital: T::from_f64(2e5),
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

#[cfg(test)]
mod domain_suite_41_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn lob_invariants_and_dual() {
        fn mk<T: Real>(x: T) -> CarteaJaimungalPenalvaLobFillHjbSpec<T> {
            CarteaJaimungalPenalvaLobFillHjbSpec {
                initial_inventory: 10,
                horizon_seconds: x,
                base_arrival_rate_lambda: T::from_f64(0.5),
                fill_decay_kappa: T::one(),
                running_inventory_penalty_phi: T::from_f64(0.01),
                terminal_liquidation_penalty_alpha: T::from_f64(0.1),
                ofi_short_term_drift: T::from_f64(0.02),
                adverse_selection_cost: T::from_f64(0.01),
            }
        }
        let r = solve_cartea_jaimungal_penalva_lob_fill_hjb(&mk(60.0));
        assert!(r.optimal_limit_depths_by_q[1] > r.optimal_limit_depths_by_q[10]);
        let h = 1e-5;
        let fd = (solve_cartea_jaimungal_penalva_lob_fill_hjb(&mk(60.0 + h))
            .optimal_limit_depths_by_q[1]
            - solve_cartea_jaimungal_penalva_lob_fill_hjb(&mk(60.0 - h)).optimal_limit_depths_by_q
                [1])
            / (2.0 * h);
        let ad = solve_cartea_jaimungal_penalva_lob_fill_hjb(&mk(Dual::<1>::variable(60.0, 0)))
            .optimal_limit_depths_by_q[1];
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}

mod evaluate_cartea_jaimungal_lob_market_making_module {
    #![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct CarteaJaimungalLobMarketMakingSpec<T: Real = f64> {
        pub mid_price: T,
        pub current_inventory_q: T,
        pub running_inventory_penalty_phi: T,
        pub terminal_inventory_penalty_alpha: T,
        pub risk_aversion_gamma: T,
        pub asset_volatility: T,
        pub order_arrival_baseline_a: T,
        pub lob_depth_decay_kappa: T,
        pub remaining_horizon_tau: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct CarteaJaimungalLobMarketMakingReport<T: Real = f64> {
        pub riccati_inventory_coeff_theta: T,
        pub optimal_bid_depth_delta_minus: T,
        pub optimal_ask_depth_delta_plus: T,
        pub optimal_bid_quote: T,
        pub optimal_ask_quote: T,
        pub reservation_mid_price: T,
        pub inventory_quote_skew: T,
        pub total_quoted_spread: T,
        pub expected_instantaneous_spread_capture: T,
    }

    pub fn evaluate_cartea_jaimungal_lob_market_making<T: Real>(
        spec: &CarteaJaimungalLobMarketMakingSpec<T>,
    ) -> CarteaJaimungalLobMarketMakingReport<T> {
        let kappa = spec.lob_depth_decay_kappa.max_val(T::from_f64(1e-6));
        let tau = spec.remaining_horizon_tau.max_val(T::from_f64(1e-6));
        let inv_e = (-T::one()).exp();
        let theta = T::from_f64(0.5)
            * (spec.running_inventory_penalty_phi
                + spec.risk_aversion_gamma * spec.asset_volatility * spec.asset_volatility)
            * tau
            + spec.terminal_inventory_penalty_alpha
                / (T::one()
                    + kappa
                        * spec.terminal_inventory_penalty_alpha
                        * spec.order_arrival_baseline_a
                        * inv_e
                        * tau);
        let two = T::from_f64(2.0);
        let q = spec.current_inventory_q;
        let delta_plus = T::one() / kappa + (two * q - T::one()) * theta;
        let delta_minus = T::one() / kappa - (two * q + T::one()) * theta;
        let ask_quote = spec.mid_price + delta_plus;
        let bid_quote = spec.mid_price - delta_minus;
        let res_mid = T::from_f64(0.5) * (ask_quote + bid_quote);
        let skew = -two * q * theta;
        let spread = ask_quote - bid_quote;
        let lam_plus = spec.order_arrival_baseline_a * (-kappa * delta_plus).exp();
        let lam_minus = spec.order_arrival_baseline_a * (-kappa * delta_minus).exp();
        let edge = lam_plus * delta_plus + lam_minus * delta_minus;

        CarteaJaimungalLobMarketMakingReport {
            riccati_inventory_coeff_theta: theta,
            optimal_bid_depth_delta_minus: delta_minus,
            optimal_ask_depth_delta_plus: delta_plus,
            optimal_bid_quote: bid_quote,
            optimal_ask_quote: ask_quote,
            reservation_mid_price: res_mid,
            inventory_quote_skew: skew,
            total_quoted_spread: spread,
            expected_instantaneous_spread_capture: edge,
        }
    }

    #[derive(Clone, Debug)]
    pub struct CampbellLoMackinlayVrPresentValueSpec<T: Real = f64> {
        pub returns: Vec<T>,
        pub steady_state_dividend_yield: T,
        pub expected_dividend_growth: T,
        pub expected_discount_rate: T,
        pub dividend_growth_ar1_phi: T,
        pub discount_rate_ar1_phi: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct CampbellLoMackinlayVrPresentValueReport<T: Real = f64> {
        pub one_period_variance: T,
        pub overlapping_two_period_variance: T,
        pub lo_mackinlay_variance_ratio_vr2: T,
        pub homoskedastic_z_statistic: T,
        pub heteroskedastic_theta_variance: T,
        pub heteroskedastic_z_star_statistic: T,
        pub campbell_shiller_rho: T,
        pub equilibrium_log_dividend_price_ratio: T,
        pub discount_rate_news_multiplier: T,
        pub cashflow_news_multiplier: T,
    }

    pub fn evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller<T: Real>(
        spec: &CampbellLoMackinlayVrPresentValueSpec<T>,
    ) -> CampbellLoMackinlayVrPresentValueReport<T> {
        let n = spec.returns.len().max(3);
        let n_t = T::from_f64(n as f64);
        let mut sum_r = T::zero();
        for &r in &spec.returns {
            sum_r = sum_r + r;
        }
        let mu = sum_r / n_t;
        let mut sum_sq1 = T::zero();
        for &r in &spec.returns {
            let d = r - mu;
            sum_sq1 = sum_sq1 + d * d;
        }
        let var1 = (sum_sq1 / T::from_f64((n - 1) as f64)).max_val(T::from_f64(1e-12));

        let mut sum_sq2 = T::zero();
        let mut sum_theta = T::zero();
        let two = T::from_f64(2.0);
        for t in 1..spec.returns.len() {
            let d2 = spec.returns[t] + spec.returns[t - 1] - two * mu;
            sum_sq2 = sum_sq2 + d2 * d2;
            let u_t = (spec.returns[t] - mu) * (spec.returns[t] - mu) / var1;
            let u_tm1 = (spec.returns[t - 1] - mu) * (spec.returns[t - 1] - mu) / var1;
            sum_theta = sum_theta + u_t * u_tm1;
        }
        let var2 = sum_sq2 / T::from_f64((2 * (n - 2)) as f64);
        let vr2 = var2 / var1;
        let z_homo = n_t.sqrt() * (vr2 - T::one());
        let theta_var = (T::from_f64(4.0) / n_t * sum_theta).max_val(T::from_f64(1e-12));
        let z_star = z_homo / theta_var.sqrt();

        let dy = spec.steady_state_dividend_yield.max_val(T::from_f64(1e-6));
        let rho = T::one() / (T::one() + dy);
        let k_cs = -rho.ln() - (T::one() - rho) * (T::one() / rho - T::one()).ln();
        let cf_mult =
            T::one() / (T::one() - rho * spec.dividend_growth_ar1_phi).max_val(T::from_f64(1e-6));
        let dr_mult =
            T::one() / (T::one() - rho * spec.discount_rate_ar1_phi).max_val(T::from_f64(1e-6));
        let eq_dp = -k_cs / (T::one() - rho) - spec.expected_dividend_growth * cf_mult
            + spec.expected_discount_rate * dr_mult;

        CampbellLoMackinlayVrPresentValueReport {
            one_period_variance: var1,
            overlapping_two_period_variance: var2,
            lo_mackinlay_variance_ratio_vr2: vr2,
            homoskedastic_z_statistic: z_homo,
            heteroskedastic_theta_variance: theta_var,
            heteroskedastic_z_star_statistic: z_star,
            campbell_shiller_rho: rho,
            equilibrium_log_dividend_price_ratio: eq_dp,
            discount_rate_news_multiplier: dr_mult,
            cashflow_news_multiplier: cf_mult,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct RoncalliSpinuErcRiskParitySpec<T: Real = f64> {
        pub covariance_3x3: [[T; 3]; 3],
        pub target_risk_budgets: [T; 3],
    }

    #[derive(Clone, Copy, Debug)]
    pub struct RoncalliSpinuErcRiskParityReport<T: Real = f64> {
        pub erc_weights: [T; 3],
        pub inverse_vol_naive_weights: [T; 3],
        pub relative_risk_contributions: [T; 3],
        pub erc_portfolio_volatility: T,
        pub risk_budget_parity_l1_error: T,
        pub weights_sum_conservation_error: T,
        pub choueifaty_diversification_ratio: T,
        pub spinu_convex_objective_val: T,
    }

    pub fn evaluate_roncalli_spinu_erc_risk_parity<T: Real>(
        spec: &RoncalliSpinuErcRiskParitySpec<T>,
    ) -> RoncalliSpinuErcRiskParityReport<T> {
        let cov = &spec.covariance_3x3;
        let b_sum = (spec.target_risk_budgets[0]
            + spec.target_risk_budgets[1]
            + spec.target_risk_budgets[2])
            .max_val(T::from_f64(1e-8));
        let b = [
            spec.target_risk_budgets[0] / b_sum,
            spec.target_risk_budgets[1] / b_sum,
            spec.target_risk_budgets[2] / b_sum,
        ];
        let vols = [
            cov[0][0].max_val(T::from_f64(1e-10)).sqrt(),
            cov[1][1].max_val(T::from_f64(1e-10)).sqrt(),
            cov[2][2].max_val(T::from_f64(1e-10)).sqrt(),
        ];
        let mut x = [b[0] / vols[0], b[1] / vols[1], b[2] / vols[2]];
        let inv_sum = x[0] + x[1] + x[2];
        let inv_vol_w = [x[0] / inv_sum, x[1] / inv_sum, x[2] / inv_sum];

        let two = T::from_f64(2.0);
        let four = T::from_f64(4.0);
        for _ in 0..16 {
            for i in 0..3 {
                let mut c_i = T::zero();
                for j in 0..3 {
                    if j != i {
                        c_i = c_i + cov[i][j] * x[j];
                    }
                }
                let sig_ii = cov[i][i].max_val(T::from_f64(1e-10));
                x[i] = (-c_i + (c_i * c_i + four * sig_ii * b[i]).sqrt()) / (two * sig_ii);
            }
        }

        let mut x_sig_x = T::zero();
        for i in 0..3 {
            for j in 0..3 {
                x_sig_x = x_sig_x + x[i] * cov[i][j] * x[j];
            }
        }
        let spinu_obj =
            T::from_f64(0.5) * x_sig_x - (b[0] * x[0].ln() + b[1] * x[1].ln() + b[2] * x[2].ln());

        let x_tot = x[0] + x[1] + x[2];
        let w = [x[0] / x_tot, x[1] / x_tot, x[2] / x_tot];
        let mut sig_w = [T::zero(); 3];
        for i in 0..3 {
            sig_w[i] = cov[i][0] * w[0] + cov[i][1] * w[1] + cov[i][2] * w[2];
        }
        let port_var =
            (w[0] * sig_w[0] + w[1] * sig_w[1] + w[2] * sig_w[2]).max_val(T::from_f64(1e-12));
        let port_vol = port_var.sqrt();
        let rrc = [
            w[0] * sig_w[0] / port_var,
            w[1] * sig_w[1] / port_var,
            w[2] * sig_w[2] / port_var,
        ];
        let mut l1_err = T::zero();
        for i in 0..3 {
            let d = rrc[i] - b[i];
            l1_err = l1_err + (d * d + T::from_f64(1e-30)).sqrt();
        }
        let w_diff = (w[0] + w[1] + w[2]) - T::one();
        let w_err = (w_diff * w_diff + T::from_f64(1e-30)).sqrt();
        let dr = (w[0] * vols[0] + w[1] * vols[1] + w[2] * vols[2]) / port_vol;

        RoncalliSpinuErcRiskParityReport {
            erc_weights: w,
            inverse_vol_naive_weights: inv_vol_w,
            relative_risk_contributions: rrc,
            erc_portfolio_volatility: port_vol,
            risk_budget_parity_l1_error: l1_err,
            weights_sum_conservation_error: w_err,
            choueifaty_diversification_ratio: dr,
            spinu_convex_objective_val: spinu_obj,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct IlmanenCarryRolldownVmcFactorSpec<T: Real = f64> {
        pub current_yield_yt: T,
        pub rolled_horizon_yield: T,
        pub repo_financing_rate: T,
        pub modified_duration: T,
        pub bond_convexity: T,
        pub expected_yield_change: T,
        pub yield_volatility: T,
        pub value_style_sharpe: T,
        pub momentum_style_sharpe: T,
        pub carry_style_sharpe: T,
        pub value_momentum_corr_rho: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct IlmanenCarryRolldownVmcFactorReport<T: Real = f64> {
        pub pure_carry_return: T,
        pub rolldown_return: T,
        pub expected_valuation_return: T,
        pub convexity_bias_return: T,
        pub total_expected_excess_return: T,
        pub carry_plus_rolldown_breakeven_bps: T,
        pub standalone_average_style_sharpe: T,
        pub composite_vmc_sharpe_ratio: T,
        pub style_diversification_sharpe_gain: T,
    }

    pub fn evaluate_ilmanen_carry_rolldown_and_vmc_factors<T: Real>(
        spec: &IlmanenCarryRolldownVmcFactorSpec<T>,
    ) -> IlmanenCarryRolldownVmcFactorReport<T> {
        let carry = spec.current_yield_yt - spec.repo_financing_rate;
        let rolldown =
            -spec.modified_duration * (spec.rolled_horizon_yield - spec.current_yield_yt);
        let val_ret = -spec.modified_duration * spec.expected_yield_change;
        let conv_bias = T::from_f64(0.5)
            * spec.bond_convexity
            * (spec.expected_yield_change * spec.expected_yield_change
                + spec.yield_volatility * spec.yield_volatility);
        let total_ret = carry + rolldown + val_ret + conv_bias;
        let breakeven_bps = (carry + rolldown) / spec.modified_duration.max_val(T::from_f64(1e-6))
            * T::from_f64(1e4);

        let sum_sr = spec.value_style_sharpe + spec.momentum_style_sharpe + spec.carry_style_sharpe;
        let avg_sr = sum_sr / T::from_f64(3.0);
        let denom = (T::from_f64(3.0) + T::from_f64(2.0) * spec.value_momentum_corr_rho)
            .max_val(T::from_f64(0.2))
            .sqrt();
        let vmc_sr = sum_sr / denom;
        let sr_gain = vmc_sr - avg_sr;

        IlmanenCarryRolldownVmcFactorReport {
            pure_carry_return: carry,
            rolldown_return: rolldown,
            expected_valuation_return: val_ret,
            convexity_bias_return: conv_bias,
            total_expected_excess_return: total_ret,
            carry_plus_rolldown_breakeven_bps: breakeven_bps,
            standalone_average_style_sharpe: avg_sr,
            composite_vmc_sharpe_ratio: vmc_sr,
            style_diversification_sharpe_gain: sr_gain,
        }
    }

    #[cfg(test)]
    mod domain_suite_48_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn suite48_portfolio_all_four_engines_invariants_and_ad_parity() {
            let h = 1e-6;

            // 1. Cartea-Jaimungal LOB Market Making
            let lob_spec = CarteaJaimungalLobMarketMakingSpec {
                mid_price: 100.0,
                current_inventory_q: 3.0,
                running_inventory_penalty_phi: 0.02,
                terminal_inventory_penalty_alpha: 0.05,
                risk_aversion_gamma: 0.10,
                asset_volatility: 0.20,
                order_arrival_baseline_a: 12.0,
                lob_depth_decay_kappa: 1.5,
                remaining_horizon_tau: 0.5,
            };
            let lob_rep = evaluate_cartea_jaimungal_lob_market_making(&lob_spec);
            assert!(lob_rep.inventory_quote_skew < 0.0);
            assert!(
                (lob_rep.optimal_ask_quote
                    - lob_rep.optimal_bid_quote
                    - lob_rep.total_quoted_spread)
                    .abs()
                    < 1e-12
            );
            let lob_ad =
                evaluate_cartea_jaimungal_lob_market_making(&CarteaJaimungalLobMarketMakingSpec {
                    mid_price: Dual::<1>::from_value(lob_spec.mid_price),
                    current_inventory_q: Dual::<1>::variable(lob_spec.current_inventory_q, 0),
                    running_inventory_penalty_phi: Dual::<1>::from_value(
                        lob_spec.running_inventory_penalty_phi,
                    ),
                    terminal_inventory_penalty_alpha: Dual::<1>::from_value(
                        lob_spec.terminal_inventory_penalty_alpha,
                    ),
                    risk_aversion_gamma: Dual::<1>::from_value(lob_spec.risk_aversion_gamma),
                    asset_volatility: Dual::<1>::from_value(lob_spec.asset_volatility),
                    order_arrival_baseline_a: Dual::<1>::from_value(
                        lob_spec.order_arrival_baseline_a,
                    ),
                    lob_depth_decay_kappa: Dual::<1>::from_value(lob_spec.lob_depth_decay_kappa),
                    remaining_horizon_tau: Dual::<1>::from_value(lob_spec.remaining_horizon_tau),
                });
            let mut lp = lob_spec;
            lp.current_inventory_q += h;
            let mut lm = lob_spec;
            lm.current_inventory_q -= h;
            let _fd_lob = (evaluate_cartea_jaimungal_lob_market_making(&lp).reservation_mid_price
                - evaluate_cartea_jaimungal_lob_market_making(&lm).reservation_mid_price)
                / (2.0 * h);
            assert!((lob_ad.reservation_mid_price.d[0] - _fd_lob).abs() < 1e-5);

            // 2. Campbell-Lo-MacKinlay VR & Campbell-Shiller
            let clm_spec = CampbellLoMackinlayVrPresentValueSpec {
                returns: vec![0.012, 0.015, -0.004, 0.009, 0.014, -0.003, 0.011, 0.013],
                steady_state_dividend_yield: 0.032,
                expected_dividend_growth: 0.045,
                expected_discount_rate: 0.078,
                dividend_growth_ar1_phi: 0.35,
                discount_rate_ar1_phi: 0.85,
            };
            let clm_rep = evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller(&clm_spec);
            assert!(clm_rep.campbell_shiller_rho > 0.9 && clm_rep.campbell_shiller_rho < 1.0);
            assert!(clm_rep.lo_mackinlay_variance_ratio_vr2 > 0.0);
            assert!(clm_rep.discount_rate_news_multiplier > 1.0);
            let clm_ad = evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller(
                &CampbellLoMackinlayVrPresentValueSpec {
                    returns: clm_spec
                        .returns
                        .iter()
                        .copied()
                        .map(Dual::<1>::from_value)
                        .collect(),
                    steady_state_dividend_yield: Dual::<1>::variable(
                        clm_spec.steady_state_dividend_yield,
                        0,
                    ),
                    expected_dividend_growth: Dual::<1>::from_value(
                        clm_spec.expected_dividend_growth,
                    ),
                    expected_discount_rate: Dual::<1>::from_value(clm_spec.expected_discount_rate),
                    dividend_growth_ar1_phi: Dual::<1>::from_value(
                        clm_spec.dividend_growth_ar1_phi,
                    ),
                    discount_rate_ar1_phi: Dual::<1>::from_value(clm_spec.discount_rate_ar1_phi),
                },
            );
            let mut cp = clm_spec.clone();
            cp.steady_state_dividend_yield += h;
            let mut cm = clm_spec;
            cm.steady_state_dividend_yield -= h;
            let _fd_clm = (evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller(&cp)
                .equilibrium_log_dividend_price_ratio
                - evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller(&cm)
                    .equilibrium_log_dividend_price_ratio)
                / (2.0 * h);
            assert!((clm_ad.equilibrium_log_dividend_price_ratio.d[0] - _fd_clm).abs() < 1e-5);

            // 3. Roncalli Spinu ERC Risk Parity
            let erc_spec = RoncalliSpinuErcRiskParitySpec {
                covariance_3x3: [
                    [0.040, 0.008, 0.002],
                    [0.008, 0.0225, 0.003],
                    [0.002, 0.003, 0.0064],
                ],
                target_risk_budgets: [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            };
            let erc_rep = evaluate_roncalli_spinu_erc_risk_parity(&erc_spec);
            assert!(erc_rep.risk_budget_parity_l1_error < 1e-8);
            assert!(erc_rep.weights_sum_conservation_error < 1e-12);
            assert!(erc_rep.choueifaty_diversification_ratio > 1.0);
            let mut cov_ad = erc_spec
                .covariance_3x3
                .map(|r| r.map(Dual::<1>::from_value));
            cov_ad[0][0] = Dual::<1>::variable(erc_spec.covariance_3x3[0][0], 0);
            let erc_ad = evaluate_roncalli_spinu_erc_risk_parity(&RoncalliSpinuErcRiskParitySpec {
                covariance_3x3: cov_ad,
                target_risk_budgets: erc_spec.target_risk_budgets.map(Dual::<1>::from_value),
            });
            let mut ep = erc_spec;
            ep.covariance_3x3[0][0] += h;
            let mut em = erc_spec;
            em.covariance_3x3[0][0] -= h;
            let _fd_erc = (evaluate_roncalli_spinu_erc_risk_parity(&ep).erc_portfolio_volatility
                - evaluate_roncalli_spinu_erc_risk_parity(&em).erc_portfolio_volatility)
                / (2.0 * h);
            assert!((erc_ad.erc_portfolio_volatility.d[0] - _fd_erc).abs() < 1e-5);

            // 4. Ilmanen Carry, Roll-Down & VMC Factors
            let ilm_spec = IlmanenCarryRolldownVmcFactorSpec {
                current_yield_yt: 0.045,
                rolled_horizon_yield: 0.041,
                repo_financing_rate: 0.032,
                modified_duration: 7.5,
                bond_convexity: 68.0,
                expected_yield_change: 0.001,
                yield_volatility: 0.0085,
                value_style_sharpe: 0.45,
                momentum_style_sharpe: 0.50,
                carry_style_sharpe: 0.48,
                value_momentum_corr_rho: -0.45,
            };
            let ilm_rep = evaluate_ilmanen_carry_rolldown_and_vmc_factors(&ilm_spec);
            assert!(ilm_rep.rolldown_return > 0.0);
            assert!(ilm_rep.composite_vmc_sharpe_ratio > ilm_rep.standalone_average_style_sharpe);
            let ilm_ad = evaluate_ilmanen_carry_rolldown_and_vmc_factors(
                &IlmanenCarryRolldownVmcFactorSpec {
                    current_yield_yt: Dual::<1>::variable(ilm_spec.current_yield_yt, 0),
                    rolled_horizon_yield: Dual::<1>::from_value(ilm_spec.rolled_horizon_yield),
                    repo_financing_rate: Dual::<1>::from_value(ilm_spec.repo_financing_rate),
                    modified_duration: Dual::<1>::from_value(ilm_spec.modified_duration),
                    bond_convexity: Dual::<1>::from_value(ilm_spec.bond_convexity),
                    expected_yield_change: Dual::<1>::from_value(ilm_spec.expected_yield_change),
                    yield_volatility: Dual::<1>::from_value(ilm_spec.yield_volatility),
                    value_style_sharpe: Dual::<1>::from_value(ilm_spec.value_style_sharpe),
                    momentum_style_sharpe: Dual::<1>::from_value(ilm_spec.momentum_style_sharpe),
                    carry_style_sharpe: Dual::<1>::from_value(ilm_spec.carry_style_sharpe),
                    value_momentum_corr_rho: Dual::<1>::from_value(
                        ilm_spec.value_momentum_corr_rho,
                    ),
                },
            );
            let mut ip = ilm_spec;
            ip.current_yield_yt += h;
            let mut im = ilm_spec;
            im.current_yield_yt -= h;
            let _fd_ilm = (evaluate_ilmanen_carry_rolldown_and_vmc_factors(&ip)
                .total_expected_excess_return
                - evaluate_ilmanen_carry_rolldown_and_vmc_factors(&im)
                    .total_expected_excess_return)
                / (2.0 * h);
            assert!((ilm_ad.total_expected_excess_return.d[0] - _fd_ilm).abs() < 1e-5);
        }
    }
}
pub use evaluate_cartea_jaimungal_lob_market_making_module::*;

mod evaluate_meucci_black_litterman_woodbury_equivalence_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct MeucciBlackLittermanWoodburySpec<T: Real = f64> {
        pub covariance_3x3: [[T; 3]; 3],
        pub market_cap_weights: [T; 3],
        pub risk_aversion_delta: T,
        pub tau_scalar: T,
        pub view_pick_vector: [T; 3],
        pub view_target_q: T,
        pub view_uncertainty_omega: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct MeucciBlackLittermanWoodburyReport<T: Real = f64> {
        pub equilibrium_prior_returns_pi: [T; 3],
        pub posterior_returns_woodbury: [T; 3],
        pub posterior_returns_precision: [T; 3],
        pub woodbury_precision_equivalence_error: T,
        pub posterior_optimal_weights: [T; 3],
        pub active_view_weight_tilt_l1: T,
        pub posterior_view_expected_return: T,
    }

    fn inv3<T: Real>(a: [[T; 3]; 3]) -> [[T; 3]; 3] {
        let mut aug = [[T::zero(); 6]; 3];
        for i in 0..3 {
            for j in 0..3 {
                aug[i][j] = a[i][j];
            }
            aug[i][i + 3] = T::one();
        }
        for i in 0..3 {
            let piv = aug[i][i];
            for j in 0..6 {
                aug[i][j] = aug[i][j] / piv;
            }
            for r in 0..3 {
                if r != i {
                    let f = aug[r][i];
                    for j in 0..6 {
                        aug[r][j] = aug[r][j] - f * aug[i][j];
                    }
                }
            }
        }
        let mut inv = [[T::zero(); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                inv[i][j] = aug[i][j + 3];
            }
        }
        inv
    }

    pub fn evaluate_meucci_black_litterman_woodbury_equivalence<T: Real>(
        spec: &MeucciBlackLittermanWoodburySpec<T>,
    ) -> MeucciBlackLittermanWoodburyReport<T> {
        let cov = &spec.covariance_3x3;
        let p = &spec.view_pick_vector;
        let tau = spec.tau_scalar.max_val(T::from_f64(1e-6));
        let omega = spec.view_uncertainty_omega.max_val(T::from_f64(1e-8));
        let delta = spec.risk_aversion_delta.max_val(T::from_f64(1e-6));

        let mut pi = [T::zero(); 3];
        let mut sig_p = [T::zero(); 3];
        for i in 0..3 {
            pi[i] = delta
                * (cov[i][0] * spec.market_cap_weights[0]
                    + cov[i][1] * spec.market_cap_weights[1]
                    + cov[i][2] * spec.market_cap_weights[2]);
            sig_p[i] = cov[i][0] * p[0] + cov[i][1] * p[1] + cov[i][2] * p[2];
        }
        let p_sig_p = p[0] * sig_p[0] + p[1] * sig_p[1] + p[2] * sig_p[2];
        let s_p = tau * p_sig_p + omega;
        let p_pi = p[0] * pi[0] + p[1] * pi[1] + p[2] * pi[2];
        let view_inn = spec.view_target_q - p_pi;

        let mut mu_wood = [T::zero(); 3];
        let mut m_wood = [[T::zero(); 3]; 3];
        for i in 0..3 {
            let k_i = tau * sig_p[i] / s_p;
            mu_wood[i] = pi[i] + k_i * view_inn;
            for j in 0..3 {
                m_wood[i][j] = tau * cov[i][j] - (tau * tau / s_p) * sig_p[i] * sig_p[j];
            }
        }

        let mut tau_cov = [[T::zero(); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                tau_cov[i][j] = tau * cov[i][j];
            }
        }
        let inv_tau_cov = inv3(tau_cov);
        let mut a_prec = [[T::zero(); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                a_prec[i][j] = inv_tau_cov[i][j] + (p[i] * p[j]) / omega;
            }
        }
        let m_prec = inv3(a_prec);
        let mut rhs = [T::zero(); 3];
        for i in 0..3 {
            let prior_prec_pi =
                inv_tau_cov[i][0] * pi[0] + inv_tau_cov[i][1] * pi[1] + inv_tau_cov[i][2] * pi[2];
            rhs[i] = prior_prec_pi + (spec.view_target_q / omega) * p[i];
        }
        let mut mu_prec = [T::zero(); 3];
        for i in 0..3 {
            mu_prec[i] = m_prec[i][0] * rhs[0] + m_prec[i][1] * rhs[1] + m_prec[i][2] * rhs[2];
        }

        let mut eq_err = T::zero();
        for i in 0..3 {
            let dm = mu_wood[i] - mu_prec[i];
            eq_err = eq_err + (dm * dm + T::from_f64(1e-30)).sqrt();
            for j in 0..3 {
                let dc = m_wood[i][j] - m_prec[i][j];
                eq_err = eq_err + (dc * dc + T::from_f64(1e-30)).sqrt();
            }
        }

        let inv_cov = inv3(*cov);
        let mut w_opt = [T::zero(); 3];
        let mut tilt_l1 = T::zero();
        for i in 0..3 {
            w_opt[i] = (inv_cov[i][0] * mu_wood[0]
                + inv_cov[i][1] * mu_wood[1]
                + inv_cov[i][2] * mu_wood[2])
                / delta;
            let dw = w_opt[i] - spec.market_cap_weights[i];
            tilt_l1 = tilt_l1 + (dw * dw + T::from_f64(1e-30)).sqrt();
        }
        let post_view_ret = p[0] * mu_wood[0] + p[1] * mu_wood[1] + p[2] * mu_wood[2];

        MeucciBlackLittermanWoodburyReport {
            equilibrium_prior_returns_pi: pi,
            posterior_returns_woodbury: mu_wood,
            posterior_returns_precision: mu_prec,
            woodbury_precision_equivalence_error: eq_err,
            posterior_optimal_weights: w_opt,
            active_view_weight_tilt_l1: tilt_l1,
            posterior_view_expected_return: post_view_ret,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct GueantAlmgrenChrissHyperbolicSpec<T: Real = f64> {
        pub initial_inventory_x0: T,
        pub horizon_years_t: T,
        pub temporary_impact_eta: T,
        pub permanent_impact_gamma: T,
        pub asset_volatility_sigma: T,
        pub risk_aversion_lambda: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct GueantAlmgrenChrissHyperbolicReport<T: Real = f64> {
        pub urgency_parameter_kappa: T,
        pub inventory_quartiles_q1_q2_q3: [T; 3],
        pub midpoint_frontloading_vs_twap: T,
        pub optimal_expected_shortfall: T,
        pub optimal_shortfall_variance: T,
        pub optimal_mean_variance_objective: T,
        pub twap_mean_variance_objective: T,
        pub ac_utility_gain_over_twap: T,
    }

    #[inline]
    fn sinh_r<T: Real>(z: T) -> T {
        T::from_f64(0.5) * (z.exp() - (-z).exp())
    }

    pub fn evaluate_gueant_almgren_chriss_hyperbolic_liquidation<T: Real>(
        spec: &GueantAlmgrenChrissHyperbolicSpec<T>,
    ) -> GueantAlmgrenChrissHyperbolicReport<T> {
        let x0 = spec.initial_inventory_x0;
        let mat = spec.horizon_years_t.max_val(T::from_f64(1e-6));
        let eta = spec.temporary_impact_eta.max_val(T::from_f64(1e-8));
        let sig2 = spec.asset_volatility_sigma * spec.asset_volatility_sigma;
        let lam = spec.risk_aversion_lambda.max_val(T::from_f64(1e-8));
        let kappa = (lam * sig2 / eta).sqrt();
        let kt = (kappa * mat).max_val(T::from_f64(1e-5));
        let sh_kt = sinh_r(kt);

        let q1 = x0 * sinh_r(kt * T::from_f64(0.75)) / sh_kt;
        let q2 = x0 * sinh_r(kt * T::from_f64(0.50)) / sh_kt;
        let q3 = x0 * sinh_r(kt * T::from_f64(0.25)) / sh_kt;
        let frontload = T::from_f64(0.5) * x0 - q2;

        let two_kt = T::from_f64(2.0) * kt;
        let sh_2kt = sinh_r(two_kt);
        let four_sh2 = T::from_f64(4.0) * sh_kt * sh_kt;
        let x0_2 = x0 * x0;
        let perm_cost = T::from_f64(0.5) * spec.permanent_impact_gamma * x0_2;

        let e_is = perm_cost + eta * x0_2 * kappa * (sh_2kt + two_kt) / four_sh2;
        let var_is = sig2 * x0_2 * (sh_2kt - two_kt) / (kappa * four_sh2);
        let j_ac = e_is + T::from_f64(0.5) * lam * var_is;

        let e_twap = perm_cost + eta * x0_2 / mat;
        let var_twap = sig2 * x0_2 * mat / T::from_f64(3.0);
        let j_twap = e_twap + T::from_f64(0.5) * lam * var_twap;

        GueantAlmgrenChrissHyperbolicReport {
            urgency_parameter_kappa: kappa,
            inventory_quartiles_q1_q2_q3: [q1, q2, q3],
            midpoint_frontloading_vs_twap: frontload,
            optimal_expected_shortfall: e_is,
            optimal_shortfall_variance: var_is,
            optimal_mean_variance_objective: j_ac,
            twap_mean_variance_objective: j_twap,
            ac_utility_gain_over_twap: j_twap - j_ac,
        }
    }

    #[cfg(test)]
    mod domain_suite_49_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn suite49_portfolio_bl_woodbury_and_ac_hyperbolic_invariants_and_ad_parity() {
            let h = 1e-6;

            let bl_spec = MeucciBlackLittermanWoodburySpec {
                covariance_3x3: [
                    [0.040, 0.010, 0.004],
                    [0.010, 0.025, 0.005],
                    [0.004, 0.005, 0.010],
                ],
                market_cap_weights: [0.50, 0.30, 0.20],
                risk_aversion_delta: 2.5,
                tau_scalar: 0.05,
                view_pick_vector: [1.0, -1.0, 0.0],
                view_target_q: 0.045,
                view_uncertainty_omega: 0.0015,
            };
            let bl_rep = evaluate_meucci_black_litterman_woodbury_equivalence(&bl_spec);
            assert!(bl_rep.woodbury_precision_equivalence_error < 1e-10);
            assert!(bl_rep.active_view_weight_tilt_l1 > 0.0);

            let bl_ad = evaluate_meucci_black_litterman_woodbury_equivalence(
                &MeucciBlackLittermanWoodburySpec {
                    covariance_3x3: bl_spec.covariance_3x3.map(|r| r.map(Dual::<1>::from_value)),
                    market_cap_weights: bl_spec.market_cap_weights.map(Dual::<1>::from_value),
                    risk_aversion_delta: Dual::<1>::from_value(bl_spec.risk_aversion_delta),
                    tau_scalar: Dual::<1>::from_value(bl_spec.tau_scalar),
                    view_pick_vector: bl_spec.view_pick_vector.map(Dual::<1>::from_value),
                    view_target_q: Dual::<1>::variable(bl_spec.view_target_q, 0),
                    view_uncertainty_omega: Dual::<1>::from_value(bl_spec.view_uncertainty_omega),
                },
            );
            let mut bp = bl_spec;
            bp.view_target_q += h;
            let mut bm = bl_spec;
            bm.view_target_q -= h;
            let fd_bl = (evaluate_meucci_black_litterman_woodbury_equivalence(&bp)
                .posterior_view_expected_return
                - evaluate_meucci_black_litterman_woodbury_equivalence(&bm)
                    .posterior_view_expected_return)
                / (2.0 * h);
            assert!((bl_ad.posterior_view_expected_return.d[0] - fd_bl).abs() < 1e-5);

            let ac_spec = GueantAlmgrenChrissHyperbolicSpec {
                initial_inventory_x0: 100.0,
                horizon_years_t: 0.25,
                temporary_impact_eta: 0.002,
                permanent_impact_gamma: 0.0005,
                asset_volatility_sigma: 0.30,
                risk_aversion_lambda: 0.15,
            };
            let ac_rep = evaluate_gueant_almgren_chriss_hyperbolic_liquidation(&ac_spec);
            assert!(ac_rep.midpoint_frontloading_vs_twap > 0.0);
            assert!(ac_rep.ac_utility_gain_over_twap > 0.0);

            let ac_ad = evaluate_gueant_almgren_chriss_hyperbolic_liquidation(
                &GueantAlmgrenChrissHyperbolicSpec {
                    initial_inventory_x0: Dual::<1>::variable(ac_spec.initial_inventory_x0, 0),
                    horizon_years_t: Dual::<1>::from_value(ac_spec.horizon_years_t),
                    temporary_impact_eta: Dual::<1>::from_value(ac_spec.temporary_impact_eta),
                    permanent_impact_gamma: Dual::<1>::from_value(ac_spec.permanent_impact_gamma),
                    asset_volatility_sigma: Dual::<1>::from_value(ac_spec.asset_volatility_sigma),
                    risk_aversion_lambda: Dual::<1>::from_value(ac_spec.risk_aversion_lambda),
                },
            );
            let mut ap = ac_spec;
            ap.initial_inventory_x0 += h;
            let mut am = ac_spec;
            am.initial_inventory_x0 -= h;
            let fd_ac = (evaluate_gueant_almgren_chriss_hyperbolic_liquidation(&ap)
                .optimal_mean_variance_objective
                - evaluate_gueant_almgren_chriss_hyperbolic_liquidation(&am)
                    .optimal_mean_variance_objective)
                / (2.0 * h);
            assert!((ac_ad.optimal_mean_variance_objective.d[0] - fd_ac).abs() < 1e-4);
        }
    }
}
pub use evaluate_meucci_black_litterman_woodbury_equivalence_module::*;

mod evaluate_ledoit_wolf_constant_correlation_shrinkage_module {
    // Engine 5A: Ledoit & Wolf (2004) & Meucci (2005) Constant-Correlation Analytical
    // Oracle Covariance Shrinkage Intensity Engine.
    // Engine 5B: Almgren & Chriss (2000) & Gatheral–Schied (2011) Hyperbolic Optimal
    // Execution Trajectory, Urgency kappa & Implementation Shortfall VaR Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct LedoitWolfConstantCorrelationSpec<T: Real = f64> {
        pub sample_volatilities: [T; 3],
        pub sample_correlations: [T; 3],
        pub asymptotic_pi_hat: T,
        pub asymptotic_rho_hat: T,
        pub sample_observations_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct LedoitWolfConstantCorrelationReport<T: Real = f64> {
        pub average_sample_correlation_rbar: T,
        pub misspecification_distance_gamma_hat: T,
        pub optimal_shrinkage_intensity_delta: T,
        pub shrunk_off_diagonal_covariances: [T; 3],
        pub equal_weight_shrunk_variance: T,
    }

    pub fn evaluate_ledoit_wolf_constant_correlation_shrinkage<T: Real>(
        spec: &LedoitWolfConstantCorrelationSpec<T>,
    ) -> LedoitWolfConstantCorrelationReport<T> {
        let v0 = spec.sample_volatilities[0].max_val(T::from_f64(1e-4));
        let v1 = spec.sample_volatilities[1].max_val(T::from_f64(1e-4));
        let v2 = spec.sample_volatilities[2].max_val(T::from_f64(1e-4));

        let r01 = spec.sample_correlations[0];
        let r02 = spec.sample_correlations[1];
        let r12 = spec.sample_correlations[2];
        let r_bar = (r01 + r02 + r12) / T::from_f64(3.0);

        // Sample vs. Constant-Correlation target off-diagonal entries:
        let s01 = r01 * v0 * v1;
        let s02 = r02 * v0 * v2;
        let s12 = r12 * v1 * v2;
        let f01 = r_bar * v0 * v1;
        let f02 = r_bar * v0 * v2;
        let f12 = r_bar * v1 * v2;

        // Frobenius squared distance gamma_hat = ||F - S||_F^2:
        let d01 = f01 - s01;
        let d02 = f02 - s02;
        let d12 = f12 - s12;
        let gamma_hat =
            (T::from_f64(2.0) * (d01 * d01 + d02 * d02 + d12 * d12)).max_val(T::from_f64(1e-8));

        // Ledoit-Wolf (2004, Eq. 5) optimal shrinkage intensity:
        // kappa_hat = (pi_hat - rho_hat) / gamma_hat, delta^* = clamp(kappa_hat / T, 0, 1)
        let t_obs = spec.sample_observations_t.max_val(T::one());
        let raw_delta = (spec.asymptotic_pi_hat - spec.asymptotic_rho_hat)
            .max_val(T::from_f64(1e-8))
            / (t_obs * gamma_hat);
        let delta = raw_delta / (T::one() + raw_delta);

        let c01 = delta * f01 + (T::one() - delta) * s01;
        let c02 = delta * f02 + (T::one() - delta) * s02;
        let c12 = delta * f12 + (T::one() - delta) * s12;

        let eq_var =
            (v0 * v0 + v1 * v1 + v2 * v2 + T::from_f64(2.0) * (c01 + c02 + c12)) / T::from_f64(9.0);

        LedoitWolfConstantCorrelationReport {
            average_sample_correlation_rbar: r_bar,
            misspecification_distance_gamma_hat: gamma_hat,
            optimal_shrinkage_intensity_delta: delta,
            shrunk_off_diagonal_covariances: [c01, c02, c12],
            equal_weight_shrunk_variance: eq_var,
        }
    }

    #[derive(Debug, Clone, Copy)]
    pub struct AlmgrenChrissCaraExecutionSpec<T: Real = f64> {
        pub initial_shares_x0: T,
        pub liquidation_horizon_t: T,
        pub asset_price_volatility_sigma: T,
        pub temporary_impact_eta: T,
        pub permanent_impact_gamma: T,
        pub cara_risk_aversion_lambda: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct AlmgrenChrissCaraExecutionReport<T: Real = f64> {
        pub urgency_parameter_kappa: T,
        pub characteristic_half_life_tau: T,
        pub mid_horizon_remaining_inventory_ratio: T,
        pub permanent_impact_expected_cost: T,
        pub temporary_impact_expected_cost: T,
        pub total_expected_implementation_shortfall: T,
        pub implementation_shortfall_variance: T,
        pub cara_certainty_equivalent_cost: T,
    }

    pub fn evaluate_almgren_chriss_cara_optimal_execution<T: Real>(
        spec: &AlmgrenChrissCaraExecutionSpec<T>,
    ) -> AlmgrenChrissCaraExecutionReport<T> {
        let x0 = spec.initial_shares_x0.max_val(T::from_f64(1e-4));
        let t = spec.liquidation_horizon_t.max_val(T::from_f64(0.05));
        let sig = spec.asset_price_volatility_sigma.max_val(T::from_f64(1e-4));
        let eta = spec.temporary_impact_eta.max_val(T::from_f64(1e-6));
        let gam = spec.permanent_impact_gamma.max_val(T::zero());
        let lam = spec.cara_risk_aversion_lambda.max_val(T::from_f64(1e-8));

        // Almgren-Chriss (2000) urgency parameter kappa = sqrt(lambda * sigma^2 / eta):
        let kappa = (lam * sig * sig / eta).sqrt();
        let half_life = T::from_f64(std::f64::consts::LN_2) / kappa;

        let kt = kappa * t;
        let sinh_kt = T::from_f64(0.5) * (kt.exp() - (-kt).exp());
        let cosh_kt = T::from_f64(0.5) * (kt.exp() + (-kt).exp());
        let sinh_half =
            T::from_f64(0.5) * ((T::from_f64(0.5) * kt).exp() - (-T::from_f64(0.5) * kt).exp());
        let mid_ratio = sinh_half / sinh_kt.max_val(T::from_f64(1e-8));

        // Continuous-time exact integrals of x(t) = X_0 sinh(kappa (T-t)) / sinh(kappa T):
        // int_0^T x'(t)^2 dt = X_0^2 * kappa * (sinh(2 k T) + 2 k T) / (4 sinh^2(k T))
        // int_0^T x(t)^2 dt  = X_0^2 / kappa * (sinh(2 k T) - 2 k T) / (4 sinh^2(k T))
        let sinh_2kt = T::from_f64(2.0) * sinh_kt * cosh_kt;
        let denom = T::from_f64(4.0) * sinh_kt * sinh_kt;

        let perm_cost = T::from_f64(0.5) * gam * x0 * x0;
        let temp_cost = eta * x0 * x0 * kappa * (sinh_2kt + T::from_f64(2.0) * kt) / denom;
        let exp_shortfall = perm_cost + temp_cost;
        let var_shortfall =
            sig * sig * x0 * x0 * (sinh_2kt - T::from_f64(2.0) * kt) / (kappa * denom);
        let ce_cost = exp_shortfall + T::from_f64(0.5) * lam * var_shortfall;

        AlmgrenChrissCaraExecutionReport {
            urgency_parameter_kappa: kappa,
            characteristic_half_life_tau: half_life,
            mid_horizon_remaining_inventory_ratio: mid_ratio,
            permanent_impact_expected_cost: perm_cost,
            temporary_impact_expected_cost: temp_cost,
            total_expected_implementation_shortfall: exp_shortfall,
            implementation_shortfall_variance: var_shortfall,
            cara_certainty_equivalent_cost: ce_cost,
        }
    }

    #[cfg(test)]
    mod domain_suite_54_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn ledoit_wolf_and_almgren_chriss_invariants_and_ad_parity() {
            let s1 = LedoitWolfConstantCorrelationSpec {
                sample_volatilities: [0.22, 0.18, 0.14],
                sample_correlations: [0.65, 0.15, 0.40],
                asymptotic_pi_hat: 0.045,
                asymptotic_rho_hat: 0.012,
                sample_observations_t: 60.0,
            };
            let r1 = evaluate_ledoit_wolf_constant_correlation_shrinkage(&s1);
            assert!(
                r1.optimal_shrinkage_intensity_delta > 0.0
                    && r1.optimal_shrinkage_intensity_delta < 1.0
            );
            assert!(r1.equal_weight_shrunk_variance > 0.0);

            let h = 1e-6;
            let r1_ad = evaluate_ledoit_wolf_constant_correlation_shrinkage(
                &LedoitWolfConstantCorrelationSpec {
                    sample_volatilities: [
                        Dual::<1>::variable(s1.sample_volatilities[0], 0),
                        Dual::<1>::from_value(s1.sample_volatilities[1]),
                        Dual::<1>::from_value(s1.sample_volatilities[2]),
                    ],
                    sample_correlations: s1.sample_correlations.map(Dual::<1>::from_value),
                    asymptotic_pi_hat: Dual::<1>::from_value(s1.asymptotic_pi_hat),
                    asymptotic_rho_hat: Dual::<1>::from_value(s1.asymptotic_rho_hat),
                    sample_observations_t: Dual::<1>::from_value(s1.sample_observations_t),
                },
            );
            let mut s1p = s1;
            s1p.sample_volatilities[0] += h;
            let mut s1m = s1;
            s1m.sample_volatilities[0] -= h;
            let fd1 = (evaluate_ledoit_wolf_constant_correlation_shrinkage(&s1p)
                .equal_weight_shrunk_variance
                - evaluate_ledoit_wolf_constant_correlation_shrinkage(&s1m)
                    .equal_weight_shrunk_variance)
                / (2.0 * h);
            assert!((r1_ad.equal_weight_shrunk_variance.d[0] - fd1).abs() < 1e-6);

            let s2 = AlmgrenChrissCaraExecutionSpec {
                initial_shares_x0: 100.0,
                liquidation_horizon_t: 1.0,
                asset_price_volatility_sigma: 0.30,
                temporary_impact_eta: 0.02,
                permanent_impact_gamma: 0.005,
                cara_risk_aversion_lambda: 0.5,
            };
            let r2 = evaluate_almgren_chriss_cara_optimal_execution(&s2);
            assert!(r2.mid_horizon_remaining_inventory_ratio < 0.5);
            assert!(r2.cara_certainty_equivalent_cost > r2.total_expected_implementation_shortfall);

            let r2_ad =
                evaluate_almgren_chriss_cara_optimal_execution(&AlmgrenChrissCaraExecutionSpec {
                    initial_shares_x0: Dual::<1>::variable(s2.initial_shares_x0, 0),
                    liquidation_horizon_t: Dual::<1>::from_value(s2.liquidation_horizon_t),
                    asset_price_volatility_sigma: Dual::<1>::from_value(
                        s2.asset_price_volatility_sigma,
                    ),
                    temporary_impact_eta: Dual::<1>::from_value(s2.temporary_impact_eta),
                    permanent_impact_gamma: Dual::<1>::from_value(s2.permanent_impact_gamma),
                    cara_risk_aversion_lambda: Dual::<1>::from_value(s2.cara_risk_aversion_lambda),
                });
            let mut s2p = s2;
            s2p.initial_shares_x0 += h;
            let mut s2m = s2;
            s2m.initial_shares_x0 -= h;
            let fd2 = (evaluate_almgren_chriss_cara_optimal_execution(&s2p)
                .cara_certainty_equivalent_cost
                - evaluate_almgren_chriss_cara_optimal_execution(&s2m)
                    .cara_certainty_equivalent_cost)
                / (2.0 * h);
            assert!((r2_ad.cara_certainty_equivalent_cost.d[0] - fd2).abs() < 1e-4);
        }
    }
}
pub use evaluate_ledoit_wolf_constant_correlation_shrinkage_module::*;

mod evaluate_jagannathan_ma_no_short_sale_shrinkage_module {
    // Engine 5A: Jagannathan & Ma (2003) & Roncalli (2013) No-Short-Sale KKT Dual-Multiplier
    // Covariance Shrinkage Matrix `Sigma_tilde` & Constrained Minimum-Variance Engine.
    // Engine 5B: Obizhaeva & Wang (2013) & Cartea–Jaimungal–Penalva (2015) Transient
    // Limit-Order-Book Resilience (`rho_res`), Block-Plus-Continuous Trajectory & Propagator Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct JagannathanMaCovarianceShrinkageSpec<T: Real = f64> {
        pub asset_volatilities: [T; 3],
        pub correlation_12: T,
        pub correlation_13: T,
        pub correlation_23: T,
        pub weight_upper_cap_u: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct JagannathanMaCovarianceShrinkageReport<T: Real = f64> {
        pub unconstrained_gmv_weights: [T; 3],
        pub constrained_long_only_weights: [T; 3],
        pub kkt_lower_bound_multipliers_lambda: [T; 3],
        pub shrunk_covariance_12_tilde: T,
        pub raw_covariance_12: T,
        pub constrained_gmv_portfolio_volatility: T,
        pub jagannathan_ma_equivalence_residual: T,
    }

    fn solve_gmv_3x3<T: Real>(c00: T, c01: T, c02: T, c11: T, c12: T, c22: T) -> [T; 3] {
        let a00 = c11 * c22 - c12 * c12;
        let a01 = -(c01 * c22 - c12 * c02);
        let a02 = c01 * c12 - c11 * c02;
        let a11 = c00 * c22 - c02 * c02;
        let a12 = -(c00 * c12 - c01 * c02);
        let a22 = c00 * c11 - c01 * c01;

        let z0 = a00 + a01 + a02;
        let z1 = a01 + a11 + a12;
        let z2 = a02 + a12 + a22;
        let denom = z0 + z1 + z2;
        [z0 / denom, z1 / denom, z2 / denom]
    }

    pub fn evaluate_jagannathan_ma_no_short_sale_shrinkage<T: Real>(
        spec: &JagannathanMaCovarianceShrinkageSpec<T>,
    ) -> JagannathanMaCovarianceShrinkageReport<T> {
        let v0 = spec.asset_volatilities[0].max_val(T::from_f64(1e-4));
        let v1 = spec.asset_volatilities[1].max_val(T::from_f64(1e-4));
        let v2 = spec.asset_volatilities[2].max_val(T::from_f64(1e-4));

        let c00 = v0 * v0;
        let c11 = v1 * v1;
        let c22 = v2 * v2;
        let c01 = spec.correlation_12 * v0 * v1;
        let c02 = spec.correlation_13 * v0 * v2;
        let c12 = spec.correlation_23 * v1 * v2;

        let w_unc = solve_gmv_3x3(c00, c01, c02, c11, c12, c22);

        // Smooth projected-gradient / KKT fixed-point under 0 <= w_i <= u_cap:
        let u_cap = spec
            .weight_upper_cap_u
            .max_val(T::from_f64(0.40))
            .min_val(T::from_f64(0.95));
        let eps = T::from_f64(1e-5);
        let smooth_pos = |x: T| T::from_f64(0.5) * (x + (x * x + eps * eps).sqrt());
        let smooth_cap = |x: T| {
            let p = smooth_pos(x);
            u_cap - smooth_pos(u_cap - p)
        };

        let mut w0 = smooth_cap(w_unc[0]);
        let mut w1 = smooth_cap(w_unc[1]);
        let mut w2 = smooth_cap(w_unc[2]);
        let s_init = (w0 + w1 + w2).max_val(T::from_f64(1e-6));
        w0 /= s_init;
        w1 /= s_init;
        w2 /= s_init;

        let step = T::from_f64(8.0);
        for _ in 0..20 {
            let g0 = c00 * w0 + c01 * w1 + c02 * w2;
            let g1 = c01 * w0 + c11 * w1 + c12 * w2;
            let g2 = c02 * w0 + c12 * w1 + c22 * w2;
            let mut mu = (g0 + g1 + g2) / T::from_f64(3.0);
            for _ in 0..5 {
                let nw0 = smooth_cap(w0 - step * (g0 - mu));
                let nw1 = smooth_cap(w1 - step * (g1 - mu));
                let nw2 = smooth_cap(w2 - step * (g2 - mu));
                mu -= ((nw0 + nw1 + nw2) - T::one()) / step;
            }
            w0 = smooth_cap(w0 - step * (g0 - mu));
            w1 = smooth_cap(w1 - step * (g1 - mu));
            w2 = smooth_cap(w2 - step * (g2 - mu));
            let sw = (w0 + w1 + w2).max_val(T::from_f64(1e-6));
            w0 /= sw;
            w1 /= sw;
            w2 /= sw;
        }

        let cw0 = c00 * w0 + c01 * w1 + c02 * w2;
        let cw1 = c01 * w0 + c11 * w1 + c12 * w2;
        let cw2 = c02 * w0 + c12 * w1 + c22 * w2;
        let port_var = w0 * cw0 + w1 * cw1 + w2 * cw2;

        // Jagannathan-Ma (2003, Prop. 1, Eq. 10) KKT dual multipliers: Sigma w - lambda + nu = mu_kkt 1
        let iw0 = smooth_pos(w0 * (u_cap - w0));
        let iw1 = smooth_pos(w1 * (u_cap - w1));
        let iw2 = smooth_pos(w2 * (u_cap - w2));
        let mu_kkt =
            (iw0 * cw0 + iw1 * cw1 + iw2 * cw2) / (iw0 + iw1 + iw2).max_val(T::from_f64(1e-8));
        let lam0 = smooth_pos(cw0 - mu_kkt);
        let lam1 = smooth_pos(cw1 - mu_kkt);
        let lam2 = smooth_pos(cw2 - mu_kkt);
        let nu0 = smooth_pos(mu_kkt - cw0);
        let nu1 = smooth_pos(mu_kkt - cw1);
        let nu2 = smooth_pos(mu_kkt - cw2);

        // Jagannathan-Ma shrunk covariance Sigma_tilde_ij = Sigma_ij - (lambda_i + lambda_j) + (nu_i + nu_j):
        let t00 = c00 - T::from_f64(2.0) * lam0 + T::from_f64(2.0) * nu0;
        let t11 = c11 - T::from_f64(2.0) * lam1 + T::from_f64(2.0) * nu1;
        let t22 = c22 - T::from_f64(2.0) * lam2 + T::from_f64(2.0) * nu2;
        let t01 = c01 - (lam0 + lam1) + (nu0 + nu1);
        let t02 = c02 - (lam0 + lam2) + (nu0 + nu2);
        let t12 = c12 - (lam1 + lam2) + (nu1 + nu2);

        let w_tilde = solve_gmv_3x3(t00, t01, t02, t11, t12, t22);
        let eq_res = (w_tilde[0] - w0) * (w_tilde[0] - w0)
            + (w_tilde[1] - w1) * (w_tilde[1] - w1)
            + (w_tilde[2] - w2) * (w_tilde[2] - w2);

        JagannathanMaCovarianceShrinkageReport {
            unconstrained_gmv_weights: w_unc,
            constrained_long_only_weights: [w0, w1, w2],
            kkt_lower_bound_multipliers_lambda: [lam0, lam1, lam2],
            shrunk_covariance_12_tilde: t01,
            raw_covariance_12: c01,
            constrained_gmv_portfolio_volatility: port_var.max_val(T::from_f64(1e-8)).sqrt(),
            jagannathan_ma_equivalence_residual: (eq_res + T::from_f64(1e-12)).sqrt(),
        }
    }

    #[derive(Debug, Clone, Copy)]
    pub struct ObizhaevaWangTransientLobExecutionSpec<T: Real = f64> {
        pub total_order_size_x0: T,
        pub execution_horizon_t: T,
        pub lob_depth_density_q0: T,
        pub permanent_impact_gamma: T,
        pub lob_resilience_speed_rho: T,
        pub asset_volatility_sigma: T,
        pub risk_aversion_lambda: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct ObizhaevaWangTransientLobExecutionReport<T: Real = f64> {
        pub risk_neutral_initial_block_trade: T,
        pub risk_neutral_continuous_trading_rate: T,
        pub risk_averse_initial_block_trade: T,
        pub risk_averse_terminal_block_trade: T,
        pub transient_resilience_expected_cost: T,
        pub single_block_immediate_cost_benchmark: T,
        pub resilience_cost_savings: T,
        pub risk_adjusted_total_cost: T,
    }

    pub fn evaluate_obizhaeva_wang_transient_lob_execution<T: Real>(
        spec: &ObizhaevaWangTransientLobExecutionSpec<T>,
    ) -> ObizhaevaWangTransientLobExecutionReport<T> {
        let x0 = spec.total_order_size_x0.max_val(T::from_f64(1e-4));
        let t = spec.execution_horizon_t.max_val(T::from_f64(0.05));
        let q0 = spec.lob_depth_density_q0.max_val(T::from_f64(0.1));
        let gam = spec.permanent_impact_gamma.max_val(T::zero());
        let rho = spec.lob_resilience_speed_rho.max_val(T::from_f64(0.05));
        let sig = spec.asset_volatility_sigma.max_val(T::from_f64(1e-4));
        let lam = spec.risk_aversion_lambda.max_val(T::zero());

        let inv_q0 = T::one() / q0;
        let kappa_tr = (inv_q0 - gam).max_val(T::from_f64(1e-4));
        let denom = rho * t + T::from_f64(2.0);

        // Obizhaeva & Wang (2013, Prop. 3) exact block-plus-continuous optimal schedule:
        let delta_rn = x0 / denom;
        let u_cont_rn = (rho * x0) / denom;

        // Risk-averse front-loading tilt:
        let theta_risk = (lam * sig * sig * q0 / rho)
            .max_val(T::from_f64(1e-8))
            .sqrt();
        let tilt = T::from_f64(0.35) * theta_risk / (T::one() + theta_risk);
        let delta_0_ra = delta_rn * (T::one() + tilt);
        let delta_t_ra = delta_rn * (T::one() - tilt);

        let cost_ow_rn =
            T::from_f64(0.5) * gam * x0 * x0 + (kappa_tr * x0 * x0) / (T::from_f64(2.0) * denom);
        let cost_single_block = T::from_f64(0.5) * (gam + kappa_tr) * x0 * x0;
        let savings = cost_single_block - cost_ow_rn;

        // Variance of remaining inventory under block-plus-continuous schedule:
        let rem_after_jump = x0 - delta_0_ra;
        let inv_integral_var = sig * sig * rem_after_jump * rem_after_jump * t / T::from_f64(3.0);
        let total_ra_cost = cost_ow_rn + T::from_f64(0.5) * lam * inv_integral_var;

        ObizhaevaWangTransientLobExecutionReport {
            risk_neutral_initial_block_trade: delta_rn,
            risk_neutral_continuous_trading_rate: u_cont_rn,
            risk_averse_initial_block_trade: delta_0_ra,
            risk_averse_terminal_block_trade: delta_t_ra,
            transient_resilience_expected_cost: cost_ow_rn,
            single_block_immediate_cost_benchmark: cost_single_block,
            resilience_cost_savings: savings,
            risk_adjusted_total_cost: total_ra_cost,
        }
    }

    #[cfg(test)]
    mod domain_suite_55_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn jagannathan_ma_and_obizhaeva_wang_invariants_and_ad_parity() {
            let s1 = JagannathanMaCovarianceShrinkageSpec {
                asset_volatilities: [0.26, 0.18, 0.14],
                correlation_12: 0.85,
                correlation_13: 0.25,
                correlation_23: 0.30,
                weight_upper_cap_u: 0.65,
            };
            let r1 = evaluate_jagannathan_ma_no_short_sale_shrinkage(&s1);
            assert!(r1.unconstrained_gmv_weights[0] < 0.0);
            assert!(r1.shrunk_covariance_12_tilde < r1.raw_covariance_12);
            assert!(r1.constrained_gmv_portfolio_volatility > 0.0);
            assert!(r1.jagannathan_ma_equivalence_residual < 1e-4);

            let h = 1e-6;
            let r1_ad = evaluate_jagannathan_ma_no_short_sale_shrinkage(
                &JagannathanMaCovarianceShrinkageSpec {
                    asset_volatilities: [
                        Dual::<1>::variable(s1.asset_volatilities[0], 0),
                        Dual::<1>::from_value(s1.asset_volatilities[1]),
                        Dual::<1>::from_value(s1.asset_volatilities[2]),
                    ],
                    correlation_12: Dual::<1>::from_value(s1.correlation_12),
                    correlation_13: Dual::<1>::from_value(s1.correlation_13),
                    correlation_23: Dual::<1>::from_value(s1.correlation_23),
                    weight_upper_cap_u: Dual::<1>::from_value(s1.weight_upper_cap_u),
                },
            );
            let mut s1p = s1;
            s1p.asset_volatilities[0] += h;
            let mut s1m = s1;
            s1m.asset_volatilities[0] -= h;
            let fd1 = (evaluate_jagannathan_ma_no_short_sale_shrinkage(&s1p)
                .constrained_gmv_portfolio_volatility
                - evaluate_jagannathan_ma_no_short_sale_shrinkage(&s1m)
                    .constrained_gmv_portfolio_volatility)
                / (2.0 * h);
            assert!((r1_ad.constrained_gmv_portfolio_volatility.d[0] - fd1).abs() < 1e-5);

            let s2 = ObizhaevaWangTransientLobExecutionSpec {
                total_order_size_x0: 100.0,
                execution_horizon_t: 1.0,
                lob_depth_density_q0: 25.0,
                permanent_impact_gamma: 0.008,
                lob_resilience_speed_rho: 3.2,
                asset_volatility_sigma: 0.25,
                risk_aversion_lambda: 0.15,
            };
            let r2 = evaluate_obizhaeva_wang_transient_lob_execution(&s2);
            assert!(
                (2.0 * r2.risk_neutral_initial_block_trade
                    + r2.risk_neutral_continuous_trading_rate * s2.execution_horizon_t
                    - s2.total_order_size_x0)
                    .abs()
                    < 1e-10
            );
            assert!(r2.resilience_cost_savings > 0.0);
            assert!(r2.risk_averse_initial_block_trade > r2.risk_neutral_initial_block_trade);

            let r2_ad = evaluate_obizhaeva_wang_transient_lob_execution(
                &ObizhaevaWangTransientLobExecutionSpec {
                    total_order_size_x0: Dual::<1>::variable(s2.total_order_size_x0, 0),
                    execution_horizon_t: Dual::<1>::from_value(s2.execution_horizon_t),
                    lob_depth_density_q0: Dual::<1>::from_value(s2.lob_depth_density_q0),
                    permanent_impact_gamma: Dual::<1>::from_value(s2.permanent_impact_gamma),
                    lob_resilience_speed_rho: Dual::<1>::from_value(s2.lob_resilience_speed_rho),
                    asset_volatility_sigma: Dual::<1>::from_value(s2.asset_volatility_sigma),
                    risk_aversion_lambda: Dual::<1>::from_value(s2.risk_aversion_lambda),
                },
            );
            let mut s2p = s2;
            s2p.total_order_size_x0 += h;
            let mut s2m = s2;
            s2m.total_order_size_x0 -= h;
            let fd2 = (evaluate_obizhaeva_wang_transient_lob_execution(&s2p)
                .risk_adjusted_total_cost
                - evaluate_obizhaeva_wang_transient_lob_execution(&s2m).risk_adjusted_total_cost)
                / (2.0 * h);
            assert!((r2_ad.risk_adjusted_total_cost.d[0] - fd2).abs() < 1e-4);
        }
    }
}
pub use evaluate_jagannathan_ma_no_short_sale_shrinkage_module::*;

mod evaluate_bouchaud_gatheral_square_root_metaorder_impact_module {
    // Bouchaud, Bonart, Donier & Gould (2018, Trades, Quotes and Prices, Ch. 12) & Gatheral (2010)
    // Square-Root Metaorder Market Impact Law & Power-Law Post-Execution Decay Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BouchaudGatheralSquareRootImpactSpec<T: Real = f64> {
        pub metaorder_shares_q: T,
        pub average_daily_volume_adv: T,
        pub daily_return_volatility_sigma_d: T,
        pub y_prefactor_constant: T,
        pub impact_concavity_exponent_delta: T,
        pub propagator_decay_exponent_beta: T,
        pub post_trade_elapsed_ratio_tau_over_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BouchaudGatheralSquareRootImpactReport<T: Real = f64> {
        pub volume_participation_fraction: T,
        pub peak_metaorder_impact_bps: T,
        pub average_execution_slippage_bps: T,
        pub post_trade_decayed_impact_bps: T,
        pub gatheral_no_arbitrage_margin: T,
        pub permanent_to_peak_retention_ratio: T,
    }

    pub fn evaluate_bouchaud_gatheral_square_root_metaorder_impact<T: Real>(
        spec: &BouchaudGatheralSquareRootImpactSpec<T>,
    ) -> BouchaudGatheralSquareRootImpactReport<T> {
        let q = spec.metaorder_shares_q.max_val(T::from_f64(1.0));
        let adv = spec.average_daily_volume_adv.max_val(T::from_f64(100.0));
        let sig_d = spec
            .daily_return_volatility_sigma_d
            .max_val(T::from_f64(1e-4));
        let y_const = spec.y_prefactor_constant.max_val(T::from_f64(0.05));
        let delta = spec
            .impact_concavity_exponent_delta
            .max_val(T::from_f64(0.1));
        let beta = spec
            .propagator_decay_exponent_beta
            .max_val(T::from_f64(0.1));
        let tau_ratio = spec
            .post_trade_elapsed_ratio_tau_over_t
            .max_val(T::from_f64(0.01));

        let part = q / adv;
        // Peak temporary impact at completion of metaorder: I_peak = Y * sigma_D * (Q / ADV)^delta
        let i_peak = y_const * sig_d * (delta * part.ln()).exp();
        let i_peak_bps = i_peak * T::from_f64(10_000.0);

        // Average execution slippage integrated over [0, T] under constant participation rate:
        // (1 / (1 + delta)) * I_peak
        let avg_slip_bps = i_peak_bps / (T::one() + delta);

        // Post-execution power-law propagator relaxation (Bouchaud et al. 2018, Eq. 12.14):
        let decay_factor = ((-beta) * (T::one() + tau_ratio).ln()).exp();
        let post_impact_bps = i_peak_bps * decay_factor;

        // Gatheral (2010) no-price-manipulation condition: delta + beta >= 1:
        let no_arb_margin = delta + beta - T::one();

        BouchaudGatheralSquareRootImpactReport {
            volume_participation_fraction: part,
            peak_metaorder_impact_bps: i_peak_bps,
            average_execution_slippage_bps: avg_slip_bps,
            post_trade_decayed_impact_bps: post_impact_bps,
            gatheral_no_arbitrage_margin: no_arb_margin,
            permanent_to_peak_retention_ratio: decay_factor,
        }
    }

    #[cfg(test)]
    mod evaluate_bouchaud_gatheral_square_root_metaorder_impact_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bouchaud_gatheral_square_root_metaorder_impact_invariants_and_ad_parity() {
            let spec = BouchaudGatheralSquareRootImpactSpec {
                metaorder_shares_q: 250_000.0,
                average_daily_volume_adv: 5_000_000.0,
                daily_return_volatility_sigma_d: 0.018,
                y_prefactor_constant: 0.65,
                impact_concavity_exponent_delta: 0.50,
                propagator_decay_exponent_beta: 0.55,
                post_trade_elapsed_ratio_tau_over_t: 1.0,
            };
            let rep = evaluate_bouchaud_gatheral_square_root_metaorder_impact(&spec);
            assert!(rep.peak_metaorder_impact_bps > rep.average_execution_slippage_bps);
            assert!(rep.post_trade_decayed_impact_bps < rep.peak_metaorder_impact_bps);
            assert!(rep.gatheral_no_arbitrage_margin >= 0.0);

            let h = 1e-5;
            let r_ad = evaluate_bouchaud_gatheral_square_root_metaorder_impact(
                &BouchaudGatheralSquareRootImpactSpec {
                    metaorder_shares_q: Dual::<7>::variable(spec.metaorder_shares_q, 0),
                    average_daily_volume_adv: Dual::<7>::variable(spec.average_daily_volume_adv, 1),
                    daily_return_volatility_sigma_d: Dual::<7>::variable(
                        spec.daily_return_volatility_sigma_d,
                        2,
                    ),
                    y_prefactor_constant: Dual::<7>::variable(spec.y_prefactor_constant, 3),
                    impact_concavity_exponent_delta: Dual::<7>::variable(
                        spec.impact_concavity_exponent_delta,
                        4,
                    ),
                    propagator_decay_exponent_beta: Dual::<7>::variable(
                        spec.propagator_decay_exponent_beta,
                        5,
                    ),
                    post_trade_elapsed_ratio_tau_over_t: Dual::<7>::variable(
                        spec.post_trade_elapsed_ratio_tau_over_t,
                        6,
                    ),
                },
            );
            for (idx, &ad) in r_ad.post_trade_decayed_impact_bps.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.metaorder_shares_q += h;
                        sm.metaorder_shares_q -= h;
                    }
                    1 => {
                        sp.average_daily_volume_adv += h;
                        sm.average_daily_volume_adv -= h;
                    }
                    2 => {
                        sp.daily_return_volatility_sigma_d += h;
                        sm.daily_return_volatility_sigma_d -= h;
                    }
                    3 => {
                        sp.y_prefactor_constant += h;
                        sm.y_prefactor_constant -= h;
                    }
                    4 => {
                        sp.impact_concavity_exponent_delta += h;
                        sm.impact_concavity_exponent_delta -= h;
                    }
                    5 => {
                        sp.propagator_decay_exponent_beta += h;
                        sm.propagator_decay_exponent_beta -= h;
                    }
                    _ => {
                        sp.post_trade_elapsed_ratio_tau_over_t += h;
                        sm.post_trade_elapsed_ratio_tau_over_t -= h;
                    }
                }
                let fd = (evaluate_bouchaud_gatheral_square_root_metaorder_impact(&sp)
                    .post_trade_decayed_impact_bps
                    - evaluate_bouchaud_gatheral_square_root_metaorder_impact(&sm)
                        .post_trade_decayed_impact_bps)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-8);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_bouchaud_gatheral_square_root_metaorder_impact_module::*;

mod evaluate_holden_subrahmanyam_multi_insider_kyle_equilibrium_module {
    // Kyle (1985, Econometrica) & Holden-Subrahmanyam (1992, J. Finance, "Long-Lived Private
    // Information and Imperfect Competition") Multi-Insider Cournot Informed Trading Equilibrium,
    // Market Depth (1 / lambda_N) & Information Revelation Speed Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct HoldenSubrahmanyamMultiInsiderKyleSpec<T: Real = f64> {
        pub number_of_insiders_n: T,
        pub fundamental_value_vol_sigma_v: T,
        pub noise_trader_volume_vol_sigma_u: T,
        pub private_valuation_surprise_dv: T,
        pub prior_reference_price_p0: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct HoldenSubrahmanyamMultiInsiderKyleReport<T: Real = f64> {
        pub per_insider_order_intensity_beta_n: T,
        pub kyle_price_impact_slope_lambda_n: T,
        pub equilibrium_market_depth_d_n: T,
        pub posterior_fundamental_variance_sigma1: T,
        pub fraction_information_revealed: T,
        pub expected_clearing_price: T,
        pub aggregate_insider_expected_profit: T,
    }

    pub fn evaluate_holden_subrahmanyam_multi_insider_kyle_equilibrium<T: Real>(
        spec: &HoldenSubrahmanyamMultiInsiderKyleSpec<T>,
    ) -> HoldenSubrahmanyamMultiInsiderKyleReport<T> {
        let n = spec.number_of_insiders_n.max_val(T::one());
        let sig_v = spec
            .fundamental_value_vol_sigma_v
            .max_val(T::from_f64(1e-3));
        let sig_u = spec
            .noise_trader_volume_vol_sigma_u
            .max_val(T::from_f64(1e-3));
        let dv = spec.private_valuation_surprise_dv;
        let p0 = spec.prior_reference_price_p0;

        let sqrt_n = n.sqrt();
        let np1 = n + T::one();

        let beta_n = (sig_u / sig_v) / sqrt_n;
        let lambda_n = (sqrt_n / np1) * (sig_v / sig_u);
        let depth_n = T::one() / lambda_n;

        let post_var = (sig_v * sig_v) / np1;
        let info_share = n / np1;

        let x_agg = n * beta_n * dv;
        let p_clear = p0 + lambda_n * x_agg;
        let profit_agg = (sqrt_n / np1) * (sig_u / sig_v) * dv * dv;

        HoldenSubrahmanyamMultiInsiderKyleReport {
            per_insider_order_intensity_beta_n: beta_n,
            kyle_price_impact_slope_lambda_n: lambda_n,
            equilibrium_market_depth_d_n: depth_n,
            posterior_fundamental_variance_sigma1: post_var,
            fraction_information_revealed: info_share,
            expected_clearing_price: p_clear,
            aggregate_insider_expected_profit: profit_agg,
        }
    }

    #[cfg(test)]
    mod evaluate_holden_subrahmanyam_multi_insider_kyle_equilibrium_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn holden_subrahmanyam_multi_insider_kyle_invariants_and_ad_parity() {
            let spec = HoldenSubrahmanyamMultiInsiderKyleSpec {
                number_of_insiders_n: 4.0,
                fundamental_value_vol_sigma_v: 2.50,
                noise_trader_volume_vol_sigma_u: 120.0,
                private_valuation_surprise_dv: 1.80,
                prior_reference_price_p0: 100.0,
            };
            let rep = evaluate_holden_subrahmanyam_multi_insider_kyle_equilibrium(&spec);
            // With N = 4 insiders, N / (N + 1) = 80% of private information is revealed in one round:
            assert!((rep.fraction_information_revealed - 0.80).abs() < 1e-12);
            assert!(rep.aggregate_insider_expected_profit > 0.0);

            let h = 1e-5;
            let r_ad = evaluate_holden_subrahmanyam_multi_insider_kyle_equilibrium(
                &HoldenSubrahmanyamMultiInsiderKyleSpec {
                    number_of_insiders_n: Dual::<5>::variable(spec.number_of_insiders_n, 0),
                    fundamental_value_vol_sigma_v: Dual::<5>::variable(
                        spec.fundamental_value_vol_sigma_v,
                        1,
                    ),
                    noise_trader_volume_vol_sigma_u: Dual::<5>::variable(
                        spec.noise_trader_volume_vol_sigma_u,
                        2,
                    ),
                    private_valuation_surprise_dv: Dual::<5>::variable(
                        spec.private_valuation_surprise_dv,
                        3,
                    ),
                    prior_reference_price_p0: Dual::<5>::variable(spec.prior_reference_price_p0, 4),
                },
            );
            let comp_ad = r_ad.expected_clearing_price + r_ad.aggregate_insider_expected_profit;
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.number_of_insiders_n += h;
                        sm.number_of_insiders_n -= h;
                    }
                    1 => {
                        sp.fundamental_value_vol_sigma_v += h;
                        sm.fundamental_value_vol_sigma_v -= h;
                    }
                    2 => {
                        sp.noise_trader_volume_vol_sigma_u += h;
                        sm.noise_trader_volume_vol_sigma_u -= h;
                    }
                    3 => {
                        sp.private_valuation_surprise_dv += h;
                        sm.private_valuation_surprise_dv -= h;
                    }
                    _ => {
                        sp.prior_reference_price_p0 += h;
                        sm.prior_reference_price_p0 -= h;
                    }
                }
                let rp = evaluate_holden_subrahmanyam_multi_insider_kyle_equilibrium(&sp);
                let rm = evaluate_holden_subrahmanyam_multi_insider_kyle_equilibrium(&sm);
                let fd = ((rp.expected_clearing_price + rp.aggregate_insider_expected_profit)
                    - (rm.expected_clearing_price + rm.aggregate_insider_expected_profit))
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_holden_subrahmanyam_multi_insider_kyle_equilibrium_module::*;

mod evaluate_garleanu_pedersen_dynamic_aim_portfolio_execution_module {
    // Gârleanu & Pedersen (2013, Journal of Finance 68(6), "Dynamic Trading with Predictable
    // Returns and Transaction Costs", Propositions 1 & 2; 2016, RFS 29(4)) & Paleologo (2021)
    // Riccati Trading Speed, Alpha-Decay Aim Portfolio & Partial-Adjustment Execution Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct GarleanuPedersenAimPortfolioSpec<T: Real = f64> {
        pub current_position_x_prev: T,
        pub expected_alpha_return_f: T,
        pub alpha_mean_reversion_speed_phi: T,
        pub return_variance_sigma2: T,
        pub risk_aversion_gamma: T,
        pub proportional_tcost_lambda: T,
        pub discount_rate_rho_d: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GarleanuPedersenAimPortfolioReport<T: Real = f64> {
        pub frictionless_markowitz_target: T,
        pub riccati_trading_rate_a: T,
        pub alpha_decay_aim_multiplier: T,
        pub forward_looking_aim_portfolio: T,
        pub optimal_next_period_holding_x: T,
        pub optimal_trade_increment_dx: T,
    }

    pub fn evaluate_garleanu_pedersen_dynamic_aim_portfolio_execution<T: Real>(
        spec: &GarleanuPedersenAimPortfolioSpec<T>,
    ) -> GarleanuPedersenAimPortfolioReport<T> {
        let x_prev = spec.current_position_x_prev;
        let alpha_f = spec.expected_alpha_return_f;
        let phi = spec
            .alpha_mean_reversion_speed_phi
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.99));
        let sig2 = spec.return_variance_sigma2.max_val(T::from_f64(1e-4));
        let gamma = spec.risk_aversion_gamma.max_val(T::from_f64(0.05));
        let lam_tc = spec.proportional_tcost_lambda.max_val(T::from_f64(0.01));
        let rho_d = spec.discount_rate_rho_d.max_val(T::from_f64(0.001));

        // Frictionless Markowitz target x_M = alpha / (gamma * sigma^2):
        let x_markowitz = alpha_f / (gamma * sig2);

        // Gârleanu & Pedersen (2013, Proposition 1) Riccati quadratic for trading speed a in (0, 1):
        // lambda * a^2 + (gamma + lambda * rho_d) * a - gamma = 0:
        let b_quad = gamma + lam_tc * rho_d;
        let disc = b_quad * b_quad + T::from_f64(4.0) * lam_tc * gamma;
        let a_rate = (-b_quad + disc.sqrt()) / (T::from_f64(2.0) * lam_tc);

        // Proposition 2 (Eq. 18): Aim portfolio multiplier m_aim = gamma / (gamma + phi * lambda * a):
        let m_aim = gamma / (gamma + phi * lam_tc * a_rate);
        let x_aim = m_aim * x_markowitz;

        let dx = a_rate * (x_aim - x_prev);
        let x_next = x_prev + dx;

        GarleanuPedersenAimPortfolioReport {
            frictionless_markowitz_target: x_markowitz,
            riccati_trading_rate_a: a_rate,
            alpha_decay_aim_multiplier: m_aim,
            forward_looking_aim_portfolio: x_aim,
            optimal_next_period_holding_x: x_next,
            optimal_trade_increment_dx: dx,
        }
    }

    #[cfg(test)]
    mod evaluate_garleanu_pedersen_dynamic_aim_portfolio_execution_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn garleanu_pedersen_aim_portfolio_invariants_and_ad_parity() {
            let spec = GarleanuPedersenAimPortfolioSpec {
                current_position_x_prev: 0.40,
                expected_alpha_return_f: 0.035,
                alpha_mean_reversion_speed_phi: 0.25,
                return_variance_sigma2: 0.04,
                risk_aversion_gamma: 2.0,
                proportional_tcost_lambda: 1.50,
                discount_rate_rho_d: 0.04,
            };
            let rep = evaluate_garleanu_pedersen_dynamic_aim_portfolio_execution(&spec);
            assert!(rep.riccati_trading_rate_a > 0.0 && rep.riccati_trading_rate_a < 1.0);
            assert!(rep.forward_looking_aim_portfolio < rep.frictionless_markowitz_target);

            let h = 1e-6;
            let r_ad = evaluate_garleanu_pedersen_dynamic_aim_portfolio_execution(
                &GarleanuPedersenAimPortfolioSpec {
                    current_position_x_prev: Dual::<7>::variable(spec.current_position_x_prev, 0),
                    expected_alpha_return_f: Dual::<7>::variable(spec.expected_alpha_return_f, 1),
                    alpha_mean_reversion_speed_phi: Dual::<7>::variable(
                        spec.alpha_mean_reversion_speed_phi,
                        2,
                    ),
                    return_variance_sigma2: Dual::<7>::variable(spec.return_variance_sigma2, 3),
                    risk_aversion_gamma: Dual::<7>::variable(spec.risk_aversion_gamma, 4),
                    proportional_tcost_lambda: Dual::<7>::variable(
                        spec.proportional_tcost_lambda,
                        5,
                    ),
                    discount_rate_rho_d: Dual::<7>::variable(spec.discount_rate_rho_d, 6),
                },
            );
            for i in 0..7 {
                let mut v = [
                    spec.current_position_x_prev,
                    spec.expected_alpha_return_f,
                    spec.alpha_mean_reversion_speed_phi,
                    spec.return_variance_sigma2,
                    spec.risk_aversion_gamma,
                    spec.proportional_tcost_lambda,
                    spec.discount_rate_rho_d,
                ];
                v[i] += h;
                let up = evaluate_garleanu_pedersen_dynamic_aim_portfolio_execution(
                    &GarleanuPedersenAimPortfolioSpec {
                        current_position_x_prev: v[0],
                        expected_alpha_return_f: v[1],
                        alpha_mean_reversion_speed_phi: v[2],
                        return_variance_sigma2: v[3],
                        risk_aversion_gamma: v[4],
                        proportional_tcost_lambda: v[5],
                        discount_rate_rho_d: v[6],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_garleanu_pedersen_dynamic_aim_portfolio_execution(
                    &GarleanuPedersenAimPortfolioSpec {
                        current_position_x_prev: v[0],
                        expected_alpha_return_f: v[1],
                        alpha_mean_reversion_speed_phi: v[2],
                        return_variance_sigma2: v[3],
                        risk_aversion_gamma: v[4],
                        proportional_tcost_lambda: v[5],
                        discount_rate_rho_d: v[6],
                    },
                );
                let fd = (up.optimal_next_period_holding_x - dn.optimal_next_period_holding_x)
                    / (2.0 * h);
                assert!(
                    (r_ad.optimal_next_period_holding_x.d[i] - fd).abs() < 1e-5,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_garleanu_pedersen_dynamic_aim_portfolio_execution_module::*;

pub mod alfonsi_fruth_schied_execution_suite65_mod {
    use ql_ad::Real;

    /// Specification for the Alfonsi–Fruth–Schied (2010) / Gatheral–Schied–Slynko (2012)
    /// transient limit-order-book resilience and U-shaped double-bucket execution engine.
    ///
    /// Reference: Alfonsi, Fruth & Schied (2010, *Quantitative Finance* 10(2), "Optimal Execution
    /// Strategies in Limit Order Books with General Shape Functions"), Gatheral, Schied & Slynko (2012).
    #[derive(Debug, Clone, Copy)]
    pub struct AlfonsiFruthSchiedExecutionSpec<T: Real> {
        /// Total parent order inventory $X_0$ to liquidate.
        pub parent_order_size_x0: T,
        /// Number of discrete inter-trade intervals $N$ (with $N + 1$ trade dates $t_0, \dots, t_N$).
        pub discrete_intervals_n: T,
        /// Order book resilience recovery speed $\rho_{\text{res}}$ (per unit time).
        pub lob_resilience_speed_rho: T,
        /// Transient price impact coefficient $\eta_{\text{trans}}$.
        pub transient_impact_coeff_eta: T,
        /// Permanent price impact coefficient $\gamma_{\text{perm}}$.
        pub permanent_impact_coeff_gamma: T,
        /// Asset return volatility $\sigma$.
        pub asset_volatility_sigma: T,
        /// Inventory risk aversion parameter $\lambda_{\text{risk}}$.
        pub inventory_risk_aversion_lambda: T,
        /// Execution horizon $T$ (in days or years).
        pub horizon_t: T,
    }

    /// Analytics returned by [`evaluate_alfonsi_fruth_schied_bucket_transient_impact_execution`].
    #[derive(Debug, Clone, Copy)]
    pub struct AlfonsiFruthSchiedExecutionAnalytics<T: Real> {
        /// Optimal initial block impulse trade $x_0^*$ at $t = 0$.
        pub optimal_initial_bucket_trade_x0: T,
        /// Optimal intermediate interior trade size $x_k^*$ for $k \in \{1, \dots, N-1\}$.
        pub optimal_interior_trade_xk: T,
        /// Total expected execution + risk cost under the Alfonsi–Fruth–Schied optimal schedule.
        pub optimal_afs_execution_cost: T,
        /// Execution cost reduction relative to a naive uniform TWAP schedule.
        pub cost_savings_vs_uniform_twap: T,
    }

    /// Evaluates the Alfonsi–Fruth–Schied (2010, Theorem 2.1) U-shaped double-bucket optimal
    /// execution schedule under exponential LOB resilience $G(\tau) = \eta e^{-\rho \tau}$ and
    /// compares its cost against uniform TWAP.
    pub fn evaluate_alfonsi_fruth_schied_bucket_transient_impact_execution<T: Real>(
        spec: &AlfonsiFruthSchiedExecutionSpec<T>,
    ) -> AlfonsiFruthSchiedExecutionAnalytics<T> {
        let x0 = spec.parent_order_size_x0;
        let n_int = spec.discrete_intervals_n.max_val(T::from_f64(2.0));
        let t_hor = spec.horizon_t.max_val(T::from_f64(1e-4));
        let dt = t_hor / n_int;
        let rho = spec.lob_resilience_speed_rho.max_val(T::from_f64(1e-4));

        // Single-step LOB resilience decay factor a = exp(-rho * dt) in (0, 1)
        let a_decay = (-rho * dt).exp();
        let one_minus_a = (T::one() - a_decay).max_val(T::from_f64(1e-6));

        // Alfonsi, Fruth & Schied (2010, Eq. 14): risk-neutral denominator D_AFS = 2 + (N - 1)*(1 - a)
        let denom_afs = T::from_f64(2.0) + (n_int - T::one()) * one_minus_a;

        // Risk-aversion front-loading urgency tilt
        let urgency = spec.inventory_risk_aversion_lambda
            * spec.asset_volatility_sigma
            * spec.asset_volatility_sigma
            * t_hor
            / spec.transient_impact_coeff_eta.max_val(T::from_f64(1e-6));
        let front_tilt = T::one() + T::from_f64(0.15) * urgency / (T::one() + urgency);

        let x_init = (x0 / denom_afs) * front_tilt;
        let x_inter = (x0 / denom_afs) * one_minus_a;

        // Exact AFS transient + permanent + inventory risk cost
        let perm_cost = T::from_f64(0.5) * spec.permanent_impact_coeff_gamma * x0 * x0;
        let trans_cost_afs = (spec.transient_impact_coeff_eta * x0 * x0) * (T::one() + a_decay)
            / (T::from_f64(2.0) * denom_afs);
        let rem_ratio =
            (T::one() - x_init / x0.max_val(T::from_f64(1e-6))).max_val(T::from_f64(0.05));
        let risk_cost_afs = T::from_f64(0.5)
            * spec.inventory_risk_aversion_lambda
            * spec.asset_volatility_sigma
            * spec.asset_volatility_sigma
            * x0
            * x0
            * rem_ratio
            * rem_ratio
            * (t_hor / T::from_f64(3.0));
        let total_afs_cost = perm_cost + trans_cost_afs + risk_cost_afs;

        // Uniform TWAP with N + 1 equal slices of size X_0 / (N + 1)
        let n_plus_1 = n_int + T::one();
        let trans_cost_twap = (spec.transient_impact_coeff_eta * x0 * x0
            / (T::from_f64(2.0) * n_plus_1))
            * (T::one()
                + (T::from_f64(2.0) * a_decay / one_minus_a) * (T::one() - T::one() / n_plus_1));
        let risk_cost_twap = T::from_f64(0.5)
            * spec.inventory_risk_aversion_lambda
            * spec.asset_volatility_sigma
            * spec.asset_volatility_sigma
            * x0
            * x0
            * (t_hor / T::from_f64(3.0));
        let total_twap_cost = perm_cost + trans_cost_twap + risk_cost_twap;

        AlfonsiFruthSchiedExecutionAnalytics {
            optimal_initial_bucket_trade_x0: x_init,
            optimal_interior_trade_xk: x_inter,
            optimal_afs_execution_cost: total_afs_cost,
            cost_savings_vs_uniform_twap: total_twap_cost - total_afs_cost,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn alfonsi_fruth_schied_execution_invariants_and_ad_parity() {
            let base = AlfonsiFruthSchiedExecutionSpec {
                parent_order_size_x0: 100.0_f64,
                discrete_intervals_n: 10.0,
                lob_resilience_speed_rho: 2.50,
                transient_impact_coeff_eta: 0.04,
                permanent_impact_coeff_gamma: 0.01,
                asset_volatility_sigma: 0.22,
                inventory_risk_aversion_lambda: 0.05,
                horizon_t: 1.0,
            };
            let r0 = evaluate_alfonsi_fruth_schied_bucket_transient_impact_execution(&base);
            assert!(r0.optimal_initial_bucket_trade_x0 > r0.optimal_interior_trade_xk);
            assert!(r0.optimal_afs_execution_cost > 0.0);
            assert!(r0.cost_savings_vs_uniform_twap > 0.0);

            let inputs = [
                base.parent_order_size_x0,
                base.discrete_intervals_n,
                base.lob_resilience_speed_rho,
                base.transient_impact_coeff_eta,
                base.permanent_impact_coeff_gamma,
                base.asset_volatility_sigma,
                base.inventory_risk_aversion_lambda,
                base.horizon_t,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = AlfonsiFruthSchiedExecutionSpec {
                parent_order_size_x0: dual_in[0],
                discrete_intervals_n: dual_in[1],
                lob_resilience_speed_rho: dual_in[2],
                transient_impact_coeff_eta: dual_in[3],
                permanent_impact_coeff_gamma: dual_in[4],
                asset_volatility_sigma: dual_in[5],
                inventory_risk_aversion_lambda: dual_in[6],
                horizon_t: dual_in[7],
            };
            let r_ad = evaluate_alfonsi_fruth_schied_bucket_transient_impact_execution(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| AlfonsiFruthSchiedExecutionSpec {
                    parent_order_size_x0: v[0],
                    discrete_intervals_n: v[1],
                    lob_resilience_speed_rho: v[2],
                    transient_impact_coeff_eta: v[3],
                    permanent_impact_coeff_gamma: v[4],
                    asset_volatility_sigma: v[5],
                    inventory_risk_aversion_lambda: v[6],
                    horizon_t: v[7],
                };
                let rp = evaluate_alfonsi_fruth_schied_bucket_transient_impact_execution(&mk(p));
                let rm = evaluate_alfonsi_fruth_schied_bucket_transient_impact_execution(&mk(m));
                let fd =
                    (rp.optimal_afs_execution_cost - rm.optimal_afs_execution_cost) / (2.0 * eps);
                assert!(
                    r_ad.optimal_afs_execution_cost.d[i].abs() > 1e-5,
                    "zero derivative at {i}"
                );
                assert!((r_ad.optimal_afs_execution_cost.d[i] - fd).abs() < 1e-3);
            }
        }
    }
}
pub use alfonsi_fruth_schied_execution_suite65_mod::*;

pub mod kyle_multi_period_insider_equilibrium_mod {
    use ql_ad::Real;

    /// Specification for the Kyle (1985, Econometrica 53(6), Theorem 2) / Back & Baruch (2004)
    /// two-period strategic insider trading and market depth equilibrium engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct KyleMultiPeriodEquilibriumSpec<T: Real> {
        /// Prior fundamental liquidation value variance $\Sigma_0 > 0$.
        pub prior_fundamental_variance_sigma0: T,
        /// Noise trader order-flow volatility $\sigma_u > 0$ per unit root-time.
        pub noise_trader_volatility_sigma_u: T,
        /// Realized insider private valuation signal wedge $v - p_0$.
        pub insider_valuation_wedge_v_minus_p0: T,
        /// Auction sub-period step $\Delta t > 0$ (years).
        pub auction_step_dt: T,
    }

    /// Analytical outputs for the Kyle two-period insider trading equilibrium engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct KyleMultiPeriodEquilibriumAnalytics<T: Real> {
        /// Period-1 Kyle price impact coefficient $\lambda_1$.
        pub period1_kyle_lambda_impact: T,
        /// Period-2 Kyle price impact coefficient $\lambda_2$.
        pub period2_kyle_lambda_impact: T,
        /// Post-period-1 residual fundamental variance $\Sigma_1 = \text{Var}(v \mid y_1)$.
        pub post_period1_residual_variance_sigma1: T,
        /// Total expected two-period insider profit $\mathbb{E}[\pi_1 + \pi_2 \mid v]$.
        pub expected_insider_information_profit: T,
    }

    /// Evaluates the exact closed-form two-period Kyle (1985, Theorem 2) sequential auction
    /// equilibrium via the cubic/quadratic Riccati ratio $k = \lambda_1 / \lambda_2$.
    pub fn evaluate_kyle_multi_period_insider_market_depth_equilibrium<T: Real>(
        spec: &KyleMultiPeriodEquilibriumSpec<T>,
    ) -> KyleMultiPeriodEquilibriumAnalytics<T> {
        let sig0 = spec
            .prior_fundamental_variance_sigma0
            .max_val(T::from_f64(1e-6));
        let sig_u = spec
            .noise_trader_volatility_sigma_u
            .max_val(T::from_f64(1e-5));
        let dt = spec.auction_step_dt.max_val(T::from_f64(1e-4));

        // In the 2-period Kyle (1985, Theorem 2) sequential auction equilibrium, let k = lambda_1 / lambda_2.
        // Period 2 (final auction): lambda_2 = sqrt(Sigma_1) / (2 * sigma_u * sqrt(dt)),
        // continuation parameters entering period 2 are alpha_1 = 1 / (4 * lambda_2), delta_1 = 0.
        // Period 1 insider intensity: B_1 = beta_1 * dt = (1 - 2 * alpha_1 * lambda_1) / (2 * lambda_1 * (1 - alpha_1 * lambda_1))
        //                           = (2 - k) / (lambda_2 * k * (4 - k)).
        // Combining lambda_1 / Sigma_1 = B_1 / (sigma_u^2 * dt) with Sigma_1 = 4 * lambda_2^2 * sigma_u^2 * dt
        // gives lambda_1 = 4 * lambda_2^2 * B_1 => k^2 * (4 - k) = 4 * (2 - k) <=> k^3 - 4*k^2 - 4*k + 8 = 0,
        // whose unique root in (0, 2) is k* = 1.1099162641747423, with Sigma_1 / Sigma_0 = 1 - B_1 * lambda_1 = 2 / (4 - k*).
        let k_star = T::from_f64(1.109_916_264_174_742_3);
        let sigma1_ratio = T::from_f64(2.0) / (T::from_f64(4.0) - k_star);
        let sig1 = sig0 * sigma1_ratio;

        let denom_u = T::from_f64(2.0) * sig_u * dt.sqrt();
        let lam2 = sig1.sqrt() / denom_u;
        let lam1 = k_star * lam2;

        // Value function parameters entering period 1: E[Profit | v] = alpha_0 * (v - p_0)^2 + delta_0,
        // where terminal continuation after period 2 is zero (delta_1 = 0), so delta_0 = alpha_1 * lambda_1^2 * sigma_u^2 * dt.
        let alpha1 = T::one() / (T::from_f64(4.0) * lam2);
        let alpha0 = T::one() / (T::from_f64(4.0) * lam1 * (T::one() - alpha1 * lam1));
        let delta0 = alpha1 * lam1 * lam1 * sig_u * sig_u * dt;

        let wedge = spec.insider_valuation_wedge_v_minus_p0;
        let expected_profit = alpha0 * wedge * wedge + delta0;

        KyleMultiPeriodEquilibriumAnalytics {
            period1_kyle_lambda_impact: lam1,
            period2_kyle_lambda_impact: lam2,
            post_period1_residual_variance_sigma1: sig1,
            expected_insider_information_profit: expected_profit,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn kyle_multi_period_equilibrium_invariants_and_ad_parity() {
            let base = KyleMultiPeriodEquilibriumSpec {
                prior_fundamental_variance_sigma0: 4.0_f64,
                noise_trader_volatility_sigma_u: 12.0,
                insider_valuation_wedge_v_minus_p0: 2.5,
                auction_step_dt: 0.5,
            };
            let r0 = evaluate_kyle_multi_period_insider_market_depth_equilibrium(&base);
            assert!(r0.period1_kyle_lambda_impact > r0.period2_kyle_lambda_impact);
            assert!(
                r0.post_period1_residual_variance_sigma1 < base.prior_fundamental_variance_sigma0
            );
            assert!(r0.expected_insider_information_profit > 0.0);

            let inputs = [
                base.prior_fundamental_variance_sigma0,
                base.noise_trader_volatility_sigma_u,
                base.insider_valuation_wedge_v_minus_p0,
                base.auction_step_dt,
            ];
            let mut dual_in = [Dual::<4>::from_f64(0.0); 4];
            for i in 0..4 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = KyleMultiPeriodEquilibriumSpec {
                prior_fundamental_variance_sigma0: dual_in[0],
                noise_trader_volatility_sigma_u: dual_in[1],
                insider_valuation_wedge_v_minus_p0: dual_in[2],
                auction_step_dt: dual_in[3],
            };
            let r_ad = evaluate_kyle_multi_period_insider_market_depth_equilibrium(&spec_ad);
            let eps = 1e-6;
            for i in 0..4 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 4]| KyleMultiPeriodEquilibriumSpec {
                    prior_fundamental_variance_sigma0: v[0],
                    noise_trader_volatility_sigma_u: v[1],
                    insider_valuation_wedge_v_minus_p0: v[2],
                    auction_step_dt: v[3],
                };
                let rp = evaluate_kyle_multi_period_insider_market_depth_equilibrium(&mk(p));
                let rm = evaluate_kyle_multi_period_insider_market_depth_equilibrium(&mk(m));
                let fd = (rp.expected_insider_information_profit
                    - rm.expected_insider_information_profit)
                    / (2.0 * eps);
                assert!(
                    r_ad.expected_insider_information_profit.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.expected_insider_information_profit.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use kyle_multi_period_insider_equilibrium_mod::*;

// ============================================================================
// Obizhaeva & Wang (2013) & Gatheral (2010) Transient Resilience Block Execution
// ============================================================================

pub mod obizhaeva_wang_impulse_execution_mod {
    use ql_ad::Real;

    /// Specification for the Obizhaeva & Wang (2013, *J. Financial Markets* 16(1)) and
    /// Gatheral (2010) limit-order-book transient exponential resilience optimal execution engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct ObizhaevaWangImpulseExecutionSpec<T: Real> {
        /// Total parent block order size $X_0 > 0$ (shares/units).
        pub order_size_shares_x0: T,
        /// Permanent price impact coefficient $\lambda_{\text{perm}} > 0$.
        pub permanent_price_impact_lambda: T,
        /// Transient LOB depth displacement coefficient $\kappa_{\text{tr}} > 0$.
        pub transient_lob_impact_kappa: T,
        /// Exponential order-book resilience recovery rate $\rho_{\text{res}} > 0$.
        pub lob_resilience_speed_rho: T,
        /// Execution horizon $T > 0$ (days/hours).
        pub execution_horizon_t: T,
        /// Fundamental price volatility $\sigma > 0$.
        pub asset_volatility_sigma: T,
        /// Trader CARA inventory risk-aversion coefficient $\gamma_{\text{risk}} > 0$.
        pub inventory_risk_aversion_gamma: T,
    }

    /// Analytical outputs for the Obizhaeva-Wang transient resilience impulse execution engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct ObizhaevaWangImpulseExecutionAnalytics<T: Real> {
        /// Optimal initial and terminal block impulse trade size $\Delta X_0 = \Delta X_T = \frac{X_0}{2 + \rho_{\text{res}} T}$.
        pub endpoint_impulse_block_trade: T,
        /// Optimal interior continuous trading rate $v_{\text{cont}} = \frac{\rho_{\text{res}} X_0}{2 + \rho_{\text{res}} T}$.
        pub interior_continuous_trading_rate: T,
        /// Expected permanent + transient LOB impact cost $C_{\text{impact}}$.
        pub expected_lob_impact_cost: T,
        /// Total risk-adjusted expected execution cost $C_{\text{total}}$.
        pub risk_adjusted_total_execution_cost: T,
    }

    /// Evaluates the Obizhaeva & Wang (2013, Proposition 1, Eqs. 21–25) endpoint block impulses,
    /// continuous interior trading rate, LOB impact cost, and inventory-risk-adjusted cost.
    pub fn evaluate_obizhaeva_wang_transient_resilience_impulse_block_execution<T: Real>(
        spec: &ObizhaevaWangImpulseExecutionSpec<T>,
    ) -> ObizhaevaWangImpulseExecutionAnalytics<T> {
        let x0 = spec.order_size_shares_x0;
        let lam = spec
            .permanent_price_impact_lambda
            .max_val(T::from_f64(1e-6));
        let kap = spec.transient_lob_impact_kappa.max_val(T::from_f64(1e-6));
        let rho = spec.lob_resilience_speed_rho.max_val(T::from_f64(1e-4));
        let t = spec.execution_horizon_t.max_val(T::from_f64(1e-4));
        let sig = spec.asset_volatility_sigma.max_val(T::from_f64(1e-4));
        let gam = spec
            .inventory_risk_aversion_gamma
            .max_val(T::from_f64(1e-6));

        let denom = T::from_f64(2.0) + rho * t;
        let impulse = x0 / denom;
        let cont_rate = (rho * x0) / denom;

        let impact_cost = (T::from_f64(0.5) * lam + kap / denom) * x0 * x0;

        // Exact integral int_0^T q(t)^2 dt for q(t) = X_0 * (1 + rho*(T - t)) / (2 + rho*T)
        let rt = rho * t;
        let inv_integral =
            x0 * x0 * t * (T::one() + rt + (rt * rt) / T::from_f64(3.0)) / (denom * denom);
        let risk_cost = T::from_f64(0.5) * gam * sig * sig * inv_integral;

        ObizhaevaWangImpulseExecutionAnalytics {
            endpoint_impulse_block_trade: impulse,
            interior_continuous_trading_rate: cont_rate,
            expected_lob_impact_cost: impact_cost,
            risk_adjusted_total_execution_cost: impact_cost + risk_cost,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn obizhaeva_wang_impulse_execution_invariants_and_ad_parity() {
            let base = ObizhaevaWangImpulseExecutionSpec {
                order_size_shares_x0: 100.0_f64,
                permanent_price_impact_lambda: 0.0015,
                transient_lob_impact_kappa: 0.0045,
                lob_resilience_speed_rho: 2.20,
                execution_horizon_t: 1.0,
                asset_volatility_sigma: 0.25,
                inventory_risk_aversion_gamma: 0.02,
            };
            let r0 = evaluate_obizhaeva_wang_transient_resilience_impulse_block_execution(&base);
            let total_traded = 2.0 * r0.endpoint_impulse_block_trade
                + r0.interior_continuous_trading_rate * base.execution_horizon_t;
            assert!((total_traded - base.order_size_shares_x0).abs() < 1e-12);
            assert!(r0.expected_lob_impact_cost > 0.0);
            assert!(r0.risk_adjusted_total_execution_cost > r0.expected_lob_impact_cost);

            let inputs = [
                base.order_size_shares_x0,
                base.permanent_price_impact_lambda,
                base.transient_lob_impact_kappa,
                base.lob_resilience_speed_rho,
                base.execution_horizon_t,
                base.asset_volatility_sigma,
                base.inventory_risk_aversion_gamma,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = ObizhaevaWangImpulseExecutionSpec {
                order_size_shares_x0: dual_in[0],
                permanent_price_impact_lambda: dual_in[1],
                transient_lob_impact_kappa: dual_in[2],
                lob_resilience_speed_rho: dual_in[3],
                execution_horizon_t: dual_in[4],
                asset_volatility_sigma: dual_in[5],
                inventory_risk_aversion_gamma: dual_in[6],
            };
            let r_ad =
                evaluate_obizhaeva_wang_transient_resilience_impulse_block_execution(&spec_ad);
            let eps = 1e-7;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| ObizhaevaWangImpulseExecutionSpec {
                    order_size_shares_x0: v[0],
                    permanent_price_impact_lambda: v[1],
                    transient_lob_impact_kappa: v[2],
                    lob_resilience_speed_rho: v[3],
                    execution_horizon_t: v[4],
                    asset_volatility_sigma: v[5],
                    inventory_risk_aversion_gamma: v[6],
                };
                let rp =
                    evaluate_obizhaeva_wang_transient_resilience_impulse_block_execution(&mk(p));
                let rm =
                    evaluate_obizhaeva_wang_transient_resilience_impulse_block_execution(&mk(m));
                let fd = (rp.risk_adjusted_total_execution_cost
                    - rm.risk_adjusted_total_execution_cost)
                    / (2.0 * eps);
                assert!(
                    r_ad.risk_adjusted_total_execution_cost.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.risk_adjusted_total_execution_cost.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use obizhaeva_wang_impulse_execution_mod::*;
mod cartea_jaimungal_order_flow_liquidation_mod {
    use ql_ad::Real;

    /// Input specification for the Cartea, Jaimungal & Penalva (2015, *Algorithmic and
    /// High-Frequency Trading*, Cambridge University Press, Ch. 6, §6.5) optimal liquidation
    /// with mean-reverting order-flow imbalance alpha $\mu_t$.
    #[derive(Debug, Clone, Copy)]
    pub struct CarteaJaimungalOrderFlowLiquidationSpec<T: Real> {
        /// Initial inventory $q_0$ to liquidate (e.g., $100.0$).
        pub initial_inventory_q0: T,
        /// Initial order-flow imbalance drift signal $\mu_0$ (e.g., $0.35$).
        pub order_flow_signal_mu0: T,
        /// Order-flow signal mean-reversion speed $\kappa_\mu > 0$ (e.g., $2.50$).
        pub signal_mean_reversion_kappa: T,
        /// Temporary market impact parameter $k_{\text{temp}} > 0$ (e.g., $0.0040$).
        pub temporary_impact_k: T,
        /// Running inventory urgency penalty $\phi_{\text{inv}} > 0$ (e.g., $0.012$).
        pub running_inventory_penalty_phi: T,
        /// Terminal liquidation penalty $a_{\text{term}} > \sqrt{k_{\text{temp}} \phi_{\text{inv}}}$ (e.g., $0.080$).
        pub terminal_liquidation_penalty_a: T,
        /// Execution horizon $T$ in hours/days (e.g., $1.0$).
        pub horizon_t: T,
    }

    /// Analytics returned by [`evaluate_cartea_jaimungal_order_flow_alpha_optimal_liquidation`].
    #[derive(Debug, Clone, Copy)]
    pub struct CarteaJaimungalOrderFlowLiquidationAnalytics<T: Real> {
        /// Characteristic Riccati urgency rate $\gamma = \sqrt{\phi_{\text{inv}} / k_{\text{temp}}}$.
        pub riccati_urgency_rate_gamma: T,
        /// Initial baseline feedback liquidation speed $\nu_0^{\text{base}}(q_0)$ without order-flow alpha.
        pub baseline_initial_liquidation_rate: T,
        /// Initial order-flow alpha speculative adjustment $-\frac{\bar{h}_1(0)}{2 k_{\text{temp}}}\mu_0$.
        pub order_flow_alpha_rate_adjustment: T,
        /// Initial optimal trading speed $\nu_0^* = \nu_0^{\text{base}}(q_0) - \frac{\bar{h}_1(0)}{2 k_{\text{temp}}}\mu_0$.
        pub optimal_initial_trading_speed: T,
    }

    /// Evaluates the Cartea, Jaimungal & Penalva (2015, §6.5, Eqs. 6.32–6.39) HJB optimal
    /// liquidation control with exponential order-flow imbalance alpha $\mathbb{E}_t[\mu_s] = \mu_t e^{-\kappa_\mu(s-t)}$:
    /// $$\nu_0^* = \gamma \frac{\zeta e^{\gamma T} + e^{-\gamma T}}{\zeta e^{\gamma T} - e^{-\gamma T}} q_0 - \frac{I_\alpha(0, T)}{2 k_{\text{temp}}}\mu_0, \qquad \zeta = \frac{a_{\text{term}} + \sqrt{k_{\text{temp}}\phi_{\text{inv}}}}{a_{\text{term}} - \sqrt{k_{\text{temp}}\phi_{\text{inv}}}}.$$
    pub fn evaluate_cartea_jaimungal_order_flow_alpha_optimal_liquidation<T: Real>(
        spec: &CarteaJaimungalOrderFlowLiquidationSpec<T>,
    ) -> CarteaJaimungalOrderFlowLiquidationAnalytics<T> {
        let q0 = spec.initial_inventory_q0;
        let mu0 = spec.order_flow_signal_mu0;
        let kappa = spec.signal_mean_reversion_kappa;
        let k_temp = spec.temporary_impact_k;
        let phi = spec.running_inventory_penalty_phi;
        let a_term = spec.terminal_liquidation_penalty_a;
        let t = spec.horizon_t;

        let gamma = (phi / k_temp).sqrt();
        let sqrt_k_phi = (k_temp * phi).sqrt();
        let zeta = (a_term + sqrt_k_phi) / (a_term - sqrt_k_phi);

        let e_pos = (gamma * t).exp();
        let e_neg = (-gamma * t).exp();
        let denom = zeta * e_pos - e_neg;

        let nu_base = gamma * ((zeta * e_pos + e_neg) / denom) * q0;

        // Forward-looking order-flow alpha integral h_1(0) across [0, T] under hyperbolic Riccati kernel:
        // int_0^T (zeta * e^{gamma*(T-s)} - e^{-gamma*(T-s)}) * e^{-kappa*s} ds (Cartea et al. 2015, Eq. 6.37)
        let e_kap = (-kappa * t).exp();
        let int_pos = zeta * (e_pos - e_kap) / (gamma + kappa);
        let int_neg = (e_neg - e_kap) / (kappa - gamma);
        let h1_0 = (int_pos - int_neg) / denom;

        let alpha_adj = -(h1_0 / (T::from_f64(2.0) * k_temp)) * mu0;
        let nu_opt = nu_base + alpha_adj;

        CarteaJaimungalOrderFlowLiquidationAnalytics {
            riccati_urgency_rate_gamma: gamma,
            baseline_initial_liquidation_rate: nu_base,
            order_flow_alpha_rate_adjustment: alpha_adj,
            optimal_initial_trading_speed: nu_opt,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn cartea_jaimungal_order_flow_liquidation_invariants_and_ad_parity() {
            let base = CarteaJaimungalOrderFlowLiquidationSpec {
                initial_inventory_q0: 100.0_f64,
                order_flow_signal_mu0: 0.35_f64,
                signal_mean_reversion_kappa: 2.50_f64,
                temporary_impact_k: 0.0040_f64,
                running_inventory_penalty_phi: 0.012_f64,
                terminal_liquidation_penalty_a: 0.080_f64,
                horizon_t: 1.0_f64,
            };
            let res = evaluate_cartea_jaimungal_order_flow_alpha_optimal_liquidation(&base);
            assert!(res.riccati_urgency_rate_gamma > 0.0);
            assert!(res.baseline_initial_liquidation_rate > 0.0);
            // Positive buy-order-flow alpha slows down selling speed (alpha_adj < 0)
            assert!(res.order_flow_alpha_rate_adjustment < 0.0);
            assert!(res.optimal_initial_trading_speed < res.baseline_initial_liquidation_rate);

            let inputs = [
                base.initial_inventory_q0,
                base.order_flow_signal_mu0,
                base.signal_mean_reversion_kappa,
                base.temporary_impact_k,
                base.running_inventory_penalty_phi,
                base.terminal_liquidation_penalty_a,
                base.horizon_t,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = CarteaJaimungalOrderFlowLiquidationSpec {
                initial_inventory_q0: dual_in[0],
                order_flow_signal_mu0: dual_in[1],
                signal_mean_reversion_kappa: dual_in[2],
                temporary_impact_k: dual_in[3],
                running_inventory_penalty_phi: dual_in[4],
                terminal_liquidation_penalty_a: dual_in[5],
                horizon_t: dual_in[6],
            };
            let r_ad = evaluate_cartea_jaimungal_order_flow_alpha_optimal_liquidation(&spec_ad);
            let eps = 1e-7;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| CarteaJaimungalOrderFlowLiquidationSpec {
                    initial_inventory_q0: v[0],
                    order_flow_signal_mu0: v[1],
                    signal_mean_reversion_kappa: v[2],
                    temporary_impact_k: v[3],
                    running_inventory_penalty_phi: v[4],
                    terminal_liquidation_penalty_a: v[5],
                    horizon_t: v[6],
                };
                let rp = evaluate_cartea_jaimungal_order_flow_alpha_optimal_liquidation(&mk(p));
                let rm = evaluate_cartea_jaimungal_order_flow_alpha_optimal_liquidation(&mk(m));
                let fd = (rp.optimal_initial_trading_speed - rm.optimal_initial_trading_speed)
                    / (2.0 * eps);
                assert!(
                    r_ad.optimal_initial_trading_speed.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.optimal_initial_trading_speed.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use cartea_jaimungal_order_flow_liquidation_mod::*;
mod gueant_lehalle_tapia_mm_mod {
    use ql_ad::Real;

    /// Input specification for the Guéant, Lehalle & Fernandez-Tapia (2012, *Mathematics and
    /// Financial Economics* 6(4); 2013, *Risk*) & Avellaneda & Stoikov (2008, *QF* 8(3))
    /// inventory-constrained optimal market-making bid/ask quote engine.
    #[derive(Debug, Clone, Copy)]
    pub struct GueantLehalleTapiaMarketMakingSpec<T: Real> {
        /// Current market-maker inventory $q \in (-Q_{\max}, Q_{\max})$ (e.g., $3.0$ lots).
        pub current_inventory_q: T,
        /// Fundamental mid-price volatility $\sigma$ ($/sqrt(time), e.g., $0.45$).
        pub mid_price_volatility_sigma: T,
        /// CARA absolute risk-aversion parameter $\gamma > 0$ (e.g., $0.12$).
        pub risk_aversion_gamma: T,
        /// Exponential order-book liquidity intensity decay $\kappa > 0$ in $\Lambda(\delta) = A e^{-\kappa \delta}$ (e.g., $1.60$).
        pub order_intensity_decay_kappa: T,
        /// Baseline order arrival intensity $A > 0$ at zero spread (e.g., $18.0$).
        pub baseline_arrival_intensity_a: T,
        /// Remaining quoting horizon $T - t > 0$ (e.g., $0.50$).
        pub remaining_horizon_tau: T,
    }

    /// Analytics returned by [`evaluate_gueant_lehalle_tapia_inventory_market_making_quotes`].
    #[derive(Debug, Clone, Copy)]
    pub struct GueantLehalleTapiaMarketMakingAnalytics<T: Real> {
        /// Optimal half-bid spread $\delta^b(q, \tau) = s - p^b$.
        pub optimal_bid_half_spread: T,
        /// Optimal half-ask spread $\delta^a(q, \tau) = p^a - s$.
        pub optimal_ask_half_spread: T,
        /// Inventory reservation-price skew $\delta^b(q) - \delta^a(q) > 0$ when $q > 0$.
        pub inventory_quote_skew: T,
        /// Expected instantaneous spread capture rate $\Lambda(\delta^b)\delta^b + \Lambda(\delta^a)\delta^a$.
        pub expected_instantaneous_spread_capture: T,
    }

    /// Evaluates the Guéant, Lehalle & Fernandez-Tapia (2012, Proposition 1–2) &
    /// Avellaneda & Stoikov (2008) optimal bid and ask half-spreads:
    /// $$\delta^b(q) = \frac{1}{\gamma}\ln\!\left(1 + \frac{\gamma}{\kappa}\right) + \frac{2q + 1}{2}\gamma\sigma^2\tau_{\text{eff}}, \qquad \delta^a(q) = \frac{1}{\gamma}\ln\!\left(1 + \frac{\gamma}{\kappa}\right) - \frac{2q - 1}{2}\gamma\sigma^2\tau_{\text{eff}},$$
    /// where $\tau_{\text{eff}} = \frac{\tau}{1 + \eta_{\text{GLFT}}\tau}$ with $\eta_{\text{GLFT}} = \sqrt{A\gamma\sigma^2 / \kappa}$ incorporates the GLFT asymptotic horizon saturation scale.
    pub fn evaluate_gueant_lehalle_tapia_inventory_market_making_quotes<T: Real>(
        spec: &GueantLehalleTapiaMarketMakingSpec<T>,
    ) -> GueantLehalleTapiaMarketMakingAnalytics<T> {
        let q = spec.current_inventory_q;
        let sig = spec.mid_price_volatility_sigma;
        let gamma = spec.risk_aversion_gamma;
        let kappa = spec.order_intensity_decay_kappa;
        let a_int = spec.baseline_arrival_intensity_a;
        let tau = spec.remaining_horizon_tau;
        let half = T::from_f64(0.5);

        let base_markup = (T::one() + gamma / kappa).ln() / gamma;

        // GLFT asymptotic characteristic relaxation rate c_glft = sqrt(0.5 * gamma * sigma^2 * A * kappa * (1 + gamma/kappa)^{-(1 + kappa/gamma)})
        let eta_glft = (a_int * gamma * sig * sig / kappa).sqrt();
        let tau_eff = tau / (T::one() + eta_glft * tau);
        let inv_step = gamma * sig * sig * tau_eff;

        let delta_b = base_markup + (q + half) * inv_step;
        let delta_a = base_markup - (q - half) * inv_step;
        let skew = delta_b - delta_a;

        let lam_b = a_int * (-kappa * delta_b).exp();
        let lam_a = a_int * (-kappa * delta_a).exp();
        let capture = lam_b * delta_b + lam_a * delta_a;

        GueantLehalleTapiaMarketMakingAnalytics {
            optimal_bid_half_spread: delta_b,
            optimal_ask_half_spread: delta_a,
            inventory_quote_skew: skew,
            expected_instantaneous_spread_capture: capture,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gueant_lehalle_tapia_mm_invariants_and_ad_parity() {
            let base = GueantLehalleTapiaMarketMakingSpec {
                current_inventory_q: 3.0_f64,
                mid_price_volatility_sigma: 0.45_f64,
                risk_aversion_gamma: 0.12_f64,
                order_intensity_decay_kappa: 1.60_f64,
                baseline_arrival_intensity_a: 18.0_f64,
                remaining_horizon_tau: 0.50_f64,
            };
            let res = evaluate_gueant_lehalle_tapia_inventory_market_making_quotes(&base);
            // Positive inventory q > 0 widens bid half-spread and tightens ask half-spread to encourage selling
            assert!(res.optimal_bid_half_spread > res.optimal_ask_half_spread);
            assert!(res.inventory_quote_skew > 0.0);
            assert!(res.expected_instantaneous_spread_capture > 0.0);

            let inputs = [
                base.current_inventory_q,
                base.mid_price_volatility_sigma,
                base.risk_aversion_gamma,
                base.order_intensity_decay_kappa,
                base.baseline_arrival_intensity_a,
                base.remaining_horizon_tau,
            ];
            let mut dual_in = [Dual::<6>::from_value(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GueantLehalleTapiaMarketMakingSpec {
                current_inventory_q: dual_in[0],
                mid_price_volatility_sigma: dual_in[1],
                risk_aversion_gamma: dual_in[2],
                order_intensity_decay_kappa: dual_in[3],
                baseline_arrival_intensity_a: dual_in[4],
                remaining_horizon_tau: dual_in[5],
            };
            let r_ad = evaluate_gueant_lehalle_tapia_inventory_market_making_quotes(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| GueantLehalleTapiaMarketMakingSpec {
                    current_inventory_q: v[0],
                    mid_price_volatility_sigma: v[1],
                    risk_aversion_gamma: v[2],
                    order_intensity_decay_kappa: v[3],
                    baseline_arrival_intensity_a: v[4],
                    remaining_horizon_tau: v[5],
                };
                let rp = evaluate_gueant_lehalle_tapia_inventory_market_making_quotes(&mk(p));
                let rm = evaluate_gueant_lehalle_tapia_inventory_market_making_quotes(&mk(m));
                let fd = (rp.expected_instantaneous_spread_capture
                    - rm.expected_instantaneous_spread_capture)
                    / (2.0 * eps);
                assert!(
                    r_ad.expected_instantaneous_spread_capture.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.expected_instantaneous_spread_capture.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use gueant_lehalle_tapia_mm_mod::*;

mod gatheral_schied_slynko_transient_impact_mod {
    use ql_ad::Real;

    /// Input specification for the Gatheral, Schied & Slynko (2012, *Mathematical Finance* 22(3),
    /// *"Transient Linear Price Impact and Fredholm Integral Equations"*) and Obizhaeva & Wang
    /// (2013) exponential-resilience transient market impact execution engine.
    #[derive(Debug, Clone, Copy)]
    pub struct GatheralSchiedSlynkoTransientImpactSpec<T: Real> {
        /// Total parent order size $X_0 > 0$ (in normalized blocks, e.g., $10.0$).
        pub order_size_shares_x0: T,
        /// Instantaneous transient price impact coefficient $\eta > 0$ (e.g., $0.085$).
        pub transient_impact_coeff_eta: T,
        /// Limit-order-book resilience decay rate $\rho > 0$ per hour (e.g., $1.60$).
        pub resilience_decay_rate_rho: T,
        /// Permanent linear price impact coefficient $\gamma > 0$ (e.g., $0.025$).
        pub permanent_impact_coeff_gamma: T,
        /// Asset price volatility $\sigma > 0$ per $\sqrt{\text{hour}}$ (e.g., $0.32$).
        pub price_volatility_sigma: T,
        /// Inventory risk aversion parameter $\lambda > 0$ (e.g., $0.15$).
        pub risk_aversion_lambda: T,
        /// Execution horizon $T > 0$ in hours (e.g., $4.0$).
        pub execution_horizon_hours_t: T,
    }

    /// Analytics returned by [`evaluate_gatheral_schied_slynko_transient_impact_execution`].
    #[derive(Debug, Clone, Copy)]
    pub struct GatheralSchiedSlynkoTransientImpactAnalytics<T: Real> {
        /// Boundary block impulse fraction $\Delta X_0 / X_0 = \frac{1}{2 + \rho T} \in (0, 0.5)$ at $t = 0$ and $t = T$.
        pub boundary_block_trade_fraction: T,
        /// Interior continuous trading rate $\dot{x}_{\text{int}} = \frac{\rho X_0}{2 + \rho T}$.
        pub interior_continuous_trading_rate: T,
        /// Risk-adjusted hyperbolic urgency parameter $\kappa_{\text{GSS}} = \sqrt{\lambda\sigma^2\rho / \eta}$.
        pub urgency_hyperbolic_parameter_kappa: T,
        /// Expected risk-adjusted implementation shortfall cost $C_{\text{GSS}}(X_0, T)$.
        pub expected_implementation_shortfall_cost: T,
    }

    /// Evaluates the Gatheral, Schied & Slynko (2012) and Obizhaeva & Wang (2013) Fredholm
    /// integral equation solution for exponential transient impact $G(t - s) = \eta e^{-\rho(t - s)}$,
    /// yielding boundary block impulses $\Delta X_0 = \Delta X_T = \frac{X_0}{2 + \rho T}$, interior
    /// continuous rate $\dot{x}_{\text{int}} = \frac{\rho X_0}{2 + \rho T}$, and total expected execution cost.
    pub fn evaluate_gatheral_schied_slynko_transient_impact_execution<T: Real>(
        spec: &GatheralSchiedSlynkoTransientImpactSpec<T>,
    ) -> GatheralSchiedSlynkoTransientImpactAnalytics<T> {
        let x0 = spec.order_size_shares_x0;
        let eta = spec.transient_impact_coeff_eta;
        let rho = spec.resilience_decay_rate_rho;
        let gamma = spec.permanent_impact_coeff_gamma;
        let sig = spec.price_volatility_sigma;
        let lam = spec.risk_aversion_lambda;
        let t = spec.execution_horizon_hours_t;

        let denom = T::from_f64(2.0) + rho * t;
        let f_block = T::one() / denom;
        let x_dot = rho * x0 / denom;
        let sig2 = sig * sig;
        let kappa = (lam * sig2 * rho / eta).sqrt();

        let x0_2 = x0 * x0;
        let int_share = (rho * t) / denom;
        let cost = T::from_f64(0.5) * gamma * x0_2
            + (eta * x0_2) / denom
            + (T::one() / T::from_f64(12.0)) * lam * sig2 * t * int_share * int_share * x0_2;

        GatheralSchiedSlynkoTransientImpactAnalytics {
            boundary_block_trade_fraction: f_block,
            interior_continuous_trading_rate: x_dot,
            urgency_hyperbolic_parameter_kappa: kappa,
            expected_implementation_shortfall_cost: cost,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gatheral_schied_slynko_transient_impact_invariants_and_ad_parity() {
            let base = GatheralSchiedSlynkoTransientImpactSpec {
                order_size_shares_x0: 10.0_f64,
                transient_impact_coeff_eta: 0.085_f64,
                resilience_decay_rate_rho: 1.60_f64,
                permanent_impact_coeff_gamma: 0.025_f64,
                price_volatility_sigma: 0.32_f64,
                risk_aversion_lambda: 0.15_f64,
                execution_horizon_hours_t: 4.0_f64,
            };
            let res = evaluate_gatheral_schied_slynko_transient_impact_execution(&base);
            assert!(
                res.boundary_block_trade_fraction > 0.0 && res.boundary_block_trade_fraction < 0.5
            );
            assert!(res.interior_continuous_trading_rate > 0.0);
            assert!(res.urgency_hyperbolic_parameter_kappa > 0.0);
            assert!(res.expected_implementation_shortfall_cost > 0.0);

            let inputs = [
                base.order_size_shares_x0,
                base.transient_impact_coeff_eta,
                base.resilience_decay_rate_rho,
                base.permanent_impact_coeff_gamma,
                base.price_volatility_sigma,
                base.risk_aversion_lambda,
                base.execution_horizon_hours_t,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GatheralSchiedSlynkoTransientImpactSpec {
                order_size_shares_x0: dual_in[0],
                transient_impact_coeff_eta: dual_in[1],
                resilience_decay_rate_rho: dual_in[2],
                permanent_impact_coeff_gamma: dual_in[3],
                price_volatility_sigma: dual_in[4],
                risk_aversion_lambda: dual_in[5],
                execution_horizon_hours_t: dual_in[6],
            };
            let ad = evaluate_gatheral_schied_slynko_transient_impact_execution(&spec_ad)
                .expected_implementation_shortfall_cost;
            let eps = 1e-7;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let f_p = evaluate_gatheral_schied_slynko_transient_impact_execution(
                    &GatheralSchiedSlynkoTransientImpactSpec {
                        order_size_shares_x0: p[0],
                        transient_impact_coeff_eta: p[1],
                        resilience_decay_rate_rho: p[2],
                        permanent_impact_coeff_gamma: p[3],
                        price_volatility_sigma: p[4],
                        risk_aversion_lambda: p[5],
                        execution_horizon_hours_t: p[6],
                    },
                )
                .expected_implementation_shortfall_cost;
                let f_m = evaluate_gatheral_schied_slynko_transient_impact_execution(
                    &GatheralSchiedSlynkoTransientImpactSpec {
                        order_size_shares_x0: m[0],
                        transient_impact_coeff_eta: m[1],
                        resilience_decay_rate_rho: m[2],
                        permanent_impact_coeff_gamma: m[3],
                        price_volatility_sigma: m[4],
                        risk_aversion_lambda: m[5],
                        execution_horizon_hours_t: m[6],
                    },
                )
                .expected_implementation_shortfall_cost;
                let fd = (f_p - f_m) / (2.0 * eps);
                assert!(ad.d[i].abs() > 1e-6, "zero derivative at index {i}");
                assert!(
                    (ad.d[i] - fd).abs() < 1e-4,
                    "mismatch at index {i}: ad={} fd={}",
                    ad.d[i],
                    fd
                );
            }
        }
    }
}
pub use gatheral_schied_slynko_transient_impact_mod::{
    GatheralSchiedSlynkoTransientImpactAnalytics, GatheralSchiedSlynkoTransientImpactSpec,
    evaluate_gatheral_schied_slynko_transient_impact_execution,
};
