use crate::CommodityFuturesCurve;
use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct WarinBaeuerleStorageRollingHedgeSpec<T: Real = f64> {
    pub forward_curve: Vec<T>,
    pub decision_times: Vec<T>,
    pub max_capacity: T,
    pub initial_inventory: T,
    pub max_injection_step: T,
    pub max_withdrawal_step: T,
    pub injection_loss_fraction: T,
    pub spot_ou_volatility: T,
    pub spot_ou_mean_reversion: T,
    pub discount_rate: T,
}
#[derive(Clone, Debug)]
pub struct WarinBaeuerleStorageRollingHedgeReport<T: Real = f64> {
    pub static_intrinsic_value: T,
    pub stochastic_total_storage_value: T,
    pub extrinsic_optionality_value: T,
    pub forward_pillar_deltas: Vec<T>,
    pub unhedged_pnl_variance: T,
    pub warin_delta_hedged_pnl_variance: T,
    pub hedge_variance_reduction_ratio: T,
}
pub fn evaluate_warin_baeuerle_storage_rolling_hedge<T: Real>(
    s: &WarinBaeuerleStorageRollingHedgeSpec<T>,
) -> WarinBaeuerleStorageRollingHedgeReport<T> {
    const Z_GH5: [f64; 5] = [
        -2.8569700139,
        -1.3556261799,
        0.0,
        1.3556261799,
        2.8569700139,
    ];
    const W_GH5: [f64; 5] = [
        0.0112574113,
        0.2220759220,
        0.5333333333,
        0.2220759220,
        0.0112574113,
    ];
    let n = s.forward_curve.len();
    let cap = s.max_capacity.max_val(T::from_f64(1e-8));
    let inv_step = cap / T::from_f64(10.0);
    let init_idx = (s.initial_inventory / inv_step)
        .val()
        .round()
        .clamp(0.0, 10.0) as usize;
    let tau = T::from_f64(0.05);
    let softmax = |vals: &[T]| -> T {
        let m = vals
            .iter()
            .map(|v| v.val())
            .fold(f64::NEG_INFINITY, f64::max);
        let mm = T::from_f64(m);
        let sum = vals
            .iter()
            .fold(T::zero(), |acc, &v| acc + ((v - mm) / tau).exp());
        mm + tau * sum.ln()
    };
    let run_dp = |spots: &[T]| -> (T, Vec<T>) {
        let mut next = vec![T::zero(); 11];
        let mut initial_path = vec![T::zero(); n];
        for i in (0..n).rev() {
            let mut cur = vec![T::zero(); 11];
            let dt = s.decision_times.get(i).copied().unwrap_or(T::zero());
            let df = (-s.discount_rate * dt).exp();
            for g in 0..11 {
                let inv = T::from_f64(g as f64) * inv_step;
                let inj = s.max_injection_step.min_val((cap - inv).max_val(T::zero()));
                let wd = s.max_withdrawal_step.min_val(inv.max_val(T::zero()));
                let gi = ((inv + inj) / inv_step).val().round().clamp(0.0, 10.0) as usize;
                let gw = ((inv - wd) / inv_step).val().round().clamp(0.0, 10.0) as usize;
                let hold = df * next[g];
                let inject =
                    df * (next[gi] - spots[i] * inj * (T::one() + s.injection_loss_fraction));
                let withdraw =
                    df * (next[gw] + spots[i] * wd * (T::one() - s.injection_loss_fraction));
                cur[g] = softmax(&[hold, inject, withdraw]);
            }
            initial_path[i] = cur[init_idx];
            next = cur;
        }
        (next[init_idx], initial_path)
    };
    let static_spots = s.forward_curve.clone();
    let (static_value, _) = run_dp(&static_spots);
    let mut sigma_i = Vec::with_capacity(n);
    let mut spot_nodes = Vec::with_capacity(n);
    let mut value_nodes = vec![vec![T::zero(); 5]; n];
    for i in 0..n {
        let ti = s
            .decision_times
            .get(i)
            .copied()
            .unwrap_or(T::zero())
            .max_val(T::zero());
        let kappa = s.spot_ou_mean_reversion.abs().max_val(T::from_f64(1e-8));
        let si = s.spot_ou_volatility.abs()
            * ((T::one() - (-T::from_f64(2.0) * kappa * ti).exp()) / (T::from_f64(2.0) * kappa))
                .max_val(T::from_f64(1e-16))
                .sqrt();
        sigma_i.push(si);
        spot_nodes.push(
            Z_GH5
                .iter()
                .map(|&z| {
                    s.forward_curve[i] * (si * T::from_f64(z) - T::from_f64(0.5) * si * si).exp()
                })
                .collect::<Vec<_>>(),
        );
    }
    for q in 0..5 {
        let spots = (0..n).map(|i| spot_nodes[i][q]).collect::<Vec<_>>();
        let (_, path) = run_dp(&spots);
        for i in 0..n {
            value_nodes[i][q] = path[i];
        }
    }
    let stochastic_total = if n == 0 {
        static_value
    } else {
        (0..5).fold(T::zero(), |acc, q| {
            acc + T::from_f64(W_GH5[q]) * value_nodes[0][q]
        })
    };
    let mut deltas = vec![T::zero(); n];
    let mut unhedged = T::zero();
    let mut hedged = T::zero();
    for i in 0..n {
        let mean_v = (0..5).fold(T::zero(), |a, q| {
            a + T::from_f64(W_GH5[q]) * value_nodes[i][q]
        });
        let mean_s = (0..5).fold(T::zero(), |a, q| {
            a + T::from_f64(W_GH5[q]) * spot_nodes[i][q]
        });
        let cov = (0..5).fold(T::zero(), |a, q| {
            a + T::from_f64(W_GH5[q]) * (value_nodes[i][q] - mean_v) * (spot_nodes[i][q] - mean_s)
        });
        let var_s = (0..5).fold(T::zero(), |a, q| {
            a + T::from_f64(W_GH5[q]) * (spot_nodes[i][q] - mean_s) * (spot_nodes[i][q] - mean_s)
        });
        deltas[i] = cov / var_s.max_val(T::from_f64(1e-12));
        let var_v = (0..5).fold(T::zero(), |a, q| {
            a + T::from_f64(W_GH5[q]) * (value_nodes[i][q] - mean_v) * (value_nodes[i][q] - mean_v)
        });
        let var_h = (0..5).fold(T::zero(), |a, q| {
            let r = value_nodes[i][q] - deltas[i] * (spot_nodes[i][q] - s.forward_curve[i]);
            a + T::from_f64(W_GH5[q]) * (r - mean_v) * (r - mean_v)
        });
        unhedged += var_v;
        hedged += var_h;
    }
    let extr = (stochastic_total - static_value).max_val(T::zero());
    WarinBaeuerleStorageRollingHedgeReport {
        static_intrinsic_value: static_value,
        stochastic_total_storage_value: stochastic_total.max_val(static_value),
        extrinsic_optionality_value: extr,
        forward_pillar_deltas: deltas,
        unhedged_pnl_variance: unhedged,
        warin_delta_hedged_pnl_variance: hedged,
        hedge_variance_reduction_ratio: T::one() - hedged / unhedged.max_val(T::from_f64(1e-12)),
    }
}

#[derive(Clone, Debug)]
pub struct SwingOptionSpec<T: Real> {
    pub exercise_times: Vec<f64>,
    pub strike: T,
    pub unit_volume: T,
    pub min_exercises: usize,
    pub max_exercises: usize,
    pub is_call: bool,
}

#[derive(Clone, Debug)]
pub struct GasStorageSpec<T: Real> {
    pub decision_times: Vec<f64>,
    pub max_capacity: f64,
    pub initial_inventory: f64,
    pub max_injection_rate: f64,
    pub max_withdrawal_rate: f64,
    pub injection_loss: T,
    pub withdrawal_loss: T,
    pub injection_cost: T,
    pub withdrawal_cost: T,
    pub inventory_grid_points: usize,
}

#[derive(Clone, Debug)]
struct OuTree<T: Real> {
    times: Vec<f64>,
    dx: T,
    states: Vec<Vec<T>>,
    probs: Vec<Vec<T>>,
    j_bounds: Vec<isize>,
    kappa: T,
    sigma: T,
}

fn ou_var<T: Real>(kappa: T, sigma: T, dt: f64) -> T {
    let dt_t = T::from_f64(dt);
    let u = kappa * dt_t;
    if u.val().abs() < 1e-6 {
        sigma
            * sigma
            * (dt_t - kappa * dt_t * dt_t + kappa * kappa * dt_t * dt_t * dt_t * (2.0 / 3.0))
    } else {
        sigma * sigma * (T::one() - (-kappa * 2.0 * dt_t).exp()) / (kappa * 2.0)
    }
}

fn make_times(raw: &[f64], steps_per_interval: usize) -> Vec<f64> {
    if raw.is_empty() {
        return vec![0.0];
    }
    let mut sorted = raw.to_vec();
    sorted.sort_by(f64::total_cmp);
    sorted.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    let mut times = vec![0.0];
    let mut prev = 0.0;
    let n = steps_per_interval.max(1);
    for &end in &sorted {
        if end <= prev + 1e-14 {
            continue;
        }
        for i in 1..=n {
            times.push(prev + (end - prev) * i as f64 / n as f64);
        }
        prev = end;
    }
    times
}

fn build_tree<T: Real>(times: Vec<f64>, kappa: T, sigma: T) -> OuTree<T> {
    let max_dt = times
        .windows(2)
        .map(|w| (w[1] - w[0]).max(0.0))
        .fold(0.0, f64::max)
        .max(1e-12);
    let dx = (ou_var(kappa, sigma, max_dt) * 3.0)
        .max_val(T::from_f64(1e-16))
        .sqrt()
        .max_val(T::from_f64(1e-8));
    let j_max = if kappa.val() * max_dt > 1e-4 {
        (0.184 / (kappa.val() * max_dt)).ceil() as isize
    } else {
        (times.len() - 1).max(1) as isize
    };
    let mut j_bounds = Vec::with_capacity(times.len());
    let mut states: Vec<Vec<T>> = Vec::with_capacity(times.len());
    let mut probs: Vec<Vec<T>> = Vec::with_capacity(times.len());
    for m in 0..times.len() {
        let jb = (m as isize).min(j_max);
        j_bounds.push(jb);
        states.push((-jb..=jb).map(|j| dx * T::from_f64(j as f64)).collect());
        probs.push(vec![T::zero(); (2 * jb + 1) as usize]);
    }
    probs[0][0] = T::one();
    let dx2 = dx * dx;
    for m in 0..times.len().saturating_sub(1) {
        let dt = (times[m + 1] - times[m]).max(1e-12);
        let v = ou_var(kappa, sigma, dt);
        let jbn = j_bounds[m + 1];
        for (idx, q) in probs[m].clone().into_iter().enumerate() {
            if q.val() == 0.0 {
                continue;
            }
            let x = states[m][idx];
            let mu = x * (-kappa * T::from_f64(dt)).exp();
            let k = ((mu / dx).val().round() as isize).clamp(-jbn + 1, jbn - 1);
            let e = mu - dx * T::from_f64(k as f64);
            let pd = (v + e * e - e * dx) / (dx2 * 2.0);
            let pm = T::one() - (v + e * e) / dx2;
            let pu = (v + e * e + e * dx) / (dx2 * 2.0);
            for (node, prob) in [(k - 1, q * pd), (k, q * pm), (k + 1, q * pu)] {
                let ii = (node + jbn).clamp(0, 2 * jbn) as usize;
                probs[m + 1][ii] += prob;
            }
        }
    }
    OuTree {
        times,
        dx,
        states,
        probs,
        j_bounds,
        kappa,
        sigma,
    }
}

fn transition<T: Real>(tree: &OuTree<T>, m: usize, j: isize) -> (isize, T, T, T, T) {
    let dt = (tree.times[m + 1] - tree.times[m]).max(1e-12);
    let v = ou_var(tree.kappa, tree.sigma, dt);
    let dx2 = tree.dx * tree.dx;
    let mu = tree.dx * T::from_f64(j as f64) * (-tree.kappa * T::from_f64(dt)).exp();
    let jbn = tree.j_bounds[m + 1];
    let k = ((mu / tree.dx).val().round() as isize).clamp(-jbn + 1, jbn - 1);
    let e = mu - tree.dx * T::from_f64(k as f64);
    let pd = (v + e * e - e * tree.dx) / (dx2 * 2.0);
    let pm = T::one() - (v + e * e) / dx2;
    let pu = (v + e * e + e * tree.dx) / (dx2 * 2.0);
    (k, pd, pm, pu, T::from_f64(dt))
}

fn node_spots<T: Real>(tree: &OuTree<T>, curve: &CommodityFuturesCurve<T>, m: usize) -> Vec<T> {
    let fwd = curve.futures_price(T::from_f64(tree.times[m].max(1e-8)));
    let den = tree.probs[m]
        .iter()
        .zip(&tree.states[m])
        .fold(T::zero(), |a, (p, x)| a + *p * x.exp());
    tree.states[m]
        .iter()
        .map(|x| fwd * x.exp() / den.max_val(T::from_f64(1e-14)))
        .collect()
}

pub fn price_swing_option<T: Real>(
    spec: &SwingOptionSpec<T>,
    curve: &CommodityFuturesCurve<T>,
    kappa: T,
    sigma: T,
    rate: T,
    steps_per_interval: usize,
) -> T {
    if spec.exercise_times.is_empty() {
        return T::zero();
    }
    let mut exercise_times = spec.exercise_times.clone();
    exercise_times.sort_by(f64::total_cmp);
    exercise_times.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    if exercise_times.is_empty() {
        return T::zero();
    }
    let times = make_times(&exercise_times, steps_per_interval);
    let tree = build_tree(times, kappa, sigma);
    let mut ex_rank = vec![None; tree.times.len()];
    for (r, t) in exercise_times.iter().enumerate() {
        if let Some(m) = tree.times.iter().position(|x| (*x - *t).abs() < 1e-10) {
            ex_rank[m] = Some(r);
        }
    }
    let n = tree.times.len().saturating_sub(1);
    let umax = spec.max_exercises.min(exercise_times.len());
    let terminal_width = (2 * tree.j_bounds[n] + 1) as usize;
    let mut values = vec![vec![T::zero(); terminal_width]; umax + 1];
    let final_spots = node_spots(&tree, curve, n);
    let final_remaining = 1usize;
    for (u, row) in values.iter_mut().enumerate().take(umax + 1) {
        for (idx, spot) in final_spots.iter().enumerate() {
            if u == umax {
                continue;
            }
            let cf = spec.unit_volume
                * if spec.is_call {
                    *spot - spec.strike
                } else {
                    spec.strike - *spot
                };
            row[idx] = if u + final_remaining <= spec.min_exercises {
                cf
            } else {
                T::zero().max_val(cf)
            };
        }
    }
    for m in (0..n).rev() {
        let jb = tree.j_bounds[m];
        let jbn = tree.j_bounds[m + 1];
        let spots = node_spots(&tree, curve, m);
        let mut next = vec![vec![T::zero(); (2 * jb + 1) as usize]; umax + 1];
        for u in 0..=umax {
            for j in -jb..=jb {
                let idx = (j + jb) as usize;
                let (k, pd, pm, pu, dt) = transition(&tree, m, j);
                let val = |node: isize| values[u][(node + jbn).clamp(0, 2 * jbn) as usize];
                next[u][idx] =
                    (-rate * dt).exp() * (pd * val(k - 1) + pm * val(k) + pu * val(k + 1));
            }
        }
        if let Some(rank) = ex_rank[m] {
            let remaining = exercise_times.len() - rank;
            for u in 0..=umax {
                for j in -jb..=jb {
                    let idx = (j + jb) as usize;
                    let cf = spec.unit_volume
                        * if spec.is_call {
                            spots[idx] - spec.strike
                        } else {
                            spec.strike - spots[idx]
                        };
                    if u == umax {
                        continue;
                    }
                    let take = cf + next[u + 1][idx];
                    next[u][idx] = if u + remaining <= spec.min_exercises {
                        take
                    } else {
                        next[u][idx].max_val(take)
                    };
                }
            }
        }
        values = next;
    }
    values[0][0]
}

fn interpolate<T: Real>(values: &[T], x: f64, dv: f64) -> T {
    if x <= 0.0 {
        return values[0];
    }
    let g = values.len() - 1;
    if x >= g as f64 * dv {
        return values[g];
    }
    let z = x / dv;
    let i = z.floor() as usize;
    let w = z - i as f64;
    values[i] * (1.0 - w) + values[i + 1] * w
}

pub fn price_gas_storage<T: Real>(
    spec: &GasStorageSpec<T>,
    curve: &CommodityFuturesCurve<T>,
    kappa: T,
    sigma: T,
    rate: T,
) -> T {
    if spec.decision_times.is_empty() || spec.max_capacity <= 0.0 {
        return T::zero();
    }
    let times = make_times(&spec.decision_times, 1);
    let tree = build_tree(times, kappa, sigma);
    let g = spec.inventory_grid_points.max(2);
    let dv = spec.max_capacity / g as f64;
    let n = tree.times.len().saturating_sub(1);
    let terminal_width = (2 * tree.j_bounds[n] + 1) as usize;
    let mut values = vec![vec![T::zero(); terminal_width]; g + 1];
    let mut decision = vec![false; tree.times.len()];
    for t in &spec.decision_times {
        if let Some(m) = tree.times.iter().position(|x| (*x - *t).abs() < 1e-10) {
            decision[m] = true;
        }
    }
    if n > 0 && decision[n] {
        let dt_n = (tree.times[n] - tree.times[n - 1]).max(1e-12);
        let final_spots = node_spots(&tree, curve, n);
        for (gi, row) in values.iter_mut().enumerate().take(g + 1) {
            let inv = gi as f64 * dv;
            let max_withdraw = inv.min(spec.max_withdrawal_rate * dt_n);
            for (idx, spot) in final_spots.iter().enumerate() {
                let withdraw_cash = ((T::one() - spec.withdrawal_loss) * *spot
                    - spec.withdrawal_cost)
                    * max_withdraw;
                row[idx] = T::zero().max_val(withdraw_cash);
            }
        }
    }
    for m in (0..tree.times.len().saturating_sub(1)).rev() {
        let jb = tree.j_bounds[m];
        let jbn = tree.j_bounds[m + 1];
        let spots = node_spots(&tree, curve, m);
        let mut next = vec![vec![T::zero(); (2 * jb + 1) as usize]; g + 1];
        let mut cont_by_j = vec![vec![T::zero(); g + 1]; (2 * jb + 1) as usize];
        for j in -jb..=jb {
            let (k, pd, pm, pu, dt) = transition(&tree, m, j);
            let row = &mut cont_by_j[(j + jb) as usize];
            for gp in 0..=g {
                let val = |node: isize| values[gp][(node + jbn).clamp(0, 2 * jbn) as usize];
                row[gp] = (-rate * dt).exp() * (pd * val(k - 1) + pm * val(k) + pu * val(k + 1));
            }
        }
        for gi in 0..=g {
            let inv = gi as f64 * dv;
            for j in -jb..=jb {
                let idx = (j + jb) as usize;
                let (_, _, _, _, dt) = transition(&tree, m, j);
                let cont = &cont_by_j[idx];
                let lo = (inv - spec.max_withdrawal_rate * dt.val()).max(0.0);
                let hi = (inv + spec.max_injection_rate * dt.val()).min(spec.max_capacity);
                let mut best = T::from_f64(-1e100);
                for (gp, cont_value) in cont.iter().enumerate().take(g + 1) {
                    let vp = gp as f64 * dv;
                    if vp + 1e-12 < lo || vp - 1e-12 > hi {
                        continue;
                    }
                    let dvv = vp - inv;
                    let cash = if dvv > 0.0 {
                        -((T::one() + spec.injection_loss) * spots[idx] + spec.injection_cost) * dvv
                    } else if dvv < 0.0 {
                        ((T::one() - spec.withdrawal_loss) * spots[idx] - spec.withdrawal_cost)
                            * (-dvv)
                    } else {
                        T::zero()
                    };
                    best = best.max_val(cash + *cont_value);
                }
                for vp in [lo, hi] {
                    let dvv = vp - inv;
                    let cash = if dvv > 0.0 {
                        -((T::one() + spec.injection_loss) * spots[idx] + spec.injection_cost) * dvv
                    } else if dvv < 0.0 {
                        ((T::one() - spec.withdrawal_loss) * spots[idx] - spec.withdrawal_cost)
                            * (-dvv)
                    } else {
                        T::zero()
                    };
                    best = best.max_val(cash + interpolate(cont, vp, dv));
                }
                next[gi][idx] = if decision[m] { best } else { cont[gi] };
            }
        }
        values = next;
    }
    let initial = spec.initial_inventory.clamp(0.0, spec.max_capacity);
    interpolate(
        &values.iter().map(|x| x[0]).collect::<Vec<_>>(),
        initial,
        dv,
    )
}

#[derive(Clone, Debug)]
pub struct TakeOrPayContractSpec<T: Real = f64> {
    pub decision_times: Vec<f64>,
    pub strike: T,
    pub daily_min_qty: f64,
    pub daily_max_qty: f64,
    pub cumulative_min_qty: f64,
    pub cumulative_max_qty: f64,
    pub shortfall_penalty_rate: T,
    pub daily_discretization_steps: usize,
}

#[derive(Copy, Clone, Debug)]
pub struct TakeOrPayValuation<T: Real> {
    pub npv: T,
    pub min_take_swap_npv: T,
    pub swing_optionality_value: T,
}

#[allow(clippy::needless_range_loop)]
pub fn price_take_or_pay_contract<T: Real>(
    spec: &TakeOrPayContractSpec<T>,
    curve: &CommodityFuturesCurve<T>,
    kappa: T,
    sigma: T,
    rate: T,
) -> TakeOrPayValuation<T> {
    if spec.decision_times.is_empty() || spec.daily_max_qty < spec.daily_min_qty {
        return TakeOrPayValuation {
            npv: T::zero(),
            min_take_swap_npv: T::zero(),
            swing_optionality_value: T::zero(),
        };
    }
    let mut decision_times = spec.decision_times.clone();
    decision_times.sort_by(f64::total_cmp);
    decision_times.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    decision_times.retain(|&t| t >= 0.0);
    if decision_times.is_empty() {
        return TakeOrPayValuation {
            npv: T::zero(),
            min_take_swap_npv: T::zero(),
            swing_optionality_value: T::zero(),
        };
    }
    let times = make_times(&decision_times, 1);
    let tree = build_tree(times, kappa, sigma);
    let mut dec_rank = vec![None; tree.times.len()];
    for (r, &t) in decision_times.iter().enumerate() {
        if let Some(m) = tree.times.iter().position(|x| (*x - t).abs() < 1e-10) {
            dec_rank[m] = Some(r);
        }
    }
    let n_dec = decision_times.len();
    let m_steps = tree.times.len();
    let dq_span = (spec.daily_max_qty - spec.daily_min_qty).max(0.0);
    let l_sub = spec.daily_discretization_steps.max(1);
    let dq = if dq_span <= 1e-12 {
        1.0
    } else {
        dq_span / (l_sub as f64)
    };
    let max_units = if dq_span <= 1e-12 { 0 } else { n_dec * l_sub };
    let n_last = m_steps.saturating_sub(1);
    let term_w = (2 * tree.j_bounds[n_last] + 1) as usize;
    let mut values = vec![vec![T::zero(); term_w]; max_units + 1];

    let final_spots = node_spots(&tree, curve, n_last);
    let prev_dec_last = n_dec.saturating_sub(1);
    let max_u_last = if dq_span <= 1e-12 {
        0
    } else {
        prev_dec_last * l_sub
    };
    for u in 0..=max_u_last {
        let q_cum_before = (prev_dec_last as f64) * spec.daily_min_qty + (u as f64) * dq;
        for (idx, &spot) in final_spots.iter().enumerate() {
            let mut best = T::from_f64(-1e100);
            let max_step = if dq_span <= 1e-12 { 0 } else { l_sub };
            for step_u in 0..=max_step {
                let q_take = spec.daily_min_qty + (step_u as f64) * dq;
                let q_tot = q_cum_before + q_take;
                if q_tot > spec.cumulative_max_qty + 1e-9 {
                    continue;
                }
                let shortfall = (spec.cumulative_min_qty - q_tot).max(0.0);
                let val = (spot - spec.strike) * q_take - spec.shortfall_penalty_rate * shortfall;
                best = best.max_val(val);
            }
            values[u][idx] = best;
        }
    }

    for m in (0..n_last).rev() {
        let jb = tree.j_bounds[m];
        let jbn = tree.j_bounds[m + 1];
        let spots = node_spots(&tree, curve, m);
        let mut next = vec![vec![T::zero(); (2 * jb + 1) as usize]; max_units + 1];
        if let Some(r) = dec_rank[m] {
            let max_u_m = if dq_span <= 1e-12 { 0 } else { r * l_sub };
            for u in 0..=max_u_m {
                let q_cum_before = (r as f64) * spec.daily_min_qty + (u as f64) * dq;
                for j in -jb..=jb {
                    let idx = (j + jb) as usize;
                    let (k, pd, pm, pu, dt) = transition(&tree, m, j);
                    let df = (-rate * dt).exp();
                    let mut best = T::from_f64(-1e100);
                    let max_step = if dq_span <= 1e-12 { 0 } else { l_sub };
                    for step_u in 0..=max_step {
                        let q_take = spec.daily_min_qty + (step_u as f64) * dq;
                        if q_cum_before + q_take > spec.cumulative_max_qty + 1e-9 {
                            continue;
                        }
                        let u_next = (u + step_u).min(max_units);
                        let val_n =
                            |node: isize| values[u_next][(node + jbn).clamp(0, 2 * jbn) as usize];
                        let cont = df * (pd * val_n(k - 1) + pm * val_n(k) + pu * val_n(k + 1));
                        let total = (spots[idx] - spec.strike) * q_take + cont;
                        best = best.max_val(total);
                    }
                    next[u][idx] = best;
                }
            }
        } else {
            for j in -jb..=jb {
                let idx = (j + jb) as usize;
                let (k, pd, pm, pu, dt) = transition(&tree, m, j);
                let df = (-rate * dt).exp();
                let val_n = |node: isize| values[0][(node + jbn).clamp(0, 2 * jbn) as usize];
                next[0][idx] = df * (pd * val_n(k - 1) + pm * val_n(k) + pu * val_n(k + 1));
            }
        }
        values = next;
    }

    let min_per_step =
        (spec.cumulative_min_qty / (n_dec as f64)).clamp(spec.daily_min_qty, spec.daily_max_qty);
    let mut swap_npv = T::zero();
    for &t in &decision_times {
        let df = (-rate * t).exp();
        let fwd = curve.futures_price(T::from_f64(t));
        swap_npv += df * (fwd - spec.strike) * min_per_step;
    }
    let npv = values[0][0];
    TakeOrPayValuation {
        npv,
        min_take_swap_npv: swap_npv,
        swing_optionality_value: (npv - swap_npv).max_val(T::zero()),
    }
}

#[derive(Clone, Debug)]
pub struct BessSpec<T: Real = f64> {
    pub decision_times: Vec<f64>,
    pub energy_capacity_mwh: f64,
    pub initial_soc_mwh: f64,
    pub max_charge_mw: f64,
    pub max_discharge_mw: f64,
    pub charge_efficiency: T,
    pub discharge_efficiency: T,
    pub cycle_degradation_cost: T,
    pub ancillary_reserve_rate: T,
    pub soc_grid_points: usize,
}

#[derive(Copy, Clone, Debug)]
pub struct BessValuation<T: Real> {
    pub total_npv: T,
    pub intrinsic_npv: T,
    pub extrinsic_optionality: T,
}

#[allow(clippy::needless_range_loop)]
fn solve_bess_tree<T: Real>(
    spec: &BessSpec<T>,
    curve: &CommodityFuturesCurve<T>,
    kappa: T,
    sigma: T,
    rate: T,
) -> T {
    if spec.decision_times.is_empty() || spec.energy_capacity_mwh <= 0.0 {
        return T::zero();
    }
    let times = make_times(&spec.decision_times, 1);
    let tree = build_tree(times, kappa, sigma);
    let mut decision = vec![false; tree.times.len()];
    for &t in &spec.decision_times {
        if let Some(m) = tree.times.iter().position(|x| (*x - t).abs() < 1e-10) {
            decision[m] = true;
        }
    }
    let g = spec.soc_grid_points.max(2);
    let de = spec.energy_capacity_mwh / (g as f64);
    let n = tree.times.len().saturating_sub(1);
    let term_w = (2 * tree.j_bounds[n] + 1) as usize;
    let mut values = vec![vec![T::zero(); term_w]; g + 1];
    let eta_c = spec
        .charge_efficiency
        .max_val(T::from_f64(0.1))
        .min_val(T::one());
    let eta_d = spec
        .discharge_efficiency
        .max_val(T::from_f64(0.1))
        .min_val(T::one());

    if n > 0 {
        let dt_n = (tree.times[n] - tree.times[n - 1]).max(1e-6);
        let final_spots = node_spots(&tree, curve, n);
        for (gi, row) in values.iter_mut().enumerate().take(g + 1) {
            let soc = (gi as f64) * de;
            let max_dis = soc.min(spec.max_discharge_mw * dt_n);
            for (idx, &spot) in final_spots.iter().enumerate() {
                let cash = (eta_d * spot - spec.cycle_degradation_cost) * max_dis;
                row[idx] = cash.max_val(T::zero());
            }
        }
    }

    for m in (0..n).rev() {
        let jb = tree.j_bounds[m];
        let jbn = tree.j_bounds[m + 1];
        let spots = node_spots(&tree, curve, m);
        let mut next = vec![vec![T::zero(); (2 * jb + 1) as usize]; g + 1];
        let mut cont_by_j = vec![vec![T::zero(); g + 1]; (2 * jb + 1) as usize];
        for j in -jb..=jb {
            let (k, pd, pm, pu, dt) = transition(&tree, m, j);
            let row = &mut cont_by_j[(j + jb) as usize];
            for gp in 0..=g {
                let val = |node: isize| values[gp][(node + jbn).clamp(0, 2 * jbn) as usize];
                row[gp] = (-rate * dt).exp() * (pd * val(k - 1) + pm * val(k) + pu * val(k + 1));
            }
        }
        for gi in 0..=g {
            let soc = (gi as f64) * de;
            let headroom = soc.min(spec.energy_capacity_mwh - soc).max(0.0);
            for j in -jb..=jb {
                let idx = (j + jb) as usize;
                if !decision[m] {
                    next[gi][idx] = cont_by_j[idx][gi];
                    continue;
                }
                let (_, _, _, _, dt) = transition(&tree, m, j);
                let dt_val = dt.val().max(1e-6);
                let lo = (soc - spec.max_discharge_mw * dt_val).max(0.0);
                let hi = (soc + spec.max_charge_mw * dt_val).min(spec.energy_capacity_mwh);
                let anc = spec.ancillary_reserve_rate * (headroom * dt_val);
                let mut best = T::from_f64(-1e100);
                for gp in 0..=g {
                    let soc_next = (gp as f64) * de;
                    if soc_next + 1e-12 < lo || soc_next - 1e-12 > hi {
                        continue;
                    }
                    let d_soc = soc_next - soc;
                    let cash = if d_soc > 0.0 {
                        -(spots[idx] / eta_c + spec.cycle_degradation_cost) * d_soc
                    } else if d_soc < 0.0 {
                        (eta_d * spots[idx] - spec.cycle_degradation_cost) * (-d_soc)
                    } else {
                        T::zero()
                    };
                    best = best.max_val(cash + anc + cont_by_j[idx][gp]);
                }
                for soc_next in [lo, hi] {
                    let d_soc = soc_next - soc;
                    let cash = if d_soc > 0.0 {
                        -(spots[idx] / eta_c + spec.cycle_degradation_cost) * d_soc
                    } else if d_soc < 0.0 {
                        (eta_d * spots[idx] - spec.cycle_degradation_cost) * (-d_soc)
                    } else {
                        T::zero()
                    };
                    let cont = interpolate(&cont_by_j[idx], soc_next, de);
                    best = best.max_val(cash + anc + cont);
                }
                next[gi][idx] = best;
            }
        }
        values = next;
    }
    let init_soc = spec.initial_soc_mwh.clamp(0.0, spec.energy_capacity_mwh);
    interpolate(
        &values.iter().map(|r| r[0]).collect::<Vec<_>>(),
        init_soc,
        de,
    )
}

pub fn price_bess_arbitrage_dp<T: Real>(
    spec: &BessSpec<T>,
    power_curve: &CommodityFuturesCurve<T>,
    kappa: T,
    sigma: T,
    rate: T,
) -> BessValuation<T> {
    let total = solve_bess_tree(spec, power_curve, kappa, sigma, rate);
    let intrinsic = solve_bess_tree(spec, power_curve, kappa, T::from_f64(1e-6), rate);
    BessValuation {
        total_npv: total,
        intrinsic_npv: intrinsic,
        extrinsic_optionality: (total - intrinsic).max_val(T::zero()),
    }
}

#[derive(Clone, Debug)]
pub struct BessAnalyticalSdpDegradationStorageSpec<T: Real = f64> {
    pub energy_capacity_mwh: T,
    pub max_power_mw: T,
    pub initial_soc_fraction: T,
    pub min_soc_fraction: T,
    pub max_soc_fraction: T,
    pub charge_efficiency: T,
    pub discharge_efficiency: T,
    pub battery_capex_usd: T,
    pub rated_cycle_life: T,
    pub peak_mean_price_mwh: T,
    pub offpeak_mean_price_mwh: T,
    pub price_volatility_mwh: T,
    pub horizon_days: usize,
    pub discount_rate_annual: T,
}

#[derive(Clone, Debug)]
pub struct BessAnalyticalSdpDegradationStorageReport<T: Real = f64> {
    pub round_trip_efficiency: T,
    pub marginal_degradation_cost_per_mwh: T,
    pub shadow_water_value_per_mwh: T,
    pub optimal_charge_threshold_price: T,
    pub optimal_discharge_threshold_price: T,
    pub hysteresis_deadband_width: T,
    pub negative_price_charge_probability: T,
    pub expected_daily_net_arbitrage_usd: T,
    pub expected_annual_equivalent_cycles: T,
    pub total_bess_npv_usd: T,
    pub gross_arbitrage_without_degradation_usd: T,
}

pub fn evaluate_bess_analytical_sdp_degradation_storage<T: Real>(
    s: &BessAnalyticalSdpDegradationStorageSpec<T>,
) -> BessAnalyticalSdpDegradationStorageReport<T> {
    let rt = s.charge_efficiency * s.discharge_efficiency;
    let delta_e =
        s.energy_capacity_mwh * (s.max_soc_fraction - s.min_soc_fraction).max_val(T::zero());
    let deg = s.battery_capex_usd
        / (T::from_f64(2.0)
            * s.rated_cycle_life.max_val(T::from_f64(1e-8))
            * delta_e.max_val(T::from_f64(1e-8)));
    let midpoint = (s.min_soc_fraction + s.max_soc_fraction) * T::from_f64(0.5);
    let water = (s.discharge_efficiency * s.peak_mean_price_mwh
        + s.offpeak_mean_price_mwh / s.charge_efficiency.max_val(T::from_f64(1e-8)))
        * T::from_f64(0.5)
        * (T::one() - T::from_f64(0.1) * (s.initial_soc_fraction - midpoint));
    let charge = s.charge_efficiency * water - deg;
    let discharge = (water + deg) / s.discharge_efficiency.max_val(T::from_f64(1e-8));
    let deadband = (discharge - charge).max_val(T::from_f64(1e-8));
    let neg_prob =
        (-s.offpeak_mean_price_mwh / s.price_volatility_mwh.max_val(T::from_f64(1e-8))).norm_cdf();
    let v_init = (s.initial_soc_fraction - s.min_soc_fraction).max_val(T::zero())
        * s.energy_capacity_mwh
        * s.discharge_efficiency
        * water;
    let daily = delta_e
        * (s.discharge_efficiency * s.peak_mean_price_mwh
            - s.offpeak_mean_price_mwh / s.charge_efficiency.max_val(T::from_f64(1e-8))
            + T::from_f64(0.4) * s.price_volatility_mwh);
    let daily_deg = T::from_f64(2.0) * delta_e * deg;
    let mut annuity = T::zero();
    for d in 1..=s.horizon_days {
        annuity +=
            (T::one() + s.discount_rate_annual / T::from_f64(365.0)).pow(T::from_f64(-(d as f64)));
    }
    let gross = v_init + annuity * daily;
    let npv = gross - annuity * daily_deg;
    BessAnalyticalSdpDegradationStorageReport {
        round_trip_efficiency: rt,
        marginal_degradation_cost_per_mwh: deg,
        shadow_water_value_per_mwh: water,
        optimal_charge_threshold_price: charge,
        optimal_discharge_threshold_price: discharge,
        hysteresis_deadband_width: deadband,
        negative_price_charge_probability: neg_prob,
        expected_daily_net_arbitrage_usd: daily,
        expected_annual_equivalent_cycles: T::from_f64(365.0)
            * (delta_e / s.energy_capacity_mwh.max_val(T::from_f64(1e-8))),
        total_bess_npv_usd: npv,
        gross_arbitrage_without_degradation_usd: gross,
    }
}

#[cfg(test)]
mod phase17_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn take_or_pay_bang_bang_and_bess_arbitrage_ad() {
        let curve =
            CommodityFuturesCurve::new(vec![0.25, 0.50, 0.75, 1.0], vec![45.0, 58.0, 42.0, 62.0])
                .unwrap();
        let top1 = TakeOrPayContractSpec {
            decision_times: vec![0.25, 0.50, 0.75, 1.0],
            strike: 50.0,
            daily_min_qty: 1.0,
            daily_max_qty: 3.0,
            cumulative_min_qty: 6.0,
            cumulative_max_qty: 10.0,
            shortfall_penalty_rate: 50.0,
            daily_discretization_steps: 1,
        };
        let mut top4 = top1.clone();
        top4.daily_discretization_steps = 4;
        let v1 = price_take_or_pay_contract(&top1, &curve, 1.5, 0.30, 0.03);
        let v4 = price_take_or_pay_contract(&top4, &curve, 1.5, 0.30, 0.03);
        assert!((v1.npv - v4.npv).abs() < 1e-10);
        assert!(v1.swing_optionality_value > 0.0);

        let mut top_frac1 = top1.clone();
        top_frac1.cumulative_min_qty = 7.0;
        top_frac1.cumulative_max_qty = 7.0;
        top_frac1.daily_discretization_steps = 1;
        let mut top_frac2 = top_frac1.clone();
        top_frac2.daily_discretization_steps = 2;
        let vf1 = price_take_or_pay_contract(&top_frac1, &curve, 1.5, 0.30, 0.03);
        let vf2 = price_take_or_pay_contract(&top_frac2, &curve, 1.5, 0.30, 0.03);
        assert!(vf2.npv > vf1.npv);

        let bess = BessSpec {
            decision_times: vec![0.25, 0.50, 0.75, 1.0],
            energy_capacity_mwh: 4.0,
            initial_soc_mwh: 2.0,
            max_charge_mw: 8.0,
            max_discharge_mw: 8.0,
            charge_efficiency: 0.92,
            discharge_efficiency: 0.92,
            cycle_degradation_cost: 1.5,
            ancillary_reserve_rate: 2.0,
            soc_grid_points: 4,
        };
        let bv = price_bess_arbitrage_dp(&bess, &curve, 1.5, 0.35, 0.03);
        assert!(bv.total_npv >= bv.intrinsic_npv - 1e-10);
        assert!(bv.extrinsic_optionality >= 0.0);

        let h = 1e-5;
        let eval_eff = |eta: f64| {
            let mut b2 = bess.clone();
            b2.discharge_efficiency = eta;
            price_bess_arbitrage_dp(&b2, &curve, 1.5, 0.35, 0.03).total_npv
        };
        let fd = (eval_eff(0.92 + h) - eval_eff(0.92 - h)) / (2.0 * h);
        type D = Dual<1>;
        let curve_d = CommodityFuturesCurve::new(
            vec![0.25, 0.50, 0.75, 1.0],
            vec![
                D::from_value(45.0),
                D::from_value(58.0),
                D::from_value(42.0),
                D::from_value(62.0),
            ],
        )
        .unwrap();
        let bess_d = BessSpec {
            decision_times: vec![0.25, 0.50, 0.75, 1.0],
            energy_capacity_mwh: 4.0,
            initial_soc_mwh: 2.0,
            max_charge_mw: 8.0,
            max_discharge_mw: 8.0,
            charge_efficiency: D::from_value(0.92),
            discharge_efficiency: D::variable(0.92, 0),
            cycle_degradation_cost: D::from_value(1.5),
            ancillary_reserve_rate: D::from_value(2.0),
            soc_grid_points: 4,
        };
        let ad = price_bess_arbitrage_dp(
            &bess_d,
            &curve_d,
            D::from_value(1.5),
            D::from_value(0.35),
            D::from_value(0.03),
        )
        .total_npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase27_tests {
    use super::*;
    use ql_ad::Dual;
    fn make<T: Real>(x: T) -> WarinBaeuerleStorageRollingHedgeSpec<T> {
        WarinBaeuerleStorageRollingHedgeSpec {
            forward_curve: vec![x, T::from_f64(3.1), T::from_f64(3.2)],
            decision_times: vec![T::from_f64(0.5), T::one(), T::from_f64(1.5)],
            max_capacity: T::from_f64(10.0),
            initial_inventory: T::from_f64(2.0),
            max_injection_step: T::one(),
            max_withdrawal_step: T::one(),
            injection_loss_fraction: T::from_f64(0.02),
            spot_ou_volatility: T::from_f64(0.3),
            spot_ou_mean_reversion: T::one(),
            discount_rate: T::from_f64(0.03),
        }
    }
    #[test]
    fn dual_fd_parity() {
        let h = 1e-5;
        let fd = (evaluate_warin_baeuerle_storage_rolling_hedge(&make(3.0 + h))
            .stochastic_total_storage_value
            - evaluate_warin_baeuerle_storage_rolling_hedge(&make(3.0 - h))
                .stochastic_total_storage_value)
            / (2.0 * h);
        let d = evaluate_warin_baeuerle_storage_rolling_hedge(&make(Dual::<1>::variable(3.0, 0)))
            .stochastic_total_storage_value;
        assert!(d.d[0].abs() > 1e-6 && (d.d[0] - fd).abs() < 1e-5);
    }

    #[test]
    fn phase27_invariants() {
        let r = evaluate_warin_baeuerle_storage_rolling_hedge(&make(3.0));
        assert!(r.stochastic_total_storage_value >= r.static_intrinsic_value);
        assert!(r.extrinsic_optionality_value > 0.0);
        assert!(r.hedge_variance_reduction_ratio > 0.5);
    }
}

#[cfg(test)]
mod phase28_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(x: T) -> BessAnalyticalSdpDegradationStorageSpec<T> {
        BessAnalyticalSdpDegradationStorageSpec {
            energy_capacity_mwh: T::from_f64(100.0) + x,
            max_power_mw: T::from_f64(50.0),
            initial_soc_fraction: T::from_f64(0.5),
            min_soc_fraction: T::from_f64(0.1),
            max_soc_fraction: T::from_f64(0.9),
            charge_efficiency: T::from_f64(0.9),
            discharge_efficiency: T::from_f64(0.9),
            battery_capex_usd: T::from_f64(100000.0),
            rated_cycle_life: T::from_f64(4000.0),
            peak_mean_price_mwh: T::from_f64(100.0),
            offpeak_mean_price_mwh: T::from_f64(20.0),
            price_volatility_mwh: T::from_f64(30.0),
            horizon_days: 365,
            discount_rate_annual: T::from_f64(0.08),
        }
    }
    #[test]
    fn bess_invariants_and_dual() {
        let r = evaluate_bess_analytical_sdp_degradation_storage(&mk(0.0));
        assert!(r.optimal_discharge_threshold_price > r.optimal_charge_threshold_price);
        assert!(r.hysteresis_deadband_width > 0.0);
        assert!(r.gross_arbitrage_without_degradation_usd >= r.total_bess_npv_usd);
        let h = 1e-5;
        let fd = (evaluate_bess_analytical_sdp_degradation_storage(&mk(h)).total_bess_npv_usd
            - evaluate_bess_analytical_sdp_degradation_storage(&mk(-h)).total_bess_npv_usd)
            / (2.0 * h);
        let ad = evaluate_bess_analytical_sdp_degradation_storage(&mk(Dual::<1>::variable(0.0, 0)))
            .total_bess_npv_usd;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
