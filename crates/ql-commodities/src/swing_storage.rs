use crate::CommodityFuturesCurve;
use ql_ad::Real;

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
