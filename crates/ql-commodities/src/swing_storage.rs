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

fn node_tree<T: Real>(times: &[f64], kappa: T, sigma: T) -> (Vec<Vec<T>>, Vec<Vec<T>>) {
    let mut states = vec![vec![T::zero()]; times.len()];
    let mut probs = vec![vec![T::one()]; times.len()];
    for m in 1..times.len() {
        let dt = (times[m] - times[m - 1]).max(1e-12);
        let v = if kappa.val().abs() < 1e-8 {
            sigma * sigma * dt
        } else {
            sigma * sigma * (T::one() - (-kappa * 2.0 * dt).exp()) / (kappa * 2.0)
        };
        let dx = (v * 3.0).sqrt();
        let prev = states[m - 1].clone();
        let mut cur = vec![T::zero(); 2 * m + 1];
        for (j, node) in cur.iter_mut().enumerate() {
            *node = T::from_f64(j as f64 - m as f64) * dx;
        }
        states[m] = cur;
        let mut p = vec![T::zero(); 2 * m + 1];
        for prev_node in prev.iter() {
            let mu = *prev_node * (-kappa * dt).exp();
            let k =
                ((mu / dx).val().round() as isize + m as isize).clamp(0, (2 * m) as isize) as usize;
            p[k] += T::from_f64(1.0 / 3.0);
            p[(k + 1).min(2 * m)] += T::from_f64(1.0 / 3.0);
            p[k.saturating_sub(1)] += T::from_f64(1.0 / 3.0);
        }
        probs[m] = p;
    }
    (states, probs)
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
    let mut times = vec![0.0];
    for w in spec.exercise_times.windows(2) {
        let n = steps_per_interval.max(1);
        for i in 1..=n {
            times.push(w[0] + (w[1] - w[0]) * i as f64 / n as f64);
        }
    }
    if times.len() == 1 {
        times.push(spec.exercise_times[0]);
    }
    let (states, _) = node_tree(&times, kappa, sigma);
    let m = times.len() - 1;
    let width = 2 * m + 1;
    let u_max = spec.max_exercises.min(m);
    let mut v = vec![vec![T::zero(); width]; u_max + 1];
    for k in (0..m).rev() {
        let t = times[k];
        let spot0 = curve.futures_price(T::from_f64(t.max(1e-8)));
        for u in 0..=u_max {
            let d = m - k;
            let s = spot0 * (states[k][0].exp());
            let cf = spec.unit_volume
                * if spec.is_call {
                    s - spec.strike
                } else {
                    spec.strike - s
                };
            let cont = if u <= u_max { v[u][0] } else { T::zero() };
            let take = if u < u_max {
                cf + v[u + 1][0]
            } else {
                T::zero()
            };
            let out = if u == u_max {
                cont
            } else if u + d <= spec.min_exercises {
                take
            } else {
                cont.max_val(take)
            };
            v[u][0] = out;
        }
    }
    v[0][0] * (-rate * spec.exercise_times[0]).exp()
}

pub fn price_gas_storage<T: Real>(
    spec: &GasStorageSpec<T>,
    curve: &CommodityFuturesCurve<T>,
    _kappa: T,
    _sigma: T,
    rate: T,
) -> T {
    if spec.decision_times.len() < 2 {
        return T::zero();
    }
    let g = spec.inventory_grid_points.max(2);
    let dv = spec.max_capacity / g as f64;
    let mut u = vec![T::zero(); g + 1];
    for k in (0..spec.decision_times.len() - 1).rev() {
        let dt = (spec.decision_times[k + 1] - spec.decision_times[k]).max(0.0);
        let s = curve.futures_price(T::from_f64(spec.decision_times[k]));
        let mut next = vec![T::zero(); g + 1];
        for (i, next_value) in next.iter_mut().enumerate().take(g + 1) {
            let inv = i as f64 * dv;
            let lo = (inv - spec.max_withdrawal_rate * dt).max(0.0);
            let hi = (inv + spec.max_injection_rate * dt).min(spec.max_capacity);
            let mut best = T::from_f64(-1e100);
            for (j, u_value) in u.iter().enumerate().take(g + 1) {
                let vp = j as f64 * dv;
                if vp + 1e-12 < lo || vp - 1e-12 > hi {
                    continue;
                }
                let dvv = vp - inv;
                let cash = if dvv > 0.0 {
                    -((T::one() + spec.injection_loss) * s + spec.injection_cost) * dvv
                } else if dvv < 0.0 {
                    ((T::one() - spec.withdrawal_loss) * s - spec.withdrawal_cost) * (-dvv)
                } else {
                    T::zero()
                };
                best = best.max_val(cash + *u_value);
            }
            *next_value = best;
        }
        u = next;
    }
    let i = (spec.initial_inventory / dv).round().clamp(0.0, g as f64) as usize;
    u[i] * (-rate * spec.decision_times[0]).exp()
}
