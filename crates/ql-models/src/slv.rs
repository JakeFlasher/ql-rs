use crate::HestonParams;
use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct SlvLeverageSlice<T: Real> {
    pub time: f64,
    pub log_moneyness_grid: Vec<f64>,
    pub leverage: Vec<T>,
}
#[derive(Clone, Debug)]
pub struct CalibratedSlvModel<T: Real> {
    pub heston: HestonParams<T>,
    pub mixing_fraction: T,
    pub slices: Vec<SlvLeverageSlice<T>>,
}
impl<T: Real> CalibratedSlvModel<T> {
    pub fn leverage_at(&self, step_idx: usize, log_moneyness: T) -> T {
        let s = &self.slices[step_idx.min(self.slices.len().saturating_sub(1))];
        if s.leverage.len() == 1 {
            return s.leverage[0];
        }
        let x = log_moneyness.val();
        if x <= s.log_moneyness_grid[0] {
            return s.leverage[0];
        }
        if x >= *s.log_moneyness_grid.last().unwrap() {
            return *s.leverage.last().unwrap();
        }
        let i = s
            .log_moneyness_grid
            .partition_point(|z| *z <= x)
            .saturating_sub(1);
        let x0 = s.log_moneyness_grid[i];
        let x1 = s.log_moneyness_grid[i + 1];
        let w = (log_moneyness - T::from_f64(x0)) / T::from_f64((x1 - x0).max(1e-12));
        s.leverage[i] * (T::one() - w) + s.leverage[i + 1] * w
    }
}
fn normal(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let u = ((*seed >> 11) as f64 + 0.5) / 9007199254740992.0;
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let v = ((*seed >> 11) as f64 + 0.5) / 9007199254740992.0;
    (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()
}
#[allow(clippy::too_many_arguments)]
pub fn calibrate_slv_particle_leverage<T, LV>(
    spot: T,
    rate: T,
    div: T,
    expiry: f64,
    heston: &HestonParams<T>,
    mixing_fraction: T,
    dupire_local_vol: LV,
    n_steps: usize,
    n_particles: usize,
    n_grid: usize,
    seed: u64,
) -> CalibratedSlvModel<T>
where
    T: Real,
    LV: Fn(f64, T) -> T,
{
    let steps = n_steps.max(1);
    let np = n_particles.max(2);
    let ng = n_grid.max(3);
    let dt = expiry / steps as f64;
    let mut x = vec![T::zero(); np];
    let mut v = vec![heston.v0; np];
    let mut slices = Vec::with_capacity(steps);
    for m in 0..steps {
        let t = (m as f64 * dt).max(1e-8);
        let width = (heston.v0.val().max(1e-8).sqrt() * t.sqrt() * 3.5).max(0.1);
        let grid = (0..ng)
            .map(|i| -width + 2.0 * width * i as f64 / (ng - 1) as f64)
            .collect::<Vec<_>>();
        let mean_v = v.iter().copied().fold(T::zero(), |a, z| a + z) / np as f64;
        let sx = x.iter().map(|z| z.val()).sum::<f64>() / np as f64;
        let bw = (1.5 * x.iter().map(|z| (z.val() - sx).powi(2)).sum::<f64>().sqrt()
            / (np as f64).sqrt()
            * (np as f64).powf(-0.2))
        .max(0.02);
        let lev = grid
            .iter()
            .map(|&z| {
                let mut num = T::from_f64(1e-4) * mean_v;
                let mut den = T::from_f64(1e-4);
                for i in 0..np {
                    let u = (x[i] - T::from_f64(z)) / T::from_f64(bw);
                    let kw = (-u * u * 0.5).exp();
                    num += v[i] * kw;
                    den += kw;
                }
                let v_cond = (num / den).max_val(T::from_f64(1e-8));
                let dup = dupire_local_vol(
                    t,
                    spot * ((rate - div) * T::from_f64(t)).exp() * T::from_f64(z).exp(),
                );
                dup / v_cond.sqrt()
            })
            .map(|z| z.max_val(T::from_f64(0.1)).min_val(T::from_f64(5.0)))
            .collect();
        slices.push(SlvLeverageSlice {
            time: m as f64 * dt,
            log_moneyness_grid: grid,
            leverage: lev,
        });
        let mut next_x = vec![T::zero(); np];
        let mut next_v = vec![T::zero(); np];
        let sqdt = dt.sqrt();
        for i in 0..np {
            let z1 = normal(&mut (seed.wrapping_add((m * np + i) as u64)));
            let z2 = normal(&mut (seed.wrapping_add((m * np + i + np) as u64)));
            let zv = heston.rho * z1
                + (T::one() - heston.rho * heston.rho)
                    .max_val(T::zero())
                    .sqrt()
                    * z2;
            let vv = v[i].max_val(T::zero());
            next_v[i] = (vv
                + heston.kappa * (heston.theta - vv) * dt
                + mixing_fraction * heston.sigma_v * vv.sqrt() * T::from_f64(sqdt) * zv)
                .max_val(T::zero());
            let l = {
                let sl = &slices[m];
                let xx = x[i].val();
                if xx <= sl.log_moneyness_grid[0] {
                    sl.leverage[0]
                } else if xx >= *sl.log_moneyness_grid.last().unwrap() {
                    *sl.leverage.last().unwrap()
                } else {
                    let k = sl
                        .log_moneyness_grid
                        .partition_point(|z| *z <= xx)
                        .saturating_sub(1);
                    let w = (x[i] - T::from_f64(sl.log_moneyness_grid[k]))
                        / T::from_f64(sl.log_moneyness_grid[k + 1] - sl.log_moneyness_grid[k]);
                    sl.leverage[k] * (T::one() - w) + sl.leverage[k + 1] * w
                }
            };
            next_x[i] = x[i] - l * l * vv * (dt * 0.5)
                + l * vv.sqrt() * T::from_f64(sqdt) * T::from_f64(z1);
        }
        let mean_exp = next_x
            .iter()
            .copied()
            .map(|z| z.exp())
            .fold(T::zero(), |a, z| a + z)
            / np as f64;
        let shift = mean_exp.ln();
        for i in 0..np {
            x[i] = next_x[i] - shift;
            v[i] = next_v[i];
        }
    }
    let _ = spot;
    CalibratedSlvModel {
        heston: *heston,
        mixing_fraction,
        slices,
    }
}
#[allow(clippy::too_many_arguments)]
pub fn price_slv_european_mc<T: Real>(
    model: &CalibratedSlvModel<T>,
    spot: T,
    strike: T,
    rate: T,
    div: T,
    expiry: f64,
    n_paths: usize,
    is_call: bool,
    seed: u64,
) -> T {
    let n = n_paths.max(2);
    let steps = model.slices.len().max(1);
    let dt = expiry / steps as f64;
    let sqdt = dt.sqrt();
    let rc = (T::one() - model.heston.rho * model.heston.rho)
        .max_val(T::zero())
        .sqrt();
    let mut sum = T::zero();
    for p in 0..n {
        let sign = if p % 2 == 0 { 1. } else { -1. };
        let mut rng = seed.wrapping_add((p / 2) as u64);
        let mut x = T::zero();
        let mut v = model.heston.v0;
        for m in 0..steps {
            let z1 = normal(&mut rng) * sign;
            let z2 = normal(&mut rng) * sign;
            let zv = model.heston.rho * z1 + rc * z2;
            let vv = v.max_val(T::zero());
            let l = model.leverage_at(m, x);
            x += -l * l * vv * (dt * 0.5) + l * vv.sqrt() * T::from_f64(sqdt * z1);
            v = (vv
                + model.heston.kappa * (model.heston.theta - vv) * dt
                + model.mixing_fraction
                    * model.heston.sigma_v
                    * vv.sqrt()
                    * T::from_f64(sqdt)
                    * zv)
                .max_val(T::zero());
        }
        let st = spot * ((rate - div) * T::from_f64(expiry) + x).exp();
        sum += if is_call {
            (st - strike).max_val(T::zero())
        } else {
            (strike - st).max_val(T::zero())
        };
    }
    sum / n as f64 * (-rate * T::from_f64(expiry)).exp()
}
