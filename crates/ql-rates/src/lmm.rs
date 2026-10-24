use ql_ad::Real;
use ql_mc::PhiloxRng;

#[derive(Clone, Debug)]
pub struct ForwardMarketModel<T: Real> {
    pub tenor_times: Vec<f64>,
    pub initial_forwards: Vec<T>,
    pub initial_df_t0: T,
    pub displacements: Vec<T>,
    pub volatilities: Vec<T>,
    pub corr_decay: f64,
    pub corr_floor: f64,
    pub backward_looking_rfr: bool,
}
#[derive(Clone, Debug)]
pub struct TarnSpec<T: Real> {
    pub notional: T,
    pub target_coupon: T,
    pub fixed_rate: T,
    pub gearing: T,
    pub coupon_floor: T,
    pub coupon_cap: T,
    pub smoothing_eps: f64,
}
#[derive(Clone, Debug)]
pub struct CallableRangeAccrualSpec<T: Real> {
    pub notional: T,
    pub coupon_rate: T,
    pub lower_bound: T,
    pub upper_bound: T,
    pub call_period_indices: Vec<usize>,
    pub smoothing_eps: f64,
}
impl<T: Real> ForwardMarketModel<T> {
    pub fn discount_factor_0(&self, period_idx: usize) -> T {
        let i = period_idx.min(self.initial_forwards.len());
        let mut df = self.initial_df_t0;
        for j in 0..i {
            let tau = self.tenor_times[j + 1] - self.tenor_times[j];
            df /= T::one() + self.initial_forwards[j] * tau;
        }
        df
    }

    fn corr_chol(&self) -> Vec<Vec<f64>> {
        let m = self.initial_forwards.len();
        let floor = self.corr_floor.clamp(0.0, 0.999);
        let decay = self.corr_decay.max(0.0);
        let mut l = vec![vec![0.0; m]; m];
        for i in 0..m {
            for j in 0..=i {
                let ti = self.tenor_times.get(i).copied().unwrap_or(i as f64);
                let tj = self.tenor_times.get(j).copied().unwrap_or(j as f64);
                let c = if i == j {
                    1.0
                } else {
                    floor + (1.0 - floor) * (-decay * (ti - tj).abs()).exp()
                };
                let mut s = c;
                for k in 0..j {
                    s -= l[i][k] * l[j][k];
                }
                l[i][j] = if i == j {
                    s.max(1e-12).sqrt()
                } else {
                    s / l[j][j].max(1e-12)
                };
            }
        }
        l
    }

    fn simulate_terminal_paths(&self, n_paths: usize, seed: u64) -> Vec<(Vec<T>, Vec<T>)> {
        let m = self.initial_forwards.len();
        let l = self.corr_chol();
        let rng = PhiloxRng::new(seed);
        let df_terminal = self.discount_factor_0(m);
        let rfr_scale = if self.backward_looking_rfr {
            (1.0 / 3.0_f64).sqrt()
        } else {
            1.0
        };
        let mut out = Vec::with_capacity(n_paths.max(1));
        for p in 0..n_paths.max(1) {
            let mut f = self.initial_forwards.clone();
            let mut prev_t = 0.0;
            for step in 0..m {
                let t = self.tenor_times[step].max(prev_t);
                let dt = t - prev_t;
                prev_t = t;
                if dt <= 1e-12 {
                    continue;
                }
                let z_raw: Vec<f64> = (0..m)
                    .map(|j| rng.normal4(p as u64, (step * m + j) as u32)[0])
                    .collect();
                let mut z = vec![0.0; m];
                for i in step..m {
                    for j in 0..=i {
                        z[i] += l[i][j] * z_raw[j];
                    }
                }
                let mut f_pred = f.clone();
                let mut drift0 = vec![T::zero(); m];
                for i in step..m {
                    let si = self.displacements.get(i).copied().unwrap_or(T::zero());
                    let vi = self.volatilities.get(i).copied().unwrap_or(T::zero()) * rfr_scale;
                    let mut mu = T::zero();
                    for j in (i + 1)..m {
                        let tau_j = self.tenor_times[j + 1] - self.tenor_times[j];
                        let sj = self.displacements.get(j).copied().unwrap_or(T::zero());
                        let vj = self.volatilities.get(j).copied().unwrap_or(T::zero()) * rfr_scale;
                        let rho_ij: f64 = (0..=i.min(j)).map(|k| l[i][k] * l[j][k]).sum();
                        mu -= (f[j] + sj) * vj * (tau_j * rho_ij)
                            / (T::one() + f[j] * tau_j).max_val(T::from_f64(1e-6));
                    }
                    drift0[i] = mu * vi;
                    let fi_s = (f[i] + si).max_val(T::from_f64(1e-8));
                    f_pred[i] = fi_s
                        * ((drift0[i] - vi * vi * 0.5) * dt + vi * (dt.sqrt() * z[i])).exp()
                        - si;
                }
                for i in step..m {
                    let si = self.displacements.get(i).copied().unwrap_or(T::zero());
                    let vi = self.volatilities.get(i).copied().unwrap_or(T::zero()) * rfr_scale;
                    let mut mu1 = T::zero();
                    for j in (i + 1)..m {
                        let tau_j = self.tenor_times[j + 1] - self.tenor_times[j];
                        let sj = self.displacements.get(j).copied().unwrap_or(T::zero());
                        let vj = self.volatilities.get(j).copied().unwrap_or(T::zero()) * rfr_scale;
                        let rho_ij: f64 = (0..=i.min(j)).map(|k| l[i][k] * l[j][k]).sum();
                        mu1 -= (f_pred[j] + sj) * vj * (tau_j * rho_ij)
                            / (T::one() + f_pred[j] * tau_j).max_val(T::from_f64(1e-6));
                    }
                    let drift = (drift0[i] + mu1 * vi) * 0.5;
                    let fi_s = (f[i] + si).max_val(T::from_f64(1e-8));
                    f[i] =
                        fi_s * ((drift - vi * vi * 0.5) * dt + vi * (dt.sqrt() * z[i])).exp() - si;
                }
            }
            let mut numeraire_dfs = vec![df_terminal; m];
            let mut acc = df_terminal;
            for i in (0..m).rev() {
                numeraire_dfs[i] = acc;
                let tau = self.tenor_times[i + 1] - self.tenor_times[i];
                acc *= (T::one() + f[i] * tau).max_val(T::from_f64(1e-6));
            }
            out.push((f, numeraire_dfs));
        }
        out
    }

    pub fn price_caplet_mc(&self, period_idx: usize, strike: T, n_paths: usize, seed: u64) -> T {
        if period_idx >= self.initial_forwards.len() || self.tenor_times.len() < period_idx + 2 {
            return T::zero();
        }
        let tau = self.tenor_times[period_idx + 1] - self.tenor_times[period_idx];
        let paths = self.simulate_terminal_paths(n_paths, seed);
        let mut sum = T::zero();
        for (f, dfs) in &paths {
            sum += (f[period_idx] - strike).max_val(T::zero()) * tau * dfs[period_idx];
        }
        sum / paths.len() as f64
    }

    pub fn price_tarn_mc(&self, spec: &TarnSpec<T>, n_paths: usize, seed: u64) -> T {
        let m = self.initial_forwards.len();
        let paths = self.simulate_terminal_paths(n_paths, seed);
        let eps = spec.smoothing_eps.max(1e-6);
        let mut total = T::zero();
        for (f, dfs) in &paths {
            let mut accrued = T::zero();
            let mut pv = T::zero();
            for i in 0..m {
                let tau = self.tenor_times[i + 1] - self.tenor_times[i];
                let c = (spec.fixed_rate - spec.gearing * f[i])
                    .max_val(spec.coupon_floor)
                    .min_val(spec.coupon_cap);
                let gate = spec.target_coupon.smooth_heaviside_c2(accrued, eps);
                let remaining = (spec.target_coupon - accrued).max_val(T::zero());
                let period_cpn = c * tau;
                let paid = period_cpn.min_val(remaining) * gate;
                pv += paid * dfs[i];
                accrued += period_cpn;
            }
            total += pv * spec.notional;
        }
        total / paths.len() as f64
    }

    pub fn price_callable_range_accrual_mc(
        &self,
        spec: &CallableRangeAccrualSpec<T>,
        n_paths: usize,
        seed: u64,
    ) -> T {
        let m = self.initial_forwards.len();
        let paths = self.simulate_terminal_paths(n_paths, seed);
        let np = paths.len();
        let eps = spec.smoothing_eps.max(1e-6);
        let mut cf = vec![vec![T::zero(); m]; np];
        for (p, (f, dfs)) in paths.iter().enumerate() {
            for i in 0..m {
                let ind = f[i].smooth_heaviside_c2(spec.lower_bound, eps)
                    * (T::one() - f[i].smooth_heaviside_c2(spec.upper_bound, eps));
                let tau = self.tenor_times[i + 1] - self.tenor_times[i];
                cf[p][i] = spec.notional * spec.coupon_rate * tau * ind * dfs[i];
            }
        }
        let mut cont: Vec<T> = (0..np).map(|p| cf[p][m - 1]).collect();
        for i in (0..(m - 1)).rev() {
            if spec.call_period_indices.contains(&(i + 1)) {
                let f_mean = paths.iter().map(|(f, _)| f[i].val()).sum::<f64>() / np as f64;
                let y_mean = cont.iter().map(|v| v.val()).sum::<f64>() / np as f64;
                let mut cov = 0.0;
                let mut var = 1e-12;
                for p in 0..np {
                    let dx = paths[p].0[i].val() - f_mean;
                    cov += dx * (cont[p].val() - y_mean);
                    var += dx * dx;
                }
                let beta = cov / var;
                let rem_tau = self.tenor_times[m] - self.tenor_times[i + 1];
                let call_threshold =
                    (spec.notional * spec.coupon_rate * T::from_f64(0.5 * rem_tau))
                        * self.discount_factor_0(i + 1);
                for p in 0..np {
                    let expected_cont = y_mean + beta * (paths[p].0[i].val() - f_mean);
                    if expected_cont > call_threshold.val() {
                        cont[p] = call_threshold;
                    }
                }
            }
            for p in 0..np {
                cont[p] += cf[p][i];
            }
        }
        cont.into_iter().fold(T::zero(), |a, z| a + z) / np as f64
    }
}

pub fn price_tarn_lmm<T: Real>(
    m: &ForwardMarketModel<T>,
    s: &TarnSpec<T>,
    np: usize,
    seed: u64,
) -> T {
    m.price_tarn_mc(s, np, seed)
}

pub fn price_callable_range_accrual_lmm<T: Real>(
    m: &ForwardMarketModel<T>,
    s: &CallableRangeAccrualSpec<T>,
    np: usize,
    seed: u64,
) -> T {
    m.price_callable_range_accrual_mc(s, np, seed)
}
