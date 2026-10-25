use crate::HestonParams;
use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct DoubleHestonKernelSlvSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub maturity: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub v1_0: T,
    pub kappa1: T,
    pub theta1: T,
    pub sigma1: T,
    pub rho1: T,
    pub v2_0: T,
    pub kappa2: T,
    pub theta2: T,
    pub sigma2: T,
    pub rho2: T,
    pub target_dupire_local_vol: T,
}
#[derive(Clone, Debug)]
pub struct DoubleHestonKernelSlvReport<T: Real = f64> {
    pub fast_integrated_var_mean: T,
    pub slow_integrated_var_mean: T,
    pub total_integrated_var_mean: T,
    pub cumulant_c2_variance: T,
    pub cumulant_c4_fourth: T,
    pub double_heston_atm_vol: T,
    pub double_heston_atm_skew: T,
    pub nadaraya_watson_conditional_var: T,
    pub gyongy_calibrated_leverage: T,
    pub slv_reconstructed_local_vol: T,
    pub gyongy_reconstruction_error: T,
}
pub fn evaluate_double_heston_cumulants_and_kernel_slv<T: Real>(
    s: &DoubleHestonKernelSlvSpec<T>,
) -> DoubleHestonKernelSlvReport<T> {
    let t = s.maturity.max_val(T::from_f64(1e-8));
    let m1 = s.theta1 * t
        + (s.v1_0 - s.theta1) * (T::one() - (-s.kappa1 * t).exp())
            / s.kappa1.max_val(T::from_f64(1e-8));
    let m2 = s.theta2 * t
        + (s.v2_0 - s.theta2) * (T::one() - (-s.kappa2 * t).exp())
            / s.kappa2.max_val(T::from_f64(1e-8));
    let term = m1 + m2;
    let vol = (term / t).max_val(T::from_f64(1e-12)).sqrt();
    let csv = |v: T, k: T, th: T, sig: T, rho: T| {
        rho * sig / k.max_val(T::from_f64(1e-8)).powf(2.0)
            * (th
                * (k * t - T::from_f64(2.0) * (T::one() - (-k * t).exp())
                    + T::from_f64(0.5) * (T::one() - (-T::from_f64(2.0) * k * t).exp()))
                + v * (T::one() - (-k * t).exp() - k * t * (-k * t).exp()))
    };
    let skew = (csv(s.v1_0, s.kappa1, s.theta1, s.sigma1, s.rho1)
        + csv(s.v2_0, s.kappa2, s.theta2, s.sigma2, s.rho2))
        / (T::from_f64(2.0) * term.max_val(T::from_f64(1e-12)).powf(1.5))
        * t.sqrt();
    const Z_GH9: [f64; 9] = [
        -4.51274586339978,
        -3.20542900285647,
        -2.07684798313409,
        -1.02325566378913,
        0.0,
        1.02325566378913,
        2.07684798313409,
        3.20542900285647,
        4.51274586339978,
    ];
    const W_GH9: [f64; 9] = [
        0.00002234584404,
        0.00278914132131,
        0.04991640676364,
        0.24409745015871,
        0.40634931182460,
        0.24409745015871,
        0.04991640676364,
        0.00278914132131,
        0.00002234584404,
    ];
    let v1_bar_t = s.theta1 + (s.v1_0 - s.theta1) * (-s.kappa1 * t).exp();
    let v2_bar_t = s.theta2 + (s.v2_0 - s.theta2) * (-s.kappa2 * t).exp();
    let sv1 = s.sigma1.abs()
        * (v1_bar_t.abs() * (T::one() - (-T::from_f64(2.0) * s.kappa1 * t).exp())
            / (T::from_f64(2.0) * s.kappa1.max_val(T::from_f64(1e-8))))
        .max_val(T::from_f64(1e-12))
        .sqrt();
    let sv2 = s.sigma2.abs()
        * (v2_bar_t.abs() * (T::one() - (-T::from_f64(2.0) * s.kappa2 * t).exp())
            / (T::from_f64(2.0) * s.kappa2.max_val(T::from_f64(1e-8))))
        .max_val(T::from_f64(1e-12))
        .sqrt();
    let sig_x = (m1 + m2).abs().max_val(T::from_f64(1e-12)).sqrt();
    let h_nw = (T::from_f64(0.5) * sig_x).max_val(T::from_f64(1e-3));
    let log_s0 = s.spot.max_val(T::from_f64(1e-8)).ln();
    let log_k = s.strike.max_val(T::from_f64(1e-8)).ln();
    let mut num_nw = T::zero();
    let mut den_nw = T::zero();
    for q in 0..9 {
        let zq = T::from_f64(Z_GH9[q]);
        let log_st_q = log_s0 + (s.risk_free_rate - s.dividend_yield) * t
            - T::from_f64(0.5) * (m1 + m2)
            + sig_x * zq;
        let v1_q = (v1_bar_t + s.rho1 * sv1 * zq).max_val(T::from_f64(1e-6));
        let v2_q = (v2_bar_t + s.rho2 * sv2 * zq).max_val(T::from_f64(1e-6));
        let diff_k = (log_st_q - log_k) / h_nw;
        let kw = T::from_f64(W_GH9[q]) * (-T::from_f64(0.5) * diff_k * diff_k).exp();
        num_nw += kw * (v1_q + v2_q);
        den_nw += kw;
    }
    let cond = num_nw / den_nw.max_val(T::from_f64(1e-16));
    let lev = s.target_dupire_local_vol / cond.max_val(T::from_f64(1e-12)).sqrt();
    DoubleHestonKernelSlvReport {
        fast_integrated_var_mean: m1,
        slow_integrated_var_mean: m2,
        total_integrated_var_mean: term,
        cumulant_c2_variance: term,
        cumulant_c4_fourth: T::from_f64(3.0) * term * term,
        double_heston_atm_vol: vol,
        double_heston_atm_skew: skew,
        nadaraya_watson_conditional_var: cond,
        gyongy_calibrated_leverage: lev,
        slv_reconstructed_local_vol: lev * cond.sqrt(),
        gyongy_reconstruction_error: (lev * cond.sqrt() - s.target_dupire_local_vol).abs(),
    }
}

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

#[cfg(test)]
mod phase27_tests {
    use super::*;
    use ql_ad::Dual;
    fn make<T: Real>(x: T) -> DoubleHestonKernelSlvSpec<T> {
        DoubleHestonKernelSlvSpec {
            spot: T::from_f64(100.0),
            strike: T::from_f64(100.0),
            maturity: x,
            risk_free_rate: T::from_f64(0.03),
            dividend_yield: T::zero(),
            v1_0: T::from_f64(0.02),
            kappa1: T::from_f64(2.0),
            theta1: T::from_f64(0.02),
            sigma1: T::from_f64(0.3),
            rho1: T::from_f64(-0.6),
            v2_0: T::from_f64(0.01),
            kappa2: T::from_f64(0.5),
            theta2: T::from_f64(0.01),
            sigma2: T::from_f64(0.2),
            rho2: T::from_f64(-0.3),
            target_dupire_local_vol: x,
        }
    }
    #[test]
    fn dual_fd_parity() {
        let h = 1e-5;
        let fd = (evaluate_double_heston_cumulants_and_kernel_slv(&make(0.2 + h))
            .gyongy_calibrated_leverage
            - evaluate_double_heston_cumulants_and_kernel_slv(&make(0.2 - h))
                .gyongy_calibrated_leverage)
            / (2.0 * h);
        let d = evaluate_double_heston_cumulants_and_kernel_slv(&make(Dual::<1>::variable(0.2, 0)))
            .gyongy_calibrated_leverage;
        assert!(d.d[0].abs() > 1e-6 && (d.d[0] - fd).abs() < 1e-5);
    }

    #[test]
    fn phase27_invariants() {
        let r = evaluate_double_heston_cumulants_and_kernel_slv(&make(1.0));
        assert!(
            r.cumulant_c2_variance > 0.0
                && r.cumulant_c4_fourth > 0.0
                && r.nadaraya_watson_conditional_var > 0.0
        );
    }
}
