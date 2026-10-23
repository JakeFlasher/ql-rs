use core::simd::Simd;
use core::simd::num::SimdFloat;
use ql_ad::{Real, SimdReal};
use ql_models::HestonParams;
use ql_vol::black_scholes_price;

#[derive(Clone, Debug)]
pub struct GbmSoaSimulator {
    pub spot: f64,
    pub vol: f64,
    pub rate: f64,
    pub div: f64,
    pub expiry: f64,
    pub n_paths: usize,
    pub n_steps: usize,
    pub seed: u64,
}

impl GbmSoaSimulator {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        spot: f64,
        vol: f64,
        rate: f64,
        div: f64,
        expiry: f64,
        n_paths: usize,
        n_steps: usize,
        seed: u64,
    ) -> Self {
        Self {
            spot,
            vol,
            rate,
            div,
            expiry,
            n_paths,
            n_steps,
            seed,
        }
    }

    pub fn price_european(&self, strike: f64, is_call: bool) -> (f64, f64) {
        let rng = crate::PhiloxRng::new(self.seed);
        let steps = self.n_steps.max(1);
        let dt = self.expiry / steps as f64;
        let drift = Simd::<f64, 8>::splat((self.rate - self.div - 0.5 * self.vol * self.vol) * dt);
        let vol_sqrt_dt = Simd::<f64, 8>::splat(self.vol * dt.sqrt());
        let k_vec = Simd::<f64, 8>::splat(strike);
        let zero = Simd::<f64, 8>::splat(0.0);
        let df = (-self.rate * self.expiry).exp();

        let n_batches = self.n_paths.max(1).div_ceil(16);
        let mut sum = 0.0;
        let mut ss = 0.0;
        for b in 0..n_batches {
            let mut log_s_pos = Simd::<f64, 8>::splat(self.spot.ln());
            let mut log_s_neg = log_s_pos;
            let base_path = (b as u64) * 2;
            for j in 0..steps {
                let z = rng.normal_simd8(base_path, j as u32);
                log_s_pos += drift + vol_sqrt_dt * z;
                log_s_neg += drift - vol_sqrt_dt * z;
            }
            let s_pos = SimdReal::exp(log_s_pos);
            let s_neg = SimdReal::exp(log_s_neg);
            let pay_pos = if is_call {
                (s_pos - k_vec).simd_max(zero)
            } else {
                (k_vec - s_pos).simd_max(zero)
            };
            let pay_neg = if is_call {
                (s_neg - k_vec).simd_max(zero)
            } else {
                (k_vec - s_neg).simd_max(zero)
            };
            let pos = pay_pos.to_array();
            let neg = pay_neg.to_array();
            let base = b * 16;
            for lane in 0..8 {
                if base + lane < self.n_paths {
                    let pay = pos[lane] * df;
                    sum += pay;
                    ss += pay * pay;
                }
                if base + lane + 8 < self.n_paths {
                    let pay = neg[lane] * df;
                    sum += pay;
                    ss += pay * pay;
                }
            }
        }
        let n_eff = self.n_paths.max(1) as f64;
        let mean = sum / n_eff;
        let se = ((ss / n_eff - mean * mean).max(0.0) / n_eff).sqrt();
        (mean, se)
    }

    pub fn analytical(&self, strike: f64, is_call: bool) -> f64 {
        black_scholes_price(
            self.spot,
            strike,
            self.vol,
            self.expiry,
            self.rate,
            self.div,
            is_call,
        )
    }

    pub fn price_european_qmc(&self, strike: f64, is_call: bool) -> f64 {
        let steps = self.n_steps.max(1);
        let mut sob = crate::SobolGenerator::new(steps);
        let times = (1..=steps)
            .map(|i| self.expiry * i as f64 / steps as f64)
            .collect::<Vec<_>>();
        let bridge = crate::BrownianBridge::new(times);
        let mut z = vec![0.0; steps];
        let mut dw = vec![0.0; steps];
        let dt = self.expiry / steps as f64;
        let drift = (self.rate - self.div - 0.5 * self.vol * self.vol) * dt;
        let mut sum = 0.0;
        let n = self.n_paths.max(1);
        for _ in 0..n {
            sob.next_normal_point(&mut z);
            bridge.transform(&z, &mut dw);
            let mut log_s_pos = self.spot.ln();
            let mut log_s_neg = log_s_pos;
            for &w in &dw {
                log_s_pos += drift + self.vol * w;
                log_s_neg += drift - self.vol * w;
            }
            let s_pos = log_s_pos.exp();
            let s_neg = log_s_neg.exp();
            let pay_pos = if is_call {
                (s_pos - strike).max(0.0)
            } else {
                (strike - s_pos).max(0.0)
            };
            let pay_neg = if is_call {
                (s_neg - strike).max(0.0)
            } else {
                (strike - s_neg).max(0.0)
            };
            sum += 0.5 * (pay_pos + pay_neg);
        }
        sum / n as f64 * (-self.rate * self.expiry).exp()
    }
}

/// Randomized Owen-scrambled Sobol Asian arithmetic-average pricer using the
/// same Brownian bridge construction as the QMC simulator.
#[allow(clippy::too_many_arguments)]
pub fn price_asian_arithmetic_rqmc(
    spot: f64,
    strike: f64,
    rate: f64,
    div: f64,
    vol: f64,
    fixing_times: &[f64],
    expiry: f64,
    is_call: bool,
    points_per_scramble: u32,
    scrambles: u32,
    seed: u32,
) -> (f64, f64) {
    assert!(
        fixing_times.len() <= 64,
        "RQMC Asian supports at most 64 fixings"
    );
    let times = fixing_times.to_vec();
    let bridge = crate::BrownianBridge::new(times.clone());
    crate::rqmc_estimate_with_stderr(
        times.len(),
        points_per_scramble,
        scrambles,
        seed,
        |normals| {
            let mut increments = [0.0f64; 64];
            bridge.transform(normals, &mut increments[..times.len()]);
            let mut log_pos = spot.ln();
            let mut log_neg = spot.ln();
            let mut sum_pos = 0.0;
            let mut sum_neg = 0.0;
            let mut prev = 0.0;
            for i in 0..times.len() {
                let t = times[i];
                let dt = t - prev;
                let drift = (rate - div - 0.5 * vol * vol) * dt;
                log_pos += drift + vol * increments[i];
                log_neg += drift - vol * increments[i];
                sum_pos += log_pos.exp();
                sum_neg += log_neg.exp();
                prev = t;
            }
            let avg_pos = sum_pos / times.len().max(1) as f64;
            let avg_neg = sum_neg / times.len().max(1) as f64;
            let payoff = |avg: f64| {
                if is_call {
                    (avg - strike).max(0.0)
                } else {
                    (strike - avg).max(0.0)
                }
            };
            0.5 * (payoff(avg_pos) + payoff(avg_neg)) * (-rate * expiry).exp()
        },
    )
}

#[derive(Clone, Debug)]
pub struct HestonQeSimulator {
    pub params: HestonParams<f64>,
    pub spot: f64,
    pub rate: f64,
    pub div: f64,
    pub expiry: f64,
    pub n_paths: usize,
    pub n_steps: usize,
    pub seed: u64,
}

impl HestonQeSimulator {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        params: HestonParams<f64>,
        spot: f64,
        rate: f64,
        div: f64,
        expiry: f64,
        n_paths: usize,
        n_steps: usize,
        seed: u64,
    ) -> Self {
        Self {
            params,
            spot,
            rate,
            div,
            expiry,
            n_paths,
            n_steps,
            seed,
        }
    }

    fn qe_step(&self, s: f64, v: f64, dt: f64, z_x: f64, z_v: f64) -> (f64, f64) {
        let kappa = self.params.kappa;
        let theta = self.params.theta;
        let sigma_v = self.params.sigma_v;
        let rho = self.params.rho;

        let e = (-kappa * dt).exp();
        let m = theta + (v - theta) * e;
        let s2 = if sigma_v.abs() < 1e-14 {
            0.0
        } else if kappa.abs() < 1e-14 {
            v * sigma_v * sigma_v * dt
        } else {
            (v * sigma_v * sigma_v * e / kappa * (1.0 - e)
                + theta * sigma_v * sigma_v / (2.0 * kappa) * (1.0 - e) * (1.0 - e))
                .max(0.0)
        };
        if sigma_v.abs() < 1e-14 {
            let v_int = if kappa.abs() < 1e-12 {
                v * dt
            } else {
                theta * dt + (v - theta) * (1.0 - e) / kappa
            }.max(0.0);
            let log_next = s + (self.rate - self.div - 0.5 * v_int / dt) * dt
                + v_int.sqrt() * z_x;
            return (log_next, m.max(0.0));
        }
        let s2 = s2.max(1e-16);
        let psi = s2 / (m * m).max(1e-16);
        let vv = if psi <= 1.5 {
            let inv_psi = 2.0 / psi;
            let b2 = inv_psi - 1.0 + (inv_psi * (inv_psi - 1.0)).sqrt();
            let b = b2.sqrt();
            let a = m / (1.0 + b2);
            a * (b + z_v) * (b + z_v)
        } else {
            let pp = (psi - 1.0) / (psi + 1.0);
            let beta = (1.0 - pp) / m.max(1e-16);
            let uu = ql_ad::Real::norm_cdf(z_v);
            if uu <= pp {
                0.0
            } else {
                ((1.0 - pp) / (1.0 - uu)).ln() / beta
            }
        }
        .max(0.0);

        // Andersen QE-M martingale correction (Eqs. 38--40).  The correction
        // is evaluated from the same conditional variance law as the state
        // transition, preserving E[S_{t+dt}|F_t] = S_t exp((r-q)dt).
        let g1 = 0.5; let g2 = 0.5;
        let k1 = g1 * dt * (kappa * rho / sigma_v - 0.5) - rho / sigma_v;
        let k2 = g2 * dt * (kappa * rho / sigma_v - 0.5) + rho / sigma_v;
        let k3 = g1 * dt * (1.0 - rho * rho);
        let k4 = g2 * dt * (1.0 - rho * rho);
        let a_corr = k2 + 0.5 * k4;
        let mut k0 = 0.0;
        if psi <= 1.5 {
            let inv_psi = 2.0 / psi;
            let b2 = inv_psi - 1.0 + (inv_psi * (inv_psi - 1.0)).sqrt();
            let b = b2.sqrt();
            let aa = m / (1.0 + b2);
            let den = 1.0 - 2.0 * a_corr * aa;
            if den > 1e-12 {
                k0 = -a_corr * b * b * aa / den + 0.5 * den.ln()
                    - (k1 + 0.5 * k3) * v;
            }
        } else {
            let pp = (psi - 1.0) / (psi + 1.0);
            let beta = (1.0 - pp) / m.max(1e-16);
            if beta > a_corr + 1e-12 {
                k0 = -(pp + beta * (1.0 - pp) / (beta - a_corr)).ln()
                    - (k1 + 0.5 * k3) * v;
            }
        }
        let s_next = s + (self.rate - self.div) * dt + k0 + k1 * v + k2 * vv
            + ((k3 * v + k4 * vv).max(0.0)).sqrt() * z_x;
        (s_next, vv)
    }

    pub fn price_european(&self, strike: f64, is_call: bool) -> (f64, f64) {
        let rng = crate::PhiloxRng::new(self.seed);
        let steps = self.n_steps.max(1);
        let dt = self.expiry / steps as f64;
        let df = (-self.rate * self.expiry).exp();
        let n = self.n_paths.max(1);
        let mut sum = 0.0;
        let mut ss = 0.0;
        for p in 0..n {
            let mut s_pos = self.spot.ln();
            let mut v_pos = self.params.v0;
            let mut s_neg = s_pos;
            let mut v_neg = v_pos;
            for j in 0..steps {
                let z = rng.normal4(p as u64, j as u32);
                (s_pos, v_pos) = self.qe_step(s_pos, v_pos, dt, z[0], z[1]);
                (s_neg, v_neg) = self.qe_step(s_neg, v_neg, dt, -z[0], -z[1]);
            }
            let st_pos = s_pos.exp();
            let st_neg = s_neg.exp();
            let pay_pos = if is_call {
                (st_pos - strike).max(0.0)
            } else {
                (strike - st_pos).max(0.0)
            };
            let pay_neg = if is_call {
                (st_neg - strike).max(0.0)
            } else {
                (strike - st_neg).max(0.0)
            };
            let pay = 0.5 * (pay_pos + pay_neg) * df;
            sum += pay;
            ss += pay * pay;
        }
        let mean = sum / n as f64;
        let se = ((ss / n as f64 - mean * mean).max(0.0) / n as f64).sqrt();
        (mean, se)
    }
}

#[allow(clippy::too_many_arguments)]
pub fn mc_gbm_pathwise_tape(
    spot: f64,
    vol: f64,
    rate: f64,
    div: f64,
    strike: f64,
    expiry: f64,
    n_paths: usize,
    seed: u64,
) -> (f64, f64, f64, f64, usize) {
    ql_ad::TapeScope::start_recording();
    let s = ql_ad::TapeScope::with_active(|t| t.input(spot));
    let v = ql_ad::TapeScope::with_active(|t| t.input(vol));
    let r = ql_ad::TapeScope::with_active(|t| t.input(rate));
    let d = ql_ad::TapeScope::with_active(|t| t.input(div));
    let k = ql_ad::TapeScope::with_active(|t| t.input(strike));
    let t = ql_ad::TapeScope::with_active(|t| t.input(expiry));
    let drift = (r - d - 0.5 * v * v) * t;
    let vol_sqrt_t = v * t.sqrt();
    let df = (-r * t).exp();
    ql_ad::TapeScope::set_mark();
    let mark = ql_ad::TapeScope::with_active(|tape| tape.mark());
    let mut mark_adj = vec![0.0; mark];
    let mut price = 0.0;
    let rng = crate::PhiloxRng::new(seed);
    for p in 0..n_paths {
        let z = rng.normal4(p as u64, 0)[0];
        let st_pos = s * (drift + vol_sqrt_t * z).exp();
        let st_neg = s * (drift - vol_sqrt_t * z).exp();
        let pay_pos = (st_pos - k).max_val(ql_ad::TapeVar::from_f64(0.0));
        let pay_neg = (st_neg - k).max_val(ql_ad::TapeVar::from_f64(0.0));
        let pay = (pay_pos + pay_neg) * 0.5 * df;
        price += pay.val;
        ql_ad::TapeScope::with_active(|tape| {
            let mut adj = vec![0.0; tape.len()];
            adj[pay.idx as usize] = 1.0;
            tape.propagate_to_mark(&mut adj);
            for i in 0..mark {
                mark_adj[i] += adj[i];
            }
            tape.rewind_to_mark();
        });
    }
    let peak = ql_ad::TapeScope::with_active(|tape| {
        tape.propagate_mark_to_start(&mut mark_adj);
        tape.len()
    });
    ql_ad::TapeScope::stop_recording();
    let n = n_paths.max(1) as f64;
    (
        price / n,
        mark_adj.get(s.idx as usize).copied().unwrap_or(0.0) / n,
        mark_adj.get(v.idx as usize).copied().unwrap_or(0.0) / n,
        mark_adj.get(r.idx as usize).copied().unwrap_or(0.0) / n,
        peak,
    )
}
