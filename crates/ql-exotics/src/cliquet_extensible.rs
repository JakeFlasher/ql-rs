use ql_ad::Real;
use ql_math::bivariate_normal_cdf;
use ql_vol::black_scholes_price;

pub struct ForwardStartOptionEngine;
impl ForwardStartOptionEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        spot: T,
        moneyness_alpha: T,
        rate: T,
        div: T,
        vol: T,
        reset_time: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        let tau = (expiry - reset_time).max_val(T::zero());
        let f = spot * ((rate - div) * reset_time).exp();
        black_scholes_price(f, moneyness_alpha * f, vol, tau, rate, div, is_call)
            * (-rate * reset_time).exp()
    }
}
#[derive(Clone, Debug)]
pub struct CliquetOptionSpec<T: Real = f64> {
    pub spot: T,
    pub rate: T,
    pub div: T,
    pub vol: T,
    pub reset_times: Vec<f64>,
    pub local_floor: T,
    pub local_cap: T,
    pub notional: T,
    pub global_floor: T,
    pub global_cap: T,
}
pub struct AnalyticCliquetEngine;
impl AnalyticCliquetEngine {
    pub fn price<T: Real>(s: &CliquetOptionSpec<T>, is_call: bool) -> T {
        let mut v = T::zero();
        for w in s.reset_times.windows(2) {
            let dt = T::from_f64(w[1] - w[0]);
            if dt.val() <= 0.0 {
                continue;
            }
            let forward_ratio = ((s.rate - s.div) * dt).exp();
            let df = (-s.rate * T::from_f64(w[1])).exp();
            if is_call {
                let floor_k = T::one() + s.local_floor;
                let cap_k = T::one() + s.local_cap;
                let floor = s.local_floor
                    + black_scholes_price(
                        forward_ratio,
                        floor_k,
                        s.vol,
                        dt,
                        T::zero(),
                        T::zero(),
                        true,
                    )
                    - black_scholes_price(
                        forward_ratio,
                        cap_k,
                        s.vol,
                        dt,
                        T::zero(),
                        T::zero(),
                        true,
                    );
                v += df * floor;
            } else {
                let floor_k = T::one() - s.local_floor;
                let cap_k = T::one() - s.local_cap;
                let floor = s.local_floor
                    + black_scholes_price(
                        forward_ratio,
                        floor_k,
                        s.vol,
                        dt,
                        T::zero(),
                        T::zero(),
                        false,
                    )
                    - black_scholes_price(
                        forward_ratio,
                        cap_k,
                        s.vol,
                        dt,
                        T::zero(),
                        T::zero(),
                        false,
                    );
                v += df * floor;
            }
        }
        s.notional * v
    }

    /// Monte-Carlo fallback for globally floored/capped cliquets.  The local
    /// analytic engine above is exact only when resets are independently
    /// aggregated; this path simulation preserves the global barriers.
    pub fn price_mc_global<T: Real>(s: &CliquetOptionSpec<T>, n_paths: usize, seed: u64) -> T {
        let rng = ql_mc::PhiloxRng::new(seed);
        let mut sum = T::zero();
        for p in 0..n_paths.max(1) {
            let mut spot = s.spot;
            let mut total = T::zero();
            for (j, w) in s.reset_times.windows(2).enumerate() {
                let dt = T::from_f64(w[1] - w[0]);
                if dt.val() <= 0.0 {
                    continue;
                }
                let z = rng.normal4(p as u64, j as u32)[0];
                let drift = (s.rate - s.div - s.vol * s.vol * 0.5) * dt;
                let diffusion = s.vol * dt.sqrt() * T::from_f64(z);
                let next = spot * (drift + diffusion).exp();
                let r = next / spot - T::one();
                let local = r.max_val(s.local_floor).min_val(s.local_cap);
                spot = next;
                total += local;
            }
            let payoff = total.max_val(s.global_floor).min_val(s.global_cap);
            let maturity = T::from_f64(*s.reset_times.last().unwrap_or(&0.0));
            let pv = s.notional * (-s.rate * maturity).exp() * payoff;
            sum += pv;
        }
        sum / T::from_f64(n_paths.max(1) as f64)
    }
}
pub struct WriterExtensibleOptionEngine;
impl WriterExtensibleOptionEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        spot: T,
        k1: T,
        k2: T,
        rate: T,
        div: T,
        vol: T,
        t1: T,
        t2: T,
        is_call: bool,
    ) -> T {
        let st1 = vol * t1.sqrt();
        let st2 = vol * t2.sqrt();
        let d11 = ((spot / k1).ln() + (rate - div + vol * vol * 0.5) * t1) / st1;
        let d21 = d11 - st1;
        let d12 = ((spot / k2).ln() + (rate - div + vol * vol * 0.5) * t2) / st2;
        let d22 = d12 - st2;
        let rho = -(t1 / t2).sqrt();
        if is_call {
            black_scholes_price(spot, k1, vol, t1, rate, div, true)
                + spot * (-div * t2).exp() * bivariate_normal_cdf(-d11, d12, rho)
                - k2 * (-rate * t2).exp() * bivariate_normal_cdf(-d21, d22, rho)
        } else {
            black_scholes_price(spot, k1, vol, t1, rate, div, false)
                + k2 * (-rate * t2).exp() * bivariate_normal_cdf(d21, -d22, rho)
                - spot * (-div * t2).exp() * bivariate_normal_cdf(d11, -d12, rho)
        }
    }
}

#[derive(Clone, Debug)]
pub struct HimalayaOptionSpec<T: Real = f64> {
    pub spots: Vec<T>,
    pub vols: Vec<T>,
    pub div_yields: Vec<T>,
    pub correlation: T,
    pub rate: T,
    pub observation_times: Vec<f64>,
    pub strike: T,
    pub notional: T,
    pub temperature: f64,
}

pub fn price_himalaya_option_mc_ad<T: Real>(
    s: &HimalayaOptionSpec<T>,
    n_paths: usize,
    seed: u64,
) -> T {
    let m = s
        .spots
        .len()
        .min(s.vols.len())
        .min(s.div_yields.len())
        .min(s.observation_times.len());
    if m == 0 {
        return T::zero();
    }
    let paths = n_paths.max(2);
    let rho_min = -1.0 / (m.max(2) as f64 - 1.0) + 1e-4;
    let rho = s
        .correlation
        .max_val(T::from_f64(rho_min))
        .min_val(T::from_f64(0.95));
    let b_orth = (T::one() - rho).max_val(T::from_f64(1e-12)).sqrt();
    let a_comm = (T::one() + rho * (m as f64 - 1.0))
        .max_val(T::from_f64(1e-12))
        .sqrt();
    let temp = s.temperature.max(1e-3);
    let t_end = s.observation_times[m - 1].max(0.0);
    let df = (-s.rate * t_end).exp();
    let rng = ql_mc::PhiloxRng::new(seed);
    let mut sum = T::zero();
    for pth in 0..paths {
        let anti = pth % 2 == 1;
        let id = (pth / 2) as u64;
        let mut cur = s.spots[..m].to_vec();
        let mut avail = vec![T::one(); m];
        let mut locked_sum = T::zero();
        let mut prev_t = 0.0;
        for step in 0..m {
            let t_now = s.observation_times[step].max(prev_t);
            let dt = (t_now - prev_t).max(1e-8);
            let sqrt_dt = dt.sqrt();
            let sgn = if anti { -1.0 } else { 1.0 };
            let mut z_raw = vec![0.0; m];
            let mut z_mean = 0.0;
            for (i, z_slot) in z_raw.iter_mut().enumerate().take(m) {
                let zi = rng.normal4(id, (step * m + i) as u32)[0] * sgn;
                *z_slot = zi;
                z_mean += zi;
            }
            z_mean /= m as f64;
            for (i, &zi) in z_raw.iter().enumerate().take(m) {
                let dw = b_orth * (zi - z_mean) + a_comm * z_mean;
                let v = s.vols[i];
                cur[i] *= ((s.rate - s.div_yields[i] - v * v * 0.5) * dt + v * sqrt_dt * dw).exp();
            }
            let rets: Vec<T> = (0..m).map(|i| cur[i] / s.spots[i] - T::one()).collect();
            let r_max = rets
                .iter()
                .map(|r| r.val())
                .fold(f64::NEG_INFINITY, f64::max);
            let mut weights = vec![T::zero(); m];
            let mut w_sum = T::zero();
            for i in 0..m {
                let w = avail[i] * ((rets[i] - r_max) / temp).exp();
                weights[i] = w;
                w_sum += w;
            }
            let w_denom = w_sum.max_val(T::from_f64(1e-14));
            let mut r_star = T::zero();
            for i in 0..m {
                let p_i = weights[i] / w_denom;
                r_star += p_i * rets[i];
                avail[i] = (avail[i] * (T::one() - p_i)).max_val(T::from_f64(1e-8));
            }
            locked_sum += r_star;
            prev_t = t_now;
        }
        let avg_ret = locked_sum / m as f64;
        let diff = avg_ret - s.strike;
        let payoff = (diff + (diff * diff + 1e-8).sqrt()) * 0.5;
        sum += s.notional * df * payoff;
    }
    sum / paths as f64
}

#[derive(Clone, Debug)]
pub struct NapoleonOptionSpec<T: Real = f64> {
    pub spot: T,
    pub rate: T,
    pub div_yield: T,
    pub vol: T,
    pub reset_times: Vec<f64>,
    pub fixed_coupon: T,
    pub local_floor: T,
    pub local_cap: T,
    pub global_floor: T,
    pub notional: T,
    pub temperature: f64,
}

pub fn price_napoleon_option_mc_ad<T: Real>(
    s: &NapoleonOptionSpec<T>,
    n_paths: usize,
    seed: u64,
) -> T {
    let m = s.reset_times.len().saturating_sub(1);
    if m == 0 {
        return T::zero();
    }
    let paths = n_paths.max(2);
    let temp = s.temperature.max(1e-3);
    let t_end = s.reset_times[m].max(0.0);
    let df = (-s.rate * t_end).exp();
    let rng = ql_mc::PhiloxRng::new(seed);
    let mut sum = T::zero();
    for pth in 0..paths {
        let anti = pth % 2 == 1;
        let id = (pth / 2) as u64;
        let sgn = if anti { -1.0 } else { 1.0 };
        let mut r_vec = Vec::with_capacity(m);
        for j in 0..m {
            let dt = (s.reset_times[j + 1] - s.reset_times[j]).max(1e-8);
            let z = rng.normal4(id, j as u32)[0] * sgn;
            let r_j = ((s.rate - s.div_yield - s.vol * s.vol * 0.5) * dt + s.vol * dt.sqrt() * z)
                .exp()
                - T::one();
            let clipped = r_j.max_val(s.local_floor).min_val(s.local_cap);
            r_vec.push(clipped);
        }
        let r_min_val = r_vec.iter().map(|r| r.val()).fold(f64::INFINITY, f64::min);
        let mut lse = T::zero();
        for &r in &r_vec {
            lse += (-(r - r_min_val) / temp).exp();
        }
        let soft_min = T::from_f64(r_min_val) - (lse / m as f64).ln() * temp;
        let raw_coupon = s.fixed_coupon + soft_min;
        let diff = raw_coupon - s.global_floor;
        let coupon = s.global_floor + (diff + (diff * diff + 1e-8).sqrt()) * 0.5;
        sum += s.notional * df * coupon;
    }
    sum / paths as f64
}

#[derive(Clone, Debug)]
pub struct EverestOptionSpec<T: Real = f64> {
    pub spots: Vec<T>,
    pub vols: Vec<T>,
    pub correlation: T,
    pub rate: T,
    pub expiry: T,
    pub strike: T,
    pub notional: T,
    pub tau_temp: T,
}

pub fn price_everest_option_mc_ad<T: Real>(s: &EverestOptionSpec<T>, paths: usize) -> T {
    let m = s.spots.len().min(s.vols.len());
    if m == 0 {
        return T::zero();
    }
    let n_paths = paths.max(2);
    let t = s.expiry.max_val(T::from_f64(1e-8));
    let sqrt_t = t.sqrt();
    let df = (-s.rate * t).exp();
    let rho_min = -1.0 / (m.max(2) as f64 - 1.0) + 1e-4;
    let rho = s
        .correlation
        .max_val(T::from_f64(rho_min))
        .min_val(T::from_f64(0.99));
    let b_orth = (T::one() - rho).max_val(T::from_f64(1e-12)).sqrt();
    let a_comm = (T::one() + rho * (m as f64 - 1.0))
        .max_val(T::from_f64(1e-12))
        .sqrt();
    let tau = s.tau_temp.max_val(T::from_f64(1e-3));
    let rng = ql_mc::PhiloxRng::new(0xE7E2_E570_u64);
    let mut sum = T::zero();

    for pth in 0..n_paths {
        let anti = pth % 2 == 1;
        let id = (pth / 2) as u64;
        let sgn = if anti { -1.0 } else { 1.0 };
        let mut z_raw = vec![0.0; m];
        let mut z_mean = 0.0;
        for (i, slot) in z_raw.iter_mut().enumerate().take(m) {
            let zi = rng.normal4(id, i as u32)[0] * sgn;
            *slot = zi;
            z_mean += zi;
        }
        z_mean /= m as f64;
        let mut rets = Vec::with_capacity(m);
        let mut r_min = f64::INFINITY;
        for (i, &zi) in z_raw.iter().enumerate().take(m) {
            let dw = b_orth * (zi - z_mean) + a_comm * z_mean;
            let v = s.vols[i];
            let r_i = (s.spots[i] / s.spots[0].max_val(T::from_f64(1e-12)))
                * ((s.rate - v * v * 0.5) * t + v * sqrt_t * dw).exp();
            r_min = r_min.min(r_i.val());
            rets.push(r_i);
        }
        let mut w_sum = T::zero();
        let mut rw_sum = T::zero();
        for &r_i in &rets {
            let w_i = (-(r_i - r_min) / tau).exp();
            w_sum += w_i;
            rw_sum += w_i * r_i;
        }
        let r_worst = rw_sum / w_sum.max_val(T::from_f64(1e-12));
        let diff = r_worst - s.strike;
        let pay = (diff + (diff * diff + T::from_f64(1e-8)).sqrt()) * 0.5;
        sum += s.notional * df * pay;
    }
    sum / (n_paths as f64)
}

#[cfg(test)]
mod phase16_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn himalaya_and_napoleon_ad_and_short_vega() {
        let nap = NapoleonOptionSpec {
            spot: 100.0,
            rate: 0.02,
            div_yield: 0.01,
            vol: 0.15,
            reset_times: vec![0.0, 0.25, 0.5, 0.75, 1.0],
            fixed_coupon: 0.08,
            local_floor: -0.15,
            local_cap: 0.10,
            global_floor: 0.0,
            notional: 100.0,
            temperature: 0.01,
        };
        let p_low = price_napoleon_option_mc_ad(&nap, 128, 19);
        let mut nap_hi = nap.clone();
        nap_hi.vol = 0.35;
        let p_hi = price_napoleon_option_mc_ad(&nap_hi, 128, 19);
        assert!(p_low > p_hi);
        let him = HimalayaOptionSpec {
            spots: vec![100.0, 100.0, 100.0],
            vols: vec![0.2, 0.22, 0.25],
            div_yields: vec![0.01, 0.01, 0.01],
            correlation: 0.3,
            rate: 0.02,
            observation_times: vec![0.33, 0.67, 1.0],
            strike: 0.0,
            notional: 100.0,
            temperature: 0.02,
        };
        let h = 1e-5;
        let eval_h = |v0: f64| {
            let mut sp = him.clone();
            sp.vols[0] = v0;
            price_himalaya_option_mc_ad(&sp, 64, 11)
        };
        let fd = (eval_h(0.2 + h) - eval_h(0.2 - h)) / (2.0 * h);
        type D = Dual<1>;
        let hd = HimalayaOptionSpec {
            spots: vec![D::from_value(100.0); 3],
            vols: vec![
                D::variable(0.2, 0),
                D::from_value(0.22),
                D::from_value(0.25),
            ],
            div_yields: vec![D::from_value(0.01); 3],
            correlation: D::from_value(0.3),
            rate: D::from_value(0.02),
            observation_times: vec![0.33, 0.67, 1.0],
            strike: D::zero(),
            notional: D::from_value(100.0),
            temperature: 0.02,
        };
        let ad = price_himalaya_option_mc_ad(&hd, 64, 11);
        assert!((ad.d[0] - fd).abs() < 1e-4);
    }
}

#[derive(Clone, Debug)]
pub struct AtlasOptionSpec<T: Real = f64> {
    pub spots: Vec<T>,
    pub vols: Vec<T>,
    pub div_yields: Vec<T>,
    pub correlation: T,
    pub rate: T,
    pub expiry: f64,
    pub strike: T,
    pub notional: T,
    pub remove_best: usize,
    pub remove_worst: usize,
}

pub fn price_atlas_option_mc_ad<T: Real>(
    s: &AtlasOptionSpec<T>,
    n_paths: usize,
    seed: u64,
) -> T {
    let m = s.spots.len().min(s.vols.len()).min(s.div_yields.len());
    if m == 0 || s.remove_best + s.remove_worst >= m {
        return T::zero();
    }
    let kept = m - s.remove_best - s.remove_worst;
    let paths = n_paths.max(2);
    let t = s.expiry.max(1e-8);
    let sqrt_t = t.sqrt();
    let df = (-s.rate * t).exp();
    let rho_min = -1.0 / (m.max(2) as f64 - 1.0) + 1e-4;
    let rho = s
        .correlation
        .max_val(T::from_f64(rho_min))
        .min_val(T::from_f64(0.95));
    let b_orth = (T::one() - rho).max_val(T::from_f64(1e-12)).sqrt();
    let a_comm = (T::one() + rho * (m as f64 - 1.0))
        .max_val(T::from_f64(1e-12))
        .sqrt();

    let rng = ql_mc::PhiloxRng::new(seed);
    let mut sum = T::zero();
    for pth in 0..paths {
        let anti = pth % 2 == 1;
        let id = (pth / 2) as u64;
        let sgn = if anti { -1.0 } else { 1.0 };
        let mut z_raw = vec![0.0; m];
        let mut z_mean = 0.0;
        for (i, slot) in z_raw.iter_mut().enumerate().take(m) {
            let zi = rng.normal4(id, i as u32)[0] * sgn;
            *slot = zi;
            z_mean += zi;
        }
        z_mean /= m as f64;
        let mut perf: Vec<(f64, T)> = Vec::with_capacity(m);
        for (i, &zi) in z_raw.iter().enumerate().take(m) {
            let dw = b_orth * (zi - z_mean) + a_comm * z_mean;
            let v = s.vols[i];
            let ret_i = ((s.rate - s.div_yields[i] - v * v * 0.5) * t + v * sqrt_t * dw).exp()
                - T::one();
            perf.push((ret_i.val(), ret_i));
        }
        perf.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut avg_kept = T::zero();
        for &(_, r_t) in perf.iter().skip(s.remove_worst).take(kept) {
            avg_kept += r_t;
        }
        avg_kept /= T::from_f64(kept as f64);
        let diff = avg_kept - s.strike;
        let payoff = (diff + (diff * diff + 1e-8).sqrt()) * 0.5;
        sum += s.notional * df * payoff;
    }
    sum / paths as f64
}

#[derive(Clone, Debug)]
pub struct AltiplanoOptionSpec<T: Real = f64> {
    pub spots: Vec<T>,
    pub vols: Vec<T>,
    pub div_yields: Vec<T>,
    pub correlation: T,
    pub rate: T,
    pub observation_times: Vec<f64>,
    pub coupon_barrier_ratio: T,
    pub fixed_coupon: T,
    pub vanilla_strike: T,
    pub notional: T,
    pub smoothing_eps: f64,
}

pub fn price_altiplano_option_mc_ad<T: Real>(
    s: &AltiplanoOptionSpec<T>,
    n_paths: usize,
    seed: u64,
) -> T {
    let m = s.spots.len().min(s.vols.len()).min(s.div_yields.len());
    if m == 0 || s.observation_times.is_empty() {
        return T::zero();
    }
    let paths = n_paths.max(2);
    let rho_min = -1.0 / (m.max(2) as f64 - 1.0) + 1e-4;
    let rho = s
        .correlation
        .max_val(T::from_f64(rho_min))
        .min_val(T::from_f64(0.95));
    let b_orth = (T::one() - rho).max_val(T::from_f64(1e-12)).sqrt();
    let a_comm = (T::one() + rho * (m as f64 - 1.0))
        .max_val(T::from_f64(1e-12))
        .sqrt();
    let eps = T::from_f64(s.smoothing_eps.max(1e-3));
    let t_end = s.observation_times.last().copied().unwrap_or(1.0).max(0.0);
    let df = (-s.rate * t_end).exp();
    let rng = ql_mc::PhiloxRng::new(seed);
    let mut sum = T::zero();

    for pth in 0..paths {
        let anti = pth % 2 == 1;
        let id = (pth / 2) as u64;
        let sgn = if anti { -1.0 } else { 1.0 };
        let mut ratios = vec![T::one(); m];
        let mut no_breach = T::one();
        let mut prev_t = 0.0;
        for (step, &t_now) in s.observation_times.iter().enumerate() {
            let dt = (t_now - prev_t).max(1e-8);
            let sqrt_dt = dt.sqrt();
            let mut z_raw = vec![0.0; m];
            let mut z_mean = 0.0;
            for (i, slot) in z_raw.iter_mut().enumerate().take(m) {
                let zi = rng.normal4(id, (step * m + i) as u32)[0] * sgn;
                *slot = zi;
                z_mean += zi;
            }
            z_mean /= m as f64;
            for (i, &zi) in z_raw.iter().enumerate().take(m) {
                let dw = b_orth * (zi - z_mean) + a_comm * z_mean;
                let v = s.vols[i];
                ratios[i] *= ((s.rate - s.div_yields[i] - v * v * 0.5) * dt + v * sqrt_dt * dw).exp();
                let dist = ratios[i] - s.coupon_barrier_ratio;
                let w_ok = T::from_f64(0.5)
                    * (T::one() + dist / (dist * dist + eps * eps).sqrt());
                no_breach *= w_ok;
            }
            prev_t = t_now.max(prev_t);
        }
        let avg_ratio = ratios.iter().copied().fold(T::zero(), |a, b| a + b) / (m as f64);
        let diff = avg_ratio - s.vanilla_strike;
        let call_payoff = (diff + (diff * diff + 1e-8).sqrt()) * 0.5;
        let payoff = no_breach * s.fixed_coupon + (T::one() - no_breach) * call_payoff;
        sum += s.notional * df * payoff;
    }
    sum / paths as f64
}

#[cfg(test)]
mod phase17_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn atlas_and_altiplano_monotonicity_and_ad() {
        let mut atlas = AtlasOptionSpec {
            spots: vec![100.0, 100.0, 100.0, 100.0],
            vols: vec![0.20, 0.22, 0.25, 0.28],
            div_yields: vec![0.01; 4],
            correlation: 0.25,
            rate: 0.03,
            expiry: 1.0,
            strike: 0.0,
            notional: 100.0,
            remove_best: 1,
            remove_worst: 0,
        };
        let p_drop_best = price_atlas_option_mc_ad(&atlas, 64, 23);
        atlas.remove_best = 0;
        atlas.remove_worst = 1;
        let p_drop_worst = price_atlas_option_mc_ad(&atlas, 64, 23);
        assert!(p_drop_worst > p_drop_best);

        let alt = AltiplanoOptionSpec {
            spots: vec![100.0, 100.0, 100.0],
            vols: vec![0.18, 0.20, 0.22],
            div_yields: vec![0.01; 3],
            correlation: 0.30,
            rate: 0.02,
            observation_times: vec![0.5, 1.0],
            coupon_barrier_ratio: 0.85,
            fixed_coupon: 0.12,
            vanilla_strike: 1.0,
            notional: 100.0,
            smoothing_eps: 0.02,
        };
        let h = 1e-5;
        let eval_b = |b_rat: f64| {
            let mut sp = alt.clone();
            sp.coupon_barrier_ratio = b_rat;
            price_altiplano_option_mc_ad(&sp, 64, 31)
        };
        let fd = (eval_b(0.85 + h) - eval_b(0.85 - h)) / (2.0 * h);
        type D = Dual<1>;
        let alt_d = AltiplanoOptionSpec {
            spots: vec![D::from_value(100.0); 3],
            vols: vec![
                D::from_value(0.18),
                D::from_value(0.20),
                D::from_value(0.22),
            ],
            div_yields: vec![D::from_value(0.01); 3],
            correlation: D::from_value(0.30),
            rate: D::from_value(0.02),
            observation_times: vec![0.5, 1.0],
            coupon_barrier_ratio: D::variable(0.85, 0),
            fixed_coupon: D::from_value(0.12),
            vanilla_strike: D::from_value(1.0),
            notional: D::from_value(100.0),
            smoothing_eps: 0.02,
        };
        let ad = price_altiplano_option_mc_ad(&alt_d, 64, 31);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn everest_correlation_monotonicity_and_ad() {
        let mut s = EverestOptionSpec {
            spots: vec![100.0, 100.0, 100.0],
            vols: vec![0.20, 0.25, 0.22],
            correlation: 0.1,
            rate: 0.03,
            expiry: 1.0,
            strike: 0.90,
            notional: 100.0,
            tau_temp: 0.02,
        };
        let p_low_corr = price_everest_option_mc_ad(&s, 64);
        s.correlation = 0.85;
        let p_high_corr = price_everest_option_mc_ad(&s, 64);
        assert!(p_high_corr > p_low_corr);

        let h = 1e-5;
        let eval_rho = |r: f64| {
            let mut sp = s.clone();
            sp.correlation = r;
            price_everest_option_mc_ad(&sp, 64)
        };
        let fd = (eval_rho(0.5 + h) - eval_rho(0.5 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sd = EverestOptionSpec {
            spots: vec![D::from_value(100.0); 3],
            vols: vec![D::from_value(0.20), D::from_value(0.25), D::from_value(0.22)],
            correlation: D::variable(0.5, 0),
            rate: D::from_value(0.03),
            expiry: D::one(),
            strike: D::from_value(0.90),
            notional: D::from_value(100.0),
            tau_temp: D::from_value(0.02),
        };
        let ad = price_everest_option_mc_ad(&sd, 64);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[derive(Clone, Debug)]
pub struct BouzoubaaReverseAndGlobalCliquetSpec<T: Real = f64>{pub notional:T,pub rate:T,pub div_yield:T,pub reset_times:Vec<f64>,pub forward_vols:Vec<T>,pub forward_skews:Vec<T>,pub initial_coupon_cmax:T,pub local_floor:T,pub local_cap:T,pub global_floor:T,pub global_cap:T,pub is_reverse_cliquet:bool}
#[derive(Copy, Clone, Debug)]
pub struct BouzoubaaCliquetValuation<T: Real = f64>{pub npv:T,pub sum_of_local_caps_npv:T,pub global_floor_cap_optionality:T,pub total_return_mean:T,pub total_return_std:T,pub total_return_skewness:T,pub total_return_excess_kurtosis:T,pub forward_skew_adjustment_pv:T}
pub fn price_bouzoubaa_reverse_and_global_cliquet<T: Real>(s:&BouzoubaaReverseAndGlobalCliquetSpec<T>)->BouzoubaaCliquetValuation<T>{
    let nodes=[-2.532731674, -1.756683649, -1.036610829, -0.342901327,0.342901327,1.036610829,1.756683649,2.532731674];let weights=[0.011257411,0.222075923,0.625704404,0.813858273,0.813858273,0.625704404,0.222075923,0.011257411];let m=s.reset_times.len();let mut means=Vec::new();let mut vars=Vec::new();let mut k3=T::zero();let mut k4=T::zero();let mut mean=T::zero();let mut localpv=T::zero();let mut skewpv=T::zero();let mut prev=0.0;for i in 0..m{let dt=(s.reset_times[i]-prev).max(0.0);prev=s.reset_times[i];let vol=s.forward_vols.get(i).copied().unwrap_or(T::zero());let sk=s.forward_skews.get(i).copied().unwrap_or(T::zero());let mu=(s.rate-s.div_yield-vol*vol*0.5)*dt as f64;let vv=vol*dt.sqrt();let mut ys=Vec::new();for (&z,&w) in nodes.iter().zip(weights.iter()){let r=(mu+vv*z-sk*vv*vv*0.5*(z*z-1.0)).exp()-T::one();let y=if s.is_reverse_cliquet{T::zero().min_val(r.max_val(s.local_floor))}else{r.max_val(s.local_floor).min_val(s.local_cap)};ys.push(y);mean+=T::from_f64(w/2.0_f64.sqrt())*y;}let mm=ys.iter().copied().fold(T::zero(),|a,b|a+b)/ys.len() as f64;let mut v=T::zero();for y in &ys{let d=*y-mm;v+=d*d;}v*=T::from_f64(1.0/ys.len() as f64);vars.push(v);means.push(mm);for y in &ys{let d=*y-mm;k3+=d*d*d/ys.len() as f64;k4+=d*d*d*d/ys.len() as f64;}localpv+=mm;skewpv+=sk*vv;}
    let mu=s.initial_coupon_cmax+means.iter().copied().fold(T::zero(),|a,b|a+b);let var=vars.iter().copied().fold(T::zero(),|a,b|a+b);let sd=var.max_val(T::from_f64(1e-12)).sqrt();let skew=k3/sd.powf(3.0);let kurt=k4/sd.powf(4.0)-3.0;let call=|k:T|{let d=(mu-k)/sd; (mu-k)*d.norm_cdf()+sd*d.norm_pdf()*(T::one()+skew*d/6.0+kurt*(d*d-T::one())/24.0)};let val=(s.global_floor+call(s.global_floor)-call(s.global_cap)).max_val(s.global_floor).min_val(s.global_cap);let df=(-s.rate*s.reset_times.last().copied().unwrap_or(0.0)).exp();let npv=s.notional*df*val;let local=s.notional*df*mu;BouzoubaaCliquetValuation{npv,sum_of_local_caps_npv:local,global_floor_cap_optionality:npv-local,total_return_mean:mu,total_return_std:sd,total_return_skewness:skew,total_return_excess_kurtosis:kurt,forward_skew_adjustment_pv:s.notional*df*skewpv}
}
