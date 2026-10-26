use crate::barrier::{BarrierOptionSpec, MonitoringFrequency, price_barrier_option};
use ql_ad::Real;
use ql_math::{bivariate_normal_cdf, norm_cdf};
use ql_vol::black_scholes_price;

#[derive(Clone, Copy, Debug)]
pub struct ChesneyYorCumulativeParisianRainbowSpec<T: Real = f64> {
    pub spot_1: T,
    pub spot_2: T,
    pub strike: T,
    pub par_barrier: T,
    pub excursion_window_d: T,
    pub vol_1: T,
    pub vol_2: T,
    pub div_1: T,
    pub div_2: T,
    pub correlation: T,
    pub risk_free_rate: T,
    pub expiry: T,
    pub is_up_and_out: bool,
}
#[derive(Clone, Debug)]
pub struct ChesneyYorCumulativeParisianRainbowReport<T: Real = f64> {
    pub vanilla_call_1: T,
    pub vanilla_call_2: T,
    pub instant_ko_barrier_price: T,
    pub consecutive_parisian_ko_price: T,
    pub cumulative_parisian_ko_price: T,
    pub parisian_protection_premium: T,
    pub rainbow_call_on_max: T,
    pub rainbow_call_on_min: T,
    pub rainbow_put_on_min: T,
    pub stulz_rainbow_parity_error: T,
}
pub fn price_chesney_yor_cumulative_parisian_and_rainbow<T: Real>(
    s: &ChesneyYorCumulativeParisianRainbowSpec<T>,
) -> ChesneyYorCumulativeParisianRainbowReport<T> {
    let t = s.expiry.max_val(T::from_f64(1e-8));
    let c1 = black_scholes_price(
        s.spot_1,
        s.strike,
        s.vol_1,
        t,
        s.risk_free_rate,
        s.div_1,
        true,
    );
    let c2 = black_scholes_price(
        s.spot_2,
        s.strike,
        s.vol_2,
        t,
        s.risk_free_rate,
        s.div_2,
        true,
    );
    let bs = BarrierOptionSpec {
        option_type: ql_instruments::OptionType::Call,
        barrier_type: if s.is_up_and_out {
            ql_instruments::BarrierType::UpOut
        } else {
            ql_instruments::BarrierType::DownOut
        },
        spot: s.spot_1,
        strike: s.strike,
        barrier: s.par_barrier,
        rebate: T::zero(),
        rate: s.risk_free_rate,
        div_yield: s.div_1,
        vol: s.vol_1,
        expiry: t,
        monitoring: MonitoringFrequency::Continuous,
    };
    let ko = price_barrier_option(&bs);
    let win = (T::one() - (s.excursion_window_d / t).powf(2.0)).max_val(T::zero());
    let kc = T::from_f64(std::f64::consts::PI / 2.0).sqrt() * s.vol_1 * s.excursion_window_d.sqrt();
    let km = T::from_f64(0.5).sqrt() * s.vol_1 * s.excursion_window_d.sqrt();
    let shifted_c = if s.is_up_and_out {
        s.par_barrier * kc.exp()
    } else {
        s.par_barrier * (-kc).exp()
    };
    let shifted_m = if s.is_up_and_out {
        s.par_barrier * km.exp()
    } else {
        s.par_barrier * (-km).exp()
    };
    let bsc = BarrierOptionSpec {
        barrier: shifted_c,
        ..bs
    };
    let bsm = BarrierOptionSpec {
        barrier: shifted_m,
        ..bs
    };
    let consec = c1 - win * (c1 - price_barrier_option(&bsc));
    let cumul = c1 - win * (c1 - price_barrier_option(&bsm));
    let sx = (s.vol_1 * s.vol_1 + s.vol_2 * s.vol_2
        - T::from_f64(2.0) * s.correlation * s.vol_1 * s.vol_2)
        .max_val(T::from_f64(1e-10))
        .sqrt();
    let r1 = (s.vol_1 - s.correlation * s.vol_2) / sx;
    let r2 = (s.vol_2 - s.correlation * s.vol_1) / sx;
    let g1 = ((s.spot_1 / s.strike).ln()
        + (s.risk_free_rate - s.div_1 + T::from_f64(0.5) * s.vol_1 * s.vol_1) * t)
        / (s.vol_1 * t.sqrt());
    let g2 = ((s.spot_2 / s.strike).ln()
        + (s.risk_free_rate - s.div_2 + T::from_f64(0.5) * s.vol_2 * s.vol_2) * t)
        / (s.vol_2 * t.sqrt());
    let g12 = ((s.spot_1 / s.spot_2).ln() + (s.div_2 - s.div_1 - T::from_f64(0.5) * sx * sx) * t)
        / (sx * t.sqrt());
    let g21 = ((s.spot_2 / s.spot_1).ln() + (s.div_1 - s.div_2 - T::from_f64(0.5) * sx * sx) * t)
        / (sx * t.sqrt());
    let minc = s.spot_1 * (-s.div_1 * t).exp() * bivariate_normal_cdf(g1, g21, -r1)
        + s.spot_2 * (-s.div_2 * t).exp() * bivariate_normal_cdf(g2, g12, -r2)
        - s.strike
            * (-s.risk_free_rate * t).exp()
            * bivariate_normal_cdf(
                g1 - s.vol_1 * t.sqrt(),
                g2 - s.vol_2 * t.sqrt(),
                s.correlation,
            );
    let maxc_stulz = s.spot_1 * (-s.div_1 * t).exp() * bivariate_normal_cdf(g1, -g21, r1)
        + s.spot_2 * (-s.div_2 * t).exp() * bivariate_normal_cdf(g2, -g12, r2)
        - s.strike
            * (-s.risk_free_rate * t).exp()
            * (T::one()
                - bivariate_normal_cdf(
                    -g1 + s.vol_1 * t.sqrt(),
                    -g2 + s.vol_2 * t.sqrt(),
                    s.correlation,
                ));
    let maxc = maxc_stulz;
    let f_min = s.spot_1 * (-s.div_1 * t).exp() * g21.norm_cdf()
        + s.spot_2 * (-s.div_2 * t).exp() * g12.norm_cdf();
    let putmin = (minc + s.strike * (-s.risk_free_rate * t).exp() - f_min).max_val(T::zero());
    ChesneyYorCumulativeParisianRainbowReport {
        vanilla_call_1: c1,
        vanilla_call_2: c2,
        instant_ko_barrier_price: ko,
        consecutive_parisian_ko_price: consec.max_val(ko).min_val(c1),
        cumulative_parisian_ko_price: cumul.max_val(ko).min_val(consec).min_val(c1),
        parisian_protection_premium: (consec - ko).max_val(T::zero()),
        rainbow_call_on_max: maxc.max_val(T::zero()),
        rainbow_call_on_min: minc.max_val(T::zero()),
        rainbow_put_on_min: putmin,
        stulz_rainbow_parity_error: (maxc + minc - c1 - c2).abs(),
    }
}

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
mod domain_suite_28_tests {
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

pub fn price_atlas_option_mc_ad<T: Real>(s: &AtlasOptionSpec<T>, n_paths: usize, seed: u64) -> T {
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
            let ret_i =
                ((s.rate - s.div_yields[i] - v * v * 0.5) * t + v * sqrt_t * dw).exp() - T::one();
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
                ratios[i] *=
                    ((s.rate - s.div_yields[i] - v * v * 0.5) * dt + v * sqrt_dt * dw).exp();
                let dist = ratios[i] - s.coupon_barrier_ratio;
                let w_ok = T::from_f64(0.5) * (T::one() + dist / (dist * dist + eps * eps).sqrt());
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
mod domain_suite_29_tests {
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

#[derive(Clone, Debug)]
pub struct BouzoubaaReverseAndGlobalCliquetSpec<T: Real = f64> {
    pub notional: T,
    pub rate: T,
    pub div_yield: T,
    pub reset_times: Vec<f64>,
    pub forward_vols: Vec<T>,
    pub forward_skews: Vec<T>,
    pub initial_coupon_cmax: T,
    pub local_floor: T,
    pub local_cap: T,
    pub global_floor: T,
    pub global_cap: T,
    pub is_reverse_cliquet: bool,
}
#[derive(Copy, Clone, Debug)]
pub struct BouzoubaaCliquetValuation<T: Real = f64> {
    pub npv: T,
    pub sum_of_local_caps_npv: T,
    pub global_floor_cap_optionality: T,
    pub total_return_mean: T,
    pub total_return_std: T,
    pub total_return_skewness: T,
    pub total_return_excess_kurtosis: T,
    pub forward_skew_adjustment_pv: T,
}
pub fn price_bouzoubaa_reverse_and_global_cliquet<T: Real>(
    s: &BouzoubaaReverseAndGlobalCliquetSpec<T>,
) -> BouzoubaaCliquetValuation<T> {
    let nodes = [
        -2.532731674,
        -1.756683649,
        -1.036610829,
        -0.342901327,
        0.342901327,
        1.036610829,
        1.756683649,
        2.532731674,
    ];
    let weights = [
        0.011257411,
        0.222075923,
        0.625704404,
        0.813858273,
        0.813858273,
        0.625704404,
        0.222075923,
        0.011257411,
    ];
    let eval_cliquet_core = |with_skew: bool| {
        let mut mean = T::zero();
        let mut variance = T::zero();
        let mut k3 = T::zero();
        let mut k4 = T::zero();
        let mut previous = 0.0;
        for i in 0..s.reset_times.len() {
            let dt = (s.reset_times[i] - previous).max(0.0);
            previous = s.reset_times[i];
            let vol = s.forward_vols.get(i).copied().unwrap_or(T::zero());
            let skew = if with_skew {
                s.forward_skews.get(i).copied().unwrap_or(T::zero())
            } else {
                T::zero()
            };
            let mu = (s.rate - s.div_yield - vol * vol * 0.5) * dt;
            let sd = vol * T::from_f64(dt.sqrt());
            let mut ys = [T::zero(); 8];
            let mut period_mean = T::zero();
            for (q, (&x, &weight)) in nodes.iter().zip(weights.iter()).enumerate() {
                let z = T::from_f64(2.0_f64.sqrt() * x);
                let y = (mu + sd * z - skew * sd * sd * 0.5 * (z * z - T::one())).exp() - T::one();
                let clipped = if s.is_reverse_cliquet {
                    T::zero().min_val(y.max_val(s.local_floor))
                } else {
                    y.max_val(s.local_floor).min_val(s.local_cap)
                };
                ys[q] = clipped;
                period_mean += T::from_f64(weight / std::f64::consts::PI.sqrt()) * clipped;
            }
            let mut period_var = T::zero();
            let mut period_k3 = T::zero();
            let mut period_k4 = T::zero();
            for (q, _) in nodes.iter().enumerate() {
                let p = T::from_f64(weights[q] / std::f64::consts::PI.sqrt());
                let d = ys[q] - period_mean;
                period_var += p * d * d;
                period_k3 += p * d * d * d;
                period_k4 += p * d * d * d * d;
            }
            mean += period_mean;
            variance += period_var;
            k3 += period_k3;
            k4 += period_k4 - T::from_f64(3.0) * period_var * period_var;
        }
        let total_mean = s.initial_coupon_cmax + mean;
        let total_var = variance.max_val(T::from_f64(1e-14));
        let total_sd = total_var.sqrt();
        let total_skew = k3 / total_sd.powf(3.0);
        let total_kurt = k4 / (total_var * total_var);
        let call = |k: T| {
            let d = (total_mean - k) / total_sd;
            (total_mean - k) * d.norm_cdf()
                + total_sd
                    * d.norm_pdf()
                    * (T::one() + total_skew * d / 6.0 + total_kurt * (d * d - T::one()) / 24.0)
        };
        let value = (s.global_floor + call(s.global_floor) - call(s.global_cap))
            .max_val(s.global_floor)
            .min_val(s.global_cap);
        let df = (-s.rate * s.reset_times.last().copied().unwrap_or(0.0)).exp();
        let npv = s.notional * df * value;
        (npv, total_mean, total_sd, total_skew, total_kurt)
    };
    let (npv, mean, sd, skew, kurt) = eval_cliquet_core(true);
    let zero_skew_npv = eval_cliquet_core(false).0;
    let df = (-s.rate * s.reset_times.last().copied().unwrap_or(0.0)).exp();
    BouzoubaaCliquetValuation {
        npv,
        sum_of_local_caps_npv: s.notional * df * mean,
        global_floor_cap_optionality: npv - s.notional * df * mean,
        total_return_mean: mean,
        total_return_std: sd,
        total_return_skewness: skew,
        total_return_excess_kurtosis: kurt,
        forward_skew_adjustment_pv: npv - zero_skew_npv,
    }
}

#[cfg(test)]
mod domain_suite_30_tests {
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
            vols: vec![
                D::from_value(0.20),
                D::from_value(0.25),
                D::from_value(0.22),
            ],
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

#[cfg(test)]
mod domain_suite_33_tests {
    use super::*;

    #[test]
    fn bouzoubaa_cumulants_and_skew_adjustment() {
        let s = BouzoubaaReverseAndGlobalCliquetSpec {
            notional: 100.0,
            rate: 0.03,
            div_yield: 0.0,
            reset_times: vec![0.5, 1.0],
            forward_vols: vec![0.2, 0.2],
            forward_skews: vec![0.1, 0.1],
            initial_coupon_cmax: 0.0,
            local_floor: -0.1,
            local_cap: 0.1,
            global_floor: -0.2,
            global_cap: 0.2,
            is_reverse_cliquet: false,
        };
        let r = price_bouzoubaa_reverse_and_global_cliquet(&s);
        assert!(r.total_return_std.val() >= 0.0);
        assert!(r.forward_skew_adjustment_pv.val().is_finite());
    }
}

#[cfg(test)]
mod domain_suite_35_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn parisian_order() {
        let s = ChesneyYorCumulativeParisianRainbowSpec {
            spot_1: 100.0,
            spot_2: 102.0,
            strike: 100.0,
            par_barrier: 120.0,
            excursion_window_d: 0.1,
            vol_1: 0.2,
            vol_2: 0.2,
            div_1: 0.0,
            div_2: 0.0,
            correlation: 0.5,
            risk_free_rate: 0.02,
            expiry: 1.0,
            is_up_and_out: true,
        };
        let r = price_chesney_yor_cumulative_parisian_and_rainbow(&s);
        assert!(r.instant_ko_barrier_price <= r.consecutive_parisian_ko_price);
        assert!(r.stulz_rainbow_parity_error < 1e-12);
        assert!(r.rainbow_put_on_min > 0.0);
    }

    #[test]
    fn dual_fd_parity() {
        fn mk<T: Real>(x: T) -> ChesneyYorCumulativeParisianRainbowSpec<T> {
            ChesneyYorCumulativeParisianRainbowSpec {
                spot_1: T::from_f64(100.0),
                spot_2: x,
                strike: T::from_f64(100.0),
                par_barrier: T::from_f64(120.0),
                excursion_window_d: T::from_f64(0.1),
                vol_1: T::from_f64(0.2),
                vol_2: T::from_f64(0.2),
                div_1: T::zero(),
                div_2: T::zero(),
                correlation: T::from_f64(0.5),
                risk_free_rate: T::from_f64(0.02),
                expiry: T::one(),
                is_up_and_out: true,
            }
        }
        let h = 1e-5;
        let fd = (price_chesney_yor_cumulative_parisian_and_rainbow(&mk(102.0 + h))
            .rainbow_call_on_max
            .val()
            - price_chesney_yor_cumulative_parisian_and_rainbow(&mk(102.0 - h))
                .rainbow_call_on_max
                .val())
            / (2.0 * h);
        let ad =
            price_chesney_yor_cumulative_parisian_and_rainbow(&mk(Dual::<1>::variable(102.0, 0)))
                .rainbow_call_on_max;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

mod evaluate_bouzoubaa_napoleon_and_reverse_cliquet_module {
    use ql_ad::Real;

    #[inline]
    fn smin<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b - (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[inline]
    fn smax<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[derive(Clone, Debug)]
    pub struct BouzoubaaNapoleonReverseCliquetSpec<T: Real = f64> {
        pub period_returns: Vec<T>,
        pub napoleon_fixed_coupon: T,
        pub global_floor: T,
        pub local_cap: T,
        pub local_floor: T,
        pub forward_volatility: T,
        pub forward_skew: T,
        pub risk_free_rate: T,
        pub maturity_years: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct BouzoubaaNapoleonReverseCliquetReport<T: Real = f64> {
        pub worst_negative_period_return: T,
        pub cumulative_negative_return_sum: T,
        pub napoleon_coupon_payoff: T,
        pub reverse_cliquet_coupon_payoff: T,
        pub standard_capped_cliquet_payoff: T,
        pub forward_vol_cliquet_convexity_cost: T,
        pub napoleon_capital_protected_note_pv: T,
    }

    pub fn evaluate_bouzoubaa_napoleon_and_reverse_cliquet<T: Real>(
        s: &BouzoubaaNapoleonReverseCliquetSpec<T>,
    ) -> BouzoubaaNapoleonReverseCliquetReport<T> {
        let mut w_min = T::zero();
        let mut neg_sum = T::zero();
        let mut std_sum = T::zero();
        for &rm in &s.period_returns {
            let neg_m = smin(rm, T::zero());
            w_min = smin(w_min, neg_m);
            neg_sum += neg_m;
            std_sum += smax(s.local_floor, smin(s.local_cap, rm));
        }
        let c_nap = smax(s.global_floor, s.napoleon_fixed_coupon + w_min);
        let c_rev = smax(s.global_floor, s.napoleon_fixed_coupon + neg_sum);
        let c_std = smax(s.global_floor, std_sum);
        let m = (s.period_returns.len().max(1)) as f64;
        let cost = T::from_f64(0.5)
            * s.forward_volatility
            * (s.maturity_years / T::from_f64(m)).sqrt()
            * T::from_f64((m + 1.0).ln())
            * (T::one() - s.forward_skew);
        let pv = (-s.risk_free_rate * s.maturity_years).exp() * (T::one() + c_nap - cost);
        BouzoubaaNapoleonReverseCliquetReport {
            worst_negative_period_return: w_min,
            cumulative_negative_return_sum: neg_sum,
            napoleon_coupon_payoff: c_nap,
            reverse_cliquet_coupon_payoff: c_rev,
            standard_capped_cliquet_payoff: c_std,
            forward_vol_cliquet_convexity_cost: cost,
            napoleon_capital_protected_note_pv: pv,
        }
    }

    #[cfg(test)]
    mod domain_suite_47_tests {
        use super::*;
        use ql_ad::Dual;
        fn spec<T: Real>(x: T) -> BouzoubaaNapoleonReverseCliquetSpec<T> {
            BouzoubaaNapoleonReverseCliquetSpec {
                period_returns: vec![
                    T::from_f64(0.03),
                    T::from_f64(-0.04) + x,
                    T::from_f64(0.02),
                    T::from_f64(-0.025),
                ],
                napoleon_fixed_coupon: T::from_f64(0.08),
                global_floor: T::zero(),
                local_cap: T::from_f64(0.04),
                local_floor: T::from_f64(-0.04),
                forward_volatility: T::from_f64(0.20),
                forward_skew: T::from_f64(-0.15),
                risk_free_rate: T::from_f64(0.03),
                maturity_years: T::from_f64(2.0),
            }
        }
        #[test]
        fn suite47_exotics_invariants_and_dual_fd() {
            let r = evaluate_bouzoubaa_napoleon_and_reverse_cliquet(&spec(0.0));
            assert!(r.napoleon_coupon_payoff >= 0.0);
            assert!(r.napoleon_coupon_payoff >= r.reverse_cliquet_coupon_payoff - 1e-8);
            let h = 1e-6;
            let fd = (evaluate_bouzoubaa_napoleon_and_reverse_cliquet(&spec(h))
                .napoleon_capital_protected_note_pv
                - evaluate_bouzoubaa_napoleon_and_reverse_cliquet(&spec(-h))
                    .napoleon_capital_protected_note_pv)
                / (2.0 * h);
            let ad =
                evaluate_bouzoubaa_napoleon_and_reverse_cliquet(&spec(Dual::<1>::variable(0.0, 0)))
                    .napoleon_capital_protected_note_pv;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_bouzoubaa_napoleon_and_reverse_cliquet_module::*;

mod evaluate_haug_eso_vesting_barrier_dilution_module {
    // Haug (2007) & Hull–White (2004) Executive Stock Option (`ESO`) Vesting-Forfeiture,
    // Suboptimal Exercise Multiple ($M$-Barrier) & Warrant Dilution Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct HaugEsoVestingBarrierDilutionSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub strike_price_k: T,
        pub risk_free_rate_r: T,
        pub dividend_yield_q: T,
        pub volatility_sigma: T,
        pub contractual_maturity_t: T,
        pub vesting_period_tv: T,
        pub forfeiture_hazard_lambda_f: T,
        pub suboptimal_exercise_multiple_m: T,
        pub grant_to_shares_ratio_gamma: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct HaugEsoVestingBarrierDilutionReport<T: Real = f64> {
        pub vesting_survival_probability: T,
        pub effective_expected_life_teff: T,
        pub plain_black_scholes_call: T,
        pub forfeiture_adjusted_eso_value: T,
        pub diluted_per_option_grant_cost: T,
        pub shareholder_dilution_discount_factor: T,
    }

    pub fn evaluate_haug_eso_vesting_barrier_dilution<T: Real>(
        spec: &HaugEsoVestingBarrierDilutionSpec<T>,
    ) -> HaugEsoVestingBarrierDilutionReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
        let k = spec.strike_price_k.max_val(T::from_f64(1e-4));
        let r = spec.risk_free_rate_r;
        let q = spec.dividend_yield_q;
        let v = spec.volatility_sigma.max_val(T::from_f64(1e-4));
        let t = spec.contractual_maturity_t.max_val(T::from_f64(0.5));
        let tv = spec.vesting_period_tv.max_val(T::from_f64(0.1)).min_val(t);
        let lam_f = spec.forfeiture_hazard_lambda_f.max_val(T::zero());
        let m_mult = spec
            .suboptimal_exercise_multiple_m
            .max_val(T::from_f64(1.05));

        // Pre-vesting survival probability exp(-lambda_f * T_v):
        let surv_prob = (-lam_f * tv).exp();

        // Hull-White (2004) / FAS 123R effective life reduction from voluntary exercise barrier S^* = M * K:
        let log_m = m_mult.ln();
        let barrier_speed = (v * v
            + T::from_f64(2.0) * log_m / (t - tv).max_val(T::from_f64(0.25)))
        .max_val(T::from_f64(0.05));
        let t_eff = tv + (t - tv) / (T::one() + T::from_f64(0.5) * log_m * barrier_speed);

        // Plain Black-Scholes at contractual maturity T vs. ESO at effective life T_eff:
        let std_t = v * t.sqrt();
        let d1_t = ((s0 / k).ln() + (r - q + T::from_f64(0.5) * v * v) * t) / std_t;
        let d2_t = d1_t - std_t;
        let bs_full = s0 * (-q * t).exp() * norm_cdf(d1_t) - k * (-r * t).exp() * norm_cdf(d2_t);

        let std_e = v * t_eff.sqrt();
        let d1_e = ((s0 / k).ln() + (r - q + T::from_f64(0.5) * v * v) * t_eff) / std_e;
        let d2_e = d1_e - std_e;
        let bs_eff =
            s0 * (-q * t_eff).exp() * norm_cdf(d1_e) - k * (-r * t_eff).exp() * norm_cdf(d2_e);

        let eso_undiluted = surv_prob * bs_eff;
        let dil_factor =
            T::one() / (T::one() + spec.grant_to_shares_ratio_gamma.max_val(T::zero()));
        let eso_diluted = eso_undiluted * dil_factor;

        HaugEsoVestingBarrierDilutionReport {
            vesting_survival_probability: surv_prob,
            effective_expected_life_teff: t_eff,
            plain_black_scholes_call: bs_full,
            forfeiture_adjusted_eso_value: eso_undiluted,
            diluted_per_option_grant_cost: eso_diluted,
            shareholder_dilution_discount_factor: dil_factor,
        }
    }

    #[cfg(test)]
    mod domain_suite_53_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn haug_eso_vesting_barrier_dilution_invariants_and_ad_parity() {
            let spec = HaugEsoVestingBarrierDilutionSpec {
                spot_price_s0: 100.0,
                strike_price_k: 100.0,
                risk_free_rate_r: 0.04,
                dividend_yield_q: 0.015,
                volatility_sigma: 0.30,
                contractual_maturity_t: 10.0,
                vesting_period_tv: 3.0,
                forfeiture_hazard_lambda_f: 0.06,
                suboptimal_exercise_multiple_m: 1.8,
                grant_to_shares_ratio_gamma: 0.05,
            };
            let rep = evaluate_haug_eso_vesting_barrier_dilution(&spec);
            assert!(rep.effective_expected_life_teff < spec.contractual_maturity_t);
            assert!(rep.forfeiture_adjusted_eso_value < rep.plain_black_scholes_call);
            assert!(rep.diluted_per_option_grant_cost < rep.forfeiture_adjusted_eso_value);

            let h = 1e-6;
            let r_ad =
                evaluate_haug_eso_vesting_barrier_dilution(&HaugEsoVestingBarrierDilutionSpec {
                    spot_price_s0: Dual::<1>::variable(spec.spot_price_s0, 0),
                    strike_price_k: Dual::<1>::from_value(spec.strike_price_k),
                    risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
                    dividend_yield_q: Dual::<1>::from_value(spec.dividend_yield_q),
                    volatility_sigma: Dual::<1>::from_value(spec.volatility_sigma),
                    contractual_maturity_t: Dual::<1>::from_value(spec.contractual_maturity_t),
                    vesting_period_tv: Dual::<1>::from_value(spec.vesting_period_tv),
                    forfeiture_hazard_lambda_f: Dual::<1>::from_value(
                        spec.forfeiture_hazard_lambda_f,
                    ),
                    suboptimal_exercise_multiple_m: Dual::<1>::from_value(
                        spec.suboptimal_exercise_multiple_m,
                    ),
                    grant_to_shares_ratio_gamma: Dual::<1>::from_value(
                        spec.grant_to_shares_ratio_gamma,
                    ),
                });
            let mut sp = spec;
            sp.spot_price_s0 += h;
            let mut sm = spec;
            sm.spot_price_s0 -= h;
            let fd = (evaluate_haug_eso_vesting_barrier_dilution(&sp)
                .diluted_per_option_grant_cost
                - evaluate_haug_eso_vesting_barrier_dilution(&sm).diluted_per_option_grant_cost)
                / (2.0 * h);
            assert!((r_ad.diluted_per_option_grant_cost.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_haug_eso_vesting_barrier_dilution_module::*;

mod evaluate_longstaff_extendible_option_module {
    // Haug (2007, Ch. 4.9) & Longstaff (1990) Writer-Extendible & Holder-Extendible Option
    // Bivariate Normal Valuation Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[inline]
    fn smax<T: Real>(a: T, b: T, eps: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b + (d * d + eps * eps).sqrt())
    }

    #[derive(Debug, Clone, Copy)]
    pub struct LongstaffExtendibleOptionSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub initial_strike_k1: T,
        pub extended_strike_k2: T,
        pub initial_maturity_t1: T,
        pub extended_maturity_t2: T,
        pub holder_extension_fee_a: T,
        pub risk_free_rate_r: T,
        pub dividend_yield_q: T,
        pub volatility_sigma: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct LongstaffExtendibleOptionReport<T: Real = f64> {
        pub base_vanilla_call_t1: T,
        pub writer_extension_probability: T,
        pub writer_extension_value_addon: T,
        pub writer_extendible_call_value: T,
        pub holder_extendible_call_value: T,
        pub holder_over_writer_optionality_spread: T,
    }

    fn bvn_cdf_gl8<T: Real>(a: T, b: T, rho: T) -> T {
        let nodes = [
            -0.960_289_856_497_536_3,
            -0.796_666_477_413_626_7,
            -0.525_532_409_916_329,
            -0.183_434_642_495_649_8,
            0.183_434_642_495_649_8,
            0.525_532_409_916_329,
            0.796_666_477_413_626_7,
            0.960_289_856_497_536_3,
        ];
        let weights = [
            0.101_228_536_290_376_3,
            0.222_381_034_453_374_5,
            0.313_706_645_877_887_3,
            0.362_683_783_378_362,
            0.362_683_783_378_362,
            0.313_706_645_877_887_3,
            0.222_381_034_453_374_5,
            0.101_228_536_290_376_3,
        ];
        let r = rho.max_val(T::from_f64(-0.995)).min_val(T::from_f64(0.995));
        let half_r = T::from_f64(0.5) * r;
        let inv_two_pi = T::from_f64(0.159_154_943_091_895);
        let mut integral = T::zero();
        for (&node, &weight) in nodes.iter().zip(weights.iter()) {
            let xi = half_r * (T::from_f64(node) + T::one());
            let omx2 = (T::one() - xi * xi).max_val(T::from_f64(1e-8));
            let num = -(a * a - T::from_f64(2.0) * xi * a * b + b * b) / (T::from_f64(2.0) * omx2);
            integral += T::from_f64(weight) * num.exp() / omx2.sqrt();
        }
        norm_cdf(a) * norm_cdf(b) + half_r * inv_two_pi * integral
    }

    pub fn evaluate_longstaff_extendible_option<T: Real>(
        spec: &LongstaffExtendibleOptionSpec<T>,
    ) -> LongstaffExtendibleOptionReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
        let k1 = spec.initial_strike_k1.max_val(T::from_f64(1e-4));
        let k2 = spec.extended_strike_k2.max_val(T::from_f64(1e-4));
        let t1 = spec.initial_maturity_t1.max_val(T::from_f64(0.05));
        let t2 = spec.extended_maturity_t2.max_val(t1 + T::from_f64(0.10));
        let fee_a = spec.holder_extension_fee_a.max_val(T::zero());
        let r = spec.risk_free_rate_r;
        let q = spec.dividend_yield_q;
        let sig = spec.volatility_sigma.max_val(T::from_f64(1e-4));

        let sqrt_t1 = t1.sqrt();
        let sqrt_t2 = t2.sqrt();
        let rho = (t1 / t2).sqrt();

        let z1 = ((s0 / k2).ln() + (r - q + T::from_f64(0.5) * sig * sig) * t2) / (sig * sqrt_t2);
        let z2 = ((s0 / k1).ln() + (r - q + T::from_f64(0.5) * sig * sig) * t1) / (sig * sqrt_t1);

        // Base vanilla call C(S0, K1, t1):
        let base_call = s0 * (-q * t1).exp() * norm_cdf(z2)
            - k1 * (-r * t1).exp() * norm_cdf(z2 - sig * sqrt_t1);

        // Longstaff (1990) / Haug (2007, Sec. 4.9.1) Writer-Extendible Call:
        // Extends to (K2, t2) automatically if S(t1) < K1 (out of the money at t1)
        let p_ext_writer = norm_cdf(-(z2 - sig * sqrt_t1));
        let writer_addon = s0 * (-q * t2).exp() * bvn_cdf_gl8(z1, -z2, -rho)
            - k2 * (-r * t2).exp() * bvn_cdf_gl8(z1 - sig * sqrt_t2, -(z2 - sig * sqrt_t1), -rho);
        let writer_ext_call = base_call + writer_addon;

        // Longstaff (1990) / Haug (2007, Sec. 4.9.2) Holder-Extendible Call:
        // At t1, holder chooses max(S(t1) - K1, C_BS(S(t1), K2, t2 - t1) - A, 0)
        // Evaluated via 8-point Gauss-Legendre quadrature over standard normal z in [-4.5, 4.5]:
        let gl_nodes = [
            -0.960_289_856_497_536_3,
            -0.796_666_477_413_626_7,
            -0.525_532_409_916_329,
            -0.183_434_642_495_649_8,
            0.183_434_642_495_649_8,
            0.525_532_409_916_329,
            0.796_666_477_413_626_7,
            0.960_289_856_497_536_3,
        ];
        let gl_weights = [
            0.101_228_536_290_376_3,
            0.222_381_034_453_374_5,
            0.313_706_645_877_887_3,
            0.362_683_783_378_362,
            0.362_683_783_378_362,
            0.313_706_645_877_887_3,
            0.222_381_034_453_374_5,
            0.101_228_536_290_376_3,
        ];
        let dt = t2 - t1;
        let sqrt_dt = dt.sqrt();
        let half_span = T::from_f64(4.5);
        let eps = T::from_f64(1e-4);
        let drift_t1 = (r - q - T::from_f64(0.5) * sig * sig) * t1;
        let mut holder_expectation = T::zero();
        for (&node, &weight) in gl_nodes.iter().zip(gl_weights.iter()) {
            let z = half_span * T::from_f64(node);
            let s_t1 = s0 * (drift_t1 + sig * sqrt_t1 * z).exp();
            let d1_rem =
                ((s_t1 / k2).ln() + (r - q + T::from_f64(0.5) * sig * sig) * dt) / (sig * sqrt_dt);
            let d2_rem = d1_rem - sig * sqrt_dt;
            let c_rem =
                s_t1 * (-q * dt).exp() * norm_cdf(d1_rem) - k2 * (-r * dt).exp() * norm_cdf(d2_rem);
            let imm_ex = s_t1 - k1;
            let ext_val = c_rem - fee_a;
            let best_active = smax(imm_ex, ext_val, eps);
            let payoff_t1 = smax(best_active, T::zero(), eps);
            holder_expectation += T::from_f64(weight) * half_span * norm_pdf(z) * payoff_t1;
        }
        let holder_ext_call = (-r * t1).exp() * holder_expectation;
        let holder_over_writer = holder_ext_call - writer_ext_call;

        LongstaffExtendibleOptionReport {
            base_vanilla_call_t1: base_call,
            writer_extension_probability: p_ext_writer,
            writer_extension_value_addon: writer_addon,
            writer_extendible_call_value: writer_ext_call,
            holder_extendible_call_value: holder_ext_call,
            holder_over_writer_optionality_spread: holder_over_writer,
        }
    }

    #[cfg(test)]
    mod evaluate_longstaff_extendible_option_module_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn longstaff_extendible_option_invariants_and_ad_parity() {
            let spec = LongstaffExtendibleOptionSpec {
                spot_price_s0: 100.0,
                initial_strike_k1: 102.0,
                extended_strike_k2: 105.0,
                initial_maturity_t1: 0.50,
                extended_maturity_t2: 1.25,
                holder_extension_fee_a: 1.50,
                risk_free_rate_r: 0.045,
                dividend_yield_q: 0.015,
                volatility_sigma: 0.25,
            };
            let rep = evaluate_longstaff_extendible_option(&spec);
            assert!(rep.writer_extendible_call_value > rep.base_vanilla_call_t1);
            assert!(rep.holder_extendible_call_value > rep.base_vanilla_call_t1);
            let mut spec_zero_fee = spec;
            spec_zero_fee.holder_extension_fee_a = 0.0;
            let rep_zero_fee = evaluate_longstaff_extendible_option(&spec_zero_fee);
            assert!(
                rep_zero_fee.holder_extendible_call_value
                    > rep_zero_fee.writer_extendible_call_value
            );

            let h = 1e-6;
            let r_ad = evaluate_longstaff_extendible_option(&LongstaffExtendibleOptionSpec {
                spot_price_s0: Dual::<1>::variable(spec.spot_price_s0, 0),
                initial_strike_k1: Dual::<1>::from_value(spec.initial_strike_k1),
                extended_strike_k2: Dual::<1>::from_value(spec.extended_strike_k2),
                initial_maturity_t1: Dual::<1>::from_value(spec.initial_maturity_t1),
                extended_maturity_t2: Dual::<1>::from_value(spec.extended_maturity_t2),
                holder_extension_fee_a: Dual::<1>::from_value(spec.holder_extension_fee_a),
                risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
                dividend_yield_q: Dual::<1>::from_value(spec.dividend_yield_q),
                volatility_sigma: Dual::<1>::from_value(spec.volatility_sigma),
            });
            let mut sp = spec;
            sp.spot_price_s0 += h;
            let mut sm = spec;
            sm.spot_price_s0 -= h;
            let fd = (evaluate_longstaff_extendible_option(&sp).writer_extendible_call_value
                - evaluate_longstaff_extendible_option(&sm).writer_extendible_call_value)
                / (2.0 * h);
            assert!((r_ad.writer_extendible_call_value.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_longstaff_extendible_option_module::*;

mod evaluate_hyer_lipton_pugachevsky_passport_option_module {
    // Hyer, Lipton & Pugachevsky (1997, Math. Finance) & Haug (2007, Ch. 4.16)
    // Passport Option (Option on a Trading Account) & Optimal Bang-Bang Leverage Strategy Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct HyerLiptonPugachevskyPassportOptionSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub initial_trading_account_w0: T,
        pub max_position_limit_q_max: T,
        pub risk_free_rate_r: T,
        pub cost_of_carry_b: T,
        pub volatility_sigma: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct HyerLiptonPugachevskyPassportOptionReport<T: Real = f64> {
        pub zero_carry_passport_value: T,
        pub carry_adjusted_passport_value: T,
        pub benchmark_atm_vanilla_call: T,
        pub passport_to_vanilla_ratio: T,
        pub optimal_switching_local_time_value: T,
    }

    pub fn evaluate_hyer_lipton_pugachevsky_passport_option<T: Real>(
        spec: &HyerLiptonPugachevskyPassportOptionSpec<T>,
    ) -> HyerLiptonPugachevskyPassportOptionReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1.0));
        let w0 = spec.initial_trading_account_w0;
        let q_max = spec.max_position_limit_q_max.max_val(T::from_f64(0.1));
        let r = spec.risk_free_rate_r;
        let b = spec.cost_of_carry_b;
        let sig = spec.volatility_sigma.max_val(T::from_f64(0.05));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let std = sig * t.sqrt();
        let d_half = T::from_f64(0.5) * std;

        // Zero-carry (b = 0) exact Hyer-Lipton-Pugachevsky / Haug (2007, Ch. 4.16) base value
        // plus smooth initial wealth offset w0:
        let scale_s = q_max * s0 * (-r * t).exp();
        let base_zero_carry = scale_s * (norm_cdf(d_half) - norm_cdf(-d_half));

        let w_norm = w0 / (q_max * s0).max_val(T::from_f64(1.0));
        let c_w = T::from_f64(0.20) * scale_s;
        let wealth_adj = T::from_f64(0.5) * w0 + T::from_f64(0.5) * (w0 * w0 + c_w * c_w).sqrt()
            - T::from_f64(0.5) * c_w;
        let v_zero_carry = base_zero_carry + wealth_adj;

        // Carry adjustment when b != 0 (holder tilts toward positive carry sign(b)):
        let carry_diff = (b * t).exp() - T::one();
        let smooth_abs_carry =
            (carry_diff * carry_diff + T::from_f64(1e-12)).sqrt() - T::from_f64(1e-6);
        let carry_gain = scale_s * smooth_abs_carry * norm_cdf(d_half)
            + T::from_f64(0.25) * scale_s * b * t * (T::one() + w_norm);
        let v_passport = (v_zero_carry + carry_gain).max_val(T::from_f64(1e-4));

        // Benchmark ATM Black-Scholes call on q_max units:
        let d1_van = (b * t + T::from_f64(0.5) * std * std) / std;
        let d2_van = d1_van - std;
        let van_call = q_max
            * s0
            * (((b - r) * t).exp() * norm_cdf(d1_van) - (-r * t).exp() * norm_cdf(d2_van));

        let ratio = v_passport / van_call.max_val(T::from_f64(1e-4));
        let switch_val = (v_passport - T::from_f64(0.5) * van_call).max_val(T::zero());

        HyerLiptonPugachevskyPassportOptionReport {
            zero_carry_passport_value: v_zero_carry,
            carry_adjusted_passport_value: v_passport,
            benchmark_atm_vanilla_call: van_call,
            passport_to_vanilla_ratio: ratio,
            optimal_switching_local_time_value: switch_val,
        }
    }

    #[cfg(test)]
    mod evaluate_hyer_lipton_pugachevsky_passport_option_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hyer_lipton_pugachevsky_passport_option_invariants_and_ad_parity() {
            let spec = HyerLiptonPugachevskyPassportOptionSpec {
                spot_price_s0: 100.0,
                initial_trading_account_w0: 2.0,
                max_position_limit_q_max: 1.0,
                risk_free_rate_r: 0.04,
                cost_of_carry_b: 0.02,
                volatility_sigma: 0.25,
                maturity_years_t: 1.0,
            };
            let rep = evaluate_hyer_lipton_pugachevsky_passport_option(&spec);
            assert!(rep.carry_adjusted_passport_value > 0.0);
            assert!(rep.passport_to_vanilla_ratio > 0.5);

            let mut spec_zero = spec;
            spec_zero.initial_trading_account_w0 = 0.0;
            spec_zero.cost_of_carry_b = 0.0;
            let rep_zero = evaluate_hyer_lipton_pugachevsky_passport_option(&spec_zero);
            assert!(
                (rep_zero.zero_carry_passport_value - rep_zero.benchmark_atm_vanilla_call).abs()
                    < 1e-12
            );

            let h = 1e-5;
            let r_ad = evaluate_hyer_lipton_pugachevsky_passport_option(
                &HyerLiptonPugachevskyPassportOptionSpec {
                    spot_price_s0: Dual::<7>::variable(spec.spot_price_s0, 0),
                    initial_trading_account_w0: Dual::<7>::variable(
                        spec.initial_trading_account_w0,
                        1,
                    ),
                    max_position_limit_q_max: Dual::<7>::variable(spec.max_position_limit_q_max, 2),
                    risk_free_rate_r: Dual::<7>::variable(spec.risk_free_rate_r, 3),
                    cost_of_carry_b: Dual::<7>::variable(spec.cost_of_carry_b, 4),
                    volatility_sigma: Dual::<7>::variable(spec.volatility_sigma, 5),
                    maturity_years_t: Dual::<7>::variable(spec.maturity_years_t, 6),
                },
            );
            for (idx, &ad) in r_ad.carry_adjusted_passport_value.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_price_s0 += h;
                        sm.spot_price_s0 -= h;
                    }
                    1 => {
                        sp.initial_trading_account_w0 += h;
                        sm.initial_trading_account_w0 -= h;
                    }
                    2 => {
                        sp.max_position_limit_q_max += h;
                        sm.max_position_limit_q_max -= h;
                    }
                    3 => {
                        sp.risk_free_rate_r += h;
                        sm.risk_free_rate_r -= h;
                    }
                    4 => {
                        sp.cost_of_carry_b += h;
                        sm.cost_of_carry_b -= h;
                    }
                    5 => {
                        sp.volatility_sigma += h;
                        sm.volatility_sigma -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_hyer_lipton_pugachevsky_passport_option(&sp)
                    .carry_adjusted_passport_value
                    - evaluate_hyer_lipton_pugachevsky_passport_option(&sm)
                        .carry_adjusted_passport_value)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_hyer_lipton_pugachevsky_passport_option_module::*;

mod evaluate_rubinstein_margrabe_exchange_on_exchange_option_module {
    // Rubinstein (1991, Risk, "Somewhere Over the Rainbow") & Haug (2007, The Complete Guide
    // to Option Pricing Formulas, Ch. 5.6, p. 216) Two-Asset Exchange-Option-on-an-Exchange-Option
    // & Margrabe Outperformance Compound Convexity Engine.

    use ql_ad::Real;
    use ql_math::{bivariate_normal_cdf, norm_cdf};

    #[derive(Debug, Clone, Copy)]
    pub struct RubinsteinMargrabeExchangeOnExchangeSpec<T: Real = f64> {
        pub spot_asset1_s1: T,
        pub spot_asset2_s2: T,
        pub compound_strike_fraction_q: T,
        pub dividend_yield_q1: T,
        pub dividend_yield_q2: T,
        pub volatility_sigma1: T,
        pub volatility_sigma2: T,
        pub asset_correlation_rho: T,
        pub underlying_expiry_years_t2: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct RubinsteinMargrabeExchangeOnExchangeReport<T: Real = f64> {
        pub margrabe_relative_volatility: T,
        pub underlying_margrabe_exchange_value: T,
        pub critical_exercise_price_ratio_u_star: T,
        pub temporal_brownian_correlation_rho12: T,
        pub exchange_on_exchange_call_value: T,
        pub compound_leverage_ratio: T,
    }

    pub fn evaluate_rubinstein_margrabe_exchange_on_exchange_option<T: Real>(
        spec: &RubinsteinMargrabeExchangeOnExchangeSpec<T>,
    ) -> RubinsteinMargrabeExchangeOnExchangeReport<T> {
        let s1 = spec.spot_asset1_s1.max_val(T::from_f64(1.0));
        let s2 = spec.spot_asset2_s2.max_val(T::from_f64(1.0));
        let q_strike = spec
            .compound_strike_fraction_q
            .max_val(T::from_f64(0.005))
            .min_val(T::from_f64(0.50));
        let q1 = spec.dividend_yield_q1;
        let q2 = spec.dividend_yield_q2;
        let sig1 = spec.volatility_sigma1.max_val(T::from_f64(0.02));
        let sig2 = spec.volatility_sigma2.max_val(T::from_f64(0.02));
        let rho = spec
            .asset_correlation_rho
            .max_val(T::from_f64(-0.95))
            .min_val(T::from_f64(0.95));
        let t2 = spec.underlying_expiry_years_t2.max_val(T::from_f64(0.20));
        let t1 = T::from_f64(0.40) * t2;
        let tau = t2 - t1;

        let sig_m = (sig1 * sig1 + sig2 * sig2 - T::from_f64(2.0) * rho * sig1 * sig2)
            .max_val(T::from_f64(1e-4))
            .sqrt();

        let std2 = sig_m * t2.sqrt();
        let d1_m = ((s1 / s2).ln() + (q2 - q1 + T::from_f64(0.5) * sig_m * sig_m) * t2) / std2;
        let d2_m = d1_m - std2;
        let v_margrabe =
            s1 * (-q1 * t2).exp() * norm_cdf(d1_m) - s2 * (-q2 * t2).exp() * norm_cdf(d2_m);

        // Solve critical ratio u* = S_1(T_1) / S_2(T_1) where normalized Margrabe value m(u*) == q_strike
        // via 6-step Newton-Raphson (preserving exact implicit-function Dual/Dual2 derivatives):
        let std_tau = sig_m * tau.sqrt();
        let df1_tau = (-q1 * tau).exp();
        let df2_tau = (-q2 * tau).exp();
        let mut u = T::one() + q_strike;
        for _ in 0..6 {
            let d1_u = (u.ln() + (q2 - q1 + T::from_f64(0.5) * sig_m * sig_m) * tau) / std_tau;
            let d2_u = d1_u - std_tau;
            let m_u = u * df1_tau * norm_cdf(d1_u) - df2_tau * norm_cdf(d2_u);
            let mp_u = (df1_tau * norm_cdf(d1_u)).max_val(T::from_f64(1e-4));
            u = (u - (m_u - q_strike) / mp_u).max_val(T::from_f64(0.20));
        }

        let std1 = sig_m * t1.sqrt();
        let y1 = ((s1 / (u * s2)).ln() + (q2 - q1 + T::from_f64(0.5) * sig_m * sig_m) * t1) / std1;
        let y2 = y1 - std1;
        let rho12 = (t1 / t2).sqrt();

        let v_exonex = s1 * (-q1 * t2).exp() * bivariate_normal_cdf(y1, d1_m, rho12)
            - s2 * (-q2 * t2).exp() * bivariate_normal_cdf(y2, d2_m, rho12)
            - q_strike * s2 * (-q2 * t1).exp() * norm_cdf(y2);
        let lev = v_margrabe / v_exonex.max_val(T::from_f64(1e-4));

        RubinsteinMargrabeExchangeOnExchangeReport {
            margrabe_relative_volatility: sig_m,
            underlying_margrabe_exchange_value: v_margrabe,
            critical_exercise_price_ratio_u_star: u,
            temporal_brownian_correlation_rho12: rho12,
            exchange_on_exchange_call_value: v_exonex,
            compound_leverage_ratio: lev,
        }
    }

    #[cfg(test)]
    mod evaluate_rubinstein_margrabe_exchange_on_exchange_option_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn rubinstein_margrabe_exchange_on_exchange_invariants_and_ad_parity() {
            let spec = RubinsteinMargrabeExchangeOnExchangeSpec {
                spot_asset1_s1: 105.0,
                spot_asset2_s2: 100.0,
                compound_strike_fraction_q: 0.04,
                dividend_yield_q1: 0.015,
                dividend_yield_q2: 0.025,
                volatility_sigma1: 0.24,
                volatility_sigma2: 0.20,
                asset_correlation_rho: 0.55,
                underlying_expiry_years_t2: 1.50,
            };
            let rep = evaluate_rubinstein_margrabe_exchange_on_exchange_option(&spec);
            assert!(rep.exchange_on_exchange_call_value > 0.0);
            assert!(rep.exchange_on_exchange_call_value < rep.underlying_margrabe_exchange_value);
            assert!(rep.critical_exercise_price_ratio_u_star > 0.8);

            let h = 1e-5;
            let r_ad = evaluate_rubinstein_margrabe_exchange_on_exchange_option(
                &RubinsteinMargrabeExchangeOnExchangeSpec {
                    spot_asset1_s1: Dual::<9>::variable(spec.spot_asset1_s1, 0),
                    spot_asset2_s2: Dual::<9>::variable(spec.spot_asset2_s2, 1),
                    compound_strike_fraction_q: Dual::<9>::variable(
                        spec.compound_strike_fraction_q,
                        2,
                    ),
                    dividend_yield_q1: Dual::<9>::variable(spec.dividend_yield_q1, 3),
                    dividend_yield_q2: Dual::<9>::variable(spec.dividend_yield_q2, 4),
                    volatility_sigma1: Dual::<9>::variable(spec.volatility_sigma1, 5),
                    volatility_sigma2: Dual::<9>::variable(spec.volatility_sigma2, 6),
                    asset_correlation_rho: Dual::<9>::variable(spec.asset_correlation_rho, 7),
                    underlying_expiry_years_t2: Dual::<9>::variable(
                        spec.underlying_expiry_years_t2,
                        8,
                    ),
                },
            );
            for (idx, &ad) in r_ad.exchange_on_exchange_call_value.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_asset1_s1 += h;
                        sm.spot_asset1_s1 -= h;
                    }
                    1 => {
                        sp.spot_asset2_s2 += h;
                        sm.spot_asset2_s2 -= h;
                    }
                    2 => {
                        sp.compound_strike_fraction_q += h;
                        sm.compound_strike_fraction_q -= h;
                    }
                    3 => {
                        sp.dividend_yield_q1 += h;
                        sm.dividend_yield_q1 -= h;
                    }
                    4 => {
                        sp.dividend_yield_q2 += h;
                        sm.dividend_yield_q2 -= h;
                    }
                    5 => {
                        sp.volatility_sigma1 += h;
                        sm.volatility_sigma1 -= h;
                    }
                    6 => {
                        sp.volatility_sigma2 += h;
                        sm.volatility_sigma2 -= h;
                    }
                    7 => {
                        sp.asset_correlation_rho += h;
                        sm.asset_correlation_rho -= h;
                    }
                    _ => {
                        sp.underlying_expiry_years_t2 += h;
                        sm.underlying_expiry_years_t2 -= h;
                    }
                }
                let fd = (evaluate_rubinstein_margrabe_exchange_on_exchange_option(&sp)
                    .exchange_on_exchange_call_value
                    - evaluate_rubinstein_margrabe_exchange_on_exchange_option(&sm)
                        .exchange_on_exchange_call_value)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_rubinstein_margrabe_exchange_on_exchange_option_module::*;

mod evaluate_turnbull_wakeman_edgeworth_arithmetic_asian_option_module {
    // Turnbull & Wakeman (1991, Journal of Financial and Quantitative Analysis 26(3),
    // "A Quick Algorithm for Pricing European Average Options"), Levy (1992, JIMF 11(5)),
    // & Kemna & Vorst (1990, JBF 14(1))
    // Moment-Matched Arithmetic Asian Call vs. Geometric Asian Dispersion Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct TurnbullWakemanArithmeticAsianSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub strike_price_k: T,
        pub risk_free_rate_r: T,
        pub dividend_yield_q: T,
        pub volatility_sigma: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct TurnbullWakemanArithmeticAsianReport<T: Real = f64> {
        pub arithmetic_average_forward_m1: T,
        pub turnbull_wakeman_effective_volatility: T,
        pub turnbull_wakeman_arithmetic_call_value: T,
        pub kemna_vorst_geometric_call_value: T,
        pub arithmetic_over_geometric_dispersion_spread: T,
        pub asian_call_spot_delta: T,
    }

    pub fn evaluate_turnbull_wakeman_edgeworth_arithmetic_asian_option<T: Real>(
        spec: &TurnbullWakemanArithmeticAsianSpec<T>,
    ) -> TurnbullWakemanArithmeticAsianReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1.0));
        let k = spec.strike_price_k.max_val(T::from_f64(1.0));
        let r = spec.risk_free_rate_r;
        let q = spec.dividend_yield_q;
        let sig = spec.volatility_sigma.max_val(T::from_f64(0.02));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let b = (r - q).max_val(T::from_f64(0.005));
        let var = sig * sig;
        let bt = b * t;
        let e_bt = bt.exp();

        // Exact continuous first and second moments of (1/T) * int_0^T (S_u / S_0) du:
        let m1_ratio = (e_bt - T::one()) / bt;
        let bpv = b + var;
        let twobpv = T::from_f64(2.0) * b + var;
        let t2 = t * t;
        let m2_ratio = (T::from_f64(2.0) * (twobpv * t).exp()) / (bpv * twobpv * t2)
            + (T::from_f64(2.0) / (b * t2)) * (T::one() / twobpv - e_bt / bpv);

        let f_arith = s0 * m1_ratio;
        let var_tw = ((m2_ratio / (m1_ratio * m1_ratio)).max_val(T::from_f64(1.0001))).ln() / t;
        let sig_tw = var_tw.sqrt();

        let sqrt_t = t.sqrt();
        let df_r = (-r * t).exp();

        // Exact continuous third moment M3 of (1/T) * int_0^T (S_u / S_0) du:
        let a1 = b + T::from_f64(2.0) * var;
        let a2 = bpv;
        let a12 = a1 + a2;
        let a123 = a12 + b;
        let a23 = twobpv;
        let int_123 = (((a123 * t).exp() - T::one()) / a123 - (e_bt - T::one()) / b) / a12;
        let int_23 = (((a23 * t).exp() - T::one()) / a23 - (e_bt - T::one()) / b) / a2;
        let m3_ratio = (T::from_f64(6.0) / (t2 * t * a1)) * (int_123 - int_23);
        let m3_ln_ratio = (m2_ratio * m2_ratio * m2_ratio) / (m1_ratio * m1_ratio * m1_ratio);
        let dk3 = s0 * s0 * s0 * (m3_ratio - m3_ln_ratio);

        let std_tw = sig_tw * sqrt_t;
        let d1_tw = ((f_arith / k).ln() + T::from_f64(0.5) * var_tw * t) / std_tw;
        let d2_tw = d1_tw - std_tw;
        let c_ln = df_r * (f_arith * norm_cdf(d1_tw) - k * norm_cdf(d2_tw));
        let edgeworth_adj = df_r
            * (dk3 / T::from_f64(6.0))
            * (norm_pdf(d2_tw) / (k * k * std_tw))
            * (T::one() - d2_tw / std_tw);
        let c_tw = c_ln + edgeworth_adj;

        // Kemna & Vorst (1990) exact continuous geometric Asian call benchmark:
        let sig_g = sig / T::from_f64(3.0).sqrt();
        let b_g = T::from_f64(0.5) * (b - var / T::from_f64(6.0));
        let f_geom = s0 * (b_g * t).exp();
        let std_g = sig_g * sqrt_t;
        let d1_g = ((f_geom / k).ln() + T::from_f64(0.5) * sig_g * sig_g * t) / std_g;
        let d2_g = d1_g - std_g;
        let c_kv = df_r * (f_geom * norm_cdf(d1_g) - k * norm_cdf(d2_g));

        let spread = c_tw - c_kv;
        let d2_over_std = d2_tw / std_tw;
        let d_edgeworth_ds0 = df_r
            * (s0 * s0 * (m3_ratio - m3_ln_ratio) / T::from_f64(6.0))
            * (norm_pdf(d2_tw) / (k * k * std_tw))
            * ((T::from_f64(3.0) - d2_over_std) * (T::one() - d2_over_std)
                - T::one() / (std_tw * std_tw));
        let delta = df_r * m1_ratio * norm_cdf(d1_tw) + d_edgeworth_ds0;

        TurnbullWakemanArithmeticAsianReport {
            arithmetic_average_forward_m1: f_arith,
            turnbull_wakeman_effective_volatility: sig_tw,
            turnbull_wakeman_arithmetic_call_value: c_tw,
            kemna_vorst_geometric_call_value: c_kv,
            arithmetic_over_geometric_dispersion_spread: spread,
            asian_call_spot_delta: delta,
        }
    }

    #[cfg(test)]
    mod evaluate_turnbull_wakeman_edgeworth_arithmetic_asian_option_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn turnbull_wakeman_arithmetic_asian_invariants_and_ad_parity() {
            let spec = TurnbullWakemanArithmeticAsianSpec {
                spot_price_s0: 100.0,
                strike_price_k: 100.0,
                risk_free_rate_r: 0.05,
                dividend_yield_q: 0.015,
                volatility_sigma: 0.25,
                maturity_years_t: 1.0,
            };
            let rep = evaluate_turnbull_wakeman_edgeworth_arithmetic_asian_option(&spec);
            // By Jensen's inequality (AM >= GM), arithmetic Asian call strictly exceeds geometric Asian call:
            assert!(
                rep.turnbull_wakeman_arithmetic_call_value > rep.kemna_vorst_geometric_call_value
            );
            assert!(rep.arithmetic_over_geometric_dispersion_spread > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_turnbull_wakeman_edgeworth_arithmetic_asian_option(
                &TurnbullWakemanArithmeticAsianSpec {
                    spot_price_s0: Dual::<6>::variable(spec.spot_price_s0, 0),
                    strike_price_k: Dual::<6>::variable(spec.strike_price_k, 1),
                    risk_free_rate_r: Dual::<6>::variable(spec.risk_free_rate_r, 2),
                    dividend_yield_q: Dual::<6>::variable(spec.dividend_yield_q, 3),
                    volatility_sigma: Dual::<6>::variable(spec.volatility_sigma, 4),
                    maturity_years_t: Dual::<6>::variable(spec.maturity_years_t, 5),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.spot_price_s0,
                    spec.strike_price_k,
                    spec.risk_free_rate_r,
                    spec.dividend_yield_q,
                    spec.volatility_sigma,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_turnbull_wakeman_edgeworth_arithmetic_asian_option(
                    &TurnbullWakemanArithmeticAsianSpec {
                        spot_price_s0: v[0],
                        strike_price_k: v[1],
                        risk_free_rate_r: v[2],
                        dividend_yield_q: v[3],
                        volatility_sigma: v[4],
                        maturity_years_t: v[5],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_turnbull_wakeman_edgeworth_arithmetic_asian_option(
                    &TurnbullWakemanArithmeticAsianSpec {
                        spot_price_s0: v[0],
                        strike_price_k: v[1],
                        risk_free_rate_r: v[2],
                        dividend_yield_q: v[3],
                        volatility_sigma: v[4],
                        maturity_years_t: v[5],
                    },
                );
                let fd = (up.turnbull_wakeman_arithmetic_call_value
                    - dn.turnbull_wakeman_arithmetic_call_value)
                    / (2.0 * h);
                assert!(
                    (r_ad.turnbull_wakeman_arithmetic_call_value.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_turnbull_wakeman_edgeworth_arithmetic_asian_option_module::*;

/// Specification for the Wilmott (2002) / Gatheral (2006) globally floored, locally capped/floored
/// multi-reset cliquet option with forward volatility skew sensitivity.
///
/// Reference: Wilmott (2002, *IJTAF* 5(6), "Cliquet Options and Volatility Models"),
/// Gatheral (2006, *The Volatility Surface*, Ch. 10), Carr & Lee (2009).
#[derive(Debug, Clone, Copy)]
pub struct WilmottCarrCliquetSpec<T: Real> {
    /// Contract notional $N$.
    pub notional_n: T,
    /// Local period return floor $F_{\text{loc}}$ (e.g. $0.0$).
    pub local_return_floor_f_loc: T,
    /// Local period return cap $C_{\text{loc}}$ (e.g. $0.08$).
    pub local_return_cap_c_loc: T,
    /// Global cumulative return floor $F_{\text{glob}}$ across all reset periods.
    pub global_cumulative_floor_f_glob: T,
    /// Risk-free interest rate $r$.
    pub risk_free_rate_r: T,
    /// Dividend yield $q$.
    pub dividend_yield_q: T,
    /// Forward-starting ATM volatility $\sigma_{\text{ATM}}$.
    pub forward_atm_volatility_sigma: T,
    /// Forward smile skew slope $\eta_{\text{skew}}$ per unit moneyness.
    pub forward_skew_slope_eta: T,
    /// Total maturity $T$ across $N_{\text{reset}} = 4$ reset periods.
    pub maturity_years_t: T,
}

/// Analytics returned by [`evaluate_wilmott_carr_globally_floored_locally_capped_cliquet`].
#[derive(Debug, Clone, Copy)]
pub struct WilmottCarrCliquetAnalytics<T: Real> {
    /// Expected locally capped/floored single-period return $\mathbb{E}^{\mathbb{Q}}[Y_i]$.
    pub expected_single_period_collared_return: T,
    /// Discounted value of the global floor protection $\mathbb{E}^{\mathbb{Q}}[(F_{\text{glob}} - \sum Y_i)^+]$.
    pub embedded_global_floor_option_value: T,
    /// Forward-skew convexity adjustment relative to flat Black–Scholes valuation.
    pub forward_skew_valuation_wedge: T,
    /// Present value of the globally floored, locally capped cliquet structure.
    pub total_cliquet_present_value: T,
}

/// Evaluates a 4-period globally floored, locally capped/floored cliquet via exact truncated lognormal
/// local collar moments and Bachelier/Gaussian global floor aggregation with forward-skew adjustment.
pub fn evaluate_wilmott_carr_globally_floored_locally_capped_cliquet<T: Real>(
    spec: &WilmottCarrCliquetSpec<T>,
) -> WilmottCarrCliquetAnalytics<T> {
    const N_PERIODS: usize = 4;
    let n_f = T::from_f64(N_PERIODS as f64);
    let t_mat = spec.maturity_years_t.max_val(T::from_f64(1e-4));
    let dt = t_mat / n_f;
    let sqrt_dt = dt.sqrt();
    let carry = spec.risk_free_rate_r - spec.dividend_yield_q;
    let fwd_gross = (carry * dt).exp();

    let k_floor = (T::one() + spec.local_return_floor_f_loc).max_val(T::from_f64(0.1));
    let k_cap = (T::one() + spec.local_return_cap_c_loc).max_val(k_floor + T::from_f64(1e-4));

    // Helper: undiscounted expected call E[(S_dt / S_0 - K)^+] and second moment contribution under vol_k
    let bs_undisc_call = |k_gross: T, vol: T| -> T {
        let v = vol.max_val(T::from_f64(1e-4));
        let v_sqrt = v * sqrt_dt;
        let d1 = ((fwd_gross / k_gross).ln() + T::from_f64(0.5) * v * v * dt) / v_sqrt;
        let d2 = d1 - v_sqrt;
        fwd_gross * norm_cdf(d1) - k_gross * norm_cdf(d2)
    };

    let sig_atm = spec.forward_atm_volatility_sigma.max_val(T::from_f64(1e-4));
    let sig_floor = sig_atm + spec.forward_skew_slope_eta * (k_floor - T::one());
    let sig_cap = sig_atm + spec.forward_skew_slope_eta * (k_cap - T::one());

    // Locally collared return Y_i = F_loc + (1 + R_i - k_floor)^+ - (1 + R_i - k_cap)^+
    let call_floor_skew = bs_undisc_call(k_floor, sig_floor);
    let call_cap_skew = bs_undisc_call(k_cap, sig_cap);
    let mu_y_skew = spec.local_return_floor_f_loc + call_floor_skew - call_cap_skew;

    let call_floor_flat = bs_undisc_call(k_floor, sig_atm);
    let call_cap_flat = bs_undisc_call(k_cap, sig_atm);
    let mu_y_flat = spec.local_return_floor_f_loc + call_floor_flat - call_cap_flat;

    // Exact second moment E[(Y_i - F_loc)^2] of the locally collared return via
    // Z^2 = ((X - k_floor)^+)^2 - ((X - k_cap)^+)^2 - 2*(k_cap - k_floor)*(X - k_cap)^+
    let bs_undisc_sq_call = |k_gross: T, vol: T| -> T {
        let v = vol.max_val(T::from_f64(1e-4));
        let v_sqrt = v * sqrt_dt;
        let d1 = ((fwd_gross / k_gross).ln() + T::from_f64(0.5) * v * v * dt) / v_sqrt;
        let d2 = d1 - v_sqrt;
        let d3 = d1 + v_sqrt;
        fwd_gross * fwd_gross * (v * v * dt).exp() * norm_cdf(d3)
            - T::from_f64(2.0) * k_gross * fwd_gross * norm_cdf(d1)
            + k_gross * k_gross * norm_cdf(d2)
    };
    let width = k_cap - k_floor;
    let ez = call_floor_flat - call_cap_flat;
    let ez2 = bs_undisc_sq_call(k_floor, sig_atm)
        - bs_undisc_sq_call(k_cap, sig_atm)
        - T::from_f64(2.0) * width * call_cap_flat;
    let var_y = (ez2 - ez * ez).max_val(T::from_f64(1e-8));
    let std_y = var_y.sqrt();

    // Sum S_N = sum_{i=1}^N Y_i has mean M_sum = N * mu_y_skew and std Sigma_sum = sqrt(N) * std_y
    let m_sum = n_f * mu_y_skew;
    let sig_sum = n_f.sqrt() * std_y;
    let z_glob = (spec.global_cumulative_floor_f_glob - m_sum) / sig_sum;
    let inv_sqrt_2pi = T::from_f64(0.3989422804014327);
    let phi_z = inv_sqrt_2pi * (-T::from_f64(0.5) * z_glob * z_glob).exp();

    let df_t = (-spec.risk_free_rate_r * t_mat).exp();
    let global_floor_put = spec.notional_n * df_t * sig_sum * (z_glob * norm_cdf(z_glob) + phi_z);
    let skew_wedge = spec.notional_n * df_t * n_f * (mu_y_skew - mu_y_flat);
    let total_pv = spec.notional_n * df_t * m_sum + global_floor_put;

    WilmottCarrCliquetAnalytics {
        expected_single_period_collared_return: mu_y_skew,
        embedded_global_floor_option_value: global_floor_put,
        forward_skew_valuation_wedge: skew_wedge,
        total_cliquet_present_value: total_pv,
    }
}

#[cfg(test)]
mod wilmott_carr_cliquet_suite64_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn wilmott_carr_cliquet_invariants_and_ad_parity() {
        let base = WilmottCarrCliquetSpec {
            notional_n: 100.0_f64,
            local_return_floor_f_loc: 0.00,
            local_return_cap_c_loc: 0.08,
            global_cumulative_floor_f_glob: 0.10,
            risk_free_rate_r: 0.04,
            dividend_yield_q: 0.015,
            forward_atm_volatility_sigma: 0.22,
            forward_skew_slope_eta: -0.25,
            maturity_years_t: 2.0,
        };
        let r0 = evaluate_wilmott_carr_globally_floored_locally_capped_cliquet(&base);
        assert!(
            r0.expected_single_period_collared_return > 0.0
                && r0.expected_single_period_collared_return < 0.08
        );
        assert!(r0.embedded_global_floor_option_value > 0.0);
        assert!(r0.forward_skew_valuation_wedge > 0.0);
        assert!(r0.total_cliquet_present_value > r0.embedded_global_floor_option_value);

        let inputs = [
            base.notional_n,
            base.local_return_floor_f_loc,
            base.local_return_cap_c_loc,
            base.global_cumulative_floor_f_glob,
            base.risk_free_rate_r,
            base.dividend_yield_q,
            base.forward_atm_volatility_sigma,
            base.forward_skew_slope_eta,
            base.maturity_years_t,
        ];
        let mut dual_in = [Dual::<9>::from_f64(0.0); 9];
        for i in 0..9 {
            dual_in[i] = Dual::variable(inputs[i], i);
        }
        let spec_ad = WilmottCarrCliquetSpec {
            notional_n: dual_in[0],
            local_return_floor_f_loc: dual_in[1],
            local_return_cap_c_loc: dual_in[2],
            global_cumulative_floor_f_glob: dual_in[3],
            risk_free_rate_r: dual_in[4],
            dividend_yield_q: dual_in[5],
            forward_atm_volatility_sigma: dual_in[6],
            forward_skew_slope_eta: dual_in[7],
            maturity_years_t: dual_in[8],
        };
        let r_ad = evaluate_wilmott_carr_globally_floored_locally_capped_cliquet(&spec_ad);
        let eps = 1e-6;
        for i in 0..9 {
            let mut p = inputs;
            let mut m = inputs;
            p[i] += eps;
            m[i] -= eps;
            let mk = |v: [f64; 9]| WilmottCarrCliquetSpec {
                notional_n: v[0],
                local_return_floor_f_loc: v[1],
                local_return_cap_c_loc: v[2],
                global_cumulative_floor_f_glob: v[3],
                risk_free_rate_r: v[4],
                dividend_yield_q: v[5],
                forward_atm_volatility_sigma: v[6],
                forward_skew_slope_eta: v[7],
                maturity_years_t: v[8],
            };
            let rp = evaluate_wilmott_carr_globally_floored_locally_capped_cliquet(&mk(p));
            let rm = evaluate_wilmott_carr_globally_floored_locally_capped_cliquet(&mk(m));
            let fd =
                (rp.total_cliquet_present_value - rm.total_cliquet_present_value) / (2.0 * eps);
            assert!(
                r_ad.total_cliquet_present_value.d[i].abs() > 1e-6,
                "zero derivative at {i}"
            );
            assert!((r_ad.total_cliquet_present_value.d[i] - fd).abs() < 1e-4);
        }
    }
}

pub mod pechtl_haug_time_switch_option_mod {
    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    /// Specification for the Pechtl (1995) / Haug (2007, Ch. 4.16) discrete time-switch
    /// occupation-time accumulator option with smile skew correction.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct PechtlHaugTimeSwitchSpec<T: Real> {
        /// Underlying spot price $S_0 > 0$.
        pub spot_price_s0: T,
        /// Threshold strike price $K > 0$ for in-the-money occupation accrual.
        pub strike_price_k: T,
        /// Annualized accrual payout rate $A > 0$ per unit occupation time.
        pub accrual_rate_a: T,
        /// Already-accumulated in-the-money occupation time $t_{\text{acc}} \ge 0$ (years).
        pub accumulated_occupation_time_t_acc: T,
        /// Continuous risk-free rate $r$.
        pub risk_free_rate_r: T,
        /// Continuous dividend yield $q$.
        pub dividend_yield_q: T,
        /// Black-Scholes volatility $\sigma > 0$.
        pub volatility_sigma: T,
        /// Risk-reversal smile skew slope $\nu_{\text{skew}} = \partial \sigma(K)/\partial \ln K$.
        pub smile_skew_slope_nu: T,
        /// Remaining maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Pechtl-Haug discrete time-switch accumulator option engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct PechtlHaugTimeSwitchAnalytics<T: Real> {
        /// Expected future in-the-money occupation time $\sum_{m=1}^M \Delta t \, \mathbb{Q}(S(t_m) > K)$.
        pub expected_in_the_money_occupation_time: T,
        /// Present value of already-accumulated occupation time $A \cdot t_{\text{acc}} e^{-r T}$.
        pub already_accumulated_present_value: T,
        /// Smile skew digital probability adjustment wedge.
        pub smile_skew_occupation_adjustment: T,
        /// Total discounted present value of the discrete time-switch option.
        pub time_switch_option_present_value: T,
    }

    /// Evaluates the Pechtl (1995) / Haug (2007, Ch. 4.16) discrete time-switch accumulator option
    /// across $M = 8$ monitoring dates with digital smile skew correction.
    pub fn evaluate_pechtl_haug_discrete_time_switch_accumulator_option<T: Real>(
        spec: &PechtlHaugTimeSwitchSpec<T>,
    ) -> PechtlHaugTimeSwitchAnalytics<T> {
        const MON_DATES: usize = 8;
        let t_mat = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let dt = t_mat / T::from_f64(MON_DATES as f64);
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
        let k = spec.strike_price_k.max_val(T::from_f64(1e-4));
        let sig = spec.volatility_sigma.max_val(T::from_f64(1e-4));
        let log_sk = (s0 / k).ln();
        let drift = spec.risk_free_rate_r - spec.dividend_yield_q - T::from_f64(0.5) * sig * sig;

        let mut occ_time = T::zero();
        let mut skew_wedge_pv = T::zero();
        let mut future_accrual_pv = T::zero();
        let df_t = (-spec.risk_free_rate_r * t_mat).exp();

        for m in 1..=MON_DATES {
            let tm = T::from_f64(m as f64) * dt;
            let v_sqrt = sig * tm.sqrt();
            let d2 = (log_sk + drift * tm) / v_sqrt;

            // Digital call probability under smile skew: Phi(d2) - (d sigma / d ln K) * sqrt(t_m) * phi(d2)
            let p_bs = norm_cdf(d2);
            let skew_corr = -spec.smile_skew_slope_nu * tm.sqrt() * norm_pdf(d2);
            let p_adj = p_bs + skew_corr;

            occ_time += dt * p_adj;
            skew_wedge_pv += df_t * spec.accrual_rate_a * dt * skew_corr;
            future_accrual_pv += df_t * spec.accrual_rate_a * dt * p_adj;
        }

        let acc_pv = df_t * spec.accrual_rate_a * spec.accumulated_occupation_time_t_acc;

        PechtlHaugTimeSwitchAnalytics {
            expected_in_the_money_occupation_time: occ_time,
            already_accumulated_present_value: acc_pv,
            smile_skew_occupation_adjustment: skew_wedge_pv,
            time_switch_option_present_value: future_accrual_pv + acc_pv,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn pechtl_haug_time_switch_invariants_and_ad_parity() {
            let base = PechtlHaugTimeSwitchSpec {
                spot_price_s0: 104.0_f64,
                strike_price_k: 100.0,
                accrual_rate_a: 25.0,
                accumulated_occupation_time_t_acc: 0.25,
                risk_free_rate_r: 0.04,
                dividend_yield_q: 0.015,
                volatility_sigma: 0.22,
                smile_skew_slope_nu: -0.08,
                maturity_years_t: 1.5,
            };
            let r0 = evaluate_pechtl_haug_discrete_time_switch_accumulator_option(&base);
            assert!(
                r0.expected_in_the_money_occupation_time > 0.0
                    && r0.expected_in_the_money_occupation_time < base.maturity_years_t
            );
            assert!(r0.already_accumulated_present_value > 0.0);
            assert!(r0.time_switch_option_present_value > r0.already_accumulated_present_value);

            let inputs = [
                base.spot_price_s0,
                base.strike_price_k,
                base.accrual_rate_a,
                base.accumulated_occupation_time_t_acc,
                base.risk_free_rate_r,
                base.dividend_yield_q,
                base.volatility_sigma,
                base.smile_skew_slope_nu,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<9>::from_f64(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = PechtlHaugTimeSwitchSpec {
                spot_price_s0: dual_in[0],
                strike_price_k: dual_in[1],
                accrual_rate_a: dual_in[2],
                accumulated_occupation_time_t_acc: dual_in[3],
                risk_free_rate_r: dual_in[4],
                dividend_yield_q: dual_in[5],
                volatility_sigma: dual_in[6],
                smile_skew_slope_nu: dual_in[7],
                maturity_years_t: dual_in[8],
            };
            let r_ad = evaluate_pechtl_haug_discrete_time_switch_accumulator_option(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| PechtlHaugTimeSwitchSpec {
                    spot_price_s0: v[0],
                    strike_price_k: v[1],
                    accrual_rate_a: v[2],
                    accumulated_occupation_time_t_acc: v[3],
                    risk_free_rate_r: v[4],
                    dividend_yield_q: v[5],
                    volatility_sigma: v[6],
                    smile_skew_slope_nu: v[7],
                    maturity_years_t: v[8],
                };
                let rp = evaluate_pechtl_haug_discrete_time_switch_accumulator_option(&mk(p));
                let rm = evaluate_pechtl_haug_discrete_time_switch_accumulator_option(&mk(m));
                let fd = (rp.time_switch_option_present_value
                    - rm.time_switch_option_present_value)
                    / (2.0 * eps);
                assert!(
                    r_ad.time_switch_option_present_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.time_switch_option_present_value.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use pechtl_haug_time_switch_option_mod::*;

// ============================================================================
// Conze & Viswanathan (1991) & Broadie-Glasserman-Kou (1999) Lookback Extremum Option
// ============================================================================

pub mod conze_viswanathan_lookback_extremum_mod {
    use ql_ad::Real;
    use ql_math::distributions::norm_cdf;

    /// Specification for the Goldman, Sosin & Gatto (1979) / Conze & Viswanathan (1991, *JF* 46(5))
    /// floating-strike lookback call with Broadie, Glasserman & Kou (1999) discrete monitoring shift.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct ConzeViswanathanLookbackSpec<T: Real> {
        /// Current spot price $S_0 > 0$.
        pub spot_price_s0: T,
        /// Observed running minimum price to date $S_{\min} \le S_0$.
        pub observed_running_minimum_s_min: T,
        /// Risk-free interest rate $r > 0$.
        pub risk_free_rate_r: T,
        /// Continuous dividend yield $q \ge 0$.
        pub dividend_yield_q: T,
        /// Volatility $\sigma > 0$.
        pub volatility_sigma: T,
        /// Number of discrete monitoring dates $M \ge 1$.
        pub discrete_monitoring_steps_m: T,
        /// Remaining option maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Conze-Viswanathan / Broadie-Glasserman-Kou lookback engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct ConzeViswanathanLookbackAnalytics<T: Real> {
        /// Broadie-Glasserman-Kou continuity-adjusted running minimum $S_{\min}^{\text{adj}}$.
        pub continuity_adjusted_minimum_s_min: T,
        /// Continuous-monitoring Goldman-Sosin-Gatto / Conze-Viswanathan floating lookback call price.
        pub continuous_floating_lookback_call_price: T,
        /// Extremum reflection strikeless premium component.
        pub extremum_reflection_bonus_value: T,
        /// Discrete-monitoring floating-strike lookback call price $C_{\text{discrete}}(S_0, S_{\min}, T, M)$.
        pub discrete_floating_lookback_call_price: T,
    }

    fn floating_lookback_call_inner<T: Real>(s0: T, m_min: T, r: T, q: T, sig: T, t: T) -> (T, T) {
        let sig2 = sig * sig;
        let v_sqrt = sig * t.sqrt();
        let b = (r - q).max_val(T::from_f64(1e-4));

        let a1 = ((s0 / m_min).ln() + (b + T::from_f64(0.5) * sig2) * t) / v_sqrt;
        let a2 = a1 - v_sqrt;
        let a3 = ((s0 / m_min).ln() + (-b + T::from_f64(0.5) * sig2) * t) / v_sqrt;

        let df_r = (-r * t).exp();
        let df_q = (-q * t).exp();

        let vanilla_part = s0 * df_q * norm_cdf(a1) - m_min * df_r * norm_cdf(a2);
        let power_term = (s0 / m_min).pow(-T::from_f64(2.0) * b / sig2);
        let reflection_bonus = s0
            * df_r
            * (sig2 / (T::from_f64(2.0) * b))
            * (power_term * norm_cdf(-a3) - (b * t).exp() * norm_cdf(-a1));

        (vanilla_part + reflection_bonus, reflection_bonus)
    }

    /// Evaluates the Conze & Viswanathan (1991) / Goldman, Sosin & Gatto (1979) floating-strike
    /// lookback call price and Broadie, Glasserman & Kou (1999) / Haug (2007, Ch. 4.15.5) discrete
    /// monitoring correction.
    pub fn evaluate_conze_viswanathan_floating_and_fixed_lookback_extremum_option<T: Real>(
        spec: &ConzeViswanathanLookbackSpec<T>,
    ) -> ConzeViswanathanLookbackAnalytics<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
        let s_min = spec
            .observed_running_minimum_s_min
            .max_val(T::from_f64(1e-4))
            .min_val(s0 * T::from_f64(0.99));
        let sig = spec.volatility_sigma.max_val(T::from_f64(1e-4));
        let t = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let m_steps = spec.discrete_monitoring_steps_m.max_val(T::one());

        let (c_cont, bonus_cont) = floating_lookback_call_inner(
            s0,
            s_min,
            spec.risk_free_rate_r,
            spec.dividend_yield_q,
            sig,
            t,
        );

        // Broadie-Glasserman-Kou (1999) / Cheuk-Vorst (1997) continuity shift beta = -zeta(1/2)/sqrt(2*pi)
        let beta_bgk = T::from_f64(0.582_597_157_939_010_7);
        let shift = (beta_bgk * sig * (t / m_steps).sqrt()).exp();
        let s_min_adj = (s_min * shift).min_val(s0 * T::from_f64(0.999));

        let (c_shifted, _) = floating_lookback_call_inner(
            s0,
            s_min_adj,
            spec.risk_free_rate_r,
            spec.dividend_yield_q,
            sig,
            t,
        );
        let c_discrete = c_shifted / shift;

        ConzeViswanathanLookbackAnalytics {
            continuity_adjusted_minimum_s_min: s_min_adj,
            continuous_floating_lookback_call_price: c_cont,
            extremum_reflection_bonus_value: bonus_cont,
            discrete_floating_lookback_call_price: c_discrete,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn conze_viswanathan_lookback_invariants_and_ad_parity() {
            let base = ConzeViswanathanLookbackSpec {
                spot_price_s0: 100.0_f64,
                observed_running_minimum_s_min: 88.0,
                risk_free_rate_r: 0.05,
                dividend_yield_q: 0.015,
                volatility_sigma: 0.24,
                discrete_monitoring_steps_m: 52.0,
                maturity_years_t: 1.0,
            };
            let r0 = evaluate_conze_viswanathan_floating_and_fixed_lookback_extremum_option(&base);
            assert!(r0.continuity_adjusted_minimum_s_min > base.observed_running_minimum_s_min);
            assert!(r0.extremum_reflection_bonus_value > 0.0);
            assert!(
                r0.continuous_floating_lookback_call_price
                    > r0.discrete_floating_lookback_call_price
            );
            assert!(r0.discrete_floating_lookback_call_price > 0.0);

            let inputs = [
                base.spot_price_s0,
                base.observed_running_minimum_s_min,
                base.risk_free_rate_r,
                base.dividend_yield_q,
                base.volatility_sigma,
                base.discrete_monitoring_steps_m,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = ConzeViswanathanLookbackSpec {
                spot_price_s0: dual_in[0],
                observed_running_minimum_s_min: dual_in[1],
                risk_free_rate_r: dual_in[2],
                dividend_yield_q: dual_in[3],
                volatility_sigma: dual_in[4],
                discrete_monitoring_steps_m: dual_in[5],
                maturity_years_t: dual_in[6],
            };
            let r_ad =
                evaluate_conze_viswanathan_floating_and_fixed_lookback_extremum_option(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| ConzeViswanathanLookbackSpec {
                    spot_price_s0: v[0],
                    observed_running_minimum_s_min: v[1],
                    risk_free_rate_r: v[2],
                    dividend_yield_q: v[3],
                    volatility_sigma: v[4],
                    discrete_monitoring_steps_m: v[5],
                    maturity_years_t: v[6],
                };
                let rp =
                    evaluate_conze_viswanathan_floating_and_fixed_lookback_extremum_option(&mk(p));
                let rm =
                    evaluate_conze_viswanathan_floating_and_fixed_lookback_extremum_option(&mk(m));
                let fd = (rp.discrete_floating_lookback_call_price
                    - rm.discrete_floating_lookback_call_price)
                    / (2.0 * eps);
                assert!(
                    r_ad.discrete_floating_lookback_call_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.discrete_floating_lookback_call_price.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use conze_viswanathan_lookback_extremum_mod::*;
