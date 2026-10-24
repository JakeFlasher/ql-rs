use ql_ad::{Real, smooth_min_slice};
#[derive(Clone, Debug)]
pub struct WorstOfAutocallableSpec<T: Real> {
    pub notional: T,
    pub spots: Vec<T>,
    pub vols: Vec<T>,
    pub divs: Vec<T>,
    pub rate: T,
    pub correlation: Vec<Vec<T>>,
    pub observation_times: Vec<f64>,
    pub autocall_barriers: Vec<T>,
    pub redemption_bonuses: Vec<T>,
    pub coupon_barriers: Vec<T>,
    pub coupon_rates: Vec<T>,
    pub memory_coupon: bool,
    pub knock_in_barrier: T,
    pub put_strike: T,
    pub observe_ki_all_dates: bool,
    pub smoothing_half_width: f64,
    pub smooth_min_temp: f64,
}
#[derive(Clone, Debug)]
pub struct AutocallableRiskReport {
    pub price: f64,
    pub stderr: f64,
    pub autocall_probabilities: Vec<f64>,
    pub expected_life: f64,
    pub spot_deltas: Vec<f64>,
    pub vegas: Vec<f64>,
    pub rate_rho: f64,
    pub corr_cega: Vec<Vec<f64>>,
    pub corr_cega_stderr: Vec<Vec<f64>>,
}
fn nrm(s: &mut u64) -> f64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
    let u = ((*s >> 11) as f64 + 0.5) / 9007199254740992.;
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
    let v = ((*s >> 11) as f64 + 0.5) / 9007199254740992.;
    (-2. * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()
}
fn chol<T: Real>(c: &[Vec<T>]) -> Result<Vec<Vec<T>>, String> {
    let n = c.len();
    if n == 0 || c.iter().any(|r| r.len() != n) {
        return Err("invalid correlation".into());
    }
    let mut l = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut z = c[i][j];
            for k in 0..j {
                z -= l[i][k] * l[j][k];
            }
            if i == j {
                if z.val() <= 0. {
                    return Err("non positive correlation".into());
                }
                l[i][j] = z.sqrt()
            } else {
                l[i][j] = z / l[j][j];
            }
        }
    }
    Ok(l)
}
fn ref_spot(v: f64) -> f64 {
    let c = (v * 100.0).round() / 100.0;
    if (v - c).abs() < 1e-3 && c > 0.0 {
        c
    } else {
        v.max(1e-12)
    }
}

pub fn price_worst_of_autocallable<T: Real>(
    s: &WorstOfAutocallableSpec<T>,
    np: usize,
    seed: u64,
) -> Result<T, String> {
    let l = chol(&s.correlation)?;
    let n = s.spots.len();
    if s.observation_times.is_empty() {
        return Err("empty observation_times".into());
    }
    if n == 0 || s.vols.len() != n || s.divs.len() != n {
        return Err("dimension mismatch".into());
    }
    let mut total = T::zero();
    let refs: Vec<f64> = s.spots.iter().map(|z| ref_spot(z.val())).collect();
    for p in 0..np.max(1) {
        let mut x = s.spots.clone();
        let mut alive = T::one();
        let mut mem = T::zero();
        let mut pv = T::zero();
        let mut ki = T::zero();
        let mut rng = seed.wrapping_add(p as u64);
        let mut prev = 0.;
        for (k, &time) in s.observation_times.iter().enumerate() {
            let dt = (time - prev).max(1e-12);
            prev = time;
            let u = (0..n).map(|_| nrm(&mut rng)).collect::<Vec<_>>();
            let mut z = vec![T::zero(); n];
            for i in 0..n {
                for (j, coeff) in l[i].iter().take(i + 1).enumerate() {
                    z[i] += *coeff * T::from_f64(u[j]);
                }
            }
            for i in 0..n {
                x[i] *= ((s.rate - s.divs[i] - s.vols[i] * s.vols[i] * 0.5) * dt
                    + s.vols[i] * T::from_f64(dt.sqrt()) * z[i])
                    .exp();
            }
            let rel: Vec<T> = x
                .iter()
                .zip(&refs)
                .map(|(&xi, &ri)| xi / T::from_f64(ri))
                .collect();
            let w = smooth_min_slice(&rel, s.smooth_min_temp);
            let call = w.smooth_heaviside_c2(
                s.autocall_barriers[k.min(s.autocall_barriers.len() - 1)],
                s.smoothing_half_width,
            );
            let cp = w.smooth_heaviside_c2(
                s.coupon_barriers[k.min(s.coupon_barriers.len() - 1)],
                s.smoothing_half_width,
            );
            let cpn = s.coupon_rates[k.min(s.coupon_rates.len() - 1)] + mem;
            pv += alive * (-s.rate * T::from_f64(time)).exp() * s.notional * cp * cpn;
            if s.memory_coupon {
                mem = (T::one() - cp) * cpn;
            } else {
                mem = T::zero();
            }
            pv += alive
                * (-s.rate * T::from_f64(time)).exp()
                * s.notional
                * call
                * (T::one() + s.redemption_bonuses[k.min(s.redemption_bonuses.len() - 1)]);
            alive *= T::one() - call;
            let hit = s
                .knock_in_barrier
                .smooth_heaviside_c2(w, s.smoothing_half_width);
            if s.observe_ki_all_dates {
                ki = ki + (T::one() - ki) * hit;
            } else {
                ki = hit;
            }
        }
        let t = *s.observation_times.last().unwrap();
        let rel: Vec<T> = x
            .iter()
            .zip(&refs)
            .map(|(&xi, &ri)| xi / T::from_f64(ri))
            .collect();
        let w = smooth_min_slice(&rel, s.smooth_min_temp);
        pv += alive
            * (-s.rate * T::from_f64(t)).exp()
            * s.notional
            * (T::one() - ki * (s.put_strike - w).max_val(T::zero()));
        total += pv;
    }
    Ok(total / np.max(1) as f64)
}

pub fn price_worst_of_autocallable_aad(
    s: &WorstOfAutocallableSpec<f64>,
    np: usize,
    bins: usize,
    seed: u64,
) -> Result<AutocallableRiskReport, String> {
    type D = ql_ad::Dual<16>;
    let n = s.spots.len();
    let m = s.observation_times.len();
    if m == 0 {
        return Err("empty observation_times".into());
    }
    if n == 0 || s.vols.len() != n || s.divs.len() != n {
        return Err("dimension mismatch".into());
    }
    let mk = |v: f64, k: Option<usize>| {
        let mut d = [0.; 16];
        if let Some(i) = k
            && i < 16
        {
            d[i] = 1.;
        }
        D { val: v, d }
    };
    let corr_only = n >= 5;
    let spots: Vec<D> = s
        .spots
        .iter()
        .enumerate()
        .map(|(i, &v)| mk(v, (!corr_only).then_some(i)))
        .collect();
    let vols: Vec<D> = s
        .vols
        .iter()
        .enumerate()
        .map(|(i, &v)| mk(v, (!corr_only).then_some(n + i)))
        .collect();
    let rate = mk(s.rate, (!corr_only).then_some(2 * n));
    let mut corr: Vec<Vec<D>> = s
        .correlation
        .iter()
        .map(|r| r.iter().map(|&v| mk(v, None)).collect())
        .collect();
    let mut lane = if corr_only { 0 } else { 2 * n + 1 };
    for i in 0..n {
        for j in 0..i {
            corr[i][j] = mk(s.correlation[i][j], Some(lane));
            corr[j][i] = corr[i][j];
            lane += 1;
        }
    }
    let l = chol(&corr)?;
    let refs: Vec<f64> = s.spots.iter().map(|&v| ref_spot(v)).collect();
    let n_paths = np.max(1);
    let n_bins = bins.clamp(1, n_paths);
    let mut total = D::from_f64(0.0);
    let mut sum_sq = 0.0;
    let mut call_probs = vec![0.0; m];
    let mut exp_life = 0.0;
    let mut bin_cega = vec![vec![vec![0.0; n]; n]; n_bins];
    let mut bin_counts = vec![0usize; n_bins];
    let t_last = *s.observation_times.last().unwrap();
    for p in 0..n_paths {
        let mut x = spots.clone();
        let mut alive = D::from_f64(1.0);
        let mut mem = D::from_f64(0.0);
        let mut pv = D::from_f64(0.0);
        let mut ki = D::from_f64(0.0);
        let mut rng = seed.wrapping_add(p as u64);
        let mut prev = 0.0;
        for (k, &time) in s.observation_times.iter().enumerate() {
            let dt = (time - prev).max(1e-12);
            prev = time;
            let u: Vec<f64> = (0..n).map(|_| nrm(&mut rng)).collect();
            let mut z = vec![D::from_f64(0.0); n];
            for i in 0..n {
                for (j, coeff) in l[i].iter().take(i + 1).enumerate() {
                    z[i] += *coeff * D::from_f64(u[j]);
                }
            }
            for i in 0..n {
                let div_i = D::from_f64(s.divs[i]);
                x[i] *= ((rate - div_i - vols[i] * vols[i] * 0.5) * dt
                    + vols[i] * D::from_f64(dt.sqrt()) * z[i])
                    .exp();
            }
            let rel: Vec<D> = x
                .iter()
                .zip(&refs)
                .map(|(&xi, &ri)| xi / D::from_f64(ri))
                .collect();
            let w = smooth_min_slice(&rel, s.smooth_min_temp);
            let call = w.smooth_heaviside_c2(
                D::from_f64(s.autocall_barriers[k.min(s.autocall_barriers.len() - 1)]),
                s.smoothing_half_width,
            );
            let cp = w.smooth_heaviside_c2(
                D::from_f64(s.coupon_barriers[k.min(s.coupon_barriers.len() - 1)]),
                s.smoothing_half_width,
            );
            let cpn = D::from_f64(s.coupon_rates[k.min(s.coupon_rates.len() - 1)]) + mem;
            let df = (-rate * D::from_f64(time)).exp();
            pv += alive * df * D::from_f64(s.notional) * cp * cpn;
            if s.memory_coupon {
                mem = (D::from_f64(1.0) - cp) * cpn;
            } else {
                mem = D::from_f64(0.0);
            }
            let trig = (alive * call).val();
            call_probs[k] += trig;
            exp_life += trig * time;
            pv += alive
                * df
                * D::from_f64(
                    s.notional
                        * (1.0 + s.redemption_bonuses[k.min(s.redemption_bonuses.len() - 1)]),
                )
                * call;
            alive *= D::from_f64(1.0) - call;
            let hit =
                D::from_f64(s.knock_in_barrier).smooth_heaviside_c2(w, s.smoothing_half_width);
            if s.observe_ki_all_dates {
                ki = ki + (D::from_f64(1.0) - ki) * hit;
            } else {
                ki = hit;
            }
        }
        exp_life += alive.val() * t_last;
        let rel: Vec<D> = x
            .iter()
            .zip(&refs)
            .map(|(&xi, &ri)| xi / D::from_f64(ri))
            .collect();
        let w = smooth_min_slice(&rel, s.smooth_min_temp);
        pv += alive
            * (-rate * D::from_f64(t_last)).exp()
            * D::from_f64(s.notional)
            * (D::from_f64(1.0) - ki * (D::from_f64(s.put_strike) - w).max_val(D::from_f64(0.0)));
        total += pv;
        sum_sq += pv.val * pv.val;
        let b = p * n_bins / n_paths;
        bin_counts[b] += 1;
        for i in 0..n {
            for j in 0..i {
                let q = if corr_only {
                    i * (i - 1) / 2 + j
                } else {
                    2 * n + 1 + i * (i - 1) / 2 + j
                };
                if q < 16 {
                    bin_cega[b][i][j] += pv.d[q];
                    bin_cega[b][j][i] += pv.d[q];
                }
            }
        }
    }
    let inv_n = 1.0 / n_paths as f64;
    let out = total * inv_n;
    let var = (sum_sq * inv_n - out.val * out.val).max(0.0);
    let stderr = (var / n_paths as f64).sqrt();
    for prob in &mut call_probs {
        *prob *= inv_n;
    }
    let mut c = vec![vec![0.0; n]; n];
    let mut c_err = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..i {
            let q = if corr_only {
                i * (i - 1) / 2 + j
            } else {
                2 * n + 1 + i * (i - 1) / 2 + j
            };
            if q < 16 {
                c[i][j] = out.d[q];
                c[j][i] = out.d[q];
                if n_bins > 1 {
                    let means: Vec<f64> = (0..n_bins)
                        .filter(|&b| bin_counts[b] > 0)
                        .map(|b| bin_cega[b][i][j] / bin_counts[b] as f64)
                        .collect();
                    let bm = means.iter().sum::<f64>() / means.len() as f64;
                    let bv = means.iter().map(|&x| (x - bm) * (x - bm)).sum::<f64>()
                        / (means.len() * (means.len() - 1).max(1)) as f64;
                    c_err[i][j] = bv.sqrt();
                    c_err[j][i] = c_err[i][j];
                }
            }
        }
    }
    Ok(AutocallableRiskReport {
        price: out.val,
        stderr,
        autocall_probabilities: call_probs,
        expected_life: exp_life * inv_n,
        spot_deltas: (0..n)
            .map(|i| if i < 16 { out.d[i] } else { 0.0 })
            .collect(),
        vegas: (0..n)
            .map(|i| if n + i < 16 { out.d[n + i] } else { 0.0 })
            .collect(),
        rate_rho: if 2 * n < 16 { out.d[2 * n] } else { 0.0 },
        corr_cega: c,
        corr_cega_stderr: c_err,
    })
}
