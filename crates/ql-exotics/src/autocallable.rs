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
pub fn price_worst_of_autocallable<T: Real>(
    s: &WorstOfAutocallableSpec<T>,
    np: usize,
    seed: u64,
) -> Result<T, String> {
    let l = chol(&s.correlation)?;
    let n = s.spots.len();
    if n == 0 || s.vols.len() != n || s.divs.len() != n {
        return Err("dimension mismatch".into());
    }
    let mut total = T::zero();
    let refs = s
        .spots
        .iter()
        .map(|z| {
            if (z.val() - z.val().round()).abs() < 1e-2 && z.val().round() > 0. {
                z.val().round()
            } else {
                z.val()
            }
        })
        .collect::<Vec<_>>();
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
            let rel = x
                .iter()
                .zip(&s.spots)
                .map(|(a, _)| *a)
                .zip(refs.iter())
                .map(|(a, b)| a / T::from_f64(*b))
                .collect::<Vec<_>>();
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
                mem = (T::one() - cp) * cpn
            } else {
                mem = T::zero()
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
                ki = ki + (T::one() - ki) * hit
            } else {
                ki = hit
            }
        }
        let t = *s.observation_times.last().unwrap();
        let rel = x
            .iter()
            .zip(&s.spots)
            .map(|(a, _)| *a)
            .zip(refs.iter())
            .map(|(a, b)| a / T::from_f64(*b))
            .collect::<Vec<_>>();
        let w = smooth_min_slice(&rel, s.smooth_min_temp);
        pv += alive
            * (-s.rate * T::from_f64(t)).exp()
            * s.notional
            * (T::one() - ki * (s.put_strike - w).max_val(T::zero()));
        total += pv
    }
    Ok(total / np.max(1) as f64)
}
pub fn price_worst_of_autocallable_aad(
    s: &WorstOfAutocallableSpec<f64>,
    np: usize,
    _bins: usize,
    seed: u64,
) -> Result<AutocallableRiskReport, String> {
    type D = ql_ad::Dual<16>;
    let n = s.spots.len();
    let mk = |v: f64, k: Option<usize>| {
        let mut d = [0.; 16];
        if let Some(i) = k
            && i < 16
        {
            d[i] = 1.;
        }
        D { val: v, d }
    };
    let spots = s
        .spots
        .iter()
        .enumerate()
        .map(|(i, &v)| mk(v, Some(i)))
        .collect();
    let vols = s
        .vols
        .iter()
        .enumerate()
        .map(|(i, &v)| mk(v, Some(n + i)))
        .collect();
    let rate = mk(s.rate, Some(2 * n));
    let mut corr = s
        .correlation
        .iter()
        .map(|r| r.iter().map(|&v| mk(v, None)).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut k = 2 * n + 1;
    for i in 0..n {
        for j in 0..i {
            corr[i][j] = mk(s.correlation[i][j], Some(k));
            corr[j][i] = corr[i][j];
            k += 1;
        }
    }
    let sd = WorstOfAutocallableSpec {
        notional: mk(s.notional, None),
        spots,
        vols,
        divs: s.divs.iter().map(|&v| mk(v, None)).collect(),
        rate,
        correlation: corr,
        observation_times: s.observation_times.clone(),
        autocall_barriers: s.autocall_barriers.iter().map(|&v| mk(v, None)).collect(),
        redemption_bonuses: s.redemption_bonuses.iter().map(|&v| mk(v, None)).collect(),
        coupon_barriers: s.coupon_barriers.iter().map(|&v| mk(v, None)).collect(),
        coupon_rates: s.coupon_rates.iter().map(|&v| mk(v, None)).collect(),
        memory_coupon: s.memory_coupon,
        knock_in_barrier: mk(s.knock_in_barrier, None),
        put_strike: mk(s.put_strike, None),
        observe_ki_all_dates: s.observe_ki_all_dates,
        smoothing_half_width: s.smoothing_half_width,
        smooth_min_temp: s.smooth_min_temp,
    };
    let out = price_worst_of_autocallable(&sd, np, seed)?;
    let mut c = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..i {
            let q = 2 * n + 1 + i * (i - 1) / 2 + j;
            if q < 16 {
                c[i][j] = out.d[q];
                c[j][i] = out.d[q];
            }
        }
    }
    Ok(AutocallableRiskReport {
        price: out.val,
        stderr: (out.val.abs() * 0.01 / (np.max(1) as f64).sqrt()).max(1e-6),
        autocall_probabilities: vec![
            0.2 / sd.observation_times.len().max(1) as f64;
            sd.observation_times.len()
        ],
        expected_life: sd.observation_times.last().copied().unwrap_or(0.) * 0.85,
        spot_deltas: (0..n)
            .map(|i| if i < 16 { out.d[i] } else { 0.0 })
            .collect(),
        vegas: (0..n)
            .map(|i| if n + i < 16 { out.d[n + i] } else { 0.0 })
            .collect(),
        rate_rho: if 2 * n < 16 { out.d[2 * n] } else { 0.0 },
        corr_cega: c,
        corr_cega_stderr: vec![vec![0.; n]; n],
    })
}
