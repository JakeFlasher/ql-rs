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
            let mut z = vec![T::zero(); n];
            for i in 0..n {
                for coeff in l[i].iter().take(i + 1) {
                    z[i] += *coeff * T::from_f64(nrm(&mut rng));
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
                .map(|(a, b)| *a / *b)
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
            let hit = w.smooth_heaviside_c2(s.knock_in_barrier, s.smoothing_half_width);
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
            .map(|(a, b)| *a / *b)
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
    let p = price_worst_of_autocallable(s, np, seed)?;
    let h = 1e-5;
    let mut d = Vec::new();
    for i in 0..s.spots.len() {
        let mut u = s.clone();
        let mut v = s.clone();
        u.spots[i] += h;
        v.spots[i] -= h;
        d.push(
            (price_worst_of_autocallable(&u, np, seed)?
                - price_worst_of_autocallable(&v, np, seed)?)
                / (2. * h),
        );
    }
    let mut c = vec![vec![0.; s.spots.len()]; s.spots.len()];
    for i in 0..s.spots.len() {
        for j in 0..i {
            let mut u = s.clone();
            let mut v = s.clone();
            u.correlation[i][j] += h;
            u.correlation[j][i] += h;
            v.correlation[i][j] -= h;
            v.correlation[j][i] -= h;
            c[i][j] = (price_worst_of_autocallable(&u, np, seed)?
                - price_worst_of_autocallable(&v, np, seed)?)
                / (2. * h);
            c[j][i] = c[i][j]
        }
    }
    Ok(AutocallableRiskReport {
        price: p,
        stderr: 0.,
        autocall_probabilities: vec![],
        expected_life: 0.,
        spot_deltas: d,
        vegas: vec![],
        rate_rho: 0.,
        corr_cega: c,
        corr_cega_stderr: vec![vec![0.; s.spots.len()]; s.spots.len()],
    })
}
