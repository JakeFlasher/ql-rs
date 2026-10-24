use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_mc::PhiloxRng;
use ql_models::HullWhite1F;
#[derive(Copy, Clone, Debug)]
pub struct Cheyette1F<T: Real> {
    pub mu: T,
}
impl<T: Real> Cheyette1F<T> {
    pub fn g_fn(&self, t: f64, big_t: f64) -> T {
        let tau = (big_t - t).max(0.);
        Self::g(self.mu, tau)
    }
    fn g(mu: T, tau: f64) -> T {
        if (mu.val() * tau).abs() < 1e-6 {
            T::from_f64(tau) - mu * tau * tau * 0.5 + mu * mu * tau * tau * tau / 6.
        } else {
            (T::one() - (-mu * tau).exp()) / mu
        }
    }
    pub fn discount_bond<C: YieldTermStructure<T>>(
        &self,
        curve: &C,
        t: f64,
        big_t: f64,
        x: T,
        y: T,
    ) -> T {
        let g = self.g_fn(t, big_t);
        curve.discount_factor(big_t) / curve.discount_factor(t)
            * (-g * x - T::from_f64(0.5) * g * g * y).exp()
    }
    pub fn explicit_local_variance_first_order(
        mu: T,
        k: T,
        w: T,
        dw_dt: T,
        dw_dk: T,
        d2w_dk2: T,
    ) -> T {
        let den = (T::one() - k * dw_dk / (T::from_f64(2.) * w)).powf(2.)
            + T::from_f64(0.5) * (d2w_dk2 - dw_dk * dw_dk / (T::from_f64(2.) * w));
        (dw_dt + mu * (T::from_f64(2.) * w - k * dw_dk) + w * dw_dk) / den
    }
    pub fn explicit_local_variance_third_order(
        mu: T,
        k: T,
        w: T,
        dw_dt: T,
        dw_dk: T,
        d2w_dk2: T,
    ) -> T {
        let base = Self::explicit_local_variance_first_order(mu, k, w, dw_dt, dw_dk, d2w_dk2);
        base + dw_dk * dw_dk * dw_dk / (w * 24.0).max_val(T::from_f64(1e-8))
    }
    pub fn simulate_short_rate_call(
        &self,
        expiry: f64,
        strike: f64,
        n_paths: usize,
        n_steps: usize,
        seed: u64,
        total_var_fn: impl Fn(f64, f64) -> (f64, f64, f64, f64),
    ) -> f64 {
        let rng = PhiloxRng::new(seed);
        let dt = expiry / n_steps.max(1) as f64;
        let mut sum = 0.;
        for p in 0..n_paths.max(1) {
            let sign = if p % 2 == 0 { 1. } else { -1. };
            let id = (p / 2) as u64;
            let mut x = 0.;
            let mut y = 0.;
            let mut int_x = 0.;
            for j in 0..n_steps.max(1) {
                let t = (j + 1) as f64 * dt;
                let (wv, dwv, dkv, d2v) = total_var_fn(t, x);
                let w = T::from_f64(wv.max(1e-10));
                let var = Self::explicit_local_variance_third_order(
                    self.mu,
                    T::from_f64(x),
                    w,
                    T::from_f64(dwv),
                    T::from_f64(dkv),
                    T::from_f64(d2v),
                )
                .val()
                .max(1e-10);
                let z = rng.normal4(id, j as u32)[0] * sign;
                int_x += x * dt;
                x += (y - self.mu.val() * x) * dt + var.sqrt() * dt.sqrt() * z;
                y += (var - 2. * self.mu.val() * y) * dt;
            }
            sum += (-int_x).exp() * (x - strike).max(0.);
        }
        sum / n_paths.max(1) as f64
    }
}

#[derive(Clone, Debug)]
pub struct CalibratedCoterminalCascade<T: Real = f64> {
    pub mu: T,
    pub exercise_times: Vec<f64>,
    pub piecewise_sigmas: Vec<T>,
    pub model_normal_vols: Vec<T>,
}
pub fn calibrate_cheyette_coterminal_ad<T: Real, C: YieldTermStructure<T>>(
    curve: &C,
    mu: T,
    swaptions: &[crate::EuropeanSwaption],
    market_normal_vols: &[T],
) -> Result<CalibratedCoterminalCascade<T>, String> {
    if swaptions.len() != market_normal_vols.len() {
        return Err("dimension mismatch".into());
    }
    if swaptions.windows(2).any(|w| w[1].expiry <= w[0].expiry)
        || swaptions.iter().any(|s| {
            !s.expiry.is_finite()
                || s.expiry <= 0.0
                || s.swap_times.len() < 2
                || s.swap_times.windows(2).any(|w| w[1] <= w[0])
        })
        || market_normal_vols
            .iter()
            .any(|v| !v.val().is_finite() || v.val() < 0.0)
    {
        return Err("invalid coterminal inputs".into());
    }
    let m = swaptions.len();
    let mu0 = mu.val();
    let mut piecewise_sigmas = Vec::with_capacity(m);
    let mut model_normal_vols = Vec::with_capacity(m);
    let mut prior_vx = T::zero();
    let mut prev_t = 0.0;
    // Bond loading and swap-rate loading at x=0.  This is the first-order
    // Cheyette normal-volatility map used by coterminal cascade calibration.
    let swap_loading = |s: &crate::EuropeanSwaption| {
        let t = s.expiry;
        let bond = |u: f64, x: f64| {
            if u <= t {
                1.0
            } else {
                let g = if mu0.abs() * (u - t) < 1e-8 {
                    u - t
                } else {
                    (1.0 - (-mu0 * (u - t)).exp()) / mu0
                };
                curve.discount_factor(u).val() / curve.discount_factor(t).val() * (-g * x).exp()
            }
        };
        let swap_rate = |x: f64| {
            let mut ann = 0.0;
            let mut floating = 0.0;
            for w in s.swap_times.windows(2) {
                let tau = w[1] - w[0];
                let p0 = bond(w[0], x);
                let p1 = bond(w[1], x);
                ann += tau * p1;
                floating += (p0 / p1 - 1.0) * p1;
            }
            if ann.abs() > 1e-14 {
                floating / ann
            } else {
                0.0
            }
        };
        let h = 1e-5;
        (swap_rate(h) - swap_rate(-h)) / (2.0 * h)
    };
    for (s, market) in swaptions.iter().zip(market_normal_vols.iter()) {
        let t = s.expiry;
        let loading = swap_loading(s).abs();
        if !loading.is_finite() || loading <= 1e-12 {
            return Err("singular swap-rate loading".into());
        }
        let target_vx = *market * *market * T::from_f64(t / (loading * loading));
        let dt_seg = t - prev_t;
        let dt_ad = T::from_f64(dt_seg);
        let seg = if mu0.abs() < 1e-8 {
            dt_ad
        } else {
            (T::one() - (-T::from_f64(2.0) * mu * dt_ad).exp()) / (T::from_f64(2.0) * mu)
        };
        if !seg.val().is_finite() || seg.val() <= 1e-14 {
            return Err("invalid Cheyette variance interval".into());
        }
        let decay = (-T::from_f64(2.0) * mu * dt_ad).exp();
        let decayed_prior = prior_vx * decay;
        let sigma2 = ((target_vx - decayed_prior) / seg).max_val(T::zero());
        let sigma = sigma2.sqrt();
        piecewise_sigmas.push(sigma);
        prior_vx = decayed_prior + sigma2 * seg;
        model_normal_vols
            .push(T::from_f64(loading) * (prior_vx / T::from_f64(t)).max_val(T::zero()).sqrt());
        prev_t = t;
    }
    Ok(CalibratedCoterminalCascade {
        mu,
        exercise_times: swaptions.iter().map(|s| s.expiry).collect(),
        piecewise_sigmas,
        model_normal_vols,
    })
}
pub fn price_callable_range_accrual_hw1f<
    T: Real,
    D: YieldTermStructure<T>,
    F: YieldTermStructure<T>,
>(
    _model: &HullWhite1F<T>,
    ois: &D,
    fwd: &F,
    spec: &crate::CallableRangeAccrualSpec<T>,
    steps: usize,
) -> (T, T, T) {
    if spec.call_period_indices.is_empty() {
        return (T::zero(), T::zero(), T::zero());
    }
    // The schedule stores actual period indices; using the loop ordinal here
    // silently shifted all callable dates whenever indices were sparse.  A
    // modest sub-grid average makes the corridor indicator less sensitive to
    // the arbitrary midpoint chosen for a long period while retaining AD.
    let sub = steps.max(1);
    let mut non = T::zero();
    for &period_idx in &spec.call_period_indices {
        let start = period_idx as f64;
        let end = start + 1.0;
        let mut accrual_gate = T::zero();
        for q in 0..sub {
            let t = start + (q as f64 + 0.5) / sub as f64;
            let dt = 1.0 / sub as f64;
            let f = (fwd.discount_factor(t - 0.5) / fwd.discount_factor(t) - T::one()) * 2.0;
            let gate = (f - spec.lower_bound)
                .smooth_heaviside_c2(T::zero(), spec.smoothing_eps.max(1e-6))
                * (T::one()
                    - (f - spec.upper_bound)
                        .smooth_heaviside_c2(T::zero(), spec.smoothing_eps.max(1e-6)));
            accrual_gate += gate * dt;
        }
        non += spec.notional * spec.coupon_rate * accrual_gate * ois.discount_factor(end);
    }
    let max_coupon = spec.notional * spec.coupon_rate * spec.call_period_indices.len() as f64;
    let callable = non.min_val(max_coupon.max_val(T::zero()));
    (callable, non, non - callable)
}
