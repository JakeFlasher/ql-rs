use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_mc::PhiloxRng;
use ql_models::HullWhite1F;
use ql_vol::bachelier_price;

#[derive(Clone, Debug)]
pub struct PiterbargQuasiGaussianAmortizingSwaptionSpec<T: Real = f64> {
    pub option_expiry: f64,
    pub second_exercise_time: Option<f64>,
    pub swap_times: Vec<f64>,
    pub period_notionals: Vec<T>,
    pub strike: T,
    pub mean_reversion: T,
    pub sigma0: T,
    pub local_vol_skew_b: T,
    pub cms_pay_delay: f64,
    pub is_payer: bool,
}
#[derive(Copy, Clone, Debug)]
pub struct PiterbargQuasiGaussianSwaptionValuation<T: Real = f64> {
    pub european_amortizing_npv: T,
    pub bermudan_amortizing_npv: T,
    pub amortizing_forward_swap_rate: T,
    pub amortizing_annuity: T,
    pub effective_normal_vol: T,
    pub multi_curve_cms_rate: T,
    pub cms_convexity_adjustment: T,
}

pub fn price_piterbarg_quasi_gaussian_amortizing_swaption<
    T: Real,
    D: YieldTermStructure<T>,
    F: YieldTermStructure<T>,
>(
    spec: &PiterbargQuasiGaussianAmortizingSwaptionSpec<T>,
    ois_discount: &D,
    projection_fwd: &F,
) -> PiterbargQuasiGaussianSwaptionValuation<T> {
    let m = spec
        .period_notionals
        .len()
        .min(spec.swap_times.len().saturating_sub(1));
    let te = spec.option_expiry.max(1e-10);
    let kappa = spec.mean_reversion.max_val(T::from_f64(1e-8));
    let mut ann = T::zero();
    let mut flt = T::zero();
    for j in 0..m {
        let t0 = spec.swap_times[j];
        let t1 = spec.swap_times[j + 1];
        let tau = t1 - t0;
        let n = spec.period_notionals[j];
        let pd1 = ois_discount.discount_factor(t1);
        let pf0 = projection_fwd.discount_factor(t0);
        let pf1 = projection_fwd.discount_factor(t1);
        let l = (pf0 / pf1 - T::one()) / tau;
        ann += n * tau * pd1;
        flt += n * tau * pd1 * l;
    }
    let s = flt / ann.max_val(T::from_f64(1e-12));
    let mut qsum = T::zero();
    let mut q2 = T::zero();
    for j in 0..m {
        let t0 = spec.swap_times[j];
        let t1 = spec.swap_times[j + 1];
        let tau = t1 - t0;
        let n = spec.period_notionals[j];
        let pd0 = ois_discount.discount_factor(t0);
        let pd1 = ois_discount.discount_factor(t1);
        let pf0 = projection_fwd.discount_factor(t0);
        let pf1 = projection_fwd.discount_factor(t1);
        let beta = (pf0 / pf1) / (pd0 / pd1);
        let g0 = (T::one() - (-kappa * (t0 - te).max(0.0)).exp()) / kappa;
        let g1 = (T::one() - (-kappa * (t1 - te).max(0.0)).exp()) / kappa;
        qsum += n * (beta * pd0 * g0 - pd1 * g1 + T::from_f64(tau) * s * pd1 * g1);
        q2 += n * (-beta * pd0 * g0 * g0 + pd1 * g1 * g1 * (T::one() - T::from_f64(tau) * s));
    }
    let q = qsum / ann.max_val(T::from_f64(1e-12));
    let mut g_ann = T::zero();
    for j in 0..m {
        let t1 = spec.swap_times[j + 1];
        let tau = spec.swap_times[j + 1] - spec.swap_times[j];
        let n = spec.period_notionals[j];
        let pd = ois_discount.discount_factor(t1);
        let g = (T::one() - (-kappa * (t1 - te).max(0.0)).exp()) / kappa;
        g_ann += n * T::from_f64(tau) * pd * g;
    }
    let g_ann = g_ann / ann.max_val(T::from_f64(1e-12));
    let qpp = q2 / ann.max_val(T::from_f64(1e-12)) + q * g_ann * 2.0;
    let varfac = ((T::one() - (-kappa * te * 2.0).exp()) / (kappa * te * 2.0))
        .max_val(T::from_f64(1e-12))
        .sqrt();
    let vol = (q.abs()
        * spec.sigma0
        * varfac
        * (T::one()
            + (spec.local_vol_skew_b + qpp / (q.abs().max_val(T::from_f64(1e-10)) * 2.0))
                * (spec.strike - s)
                * 0.5))
        .max_val(T::from_f64(1e-10));
    let npv = ann
        * bachelier_price(
            s,
            spec.strike,
            vol,
            T::from_f64(te),
            T::one(),
            spec.is_payer,
        );
    let pay = spec.swap_times.first().copied().unwrap_or(te) + spec.cms_pay_delay;
    let gpay = (T::one() - (-kappa * (pay - te).max(0.0)).exp()) / kappa;
    let mut alpha_num = T::zero();
    for j in 0..m {
        let t1 = spec.swap_times[j + 1];
        let tau = spec.swap_times[j + 1] - spec.swap_times[j];
        let n = spec.period_notionals[j];
        let pd = ois_discount.discount_factor(t1);
        let gj = (T::one() - (-kappa * (t1 - te).max(0.0)).exp()) / kappa;
        alpha_num += n * tau * pd * (gj - gpay);
    }
    let alpha = alpha_num / (ann.max_val(T::from_f64(1e-12)) * q.abs().max_val(T::from_f64(1e-10)));
    let cms_adj = alpha * vol * vol * te;
    let cms = s + cms_adj;
    let berm = if let Some(t2) = spec.second_exercise_time {
        if t2 <= te {
            npv
        } else {
            let tail: Vec<usize> = (0..m).filter(|&j| spec.swap_times[j + 1] > t2).collect();
            if tail.is_empty() {
                npv
            } else {
                let mut a2 = T::zero();
                let mut f2 = T::zero();
                for &j in &tail {
                    let tau = spec.swap_times[j + 1] - spec.swap_times[j];
                    let n = spec.period_notionals[j];
                    let pd = ois_discount.discount_factor(spec.swap_times[j + 1]);
                    let pf0 = projection_fwd.discount_factor(spec.swap_times[j]);
                    let pf1 = projection_fwd.discount_factor(spec.swap_times[j + 1]);
                    let l = (pf0 / pf1 - T::one()) / tau;
                    a2 += n * T::from_f64(tau) * pd;
                    f2 += n * T::from_f64(tau) * pd * l;
                }
                let s2 = f2 / a2.max_val(T::from_f64(1e-12));
                let sigma_n2 = vol;
                let rho12 = (-kappa * (t2 - te)).exp();
                let rem_t = (t2 - te).max(1e-12);
                let rem_var =
                    (sigma_n2 * sigma_n2 * (T::from_f64(t2) - rho12 * rho12 * T::from_f64(te)))
                        .max_val(T::from_f64(1e-12));
                let rem_vol = rem_var.sqrt() / T::from_f64(rem_t).sqrt();
                let nodes = [
                    -3.1909932018,
                    -2.2665805845,
                    -1.4685532892,
                    -0.7235510188,
                    0.0,
                    0.7235510188,
                    1.4685532892,
                    2.2665805845,
                    3.1909932018,
                ];
                let weights = [
                    0.0003960698,
                    0.0149487324,
                    0.1187461380,
                    0.3740123450,
                    0.4929700840,
                    0.3740123450,
                    0.1187461380,
                    0.0149487324,
                    0.0003960698,
                ];
                let weight_norm = weights.iter().sum::<f64>();
                let gh_sum = nodes
                    .iter()
                    .zip(weights.iter())
                    .fold(T::zero(), |acc, (&z, &w)| {
                        let zt = T::from_f64(z);
                        let s1 = s + vol * T::from_f64(te).sqrt() * zt;
                        let immediate = ann
                            * (if spec.is_payer {
                                s1 - spec.strike
                            } else {
                                spec.strike - s1
                            })
                            .max_val(T::zero());
                        let s2z = s2 + rho12 * sigma_n2 * T::from_f64(te).sqrt() * zt;
                        let continuation = a2
                            * bachelier_price(
                                s2z,
                                spec.strike,
                                rem_vol,
                                T::from_f64(rem_t),
                                T::one(),
                                spec.is_payer,
                            );
                        let diff = immediate - continuation;
                        let smooth_max = continuation
                            + (diff * diff + T::from_f64(1e-10)).sqrt() * 0.5
                            + diff * 0.5;
                        acc + T::from_f64(w / weight_norm) * smooth_max
                    });
                gh_sum.max_val(npv)
            }
        }
    } else {
        npv
    };
    PiterbargQuasiGaussianSwaptionValuation {
        european_amortizing_npv: npv,
        bermudan_amortizing_npv: berm,
        amortizing_forward_swap_rate: s,
        amortizing_annuity: ann,
        effective_normal_vol: vol,
        multi_curve_cms_rate: cms,
        cms_convexity_adjustment: cms_adj,
    }
}
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
    // Bond and swap-rate loading at x=0.  Keep the whole finite-difference
    // stencil in `T`: the curve may itself be a Dual/Tape curve and detaching
    // its discount factors here would silently erase quote sensitivities.
    let swap_loading = |s: &crate::EuropeanSwaption| -> T {
        let t = s.expiry;
        let model = Cheyette1F { mu };
        let bond = |u: f64, x: T| -> T {
            if u <= t {
                T::one()
            } else {
                let g = model.g_fn(t, u);
                curve.discount_factor(u) / curve.discount_factor(t) * (-g * x).exp()
            }
        };
        let swap_rate = |x: T| -> T {
            let mut ann = T::zero();
            let mut floating = T::zero();
            for w in s.swap_times.windows(2) {
                let tau = w[1] - w[0];
                let p0 = bond(w[0], x);
                let p1 = bond(w[1], x);
                ann += p1 * tau;
                floating += p0 - p1;
            }
            if ann.val().abs() > 1e-14 {
                floating / ann
            } else {
                T::zero()
            }
        };
        let h = 1e-5;
        let hv = T::from_f64(h);
        (swap_rate(hv) - swap_rate(-hv)) / (2.0 * h)
    };
    for (s, market) in swaptions.iter().zip(market_normal_vols.iter()) {
        let t = s.expiry;
        let loading = swap_loading(s).abs();
        if !loading.val().is_finite() || loading.val() <= 1e-12 {
            return Err("singular swap-rate loading".into());
        }
        let target_vx = *market * *market * T::from_f64(t) / (loading * loading);
        let dt_seg = t - prev_t;
        let dt_ad = T::from_f64(dt_seg);
        let seg = if (mu0 * dt_seg).abs() < 1e-6 {
            // (1-exp(-2 mu dt))/(2 mu) = dt-mu dt²+2/3 mu²dt³+...
            dt_ad * (T::one() - mu * dt_ad)
                + T::from_f64(2.0 / 3.0) * mu * mu * dt_ad * dt_ad * dt_ad
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
        model_normal_vols.push(loading * (prior_vx / T::from_f64(t)).max_val(T::zero()).sqrt());
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
    model: &HullWhite1F<T>,
    ois: &D,
    fwd: &F,
    spec: &crate::CallableRangeAccrualSpec<T>,
    steps: usize,
) -> (T, T, T) {
    if spec.call_period_indices.is_empty() {
        return (T::zero(), T::zero(), T::zero());
    }
    let sub = steps.max(1);
    let a = model.a.abs().max_val(T::from_f64(1e-6));
    let mut non = T::zero();
    let mut callable = T::zero();
    for &period_idx in &spec.call_period_indices {
        let start = period_idx as f64;
        let end = start + 1.0;
        let std_x = model.sigma.abs()
            * ((T::one() - (-a * (start + 0.5) * 2.0).exp()) / (a * 2.0))
                .max_val(T::zero())
                .sqrt();
        let df_end = ois.discount_factor(end);
        let call_strike = spec.notional * spec.coupon_rate * df_end * 0.5;
        let mut node_non = T::zero();
        let mut node_callable = T::zero();
        for (z, w) in [
            (-1.224744871391589, 1.0 / 6.0),
            (0.0, 4.0 / 6.0),
            (1.224744871391589, 1.0 / 6.0),
        ] {
            let mut accrual_gate = T::zero();
            for q in 0..sub {
                let t = start + (q as f64 + 0.5) / sub as f64;
                let dt = 1.0 / sub as f64;
                let fwd_rate =
                    (fwd.discount_factor(t - 0.5) / fwd.discount_factor(t) - T::one()) * 2.0;
                let f = fwd_rate + std_x * z;
                let gate = (f - spec.lower_bound)
                    .smooth_heaviside_c2(T::zero(), spec.smoothing_eps.max(1e-6))
                    * (T::one()
                        - (f - spec.upper_bound)
                            .smooth_heaviside_c2(T::zero(), spec.smoothing_eps.max(1e-6)));
                accrual_gate += gate * dt;
            }
            // Normalize the conditional Gaussian bond factor so its
            // expectation remains the deterministic OIS discount factor.
            // This is the exponential-martingale convexity term and keeps
            // the short-rate vega live at the central quadrature node.
            let state_df = (-std_x * z - T::from_f64(0.5) * std_x * std_x).exp();
            let cpn = spec.notional * spec.coupon_rate * accrual_gate * df_end * state_df;
            node_non += cpn * w;
            node_callable += cpn.min_val(call_strike) * w;
        }
        non += node_non;
        callable += node_callable;
    }
    let embedded = (non - callable).max_val(T::zero());
    (callable, non, embedded)
}

#[cfg(test)]
mod phase13_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_models::HullWhite1F;
    use ql_time::Date;
    #[test]
    fn cheyette_callable_has_positive_embedded_option_and_sigma_vega() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let c = FlatForwardCurve::new(d, Dual::<1>::from_value(0.03));
        let s = crate::CallableRangeAccrualSpec {
            notional: Dual::<1>::from_value(100.),
            coupon_rate: Dual::<1>::from_value(0.05),
            lower_bound: Dual::<1>::from_value(0.01),
            upper_bound: Dual::<1>::from_value(0.05),
            call_period_indices: vec![1, 2],
            smoothing_eps: 1e-3,
        };
        let x = price_callable_range_accrual_hw1f(
            &HullWhite1F {
                a: Dual::<1>::from_value(0.05),
                sigma: Dual::<1>::variable(0.015, 0),
            },
            &c,
            &c,
            &s,
            4,
        );
        assert!(x.0.val().is_finite() && x.1.val().is_finite());
        assert!(x.0.val() < x.1.val());
        assert!(x.2.val() > 1e-2);
        assert!(x.2.d[0].abs() > 1e-4);
    }
}

#[derive(Clone, Debug)]
pub struct CheyetteSwaptionSpec<T: Real = f64> {
    pub notional: T,
    pub strike: T,
    pub option_expiry: f64,
    pub payment_times: Vec<f64>,
    pub mean_reversion: T,
    pub sigma0: T,
    pub skew_b: T,
    pub is_payer: bool,
}

#[derive(Copy, Clone, Debug)]
pub struct CheyetteSwaptionValuation<T: Real> {
    pub npv: T,
    pub forward_swap_rate: T,
    pub annuity: T,
    pub effective_normal_vol: T,
    pub effective_skew: T,
    pub implied_normal_vol: T,
}

pub fn price_cheyette_swaption_markov_projection<T: Real, C: YieldTermStructure<T>>(
    s: &CheyetteSwaptionSpec<T>,
    curve: &C,
) -> CheyetteSwaptionValuation<T> {
    let te = s.option_expiry.max(0.0);
    if te <= 0.0 || s.payment_times.is_empty() {
        return CheyetteSwaptionValuation {
            npv: T::zero(),
            forward_swap_rate: T::zero(),
            annuity: T::zero(),
            effective_normal_vol: T::zero(),
            effective_skew: T::zero(),
            implied_normal_vol: T::zero(),
        };
    }
    let chey = Cheyette1F {
        mu: s.mean_reversion,
    };
    let mut ann = T::zero();
    let mut prev_t = te;
    let mut last_t = te;
    for &t in &s.payment_times {
        if t <= te {
            continue;
        }
        let tau = (t - prev_t).max(0.0);
        ann += curve.discount_factor(t) * tau;
        prev_t = t;
        last_t = t;
    }
    if ann.val() <= 1e-14 || last_t <= te {
        return CheyetteSwaptionValuation {
            npv: T::zero(),
            forward_swap_rate: T::zero(),
            annuity: T::zero(),
            effective_normal_vol: T::zero(),
            effective_skew: T::zero(),
            implied_normal_vol: T::zero(),
        };
    }
    let p_te = curve.discount_factor(te);
    let p_tm = curve.discount_factor(last_t);
    let s0 = (p_te - p_tm) / ann;
    let mut wg_sum = T::zero();
    let mut wg2_sum = T::zero();
    prev_t = te;
    for &t in &s.payment_times {
        if t <= te {
            continue;
        }
        let tau = (t - prev_t).max(0.0);
        let wj = curve.discount_factor(t) * tau / ann;
        let gj = chey.g_fn(te, t);
        wg_sum += wj * gj;
        wg2_sum += wj * gj * gj;
        prev_t = t;
    }
    let g_tm = chey.g_fn(te, last_t);
    let q0 = (p_tm / ann * g_tm + s0 * wg_sum).max_val(T::from_f64(1e-10));
    let k2_te = Cheyette1F::g(s.mean_reversion * 2.0, te);
    let var_n = s.sigma0 * s.sigma0 * q0 * q0 * k2_te;
    let eff_vol = (var_n / te).max_val(T::from_f64(1e-16)).sqrt();
    let kt = s.mean_reversion * te;
    let time_weight = if kt.val().abs() < 1e-2 {
        T::from_f64(0.5) + kt / 6.0 - kt * kt / 24.0 - kt * kt * kt / 72.0
    } else {
        let e1 = (-kt).exp();
        let e2 = e1 * e1;
        let e3 = e2 * e1;
        let num = (T::one() - e3) / 3.0 - e2 * (T::one() - e1);
        let den = (T::one() - e2) * (T::one() - e2);
        num * 2.0 / den.max_val(T::from_f64(1e-12))
    };
    let gamma_bond = ((p_tm / ann * g_tm * g_tm - s0 * wg2_sum) / q0 + wg_sum * 2.0) / q0;
    let eff_skew = -s.skew_b * time_weight / q0 + gamma_bond * 0.25;
    let dk = s.strike - s0;
    let imp_vol = (eff_vol
        * (T::one() + eff_skew * dk * 0.5 + eff_skew * eff_skew * dk * dk / 12.0))
        .max_val(T::from_f64(1e-8));
    let st = imp_vol * te.sqrt();
    let omega = if s.is_payer { 1.0 } else { -1.0 };
    let d = (s0 - s.strike) / st;
    let unit = T::from_f64(omega) * (s0 - s.strike) * (T::from_f64(omega) * d).norm_cdf()
        + st * d.norm_pdf();
    CheyetteSwaptionValuation {
        npv: s.notional * ann * unit,
        forward_swap_rate: s0,
        annuity: ann,
        effective_normal_vol: eff_vol,
        effective_skew: eff_skew,
        implied_normal_vol: imp_vol,
    }
}

#[cfg(test)]
mod phase16_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;
    #[test]
    fn cheyette_swaption_parity_and_ad() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let c = FlatForwardCurve::new(d, 0.03);
        let mut s = CheyetteSwaptionSpec {
            notional: 100.0,
            strike: 0.032,
            option_expiry: 1.0,
            payment_times: vec![2.0, 3.0, 4.0, 5.0],
            mean_reversion: 0.05,
            sigma0: 0.008,
            skew_b: -1.5,
            is_payer: true,
        };
        let pay = price_cheyette_swaption_markov_projection(&s, &c);
        s.is_payer = false;
        let rec = price_cheyette_swaption_markov_projection(&s, &c);
        let intrinsic = s.notional * pay.annuity * (pay.forward_swap_rate - s.strike);
        assert!((pay.npv - rec.npv - intrinsic).abs() < 1e-10);
        let h = 1e-5;
        let eval_s = |sig: f64| {
            let mut sp = s.clone();
            sp.is_payer = true;
            sp.sigma0 = sig;
            price_cheyette_swaption_markov_projection(&sp, &c).npv
        };
        let fd = (eval_s(0.008 + h) - eval_s(0.008 - h)) / (2.0 * h);
        type D = Dual<1>;
        let cd = FlatForwardCurve::new(d, D::from_value(0.03));
        let sd = CheyetteSwaptionSpec {
            notional: D::from_value(100.0),
            strike: D::from_value(0.032),
            option_expiry: 1.0,
            payment_times: vec![2.0, 3.0, 4.0, 5.0],
            mean_reversion: D::from_value(0.05),
            sigma0: D::variable(0.008, 0),
            skew_b: D::from_value(-1.5),
            is_payer: true,
        };
        let ad = price_cheyette_swaption_markov_projection(&sd, &cd).npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase21_tests {
    use super::*;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;

    #[test]
    fn amortizing_bermudan_collapse_and_premium() {
        let curve = FlatForwardCurve::new(Date(0), 0.03);
        let base = PiterbargQuasiGaussianAmortizingSwaptionSpec {
            option_expiry: 1.0,
            second_exercise_time: None,
            swap_times: vec![1.0, 2.0, 3.0],
            period_notionals: vec![1.0, 0.9],
            strike: 0.03,
            mean_reversion: 0.05,
            sigma0: 0.01,
            local_vol_skew_b: 0.1,
            cms_pay_delay: 0.0,
            is_payer: true,
        };
        let european = price_piterbarg_quasi_gaussian_amortizing_swaption(&base, &curve, &curve);
        assert!(
            (european.bermudan_amortizing_npv - european.european_amortizing_npv).abs() < 1e-12
        );
        let mut berm_spec = base.clone();
        berm_spec.second_exercise_time = Some(2.0);
        let berm = price_piterbarg_quasi_gaussian_amortizing_swaption(&berm_spec, &curve, &curve);
        assert!(berm.bermudan_amortizing_npv >= berm.european_amortizing_npv - 1e-12);
    }
}
