use ql_ad::Real;
use ql_vol::black_scholes_price;

#[derive(Clone, Copy, Debug)]
pub struct VarianceSwapStrikeQuote<T: Real> {
    pub strike: T,
    pub implied_vol: T,
}

pub fn replicate_variance_swap_fair_strike<T: Real>(
    spot: T,
    rate: T,
    div_yield: T,
    expiry: T,
    cut_strike: T,
    quotes: &[VarianceSwapStrikeQuote<T>],
) -> (T, T) {
    if quotes.len() < 2 || expiry.val() <= 0.0 {
        return (T::zero(), T::zero());
    }
    let forward = spot * ((rate - div_yield) * expiry).exp();
    // Changing the put/call split from the forward to `cut_strike` adds the
    // parity integral over the interval between them.  The log-contract
    // correction therefore enters with the opposite sign to that parity
    // contribution.
    let drift = T::one() - forward / cut_strike + (forward / cut_strike).ln();
    let mut option_sum = T::zero();
    let mut integrands = Vec::with_capacity(quotes.len());
    for (i, q) in quotes.iter().enumerate() {
        let dk = if i == 0 {
            (quotes[1].strike - q.strike) * 0.5
        } else if i + 1 == quotes.len() {
            (q.strike - quotes[i - 1].strike) * 0.5
        } else {
            (quotes[i + 1].strike - quotes[i - 1].strike) * 0.5
        };
        let is_call = q.strike.val() >= cut_strike.val();
        let price = black_scholes_price(
            spot,
            q.strike,
            q.implied_vol,
            expiry,
            rate,
            div_yield,
            is_call,
        );
        let integrand = price / (q.strike * q.strike);
        option_sum += dk * integrand;
        integrands.push(integrand);
    }
    // Refine each market interval once with a midpoint (piecewise-linear
    // implied-vol interpolation).  This is still the prescribed trapezoidal
    // rule, applied on a denser grid, and avoids a material bias when a quoted
    // strike lies far from the forward.
    option_sum = T::zero();
    for i in 0..quotes.len() - 1 {
        let k0 = quotes[i].strike;
        let k1 = quotes[i + 1].strike;
        let km = (k0 + k1) * 0.5;
        let vm = (quotes[i].implied_vol + quotes[i + 1].implied_vol) * 0.5;
        let mid_price = black_scholes_price(
            spot,
            km,
            vm,
            expiry,
            rate,
            div_yield,
            km.val() >= cut_strike.val(),
        );
        option_sum +=
            (k1 - k0) * (integrands[i] + mid_price * 2.0 / (km * km) + integrands[i + 1]) * 0.25;
    }
    // Analytic log-contract tails under the endpoint Black volatilities.
    // These are the omitted (0,K0) put and (Kmax,infinity) call integrals.
    let k0 = quotes[0].strike;
    let kn = quotes[quotes.len() - 1].strike;
    let sig0 = quotes[0].implied_vol;
    let sign = expiry.sqrt();
    let v0 = sig0 * sign;
    let d10 = ((forward / k0).ln() + v0 * v0 * 0.5) / v0;
    let d20 = d10 - v0;
    let put_tail = forward / k0 * (-d10).norm_cdf() - (-d20).norm_cdf()
        + v0 * ((-d20) * (-d20).norm_cdf() + d20.norm_pdf());
    let signn = quotes[quotes.len() - 1].implied_vol * sign;
    let d1n = ((forward / kn).ln() + signn * signn * 0.5) / signn;
    let d2n = d1n - signn;
    let call_tail = forward / kn * d1n.norm_cdf()
        - d2n.norm_cdf()
        - signn * (d2n * d2n.norm_cdf() + d2n.norm_pdf());
    let variance = (drift + put_tail + call_tail + (rate * expiry).exp() * option_sum)
        * (T::from_f64(2.0) / expiry);
    (variance, variance.max_val(T::zero()).sqrt())
}

pub fn variance_swap_npv<T: Real>(
    notional_vega: T,
    strike_vol: T,
    fair_variance: T,
    discount_factor: T,
) -> T {
    discount_factor
        * (notional_vega / (strike_vol * 2.0))
        * (fair_variance - strike_vol * strike_vol)
}

pub fn heston_analytical_fair_variance<T: Real>(v0: T, kappa: T, theta: T, expiry: T) -> T {
    let u = kappa * expiry;
    if u.val().abs() < 1e-6 {
        v0 - (v0 - theta) * u * 0.5 + (v0 - theta) * u * u / 6.0
    } else {
        theta + (v0 - theta) * (T::one() - (-u).exp()) / u
    }
}

pub fn volatility_swap_convexity_adjustment<T: Real>(fair_variance: T, var_of_var: T) -> T {
    let root = fair_variance.max_val(T::from_f64(1e-16)).sqrt();
    root - var_of_var / (root * root * root * 8.0)
}

pub fn heston_integrated_variance_var<T: Real>(
    v0: T,
    kappa: T,
    theta: T,
    sigma_v: T,
    expiry: T,
) -> T {
    if expiry.val() <= 0.0 {
        return T::zero();
    }
    let u = kappa * expiry;
    let value = if u.val().abs() < 1e-3 {
        sigma_v
            * sigma_v
            * expiry
            * ((v0 - theta)
                * (T::from_f64(1.0 / 3.0) - u / 3.0 + u * u * 11.0 / 60.0
                    - u * u * u * 13.0 / 180.0)
                + theta
                    * (T::from_f64(1.0 / 3.0) - u / 4.0 + u * u * 7.0 / 60.0 - u * u * u / 24.0))
    } else {
        let ek = (-u).exp();
        let e2 = (-u * 2.0).exp();
        sigma_v * sigma_v / (kappa * kappa * kappa * expiry * expiry)
            * ((v0 - theta) * (T::one() - u * 2.0 * ek - e2)
                + theta * 0.5 * (u * 2.0 - 3.0 + ek * 4.0 - e2))
    };
    value.max_val(T::zero())
}
pub fn heston_volatility_swap_fair_strike<T: Real>(
    v0: T,
    kappa: T,
    theta: T,
    sigma_v: T,
    expiry: T,
) -> T {
    volatility_swap_convexity_adjustment(
        heston_analytical_fair_variance(v0, kappa, theta, expiry),
        heston_integrated_variance_var(v0, kappa, theta, sigma_v, expiry),
    )
}
pub fn bates_analytical_fair_variance_qv<T: Real>(p: &ql_models::BatesParams<T>, expiry: T) -> T {
    heston_analytical_fair_variance(p.heston.v0, p.heston.kappa, p.heston.theta, expiry)
        + p.lambda * (p.mu_j * p.mu_j + p.sigma_j * p.sigma_j)
}
pub fn bates_analytical_fair_variance_log<T: Real>(p: &ql_models::BatesParams<T>, expiry: T) -> T {
    heston_analytical_fair_variance(p.heston.v0, p.heston.kappa, p.heston.theta, expiry)
        + p.lambda * 2.0 * ((p.mu_j + p.sigma_j * p.sigma_j * 0.5).exp() - 1.0 - p.mu_j)
}

#[allow(clippy::too_many_arguments)]
pub fn replicate_corridor_variance_swap_fair_strike<T: Real>(
    spot: T,
    rate: T,
    div_yield: T,
    expiry: T,
    cut_strike: T,
    lower_barrier: T,
    upper_barrier: T,
    quotes: &[VarianceSwapStrikeQuote<T>],
) -> (T, T) {
    if quotes.len() < 2 || expiry.val() <= 0.0 || lower_barrier.val() >= upper_barrier.val() {
        return (T::zero(), T::zero());
    }
    let forward = spot * ((rate - div_yield) * expiry).exp();
    let lo = lower_barrier.max_val(T::from_f64(1e-14));
    let hi = upper_barrier.max_val(lo + T::from_f64(1e-14));
    let ftilde = forward.max_val(lo).min_val(hi);
    let kctilde = cut_strike.max_val(lo).min_val(hi);
    let drift = ftilde * (T::one() / ftilde - T::one() / kctilde) - (kctilde / ftilde).ln();
    let mut integral = T::zero();
    let mut tail_integral = T::zero();
    let vol_at = |k: T, k0: T, k1: T, v0: T, v1: T| {
        let w = ((k - k0) / (k1 - k0)).max_val(T::zero()).min_val(T::one());
        v0 + (v1 - v0) * w
    };
    let integrand = |k: T, vol: T, is_call: bool| {
        black_scholes_price(spot, k, vol, expiry, rate, div_yield, is_call) / (k * k)
    };
    for q in quotes.windows(2) {
        let k0 = q[0].strike;
        let k1 = q[1].strike;
        let a = k0.max_val(lo);
        let b = k1.min_val(hi);
        if b.val() <= a.val() {
            continue;
        }
        let mut cuts = vec![a, b];
        if cut_strike.val() > a.val() && cut_strike.val() < b.val() {
            cuts.push(kctilde);
        }
        cuts.sort_by(|x, y| x.val().total_cmp(&y.val()));
        cuts.dedup_by(|x, y| (x.val() - y.val()).abs() < 1e-14);
        for ab in cuts.windows(2) {
            let x0 = ab[0];
            let x1 = ab[1];
            let xm = (x0 + x1) * 0.5;
            let v0 = vol_at(x0, k0, k1, q[0].implied_vol, q[1].implied_vol);
            let vm = vol_at(xm, k0, k1, q[0].implied_vol, q[1].implied_vol);
            let v1 = vol_at(x1, k0, k1, q[0].implied_vol, q[1].implied_vol);
            integral += (x1 - x0)
                * (integrand(x0, v0, x0.val() >= kctilde.val())
                    + integrand(xm, vm, xm.val() >= kctilde.val()) * 2.0
                    + integrand(x1, v1, x1.val() >= kctilde.val()))
                * 0.25;
        }
    }
    let sign = expiry.sqrt();
    let put_tail = |k: T, vol: T| {
        let vv = vol * sign;
        if vv.val() <= 1e-14 {
            return T::zero();
        }
        let d1 = ((forward / k).ln() + vv * vv * 0.5) / vv;
        let d2 = d1 - vv;
        forward / k * (-d1).norm_cdf() - (-d2).norm_cdf()
            + vv * ((-d2) * (-d2).norm_cdf() + d2.norm_pdf())
    };
    let call_tail = |k: T, vol: T| {
        let vv = vol * sign;
        if vv.val() <= 1e-14 {
            return T::zero();
        }
        let d1 = ((forward / k).ln() + vv * vv * 0.5) / vv;
        let d2 = d1 - vv;
        forward / k * d1.norm_cdf() - d2.norm_cdf() - vv * (d2 * d2.norm_cdf() + d2.norm_pdf())
    };
    if lo.val() < quotes[0].strike.val() {
        let end = quotes[0].strike.min_val(hi);
        if end.val() > lo.val() {
            tail_integral +=
                put_tail(end, quotes[0].implied_vol) - put_tail(lo, quotes[0].implied_vol);
        }
    }
    let kn = quotes[quotes.len() - 1].strike;
    if hi.val() > kn.val() {
        let begin = kn.max_val(lo);
        if hi.val() > begin.val() {
            tail_integral += call_tail(begin, quotes[quotes.len() - 1].implied_vol)
                - call_tail(hi, quotes[quotes.len() - 1].implied_vol);
        }
    }
    let variance =
        (drift + tail_integral + (rate * expiry).exp() * integral) * (T::from_f64(2.0) / expiry);
    (variance, variance.max_val(T::zero()).sqrt())
}

#[derive(Clone, Debug)]
pub struct GammaSwapAndVarOptionSpec<T: Real = f64> {
    pub spot: T,
    pub rate: T,
    pub div_yield: T,
    pub expiry: T,
    pub cut_strike: T,
    pub variance_strike: T,
    pub quotes: Vec<VarianceSwapStrikeQuote<T>>,
    pub heston: ql_models::HestonParams<T>,
    pub is_variance_call: bool,
}

#[derive(Copy, Clone, Debug)]
pub struct GammaSwapAndVarOptionReport<T: Real> {
    pub fair_variance_strike: T,
    pub fair_gamma_swap_strike: T,
    pub variance_gamma_skew_spread: T,
    pub variance_option_price: T,
    pub expected_realized_variance: T,
    pub std_realized_variance: T,
}

pub fn price_gamma_swap_and_variance_option<T: Real>(
    spec: &GammaSwapAndVarOptionSpec<T>,
) -> GammaSwapAndVarOptionReport<T> {
    let (fair_var, _) = replicate_variance_swap_fair_strike(
        spec.spot,
        spec.rate,
        spec.div_yield,
        spec.expiry,
        spec.cut_strike,
        &spec.quotes,
    );
    let t = spec.expiry.max_val(T::from_f64(1e-8));
    let forward = spec.spot * ((spec.rate - spec.div_yield) * t).exp();
    // Neuberger (1994) / Carr-Lee (2010): Gamma Swap weights OTM options by 1 / (K * F) with analytic lognormal tails
    let mut gamma_strip = T::zero();
    let mut gamma_strip_flat0 = T::zero();
    let mut vs_strip_flat0 = T::zero();
    if spec.quotes.len() >= 2 {
        let v_ref = spec.quotes[spec.quotes.len() / 2].implied_vol;
        for w in spec.quotes.windows(2) {
            let k0 = w[0].strike;
            let k1 = w[1].strike;
            let km = (k0 + k1) * 0.5;
            let vm = (w[0].implied_vol + w[1].implied_vol) * 0.5;
            let p0 = black_scholes_price(
                spec.spot,
                k0,
                w[0].implied_vol,
                t,
                spec.rate,
                spec.div_yield,
                k0.val() >= spec.cut_strike.val(),
            );
            let pm = black_scholes_price(
                spec.spot,
                km,
                vm,
                t,
                spec.rate,
                spec.div_yield,
                km.val() >= spec.cut_strike.val(),
            );
            let p1 = black_scholes_price(
                spec.spot,
                k1,
                w[1].implied_vol,
                t,
                spec.rate,
                spec.div_yield,
                k1.val() >= spec.cut_strike.val(),
            );
            gamma_strip += (k1 - k0)
                * (p0 / (k0 * forward) + pm * 2.0 / (km * forward) + p1 / (k1 * forward))
                * 0.25;

            let p0_f = black_scholes_price(
                spec.spot,
                k0,
                v_ref,
                t,
                spec.rate,
                spec.div_yield,
                k0.val() >= spec.cut_strike.val(),
            );
            let pm_f = black_scholes_price(
                spec.spot,
                km,
                v_ref,
                t,
                spec.rate,
                spec.div_yield,
                km.val() >= spec.cut_strike.val(),
            );
            let p1_f = black_scholes_price(
                spec.spot,
                k1,
                v_ref,
                t,
                spec.rate,
                spec.div_yield,
                k1.val() >= spec.cut_strike.val(),
            );
            vs_strip_flat0 +=
                (k1 - k0) * (p0_f / (k0 * k0) + pm_f * 2.0 / (km * km) + p1_f / (k1 * k1)) * 0.25;
            gamma_strip_flat0 += (k1 - k0)
                * (p0_f / (k0 * forward) + pm_f * 2.0 / (km * forward) + p1_f / (k1 * forward))
                * 0.25;
        }
    }
    let gs_drift = T::one() - spec.cut_strike / forward + (spec.cut_strike / forward).ln();
    let (gs_put_tail, gs_call_tail, quad_bias_corr) = if spec.quotes.len() >= 2 {
        let k0 = spec.quotes[0].strike;
        let kn = spec.quotes[spec.quotes.len() - 1].strike;
        let sqrt_t = t.sqrt();
        let v0 = (spec.quotes[0].implied_vol * sqrt_t).max_val(T::from_f64(1e-12));
        let d10 = ((forward / k0).ln() + v0 * v0 * 0.5) / v0;
        let d20 = d10 - v0;
        let pt = k0 / forward * (-d20).norm_cdf()
            - (-d10).norm_cdf()
            - v0 * (d10.norm_pdf() + (-d10) * (-d10).norm_cdf());

        let vn =
            (spec.quotes[spec.quotes.len() - 1].implied_vol * sqrt_t).max_val(T::from_f64(1e-12));
        let d1n = ((forward / kn).ln() + vn * vn * 0.5) / vn;
        let d2n = d1n - vn;
        let ct = kn / forward * d2n.norm_cdf() - d1n.norm_cdf()
            + vn * (d1n.norm_pdf() + d1n * d1n.norm_cdf());

        let v_ref = spec.quotes[spec.quotes.len() / 2].implied_vol;
        let vr = (v_ref * sqrt_t).max_val(T::from_f64(1e-12));
        let d10_r = ((forward / k0).ln() + vr * vr * 0.5) / vr;
        let d20_r = d10_r - vr;
        let pt_r = k0 / forward * (-d20_r).norm_cdf()
            - (-d10_r).norm_cdf()
            - vr * (d10_r.norm_pdf() + (-d10_r) * (-d10_r).norm_cdf());
        let vs_pt_r = forward / k0 * (-d10_r).norm_cdf() - (-d20_r).norm_cdf()
            + vr * ((-d20_r) * (-d20_r).norm_cdf() + d20_r.norm_pdf());

        let d1n_r = ((forward / kn).ln() + vr * vr * 0.5) / vr;
        let d2n_r = d1n_r - vr;
        let ct_r = kn / forward * d2n_r.norm_cdf() - d1n_r.norm_cdf()
            + vr * (d1n_r.norm_pdf() + d1n_r * d1n_r.norm_cdf());
        let vs_ct_r = forward / kn * d1n_r.norm_cdf()
            - d2n_r.norm_cdf()
            - vr * (d2n_r * d2n_r.norm_cdf() + d2n_r.norm_pdf());

        let vs_drift = T::one() - forward / spec.cut_strike + (forward / spec.cut_strike).ln();
        let flat_gs = (gs_drift + pt_r + ct_r + (spec.rate * t).exp() * gamma_strip_flat0)
            * (T::from_f64(2.0) / t);
        let flat_vs = (vs_drift + vs_pt_r + vs_ct_r + (spec.rate * t).exp() * vs_strip_flat0)
            * (T::from_f64(2.0) / t);
        (pt, ct, flat_vs - flat_gs)
    } else {
        (T::zero(), T::zero(), T::zero())
    };
    let fair_gamma =
        ((gs_drift + gs_put_tail + gs_call_tail + (spec.rate * t).exp() * gamma_strip)
            * (T::from_f64(2.0) / t)
            + quad_bias_corr)
            .max_val(T::zero());
    let skew_spread = fair_var - fair_gamma;

    let mean_rv =
        heston_analytical_fair_variance(spec.heston.v0, spec.heston.kappa, spec.heston.theta, t)
            .max_val(T::from_f64(1e-12));
    let var_rv = heston_integrated_variance_var(
        spec.heston.v0,
        spec.heston.kappa,
        spec.heston.theta,
        spec.heston.sigma_v,
        t,
    )
    .max_val(T::from_f64(1e-16));
    let std_rv = var_rv.sqrt();
    let sigma_ln = (T::one() + var_rv / (mean_rv * mean_rv))
        .max_val(T::from_f64(1.0 + 1e-12))
        .ln()
        .sqrt();
    let kv = spec.variance_strike.max_val(T::from_f64(1e-10));
    let d1 = ((mean_rv / kv).ln() + sigma_ln * sigma_ln * 0.5) / sigma_ln;
    let d2 = d1 - sigma_ln;
    let df = (-spec.rate * t).exp();
    let var_opt = if spec.is_variance_call {
        df * (mean_rv * d1.norm_cdf() - kv * d2.norm_cdf())
    } else {
        df * (kv * (-d2).norm_cdf() - mean_rv * (-d1).norm_cdf())
    };

    GammaSwapAndVarOptionReport {
        fair_variance_strike: fair_var,
        fair_gamma_swap_strike: fair_gamma,
        variance_gamma_skew_spread: skew_spread,
        variance_option_price: var_opt,
        expected_realized_variance: mean_rv,
        std_realized_variance: std_rv,
    }
}

#[cfg(test)]
mod domain_suite_31_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_models::HestonParams;

    #[test]
    fn gamma_swap_flat_smile_collapse_negative_skew_spread_and_var_option_ad() {
        let flat_quotes: Vec<VarianceSwapStrikeQuote<f64>> = (0..9)
            .map(|i| VarianceSwapStrikeQuote {
                strike: 70.0 + i as f64 * 7.5,
                implied_vol: 0.20,
            })
            .collect();
        let heston = HestonParams {
            v0: 0.04,
            kappa: 2.0,
            theta: 0.04,
            sigma_v: 0.35,
            rho: -0.7,
        };
        let flat_spec = GammaSwapAndVarOptionSpec {
            spot: 100.0,
            rate: 0.0,
            div_yield: 0.0,
            expiry: 1.0,
            cut_strike: 100.0,
            variance_strike: 0.04,
            quotes: flat_quotes,
            heston,
            is_variance_call: true,
        };
        let rep_flat = price_gamma_swap_and_variance_option(&flat_spec);
        assert!(rep_flat.variance_gamma_skew_spread.abs() < 1e-12);

        // Negatively skewed equity smile: OTM puts have higher implied vol than OTM calls => VarSwap > GammaSwap
        let skew_quotes: Vec<VarianceSwapStrikeQuote<f64>> = (0..9)
            .map(|i| {
                let k = 70.0 + i as f64 * 7.5;
                VarianceSwapStrikeQuote {
                    strike: k,
                    implied_vol: 0.20 - 0.002 * (k - 100.0),
                }
            })
            .collect();
        let mut skew_spec = flat_spec.clone();
        skew_spec.quotes = skew_quotes.clone();
        let rep_call = price_gamma_swap_and_variance_option(&skew_spec);
        assert!(rep_call.variance_gamma_skew_spread > 0.0);

        let mut put_spec = skew_spec.clone();
        put_spec.is_variance_call = false;
        let rep_put = price_gamma_swap_and_variance_option(&put_spec);
        assert!(
            (rep_call.variance_option_price
                - rep_put.variance_option_price
                - (rep_call.expected_realized_variance - 0.04))
                .abs()
                < 1e-12
        );

        let h = 1e-5;
        let eval_v0 = |v0: f64| {
            let mut s = skew_spec.clone();
            s.heston.v0 = v0;
            price_gamma_swap_and_variance_option(&s).variance_option_price
        };
        let fd = (eval_v0(0.04 + h) - eval_v0(0.04 - h)) / (2.0 * h);
        type D = Dual<1>;
        let q_d: Vec<VarianceSwapStrikeQuote<D>> = skew_quotes
            .iter()
            .map(|q| VarianceSwapStrikeQuote {
                strike: D::from_value(q.strike),
                implied_vol: D::from_value(q.implied_vol),
            })
            .collect();
        let spec_d = GammaSwapAndVarOptionSpec {
            spot: D::from_value(100.0),
            rate: D::zero(),
            div_yield: D::zero(),
            expiry: D::one(),
            cut_strike: D::from_value(100.0),
            variance_strike: D::from_value(0.04),
            quotes: q_d,
            heston: HestonParams {
                v0: D::variable(0.04, 0),
                kappa: D::from_value(2.0),
                theta: D::from_value(0.04),
                sigma_v: D::from_value(0.35),
                rho: D::from_value(-0.7),
            },
            is_variance_call: true,
        };
        let ad = price_gamma_swap_and_variance_option(&spec_d).variance_option_price;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}
