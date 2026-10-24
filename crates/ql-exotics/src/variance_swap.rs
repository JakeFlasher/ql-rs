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
    // In the practically unbounded-corridor limit retain the exact primary
    // replication (including its endpoint Black tails) rather than introducing
    // a second quadrature discretisation.
    if lower_barrier.val() <= 1.0e-4 && upper_barrier.val() >= 1.0e4 {
        return replicate_variance_swap_fair_strike(
            spot, rate, div_yield, expiry, cut_strike, quotes,
        );
    }
    let forward = spot * ((rate - div_yield) * expiry).exp();
    let lo = lower_barrier.max_val(T::from_f64(1e-14));
    let hi = upper_barrier.max_val(lo + T::from_f64(1e-14));
    let ftilde = forward.max_val(lo).min_val(hi);
    let kctilde = cut_strike.max_val(lo).min_val(hi);
    let drift = ftilde * (T::one() / ftilde - T::one() / kctilde) - (kctilde / ftilde).ln();
    let mut integral = T::zero();
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
        if kc_val(cut_strike, a, b) {
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
            let is_call = xm.val() >= kctilde.val();
            integral += (x1 - x0)
                * (integrand(x0, v0, is_call)
                    + integrand(xm, vm, is_call) * 2.0
                    + integrand(x1, v1, is_call))
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
            integral += put_tail(end, quotes[0].implied_vol) - put_tail(lo, quotes[0].implied_vol);
        }
    }
    let kn = quotes[quotes.len() - 1].strike;
    if hi.val() > kn.val() {
        let begin = kn.max_val(lo);
        if hi.val() > begin.val() {
            integral += call_tail(begin, quotes[quotes.len() - 1].implied_vol)
                - call_tail(hi, quotes[quotes.len() - 1].implied_vol);
        }
    }
    let variance = (drift + (rate * expiry).exp() * integral) * (T::from_f64(2.0) / expiry);
    (variance, variance.max_val(T::zero()).sqrt())
}

fn kc_val<T: Real>(cut: T, a: T, b: T) -> bool {
    cut.val() > a.val() && cut.val() < b.val()
}
