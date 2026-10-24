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
    if kappa.val().abs() < 1e-12 {
        v0
    } else {
        theta + (v0 - theta) * (T::one() - (-kappa * expiry).exp()) / (kappa * expiry)
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
    let f = spot * ((rate - div_yield) * expiry).exp();
    let lo = lower_barrier.max_val(T::zero());
    let hi = upper_barrier;
    let fc = f.max_val(lo).min_val(hi);
    let kc = cut_strike.max_val(lo).min_val(hi);
    let mut sum = T::zero();
    for w in quotes.windows(2) {
        let k0 = w[0].strike;
        let k1 = w[1].strike;
        let km = (k0 + k1) * 0.5;
        if km.val() < lo.val() || km.val() > hi.val() {
            continue;
        }
        let v = (w[0].implied_vol + w[1].implied_vol) * 0.5;
        let c = black_scholes_price(spot, km, v, expiry, rate, div_yield, km.val() >= kc.val());
        sum += (k1 - k0) * c / (km * km);
    }
    let drift = fc * (T::one() / fc - T::one() / kc) - (kc / fc).ln();
    let variance = (drift + (rate * expiry).exp() * sum) * (T::from_f64(2.0) / expiry);
    (variance, variance.max_val(T::zero()).sqrt())
}
