//! Continuous and discretely monitored lookback options.
//! Formulas follow the Goldman--Sosin--Gatto and Conze--Viswanathan
//! implementations used by QuantLib's analytic lookback engines.

use ql_ad::Real;
use ql_instruments::OptionType;

/// Floating-strike continuous lookback. `extremum` is the running minimum
/// for a call and the running maximum for a put.
#[allow(clippy::too_many_arguments)]
pub fn price_floating_lookback_continuous<T: Real>(
    spot: T,
    extremum: T,
    rate: T,
    dividend: T,
    vol: T,
    expiry: T,
    option_type: OptionType,
) -> T {
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        return if option_type.is_call() {
            (spot - extremum).max_val(T::zero())
        } else {
            (extremum - spot).max_val(T::zero())
        };
    }
    let eta = if option_type.is_call() { 1.0 } else { -1.0 };
    let std = vol * expiry.sqrt();
    let lambda = T::from_f64(2.0) * (rate - dividend) / (vol * vol);
    let s = spot / extremum;
    let d1 = s.ln() / std + T::from_f64(0.5) * (lambda + 1.0) * std;
    let et = T::from_f64(eta);
    let n1 = (et * d1).norm_cdf();
    let n2 = (et * (d1 - std)).norm_cdf();
    let n3 = (et * (-d1 + lambda * std)).norm_cdf();
    let n4 = (-et * d1).norm_cdf();
    let dr = (-rate * expiry).exp();
    let dq = (-dividend * expiry).exp();
    let ratio = if lambda.val().abs() < 1e-7 {
        // Expand the 0/0 quotient through first order in lambda.  Using only
        // the value limit would drop the carry (rate/dividend) sensitivity
        // when lambda is represented by a Dual at exactly zero.
        let d0 = d1 - T::from_f64(0.5) * std * lambda;
        let n0 = d0.norm_pdf();
        let c0 = (-et * d0).norm_cdf();
        let lns = s.ln();
        let aa = T::from_f64(0.5) * std * std;
        let r0 = et * std * n0 - (lns + aa) * c0;
        let r1 = T::from_f64(0.5) * (lns * lns - aa * aa) * c0
            + et * std * T::from_f64(0.5) * n0 * (aa - lns);
        r0 + lambda * r1
    } else {
        (s.pow(-lambda) * n3 - dq / dr * n4) / lambda
    };
    T::from_f64(eta) * (spot * dq * n1 - extremum * dr * n2 + spot * dr * ratio)
}

/// Fixed-strike continuous lookback. `extremum` is the running maximum for a
/// call and running minimum for a put (the convention used by QuantLib).
#[allow(clippy::too_many_arguments)]
pub fn price_fixed_lookback_continuous<T: Real>(
    spot: T,
    extremum: T,
    strike: T,
    rate: T,
    dividend: T,
    vol: T,
    expiry: T,
    option_type: OptionType,
) -> T {
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        return if option_type.is_call() {
            (extremum - strike).max_val(T::zero())
        } else {
            (strike - extremum).max_val(T::zero())
        };
    }
    let eta = if option_type.is_call() { 1.0 } else { -1.0 };
    let std = vol * expiry.sqrt();
    let lambda = T::from_f64(2.0) * (rate - dividend) / (vol * vol);
    let a_term = |m: T| {
        let ss = spot / m;
        let d1 = ss.ln() / std + T::from_f64(0.5) * (lambda + 1.0) * std;
        let et = T::from_f64(eta);
        let n1 = (et * d1).norm_cdf();
        let n2 = (et * (d1 - std)).norm_cdf();
        let n3 = (et * (d1 - lambda * std)).norm_cdf();
        let n4 = (et * d1).norm_cdf();
        let ratio = if lambda.val().abs() < 1e-7 {
            // First-order expansion of the quotient, retaining the exact
            // lambda=0 value and its AD carry sensitivity.
            let d0 = d1 - T::from_f64(0.5) * std * lambda;
            let n0 = d0.norm_pdf();
            let c0 = (et * d0).norm_cdf();
            let lns = ss.ln();
            let aa = T::from_f64(0.5) * std * std;
            let r0 = -(lns + aa) * c0 - et * std * n0;
            let r1 = T::from_f64(0.5) * (lns * lns - aa * aa) * c0
                + et * std * T::from_f64(0.5) * n0 * (lns - aa);
            r0 + lambda * r1
        } else {
            (ss.pow(-lambda) * n3 - (-dividend * expiry).exp() / (-rate * expiry).exp() * n4)
                / lambda
        };
        T::from_f64(eta)
            * (spot * (-dividend * expiry).exp() * n1
                - m * (-rate * expiry).exp() * n2
                - spot * (-rate * expiry).exp() * ratio)
    };
    let b_term = a_term(strike);
    if option_type.is_call() {
        if strike.val() <= extremum.val() {
            a_term(extremum) + T::from_f64(eta) * (-rate * expiry).exp() * (extremum - strike)
        } else {
            b_term
        }
    } else if strike.val() >= extremum.val() {
        a_term(extremum) + T::from_f64(eta) * (-rate * expiry).exp() * (extremum - strike)
    } else {
        b_term
    }
}

/// BGK continuity-corrected discrete floating lookback.  `observations` is
/// the number of monitoring dates over the option life.
#[allow(clippy::too_many_arguments)]
pub fn price_floating_lookback_discrete<T: Real>(
    spot: T,
    extremum: T,
    rate: T,
    dividend: T,
    vol: T,
    expiry: T,
    option_type: OptionType,
    observations: usize,
) -> T {
    let dt = expiry / observations.max(1) as f64;
    let beta = 0.5825971579390107;
    let adjusted = if option_type.is_call() {
        // A discretely sampled minimum is above the continuous minimum.
        extremum * (vol * dt.sqrt() * beta).exp()
    } else {
        // A discretely sampled maximum is below the continuous maximum.
        extremum * (-(vol * dt.sqrt() * beta)).exp()
    };
    price_floating_lookback_continuous(spot, adjusted, rate, dividend, vol, expiry, option_type)
}

/// BGK continuity-corrected discrete fixed lookback.
#[allow(clippy::too_many_arguments)]
pub fn price_fixed_lookback_discrete<T: Real>(
    spot: T,
    extremum: T,
    strike: T,
    rate: T,
    dividend: T,
    vol: T,
    expiry: T,
    option_type: OptionType,
    observations: usize,
) -> T {
    let dt = expiry / observations.max(1) as f64;
    let beta = 0.5825971579390107;
    let adjusted = if option_type.is_call() {
        // A discretely sampled maximum is below the continuous maximum.
        extremum * (-(vol * dt.sqrt() * beta)).exp()
    } else {
        // A discretely sampled minimum is above the continuous minimum.
        extremum * (vol * dt.sqrt() * beta).exp()
    };
    price_fixed_lookback_continuous(
        spot,
        adjusted,
        strike,
        rate,
        dividend,
        vol,
        expiry,
        option_type,
    )
}
