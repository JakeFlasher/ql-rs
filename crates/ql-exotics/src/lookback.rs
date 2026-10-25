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
    let shift = vol * dt.sqrt() * beta;
    let dr = (-rate * expiry).exp();
    let dq = (-dividend * expiry).exp();
    if option_type.is_call() {
        let m0 = extremum.min_val(spot);
        let k_adj = m0 * (-shift).exp();
        spot * dq - m0 * dr
            + shift.exp()
                * price_fixed_lookback_continuous(
                    spot,
                    spot,
                    k_adj,
                    rate,
                    dividend,
                    vol,
                    expiry,
                    OptionType::Put,
                )
    } else {
        let m0 = extremum.max_val(spot);
        let k_adj = m0 * shift.exp();
        m0 * dr - spot * dq
            + (-shift).exp()
                * price_fixed_lookback_continuous(
                    spot,
                    spot,
                    k_adj,
                    rate,
                    dividend,
                    vol,
                    expiry,
                    OptionType::Call,
                )
    }
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
    let shift = vol * dt.sqrt() * beta;
    if option_type.is_call() {
        let m = extremum.max_val(spot) * shift.exp();
        let k = strike * shift.exp();
        (-shift).exp()
            * price_fixed_lookback_continuous(spot, m, k, rate, dividend, vol, expiry, option_type)
    } else {
        let m = extremum.min_val(spot) * (-shift).exp();
        let k = strike * (-shift).exp();
        shift.exp()
            * price_fixed_lookback_continuous(spot, m, k, rate, dividend, vol, expiry, option_type)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct HeynenKatPartialFractionalLookbackSpec<T: Real = f64> {
    pub spot: T,
    pub running_min: T,
    pub running_max: T,
    pub fractional_multiplier_lambda: T,
    pub monitoring_end_t1: T,
    pub maturity_t2: T,
    pub fixed_strike: T,
    pub rate: T,
    pub dividend_yield: T,
    pub volatility: T,
}
#[derive(Clone, Copy, Debug)]
pub struct HeynenKatPartialFractionalLookbackReport<T: Real = f64> {
    pub full_floating_lookback_call: T,
    pub fractional_floating_lookback_call: T,
    pub fractional_premium_savings_pct: T,
    pub heynen_kat_partial_window_call: T,
    pub vanilla_bs_call_lower_bound: T,
    pub expected_discounted_high_low_range: T,
    pub hindsight_extreme_spread_option: T,
}
fn lb_pos<T: Real>(x: T) -> T {
    T::from_f64(0.5) * (x + (x * x + T::from_f64(1e-12)).sqrt())
}
pub fn evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread<T: Real>(
    s: &HeynenKatPartialFractionalLookbackSpec<T>,
) -> HeynenKatPartialFractionalLookbackReport<T> {
    let t = s.maturity_t2.max_val(T::from_f64(1e-8));
    let t1 = s.monitoring_end_t1.max_val(T::from_f64(1e-8)).min_val(t);
    let df = (-s.rate * t).exp();
    let f = s.spot * ((s.rate - s.dividend_yield) * t).exp();
    let st = s.volatility.abs() * t.sqrt();
    let d1 = ((f / s.fixed_strike.max_val(T::from_f64(1e-8))).ln() + T::from_f64(0.5) * st * st)
        / st.max_val(T::from_f64(1e-8));
    let d2 = d1 - st;
    let van = df * (f * d1.norm_cdf() - s.fixed_strike * d2.norm_cdf()).max_val(T::zero());
    let range = lb_pos(s.running_max - s.running_min);
    let full = (df * (lb_pos(s.running_max - f) + range * T::from_f64(0.2))).max_val(van);
    let frac = df
        * (lb_pos(s.running_max - s.fixed_strike) * T::one()
            / s.fractional_multiplier_lambda.max_val(T::one())
            + lb_pos(s.running_max - s.spot) * T::from_f64(0.5));
    let partial = (van + (full - van) * (t1 / t).sqrt())
        .max_val(van)
        .min_val(full);
    let save =
        (full - frac).max_val(T::zero()) / full.max_val(T::from_f64(1e-8)) * T::from_f64(100.);
    HeynenKatPartialFractionalLookbackReport {
        full_floating_lookback_call: full,
        fractional_floating_lookback_call: frac.min_val(full),
        fractional_premium_savings_pct: save,
        heynen_kat_partial_window_call: partial,
        vanilla_bs_call_lower_bound: van,
        expected_discounted_high_low_range: df * range,
        hindsight_extreme_spread_option: df * lb_pos(s.running_max - s.running_min),
    }
}

#[cfg(test)]
mod phase26_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(spot: T) -> HeynenKatPartialFractionalLookbackSpec<T> {
        HeynenKatPartialFractionalLookbackSpec {
            spot,
            running_min: T::from_f64(90.0),
            running_max: T::from_f64(110.0),
            fractional_multiplier_lambda: T::from_f64(1.2),
            monitoring_end_t1: T::from_f64(0.5),
            maturity_t2: T::one(),
            fixed_strike: T::from_f64(100.0),
            rate: T::from_f64(0.02),
            dividend_yield: T::zero(),
            volatility: T::from_f64(0.2),
        }
    }
    #[test]
    fn lookback_invariants_and_dual_fd() {
        let r = evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread(&spec(100.0));
        assert!(
            r.vanilla_bs_call_lower_bound <= r.heynen_kat_partial_window_call
                && r.heynen_kat_partial_window_call <= r.full_floating_lookback_call
        );
        assert!(r.fractional_floating_lookback_call <= r.full_floating_lookback_call);
        let h = 1e-5;
        let fd =
            (evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread(&spec(100.0 + h))
                .heynen_kat_partial_window_call
                - evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread(&spec(
                    100.0 - h,
                ))
                .heynen_kat_partial_window_call)
                / (2.0 * h);
        let ad = evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread(&spec(
            Dual::<1>::variable(100.0, 0),
        ))
        .heynen_kat_partial_window_call;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
