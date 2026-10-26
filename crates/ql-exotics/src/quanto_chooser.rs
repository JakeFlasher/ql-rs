//! Simple chooser and quanto Black--Scholes analytics.

use ql_ad::Real;
use ql_ad::ift_lift_1d;
use ql_instruments::OptionType;
use ql_math::{bivariate_normal_cdf, brent_solve};
use ql_vol::black_scholes_price;

#[derive(Copy, Clone, Debug)]
pub struct ComplexChooserSpec<T: Real> {
    pub spot: T,
    pub strike_call: T,
    pub strike_put: T,
    pub choosing_time: T,
    pub call_expiry: T,
    pub put_expiry: T,
    pub rate: T,
    pub div_yield: T,
    pub vol: T,
}

/// Rubinstein simple chooser with decision date `chooser_time`.  The holder
/// chooses a call or put with common maturity `expiry`; the decomposition is
/// a vanilla call plus a put with adjusted strike and expiry at the decision
/// date.
#[allow(clippy::too_many_arguments)]
pub fn price_simple_chooser<T: Real>(
    spot: T,
    strike: T,
    rate: T,
    dividend: T,
    vol: T,
    expiry: T,
    chooser_time: T,
) -> T {
    let tc = chooser_time
        .max_val(T::zero())
        .min_val(expiry.max_val(T::zero()));
    let rem = expiry - tc;
    let adjusted_strike = strike * (-(rate - dividend) * rem).exp();
    black_scholes_price(spot, strike, vol, expiry, rate, dividend, true)
        + (-(dividend) * rem).exp()
            * black_scholes_price(spot, adjusted_strike, vol, tc, rate, dividend, false)
}

/// Quanto European option.  Under the domestic measure the foreign asset's
/// effective dividend yield is `r_for + rho*sigma_s*sigma_fx`.
#[allow(clippy::too_many_arguments)]
pub fn price_quanto_option<T: Real>(
    spot: T,
    strike: T,
    vol_asset: T,
    expiry: T,
    domestic_rate: T,
    foreign_rate: T,
    fx_vol: T,
    correlation: T,
    is_call: bool,
) -> T {
    let effective_div = foreign_rate + correlation * vol_asset * fx_vol;
    black_scholes_price(
        spot,
        strike,
        vol_asset,
        expiry,
        domestic_rate,
        effective_div,
        is_call,
    )
}

pub fn price_complex_chooser<T: Real>(s: &ComplexChooserSpec<T>) -> T {
    let tc = s.choosing_time;
    let call_tau = (s.call_expiry - tc).max_val(T::from_f64(1e-12));
    let put_tau = (s.put_expiry - tc).max_val(T::from_f64(1e-12));
    let residual_f64 = |x: f64| {
        black_scholes_price(
            x,
            s.strike_call.val(),
            s.vol.val(),
            call_tau.val(),
            s.rate.val(),
            s.div_yield.val(),
            true,
        )
        .val()
            - black_scholes_price(
                x,
                s.strike_put.val(),
                s.vol.val(),
                put_tau.val(),
                s.rate.val(),
                s.div_yield.val(),
                false,
            )
            .val()
    };
    let residual = |x: T| {
        black_scholes_price(x, s.strike_call, s.vol, call_tau, s.rate, s.div_yield, true)
            - black_scholes_price(x, s.strike_put, s.vol, put_tau, s.rate, s.div_yield, false)
    };
    let root = match brent_solve(
        residual_f64,
        1e-8,
        s.spot.val() * 20.0 + s.strike_call.val() * 10.0 + s.strike_put.val() * 10.0 + 1.0,
        1e-12,
        200,
    ) {
        Ok(x) => x,
        Err(_) => {
            // Safeguarded Newton fallback using the same exact analytical
            // delta difference used by the AD lift below.
            let mut x = s.spot.val().max(1e-8);
            for _ in 0..64 {
                let xt = T::from_f64(x);
                let dc = ((xt / s.strike_call).ln()
                    + (s.rate - s.div_yield + s.vol * s.vol * 0.5) * call_tau)
                    / (s.vol * call_tau.sqrt());
                let dp = ((xt / s.strike_put).ln()
                    + (s.rate - s.div_yield + s.vol * s.vol * 0.5) * put_tau)
                    / (s.vol * put_tau.sqrt());
                let d = (-s.div_yield * call_tau).exp() * dc.norm_cdf()
                    + (-s.div_yield * put_tau).exp() * (-dp).norm_cdf();
                if d.val().abs() <= 1e-14 {
                    break;
                }
                let xn = (x - residual_f64(x) / d.val()).max(1e-10);
                if (xn - x).abs() < 1e-12 {
                    x = xn;
                    break;
                }
                x = xn;
            }
            x
        }
    };
    let root_t = T::from_f64(root);
    let d1c = ((root_t / s.strike_call).ln()
        + (s.rate - s.div_yield + s.vol * s.vol * 0.5) * call_tau)
        / (s.vol * call_tau.sqrt());
    let d1p = ((root_t / s.strike_put).ln()
        + (s.rate - s.div_yield + s.vol * s.vol * 0.5) * put_tau)
        / (s.vol * put_tau.sqrt());
    let fx = (-s.div_yield * call_tau).exp() * d1c.norm_cdf()
        + (-s.div_yield * put_tau).exp() * (-d1p).norm_cdf();
    let fxx = (-s.div_yield.val() * call_tau.val()).exp() * d1c.val().norm_pdf()
        / (root * s.vol.val() * call_tau.val().sqrt())
        - (-s.div_yield.val() * put_tau.val()).exp() * d1p.val().norm_pdf()
            / (root * s.vol.val() * put_tau.val().sqrt());
    let i_star = if fx.val().abs() > 1e-14 {
        ift_lift_1d(root, residual, fx, fxx)
    } else {
        T::from_f64(root)
    };
    let st = s.vol * tc.sqrt();
    let d1 = ((s.spot / i_star).ln() + (s.rate - s.div_yield + s.vol * s.vol * 0.5) * tc) / st;
    let d2 = d1 - st;
    let yc = ((s.spot / s.strike_call).ln()
        + (s.rate - s.div_yield + s.vol * s.vol * 0.5) * s.call_expiry)
        / (s.vol * s.call_expiry.sqrt());
    let yp = ((s.spot / s.strike_put).ln()
        + (s.rate - s.div_yield + s.vol * s.vol * 0.5) * s.put_expiry)
        / (s.vol * s.put_expiry.sqrt());
    let r1 = (tc / s.call_expiry).sqrt();
    let r2 = (tc / s.put_expiry).sqrt();
    s.spot * (-s.div_yield * s.call_expiry).exp() * bivariate_normal_cdf(d1, yc, r1)
        - s.strike_call
            * (-s.rate * s.call_expiry).exp()
            * bivariate_normal_cdf(d2, yc - s.vol * s.call_expiry.sqrt(), r1)
        - s.spot * (-s.div_yield * s.put_expiry).exp() * bivariate_normal_cdf(-d1, -yp, r2)
        + s.strike_put
            * (-s.rate * s.put_expiry).exp()
            * bivariate_normal_cdf(-d2, -yp + s.vol * s.put_expiry.sqrt(), r2)
}

/// Convenience alias accepting an instrument option type.
#[allow(clippy::too_many_arguments)]
pub fn price_quanto_option_type<T: Real>(
    spot: T,
    strike: T,
    vol_asset: T,
    expiry: T,
    domestic_rate: T,
    foreign_rate: T,
    fx_vol: T,
    correlation: T,
    option_type: OptionType,
) -> T {
    price_quanto_option(
        spot,
        strike,
        vol_asset,
        expiry,
        domestic_rate,
        foreign_rate,
        fx_vol,
        correlation,
        option_type.is_call(),
    )
}

mod evaluate_rubinstein_simple_chooser_option_module {
    // Bouzoubaa & Osseiran (2010) & Haug (2007) Rubinstein (1991) Chooser Option
    // Put–Call Parity Decomposition & Decision-Flexibility Premium Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct RubinsteinChooserOptionSpec<T: Real = f64> {
        pub spot_price: T,
        pub strike_price: T,
        pub risk_free_rate: T,
        pub dividend_yield: T,
        pub volatility: T,
        pub choice_date_tc: T,
        pub final_maturity_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct RubinsteinChooserOptionReport<T: Real = f64> {
        pub underlying_call_t: T,
        pub underlying_put_t: T,
        pub choice_date_effective_put_tc: T,
        pub rubinstein_chooser_price: T,
        pub straddle_upper_bound: T,
        pub decision_flexibility_premium: T,
    }

    pub fn evaluate_rubinstein_simple_chooser_option<T: Real>(
        spec: &RubinsteinChooserOptionSpec<T>,
    ) -> RubinsteinChooserOptionReport<T> {
        let s = spec.spot_price.max_val(T::from_f64(1e-4));
        let k = spec.strike_price.max_val(T::from_f64(1e-4));
        let r = spec.risk_free_rate;
        let q = spec.dividend_yield;
        let v = spec.volatility.max_val(T::from_f64(1e-4));
        let tc = spec.choice_date_tc.max_val(T::from_f64(0.01));
        let t = spec.final_maturity_t.max_val(tc + T::from_f64(0.01));

        // Standard Black-Scholes call and put for maturity T:
        let std_t = v * t.sqrt();
        let d1 = ((s / k).ln() + (r - q + T::from_f64(0.5) * v * v) * t) / std_t;
        let d2 = d1 - std_t;
        let df_q_t = (-q * t).exp();
        let df_r_t = (-r * t).exp();
        let call_t = s * df_q_t * norm_cdf(d1) - k * df_r_t * norm_cdf(d2);
        let put_t = k * df_r_t * norm_cdf(-d2) - s * df_q_t * norm_cdf(-d1);

        // By Rubinstein (1991) put-call parity at choice date T_c:
        // max(C(S_Tc, K, T-Tc), P(S_Tc, K, T-Tc)) = C(S_Tc, K, T-Tc) + e^{-q(T-Tc)} max(0, K e^{-(r-q)(T-Tc)} - S_Tc)
        let std_tc = v * tc.sqrt();
        let y1 = ((s / k).ln() + (r - q) * t + T::from_f64(0.5) * v * v * tc) / std_tc;
        let y2 = y1 - std_tc;
        let put_tc = k * df_r_t * norm_cdf(-y2) - s * df_q_t * norm_cdf(-y1);

        let chooser = call_t + put_tc;
        let straddle = call_t + put_t;
        let max_cp = if call_t > put_t { call_t } else { put_t };
        let flex_prem = chooser - max_cp;

        RubinsteinChooserOptionReport {
            underlying_call_t: call_t,
            underlying_put_t: put_t,
            choice_date_effective_put_tc: put_tc,
            rubinstein_chooser_price: chooser,
            straddle_upper_bound: straddle,
            decision_flexibility_premium: flex_prem,
        }
    }

    #[cfg(test)]
    mod domain_suite_51_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn rubinstein_chooser_option_invariants_and_ad_parity() {
            let spec = RubinsteinChooserOptionSpec {
                spot_price: 100.0,
                strike_price: 100.0,
                risk_free_rate: 0.05,
                dividend_yield: 0.02,
                volatility: 0.24,
                choice_date_tc: 0.4,
                final_maturity_t: 1.0,
            };
            let rep = evaluate_rubinstein_simple_chooser_option(&spec);
            // No-arbitrage Sandwich: max(C, P) < Chooser < C + P (Straddle)
            assert!(rep.decision_flexibility_premium > 0.0);
            assert!(rep.rubinstein_chooser_price < rep.straddle_upper_bound);

            let h = 1e-6;
            let r_ad = evaluate_rubinstein_simple_chooser_option(&RubinsteinChooserOptionSpec {
                spot_price: Dual::<1>::variable(spec.spot_price, 0),
                strike_price: Dual::<1>::from_value(spec.strike_price),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                dividend_yield: Dual::<1>::from_value(spec.dividend_yield),
                volatility: Dual::<1>::from_value(spec.volatility),
                choice_date_tc: Dual::<1>::from_value(spec.choice_date_tc),
                final_maturity_t: Dual::<1>::from_value(spec.final_maturity_t),
            });
            let mut sp = spec;
            sp.spot_price += h;
            let mut sm = spec;
            sm.spot_price -= h;
            let fd = (evaluate_rubinstein_simple_chooser_option(&sp).rubinstein_chooser_price
                - evaluate_rubinstein_simple_chooser_option(&sm).rubinstein_chooser_price)
                / (2.0 * h);
            assert!((r_ad.rubinstein_chooser_price.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_rubinstein_simple_chooser_option_module::*;
