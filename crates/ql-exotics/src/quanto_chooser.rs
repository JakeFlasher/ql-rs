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
