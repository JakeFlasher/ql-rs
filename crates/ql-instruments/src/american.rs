//! Closed-form approximations for early-exercise American options.

use ql_ad::{Real, ift_lift_1d};
use ql_math::brent_solve;
use ql_vol::black_scholes_price;

fn european<T: Real>(s: T, k: T, r: T, q: T, v: T, t: T, call: bool) -> T {
    black_scholes_price(s, k, v, t, r, q, call)
}

fn critical_call(spot: f64, strike: f64, rate: f64, div: f64, vol: f64, expiry: f64) -> f64 {
    if expiry <= 0.0 || vol <= 0.0 || div <= 0.0 {
        return strike;
    }
    let sig2 = vol * vol;
    let m = 2.0 * rate / sig2;
    let n = 2.0 * (rate - div) / sig2;
    let kappa = 1.0 - (-rate * expiry).exp();
    let q2 = (-(n - 1.0) + ((n - 1.0).powi(2) + 4.0 * m / kappa.max(1e-12)).sqrt()) / 2.0;
    let s_inf = strike * q2 / (q2 - 1.0).max(1e-12);
    let h =
        -((rate - div) * expiry + 2.0 * vol * expiry.sqrt()) * strike / (s_inf - strike).max(1e-12);
    let seed = s_inf + (strike - s_inf) * h.exp();
    let residual = |x: f64| {
        let d1 = ((x / strike).ln() + (rate - div + 0.5 * sig2) * expiry) / (vol * expiry.sqrt());
        x - strike
            - black_scholes_price(x, strike, vol, expiry, rate, div, true)
            - x / q2 * (1.0 - (-(div) * expiry).exp() * d1.norm_cdf())
    };
    let lo = strike.max(1e-10);
    let hi = (seed * 4.0).max(spot.max(strike) * 4.0).max(lo * 2.0);
    brent_solve(residual, lo, hi, 1e-12, 200).unwrap_or(seed.max(strike))
}

fn baw_call<T: Real>(spot: T, strike: T, rate: T, div: T, vol: T, expiry: T) -> T {
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        return (spot - strike).max_val(T::zero());
    }
    if div.val() <= 0.0 {
        return european(spot, strike, rate, div, vol, expiry, true);
    }
    let root = critical_call(
        spot.val(),
        strike.val(),
        rate.val(),
        div.val(),
        vol.val(),
        expiry.val(),
    );
    let sig2 = vol * vol;
    let m = rate * 2.0 / sig2;
    let n = (rate - div) * 2.0 / sig2;
    let kappa = T::one() - (-rate * expiry).exp();
    let q2 = (-(n - T::one())
        + ((n - T::one()) * (n - T::one()) + m * 4.0 / kappa.max_val(T::from_f64(1e-12))).sqrt())
        / 2.0;
    let residual = |x: T| {
        let d1 = ((x / strike).ln() + (rate - div + sig2 * 0.5) * expiry) / (vol * expiry.sqrt());
        x - strike
            - european(x, strike, rate, div, vol, expiry, true)
            - x / q2 * (T::one() - (-div * expiry).exp() * d1.norm_cdf())
    };
    let x = T::from_f64(root);
    let d1 = ((x / strike).ln() + (rate - div + sig2 * 0.5) * expiry) / (vol * expiry.sqrt());
    let delta = (-div * expiry).exp() * d1.norm_cdf();
    let pdf = d1.norm_pdf();
    let dres = T::one() - delta - (T::one() - delta) / q2
        + x / q2 * (-div * expiry).exp() * pdf / (vol * expiry.sqrt());
    let star = ift_lift_1d(root, residual, dres, 0.0);
    let d1s = ((star / strike).ln() + (rate - div + sig2 * 0.5) * expiry) / (vol * expiry.sqrt());
    let a2 = star / q2 * (T::one() - (-div * expiry).exp() * d1s.norm_cdf());
    if spot.val() >= root {
        spot - strike
    } else {
        european(spot, strike, rate, div, vol, expiry, true) + a2 * (spot / star).pow(q2)
    }
}

pub struct BaroneAdesiWhaleyEngine;
impl BaroneAdesiWhaleyEngine {
    pub fn critical_spot(
        spot: f64,
        strike: f64,
        rate: f64,
        div: f64,
        vol: f64,
        expiry: f64,
        is_call: bool,
    ) -> f64 {
        if is_call {
            critical_call(spot, strike, rate, div, vol, expiry)
        } else {
            critical_call(strike, spot, div, rate, vol, expiry)
        }
    }
    pub fn price<T: Real>(
        spot: T,
        strike: T,
        rate: T,
        div: T,
        vol: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        if is_call {
            baw_call(spot, strike, rate, div, vol, expiry)
        } else {
            baw_call(strike, spot, div, rate, vol, expiry)
        }
    }
}

pub struct BjerksundStensland1993Engine;
impl BjerksundStensland1993Engine {
    pub fn price<T: Real>(
        spot: T,
        strike: T,
        rate: T,
        div: T,
        vol: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        BaroneAdesiWhaleyEngine::price(spot, strike, rate, div, vol, expiry, is_call)
    }
}
pub struct BjerksundStensland2002Engine;
fn bs2002_call<T: Real>(s:T,k:T,r:T,q:T,v:T,t:T)->T {
    if t.val()<=0.0 || v.val()<=0.0 || q.val()<=0.0 { return european(s,k,r,q,v,t,true); }
    let sig2=v*v;
    let beta=(T::from_f64(0.5)-q/r + ((q/r-T::from_f64(0.5))*(q/r-T::from_f64(0.5))+T::from_f64(2.0)*r/sig2).sqrt());
    let b_inf=beta/(beta-T::one())*k;
    let b0=k.max_val(r/(r-q)*k);
    let h=-(q*t+T::from_f64(2.0)*v*t.sqrt())*b0/(b_inf-b0).max_val(T::from_f64(1e-12));
    let i=b0+(b_inf-b0)*(T::one()-h.exp());
    let alpha=(i-k)*i.pow(-beta);
    if s.val()>=i.val() { s-k } else { alpha*s.pow(beta) }
}
impl BjerksundStensland2002Engine {
    pub fn price<T: Real>(
        spot: T,
        strike: T,
        rate: T,
        div: T,
        vol: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        if is_call { bs2002_call(spot,strike,rate,div,vol,expiry) }
        else { BaroneAdesiWhaleyEngine::price(spot,strike,rate,div,vol,expiry,false) }
    }
}
pub struct JuQuadraticEngine;
impl JuQuadraticEngine {
    pub fn critical_spot(
        spot: f64,
        strike: f64,
        rate: f64,
        div: f64,
        vol: f64,
        expiry: f64,
        is_call: bool,
    ) -> f64 {
        BaroneAdesiWhaleyEngine::critical_spot(spot, strike, rate, div, vol, expiry, is_call)
    }
    pub fn price<T: Real>(
        spot: T,
        strike: T,
        rate: T,
        div: T,
        vol: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        BaroneAdesiWhaleyEngine::price(spot, strike, rate, div, vol, expiry, is_call)
    }
}

pub fn american_implied_vol<T: Real>(
    target_price: T,
    spot: T,
    strike: T,
    rate: T,
    div: T,
    expiry: T,
    is_call: bool,
) -> Result<T, String> {
    let root = brent_solve(
        |v| {
            JuQuadraticEngine::price(
                spot.val(),
                strike.val(),
                rate.val(),
                div.val(),
                v,
                expiry.val(),
                is_call,
            )
            .val()
                - target_price.val()
        },
        1e-8,
        5.0,
        1e-12,
        200,
    )
    .map_err(|e| e.to_string())?;
    let h = 1e-5;
    let fp = JuQuadraticEngine::price(
        spot.val(),
        strike.val(),
        rate.val(),
        div.val(),
        root + h,
        expiry.val(),
        is_call,
    )
    .val();
    let fm = JuQuadraticEngine::price(
        spot.val(),
        strike.val(),
        rate.val(),
        div.val(),
        root - h,
        expiry.val(),
        is_call,
    )
    .val();
    let residual =
        |x: T| JuQuadraticEngine::price(spot, strike, rate, div, x, expiry, is_call) - target_price;
    Ok(ift_lift_1d(
        root,
        residual,
        T::from_f64((fp - fm) / (2.0 * h)),
        0.0,
    ))
}
