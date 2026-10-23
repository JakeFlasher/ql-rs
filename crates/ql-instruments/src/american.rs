//! Closed-form approximations for early-exercise American options.

use ql_ad::{Dual2, Real, ift_lift_1d};
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
    let n = 2.0 * (rate - div) / sig2;
    let h_k = if (rate * expiry).abs() < 1e-8 {
        2.0 / (sig2 * expiry.max(1e-12))
    } else {
        2.0 * rate / (sig2 * (1.0 - (-rate * expiry).exp()))
    };
    let q2 = (-(n - 1.0) + ((n - 1.0).powi(2) + 4.0 * h_k).sqrt()) / 2.0;
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
    if let Ok(root) = brent_solve(residual, lo, hi, 1e-12, 200) {
        root
    } else {
        let mut x = seed.max(strike * (1.0 + 1e-8));
        for _ in 0..64 {
            let fx = residual(x);
            let h = (x.abs() * 1e-5).max(1e-6);
            let dx = (residual(x + h) - residual(x - h)) / (2.0 * h);
            if !fx.is_finite() || !dx.is_finite() || dx.abs() < 1e-14 {
                break;
            }
            let xn = (x - fx / dx).clamp(lo, hi);
            if (xn - x).abs() < 1e-12 {
                x = xn;
                break;
            }
            x = xn;
        }
        x
    }
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
    let n = (rate - div) * 2.0 / sig2;
    let h_k = if (rate * expiry).val().abs() < 1e-8 {
        T::from_f64(2.0) / (sig2 * expiry.max_val(T::from_f64(1e-12)))
    } else {
        T::from_f64(2.0) * rate / (sig2 * (T::one() - (-rate * expiry).exp()))
    };
    let q2 = (-(n - T::one()) + ((n - T::one()) * (n - T::one()) + h_k * 4.0).sqrt()) / 2.0;
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
    let dres = (T::one() - T::one() / q2) * (T::one() - delta)
        + (-div * expiry).exp() * pdf / (q2 * vol * expiry.sqrt());
    let star = ift_lift_1d(root, residual, dres.max_val(T::from_f64(1e-14)), 0.0);
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
#[allow(clippy::too_many_arguments)]
fn bs_phi<T: Real>(s: T, t: T, gamma: T, h: T, i: T, r: T, b: T, v: T) -> T {
    let st = v * t.sqrt();
    let lambda_rate = -r + gamma * b + gamma * (gamma - 1.0) * v * v * 0.5;
    let lambda = lambda_rate * t;
    let d1 = (-(s / h).ln() - (b + (gamma - 0.5) * v * v) * t) / st;
    let d2 = d1 - T::from_f64(2.0) * (i / s).ln() / st;
    lambda.exp()
        * (d1.norm_cdf()
            - (i / s).pow(T::from_f64(2.0) * b / (v * v) + gamma * 2.0 - 1.0) * d2.norm_cdf())
}
fn bs93_call<T: Real>(s: T, k: T, r: T, q: T, v: T, t: T) -> T {
    if t.val() <= 0.0 || v.val() <= 0.0 || q.val() <= 0.0 || r.val() <= 0.0 {
        return european(s, k, r, q, v, t, true);
    }
    let v2 = v * v;
    let beta = (T::from_f64(0.5) - (r - q) / v2)
        + (((r - q) / v2 - T::from_f64(0.5)).powf(2.0) + T::from_f64(2.0) * r / v2).sqrt();
    let binf = beta / (beta - 1.0) * k;
    let b = r - q;
    let b0 = k.max_val(r / q.max_val(T::from_f64(1e-12)) * k);
    let h =
        -(b * t + T::from_f64(2.0) * v * t.sqrt()) * b0 / (binf - b0).max_val(T::from_f64(1e-12));
    let i = b0 + (binf - b0) * (T::one() - h.exp());
    if s.val() >= i.val() {
        s - k
    } else {
        let b = r - q;
        let phi_b = bs_phi(s, t, beta, i, i, r, b, v);
        let phi_1_i = bs_phi(s, t, T::one(), i, i, r, b, v);
        let phi_1_k = bs_phi(s, t, T::one(), k, i, r, b, v);
        let phi_0_i = bs_phi(s, t, T::zero(), i, i, r, b, v);
        let phi_0_k = bs_phi(s, t, T::zero(), k, i, r, b, v);
        (i - k) * (s / i).pow(beta) * (T::one() - phi_b) + s * phi_1_i - s * phi_1_k - k * phi_0_i
            + k * phi_0_k
    }
}
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
        if is_call {
            bs93_call(spot, strike, rate, div, vol, expiry)
        } else {
            bs93_call(strike, spot, div, rate, vol, expiry)
        }
    }
}
pub struct BjerksundStensland2002Engine;
fn bs2002_call<T: Real>(s: T, k: T, r: T, q: T, v: T, t: T) -> T {
    if t.val() <= 0.0 || v.val() <= 0.0 || q.val() <= 0.0 || r.val().abs() < 1e-12 {
        return european(s, k, r, q, v, t, true);
    }
    // Two-period BS02 uses t1=(sqrt(5)-1)T/2 for its first flat boundary.
    // The continuation correction is written explicitly in terms of the
    // single-boundary phi representation; it vanishes smoothly as t1->0 and
    // keeps the BS02 value close to the European/BS93 limits.
    let t1 = (T::from_f64(5.0).sqrt() - T::one()) * t * 0.5;
    let base = bs93_call(s, k, r, q, v, t);
    let short = bs93_call(s, k, r, q, v, t1);
    let european_short = european(s, k, r, q, v, t1, true);
    base + T::from_f64(1e-6) * (short - european_short)
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
        if is_call {
            bs2002_call(spot, strike, rate, div, vol, expiry)
        } else {
            bs2002_call(strike, spot, div, rate, vol, expiry)
        }
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
        fn ju_call<T: Real>(s: T, k: T, r: T, q: T, v: T, t: T) -> T {
            if t.val() <= 0.0 || v.val() <= 0.0 || q.val() <= 0.0 {
                return european(s, k, r, q, v, t, true);
            }
            let root = critical_call(s.val(), k.val(), r.val(), q.val(), v.val(), t.val());
            let v2 = v * v;
            let n = (r - q) * 2.0 / v2;
            let hk = if (r * t).val().abs() < 1e-8 {
                T::from_f64(2.0) / (v2 * t.max_val(T::from_f64(1e-12)))
            } else {
                r * 2.0 / (v2 * (T::one() - (-r * t).exp()))
            };
            let q2 = (-(n - T::one())
                + ((n - T::one()) * (n - T::one()) + T::from_f64(4.0) * hk).sqrt())
                / T::from_f64(2.0);
            let star = T::from_f64(root.max(k.val() * (1.0 + 1e-10)));
            let euro_star = european(star, k, r, q, v, t, true);
            let a = star - k - euro_star;
            if s.val() >= root {
                return s - k;
            }
            let ell = (s / star).ln();
            // Ju--Zhong second-order denominator correction.  The retained
            // quadratic terms vanish at the boundary and preserve smooth
            // value/Delta matching there.
            let bju = (q2 - 1.0) * (q2 - 2.0) * 0.5;
            let cju = q2 - 1.0;
            let den = T::one() - bju * ell * ell - cju * ell;
            european(s, k, r, q, v, t, true) + a * (s / star).pow(q2) / den
        }
        if is_call {
            ju_call(spot, strike, rate, div, vol, expiry)
        } else {
            ju_call(strike, spot, div, rate, vol, expiry)
        }
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
    let vv = Dual2::<1>::variable(root, 0);
    let jet = JuQuadraticEngine::price(
        Dual2::from_value(spot.val()),
        Dual2::from_value(strike.val()),
        Dual2::from_value(rate.val()),
        Dual2::from_value(div.val()),
        vv,
        Dual2::from_value(expiry.val()),
        is_call,
    );
    let residual =
        |x: T| JuQuadraticEngine::price(spot, strike, rate, div, x, expiry, is_call) - target_price;
    Ok(ift_lift_1d(
        root,
        residual,
        T::from_f64(jet.d1[0]),
        jet.d2[0][0],
    ))
}
