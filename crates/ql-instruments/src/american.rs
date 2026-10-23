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
    let gamma_e = (-div * expiry).exp() * pdf / (x * vol * expiry.sqrt());
    let a = (-div * expiry).exp() * pdf / (vol * expiry.sqrt());
    let d2res = -gamma_e - (-a / x + a * d1 / (x * vol * expiry.sqrt())) / q2;
    let star = ift_lift_1d(
        root,
        residual,
        dres.max_val(T::from_f64(1e-14)),
        d2res.val(),
    );
    let d1s = ((star / strike).ln() + (rate - div + sig2 * 0.5) * expiry) / (vol * expiry.sqrt());
    let a2 = star / q2 * (T::one() - (-div * expiry).exp() * d1s.norm_cdf());
    if spot.val() >= root {
        spot - strike
    } else {
        european(spot, strike, rate, div, vol, expiry, true) + a2 * (spot / star).pow(q2)
    }
}

fn critical_put(_spot: f64, strike: f64, rate: f64, div: f64, vol: f64, expiry: f64) -> f64 {
    if expiry <= 0.0 || vol <= 0.0 {
        return strike;
    }
    let v2 = vol * vol;
    let n = 2.0 * (rate - div) / v2;
    let hk = if (rate * expiry).abs() < 1e-8 {
        2.0 / (v2 * expiry.max(1e-12))
    } else {
        2.0 * rate / (v2 * (1.0 - (-rate * expiry).exp()))
    };
    let q1 = (-(n - 1.0) - ((n - 1.0).powi(2) + 4.0 * hk).sqrt()) / 2.0;
    let dres = |x: f64| {
        let d1 = ((x / strike).ln() + (rate - div + 0.5 * v2) * expiry) / (vol * expiry.sqrt());
        let pe = european(x, strike, rate, div, vol, expiry, false).val();
        strike - x - pe + x / q1 * (1.0 - (-div * expiry).exp() * (-d1).norm_cdf())
    };
    let hi = (strike * (1.0 - 1e-10)).max(1e-10);
    brent_solve(dres, 1e-10, hi, 1e-12, 200).unwrap_or_else(|_| {
        let mut x = (strike * 0.5).max(1e-10);
        for _ in 0..80 {
            let f = dres(x);
            let h = (x.abs() * 1e-5).max(1e-6);
            let df = (dres(x + h) - dres(x - h)) / (2.0 * h);
            if !df.is_finite() || df.abs() < 1e-14 {
                break;
            }
            x = (x - f / df).clamp(1e-10, hi);
        }
        x
    })
}

fn baw_put<T: Real>(spot: T, strike: T, rate: T, div: T, vol: T, expiry: T) -> T {
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        return (strike - spot).max_val(T::zero());
    }
    let v2 = vol * vol;
    let n = (rate - div) * 2.0 / v2;
    let hk = if (rate * expiry).val().abs() < 1e-8 {
        T::from_f64(2.0) / (v2 * expiry.max_val(T::from_f64(1e-12)))
    } else {
        rate * 2.0 / (v2 * (T::one() - (-rate * expiry).exp()))
    };
    let q1 = (-(n - T::one()) - ((n - T::one()) * (n - T::one()) + hk * 4.0).sqrt()) / 2.0;
    let root = critical_put(
        spot.val(),
        strike.val(),
        rate.val(),
        div.val(),
        vol.val(),
        expiry.val(),
    );
    let residual = |x: T| {
        let d1 = ((x / strike).ln() + (rate - div + v2 * 0.5) * expiry) / (vol * expiry.sqrt());
        strike - x - european(x, strike, rate, div, vol, expiry, false)
            + x / q1 * (T::one() - (-div * expiry).exp() * (-d1).norm_cdf())
    };
    let x = T::from_f64(root);
    let d1 = ((x / strike).ln() + (rate - div + v2 * 0.5) * expiry) / (vol * expiry.sqrt());
    let d = (-div * expiry).exp() * (-d1).norm_cdf();
    let pdf = d1.norm_pdf();
    let fx = -(-div * expiry).exp() * d1.norm_cdf()
        + T::one() / q1 * (T::one() - d)
        + (-div * expiry).exp() * pdf / (q1 * vol * expiry.sqrt());
    // second derivative of the residual, assembled from the BS gamma and the
    // analytic derivative of the boundary term.
    let gamma_e = (-div * expiry).exp() * pdf / (x * vol * expiry.sqrt());
    let aa = (-div * expiry).exp() * pdf / (vol * expiry.sqrt());
    let fxx = -gamma_e + (-aa / x + aa * d1 / (x * vol * expiry.sqrt())) / q1;
    let star = ift_lift_1d(root, residual, fx, fxx.val());
    let ds = ((star / strike).ln() + (rate - div + v2 * 0.5) * expiry) / (vol * expiry.sqrt());
    let acoef = -star / q1 * (T::one() - (-div * expiry).exp() * (-ds).norm_cdf());
    if spot.val() <= root {
        strike - spot
    } else {
        european(spot, strike, rate, div, vol, expiry, false) + acoef * (spot / star).pow(q1)
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
            critical_put(spot, strike, rate, div, vol, expiry)
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
            baw_put(spot, strike, rate, div, vol, expiry)
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
    // Early exercise is unnecessary only when b=r-q >= r (q<=0).  In
    // particular r=0 with q>0 still has an optimal finite exercise boundary.
    if t.val() <= 0.0 || v.val() <= 0.0 || q.val() <= 0.0 {
        return european(s, k, r, q, v, t, true);
    }
    let v2 = v * v;
    let beta = (T::from_f64(0.5) - (r - q) / v2)
        + (((r - q) / v2 - T::from_f64(0.5)).powf(2.0) + T::from_f64(2.0) * r / v2).sqrt();
    let binf = beta / (beta - 1.0) * k;
    let b = r - q;
    let b0 = if r.val() <= 0.0 {
        k
    } else {
        k.max_val(r / q.max_val(T::from_f64(1e-12)) * k)
    };
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
#[allow(clippy::too_many_arguments)]
fn bs_psi<T: Real>(s: T, t: T, gamma: T, h: T, i2: T, i1: T, t1: T, r: T, b: T, v: T) -> T {
    // Two-time correction from Bjerksund--Stensland (2002).  The four
    // bivariate-normal terms are retained in AD form; no finite differences or
    // primal extraction is used in the pricing graph.
    let st = v * t.sqrt();
    let st1 = v * t1.sqrt();
    let rho = (t1 / t).sqrt();
    let kappa = T::from_f64(2.0) * b / (v * v) + T::from_f64(2.0) * gamma - T::one();
    let lambda = -r + gamma * b + gamma * (gamma - T::one()) * v * v * 0.5;
    let e1 = (-(s / h).ln() - (b + (gamma - 0.5) * v * v) * t) / st;
    let e2 = e1 - T::from_f64(2.0) * (i2 / s).ln() / st;
    let e3 = (-(s / i1).ln() - (b + (gamma - 0.5) * v * v) * t) / st;
    let e4 = e3 - T::from_f64(2.0) * (i1 / s).ln() / st;
    let f1 = (-(s / h).ln() - (b + (gamma - 0.5) * v * v) * t1) / st1;
    let f2 = f1 - T::from_f64(2.0) * (i2 / s).ln() / st1;
    let f3 = (-(s / i1).ln() - (b + (gamma - 0.5) * v * v) * t1) / st1;
    let f4 = f3 - T::from_f64(2.0) * (i1 / s).ln() / st1;
    let p11 = ql_math::bivariate_normal_cdf(-e1, -f1, rho);
    let p12 = ql_math::bivariate_normal_cdf(-e2, -f2, rho);
    let p13 = ql_math::bivariate_normal_cdf(-e3, -f3, rho);
    let p14 = ql_math::bivariate_normal_cdf(-e4, -f4, rho);
    lambda.exp() * (p11 - (i2 / s).pow(kappa) * p12 - (i1 / s).pow(kappa) * (p13 - p14))
}
fn bs2002_call<T: Real>(s: T, k: T, r: T, q: T, v: T, t: T) -> T {
    if t.val() <= 0.0 || v.val() <= 0.0 || q.val() <= 0.0 {
        return european(s, k, r, q, v, t, true);
    }
    let t1 = (T::from_f64(5.0).sqrt() - T::one()) * t * 0.5;
    let v2 = v * v;
    let b = r - q;
    let beta = (T::from_f64(0.5) - b / v2)
        + ((b / v2 - T::from_f64(0.5)).powf(2.0) + T::from_f64(2.0) * r / v2).sqrt();
    let binf = beta / (beta - T::one()) * k;
    let b0 = if r.val() <= 0.0 {
        k
    } else {
        k.max_val(r / q.max_val(T::from_f64(1e-12)) * k)
    };
    let denom = (binf - b0).max_val(T::from_f64(1e-12)) * b0;
    let h1 = -(b * t1 + T::from_f64(2.0) * v * t1.sqrt()) * k * k / denom;
    let h2 = -(b * t + T::from_f64(2.0) * v * t.sqrt()) * k * k / denom;
    let i1 = b0 + (binf - b0) * (T::one() - h1.exp());
    let i2 = b0 + (binf - b0) * (T::one() - h2.exp());
    let a1 = (i1 - k) / i1.pow(beta);
    let a2 = (i2 - k) / i2.pow(beta);
    if s.val() >= i2.val() {
        return s - k;
    }
    let phi = |gamma: T, h: T, i: T| bs_phi(s, t1, gamma, h, i, r, b, v);
    a2 * s.pow(beta) - a2 * phi(beta, i2, i2) + phi(T::one(), i2, i2)
        - phi(T::one(), i1, i2)
        - k * phi(T::zero(), i2, i2)
        + k * phi(T::zero(), i1, i2)
        + a1 * phi(beta, i1, i2)
        - a1 * bs_psi(s, t, beta, i2, i2, i1, t1, r, b, v)
        + bs_psi(s, t, T::one(), i2, i2, i1, t1, r, b, v)
        - bs_psi(s, t, T::one(), k, i2, i1, t1, r, b, v)
        - k * bs_psi(s, t, T::zero(), i2, i2, i1, t1, r, b, v)
        + k * bs_psi(s, t, T::zero(), k, i2, i1, t1, r, b, v)
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
            let alpha = r * 2.0 / v2;
            let beta = n;
            let h = T::one() - (-r * t).exp();
            let alpha_over_h = if h.val().abs() < 1e-10 {
                T::from_f64(2.0) / v2
            } else {
                alpha / h
            };
            let disc = ((beta - T::one()) * (beta - T::one()) + alpha_over_h * 4.0).sqrt();
            let q2 = (-(beta - T::one()) + disc) * 0.5;
            let lambda_prime = if h.val().abs() < 1e-10 {
                T::zero()
            } else {
                -alpha / (h * h * disc)
            };
            let residual = |x: T| {
                let d1 = ((x / k).ln() + (r - q + v2 * 0.5) * t) / (v * t.sqrt());
                x - k
                    - european(x, k, r, q, v, t, true)
                    - x / q2 * (T::one() - (-q * t).exp() * d1.norm_cdf())
            };
            let star0 = root.max(k.val() * (1.0 + 1e-10));
            let star_t = T::from_f64(star0);
            let d1s0 = ((star_t / k).ln() + (r - q + v2 * 0.5) * t) / (v * t.sqrt());
            let dlt = (-q * t).exp() * d1s0.norm_cdf();
            let pdf = d1s0.norm_pdf();
            let f_x = (T::one() - T::one() / q2) * (T::one() - dlt)
                + (-q * t).exp() * pdf / (q2 * v * t.sqrt());
            let gam = (-q * t).exp() * pdf / (star_t * v * t.sqrt());
            let aa = (-q * t).exp() * pdf / (v * t.sqrt());
            let f_xx = -gam - (-aa / star_t + aa * d1s0 / (star_t * v * t.sqrt())) / q2;
            let star = ift_lift_1d(star0, residual, f_x, f_xx.val());
            let euro_star = european(star, k, r, q, v, t, true);
            let h_a = star - k - euro_star;
            if s.val() >= root {
                return s - k;
            }
            let ell = (s / star).ln();
            let d_ve_dt = -q * star * (-q * t).exp() * d1s0.norm_cdf()
                + r * k * (-r * t).exp() * (d1s0 - v * t.sqrt()).norm_cdf()
                + star * (-q * t).exp() * pdf * v / (t.sqrt() * 2.0);
            let ve_h = if r.val().abs() < 1e-10 {
                T::zero()
            } else {
                d_ve_dt / (r * (-r * t).exp())
            };
            let denom_l = q2 * 2.0 + beta - 1.0;
            let bju = (T::one() - h) * alpha * lambda_prime / (denom_l * 2.0);
            let cju = -(T::one() - h) / denom_l
                * (alpha * ve_h / h_a + alpha_over_h + alpha * lambda_prime / (denom_l * 2.0));
            let den = T::one() - bju * ell * ell - cju * ell;
            european(s, k, r, q, v, t, true) + h_a * (s / star).pow(q2) / den
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
