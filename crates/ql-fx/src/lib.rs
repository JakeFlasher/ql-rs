use ql_ad::Real;
use ql_math::inv_norm_cdf;
use thiserror::Error;

pub mod tarf;
pub mod xccs;
pub use tarf::*;
pub use xccs::*;

#[derive(Debug, Error, Clone)]
pub enum FxError {
    #[error("invalid FX smile quote")]
    InvalidQuote,
}
#[inline]
pub fn garman_kohlhagen_price<T: Real>(
    spot: T,
    strike: T,
    vol: T,
    expiry: T,
    r_dom: T,
    r_for: T,
    is_call: bool,
) -> T {
    if expiry.val() <= 0. {
        return (if is_call {
            spot - strike
        } else {
            strike - spot
        })
        .max_val(T::zero());
    }
    if vol.val() <= 0. {
        // With zero volatility the terminal spot is deterministic.  The
        // payoff is therefore evaluated at the forward and discounted at the
        // domestic rate; using max(spot - strike) * DF_for would incorrectly
        // ignore the carry over a positive expiry.
        let f = spot * ((r_dom - r_for) * expiry).exp();
        let df = (-r_dom * expiry).exp();
        return df * (if is_call { f - strike } else { strike - f }).max_val(T::zero());
    }
    let f = spot * ((r_dom - r_for) * expiry).exp();
    let df = (-r_dom * expiry).exp();
    let s = vol * expiry.sqrt();
    let d1 = ((f / strike).ln() + vol * vol * expiry * 0.5) / s;
    let d2 = d1 - s;
    let sign = if is_call { 1. } else { -1. };
    // The option orientation multiplies the *whole* Black payoff.  Putting
    // the sign on both notionals (as the old expression did) reverses the
    // put value and violates put-call parity.
    df * T::from_f64(sign)
        * (f * (T::from_f64(sign) * d1).norm_cdf() - strike * (T::from_f64(sign) * d2).norm_cdf())
}
#[inline]
pub fn garman_kohlhagen_spot_delta<T: Real>(
    spot: T,
    strike: T,
    vol: T,
    expiry: T,
    r_dom: T,
    r_for: T,
    is_call: bool,
) -> T {
    let f = spot * ((r_dom - r_for) * expiry).exp();
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        let sign = if is_call { 1.0 } else { -1.0 };
        // Deterministic terminal payoff delta; at the kink use the symmetric
        // one-sided convention.  Avoid evaluating d1/vol when vol is zero.
        let indicator = if f.val() > strike.val() {
            1.0
        } else if f.val() < strike.val() {
            0.0
        } else {
            0.5
        };
        return T::from_f64(sign * indicator) * (-r_for * expiry).exp();
    }
    let s = vol * expiry.sqrt();
    let d1 = ((f / strike).ln() + vol * vol * expiry * 0.5) / s;
    let sign = if is_call { 1. } else { -1. };
    T::from_f64(sign) * (-r_for * expiry).exp() * (T::from_f64(sign) * d1).norm_cdf()
}
#[inline]
pub fn garman_kohlhagen_vega<T: Real>(
    spot: T,
    strike: T,
    vol: T,
    expiry: T,
    r_dom: T,
    r_for: T,
) -> T {
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        return T::zero();
    }
    let f = spot * ((r_dom - r_for) * expiry).exp();
    let s = vol * expiry.sqrt();
    let d1 = ((f / strike).ln() + vol * vol * expiry * 0.5) / s;
    (-r_dom * expiry).exp() * f * d1.norm_pdf() * expiry.sqrt()
}
pub fn inv_norm_cdf_real<T: Real>(p: T) -> T {
    let mut z = T::from_f64(inv_norm_cdf(p.val()));
    for _ in 0..2 {
        z -= (z.norm_cdf() - p) / z.norm_pdf();
    }
    z
}
#[derive(Copy, Clone, Debug)]
pub struct FxSmileQuote<T: Real> {
    pub atm_vol: T,
    pub rr_25d: T,
    pub bf_25d: T,
}
impl<T: Real> FxSmileQuote<T> {
    pub fn vol_25d_put(&self) -> T {
        self.atm_vol + self.bf_25d - self.rr_25d * 0.5
    }
    pub fn vol_25d_call(&self) -> T {
        self.atm_vol + self.bf_25d + self.rr_25d * 0.5
    }
}
#[derive(Clone, Debug)]
pub struct FxVannaVolgaSmile<T: Real> {
    pub spot: T,
    pub r_dom: T,
    pub r_for: T,
    pub expiry: T,
    pub quote: FxSmileQuote<T>,
    pub k1: T,
    pub k2: T,
    pub k3: T,
    pub sigma1: T,
    pub sigma2: T,
    pub sigma3: T,
}
impl<T: Real> FxVannaVolgaSmile<T> {
    pub fn new(spot: T, r_dom: T, r_for: T, expiry: T, quote: FxSmileQuote<T>) -> Self {
        let f = spot * ((r_dom - r_for) * expiry).exp();
        let root = inv_norm_cdf_real(T::from_f64(0.25) * (r_for * expiry).exp());
        let s1 = quote.vol_25d_put();
        let s2 = quote.atm_vol;
        let s3 = quote.vol_25d_call();
        let k2 = f * (s2 * s2 * expiry * 0.5).exp();
        let k1 = f * (s1 * expiry.sqrt() * root + s1 * s1 * expiry * 0.5).exp();
        let k3 = f * (-s3 * expiry.sqrt() * root + s3 * s3 * expiry * 0.5).exp();
        Self {
            spot,
            r_dom,
            r_for,
            expiry,
            quote,
            k1,
            k2,
            k3,
            sigma1: s1,
            sigma2: s2,
            sigma3: s3,
        }
    }
    fn weights(&self, k: T) -> (T, T, T) {
        let x1 = self.k1.ln();
        let x2 = self.k2.ln();
        let x3 = self.k3.ln();
        let x = k.ln();
        let l1 = (x2 - x) * (x3 - x) / ((x2 - x1) * (x3 - x1));
        let l2 = (x - x1) * (x3 - x) / ((x2 - x1) * (x3 - x2));
        let l3 = (x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2));
        (l1, l2, l3)
    }
    pub fn implied_vol_first_order(&self, k: T) -> T {
        let (a, b, c) = self.weights(k);
        a * self.sigma1 + b * self.sigma2 + c * self.sigma3
    }
    pub fn implied_vol(&self, k: T) -> T {
        let (l1, l2, l3) = self.weights(k);
        let first = l1 * self.sigma1 + l2 * self.sigma2 + l3 * self.sigma3;
        let forward = self.spot * ((self.r_dom - self.r_for) * self.expiry).exp();
        let d12 = |strike: T, sigma: T| {
            let s = sigma * self.expiry.sqrt();
            let d1 = ((forward / strike).ln() + sigma * sigma * self.expiry * 0.5) / s;
            (d1, d1 - s)
        };
        let (d1, d2) = d12(k, self.sigma2);
        // Castagna--Mercurio's second-order correction is expressed in
        // vega/volga units at the ATM volatility.  Using the pillar
        // volatilities here breaks the exact interpolation property and
        // requires non-differentiable special cases at the pillars.
        let (d1p, d2p) = d12(self.k1, self.sigma2);
        let (d1c, d2c) = d12(self.k3, self.sigma2);
        let dprod = d1 * d2;
        let d2corr = l1 * d1p * d2p * (self.sigma1 - self.sigma2).powf(2.0)
            + l3 * d1c * d2c * (self.sigma3 - self.sigma2).powf(2.0);
        let num = self.sigma2 * (first - self.sigma2) * 2.0 + d2corr;
        // Rationalising the Castagna--Mercurio root removes the removable
        // singularity at d1*d2=0 while retaining all AD lanes.
        let rad = (self.sigma2 * self.sigma2 + dprod * num)
            .max_val(T::from_f64(1e-16))
            .sqrt();
        self.sigma2 + num / (self.sigma2 + rad)
    }
    pub fn price_vanna_volga(&self, k: T, is_call: bool) -> T {
        let base = garman_kohlhagen_price(
            self.spot,
            k,
            self.sigma2,
            self.expiry,
            self.r_dom,
            self.r_for,
            is_call,
        );
        let (a, _, c) = self.weights(k);
        let v = garman_kohlhagen_vega(
            self.spot,
            k,
            self.sigma2,
            self.expiry,
            self.r_dom,
            self.r_for,
        );
        let v1 = garman_kohlhagen_vega(
            self.spot,
            self.k1,
            self.sigma2,
            self.expiry,
            self.r_dom,
            self.r_for,
        );
        let v3 = garman_kohlhagen_vega(
            self.spot,
            self.k3,
            self.sigma2,
            self.expiry,
            self.r_dom,
            self.r_for,
        );
        base + v / v1
            * a
            * (garman_kohlhagen_price(
                self.spot,
                self.k1,
                self.sigma1,
                self.expiry,
                self.r_dom,
                self.r_for,
                is_call,
            ) - garman_kohlhagen_price(
                self.spot,
                self.k1,
                self.sigma2,
                self.expiry,
                self.r_dom,
                self.r_for,
                is_call,
            ))
            + v / v3
                * c
                * (garman_kohlhagen_price(
                    self.spot,
                    self.k3,
                    self.sigma3,
                    self.expiry,
                    self.r_dom,
                    self.r_for,
                    is_call,
                ) - garman_kohlhagen_price(
                    self.spot,
                    self.k3,
                    self.sigma2,
                    self.expiry,
                    self.r_dom,
                    self.r_for,
                    is_call,
                ))
    }
}
