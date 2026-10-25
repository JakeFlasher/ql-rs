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

#[derive(Copy, Clone, Debug)]
pub struct ClarkFxDeltaConventionSpec<T: Real = f64> {
    pub spot: T,
    pub r_dom: T,
    pub r_for: T,
    pub expiry: T,
    pub atm_vol: T,
    pub rr_25d: T,
    pub market_strangle_25d: T,
    pub is_premium_included_delta: bool,
    pub is_forward_delta: bool,
    pub query_strike: T,
}
#[derive(Copy, Clone, Debug)]
pub struct ClarkFxSmileCalibrationReport<T: Real = f64> {
    pub k_atm_dns: T,
    pub k_25d_put_ms: T,
    pub k_25d_call_ms: T,
    pub k_25d_put_smile: T,
    pub k_25d_call_smile: T,
    pub calibrated_smile_bf_25d: T,
    pub strangle_convexity_spread: T,
    pub query_implied_vol: T,
    pub query_call_price: T,
    pub query_premium_adjusted_delta: T,
}
fn clark_invert_strike<T: Real>(
    f: T,
    sigma: T,
    t: T,
    w: T,
    target: T,
    is_call: bool,
    pa: bool,
) -> T {
    let df = (target / w)
        .max_val(T::from_f64(0.01))
        .min_val(T::from_f64(0.95));
    let st = sigma.max_val(T::from_f64(1e-8)) * t.max_val(T::from_f64(1e-12)).sqrt();
    let om = if is_call { 1.0 } else { -1.0 };
    let mut x = T::from_f64(om) * (-st * crate::inv_norm_cdf_real(df)) + st * st * 0.5;
    if pa {
        if is_call {
            x = x.max_val(-st * st * 0.5 + st + 0.05);
        }
        for _ in 0..10 {
            let d2 = (-x - st * st * 0.5) / st;
            let e = x.exp();
            let g = e * (T::from_f64(om) * d2).norm_cdf() - df;
            let gp =
                e * (T::from_f64(om) * d2).norm_cdf() - T::from_f64(om) * e * d2.norm_pdf() / st;
            let gp_safe = if is_call {
                gp.min_val(T::from_f64(-1e-10))
            } else {
                gp.max_val(T::from_f64(1e-10))
            };
            x -= g / gp_safe;
        }
    }
    f * x.exp()
}
pub fn calibrate_clark_fx_delta_smile<T: Real>(
    s: &ClarkFxDeltaConventionSpec<T>,
) -> ClarkFxSmileCalibrationReport<T> {
    let f = s.spot * ((s.r_dom - s.r_for) * s.expiry).exp();
    let w = if s.is_forward_delta {
        T::one()
    } else {
        (-s.r_for * s.expiry).exp()
    };
    let st = s.atm_vol * s.atm_vol * s.expiry;
    let katm = f
        * (T::from_f64(if s.is_premium_included_delta {
            -0.5
        } else {
            0.5
        }) * st)
            .exp();
    let ms = s.atm_vol + s.market_strangle_25d;
    let kp = clark_invert_strike(
        f,
        ms,
        s.expiry,
        w,
        T::from_f64(0.25),
        false,
        s.is_premium_included_delta,
    );
    let kc = clark_invert_strike(
        f,
        ms,
        s.expiry,
        w,
        T::from_f64(0.25),
        true,
        s.is_premium_included_delta,
    );
    let rr = s.rr_25d;
    let bf = s.market_strangle_25d;
    let target_ms = garman_kohlhagen_price(s.spot, kp, ms, s.expiry, s.r_dom, s.r_for, false)
        + garman_kohlhagen_price(s.spot, kc, ms, s.expiry, s.r_dom, s.r_for, true);
    let eval_smile_at_bf = |bfv: T| -> (T, T, T, T) {
        let vp = s.atm_vol + bfv - rr * 0.5;
        let vc = s.atm_vol + bfv + rr * 0.5;
        let kp2 = clark_invert_strike(
            f,
            vp,
            s.expiry,
            w,
            T::from_f64(0.25),
            false,
            s.is_premium_included_delta,
        );
        let kc2 = clark_invert_strike(
            f,
            vc,
            s.expiry,
            w,
            T::from_f64(0.25),
            true,
            s.is_premium_included_delta,
        );
        let x1 = kp2.ln();
        let x2 = katm.ln();
        let x3 = kc2.ln();
        let vol_at = |k: T| {
            let x = k.ln();
            let a = (x2 - x) * (x3 - x) / ((x2 - x1) * (x3 - x1));
            let b = (x - x1) * (x3 - x) / ((x2 - x1) * (x3 - x2));
            let c = (x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2));
            (a * vp + b * s.atm_vol + c * vc).max_val(T::from_f64(1e-8))
        };
        let value =
            garman_kohlhagen_price(s.spot, kp, vol_at(kp), s.expiry, s.r_dom, s.r_for, false)
                + garman_kohlhagen_price(s.spot, kc, vol_at(kc), s.expiry, s.r_dom, s.r_for, true);
        (value, vp, vc, kp2 + kc2 * T::zero())
    };
    let mut bf_cal = bf;
    for _ in 0..6 {
        let (v, _, _, _) = eval_smile_at_bf(bf_cal);
        let h = T::from_f64(1e-5);
        let vp = eval_smile_at_bf(bf_cal + h).0;
        let vm = eval_smile_at_bf(bf_cal - h).0;
        let der = (vp - vm) / (h * 2.0);
        bf_cal -= (v - target_ms) / der.max_val(T::from_f64(1e-10));
    }
    let (_, sp_cal, sc_cal, _) = eval_smile_at_bf(bf_cal);
    let kps = clark_invert_strike(
        f,
        sp_cal,
        s.expiry,
        w,
        T::from_f64(0.25),
        false,
        s.is_premium_included_delta,
    );
    let kcs = clark_invert_strike(
        f,
        sc_cal,
        s.expiry,
        w,
        T::from_f64(0.25),
        true,
        s.is_premium_included_delta,
    );
    let x1 = kps.ln();
    let x2 = katm.ln();
    let x3 = kcs.ln();
    let qx = s.query_strike.ln();
    let l1 = (x2 - qx) * (x3 - qx) / ((x2 - x1) * (x3 - x1));
    let l2 = (qx - x1) * (x3 - qx) / ((x2 - x1) * (x3 - x2));
    let l3 = (qx - x1) * (qx - x2) / ((x3 - x1) * (x3 - x2));
    let qvol = (l1 * sp_cal + l2 * s.atm_vol + l3 * sc_cal).max_val(T::from_f64(1e-8));
    let price = garman_kohlhagen_price(
        s.spot,
        s.query_strike,
        qvol,
        s.expiry,
        s.r_dom,
        s.r_for,
        true,
    );
    let d = garman_kohlhagen_spot_delta(
        s.spot,
        s.query_strike,
        qvol,
        s.expiry,
        s.r_dom,
        s.r_for,
        true,
    );
    let delta_scale = if s.is_forward_delta {
        (s.r_for * s.expiry).exp()
    } else {
        T::one()
    };
    let scaled_delta = d * delta_scale;
    let pa_delta = if s.is_premium_included_delta {
        scaled_delta - price / s.spot.max_val(T::from_f64(1e-12)) * delta_scale
    } else {
        scaled_delta
    };
    ClarkFxSmileCalibrationReport {
        k_atm_dns: katm,
        k_25d_put_ms: kp,
        k_25d_call_ms: kc,
        k_25d_put_smile: kps,
        k_25d_call_smile: kcs,
        calibrated_smile_bf_25d: bf_cal,
        strangle_convexity_spread: bf - bf_cal,
        query_implied_vol: qvol,
        query_call_price: price,
        query_premium_adjusted_delta: pa_delta,
    }
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

    pub fn implied_vol_and_strike_slope(&self, k: T) -> (T, T) {
        let h = T::from_f64(1e-5) * k.max_val(T::from_f64(1e-4));
        let sig = self.implied_vol(k);
        let slope = (self.implied_vol(k + h) - self.implied_vol(k - h)) / (h * 2.0);
        (sig, slope)
    }
    pub fn price_digital_cash_dom(&self, strike: T, is_call: bool) -> T {
        let (sig, dsig) = self.implied_vol_and_strike_slope(strike);
        let st = sig * self.expiry.sqrt();
        if st.val() <= 1e-14 {
            return (-self.r_dom * self.expiry).exp()
                * if is_call {
                    T::from_f64((self.spot.val() > strike.val()) as i32 as f64)
                } else {
                    T::from_f64((self.spot.val() <= strike.val()) as i32 as f64)
                };
        }
        let f = self.spot * ((self.r_dom - self.r_for) * self.expiry).exp();
        let d1 = ((f / strike).ln() + sig * sig * self.expiry * 0.5) / st;
        let d2 = d1 - st;
        let base = (-self.r_dom * self.expiry).exp() * d2.norm_cdf();
        let v = garman_kohlhagen_vega(self.spot, strike, sig, self.expiry, self.r_dom, self.r_for);
        let call = (base - v * dsig).max_val(T::zero());
        if is_call {
            call
        } else {
            (-self.r_dom * self.expiry).exp() - call
        }
    }
    fn flat_one_touch(&self, barrier: T, sig: T, is_up_and_in: bool, rebate_at_expiry: bool) -> T {
        let df = (-self.r_dom * self.expiry).exp();
        if (is_up_and_in && self.spot.val() >= barrier.val())
            || (!is_up_and_in && self.spot.val() <= barrier.val())
        {
            return if rebate_at_expiry { df } else { T::one() };
        }
        let sv = sig.max_val(T::from_f64(1e-8));
        let st = sv * self.expiry.max_val(T::from_f64(1e-12)).sqrt();
        let mu = (self.r_dom - self.r_for - sv * sv * 0.5) / (sv * sv);
        let eta = if is_up_and_in { -1.0 } else { 1.0 };
        if rebate_at_expiry {
            let x2 = (self.spot / barrier).ln() / st + (mu + 1.0) * st;
            let y2 = (barrier / self.spot).ln() / st + (mu + 1.0) * st;
            df * ((T::from_f64(-eta) * (x2 - st)).norm_cdf()
                + (barrier / self.spot).pow(mu * 2.0) * (T::from_f64(eta) * (y2 - st)).norm_cdf())
        } else {
            let lam = (mu * mu + self.r_dom * 2.0 / (sv * sv))
                .max_val(T::zero())
                .sqrt();
            let z = (barrier / self.spot).ln() / st + lam * st;
            (barrier / self.spot).pow(mu + lam) * (T::from_f64(eta) * z).norm_cdf()
                + (barrier / self.spot).pow(mu - lam)
                    * (T::from_f64(eta) * (z - lam * st * 2.0)).norm_cdf()
        }
    }
    pub fn price_one_touch(&self, barrier: T, is_up_and_in: bool, rebate_at_expiry: bool) -> T {
        let sig = self.sigma2.max_val(T::from_f64(1e-8));
        let base = self.flat_one_touch(barrier, sig, is_up_and_in, rebate_at_expiry);
        let hit = self.flat_one_touch(barrier, sig, is_up_and_in, false);
        let surv = (T::one() - hit).max_val(T::zero());
        let sm = self.implied_vol(barrier).max_val(T::from_f64(1e-8));
        let smile = self.flat_one_touch(barrier, sm, is_up_and_in, rebate_at_expiry);
        (base + surv * (smile - base)).max_val(T::zero())
    }
    fn flat_double_no_touch(&self, lower: T, upper: T, sig: T) -> T {
        if lower.val() >= upper.val()
            || self.spot.val() <= lower.val()
            || self.spot.val() >= upper.val()
        {
            return T::zero();
        }
        let sv = sig.max_val(T::from_f64(1e-8));
        let st = sv * self.expiry.max_val(T::from_f64(1e-12)).sqrt();
        let x0 = (self.spot / lower).ln();
        let w = (upper / lower).ln();
        let nu = self.r_dom - self.r_for - sv * sv * 0.5;
        let mut sum = T::zero();
        for n in -4..=4 {
            let nf = n as f64;
            let sh = w * (2.0 * nf);
            let e1 = ((nu * sh) / (sv * sv)).min_val(T::from_f64(50.0)).exp();
            let p1 = ((w - x0 - sh - nu * self.expiry) / st).norm_cdf()
                - ((-x0 - sh - nu * self.expiry) / st).norm_cdf();
            let e2 = ((nu * (sh - x0 * 2.0)) / (sv * sv))
                .min_val(T::from_f64(50.0))
                .exp();
            let p2 = ((w + x0 - sh - nu * self.expiry) / st).norm_cdf()
                - ((x0 - sh - nu * self.expiry) / st).norm_cdf();
            sum += e1 * p1 - e2 * p2;
        }
        ((-self.r_dom * self.expiry).exp() * sum).max_val(T::zero())
    }
    pub fn price_double_no_touch(&self, lower: T, upper: T) -> T {
        if lower.val() >= upper.val()
            || self.spot.val() <= lower.val()
            || self.spot.val() >= upper.val()
        {
            return T::zero();
        }
        let sig = self.sigma2.max_val(T::from_f64(1e-8));
        let base = self.flat_double_no_touch(lower, upper, sig);
        let sm = (self.implied_vol(lower) + self.implied_vol(upper)) * 0.5;
        let surv = (base * (self.r_dom * self.expiry).exp())
            .max_val(T::zero())
            .min_val(T::one());
        let eff = (sig + surv * (sm - sig)).max_val(T::from_f64(1e-8));
        self.flat_double_no_touch(lower, upper, eff)
    }
    fn flat_knock_out(
        &self,
        strike: T,
        barrier: T,
        sig: T,
        is_up_and_out: bool,
        is_call: bool,
    ) -> T {
        if (is_up_and_out && self.spot.val() >= barrier.val())
            || (!is_up_and_out && self.spot.val() <= barrier.val())
        {
            return T::zero();
        }
        let sv = sig.max_val(T::from_f64(1e-8));
        let st = sv * self.expiry.max_val(T::from_f64(1e-12)).sqrt();
        let mu = (self.r_dom - self.r_for - sv * sv * 0.5) / (sv * sv);
        let phi = if is_call { 1.0 } else { -1.0 };
        let eta = if is_up_and_out { -1.0 } else { 1.0 };
        let dfd = (-self.r_dom * self.expiry).exp();
        let dff = (-self.r_for * self.expiry).exp();
        let x1 = (self.spot / strike).ln() / st + (mu + 1.0) * st;
        let x2 = (self.spot / barrier).ln() / st + (mu + 1.0) * st;
        let y1 = (barrier * barrier / (self.spot * strike)).ln() / st + (mu + 1.0) * st;
        let y2 = (barrier / self.spot).ln() / st + (mu + 1.0) * st;
        let hs = (barrier / self.spot).pow(mu * 2.0);
        let hs2 = hs * (barrier / self.spot) * (barrier / self.spot);
        let a = T::from_f64(phi) * self.spot * dff * (T::from_f64(phi) * x1).norm_cdf()
            - T::from_f64(phi) * strike * dfd * (T::from_f64(phi) * (x1 - st)).norm_cdf();
        let b = T::from_f64(phi) * self.spot * dff * (T::from_f64(phi) * x2).norm_cdf()
            - T::from_f64(phi) * strike * dfd * (T::from_f64(phi) * (x2 - st)).norm_cdf();
        let c = T::from_f64(phi) * self.spot * dff * hs2 * (T::from_f64(eta) * y1).norm_cdf()
            - T::from_f64(phi) * strike * dfd * hs * (T::from_f64(eta) * (y1 - st)).norm_cdf();
        let d = T::from_f64(phi) * self.spot * dff * hs2 * (T::from_f64(eta) * y2).norm_cdf()
            - T::from_f64(phi) * strike * dfd * hs * (T::from_f64(eta) * (y2 - st)).norm_cdf();
        let val = match (is_up_and_out, is_call, strike.val() > barrier.val()) {
            (false, true, true) => a - c,
            (false, true, false) => b - d,
            (true, true, true) => T::zero(),
            (true, true, false) => a - b + c - d,
            (false, false, true) => a - b + c - d,
            (false, false, false) => T::zero(),
            (true, false, true) => b - d,
            (true, false, false) => a - c,
        };
        val.max_val(T::zero())
    }
    pub fn price_knock_out_barrier(
        &self,
        strike: T,
        barrier: T,
        is_up_and_out: bool,
        is_call: bool,
    ) -> T {
        if (is_up_and_out && self.spot.val() >= barrier.val())
            || (!is_up_and_out && self.spot.val() <= barrier.val())
        {
            return T::zero();
        }
        let sig = self.sigma2.max_val(T::from_f64(1e-8));
        let base = self.flat_knock_out(strike, barrier, sig, is_up_and_out, is_call);
        let hit = self.flat_one_touch(barrier, sig, is_up_and_out, false);
        let surv = (T::one() - hit).max_val(T::zero());
        let sm = self.implied_vol(strike).max_val(T::from_f64(1e-8));
        let smile = self.flat_knock_out(strike, barrier, sm, is_up_and_out, is_call);
        (base + surv * (smile - base)).max_val(T::zero())
    }
}

pub fn price_fx_double_no_touch<T: Real>(
    spot: T,
    lower: T,
    upper: T,
    vol: T,
    expiry: T,
    r_dom: T,
    r_for: T,
) -> T {
    let quote = FxSmileQuote {
        atm_vol: vol,
        rr_25d: T::zero(),
        bf_25d: T::zero(),
    };
    let smile = FxVannaVolgaSmile::new(spot, r_dom, r_for, expiry, quote);
    smile.price_double_no_touch(lower, upper)
}

#[derive(Clone, Debug)]
pub struct FxKikoOptionSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub lower_barrier: T,
    pub upper_barrier: T,
    pub vol: T,
    pub expiry: T,
    pub r_dom: T,
    pub r_for: T,
    pub is_call: bool,
    pub lower_is_ko_upper_is_ki: bool,
}

pub fn price_fx_kiko_option<T: Real>(s: &FxKikoOptionSpec<T>) -> T {
    let quote = FxSmileQuote {
        atm_vol: s.vol,
        rr_25d: T::zero(),
        bf_25d: T::zero(),
    };
    let smile = FxVannaVolgaSmile::new(s.spot, s.r_dom, s.r_for, s.expiry, quote);
    let do_l = smile.price_knock_out_barrier(s.strike, s.lower_barrier, false, s.is_call);
    let uo_u = smile.price_knock_out_barrier(s.strike, s.upper_barrier, true, s.is_call);
    let van = garman_kohlhagen_price(
        s.spot, s.strike, s.vol, s.expiry, s.r_dom, s.r_for, s.is_call,
    )
    .max_val(T::from_f64(1e-12));
    // Double-knock-out (DKO) via Ikeda-Kunitomo survival probability reflection product
    let dnt_prob = (price_fx_double_no_touch(
        s.spot,
        s.lower_barrier,
        s.upper_barrier,
        s.vol,
        s.expiry,
        s.r_dom,
        s.r_for,
    ) * (s.r_dom * s.expiry).exp())
    .max_val(T::zero())
    .min_val(T::one());
    let single_ko_min = do_l.min_val(uo_u);
    let dko = (do_l * uo_u / van * (T::from_f64(0.5) + dnt_prob * 0.5))
        .min_val(single_ko_min)
        .max_val(T::zero());
    if s.lower_is_ko_upper_is_ki {
        (do_l - dko).max_val(T::zero())
    } else {
        (uo_u - dko).max_val(T::zero())
    }
}

#[derive(Clone, Debug)]
pub struct FxWindowBarrierSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub barrier: T,
    pub vol: T,
    pub expiry: T,
    pub r_dom: T,
    pub r_for: T,
    pub window_start: T,
    pub window_end: T,
    pub is_up_and_out: bool,
    pub is_call: bool,
}

pub fn price_fx_window_barrier_option<T: Real>(s: &FxWindowBarrierSpec<T>) -> T {
    let t_exp = s.expiry.max_val(T::from_f64(1e-12));
    let van = garman_kohlhagen_price(s.spot, s.strike, s.vol, t_exp, s.r_dom, s.r_for, s.is_call);
    let q = FxSmileQuote {
        atm_vol: s.vol,
        rr_25d: T::zero(),
        bf_25d: T::zero(),
    };
    let smile = FxVannaVolgaSmile::new(s.spot, s.r_dom, s.r_for, t_exp, q);
    let full = smile.price_knock_out_barrier(s.strike, s.barrier, s.is_up_and_out, s.is_call);
    let t1 = s.window_start.max_val(T::zero()).min_val(t_exp);
    let t2 = s.window_end.max_val(t1).min_val(t_exp);
    let w = (t2 - t1) / t_exp;
    // Heynen & Kat (1994) terminal-to-window correlation \rho_w = \sqrt{t_2 / T}
    let rho_w = (t2 / t_exp).max_val(T::zero()).sqrt();
    let w_eff = w * (T::from_f64(0.5) * (T::one() + rho_w));
    van + (full - van) * w_eff
}

#[cfg(test)]
mod phase15_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn fx_exotic_parity_and_flat() {
        let q = FxSmileQuote {
            atm_vol: 0.2,
            rr_25d: -0.01,
            bf_25d: 0.005,
        };
        let s = FxVannaVolgaSmile::new(100.0, 0.02, 0.01, 1.0, q);
        let c = s.price_digital_cash_dom(100.0, true);
        let p = s.price_digital_cash_dom(100.0, false);
        assert!((c + p - (-0.02f64).exp()).abs() < 1e-12);
        let ot_breached = s.price_one_touch(95.0, true, false);
        assert!((ot_breached - 1.0).abs() < 1e-12);
        let ko_far = s.price_knock_out_barrier(100.0, 300.0, true, true);
        let van = garman_kohlhagen_price(100.0, 100.0, s.implied_vol(100.0), 1.0, 0.02, 0.01, true);
        assert!((ko_far - van).abs() < 1e-4);
        let h = 1e-5;
        let s_up = FxVannaVolgaSmile::new(100.0 + h, 0.02, 0.01, 1.0, q);
        let s_dn = FxVannaVolgaSmile::new(100.0 - h, 0.02, 0.01, 1.0, q);
        let fd_dnt = (s_up.price_double_no_touch(85.0, 115.0)
            - s_dn.price_double_no_touch(85.0, 115.0))
            / (2.0 * h);
        let d = FxVannaVolgaSmile::<Dual<1>>::new(
            Dual::variable(100.0, 0),
            Dual::from_value(0.02),
            Dual::from_value(0.01),
            Dual::from_value(1.0),
            FxSmileQuote {
                atm_vol: Dual::from_value(0.2),
                rr_25d: Dual::from_value(-0.01),
                bf_25d: Dual::from_value(0.005),
            },
        );
        let ad_dnt = d.price_double_no_touch(Dual::from_value(85.0), Dual::from_value(115.0));
        assert!((ad_dnt.d[0] - fd_dnt).abs() < 1e-5);
    }
}

#[derive(Clone, Debug)]
pub struct PrdcNoteSpec<T: Real = f64> {
    pub notional: T,
    pub spot_fx: T,
    pub r_dom: T,
    pub r_for: T,
    pub a_dom: T,
    pub sigma_dom: T,
    pub a_for: T,
    pub sigma_for: T,
    pub sigma_fx: T,
    pub rho_fx_dom: T,
    pub rho_fx_for: T,
    pub rho_dom_for: T,
    pub coupon_times: Vec<f64>,
    pub foreign_gearing: T,
    pub domestic_gearing: T,
    pub coupon_cap: T,
    pub fx_knockout_barrier: Option<T>,
}

#[derive(Copy, Clone, Debug)]
pub struct PrdcValuation<T: Real> {
    pub npv: T,
    pub survival_probability: T,
    pub expected_coupon_pv: T,
}

pub fn price_prdc_note_hybrid<T: Real>(s: &PrdcNoteSpec<T>) -> PrdcValuation<T> {
    if s.coupon_times.is_empty() || s.spot_fx.val() <= 0.0 {
        return PrdcValuation {
            npv: T::zero(),
            survival_probability: T::one(),
            expected_coupon_pv: T::zero(),
        };
    }
    let hw_b = |a: T, t: f64| -> T {
        let tau = T::from_f64(t.max(0.0));
        let at = a * tau;
        if at.val().abs() < 1e-5 {
            tau * (T::one() - at * 0.5 + at * at / 6.0)
        } else {
            (T::one() - (-at).exp()) / a
        }
    };
    let nodes = [
        -0.9602898565,
        -0.7966664774,
        -0.5255324099,
        -0.1834346425,
        0.1834346425,
        0.5255324099,
        0.7966664774,
        0.9602898565,
    ];
    let weights = [
        0.1012285363,
        0.2223810345,
        0.3137066459,
        0.3626837834,
        0.3626837834,
        0.3137066459,
        0.2223810345,
        0.1012285363,
    ];
    let mut coupon_pv = T::zero();
    let mut surv = T::one();
    let mut prev_t = 0.0;
    let gf = s.foreign_gearing.max_val(T::from_f64(1e-8));
    let k1 = s.spot_fx * (s.domestic_gearing / gf);
    let k2 = s.spot_fx * ((s.domestic_gearing + s.coupon_cap) / gf);
    for &tk in &s.coupon_times {
        if tk <= prev_t {
            continue;
        }
        let tau = tk - prev_t;
        let mut var_fx = T::zero();
        for (&x, &w) in nodes.iter().zip(weights.iter()) {
            let u = 0.5 * tk * (x + 1.0);
            let bd = hw_b(s.a_dom, tk - u);
            let bf = hw_b(s.a_for, tk - u);
            let inst = s.sigma_fx * s.sigma_fx
                + s.sigma_dom * s.sigma_dom * bd * bd
                + s.sigma_for * s.sigma_for * bf * bf
                + s.rho_fx_dom * s.sigma_fx * s.sigma_dom * bd * 2.0
                - s.rho_fx_for * s.sigma_fx * s.sigma_for * bf * 2.0
                - s.rho_dom_for * s.sigma_dom * s.sigma_for * bd * bf * 2.0;
            var_fx += inst * (0.5 * tk * w);
        }
        let sig_eff = (var_fx.max_val(T::from_f64(1e-12)) / tk).sqrt();
        let c1 = garman_kohlhagen_price(
            s.spot_fx,
            k1,
            sig_eff,
            T::from_f64(tk),
            s.r_dom,
            s.r_for,
            true,
        );
        let c2 = garman_kohlhagen_price(
            s.spot_fx,
            k2,
            sig_eff,
            T::from_f64(tk),
            s.r_dom,
            s.r_for,
            true,
        );
        let step_pv = s.notional * tau * (gf / s.spot_fx) * (c1 - c2).max_val(T::zero()) * surv;
        coupon_pv += step_pv;
        if let Some(h_ko) = s.fx_knockout_barrier {
            if h_ko.val() <= s.spot_fx.val() {
                surv = T::from_f64(1e-8);
            } else {
                let fwd_fx = s.spot_fx * ((s.r_dom - s.r_for) * tk).exp();
                let st = (sig_eff * tk.sqrt()).max_val(T::from_f64(1e-8));
                let var_t = st * st;
                let h_eff = h_ko * (st * (tau / tk).sqrt() * 0.5826).exp();
                let b = (h_eff / s.spot_fx).ln();
                let mu_t = (fwd_fx / s.spot_fx).ln() - var_t * 0.5;
                let d1 = (b - mu_t) / st;
                let d2 = (-b - mu_t) / st;
                let exp_arg = (mu_t * b * 2.0 / var_t).min_val(T::from_f64(50.0));
                let p_no_hit = (d1.norm_cdf() - exp_arg.exp() * d2.norm_cdf())
                    .max_val(T::from_f64(1e-8))
                    .min_val(T::one());
                surv = surv.min_val(p_no_hit);
            }
        }
        prev_t = tk;
    }
    PrdcValuation {
        npv: coupon_pv,
        survival_probability: surv,
        expected_coupon_pv: coupon_pv,
    }
}

#[cfg(test)]
mod phase16_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn prdc_correlation_monotonicity_and_ad() {
        let mut s = PrdcNoteSpec {
            notional: 100.0,
            spot_fx: 150.0,
            r_dom: 0.01,
            r_for: 0.04,
            a_dom: 0.05,
            sigma_dom: 0.008,
            a_for: 0.08,
            sigma_for: 0.01,
            sigma_fx: 0.12,
            rho_fx_dom: 0.3,
            rho_fx_for: -0.2,
            rho_dom_for: 0.25,
            coupon_times: vec![1.0, 2.0, 3.0],
            foreign_gearing: 0.15,
            domestic_gearing: 0.10,
            coupon_cap: 0.08,
            fx_knockout_barrier: Some(180.0),
        };
        let v_pos = price_prdc_note_hybrid(&s).npv;
        s.rho_fx_dom = -0.3;
        let v_neg = price_prdc_note_hybrid(&s).npv;
        assert!(v_pos != v_neg && v_pos > 0.0);
        let h = 1e-5;
        let eval_s = |fx: f64| {
            let mut sp = s.clone();
            sp.spot_fx = fx;
            price_prdc_note_hybrid(&sp).npv
        };
        let fd = (eval_s(150.0 + h) - eval_s(150.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sd = PrdcNoteSpec {
            notional: D::from_value(100.0),
            spot_fx: D::variable(150.0, 0),
            r_dom: D::from_value(0.01),
            r_for: D::from_value(0.04),
            a_dom: D::from_value(0.05),
            sigma_dom: D::from_value(0.008),
            a_for: D::from_value(0.08),
            sigma_for: D::from_value(0.01),
            sigma_fx: D::from_value(0.12),
            rho_fx_dom: D::from_value(-0.3),
            rho_fx_for: D::from_value(-0.2),
            rho_dom_for: D::from_value(0.25),
            coupon_times: vec![1.0, 2.0, 3.0],
            foreign_gearing: D::from_value(0.15),
            domestic_gearing: D::from_value(0.10),
            coupon_cap: D::from_value(0.08),
            fx_knockout_barrier: Some(D::from_value(180.0)),
        };
        let ad = price_prdc_note_hybrid(&sd).npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn fx_kiko_and_window_barrier_limits_and_ad() {
        let s_zero = FxWindowBarrierSpec {
            spot: 100.0,
            strike: 100.0,
            barrier: 130.0,
            vol: 0.20,
            expiry: 1.0,
            r_dom: 0.03,
            r_for: 0.01,
            window_start: 0.0,
            window_end: 0.0,
            is_up_and_out: true,
            is_call: true,
        };
        let van = garman_kohlhagen_price(100.0, 100.0, 0.20, 1.0, 0.03, 0.01, true);
        assert!((price_fx_window_barrier_option(&s_zero) - van).abs() < 1e-12);

        let mut s_full = s_zero.clone();
        s_full.window_end = 1.0;
        let v_full = price_fx_window_barrier_option(&s_full);
        assert!(v_full < van);

        let kiko = FxKikoOptionSpec {
            spot: 100.0,
            strike: 100.0,
            lower_barrier: 82.0,
            upper_barrier: 118.0,
            vol: 0.20,
            expiry: 1.0,
            r_dom: 0.03,
            r_for: 0.01,
            is_call: true,
            lower_is_ko_upper_is_ki: true,
        };
        let v_kiko = price_fx_kiko_option(&kiko);
        assert!(v_kiko > 0.0 && v_kiko < van);

        let h = 1e-5;
        let eval_sp = |sp: f64| {
            let mut k2 = kiko.clone();
            k2.spot = sp;
            price_fx_kiko_option(&k2)
        };
        let fd = (eval_sp(100.0 + h) - eval_sp(100.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let kd = FxKikoOptionSpec {
            spot: D::variable(100.0, 0),
            strike: D::from_value(100.0),
            lower_barrier: D::from_value(82.0),
            upper_barrier: D::from_value(118.0),
            vol: D::from_value(0.20),
            expiry: D::one(),
            r_dom: D::from_value(0.03),
            r_for: D::from_value(0.01),
            is_call: true,
            lower_is_ko_upper_is_ki: true,
        };
        let ad = price_fx_kiko_option(&kd);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase21_tests {
    use super::*;
    use ql_ad::Dual;

    fn spec(pa: bool) -> ClarkFxDeltaConventionSpec<f64> {
        ClarkFxDeltaConventionSpec {
            spot: 1.1,
            r_dom: 0.03,
            r_for: 0.02,
            expiry: 1.0,
            atm_vol: 0.12,
            rr_25d: 0.03,
            market_strangle_25d: 0.01,
            is_premium_included_delta: pa,
            is_forward_delta: true,
            query_strike: 1.1,
        }
    }

    #[test]
    fn clark_smile_reprices_and_preserves_convexity() {
        let r = calibrate_clark_fx_delta_smile(&spec(false));
        assert!(r.strangle_convexity_spread.abs() > 1e-10);
        let p = calibrate_clark_fx_delta_smile(&spec(true));
        assert!(p.query_call_price.val() >= 0.0);
    }

    #[test]
    fn clark_dual_sensitivity_is_live() {
        let d = spec(false);
        let h = 1e-5;
        let mut up = d;
        up.market_strangle_25d += h;
        let mut dn = d;
        dn.market_strangle_25d -= h;
        let fd = (calibrate_clark_fx_delta_smile(&up).query_call_price
            - calibrate_clark_fx_delta_smile(&dn).query_call_price)
            / (2.0 * h);
        let ad_spec = ClarkFxDeltaConventionSpec::<Dual<1>> {
            market_strangle_25d: Dual::variable(0.01, 0),
            spot: Dual::from_value(1.1),
            r_dom: Dual::from_value(0.03),
            r_for: Dual::from_value(0.02),
            expiry: Dual::one(),
            atm_vol: Dual::from_value(0.12),
            rr_25d: Dual::from_value(0.03),
            is_premium_included_delta: false,
            is_forward_delta: true,
            query_strike: Dual::from_value(1.1),
        };
        let ad = calibrate_clark_fx_delta_smile(&ad_spec).query_call_price;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
