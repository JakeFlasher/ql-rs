use ql_ad::Real;
use ql_math::bivariate_normal_cdf;
use ql_vol::black_scholes_price;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum RainbowKind {
    CallOnMin,
    CallOnMax,
    PutOnMin,
    PutOnMax,
}

pub struct StulzRainbowEngine;
impl StulzRainbowEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        s1: T,
        s2: T,
        k: T,
        r: T,
        q1: T,
        q2: T,
        v1: T,
        v2: T,
        rho: T,
        t: T,
        kind: RainbowKind,
    ) -> T {
        if k.val() <= 0.0 {
            return match kind {
                RainbowKind::CallOnMin => Self::expected_min(s1, s2, q1, q2, v1, v2, rho, t),
                RainbowKind::CallOnMax => s1 * (-q1 * t).exp() + s2 * (-q2 * t).exp(),
                RainbowKind::PutOnMin | RainbowKind::PutOnMax => T::zero(),
            };
        }
        // Stulz (1982) image formula for a call on the minimum.
        let variance = (v1 * v1 + v2 * v2 - rho * v1 * v2 * 2.0).max_val(T::from_f64(1e-16));
        let sig = variance.sqrt();
        let sig_t = (variance * t).sqrt().max_val(T::from_f64(1e-8));
        let d = ((s1 / s2).ln() + (q2 - q1 + variance * 0.5) * t) / sig_t;
        let st1 = v1 * t.sqrt();
        let st2 = v2 * t.sqrt();
        let y1 = ((s1 / k).ln() + (r - q1 + v1 * v1 * 0.5) * t) / st1;
        let y2 = ((s2 / k).ln() + (r - q2 + v2 * v2 * 0.5) * t) / st2;
        let rho1 = (v1 - rho * v2) / sig;
        let rho2 = (v2 - rho * v1) / sig;
        let cmin = s1 * (-q1 * t).exp() * bivariate_normal_cdf(y1, -d, -rho1)
            + s2 * (-q2 * t).exp() * bivariate_normal_cdf(y2, d - sig_t, -rho2)
            - k * (-r * t).exp() * bivariate_normal_cdf(y1 - st1, y2 - st2, rho);
        let c1 = black_scholes_price(s1, k, v1, t, r, q1, true);
        let c2 = black_scholes_price(s2, k, v2, t, r, q2, true);
        let cmax = c1 + c2 - cmin;
        let emin = Self::expected_min(s1, s2, q1, q2, v1, v2, rho, t);
        let emax = s1 * (-q1 * t).exp() + s2 * (-q2 * t).exp() - emin;
        match kind {
            RainbowKind::CallOnMin => cmin,
            RainbowKind::CallOnMax => cmax,
            RainbowKind::PutOnMin => cmin - emin + k * (-r * t).exp(),
            RainbowKind::PutOnMax => cmax - emax + k * (-r * t).exp(),
        }
    }
    #[allow(clippy::too_many_arguments)]
    fn expected_min<T: Real>(s1: T, s2: T, q1: T, q2: T, v1: T, v2: T, rho: T, t: T) -> T {
        let var = (v1 * v1 + v2 * v2 - rho * 2.0 * v1 * v2).max_val(T::from_f64(1e-14));
        let st = var.sqrt() * t.sqrt();
        if st.val() < 1e-8 {
            return s1.min_val(s2);
        }
        let d = ((s1 / s2).ln() + (q2 - q1 + var * 0.5) * t) / st;
        s1 * (-q1 * t).exp() * d.norm_cdf() + s2 * (-q2 * t).exp() * (-d + st).norm_cdf()
    }
}

pub struct TwoAssetCorrelationEngine;
impl TwoAssetCorrelationEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        s1: T,
        s2: T,
        k1: T,
        k2: T,
        r: T,
        q1: T,
        q2: T,
        v1: T,
        v2: T,
        rho: T,
        t: T,
        is_call: bool,
    ) -> T {
        let st1 = v1 * t.sqrt();
        let st2 = v2 * t.sqrt();
        let y1 = ((s1 / k1).ln() + (r - q1 - v1 * v1 * 0.5) * t) / st1;
        let y2 = ((s2 / k2).ln() + (r - q2 - v2 * v2 * 0.5) * t) / st2;
        if is_call {
            (s2 * (-q2 * t).exp() * bivariate_normal_cdf(y2 + st2, y1 + rho * st2, rho)
                - k2 * (-r * t).exp() * bivariate_normal_cdf(y2, y1, rho))
            .max_val(T::zero())
        } else {
            (k2 * (-r * t).exp() * bivariate_normal_cdf(-y2, -y1, rho)
                - s2 * (-q2 * t).exp() * bivariate_normal_cdf(-y2 - st2, -y1 - rho * st2, rho))
            .max_val(T::zero())
        }
    }
}

pub struct DengLiZhouSpreadBasketEngine;
impl DengLiZhouSpreadBasketEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price_spread<T: Real>(
        s1: T,
        s2: T,
        k: T,
        r: T,
        q1: T,
        q2: T,
        v1: T,
        v2: T,
        rho: T,
        t: T,
        is_call: bool,
    ) -> T {
        if k.val().abs() < 1e-14 {
            let vv = (v1 * v1 + v2 * v2 - rho * 2.0 * v1 * v2).max_val(T::from_f64(1e-14));
            let st = vv.sqrt() * t.sqrt();
            let d1 = ((s1 / s2).ln() + (q2 - q1 + vv * 0.5) * t) / st;
            let d2 = d1 - st;
            let call = s1 * (-q1 * t).exp() * d1.norm_cdf() - s2 * (-q2 * t).exp() * d2.norm_cdf();
            return if is_call {
                call
            } else {
                call - s1 * (-q1 * t).exp() + s2 * (-q2 * t).exp()
            };
        }
        let f1 = s1 * ((r - q1) * t).exp();
        let f2 = s2 * ((r - q2) * t).exp();
        // Kirk's first-order boundary approximation, with forward prices
        // discounted exactly once.  The spread variance has the negative
        // covariance term required by Var(S1-S2).
        let f = f1 / (f1 + f2).max_val(T::from_f64(1e-14));
        let v = (v1 * v1 * f * f + v2 * v2 * (T::one() - f) * (T::one() - f)
            - rho * v1 * v2 * f * (T::one() - f) * 2.0)
            .max_val(T::from_f64(1e-16))
            .sqrt();
        ql_vol::black76_price(f1, f2 + k, v, t, (-r * t).exp(), is_call)
    }
}

pub struct LevyBasketEngine;
impl LevyBasketEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price_basket<T: Real>(
        spots: &[T],
        weights: &[T],
        vols: &[T],
        divs: &[T],
        corr: &[Vec<T>],
        strike: T,
        rate: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        assert_eq!(spots.len(), weights.len());
        let mut forward = T::zero();
        let mut var = T::zero();
        for i in 0..spots.len() {
            forward += weights[i] * spots[i] * ((rate - divs[i]) * expiry).exp();
        }
        for i in 0..spots.len() {
            for j in 0..spots.len() {
                var += weights[i]
                    * spots[i]
                    * ((rate - divs[i]) * expiry).exp()
                    * weights[j]
                    * spots[j]
                    * ((rate - divs[j]) * expiry).exp()
                    * ((corr[i][j] * vols[i] * vols[j] * expiry).exp() - T::one());
            }
        }
        let ratio = (T::one()
            + var.max_val(T::zero()) / (forward * forward).max_val(T::from_f64(1e-24)))
        .ln();
        let vol = (ratio / expiry.max_val(T::from_f64(1e-12)))
            .max_val(T::zero())
            .sqrt();
        ql_vol::black76_price(
            forward,
            strike,
            vol,
            expiry,
            (-rate * expiry).exp(),
            is_call,
        )
    }
}
