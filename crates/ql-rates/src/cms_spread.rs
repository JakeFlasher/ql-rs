use ql_ad::Real;
use ql_vol::bachelier_price;

#[derive(Copy, Clone, Debug)]
pub struct CmsSpreadSpec<T: Real> {
    pub fixing_time: T,
    pub payment_time: T,
    pub accrual_factor: T,
    pub discount_df: T,
    pub gearing1: T,
    pub gearing2: T,
    pub spread: T,
    pub swap_rate1: T,
    pub swap_rate2: T,
    pub vol1: T,
    pub vol2: T,
    pub rho: T,
    pub shift1: T,
    pub shift2: T,
    pub mean_reversion: T,
    pub annuity1: T,
    pub annuity2: T,
    pub swap_tenor1: f64,
    pub swap_tenor2: f64,
}

pub struct LognormalCmsSpreadPricer;
impl LognormalCmsSpreadPricer {
    fn convexity_factor<T: Real>(a: T, tenor: f64, annuity: T, dt: T) -> T {
        let u = T::from_f64(tenor.max(0.0));
        let mr = (a * u).abs();
        let b = if mr.val() < 1e-8 {
            u - a * u * u * 0.5 + a * a * u * u * u / 6.0
        } else {
            (T::one() - (-a * u).exp()) / a
        };
        // Linear-TSR annuity convexity factor.  The annuity ratio keeps the
        // payment leg and the underlying swap numeraire in the adjustment.
        b * annuity / (T::one() + annuity * u) * dt.max_val(T::zero())
    }
    fn adjusted_rates<T: Real>(s: &CmsSpreadSpec<T>) -> (T, T) {
        let tf = s.fixing_time.max_val(T::zero());
        let dt = (s.payment_time - s.fixing_time).max_val(T::zero());
        let c1 = Self::convexity_factor(s.mean_reversion, s.swap_tenor1, s.annuity1, dt);
        let c2 = Self::convexity_factor(s.mean_reversion, s.swap_tenor2, s.annuity2, dt);
        (
            s.swap_rate1 + s.swap_rate1 * s.swap_rate1 * s.vol1 * s.vol1 * tf * c1 * 0.5,
            s.swap_rate2 + s.swap_rate2 * s.swap_rate2 * s.vol2 * s.vol2 * tf * c2 * 0.5,
        )
    }
    fn forward<T: Real>(s: &CmsSpreadSpec<T>) -> T {
        let (m1, m2) = Self::adjusted_rates(s);
        s.gearing1 * m1 - s.gearing2 * m2 + s.spread
    }
    fn normal_vol<T: Real>(s: &CmsSpreadSpec<T>) -> T {
        let t = s.fixing_time.max_val(T::zero());
        (s.gearing1 * s.gearing1 * s.vol1 * s.vol1 + s.gearing2 * s.gearing2 * s.vol2 * s.vol2
            - s.gearing1 * s.gearing2 * s.rho * s.vol1 * s.vol2 * 2.0)
            .max_val(T::zero())
            .sqrt()
            * t.sqrt()
    }
    pub fn coupon_rate_and_pv<T: Real>(s: &CmsSpreadSpec<T>) -> (T, T) {
        let x = Self::forward(s);
        (x, s.discount_df * s.accrual_factor * x)
    }
    pub fn caplet_pv<T: Real>(s: &CmsSpreadSpec<T>, strike: T) -> T {
        Self::lognormal_caplet_pv(s, strike)
    }
    pub fn floorlet_pv<T: Real>(s: &CmsSpreadSpec<T>, strike: T) -> T {
        // Put-call parity is exact for the same shifted-lognormal law and
        // avoids introducing a second, independently approximated quadrature.
        Self::caplet_pv(s, strike) - s.discount_df * s.accrual_factor * (Self::forward(s) - strike)
    }
    pub fn collar_pv<T: Real>(s: &CmsSpreadSpec<T>, cap_strike: T, floor_strike: T) -> T {
        Self::caplet_pv(s, cap_strike) - Self::floorlet_pv(s, floor_strike)
    }
    pub fn bachelier_spread_option_pv<T: Real>(s: &CmsSpreadSpec<T>, strike: T, is_cap: bool) -> T {
        s.discount_df
            * s.accrual_factor
            * bachelier_price(
                Self::forward(s),
                strike,
                Self::normal_vol(s),
                T::one(),
                T::one(),
                is_cap,
            )
    }
    pub fn lognormal_caplet_pv<T: Real>(s: &CmsSpreadSpec<T>, strike: T) -> T {
        // Brigo--Mercurio conditional shifted-lognormal quadrature.  The
        // 8-point Gauss--Hermite rule is used symmetrically (16 nodes); this
        // keeps all arithmetic generic in T for exact AD Greeks.
        const X: [f64; 8] = [
            -2.930_637_420_257_244,
            -1.981_656_756_695_843,
            -1.157_193_712_446_78,
            -0.381_186_990_207_322_1,
            0.381_186_990_207_322_1,
            1.157_193_712_446_78,
            1.981_656_756_695_843,
            2.930_637_420_257_244,
        ];
        const W: [f64; 8] = [
            0.000_199_604_072_211_367_6,
            0.017_077_983_007_413_48,
            0.207_802_325_814_891_9,
            0.661_147_012_558_241_3,
            0.661_147_012_558_241_3,
            0.207_802_325_814_891_9,
            0.017_077_983_007_413_48,
            0.000_199_604_072_211_367_6,
        ];
        let tf = s.fixing_time.max_val(T::zero());
        let sqt = tf.sqrt();
        let (m1, m2) = Self::adjusted_rates(s);
        let z1 = (m1 + s.shift1).max_val(T::from_f64(1e-12));
        let z2 = (m2 + s.shift2).max_val(T::from_f64(1e-12));
        let mut value = T::zero();
        for i in 0..X.len() {
            let z = T::from_f64(X[i] * 2.0_f64.sqrt());
            let e1 = (s.vol1 * sqt * z - s.vol1 * s.vol1 * tf * 0.5).exp();
            let e2 = (s.vol2 * sqt * (s.rho * z) - s.vol2 * s.vol2 * tf * 0.5).exp();
            let s1 = z1 * e1 - s.shift1;
            let s2 = z2 * e2 - s.shift2;
            let spread = s.gearing1 * s1 - s.gearing2 * s2 + s.spread;
            value += T::from_f64(W[i] / std::f64::consts::PI.sqrt())
                * (spread - strike).max_val(T::zero());
        }
        s.discount_df * s.accrual_factor * value
    }
}
