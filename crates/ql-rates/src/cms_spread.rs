use ql_ad::Real;
use ql_vol::{bachelier_price, black76_price};

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
    fn forward<T: Real>(s: &CmsSpreadSpec<T>) -> T {
        s.gearing1 * s.swap_rate1 - s.gearing2 * s.swap_rate2 + s.spread
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
        s.discount_df
            * s.accrual_factor
            * bachelier_price(
                Self::forward(s),
                strike,
                Self::normal_vol(s),
                T::one(),
                T::one(),
                true,
            )
    }
    pub fn floorlet_pv<T: Real>(s: &CmsSpreadSpec<T>, strike: T) -> T {
        s.discount_df
            * s.accrual_factor
            * bachelier_price(
                Self::forward(s),
                strike,
                Self::normal_vol(s),
                T::one(),
                T::one(),
                false,
            )
    }
    pub fn collar_pv<T: Real>(s: &CmsSpreadSpec<T>, cap_strike: T, floor_strike: T) -> T {
        Self::floorlet_pv(s, floor_strike) - Self::caplet_pv(s, cap_strike)
    }
    pub fn bachelier_spread_option_pv<T: Real>(s: &CmsSpreadSpec<T>, strike: T, is_cap: bool) -> T {
        if is_cap {
            Self::caplet_pv(s, strike)
        } else {
            Self::floorlet_pv(s, strike)
        }
    }
    pub fn lognormal_caplet_pv<T: Real>(s: &CmsSpreadSpec<T>, strike: T) -> T {
        let x = Self::forward(s);
        let v = Self::normal_vol(s) / x.abs().max_val(T::from_f64(1e-12));
        s.discount_df * s.accrual_factor * black76_price(x, strike, v, T::one(), T::one(), true)
    }
}
