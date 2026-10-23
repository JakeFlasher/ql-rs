use ql_ad::Real;
use ql_vol::black_scholes_price;

pub struct ForwardStartOptionEngine;
impl ForwardStartOptionEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        spot: T,
        moneyness_alpha: T,
        rate: T,
        div: T,
        vol: T,
        reset_time: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        let tau = (expiry - reset_time).max_val(T::zero());
        let f = spot * ((rate - div) * reset_time).exp();
        black_scholes_price(f, moneyness_alpha * f, vol, tau, rate, div, is_call)
            * (-rate * reset_time).exp()
    }
}
#[derive(Clone, Debug)]
pub struct CliquetOptionSpec<T: Real = f64> {
    pub spot: T,
    pub rate: T,
    pub div: T,
    pub vol: T,
    pub reset_times: Vec<f64>,
    pub local_floor: T,
    pub local_cap: T,
    pub notional: T,
}
pub struct AnalyticCliquetEngine;
impl AnalyticCliquetEngine {
    pub fn price<T: Real>(s: &CliquetOptionSpec<T>, is_call: bool) -> T {
        let mut v = T::zero();
        for w in s.reset_times.windows(2) {
            let dt = T::from_f64(w[1] - w[0]);
            let f = s.spot * ((s.rate - s.div) * T::from_f64(w[0])).exp();
            let c = black_scholes_price(f, f, s.vol, dt, s.rate, s.div, true);
            let clipped = c.max_val(s.local_floor).min_val(s.local_cap);
            v += clipped * (-s.rate * T::from_f64(w[1])).exp();
        }
        let _ = is_call;
        s.notional * v
    }
}
pub struct WriterExtensibleOptionEngine;
impl WriterExtensibleOptionEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        spot: T,
        k1: T,
        k2: T,
        rate: T,
        div: T,
        vol: T,
        t1: T,
        t2: T,
        is_call: bool,
    ) -> T {
        black_scholes_price(spot, k1, vol, t1, rate, div, is_call)
            + black_scholes_price(spot, k2, vol, t2, rate, div, is_call)
    }
}
