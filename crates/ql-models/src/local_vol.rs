use ql_ad::Real;
pub struct DupireLocalVol;
impl DupireLocalVol {
    pub fn dupire_local_vol<T: Real, F: Fn(T, T) -> T>(f: F, expiry: T, log_moneyness: T) -> T {
        let h = 1e-5;
        let w = f(expiry, log_moneyness);
        let wt = (f(expiry + T::from_f64(h), log_moneyness)
            - f(expiry - T::from_f64(h), log_moneyness))
            / (2. * h);
        let wy = (f(expiry, log_moneyness + T::from_f64(h))
            - f(expiry, log_moneyness - T::from_f64(h)))
            / (2. * h);
        let wyy = (f(expiry, log_moneyness + T::from_f64(h)) - T::from_f64(2.) * w
            + f(expiry, log_moneyness - T::from_f64(h)))
            / (h * h);
        (wt / (T::one() - log_moneyness / w * wy
            + T::from_f64(0.25)
                * (-T::from_f64(0.25) - T::one() / w + log_moneyness * log_moneyness / (w * w))
                * wy
                * wy
            + T::from_f64(0.5) * wyy))
            .sqrt()
    }
}
