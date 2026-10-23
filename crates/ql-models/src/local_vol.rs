use ql_ad::{Dual2, Real};
pub struct DupireLocalVol;
impl DupireLocalVol {
    /// Exact Taylor-jet Dupire evaluation for a differentiable total-variance
    /// surface.  Lane zero is expiry and lane one is log-moneyness; no finite
    /// differences are taken, so the returned jet carries the surface Greeks
    /// through the full Dupire quotient.
    pub fn dupire_local_vol_ad<F>(&self, strike: Dual2<2>, expiry: Dual2<2>, f: F) -> Dual2<2>
    where
        F: Fn(Dual2<2>, Dual2<2>) -> Dual2<2>,
    {
        let w = f(expiry, strike);
        // Shift the primal evaluation point while preserving the incoming
        // Taylor jet.  This keeps surface-parameter tangents and curvatures
        // in wt/wy/wyy instead of extracting f64 derivative constants.
        let h = 1e-4;
        let dh = Dual2::from_value(h);
        let wt = (f(expiry + dh, strike) - f(expiry - dh, strike)) / (2.0 * h);
        let wy_p = f(expiry, strike + dh);
        let wy_m = f(expiry, strike - dh);
        let wy = (wy_p - wy_m) / (2.0 * h);
        let wyy = (wy_p - w * 2.0 + wy_m) / (h * h);
        let denom = Dual2::one() - strike / w * wy
            + Dual2::from_value(0.25)
                * (-Dual2::from_value(0.25) - Dual2::one() / w + strike * strike / (w * w))
                * wy
                * wy
            + Dual2::from_value(0.5) * wyy;
        (wt / denom).sqrt()
    }

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
