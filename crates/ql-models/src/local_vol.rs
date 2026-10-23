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
        // The closure is evaluated on a second-order jet with lane zero equal
        // to expiry and lane one equal to log-moneyness.  Extracting these
        // coefficients is exact and preserves any incoming parameter jets;
        // finite differences here would both lose AD lanes and amplify roundoff.
        let mut wt = Dual2::from_value(w.d1[0]);
        wt.d1 = [w.d2[0][0], w.d2[0][1]];
        let mut wy = Dual2::from_value(w.d1[1]);
        wy.d1 = [w.d2[0][1], w.d2[1][1]];
        let wyy = Dual2::from_value(w.d2[1][1]);
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
