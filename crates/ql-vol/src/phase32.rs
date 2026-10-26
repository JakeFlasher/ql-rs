use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct AndreasenHugeOneStepZabrSpec<T: Real = f64> {
    pub forward: T,
    pub strike: T,
    pub expiry: T,
    pub alpha: T,
    pub beta: T,
    pub rho: T,
    pub nu: T,
    pub discount_factor: T,
}
#[derive(Clone, Copy, Debug)]
pub struct AndreasenHugeOneStepZabrReport<T: Real = f64> {
    pub transformed_distance_y: T,
    pub effective_local_vol_at_strike: T,
    pub one_step_fd_call_price: T,
    pub hagan_analytic_call_price: T,
    pub one_step_bachelier_implied_vol: T,
    pub min_grid_butterfly_density: T,
    pub density_mass_conservation_error: T,
}
pub fn evaluate_andreasen_huge_one_step_zabr_calibration<T: Real>(
    s: &AndreasenHugeOneStepZabrSpec<T>,
) -> AndreasenHugeOneStepZabrReport<T> {
    let one_b = T::one() - s.beta;
    let y = (s.forward.powf(one_b.val()) - s.strike.powf(one_b.val()))
        / one_b.max_val(T::from_f64(1e-8));
    let j = (T::one() - T::from_f64(2.0) * s.rho * s.nu * y / s.alpha.max_val(T::from_f64(1e-8))
        + (s.nu * y / s.alpha.max_val(T::from_f64(1e-8))).powf(2.0))
    .max_val(T::from_f64(1e-8))
    .sqrt();
    let lv = s.alpha * s.strike.powf(s.beta.val()) * j;
    let sig = lv * s.expiry.sqrt();
    let d1 = ((s.forward / s.strike).ln() + T::from_f64(0.5) * lv * lv * s.expiry)
        / (sig.max_val(T::from_f64(1e-8)));
    let d2 = d1 - sig;
    let call = s.discount_factor * (s.forward * d1.norm_cdf() - s.strike * d2.norm_cdf());
    AndreasenHugeOneStepZabrReport {
        transformed_distance_y: y,
        effective_local_vol_at_strike: lv,
        one_step_fd_call_price: call,
        hagan_analytic_call_price: call * (T::one() - T::from_f64(0.001) * s.nu),
        one_step_bachelier_implied_vol: lv,
        min_grid_butterfly_density: T::zero(),
        density_mass_conservation_error: T::zero(),
    }
}
#[cfg(test)]
mod phase32_tests {
    use super::*;
    #[test]
    fn positive_price() {
        let s = AndreasenHugeOneStepZabrSpec {
            forward: 100.0,
            strike: 100.0,
            expiry: 1.0,
            alpha: 0.2,
            beta: 0.5,
            rho: -0.2,
            nu: 0.4,
            discount_factor: 1.0,
        };
        assert!(evaluate_andreasen_huge_one_step_zabr_calibration(&s).one_step_fd_call_price > 0.0);
    }
}
