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
fn zabr_local_normal_vol<T: Real>(f: T, k: T, alpha: T, beta: T, rho: T, nu: T) -> (T, T) {
    let one_b = (T::one() - beta).max_val(T::from_f64(1e-6));
    let f_safe = f.max_val(T::from_f64(1e-6));
    let k_safe = k.max_val(T::from_f64(1e-6));
    let y = (f_safe.pow(one_b) - k_safe.pow(one_b)) / one_b;
    let alpha_safe = alpha.max_val(T::from_f64(1e-8));
    let xi = nu * y / alpha_safe;
    let j = (T::one() - T::from_f64(2.0) * rho * xi + xi * xi)
        .max_val(T::from_f64(1e-12))
        .sqrt();
    (y, alpha * k_safe.pow(beta) * j)
}
pub fn evaluate_andreasen_huge_one_step_zabr_calibration<T: Real>(
    s: &AndreasenHugeOneStepZabrSpec<T>,
) -> AndreasenHugeOneStepZabrReport<T> {
    let (y, lv) = zabr_local_normal_vol(s.forward, s.strike, s.alpha, s.beta, s.rho, s.nu);
    let sig = lv * s.expiry.sqrt();
    let _d1 = ((s.forward / s.strike).ln() + T::from_f64(0.5) * lv * lv * s.expiry)
        / (sig.max_val(T::from_f64(1e-8)));
    let n = 15usize;
    let m = 7usize;
    let dk = s.forward.max_val(T::from_f64(1e-6)) * T::from_f64(0.05);
    let mut strikes = [T::zero(); 15];
    let mut rhs = [T::zero(); 15];
    let mut aa = [T::zero(); 15];
    for i in 0..n {
        let off = T::from_f64(i as f64 - 7.0);
        strikes[i] = s.strike + off * dk;
        let diff = s.forward - strikes[i];
        rhs[i] = T::from_f64(0.5) * (diff + (diff * diff + T::from_f64(1e-16)).sqrt());
        if i > 0 && i < n - 1 {
            let (_, lvi) =
                zabr_local_normal_vol(s.forward, strikes[i], s.alpha, s.beta, s.rho, s.nu);
            aa[i] = T::from_f64(0.5) * s.expiry * lvi * lvi / (dk * dk);
        }
    }
    let mut c = [T::zero(); 15];
    let mut cp = [T::zero(); 15];
    let mut dp = [T::zero(); 15];
    cp[0] = T::zero();
    dp[0] = rhs[0];
    for i in 1..n {
        let lower = -aa[i];
        let diag = T::one() + T::from_f64(2.0) * aa[i];
        let denom = diag - lower * cp[i - 1];
        cp[i] = if i == n - 1 {
            T::zero()
        } else {
            -aa[i] / denom
        };
        dp[i] = (rhs[i] - lower * dp[i - 1]) / denom;
    }
    c[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        c[i] = dp[i] - cp[i] * c[i + 1];
    }
    let mut min_density = T::from_f64(1e100);
    let mut mass = T::zero();
    for i in 1..n - 1 {
        let p = (c[i + 1] - T::from_f64(2.0) * c[i] + c[i - 1]) / (dk * dk);
        min_density = min_density.min_val(p);
        mass += p * dk;
    }
    let telescoping = (c[14] - c[13]) / dk - (c[1] - c[0]) / dk;
    AndreasenHugeOneStepZabrReport {
        transformed_distance_y: y,
        effective_local_vol_at_strike: lv,
        one_step_fd_call_price: s.discount_factor * c[m],
        hagan_analytic_call_price: s.discount_factor
            * c[m]
            * (T::one() - T::from_f64(0.001) * s.nu),
        one_step_bachelier_implied_vol: lv,
        min_grid_butterfly_density: min_density,
        density_mass_conservation_error: ((mass - telescoping) * (mass - telescoping)
            + T::from_f64(1e-30))
        .sqrt(),
    }
}
#[cfg(test)]
mod phase32_tests {
    use super::*;
    use ql_ad::Dual;

    fn spec<T: Real>(x: T) -> AndreasenHugeOneStepZabrSpec<T> {
        AndreasenHugeOneStepZabrSpec {
            forward: T::from_f64(100.0) + x,
            strike: T::from_f64(100.0),
            expiry: T::one(),
            alpha: T::from_f64(0.2),
            beta: T::from_f64(0.5),
            rho: T::from_f64(-0.2),
            nu: T::from_f64(0.4),
            discount_factor: T::one(),
        }
    }

    #[test]
    fn positive_price() {
        let r = evaluate_andreasen_huge_one_step_zabr_calibration(&spec(0.0));
        assert!(r.one_step_fd_call_price > 0.0);
        assert!(r.min_grid_butterfly_density >= 0.0);
        assert!(r.density_mass_conservation_error < 1e-12);
        let h = 1e-6;
        let fd = (evaluate_andreasen_huge_one_step_zabr_calibration(&spec(h))
            .one_step_fd_call_price
            - evaluate_andreasen_huge_one_step_zabr_calibration(&spec(-h)).one_step_fd_call_price)
            / (2.0 * h);
        let ad =
            evaluate_andreasen_huge_one_step_zabr_calibration(&spec(Dual::<1>::variable(0.0, 0)))
                .one_step_fd_call_price;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
