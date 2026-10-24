use ql_ad::Real;

/// Hagan's normal SABR approximation, generic over the AD scalar.
pub fn normal_sabr_vol<T: Real>(f: T, k: T, expiry: T, alpha: T, beta: T, rho: T, nu: T) -> T {
    let fk = (f * k).max_val(T::from_f64(1e-14));
    let omb = T::one() - beta;
    let mid = fk.pow(omb * 0.5);
    let z = nu / alpha.max_val(T::from_f64(1e-12)) * mid * (f - k);
    let ratio = if z.val().abs() < 1e-4 {
        T::one() - rho * z * 0.5
            + (T::from_f64(2.0) - T::from_f64(3.0) * rho * rho) * z * z / 12.0
            + rho * (T::from_f64(5.0) * rho * rho - T::from_f64(3.0)) * z * z * z / 24.0
    } else {
        let d = (T::one() - T::from_f64(2.0) * rho * z + z * z).sqrt();
        z / ((d + z - rho) / (T::one() - rho)).ln()
    };
    alpha
        * mid
        * ratio
        * (T::one()
            + expiry
                * (-(beta * (T::from_f64(2.0) - beta)) * alpha * alpha / (fk.pow(omb) * 24.0)
                    + rho * beta * nu * alpha / (mid * 4.0)
                    + (T::from_f64(2.0) - T::from_f64(3.0) * rho * rho) * nu * nu / 24.0))
}
