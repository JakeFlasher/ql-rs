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

#[derive(Copy, Clone, Debug)]
pub struct FreeBoundarySabrParams<T: Real = f64> {
    pub forward: T,
    pub expiry: T,
    pub alpha: T,
    pub beta: T,
    pub rho: T,
    pub nu: T,
}

#[derive(Copy, Clone, Debug)]
pub struct FreeBoundarySabrQuote<T: Real> {
    pub normal_vol: T,
    pub call_price: T,
    pub put_price: T,
    pub zero_strike_absorption_prob: T,
    pub risk_neutral_density: T,
}

fn signed_cev_pow<T: Real>(x: T, omb: T) -> T {
    let eps = T::from_f64(1e-8);
    let abs_x = (x * x + eps * eps).sqrt();
    let sgn = x / abs_x;
    sgn * (abs_x.ln() * omb).exp()
}

fn fb_sabr_normal_vol_and_call_undisc<T: Real>(
    p: &FreeBoundarySabrParams<T>,
    strike: T,
) -> (T, T, T) {
    let t = p.expiry.max_val(T::from_f64(1e-8));
    let omb = (T::one() - p.beta).max_val(T::from_f64(0.1));
    let eps = T::from_f64(1e-8);
    let f_abs = (p.forward * p.forward + eps * eps).sqrt();
    let k_abs = (strike * strike + eps * eps).sqrt();
    let mid_abs = ((f_abs + k_abs) * 0.5).max_val(T::from_f64(1e-6));
    let c_mid = (mid_abs.ln() * p.beta).exp();

    let df_cev = (signed_cev_pow(p.forward, omb) - signed_cev_pow(strike, omb)) / omb;
    let diff = p.forward - strike;
    // Exact CEV chord (F - K) / int_K^F |u|^{-beta} du with smooth ATM limit c_mid
    let chord_w = (diff * diff) / (diff * diff + T::from_f64(1e-8));
    let c_chord = c_mid * (T::one() - chord_w)
        + (diff / (df_cev + T::from_f64(1e-12)))
            .abs()
            .max_val(T::from_f64(1e-8))
            * chord_w;

    let y_cev = df_cev / p.alpha.max_val(T::from_f64(1e-8));
    let z = p.nu * y_cev;
    let z_over_chi = if z.val().abs() < 1e-4 {
        T::one() - p.rho * z * 0.5
            + (T::from_f64(2.0) - p.rho * p.rho * 3.0) * z * z / 12.0
    } else {
        let disc = (T::one() - p.rho * z * 2.0 + z * z)
            .max_val(T::from_f64(1e-12))
            .sqrt();
        let chi = ((disc + z - p.rho) / (T::one() - p.rho).max_val(T::from_f64(1e-8))).ln();
        z / chi
    };
    let time_corr = T::one()
        + t * (-(p.beta * (T::from_f64(2.0) - p.beta)) * p.alpha * p.alpha
            / (mid_abs * mid_abs / (c_mid * c_mid) * 24.0)
            + p.rho * p.beta * p.nu * p.alpha / (mid_abs / c_mid * 4.0)
            + (T::from_f64(2.0) - p.rho * p.rho * 3.0) * p.nu * p.nu / 24.0);
    let sigma_n = (p.alpha * c_chord * z_over_chi * time_corr).max_val(T::from_f64(1e-8));
    let std = sigma_n * t.sqrt();
    let d = diff / std;
    let call_undisc = diff * d.norm_cdf() + std * d.norm_pdf();
    (sigma_n, call_undisc, std)
}

pub fn evaluate_free_boundary_sabr<T: Real>(
    p: &FreeBoundarySabrParams<T>,
    strike: T,
    discount_factor: T,
) -> FreeBoundarySabrQuote<T> {
    let eps = T::from_f64(1e-8);
    let f_abs = (p.forward * p.forward + eps * eps).sqrt();
    let (sigma_n, call_undisc, std) = fb_sabr_normal_vol_and_call_undisc(p, strike);
    let put_undisc = call_undisc - (p.forward - strike);
    let z0 = -f_abs / std;
    let abs_prob = (z0.norm_cdf() * 2.0).max_val(T::zero()).min_val(T::one());

    // Breeden-Litzenberger / Hagan (2014) smile-consistent density d^2 C_undisc / dK^2
    let hk = 1e-4_f64;
    let hk_t = T::from_f64(hk);
    let (_, c_up, _) = fb_sabr_normal_vol_and_call_undisc(p, strike + hk_t);
    let (_, c_dn, _) = fb_sabr_normal_vol_and_call_undisc(p, strike - hk_t);
    let density = ((c_up - call_undisc * 2.0 + c_dn) / (hk * hk)).max_val(T::zero());

    FreeBoundarySabrQuote {
        normal_vol: sigma_n,
        call_price: discount_factor * call_undisc,
        put_price: discount_factor * put_undisc,
        zero_strike_absorption_prob: abs_prob,
        risk_neutral_density: density,
    }
}

#[cfg(test)]
mod phase19_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn free_boundary_sabr_negative_strikes_bachelier_parity_and_ad() {
        let p = FreeBoundarySabrParams {
            forward: -0.004,
            expiry: 1.0,
            alpha: 0.008,
            beta: 0.25,
            rho: -0.20,
            nu: 0.45,
        };
        for &k in &[-0.015, -0.004, 0.0, 0.010] {
            let q = evaluate_free_boundary_sabr(&p, k, 0.99);
            assert!(q.normal_vol > 0.0 && q.normal_vol.is_finite());
            assert!(q.risk_neutral_density > 0.0);
            assert!((q.call_price - q.put_price - 0.99 * (p.forward - k)).abs() < 1e-12);
        }

        let h = 1e-5;
        let eval_f = |f: f64| {
            let mut sp = p;
            sp.forward = f;
            evaluate_free_boundary_sabr(&sp, -0.002, 0.99).call_price
        };
        let fd = (eval_f(-0.004 + h) - eval_f(-0.004 - h)) / (2.0 * h);
        type D = Dual<1>;
        let pd = FreeBoundarySabrParams {
            forward: D::variable(-0.004, 0),
            expiry: D::one(),
            alpha: D::from_value(0.008),
            beta: D::from_value(0.25),
            rho: D::from_value(-0.20),
            nu: D::from_value(0.45),
        };
        let ad = evaluate_free_boundary_sabr(&pd, D::from_value(-0.002), D::from_value(0.99))
            .call_price;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 2e-6);
    }
}
