#![allow(clippy::excessive_precision)]
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
        T::one() - p.rho * z * 0.5 + (T::from_f64(2.0) - p.rho * p.rho * 3.0) * z * z / 12.0
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

#[derive(Copy, Clone, Debug)]
pub struct ThirdOrderSpeedVolExpansionSpec<T: Real = f64> {
    pub xi0: T,
    pub omega_vol_of_vol: T,
    pub k_mean_reversion: T,
    pub expiry: T,
    pub window_delta: T,
    pub strike: T,
    pub discount_factor: T,
}

#[derive(Copy, Clone, Debug)]
pub struct ThirdOrderSpeedVolExpansionQuote<T: Real = f64> {
    pub model_forward: T,
    pub call_price: T,
    pub put_price: T,
    pub implied_vol: T,
    pub bs_gamma: T,
    pub bs_speed: T,
    pub p0_call: T,
    pub p1_call: T,
    pub p2_call: T,
    pub p3_speed_call: T,
}

/// Bourgey, De Marco & Gobet (2022) Corollary 6 & Proposition 7 Third-Order Black-Scholes
/// Speed ($\mathrm{Speed}_{\mathrm{BS}}$) Expansion (`papers/bourgey_demarco_gobet_2022_bergomi_expansions.pdf`, pp. 9-10).
pub fn evaluate_third_order_speed_vol_expansion<T: Real>(
    spec: &ThirdOrderSpeedVolExpansionSpec<T>,
) -> ThirdOrderSpeedVolExpansionQuote<T> {
    let t = spec.expiry.max_val(T::from_f64(1e-4));
    let delta = spec.window_delta.max_val(T::from_f64(1e-4));
    let k_mr = spec.k_mean_reversion.max_val(T::from_f64(1e-4));
    let om = spec.omega_vol_of_vol.max_val(T::from_f64(1e-6));
    let xi0 = spec.xi0.max_val(T::from_f64(1e-8));
    let kappa = spec.strike.max_val(T::from_f64(1e-8));

    // Proposition 7 & Eq. (2.18): exact 1F Bergomi integrated kernel moments on [T, T + Delta]
    let kd = k_mr * delta;
    let one_m_e2kt = T::one() - (-k_mr * t * 2.0).exp();
    let om2 = om * om;
    let a_t = om2 * one_m_e2kt / (k_mr * 2.0);

    let m1 = (T::one() - (-kd).exp()) / kd;
    let m2 = (T::one() - (-kd * 2.0).exp()) / (kd * 2.0);
    let m3 = (T::one() - (-kd * 3.0).exp()) / (kd * 3.0);
    let m4 = (T::one() - (-kd * 4.0).exp()) / (kd * 4.0);

    let c_bar = (a_t * m1 * m1).max_val(T::from_f64(1e-12));
    let v_bar = a_t * m2;
    let mu_p = -v_bar * 0.5 + xi0.ln();
    let sig_p2 = c_bar;
    let sig_p = sig_p2.sqrt();

    let var_v = (a_t * a_t * (m4 - m2 * m2)).max_val(T::zero());
    let cov_cv = a_t * a_t * m1 * (m3 - m1 * m2);
    let var_c = (a_t * a_t * m1 * m1 * (m2 - m1 * m1)).max_val(T::zero());

    let gamma1 = var_v * 0.125 + (v_bar - c_bar) * 0.5;
    let gamma2 = -cov_cv * 0.5;
    let gamma3 = var_c * 0.5;

    // Corollary 6: x_P, s_P, BS Gamma, BS Speed, and basis functions P_0, P_1, P_2, P_3
    let x_p = (mu_p * 0.5 + sig_p2 * 0.125).exp();
    let s_p = (sig_p * 0.5).max_val(T::from_f64(1e-8));
    let log_xk = (x_p / kappa).ln();
    let d1 = log_xk / s_p + s_p * 0.5;
    let d2 = d1 - s_p;

    let p0_call = x_p * d1.norm_cdf() - kappa * d2.norm_cdf();
    let p1_call = x_p * d1.norm_cdf() * 0.5;
    let bs_gamma = d1.norm_pdf() / (x_p * s_p);
    let bs_speed = -(bs_gamma / x_p) * (log_xk / (s_p * s_p) + T::from_f64(1.5));
    let p2_call = p1_call * 0.5 + x_p * x_p * bs_gamma * 0.25;
    let p3_speed_call = -p1_call * 0.5 + p2_call * 1.5 + x_p * x_p * x_p * bs_speed * 0.125;

    let pf0 = x_p;
    let pf1 = x_p * 0.5;
    let pf2 = x_p * 0.25;
    let pf3 = x_p * 0.125;
    let model_forward =
        (pf0 + gamma1 * pf1 + gamma2 * pf2 + gamma3 * pf3).max_val(T::from_f64(1e-6));

    let raw_call = p0_call + gamma1 * p1_call + gamma2 * p2_call + gamma3 * p3_speed_call;
    let intrinsic = (model_forward - kappa).max_val(T::zero());
    let call_undisc = raw_call.max_val(intrinsic);
    let put_undisc = (call_undisc - (model_forward - kappa)).max_val(T::zero());
    let call_price = spec.discount_factor * call_undisc;
    let put_price = spec.discount_factor * put_undisc;

    // Invert call_undisc against Black-76 around model_forward via ift_lift_1d
    let f_val = model_forward.val();
    let k_val = kappa.val();
    let t_val = t.val();
    let target_val = call_undisc.val();
    let implied_vol = if target_val <= (f_val - k_val).max(0.0) + 1e-6 {
        s_p / t.sqrt()
    } else {
        let bs_f64 =
            |v: f64| -> f64 { crate::black_scholes_price(f_val, k_val, v, t_val, 0.0, 0.0, true) };
        let mut lo = 1e-4_f64;
        let mut hi = 5.0_f64;
        for _ in 0..36 {
            let mid = 0.5 * (lo + hi);
            if bs_f64(mid) < target_val {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let mut iv_root = 0.5 * (lo + hi);
        let h_iv = 1e-5;
        for _ in 0..4 {
            let r = bs_f64(iv_root) - target_val;
            let dr = (bs_f64(iv_root + h_iv) - bs_f64(iv_root - h_iv)) / (2.0 * h_iv);
            if dr.abs() > 1e-12 {
                iv_root = (iv_root - r / dr).clamp(1e-4, 5.0);
            }
        }
        let f_xx = (bs_f64(iv_root + h_iv) - 2.0 * bs_f64(iv_root) + bs_f64(iv_root - h_iv))
            / (h_iv * h_iv);
        let iv_t = T::from_f64(iv_root);
        let h_t = T::from_f64(h_iv);
        let c_up = crate::black_scholes_price(
            model_forward,
            kappa,
            iv_t + h_t,
            t,
            T::zero(),
            T::zero(),
            true,
        );
        let c_dn = crate::black_scholes_price(
            model_forward,
            kappa,
            iv_t - h_t,
            t,
            T::zero(),
            T::zero(),
            true,
        );
        let f_x = (c_up - c_dn) / (2.0 * h_iv);
        ql_ad::ift_lift_1d(
            iv_root,
            |v| {
                crate::black_scholes_price(model_forward, kappa, v, t, T::zero(), T::zero(), true)
                    - call_undisc
            },
            f_x,
            f_xx,
        )
    };

    ThirdOrderSpeedVolExpansionQuote {
        model_forward,
        call_price,
        put_price,
        implied_vol,
        bs_gamma,
        bs_speed,
        p0_call,
        p1_call,
        p2_call,
        p3_speed_call,
    }
}

#[cfg(test)]
mod domain_suite_32_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn third_order_speed_vol_expansion_put_call_parity_repricing_and_ad() {
        let spec = ThirdOrderSpeedVolExpansionSpec {
            xi0: 0.04,
            omega_vol_of_vol: 1.8,
            k_mean_reversion: 3.0,
            expiry: 0.25,
            window_delta: 1.0 / 12.0,
            strike: 0.20,
            discount_factor: 0.99,
        };
        let q = evaluate_third_order_speed_vol_expansion(&spec);
        assert!(q.bs_gamma > 0.0);
        assert!(q.bs_speed.is_finite());
        assert!((q.call_price - q.put_price - 0.99 * (q.model_forward - 0.20)).abs() < 1e-12);
        let repriced =
            crate::black_scholes_price(q.model_forward, 0.20, q.implied_vol, 0.25, 0.0, 0.0, true)
                * 0.99;
        assert!((repriced - q.call_price).abs() < 1e-7);

        let h = 1e-5;
        let eval_xi = |xi: f64| {
            let mut s = spec;
            s.xi0 = xi;
            evaluate_third_order_speed_vol_expansion(&s).call_price
        };
        let fd = (eval_xi(0.04 + h) - eval_xi(0.04 - h)) / (2.0 * h);
        type D = Dual<1>;
        let spec_d = ThirdOrderSpeedVolExpansionSpec {
            xi0: D::variable(0.04, 0),
            omega_vol_of_vol: D::from_value(1.8),
            k_mean_reversion: D::from_value(3.0),
            expiry: D::from_value(0.25),
            window_delta: D::from_value(1.0 / 12.0),
            strike: D::from_value(0.20),
            discount_factor: D::from_value(0.99),
        };
        let ad = evaluate_third_order_speed_vol_expansion(&spec_d).call_price;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_31_tests {
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
        let ad =
            evaluate_free_boundary_sabr(&pd, D::from_value(-0.002), D::from_value(0.99)).call_price;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 2e-6);
    }
}
