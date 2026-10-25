#![allow(clippy::excessive_precision)]
use crate::BatesParams;
use ql_ad::Real;
pub fn replicate_cboe_vix_index<T: Real>(expiry: T, rate: T, forward: T, quotes: &[(T, T)]) -> T {
    // CBOE's strip replication requires a positive expiry, forward and strike
    // grid.  Returning a zero quote for an unusable strip is preferable to
    // propagating NaNs through a larger calibration graph.
    if quotes.len() < 2
        || !expiry.val().is_finite()
        || expiry.val() <= 0.0
        || !forward.val().is_finite()
        || forward.val() <= 0.0
        || !rate.val().is_finite()
    {
        return T::zero();
    }
    let mut q = quotes.to_vec();
    q.sort_by(|a, b| a.0.val().total_cmp(&b.0.val()));
    if q.iter()
        .any(|(k, p)| !k.val().is_finite() || k.val() <= 0.0 || !p.val().is_finite())
    {
        return T::zero();
    }
    // Duplicate strikes do not define a finite integration cell.  Keep the
    // first occurrence to avoid double-counting their option quote.
    q.dedup_by(|a, b| (a.0.val() - b.0.val()).abs() <= 1e-14);
    if q.len() < 2 {
        return T::zero();
    }
    let mut sum = T::zero();
    for i in 0..q.len() {
        let k = q[i].0;
        let dk = if i == 0 {
            q[1].0 - k
        } else if i + 1 == q.len() {
            k - q[i - 1].0
        } else {
            (q[i + 1].0 - q[i - 1].0) * 0.5
        };
        if dk.val() <= 0.0 {
            return T::zero();
        }
        sum += dk / (k * k) * q[i].1.max_val(T::zero());
    }
    let k0 = q
        .iter()
        .filter(|(k, _)| k.val() <= forward.val())
        .map(|(k, _)| *k)
        .last()
        .unwrap_or(q[0].0);
    let var = (sum * T::from_f64(2.0) / expiry) * (rate * expiry).exp()
        - (forward / k0 - T::one()).powf(2.0) / expiry;
    T::from_f64(100.0) * var.max_val(T::zero()).sqrt()
}
#[derive(Copy, Clone, Debug)]
pub struct VixOptionValuation<T: Real> {
    pub vix_future: T,
    pub call_price: T,
    pub put_price: T,
}
pub fn heston_bates_vix_future_and_option_price<T: Real>(
    p: &BatesParams<T>,
    expiry: T,
    vix_window: f64,
    strike: T,
    discount_factor: T,
) -> VixOptionValuation<T> {
    if !expiry.val().is_finite()
        || expiry.val() < 0.0
        || !vix_window.is_finite()
        || vix_window <= 0.0
        || !strike.val().is_finite()
        || strike.val() <= 0.0
        || !discount_factor.val().is_finite()
        || discount_factor.val() < 0.0
    {
        return VixOptionValuation {
            vix_future: T::zero(),
            call_price: T::zero(),
            put_price: T::zero(),
        };
    }
    let k = p.heston.kappa;
    let kw = k * T::from_f64(vix_window);
    let a = if kw.val().abs() < 1e-6 {
        T::one() - kw * 0.5 + kw * kw / 6.0
    } else {
        (T::one() - (-kw).exp()) / kw
    };
    let kappa_j = (p.mu_j + p.sigma_j * p.sigma_j * 0.5).exp() - T::one();
    let jump = p.lambda * (kappa_j - p.mu_j);
    let b = p.heston.theta * (T::one() - a) + T::from_f64(2.0) * jump;
    let mean =
        p.heston.v0 * (-k * expiry).exp() + p.heston.theta * (T::one() - (-k * expiry).exp());
    let mu_var = (a * mean + b).max_val(T::from_f64(1e-12));
    let t_exp = expiry.max_val(T::from_f64(1e-8));
    let e_k = (-k * t_exp).exp();
    let kt = k * t_exp;
    let var_v = if kt.val().abs() < 1e-6 {
        p.heston.sigma_v
            * p.heston.sigma_v
            * t_exp
            * (p.heston.v0 + (p.heston.theta * 0.5 - p.heston.v0) * kt).max_val(T::zero())
    } else {
        let s2_k = p.heston.sigma_v * p.heston.sigma_v / k;
        (p.heston.v0 * s2_k * (e_k - e_k * e_k)
            + p.heston.theta * s2_k * 0.5 * (T::one() - e_k) * (T::one() - e_k))
            .max_val(T::zero())
    };
    let rel_var = (a * a * var_v) / (mu_var * mu_var);
    let total_vix_var = (T::one() + rel_var).ln() * 0.25;
    let vvol = (total_vix_var / t_exp).max_val(T::from_f64(1e-12)).sqrt();
    let f = mu_var.sqrt() * (-total_vix_var * 0.5).exp();
    if expiry.val() <= 0.0 {
        return VixOptionValuation {
            vix_future: f,
            call_price: (f - strike).max_val(T::zero()) * discount_factor,
            put_price: (strike - f).max_val(T::zero()) * discount_factor,
        };
    }
    let c = ql_vol::black_scholes_price(f, strike, vvol, t_exp, T::zero(), T::zero(), true)
        * discount_factor;
    let put = ql_vol::black_scholes_price(f, strike, vvol, t_exp, T::zero(), T::zero(), false)
        * discount_factor;
    VixOptionValuation {
        vix_future: f,
        call_price: c,
        put_price: put,
    }
}

#[derive(Copy, Clone, Debug)]
pub struct QuinticOuParams<T: Real = f64> {
    pub xi0: T,
    pub kappa: T,
    pub nu: T,
    pub rho: T,
    pub alpha0: T,
    pub alpha1: T,
    pub alpha3: T,
    pub alpha5: T,
}

#[derive(Copy, Clone, Debug)]
pub struct QuinticOuVixSpxReport<T: Real> {
    pub vix_future: T,
    pub vix_call_price: T,
    pub vix_put_price: T,
    pub vix_implied_vol: T,
    pub spx_atm_vol: T,
    pub spx_atm_skew: T,
}

pub fn evaluate_quintic_ou_spx_vix<T: Real>(
    p: &QuinticOuParams<T>,
    expiry: T,
    vix_window: f64,
    vix_strike: T,
    discount_factor: T,
) -> QuinticOuVixSpxReport<T> {
    // 8-point Gauss-Hermite nodes for standard normal Z ~ N(0, 1)
    const GH_Z: [f64; 8] = [
        -4.14454718612589,
        -2.80248586128754,
        -1.63651904243511,
        -0.53907981135138,
        0.53907981135138,
        1.63651904243511,
        2.80248586128754,
        4.14454718612589,
    ];
    const GH_W: [f64; 8] = [
        0.00011261453835,
        0.00963522012073,
        0.11723990766013,
        0.37301225768079,
        0.37301225768079,
        0.11723990766013,
        0.00963522012073,
        0.00011261453835,
    ];

    let t_exp = expiry.max_val(T::from_f64(1e-6));
    let delta = vix_window.max(1e-4);
    let poly_sq_exp = |m_z: T, s_z: T| -> T {
        let mut s = T::zero();
        for (&z_i, &w_i) in GH_Z.iter().zip(GH_W.iter()) {
            let x = m_z + s_z * z_i;
            let x2 = x * x;
            let p_x = p.alpha0 + p.alpha1 * x + p.alpha3 * x * x2 + p.alpha5 * x * x2 * x2;
            s += p_x * p_x * w_i;
        }
        s.max_val(T::from_f64(1e-12))
    };

    // Unconditional normalization E[P(Z)^2] so E[sigma_u^2] == xi0 identically
    let norm_denom = poly_sq_exp(T::zero(), T::one());

    // 4-point Gauss-Legendre over u in [T, T + delta] for conditional VIX_T^2(Z)
    const GL_X: [f64; 4] = [
        -0.8611363115940526,
        -0.3399810435848563,
        0.3399810435848563,
        0.8611363115940526,
    ];
    const GL_W: [f64; 4] = [
        0.3478548451374538,
        0.6521451548625461,
        0.6521451548625461,
        0.3478548451374538,
    ];

    let var_x_t = p.nu * p.nu * (T::one() - (-p.kappa * t_exp * 2.0).exp())
        / (p.kappa * 2.0).max_val(T::from_f64(1e-8));
    let mut vix_fut = T::zero();
    let mut vix_call_undisc = T::zero();
    let eps = T::from_f64(1e-4);

    for (&z_node, &w_node) in GH_Z.iter().zip(GH_W.iter()) {
        let mut vix2_z = T::zero();
        for (&gx, &gw) in GL_X.iter().zip(GL_W.iter()) {
            let tau = 0.5 * delta * (gx + 1.0);
            let var_inc = p.nu * p.nu * (T::one() - (-p.kappa * (2.0 * tau)).exp())
                / (p.kappa * 2.0).max_val(T::from_f64(1e-8));
            let var_u = var_x_t * (-p.kappa * (2.0 * tau)).exp() + var_inc;
            let rho_u = ((var_x_t * (-p.kappa * (2.0 * tau)).exp())
                / var_u.max_val(T::from_f64(1e-12)))
            .max_val(T::zero())
            .min_val(T::from_f64(0.999999))
            .sqrt();
            let s_cond = (T::one() - rho_u * rho_u).max_val(T::from_f64(1e-12)).sqrt();
            let m_cond = rho_u * z_node;
            let cond_exp = p.xi0 * poly_sq_exp(m_cond, s_cond) / norm_denom;
            vix2_z += cond_exp * (0.5 * gw);
        }
        let vix_z = T::from_f64(100.0) * vix2_z.max_val(T::from_f64(1e-12)).sqrt();
        vix_fut += vix_z * w_node;
        let diff = vix_z - vix_strike;
        let call_payoff = (diff + (diff * diff + eps * eps).sqrt()) * 0.5;
        vix_call_undisc += call_payoff * w_node;
    }

    let vix_call_price = discount_factor * vix_call_undisc;
    let vix_put_price = vix_call_price - discount_factor * (vix_fut - vix_strike);
    let atm_poly_deriv = p.alpha1 + p.alpha3 * 3.0 + p.alpha5 * 15.0;
    let eff_vvol = (p.nu * atm_poly_deriv.abs() / norm_denom.sqrt()).max_val(T::from_f64(0.05));

    // Invert undiscounted VIX call price via Black-76 + 1D IFT lift
    let f_val = vix_fut.val().max(1e-8);
    let k_val = vix_strike.val().max(1e-8);
    let t_val = t_exp.val().max(1e-8);
    let intrinsic_val = (f_val - k_val).max(0.0);
    let target_val = vix_call_undisc.val();
    let vix_iv = if target_val <= intrinsic_val + 1e-5 {
        eff_vvol
    } else {
        let bs_undisc_f64 = |vol: f64| -> f64 {
            ql_vol::black_scholes_price(f_val, k_val, vol, t_val, 0.0, 0.0, true)
        };
        let mut lo = 1e-4_f64;
        let mut hi = 5.0_f64;
        for _ in 0..36 {
            let mid = 0.5 * (lo + hi);
            if bs_undisc_f64(mid) < target_val {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let mut iv_root = 0.5 * (lo + hi);
        let h_iv = 1e-5;
        for _ in 0..4 {
            let r = bs_undisc_f64(iv_root) - target_val;
            let dr = (bs_undisc_f64(iv_root + h_iv) - bs_undisc_f64(iv_root - h_iv)) / (2.0 * h_iv);
            if dr.abs() > 1e-12 {
                iv_root = (iv_root - r / dr).clamp(1e-4, 5.0);
            }
        }
        let f_xx = (bs_undisc_f64(iv_root + h_iv) - 2.0 * bs_undisc_f64(iv_root) + bs_undisc_f64(iv_root - h_iv))
            / (h_iv * h_iv);
        let iv_t = T::from_f64(iv_root);
        let h_t = T::from_f64(h_iv);
        let c_up = ql_vol::black_scholes_price(vix_fut, vix_strike, iv_t + h_t, t_exp, T::zero(), T::zero(), true);
        let c_dn = ql_vol::black_scholes_price(vix_fut, vix_strike, iv_t - h_t, t_exp, T::zero(), T::zero(), true);
        let f_x = (c_up - c_dn) / (2.0 * h_iv);
        ql_ad::ift_lift_1d(
            iv_root,
            |s| ql_vol::black_scholes_price(vix_fut, vix_strike, s, t_exp, T::zero(), T::zero(), true) - vix_call_undisc,
            f_x,
            f_xx,
        )
    };
    let spx_atm_vol = p.xi0.max_val(T::from_f64(1e-12)).sqrt();
    let spx_atm_skew = p.rho * eff_vvol
        * (T::one() - (-p.kappa * t_exp).exp())
        / (p.kappa * t_exp).max_val(T::from_f64(1e-6));

    QuinticOuVixSpxReport {
        vix_future: vix_fut,
        vix_call_price,
        vix_put_price,
        vix_implied_vol: vix_iv,
        spx_atm_vol,
        spx_atm_skew,
    }
}

#[cfg(test)]
mod phase19_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn quintic_ou_flat_collapse_vix_parity_and_ad() {
        let flat = QuinticOuParams {
            xi0: 0.04,
            kappa: 2.0,
            nu: 0.8,
            rho: -0.7,
            alpha0: 1.0,
            alpha1: 0.0,
            alpha3: 0.0,
            alpha5: 0.0,
        };
        let rep_flat = evaluate_quintic_ou_spx_vix(&flat, 0.5, 1.0 / 12.0, 18.0, 0.985);
        // When alpha1 = alpha3 = alpha5 = 0, VIX_T is deterministically 100 * sqrt(0.04) = 20.0
        assert!((rep_flat.vix_future - 20.0).abs() < 1e-6);

        let p = QuinticOuParams {
            xi0: 0.04,
            kappa: 2.5,
            nu: 0.9,
            rho: -0.75,
            alpha0: 1.0,
            alpha1: 0.35,
            alpha3: 0.12,
            alpha5: 0.03,
        };
        let rep = evaluate_quintic_ou_spx_vix(&p, 0.5, 1.0 / 12.0, 20.0, 0.985);
        assert!(
            (rep.vix_call_price - rep.vix_put_price - 0.985 * (rep.vix_future - 20.0)).abs()
                < 1e-12
        );
        assert!((ql_vol::black_scholes_price(rep.vix_future, 20.0, rep.vix_implied_vol, 0.5, 0.0, 0.0, true) * 0.985 - rep.vix_call_price).abs() < 1e-6);
        assert!(rep.spx_atm_skew < 0.0);

        let h = 1e-5;
        let eval_xi = |xi: f64| {
            let mut q = p;
            q.xi0 = xi;
            evaluate_quintic_ou_spx_vix(&q, 0.5, 1.0 / 12.0, 20.0, 0.985).vix_call_price
        };
        let fd = (eval_xi(0.04 + h) - eval_xi(0.04 - h)) / (2.0 * h);
        type D = Dual<1>;
        let pd = QuinticOuParams {
            xi0: D::variable(0.04, 0),
            kappa: D::from_value(2.5),
            nu: D::from_value(0.9),
            rho: D::from_value(-0.75),
            alpha0: D::one(),
            alpha1: D::from_value(0.35),
            alpha3: D::from_value(0.12),
            alpha5: D::from_value(0.03),
        };
        let ad = evaluate_quintic_ou_spx_vix(
            &pd,
            D::from_value(0.5),
            1.0 / 12.0,
            D::from_value(20.0),
            D::from_value(0.985),
        )
        .vix_call_price;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}


#[derive(Copy, Clone, Debug)]
pub struct MixedRoughBergomiVixSpec<T: Real = f64> {
    pub xi0: T,
    pub hurst: T,
    pub mixing_weight_lambda: T,
    pub nu1: T,
    pub nu2: T,
    pub expiry: T,
    pub vix_window: f64,
    pub strike: T,
    pub discount_factor: T,
}

#[derive(Copy, Clone, Debug)]
pub struct MixedRoughBergomiVixReport<T: Real = f64> {
    pub vix_future: T,
    pub vix_call_price: T,
    pub vix_put_price: T,
    pub vix_implied_vol: T,
    pub proxy_sigma1: T,
    pub proxy_sigma2: T,
    pub vix_atm_skew_proxy: T,
}

/// Bourgey, De Marco & Gobet (2022) Theorem 3.1 & Proposition 9 Mixed Rough Bergomi
/// VIX Futures & Option Weak-Mixing Expansion (`papers/bourgey_demarco_gobet_2022_bergomi_expansions.pdf`).
pub fn evaluate_mixed_rough_bergomi_vix<T: Real>(
    spec: &MixedRoughBergomiVixSpec<T>,
) -> MixedRoughBergomiVixReport<T> {
    const GL4_X: [f64; 4] = [
        -0.861_136_311_594_052_6,
        -0.339_981_043_584_856_3,
        0.339_981_043_584_856_3,
        0.861_136_311_594_052_6,
    ];
    const GL4_W: [f64; 4] = [
        0.347_854_845_137_453_8,
        0.652_145_154_862_546_1,
        0.652_145_154_862_546_1,
        0.347_854_845_137_453_8,
    ];
    const GH9_X: [f64; 9] = [
        -3.190_993_201_781_528,
        -2.266_580_584_531_843,
        -1.468_553_289_216_668,
        -0.723_551_018_752_837_6,
        0.0,
        0.723_551_018_752_837_6,
        1.468_553_289_216_668,
        2.266_580_584_531_843,
        3.190_993_201_781_528,
    ];
    const GH9_W: [f64; 9] = [
        3.960_697_726_326_438e-5,
        0.004_943_624_275_536_947,
        0.088_474_527_394_376_57,
        0.432_651_559_002_555_8,
        0.720_235_215_606_051,
        0.432_651_559_002_555_8,
        0.088_474_527_394_376_57,
        0.004_943_624_275_536_947,
        3.960_697_726_326_438e-5,
    ];

    let t_exp = spec.expiry.max_val(T::from_f64(1e-4));
    let delta = spec.vix_window.max(1e-4);
    let h = spec
        .hurst
        .max_val(T::from_f64(0.02))
        .min_val(T::from_f64(0.49));
    let two_h = h * 2.0;
    let lam = spec
        .mixing_weight_lambda
        .max_val(T::zero())
        .min_val(T::one());

    // Fractional Riemann-Liouville variance integral int_0^T K_H(u, s)^2 ds = u^{2H} - (u - T)^{2H}
    let mut avg_k2 = T::zero();
    let mut avg_k_cov = T::zero();
    for i in 0..4 {
        let tau_i = 0.5 * delta * (1.0 + GL4_X[i]);
        let wi = 0.5 * GL4_W[i];
        let u_i = t_exp + tau_i;
        let k2_i = u_i.pow(two_h) - T::from_f64(tau_i).pow(two_h);
        avg_k2 += k2_i * wi;
        for j in 0..4 {
            let tau_j = 0.5 * delta * (1.0 + GL4_X[j]);
            let wj = 0.5 * GL4_W[j];
            let u_j = t_exp + tau_j;
            let k2_j = u_j.pow(two_h) - T::from_f64(tau_j).pow(two_h);
            // High-correlation fractional kernel overlap across the 30-day VIX window
            let dt_ij = (tau_i - tau_j).abs();
            let corr_ij = (T::one() - T::from_f64(0.15 * (dt_ij / delta)).pow(two_h))
                .max_val(T::from_f64(0.7));
            avg_k_cov += (k2_i * k2_j).sqrt() * corr_ij * (wi * wj);
        }
    }

    let sigma1 = (spec.nu1 * avg_k_cov.max_val(T::from_f64(1e-14)).sqrt()).max_val(T::from_f64(1e-6));
    let sigma2 = (spec.nu2 * avg_k_cov.max_val(T::from_f64(1e-14)).sqrt()).max_val(T::from_f64(1e-6));
    let mu1 = -spec.nu1 * spec.nu1 * avg_k2 * 0.5;
    let mu2 = -spec.nu2 * spec.nu2 * avg_k2 * 0.5;

    let res_var = (avg_k2 - avg_k_cov).max_val(T::zero());
    let gamma1_1 = -spec.nu1 * spec.nu1 * res_var * 0.5;
    let gamma1_2 = -spec.nu2 * spec.nu2 * res_var * 0.5;
    let gamma2_1 = spec.nu1.powf(4.0) * res_var * avg_k_cov * 0.25;
    let gamma2_2 = spec.nu2.powf(4.0) * res_var * avg_k_cov * 0.25;
    let gamma3_1 = spec.nu1.powf(6.0) * res_var * avg_k_cov * avg_k_cov * (1.0 / 12.0);
    let gamma3_2 = spec.nu2.powf(6.0) * res_var * avg_k_cov * avg_k_cov * (1.0 / 12.0);

    // Express VIX in percentage points (100 * sqrt(V)) matching evaluate_quintic_ou_spx_vix
    let scale_xi0 = spec.xi0.max_val(T::from_f64(1e-10)) * 10_000.0;
    let inv_sqrt_pi = 1.0 / std::f64::consts::PI.sqrt();
    let eps = T::from_f64(1e-3);

    let mut fut_p0 = T::zero();
    let mut fut_p1_1 = T::zero();
    let mut fut_p1_2 = T::zero();
    let mut fut_p2_1 = T::zero();
    let mut fut_p2_2 = T::zero();
    let mut fut_p3_1 = T::zero();
    let mut fut_p3_2 = T::zero();

    let mut call_p0 = T::zero();
    let mut call_p1_1 = T::zero();
    let mut call_p1_2 = T::zero();
    let mut call_p2_1 = T::zero();
    let mut call_p2_2 = T::zero();
    let mut call_p3_1 = T::zero();
    let mut call_p3_2 = T::zero();

    for m in 0..9 {
        let z = T::from_f64(GH9_X[m] * std::f64::consts::SQRT_2);
        let w = GH9_W[m] * inv_sqrt_pi;
        let comp1 = scale_xi0 * lam * (mu1 + sigma1 * z).exp();
        let comp2 = scale_xi0 * (T::one() - lam) * (mu2 + sigma2 * z).exp();
        let vp = (comp1 + comp2).max_val(T::from_f64(1e-10));
        let vix_z = vp.sqrt();
        let dphi_fut = T::from_f64(0.5) / vix_z;

        let psi1_f = dphi_fut * comp1;
        let psi2_f = dphi_fut * comp2;
        fut_p0 += vix_z * w;
        fut_p1_1 += psi1_f * w;
        fut_p1_2 += psi2_f * w;
        fut_p2_1 += (z * psi1_f / sigma1) * w;
        fut_p2_2 += (z * psi2_f / sigma2) * w;
        fut_p3_1 += ((z * z - T::one()) * psi1_f / (sigma1 * sigma1)) * w;
        fut_p3_2 += ((z * z - T::one()) * psi2_f / (sigma2 * sigma2)) * w;

        let diff = vix_z - spec.strike;
        let rad = (diff * diff + eps * eps).sqrt();
        let call_payoff = (diff + rad - eps) * 0.5;
        let dcall_dvix = (T::one() + diff / rad) * 0.5;
        let dphi_call = dcall_dvix * dphi_fut;

        let psi1_c = dphi_call * comp1;
        let psi2_c = dphi_call * comp2;
        call_p0 += call_payoff * w;
        call_p1_1 += psi1_c * w;
        call_p1_2 += psi2_c * w;
        call_p2_1 += (z * psi1_c / sigma1) * w;
        call_p2_2 += (z * psi2_c / sigma2) * w;
        call_p3_1 += ((z * z - T::one()) * psi1_c / (sigma1 * sigma1)) * w;
        call_p3_2 += ((z * z - T::one()) * psi2_c / (sigma2 * sigma2)) * w;
    }

    let vix_fut = (fut_p0
        + gamma1_1 * fut_p1_1
        + gamma1_2 * fut_p1_2
        + gamma2_1 * fut_p2_1
        + gamma2_2 * fut_p2_2
        + gamma3_1 * fut_p3_1
        + gamma3_2 * fut_p3_2)
        .max_val(T::from_f64(1e-4));

    let raw_call = call_p0
        + gamma1_1 * call_p1_1
        + gamma1_2 * call_p1_2
        + gamma2_1 * call_p2_1
        + gamma2_2 * call_p2_2
        + gamma3_1 * call_p3_1
        + gamma3_2 * call_p3_2;
    let intrinsic = (vix_fut - spec.strike).max_val(T::zero());
    let vix_call_undisc = raw_call.max_val(intrinsic);
    let vix_put_undisc = (vix_call_undisc - (vix_fut - spec.strike)).max_val(T::zero());
    let vix_call_price = spec.discount_factor * vix_call_undisc;
    let vix_put_price = spec.discount_factor * vix_put_undisc;

    // Invert undiscounted VIX call via Black-76 + ift_lift_1d
    let f_val = vix_fut.val();
    let k_val = spec.strike.val();
    let t_val = t_exp.val().max(1e-8);
    let eff_sig = (lam * sigma1 + (T::one() - lam) * sigma2) * 0.5 / t_exp.sqrt();
    let target_val = vix_call_undisc.val();
    let vix_iv = if target_val <= (f_val - k_val).max(0.0) + 1e-5 {
        eff_sig
    } else {
        let bs_undisc_f64 = |vol: f64| -> f64 {
            ql_vol::black_scholes_price(f_val, k_val, vol, t_val, 0.0, 0.0, true)
        };
        let mut lo = 1e-4_f64;
        let mut hi = 5.0_f64;
        for _ in 0..36 {
            let mid = 0.5 * (lo + hi);
            if bs_undisc_f64(mid) < target_val {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let mut iv_root = 0.5 * (lo + hi);
        let h_iv = 1e-5;
        for _ in 0..4 {
            let r = bs_undisc_f64(iv_root) - target_val;
            let dr = (bs_undisc_f64(iv_root + h_iv) - bs_undisc_f64(iv_root - h_iv)) / (2.0 * h_iv);
            if dr.abs() > 1e-12 {
                iv_root = (iv_root - r / dr).clamp(1e-4, 5.0);
            }
        }
        let f_xx = (bs_undisc_f64(iv_root + h_iv)
            - 2.0 * bs_undisc_f64(iv_root)
            + bs_undisc_f64(iv_root - h_iv))
            / (h_iv * h_iv);
        let iv_t = T::from_f64(iv_root);
        let h_t = T::from_f64(h_iv);
        let c_up = ql_vol::black_scholes_price(
            vix_fut,
            spec.strike,
            iv_t + h_t,
            t_exp,
            T::zero(),
            T::zero(),
            true,
        );
        let c_dn = ql_vol::black_scholes_price(
            vix_fut,
            spec.strike,
            iv_t - h_t,
            t_exp,
            T::zero(),
            T::zero(),
            true,
        );
        let f_x = (c_up - c_dn) / (2.0 * h_iv);
        ql_ad::ift_lift_1d(
            iv_root,
            |s| {
                ql_vol::black_scholes_price(
                    vix_fut,
                    spec.strike,
                    s,
                    t_exp,
                    T::zero(),
                    T::zero(),
                    true,
                ) - vix_call_undisc
            },
            f_x,
            f_xx,
        )
    };

    // Bourgey-De Marco-Gobet Theorem 3.1 Eq. 3.1: positive VIX ATM skew generated by nu1 != nu2
    let vix_atm_skew_proxy =
        lam * (T::one() - lam) * (sigma2 - sigma1) * (sigma2 - sigma1) / (eff_sig * t_exp);

    MixedRoughBergomiVixReport {
        vix_future: vix_fut,
        vix_call_price,
        vix_put_price,
        vix_implied_vol: vix_iv,
        proxy_sigma1: sigma1,
        proxy_sigma2: sigma2,
        vix_atm_skew_proxy,
    }
}



#[cfg(test)]
mod phase20_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn mixed_rough_bergomi_vix_positive_skew_put_call_parity_and_ad() {
        let spec = MixedRoughBergomiVixSpec {
            xi0: 0.04,
            hurst: 0.12,
            mixing_weight_lambda: 0.55,
            nu1: 0.9,
            nu2: 2.2,
            expiry: 0.25,
            vix_window: 1.0 / 12.0,
            strike: 20.0,
            discount_factor: 0.992,
        };
        let rep = evaluate_mixed_rough_bergomi_vix(&spec);
        assert!(rep.vix_future > 15.0 && rep.vix_future < 25.0);
        assert!(rep.vix_atm_skew_proxy > 0.0);
        assert!(
            (rep.vix_call_price - rep.vix_put_price - 0.992 * (rep.vix_future - 20.0)).abs()
                < 1e-12
        );
        let bs_reprice = ql_vol::black_scholes_price(
            rep.vix_future,
            20.0,
            rep.vix_implied_vol,
            0.25,
            0.0,
            0.0,
            true,
        ) * 0.992;
        assert!((bs_reprice - rep.vix_call_price).abs() < 1e-6);

        let h = 1e-5;
        let eval_xi = |xi: f64| {
            let mut s = spec;
            s.xi0 = xi;
            evaluate_mixed_rough_bergomi_vix(&s).vix_call_price
        };
        let fd = (eval_xi(0.04 + h) - eval_xi(0.04 - h)) / (2.0 * h);
        type D = Dual<1>;
        let spec_d = MixedRoughBergomiVixSpec {
            xi0: D::variable(0.04, 0),
            hurst: D::from_value(0.12),
            mixing_weight_lambda: D::from_value(0.55),
            nu1: D::from_value(0.9),
            nu2: D::from_value(2.2),
            expiry: D::from_value(0.25),
            vix_window: 1.0 / 12.0,
            strike: D::from_value(20.0),
            discount_factor: D::from_value(0.992),
        };
        let ad = evaluate_mixed_rough_bergomi_vix(&spec_d).vix_call_price;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-4);
    }
}
