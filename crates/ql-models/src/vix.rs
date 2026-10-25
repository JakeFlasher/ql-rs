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
    let moneyness = (vix_strike / vix_fut.max_val(T::from_f64(1e-6))).ln();
    let vix_iv = eff_vvol * (T::one() + p.alpha3 * 0.15 * moneyness + p.alpha5 * 0.08 * moneyness * moneyness);
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
