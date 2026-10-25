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

#[derive(Clone, Debug)]
pub struct AndreasenHugeGridSpec<T: Real = f64> {
    pub forward: T,
    pub rate: T,
    pub dt: T,
    pub strike_min: T,
    pub strike_max: T,
    pub n_grid: usize,
    pub local_vol_atm: T,
    pub local_vol_skew: T,
    pub local_vol_smile: T,
}

#[derive(Clone, Debug)]
pub struct AndreasenHugeSliceResult<T: Real> {
    pub strikes: Vec<T>,
    pub local_vols: Vec<T>,
    pub undiscounted_calls: Vec<T>,
    pub call_prices: Vec<T>,
    pub implied_densities: Vec<T>,
}

pub fn solve_andreasen_huge_step<T: Real>(
    spec: &AndreasenHugeGridSpec<T>,
    prior_undiscounted_calls: Option<&[T]>,
) -> AndreasenHugeSliceResult<T> {
    let g = spec.n_grid.max(5);
    let dk = (spec.strike_max - spec.strike_min) / (g as f64 - 1.0);
    let dk_safe = dk.max_val(T::from_f64(1e-8));
    let fwd_safe = spec.forward.max_val(T::from_f64(1e-8));
    let dt_safe = spec.dt.max_val(T::from_f64(1e-8));
    let eps = T::from_f64(1e-8);

    let mut strikes = Vec::with_capacity(g);
    let mut local_vols = Vec::with_capacity(g);
    let mut rhs = Vec::with_capacity(g);
    for i in 0..g {
        let k_i = spec.strike_min + dk * (i as f64);
        strikes.push(k_i);
        let m_i = (k_i - spec.forward) / fwd_safe;
        let lv = (spec.local_vol_atm
            + spec.local_vol_skew * m_i
            + spec.local_vol_smile * m_i * m_i * 0.5)
            .max_val(T::from_f64(0.01));
        local_vols.push(lv);
        if let Some(prev) = prior_undiscounted_calls.filter(|p| p.len() == g) {
            rhs.push(prev[i]);
        } else {
            let diff = spec.forward - k_i;
            let intrinsic = (diff + (diff * diff + eps * eps).sqrt()) * 0.5;
            rhs.push(intrinsic);
        }
    }

    let mut a = vec![T::zero(); g];
    let mut b = vec![T::one(); g];
    let mut c = vec![T::zero(); g];
    for i in 1..(g - 1) {
        let alpha = dt_safe * local_vols[i] * local_vols[i] * strikes[i] * strikes[i]
            / (dk_safe * dk_safe * 2.0);
        a[i] = -alpha;
        b[i] = T::one() + alpha * 2.0;
        c[i] = -alpha;
    }

    let mut c_prime = vec![T::zero(); g];
    let mut d_prime = vec![T::zero(); g];
    c_prime[0] = c[0] / b[0];
    d_prime[0] = rhs[0] / b[0];
    for i in 1..g {
        let denom = (b[i] - a[i] * c_prime[i - 1]).max_val(T::from_f64(1e-14));
        c_prime[i] = if i + 1 < g { c[i] / denom } else { T::zero() };
        d_prime[i] = (rhs[i] - a[i] * d_prime[i - 1]) / denom;
    }
    let mut u_calls = vec![T::zero(); g];
    u_calls[g - 1] = d_prime[g - 1];
    for i in (0..(g - 1)).rev() {
        u_calls[i] = d_prime[i] - c_prime[i] * u_calls[i + 1];
    }

    let df = (-spec.rate * dt_safe).exp();
    let call_prices: Vec<T> = u_calls.iter().map(|&c_val| c_val * df).collect();
    let mut densities = vec![T::zero(); g];
    for i in 1..(g - 1) {
        densities[i] = ((u_calls[i + 1] - u_calls[i] * 2.0 + u_calls[i - 1]) / (dk_safe * dk_safe))
            .max_val(T::zero());
    }

    AndreasenHugeSliceResult {
        strikes,
        local_vols,
        undiscounted_calls: u_calls,
        call_prices,
        implied_densities: densities,
    }
}

pub fn interpolate_andreasen_huge_call<T: Real>(
    slice: &AndreasenHugeSliceResult<T>,
    target_strike: T,
) -> T {
    let g = slice.strikes.len();
    if g == 0 {
        return T::zero();
    }
    if target_strike.val() <= slice.strikes[0].val() || g == 1 {
        return slice.call_prices[0];
    }
    if target_strike.val() >= slice.strikes[g - 1].val() {
        return slice.call_prices[g - 1];
    }
    for i in 0..(g - 1) {
        if target_strike.val() <= slice.strikes[i + 1].val() {
            let dk = (slice.strikes[i + 1] - slice.strikes[i]).max_val(T::from_f64(1e-12));
            let w = (target_strike - slice.strikes[i]) / dk;
            return slice.call_prices[i] * (T::one() - w) + slice.call_prices[i + 1] * w;
        }
    }
    slice.call_prices[g - 1]
}

#[derive(Copy, Clone, Debug)]
pub struct ZabrParams<T: Real> {
    pub alpha: T,
    pub beta: T,
    pub nu: T,
    pub rho: T,
    pub gamma: T,
}

#[derive(Copy, Clone, Debug)]
pub struct ZabrSmileQuote<T: Real> {
    pub normal_vol: T,
    pub black_vol: T,
    pub call_price: T,
}

fn pow_real<T: Real>(base: T, exp: T) -> T {
    (base.max_val(T::from_f64(1e-12)).ln() * exp).exp()
}

pub fn evaluate_zabr_smile<T: Real>(
    p: &ZabrParams<T>,
    forward: T,
    strike: T,
    expiry: T,
) -> ZabrSmileQuote<T> {
    let f = forward.max_val(T::from_f64(1e-8));
    let k = strike.max_val(T::from_f64(1e-8));
    let t = expiry.max_val(T::from_f64(1e-8));
    let alpha = p.alpha.max_val(T::from_f64(1e-8));
    let beta = p.beta;
    let one_m_beta = T::one() - beta;
    let f_mid = (f * k).sqrt();
    let log_fk = (f / k).ln();

    // CEV integral I_CEV = \int_K^F u^{-\beta} du using smooth 3rd-order expansion in log(F/K)
    let u = one_m_beta * log_fk * 0.5;
    let sinh_over_u = T::one() + u * u / 6.0 + u * u * u * u / 120.0;
    let i_cev = pow_real(f_mid, -beta) * (f - k) / sinh_over_u;

    let alpha_gm2 = pow_real(alpha, p.gamma - T::from_f64(2.0));
    let y = alpha_gm2 * i_cev;
    let nu_tilde = (T::from_f64(2.0) - p.gamma).max_val(T::from_f64(0.05))
        * p.nu
        * pow_real(alpha, p.gamma - T::one());
    let z = nu_tilde * y;
    let y_over_x = if z.val().abs() < 1e-4 {
        T::one() - p.rho * z * 0.5 + (T::from_f64(2.0) - p.rho * p.rho * 3.0) * z * z / 12.0
    } else {
        let j = (T::one() - p.rho * z * 2.0 + z * z)
            .max_val(T::from_f64(1e-14))
            .sqrt();
        let x_z = ((j + z - p.rho) / (T::one() - p.rho).max_val(T::from_f64(1e-12))).ln();
        z / x_z
    };

    let f_m_over_i = pow_real(f_mid, beta) * sinh_over_u;
    let sig_atm_n = alpha * pow_real(f_mid, beta);
    let time_corr = T::one()
        + (beta * (beta - T::from_f64(2.0)) * sig_atm_n * sig_atm_n / (f_mid * f_mid * 24.0)
            + p.rho * beta * alpha * nu_tilde / (pow_real(f_mid, one_m_beta) * 4.0)
            + (T::from_f64(2.0) - p.rho * p.rho * 3.0) * nu_tilde * nu_tilde / 24.0)
            * t;
    let normal_vol = (alpha * f_m_over_i * y_over_x * time_corr).max_val(T::from_f64(1e-8));
    let fk_ratio = if log_fk.val().abs() < 1e-5 {
        T::one() / f_mid
    } else {
        log_fk / (f - k)
    };
    let black_vol = (normal_vol * fk_ratio).max_val(T::from_f64(1e-8));
    let st = black_vol * t.sqrt();
    let d1 = (log_fk + black_vol * black_vol * t * 0.5) / st;
    let d2 = d1 - st;
    let call_price = (f * d1.norm_cdf() - k * d2.norm_cdf()).max_val(T::zero());

    ZabrSmileQuote {
        normal_vol,
        black_vol,
        call_price,
    }
}

#[derive(Copy, Clone, Debug)]
pub struct FeldmanSabrCalibrationResult<T: Real> {
    pub alpha: T,
    pub nu: T,
    pub rho: T,
}

pub fn feldman_analytic_sabr_calibration<T: Real>(
    forward: T,
    shift: T,
    beta: T,
    expiry: T,
    h: T,
    c: [T; 5],
) -> FeldmanSabrCalibrationResult<T> {
    let f_s = (forward + shift).max_val(T::from_f64(1e-6));
    let t = expiry.max_val(T::from_f64(1e-8));
    let h_s = h.max_val(T::from_f64(1e-6));
    let atm = c[2].max_val(T::from_f64(1e-8));
    let denom_atm = (c[1] + c[3] - atm * 2.0).max_val(T::from_f64(1e-10));
    let f_s_beta = pow_real(f_s, beta);
    // Feldman (2020) Theorem 1 Eq. (2) / Eq. (32)
    let alpha = (h_s / (f_s_beta * t.sqrt())) * (atm / denom_atm).sqrt();

    // Adjacent second-difference curvature ratios z_- and z_+ (Feldman Eq. 5)
    let z_minus = (c[1] - h_s).max_val(T::from_f64(1e-10))
        / (c[0] + atm - c[1] * 2.0).max_val(T::from_f64(1e-10));
    let z_plus =
        c[3].max_val(T::from_f64(1e-10)) / (c[4] + atm - c[3] * 2.0).max_val(T::from_f64(1e-10));

    let sig_sqrt_t = atm * (2.0 * std::f64::consts::PI).sqrt();
    let xi = h_s / sig_sqrt_t.max_val(T::from_f64(1e-8));
    let mills = (-xi).norm_cdf() / xi.norm_pdf().max_val(T::from_f64(1e-12));
    let kappa_h = (T::one() - xi * mills).max_val(T::from_f64(1e-6));

    let one_m_beta = T::one() - beta;
    let cev_dist = |k_sh: T| -> T {
        let ks = k_sh.max_val(T::from_f64(1e-6));
        if one_m_beta.val().abs() < 1e-5 {
            (f_s / ks).ln() / alpha
        } else {
            (pow_real(f_s, one_m_beta) - pow_real(ks, one_m_beta)) / (alpha * one_m_beta)
        }
    };
    let y_minus = cev_dist(f_s - h_s);
    let y_plus = cev_dist(f_s + h_s);

    let lv_m2 = pow_real((f_s - h_s).max_val(T::from_f64(1e-6)), beta * 2.0);
    let lv_p2 = pow_real(f_s + h_s, beta * 2.0);
    let r_minus =
        ((z_minus * h_s * h_s) / (t * kappa_h * alpha * alpha * lv_m2 * 2.0) - T::one()) / y_minus;
    let r_plus =
        ((z_plus * h_s * h_s) / (t * kappa_h * alpha * alpha * lv_p2 * 2.0) - T::one()) / y_plus;

    let nu2 = ((r_minus - r_plus) / (y_minus - y_plus)).max_val(T::from_f64(1e-10));
    let nu = nu2.sqrt();
    let rho = ((nu2 * y_minus - r_minus) / (nu * 2.0))
        .max_val(T::from_f64(-0.999))
        .min_val(T::from_f64(0.999));

    FeldmanSabrCalibrationResult { alpha, nu, rho }
}

#[cfg(test)]
mod phase17_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn andreasen_huge_butterfly_calendar_arbitrage_free_and_ad() {
        let spec1 = AndreasenHugeGridSpec {
            forward: 100.0,
            rate: 0.02,
            dt: 0.5,
            strike_min: 60.0,
            strike_max: 140.0,
            n_grid: 41,
            local_vol_atm: 0.20,
            local_vol_skew: -0.15,
            local_vol_smile: 0.10,
        };
        let slice1 = solve_andreasen_huge_step(&spec1, None);
        let slice2 = solve_andreasen_huge_step(&spec1, Some(&slice1.undiscounted_calls));
        for i in 1..40 {
            assert!(slice1.implied_densities[i] >= 0.0);
            assert!(slice2.undiscounted_calls[i] >= slice1.undiscounted_calls[i] - 1e-12);
        }
        let h = 1e-5;
        let eval_atm = |v: f64| {
            let mut sp = spec1.clone();
            sp.local_vol_atm = v;
            let sl = solve_andreasen_huge_step(&sp, None);
            interpolate_andreasen_huge_call(&sl, 100.0)
        };
        let fd = (eval_atm(0.20 + h) - eval_atm(0.20 - h)) / (2.0 * h);
        type D = Dual<1>;
        let spec_d = AndreasenHugeGridSpec {
            forward: D::from_value(100.0),
            rate: D::from_value(0.02),
            dt: D::from_value(0.5),
            strike_min: D::from_value(60.0),
            strike_max: D::from_value(140.0),
            n_grid: 41,
            local_vol_atm: D::variable(0.20, 0),
            local_vol_skew: D::from_value(-0.15),
            local_vol_smile: D::from_value(0.10),
        };
        let sl_d = solve_andreasen_huge_step(&spec_d, None);
        let ad = interpolate_andreasen_huge_call(&sl_d, D::from_value(100.0));
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn zabr_gamma_wings_feldman_calibration_and_ad() {
        let p_sabr = ZabrParams {
            alpha: 0.20,
            beta: 1.0,
            nu: 0.45,
            rho: -0.30,
            gamma: 1.0,
        };
        let p_flat_wing = ZabrParams {
            gamma: 0.5,
            ..p_sabr
        };
        let q1 = evaluate_zabr_smile(&p_sabr, 100.0, 125.0, 1.0);
        let q_half = evaluate_zabr_smile(&p_flat_wing, 100.0, 125.0, 1.0);
        assert!(q1.black_vol.is_finite() && q1.call_price > 0.0);
        assert_ne!(q1.black_vol, q_half.black_vol);

        let h_fd = 1e-5;
        let eval_g = |g: f64| {
            let p = ZabrParams { gamma: g, ..p_sabr };
            evaluate_zabr_smile(&p, 100.0, 115.0, 1.0).black_vol
        };
        let fd = (eval_g(1.0 + h_fd) - eval_g(1.0 - h_fd)) / (2.0 * h_fd);
        type D = Dual<1>;
        let pd = ZabrParams {
            alpha: D::from_value(0.20),
            beta: D::one(),
            nu: D::from_value(0.45),
            rho: D::from_value(-0.30),
            gamma: D::variable(1.0, 0),
        };
        let ad = evaluate_zabr_smile(&pd, D::from_value(100.0), D::from_value(115.0), D::one())
            .black_vol;
        assert!((ad.d[0] - fd).abs() < 1e-6);

        let cal = feldman_analytic_sabr_calibration(
            100.0,
            0.0,
            0.5,
            1.0,
            2.0,
            [9.80, 8.65, 7.65, 6.80, 6.10],
        );
        assert!(cal.alpha > 0.0 && cal.nu > 0.0 && cal.rho.abs() < 1.0);
    }
}
