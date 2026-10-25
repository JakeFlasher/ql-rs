use crate::{Complex, CosHestonEngine, HestonParams};
use ql_ad::Real;
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LiftedHestonParams<T: Real> {
    pub v0: T,
    pub lambda: T,
    pub theta: T,
    pub nu: T,
    pub rho: T,
    pub hurst: f64,
    pub n_factors: usize,
    pub r_n: f64,
}
pub fn lifted_heston_kernel(h: f64, n: usize, r_n: f64) -> (Vec<f64>, Vec<f64>) {
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    if (h - 0.5).abs() < 1e-8 {
        return (vec![1.0 / n as f64; n], vec![0.0; n]);
    }
    if n == 1 {
        return (vec![1.], vec![0.]);
    }
    let a = h + 0.5;
    let den = libm::tgamma(a) * libm::tgamma(2. - a);
    let mut c = Vec::with_capacity(n);
    let mut g = Vec::with_capacity(n);
    for i in 1..=n {
        c.push(
            (r_n.powf(1. - a) - 1.) * r_n.powf((a - 1.) * (1. + n as f64 / 2.)) / den
                * r_n.powf((1. - a) * i as f64),
        );
        g.push(
            (1. - a) / (2. - a) * (r_n.powf(2. - a) - 1.) / (r_n.powf(1. - a) - 1.)
                * r_n.powf(i as f64 - 1. - n as f64 / 2.),
        );
    }
    (c, g)
}

#[inline]
fn phi1<T: Real>(x: T) -> T {
    if x.val().abs() < 1e-6 {
        T::one() - x * 0.5 + x * x / 6.0 - x * x * x / 24.0
    } else {
        (T::one() - (-x).exp()) / x
    }
}

/// Characteristic function of the Markovian lifted rough-Heston model with
/// unit spot.  The factors are integrated with the exponential-midpoint
/// scheme, while the common Riccati non-linearity is evaluated at the
/// weighted aggregate of the factor states.
#[allow(clippy::too_many_arguments)]
pub fn lifted_heston_char_func_unit_spot<T: Real>(
    u: T,
    expiry: T,
    rate: T,
    div: T,
    p: &LiftedHestonParams<T>,
    steps: usize,
) -> Complex<T> {
    let z = Complex::new(T::zero(), u);

    // A one-factor lift is exactly the classical Heston model.  Keeping this
    // branch also avoids a loss of accuracy in the short integration for the
    // limiting kernel gamma=0.
    if (p.hurst - 0.5).abs() < 1e-8 {
        let h = HestonParams {
            v0: p.v0,
            kappa: p.lambda,
            theta: p.theta,
            sigma_v: p.nu,
            rho: p.rho,
        };
        return h.char_func_unit_spot(Complex::from_real(u), rate, div, expiry);
    }

    let (c, gamma) = lifted_heston_kernel(p.hurst, p.n_factors, p.r_n);
    let n = c.len();
    // A lifted Riccati system is considerably stiffer than the classical
    // one-factor limit.  Enforce a conservative minimum subdivision so small
    // caller step counts cannot drive the exponential-midpoint update into a
    // non-finite state (the requested count remains an upper-level accuracy
    // knob for counts >= 64).
    let n_steps = steps.max(64);
    let dt = expiry / T::from_f64(n_steps as f64);
    let half_exp: Vec<T> = gamma
        .iter()
        .map(|&g| (-T::from_f64(0.5 * g) * dt).exp())
        .collect();
    let half_phi: Vec<T> = gamma
        .iter()
        .map(|&g| phi1(T::from_f64(0.5 * g) * dt) * 0.5)
        .collect();
    let full_exp: Vec<T> = gamma
        .iter()
        .map(|&g| (-T::from_f64(g) * dt).exp())
        .collect();
    let full_phi: Vec<T> = gamma.iter().map(|&g| phi1(T::from_f64(g) * dt)).collect();

    let riccati = |psi: Complex<T>| -> Complex<T> {
        let half = T::from_f64(0.5);
        (z * z - z) * half
            + (z * p.rho * p.nu - p.lambda) * psi
            + (psi * psi) * (half * p.nu * p.nu)
    };
    let mut psi = vec![Complex::zero(); n];
    let mut phi_int = Complex::zero();
    for step in 0..n_steps {
        let mut aggregate = Complex::zero();
        for i in 0..n {
            aggregate = aggregate + psi[i] * T::from_f64(c[i]);
        }
        let f_k = riccati(aggregate);
        let mut mid_aggregate = Complex::zero();
        for i in 0..n {
            let mid = psi[i] * half_exp[i] + f_k * (dt * half_phi[i]);
            mid_aggregate = mid_aggregate + mid * T::from_f64(c[i]);
        }
        let f_mid = riccati(mid_aggregate);
        for i in 0..n {
            psi[i] = psi[i] * full_exp[i] + f_mid * (dt * full_phi[i]);
        }
        let tau_mid = expiry - dt * T::from_f64(step as f64 + 0.5);
        let kernel_mid = c
            .iter()
            .zip(gamma.iter())
            .map(|(&ci, &gi)| T::from_f64(ci) * phi1(T::from_f64(gi) * tau_mid))
            .fold(T::zero(), |acc, x| acc + x);
        let g_mid = p.v0 + p.lambda * p.theta * tau_mid * kernel_mid;
        phi_int = phi_int + f_mid * (dt * g_mid);
    }
    (z * (rate - div) * expiry + phi_int).exp()
}

pub fn lifted_heston_char_fn<T: Real>(
    u: T,
    expiry: T,
    spot: T,
    rate: T,
    div: T,
    p: &LiftedHestonParams<T>,
    steps: usize,
) -> Complex<T> {
    let z = Complex::new(T::zero(), u);
    lifted_heston_char_func_unit_spot(u, expiry, rate, div, p, steps) * (z * spot.ln()).exp()
}
#[allow(clippy::too_many_arguments)]
pub fn lifted_heston_cos_price<T: Real>(
    spot: T,
    strike: T,
    expiry: T,
    rate: T,
    div: T,
    p: &LiftedHestonParams<T>,
    n_terms: usize,
    l_trunc: f64,
    ode_steps: usize,
    is_call: bool,
) -> T {
    CosHestonEngine { n_terms, l_trunc }.price_impl(
        |u| lifted_heston_char_func_unit_spot(u.re, expiry, rate, div, p, ode_steps),
        p.v0,
        p.theta,
        p.nu,
        p.lambda,
        p.rho,
        spot,
        rate,
        div,
        expiry,
        std::slice::from_ref(&strike),
        is_call,
    )[0]
}

#[derive(Clone, Debug, PartialEq)]
pub struct RoughBergomiParams<T: Real> {
    pub xi0: T,
    pub eta: T,
    pub rho: T,
    pub hurst: f64,
}

#[allow(clippy::too_many_arguments)]
pub fn rbergomi_hybrid_mc_price<T: Real>(
    spot: T,
    strike: T,
    expiry: T,
    rate: T,
    div: T,
    params: &RoughBergomiParams<T>,
    n_steps: usize,
    n_paths: usize,
    is_call: bool,
    seed: u64,
) -> T {
    assert!((0.01..0.499).contains(&params.hurst));
    let n = n_steps.max(2);
    let dt = expiry / T::from_f64(n as f64);
    let h = params.hurst;
    let ch = (2. * h).sqrt() / (h + 0.5);
    let sh = (1. - ch * ch).max(0.).sqrt();
    let mut paths: Vec<(T, T)> = Vec::with_capacity(n_paths.max(2));
    fn nr(s: &mut u64) -> f64 {
        *s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
        let u = ((*s >> 11) as f64 + 0.5) / 9007199254740992.;
        *s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
        let v = ((*s >> 11) as f64 + 0.5) / 9007199254740992.;
        (-2. * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()
    }
    let mut weights = vec![0.; n + 1];
    for (m, weight) in weights.iter_mut().enumerate().take(n + 1).skip(2) {
        *weight = ((m as f64).powf(h + 0.5) - ((m - 1) as f64).powf(h + 0.5)) / (h + 0.5);
    }
    for p in 0..n_paths.max(2) {
        let sign = if p % 2 == 0 { 1. } else { -1. };
        let mut rng = seed.wrapping_add((p / 2) as u64);
        let mut dw = vec![T::zero(); n];
        let mut wh = vec![T::zero(); n];
        for j in 0..n {
            let u1 = nr(&mut rng) * sign;
            let u2 = nr(&mut rng) * sign;
            dw[j] = dt.sqrt() * T::from_f64(u1);
            wh[j] = dt.powf(h) * T::from_f64(ch * u1 + sh * u2);
        }
        let mut s_par = spot * ((rate - div) * expiry).exp();
        let mut qv = T::zero();
        let mut v = params.xi0;
        for k in 0..n {
            s_par *= (params.rho * v.sqrt() * dw[k] - params.rho * params.rho * v * dt * 0.5).exp();
            qv += (T::one() - params.rho * params.rho).max_val(T::zero()) * v * dt;
            let mut conv = wh[k];
            for j in 0..k {
                conv += T::from_f64((2. * h).sqrt())
                    * (T::from_f64(weights[k + 1 - j]) * dt.powf(h - 0.5))
                    * dw[j];
            }
            let t = T::from_f64((k + 1) as f64) * dt;
            v = params.xi0
                * (params.eta * conv - params.eta * params.eta * t.powf(2. * h) * 0.5).exp();
        }
        let vol = (qv / expiry.max_val(T::from_f64(1e-12)))
            .max_val(T::from_f64(1e-16))
            .sqrt();
        paths.push((s_par, vol));
    }
    let target_forward = spot * ((rate - div) * expiry).exp();
    let mean_spot = paths
        .iter()
        .map(|(s_path, _)| *s_path)
        .fold(T::zero(), |a, z| a + z)
        / paths.len() as f64;
    let martingale_scale = target_forward / mean_spot.max_val(T::from_f64(1e-12));
    let sum = paths
        .into_iter()
        .map(|(s_path, vol)| {
            ql_vol::black_scholes_price(
                s_path * martingale_scale,
                strike,
                vol,
                expiry,
                T::zero(),
                T::zero(),
                is_call,
            ) * (-rate * expiry).exp()
        })
        .fold(T::zero(), |a, z| a + z);
    sum / n_paths.max(2) as f64
}

#[derive(Copy, Clone, Debug)]
pub struct QuadraticRoughHestonParams<T: Real = f64> {
    pub hurst: f64,
    pub nu: T,
    pub lambda: T,
    pub a: T,
    pub b: T,
    pub c: T,
    pub z0: T,
}

pub fn qrh_implied_vol_expansion<T: Real>(
    p: &QuadraticRoughHestonParams<T>,
    forward: T,
    strike: T,
    expiry: T,
) -> T {
    let h = p.hurst.clamp(0.02, 0.49);
    let t = expiry.max_val(T::from_f64(1e-8));
    let zb = p.z0 - p.b;
    let v0 = (p.a * zb * zb + p.c).max_val(T::from_f64(1e-10));
    let gh = libm::tgamma(h + 1.5);
    let th = t.powf(h);
    let eta_h = p.nu * th / (T::from_f64(gh) + p.lambda * t.powf(h + 0.5));
    let atm_var = v0 + p.a * eta_h * eta_h * v0 * (2.0 * h / (2.0 * h + 1.0)) * 0.5;
    let atm = atm_var.max_val(T::from_f64(1e-10)).sqrt();
    let k = (strike / forward).ln();
    let skew = p.a * zb * eta_h / (atm * atm * t.sqrt());
    let curv = eta_h * eta_h / (atm * atm * atm * t)
        * (p.a + p.a * p.a * zb * zb * 2.0 / (atm * atm) * ((1.0 - 2.0 * h) / (1.0 + 2.0 * h)));
    let total =
        (T::one() + skew * k * 2.0 + (curv + skew * skew) * k * k).max_val(T::from_f64(0.04));
    atm * total.sqrt()
}

#[allow(clippy::too_many_arguments)]
pub fn qrh_mc_price<T: Real>(
    p: &QuadraticRoughHestonParams<T>,
    spot: T,
    strike: T,
    rate: T,
    div: T,
    expiry: T,
    n_steps: usize,
    n_paths: usize,
    is_call: bool,
    seed: u64,
) -> T {
    let h = p.hurst.clamp(0.02, 0.49);
    let steps = n_steps.max(4);
    let paths = n_paths.max(2);
    let dt = expiry / steps as f64;
    let sqrt_dt = dt.sqrt();
    let gh = libm::tgamma(h + 1.5);
    let dt_h = dt.powf(h + 0.5);
    let mut w_k = Vec::with_capacity(steps + 1);
    w_k.push(T::zero());
    for k in 1..=steps {
        let diff = (k as f64).powf(h + 0.5) - ((k - 1) as f64).powf(h + 0.5);
        w_k.push(dt_h * (diff / gh));
    }
    let rng = ql_math::PhiloxRng::new(seed);
    let mut term_spots = Vec::with_capacity(paths);
    let mut mean_spot = T::zero();
    for path in 0..paths {
        let anti = path % 2 == 1;
        let id = (path / 2) as u64;
        let mut log_s = spot.ln();
        let mut z = p.z0;
        let mut dz_hist = vec![T::zero(); steps];
        for i in 0..steps {
            let zb = z - p.b;
            let v = (p.a * zb * zb + p.c).max_val(T::from_f64(1e-10));
            let eps = rng.normal4(id, i as u32)[0] * if anti { -1.0 } else { 1.0 };
            let dw = sqrt_dt * eps;
            log_s += (rate - div - v * 0.5) * dt + v.sqrt() * dw;
            dz_hist[i] = -p.lambda * z * dt - p.nu * v.sqrt() * dw;
            let mut conv = T::zero();
            for (j, &dz_j) in dz_hist.iter().enumerate().take(i + 1) {
                conv += w_k[i + 1 - j] * dz_j / dt;
            }
            z = p.z0 + conv;
        }
        let s_t = log_s.exp();
        mean_spot += s_t;
        term_spots.push(s_t);
    }
    mean_spot = mean_spot / paths as f64;
    let target_fwd = spot * ((rate - div) * expiry).exp();
    let scale = target_fwd / mean_spot.max_val(T::from_f64(1e-12));
    let df = (-rate * expiry).exp();
    let mut sum = T::zero();
    for s_raw in term_spots {
        let s_t = s_raw * scale;
        let diff = if is_call { s_t - strike } else { strike - s_t };
        let payoff = (diff + (diff * diff + 1e-8).sqrt()) * 0.5;
        sum += payoff * df;
    }
    sum / paths as f64
}

#[cfg(test)]
mod phase16_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn qrh_flat_limit_and_ad() {
        let p_flat = QuadraticRoughHestonParams {
            hurst: 0.1,
            nu: 0.3,
            lambda: 1.0,
            a: 0.0,
            b: 0.2,
            c: 0.04,
            z0: 0.0,
        };
        let iv = qrh_implied_vol_expansion(&p_flat, 100.0, 95.0, 1.0);
        assert!((iv - 0.2).abs() < 1e-12);
        let p = QuadraticRoughHestonParams {
            hurst: 0.1,
            nu: 0.4,
            lambda: 1.2,
            a: 0.5,
            b: 0.25,
            c: 0.02,
            z0: 0.0,
        };
        let iv_put = qrh_implied_vol_expansion(&p, 100.0, 90.0, 1.0);
        let iv_call = qrh_implied_vol_expansion(&p, 100.0, 110.0, 1.0);
        assert!(iv_put > iv_call);
        let mc_c = qrh_mc_price(&p, 100.0, 100.0, 0.02, 0.01, 1.0, 16, 128, true, 7);
        let mc_p = qrh_mc_price(&p, 100.0, 100.0, 0.02, 0.01, 1.0, 16, 128, false, 7);
        let parity = 100.0 * (-0.01f64).exp() - 100.0 * (-0.02f64).exp();
        assert!((mc_c - mc_p - parity).abs() < 0.35);
        let h = 1e-5;
        let eval_b = |b_val: f64| {
            let mut p2 = p;
            p2.b = b_val;
            qrh_implied_vol_expansion(&p2, 100.0, 92.0, 1.0)
        };
        let fd = (eval_b(0.25 + h) - eval_b(0.25 - h)) / (2.0 * h);
        type D = Dual<1>;
        let pd = QuadraticRoughHestonParams {
            hurst: 0.1,
            nu: D::from_value(0.4),
            lambda: D::from_value(1.2),
            a: D::from_value(0.5),
            b: D::variable(0.25, 0),
            c: D::from_value(0.02),
            z0: D::zero(),
        };
        let ad = qrh_implied_vol_expansion(
            &pd,
            D::from_value(100.0),
            D::from_value(92.0),
            D::from_value(1.0),
        );
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}
