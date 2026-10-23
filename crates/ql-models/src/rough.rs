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
    if n <= 1 {
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
fn phi1(x: f64) -> f64 {
    if x.abs() < 1e-8 {
        1.0 + 0.5 * x + x * x / 6.0 + x * x * x / 24.0
    } else {
        (x.exp() - 1.0) / x
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
    if p.n_factors <= 1 {
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
    let n_steps = steps.max(1);
    let dt = expiry / T::from_f64(n_steps as f64);
    let dt_val = dt.val();
    let half_exp: Vec<f64> = gamma.iter().map(|&g| (-0.5 * g * dt_val).exp()).collect();
    let half_phi: Vec<f64> = gamma
        .iter()
        .map(|&g| 0.5 * phi1(-0.5 * g * dt_val))
        .collect();
    let full_exp: Vec<f64> = gamma.iter().map(|&g| (-g * dt_val).exp()).collect();
    let full_phi: Vec<f64> = gamma.iter().map(|&g| phi1(-g * dt_val)).collect();

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
            let mid = psi[i] * T::from_f64(half_exp[i]) + f_k * (dt * T::from_f64(half_phi[i]));
            mid_aggregate = mid_aggregate + mid * T::from_f64(c[i]);
        }
        let f_mid = riccati(mid_aggregate);
        for i in 0..n {
            psi[i] = psi[i] * T::from_f64(full_exp[i]) + f_mid * (dt * T::from_f64(full_phi[i]));
        }
        let tau_mid = expiry - dt * T::from_f64(step as f64 + 0.5);
        let kernel_mid: f64 = c
            .iter()
            .zip(gamma.iter())
            .map(|(&ci, &gi)| ci * phi1(-gi * tau_mid.val()))
            .sum();
        let g_mid = p.v0 + p.lambda * p.theta * tau_mid * T::from_f64(kernel_mid);
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
        spot,
        rate,
        div,
        expiry,
        std::slice::from_ref(&strike),
        is_call,
    )[0]
}
