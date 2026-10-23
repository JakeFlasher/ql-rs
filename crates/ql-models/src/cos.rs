use crate::{BatesParams, Complex, HestonParams};
use ql_ad::Real;
#[derive(Copy, Clone, Debug)]
pub struct CosHestonEngine {
    pub n_terms: usize,
    pub l_trunc: f64,
}
impl Default for CosHestonEngine {
    fn default() -> Self {
        Self {
            n_terms: 160,
            l_trunc: 12.,
        }
    }
}

/// Convenience entry point for callers that want a configured COS Heston
/// price without constructing an engine object.
#[allow(clippy::too_many_arguments)]
pub fn heston_cos_price<T: Real>(
    spot: T,
    strike: T,
    expiry: T,
    rate: T,
    div: T,
    params: &HestonParams<T>,
    n_terms: usize,
    l_trunc: f64,
    is_call: bool,
) -> T {
    CosHestonEngine { n_terms, l_trunc }.price(params, spot, strike, rate, div, expiry, is_call)
}

impl CosHestonEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        &self,
        p: &HestonParams<T>,
        spot: T,
        strike: T,
        rate: T,
        div: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        self.price_strikes(
            p,
            spot,
            rate,
            div,
            expiry,
            std::slice::from_ref(&strike),
            is_call,
        )[0]
    }
    #[allow(clippy::too_many_arguments)]
    pub fn price_bates<T: Real>(
        &self,
        p: &BatesParams<T>,
        spot: T,
        strike: T,
        rate: T,
        div: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        self.price_bates_strikes(
            p,
            spot,
            rate,
            div,
            expiry,
            std::slice::from_ref(&strike),
            is_call,
        )[0]
    }
    #[allow(clippy::too_many_arguments)]
    pub fn price_strikes<T: Real>(
        &self,
        p: &HestonParams<T>,
        spot: T,
        rate: T,
        div: T,
        expiry: T,
        strikes: &[T],
        is_call: bool,
    ) -> Vec<T> {
        self.price_impl(
            |u| p.char_func_unit_spot(u, rate, div, expiry),
            p.v0,
            p.theta,
            p.sigma_v,
            p.kappa,
            p.rho,
            spot,
            rate,
            div,
            expiry,
            strikes,
            is_call,
        )
    }
    #[allow(clippy::too_many_arguments)]
    pub fn price_bates_strikes<T: Real>(
        &self,
        p: &BatesParams<T>,
        spot: T,
        rate: T,
        div: T,
        expiry: T,
        strikes: &[T],
        is_call: bool,
    ) -> Vec<T> {
        self.price_impl_extra(
            |u| p.char_func_unit_spot(u, rate, div, expiry),
            p.heston.v0,
            p.heston.theta,
            p.heston.sigma_v,
            p.heston.kappa,
            p.heston.rho,
            spot,
            rate,
            div,
            expiry,
            strikes,
            is_call,
            p.lambda * (p.mu_j * p.mu_j + p.sigma_j * p.sigma_j) * expiry,
            p.lambda
                * (p.mu_j + T::one() - (p.mu_j + p.sigma_j * p.sigma_j * 0.5).exp())
                * expiry,
        )
    }
    #[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
    pub(crate) fn price_impl<T: Real, F: Fn(Complex<T>) -> Complex<T>>(
        &self,
        cf: F,
        v0: T,
        theta: T,
        sigma_v: T,
        kappa: T,
        rho: T,
        spot: T,
        rate: T,
        div: T,
        expiry: T,
        strikes: &[T],
        is_call: bool,
    ) -> Vec<T> {
        self.price_impl_extra(
            cf,
            v0,
            theta,
            sigma_v,
            kappa,
            rho,
            spot,
            rate,
            div,
            expiry,
            strikes,
            is_call,
            T::zero(),
            T::zero(),
        )
    }

    #[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
    pub(crate) fn price_impl_extra<T: Real, F: Fn(Complex<T>) -> Complex<T>>(
        &self,
        cf: F,
        v0: T,
        theta: T,
        sigma_v: T,
        kappa: T,
        rho: T,
        spot: T,
        rate: T,
        div: T,
        expiry: T,
        strikes: &[T],
        is_call: bool,
        jump_variance: T,
        jump_drift: T,
    ) -> Vec<T> {
        let t = expiry.max_val(T::from_f64(1e-8));
        let delta = v0 - theta;
        let mean_var = if kappa.val().abs() < 1e-8 {
            v0 * t - delta * kappa * t * t * 0.5 + delta * kappa * kappa * t * t * t / 6.0
        } else {
            theta * t + delta * (T::one() - (-kappa * t).exp()) / kappa
        };
        let j = |z: T| {
            if z.val().abs() < 1e-7 {
                t + z * t * t * 0.5 + z * z * t * t * t / 6.0
            } else {
                ((z * t).exp() - T::one()) / z
            }
        };
        let (var_i, cov_im) = if kappa.val().abs() < 1e-7 {
            (
                sigma_v * sigma_v * v0 * t * t * t / 3.0,
                rho * sigma_v * v0 * t * t * 0.5,
            )
        } else {
            let jk = j(-kappa);
            let j2k = j(kappa * -2.0);
            let ek = (-kappa * t).exp();
            let a_int = theta * 0.5 * (t - j2k - ek * (j(kappa) - jk))
                + delta * (jk - j2k - ek * (t - jk));
            let vi = sigma_v * sigma_v * 2.0 * a_int / (kappa * kappa);
            let ci = rho * sigma_v / kappa
                * (theta * (t - ek * j(kappa)) + delta * (jk - ek * t));
            (vi, ci)
        };
        let var = (mean_var + var_i * 0.25 - cov_im + jump_variance)
            .max_val(T::from_f64(1e-12));
        let c1 = (rate - div) * t - mean_var * 0.5 + jump_drift;
        let w = T::from_f64(self.l_trunc) * var.sqrt();
        let a = c1 - w;
        let b = c1 + w;
        let width = b - a;
        let mut ph = Vec::with_capacity(self.n_terms);
        for k in 0..self.n_terms {
            ph.push(cf(Complex::from_real(
                T::from_f64(k as f64 * std::f64::consts::PI) / width,
            )));
        }
        strikes
            .iter()
            .map(|&k| {
                let x0 = (spot / k).ln();
                let mut put = T::zero();
                for n in 0..self.n_terms {
                    let om = T::from_f64(n as f64 * std::f64::consts::PI) / width;
                    let psi = if n == 0 { -a } else { (-om * a).sin() / om };
                    let chi =
                        ((-om * a).cos() - a.exp() + om * (-om * a).sin()) / (T::one() + om * om);
                    let vk = T::from_f64(2.0) / width * (-chi + psi);
                    let phase = Complex::new(T::zero(), om * (x0 - a)).exp();
                    let term = (ph[n] * phase).re * vk;
                    put += if n == 0 { term * 0.5 } else { term };
                }
                let put = k * (-rate * expiry).exp() * put;
                if is_call {
                    put + spot * (-div * expiry).exp() - k * (-rate * expiry).exp()
                } else {
                    put
                }
            })
            .collect()
    }
}
