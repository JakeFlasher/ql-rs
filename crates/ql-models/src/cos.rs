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
        self.price_impl(
            |u| p.char_func_unit_spot(u, rate, div, expiry),
            p.heston.v0,
            p.heston.theta,
            p.heston.sigma_v,
            p.heston.kappa,
            spot,
            rate,
            div,
            expiry,
            strikes,
            is_call,
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
        spot: T,
        rate: T,
        div: T,
        expiry: T,
        strikes: &[T],
        is_call: bool,
    ) -> Vec<T> {
        let t = expiry.val().max(1e-8);
        let c1 = (rate.val() - div.val() - 0.5 * theta.val()) * t
            + (1. - (-kappa.val() * t).exp()) * (theta.val() - v0.val())
                / (2. * kappa.val().max(1e-8));
        let var = (0.5 * (v0.val() + theta.val()) + 0.1 * sigma_v.val() * sigma_v.val()) * t;
        let w = self.l_trunc * var.max(1e-8).sqrt();
        let a = c1 - w;
        let b = c1 + w;
        let width = b - a;
        let mut ph = Vec::with_capacity(self.n_terms);
        for k in 0..self.n_terms {
            ph.push(cf(Complex::from_real(T::from_f64(
                k as f64 * std::f64::consts::PI / width,
            ))));
        }
        strikes
            .iter()
            .map(|&k| {
                let x0 = (spot / k).ln();
                let mut put = T::zero();
                for n in 0..self.n_terms {
                    let om = n as f64 * std::f64::consts::PI / width;
                    let psi = if n == 0 { -a } else { (-om * a).sin() / om };
                    let chi = ((-om * a).cos() - a.exp() + om * (-om * a).sin()) / (1. + om * om);
                    let vk = T::from_f64(2. / width * (-chi + psi));
                    let phase =
                        Complex::new(T::zero(), T::from_f64(om) * (x0 - T::from_f64(a))).exp();
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
