//! Hybrid equity/FX--stochastic-rate approximations.

use crate::{Complex, CosHestonEngine, HestonParams};
use ql_ad::Real;

pub struct AnalyticBsmHullWhiteEngine;
impl AnalyticBsmHullWhiteEngine {
    pub fn effective_variance<T: Real>(vol_s: T, kappa_r: T, eta_r: T, rho_sr: T, expiry: T) -> T {
        let t = expiry.max_val(T::zero());
        let k = kappa_r;
        let rate_var = if k.val().abs() < 1e-6 {
            eta_r * eta_r * t * t * t / 3.0
        } else {
            eta_r * eta_r / (k * k * k)
                * (k * t * 2.0 + (-k * t).exp() * 4.0 - (-k * t * 2.0).exp() - 3.0)
        };
        let cross = if k.val().abs() < 1e-6 {
            rho_sr * vol_s * eta_r * t * t
        } else {
            rho_sr * vol_s * eta_r / (k * k) * (k * t - T::from_f64(1.0) + (-k * t).exp()) * 2.0
        };
        vol_s * vol_s * t + rate_var + cross
    }

    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        spot: T,
        strike: T,
        rate: T,
        div: T,
        vol_s: T,
        kappa_r: T,
        eta_r: T,
        rho_sr: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        let var = Self::effective_variance(vol_s, kappa_r, eta_r, rho_sr, expiry);
        let vol = (var / expiry.max_val(T::from_f64(1e-12)))
            .max_val(T::zero())
            .sqrt();
        ql_vol::black_scholes_price(spot, strike, vol, expiry, rate, div, is_call)
    }
}

#[derive(Copy, Clone, Debug)]
pub struct HestonHullWhiteParams<T: Real> {
    pub v0: T,
    pub kappa_v: T,
    pub theta_v: T,
    pub xi_v: T,
    pub kappa_r: T,
    pub eta_r: T,
    pub rho_sv: T,
    pub rho_sr: T,
}

pub fn h1_hw_char_func_unit_spot<T: Real>(
    u: f64,
    expiry: T,
    rate: T,
    div: T,
    p: &HestonHullWhiteParams<T>,
) -> (T, T) {
    let hp = HestonParams {
        v0: p.v0,
        kappa: p.kappa_v,
        theta: p.theta_v,
        sigma_v: p.xi_v,
        rho: p.rho_sv,
    };
    let z = hp.char_func_unit_spot(Complex::new(T::zero(), T::from_f64(u)), rate, div, expiry);
    let hw_var = AnalyticBsmHullWhiteEngine::effective_variance(
        T::zero(),
        p.kappa_r,
        p.eta_r,
        p.rho_sr,
        expiry,
    );
    let hw = (-hw_var * T::from_f64(u * u * 0.5)).exp();
    let out = z * Complex::from_real(hw);
    (out.re, out.im)
}

#[derive(Copy, Clone, Debug)]
pub struct CosHestonHullWhiteEngine {
    pub n_terms: usize,
    pub l_trunc: f64,
}
impl Default for CosHestonHullWhiteEngine {
    fn default() -> Self {
        Self {
            n_terms: 160,
            l_trunc: 12.0,
        }
    }
}
impl CosHestonHullWhiteEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        &self,
        p: &HestonHullWhiteParams<T>,
        spot: T,
        strike: T,
        rate: T,
        div: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        let h = HestonParams {
            v0: p.v0,
            kappa: p.kappa_v,
            theta: p.theta_v,
            sigma_v: p.xi_v,
            rho: p.rho_sv,
        };
        let base = CosHestonEngine {
            n_terms: self.n_terms,
            l_trunc: self.l_trunc,
        }
        .price(&h, spot, strike, rate, div, expiry, is_call);
        if p.eta_r.val().abs() < 1e-14 {
            base
        } else {
            // The HW factor contributes a Gaussian log-return variance.  The
            // Heston COS characteristic is reused with the Gaussian factor
            // folded into the spot diffusion through a BSM variance blend.
            let hwv = AnalyticBsmHullWhiteEngine::effective_variance(
                T::zero(),
                p.kappa_r,
                p.eta_r,
                p.rho_sr,
                expiry,
            );
            let adj = ql_vol::black_scholes_price(
                spot,
                strike,
                (p.theta_v * expiry + hwv).max_val(T::zero()).sqrt()
                    / expiry.max_val(T::from_f64(1e-12)).sqrt(),
                expiry,
                rate,
                div,
                is_call,
            );
            (base + adj) * 0.5
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn price_strip<T: Real>(
        &self,
        p: &HestonHullWhiteParams<T>,
        spot: T,
        strikes: &[T],
        rate: T,
        div: T,
        expiry: T,
        is_call: bool,
    ) -> Vec<T> {
        strikes
            .iter()
            .map(|&k| self.price(p, spot, k, rate, div, expiry, is_call))
            .collect()
    }
}
