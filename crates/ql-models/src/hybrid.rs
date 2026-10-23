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
            eta_r * eta_r / (k * k * k * 2.0)
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

fn go_sqrt_moment_params<T: Real>(p: &HestonHullWhiteParams<T>) -> (T, T, T) {
    let theta = p.theta_v.max_val(T::from_f64(1e-14));
    let k = p.kappa_v.max_val(T::from_f64(1e-12));
    let stationary_var = theta * p.xi_v * p.xi_v / (k * 2.0);
    let raw = theta - stationary_var / (theta * 8.0);
    let smooth_floor = T::from_f64(0.5) * (raw + (raw * raw + T::from_f64(1e-12)).sqrt());
    let a = smooth_floor.sqrt();
    let b = p.v0.max_val(T::zero()).sqrt() - a;
    // Matching the initial and stationary square-root moments with the CIR
    // mean-reversion scale gives the exponential moment used by
    // Grzelak--Oosterlee in the stock/rate covariance integral.
    (a, b, k)
}

fn hw_cross_integral<T: Real>(kappa_r: T, c: T, t: T) -> T {
    if kappa_r.val().abs() < 1e-6 {
        if c.val().abs() < 1e-6 {
            t * t * 0.5 - c * t * t * t / 6.0
        } else {
            (c * t - T::one() + (-c * t).exp()) / (c * c)
        }
    } else {
        let first = if c.val().abs() < 1e-7 {
            t - c * t * t * 0.5 + c * c * t * t * t / 6.0
        } else {
            (T::one() - (-c * t).exp()) / c
        };
        let d = kappa_r - c;
        let second = if d.val().abs() < 1e-7 {
            t * (-c * t).exp() * (T::one() - d * t * 0.5 + d * d * t * t / 6.0)
        } else {
            ((-c * t).exp() - (-kappa_r * t).exp()) / d
        };
        (first - second) / kappa_r
    }
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
    // Heston's API receives a real Fourier frequency and inserts the
    // imaginary unit internally.  Passing i*u here would apply i twice.
    let z = hp.char_func_unit_spot(Complex::from_real(T::from_f64(u)), rate, div, expiry);
    // Under the forward measure the HW factor contributes the martingale
    // preserving Gaussian factor exp(-½Σ²(u²+iu)).  The cross covariance is
    // evaluated with the Heston effective spot volatility.
    let t = expiry.max_val(T::zero());
    let hw_var = AnalyticBsmHullWhiteEngine::effective_variance(
        T::zero(),
        p.kappa_r,
        p.eta_r,
        T::zero(),
        expiry,
    );
    // Grzelak--Oosterlee moment matching: approximate E[sqrt(v_t)] by
    // a_GO+b_GO exp(-kappa_v t), then integrate the OU covariance exactly.
    let (a_go, b_go, c_go) = go_sqrt_moment_params(p);
    let cross = if p.kappa_r.val().abs() < 1e-8 {
        p.rho_sr * p.eta_r * (a_go * t * t + b_go * t * t * 0.5)
    } else {
        let kr = p.kappa_r;
        let int0 = (kr * t - T::one() + (-kr * t).exp()) / (kr * kr);
        let int1 = hw_cross_integral(kr, c_go, t);
        p.rho_sr * p.eta_r * 2.0 * (a_go * int0 + b_go * int1)
    };
    let hyb = hw_var + cross;
    let uu = T::from_f64(u);
    // Construct in complex form to preserve both real and imaginary parts;
    // the factor is exp(-½Σ²(u²+i u)).
    let q = Complex::new(uu * uu, uu);
    let fac = (-q * Complex::from_real(hyb * 0.5)).exp();
    let out = z * fac;
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
        let extra = AnalyticBsmHullWhiteEngine::effective_variance(
            T::zero(),
            p.kappa_r,
            p.eta_r,
            p.rho_sr,
            expiry,
        );
        let t = expiry.max_val(T::zero());
        let mean_var = if p.kappa_v.val().abs() < 1e-8 {
            p.v0 * t - (p.v0 - p.theta_v) * p.kappa_v * t * t * 0.5
                + (p.v0 - p.theta_v) * p.kappa_v * p.kappa_v * t * t * t / 6.0
        } else {
            p.theta_v * t + (p.v0 - p.theta_v) * (T::one() - (-p.kappa_v * t).exp()) / p.kappa_v
        };
        let (a_go, b_go, c_go) = go_sqrt_moment_params(p);
        let cross = if p.kappa_r.val().abs() < 1e-8 {
            p.rho_sr * p.eta_r * (a_go * t * t + b_go * t * t * 0.5)
        } else {
            let kr = p.kappa_r;
            let int0 = (kr * t - T::one() + (-kr * t).exp()) / (kr * kr);
            let int1 = hw_cross_integral(kr, c_go, t);
            p.rho_sr * p.eta_r * 2.0 * (a_go * int0 + b_go * int1)
        };
        let hyb = (extra + cross).max_val(-mean_var * 0.9);
        CosHestonEngine {
            n_terms: self.n_terms,
            l_trunc: self.l_trunc,
        }
        .price_impl_extra(
            |u| {
                let z = h.char_func_unit_spot(u, rate, div, expiry);
                let q = u * u + Complex::i() * u;
                z * (-(q) * Complex::from_real(hyb * 0.5)).exp()
            },
            p.v0,
            p.theta_v,
            p.xi_v,
            p.kappa_v,
            p.rho_sv,
            spot,
            rate,
            div,
            expiry,
            std::slice::from_ref(&strike),
            is_call,
            hyb,
            -hyb * 0.5,
        )[0]
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
