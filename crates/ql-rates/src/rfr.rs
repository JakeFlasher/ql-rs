//! Backward-looking risk-free-rate (RFR) caplets and Willems effective SABR.

use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_vol::black76_price;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum RfrRateType {
    ForwardLooking,
    BackwardLooking,
}

#[derive(Copy, Clone, Debug)]
pub struct RfrCapletSpec<T: Real> {
    pub notional: T,
    pub forward_rate: T,
    pub strike: T,
    pub tau_0: T,
    pub tau_1: T,
    pub accrual_fraction: T,
    pub discount_tau_1: T,
    pub is_caplet: bool,
    pub rate_type: RfrRateType,
}

pub fn rfr_effective_variance_time<T: Real>(
    tau_0: T,
    tau_1: T,
    q_decay: T,
    rate_type: RfrRateType,
) -> T {
    if rate_type == RfrRateType::ForwardLooking {
        return tau_0.max_val(T::zero());
    }
    let start = tau_0.max_val(T::zero());
    let delta = (tau_1 - tau_0).max_val(T::from_f64(1e-14));
    let q = q_decay.max_val(T::zero());
    start
        + (tau_1 - start).max_val(T::zero()).pow(q * 2.0 + 1.0)
            / ((q * 2.0 + 1.0) * delta.pow(q * 2.0))
}

fn bachelier_unit<T: Real>(f: T, k: T, sigma: T, expiry: T, is_call: bool) -> T {
    let st = sigma.abs() * expiry.max_val(T::zero()).sqrt();
    let sign = if is_call { 1.0 } else { -1.0 };
    if st.val() <= 1e-14 {
        return (T::from_f64(sign) * (f - k)).max_val(T::zero());
    }
    let d = (f - k) / st;
    T::from_f64(sign) * (f - k) * (T::from_f64(sign) * d).norm_cdf() + st * d.norm_pdf()
}

pub fn price_rfr_caplet_black<T: Real>(
    spec: &RfrCapletSpec<T>,
    sigma: T,
    shift: T,
    q_decay: T,
) -> T {
    let te = rfr_effective_variance_time(spec.tau_0, spec.tau_1, q_decay, spec.rate_type);
    let unit = black76_price(
        spec.forward_rate + shift,
        spec.strike + shift,
        sigma,
        te,
        spec.discount_tau_1,
        spec.is_caplet,
    );
    spec.notional * spec.accrual_fraction * unit
}

pub fn price_rfr_caplet_bachelier<T: Real>(
    spec: &RfrCapletSpec<T>,
    normal_sigma: T,
    q_decay: T,
) -> T {
    let te = rfr_effective_variance_time(spec.tau_0, spec.tau_1, q_decay, spec.rate_type);
    spec.notional
        * spec.accrual_fraction
        * spec.discount_tau_1
        * bachelier_unit(
            spec.forward_rate,
            spec.strike,
            normal_sigma,
            te,
            spec.is_caplet,
        )
}

#[derive(Copy, Clone, Debug)]
pub struct EffectiveSabrParams<T: Real> {
    pub alpha_hat: T,
    pub beta: T,
    pub rho_hat: T,
    pub nu_hat: T,
    pub exercise_time: T,
}

pub fn willems_backward_sabr_params<T: Real>(
    tau_0: T,
    tau_1: T,
    alpha: T,
    beta: T,
    rho: T,
    nu: T,
    q: T,
) -> EffectiveSabrParams<T> {
    let q = q.max_val(T::zero());
    // The two closed forms meet at the accrual boundary.  Treat an
    // infinitesimal negative time as the boundary itself to avoid numerical
    // branch noise when a date is represented by an AD perturbation.
    if tau_0.val() >= 0.0 {
        let tau = q * tau_0 * 2.0 + tau_1;
        let gamma = tau
            * (tau.powf(3.0) * 2.0
                + tau_1.powf(3.0)
                + (q * q * 4.0 - q * 2.0) * tau_0.powf(3.0)
                + q * 6.0 * tau_0 * tau_0 * tau_1)
            / ((q * 4.0 + 3.0) * (q * 2.0 + 1.0))
            + q * rho
                * rho
                * 3.0
                * (tau_1 - tau_0).powf(2.0)
                * (tau * tau * 3.0 - tau_1 * tau_1 + q * 5.0 * tau_0 * tau_0 + tau_0 * tau_1 * 4.0)
                / ((q * 4.0 + 3.0) * (q * 3.0 + 2.0).powf(2.0));
        let rho_hat = rho * (tau * tau * 3.0 + q * 2.0 * tau_0 * tau_0 + tau_1 * tau_1)
            / (gamma.sqrt() * (q * 6.0 + 4.0));
        let nu_hat_sq = nu * nu * gamma * (q * 2.0 + 1.0) / (tau.powf(3.0) * tau_1);
        let h = nu * nu * (tau * tau + q * 2.0 * tau_0 * tau_0 + tau_1 * tau_1)
            / (tau_1 * tau * (q + 1.0) * 2.0)
            - nu_hat_sq;
        let alpha_hat = alpha * (tau / ((q * 2.0 + 1.0) * tau_1)).sqrt() * (h * tau_1 * 0.25).exp();
        EffectiveSabrParams {
            alpha_hat,
            beta,
            rho_hat,
            nu_hat: nu_hat_sq.max_val(T::zero()).sqrt(),
            exercise_time: tau_1,
        }
    } else {
        let zeta = T::from_f64(3.0) / (q * 4.0 + 3.0)
            * (T::one() / (q * 2.0 + 1.0) + rho * rho * q * 2.0 / (q * 3.0 + 2.0).powf(2.0));
        let rho_hat = rho * 2.0 / (zeta.sqrt() * (q * 3.0 + 2.0));
        let nu_hat = nu * (zeta * (q * 2.0 + 1.0)).sqrt();
        let alpha_hat = alpha / (q * 2.0 + 1.0).sqrt()
            * (tau_1 / (tau_1 - tau_0)).pow(q)
            * ((nu * nu / (q + 1.0) - nu_hat * nu_hat) * tau_1 * 0.25).exp();
        EffectiveSabrParams {
            alpha_hat,
            beta,
            rho_hat,
            nu_hat,
            exercise_time: tau_1,
        }
    }
}

fn hagan_black_vol<T: Real>(f: T, k: T, expiry: T, alpha: T, beta: T, rho: T, nu: T) -> T {
    let one_minus = T::one() - beta;
    let fk = f * k;
    let log_fk = (f / k).ln();
    let z = nu / alpha * fk.pow(one_minus * 0.5) * log_fk;
    let z_over_x = if z.val().abs() < 1e-8 {
        T::one() - rho * z * 0.5 + (T::from_f64(2.0) - rho * rho * 3.0) * z * z / 12.0
    } else {
        let x = ((T::one() - rho * z * 2.0 + z * z).sqrt() + z - rho) / (T::one() - rho);
        z / x.ln()
    };
    let corr = T::one()
        + expiry
            * (one_minus * one_minus * alpha * alpha / (fk.pow(one_minus) * 24.0)
                + rho * beta * nu * alpha / (fk.pow(one_minus * 0.5) * 4.0)
                + (T::from_f64(2.0) - rho * rho * 3.0) * nu * nu / 24.0);
    alpha * z_over_x * corr / fk.pow(one_minus * 0.5)
}

#[allow(clippy::too_many_arguments)]
pub fn willems_backward_sabr_implied_vol<T: Real>(
    forward: T,
    strike: T,
    tau_0: T,
    tau_1: T,
    alpha: T,
    beta: T,
    rho: T,
    nu: T,
    q: T,
) -> T {
    let p = willems_backward_sabr_params(tau_0, tau_1, alpha, beta, rho, nu, q);
    hagan_black_vol(
        forward,
        strike,
        p.exercise_time,
        p.alpha_hat,
        p.beta,
        p.rho_hat,
        p.nu_hat,
    )
}

pub fn willems_backward_sabr_caplet_price<T: Real>(
    spec: &RfrCapletSpec<T>,
    alpha: T,
    beta: T,
    rho: T,
    nu: T,
    q: T,
) -> T {
    if spec.rate_type == RfrRateType::BackwardLooking {
        let p = willems_backward_sabr_params(spec.tau_0, spec.tau_1, alpha, beta, rho, nu, q);
        let vol = hagan_black_vol(
            spec.forward_rate,
            spec.strike,
            p.exercise_time,
            p.alpha_hat,
            p.beta,
            p.rho_hat,
            p.nu_hat,
        );
        spec.notional
            * spec.accrual_fraction
            * black76_price(
                spec.forward_rate,
                spec.strike,
                vol,
                p.exercise_time,
                spec.discount_tau_1,
                spec.is_caplet,
            )
    } else {
        let t0 = spec.tau_0.max_val(T::zero());
        let vol = hagan_black_vol(spec.forward_rate, spec.strike, t0, alpha, beta, rho, nu);
        spec.notional
            * spec.accrual_fraction
            * black76_price(
                spec.forward_rate,
                spec.strike,
                vol,
                t0,
                spec.discount_tau_1,
                spec.is_caplet,
            )
    }
}

pub fn willems_backward_sabr_caplet_strip<T: Real>(
    spec: &RfrCapletSpec<T>,
    strikes: &[T],
    alpha: T,
    beta: T,
    rho: T,
    nu: T,
    q: T,
) -> Vec<T> {
    if spec.rate_type == RfrRateType::BackwardLooking {
        let p = willems_backward_sabr_params(spec.tau_0, spec.tau_1, alpha, beta, rho, nu, q);
        strikes
            .iter()
            .map(|&k| {
                let vol = hagan_black_vol(
                    spec.forward_rate,
                    k,
                    p.exercise_time,
                    p.alpha_hat,
                    p.beta,
                    p.rho_hat,
                    p.nu_hat,
                );
                spec.notional
                    * spec.accrual_fraction
                    * black76_price(
                        spec.forward_rate,
                        k,
                        vol,
                        p.exercise_time,
                        spec.discount_tau_1,
                        spec.is_caplet,
                    )
            })
            .collect()
    } else {
        let t0 = spec.tau_0.max_val(T::zero());
        strikes
            .iter()
            .map(|&k| {
                let vol = hagan_black_vol(spec.forward_rate, k, t0, alpha, beta, rho, nu);
                spec.notional
                    * spec.accrual_fraction
                    * black76_price(
                        spec.forward_rate,
                        k,
                        vol,
                        t0,
                        spec.discount_tau_1,
                        spec.is_caplet,
                    )
            })
            .collect()
    }
}

fn hw_b<T: Real>(a: T, dt: T) -> T {
    if a.val().abs() < 1e-10 {
        dt
    } else {
        (T::one() - (-a * dt).exp()) / a
    }
}

#[allow(clippy::too_many_arguments)]
pub fn hw1f_rfr_caplet_price<T: Real>(
    curve: &impl YieldTermStructure<T>,
    a: T,
    sigma: T,
    tau_0: f64,
    tau_1: f64,
    strike: T,
    notional: T,
    is_caplet: bool,
    rate_type: RfrRateType,
) -> T {
    let t0 = T::from_f64(tau_0.max(0.0));
    let t1 = T::from_f64(tau_1);
    let dt = t1 - T::from_f64(tau_0);
    let b = hw_b(a, dt);
    let srr = if a.val().abs() < 1e-10 {
        sigma * sigma * t0
    } else {
        sigma * sigma / (a * 2.0) * (T::one() - (-a * t0 * 2.0).exp())
    };
    let szz = if a.val().abs() < 1e-10 {
        sigma * sigma * dt.powf(3.0) / 3.0
    } else {
        sigma * sigma / (a * a)
            * (dt - (T::one() - (-a * dt).exp()) * 2.0 / a
                + (T::one() - (-a * dt * 2.0).exp()) / (a * 2.0))
    };
    let variance = match rate_type {
        RfrRateType::ForwardLooking => b * b * srr,
        RfrRateType::BackwardLooking => b * b * srr + szz,
    };
    let p0 = curve.discount_factor(tau_0.max(0.0));
    let p1 = curve.discount_factor(tau_1);
    let k_inv = T::one() + dt * strike;
    let st = variance.max_val(T::zero()).sqrt();
    if st.val() < 1e-14 {
        return notional
            * (if is_caplet {
                (p0 - k_inv * p1).max_val(T::zero())
            } else {
                (k_inv * p1 - p0).max_val(T::zero())
            });
    }
    let d1 = ((p0 / (k_inv * p1)).ln() + variance * 0.5) / st;
    let d2 = d1 - st;
    notional
        * if is_caplet {
            p0 * d1.norm_cdf() - k_inv * p1 * d2.norm_cdf()
        } else {
            k_inv * p1 * (-d2).norm_cdf() - p0 * (-d1).norm_cdf()
        }
}

pub fn hw1f_sofr_futures_rate_3m<T: Real>(
    curve: &impl YieldTermStructure<T>,
    a: T,
    sigma: T,
    tau_0: f64,
    tau_1: f64,
) -> T {
    let t0 = T::from_f64(tau_0.max(0.0));
    let dt = T::from_f64(tau_1 - tau_0);
    let b = hw_b(a, dt);
    let srr = if a.val().abs() < 1e-10 {
        sigma * sigma * t0
    } else {
        sigma * sigma / (a * 2.0) * (T::one() - (-a * t0 * 2.0).exp())
    };
    let srz = if a.val().abs() < 1e-10 {
        sigma * sigma * t0 * t0 * 0.5
    } else {
        sigma * sigma / (a * a) * (T::one() - (-a * t0).exp()).powf(2.0) * 0.5
    };
    let szz = if a.val().abs() < 1e-10 {
        sigma * sigma * dt.powf(3.0) / 3.0
    } else {
        sigma * sigma / (a * a)
            * (dt - (T::one() - (-a * dt).exp()) * 2.0 / a
                + (T::one() - (-a * dt * 2.0).exp()) / (a * 2.0))
    };
    let conv = b * srz + b * b * srr * 0.5 + szz;
    ((curve.discount_factor(tau_0.max(0.0)) / curve.discount_factor(tau_1)) * conv.exp() - T::one())
        / dt
}
