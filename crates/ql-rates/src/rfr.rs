//! Backward-looking risk-free-rate (RFR) caplets and Willems effective SABR.
#![allow(clippy::excessive_precision)]

use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_vol::black76_price;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum RfrRateType {
    ForwardLooking,
    BackwardLooking,
}

impl RfrRateType {
    /// Alias used by compounded-in-arrears APIs.
    #[allow(non_upper_case_globals)]
    pub const BackwardCompounded: Self = Self::BackwardLooking;
}

/// Internal adapter allowing legacy and compounded-call argument ordering.
pub enum HwRfrArg { Bool(bool), Rate(RfrRateType) }
impl From<bool> for HwRfrArg { fn from(v: bool) -> Self { Self::Bool(v) } }
impl From<RfrRateType> for HwRfrArg { fn from(v: RfrRateType) -> Self { Self::Rate(v) } }

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
    arg7: impl Into<HwRfrArg>,
    arg8: impl Into<HwRfrArg>,
) -> T {
    let (is_caplet, rate_type) = match (arg7.into(), arg8.into()) {
        (HwRfrArg::Bool(b), HwRfrArg::Rate(r)) => (b, r),
        (HwRfrArg::Rate(r), HwRfrArg::Bool(b)) => (b, r),
        _ => (true, RfrRateType::BackwardLooking),
    };
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

#[derive(Clone, Debug)]
pub struct TenorBasisSwapSpec<T: Real = f64>{pub notional:T,pub short_leg_times:Vec<f64>,pub long_leg_times:Vec<f64>,pub short_leg_basis_spread:T,pub pay_short_leg:bool}
#[derive(Copy,Clone,Debug)]
pub struct TenorBasisValuation<T: Real>{pub npv:T,pub fair_short_basis_spread:T,pub short_leg_bpv:T,pub long_leg_pv:T}
pub fn price_tenor_basis_swap<T: Real>(s:&TenorBasisSwapSpec<T>, d:&impl YieldTermStructure<T>, sh:&impl YieldTermStructure<T>, lo:&impl YieldTermStructure<T>)->TenorBasisValuation<T>{
 let mut sf=T::zero(); let mut bpv=T::zero(); let mut prev=0.0; for &t in &s.short_leg_times{if t<=prev+1e-8{continue} let tau=t-prev; let p0=if prev<=1e-12{T::one()}else{sh.discount_factor(prev)}; let p1=sh.discount_factor(t).max_val(T::from_f64(1e-12)); let df=d.discount_factor(t); sf += s.notional*(p0/p1-T::one())*df; bpv += s.notional*df*tau; prev=t;} let mut lp=T::zero(); prev=0.0; for &t in &s.long_leg_times{if t<=prev+1e-8{continue} let p0=if prev<=1e-12{T::one()}else{lo.discount_factor(prev)}; let p1=lo.discount_factor(t).max_val(T::from_f64(1e-12)); lp += s.notional*(p0/p1-T::one())*d.discount_factor(t); prev=t;} let fair=(lp-sf)/bpv.max_val(T::from_f64(1e-12)); let st=sf+s.short_leg_basis_spread*bpv; let npv=if s.pay_short_leg{lp-st}else{st-lp}; TenorBasisValuation{npv,fair_short_basis_spread:fair,short_leg_bpv:bpv,long_leg_pv:lp}
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;

    #[test]
    fn tenor_basis_fair_spread_and_ad() {
        let d = Date(0);
        let ois = FlatForwardCurve::new(d, 0.03);
        let short_c = FlatForwardCurve::new(d, 0.032);
        let long_c = FlatForwardCurve::new(d, 0.035);
        let mut s = TenorBasisSwapSpec {
            notional: 1_000_000.0,
            short_leg_times: vec![0.25, 0.50, 0.75, 1.0],
            long_leg_times: vec![0.50, 1.0],
            short_leg_basis_spread: 0.0,
            pay_short_leg: true,
        };
        let v0 = price_tenor_basis_swap(&s, &ois, &short_c, &long_c);
        assert!(v0.fair_short_basis_spread > 0.0);
        s.short_leg_basis_spread = v0.fair_short_basis_spread;
        let v_par = price_tenor_basis_swap(&s, &ois, &short_c, &long_c);
        assert!(v_par.npv.abs() < 1e-9);

        let h = 1e-5;
        let eval_r = |r_long: f64| {
            let lc = FlatForwardCurve::new(d, r_long);
            price_tenor_basis_swap(&s, &ois, &short_c, &lc).fair_short_basis_spread
        };
        let fd = (eval_r(0.035 + h) - eval_r(0.035 - h)) / (2.0 * h);
        type D = Dual<1>;
        let ois_d = FlatForwardCurve::new(d, D::from_value(0.03));
        let sh_d = FlatForwardCurve::new(d, D::from_value(0.032));
        let lo_d = FlatForwardCurve::new(d, D::variable(0.035, 0));
        let sd = TenorBasisSwapSpec {
            notional: D::from_value(1_000_000.0),
            short_leg_times: vec![0.25, 0.50, 0.75, 1.0],
            long_leg_times: vec![0.50, 1.0],
            short_leg_basis_spread: D::zero(),
            pay_short_leg: true,
        };
        let ad = price_tenor_basis_swap(&sd, &ois_d, &sh_d, &lo_d).fair_short_basis_spread;
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}


#[derive(Clone, Debug)]
pub struct SinhOuRfrCapletSpec<T: Real = f64> {
    pub notional: T,
    pub start_time: f64,
    pub end_time: f64,
    pub strike: T,
    pub alpha_mean_reversion: T,
    pub sigma_short_rate: T,
    pub gamma_smile: T,
    pub y_star_skew: T,
    pub is_caplet: bool,
    pub is_backward_compounded: bool,
}

#[derive(Copy, Clone, Debug)]
pub struct SinhOuRfrCapletValuation<T: Real = f64> {
    pub npv: T,
    pub base_variance: T,
    pub effective_variance: T,
    pub smile_skew_variance_adj: T,
    pub implied_hw_short_rate_vol: T,
    pub c1_shift: T,
    pub c2_skew: T,
    pub c3_smile: T,
}

fn sinh_ou_phi_r<T: Real>(alpha: T, t: f64, v: f64) -> T {
    let tau = (v - t).max(0.0);
    (-alpha * tau).exp()
}

fn sinh_ou_sigma_rr<T: Real>(alpha: T, sigma: T, t: f64, v: f64) -> T {
    let tau = (v - t).max(0.0);
    let s2 = sigma * sigma;
    if alpha.val().abs() < 1e-8 {
        s2 * tau * (T::one() - alpha * tau)
    } else {
        s2 * (T::one() - (-alpha * (2.0 * tau)).exp()) / (alpha * 2.0)
    }
}

fn sinh_ou_psi_r<T: Real>(alpha: T, sigma: T, gamma: T, t: f64, v: f64) -> T {
    (gamma * gamma * sinh_ou_sigma_rr(alpha, sigma, t, v) * 0.5).exp()
}

const GL4_NODES: [f64; 4] = [
    -0.861_136_311_594_052_6,
    -0.339_981_043_584_856_3,
    0.339_981_043_584_856_3,
    0.861_136_311_594_052_6,
];
const GL4_WEIGHTS: [f64; 4] = [
    0.347_854_845_137_453_8,
    0.652_145_154_862_546_1,
    0.652_145_154_862_546_1,
    0.347_854_845_137_453_8,
];

fn sinh_ou_sigma_rz<T: Real>(alpha: T, sigma: T, gamma: T, t: f64, v: f64) -> T {
    if v <= t + 1e-12 {
        return T::zero();
    }
    let mid = 0.5 * (t + v);
    let half = 0.5 * (v - t);
    let mut sum = T::zero();
    for m in 0..4 {
        let u = mid + half * GL4_NODES[m];
        let integrand = sinh_ou_psi_r(alpha, sigma, gamma, t, u)
            * sinh_ou_phi_r(alpha, u, v)
            * sinh_ou_sigma_rr(alpha, sigma, t, u);
        sum += integrand * (half * GL4_WEIGHTS[m]);
    }
    sum
}

fn sinh_ou_sigma_zz<T: Real>(alpha: T, sigma: T, gamma: T, t: f64, v: f64) -> T {
    if v <= t + 1e-12 {
        return T::zero();
    }
    let mid = 0.5 * (t + v);
    let half = 0.5 * (v - t);
    let mut sum = T::zero();
    for m in 0..4 {
        let u = mid + half * GL4_NODES[m];
        let integrand =
            sinh_ou_psi_r(alpha, sigma, gamma, t, u) * sinh_ou_sigma_rz(alpha, sigma, gamma, t, u);
        sum += integrand * (half * GL4_WEIGHTS[m]);
    }
    sum * 2.0
}

fn sinh_ou_b_star<T: Real>(alpha: T, sigma: T, gamma: T, t: f64, v: f64) -> T {
    if v <= t + 1e-12 {
        return T::zero();
    }
    let mid = 0.5 * (t + v);
    let half = 0.5 * (v - t);
    let mut sum = T::zero();
    for m in 0..4 {
        let u = mid + half * GL4_NODES[m];
        let integrand = sinh_ou_psi_r(alpha, sigma, gamma, t, u) * sinh_ou_phi_r(alpha, t, u);
        sum += integrand * (half * GL4_WEIGHTS[m]);
    }
    sum
}

fn sinh_ou_b_plus<T: Real>(alpha: T, sigma: T, gamma: T, t: f64, t1: f64, v: f64) -> T {
    if v <= t1 + 1e-12 {
        return T::zero();
    }
    let mid = 0.5 * (t1 + v);
    let half = 0.5 * (v - t1);
    let mut sum = T::zero();
    for m in 0..4 {
        let u = mid + half * GL4_NODES[m];
        let integrand = sinh_ou_psi_r(alpha, sigma, gamma, t, u) * sinh_ou_phi_r(alpha, t1, u);
        sum += integrand * (half * GL4_WEIGHTS[m]);
    }
    sum
}

/// Turfus & Romero-Bermudez (2023) analytic backward-looking compounded RFR & forward-looking
/// term caplet/floorlet pricing with explicit smile ($\gamma$) and skew ($y^*$).
pub fn price_sinh_ou_rfr_caplet<T: Real>(
    spec: &SinhOuRfrCapletSpec<T>,
    curve: &impl YieldTermStructure<T>,
) -> SinhOuRfrCapletValuation<T> {
    let t1 = spec.start_time.max(0.0);
    let t2 = spec.end_time.max(t1 + 1e-8);
    let delta = t2 - t1;
    let alpha = spec.alpha_mean_reversion;
    let sigma = spec.sigma_short_rate;
    let gamma = spec.gamma_smile;
    let y_star = spec.y_star_skew;

    let d1_df = curve.discount_factor(t1).max_val(T::from_f64(1e-12));
    let d2_df = curve.discount_factor(t2).max_val(T::from_f64(1e-12));
    let k_inv = T::one() + spec.strike * delta;

    let b_star_12 = sinh_ou_b_star(alpha, sigma, gamma, t1, t2);
    let sig_rr_01 = sinh_ou_sigma_rr(alpha, sigma, 0.0, t1);
    let sig_zz_12 = sinh_ou_sigma_zz(alpha, sigma, gamma, t1, t2);
    let v_comp = b_star_12 * b_star_12 * sig_rr_01 + sig_zz_12;
    let v_term = (b_star_12 * b_star_12 * sig_rr_01).max_val(T::from_f64(1e-14));
    let v_base = if spec.is_backward_compounded {
        v_comp.max_val(T::from_f64(1e-14))
    } else {
        v_term
    };
    let sqrt_v = v_base.sqrt();

    let mid = 0.5 * (t1 + t2);
    let half = 0.5 * delta;
    let mut c1 = T::zero();
    let mut c2 = T::zero();
    let mut c3 = T::zero();
    for m in 0..4 {
        let u = mid + half * GL4_NODES[m];
        let w = half * GL4_WEIGHTS[m];
        let y_st = gamma
            * (y_star
                - sinh_ou_b_plus(alpha, sigma, gamma, 0.0, u, t2)
                    * sinh_ou_sigma_rr(alpha, sigma, 0.0, u)
                - sinh_ou_sigma_rz(alpha, sigma, gamma, 0.0, u));
        let psi_val = if spec.is_backward_compounded {
            (b_star_12 * sinh_ou_phi_r(alpha, t1, u) * sig_rr_01
                + sinh_ou_sigma_rz(alpha, sigma, gamma, t1, u)
                + sinh_ou_b_plus(alpha, sigma, gamma, t1, u, t2)
                    * sinh_ou_sigma_rr(alpha, sigma, t1, u))
                / sqrt_v
        } else {
            sinh_ou_phi_r(alpha, t1, u) * sig_rr_01.max_val(T::from_f64(1e-16)).sqrt() / sqrt_v
        };
        let psi0 = sinh_ou_psi_r(alpha, sigma, gamma, 0.0, u);
        let psi1 = sinh_ou_psi_r(alpha, sigma, gamma, t1, u);
        let ch = (y_st.exp() + (-y_st).exp()) * 0.5;
        let sh = (y_st.exp() - (-y_st).exp()) * 0.5;
        c1 += (psi0 * ch - psi1) * psi_val * w;
        c2 += psi0 * sh * gamma * psi_val * psi_val * (0.5 * w);
        c3 += psi0 * ch * gamma * gamma * psi_val * psi_val * psi_val * (w / 6.0);
    }

    let dz_star = (k_inv * d2_df / d1_df).max_val(T::from_f64(1e-12)).ln();
    let d2_0 = (-dz_star - v_base * 0.5) / sqrt_v;
    let v_adj = sqrt_v * 2.0 * (c1 - c2 * d2_0 + c3 * (d2_0 * d2_0 - T::one()));
    let v_eff = (v_base + v_adj).max_val(T::from_f64(1e-14));
    let st_eff = v_eff.sqrt();
    let d1 = (-dz_star + v_eff * 0.5) / st_eff;
    let d2 = d1 - st_eff;
    let npv = spec.notional
        * if spec.is_caplet {
            d1_df * d1.norm_cdf() - k_inv * d2_df * d2.norm_cdf()
        } else {
            k_inv * d2_df * (-d2).norm_cdf() - d1_df * (-d1).norm_cdf()
        };

    // Baseline Hull-White variance (gamma = 0) for implied HW short-rate volatility
    let b_hw = sinh_ou_b_star(alpha, sigma, T::zero(), t1, t2);
    let v_hw = if spec.is_backward_compounded {
        (b_hw * b_hw * sig_rr_01 + sinh_ou_sigma_zz(alpha, sigma, T::zero(), t1, t2))
            .max_val(T::from_f64(1e-14))
    } else {
        (b_hw * b_hw * sig_rr_01).max_val(T::from_f64(1e-14))
    };
    let implied_hw_vol = sigma * (v_eff / v_hw).sqrt();

    SinhOuRfrCapletValuation {
        npv,
        base_variance: v_base,
        effective_variance: v_eff,
        smile_skew_variance_adj: v_adj,
        implied_hw_short_rate_vol: implied_hw_vol,
        c1_shift: c1,
        c2_skew: c2,
        c3_smile: c3,
    }
}

#[derive(Clone, Debug)]
pub struct CirppSofrFuturesSpec<T: Real = f64> {
    pub start_time: f64,
    pub end_time: f64,
    pub x0: T,
    pub b_drift: T,
    pub beta_mean_reversion: T,
    pub sigma: T,
}

#[derive(Copy, Clone, Debug)]
pub struct CirppSofrFuturesValuation<T: Real = f64> {
    pub ois_forward_rate: T,
    pub futures_rate_3m_compounded: T,
    pub convexity_adjustment_3m: T,
    pub futures_rate_1m_arithmetic: T,
    pub convexity_adjustment_1m: T,
    pub shift_integral_l: T,
}

fn cir_bond_phi_psi<T: Real>(b: T, beta: T, sigma: T, tau: f64) -> (T, T) {
    if tau <= 1e-12 {
        return (T::zero(), T::zero());
    }
    let s2 = (sigma * sigma).max_val(T::from_f64(1e-14));
    let th1 = (beta * beta + s2 * 2.0).sqrt();
    let e_th = (th1 * tau).exp();
    let denom = th1 * 2.0 + (beta + th1) * (e_th - T::one());
    let num = th1 * 2.0 * ((th1 + beta) * (0.5 * tau)).exp();
    let phi = (b * 2.0 / s2) * (num / denom).max_val(T::from_f64(1e-14)).ln();
    let psi = -(e_th - T::one()) * 2.0 / denom;
    (phi, psi)
}

fn cir_compounded_phi_psi<T: Real>(b: T, beta: T, sigma: T, tau: f64) -> (T, T) {
    if tau <= 1e-12 {
        return (T::zero(), T::zero());
    }
    let s2 = (sigma * sigma).max_val(T::from_f64(1e-14));
    let disc = beta * beta - s2 * 2.0;
    if disc.val() > 1e-10 {
        let th = disc.sqrt();
        let e_th = (th * tau).exp();
        let denom = th * 2.0 + (beta + th) * (e_th - T::one());
        let num = th * 2.0 * ((th + beta) * (0.5 * tau)).exp();
        let phi = (b * 2.0 / s2) * (num / denom).max_val(T::from_f64(1e-14)).ln();
        let psi = (e_th - T::one()) * 2.0 / denom;
        (phi, psi)
    } else if disc.val() < -1e-10 {
        let om = (-disc).sqrt();
        let half_om = om * (0.5 * tau);
        let denom = om * half_om.cos() + beta * half_om.sin();
        let num = om * (beta * (0.5 * tau)).exp();
        let phi = (b * 2.0 / s2) * (num / denom).max_val(T::from_f64(1e-14)).ln();
        let psi = half_om.sin() * 2.0 / denom;
        (phi, psi)
    } else {
        let denom = T::from_f64(2.0) + beta * tau;
        let phi = (b * 2.0 / s2)
            * (T::from_f64(2.0) * (beta * (0.5 * tau)).exp() / denom)
                .max_val(T::from_f64(1e-14))
                .ln();
        let psi = T::from_f64(2.0 * tau) / denom;
        (phi, psi)
    }
}

/// Fontana (2023) CIR++ exact affine 1M (arithmetic average) & 3M (compounded in-arrears)
/// SOFR futures rates and futures-forward convexity adjustments (Proposition 2.5, Eq. 6.1, Prop. 6.2-6.3).
pub fn price_cirpp_sofr_futures<T: Real>(
    spec: &CirppSofrFuturesSpec<T>,
    curve: &impl YieldTermStructure<T>,
) -> CirppSofrFuturesValuation<T> {
    let s = spec.start_time.max(0.0);
    let t = spec.end_time.max(s + 1e-8);
    let dt = t - s;
    let b = spec.b_drift;
    let beta = spec.beta_mean_reversion.max_val(T::from_f64(1e-8));
    let sigma = spec.sigma;
    let x0 = spec.x0;

    let ps = curve.discount_factor(s).max_val(T::from_f64(1e-12));
    let pt = curve.discount_factor(t).max_val(T::from_f64(1e-12));
    let ois_fwd = (ps / pt - T::one()) / dt;
    let log_bond_ratio = (ps / pt).ln();

    if sigma.val().abs() < 1e-10 {
        let cont_fwd = log_bond_ratio / dt;
        return CirppSofrFuturesValuation {
            ois_forward_rate: ois_fwd,
            futures_rate_3m_compounded: ois_fwd,
            convexity_adjustment_3m: T::zero(),
            futures_rate_1m_arithmetic: cont_fwd,
            convexity_adjustment_1m: T::zero(),
            shift_integral_l: log_bond_ratio,
        };
    }

    let (phi_s_neg, psi_s_neg) = cir_bond_phi_psi(b, beta, sigma, s);
    let (phi_t_neg, psi_t_neg) = cir_bond_phi_psi(b, beta, sigma, t);
    let l_st = log_bond_ratio + (phi_t_neg - phi_s_neg) + (psi_t_neg - psi_s_neg) * x0;

    let (phi_tau_pos, u_s) = cir_compounded_phi_psi(b, beta, sigma, dt);
    let s2 = sigma * sigma;
    let one_minus_eb = T::one() - (-beta * s).exp();
    let denom_s = (T::one() - u_s * s2 * one_minus_eb / (beta * 2.0)).max_val(T::from_f64(1e-12));
    let phi_s_u = -(b * 2.0 / s2) * denom_s.ln();
    let psi_s_u = u_s * (-beta * s).exp() / denom_s;

    let log_mgf_3m = l_st + phi_tau_pos + phi_s_u + psi_s_u * x0;
    let fut_3m = (log_mgf_3m.exp() - T::one()) / dt;
    let conv_3m = fut_3m - ois_fwd;

    let ey_st = (b / beta) * dt + (x0 - b / beta) * ((-beta * s).exp() - (-beta * t).exp()) / beta;
    let fut_1m = (l_st + ey_st) / dt;
    let cont_fwd = log_bond_ratio / dt;
    let conv_1m = fut_1m - cont_fwd;

    CirppSofrFuturesValuation {
        ois_forward_rate: ois_fwd,
        futures_rate_3m_compounded: fut_3m,
        convexity_adjustment_3m: conv_3m,
        futures_rate_1m_arithmetic: fut_1m,
        convexity_adjustment_1m: conv_1m,
        shift_integral_l: l_st,
    }
}



#[cfg(test)]
mod phase20_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;

    #[test]
    fn sinh_ou_rfr_caplet_hw_collapse_smile_skew_and_cirpp_sofr_futures_ad() {
        let d = Date(0);
        let curve = FlatForwardCurve::new(d, 0.03);
        // When gamma_smile == 0.0, Sinh-OU collapses identically to Hull-White 1F
        let hw_spec = SinhOuRfrCapletSpec {
            notional: 1_000_000.0,
            start_time: 0.5,
            end_time: 0.75,
            strike: 0.03,
            alpha_mean_reversion: 0.05,
            sigma_short_rate: 0.01,
            gamma_smile: 0.0,
            y_star_skew: -0.2,
            is_caplet: true,
            is_backward_compounded: true,
        };
        let v_hw_collapse = price_sinh_ou_rfr_caplet(&hw_spec, &curve);
        let hw_exact = hw1f_rfr_caplet_price(
            &curve,
            0.05,
            0.01,
            0.5,
            0.75,
            0.03,
            1_000_000.0,
            true,
            RfrRateType::BackwardCompounded,
        );
        assert!((v_hw_collapse.npv - hw_exact).abs() < 1e-8);
        assert!((v_hw_collapse.implied_hw_short_rate_vol - 0.01).abs() < 1e-12);

        // Positive gamma_smile produces smile curvature C3 > 0 and put-call parity holds
        let smile_spec = SinhOuRfrCapletSpec {
            gamma_smile: 1.5,
            ..hw_spec.clone()
        };
        let cap = price_sinh_ou_rfr_caplet(&smile_spec, &curve);
        let flr = price_sinh_ou_rfr_caplet(
            &SinhOuRfrCapletSpec {
                is_caplet: false,
                ..smile_spec.clone()
            },
            &curve,
        );
        assert!(cap.c3_smile > 0.0);
        let p1 = curve.discount_factor(0.5);
        let p2 = curve.discount_factor(0.75);
        let parity = 1_000_000.0 * (p1 - (1.0 + 0.03 * 0.25) * p2);
        assert!((cap.npv - flr.npv - parity).abs() < 1e-8);

        // Fontana (2023) CIR++ SOFR futures: positive convexity adjustments for sigma > 0
        let cir_spec = CirppSofrFuturesSpec {
            start_time: 0.5,
            end_time: 0.75,
            x0: 0.03,
            b_drift: 0.006,
            beta_mean_reversion: 0.20,
            sigma: 0.06,
        };
        let fut = price_cirpp_sofr_futures(&cir_spec, &curve);
        assert!(fut.convexity_adjustment_3m > 0.0);
        assert!(fut.convexity_adjustment_1m > 0.0);
        assert!(fut.futures_rate_3m_compounded > fut.ois_forward_rate);

        let h = 1e-5;
        let eval_g = |g: f64| {
            let mut s = smile_spec.clone();
            s.gamma_smile = g;
            price_sinh_ou_rfr_caplet(&s, &curve).npv
        };
        let fd = (eval_g(1.5 + h) - eval_g(1.5 - h)) / (2.0 * h);
        type D = Dual<1>;
        let curve_d = FlatForwardCurve::new(d, D::from_value(0.03));
        let spec_d = SinhOuRfrCapletSpec {
            notional: D::from_value(1_000_000.0),
            start_time: 0.5,
            end_time: 0.75,
            strike: D::from_value(0.03),
            alpha_mean_reversion: D::from_value(0.05),
            sigma_short_rate: D::from_value(0.01),
            gamma_smile: D::variable(1.5, 0),
            y_star_skew: D::from_value(-0.2),
            is_caplet: true,
            is_backward_compounded: true,
        };
        let ad = price_sinh_ou_rfr_caplet(&spec_d, &curve_d).npv;
        assert!((ad.d[0] - fd).abs() < 1e-3);
    }
}
