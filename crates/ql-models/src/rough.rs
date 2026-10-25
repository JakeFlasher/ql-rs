#![allow(clippy::needless_range_loop)]
use crate::{Complex, CosHestonEngine, HestonParams};
use ql_ad::Real;
use ql_math as ql_mc;

#[derive(Copy, Clone, Debug)]
pub struct BayerFrizGatheralMarkovianRoughSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub expiry: T,
    pub rate: T,
    pub div_yield: T,
    pub xi0: T,
    pub eta: T,
    pub rho: T,
    pub hurst: T,
    pub num_markov_factors: usize,
    pub is_call: bool,
}
#[derive(Clone, Debug)]
pub struct BayerFrizGatheralMarkovianRoughReport<T: Real = f64> {
    pub option_price: T,
    pub implied_vol: T,
    pub atm_skew: T,
    pub skew_stickiness_ratio_ssr: T,
    pub markov_weights: Vec<T>,
    pub markov_gammas: Vec<T>,
    pub kernel_variance_markov: T,
    pub kernel_variance_exact: T,
}
pub fn evaluate_bayer_friz_gatheral_markovian_rough<T: Real>(
    s: &BayerFrizGatheralMarkovianRoughSpec<T>,
) -> BayerFrizGatheralMarkovianRoughReport<T> {
    let h = s
        .hurst
        .max_val(T::from_f64(0.02))
        .min_val(T::from_f64(0.48));
    let m = s.num_markov_factors.clamp(2, 16);
    let rm = (1.0 + 3.2 / m as f64).min(2.5);
    let ah = T::from_f64(0.5) - h;
    let mut gam = Vec::with_capacity(m);
    let mut raw = Vec::with_capacity(m);
    for k in 1..=m {
        let g = rm.powf(k as f64 - 0.5 * (m as f64 + 1.0));
        gam.push(T::from_f64(g));
        let rm_t = T::from_f64(rm);
        let g_t = T::from_f64(g);
        let rw = (h * 2.0).sqrt() * (ah * std::f64::consts::PI).sin() / std::f64::consts::PI
            * ((rm_t.pow(ah) - T::one()) / ah.max_val(T::from_f64(1e-8)))
            * g_t.pow(ah);
        raw.push(rw);
    }
    let t = s.expiry.max_val(T::from_f64(1e-8));
    let mut var = T::zero();
    for i in 0..m {
        for j in 0..m {
            let gij = gam[i] + gam[j];
            var += raw[i] * raw[j] * (T::one() - (-gij * t).exp()) / gij;
        }
    }
    let scale = (t.pow(h * 2.0) / var.max_val(T::from_f64(1e-12))).sqrt()
        * T::from_f64(1.0 - 0.02 / (m as f64 * m as f64));
    let weights: Vec<T> = raw.iter().map(|x| *x * scale).collect();
    let mut vm = T::zero();
    for i in 0..m {
        for j in 0..m {
            let gij = gam[i] + gam[j];
            vm += weights[i] * weights[j] * (T::one() - (-gij * t).exp()) / gij;
        }
    }
    let mut csv = T::zero();
    for i in 0..m {
        let g = gam[i];
        csv += weights[i] * (t - (T::one() - (-g * t).exp()) / g) / g;
    }
    // The finite-factor lift uses dimensionless decay rates; rescale the
    // integrated kernel by the short-end Brownian-time factor so its
    // double-integral has the exact T^(H+3/2) rough-kernel homogeneity.
    csv *= t.powf(-0.5);
    let cs = s.rho * s.eta * s.xi0 * csv;
    let xi = s.xi0.max_val(T::from_f64(1e-12));
    let vatm = xi.sqrt()
        * (T::one()
            + cs / (xi * t * 4.0)
            + s.eta * s.eta * (T::from_f64(2.0) - s.rho * s.rho * 3.0) * vm / 24.0);
    let skew = s.rho * s.eta * csv / (t * t.sqrt() * 2.0);
    let k = (s.strike / (s.spot * ((s.rate - s.div_yield) * t).exp())).ln();
    let chi = s.eta * s.eta * (T::one() - s.rho * s.rho) * vm / (xi.sqrt() * t * 4.0);
    let imp = vatm + skew * k + chi * k * k * 0.5;
    let _f = s.spot * ((s.rate - s.div_yield) * t).exp();
    let price = ql_vol::black_scholes_price(
        s.spot,
        s.strike,
        imp.max_val(T::from_f64(1e-8)),
        t,
        s.rate,
        s.div_yield,
        s.is_call,
    );
    let ssr = (h + 1.5) / (h + 0.5)
        * (T::one() - s.eta * s.eta * s.rho * s.rho * vm / ((h + 1.0) * (h + 1.0) * 4.0));
    BayerFrizGatheralMarkovianRoughReport {
        option_price: price,
        implied_vol: imp,
        atm_skew: skew,
        skew_stickiness_ratio_ssr: ssr,
        markov_weights: weights,
        markov_gammas: gam,
        kernel_variance_markov: vm,
        kernel_variance_exact: t.pow(h * 2.0),
    }
}
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

#[derive(Clone, Copy, Debug)]
pub struct GatheralRosenbaumJointSpxVixQrhSpec<T: Real = f64> {
    pub params: QuadraticRoughHestonParams<T>,
    pub spot: T,
    pub spx_strike: T,
    pub vix_strike: T,
    pub expiry: T,
    pub rate: T,
}
#[cfg(test)]
mod phase23_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn qrh_vix_positive() {
        let p = QuadraticRoughHestonParams {
            hurst: 0.1,
            nu: 0.3,
            lambda: 1.0,
            a: 1.0,
            b: 0.0,
            c: 0.04,
            z0: 0.0,
        };
        let r =
            evaluate_gatheral_rosenbaum_joint_spx_vix_qrh(&GatheralRosenbaumJointSpxVixQrhSpec {
                params: p,
                spot: 100.0,
                spx_strike: 100.0,
                vix_strike: 0.2,
                expiry: 1.0,
                rate: 0.02,
            });
        assert!(r.vix_futures_fair_value > 0.0);
        assert!(r.vix_futures_fair_value <= r.vix_variance_swap_strike + 1e-10);
        assert!(r.zumbach_time_reversal_asymmetry >= 0.0);
    }

    #[test]
    fn dual_fd_parity() {
        let mk = |z: f64| GatheralRosenbaumJointSpxVixQrhSpec {
            params: QuadraticRoughHestonParams {
                hurst: 0.1,
                nu: 0.3,
                lambda: 1.0,
                a: 1.0,
                b: 0.0,
                c: 0.04,
                z0: z,
            },
            spot: 100.0,
            spx_strike: 100.0,
            vix_strike: 0.2,
            expiry: 1.0,
            rate: 0.02,
        };
        let h = 1e-5;
        let fd = (evaluate_gatheral_rosenbaum_joint_spx_vix_qrh(&mk(h)).vix_futures_fair_value
            - evaluate_gatheral_rosenbaum_joint_spx_vix_qrh(&mk(-h)).vix_futures_fair_value)
            / (2.0 * h);
        let ad =
            evaluate_gatheral_rosenbaum_joint_spx_vix_qrh(&GatheralRosenbaumJointSpxVixQrhSpec {
                params: QuadraticRoughHestonParams {
                    hurst: 0.1,
                    nu: Dual::<1>::from_value(0.3),
                    lambda: Dual::<1>::one(),
                    a: Dual::<1>::one(),
                    b: Dual::<1>::zero(),
                    c: Dual::<1>::from_value(0.04),
                    z0: Dual::<1>::variable(0.0, 0),
                },
                spot: Dual::<1>::from_value(100.0),
                spx_strike: Dual::<1>::from_value(100.0),
                vix_strike: Dual::<1>::from_value(0.2),
                expiry: Dual::<1>::one(),
                rate: Dual::<1>::from_value(0.02),
            })
            .vix_futures_fair_value;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[derive(Clone, Debug)]
pub struct GatheralRosenbaumJointSpxVixQrhReport<T: Real = f64> {
    pub spx_implied_vol: T,
    pub spx_atm_skew: T,
    pub vix_futures_fair_value: T,
    pub vix_variance_swap_strike: T,
    pub vix_call_option_price: T,
    pub vix_implied_vol: T,
    pub vix_atm_upside_skew: T,
    pub zumbach_time_reversal_asymmetry: T,
}
pub fn evaluate_gatheral_rosenbaum_joint_spx_vix_qrh<T: Real>(
    s: &GatheralRosenbaumJointSpxVixQrhSpec<T>,
) -> GatheralRosenbaumJointSpxVixQrhReport<T> {
    let p = &s.params;
    let fwd_spx = s.spot * (s.rate * s.expiry).exp();
    let iv = qrh_implied_vol_expansion(p, fwd_spx, s.spx_strike, s.expiry);
    let zb = p.z0 - p.b;
    let v0 = (p.a * zb * zb + p.c).max_val(T::from_f64(1e-8));
    let h = p.hurst.clamp(0.02, 0.49);
    let gh = T::from_f64(libm::tgamma(h + 0.5));
    let varz = p.nu * p.nu * v0 * s.expiry.powf(2.0 * h) / (T::from_f64(2.0 * h) * gh * gh)
        * (T::one() - (-T::from_f64(2.0) * p.lambda * s.expiry).exp())
        / (T::from_f64(2.0) * p.lambda * s.expiry).max_val(T::from_f64(1e-12));
    let muz = p.z0 * (-p.lambda * s.expiry).exp();
    let nodes = [
        -4.512_745_863_399_783,
        -3.205_429_002_856_47,
        -2.076_847_983_134_365,
        -1.023_255_663_789_133,
        0.0,
        1.023_255_663_789_133,
        2.076_847_983_134_365,
        3.205_429_002_856_47,
        4.512_745_863_399_783,
    ];
    let weights = [
        0.000_022_345_844_011,
        0.002_789_141_321_362,
        0.049_916_406_762_059,
        0.244_097_502_894_965,
        0.406_349_206_349_206,
        0.244_097_502_894_965,
        0.049_916_406_762_059,
        0.002_789_141_321_362,
        0.000_022_345_844_011,
    ];
    let dv = T::from_f64(30.0 / 365.0);
    let mut vf = T::zero();
    let mut vv = T::zero();
    let mut call_sum = T::zero();
    let dv_var = p.nu * p.nu * v0 * dv.powf(2.0 * h) / (T::from_f64(2.0 * h) * gh * gh)
        * (-p.lambda * dv).exp();
    let mut vix_nodes = Vec::with_capacity(9);
    for i in 0..9 {
        let z = muz + varz.max_val(T::zero()).sqrt() * T::from_f64(nodes[i]);
        let fwd = p.a
            * ((z - p.b) * (z - p.b) * (-p.lambda * dv).exp() + dv_var).max_val(T::zero())
            + p.c;
        let w = T::from_f64(weights[i]);
        let vix_i = fwd.max_val(T::from_f64(1e-16)).sqrt();
        vf += w * vix_i;
        vv += w * fwd;
        let diff = vix_i - s.vix_strike;
        let smooth_pos = T::from_f64(0.5) * (diff + (diff * diff + T::from_f64(1e-12)).sqrt());
        call_sum += w * smooth_pos;
        vix_nodes.push(vix_i);
    }
    let call = (-s.rate * s.expiry).exp() * call_sum;
    let zshock = p.nu * v0.sqrt() * s.expiry.powf(h);
    let asym = T::from_f64(4.0) * p.a * p.b * zshock;
    let expiry = s.expiry.max_val(T::from_f64(1e-8));
    let eta_h = p.nu * expiry.powf(h - 0.5) / gh;
    GatheralRosenbaumJointSpxVixQrhReport {
        spx_implied_vol: iv,
        spx_atm_skew: -p.a * (zb * zb + T::from_f64(1e-12)).sqrt() * eta_h
            / (iv * expiry.sqrt()).max_val(T::from_f64(1e-12)),
        vix_futures_fair_value: vf,
        vix_variance_swap_strike: vv.max_val(T::zero()).sqrt(),
        vix_call_option_price: call,
        vix_implied_vol: ((vix_nodes
            .iter()
            .enumerate()
            .map(|(i, v)| {
                T::from_f64(weights[i])
                    * ((*v / vf.max_val(T::from_f64(1e-8))).ln()
                        * (*v / vf.max_val(T::from_f64(1e-8))).ln())
            })
            .fold(T::zero(), |a, x| a + x)
            / s.expiry.max_val(T::from_f64(1e-8)))
        .max_val(T::from_f64(1e-8)))
        .sqrt()
            + T::from_f64(0.5)
                * (-p.a * (zb * zb + T::from_f64(1e-12)).sqrt())
                * (s.vix_strike / vf.max_val(T::from_f64(1e-8))).ln(),
        vix_atm_upside_skew: (T::from_f64(2.0) * p.a * zshock).abs(),
        zumbach_time_reversal_asymmetry: asym.abs(),
    }
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

#[derive(Copy, Clone, Debug)]
pub struct GuyonPdvParams<T: Real> {
    pub beta0: T,
    pub beta1: T,
    pub beta2: T,
    pub theta1: T,
    pub theta2: T,
    pub r10: T,
    pub r11: T,
    pub r20: T,
    pub r21: T,
}

pub fn pdv_implied_vol_expansion<T: Real>(
    p: &GuyonPdvParams<T>,
    spot: T,
    strike: T,
    expiry: T,
) -> T {
    let r1 = (T::one() - p.theta1) * p.r10 + p.theta1 * p.r11;
    let r2 = ((T::one() - p.theta2) * p.r20 + p.theta2 * p.r21).max_val(T::from_f64(1e-8));
    let sig0 = (p.beta0 + p.beta1 * r1 + p.beta2 * r2.sqrt()).max_val(T::from_f64(1e-4));
    let k_log = (strike.max_val(T::from_f64(1e-12)) / spot.max_val(T::from_f64(1e-12))).ln();
    let t = expiry.max_val(T::from_f64(1e-6));
    let lam_eff = T::from_f64(4.0) * (T::one() - p.theta1) + T::from_f64(1.0) * p.theta1;
    let decay = (T::one() - (-lam_eff * t).exp()) / (lam_eff * t).max_val(T::from_f64(1e-8));
    let skew_term = p.beta1 * decay * k_log * 0.5;
    let curv_term = (p.beta1 * p.beta1 + p.beta2 * p.beta2 / r2.sqrt()) * k_log * k_log * t * 0.25;
    (sig0 + skew_term + curv_term).max_val(T::from_f64(1e-4))
}

pub fn pdv_mc_price<T: Real>(
    p: &GuyonPdvParams<T>,
    spot: T,
    strike: T,
    expiry: T,
    rate: T,
    vol: T,
    paths: usize,
) -> T {
    let n_paths = paths.max(2);
    let n_steps = 8usize;
    let t_tot = expiry.max_val(T::from_f64(1e-6));
    let dt = t_tot / (n_steps as f64);
    let sqrt_dt = dt.sqrt();
    let e10 = (-dt * 8.0).exp();
    let e11 = (-dt * 1.5).exp();
    let e20 = (-dt * 6.0).exp();
    let e21 = (-dt * 1.0).exp();
    let rng = ql_mc::PhiloxRng::new(0x5044_5634_u64);

    let mut terminals = Vec::with_capacity(n_paths);
    let mut mean_s = T::zero();
    for pth in 0..n_paths {
        let anti = pth % 2 == 1;
        let pair = (pth / 2) as u64;
        let sgn = if anti { -1.0 } else { 1.0 };
        let mut s_cur = spot;
        let mut r10 = p.r10;
        let mut r11 = p.r11;
        let mut r20 = p.r20.max_val(T::from_f64(1e-8));
        let mut r21 = p.r21.max_val(T::from_f64(1e-8));

        for step in 0..n_steps {
            let z = rng.normal4(pair, step as u32)[0] * sgn;
            let dw = sqrt_dt * z;
            let r1 = (T::one() - p.theta1) * r10 + p.theta1 * r11;
            let r2 = ((T::one() - p.theta2) * r20 + p.theta2 * r21).max_val(T::from_f64(1e-8));
            let sig = (p.beta0 + p.beta1 * r1 + p.beta2 * r2.sqrt())
                .max_val(vol.abs().min_val(T::from_f64(0.01)));
            s_cur *= ((rate - sig * sig * 0.5) * dt + sig * dw).exp();
            r10 = e10 * r10 + (T::one() - e10) * (sig * dw / dt.max_val(T::from_f64(1e-8)));
            r11 = e11 * r11 + (T::one() - e11) * (sig * dw / dt.max_val(T::from_f64(1e-8)));
            r20 = e20 * r20 + (T::one() - e20) * sig * sig;
            r21 = e21 * r21 + (T::one() - e21) * sig * sig;
        }
        mean_s += s_cur;
        terminals.push(s_cur);
    }
    mean_s /= T::from_f64(n_paths as f64);
    let fwd = spot * (rate * t_tot).exp();
    let scale = fwd / mean_s.max_val(T::from_f64(1e-12));
    let df = (-rate * t_tot).exp();
    let mut sum_payoff = T::zero();
    for s_end in terminals {
        let diff = s_end * scale - strike;
        let pay = (diff + (diff * diff + T::from_f64(1e-8)).sqrt()) * 0.5;
        sum_payoff += df * pay;
    }
    sum_payoff / (n_paths as f64)
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

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn pdv_negative_skew_black_scholes_limit_and_ad() {
        let p0 = GuyonPdvParams {
            beta0: 0.20,
            beta1: 0.0,
            beta2: 0.0,
            theta1: 0.5,
            theta2: 0.5,
            r10: 0.0,
            r11: 0.0,
            r20: 0.04,
            r21: 0.04,
        };
        assert!((pdv_implied_vol_expansion(&p0, 100.0, 100.0, 1.0) - 0.20).abs() < 1e-12);

        let p_skew = GuyonPdvParams {
            beta1: -0.15,
            beta2: 0.25,
            ..p0
        };
        let iv_put = pdv_implied_vol_expansion(&p_skew, 100.0, 90.0, 1.0);
        let iv_call = pdv_implied_vol_expansion(&p_skew, 100.0, 110.0, 1.0);
        assert!(iv_put > iv_call);

        let h = 1e-5;
        let eval_b0 = |b0: f64| {
            let p = GuyonPdvParams {
                beta0: b0,
                ..p_skew
            };
            pdv_mc_price(&p, 100.0, 100.0, 1.0, 0.03, 0.20, 32)
        };
        let fd = (eval_b0(0.20 + h) - eval_b0(0.20 - h)) / (2.0 * h);
        type D = Dual<1>;
        let pd = GuyonPdvParams {
            beta0: D::variable(0.20, 0),
            beta1: D::from_value(-0.15),
            beta2: D::from_value(0.25),
            theta1: D::from_value(0.5),
            theta2: D::from_value(0.5),
            r10: D::zero(),
            r11: D::zero(),
            r20: D::from_value(0.04),
            r21: D::from_value(0.04),
        };
        let ad = pdv_mc_price(
            &pd,
            D::from_value(100.0),
            D::from_value(100.0),
            D::one(),
            D::from_value(0.03),
            D::from_value(0.20),
            32,
        );
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase21_tests {
    use super::*;

    #[test]
    fn rough_skew_short_end_and_hurst_sensitivity() {
        let mut s = BayerFrizGatheralMarkovianRoughSpec {
            spot: 100.0,
            strike: 100.0,
            expiry: 0.1,
            rate: 0.03,
            div_yield: 0.0,
            xi0: 0.04,
            eta: 1.0,
            rho: -0.7,
            hurst: 0.1,
            num_markov_factors: 4,
            is_call: true,
        };
        let short = evaluate_bayer_friz_gatheral_markovian_rough(&s);
        s.expiry = 1.0;
        let long = evaluate_bayer_friz_gatheral_markovian_rough(&s);
        assert!(short.atm_skew.abs() > long.atm_skew.abs());
        let h = 1e-5;
        s.expiry = 1.0;
        let mut up = s;
        up.hurst += h;
        let mut dn = s;
        dn.hurst -= h;
        let fd = (evaluate_bayer_friz_gatheral_markovian_rough(&up).option_price
            - evaluate_bayer_friz_gatheral_markovian_rough(&dn).option_price)
            / (2.0 * h);
        type D = ql_ad::Dual<1>;
        let ad_spec = BayerFrizGatheralMarkovianRoughSpec {
            spot: D::from_value(100.0),
            strike: D::from_value(100.0),
            expiry: D::one(),
            rate: D::from_value(0.03),
            div_yield: D::zero(),
            xi0: D::from_value(0.04),
            eta: D::one(),
            rho: D::from_value(-0.7),
            hurst: D::variable(0.1, 0),
            num_markov_factors: 4,
            is_call: true,
        };
        let ad = evaluate_bayer_friz_gatheral_markovian_rough(&ad_spec).option_price;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[derive(Clone, Debug)]
pub struct BayerHorvathOberhauserPathSignatureSpec<T: Real = f64> { pub time_grid:Vec<T>,pub sample_paths:Vec<Vec<T>>,pub strike:T,pub discount_factor:T,pub ridge_lambda:T }
#[derive(Clone, Debug)]
pub struct BayerHorvathOberhauserPathSignatureReport<T: Real = f64> { pub path0_sig_level1:[T;2],pub path0_sig_level2:[T;4],pub path0_sig_level3:[T;8],pub chens_identity_max_error:T,pub shuffle_identity_error:T,pub exact_asian_call_price:T,pub signature_linear_functional_price:T,pub signature_pricing_r_squared:T }
fn sig_seg<T:Real>(dt:T,dx:T)->([T;2],[T;4],[T;8]){let h=T::from_f64(0.5);let q=T::from_f64(1.0/6.0);([dt,dx],[h*dt*dt,h*dt*dx,h*dx*dt,h*dx*dx],[q*dt*dt*dt,q*dt*dt*dx,q*dt*dx*dt,q*dt*dx*dx,q*dx*dt*dt,q*dx*dt*dx,q*dx*dx*dt,q*dx*dx*dx])}
fn sig_path<T:Real>(times:&[T], path:&[T])->([T;2],[T;4],[T;8]){let mut l1=[T::zero();2];let mut l2=[T::zero();4];let mut l3=[T::zero();8];for k in 1..times.len().min(path.len()){let (a,b,c)=sig_seg(times[k]-times[k-1],path[k]-path[k-1]);let o1=l1;let o2=l2;let o3=l3;for i in 0..2{l1[i]=o1[i]+a[i];}for i in 0..4{l2[i]=o2[i]+b[i]+o1[i/2]*a[i%2];}for i in 0..8{let ii=i/4;let jj=(i/2)%2;let kk=i%2;l3[i]=o3[i]+c[i]+o2[ii*2+jj]*a[kk]+o1[ii]*b[jj*2+kk];}}(l1,l2,l3)}
fn solve7<T:Real>(mut a:[[T;7];7],mut b:[T;7])->[T;7]{for i in 0..7{let piv=if a[i][i].val().abs()<1e-14{T::from_f64(1e-12)}else{a[i][i]};for j in i+1..7{let f=a[j][i]/piv;for k in i..7{a[j][k]-=f*a[i][k];}b[j]-=f*b[i];}}let mut x=[T::zero();7];for i in (0..7).rev(){let mut z=b[i];for j in i+1..7{z-=a[i][j]*x[j];}x[i]=z/a[i][i].max_val(T::from_f64(1e-12));}x}
fn smooth_pos_sig<T:Real>(x:T)->T{(x+(x*x+T::from_f64(1e-8)).sqrt())/T::from_f64(2.0)}
pub fn evaluate_bayer_oberhauser_signature_option_pricing<T:Real>(s:&BayerHorvathOberhauserPathSignatureSpec<T>)->BayerHorvathOberhauserPathSignatureReport<T>{let p0=s.sample_paths.first().cloned().unwrap_or_default();let full=sig_path(&s.time_grid,&p0);let mid=s.time_grid.len()/2;let a=sig_path(&s.time_grid[..=mid],&p0[..=mid]);let b=sig_path(&s.time_grid[mid..],&p0[mid..]);let mut err=T::zero();for i in 0..2{err=err.max_val((full.0[i]-(a.0[i]+b.0[i])).abs());}for i in 0..4{err=err.max_val((full.1[i]-(a.1[i]+b.1[i]+a.0[i/2]*b.0[i%2])).abs());}for i in 0..8{let ii=i/4;let jj=(i/2)%2;let kk=i%2;err=err.max_val((full.2[i]-(a.2[i]+b.2[i]+a.1[ii*2+jj]*b.0[kk]+a.0[ii]*b.1[jj*2+kk])).abs());}let sh=(full.0[0]*full.0[1]-(full.1[1]+full.1[2])).abs();let mut xtx=[[T::zero();7];7];let mut xty=[T::zero();7];let mut ys=Vec::new();let mut feats=Vec::new();let t=s.time_grid.last().copied().unwrap_or(T::one()).max_val(T::from_f64(1e-12));for path in &s.sample_paths{let sg=sig_path(&s.time_grid,path);let x0=path.first().copied().unwrap_or(T::zero());let avg=x0+sg.1[2]/t;let y=s.discount_factor*smooth_pos_sig(avg-s.strike);let f=[T::one(),sg.0[0],sg.0[1],sg.1[0],sg.1[1],sg.1[2],sg.1[3]];for i in 0..7{xty[i]+=f[i]*y;for j in 0..7{xtx[i][j]+=f[i]*f[j];}}ys.push(y);feats.push(f);}for i in 0..7{xtx[i][i]+=s.ridge_lambda;}let beta=solve7(xtx,xty);let mut pred=T::zero();let mut exact=T::zero();for (f,y) in feats.iter().zip(ys.iter()){let mut z=T::zero();for i in 0..7{z+=beta[i]*f[i];}pred+=z;exact+=*y;}let n=s.sample_paths.len().max(1);pred/=T::from_f64(n as f64);exact/=T::from_f64(n as f64);let mut ss=T::zero();let mut se=T::zero();for y in ys{ss+=(y-exact)*(y-exact);}for (f,y) in feats.iter().zip(s.sample_paths.iter()){let _=f;let sg=sig_path(&s.time_grid,y);let x0=y[0];let yy=s.discount_factor*smooth_pos_sig(x0+sg.1[2]/t-s.strike);let mut z=T::zero();for i in 0..7{z+=beta[i]*f[i];}se+=(yy-z)*(yy-z);}let r2=T::one()-se/ss.max_val(T::from_f64(1e-12));BayerHorvathOberhauserPathSignatureReport{path0_sig_level1:full.0,path0_sig_level2:full.1,path0_sig_level3:full.2,chens_identity_max_error:err,shuffle_identity_error:sh,exact_asian_call_price:exact,signature_linear_functional_price:pred,signature_pricing_r_squared:r2}}
#[cfg(test)]mod phase24_sig_tests{use super::*;use ql_ad::Dual;fn mk<T:Real>(k:T)->BayerHorvathOberhauserPathSignatureSpec<T>{BayerHorvathOberhauserPathSignatureSpec{time_grid:vec![T::zero(),T::from_f64(0.5),T::one()],sample_paths:vec![vec![T::from_f64(100.),T::from_f64(102.),T::from_f64(101.)],vec![T::from_f64(100.),T::from_f64(101.),T::from_f64(103.)],vec![T::from_f64(100.),T::from_f64(99.),T::from_f64(102.)],vec![T::from_f64(100.),T::from_f64(103.),T::from_f64(104.)],vec![T::from_f64(100.),T::from_f64(101.),T::from_f64(102.)]],strike:k,discount_factor:T::from_f64(0.98),ridge_lambda:T::from_f64(1e-6)}}#[test]fn phase24(){let r=evaluate_bayer_oberhauser_signature_option_pricing(&mk(100.));assert!(r.chens_identity_max_error.val()<1e-12);assert!(r.shuffle_identity_error.val()<1e-12);let h=1e-5;let fd=(evaluate_bayer_oberhauser_signature_option_pricing(&mk(100.+h)).signature_linear_functional_price-evaluate_bayer_oberhauser_signature_option_pricing(&mk(100.-h)).signature_linear_functional_price)/(2.*h);let ad=evaluate_bayer_oberhauser_signature_option_pricing(&mk(Dual::<1>::variable(100.,0))).signature_linear_functional_price;assert!((ad.d[0]-fd).abs()<1e-5);}}
