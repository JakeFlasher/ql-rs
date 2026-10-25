use ql_ad::Real;

#[derive(Copy, Clone, Debug)]
pub struct Bergomi2FParams<T: Real = f64> {
    pub xi0: T,
    pub omega: T,
    pub theta: T,
    pub k1: T,
    pub k2: T,
    pub rho12: T,
    pub rho_s1: T,
    pub rho_s2: T,
}
fn bg_j<T: Real>(a: T, t: T) -> T {
    let at = a * t;
    if at.val().abs() < 1e-5 {
        t - at * t * 0.5 + at * at * t / 6.0
    } else {
        (T::one() - (-at).exp()) / a
    }
}
fn bg_i2<T: Real>(a: T, t: T) -> T {
    let at = a * t;
    if at.val().abs() < 1e-5 {
        t * t * 0.5 - at * t * t / 6.0 + at * at * t * t / 24.0
    } else {
        (a * t - (T::one() - (-at).exp())) / (a * a)
    }
}
fn bg_k<T: Real>(a: T, b: T, t: T) -> T {
    let nodes = [
        -0.9602898565,
        -0.7966664774,
        -0.5255324099,
        -0.1834346425,
        0.1834346425,
        0.5255324099,
        0.7966664774,
        0.9602898565,
    ];
    let weights = [
        0.1012285363,
        0.2223810345,
        0.3137066459,
        0.3626837834,
        0.3626837834,
        0.3137066459,
        0.2223810345,
        0.1012285363,
    ];
    let mut z = T::zero();
    for (&x, &w) in nodes.iter().zip(weights.iter()) {
        let s = t * (x + 1.0) * 0.5;
        z += bg_j(a, s) * bg_j(b, s) * (t * 0.5 * w);
    }
    z
}
fn bg_cross<T: Real>(a: T, b: T, t: T) -> T {
    let nodes = [
        -0.9602898565,
        -0.7966664774,
        -0.5255324099,
        -0.1834346425,
        0.1834346425,
        0.5255324099,
        0.7966664774,
        0.9602898565,
    ];
    let weights = [
        0.1012285363,
        0.2223810345,
        0.3137066459,
        0.3626837834,
        0.3626837834,
        0.3137066459,
        0.2223810345,
        0.2223810345 * 0.0 + 0.1012285363,
    ];
    let mut z = T::zero();
    for (&x, &w) in nodes.iter().zip(weights.iter()) {
        let s = t * (x + 1.0) * 0.5;
        z += bg_j(b, t - s) * bg_j(a, s) * (t * 0.5 * w);
    }
    z
}

pub fn bergomi_2f_smile_moments_bg<T: Real>(p: &Bergomi2FParams<T>, expiry: T) -> (T, T, T) {
    let n = ((T::one() - p.theta) * (T::one() - p.theta)
        + p.theta * p.theta
        + p.rho12 * p.theta * (T::one() - p.theta) * 2.0)
        .max_val(T::from_f64(1e-12))
        .sqrt();
    let w1 = (T::one() - p.theta) / n;
    let w2 = p.theta / n;
    let w0 = p.xi0 * expiry;
    let s0 = p.xi0.max_val(T::from_f64(1e-14)).sqrt();
    if w0.val() <= 1e-14 {
        return (s0, T::zero(), T::zero());
    }
    let cxi = p.omega
        * p.xi0.powf(1.5)
        * (w1 * p.rho_s1 * bg_i2(p.k1, expiry) + w2 * p.rho_s2 * bg_i2(p.k2, expiry));
    let cxx = p.omega
        * p.omega
        * p.xi0
        * p.xi0
        * (w1 * w1 * bg_k(p.k1, p.k1, expiry)
            + w2 * w2 * bg_k(p.k2, p.k2, expiry)
            + w1 * w2 * p.rho12 * bg_k(p.k1, p.k2, expiry) * 2.0);
    let cmu = p.omega
        * p.omega
        * p.xi0
        * p.xi0
        * (w1
            * w1
            * p.rho_s1
            * (bg_cross(p.k1, p.k1, expiry) + p.rho_s1 * bg_k(p.k1, p.k1, expiry) * 0.5)
            + w2 * w2
                * p.rho_s2
                * (bg_cross(p.k2, p.k2, expiry) + p.rho_s2 * bg_k(p.k2, p.k2, expiry) * 0.5)
            + w1 * w2
                * (p.rho_s1 * p.rho12 * bg_cross(p.k1, p.k2, expiry)
                    + p.rho_s2 * p.rho12 * bg_cross(p.k2, p.k1, expiry)
                    + (p.rho_s1 * p.rho_s2) * bg_k(p.k1, p.k2, expiry)));
    let atm = s0
        * (T::one()
            + cxi / (w0 * 4.0)
            + (cxi * cxi * 12.0 - w0 * (w0 + 4.0) * cxx + w0 * (w0 - 4.0) * cmu * 4.0)
                / (w0 * w0 * w0 * 32.0));
    let skew =
        s0 * (cxi / (w0 * w0 * 2.0) + (w0 * cmu * 4.0 - cxi * cxi * 3.0) / (w0 * w0 * w0 * 8.0));
    let curv = s0 * (w0 * cmu * 4.0 + w0 * cxx - cxi * cxi * 6.0) / (w0 * w0 * w0 * w0 * 4.0);
    (atm, skew, curv)
}
pub fn bergomi_2f_implied_vol_bg<T: Real>(
    p: &Bergomi2FParams<T>,
    forward: T,
    strike: T,
    expiry: T,
) -> T {
    let (a, s, c) = bergomi_2f_smile_moments_bg(p, expiry);
    let k = (strike / forward).ln();
    (a + s * k + c * k * k * 0.5).max_val(T::from_f64(1e-8))
}
#[allow(clippy::too_many_arguments)]
pub fn bergomi_2f_mc_price<T: Real>(
    p: &Bergomi2FParams<T>,
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
    let dt = expiry / n_steps.max(1) as f64;
    let mut total = T::zero();
    for path in 0..n_paths.max(1) {
        let mut x = spot;
        let mut xi = p.xi0;
        let zsign = if path % 2 == 0 { 1.0 } else { -1.0 };
        for j in 0..n_steps.max(1) {
            let phase = (path as f64 * 12.9898 + j as f64 * 78.233 + seed as f64 * 0.001).sin();
            let z = zsign * phase;
            let z2 = zsign * (phase * 1.37 + 0.3).sin();
            let v = xi.max_val(T::zero()).sqrt();
            x *= ((rate - div - v * v * 0.5) * dt + v * dt.sqrt() * z).exp();
            let f1 = (-p.k1 * dt).exp();
            let f2 = (-p.k2 * dt).exp();
            xi = p.xi0
                + p.omega
                    * (wgt(p.theta, p.rho12) * f1 * z * dt.sqrt()
                        + wgt2(p.theta, p.rho12) * f2 * z2 * dt.sqrt());
        }
        let payoff = if is_call {
            (x - strike).max_val(T::zero())
        } else {
            (strike - x).max_val(T::zero())
        };
        total += payoff * (-rate * expiry).exp();
    }
    total / n_paths.max(1) as f64
}
fn wgt<T: Real>(theta: T, rho: T) -> T {
    (T::one() - theta).max_val(T::zero())
        / ((T::one() - theta) * (T::one() - theta)
            + theta * theta
            + rho * theta * (T::one() - theta) * 2.0)
            .max_val(T::from_f64(1e-12))
            .sqrt()
}
fn wgt2<T: Real>(theta: T, rho: T) -> T {
    theta.max_val(T::zero())
        / ((T::one() - theta) * (T::one() - theta)
            + theta * theta
            + rho * theta * (T::one() - theta) * 2.0)
            .max_val(T::from_f64(1e-12))
            .sqrt()
}

#[cfg(test)]
mod phase15_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn bergomi_limits() {
        let p = Bergomi2FParams {
            xi0: 0.04,
            omega: 0.0,
            theta: 0.5,
            k1: 1.0,
            k2: 2.0,
            rho12: 0.2,
            rho_s1: -0.5,
            rho_s2: -0.3,
        };
        let (a, s, c) = bergomi_2f_smile_moments_bg(&p, 1.0);
        assert!((a - 0.2).abs() < 1e-12 && s.abs() < 1e-12 && c.abs() < 1e-12);
        let pd = Bergomi2FParams::<Dual<1>> {
            xi0: Dual::from_value(0.04),
            omega: Dual::variable(0.2, 0),
            theta: Dual::from_value(0.5),
            k1: Dual::from_value(1.0),
            k2: Dual::from_value(2.0),
            rho12: Dual::from_value(0.2),
            rho_s1: Dual::from_value(-0.5),
            rho_s2: Dual::from_value(-0.3),
        };
        assert!(
            bergomi_2f_implied_vol_bg(
                &pd,
                Dual::from_value(100.0),
                Dual::from_value(100.0),
                Dual::from_value(1.0)
            )
            .d[0]
                .is_finite()
        );
    }
}
