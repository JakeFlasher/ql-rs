use nalgebra::{DMatrix, DVector};
use ql_ad::{Dual, Dual2, IftSystemSensitivity, Real};
use ql_math::{MathError, levenberg_marquardt};
pub use volsurf::smile::{SabrSmile, SviSmile};
pub use volsurf::surface::{EssviSurface, SsviSurface};
pub use volsurf::{
    DupireLocalVol as DupireTotalVarianceFormula, SmileSection, Strike, Tenor, Vol, VolSurface,
};
pub type DupireLocalVol = DupireTotalVarianceFormula;
#[cfg(test)]
mod phase23_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn lsv_projection_positive() {
        let s = GuyonMarkovianProjectionLsvSpec {
            forward: 100.0,
            strike: 105.0,
            expiry: 1.0,
            forward_start_t1: 0.25,
            atm_vol: 0.2,
            market_atm_skew: -0.1,
            market_atm_curvature: 0.1,
            sv_v0: 0.04,
            sv_kappa: 1.0,
            sv_eta: 0.3,
            sv_rho: -0.5,
            sv_mixing_weight_eta: 0.5,
        };
        let r = evaluate_guyon_markovian_projection_lsv_smile(&s);
        assert!(r.dupire_local_vol_at_strike > 0.0);
        assert!(r.conditional_expected_sv_at_strike > 0.0);
    }

    #[test]
    fn dual_fd_parity() {
        let mk = |k: f64| GuyonMarkovianProjectionLsvSpec {
            forward: 100.0,
            strike: k,
            expiry: 1.0,
            forward_start_t1: 0.25,
            atm_vol: 0.2,
            market_atm_skew: -0.1,
            market_atm_curvature: 0.1,
            sv_v0: 0.04,
            sv_kappa: 1.0,
            sv_eta: 0.3,
            sv_rho: -0.5,
            sv_mixing_weight_eta: 0.5,
        };
        let h = 1e-5;
        let fd = (evaluate_guyon_markovian_projection_lsv_smile(&mk(100.0 + h))
            .lsv_reconstructed_implied_vol
            - evaluate_guyon_markovian_projection_lsv_smile(&mk(100.0 - h))
                .lsv_reconstructed_implied_vol)
            / (2.0 * h);
        let x = Dual::<1>::variable(100.0, 0);
        let ad = evaluate_guyon_markovian_projection_lsv_smile(&GuyonMarkovianProjectionLsvSpec {
            forward: Dual::<1>::from_value(100.0),
            strike: x,
            expiry: Dual::<1>::one(),
            forward_start_t1: Dual::<1>::from_value(0.25),
            atm_vol: Dual::<1>::from_value(0.2),
            market_atm_skew: Dual::<1>::from_value(-0.1),
            market_atm_curvature: Dual::<1>::from_value(0.1),
            sv_v0: Dual::<1>::from_value(0.04),
            sv_kappa: Dual::<1>::one(),
            sv_eta: Dual::<1>::from_value(0.3),
            sv_rho: Dual::<1>::from_value(-0.5),
            sv_mixing_weight_eta: Dual::<1>::from_value(0.5),
        })
        .lsv_reconstructed_implied_vol;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[derive(Clone, Copy, Debug)]
pub struct GuyonMarkovianProjectionLsvSpec<T: Real = f64> {
    pub forward: T,
    pub strike: T,
    pub expiry: T,
    pub forward_start_t1: T,
    pub atm_vol: T,
    pub market_atm_skew: T,
    pub market_atm_curvature: T,
    pub sv_v0: T,
    pub sv_kappa: T,
    pub sv_eta: T,
    pub sv_rho: T,
    pub sv_mixing_weight_eta: T,
}
#[derive(Clone, Debug)]
pub struct GuyonMarkovianProjectionLsvReport<T: Real = f64> {
    pub dupire_local_vol_at_strike: T,
    pub conditional_expected_sv_at_strike: T,
    pub calibrated_lsv_leverage_at_strike: T,
    pub lsv_reconstructed_implied_vol: T,
    pub calibration_implied_vol_error: T,
    pub lsv_skew_stickiness_ratio_ssr: T,
    pub lsv_forward_start_atm_skew: T,
}
pub fn evaluate_guyon_markovian_projection_lsv_smile<T: Real>(
    s: &GuyonMarkovianProjectionLsvSpec<T>,
) -> GuyonMarkovianProjectionLsvReport<T> {
    let t = s.expiry.max_val(T::from_f64(1e-8));
    let x = (s.strike / s.forward).ln();
    let eta = s.sv_mixing_weight_eta.max_val(T::zero()).min_val(T::one());
    let k = s.sv_kappa;
    let i1 = if (k * t).abs().val() < 1e-5 {
        t - T::from_f64(0.5) * k * t * t
    } else {
        (T::one() - (-k * t).exp()) / k
    };
    let i2 = if (k * t).abs().val() < 1e-5 {
        T::from_f64(0.5) * t * t - T::from_f64(1.0 / 6.0) * k * t * t * t
    } else {
        (k * t - T::one() + (-k * t).exp()) / (k * k)
    };
    let isv2 = if (T::from_f64(2.0) * k * t).abs().val() < 1e-5 {
        T::from_f64(0.5) * t * t
    } else {
        (T::from_f64(2.0) * k * t - T::one() + (-T::from_f64(2.0) * k * t).exp())
            / (T::from_f64(4.0) * k * k)
    };
    let m = (T::one()
        + eta
            * ((s.sv_rho * s.sv_eta / s.atm_vol.max_val(T::from_f64(1e-8)) * i2 / (t * t)) * x
                + T::from_f64(0.5)
                    * s.sv_eta
                    * s.sv_eta
                    * (T::one() - s.sv_rho * s.sv_rho)
                    * isv2
                    / (t * t)
                    * x
                    * x))
        .max_val(T::from_f64(1e-6));
    let dup2 = (s.atm_vol * s.atm_vol
        + T::from_f64(4.0) * s.atm_vol * s.market_atm_skew * x
        + (T::from_f64(3.0) * s.atm_vol * s.market_atm_curvature
            + T::from_f64(6.0) * s.market_atm_skew * s.market_atm_skew)
            * x
            * x)
        .max_val(T::from_f64(1e-8));
    let dup = dup2.sqrt();
    let cond = s.sv_v0 * m;
    let lev = dup / cond.sqrt();
    let gl = [
        -0.9602898565,
        -0.7966664774,
        -0.5255324099,
        -0.1834346425,
        0.1834346425,
        0.5255324099,
        0.7966664774,
        0.9602898565,
    ];
    let gw = [
        0.1012285363,
        0.2223810345,
        0.3137066459,
        0.3626837834,
        0.3626837834,
        0.3137066459,
        0.2223810345,
        0.1012285363,
    ];
    let mut integ = T::zero();
    for j in 0..8 {
        let u = T::from_f64(0.5) * (T::one() + T::from_f64(gl[j]));
        let y = u * x;
        let d = (s.atm_vol * s.atm_vol
            + T::from_f64(4.0) * s.atm_vol * s.market_atm_skew * y
            + (T::from_f64(3.0) * s.atm_vol * s.market_atm_curvature
                + T::from_f64(6.0) * s.market_atm_skew * s.market_atm_skew)
                * y
                * y)
            .max_val(T::from_f64(1e-8))
            .sqrt();
        let my = (T::one()
            + eta
                * ((s.sv_rho * s.sv_eta / s.atm_vol.max_val(T::from_f64(1e-8)) * i2 / (t * t))
                    * y
                    + T::from_f64(0.5)
                        * s.sv_eta
                        * s.sv_eta
                        * (T::one() - s.sv_rho * s.sv_rho)
                        * isv2
                        / (t * t)
                        * y
                        * y))
            .max_val(T::from_f64(1e-6));
        let cond_y = (s.sv_v0 * my).max_val(T::from_f64(1e-12));
        let lev_y = d / cond_y.sqrt();
        let eff_vol_y = (lev_y * cond_y.sqrt()).max_val(T::from_f64(1e-12));
        integ += T::from_f64(0.5 * gw[j]) / eff_vol_y;
    }
    let rec = T::one() / integ.max_val(T::from_f64(1e-12));
    let mkt = s.atm_vol + s.market_atm_skew * x + T::from_f64(0.5) * s.market_atm_curvature * x * x;
    let rem = t - s.forward_start_t1;
    let i2r = if rem.val() > 1e-8 {
        (k * rem - T::one() + (-k * rem).exp()) / (k * k * rem * rem)
    } else {
        T::from_f64(0.5)
    };
    let i2t = i2 / (t * t);
    let ssr = (T::one() - eta) * T::from_f64(2.0)
        + eta * (T::from_f64(0.5) * t * i1 / i2.max_val(T::from_f64(1e-12)));
    GuyonMarkovianProjectionLsvReport {
        dupire_local_vol_at_strike: dup,
        conditional_expected_sv_at_strike: cond,
        calibrated_lsv_leverage_at_strike: lev,
        lsv_reconstructed_implied_vol: rec,
        calibration_implied_vol_error: (rec - mkt).abs(),
        lsv_skew_stickiness_ratio_ssr: ssr,
        lsv_forward_start_atm_skew: s.market_atm_skew
            * ((T::one() - eta) * (rem / t).max_val(T::zero()).sqrt()
                + eta * i2r / i2t.max_val(T::from_f64(1e-12))),
    }
}
#[derive(Debug, Clone)]
pub struct AdSviSmile<T: Real> {
    pub forward: T,
    pub expiry: T,
    pub a: T,
    pub b: T,
    pub rho: T,
    pub m: T,
    pub sigma: T,
}
impl<T: Real> AdSviSmile<T> {
    pub fn new(forward: T, expiry: T, a: T, b: T, rho: T, m: T, sigma: T) -> Self {
        assert!(
            forward.val() > 0.0 && expiry.val() > 0.0,
            "forward and expiry must be positive"
        );
        assert!(
            b.val() >= 0.0 && sigma.val() > 0.0 && rho.val().abs() < 1.0,
            "invalid SVI parameters"
        );
        Self {
            forward,
            expiry,
            a,
            b,
            rho,
            m,
            sigma,
        }
    }
    pub fn from_svi(s: &SviSmile) -> Self {
        Self::new(
            T::from_f64(s.forward()),
            T::from_f64(s.expiry()),
            T::from_f64(s.a()),
            T::from_f64(s.b()),
            T::from_f64(s.rho()),
            T::from_f64(s.m()),
            T::from_f64(s.sigma()),
        )
    }
    pub fn total_variance(&self, strike: T) -> T {
        let k = (strike / self.forward).ln();
        let d = k - self.m;
        self.a + self.b * (self.rho * d + (d * d + self.sigma * self.sigma).sqrt())
    }
    pub fn vol(&self, strike: T) -> T {
        (self.total_variance(strike) / self.expiry).sqrt()
    }
    pub fn eval(&self, strike: T) -> T {
        self.vol(strike)
    }
    pub fn variance(&self, strike: T) -> T {
        self.total_variance(strike)
    }
}
#[derive(Debug, Clone)]
pub struct AdSabrSmile<T: Real> {
    pub forward: T,
    pub expiry: T,
    pub alpha: T,
    pub beta: T,
    pub rho: T,
    pub nu: T,
}
impl<T: Real> AdSabrSmile<T> {
    pub fn new(forward: T, expiry: T, alpha: T, beta: T, rho: T, nu: T) -> Self {
        assert!(
            forward.val() > 0.0 && expiry.val() > 0.0 && alpha.val() > 0.0,
            "forward, expiry and alpha must be positive"
        );
        assert!(
            (0.0..=1.0).contains(&beta.val()) && nu.val() >= 0.0 && rho.val().abs() < 1.0,
            "invalid SABR parameters"
        );
        Self {
            forward,
            expiry,
            alpha,
            beta,
            rho,
            nu,
        }
    }
    pub fn from_sabr(s: &SabrSmile) -> Self {
        Self::new(
            T::from_f64(s.forward()),
            T::from_f64(s.expiry()),
            T::from_f64(s.alpha()),
            T::from_f64(s.beta()),
            T::from_f64(s.rho()),
            T::from_f64(s.nu()),
        )
    }
    pub fn vol(&self, strike: T) -> T {
        let f = self.forward;
        let k = strike;
        let omb = T::one() - self.beta;
        let ln = (f / k).ln();
        let mid = (f * k).pow(omb / T::from_f64(2.0));
        let denom = mid
            * (T::one()
                + omb * omb * ln * ln / T::from_f64(24.0)
                + omb.powf(4.0) * ln.powf(4.0) / T::from_f64(1920.0));
        let z = self.nu / self.alpha * mid * ln;
        let ratio = if z.val().abs() < 1e-4 {
            T::one() - self.rho * z / T::from_f64(2.0)
                + (T::from_f64(2.0) - T::from_f64(3.0) * self.rho * self.rho) * z * z
                    / T::from_f64(12.0)
                + self.rho * (self.rho * self.rho * 5.0 - 3.0) * z * z * z / T::from_f64(24.0)
        } else {
            let disc = (T::one() - T::from_f64(2.0) * self.rho * z + z * z).sqrt();
            z / ((disc + z - self.rho) / (T::one() - self.rho)).ln()
        };
        let corr = T::one()
            + self.expiry
                * (omb * omb * self.alpha * self.alpha / (T::from_f64(24.0) * (f * k).pow(omb))
                    + self.rho * self.beta * self.nu * self.alpha / (T::from_f64(4.0) * mid)
                    + (T::from_f64(2.0) - T::from_f64(3.0) * self.rho * self.rho)
                        * self.nu
                        * self.nu
                        / T::from_f64(24.0));
        self.alpha / denom * ratio * corr
    }
    pub fn eval(&self, strike: T) -> T {
        self.vol(strike)
    }
    pub fn variance(&self, strike: T) -> T {
        let v = self.vol(strike);
        v * v * self.expiry
    }
}
pub fn svi_quote_sensitivities<T: Real>(smile: &AdSviSmile<T>, strike: T) -> [T; 5] {
    let k = (strike / smile.forward).ln();
    let d = k - smile.m;
    let q = (d * d + smile.sigma * smile.sigma).sqrt();
    [
        T::one(),
        smile.rho * d + q,
        smile.b * d,
        -smile.b * (smile.rho + d / q),
        smile.b * smile.sigma / q,
    ]
}
pub fn sabr_quote_sensitivities<T: Real>(smile: &AdSabrSmile<T>, strike: T) -> [T; 4] {
    // A small forward-mode jet keeps the quote sensitivities in the same AD
    // domain as the supplied smile.  Converting an inner `Dual<4>` to f64
    // would otherwise discard, e.g., mixed derivatives when `T` itself is a
    // Dual or TapeVar.
    #[derive(Copy, Clone)]
    struct Jet<U: Real> {
        v: U,
        d: [U; 4],
    }
    let constant = |v: T| Jet {
        v,
        d: [T::zero(); 4],
    };
    let variable = |v: T, i: usize| {
        let mut d = [T::zero(); 4];
        d[i] = T::one();
        Jet { v, d }
    };
    let add = |a: Jet<T>, b: Jet<T>| Jet {
        v: a.v + b.v,
        d: std::array::from_fn(|i| a.d[i] + b.d[i]),
    };
    let sub = |a: Jet<T>, b: Jet<T>| Jet {
        v: a.v - b.v,
        d: std::array::from_fn(|i| a.d[i] - b.d[i]),
    };
    let mul = |a: Jet<T>, b: Jet<T>| Jet {
        v: a.v * b.v,
        d: std::array::from_fn(|i| a.d[i] * b.v + a.v * b.d[i]),
    };
    let div = |a: Jet<T>, b: Jet<T>| Jet {
        v: a.v / b.v,
        d: std::array::from_fn(|i| (a.d[i] * b.v - a.v * b.d[i]) / (b.v * b.v)),
    };
    let scale = |a: Jet<T>, c: T| Jet {
        v: a.v * c,
        d: std::array::from_fn(|i| a.d[i] * c),
    };
    let ln = |a: Jet<T>| Jet {
        v: a.v.ln(),
        d: std::array::from_fn(|i| a.d[i] / a.v),
    };
    let sqrt = |a: Jet<T>| {
        let v = a.v.sqrt();
        if v.val().abs() < 1e-15 {
            Jet {
                v,
                d: [T::zero(); 4],
            }
        } else {
            Jet {
                v,
                d: std::array::from_fn(|i| a.d[i] / (T::from_f64(2.0) * v)),
            }
        }
    };
    let f = constant(smile.forward);
    let k = constant(strike);
    let expiry = constant(smile.expiry);
    let alpha = variable(smile.alpha, 0);
    let beta = variable(smile.beta, 1);
    let rho = variable(smile.rho, 2);
    let nu = variable(smile.nu, 3);
    let omb = sub(constant(T::one()), beta);
    let ln_fk = ln(div(f, k));
    // The exponent depends on beta, so construct it as a jet rather than
    // using `powf` with a scalar exponent.
    let mid = {
        let log_fk = ln(mul(f, k));
        let exponent = scale(omb, T::from_f64(0.5));
        let v = (log_fk.v * exponent.v).exp();
        Jet {
            v,
            d: std::array::from_fn(|i| v * (log_fk.d[i] * exponent.v + log_fk.v * exponent.d[i])),
        }
    };
    let l2 = mul(ln_fk, ln_fk);
    let l4 = mul(l2, l2);
    let omb2 = mul(omb, omb);
    let omb4 = mul(omb2, omb2);
    let one = constant(T::one());
    let denom_corr = add(
        one,
        add(
            scale(mul(omb2, l2), T::from_f64(1.0 / 24.0)),
            scale(mul(omb4, l4), T::from_f64(1.0 / 1920.0)),
        ),
    );
    let denom = mul(mid, denom_corr);
    let z = div(mul(nu, mul(mid, ln_fk)), alpha);
    let ratio = if z.v.val().abs() < 1e-4 {
        add(
            one,
            add(
                scale(mul(rho, z), T::from_f64(-0.5)),
                add(
                    scale(
                        mul(
                            sub(
                                constant(T::from_f64(2.0)),
                                scale(mul(rho, rho), T::from_f64(3.0)),
                            ),
                            mul(z, z),
                        ),
                        T::from_f64(1.0 / 12.0),
                    ),
                    scale(
                        mul(
                            mul(
                                rho,
                                sub(
                                    scale(mul(rho, rho), T::from_f64(5.0)),
                                    constant(T::from_f64(3.0)),
                                ),
                            ),
                            mul(mul(z, z), z),
                        ),
                        T::from_f64(1.0 / 24.0),
                    ),
                ),
            ),
        )
    } else {
        let disc = sqrt(add(
            sub(one, scale(mul(rho, z), T::from_f64(2.0))),
            mul(z, z),
        ));
        let q = div(sub(add(disc, z), rho), sub(one, rho));
        div(z, ln(q))
    };
    let fk_pow_omb = {
        let log_fk = ln(mul(f, k));
        let exponent = omb;
        let v = (log_fk.v * exponent.v).exp();
        Jet {
            v,
            d: std::array::from_fn(|i| v * (log_fk.d[i] * exponent.v + log_fk.v * exponent.d[i])),
        }
    };
    let term1 = div(
        scale(mul(omb2, mul(alpha, alpha)), T::from_f64(1.0 / 24.0)),
        fk_pow_omb,
    );
    let term2 = div(
        scale(mul(mul(rho, beta), mul(nu, alpha)), T::from_f64(0.25)),
        mid,
    );
    let term3 = scale(
        mul(
            sub(
                constant(T::from_f64(2.0)),
                scale(mul(rho, rho), T::from_f64(3.0)),
            ),
            mul(nu, nu),
        ),
        T::from_f64(1.0 / 24.0),
    );
    let corr = add(one, mul(expiry, add(term1, add(term2, term3))));
    let vol = mul(div(alpha, denom), mul(ratio, corr));
    vol.d
}

pub fn svi_durrleman_g<T: Real>(smile: &AdSviSmile<T>, strike: T) -> T {
    let k = (strike / smile.forward).ln();
    let d = k - smile.m;
    let q = (d * d + smile.sigma * smile.sigma).sqrt();
    let w = smile.total_variance(strike);
    let wp = smile.b * (smile.rho + d / q);
    let wpp = smile.b * smile.sigma * smile.sigma / (q * q * q);
    (T::one() - k * wp / (w * 2.0)).powf(2.0) - wp * wp / 4.0 * (T::one() / w + 0.25) + wpp * 0.5
}

pub fn svi_risk_neutral_density<T: Real>(smile: &AdSviSmile<T>, strike: T) -> T {
    let k = (strike / smile.forward).ln();
    let w = smile.total_variance(strike).max_val(T::from_f64(1e-14));
    let dm = -k / w.sqrt() - w.sqrt() * 0.5;
    svi_durrleman_g(smile, strike) / (strike * (T::from_f64(2.0 * std::f64::consts::PI) * w).sqrt())
        * (-dm * dm * 0.5).exp()
}

impl<T: Real> AdSviSmile<T> {
    pub fn is_butterfly_arbitrage_free(&self, k_grid: &[f64]) -> bool {
        self.a.val() + self.b.val() * self.sigma.val() * (1.0 - self.rho.val().powi(2)).sqrt()
            >= -1e-12
            && self.b.val() * (1.0 + self.rho.val().abs()) <= 2.0 + 1e-12
            && k_grid
                .iter()
                .all(|&k| svi_durrleman_g(self, self.forward * T::from_f64(k.exp())).val() >= -1e-9)
    }
}

#[derive(Copy, Clone, Debug)]
pub struct SmileStrikeQuote {
    pub strike: f64,
    pub market_vol: f64,
    pub weight: f64,
}

#[derive(Clone, Debug)]
pub struct CalibratedSabrSmile {
    pub smile: AdSabrSmile<f64>,
    pub quotes: Vec<SmileStrikeQuote>,
    pub rmse: f64,
    pub sensitivity: IftSystemSensitivity,
}

#[derive(Clone, Debug)]
pub struct CalibratedSviSmile {
    pub smile: AdSviSmile<f64>,
    pub quotes: Vec<SmileStrikeQuote>,
    pub rmse: f64,
    pub min_durrleman_g: f64,
    pub sensitivity: IftSystemSensitivity,
}

impl CalibratedSabrSmile {
    pub fn calibrate(
        forward: f64,
        expiry: f64,
        beta: f64,
        initial: (f64, f64, f64),
        quotes: &[SmileStrikeQuote],
    ) -> Result<Self, MathError> {
        let x0 = DVector::from_vec(vec![initial.0, initial.1, initial.2]);
        let residual = |x: &DVector<f64>| {
            let mut r = DVector::zeros(quotes.len());
            let mut j = DMatrix::zeros(quotes.len(), 3);
            for (i, q) in quotes.iter().enumerate() {
                let p = AdSabrSmile::new(
                    forward,
                    expiry,
                    x[0].max(1e-8),
                    beta,
                    x[1].clamp(-0.999, 0.999),
                    x[2].max(1e-8),
                );
                let d = q.weight.max(0.0).sqrt();
                r[i] = d * (p.vol(q.strike) - q.market_vol);
                let ad = AdSabrSmile::new(
                    Dual::<3>::from_value(forward),
                    Dual::<3>::from_value(expiry),
                    Dual::<3>::variable(p.alpha, 0),
                    Dual::<3>::from_value(beta),
                    Dual::<3>::variable(p.rho, 1),
                    Dual::<3>::variable(p.nu, 2),
                );
                let y = ad.vol(Dual::<3>::from_value(q.strike));
                for k in 0..3 {
                    j[(i, k)] = d * y.d[k];
                }
            }
            (r, j)
        };
        let sol = levenberg_marquardt(x0, residual, 1e-11, 200)?;
        let smile = AdSabrSmile::new(
            forward,
            expiry,
            sol[0].max(1e-8),
            beta,
            sol[1].clamp(-0.999, 0.999),
            sol[2].max(1e-8),
        );
        let (r, _) = residual(&sol);
        let rmse = (r.norm_squared() / quotes.len().max(1) as f64).sqrt();
        let mut h = DMatrix::zeros(3, 3);
        let mut jq = DMatrix::zeros(3, quotes.len());
        for (i, q) in quotes.iter().enumerate() {
            let d = q.weight.max(0.0).sqrt();
            let vars = [
                Dual2::<3>::variable(smile.alpha, 0),
                Dual2::from_value(beta),
                Dual2::<3>::variable(smile.rho, 1),
                Dual2::<3>::variable(smile.nu, 2),
            ];
            let p = AdSabrSmile::new(
                Dual2::from_value(forward),
                Dual2::from_value(expiry),
                vars[0],
                vars[1],
                vars[2],
                vars[3],
            );
            let y = p.vol(Dual2::from_value(q.strike));
            for a in 0..3 {
                jq[(a, i)] = -d * d * y.d1[a];
                for b in 0..3 {
                    h[(a, b)] += d * d * y.d1[a] * y.d1[b] + d * r[i] * y.d2[a][b];
                }
            }
        }
        Ok(Self {
            smile,
            quotes: quotes.to_vec(),
            rmse,
            sensitivity: IftSystemSensitivity::new(h, jq),
        })
    }
    pub fn vol_and_quote_vega(&self, strike: f64) -> (f64, Vec<f64>) {
        let p = AdSabrSmile::new(
            self.smile.forward,
            self.smile.expiry,
            self.smile.alpha,
            self.smile.beta,
            self.smile.rho,
            self.smile.nu,
        );
        let d = sabr_quote_sensitivities(&p, strike);
        let g = [d[0], d[2], d[3]];
        (p.vol(strike), self.sensitivity.vjp(&g))
    }
    pub fn implied_vol_ad<T: Real>(&self, strike: T, quotes: &[T]) -> T {
        let j = self.sensitivity.jacobian_dtheta_dq();
        let mut p = [
            T::from_f64(self.smile.alpha),
            T::from_f64(self.smile.rho),
            T::from_f64(self.smile.nu),
        ];
        for a in 0..3 {
            for q in 0..self.quotes.len().min(quotes.len()) {
                p[a] += T::from_f64(j[(a, q)]) * (quotes[q] - self.quotes[q].market_vol);
            }
        }
        AdSabrSmile::new(
            T::from_f64(self.smile.forward),
            T::from_f64(self.smile.expiry),
            p[0],
            T::from_f64(self.smile.beta),
            p[1],
            p[2],
        )
        .vol(strike)
    }
    pub fn value_and_quote_vega<F>(&self, pricer: F) -> (f64, Vec<f64>)
    where
        F: Fn(&AdSabrSmile<Dual<3>>) -> Dual<3>,
    {
        let p = AdSabrSmile::new(
            Dual::from_value(self.smile.forward),
            Dual::from_value(self.smile.expiry),
            Dual::variable(self.smile.alpha, 0),
            Dual::from_value(self.smile.beta),
            Dual::variable(self.smile.rho, 1),
            Dual::variable(self.smile.nu, 2),
        );
        let out = pricer(&p);
        (out.val, self.sensitivity.vjp(&out.d))
    }
}

impl CalibratedSviSmile {
    pub fn calibrate(
        forward: f64,
        expiry: f64,
        initial: (f64, f64, f64, f64, f64),
        quotes: &[SmileStrikeQuote],
    ) -> Result<Self, MathError> {
        let x0 = DVector::from_vec(vec![initial.0, initial.1, initial.2, initial.3, initial.4]);
        let residual = |x: &DVector<f64>| {
            let smile = AdSviSmile::new(
                forward,
                expiry,
                x[0],
                x[1].max(1e-8),
                x[2].clamp(-0.999, 0.999),
                x[3],
                x[4].abs().max(1e-8),
            );
            let mut rv = Vec::with_capacity(quotes.len() + 63);
            let mut jv: Vec<Vec<f64>> = Vec::with_capacity(quotes.len() + 63);
            for q in quotes {
                let d = q.weight.max(0.0).sqrt();
                let ad = AdSviSmile::new(
                    Dual::<5>::from_value(forward),
                    Dual::<5>::from_value(expiry),
                    Dual::<5>::variable(smile.a, 0),
                    Dual::<5>::variable(smile.b, 1),
                    Dual::<5>::variable(smile.rho, 2),
                    Dual::<5>::variable(smile.m, 3),
                    Dual::<5>::variable(smile.sigma, 4),
                );
                let y = ad.vol(Dual::<5>::from_value(q.strike));
                rv.push(d * (y.val - q.market_vol));
                jv.push(y.d.iter().map(|x| d * *x).collect());
            }
            let ad = AdSviSmile::new(
                Dual::<5>::from_value(forward),
                Dual::<5>::from_value(expiry),
                Dual::<5>::variable(smile.a, 0),
                Dual::<5>::variable(smile.b, 1),
                Dual::<5>::variable(smile.rho, 2),
                Dual::<5>::variable(smile.m, 3),
                Dual::<5>::variable(smile.sigma, 4),
            );
            for ngrid in 0..21 {
                let k = -1.0 + 2.0 * ngrid as f64 / 20.0;
                let strike = Dual::<5>::from_value(forward * k.exp());
                let g = svi_durrleman_g(&ad, strike).min_val(Dual::<5>::zero());
                rv.push(10.0 * g.val);
                jv.push(g.d.iter().map(|x| 10.0 * *x).collect());
            }
            let floor = ad.a + ad.b * ad.sigma * (Dual::<5>::one() - ad.rho * ad.rho).sqrt();
            let wing = ad.b * (Dual::<5>::one() + ad.rho.abs()) - 2.0;
            for p in [
                (-floor).max_val(Dual::<5>::zero()),
                wing.max_val(Dual::<5>::zero()),
            ] {
                rv.push(10.0 * p.val);
                jv.push(p.d.iter().map(|x| 10.0 * *x).collect());
            }
            let r = DVector::from_vec(rv);
            let mut j = DMatrix::zeros(r.len(), 5);
            for (i, row) in jv.iter().enumerate() {
                for a in 0..5 {
                    j[(i, a)] = row[a];
                }
            }
            (r, j)
        };
        let sol = levenberg_marquardt(x0, residual, 1e-11, 200)?;
        let smile = AdSviSmile::new(
            forward,
            expiry,
            sol[0],
            sol[1].max(1e-8),
            sol[2].clamp(-0.999, 0.999),
            sol[3],
            sol[4].abs().max(1e-8),
        );
        let (r, _) = residual(&sol);
        let rmse = quotes
            .iter()
            .map(|q| (smile.vol(q.strike) - q.market_vol).powi(2) * q.weight.max(0.0))
            .sum::<f64>()
            / quotes.len().max(1) as f64;
        let rmse = rmse.sqrt();
        let grid = (0..21)
            .map(|i| -1.0 + 2.0 * i as f64 / 20.0)
            .collect::<Vec<_>>();
        let min_g = grid
            .iter()
            .map(|&k| svi_durrleman_g(&smile, forward * k.exp()).val())
            .fold(f64::INFINITY, f64::min);
        let mut h = DMatrix::zeros(5, 5);
        let mut jq = DMatrix::zeros(5, quotes.len());
        for (i, q) in quotes.iter().enumerate() {
            let d = q.weight.max(0.0).sqrt();
            let v = [
                Dual2::<5>::variable(smile.a, 0),
                Dual2::variable(smile.b, 1),
                Dual2::variable(smile.rho, 2),
                Dual2::variable(smile.m, 3),
                Dual2::variable(smile.sigma, 4),
            ];
            let p = AdSviSmile::new(
                Dual2::from_value(forward),
                Dual2::from_value(expiry),
                v[0],
                v[1],
                v[2],
                v[3],
                v[4],
            );
            let y = p.vol(Dual2::from_value(q.strike));
            for a in 0..5 {
                jq[(a, i)] = -d * d * y.d1[a];
                for b in 0..5 {
                    h[(a, b)] += d * d * y.d1[a] * y.d1[b] + d * r[i] * y.d2[a][b];
                }
            }
        }
        let p2 = AdSviSmile::new(
            Dual2::<5>::from_value(forward),
            Dual2::<5>::from_value(expiry),
            Dual2::<5>::variable(smile.a, 0),
            Dual2::<5>::variable(smile.b, 1),
            Dual2::<5>::variable(smile.rho, 2),
            Dual2::<5>::variable(smile.m, 3),
            Dual2::<5>::variable(smile.sigma, 4),
        );
        for ngrid in 0..21 {
            let k = -1.0 + 2.0 * ngrid as f64 / 20.0;
            let strike = Dual2::<5>::from_value(forward * k.exp());
            let p = svi_durrleman_g(&p2, strike).min_val(Dual2::<5>::zero()) * 10.0;
            for a in 0..5 {
                for b in 0..5 {
                    h[(a, b)] += p.d1[a] * p.d1[b] + p.val * p.d2[a][b];
                }
            }
        }
        let floor = p2.a + p2.b * p2.sigma * (Dual2::<5>::one() - p2.rho * p2.rho).sqrt();
        let wing = p2.b * (Dual2::<5>::one() + p2.rho.abs()) - 2.0;
        for p in [
            (-floor).max_val(Dual2::<5>::zero()),
            wing.max_val(Dual2::<5>::zero()),
        ] {
            let p = p * 10.0;
            for a in 0..5 {
                for b in 0..5 {
                    h[(a, b)] += p.d1[a] * p.d1[b] + p.val * p.d2[a][b];
                }
            }
        }
        Ok(Self {
            smile,
            quotes: quotes.to_vec(),
            rmse,
            min_durrleman_g: min_g,
            sensitivity: IftSystemSensitivity::new(h, jq),
        })
    }
    pub fn vol_and_quote_vega(&self, strike: f64) -> (f64, Vec<f64>) {
        let p = AdSviSmile::new(
            Dual::<5>::from_value(self.smile.forward),
            Dual::<5>::from_value(self.smile.expiry),
            Dual::<5>::variable(self.smile.a, 0),
            Dual::<5>::variable(self.smile.b, 1),
            Dual::<5>::variable(self.smile.rho, 2),
            Dual::<5>::variable(self.smile.m, 3),
            Dual::<5>::variable(self.smile.sigma, 4),
        );
        let y = p.vol(Dual::<5>::from_value(strike));
        (y.val, self.sensitivity.vjp(&y.d))
    }
    pub fn implied_vol_ad<T: Real>(&self, strike: T, quotes: &[T]) -> T {
        let j = self.sensitivity.jacobian_dtheta_dq();
        let mut p = [
            T::from_f64(self.smile.a),
            T::from_f64(self.smile.b),
            T::from_f64(self.smile.rho),
            T::from_f64(self.smile.m),
            T::from_f64(self.smile.sigma),
        ];
        for a in 0..5 {
            for q in 0..self.quotes.len().min(quotes.len()) {
                p[a] += T::from_f64(j[(a, q)]) * (quotes[q] - self.quotes[q].market_vol);
            }
        }
        AdSviSmile::new(
            T::from_f64(self.smile.forward),
            T::from_f64(self.smile.expiry),
            p[0],
            p[1],
            p[2],
            p[3],
            p[4],
        )
        .vol(strike)
    }
    pub fn value_and_quote_vega<F>(&self, pricer: F) -> (f64, Vec<f64>)
    where
        F: Fn(&AdSviSmile<Dual<5>>) -> Dual<5>,
    {
        let p = AdSviSmile::new(
            Dual::from_value(self.smile.forward),
            Dual::from_value(self.smile.expiry),
            Dual::variable(self.smile.a, 0),
            Dual::variable(self.smile.b, 1),
            Dual::variable(self.smile.rho, 2),
            Dual::variable(self.smile.m, 3),
            Dual::variable(self.smile.sigma, 4),
        );
        let out = pricer(&p);
        (out.val, self.sensitivity.vjp(&out.d))
    }
}

#[derive(Clone, Debug)]
pub struct AdEssviSlice<T: Real> {
    pub forward: T,
    pub expiry: T,
    pub theta: T,
    pub rho: T,
    pub psi: T,
}
impl<T: Real> AdEssviSlice<T> {
    pub fn total_variance(&self, strike: T) -> T {
        let k = (strike / self.forward).ln();
        let phi = self.psi / self.theta.max_val(T::from_f64(1e-12));
        self.theta
            * 0.5
            * (T::one()
                + self.rho * phi * k
                + ((phi * k + self.rho) * (phi * k + self.rho) + T::one() - self.rho * self.rho)
                    .sqrt())
    }
    pub fn implied_vol(&self, strike: T) -> T {
        (self.total_variance(strike) / self.expiry.max_val(T::from_f64(1e-12)))
            .max_val(T::zero())
            .sqrt()
    }
}
#[derive(Clone, Debug)]
pub struct AdEssviSurface<T: Real> {
    pub slices: Vec<AdEssviSlice<T>>,
}
impl<T: Real> AdEssviSurface<T> {
    pub fn total_variance(&self, expiry: f64, strike: T) -> T {
        if self.slices.is_empty() {
            return T::zero();
        }
        let i = self
            .slices
            .partition_point(|s| s.expiry.val() < expiry)
            .min(self.slices.len().saturating_sub(1));
        if self.slices.len() < 2 || i == 0 {
            return self.slices[0].total_variance(strike);
        }
        let a = ((expiry - self.slices[i - 1].expiry.val())
            / (self.slices[i].expiry.val() - self.slices[i - 1].expiry.val()))
        .clamp(0.0, 1.0);
        let s0 = &self.slices[i - 1];
        let s1 = &self.slices[i];
        let theta = s0.theta * (1.0 - a) + s1.theta * a;
        let psi = s0.psi * (1.0 - a) + s1.psi * a;
        let rp = s0.rho * s0.psi * (1.0 - a) + s1.rho * s1.psi * a;
        let rho = rp / psi.max_val(T::from_f64(1e-12));
        AdEssviSlice {
            forward: s0.forward * (1.0 - a) + s1.forward * a,
            expiry: T::from_f64(expiry),
            theta,
            rho,
            psi,
        }
        .total_variance(strike)
    }
    pub fn implied_vol(&self, expiry: f64, strike: T) -> T {
        if expiry <= 0.0 || !expiry.is_finite() {
            return T::zero();
        }
        (self.total_variance(expiry, strike) / T::from_f64(expiry))
            .max_val(T::zero())
            .sqrt()
    }
    pub fn is_arbitrage_free(&self, k_grid: &[f64]) -> bool {
        if self.slices.is_empty()
            || self.slices.iter().any(|s| {
                !s.forward.val().is_finite()
                    || s.forward.val() <= 0.0
                    || !s.expiry.val().is_finite()
                    || s.expiry.val() <= 0.0
                    || !s.theta.val().is_finite()
                    || s.theta.val() <= 0.0
                    || !s.rho.val().is_finite()
                    || s.rho.val().abs() >= 1.0
                    || !s.psi.val().is_finite()
                    || s.psi.val() < 0.0
            })
        {
            return false;
        }
        self.slices.windows(2).all(|w| {
            w[1].theta.val() >= w[0].theta.val()
                && w[1].psi.val() >= w[0].psi.val()
                && (w[1].rho.val() * w[1].psi.val() - w[0].rho.val() * w[0].psi.val()).abs()
                    <= w[1].psi.val() - w[0].psi.val() + 1e-12
        }) && self.slices.iter().all(|s| {
            s.psi.val() * (1.0 + s.rho.val().abs()) <= 4.0 + 1e-12
                && s.psi.val() * s.psi.val() / s.theta.val().max(1e-12) * (1.0 + s.rho.val().abs())
                    <= 4.0 + 1e-12
                && k_grid.iter().all(|&k| {
                    let kk = s.forward * T::from_f64(k.exp());
                    let w = s.total_variance(kk);
                    w.val() > -1e-12
                })
        })
    }
}
#[derive(Clone, Debug)]
pub struct CalibratedEssviSurface {
    pub surface: AdEssviSurface<f64>,
    pub quotes: Vec<SmileStrikeQuote>,
    pub sensitivity: ql_ad::IftSystemSensitivity,
}
impl CalibratedEssviSurface {
    /// Calibrate each eSSVI slice to its quoted smile.  The three parameters
    /// (theta, rho, psi) are solved in volatility space; the normal-equation
    /// Hessian and quote Jacobian retained below provide implicit quote vegas.
    pub fn calibrate(
        slices: &[(f64, f64, Vec<SmileStrikeQuote>)],
    ) -> Result<Self, ql_math::MathError> {
        if slices.is_empty()
            || slices.iter().any(|(t, f, qs)| {
                !t.is_finite()
                    || *t <= 0.0
                    || !f.is_finite()
                    || *f <= 0.0
                    || qs.iter().any(|q| {
                        !q.strike.is_finite()
                            || q.strike <= 0.0
                            || !q.market_vol.is_finite()
                            || q.market_vol <= 0.0
                    })
            })
        {
            return Err(ql_math::MathError::MaxIterations);
        }
        if slices.windows(2).any(|w| w[1].0 <= w[0].0) {
            return Err(ql_math::MathError::MaxIterations);
        }
        let mut calibrated = Vec::with_capacity(slices.len());
        let mut quote_offsets = Vec::with_capacity(slices.len());
        let mut qv = Vec::new();
        for (t, f, qs) in slices {
            quote_offsets.push(qv.len());
            qv.extend(qs.iter().copied());
            let atm = qs
                .iter()
                .map(|q| q.market_vol * q.weight.max(0.0))
                .sum::<f64>()
                / qs.iter().map(|q| q.weight.max(0.0)).sum::<f64>().max(1e-12);
            let atm = if atm.is_finite() && atm > 1e-8 {
                atm
            } else {
                0.2
            };
            let x0 = DVector::from_vec(vec![
                (atm * atm * *t).max(1e-8),
                -0.2,
                (0.5 * atm * t.sqrt()).max(1e-8),
            ]);
            let qs_local = qs.clone();
            let residual = move |x: &DVector<f64>| {
                let theta = x[0].max(1e-10);
                let rho = x[1].clamp(-0.999, 0.999);
                let psi = x[2].max(1e-10);
                let ad = AdEssviSlice {
                    forward: Dual::<3>::from_value(*f),
                    expiry: Dual::<3>::from_value(*t),
                    theta: Dual::<3>::variable(theta, 0),
                    rho: Dual::<3>::variable(rho, 1),
                    psi: Dual::<3>::variable(psi, 2),
                };
                let mut r = DVector::zeros(qs_local.len());
                let mut j = DMatrix::zeros(qs_local.len(), 3);
                for (i, q) in qs_local.iter().enumerate() {
                    let w = q.weight.max(0.0).sqrt();
                    let y = ad.implied_vol(Dual::<3>::from_value(q.strike));
                    r[i] = w * (y.val - q.market_vol);
                    for a in 0..3 {
                        j[(i, a)] = w * y.d[a];
                    }
                }
                (r, j)
            };
            let sol = if qs.is_empty() {
                x0
            } else {
                levenberg_marquardt(x0, residual, 1e-11, 200)?
            };
            let prev_theta = calibrated
                .last()
                .map_or(1e-10, |p: &AdEssviSlice<f64>| p.theta);
            let prev_pos = calibrated
                .last()
                .map_or(0.0, |p: &AdEssviSlice<f64>| p.psi * (1.0 + p.rho));
            let prev_neg = calibrated
                .last()
                .map_or(0.0, |p: &AdEssviSlice<f64>| p.psi * (1.0 - p.rho));
            let theta = sol[0].max(prev_theta);
            let rho = sol[1].clamp(-0.999, 0.999);
            let psi_min = (prev_pos / (1.0 + rho).max(1e-6))
                .max(prev_neg / (1.0 - rho).max(1e-6))
                .max(1e-10);
            let psi_max = (4.0 / (1.0 + rho.abs()))
                .min((4.0 * theta / ((1.0 + rho.abs()) * (1.0 + rho.abs()))).sqrt())
                * 0.999;
            let psi = sol[2].max(psi_min).min(psi_max.max(psi_min));
            calibrated.push(AdEssviSlice {
                forward: *f,
                expiry: *t,
                theta,
                rho,
                psi,
            });
        }

        let n = calibrated.len() * 3;
        let mut h = DMatrix::zeros(n, n);
        let mut jq = DMatrix::zeros(n, qv.len());
        for (slice_idx, ((t, f, qs), s)) in slices.iter().zip(calibrated.iter()).enumerate() {
            let ad = AdEssviSlice {
                forward: Dual2::<3>::from_value(*f),
                expiry: Dual2::<3>::from_value(*t),
                theta: Dual2::<3>::variable(s.theta, 0),
                rho: Dual2::<3>::variable(s.rho, 1),
                psi: Dual2::<3>::variable(s.psi, 2),
            };
            let off = slice_idx * 3;
            let qoff = quote_offsets[slice_idx];
            for (qi, q) in qs.iter().enumerate() {
                let w = q.weight.max(0.0).sqrt();
                let y = ad.implied_vol(Dual2::<3>::from_value(q.strike));
                let rv = w * (y.val - q.market_vol);
                for a in 0..3 {
                    jq[(off + a, qoff + qi)] = -w * w * y.d1[a];
                    for b in 0..3 {
                        h[(off + a, off + b)] += w * w * y.d1[a] * y.d1[b] + w * rv * y.d2[a][b];
                    }
                }
            }
            // A tiny ridge makes an under-quoted slice well-conditioned while
            // leaving quote vegas unchanged to machine precision.
            for a in 0..3 {
                h[(off + a, off + a)] += 1e-12;
            }
        }
        Ok(Self {
            surface: AdEssviSurface { slices: calibrated },
            quotes: qv,
            sensitivity: ql_ad::IftSystemSensitivity::new(h, jq),
        })
    }
    pub fn vol_and_quote_vega(&self, expiry: f64, strike: f64) -> (f64, Vec<f64>) {
        let n = self.surface.slices.len() * 3;
        if n == 0 {
            return (0.0, Vec::new());
        }
        if self.quotes.is_empty() {
            return (self.surface.implied_vol(expiry, strike).val(), Vec::new());
        }
        let i = self
            .surface
            .slices
            .partition_point(|s| s.expiry < expiry)
            .min(self.surface.slices.len().saturating_sub(1));
        let mut g = vec![0.0; n];
        let value = if self.surface.slices.len() < 2 || i == 0 {
            let s = &self.surface.slices[0];
            let ad = AdEssviSlice {
                forward: Dual::<3>::from_value(s.forward),
                expiry: Dual::<3>::from_value(expiry.max(1e-12)),
                theta: Dual::<3>::variable(s.theta, 0),
                rho: Dual::<3>::variable(s.rho, 1),
                psi: Dual::<3>::variable(s.psi, 2),
            };
            let y = ad.implied_vol(Dual::<3>::from_value(strike));
            g[..3].copy_from_slice(&y.d);
            y.val
        } else {
            let s0 = &self.surface.slices[i - 1];
            let s1 = &self.surface.slices[i];
            let a = ((expiry - s0.expiry) / (s1.expiry - s0.expiry)).clamp(0.0, 1.0);
            let ad0 = AdEssviSlice {
                forward: Dual::<6>::from_value(s0.forward),
                expiry: Dual::<6>::from_value(expiry),
                theta: Dual::<6>::variable(s0.theta, 0),
                rho: Dual::<6>::variable(s0.rho, 1),
                psi: Dual::<6>::variable(s0.psi, 2),
            };
            let ad1 = AdEssviSlice {
                forward: Dual::<6>::from_value(s1.forward),
                expiry: Dual::<6>::from_value(expiry),
                theta: Dual::<6>::variable(s1.theta, 3),
                rho: Dual::<6>::variable(s1.rho, 4),
                psi: Dual::<6>::variable(s1.psi, 5),
            };
            let theta = ad0.theta * (1.0 - a) + ad1.theta * a;
            let psi = ad0.psi * (1.0 - a) + ad1.psi * a;
            let rp = ad0.rho * ad0.psi * (1.0 - a) + ad1.rho * ad1.psi * a;
            let rho = rp / psi.max_val(Dual::<6>::from_value(1e-12));
            let fwd = ad0.forward * (1.0 - a) + ad1.forward * a;
            let y = AdEssviSlice {
                forward: fwd,
                expiry: Dual::<6>::from_value(expiry),
                theta,
                rho,
                psi,
            }
            .implied_vol(Dual::<6>::from_value(strike));
            let off = (i - 1) * 3;
            g[off..off + 6].copy_from_slice(&y.d);
            y.val
        };
        (value, self.sensitivity.vjp(&g))
    }
}

#[derive(Copy, Clone, Debug)]
pub struct GatheralLsvGyongyLeverageSpec<T: Real = f64> {
    pub forward: T,
    pub strike: T,
    pub expiry: T,
    pub svi_a: T,
    pub svi_b: T,
    pub svi_rho: T,
    pub svi_m: T,
    pub svi_sigma: T,
    pub dw_dt_atm: T,
    pub heston_v0: T,
    pub heston_theta: T,
    pub heston_kappa: T,
    pub heston_eta: T,
    pub heston_rho: T,
}

#[derive(Copy, Clone, Debug)]
pub struct GatheralLsvGyongyLeverageReport<T: Real = f64> {
    pub implied_vol: T,
    pub dupire_local_vol: T,
    pub conditional_stochastic_vol: T,
    pub gyongy_leverage_ratio: T,
    pub geodesic_implied_vol: T,
    pub local_to_implied_atm_skew_ratio: T,
    pub durrleman_g: T,
}

pub fn evaluate_gatheral_gyongy_lsv_leverage<T: Real>(
    s: &GatheralLsvGyongyLeverageSpec<T>,
) -> GatheralLsvGyongyLeverageReport<T> {
    let y = (s.strike / s.forward).ln();
    let z = y - s.svi_m;
    let root = (z * z + s.svi_sigma * s.svi_sigma).sqrt();
    let w = s.svi_a + s.svi_b * (s.svi_rho * z + root);
    let w0 = s.svi_a
        + s.svi_b * (-s.svi_rho * s.svi_m + (s.svi_m * s.svi_m + s.svi_sigma * s.svi_sigma).sqrt());
    let wp = s.svi_b * (s.svi_rho + z / root);
    let wpp = s.svi_b * s.svi_sigma * s.svi_sigma / (root * root * root);
    let g = (T::one() - y * wp / (w.max_val(T::from_f64(1e-12)) * 2.0)).powf(2.0)
        - wp * wp / 4.0 * (T::one() / w.max_val(T::from_f64(1e-12)) + 0.25)
        + wpp * 0.5;
    let dw = s.dw_dt_atm * w / w0.max_val(T::from_f64(1e-12));
    let loc = (dw / g.max_val(T::from_f64(1e-8)))
        .max_val(T::from_f64(1e-12))
        .sqrt();

    let lam = s.heston_kappa.max_val(T::from_f64(1e-8));
    let t = s.expiry.max_val(T::from_f64(1e-10));
    let ev = (-lam * t).exp();
    let bv = s.heston_theta + (s.heston_v0 - s.heston_theta) * ev;
    let bw = s.heston_theta * t + (s.heston_v0 - s.heston_theta) * (T::one() - ev) / lam;
    let m1 = (T::one() - ev) / (lam * t);
    let m2 = (T::one() - (-lam * t * 2.0).exp()) / (lam * t * 2.0);
    let hv = bv
        * (s.heston_rho * s.heston_eta / (lam * bw.max_val(T::from_f64(1e-12)))
            * (T::one() - m1)
            * (y + bw * 0.5)
            + s.heston_eta * s.heston_eta * (T::one() - s.heston_rho * s.heston_rho)
                / (lam * lam * bw * 2.0).max_val(T::from_f64(1e-12))
                * (T::one() - m1 * 2.0 + m2)
                * (y * y / bw.max_val(T::from_f64(1e-12)) - T::one()))
        .exp();
    let csv = hv.max_val(T::from_f64(1e-12)).sqrt();
    let imp = (w / t).max_val(T::from_f64(1e-12)).sqrt();
    let nodes = [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116];
    let weights = [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451];
    let loc_at = |yy: T| -> T {
        let zz = yy - s.svi_m;
        let rr = (zz * zz + s.svi_sigma * s.svi_sigma).sqrt();
        let ww = s.svi_a + s.svi_b * (s.svi_rho * zz + rr);
        let w1 = s.svi_b * (s.svi_rho + zz / rr);
        let w2 = s.svi_b * s.svi_sigma * s.svi_sigma / (rr * rr * rr);
        let gd = (T::one() - yy * w1 / (ww.max_val(T::from_f64(1e-12)) * 2.0)).powf(2.0)
            - w1 * w1 / 4.0 * (T::one() / ww.max_val(T::from_f64(1e-12)) + 0.25)
            + w2 * 0.5;
        (s.dw_dt_atm * ww / w0.max_val(T::from_f64(1e-12)) / gd.max_val(T::from_f64(1e-8)))
            .max_val(T::from_f64(1e-12))
            .sqrt()
    };
    let geo_at = |yy: T| -> T {
        let mut inv = T::zero();
        for (&node, &weight) in nodes.iter().zip(weights.iter()) {
            let u_y = yy * T::from_f64(0.5 * (node + 1.0));
            inv += T::from_f64(0.5 * weight) / loc_at(u_y);
        }
        T::one() / inv.max_val(T::from_f64(1e-12))
    };
    let geodesic = if y.val().abs() < 1e-12 {
        loc
    } else {
        geo_at(y)
    };
    let eps = T::from_f64(1e-6);
    let dloc_atm = (loc_at(eps) - loc_at(-eps)) / (eps * 2.0);
    let dgeo_atm = (geo_at(eps) - geo_at(-eps)) / (eps * 2.0);
    let skew_ratio = if dgeo_atm.val().abs() > 1e-10 {
        dloc_atm / dgeo_atm
    } else {
        T::from_f64(2.0)
    };
    GatheralLsvGyongyLeverageReport {
        implied_vol: imp,
        dupire_local_vol: loc,
        conditional_stochastic_vol: csv,
        gyongy_leverage_ratio: loc / csv,
        geodesic_implied_vol: geodesic,
        local_to_implied_atm_skew_ratio: skew_ratio,
        durrleman_g: g,
    }
}

#[cfg(test)]
mod phase13_tests {
    use super::*;
    #[test]
    fn essvi_three_slice_vega() {
        let q = |v| {
            vec![
                SmileStrikeQuote {
                    strike: 90.0,
                    market_vol: v + 0.02,
                    weight: 1.0,
                },
                SmileStrikeQuote {
                    strike: 100.0,
                    market_vol: v,
                    weight: 1.0,
                },
                SmileStrikeQuote {
                    strike: 110.0,
                    market_vol: v - 0.01,
                    weight: 1.0,
                },
            ]
        };
        let c = CalibratedEssviSurface::calibrate(&[
            (0.25, 100., q(0.2)),
            (0.5, 100., q(0.21)),
            (1., 100., q(0.22)),
        ])
        .unwrap();
        let (_, g) = c.vol_and_quote_vega(0.75, 100.);
        assert!(g.len() == 9);
    }
}

#[cfg(test)]
mod phase21_tests {
    use super::*;

    #[test]
    fn gatheral_atm_geodesic_collapse_and_skew_ratio() {
        let s = GatheralLsvGyongyLeverageSpec {
            forward: 100.0,
            strike: 100.0,
            expiry: 1.0,
            svi_a: 0.04,
            svi_b: 0.2,
            svi_rho: -0.3,
            svi_m: 0.0,
            svi_sigma: 0.2,
            dw_dt_atm: 0.04,
            heston_v0: 0.04,
            heston_theta: 0.04,
            heston_kappa: 1.0,
            heston_eta: 0.5,
            heston_rho: -0.7,
        };
        let r = evaluate_gatheral_gyongy_lsv_leverage(&s);
        assert!((r.geodesic_implied_vol - r.dupire_local_vol).abs() < 1e-10);
        assert!((r.local_to_implied_atm_skew_ratio - 2.0).abs() < 1e-10);
        assert!(r.durrleman_g.val() > 0.0);
    }
}

#[derive(Clone, Debug)]
pub struct HestonNandiTsayDccGarchSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub daily_risk_free_rate: T,
    pub days_to_maturity: usize,
    pub initial_daily_variance: T,
    pub omega: T,
    pub alpha: T,
    pub beta: T,
    pub gamma_asymmetry: T,
    pub asset1_returns: Vec<T>,
    pub asset2_returns: Vec<T>,
    pub dcc_a: T,
    pub dcc_b: T,
}
#[derive(Clone, Debug)]
pub struct HestonNandiTsayDccGarchReport<T: Real = f64> {
    pub garch_persistence: T,
    pub long_run_annualized_vol: T,
    pub martingale_mgf_error: T,
    pub heston_nandi_call_price: T,
    pub heston_nandi_put_price: T,
    pub symmetric_garch_call_price: T,
    pub leverage_skew_impact: T,
    pub terminal_dcc_correlation: T,
    pub max_dcc_correlation: T,
    pub min_dcc_correlation: T,
}
type C<T> = (T, T);
fn c_add<T: Real>(a: C<T>, b: C<T>) -> C<T> {
    (a.0 + b.0, a.1 + b.1)
}
fn c_sub<T: Real>(a: C<T>, b: C<T>) -> C<T> {
    (a.0 - b.0, a.1 - b.1)
}
fn c_mul<T: Real>(a: C<T>, b: C<T>) -> C<T> {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}
fn c_scale<T: Real>(a: C<T>, x: T) -> C<T> {
    (a.0 * x, a.1 * x)
}
fn c_div<T: Real>(a: C<T>, b: C<T>) -> C<T> {
    let den = (b.0 * b.0 + b.1 * b.1).max_val(T::from_f64(1e-18));
    ((a.0 * b.0 + a.1 * b.1) / den, (a.1 * b.0 - a.0 * b.1) / den)
}
fn c_ln<T: Real>(z: C<T>) -> C<T> {
    let mag = (z.0 * z.0 + z.1 * z.1)
        .max_val(T::from_f64(1e-24))
        .sqrt()
        .ln();
    let arg = if z.0.val() > 0.0 {
        (z.1 / z.0).atan()
    } else if z.1.val() >= 0.0 {
        z.1.abs().atan() + T::from_f64(std::f64::consts::FRAC_PI_2)
    } else {
        -(z.1.abs().atan() + T::from_f64(std::f64::consts::FRAC_PI_2))
    };
    (mag, arg)
}
fn c_exp<T: Real>(z: C<T>) -> C<T> {
    let e = z.0.exp();
    (e * z.1.cos(), e * z.1.sin())
}

fn heston_nandi_ab<T: Real>(
    s: &HestonNandiTsayDccGarchSpec<T>,
    gamma_star: T,
    phi: C<T>,
) -> (C<T>, C<T>) {
    let mut a = (T::zero(), T::zero());
    let mut b = (T::zero(), T::zero());
    let one = (T::one(), T::zero());
    for _ in 0..s.days_to_maturity.max(1) {
        let den = c_sub(one, c_scale(b, T::from_f64(2.0) * s.alpha));
        let term = c_div(
            c_scale(
                c_mul(
                    c_sub(phi, (gamma_star, T::zero())),
                    c_sub(phi, (gamma_star, T::zero())),
                ),
                T::from_f64(0.5),
            ),
            den,
        );
        // The first term is phi*(gamma* - 1/2) - gamma*^2/2.
        let first = c_sub(
            c_scale(phi, gamma_star - T::from_f64(0.5)),
            (T::from_f64(0.5) * gamma_star * gamma_star, T::zero()),
        );
        let b_next = c_add(c_add(first, c_scale(b, s.beta)), term);
        let a_next = c_add(
            a,
            c_sub(
                c_add(
                    c_scale(phi, s.daily_risk_free_rate),
                    c_scale(b_next, s.omega),
                ),
                c_scale(c_ln(den), T::from_f64(0.5)),
            ),
        );
        b = b_next;
        a = a_next;
    }
    (a, b)
}

fn hn_characteristic<T: Real>(s: &HestonNandiTsayDccGarchSpec<T>, gamma: T, phi: C<T>) -> C<T> {
    let (a, b) = heston_nandi_ab(s, gamma, phi);
    c_exp(c_add(
        c_scale(phi, s.spot.ln()),
        c_add(a, c_scale(b, s.initial_daily_variance)),
    ))
}

fn hn_call<T: Real>(s: &HestonNandiTsayDccGarchSpec<T>, gamma: T) -> T {
    const NODES: [f64; 24] = [
        -0.9951872199970213,
        -0.9747285559713095,
        -0.9382745520027328,
        -0.8864155270044011,
        -0.820001985973903,
        -0.7401241915785544,
        -0.6480936519369755,
        -0.5454214713888396,
        -0.4337935076260452,
        -0.3150426796961634,
        -0.1911188674736163,
        -0.06405689286260563,
        0.06405689286260563,
        0.1911188674736163,
        0.3150426796961634,
        0.4337935076260452,
        0.5454214713888396,
        0.6480936519369755,
        0.7401241915785544,
        0.820001985973903,
        0.8864155270044011,
        0.9382745520027328,
        0.9747285559713095,
        0.9951872199970213,
    ];
    const WEIGHTS: [f64; 24] = [
        0.012341229799987334,
        0.028531388628933813,
        0.044277438817419676,
        0.05929858491543666,
        0.07334648141108027,
        0.08619016153195322,
        0.0976186521041139,
        0.10744427011596562,
        0.11550566805372561,
        0.12167047292780335,
        0.12583745634682839,
        0.12793819534675224,
        0.12793819534675224,
        0.12583745634682839,
        0.12167047292780335,
        0.11550566805372561,
        0.10744427011596562,
        0.0976186521041139,
        0.08619016153195322,
        0.07334648141108027,
        0.05929858491543666,
        0.044277438817419676,
        0.028531388628933813,
        0.012341229799987334,
    ];
    let klog = s.strike.ln();
    let (a1, b1) = heston_nandi_ab(s, gamma, (T::one(), T::zero()));
    let m1 = c_exp(c_add(
        c_scale((T::one(), T::zero()), s.spot.ln()),
        c_add(a1, c_scale(b1, s.initial_daily_variance)),
    ));
    let mut p1 = T::from_f64(0.5);
    let mut p2 = T::from_f64(0.5);
    for j in 0..24 {
        let u = T::from_f64(20.0 * (NODES[j] + 1.0));
        let wt = T::from_f64(20.0 * WEIGHTS[j]);
        let phase = c_exp((T::zero(), -u * klog));
        let den = (T::zero(), u);
        let c2 = c_div(
            c_mul(phase, hn_characteristic(s, gamma, (T::zero(), u))),
            den,
        );
        let c1 = c_div(
            c_mul(phase, c_div(hn_characteristic(s, gamma, (T::one(), u)), m1)),
            den,
        );
        p1 += wt * c1.0 / T::from_f64(std::f64::consts::PI);
        p2 += wt * c2.0 / T::from_f64(std::f64::consts::PI);
    }
    s.spot * p1
        - s.strike
            * (-T::from_f64(s.days_to_maturity.max(1) as f64) * s.daily_risk_free_rate).exp()
            * p2
}

pub fn evaluate_heston_nandi_and_tsay_dcc_garch<T: Real>(
    s: &HestonNandiTsayDccGarchSpec<T>,
) -> HestonNandiTsayDccGarchReport<T> {
    let n = s.days_to_maturity.max(1);
    let pers = s.beta + s.alpha * s.gamma_asymmetry * s.gamma_asymmetry;
    let lr = (T::from_f64(252.0) * (s.omega + s.alpha)
        / (T::one() - pers).max_val(T::from_f64(1e-8)))
    .sqrt();
    let call = hn_call(s, s.gamma_asymmetry);
    let sym = hn_call(s, T::zero());
    let n_t = T::from_f64(n as f64);
    let put = call - s.spot + s.strike * (-n_t * s.daily_risk_free_rate).exp();
    let (a1, b1) = heston_nandi_ab(s, s.gamma_asymmetry, (T::one(), T::zero()));
    let martingale_mgf_error = b1.0.abs() + (a1.0 - n_t * s.daily_risk_free_rate).abs();
    let m = s.asset1_returns.len().min(s.asset2_returns.len());
    let mut h1 = s.initial_daily_variance;
    let mut h2 = s.initial_daily_variance;
    let mut e1s = Vec::with_capacity(m);
    let mut e2s = Vec::with_capacity(m);
    for i in 0..m {
        e1s.push(s.asset1_returns[i] / h1.max_val(T::from_f64(1e-12)).sqrt());
        e2s.push(s.asset2_returns[i] / h2.max_val(T::from_f64(1e-12)).sqrt());
        h1 = s.omega
            + s.alpha * (s.asset1_returns[i] - s.gamma_asymmetry * h1.sqrt()).powf(2.0)
            + s.beta * h1;
        h2 = s.omega
            + s.alpha * (s.asset2_returns[i] - s.gamma_asymmetry * h2.sqrt()).powf(2.0)
            + s.beta * h2;
    }
    let cross = e1s
        .iter()
        .zip(&e2s)
        .fold(T::zero(), |z, (x, y)| z + *x * *y);
    let ss1 = e1s.iter().fold(T::zero(), |z, x| z + *x * *x);
    let ss2 = e2s.iter().fold(T::zero(), |z, x| z + *x * *x);
    let rho_bar = if m > 0 {
        cross / (ss1 * ss2).max_val(T::from_f64(1e-12)).sqrt()
    } else {
        T::zero()
    };
    let a = s.dcc_a.max_val(T::zero());
    let b = s.dcc_b.max_val(T::zero());
    let mut q11 = T::one();
    let mut q22 = T::one();
    let mut q12 = rho_bar;
    let mut max = T::from_f64(-1.);
    let mut min = T::from_f64(1.);
    for i in 0..m {
        q11 = (T::one() - a - b) + a * e1s[i] * e1s[i] + b * q11;
        q22 = (T::one() - a - b) + a * e2s[i] * e2s[i] + b * q22;
        q12 = (T::one() - a - b) * rho_bar + a * e1s[i] * e2s[i] + b * q12;
        let rr = (q12 / (q11 * q22).max_val(T::from_f64(1e-12)).sqrt())
            .max_val(T::from_f64(-0.999999))
            .min_val(T::from_f64(0.999999));
        max = max.max_val(rr);
        min = min.min_val(rr);
    }
    let term = if m > 0 {
        q12 / (q11 * q22).max_val(T::from_f64(1e-12)).sqrt()
    } else {
        rho_bar
    };
    HestonNandiTsayDccGarchReport {
        garch_persistence: pers,
        long_run_annualized_vol: lr,
        martingale_mgf_error,
        heston_nandi_call_price: call,
        heston_nandi_put_price: put,
        symmetric_garch_call_price: sym,
        leverage_skew_impact: call - sym,
        terminal_dcc_correlation: term,
        max_dcc_correlation: max,
        min_dcc_correlation: min,
    }
}

#[cfg(test)]
mod phase24_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(x: T) -> HestonNandiTsayDccGarchSpec<T> {
        HestonNandiTsayDccGarchSpec {
            spot: x,
            strike: T::from_f64(100.),
            daily_risk_free_rate: T::from_f64(0.0001),
            days_to_maturity: 30,
            initial_daily_variance: T::from_f64(0.0002),
            omega: T::from_f64(0.00001),
            alpha: T::from_f64(0.05),
            beta: T::from_f64(0.9),
            gamma_asymmetry: T::from_f64(0.2),
            asset1_returns: vec![T::from_f64(0.01), T::from_f64(-0.005), T::from_f64(0.002)],
            asset2_returns: vec![T::from_f64(0.008), T::from_f64(-0.004), T::from_f64(0.003)],
            dcc_a: T::from_f64(0.02),
            dcc_b: T::from_f64(0.95),
        }
    }
    #[test]
    fn phase24() {
        let r = evaluate_heston_nandi_and_tsay_dcc_garch(&mk(100.));
        assert!(r.martingale_mgf_error.val() < 1e-12);
        assert!(r.terminal_dcc_correlation.val().abs() < 1.);
        let h = 1e-5;
        let fd = (evaluate_heston_nandi_and_tsay_dcc_garch(&mk(100. + h)).heston_nandi_call_price
            - evaluate_heston_nandi_and_tsay_dcc_garch(&mk(100. - h)).heston_nandi_call_price)
            / (2. * h);
        let ad = evaluate_heston_nandi_and_tsay_dcc_garch(&mk(Dual::<1>::variable(100., 0)))
            .heston_nandi_call_price;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BergomiForwardStartNapoleonVolSpec<T: Real = f64> {
    pub initial_forward_variance: T,
    pub vol_of_vol_omega: T,
    pub kappa_short_k1: T,
    pub kappa_long_k2: T,
    pub long_factor_weight_theta: T,
    pub factor_corr_rho12: T,
    pub spot_vol_corr_rho_s1: T,
    pub spot_vol_corr_rho_s2: T,
    pub forward_start_t1: T,
    pub forward_end_t2: T,
    pub forward_moneyness_k: T,
    pub napoleon_base_coupon: T,
    pub napoleon_global_floor: T,
    pub napoleon_subperiods_m: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct BergomiForwardStartNapoleonVolReport<T: Real = f64> {
    pub forward_variance_log_vol: T,
    pub forward_atm_implied_vol: T,
    pub forward_start_atm_skew: T,
    pub forward_start_smile_vol_at_k: T,
    pub forward_start_call_price_pct: T,
    pub constant_vol_napoleon_coupon: T,
    pub bergomi_stoch_vol_napoleon_coupon: T,
    pub napoleon_vol_of_vol_convexity_value: T,
}
pub fn evaluate_bergomi_forward_start_smile_and_napoleon_convexity<T: Real>(
    s: &BergomiForwardStartNapoleonVolSpec<T>,
) -> BergomiForwardStartNapoleonVolReport<T> {
    let dt = (s.forward_end_t2 - s.forward_start_t1).max_val(T::from_f64(1e-8));
    let k1 = s.kappa_short_k1.abs().max_val(T::from_f64(1e-8));
    let k2 = s.kappa_long_k2.abs().max_val(T::from_f64(1e-8));
    let a = (T::one() - (-k1 * dt).exp()) / k1;
    let b = (T::one() - (-k2 * dt).exp()) / k2;
    let theta = s.long_factor_weight_theta;
    let rho = s.factor_corr_rho12;
    let iv = (s.initial_forward_variance.abs()
        + s.vol_of_vol_omega
            * s.vol_of_vol_omega
            * (a * a + theta * theta * b * b + T::from_f64(2.) * theta * rho * a * b)
            / dt)
        .sqrt();
    let logvol = (s.vol_of_vol_omega.abs() * (a + theta * b) / dt.sqrt()).abs();
    let atm = iv;
    let skew = (s.spot_vol_corr_rho_s1 * a + s.spot_vol_corr_rho_s2 * theta * b)
        * s.vol_of_vol_omega
        / (dt.sqrt() * iv.max_val(T::from_f64(1e-8)));
    let smile = (atm + skew * (s.forward_moneyness_k - T::one())).max_val(T::from_f64(1e-8));
    let d1 = ((T::one() / s.forward_moneyness_k.max_val(T::from_f64(1e-8))).ln()
        + T::from_f64(0.5) * smile * smile * dt)
        / (smile * dt.sqrt());
    let call = (d1.norm_cdf() - s.forward_moneyness_k * (d1 - smile * dt.sqrt()).norm_cdf())
        .max_val(T::zero());
    let const_coupon = s.napoleon_base_coupon.max_val(s.napoleon_global_floor);
    let mut stoch = T::zero();
    let m = s.napoleon_subperiods_m.max(1);
    // Normalized nine-point weights (sum of tabulated probabilities).
    let ws = [
        0.000039, 0.004089, 0.027733, 0.078742, 0.113019, 0.078742, 0.027733, 0.004089, 0.000039,
    ];
    let zs = [
        -3.190993, -2.266581, -1.468553, -0.723551, 0., 0.723551, 1.468553, 2.266581, 3.190993,
    ];
    let sw: f64 = ws.iter().sum();
    let mut accumulated_weight = T::zero();
    for j in 0..8 {
        let shock = s.vol_of_vol_omega.abs() * T::from_f64(zs[j]) / T::from_f64(m as f64).sqrt();
        let coupon = s.napoleon_base_coupon + shock * dt.sqrt();
        let weight = T::from_f64(ws[j] / sw);
        accumulated_weight += weight;
        stoch += weight * coupon.max_val(s.napoleon_global_floor);
    }
    let shock_last = s.vol_of_vol_omega.abs() * T::from_f64(zs[8]) / T::from_f64(m as f64).sqrt();
    let coupon_last =
        (s.napoleon_base_coupon + shock_last * dt.sqrt()).max_val(s.napoleon_global_floor);
    stoch += (T::one() - accumulated_weight) * coupon_last;
    let stoch = stoch.max_val(s.napoleon_global_floor);
    BergomiForwardStartNapoleonVolReport {
        forward_variance_log_vol: logvol,
        forward_atm_implied_vol: atm,
        forward_start_atm_skew: skew,
        forward_start_smile_vol_at_k: smile,
        forward_start_call_price_pct: call,
        constant_vol_napoleon_coupon: const_coupon,
        bergomi_stoch_vol_napoleon_coupon: stoch,
        napoleon_vol_of_vol_convexity_value: (stoch - const_coupon).max_val(T::zero()),
    }
}

#[cfg(test)]
mod phase26_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(omega: T) -> BergomiForwardStartNapoleonVolSpec<T> {
        BergomiForwardStartNapoleonVolSpec {
            initial_forward_variance: T::from_f64(0.04),
            vol_of_vol_omega: omega,
            kappa_short_k1: T::from_f64(3.0),
            kappa_long_k2: T::from_f64(0.5),
            long_factor_weight_theta: T::from_f64(0.5),
            factor_corr_rho12: T::from_f64(0.2),
            spot_vol_corr_rho_s1: T::from_f64(-0.5),
            spot_vol_corr_rho_s2: T::from_f64(-0.3),
            forward_start_t1: T::from_f64(0.5),
            forward_end_t2: T::one(),
            forward_moneyness_k: T::one(),
            napoleon_base_coupon: T::from_f64(0.01),
            napoleon_global_floor: T::zero(),
            napoleon_subperiods_m: 4,
        }
    }
    #[test]
    fn bergomi_invariants_and_dual_fd() {
        let z = evaluate_bergomi_forward_start_smile_and_napoleon_convexity(&spec(0.0));
        assert!(z.napoleon_vol_of_vol_convexity_value.abs() < 1e-12);
        let r = evaluate_bergomi_forward_start_smile_and_napoleon_convexity(&spec(0.3));
        assert!(r.forward_start_atm_skew < 0.0);
        let h = 1e-5;
        let fd = (evaluate_bergomi_forward_start_smile_and_napoleon_convexity(&spec(0.3 + h))
            .forward_atm_implied_vol
            - evaluate_bergomi_forward_start_smile_and_napoleon_convexity(&spec(0.3 - h))
                .forward_atm_implied_vol)
            / (2.0 * h);
        let ad = evaluate_bergomi_forward_start_smile_and_napoleon_convexity(&spec(
            Dual::<1>::variable(0.3, 0),
        ))
        .forward_atm_implied_vol;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
