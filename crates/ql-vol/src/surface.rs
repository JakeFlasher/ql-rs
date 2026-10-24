use nalgebra::{DMatrix, DVector};
use ql_ad::{Dual2, IftSystemSensitivity, Real};
use ql_math::{MathError, levenberg_marquardt};
pub use volsurf::smile::{SabrSmile, SviSmile};
pub use volsurf::surface::{EssviSurface, SsviSurface};
pub use volsurf::{
    DupireLocalVol as DupireTotalVarianceFormula, SmileSection, Strike, Tenor, Vol, VolSurface,
};
pub type DupireLocalVol = DupireTotalVarianceFormula;
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
                .all(|&k| svi_durrleman_g(self, T::from_f64(k)).val() >= -1e-9)
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
                for k in 0..3 {
                    let mut v = [p.alpha, p.rho, p.nu];
                    let h = 1e-5;
                    v[k] += h;
                    let up =
                        AdSabrSmile::new(forward, expiry, v[0], beta, v[1], v[2]).vol(q.strike);
                    v[k] -= 2. * h;
                    let dn =
                        AdSabrSmile::new(forward, expiry, v[0], beta, v[1], v[2]).vol(q.strike);
                    j[(i, k)] = d * (up - dn) / (2. * h);
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
                jq[(a, i)] = -d * y.d1[a];
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
        (p.vol(strike), self.sensitivity.vjp(&d[..3]))
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
            let mut r = DVector::zeros(quotes.len());
            let mut j = DMatrix::zeros(quotes.len(), 5);
            for (i, q) in quotes.iter().enumerate() {
                let d = q.weight.max(0.0).sqrt();
                r[i] = d * (smile.vol(q.strike) - q.market_vol);
                let h = 1e-5;
                let base = [x[0], x[1], x[2], x[3], x[4]];
                for k in 0..5 {
                    let mut up = base;
                    let mut dn = base;
                    up[k] += h;
                    dn[k] -= h;
                    let pu = AdSviSmile::new(
                        forward,
                        expiry,
                        up[0],
                        up[1].max(1e-8),
                        up[2].clamp(-0.999, 0.999),
                        up[3],
                        up[4].abs().max(1e-8),
                    )
                    .vol(q.strike);
                    let pd = AdSviSmile::new(
                        forward,
                        expiry,
                        dn[0],
                        dn[1].max(1e-8),
                        dn[2].clamp(-0.999, 0.999),
                        dn[3],
                        dn[4].abs().max(1e-8),
                    )
                    .vol(q.strike);
                    j[(i, k)] = d * (pu - pd) / (2. * h);
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
        let rmse = (r.norm_squared() / quotes.len().max(1) as f64).sqrt();
        let grid: Vec<f64> = quotes.iter().map(|q| (q.strike / forward).ln()).collect();
        let min_g = grid
            .iter()
            .map(|&k| svi_durrleman_g(&smile, forward * k.exp()))
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
                jq[(a, i)] = -d * y.d1[a];
                for b in 0..5 {
                    h[(a, b)] += d * d * y.d1[a] * y.d1[b] + d * r[i] * y.d2[a][b];
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
        let d = svi_quote_sensitivities(&self.smile, strike);
        (self.smile.vol(strike), self.sensitivity.vjp(&d))
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
}
