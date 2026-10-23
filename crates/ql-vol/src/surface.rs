use ql_ad::Real;
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
        let ratio = if z.val().abs() < 1e-6 {
            T::one() - self.rho * z / T::from_f64(2.0)
                + (T::from_f64(2.0) - T::from_f64(3.0) * self.rho * self.rho) * z * z
                    / T::from_f64(12.0)
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
    let ratio = if z.v.val().abs() < 1e-6 {
        add(
            one,
            add(
                scale(mul(rho, z), T::from_f64(-0.5)),
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
