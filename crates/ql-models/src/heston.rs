use crate::Complex;
use ql_ad::Real;
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct HestonParams<T: Real> {
    pub v0: T,
    pub kappa: T,
    pub theta: T,
    pub sigma_v: T,
    pub rho: T,
}
impl<T: Real> HestonParams<T> {
    pub fn feller_ratio(&self) -> T {
        if self.sigma_v.val().abs() < 1e-14 {
            T::from_f64(f64::INFINITY)
        } else {
            self.kappa * self.theta * 2. / (self.sigma_v * self.sigma_v)
        }
    }
    pub fn satisfies_feller(&self) -> bool {
        self.feller_ratio().val() >= 1.
    }
    pub fn char_func_unit_spot(&self, u: Complex<T>, rate: T, div: T, tau: T) -> Complex<T> {
        let i = Complex::i();
        if self.sigma_v.val().abs() < 1e-12 {
            let kappa_tau = self.kappa * tau;
            let integrated_variance = if self.kappa.val().abs() < 1e-12 {
                self.v0 * tau
            } else {
                self.theta * tau
                    + (self.v0 - self.theta) * (T::one() - (-kappa_tau).exp()) / self.kappa
            };
            let mean = (rate - div) * tau - T::from_f64(0.5) * integrated_variance;
            return (i * u * mean
                - Complex::from_real(T::from_f64(0.5) * integrated_variance) * u * u)
                .exp();
        }
        // Andersen--Lake cancellation-free form with the Little Heston Trap
        // logarithm.  Rationalising xi-d avoids catastrophic cancellation for
        // large |u| and the branch-stable y representation avoids 1-g poles.
        let xi = Complex::from_real(self.kappa) - i * self.rho * self.sigma_v * u;
        let sig2 = self.sigma_v * self.sigma_v;
        let u2iu = u * u + i * u;
        let d = (xi * xi + Complex::from_real(sig2) * u2iu).sqrt();
        let dot = xi.re.val() * d.re.val() + xi.im.val() * d.im.val();
        let r = if dot > 0.0 {
            (Complex::from_real(-sig2) * u2iu) / (xi + d)
        } else {
            xi - d
        };
        let ed = (-d * tau).exp();
        let y = (ed - Complex::one()) / (d * T::from_f64(2.0));
        let one_minus_ry = Complex::one() - r * y;
        let a = Complex::from_real(self.kappa * self.theta / sig2)
            * (r * tau - one_minus_ry.ln() * T::from_f64(2.0));
        let b = (u2iu * y) / one_minus_ry;
        (i * u * (rate - div) * tau + a + b * self.v0).exp()
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct DoubleHestonParams<T: Real> {
    pub factor1: HestonParams<T>,
    pub factor2: HestonParams<T>,
}
impl<T: Real> DoubleHestonParams<T> {
    pub fn satisfies_feller(&self) -> bool {
        self.factor1.satisfies_feller() && self.factor2.satisfies_feller()
    }
    pub fn char_func_unit_spot(&self, u: Complex<T>, rate: T, div: T, tau: T) -> Complex<T> {
        let i = Complex::i();
        (i * u * (rate - div) * tau).exp()
            * self
                .factor1
                .char_func_unit_spot(u, T::zero(), T::zero(), tau)
            * self
                .factor2
                .char_func_unit_spot(u, T::zero(), T::zero(), tau)
    }
}
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BatesParams<T: Real> {
    pub heston: HestonParams<T>,
    pub lambda: T,
    pub mu_j: T,
    pub sigma_j: T,
}
impl<T: Real> BatesParams<T> {
    pub fn char_func_unit_spot(&self, u: Complex<T>, rate: T, div: T, tau: T) -> Complex<T> {
        let i = Complex::i();
        let jump = (i * u * self.mu_j
            - (Complex::from_real(self.sigma_j * self.sigma_j) * u * u) * T::from_f64(0.5))
        .exp();
        let comp = (self.mu_j + self.sigma_j * self.sigma_j * 0.5).exp() - 1.;
        self.heston.char_func_unit_spot(u, rate, div, tau)
            * (Complex::from_real(self.lambda * tau) * (jump - Complex::one() - i * u * comp)).exp()
    }
}
