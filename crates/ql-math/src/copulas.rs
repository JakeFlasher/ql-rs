use crate::{bivariate_normal_cdf, inv_norm_cdf};
use ql_ad::Real;

fn clamp_unit<T: Real>(x: T) -> T {
    x.max_val(T::from_f64(1e-14))
        .min_val(T::from_f64(1.0 - 1e-14))
}

#[derive(Copy, Clone, Debug)]
pub struct GaussianCopula<T: Real> {
    pub rho: T,
}
impl<T: Real> GaussianCopula<T> {
    pub fn cdf(&self, u: T, v: T) -> T {
        bivariate_normal_cdf(
            inv_norm_cdf(clamp_unit(u)),
            inv_norm_cdf(clamp_unit(v)),
            self.rho,
        )
    }
    pub fn conditional_cdf_v_given_u(&self, u: T, v: T) -> T {
        let x = inv_norm_cdf(clamp_unit(u));
        let y = inv_norm_cdf(clamp_unit(v));
        ((y - self.rho * x) / (T::one() - self.rho * self.rho).sqrt()).norm_cdf()
    }
    pub fn kendall_tau(&self) -> T {
        T::from_f64(2.0 / std::f64::consts::PI)
            * (self.rho / (T::one() - self.rho * self.rho).sqrt()).atan()
    }
    pub fn spearman_rho(&self) -> T {
        let x = self.rho * 0.5;
        T::from_f64(6.0 / std::f64::consts::PI) * (x / (T::one() - x * x).sqrt()).atan()
    }
}

#[derive(Copy, Clone, Debug)]
pub struct ClaytonCopula<T: Real> {
    pub theta: T,
}
impl<T: Real> ClaytonCopula<T> {
    pub fn cdf(&self, u: T, v: T) -> T {
        let u = clamp_unit(u);
        let v = clamp_unit(v);
        (u.pow(-self.theta) + v.pow(-self.theta) - 1.0).pow(-T::one() / self.theta)
    }
    pub fn conditional_cdf_v_given_u(&self, u: T, v: T) -> T {
        let u = clamp_unit(u);
        let v = clamp_unit(v);
        u.pow(-(self.theta + 1.0))
            * (u.pow(-self.theta) + v.pow(-self.theta) - 1.0)
                .pow(-(T::one() + T::one() / self.theta))
    }
    pub fn kendall_tau(&self) -> T {
        self.theta / (self.theta + 2.0)
    }
    pub fn lower_tail_dependence(&self) -> T {
        T::from_f64(2.0).pow(-T::one() / self.theta)
    }
}

#[derive(Copy, Clone, Debug)]
pub struct FrankCopula<T: Real> {
    pub theta: T,
}
impl<T: Real> FrankCopula<T> {
    pub fn cdf(&self, u: T, v: T) -> T {
        let u = clamp_unit(u);
        let v = clamp_unit(v);
        let et = (-self.theta).exp();
        let eu = (-self.theta * u).exp();
        let ev = (-self.theta * v).exp();
        (T::one() + (eu - 1.0) * (ev - 1.0) / (et - 1.0)).ln() / (-self.theta)
    }
    pub fn conditional_cdf_v_given_u(&self, u: T, v: T) -> T {
        let u = clamp_unit(u);
        let v = clamp_unit(v);
        let eu = (-self.theta * u).exp();
        let ev = (-self.theta * v).exp();
        eu * (ev - 1.0) / ((-self.theta).exp() - 1.0 + (eu - 1.0) * (ev - 1.0))
    }
}

#[derive(Copy, Clone, Debug)]
pub struct GumbelCopula<T: Real> {
    pub theta: T,
}
impl<T: Real> GumbelCopula<T> {
    pub fn cdf(&self, u: T, v: T) -> T {
        let u = clamp_unit(u);
        let v = clamp_unit(v);
        (-((-u.ln()).pow(self.theta) + (-v.ln()).pow(self.theta)).pow(T::one() / self.theta)).exp()
    }
    pub fn kendall_tau(&self) -> T {
        T::one() - T::one() / self.theta
    }
    pub fn upper_tail_dependence(&self) -> T {
        T::from_f64(2.0) - T::from_f64(2.0).pow(T::one() / self.theta)
    }
}

#[derive(Copy, Clone, Debug)]
pub struct FarlieGumbelMorgensternCopula<T: Real> {
    pub theta: T,
}
impl<T: Real> FarlieGumbelMorgensternCopula<T> {
    pub fn cdf(&self, u: T, v: T) -> T {
        u * v + self.theta * u * v * (T::one() - u) * (T::one() - v)
    }
}
