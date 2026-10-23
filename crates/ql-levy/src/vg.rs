use crate::LevyModel;
use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct VarianceGamma<T: Real> {
    pub sigma: T,
    pub nu: T,
    pub theta: T,
}
impl<T: Real> LevyModel<T> for VarianceGamma<T> {
    fn martingale_compensator(&self) -> T {
        if self.nu.val().abs() < 1e-10 {
            // Brownian limit of the VG compensator; avoid a 0/0 when callers
            // intentionally use nu=0 as a degeneracy check.
            -self.theta - self.sigma * self.sigma * 0.5
        } else {
            (T::one() - self.theta * self.nu - self.sigma * self.sigma * self.nu * 0.5).ln()
                / self.nu
        }
    }
    fn char_exponent(&self, u: T) -> (T, T) {
        if self.nu.val().abs() < 1e-10 {
            return (-self.sigma * self.sigma * u * u * 0.5, self.theta * u);
        }
        let a = T::one() + self.sigma * self.sigma * self.nu * u * u * 0.5;
        let b = -self.theta * self.nu * u;
        let rr = (a * a + b * b).sqrt();
        let ph = (b / a).atan();
        (-rr.ln() / self.nu, -ph / self.nu)
    }
    fn cumulants(&self, r: T, q: T, t: T) -> (T, T, T) {
        let om = self.martingale_compensator();
        (
            (r - q + om + self.theta) * t,
            (self.sigma * self.sigma + self.nu * self.theta * self.theta) * t,
            (self.sigma.powf(4.0) * self.nu
                + self.theta.powf(4.0) * self.nu.powf(3.0) * 2.0
                + self.sigma * self.sigma * self.theta * self.theta * self.nu * self.nu * 4.0)
                * 3.0
                * t,
        )
    }
}
