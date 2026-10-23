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
            // Taylor expansion of log(1-nu*theta-nu*sigma²/2)/nu;
            // retaining nu terms keeps the VG->Black limit differentiable.
            let a = self.theta + self.sigma * self.sigma * 0.5;
            -a - self.nu * a * a * 0.5 - self.nu * self.nu * a * a * a / 3.0
        } else {
            (T::one() - self.theta * self.nu - self.sigma * self.sigma * self.nu * 0.5).ln()
                / self.nu
        }
    }
    fn char_exponent(&self, u: T) -> (T, T) {
        if self.nu.val().abs() < 1e-10 {
            // Polar-form VG exponent expanded around nu=0.  This is the
            // Brownian limit plus its first two variance-gamma corrections.
            let cr = self.sigma * self.sigma * u * u * 0.5;
            let ci = -self.theta * u;
            let c2r = cr * cr - ci * ci;
            let c2i = cr * ci * 2.0;
            let c3r = c2r * cr - c2i * ci;
            let c3i = c2r * ci + c2i * cr;
            return (
                -cr + self.nu * c2r * 0.5 - self.nu * self.nu * c3r / 3.0,
                -ci + self.nu * c2i * 0.5 - self.nu * self.nu * c3i / 3.0,
            );
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
