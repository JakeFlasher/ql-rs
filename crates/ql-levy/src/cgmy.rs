use crate::LevyModel;
use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct CgmyProcess<T: Real> {
    pub sigma: T,
    pub c: T,
    pub g: T,
    pub m: T,
    pub y: T,
}
impl<T: Real> CgmyProcess<T> {
    fn gamma_neg_y(&self) -> T {
        let y0 = self.y.val();
        let scalar = |y: f64| libm::tgamma(2.0 - y) / (y * (y - 1.0));
        let base = scalar(y0);
        let h = 1e-5_f64.max(y0.abs() * 1e-5);
        let slope = (scalar(y0 + h) - scalar(y0 - h)) / (2.0 * h);
        T::from_f64(base) + (self.y - y0) * slope
    }
}
impl<T: Real> LevyModel<T> for CgmyProcess<T> {
    fn martingale_compensator(&self) -> T {
        -self.sigma * self.sigma * 0.5
            - self.c
                * self.gamma_neg_y()
                * ((self.m - 1.0).pow(self.y) - self.m.pow(self.y) + (self.g + 1.0).pow(self.y)
                    - self.g.pow(self.y))
    }
    fn char_exponent(&self, u: T) -> (T, T) {
        let rm = (self.m * self.m + u * u).sqrt();
        let rg = (self.g * self.g + u * u).sqrt();
        let am = (self.y * (u / self.m).atan()).cos();
        let sm = (self.y * (u / self.m).atan()).sin();
        let ag = (self.y * (u / self.g).atan()).cos();
        let sg = (self.y * (u / self.g).atan()).sin();
        let mm = rm.pow(self.y);
        let gg = rg.pow(self.y);
        let re = -self.sigma * self.sigma * u * u * 0.5
            + self.c
                * self.gamma_neg_y()
                * (mm * am - self.m.pow(self.y) + gg * ag - self.g.pow(self.y));
        let im = self.c * self.gamma_neg_y() * (-mm * sm + gg * sg);
        (re, im)
    }
    fn cumulants(&self, r: T, q: T, t: T) -> (T, T, T) {
        let om = self.martingale_compensator();
        let jump_mean = -self.c
            * self.gamma_neg_y()
            * self.y
            * (self.m.pow(self.y - 1.0) - self.g.pow(self.y - 1.0));
        let c1 = (r - q + om + jump_mean) * t;
        let c2 = (self.sigma * self.sigma
            + self.c
                * self.gamma_neg_y()
                * self.y
                * (self.y - 1.0)
                * (self.m.pow(self.y - 2.0) + self.g.pow(self.y - 2.0)))
            * t;
        let c4 = (self.c
            * self.gamma_neg_y()
            * self.y
            * (self.y - 1.0)
            * (self.y - 2.0)
            * (self.y - 3.0)
            * (self.m.pow(self.y - 4.0) + self.g.pow(self.y - 4.0)))
            * t;
        (c1, c2, c4)
    }
}
