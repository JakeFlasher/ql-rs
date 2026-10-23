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
    /// Lanczos gamma approximation evaluated in the AD scalar type.  CGMY
    /// uses Gamma(2-Y)/(Y(Y-1)); keeping the approximation in `T` preserves
    /// curvature Greeks with respect to `Y` instead of lifting a finite
    /// difference constant from f64.
    fn gamma_lanczos(z: T) -> T {
        const C: [f64; 9] = [
            0.9999999999998099,
            676.5203681218851,
            -1259.1392167224028,
            771.3234287776531,
            -176.6150291621406,
            12.507343278686905,
            -0.13857109526572012,
            9.984369578019572e-6,
            1.5056327351493116e-7,
        ];
        // All production CGMY calibrations have 0 < Y < 2, hence 2-Y > 0;
        // the reflection branch is retained for a well-defined generic API.
        if z.val() < 0.5 {
            let pi = T::from_f64(std::f64::consts::PI);
            pi / ((pi * z).sin() * Self::gamma_lanczos(T::one() - z))
        } else {
            let zz = z - 1.0;
            let mut x = T::from_f64(C[0]);
            for (i, &ci) in C.iter().enumerate().skip(1) {
                x += T::from_f64(ci) / (zz + i as f64);
            }
            let t = zz + 7.5;
            T::from_f64((2.0 * std::f64::consts::PI).sqrt()) * t.pow(zz + 0.5) * (-t).exp() * x
        }
    }
    fn gamma_neg_y(&self) -> T {
        Self::gamma_lanczos(T::from_f64(2.0) - self.y) / (self.y * (self.y - 1.0))
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
