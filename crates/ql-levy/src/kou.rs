use crate::LevyModel;
use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct KouJumpDiffusion<T: Real> {
    pub sigma: T,
    pub lambda_j: T,
    pub p_up: T,
    pub eta1: T,
    pub eta2: T,
}
impl<T: Real> LevyModel<T> for KouJumpDiffusion<T> {
    fn martingale_compensator(&self) -> T {
        -self.sigma * self.sigma * 0.5
            - self.lambda_j
                * (self.p_up * self.eta1 / (self.eta1 - 1.0)
                    + (T::one() - self.p_up) * self.eta2 / (self.eta2 + 1.0)
                    - 1.0)
    }
    fn char_exponent(&self, u: T) -> (T, T) {
        let u2 = u * u;
        let r1 = self.eta1 * self.eta1 / (self.eta1 * self.eta1 + u2);
        let r2 = self.eta2 * self.eta2 / (self.eta2 * self.eta2 + u2);
        let re = -self.sigma * self.sigma * u2 * 0.5
            + self.lambda_j * (self.p_up * r1 + (T::one() - self.p_up) * r2 - 1.0);
        let im = self.lambda_j
            * u
            * (self.p_up * self.eta1 / (self.eta1 * self.eta1 + u2)
                - (T::one() - self.p_up) * self.eta2 / (self.eta2 * self.eta2 + u2));
        (re, im)
    }
    fn char_exponent_shifted(&self, u: T, beta: T) -> (T, T) {
        let u2 = u * u;
        let e1 = self.eta1 - beta;
        let e2 = self.eta2 + beta;
        let d1 = e1 * e1 + u2;
        let d2 = e2 * e2 + u2;
        let b2_u2 = beta * beta - u2;
        let re = self.sigma * self.sigma * b2_u2 * 0.5
            + self.lambda_j
                * (self.p_up * self.eta1 * e1 / d1 + (T::one() - self.p_up) * self.eta2 * e2 / d2
                    - 1.0);
        let im = self.sigma * self.sigma * u * beta
            + self.lambda_j
                * u
                * (self.p_up * self.eta1 / d1 - (T::one() - self.p_up) * self.eta2 / d2);
        (re, im)
    }
    fn cumulants(&self, r: T, q: T, t: T) -> (T, T, T) {
        let om = self.martingale_compensator();
        let c1 = (r - q
            + om
            + self.lambda_j * (self.p_up / self.eta1 - (T::one() - self.p_up) / self.eta2))
            * t;
        let c2 = (self.sigma * self.sigma
            + self.lambda_j
                * 2.0
                * (self.p_up / (self.eta1 * self.eta1)
                    + (T::one() - self.p_up) / (self.eta2 * self.eta2)))
            * t;
        let c4 = self.lambda_j
            * 24.0
            * (self.p_up / self.eta1.powf(4.0) + (T::one() - self.p_up) / self.eta2.powf(4.0))
            * t;
        (c1, c2, c4)
    }
}
