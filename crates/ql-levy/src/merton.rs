use crate::LevyModel;
use ql_ad::Real;
use ql_instruments::OptionType;
use ql_vol::black_scholes_price;
#[derive(Clone, Copy, Debug)]
pub struct MertonJumpDiffusion<T: Real> {
    pub sigma: T,
    pub lambda_j: T,
    pub mu_j: T,
    pub delta_j: T,
}
impl<T: Real> LevyModel<T> for MertonJumpDiffusion<T> {
    fn martingale_compensator(&self) -> T {
        -self.sigma * self.sigma * 0.5
            - self.lambda_j * ((self.mu_j + self.delta_j * self.delta_j * 0.5).exp() - 1.0)
    }
    fn char_exponent(&self, u: T) -> (T, T) {
        let q = self.delta_j * self.delta_j;
        let e = (-q * u * u * 0.5).exp();
        let ang = u * self.mu_j;
        let re = -self.sigma * self.sigma * u * u * 0.5 + self.lambda_j * (e * ang.cos() - 1.0);
        let im = self.lambda_j * e * ang.sin();
        (re, im)
    }
    fn char_exponent_shifted(&self, u: T, beta: T) -> (T, T) {
        let q = self.delta_j * self.delta_j;
        let b2_u2 = beta * beta - u * u;
        let e = (beta * self.mu_j + q * b2_u2 * 0.5).exp();
        let ang = u * (self.mu_j + q * beta);
        let re = self.sigma * self.sigma * b2_u2 * 0.5 + self.lambda_j * (e * ang.cos() - 1.0);
        let im = self.sigma * self.sigma * u * beta + self.lambda_j * e * ang.sin();
        (re, im)
    }
    fn cumulants(&self, r: T, q: T, t: T) -> (T, T, T) {
        let om = self.martingale_compensator();
        let c1 = (r - q + om + self.lambda_j * self.mu_j) * t;
        let c2 = (self.sigma * self.sigma
            + self.lambda_j * (self.mu_j * self.mu_j + self.delta_j * self.delta_j))
            * t;
        let c4 = self.lambda_j
            * (self.mu_j.powf(4.0)
                + self.mu_j * self.mu_j * 6.0 * self.delta_j * self.delta_j
                + self.delta_j.powf(4.0) * 3.0)
            * t;
        (c1, c2, c4)
    }
}
impl<T: Real> MertonJumpDiffusion<T> {
    #[allow(clippy::too_many_arguments)]
    pub fn merton_series_price(
        &self,
        spot: T,
        strike: T,
        rate: T,
        div: T,
        expiry: T,
        kind: OptionType,
        max_terms: usize,
    ) -> T {
        let kappa = (self.mu_j + self.delta_j * self.delta_j * 0.5).exp() - 1.0;
        let lam = self.lambda_j;
        let mut w = (-lam * expiry).exp();
        let mut sum = T::zero();
        for n in 0..max_terms {
            let var = self.sigma * self.sigma
                + self.delta_j * self.delta_j * (n as f64) / expiry.max_val(T::from_f64(1e-12));
            let sig = var.sqrt();
            let rn = rate - self.lambda_j * kappa
                + (self.mu_j + self.delta_j * self.delta_j * 0.5) * (n as f64)
                    / expiry.max_val(T::from_f64(1e-12));
            let conditional =
                black_scholes_price(spot, strike, sig, expiry, rn, div, kind.is_call());
            // `black_scholes_price` discounts at r-lambda*kappa.  Convert it
            // to the risk-neutral discount r before applying the Poisson
            // probability.  The equivalent representation uses rn above.
            sum += w * ((rn - rate) * expiry).exp() * conditional;
            w *= lam * expiry / (n as f64 + 1.0);
        }
        sum
    }
}

/// Free-function spelling convenient for benchmark and scripting callers.
#[allow(clippy::too_many_arguments)]
pub fn merton_series_price<T: Real>(
    model: &MertonJumpDiffusion<T>,
    spot: T,
    strike: T,
    rate: T,
    div_yield: T,
    expiry: T,
    option_type: OptionType,
    max_terms: usize,
) -> T {
    model.merton_series_price(
        spot,
        strike,
        rate,
        div_yield,
        expiry,
        option_type,
        max_terms,
    )
}
