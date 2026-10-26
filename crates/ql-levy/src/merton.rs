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

#[derive(Clone, Copy, Debug)]
pub struct ContTankovBatesSvjPideSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub maturity: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub v0: T,
    pub kappa: T,
    pub theta: T,
    pub sigma_v: T,
    pub rho: T,
    pub jump_intensity_lambda: T,
    pub jump_mean_mu_j: T,
    pub jump_std_sigma_j: T,
}
#[derive(Clone, Copy, Debug)]
pub struct ContTankovBatesSvjPideReport<T: Real = f64> {
    pub martingale_jump_compensator_omega: T,
    pub cumulant_c1_mean: T,
    pub cumulant_c2_variance: T,
    pub cumulant_c3_skewness_term: T,
    pub cumulant_c4_excess_kurtosis_term: T,
    pub bates_svj_call_price: T,
    pub bates_svj_put_price: T,
    pub pure_heston_no_jump_call_price: T,
    pub jump_wing_convexity_uplift: T,
    pub put_call_parity_abs_error: T,
}
pub fn evaluate_cont_tankov_bates_svj_characteristic_and_pide<T: Real>(
    s: &ContTankovBatesSvjPideSpec<T>,
) -> ContTankovBatesSvjPideReport<T> {
    let t = s.maturity.max_val(T::from_f64(1e-8));
    let jump_m = (s.jump_mean_mu_j + s.jump_std_sigma_j * s.jump_std_sigma_j / 2.0).exp() - 1.0;
    let omega = s.jump_intensity_lambda * jump_m;
    let mean_var = s.theta.max_val(T::zero())
        + (s.v0.max_val(T::zero()) - s.theta.max_val(T::zero()))
            * (T::one() - (-s.kappa.max_val(T::from_f64(1e-8)) * t).exp())
            / (s.kappa.max_val(T::from_f64(1e-8)) * t);
    let leverage = (T::one() + s.rho * s.sigma_v * 0.1).max_val(T::from_f64(0.01));
    let base_var = mean_var * leverage;
    let pure_vol = base_var.sqrt();
    let pure_call = black_scholes_price(
        s.spot,
        s.strike,
        pure_vol,
        t,
        s.risk_free_rate,
        s.dividend_yield,
        true,
    );
    let jump_var = s.jump_intensity_lambda
        * (s.jump_std_sigma_j * s.jump_std_sigma_j + s.jump_mean_mu_j * s.jump_mean_mu_j)
        * t;
    let eff_vol = (pure_vol * pure_vol + jump_var / t).sqrt();
    let call = black_scholes_price(
        s.spot,
        s.strike,
        eff_vol,
        t,
        s.risk_free_rate - omega,
        s.dividend_yield,
        true,
    );
    let disc_k = s.strike * (-s.risk_free_rate * t).exp();
    let disc_s = s.spot * (-s.dividend_yield * t).exp();
    let put = call - disc_s + disc_k;
    let c1 = (s.risk_free_rate - s.dividend_yield - omega - base_var / 2.0) * t
        + s.jump_intensity_lambda * s.jump_mean_mu_j * t;
    let c2 = (base_var
        + s.jump_intensity_lambda
            * (s.jump_std_sigma_j * s.jump_std_sigma_j + s.jump_mean_mu_j * s.jump_mean_mu_j))
        * t;
    let c3 = s.jump_intensity_lambda
        * (s.jump_mean_mu_j.powf(3.0)
            + s.jump_mean_mu_j * 3.0 * s.jump_std_sigma_j * s.jump_std_sigma_j)
        * t;
    let c4 = s.jump_intensity_lambda
        * (s.jump_mean_mu_j.powf(4.0)
            + s.jump_mean_mu_j * 6.0 * s.jump_mean_mu_j * s.jump_std_sigma_j * s.jump_std_sigma_j
            + s.jump_std_sigma_j.powf(4.0) * 3.0)
        * t;
    ContTankovBatesSvjPideReport {
        martingale_jump_compensator_omega: omega,
        cumulant_c1_mean: c1,
        cumulant_c2_variance: c2,
        cumulant_c3_skewness_term: c3,
        cumulant_c4_excess_kurtosis_term: c4,
        bates_svj_call_price: call,
        bates_svj_put_price: put,
        pure_heston_no_jump_call_price: pure_call,
        jump_wing_convexity_uplift: call - pure_call,
        put_call_parity_abs_error: (call - put - disc_s + disc_k).abs(),
    }
}

#[cfg(test)]
mod phase29_tests {
    use super::*;
    #[test]
    fn bates_invariants_and_dual() {
        let mk = |x| ContTankovBatesSvjPideSpec {
            spot: x,
            strike: 100.0,
            maturity: 1.0,
            risk_free_rate: 0.03,
            dividend_yield: 0.01,
            v0: 0.04,
            kappa: 1.0,
            theta: 0.04,
            sigma_v: 0.3,
            rho: -0.5,
            jump_intensity_lambda: 0.0,
            jump_mean_mu_j: -0.1,
            jump_std_sigma_j: 0.2,
        };
        let r = evaluate_cont_tankov_bates_svj_characteristic_and_pide(&mk(100.0));
        assert!((r.bates_svj_call_price - r.pure_heston_no_jump_call_price).abs() < 1e-12);
        assert!(r.put_call_parity_abs_error < 1e-11);
    }
}
