//! Geman (2005) & Eydeland–Wolyniec (2003) Gibson–Schwartz (1990) Two-Factor Spot &
//! Stochastic Convenience-Yield Commodity Futures & Samuelson Volatility Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct GemanGibsonSchwartzSpec<T: Real = f64> {
    pub spot_price_s0: T,
    pub convenience_yield_delta0: T,
    pub long_run_convenience_alpha: T,
    pub mean_reversion_kappa: T,
    pub spot_volatility_sigma1: T,
    pub convenience_volatility_sigma2: T,
    pub spot_convenience_corr_rho: T,
    pub risk_free_rate: T,
    pub maturity_years_tau: T,
}

#[derive(Debug, Clone, Copy)]
pub struct GemanGibsonSchwartzReport<T: Real = f64> {
    pub convenience_loading_b_tau: T,
    pub affine_convexity_a_tau: T,
    pub gibson_schwartz_futures_price: T,
    pub instantaneous_futures_vol: T,
    pub samuelson_vol_attenuation_ratio: T,
    pub annualized_roll_yield: T,
}

pub fn evaluate_geman_gibson_schwartz_commodity_futures<T: Real>(
    spec: &GemanGibsonSchwartzSpec<T>,
) -> GemanGibsonSchwartzReport<T> {
    let tau = spec.maturity_years_tau.max_val(T::from_f64(0.05));
    let k = spec.mean_reversion_kappa.max_val(T::from_f64(1e-3));
    let k2 = k * k;
    let k3 = k2 * k;
    let s1 = spec.spot_volatility_sigma1;
    let s2 = spec.convenience_volatility_sigma2;
    let rho = spec.spot_convenience_corr_rho;

    // B(tau) = (1 - exp(-kappa * tau)) / kappa
    let exp_kt = (-k * tau).exp();
    let b_tau = (T::one() - exp_kt) / k;

    // Exact Gibson-Schwartz (1990) affine convexity term A(tau):
    let term1 = (spec.risk_free_rate - spec.long_run_convenience_alpha
        + (s2 * s2) / (T::from_f64(2.0) * k2)
        - (rho * s1 * s2) / k)
        * tau;
    let term2 =
        (s2 * s2) * (T::one() - (-T::from_f64(2.0) * k * tau).exp()) / (T::from_f64(4.0) * k3);
    let term3 = (spec.long_run_convenience_alpha * k + rho * s1 * s2 - (s2 * s2) / k)
        * (T::one() - exp_kt)
        / k2;
    let a_tau = term1 + term2 + term3;

    let fut_price = spec.spot_price_s0 * (-spec.convenience_yield_delta0 * b_tau + a_tau).exp();

    // Instantaneous futures volatility sigma_F(tau) = sqrt(sigma_1^2 + sigma_2^2 B(tau)^2 - 2 rho sigma_1 sigma_2 B(tau))
    let var_f = (s1 * s1 + s2 * s2 * b_tau * b_tau - T::from_f64(2.0) * rho * s1 * s2 * b_tau)
        .max_val(T::from_f64(1e-8));
    let sig_f = var_f.sqrt();
    let samuelson_ratio = sig_f / s1.max_val(T::from_f64(1e-6));
    let roll_yield = (spec.spot_price_s0 / fut_price).ln() / tau;

    GemanGibsonSchwartzReport {
        convenience_loading_b_tau: b_tau,
        affine_convexity_a_tau: a_tau,
        gibson_schwartz_futures_price: fut_price,
        instantaneous_futures_vol: sig_f,
        samuelson_vol_attenuation_ratio: samuelson_ratio,
        annualized_roll_yield: roll_yield,
    }
}

#[cfg(test)]
mod phase39_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn geman_gibson_schwartz_invariants_and_ad_parity() {
        let spec = GemanGibsonSchwartzSpec {
            spot_price_s0: 82.0,
            convenience_yield_delta0: 0.14,
            long_run_convenience_alpha: 0.08,
            mean_reversion_kappa: 1.4,
            spot_volatility_sigma1: 0.34,
            convenience_volatility_sigma2: 0.26,
            spot_convenience_corr_rho: 0.65,
            risk_free_rate: 0.04,
            maturity_years_tau: 1.5,
        };
        let rep = evaluate_geman_gibson_schwartz_commodity_futures(&spec);
        // High positive convenience yield induces backwardation (positive roll yield):
        assert!(rep.annualized_roll_yield > 0.0);
        // Positive correlation between spot and convenience yield dampens far-dated futures vol (Samuelson effect):
        assert!(rep.samuelson_vol_attenuation_ratio < 1.0);

        let h = 1e-6;
        let r_ad = evaluate_geman_gibson_schwartz_commodity_futures(&GemanGibsonSchwartzSpec {
            spot_price_s0: Dual::<1>::from_value(spec.spot_price_s0),
            convenience_yield_delta0: Dual::<1>::variable(spec.convenience_yield_delta0, 0),
            long_run_convenience_alpha: Dual::<1>::from_value(spec.long_run_convenience_alpha),
            mean_reversion_kappa: Dual::<1>::from_value(spec.mean_reversion_kappa),
            spot_volatility_sigma1: Dual::<1>::from_value(spec.spot_volatility_sigma1),
            convenience_volatility_sigma2: Dual::<1>::from_value(
                spec.convenience_volatility_sigma2,
            ),
            spot_convenience_corr_rho: Dual::<1>::from_value(spec.spot_convenience_corr_rho),
            risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
            maturity_years_tau: Dual::<1>::from_value(spec.maturity_years_tau),
        });
        let mut sp = spec;
        sp.convenience_yield_delta0 += h;
        let mut sm = spec;
        sm.convenience_yield_delta0 -= h;
        let fd = (evaluate_geman_gibson_schwartz_commodity_futures(&sp)
            .gibson_schwartz_futures_price
            - evaluate_geman_gibson_schwartz_commodity_futures(&sm).gibson_schwartz_futures_price)
            / (2.0 * h);
        assert!((r_ad.gibson_schwartz_futures_price.d[0] - fd).abs() < 1e-5);
    }
}
