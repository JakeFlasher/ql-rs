//! Adrian, Crump & Moench (2013, `ACM`) & Cochrane–Piazzesi (2005) Tent-Factor
//! Return-Forecasting & 10Y Affine Term Premium Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct AcmCochranePiazzesiTermPremiumSpec<T: Real = f64> {
    pub forward_rates_1y_to_5y: [T; 5],
    pub tent_weights_gamma: [T; 5],
    pub long_run_neutral_short_rate_r_star: T,
    pub physical_var_mean_reversion_phi: T,
    pub market_price_of_level_risk_lambda0: T,
    pub market_price_of_slope_risk_lambda1: T,
    pub horizon_years_n: T,
}

#[derive(Debug, Clone, Copy)]
pub struct AcmCochranePiazzesiTermPremiumReport<T: Real = f64> {
    pub cochrane_piazzesi_tent_factor_cp: T,
    pub expectations_hypothesis_10y_yield: T,
    pub acm_10y_term_premium: T,
    pub fitted_10y_nominal_yield: T,
    pub expected_1y_excess_bond_return: T,
}

pub fn evaluate_acm_cochrane_piazzesi_affine_term_premium<T: Real>(
    spec: &AcmCochranePiazzesiTermPremiumSpec<T>,
) -> AcmCochranePiazzesiTermPremiumReport<T> {
    let f = spec.forward_rates_1y_to_5y;
    let g = spec.tent_weights_gamma;
    let n = spec.horizon_years_n.max_val(T::one());
    let phi = spec
        .physical_var_mean_reversion_phi
        .max_val(T::from_f64(0.05));

    // Cochrane-Piazzesi (2005) tent-shaped linear combination of 1Y..5Y forward rates:
    let mut cp_factor = T::zero();
    for i in 0..5 {
        cp_factor += g[i] * f[i];
    }

    // Physical P-measure Expectations Hypothesis average short rate over [0, N]:
    let r0 = f[0];
    let r_star = spec.long_run_neutral_short_rate_r_star;
    let p_decay = (T::one() - (-phi * n).exp()) / (phi * n);
    let eh_yield = r_star + (r0 - r_star) * p_decay;

    // Slope factor (5Y forward minus 1Y spot) + CP tent factor drive the ACM term premium:
    let slope = f[4] - f[0];
    let term_premium = spec.market_price_of_level_risk_lambda0 * cp_factor
        + spec.market_price_of_slope_risk_lambda1 * slope * (T::one() - p_decay);

    let y10 = eh_yield + term_premium;
    let exp_excess_ret = cp_factor + T::from_f64(0.5) * term_premium * n;

    AcmCochranePiazzesiTermPremiumReport {
        cochrane_piazzesi_tent_factor_cp: cp_factor,
        expectations_hypothesis_10y_yield: eh_yield,
        acm_10y_term_premium: term_premium,
        fitted_10y_nominal_yield: y10,
        expected_1y_excess_bond_return: exp_excess_ret,
    }
}

#[cfg(test)]
mod phase42_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn acm_cochrane_piazzesi_affine_term_premium_invariants_and_ad_parity() {
        let spec = AcmCochranePiazzesiTermPremiumSpec {
            forward_rates_1y_to_5y: [0.042, 0.045, 0.049, 0.047, 0.046],
            tent_weights_gamma: [-1.2, 0.4, 2.1, 0.5, -1.5],
            long_run_neutral_short_rate_r_star: 0.035,
            physical_var_mean_reversion_phi: 0.22,
            market_price_of_level_risk_lambda0: 0.35,
            market_price_of_slope_risk_lambda1: 0.80,
            horizon_years_n: 10.0,
        };
        let rep = evaluate_acm_cochrane_piazzesi_affine_term_premium(&spec);
        assert!(rep.cochrane_piazzesi_tent_factor_cp > 0.0);
        assert!(rep.acm_10y_term_premium > 0.0);
        assert!(rep.fitted_10y_nominal_yield > rep.expectations_hypothesis_10y_yield);

        let h = 1e-6;
        let r_ad = evaluate_acm_cochrane_piazzesi_affine_term_premium(
            &AcmCochranePiazzesiTermPremiumSpec {
                forward_rates_1y_to_5y: [
                    Dual::<1>::from_value(spec.forward_rates_1y_to_5y[0]),
                    Dual::<1>::from_value(spec.forward_rates_1y_to_5y[1]),
                    Dual::<1>::variable(spec.forward_rates_1y_to_5y[2], 0),
                    Dual::<1>::from_value(spec.forward_rates_1y_to_5y[3]),
                    Dual::<1>::from_value(spec.forward_rates_1y_to_5y[4]),
                ],
                tent_weights_gamma: spec.tent_weights_gamma.map(Dual::<1>::from_value),
                long_run_neutral_short_rate_r_star: Dual::<1>::from_value(
                    spec.long_run_neutral_short_rate_r_star,
                ),
                physical_var_mean_reversion_phi: Dual::<1>::from_value(
                    spec.physical_var_mean_reversion_phi,
                ),
                market_price_of_level_risk_lambda0: Dual::<1>::from_value(
                    spec.market_price_of_level_risk_lambda0,
                ),
                market_price_of_slope_risk_lambda1: Dual::<1>::from_value(
                    spec.market_price_of_slope_risk_lambda1,
                ),
                horizon_years_n: Dual::<1>::from_value(spec.horizon_years_n),
            },
        );
        let mut sp = spec;
        sp.forward_rates_1y_to_5y[2] += h;
        let mut sm = spec;
        sm.forward_rates_1y_to_5y[2] -= h;
        let fd = (evaluate_acm_cochrane_piazzesi_affine_term_premium(&sp).fitted_10y_nominal_yield
            - evaluate_acm_cochrane_piazzesi_affine_term_premium(&sm).fitted_10y_nominal_yield)
            / (2.0 * h);
        assert!((r_ad.fitted_10y_nominal_yield.d[0] - fd).abs() < 1e-6);
    }
}
