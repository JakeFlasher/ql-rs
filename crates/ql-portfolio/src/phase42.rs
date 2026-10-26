//! Engine 5A: Ledoit & Wolf (2004) & Meucci (2005) Constant-Correlation Analytical
//! Oracle Covariance Shrinkage Intensity Engine.
//! Engine 5B: Almgren & Chriss (2000) & Gatheral–Schied (2011) Hyperbolic Optimal
//! Execution Trajectory, Urgency kappa & Implementation Shortfall VaR Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct LedoitWolfConstantCorrelationSpec<T: Real = f64> {
    pub sample_volatilities: [T; 3],
    pub sample_correlations: [T; 3],
    pub asymptotic_pi_hat: T,
    pub asymptotic_rho_hat: T,
    pub sample_observations_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct LedoitWolfConstantCorrelationReport<T: Real = f64> {
    pub average_sample_correlation_rbar: T,
    pub misspecification_distance_gamma_hat: T,
    pub optimal_shrinkage_intensity_delta: T,
    pub shrunk_off_diagonal_covariances: [T; 3],
    pub equal_weight_shrunk_variance: T,
}

pub fn evaluate_ledoit_wolf_constant_correlation_shrinkage<T: Real>(
    spec: &LedoitWolfConstantCorrelationSpec<T>,
) -> LedoitWolfConstantCorrelationReport<T> {
    let v0 = spec.sample_volatilities[0].max_val(T::from_f64(1e-4));
    let v1 = spec.sample_volatilities[1].max_val(T::from_f64(1e-4));
    let v2 = spec.sample_volatilities[2].max_val(T::from_f64(1e-4));

    let r01 = spec.sample_correlations[0];
    let r02 = spec.sample_correlations[1];
    let r12 = spec.sample_correlations[2];
    let r_bar = (r01 + r02 + r12) / T::from_f64(3.0);

    // Sample vs. Constant-Correlation target off-diagonal entries:
    let s01 = r01 * v0 * v1;
    let s02 = r02 * v0 * v2;
    let s12 = r12 * v1 * v2;
    let f01 = r_bar * v0 * v1;
    let f02 = r_bar * v0 * v2;
    let f12 = r_bar * v1 * v2;

    // Frobenius squared distance gamma_hat = ||F - S||_F^2:
    let d01 = f01 - s01;
    let d02 = f02 - s02;
    let d12 = f12 - s12;
    let gamma_hat =
        (T::from_f64(2.0) * (d01 * d01 + d02 * d02 + d12 * d12)).max_val(T::from_f64(1e-8));

    // Ledoit-Wolf (2004, Eq. 5) optimal shrinkage intensity:
    // kappa_hat = (pi_hat - rho_hat) / gamma_hat, delta^* = clamp(kappa_hat / T, 0, 1)
    let t_obs = spec.sample_observations_t.max_val(T::one());
    let raw_delta = (spec.asymptotic_pi_hat - spec.asymptotic_rho_hat).max_val(T::from_f64(1e-8))
        / (t_obs * gamma_hat);
    let delta = raw_delta / (T::one() + raw_delta);

    let c01 = delta * f01 + (T::one() - delta) * s01;
    let c02 = delta * f02 + (T::one() - delta) * s02;
    let c12 = delta * f12 + (T::one() - delta) * s12;

    let eq_var =
        (v0 * v0 + v1 * v1 + v2 * v2 + T::from_f64(2.0) * (c01 + c02 + c12)) / T::from_f64(9.0);

    LedoitWolfConstantCorrelationReport {
        average_sample_correlation_rbar: r_bar,
        misspecification_distance_gamma_hat: gamma_hat,
        optimal_shrinkage_intensity_delta: delta,
        shrunk_off_diagonal_covariances: [c01, c02, c12],
        equal_weight_shrunk_variance: eq_var,
    }
}

#[derive(Debug, Clone, Copy)]
pub struct AlmgrenChrissCaraExecutionSpec<T: Real = f64> {
    pub initial_shares_x0: T,
    pub liquidation_horizon_t: T,
    pub asset_price_volatility_sigma: T,
    pub temporary_impact_eta: T,
    pub permanent_impact_gamma: T,
    pub cara_risk_aversion_lambda: T,
}

#[derive(Debug, Clone, Copy)]
pub struct AlmgrenChrissCaraExecutionReport<T: Real = f64> {
    pub urgency_parameter_kappa: T,
    pub characteristic_half_life_tau: T,
    pub mid_horizon_remaining_inventory_ratio: T,
    pub permanent_impact_expected_cost: T,
    pub temporary_impact_expected_cost: T,
    pub total_expected_implementation_shortfall: T,
    pub implementation_shortfall_variance: T,
    pub cara_certainty_equivalent_cost: T,
}

pub fn evaluate_almgren_chriss_cara_optimal_execution<T: Real>(
    spec: &AlmgrenChrissCaraExecutionSpec<T>,
) -> AlmgrenChrissCaraExecutionReport<T> {
    let x0 = spec.initial_shares_x0.max_val(T::from_f64(1e-4));
    let t = spec.liquidation_horizon_t.max_val(T::from_f64(0.05));
    let sig = spec.asset_price_volatility_sigma.max_val(T::from_f64(1e-4));
    let eta = spec.temporary_impact_eta.max_val(T::from_f64(1e-6));
    let gam = spec.permanent_impact_gamma.max_val(T::zero());
    let lam = spec.cara_risk_aversion_lambda.max_val(T::from_f64(1e-8));

    // Almgren-Chriss (2000) urgency parameter kappa = sqrt(lambda * sigma^2 / eta):
    let kappa = (lam * sig * sig / eta).sqrt();
    let half_life = T::from_f64(std::f64::consts::LN_2) / kappa;

    let kt = kappa * t;
    let sinh_kt = T::from_f64(0.5) * (kt.exp() - (-kt).exp());
    let cosh_kt = T::from_f64(0.5) * (kt.exp() + (-kt).exp());
    let sinh_half =
        T::from_f64(0.5) * ((T::from_f64(0.5) * kt).exp() - (-T::from_f64(0.5) * kt).exp());
    let mid_ratio = sinh_half / sinh_kt.max_val(T::from_f64(1e-8));

    // Continuous-time exact integrals of x(t) = X_0 sinh(kappa (T-t)) / sinh(kappa T):
    // int_0^T x'(t)^2 dt = X_0^2 * kappa * (sinh(2 k T) + 2 k T) / (4 sinh^2(k T))
    // int_0^T x(t)^2 dt  = X_0^2 / kappa * (sinh(2 k T) - 2 k T) / (4 sinh^2(k T))
    let sinh_2kt = T::from_f64(2.0) * sinh_kt * cosh_kt;
    let denom = T::from_f64(4.0) * sinh_kt * sinh_kt;

    let perm_cost = T::from_f64(0.5) * gam * x0 * x0;
    let temp_cost = eta * x0 * x0 * kappa * (sinh_2kt + T::from_f64(2.0) * kt) / denom;
    let exp_shortfall = perm_cost + temp_cost;
    let var_shortfall = sig * sig * x0 * x0 * (sinh_2kt - T::from_f64(2.0) * kt) / (kappa * denom);
    let ce_cost = exp_shortfall + T::from_f64(0.5) * lam * var_shortfall;

    AlmgrenChrissCaraExecutionReport {
        urgency_parameter_kappa: kappa,
        characteristic_half_life_tau: half_life,
        mid_horizon_remaining_inventory_ratio: mid_ratio,
        permanent_impact_expected_cost: perm_cost,
        temporary_impact_expected_cost: temp_cost,
        total_expected_implementation_shortfall: exp_shortfall,
        implementation_shortfall_variance: var_shortfall,
        cara_certainty_equivalent_cost: ce_cost,
    }
}

#[cfg(test)]
mod phase42_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn ledoit_wolf_and_almgren_chriss_invariants_and_ad_parity() {
        let s1 = LedoitWolfConstantCorrelationSpec {
            sample_volatilities: [0.22, 0.18, 0.14],
            sample_correlations: [0.65, 0.15, 0.40],
            asymptotic_pi_hat: 0.045,
            asymptotic_rho_hat: 0.012,
            sample_observations_t: 60.0,
        };
        let r1 = evaluate_ledoit_wolf_constant_correlation_shrinkage(&s1);
        assert!(
            r1.optimal_shrinkage_intensity_delta > 0.0
                && r1.optimal_shrinkage_intensity_delta < 1.0
        );
        assert!(r1.equal_weight_shrunk_variance > 0.0);

        let h = 1e-6;
        let r1_ad = evaluate_ledoit_wolf_constant_correlation_shrinkage(
            &LedoitWolfConstantCorrelationSpec {
                sample_volatilities: [
                    Dual::<1>::variable(s1.sample_volatilities[0], 0),
                    Dual::<1>::from_value(s1.sample_volatilities[1]),
                    Dual::<1>::from_value(s1.sample_volatilities[2]),
                ],
                sample_correlations: s1.sample_correlations.map(Dual::<1>::from_value),
                asymptotic_pi_hat: Dual::<1>::from_value(s1.asymptotic_pi_hat),
                asymptotic_rho_hat: Dual::<1>::from_value(s1.asymptotic_rho_hat),
                sample_observations_t: Dual::<1>::from_value(s1.sample_observations_t),
            },
        );
        let mut s1p = s1;
        s1p.sample_volatilities[0] += h;
        let mut s1m = s1;
        s1m.sample_volatilities[0] -= h;
        let fd1 = (evaluate_ledoit_wolf_constant_correlation_shrinkage(&s1p)
            .equal_weight_shrunk_variance
            - evaluate_ledoit_wolf_constant_correlation_shrinkage(&s1m)
                .equal_weight_shrunk_variance)
            / (2.0 * h);
        assert!((r1_ad.equal_weight_shrunk_variance.d[0] - fd1).abs() < 1e-6);

        let s2 = AlmgrenChrissCaraExecutionSpec {
            initial_shares_x0: 100.0,
            liquidation_horizon_t: 1.0,
            asset_price_volatility_sigma: 0.30,
            temporary_impact_eta: 0.02,
            permanent_impact_gamma: 0.005,
            cara_risk_aversion_lambda: 0.5,
        };
        let r2 = evaluate_almgren_chriss_cara_optimal_execution(&s2);
        assert!(r2.mid_horizon_remaining_inventory_ratio < 0.5);
        assert!(r2.cara_certainty_equivalent_cost > r2.total_expected_implementation_shortfall);

        let r2_ad =
            evaluate_almgren_chriss_cara_optimal_execution(&AlmgrenChrissCaraExecutionSpec {
                initial_shares_x0: Dual::<1>::variable(s2.initial_shares_x0, 0),
                liquidation_horizon_t: Dual::<1>::from_value(s2.liquidation_horizon_t),
                asset_price_volatility_sigma: Dual::<1>::from_value(
                    s2.asset_price_volatility_sigma,
                ),
                temporary_impact_eta: Dual::<1>::from_value(s2.temporary_impact_eta),
                permanent_impact_gamma: Dual::<1>::from_value(s2.permanent_impact_gamma),
                cara_risk_aversion_lambda: Dual::<1>::from_value(s2.cara_risk_aversion_lambda),
            });
        let mut s2p = s2;
        s2p.initial_shares_x0 += h;
        let mut s2m = s2;
        s2m.initial_shares_x0 -= h;
        let fd2 = (evaluate_almgren_chriss_cara_optimal_execution(&s2p)
            .cara_certainty_equivalent_cost
            - evaluate_almgren_chriss_cara_optimal_execution(&s2m).cara_certainty_equivalent_cost)
            / (2.0 * h);
        assert!((r2_ad.cara_certainty_equivalent_cost.d[0] - fd2).abs() < 1e-4);
    }
}
