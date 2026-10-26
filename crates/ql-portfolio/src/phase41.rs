//! Engine 5A: Qian (2005) & Spinu (2013) Equal Risk Contribution (`ERC` / Risk Parity)
//! Convex Log-Barrier Solver & Diversification Ratio Engine.
//! Engine 5B: Rockafellar & Uryasev (2000) & Meucci (2005) Auxiliary-Threshold CVaR
//! Minimization & KKT Dual Tail-Probability Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct QianSpinuRiskParityLogBarrierSpec<T: Real = f64> {
    pub asset_volatilities: [T; 3],
    pub risk_budgets_b: [T; 3],
    pub corr_01: T,
    pub corr_02: T,
    pub corr_12: T,
}

#[derive(Debug, Clone, Copy)]
pub struct QianSpinuRiskParityLogBarrierReport<T: Real = f64> {
    pub optimal_erc_weights: [T; 3],
    pub risk_contribution_shares: [T; 3],
    pub portfolio_volatility: T,
    pub diversification_ratio: T,
    pub max_risk_budget_parity_error: T,
}

pub fn evaluate_qian_risk_parity_spinu_log_barrier<T: Real>(
    spec: &QianSpinuRiskParityLogBarrierSpec<T>,
) -> QianSpinuRiskParityLogBarrierReport<T> {
    let v0 = spec.asset_volatilities[0].max_val(T::from_f64(1e-4));
    let v1 = spec.asset_volatilities[1].max_val(T::from_f64(1e-4));
    let v2 = spec.asset_volatilities[2].max_val(T::from_f64(1e-4));
    let b_sum = (spec.risk_budgets_b[0] + spec.risk_budgets_b[1] + spec.risk_budgets_b[2])
        .max_val(T::from_f64(1e-6));
    let b = [
        spec.risk_budgets_b[0] / b_sum,
        spec.risk_budgets_b[1] / b_sum,
        spec.risk_budgets_b[2] / b_sum,
    ];

    let c00 = v0 * v0;
    let c11 = v1 * v1;
    let c22 = v2 * v2;
    let c01 = spec.corr_01 * v0 * v1;
    let c02 = spec.corr_02 * v0 * v2;
    let c12 = spec.corr_12 * v1 * v2;

    // Spinu (2013) coordinatewise exact quadratic updates on 0.5 y^T Sigma y - sum b_i ln y_i:
    // For asset i, c_{ii} y_i^2 + (sum_{j != i} c_{ij} y_j) y_i - b_i = 0:
    let mut y = [b[0].sqrt() / v0, b[1].sqrt() / v1, b[2].sqrt() / v2];
    for _ in 0..8 {
        let s0 = c01 * y[1] + c02 * y[2];
        y[0] = (-s0 + (s0 * s0 + T::from_f64(4.0) * c00 * b[0]).sqrt()) / (T::from_f64(2.0) * c00);
        let s1 = c01 * y[0] + c12 * y[2];
        y[1] = (-s1 + (s1 * s1 + T::from_f64(4.0) * c11 * b[1]).sqrt()) / (T::from_f64(2.0) * c11);
        let s2 = c02 * y[0] + c12 * y[1];
        y[2] = (-s2 + (s2 * s2 + T::from_f64(4.0) * c22 * b[2]).sqrt()) / (T::from_f64(2.0) * c22);
    }

    let y_sum = (y[0] + y[1] + y[2]).max_val(T::from_f64(1e-8));
    let w = [y[0] / y_sum, y[1] / y_sum, y[2] / y_sum];

    let mw0 = c00 * w[0] + c01 * w[1] + c02 * w[2];
    let mw1 = c01 * w[0] + c11 * w[1] + c12 * w[2];
    let mw2 = c02 * w[0] + c12 * w[1] + c22 * w[2];
    let port_var = (w[0] * mw0 + w[1] * mw1 + w[2] * mw2).max_val(T::from_f64(1e-10));
    let port_vol = port_var.sqrt();

    let rc = [
        (w[0] * mw0) / port_var,
        (w[1] * mw1) / port_var,
        (w[2] * mw2) / port_var,
    ];
    let err0 = (rc[0] - b[0]).abs();
    let err1 = (rc[1] - b[1]).abs();
    let err2 = (rc[2] - b[2]).abs();
    let max_err = err0.max_val(err1).max_val(err2);
    let weighted_vol = w[0] * v0 + w[1] * v1 + w[2] * v2;
    let dr = weighted_vol / port_vol;

    QianSpinuRiskParityLogBarrierReport {
        optimal_erc_weights: w,
        risk_contribution_shares: rc,
        portfolio_volatility: port_vol,
        diversification_ratio: dr,
        max_risk_budget_parity_error: max_err,
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RockafellarUryasevCvarAuxiliarySpec<T: Real = f64> {
    pub portfolio_weights: [T; 2],
    pub scenario_losses_asset0: [T; 4],
    pub scenario_losses_asset1: [T; 4],
    pub confidence_alpha: T,
    pub smoothing_epsilon: T,
}

#[derive(Debug, Clone, Copy)]
pub struct RockafellarUryasevCvarAuxiliaryReport<T: Real = f64> {
    pub optimal_var_threshold_zeta: T,
    pub minimized_auxiliary_cvar: T,
    pub expected_portfolio_loss: T,
    pub tail_excess_loss_integral: T,
    pub dual_tail_weight_sum: T,
}

pub fn evaluate_rockafellar_uryasev_cvar_minimization<T: Real>(
    spec: &RockafellarUryasevCvarAuxiliarySpec<T>,
) -> RockafellarUryasevCvarAuxiliaryReport<T> {
    let w0 = spec.portfolio_weights[0];
    let w1 = spec.portfolio_weights[1];
    let one_minus_alpha = (T::one() - spec.confidence_alpha).max_val(T::from_f64(0.05));
    let eps = spec.smoothing_epsilon.max_val(T::from_f64(1e-4));

    let mut losses = [T::zero(); 4];
    let mut mean_loss = T::zero();
    for s in 0..4 {
        losses[s] = w0 * spec.scenario_losses_asset0[s] + w1 * spec.scenario_losses_asset1[s];
        mean_loss += T::from_f64(0.25) * losses[s];
    }

    // Smooth Newton minimization of Rockafellar-Uryasev auxiliary function F_alpha(w, zeta):
    // F_alpha(w, zeta) = zeta + (1 / (1 - alpha)) * (1/4) sum_{s=1}^4 pos_{eff_eps}(L_s - zeta)
    // dF/dzeta = 1 - (1 / (1 - alpha)) * (1/4) sum_{s=1}^4 psi_{eff_eps}(L_s - zeta) = 0
    let span = (losses[3] - losses[2]).abs().max_val(T::from_f64(0.01));
    let eff_eps = eps + T::from_f64(0.25) * span;
    let mut zeta = T::from_f64(0.5) * (losses[2] + losses[3]);
    for _ in 0..8 {
        let mut prob_sum = T::zero();
        let mut dens_sum = T::zero();
        for &ls in &losses {
            let diff = ls - zeta;
            let r = (diff * diff + eff_eps * eff_eps).sqrt();
            let psi = T::from_f64(0.5) * (T::one() + diff / r);
            let dpsi = T::from_f64(0.5) * (eff_eps * eff_eps) / (r * r * r);
            prob_sum += T::from_f64(0.25) * psi;
            dens_sum += T::from_f64(0.25) * dpsi;
        }
        let grad = T::one() - prob_sum / one_minus_alpha;
        let hess = (dens_sum / one_minus_alpha).max_val(T::from_f64(1e-4));
        zeta -= grad / hess;
    }

    let mut excess_sum = T::zero();
    let mut dual_sum = T::zero();
    for &ls in &losses {
        let diff = ls - zeta;
        let r = (diff * diff + eff_eps * eff_eps).sqrt();
        let pos = T::from_f64(0.5) * (diff + r);
        excess_sum += T::from_f64(0.25) * pos;
        let psi = T::from_f64(0.5) * (T::one() + diff / r);
        dual_sum += T::from_f64(0.25) * psi / one_minus_alpha;
    }
    let cvar = zeta + excess_sum / one_minus_alpha;

    RockafellarUryasevCvarAuxiliaryReport {
        optimal_var_threshold_zeta: zeta,
        minimized_auxiliary_cvar: cvar,
        expected_portfolio_loss: mean_loss,
        tail_excess_loss_integral: excess_sum,
        dual_tail_weight_sum: dual_sum,
    }
}

#[cfg(test)]
mod phase41_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn qian_spinu_erc_and_rockafellar_uryasev_cvar_invariants_and_ad_parity() {
        let s1 = QianSpinuRiskParityLogBarrierSpec {
            asset_volatilities: [0.18, 0.12, 0.08],
            risk_budgets_b: [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            corr_01: 0.25,
            corr_02: -0.10,
            corr_12: 0.15,
        };
        let r1 = evaluate_qian_risk_parity_spinu_log_barrier(&s1);
        assert!(r1.max_risk_budget_parity_error < 1e-10);
        assert!(r1.diversification_ratio > 1.0);
        assert!(r1.optimal_erc_weights[2] > r1.optimal_erc_weights[0]);

        let h = 1e-6;
        let r1_ad =
            evaluate_qian_risk_parity_spinu_log_barrier(&QianSpinuRiskParityLogBarrierSpec {
                asset_volatilities: [
                    Dual::<1>::variable(s1.asset_volatilities[0], 0),
                    Dual::<1>::from_value(s1.asset_volatilities[1]),
                    Dual::<1>::from_value(s1.asset_volatilities[2]),
                ],
                risk_budgets_b: s1.risk_budgets_b.map(Dual::<1>::from_value),
                corr_01: Dual::<1>::from_value(s1.corr_01),
                corr_02: Dual::<1>::from_value(s1.corr_02),
                corr_12: Dual::<1>::from_value(s1.corr_12),
            });
        let mut s1p = s1;
        s1p.asset_volatilities[0] += h;
        let mut s1m = s1;
        s1m.asset_volatilities[0] -= h;
        let fd1 = (evaluate_qian_risk_parity_spinu_log_barrier(&s1p).portfolio_volatility
            - evaluate_qian_risk_parity_spinu_log_barrier(&s1m).portfolio_volatility)
            / (2.0 * h);
        assert!((r1_ad.portfolio_volatility.d[0] - fd1).abs() < 1e-6);

        let s2 = RockafellarUryasevCvarAuxiliarySpec {
            portfolio_weights: [0.6, 0.4],
            scenario_losses_asset0: [-0.02, 0.01, 0.04, 0.09],
            scenario_losses_asset1: [-0.01, 0.005, 0.02, 0.05],
            confidence_alpha: 0.75,
            smoothing_epsilon: 0.002,
        };
        let r2 = evaluate_rockafellar_uryasev_cvar_minimization(&s2);
        assert!(r2.minimized_auxiliary_cvar > r2.expected_portfolio_loss);
        assert!((r2.dual_tail_weight_sum - 1.0).abs() < 1e-5);

        let r2_ad =
            evaluate_rockafellar_uryasev_cvar_minimization(&RockafellarUryasevCvarAuxiliarySpec {
                portfolio_weights: [
                    Dual::<1>::variable(s2.portfolio_weights[0], 0),
                    Dual::<1>::from_value(s2.portfolio_weights[1]),
                ],
                scenario_losses_asset0: s2.scenario_losses_asset0.map(Dual::<1>::from_value),
                scenario_losses_asset1: s2.scenario_losses_asset1.map(Dual::<1>::from_value),
                confidence_alpha: Dual::<1>::from_value(s2.confidence_alpha),
                smoothing_epsilon: Dual::<1>::from_value(s2.smoothing_epsilon),
            });
        let mut s2p = s2;
        s2p.portfolio_weights[0] += h;
        let mut s2m = s2;
        s2m.portfolio_weights[0] -= h;
        let fd2 = (evaluate_rockafellar_uryasev_cvar_minimization(&s2p).minimized_auxiliary_cvar
            - evaluate_rockafellar_uryasev_cvar_minimization(&s2m).minimized_auxiliary_cvar)
            / (2.0 * h);
        assert!((r2_ad.minimized_auxiliary_cvar.d[0] - fd2).abs() < 1e-5);
    }
}
