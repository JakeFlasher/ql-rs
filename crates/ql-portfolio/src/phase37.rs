#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct MeucciBlackLittermanWoodburySpec<T: Real = f64> {
    pub covariance_3x3: [[T; 3]; 3],
    pub market_cap_weights: [T; 3],
    pub risk_aversion_delta: T,
    pub tau_scalar: T,
    pub view_pick_vector: [T; 3],
    pub view_target_q: T,
    pub view_uncertainty_omega: T,
}

#[derive(Clone, Copy, Debug)]
pub struct MeucciBlackLittermanWoodburyReport<T: Real = f64> {
    pub equilibrium_prior_returns_pi: [T; 3],
    pub posterior_returns_woodbury: [T; 3],
    pub posterior_returns_precision: [T; 3],
    pub woodbury_precision_equivalence_error: T,
    pub posterior_optimal_weights: [T; 3],
    pub active_view_weight_tilt_l1: T,
    pub posterior_view_expected_return: T,
}

fn inv3<T: Real>(a: [[T; 3]; 3]) -> [[T; 3]; 3] {
    let mut aug = [[T::zero(); 6]; 3];
    for i in 0..3 {
        for j in 0..3 {
            aug[i][j] = a[i][j];
        }
        aug[i][i + 3] = T::one();
    }
    for i in 0..3 {
        let piv = aug[i][i];
        for j in 0..6 {
            aug[i][j] = aug[i][j] / piv;
        }
        for r in 0..3 {
            if r != i {
                let f = aug[r][i];
                for j in 0..6 {
                    aug[r][j] = aug[r][j] - f * aug[i][j];
                }
            }
        }
    }
    let mut inv = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            inv[i][j] = aug[i][j + 3];
        }
    }
    inv
}

pub fn evaluate_meucci_black_litterman_woodbury_equivalence<T: Real>(
    spec: &MeucciBlackLittermanWoodburySpec<T>,
) -> MeucciBlackLittermanWoodburyReport<T> {
    let cov = &spec.covariance_3x3;
    let p = &spec.view_pick_vector;
    let tau = spec.tau_scalar.max_val(T::from_f64(1e-6));
    let omega = spec.view_uncertainty_omega.max_val(T::from_f64(1e-8));
    let delta = spec.risk_aversion_delta.max_val(T::from_f64(1e-6));

    let mut pi = [T::zero(); 3];
    let mut sig_p = [T::zero(); 3];
    for i in 0..3 {
        pi[i] = delta
            * (cov[i][0] * spec.market_cap_weights[0]
                + cov[i][1] * spec.market_cap_weights[1]
                + cov[i][2] * spec.market_cap_weights[2]);
        sig_p[i] = cov[i][0] * p[0] + cov[i][1] * p[1] + cov[i][2] * p[2];
    }
    let p_sig_p = p[0] * sig_p[0] + p[1] * sig_p[1] + p[2] * sig_p[2];
    let s_p = tau * p_sig_p + omega;
    let p_pi = p[0] * pi[0] + p[1] * pi[1] + p[2] * pi[2];
    let view_inn = spec.view_target_q - p_pi;

    let mut mu_wood = [T::zero(); 3];
    let mut m_wood = [[T::zero(); 3]; 3];
    for i in 0..3 {
        let k_i = tau * sig_p[i] / s_p;
        mu_wood[i] = pi[i] + k_i * view_inn;
        for j in 0..3 {
            m_wood[i][j] = tau * cov[i][j] - (tau * tau / s_p) * sig_p[i] * sig_p[j];
        }
    }

    let mut tau_cov = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            tau_cov[i][j] = tau * cov[i][j];
        }
    }
    let inv_tau_cov = inv3(tau_cov);
    let mut a_prec = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            a_prec[i][j] = inv_tau_cov[i][j] + (p[i] * p[j]) / omega;
        }
    }
    let m_prec = inv3(a_prec);
    let mut rhs = [T::zero(); 3];
    for i in 0..3 {
        let prior_prec_pi =
            inv_tau_cov[i][0] * pi[0] + inv_tau_cov[i][1] * pi[1] + inv_tau_cov[i][2] * pi[2];
        rhs[i] = prior_prec_pi + (spec.view_target_q / omega) * p[i];
    }
    let mut mu_prec = [T::zero(); 3];
    for i in 0..3 {
        mu_prec[i] = m_prec[i][0] * rhs[0] + m_prec[i][1] * rhs[1] + m_prec[i][2] * rhs[2];
    }

    let mut eq_err = T::zero();
    for i in 0..3 {
        let dm = mu_wood[i] - mu_prec[i];
        eq_err = eq_err + (dm * dm + T::from_f64(1e-30)).sqrt();
        for j in 0..3 {
            let dc = m_wood[i][j] - m_prec[i][j];
            eq_err = eq_err + (dc * dc + T::from_f64(1e-30)).sqrt();
        }
    }

    let inv_cov = inv3(*cov);
    let mut w_opt = [T::zero(); 3];
    let mut tilt_l1 = T::zero();
    for i in 0..3 {
        w_opt[i] =
            (inv_cov[i][0] * mu_wood[0] + inv_cov[i][1] * mu_wood[1] + inv_cov[i][2] * mu_wood[2])
                / delta;
        let dw = w_opt[i] - spec.market_cap_weights[i];
        tilt_l1 = tilt_l1 + (dw * dw + T::from_f64(1e-30)).sqrt();
    }
    let post_view_ret = p[0] * mu_wood[0] + p[1] * mu_wood[1] + p[2] * mu_wood[2];

    MeucciBlackLittermanWoodburyReport {
        equilibrium_prior_returns_pi: pi,
        posterior_returns_woodbury: mu_wood,
        posterior_returns_precision: mu_prec,
        woodbury_precision_equivalence_error: eq_err,
        posterior_optimal_weights: w_opt,
        active_view_weight_tilt_l1: tilt_l1,
        posterior_view_expected_return: post_view_ret,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct GueantAlmgrenChrissHyperbolicSpec<T: Real = f64> {
    pub initial_inventory_x0: T,
    pub horizon_years_t: T,
    pub temporary_impact_eta: T,
    pub permanent_impact_gamma: T,
    pub asset_volatility_sigma: T,
    pub risk_aversion_lambda: T,
}

#[derive(Clone, Copy, Debug)]
pub struct GueantAlmgrenChrissHyperbolicReport<T: Real = f64> {
    pub urgency_parameter_kappa: T,
    pub inventory_quartiles_q1_q2_q3: [T; 3],
    pub midpoint_frontloading_vs_twap: T,
    pub optimal_expected_shortfall: T,
    pub optimal_shortfall_variance: T,
    pub optimal_mean_variance_objective: T,
    pub twap_mean_variance_objective: T,
    pub ac_utility_gain_over_twap: T,
}

#[inline]
fn sinh_r<T: Real>(z: T) -> T {
    T::from_f64(0.5) * (z.exp() - (-z).exp())
}

pub fn evaluate_gueant_almgren_chriss_hyperbolic_liquidation<T: Real>(
    spec: &GueantAlmgrenChrissHyperbolicSpec<T>,
) -> GueantAlmgrenChrissHyperbolicReport<T> {
    let x0 = spec.initial_inventory_x0;
    let mat = spec.horizon_years_t.max_val(T::from_f64(1e-6));
    let eta = spec.temporary_impact_eta.max_val(T::from_f64(1e-8));
    let sig2 = spec.asset_volatility_sigma * spec.asset_volatility_sigma;
    let lam = spec.risk_aversion_lambda.max_val(T::from_f64(1e-8));
    let kappa = (lam * sig2 / eta).sqrt();
    let kt = (kappa * mat).max_val(T::from_f64(1e-5));
    let sh_kt = sinh_r(kt);

    let q1 = x0 * sinh_r(kt * T::from_f64(0.75)) / sh_kt;
    let q2 = x0 * sinh_r(kt * T::from_f64(0.50)) / sh_kt;
    let q3 = x0 * sinh_r(kt * T::from_f64(0.25)) / sh_kt;
    let frontload = T::from_f64(0.5) * x0 - q2;

    let two_kt = T::from_f64(2.0) * kt;
    let sh_2kt = sinh_r(two_kt);
    let four_sh2 = T::from_f64(4.0) * sh_kt * sh_kt;
    let x0_2 = x0 * x0;
    let perm_cost = T::from_f64(0.5) * spec.permanent_impact_gamma * x0_2;

    let e_is = perm_cost + eta * x0_2 * kappa * (sh_2kt + two_kt) / four_sh2;
    let var_is = sig2 * x0_2 * (sh_2kt - two_kt) / (kappa * four_sh2);
    let j_ac = e_is + T::from_f64(0.5) * lam * var_is;

    let e_twap = perm_cost + eta * x0_2 / mat;
    let var_twap = sig2 * x0_2 * mat / T::from_f64(3.0);
    let j_twap = e_twap + T::from_f64(0.5) * lam * var_twap;

    GueantAlmgrenChrissHyperbolicReport {
        urgency_parameter_kappa: kappa,
        inventory_quartiles_q1_q2_q3: [q1, q2, q3],
        midpoint_frontloading_vs_twap: frontload,
        optimal_expected_shortfall: e_is,
        optimal_shortfall_variance: var_is,
        optimal_mean_variance_objective: j_ac,
        twap_mean_variance_objective: j_twap,
        ac_utility_gain_over_twap: j_twap - j_ac,
    }
}

#[cfg(test)]
mod phase37_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn phase37_portfolio_bl_woodbury_and_ac_hyperbolic_invariants_and_ad_parity() {
        let h = 1e-6;

        let bl_spec = MeucciBlackLittermanWoodburySpec {
            covariance_3x3: [
                [0.040, 0.010, 0.004],
                [0.010, 0.025, 0.005],
                [0.004, 0.005, 0.010],
            ],
            market_cap_weights: [0.50, 0.30, 0.20],
            risk_aversion_delta: 2.5,
            tau_scalar: 0.05,
            view_pick_vector: [1.0, -1.0, 0.0],
            view_target_q: 0.045,
            view_uncertainty_omega: 0.0015,
        };
        let bl_rep = evaluate_meucci_black_litterman_woodbury_equivalence(&bl_spec);
        assert!(bl_rep.woodbury_precision_equivalence_error < 1e-10);
        assert!(bl_rep.active_view_weight_tilt_l1 > 0.0);

        let bl_ad = evaluate_meucci_black_litterman_woodbury_equivalence(
            &MeucciBlackLittermanWoodburySpec {
                covariance_3x3: bl_spec.covariance_3x3.map(|r| r.map(Dual::<1>::from_value)),
                market_cap_weights: bl_spec.market_cap_weights.map(Dual::<1>::from_value),
                risk_aversion_delta: Dual::<1>::from_value(bl_spec.risk_aversion_delta),
                tau_scalar: Dual::<1>::from_value(bl_spec.tau_scalar),
                view_pick_vector: bl_spec.view_pick_vector.map(Dual::<1>::from_value),
                view_target_q: Dual::<1>::variable(bl_spec.view_target_q, 0),
                view_uncertainty_omega: Dual::<1>::from_value(bl_spec.view_uncertainty_omega),
            },
        );
        let mut bp = bl_spec;
        bp.view_target_q += h;
        let mut bm = bl_spec;
        bm.view_target_q -= h;
        let fd_bl = (evaluate_meucci_black_litterman_woodbury_equivalence(&bp)
            .posterior_view_expected_return
            - evaluate_meucci_black_litterman_woodbury_equivalence(&bm)
                .posterior_view_expected_return)
            / (2.0 * h);
        assert!((bl_ad.posterior_view_expected_return.d[0] - fd_bl).abs() < 1e-5);

        let ac_spec = GueantAlmgrenChrissHyperbolicSpec {
            initial_inventory_x0: 100.0,
            horizon_years_t: 0.25,
            temporary_impact_eta: 0.002,
            permanent_impact_gamma: 0.0005,
            asset_volatility_sigma: 0.30,
            risk_aversion_lambda: 0.15,
        };
        let ac_rep = evaluate_gueant_almgren_chriss_hyperbolic_liquidation(&ac_spec);
        assert!(ac_rep.midpoint_frontloading_vs_twap > 0.0);
        assert!(ac_rep.ac_utility_gain_over_twap > 0.0);

        let ac_ad = evaluate_gueant_almgren_chriss_hyperbolic_liquidation(
            &GueantAlmgrenChrissHyperbolicSpec {
                initial_inventory_x0: Dual::<1>::variable(ac_spec.initial_inventory_x0, 0),
                horizon_years_t: Dual::<1>::from_value(ac_spec.horizon_years_t),
                temporary_impact_eta: Dual::<1>::from_value(ac_spec.temporary_impact_eta),
                permanent_impact_gamma: Dual::<1>::from_value(ac_spec.permanent_impact_gamma),
                asset_volatility_sigma: Dual::<1>::from_value(ac_spec.asset_volatility_sigma),
                risk_aversion_lambda: Dual::<1>::from_value(ac_spec.risk_aversion_lambda),
            },
        );
        let mut ap = ac_spec;
        ap.initial_inventory_x0 += h;
        let mut am = ac_spec;
        am.initial_inventory_x0 -= h;
        let fd_ac = (evaluate_gueant_almgren_chriss_hyperbolic_liquidation(&ap)
            .optimal_mean_variance_objective
            - evaluate_gueant_almgren_chriss_hyperbolic_liquidation(&am)
                .optimal_mean_variance_objective)
            / (2.0 * h);
        assert!((ac_ad.optimal_mean_variance_objective.d[0] - fd_ac).abs() < 1e-4);
    }
}
