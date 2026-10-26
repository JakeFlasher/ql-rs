#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct LjungqvistNewKeynesianReeSpec<T: Real = f64> {
    pub risk_aversion_sigma: T,
    pub discount_factor_beta: T,
    pub calvo_phillips_kappa: T,
    pub taylor_inflation_phi_pi: T,
    pub taylor_output_phi_y: T,
    pub natural_rate_persistence_rho: T,
    pub natural_rate_shock_rn: T,
}

#[derive(Clone, Copy, Debug)]
pub struct LjungqvistNewKeynesianReeReport<T: Real = f64> {
    pub blanchard_kahn_determinacy_margin: T,
    pub output_gap_loading_psi_y: T,
    pub inflation_loading_psi_pi: T,
    pub equilibrium_output_gap: T,
    pub equilibrium_inflation: T,
    pub equilibrium_policy_rate: T,
    pub ex_ante_real_interest_rate: T,
    pub ree_is_nkpc_conservation_error: T,
}

pub fn evaluate_ljungqvist_new_keynesian_three_equation_ree<T: Real>(
    spec: &LjungqvistNewKeynesianReeSpec<T>,
) -> LjungqvistNewKeynesianReeReport<T> {
    let sig = spec.risk_aversion_sigma.max_val(T::from_f64(1e-6));
    let beta = spec.discount_factor_beta;
    let kappa = spec.calvo_phillips_kappa;
    let phi_pi = spec.taylor_inflation_phi_pi;
    let phi_y = spec.taylor_output_phi_y;
    let rho = spec.natural_rate_persistence_rho;
    let rn = spec.natural_rate_shock_rn;

    let bk_margin = kappa * (phi_pi - T::one()) + (T::one() - beta) * phi_y;
    let one_m_brho = T::one() - beta * rho;
    let lam_nk = ((sig * (T::one() - rho) + phi_y) * one_m_brho + kappa * (phi_pi - rho))
        .max_val(T::from_f64(1e-8));
    let psi_y = one_m_brho / lam_nk;
    let psi_pi = kappa / lam_nk;

    let y_t = psi_y * rn;
    let pi_t = psi_pi * rn;
    let i_t = phi_pi * pi_t + phi_y * y_t;
    let ey_tp1 = rho * y_t;
    let epi_tp1 = rho * pi_t;
    let real_r = i_t - epi_tp1;

    let res_is = y_t - (ey_tp1 - (real_r - rn) / sig);
    let res_nkpc = pi_t - (beta * epi_tp1 + kappa * y_t);
    let err = (res_is * res_is + res_nkpc * res_nkpc + T::from_f64(1e-30)).sqrt();

    LjungqvistNewKeynesianReeReport {
        blanchard_kahn_determinacy_margin: bk_margin,
        output_gap_loading_psi_y: psi_y,
        inflation_loading_psi_pi: psi_pi,
        equilibrium_output_gap: y_t,
        equilibrium_inflation: pi_t,
        equilibrium_policy_rate: i_t,
        ex_ante_real_interest_rate: real_r,
        ree_is_nkpc_conservation_error: err,
    }
}

#[cfg(test)]
mod phase37_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn ljungqvist_new_keynesian_ree_invariants_and_ad_parity() {
        let spec = LjungqvistNewKeynesianReeSpec {
            risk_aversion_sigma: 1.5,
            discount_factor_beta: 0.99,
            calvo_phillips_kappa: 0.15,
            taylor_inflation_phi_pi: 1.50,
            taylor_output_phi_y: 0.125,
            natural_rate_persistence_rho: 0.80,
            natural_rate_shock_rn: 0.015,
        };
        let rep = evaluate_ljungqvist_new_keynesian_three_equation_ree(&spec);
        assert!(rep.blanchard_kahn_determinacy_margin > 0.0);
        assert!(rep.ree_is_nkpc_conservation_error < 1e-12);

        let h = 1e-6;
        let r_ad =
            evaluate_ljungqvist_new_keynesian_three_equation_ree(&LjungqvistNewKeynesianReeSpec {
                risk_aversion_sigma: Dual::<1>::from_value(spec.risk_aversion_sigma),
                discount_factor_beta: Dual::<1>::from_value(spec.discount_factor_beta),
                calvo_phillips_kappa: Dual::<1>::from_value(spec.calvo_phillips_kappa),
                taylor_inflation_phi_pi: Dual::<1>::variable(spec.taylor_inflation_phi_pi, 0),
                taylor_output_phi_y: Dual::<1>::from_value(spec.taylor_output_phi_y),
                natural_rate_persistence_rho: Dual::<1>::from_value(
                    spec.natural_rate_persistence_rho,
                ),
                natural_rate_shock_rn: Dual::<1>::from_value(spec.natural_rate_shock_rn),
            });
        let mut sp = spec;
        sp.taylor_inflation_phi_pi += h;
        let mut sm = spec;
        sm.taylor_inflation_phi_pi -= h;
        let fd = (evaluate_ljungqvist_new_keynesian_three_equation_ree(&sp)
            .equilibrium_policy_rate
            - evaluate_ljungqvist_new_keynesian_three_equation_ree(&sm).equilibrium_policy_rate)
            / (2.0 * h);
        assert!((r_ad.equilibrium_policy_rate.d[0] - fd).abs() < 1e-5);
    }
}
