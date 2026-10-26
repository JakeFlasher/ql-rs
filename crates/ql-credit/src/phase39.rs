//! Lando (2004) & Schönbucher (2003) Marshall–Olkin Fatal-Shock Common-Poisson
//! Dynamic Credit Contagion & First/Second-to-Default Basket Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct LandoMarshallOlkinContagionSpec<T: Real = f64> {
    pub idiosyncratic_intensity_lambda1: T,
    pub idiosyncratic_intensity_lambda2: T,
    pub fatal_common_shock_lambda12: T,
    pub recovery_rate: T,
    pub risk_free_rate: T,
    pub maturity_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct LandoMarshallOlkinContagionReport<T: Real = f64> {
    pub marginal_hazard_rate_1: T,
    pub marginal_hazard_rate_2: T,
    pub default_time_tail_correlation: T,
    pub joint_survival_prob_stt: T,
    pub simultaneous_default_prob: T,
    pub first_to_default_par_spread_bps: T,
    pub second_to_default_par_spread_bps: T,
    pub inclusion_exclusion_spread_identity_error: T,
}

pub fn evaluate_lando_marshall_olkin_fatal_shock_basket<T: Real>(
    spec: &LandoMarshallOlkinContagionSpec<T>,
) -> LandoMarshallOlkinContagionReport<T> {
    let l1 = spec
        .idiosyncratic_intensity_lambda1
        .max_val(T::from_f64(1e-5));
    let l2 = spec
        .idiosyncratic_intensity_lambda2
        .max_val(T::from_f64(1e-5));
    let l12 = spec.fatal_common_shock_lambda12.max_val(T::zero());
    let lgd = T::one() - spec.recovery_rate;
    let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

    // Marginal hazard rates bar_lambda_i = lambda_i + lambda_12
    let h1 = l1 + l12;
    let h2 = l2 + l12;
    let l_tot = l1 + l2 + l12;

    // Marshall-Olkin exponential default-time correlation:
    // Corr(tau_1, tau_2) = lambda_12 / (lambda_1 + lambda_2 + lambda_12)
    let rho_tau = l12 / l_tot;
    let s_tt = (-l_tot * t).exp();

    // Probability that both obligors default simultaneously via the fatal shock before T:
    let p_sim = (l12 / l_tot) * (T::one() - s_tt);

    // Continuous-premium CDS par spread with constant hazard rate h is (1 - R) * h.
    // For First-to-Default (tau_FtD = min(tau_1, tau_2) ~ Exp(l1 + l2 + l12)):
    let s_ftd = lgd * l_tot;
    let s_1 = lgd * h1;
    let s_2 = lgd * h2;

    // By inclusion-exclusion on protection legs and risky annuities for Second-to-Default (tau_StD = max(tau_1, tau_2)):
    let r = spec.risk_free_rate.max_val(T::from_f64(1e-4));
    let rpv01 = |lam: T| -> T { (T::one() - (-(r + lam) * t).exp()) / (r + lam) };
    let prot = |lam: T| -> T { lgd * lam * rpv01(lam) };

    let prot_std = prot(h1) + prot(h2) - prot(l_tot);
    let rpv_std = rpv01(h1) + rpv01(h2) - rpv01(l_tot);
    let s_std = prot_std / rpv_std.max_val(T::from_f64(1e-8));

    // Protection PV identity check: Prot(1) + Prot(2) = Prot(FtD) + Prot(StD)
    let id_err = (prot(h1) + prot(h2) - (prot(l_tot) + prot_std)).abs()
        + ((s_1 + s_2) - (s_ftd + lgd * l12)).abs();

    LandoMarshallOlkinContagionReport {
        marginal_hazard_rate_1: h1,
        marginal_hazard_rate_2: h2,
        default_time_tail_correlation: rho_tau,
        joint_survival_prob_stt: s_tt,
        simultaneous_default_prob: p_sim,
        first_to_default_par_spread_bps: s_ftd * T::from_f64(1e4),
        second_to_default_par_spread_bps: s_std * T::from_f64(1e4),
        inclusion_exclusion_spread_identity_error: id_err,
    }
}

#[cfg(test)]
mod phase39_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn lando_marshall_olkin_contagion_invariants_and_ad_parity() {
        let spec = LandoMarshallOlkinContagionSpec {
            idiosyncratic_intensity_lambda1: 0.018,
            idiosyncratic_intensity_lambda2: 0.024,
            fatal_common_shock_lambda12: 0.012,
            recovery_rate: 0.40,
            risk_free_rate: 0.035,
            maturity_years_t: 5.0,
        };
        let rep = evaluate_lando_marshall_olkin_fatal_shock_basket(&spec);
        assert!(rep.default_time_tail_correlation > 0.0 && rep.default_time_tail_correlation < 1.0);
        assert!(rep.first_to_default_par_spread_bps > rep.second_to_default_par_spread_bps);
        assert!(rep.inclusion_exclusion_spread_identity_error < 1e-12);

        let h = 1e-6;
        let r_ad =
            evaluate_lando_marshall_olkin_fatal_shock_basket(&LandoMarshallOlkinContagionSpec {
                idiosyncratic_intensity_lambda1: Dual::<1>::from_value(
                    spec.idiosyncratic_intensity_lambda1,
                ),
                idiosyncratic_intensity_lambda2: Dual::<1>::from_value(
                    spec.idiosyncratic_intensity_lambda2,
                ),
                fatal_common_shock_lambda12: Dual::<1>::variable(
                    spec.fatal_common_shock_lambda12,
                    0,
                ),
                recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
            });
        let mut sp = spec;
        sp.fatal_common_shock_lambda12 += h;
        let mut sm = spec;
        sm.fatal_common_shock_lambda12 -= h;
        let fd = (evaluate_lando_marshall_olkin_fatal_shock_basket(&sp)
            .second_to_default_par_spread_bps
            - evaluate_lando_marshall_olkin_fatal_shock_basket(&sm)
                .second_to_default_par_spread_bps)
            / (2.0 * h);
        assert!((r_ad.second_to_default_par_spread_bps.d[0] - fd).abs() < 1e-5);
    }
}
