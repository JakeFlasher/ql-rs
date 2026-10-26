use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct RomerVeronesiTaylorRuleMacroCurveSpec<T: Real = f64> {
    pub natural_real_rate_r_star: T,
    pub inflation_target_pi_star: T,
    pub current_inflation: T,
    pub current_output_gap: T,
    pub taylor_inflation_weight_phi_pi: T,
    pub taylor_output_weight_phi_y: T,
    pub inflation_mean_reversion: T,
    pub output_gap_mean_reversion: T,
    pub inflation_vol: T,
    pub output_gap_vol: T,
    pub market_price_of_inflation_risk: T,
    pub zlb_floor_rate: T,
    pub maturities: Vec<T>,
}
#[derive(Clone, Debug)]
pub struct RomerVeronesiTaylorRuleMacroCurveReport<T: Real = f64> {
    pub unconstrained_taylor_shadow_rate: T,
    pub zlb_constrained_policy_rate: T,
    pub macro_affine_zero_yields: Vec<T>,
    pub expectations_hypothesis_yields: Vec<T>,
    pub macro_term_premia_bps: Vec<T>,
    pub yield_inflation_beta_10y: T,
    pub yield_output_gap_beta_10y: T,
}
pub fn evaluate_romer_veronesi_taylor_rule_macro_yield_curve<T: Real>(
    s: &RomerVeronesiTaylorRuleMacroCurveSpec<T>,
) -> RomerVeronesiTaylorRuleMacroCurveReport<T> {
    let shadow = s.natural_real_rate_r_star
        + s.inflation_target_pi_star
        + s.taylor_inflation_weight_phi_pi * (s.current_inflation - s.inflation_target_pi_star)
        + s.taylor_output_weight_phi_y * s.current_output_gap;
    let policy = shadow.max_val(s.zlb_floor_rate);
    let mut affine = Vec::with_capacity(s.maturities.len());
    let mut eh = Vec::with_capacity(s.maturities.len());
    let mut premia = Vec::with_capacity(s.maturities.len());
    for &t in &s.maturities {
        let ti = (-s.inflation_mean_reversion * t).exp();
        let ty = (-s.output_gap_mean_reversion * t).exp();
        let inflation_state =
            s.inflation_target_pi_star + (s.current_inflation - s.inflation_target_pi_star) * ti;
        let output_state = s.current_output_gap * ty;
        let expected = s.natural_real_rate_r_star
            + inflation_state
            + s.taylor_inflation_weight_phi_pi * (inflation_state - s.inflation_target_pi_star)
            + s.taylor_output_weight_phi_y * output_state;
        let risk =
            s.market_price_of_inflation_risk * s.inflation_vol * s.inflation_vol * (T::one() - ti)
                + s.output_gap_vol * s.output_gap_vol * (T::one() - ty) * T::from_f64(0.5);
        let y = expected + risk;
        affine.push(y);
        eh.push(expected);
        premia.push(risk * T::from_f64(10_000.0));
    }
    let horizon = T::from_f64(10.0);
    let bi = s.taylor_inflation_weight_phi_pi * (-s.inflation_mean_reversion * horizon).exp();
    let by = s.taylor_output_weight_phi_y * (-s.output_gap_mean_reversion * horizon).exp();
    RomerVeronesiTaylorRuleMacroCurveReport {
        unconstrained_taylor_shadow_rate: shadow,
        zlb_constrained_policy_rate: policy,
        macro_affine_zero_yields: affine,
        expectations_hypothesis_yields: eh,
        macro_term_premia_bps: premia,
        yield_inflation_beta_10y: bi,
        yield_output_gap_beta_10y: by,
    }
}

#[cfg(test)]
mod phase30_tests {
    use super::*;
    use ql_ad::Dual;
    fn make<T: Real>(x: T) -> RomerVeronesiTaylorRuleMacroCurveSpec<T> {
        RomerVeronesiTaylorRuleMacroCurveSpec {
            natural_real_rate_r_star: T::from_f64(0.02) + x,
            inflation_target_pi_star: T::from_f64(0.02),
            current_inflation: T::from_f64(0.025),
            current_output_gap: T::from_f64(0.01),
            taylor_inflation_weight_phi_pi: T::from_f64(1.5),
            taylor_output_weight_phi_y: T::from_f64(0.5),
            inflation_mean_reversion: T::from_f64(0.4),
            output_gap_mean_reversion: T::from_f64(0.3),
            inflation_vol: T::from_f64(0.2),
            output_gap_vol: T::from_f64(0.1),
            market_price_of_inflation_risk: T::from_f64(0.2),
            zlb_floor_rate: T::zero(),
            maturities: vec![T::one(), T::from_f64(10.0)],
        }
    }
    #[test]
    fn phase30_taylor_invariants_and_dual_fd() {
        let r = evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&make(0.0));
        assert!(r.zlb_constrained_policy_rate >= 0.0);
        let h = 1e-6;
        let fd = (evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&make(h))
            .unconstrained_taylor_shadow_rate
            - evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&make(-h))
                .unconstrained_taylor_shadow_rate)
            / (2.0 * h);
        let ad = evaluate_romer_veronesi_taylor_rule_macro_yield_curve(&make(Dual::<1>::variable(
            0.0, 0,
        )))
        .unconstrained_taylor_shadow_rate;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
