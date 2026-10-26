//! Haug (2007) & Hull–White (2004) Executive Stock Option (`ESO`) Vesting-Forfeiture,
//! Suboptimal Exercise Multiple ($M$-Barrier) & Warrant Dilution Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct HaugEsoVestingBarrierDilutionSpec<T: Real = f64> {
    pub spot_price_s0: T,
    pub strike_price_k: T,
    pub risk_free_rate_r: T,
    pub dividend_yield_q: T,
    pub volatility_sigma: T,
    pub contractual_maturity_t: T,
    pub vesting_period_tv: T,
    pub forfeiture_hazard_lambda_f: T,
    pub suboptimal_exercise_multiple_m: T,
    pub grant_to_shares_ratio_gamma: T,
}

#[derive(Debug, Clone, Copy)]
pub struct HaugEsoVestingBarrierDilutionReport<T: Real = f64> {
    pub vesting_survival_probability: T,
    pub effective_expected_life_teff: T,
    pub plain_black_scholes_call: T,
    pub forfeiture_adjusted_eso_value: T,
    pub diluted_per_option_grant_cost: T,
    pub shareholder_dilution_discount_factor: T,
}

pub fn evaluate_haug_eso_vesting_barrier_dilution<T: Real>(
    spec: &HaugEsoVestingBarrierDilutionSpec<T>,
) -> HaugEsoVestingBarrierDilutionReport<T> {
    let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
    let k = spec.strike_price_k.max_val(T::from_f64(1e-4));
    let r = spec.risk_free_rate_r;
    let q = spec.dividend_yield_q;
    let v = spec.volatility_sigma.max_val(T::from_f64(1e-4));
    let t = spec.contractual_maturity_t.max_val(T::from_f64(0.5));
    let tv = spec.vesting_period_tv.max_val(T::from_f64(0.1)).min_val(t);
    let lam_f = spec.forfeiture_hazard_lambda_f.max_val(T::zero());
    let m_mult = spec
        .suboptimal_exercise_multiple_m
        .max_val(T::from_f64(1.05));

    // Pre-vesting survival probability exp(-lambda_f * T_v):
    let surv_prob = (-lam_f * tv).exp();

    // Hull-White (2004) / FAS 123R effective life reduction from voluntary exercise barrier S^* = M * K:
    let log_m = m_mult.ln();
    let barrier_speed = (v * v + T::from_f64(2.0) * log_m / (t - tv).max_val(T::from_f64(0.25)))
        .max_val(T::from_f64(0.05));
    let t_eff = tv + (t - tv) / (T::one() + T::from_f64(0.5) * log_m * barrier_speed);

    // Plain Black-Scholes at contractual maturity T vs. ESO at effective life T_eff:
    let std_t = v * t.sqrt();
    let d1_t = ((s0 / k).ln() + (r - q + T::from_f64(0.5) * v * v) * t) / std_t;
    let d2_t = d1_t - std_t;
    let bs_full = s0 * (-q * t).exp() * norm_cdf(d1_t) - k * (-r * t).exp() * norm_cdf(d2_t);

    let std_e = v * t_eff.sqrt();
    let d1_e = ((s0 / k).ln() + (r - q + T::from_f64(0.5) * v * v) * t_eff) / std_e;
    let d2_e = d1_e - std_e;
    let bs_eff = s0 * (-q * t_eff).exp() * norm_cdf(d1_e) - k * (-r * t_eff).exp() * norm_cdf(d2_e);

    let eso_undiluted = surv_prob * bs_eff;
    let dil_factor = T::one() / (T::one() + spec.grant_to_shares_ratio_gamma.max_val(T::zero()));
    let eso_diluted = eso_undiluted * dil_factor;

    HaugEsoVestingBarrierDilutionReport {
        vesting_survival_probability: surv_prob,
        effective_expected_life_teff: t_eff,
        plain_black_scholes_call: bs_full,
        forfeiture_adjusted_eso_value: eso_undiluted,
        diluted_per_option_grant_cost: eso_diluted,
        shareholder_dilution_discount_factor: dil_factor,
    }
}

#[cfg(test)]
mod phase41_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn haug_eso_vesting_barrier_dilution_invariants_and_ad_parity() {
        let spec = HaugEsoVestingBarrierDilutionSpec {
            spot_price_s0: 100.0,
            strike_price_k: 100.0,
            risk_free_rate_r: 0.04,
            dividend_yield_q: 0.015,
            volatility_sigma: 0.30,
            contractual_maturity_t: 10.0,
            vesting_period_tv: 3.0,
            forfeiture_hazard_lambda_f: 0.06,
            suboptimal_exercise_multiple_m: 1.8,
            grant_to_shares_ratio_gamma: 0.05,
        };
        let rep = evaluate_haug_eso_vesting_barrier_dilution(&spec);
        assert!(rep.effective_expected_life_teff < spec.contractual_maturity_t);
        assert!(rep.forfeiture_adjusted_eso_value < rep.plain_black_scholes_call);
        assert!(rep.diluted_per_option_grant_cost < rep.forfeiture_adjusted_eso_value);

        let h = 1e-6;
        let r_ad = evaluate_haug_eso_vesting_barrier_dilution(&HaugEsoVestingBarrierDilutionSpec {
            spot_price_s0: Dual::<1>::variable(spec.spot_price_s0, 0),
            strike_price_k: Dual::<1>::from_value(spec.strike_price_k),
            risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
            dividend_yield_q: Dual::<1>::from_value(spec.dividend_yield_q),
            volatility_sigma: Dual::<1>::from_value(spec.volatility_sigma),
            contractual_maturity_t: Dual::<1>::from_value(spec.contractual_maturity_t),
            vesting_period_tv: Dual::<1>::from_value(spec.vesting_period_tv),
            forfeiture_hazard_lambda_f: Dual::<1>::from_value(spec.forfeiture_hazard_lambda_f),
            suboptimal_exercise_multiple_m: Dual::<1>::from_value(
                spec.suboptimal_exercise_multiple_m,
            ),
            grant_to_shares_ratio_gamma: Dual::<1>::from_value(spec.grant_to_shares_ratio_gamma),
        });
        let mut sp = spec;
        sp.spot_price_s0 += h;
        let mut sm = spec;
        sm.spot_price_s0 -= h;
        let fd = (evaluate_haug_eso_vesting_barrier_dilution(&sp).diluted_per_option_grant_cost
            - evaluate_haug_eso_vesting_barrier_dilution(&sm).diluted_per_option_grant_cost)
            / (2.0 * h);
        assert!((r_ad.diluted_per_option_grant_cost.d[0] - fd).abs() < 1e-5);
    }
}
