#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct SchonbucherLandoRecoveryConventionSpec<T: Real = f64> {
    pub par_amount: T,
    pub risk_free_rate: T,
    pub hazard_rate_lambda: T,
    pub recovery_rate: T,
    pub coupon_rate: T,
    pub maturity_years: T,
}

#[derive(Clone, Copy, Debug)]
pub struct SchonbucherLandoRecoveryConventionReport<T: Real = f64> {
    pub rmv_zero_bond_price: T,
    pub rt_zero_bond_price: T,
    pub rfv_zero_bond_price: T,
    pub rmv_implied_spread_bps: T,
    pub rt_implied_spread_bps: T,
    pub rfv_implied_spread_bps: T,
    pub rt_minus_rmv_jensen_wedge: T,
    pub rfv_minus_rt_discounting_wedge: T,
    pub rfv_coupon_bond_fair_value: T,
}

pub fn evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle<T: Real>(
    spec: &SchonbucherLandoRecoveryConventionSpec<T>,
) -> SchonbucherLandoRecoveryConventionReport<T> {
    let f = spec.par_amount;
    let r = spec.risk_free_rate;
    let lam = spec.hazard_rate_lambda;
    let rec = spec.recovery_rate;
    let mat = spec.maturity_years.max_val(T::from_f64(1e-6));

    let s_rmv = (T::one() - rec) * lam;
    let b_rmv = f * (-(r + s_rmv) * mat).exp();

    let surv = (-lam * mat).exp();
    let surv_rt = surv + rec * (T::one() - surv);
    let b_rt = f * (-r * mat).exp() * surv_rt;

    let r_plus_l = (r + lam).max_val(T::from_f64(1e-8));
    let df_rl = (-r_plus_l * mat).exp();
    let b_rfv = f * df_rl + rec * f * (lam / r_plus_l) * (T::one() - df_rl);

    let bps = T::from_f64(1e4);
    let s_rmv_bps = s_rmv * bps;
    let s_rt_bps = -surv_rt.ln() / mat * bps;
    let s_rfv_bps = -(b_rfv / (f * (-r * mat).exp())).ln() / mat * bps;

    let cpn_pv = f * spec.coupon_rate * (T::one() - df_rl) / r_plus_l;

    SchonbucherLandoRecoveryConventionReport {
        rmv_zero_bond_price: b_rmv,
        rt_zero_bond_price: b_rt,
        rfv_zero_bond_price: b_rfv,
        rmv_implied_spread_bps: s_rmv_bps,
        rt_implied_spread_bps: s_rt_bps,
        rfv_implied_spread_bps: s_rfv_bps,
        rt_minus_rmv_jensen_wedge: b_rt - b_rmv,
        rfv_minus_rt_discounting_wedge: b_rfv - b_rt,
        rfv_coupon_bond_fair_value: b_rfv + cpn_pv,
    }
}

#[cfg(test)]
mod phase37_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn schonbucher_lando_rmv_rt_rfv_invariants_and_ad_parity() {
        let spec = SchonbucherLandoRecoveryConventionSpec {
            par_amount: 100.0,
            risk_free_rate: 0.04,
            hazard_rate_lambda: 0.035,
            recovery_rate: 0.40,
            coupon_rate: 0.05,
            maturity_years: 5.0,
        };
        let rep = evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle(&spec);
        assert!(rep.rfv_zero_bond_price > rep.rt_zero_bond_price);
        assert!(rep.rt_zero_bond_price > rep.rmv_zero_bond_price);
        assert!(rep.rt_minus_rmv_jensen_wedge > 0.0);
        assert!(rep.rfv_minus_rt_discounting_wedge > 0.0);

        let h = 1e-6;
        let r_ad = evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle(
            &SchonbucherLandoRecoveryConventionSpec {
                par_amount: Dual::<1>::from_value(spec.par_amount),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                hazard_rate_lambda: Dual::<1>::variable(spec.hazard_rate_lambda, 0),
                recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
                coupon_rate: Dual::<1>::from_value(spec.coupon_rate),
                maturity_years: Dual::<1>::from_value(spec.maturity_years),
            },
        );
        let mut sp = spec;
        sp.hazard_rate_lambda += h;
        let mut sm = spec;
        sm.hazard_rate_lambda -= h;
        let fd = (evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle(&sp)
            .rfv_coupon_bond_fair_value
            - evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle(&sm)
                .rfv_coupon_bond_fair_value)
            / (2.0 * h);
        assert!((r_ad.rfv_coupon_bond_fair_value.d[0] - fd).abs() < 1e-5);
    }
}
