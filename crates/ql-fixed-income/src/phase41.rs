//! Tuckman & Serrat (2011) & Fabozzi (2021) Special Repo Rate (`GC` vs. `Special`),
//! TMPG Fail-Penalty Floor & Cheapest-to-Deliver (`CTD`) Implied-Repo Switch Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct TuckmanSpecialRepoFailPenaltyCtdSpec<T: Real = f64> {
    pub general_collateral_rate_gc: T,
    pub tmpg_fail_penalty_rate_p: T,
    pub issue_short_demand_ratio_d: T,
    pub lendable_supply_elasticity_k: T,
    pub bond_clean_price: T,
    pub accrued_interest_at_spot: T,
    pub coupon_rate_annual: T,
    pub futures_invoice_price_at_delivery: T,
    pub holding_period_year_frac_tau: T,
}

#[derive(Debug, Clone, Copy)]
pub struct TuckmanSpecialRepoFailPenaltyCtdReport<T: Real = f64> {
    pub equilibrium_special_repo_rate: T,
    pub repo_specialness_spread_bps: T,
    pub tmpg_fail_floor_rate: T,
    pub financing_advantage_pv_per_100: T,
    pub implied_repo_rate_irr: T,
    pub net_basis_after_specialness: T,
    pub fail_probability_index: T,
}

pub fn evaluate_tuckman_special_repo_fail_penalty_ctd<T: Real>(
    spec: &TuckmanSpecialRepoFailPenaltyCtdSpec<T>,
) -> TuckmanSpecialRepoFailPenaltyCtdReport<T> {
    let gc = spec.general_collateral_rate_gc;
    let pen = spec.tmpg_fail_penalty_rate_p.max_val(T::from_f64(0.01));
    let tau = spec.holding_period_year_frac_tau.max_val(T::from_f64(0.01));

    // TMPG fail floor rate = -penalty (or GC - penalty in scarcity equilibrium):
    let fail_floor = gc - pen;
    let scarcity = (spec.issue_short_demand_ratio_d
        * spec.lendable_supply_elasticity_k.max_val(T::from_f64(0.01)))
    .max_val(T::zero());
    let fail_prob = T::one() - (-scarcity).exp();

    // Equilibrium special repo rate interpolates smoothly between GC and the TMPG fail floor:
    let special_rate = gc - pen * fail_prob;
    let specialness_bps = (gc - special_rate) * T::from_f64(1e4);

    let dirty_spot = (spec.bond_clean_price + spec.accrued_interest_at_spot).max_val(T::one());
    let fin_adv_pv = dirty_spot * (gc - special_rate) * tau;

    // Implied repo rate (IRR) from cash-and-carry to futures delivery:
    let carry_coupon = T::from_f64(100.0) * spec.coupon_rate_annual * tau;
    let proceeds = spec.futures_invoice_price_at_delivery + carry_coupon;
    let irr = ((proceeds - dirty_spot) / dirty_spot) / tau;

    // Net basis = DirtySpot * (1 + r_special * tau) - CarryCoupon - InvoicePrice:
    let net_basis = dirty_spot * (T::one() + special_rate * tau)
        - carry_coupon
        - spec.futures_invoice_price_at_delivery;

    TuckmanSpecialRepoFailPenaltyCtdReport {
        equilibrium_special_repo_rate: special_rate,
        repo_specialness_spread_bps: specialness_bps,
        tmpg_fail_floor_rate: fail_floor,
        financing_advantage_pv_per_100: fin_adv_pv,
        implied_repo_rate_irr: irr,
        net_basis_after_specialness: net_basis,
        fail_probability_index: fail_prob,
    }
}

#[cfg(test)]
mod phase41_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn tuckman_special_repo_fail_penalty_ctd_invariants_and_ad_parity() {
        let spec = TuckmanSpecialRepoFailPenaltyCtdSpec {
            general_collateral_rate_gc: 0.048,
            tmpg_fail_penalty_rate_p: 0.030,
            issue_short_demand_ratio_d: 1.4,
            lendable_supply_elasticity_k: 0.55,
            bond_clean_price: 98.5,
            accrued_interest_at_spot: 0.85,
            coupon_rate_annual: 0.040,
            futures_invoice_price_at_delivery: 99.65,
            holding_period_year_frac_tau: 0.25,
        };
        let rep = evaluate_tuckman_special_repo_fail_penalty_ctd(&spec);
        assert!(rep.equilibrium_special_repo_rate < spec.general_collateral_rate_gc);
        assert!(rep.equilibrium_special_repo_rate > rep.tmpg_fail_floor_rate);
        assert!(rep.repo_specialness_spread_bps > 0.0);

        let h = 1e-6;
        let r_ad =
            evaluate_tuckman_special_repo_fail_penalty_ctd(&TuckmanSpecialRepoFailPenaltyCtdSpec {
                general_collateral_rate_gc: Dual::<1>::variable(spec.general_collateral_rate_gc, 0),
                tmpg_fail_penalty_rate_p: Dual::<1>::from_value(spec.tmpg_fail_penalty_rate_p),
                issue_short_demand_ratio_d: Dual::<1>::from_value(spec.issue_short_demand_ratio_d),
                lendable_supply_elasticity_k: Dual::<1>::from_value(
                    spec.lendable_supply_elasticity_k,
                ),
                bond_clean_price: Dual::<1>::from_value(spec.bond_clean_price),
                accrued_interest_at_spot: Dual::<1>::from_value(spec.accrued_interest_at_spot),
                coupon_rate_annual: Dual::<1>::from_value(spec.coupon_rate_annual),
                futures_invoice_price_at_delivery: Dual::<1>::from_value(
                    spec.futures_invoice_price_at_delivery,
                ),
                holding_period_year_frac_tau: Dual::<1>::from_value(
                    spec.holding_period_year_frac_tau,
                ),
            });
        let mut sp = spec;
        sp.general_collateral_rate_gc += h;
        let mut sm = spec;
        sm.general_collateral_rate_gc -= h;
        let fd = (evaluate_tuckman_special_repo_fail_penalty_ctd(&sp).net_basis_after_specialness
            - evaluate_tuckman_special_repo_fail_penalty_ctd(&sm).net_basis_after_specialness)
            / (2.0 * h);
        assert!((r_ad.net_basis_after_specialness.d[0] - fd).abs() < 1e-6);
    }
}
