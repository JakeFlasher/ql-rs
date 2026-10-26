use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct TuckmanRepoSpecialnessBasisSpec<T: Real = f64> {
    pub otr_dirty_price: T,
    pub otr_yield: T,
    pub oftr_yield: T,
    pub modified_duration: T,
    pub gc_repo_rate: T,
    pub special_repo_rate: T,
    pub specialness_horizon_years: T,
    pub repo_haircut_fraction: T,
    pub equity_funding_rate: T,
}
#[derive(Clone, Copy, Debug)]
pub struct TuckmanRepoSpecialnessBasisReport<T: Real = f64> {
    pub repo_specialness_spread_bps: T,
    pub repo_specialness_pv_per_100: T,
    pub total_oftr_otr_yield_spread_bps: T,
    pub repo_specialness_yield_equiv_bps: T,
    pub structural_liquidity_yield_spread_bps: T,
    pub spread_decomposition_conservation_error: T,
    pub levered_basis_net_carry_bps: T,
    pub haircut_funded_levered_roe: T,
}
pub fn evaluate_tuckman_repo_specialness_and_oftr_basis<T: Real>(
    s: &TuckmanRepoSpecialnessBasisSpec<T>,
) -> TuckmanRepoSpecialnessBasisReport<T> {
    let d = s.gc_repo_rate - s.special_repo_rate;
    let pv = s.otr_dirty_price
        * d
        * s.specialness_horizon_years
        * (-s.gc_repo_rate * s.specialness_horizon_years).exp();
    let total = (s.oftr_yield - s.otr_yield) * T::from_f64(1e4);
    let repo = pv / (s.otr_dirty_price * s.modified_duration) * T::from_f64(1e4);
    let liq = total - repo;
    let carry = (s.oftr_yield - s.otr_yield) + (s.gc_repo_rate - s.special_repo_rate)
        - T::from_f64(2.0) * s.repo_haircut_fraction * (s.equity_funding_rate - s.gc_repo_rate);
    let roe = carry / (T::from_f64(2.0) * s.repo_haircut_fraction.max_val(T::from_f64(1e-12)));
    TuckmanRepoSpecialnessBasisReport {
        repo_specialness_spread_bps: d * T::from_f64(1e4),
        repo_specialness_pv_per_100: pv * T::from_f64(100.0) / s.otr_dirty_price,
        total_oftr_otr_yield_spread_bps: total,
        repo_specialness_yield_equiv_bps: repo,
        structural_liquidity_yield_spread_bps: liq,
        spread_decomposition_conservation_error: {
            let diff = total - repo - liq;
            (diff * diff + T::from_f64(1e-30)).sqrt()
        },
        levered_basis_net_carry_bps: carry * T::from_f64(1e4),
        haircut_funded_levered_roe: roe,
    }
}

#[cfg(test)]
mod phase34_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> TuckmanRepoSpecialnessBasisSpec<T> {
        TuckmanRepoSpecialnessBasisSpec {
            otr_dirty_price: T::from_f64(100.0) + x,
            otr_yield: T::from_f64(0.04),
            oftr_yield: T::from_f64(0.0415),
            modified_duration: T::from_f64(7.5),
            gc_repo_rate: T::from_f64(0.05) + x * T::from_f64(0.01),
            special_repo_rate: T::from_f64(0.035),
            specialness_horizon_years: T::from_f64(0.25),
            repo_haircut_fraction: T::from_f64(0.02),
            equity_funding_rate: T::from_f64(0.06),
        }
    }
    #[test]
    fn phase34_curves_invariants_and_dual_fd() {
        let r = evaluate_tuckman_repo_specialness_and_oftr_basis(&spec(0.0));
        assert!(
            r.spread_decomposition_conservation_error < 1e-12
                && r.repo_specialness_spread_bps > 0.0
        );
        let h = 1e-6;
        let fd = (evaluate_tuckman_repo_specialness_and_oftr_basis(&spec(h))
            .haircut_funded_levered_roe
            - evaluate_tuckman_repo_specialness_and_oftr_basis(&spec(-h))
                .haircut_funded_levered_roe)
            / (2.0 * h);
        let ad =
            evaluate_tuckman_repo_specialness_and_oftr_basis(&spec(Dual::<1>::variable(0.0, 0)))
                .haircut_funded_levered_roe;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
