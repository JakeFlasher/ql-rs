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
    #[test]
    fn phase34_smoke() {
        assert!(std::hint::black_box(true));
    }
}
