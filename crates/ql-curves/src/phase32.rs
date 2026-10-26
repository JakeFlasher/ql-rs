use ql_ad::Real;
#[derive(Clone, Debug)]
pub struct AdbRbiNdirsAssetSwapSpec<T: Real = f64> {
    pub payment_tenors: Vec<T>,
    pub onshore_forward_rates: Vec<T>,
    pub usd_ois_discount_factors: Vec<T>,
    pub local_rate_vol: T,
    pub ndf_fx_vol: T,
    pub rate_fx_corr: T,
    pub sovereign_bond_yield: T,
    pub usd_sofr_yield: T,
    pub ndf_ccs_basis_spread: T,
    pub withholding_tax_rate: T,
}
#[derive(Clone, Copy, Debug)]
pub struct AdbRbiNdirsAssetSwapReport<T: Real = f64> {
    pub onshore_par_irs_rate: T,
    pub offshore_par_ndirs_rate: T,
    pub ndirs_quanto_convexity_bps: T,
    pub after_tax_sovereign_yield: T,
    pub fx_hedged_usd_asset_swap_yield: T,
    pub usd_asset_swap_pickup_bps: T,
}
pub fn evaluate_adb_rbi_ndirs_and_asset_swap_basis<T: Real>(
    s: &AdbRbiNdirsAssetSwapSpec<T>,
) -> AdbRbiNdirsAssetSwapReport<T> {
    let mut w = T::zero();
    let mut on = T::zero();
    let mut off = T::zero();
    for i in 0..s
        .payment_tenors
        .len()
        .min(s.onshore_forward_rates.len())
        .min(s.usd_ois_discount_factors.len())
    {
        let dt = if i == 0 {
            s.payment_tenors[0]
        } else {
            s.payment_tenors[i] - s.payment_tenors[i - 1]
        };
        let wt = dt.max_val(T::from_f64(1e-8)) * s.usd_ois_discount_factors[i];
        let cv = (s.ndf_fx_vol * s.ndf_fx_vol - s.rate_fx_corr * s.local_rate_vol * s.ndf_fx_vol)
            * s.payment_tenors[i];
        w += wt;
        on += wt * s.onshore_forward_rates[i];
        off += wt * (s.onshore_forward_rates[i] + cv);
    }
    let onr = on / w.max_val(T::from_f64(1e-8));
    let offr = off / w.max_val(T::from_f64(1e-8));
    let pickup = s.sovereign_bond_yield * (T::one() - s.withholding_tax_rate) - s.usd_sofr_yield
        + s.ndf_ccs_basis_spread;
    AdbRbiNdirsAssetSwapReport {
        onshore_par_irs_rate: onr,
        offshore_par_ndirs_rate: offr,
        ndirs_quanto_convexity_bps: (offr - onr) * T::from_f64(1e4),
        after_tax_sovereign_yield: s.sovereign_bond_yield * (T::one() - s.withholding_tax_rate),
        fx_hedged_usd_asset_swap_yield: s.usd_sofr_yield + pickup,
        usd_asset_swap_pickup_bps: pickup * T::from_f64(1e4),
    }
}
#[cfg(test)]
mod phase32_tests {
    use super::*;
    #[test]
    fn convexity_zero_fx() {
        let s = AdbRbiNdirsAssetSwapSpec {
            payment_tenors: vec![1.0],
            onshore_forward_rates: vec![0.03],
            usd_ois_discount_factors: vec![1.0],
            local_rate_vol: 0.1,
            ndf_fx_vol: 0.0,
            rate_fx_corr: 0.0,
            sovereign_bond_yield: 0.04,
            usd_sofr_yield: 0.03,
            ndf_ccs_basis_spread: 0.0,
            withholding_tax_rate: 0.0,
        };
        assert_eq!(
            evaluate_adb_rbi_ndirs_and_asset_swap_basis(&s).ndirs_quanto_convexity_bps,
            0.0
        );
    }
}
