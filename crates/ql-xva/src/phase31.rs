use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct GregoryCsaThresholdMtaCtdColvaSpec<T: Real = f64> {
    pub time_grid: Vec<T>,
    pub expected_mtm_profile: Vec<T>,
    pub exposure_vol: T,
    pub mpor_years: T,
    pub csa_threshold: T,
    pub minimum_transfer_amount: T,
    pub independent_amount: T,
    pub hazard_rate: T,
    pub lgd: T,
    pub risk_free_rate: T,
    pub ccy1_collateral_spread: T,
    pub ccy2_collateral_spread: T,
    pub collateral_spread_vol: T,
}

#[derive(Clone, Copy, Debug)]
pub struct GregoryCsaThresholdMtaCtdColvaReport<T: Real = f64> {
    pub effective_csa_threshold: T,
    pub uncollateralized_cva: T,
    pub threshold_mta_collateralized_cva: T,
    pub csa_cva_reduction_benefit: T,
    pub single_currency_colva: T,
    pub ctd_multi_currency_colva: T,
    pub ctd_collateral_switch_option_value: T,
}

pub fn evaluate_gregory_csa_threshold_mta_and_ctd_colva<T: Real>(
    s: &GregoryCsaThresholdMtaCtdColvaSpec<T>,
) -> GregoryCsaThresholdMtaCtdColvaReport<T> {
    let threshold = (s.csa_threshold + s.minimum_transfer_amount / T::from_f64(2.0)
        - s.independent_amount)
        .max_val(T::zero());
    let sigma = (s.exposure_vol * s.mpor_years.sqrt()).max_val(T::from_f64(1e-10));
    let mut ucva = T::zero();
    let mut ccva = T::zero();
    let mut single = T::zero();
    let mut ctd = T::zero();
    for i in 0..s.time_grid.len().min(s.expected_mtm_profile.len()) {
        let t = s.time_grid[i];
        let mu = s.expected_mtm_profile[i];
        let z0 = mu / sigma;
        let zt = (mu - threshold) / sigma;
        let epe0 = mu * z0.norm_cdf() + sigma * z0.norm_pdf();
        let epet = (mu - threshold) * zt.norm_cdf() + sigma * zt.norm_pdf();
        let disc = (-s.risk_free_rate * t).exp();
        let dp = s.hazard_rate * (-s.hazard_rate * t).exp() * disc * s.lgd;
        ucva += dp * epe0;
        ccva += dp * epet;
        let ds = s.collateral_spread_vol * t.sqrt();
        let ds_safe = ds.max_val(T::from_f64(1e-10));
        let zeta = (s.ccy1_collateral_spread - s.ccy2_collateral_spread) / ds_safe;
        let max_spread = s.ccy1_collateral_spread * zeta.norm_cdf()
            + s.ccy2_collateral_spread * (-zeta).norm_cdf()
            + ds * zeta.norm_pdf();
        single += disc * s.ccy1_collateral_spread;
        ctd += disc * max_spread;
    }
    let option = (ctd - single).max_val(T::zero());
    GregoryCsaThresholdMtaCtdColvaReport {
        effective_csa_threshold: threshold,
        uncollateralized_cva: ucva,
        threshold_mta_collateralized_cva: ccva,
        csa_cva_reduction_benefit: ucva - ccva,
        single_currency_colva: single,
        ctd_multi_currency_colva: ctd,
        ctd_collateral_switch_option_value: option,
    }
}

#[cfg(test)]
mod phase31_tests {
    use super::*;
    use ql_ad::Dual;

    fn spec<T: Real>(x: T) -> GregoryCsaThresholdMtaCtdColvaSpec<T> {
        GregoryCsaThresholdMtaCtdColvaSpec {
            time_grid: vec![T::one(), T::from_f64(2.0), T::from_f64(3.0)],
            expected_mtm_profile: vec![T::from_f64(5.0) + x, T::from_f64(8.0), T::from_f64(10.0)],
            exposure_vol: T::from_f64(3.0),
            mpor_years: T::from_f64(10.0 / 252.0),
            csa_threshold: T::from_f64(2.0),
            minimum_transfer_amount: T::one(),
            independent_amount: T::from_f64(0.5),
            hazard_rate: T::from_f64(0.02),
            lgd: T::from_f64(0.6),
            risk_free_rate: T::from_f64(0.03),
            ccy1_collateral_spread: T::from_f64(0.01),
            ccy2_collateral_spread: T::from_f64(0.008),
            collateral_spread_vol: T::from_f64(0.01),
        }
    }

    #[test]
    fn phase31_invariants_and_dual_fd() {
        let r = evaluate_gregory_csa_threshold_mta_and_ctd_colva(&spec(0.0));
        assert!(r.threshold_mta_collateralized_cva <= r.uncollateralized_cva);
        assert!(r.ctd_multi_currency_colva >= r.single_currency_colva);
        assert!(r.ctd_collateral_switch_option_value > 0.0);
        let h = 1e-6;
        let fd = (evaluate_gregory_csa_threshold_mta_and_ctd_colva(&spec(h)).uncollateralized_cva
            - evaluate_gregory_csa_threshold_mta_and_ctd_colva(&spec(-h)).uncollateralized_cva)
            / (2.0 * h);
        let ad =
            evaluate_gregory_csa_threshold_mta_and_ctd_colva(&spec(Dual::<1>::variable(0.0, 0)))
                .uncollateralized_cva;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
