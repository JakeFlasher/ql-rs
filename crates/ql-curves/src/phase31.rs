use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct FabozziKeyRateLscImmunizationSpec<T: Real = f64> {
    pub cashflow_tenors: Vec<T>,
    pub cashflows: Vec<T>,
    pub zero_rates: Vec<T>,
    pub key_rate_tenors: Vec<T>,
    pub ns_decay_lambda: T,
    pub liability_horizon_years: T,
    pub barbell_short_tenor: T,
    pub barbell_long_tenor: T,
    pub twist_slope_shock_bps: T,
    pub butterfly_curvature_shock_bps: T,
}

#[derive(Clone, Debug)]
pub struct FabozziKeyRateLscImmunizationReport<T: Real = f64> {
    pub bond_present_value: T,
    pub fisher_weil_duration: T,
    pub key_rate_durations: Vec<T>,
    pub krd_partition_unity_error: T,
    pub lsc_level_duration: T,
    pub lsc_slope_duration: T,
    pub lsc_curvature_duration: T,
    pub barbell_short_weight: T,
    pub barbell_convexity_advantage: T,
    pub non_parallel_shock_pnl_pct: T,
}

pub fn evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization<T: Real>(
    s: &FabozziKeyRateLscImmunizationSpec<T>,
) -> FabozziKeyRateLscImmunizationReport<T> {
    let n = s
        .cashflow_tenors
        .len()
        .min(s.cashflows.len())
        .min(s.zero_rates.len());
    assert!(n > 0 && !s.key_rate_tenors.is_empty());
    let mut pv = T::zero();
    let mut dur_num = T::zero();
    let mut conv_num = T::zero();
    let mut lsc_slope = T::zero();
    let mut lsc_curv = T::zero();
    let lambda = s.ns_decay_lambda.max_val(T::from_f64(1e-12));
    let mut krd_num = vec![T::zero(); s.key_rate_tenors.len()];
    for i in 0..n {
        let tau = s.cashflow_tenors[i];
        let pv_cf = s.cashflows[i] * (-s.zero_rates[i] * tau).exp();
        pv += pv_cf;
        dur_num += tau * pv_cf;
        conv_num += tau * tau * pv_cf;
        let x = lambda * tau;
        let e = (-x).exp();
        let l2 = (T::one() - e) / x.max_val(T::from_f64(1e-12));
        lsc_slope += tau * l2 * pv_cf;
        lsc_curv += tau * (l2 - e) * pv_cf;
        for k in 0..s.key_rate_tenors.len() {
            let tk = s.key_rate_tenors[k];
            let w = if (k == 0 && tau <= tk) || (k + 1 == s.key_rate_tenors.len() && tau >= tk) {
                T::one()
            } else if k > 0 && tau >= s.key_rate_tenors[k - 1] && tau <= tk {
                (tau - s.key_rate_tenors[k - 1])
                    / (tk - s.key_rate_tenors[k - 1]).max_val(T::from_f64(1e-12))
            } else if k + 1 < s.key_rate_tenors.len()
                && tau >= tk
                && tau <= s.key_rate_tenors[k + 1]
            {
                (s.key_rate_tenors[k + 1] - tau)
                    / (s.key_rate_tenors[k + 1] - tk).max_val(T::from_f64(1e-12))
            } else {
                T::zero()
            };
            krd_num[k] += tau * w * pv_cf;
        }
    }
    let den = pv.abs().max_val(T::from_f64(1e-12));
    let duration = dur_num / den;
    let krd: Vec<T> = krd_num.into_iter().map(|v| v / den).collect();
    let sum_krd = krd.iter().fold(T::zero(), |a, &b| a + b);
    let short = s.barbell_short_tenor;
    let long = s.barbell_long_tenor;
    let horizon = s.liability_horizon_years;
    let ws = (long - horizon) / (long - short).max_val(T::from_f64(1e-12));
    let barbell_conv = ws * short * short + (T::one() - ws) * long * long;
    let bullet_conv = horizon * horizon;
    let non_parallel = -(lsc_slope / den * s.twist_slope_shock_bps * T::from_f64(1e-4)
        + lsc_curv / den * s.butterfly_curvature_shock_bps * T::from_f64(1e-4));
    FabozziKeyRateLscImmunizationReport {
        bond_present_value: pv,
        fisher_weil_duration: duration,
        key_rate_durations: krd,
        krd_partition_unity_error: (duration - sum_krd).abs(),
        lsc_level_duration: duration,
        lsc_slope_duration: lsc_slope / den,
        lsc_curvature_duration: lsc_curv / den,
        barbell_short_weight: ws,
        barbell_convexity_advantage: barbell_conv - bullet_conv,
        non_parallel_shock_pnl_pct: non_parallel,
    }
}

#[cfg(test)]
mod phase31_tests {
    use super::*;
    use ql_ad::Dual;

    fn spec<T: Real>(x: T) -> FabozziKeyRateLscImmunizationSpec<T> {
        FabozziKeyRateLscImmunizationSpec {
            cashflow_tenors: vec![T::one(), T::from_f64(5.0), T::from_f64(10.0)],
            cashflows: vec![T::from_f64(2.0) + x, T::from_f64(2.0), T::from_f64(102.0)],
            zero_rates: vec![T::from_f64(0.03), T::from_f64(0.035), T::from_f64(0.04)],
            key_rate_tenors: vec![T::one(), T::from_f64(5.0), T::from_f64(10.0)],
            ns_decay_lambda: T::from_f64(0.5),
            liability_horizon_years: T::from_f64(5.0),
            barbell_short_tenor: T::one(),
            barbell_long_tenor: T::from_f64(10.0),
            twist_slope_shock_bps: T::from_f64(20.0),
            butterfly_curvature_shock_bps: T::from_f64(10.0),
        }
    }

    #[test]
    fn invariants_and_dual_fd() {
        let r = evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization(&spec(0.0));
        assert!(r.krd_partition_unity_error < 1e-12);
        assert!(r.barbell_convexity_advantage > 0.0);
        let h = 1e-6;
        let fd = (evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization(&spec(h))
            .bond_present_value
            - evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization(&spec(-h)).bond_present_value)
            / (2.0 * h);
        let ad = evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization(&spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .bond_present_value;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
