//! Andersen & Piterbarg (2010) & Brigo–Mercurio (2006) Constant Maturity Swap (`CMS`)
//! Cash-vs-Physical Annuity Ratio & Swaption-Replication Convexity Adjustment Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct AndersenPiterbargCmsReplicationSpec<T: Real = f64> {
    pub forward_swap_rate_s0: T,
    pub swap_tenor_years_n: T,
    pub payment_frequency_m: T,
    pub fixing_expiry_te: T,
    pub normal_swaption_vol_sigma_n: T,
    pub ois_spread_delta_r: T,
    pub mean_reversion_speed_a: T,
}

#[derive(Debug, Clone, Copy)]
pub struct AndersenPiterbargCmsReplicationReport<T: Real = f64> {
    pub cash_settled_annuity_gc: T,
    pub physical_ois_annuity_ap: T,
    pub annuity_derivative_ratio_g1: T,
    pub annuity_second_deriv_ratio_g2: T,
    pub replication_cms_convexity_adj_bps: T,
    pub convexity_adjusted_cms_rate: T,
    pub in_arrears_swap_rate_convexity_bps: T,
}

fn cash_annuity<T: Real>(s: T, n_years: T, m_freq: T) -> T {
    let s_pos = s.max_val(T::from_f64(1e-4));
    let m = m_freq.max_val(T::one());
    let total_periods = m * n_years.max_val(T::one());
    let base = T::one() + s_pos / m;
    (T::one() - base.pow(-total_periods)) / s_pos
}

pub fn evaluate_andersen_piterbarg_cms_replication_convexity<T: Real>(
    spec: &AndersenPiterbargCmsReplicationSpec<T>,
) -> AndersenPiterbargCmsReplicationReport<T> {
    let s0 = spec.forward_swap_rate_s0.max_val(T::from_f64(1e-4));
    let n = spec.swap_tenor_years_n.max_val(T::one());
    let m = spec.payment_frequency_m.max_val(T::one());
    let te = spec.fixing_expiry_te.max_val(T::from_f64(0.05));
    let sig_n = spec.normal_swaption_vol_sigma_n.max_val(T::from_f64(1e-5));

    let gc = cash_annuity(s0, n, m);
    let s_ois = (s0 - spec.ois_spread_delta_r).max_val(T::from_f64(1e-4));
    let ap = cash_annuity(s_ois, n, m);

    // Analytic finite-step derivatives of cash annuity G(S) in T: Real:
    let ds = T::from_f64(1e-4);
    let g_up = cash_annuity(s0 + ds, n, m);
    let g_dn = cash_annuity(s0 - ds, n, m);
    let g_prime = (g_up - g_dn) / (T::from_f64(2.0) * ds);
    let g_second = (g_up - T::from_f64(2.0) * gc + g_dn) / (ds * ds);

    // Andersen-Piterbarg (2010, Vol. I Ch. 16) annuity-slope ratio -G'(S_0)/G(S_0) > 0:
    let g1_ratio = (-g_prime / gc).max_val(T::from_f64(1e-4));
    let g2_ratio = (g_second / gc).max_val(T::from_f64(1e-4));

    // Mean-reversion damping of swap-rate variance over [0, T_e]:
    let a = spec.mean_reversion_speed_a.max_val(T::from_f64(1e-3));
    let mr_factor = (T::one() - (T::from_f64(-2.0) * a * te).exp()) / (T::from_f64(2.0) * a);
    let var_s = sig_n * sig_n * mr_factor * (ap / gc);

    // CMS convexity adjustment = (-G'(S_0) / G(S_0)) * Var(S_{T_e}):
    let cms_conv = g1_ratio * var_s;
    let cms_rate = s0 + cms_conv;
    let in_arrears_bps = (cms_conv + T::from_f64(0.5) * g2_ratio * s0 * var_s) * T::from_f64(1e4);

    AndersenPiterbargCmsReplicationReport {
        cash_settled_annuity_gc: gc,
        physical_ois_annuity_ap: ap,
        annuity_derivative_ratio_g1: g1_ratio,
        annuity_second_deriv_ratio_g2: g2_ratio,
        replication_cms_convexity_adj_bps: cms_conv * T::from_f64(1e4),
        convexity_adjusted_cms_rate: cms_rate,
        in_arrears_swap_rate_convexity_bps: in_arrears_bps,
    }
}

#[cfg(test)]
mod phase41_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn andersen_piterbarg_cms_replication_invariants_and_ad_parity() {
        let spec = AndersenPiterbargCmsReplicationSpec {
            forward_swap_rate_s0: 0.038,
            swap_tenor_years_n: 10.0,
            payment_frequency_m: 2.0,
            fixing_expiry_te: 5.0,
            normal_swaption_vol_sigma_n: 0.0085,
            ois_spread_delta_r: 0.0015,
            mean_reversion_speed_a: 0.04,
        };
        let rep = evaluate_andersen_piterbarg_cms_replication_convexity(&spec);
        assert!(rep.cash_settled_annuity_gc > 0.0);
        assert!(rep.replication_cms_convexity_adj_bps > 0.0);
        assert!(rep.convexity_adjusted_cms_rate > spec.forward_swap_rate_s0);

        let h = 1e-6;
        let r_ad = evaluate_andersen_piterbarg_cms_replication_convexity(
            &AndersenPiterbargCmsReplicationSpec {
                forward_swap_rate_s0: Dual::<1>::variable(spec.forward_swap_rate_s0, 0),
                swap_tenor_years_n: Dual::<1>::from_value(spec.swap_tenor_years_n),
                payment_frequency_m: Dual::<1>::from_value(spec.payment_frequency_m),
                fixing_expiry_te: Dual::<1>::from_value(spec.fixing_expiry_te),
                normal_swaption_vol_sigma_n: Dual::<1>::from_value(
                    spec.normal_swaption_vol_sigma_n,
                ),
                ois_spread_delta_r: Dual::<1>::from_value(spec.ois_spread_delta_r),
                mean_reversion_speed_a: Dual::<1>::from_value(spec.mean_reversion_speed_a),
            },
        );
        let mut sp = spec;
        sp.forward_swap_rate_s0 += h;
        let mut sm = spec;
        sm.forward_swap_rate_s0 -= h;
        let fd = (evaluate_andersen_piterbarg_cms_replication_convexity(&sp)
            .convexity_adjusted_cms_rate
            - evaluate_andersen_piterbarg_cms_replication_convexity(&sm)
                .convexity_adjusted_cms_rate)
            / (2.0 * h);
        assert!((r_ad.convexity_adjusted_cms_rate.d[0] - fd).abs() < 1e-6);
    }
}
