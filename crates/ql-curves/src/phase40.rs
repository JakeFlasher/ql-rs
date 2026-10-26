//! Brigo & Mercurio (2006) & Andersen–Piterbarg (2010) IBOR-to-OIS Fallback Spread,
//! Compounded-in-Arrears vs. In-Advance Convexity & Seasoned Accrual Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct BrigoOisInArrearsFallbackSpec<T: Real = f64> {
    pub ois_df_start_p0t1: T,
    pub ois_df_end_p0t2: T,
    pub accrual_year_fraction_delta: T,
    pub short_rate_vol_sigma_r: T,
    pub mean_reversion_a: T,
    pub isda_historical_fallback_spread: T,
    pub realized_accrued_rfr_rate: T,
    pub elapsed_accrual_fraction_w: T,
}

#[derive(Debug, Clone, Copy)]
pub struct BrigoOisInArrearsFallbackReport<T: Real = f64> {
    pub forward_compounded_in_arrears_rate: T,
    pub in_advance_timing_convexity_adj_bps: T,
    pub convexity_adjusted_in_advance_rate: T,
    pub isda_fallback_ibor_all_in_rate: T,
    pub seasoned_period_blended_rfr_rate: T,
    pub martingale_arrears_parity_error: T,
}

pub fn evaluate_brigo_ois_compounded_in_arrears_and_fallback<T: Real>(
    spec: &BrigoOisInArrearsFallbackSpec<T>,
) -> BrigoOisInArrearsFallbackReport<T> {
    let delta = spec.accrual_year_fraction_delta.max_val(T::from_f64(1e-4));
    let p1 = spec.ois_df_start_p0t1.max_val(T::from_f64(1e-6));
    let p2 = spec.ois_df_end_p0t2.max_val(T::from_f64(1e-6));

    // Under T_2-forward measure, daily compounded-in-arrears RFR paid at T_2 has
    // exact martingale expectation F_0(T_1, T_2) = (P(0, T_1) / P(0, T_2) - 1) / delta:
    let ratio = p1 / p2;
    let f_arrears = (ratio - T::one()) / delta;
    let parity_err = ((T::one() + delta * f_arrears) * p2 - p1).abs();

    // When a term RFR rate is paid in-advance at T_1 (or fixed in-arrears without compounding),
    // Hull-White timing convexity adjustment over [T_1, T_2] is:
    // var_int = sigma_r^2 * delta^3 / 3 * (1 - 0.5 * a * delta)
    let a_damp =
        (T::one() - T::from_f64(0.5) * spec.mean_reversion_a * delta).max_val(T::from_f64(0.2));
    let var_int = spec.short_rate_vol_sigma_r * spec.short_rate_vol_sigma_r * delta * delta * delta
        / T::from_f64(3.0)
        * a_damp;
    let f_in_advance = (ratio * (-var_int).exp() - T::one()) / delta;
    let conv_adj_bps = (f_arrears - f_in_advance) * T::from_f64(1e4);

    // ISDA IBOR Fallback all-in rate = Compounded-in-arrears RFR + 5Y median spread adjustment:
    let fallback_ibor = f_arrears + spec.isda_historical_fallback_spread;

    // Seasoned mid-period valuation with elapsed weight w in [0, 1]:
    // (1 + delta * R_blended) = (1 + w * delta * R_realized) * (1 + (1 - w) * delta * F_arrears)
    let w = spec
        .elapsed_accrual_fraction_w
        .max_val(T::zero())
        .min_val(T::one());
    let gross_realized = T::one() + w * delta * spec.realized_accrued_rfr_rate;
    let gross_remaining = T::one() + (T::one() - w) * delta * f_arrears;
    let blended_rate = (gross_realized * gross_remaining - T::one()) / delta;

    BrigoOisInArrearsFallbackReport {
        forward_compounded_in_arrears_rate: f_arrears,
        in_advance_timing_convexity_adj_bps: conv_adj_bps,
        convexity_adjusted_in_advance_rate: f_in_advance,
        isda_fallback_ibor_all_in_rate: fallback_ibor,
        seasoned_period_blended_rfr_rate: blended_rate,
        martingale_arrears_parity_error: parity_err,
    }
}

#[cfg(test)]
mod phase40_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn brigo_ois_in_arrears_fallback_invariants_and_ad_parity() {
        let spec = BrigoOisInArrearsFallbackSpec {
            ois_df_start_p0t1: 0.960,
            ois_df_end_p0t2: 0.938,
            accrual_year_fraction_delta: 0.5,
            short_rate_vol_sigma_r: 0.012,
            mean_reversion_a: 0.15,
            isda_historical_fallback_spread: 0.0026161,
            realized_accrued_rfr_rate: 0.049,
            elapsed_accrual_fraction_w: 0.4,
        };
        let rep = evaluate_brigo_ois_compounded_in_arrears_and_fallback(&spec);
        assert!(rep.martingale_arrears_parity_error < 1e-12);
        assert!(rep.in_advance_timing_convexity_adj_bps > 0.0);
        assert!(rep.isda_fallback_ibor_all_in_rate > rep.forward_compounded_in_arrears_rate);

        let h = 1e-6;
        let r_ad =
            evaluate_brigo_ois_compounded_in_arrears_and_fallback(&BrigoOisInArrearsFallbackSpec {
                ois_df_start_p0t1: Dual::<1>::from_value(spec.ois_df_start_p0t1),
                ois_df_end_p0t2: Dual::<1>::variable(spec.ois_df_end_p0t2, 0),
                accrual_year_fraction_delta: Dual::<1>::from_value(
                    spec.accrual_year_fraction_delta,
                ),
                short_rate_vol_sigma_r: Dual::<1>::from_value(spec.short_rate_vol_sigma_r),
                mean_reversion_a: Dual::<1>::from_value(spec.mean_reversion_a),
                isda_historical_fallback_spread: Dual::<1>::from_value(
                    spec.isda_historical_fallback_spread,
                ),
                realized_accrued_rfr_rate: Dual::<1>::from_value(spec.realized_accrued_rfr_rate),
                elapsed_accrual_fraction_w: Dual::<1>::from_value(spec.elapsed_accrual_fraction_w),
            });
        let mut sp = spec;
        sp.ois_df_end_p0t2 += h;
        let mut sm = spec;
        sm.ois_df_end_p0t2 -= h;
        let fd = (evaluate_brigo_ois_compounded_in_arrears_and_fallback(&sp)
            .seasoned_period_blended_rfr_rate
            - evaluate_brigo_ois_compounded_in_arrears_and_fallback(&sm)
                .seasoned_period_blended_rfr_rate)
            / (2.0 * h);
        assert!((r_ad.seasoned_period_blended_rfr_rate.d[0] - fd).abs() < 1e-7);
    }
}
