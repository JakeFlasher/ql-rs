//! O'Kane (2008) & Brigo–Morini (2008) Credit Index Intrinsic vs. Market Skew (`CDX`/`iTraxx`),
//! Constituent Dispersion Convexity & Index Arbitrage Basis Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct OKaneCdsIndexIntrinsicSkewSpec<T: Real = f64> {
    pub constituent_spreads_s: [T; 3],
    pub market_quoted_index_spread: T,
    pub standard_index_coupon_c: T,
    pub recovery_rate_r: T,
    pub ois_discount_rate_r: T,
    pub maturity_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct OKaneCdsIndexIntrinsicSkewReport<T: Real = f64> {
    pub equal_weight_average_spread: T,
    pub rpv01_weighted_intrinsic_spread: T,
    pub constituent_dispersion_wedge_bps: T,
    pub market_minus_intrinsic_skew_bps: T,
    pub intrinsic_portfolio_upfront_pct: T,
    pub flat_curve_index_upfront_pct: T,
    pub index_vs_constituents_basis_pv_bps: T,
}

fn rpv01_flat<T: Real>(r: T, lam: T, t: T) -> T {
    let denom = (r + lam).max_val(T::from_f64(1e-6));
    (T::one() - (-(r + lam) * t).exp()) / denom
}

pub fn evaluate_okane_cds_index_intrinsic_skew_basis<T: Real>(
    spec: &OKaneCdsIndexIntrinsicSkewSpec<T>,
) -> OKaneCdsIndexIntrinsicSkewReport<T> {
    let lgd = (T::one() - spec.recovery_rate_r).max_val(T::from_f64(0.05));
    let r = spec.ois_discount_rate_r;
    let t = spec.maturity_years_t.max_val(T::from_f64(0.25));
    let c = spec.standard_index_coupon_c;
    let third = T::from_f64(1.0 / 3.0);

    let mut avg_s = T::zero();
    let mut num_s = T::zero();
    let mut den_rpv = T::zero();
    let mut intr_upfront = T::zero();

    for i in 0..3 {
        let si = spec.constituent_spreads_s[i].max_val(T::from_f64(1e-5));
        let lami = si / lgd;
        let rpvi = rpv01_flat(r, lami, t);
        avg_s += third * si;
        num_s += third * si * rpvi;
        den_rpv += third * rpvi;
        intr_upfront += third * (si - c) * rpvi;
    }

    // Because RPV01(S_i) is decreasing and convex in S_i, RPV01-weighted intrinsic spread
    // is strictly below the simple average spread when constituents have dispersion:
    let s_intr = num_s / den_rpv.max_val(T::from_f64(1e-8));
    let disp_wedge_bps = (avg_s - s_intr) * T::from_f64(1e4);
    let mkt_skew_bps = (spec.market_quoted_index_spread - s_intr) * T::from_f64(1e4);

    let lam_mkt = spec.market_quoted_index_spread.max_val(T::from_f64(1e-5)) / lgd;
    let rpv_mkt = rpv01_flat(r, lam_mkt, t);
    let mkt_upfront = (spec.market_quoted_index_spread - c) * rpv_mkt;
    let basis_pv_bps = (mkt_upfront - intr_upfront) * T::from_f64(1e4);

    OKaneCdsIndexIntrinsicSkewReport {
        equal_weight_average_spread: avg_s,
        rpv01_weighted_intrinsic_spread: s_intr,
        constituent_dispersion_wedge_bps: disp_wedge_bps,
        market_minus_intrinsic_skew_bps: mkt_skew_bps,
        intrinsic_portfolio_upfront_pct: intr_upfront * T::from_f64(100.0),
        flat_curve_index_upfront_pct: mkt_upfront * T::from_f64(100.0),
        index_vs_constituents_basis_pv_bps: basis_pv_bps,
    }
}

#[cfg(test)]
mod phase42_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn okane_cds_index_intrinsic_skew_basis_invariants_and_ad_parity() {
        let spec = OKaneCdsIndexIntrinsicSkewSpec {
            constituent_spreads_s: [0.0060, 0.0120, 0.0300],
            market_quoted_index_spread: 0.0155,
            standard_index_coupon_c: 0.0100,
            recovery_rate_r: 0.40,
            ois_discount_rate_r: 0.038,
            maturity_years_t: 5.0,
        };
        let rep = evaluate_okane_cds_index_intrinsic_skew_basis(&spec);
        // Dispersion convexity theorem: S_intrinsic < S_average for heterogeneous spreads:
        assert!(rep.rpv01_weighted_intrinsic_spread < rep.equal_weight_average_spread);
        assert!(rep.constituent_dispersion_wedge_bps > 0.0);

        let h = 1e-6;
        let r_ad = evaluate_okane_cds_index_intrinsic_skew_basis(&OKaneCdsIndexIntrinsicSkewSpec {
            constituent_spreads_s: [
                Dual::<1>::variable(spec.constituent_spreads_s[0], 0),
                Dual::<1>::from_value(spec.constituent_spreads_s[1]),
                Dual::<1>::from_value(spec.constituent_spreads_s[2]),
            ],
            market_quoted_index_spread: Dual::<1>::from_value(spec.market_quoted_index_spread),
            standard_index_coupon_c: Dual::<1>::from_value(spec.standard_index_coupon_c),
            recovery_rate_r: Dual::<1>::from_value(spec.recovery_rate_r),
            ois_discount_rate_r: Dual::<1>::from_value(spec.ois_discount_rate_r),
            maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
        });
        let mut sp = spec;
        sp.constituent_spreads_s[0] += h;
        let mut sm = spec;
        sm.constituent_spreads_s[0] -= h;
        let fd = (evaluate_okane_cds_index_intrinsic_skew_basis(&sp)
            .rpv01_weighted_intrinsic_spread
            - evaluate_okane_cds_index_intrinsic_skew_basis(&sm).rpv01_weighted_intrinsic_spread)
            / (2.0 * h);
        assert!((r_ad.rpv01_weighted_intrinsic_spread.d[0] - fd).abs() < 1e-7);
    }
}
