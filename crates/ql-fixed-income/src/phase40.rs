//! Fabozzi (2021) & Tuckman (2011) Collared Floating-Rate Note (`FRN`) Discount Margin (`DM`),
//! Spread Duration vs. Rate Duration & Option-Adjusted Discount Margin (`OADM`) Engine.

use ql_ad::Real;
use ql_math::{norm_cdf, norm_pdf};

#[derive(Debug, Clone, Copy)]
pub struct FabozziCollaredFrnDiscountMarginSpec<T: Real = f64> {
    pub par_face: T,
    pub index_rate_r: T,
    pub quoted_margin_qm: T,
    pub discount_margin_dm: T,
    pub cap_strike: T,
    pub floor_strike: T,
    pub rate_normal_vol_sigma: T,
    pub time_to_next_reset_tau1: T,
}

#[derive(Debug, Clone, Copy)]
pub struct FabozziCollaredFrnDiscountMarginReport<T: Real = f64> {
    pub straight_frn_price: T,
    pub par_reset_identity_error: T,
    pub embedded_cap_pv: T,
    pub embedded_floor_pv: T,
    pub collared_frn_price: T,
    pub spread_duration_dm: T,
    pub interest_rate_duration: T,
    pub option_adjusted_discount_margin_oadm: T,
}

fn bachelier_caplet<T: Real>(fwd: T, strike: T, std: T, is_cap: bool) -> T {
    let s = std.max_val(T::from_f64(1e-8));
    let diff = if is_cap { fwd - strike } else { strike - fwd };
    let d = diff / s;
    diff * norm_cdf(d) + s * norm_pdf(d)
}

fn straight_frn_pv<T: Real>(par: T, r: T, qm: T, dm: T, tau1: T) -> T {
    let n = 4usize;
    let mut pv = T::zero();
    let c_rate = r + qm;
    let disc_rate = (r + dm).max_val(T::from_f64(1e-4));
    for k in 1..=n {
        let tk = tau1 + T::from_f64((k - 1) as f64);
        let df = (T::one() + disc_rate).pow(-tk);
        let accrual = if k == 1 { tau1 } else { T::one() };
        let cf = if k == n {
            par * (T::one() + c_rate * accrual)
        } else {
            par * c_rate * accrual
        };
        pv += cf * df;
    }
    pv
}

pub fn evaluate_fabozzi_collared_frn_discount_margin<T: Real>(
    spec: &FabozziCollaredFrnDiscountMarginSpec<T>,
) -> FabozziCollaredFrnDiscountMarginReport<T> {
    let par = spec.par_face;
    let r = spec.index_rate_r;
    let qm = spec.quoted_margin_qm;
    let dm = spec.discount_margin_dm;
    let tau1 = spec.time_to_next_reset_tau1.max_val(T::from_f64(0.01));

    let v_straight = straight_frn_pv(par, r, qm, dm, tau1);
    // Exact Fabozzi/Tuckman FRN par-reset identity: at reset (tau1 = 1) with DM == QM, V == Par!
    let v_par_check = straight_frn_pv(par, r, qm, qm, T::one());
    let par_err = (v_par_check - par).abs();

    let fwd_coupon = r + qm;
    let disc_rate = (r + dm).max_val(T::from_f64(1e-4));
    let mut cap_pv = T::zero();
    let mut floor_pv = T::zero();
    let mut risky_dv01 = T::zero();

    for k in 1..=4usize {
        let tk = tau1 + T::from_f64((k - 1) as f64);
        let df = (T::one() + disc_rate).pow(-tk);
        let std = spec.rate_normal_vol_sigma * tk.sqrt();
        cap_pv += par * df * bachelier_caplet(fwd_coupon, spec.cap_strike, std, true);
        floor_pv += par * df * bachelier_caplet(fwd_coupon, spec.floor_strike, std, false);
        risky_dv01 += par * df;
    }

    // Collared FRN = Straight FRN + Long Floor - Short Cap:
    let v_collared = v_straight + floor_pv - cap_pv;

    // Analytical spread duration (w.r.t. DM) vs. interest-rate duration (w.r.t. parallel index shift):
    let bump = T::from_f64(1e-4);
    let sd_dm = (straight_frn_pv(par, r, qm, dm - bump, tau1)
        - straight_frn_pv(par, r, qm, dm + bump, tau1))
        / (T::from_f64(2.0) * bump * v_straight);
    let ir_dur = (straight_frn_pv(par, r - bump, qm, dm, tau1)
        - straight_frn_pv(par, r + bump, qm, dm, tau1))
        / (T::from_f64(2.0) * bump * v_straight);

    // Option-Adjusted Discount Margin (OADM):
    let oadm = dm + (floor_pv - cap_pv) / risky_dv01.max_val(T::from_f64(1e-6));

    FabozziCollaredFrnDiscountMarginReport {
        straight_frn_price: v_straight,
        par_reset_identity_error: par_err,
        embedded_cap_pv: cap_pv,
        embedded_floor_pv: floor_pv,
        collared_frn_price: v_collared,
        spread_duration_dm: sd_dm,
        interest_rate_duration: ir_dur.abs(),
        option_adjusted_discount_margin_oadm: oadm,
    }
}

#[cfg(test)]
mod phase40_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn fabozzi_collared_frn_discount_margin_invariants_and_ad_parity() {
        let spec = FabozziCollaredFrnDiscountMarginSpec {
            par_face: 100.0,
            index_rate_r: 0.042,
            quoted_margin_qm: 0.012,
            discount_margin_dm: 0.015,
            cap_strike: 0.065,
            floor_strike: 0.035,
            rate_normal_vol_sigma: 0.011,
            time_to_next_reset_tau1: 0.25,
        };
        let rep = evaluate_fabozzi_collared_frn_discount_margin(&spec);
        // Par-reset theorem & DM > QM => discount price < 100:
        assert!(rep.par_reset_identity_error < 1e-11);
        assert!(rep.straight_frn_price < 100.0);
        // Spread duration of an FRN is much larger than its interest-rate duration:
        assert!(rep.spread_duration_dm > rep.interest_rate_duration * 4.0);

        let h = 1e-6;
        let r_ad =
            evaluate_fabozzi_collared_frn_discount_margin(&FabozziCollaredFrnDiscountMarginSpec {
                par_face: Dual::<1>::from_value(spec.par_face),
                index_rate_r: Dual::<1>::from_value(spec.index_rate_r),
                quoted_margin_qm: Dual::<1>::from_value(spec.quoted_margin_qm),
                discount_margin_dm: Dual::<1>::variable(spec.discount_margin_dm, 0),
                cap_strike: Dual::<1>::from_value(spec.cap_strike),
                floor_strike: Dual::<1>::from_value(spec.floor_strike),
                rate_normal_vol_sigma: Dual::<1>::from_value(spec.rate_normal_vol_sigma),
                time_to_next_reset_tau1: Dual::<1>::from_value(spec.time_to_next_reset_tau1),
            });
        let mut sp = spec;
        sp.discount_margin_dm += h;
        let mut sm = spec;
        sm.discount_margin_dm -= h;
        let fd = (evaluate_fabozzi_collared_frn_discount_margin(&sp).collared_frn_price
            - evaluate_fabozzi_collared_frn_discount_margin(&sm).collared_frn_price)
            / (2.0 * h);
        assert!((r_ad.collared_frn_price.d[0] - fd).abs() < 1e-5);
    }
}
