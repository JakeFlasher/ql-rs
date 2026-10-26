//! De Spiegeleer, Schoutens & Van Hulle (2014) Convertible Bond Soft-Call Parity Trigger,
//! Make-Whole Additional Shares Table & Stock-Borrow Fee Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct DeSpiegeleerCbSoftCallMakeWholeSpec<T: Real = f64> {
    pub spot_price_s0: T,
    pub par_value: T,
    pub conversion_ratio_cr: T,
    pub soft_call_trigger_multiple: T,
    pub make_whole_max_bonus_shares: T,
    pub coupon_rate_c: T,
    pub risky_bond_yield_y: T,
    pub risk_free_rate_r: T,
    pub dividend_and_borrow_cost_q: T,
    pub equity_volatility_sigma: T,
    pub maturity_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct DeSpiegeleerCbSoftCallMakeWholeReport<T: Real = f64> {
    pub conversion_price_k_conv: T,
    pub soft_call_barrier_price_b: T,
    pub conversion_parity_value: T,
    pub straight_bond_floor: T,
    pub embedded_conversion_call_value: T,
    pub soft_call_issuer_cap_reduction: T,
    pub make_whole_table_protection_pv: T,
    pub convertible_bond_fair_value: T,
}

fn bs_call<T: Real>(s: T, k: T, r: T, q: T, sig: T, t: T) -> T {
    let std = sig * t.sqrt();
    let d1 = ((s / k).ln() + (r - q + T::from_f64(0.5) * sig * sig) * t) / std;
    let d2 = d1 - std;
    s * (-q * t).exp() * norm_cdf(d1) - k * (-r * t).exp() * norm_cdf(d2)
}

pub fn evaluate_despiegeleer_cb_soft_call_make_whole<T: Real>(
    spec: &DeSpiegeleerCbSoftCallMakeWholeSpec<T>,
) -> DeSpiegeleerCbSoftCallMakeWholeReport<T> {
    let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
    let par = spec.par_value.max_val(T::one());
    let cr = spec.conversion_ratio_cr.max_val(T::from_f64(1e-4));
    let trig_mult = spec.soft_call_trigger_multiple.max_val(T::from_f64(1.05));
    let t = spec.maturity_years_t.max_val(T::from_f64(0.25));
    let sig = spec.equity_volatility_sigma.max_val(T::from_f64(1e-4));
    let y = spec.risky_bond_yield_y.max_val(T::from_f64(1e-4));

    let k_conv = par / cr;
    let b_soft = k_conv * trig_mult;
    let parity = cr * s0;

    // Straight risky bond floor:
    let annuity = (T::one() - (-y * t).exp()) / y;
    let bond_floor = par * spec.coupon_rate_c * annuity + par * (-y * t).exp();

    // Embedded conversion call struck at effective floor strike K_eff = BondFloor / CR:
    let k_eff = (bond_floor / cr).max_val(T::from_f64(1e-4));
    let conv_call = cr
        * bs_call(
            s0,
            k_eff,
            spec.risk_free_rate_r,
            spec.dividend_and_borrow_cost_q,
            sig,
            t,
        );

    // Provisional soft-call cap above B_soft (weighted by M-of-N excursion damping):
    let soft_cap_raw = cr
        * bs_call(
            s0,
            b_soft,
            spec.risk_free_rate_r,
            spec.dividend_and_borrow_cost_q,
            sig,
            t,
        );
    let excursion_prob = (-T::from_f64(0.35) * sig * t.sqrt()).exp();
    let soft_call_reduction = T::from_f64(0.5) * excursion_prob * soft_cap_raw;

    // Make-whole bonus shares protection PV when called or acquired above K_conv:
    let mw_pv = spec.make_whole_max_bonus_shares
        * (bs_call(
            s0,
            k_conv,
            spec.risk_free_rate_r,
            spec.dividend_and_borrow_cost_q,
            sig,
            t,
        ) - bs_call(
            s0,
            b_soft,
            spec.risk_free_rate_r,
            spec.dividend_and_borrow_cost_q,
            sig,
            t,
        ))
        .max_val(T::zero());

    let cb_val = bond_floor + conv_call - soft_call_reduction + mw_pv;

    DeSpiegeleerCbSoftCallMakeWholeReport {
        conversion_price_k_conv: k_conv,
        soft_call_barrier_price_b: b_soft,
        conversion_parity_value: parity,
        straight_bond_floor: bond_floor,
        embedded_conversion_call_value: conv_call,
        soft_call_issuer_cap_reduction: soft_call_reduction,
        make_whole_table_protection_pv: mw_pv,
        convertible_bond_fair_value: cb_val,
    }
}

#[cfg(test)]
mod phase42_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn despiegeleer_cb_soft_call_make_whole_invariants_and_ad_parity() {
        let spec = DeSpiegeleerCbSoftCallMakeWholeSpec {
            spot_price_s0: 105.0,
            par_value: 100.0,
            conversion_ratio_cr: 1.0,
            soft_call_trigger_multiple: 1.30,
            make_whole_max_bonus_shares: 0.18,
            coupon_rate_c: 0.025,
            risky_bond_yield_y: 0.055,
            risk_free_rate_r: 0.040,
            dividend_and_borrow_cost_q: 0.018,
            equity_volatility_sigma: 0.28,
            maturity_years_t: 4.0,
        };
        let rep = evaluate_despiegeleer_cb_soft_call_make_whole(&spec);
        assert!(rep.convertible_bond_fair_value > rep.straight_bond_floor);
        assert!(rep.soft_call_barrier_price_b > rep.conversion_price_k_conv);
        assert!(rep.make_whole_table_protection_pv > 0.0);

        let h = 1e-6;
        let r_ad =
            evaluate_despiegeleer_cb_soft_call_make_whole(&DeSpiegeleerCbSoftCallMakeWholeSpec {
                spot_price_s0: Dual::<1>::variable(spec.spot_price_s0, 0),
                par_value: Dual::<1>::from_value(spec.par_value),
                conversion_ratio_cr: Dual::<1>::from_value(spec.conversion_ratio_cr),
                soft_call_trigger_multiple: Dual::<1>::from_value(spec.soft_call_trigger_multiple),
                make_whole_max_bonus_shares: Dual::<1>::from_value(
                    spec.make_whole_max_bonus_shares,
                ),
                coupon_rate_c: Dual::<1>::from_value(spec.coupon_rate_c),
                risky_bond_yield_y: Dual::<1>::from_value(spec.risky_bond_yield_y),
                risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
                dividend_and_borrow_cost_q: Dual::<1>::from_value(spec.dividend_and_borrow_cost_q),
                equity_volatility_sigma: Dual::<1>::from_value(spec.equity_volatility_sigma),
                maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
            });
        let mut sp = spec;
        sp.spot_price_s0 += h;
        let mut sm = spec;
        sm.spot_price_s0 -= h;
        let fd = (evaluate_despiegeleer_cb_soft_call_make_whole(&sp).convertible_bond_fair_value
            - evaluate_despiegeleer_cb_soft_call_make_whole(&sm).convertible_bond_fair_value)
            / (2.0 * h);
        assert!((r_ad.convertible_bond_fair_value.d[0] - fd).abs() < 1e-5);
    }
}
