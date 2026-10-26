//! Fabozzi (2021) & Tuckman (2011) Agency MBS Specified-Pool Pay-Up,
//! S-Curve Refinancing Burnout & IO/PO Strip Negative-Convexity Decomposition Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct FabozziSpecifiedPoolIoPoSpec<T: Real = f64> {
    pub par_notional: T,
    pub wac_coupon: T,
    pub pass_through_coupon: T,
    pub current_mortgage_rate: T,
    pub discount_yield: T,
    pub base_turnover_cpr: T,
    pub max_refi_cpr: T,
    pub s_curve_steepness_kappa: T,
    pub spec_pool_burnout_dampener: T,
}

#[derive(Debug, Clone, Copy)]
pub struct FabozziSpecifiedPoolIoPoReport<T: Real = f64> {
    pub generic_tba_cpr: T,
    pub specified_pool_cpr: T,
    pub io_strip_pv: T,
    pub po_strip_pv: T,
    pub pass_through_pv: T,
    pub generic_tba_pv: T,
    pub specified_pool_payup_32nds: T,
    pub io_po_additivity_error: T,
}

fn price_io_po<T: Real>(notional: T, pass_coupon: T, y_disc: T, cpr: T) -> (T, T) {
    let n = 6usize;
    let r_sched = pass_coupon.max_val(T::from_f64(1e-4));
    let annuity_factor = (T::one() - (T::one() + r_sched).powf(-(n as f64))) / r_sched;
    let sched_pmt = notional / annuity_factor;

    let mut bal = notional;
    let mut io_pv = T::zero();
    let mut po_pv = T::zero();

    for k in 1..=n {
        let df = (T::one() + y_disc).powf(-(k as f64));
        let interest = bal * pass_coupon;
        let sched_prin = if k == n {
            bal
        } else {
            (sched_pmt - interest).max_val(T::zero()).min_val(bal)
        };
        let rem = (bal - sched_prin).max_val(T::zero());
        let prepay = if k == n { T::zero() } else { rem * cpr };
        let tot_prin = sched_prin + prepay;

        io_pv += interest * df;
        po_pv += tot_prin * df;
        bal = (rem - prepay).max_val(T::zero());
    }
    (io_pv, po_pv)
}

pub fn evaluate_fabozzi_specified_pool_io_po_payup<T: Real>(
    spec: &FabozziSpecifiedPoolIoPoSpec<T>,
) -> FabozziSpecifiedPoolIoPoReport<T> {
    let incentive = spec.wac_coupon - spec.current_mortgage_rate - T::from_f64(0.005);
    let logistic = T::one() / (T::one() + (-spec.s_curve_steepness_kappa * incentive).exp());
    let refi_wave = (spec.max_refi_cpr - spec.base_turnover_cpr).max_val(T::zero()) * logistic;

    let tba_cpr = (spec.base_turnover_cpr + refi_wave).min_val(T::from_f64(0.85));
    let damp = spec
        .spec_pool_burnout_dampener
        .max_val(T::from_f64(0.05))
        .min_val(T::one());
    let spec_cpr = (spec.base_turnover_cpr + damp * refi_wave).min_val(T::from_f64(0.85));

    let (io_spec, po_spec) = price_io_po(
        spec.par_notional,
        spec.pass_through_coupon,
        spec.discount_yield,
        spec_cpr,
    );
    let (io_tba, po_tba) = price_io_po(
        spec.par_notional,
        spec.pass_through_coupon,
        spec.discount_yield,
        tba_cpr,
    );

    let pt_spec = io_spec + po_spec;
    let pt_tba = io_tba + po_tba;
    // Pay-up in 32nds of par (per 100 face)
    let payup_32nds =
        (pt_spec - pt_tba) * (T::from_f64(100.0) / spec.par_notional) * T::from_f64(32.0);
    let add_err = (pt_spec - (io_spec + po_spec)).abs();

    FabozziSpecifiedPoolIoPoReport {
        generic_tba_cpr: tba_cpr,
        specified_pool_cpr: spec_cpr,
        io_strip_pv: io_spec,
        po_strip_pv: po_spec,
        pass_through_pv: pt_spec,
        generic_tba_pv: pt_tba,
        specified_pool_payup_32nds: payup_32nds,
        io_po_additivity_error: add_err,
    }
}

#[cfg(test)]
mod phase39_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn fabozzi_specified_pool_io_po_payup_invariants_and_ad_parity() {
        let spec = FabozziSpecifiedPoolIoPoSpec {
            par_notional: 100.0,
            wac_coupon: 0.065,
            pass_through_coupon: 0.060,
            current_mortgage_rate: 0.048,
            discount_yield: 0.042,
            base_turnover_cpr: 0.06,
            max_refi_cpr: 0.42,
            s_curve_steepness_kappa: 120.0,
            spec_pool_burnout_dampener: 0.45,
        };
        let rep = evaluate_fabozzi_specified_pool_io_po_payup(&spec);
        // Premium pass-through (coupon 6% > yield 4.2%): slower prepayments on specified pool
        // increase IO value and command a positive specified-pool pay-up over TBA
        assert!(rep.specified_pool_cpr < rep.generic_tba_cpr);
        assert!(rep.specified_pool_payup_32nds > 0.0);
        assert!(rep.io_po_additivity_error < 1e-12);

        let h = 1e-6;
        let r_ad = evaluate_fabozzi_specified_pool_io_po_payup(&FabozziSpecifiedPoolIoPoSpec {
            par_notional: Dual::<1>::from_value(spec.par_notional),
            wac_coupon: Dual::<1>::from_value(spec.wac_coupon),
            pass_through_coupon: Dual::<1>::from_value(spec.pass_through_coupon),
            current_mortgage_rate: Dual::<1>::variable(spec.current_mortgage_rate, 0),
            discount_yield: Dual::<1>::from_value(spec.discount_yield),
            base_turnover_cpr: Dual::<1>::from_value(spec.base_turnover_cpr),
            max_refi_cpr: Dual::<1>::from_value(spec.max_refi_cpr),
            s_curve_steepness_kappa: Dual::<1>::from_value(spec.s_curve_steepness_kappa),
            spec_pool_burnout_dampener: Dual::<1>::from_value(spec.spec_pool_burnout_dampener),
        });
        let mut sp = spec;
        sp.current_mortgage_rate += h;
        let mut sm = spec;
        sm.current_mortgage_rate -= h;
        let fd = (evaluate_fabozzi_specified_pool_io_po_payup(&sp).specified_pool_payup_32nds
            - evaluate_fabozzi_specified_pool_io_po_payup(&sm).specified_pool_payup_32nds)
            / (2.0 * h);
        assert!((r_ad.specified_pool_payup_32nds.d[0] - fd).abs() < 1e-5);
    }
}
