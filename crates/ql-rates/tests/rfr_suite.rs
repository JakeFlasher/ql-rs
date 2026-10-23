use ql_ad::{Dual, Real};
use ql_curves::FlatForwardCurve;
use ql_rates::*;
use ql_time::Date;

#[test]
fn willems_effective_params_and_limits() {
    let p = willems_backward_sabr_params(0.5, 1.0, 0.10, 1.0, -0.50, 0.50, 1.0);
    assert!((p.alpha_hat - 0.08203).abs() < 5e-4);
    assert!((p.rho_hat + 0.5030).abs() < 5e-4);
    assert!((p.nu_hat - 0.4112).abs() < 5e-4);
    let left = willems_backward_sabr_params(-1e-12, 1.0, 0.10, 0.5, -0.3, 0.4, 1.0);
    let right = willems_backward_sabr_params(0.0, 1.0, 0.10, 0.5, -0.3, 0.4, 1.0);
    assert!((left.alpha_hat - right.alpha_hat).abs() < 1e-11);
    assert!((left.rho_hat - right.rho_hat).abs() < 1e-11);
    assert!((left.nu_hat - right.nu_hat).abs() < 1e-11);
}

#[test]
fn willems_forward_looking_uses_original_sabr_and_limits() {
    let base = RfrCapletSpec {
        notional: 1.0,
        forward_rate: 0.05,
        strike: 0.05,
        tau_0: 0.5,
        tau_1: 1.0,
        accrual_fraction: 0.5,
        discount_tau_1: 0.98,
        is_caplet: true,
        rate_type: RfrRateType::ForwardLooking,
    };
    let backward = RfrCapletSpec {
        rate_type: RfrRateType::BackwardLooking,
        ..base
    };
    let forward = willems_backward_sabr_caplet_price(&base, 0.10, 1.0, -0.50, 0.50, 1.0);
    let backward_value = willems_backward_sabr_caplet_price(&backward, 0.10, 1.0, -0.50, 0.50, 1.0);
    assert!(backward_value > forward);

    let backward_large_q =
        willems_backward_sabr_caplet_price(&backward, 0.10, 1.0, -0.50, 0.50, 100_000.0);
    assert!((backward_large_q - forward).abs() < 1e-5);

    let strikes = [0.03, 0.04, 0.05, 0.06, 0.07];
    let strip = willems_backward_sabr_caplet_strip(&base, &strikes, 0.10, 1.0, -0.50, 0.50, 1.0);
    for (&strike, &value) in strikes.iter().zip(&strip) {
        let mut s = base;
        s.strike = strike;
        let direct = willems_backward_sabr_caplet_price(&s, 0.10, 1.0, -0.50, 0.50, 1.0);
        assert!((value - direct).abs() < 1e-12);
    }
}

#[test]
fn rfr_caplet_parity_and_ad() {
    let spec = RfrCapletSpec {
        notional: 1_000_000.0,
        forward_rate: 0.045,
        strike: 0.04,
        tau_0: 0.5,
        tau_1: 1.0,
        accrual_fraction: 0.5,
        discount_tau_1: 0.98,
        is_caplet: true,
        rate_type: RfrRateType::BackwardLooking,
    };
    let cap = price_rfr_caplet_black(&spec, 0.3, 0.0, 1.0);
    let mut floor = spec;
    floor.is_caplet = false;
    let flr = price_rfr_caplet_black(&floor, 0.3, 0.0, 1.0);
    let diff = cap - flr;
    let expected = spec.notional
        * spec.accrual_fraction
        * spec.discount_tau_1
        * (spec.forward_rate - spec.strike);
    assert!((diff - expected).abs() < 1e-8);
    let dspec = RfrCapletSpec {
        notional: Dual::<4>::variable(1_000_000.0, 0),
        forward_rate: Dual::variable(0.045, 1),
        strike: Dual::variable(0.04, 2),
        tau_0: Dual::from_f64(0.5),
        tau_1: Dual::from_f64(1.0),
        accrual_fraction: Dual::from_f64(0.5),
        discount_tau_1: Dual::variable(0.98, 3),
        is_caplet: true,
        rate_type: RfrRateType::BackwardLooking,
    };
    let x = price_rfr_caplet_bachelier(&dspec, Dual::from_f64(0.01), Dual::from_f64(1.0));
    assert!(x.val().is_finite() && x.d.iter().all(|v| v.is_finite()));
}

#[test]
fn hull_white_rfr_and_futures_are_finite() {
    let d = Date::from_ymd(2024, 1, 1).unwrap();
    let c = FlatForwardCurve::new(d, 0.03);
    let x = hw1f_rfr_caplet_price(
        &c,
        0.1,
        0.01,
        0.5,
        1.0,
        0.04,
        1_000_000.0,
        true,
        RfrRateType::BackwardLooking,
    );
    assert!(x.is_finite() && x > 0.0);
    let f = hw1f_sofr_futures_rate_3m(&c, 0.1, 0.01, 0.5, 0.75);
    assert!(f.is_finite());
}
