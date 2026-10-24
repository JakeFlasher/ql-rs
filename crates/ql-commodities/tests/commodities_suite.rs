use ql_ad::Dual;
use ql_commodities::*;
use ql_instruments::OptionType;
#[test]
fn seasonal_curve_and_ift_delta() {
    let s = SeasonalParams {
        a1: 0.01,
        b1: 0.02,
        a2: 0.003,
        b2: 0.004,
    };
    let c = CommodityFuturesCurve::from_futures_pillars(vec![1., 2., 3.], &[80., 82., 85.], s);
    for (i, &t) in [1., 2., 3.].iter().enumerate() {
        assert!((c.futures_price(t) - [80., 82., 85.][i]).abs() < 1e-12);
    }
    let cal = CalibratedCommodityCurve::new(vec![1., 2., 3.], vec![80., 82., 85.], s);
    let (base, d) = cal.dp_dfutures(|curve| curve.futures_price(Dual::from_value(2.5)));
    assert!((base - cal.curve().futures_price(2.5)).abs() < 1e-12);
    for j in 0..3 {
        let h = 1e-5;
        let mut up = [80., 82., 85.];
        let mut dn = up;
        up[j] += h;
        dn[j] -= h;
        let vp = CalibratedCommodityCurve::new(vec![1., 2., 3.], up.to_vec(), s)
            .curve()
            .futures_price(2.5);
        let vm = CalibratedCommodityCurve::new(vec![1., 2., 3.], dn.to_vec(), s)
            .curve()
            .futures_price(2.5);
        assert!((d[j] - (vp - vm) / (2. * h)).abs() < 1e-6);
    }
}
#[test]
fn schwartz_and_spread_parity() {
    let p = SchwartzSmith2F {
        chi0: 0.,
        xi0: 4.,
        kappa: 1.,
        sigma_chi: 0.3,
        mu_xi: 0.,
        sigma_xi: 0.2,
        rho: -0.3,
    };
    let f = p.futures_price(1.);
    assert!(f > 0.);
    assert!((p.a_term(1.) - p.log_futures_variance(1., 1.) * 0.5).abs() < 1e-14);
    let v = p.price_futures_option(1., 1.5, f, 0.03, OptionType::Call);
    assert!(v >= 0.);
    let m = margrabe_exchange_option(100., 90., 0.2, 0.25, 0.2, 0.03, 1., 1., 1.);
    let k = kirk_spread_option(100., 90., 0., 0.2, 0.25, 0.2, 0.03, 1., OptionType::Call);
    assert!((m - k).abs() < 1e-12);
}

#[test]
fn swing_forward_strip_and_storage_dp_are_tree_based() {
    let curve = CommodityFuturesCurve::from_futures_pillars(
        vec![1.0, 2.0, 3.0],
        &[80.0, 81.0, 82.0],
        SeasonalParams {
            a1: 0.0,
            b1: 0.0,
            a2: 0.0,
            b2: 0.0,
        },
    );
    let swing = SwingOptionSpec {
        exercise_times: vec![1.0, 2.0, 3.0],
        strike: 79.0,
        unit_volume: 1.0,
        min_exercises: 3,
        max_exercises: 3,
        is_call: true,
    };
    let value = price_swing_option(&swing, &curve, 1.0, 0.3, 0.02, 3);
    let strip = (1..=3)
        .map(|i| (-0.02 * i as f64).exp() * (curve.futures_price(i as f64) - 79.0))
        .sum::<f64>();
    assert!((value - strip).abs() < 1e-10);
    let storage = GasStorageSpec {
        decision_times: vec![1.0, 2.0, 3.0],
        max_capacity: 10.0,
        initial_inventory: 5.0,
        max_injection_rate: 5.0,
        max_withdrawal_rate: 5.0,
        injection_loss: 0.0,
        withdrawal_loss: 0.0,
        injection_cost: 0.0,
        withdrawal_cost: 0.0,
        inventory_grid_points: 10,
    };
    let low = price_gas_storage(&storage, &curve, 1.0, 0.1, 0.02);
    let high = price_gas_storage(&storage, &curve, 1.0, 0.5, 0.02);
    assert!(low.is_finite() && high.is_finite());
}

#[test]
fn gas_storage_terminal_decisions_have_value() {
    let curve = CommodityFuturesCurve::from_futures_pillars(
        vec![1.0, 2.0],
        &[50.0, 100.0],
        SeasonalParams {
            a1: 0.0,
            b1: 0.0,
            a2: 0.0,
            b2: 0.0,
        },
    );
    let base = GasStorageSpec {
        decision_times: vec![1.0, 2.0],
        max_capacity: 10.0,
        initial_inventory: 0.0,
        max_injection_rate: 10.0,
        max_withdrawal_rate: 10.0,
        injection_loss: 0.0,
        withdrawal_loss: 0.0,
        injection_cost: 0.0,
        withdrawal_cost: 0.0,
        inventory_grid_points: 10,
    };
    assert!(price_gas_storage(&base, &curve, 1.0, 0.1, 0.02) > 0.0);
    let one = GasStorageSpec {
        decision_times: vec![1.0],
        initial_inventory: 5.0,
        ..base
    };
    assert!(price_gas_storage(&one, &curve, 1.0, 0.1, 0.02) > 0.0);
}
