use ql_exotics::{DengLiZhouSpreadBasketEngine, RainbowKind, StulzRainbowEngine};
use ql_vol::black_scholes_price;

#[test]
fn rainbow_min_max_and_spread_limits_are_finite() {
    let min = StulzRainbowEngine::price(
        100.0,
        95.0,
        100.0,
        0.03,
        0.01,
        0.02,
        0.2,
        0.25,
        0.4,
        1.0,
        RainbowKind::CallOnMin,
    );
    let max = StulzRainbowEngine::price(
        100.0,
        95.0,
        100.0,
        0.03,
        0.01,
        0.02,
        0.2,
        0.25,
        0.4,
        1.0,
        RainbowKind::CallOnMax,
    );
    let c1 = black_scholes_price(100.0, 100.0, 0.2, 1.0, 0.03, 0.01, true);
    let c2 = black_scholes_price(95.0, 100.0, 0.25, 1.0, 0.03, 0.02, true);
    assert!((min + max - c1 - c2).abs() < 1e-10);
    assert!(
        DengLiZhouSpreadBasketEngine::price_spread(
            100.0, 95.0, 0.0, 0.03, 0.01, 0.02, 0.2, 0.25, 0.4, 1.0, true
        )
        .is_finite()
    );
}
