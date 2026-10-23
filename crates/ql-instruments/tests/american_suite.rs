use ql_instruments::{
    AmericanOption, BaroneAdesiWhaleyEngine, JuQuadraticEngine, OptionType, american_implied_vol,
};
use ql_vol::black_scholes_price;

#[test]
fn american_bounds_and_zero_dividend_call_reduction() {
    let am = BaroneAdesiWhaleyEngine::price(100.0, 100.0, 0.03, 0.01, 0.2, 1.0, true);
    let eu = black_scholes_price(100.0, 100.0, 0.2, 1.0, 0.03, 0.01, true);
    assert!(am >= eu && am >= 0.0);
    let am_q0 = BaroneAdesiWhaleyEngine::price(100.0, 100.0, 0.03, 0.0, 0.2, 1.0, true);
    let eu_q0 = black_scholes_price(100.0, 100.0, 0.2, 1.0, 0.03, 0.0, true);
    assert!((am_q0 - eu_q0).abs() < 1e-12);
}

#[test]
fn american_implied_vol_round_trip() {
    let target = JuQuadraticEngine::price(100.0, 95.0, 0.03, 0.01, 0.24, 1.5, false);
    let implied = american_implied_vol(target, 100.0, 95.0, 0.03, 0.01, 1.5, false).unwrap();
    assert!((implied - 0.24).abs() < 1e-8);
    let option = AmericanOption::new(100.0, 95.0, 1.5, 0.24, 0.03, OptionType::Put);
    assert!(option.price::<f64>() >= 0.0);
}
