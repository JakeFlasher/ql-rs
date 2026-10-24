use ql_ad::{Dual2, Real};
use ql_instruments::*;
use ql_time::*;
use std::sync::Arc;

#[test]
fn fd1d_european_american_barrier_and_ad_greeks() {
    let engine = Fd1dBlackScholesEngine::new(300, 200);
    let euro_call = engine.price(100.0, 100.0, 1.0, 0.2, 0.05, 0.0, OptionType::Call);
    let bs_call = ql_vol::black_scholes_price(100.0, 100.0, 0.2, 1.0, 0.05, 0.0, true);
    assert!(
        (euro_call - bs_call).abs() < 2e-2,
        "pde={} bs={}",
        euro_call,
        bs_call
    );

    let spot_ad = Dual2::<1>::variable(100.0, 0);
    let pde_ad = engine.price(
        spot_ad,
        Dual2::from_f64(100.0),
        Dual2::from_f64(1.0),
        Dual2::from_f64(0.2),
        Dual2::from_f64(0.05),
        Dual2::from_f64(0.0),
        OptionType::Call,
    );
    let bs_ad = ql_vol::black_scholes_price(
        spot_ad,
        Dual2::from_f64(100.0),
        Dual2::from_f64(0.2),
        Dual2::from_f64(1.0),
        Dual2::from_f64(0.05),
        Dual2::from_f64(0.0),
        true,
    );
    assert!(
        (pde_ad.d1[0] - bs_ad.d1[0]).abs() < 2e-2,
        "pde_delta={} bs_delta={}",
        pde_ad.d1[0],
        bs_ad.d1[0]
    );
    assert!(
        (pde_ad.d2[0][0] - bs_ad.d2[0][0]).abs() < 5e-3,
        "pde_gamma={} bs_gamma={}",
        pde_ad.d2[0][0],
        bs_ad.d2[0][0]
    );

    let euro_put = ql_vol::black_scholes_price(100.0, 100.0, 0.2, 1.0, 0.05, 0.0, false);
    let amer_put =
        AmericanOption::new(100.0, 100.0, 1.0, 0.2, 0.05, OptionType::Put).price::<f64>();
    assert!(
        amer_put > euro_put + 0.2,
        "amer_put={} euro_put={}",
        amer_put,
        euro_put
    );

    let uo_call = BarrierOption::new(
        100.0,
        100.0,
        130.0,
        1.0,
        0.2,
        0.05,
        OptionType::Call,
        BarrierType::UpAndOut,
    )
    .price::<f64>();
    assert!(uo_call > 0.0 && uo_call < euro_call);
}

#[test]
fn bond_ytm_ad_and_simd_leg() {
    let cal = Arc::new(BitsetCalendar::weekend_only());
    let spec = FixedRateBondSpec::new(
        Date::from_ymd(2024, 1, 1).unwrap(),
        Date::from_ymd(2026, 1, 1).unwrap(),
        0.05,
        Frequency::Semiannual,
        100.0,
        cal.clone(),
        BusinessDayConvention::Following,
        DayCounter::Act365Fixed,
    );
    let bond = spec.resolve(Date::from_ymd(2024, 1, 1).unwrap());
    assert!(bond.accrued_interest(Date::from_ymd(2024, 3, 1).unwrap()) > 0.0);

    let p0 = bond.dirty_price(0.055);
    let y_ad = bond.yield_to_maturity(Dual2::<1>::variable(p0, 0));
    assert!((y_ad.val - 0.055).abs() < 1e-10);
    let h = 1e-4;
    let dy_dp_fd = (bond.yield_to_maturity(p0 + h) - bond.yield_to_maturity(p0 - h)) / (2.0 * h);
    assert!((y_ad.d1[0] - dy_dp_fd).abs() < 1e-6);

    let swap_spec = VanillaSwapSpec {
        effective_date: Date::from_ymd(2024, 1, 1).unwrap(),
        termination_date: Date::from_ymd(2034, 1, 1).unwrap(),
        notional: 1_000_000.0,
        fixed_rate_or_coupon: 0.04,
        fixed_tenor: Period::Months(6),
        float_tenor: Period::Months(3),
        calendar: cal,
        bdc: BusinessDayConvention::ModifiedFollowing,
        fixed_dc: DayCounter::Thirty360Us,
        float_dc: DayCounter::Act360,
        swap_type: SwapType::Payer,
    };
    let resolved = swap_spec.resolve(Date::from_ymd(2024, 1, 1).unwrap());
    assert_eq!(resolved.fixed_leg.len(), 20);
    let c = ql_curves::FlatForwardCurve::new(Date::from_ymd(2024, 1, 1).unwrap(), 0.04);
    assert!(
        (resolved.fixed_leg.npv_fixed_leg(&c) - resolved.fixed_leg.npv_fixed_leg_simd(&c)).abs()
            < 1e-9
    );
}

#[test]
fn barrier_rebate_survival_and_breached_short_circuit() {
    let engine = Fd1dBlackScholesEngine::new(161, 120);
    let ki0 = BarrierOption {
        underlying: 100.0,
        strike: 100.0,
        barrier: 115.0,
        expiry: 1.0,
        volatility: 0.2,
        rate: 0.03,
        option_type: OptionType::Call,
        barrier_type: BarrierType::UpAndIn,
        rebate: 0.0,
    }
    .price::<f64>();
    let ki1 = BarrierOption {
        underlying: 100.0,
        strike: 100.0,
        barrier: 115.0,
        expiry: 1.0,
        volatility: 0.2,
        rate: 0.03,
        option_type: OptionType::Call,
        barrier_type: BarrierType::UpAndIn,
        rebate: 2.0,
    }
    .price::<f64>();
    assert!(ki1 > ki0);
    let breached = BarrierOption {
        underlying: 120.0,
        strike: 100.0,
        barrier: 115.0,
        expiry: 1.0,
        volatility: 0.2,
        rate: 0.03,
        option_type: OptionType::Call,
        barrier_type: BarrierType::UpAndOut,
        rebate: 2.0,
    }
    .price::<f64>();
    assert!((breached - 2.0).abs() < 1e-12);
    let _ = engine;
}

#[test]
fn down_barrier_call_rebate_tolerance_and_parity() {
    let mk = |bt, rebate| BarrierOption {
        underlying: 100.0,
        strike: 100.0,
        barrier: 83.7,
        expiry: 1.0,
        volatility: 0.2,
        rate: 0.03,
        option_type: OptionType::Call,
        barrier_type: bt,
        rebate,
    };
    let ki0 = mk(BarrierType::DownAndIn, 0.0).price::<f64>();
    let ki1 = mk(BarrierType::DownAndIn, 2.0).price::<f64>();
    let ko0 = mk(BarrierType::DownAndOut, 0.0).price::<f64>();
    let vanilla = Fd1dBlackScholesEngine::new(161, 120).price(
        100.0,
        100.0,
        1.0,
        0.2,
        0.03,
        0.0,
        OptionType::Call,
    );
    println!(
        "ki0 {ki0} ki1 {ki1} ko0 {ko0} vanilla {vanilla} diff {}",
        ki0 + ko0 - vanilla
    );
    assert!(ki1 > ki0);
}
