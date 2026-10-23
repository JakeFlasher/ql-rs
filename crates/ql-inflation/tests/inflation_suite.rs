use ql_ad::{Dual, Real};
use ql_curves::{DiscountCurve, FlatForwardCurve};
use ql_inflation::*;
use ql_time::Date;

fn setup() -> (
    Date,
    DiscountCurve<f64>,
    MultiplicativeSeasonality,
    Vec<(f64, f64)>,
) {
    let base = Date::from_ymd(2024, 1, 1).unwrap();
    let times = vec![1.0, 2.0, 5.0, 10.0];
    let discount = DiscountCurve::new(
        base,
        times.clone(),
        times.iter().map(|t| (-0.02 * t).exp()).collect(),
    );
    let mut factors = [1.0; 12];
    for (i, f) in factors.iter_mut().enumerate().take(11) {
        *f = 1.0 + 0.005 * ((i as f64) - 5.0) / 5.0;
    }
    factors[11] = 1.0 / factors[..11].iter().product::<f64>();
    let seasonality = MultiplicativeSeasonality::new(factors).unwrap();
    let quotes = vec![(1.0, 0.02), (2.0, 0.022), (5.0, 0.025), (10.0, 0.027)];
    (base, discount, seasonality, quotes)
}

#[test]
fn zciis_bootstrap_seasonality_and_daily_cpi_interpolation() {
    let (base, discount, seasonality, quotes) = setup();
    let calibrated = CalibratedInflationCurve::bootstrap_with_seasonality(
        base,
        100.0,
        discount.clone(),
        &quotes,
        Some(seasonality.clone()),
    )
    .unwrap();
    for &(t, k) in &quotes {
        let swap = ZeroCouponInflationSwap {
            notional: 1.0,
            fixed_rate: k,
            maturity_years: t,
            target_month: None,
            pay_fixed: true,
        };
        assert!(swap.npv(&discount, &calibrated.curve).abs() < 1e-12);
        assert!((swap.fair_rate(&calibrated.curve) - k).abs() < 1e-12);
    }
    let jan = calibrated.curve.forward_cpi(5.0, Some(1));
    let jul = calibrated.curve.forward_cpi(5.0, Some(7));
    let ratio = jul / jan;
    assert!(
        (ratio.val() - seasonality.factor_for_month(7) / seasonality.factor_for_month(1)).abs()
            < 1e-12
    );
    let x = ZeroInflationCurve::<f64>::interpolated_cpi(100.0, 102.0, 16, 30);
    assert!((x - 101.0).abs() < 1e-12);
}

#[test]
fn ift_inflation_delta_matches_rebootstrap_fd() {
    let (base, discount, seasonality, quotes) = setup();
    let calibrated = CalibratedInflationCurve::bootstrap_with_seasonality(
        base,
        100.0,
        discount.clone(),
        &quotes,
        Some(seasonality.clone()),
    )
    .unwrap();
    let (value, delta) = calibrated.value_and_par_inflation_delta(|curve| {
        let d = discount.to_ad_curve(
            &discount
                .discounts
                .iter()
                .map(|&x| ql_ad::TapeVar::from_f64(x))
                .collect::<Vec<_>>(),
        );
        let z = ZeroCouponInflationSwap {
            notional: 1_000_000.0,
            fixed_rate: 0.021,
            maturity_years: 4.0,
            target_month: Some(9),
            pay_fixed: true,
        };
        let y = YearOnYearInflationSwap {
            notional: 1_000_000.0,
            fixed_rate: 0.024,
            maturity_years: 5,
            pay_fixed: true,
        };
        z.npv(&d, curve) + y.npv(&d, curve)
    });
    assert_eq!(delta.len(), quotes.len());
    assert!(value.abs() < 1e9);
    let h = 1e-6;
    for j in 0..quotes.len() {
        let mut up = quotes.clone();
        let mut dn = quotes.clone();
        up[j].1 += h;
        dn[j].1 -= h;
        let cp = CalibratedInflationCurve::bootstrap_with_seasonality(
            base,
            100.0,
            discount.clone(),
            &up,
            Some(seasonality.clone()),
        )
        .unwrap();
        let cm = CalibratedInflationCurve::bootstrap_with_seasonality(
            base,
            100.0,
            discount.clone(),
            &dn,
            Some(seasonality.clone()),
        )
        .unwrap();
        let z = ZeroCouponInflationSwap {
            notional: 1_000_000.0,
            fixed_rate: 0.021,
            maturity_years: 4.0,
            target_month: Some(9),
            pay_fixed: true,
        };
        let y = YearOnYearInflationSwap {
            notional: 1_000_000.0,
            fixed_rate: 0.024,
            maturity_years: 5,
            pay_fixed: true,
        };
        let vp = z.npv(&discount, &cp.curve) + y.npv(&discount, &cp.curve);
        let vm = z.npv(&discount, &cm.curve) + y.npv(&discount, &cm.curve);
        let fd = (vp - vm) / (2.0 * h);
        assert!((delta[j] - fd).abs() < 1e-6 * fd.abs().max(1.0), "j={j}");
    }
}

#[test]
fn legacy_flat_curve_constructor_remains_valid() {
    let d = Date::from_ymd(2024, 1, 1).unwrap();
    let c = ZeroInflationCurve::new(d, 100.0, vec![1.0, 5.0], vec![0.02, 0.025], None);
    let dc = FlatForwardCurve::new(d, 0.02);
    let z = ZeroCouponInflationSwap {
        notional: 1e6,
        fixed_rate: 0.02,
        maturity_years: 1.,
        target_month: None,
        pay_fixed: true,
    };
    let _ = z.npv(&dc, &c);
}

#[test]
fn yoy_inflation_cap_floor_collar_parity_and_dual_greeks() {
    let base = Date::from_ymd(2024, 1, 1).unwrap();
    let discount = FlatForwardCurve::new(base, 0.02);
    let inflation = ZeroInflationCurve::new(
        base,
        100.0,
        vec![1.0, 2.0, 3.0],
        vec![0.025, 0.027, 0.029],
        None,
    );
    let maturity = 3;
    let strike = 0.02;
    let common = |kind| YoYInflationCapFloor {
        notional: 1_000_000.0,
        maturity_years: maturity,
        cap_rate: Some(strike),
        floor_rate: Some(strike),
        gearing: 1.0,
        spread: 0.0,
        cap_floor_type: kind,
    };
    let yoy = YearOnYearInflationSwap {
        notional: 1_000_000.0,
        fixed_rate: strike,
        maturity_years: maturity,
        pay_fixed: true,
    };
    for model in [
        YoYOptionletModel::Black,
        YoYOptionletModel::UnitDisplacedBlack,
        YoYOptionletModel::Bachelier,
    ] {
        let cap = price_yoy_inflation_capfloor(
            &common(YoYInflationCapFloorType::Cap),
            &discount,
            &inflation,
            0.10,
            model,
        );
        let floor = price_yoy_inflation_capfloor(
            &common(YoYInflationCapFloorType::Floor),
            &discount,
            &inflation,
            0.10,
            model,
        );
        let swap = yoy.npv(&discount, &inflation);
        assert!((cap - floor - swap).abs() < 1e-10, "model {model:?}");

        let collar = YoYInflationCapFloor {
            cap_rate: Some(0.035),
            floor_rate: Some(0.015),
            cap_floor_type: YoYInflationCapFloorType::Collar,
            ..common(YoYInflationCapFloorType::Collar)
        };
        let cap_hi = YoYInflationCapFloor {
            cap_rate: Some(0.035),
            cap_floor_type: YoYInflationCapFloorType::Cap,
            ..common(YoYInflationCapFloorType::Cap)
        };
        let floor_lo = YoYInflationCapFloor {
            floor_rate: Some(0.015),
            cap_floor_type: YoYInflationCapFloorType::Floor,
            ..common(YoYInflationCapFloorType::Floor)
        };
        let collar_value =
            price_yoy_inflation_capfloor(&collar, &discount, &inflation, 0.10, model);
        let decomposition =
            price_yoy_inflation_capfloor(&cap_hi, &discount, &inflation, 0.10, model)
                - price_yoy_inflation_capfloor(&floor_lo, &discount, &inflation, 0.10, model);
        assert!((collar_value - decomposition).abs() < 1e-10);
    }

    let ad_inflation = ZeroInflationCurve::new(
        base,
        100.0,
        vec![1.0, 2.0, 3.0],
        vec![
            Dual::<3>::variable(0.025, 1),
            Dual::<3>::variable(0.027, 2),
            Dual::from_value(0.029),
        ],
        None,
    );
    let ad_discount = FlatForwardCurve::new(base, Dual::<3>::from_value(0.02));
    let ad_instrument = YoYInflationCapFloor {
        notional: Dual::from_value(1.0),
        maturity_years: maturity,
        cap_rate: Some(Dual::from_value(strike)),
        floor_rate: Some(Dual::from_value(strike)),
        gearing: Dual::from_value(1.0),
        spread: Dual::from_value(0.0),
        cap_floor_type: YoYInflationCapFloorType::Cap,
    };
    let ad_value = price_yoy_inflation_capfloor(
        &ad_instrument,
        &ad_discount,
        &ad_inflation,
        Dual::<3>::variable(0.10, 0),
        YoYOptionletModel::Bachelier,
    );
    let h = 1e-6;
    let unit = |kind| YoYInflationCapFloor {
        notional: 1.0,
        maturity_years: maturity,
        cap_rate: Some(strike),
        floor_rate: Some(strike),
        gearing: 1.0,
        spread: 0.0,
        cap_floor_type: kind,
    };
    for lane in 0..3 {
        let mut up = inflation.clone();
        let mut down = inflation.clone();
        let mut vol_up = 0.10;
        let mut vol_down = 0.10;
        if lane == 0 {
            vol_up += h;
            vol_down -= h;
        } else {
            up.zero_rates[lane - 1] += h;
            down.zero_rates[lane - 1] -= h;
        }
        let up_value = price_yoy_inflation_capfloor(
            &unit(YoYInflationCapFloorType::Cap),
            &discount,
            &up,
            vol_up,
            YoYOptionletModel::Bachelier,
        );
        let down_value = price_yoy_inflation_capfloor(
            &unit(YoYInflationCapFloorType::Cap),
            &discount,
            &down,
            vol_down,
            YoYOptionletModel::Bachelier,
        );
        let fd = (up_value - down_value) / (2.0 * h);
        assert!(
            (ad_value.d[lane] - fd).abs() < 1e-8,
            "lane {lane}: {} vs {fd}",
            ad_value.d[lane]
        );
    }
}
