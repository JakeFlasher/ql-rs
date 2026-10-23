use ql_ad::{Dual, Real};
use ql_curves::{DiscountCurve, FlatForwardCurve};
use ql_fixed_income::*;
use ql_inflation::ZeroInflationCurve;
use ql_time::Date;

#[test]
fn nss_short_end_and_ols_recovery() {
    let curve = NelsonSiegelSvensson {
        beta0: 0.03,
        beta1: -0.01,
        beta2: 0.02,
        beta3: -0.005,
        tau1: 1.5,
        tau2: 5.0,
    };
    assert!((curve.zero_rate(0.0) - 0.02).abs() < 1e-14);
    assert!((curve.forward_rate(0.0) - 0.02).abs() < 1e-14);
    let times: Vec<f64> = (1..=12).map(|i| i as f64 / 2.0).collect();
    let rates: Vec<f64> = times.iter().map(|&t| curve.zero_rate(t)).collect();
    let fit = NelsonSiegelSvensson::<f64>::fit_ols_for_taus(&times, &rates, 1.5, 5.0);
    for (a, b) in [
        (fit.beta0, curve.beta0),
        (fit.beta1, curve.beta1),
        (fit.beta2, curve.beta2),
        (fit.beta3, curve.beta3),
    ] {
        assert!((a - b).abs() < 1e-10, "{a} {b}");
    }
}

fn bond() -> CallableFixedRateBond<f64> {
    CallableFixedRateBond {
        face_amount: 100.0,
        coupon_rate: 0.05,
        coupon_times: (1..=5).map(|x| x as f64).collect(),
        maturity: 5.0,
        schedule: vec![BondCallability {
            exercise_time: 3.0,
            clean_strike: 100.0,
            is_call: true,
        }],
    }
}

#[test]
fn callable_zero_vol_and_oas() {
    let b = bond();
    let black = b.black_callable_bond_price(0.03, 0.0);
    let tree = b.tree_callable_bond_price(0.03, 0.0, 0.0, 200);
    assert!((black - tree).abs() < 1e-8);
    let market = b.tree_callable_bond_price(0.03, 0.0, 0.01, 200);
    let (oas, _, _) = b.solve_callable_bond_oas(market, 0.03, 0.0, 200);
    assert!((oas - 0.01).abs() < 1e-7, "{oas}");
}

#[test]
fn convertible_limits_and_dual_greeks() {
    let spec = ConvertibleBondSpec {
        spot: 103.0,
        conversion_ratio: 0.0,
        face_value: 100.0,
        coupon_rate: 0.04,
        coupon_times: (1..=3).map(|x| x as f64).collect(),
        maturity: 3.0,
        rate: 0.03,
        div_yield: 0.0,
        credit_spread: 0.02,
        vol: 0.2,
        call_price: None,
        put_price: None,
    };
    let (value, floor) = price_convertible_bond_tf(&spec, 120);
    assert!((value - floor).abs() < 1e-10);
    let mut conv = spec.clone();
    conv.conversion_ratio = 1.0;
    let (value, floor) = price_convertible_bond_tf(&conv, 120);
    assert!(value + 1e-10 >= floor);
    assert!(value + 1e-10 >= conv.spot * conv.conversion_ratio);
    let (ad_value, greeks) = convertible_bond_greeks_dual4(&conv, 120);
    assert!((ad_value - value).abs() < 1e-10);
    let h = 1e-5;
    let mut up = conv.clone();
    let mut dn = conv.clone();
    up.spot += h;
    dn.spot -= h;
    let fd =
        (price_convertible_bond_tf(&up, 120).0 - price_convertible_bond_tf(&dn, 120).0) / (2.0 * h);
    assert!((greeks[0] - fd).abs() < 1e-5);
    let _ = Dual::<4>::from_f64(greeks[0]);
}

#[test]
fn amortizing_cpi_and_agency_mbs_prepayment_oas() {
    let reference = Date::from_ymd(2024, 1, 2).unwrap();
    let periods = 5;
    let coupon = 0.05;
    let notionals = sinking_notionals(periods, coupon, 1, 100.0);
    let payment = notionals[0] * coupon + (notionals[0] - notionals[1]);
    for i in 0..periods {
        let debt_service = notionals[i] * coupon + (notionals[i] - notionals[i + 1]);
        assert!((debt_service - payment).abs() < 1e-12);
    }
    let times: Vec<f64> = (1..=periods).map(|i| i as f64).collect();
    let effective_rate = (1.0 + coupon).ln();
    let curve = DiscountCurve::new(
        reference,
        times.clone(),
        times.iter().map(|t| (-effective_rate * t).exp()).collect(),
    );
    let bond = AmortizingFixedRateBond {
        initial_notional: 100.0,
        coupon_rate: coupon,
        frequency: 1,
        payment_times: times,
    };
    assert!((price_amortizing_bond(&bond, &curve) - 100.0).abs() < 1e-10);

    let inflation =
        ZeroInflationCurve::new(reference, 100.0, vec![1.0, 5.0], vec![-0.03, -0.02], None);
    let cpi_spec = CpiBondSpec {
        face_amount: 100.0,
        base_cpi: 100.0,
        real_coupon_rate: 0.02,
        payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
        accrual_fractions: vec![1.0; 5],
        deflation_floor: true,
    };
    let nominal = FlatForwardCurve::new(reference, 0.02);
    let floored = price_cpi_bond(&cpi_spec, &nominal, &inflation);
    let unfloored = price_cpi_bond(
        &CpiBondSpec {
            deflation_floor: false,
            ..cpi_spec.clone()
        },
        &nominal,
        &inflation,
    );
    assert!(floored > unfloored);
    assert!(breakeven_inflation_rate(&cpi_spec, floored, &nominal).is_finite());

    let mbs_base = AgencyMbsSpec {
        original_balance: 1_000_000.0,
        wac: 0.05,
        pass_through_rate: 0.045,
        term_months: 120,
        seasoning_months: 0,
        prepayment: PrepaymentModel::Psa {
            psa_multiplier: 0.0,
        },
    };
    let (zero_psa_price, zero_psa_wal) = price_agency_mbs(&mbs_base, &nominal, 0.0);
    let (_, one_psa_wal) = price_agency_mbs(
        &AgencyMbsSpec {
            prepayment: PrepaymentModel::Psa {
                psa_multiplier: 1.0,
            },
            ..mbs_base.clone()
        },
        &nominal,
        0.0,
    );
    let (_, three_psa_wal) = price_agency_mbs(
        &AgencyMbsSpec {
            prepayment: PrepaymentModel::Psa {
                psa_multiplier: 3.0,
            },
            ..mbs_base.clone()
        },
        &nominal,
        0.0,
    );
    assert!(zero_psa_wal > one_psa_wal && one_psa_wal > three_psa_wal);
    let oas = 0.0075;
    let market = price_agency_mbs(&mbs_base, &nominal, oas).0;
    let solved = solve_mbs_oas(&mbs_base, &nominal, market).unwrap();
    assert!((solved - oas).abs() < 1e-11);
    assert!(zero_psa_price.is_finite());
}

#[test]
fn cat_bond_loss_free_and_zero_intensity_limits() {
    let spec = CatBondSpec {
        face_amount: 100.0,
        coupon_rate: 0.05,
        payment_times: vec![1.0, 2.0, 3.0],
        accrual_fractions: vec![1.0; 3],
        discount_rate: 0.03,
        risk: NotionalRisk::Digital { threshold: 20.0 },
        peril: CatPerilSpec::PoissonExponential {
            lambda: 0.0,
            mean_severity: 10.0,
        },
    };
    let valuation = spec.price_analytical();
    assert!((valuation.price - spec.loss_free_price()).abs() < 1e-12);
    assert!(valuation.expected_loss.abs() < 1e-12);
}
