use ql_ad::{Real, TapeVar};
use ql_curves::{DiscountCurve, YieldTermStructure};
use ql_fx::*;
use ql_time::Date;

fn curves() -> (Date, DiscountCurve<f64>, DiscountCurve<f64>) {
    let reference = Date::from_ymd(2024, 1, 2).unwrap();
    let times: Vec<f64> = (1..=5).map(|i| i as f64).collect();
    let domestic = DiscountCurve::new(
        reference,
        times.clone(),
        times.iter().map(|t| (-0.025 * t).exp()).collect(),
    );
    let foreign = DiscountCurve::new(
        reference,
        times.clone(),
        times.iter().map(|t| (-0.015 * t).exp()).collect(),
    );
    (reference, domestic, foreign)
}

fn quotes() -> Vec<XccsQuote> {
    (1..=5)
        .map(|year| XccsQuote {
            maturity_years: year as f64,
            basis_spread: 0.0005 + year as f64 * 0.0001,
        })
        .collect()
}

fn calibration_instrument(quote: &XccsQuote, spot: f64) -> ConstNotionalXccs {
    let n = quote.maturity_years.round() as usize;
    ConstNotionalXccs {
        dom_notional: spot,
        for_notional: 1.0,
        spot_fx: spot,
        payment_times: (1..=n).map(|i| i as f64).collect(),
        accrual_fractions: vec![1.0; n],
        for_basis_spread: quote.basis_spread,
        dom_spread: 0.0,
        receive_foreign: true,
    }
}

#[test]
fn fx_forward_points_and_csa_df_round_trip() {
    let spot = 1.10;
    let dom = 0.96;
    let foreign = 0.98;
    let forward = fx_forward_rate(spot, dom, foreign);
    assert!((fx_swap_points(spot, dom, foreign) - (forward - spot)).abs() < 1e-14);
    assert!((implied_foreign_csa_df_from_fx_forward(spot, forward, dom) - foreign).abs() < 1e-14);
}

#[test]
fn xccs_bootstrap_reprices_and_zero_basis_degenerates_to_ois() {
    let (reference, domestic, foreign) = curves();
    let spot = 1.10;
    let quotes = quotes();
    let calibrated =
        CalibratedXccsCurve::bootstrap(reference, spot, domestic.clone(), foreign.clone(), &quotes)
            .unwrap();
    for quote in &quotes {
        let instrument = calibration_instrument(quote, spot);
        let npv = instrument.npv(&domestic, &foreign, &calibrated.for_csa_curve);
        assert!(npv.abs() < 1e-12, "{}Y NPV={npv}", quote.maturity_years);
        let fair =
            instrument.fair_foreign_basis_spread(&domestic, &foreign, &calibrated.for_csa_curve);
        assert!((fair - quote.basis_spread).abs() < 1e-12);
    }

    let zero_quotes: Vec<XccsQuote> = quotes
        .iter()
        .map(|q| XccsQuote {
            maturity_years: q.maturity_years,
            basis_spread: 0.0,
        })
        .collect();
    let zero = CalibratedXccsCurve::bootstrap(
        reference,
        spot,
        domestic.clone(),
        foreign.clone(),
        &zero_quotes,
    )
    .unwrap();
    for &t in &zero.for_csa_curve.times {
        assert!((zero.for_csa_curve.discount_factor(t) - foreign.discount_factor(t)).abs() < 1e-12);
    }
}

#[test]
fn xccs_ift_basis_dv01_matches_bump_and_rebootstrap() {
    let (reference, domestic, foreign) = curves();
    let spot = 1.10;
    let quotes = quotes();
    let calibrated =
        CalibratedXccsCurve::bootstrap(reference, spot, domestic.clone(), foreign.clone(), &quotes)
            .unwrap();
    let portfolio = ConstNotionalXccs {
        dom_notional: 1_000.0,
        for_notional: 900.0,
        spot_fx: spot,
        payment_times: (1..=5).map(|i| i as f64).collect(),
        accrual_fractions: vec![1.0; 5],
        for_basis_spread: 0.0021,
        dom_spread: 0.0003,
        receive_foreign: true,
    };
    let (value, sensitivities) = calibrated.value_and_par_basis_dv01(|foreign_csa| {
        let domestic_ad = domestic.to_ad_curve(
            &domestic
                .discounts
                .iter()
                .map(|&x| TapeVar::from_f64(x))
                .collect::<Vec<_>>(),
        );
        let foreign_ad = foreign.to_ad_curve(
            &foreign
                .discounts
                .iter()
                .map(|&x| TapeVar::from_f64(x))
                .collect::<Vec<_>>(),
        );
        portfolio.npv(&domestic_ad, &foreign_ad, foreign_csa)
    });
    assert!(value.is_finite());
    assert_eq!(sensitivities.len(), quotes.len());

    let h = 1e-6;
    for (j, sensitivity) in sensitivities.iter().enumerate() {
        let mut up = quotes.clone();
        let mut down = quotes.clone();
        up[j].basis_spread += h;
        down[j].basis_spread -= h;
        let cup =
            CalibratedXccsCurve::bootstrap(reference, spot, domestic.clone(), foreign.clone(), &up)
                .unwrap();
        let cdown = CalibratedXccsCurve::bootstrap(
            reference,
            spot,
            domestic.clone(),
            foreign.clone(),
            &down,
        )
        .unwrap();
        let vp = portfolio.npv(&domestic, &foreign, &cup.for_csa_curve);
        let vm = portfolio.npv(&domestic, &foreign, &cdown.for_csa_curve);
        let fd = (vp - vm) / (2.0 * h);
        assert!(
            (sensitivity - fd).abs() < 1e-6,
            "pillar {j}: {sensitivity} vs {fd}"
        );
    }
}
