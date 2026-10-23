use ql_ad::Real;
use ql_credit::*;
use ql_curves::DiscountCurve;
use ql_time::Date;

fn strip() -> (Date, DiscountCurve<f64>, Vec<CdsQuote>) {
    let d = Date::from_ymd(2024, 1, 1).unwrap();
    let maturities = vec![1.0, 3.0, 5.0, 7.0, 10.0];
    let discount = DiscountCurve::new(
        d,
        maturities.clone(),
        maturities.iter().map(|t| (-0.02 * t).exp()).collect(),
    );
    let quotes = maturities
        .iter()
        .enumerate()
        .map(|(i, &maturity_years)| CdsQuote {
            maturity_years,
            par_spread: 0.01 + 0.001 * i as f64,
        })
        .collect();
    (d, discount, quotes)
}

fn benchmark(q: CdsQuote) -> CreditDefaultSwap {
    let n = (q.maturity_years * 4.0) as usize;
    CreditDefaultSwap {
        side: ProtectionSide::Buyer,
        notional: 1.0,
        running_spread: q.par_spread,
        recovery_rate: 0.4,
        payment_times: (1..=n).map(|i| i as f64 / 4.0).collect(),
        accrual_fractions: vec![0.25; n],
    }
}

#[test]
fn isda_cds_bootstrap_reprices_benchmarks_and_hazard_intervals() {
    let (d, discount, quotes) = strip();
    let calibrated = CalibratedCreditCurve::bootstrap(d, discount.clone(), 0.4, &quotes).unwrap();
    for q in &quotes {
        let v = IsdaCdsEngine::price(&benchmark(*q), &discount, &calibrated.curve);
        assert!(v.npv.abs() < 1e-12, "{}", v.npv);
        assert!((v.fair_spread - q.par_spread).abs() < 1e-12);
    }
    let h = &calibrated.curve.hazard_rates;
    assert_eq!(calibrated.curve.hazard_rate(0.5), h[0]);
    assert_eq!(calibrated.curve.hazard_rate(1.0), h[0]);
    assert_eq!(calibrated.curve.hazard_rate(1.000001), h[1]);
    assert_eq!(calibrated.curve.hazard_rate(2.9), h[1]);
    assert_eq!(calibrated.curve.hazard_rate(3.000001), h[2]);
}

#[test]
fn ift_par_cs01_matches_rebootstrap_central_fd() {
    let (d, discount, quotes) = strip();
    let calibrated = CalibratedCreditCurve::bootstrap(d, discount.clone(), 0.4, &quotes).unwrap();
    let portfolio: Vec<CreditDefaultSwap> = (1..=8)
        .map(|i| CreditDefaultSwap {
            side: if i % 2 == 0 {
                ProtectionSide::Buyer
            } else {
                ProtectionSide::Seller
            },
            notional: 1_000_000.0,
            running_spread: 0.009 + 0.0007 * i as f64,
            recovery_rate: 0.4,
            payment_times: (1..=(i * 4)).map(|j| j as f64 / 4.0).collect(),
            accrual_fractions: vec![0.25; i * 4],
        })
        .collect();
    let (value, cs01) = calibrated.value_and_par_cs01(|curve| {
        let discount_ad = discount.to_ad_curve(
            &discount
                .discounts
                .iter()
                .map(|&x| ql_ad::TapeVar::from_f64(x))
                .collect::<Vec<_>>(),
        );
        portfolio
            .iter()
            .map(|cds| IsdaCdsEngine::price(cds, &discount_ad, curve).npv)
            .fold(ql_ad::TapeVar::from_f64(0.0), |a, x| a + x)
    });
    let direct: f64 = portfolio
        .iter()
        .map(|cds| IsdaCdsEngine::price(cds, &discount, &calibrated.curve).npv)
        .sum();
    assert!((value - direct).abs() < 1e-8);
    let h = 1e-6;
    for j in 0..quotes.len() {
        let mut up = quotes.clone();
        let mut dn = quotes.clone();
        up[j].par_spread += h;
        dn[j].par_spread -= h;
        let cp = SurvivalProbabilityCurve::bootstrap_quotes(d, discount.clone(), 0.4, &up).unwrap();
        let cm = SurvivalProbabilityCurve::bootstrap_quotes(d, discount.clone(), 0.4, &dn).unwrap();
        let vp: f64 = portfolio
            .iter()
            .map(|x| IsdaCdsEngine::price(x, &discount, &cp).npv)
            .sum();
        let vm: f64 = portfolio
            .iter()
            .map(|x| IsdaCdsEngine::price(x, &discount, &cm).npv)
            .sum();
        let fd = (vp - vm) / (2.0 * h);
        assert!(
            (cs01[j] - fd).abs() < 1e-5 * fd.abs().max(1.0),
            "j={j} ad={} fd={fd}",
            cs01[j]
        );
    }
}
