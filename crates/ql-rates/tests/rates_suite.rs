use ql_ad::{Dual, Real};
use ql_curves::{DiscountCurve, YieldTermStructure};
use ql_rates::*;

fn curve() -> (DiscountCurve<f64>, DiscountCurve<f64>) {
    let d = ql_time::Date::from_ymd(2024, 1, 1).unwrap();
    let t = vec![0.25, 0.5, 1., 2., 5., 10.];
    let df = t.iter().map(|x| (-0.02 * x).exp()).collect();
    let f = t.iter().map(|x| (-0.025 * x).exp()).collect();
    (
        DiscountCurve::new(d, t.clone(), df),
        DiscountCurve::new(d, t, f),
    )
}
#[test]
fn caplet_price_and_strip_reprice() {
    let (d, f) = curve();
    let cap = IborCapFloor::from_schedule(1e6, 0.02, &[0.0, 0.25, 0.5, 0.75, 1.0], true);
    let p = cap.price_bachelier(&d, &f, &[0.01; 4]);
    assert!(p > 0.);
    let quotes = (1..=4)
        .map(|m| {
            (
                IborCapFloor::from_schedule(
                    1e6,
                    0.02,
                    &(0..=m).map(|i| i as f64 * 0.25).collect::<Vec<_>>(),
                    true,
                ),
                0.01,
            )
        })
        .collect::<Vec<_>>();
    let strip = CalibratedCapletVolCurve::bootstrap_bachelier(&d, &f, &quotes).unwrap();
    for (c, q) in &quotes {
        let lhs = c.price_bachelier(&d, &f, &strip.caplet_vols);
        let rhs = c.price_bachelier(&d, &f, &vec![*q; c.caplets.len()]);
        assert!((lhs - rhs).abs() < 1e-6);
    }

    // The IFT strip must map a perturbation of a flat par-cap quote to a
    // *positive* change in the corresponding caplet volatility.  This also
    // guards the sign in the lower-triangular Jacobian solve.
    let (_, cs01) = strip.value_and_par_cap_vega(|vols| {
        let cap = &strip.quotes[3].0;
        cap.price_bachelier(
            &d.to_ad_curve(
                &d.discounts
                    .iter()
                    .map(|&x| ql_ad::TapeVar::from_f64(x))
                    .collect::<Vec<_>>(),
            ),
            &f.to_ad_curve(
                &f.discounts
                    .iter()
                    .map(|&x| ql_ad::TapeVar::from_f64(x))
                    .collect::<Vec<_>>(),
            ),
            vols,
        )
    });
    assert!(cs01.iter().all(|x| *x >= -1e-8));
    let h = 1e-7;
    for j in 0..strip.quotes.len() {
        let mut up = strip.quotes.clone();
        let mut down = strip.quotes.clone();
        up[j].1 += h;
        down[j].1 -= h;
        let up_curve = CalibratedCapletVolCurve::bootstrap_bachelier(&d, &f, &up).unwrap();
        let down_curve = CalibratedCapletVolCurve::bootstrap_bachelier(&d, &f, &down).unwrap();
        let cap = &strip.quotes[3].0;
        let vp = cap.price_bachelier(&d, &f, &up_curve.caplet_vols);
        let vm = cap.price_bachelier(&d, &f, &down_curve.caplet_vols);
        let fd = (vp - vm) / (2.0 * h);
        assert!((cs01[j] - fd).abs() < 1e-4, "j={j} ift={} fd={fd}", cs01[j]);
    }
}
#[test]
fn swaption_settlement_and_sabr_ad() {
    let (d, f) = curve();
    let s = EuropeanSwaption {
        notional: 1e6,
        strike: 0.025,
        expiry: 1.,
        swap_times: vec![1., 2., 3., 4., 5.],
        is_payer: true,
        settlement: SwaptionSettlement::Physical,
    };
    let p = s.price_bachelier(&d, &f, 0.01);
    let c = EuropeanSwaption {
        settlement: SwaptionSettlement::CashParYield,
        ..s.clone()
    };
    assert!(p > 0. && c.price_bachelier(&d, &f, 0.01) > 0.);
    let da = d.to_ad_curve(
        &d.discounts
            .iter()
            .map(|&x| Dual::<3>::from_value(x))
            .collect::<Vec<_>>(),
    );
    let fa = f.to_ad_curve(
        &f.discounts
            .iter()
            .map(|&x| Dual::<3>::from_value(x))
            .collect::<Vec<_>>(),
    );
    let a = s.price_sabr_bachelier(
        &da,
        &fa,
        Dual::variable(0.01, 0),
        0.5,
        Dual::variable(-0.2, 1),
        Dual::variable(0.5, 2),
    );
    assert!(a.val() > 0. && a.d.iter().any(|x| x.abs() > 0.));
}
#[test]
fn g2_and_cheyette_and_cms() {
    let (d, _) = curve();
    let g = G2PlusPlus {
        a: 0.1,
        sigma: 0.01,
        b: 0.3,
        eta: 0.015,
        rho: -0.7,
    };
    let p = g.discount_bond_option(&d, true, 0.95, 1., 5.);
    assert!(p >= 0.);
    let sw = EuropeanSwaption {
        notional: 1e6,
        strike: 0.02,
        expiry: 1.,
        swap_times: (1..=6).map(|x| x as f64).collect(),
        is_payer: true,
        settlement: SwaptionSettlement::Physical,
    };
    assert!(g.swaption_price(&d, &sw, 16).val() > 0.);
    let v = Cheyette1F::<f64>::explicit_local_variance_third_order(0., 0., 0.04, 0.04, 0., 0.);
    assert!((v - 0.04).abs() < 1e-12);
    let cms = LinearTsrCmsPricer { mean_reversion: 0. };
    let r = cms.cms_swaplet_rate(&d, 1., 1., &[1., 2., 3., 4., 5.], |_| 0.01);
    let ann = (2..=5).map(|i| d.discount_factor(i as f64)).sum::<f64>();
    let forward = (d.discount_factor(1.) - d.discount_factor(5.)) / ann;
    assert!(r > forward);
    let mr = Dual::<1>::variable(0.0, 0);
    let d_ad = d.to_ad_curve(
        &d.discounts
            .iter()
            .map(|&x| Dual::<1>::from_value(x))
            .collect::<Vec<_>>(),
    );
    let ad_rate = LinearTsrCmsPricer { mean_reversion: mr }.cms_swaplet_rate(
        &d_ad,
        1.,
        1.,
        &[1., 2., 3., 4., 5.],
        |_| Dual::from_value(0.01),
    );
    let h = 1e-5;
    let up = LinearTsrCmsPricer { mean_reversion: h }.cms_swaplet_rate(
        &d,
        1.,
        1.,
        &[1., 2., 3., 4., 5.],
        |_| 0.01,
    );
    let down = LinearTsrCmsPricer { mean_reversion: -h }.cms_swaplet_rate(
        &d,
        1.,
        1.,
        &[1., 2., 3., 4., 5.],
        |_| 0.01,
    );
    assert!((ad_rate.d[0] - (up - down) / (2. * h)).abs() < 1e-6);
}

#[test]
fn g2_swaption_dual_parameters_are_live() {
    let (d, _) = curve();
    let sw = EuropeanSwaption {
        notional: 1e6,
        strike: 0.02,
        expiry: 1.0,
        swap_times: (1..=6).map(|x| x as f64).collect(),
        is_payer: true,
        settlement: SwaptionSettlement::Physical,
    };
    let p = G2PlusPlus {
        a: Dual::variable(0.1, 0),
        sigma: Dual::variable(0.01, 1),
        b: Dual::variable(0.3, 2),
        eta: Dual::variable(0.015, 3),
        rho: Dual::variable(-0.7, 4),
    };
    let da = d.to_ad_curve(
        &d.discounts
            .iter()
            .map(|&x| Dual::<5>::from_value(x))
            .collect::<Vec<_>>(),
    );
    let ad = p.swaption_price(&da, &sw, 24);
    assert!(ad.d.iter().all(|x| x.abs() > 1e-10));
    let base = [0.1, 0.01, 0.3, 0.015, -0.7];
    for j in 0..5 {
        let mut up = base;
        let mut down = base;
        let h = 1e-6;
        up[j] += h;
        down[j] -= h;
        let mk = |x: [f64; 5]| G2PlusPlus {
            a: x[0],
            sigma: x[1],
            b: x[2],
            eta: x[3],
            rho: x[4],
        };
        let fd =
            (mk(up).swaption_price(&d, &sw, 24) - mk(down).swaption_price(&d, &sw, 24)) / (2.0 * h);
        assert!(
            (ad.d[j] - fd).abs() / fd.abs().max(1.0) < 1e-6,
            "j={j} ad={} fd={fd}",
            ad.d[j]
        );
    }
}
