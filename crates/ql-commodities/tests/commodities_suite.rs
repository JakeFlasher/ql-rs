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
