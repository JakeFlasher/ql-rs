use ql_ad::*;
use ql_curves::DiscountCurve;
use ql_models::*;
use ql_time::Date;

fn complex_fn<T: Real>(re: T, im: T) -> Complex<T> {
    let z = Complex::new(re, im);
    (z.exp() + z.ln()) * z.sqrt() / (z + Complex::from_f64(0.5, -0.2))
}

#[test]
fn complex_ad_first_and_second_derivatives_match_fd() {
    let x0 = 1.2;
    let y0 = 0.4;
    let w = complex_fn(Dual2::<2>::variable(x0, 0), Dual2::<2>::variable(y0, 1));
    let h = 1e-5;

    let fx_re = (complex_fn(x0 + h, y0).re - complex_fn(x0 - h, y0).re) / (2.0 * h);
    let fy_re = (complex_fn(x0, y0 + h).re - complex_fn(x0, y0 - h).re) / (2.0 * h);
    let fxx_re = (complex_fn(x0 + h, y0).re - 2.0 * complex_fn(x0, y0).re
        + complex_fn(x0 - h, y0).re)
        / (h * h);
    let fxy_re = (complex_fn(x0 + h, y0 + h).re
        - complex_fn(x0 + h, y0 - h).re
        - complex_fn(x0 - h, y0 + h).re
        + complex_fn(x0 - h, y0 - h).re)
        / (4.0 * h * h);

    assert!((w.re.d1[0] - fx_re).abs() < 1e-8);
    assert!((w.re.d1[1] - fy_re).abs() < 1e-8);
    assert!((w.re.d2[0][0] - fxx_re).abs() < 1e-5);
    assert!((w.re.d2[0][1] - fxy_re).abs() < 1e-5);
}

#[test]
fn cos_heston_parity_and_dual5_sensitivities_match_fd() {
    let base = [0.04, 1.5, 0.04, 0.3, -0.7];
    let p_f64 = HestonParams {
        v0: base[0],
        kappa: base[1],
        theta: base[2],
        sigma_v: base[3],
        rho: base[4],
    };
    let e = CosHestonEngine::default();
    let c = e.price(&p_f64, 100.0, 100.0, 0.02, 0.01, 1.0, true);
    let put = e.price(&p_f64, 100.0, 100.0, 0.02, 0.01, 1.0, false);
    let parity = 100.0 * (-0.01f64).exp() - 100.0 * (-0.02f64).exp();
    assert!((c - put - parity).abs() < 1e-12);

    let p_dual = HestonParams {
        v0: Dual::<5>::variable(base[0], 0),
        kappa: Dual::<5>::variable(base[1], 1),
        theta: Dual::<5>::variable(base[2], 2),
        sigma_v: Dual::<5>::variable(base[3], 3),
        rho: Dual::<5>::variable(base[4], 4),
    };
    let c_ad = e.price(
        &p_dual,
        Dual::from_value(100.0),
        Dual::from_value(100.0),
        Dual::from_value(0.02),
        Dual::from_value(0.01),
        Dual::from_value(1.0),
        true,
    );
    assert!((c_ad.val - c).abs() < 1e-12);

    let h = 1e-6;
    for idx in 0..5 {
        let mut p_plus = base;
        let mut p_minus = base;
        p_plus[idx] += h;
        p_minus[idx] -= h;
        let c_p = e.price(
            &HestonParams {
                v0: p_plus[0],
                kappa: p_plus[1],
                theta: p_plus[2],
                sigma_v: p_plus[3],
                rho: p_plus[4],
            },
            100.0,
            100.0,
            0.02,
            0.01,
            1.0,
            true,
        );
        let c_m = e.price(
            &HestonParams {
                v0: p_minus[0],
                kappa: p_minus[1],
                theta: p_minus[2],
                sigma_v: p_minus[3],
                rho: p_minus[4],
            },
            100.0,
            100.0,
            0.02,
            0.01,
            1.0,
            true,
        );
        let fd = (c_p - c_m) / (2.0 * h);
        assert!(
            (c_ad.d[idx] - fd).abs() < 1e-5,
            "idx {idx}: ad={} fd={}",
            c_ad.d[idx],
            fd
        );
    }
}

#[test]
fn calibrate_heston_recovers_synthetic_surface() {
    let true_params = HestonParams {
        v0: 0.04,
        kappa: 1.8,
        theta: 0.045,
        sigma_v: 0.35,
        rho: -0.65,
    };
    let engine = CosHestonEngine::default();
    let mut quotes = Vec::new();
    for &expiry in &[0.5, 1.0] {
        for &strike in &[85.0, 95.0, 100.0, 105.0, 115.0] {
            let price = engine.price(&true_params, 100.0, strike, 0.03, 0.0, expiry, true);
            quotes.push(HestonMarketQuote {
                strike,
                expiry,
                market_price: price,
                is_call: true,
            });
        }
    }
    let guess = HestonParams {
        v0: 0.035,
        kappa: 1.4,
        theta: 0.04,
        sigma_v: 0.30,
        rho: -0.55,
    };
    let cal = calibrate_heston(100.0, 0.03, 0.0, guess, &quotes).unwrap();
    assert!(cal.rmse < 1e-4, "rmse={}", cal.rmse);
}

#[test]
fn hull_white_swaption_jamshidian_ad_matches_fd() {
    let ref_date = Date::from_ymd(2024, 1, 1).unwrap();
    let curve_f64 = DiscountCurve::new(
        ref_date,
        vec![1.0, 2.0, 3.0, 4.0],
        vec![
            (-0.03f64).exp(),
            (-0.06f64).exp(),
            (-0.09f64).exp(),
            (-0.12f64).exp(),
        ],
    );
    let curve_ad = DiscountCurve::new(
        ref_date,
        vec![1.0, 2.0, 3.0, 4.0],
        vec![
            Dual::<2>::from_value((-0.03f64).exp()),
            Dual::<2>::from_value((-0.06f64).exp()),
            Dual::<2>::from_value((-0.09f64).exp()),
            Dual::<2>::from_value((-0.12f64).exp()),
        ],
    );

    let a0 = 0.08;
    let sig0 = 0.012;
    let hw_ad = HullWhite1F {
        a: Dual::<2>::variable(a0, 0),
        sigma: Dual::<2>::variable(sig0, 1),
    };
    let pay_times = [2.0, 3.0, 4.0];
    let accruals = [1.0, 1.0, 1.0];
    let swaption_ad = hw_ad
        .swaption_jamshidian(
            &curve_ad,
            1.0,
            &pay_times,
            &accruals,
            Dual::from_value(0.03),
            true,
        )
        .unwrap();

    let h = 1e-6;
    let price_fn = |a: f64, sigma: f64| {
        HullWhite1F { a, sigma }
            .swaption_jamshidian(&curve_f64, 1.0, &pay_times, &accruals, 0.03, true)
            .unwrap()
    };
    let da_fd = (price_fn(a0 + h, sig0) - price_fn(a0 - h, sig0)) / (2.0 * h);
    let dsig_fd = (price_fn(a0, sig0 + h) - price_fn(a0, sig0 - h)) / (2.0 * h);

    assert!((swaption_ad.d[0] - da_fd).abs() < 1e-6);
    assert!((swaption_ad.d[1] - dsig_fd).abs() < 1e-6);
}

#[test]
fn lifted_heston_n1_matches_classical() {
    let p = HestonParams {
        v0: 0.04,
        kappa: 1.5,
        theta: 0.04,
        sigma_v: 0.3,
        rho: -0.7,
    };
    let l = LiftedHestonParams {
        v0: p.v0,
        lambda: p.kappa,
        theta: p.theta,
        nu: p.sigma_v,
        rho: p.rho,
        hurst: 0.1,
        n_factors: 1,
        r_n: 2.5,
    };
    let a = CosHestonEngine::default().price(&p, 100., 100., 0.02, 0., 1., true);
    let b = lifted_heston_cos_price(100., 100., 1., 0.02, 0., &l, 160, 12., 128, true);
    assert!((a - b).abs() < 1e-4);
}

#[test]
fn rough_heston_dual5_sensitivities_match_fd_and_short_dated_skew() {
    let base = [0.04, 1.5, 0.04, 0.3, -0.7];
    let make = |x: [f64; 5]| LiftedHestonParams {
        v0: x[0],
        lambda: x[1],
        theta: x[2],
        nu: x[3],
        rho: x[4],
        hurst: 0.1,
        n_factors: 8,
        r_n: 2.5,
    };
    let p = make(base);
    let call = lifted_heston_cos_price(100.0, 100.0, 1.0, 0.02, 0.0, &p, 128, 10.0, 64, true);
    let put = lifted_heston_cos_price(100.0, 100.0, 1.0, 0.02, 0.0, &p, 128, 10.0, 64, false);
    assert!(call > 0.0 && put > 0.0);
    assert!((call - put - (100.0 * (-0.0f64).exp() - 100.0 * (-0.02f64).exp())).abs() < 1e-10);
    let p1 = LiftedHestonParams { n_factors: 1, ..p };
    let classical = lifted_heston_cos_price(100.0, 100.0, 1.0, 0.02, 0.0, &p1, 128, 10.0, 64, true);
    assert!((call - classical).abs() > 1e-3);
    let pd = LiftedHestonParams {
        v0: Dual::<5>::variable(base[0], 0),
        lambda: Dual::<5>::variable(base[1], 1),
        theta: Dual::<5>::variable(base[2], 2),
        nu: Dual::<5>::variable(base[3], 3),
        rho: Dual::<5>::variable(base[4], 4),
        hurst: 0.1,
        n_factors: 8,
        r_n: 2.5,
    };
    let ad = lifted_heston_cos_price(
        Dual::from_value(100.0),
        Dual::from_value(100.0),
        Dual::from_value(1.0),
        Dual::from_value(0.02),
        Dual::from_value(0.0),
        &pd,
        128,
        10.0,
        64,
        true,
    );
    let h = 1e-6;
    for j in 0..5 {
        let mut xp = base;
        let mut xm = base;
        xp[j] += h;
        xm[j] -= h;
        let fp =
            lifted_heston_cos_price(100.0, 100.0, 1.0, 0.02, 0.0, &make(xp), 128, 10.0, 64, true);
        let fm =
            lifted_heston_cos_price(100.0, 100.0, 1.0, 0.02, 0.0, &make(xm), 128, 10.0, 64, true);
        assert!((ad.d[j] - (fp - fm) / (2.0 * h)).abs() < 1e-4, "j={j}");
    }
}

#[test]
fn phase12_rbergomi_and_slv_particle_invariants() {
    let p = ql_models::RoughBergomiParams {
        xi0: 0.04,
        eta: 1.2,
        rho: -0.7,
        hurst: 0.1,
    };
    let call =
        ql_models::rbergomi_hybrid_mc_price(100.0, 100.0, 1.0, 0.02, 0.0, &p, 16, 128, true, 42);
    let put =
        ql_models::rbergomi_hybrid_mc_price(100.0, 100.0, 1.0, 0.02, 0.0, &p, 16, 128, false, 42);
    let parity = call - put - (100.0 - 100.0 * (-0.02_f64).exp());
    assert!(call > 0.0 && put > 0.0 && parity.abs() < 1e-10);

    let heston = ql_models::HestonParams {
        v0: 0.04,
        kappa: 2.0,
        theta: 0.04,
        sigma_v: 0.3,
        rho: -0.5,
    };
    let slv = ql_models::calibrate_slv_particle_leverage(
        100.0,
        0.02,
        0.0,
        1.0,
        &heston,
        0.5,
        |_, _| 0.2,
        16,
        64,
        3,
        7,
    );
    let slv_px =
        ql_models::price_slv_european_mc(&slv, 100.0, 100.0, 0.02, 0.0, 1.0, 128, true, 11);
    assert!(
        slv.slices
            .iter()
            .flat_map(|slice| slice.leverage.iter())
            .all(|&l| l > 0.2 && l < 5.0)
    );
    assert!(slv_px > 5.0 && slv_px < 20.0);
}
