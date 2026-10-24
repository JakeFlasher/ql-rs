use ql_ad::{Dual, Real};
use ql_curves::FlatForwardCurve;
use ql_models::HullWhite1F;
use ql_time::Date;
use ql_xva::{
    ChebyshevInterpolant1D, ChebyshevTensor2D, CsaSpec, NettingSetSwap, NpvCube, SimmBucket,
    compute_bruteforce_dim_and_mva, compute_chebyshev_dim_and_mva, compute_isda_simm, compute_xva,
    simulate_hw1f_exposure_cube, simulate_hw1f_multi_netting_set_cube,
};

fn setup() -> (
    FlatForwardCurve<f64>,
    HullWhite1F<f64>,
    NettingSetSwap,
    NpvCube<f64>,
) {
    let ref_date = Date::from_ymd(2024, 1, 1).unwrap();
    let curve = FlatForwardCurve::new(ref_date, 0.02);
    let model = HullWhite1F {
        a: 0.1,
        sigma: 0.01,
    };
    let swap = NettingSetSwap {
        notional: 100_000.0,
        fixed_rate: 0.0202,
        payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
        accruals: vec![1.0; 5],
        is_payer: true,
    };
    let cube = simulate_hw1f_exposure_cube(
        &model,
        &curve,
        std::slice::from_ref(&swap),
        &[1.0, 2.0, 3.0, 4.0],
        1_024,
        42,
    );
    (curve, model, swap, cube)
}

#[test]
fn multi_netting_set_cube_keeps_legal_sets_separate() {
    let (curve, model, swap, _) = setup();
    let other = NettingSetSwap {
        notional: -50_000.0,
        ..swap.clone()
    };
    let cube = simulate_hw1f_multi_netting_set_cube(
        &model,
        &curve,
        &[vec![swap], vec![other]],
        &[1.0, 2.0],
        64,
        7,
    );
    assert_eq!(cube.n_netting_sets, 2);
    assert_ne!(cube.get(0, 0, 0), cube.get(1, 0, 0));
}

fn lift_cube(cube: &NpvCube<f64>) -> NpvCube<Dual<3>> {
    NpvCube {
        n_netting_sets: cube.n_netting_sets,
        times: cube.times.clone(),
        n_paths: cube.n_paths,
        data: cube.data.iter().map(|x| Dual::from_value(*x)).collect(),
    }
}

fn xva_f64(
    cube: &NpvCube<f64>,
    curve: &FlatForwardCurve<f64>,
    hc: f64,
    hs: f64,
    fs: f64,
) -> ql_xva::XvaReport<f64> {
    compute_xva(cube, curve, None, hc, 0.6, hs, 0.6, fs)
}

#[test]
fn hw1f_exposure_cube_csa_and_xva_sensitivities() {
    let (curve, _model, _swap, cube) = setup();
    let raw = compute_xva(&cube, &curve, None, 0.02, 0.6, 0.01, 0.6, 0.005);
    for k in 0..cube.times.len() {
        assert!(raw.ee[k] > 0.0);
        assert!(raw.ene[k] < 0.0);
        assert!(raw.pfe_95[k] > raw.ee[k]);
    }
    let csa = CsaSpec {
        threshold_cpty: 500.0,
        threshold_self: 500.0,
        mta: 50.0,
    };
    let collateral = compute_xva(&cube, &curve, Some(&csa), 0.02, 0.6, 0.01, 0.6, 0.005);
    assert!(collateral.epe < raw.epe);
    assert!(collateral.cva < raw.cva);

    let dual_cube = lift_cube(&cube);
    let dual_curve = FlatForwardCurve::new(curve.reference_date, Dual::from_value(0.02));
    let dual = compute_xva(
        &dual_cube,
        &dual_curve,
        None,
        Dual::variable(0.02, 0),
        Dual::from_value(0.6),
        Dual::variable(0.01, 1),
        Dual::from_value(0.6),
        Dual::variable(0.005, 2),
    );
    let h = 1e-6;
    let hc_fd = (xva_f64(&cube, &curve, 0.02 + h, 0.01, 0.005).cva
        - xva_f64(&cube, &curve, 0.02 - h, 0.01, 0.005).cva)
        / (2.0 * h);
    let hs_fd = (xva_f64(&cube, &curve, 0.02, 0.01 + h, 0.005).cva
        - xva_f64(&cube, &curve, 0.02, 0.01 - h, 0.005).cva)
        / (2.0 * h);
    let dva_hs_fd = (xva_f64(&cube, &curve, 0.02, 0.01 + h, 0.005).dva
        - xva_f64(&cube, &curve, 0.02, 0.01 - h, 0.005).dva)
        / (2.0 * h);
    let fva_fs_fd = (xva_f64(&cube, &curve, 0.02, 0.01, 0.005 + h).fva
        - xva_f64(&cube, &curve, 0.02, 0.01, 0.005 - h).fva)
        / (2.0 * h);
    assert!((dual.cva.d[0] - hc_fd).abs() < 1e-6);
    assert!((dual.cva.d[1] - hs_fd).abs() < 1e-6);
    assert!((dual.dva.d[1] - dva_hs_fd).abs() < 1e-6);
    assert!((dual.fva.d[2] - fva_fs_fd).abs() < 1e-6);

    // Also verify end-to-end AD sensitivities through simulate_hw1f_exposure_cube wrt (a, sigma)
    let model_dual = HullWhite1F {
        a: Dual::<2>::variable(0.1, 0),
        sigma: Dual::<2>::variable(0.01, 1),
    };
    let curve_dual2 = FlatForwardCurve::new(curve.reference_date, Dual::<2>::from_value(0.02));
    let cube_dual2 = simulate_hw1f_exposure_cube(
        &model_dual,
        &curve_dual2,
        std::slice::from_ref(&_swap),
        &[1.0, 2.0, 3.0, 4.0],
        256,
        42,
    );
    let xva_dual2 = compute_xva(
        &cube_dual2,
        &curve_dual2,
        None,
        Dual::from_value(0.02),
        Dual::from_value(0.6),
        Dual::from_value(0.01),
        Dual::from_value(0.6),
        Dual::from_value(0.005),
    );
    let eval_hw = |a_val: f64, sig_val: f64| {
        let m = HullWhite1F {
            a: a_val,
            sigma: sig_val,
        };
        let c = simulate_hw1f_exposure_cube(
            &m,
            &curve,
            std::slice::from_ref(&_swap),
            &[1.0, 2.0, 3.0, 4.0],
            256,
            42,
        );
        compute_xva(&c, &curve, None, 0.02, 0.6, 0.01, 0.6, 0.005).cva
    };
    let cva_da_fd = (eval_hw(0.1 + h, 0.01) - eval_hw(0.1 - h, 0.01)) / (2.0 * h);
    let cva_dsig_fd = (eval_hw(0.1, 0.01 + h) - eval_hw(0.1, 0.01 - h)) / (2.0 * h);
    assert!((xva_dual2.cva.d[0] - cva_da_fd).abs() < 1e-4);
    assert!((xva_dual2.cva.d[1] - cva_dsig_fd).abs() < 1e-3);
}

#[test]
fn isda_simm_dual_sensitivities_match_fd() {
    let b0 = SimmBucket {
        sensitivities: vec![Dual::<4>::variable(1.0, 0), Dual::<4>::variable(-2.0, 1)],
        risk_weights: vec![0.5, 0.4],
        concentration_threshold: 1.0,
        intra_corr: vec![vec![1.0, 0.2], vec![0.2, 1.0]],
    };
    let b1 = SimmBucket {
        sensitivities: vec![Dual::<4>::variable(0.5, 2), Dual::<4>::variable(-1.5, 3)],
        risk_weights: vec![0.3, 0.6],
        concentration_threshold: 1.0,
        intra_corr: vec![vec![1.0, 0.1], vec![0.1, 1.0]],
    };
    let value = compute_isda_simm(&[b0, b1], &[vec![1.0, 0.25], vec![0.25, 1.0]]);
    let h = 1e-6;
    let base = [1.0, -2.0, 0.5, -1.5];
    let fd = |idx: usize, sign: f64| {
        let mut x = base;
        x[idx] += sign * h;
        let a = SimmBucket {
            sensitivities: vec![x[0], x[1]],
            risk_weights: vec![0.5, 0.4],
            concentration_threshold: 1.0,
            intra_corr: vec![vec![1.0, 0.2], vec![0.2, 1.0]],
        };
        let c = SimmBucket {
            sensitivities: vec![x[2], x[3]],
            risk_weights: vec![0.3, 0.6],
            concentration_threshold: 1.0,
            intra_corr: vec![vec![1.0, 0.1], vec![0.1, 1.0]],
        };
        compute_isda_simm(&[a, c], &[vec![1.0, 0.25], vec![0.25, 1.0]])
    };
    for (i, ad) in value.d.iter().enumerate() {
        let d = (fd(i, 1.0) - fd(i, -1.0)) / (2.0 * h);
        assert!((*ad - d).abs() < 1e-7, "i={i} ad={} fd={d}", ad);
    }
}

#[test]
fn chebyshev_interpolation_and_dim_parity() {
    let c = ChebyshevInterpolant1D::from_fn(-1., 1., 8, |x| x * x + 2. * x + 1.);
    assert!((c.eval_barycentric(0.3) - 1.69).abs() < 1e-10);
    let a = ql_ad::Dual::<1>::variable(0.2, 0);
    let v = c.eval_clenshaw(a);
    assert!((v.val() - 1.44).abs() < 1e-8);
    assert!((v.d[0] - 2.4).abs() < 1e-7);
    let t = ChebyshevTensor2D::from_fn(-1., 1., -1., 1., 6, |x, y| x * x + y * y);
    assert!((t.eval_barycentric(0.2, -0.3) - 0.13).abs() < 1e-8);
    let tx = ql_ad::Dual::<2>::variable(0.2, 0);
    let ty = ql_ad::Dual::<2>::variable(-0.3, 1);
    let tv = ChebyshevTensor2D::from_fn(-1., 1., -1., 1., 6, |x, y| 1. + x * x + 2. * y * y)
        .eval_clenshaw(tx, ty);
    assert!((tv.val() - 1.22).abs() < 1e-8);
    assert!((tv.d[0] - 0.4).abs() < 1e-8);
    assert!((tv.d[1] + 1.2).abs() < 1e-8);
    let d = Date::from_ymd(2024, 1, 1).unwrap();
    let curve = FlatForwardCurve::new(d, 0.02);
    let model = HullWhite1F {
        a: 0.05,
        sigma: 0.01,
    };
    let swap = NettingSetSwap {
        notional: 1e6,
        fixed_rate: 0.02,
        payment_times: (1..=5).map(|x| x as f64).collect(),
        accruals: vec![1.; 5],
        is_payer: true,
    };
    let times = (1..=4).map(|x| x as f64).collect::<Vec<_>>();
    let brute = compute_bruteforce_dim_and_mva(
        &model,
        &curve,
        std::slice::from_ref(&swap),
        &times,
        32,
        6,
        0.01,
        11,
    );
    let fast = compute_chebyshev_dim_and_mva(
        &model,
        &curve,
        std::slice::from_ref(&swap),
        &times,
        32,
        6,
        0.01,
        11,
    );
    assert!(brute.mva > 0.0 && fast.mva > 0.0);
    assert!((brute.mva - fast.mva).abs() / brute.mva < 1e-4);
}

#[test]
fn sa_ccr_empty_classes_and_forward_starting_ir() {
    use ql_ad::Dual;
    use ql_xva::{SaCcrAssetClass, SaCcrNettingSet, SaCcrTrade, compute_sa_ccr_ead};
    let ns = SaCcrNettingSet {
        trades: vec![SaCcrTrade {
            mtm: Dual::<2>::variable(10.0, 0),
            notional: Dual::from_value(1_000_000.0),
            maturity: 5.0,
            asset_class: SaCcrAssetClass::InterestRate {
                ccy: "USD".into(),
                start_time: 2.0,
                end_time: 7.0,
            },
            is_long_primary: true,
            option_spec: None,
        }],
        collateral: Dual::from_value(0.0),
        is_margined: false,
        threshold: Dual::from_value(0.0),
        mta: Dual::from_value(0.0),
        nica: Dual::from_value(0.0),
        mpor_years: 0.5,
    };
    let report = compute_sa_ccr_ead(&ns);
    assert!(report.ead.val.is_finite() && report.ead.d.iter().all(|x| x.is_finite()));
    let empty = SaCcrNettingSet {
        trades: Vec::new(),
        ..ns
    };
    let r0 = compute_sa_ccr_ead(&empty);
    assert!(r0.ead.val.is_finite() && r0.ead.d.iter().all(|x| x.is_finite()));
}

#[test]
fn wwr_zero_preserves_ad_and_epe() {
    use ql_ad::Dual;
    let mut cube = NpvCube::<Dual<2>>::new(1, vec![0.0, 1.0, 2.0], 2);
    for (t, vals) in [(0, [0.0, 0.0]), (1, [1.0, 2.0]), (2, [2.0, 3.0])] {
        for (p, v) in vals.into_iter().enumerate() {
            cube.set(0, t, p, Dual::from_value(v));
        }
    }
    let curve = FlatForwardCurve::new(Date(0), Dual::from_value(0.02));
    let report = ql_xva::compute_xva_with_wwr(
        &cube,
        &curve,
        None,
        Dual::from_value(0.02),
        Dual::from_value(0.6),
        Dual::from_value(0.01),
        Dual::from_value(0.5),
        Dual::from_value(0.01),
        Dual::variable(0.0, 0),
    );
    assert!(report.epe.val > 0.0 && report.cva.val.is_finite() && report.cva.d[0].is_finite());
}
