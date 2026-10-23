use ql_ad::{Dual, Real};
use ql_instruments::OptionType;
use ql_levy::{
    CgmyProcess, KouJumpDiffusion, LevyModel, MertonJumpDiffusion, VarianceGamma, cos_levy_price,
};
use ql_vol::black_scholes_price;

#[test]
fn merton_cos_matches_poisson_series() {
    let m = MertonJumpDiffusion {
        sigma: 0.2,
        lambda_j: 0.3,
        mu_j: -0.1,
        delta_j: 0.2,
    };
    for &k in &[80.0, 100.0, 120.0] {
        for kind in [OptionType::Call, OptionType::Put] {
            let series = m.merton_series_price(100.0, k, 0.03, 0.0, 1.0, kind, 40);
            let cos = cos_levy_price(&m, 100.0, k, 0.03, 0.0, 1.0, kind, 256);
            assert!((series - cos).abs() < 1e-8, "{k} {kind:?}: {series} {cos}");
        }
    }
}

fn assert_cumulants_from_characteristic<M: LevyModel<f64>>(model: M) {
    let r = 0.03;
    let q = 0.01;
    let t = 1.7;
    let h = 1e-4;
    let (rp, ip) = model.char_exponent(h);
    let (rm, im) = model.char_exponent(-h);
    let (r0, _) = model.char_exponent(0.0);
    let (c1, c2, _) = model.cumulants(r, q, t);
    let c1_fd = (r - q + model.martingale_compensator() + (ip - im) / (2.0 * h)) * t;
    let c2_fd = -(rp - 2.0 * r0 + rm) / (h * h) * t;
    assert!((c1 - c1_fd).abs() < 1e-6, "c1 {c1} {c1_fd}");
    assert!((c2 - c2_fd).abs() < 1e-5, "c2 {c2} {c2_fd}");
}

#[test]
fn cumulants_match_char_exponent_derivatives() {
    assert_cumulants_from_characteristic(MertonJumpDiffusion {
        sigma: 0.2,
        lambda_j: 0.3,
        mu_j: -0.1,
        delta_j: 0.2,
    });
    assert_cumulants_from_characteristic(KouJumpDiffusion {
        sigma: 0.2,
        lambda_j: 0.2,
        p_up: 0.4,
        eta1: 5.0,
        eta2: 8.0,
    });
    assert_cumulants_from_characteristic(VarianceGamma {
        sigma: 0.18,
        nu: 0.2,
        theta: -0.1,
    });
    assert_cumulants_from_characteristic(CgmyProcess {
        sigma: 0.1,
        c: 0.5,
        g: 5.0,
        m: 10.0,
        y: 0.7,
    });
}

#[test]
fn zero_jump_models_reduce_to_black_scholes() {
    let reference = black_scholes_price(100.0, 100.0, 0.2, 1.0, 0.03, 0.0, true);
    let m = MertonJumpDiffusion {
        sigma: 0.2,
        lambda_j: 0.0,
        mu_j: -0.1,
        delta_j: 0.2,
    };
    let k = KouJumpDiffusion {
        sigma: 0.2,
        lambda_j: 0.0,
        p_up: 0.4,
        eta1: 5.0,
        eta2: 8.0,
    };
    let vg = VarianceGamma {
        sigma: 0.2,
        nu: 1e-6,
        theta: 0.0,
    };
    for value in [
        cos_levy_price(&m, 100.0, 100.0, 0.03, 0.0, 1.0, OptionType::Call, 256),
        cos_levy_price(&k, 100.0, 100.0, 0.03, 0.0, 1.0, OptionType::Call, 256),
        cos_levy_price(&vg, 100.0, 100.0, 0.03, 0.0, 1.0, OptionType::Call, 256),
    ] {
        assert!((value - reference).abs() < 1e-6, "{value} vs {reference}");
    }
}

fn fd<M: LevyModel<f64> + Copy>(mk: impl Fn(f64) -> M, base: f64) -> f64 {
    let h = 1e-5;
    let up = cos_levy_price(
        &mk(base + h),
        100.0,
        100.0,
        0.03,
        0.0,
        1.0,
        OptionType::Call,
        192,
    );
    let dn = cos_levy_price(
        &mk(base - h),
        100.0,
        100.0,
        0.03,
        0.0,
        1.0,
        OptionType::Call,
        192,
    );
    (up - dn) / (2.0 * h)
}

#[test]
fn levy_dual_sensitivities_match_fd() {
    let vg = VarianceGamma {
        sigma: Dual::<5>::variable(0.18, 0),
        nu: Dual::<5>::variable(0.2, 1),
        theta: Dual::<5>::variable(-0.1, 2),
    };
    let x = cos_levy_price(
        &vg,
        Dual::<5>::variable(100.0, 3),
        Dual::from_f64(100.0),
        Dual::<5>::variable(0.03, 4),
        Dual::zero(),
        Dual::one(),
        OptionType::Call,
        192,
    );
    let vg0 = VarianceGamma {
        sigma: 0.18,
        nu: 0.2,
        theta: -0.1,
    };
    let fds = [
        fd(|v| VarianceGamma { sigma: v, ..vg0 }, 0.18),
        fd(|v| VarianceGamma { nu: v, ..vg0 }, 0.2),
        fd(|v| VarianceGamma { theta: v, ..vg0 }, -0.1),
    ];
    for i in 0..3 {
        assert!((x.d[i] - fds[i]).abs() < 1e-5, "{} {}", x.d[i], fds[i]);
    }
    let kd = KouJumpDiffusion {
        sigma: Dual::<5>::variable(0.18, 0),
        lambda_j: Dual::<5>::variable(0.2, 1),
        p_up: Dual::<5>::variable(0.4, 2),
        eta1: Dual::<5>::variable(5.0, 3),
        eta2: Dual::<5>::variable(8.0, 4),
    };
    let kprice = cos_levy_price(
        &kd,
        Dual::from_f64(100.0),
        Dual::from_f64(100.0),
        Dual::from_f64(0.03),
        Dual::zero(),
        Dual::one(),
        OptionType::Call,
        128,
    );
    let kv = [0.18, 0.2, 0.4, 5.0, 8.0];
    for j in 0..5 {
        let h = 1e-5;
        let mut up = kv;
        let mut dn = kv;
        up[j] += h;
        dn[j] -= h;
        let mup = KouJumpDiffusion {
            sigma: up[0],
            lambda_j: up[1],
            p_up: up[2],
            eta1: up[3],
            eta2: up[4],
        };
        let mdn = KouJumpDiffusion {
            sigma: dn[0],
            lambda_j: dn[1],
            p_up: dn[2],
            eta1: dn[3],
            eta2: dn[4],
        };
        let fd = (cos_levy_price(&mup, 100.0, 100.0, 0.03, 0.0, 1.0, OptionType::Call, 128)
            - cos_levy_price(&mdn, 100.0, 100.0, 0.03, 0.0, 1.0, OptionType::Call, 128))
            / (2.0 * h);
        assert!(
            (kprice.d[j] - fd).abs() < 1e-5,
            "Kou {j}: {} {fd}",
            kprice.d[j]
        );
    }
    let cg = CgmyProcess {
        sigma: Dual::<5>::variable(0.1, 0),
        c: Dual::<5>::variable(0.5, 1),
        g: Dual::<5>::variable(5.0, 2),
        m: Dual::<5>::variable(10.0, 3),
        y: Dual::<5>::variable(0.7, 4),
    };
    let cprice = cos_levy_price(
        &cg,
        Dual::from_f64(100.0),
        Dual::from_f64(100.0),
        Dual::from_f64(0.03),
        Dual::zero(),
        Dual::one(),
        OptionType::Call,
        128,
    );
    assert!(cprice.val().is_finite());
    assert!(cprice.d.iter().all(|x| x.is_finite()));
    let base = CgmyProcess {
        sigma: 0.1,
        c: 0.5,
        g: 5.0,
        m: 10.0,
        y: 0.7,
    };
    let vals = [base.sigma, base.c, base.g, base.m, base.y];
    for j in 0..5 {
        let h = 1e-5;
        let mut up = vals;
        let mut dn = vals;
        up[j] += h;
        dn[j] -= h;
        let mup = CgmyProcess {
            sigma: up[0],
            c: up[1],
            g: up[2],
            m: up[3],
            y: up[4],
        };
        let mdn = CgmyProcess {
            sigma: dn[0],
            c: dn[1],
            g: dn[2],
            m: dn[3],
            y: dn[4],
        };
        let fd = (cos_levy_price(&mup, 100.0, 100.0, 0.03, 0.0, 1.0, OptionType::Call, 128)
            - cos_levy_price(&mdn, 100.0, 100.0, 0.03, 0.0, 1.0, OptionType::Call, 128))
            / (2.0 * h);
        assert!(
            (cprice.d[j] - fd).abs() < 1e-4,
            "CGMY {j}: {} {fd}",
            cprice.d[j]
        );
    }
}
