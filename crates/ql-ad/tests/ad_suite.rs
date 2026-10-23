#![feature(portable_simd)]

use core::simd::Simd;
use ql_ad::*;

#[test]
fn tail_safe_erfc_and_normal_cdf_match_scalar() {
    let erfc_scalar = Real::erfc(6.0_f64);
    let erfc_dual = Dual::<1>::from_value(6.0).erfc();
    assert!(erfc_scalar > 2e-17);
    assert!((erfc_scalar - erfc_dual.val).abs() < 1e-30);

    let cdf_scalar = (-8.0_f64).norm_cdf();
    let cdf_dual = Dual::<1>::from_value(-8.0).norm_cdf();
    assert!(cdf_scalar > 6e-16);
    assert!((cdf_scalar - cdf_dual.val).abs() < 1e-30);
}

fn composite_fn<T: Real>(x: T) -> T {
    x.exp()
        + x.ln()
        + x.sin()
        + x.cos()
        + x.tan()
        + x.atan()
        + x.sinh()
        + x.cosh()
        + x.tanh()
        + x.powf(1.7)
        + x.erf()
        + x.norm_cdf()
        + x.norm_pdf()
}

#[test]
fn dual_and_tape_match_central_fd() {
    let x0 = 1.2;
    let x_dual = Dual::<4>::variable(x0, 0);
    let y_dual = composite_fn(x_dual);

    TapeScope::start_recording();
    let x_tape = TapeScope::with_active(|t| t.input(x0));
    let y_tape = composite_fn(x_tape);
    let grad = TapeScope::gradient_from(y_tape);
    let d_tape = grad.gradient(x_tape);
    TapeScope::stop_recording();

    let h = 1e-6;
    let fd = (composite_fn(x0 + h) - composite_fn(x0 - h)) / (2.0 * h);

    assert!((y_dual.val - y_tape.val).abs() < 1e-14);
    assert!((y_dual.d[0] - d_tape).abs() < 1e-12);
    assert!((y_dual.d[0] - fd).abs() < 1e-6);

    let s = Dual::<1>::variable(0.2, 0).smooth_heaviside(Dual::<1>::zero(), 1.0);
    assert!((s.d[0] - 0.72).abs() < 1e-12);
}

#[test]
fn dual2_hessian_matches_fd() {
    fn f3<T: Real>(x: T, y: T, z: T) -> T {
        (x * y).exp() * (y * z).sin() + (x - z).norm_cdf()
    }
    let pt = [1.2, 0.7, -0.4];
    let x = Dual2::<3>::variable(pt[0], 0);
    let y = Dual2::<3>::variable(pt[1], 1);
    let z = Dual2::<3>::variable(pt[2], 2);
    let out = f3(x, y, z);

    let h = 1e-5;
    for i in 0..3 {
        let mut p = pt;
        let mut m = pt;
        p[i] += h;
        m[i] -= h;
        let g_fd = (f3(p[0], p[1], p[2]) - f3(m[0], m[1], m[2])) / (2.0 * h);
        assert!((out.d1[i] - g_fd).abs() < 1e-6);
        for j in 0..3 {
            let mut pp = pt;
            let mut pm = pt;
            let mut mp = pt;
            let mut mm = pt;
            pp[i] += h;
            pp[j] += h;
            pm[i] += h;
            pm[j] -= h;
            mp[i] -= h;
            mp[j] += h;
            mm[i] -= h;
            mm[j] -= h;
            let h_fd =
                (f3(pp[0], pp[1], pp[2]) - f3(pm[0], pm[1], pm[2]) - f3(mp[0], mp[1], mp[2])
                    + f3(mm[0], mm[1], mm[2]))
                    / (4.0 * h * h);
            assert!((out.d2[i][j] - h_fd).abs() < 1e-4);
        }
    }
}

#[test]
fn tape_checkpoint_mark_and_rewind() {
    let mut tape = ArenaTape::<f64>::new();
    let r = tape.input(0.05);
    let neg_one = tape.constant(-1.0);
    let neg_r = tape.mul(r, neg_one);
    let df = tape.exp(neg_r);
    tape.set_mark();
    let mark = tape.mark();

    let mut mark_adjoints = vec![0.0; mark];
    for path in 1..=4 {
        let c = tape.constant(path as f64 * 10.0);
        let payoff = tape.mul(df, c);
        let mut adj = vec![0.0; tape.len()];
        adj[payoff.idx as usize] = 1.0;
        tape.propagate_to_mark(&mut adj);
        for i in 0..mark {
            mark_adjoints[i] += adj[i];
        }
        tape.rewind_to_mark();
        assert_eq!(tape.len(), mark);
    }
    tape.propagate_mark_to_start(&mut mark_adjoints);
    let expected = -(1.0 + 2.0 + 3.0 + 4.0) * 10.0 * (-0.05_f64).exp();
    assert!((mark_adjoints[r.idx as usize] - expected).abs() < 1e-12);
}

#[test]
fn simd_exp_matches_scalar() {
    let x = Simd::<f64, 4>::from_array([-2.0, -0.1, 1.0, 5.0]);
    let y = <Simd<f64, 4> as SimdReal>::exp(x).to_array();
    for (a, b) in y.into_iter().zip([-2.0, -0.1, 1.0, 5.0]) {
        assert!((a - b.exp()).abs() < 1e-12);
    }
}
