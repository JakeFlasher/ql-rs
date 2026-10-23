#![feature(portable_simd)]
use ql_ad::Dual2;
use ql_math::*;

#[test]
fn brent_solve_ad_second_order_cross_derivative() {
    let solve_for = |t1: Dual2<2>, t2: Dual2<2>| {
        brent_solve_ad::<Dual2<2>, _>(|z| z * z * z + t1 * z - t2, 0.0, 3.0, 1e-13, 200).unwrap()
    };
    let x = solve_for(Dual2::variable(2.0, 0), Dual2::variable(5.0, 1));
    let x0 = x.val;
    assert!((x0 * x0 * x0 + 2.0 * x0 - 5.0).abs() < 1e-12);

    let h = 1e-5;
    let f_scalar =
        |a: f64, b: f64| brent_solve(|z| z * z * z + a * z - b, 0.0, 3.0, 1e-14, 200).unwrap();
    let d1_fd = (f_scalar(2.0 + h, 5.0) - f_scalar(2.0 - h, 5.0)) / (2.0 * h);
    let d2_fd = (f_scalar(2.0, 5.0 + h) - f_scalar(2.0, 5.0 - h)) / (2.0 * h);
    let cross_fd =
        (f_scalar(2.0 + h, 5.0 + h) - f_scalar(2.0 + h, 5.0 - h) - f_scalar(2.0 - h, 5.0 + h)
            + f_scalar(2.0 - h, 5.0 - h))
            / (4.0 * h * h);

    assert!((x.d1[0] - d1_fd).abs() < 1e-6);
    assert!((x.d1[1] - d2_fd).abs() < 1e-6);
    assert!((x.d2[0][1] - cross_fd).abs() < 1e-5);
}

#[test]
fn interpolators_and_distributions() {
    let x = vec![0.0, 1.0, 2.0, 3.0];
    let y = vec![0.0, 1.0, 4.0, 9.0];
    let m = MonotoneCubicInterpolator::new(x.clone(), y.clone());
    let n = NaturalCubicSpline::new(x, y);
    assert!((m.eval(1.5) - 2.25).abs() < 0.2);
    assert!((n.eval(2.0) - 4.0).abs() < 1e-12);
    assert!((inv_norm_cdf(0.975) - 1.959963984540054).abs() < 1e-12);
}

#[test]
fn simd_thomas_matches_8_scalar_systems() {
    let n = 16;
    let mut lo_simd = vec![std::simd::Simd::<f64, 8>::splat(0.0); n - 1];
    let mut di_simd = vec![std::simd::Simd::<f64, 8>::splat(0.0); n];
    let mut up_simd = vec![std::simd::Simd::<f64, 8>::splat(0.0); n - 1];
    let mut rhs_simd = vec![std::simd::Simd::<f64, 8>::splat(0.0); n];

    for i in 0..n {
        let mut d_arr = [0.0; 8];
        let mut r_arr = [0.0; 8];
        for lane in 0..8 {
            d_arr[lane] = 2.5 + 0.1 * lane as f64 + 0.02 * i as f64;
            r_arr[lane] = (i as f64 + 1.0) * (lane as f64 + 0.5);
        }
        di_simd[i] = std::simd::Simd::from_array(d_arr);
        rhs_simd[i] = std::simd::Simd::from_array(r_arr);
        if i + 1 < n {
            let mut l_arr = [0.0; 8];
            let mut u_arr = [0.0; 8];
            for lane in 0..8 {
                l_arr[lane] = -0.8 - 0.01 * lane as f64;
                u_arr[lane] = -0.7 - 0.02 * lane as f64;
            }
            lo_simd[i] = std::simd::Simd::from_array(l_arr);
            up_simd[i] = std::simd::Simd::from_array(u_arr);
        }
    }

    let mut out_simd = vec![std::simd::Simd::<f64, 8>::splat(0.0); n];
    let mut scratch_simd = vec![std::simd::Simd::<f64, 8>::splat(0.0); n];
    SimdThomasSolver::<8>::new(n).solve_in_place(
        &lo_simd,
        &di_simd,
        &up_simd,
        &rhs_simd,
        &mut out_simd,
        &mut scratch_simd,
    );

    for lane in 0..8 {
        let lo: Vec<f64> = lo_simd.iter().map(|v| v.to_array()[lane]).collect();
        let di: Vec<f64> = di_simd.iter().map(|v| v.to_array()[lane]).collect();
        let up: Vec<f64> = up_simd.iter().map(|v| v.to_array()[lane]).collect();
        let rhs: Vec<f64> = rhs_simd.iter().map(|v| v.to_array()[lane]).collect();
        let mut out = vec![0.0; n];
        let mut scratch = vec![0.0; n];
        solve_tridiagonal(&lo, &di, &up, &rhs, &mut out, &mut scratch);
        for i in 0..n {
            assert!((out_simd[i].to_array()[lane] - out[i]).abs() < 1e-13);
        }
    }
}
