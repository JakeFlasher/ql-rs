use ql_ad::*;
use ql_vol::*;

#[test]
fn implied_vol_first_and_second_derivatives_match_fd() {
    let f0 = 100.0;
    let k0 = 100.0;
    let t0 = 1.0;
    let p0 = black76_undiscounted(f0, k0, 0.2, t0, true);

    let f = Dual2::<3>::variable(f0, 0);
    let k = Dual2::<3>::variable(k0, 1);
    let p = Dual2::<3>::variable(p0, 2);
    let iv = implied_black_vol(f, k, Dual2::from_value(t0), p, true).unwrap();
    assert!((iv.val - 0.2).abs() < 1e-12);

    let h = 1e-4;
    let iv_p_plus = implied_black_vol(f0, k0, t0, p0 + h, true).unwrap();
    let iv_p_minus = implied_black_vol(f0, k0, t0, p0 - h, true).unwrap();
    let d_iv_dp = (iv_p_plus - iv_p_minus) / (2.0 * h);
    let d2_iv_dp2 = (iv_p_plus - 2.0 * iv.val + iv_p_minus) / (h * h);

    assert!((iv.d1[2] - d_iv_dp).abs() < 1e-6);
    assert!((iv.d2[2][2] - d2_iv_dp2).abs() < 1e-5);
}

#[test]
fn svi_and_sabr_sensitivities_match_fd() {
    let s = AdSviSmile::<f64>::new(100.0, 1.0, 0.04, 0.2, -0.3, 0.05, 0.4);
    let g = svi_quote_sensitivities(&s, 105.0);
    let h = 1e-6;
    let dm_fd = (AdSviSmile::<f64>::new(100.0, 1.0, 0.04, 0.2, -0.3, 0.05 + h, 0.4)
        .total_variance(105.0)
        - AdSviSmile::<f64>::new(100.0, 1.0, 0.04, 0.2, -0.3, 0.05 - h, 0.4).total_variance(105.0))
        / (2.0 * h);
    assert!((g[3] - dm_fd).abs() < 1e-7);

    let a = AdSabrSmile::<f64>::new(100.0, 1.0, 0.2, 0.5, -0.2, 0.4);
    let sg = sabr_quote_sensitivities(&a, 105.0);
    let d_alpha_fd = (AdSabrSmile::<f64>::new(100.0, 1.0, 0.2 + h, 0.5, -0.2, 0.4).vol(105.0)
        - AdSabrSmile::<f64>::new(100.0, 1.0, 0.2 - h, 0.5, -0.2, 0.4).vol(105.0))
        / (2.0 * h);
    assert!((sg[0] - d_alpha_fd).abs() < 1e-7);
}
