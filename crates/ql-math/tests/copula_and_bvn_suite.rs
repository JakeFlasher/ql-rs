use ql_ad::{Dual, Real};
use ql_math::*;

#[test]
fn bvn_identities_and_bounds() {
    for &rho in &[-0.99, -0.5, 0.0, 0.5, 0.99] {
        let a = -0.4;
        let b = 0.7;
        let c = bivariate_normal_cdf(a, b, rho);
        let pa = a.norm_cdf();
        let pb = b.norm_cdf();
        assert!(c + 1e-12 >= (pa + pb - 1.0).max(0.0));
        assert!(c <= pa.min(pb) + 1e-12);
        assert!((c - bivariate_normal_cdf(b, a, rho)).abs() < 1e-12);
        let reflect = c + bivariate_normal_cdf(a, -b, -rho);
        assert!((reflect - pa).abs() < 1e-10, "{rho}: {reflect} {pa}");
        let origin = bivariate_normal_cdf(0.0, 0.0, rho);
        let expected = 0.25 + rho.asin() / (2.0 * std::f64::consts::PI);
        assert!(
            (origin - expected).abs() < 1e-10,
            "{rho}: {origin} {expected}"
        );
    }
}

#[test]
fn bvn_dual_derivatives_match_identities() {
    let a = Dual::<3>::variable(-0.3, 0);
    let b = Dual::<3>::variable(0.6, 1);
    let rho = Dual::<3>::variable(0.4, 2);
    let x = bivariate_normal_cdf(a, b, rho);
    let h = 1e-5;
    let fa = (bivariate_normal_cdf(-0.3 + h, 0.6, 0.4) - bivariate_normal_cdf(-0.3 - h, 0.6, 0.4))
        / (2.0 * h);
    let fb = (bivariate_normal_cdf(-0.3, 0.6 + h, 0.4) - bivariate_normal_cdf(-0.3, 0.6 - h, 0.4))
        / (2.0 * h);
    let fr = (bivariate_normal_cdf(-0.3, 0.6, 0.4 + h) - bivariate_normal_cdf(-0.3, 0.6, 0.4 - h))
        / (2.0 * h);
    let analytic_a =
        (-0.3_f64).norm_pdf() * ((0.6 - 0.4 * -0.3) / (1.0 - 0.4 * 0.4).sqrt()).norm_cdf();
    let analytic_b =
        0.6_f64.norm_pdf() * ((-0.3 - 0.4 * 0.6) / (1.0 - 0.4 * 0.4).sqrt()).norm_cdf();
    assert!((x.d[0] - fa).abs() < 1e-8);
    assert!((x.d[1] - fb).abs() < 1e-8);
    assert!((x.d[2] - fr).abs() < 1e-8);
    assert!((x.d[0] - analytic_a).abs() < 1e-10);
    assert!((x.d[1] - analytic_b).abs() < 1e-10);
    assert!((x.d[2] - bivariate_normal_pdf(-0.3, 0.6, 0.4)).abs() < 1e-10);
}

#[test]
fn copula_boundaries_rectangles_and_conditionals() {
    let copulas: Vec<Box<dyn Fn(f64, f64) -> f64>> = vec![
        Box::new(|u, v| GaussianCopula { rho: 0.4 }.cdf(u, v)),
        Box::new(|u, v| ClaytonCopula { theta: 1.5 }.cdf(u, v)),
        Box::new(|u, v| FrankCopula { theta: 2.0 }.cdf(u, v)),
        Box::new(|u, v| GumbelCopula { theta: 1.5 }.cdf(u, v)),
        Box::new(|u, v| FarlieGumbelMorgensternCopula { theta: 0.5 }.cdf(u, v)),
    ];
    for cdf in copulas {
        assert!(cdf(0.4, 0.0).abs() < 1e-12);
        assert!((cdf(0.4, 1.0) - 0.4).abs() < 1e-10);
        assert!((cdf(1.0, 0.4) - 0.4).abs() < 1e-10);
        let rect = cdf(0.7, 0.8) - cdf(0.3, 0.8) - cdf(0.7, 0.2) + cdf(0.3, 0.2);
        assert!(rect >= -1e-10);
    }
    let g = GaussianCopula {
        rho: Dual::<2>::variable(0.4, 0),
    };
    let u = Dual::<2>::variable(0.4, 1);
    let v = Dual::<2>::from_f64(0.6);
    let cond = g.conditional_cdf_v_given_u(u, v);
    assert!(cond.d[1].is_finite());
}

#[test]
fn bvn_high_correlation_diagonal_dual_is_finite_and_exact() {
    for &(av, bv, rv) in &[(0.5, 0.5, 0.995), (-0.3, 0.6, -0.96)] {
        let a = Dual::<3>::variable(av, 0);
        let b = Dual::<3>::variable(bv, 1);
        let rho = Dual::<3>::variable(rv, 2);
        let value = bivariate_normal_cdf(a, b, rho);
        assert!(value.val.is_finite());
        assert!(value.d.iter().all(|x| x.is_finite()));
        let den = (1.0 - rv * rv).sqrt();
        let analytic_a = av.norm_pdf() * ((bv - rv * av) / den).norm_cdf();
        let analytic_b = bv.norm_pdf() * ((av - rv * bv) / den).norm_cdf();
        let analytic_rho = bivariate_normal_pdf(av, bv, rv);
        assert!(
            (value.d[0] - analytic_a).abs() < 1e-10,
            "a={av} b={bv} rho={rv}: {} {}",
            value.d[0],
            analytic_a
        );
        assert!(
            (value.d[1] - analytic_b).abs() < 1e-10,
            "a={av} b={bv} rho={rv}: {} {}",
            value.d[1],
            analytic_b
        );
        assert!(
            (value.d[2] - analytic_rho).abs() < 1e-10,
            "a={av} b={bv} rho={rv}: {} {}",
            value.d[2],
            analytic_rho
        );
    }
}
