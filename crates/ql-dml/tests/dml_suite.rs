use ql_dml::*;
#[test]
fn pca_and_regression() {
    let states = (0..20)
        .map(|i| vec![i as f64, (i as f64) * 2.])
        .collect::<Vec<_>>();
    let p = fit_classical_pca(&states);
    assert!(p.eigenvalues[0] > p.eigenvalues[1]);
    let d = states.iter().map(|x| vec![-x[0], x[1]]).collect::<Vec<_>>();
    let q = fit_differential_pca_non_central(&d);
    assert!(q.explained_ratio(1) > 0.99);
    let y = states
        .iter()
        .map(|x| x[0] * x[0] + x[1])
        .collect::<Vec<_>>();
    let r = DiffPolynomialRegressor::fit(&states, &y, &d, 2, 0., 1e-8);
    assert!((r.predict_value(&[2., 4.]) - 8.).abs() < 1.);
}

#[test]
fn differential_pca_finds_spread_risk_axis() {
    let states = (0..512)
        .map(|i| {
            let z = (i as f64 - 255.5) / 50.0;
            vec![z + 0.01 * (i as f64).sin(), z + 0.01 * (i as f64).cos()]
        })
        .collect::<Vec<_>>();
    let classical = fit_classical_pca(&states);
    let differentials = (0..512)
        .map(|i| {
            let scale = 0.5 + 0.1 * (i as f64).sin().abs();
            vec![-scale, scale]
        })
        .collect::<Vec<_>>();
    let differential = fit_differential_pca_non_central(&differentials);
    let common = [1.0 / 2.0_f64.sqrt(), 1.0 / 2.0_f64.sqrt()];
    let spread = [-1.0 / 2.0_f64.sqrt(), 1.0 / 2.0_f64.sqrt()];
    let common_alignment = (0..2)
        .map(|i| classical.eigenvectors[i][0] * common[i])
        .sum::<f64>()
        .abs();
    let spread_alignment = (0..2)
        .map(|i| differential.eigenvectors[i][0] * spread[i])
        .sum::<f64>()
        .abs();
    assert!(common_alignment > 0.999);
    assert!(spread_alignment > 0.999);
}
#[test]
fn cholesky_ad_matches_fd() {
    let c = vec![vec![1., 0.3], vec![0.3, 1.]];
    let l = cholesky_lower(&c).unwrap();
    let mut lb = vec![vec![0.; 2]; 2];
    lb[1][0] = 1.;
    let a = cholesky_adjoint(&c, &lb).unwrap();
    let h = 1e-5;
    let mut up = c.clone();
    let mut dn = c.clone();
    up[1][0] += h;
    up[0][1] += h;
    dn[1][0] -= h;
    dn[0][1] -= h;
    let lp = cholesky_lower(&up).unwrap();
    let lm = cholesky_lower(&dn).unwrap();
    let fd = (lp[1][0] - lm[1][0]) / (2. * h);
    assert!((fd - a[1][0]).abs() < 1e-5);
    assert!(l[0][0] > 0.);
}
#[test]
fn basket_corr_report() {
    let c = vec![vec![1., 0.2], vec![0.2, 1.]];
    let r = price_correlated_basket_call_aad(
        &[100., 100.],
        &[0.2, 0.2],
        &[0.5, 0.5],
        &c,
        0.02,
        1.,
        100.,
        2048,
        8,
        1e-6,
        3,
    )
    .unwrap();
    assert!(r.price > 0. && r.spot_deltas.len() == 2);
    assert!(r.corr_stderr[0][1] > 0.0);
    let h = 1e-4;
    let mut up = c.clone();
    let mut dn = c.clone();
    up[0][1] += h;
    up[1][0] += h;
    dn[0][1] -= h;
    dn[1][0] -= h;
    let pu = price_correlated_basket_call_aad(
        &[100., 100.],
        &[0.2, 0.2],
        &[0.5, 0.5],
        &up,
        0.02,
        1.,
        100.,
        2048,
        8,
        1e-6,
        3,
    )
    .unwrap()
    .price;
    let pd = price_correlated_basket_call_aad(
        &[100., 100.],
        &[0.2, 0.2],
        &[0.5, 0.5],
        &dn,
        0.02,
        1.,
        100.,
        2048,
        8,
        1e-6,
        3,
    )
    .unwrap()
    .price;
    assert!((r.corr_sensitivities[0][1] - (pu - pd) / (2. * h)).abs() < 1e-5);
}
