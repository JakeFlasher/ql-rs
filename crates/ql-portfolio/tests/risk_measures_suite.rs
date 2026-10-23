use ql_portfolio::*;

#[test]
fn acerbi_tasche_properties() {
    let x = vec![-5.0, -2.0, -1.0, 1.0, 2.0, 3.0];
    let y = vec![-3.0, -1.0, 0.0, 1.0, 2.0, 4.0];
    let xy: Vec<f64> = x.iter().zip(&y).map(|(a, b)| a + b).collect();
    let confidence = 0.75;
    let es_x = empirical_expected_shortfall(&x, confidence);
    let es_y = empirical_expected_shortfall(&y, confidence);
    let es_xy = empirical_expected_shortfall(&xy, confidence);
    assert!(es_xy <= es_x + es_y + 1e-12);
    let scaled: Vec<f64> = x.iter().map(|v| 2.0 * v).collect();
    assert!((empirical_expected_shortfall(&scaled, confidence) - 2.0 * es_x).abs() < 1e-12);
    let shifted: Vec<f64> = x.iter().map(|v| v + 3.0).collect();
    assert!((empirical_expected_shortfall(&shifted, confidence) - (es_x - 3.0)).abs() < 1e-12);
    assert!(es_x >= empirical_var(&x, confidence));
}

#[test]
fn gaussian_and_euler_allocations_match_fd() {
    let means = vec![0.01, 0.02, -0.005];
    let covariance = vec![
        vec![0.04, 0.01, 0.0],
        vec![0.01, 0.09, 0.005],
        vec![0.0, 0.005, 0.025],
    ];
    let weights = vec![0.4, 0.35, 0.25];
    let report = euler_gaussian_es_allocation_aad(&means, &covariance, &weights, 0.975);
    assert!((report.component_es.iter().sum::<f64>() - report.es).abs() < 1e-12);
    let h = 1e-6;
    for i in 0..weights.len() {
        let mut up = weights.clone();
        let mut down = weights.clone();
        up[i] += h;
        down[i] -= h;
        let fd = euler_gaussian_es_allocation_aad(&means, &covariance, &up, 0.975).es;
        let fm = euler_gaussian_es_allocation_aad(&means, &covariance, &down, 0.975).es;
        assert!((report.marginal_es[i] - (fd - fm) / (2.0 * h)).abs() < 1e-8);
    }
    let scenarios = vec![
        vec![-0.10, -0.05, 0.01],
        vec![-0.04, 0.02, -0.02],
        vec![0.01, -0.03, 0.02],
        vec![0.03, 0.04, -0.01],
        vec![0.06, 0.01, 0.03],
    ];
    let scenario = euler_scenario_es_allocation(&scenarios, &weights, 0.8);
    assert!((scenario.component_es.iter().sum::<f64>() - scenario.es).abs() < 1e-12);
    assert!(gaussian_expected_shortfall(0.0, 1.0, 0.975) > gaussian_var(0.0, 1.0, 0.975));
}
