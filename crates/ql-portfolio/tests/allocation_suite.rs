use ql_portfolio::{
    equal_risk_contribution_weights, euler_risk_contributions, risk_budgeting_weights,
};

#[test]
fn risk_parity_matches_budget_contributions() {
    let cov = vec![vec![0.04, 0.01], vec![0.01, 0.09]];
    let w = equal_risk_contribution_weights(&cov);
    let (sigma, rc) = euler_risk_contributions(&cov, &w);
    assert!(sigma > 0.0);
    assert!((rc[0] - rc[1]).abs() < 1e-12);
    let custom = risk_budgeting_weights(&cov, &[0.25, 0.75]);
    let (_, rc2) = euler_risk_contributions(&cov, &custom);
    assert!((rc2[0] / (rc2[0] + rc2[1]) - 0.25).abs() < 1e-10);
}

#[test]
fn box_markowitz_budget_bounds_and_errors() {
    use ql_portfolio::{
        AllocationError, BoxConstrainedMarkowitzSpec, box_constrained_markowitz_weights,
        try_box_constrained_markowitz_weights,
    };
    let spec = BoxConstrainedMarkowitzSpec {
        covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
        expected_returns: vec![0.08, 0.04],
        risk_aversion: 3.0,
        lower_bounds: vec![0.1, 0.1],
        upper_bounds: vec![0.8, 0.8],
    };
    let w = box_constrained_markowitz_weights(&spec);
    assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    assert!(
        w.iter()
            .zip(&spec.lower_bounds)
            .all(|(x, l)| *x >= *l - 1e-12)
    );
    assert!(
        w.iter()
            .zip(&spec.upper_bounds)
            .all(|(x, u)| *x <= *u + 1e-12)
    );
    let bad = BoxConstrainedMarkowitzSpec {
        lower_bounds: vec![0.8, 0.8],
        upper_bounds: vec![1.0, 1.0],
        ..spec
    };
    assert_eq!(
        try_box_constrained_markowitz_weights(&bad),
        Err(AllocationError::InfeasibleBounds)
    );
}
