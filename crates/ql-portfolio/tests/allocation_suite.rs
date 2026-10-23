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
