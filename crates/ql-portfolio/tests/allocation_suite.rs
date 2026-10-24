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

#[test]
fn long_only_markowitz_wrapper_is_budgeted() {
    let cov = vec![vec![0.04, 0.01], vec![0.01, 0.09]];
    let w = ql_portfolio::long_only_markowitz_weights(&cov, &[0.06, 0.03], 3.0);
    assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    assert!(w.iter().all(|x| *x >= -1e-12));
}

#[test]
fn phase12_mean_cvar_and_execution_invariants() {
    let cvar_spec = ql_portfolio::MeanCvarSpec {
        scenario_returns: vec![
            vec![-0.08, -0.01],
            vec![-0.06, -0.02],
            vec![0.02, 0.01],
            vec![0.05, 0.02],
            vec![0.06, 0.03],
        ],
        expected_returns: vec![0.07, 0.04],
        risk_aversion: 4.0,
        confidence: 0.6,
        smoothing_eps: 1e-3,
        lower_bounds: vec![0.0, 0.0],
        upper_bounds: vec![1.0, 1.0],
    };
    let res1 = ql_portfolio::mean_cvar_weights(&cvar_spec);
    let mut cvar_spec2 = cvar_spec.clone();
    cvar_spec2.smoothing_eps = 5e-2;
    let res2 = ql_portfolio::mean_cvar_weights(&cvar_spec2);
    assert!((res1.weights.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    assert!((res1.cvar - res2.cvar).abs() > 1e-4);

    let ow = ql_portfolio::obizhaeva_wang_execution(&ql_portfolio::ObizhaevaWangSpec {
        initial_inventory: 100.0,
        trade_times: vec![0.0, 0.25, 0.5, 0.75, 1.0],
        volatility: 0.2,
        transient_impact: 0.02,
        resilience_speed: 1.5,
        temporary_impact: 0.001,
        risk_aversion: 0.0,
    });
    assert!((ow.trade_schedule.iter().sum::<f64>() - 100.0).abs() < 1e-10);
    assert!(ow.inventory_schedule.last().unwrap().abs() < 1e-10);
    assert!(ow.trade_schedule[0] > ow.trade_schedule[1]);
    assert!(ow.trade_schedule[3] > ow.trade_schedule[1]);
}
