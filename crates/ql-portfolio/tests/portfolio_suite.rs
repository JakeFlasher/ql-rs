use ql_ad::TapeScope;
use ql_models::HestonParams;
use ql_portfolio::{MarketSnapshot, PortfolioDag, PortfolioTrade};
use ql_scripting::{BytecodeProgram, ScriptedProduct};
use ql_time::Date;

fn snapshot() -> MarketSnapshot {
    MarketSnapshot {
        ref_date: Date::from_ymd(2024, 1, 1).unwrap(),
        ois_quotes: [0.02, 0.021, 0.022],
        fwd_quotes: [0.025, 0.026, 0.027],
        inflation_quotes: vec![("CPI".to_string(), 0.023)],
        heston_params: HestonParams {
            v0: 0.04,
            kappa: 1.5,
            theta: 0.04,
            sigma_v: 0.3,
            rho: -0.7,
        },
        eq_spot: 100.0,
        eq_vol: 0.2,
    }
}

fn trades() -> Vec<PortfolioTrade> {
    let note = BytecodeProgram::compile(&ScriptedProduct::asian_call(vec![0.5, 1.0], 100.0, 1.0));
    vec![
        PortfolioTrade::IrsSwap {
            id: 1,
            notional: 1e6,
            fixed_rate: 0.02,
            maturity_years: 5,
        },
        PortfolioTrade::IrsSwap {
            id: 2,
            notional: 2e6,
            fixed_rate: 0.025,
            maturity_years: 2,
        },
        PortfolioTrade::HestonEuropeanOption {
            id: 3,
            strike: 100.0,
            expiry: 1.0,
            is_call: true,
            notional: 1.0,
        },
        PortfolioTrade::HestonEuropeanOption {
            id: 4,
            strike: 110.0,
            expiry: 2.0,
            is_call: false,
            notional: 1.0,
        },
        PortfolioTrade::ScriptedNote {
            id: 5,
            program: note,
            expiry: 1.0,
        },
    ]
}

#[test]
fn incremental_cache_skips_clean_subgraphs_and_aggregates_risk() {
    let mut s = snapshot();
    let ts = trades();
    let mut dag = PortfolioDag::new();

    let first = dag.evaluate_portfolio(&s, &ts);
    assert_eq!(first.curve_recomputes, 1);
    assert_eq!(first.trade_recomputes, 5);
    assert_eq!(first.total_par_dv01.len(), 6);
    assert!(first.total_eq_delta > 0.0);

    let second = dag.evaluate_portfolio(&s, &ts);
    assert_eq!(second.curve_recomputes, 0);
    assert_eq!(second.trade_recomputes, 0);

    s.eq_spot = 102.0;
    let third = dag.evaluate_portfolio(&s, &ts);
    assert_eq!(third.curve_recomputes, 0);
    assert_eq!(third.trade_recomputes, 3);

    s.fwd_quotes[2] += 1e-4;
    let fourth = dag.evaluate_portfolio(&s, &ts);
    assert_eq!(fourth.curve_recomputes, 1);
    assert_eq!(fourth.trade_recomputes, 5);
    let predicted = third.total_par_dv01[5] * 1e-4;
    let observed = fourth.total_npv - third.total_npv;
    assert!((observed - predicted).abs() < 1e-2);
}

#[test]
fn market_context_typed_named_and_ad_views_share_curves() {
    let ctx = snapshot().to_market_context();
    let discount = ctx.discount("UsdSofr").expect("discount alias");
    let forward = ctx.forward("UsdSofr").expect("forward alias");
    assert!((discount.discount_factor(1.0) - forward.discount_factor(1.0)).abs() > 1e-8);
    let guard = TapeScope::guard();
    let ad = ctx.to_ad();
    let out = ad.discount("UsdSofr").unwrap().discount_factor(1.0);
    let adj = TapeScope::gradient_from(out);
    assert!(adj.adjoints.iter().any(|x| x.abs() > 0.0));
    let _ = guard.stop();
}

#[test]
fn bermudan_swaption_is_curve_dirty_and_has_par_dv01() {
    let m = snapshot();
    let trade = PortfolioTrade::BermudanSwaption {
        id: 99,
        notional: 1_000_000.0,
        fixed_rate: 0.02,
        exercise_years: vec![1, 2],
        maturity_years: 5,
        is_payer: true,
        hw_a: 0.05,
        hw_sigma: 0.01,
    };
    let mut dag = PortfolioDag::new();
    let first = dag.evaluate_portfolio(&m, std::slice::from_ref(&trade));
    assert!(first.total_par_dv01.iter().any(|x| x.abs() > 0.0));
    let mut eq = m.clone();
    eq.eq_spot += 1.0;
    let second = dag.evaluate_portfolio(&eq, std::slice::from_ref(&trade));
    assert_eq!(second.last_trade_recomputes, 0);
    let mut curve = m;
    curve.ois_quotes[0] += 1e-4;
    let third = dag.evaluate_portfolio(&curve, std::slice::from_ref(&trade));
    assert!(third.last_trade_recomputes > 0);
}
