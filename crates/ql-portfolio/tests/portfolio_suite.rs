use ql_models::HestonParams;
use ql_portfolio::{MarketSnapshot, PortfolioDag, PortfolioTrade};
use ql_scripting::{BytecodeProgram, ScriptedProduct};
use ql_time::Date;

fn snapshot() -> MarketSnapshot {
    MarketSnapshot {
        ref_date: Date::from_ymd(2024, 1, 1).unwrap(),
        ois_quotes: [0.02, 0.021, 0.022],
        fwd_quotes: [0.025, 0.026, 0.027],
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
    assert_eq!(second.curve_recomputes, 1);
    assert_eq!(second.trade_recomputes, 5);

    s.eq_spot = 102.0;
    let third = dag.evaluate_portfolio(&s, &ts);
    assert_eq!(third.curve_recomputes, 1);
    assert_eq!(third.trade_recomputes, 8);

    s.fwd_quotes[2] += 1e-4;
    let fourth = dag.evaluate_portfolio(&s, &ts);
    assert_eq!(fourth.curve_recomputes, 2);
    assert_eq!(fourth.trade_recomputes, 13);
    let predicted = third.total_par_dv01[5] * 1e-4;
    let observed = fourth.total_npv - third.total_npv;
    assert!((observed - predicted).abs() < 1e-2);
}
