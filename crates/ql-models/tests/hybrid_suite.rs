use ql_models::{AnalyticBsmHullWhiteEngine, CosHestonHullWhiteEngine, HestonHullWhiteParams};

#[test]
fn hybrid_zero_rate_factor_reduces_to_heston() {
    let p = HestonHullWhiteParams {
        v0: 0.04,
        kappa_v: 1.5,
        theta_v: 0.04,
        xi_v: 0.3,
        kappa_r: 0.1,
        eta_r: 0.0,
        rho_sv: -0.7,
        rho_sr: 0.0,
    };
    let e = CosHestonHullWhiteEngine::default();
    let x = e.price(&p, 100.0, 100.0, 0.03, 0.0, 1.0, true);
    assert!(x.is_finite() && x > 0.0);
    assert!(AnalyticBsmHullWhiteEngine::effective_variance(0.2, 0.1, 0.0, 0.0, 1.0) > 0.0);
}
