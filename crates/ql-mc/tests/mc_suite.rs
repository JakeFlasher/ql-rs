use ql_instruments::{Fd1dBlackScholesEngine, OptionType};
use ql_mc::*;
use ql_models::*;
use ql_vol::black76_greeks;

#[test]
fn qmc_sobol_bridge_and_simd_gbm_match_black_scholes() {
    let sim = GbmSoaSimulator::new(100.0, 0.2, 0.03, 0.01, 1.0, 8_192, 16, 42);
    let exact = sim.analytical(100.0, true);
    let qmc = sim.price_european_qmc(100.0, true);
    let (mc, se) = sim.price_european(100.0, true);

    assert!((qmc - exact).abs() < 1e-2, "qmc={qmc} exact={exact}");
    assert!(
        (mc - exact).abs() < 3.5 * se + 2e-2,
        "mc={mc} exact={exact} se={se}"
    );
}

#[test]
fn heston_qe_matches_cos_in_feller_violating_regime() {
    let p = HestonParams {
        v0: 0.04,
        kappa: 1.0,
        theta: 0.04,
        sigma_v: 0.5,
        rho: -0.7,
    };
    assert!(!p.satisfies_feller());
    let sim = HestonQeSimulator::new(p, 100.0, 0.02, 0.0, 1.0, 16_384, 32, 7);
    let (mc_price, se) = sim.price_european(100.0, true);
    let cos_price = CosHestonEngine::default().price(&p, 100.0, 100.0, 0.02, 0.0, 1.0, true);
    assert!(
        (mc_price - cos_price).abs() < 3.5 * se + 2e-2,
        "mc={mc_price} cos={cos_price} se={se}"
    );
}

#[test]
fn tape_checkpoint_pathwise_greeks_and_constant_memory() {
    let spot: f64 = 100.0;
    let vol: f64 = 0.2;
    let rate: f64 = 0.03;
    let div: f64 = 0.0;
    let strike: f64 = 100.0;
    let expiry: f64 = 1.0;

    let (p, delta, vega, _rho, len) =
        mc_gbm_pathwise_tape(spot, vol, rate, div, strike, expiry, 16_384, 42);
    let fwd = spot * ((rate - div) * expiry).exp();
    let df = (-rate * expiry).exp();
    let (bs_p, bs_d_fwd, _, bs_vega, _, _, _, _) =
        black76_greeks(fwd, strike, vol, expiry, df, true);
    let bs_delta = bs_d_fwd * ((rate - div) * expiry).exp();

    assert!((p - bs_p).abs() < 1e-1, "p={p} bs_p={bs_p}");
    assert!(
        (delta - bs_delta).abs() < 2e-2,
        "delta={delta} bs_delta={bs_delta}"
    );
    assert!(
        (vega - bs_vega).abs() < 6e-1,
        "vega={vega} bs_vega={bs_vega}"
    );
    assert!(len <= 20, "peak tape len={len}");
}

#[test]
fn two_pass_amc_matches_1d_pde_american_put() {
    let amc = AmcEngine {
        n_calibration_paths: 8_192,
        n_valuation_paths: 16_384,
        n_steps: 32,
        basis_degree: 3,
        seed: 123,
    };
    let (amc_put, se) = amc.price_american(105.0, 100.0, 0.20, 0.05, 1.0, false);
    let pde_put = Fd1dBlackScholesEngine::new(200, 200).price_american_ad(
        100.0,
        105.0,
        1.0,
        0.20,
        0.05,
        OptionType::Put,
    );
    let eur_put = ql_vol::black_scholes_price(100.0, 105.0, 0.20, 1.0, 0.05, 0.0, false);

    assert!(amc_put > eur_put + 0.3, "amc={amc_put} eur={eur_put}");
    assert!(
        (amc_put - pde_put).abs() < 3.5 * se + 8e-2,
        "amc={amc_put} pde={pde_put} se={se}"
    );
}

#[test]
fn burley_owen_points_are_stratified() {
    let g = BurleyOwenSobol::new(3, 17);
    let mut bins = vec![vec![0u32; 1024]; 3];
    let mut point = vec![0.0; 3];
    for i in 0..1024u32 {
        g.point_at(i, &mut point);
        for d in 0..3 {
            let b = (point[d] * 1024.0).floor() as usize;
            bins[d][b.min(1023)] += 1;
        }
    }
    assert!(bins.iter().all(|v| v.iter().all(|&x| x == 1)));
}

#[test]
fn rqmc_asian_returns_finite_stderr() {
    let (price, se) = price_asian_arithmetic_rqmc(
        100.0,
        100.0,
        0.03,
        0.0,
        0.2,
        &(1..=16).map(|i| i as f64 / 16.0).collect::<Vec<_>>(),
        1.0,
        true,
        128,
        4,
        19,
    );
    assert!(price > 0.0 && se > 0.0 && se < price);
    let mut samples = Vec::with_capacity(2048);
    let rng = PhiloxRng::new(19);
    for p in 0..2048u64 {
        let mut log_s = 100.0f64.ln();
        let mut prev = 0.0;
        let mut sum = 0.0;
        for (j, tm) in (1..=16).map(|i| i as f64 / 16.0).enumerate() {
            let dt = tm - prev;
            log_s += (0.03 - 0.5 * 0.2 * 0.2) * dt + 0.2 * dt.sqrt() * rng.normal4(p, j as u32)[0];
            sum += log_s.exp();
            prev = tm;
        }
        samples.push(((sum / 16.0) - 100.0).max(0.0) * (-0.03f64).exp());
    }
    let mean = samples.iter().sum::<f64>() / samples.len() as f64;
    let mc_se = (samples.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
        / (samples.len() as f64 * (samples.len() - 1) as f64))
        .sqrt();
    assert!(mc_se > se * 5.0, "rqmc_se={se} mc_se={mc_se}");
    assert!(price > 0.0, "price={price}");
}
