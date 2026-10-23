use ql_benches::run_all_benchmarks;

#[test]
fn benchmarks_verify_parity_and_sota_speedups() {
    let rows = run_all_benchmarks(2);
    assert_eq!(rows.len(), 20);
    let tolerances = [
        ("01_calendar_popcnt", 0.0),
        ("02_leg_simd_soa", 1e-6),
        ("03_multicurve_ift_dv01", 5e-2),
        ("04_implied_vol_ift_dual2", 1e-3),
        ("05_heston_cos_dual5", 1e-3),
        ("06_gbm_mc_simd8", 1e-9),
        ("07_scripted_autocall_fuzzy_aad", 1e-2),
        ("08_portfolio_dag_incremental", 1e-9),
        ("09_credit_cds_ift_cs01", 1e-2),
        ("10_fx_vanna_volga_dual3", 1e-3),
        ("11_capriotti_giles_cholesky_corr_aad", 1e-5),
        ("12_chebyshev_dynamic_simm_mva", 1e-3),
        ("13_reiner_rubinstein_barrier_dual5", 1e-6),
        ("14_burley_owen_rqmc_variance_reduction", 1.5e-2),
        ("15_levy_cos_vg_dual5", 1e-4),
        ("16_asb_cdo_tranche_recursion_dual", 1e-6),
        ("17_rfr_sabr_caplet_strip_dual4", 1e-4),
        ("18_xccs_multicurrency_ift_dv01", 1e-4),
        ("19_american_baw_ju_ift_dual5", 1e-6),
        ("20_risk_parity_ift_allocation_aad", 1e-8),
    ];
    for row in &rows {
        let tolerance = tolerances
            .iter()
            .find(|(id, _)| *id == row.id)
            .map(|(_, tol)| *tol)
            .expect("tolerance for benchmark row");
        assert!(row.max_abs_diff <= tolerance);
    }
    let log_mean = rows
        .iter()
        .map(|row| row.speedup.max(1e-12).ln())
        .sum::<f64>()
        / rows.len() as f64;
    let geometric_mean = log_mean.exp();
    assert!(
        geometric_mean > 5.0,
        "geometric mean speedup={geometric_mean}"
    );
    for id in [
        "01_calendar_popcnt",
        "02_leg_simd_soa",
        "03_multicurve_ift_dv01",
        "08_portfolio_dag_incremental",
        "09_credit_cds_ift_cs01",
        "15_levy_cos_vg_dual5",
        "16_asb_cdo_tranche_recursion_dual",
    ] {
        let row = rows.iter().find(|row| row.id == id).expect("benchmark row");
        assert!(row.speedup > 3.5, "{id} speedup={}", row.speedup);
    }
    for row in &rows {
        let minimum = match row.id {
            "04_implied_vol_ift_dual2" | "10_fx_vanna_volga_dual3" => 2.0,
            "12_chebyshev_dynamic_simm_mva" => 6.0,
            "14_burley_owen_rqmc_variance_reduction" => 0.1,
            "15_levy_cos_vg_dual5" => 2.5,
            "16_asb_cdo_tranche_recursion_dual" => 2.5,
            _ => 2.5,
        };
        assert!(row.speedup > minimum, "{} speedup={}", row.id, row.speedup);
    }
}
