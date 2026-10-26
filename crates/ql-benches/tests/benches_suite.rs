use ql_benches::harness::with_config_override;
use ql_benches::run_all_benchmark_suites_with_config;

#[test]
fn benchmarks_verify_parity_and_optional_sota_speedups() {
    // A short deterministic run is enough to exercise all 20 numerical
    // parity paths during normal `cargo test`; timing gates are inherently
    // host-dependent and opt in through QL_RUN_TIMING_GATES=1.
    let rows = run_all_benchmark_suites_with_config(1, 3);
    assert_eq!(rows.len(), 47);
    for row in &rows {
        assert!(!row.id.is_empty());
        assert!(!row.description.is_empty());
        assert!(!row.category.is_empty());
        assert!(row.baseline_median_ns.is_finite() && row.baseline_median_ns > 0.0);
        assert!(row.sota_median_ns.is_finite() && row.sota_median_ns > 0.0);
        assert!(row.speedup.is_finite() && row.speedup > 0.0);
        assert!(row.parity_ok, "{} parity flag is false", row.id);
        assert!(row.max_abs_diff.is_finite());
        assert!(
            row.max_abs_diff <= row.parity_tol,
            "{} parity diff={} tol={}",
            row.id,
            row.max_abs_diff,
            row.parity_tol
        );
    }

    if std::env::var("QL_RUN_TIMING_GATES").is_ok() {
        let log_mean = rows
            .iter()
            .map(|row| row.speedup.max(1e-12).ln())
            .sum::<f64>()
            / rows.len() as f64;
        assert!(log_mean.exp() >= 1.50, "geomean speedup={}", log_mean.exp());
        for row in &rows {
            assert!(row.speedup >= 1.05, "{} speedup={}", row.id, row.speedup);
        }
    }
}

#[test]
fn phase22_suite34_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("34_"))
        .expect("suite 34");
    assert!(row.parity_ok);
}
#[test]
fn phase23_suite35_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("35_"))
        .expect("suite 35");
    assert!(row.parity_ok);
}

#[test]
fn phase24_suite36_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("36_"))
        .expect("suite 36");
    assert!(row.parity_ok);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase25_suite37_workload_is_registered() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_phase25_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin(1)
    });
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase26_suite38_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("38_"))
        .expect("suite 38");
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn phase27_suite39_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("39_"))
        .expect("suite 39");
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase28_suite40_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("40_"))
        .expect("suite 40");
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase29_suite41_workload_is_registered() {
    let rows = ql_benches::run_all_benchmarks(1);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("41_"))
        .expect("suite 41");
    assert!(row.max_abs_diff < 1e-4);
}

#[test]
fn phase30_suite42_workload_is_registered() {
    let rows = ql_benches::run_all_benchmarks(1);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("42_"))
        .expect("suite 42");
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
}

#[test]
fn phase31_suite43_workload_is_registered() {
    let rows = ql_benches::run_all_benchmarks(1);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("43_"))
        .expect("suite 43");
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
}

#[test]
fn phase32_suite44_workload_is_registered() {
    let rows = ql_benches::run_all_benchmarks(1);
    assert!(rows.iter().any(|r| r.id.starts_with("44_")));
    let row = ql_benches::benchmark_phase32_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase33_suite45_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("45_")));
    let row = ql_benches::benchmark_phase33_rebonato_cms_lm_all_engines(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase34_suite46_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("46_")));
    let row = ql_benches::benchmark_phase34_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase35_suite47_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("47_")));
    let row = ql_benches::benchmark_phase35_hagan_west_ctd_afv_kth_sacva_gs_kirk_fx_vv_napoleon_tc_meucci_krr_eva(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}
