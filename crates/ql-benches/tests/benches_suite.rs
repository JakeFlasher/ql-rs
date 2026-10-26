use ql_benches::harness::with_config_override;
use ql_benches::run_all_benchmark_suites_with_config;

#[test]
fn benchmarks_verify_parity_and_optional_sota_speedups() {
    // A short deterministic run is enough to exercise all 20 numerical
    // parity paths during normal `cargo test`; timing gates are inherently
    // host-dependent and opt in through QL_RUN_TIMING_GATES=1.
    let rows = run_all_benchmark_suites_with_config(1, 3);
    assert_eq!(rows.len(), 53);
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

#[test]
fn phase36_suite48_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("48_")));
    let row = ql_benches::benchmark_phase36_brigo_csa_ctd_fabozzi_callable_oas_lando_jlt_migration_gatheral_corridor_varswap_glasserman_lsm_rogers_dual_cartea_lob_mm_lo_mackinlay_vr_lopez_dsr_purged_cv_roncalli_spinu_erc_ilmanen_vmc_penman_reoi_mccall_bellman(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase37_suite49_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("49_")));
    let row = ql_benches::benchmark_phase37_brigo_nss_krd_tuckman_tips_floor_schonbucher_rmv_rt_rfv_gatheral_sticky_lv_gregory_simm_ccp_mva_schwartz_smith_kalman_wystup_fx_dnt_bouzoubaa_timer_meucci_bl_woodbury_gueant_ac_sinh_zou_elastic_net_nk_ree(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase38_suite50_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("50_")));
    let row = ql_benches::benchmark_phase38_hagan_west_tension_duffie_lando_pac_cmo_fama_macbeth_shanken_gatheral_asvi_bgk_discrete_barrier_bouchaud_cross_impact_karnosky_singer_merton_jump_campbell_viceira_ohlson_lid_griewank_revolve(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase39_suite51_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("51_")));
    let row = ql_benches::benchmark_phase39_brigo_xccy_mtm_fabozzi_spec_pool_io_po_lando_marshall_olkin_gatheral_vix_vvix_gregory_wwr_cva_geman_gibson_schwartz_clark_fx_fva_rubinstein_chooser_grinold_flam_tc_meucci_es_euler_james_stein_sure_campbell_cochrane_habit(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase40_suite52_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("52_")));
    let row = ql_benches::benchmark_phase40_brigo_ois_fallback_fabozzi_collared_frn_schonbucher_cds_swaption_fep_bergomi_2f_skew_gregory_sa_ccr_ead_eydeland_spark_spread_wystup_prdc_geske_compound_qian_alpha_decay_roncalli_enb_johansen_vecm_bansal_yaron_lrr(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn phase41_suite53_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("53_")));
    let row = ql_benches::benchmark_phase41_andersen_cms_replication_tuckman_special_repo_okane_cds_upfront_gatheral_ssr_tanaka_gregory_saccr_kva_eydeland_crack_spread_wystup_prdc_cap_haug_eso_dilution_qian_spinu_erc_rockafellar_cvar_htf_pls_krylov_piazzesi_fomc(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}
