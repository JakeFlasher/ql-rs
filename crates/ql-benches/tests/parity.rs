use ql_benches::harness::with_config_override;
use ql_benches::{
    benchmark_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad,
    benchmark_phase25_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin,
    run_all_benchmarks,
    suites::{
        benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad,
        benchmark_phase24_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp,
        benchmark_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml, phase23_workloads,
    },
};

#[test]
fn phase21_suite33_and_invariants() {
    let rows = with_config_override(1, 3, || run_all_benchmarks(1));
    assert_eq!(rows.len(), 39);
    let row = with_config_override(1, 3, || {
        benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad(1)
    });
    assert!(row.max_abs_diff < 1e-4);
    assert!(row.parity_ok);
}

#[test]
fn phase23_suite35_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml(1)
    });
    assert_eq!(
        row.id,
        "35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml"
    );
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
    let values = phase23_workloads([0.2_f64, -0.1, 0.3, -0.2]);
    assert_eq!(values.len(), 12);
    assert!(values.iter().all(|v| v.is_finite()));
}

#[test]
fn phase22_suite34_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad(1)
    });
    assert_eq!(row.id, "34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad");
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
}

#[test]
fn phase24_suite36_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_phase24_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp(1)
    });
    assert!(row.parity_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn phase25_suite37_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_phase25_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin(1)
    });
    assert!(row.parity_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn phase26_suite38_and_invariants() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_phase26_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc(1)
    });
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn phase27_suite39_and_invariants() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_phase27_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt(1)
    });
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}
