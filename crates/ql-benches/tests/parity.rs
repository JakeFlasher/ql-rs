use ql_benches::harness::with_config_override;
use ql_benches::{
    benchmark_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad, run_all_benchmarks,
    suites::{
        benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad,
        benchmark_phase24_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp,
        benchmark_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml, phase23_workloads,
    },
};

#[test]
fn phase21_suite33_and_invariants() {
    let rows = with_config_override(1, 3, || run_all_benchmarks(1));
    assert_eq!(rows.len(), 36);
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
