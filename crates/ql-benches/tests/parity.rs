use ql_benches::harness::with_config_override;
use ql_benches::{
    benchmark_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad, run_all_benchmarks,
    suites::benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad,
};

#[test]
fn phase21_suite33_and_invariants() {
    let rows = with_config_override(1, 3, || run_all_benchmarks(1));
    assert_eq!(rows.len(), 34);
    let row = with_config_override(1, 3, || {
        benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad(1)
    });
    assert!(row.max_abs_diff < 1e-4);
    assert!(row.parity_ok);
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
