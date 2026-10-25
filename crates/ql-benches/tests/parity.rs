use ql_benches::{run_all_benchmarks, suites::benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad};
use ql_benches::harness::with_config_override;

#[test]
fn phase21_suite33_and_invariants() {
    let rows = with_config_override(1, 3, || run_all_benchmarks(1));
    assert_eq!(rows.len(), 33);
    let row = with_config_override(1, 3, || benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad(1));
    assert!(row.max_abs_diff < 1e-4);
    assert!(row.parity_ok);
}
