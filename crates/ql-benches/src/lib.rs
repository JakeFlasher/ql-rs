pub mod harness;
pub mod suites;
/// Number of registered benchmark suites in the current release.
pub const SUITE_COUNT: usize = 33;
pub use harness::{BenchConfig, ComparisonRow, format_markdown_table, run_comparison};
pub use suites::{
    run_all_benchmark_suites_with_Config, run_all_benchmark_suites_with_config, run_all_benchmarks,
    run_all_suites,
};
