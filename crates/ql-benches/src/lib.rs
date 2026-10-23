pub mod harness;
pub mod suites;
pub use harness::{BenchConfig, ComparisonRow, format_markdown_table, run_comparison};
pub use suites::run_all_benchmarks;
