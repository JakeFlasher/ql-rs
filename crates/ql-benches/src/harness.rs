use std::hint::black_box;
use std::time::Instant;

#[derive(Copy, Clone, Debug)]
pub struct BenchConfig {
    pub warmup_iters: usize,
    pub samples: usize,
    pub iters_per_sample: usize,
}

#[derive(Clone, Debug)]
pub struct ComparisonRow {
    pub id: &'static str,
    pub description: &'static str,
    pub baseline_median_ns: f64,
    pub baseline_iqr_ns: f64,
    pub sota_median_ns: f64,
    pub sota_iqr_ns: f64,
    pub speedup: f64,
    pub max_abs_diff: f64,
}

pub fn measure_ns<F: FnMut() -> f64>(cfg: BenchConfig, mut f: F) -> (f64, f64, f64) {
    for _ in 0..cfg.warmup_iters {
        black_box(f());
    }
    if cfg.samples == 0 || cfg.iters_per_sample == 0 {
        return (0.0, 0.0, 0.0);
    }
    let mut samples = Vec::with_capacity(cfg.samples);
    let mut last = 0.0;
    for _ in 0..cfg.samples {
        let start = Instant::now();
        for _ in 0..cfg.iters_per_sample {
            last = black_box(f());
        }
        let elapsed = start.elapsed();
        samples.push(elapsed.as_nanos() as f64 / cfg.iters_per_sample as f64);
    }
    samples.sort_by(f64::total_cmp);
    let n = samples.len();
    let median = if n % 2 == 0 {
        0.5 * (samples[n / 2 - 1] + samples[n / 2])
    } else {
        samples[n / 2]
    };
    let iqr = samples[3 * n / 4] - samples[n / 4];
    (median, iqr, last)
}

pub fn run_comparison<B, S>(
    id: &'static str,
    description: &'static str,
    cfg: BenchConfig,
    tol: f64,
    baseline: B,
    sota: S,
) -> ComparisonRow
where
    B: FnMut() -> f64,
    S: FnMut() -> f64,
{
    let (base_median, base_iqr, base_value) = measure_ns(cfg, baseline);
    let (sota_median, sota_iqr, sota_value) = measure_ns(cfg, sota);
    let diff = (base_value - sota_value).abs();
    assert!(
        diff <= tol,
        "benchmark {id} parity failed: baseline={base_value:.17e}, sota={sota_value:.17e}, diff={diff:.3e}, tol={tol:.3e}"
    );
    ComparisonRow {
        id,
        description,
        baseline_median_ns: base_median,
        baseline_iqr_ns: base_iqr,
        sota_median_ns: sota_median,
        sota_iqr_ns: sota_iqr,
        speedup: base_median / sota_median.max(1e-9),
        max_abs_diff: diff,
    }
}

pub fn format_markdown_table(rows: &[ComparisonRow]) -> String {
    let mut out = String::from(
        "| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |\n|---|---|---:|---:|---:|---:|\n",
    );
    for row in rows {
        out.push_str(&format!(
            "| {} | {} | {:.2} (IQR {:.2}) | {:.2} (IQR {:.2}) | {:.2}x | {:.3e} |\n",
            row.id,
            row.description,
            row.baseline_median_ns,
            row.baseline_iqr_ns,
            row.sota_median_ns,
            row.sota_iqr_ns,
            row.speedup,
            row.max_abs_diff,
        ));
    }
    out
}
