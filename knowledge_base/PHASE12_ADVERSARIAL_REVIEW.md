# Whole-Project Adversarial Review Log (Phases 1–12, All 23 Crates at `f1340fe`)

This document records the multi-round whole-project adversarial audit across all 23 crates (`crates/*`, Phases 1 through 12):
- **Big-Picture & Workflow Architect** (`d7266665-6096-4a67-b5c8-cc65aa5bd569`)
- **Rust Systems & AD/SIMD Auditor** (`fdb49aaf-a928-4e2c-94a2-37df1b042e24`)
- **Mathematical & Pricing Engine Auditor** (`1e332c19-68a5-4c2c-b2c0-f48b78083eee`)
- **Risk, Portfolio, Credit, XVA & SOTA Literature Auditor** (`f94752af-fc1c-4eff-92e6-1872ed9bc4d1`)

---

## 1. Whole-Project Pass 1 (`e6fb004..f1e93e5`) & Pass 2 (`f1e93e5..f1340fe`) Verification

All **21 core quantitative library crates (`ql-ad` through `ql-xva`)** have reached **`0 CRITICAL / 0 HIGH / 0 MEDIUM` convergence** across all 4 specialist audit domains:
- **`Mathematical & Pricing Engine Auditor` (`1e332c19`)**: **CONVERGED (`0 CRITICAL / 0 HIGH / 0 MEDIUM`)** across `ql-models`, `ql-rates`, `ql-exotics`, `ql-vol`, `ql-curves`, `ql-fixed-income`, `ql-mc`, `ql-pde`, `ql-instruments`, `ql-math`, `ql-time`, and `ql-scripting`.
- **`Risk, Portfolio, Credit, XVA & SOTA Literature Auditor` (`f94752af`)**: **CONVERGED (`0 CRITICAL / 0 HIGH / 0 MEDIUM`)** across `ql-portfolio`, `ql-xva`, `ql-credit`, `ql-commodities`, `ql-fx`, `ql-inflation`, `ql-levy`, and `ql-dml`.
- **`Rust Systems & AD/SIMD Auditor` (`fdb49aaf`)**: **CONVERGED (`0 CRITICAL / 0 HIGH / 0 MEDIUM`)** across all 23 crates (`ql-ad` through `ql-terminal`).
- **`Big-Picture & Workflow Architect` (`d7266665`)**: Confirmed **`[WP-01]` through `[WP-06]` are 100% RESOLVED & CONVERGED**, leaving only the final **`[WP-07]` (`ql-benches`, `ql-terminal`, `README.md`)** benchmark/workstation/documentation synchronization in **Whole-Project Pass 3**.

---

## 2. Whole-Project Pass 3 Final Polish Scope (`CODEX_WHOLE_PROJECT_PASS3_PROMPT.md`)

| ID | Severity | Location | Remediation in Pass 3 |
| :--- | :--- | :--- | :--- |
| **`[WP-07a]`** | **MEDIUM** | `crates/ql-benches/src/suites.rs:497-513, 1123-1174, 1573-1654, 1932-2109` | **1. Suite 06 (`scalar_gbm_price`)**: Accumulate log-returns across steps and call `libm::exp` once at maturity (matching `GbmSoaSimulator::price_european`).<br>**2. Suite 15 (`benchmark_levy_cos`)**: Use `cos_levy_price_strip` in the baseline closure (matching the SOTA closure).<br>**3. Suite 20 (`benchmark_risk_parity`)**: Use 2-point central FD (`h = 1e-6`, `tol = 1e-5`) instead of 4-solve Richardson extrapolation in the baseline closure, and remove the duplicate `cov_ad` pass in the SOTA closure.<br>**4. Suites 23 & 24 (`benchmark_slv_autocall` & `benchmark_frtb_cvar_execution`)**: Include `rbergomi_hybrid_mc_price` (`Dual<3>`) alongside `price_worst_of_autocallable_aad` in Suite 23, and include `compute_frtb_sa` (`Dual<4>`) alongside `almgren_chriss_execution` (`Dual<4>`) in Suite 24 so their workloads match their suite IDs/labels. |
| **`[WP-07b]`** | **MEDIUM** | `crates/ql-terminal/src/lib.rs:1268-1296, 1617, 1641`, `crates/ql-terminal/src/main.rs:35-41`, `crates/ql-terminal/tests/terminal_suite.rs:78` | **1. `render_dashboard`**: Use `self.chart_mode.label()` at lines 1617 & 1641 so `set_chart_mode` updates both the Braille canvas and the header/footer labels.<br>**2. `step_tick`**: Include `price_worst_of_autocallable` in `rough_metric` (`slv_rough_worst_of_metric`) and `compute_frtb_sa` in `execution_metric` (`frtb_cvar_execution_metric`).<br>**3. `main.rs`**: Include `ScenarioKind::StagflationPlus150Bp` in the `--snapshot` loop. |
| **`[WP-07c]`** | **LOW** | `README.md:38, 53-75`, `crates/ql-benches/BENCHMARK_RESULTS.md`, `knowledge_base/08_BENCHMARK_RESULTS.md` | Update `README.md:38` (`20-suite` -> `24-suite`) and sync the 24-row benchmark table across `BENCHMARK_RESULTS.md`, `knowledge_base/08_BENCHMARK_RESULTS.md`, and `README.md`. |

## 5. Whole-Project Pass 3 Polish (`f1340fe..working tree`)

Final WP-07 polish aligned scalar GBM, Lévy COS, risk-parity, SLV/rBergomi/autocallable, and FRTB/execution benchmark workloads; synchronized live terminal chart/scenario/metric wiring; and refreshed the 24-suite benchmark tables in the reports and README.

## 6. Whole-Project Pass 3 Polish (`f1340fe..working tree`)

WP-07 benchmark workload parity, terminal chart/scenario/metric synchronization, and README/report consistency are complete. Release regeneration now records all 24 suites, including aligned SLV/rBergomi/autocallable and FRTB/execution workloads.
