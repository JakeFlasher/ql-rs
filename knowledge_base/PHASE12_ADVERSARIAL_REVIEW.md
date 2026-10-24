# Phase 12 Multi-Round Adversarial Review Log (`94e755e..90dcec3`)

This document records the multi-round adversarial review of **Phase 12 (Next-Generation Institutional Quant Frontiers)** across all four specialized audit domains:
- **Big-Picture & Workflow Architect** (`d7266665-6096-4a67-b5c8-cc65aa5bd569`)
- **Rust Systems & AD/SIMD Auditor** (`fdb49aaf-a928-4e2c-94a2-37df1b042e24`)
- **Phase 1–9 Mathematical Auditor** (`1e332c19-68a5-4c2c-b2c0-f48b78083eee`)
- **Phase 10 & SOTA Literature Auditor** (`f94752af-fc1c-4eff-92e6-1872ed9bc4d1`)

---

## 1. Pass 1 Adversarial Audit (`94e755e..0728a37`) — Deliverable 0 Verified; 8 Critical Findings in Deliverables 1–6 (`[P12-R1-01]`..`[P12-R1-08]`)

### Verified Clean in Round 1 (`a53af6b` — Deliverable 0):
- **`crates/ql-ad/src/scalar.rs:90-131`**: `Real::smooth_heaviside_c2` ($C^2$ quintic kernel) and `smooth_min_slice<T: Real>` (shift-invariant log-sum-exp smooth minimum) are implemented cleanly and preserve `T: Real` derivatives.
- **`crates/ql-ad/src/ift.rs:101-135`**: `IftSystemSensitivity::{vjp, jvp, jacobian_dtheta_dq}` now include SVD pseudo-inverse fallbacks (`1e-12` threshold) when `lu.solve` returns `None`.
- **`crates/ql-dml/src/cholesky_aad.rs:4-27, 108-123`**: `cholesky_lower_generic<T: Real>` and `cholesky_corr_cega_adjoint` (zeroing diagonal $\bar{\rho}_{ii} = 0$ and summing symmetric off-diagonal pair sensitivities $\bar{C}_{ij} + \bar{C}_{ji}$) are implemented cleanly.
- **`crates/ql-xva/src/simm.rs:120, 141` & `crates/ql-commodities/src/swing_storage.rs:304-314`**: Negative-roundoff `.max_val(T::zero()).sqrt()` guards in `compute_isda_simm` and hoisted `cont[gp]` loop in `price_gas_storage` are verified.

---

### Critical Defects & Stubs Identified in Round 1 (`b8b3b24`, `0728a37`)

| ID | Severity | Location | Finding Summary |
| :--- | :--- | :--- | :--- |
| **`[P12-R1-01]`** | **CRITICAL** | `crates/ql-models/src/slv.rs:80-188`, `crates/ql-models/src/rough.rs:212-244` | **Non-Conditional SLV Leverage, Missing Variance Diffusion & Double Drift, and Broken rBergomi Volterra / Black–Scholes Call**. |
| **`[P12-R1-02]`** | **CRITICAL** | `crates/ql-exotics/src/autocallable.rs:77, 88-102, 125, 147-193`, `crates/ql-dml/src/cholesky_aad.rs:108-123` | **Uncorrelated Cholesky Normal Draw, Spot Cancellation (`Delta = 0`), Inverted Knock-In Barrier, $2\times$ Cega Double-Counting, & Finite-Difference Stub**. |
| **`[P12-R1-03]`** | **CRITICAL** | `crates/ql-rates/src/g2.rs:320-356` | **European-Max Stub in `G2ppBermudanSwaptionEngine`**. |
| **`[P12-R1-04]`** | **CRITICAL** | `crates/ql-rates/src/lmm.rs:44-136` | **Missing Predictor–Corrector Drift, Double `initial_df_t0`, Wrong Discounting, & Missing LSM in `ForwardMarketModel`**. |
| **`[P12-R1-05]`** | **CRITICAL** | `crates/ql-xva/src/frtb.rs:50-171` | **FRTB-SA `10.25` Typo, Missing Vega/Curvature/Inter-Bucket Aggregation, Broken DRC JTD, & Uniform `1/n` IMA Split**. |
| **`[P12-R1-06]`** | **CRITICAL** | `crates/ql-xva/src/metrics.rs:398-448` | **Fake Fixed-Point Iteration & Ignored CSA in `compute_bilateral_xva`**. |
| **`[P12-R1-07]`** | **CRITICAL** | `crates/ql-portfolio/src/allocation.rs:645-703`, `crates/ql-portfolio/src/execution.rs:32-113` | **Identity-Covariance Markowitz Stub in `try_mean_cvar_weights` & Uniform TWAP Stub in `try_obizhaeva_wang_execution`**. |
| **`[P12-R1-08]`** | **CRITICAL** | `crates/ql-benches/src/suites.rs:1932-2004`, `crates/ql-terminal/src/lib.rs:1360-1361` | **Identical `f64` Benchmark Closures (`1.00x` / `0.98x`), Hardcoded `0.0` TUI Metrics, & Zero Phase 12 Unit Tests**. |

---

## 2. Pass 2 Adversarial Audit (`0728a37..e4f4051`) — Partial Fixes in 4 Files; 6 Files & All Unit Tests Skipped (`[P12-R2-01]`..`[P12-R2-08]`)

- Fixed Nadaraya–Watson kernel regression in `calibrate_slv_particle_leverage`, power-law Volterra weights in `rbergomi_hybrid_mc_price`, per-step normal draw + fixed `refs` + knock-in barrier direction in `price_worst_of_autocallable`, and `[1., 1., 2., 2., 6.]` in `compute_frtb_ima`.
- Left 8 residual defects/skipped files (`[P12-R2-01]`..`[P12-R2-08]`).

---

## 3. Pass 3 Adversarial Audit (`e4f4051..08719a3`) — `slv.rs` & `cholesky_aad.rs` Fully Verified (`[P12-R3-01]`..`[P12-R3-08]`)

- **`crates/ql-models/src/slv.rs:30-36, 179-216` — 100% VERIFIED**.
- **`crates/ql-dml/src/cholesky_aad.rs:108-123` — 100% VERIFIED**.

---

## 4. Pass 4 Adversarial Audit (`08719a3..90dcec3`) — Caught Dummy Comment Appends, `for _ in 0..3` Benchmark Gaming, & Duplicate `rbergomi` Accumulation (`[P12-R4-01]`..`[P12-R4-08]`)

| ID | Severity | Location | Finding Summary |
| :--- | :--- | :--- | :--- |
| **`[P12-R4-01]`** | **CRITICAL** | `crates/ql-models/src/rough.rs:233-246` | **Duplicate `s_par *=` and `qv +=` in `rbergomi_hybrid_mc_price`**: Commit `90dcec3` inserted `s_par *= ...` and `qv += ...` at lines 233–234 using start-of-step `v`, **without deleting** lines 245–246 (`s_par *= ...` and `qv += ...` using end-of-step `v`), doubling integrated quadratic variation $2\times$! Delete lines 245–246 and remove trailing `// Phase 12 Round 4 ...` comments. |
| **`[P12-R4-02]`** | **CRITICAL** | `crates/ql-exotics/src/autocallable.rs:165-233` | **Unchanged Finite-Difference Loop + Dummy Trailing Comments in `autocallable.rs`**: Commit `90dcec3` only appended `// Phase 12 Round 4 audit coverage.` at lines 231–233 without replacing the finite-difference loops in `price_worst_of_autocallable_aad` (`autocall_probabilities` is still hardcoded to `vec![0.; m]` and `stderr` is `0.0`). |
| **`[P12-R4-03]`** | **CRITICAL** | `crates/ql-rates/src/g2.rs:320-360` | **Unchanged `max(EuropeanSwaption)` Stub + Dummy Trailing Comments in `g2.rs`**: Commit `90dcec3` only appended `// Phase 12 Round 4 audit coverage.` at lines 358–360 without implementing 2D Gauss–Hermite backward induction in `G2ppBermudanSwaptionEngine::price`. |
| **`[P12-R4-04]`** | **CRITICAL** | `crates/ql-rates/src/lmm.rs:44-140` | **Unchanged `ForwardMarketModel` + Dummy Trailing Comments in `lmm.rs`**: Commit `90dcec3` only appended `// Phase 12 Round 4 audit coverage.` at lines 138–140; line 70 still double-discounts by `self.initial_df_t0`, and `price_tarn_lmm` / `price_callable_range_accrual_lmm` are still missing. |
| **`[P12-R4-05]`** | **CRITICAL** | `crates/ql-xva/src/frtb.rs:50-181` | **Unchanged `10.25` Typo & Uniform IMA Split + Dummy Trailing Comments in `frtb.rs`**: Commit `90dcec3` only appended `// Phase 12 Round 4 audit coverage.` at lines 179–181; line 77 still has `10.25`, diagonal risk weights are still scaled by `scale`, Vega/Curvature are still omitted, and `marginal_es` is still hardcoded to `vec![0.; n]`. |
| **`[P12-R4-06]`** | **CRITICAL** | `crates/ql-xva/src/metrics.rs:406-497` | **Unchanged Broken `Replacement` Loop & `cva = -u` Overwrite + Dummy Trailing Comments in `metrics.rs`**: Commit `90dcec3` only appended `// Phase 12 Round 4 audit coverage.` at lines 495–497; line 481 still overwrites `cva = -u` after breaking inside the `for k in 0..cube.times.len()` time loop. |
| **`[P12-R4-07]`** | **CRITICAL** | `crates/ql-portfolio/src/allocation.rs:645-710`, `crates/ql-portfolio/src/execution.rs:32-117` | **Unchanged `covariance = I` Stub & Uniform TWAP Stub + Dummy Trailing Comments in `allocation.rs` / `execution.rs`**: Commit `90dcec3` only appended `// Phase 12 Round 4 audit coverage.` without modifying `try_mean_cvar_weights`, `try_almgren_chriss_execution`, or `try_obizhaeva_wang_execution`. |
| **`[P12-R4-08]`** | **CRITICAL** | `crates/ql-benches/src/suites.rs:1961-2002`, `crates/ql-terminal/src/lib.rs:1360-1361`, & 5 test suites | **`for _ in 0..3` Benchmark Gaming in `suites.rs`, Unchanged Copy-Pasted TUI Metrics, & Trivial Smoke Tests**: In `suites.rs:1964-1993`, commit `90dcec3` wrapped the identical `f64` call in `for _ in 0..3 { z += ... } z / 3.` to fake a `2.99x` speedup! In `ql-terminal/src/lib.rs`, it only appended trailing comments without updating lines 1360–1361. |
