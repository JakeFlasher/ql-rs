# Phase 11 Multi-Round Adversarial Review Log (`2ed957f..3eb5c00`)

This document records the complete multi-round adversarial review of **Phase 11 (Production Institutional Closure)** across all four specialized audit domains:
- **Big-Picture & Workflow Architect** (`d7266665-6096-4a67-b5c8-cc65aa5bd569`)
- **Rust Systems & AD/SIMD Auditor** (`fdb49aaf-a928-4e2c-94a2-37df1b042e24`)
- **Phase 1–9 Mathematical Auditor** (`1e332c19-68a5-4c2c-b2c0-f48b78083eee`)
- **Phase 10 & SOTA Literature Auditor** (`f94752af-fc1c-4eff-92e6-1872ed9bc4d1`)

---

## 1. Pass 1 Adversarial Audit (`2ed957f..9774fcf`) — 8 Findings (`[P11-R1-01]`..`[P11-R1-08]`)

| ID | Severity | Location | Summary | Status |
| :--- | :--- | :--- | :--- | :--- |
| **`[P11-R1-01]`** | **CRITICAL** | `crates/ql-instruments/src/pde1d.rs`, `crates/ql-pde/src/adi.rs` | Deliverable 0 skipped in Round 1: `Fd1dBlackScholesEngine::price_barrier_with_rebate` subtracted the Knock-Out hit rebate for `UpAndIn \| DownAndIn`, and `Fd2dHestonEngine::solve` omitted `\|\| grid_center.is_some()`. | Addressed in R2/R3 (1 residual edge case `[P11-R3-01]` in `pde1d.rs:146-220`) |
| **`[P11-R1-02]`** | **CRITICAL** | `crates/ql-rates/src/swaption.rs`, `crates/ql-portfolio/src/dag.rs` | Receiver `bermudan_payoff` sign inversion; uncalibrated Hull–White tree without OU mean-reversion or Arrow–Debreu shifts; fake 1D PDE wrapper; `PortfolioTrade::BermudanSwaption` hash/schedule defects. | **RESOLVED** in R2/R3 |
| **`[P11-R1-03]`** | **CRITICAL** | `crates/ql-vol/src/surface.rs` | Wrong SABR VJP slice (`&d[..3]` instead of `&[d[0], d[2], d[3]]`), SVI `NaN` log-moneyness/strike confusion in `is_butterfly_arbitrage_free`, and missing `d * d` weight scaling in `jq`. | **RESOLVED** in R2/R3 |
| **`[P11-R1-04]`** | **CRITICAL** | `crates/ql-xva/src/saccr.rs`, `crates/ql-xva/src/metrics.rs` | Unprotected `0.0.sqrt()` in SA-CCR producing `NaN` AD tangents; `sd(0.0, *end_time)` on forward-starting IR trades; WWR tangent stripping via `T::from_f64(zv)` and `epe: T::zero()`. | **RESOLVED** in R2/R3 |
| **`[P11-R1-05]`** | **CRITICAL** | `crates/ql-commodities/src/swing_storage.rs` | Non-tree stubs in `price_swing_option` and `price_gas_storage` (ignored transition probabilities, `_kappa`, `_sigma`, and step discounting). | **RESOLVED** in R2/R3 |
| **`[P11-R1-06]`** | **HIGH** | `crates/ql-portfolio/src/allocation.rs` | Active-set QP clamped all bounds simultaneously without checking KKT multipliers or re-solving free variables (violating $\sum_i w_i = 1$). | **RESOLVED** in R2/R3 |
| **`[P11-R1-07]`** | **HIGH** | `crates/ql-exotics/src/variance_swap.rs` | Corridor variance swap midpoint indicator stripped barrier AD derivatives and omitted analytic log-contract tails. | **RESOLVED** in R2/R3 |
| **`[P11-R1-08]`** | **CRITICAL** | `crates/ql-benches/src/suites.rs`, `crates/ql-terminal/src/lib.rs` | Hardcoded `0.123456789` benchmark stub in `22_sabr_svi_ift_and_qp_kkt_aad`, identical `f64` closures in `21_bermudan_swaption_hw1f_aad`, unwired TUI tabs 19/20, and missing Phase 11 unit tests. | **RESOLVED** in R2/R3 |

---

## 2. Pass 2 Adversarial Audit (`9774fcf..e2695bc`) — 7 Findings (`[P11-R2-01]`..`[P11-R2-07]`)

| ID | Severity | Location | Summary | Pass 3 Verdict (`e2695bc..3eb5c00`) |
| :--- | :--- | :--- | :--- | :--- |
| **`[P11-R2-01]`** | **CRITICAL** | `crates/ql-instruments/src/pde1d.rs`, `crates/ql-pde/src/adi.rs` | Barrier survival rebates, origin breach handling, boundary tolerances, and ADI parity. | **RESOLVED** in Round 4 |
| **`[P11-R2-02]`** | **CRITICAL** | `crates/ql-vol/src/surface.rs` | `CalibratedSviSmile::vol_and_quote_vega` passed total-variance gradient $\nabla_\theta w$ instead of implied-volatility gradient $\nabla_\theta \sigma_{\text{BS}}$; `floor`/`wing` penalties duplicated $21\times$ inside `ngrid` loop. | **RESOLVED** (`surface.rs:548-563, 634-677`, `vol_suite.rs:46-61`) |
| **`[P11-R2-03]`** | **HIGH** | `crates/ql-commodities/src/swing_storage.rs` | `price_gas_storage` initialized `values` at `m = n` to `0.0` and only looped `0..n`, ignoring the final decision date $t_M$; `price_swing_option` did not sort/deduplicate `exercise_times`. | **RESOLVED** (`swing_storage.rs:168-174, 282-298`, `commodities_suite.rs:97-128`) |
| **`[P11-R2-04]`** | **HIGH** | `crates/ql-exotics/src/variance_swap.rs` | `replicate_corridor_variance_swap_fair_strike` double-compounded `put_tail`/`call_tail` by $e^{rT}$ and used a `1e-4 / 1e4` delegation bypass. | **RESOLVED** (`variance_swap.rs:181-272`, `exotics_suite.rs:144-171`) |
| **`[P11-R2-05]`** | **HIGH** | `crates/ql-portfolio/src/allocation.rs`, `crates/ql-portfolio/src/dag.rs` | Missing `try_long_only_markowitz_weights` / `long_only_markowitz_weights`; `lambda = 0.0` when `free.is_empty()`; XOR hash folding of `curve_hash ^ equity_hash` in `trade_dependency_hash`. | **RESOLVED** (`allocation.rs:65-86, 230-252`, `dag.rs:145`, `allocation_suite.rs:53-59`, `portfolio_suite.rs:105-130`) |
| **`[P11-R2-06]`** | **MEDIUM** | `crates/ql-xva/src/metrics.rs`, `crates/ql-ad/src/dual.rs` | `compute_xva_with_wwr` early-returned at `wwr_sensitivity == 0`, zeroing out $\left.\frac{\partial \text{CVA}}{\partial b}\right|_{b=0}$; `Dual2::sqrt` lacked a `g == 0.0` guard. | **RESOLVED** (`metrics.rs:233-244`, `dual.rs:605-612`, `xva_suite.rs:299-323`) |
| **`[P11-R2-07]`** | **HIGH** | `crates/ql-benches/src/suites.rs`, `crates/ql-terminal/src/lib.rs` | `benchmark_sabr_qp` had zero active box constraints and bumped only `mu` in FD while SOTA seeded `mu` + `upper`; `DeskWorkstation::step_tick` passed 1 quote to `CalibratedSabrSmile::calibrate` and did not render lines 19 & 20 in Quadrant 3. | **RESOLVED** (`suites.rs:1775-1928`, `ql-terminal/src/lib.rs:1218-1260, 1560-1565`) |

---

## 3. Pass 3 Adversarial Audit (`e2695bc..3eb5c00`) — 1 Remaining Finding (`[P11-R3-01]`)

Three of the four audit domains (`Big-Picture & Workflow Architect`, `Rust Systems & AD/SIMD Auditor`, and `Phase 10 & SOTA Literature Auditor`) reported **zero remaining findings** across `e2695bc..3eb5c00`. The `Phase 1–9 Mathematical Auditor` identified **one remaining MEDIUM defect** in `crates/ql-instruments/src/pde1d.rs`:

### `[P11-R3-01]` (MEDIUM) — **RESOLVED in Round 4** — Inconsistent `1e-12` Barrier Tolerance in Tridiagonal System Assembly & Lower-Boundary `lowb` Override for Call Options (`crates/ql-instruments/src/pde1d.rs:146-220`)

1. **Missing `eps` in Tridiagonal Row Override (`crates/ql-instruments/src/pde1d.rs:204-216` and `147, 167, 186`)**:
   - While `let eps = 1e-12 * b.val().abs().max(1.0);` was added to the terminal condition (`lines 96-98`), the post-step projection (`lines 237-239`), and the $t=0$ short-circuit (`line 295`), the implicit tridiagonal system row override at `lines 204-216` (as well as the boundary checks at `lines 147, 167, 186`) still uses un-toleranced inequalities (`sj >= barrier.val()` / `sj <= barrier.val()`):
     ```rust
     if let Some((barrier, bt, rb)) = barrier {
         for j in 1..n {
             let sj = (xmin + j as f64 * dx).exp();
             let knocked = (matches!(bt, BarrierType::UpAndOut) && sj >= barrier.val())
                 || (matches!(bt, BarrierType::DownAndOut) && sj <= barrier.val());
     ```
   - Because `xmin = barrier_level.val().ln() - j * dx` (`line 81`), floating-point rounding in `(xmin + j * dx).exp()` can land 1 ULP inside the barrier (e.g., `sj = b - 1.4e-14` for `UpAndOut` or `sj = b + 1.4e-14` for `DownAndOut`). When that occurs, `knocked` is `false` at the barrier node $j_{\text{bar}}$ during `solve_tridiagonal`, solving row $j_{\text{bar}} - 1$ as an interior PDE equation instead of the Dirichlet condition $v_{j_{\text{bar}}} = \text{rb}$ and leaking diffused value into the adjacent alive node $j_{\text{bar}} - 1$ before the post-step projection at `lines 237-242` overwrites `v[j_bar] = rb`.

2. **Asymmetric Lower Boundary (`lowb`) for `OptionType::Call` with `DownAndOut` Rebate or `terminal_rebate` (`crates/ql-instruments/src/pde1d.rs:166-203, 219`)**:
   - At the upper boundary (`j = n`), `highb` (`lines 185-201`) properly prioritizes:
     1. `UpAndOut` barrier rebate `rb` when `smax_actual >= b.val() - eps`,
     2. `terminal_rebate` discounted value `tr * (-rate * tnew).exp()`,
     3. `callb` if `kind.is_call()` else `T::zero()`.
   - However, at the lower boundary (`j = 0`), `lines 202` and `219` write:
     ```rust
     rhs[0] -= lower[0] * if kind.is_call() { T::zero() } else { putb };
     ...
     v[0] = if kind.is_call() { T::zero() } else { putb };
     ```
   - Because `if kind.is_call() { T::zero() } else { putb }` is applied *outside* `putb`, whenever `kind.is_call()` is true, the lower boundary at $S_{\min}$ (`j = 0`) is forced to `0.0` during `solve_tridiagonal` even when:
     - `barrier = Some((b, BarrierType::DownAndOut, rb))` with `rb > 0` (where $S_{\min} \le B$ is knocked out and has Dirichlet value `rb`), or
     - `terminal_rebate = Some(tr)` (e.g., the survival rebate leg of an `UpAndIn` Call at `lines 340-352`, where $S_{\min}$ is in the alive region and has asymptotic value `tr * (-rate * tnew).exp()`).
   - Defining `lowb` symmetrically with `highb` and using `rhs[0] -= lower[0] * lowb;` and `v[0] = lowb;` resolves both issues cleanly.


## 4. Pass 4 Surgical Remediation (`3eb5c00` -> Round 4)

`[P11-R3-01]` is **RESOLVED**. `Fd1dBlackScholesEngine::fd` now applies scale-aware `1e-12` barrier tolerances to all boundary, tridiagonal row, and projection comparisons; constructs a direction-aware symmetric `lowb` boundary (including knock-out rebates and terminal survival rebates); and preserves lower-boundary values for calls. Regression coverage includes non-power-of-two Down barriers, alive-region rebate monotonicity, in/out parity, and breached-at-origin short-circuiting.

All Phase 11 adversarial findings are now resolved.
