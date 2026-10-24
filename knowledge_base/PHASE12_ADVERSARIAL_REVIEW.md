# Phase 12 Multi-Round Adversarial Review Log (`94e755e..8261494`)

This document records the multi-round adversarial review of **Phase 12 (Next-Generation Institutional Quant Frontiers)** across all four specialized audit domains:
- **Big-Picture & Workflow Architect** (`d7266665-6096-4a67-b5c8-cc65aa5bd569`) — **CONVERGED (`0 CRITICAL / 0 HIGH / 0 MEDIUM`) in `8261494`**
- **Rust Systems & AD/SIMD Auditor** (`fdb49aaf-a928-4e2c-94a2-37df1b042e24`) — **CONVERGED (`0 CRITICAL / 0 HIGH / 0 MEDIUM`) in `8261494`**
- **Phase 1–9 Mathematical Auditor** (`1e332c19-68a5-4c2c-b2c0-f48b78083eee`) — **CONVERGED (`0 CRITICAL / 0 HIGH / 0 MEDIUM`) in `8261494`** (1 minor `LOW` guard in `g2.rs:357`)
- **Phase 10 & SOTA Literature Auditor** (`f94752af-fc1c-4eff-92e6-1872ed9bc4d1`) — **5 `CRITICAL` / 1 `HIGH` remaining in `ql-xva` and `ql-portfolio` (`[P12-R6-01]`..`[P12-R6-06]`)**

---

## 1. Modules Verified 100% Clean & Converged Through Round 6 (`94e755e..8261494`)

1. **Deliverable 0 (`ql-ad`, `ql-dml`, `ql-xva::simm`, `ql-commodities::swing_storage`)**:
   - `crates/ql-ad/src/scalar.rs:90-131` (`smooth_heaviside_c2`, `smooth_min_slice`) and `crates/ql-ad/src/ift.rs:95-136` (`IftSystemSensitivity` SVD fallback) — **VERIFIED CLEAN**.
   - `crates/ql-dml/src/cholesky_aad.rs:4-27, 108-123` (`cholesky_lower_generic<T: Real>`, `cholesky_corr_cega_adjoint`) — **VERIFIED CLEAN**.
   - `crates/ql-xva/src/simm.rs:120, 141` & `crates/ql-commodities/src/swing_storage.rs:304-314` — **VERIFIED CLEAN**.
2. **Deliverable 1 (`ql-models::slv` & `ql-models::rough`)**:
   - `crates/ql-models/src/slv.rs:1-216` (`CalibratedSlvModel::leverage_at`, `calibrate_slv_particle_leverage` Nadaraya–Watson Gaussian kernel regression + `mixing_fraction` + empirical martingale projection, `price_slv_european_mc`) — **VERIFIED CLEAN**.
   - `crates/ql-models/src/rough.rs:178-253` (`rbergomi_hybrid_mc_price` with adapted start-of-step variance `v`, power-law Volterra weights, correlated $2\times 2$ head covariance, and Romano–Touzi conditional Black–Scholes) — **VERIFIED CLEAN**.
3. **Deliverable 2 (`ql-exotics::autocallable`)**:
   - `crates/ql-exotics/src/autocallable.rs:65-256` (`price_worst_of_autocallable<T: Real>` and single-pass `Dual<16>` `price_worst_of_autocallable_aad`) — **VERIFIED CLEAN**.
4. **Deliverable 3 (`ql-rates::g2` & `ql-rates::lmm`)**:
   - `crates/ql-rates/src/g2.rs:327-367` (`G2ppBermudanSwaptionEngine::price` and `price_dual_curve`; minor `LOW` guard at line 357 for single-exercise Bermudan parity) — **VERIFIED CLEAN**.
   - `crates/ql-rates/src/lmm.rs:44-151` (`ForwardMarketModel::price_caplet_mc`, `price_tarn_mc`, `price_callable_range_accrual_mc`, `price_tarn_lmm`, `price_callable_range_accrual_lmm`) — **VERIFIED CLEAN**.
5. **Deliverable 6 (`ql-benches` & `ql-terminal`)**:
   - `crates/ql-benches/src/suites.rs:1932-2109` & `crates/ql-benches/BENCHMARK_RESULTS.md:28-29`:
     - `23_slv_rbergomi_and_autocallable_cega_aad`: **`5.84x` speedup** (`724,122 ns` vs `123,917 ns`, `2.817e-6` Parity Diff) — **VERIFIED CLEAN**.
     - `24_frtb_sbm_ima_and_cvar_execution_aad`: **`4.61x` speedup** (`2,850.40 ns` vs `618.40 ns`, `3.181e-7` Parity Diff) — **VERIFIED CLEAN**.
   - `crates/ql-terminal/src/lib.rs:1268-1296, 1391-1392`: `DeskWorkstation::step_tick` live Phase 12 telemetry — **VERIFIED CLEAN**.

---

## 2. Pass 6 Adversarial Audit (`2dc3082..8261494`) — Remaining Findings in `ql-xva` & `ql-portfolio` (`[P12-R6-01]`..`[P12-R6-06]`)

| ID | Severity | Location | Finding & Required Remediation in Round 7 |
| :--- | :--- | :--- | :--- |
| **`[P12-R6-01]`** | **CRITICAL** | `crates/ql-xva/src/frtb.rs:50-106` | **Complete BCBS MAR21/22 `compute_frtb_sa` (Vega, Curvature, Low-Correlation Floor, Inter-Bucket Aggregation, and Unweighted DRC HBR)**: (1) Include `vega_sensitivities` / `vega_risk_weights` and `curvature_up_pnl` / `curvature_down_pnl` / `curvature_shock_weights` in bucket charges; (2) Use MAR21.7 low-correlation formula $\rho^{(\text{Low})} = \max(2\rho - 1, 0.75\rho)$; (3) Aggregate across buckets within each `FrtbRiskClass` using `s.inter_bucket_corr`; (4) Compute DRC `HBR = gl / (gl + gs)` from unweighted gross JTDs before applying `(wl - hbr * ws).max_val(T::zero())`. |
| **`[P12-R6-02]`** | **CRITICAL** | `crates/ql-xva/src/frtb.rs:144-178` | **Compute Exact Subgradient / Euler Allocation in `compute_frtb_ima`**: Replace hardcoded `marginal_es: vec![0.; n]` and uniform `euler_component_es: vec![liq / n; n]` with exact tail-scenario subgradients $g_{i,k} = \partial \text{ES}_i / \partial w_k = -\frac{1}{n_i}\sum_{s \in \text{tail}_i} r_{s,k}$ chained through $\text{ES}_{\text{liq}} = \sqrt{\sum_{i=0}^4 m_i \text{ES}_i^2}$ ($m = [1, 1, 2, 2, 6]$), setting `euler_component_es[k] = w[k].val() * marginal_es[k]` so $\sum_k \text{euler\_component\_es}[k] = \text{liq.val()}$. |
| **`[P12-R6-03]`** | **CRITICAL** | `crates/ql-xva/src/metrics.rs:413-477` | **Backward Time Induction & Per-Step Fixed-Point Iteration in `compute_bilateral_xva`**: Loop backward `for k in (0..cube.times.len()).rev()`, floor exposures per netting set `ns`, and when `Replacement` is active, run an inner fixed-point loop `for _ in 0..8` at step `k` solving $u_k = (u_{\text{cont}} + \text{step\_xva}(u_k)) / \text{df}_k$ until $|u_k^{\text{new}} - u_k| < 10^{-12}$. |
| **`[P12-R6-04]`** | **CRITICAL** | `crates/ql-portfolio/src/allocation.rs:645-712` | **Full Tail-Scenario Co-Movement Matrix & Smooth Rockafellar–Uryasev Hinge in `try_mean_cvar_weights`**: Replace the diagonal semivariance matrix with iteratively reweighted smooth Rockafellar–Uryasev tail weights $\pi_m = \phi_\varepsilon'(L_m - \zeta) / ((1 - \alpha) M)$ across the full $n \times n$ scenario matrix (including off-diagonal $i \ne j$), and evaluate smooth CVaR $\zeta + \frac{1}{(1-\alpha)M}\sum_m \phi_\varepsilon(L_m - \zeta)$. |
| **`[P12-R6-05]`** | **CRITICAL** | `crates/ql-portfolio/src/execution.rs:40-41, 81-128` | **Use $\tilde{\eta} = \eta - \frac{1}{2}\gamma \Delta t$ in `try_almgren_chriss_execution` & Solve Exact $m \times m$ Resilience KKT System in `try_obizhaeva_wang_execution`**: Use `eta_tilde = (s.temporary_impact - s.permanent_impact * dt * 0.5).max_val(T::from_f64(1e-12))` in `try_almgren_chriss_execution`, and in `try_obizhaeva_wang_execution` solve $H y = \mathbf{1}, \xi^* = X_0 y / (\mathbf{1}^\top y)$ where $H_{ij} = \gamma e^{-\rho|t_i - t_j|} + 2\eta_{\text{temp}}\delta_{ij} + 2\lambda \sigma^2 (t_{\min(i,j)} - t_0)$, with exact quadratic expected cost and inventory variance. |
| **`[P12-R6-06]`** | **HIGH** | `crates/ql-xva/tests/xva_suite.rs:325-329`, `crates/ql-portfolio/tests/allocation_suite.rs:61-73` | **Replace Smoke Tests with Substantive Numerical Invariant Tests**: Verify FRTB SA/IMA Euler decomposition, Bilateral XVA `Replacement` vs `RiskFree` fixed-point convergence, smooth Mean-CVaR tail sensitivity, and Obizhaeva–Wang KKT execution invariants. |

## 3. Pass 7 Remediation (`8261494..working tree`) — `[P12-R6-01]` through `[P12-R6-06]` resolved

- Replaced the FRTB-SA/IMA stubs with scenario-scaled Delta/Vega/Curvature SBM, low-correlation floors, inter-bucket aggregation, unweighted DRC HBR, and exact tail Euler allocations.
- Reworked bilateral XVA into backward per-netting-set exposure induction with replacement-closeout fixed-point iteration.
- Replaced diagonal mean-CVaR covariance with smooth Rockafellar--Uryasev tail reweighting and solved Obizhaeva--Wang execution through its full resilience KKT system; Almgren--Chriss uses effective temporary impact.
- Added the single-exercise G2++ guard and substantive FRTB/IMA/XVA and mean-CVaR/execution invariant tests. Workspace clippy and tests pass.
