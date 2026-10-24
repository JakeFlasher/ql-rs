# Phase 13 Multi-Agent Adversarial Audit & Remediation Log

## Overview
- **Audit Scope**: Phase 13 Institutional Expansion (`9221898..1cb863a`) across all 6 Phase 13 pillars (`ql-curves`, `ql-vol`, `ql-models`, `ql-exotics`, `ql-mc`, `ql-rates`, `ql-inflation`, `ql-credit`, `ql-xva`, `ql-portfolio`), Benchmark Suite 25 (`ql-benches`), Terminal Tab 23 (`ql-terminal`), and whole-project changes (`ql-dml`, `ql-fixed-income`, `README.md`, `knowledge_base/*`).
- **Specialist Audit Panel**:
  1. **Big-Picture & Workflow Architect (`d7266665`)**: Public API completeness, cross-crate architecture, Benchmark Suite 25, Terminal Tab 23, and test coverage.
  2. **Rust Systems & AD/SIMD Auditor (`fdb49aaf`)**: Generic `T: Real` AD/IFT correctness (`Dual`, `Dual2`, `AadVar`), KKT bordered Hessians, and Benchmark Suite 25 parity diagnostics.
  3. **Mathematical & Pricing Engine Auditor (`1e332c19`)**: Hagan–West Monotone Convex & step-jump OIS curves, global eSSVI & VIX strip/Heston–Bates quadrature, Parisian/Parasian 2D Markov lattice & Andersen–Broadie primal–dual bounds, and Cheyette/HW1F/LPI inflation engines.
  4. **Risk, Portfolio, Credit, XVA & SOTA Literature Auditor (`f94752af`)**: Heterogeneous $k$-th-to-default baskets, MPoR cashflow-aware XVA & Contingent CDS, Ledoit–Wolf shrinkage, SOCP Robust Ellipsoidal Markowitz, and primary literature verification (`papers/60_*.pdf` through `papers/66_*.pdf`).

---

## Round 1 Audit Verdict (`9221898..1cb863a`): **NOT CONVERGED**
- **Summary of Findings**: **5 CRITICAL / 6 HIGH / 1 MEDIUM / 1 LOW** (13 distinct items).
- **Verified Clean in Round 1**:
  - Whole-project remediation edits in `crates/ql-dml/src/{cholesky_aad.rs,pca.rs,regression.rs}`, `crates/ql-fixed-income/src/{cat_bond.rs,convertible_bond.rs,mbs_and_cpi.rs,nss.rs}`, and `crates/ql-mc/src/paths.rs`.
  - Model-free CBOE VIX replication strip `replicate_cboe_vix_index` (`crates/ql-models/src/vix.rs:14-58`).
  - 2D Markov Parisian/Parasian lattice `price_parisian_option_lattice` (`crates/ql-exotics/src/parisian.rs:146-215`).

### Round 1 Findings Table

| ID | Severity | Crate / Location | Finding Summary |
|---|---|---|---|
| **P13-R1-C1** | **CRITICAL** | `crates/ql-mc/src/andersen_broadie.rs:16-129` | `AmcEngine::andersen_broadie_primal_dual_bounds` omits the Doob–Meyer martingale subtraction $\pi_k$, ignores `n_inner`, and prematurely `break`s out of the outer path loop upon primal exercise $\hat{\tau}$, computing a truncated lookback payoff instead of the Rogers (2002) / Andersen–Broadie (2004) dual upper bound. |
| **P13-R1-C2** | **CRITICAL** | `crates/ql-models/src/vix.rs:60-110` | `heston_bates_vix_future_and_option_price` uses Jensen's upper bound $\sqrt{\mathbb{E}[\text{VIX}_T^2]}$ (zero vol-of-vol convexity adjustment and zero $\sigma_v$ sensitivity on `vix_future`), passes `vvol = 0.5 * sigma_v * expiry.sqrt()` to `black_scholes_price` (which multiplies by `expiry.sqrt()` a second time), and linearizes the Bates jump compensator. |
| **P13-R1-C3** | **CRITICAL** | `crates/ql-xva/src/metrics.rs:513-532` | `compute_mpor_xva` is a 1-line pass-through stub that ignores `_period_cashflows`, `s.mpor_years`, and `s.unpaid_cashflow_weight`, delegating directly to zero-MPoR `compute_xva`. |
| **P13-R1-C4** | **CRITICAL** | `crates/ql-rates/src/cheyette.rs:201-241` | `price_callable_range_accrual_hw1f` ignores `_model: &HullWhite1F<T>`, evaluates deterministic forwards, and sets `callable = non.min_val(max_coupon)`, yielding `embedded_call_value == 0` identically for all positive curves and all $(a, \sigma)$. |
| **P13-R1-C5** | **CRITICAL** | `crates/ql-*/tests/*` | Zero (`0`) unit or integration tests were added across all 10 domain crates (`ql-curves`, `ql-vol`, `ql-models`, `ql-exotics`, `ql-mc`, `ql-rates`, `ql-inflation`, `ql-credit`, `ql-xva`, `ql-portfolio`) for Phase 13 Pillars 1–6. |
| **P13-R1-H1** | **HIGH** | `crates/ql-curves/src/curves.rs:222-370, 433-445` & `bootstrap.rs:893-967` | `MonotoneConvexZeroCurve<T>` evaluates a single unconstrained quadratic instead of the Hagan–West (2006) monotonicity clamping and 4-region piecewise-quadratic spline (allowing negative forwards even when `enforce_positivity == true`); `MeetingDateJumpOisCurve::instantaneous_forward` applies `meeting_jumps[0]` on $t \in [0, m_0)$ due to `0usize.saturating_sub(1) == 0`; `CalibratedMeetingJumpOisCurve` uses FD instead of `Dual<1>` for `jt` and full matrix inversion instead of `vjp`. |
| **P13-R1-H2** | **HIGH** | `crates/ql-vol/src/surface.rs:857-911, 949-1010` | `CalibratedEssviSurface::vol_and_quote_vega` writes into `g[..6]` instead of `g[(i-1)*3..(i+1)*3]` (corrupting quote vegas for any surface with $\ge 3$ slices when $i \ge 2$), uses `s.expiry` instead of `expiry` at `i == 0`, and `calibrate` omits Hendriks–Martini / Mingone butterfly and calendar no-arbitrage constraints. |
| **P13-R1-H3** | **HIGH** | `crates/ql-exotics/src/parisian.rs:280-284` & `crates/ql-benches/src/suites.rs:2261-2500` | `price_parisian_option_mc_ad` passes `s.window_duration.val()` as `half_width` to `clock.smooth_heaviside_c2`, widening the gate to $[0, 2D]$ and creating a `3.440e-1` AD vs FD parity error in Benchmark Suite 25 (masked by `parity_tol: 5e-1` and leftover `QL_DEBUG_PHASE13` `eprintln!` blocks). |
| **P13-R1-H4** | **HIGH** | `crates/ql-rates/src/cheyette.rs:136-193` | `calibrate_cheyette_coterminal_ad` omits Ornstein–Uhlenbeck variance decay `exp(-2 * mu * dt)` when solving for `sigma2` (`line 186` vs `line 189`), causing `model_normal_vols` to miss `market_normal_vols` for $i \ge 1$ when $\mu \ne 0$, and strips `T: Real` derivatives in `swap_loading`. |
| **P13-R1-H5** | **HIGH** | `crates/ql-credit/src/cdo.rs:407-437` & `crates/ql-xva/src/metrics.rs:533-575` | `price_kth_to_default_basket` replaces entity-specific recovery rates $R_i$ with unweighted average LGD instead of leave-one-out cavity conditioning $P^{(k-1,-i)}(t | Z)$ (violating conservation of basket protection under heterogeneous $R_i$) and omits half-period default accrual; `price_contingent_cds` omits survival probability and `running_spread` in `prem`. |
| **P13-R1-H6** | **HIGH** | `crates/ql-portfolio/src/allocation.rs:905-1100` | `ledoit_wolf_shrinkage` omits the sample-size factor $1/n$ in $\hat{\delta}^* = \text{clamp}(\frac{\hat{\pi}-\hat{\rho}}{n\hat{\gamma}}, 0, 1)$, omits $\hat{\rho}$ for `ConstantCorrelation`, and omits `DiagonalVariance`; `try_robust_markowitz_weights` omits the rank-1 SOCP Hessian term $-\frac{\kappa}{\sigma_{\text{est}}^3}(\Omega_\mu w)(\Omega_\mu w)^\top$ from the KKT IFT system. |
| **P13-R1-M1** | **MEDIUM** | `crates/ql-inflation/src/lib.rs:450-465` | `lpi_inflation_swap_npv` ignores `vol` and `model` when `is_compounded_zero_coupon == true`, giving zero inflation vega on compounded zero-coupon LPI swaps. |
| **P13-R1-L1** | **LOW** | `crates/ql-terminal/src/lib.rs:21-27` & `knowledge_base/*` | `ql_terminal::prelude` omits Phase 12–13 `ql_models` modules (`vix`, `rough`, `slv`), and `knowledge_base/{05,07}*.md` have minor stale counts (`24-suite` / `22-tab` instead of `25-suite` / `23-tab`). |
