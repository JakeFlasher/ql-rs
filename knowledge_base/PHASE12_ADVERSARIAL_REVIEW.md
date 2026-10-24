# Whole-Project Adversarial Review Log (Phases 1–12, All 23 Crates at `e6fb004`)

Following convergence of the Phase 12 Round 7 diff (`8261494..e6fb004`), the adversarial audit scope was expanded per user directive to a **comprehensive, whole-project audit across all 23 crates (`crates/*`) spanning Phases 1 through 12**:
- **Big-Picture & Workflow Architect** (`d7266665-6096-4a67-b5c8-cc65aa5bd569`): Audited all 23 crates, 24 benchmark suites (`ql-benches`), 22 terminal workstation modes/tabs (`ql-terminal`), integration suites, and documentation (`README.md`, `knowledge_base/*`).
- **Rust Systems & AD/SIMD Auditor** (`fdb49aaf-a928-4e2c-94a2-37df1b042e24`): Audited all 23 crates for AD (`Dual<N>`, `Dual2<N>`, `Tape`, `IftSystemSensitivity`) propagation, SIMD (`f64x4`), arena tape hygiene, and numerical guards.
- **Mathematical & Pricing Engine Auditor** (`1e332c19-68a5-4c2c-b2c0-f48b78083eee`): Audited all 12 pricing, model, calibration, curve, volatility, Monte Carlo, PDE, fixed-income, instrument, math, time, and scripting crates (`ql-models`, `ql-exotics`, `ql-rates`, `ql-pde`, `ql-mc`, `ql-vol`, `ql-curves`, `ql-fixed-income`, `ql-instruments`, `ql-math`, `ql-time`, `ql-scripting`).
- **Risk, Portfolio, Credit, XVA & SOTA Literature Auditor** (`f94752af-fc1c-4eff-92e6-1872ed9bc4d1`): Audited all 45 source and test files across `ql-xva`, `ql-portfolio`, `ql-commodities`, `ql-credit`, `ql-fx`, `ql-inflation`, `ql-levy`, and `ql-dml`.

---

## 1. Whole-Project Crate Audit Matrix (`e6fb004`)

| Crate | Whole-Project Status | Key Findings Across Phases 1–12 |
| :--- | :--- | :--- |
| `ql-ad` | **VERIFIED CLEAN** | `Dual<N>`, `Dual2<N>`, `ArenaTape`, `ift_lift_1d`, `IftSystemSensitivity` (with SVD fallback), `smooth_heaviside_c2`, `smooth_min_slice`, and SIMD `f64x4` verified clean. |
| `ql-time` | **VERIFIED CLEAN** | 64-bit civil `Date`, `BitsetCalendar` (POPCNT `business_days_between`), 9 day-count conventions, and `ScheduleBuilder` verified clean. |
| `ql-mc` | **VERIFIED CLEAN** | Philox-4x32-10 RNG, SIMD Box–Muller, Joe–Kuo Sobol + Burley–Owen scrambling, Brownian bridge, Andersen QE-M Heston, and two-pass Chebyshev LSM AMC verified clean. |
| `ql-scripting` | **VERIFIED CLEAN** | AST + bytecode VM with fuzzy pathwise AAD verified clean. |
| `ql-levy` | **VERIFIED CLEAN** | VG, CGMY, Merton, Kou, and Fang–Oosterlee COS engines verified clean (minor test addition for `cos_levy_price_strip`). |
| `ql-models` | **`[WP-01]` (`CRITICAL`)** | `rbergomi_hybrid_mc_price` (`rough.rs:216, 220`) applies an extra `.powf(1. / (h - 0.5))` that inverts the power-law Volterra decay and seeds antithetic paths with `p` instead of `p / 2`; `models_suite.rs` lacks SLV and rBergomi convergence tests. |
| `ql-exotics` | **`[WP-02]` (`CRITICAL`)** | `price_worst_of_autocallable` (`autocallable.rs:76-86`) only rounds integer spots (`(z - z.round()).abs() < 1e-2`), zeroing out FD spot deltas at non-integer spots (e.g. `100.25`), and `price_worst_of_autocallable_aad` (`lines 240-254`) hardcodes `stderr`, `autocall_probabilities`, `expected_life`, and `corr_cega_stderr` while truncating $N \ge 5$ baskets at `Dual<16>`. |
| `ql-rates` | **`[WP-03]` (`CRITICAL`)** | (1) `G2ppBermudanSwaptionEngine::price` (`g2.rs:327-371`) uses `0.75 * ||v||_2 + 0.30 * best` instead of 2D $(x, y)$ grid backward induction and ignores `_ibor_forward`; (2) `ForwardMarketModel` (`lmm.rs:72-132`) omits cross-rate correlation, Hunter–Jäckel–Joshi predictor–corrector drift, `displacements`, `backward_looking_rfr`, `smoothing_eps`, and LSM continuation regression; (3) `Cheyette1F` (`cheyette.rs:44-93`) omits discounting in `simulate_short_rate_call`. |
| `ql-portfolio` | **`[WP-04]` (`HIGH`)** | `try_mean_cvar_weights` (`allocation.rs:676-732`) solves a heuristic second-moment Markowitz surrogate instead of the exact Rockafellar–Uryasev (2000) smooth objective; `try_obizhaeva_wang_execution` (`execution.rs:102, 168`) has an off-by-one post-trade inventory index (`inv[k]` vs `inv[k + 1]`). |
| `ql-xva` | **`[WP-05]` (`MEDIUM`)** | `bucket_curv` (`frtb.rs:95-124`) mixes upward and downward curvature shocks per risk factor instead of bucket-level $\max(K_b^+, K_b^-)$, `sbm_scenario_charge` (`lines 165-189`) omits MAR21.102 $\psi(S_b, S_c)$ and MAR21.4(5)(b) fallback, and `compute_frtb_ima` (`line 279`) uses integer `.ceil()` tail count instead of fractional Acerbi–Tasche tail weights. |
| `ql-fixed-income`, `ql-math`, `ql-commodities`, `ql-inflation`, `ql-curves`, `ql-credit`, `ql-dml` | **`[WP-06]` (`HIGH`/`MEDIUM`)** | Cross-crate AD/numerical hygiene: `cat_bond_z_spread` (`cat_bond.rs:280-339`) and `brent_solve_ad` (`solvers.rs:142`) record scalar Brent iterations onto active `ArenaTape`; `breakeven_inflation_rate` (`mbs_and_cpi.rs:108-116`) ignores real coupons; `price_swing_option` (`swing_storage.rs:19-74`) ignores `min_exercises`; `bootstrap_with_seasonality` (`ql-inflation/src/lib.rs:172-253`) passes `target_month: None`; `PiecewiseCubicZeroCurve` (`curves.rs:178`) aliases linear `ZeroCurve`; `cdo.rs:229` and `cholesky_aad.rs:233` minor AD polish. |
| `ql-benches` & `ql-terminal` | **`[WP-07]` (`HIGH`)** | Suites 06, 15, 20 baseline parity & Suites 21–24 subsystem coverage (`suites.rs`); `ql-terminal` `TabId` / `chart_mode` header sync and full Phase 12 metric wiring (`lib.rs`, `main.rs`). |

---

## 2. Two-Pass Whole-Project Remediation Roadmap

- **Pass 1 (`CODEX_WHOLE_PROJECT_PASS1_PROMPT.md`)**: Remediates all `CRITICAL` and `HIGH` mathematical/pricing/risk findings in `ql-models`, `ql-exotics`, `ql-rates`, `ql-portfolio`, and `ql-xva` (`[WP-01]` through `[WP-05]`) along with substantive integration tests in `models_suite.rs`, `exotics_suite.rs`, and `rates_suite.rs`.
- **Pass 2 (`CODEX_WHOLE_PROJECT_PASS2_PROMPT.md`)**: Remediates the cross-crate AD/numerical hygiene items (`[WP-06]`) and benchmark/terminal/documentation alignment (`[WP-07]`).

## 3. Whole-Project Pass 1 Remediation (`e6fb004..working tree`)

The five Pass 1 findings were remediated: rough Bergomi now uses decaying Volterra weights, true antithetic seeds, and a discrete martingale correction; the autocallable engine anchors decimal reference spots, reports path/bin diagnostics, and handles larger correlation baskets; G2++ Bermudan pricing uses a dual-curve two-factor grid, the LMM uses correlated predictor/corrector terminal paths with LSM callable continuation, and Cheyette paths are discounted; portfolio Mean-CVaR and transient execution use their specified smooth/KKT objectives; and FRTB curvature direction/fallback and fractional Acerbi--Tasche tails are implemented. Substantive model, exotic, and rates invariants pass with workspace clippy/tests.
