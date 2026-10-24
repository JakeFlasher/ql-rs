# Whole-Project Adversarial Review Log (Phases 1–12, All 23 Crates at `f1e93e5`)

This document records the multi-round whole-project adversarial audit across all 23 crates (`crates/*`, Phases 1 through 12):
- **Big-Picture & Workflow Architect** (`d7266665-6096-4a67-b5c8-cc65aa5bd569`)
- **Rust Systems & AD/SIMD Auditor** (`fdb49aaf-a928-4e2c-94a2-37df1b042e24`)
- **Mathematical & Pricing Engine Auditor** (`1e332c19-68a5-4c2c-b2c0-f48b78083eee`)
- **Risk, Portfolio, Credit, XVA & SOTA Literature Auditor** (`f94752af-fc1c-4eff-92e6-1872ed9bc4d1`)

---

## 1. Whole-Project Pass 1 Verification (`e6fb004..f1e93e5`)

- **`[WP-01]` (`ql-models` `rough.rs` & `models_suite.rs`)**: **100% RESOLVED** — Decaying power-law Volterra weights, paired antithetic RNG seeds, empirical martingale control (`< 1e-10` put–call parity), and SLV particle leverage tests verified by all 4 auditors.
- **`[WP-05]` (`ql-xva` `frtb.rs`)**: **100% RESOLVED** — Separate bucket-level upward/downward curvature (`max(k_up, k_dn)`), MAR21.102 $\psi(S_u, S_v)$ indicator, MAR21.4(5)(b) capped-$\hat{S}_b$ fallback, and coherent Acerbi–Tasche fractional tail ES verified by all 4 auditors.
- **`[WP-04]` (`ql-portfolio` `execution.rs` & `allocation.rs`)**: `execution.rs` is **100% RESOLVED**; `allocation.rs` (`project_capped_simplex` and `try_mean_cvar_weights`) has a 2-item active-set KKT / normalized-step polish queued in Pass 2.
- **`[WP-02]` (`ql-exotics` `autocallable.rs`) & `[WP-03]` (`ql-rates` `g2.rs`, `lmm.rs`, `cheyette.rs`)**: `cheyette.rs` is **100% RESOLVED**; `autocallable.rs` (`Dual<64>` to eliminate $N \ge 5$ lane collision + 4-decimal `ref_spot`), `g2.rs` (exact 5-point Gauss–Hermite normal weights + $\mathbb{Q}^{T_N}$ terminal bond deflation), and `lmm.rs:230` (`cont[p] = call_threshold;`) are queued in Pass 2 alongside `[WP-06]`.

---

## 2. Whole-Project Pass 2 Remediation Scope (`CODEX_WHOLE_PROJECT_PASS2_PROMPT.md`)

| ID | Severity | Location | Remediation in Pass 2 |
| :--- | :--- | :--- | :--- |
| **`[WP-02a/b]`** | **HIGH** | `crates/ql-exotics/src/autocallable.rs:65-72, 170-225, 358-364`, `exotics_suite.rs` | Use `Dual<64>` in `price_worst_of_autocallable_aad` (eliminating the `corr_only` lane collision for $N \ge 5$ baskets while keeping 1-pass evaluation up to $N = 9$), round `ref_spot` to `1e-4`, and add a 5-asset test in `exotics_suite.rs`. |
| **`[WP-03a/b]`** | **HIGH** | `crates/ql-rates/src/g2.rs:385-496`, `crates/ql-rates/src/lmm.rs:230` | Use exact 5-point Gauss–Hermite normal nodes/weights and $\mathbb{Q}^{T_N}$ terminal bond deflation in `G2ppBermudanSwaptionEngine::price_dual_curve`; replace `cont[p].min_val(call_threshold)` with `call_threshold` in `lmm.rs:230`. |
| **`[WP-04a/b]`** | **MEDIUM** | `crates/ql-portfolio/src/allocation.rs:680-765` | Compute `shift` in `T: Real` over the free set in `project_capped_simplex` and normalize gradient steps with backtracking line search on the smooth Rockafellar–Uryasev objective in `try_mean_cvar_weights`. |
| **`[WP-06a]`** | **HIGH** | `crates/ql-math/src/solvers.rs:135-150`, `crates/ql-fixed-income/src/cat_bond.rs:279-340` | Rewind temporary `ArenaTape` operations after each scalar `brent_solve` iteration in `brent_solve_ad`, and delegate `cat_bond_z_spread` to `brent_solve_ad` so both `PoissonExponential` and `PoissonPareto` have exact IFT derivatives and zero tape bloat. |
| **`[WP-06b]`** | **MEDIUM** | `crates/ql-fixed-income/src/mbs_and_cpi.rs:108-116`, `crates/ql-inflation/src/lib.rs:85-88` | Solve coupon-bearing `breakeven_inflation_rate` via 4 Newton steps in `T: Real`; infer calendar target month from `base_date` when `target_month` is `None` in `ZeroInflationCurve::forward_cpi`. |
| **`[WP-06c]`** | **MEDIUM** | `crates/ql-curves/src/curves.rs:178`, `crates/ql-curves/src/bootstrap.rs:495-523`, `crates/ql-dml/src/cholesky_aad.rs:233`, `crates/ql-pde/src/grid2d.rs:53` | Alias `PiecewiseCubicZeroCurve` to `MonotoneCubicZeroCurve`; use `self.sensitivity.vjp(&g)` in `value_par_dv01_and_gamma`; use `cholesky_corr_cega_adjoint` in `price_correlated_basket_call_aad`; stabilize `HestonGrid2D::new_centered` variance focus at `0.04` for $v_0 \le 0.1$. |

## 4. Whole-Project Pass 2 Remediation (`f1e93e5..working tree`)

Pass 2 completed the remaining refinement set: chunked 16-lane autocallable AAD with four-decimal spot anchoring, terminal-measure G2++ Gauss--Hermite deflation and LMM callable exercise assignment, active-set Mean-CVaR projection/line search, and cross-crate AD hygiene for Brent roots, CAT-bond z-spreads, CPI breakeven coupons, seasonal CPI months, cubic zero curves, bootstrap SVD fallbacks, Cholesky CEGA diagonals, and Heston variance grids.
