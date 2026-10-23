# CODEX MASTER PROMPT — Phase 10 Comprehensive Adversarial Remediation (Architecture, Rust/AD Core, Phases 1–9 & Phase 10 SOTA Parity)

You are operating as a Principal Quantitative Architect and Rust Systems Lead on `/home/jakeshea/rust_quantlib` (`ql-rs`).

We have completed a multi-pass, workspace-wide adversarial review across all 23 crates (`ql-ad`, `ql-time`, `ql-math`, `ql-curves`, `ql-instruments`, `ql-vol`, `ql-models`, `ql-mc`, `ql-pde`, `ql-rates`, `ql-fixed-income`, `ql-credit`, `ql-inflation`, `ql-fx`, `ql-commodities`, `ql-levy`, `ql-exotics`, `ql-dml`, `ql-scripting`, `ql-portfolio`, `ql-xva`, `ql-benches`, `ql-terminal`), cross-checked against the 45 papers in `/home/jakeshea/rust_quantlib/papers/` and the 12 reference repositories in `/home/jakeshea/rust_quantlib/repos/` (`QuantLib`, `Strata`, `ORE`, `rateslib`, `quantsupport`, etc.).

The detailed findings are recorded in:
- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_ARCHITECTURE.md`
- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_RUST_ARCH.md`
- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_MATH.md`
- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_SECOND_PASS.md`

---

## STRICT ANTI-SHORTCUT & ENGINEERING RULES (ZERO TOLERANCE)

1. **NO ALIASING DISTINCT ENGINES**: Never implement a named pricing engine (`BjerksundStensland1993Engine`, `BjerksundStensland2002Engine`, `JuQuadraticEngine`, `CatBondSpec::price_mc`, `CosHestonHullWhiteEngine`, `DengLiZhouSpreadBasketEngine`) by calling another engine or averaging two prices.
2. **NO SILENT `.val()` DERIVATIVE STRIPPING IN GENERIC `<T: Real>` FUNCTIONS**: Any function generic over `<T: Real>` must propagate `Dual<N>`, `Dual2<N>`, and `TapeVar<f64>` derivatives with respect to **all** `T` arguments (either via direct `<T: Real>` arithmetic or exact 1st + 2nd order Implicit Function Theorem lifting). Never convert `T` inputs to `f64` via `.val()`, solve in `f64`, and wrap constants back in `T::from_f64(...)`.
3. **NO SILENT `.unwrap_or(...)` SOLVER FALLBACKS**: Never hide a failed root solve via `.unwrap_or(guess)`, `.unwrap_or(0.0)`, or `.unwrap_or(spot)` and feed the unconverged value into `ift_lift_1d`. Use analytical bracket expansion / safeguarded Newton-bisection and return a typed `Result` when inputs are out of domain.
4. **NO UNUSED STRUCT FIELDS OR `let _ = is_call`**: Every field of a specification struct and every option orientation flag (`is_call`, `option_type`) must enter the mathematical valuation formula faithfully.
5. **Preserve Existing Public APIs Where Tested**: Keep existing public function/method names working so existing test suites continue to compile, while adding `try_*` / `Result`-returning constructors/solvers and fixing all broken formulas.

---

## WORKSTREAM I — Big-Picture Domain Architecture, Workflow Framework, `ql-portfolio` & `ql-xva`

### 1.1 Promote `MarketContext<T>` & Canonical Valuation Pipeline (`ql-curves`, `ql-instruments`, `ql-fixed-income`, `ql-portfolio`)
- **File**: `crates/ql-curves/src/market_context.rs`
  - Expand `MarketContext<T>` (and alias `pub type MarketEnvironment<T> = MarketContext<T>;`) so it is no longer dead code:
    - Store `reference_date: Date`, `discount_curves: HashMap<String, Arc<dyn YieldTermStructure<T>>>`, `forward_curves: HashMap<String, Arc<dyn YieldTermStructure<T>>>`, `spots: HashMap<String, T>`, and `vols: HashMap<String, T>`.
    - Add convenience builders and lookup methods (`with_discount`, `with_forward`, `with_spot`, `with_vol`, `discount_curve`, `forward_curve`, `spot`, `vol`).
- **File**: `crates/ql-curves/src/curves.rs:10-15` & `crates/ql-fixed-income/src/nss.rs:3-36`
  - In `YieldTermStructure::forward_rate` (`curves.rs:10-15`), compute `t1 = dc.year_fraction(self.reference_date(), d1, None, None)` and `t2 = dc.year_fraction(self.reference_date(), d2, None, None)` instead of hardcoding `(d - self.reference_date()) as f64 / 365.0`.
  - In `crates/ql-fixed-income/src/nss.rs`, add `pub reference_date: Date` (defaulting to `Date::from_ymd(2025, 1, 2).unwrap()` in existing constructors or via `with_reference_date`) and implement `YieldTermStructure<T>` for `NelsonSiegelSvensson<T>`.
- **Unify Duplicate Domain Types Across Crates**:
  - Re-export `ql_instruments::BarrierType` in `crates/ql-exotics/src/barrier.rs` (`pub use ql_instruments::BarrierType;`) instead of defining a second colliding `enum BarrierType`.
  - In `crates/ql-portfolio/src/dag.rs:293-328`, allow `PortfolioTrade::IrsSwap` to use `ResolvedVanillaSwap` / fractional accrual schedules (and ensure `PortfolioDag` can evaluate a `ResolvedVanillaSwap` against `MarketContext<T>`).

### 1.2 Fix `PortfolioDag` Cache Eviction, Telemetry Counters, Float Hashing & Pillar Labels (`crates/ql-portfolio/src/dag.rs`, `crates/ql-terminal/src/lib.rs`)
- **File**: `crates/ql-portfolio/src/dag.rs`
  1. **Canonical float hashing (`lines 10-17, 49-58, 94`)**:
     Define a single helper:
     ```rust
     #[inline]
     fn canonical_f64_bits(x: f64) -> u64 {
         if x == 0.0 {
             0.0_f64.to_bits()
         } else if x.is_nan() {
             0x7ff8_0000_0000_0000_u64
         } else {
             x.to_bits()
         }
     }
     ```
     Use `canonical_f64_bits` in `hash_f64_slice`, `hash_program`, and `trade_market_hash`.
  2. **Evict removed trade IDs (`lines 221-240`)**:
     In `evaluate_portfolio`, after validating unique trade IDs in `seen`, call `self.cache.retain(|id, _| seen.contains(id));` so removed/matured trades do not leak memory in long-running services. Also add `pub fn clear_cache(&mut self)` and `pub fn cached_trade_count(&self) -> usize`.
  3. **Fix telemetry counter semantics (`lines 206-265`)**:
     Reset `self.trade_recomputes = dirty.len();` per `evaluate_portfolio` call (matching `self.curve_recomputes = 0/1`) and add `pub total_trade_recomputes: usize` and `pub total_curve_recomputes: usize` for cumulative lifetime telemetry.
  4. **Fix duplicate ID panic (`lines 223-227`)**:
     Add `pub fn try_evaluate_portfolio(&mut self, trades: &[PortfolioTrade], market: &MarketSnapshot) -> Result<PortfolioRiskReport, String>` returning `Err` on duplicate trade IDs or curve calibration failure.
  5. **Fix TUI pillar labels (`crates/ql-terminal/src/lib.rs:894-907`)**:
     Label the 6 entries of `par_rate_dv01` accurately as `["OIS 1Y", "OIS 2Y", "OIS 5Y", "FWD 1Y", "FWD 2Y", "FWD 5Y"]` instead of mislabeling indices `3..6` as `["10Y", "20Y", "30Y"]`.

### 1.3 Fix `ql-xva` CSA Exposure Formula, Netting-Set Boundaries, Time-Integrated EPE & ISDA SIMM v2.6 (`crates/ql-xva/src/metrics.rs`, `simm.rs`, `chebyshev.rs`)
- **File**: `crates/ql-xva/src/metrics.rs`
  1. **Fix inverted `collateralized_exposure` (`lines 11-19`)**:
     Under a CSA with counterparty threshold $H_{\text{cpty}} \ge 0$ and minimum transfer amount $\text{MTA} \ge 0$, when uncollateralized exposure $E = \max(V, 0)$ exceeds $H_{\text{cpty}} + \text{MTA}$, variation margin $C = E - H_{\text{cpty}}$ is posted, leaving **remaining uncollateralized exposure** $E - C = H_{\text{cpty}}$ (or with partial/continuous AD tangent propagation, `if excess.val() > mta { e.min_val(T::from_f64(threshold_cpty)) } else { e }`) — NOT `excess`!
     Verify that for `threshold_cpty = 0.0, mta = 0.0` (zero-threshold CSA), `collateralized_exposure(mtm, 0.0, 0.0)` returns `0.0`, and for `threshold_cpty = 500.0, mta = 50.0`, `collateral.epe < raw.epe` and `collateral.cva < raw.cva`.
  2. **Respect legal netting-set boundaries (`lines 34-99`)**:
     In `expected_positive_exposure`, `potential_future_exposure`, and `compute_xva`, apply `.max_val(T::zero())` (and `collateralized_exposure` / negative exposure floor `.min_val(T::zero())`) **per netting set `ns` first**, and then sum across `ns in 0..cube.n_netting_sets`:
     ```rust
     let mut pos = T::zero();
     let mut neg = T::zero();
     for ns in 0..cube.n_netting_sets {
         let v = cube.get(p, t, ns);
         let e_pos = match csa {
             Some(c) => collateralized_exposure(v, c.threshold_cpty, c.mta),
             None => v.max_val(T::zero()),
         };
         let e_neg = match csa {
             Some(c) => -collateralized_exposure(-v, c.threshold_self, c.mta),
             None => v.min_val(T::zero()),
         };
         pos += e_pos;
         neg += e_neg;
     }
     ```
  3. **Time-integrated EPE (`line 107`)**:
     Compute `epe` using trapezoidal time integration over `cube.times` (starting from $t_0 = 0, EE(0) = 0$ or across `cube.times` weighted by $\Delta t_k / T_{\max}$) when `cube.times.last() > 0`.
- **File**: `crates/ql-xva/src/simm.rs:52-63` & `crates/ql-xva/src/chebyshev.rs:306-319`
  1. In `compute_bucket_margin` (`simm.rs:52-63`), include the ISDA SIMM v2.6 concentration correlation factor on off-diagonal terms ($i \ne j$):
     $$f_{ij} = \frac{\min(CR_i, CR_j)}{\max(CR_i, CR_j)}$$
     Wait: check `crates/ql-xva/tests/xva_suite.rs:154-191` (`isda_simm_dual_sensitivities_match_fd`), where `concentration_threshold = 1.0` and sensitivities are `[1.0, -2.0]`, `[0.5, -1.5]`; because that test compares AD vs. FD of `compute_isda_simm` itself, including $f_{ij} = \min(CR_i, CR_j) / \max(CR_i, CR_j)$ (computed in `T: Real` via `cr[i].min_val(cr[j]) / cr[i].max_val(cr[j])`) preserves exact AD==FD parity while making `compute_isda_simm` 100% ISDA SIMM v2.6 compliant!
  2. Add `SimmRiskClass` (`InterestRate`, `CreditQualifying`, `CreditNonQualifying`, `Equity`, `Commodity`, `Fx`) and `pub fn compute_isda_simm_6_risk_classes<T: Real>(class_margins: &[T; 6]) -> T` using the official ISDA SIMM v2.6 $6 \times 6$ cross-risk-class correlation matrix $\psi_{r,s}$.
  3. In `simm_from_iter` (`crates/ql-xva/src/chebyshev.rs:306-319`), net sensitivities across trades by pillar/factor before squaring if multiple sensitivities share the same factor, or provide a clear netted sensitivity aggregation.

### 1.4 Fix `cornish_fisher_var` Skewness Sign (`crates/ql-portfolio/src/risk_measures.rs:103-117`)
- **File**: `crates/ql-portfolio/src/risk_measures.rs:110-117`
  - Because `mean_pnl`, `std_pnl`, `skew`, `excess_kurtosis` are moments of **PnL** (where negative PnL is loss $L = -\text{PnL}$, so $S_L = -S_{\text{PnL}} = -\text{skew}$), negative PnL skewness (`skew < 0`) must **increase** VaR:
    ```rust
    let z = T::from_f64(ql_math::inv_norm_cdf(confidence));
    let z2 = z * z;
    let z3 = z2 * z;
    let adjusted = z - (z2 - 1.0) * skew / 6.0 + (z3 - z * 3.0) * excess_kurtosis / 24.0
        - (z3 * 2.0 - z * 5.0) * skew * skew / 36.0;
    mean_pnl.neg() + std_pnl * adjusted
    ```

### 1.5 Decouple `ql-terminal` & Fast Deterministic `cargo test` (`crates/ql-terminal`, `crates/ql-benches`)
- **File**: `crates/ql-terminal/src/lib.rs:3-28, 1411-1438`
  - Remove `pub use ql_benches::*;` from `ql_terminal::prelude`.
  - Fix `right_bottom` in `crates/ql-terminal/src/lib.rs:1411-1438` so the Phase 10 telemetry lines for Tabs 17 and 18 are visible within the 20 rendered rows (`right_bottom[0..20]`) instead of being pushed at indices 21–22 and truncated!
- **File**: `crates/ql-benches/tests/benches_suite.rs`, `crates/ql-benches/src/main.rs`, `crates/ql-benches/benches/quant_benchmarks.rs`
  - In `crates/ql-benches/tests/benches_suite.rs`, check `std::env::var("QL_RUN_TIMING_GATES").is_ok()`: by default during `cargo test`, run `run_all_benchmark_suites_with_Config(1, 3)` (1 warmup, 3 samples — takes $< 2\text{s}$!) and assert all 20 suite names, categories, finite positive timings, and **strict numerical parity (`r.parity_ok`, `r.max_abs_diff <= r.parity_tol`)**, while only enforcing hardware-dependent wall-clock `speedup` thresholds when `QL_RUN_TIMING_GATES=1` is set!
  - In `crates/ql-benches/src/main.rs` and `crates/ql-benches/benches/quant_benchmarks.rs`, only overwrite `BENCHMARK_RESULTS.md` and `knowledge_base/08_BENCHMARK_RESULTS.md` when `--write-report` or `QL_WRITE_BENCH_REPORT=1` is set (note: after fixing Suite 19 & Suite 20 in Workstream IV below, regenerate and commit the updated `crates/ql-benches/BENCHMARK_RESULTS.md` and `knowledge_base/08_BENCHMARK_RESULTS.md` once!).

---

## WORKSTREAM II — Rust Systems, AD/IFT Core (`ql-ad`), Linear Algebra, SIMD & Date Safety

### 2.1 Fix `0.0 * Inf = NaN` Inactive-Lane Poisoning & `powf`/`pow`/`div` in `Dual<N>`, `Dual2<N>` & `ArenaTape` (`crates/ql-ad/src/dual.rs`, `tape.rs`)
- **File**: `crates/ql-ad/src/dual.rs`
  1. **`Dual::scale` (`lines 58-64`)**:
     ```rust
     #[inline]
     fn scale(&self, v: f64, s: f64) -> Self {
         let mut d = [0.0; N];
         for (o, &x) in d.iter_mut().zip(self.d.iter()) {
             *o = if x == 0.0 { 0.0 } else { x * s };
         }
         Self { val: v, d }
     }
     ```
  2. **`Dual2::unary` (`lines 334-344`)**:
     Guard zero `d1[i]`, `d2[i][j]`, and `d1[i] * d1[j]` before multiplying by `f1` or `f2` so `0.0 * INFINITY` never produces `NaN` in inactive lanes.
  3. **`Dual::powf` (`lines 239-241`) and `Dual2::powf` (`lines 515-525`)**:
     Handle special exponents cleanly without evaluating negative powers of `0.0`:
     - If `n == 0.0`: return `Self::one()`.
     - If `n == 1.0`: return `self`.
     - If `n == 2.0`: return `self * self`.
     - Otherwise evaluate `f1 = if self.val == 0.0 && n > 1.0 { 0.0 } else { n * self.val.powf(n - 1.0) }` and `f2 = if self.val == 0.0 && n > 2.0 { 0.0 } else { n * (n - 1.0) * self.val.powf(n - 2.0) }`.
  4. **`Dual2::div` (`lines 487-489`)**:
     Replace `self * rhs.powf(-1.0)` with direct reciprocal `unary`:
     ```rust
     let inv = 1.0 / r.val;
     let inv2 = inv * inv;
     self * r.unary(inv, -inv2, 2.0 * inv2 * inv)
     ```
  5. **`PartialEq` / `PartialOrd` consistency (`lines 17, 30-35, 280, 295-300`)**:
     Remove `#[derive(PartialEq)]` and `impl Eq` on `Dual<N>` and `Dual2<N>`; implement `PartialEq` by comparing `self.val == other.val` so `a == b` iff `a.partial_cmp(&b) == Some(Ordering::Equal)`.
- **File**: `crates/ql-ad/src/tape.rs`
  1. In `ArenaTape::gradient_from` (`line 244`), guard `if bar == 0.0 { continue; }` and only add `bar * op.w_lhs` if `op.w_lhs != 0.0` (and similarly for `rhs`).
  2. Upgrade `NEXT_TAPE_ID` (`line 8`), `TapeVar::tape_id` (`line 13`), `ArenaTape::tape_id` (`line 46`), and `AdjointResult::tape_id` from `u16` to `u32` (`AtomicU32`). Note: `TapeVar<f64>` (`val: f64, idx: u32, tape_id: u32`) remains **exactly 16 bytes**!
  3. In `AdjointResult::gradient` (`lines 29-35`), check `assert!(var.idx == CONST_IDX || var.tape_id == self.tape_id, "cross-tape variable")` and add `pub fn try_gradient(&self, var: TapeVar<f64>) -> Result<f64, &'static str>`.
  4. In `TapeScope` (`lines 256-312`):
     - Fix `with_active`: do **not** silently start and leak an active tape when `had == false`; if no tape is active, panic with a clear message or require `start_recording`.
     - Add an RAII `pub struct TapeGuard` returned by `TapeScope::guard()` / `TapeScope::guard_with_capacity(cap)` whose `Drop` implementation restores the previous tape (or clears `ACTIVE_TAPE`), ensuring unwind/early-return safety.

### 2.2 Fix `ift_lift_1d` Primal Residual Contamination (`crates/ql-ad/src/ift.rs:4-11`)
- **File**: `crates/ql-ad/src/ift.rs:4-11`
  - Currently:
    ```rust
    pub fn ift_lift_1d<T: Real, F: Fn(T) -> T>(x_star: f64, residual: F, f_x: T, f_xx: f64) -> T {
        let x = T::from_f64(x_star);
        let d = -residual(x) / f_x;
        x + d - T::from_f64(0.5 * f_xx) * d * d / f_x
    }
    ```
    When `residual(x).val() != 0.0`, `d.val() != 0.0`, so `d * d` injects a spurious first-order derivative $2 d_0 d_1$ into `Dual<N>`!
  - Fix `ift_lift_1d` to:
    1. Guard `if f_x.val().abs() <= 1e-14 || !f_x.val().is_finite() { return T::from_f64(x_star); }`.
    2. Evaluate `let r_raw = residual(T::from_f64(x_star));` and perform one safeguarded primal Newton polish step if `r_raw.val().abs() > 1e-15`:
       ```rust
       let fx0 = f_x.val();
       if !fx0.is_finite() || fx0.abs() <= 1e-14 {
           return T::from_f64(x_star);
       }
       let r0 = residual(T::from_f64(x_star)).val();
       let x_polished = if r0.is_finite() && r0.abs() > 1e-15 {
           let step = (-r0 / fx0).clamp(-0.5 * x_star.abs().max(1.0), 0.5 * x_star.abs().max(1.0));
           let cand = x_star + step;
           let r_cand = residual(T::from_f64(cand)).val();
           if r_cand.is_finite() && r_cand.abs() < r0.abs() { cand } else { x_star }
       } else {
           x_star
       };
       let x = T::from_f64(x_polished);
       let r_eval = residual(x);
       let r_diff = r_eval - T::from_f64(r_eval.val()); // exact zero primal: r_diff.val() == 0.0
       let d = -r_diff / f_x;
       x + d - T::from_f64(0.5 * f_xx) * d * d / f_x
       ```
    Look what this achieves: `r_diff.val() == 0.0` identically, so `d.val() == 0.0` identically, `d * d` has zero first derivative (`Dual`) and exact second derivative (`Dual2`), and primal residual noise never contaminates AD sensitivities!

### 2.3 Replace `SimdReal::norm_cdf` 5-Term `erf` with High-Precision Vectorized `simd_erfc` (`crates/ql-ad/src/scalar.rs:228-293`, `crates/ql-math/src/distributions.rs:213-218`)
- **File**: `crates/ql-ad/src/scalar.rs:228-293`
  - Implement a double-precision vectorized `simd_erfc4(z: Simd<f64, 4>) -> Simd<f64, 4>` and `simd_erfc8(z: Simd<f64, 8>) -> Simd<f64, 8>` using Cody's (1969) / Cephes rational minimax approximation on $|z| \le 0.46875$, $0.46875 < |z| \le 4.0$, and $|z| > 4.0$ (where $\operatorname{erfc}(|z|) = \frac{e^{-z^2}}{|z|\sqrt{\pi}} R(1/z^2)$ avoids left-tail subtraction cancellation completely, and for $z < 0$, $\operatorname{erfc}(z) = 2.0 - \operatorname{erfc}(|z|)$), and define:
    $$\Phi(x) = \frac{1}{2}\operatorname{erfc}\!\left(-\frac{x}{\sqrt{2}}\right)$$
    so `Simd<f64, 4>::norm_cdf` and `Simd<f64, 8>::norm_cdf` match `f64::norm_cdf` to $< 10^{-14}$ relative error even in the deep left tail $x \in [-30, -5]$!
- **File**: `crates/ql-math/src/distributions.rs:147-218`
  - In `norm_cdf_f64` (`line 214`), replace `0.5 * (1.0 + erf_f64(x * FRAC_1_SQRT_2))` with `0.5 * libm::erfc(-x * std::f64::consts::FRAC_1_SQRT_2)` so Newton refinement in `inv_norm_cdf_f64` never suffers left-tail cancellation for $p < 10^{-15}$.

### 2.4 Fix `brent_solve`, `Date` Civil Overflow, `BitsetCalendar::adjust`, and `PcaDecomposition` Mean-Centering
- **File**: `crates/ql-math/src/solvers.rs:28`
  - Replace `if fa * fb > 0.0` with:
    ```rust
    if !fa.is_finite() || !fb.is_finite() || (fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0) {
        return Err(MathError::NotBracketed);
    }
    ```
- **File**: `crates/ql-time/src/date.rs:205-228` & `crates/ql-time/src/calendar.rs:209-211`
  - Perform internal arithmetic in `civil_from_days` and `days_from_civil` in `i64` so extreme `Date` values never overflow `i32`.
  - In `BitsetCalendar::adjust` (`calendar.rs:205-220`), guard `d.0 > -719_468` / `d.0 < 2_932_896` in the `while !self.is_business_day(d)` loops so `Date(i32::MIN)` or `Date(i32::MAX)` cannot loop forever.
- **File**: `crates/ql-dml/src/pca.rs:1-67`
  - Add `pub mean: Vec<f64>` to `PcaDecomposition`.
  - In `encode(&self, x: &[f64], p: usize)`, subtract `self.mean[i]` (`(x[i] - self.mean[i])`), and in `decode(&self, l: &[f64])`, add `self.mean[i]` back (`self.mean[i] + sum`), so `decode(&encode(x, self.dim))` reconstructs `x` exactly when `center == true` (`fit_classical_pca`, `fit_differential_pca_central`), while `center == false` (`fit_differential_pca_non_central`) has `self.mean = vec![0.0; d]` and is unchanged!

---

## WORKSTREAM III — Phase 1–9 & Part I (`A1–A12`) Mathematical & Algorithmic Remediations

### 3.1 Fix `A2` Lifted Rough Heston Double-Negation Sign Bug & $H = 0.5$ Singularity (`crates/ql-models/src/rough.rs`)
- **File**: `crates/ql-models/src/rough.rs:14-33, 62-124`
  1. **Fix double-negation in `phi1` calls (`lines 83, 89, 118`)**:
     `phi1(x)` at line 40 already computes $\frac{1 - e^{-x}}{x}$. Passing `-0.5 * g * dt`, `-g * dt`, and `-gi * tau_mid` evaluated $\frac{e^{+\gamma_i t} - 1}{\gamma_i t}$, overflowing to `NaN` for $n = 20$!
     Change lines 83, 89, 118 to pass **positive** arguments:
     ```rust
     let half_phi: Vec<T> = gamma
         .iter()
         .map(|&g| phi1(T::from_f64(0.5 * g) * dt) * 0.5)
         .collect();
     let full_phi: Vec<T> = gamma
         .iter()
         .map(|&g| phi1(T::from_f64(g) * dt))
         .collect();
     ...
     let kernel_mid = c
         .iter()
         .zip(gamma.iter())
         .map(|(&ci, &gi)| T::from_f64(ci) * phi1(T::from_f64(gi) * tau_mid))
         .fold(T::zero(), |acc, x| acc + x);
     ```
  2. **Fix $H = 0.5$ ($0/0$) singularity in `lifted_heston_kernel` (`lines 18-32`)**:
     When $H \to 0.5$ ($\alpha = H + 0.5 \to 1.0$, $\epsilon = 1 - \alpha \to 0$), $\Gamma(\alpha)\Gamma(2-\alpha) \to 1$, $\gamma_i \to 0$, and $c_i \to \frac{1}{n}$ (with $\sum_{i=1}^n c_i = 1$). When $(0.5 - h).abs() < 1e-8$, return `(vec![1.0 / n as f64; n], vec![0.0; n])` (and remove the `if p.n_factors <= 1` special-case bypass at line 62 so all $n \ge 1$ flow through the exact solver!).

### 3.2 Fix `A10` `DupireLocalVol::dupire_local_vol_ad` Constant-Derivative Bug (`crates/ql-models/src/local_vol.rs`)
- **File**: `crates/ql-models/src/local_vol.rs:8-45`
  - Also re-export `DupireLocalVol` cleanly or keep methods compatible.
  - In `dupire_local_vol_ad`:
    Currently, `let wt = w.d1[0]; let wy = w.d1[1]; let wyy = w.d2[1][1];` extracts `f64` constants and loses all derivatives of $w_T, w_y, w_{yy}$ in the returned `Dual2<2>`!
    Instead, compute `(wt, wy, wyy)` as `Dual2<2>` values by applying central 4th-order / step-`h = 1e-4` shifts in the $(T, y)$ arguments while preserving the incoming `Dual2<2>` tangents/Hessians of `(strike, expiry)` (or any seeded parameter lanes in `f`):
    ```rust
    pub fn dupire_local_vol_ad<F>(&self, strike: Dual2<2>, expiry: Dual2<2>, f: F) -> Dual2<2>
    where
        F: Fn(Dual2<2>, Dual2<2>) -> Dual2<2>,
    {
        let h = 1e-4;
        let dh = Dual2::from_value(h);
        let w = f(expiry, strike);
        let wt = (f(expiry + dh, strike) - f(expiry - dh, strike)) / (2.0 * h);
        let wy_p = f(expiry, strike + dh);
        let wy_m = f(expiry, strike - dh);
        let wy = (wy_p - wy_m) / (2.0 * h);
        let wyy = (wy_p - w * 2.0 + wy_m) / (h * h);
        let denom = Dual2::one() - strike / w * wy
            + Dual2::from_value(0.25)
                * (-Dual2::from_value(0.25) - Dual2::one() / w + strike * strike / (w * w))
                * wy
                * wy
            + Dual2::from_value(0.5) * wyy;
        (wt / denom).sqrt()
    }
    ```
    Why does this work for BOTH $(T, y)$ sensitivities and surface parameter sensitivities? Because `expiry + dh` shifts the primal value by `h` while keeping the `Dual2<2>` first and second derivative lanes intact, so `wt`, `wy`, and `wyy` are full `Dual2<2>` numbers carrying exact 1st and 2nd derivatives of $\partial_T w, \partial_y w, \partial_{yy} w$!

### 3.3 Fix `HestonParams::char_func_unit_spot`, `CosHestonEngine` Cumulants (`c1, c2, c4`), & `HestonQeSimulator` (`crates/ql-models/src/heston.rs`, `cos.rs`, `crates/ql-mc/src/paths.rs`)
- **File**: `crates/ql-models/src/heston.rs:37-49`
  - Implement Andersen–Lake cancellation-free `r = xi - d` and Little Heston Trap log stabilization (matching QuantLib `AnalyticHestonEngine::lnChF`, `repos/QuantLib/ql/pricingengines/vanilla/analytichestonengine.cpp:551-576`):
    ```rust
    let xi = Complex::from_real(self.kappa) - i * self.rho * self.sigma_v * u;
    let sig2 = self.sigma_v * self.sigma_v;
    let u2_plus_iu = u * u + i * u;
    let d = (xi * xi + Complex::from_real(sig2) * u2_plus_iu).sqrt();
    // Andersen-Lake cancellation avoidance: when Re(xi)*Re(d) + Im(xi)*Im(d) > 0,
    // rationalize r = xi - d = -sigma_v^2 * (u^2 + i*u) / (xi + d)
    let r = if (xi.re * d.re + xi.im * d.im).val() > 0.0 {
        (Complex::from_real(-sig2) * u2_plus_iu) / (xi + d)
    } else {
        xi - d
    };
    let ed = (-d * tau).exp();
    // y = (exp(-d*tau) - 1) / (2*d) = -0.5 * tau * phi1(d * tau)
    let y = (ed - Complex::one()) / (d * T::from_f64(2.0));
    let one_minus_ry = Complex::one() - r * y;
    let a = Complex::from_real(self.kappa * self.theta / sig2)
        * (r * tau - one_minus_ry.ln() * T::from_f64(2.0));
    let b = (u2_plus_iu * y) / one_minus_ry;
    (i * u * (rate - div) * tau + a + b * self.v0).exp()
    ```
    Notice how `r * y = (xi - d)(e^{-d\tau} - 1)/(2d) = \frac{\xi - d}{\xi + d} \frac{e^{-d\tau} - 1}{1 - g} \implies 1 - r y = \frac{1 - g e^{-d\tau}}{1 - g}$, AND in `b = (u^2 + i u) y / (1 - r y)`, $\sigma_v^2$ cancels out of the denominator completely!
- **File**: `crates/ql-models/src/cos.rs:100-195`
  - Pass `p.rho` into `price_impl` / `price_impl_extra`, and for Bates (`price_bates`, `cos.rs:110-126`), pass BOTH the jump drift adjustment `jump_c1 = p.lambda * (p.mu_j - ((p.mu_j + p.sigma_j * p.sigma_j * 0.5).exp() - 1.0)) * expiry` AND `jump_c2 = p.lambda * (p.mu_j * p.mu_j + p.sigma_j * p.sigma_j) * expiry` AND `jump_c4 = p.lambda * (p.mu_j.powf(4.0) + p.mu_j * p.mu_j * p.sigma_j * p.sigma_j * 6.0 + p.sigma_j.powf(4.0) * 3.0) * expiry`.
  - Compute the Fang–Oosterlee (2008, Table A.1) Heston cumulants $c_1, c_2$ (including the full $\sigma_v, \rho$ terms in $c_2$) and set the truncation half-width to `w = self.l_trunc * (c2.abs() + c4.max_val(T::zero()).sqrt()).max_val(T::from_f64(1e-12)).sqrt()`.
- **File**: `crates/ql-mc/src/paths.rs:254-300`
  1. In `qe_step` when `sigma_v.abs() < 1e-14` (`lines 263-267`):
     Compute exact integrated variance $v_{\text{int}} = \theta \Delta t + (v - \theta)\frac{1 - e^{-\kappa \Delta t}}{\kappa}$ (with $\Delta t$ limit when $\kappa \to 0$) and include the Brownian diffusion term:
     ```rust
     if sigma_v.abs() < 1e-14 {
         let v_int = if kappa.abs() < 1e-12 {
             v * dt
         } else {
             theta * dt + (v - theta) * (1.0 - e) / kappa
         }.max(0.0);
         let log_next = s + (self.rate - self.div - 0.5 * v_int / dt) * dt + v_int.sqrt() * z_x;
         return (log_next, m.max(0.0));
     }
     ```
  2. In `qe_step` for $\sigma_v > 0$ (`lines 268-298`), implement Andersen's (2008, Section 3.2.3, Eqs. 38–40) QE-M martingale-corrected drift $K_0^*$:
     Let $\gamma_1 = 0.5, \gamma_2 = 0.5$, $K_1 = \gamma_1 \Delta t \left(\frac{\kappa\rho}{\sigma_v} - \frac{1}{2}\right) - \frac{\rho}{\sigma_v}$, $K_2 = \gamma_2 \Delta t \left(\frac{\kappa\rho}{\sigma_v} - \frac{1}{2}\right) + \frac{\rho}{\sigma_v}$, $K_3 = \gamma_1 \Delta t (1 - \rho^2)$, $K_4 = \gamma_2 \Delta t (1 - \rho^2)$, and $A = K_2 + \frac{1}{2}K_4$.
     - In the quadratic branch ($\psi \le \psi_c$), if $1 - 2 A a > 10^{-12}$, set $K_0^* = -\frac{A b^2 a}{1 - 2 A a} + \frac{1}{2}\ln(1 - 2 A a) - (K_1 + \frac{1}{2}K_3) v$;
     - In the exponential branch ($\psi > \psi_c$), if $\beta > A + 10^{-12}$, set $K_0^* = -\ln\!\left(p + \frac{\beta(1 - p)}{\beta - A}\right) - (K_1 + \frac{1}{2}K_3) v$;
     - Update `s_next = s + (self.rate - self.div) * dt + k0_star + k1 * v + k2 * v_next + ((k3 * v + k4 * v_next).max(0.0)).sqrt() * z_x`.

### 3.4 Fix `ql-curves` Bootstrap Helpers, `BasisSwapInstrument`, and `FraInstrument` IFT Hessian (`crates/ql-curves/src/bootstrap.rs`)
- **File**: `crates/ql-curves/src/bootstrap.rs`
  1. In `RateInstrumentKind` (`lines 185-227`), add `Fra(FraInstrument)` and `BasisSwap(BasisSwapInstrument)` and match them in `from_single_curve_helpers` (`lines 266-271`) and `from_multicurve_helpers` (`lines 343-349`) so `self.instruments.len() == k` holds when FRAs or basis swaps are included!
  2. In `BasisSwapInstrument::residual_two_curves` (`lines 606-624`), accumulate the discount annuity `annuity += p * 0.25` alongside `x` and return `x / annuity - T::from_f64(self.spread)`.
  3. In `bootstrap_piecewise_curve` and `bootstrap_global_multicurve` (`lines 146-171`), solve sequential par-swap / deposit residuals (where a pillar at $t_k \le 1.0$ is treated as a simple deposit $P(0, t_k) = \frac{1}{1 + q_k t_k}$ and pillars $t_k > 1.0$ solve the annual par-swap condition $q_k \sum_{j=1}^k \tau_j P(0, t_j) + P(0, t_k) = 1$) so par swap rates at the pillars reproduce the input quotes!

### 3.5 Fix `normal_sabr_vol`, Discrete Lookback & Barrier, Inflation & Commodity Extrapolation, Callable/Convertible Bonds, PDE Barriers, ASB CDO, and FX Vanna-Volga
- **File**: `crates/ql-rates/src/swaption.rs:14-20`
  - In `normal_sabr_vol`, replace the lognormal term `omb * omb * alpha * alpha / (fk.pow(omb) * 24.0)` with Hagan's Bachelier term:
    `-T::from_f64(beta * (2.0 - beta)) * alpha * alpha / (fk.pow(omb) * 24.0)`.
- **File**: `crates/ql-exotics/src/lookback.rs:127-181` & `crates/ql-exotics/src/barrier.rs:151-155`
  1. In `price_fixed_lookback_discrete` (`lookback.rs:151-181`), apply Broadie, Glasserman & Kou (1999, Section 4):
     Shift the **strike** $K_{\text{adj}} = K e^{+\eta \beta \sigma \sqrt{\Delta t}}$ (where $\eta = +1$ for a call on the maximum and $\eta = -1$ for a put on the minimum) and scale the continuous fixed-lookback result by $e^{-\eta \beta \sigma \sqrt{\Delta t}}$ (with adjusted running extremum $M_{\text{adj}} = \max(M_0, S_0) e^{+\eta \beta \sigma \sqrt{\Delta t}}$ for call / $\min(m_0, S_0) e^{+\eta \beta \sigma \sqrt{\Delta t}}$ for put):
     - Specifically, by BGK (1999, Theorem 2 & Section 4), for a discrete fixed-strike lookback call on the maximum:
       $$\mathbb{E}\!\left[e^{-rT}(\max_{0 \le k \le N} S_{t_k} - K)^+\right] \approx e^{-\beta \sigma \sqrt{\Delta t}} \, V_{\text{fixed, cont}}^{\text{call}}\!\left(S_0, \max(M_0, S_0) e^{+\beta\sigma\sqrt{\Delta t}}, K e^{+\beta\sigma\sqrt{\Delta t}}\right)$$
       and for a discrete fixed-strike lookback put on the minimum:
       $$\mathbb{E}\!\left[e^{-rT}(K - \min_{0 \le k \le N} S_{t_k})^+\right] \approx e^{+\beta \sigma \sqrt{\Delta t}} \, V_{\text{fixed, cont}}^{\text{put}}\!\left(S_0, \min(m_0, S_0) e^{-\beta\sigma\sqrt{\Delta t}}, K e^{-\beta\sigma\sqrt{\Delta t}}\right)$$
       Notice that now `K` is shifted to $K e^{\pm \beta \sigma \sqrt{\Delta t}}$, so `price_fixed_lookback_discrete` is **no longer a no-op** when $K > M_0$, converges monotonically to `price_fixed_lookback_continuous` as `observations -> inf`, and satisfies `discrete <= continuous`!
  2. In `price_floating_lookback_discrete` (`lookback.rs:127-147`), apply BGK (1999, Proposition 5):
     Using the duality between floating-strike and fixed-strike lookbacks ($S_T - m_N = (S_T - m_0) + (m_0 - m_N)^+$ for a call, $M_N - S_T = (M_0 - S_T) + (M_N - M_0)^+$ for a put):
     - For a call (minimum $m_0 \le S_0$):
       `let m0 = extremum.min_val(spot);`
       `let shift = (vol * dt.sqrt() * beta).exp();`
       Evaluate `price_floating_lookback_continuous(spot, (m0 * shift).min_val(spot), ...)` when `m0 * shift <= spot`, and when `m0 * shift > spot` (including $m_0 = S_0$), use the exact BGK Proposition 5 formula:
       $$V_{\text{float, disc}}^{\text{call}}(S_0, m_0) = e^{-\beta\sigma\sqrt{\Delta t}} V_{\text{float, cont}}^{\text{call}}(S_0, m_0 e^{\beta\sigma\sqrt{\Delta t}}) \quad \text{via } S_0 e^{-qT} - m_0 e^{-rT} + e^{+\beta\sigma\sqrt{\Delta t}} V_{\text{fixed, cont}}^{\text{put}}(S_0, S_0, m_0 e^{-\beta\sigma\sqrt{\Delta t}})$$
       Checking the identity: at $t=0$, $S_T - m_N = S_T - m_0 + (m_0 - m_N)^+$, so the floating lookback call equals the forward $S_0 e^{-qT} - m_0 e^{-rT}$ plus a **fixed-strike lookback put** with strike $m_0$! And for a floating lookback put ($M_N - S_T = M_0 - S_T + (M_N - M_0)^+$), it equals $M_0 e^{-rT} - S_0 e^{-qT}$ plus a **fixed-strike lookback call** with strike $M_0$!
       That duality identity holds for **all** $m_0 \le S_0$ (including $m_0 = S_0$), never violates the domain $m_0 \le S_0$, satisfies `discrete <= continuous`, and converges to `price_floating_lookback_continuous` as `observations -> inf`!
  3. In `crates/ql-exotics/src/barrier.rs:151-155`, check whether `spec.spot` has breached `spec.barrier` (the contractual barrier, NOT `adjusted`) before applying the BGK shift, and clamp `adjusted` so it cannot cross `spec.spot` if `spec.spot` has not breached `spec.barrier`.
- **File**: `crates/ql-inflation/src/lib.rs:55-69, 326-334` & `crates/ql-commodities/src/futures_curve.rs:44-63`
  - Fix left-extrapolation in `ZeroInflationCurve::z` (`ql-inflation/src/lib.rs:55-69`) and `CommodityFuturesCurve::futures_price` (`futures_curve.rs:44-63`):
    If `t <= times[0]`, return the first pillar value (`self.zero_rates[0]` / `(self.log_base_forwards[0] + self.seasonal.factor(t)).exp()`) directly without linearly extrapolating with `w < 0`.
  - In `price_yoy_inflation_capfloor` (`ql-inflation/src/lib.rs:326-334`), for `YoyVolType::Black`, the payoff is $\tau_i \cdot g \cdot \max(\pm (I_t / I_{t-1} - 1 - (K - s)/g), 0)$, so price on the shifted YoY rate or scale by `gearing` cleanly.
- **File**: `crates/ql-fx/src/lib.rs:187-196`
  - Replace the `if dprod.val().abs() < 1e-8 { first } else { ... }` branch in `FxVannaVolgaSmile::implied_vol` with the rationalized, branch-free $C^\infty$ identity:
    ```rust
    let dprod = d1 * d2;
    let d2corr = l1 * d1p * d2p * (self.sigma1 - self.sigma2).powf(2.0)
        + l3 * d1c * d2c * (self.sigma3 - self.sigma2).powf(2.0);
    let num = self.sigma2 * (first - self.sigma2) * 2.0 + d2corr;
    let rad = (self.sigma2 * self.sigma2 + dprod * num).max_val(T::from_f64(1e-16)).sqrt();
    self.sigma2 + num / (self.sigma2 + rad)
    ```
- **File**: `crates/ql-fixed-income/src/callable_bond.rs:104-179` & `crates/ql-fixed-income/src/convertible_bond.rs:17-82`
  - In `black_callable_bond_price` (`callable_bond.rs:106-118`), take the maximum embedded European call option value across call dates (and maximum put across put dates) instead of summing all European calls across every call date.
  - In `tree_callable_bond_price` (`callable_bond.rs:120-179`) and `price_convertible_bond_tf` (`convertible_bond.rs:28-30`), replace the half-closed/closed interval `<= dt * 0.5 + 1e-12` coupon lookup with exact step assignment (`let step_idx = ((ct / dt).val().round() as usize).min(n);` checked against `step_idx == j`) so a coupon is assigned to **at most one** tree step and never double-counted; and in `tree_callable_bond_price`, include the Ho–Lee convexity drift adjustment $\frac{1}{2}\sigma^2 t_j \Delta t$ (or normalize the step discount factor so an option-free zero-coupon bond reproduces $e^{-(r + \text{oas})T}$).
- **File**: `crates/ql-pde/src/adi.rs:167, 385-512` & `crates/ql-instruments/src/pde1d.rs:58-86, 134-160`
  - In `Fd2dHestonEngine::price_barrier` (`adi.rs:167`), pass `T::from_f64(o.rebate)` instead of `T::zero()`.
  - In both `Fd2dHestonEngine::solve` (`adi.rs`) and `Fd1dBlackScholesEngine::fd` (`pde1d.rs`), enforce the absorbing barrier rows (`a[i] = 0, b[i] = 1, c[i] = 0, rhs[i] = rebate`) **inside** the tridiagonal matrix before solving `solve_tridiagonal` / `sweep_x` so value never diffuses across knocked-out nodes during the implicit step!
- **File**: `crates/ql-mc/src/amc.rs:202-241`
  - Add `pub fn price_with_policy_ad_and_div<T: Real>(&self, ..., is_call: bool) -> T` (and keep `price_with_policy_ad` delegating to `is_call = false`) so both American calls and puts are supported with AD.
- **File**: `crates/ql-credit/src/cdo.rs:173-188`
  - Remove the `out * target_total / quadrature_total` rescaling in `asb_recursive_expected_tranche_loss` and clamp the returned expected tranche loss to `[0.0, d - a]`.

---

## WORKSTREAM IV — Phase 10 Round 1 Complete Mathematical & Algorithmic Remediation

### 4.1 `crates/ql-instruments/src/american.rs` — Distinct BAW (1987), BS (1993), BS (2002), and Ju–Zhong (1999) Engines + AD Implied Vol
1. **`BaroneAdesiWhaleyEngine`**:
   - Fix $r = 0$ division by zero (`m = 2.0 * rate / sig2`, `n = 2.0 * b / sig2`): notice that $M / (1 - e^{-rT}) = \frac{2 r / \sigma^2}{1 - e^{-rT}} \xrightarrow{r \to 0} \frac{2}{\sigma^2 T}$ is completely nonsingular! Define:
     $$h_K = \text{if } |r T| < 10^{-8} \text{ then } \frac{2}{\sigma^2 T} \text{ else } \frac{2 r}{\sigma^2 (1 - e^{-r T})}$$
   - Implement both the **Call** free boundary ($q_2 = \frac{-(n - 1) + \sqrt{(n - 1)^2 + 4 h_K}}{2} > 1$, $S^* > K$) AND the native **Put** free boundary ($q_1 = \frac{-(n - 1) - \sqrt{(n - 1)^2 + 4 h_K}}{2} < 0$, $S^{**} < K$ solving $K - S^{**} = P_{\text{euro}}(S^{**}) - \frac{S^{**}}{q_1}(1 - e^{(b-r)T}\Phi(-d_1(S^{**})))$).
   - Fix `dres` (`line 68`): for the call residual $R(S^*) = S^* - K - C_E(S^*) - \frac{S^*}{q_2}(1 - e^{-q T}\Phi(d_1(S^*)))$, since $\frac{\partial C_E}{\partial S^*} = e^{-qT}\Phi(d_1)$ and $\frac{\partial d_1}{\partial S^*} = \frac{1}{S^* \sigma \sqrt{T}}$, the factors of $S^*$ **cancel** in $\frac{S^*}{q_2} \frac{\partial \Phi(d_1)}{\partial S^*} = \frac{\phi(d_1)}{q_2 \sigma \sqrt{T}}$! Thus:
     $$\frac{\partial R_{\text{call}}}{\partial S^*} = \left(1 - \frac{1}{q_2}\right)\left(1 - e^{-qT}\Phi(d_1)\right) + \frac{e^{-qT}\phi(d_1)}{q_2 \sigma \sqrt{T}}$$
     (and for the put residual $R_{\text{put}}(S^{**}) = K - S^{**} - P_E(S^{**}) + \frac{S^{**}}{q_1}(1 - e^{-qT}\Phi(-d_1(S^{**})))$, $\frac{\partial R_{\text{put}}}{\partial S^{**}} = -\left(1 - \frac{1}{q_1}\right)\left(1 - e^{-qT}\Phi(-d_1)\right) + \frac{e^{-qT}\phi(d_1)}{q_1 \sigma \sqrt{T}}$). Also pass the second derivative `f_xx` into `ift_lift_1d` so `Dual2` Gamma/Vanna/Volga are exact!
2. **`BjerksundStensland1993Engine`**:
   - Implement the genuine Bjerksund & Stensland (1993) single flat boundary $I = B_0 + (B_\infty - B_0)(1 - e^{h(T)})$ closed-form formula with $\alpha = (I - K) I^{-\beta}$, $\beta = \left(\frac{1}{2} - \frac{b}{\sigma^2}\right) + \sqrt{\left(\frac{b}{\sigma^2} - \frac{1}{2}\right)^2 + \frac{2r}{\sigma^2}}$, and $\phi(S, T, \gamma, H, I)$ (using Bjerksund–Stensland put-call transformation $P(S, K, T, r, b, \sigma) = C(K, S, T, r - b, -b, \sigma)$ for puts).
3. **`BjerksundStensland2002Engine`**:
   - Implement the genuine Bjerksund & Stensland (2002) two-step boundary $(t_1 = \frac{1}{2}(\sqrt{5}-1)T, I_1, I_2)$ formula using `bivariate_normal_cdf` ($\Psi$) and $\phi$ functions (matching `repos/QuantLib/ql/pricingengines/vanilla/bjerksundstenslandengine.cpp`).
4. **`JuQuadraticEngine`**:
   - Implement the genuine Ju (1999) / Ju & Zhong (1999, Paper 41) quadratic approximation with the second-order correction $\frac{A (S/S^*)^{q_2}}{1 - b_{\text{Ju}} (\ln(S/S^*))^2 - c_{\text{Ju}} \ln(S/S^*)}$ where $(b_{\text{Ju}}, c_{\text{Ju}})$ are given by Ju (1999, Eqs. 23–26) (matching `repos/QuantLib/ql/pricingengines/vanilla/juquadraticengine.cpp`).
5. **`american_implied_vol`**:
   - Compute `f_x` (Vega) and `f_xx` (Volga) at `root` using `Dual2<1>::variable(root, 0)` instead of finite differences!

### 4.2 `crates/ql-models/src/hybrid.rs` — Exact BSM–Hull–White & Grzelak–Oosterlee (2011) H1-HW COS Engine
1. **`AnalyticBsmHullWhiteEngine::effective_variance` (`lines 8-23`)**:
   - Fix line 14–15 by adding `/ 2.0`:
     ```rust
     let rate_var = if k.val().abs() < 1e-6 {
         eta_r * eta_r * t * t * t / 3.0 - eta_r * eta_r * k * t.powf(4.0) * 0.25
     } else {
         eta_r * eta_r / (k * k * k * 2.0)
             * (k * t * 2.0 + (-k * t).exp() * 4.0 - (-k * t * 2.0).exp() - 3.0)
     };
     let cross = if k.val().abs() < 1e-6 {
         rho_sr * vol_s * eta_r * t * t - rho_sr * vol_s * eta_r * k * t * t * t / 3.0
     } else {
         rho_sr * vol_s * eta_r / (k * k) * (k * t - T::one() + (-k * t).exp()) * 2.0
     };
     ```
2. **`h1_hw_char_func_unit_spot` & `CosHestonHullWhiteEngine` (`lines 58-165`)**:
   - Implement Grzelak & Oosterlee (2011, Paper 35, Section 3.1) / QuantLib `AnalyticH1HWEngine` (`repos/QuantLib/ql/pricingengines/vanilla/analytich1hwengine.cpp`):
     - Under the $T$-forward measure for log-spot $x_T = \ln(S_T / S_0)$, the Hull–White integrated rate variance is:
       $$\sigma_{\text{HW}}^2(T) = \frac{\eta_r^2}{2 \kappa_r^3}\left(2\kappa_r T + 4 e^{-\kappa_r T} - e^{-2\kappa_r T} - 3\right)$$
     - The equity-rate correlation term $I_4(u, T)$ using Grzelak–Oosterlee's Gaussian moment-matching of $\mathbb{E}[\sqrt{v_t}] \approx c + d e^{-e t}$ (with exact mean $\bar{v}_{\text{sqrt}}(t) = \sqrt{\theta_v + (v_0 - \theta_v)e^{-\kappa_v t} - \frac{\xi_v^2}{8\kappa_v}(1 - e^{-\kappa_v t})}$ or effective spot volatility $\bar{\sigma}_S = \sqrt{\theta_v + (v_0 - \theta_v)\frac{1 - e^{-\kappa_v T}}{\kappa_v T}}$) adds the coupled equity-rate covariance:
       $$C_{sr}(T) = \frac{2 \rho_{sr} \eta_r \bar{\sigma}_S}{\kappa_r^2}\left(\kappa_r T - 1 + e^{-\kappa_r T}\right)$$
     - Thus the total hybrid variance added to the Heston characteristic exponent is $\Sigma_{\text{hyb}}^2(T) = \sigma_{\text{HW}}^2(T) + C_{sr}(T)$, with martingale-preserving characteristic factor:
       $$\phi_{\text{hyb}}(u) = \phi_{\text{Heston}}(u) \cdot \exp\!\left(-\frac{1}{2}\Sigma_{\text{hyb}}^2(T)\,(u^2 + i u)\right)$$
       Notice that at $u = -i$, $u^2 + i u = (-i)^2 + i(-i) = -1 + 1 = 0$, so $\phi_{\text{hyb}}(-i) = e^{(r-q)T}$ (**exact martingale conservation**!), $\eta_r \to 0$ recovers `CosHestonEngine` continuously, $\xi_v \to 0, v_0 = \theta_v = \sigma_S^2$ recovers `AnalyticBsmHullWhiteEngine` **to machine precision**, and $\rho_{sr}$ has the exact Grzelak–Oosterlee cross-covariance effect!
     - Wire `CosHestonHullWhiteEngine::price` to call `CosHestonEngine { n_terms: self.n_terms, l_trunc: self.l_trunc }.price_impl_extra(|u| h1_hw_char_func_unit_spot(u, p, rate, div, expiry), ...)` with `extra_variance = Sigma_hyb^2(T)`.

### 4.3 `crates/ql-exotics/src/basket_rainbow.rs` & `cliquet_extensible.rs` — Stulz (1982), Zhang (1995), Li–Deng–Zhou (2008), Levy (1992), Cliquet & Longstaff (1990)
1. **`StulzRainbowEngine` (`basket_rainbow.rs:4-69`)**:
   - Transcribe the exact Stulz (1982) / QuantLib `StulzEngine` (`repos/QuantLib/ql/pricingengines/basket/stulzengine.cpp`) formulas for `euro_two_asset_max_call`, `euro_two_asset_min_call`, `euro_two_asset_max_put`, and `euro_two_asset_min_put`:
     - Let $\sigma = \sqrt{\sigma_1^2 + \sigma_2^2 - 2\rho\sigma_1\sigma_2}\sqrt{T}$, $\rho_1 = \frac{\sigma_1 - \rho\sigma_2}{\sigma/\sqrt{T}}$, $\rho_2 = \frac{\sigma_2 - \rho\sigma_1}{\sigma/\sqrt{T}}$, $D_1 = \frac{\ln(F_1/F_2) + \frac{1}{2}\sigma^2}{\sigma}$, $D_2 = D_1 - \sigma$, $d_{i,1} = \frac{\ln(F_i/K) + \frac{1}{2}\sigma_i^2 T}{\sigma_i\sqrt{T}}$, $d_{i,2} = d_{i,1} - \sigma_i\sqrt{T}$:
       - $C_{\min}(S_1, S_2, K) = e^{-rT}\left[F_1 M(d_{1,1}, -D_1; -\rho_1) + F_2 M(d_{2,1}, D_2; -\rho_2) - K M(d_{1,2}, d_{2,2}; \rho)\right]$
       - $C_{\max}(S_1, S_2, K) = C_{\text{BS}}(S_1, K) + C_{\text{BS}}(S_2, K) - C_{\min}(S_1, S_2, K)$
       - Forward values of min/max at $K \to 0$:
         $$PV(\min(S_1, S_2)) = e^{-rT}\left[F_1 \Phi(-D_1) + F_2 \Phi(D_2)\right], \qquad PV(\max(S_1, S_2)) = e^{-rT}(F_1 + F_2) - PV(\min(S_1, S_2))$$
       - Put-call parity on min and max:
         $$P_{\min}(S_1, S_2, K) = K e^{-rT} - PV(\min(S_1, S_2)) + C_{\min}(S_1, S_2, K)$$
         $$P_{\max}(S_1, S_2, K) = K e^{-rT} - PV(\max(S_1, S_2)) + C_{\max}(S_1, S_2, K)$$
2. **`TwoAssetCorrelationEngine` (`basket_rainbow.rs:71-101`)**:
   - Implement Zhang (1995) / Haug (2007, Section 4.11.1):
     - Let $y_1 = \frac{\ln(S_1/K_1) + (r - q_1 - \frac{1}{2}\sigma_1^2)T}{\sigma_1\sqrt{T}}$ and $y_2 = \frac{\ln(S_2/K_2) + (r - q_2 - \frac{1}{2}\sigma_2^2)T}{\sigma_2\sqrt{T}}$:
       - For a Call ($S_2 > K_2$ if $S_1 > K_1$):
         $$C = S_2 e^{-q_2 T} M(y_2 + \sigma_2\sqrt{T},\, y_1 + \rho\sigma_2\sqrt{T};\, \rho) - K_2 e^{-rT} M(y_2,\, y_1;\, \rho)$$
       - For a Put ($K_2 - S_2$ if $S_1 < K_1$):
         $$P = K_2 e^{-rT} M(-y_2,\, -y_1;\, \rho) - S_2 e^{-q_2 T} M(-y_2 - \sigma_2\sqrt{T},\, -y_1 - \rho\sigma_2\sqrt{T};\, \rho)$$
3. **`DengLiZhouSpreadBasketEngine` (`basket_rainbow.rs:103-138`)**:
   - Fix the double-compounding bug (`f1, f2` are already forwards, so discount by `df = (-r * t).exp()`, NEVER pass `f1` into `black_scholes_price(..., r, 0.0)`!) and implement the Li, Deng & Zhou (2008, Paper 36, Eqs. 15–24) second-order boundary approximation $J_0 + J_1 + J_2$ (and for $K = 0$, exact Margrabe formula).
4. **`LevyBasketEngine` (`basket_rainbow.rs:140-171`)**:
   - Fix the double-compounding bug and implement Levy (1992):
     $$M_1 = \sum_{i=1}^n w_i F_i, \qquad M_2 = \sum_{i=1}^n \sum_{j=1}^n w_i w_j F_i F_j e^{\rho_{ij}\sigma_i\sigma_j T}, \qquad \sigma_{\text{Levy}} = \sqrt{\frac{\ln(M_2 / M_1^2)}{T}}$$
     and return `black76_price(m1, strike, sigma_levy, expiry, (-rate * expiry).exp(), is_call)`.
5. **`AnalyticCliquetEngine` (`cliquet_extensible.rs:33-47`)**:
   - For each reset interval $\Delta t_i = t_i - t_{i-1} > 0$, under Black–Scholes with forward return ratio $m_i = e^{(r - q)\Delta t_i}$ and $1 + R_i = \frac{S_{t_i}}{S_{t_{i-1}}}$, the locally floored & capped return $\min(\max(R_i, F_{\text{loc}}), C_{\text{loc}})$ for a call (`is_call = true`) or put (`is_call = false`) decomposes into exact Black76 options on $1 + R_i$ with forward $m_i$, volatility $\sigma$, and maturity $\Delta t_i$:
     - For `is_call == true`, $\max(\min(R_i, C_{\text{loc}}), F_{\text{loc}}) = F_{\text{loc}} + \text{Call}_{76}(m_i, 1 + F_{\text{loc}}, \sigma, \Delta t_i) - \text{Call}_{76}(m_i, 1 + C_{\text{loc}}, \sigma, \Delta t_i)$, discounted by $e^{-r t_i}$ and multiplied by `notional`.
     - For `is_call == false`, $\max(\min(-R_i, C_{\text{loc}}), F_{\text{loc}}) = F_{\text{loc}} + \text{Put}_{76}(m_i, 1 - F_{\text{loc}}, \sigma, \Delta t_i) - \text{Put}_{76}(m_i, 1 - C_{\text{loc}}, \sigma, \Delta t_i)$, discounted by $e^{-r t_i}$ and multiplied by `notional`.
   - Also add `pub fn price_mc_global<T: Real>(s: &CliquetOptionSpec<T>, n_paths: usize, seed: u64) -> T` to price the path-dependent `global_floor` and `global_cap` exactly.
6. **`WriterExtensibleOptionEngine` (`cliquet_extensible.rs:49-65`)**:
   - Implement Longstaff (1990) / Haug (2007, Section 4.9.2):
     At $t_1$, if the option is OTM ($S_{t_1} < K_1$ for call, $S_{t_1} > K_1$ for put), it extends automatically to $t_2$ with strike $K_2$:
     - Call: $V_{\text{call}} = C_{\text{BS}}(S, K_1, t_1) + S e^{-q t_2} M(-d_1(K_1, t_1),\, d_1(K_2, t_2);\, -\sqrt{t_1/t_2}) - K_2 e^{-r t_2} M(-d_2(K_1, t_1),\, d_2(K_2, t_2);\, -\sqrt{t_1/t_2})$
     - Put: $V_{\text{put}} = P_{\text{BS}}(S, K_1, t_1) + K_2 e^{-r t_2} M(d_2(K_1, t_1),\, -d_2(K_2, t_2);\, -\sqrt{t_1/t_2}) - S e^{-q t_2} M(d_1(K_1, t_1),\, -d_1(K_2, t_2);\, -\sqrt{t_1/t_2})$.

### 4.4 `crates/ql-rates/src/cms_spread.rs` — Shifted-Lognormal Gauss–Hermite & Linear-TSR CMS Spread Pricer
- **File**: `crates/ql-rates/src/cms_spread.rs`
  - Compute the Linear-TSR convexity-adjusted forward swap rates $\mu_i = S_i + \Phi_{\text{TSR}}(S_i, \sigma_i, a, T_{\text{fix}}, T_{\text{pay}}, \text{tenor}_i)$ using `mean_reversion`, `payment_time`, `fixing_time`, `swap_tenor1`, `swap_tenor2`, `annuity1`, `annuity2`.
  - In `LognormalCmsSpreadPricer::lognormal_caplet_pv` (and `lognormal_floorlet_pv`), implement Brigo & Mercurio (2006, Section 13.14.2) / Cech (2014, Paper 42) / QuantLib `LognormalCmsSpreadPricer` (`repos/QuantLib/ql/experimental/coupons/lognormalcmsspreadpricer.cpp`) 16-node Gauss–Hermite quadrature over the shifted-lognormal state $S_2(u) = -\text{shift}_2 + (\mu_2 + \text{shift}_2)\exp(\sigma_2\sqrt{T_{\text{fix}}} u - \frac{1}{2}\sigma_2^2 T_{\text{fix}})$, integrating the conditional Black76 price on $g_1(S_1 + \text{shift}_1)$ given $u$.
  - Make `collar_pv` return `Self::caplet_pv(s, cap_strike) - Self::floorlet_pv(s, floor_strike)` (consistent with `ql-inflation`).

### 4.5 `crates/ql-fixed-income/src/cat_bond.rs` — Fix `erlang_tail`, `Proportional` Trigger, Genuine Compound-Poisson MC (`PoissonPareto` & `PoissonExponential`), and `cat_bond_z_spread`
- **File**: `crates/ql-fixed-income/src/cat_bond.rs`
  1. **Fix `erlang_tail` (`lines 45-58`) and `pa, pe` (`lines 108-131`)**:
     - `erlang_tail(x, mu, n)` is currently returning $1 - e^{-\mu x}\sum_{k=0}^{n-1} \frac{(\mu x)^k}{k!} = F_n(x) = \mathbb{P}(S_t \le x \mid N_t = n)$ (the **CDF**!) for $n \ge 1$, and returning $0$ for $n = 0, x > 0$!
     - Define `erlang_cdf<T: Real>(x: T, mu: T, n: usize) -> T`:
       For $n == 0$, $\mathbb{P}(S_t \le x \mid N_t = 0) = 1.0$ (since $S_t = 0 \le x$ for $x \ge 0$).
       For $n \ge 1$, $F_n(x) = 1 - e^{-\mu x}\sum_{k=0}^{n-1} \frac{(\mu x)^k}{k!}$.
       And the tail probability (exceedance probability) is $\bar{F}_n(x) = \mathbb{P}(S_t > x \mid N_t = n) = 1 - F_n(x)$ (`0.0` for $n = 0$).
     - Then at maturity/coupon time $t$:
       - `attachment_prob` accumulates $w_n \bar{F}_n(A)$, which is **`0.0` when $\lambda = 0$**!
       - `exhaustion_prob` accumulates $w_n \bar{F}_n(E)$, which is **`0.0` when $\lambda = 0$**!
       - For `NotionalRisk::Digital { threshold }`, surviving notional fraction is $\mathbb{P}(S_t \le \text{threshold}) = \sum_{n=0}^\infty w_n F_n(\text{threshold})$.
       - For `NotionalRisk::Proportional { attachment: A, exhaustion: E }`, the conditional expected surviving notional fraction given $N_t = n$ is:
         - For $n = 0$: $1.0$.
         - For $n \ge 1$:
           $$\mathbb{E}\!\left[\text{rem}(S_t) \mid N_t = n\right] = 1 - \frac{\mathbb{E}[(S_t - A)^+ \mid N_t = n] - \mathbb{E}[(S_t - E)^+ \mid N_t = n]}{E - A}$$
           where for $S_t \sim \operatorname{Erlang}(n, \mu)$, $\mathbb{E}[(S_t - x)^+ \mid N_t = n] = \frac{n}{\mu}\bar{F}_{n+1}(x) - x \bar{F}_n(x)$!
           Look how clean and exact that formula is: when $S_t \ll A$, $\bar{F}_n(A) \to 0$ and $\bar{F}_{n+1}(A) \to 0$, so $\mathbb{E}[\text{rem}(S_t)] \to 1.0$!
  2. **Implement genuine compound-Poisson Monte Carlo in `price_mc` (`lines 155-157`)**:
     - Using `PhiloxRng::new(seed)` (from `ql_mc`), simulate `n_paths` compound-Poisson event trajectories up to $T_{\max}$ (sampling Poisson arrival counts or exponential inter-arrival times $-\ln(U_1)/\lambda$, and severities $X_j = -\mu_{\text{sev}}\ln(U_2)$ for `PoissonExponential` or $X_j = x_m U_2^{-1/\alpha}$ for `PoissonPareto`).
     - Also use a 64-path deterministic/QMC or compound-Poisson truncation in `price_analytical` when `peril` is `PoissonPareto` so `price_analytical` does not silently return `surv = 1.0`!
  3. **Fix `cat_bond_z_spread` (`lines 165-183`)**:
     - Discount each coupon/redemption cashflow at $t_i$ by $e^{-(r + z) t_i}$ inside the residual $R(z) = P(r + z) - P_{\text{mkt}}$, compute `f_x` and `f_xx` using `Dual2<1>::variable(root, 0)`, and lift via `ift_lift_1d(root, residual, f_x, f_xx)`.

### 4.6 `crates/ql-portfolio/src/allocation.rs` — Generic `<T: Real>` Linear Solver, Markowitz, Black–Litterman Posterior Covariance & Risk-Budgeting IFT
- **File**: `crates/ql-portfolio/src/allocation.rs`
  1. **Delete `solve_f64` (`lines 42-69`) and implement a generic `<T: Real>` linear solver**:
     ```rust
     pub fn solve_linear_generic<T: Real>(a: &[Vec<T>], b: &[T]) -> Result<Vec<T>, AllocationError>
     ```
     using Gaussian elimination with partial pivoting on `a[row][col].val().abs()`, returning `Err(AllocationError::SingularCovariance)` if `!pivot.val().is_finite() || pivot.val().abs() <= 1e-14` (NEVER clamping pivots with `.max(1e-14)`!).
  2. **Make `markowitz_minimum_variance_weights` and `markowitz_tangency_weights` (`lines 136-165`) genuinely generic in `<T: Real>`**:
     Call `solve_linear_generic(cov, &ones)` and `solve_linear_generic(cov, &excess)` directly in `T` so all `Dual`, `Dual2`, and `TapeVar` sensitivities with respect to `cov`, `expected_returns`, and `risk_free` propagate untouched!
  3. **Fix `black_litterman_posterior` (`lines 77-133`)**:
     - Perform all matrix operations in generic `<T: Real>` using `solve_linear_generic`!
     - Include `* tau` on the view update and compute the **updated Black–Litterman posterior covariance** (He & Litterman 1999, Paper 40):
       Let $\Sigma_\tau = \tau \Sigma$, $M = P \Sigma_\tau P^\top + \Omega$. Solve $M Y = P \Sigma_\tau$ and $M v = Q - P \Pi$:
       $$\mu_{\text{BL}} = \Pi + \Sigma_\tau P^\top v, \qquad \Sigma_{\text{BL}} = \Sigma + \Sigma_\tau - \Sigma_\tau P^\top Y$$
       (or posterior mean parameter covariance $\Sigma_{\text{BL}}^{\mu} = \Sigma_\tau - \Sigma_\tau P^\top Y$, with total return covariance $\Sigma_{\text{BL}} = \Sigma + \Sigma_{\text{BL}}^{\mu}$).
     - Solve $\Sigma_{\text{BL}} \tilde{w} = \frac{1}{\delta}\mu_{\text{BL}}$ (or normalized weights $\sum_i w_i = 1$) in generic `<T: Real>` so `posterior_returns`, `posterior_cov`, and `optimal_weights` carry exact `Dual`/`Dual2`/`TapeVar` derivatives!
  4. **Fix `risk_budgeting_weights` (`lines 201-278`)**:
     - After solving the primal positive risk-budgeting root $y^* \in \mathbb{R}_{++}^n$ via Spinu (2013) CCD on $F(y) = \frac{1}{2} y^\top \Sigma y - \sum_{i=1}^n b_i \ln y_i$ where $\nabla_y F(y^*) = \Sigma y^* - b / y^* = 0$ and $H_{ij} = \Sigma_{ij} + \delta_{ij} b_i / (y_i^*)^2$:
     - Lift $y(\Sigma, b)$ into generic `T: Real` via the Implicit Function Theorem on $\nabla_y F(y, \Sigma, b) = \Sigma y - b / y = 0$:
       Evaluate the residual at frozen $y^*$:
       $$r_i = \sum_{j=1}^n \Sigma_{ij} y_j^* - \frac{b_i}{y_i^*} \in T$$
       Solve $H \Delta y = -(r - \text{primal}(r))$ (using the Cholesky/LU factorization of $H$), set $y_i = T::\text{from\_f64}(y_i^*) + \Delta y_i$, and normalize $w_i = y_i / \sum_{k=1}^n y_k$ in `T`!
       Let's verify this mathematically!
       Because $y^*$ satisfies $\Sigma y^* - b / y^* = 0$, the exact derivative of $y^*(\Sigma, b)$ with respect to ANY parameter $\theta$ (whether in $\Sigma$ or in $b$) is:
       $$H \frac{\partial y^*}{\partial \theta} = -\left(\frac{\partial \Sigma}{\partial \theta} y^* - \frac{1}{y^*}\frac{\partial b}{\partial \theta}\right) = -\frac{\partial r(y^*; \Sigma, b)}{\partial \theta}$$
       And normalizing $w_i = \frac{y_i}{\sum_k y_k}$ in `T` automatically applies the quotient rule $\frac{\partial w}{\partial y}$!
       Let's double check: does $y^*$ depend on the scale of $\Sigma$ or $b$? No, $F(y)$ is strictly convex on $\mathbb{R}_{++}^n$ with **non-singular SPD diagonal-plus-SPD Hessian** $H = \Sigma + \operatorname{diag}(b_i / (y_i^*)^2) \succ 0$! Unlike the normalized KKT residual on $w$ (whose Jacobian is non-symmetric), the Spinu log-barrier optimality condition $\Sigma y - b / y = 0$ has a **strictly positive-definite symmetric Hessian $H = \Sigma + \operatorname{diag}(b_i / (y_i^*)^2)$** and gives the **100% exact** IFT tangent and adjoint for all $\Sigma_{ij}$ and all $b_i$ in $O(n^3)$!

### 4.7 Fix Benchmark Suite 19 & Suite 20 Integrity (`crates/ql-benches/src/suites.rs:1500-1630`)
- **File**: `crates/ql-benches/src/suites.rs`
  1. **Suite 19 (`benchmark_american_analytics`)**:
     - Baseline: 256-step CRR binomial tree with central finite-difference Delta, Gamma, Vega across a 64-option strip (32 calls, 32 puts).
     - SOTA: `JuQuadraticEngine::price` on `Dual2<2>` (seeding `spot` in lane 0 and `vol` in lane 1 to extract Price, Delta, Gamma, Vanna, Volga in a single pass!).
     - Both return a weighted checksum of `(price, delta, gamma, vega)` across all 64 options with verified parity (`tolerance = 0.05`).
  2. **Suite 20 (`benchmark_risk_parity`)**:
     - Baseline: Bump-and-re-solve Spinu CCD for all $n$ risk budgets $b_k \pm h$ and all $n(n+1)/2$ upper-triangle covariance entries $\Sigma_{ij} \pm h$ ($n = 12$, so $12 + 78 = 90$ bumped CCD solves), accumulating a weighted checksum of the weights $w_i$ AND the full Jacobian $\frac{\partial w_i}{\partial b_k}$ and $\frac{\partial w_i}{\partial \Sigma_{ij}}$.
     - SOTA: Single Spinu CCD primal solve + cached Cholesky/LU of $H = \Sigma + \operatorname{diag}(b_i / (y_i^*)^2)$ evaluating the exact IFT sensitivity matrix $\frac{\partial w}{\partial b}$ and $\frac{\partial w}{\partial \Sigma}$ (via `risk_budgeting_weights` on `Dual` or `TapeVar`), returning the **exact same weighted checksum** of weights + sensitivities with verified parity `< 1e-5`!

---

## MANDATORY VERIFICATION GATES BEFORE COMPLETION

Run and pass all of the following commands with **zero warnings and zero errors**:
1. `cargo fmt --all -- --check`
2. `cargo clippy --workspace --all-targets -- -D warnings`
3. `cargo test --workspace --all-targets` (must complete in $< 30\text{s}$ and pass all existing + new regression suites!)
4. `QL_WRITE_BENCH_REPORT=1 cargo run --release -p ql-benches` (verify all 20 benchmark suites pass parity and update `crates/ql-benches/BENCHMARK_RESULTS.md` and `knowledge_base/08_BENCHMARK_RESULTS.md`).
