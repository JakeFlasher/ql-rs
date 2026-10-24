# Phase 10 Adversarial Remediation Prompt — Round 4 (Pass 4 Verification of `9a72cde..e3a335a`)

## Context & Audit Reference Documents
Read the updated Pass 4 adversarial review documents in `knowledge_base/` before starting:
1. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_SECOND_PASS.md` — Pass 4 Executive Summary & Resolution Matrix (`P3-01`..`P3-11` -> `P4-01`..`P4-05`)
2. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_MATH.md` — Mathematical & Numerical Audit (Phases 1–10 vs. QuantLib C++ & Primary Literature)
3. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_RUST_ARCH.md` — Rust Systems, AD/IFT, SIMD & Numerical Architecture Audit
4. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_ARCHITECTURE.md` — Big-Picture Architecture, Pricing Engine & Workflow Audit

Round 3 (`9a72cde..e3a335a`) cleanly resolved 8 of the 11 Pass 3 items (`P3-02` through `P3-06`, `P3-08`, `P3-10`, and `P3-11`), including `BjerksundStensland2002Engine`, `BaroneAdesiWhaleyEngine::baw_put`, `h1_hw_char_func_unit_spot`, `LognormalCmsSpreadPricer`, `StulzRainbowEngine`, `CatBondSpec::price_mc`, `MarketContext`, `CosHestonEngine`, `dupire_local_vol_ad`, `simd_erfc`, `DayCounter`, `DeskWorkstation`, and `PortfolioDag`.

Pass 4 identified **5 remaining / newly introduced defects (`P4-01` through `P4-05`)**. Resolve **all 5 items** below in full, add the specified regression tests, and ensure zero compiler/clippy warnings across the workspace.

---

## Mandatory Remediation Tasks (`P4-01` through `P4-05`)

### `[P4-01]` Fix `JuQuadraticEngine::price` Formula Regression (`denom_l`, `alpha_ve_h`, `alpha_over_h`, Direct Put $\phi = -1$) & `american_implied_vol` IFT Lift (`crates/ql-instruments/src/american.rs:398–450, 483–500`)

1. **`JuQuadraticEngine::price` (`crates/ql-instruments/src/american.rs:388–450`)**:
   Compare against Ju (1999, Eqs. 14–18) and QuantLib `repos/QuantLib/ql/pricingengines/vanilla/juquadraticengine.cpp:136–192`:
   - **`alpha_over_h` at `|h| < 1e-10` (`line 399`)**:
     Since $\alpha = \frac{2r}{\sigma^2}$ and $h = 1 - e^{-rT} = rT + O(r^2 T^2)$, the limit as $r \to 0$ ($h \to 0$) is $\frac{2}{\sigma^2 T} = \text{T::from\_f64(2.0) / (v2 * t)}$, **not** `T::from_f64(2.0) / v2`:
     ```rust
     let alpha_over_h = if h.val().abs() < 1e-10 {
         T::from_f64(2.0) / (v2 * t)
     } else {
         alpha / h
     };
     ```
   - **`alpha_ve_h` and `denom_l` (`lines 435–442`)**:
     - Because $\alpha = \frac{2r}{\sigma^2}$ and $V_{E,h} = \frac{1}{r e^{-rT}} \frac{\partial V_E(S^*, T)}{\partial T}$, their product cancels $r$ identically for all $r$:
       $$\alpha V_{E,h} = \frac{2}{\sigma^2 e^{-rT}} \frac{\partial V_E(S^*, T)}{\partial T} = \text{T::from\_f64(2.0) / (v2 * disc0) * d\_ve\_dt}.$$
       Remove the spurious `+ (alpha_over_h - T::from_f64(2.0) / (v2 * t)) * euro_star / disc0` term that was added on line 436 in `1b98410`.
     - In Ju (1999, Eqs. 17–18) and QuantLib (`juquadraticengine.cpp:183–186`), the denominator in $b$ and $c$ is $2\lambda + \beta - 1 = \phi\sqrt{(\beta - 1)^2 + 4\alpha/h} = \phi\,\text{disc}$, **not** $\lambda - 1$:
       ```rust
       let disc0 = (-r * t).exp();
       let alpha_lambda_prime = -alpha_over_h * alpha_over_h / disc;
       let alpha_ve_h = T::from_f64(2.0) / (v2 * disc0) * d_ve_dt;
       let denom_l = phi_t * disc; // = 2.0 * lambda + beta - 1.0
       let bju = (T::one() - h) * alpha_lambda_prime / (denom_l * 2.0);
       let cju = -((T::one() - h) / denom_l)
           * (alpha_ve_h / h_a + alpha_over_h + alpha_lambda_prime / denom_l);
       ```
   - **Direct Put Formula via $\phi = \pm 1$ (`lines 388–450`)**:
     Unlike Bjerksund–Stensland, Ju (1999) and BAW define $h = 1 - e^{-rT}$; pricing a put via `ju_call(strike, spot, div, rate, ...)` replaces $h = 1 - e^{-rT}$ with $1 - e^{-qT}$ whenever $r \ne q$ (and hits $h = 0$ when $q = 0$). Unify `JuQuadraticEngine::price` for both calls ($\phi = +1$, `critical_call`, $\lambda = q_2 = \frac{1 - \beta + \text{disc}}{2}$, early exercise if $S \ge S^*$) and puts ($\phi = -1$, `critical_put`, $\lambda = q_1 = \frac{1 - \beta - \text{disc}}{2}$, early exercise if $S \le S^*$) directly as in `juquadraticengine.cpp:136–192`:
     - Let `phi = if is_call { 1.0 } else { -1.0 }; let phi_t = T::from_f64(phi);`
     - `let lambda = ((T::one() - beta) + phi_t * disc) * 0.5;`
     - `let root = if is_call { Self::critical_call(...) } else { Self::critical_put(...) };`
     - In `residual`: `phi_t * (x - k) - euro - h_a_x`, where `h_a_x = phi_t * (x - k) - euro` at $x = S^*$ and inside `residual(x)` use `let d1 = ((x / k).ln() + (r - q + v2 * 0.5) * t) / vt; let fx = phi_t * (T::one() - (-q * t).exp() * (phi_t * d1).norm_cdf()) - lambda / x * (phi_t * (x - k) - black_scholes_price(x, k, v, t, r, q, is_call));`.
     - At the lifted `star`: `let h_a = phi_t * (star - k) - euro_star;`
     - `let d_ve_dt = star * (-q * t).exp() * pdf_s * sigma / (T::from_f64(2.0) * st) - phi_t * q * star * (-q * t).exp() * (phi_t * d1s).norm_cdf() + phi_t * r * k * (-r * t).exp() * (phi_t * d2s).norm_cdf();`
     - Early exercise check: `if phi * (s.val() - star.val()) >= 0.0 { return (phi_t * (s - k)).max_val(euro_s); }`.
2. **`american_implied_vol` (`crates/ql-instruments/src/american.rs:483–500`)**:
   - Line 498 currently passes `T::from_f64(jet.d1[0])` as `f_x: T` to `ift_lift_1d`, which has zero first derivatives with respect to `(spot, strike, rate, div, expiry)` and drops the 2nd-order cross-derivative term $-2 F_{\sigma\theta}\sigma^*_\theta / F_\sigma$ under `Dual2`.
   - Compute `f_x: T` in `T` at `root` (without duplicating `TapeVar` subgraphs when `T = TapeVar<f64>`): e.g., when `T::from_f64(1.0)` is a forward jet (or via `Dual<1>` / central difference in volatility when `TapeScope::try_with_active(|_| ()).is_none()`):
     ```rust
     let fx = if ql_ad::TapeScope::try_with_active(|_| ()).is_some() {
         T::from_f64(jet.d1[0])
     } else {
         let h = 1e-5;
         (JuQuadraticEngine::price(spot, strike, rate, div, T::from_f64(root + h), expiry, is_call)
             - JuQuadraticEngine::price(spot, strike, rate, div, T::from_f64(root - h), expiry, is_call))
             / (2.0 * h)
     };
     let fxx = jet.d2[0][0];
     Ok(ift_lift_1d(root, |x| price_t(x) - target_price, fx, fxx))
     ```

---

### `[P4-02]` Fix `Fd2dHestonEngine::price_knockout_rebate` Knocked-Out Boundary & Same-Grid Knock-In Parity (`crates/ql-pde/src/adi.rs:87–112, 192–240, 462–527`)

1. **`price_knockout_rebate` Knocked-Out Dirichlet Value (`crates/ql-pde/src/adi.rs:87–112, 512–527`)**:
   - In `price_knockout_rebate` (`lines 98–111`), the 6th argument to `self.solve(...)` is the Dirichlet value on the **knocked-out region**, while the 12th argument `Some(rebate)` is the terminal survival payoff on the **un-knocked-out region**. Passing `rebate` for both causes `price_knockout_rebate` to pay `rebate` whether or not the barrier is hit, double-counting the barrier-hit rebate in `price_barrier`!
   - Pass `T::zero()` as the 6th argument (`rebate` on the knocked-out boundary) in `price_knockout_rebate` (`line 104`), keeping `Some(rebate)` as the 12th argument (`terminal_rebate`).
   - In `solve` (`lines 512–527`), when `let Some(constant) = terminal_rebate`:
     ```rust
     let bc_low = if lower_barrier.is_some() {
         rebate
     } else if let Some(constant) = terminal_rebate {
         constant * df_r
     } else if is_call { ... };
     let bc_high = if upper_barrier.is_some() {
         rebate
     } else if let Some(constant) = terminal_rebate {
         constant * df_r
     } else if is_call { ... };
     ```
2. **Avoid Unused 4th Solve & Align `vanilla` Grid in `price_barrier` (`crates/ql-pde/src/adi.rs:192–240, 462–475`)**:
   - Move `let out = self.price_knockout_barrier(...)` inside the `UpAndOut | DownAndOut` match arm so `UpAndIn | DownAndIn` does not waste an unused 2D ADI PDE solve.
   - Pass an optional `grid_center: Option<f64>` (or `align_barrier: Option<f64>`) into `solve` so that in `price_barrier` for `UpAndIn | DownAndIn`, `vanilla` is solved on the **exact same `(nx | 1)` grid centered at `option.barrier.ln()`** as `out_no_rebate`, guaranteeing exact discrete In–Out parity (`in + out == vanilla` when `rebate == 0`, and `in + out == vanilla + discounted_rebate` when `rebate > 0`).

---

### `[P4-03]` Fix `Fd1dBlackScholesEngine` Shifted Boundary Payoffs, `barrier: T` AD Sensitivity & Same-Grid Knock-In Parity (`crates/ql-instruments/src/pde1d.rs:64–225`)

`crates/ql-instruments/src/pde1d.rs` was skipped in Round 3. Implement all three fixes in `Fd1dBlackScholesEngine`:
1. **Use Shifted Boundary Spots & Fix European `putb` (`pde1d.rs:74–77, 126–137`)**:
   - After shifting `xmin` so `barrier_level` lies on node `j` (`lines 74–77`), compute:
     ```rust
     let smin_actual = xmin.exp();
     let smax_actual = (xmin + n as f64 * dx).exp();
     ```
   - Use `smax_actual` in `callb` and `smin_actual` in `putb`, and include `- T::from_f64(smin_actual) * (-dividend * tnew).exp()` in the European `putb`:
     ```rust
     let putb = if let Some((b, rb)) = barrier && smin_actual <= b.val() {
         rb
     } else if is_call {
         T::zero()
     } else if american {
         (strike - T::from_f64(smin_actual)).max_val(T::zero())
     } else {
         (strike * (-rate * tnew).exp()
             - T::from_f64(smin_actual) * (-dividend * tnew).exp())
             .max_val(T::zero())
     };
     ```
2. **Preserve AD Sensitivity $\partial V / \partial B$ wrt `barrier: T` (`pde1d.rs:74–77, 182–191`)**:
   - Define `xmin_t: T`:
     ```rust
     let xmin_t = if let Some((barrier_level, _)) = barrier {
         let target = barrier_level.val().ln();
         let j = (((target - xmin) / dx).round() as usize).clamp(1, n - 1);
         xmin = target - j as f64 * dx;
         barrier_level.ln() - T::from_f64(j as f64 * dx)
     } else {
         T::from_f64(xmin)
     };
     ```
   - In the final quadratic Lagrange interpolation (`lines 185–187`), use `xmin_t`:
     ```rust
     let xm = xmin_t + T::from_f64((i - 1) as f64 * dx);
     let x0 = xmin_t + T::from_f64(i as f64 * dx);
     let xp = xmin_t + T::from_f64((i + 1) as f64 * dx);
     ```
     Now `x - x0 = spot.ln() - barrier_level.ln() + const`, so differentiating `price_barrier` with respect to `barrier: Dual<1>` produces a non-zero barrier Greek $\partial V / \partial B$ matching central finite differences!
3. **Same-Grid Knock-In Parity & Knock-In Rebate (`pde1d.rs:205–225`)**:
   - For `UpAndIn | DownAndIn`, solve `vanilla` on the same barrier-shifted grid `xmin` as `out_no_rebate` (`rebate = T::zero()`), and when `rebate.val().abs() > 0.0` add the survival rebate leg (terminal payoff `rebate`, knocked-out boundary `0.0`) so `in + out == vanilla` when `rebate == 0`.

---

### `[P4-04]` Eliminate Dead `TapeVar` Graph Recording in `price_complex_chooser` & Use Analytical BS Gamma `fxx` (`crates/ql-exotics/src/quanto_chooser.rs:76–130`, `crates/ql-exotics/src/compound.rs:84–101`)

1. **`price_complex_chooser` (`crates/ql-exotics/src/quanto_chooser.rs:76–130`)**:
   - Currently, `brent_solve` (`line 81`) calls `residual(T::from_f64(x)).val()` on every iteration, and `fxx` (`lines 125–128`) calls `residual(T::from_f64(...)).val()` 3 more times. Because `residual` captures `s.strike_call: T, s.vol: T, ...`, when `T = TapeVar<f64>` this records dozens of dead Black–Scholes subgraphs onto the active `ArenaTape`!
   - Define a scalar `f64` closure `residual_f64` for `brent_solve` (and the Newton fallback) using `.val()` on all parameters so zero `TapeVar` nodes are recorded during root-finding.
   - Replace the 3-point finite difference `fxx` (`lines 124–128`) with the exact analytical Black–Scholes call minus put gamma at `root` (using `d1c` and `d1p` already computed on lines 116–121):
     ```rust
     let fxx = ((-s.div_yield.val() * call_tau.val()).exp() * d1c.val().norm_pdf()
         / (root * s.vol.val() * call_tau.val().sqrt()))
         - ((-s.div_yield.val() * put_tau.val()).exp() * d1p.val().norm_pdf()
             / (root * s.vol.val() * put_tau.val().sqrt()));
     ```
2. **`price_compound_option` (`crates/ql-exotics/src/compound.rs:84–85`)**:
   - Replace the 3-point finite difference `fxx` (`lines 84–85`) with the exact analytical daughter Black–Scholes gamma at `root`:
     ```rust
     let fxx = ql_vol::black_scholes_gamma(
         root,
         k2.val(),
         spec.vol.val(),
         tau12.val(),
         spec.rate.val(),
         spec.div_yield.val(),
     );
     ```

---

### `[P4-05]` Fix `CalibratedCurveSystem::from_multicurve_helpers` DayCounter & `k > 64` Fallback (`crates/ql-curves/src/bootstrap.rs:334, 491`) and `hybrid.rs:126, 191` Small-$\kappa_r$ Cross-Integral Branch

1. **`crates/ql-curves/src/bootstrap.rs:332–335, 491–515`**:
   - On line 334, replace `DayCounter::Act365Fixed.year_fraction(reference_date, h.end, None, None)` with `h.day_counter.year_fraction(reference_date, h.end, None, None)` so `fwd_helpers` pillar times use each instrument's configured day counter (matching `IrsSwapInstrument::residual_with_discount` on lines 673, 679).
   - On line 491 (`if self.instruments.len() == k && k <= 64`), add an `else if self.instruments.len() == k` fallback branch that evaluates pairwise `Dual2<2>` entries `(i, jcol)` so systems with `k > 64` pillars also include the exact IFT constraint Hessian.
2. **`crates/ql-models/src/hybrid.rs:126–127, 191–192`**:
   - In `h1_hw_char_func_unit_spot` (`line 126`) and `CosHestonHullWhiteEngine::price` (`line 191`), `hw_cross_integral(p.kappa_r, c_go, t)` (`lines 73–78`) already implements the exact $\kappa_r \to 0$ Taylor limit `(c * t - 1 + (-c * t).exp()) / (c * c)`. Update the `p.kappa_r.val().abs() < 1e-8` branch on lines 127 and 192 to `p.rho_sr * p.eta_r * (a_go * t * t + b_go * hw_cross_integral(p.kappa_r, c_go, t) * 2.0)`.

---

## Verification & Test Checklist
1. **Unit & Integration Tests**:
   - Add regression tests verifying:
     - `JuQuadraticEngine::price` call and put values and `Dual<1>` / `Dual2<2>` sensitivities (including at $r \to 0$ and $r \neq q$), and `american_implied_vol` second-order `Dual2<1>` sensitivities against central finite differences.
     - `Fd2dHestonEngine::price_barrier` and `Fd1dBlackScholesEngine::price_barrier` exact discrete In–Out parity (`UpAndIn + UpAndOut == vanilla` at `rebate == 0`, and Knock-In rebate behavior when `rebate > 0`), plus non-zero `Dual<1>` barrier sensitivity $\partial V / \partial B$ in `Fd1dBlackScholesEngine::price_barrier`.
     - `price_complex_chooser` on `TapeVar<f64>` records only a single residual evaluation on `ArenaTape` and matches `Dual2<1>` second derivatives.
     - `CalibratedCurveSystem::from_multicurve_helpers` with `Act360` / `Thirty360` `IrsSwapInstrument` helpers reprices all calibration instruments to `< 1e-12`.
2. **Quality Gates**:
   - `cargo clippy --workspace --all-targets -- -D warnings` must pass with **0 warnings**.
   - `cargo test --workspace` must pass 100% of unit, integration, and parity tests.
