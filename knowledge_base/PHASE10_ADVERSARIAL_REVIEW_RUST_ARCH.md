# Phase 10 Comprehensive Adversarial Review — Rust Programming Practices, AD/IFT Safety, Linear Algebra & Numerical Portability

**Scope**: Line-by-line audit of `ql-ad`, `ql-math`, `ql-time`, `ql-curves`, `ql-portfolio`, `ql-xva`, `ql-pde`, `ql-mc`, `ql-scripting`, `ql-dml`, and AD-generic kernels across all 23 crates.

---

## 1. Forward & Second-Order AD (`Dual<N>`, `Dual2<N>`) & `Real` Trait Contract

### RUST-01 — **CRITICAL: Inactive Lane `0.0 * Inf = NaN` Poisoning in `Dual<N>`, `Dual2<N>`, and `TapeVar`**
- **Evidence**:
  - `crates/ql-ad/src/dual.rs`:
    - `Dual::scale` (`lines 58-64`): multiplies `self.d[i] * s`. When `self.d[i] == 0.0` (an unseeded/inactive directional derivative lane) and `s` is `INFINITY` or `-INFINITY` (e.g., `0.0.ln()`, `0.0.powf(0.5)`, `0.0.sqrt()`), IEEE-754 evaluates `0.0 * INFINITY = NaN`, poisoning all inactive derivative lanes.
    - `Dual2::unary` (`lines 334-344`): evaluates `self.d1[i] * f1` and `self.d2[idx] * f1 + self.d1[i] * self.d1[j] * f2` without guarding zero `d1`/`d2` lanes against infinite `f1`/`f2`.
    - `Dual2::powf` (`lines 515-525`):
      ```rust
      fn powf(self, n: f64) -> Self {
          self.unary(
              self.val.powf(n),
              n * self.val.powf(n - 1.0),
              n * (n - 1.0) * self.val.powf(n - 2.0),
          )
      }
      ```
      When `self.val == 0.0` and `n == 1.0`, `n * (n - 1.0) * self.val.powf(n - 2.0)` evaluates `1.0 * 0.0 * 0.0_f64.powf(-1.0) = 0.0 * inf = NaN`! And when `n == 2.0` and `self.val == 0.0`, `n == 0.0` evaluates `0.0 * (-1.0) * 0.0.powf(-2.0) = NaN`!
    - `Dual::pow` (`line 244`), `Dual2::pow` (`lines 526-539`), and `TapeVar::pow` (`tape.rs:523-529`): unconditionally evaluate `self.ln()` (`(-inf)` at `self.val == 0.0` and `NaN` for `self.val < 0.0`) and `1.0 / self.val`.
    - `ArenaTape::gradient_from` (`crates/ql-ad/src/tape.rs:244`): `adj[ix] += bar * w` evaluates `0.0 * inf = NaN` when `bar == 0.0` (an unreached branch on the tape) and local weight `w` is infinite.
- **Required Remediation**:
  - In `Dual::scale`, `Dual2::unary`, and `ArenaTape::gradient_from`, guard zero derivative lanes (`if d == 0.0 { 0.0 } else { d * s }`, `if bar == 0.0 { continue; }`).
  - In `Dual::powf` and `Dual2::powf`, handle `n == 0.0`, `n == 1.0`, `n == 2.0` explicitly without evaluating negative powers of `0.0`.
  - In `Dual2::div` (`dual.rs:487-489`), replace `self * rhs.powf(-1.0)` (which calls `powf` 3 times per division!) with direct reciprocal derivatives `inv = 1.0 / rhs.val`, `f1 = -inv * inv`, `f2 = 2.0 * inv * inv * inv`.

### RUST-02 — **HIGH: Violated `PartialEq` / `PartialOrd` Axiom & Unsound `Eq` on `Dual<N>` and `Dual2<N>`**
- **Evidence (`crates/ql-ad/src/dual.rs:17, 30-35, 280, 295-300`)**:
  - `Dual<N>` and `Dual2<N>` derive `PartialEq` (comparing `val` AND all derivative arrays `d`/`d1`/`d2`) and implement `Eq` (unsound because `f64::NAN != f64::NAN`), while implementing `PartialOrd` by comparing ONLY `self.val.partial_cmp(&other.val)`.
  - Consequently, for `a = Dual::variable(1.0, 0)` and `b = Dual::from_value(1.0)`, `a.partial_cmp(&b) == Some(Ordering::Equal)` while `a == b` is `false`! This violates the Rust standard library `PartialOrd` contract (`a == b` iff `partial_cmp(a, b) == Some(Equal)`).
- **Required Remediation**:
  - Remove `Eq` from `Dual<N>` and `Dual2<N>`, and implement `PartialEq` to compare `self.val == other.val` so `PartialEq` and `PartialOrd` are mutually consistent (matching `TapeVar` in `tape.rs:19-26`).

---

## 2. Reverse-Mode Tape (`ArenaTape`, `TapeVar`, `TapeScope`) & Implicit Function Theorem (`ift.rs`)

### RUST-03 — **CRITICAL: `ArenaTape` 16-Bit ID Wraparound, Missing `tape_id` Check in `AdjointResult::gradient`, & `TapeScope` Memory Leak**
- **Evidence (`crates/ql-ad/src/tape.rs`)**:
  - `NEXT_TAPE_ID: AtomicU16` (`line 8`) wraps around after 65,535 tapes (reached in seconds during Monte Carlo or row-wise Jacobian evaluations). Upgrading `tape_id` from `u16` to `u32` (`NEXT_TAPE_ID: AtomicU32`) keeps `TapeVar<f64>` (`val: f64, idx: u32, tape_id: u32`) at **16 bytes** with zero padding while raising the wrap threshold by $65,536\times$ to $4.29 \times 10^9$!
  - `AdjointResult::gradient` (`lines 31-35`) indexes `self.adjoints[var.idx]` without asserting or checking `var.tape_id == self.tape_id`. Store `tape_id: u32` in `AdjointResult` and assert `var.tape_id == self.tape_id` (plus provide `try_gradient(&self, var: TapeVar) -> Result<f64, AdError>`).
  - `TapeScope::with_active` (`lines 268-274`): when called with no active tape (`had == false`), it creates `ArenaTape::with_capacity(1024)` and stores it in thread-local `ACTIVE_TAPE` **forever** without ever stopping it, causing every subsequent `TapeVar` operation on that thread to record into an unbounded zombie tape!
  - Add an RAII `TapeGuard` (`pub struct TapeGuard { prev: Option<ArenaTape> }` implementing `Drop`) so panics or early `Result` returns always restore/clean up the thread-local tape.

### RUST-04 — **CRITICAL: `ift_lift_1d` Injects Spurious Second-Order Residual Contamination & Divides by Unchecked Zero Jacobian**
- **Evidence (`crates/ql-ad/src/ift.rs:4-11`)**:
  ```rust
  pub fn ift_lift_1d<T: Real, F: Fn(T) -> T>(x_star: f64, residual: F) -> T {
      let x = T::from_f64(x_star);
      let r = residual(x);
      let f_x = residual(Dual::<1>::variable(x_star, 0)).d[0];
      let f_xx = residual(Dual2::<1>::variable(x_star, 0)).d2[0];
      let d = -r / f_x;
      x + d - d * d * (0.5 * f_xx / f_x)
  }
  ```
  - When `x_star` has a non-zero primal residual $r_0 = \text{residual}(x^*).\text{val()} \ne 0$ (even $10^{-12}$ from Brent, or $O(1)$ when a caller used `.unwrap_or(guess)`), `d = -r / f_x` has a non-zero primal value $d_0 = -r_0 / f_x$. Then when `T = Dual<N>`, the quadratic term `d * d` injects a spurious first derivative $2 d_0 d_1 = 2 \frac{r_0}{f_x^2} \frac{\partial r}{\partial \theta}$ into the first-order sensitivity!
  - Furthermore, if `f_x == 0.0` (e.g., deep OTM/ITM zero vega), `-r / f_x` silently produces `NaN`/`Inf`.
- **Required Remediation**:
  - Subtract the primal residual before forming `d`:
    ```rust
    let r_full = residual(x);
    let r0 = r_full.val();
    let r = r_full - T::from_f64(r0); // pure differential jet (r.val() == 0.0)
    ```
    Perform at most one safeguarded primal Newton polish step `x_refined = x_star - r0 / f_x` if `r0.abs() > 1e-14 && f_x.abs() > 1e-14`, then evaluate `r = residual(T::from_f64(x_refined)) - T::from_f64(residual(T::from_f64(x_refined)).val())`. Because `r.val() == 0.0`, `d = -r / f_x` has exact zero primal (`d.val() == 0.0`), `d * d` is purely second-order (`Dual2`), and first-order `Dual` derivatives are never contaminated by $r_0$!
  - Add `try_ift_lift_1d` and `try_ift_lift_nd` returning `Result<..., AdError>` when `!f_x.is_finite() || f_x.abs() <= 1e-14` or LU fails.

---

## 3. Silent `.val()` Derivative Stripping Across Quantitative Crates

### RUST-05 — **HIGH: Functions Claiming Generic `<T: Real>` Strip AD Derivatives via `.val()`**
- **Evidence**:
  1. `crates/ql-portfolio/src/allocation.rs:42-165`:
     - `solve_f64` (`lines 42-69`): hand-rolled `f64` Gaussian elimination that clamps negative/zero pivots via `.max(1e-14)` (corrupting non-SPD systems) and uses `partial_cmp().unwrap()` (panics on `NaN`).
     - `markowitz_minimum_variance_weights` (`lines 136-148`) and `markowitz_tangency_weights` (`lines 150-165`): strip `cov` and `expected_returns` to `f64` via `.val()`, call `solve_f64`, and wrap constants back via `T::from_f64(...)` — returning **zero AD derivatives** with respect to every return and covariance input!
     - `black_litterman_posterior` (`lines 89-133`): strips `prior_cov`, `p`, `q`, `omega` to `f64` via `.val()`, computes in `f64` `nalgebra`, and wraps constants back in `T::from_f64(...)` — returning **zero AD derivatives** for posterior returns and weights!
     - `risk_budgeting_weights` (`lines 259-270`): evaluates `grad` at frozen `f64` `cov_f` and `b_f` instead of evaluating the KKT/Focal residual $r_i(w(\theta), \Sigma(\theta), b(\theta)) = w_i (\Sigma w)_i - b_i (w^\top \Sigma w)$ in generic `T` and lifting via $-H^{-1} r_T$.
  2. `crates/ql-exotics/src/variance_swap.rs:21-27`:
     - `replicate_variance_swap_fair_strike` checks `if quotes.windows(2).all(|w| (w[0].vol.val() - w[1].vol.val()).abs() < 1e-12) { return quotes[0].vol * quotes[0].vol; }`. When a flat vol surface is seeded with per-strike `Dual` or `TapeVar` variables to compute strike-by-strike Vega buckets, this shortcut attributes 100% of Vega to `quotes[0]` and **0.0** to `quotes[1..]`! Always execute the Carr–Madan log-contract replication loop (using analytic log-contract adjustment at the boundary).
  3. `crates/ql-levy/src/cgmy.rs:12-19`:
     - `CgmyProcess::gamma_neg_y` evaluates `gamma_lanczos(-y.val())` and lifts only a first-order finite difference via `T::from_f64(g0) + (y - T::from_f64(yv)) * (-dg)`, dropping all second-order `Dual2` curvature in $Y$. Make `gamma_lanczos<T: Real>(z: T) -> T` generic over `T: Real`!
  4. `crates/ql-ad/src/complex.rs:64-68`:
     - `Complex::arg` returns `T::from_f64(PI * 0.5)` (constant, zero derivative!) when `self.re.val() == 0.0` and `self.im.val() > 0.0`, dropping $-\frac{\text{d}(\text{re})}{\text{im}}$. Use smooth `atan2` derivative $\frac{\text{im}\cdot d(\text{re}) - \text{re}\cdot d(\text{im})}{\text{re}^2 + \text{im}^2}$.
  5. **Zero-limit branch derivative dropping (`if param.val().abs() < tol`)**:
     - `HullWhite1F::b` (`crates/ql-models/src/hull_white.rs:13-17`): returns `dt` when `|a| < 1e-10`, dropping $-\frac{1}{2} a\,\Delta t^2 + \frac{1}{6} a^2 \Delta t^3$.
     - `SchwartzSmith2F::b` (`crates/ql-commodities/src/schwartz_smith.rs:16-22`): returns `t` when `|z| < 1e-12`, dropping $-\frac{1}{2} z t^2 + \frac{1}{6} z^2 t^3$.
     - `G2ppModel::b` (`crates/ql-rates/src/g2.rs:29-35`): returns `dt` when `|k| < 1e-10`, dropping $-\frac{1}{2} k \Delta t^2 + \frac{1}{6} k^2 \Delta t^3$.
     - `VarianceGammaProcess::char_func_unit_spot` (`crates/ql-levy/src/vg.rs:28-35`): drops `nu` derivatives when `|nu| < 1e-12`.
     - `FxVannaVolgaSmile::implied_vol` (`crates/ql-fx/src/lib.rs:188-196`): branches on `if dprod.val().abs() < 1e-8 { first }`, dropping the entire second-order correction and creating a discontinuity at $d_1 = 0$ and $d_2 = 0$; replace with the rationalized identity $\sigma_2 + \frac{A}{\sigma_2 + \sqrt{\sigma_2^2 + d_1 d_2 A}}$.
  6. `crates/ql-curves/src/bootstrap.rs:265-279, 343-355, 450`:
     - `CalibratedCurveSystem::from_single_curve_helpers` and `from_multicurve_helpers` only match `DepositInstrument` and `OisSwapInstrument`/`IrsSwapInstrument`, silently ignoring `FraInstrument`!
     - `value_par_dv01_and_gamma` (`line 450`) uses fixed `Dual2::<64>` and silently drops sensitivities for pillars $i \ge 64$.

---

## 4. `SimdReal` Accuracy, Root-Finder Robustness, `Date` Arithmetic & Hot-Loop Allocations

### RUST-06 — **HIGH: `SimdReal::norm_cdf` Catastrophic Left-Tail Cancellation (`~1.5e-7` Error & Underflow at `x = -8.5`)**
- **Evidence (`crates/ql-ad/src/scalar.rs:228-293`)**:
  - `simd_norm_cdf` computes `0.5 * (1.0 + simd_erf(x / sqrt(2)))` using Abramowitz & Stegun 7.1.26 (5-term polynomial, max absolute error $1.5 \times 10^{-7}$).
  - For $x < -5$, `1.0 + simd_erf(x / sqrt(2))` subtracts two nearly equal numbers, losing all significant digits and hitting exact `0.0` at $x \approx -8.5$ (where true $\Phi(-8.5) \approx 9.48 \times 10^{-18}$).
  - Meanwhile, `f64::norm_cdf` (`line 111`) uses `0.5 * libm::erfc(-x / sqrt(2))`.
- **Required Remediation**:
  - Replace `simd_erf` in `simd_norm_cdf` with a vectorized double-precision rational `simd_erfc(-x / sqrt(2)) * 0.5` (Cody 1969 / Acklam / Cephes rational minimax approximation on $|z| \in [0, 0.46875], [0.46875, 4.0], [4.0, \infty)$ using $\text{erfc}(z) = \frac{e^{-z^2}}{z \sqrt{\pi}} R(1/z^2)$ for large $z>0$ and $2 - \text{erfc}(|z|)$ for $z<0$), achieving $< 10^{-15}$ relative accuracy across $x \in [-37, 8]$.

### RUST-07 — **HIGH: `brent_solve` Product Underflow & `NaN` Bypass (`crates/ql-math/src/solvers.rs:28`)**
- **Evidence**:
  - `if fa * fb > 0.0 { return Err(MathError::NotBracketed); }`
  - If `fa` or `fb` is `NaN`, `fa * fb > 0.0` is `false`, so `brent_solve` proceeds and returns a bogus root!
  - If `|fa| = 1e-160` and `|fb| = 1e-160` with the same sign, `fa * fb` underflows to `0.0` in IEEE-754 float64, bypassing the bracket check!
- **Required Remediation**:
  - Check `if !fa.is_finite() || !fb.is_finite() || (fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0) { return Err(MathError::NotBracketed); }`.

### RUST-08 — **MEDIUM: `BitsetCalendar::adjust` Infinite Loop at `Date(i32::MIN)` & Unchecked `Date` Civil Overflow**
- **Evidence**:
  - `crates/ql-time/src/calendar.rs:209-211`: `while !self.is_business_day(d) { d = Date(d.0 - 1); }` wraps or loops forever near `Date::MIN`.
  - `crates/ql-time/src/date.rs:205-228`: `civil_from_days` and `days_from_civil` perform unchecked `i32` arithmetic (`z + 719468`, `era * 146097`) that overflows on large `i32` inputs; upgrade internal conversion math to `i64` and bound valid `Date` range.

### RUST-09 — **MEDIUM: Inner-Loop Heap Allocations in Performance-Critical Engines**
- **Evidence**:
  - `mc_gbm_pathwise_tape` (`crates/ql-mc/src/paths.rs:375`): allocates a new `ArenaTape` with `8192` capacity on every single Monte Carlo path instead of `tape.clear()` reuse.
  - `Fd2dHestonEngine::solve` (`crates/ql-pde/src/adi.rs:189-499`): allocates 14 temporary `Vec<T>` per time step inside `for _ in 0..n_t`.
  - `AmcEngine::price_bermudan_lsm` (`crates/ql-scripting/src/amc.rs:85-113`): clones `Vec<VMValue>` and allocates `HashMap` per path per exercise date; replace string-keyed `HashMap` in `Env` with slice/indexed lookup.


---

# Pass 2 Post-Remediation Audit (`3559aac..19aa6ad`) — Remaining & Newly Introduced Rust Systems / AD / SIMD Defects

Following Codex's 9-commit remediation pass (`3559aac..19aa6ad`), a full line-by-line systems audit of `ql-ad`, `ql-time`, `ql-math`, `ql-models`, `ql-exotics`, `ql-fixed-income`, `ql-portfolio`, `ql-mc`, `ql-pde`, `ql-dml`, and `ql-xva` identified the following critical regressions, shortcuts, and unfixed items:

## `[RUST-P2-01]` Algebraic Self-Cancellation in `replicate_variance_swap_fair_strike` (`crates/ql-exotics/src/variance_swap.rs:53-60`)
- **Current Defect**:
  ```rust
  if flat_surface {
      let mean_sq = quotes
          .iter()
          .map(|q| q.implied_vol * q.implied_vol)
          .fold(T::zero(), |a, b| a + b)
          / quotes.len() as f64;
      variance += mean_sq - variance;
  }
  ```
  - `variance += mean_sq - variance` is algebraically `variance = mean_sq`! It cancels out 100% of the Carr–Madan strip integral `variance`, wiping out all AD derivatives with respect to `spot`, `rate`, `div_yield`, `expiry`, `cut_strike`, and `q.strike`, and replacing the true $\frac{\Delta K_i}{K_i^2}\mathcal{V}(K_i)$ vega weights with uniform $\frac{1}{N}$ weights!
- **Required Fix**:
  - Delete the `if flat_surface` overwrite completely.
  - Instead, use proper trapezoidal strip weights ($\Delta K_0 = \frac{1}{2}(K_1 - K_0)$, $\Delta K_{n-1} = \frac{1}{2}(K_{n-1} - K_{n-2})$) on $[K_0, K_{n-1}]$ and add the exact analytical Black–Scholes log-contract left put tail on $(0, K_0)$ (using $\sigma_0 = \text{quotes}[0].\text{implied\_vol}$, $v_0 = \sigma_0\sqrt{T}$):
    $$I_{\text{put}}(K_0, \sigma_0) = \frac{F}{K_0}\Phi(-d_1(K_0)) - \Phi(-d_2(K_0)) + v_0\left(-d_2(K_0)\Phi(-d_2(K_0)) + \phi(d_2(K_0))\right)$$
    and right call tail on $(K_{n-1}, \infty)$ (using $\sigma_{n-1} = \text{quotes}[n-1].\text{implied\_vol}$, $v_m = \sigma_{n-1}\sqrt{T}$):
    $$I_{\text{call}}(K_{n-1}, \sigma_{n-1}) = \frac{F}{K_{n-1}}\Phi(d_1(K_{n-1})) - \Phi(d_2(K_{n-1})) - v_m\left(d_2(K_{n-1})\Phi(d_2(K_{n-1})) + \phi(d_2(K_{n-1}))\right)$$
    multiplied by $\frac{2}{T}$. This is exact for any flat surface AND preserves true Carr–Madan vega/spot/rate sensitivities across all quotes!

## `[RUST-P2-02]` `DupireLocalVol::dupire_local_vol_ad` Replaced `Dual2<2>` Taylor Jet Extraction with Finite Differences (`crates/ql-models/src/local_vol.rs:8-30`)
- **Current Defect**:
  - `dupire_local_vol_ad` (`lines 16-22`) replaced exact `Dual2<2>` derivative extraction (`w.d1[0]`, `w.d1[1]`, `w.d2[1][1]`) with `h = 1e-4` central finite differences (`f(expiry + dh, strike)`, etc.), directly contradicting its own doc comment (*"no finite differences are taken"*) and introducing $O(h^2)$ truncation error into the primal local volatility.
- **Required Fix**:
  - Extract $\partial_T w = w.\text{d1}[0]$, $\partial_y w = w.\text{d1}[1]$, $\partial_{yy}^2 w = w.\text{d2}[1][1]$ from the exact `Dual2<2>` jet when `expiry` and `strike` are seeded as `(0, 1)` variables, and provide a generic `dupire_local_vol_surface_ad<T: Real, F: Fn(Dual2<2>, Dual2<2>) -> Dual2<2>>` or nested jet evaluation when differentiating wrt surface parameters.

## `[RUST-P2-03]` Wrong `4.0 / 3.0` Taylor Coefficient in `HullWhite1F` (`crates/ql-models/src/hull_white.rs:26, 47`) & Wrong $\kappa_v \to 0$ Limit in `ql-models::hybrid` (`crates/ql-models/src/hybrid.rs:78, 156`)
- **Current Defect**:
  1. `crates/ql-models/src/hull_white.rs:26, 47`: computes `tt - self.a * tt * tt + self.a * self.a * tt * tt * tt * 4.0 / 3.0` for $\frac{1 - e^{-2at}}{2a}$, where the cubic coefficient is $\frac{(2a)^3 t^3}{6 \cdot 2a} = \frac{2}{3} a^2 t^3$, **NOT** $\frac{4}{3} a^2 t^3$!
  2. `crates/ql-models/src/hybrid.rs:78, 156`: returns `p.theta_v` instead of `p.v0 - (p.v0 - p.theta_v) * p.kappa_v * t * 0.5` when `|kappa_v| < 1e-10`! As $\kappa_v \to 0$, $\bar{v}(t) = \theta_v + (v_0 - \theta_v)\frac{1 - e^{-\kappa_v t}}{\kappa_v t} \to v_0$, NOT $\theta_v$!

## `[RUST-P2-04]` Remaining `ql-ad` Defects: `ift_lift_1d` Duplicate Tape Recording, Unused `TapeGuard`, and `0.0 * Inf = NaN` in `propagate_range` / `propagate_to_mark` / `pow`
- **Current Defect**:
  1. `ift_lift_1d` (`crates/ql-ad/src/ift.rs:34-49`): evaluates `residual(T::from_f64(x_star)).val()` and `residual(T::from_f64(candidate)).val()` before `residual(x)`. When `T = TapeVar<f64>`, `residual` captures active `TapeVar<f64>` parameters and records the entire pricing graph onto `ArenaTape` **2 or 3 times**! Moreover, if `x_polished != x_star`, `f_x` and `f_xx` are still evaluated at `x_star`.
  2. `TapeScope::guard()` (`crates/ql-ad/src/tape.rs:281-342`) has **zero callers** across the workspace; all 10 production call sites still use `start_recording()` / `stop_recording()`.
  3. `ArenaTape::propagate_range` (`tape.rs:185-191`) and `propagate_to_mark` (`tape.rs:204-210`) lack the `if a == T::default() { continue; }` and `op.w_lhs != T::default()` guards present in `gradient_from`, propagating `0.0 * Inf = NaN` during checkpointed AAD (`mc_gbm_pathwise_tape`).
  4. `Dual::pow` (`dual.rs:301-312`), `Dual2::pow` (`dual.rs:672`), `TapeVar::pow` (`tape.rs:608`), and `TapeVar::powf` (`tape.rs:602`) still produce `NaN` at `self.val == 0.0`.

## `[RUST-P2-05]` Scalarized `simd_erfc4` / `simd_erfc8` (`crates/ql-ad/src/scalar.rs:234-240`) and Unfixed Hot-Loop Allocations (`RUST-09`)
- **Current Defect**:
  1. `simd_erfc4` and `simd_erfc8` (`scalar.rs:234-240`) call `Simd::from_array(x.to_array().map(libm::erfc))`, completely scalarizing `SimdReal::norm_cdf` into 4 or 8 sequential scalar function calls! Replace with a branchless vectorized minimax rational approximation (`Simd::simd_lt` / `.select(...)` over central $|z| \le 4$ and tail $|z| > 4$ multiplied by `simd_exp(-z * z)`), and in `simd_exp` (`scalar.rs:211`) extract exponent bits directly via `(u.to_bits() + Simd::splat(1023_u64)) << Simd::splat(52_u64)` without `k.cast::<i64>()`.
  2. None of the hot-loop heap allocations in `RUST-09` were fixed (`mc_gbm_pathwise_tape` at `paths.rs:404`, `Fd2dHestonEngine::solve` at `adi.rs:189-500`, `AmcEngine` at `amc.rs:17-27`, `cholesky_aad` at `cholesky_aad.rs:133-144`, `ChebyshevTensor2D::eval_barycentric` at `chebyshev.rs:158-180`, `DiffPolynomialRegressor::basis` at `regression.rs:26-55`).
  3. `Date::add_months` (`date.rs:162, 170`), `Period::months` (`date.rs:84-85`), `Date::add_period` (`date.rs:186, 189`), `day_count.rs:25, 79, 131`, and `inv_norm_cdf_f64` (`distributions.rs:198-202` where `pdf == 0.0` causes `0.0 / 0.0 = NaN` for `p < 1e-300`) remain unfixed.


---

# Pass 3 Rust Systems, AD/IFT & SIMD Audit (`19aa6ad..9a72cde`)

See [`PHASE10_ADVERSARIAL_REVIEW_SECOND_PASS.md`](./PHASE10_ADVERSARIAL_REVIEW_SECOND_PASS.md#pass-3-adversarial-audit-report-19aa6ad9a72cde) for the unified Pass 3 verification matrix. Key systems, AD/IFT, and SIMD findings in `19aa6ad..9a72cde`:

1. **`[RUST-P3-01]` (`crates/ql-instruments/src/american.rs:287-346`)**:
   - `bs_psi` (`lines 295, 308`) omits `* t` in `lambda.exp()`, swaps $t_1/T$ barriers, omits $H$ in $f_3, f_4$, uses $+\rho$ instead of $-\rho$ in $p_{13}, p_{14}$, and `bs2002_call` (`lines 336-343`) omits `s.pow(beta)` and `s` factors.
2. **`[RUST-P3-02]` (`crates/ql-ad/src/ift.rs:34-39`, `crates/ql-curves/src/bootstrap.rs:385-394, 498-512`)**:
   - `ift_lift_1d` divides unpolished `r_eval` by `f_x: T`, adding a spurious $+r_0 (\partial_\theta f_x)/f_x^2$ term when $r_0 \neq 0$ and $\partial_\theta f_x \neq 0$. Applying the primal polish in `f64` (`x_polished = x_star - r0 / f_x.val()`) and using `d_tan = -(r_eval - T::from_f64(r0)) / f_x` eliminates this contamination.
   - `bootstrap.rs` skips `BasisSwapInstrument` in `from_multicurve_helpers` (`line 393`) and assumes `k <= 64` in `Dual2::<64>` (`line 502`).
3. **`[RUST-P3-03]` (`crates/ql-models/src/hybrid.rs:58-70`, `crates/ql-models/src/local_vol.rs:17-19`)**:
   - `go_sqrt_moment_params` clamps `a` with `.max_val(1e-14)`, zeroing AD tangents when Feller is violated.
   - `dupire_local_vol_dual2_log_moneyness` wraps `w.d1[0], w.d1[1], w.d2[1][1]` in `Dual2::from_value(...)`, zeroing their `.d1` tangents $\nabla w_t, \nabla w_y$.
4. **`[RUST-P3-04]` (`crates/ql-ad/src/scalar.rs:345-351`, `crates/ql-time/src/day_count.rs:25, 107, 118, 129, 134`)**:
   - `simd_erfc` does not underflow to `0.0` for `ax >= 28.0` because `simd_exp` clamps at `-708.0`, and returns `NaN` for `x = ±Inf`.
   - `DayCounter::day_count`, `thirty_*`, and `Sub for Date` still use unchecked `i32` arithmetic.
5. **`[RUST-P3-05]` (`crates/ql-instruments/src/pde1d.rs:74-77, 126-139, 205-225`, `crates/ql-pde/src/adi.rs:174-197, 420-425`)**:
   - `Fd1dBlackScholesEngine` strips `barrier`'s AD tangent via `.val().ln()` and uses unshifted `smin`/`smax` in `callb`/`putb`.
   - `Fd2dHestonEngine::price_barrier` adds a Knock-Out rebate `+(out - out_no_rebate)` to Knock-In options and `HestonGrid2D::new_centered` with even `nx = 80` places `ln(B)` between nodes `39` and `40`.


---

## Pass 4 Systems, AD/IFT & SIMD Verification (`9a72cde..e3a335a`)

All Pass 3 systems/AD/SIMD defects in `ql-ad` (`ift_lift_1d` single residual evaluation, `simd_erfc` tail masking), `ql-models` (`dupire_local_vol_ad` tangents, `cos.rs` overflow-free $c_2$), `ql-curves` (`MarketContext::to_ad` active tape registration, `Dual2<64>` IFT Hessian), and `ql-time` (`day_count.rs` `i64` arithmetic) were verified as **RESOLVED** in `9a72cde..e3a335a`. Remaining items for Round 4:

1. **`P4-01` (`crates/ql-instruments/src/american.rs:398–402, 434–450, 483–500`)**: Fix `JuQuadraticEngine::price` `denom_l` (`2*q2 + beta - 1`), `alpha_ve_h`, `alpha_over_h` at `h -> 0`, direct put $\phi = -1$, and lift `f_x: T` in `american_implied_vol`.
2. **`P4-02` & `P4-03` (`crates/ql-pde/src/adi.rs` & `crates/ql-instruments/src/pde1d.rs`)**: Fix `price_knockout_rebate` knocked-out boundary (`T::zero()`), same-grid `vanilla` leg in 1D/2D Knock-In parity, shifted `smin_actual`/`smax_actual` boundary payoffs, and `barrier_level: T` AD sensitivity in `pde1d.rs`.
3. **`P4-04` (`crates/ql-exotics/src/quanto_chooser.rs:76–130` & `crates/ql-exotics/src/compound.rs:84–101`)**: Use a scalar `f64` closure inside `brent_solve` (`quanto_chooser.rs:81`) so root-finding does not record dead `TapeVar` subgraphs onto the active `ArenaTape`, and use exact analytical Black–Scholes gamma for `fxx` in both `quanto_chooser.rs` and `compound.rs`.
4. **`P4-05` (`crates/ql-curves/src/bootstrap.rs:334, 491` & `crates/ql-models/src/hybrid.rs:126, 191`)**: Use `h.day_counter` for `fwd_helpers` pillar times (`bootstrap.rs:334`), add a pairwise `Dual2<2>` fallback when `k > 64` (`bootstrap.rs:491`), and call `hw_cross_integral` in the `|kappa_r| < 1e-8` branch (`hybrid.rs:126, 191`).


---

## Pass 5 Systems, AD/IFT & Numerical Verification (`e3a335a..c2bef1a`)

`P4-04` (`quanto_chooser.rs` & `compound.rs` `TapeVar` isolation and analytical `fxx`), `P4-05` (`bootstrap.rs` `DayCounter` and `k > 64` pairwise `Dual2<2>` fallback, `hybrid.rs` `hw_cross_integral`), and `american_implied_vol` (`american.rs:547–568`) were verified as **RESOLVED** in `e3a335a..c2bef1a`. Remaining items for Round 5 are `P5-01` (`ju_put_direct` signs in `american.rs:481–490`), `P5-02` (`pde1d.rs` `UpAndOut` `putb` boundary, `smax`, and Knock-In rebate; `adi.rs` `grid_center.is_some()` odd-`nx` alignment and lazy `out` evaluation), and `P5-03` (`convertible_bond.rs` clean call/put order and step coupon accumulation).
