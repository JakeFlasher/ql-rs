# Phase 10 Adversarial Remediation Prompt — Round 5 (Pass 5 Verification of `e3a335a..c2bef1a`)

## Context & Audit Reference Documents
Read the updated Pass 5 adversarial review documents in `knowledge_base/` before starting:
1. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_SECOND_PASS.md` — Pass 5 Executive Summary & Resolution Matrix (`P4-01`..`P4-05` -> `P5-01`..`P5-03`)
2. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_MATH.md` — Mathematical & Numerical Audit (Phases 1–10 vs. QuantLib C++ & Primary Literature)
3. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_RUST_ARCH.md` — Rust Systems, AD/IFT, SIMD & Numerical Architecture Audit
4. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_ARCHITECTURE.md` — Big-Picture Architecture, Pricing Engine & Workflow Audit

Round 4 (`e3a335a..c2bef1a`) resolved `JuQuadraticEngine::ju_call`, `american_implied_vol`, `price_complex_chooser`, `price_compound_option`, `CalibratedCurveSystem::from_multicurve_helpers` & `value_par_dv01_and_gamma` (`k > 64`), `hybrid.rs` small-$\kappa_r$ `hw_cross_integral`, `Fd2dHestonEngine::price_knockout_rebate`, and `Fd1dBlackScholesEngine` `xmin_t` AD sensitivity.

Pass 5 identified **3 remaining / newly introduced items (`P5-01`, `P5-02`, `P5-03`)**. Resolve **all 3 items** below in full, add the specified regression tests, and ensure zero compiler/clippy warnings across the workspace.

---

## Mandatory Remediation Tasks (`P5-01` through `P5-03`)

### `[P5-01]` Fix `JuQuadraticEngine::ju_put_direct` Sign Errors in `theta` and `alpha_lambda_prime` & Add Non-Positive Rate Guard (`crates/ql-instruments/src/american.rs:443–494`)

Compare `ju_put_direct` (`crates/ql-instruments/src/american.rs:443–494`) against Ju (1999, Eqs. 14–18) and QuantLib `repos/QuantLib/ql/pricingengines/vanilla/juquadraticengine.cpp:86–87, 163–186` with $\phi = -1$:
1. **Non-Positive Rate Early-Exercise Guard (`line 444`)**:
   In QuantLib `juquadraticengine.cpp:86–87` (`if (riskFreeDiscount >= 1.0 && riskFreeDiscount >= dividendDiscount) return black.value();`), when $r \le 0$ and $r \le q$, an American put is never optimal to exercise early ($S^* = 0$):
   ```rust
   if t.val() <= 0.0 || v.val() <= 0.0 || (r.val() <= 0.0 && r.val() <= q.val()) {
       return european(s, k, r, q, v, t, false);
   }
   ```
2. **European Put Maturity Derivative `theta` ($\frac{\partial P_E}{\partial T}$ at $S^*$, `lines 481–483`)**:
   For a European put $P_E(S^*, T) = K e^{-rT} N(-d_2) - S^* e^{-qT} N(-d_1)$, differentiating with respect to $T$ (or setting `phi = -1` in `juquadraticengine.cpp:176–181`) gives $+q$ and $-r$:
   $$\frac{\partial P_E}{\partial T}(S^*, T) = +q S^* e^{-qT} N(-d_1(S^*)) - r K e^{-rT} N(-d_2(S^*)) + S^* e^{-qT}\phi(d_1(S^*))\frac{\sigma}{2\sqrt{T}}.$$
   Lines 481–482 copied the call signs (`-q * star * ... + r * k * ...`) instead of the put signs (`+q * star * ... - r * k * ...`)!
3. **Sign of `alpha_lambda_prime` for Put (`line 485`)**:
   For the put root $\lambda = \frac{1 - \beta - \text{disc}}{2}$, $\lambda'(h) = +\frac{\alpha}{h^2\,\text{disc}} > 0$, so in QuantLib (`juquadraticengine.cpp:163–164` with `phi = -1`):
   $$\alpha\lambda'(h) = -\phi \frac{(\alpha/h)^2}{\text{disc}} = +\frac{(\alpha/h)^2}{\text{disc}} > 0, \qquad \text{denom} = 2\lambda + \beta - 1 = -\text{disc} < 0.$$
   Setting both `alpha_lambda_prime = -alpha_over_h * alpha_over_h / disc` and `denom = -disc` cancelled the two minus signs and made `bju > 0` instead of `bju < 0`!
   Replace lines 481–490 with:
   ```rust
   let theta = q * star * (-q * t).exp() * (-d1s).norm_cdf()
       - r * k * (-r * t).exp() * (-d1s + v * t.sqrt()).norm_cdf()
       + star * (-q * t).exp() * pdf_s * v / (t.sqrt() * 2.0);
   let disc0 = (-r * t).exp();
   let alpha_lambda_prime = alpha_over_h * alpha_over_h / disc;
   let alpha_ve_h = T::from_f64(2.0) / (v2 * disc0) * theta;
   let denom = -disc;
   let bju = (T::one() - h) * alpha_lambda_prime / (denom * 2.0);
   let cju = -((T::one() - h) / denom)
       * (alpha_ve_h / h_a + alpha_over_h + alpha_lambda_prime / denom);
   ```

---

### `[P5-02]` Fix `Fd1dBlackScholesEngine` `UpAndOut` Put Boundary, `smax` & Knock-In Rebate (`crates/ql-instruments/src/pde1d.rs:62–72, 138–163, 254–302`) and `Fd2dHestonEngine` `nx` Alignment (`crates/ql-pde/src/adi.rs:196–250, 479–487`)

1. **`Fd1dBlackScholesEngine::fd` Boundary Conditions (`crates/ql-instruments/src/pde1d.rs:138–163, 179–180`)**:
   - On lines 145–147, `let putb = if let Some((b, _, rb)) = barrier { if smin_actual <= b.val() { rb } ...` ignores the barrier type `bt`! For **every `UpAndOut` put** ($S_{\min} < S_0 < B$), `smin_actual <= b.val()` is true, clobbering the lower boundary at $S = S_{\min}$ to `rb` (`0.0`) instead of $(K e^{-r\tau} - S_{\min} e^{-q\tau})^+$!
   - Check the barrier direction on both boundaries `j = 0` (`low_bc`) and `j = n` (`high_bc`), and support an optional `terminal_rebate: Option<T>` (analogous to `Fd2dHestonEngine::solve`):
     - At `j = 0` (`smin_actual`): if `let Some((b, BarrierType::DownAndOut, rb)) = barrier && smin_actual <= b.val()`, value is `rb`; else if `let Some(tr) = terminal_rebate`, value is `tr * (-rate * tnew).exp()`; else if `kind.is_call()`, `T::zero()`; else if `american`, `(strike - T::from_f64(smin_actual)).max_val(T::zero())`; else `(strike * (-rate * tnew).exp() - T::from_f64(smin_actual) * (-dividend * tnew).exp()).max_val(T::zero())`.
     - At `j = n` (`smax_actual`): if `let Some((b, BarrierType::UpAndOut, rb)) = barrier && smax_actual >= b.val()`, value is `rb`; else if `let Some(tr) = terminal_rebate`, value is `tr * (-rate * tnew).exp()`; else if `kind.is_call()`, `(T::from_f64(smax_actual) * (-dividend * tnew).exp() - strike * (-rate * tnew).exp()).max_val(T::zero())`; else `T::zero()`.
2. **`Fd1dBlackScholesEngine::fd` `smax` & `price_barrier_with_rebate` Knock-In Rebate (`crates/ql-instruments/src/pde1d.rs:62–72, 254–302`)**:
   - In `fd` (`lines 62–72`), pass `grid_anchor: Option<(T, BarrierType)>` so when `barrier` is `None` and `grid_anchor = Some((b, bt))`, `smax` uses the exact same `2.0 * b.val()` / `4.0 * b.val()` bound as `out`, keeping `dx` identical between `vanilla` and `out`.
   - In `price_barrier_with_rebate` (`lines 254–288`), for `UpAndIn | DownAndIn`:
     - Compute `out_no_rebate` with `T::zero()` rebate (`Some((barrier, out_bt, T::zero()))`), not `rebate`.
     - When `rebate.val().abs() > 0.0`, add `rebate_survival` (solved via `self.fd` with `barrier = Some((barrier, out_bt, T::zero()))` and `terminal_rebate = Some(rebate)` so terminal payoff is `rebate` on the un-knocked-out domain and `0` on the knocked-out domain), returning `vanilla - out_no_rebate + rebate_survival`.
3. **`Fd2dHestonEngine::solve` & `price_barrier` (`crates/ql-pde/src/adi.rs:196–250, 479–487`)**:
   - On line 479 of `crates/ql-pde/src/adi.rs`, change `let nx = if lower_barrier.is_some() || upper_barrier.is_some()` to `let nx = if lower_barrier.is_some() || upper_barrier.is_some() || grid_center.is_some()` so `vanilla` in `price_barrier` (`lines 210–224`) uses the exact same `(nx | 1)` node count (`81`) as `out_no_rebate`.
   - Move `let out = self.price_knockout_barrier(...)` (`lines 196–207`) inside the `UpAndOut | DownAndOut` match arm (`line 249`) so `UpAndIn | DownAndIn` does not execute an unused 4th 2D ADI solve.

---

### `[P5-03]` Fix `price_convertible_bond_tf` Clean vs. Dirty Call/Put Strike & Step Coupon Accumulation (`crates/ql-fixed-income/src/convertible_bond.rs:25–82`)

Compare `price_convertible_bond_tf` (`crates/ql-fixed-income/src/convertible_bond.rs:25–82`) against Tsiveriotis & Fernandes (1998) and QuantLib `DiscretizedConvertible::postAdjustValuesImpl()` (`repos/QuantLib/ql/pricingengines/bond/discretizedconvertible.cpp:151–163`, and `crates/ql-fixed-income/src/callable_bond.rs:173–193`):
1. **Clean (Ex-Coupon) Call/Put Strike Order (`lines 58–78`)**:
   - Currently, `coupon` is added to `val` and `bb` (`lines 59–60`) **before** applying `s.call_price` (`lines 61–66`) and `s.put_price` (`lines 67–72`). This compares the dirty continuation value against a clean call/put price and wipes out the coupon when called on a coupon date!
   - Apply `s.call_price` and `s.put_price` to the ex-coupon continuation value `(val, bb)` **first**, then add `coupon` (`val += coupon; bb += coupon;`), and finally apply the holder's conversion option (`if conv.val() > val.val() { val = conv; bb = T::zero(); }`):
     ```rust
     if let Some(c) = s.call_price
         && val.val() > c.val()
     {
         val = c;
         bb = c;
     }
     if let Some(q) = s.put_price
         && val.val() < q.val()
     {
         val = q;
         bb = q;
     }
     let coupon = step_coupons[k];
     val += coupon;
     bb += coupon;
     let spot = s.spot * u.powf((k - i) as f64) * d.powf(i as f64);
     let conv = s.conversion_ratio * spot;
     if conv.val() > val.val() {
         val = conv;
         bb = T::zero();
     }
     ```
2. **Pre-Accumulate Step Coupons Without `break` (`lines 25–39`)**:
   - Build `let mut step_coupons = vec![T::zero(); n + 1];` before the tree loop by iterating over `&s.coupon_times`:
     ```rust
     let mut step_coupons = vec![T::zero(); n + 1];
     let mut previous = T::zero();
     for &ct in &s.coupon_times {
         if ct.val() > previous.val() {
             let step_idx = ((ct / dt).val().round() as usize).clamp(1, n);
             step_coupons[step_idx] += s.face_value * s.coupon_rate * (ct - previous);
             previous = ct;
         }
     }
     let maturity_coupon = step_coupons[n];
     ```
     And at tree step `k` (`k in (0..n).rev()`), only add `step_coupons[k]` when `k > 0` (since coupons at $t = 0$ are ex-dividend).

---

## Verification & Test Checklist
1. **Unit & Integration Tests**:
   - Add regression tests verifying:
     - `JuQuadraticEngine::price` American Put prices (`is_call = false`) match high-resolution `Fd1dBlackScholesEngine` / `BaroneAdesiWhaleyEngine` across $r > q$, $r < q$, and $r = 0, q > 0$, and satisfy $P_{\text{Amer}}(r \le 0, q = 0) = P_{\text{Euro}}$.
     - `Fd1dBlackScholesEngine::price_barrier_option` for `UpAndOut` Put matches analytical barrier pricing, satisfies exact discrete In–Out parity (`UpAndIn + UpAndOut == vanilla` at `rebate == 0` for both 1D BS and 2D Heston ADI), and increases monotonically with `rebate > 0` for Knock-In options.
     - `price_convertible_bond_tf` with a clean `call_price` on a coupon date pays `call_price + coupon` (when not converted) rather than overwriting the coupon.
2. **Quality Gates**:
   - `cargo clippy --workspace --all-targets -- -D warnings` must pass with **0 warnings**.
   - `cargo test --workspace` must pass 100% of unit, integration, and parity tests.
