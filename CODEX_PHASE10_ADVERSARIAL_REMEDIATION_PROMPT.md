# Phase 10 Round 3 Adversarial Remediation Prompt for Codex

## Context & Audit Status
Following your Round 2 remediation (`19aa6ad..9a72cde`, 16 commits, 57 files changed), our 4-agent adversarial audit loop (`Big-Picture & Workflow Architect`, `Rust Systems & AD/SIMD Auditor`, `Phase 1–9 Mathematical Auditor`, and `Phase 10 & SOTA Literature Auditor`) conducted **Pass 3** across the entire 23-crate workspace against QuantLib C++ (`repos/QuantLib`) and primary literature (`papers/01`–`papers/45`).

Full audit reports are recorded in:
- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_SECOND_PASS.md`
- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_MATH.md`
- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_RUST_ARCH.md`
- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_ARCHITECTURE.md`

While Round 2 resolved the majority of Pass 2 items (`BaroneAdesiWhaleyEngine` direct put boundary $q_1 < 0$ & $F_{xx}$, `replicate_variance_swap_fair_strike` log-contract tails, `G2PlusPlus::swaption_price` 1D quadrature with AD Newton polish, `asb_recursive_expected_tranche_loss` 32-node Gauss–Hermite + high-correlation Legendre rule, `NelsonSiegelSvensson<T>` reference date, `AllocationError`, `PortfolioRiskSummary` per-call vs cumulative counters, and `simd_erfc` vectorized fdlibm approximation), **11 concrete mathematical, AD/IFT, SIMD, and architectural bugs remain or were newly introduced in Round 2 (`P3-01` through `P3-11`)**.

Execute **all 11 remediation items (`P3-01` through `P3-11`)** below completely, without shortcuts, placeholders, or tolerance loosening.

---

## Mandatory Remediation Checklist (`P3-01` – `P3-11`)

### `[P3-01]` CRITICAL — Fix `BjerksundStensland2002Engine`, `JuQuadraticEngine`, `BaroneAdesiWhaleyEngine::baw_put`, and `american_implied_vol` (`crates/ql-instruments/src/american.rs`)

1. **`BjerksundStensland2002Engine::bs_psi` and `bs2002_call` (`crates/ql-instruments/src/american.rs:287-346`)**:
   - **Bug A (Missing `* t` in `bs_psi` exponential prefactor, `lines 295, 308`)**:
     In `bs_psi`, `let lambda = -r + gamma * b + gamma * (gamma - T::one()) * v * v * 0.5;` is an annualized rate, and line 308 computes `lambda.exp()` **without multiplying by `t`** (compare `bs_phi` at lines 228–229: `let lambda = lambda_rate * t; lambda.exp()`).
   - **Bug B (Missing `s.pow(beta)` and `s` factors in `bs2002_call`, `lines 330-345`)**:
     Neither `bs_phi` (`lines 226-235`) nor `bs_psi` (`lines 287-309`) multiplies by $S^\gamma$ (notice `bs93_call` at `line 264` explicitly multiplies `phi_b` by `(s / i).pow(beta)` and `phi_1_*` by `s`). In `bs2002_call` (`lines 330-343`), `a1 = (i1 - k) / i1.pow(beta)` and `a2 = (i2 - k) / i2.pow(beta)` already divide by $I_1^\beta, I_2^\beta$, so all $\gamma = \beta$ terms in `bs2002_call` must be multiplied by `s.pow(beta)` and all $\gamma = 1$ terms must be multiplied by `s`!
   - **Bug C (Swapped $t_1/T$ barriers, missing $H$ in $f_3, f_4$, correlation sign, and reflection power in `bs_psi`, `lines 296-308`)**:
     In Bjerksund & Stensland (2002, Eqs. 6–14), $e_1..e_4$ are evaluated at $t_1$ (`st1 = v * t1.sqrt()`) against $I_1$ and $I_2^2 / I_1$, while $f_1..f_4$ are evaluated at $T$ (`st = v * t.sqrt()`) against $H$:
     ```rust
     let drift = b + (gamma - 0.5) * v * v;
     let e1 = -((s / i1).ln() + drift * t1) / st1;
     let e2 = -((i2 * i2 / (s * i1)).ln() + drift * t1) / st1;
     let e3 = -((s / i1).ln() - drift * t1) / st1;
     let e4 = -((i2 * i2 / (s * i1)).ln() - drift * t1) / st1;
     let f1 = -((s / h).ln() + drift * t) / st;
     let f2 = -((i2 * i2 / (s * h)).ln() + drift * t) / st;
     let f3 = -((i1 * i1 / (s * h)).ln() + drift * t) / st;
     let f4 = -((s * i1 * i1 / (h * i2 * i2)).ln() + drift * t) / st;
     let p11 = ql_math::bivariate_normal_cdf(e1, f1, rho);
     let p12 = ql_math::bivariate_normal_cdf(e2, f2, rho);
     let p13 = ql_math::bivariate_normal_cdf(e3, f3, -rho);
     let p14 = ql_math::bivariate_normal_cdf(e4, f4, -rho);
     (lambda * t).exp()
         * (p11 - (i2 / s).pow(kappa) * p12 - (i1 / s).pow(kappa) * p13
             + (i1 / i2).pow(kappa) * p14)
     ```
     And in `bs2002_call` (`lines 332-345`):
     ```rust
     if s.val() >= i2.val() {
         return (s - k).max_val(T::zero());
     }
     let phi = |gamma: T, h: T, i: T| Self::bs_phi(s, t1, gamma, h, i, r, b, v);
     let sb = s.pow(beta);
     let raw = a2 * sb * (T::one() - phi(beta, i2, i2))
         + s * (phi(T::one(), i2, i2) - phi(T::one(), i1, i2))
         - k * (phi(T::zero(), i2, i2) - phi(T::zero(), i1, i2))
         + a1 * sb * (phi(beta, i1, i2) - Self::bs_psi(s, t, beta, i2, i2, i1, t1, r, b, v))
         + s * (Self::bs_psi(s, t, T::one(), i2, i2, i1, t1, r, b, v)
             - Self::bs_psi(s, t, T::one(), k, i2, i1, t1, r, b, v))
         - k * (Self::bs_psi(s, t, T::zero(), i2, i2, i1, t1, r, b, v)
             - Self::bs_psi(s, t, T::zero(), k, i2, i1, t1, r, b, v));
     raw.max_val(euro).max_val((s - k).max_val(T::zero()))
     ```

2. **`JuQuadraticEngine::price` (`crates/ql-instruments/src/american.rs:386-451`)**:
   - **Bug A (Unlifted `d1s0` and `pdf` used in `d_ve_dt`, `lines 415-418, 431-433`)**:
     Lines 415–418 evaluate `d1s0` and `pdf` at the detached constant `star_t = T::from_f64(star0)` *before* `ift_lift_1d`, and then lines 431–433 reuse `d1s0` and `pdf` inside `d_ve_dt`, dropping $\frac{\partial d_1(S^*)}{\partial S^*} \frac{\partial S^*}{\partial \theta}$ from the AD graph (which forced Suite 19 in `crates/ql-benches/src/suites.rs:1495` to use `tol = 1e-5`!).
     **Fix**: After `let star = ift_lift_1d(...);` (`line 424`), recompute `d1s = ((star / k).ln() + (r - q + v2 * 0.5) * t) / (v * t.sqrt())` and `pdf_s = d1s.norm_pdf()` from the IFT-lifted `star`, and use `d1s` and `pdf_s` in `d_ve_dt`, `ve`, and `h_a`!
   - **Bug B (Zero-rate limit in `ju_call` breaks zero-dividend American puts, `lines 403-407, 434-442`)**:
     When pricing a zero-dividend American put (`div == 0.0`) via put-call duality on line 449 (`Self::ju_call(strike, spot, div, rate, vol, expiry)`), the dual call has `r_dual = div = 0.0`. As $r \to 0$, $\alpha = 2r/\sigma^2 \to 0$ and $h = 1 - e^{-rT} \to 0$, so `alpha_over_h` $= \frac{\alpha}{h} \to \frac{2}{\sigma^2 T} \neq 0$ (`line 395`). In Ju (1999) (see QuantLib `juquadraticengine.cpp:163-182`), the combinations appearing in $b_{\text{Ju}}$ and $c_{\text{Ju}}$ are:
     $$\alpha \lambda'(h) = -\frac{(\alpha / h)^2}{\text{disc}}, \qquad \alpha V_{E,h} = \frac{2}{\sigma^2 e^{-rT}} \frac{\partial V_E}{\partial T} + \left(\frac{\alpha}{h} - \frac{2}{\sigma^2 T}\right) \frac{V_E}{e^{-rT}}$$
     BOTH of which have finite, non-zero limits at $r = 0$! Currently lines 403–404 and 434–435 zero out `lambda_prime` and `ve_h` when $|r| < 10^{-10}$, forcing `bju = 0` and `cju = 0` for all zero-dividend American puts!
     **Fix**: Compute `alpha_lambda_prime = -alpha_over_h * alpha_over_h / disc;` and `alpha_ve_h = T::from_f64(2.0) / (v2 * (-r * t).exp()) * d_ve_dt + (alpha_over_h - T::from_f64(2.0) / (v2 * t)) * ve / (-r * t).exp();` directly (which is smooth and finite for all $r$, including $r = 0$), and use `alpha_lambda_prime` and `alpha_ve_h` in `bju = (T::one() - h) * alpha_lambda_prime / (T::from_f64(2.0) * (lambda - T::one()))` and `cju = -((T::one() - h) / (lambda - T::one())) * (alpha_ve_h / h_a + lambda_h + alpha_lambda_prime / (T::from_f64(2.0) * (lambda - T::one())))`! Then tighten Suite 19 (`crates/ql-benches/src/suites.rs:1495`) tolerance from `1e-5` to `1e-7`.

3. **`BaroneAdesiWhaleyEngine::baw_put` (`crates/ql-instruments/src/american.rs:172-174`)**:
   - The derivative of $K - x - P_E(x)$ wrt $x$ is $-1 - (-e^{-qT} N(-d_1)) = -(1 - d)$ where `d = (-div * expiry).exp() * (-d1).norm_cdf()`. Line 172 wrote `-(-div * expiry).exp() * d1.norm_cdf()`, which equals $-(1 - d)$ only when `div == 0`.
   - **Fix (`american.rs:172-174`)**:
     ```rust
     let fx = -(T::one() - T::one() / q1) * (T::one() - d)
         + (-div * expiry).exp() * pdf / (q1 * vol * expiry.sqrt());
     ```

4. **`american_implied_vol` (`crates/ql-instruments/src/american.rs:483-500`)**:
   - Keep `f_x` in `T` (preserving $\partial_\theta F_\sigma$ for 2nd-order Vanna/Volga sensitivities) instead of detaching `T::from_f64(jet.d1[0])`.

---

### `[P3-02]` CRITICAL — Fix `h1_hw_char_func_unit_spot`, Grzelak–Oosterlee `int1`, `CosHestonHullWhiteEngine::price`, and `go_sqrt_moment_params` (`crates/ql-models/src/hybrid.rs`)

1. **`h1_hw_char_func_unit_spot` (`crates/ql-models/src/hybrid.rs:120-122`)**:
   - Lines 120–122 still write `let iz = Complex::i() * uu; let q = iz * iz + iz; let fac = (-q * Complex::from_real(hyb * 0.5)).exp();`. Because `iz * iz = -uu * uu`, `q = -uu^2 + i * uu` and `-q = +uu^2 - i * uu`, so `fac` has a **positive** real exponent $+\frac{1}{2}\text{hyb}\,u^2$, exploding as $|u| \to \infty$!
   - **Fix (`hybrid.rs:120-122`)** (matching line 195):
     ```rust
     let q = Complex::new(uu * uu, uu);
     let fac = (-q * Complex::from_real(hyb * 0.5)).exp();
     ```

2. **Grzelak–Oosterlee Cross-Integral `int1` (`crates/ql-models/src/hybrid.rs:103-115, 174-186`)**:
   - The exact integral $\int_0^T e^{-c t}\frac{1 - e^{-\kappa_r(T-t)}}{\kappa_r}\,dt$ (QuantLib `analytich1hwengine.cpp:130-132`, Grzelak & Oosterlee 2011) is:
     $$\text{int1} = \frac{1}{\kappa_r}\left[\frac{1 - e^{-c T}}{c} - \frac{e^{-c T} - e^{-\kappa_r T}}{\kappa_r - c}\right]$$
     with smooth Taylor limits:
     - when $|\kappa_r T| < 10^{-6}$: $\text{int1} = \frac{c T - 1 + e^{-c T}}{c^2}$ (if $|c T| < 10^{-6}$, $\frac{T^2}{2} - \frac{c T^3}{6}$),
     - when $|(\kappa_r - c) T| < 10^{-6}$: replace $\frac{e^{-c T} - e^{-\kappa_r T}}{\kappa_r - c}$ with $T e^{-c T}\left(1 - \frac{1}{2}(\kappa_r - c)T + \frac{1}{6}((\kappa_r - c)T)^2\right)$.
   - Replace the broken `ksum = kr + c_go` formula at both `lines 103-114` and `lines 174-185` with a shared helper implementing this exact integral.

3. **`CosHestonHullWhiteEngine::price` (`crates/ql-models/src/hybrid.rs:159-213`)**:
   - At `line 187`, `extra` is only the pure Hull–White rate variance ($\sim 10^{-4}$), whereas the Heston equity variance ($\sim 0.04$) is already inside `phi_heston`. When $\rho_{sr} < 0$ and $|\text{cross}| > \text{extra}$ (any $\rho_{sr} < -0.02$), `let hyb = (extra + cross).max_val(T::zero());` **clamps `hyb` to `0.0`**, erasing the entire Hull–White rate variance and equity-rate correlation!
   - **Fix**: Compute `mean_var` (the Heston integrated variance) before `hyb`, clamp `let hyb = (extra + cross).max_val(-mean_var * 0.9);` so total log-spot variance remains positive without zeroing out negative-correlation hybrid terms, and pass `-hyb * 0.5` as `jump_drift` to `price_impl_extra` at `line 210`.

4. **`go_sqrt_moment_params` (`crates/ql-models/src/hybrid.rs:58-70`)**:
   - Replace `.max_val(T::from_f64(1e-14))` on `raw = theta - stationary_var / (theta * 8.0)` (`line 65`, which zeros AD derivatives when the Feller condition $\xi_v^2 < 16 \kappa_v \theta_v$ is violated) with a smooth positive floor `let a = (T::from_f64(0.5) * (raw + (raw * raw + T::from_f64(1e-12)).sqrt())).sqrt();`.

---

### `[P3-03]` CRITICAL — Fix `LognormalCmsSpreadPricer::lognormal_caplet_pv` and `convexity_factor` (`crates/ql-rates/src/cms_spread.rs`)

1. **`LognormalCmsSpreadPricer::lognormal_caplet_pv` (`crates/ql-rates/src/cms_spread.rs:114-131`)**:
   - Line 123 sets `e2 = (s.vol2 * sqt * (s.rho * z) - s.vol2 * s.vol2 * tf * 0.5).exp()`, completely omitting the independent normal shock $\sqrt{1 - \rho^2} Z_2$ of $S_2 \mid Z_1 = z$! Consequently, $\mathbb{E}[e_2] = \exp(-\frac{1}{2}(1 - \rho^2)\sigma_2^2 T_f) < 1$ whenever $|\rho| < 1$, violating $\mathbb{E}[S_2(T_f)] = m_2$, giving $S_2$ zero variance at $\rho = 0$, and causing `floorlet_pv` (`line 73`, computed via put-call parity from `Self::forward(s)`) to return **negative prices** for OTM floorlets!
   - **Fix**: Condition on $Z_2 = z$ so $S_2(z) = z_2 \exp(\sigma_2 \sqrt{T_f} z - \frac{1}{2}\sigma_2^2 T_f) - \text{shift}_2$ has exact mean $m_2$, and evaluate the inner expectation over $Z_1 \mid Z_2 = z$ analytically using the Black-76 call formula (matching Brigo & Mercurio Section 13.16.2 / QuantLib `lognormalcmsspreadpricer.cpp:83-109` and `basket_rainbow.rs:159-173`):
     - Clamp $\rho \in (-0.999999, 0.999999)$ and let $\sigma_{1,\text{cond}} = \sigma_1 \sqrt{(1 - \rho^2) T_f}$.
     - At each Gauss–Hermite node $z = \sqrt{2}\,X_i$:
       - $S_2(z) = z_2 \exp(\sigma_2 \sqrt{T_f} z - \frac{1}{2}\sigma_2^2 T_f) - \text{shift}_2$
       - Conditional forward of $g_1(S_1 + \text{shift}_1)$: $F_{1\mid z} = g_1 z_1 \exp(\rho \sigma_1 \sqrt{T_f} z - \frac{1}{2}\rho^2 \sigma_1^2 T_f)$
       - Effective strike: $K_{\text{eff}}(z) = \text{strike} - \text{spread} + g_1 \text{shift}_1 + g_2 S_2(z)$
       - If $K_{\text{eff}}(z) \le 0$: inner expectation is $F_{1\mid z} - K_{\text{eff}}(z)$.
       - Else if $\sigma_{1,\text{cond}} \le 10^{-12}$: inner expectation is $(F_{1\mid z} - K_{\text{eff}}(z))^+$.
       - Else: $d_1 = \frac{\ln(F_{1\mid z} / K_{\text{eff}}(z))}{\sigma_{1,\text{cond}}} + \frac{1}{2}\sigma_{1,\text{cond}}$, $d_2 = d_1 - \sigma_{1,\text{cond}}$, and inner expectation is $F_{1\mid z} N(d_1) - K_{\text{eff}}(z) N(d_2)$.

2. **`convexity_factor` (`crates/ql-rates/src/cms_spread.rs:29-40`)**:
   - Line 39 multiplies the Linear-TSR annuity convexity factor by `dt = (s.payment_time - s.fixing_time).max_val(T::zero())`, which forces the CMS convexity adjustment to **zero** whenever `payment_time == fixing_time` (including in `crates/ql-rates/tests/cms_spread_suite.rs`!).
   - **Fix**: Replace line 39 with `b * (annuity / (T::one() + annuity * rate.max_val(T::from_f64(1e-4)))) * (T::one() + dt / u)` so the TSR swap-annuity convexity adjustment is non-zero whenever $T_f > 0$ and $\sigma > 0$, even when $T_p = T_f$. Add tests in `cms_spread_suite.rs` verifying that `floorlet_pv > 0` for deep OTM strikes across $\rho \in \{-0.5, 0.0, 0.5\}$ and that `adjusted_rates` exceeds `(swap_rate1, swap_rate2)` when `payment_time == fixing_time`.

---

### `[P3-04]` CRITICAL — Fix `StulzRainbowEngine::expected_min` and `DengLiZhouSpreadBasketEngine` (`crates/ql-exotics/src/basket_rainbow.rs`)

1. **`StulzRainbowEngine::expected_min` (`crates/ql-exotics/src/basket_rainbow.rs:29-35, 63-71`)**:
   - Line 70 computes `s1 * (-q1 * t).exp() * d.norm_cdf() + s2 * (-q2 * t).exp() * (-d + st).norm_cdf()`, which is $\mathbb{E}[e^{-rT}\max(S_1, S_2)]$, **not** $\mathbb{E}[e^{-rT}\min(S_1, S_2)]$! Compare QuantLib `stulzengine.cpp:70-71`:
     $$\mathbb{E}[e^{-rT}\min(S_1(T), S_2(T))] = S_1 e^{-q_1 T} N(-d) + S_2 e^{-q_2 T} N(d - \sigma\sqrt{T})$$
   - Because `expected_min` returned $\mathbb{E}[\max]$, `emin` (`line 53`) and `emax` (`line 54`) were swapped, corrupting `PutOnMin`, `PutOnMax`, and $K = 0$ `CallOnMin`!
   - **Fix (`basket_rainbow.rs:29-35, 66-71`)**:
     In `expected_min`:
     ```rust
     if st.val() < 1e-8 {
         return (s1 * (-q1 * t).exp()).min_val(s2 * (-q2 * t).exp());
     }
     let d = ((s1 / s2).ln() + (q2 - q1 + var * 0.5) * t) / st;
     s1 * (-q1 * t).exp() * (-d).norm_cdf() + s2 * (-q2 * t).exp() * (d - st).norm_cdf()
     ```
     And for `k.val() <= 0.0` (`line 32`), `RainbowKind::CallOnMax` must return `s1 * (-q1 * t).exp() + s2 * (-q2 * t).exp() - Self::expected_min(s1, s2, q1, q2, v1, v2, rho, t)`.

2. **`DengLiZhouSpreadBasketEngine::price_spread` (`crates/ql-exotics/src/basket_rainbow.rs:163-173`)**:
   - When $K < 0$, `strike_cond = s2_t + k` can be $\le 0$, causing `(f1_cond / strike_cond).ln()` on line 168 to evaluate `ln(negative) = NaN`.
   - **Fix**: Guard `if strike_cond.val() <= 0.0 { sum += w * (f1_cond - strike_cond); continue; }` before taking `.ln()`.

---

### `[P3-05]` CRITICAL — Fix `CatBondSpec::price_mc` Poisson Arrival Boundary Bug, `expected_loss` Discounting, and `cat_bond_z_spread` (`crates/ql-fixed-income/src/cat_bond.rs`)

1. **`CatBondSpec::price_mc` (`crates/ql-fixed-income/src/cat_bond.rs:179-219`)**:
   - `let mut te = T::zero();` is initialized outside the coupon loop (`line 179`). Inside `for (i, &t) in self.payment_times.iter().enumerate()` (`line 181`), when an arrival time `te` exceeds `t` (`if te.val() > t { break; }` at `lines 193-195`), it breaks out of the inner `loop` *before* drawing the Pareto severity `u2`, and on the next coupon period `i + 1`, line 192 immediately adds *another* `-u.ln() / lambda` to `te`, **discarding the catastrophe that arrived at `te`**!
   - **Fix**: By the memoryless property of the homogeneous Poisson process, reset `te` at the start of each coupon period `i`:
     ```rust
     for (i, &t) in self.payment_times.iter().enumerate() {
         let prev_t = if i == 0 { 0.0 } else { self.payment_times[i - 1] };
         let mut te = T::from_f64(prev_t);
         loop {
             ...
     ```
2. **`expected_loss` Discounting Alignment (`crates/ql-fixed-income/src/cat_bond.rs:156` vs `253`)**:
   - In `price_analytical` (`line 156`), `expected_loss` sums incremental losses discounted at each coupon date `df(t_i)`, whereas in `price_mc` (`line 253`), `expected_loss` discounts cumulative loss at maturity `df(t_final)`. Align both so `price_mc` also accumulates `loss_pv += self.face_amount * (prev_surv - surv) * df;` at each coupon date `t_i` (matching `price_analytical`).
3. **`cat_bond_z_spread` (`crates/ql-fixed-income/src/cat_bond.rs:299-304`)**:
   - For `CatPerilSpec::PoissonPareto { lambda, xm, alpha, exhaustion }`, instead of hardcoding `surv = T::one()`, compute the effective survival factor at `t` (for exhaustion trigger: `let p_ex = if self.attachment.val() > xm.val() { (xm / self.attachment).pow(alpha) } else { T::one() }; (-lambda * p_ex * t).exp()`).

---

### `[P3-06]` CRITICAL — Fix `MarketContext` Shared `curve_store` Collision, `MarketContext::to_ad`, and `MarketSnapshot::to_market_context` (`crates/ql-curves/src/market_context.rs`, `crates/ql-portfolio/src/dag.rs`)

1. **Discount/Forward Curve Collision (`crates/ql-curves/src/market_context.rs:34, 93-106, 113-178`)**:
   - Both `insert_discount_curve` and `insert_forward_curve` write `format!("{id:?}")` (`"UsdSofr"`, `"GbpSonia"`) into the **same** `self.curve_store`, and `discount(name)` (`line 143`) checks `self.curve_store.get(name)` first! Thus, after `with_discount_curve(DiscountCurveId::USD_SOFR, ois).with_forward_curve(ForwardIndexId::USD_SOFR, fwd)`, calling `ctx.discount("UsdSofr")` returns the **forward curve `fwd`** instead of the OIS discount curve `ois`!
   - **Fix**: Never look up `self.curve_store` in `discount(name)` before `self.named_discount_curves` / `self.discount_curves`, and never look up `self.curve_store` in `forward(name)` before `self.named_forward_curves` / `self.forward_curves`! Specifically:
     - In `discount(name)`: check `self.named_discount_curves.get(name)` first, then match normalized `name` (`"USD_SOFR" | "USD-SOFR" | "USDSOFR"`) against `self.discount_curves`.
     - In `forward(name)`: check `self.named_forward_curves.get(name)` first, then match normalized `name` against `self.forward_curves`.
2. **`MarketContext::<f64>::to_ad` (`crates/ql-curves/src/market_context.rs:192-224`)**:
   - Currently `to_ad` uses `TapeVar::from_f64` (`idx: CONST_IDX`), which produces **identically zero AD gradients** when differentiated on a tape!
   - **Fix**: Have `to_ad` check if an active tape exists via `ql_ad::TapeScope`: if an active tape is recording, register each curve discount node, spot, and vol as a tape input via `ql_ad::TapeScope::with_active(|t| t.input(x))` (falling back to `TapeVar::from_f64` only if no tape is active).
3. **`MarketSnapshot::to_market_context` (`crates/ql-portfolio/src/dag.rs:164-190`)**:
   - Populate the OIS and forward curves in `to_market_context` from `try_build_curve_system(self)` (falling back to `(-r * t).exp()` only if bootstrap fails), and add unit tests in `crates/ql-portfolio/tests/portfolio_suite.rs` verifying `ctx.discount("UsdSofr")`, `ctx.forward("UsdSofr")`, and `ctx.to_ad()` on an active `TapeScope::guard()`.

---

### `[P3-07]` HIGH — Fix Knock-In Rebate Parity & Barrier Grid Node Alignment in `Fd2dHestonEngine` and `Fd1dBlackScholesEngine` (`crates/ql-pde/src/adi.rs`, `crates/ql-instruments/src/pde1d.rs`)

1. **`Fd2dHestonEngine::price_barrier` Knock-In Rebate (`crates/ql-pde/src/adi.rs:174-197`)**:
   - Line 196 computes `self.price_european(...) - out_no_rebate + (out - out_no_rebate)`.
   - `out - out_no_rebate` is the PV of receiving `rebate` **when the barrier IS hit** (a Knock-Out rebate)! A Knock-In option pays its rebate at expiry $T$ **if the barrier is NOT hit** (see QuantLib `fdhestonbarrierengine.cpp:196` and `fdhestonrebateengine.cpp:110`, where `rebateOption` has terminal payoff `rebate` inside the domain AND `rebate` on the barrier, so `rebateOption - out` cancels the barrier-hit rebate and leaves only the terminal rebate on paths that never hit the barrier!).
   - **Fix**: When `rebate != 0` for a Knock-In option (`UpAndIn` or `DownAndIn`), compute the Knock-In rebate by solving the knockout PDE with terminal condition `u(x, v, 0) = rebate` inside the survival domain and Dirichlet boundary `0` on the barrier (or equivalently `rebateOption(terminal=rebate, boundary=rebate) - out`), so Knock-In rebate is paid at expiry iff the barrier is not hit!
2. **`Fd2dHestonEngine::solve_knockout_barrier` Grid Alignment (`crates/ql-pde/src/adi.rs:420-425`)**:
   - When `self.nx` is even (default `nx = 80`), `dx = (x_max - x_min) / (nx - 1)` divides by odd `nx - 1 = 79`, so `focus = option.barrier.ln()` falls halfway between nodes `39` and `40`!
   - **Fix**: Shift `grid.x_min` and `grid.x_max` (just as `pde1d.rs:74-77` does) so `option.barrier.ln()` lands on an exact integer node `i_b = ((option.barrier.val().ln() - grid.x_min) / grid.dx).round()`.
3. **`Fd1dBlackScholesEngine` (`crates/ql-instruments/src/pde1d.rs:74-77, 126-139, 182-191, 205-225`)**:
   - Use the shifted grid endpoints `smin_actual = xmin.exp()` and `smax_actual = (xmin + n as f64 * dx).exp()` in `callb` and `putb` (`lines 126-137`), include `- T::from_f64(smin_actual) * (-dividend * tnew).exp()` in the European `putb`, and preserve `barrier_level`'s AD tangent in `xmin_t = barrier_level.ln() - T::from_f64(j * dx)` during final interpolation (`lines 182-191`) so `price_barrier` carries non-zero AD sensitivity $\partial V / \partial B$.
   - In `price_barrier` (`lines 205-225`), evaluate `vanilla` on the same barrier-shifted grid `xmin` as `out` so In-Out parity is exact on the discrete mesh.

---

### `[P3-08]` HIGH — Fix `CosHestonEngine` Cumulant Overflow (`0 * Inf = NaN`) & $c_4$ Truncation Width and `dupire_local_vol_ad` Tangents (`crates/ql-models/src/cos.rs`, `crates/ql-models/src/local_vol.rs`)

1. **`CosHestonEngine::price_impl_extra` (`crates/ql-models/src/cos.rs:190-217`)**:
   - Line 195 defines `j(a) = ((a * t).exp() - 1.0) / a` and lines 207, 209 evaluate `ek * j(kappa)` where `ek = (-kappa * t).exp()`. Algebraically $e^{-\kappa t}\frac{e^{\kappa t}-1}{\kappa} \equiv \frac{1 - e^{-\kappa t}}{\kappa} = j_k$, whereas computing `(kappa * t).exp()` first overflows to `+Inf` for large $\kappa t$, producing `0.0 * Inf = NaN`! Furthermore, `ek * (j(kappa) - jk)` in `a_int` (`line 207`) cancels the linear $t$ term of the stationary CIR variance!
   - **Fix**: Use the exact Fang–Oosterlee (2008, Table A.1) / QuantLib (`coshestonengine.cpp`) closed-form $c_2$ expressed purely in terms of $e^{-\kappa t}$ and $e^{-2\kappa t}$ (with smooth Taylor expansion when $|\kappa t| < 10^{-4}$):
     $$c_2(T) = \frac{1}{8\kappa^3}\Big[\sigma_v \kappa t e^{-\kappa t}(v_0 - \theta)(8\kappa\rho - 4\sigma_v) + \kappa\rho\sigma_v(1 - e^{-\kappa t})(16\theta - 8v_0) + 2\theta\kappa t(-4\kappa\rho\sigma_v + \sigma_v^2 + 4\kappa^2) + \sigma_v^2\big((\theta - 2v_0)e^{-2\kappa t} + \theta(6e^{-\kappa t} - 7) + 2v_0\big) + 8\kappa^2(v_0 - \theta)(1 - e^{-\kappa t})\Big]$$
     and include the fourth cumulant $c_4$ approximation in `width = l_trunc * (c2 + c4.max_val(T::zero()).sqrt()).sqrt()`.
2. **`dupire_local_vol_dual2_log_moneyness` (`crates/ql-models/src/local_vol.rs:14-31`)**:
   - Lines 17–19 construct `wt = Dual2::from_value(w.d1[0])` and `wy = Dual2::from_value(w.d1[1])`, which zeroes out their `.d1` tangents $\nabla w_t = (w_{tt}, w_{ty}) = (\text{w.d2}[0][0], \text{w.d2}[0][1])$ and $\nabla w_y = (w_{ty}, w_{yy}) = (\text{w.d2}[0][1], \text{w.d2}[1][1])$.
   - **Fix**: Populate `wt.d1 = [w.d2[0][0], w.d2[0][1]]` and `wy.d1 = [w.d2[0][1], w.d2[1][1]]` so first derivatives of `wt` and `wy` propagate through `(wt / denom).sqrt()`.

---

### `[P3-09]` HIGH — Fix `ift_lift_1d` Primal Residual Tangent Contamination, `bootstrap.rs` `BasisSwapInstrument` & `Dual2<64>` Guard, and `G2PlusPlus::swaption_price` Small-Rate Drift (`crates/ql-ad/src/ift.rs`, `crates/ql-curves/src/bootstrap.rs`, `crates/ql-rates/src/g2.rs`)

1. **`ift_lift_1d` (`crates/ql-ad/src/ift.rs:34-39`)**:
   - When `r0 = r_eval.val() != 0` (within solver tolerance) and $\partial_\theta f_x \neq 0$, `let d = -r_eval / f_x;` differentiates `1 / f_x` against the non-zero primal residual `r0`, adding a spurious $+r_0 (\partial_\theta f_x) / f_x^2$ error to the derivative.
   - **Fix (`ift.rs:34-39`)**:
     ```rust
     let x = T::from_f64(x_star);
     let r_eval = residual(x);
     let r0 = r_eval.val();
     let x_polished = T::from_f64(x_star - r0 / f_x.val());
     let d_tan = -(r_eval - T::from_f64(r0)) / f_x;
     x_polished + d_tan - T::from_f64(0.5 * f_xx) * d_tan * d_tan / f_x
     ```
   - Also pass non-zero `f_xx` to `ift_lift_1d` in `crates/ql-exotics/src/compound.rs:98` and `crates/ql-exotics/src/quanto_chooser.rs:125`.
2. **`CalibratedCurveSystem` (`crates/ql-curves/src/bootstrap.rs:385-394, 489-513, 557-725`)**:
   - In `from_multicurve_helpers` (`lines 385-394`), include `BasisSwapInstrument` (matching `from_single_curve_helpers` lines 311–312) instead of `else { continue; }`.
   - In `value_par_dv01_and_gamma` (`lines 489-513`), guard `if self.instruments.len() == k && k <= 64` (falling back to pairwise `Dual2<2>` via `ql_ad::curve_hessian_1d` if `k > 64`).
   - In `SingleCurveHelper` and `IrsSwapInstrument::residual_with_discount` (`lines 557, 574, 580-603, 615, 672, 680, 725`), use `curve_time(&self.day_counter, ...)` instead of `DayCounter::Act365Fixed`.
3. **`G2PlusPlus::swaption_price` (`crates/ql-rates/src/g2.rs:222-239`)**:
   - Replace `Self::nonzero(self.a)` in `mux` and `muy` with the exact cancellation-free identities:
     $$\frac{1 - e^{-aT}}{a^2} - \frac{1 - e^{-2aT}}{2a^2} = \frac{1}{2} B(a, T)^2, \qquad \frac{1 - e^{-aT}}{ab} - \frac{1 - e^{-(a+b)T}}{b(a+b)} = \frac{B(a, T) - B(a+b, T)}{b}$$
     using a smooth Taylor helper for $(B(a, T) - B(a+b, T))/b$ when $|b T| < 10^{-5}$: $\frac{T^2}{2} - \frac{(2a + b)T^3}{6} + \frac{(3a^2 + 3ab + b^2)T^4}{24}$.

---

### `[P3-10]` MEDIUM — Fix `simd_erfc` Extreme-Tail / `±Inf` Clamp and `DayCounter` `i64` Arithmetic (`crates/ql-ad/src/scalar.rs`, `crates/ql-time/src/day_count.rs`)

1. **`simd_erfc` (`crates/ql-ad/src/scalar.rs:345-351`)**:
   - Because `simd_exp` clamps its input at `-708.0` (`line 193`), for `ax >= 28.0` (`z^2 > 708`) `simd_exp(-z * z - 0.5625)` does not underflow to `0.0`, and for `x = ±Inf`, `z - ax = Inf - Inf = NaN`.
   - **Fix**: Mask `let r = ax.simd_ge(Simd::splat(28.0)).select(Simd::splat(0.0), r);` before `let res = ax.simd_lt(...)` at `line 348`, so `simd_erfc` and `SimdReal::norm_cdf` are exact across all $|x| \ge 28$ and $x = \pm\infty$. Add a unit test in `crates/ql-ad/tests/ad_suite.rs` testing `SimdReal::erfc` and `SimdReal::norm_cdf` at `±30.0`, `±40.0`, and `±f64::INFINITY`.
2. **`DayCounter` (`crates/ql-time/src/day_count.rs:25, 107, 118, 129, 134`)**:
   - Use `saturating_sub` in `DayCounter::day_count` (`line 25`) and `Sub for Date` (`line 134`), and compute `thirty_us`, `thirty_euro`, `thirty_isda` (`lines 92-130`) in `i64` so `(i64::from(b.year()) - i64::from(a.year())) * 360` never overflows `i32`.

---

### `[P3-11]` MEDIUM — Fix `DeskWorkstation` SIMM Inflation Bucket & Blotter Row Truncation, `PortfolioDag` Date/Pillar Consistency, and `CalibratedCreditCurve` Singular Jacobian (`crates/ql-terminal/src/lib.rs`, `crates/ql-portfolio/src/dag.rs`, `crates/ql-credit/src/lib.rs`)

1. **`DeskWorkstation::step_tick` and `render_dashboard` (`crates/ql-terminal/src/lib.rs:828-878, 926-933, 1405-1435`)**:
   - In `step_tick` (`lines 861-878`), `inflation` returned by `inflation_metrics` (`line 686`) is **already** scaled by `1e-4`, and `buckets[4]` (`line 861`) is **never read** in `class_margins` (`lines 870-877`)! Combine `inflation` (without multiplying by `1e-4` a second time) into the Interest Rate SIMM risk class (`class_margins[0]`, e.g., as a second bucket in `[buckets[0], inflation_bucket]` with inter-bucket correlation or inside `buckets[0]`).
   - In `step_tick` (`lines 926-933`), replace the `for m in 1..=5` integer loop with the `ResolvedVanillaSwap` valuation (matching `PortfolioDag::value_trade`).
   - In `render_dashboard` (`lines 1405-1435`), `right_bottom` prepends 8 telemetry lines before the 15 `blotter_rows` and truncates at 20 lines (`line 1432`), cutting off the last 3 blotter rows (`"SOFR RFR SABR 1Y"`, `"EURUSD XCCS 5Y"`, `"Portfolio ES 97.5%"`), while `left_bottom` has 9 blank lines (`lines 1405-1407`)! Move the SOTA telemetry lines into the bottom of `left_bottom` (or condense `right_bottom` so all 15 blotter rows AND the SOTA/American/ERC telemetry lines fit within the 20 rows without truncating any blotter row).
   - Add a unit test for `simulate_hw1f_multi_netting_set_cube` (`crates/ql-xva/src/cube.rs:114`).
2. **`PortfolioDag` (`crates/ql-portfolio/src/dag.rs:254, 264, 495, 553`)**:
   - In `try_build_curve_system` (`lines 254, 264`), use `reference_date.add_months(12 * years, false)` instead of `reference_date.add_days(365 * years)` so helper maturities match `value_trade` (`line 463`).
   - At lines 495 and 553, slice `&params[..curves.ois_pillars]` instead of `&params[..3]`.
3. **`CalibratedCreditCurve::bootstrap` (`crates/ql-credit/src/lib.rs:333-337`)**:
   - Replace `.unwrap_or_else(|| DMatrix::zeros(m, m))` with `.ok_or_else(|| CreditError::Solver("singular credit IFT Jacobian".into()))?`.

---

## Verification Gate
After completing all 11 items (`P3-01` through `P3-11`):
1. Run `cargo fmt --all -- --check`
2. Run `cargo clippy --workspace --all-targets -- -D warnings`
3. Run `cargo test --workspace`
4. Run `cargo run --release -p ql-benches` and refresh `crates/ql-benches/BENCHMARK_RESULTS.md` and `knowledge_base/08_BENCHMARK_RESULTS.md` if any benchmark checksums/tables change.
