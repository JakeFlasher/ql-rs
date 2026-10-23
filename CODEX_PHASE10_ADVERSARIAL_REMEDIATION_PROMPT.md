# CODEX PROMPT: Phase 10 & Workspace-Wide Adversarial Remediation — Round 2 (Post-`19aa6ad` Deep Audit)

You are operating on `/home/jakeshea/rust_quantlib` (a 23-crate Rust 2024 nightly quantitative finance workspace, HEAD at commit `19aa6ad`).
A 4-agent adversarial audit of your 9-commit remediation pass (`3559aac..19aa6ad`) verified that several Pass 1 items were genuinely resolved, **but uncovered critical mathematical shortcuts, fabricated perturbation formulas, newly introduced AD regressions, tolerance/checksum hacks in `ql-benches`, untouched files (`cms_spread.rs`, `adi.rs`, `cos.rs`, `basket_rainbow.rs`), and unwired dead code (`MarketContext`, `TapeGuard`, `compute_isda_simm_6_risk_classes`)**.

Before proceeding to the next phase, you must resolve **every single defect listed below** across all 4 workstreams with zero shortcuts, zero fake coefficients, and 100% mathematical, AD/IFT, and architectural fidelity.

Read the updated Pass 2 sections in the 4 adversarial review reports in `knowledge_base/` before writing code:
1. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_ARCHITECTURE.md` (Section: *Pass 2 Post-Remediation Audit*)
2. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_RUST_ARCH.md` (Section: *Pass 2 Post-Remediation Audit*)
3. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_MATH.md` (Section: *Pass 2 Post-Remediation Mathematical Audit*)
4. `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_SECOND_PASS.md` (Section: *Pass 2 Synthesis*)

---

## Strict Zero-Tolerance Anti-Shortcut Rules

1. **NO Fabricated Perturbations or Fake Coefficients**:
   - Never write `bs93_call(...) + 1e-6 * (short - european_short)` (`american.rs:200-204`), `bju = (q2 - 1.0)*(q2 - 2.0)*0.5; cju = q2 - 1.0;` (`american.rs:270-271`), `variance += mean_sq - variance` (`variance_swap.rs:59`), or `with_reference_date(self, _reference_date: Date) -> Self { self }` (`nss.rs:79-81`).
2. **NO Tolerance Scaling Hacks or Zero-Sum Checksums in `ql-benches`**:
   - Never divide benchmark checksums by `1_000.0` with `tol = 5e-4` to mask AD/FD mismatches (`suites.rs:1495, 1520, 1544`), and never sum unweighted normalized portfolio weight derivatives $\sum_k \frac{\partial w_k}{\partial \theta_\ell} \equiv 0$ (`suites.rs:1587-1651`).
3. **NO Finite Differences Inside Exact `Dual` / `Dual2` AD Kernels**:
   - Never replace exact Taylor-jet extraction with finite differences (`local_vol.rs:16-22`, `solvers.rs:143-149`).
4. **NO `.val()` Stripping in Generic `<T: Real>` Pricers**:
   - Every generic `<T: Real>` function (`JuQuadraticEngine::price`, `CatBondSpec::price_mc`, `cat_bond_z_spread`, `AnalyticCliquetEngine::price_mc_global`, `HullWhite1F`, `G2PlusPlus`, `rfr.rs`) must propagate exact AD derivatives (`Dual<N>`, `Dual2<N>`, `TapeVar<f64>`).
5. **NO Unwired Dead Code**:
   - Every public type/helper introduced to satisfy an architectural requirement (`MarketContext<T>`, `TapeGuard`, `compute_isda_simm_6_risk_classes`) must be wired into its production callers and verified by tests.

---

## Workstream 1: Phase 10 Mathematical SOTA & AD/IFT Fidelity (`MATH-P2-01` – `MATH-P2-04`)

### 1.1 `crates/ql-instruments/src/american.rs` — Fix Zero-Dividend Put Bug, Implement True Bjerksund–Stensland (2002), and Fix `JuQuadraticEngine` + `BaroneAdesiWhaleyEngine`
- **File**: `crates/ql-instruments/src/american.rs` (cross-reference `repos/QuantLib/ql/pricingengines/vanilla/baroneadesiwhaleyengine.cpp`, `bjerksundstenslandengine.cpp`, `juquadraticengine.cpp`)
1. **Fix Zero-Dividend American Put Early-Exercise Bug in `BS1993` & `BS2002` (`american.rs:149, 157, 179-188, 193, 211-220`)**:
   - When pricing a Put via McDonald–Schroder duality `bs93_call(k, s, q, r, v, t)` (where dual rate $r_{\text{dual}} = q$ and dual dividend $q_{\text{dual}} = r$), a non-dividend-paying American Put ($q = 0, r > 0$) has $r_{\text{dual}} = 0$ and $q_{\text{dual}} = r > 0$ (so dual cost-of-carry $b_{\text{dual}} = r_{\text{dual}} - q_{\text{dual}} = -r < 0$, meaning $b_{\text{dual}} < r_{\text{dual}}$!).
   - Currently, line 149 (`if q.val() <= 0.0 || t.val() <= 0.0 || r.val() <= 0.0`) and line 193 (`r.val().abs() < 1e-12`) immediately return the **European price** whenever $r_{\text{dual}} = q = 0$!
   - **Fix**: Only fall back to European when $b = r - q \ge r$ (i.e. `q.val() <= 0.0`) or `t.val() <= 0.0`! When $r = 0$ and $q > 0$ ($b = -q < 0 = r$), early exercise IS optimal:
     - Set $\beta = \left(\frac{1}{2} - \frac{b}{\sigma^2}\right) + \sqrt{\left(\frac{b}{\sigma^2} - \frac{1}{2}\right)^2 + \frac{2 r}{\sigma^2}}$ (which equals $1 - \frac{2b}{\sigma^2} = 1 + \frac{2q}{\sigma^2} > 1$ when $r = 0, b = -q < 0$).
     - Set $B_\infty = \frac{\beta}{\beta - 1} K$, and following QuantLib (`bjerksundstenslandengine.cpp:191-196`), when $r \le 0$ set $B_0 = K$, else $B_0 = \max\!\left(K, \frac{r}{q} K\right)$.
2. **Replace Fake `bs2002_call` (`american.rs:192-205`) with Genuine Two-Stage Bjerksund & Stensland (2002)**:
   - Delete the fake `base + T::from_f64(1e-6) * (short - european_short)` at lines 200–204.
   - Implement the genuine Bjerksund & Stensland (2002) two-stage flat boundary $(I_1, I_2)$ at $t_1 = \frac{1}{2}(\sqrt{5}-1)T$ and $T$:
     - With $b = r - q$, $\beta = \left(\frac{1}{2} - \frac{b}{\sigma^2}\right) + \sqrt{\left(\frac{b}{\sigma^2} - \frac{1}{2}\right)^2 + \frac{2 r}{\sigma^2}}$, $B_\infty = \frac{\beta}{\beta - 1}K$, $B_0 = \text{if } r \le 0 \{ K \} \text{ else } \{ \max(K, \frac{r}{q}K) \}$:
     - $h_1 = -(b t_1 + 2\sigma\sqrt{t_1})\frac{K^2}{(B_\infty - B_0)B_0}$, $I_1 = B_0 + (B_\infty - B_0)(1 - e^{h_1})$,
     - $h_2 = -(b T + 2\sigma\sqrt{T})\frac{K^2}{(B_\infty - B_0)B_0}$, $I_2 = B_0 + (B_\infty - B_0)(1 - e^{h_2})$,
     - $\alpha_1 = (I_1 - K) I_1^{-\beta}$, $\alpha_2 = (I_2 - K) I_2^{-\beta}$,
     - If $S \ge I_2$, return $S - K$; otherwise evaluate:
       $$C_{\text{BS02}}(S) = \alpha_2 S^\beta - \alpha_2 \phi(S, t_1, \beta, I_2, I_2) + \phi(S, t_1, 1, I_2, I_2) - \phi(S, t_1, 1, I_1, I_2) - K\phi(S, t_1, 0, I_2, I_2) + K\phi(S, t_1, 0, I_1, I_2) + \alpha_1 \phi(S, t_1, \beta, I_1, I_2) - \alpha_1 \psi(S, T, \beta, I_2, I_2, I_1, t_1) + \psi(S, T, 1, I_2, I_2, I_1, t_1) - \psi(S, T, 1, K, I_2, I_1, t_1) - K\psi(S, T, 0, I_2, I_2, I_1, t_1) + K\psi(S, T, 0, K, I_2, I_1, t_1)$$
       where $\psi(S, T, \gamma, H, I_2, I_1, t_1)$ uses `ql_math::bivariate_normal_cdf` with correlation $\rho = \sqrt{t_1 / T}$ on the four $(e_i, f_i)$ limits exactly as in Bjerksund & Stensland (2002, Eqs. 11–18).
3. **Fix `JuQuadraticEngine::price` (`american.rs:245-289`) — True Ju (1999) Coefficients AND `ift_lift_1d` on $S^*$**:
   - Delete the fabricated `let bju = (q2 - 1.0) * (q2 - 2.0) * 0.5; let cju = q2 - 1.0;` (`lines 270-271`) and the unlifted `let star = T::from_f64(root.max(...));` (`line 260`).
   - Lift `star` via `ift_lift_1d` with exact first and second derivatives `f_x` (in `T`) and `f_xx` (in `f64`) of the critical boundary residual!
   - Compute the true Ju (1999, Eqs. 12–16; `repos/QuantLib/ql/pricingengines/vanilla/juquadraticengine.cpp:159-186`) coefficients in `T`:
     - Let $\alpha = \frac{2r}{\sigma^2}$, $\beta = \frac{2(r - q)}{\sigma^2}$, $h = 1 - e^{-rT}$, $\text{disc} = \sqrt{(\beta - 1)^2 + \frac{4\alpha}{h}}$, $\lambda = \frac{-(\beta - 1) + \text{disc}}{2} = q_2$, and $\lambda' = \frac{\partial \lambda}{\partial h} = -\frac{\alpha}{h^2\,\text{disc}}$.
     - At $S^*$, let $V_E(S^*)$ be the European call price, $h_A = S^* - K - V_E(S^*)$, and let $V_{E,h}(S^*) = \frac{1}{r e^{-rT}}\frac{\partial V_E(S^*, T)}{\partial T}$ where $\frac{\partial V_E(S^*, T)}{\partial T} = -q S^* e^{-qT}\Phi(d_1(S^*)) + r K e^{-rT}\Phi(d_2(S^*)) + \frac{S^* e^{-qT}\phi(d_1(S^*))\sigma}{2\sqrt{T}}$.
     - Then:
       $$b_{\text{Ju}} = \frac{(1 - h)\alpha \lambda'}{2(2\lambda + \beta - 1)}, \qquad c_{\text{Ju}} = -\frac{1 - h}{2\lambda + \beta - 1}\left(\frac{\alpha V_{E,h}(S^*)}{h_A} + \frac{\alpha}{h} + \frac{\alpha \lambda'}{2\lambda + \beta - 1}\right)$$
       $$\chi(S) = b_{\text{Ju}}\left(\ln\frac{S}{S^*}\right)^2 + c_{\text{Ju}}\ln\frac{S}{S^*}, \qquad V_{\text{Ju}}(S) = V_E(S) + \frac{h_A (S / S^*)^\lambda}{1 - \chi(S)}.$$
     - Notice that at $S = S^*$, $\chi(S^*) = 0$ and $\frac{\partial V_{\text{Ju}}}{\partial S^*}\big|_{S < S^*} \neq 0$, so lifting $S^*$ with `ift_lift_1d` makes `Dual<5>` and `Dual2<2>` sensitivities match central finite differences to $< 10^{-8}$!
4. **Fix `BaroneAdesiWhaleyEngine` Put Critical Boundary & Second Derivative `f_xx` (`american.rs:78-134`)**:
   - In `baw_call` (`line 92`), pass the analytical second derivative `f_xx` of the BAW residual wrt $S^*$ into `ift_lift_1d` so `Dual2` Gamma/Vanna/Volga are exact.
   - Implement direct BAW put critical boundary $S^{**} < K$ using the negative root $q_1 = \frac{1}{2}\left(-(n - 1) - \sqrt{(n - 1)^2 + \frac{4m}{1 - e^{-rT}}}\right) < 0$ and $K - S^{**} - P_E(S^{**}) = -\frac{S^{**}}{q_1}(1 - e^{-qT}\Phi(-d_1(S^{**})))$ so `BaroneAdesiWhaleyEngine::critical_spot(..., false)` returns the actual put exercise boundary $S^{**} < K$.
   - In `american_implied_vol` (`american.rs:312-329`), compute `vega` in `T` (by evaluating `price` at `vol = T::from_f64(root)` and differentiating wrt `vol` while retaining `T` tangents of `(spot, strike, rate, div, expiry)`) so `Dual2` mixed second derivatives of `american_implied_vol` are exact.

---

### 1.2 `crates/ql-models/src/hybrid.rs` — Fix Exponential Blowup Sign Error, Double-$i$ Bug, $\kappa_v \to 0$ Limit, and Grzelak–Oosterlee $I_4(u, T)$
- **File**: `crates/ql-models/src/hybrid.rs` (cross-reference `repos/QuantLib/ql/pricingengines/vanilla/analytich1hwengine.cpp`)
1. **Fix Double-$i$ and Sign Error in `h1_hw_char_func_unit_spot` (`hybrid.rs:72, 101-105`) and `CosHestonHullWhiteEngine::price` (`hybrid.rs:140-185`)**:
   - `HestonParams::char_func_unit_spot(u: Complex<T>, ...)` (`heston.rs:22-57`) takes a frequency $u$ (where `u` is real on the Fourier axis) and multiplies by $i$ internally (`i * u`).
   - In `h1_hw_char_func_unit_spot` (`line 72`), pass `Complex::from_real(T::from_f64(u))` (NOT `Complex::new(0, u)`!).
   - For real frequency $u$, the Gaussian hybrid variance/martingale factor $\exp\!\left(-\frac{1}{2}\Sigma_{\text{hyb}}^2 (u^2 + i u)\right)$ has exponent $-\frac{1}{2}\Sigma_{\text{hyb}}^2 u^2 - \frac{1}{2}i u \Sigma_{\text{hyb}}^2$ (whose real part is **NEGATIVE** $-\frac{1}{2}\Sigma_{\text{hyb}}^2 u^2 \le 0$ and satisfies $\phi(-i) = 1$!).
   - In `CosHestonHullWhiteEngine::price` (`line 152-185`), do NOT duplicate the characteristic function inside the closure! Precompute the time-dependent hybrid variance `hyb` **once** before the COS loop, pass `hyb` (NOT `T::zero()`!) as `jump_variance` to `price_impl_extra` (`line 183`) so the COS interval $[a, b]$ widens properly, and inside the closure evaluate `h.char_func_unit_spot(u, rate, div, expiry) * (-(u * u + Complex::i() * u) * Complex::from_real(hyb * 0.5)).exp()`.
2. **Fix $\kappa_v \to 0$ Limit and Grzelak–Oosterlee (2011) $\mathbb{E}[\sqrt{v_t}]$ Moment-Matching (`hybrid.rs:77-99`)**:
   - When `|kappa_v| < 1e-8`, $\bar{v}(t) = \theta_v + (v_0 - \theta_v)\frac{1 - e^{-\kappa_v t}}{\kappa_v t} \to v_0 - \frac{1}{2}(v_0 - \theta_v)\kappa_v t + \frac{1}{6}(v_0 - \theta_v)\kappa_v^2 t^2$ (NOT `p.theta_v`!).
   - Implement Grzelak & Oosterlee (2011, Lemma 3.6; `analytich1hwengine.cpp:68-135`) moment-matched exponential approximation $\mathbb{E}[\sqrt{v_t}] \approx a_{\text{GO}} + b_{\text{GO}} e^{-c_{\text{GO}} t}$ for the stock-rate covariance integral $\int_0^T \rho_{sr}\eta_r \mathbb{E}[\sqrt{v_t}] B_r(t, T)\,dt$.

---

### 1.3 `crates/ql-exotics/src/basket_rainbow.rs`, `cliquet_extensible.rs`, & `variance_swap.rs` — Fix Stulz (1982), Li–Deng–Zhou (2008), Levy (1992), Cliquet MC AD, and Variance Swap Self-Cancellation
1. **`StulzRainbowEngine::price` (`crates/ql-exotics/src/basket_rainbow.rs:20-72`)**:
   - `cmin` at lines 44–46 is still completely wrong (`CallOnMin -> S1*e^{-q1 T} + S2*e^{-q2 T}` as $K \to 0$!).
   - Implement the exact Stulz (1982, Eqs. 6–12; `repos/QuantLib/ql/pricingengines/basket/stulzengine.cpp:46-134`) formula:
     - Let $\sigma^2 = \sigma_1^2 + \sigma_2^2 - 2\rho\sigma_1\sigma_2$, $\sigma\sqrt{T} = \sqrt{\max(\sigma^2 T, 10^{-16})}$,
     - $d = \frac{\ln(S_1/S_2) + (q_2 - q_1 + \frac{1}{2}\sigma^2)T}{\sigma\sqrt{T}}$,
     - $y_1 = \frac{\ln(S_1/K) + (r - q_1 + \frac{1}{2}\sigma_1^2)T}{\sigma_1\sqrt{T}}$, $y_2 = \frac{\ln(S_2/K) + (r - q_2 + \frac{1}{2}\sigma_2^2)T}{\sigma_2\sqrt{T}}$,
     - $\rho_1 = \frac{\sigma_1 - \rho\sigma_2}{\sigma}$, $\rho_2 = \frac{\sigma_2 - \rho\sigma_1}{\sigma}$,
     - Then the European **Call on Minimum** is:
       $$C_{\min}(S_1, S_2, K) = S_1 e^{-q_1 T} \Phi_2(y_1, -d; -\rho_1) + S_2 e^{-q_2 T} \Phi_2(y_2, d - \sigma\sqrt{T}; -\rho_2) - K e^{-rT} \Phi_2(y_1 - \sigma_1\sqrt{T}, y_2 - \sigma_2\sqrt{T}; \rho)$$
     - And by Stulz (1982) parity:
       - $C_{\max}(S_1, S_2, K) = C_{\text{BS}}(S_1, K) + C_{\text{BS}}(S_2, K) - C_{\min}(S_1, S_2, K)$,
       - $P_{\min}(S_1, S_2, K) = K e^{-rT} - C_{\min}(S_1, S_2, 0) + C_{\min}(S_1, S_2, K)$ where $C_{\min}(S_1, S_2, 0) = S_1 e^{-q_1 T}\Phi(-d) + S_2 e^{-q_2 T}\Phi(d - \sigma\sqrt{T})$ is Margrabe's minimum-of-two-assets forward value,
       - $P_{\max}(S_1, S_2, K) = K e^{-rT} - C_{\max}(S_1, S_2, 0) + C_{\max}(S_1, S_2, K)$ where $C_{\max}(S_1, S_2, 0) = S_1 e^{-q_1 T} + S_2 e^{-q_2 T} - C_{\min}(S_1, S_2, 0)$.
2. **`DengLiZhouSpreadBasketEngine` (`crates/ql-exotics/src/basket_rainbow.rs:104-145`)**:
   - Currently, `price_spread` (`lines 132-144`) computes `f1 = s1 * ((r - q1)*t).exp()` and passes `f1` as `spot` into `black_scholes_price(f1, ..., r, 0.0)` (**double-compounding `f1` by `e^{rT}`!**) and uses `+ 2.0 * rho * v1 * v2` instead of `- 2.0 * rho * v1 * v2` for a spread option!
   - Implement the genuine Li, Deng & Zhou (2008) second-order log-cumulant boundary expansion $J_0 + J_1 + J_2$ (`repos/QuantLib/ql/pricingengines/basket/denglizhoubasketengine.cpp`) for both 2-asset spread options $(S_1(T) - S_2(T) - K)^+$ and $N$-asset basket options,discounted once by $e^{-rT}$ on forward prices $F_i = S_i e^{(r - q_i)T}$.
3. **`LevyBasketEngine::price_basket` (`crates/ql-exotics/src/basket_rainbow.rs:167-180`)**:
   - Replace `let vol = (var.sqrt() / forward) / expiry.sqrt();` (`lines 178-179`) with Levy's (1992) exact log-moment formula:
     $$\sigma_{\text{Levy}} = \sqrt{\frac{\ln(M_2 / M_1^2)}{T}} = \sqrt{\frac{\ln(1 + \text{var} / (F^2))}{T}}.$$
4. **`AnalyticCliquetEngine::price_mc_global` (`crates/ql-exotics/src/cliquet_extensible.rs:102-135`)**:
   - Do NOT strip inputs via `.val()`! Run the pathwise simulation directly in `T: Real` (using `s.spot`, `s.rate`, `s.div`, `s.vol`, `s.local_floor`, `s.local_cap`, `s.global_floor`, `s.global_cap`, and `max_val` / `min_val`) so `price_mc_global` preserves exact pathwise AD sensitivities (`Dual`, `Dual2`, `TapeVar`).
5. **`replicate_variance_swap_fair_strike` (`crates/ql-exotics/src/variance_swap.rs:21-62`)**:
   - Delete the `if flat_surface { ... variance += mean_sq - variance; }` block (`lines 21-23, 53-60`), which algebraically cancels `variance`!
   - Use proper trapezoidal weights $\Delta K_0 = \frac{1}{2}(K_1 - K_0)$, $\Delta K_i = \frac{1}{2}(K_{i+1} - K_{i-1})$, $\Delta K_{n-1} = \frac{1}{2}(K_{n-1} - K_{n-2})$ on $[K_0, K_{n-1}]$, and add the exact analytical Black–Scholes log-contract left put tail on $(0, K_0)$ (with $\sigma_0 = \text{quotes}[0].\text{implied\_vol}$, $v_0 = \sigma_0\sqrt{T}$) and right call tail on $(K_{n-1}, \infty)$ (with $\sigma_{n-1} = \text{quotes}[n-1].\text{implied\_vol}$, $v_m = \sigma_{n-1}\sqrt{T}$):
     $$I_{\text{put}}(K_0, \sigma_0) = \frac{F}{K_0}\Phi(-d_1(K_0)) - \Phi(-d_2(K_0)) + v_0\left(-d_2(K_0)\Phi(-d_2(K_0)) + \phi(d_2(K_0))\right),$$
     $$I_{\text{call}}(K_{n-1}, \sigma_{n-1}) = \frac{F}{K_{n-1}}\Phi(d_1(K_{n-1})) - \Phi(d_2(K_{n-1})) - v_m\left(d_2(K_{n-1})\Phi(d_2(K_{n-1})) + \phi(d_2(K_{n-1}))\right).$$
     Then $\text{variance} = \frac{2}{T}\left(\text{drift} + I_{\text{put}}(K_0, \sigma_0) + e^{rT}\sum_{i=0}^{n-1} \frac{\Delta K_i}{K_i^2} Q(K_i) + I_{\text{call}}(K_{n-1}, \sigma_{n-1})\right)$.

---

### 1.4 `crates/ql-rates/src/cms_spread.rs`, `crates/ql-fixed-income/src/cat_bond.rs`, & `crates/ql-portfolio/src/allocation.rs`
1. **`LognormalCmsSpreadPricer` (`crates/ql-rates/src/cms_spread.rs:1-83`) — Must Be Implemented (Was Untouched in Round 1!)**:
   - Implement Linear-TSR / Hagan CMS convexity adjustments using `mean_reversion`, `swap_tenor1`, `swap_tenor2`, `annuity1`, `annuity2`, `fixing_time`, and `payment_time`:
     $$\mu_i = S_i + \frac{1}{2} S_i^2 \sigma_i^2 T_f \, C_i(a, \tau_i, T_p - T_f)$$
     so every field of `CmsSpreadSpec<T>` is genuinely used.
   - Fix `collar_pv` (`line 69`) to return `Self::caplet_pv(s, cap_strike) - Self::floorlet_pv(s, floor_strike)`.
   - Replace `lognormal_caplet_pv` (`lines 78-82`) with Brigo–Mercurio / QuantLib (`repos/QuantLib/ql/experimental/coupons/lognormalcmsspreadpricer.cpp`) **16-node Gauss–Hermite quadrature** over the conditional shifted-lognormal law of $(S_1 + \text{shift}_1, S_2 + \text{shift}_2)$ in `T: Real`.
2. **`CatBondSpec` & `cat_bond_z_spread` (`crates/ql-fixed-income/src/cat_bond.rs:94-296`)**:
   - In `cat_bond_z_spread` (`lines 270-295`), `fx` and `fxx` MUST include the catastrophe survival factor `surv_i` at each payment time $t_i$ (matching `price_analytical()`), NOT loss-free `cf`! Compute the exact `fx` (in `T`) and `fxx` (in `f64`) by differentiating the survival-weighted cashflows $\text{CF}_i^{\text{surv}} = \text{face\_amount} \cdot (\text{coupon\_rate} \cdot a_i \cdot \text{surv}_i + \text{red}_i \cdot \text{surv}_i)$ with respect to $z$.
   - In `price_analytical` (`lines 151-161`) and `price_mc` (`lines 236-249`):
     - Do NOT sum the cumulative principal loss `face_amount * (1 - surv_i) * df_i` across every intermediate coupon date! Principal loss is incurred once (or incrementally $\Delta\text{surv}_i = \text{surv}_{i-1} - \text{surv}_i$ at period $i$, or at maturity $T_{\text{final}}$: `expected_loss = self.face_amount * (T::one() - surv_final) * df_final` or $\sum_i \text{face\_amount} \cdot (\text{surv}_{i-1} - \text{surv}_i) \cdot \text{df}_i$).
     - Set `attachment_prob` and `exhaustion_prob` to the horizon probabilities at final maturity $T_{\text{final}}$ (not the unweighted average across intermediate coupon dates).
   - In `price_mc` (`lines 163-249`), keep `face_amount`, `coupon_rate`, `discount_rate`, and pathwise severity/proportional-exhaustion arithmetic in `T: Real` so `price_mc` and `PoissonPareto` carry non-zero AD sensitivities!
3. **`black_litterman_posterior` & `risk_budgeting_weights` (`crates/ql-portfolio/src/allocation.rs:137-342`)**:
   - Delete `solve_f64` (`lines 51-79`) completely and use `solve_linear_generic` everywhere.
   - Do not swallow `AllocationError` via `.unwrap_or_else(...)` (`lines 174, 187, 219, 234, 242`); provide `try_black_litterman_posterior`, `try_markowitz_minimum_variance_weights`, `try_markowitz_tangency_weights`, and `try_risk_budgeting_weights` returning `Result<..., AllocationError>` (with `.expect(...)` in the infallible wrappers).
   - In `black_litterman_posterior` (`lines 215-229`), ensure consistency with `implied_equilibrium_returns` ($\Pi = \delta \Sigma w_{\text{mkt}}$): solve $(\delta \Sigma) w_{\text{raw}} = \mu_{\text{BL}}$ so that when investor views equal equilibrium ($Q = P\Pi \implies \mu_{\text{BL}} = \Pi$) or $\tau \to 0$, $w$ recovers $w_{\text{mkt}}$!
   - Expand `crates/ql-portfolio/tests/allocation_suite.rs` to test `black_litterman_posterior` equilibrium recovery, `markowitz_minimum_variance_weights`, `markowitz_tangency_weights`, and `Dual` Jacobian parity against central finite differences.

---

## Workstream 2: Phase 1–9 Mathematical & Numerical Remediation (`MATH-P2-05`)

1. **`crates/ql-models/src/cos.rs:160-186` (`CosHestonEngine::price_impl_extra`)**:
   - Pass `rho: T` into `price_impl` / `price_impl_extra`, use `sigma_v` and `rho` to compute the Fang & Oosterlee (2008, Table A.1) second cumulant $c_2^{\text{Heston}}$ (with $|\kappa T| < 10^{-4}$ Taylor branch), and for Bates (`line 113-126`) add the jump drift $+\lambda(\mu_J - ((e^{\mu_J + \frac{1}{2}\sigma_J^2} - 1)))T$ to $c_1$ rather than subtracting `0.5 * jump_variance` (`line 182`).
2. **`crates/ql-pde/src/adi.rs:144-187, 254-368, 385, 501-512` & `crates/ql-instruments/src/pde1d.rs:58-71, 276-287`**:
   - In `Fd2dHestonEngine::price_barrier` (`adi.rs:167`), pass `T::from_f64(option.rebate)` instead of `T::zero()`, and for Knock-In with non-zero rebate subtract the knock-out rebate value correctly.
   - In `Fd2dHestonEngine::solve` (`adi.rs:254-368`), restrict `sweep_x` and `sweep_v` to active interior nodes inside `(lower_barrier, upper_barrier)` with Dirichlet boundary value `rebate` on the barrier rows/columns.
   - In both `Fd1dBlackScholesEngine` (`pde1d.rs:58-71`) and `Fd2dHestonEngine` (`adi.rs:385`), shift/stretch the spatial log-spot grid in `T: Real` so a grid node aligns exactly with `barrier.ln()` (`A12`), eliminating $O(\Delta x)$ barrier quantization error and preserving non-zero AD barrier sensitivity $\partial V / \partial B$.
3. **`crates/ql-models/src/local_vol.rs:8-30` (`DupireLocalVol::dupire_local_vol_ad`)**:
   - Restore exact `Dual2<2>` Taylor-jet derivative extraction (`wt = w.d1[0]`, `wy = w.d1[1]`, `wyy = w.d2[1][1]`) when `expiry` is seeded on lane 0 and `strike` is seeded on lane 1 (matching the doc comment *"no finite differences are taken"*), and if incoming tangents are present on `strike`/`expiry`, propagate them analytically without dividing `Dual2.d2` differences by $h^2 = 10^{-8}$.
4. **`crates/ql-models/src/hull_white.rs:24-27, 45-48`**:
   - Fix the Taylor expansion of $\frac{1 - e^{-2at}}{2a}$: change `4.0 / 3.0` to `2.0 / 3.0` at both lines 26 and 47.
5. **`crates/ql-models/src/rough.rs:18-20, 68-77` & `crates/ql-mc/src/paths.rs:303-319`**:
   - In `rough.rs:68-77`, only short-circuit to classical Heston when `(p.hurst - 0.5).abs() < 1e-8`, NOT whenever `p.n_factors == 1` with $H < 0.5$.
   - In `paths.rs:303-319`, initialize `let mut k0 = -rho * kappa * theta * dt / sigma_v;` so when the QE-M regularity condition fails (`den <= 1e-12` or `beta <= a_corr + 1e-12`), `qe_step` falls back to regular Andersen QE drift instead of `k0 = 0.0`.
6. **`crates/ql-fixed-income/src/callable_bond.rs:150-156` & `crates/ql-credit/src/cdo.rs:184-191`**:
   - In `tree_callable_bond_price` (`callable_bond.rs:150-156`), use a single calibrated node short rate $r_{j,i} = r + \text{oas} + \theta_j + \sigma\sqrt{\Delta t}(2i - j)$ at node $(j, i)$ (where $\theta_j$ fits the flat/initial discount factor $e^{-(r+\text{oas})t_j}$).
   - In `asb_expected_tranche_loss_quadrature` (`cdo.rs:184-191`), remove the nonlinear `total / quadrature_total` rescaling before clamping.

---

## Workstream 3: Big-Picture Architecture, `MarketContext`, `ql-xva`, `ql-terminal`, & `ql-benches` (`ARCH-P2-01` – `ARCH-P2-05`)

1. **`NelsonSiegelSvensson<T>` (`crates/ql-fixed-income/src/nss.rs:5-110`)**:
   - Add `pub reference_date: Date` to `NelsonSiegelSvensson<T>`, make `with_reference_date(mut self, reference_date: Date) -> Self` actually set `self.reference_date = reference_date`, and use `self.reference_date` in `YieldTermStructure::reference_date` and `forward_rate` (replacing hardcoded `Date(0)`).
2. **`MarketContext<T>` (`crates/ql-curves/src/market_context.rs:27-140`) & `PortfolioDag` (`crates/ql-portfolio/src/dag.rs`)**:
   - Unify `MarketContext<T>` curve maps so enum IDs (`DiscountCurveId`, `ForwardCurveId`) and string keys share the same backing store, add `MarketContext<f64>::to_ad`, and add `MarketSnapshot::to_market_context(&self) -> MarketContext<f64>` so `MarketContext` is actively exercised by `ql-portfolio` and tests.
   - In `PortfolioDag` (`crates/ql-portfolio/src/dag.rs:188-201, 265-381`): make `curve_recomputes` and `trade_recomputes` consistently report **per-call** counts in both `evaluate_portfolio` and `try_evaluate_portfolio` (while `total_curve_recomputes` and `total_trade_recomputes` report cumulative lifetime counts), and update `DeskWorkstation::step_tick` (`crates/ql-terminal/src/lib.rs:1163, 1165`) to read `summary.last_trade_recomputes` and `summary.last_curve_recomputes > 0`.
   - In `PortfolioDag::value_trade` (`dag.rs:409-444`), resolve `PortfolioTrade::IrsSwap` through `VanillaSwapSpec::resolve -> ResolvedVanillaSwap` (using `DayCounter::Act365Fixed`) and slice `params` at `curves.ois_pillars` instead of hardcoding `3`.
   - Replace `(d - ref_date) as f64 / 365.0` with `day_counter.year_fraction(ref_date, d, None, None)` in `crates/ql-curves/src/bootstrap.rs`, `crates/ql-instruments/src/swap.rs:54-56`, and `crates/ql-instruments/src/bond.rs:71`.
   - Rename `ql_vol::DupireLocalVol` (`crates/ql-vol/src/surface.rs:4`) to `DupireTotalVarianceFormula` (with a type alias in `ql-vol`) so `ql_terminal::prelude` can remove `#[allow(ambiguous_glob_reexports)]`.
3. **`ql-xva` Netting-Set Grouping (`crates/ql-xva/src/cube.rs:70-106`) & `ql-terminal` SIMM / Scenario Wiring (`crates/ql-terminal/src/lib.rs:748-863`)**:
   - In `simulate_hw1f_exposure_cube` (`crates/ql-xva/src/cube.rs:70-106`), `swaps: &[NettingSetSwap]` represents a **single legal netting set** of swaps! Sum the MTM of all `swaps` in the netting set into netting set `0` (`n_netting_sets = 1`), and provide `simulate_hw1f_multi_netting_set_cube(model, curve, netting_sets: &[Vec<NettingSetSwap>], times, n_paths, seed)` for multiple distinct legal netting sets.
   - In `DeskWorkstation::step_tick` (`crates/ql-terminal/src/lib.rs:816-862`):
     1. Scale `summary.total_par_dv01` by `* 1e-4` (1bp DV01, matching `credit_raw` and `inflation`) so IR SIMM margin is not overstated by $10{,}000\times$!
     2. Compute each risk class margin via `compute_isda_simm` and aggregate the 6 risk classes via `compute_isda_simm_6_risk_classes` (`crates/ql-xva/src/simm.rs:18-37`) so `compute_isda_simm_6_risk_classes` is used in production.
   - In `DeskWorkstation::xva_metrics` (`crates/ql-terminal/src/lib.rs:748-766`), use `self.scenario_market()` and bump hazard rate `hc` under `ScenarioKind::CreditWideningPlus100Bp` so CVA/DVA/FVA dynamically respond to scenario changes.
4. **Fix `G2PlusPlus::swaption_price` (`crates/ql-rates/src/g2.rs:90-107, 271-284`) and `ql-credit` (`crates/ql-credit/src/lib.rs:333-337`)**:
   - In `G2PlusPlus::swaption_price` (`g2.rs:271-284`), do NOT `continue` when `brent_solve` returns `Err(NotBracketed)`! Expand `bound` (e.g. `while f_scalar(-bound) > 0.0 || f_scalar(bound) < 0.0 { bound *= 2.0; }`) so deep ITM/OTM quadrature nodes are solved rather than dropped to zero. Also fix `g2.rs:102-105` so small-$a$ and small-$b$ Taylor branches are independent (`|a| < 1e-6` vs `|b| < 1e-6`).
5. **Fix Benchmark Suites 19 & 20 in `crates/ql-benches/src/suites.rs` (`lines 1489-1655`)**:
   - **Suite 19 (`19_american_baw_ju_ift_dual5`, `suites.rs:1489-1547`)**: Remove the `/ 1_000.0` division hack (`lines 1520, 1544`) and tighten `tol` to `1e-5` once `JuQuadraticEngine::price` lifts `star` via `ift_lift_1d`.
   - **Suite 20 (`20_risk_parity_ift_allocation_aad`, `suites.rs:1549-1655`)**: Weight each asset $k \in 0..n$ by `(k + 1) as f64` in both the finite-difference checksum (`lines 1576, 1587-1591, 1605-1609`) and the `Dual<90>` checksum (`lines 1648-1651`):
     ```rust
     checksum += wu
         .iter()
         .zip(&wd)
         .enumerate()
         .map(|(k, (u, d))| (k + 1) as f64 * (u - d) / (2.0 * h))
         .sum::<f64>();
     ```
     and:
     ```rust
     let mut checksum = w.iter().enumerate().map(|(k, x)| (k + 1) as f64 * x.val).sum::<f64>();
     for lane in 0..LANES {
         checksum += w.iter().enumerate().map(|(k, x)| (k + 1) as f64 * x.d[lane]).sum::<f64>();
     }
     ```
     so the 90 Jacobian lanes do not cancel to $\frac{\partial}{\partial\theta}\sum_k w_k \equiv 0$!

---

## Workstream 4: Rust Systems, AD/IFT, Vectorized SIMD & Hot-Loop Allocation Remediation (`RUST-P2-04` – `RUST-P2-05`)

1. **`ift_lift_1d` (`crates/ql-ad/src/ift.rs:22-56`)**:
   - Do NOT call `residual(T::from_f64(x_star))` or `residual(T::from_f64(candidate))` inside `ift_lift_1d`! When `T = TapeVar<f64>`, calling `residual` 2 or 3 times records 2 or 3 copies of the entire pricing graph onto the active `ArenaTape`!
   - Evaluate `let r_eval = residual(x);` **exactly once** inside `ift_lift_1d`:
     ```rust
     pub fn ift_lift_1d<T: Real, F: FnOnce(T) -> T>(x_star: f64, residual: F, f_x: T, f_xx: f64) -> T {
         let fx0 = f_x.val();
         if !fx0.is_finite() || fx0.abs() <= 1e-14 {
             return T::from_f64(x_star);
         }
         let x = T::from_f64(x_star);
         let r_eval = residual(x);
         let d = -r_eval / f_x;
         let d_tan = -(r_eval - T::from_f64(r_eval.val())) / f_x;
         x + d - T::from_f64(0.5 * f_xx) * d_tan * d_tan / f_x
     }
     ```
     Notice how `d = -r_eval / f_x` automatically performs the exact 1-step scalar Newton polish `-r_eval.val() / fx0` on the primal value AND computes the exact first-order tangent in a **single** evaluation of `residual(x)`, while `d_tan` has zero primal so `d_tan * d_tan` is purely second-order (`Dual2`)!
2. **Adopt `TapeScope::guard()` Across All Call Sites & Harden `ArenaTape` (`crates/ql-ad/src/tape.rs`)**:
   - Migrate all `TapeScope::start_recording()` / `stop_recording()` call sites (`bootstrap.rs`, `paths.rs`, `risk_measures.rs`, `ql-credit`, `ql-inflation`, `ql-fx`, `ql-rates`, `ql-commodities`) to `TapeScope::guard()` / `guard_with_capacity(...)`.
   - In `ArenaTape::propagate_range` (`tape.rs:185-191`) and `propagate_to_mark` (`tape.rs:204-210`), skip zero adjoints/weights (`if a == T::default() { continue; }`) to prevent `0.0 * Inf = NaN` in checkpointed AAD.
   - Guard `self.val == 0.0` in `TapeVar::powf` (`tape.rs:602`), `Dual::pow` (`dual.rs:301-312`), `Dual2::pow` (`dual.rs:672`), and `TapeVar::pow` (`tape.rs:608`).
3. **True Vectorized SIMD `simd_erfc` (`crates/ql-ad/src/scalar.rs:200-245`)**:
   - Replace the scalarized `Simd::from_array(x.to_array().map(libm::erfc))` in `simd_erfc4` and `simd_erfc8` (`scalar.rs:234-240`) with a branchless vectorized rational approximation in `Simd<f64, N>` (e.g., Cody's minimax rational approximation on $|z| \le 4$ and asymptotic continued-fraction/rational expansion $\frac{e^{-z^2}}{z\sqrt{\pi}} R(1/z^2)$ on $|z| > 4$ combined via `Simd::simd_le` / `.select(...)`, using `simd_exp(-z * z)`), preserving both `< 1e-12` relative accuracy across $x \in [-8, 8]$ and true 4-lane/8-lane SIMD execution!
4. **Fix Remaining `Date` / `inv_norm_cdf` Edge Cases & Hot-Loop Allocations**:
   - In `crates/ql-time/src/date.rs:160-171`, compute `add_months` in `i64` (`let z = (y as i64) * 12 + (m as i64 - 1) + months as i64;`) and saturate to `Date::MIN` / `Date(i32::MAX)` instead of panicking on `.expect("valid")`. Use `saturating_mul` in `Period::months` (`date.rs:84-85`), `Date::add_period` (`date.rs:186, 189`), and `scale_period` (`schedule.rs:238-242`), and use `i64` day differences in `day_count.rs:25, 131`.
   - In `crates/ql-math/src/distributions.rs:198-202`, guard `if pdf > 0.0 { x -= e / pdf; }` so extreme tail probabilities (`p < 1e-300`) do not divide `0.0 / 0.0 = NaN`.
   - Hoist per-path / per-step scratch buffers out of the hot loops in `mc_gbm_pathwise_tape` (`crates/ql-mc/src/paths.rs:404`), `Fd2dHestonEngine::solve` (`crates/ql-pde/src/adi.rs:189-500`), `cholesky_aad` (`crates/ql-dml/src/cholesky_aad.rs:133-144`), and `ChebyshevTensor2D::eval_barycentric` (`crates/ql-xva/src/chebyshev.rs:158-180`).

---

## Verification Checklist Before Completing

1. `cargo fmt --all -- --check` passes with zero formatting diffs.
2. `cargo clippy --workspace --all-targets -- -D warnings` passes with zero warnings.
3. `cargo test --workspace` passes 100% of unit, integration, and parity tests.
4. `cargo run --release -p ql-benches` executes all 20 benchmark suites cleanly with:
   - Suite 19 (`19_american_baw_ju_ift_dual5`) unscaled (`no / 1_000.0`) and `tol <= 1e-5`,
   - Suite 20 (`20_risk_parity_ift_allocation_aad`) index-weighted `(k + 1) as f64` and `tol <= 1e-8`,
   - Working tree clean (`git status --short` empty).
