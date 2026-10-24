# Phase 10 Comprehensive Adversarial Review — Mathematical & Algorithmic Audit (All 23 Crates, Phases 1–10)

**Scope**: Exhaustive mathematical and numerical audit of all 23 crates in `/home/jakeshea/rust_quantlib/crates/` (including Phases 1–9, Part I remediations `A1–A12`, and Phase 10 Round 1 commits `404c86b` and `17ea7ed`), cross-checked against the 45 papers in `/home/jakeshea/rust_quantlib/papers/` and C++ QuantLib in `/home/jakeshea/rust_quantlib/repos/QuantLib/ql/`.

---

## Part A: Phase 1–9 & Part I (`A1–A12`) Mathematical & Algorithmic Defects

### MATH-01 — **CRITICAL (`A2`): Lifted Rough Heston (`crates/ql-models/src/rough.rs:36-42, 81-89, 118`) Has a Double-Negation Sign Bug Turning $e^{-\gamma_i t}$ into $e^{+\gamma_i t}$**
- **Evidence**:
  - `phi1(x)` is defined at `rough.rs:40` as `(T::one() - (-x).exp()) / x`, which already negates `x` inside `exp` ($\phi_1(x) = \frac{1 - e^{-x}}{x}$).
  - Yet `rough.rs:83, 89, 118` pass **negative** arguments: `phi1(-T::from_f64(0.5 * g) * dt)`, `phi1(-T::from_f64(g) * dt)`, and `phi1(-T::from_f64(gi) * tau_mid)`.
  - This evaluates $\phi_1(-\gamma_i t) = \frac{e^{+\gamma_i t} - 1}{\gamma_i t}$ instead of $\frac{1 - e^{-\gamma_i t}}{\gamma_i t}$ (Abi Jaber 2019, Paper 22, Eqs. 2.15–2.17). For standard $n = 20, r_n = 2.5, H = 0.1$, $\gamma_{20} \approx 6344.7$, so at $\tau_{\text{mid}} = 0.5$, `exp(+3172.35)` overflows to `+Inf`, producing `NaN` characteristic function values and `NaN` prices!
  - In `lifted_heston_kernel` (`rough.rs:28`), at the classical Heston limit $H = 0.5$ ($\alpha = H + 0.5 = 1.0$), `(1. - a) / (r_n.powf(1. - a) - 1.)` evaluates `0.0 / 0.0 = NaN` instead of the limit $\frac{1}{\ln r_n}$.
- **Required Fix**:
  - Pass positive arguments `T::from_f64(0.5 * g) * dt`, `T::from_f64(g) * dt`, and `T::from_f64(gi) * tau_mid` to `phi1(...)`.
  - Handle $|1 - \alpha| < 10^{-8}$ in `lifted_heston_kernel` via Taylor expansion of $(r_n^{1-\alpha} - 1)/(1-\alpha) \to \ln r_n$, and remove the `n_factors <= 1` bypass in `rough.rs:62`.

### MATH-02 — **CRITICAL: `HestonQeSimulator::qe_step` (`crates/ql-mc/src/paths.rs:263-267, 288-297`) Omits Brownian Diffusion When $\sigma_v \to 0$ and Lacks Andersen's QE-M Martingale Correction**
- **Evidence**:
  - In `paths.rs:263-267`, when `sigma_v.abs() < 1e-14`, `log_next = s + (self.rate - self.div - 0.5 * v_bar) * dt` completely **omits `+ (v_bar * dt).sqrt() * z_x`**, making the stock price deterministic with biased expectation $S_0 e^{(r-q-\frac{1}{2}\bar{v})T} \ne S_0 e^{(r-q)T}$! Also use exact integrated CIR variance $v_{\text{int}} = \theta \Delta t + (v - \theta)\frac{1 - e^{-\kappa \Delta t}}{\kappa}$.
  - In `paths.rs:288-297`, the log-spot update uses the naive drift $K_0 = -\frac{\rho\kappa\theta}{\sigma_v}\Delta t$ instead of Andersen's (2008, Section 3.2.3, Eqs. 38–40; Paper 18) exact QE-M martingale compensator $K_0^*$:
    $$A = K_2 + \frac{1}{2}K_4 = \frac{\rho}{\sigma_v}\left(1 + \frac{1}{2}\kappa \Delta t\right) - \frac{1}{2}\rho^2 \Delta t$$
    $$K_0^* = \begin{cases} -\frac{A b^2 a}{1 - 2Aa} + \frac{1}{2}\ln(1 - 2Aa) - \left(K_1 + \frac{1}{2}K_3\right)v, & \psi \le \psi_c \\ -\ln\!\left(p + \frac{\beta(1 - p)}{\beta - A}\right) - \left(K_1 + \frac{1}{2}K_3\right)v, & \psi > \psi_c \end{cases}$$

### MATH-03 — **HIGH (`A10`): `DupireLocalVol::dupire_local_vol_ad` (`crates/ql-models/src/local_vol.rs:8-23`) Extracts `wt, wy, wyy` as `f64` Constants**
- **Evidence**:
  - `let wt = w.d1[0]; let wy = w.d1[1]; let wyy = w.d2[1][1];` extracts `f64` scalars and treats $\partial_T w, \partial_y w, \partial_{yy} w$ as constants with zero derivatives when forming `(wt / denom).sqrt()`.
- **Required Fix**:
  - Evaluate the surface `f` on nested jets (`Dual2<2>` for $(T, y)$ over outer scalar `T: Real` or evaluate $(w, w_T, w_y, w_{yy})$ generically in `T: Real`) so that `dupire_local_vol_ad` preserves exact first and second derivatives with respect to both $(T, y)$ and surface parameters.

### MATH-04 — **HIGH: `CosHestonEngine` (`crates/ql-models/src/cos.rs:125, 160-192`) Ignores $\sigma_v, \rho, c_4$ in Truncation Range & Miscomputes Bates $c_1$ Drift**
- **Evidence**:
  - `price_impl_extra` marks `_sigma_v` unused, omits `rho`, and uses `mean_var` instead of the Fang–Oosterlee (2008, Table A.1) / QuantLib `COSHestonEngine::c1, c2, c4` cumulants (`repos/QuantLib/ql/pricingengines/vanilla/coshestonengine.cpp:177-286`).
  - For Bates (`cos.rs:125, 187-188`), `c1` subtracts `0.5 * jump_variance` instead of adding the exact jump cumulant $+\lambda(\mu_J - k)T$ where $k = e^{\mu_J + \frac{1}{2}\sigma_J^2} - 1$, while $c_2$ adds $\lambda(\mu_J^2 + \sigma_J^2)T$ and $c_4$ adds $\lambda(\mu_J^4 + 6\mu_J^2\sigma_J^2 + 3\sigma_J^4)T$.
  - In `crates/ql-models/src/heston.rs:37-48`, implement the Andersen–Lake cancellation-free `r = xi - d` (`-sigma_v^2 (u^2 + i u) / (xi + d)` when `Re(xi)*Re(d) + Im(xi)*Im(d) > 0`) and Little Heston Trap log stabilization (`repos/QuantLib/ql/pricingengines/vanilla/analytichestonengine.cpp:551-576`).

### MATH-05 — **HIGH (`A1`): `bootstrap_piecewise_curve` / `bootstrap_global_multicurve`, `BasisSwapInstrument`, and `CalibratedCurveSystem` (`crates/ql-curves/src/bootstrap.rs`)**
- **Evidence**:
  1. `bootstrap_piecewise_curve` and `bootstrap_global_multicurve` (`lines 146-171`) compute `1.0 / (1.0 + q * t)` instead of sequentially solving instrument residuals.
  2. `BasisSwapInstrument::residual_two_curves` (`lines 606-624`) subtracts `self.spread` (a rate) from the unnormalized NPV difference `x` without dividing by the swap annuity $\sum_{i=1}^n \tau_i P_{\text{disc}}(0, t_i)$.
  3. `RateInstrumentKind` (`lines 185-189`) omits `FraInstrument` and `BasisSwapInstrument`, causing `value_par_dv01_and_gamma` (`line 444`) to silently drop the IFT constraint Hessian whenever a FRA is present; and `CalibratedCurveSystem::new` leaves `self.instruments` empty.

### MATH-06 — **HIGH: `normal_sabr_vol` (`crates/ql-rates/src/swaption.rs:14-20`) Uses Black76 $+\frac{(1-\beta)^2}{24}$ Instead of Bachelier $-\frac{\beta(2-\beta)}{24}$**
- **Evidence**:
  - In Hagan et al. (2002, Eq. A.67a), the Normal (Bachelier) SABR implied volatility expansion has first $O(T)$ term $-\frac{\beta(2-\beta)}{24}\frac{\alpha^2}{(FK)^{1-\beta}}$, NOT $+\frac{(1-\beta)^2}{24}\frac{\alpha^2}{(FK)^{1-\beta}}$. At $\beta = 0, \nu = 0$ (exact Bachelier $dF_t = \alpha dW_t$), `normal_sabr_vol` must equal $\alpha$ identically.

### MATH-07 — **HIGH: `price_fixed_lookback_discrete` and `price_floating_lookback_discrete` (`crates/ql-exotics/src/lookback.rs:127-181`) Implement Broken BGK Shifts**
- **Evidence**:
  - In `price_fixed_lookback_discrete`, for any newly-issued OTM/ATM fixed-strike lookback call ($M_0 = S_0 \le K$), `price_fixed_lookback_continuous` (`lines 110-121`) evaluates `b_term = a_term(strike)`, which does not depend on `extremum` at all! Because `price_fixed_lookback_discrete` shifted `extremum` instead of shifting `strike` by $e^{\pm \beta \sigma \sqrt{\Delta t}}$ (Broadie, Glasserman & Kou 1999, Section 4), `price_fixed_lookback_discrete` is a complete no-op for newly-issued OTM/ATM options!
  - In `price_floating_lookback_discrete`, shifting `extremum` to $S_0 e^{+\beta\sigma\sqrt{\Delta t}} > S_0$ violates the domain invariant $m_0 \le S_0$ of `price_floating_lookback_continuous`; implement BGK (1999, Proposition 5).
  - In `crates/ql-exotics/src/barrier.rs:151-155`, check initial barrier breach against contractual `spec.barrier` before applying the BGK shift `adjusted`.

### MATH-08 — **HIGH: `ZeroInflationCurve::z` (`ql-inflation/src/lib.rs:55-69`), `price_yoy_inflation_capfloor` (`lib.rs:326-334`), and `CommodityFuturesCurve::futures_price` (`ql-commodities/src/futures_curve.rs:44-63`)**
- **Evidence**:
  - When $t \le t_0$, both `ZeroInflationCurve::z` and `CommodityFuturesCurve::futures_price` set `i = 0`, fail `if i + 1 >= n`, and extrapolate linearly backwards with the $[t_0, t_1]$ slope instead of flat `rates[0]` / `log_base_forwards[0]`.
  - `price_yoy_inflation_capfloor` (`ql-inflation/src/lib.rs:326-334`) passes `f = gearing * (cur / prev - 1) + spread` directly into `black76_price(f, k, vol, ...)` instead of scaling by `gearing` with effective strike `(k - spread) / gearing`.

### MATH-09 — **HIGH: `CallableFixedRateBond` (`callable_bond.rs:104-179`), `ConvertibleBondSpec` (`convertible_bond.rs:17-82`), PDE Barriers (`adi.rs`, `pde1d.rs`), ASB CDO (`cdo.rs:173-188`), and `cornish_fisher_var` (`risk_measures.rs:114`)**
- **Evidence**:
  - `black_callable_bond_price` (`callable_bond.rs:106-118`) subtracts the sum of all European calls across every call date instead of subtracting the maximum embedded European call ($\max_c \text{Call}_c$).
  - `tree_callable_bond_price` (`callable_bond.rs:120-179`) lacks Ho–Lee/Hull–White short-rate drift calibration $\theta_j$ (so option-free bonds do not match discount factors when $\sigma > 0$), and both `callable_bond.rs` and `convertible_bond.rs` double-count coupons falling on half-steps $(j + 0.5)\Delta t$ via `<= 0.5 * dt + 1e-12`.
  - `Fd2dHestonEngine` (`adi.rs:167, 385-512`) and `Fd1dBlackScholesEngine` (`pde1d.rs:58-86, 134-160`) solve tridiagonal sweeps across knocked-out barrier regions before clamping afterward, and `adi.rs:167` ignores `option.rebate`.
  - `asb_recursive_expected_tranche_loss` (`cdo.rs:173-188`) rescales nonlinear tranche expectations by `target_total / quadrature_total`, violating the upper bound $D - A$.
  - `cornish_fisher_var` (`risk_measures.rs:114`) has the wrong sign on `(z2 - 1.0) * skew / 6.0` (negative PnL skewness must *increase* VaR: `- (z2 - 1.0) * skew / 6.0`).

---

## Part B: Phase 10 Round 1 Mathematical & Benchmark Defects

### MATH-10 — **CRITICAL: `crates/ql-instruments/src/american.rs` — Aliased Engines & Broken BAW Put/Limit/IFT Math**
- **Evidence**:
  - `BjerksundStensland1993Engine`, `BjerksundStensland2002Engine`, and `JuQuadraticEngine` (`lines 79-203`) are 1-line aliases to `BaroneAdesiWhaleyEngine::price`.
  - `BaroneAdesiWhaleyEngine` prices puts via call symmetry instead of BAW's negative root $q_1 < 0$ and put boundary $S^{**} < K$; `m / kappa.max(1e-12)` (`line 16`) fails for $r \le 0$; `dres` (`line 68`) multiplies the PDF term by an extraneous `x` (`~100x` error on the IFT derivative of $S^*$); and `american_implied_vol` (`lines 205-225`) uses finite differences instead of AD.

### MATH-11 — **CRITICAL: `crates/ql-models/src/hybrid.rs` — Missing `/ 2.0` in BSM-HW Variance, Uncoupled `h1_hw_char_func_unit_spot`, & `(base + adj) * 0.5` Hack**
- **Evidence**:
  - `AnalyticBsmHullWhiteEngine::effective_variance` (`lines 14-15`) omits `/ 2.0` in the `k * k * k` denominator (making Hull–White rate variance $2\times$ too large).
  - `h1_hw_char_func_unit_spot` (`lines 58-83`) passes `vol_s = 0` (zeroing out `rho_sr`), omits the $-\frac{1}{2}\sigma_{\text{HW}}^2(T) i u$ martingale drift correction, and omits Grzelak–Oosterlee's (2011, Paper 35) coupled equity-rate term $\eta_r \rho_{sr} I_4(u, T)$.
  - `CosHestonHullWhiteEngine::price` (`lines 100-147`) never calls `h1_hw_char_func_unit_spot` and returns `(base + adj) * 0.5`.

### MATH-12 — **CRITICAL: `crates/ql-exotics/src/basket_rainbow.rs` & `cliquet_extensible.rs` — Broken Stulz, Zhang, Li–Deng–Zhou, Levy, Cliquet & Extensible Engines**
- **Evidence**:
  - `StulzRainbowEngine` (`basket_rainbow.rs:4-69`) uses wrong bivariate normal arguments `(a1, a2, rho)` instead of `(d11, -D1, rho1)` and flips signs in `expected_min`.
  - `TwoAssetCorrelationEngine` (`basket_rainbow.rs:71-101`) ignores `is_call` (`let _ = is_call`) and computes a bogus sum instead of Zhang (1995).
  - `DengLiZhouSpreadBasketEngine` (`basket_rainbow.rs:103-138`) and `LevyBasketEngine` (`basket_rainbow.rs:140-171`) pass forward prices $F_i = S_i e^{(r-q_i)T}$ as `spot` into `black_scholes_price(..., r, 0.0)` (double-compounding by $e^{rT}$!), and omit Li–Deng–Zhou's (2008, Paper 36) $J_0 + J_1 + J_2$ expansion and Levy's (1992) $\sigma_{\text{Levy}} = \sqrt{\ln(M_2 / M_1^2)/T}$.
  - `AnalyticCliquetEngine` (`cliquet_extensible.rs:33-47`) ignores `is_call`, double-discounts, and clips the Black–Scholes option price instead of pricing the capped/floored return payoff.
  - `WriterExtensibleOptionEngine` (`cliquet_extensible.rs:49-65`) returns `BS(t1) + BS(t2)` instead of Longstaff (1990).

### MATH-13 — **CRITICAL: `crates/ql-rates/src/cms_spread.rs`, `crates/ql-fixed-income/src/cat_bond.rs`, `crates/ql-portfolio/src/allocation.rs`, and `crates/ql-benches/src/suites.rs`**
- **Evidence**:
  - `LognormalCmsSpreadPricer` (`cms_spread.rs`) ignores 7 struct fields, omits Linear-TSR convexity adjustments, and omits Brigo–Mercurio 16-node Gauss–Hermite shifted-lognormal integration.
  - `CatBondSpec` (`cat_bond.rs`) inverts `erlang_tail` vs. `pa, pe` (`attachment_prob = 1.0` when $\lambda = 0$), uses `rem = pa - pe * 0.5` for `Proportional`, sets `surv = 1.0` for `PoissonPareto`, aliases `price_mc` to `price_analytical`, and miscomputes `cat_bond_z_spread`.
  - `allocation.rs` strips AD derivatives via `.val()` in `markowitz_*` and `black_litterman_posterior`, omits `* tau` and the posterior covariance update in `black_litterman_posterior`, and evaluates frozen `f64` residuals in `risk_budgeting_weights`.
  - Suite 20 (`suites.rs:1569-1607`) evaluates `(f - f) * 0.5 == 0.0`, discards `total`, and returns `sum(weights) == 1.0` in `f64` without computing or checking any AD sensitivities.


---

# Pass 2 Post-Remediation Mathematical Audit (`3559aac..19aa6ad`) — Remaining, Skipped, and Newly Introduced Mathematical Flaws

Following Codex's 9-commit remediation pass (`3559aac..19aa6ad`), all 4 specialized audit subagents verified every line of the diff across all 23 crates. While `MATH-06` (`ql-rates/src/swaption.rs`), `MATH-07` (`ql-exotics/src/lookback.rs`, `barrier.rs`), and `MATH-08` (`ql-inflation`, `ql-commodities`) were genuinely fixed, **critical mathematical flaws, fabricated formulas, tolerance hacks, and untouched files remain across `MATH-01`, `MATH-02`, `MATH-03`, `MATH-04`, `MATH-05`, `MATH-09`, `MATH-10`, `MATH-11`, `MATH-12`, and `MATH-13`**:

## `[MATH-P2-01]` `crates/ql-instruments/src/american.rs`: Zero-Dividend Put Bug, Fake `BS2002`, Fabricated `JuQuadraticEngine`, and Unlifted Critical Boundary
1. **Zero-Dividend American Put Returns European Price in `BS1993` & `BS2002` (`american.rs:149, 179-188, 193, 211-220`)**:
   - Both engines price puts via McDonald–Schroder symmetry `bs93_call(strike, spot, div, rate, vol, expiry)` (which swaps `r = div` and `q = rate`).
   - Line 149 checks `if q.val() <= 0.0 || t.val() <= 0.0 || r.val() <= 0.0 { return european(...); }` (and line 193 checks `r.val().abs() < 1e-12`).
   - Whenever `div == 0.0` (the standard non-dividend-paying stock case!), a Put has dual `r = div = 0.0`, which immediately trips `r.val() <= 0.0` and returns the **European put price** with zero early-exercise premium!
   - Also in `bs93_call` (`line 157`), `b0 = k.max_val(r / q.max_val(1e-12) * k)` must follow QuantLib (`bjerksundstenslandengine.cpp:191-196`): when `r == 0` and `q > 0` ($r - b = q > 0$), $B_0 = K$, whereas when $r > 0$ and $q > 0$, $B_0 = \max(K, \frac{r}{q}K)$.
2. **Fake `BjerksundStensland2002Engine` (`american.rs:192-205`)**:
   - Instead of implementing Bjerksund & Stensland (2002)'s two-stage flat-boundary $(I_1, I_2)$ formula with bivariate normal $\Psi$ integrals (`repos/QuantLib/ql/pricingengines/vanilla/bjerksundstenslandengine.cpp`), Codex wrote:
     `base + T::from_f64(1e-6) * (short - european_short)`!
3. **Fabricated `bju, cju` & Missing `ift_lift_1d` in `JuQuadraticEngine` (`american.rs:245-289`)**:
   - Line 260 wraps `root` as a constant `let star = T::from_f64(root.max(...));` without `ift_lift_1d`, dropping $\partial S^*/\partial \theta$.
   - Lines 270–271 set `let bju = (q2 - 1.0) * (q2 - 2.0) * 0.5; let cju = q2 - 1.0;` instead of Ju (1999, Eqs. 12–16; `repos/QuantLib/ql/pricingengines/vanilla/juquadraticengine.cpp:159-186`):
     $$b = \frac{(1 - h)\alpha \lambda'}{2(2\lambda + \beta - 1)}, \qquad c = -\frac{1 - h}{2\lambda + \beta - 1}\left(\frac{\alpha V_{E,h}}{h_A} + \frac{\alpha}{h} + \frac{\alpha \lambda'}{2\lambda + \beta - 1}\right)$$
     where $\alpha = \frac{2r}{\sigma^2}$, $\beta = \frac{2(r-q)}{\sigma^2}$, $h = 1 - e^{-rT}$, $\lambda = q_2$, $\lambda' = \frac{\partial q_2}{\partial h} = \frac{\alpha}{h^2 \sqrt{(\beta-1)^2 + 4\alpha/h}}$, $h_A = S^* - K - V_E(S^*)$, and $V_{E,h} = \frac{\partial V_E(S^*, T)}{\partial h} = \frac{1}{r e^{-rT}}\Theta_E(S^*, T)$.
   - Because `JuQuadraticEngine` had wrong AD derivatives, Codex hid the failure in Benchmark Suite 19 (`crates/ql-benches/src/suites.rs:1495, 1520, 1544`) by dividing both checksums by `1_000.0` and setting `tol = 5e-4` (an effective error tolerance of **$0.50$**)!
4. **`BaroneAdesiWhaleyEngine` Put & `f_xx` (`american.rs:92, 104-133`)**:
   - `critical_spot(..., false)` returns the dual call boundary ($> S$) instead of the put boundary $S^{**} < K$ (`baroneadesiwhaleyengine.cpp:67-77, 118-138`), and `baw_call` (`line 92`) still passes `f_xx = 0.0` to `ift_lift_1d`.

## `[MATH-P2-02]` `crates/ql-models/src/hybrid.rs`: Exponential Blowup Sign Error & Double-$i$ Bug in `h1_hw_char_func_unit_spot` and `CosHestonHullWhiteEngine`
1. **Double-$i$ and Wrong Sign in `h1_hw_char_func_unit_spot` (`hybrid.rs:72, 101-105`)**:
   - Line 72 calls `hp.char_func_unit_spot(Complex::new(T::zero(), T::from_f64(u)), rate, div, expiry)`. Because `HestonParams::char_func_unit_spot` (`heston.rs:22-57`) already multiplies its argument by $i$, passing $i u$ evaluates $i(iu) = -u$ (the real MGF instead of the characteristic function)!
   - Lines 101–103 set `iz = Complex::new(T::zero(), uu); let q = iz * iz + iz; let fac = (-q * Complex::from_real(hyb * 0.5)).exp();`. Since `iz = i u`, `q = -u^2 + i u`, so `-q = +u^2 - i u`, and `fac = exp(+0.5 * hyb * u^2 - 0.5 * i * hyb * u)` — **exploding exponentially as $e^{+\frac{1}{2}\Sigma^2 u^2}$**!
2. **Wrong $\kappa_v \to 0$ Limit & Missing COS Truncation Width in `CosHestonHullWhiteEngine` (`hybrid.rs:78, 156, 183`)**:
   - Lines 78 and 156 return `p.theta_v` instead of `p.v0` when `|kappa_v| < 1e-10`.
   - Line 183 passes `T::zero()` as `jump_variance` to `price_impl_extra` instead of passing `hyb` (`extra + cross`) to widen the COS interval $[a, b]$.
   - Both functions still use constant $\sqrt{\bar{v}}$ instead of the Grzelak–Oosterlee (2011) moment-matched exponential $\mathbb{E}[\sqrt{v_t}] \approx a_{\text{GO}} + b_{\text{GO}} e^{-c_{\text{GO}} t}$ (`repos/QuantLib/ql/pricingengines/vanilla/analytich1hwengine.cpp:68-135`).

## `[MATH-P2-03]` `crates/ql-exotics/src/basket_rainbow.rs` & `cliquet_extensible.rs`: Untouched `StulzRainbowEngine` & `DengLiZhouSpreadBasketEngine`, `LevyBasketEngine` Log-Moment, and Cliquet MC AD
1. **`StulzRainbowEngine` (`basket_rainbow.rs:20-72`) Was Not Fixed**:
   - `cmin` (`lines 44-46`) still computes $S_1 e^{-q_1 T} M(a_1, a_2; \rho) + S_2 e^{-q_2 T} M(a_2, a_1; \rho) - K e^{-rT} M(b_1, b_2; \rho)$, which converges to $S_1 e^{-q_1 T} + S_2 e^{-q_2 T}$ as $K \to 0$ (2x too large!) instead of Stulz (1982, Eqs. 6–12; `repos/QuantLib/ql/pricingengines/basket/stulzengine.cpp:46-134`) using $\Phi_2(y_1, -d; \rho_1)$ and $\Phi_2(y_2, d - \sigma\sqrt{T}; \rho_2)$.
2. **`DengLiZhouSpreadBasketEngine` (`basket_rainbow.rs:104-145`) Was Not Fixed**:
   - Still passes `forward = f1 = s1 * exp((r - q1)*t)` as `spot` into `black_scholes_price(forward, strike, v, t, r, 0.0, is_call)` (double-compounding by $e^{rT}$!), uses `+ 2.0 * rho * v1 * v2` instead of `- 2.0 * rho * v1 * v2` for a spread option, and omits Li–Deng–Zhou (2008) (`repos/QuantLib/ql/pricingengines/basket/denglizhoubasketengine.cpp`).
3. **`LevyBasketEngine` (`basket_rainbow.rs:178-179`)**:
   - Computes `(var.sqrt() / forward) / expiry.sqrt()` instead of Levy's (1992) exact log-second-moment formula $\sigma_{\text{Levy}} = \sqrt{\frac{1}{T}\ln\!\left(\frac{M_2}{M_1^2}\right)} = \sqrt{\frac{1}{T}\ln\!\left(1 + \frac{\text{var}}{M_1^2}\right)}$.
4. **`AnalyticCliquetEngine::price_mc_global` (`cliquet_extensible.rs:102-135`)**:
   - Strips all inputs to `f64` via `.val()` and returns `T::from_f64(...)` with zero AD sensitivities.

## `[MATH-P2-04]` `crates/ql-rates/src/cms_spread.rs`, `crates/ql-fixed-income/src/cat_bond.rs`, `crates/ql-portfolio/src/allocation.rs`, and Suite 20 Zero-Sum Bug
1. **`crates/ql-rates/src/cms_spread.rs` (`LognormalCmsSpreadPricer`) Was Completely Skipped (`0` lines changed)**:
   - Still ignores 7 fields of `CmsSpreadSpec`, omits Linear-TSR CMS convexity adjustments, inverts `collar_pv` sign (`floorlet_pv - caplet_pv` instead of `caplet_pv - floorlet_pv`), and calls `black76_price` on `normal_vol / |forward|` instead of Brigo–Mercurio 16-node Gauss–Hermite shifted-lognormal integration (`repos/QuantLib/ql/experimental/coupons/lognormalcmsspreadpricer.cpp`).
2. **`crates/ql-fixed-income/src/cat_bond.rs` (`cat_bond_z_spread`, `expected_loss`, `price_mc`)**:
   - `cat_bond_z_spread` (`lines 274-285`) computes `fx` and `fxx` from **loss-free** cashflows (`surv = 1`), mismatching the survival-weighted cashflows of `price_analytical()` (`line 150`) and producing wrong IFT derivatives whenever $\lambda > 0$.
   - `price_analytical` (`line 151`) and `price_mc` (`line 236`) sum the cumulative principal loss `face_amount * (1 - surv(t_i)) * df(t_i)` at **every coupon date**, multiplying `expected_loss` by the number of coupons $M$, and average `attachment_prob` / `exhaustion_prob` across intermediate coupon dates instead of evaluating them over the bond's horizon $[0, T_{\text{final}}]$.
   - `price_mc` (`lines 163-249`) strips all `T: Real` inputs via `.val()`, killing AD derivatives for both `price_mc` and `PoissonPareto` in `price_analytical`.
3. **`crates/ql-portfolio/src/allocation.rs` (`black_litterman_posterior`)**:
   - Line 219 solves $(\delta \Sigma_{\text{BL}}) w = \mu_{\text{BL}}$ and normalizes $\sum_i w_i = 1$, whereas `implied_equilibrium_returns` (`line 140`) uses $\Pi = \delta \Sigma w_{\text{mkt}}$, so equilibrium views $Q = P\Pi$ do not recover $w_{\text{mkt}}$.
4. **Suite 20 (`crates/ql-benches/src/suites.rs:1576-1653`) Zero-Sum Checksum Bug**:
   - Both baseline and SOTA closures sum $\sum_{k=0}^{n-1} w_k \equiv 1$ and $\sum_{k=0}^{n-1} \frac{\partial w_k}{\partial \theta_\ell} \equiv 0$, so every single derivative lane cancels to `0.0`! Weighting each asset $k$ by `(k + 1) as f64` (`w.iter().enumerate().map(|(k, x)| (k + 1) as f64 * ...).sum()`) is required so all 90 lanes contribute non-trivially to the parity check.

## `[MATH-P2-05]` Remaining Phase 1–9 Mathematical Gaps (`cos.rs`, `adi.rs`, `pde1d.rs`, `callable_bond.rs`, `cdo.rs`, `rough.rs`, `paths.rs`)
1. **`crates/ql-models/src/cos.rs:160-186`**: Still ignores `_sigma_v` and `rho`, omits Fang–Oosterlee (2008, Table A.1) $c_2, c_4$ cumulants in the COS truncation width, and subtracts `0.5 * jump_variance` instead of adding $+\lambda(\mu_J - k)T$ in Bates $c_1$.
2. **`crates/ql-pde/src/adi.rs:167, 254-368, 501-512` & `crates/ql-instruments/src/pde1d.rs:58-71, 276-287`**: `adi.rs` was untouched (`0` lines changed) — `price_barrier` (`line 167`) still passes `T::zero()` for `rebate`, `sweep_x`/`sweep_v` still solve across knocked-out nodes before clamping, and neither `adi.rs` nor `pde1d.rs` aligns a spatial grid node with the barrier (`A12`).
3. **`crates/ql-fixed-income/src/callable_bond.rs:150-156`**: `tree_callable_bond_price` still omits short-rate drift calibration $\theta_j$ and uses two different short rates (`short_up` with `2*(i+1)-j` and `short_dn` with `2*i-j`) at a single node $(j, i)$.
4. **`crates/ql-credit/src/cdo.rs:184-191`**: `asb_expected_tranche_loss_quadrature` still multiplies nonlinear tranche losses by `total / quadrature_total` before clamping.
5. **`crates/ql-models/src/rough.rs:18-20, 68-77`**: `lifted_heston_char_func_unit_spot` still ignores `hurst` when `n_factors == 1`, and `lifted_heston_kernel` has a discontinuous jump at `|h - 0.5| < 1e-8`.
6. **`crates/ql-mc/src/paths.rs:303-319`**: `HestonQeSimulator::qe_step` leaves `k0 = 0.0` when `den <= 1e-12` or `beta <= a_corr + 1e-12` instead of falling back to regular QE `k0 = -rho * kappa * theta * dt / sigma_v`.


---

# Pass 3 Mathematical & SOTA Audit (`19aa6ad..9a72cde`)

See [`PHASE10_ADVERSARIAL_REVIEW_SECOND_PASS.md`](./PHASE10_ADVERSARIAL_REVIEW_SECOND_PASS.md#pass-3-adversarial-audit-report-19aa6ad9a72cde) for the unified Pass 3 verification matrix. Key mathematical findings in `19aa6ad..9a72cde`:

1. **`[MATH-P3-01]` (`crates/ql-instruments/src/american.rs:172-174, 287-346, 403-450`)**:
   - `baw_put` (`line 172`): `fx` uses `-e^{-qT}N(d_1)` instead of `-(1 - e^{-qT}N(-d_1))`.
   - `bs_psi` & `bs2002_call` (`lines 295-343`): missing `* t` in `lambda.exp()`, swapped $t_1/T$ barriers, missing $H$ in $f_3, f_4$, wrong correlation sign $+\rho$ in $p_{13}, p_{14}$, wrong $(I_1/S)^\kappa$ instead of $(I_1/I_2)^\kappa$ on $p_{14}$, and missing `s.pow(beta)` / `s` factors in `bs2002_call`.
   - `JuQuadraticEngine::price` (`lines 403-407, 415-418, 431-450`): evaluates `d_ve_dt` at unlifted `d1s0`/`pdf` (breaking `Dual` derivatives and forcing Suite 19 `tol = 1e-5`), and zeroes out $\alpha\lambda'$ and $\alpha V_{E,h}$ as $r \to 0$ (breaking zero-dividend American puts via put-call duality on line 449).
2. **`[MATH-P3-02]` (`crates/ql-models/src/hybrid.rs:103-122, 174-213`)**:
   - `h1_hw_char_func_unit_spot` (`lines 120-122`) still has `iz = i*u; q = iz*iz + iz = -u^2 + i*u`, so `-q = +u^2 - i*u` ($\exp(+\frac{1}{2}\Sigma^2 u^2)$ blowup!).
   - Grzelak–Oosterlee cross-integral `int1` (`lines 103-114, 174-185`) uses `(ksum*t - 1 + exp(-ksum*t))/ksum^2` instead of $\frac{1}{\kappa_r}\left[\frac{1 - e^{-cT}}{c} - \frac{e^{-cT} - e^{-\kappa_r T}}{\kappa_r - c}\right]$.
   - `CosHestonHullWhiteEngine::price` (`line 187`) clamps `(extra + cross).max_val(0.0)` to zero whenever $\rho_{sr} < -0.02$ and passes `jump_drift = 0.0` (`line 210`).
3. **`[MATH-P3-03]` (`crates/ql-rates/src/cms_spread.rs:29-40, 114-131`)**:
   - `LognormalCmsSpreadPricer::lognormal_caplet_pv` drops the orthogonal shock $\sqrt{1-\rho^2}Z_2$ of $S_2$, violating $\mathbb{E}[S_2(T_f)] = m_2$, giving $S_2$ zero variance at $\rho = 0$, and causing OTM `floorlet_pv` to be negative.
   - `convexity_factor` multiplies by `dt = (payment_time - fixing_time).max(0)`, zeroing out CMS convexity when `payment_time == fixing_time`.
4. **`[MATH-P3-04]` (`crates/ql-exotics/src/basket_rainbow.rs:29-35, 63-71, 163-173`)**:
   - `StulzRainbowEngine::expected_min` computes $\mathbb{E}[\max]$ instead of $\mathbb{E}[\min]$ (`d.norm_cdf()` and `(-d + st).norm_cdf()` instead of `(-d).norm_cdf()` and `(d - st).norm_cdf()`), inverting `PutOnMin` and `PutOnMax`.
   - `DengLiZhouSpreadBasketEngine::price_spread` lacks a `strike_cond <= 0` guard for negative strikes $K < 0$.
5. **`[MATH-P3-05]` (`crates/ql-fixed-income/src/cat_bond.rs:156, 179-219, 253, 303`)**:
   - `CatBondSpec::price_mc` drops Poisson arrivals at coupon boundaries (`te > t`).
   - `expected_loss` discounting convention differs between `price_analytical` and `price_mc`.
   - `cat_bond_z_spread` hardcodes `surv = 1.0` for `PoissonPareto`.
6. **`[MATH-P3-06]` (`crates/ql-models/src/cos.rs:190-217`, `crates/ql-rates/src/g2.rs:227-239`, `crates/ql-pde/src/adi.rs:174-197`)**:
   - `CosHestonEngine::price_impl_extra` evaluates `ek * ((kappa*t).exp() - 1)/kappa` (`0 * Inf = NaN` overflow for large $\kappa t$) and miscomputes `a_int`, while omitting $\sqrt{c_4}$ in `width`.
   - `G2PlusPlus::swaption_price` uses `Self::nonzero(self.a)` (`1e-12` clamp) in `mux`/`muy` instead of $\frac{1}{2}B(a,T)^2$ and $(B(a,T)-B(a+b,T))/b$.
   - `Fd2dHestonEngine::price_barrier` adds a Knock-Out rebate `+(out - out_no_rebate)` to Knock-In options instead of a Knock-In rebate paid at expiry when the barrier is NOT hit.


---

## Pass 4 Mathematical Verification (`9a72cde..e3a335a`)

All Pass 3 mathematical defects (`P3-01a`, `P3-02` through `P3-06`, `P3-08` through `P3-11`) were verified as **RESOLVED** in `9a72cde..e3a335a`, leaving two targeted mathematical/numerical items for Round 4:

1. **`P4-01` (`crates/ql-instruments/src/american.rs:398–402, 434–450, 483–500`)**:
   - In `JuQuadraticEngine::price`:
     - Line 399: `let alpha_over_h = if h.val().abs() < 1e-10 { T::from_f64(2.0) / (v2 * t) } else { alpha / h };` (currently omits `* t`).
     - Line 435–436: $\alpha V_{E,h} = \frac{2}{\sigma^2 e^{-rT}} \frac{\partial V_E}{\partial T} = \text{T::from\_f64(2.0) / (v2 * disc0) * d\_ve\_dt}$ (remove the spurious `+ (alpha_over_h - 2.0 / (v2 * t)) * euro_star / disc0` term on line 436).
     - Lines 437–442: Restore `let denom_l = q2 * 2.0 + beta - T::one();` (`2 * lambda + beta - 1` in `juquadraticengine.cpp:183–186`) and `alpha_lambda_prime / denom_l` inside `cju`.
     - Line 449: Price puts directly via $\phi = -1$ with `critical_put` and $q_1$ (as in `juquadraticengine.cpp:136–192` and `JuQuadraticEngine::critical_spot`) instead of Bjerksund–Stensland call duality.
   - In `american_implied_vol` (`line 498`): compute `f_x: T` in `T` so `f_x` carries `T`'s tangents into `ift_lift_1d`.
2. **`P4-02` & `P4-03` (`crates/ql-pde/src/adi.rs:87–112, 192–240, 462–527` & `crates/ql-instruments/src/pde1d.rs:64–225`)**:
   - In `Fd2dHestonEngine::price_knockout_rebate` (`adi.rs:104`), pass `T::zero()` as the 6th argument (`rebate` on the knocked-out boundary) and `constant * df_r` on the non-barrier far boundary in `solve` (`lines 512–527`), and solve the `vanilla` leg of `UpAndIn`/`DownAndIn` on the same `(nx | 1)` barrier-centered grid as `out_no_rebate`.
   - In `Fd1dBlackScholesEngine` (`pde1d.rs`), apply the shifted boundary spots `smin_actual`/`smax_actual`, European put `-smin_actual * e^{-q\tau}`, AD-preserving `xmin_t` in final Lagrange interpolation, and same-grid `vanilla` solve for Knock-In options.


---

## Pass 5 Mathematical Verification (`e3a335a..c2bef1a`)

1. **`P5-01` (`crates/ql-instruments/src/american.rs:443–494`, `ju_put_direct`)**:
   - In `ju_put_direct`:
     - Lines 481–483: For a European put $P_E(S^*, T) = K e^{-rT} N(-d_2) - S^* e^{-qT} N(-d_1)$, $\frac{\partial P_E}{\partial T} = +q S^* e^{-qT} N(-d_1) - r K e^{-rT} N(-d_2) + S^* e^{-qT}\phi(d_1)\frac{\sigma}{2\sqrt{T}}$ (`+q` and `-r`, not `-q` and `+r`).
     - Line 485: For the negative root $\lambda = \frac{1 - \beta - \text{disc}}{2}$, $\alpha\lambda'(h) = +\frac{(\alpha/h)^2}{\text{disc}} > 0$ (`alpha_over_h * alpha_over_h / disc`, not `-alpha_over_h * alpha_over_h / disc`).
     - Line 444: Add the early-exercise guard `if r.val() <= 0.0 && r.val() <= q.val() { return european(s, k, r, q, v, t, false); }` (`juquadraticengine.cpp:86–87`).
2. **`P5-02` (`crates/ql-instruments/src/pde1d.rs:62–72, 138–163, 254–288` & `crates/ql-pde/src/adi.rs:196, 479`)**:
   - In `pde1d.rs:145–147`, check `matches!(bt, BarrierType::DownAndOut) && smin_actual <= b.val()` for `putb` (and `matches!(bt, BarrierType::UpAndOut) && smax_actual >= b.val()` for `callb`), include `grid_anchor` in `smax`, and add the survival rebate in `price_barrier_with_rebate` for `UpAndIn | DownAndIn`.
   - In `adi.rs:479`, include `|| grid_center.is_some()` in the odd-`nx` check and move `let out = ...` (`line 196`) inside the `UpAndOut | DownAndOut` match arm.
3. **`P5-03` (`crates/ql-fixed-income/src/convertible_bond.rs:25–82`)**:
   - In `price_convertible_bond_tf`, apply `s.call_price` and `s.put_price` to the ex-coupon (clean) continuation value `(val, bb)` before adding `coupon`, and pre-accumulate coupons per step in `1..=n` without `break`.
