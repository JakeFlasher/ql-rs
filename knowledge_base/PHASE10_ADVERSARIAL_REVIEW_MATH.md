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
