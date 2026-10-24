# CODEX PHASE 11 IMPLEMENTATION PROMPT — SOTA EXOTIC RATES, SECOND-ORDER IFT SMILE CALIBRATION, BASEL III SA-CCR / KVA / WWR, COMMODITY SWING & GAS STORAGE DP, AND DIFFERENTIABLE CONSTRAINED PORTFOLIO QP

You are implementing **Phase 11** of the `ql-rs` Rust quantitative finance workspace (`/home/jakeshea/rust_quantlib`).

Before writing code, read:
- `knowledge_base/00_INDEX_AND_VERDICT.md`
- `knowledge_base/04_SOTA_QUANT_COMPUTING_PARADIGMS.md` (especially §7 on Second-Order Optimality IFT, Active-Set KKT Differentiation, and Multi-Level Dynamic Programming)
- `knowledge_base/05_BLUEPRINT_MODERN_RUST_QUANTLIB.md`
- `knowledge_base/06_LITERATURE_AND_REFERENCES.md` (papers `43` through `51` in `/home/jakeshea/rust_quantlib/papers/`)
- `knowledge_base/07_BENCHMARK_METHODOLOGY_AND_STUDY.md`

---

## Non-Negotiable Engineering & Architectural Rules

1. **Zero Blanket Lint Suppressions**: Never add `#![allow(...)]` or module-wide `#[allow(dead_code, unused_variables, ...)]`. Every line must pass `cargo +nightly clippy --workspace --all-targets -- -D warnings`.
2. **AD-Generic Kernels (`T: Real`)**: All pricing, risk, and dynamic programming kernels must be generic over `T: Real` (`f64`, `Dual<N>`, `Dual2<N>`, `TapeVar<f64>`). Branch on `.val()` only for continuous/piecewise-differentiable regime selection (or envelope-theorem optimal stopping/active-set selection), while keeping all arithmetic on the selected path in `T`.
3. **Exact IFT at Calibration & Optimization Boundaries**: Never unroll iterative root-finders, Levenberg–Marquardt loops, or projected-gradient iterations onto `ArenaTape`. Solve primal systems in `f64`, compute the exact analytical / `Dual2` optimality Jacobian via the Implicit Function Theorem (Blondel et al. 2021 `papers/43`, Agrawal et al. 2019 `papers/45`), and lift tangents/adjoints in $O(1)$ linear solves.

---

## Deliverable 0 — Two Residual Pass-5 PDE Polish Items (`ql-instruments`, `ql-pde`)

1. **`crates/ql-instruments/src/pde1d.rs` (`Fd1dBlackScholesEngine::price_barrier_with_rebate`)**:
   - Currently, for `BarrierType::UpAndIn | BarrierType::DownAndIn`, `price_barrier_with_rebate` calls `self.fd(..., rebate)` for `out` and returns `vanilla - out`. When `rebate.val().abs() > 0.0`, subtracting `out(rebate)` subtracts the Knock-Out barrier-hit rebate instead of paying a Knock-In expiration rebate when the barrier is **not** hit.
   - Fix: For `UpAndIn | DownAndIn`, compute `out_no_rebate` using `T::zero()` as the rebate. When `rebate.val().abs() > 0.0`, add a `rebate_survival` solve on the identical barrier grid `(smin, smax)` with terminal condition `u(S, T) = rebate` in the un-knocked-out interior and `0` on the knocked-out boundary (and `0` on the far boundary or discounted flat rebate), so `In(rebate) + Out(0) = Vanilla + Rebate * P(tau > T)` and `In(0) + Out(rebate) = Vanilla + Rebate_hit`.
   - Add a regression test in `crates/ql-instruments/tests/` verifying Knock-In + Knock-Out parity with non-zero rebate and positive Knock-In sensitivity $\partial V_{\text{KI}} / \partial \text{rebate} > 0$.

2. **`crates/ql-pde/src/adi.rs` (`Fd2dHestonEngine::price_barrier` & `solve`)**:
   - At line 479 (`Fd2dHestonEngine::solve`), change `let nx = if lower_barrier.is_some() || upper_barrier.is_some()` to `let nx = if lower_barrier.is_some() || upper_barrier.is_some() || grid_center.is_some()` so that when `self.nx` is even, the `vanilla` leg in `price_barrier` (which passes `grid_center: Some(...)` with `None` barriers) uses the exact same odd `(self.nx | 1)` spatial node count as `out_no_rebate`.
   - At line 196 (`Fd2dHestonEngine::price_barrier`), move `let out = self.price_knockout_barrier(params, option);` inside the `BarrierType::UpAndOut | BarrierType::DownAndOut` match arm, and in the `UpAndIn | DownAndIn` arm only solve `rebate_survival` when `option.rebate.abs() > 0.0` (evaluating `0.0` otherwise), cutting `UpAndIn | DownAndIn` from 4 PDE solves to 2 when `rebate == 0.0`.

---

## Deliverable 1 — Bermudan Swaptions under 1F Hull–White Trinomial Tree & 1D FD PDE (`ql-rates`, `ql-portfolio`)

Reference: `repos/QuantLib/ql/pricingengines/swaption/treeswaptionengine.cpp`, `fdhullwhiteswaptionengine.cpp`, `fdmaffinemodelswapinnervalue.hpp`, and Yamakami & Takeuchi (2022) (`papers/50_Yamakami_Takeuchi_2022_Bermudan_Swaption_Two_Factor_Hull_White.pdf`).

### 1.1 Types & Signatures in `crates/ql-rates/src/swaption.rs`
```rust
#[derive(Clone, Debug)]
pub struct BermudanSwaption {
    pub notional: f64,
    pub strike: f64,
    pub exercise_times: Vec<f64>,
    pub swap_times: Vec<f64>,
    pub is_payer: bool,
    pub settlement: SwaptionSettlement,
}

impl BermudanSwaption {
    pub fn from_european(swaption: &EuropeanSwaption) -> Self;
    pub fn from_resolved_swap(
        swap: &ql_instruments::ResolvedVanillaSwap,
        exercise_times: Vec<f64>,
        settlement: SwaptionSettlement,
    ) -> Self;
}

#[derive(Copy, Clone, Debug)]
pub struct HullWhiteBermudanTreeEngine {
    pub steps: usize,
}

#[derive(Copy, Clone, Debug)]
pub struct Fd1dHullWhiteBermudanEngine {
    pub time_steps: usize,
    pub state_points: usize,
}
```
Expose on both `HullWhiteBermudanTreeEngine` and `Fd1dHullWhiteBermudanEngine`:
- `pub fn price<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(&self, model: &ql_models::HullWhite1F<T>, ois: &D, fwd: &F, swaption: &BermudanSwaption) -> T`
- `pub fn price_single_curve<T: Real, C: YieldTermStructure<T>>(&self, model: &ql_models::HullWhite1F<T>, curve: &C, swaption: &BermudanSwaption) -> T` (delegating to `self.price(model, curve, curve, swaption)`)
- `pub fn price_from_context<T: Real + 'static>(&self, model: &ql_models::HullWhite1F<T>, ctx: &ql_curves::MarketContext<T>, discount_id: ql_curves::DiscountCurveId, forward_id: ql_curves::ForwardIndexId, swaption: &BermudanSwaption) -> Option<T>`

### 1.2 Mathematical Specification & Adversarial Guards
1. **Zero-Mean OU State $x_t$ & Exact Bond Reconstruction Free of Finite-Difference $f^M(0, t)$**:
   - Under `HullWhite1F<T>` ($dr_t = (\theta(t) - a r_t)\,dt + \sigma\,dW_t$), write $r_t = x_t + \phi(t)$ where $dx_t = -a x_t\,dt + \sigma\,dW_t$ ($x_0 = 0$).
   - **Never call `curve.instantaneous_forward(t)`** inside the tree or PDE. Instead, reconstruct conditional discount factors at any time $t \ge 0$ and state $x_t$ directly from `ois.discount_factor` (and `fwd.discount_factor`) using the exact affine formula:
     $$P_C(t, T \mid x_t) = \frac{P_C(0, T)}{P_C(0, t)} \exp\!\left(-\frac{1}{2}\sigma^2 B(t, T)\left[B(t, T)\operatorname{Var}_r(t) + B(0, t)^2\right] - B(t, T) x_t\right),$$
     where $B(t, T) = \text{model.b\_factor}(t, T)$ and $\operatorname{Var}_r(t) = \frac{1 - e^{-2at}}{2a}$ (using the Taylor expansion $t - a t^2 + \frac{2}{3}a^2 t^3$ when $|at| < 10^{-6}$).
2. **Dual-Curve Remaining Swap Exercise Payoff at $t = T_{e_k}$**:
   - At exercise date $t = T_{e_k}$, find the first swap period $[T_{i-1}, T_i]$ (`swap_times.windows(2)`) with $T_{i-1} \ge T_{e_k} - 10^{-10}$.
   - For each remaining period $[T_{i-1}, T_i]$ ($\tau_i = T_i - T_{i-1}$), let $\beta_i = \frac{P_{\text{FWD}}(0, T_{i-1}) / P_{\text{OIS}}(0, T_{i-1})}{P_{\text{FWD}}(0, T_i) / P_{\text{OIS}}(0, T_i)}$ (or evaluate the forward rate $F_i(t, x_t) = \frac{1}{\tau_i}\left(\frac{P_{\text{FWD}}(t, T_{i-1} \mid x_t)}{P_{\text{FWD}}(t, T_i \mid x_t)} - 1\right)$).
   - For `SwaptionSettlement::Physical`:
     $$\text{SwapNPV}(t, x_t) = \omega N_{\text{nom}} \sum_{i=i_0}^{n} \left[ \beta_i P_{\text{OIS}}(t, T_{i-1} \mid x_t) - (1 + K \tau_i) P_{\text{OIS}}(t, T_i \mid x_t) \right],$$
     where $\omega = +1$ if `is_payer` else $-1$ (when `fwd == ois`, $\beta_i = 1$ and the floating leg telescopes to $P_{\text{OIS}}(t, T_{\text{start}} \mid x_t) - P_{\text{OIS}}(t, T_n \mid x_t)$).
   - For `SwaptionSettlement::CashParYield`, compute the remaining forward par swap rate $S(t, x_t) = \frac{\sum_{i=i_0}^n (\beta_i P_{\text{OIS}}(t, T_{i-1} \mid x_t) - P_{\text{OIS}}(t, T_i \mid x_t))}{\sum_{i=i_0}^n \tau_i P_{\text{OIS}}(t, T_i \mid x_t)}$ and cash-settled annuity $A_{\text{cash}}(t, x_t) = P_{\text{OIS}}(t, T_{\text{start}} \mid x_t) \sum_{m=1}^{n - i_0 + 1} \frac{\tau_{i_0 + m - 1}}{(1 + \tau_{i_0 + m - 1} S(t, x_t))^m}$, and set $\text{SwapNPV}(t, x_t) = \omega N_{\text{nom}} A_{\text{cash}}(t, x_t) (S(t, x_t) - K)$.
3. **Mandatory Exercise-Date Time Grid Alignment**:
   - Build the time grid on $[0, T_{e_M}]$ (up to the last valid exercise date $T_{e_M} > 0$) by merging $\{0.0\} \cup \{T_{e_k}\}$ with uniform substeps so that **every exercise date $T_{e_k}$ is an exact grid time $t_m$** (never skip an exercise date when `steps` is not a multiple of the exercise frequency!).
4. **`HullWhiteBermudanTreeEngine` (Recombining Trinomial Tree with Analytical Arrow–Debreu Shifts)**:
   - On time grid $0 = t_0 < t_1 < \dots < t_N = T_{e_M}$ with $\Delta t_m = t_{m+1} - t_m$ and $\Delta t_{\max} = \max_m \Delta t_m$:
   - Spatial step $\Delta x = \sqrt{3 V(\Delta t_{\max})}$ where $V(\Delta t) = \sigma^2 \frac{1 - e^{-2a\Delta t}}{2a}$ (or $\sigma^2 \Delta t$ for $|a\Delta t| < 10^{-8}$), and max node index $j_{\max} = \lceil 0.184 / (a\Delta t_{\max}) \rceil$ when $a\Delta t_{\max} > 10^{-4}$ (else $N$).
   - At step $m$ and node $j$ ($x_j = j\Delta x$), conditional mean $\mu_{m,j} = x_j e^{-a\Delta t_m}$, central next-step node $k = \operatorname{round}(\mu_{m,j} / \Delta x)$ clamped to $[-(j_{\max}-1),\, j_{\max}-1]$, residual $e = \mu_{m,j} - k\Delta x$, and moment-matching probabilities in `T`:
     $$p_d = \frac{V(\Delta t_m) + e^2 - e\Delta x}{2\Delta x^2}, \qquad p_m = 1 - \frac{V(\Delta t_m) + e^2}{\Delta x^2}, \qquad p_u = \frac{V(\Delta t_m) + e^2 + e\Delta x}{2\Delta x^2}.$$
   - Forward-propagate Arrow–Debreu prices $Q_{m, j} \in T$ from $Q_{0, 0} = 1$ with the exact analytical shift $\alpha_m \in T$:
     $$\alpha_m = \frac{1}{\Delta t_m}\ln\!\left(\frac{\sum_j Q_{m,j} e^{-x_j \Delta t_m}}{P_{\text{OIS}}(0, t_{m+1})}\right).$$
   - Backward induction from $V_{N, j} = \text{SwapNPV}(T_{e_M}, x_j).\text{max\_val}(T::\text{zero}())$:
     $$V_{m, j}^{\text{cont}} = e^{-(\alpha_m + x_j)\Delta t_m}\left(p_d V_{m+1, k-1} + p_m V_{m+1, k} + p_u V_{m+1, k+1}\right),$$
     and at each exercise step $t_m = T_{e_k}$, $V_{m, j} = V_{m, j}^{\text{cont}}.\text{max\_val}(\text{SwapNPV}(T_{e_k}, x_j))$.
5. **`Fd1dHullWhiteBermudanEngine` (1D Crank–Nicolson + Rannacher FD PDE)**:
   - On spatial grid $x \in [-x_{\max}, x_{\max}]$ with $x_{\max} = 5.5\sqrt{V(T_{e_M})}$, discretize $-\partial_t V = -a x \partial_x V + \frac{1}{2}\sigma^2 \partial_{xx} V - (x + \bar{\phi}_m) V$ with step-integrated shift $\bar{\phi}_m = -\frac{1}{\Delta t_m}\ln\!\left(\frac{P_{\text{OIS}}(0, t_{m+1})}{P_{\text{OIS}}(0, t_m)}\right) + \frac{1}{4}\sigma^2(B(0, t_m)^2 + B(0, t_{m+1})^2)$, linear/zero-gamma boundary conditions at $\pm x_{\max}$, Rannacher startup steps after terminal/exercise payoff kinks, and `ql_math::solve_tridiagonal`.
6. **Integration into `PortfolioDag` (`crates/ql-portfolio/src/dag.rs`)**:
   - Add `PortfolioTrade::BermudanSwaption { id: usize, notional: f64, fixed_rate: f64, exercise_years: Vec<usize>, maturity_years: usize, is_payer: bool, hw_a: f64, hw_sigma: f64 }`.
   - Update `trade_dependency_hash`, trade ID extraction, dirty-trade invalidation (depends on `curve_hash`), and `value_trade` (using `curves.value_and_par_dv01` with `HullWhiteBermudanTreeEngine { steps: 60 }` so Bermudan swaptions contribute exact IFT par-swap DV01s to `PortfolioRiskSummary`).

---

## Deliverable 2 — Differentiable SABR & Arbitrage-Free SVI Smile Calibration with `Dual2<N>` Second-Order Optimality IFT (`ql-vol`)

Reference: Blondel et al. (2021) (`papers/43_Blondel_et_al_2021_Efficient_Modular_Implicit_Differentiation.pdf`), Gatheral & Jacquier (2014) (`papers/46_Gatheral_Jacquier_2014_Arbitrage_Free_SVI_Volatility_Surfaces.pdf`), and `crates/ql-models/src/calibration.rs`.

### 2.1 Near-ATM `Dual2` Taylor Guard in `AdSabrSmile::vol` (`crates/ql-vol/src/surface.rs`)
- In `AdSabrSmile::vol` (around line 113) and `sabr_quote_sensitivities` (around line 242), raise the small-$z$ Taylor threshold from `1e-6` to `1e-4` and include the cubic term of $z / \chi(z)$:
  $$\frac{z}{\chi(z)} = 1 - \frac{1}{2}\rho z + \frac{2 - 3\rho^2}{12} z^2 + \frac{\rho(5\rho^2 - 3)}{24} z^3$$
  so `Dual2<3>` second derivatives with respect to $(\alpha, \rho, \nu)$ are smooth and free of $0/0$ cancellation near ATM ($K \approx F$).

### 2.2 Gatheral–Jacquier (2014) Durrleman Butterfly Function & Density (`crates/ql-vol/src/surface.rs`)
Implement on `AdSviSmile<T>` and as free functions:
- `pub fn svi_durrleman_g<T: Real>(smile: &AdSviSmile<T>, strike: T) -> T`:
  For $k = \ln(K / F)$, $w(k) = \text{smile.total\_variance}(K)$, $d = k - m$, $q = \sqrt{d^2 + \sigma^2}$, $w'(k) = b(\rho + d/q)$, and $w''(k) = b\sigma^2 / q^3$:
  $$g(k) = \left(1 - \frac{k\,w'(k)}{2\,w(k)}\right)^2 - \frac{w'(k)^2}{4}\left(\frac{1}{w(k)} + \frac{1}{4}\right) + \frac{w''(k)}{2}.$$
- `pub fn svi_risk_neutral_density<T: Real>(smile: &AdSviSmile<T>, strike: T) -> T`:
  $$p(K) = \frac{g(k)}{K \sqrt{2\pi w(k)}} \exp\!\left(-\frac{d_-(k)^2}{2}\right), \qquad d_-(k) = -\frac{k}{\sqrt{w(k)}} - \frac{\sqrt{w(k)}}{2}.$$
- `pub fn is_butterfly_arbitrage_free(&self, k_grid: &[f64]) -> bool`: checks $a + b\sigma\sqrt{1-\rho^2} \ge -10^{-12}$, Roger Lee's wing bound $b(1 + |\rho|) \le 2.0 + 10^{-12}$, and $g(k) \ge -10^{-9}$ for all $k \in \text{k\_grid}$.

### 2.3 `CalibratedSabrSmile` & `CalibratedSviSmile` with `Dual2` Second-Order Optimality IFT
```rust
#[derive(Copy, Clone, Debug)]
pub struct SmileStrikeQuote {
    pub strike: f64,
    pub market_vol: f64,
    pub weight: f64,
}

#[derive(Clone, Debug)]
pub struct CalibratedSabrSmile {
    pub smile: AdSabrSmile<f64>,
    pub quotes: Vec<SmileStrikeQuote>,
    pub rmse: f64,
    pub sensitivity: ql_ad::IftSystemSensitivity,
}

#[derive(Clone, Debug)]
pub struct CalibratedSviSmile {
    pub smile: AdSviSmile<f64>,
    pub quotes: Vec<SmileStrikeQuote>,
    pub rmse: f64,
    pub min_durrleman_g: f64,
    pub sensitivity: ql_ad::IftSystemSensitivity,
}
```
1. **`CalibratedSabrSmile::calibrate(forward: f64, expiry: f64, beta: f64, initial: (f64, f64, f64), quotes: &[SmileStrikeQuote]) -> Result<Self, ql_math::MathError>`**:
   - Fit $\theta = (\alpha, \rho, \nu) \in \mathbb{R}^3$ to weighted implied-vol residuals $r_i(\theta, q) = \sqrt{w_i}\left(\hat{\sigma}_{\text{SABR}}(K_i; \theta) - q_i\right)$ via `ql_math::levenberg_marquardt`.
   - At the converged interior optimum $\theta^* = (\alpha^*, \rho^*, \nu^*)$, evaluate `AdSabrSmile<Dual2<3>>` at every quote $K_i$ to form the exact second-order optimality Hessian (Blondel et al. 2021 `papers/43`):
     $$H_{\theta\theta} = \sum_{i=1}^M w_i \left(\nabla_\theta \hat{\sigma}_i (\nabla_\theta \hat{\sigma}_i)^\top + (\hat{\sigma}_i(\theta^*) - q_i) \nabla^2_{\theta\theta} \hat{\sigma}_i(\theta^*)\right), \qquad J_{\theta q}(k, i) = -w_i \frac{\partial \hat{\sigma}_i(\theta^*)}{\partial \theta_k},$$
     and store `IftSystemSensitivity::new(h_opt, j_theta_q)`.
   - Provide:
     - `pub fn implied_vol_ad<T: Real>(&self, strike: T, quote_vols: &[T]) -> T`: lifts $(\alpha, \rho, \nu)$ into `T` via $\theta_k(q) = \theta_k^* + \sum_{j=1}^M \left(\frac{\partial \theta^*}{\partial q}\right)_{k,j} (q_j - q_j^*)$ and evaluates `AdSabrSmile<T>::vol(strike)`.
     - `pub fn vol_and_quote_vega(&self, strike: f64) -> (f64, Vec<f64>)`
     - `pub fn value_and_quote_vega<F>(&self, pricer: F) -> (f64, Vec<f64>) where F: Fn(&AdSabrSmile<ql_ad::Dual<3>>) -> ql_ad::Dual<3>`

2. **`CalibratedSviSmile::calibrate(forward: f64, expiry: f64, initial: (f64, f64, f64, f64, f64), quotes: &[SmileStrikeQuote]) -> Result<Self, ql_math::MathError>`**:
   - Fit $\theta = (a, b, \rho, m, \sigma) \in \mathbb{R}^5$ to weighted implied-vol residuals $r_i(\theta, q) = \sqrt{w_i}\left(\hat{\sigma}_{\text{SVI}}(K_i; \theta) - q_i\right)$ augmented with smooth Gatheral–Jacquier butterfly/wing penalties if any grid point violates $g(k) \ge 0$ or $b(1+|\rho|) \le 2$.
   - Evaluate `AdSviSmile<Dual2<5>>` at $\theta^*$ to assemble the exact $5 \times 5$ optimality Hessian $H_{\theta\theta}$ (including any active penalty residuals so $H_{\theta\theta} = \partial^2 \mathcal{L} / \partial \theta^2$ matches the LM objective identically!) and $5 \times M$ mixed Jacobian $J_{\theta q}$, and store `IftSystemSensitivity::new(h_opt, j_theta_q)`.
   - Provide `implied_vol_ad<T: Real>(&self, strike: T, quote_vols: &[T]) -> T`, `vol_and_quote_vega(&self, strike: f64) -> (f64, Vec<f64>)`, and `value_and_quote_vega<F>(&self, pricer: F) -> (f64, Vec<f64>) where F: Fn(&AdSviSmile<ql_ad::Dual<5>>) -> ql_ad::Dual<5>`.

---

## Deliverable 3 — Basel III SA-CCR Regulatory EAD, Green–Kenyon–Dennis KVA, and Wrong-Way Risk XVA (`ql-xva`)

Reference: BCBS 279 (`repos/ORE/Docs/UserGuide/methods/capital/saccr.tex`), Green, Kenyon & Dennis (2014) (`papers/47_Green_Kenyon_Dennis_2014_KVA_Capital_Valuation_Adjustment.pdf`), Brigo & Vrins (2018) (`papers/48_Brigo_Vrins_2018_Disentangling_Wrong_Way_Risk_CVA.pdf`), and Hull & White (2012).

### 3.1 Basel III SA-CCR (`crates/ql-xva/src/saccr.rs` re-exported in `ql-xva/src/lib.rs`)
```rust
#[derive(Clone, Debug, PartialEq)]
pub enum SaCcrAssetClass {
    InterestRate { ccy: String, start_time: f64, end_time: f64 },
    Fx { pair: String },
    Credit { entity: String, is_index: bool, supervisory_factor: f64 },
    Equity { entity: String, is_index: bool },
    Commodity { hedging_set: String, sub_class: String, supervisory_factor: f64 },
}

#[derive(Copy, Clone, Debug)]
pub struct SaCcrOptionSpec<T: Real> {
    pub underlying_price: T,
    pub strike: T,
    pub expiry: f64,
    pub supervisory_vol: f64,
    pub shift: f64,
    pub is_call: bool,
    pub is_bought: bool,
}

#[derive(Clone, Debug)]
pub struct SaCcrTrade<T: Real> {
    pub mtm: T,
    pub notional: T,
    pub maturity: f64,
    pub asset_class: SaCcrAssetClass,
    pub is_long_primary: bool,
    pub option_spec: Option<SaCcrOptionSpec<T>>,
}

#[derive(Clone, Debug)]
pub struct SaCcrNettingSet<T: Real> {
    pub trades: Vec<SaCcrTrade<T>>,
    pub collateral: T,
    pub is_margined: bool,
    pub threshold: T,
    pub mta: T,
    pub nica: T,
    pub mpor_years: f64,
}

#[derive(Clone, Debug)]
pub struct SaCcrReport<T: Real> {
    pub rc: T,
    pub add_on_ir: T,
    pub add_on_fx: T,
    pub add_on_credit: T,
    pub add_on_equity: T,
    pub add_on_commodity: T,
    pub add_on_aggregate: T,
    pub multiplier: T,
    pub pfe: T,
    pub ead: T,
}

pub fn compute_sa_ccr_ead<T: Real>(ns: &SaCcrNettingSet<T>) -> SaCcrReport<T>;
```
Implement the exact BCBS 279 equations in `T: Real`:
- **Replacement Cost ($\text{RC}$)**:
  - Unmargined: $\text{RC} = (V - C).\text{max\_val}(T::\text{zero}())$ where $V = \sum_i \text{mtm}_i$ and $C = \text{ns.collateral}$.
  - Margined: $\text{RC} = (V - C).\text{max\_val}(\text{ns.threshold} + \text{ns.mta} - \text{ns.nica}).\text{max\_val}(T::\text{zero}())$.
- **Supervisory Duration ($\text{SD}_i$) & Adjusted Notional ($d_i$)**:
  - For `InterestRate { start_time, end_time }` and `Credit`: $S_i = \text{start\_time}.\text{max}(0.0)$, $E_i = \text{end\_time}.\text{max}(S_i + 10.0 / 250.0)$, $\text{SD}_i = \frac{\exp(-0.05 S_i) - \exp(-0.05 E_i)}{0.05}$, and $d_i = \text{notional}_i \times \text{SD}_i$.
  - For `Fx`, `Equity`, `Commodity`: $d_i = \text{notional}_i$.
- **Supervisory Delta ($\delta_i$)**:
  - Linear trade (`option_spec == None`): $+1$ if `is_long_primary` else $-1$.
  - Option trade (`Some(opt)`): let $\omega = +1$ if `opt.is_bought` else $-1$, and $\phi = +1$ if `opt.is_call` else $-1$. With $F' = \text{opt.underlying\_price} + \text{opt.shift}$, $K' = \text{opt.strike} + \text{opt.shift}$, and $\sigma\sqrt{T} = \text{opt.supervisory\_vol}\sqrt{\text{opt.expiry}.\text{max}(10.0/250.0)}$:
    $$\delta_i = \omega \cdot \phi \cdot \left(\phi \cdot \frac{\ln(F'/K') + \frac{1}{2}\sigma^2 T}{\sigma\sqrt{T}}\right).\text{norm\_cdf}().$$
- **Maturity Factor ($\text{MF}_i$)**:
  - Unmargined: $\text{MF}_i = \sqrt{\frac{\min(\max(M_i, 10.0 / 250.0), 1.0)}{1.0}}$.
  - Margined: $\text{MF}_i = 1.5\sqrt{\max(\text{mpor\_years}, 10.0 / 250.0)}$.
- **Asset-Class AddOns**:
  - **IR**: Group by currency `ccy` and maturity bucket $k \in \{1: E_i < 1,\; 2: 1 \le E_i \le 5,\; 3: E_i > 5\}$. Sum $D_{j,k} = \sum_{i \in (j,k)} \delta_i d_i \text{MF}_i$. Aggregate across the 3 buckets using the BCBS 279 quadratic form:
    $$Q_j = D_{j,1}^2 + D_{j,2}^2 + D_{j,3}^2 + 1.4 D_{j,1} D_{j,2} + 1.4 D_{j,2} D_{j,3} + 0.6 D_{j,1} D_{j,3},$$
    $$D_j^{\text{IR}} = \text{if } Q_j.\text{val}() > 0.0 \ \{ Q_j.\text{sqrt}() \} \ \text{else} \ \{ T::\text{zero}() \}, \qquad \text{AddOn}^{\text{IR}} = \sum_j 0.005 \times D_j^{\text{IR}}.$$
  - **FX**: Group by `pair`, $D_j^{\text{FX}} = \sum_{i \in j} \delta_i d_i \text{MF}_i$, $\text{AddOn}^{\text{FX}} = \sum_j 0.04 \times |D_j^{\text{FX}}|$.
  - **Equity**: Group by `entity`, $D_k^{\text{EQ}} = \sum_{i \in k} \delta_i d_i \text{MF}_i$, $\text{AddOn}_k = \text{SF}_k D_k^{\text{EQ}}$ ($\text{SF}_k = 0.20, \rho_k = 0.80$ if `is_index` else $\text{SF}_k = 0.32, \rho_k = 0.50$), and single-factor correlation aggregation $Q^{\text{EQ}} = \left(\sum_k \rho_k \text{AddOn}_k\right)^2 + \sum_k (1 - \rho_k^2) \text{AddOn}_k^2$, $\text{AddOn}^{\text{EQ}} = \sqrt{Q^{\text{EQ}}}$ (guarded at $Q^{\text{EQ}}.\text{val}() > 0$).
  - **Credit**: Same single-factor correlation formula as Equity using `supervisory_factor` and $\rho_k = 0.80$ if `is_index` else $0.50$.
  - **Commodity**: Group by `(hedging_set, sub_class)`, $\text{AddOn}_{j,k} = \text{SF}_k \sum_{i \in (j,k)} \delta_i d_i \text{MF}_i$, aggregate within hedging set $j$ with $\rho = 0.40$: $\text{AddOn}_j^{\text{COM}} = \sqrt{\left(0.40 \sum_k \text{AddOn}_{j,k}\right)^2 + 0.84 \sum_k \text{AddOn}_{j,k}^2}$.
- **Multiplier, PFE, and EAD**:
  - $\text{AddOn}^{\text{agg}} = \text{AddOn}^{\text{IR}} + \text{AddOn}^{\text{FX}} + \text{AddOn}^{\text{CR}} + \text{AddOn}^{\text{EQ}} + \text{AddOn}^{\text{COM}}$.
  - If $(V - C).\text{val}() \ge 0.0$ or $\text{AddOn}^{\text{agg}}.\text{val}() \le 10^{-15}$, $\text{multiplier} = T::\text{one}()$; else $\text{multiplier} = T::\text{one}().\text{min\_val}\!\left(T::\text{from\_f64}(0.05) + T::\text{from\_f64}(0.95) \times \left(\frac{V - C}{T::\text{from\_f64}(1.9) \times \text{AddOn}^{\text{agg}}}\right).\text{exp}()\right)$.
  - $\text{PFE} = \text{multiplier} \times \text{AddOn}^{\text{agg}}$ and $\text{EAD} = T::\text{from\_f64}(1.4) \times (\text{RC} + \text{PFE})$.

### 3.2 Green–Kenyon–Dennis (2014) KVA (`compute_kva` in `crates/ql-xva/src/metrics.rs`)
```rust
#[allow(clippy::too_many_arguments)]
pub fn compute_kva<T: Real>(
    cube: &NpvCube<T>,
    curve: &impl YieldTermStructure<T>,
    csa: Option<&CsaSpec>,
    hc: T,
    hs: T,
    cpty_risk_weight: T,
    capital_ratio: T,
    hurdle_rate: T,
    add_on_profile: Option<&[T]>,
) -> T
```
- At each cube time $t_k$ ($\Delta t_k = t_k - t_{k-1}$), compute the expected collateralized positive exposure $\text{EE}(t_k)$ as in `compute_xva`.
- Let $\text{AddOn}_k = \text{add\_on\_profile.and\_then}(|s| s.\text{get}(k).\text{copied}()).\text{unwrap\_or}(T::\text{zero}())$.
- Compute the expected SA-CCR EAD at step $k$: $\text{EAD}_k = T::\text{from\_f64}(1.4) \times (\text{EE}(t_k) + \text{AddOn}_k)$, regulatory capital $K_{\text{CCR}}(t_k) = T::\text{from\_f64}(12.5) \times \text{cpty\_risk\_weight} \times \text{capital\_ratio} \times \text{EAD}_k$, and integrate:
  $$\text{KVA} = \sum_{k=1}^M \text{hurdle\_rate} \times D(0, t_k) \times e^{-(h_c + h_s) t_k} \times K_{\text{CCR}}(t_k) \times \Delta t_k.$$

### 3.3 Wrong-Way Risk XVA (`compute_xva_with_wwr` in `crates/ql-xva/src/metrics.rs`)
```rust
#[allow(clippy::too_many_arguments)]
pub fn compute_xva_with_wwr<T: Real>(
    cube: &NpvCube<T>,
    curve: &impl YieldTermStructure<T>,
    csa: Option<&CsaSpec>,
    hc: T,
    lc: T,
    hs: T,
    ls: T,
    fs: T,
    wwr_sensitivity: T,
) -> XvaReport<T>
```
- Implement the Hull–White (2012) / ORE dynamic credit exposure-coupled intensity model (`repos/ORE/OREAnalytics/orea/aggregation/dynamiccreditxvacalculator.cpp` & `papers/48`):
  - Track pathwise counterparty survival $S_c(t_{k-1}, p) \in T$ (initialized to $T::\text{one}()$ at $t_0 = 0$) and self survival $S_b(t_k) = (-h_s t_k).\text{exp}()$.
  - At step $k$ ($\Delta t_k = t_k - t_{k-1}$), let $w_{k,p} = \text{pos\_vals}[p]$ be the positive collateralized exposure on path $p$ and $\bar{w}_k = \frac{1}{P}\sum_{p=0}^{P-1} w_{k,p}$ be the path mean.
  - Define centered weights $\beta_{k,p} = \Delta t_k \exp\!\big(\text{wwr\_sensitivity} \times (w_{k,p} - \bar{w}_k)\big)$ (clamping the exponent primal value to $[-20.0, 20.0]$ for overflow safety while keeping arithmetic in `T`).
  - If $\text{wwr\_sensitivity}.\text{val}().\text{abs}() \le 10^{-15}$: set $S_c(t_k, p) = (-h_c t_k).\text{exp}()$ for all $p$ (so `wwr_sensitivity == 0` reproduces `compute_xva` to **0 ULP**!).
  - Otherwise, solve for $u_k = \exp(\tilde{a}_k) > 0$ in
    $$g_k(u_k) = \frac{1}{P}\sum_{p=0}^{P-1} S_c(t_{k-1}, p) \exp(-u_k \beta_{k,p}) - \exp(-h_c t_k) = 0$$
    starting from $u_k = h_c / \left(\frac{1}{P}\sum_p \frac{S_c(t_{k-1}, p)}{\exp(-h_c t_{k-1})} \frac{\beta_{k,p}}{\Delta t_k}\right)$ via 8 Newton iterations $u_k \leftarrow (u_k - g_k(u_k)/g_k'(u_k)).\text{max\_val}(T::\text{from\_f64}(10^{-12}))$ directly in `T: Real`.
  - Update $S_c(t_k, p) = S_c(t_{k-1}, p)\exp(-u_k \beta_{k,p})$ (satisfying $\frac{1}{P}\sum_p S_c(t_k, p) = \exp(-h_c t_k)$ to machine precision!) and accumulate pathwise CVA, DVA, and FVA:
    $$\text{CVA} \mathrel{+}= l_c D(0, t_k) S_b(t_{k-1}) \frac{1}{P}\sum_{p=0}^{P-1} E^+(t_k, p)\big(S_c(t_{k-1}, p) - S_c(t_k, p)\big),$$
    $$\text{DVA} \mathrel{-}= l_s D(0, t_k) (S_b(t_{k-1}) - S_b(t_k)) \frac{1}{P}\sum_{p=0}^{P-1} E^-(t_k, p) S_c(t_{k-1}, p),$$
    $$\text{FVA} \mathrel{+}= f_s D(0, t_k) \Delta t_k S_b(t_k) \frac{1}{P}\sum_{p=0}^{P-1} \big(E^+(t_k, p) + E^-(t_k, p)\big) S_c(t_k, p).$$

---

## Deliverable 4 — Commodity Swing Options & Natural Gas Storage via 1F OU "Forest of Trees" Dynamic Programming (`ql-commodities`)

Reference: Warin (2012) (`papers/49_Warin_2012_Hedging_Swing_Contract_on_Gas_Market.pdf`), Bäuerle & Rieß (2015) (`papers/51_Baeuerle_Riess_2015_Gas_Storage_Valuation_Multinomial_Tree_LSM.pdf`), and `repos/QuantLib/ql/methods/finitedifferences/stepconditions/fdmsimpleswingcondition.cpp` / `fdmsimplestoragecondition.cpp`.

### 4.1 Types & Signatures in `crates/ql-commodities/src/swing_storage.rs` (re-exported in `ql-commodities/src/lib.rs`)
```rust
#[derive(Clone, Debug)]
pub struct SwingOptionSpec<T: Real> {
    pub exercise_times: Vec<f64>,
    pub strike: T,
    pub unit_volume: T,
    pub min_exercises: usize,
    pub max_exercises: usize,
    pub is_call: bool,
}

#[derive(Clone, Debug)]
pub struct GasStorageSpec<T: Real> {
    pub decision_times: Vec<f64>,
    pub max_capacity: f64,
    pub initial_inventory: f64,
    pub max_injection_rate: f64,
    pub max_withdrawal_rate: f64,
    pub injection_loss: T,
    pub withdrawal_loss: T,
    pub injection_cost: T,
    pub withdrawal_cost: T,
    pub inventory_grid_points: usize,
}

pub fn price_swing_option<T: Real>(
    spec: &SwingOptionSpec<T>,
    curve: &CommodityFuturesCurve<T>,
    kappa: T,
    sigma: T,
    rate: T,
    steps_per_interval: usize,
) -> T;

pub fn price_gas_storage<T: Real>(
    spec: &GasStorageSpec<T>,
    curve: &CommodityFuturesCurve<T>,
    kappa: T,
    sigma: T,
    rate: T,
) -> T;
```

### 4.2 Mathematical Specification
1. **Futures-Calibrated 1F Ornstein–Uhlenbeck Trinomial Spot Tree**:
   - Build a recombining trinomial tree for $dx_t = -\kappa x_t\,dt + \sigma\,dW_t$ ($x_0 = 0$) across the mandatory exercise/decision dates.
   - Forward-propagate the unconditional node probabilities $P_{m, j} \in T$ from $P_{0, 0} = 1$, and define the spot price at step $m$, node $j$ as:
     $$S_{m, j} = F(0, t_m) \frac{e^{x_j}}{\sum_l P_{m, l} e^{x_l}}, \qquad F(0, t_m) = \text{curve.futures\_price}(T::\text{from\_f64}(t_m)).$$
     This guarantees $\mathbb{E}[S_{t_m}] = F(0, t_m)$ to machine precision for all $(\kappa, \sigma)$ and propagates exact `T: Real` AD sensitivities to `CommodityFuturesCurve<T>` pillars.
2. **Swing Option 2D Backward Dynamic Programming (`price_swing_option`)**:
   - State $(j, u)$ where $u \in \{0, \dots, N_{\max}\}$ is the number of rights **already exercised** ($N_{\min} = \text{min\_exercises}, N_{\max} = \text{max\_exercises}$).
   - At exercise date $T_{e_k}$ ($k \in \{0, \dots, M-1\}$, with $d = M - k$ remaining exercise opportunities including $T_{e_k}$), the signed cashflow of exercising one right is:
     $$\text{cf}(S_{m, j}) = \text{if spec.is\_call } \{ \text{unit\_volume} \times (S_{m, j} - \text{strike}) \} \text{ else } \{ \text{unit\_volume} \times (\text{strike} - S_{m, j}) \}.$$
   - For each $u \in \{0, \dots, N_{\max}\}$:
     - If $u == N_{\max}$: $V_m(j, u) = V_m^{\text{cont}}(j, u)$.
     - Else if $u + d \le N_{\min}$ (**mandatory take-or-pay exercise**): $V_m(j, u) = \text{cf}(S_{m, j}) + V_m^{\text{cont}}(j, u + 1)$.
     - Else (**optional exercise**): $V_m(j, u) = V_m^{\text{cont}}(j, u).\text{max\_val}\!\big(\text{cf}(S_{m, j}) + V_m^{\text{cont}}(j, u + 1)\big)$.
   - Return $V_0(0, 0)$.
   - **Exact Invariants**:
     - When $N_{\min} = N_{\max} = M$, `price_swing_option` must equal the discounted forward strip $\sum_{k=1}^M e^{-r T_{e_k}} (\pm q)(F(0, T_{e_k}) - K)$ to machine precision ($< 10^{-11}$)!
     - Monotonicity: $V(N_{\min}, N_{\max} + 1) \ge V(N_{\min}, N_{\max})$ and $V(N_{\min} + 1, N_{\max}) \le V(N_{\min}, N_{\max})$.
3. **Natural Gas Storage 2D Backward Dynamic Programming (`price_gas_storage`, Bäuerle & Rieß 2015 `papers/51`)**:
   - Discretize inventory $v \in [0, V_{\max}]$ on $G + 1$ uniform levels $v_g = g \Delta v$ ($G = \text{inventory\_grid\_points}.\text{max}(2)$).
   - At decision step $t_m$ ($\Delta t_m = t_{m+1} - t_m$), from inventory level $v_g$, admissible target inventory is $v' \in [v_{\min}(v_g),\, v_{\max}(v_g)]$ with $v_{\min}(v_g) = (v_g - r_{\text{with}}\Delta t_m).\text{max}(0.0)$ and $v_{\max}(v_g) = (v_g + r_{\text{inj}}\Delta t_m).\text{min}(V_{\max})$.
   - For candidate transition $\Delta v = v' - v_g$:
     - Injection ($\Delta v > 0$): cashflow $h(S_{m, j}, \Delta v) = -\big((T::\text{one}() + w_{\text{inj}}) S_{m, j} + c_{\text{inj}}\big)\Delta v$.
     - Hold ($\Delta v = 0$): $h = 0$.
     - Withdrawal ($\Delta v < 0$): cashflow $h(S_{m, j}, \Delta v) = +\big((T::\text{one}() - w_{\text{with}}) S_{m, j} - c_{\text{with}}\big)(-\Delta v)$.
   - Maximize $h(S_{m, j}, v' - v_g) + U_m^{\text{cont}}(v', S_{m, j})$ over all grid levels $v_{g'} \in [v_{\min}(v_g), v_{\max}(v_g)]$ plus the two endpoints $v_{\min}(v_g), v_{\max}(v_g)$ (using linear interpolation on $U_m^{\text{cont}}(\cdot, S_{m, j})$). Interpolate at $t = 0$ for `initial_inventory`.

---

## Deliverable 5 — Box-Constrained Markowitz Portfolio Optimization with Exact Active-Set KKT Implicit Differentiation (`ql-portfolio`)

Reference: Agrawal et al. (2019) (`papers/45_Agrawal_et_al_2019_Differentiable_Convex_Optimization_Layers.pdf`) and Blondel et al. (2021) (`papers/43_Blondel_et_al_2021_Efficient_Modular_Implicit_Differentiation.pdf`).

### 5.1 Types & Signatures in `crates/ql-portfolio/src/allocation.rs`
- Add `InfeasibleBounds` variant to `AllocationError`.
```rust
#[derive(Clone, Debug)]
pub struct BoxConstrainedMarkowitzSpec<T: Real> {
    pub covariance: Vec<Vec<T>>,
    pub expected_returns: Vec<T>,
    pub risk_aversion: T,
    pub lower_bounds: Vec<T>,
    pub upper_bounds: Vec<T>,
}

pub fn try_box_constrained_markowitz_weights<T: Real>(
    spec: &BoxConstrainedMarkowitzSpec<T>,
) -> Result<Vec<T>, AllocationError>;

pub fn box_constrained_markowitz_weights<T: Real>(
    spec: &BoxConstrainedMarkowitzSpec<T>,
) -> Vec<T>;
```
### 5.2 Active-Set Primal Solver & Reduced-KKT `T: Real` Lift
$$\min_{w \in \mathbb{R}^N} \frac{1}{2}\gamma w^\top \Sigma w - \mu^\top w \quad \text{s.t.} \quad \sum_{i=1}^N w_i = 1, \quad l_i \le w_i \le u_i.$$
1. **Validation**: Verify dimensions (`n > 0`, `covariance` is $n \times n$, vectors have length $n$), `risk_aversion.val() > 0.0`, every $l_i.\text{val}() \le u_i.\text{val}() + 10^{-12}$, $\sum_i l_i.\text{val}() \le 1.0 + 10^{-12}$, and $\sum_i u_i.\text{val}() \ge 1.0 - 10^{-12}$ (else `Err(AllocationError::InfeasibleBounds)`).
2. **Primal Active-Set KKT Solver in `f64`**:
   - Solve the convex QP in `f64` to machine precision ($< 10^{-13}$) using primal-dual active-set iteration (or high-precision projected gradient + active-set KKT polish) to determine the exact partition of $\{0, \dots, N-1\}$ into lower-active $\mathcal{B}_L = \{i : w_i^* \le l_i + 10^{-9}\}$, upper-active $\mathcal{B}_U = \{i : w_i^* \ge u_i - 10^{-9}\}$, and strictly interior free indices $\mathcal{F}$.
3. **Exact `T: Real` KKT Implicit Differentiation (`solve_linear_generic::<T>`)**:
   - Set $w_i = \text{spec.lower\_bounds}[i]$ for $i \in \mathcal{B}_L$ and $w_i = \text{spec.upper\_bounds}[i]$ for $i \in \mathcal{B}_U$ (preserving exact `T` tangents with respect to active bounds!).
   - Let $w_{\text{fixed\_sum}} = \sum_{j \in \mathcal{B}_L \cup \mathcal{B}_U} w_j$.
   - If $|\mathcal{F}| == 0$: return `Ok(w)`.
   - If $|\mathcal{F}| == 1$ ($\mathcal{F} = \{i_0\}$): set $w_{i_0} = T::\text{one}() - w_{\text{fixed\_sum}}$ directly and return `Ok(w)`.
   - If $|\mathcal{F}| \ge 2$: assemble the $(|\mathcal{F}| + 1) \times (|\mathcal{F}| + 1)$ equality KKT system in `T`:
     $$M = \begin{bmatrix} \gamma \Sigma_{\mathcal{F},\mathcal{F}} & \mathbf{1} \\ \mathbf{1}^\top & 0 \end{bmatrix}, \qquad \text{rhs} = \begin{bmatrix} \mu_{\mathcal{F}} - \gamma \Sigma_{\mathcal{F},\mathcal{B}} w_{\mathcal{B}} \\ T::\text{one}() - w_{\text{fixed\_sum}} \end{bmatrix},$$
     solve via `solve_linear_generic(&m, &rhs)?`, scatter $w_{\mathcal{F}}$ into `w`, and return `Ok(w)`.
   - When $l_i \ll w_{\text{unconstrained}, i}^* \ll u_i$ for all $i$ ($\mathcal{B}_L = \mathcal{B}_U = \emptyset$), this reduces **identically** to the unconstrained Markowitz KKT solution!

---

## Deliverable 6 — Heston / Bates Variance-of-Variance, Volatility Swap Fair Strike, & Corridor Variance Swaps (`ql-exotics`)

Reference: Broadie & Jain (2008) (*"The Effect of Jumps and Discrete Sampling on Volatility and Variance Swaps"*), Brockhaus & Long (2000), and Gatheral (2006).

Implement in `crates/ql-exotics/src/variance_swap.rs`:
1. **`pub fn heston_integrated_variance_var<T: Real>(v0: T, kappa: T, theta: T, sigma_v: T, expiry: T) -> T`**:
   - Computes $\operatorname{Var}\!\left(\frac{1}{T}\int_0^T v_t\,dt\right)$ under Heston.
   - Let $u = \kappa \times \text{expiry}$.
   - For $|u.\text{val}()| < 10^{-3}$, use the exact 4th-order Taylor series to prevent $0/0$ cancellation in $\kappa^3 T^2$:
     $$\operatorname{Var}(V_T) = \sigma_v^2 T \left[ (v_0 - \theta)\left(\frac{1}{3} - \frac{u}{3} + \frac{11}{60}u^2 - \frac{13}{180}u^3\right) + \theta\left(\frac{1}{3} - \frac{u}{4} + \frac{7}{60}u^2 - \frac{1}{24}u^3\right) \right].$$
   - Otherwise:
     $$\operatorname{Var}(V_T) = \frac{\sigma_v^2}{\kappa^3 T^2}\left[ (v_0 - \theta)\left(1 - 2u e^{-u} - e^{-2u}\right) + \frac{\theta}{2}\left(2u - 3 + 4e^{-u} - e^{-2u}\right) \right].$$
2. **`pub fn heston_volatility_swap_fair_strike<T: Real>(v0: T, kappa: T, theta: T, sigma_v: T, expiry: T) -> T`**:
   - Combines `heston_analytical_fair_variance(v0, kappa, theta, expiry)` and `heston_integrated_variance_var(v0, kappa, theta, sigma_v, expiry)` via `volatility_swap_convexity_adjustment`.
3. **`pub fn bates_analytical_fair_variance_qv<T: Real>(params: &ql_models::BatesParams<T>, expiry: T) -> T`** and **`pub fn bates_analytical_fair_variance_log<T: Real>(params: &ql_models::BatesParams<T>, expiry: T) -> T`**:
   - True quadratic variation strike (Broadie & Jain 2008 Eq. 3.10):
     $$K_{\text{var, QV}}^{\text{Bates}} = \text{heston\_analytical\_fair\_variance}(v_0, \kappa, \theta, T) + \lambda\left(\mu_J^2 + \sigma_J^2\right).$$
   - Log-contract replication strike (Broadie & Jain 2008 Eq. 3.12):
     $$K_{\text{var, log}}^{\text{Bates}} = \text{heston\_analytical\_fair\_variance}(v_0, \kappa, \theta, T) + 2\lambda\left(e^{\mu_J + \frac{1}{2}\sigma_J^2} - 1 - \mu_J\right).$$
4. **`pub fn replicate_corridor_variance_swap_fair_strike<T: Real>(spot: T, rate: T, div_yield: T, expiry: T, cut_strike: T, lower_barrier: T, upper_barrier: T, quotes: &[VarianceSwapStrikeQuote<T>]) -> (T, T)`**:
   - Clips the static OTM option replication integral (and endpoint log-contract tails when the corridor extends beyond $[K_0, K_n]$) to $[B_{\text{lo}}, B_{\text{hi}}] = [\text{lower\_barrier}, \text{upper\_barrier}]$, with clamped forward/cut drift $\Delta_{\text{drift}} = F(1/\tilde{F} - 1/\tilde{K}_{\text{cut}}) - \ln(\tilde{K}_{\text{cut}}/\tilde{F})$ where $\tilde{F} = \operatorname{clamp}(F, B_{\text{lo}}, B_{\text{hi}})$ and $\tilde{K}_{\text{cut}} = \operatorname{clamp}(K_{\text{cut}}, B_{\text{lo}}, B_{\text{hi}})$.
   - When $B_{\text{lo}} \to 0$ and $B_{\text{hi}} \to \infty$, this must match `replicate_variance_swap_fair_strike` to machine precision; when $B_{\text{lo}} == B_{\text{hi}}$, it returns `(0, 0)`; and adjacent corridors $[B_1, B_2] \cup [B_2, B_3]$ with $K_{\text{cut}} = F$ sum to $[B_1, B_3]$.

---

## Deliverable 7 — Benchmark Suites 21 & 22 (`ql-benches`) and Terminal Workstation Tabs 19 & 20 (`ql-terminal`)

1. **`crates/ql-benches/src/suites.rs` & `crates/ql-benches/tests/benches_suite.rs`**:
   - Add two new benchmark suites to `run_all_benchmarks`:
     - **`21_bermudan_hw1f_tree_ift_dv01`**:
       - *Legacy baseline*: 14 central finite-difference bump-and-rebootstrap / re-solve runs (6 OIS par swap quotes + 6 forwarding par swap quotes + 2 Hull–White $(a, \sigma)$ parameters, 28 tree evaluations + 24 multi-curve bootstraps) of a 5-exercise-date Bermudan swaption on `HullWhiteBermudanTreeEngine`.
       - *`ql-rs` SOTA*: 1 cached `CalibratedCurveSystem` IFT Jacobian + 1 `Dual<8>` (or `TapeVar`) evaluation of `HullWhiteBermudanTreeEngine::price` projecting exact par-swap DV01s and Hull–White $(a, \sigma)$ sensitivities, asserting parity $< 10^{-4}$.
     - **`22_sabr_svi_ift_and_qp_kkt_aad`**:
       - *Legacy baseline*: 16 central finite-difference strike-quote recalibrations of an 8-strike `CalibratedSabrSmile` + $2(N + N(N+1)/2)$ perturbed primal solves of an $N=10$ `BoxConstrainedMarkowitzSpec`.
       - *`ql-rs` SOTA*: 1 `CalibratedSabrSmile` with `Dual2<3>` second-order optimality IFT Jacobian + 1 `try_box_constrained_markowitz_weights` active-set reduced-KKT `TapeVar` pass, asserting parity $< 10^{-6}$.
   - **Update `crates/ql-benches/tests/benches_suite.rs:9`** from `assert_eq!(rows.len(), 20);` to `assert_eq!(rows.len(), 22);`!
   - Run `cargo run -p ql-benches --release` to regenerate `crates/ql-benches/BENCHMARK_RESULTS.md` and `knowledge_base/08_BENCHMARK_RESULTS.md`, and update the empirical table in `knowledge_base/07_BENCHMARK_METHODOLOGY_AND_STUDY.md`.

2. **`crates/ql-terminal/src/lib.rs` & `crates/ql-terminal/tests/`**:
   - Expand `TabId` and `TabId::ALL` from 18 to **20 tabs** by adding:
     - `TabId::BermudanVolIft` (`"19 Bermudan / Vol IFT"`) — renders Bermudan vs. European swaption term structure, `CalibratedSabrSmile` / `CalibratedSviSmile` Durrleman $g(k)$ butterfly density & IFT quote Vegas.
     - `TabId::SaCcrSwingQp` (`"20 SA-CCR / Swing / QP"`) — renders Basel III SA-CCR EAD / Green–Kenyon–Dennis KVA / Wrong-Way Risk CVA ladder, Commodity Swing & Gas Storage DP surface, and Box-Constrained Markowitz active-set KKT weights.
   - **Respect the 20-Row Bottom Quadrant Budget (`crates/ql-terminal/src/lib.rs:1405–1475`)**:
     - `render_dashboard` prints `for i in 0..20` for `left_bottom` and `right_bottom`. Place new summary rows in `right_bottom` (which currently has 15 rows out of 20) or condense existing rows so neither `left_bottom.len()` nor `right_bottom.len()` exceeds 20!

---

## Mandatory Verification Gates

Before finishing, run all four commands and verify zero warnings and 100% test pass rate:
```bash
cargo +nightly fmt --all -- --check
cargo +nightly clippy --workspace --all-targets -- -D warnings
cargo +nightly test --workspace --release --no-fail-fast
cargo run -p ql-benches --release
```
