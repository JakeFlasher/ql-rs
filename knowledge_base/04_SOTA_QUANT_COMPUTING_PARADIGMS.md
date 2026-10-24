# 04 — State-of-the-Art (SOTA) Quantitative Computing Paradigms (2015–2026)

A next-generation Rust quantitative finance library must not merely re-implement the 2000-era algorithms of C++ `QuantLib` in safer syntax. Over the past two decades (2005–2026), computational finance at tier-1 investment banks, hedge funds, and academic labs has undergone a paradigm shift driven by four interlocking revolutions:

1. **Adjoint Algorithmic Differentiation (AAD) & Implicit Function Theorem (IFT) Calibration**
2. **Differential Machine Learning (DML) & Differential PCA**
3. **Vectorized (SIMD/GPU) PDE & Monte Carlo Solvers** (ADI Splitting, Two-stage parallel AMC, Quasi-Monte Carlo)
4. **Compositional Contract Algebras & Bytecode Payoff Compilation with Static Fuzzy Smoothing**

All primary literature cited below is archived locally in [`/home/jakeshea/rust_quantlib/papers/`](file:///home/jakeshea/rust_quantlib/papers/).

---

## 1. Algorithmic Differentiation (AD) in Quantitative Finance

### 1.1 The Bump-and-Revalue Bottleneck in Legacy QuantLib
In C++ `QuantLib`, computing sensitivities ($\Delta, \Gamma, \text{Vega}$, or a 500-pillar multi-curve DV01 + SABR surface Vega matrix) requires **finite-difference bump-and-revalue**:
1. Perturb market quote $q_i \leftarrow q_i + \delta q_i$ via `SimpleQuote::setValue`.
2. Trigger `Observable::notifyObservers()` across all dependent curves, surfaces, and engines.
3. Re-bootstrap the yield curve or re-calibrate the model.
4. Re-price the instrument (or re-run the Monte Carlo / PDE solver).

For $N$ market inputs, bump-and-revalue requires $\mathcal{O}(N)$ full pricing passes for first-order Greeks and $\mathcal{O}(N^2)$ for full Cross-Gamma. When $N \in [100, 2000]$ (typical for a multi-curve STIR/swap book or an XVA netting set), real-time Greeks become computationally intractable.

### 1.2 Forward-Mode vs. Reverse-Mode (Adjoint) AD Complexity Crossover
As established in Giles & Glasserman (2005) ([`01_Giles_Glasserman_2005_Smoking_Adjoints_Fast_MC_Greeks.pdf`](file:///home/jakeshea/rust_quantlib/papers/01_Giles_Glasserman_2005_Smoking_Adjoints_Fast_MC_Greeks.pdf)) and Margossian (2019) ([`08_Margossian_2019_Review_of_Automatic_Differentiation.pdf`](file:///home/jakeshea/rust_quantlib/papers/08_Margossian_2019_Review_of_Automatic_Differentiation.pdf)):
- For a pricing function $f: \mathbb{R}^N \to \mathbb{R}^M$:
  - **Forward-Mode AD (JVP, $\dot{y} = J \dot{x}$)** propagates tangents alongside primals in a single forward pass with **zero tape memory**. Each seed costs $\sim 1.5\times\text{–}2.5\times$ the primal evaluation, scaling as $\mathcal{O}(N)$ passes for all $N$ inputs. Using SIMD-packed multi-tangent dual numbers (`Dual<N>` with $N \in \{4, 8, 16\}$), $N$ tangents are propagated simultaneously in registers.
  - **Reverse-Mode / Adjoint AD (VJP, $\bar{x} = J^\top \bar{y}$)** records the directed acyclic execution graph (tape) during the forward pass and propagates adjoints backward from outputs to inputs in $\mathcal{O}(M)$ passes. For scalar price output ($M = 1$), **all $N$ sensitivities are computed in a single backward pass in $\le 4\times$ the cost of one primal evaluation, independent of $N$**.

| Workload | Inputs ($N$) | Outputs ($M$) | Optimal AD Strategy |
| :--- | :---: | :---: | :--- |
| Black-Scholes / Analytic Greeks ($\Delta, \Gamma, \mathcal{V}, \Theta, \rho$) | $4\text{–}6$ | $1$ | Stack-allocated Forward `Dual<8>` or `Dual2` (Hyper-Dual) |
| Yield Curve Bootstrap Jacobian ($J = \partial r_{\text{inst}} / \partial d_{\text{node}}$) | $30\text{–}80$ | $30\text{–}80$ | Sparsity-exploiting Forward `Dual<8>` chunks or Enzyme JVP |
| Portfolio / Single-Swap Multi-Curve DV01 + Cross-Gamma | $50\text{–}250$ | $1$ | Forward-over-Reverse or `rateslib`-style `Dual2` with IFT |
| Monte Carlo Pathwise Exotic / XVA Sensitivities | $200\text{–}5,000$ | $1$ | Reverse-Mode Block-Arena Tape with `mark`/`rewind` + SIMD path batching |

---

## 2. Implicit Function Theorem (IFT) for Zero-Tape Calibration & Bootstrapping

One of the most important breakthroughs in quantitative finance AD is **Marc Henrard's (2011/2013) Implicit Function Theorem (IFT) method** ([`03_Henrard_2011_AAD_Calibration_Implicit_Function_Theorem.pdf`](file:///home/jakeshea/rust_quantlib/papers/03_Henrard_2011_AAD_Calibration_Implicit_Function_Theorem.pdf)).

### 2.1 Why Naively Taping Through Solvers Is an Anti-Pattern
Suppose market quotes $q \in \mathbb{R}^K$ determine calibrated model parameters or curve discount nodes $\theta(q) \in \mathbb{R}^K$ by solving the calibration equation:
$$F(\theta, q) = 0$$
via $30\text{–}100$ iterations of Brent, Newton-Raphson, or Levenberg-Marquardt.
If an AD tape stays active during calibration:
1. The tape records every unrolled root-finding iteration ($50\times$ tape bloat and cache pollution).
2. Line searches, step rejections, and stopping criteria introduce non-differentiable branches or inaccurate derivatives near convergence tolerances.

### 2.2 The IFT Adjoint Shortcut
By differentiating $F(\theta(q), q) = 0$ with respect to $q$ at the converged root $\theta^*$:
$$\frac{\partial F}{\partial \theta}\bigg|_{(\theta^*, q)} \frac{\mathrm{d}\theta}{\mathrm{d}q} + \frac{\partial F}{\partial q}\bigg|_{(\theta^*, q)} = 0 \implies \frac{\mathrm{d}\theta}{\mathrm{d}q} = -\left(\frac{\partial F}{\partial \theta}\right)^{-1} \frac{\partial F}{\partial q}$$

Consequently, the SOTA calibration/bootstrapping workflow is:
1. **Solve in Pure `f64`**: Run the curve bootstrapper or Levenberg-Marquardt optimizer using raw `f64` (zero AD overhead).
2. **Evaluate Residual Once at Convergence**: At the converged root $\theta^*$, evaluate $F(\theta^*, q)$ **once** with AD (or reuse the analytical/forward-AD Jacobian $\partial F / \partial \theta$ already factored by the Newton/LM solver).
3. **Propagate Adjoints via Linear Solve**: Given price adjoints $\bar{\theta} = (\partial V / \partial \theta)^\top$ from the downstream pricing engine, compute market-quote adjoints $\bar{q} = (\partial V / \partial q)^\top$ by solving a single $K \times K$ linear system:
   $$\left(\frac{\partial F}{\partial \theta}\right)^\top w = \bar{\theta}, \qquad \bar{q} = -\left(\frac{\partial F}{\partial q}\right)^\top w$$

### 2.3 Second-Order IFT Calibration Curvature (Exact Par-Rate Cross-Gamma)
Differentiating $P(q) = V(\theta(q))$ twice with respect to market quotes $q$ requires **both** the cashflow curvature and the calibration manifold curvature:
$$\frac{\mathrm{d}^2 P}{\mathrm{d} q^2} = J^\top \left(\nabla^2_\theta V(\theta^*) - \sum_{m=1}^K \lambda_m \nabla^2_\theta F_m(\theta^*, q)\right) J, \qquad J = \frac{\mathrm{d}\theta}{\mathrm{d}q} = -J_\theta^{-1} J_q, \quad \lambda = (J_\theta^\top)^{-1} \nabla_\theta V(\theta^*)$$
Omitting the second term $-\sum_{m=1}^K \lambda_m J^\top (\nabla^2_\theta F_m) J$ yields zero par-rate convexity for any payoff linear in discount factors $\theta$ (such as a zero-coupon bond $V(\theta) = \theta_k$). Similarly, for overdetermined least-squares calibration $\min_\theta \frac{1}{2}\|r(\theta, q)\|^2$ with $r_i(\theta, q) = P_i(\theta) - q_i$, differentiating the stationarity system $g(\theta, q) = J(\theta)^\top r(\theta, q) = 0$ yields the exact second-order optimality Jacobian:
$$\frac{\partial g}{\partial \theta} = J^\top J + \sum_{i=1}^m r_i(\theta^*) \nabla^2_\theta P_i(\theta^*), \qquad \frac{\partial g}{\partial q} = -J^\top$$

### 2.4 IFT Lifting of Free-Boundary American Options & Convex Risk-Parity Optimization
1. **Free-Boundary American Option Approximations (Barone–Adesi–Whaley 1987, Ju–Zhong 1999)**:
   - Analytic American approximations determine the critical early-exercise boundary $S^*$ via a nonlinear scalar equation $f(S^*; S, K, r, q, \sigma, T) = 0$, while at $S = S^*$ the smooth-pasting condition $\frac{\partial V}{\partial S^*}\big|_{S^*} = 0$ holds to first order (and for second-order or general parameter sensitivities $\theta \in \{K, r, q, \sigma, T\}$, `ift_lift_1d` lifts the converged `f64` boundary $S^*$ into `Dual<N>` / `Dual2<N>` / `TapeVar` in a single residual evaluation without recording Newton iterations).
2. **Risk Budgeting / Equal Risk Contribution (Spinu 2013, Griveau-Billion–Richard–Roncalli 2013)**:
   - Solving the strictly convex risk-budgeting optimization $\min_{y > 0} \frac{1}{2} y^\top \Sigma y - \sum_{i=1}^N b_i \ln y_i$ satisfies the first-order optimality system $F(y, \Sigma, b) = \Sigma y - b / y = 0$ with positive-definite diagonal-plus-covariance Hessian $H_y = \frac{\partial F}{\partial y} = \Sigma + \operatorname{diag}(b_i / y_i^{*2})$.
   - Solving for $y^*$ in pure `f64` via Cyclical Coordinate Descent + Newton–Cholesky and lifting $y^*$ into `T: Real` via a single linear solve $H_y \, \mathrm{d}y = -(\mathrm{d}\Sigma \, y^* - \mathrm{d}b / y^*)$ yields exact `Dual<N>` and `TapeVar` sensitivities of risk-parity portfolio weights $w^* = y^* / (\mathbf{1}^\top y^*)$ with respect to every covariance entry $\Sigma_{ij}$ and risk budget $b_i$ in $\mathcal{O}(N^3)$ time without taping through solver iterations.

> [!TIP]
> **Verified in Cloned Repos**: Both `attack68/rateslib` (`repos/rateslib/rust/src/curves/`) and `jmelo11/quantsupport` implement curve Jacobians to propagate sensitivities from discount factors/zero rates back to par swap/FRA quotes without taping the iterative solver loop.

---

## 3. Memory-Safe High-Performance AAD Tape Architecture in Rust

We inspected three generations of AAD tapes in our cloned repositories:
1. **Legacy `RustQuant_autodiff`** ([`repos/RustQuant/crates/RustQuant_autodiff/src/graph.rs`](file:///home/jakeshea/rust_quantlib/repos/RustQuant/crates/RustQuant_autodiff/src/graph.rs)): Single-threaded `RefCell<Vec<Vertex>>`, no `mark`/`rewind`, heap reallocations, disconnected from pricing traits.
2. **Savine's C++ `CompFinance`** ([`repos/CompFinance/`](file:///home/jakeshea/rust_quantlib/repos/CompFinance/)) and **`quantsupport`** ([`repos/quantsupport/src/ad/tape.rs`](file:///home/jakeshea/rust_quantlib/repos/quantsupport/src/ad/tape.rs)): Uses a slab allocator (`BlockList<TapeNode<T>>`) plus `book: Vec<NonNull<TapeNode<T>>>` with `set_mark()` and `rewind_to_mark()`.
   - *Critique of `quantsupport`'s Rust tape*: Storing `Vec<NonNull<TapeNode<T>>>` and mutating `adj` inside `TapeNode<T>` via `unsafe { (*ptr.as_ptr()).adj = ... }` (lines 117, 143 of `tape.rs`) introduces pointer indirection during the reverse sweep and violates Rust's aliasing rules unless wrapped in `UnsafeCell`.
3. **SOTA Safe Index-Based SoA Block-Arena Tape (Proposed for Rust)**:
   - **Split Primal Recording from Adjoint Storage**:
     - The forward pass appends compact 24-byte operations `(lhs_idx: u32, rhs_idx: u32, w_lhs: f64, w_rhs: f64)` into a pre-allocated segmented arena (`BumpVec<OpRecord>`).
     - Adjoints live in a separate contiguous `Vec<T>` (`adjoints: Vec<f64>` or `Vec<f64x8>`) indexed by `NodeIdx(u32)`.
   - **Zero `unsafe` Pointer Aliasing**: The backward pass is a tight bounds-checked (or slice-iterated) reverse loop over `ops[mark..end].iter().rev()` updating `adjoints[op.lhs_idx as usize] += adj * op.w_lhs`.
   - **Checkpoint / `mark` & `rewind_to_mark` for Monte Carlo**:
     - Record curve/surface construction on `[0..mark]`.
     - Inside each Rayon worker thread, clone/share the immutable pre-mark tape prefix, simulate paths on `[mark..end]`, propagate adjoints from path payoff down to `mark`, accumulate into thread-local `mark_adjoints[0..mark]`, and call `rewind_to_mark()` in $\mathcal{O}(1)$ by resetting the arena cursor.
     - After all paths finish, reduce `mark_adjoints` across threads and execute `propagate_mark_to_start()` **once**.
   - **SIMD-Vectorized Tape (`Tape<f64x8>`)**: Because the control flow and computational graph of 8 Monte Carlo paths are identical (modulo branchless `select` masks), a single tape node records **one** `[u32; 2]` index pair for 8 simultaneous paths (`w_lhs: f64x8`, `w_rhs: f64x8`), slashing tape memory and index-lookup overhead by $8\times$.

### 3.1 Discontinuity Smoothing ("Fuzzy AAD") for Barriers, Digitals, and Callables
Pathwise AAD fails when payoffs contain Heaviside step functions $\mathbf{1}_{\{S_T > K\}}$ (digitals, barriers, discrete autocall triggers) because the pathwise derivative is zero almost everywhere and a Dirac delta at the strike.
- **Savine's Static AST Smoothing** ([`repos/Scripting/`](file:///home/jakeshea/rust_quantlib/repos/Scripting/)): Instead of requiring users to manually rewrite payoffs, an AST visitor inspects the payoff script before compilation, detects discontinuous indicator predicates `x > k` that feed into continuous cashflows, and rewrites them into smooth sigmoidal or piecewise-linear call-spread approximations of width $2\varepsilon$:
  $$I_\varepsilon(x, k) = \mathrm{clamp}\!\left(\frac{x - k + \varepsilon}{2\varepsilon},\, 0,\, 1\right)$$
- **Vibrato / Likelihood Ratio Method (LRM) & Malliavin Weights**: For pure jump/digital payoffs where smoothing width $\varepsilon$ is sensitive, hybrid Pathwise + LRM (Giles & Glasserman 2005) or Malliavin integration by parts (implemented in `stochastic-rs/src/quant/pricing/malliavin.rs`) shifts the derivative onto the transition density score $\nabla_\theta \log p(S_T \mid S_{T-\Delta t})$.

---

## 4. Differential Machine Learning (DML) & Differential PCA

Huge & Savine (2020) ([`04_Huge_Savine_2020_Differential_Machine_Learning.pdf`](file:///home/jakeshea/rust_quantlib/papers/04_Huge_Savine_2020_Differential_Machine_Learning.pdf)) and Huge & Savine (2025) ([`05_Huge_Savine_2025_Axes_That_Matter_Differential_PCA.pdf`](file:///home/jakeshea/rust_quantlib/papers/05_Huge_Savine_2025_Axes_That_Matter_Differential_PCA.pdf)) demonstrated how AAD transforms statistical learning and dimensionality reduction in quantitative finance:

1. **Differential Machine Learning (DML)**:
   - Standard neural network surrogates trained solely on noisy Monte Carlo payoff labels $Y = \text{Payoff}(\omega)$ require $10^6+$ paths to learn smooth pricing functions and produce wildly inaccurate Greeks upon differentiation.
   - Because pathwise AAD extracts unbiased pathwise differentials $Z_i = \partial Y(\omega) / \partial X_i(\omega)$ alongside $Y(\omega)$ at $\le 4\times$ cost, a neural network (or polynomial regression) can be trained with a **combined Sobolev/differential loss**:
     $$\mathcal{L}(\theta) = \mathbb{E}\!\left[\left(\hat{V}_\theta(X) - Y\right)^2\right] + \lambda \sum_{j=1}^d w_j \,\mathbb{E}\!\left[\left(\frac{\partial \hat{V}_\theta(X)}{\partial X_j} - Z_j\right)^2\right]$$
   - Each training sample provides $1 + d$ constraints instead of $1$, achieving $100\times$ sample efficiency and exact Greek convergence with $1,000\text{–}4,000$ paths.
2. **Differential PCA ("Axes That Matter", 2025)**:
   - Traditional PCA decomposes the covariance of state variables $\mathrm{Cov}(X)$, which ignores whether high-variance directions actually impact the portfolio's value.
   - **Differential PCA** performs eigendecomposition on the **expected outer product of pathwise AAD gradients** $\mathbb{E}[\nabla_X Y \, (\nabla_X Y)^\top]$, directly discovering the low-dimensional linear subspace of risk factors that explain $>99.9\%$ of payoff variance and exposure.
   - In a modern Rust QuantLib, exposing `(payoff, pathwise_gradients)` as zero-copy Arrow/ndarray/Burn tensors makes the library a first-class engine for DML and Differential PCA out of the box.

---

## 5. High-Performance Numerical Solvers: PDEs & Monte Carlo

### 5.1 Multi-Dimensional ADI Finite-Difference Schemes
For 2D/3D stochastic volatility and hybrid models (Heston, Bates, Heston-Hull-White, Cheyette), explicit schemes are stability-limited ($\Delta t = \mathcal{O}(\Delta x^2)$) and fully implicit 2D/3D schemes require solving sparse banded systems with bandwidth $N_x N_v$.
As proven by Haentjens & in 't Hout (2012) ([`09_Haentjens_intHout_2012_ADI_Heston_Hull_White_PDE.pdf`](file:///home/jakeshea/rust_quantlib/papers/09_Haentjens_intHout_2012_ADI_Heston_Hull_White_PDE.pdf)) and Düring et al. (2015) ([`10_Duering_et_al_2015_High_Order_ADI_Stochastic_Vol.pdf`](file:///home/jakeshea/rust_quantlib/papers/10_Duering_et_al_2015_High_Order_ADI_Stochastic_Vol.pdf)):
- **Alternating Direction Implicit (ADI)** operator splitting decomposes the multi-dimensional differential operator $\mathcal{A} = A_0 + A_1 + A_2 + A_3$, where $A_0$ contains the mixed derivative $\rho \sigma v S \frac{\partial^2 u}{\partial S \partial v}$ (treated explicitly) and $A_1, A_2, A_3$ contain the unidirectional convection-diffusion operators along $S$, $v$, and $r$ (treated implicitly).
- **The Four Standard ADI Splitting Schemes**:
  1. **Douglas (Do)** — First-order unconditionally stable for $\theta \ge \frac{1}{2}$, simplest per step.
  2. **Craig-Sneyd (CS)** — Second-order accurate with an extra explicit mixed-derivative predictor/corrector step.
  3. **Modified Craig-Sneyd (MCS)** — Unconditionally stable in 2D/3D with arbitrary damping parameter $\theta$, damping high-frequency oscillations from non-smooth payoff kinks (Rannacher startup + MCS).
  4. **Hundsdorfer-Verwer (HV)** — Second-order accurate with strong stiff decay for high vol-of-vol $\xi$ when the Feller condition $2\kappa\theta > \xi^2$ is violated.
- **SIMD-Interleaved Thomas Algorithm**:
  - In a 2D grid of size $N_S \times N_v$, each implicit sweep along $S$ solves $N_v$ **independent** tridiagonal systems of size $N_S$.
  - While sequential Thomas algorithm along a single line has a loop-carried dependency ($u_i$ depends on $u_{i-1}$), **interleaving 8 grid lines across `f64x8` SIMD lanes** executes 8 tridiagonal solves simultaneously at $100\%$ AVX-512/NEON FMA utilization with zero matrix allocation!

### 5.2 Two-Pass Parallel American Monte Carlo (AMC) & XVA Exposure Cubes
In C++ `QuantLib` (and brokenly in `RustQuant`), Longstaff-Schwartz American Monte Carlo stores the entire $[N_{\text{paths}} \times N_{\text{steps}} \times d]$ state tensor in RAM during the forward pass, which requires gigabytes of memory and thrashes L3 cache.
Following Oancea et al. (2019) ([`11_Oancea_et_al_2019_Functional_Accelerating_American_MC.pdf`](file:///home/jakeshea/rust_quantlib/papers/11_Oancea_et_al_2019_Functional_Accelerating_American_MC.pdf)), Kenyon & Green (2014) ([`12_Kenyon_Green_2014_Efficient_XVA_Management_Regression.pdf`](file:///home/jakeshea/rust_quantlib/papers/12_Kenyon_Green_2014_Efficient_XVA_Management_Regression.pdf)), and Albanese & Crépey (2020) ([`13_Albanese_Crepey_2020_XVA_Analysis_Balance_Sheet.pdf`](file:///home/jakeshea/rust_quantlib/papers/13_Albanese_Crepey_2020_XVA_Analysis_Balance_Sheet.pdf)):
1. **Calibration stage, $N_{\text{cal}} \approx 4,096\text{–}16,384$ paths)**:
   - Simulate a moderate batch of paths (using scrambled Sobol or counter-based Philox RNG).
   - Perform backward induction using **QR / SVD-stabilized Chebyshev or Laguerre polynomial regression** (never normal equations $X^\top X$, which square the condition number $\kappa(X)^2$) to produce an immutable `ExercisePolicy` containing regression weights $\beta_t$ at each Bermudan exercise date $t_k$.
2. **Pricing and AAD stage, $N_{\text{price}} \approx 10^5\text{–}10^6$ paths, embarrassingly parallel)**:
   - Stream SIMD chunks of paths forward in time with $\mathcal{O}(1)$ memory per path, evaluating stopping decisions against the frozen `ExercisePolicy` $\beta_t$.
   - Because regression weights $\beta_t$ are near-optimal stopping boundaries, the **envelope theorem (Glasserman 2004)** guarantees that holding $\beta_t$ constant during the AAD pass yields unbiased first-order pathwise Greeks without differentiating through the SVD regression!

---

## 6. Compositional Financial Contracts & Payoff Bytecode VMs

Writing bespoke C++ `Instrument` + `PricingEngine` + `PathPricer` subclasses for every structured note (Autocallable Phoenix, Tarf, Snowball, Cliquet, CMS Steepener) is the single greatest source of code duplication in C++ `QuantLib`.

### 6.1 Algebraic Contract Combinators (Peyton Jones, Eber, Bahr)
Peyton Jones, Eber & Seward (2000) ([`15_PeytonJones_Eber_Seward_2000_Composing_Contracts.pdf`](file:///home/jakeshea/rust_quantlib/papers/15_PeytonJones_Eber_Seward_2000_Composing_Contracts.pdf)) and Bahr, Berthold & Elsman (2015) ([`16_Bahr_Berthold_Elsman_2015_Symbolic_Financial_Contracts.pdf`](file:///home/jakeshea/rust_quantlib/papers/16_Bahr_Berthold_Elsman_2015_Symbolic_Financial_Contracts.pdf)) proved that virtually any financial derivative can be expressed compositionally from a small set of algebraic primitives (`Zero`, `One(Currency)`, `Give(c)`, `And(c1, c2)`, `Or(c1, c2)`, `Cond(obs, c1, c2)`, `Scale(obs, c)`, `When(pred, c)`, `Anytime(pred, c)`).

### 6.2 Savine Scripting + Flat Bytecode VM
In Antoine Savine's `Scripting` architecture ([`repos/Scripting/`](file:///home/jakeshea/rust_quantlib/repos/Scripting/), partially ported in `quantsupport/src/scripting/`):
1. **Declarative Event Schedule**: A product is a chronologically ordered sequence of dated events:
   ```text
   EventDate(t_i):
     if Spot("SX5E") > 1.00 * InitialSpot:
         pay 1.05 * Notional at t_i + 2bd
         terminate
     else if Spot("SX5E") > 0.70 * InitialSpot:
         pay 0.05 * Notional at t_i + 2bd
   ```
2. **Static Analysis Passes over the AST**:
   - **Variable Indexing & Liveness**: Maps named variables to flat register slots `r0..rK`.
   - **Payoff Date & Underlying Extraction**: Automatically determines the exact simulation timeline $\{t_1, \dots, t_M\}$, required market observables, and discount factor maturities needed by the Monte Carlo or PDE engine.
   - **Fuzzy AAD Smoothing**: Rewrites discontinuous conditionals into smooth transitions of configurable width $\varepsilon$.
3. **Compilation to Flat SIMD Register Bytecode**:
   - Instead of walking a heap-allocated `Box<Expr>` tree on every Monte Carlo path (which causes branch mispredictions and pointer chasing), the AST compiles into a contiguous `Vec<Instruction>` operating on a small register file `[T; 16]` where `T = f64x8` or `AdVar<f64x8>`.
   - This enables **runtime-defined structured products** (loaded from JSON/DSL strings in Python or web UIs) to execute at near-native compiled Rust speed with full AAD Greeks and XVA cube integration.

---

## 7. Second-Order Optimality IFT, Active-Set KKT Differentiation & Multi-Level Dynamic Programming

### 7.1 Modular Implicit Differentiation of Nonlinear Least-Squares Calibrations (`papers/43`)
As formalized by Blondel et al. (2021) ([`43_Blondel_et_al_2021_Efficient_Modular_Implicit_Differentiation.pdf`](file:///home/jakeshea/rust_quantlib/papers/43_Blondel_et_al_2021_Efficient_Modular_Implicit_Differentiation.pdf)):
- When model parameters $\theta^* \in \mathbb{R}^P$ (such as Heston $(v_0, \kappa, \theta, \sigma_v, \rho)$, SABR $(\alpha, \rho, \nu)$, or arbitrage-free SVI $(a, b, \rho, m, \sigma)$) are calibrated via nonlinear least squares $\min_\theta \frac{1}{2}\sum_{i=1}^M w_i \left(\hat{\sigma}(K_i; \theta) - q_i\right)^2$, the stationary optimality condition is:
  $$G(\theta^*, q) = \nabla_\theta \mathcal{L}(\theta^*, q) = \sum_{i=1}^M w_i r_i(\theta^*, q) \nabla_\theta \hat{\sigma}(K_i; \theta^*) = 0.$$
- Differentiating $G(\theta^*(q), q) = 0$ with respect to market quotes $q \in \mathbb{R}^M$ yields the **exact second-order optimality IFT Jacobian**:
  $$\frac{\partial \theta^*}{\partial q} = -\left(\nabla_{\theta\theta}^2 \mathcal{L}(\theta^*, q)\right)^{-1} \nabla_{\theta q}^2 \mathcal{L}(\theta^*, q) = \left(J^\top W J + \sum_{i=1}^M w_i r_i(\theta^*, q)\,\nabla_{\theta\theta}^2 \hat{\sigma}(K_i; \theta^*)\right)^{-1} J^\top W.$$
- Crucially, omitting the second-order residual-Hessian term $\sum_i w_i r_i \nabla_{\theta\theta}^2 \hat{\sigma}(K_i; \theta^*)$ (the Gauss–Newton approximation) introduces $O(\|r\|)$ bias whenever market smiles cannot be fit with zero residual. Evaluating $\hat{\sigma}(K_i; \theta^*)$ with `Dual2<P>` supplies both $J$ and $\nabla_{\theta\theta}^2 \hat{\sigma}$ in a single pass, enabling exact quote Vegas/Vannas/Volgas in $O(1)$ without taping through Levenberg–Marquardt iterations.

### 7.2 Active-Set KKT Implicit Differentiation for Constrained Portfolio Optimization (`papers/43`, `papers/45`)
Following Agrawal et al. (2019) ([`45_Agrawal_et_al_2019_Differentiable_Convex_Optimization_Layers.pdf`](file:///home/jakeshea/rust_quantlib/papers/45_Agrawal_et_al_2019_Differentiable_Convex_Optimization_Layers.pdf)) and Blondel et al. (2021) ([`43_Blondel_et_al_2021_Efficient_Modular_Implicit_Differentiation.pdf`](file:///home/jakeshea/rust_quantlib/papers/43_Blondel_et_al_2021_Efficient_Modular_Implicit_Differentiation.pdf)):
- For box-constrained Markowitz mean–variance allocation $\min_w \frac{1}{2}\gamma w^\top \Sigma w - \mu^\top w$ subject to $\mathbf{1}^\top w = 1$ and $l \le w \le u$, unrolling projected-gradient steps onto an AD tape is slow and numerically brittle.
- Instead, solving the primal QP in `f64` identifies the exact active bound sets $\mathcal{B}_L = \{i : w_i^* = l_i\}$, $\mathcal{B}_U = \{i : w_i^* = u_i\}$, and strictly interior free set $\mathcal{F} = \{i : l_i < w_i^* < u_i\}$.
- On $\mathcal{B}_L \cup \mathcal{B}_U$, $w_i$ equals the bound $l_i$ or $u_i$ in `T: Real` directly, while on $\mathcal{F}$, $w_{\mathcal{F}}$ solves the reduced equality KKT linear system in `T: Real`:
  $$\begin{bmatrix} \gamma \Sigma_{\mathcal{F},\mathcal{F}} & \mathbf{1} \\ \mathbf{1}^\top & 0 \end{bmatrix} \begin{bmatrix} w_{\mathcal{F}} \\ \nu \end{bmatrix} = \begin{bmatrix} \mu_{\mathcal{F}} - \gamma \Sigma_{\mathcal{F},\mathcal{B}} w_{\mathcal{B}} \\ 1 - \mathbf{1}^\top w_{\mathcal{B}} \end{bmatrix},$$
  preserving both first-order (`Dual<N>`, `TapeVar`) and second-order (`Dual2<N>`) derivatives with respect to $(\Sigma, \mu, \gamma, l, u)$ in a single linear solve.

### 7.3 Analytical Arrow–Debreu Lattices & Multi-Level Dynamic Programming (`papers/49`, `papers/50`, `papers/51`)
- In short-rate Bermudan swaption trees (`HullWhiteBermudanTreeEngine`) and commodity swing/storage "Forest of Trees" engines (`Warin 2012` [`49_Warin_2012_Hedging_Swing_Contract_on_Gas_Market.pdf`](file:///home/jakeshea/rust_quantlib/papers/49_Warin_2012_Hedging_Swing_Contract_on_Gas_Market.pdf), `Bäuerle & Rieß 2015` [`51_Baeuerle_Riess_2015_Gas_Storage_Valuation_Multinomial_Tree_LSM.pdf`](file:///home/jakeshea/rust_quantlib/papers/51_Baeuerle_Riess_2015_Gas_Storage_Valuation_Multinomial_Tree_LSM.pdf)), forward Arrow–Debreu state-price shifts $\alpha_m = \frac{1}{\Delta t_m}\ln\!\left(\frac{\sum_j Q_{m,j} e^{-x_j \Delta t_m}}{P(0, t_{m+1})}\right)$ and multiplicative futures normalizations $S_{m,j} = F(0, t_m)\frac{e^{x_j}}{\sum_l P_{m,l} e^{x_l}}$ are **closed-form in `T: Real`**, eliminating iterative root-finding inside the tree and propagating exact curve/surface sensitivities through backward dynamic programming.

---

## 8. Next-Generation Institutional Frontiers: SLV Particle Calibration, Rough Bergomi Hybrid Volterra, Multi-Asset Autocallable Cega, FMM/LMM, FRTB, and Transient-Impact Execution

### 8.1 Guyon–Henry-Labordère McKean–Vlasov SLV Particle Calibration (`papers/53`, `papers/54`)
Following Guyon & Henry-Labordère (2012), Cozma, Mariapragassam & Reisinger (2018) ([`53_Cozma_Mariapragassam_Reisinger_2018_SLV_Particle_Method.pdf`](file:///home/jakeshea/rust_quantlib/papers/53_Cozma_Mariapragassam_Reisinger_2018_SLV_Particle_Method.pdf)), and Reisinger & Tsianni (2025) ([`54_Reisinger_Tsianni_2025_Particle_System_Heston_SLV.pdf`](file:///home/jakeshea/rust_quantlib/papers/54_Reisinger_Tsianni_2025_Particle_System_Heston_SLV.pdf)):
- By Gyöngy's (1986) theorem, a Heston Stochastic Local Volatility (SLV) model $dS_t/S_t = (r-q)dt + L(t, S_t)\sqrt{V_t}dW_t^S$ matches the vanilla market implied volatility surface $\sigma_{\text{Dup}}(t, S)$ if and only if the **leverage surface** satisfies the conditional expectation identity:
  $$L(t, S) = \frac{\sigma_{\text{Dup}}(t, S)}{\sqrt{\mathbb{E}[V_t \mid S_t = S]}}.$$
- Evaluating $\widehat{\mathbb{E}}[V_{t_m} \mid S_{t_m} = x_k]$ via a Nadaraya–Watson kernel estimator over $N_p$ interacting particles in `T: Real` (holding only the Silverman bandwidth $h_m = \kappa_h \hat{\sigma}_{X_m} N_p^{-1/5}$ and spatial grid knots $x_k$ fixed in `f64`) plus empirical martingale projection enables pathwise `Dual` and `TapeVar` sensitivities through the entire McKean–Vlasov calibration and exotic payoff simulation.

### 8.2 Turbocharged Hybrid Volterra Simulation for Rough Bergomi (`papers/52`)
Following Bayer, Friz & Gatheral (2016) and McCrickerd & Pakkanen (2018) ([`52_McCrickerd_Pakkanen_2018_Turbocharging_MC_Rough_Bergomi.pdf`](file:///home/jakeshea/rust_quantlib/papers/52_McCrickerd_Pakkanen_2018_Turbocharging_MC_Rough_Bergomi.pdf)):
- In the rough Bergomi (`rBergomi`) model $V_t = \xi_0(t)\exp\!\left(\eta\widetilde{W}_t^H - \frac{1}{2}\eta^2 t^{2H}\right)$ with Hurst parameter $H \in (0, \frac{1}{2})$, the Riemann–Liouville Volterra integral $\widetilde{W}_{t_k}^H = \sqrt{2H}\int_0^{t_k}(t_k - s)^{H-1/2}dW_s^{(1)}$ splits into an exact bivariate Gaussian head $(\Delta W_k^{(1)}, \widehat{W}_k)$ and optimal root-mean-square power-law weights $b_m = \left(\frac{m^{H+1/2} - (m-1)^{H+1/2}}{H + 1/2}\right)^{\frac{1}{H-1/2}}$ for the regular tail.
- Conditioning on $\mathcal{F}_T^{W^{(1)}}$ (Romano–Touzi projection) turns the terminal payoff into an analytic Black–Scholes evaluation on $(S_T^{\parallel}, \sqrt{(1-\rho^2)\int_0^T V_t dt / T})$, cutting Monte Carlo variance by $20\times$ and yielding smooth pathwise Greeks in `T: Real`.

### 8.3 Multi-Asset Worst-Of Autocallables with Fuzzy-AAD & Capriotti–Giles Adjoint Cholesky Cega (`papers/02`, `papers/16`)
- Combining a shift-invariant smooth minimum $W_m = \text{smin}_\tau(S_{1,t_m}/S_{1,0}, \dots, S_{d,t_m}/S_{d,0})$, $C^2$ quintic Fuzzy-AAD barrier/coupon gates $H_\epsilon^{\text{C2}}(W_m - B_m)$, and the Capriotti–Giles (2010) ([`02_Capriotti_Giles_2010_Fast_Correlation_Greeks_AAD.pdf`](file:///home/jakeshea/rust_quantlib/papers/02_Capriotti_Giles_2010_Fast_Correlation_Greeks_AAD.pdf)) adjoint Cholesky pullback $\bar{\Sigma} = \text{Pullback}_L(\bar{L})$ yields all $d$ spot Deltas, $d$ Vegas, and the full $d(d-1)/2$ cross-asset correlation sensitivity matrix ($\text{Cega}_{ij} = \partial V / \partial \rho_{ij}$) of Phoenix/Snowball Worst-Of Autocallables in a single reverse pass.

### 8.4 Generalized Forward Market Model (FMM / LMM) & 2-Factor `G2++` Bermudans (`papers/50`, `papers/58`)
- Extending the Brace–Gatarek–Musiela (1997) LMM to the Lyashenko–Mercurio (2019) / López-Salas et al. (2024) ([`58_LopezSalas_2024_PDEs_Generalized_Forward_Market_Model.pdf`](file:///home/jakeshea/rust_quantlib/papers/58_LopezSalas_2024_PDEs_Generalized_Forward_Market_Model.pdf)) Generalized Forward Market Model incorporates linear volatility decay $\gamma_i(t) = (T_{i+1}-t)/(T_{i+1}-T_i)$ inside each backward-looking RFR accrual period $[T_i, T_{i+1}]$, Hunter–Jäckel–Joshi (2001) predictor–corrector terminal-measure drift stepping, and pathwise AAD for Target Redemption Notes (`TARNs`) and Callable Range Accruals alongside 2D `G2++` Bermudan swaption quadrature/lattice engines ([`50_Yamakami_Takeuchi_2022_Bermudan_Swaption_Two_Factor_Hull_White.pdf`](file:///home/jakeshea/rust_quantlib/papers/50_Yamakami_Takeuchi_2022_Bermudan_Swaption_Two_Factor_Hull_White.pdf)).

### 8.5 Basel III / IV FRTB (MAR20–MAR33) & Arbitrage-Free Bilateral Funding XVA (`papers/56`, `papers/59`)
- Under FRTB ([`59_Tasche_2018_Capital_Allocation_FRTB_Expected_Shortfall.pdf`](file:///home/jakeshea/rust_quantlib/papers/59_Tasche_2018_Capital_Allocation_FRTB_Expected_Shortfall.pdf)), the Standardized Approach (Sensitivities-Based Method Delta/Vega/Curvature under Low/Medium/High correlation scenarios + DRC + RRAO) and Internal Models Approach (Liquidity-Horizon-Scaled Expected Shortfall $ES = \sqrt{ES_1^2 + \sum_{j \ge 2}(ES_j\sqrt{(LH_j - LH_{j-1})/10})^2}$) are formulated in `T: Real` so Euler AAD attributes desk and trade-level market-risk capital in one pass.
- Under asymmetric borrowing and lending rates ($r_f^+ \ne r_f^-$) and collateral rates ($r_c$), Bichuch, Capponi & Sturm (2018) ([`56_Bichuch_Capponi_Sturm_2018_Arbitrage_Free_XVA_Funding.pdf`](file:///home/jakeshea/rust_quantlib/papers/56_Bichuch_Capponi_Sturm_2018_Arbitrage_Free_XVA_Funding.pdf)) and Burgard–Kjaer (2011) formulate total bilateral XVA ($\text{CVA} - \text{DVA} + \text{FCA} - \text{FBA}$) as a contractive backward fixed-point equation whose converged solution lifts to `T: Real` via IFT.

### 8.6 Smooth Rockafellar–Uryasev Mean–CVaR & Obizhaeva–Wang Transient-Impact Execution (`papers/55`, `papers/57`)
- Smoothing the Rockafellar–Uryasev (2000) ([`57_Rockafellar_Uryasev_2000_Optimization_of_CVaR.pdf`](file:///home/jakeshea/rust_quantlib/papers/57_Rockafellar_Uryasev_2000_Optimization_of_CVaR.pdf)) excess-loss hinge $[u]_+ \approx \frac{1}{2}(u + \sqrt{u^2 + \varepsilon^2})$ enables active-set KKT implicit differentiation of constrained Mean–CVaR portfolios without unrolling LP/QP iterations onto the AD tape.
- Coupled with Almgren–Chriss (2001) and Obizhaeva–Wang (2013) / Alfonsi–Fruth–Schied (2010) ([`55_Alfonsi_Fruth_Schied_2010_Optimal_Execution_Transient_Impact.pdf`](file:///home/jakeshea/rust_quantlib/papers/55_Alfonsi_Fruth_Schied_2010_Optimal_Execution_Transient_Impact.pdf)) transient limit-order-book resilience $G(t-s) = \theta_0 e^{-\rho(t-s)}$, portfolio rebalancing and optimal execution schedules become end-to-end differentiable in `T: Real`.

---

## 9. Phase 13 Institutional Frontiers: Monotone Convex & Meeting-Date OIS Curves, Global eSSVI & VIX Derivatives, Parisian/Parasian Exotics, Andersen–Broadie Dual AMC, Coterminal Cascade / Range Accruals / LPI Inflation, $k$-th-to-Default & Cashflow-Aware MPoR XVA, and Ledoit–Wolf / Robust Markowitz

### 9.1 Hagan–West (2006) Monotone Convex Splines & Central-Bank Meeting-Date Step-Jump OIS Curves
- Following Hagan & West (2006), interpolating discrete forwards $f_i^d = (r_i t_i - r_{i-1} t_{i-1}) / \Delta_i$ via the four-region monotone convex piecewise quadratic $g(x) = f(t_{i-1} + \Delta_i x) - f_i^d$ guarantees $\int_0^1 g(x)\,dx = 0$ (exact discount factor recovery at all pillars), non-negative instantaneous forwards $f(t) \ge 0$, and localized par-rate DV01s without the Gibbs oscillations of cubic zero-rate splines.
- Over the short-end central-bank policy horizon $[0, \tau_M]$, `MeetingDateJumpOisCurve<T>` models piecewise-constant overnight forwards $r_{\text{ON}}(t) = r_0 + \sum_{k=1}^m \Delta r_k$ across scheduled FOMC/ECB/BoE meeting dates $\tau_m$, with an analytic lower-triangular Jacobian $J_{\Delta r}$ and exact Henrard IFT Meeting-Date Jump DV01s (`CalibratedMeetingJumpOisCurve`).

### 9.2 Global Arbitrage-Free eSSVI Calibration & CBOE VIX / Heston–Bates VIX Derivatives (`papers/60`, `papers/61`, `papers/62`)
- Following Hendriks & Martini (2017) and Mingone (2022) ([`60_Mingone_2022_No_Arbitrage_Global_Parametrization_eSSVI.pdf`](file:///home/jakeshea/rust_quantlib/papers/60_Mingone_2022_No_Arbitrage_Global_Parametrization_eSSVI.pdf)), `AdEssviSurface<T>` and `CalibratedEssviSurface` parameterize multi-expiry total variance slices $w(k; \theta_i, \rho_i, \psi_i)$ subject to butterfly bounds $\psi_i(1 + |\rho_i|) \le \min(4, \sqrt{4\theta_i(1+|\rho_i|)})$ and cross-maturity calendar monotonicity $\theta_{i+1} \ge \theta_i$, $\psi_{i+1} \ge \psi_i$, $|\rho_{i+1}\psi_{i+1} - \rho_i\psi_i| \le \psi_{i+1} - \psi_i$, with `Dual2` second-order optimality IFT quote Vegas.
- Following Lian & Zhu (2013), Jacquier, Martini & Muguruza (2018) ([`61_Jacquier_Martini_Muguruza_2018_VIX_Futures_Rough_Bergomi.pdf`](file:///home/jakeshea/rust_quantlib/papers/61_Jacquier_Martini_Muguruza_2018_VIX_Futures_Rough_Bergomi.pdf)), and Horvath, Jacquier & Tankov (2020) ([`62_Horvath_Jacquier_Tankov_2020_Volatility_Options_Rough_Vol.pdf`](file:///home/jakeshea/rust_quantlib/papers/62_Horvath_Jacquier_Tankov_2020_Volatility_Options_Rough_Vol.pdf)), `ql-models::vix` combines model-free CBOE VIX strip replication with exact non-central chi-square Poisson-mixture quadrature for Heston/Bates VIX futures and options (`heston_bates_vix_future_and_option_price<T>`).

### 9.3 Parisian & Cumulative Parasian Barrier Options & Andersen–Broadie Primal–Dual AMC (`papers/63`, `papers/64`)
- Following Chesney, Jeanblanc-Picqué & Yor (1997) and Markov excursion-time stopping frameworks ([`63_General_Approach_Parisian_Stopping_Times_Markov_Processes_2021.pdf`](file:///home/jakeshea/rust_quantlib/papers/63_General_Approach_Parisian_Stopping_Times_Markov_Processes_2021.pdf)), `ql-exotics::parisian` augments the lattice state space with the integer excursion clock $m \in \{0, \dots, M\}$ (resetting on exit from $\mathcal{B}$ for consecutive Parisian options, or pausing for cumulative Parasian options) alongside smoothed pathwise AAD Monte Carlo.
- Following Andersen & Broadie (2004) and Rogers (2002) ([`64_Primal_Dual_Optimal_Stopping_Signatures_2023.pdf`](file:///home/jakeshea/rust_quantlib/papers/64_Primal_Dual_Optimal_Stopping_Signatures_2023.pdf)), `ql-mc::andersen_broadie` complements the Longstaff–Schwartz primal lower bound $L_0 = \mathbb{E}[Z_{\hat{\tau}}]$ with the Doob–Meyer martingale dual upper bound $U_0 = \mathbb{E}[\max_{1 \le k \le K}(Z_k - \pi_k)]$, bracketing the true Bermudan/American price in $[L_0, U_0]$ with pathwise AAD.

### 9.4 Coterminal Swaption Volatility Cascade IFT, Dual-Curve Callable Range Accruals & Jarrow–Yildirim LPI Inflation Swaps
- In `ql-rates`, `calibrate_cheyette_coterminal_ad` sequentially bootstraps piecewise-constant short-rate volatilities $\sigma_1, \dots, \sigma_M$ to coterminal European swaptions with exact IFT lifting, while `price_callable_range_accrual_hw1f` values dual-curve (OIS discount / IBOR-RFR forward) Bermudan callable range accrual notes with $C^1$ corridor smoothing.
- In `ql-inflation`, `lpi_inflation_swap_npv` values UK/European Limited Price Indexation (`LPI(floor, cap)`) year-on-year and compounded zero-coupon inflation swaps with Jarrow–Yildirim (2003) / Mercurio (2005) convexity adjustments and collar optionlet decomposition.

### 9.5 Heterogeneous $k$-th-to-Default Baskets & Andersen–Pykhtin–Sokol Cashflow-Aware MPoR / Contingent CDS (`papers/65`)
- In `ql-credit`, `price_kth_to_default_basket` uses the $O(N^2)$ conditional default-count recursion and cavity deconvolution under a one-factor Gaussian copula to value heterogeneous $k$-th-to-default protection and premium legs in generic `T: Real`.
- In `ql-xva`, following Andersen, Pykhtin & Sokol (2017) and Silotto, Scaringi & Bianchetti (2021) ([`65_Silotto_Scaringi_Bianchetti_2021_XVA_Model_Risk_MPoR.pdf`](file:///home/jakeshea/rust_quantlib/papers/65_Silotto_Scaringi_Bianchetti_2021_XVA_Model_Risk_MPoR.pdf)), `compute_mpor_xva` captures exposure spikes caused by unpaid counterparty cashflows and lagged variation margin over the Margin Period of Risk $[t - \delta_{\text{MPoR}}, t]$, alongside `price_contingent_cds` (`CCDS`).

### 9.6 Ledoit–Wolf Optimal Covariance Shrinkage & SOCP Robust Ellipsoidal Markowitz (`papers/66`)
- Following Ledoit & Wolf (2004, 2012) ([`66_Ledoit_Wolf_2012_Nonlinear_Shrinkage_Covariance_Matrices.pdf`](file:///home/jakeshea/rust_quantlib/papers/66_Ledoit_Wolf_2012_Nonlinear_Shrinkage_Covariance_Matrices.pdf)), `ledoit_wolf_shrinkage<T>` computes the finite-sample optimal Frobenius shrinkage intensity $\delta^* \in [0, 1]$ toward scaled-identity or constant-correlation targets in generic `T: Real`.
- Following Goldfarb & Iyengar (2003), `try_robust_markowitz_weights<T>` solves the worst-case mean-variance allocation under ellipsoidal expected-return uncertainty $\mathcal{U}_\kappa = \{\mu : (\mu - \hat{\mu})^\top \Omega_\mu^{-1}(\mu - \hat{\mu}) \le \kappa^2\}$ via iteratively reweighted active-set QP in `f64` followed by an exact bordered-Hessian KKT implicit differentiation solve in `T: Real`.
