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
