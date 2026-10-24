# Codex Master Prompt — Phase 12: Next-Generation Institutional Quant Frontiers

## Context & Architectural Tenets
Phase 11 converged at commit `94e755e` with `0 CRITICAL / 0 HIGH / 0 MEDIUM` findings across all 21 crates and 22 benchmark suites. We have expanded `papers/` to **59 primary research papers** (`52_McCrickerd_Pakkanen_2018_Turbocharging_MC_Rough_Bergomi.pdf` through `59_Tasche_2018_Capital_Allocation_FRTB_Expected_Shortfall.pdf`) and updated `knowledge_base/00_INDEX_AND_VERDICT.md`, `04_SOTA_QUANT_COMPUTING_PARADIGMS.md`, and `06_LITERATURE_AND_REFERENCES.md`.

Your mission in **Phase 12** is to implement **Deliverables 0 through 6** below with zero shortcuts, zero stubs, full `T: Real` (`f64`, `Dual<N>`, `Dual2<N>`, `TapeVar<f64>`) automatic differentiation and implicit function theorem (IFT) fidelity, comprehensive mathematical integration tests, **Benchmarks 23 & 24** in `ql-benches`, and **Workstation Tabs 21 & 22** in `ql-terminal`.

> **Critical Rules**:
> 1. **Preserve `papers/52_*.pdf`..`59_*.pdf` and `knowledge_base/*.md`**: Do NOT delete or `git clean` the newly added PDFs in `papers/` or markdown files in `knowledge_base/`. Include them in your commit.
> 2. **No Stubs or Constant/Identity Shortcuts**: Every numerical engine, calibration routine, benchmark closure, and workstation tab must execute the genuine mathematical algorithm described below.
> 3. **Strict `T: Real` AD Preservation**: Never strip AD tangents via `T::from_f64(x.val())` on state variables, payoffs, discount factors, or KKT/IFT outputs. Only discrete grid/bin indices, loop counters, and kernel bandwidths may inspect `.val()`.

---

## Deliverable 0: Foundational AD, IFT, Cholesky-Cega, & SIMM/Storage Hardening

1. **`crates/ql-ad/src/scalar.rs` ($C^2$ Smooth Heaviside & Shift-Invariant Smooth Minimum)**:
   - Add a default method on `Real`:
     ```rust
     #[inline]
     fn smooth_heaviside_c2(self, threshold: Self, half_width: f64) -> Self {
         let hw = half_width.abs().max(1e-12);
         let u = (self - threshold) / Self::from_f64(hw);
         if u.val() <= -1.0 {
             Self::zero()
         } else if u.val() >= 1.0 {
             Self::one()
         } else {
             let u2 = u * u;
             let u3 = u2 * u;
             let u5 = u3 * u2;
             Self::from_f64(0.5)
                 + Self::from_f64(15.0 / 16.0) * u
                 - Self::from_f64(5.0 / 8.0) * u3
                 + Self::from_f64(3.0 / 16.0) * u5
         }
     }
     ```
   - Add a free helper function in `crates/ql-ad/src/scalar.rs` (and re-export from `crates/ql-ad/src/lib.rs`):
     ```rust
     pub fn smooth_min_slice<T: Real>(xs: &[T], temperature: f64) -> T
     ```
     which asserts `!xs.is_empty()`, returns the exact `xs.iter().copied().reduce(|a, b| a.min_val(b)).unwrap()` when `temperature <= 1e-12 || xs.len() == 1`, and otherwise computes the shift-invariant log-sum-exp smooth minimum:
     $$\text{smin}_\tau(x_1, \dots, x_d) = m - \tau \ln\!\left(\frac{1}{d}\sum_{i=1}^d \exp\!\left(-\frac{x_i - m}{\tau}\right)\right)$$
     where $m = \min_i x_i$ (as the `T: Real` minimum element so shift invariance and `T: Real` derivatives at $\tau \to 0^+$ are exact).

2. **`crates/ql-ad/src/ift.rs` (SVD Fallback in `IftSystemSensitivity`)**:
   - In `IftSystemSensitivity::vjp`, `jvp`, and `jacobian_dtheta_dq` (`crates/ql-ad/src/ift.rs:95-117`), replace `.expect("singular IFT Jacobian")` with an SVD pseudo-inverse fallback (`self.j_theta.clone().svd(true, true).solve(&rhs, 1e-12)` / `self.j_theta.transpose().svd(true, true).solve(&rhs, 1e-12)`) when `lu.solve` returns `None`.

3. **`crates/ql-dml/src/cholesky_aad.rs` (Generic `T: Real` Cholesky & Symmetrized Zero-Diagonal Correlation Cega)**:
   - Add:
     ```rust
     pub fn cholesky_lower_generic<T: Real>(c: &[Vec<T>]) -> Result<Vec<Vec<T>>, String>
     ```
     implementing lower-triangular Cholesky factorization in generic `T: Real` (guarding `sum.max_val(T::from_f64(1e-16)).sqrt()` on the diagonal).
   - Add:
     ```rust
     pub fn cholesky_corr_cega_adjoint(
         c: &[Vec<f64>],
         l_bar: &[Vec<f64>],
     ) -> Result<Vec<Vec<f64>>, String>
     ```
     which runs `cholesky_adjoint(c, l_bar)` to obtain the symmetric sensitivity matrix $\bar{C}$, sets the diagonal $\bar{C}_{ii} = 0.0$ (since $\rho_{ii} \equiv 1$ is fixed), and sets off-diagonal entries $(i \ne j)$ to $\bar{C}_{ij} + \bar{C}_{ji} = 2\bar{C}_{ij}$ (or keeps both upper/lower entries equal to the total derivative with respect to bumping the symmetric pair $(\rho_{ij}, \rho_{ji})$ simultaneously: document clearly that `cega[i][j]` for $i \ne j$ is the total sensitivity $\frac{d V}{d \rho_{ij}} = \frac{\partial V}{\partial C_{ij}} + \frac{\partial V}{\partial C_{ji}}$ when both $(i,j)$ and $(j,i)$ are perturbed together).

4. **`crates/ql-xva/src/simm.rs` & `crates/ql-commodities/src/swing_storage.rs`**:
   - In `crates/ql-xva/src/simm.rs:120, 141`, change `k2.sqrt()` to `k2.max_val(T::zero()).sqrt()` and `total.sqrt()` to `total.max_val(T::zero()).sqrt()`.
   - In `crates/ql-commodities/src/swing_storage.rs:304-314` (`price_gas_storage`), hoist the computation of `let mut cont = vec![T::zero(); g + 1];` out of the `for gi in 0..=g` loop (computing `cont[gp]` once per `(m, j)` node instead of `g + 1` times).

---

## Deliverable 1: Guyon–Henry-Labordère SLV Particle Leverage Calibration & Rough Bergomi Hybrid Volterra Engine (`crates/ql-models/src/slv.rs`, `crates/ql-models/src/rough.rs`)

### 1.1 Stochastic Local Volatility (SLV) Particle Method (`crates/ql-models/src/slv.rs`, re-exported in `crates/ql-models/src/lib.rs`)
Implement the Guyon & Henry-Labordère (2012) / Cozma, Mariapragassam & Reisinger (2018) (`papers/53`) / Reisinger & Tsianni (2025) (`papers/54`) McKean–Vlasov particle method for Heston-SLV:
$$\frac{dS_t}{S_t} = (r - q)\,dt + L(t, S_t)\sqrt{V_t}\,dW_t^S, \qquad dV_t = \kappa(\theta - V_t)\,dt + \eta_{\text{mix}}\sigma_v\sqrt{V_t}\,dW_t^V, \qquad d\langle W^S, W^V\rangle_t = \rho\,dt$$
where $\eta_{\text{mix}} \in [0, 1]$ is the mixing fraction (`0.0` = pure Dupire local volatility, `1.0` = full Heston-SLV).

- **Types & Signatures**:
  ```rust
  #[derive(Clone, Debug)]
  pub struct SlvLeverageSlice<T: Real> {
      pub time: f64,
      pub log_moneyness_grid: Vec<f64>,
      pub leverage: Vec<T>,
  }

  #[derive(Clone, Debug)]
  pub struct CalibratedSlvModel<T: Real> {
      pub heston: HestonParams<T>,
      pub mixing_fraction: T,
      pub slices: Vec<SlvLeverageSlice<T>>,
  }

  impl<T: Real> CalibratedSlvModel<T> {
      /// Interpolate L(t_step, x) where x = ln(S / F(0, t)).
      pub fn leverage_at(&self, step_idx: usize, log_moneyness: T) -> T;
  }

  #[allow(clippy::too_many_arguments)]
  pub fn calibrate_slv_particle_leverage<T, LV>(
      spot: T,
      rate: T,
      div: T,
      expiry: f64,
      heston: &HestonParams<T>,
      mixing_fraction: T,
      dupire_local_vol: LV,
      n_steps: usize,
      n_particles: usize,
      n_grid: usize,
      seed: u64,
  ) -> CalibratedSlvModel<T>
  where
      T: Real,
      LV: Fn(f64, T) -> T;

  #[allow(clippy::too_many_arguments)]
  pub fn price_slv_european_mc<T: Real>(
      model: &CalibratedSlvModel<T>,
      spot: T,
      strike: T,
      rate: T,
      div: T,
      expiry: f64,
      n_paths: usize,
      is_call: bool,
      seed: u64,
  ) -> T;
  ```
- **Algorithm & Mathematical Invariants**:
  1. At $t_0 = 0$, all $N_p$ particles start at $X_0^{(i)} = 0$, $V_0^{(i)} = v_0$, and $L(0, x_k) = \sigma_{\text{Dup}}(0, S_0 e^{x_k}) / \sqrt{\max(v_0, 10^{-8})}$.
  2. At each step $m = 0, \dots, n_{\text{steps}} - 1$ ($\Delta t = T / n_{\text{steps}}$, $t_m = m \Delta t$):
     - For $m \ge 1$, on the log-moneyness grid $x_k \in [-3.5\bar{\sigma}\sqrt{t_m}, +3.5\bar{\sigma}\sqrt{t_m}]$ ($k = 0, \dots, n_{\text{grid}}-1$), estimate $\widehat{\mathbb{E}}[V_{t_m} \mid X_{t_m} = x_k]$ using Nadaraya–Watson kernel regression in `T: Real`:
       $$\widehat{v}_{\text{cond}}(t_m, x_k) = \frac{\sum_{i=1}^{N_p} V_m^{(i)} K\!\left(\frac{X_m^{(i)} - x_k}{h_m}\right) + \varepsilon_{\text{reg}} \bar{V}_m}{\sum_{i=1}^{N_p} K\!\left(\frac{X_m^{(i)} - x_k}{h_m}\right) + \varepsilon_{\text{reg}}}, \qquad \bar{V}_m = \frac{1}{N_p}\sum_{i=1}^{N_p} V_m^{(i)}$$
       with Gaussian or quartic biweight kernel $K(u)$ evaluated in `T: Real`, Silverman bandwidth $h_m = \max(1.5 \hat{\sigma}_{X_m} N_p^{-1/5}, 0.02)$ computed from `.val()`, and $\varepsilon_{\text{reg}} = 10^{-4}$, then set:
       $$L_m(x_k) = \text{clamp}\!\left(\frac{\sigma_{\text{Dup}}(t_m, F_{0,t_m} e^{x_k})}{\sqrt{\widehat{v}_{\text{cond}}(t_m, x_k).max\_val(10^{-8})}},\; 0.1,\; 5.0\right).$$
     - Propagate each particle $(X_m^{(i)}, V_m^{(i)})$ to $t_{m+1}$ using correlated antithetic normals $(Z_1^{(i)}, Z_2^{(i)})$, full-truncation / reflected Milstein step for $V_{m+1}^{(i)}$ with effective vol-of-vol $\xi_{\text{eff}} = \eta_{\text{mix}} \sigma_v$, and log-Euler step for $\widetilde{X}_{m+1}^{(i)}$.
     - Apply **exact empirical martingale projection** in `T: Real`:
       $$X_{m+1}^{(i)} = \widetilde{X}_{m+1}^{(i)} - \ln\!\left(\frac{1}{N_p}\sum_{j=1}^{N_p} \exp\!\left(\widetilde{X}_{m+1}^{(j)}\right)\right)$$
       so that $\frac{1}{N_p}\sum_{i=1}^{N_p} S_{m+1}^{(i)} = S_0 e^{(r - q) t_{m+1}}$ holds to machine precision ($< 10^{-12}$) at every time step!

### 1.2 Rough Bergomi (rBergomi) Hybrid Volterra Simulation (`crates/ql-models/src/rough.rs`)
Implement the Bayer, Friz & Gatheral (2016) / McCrickerd & Pakkanen (2018) (`papers/52`) hybrid Volterra scheme with Romano–Touzi conditional Black–Scholes variance reduction:
- **Types & Signatures**:
  ```rust
  #[derive(Clone, Debug, PartialEq)]
  pub struct RoughBergomiParams<T: Real> {
      pub xi0: T,
      pub eta: T,
      pub rho: T,
      pub hurst: f64,
  }

  #[allow(clippy::too_many_arguments)]
  pub fn rbergomi_hybrid_mc_price<T: Real>(
      spot: T,
      strike: T,
      expiry: T,
      rate: T,
      div: T,
      params: &RoughBergomiParams<T>,
      n_steps: usize,
      n_paths: usize,
      is_call: bool,
      seed: u64,
  ) -> T;
  ```
- **Mathematical Discretization**:
  - Assert $H \in (0.01, 0.499)$ and $n_{\text{steps}} \ge 2$. Let $\Delta t = T / n_{\text{steps}}$ (preserving `T: Real` time scaling on `expiry`).
  - At each step $j = 1, \dots, n_{\text{steps}}$, generate correlated 2D Gaussian increments $(\Delta W_j^{(1)}, \widehat{W}_j)$ matching the exact $2 \times 2$ covariance:
    $$\text{Var}(\Delta W_j^{(1)}) = \Delta t, \qquad \text{Var}(\widehat{W}_j) = \Delta t^{2H}, \qquad \text{Cov}(\Delta W_j^{(1)}, \widehat{W}_j) = \frac{\sqrt{2H}}{H + 0.5}\Delta t^{H + 0.5}.$$
  - Compute the Riemann–Liouville fBM $\widetilde{W}_{t_k}^H$ ($k = 1, \dots, n_{\text{steps}}$) using the exact singular head $\widehat{W}_k$ and optimal power-law weights $b_m = \left(\frac{m^{H+0.5} - (m-1)^{H+0.5}}{H + 0.5}\right)^{\frac{1}{H - 0.5}}$ for $m = k + 1 - j \ge 2$:
    $$\widetilde{W}_{t_k}^H = \widehat{W}_k + \sqrt{2H}\sum_{j=1}^{k-1} (b_{k+1-j}\Delta t)^{H - 0.5} \Delta W_j^{(1)}, \qquad V_{t_k} = \xi_0 \exp\!\left(\eta \widetilde{W}_{t_k}^H - \frac{1}{2}\eta^2 t_k^{2H}\right).$$
  - Apply **Romano–Touzi Conditional Black–Scholes** on each path (and its antithetic twin $(-\Delta W^{(1)}, -\widehat{W})$):
    $$S_T^{\parallel} = S_0 \exp\!\left((r - q)T + \rho\sum_{k=0}^{n_{\text{steps}}-1}\sqrt{V_{t_k}}\Delta W_{k+1}^{(1)} - \frac{1}{2}\rho^2\sum_{k=0}^{n_{\text{steps}}-1} V_{t_k}\Delta t\right), \qquad Q_T = (1 - \rho^2)\sum_{k=0}^{n_{\text{steps}}-1} V_{t_k}\Delta t,$$
    and evaluate the Black–Scholes call/put formula at spot $S_T^{\parallel}$, strike $K$, rate $0$ (since $(r-q)T$ is already in $S_T^{\parallel}$, discounting the payoff by $e^{-rT}$), and total variance $Q_T$ in `T: Real`!

---

## Deliverable 2: Multi-Asset Worst-Of Autocallable / Phoenix / Snowball Engine with Fuzzy-AAD & Capriotti–Giles Adjoint Cholesky Cega (`crates/ql-exotics/src/autocallable.rs`)

Create `crates/ql-exotics/src/autocallable.rs` and re-export its public API in `crates/ql-exotics/src/lib.rs`:

- **Types & Signatures**:
  ```rust
  #[derive(Clone, Debug)]
  pub struct WorstOfAutocallableSpec<T: Real> {
      pub notional: T,
      pub spots: Vec<T>,
      pub vols: Vec<T>,
      pub divs: Vec<T>,
      pub rate: T,
      pub correlation: Vec<Vec<T>>,
      pub observation_times: Vec<f64>,
      pub autocall_barriers: Vec<T>,
      pub redemption_bonuses: Vec<T>,
      pub coupon_barriers: Vec<T>,
      pub coupon_rates: Vec<T>,
      pub memory_coupon: bool,
      pub knock_in_barrier: T,
      pub put_strike: T,
      pub observe_ki_all_dates: bool,
      pub smoothing_half_width: f64,
      pub smooth_min_temp: f64,
  }

  #[derive(Clone, Debug)]
  pub struct AutocallableRiskReport {
      pub price: f64,
      pub stderr: f64,
      pub autocall_probabilities: Vec<f64>,
      pub expected_life: f64,
      pub spot_deltas: Vec<f64>,
      pub vegas: Vec<f64>,
      pub rate_rho: f64,
      pub corr_cega: Vec<Vec<f64>>,
      pub corr_cega_stderr: Vec<Vec<f64>>,
  }

  pub fn price_worst_of_autocallable<T: Real>(
      spec: &WorstOfAutocallableSpec<T>,
      n_paths: usize,
      seed: u64,
  ) -> Result<T, String>;

  pub fn price_worst_of_autocallable_aad(
      spec: &WorstOfAutocallableSpec<f64>,
      n_paths: usize,
      n_bins: usize,
      seed: u64,
  ) -> Result<AutocallableRiskReport, String>;
  ```
- **Mathematical Specification**:
  1. Factor `spec.correlation` via `cholesky_lower_generic(&spec.correlation)?`.
  2. For each path $p = 0..n_{\text{paths}}$ (using antithetic pairs $+U$ and $-U$), simulate correlated log-Euler GBM across $0 = t_0 < t_1 < \dots < t_M = T$.
  3. At each observation date $t_m$ ($m = 1, \dots, M$):
     - Compute relative performances $R_{i,m} = S_{i,t_m} / S_{i,0}$ and smooth worst-of performance $W_m = \text{smooth\_min\_slice}(\&R_m, \text{spec.smooth\_min\_temp})$.
     - Evaluate $C^2$ Fuzzy-AAD indicators using `smooth_heaviside_c2`:
       $$p^{\text{call}}_m = W_m.\text{smooth\_heaviside\_c2}(B^{\text{call}}_m, \epsilon), \qquad p^{\text{cpn}}_m = W_m.\text{smooth\_heaviside\_c2}(B^{\text{cpn}}_m, \epsilon).$$
     - Update Phoenix memory coupon bank $M_m$:
       - If `memory_coupon`: $C_m = p^{\text{cpn}}_m (c_m + M_{m-1})$ and $M_m = (1 - p^{\text{cpn}}_m)(c_m + M_{m-1})$.
       - Else: $C_m = p^{\text{cpn}}_m c_m$.
     - Accumulate discounted cashflow:
       $$\text{PV} \mathrel{+}= \text{Alive}_{m-1} e^{-r t_m} \text{Notional} \cdot \left[ C_m + p^{\text{call}}_m (1 + r^{\text{red}}_m) \right], \qquad \text{Alive}_m = \text{Alive}_{m-1}(1 - p^{\text{call}}_m).$$
  4. At maturity $t_M = T$:
     - Knock-in indicator: if `observe_ki_all_dates`, $I_{\text{KI}} = 1 - \prod_{m=1}^M (1 - B^{\text{KI}}.\text{smooth\_heaviside\_c2}(W_m, \epsilon))$; else $I_{\text{KI}} = B^{\text{KI}}.\text{smooth\_heaviside\_c2}(W_M, \epsilon)$.
     - Terminal redemption:
       $$\text{PV} \mathrel{+}= \text{Alive}_M e^{-r T} \text{Notional} \cdot \left[ 1 - I_{\text{KI}} \cdot (K_{\text{put}} - W_M).\text{max\_val}(0) \right].$$
  5. In `price_worst_of_autocallable_aad`, accumulate the pathwise adjoint $\bar{L}$ with respect to the Cholesky factor $L$ (either via `ArenaTape` / `Dual` seeded on `spots`, `vols`, `rate`, and the lower-triangular entries of $L$, or via explicit pathwise pullback $\bar{L}_{ij} = \frac{1}{N}\sum_{p,m} \bar{Z}_{i,m}^{(p)} U_{j,m}^{(p)}$) and project $\bar{L} \mapsto \text{corr\_cega}$ via `cholesky_corr_cega_adjoint(&spec.correlation, &l_bar)?`, computing jackknife/bin standard errors across `n_bins.max(2)` bins.

---

## Deliverable 3: 2-Factor `G2++` Bermudan Swaption Engine & Generalized Forward Market Model (FMM / LMM) for TARNs and Callable Range Accruals (`crates/ql-rates/src/g2.rs`, `crates/ql-rates/src/lmm.rs`)

### 3.1 2-Factor `G2++` Bermudan Swaption Engine (`crates/ql-rates/src/g2.rs`, re-exported in `crates/ql-rates/src/lib.rs`)
Implement `G2ppBermudanSwaptionEngine` following Yamakami & Takeuchi (2022) (`papers/50`) and Brigo & Mercurio Chapter 4.2:
```rust
#[derive(Clone, Copy, Debug)]
pub struct G2ppBermudanSwaptionEngine {
    pub state_grid_points: usize, // e.g. 21 per factor (21 x 21 tensor grid)
    pub quad_order: usize,        // e.g. 7 or 9 Gauss-Hermite nodes per factor
}

impl Default for G2ppBermudanSwaptionEngine {
    fn default() -> Self {
        Self { state_grid_points: 21, quad_order: 7 }
    }
}

impl G2ppBermudanSwaptionEngine {
    pub fn new(state_grid_points: usize, quad_order: usize) -> Self;

    pub fn price<T: Real, C: YieldCurve<T>>(
        &self,
        model: &G2PlusPlus<T>,
        curve: &C,
        bermudan: &BermudanSwaption<T>,
    ) -> T;

    pub fn price_dual_curve<T: Real, C1: YieldCurve<T>, C2: YieldCurve<T>>(
        &self,
        model: &G2PlusPlus<T>,
        ois_discount: &C1,
        ibor_forward: &C2,
        bermudan: &BermudanSwaption<T>,
    ) -> T;
}
```
- **Mathematical Rollback under Terminal Bond Measure $\mathbb{Q}^{T_N}$**:
  1. Sort and deduplicate `exercise_times` inside $[T_0, T_N)$. If empty, return `T::zero()`.
  2. Work with the numeraire-deflated value $\tilde{V}(t_e, x, y) = \frac{V(t_e, x, y)}{P_{\text{ois}}(t_e, T_N, x, y)}$, where $P_{\text{ois}}(t, T, x, y)$ is evaluated via `model.discount_bond(ois_discount, t, T, x, y)`.
  3. At each exercise date $t_e$, the intrinsic exercise payoff $E(t_e, x, y)$ (for remaining swap dates $T_{k_e}, \dots, T_N$, supporting both `SwaptionSettlement::Physical` and `SwaptionSettlement::CashParYield` and dual curves `ois_discount` / `ibor_forward` exactly as in `bermudan_payoff` in `swaption.rs`) is computed in `T: Real`.
  4. For a single exercise date ($M = 1$), integrate $\max(E(t_{e_1}, x, y), 0) \frac{P_{\text{ois}}(0, T_N)}{P_{\text{ois}}(t_{e_1}, T_N, x, y)}$ directly at $(x_0, y_0) = (0, 0)$ using 2D Gauss–Hermite quadrature under $\mathbb{Q}^{T_N}$.
  5. For $M \ge 2$ exercise dates, build a 2D tensor grid $(x_i, y_j)$ over $[-5\sigma_x(t_{e_M}), +5\sigma_x(t_{e_M})] \times [-5\sigma_y(t_{e_M}), +5\sigma_y(t_{e_M})]$, initialize $\tilde{V}(t_{e_M}, x_i, y_j) = E(t_{e_M}, x_i, y_j) / P_{\text{ois}}(t_{e_M}, T_N, x_i, y_j)$, and roll backward from $t_{m+1}$ to $t_m$ ($\Delta t = t_{m+1} - t_m$, $\tau = T_N - t_{m+1}$) using the exact conditional Gaussian mean $(\mu_x^{T_N}, \mu_y^{T_N})$ and covariance $(v_x(\Delta t), v_y(\Delta t), c_{xy}(\Delta t))$ under $\mathbb{Q}^{T_N}$ via separable 1D Gauss–Hermite quadrature (first integrating along $y'$ for each grid row $i$, then along $x'$, reducing complexity from $O(N_x^2 Q^2)$ to $O(N_x^2 Q)$!).

### 3.2 Generalized Forward Market Model (`crates/ql-rates/src/lmm.rs`, re-exported in `crates/ql-rates/src/lib.rs`)
Create `crates/ql-rates/src/lmm.rs` implementing the Shifted-Lognormal BGM / Lyashenko–Mercurio (2019) / López-Salas et al. (2024) (`papers/58`) Forward Market Model with Hunter–Jäckel–Joshi (2001) predictor–corrector terminal-measure stepping:
- **Types & Signatures**:
  ```rust
  #[derive(Clone, Debug)]
  pub struct ForwardMarketModel<T: Real> {
      pub tenor_times: Vec<f64>,       // T_0 < T_1 < ... < T_M
      pub initial_forwards: Vec<T>,    // F_0(0), ..., F_{M-1}(0)
      pub initial_df_t0: T,            // P(0, T_0)
      pub displacements: Vec<T>,       // shift delta_i >= 0
      pub volatilities: Vec<T>,        // caplet/forward vol level sigma_i
      pub corr_decay: f64,             // beta in rho_{ij} = rho_inf + (1 - rho_inf) exp(-beta |T_i - T_j|)
      pub corr_floor: f64,             // rho_inf in [0, 1)
      pub backward_looking_rfr: bool,  // true = Lyashenko-Mercurio in-accrual decay gamma_i(t) = (T_{i+1}-t)/(T_{i+1}-T_i)
  }

  #[derive(Clone, Debug)]
  pub struct TarnSpec<T: Real> {
      pub notional: T,
      pub target_coupon: T,            // total target cap T_cap (fraction of notional, e.g. 0.08)
      pub fixed_rate: T,               // e.g. 0.09 in structured coupon max(fixed_rate - gearing * F_i, floor)
      pub gearing: T,                  // e.g. 2.0 (inverse floater) or -1.0 (bullish floater)
      pub coupon_floor: T,             // e.g. 0.0
      pub coupon_cap: T,               // e.g. 0.12
      pub smoothing_eps: f64,          // fuzzy-AAD target knockout width (e.g. 1e-3)
  }

  #[derive(Clone, Debug)]
  pub struct CallableRangeAccrualSpec<T: Real> {
      pub notional: T,
      pub coupon_rate: T,
      pub lower_bound: T,
      pub upper_bound: T,
      pub call_period_indices: Vec<usize>, // indices i >= 1 where issuer can call at T_i
      pub smoothing_eps: f64,              // smooth corridor indicator width (e.g. 1e-3)
  }

  impl<T: Real> ForwardMarketModel<T> {
      /// Discount factor P(0, T_k) reconstructed from initial_df_t0 and initial_forwards.
      pub fn discount_factor_0(&self, period_idx: usize) -> T;

      /// Price a caplet on F_i with strike K settling at T_{i+1} via Predictor-Corrector MC under Q^{T_M}.
      pub fn price_caplet_mc(
          &self,
          period_idx: usize,
          strike: T,
          n_paths: usize,
          seed: u64,
      ) -> T;

      /// Price a Target Redemption Note (TARN) with pathwise AAD support.
      pub fn price_tarn_mc(
          &self,
          spec: &TarnSpec<T>,
          n_paths: usize,
          seed: u64,
      ) -> T;

      /// Price a Bermudan Callable Range Accrual note using LSM continuation regression + envelope-theorem AAD.
      pub fn price_callable_range_accrual_mc(
          &self,
          spec: &CallableRangeAccrualSpec<T>,
          n_paths: usize,
          seed: u64,
      ) -> T;
  }
  ```

---

## Deliverable 4: Basel III / IV FRTB (MAR20–MAR33) SA & IMA Engines + Burgard–Kjaer Bilateral Asymmetric-Funding XVA (`crates/ql-xva/src/frtb.rs`, `crates/ql-xva/src/metrics.rs`)

### 4.1 FRTB Standardized Approach (SA) & Internal Models Approach (IMA) (`crates/ql-xva/src/frtb.rs`, re-exported in `crates/ql-xva/src/lib.rs`)
Implement `crates/ql-xva/src/frtb.rs` following BCBS MAR20–MAR33 and Tasche (2018) (`papers/59`):
- **Types & Signatures**:
  ```rust
  #[derive(Clone, Copy, Debug, PartialEq, Eq)]
  pub enum FrtbRiskClass {
      Girr,
      CsrNonSec,
      Equity,
      Commodity,
      Fx,
  }

  #[derive(Clone, Debug)]
  pub struct FrtbSensitivityBucket<T: Real> {
      pub risk_class: FrtbRiskClass,
      pub bucket_id: usize,
      pub delta_sensitivities: Vec<T>,
      pub delta_risk_weights: Vec<f64>,
      pub vega_sensitivities: Vec<T>,
      pub vega_risk_weights: Vec<f64>,
      pub curvature_up_pnl: Vec<T>,      // V^+(k) - V_0
      pub curvature_down_pnl: Vec<T>,    // V^-(k) - V_0
      pub curvature_shock_weights: Vec<f64>,
      pub intra_bucket_corr: Vec<Vec<f64>>,
  }

  #[derive(Clone, Debug)]
  pub struct FrtbDrcObligor<T: Real> {
      pub bucket_id: usize,
      pub notional: T,       // > 0 long credit exposure, < 0 short credit exposure
      pub mark_to_market_pnl: T,
      pub lgd: f64,          // e.g. 0.75 senior, 1.00 equity/subordinated, 0.25 covered bond
      pub maturity: f64,     // in years (clamped to [0.25, 1.0])
      pub default_risk_weight: f64, // MAR22 Table 2 weight (e.g. 0.005 AAA .. 0.15 BB)
  }

  #[derive(Clone, Debug)]
  pub struct FrtbSaSpec<T: Real> {
      pub buckets: Vec<FrtbSensitivityBucket<T>>,
      pub inter_bucket_corr: Vec<(FrtbRiskClass, usize, usize, f64)>,
      pub drc_obligors: Vec<FrtbDrcObligor<T>>,
      pub exotic_residual_notionals: Vec<T>, // 1.0% RRAO weight
      pub other_residual_notionals: Vec<T>,  // 0.1% RRAO weight
  }

  #[derive(Clone, Debug)]
  pub struct FrtbSaReport<T: Real> {
      pub sbm_low: T,
      pub sbm_med: T,
      pub sbm_high: T,
      pub sbm_charge: T,
      pub drc_charge: T,
      pub rrao_charge: T,
      pub total_capital: T,
  }

  pub fn compute_frtb_sa<T: Real>(spec: &FrtbSaSpec<T>) -> FrtbSaReport<T>;

  #[derive(Clone, Debug)]
  pub struct FrtbImaSpec<T: Real> {
      pub position_weights: Vec<T>,
      /// Base 10-day scenario P&L shocks per liquidity horizon cascade j in 0..5
      /// corresponding to LH = [10, 20, 40, 60, 120] days: shape [5][n_scenarios][n_positions].
      pub horizon_scenario_pnls: [Vec<Vec<T>>; 5],
      pub confidence: f64, // default 0.975 (MAR33.3)
      pub ses_addon: T,    // non-modellable risk factor (NMRF) capital add-on
  }

  #[derive(Clone, Debug)]
  pub struct FrtbImaReport<T: Real> {
      pub horizon_es: [T; 5],
      pub liquidity_adjusted_es: T,
      pub total_ima_capital: T,
      pub marginal_es: Vec<f64>,
      pub euler_component_es: Vec<f64>,
  }

  pub fn compute_frtb_ima<T: Real>(spec: &FrtbImaSpec<T>) -> FrtbImaReport<T>;
  ```
- **Exact Regulatory Formulas**:
  - **Low / Medium / High Correlation Scenarios (MAR21.6–MAR21.7)**:
    $$\rho^{(\text{Low})} = \max(2\rho - 1.0,\; 0.75\rho), \qquad \rho^{(\text{Med})} = \rho, \qquad \rho^{(\text{High})} = \min(1.25\rho,\; 1.0).$$
  - **MAR21.4(5)(b) Non-Negative Inter-Bucket Radicand Fallback**: If $\sum_b K_b^2 + \sum_{b \ne b'} \gamma_{bb'}^{(c)} S_b S_{b'} < 0$, replace $S_b$ with $\hat{S}_b = \max(\min(S_b, K_b), -K_b)$.
  - **Curvature Charge (MAR21.5 / MAR21.98–MAR21.105)**: Compute $\text{CVR}_k^+ = -(V^+(k) - V_0 - \text{RW}_k^{(\text{Curv})} s_k)$ and $\text{CVR}_k^- = -(V^-(k) - V_0 + \text{RW}_k^{(\text{Curv})} s_k)$, use squared correlations $(\rho^{(c)})^2, (\gamma^{(c)})^2$, and apply the MAR21.101 indicator $\psi(x, y) = 0$ if $x < 0 \land y < 0$ else $1$.
  - **DRC (MAR22)**: Apply maturity weight $w_i = \text{clamp}(T_i, 0.25, 1.0)$, compute obligor $\text{NetJTD}_j$, bucket Hedge Benefit Ratio $\text{HBR}_b = \frac{\sum_{j: \text{NetJTD}_j > 0} \text{NetJTD}_j}{\sum_{j: \text{NetJTD}_j > 0} \text{NetJTD}_j + \sum_{j: \text{NetJTD}_j < 0} |\text{NetJTD}_j|}$, and $\text{DRC}_b = \max\left(0, \sum_{\text{NetJTD}_j > 0} \text{RW}_j \text{NetJTD}_j - \text{HBR}_b \sum_{\text{NetJTD}_j < 0} \text{RW}_j |\text{NetJTD}_j|\right)$.
  - **IMA Liquidity-Horizon Cascade (MAR33.4)**:
    $$\text{ES}(w) = \sqrt{\text{ES}_0^2 + \sum_{j=1}^4 \text{ES}_j^2 \cdot \frac{\text{LH}_j - \text{LH}_{j-1}}{10}}, \qquad \text{LH} = [10, 20, 40, 60, 120]$$
    with exact Acerbi–Tasche 97.5% Expected Shortfall in `T: Real` and exact Euler AAD component decomposition $\sum_i \text{euler\_component\_es}[i] = \text{liquidity\_adjusted\_es.val()}$.

### 4.2 Burgard–Kjaer (2011) / Bichuch–Capponi–Sturm (2018) Bilateral Asymmetric-Funding XVA (`crates/ql-xva/src/metrics.rs`)
Add to `crates/ql-xva/src/metrics.rs` (and re-export in `crates/ql-xva/src/lib.rs`):
```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum XvaCloseOutConvention {
    RiskFree,    // Close-out on clean OIS value V
    Replacement, // Nonlinear close-out on risky/funding-adjusted value V_hat = V + U
}

#[derive(Clone, Debug)]
pub struct BilateralFundingSpec<T: Real> {
    pub lgd_counterparty: T,
    pub lgd_own: T,
    pub funding_borrow_spread: T,   // s_F^+ >= 0 (treasury borrowing cost on positive uncollateralized value)
    pub funding_lend_spread: T,     // s_F^- (treasury lending benefit on negative uncollateralized value)
    pub collateral_pay_spread: T,   // s_X^+ = r_X^+ - r_OIS on received collateral X > 0
    pub collateral_rec_spread: T,   // s_X^- = r_X^- - r_OIS on posted collateral X < 0
    pub close_out: XvaCloseOutConvention,
}

#[derive(Clone, Debug)]
pub struct BilateralXvaReport<T: Real> {
    pub cva: T,
    pub dva: T,
    pub fca: T,
    pub fba: T,
    pub colva: T,
    pub total_xva: T, // cva - dva + fca - fba + colva (cost to desk)
    pub fixed_point_iterations: usize,
}

pub fn compute_bilateral_xva<T: Real>(
    cube: &NpvCube<T>,
    csa: Option<&CsaTerms<T>>,
    ois_dfs: &[T],
    cpty_survival: &[T],
    own_survival: &[T],
    spec: &BilateralFundingSpec<T>,
) -> BilateralXvaReport<T>;
```
- Under `XvaCloseOutConvention::RiskFree`, evaluate CVA, DVA, FCA ($s_F^+ \mathbb{E}[(V_k - X_k)^+]\Delta t_k$), FBA ($s_F^- \mathbb{E}[(X_k - V_k)^+]\Delta t_k$), and COLVA directly on $V_k - X_k$.
- Under `XvaCloseOutConvention::Replacement`, solve the backward recursive contractive fixed-point equation for the valuation adjustment $U_k = -\text{XVA}_k$ at each time step $k = K-1, \dots, 0$ in `T: Real` (converging in $\le 8$ Picard/Newton iterations to $< 10^{-13}$ while preserving `T: Real` AD tangents).

---

## Deliverable 5: Smooth Rockafellar–Uryasev Mean–CVaR Optimization & Almgren–Chriss / Obizhaeva–Wang Optimal Execution (`crates/ql-portfolio/src/allocation.rs`, `crates/ql-portfolio/src/execution.rs`)

### 5.1 Smooth Rockafellar–Uryasev (2000) Mean–CVaR Portfolio Optimization (`crates/ql-portfolio/src/allocation.rs`)
Add to `crates/ql-portfolio/src/allocation.rs` (and re-export in `crates/ql-portfolio/src/lib.rs`):
```rust
#[derive(Clone, Debug)]
pub struct MeanCvarSpec<T: Real> {
    pub scenario_returns: Vec<Vec<T>>, // shape [S][N]
    pub expected_returns: Vec<T>,      // shape [N]
    pub risk_aversion: T,              // gamma > 0
    pub confidence: f64,               // alpha in (0.5, 0.999)
    pub smoothing_eps: f64,            // smooth hinge bandwidth eps (e.g. 1e-3)
    pub lower_bounds: Vec<T>,          // shape [N]
    pub upper_bounds: Vec<T>,          // shape [N]
}

#[derive(Clone, Debug)]
pub struct MeanCvarResult<T: Real> {
    pub weights: Vec<T>,
    pub var_threshold: T,
    pub cvar: T,
    pub expected_return: T,
}

pub fn try_mean_cvar_weights<T: Real>(
    spec: &MeanCvarSpec<T>,
) -> Result<MeanCvarResult<T>, AllocationError>;

pub fn mean_cvar_weights<T: Real>(spec: &MeanCvarSpec<T>) -> MeanCvarResult<T>;
```
- **Mathematical & KKT-IFT Specification**:
  1. Using the smooth Rockafellar–Uryasev hinge $\phi_\varepsilon(u) = \frac{1}{2}(u + \sqrt{u^2 + \varepsilon^2})$ with $\phi_\varepsilon'(u) = \frac{1}{2}(1 + \frac{u}{\sqrt{u^2 + \varepsilon^2}})$ and $\phi_\varepsilon''(u) = \frac{\varepsilon^2}{2(u^2 + \varepsilon^2)^{3/2}}$, define the auxiliary CVaR objective:
     $$\mathcal{J}(w, \zeta) = \gamma\!\left(\zeta + \frac{1}{(1 - \alpha) S}\sum_{s=1}^S \phi_\varepsilon\!\left(-\sum_{j=1}^N R_{s,j} w_j - \zeta\right)\right) - \sum_{i=1}^N \mu_i w_i + \frac{\tau_{\text{reg}}}{2}\sum_{i=1}^N w_i^2$$
     with $\tau_{\text{reg}} = 10^{-4}$ ensuring a strictly positive-definite Hessian on $(w, \zeta)$.
  2. **Stage 1 (Primal `f64` Active-Set damped Newton solve)**: Solve for $(w^*, \zeta^*, \lambda^*)$ and the active bound partition $(\mathcal{L}, \mathcal{U}, \mathcal{F})$ in `f64` (verifying $\sum_i w_i^* = 1$ and $l_i \le w_i^* \le u_i$).
  3. **Stage 2 (Single-Pass `T: Real` Reduced-KKT IFT Lift)**:
     - Set $w_i = l_i$ for $i \in \mathcal{L}$ and $w_i = u_i$ for $i \in \mathcal{U}$ directly in `T: Real`.
     - Assemble the reduced stationarity residual $G_T(w_{\mathcal{F}}^*, \zeta^*, \lambda^*)$ in `T: Real` and the primal KKT Hessian matrix $H_{\text{KKT}} \in \mathbb{R}^{(|\mathcal{F}|+2)\times(|\mathcal{F}|+2)}$ at $(w^*, \zeta^*, \lambda^*)$, and solve the single linear system $z_T = z^* - H_{\text{KKT}}^{-1} G_T(z^*)$ via `solve_linear_generic` so `Dual<N>` and `TapeVar<f64>` propagate exact sensitivities with respect to `scenario_returns`, `expected_returns`, `risk_aversion`, `lower_bounds`, and `upper_bounds` without taping Newton iterations!

### 5.2 Almgren–Chriss (2001) & Obizhaeva–Wang (2013) Optimal Execution (`crates/ql-portfolio/src/execution.rs`)
Create `crates/ql-portfolio/src/execution.rs` and re-export in `crates/ql-portfolio/src/lib.rs`:
```rust
#[derive(Clone, Debug)]
pub struct AlmgrenChrissSpec<T: Real> {
    pub initial_inventory: T,
    pub horizon: T,
    pub num_steps: usize,
    pub volatility: T,
    pub temporary_impact: T,
    pub permanent_impact: T,
    pub risk_aversion: T,
}

#[derive(Clone, Debug)]
pub struct ObizhaevaWangSpec<T: Real> {
    pub initial_inventory: T,
    pub trade_times: Vec<f64>,
    pub volatility: T,
    pub transient_impact: T,
    pub resilience_speed: T,
    pub temporary_impact: T,
    pub risk_aversion: T,
}

#[derive(Clone, Debug)]
pub struct OptimalExecutionReport<T: Real> {
    pub times: Vec<f64>,
    pub inventory_schedule: Vec<T>,
    pub trade_schedule: Vec<T>,
    pub expected_cost: T,
    pub cost_variance: T,
    pub mean_variance_utility: T,
}

pub fn try_almgren_chriss_execution<T: Real>(
    spec: &AlmgrenChrissSpec<T>,
) -> Result<OptimalExecutionReport<T>, AllocationError>;

pub fn almgren_chriss_execution<T: Real>(
    spec: &AlmgrenChrissSpec<T>,
) -> OptimalExecutionReport<T>;

pub fn try_obizhaeva_wang_execution<T: Real>(
    spec: &ObizhaevaWangSpec<T>,
) -> Result<OptimalExecutionReport<T>, AllocationError>;

pub fn obizhaeva_wang_execution<T: Real>(
    spec: &ObizhaevaWangSpec<T>,
) -> OptimalExecutionReport<T>;
```
- **Mathematical Invariants**:
  - **Almgren–Chriss**: With $\tau = T/N$ and $\tilde{\eta} = \eta - \frac{1}{2}\gamma_{\text{perm}}\tau > 0$, compute discrete urgency $\kappa = \frac{1}{\tau}\ln(u + \sqrt{u^2 - 1})$ where $u = 1 + \frac{\lambda_{\text{risk}}\sigma^2\tau^2}{2\tilde{\eta}}$ (using the smooth TWAP limit $x_k = X_0(1 - k/N)$ when $\frac{\lambda_{\text{risk}}\sigma^2\tau^2}{\tilde{\eta}} < 10^{-10}$), inventory $x_k = X_0 \frac{\sinh(\kappa(T - t_k))}{\sinh(\kappa T)}$, trades $n_k = x_{k-1} - x_k$, expected cost $\frac{1}{2}\gamma_{\text{perm}} X_0^2 - \frac{1}{2}\gamma_{\text{perm}}\sum_{k=1}^N n_k^2 + \frac{\eta}{\tau}\sum_{k=1}^N n_k^2$, and variance $\sigma^2 \tau \sum_{k=1}^{N-1} x_k^2$.
  - **Obizhaeva–Wang**: Solve the $(N+1)\times(N+1)$ KKT linear system in `T: Real` with transient propagator kernel $G(t_i, t_j) = \theta_0 \exp(-\rho |t_i - t_j|)$, temporary impact $2\eta \delta_{ij}$, holding-risk covariance $\lambda_{\text{risk}}\sigma^2 \sum_{k=\max(i,j)}^{N-2}(t_{k+1} - t_k)$, and budget constraint $\sum_{i=0}^{N-1} n_i = X_0$. Verify U-shaped symmetry ($n_k = n_{N-1-k}$ and $n_0 = n_{N-1} > n_1 > 0$) when $\lambda_{\text{risk}} = 0$ on a uniform grid.

---

## Deliverable 6: Benchmarks 23 & 24 (`crates/ql-benches/src/suites.rs`) & Workstation Tabs 21 & 22 (`crates/ql-terminal/src/lib.rs`)

### 6.1 Benchmarks 23 & 24 (`crates/ql-benches/src/suites.rs`)
Add Suites 23 & 24 to `run_all_benchmarks` (`crates/ql-benches/src/suites.rs`):
1. **Suite 23 (`23_slv_rbergomi_and_autocallable_cega_aad`)**:
   - **Workload**: Combine (a) a 4-asset Worst-Of Phoenix Autocallable sensitivity pass (`spots[0..4]`, `vols[0..4]`, and all 6 off-diagonal correlations $\rho_{ij}$) and (b) an rBergomi hybrid Volterra sensitivity pass (`spot`, `xi0`, `eta`, `rho`).
   - **Baseline**: Central finite-difference perturbations of all 14 Autocallable inputs + 4 rBergomi inputs (using matching random seeds and smoothing parameters so the objective is smooth).
   - **SOTA**: Single-pass `price_worst_of_autocallable_aad` (using Capriotti–Giles `cholesky_corr_cega_adjoint`) + single `Dual<4>` pass of `rbergomi_hybrid_mc_price`.
   - **Parity Gate**: Checksum of primal + all 18 sensitivities agrees within `1e-4`, achieving a substantial speedup.
2. **Suite 24 (`24_frtb_sbm_ima_and_cvar_execution_aad`)**:
   - **Workload**: Combine (a) `compute_frtb_sa` sensitivity with respect to 8 delta/curvature inputs, (b) 6-asset `mean_cvar_weights` sensitivity with respect to expected returns $\mu_0..\mu_5$, and (c) 8-step `obizhaeva_wang_execution` sensitivity with respect to `(initial_inventory, transient_impact, resilience_speed, temporary_impact)`.
   - **Baseline**: Central finite-difference bumps across all $8 + 6 + 4 = 18$ parameters.
   - **SOTA**: `Dual<8>` pass of `compute_frtb_sa` + `Dual<6>` KKT-IFT pass of `mean_cvar_weights` + `Dual<4>` KKT pass of `obizhaeva_wang_execution`.
   - **Parity Gate**: Checksum of primal + all 18 derivatives agrees within `1e-5`, with $\ge 10\times$ speedup.
3. Update `crates/ql-benches/tests/benches_suite.rs:9` from `assert_eq!(rows.len(), 22);` to `assert_eq!(rows.len(), 24);`, run `cargo run -p ql-benches --release` to regenerate `crates/ql-benches/BENCHMARK_RESULTS.md` and `knowledge_base/08_BENCHMARK_RESULTS.md`, and document Suites 23 & 24 in `knowledge_base/07_BENCHMARK_METHODOLOGY_AND_STUDY.md` and `knowledge_base/05_BLUEPRINT_MODERN_RUST_QUANTLIB.md`.

### 6.2 Workstation Tabs 21 & 22 (`crates/ql-terminal/src/lib.rs`)
1. Extend `TabId` (`crates/ql-terminal/src/lib.rs:82-151`) with:
   - `SlvRoughWorstOf` (`21`, `"21 SLV / rBergomi / WoF"`)
   - `FrtbCvarExecution` (`22`, `"22 FRTB / CVaR / Exec"`)
   updating `TabId::ALL` to `[Self; 22]`.
2. In `WorkstationFrameData` and `DeskWorkstation::step_tick`:
   - Compute live scenario-responsive metrics for Tab 21 (calibrated SLV leverage / rBergomi price + Worst-Of Autocallable price & correlation Cega) and Tab 22 (FRTB-SA total capital + Smooth Rockafellar–Uryasev Mean–CVaR + Obizhaeva–Wang optimal execution cost).
   - Add blotter rows, tab chart series, and Quadrant-3 (`left_bottom`) summary lines (`"21 SLV/rBergomi/WoF: ..."` and `"22 FRTB/CVaR/Exec: ..."`), and update `crates/ql-terminal/tests/terminal_suite.rs` to assert `TabId::ALL.len() == 22` and verify that all 22 tabs render cleanly.

---

## Verification Checklist Before Committing
1. `cargo clippy --workspace --all-targets -- -D warnings` passes with **zero warnings**.
2. `cargo test --workspace` passes with **all unit, integration, and doc tests green**, including dedicated new Phase 12 tests in:
   - `crates/ql-models/tests/models_suite.rs` (SLV particle martingale conservation, $\eta_{\text{mix}} = 0$ local-vol limit, and rBergomi hybrid MC + `Dual` sensitivities vs FD)
   - `crates/ql-exotics/tests/exotics_suite.rs` (Worst-Of Autocallable Phoenix memory vs non-memory monotonicity, step-down barrier monotonicity, and Capriotti–Giles `corr_cega` vs central FD)
   - `crates/ql-rates/tests/rates_suite.rs` (`G2ppBermudanSwaptionEngine` 1-exercise parity with `G2PlusPlus::swaption_price`, 1F degenerate limit parity with `HullWhiteBermudanTreeEngine`, and `ForwardMarketModel` caplet/discount martingale test + TARN target monotonicity + Callable Range Accrual callability discount)
   - `crates/ql-xva/tests/xva_suite.rs` (`compute_frtb_sa` Low/Med/High scenario envelope + MAR21.4(5)(b) radicand guard + DRC Hedge Benefit Ratio, `compute_frtb_ima` Euler additivity $\sum_i \text{CC}_i = \text{ES}$, and `compute_bilateral_xva` symmetric parity with `compute_xva` + replacement close-out fixed-point convergence)
   - `crates/ql-portfolio/tests/allocation_suite.rs` (`mean_cvar_weights` budget/box feasibility + `Dual` KKT-IFT sensitivities vs central FD, `almgren_chriss_execution` TWAP limit & risk-aversion front-loading, and `obizhaeva_wang_execution` U-shape symmetry & `Dual` sensitivities vs central FD)
3. `cargo run -p ql-benches --release` completes all **24 suites** with zero parity failures and updates `BENCHMARK_RESULTS.md` and `knowledge_base/08_BENCHMARK_RESULTS.md`.
4. All 59 PDFs in `papers/` and all markdown documents in `knowledge_base/` are tracked and committed cleanly.
