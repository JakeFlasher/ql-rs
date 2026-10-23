# Phase 10 Comprehensive Adversarial Review — Big-Picture Design, Software Engineering & Workflow Framework

**Scope**: End-to-end architectural and software engineering audit of all 23 crates in `ql-rs`, benchmarked against **OpenGamma Strata** (`repos/Strata`), **QuantLib** (`repos/QuantLib`), **OpenSourceRisk/Engine (`ORE`)** (`repos/ORE`), **rateslib** (`repos/rateslib`), and **quantsupport** (`repos/quantsupport`).

---

## 1. Domain Model & Trade Lifecycle Architecture (`ContractSpec -> ResolvedInstrument -> PricingEngine`)

### ARCH-01 — **CRITICAL: `MarketContext<T>` is Dead Code & No Canonical `PricingEngine` / `ResolvedInstrument` Pipeline Exists**
- **Evidence**:
  - `MarketContext<T>` in `crates/ql-curves/src/market_context.rs:21-57` is never referenced outside its own file anywhere in the 23-crate workspace.
  - `VanillaSwapSpec::resolve -> ResolvedVanillaSwap` (`crates/ql-instruments/src/swap.rs:62-98`) implements calendar-adjusted schedules and `DayCounter` year fractions, yet `PortfolioDag::value_trade` (`crates/ql-portfolio/src/dag.rs:293-328`) completely bypasses it and prices swaps via a hardcoded integer loop `for m in 1..=n { let t = m as f64; ... }` with unit accrual `1.0`.
  - `bootstrap.rs` (`crates/ql-curves/src/bootstrap.rs:36-631`) and `YieldTermStructure::forward_rate` (`crates/ql-curves/src/curves.rs:10-15`) hardcode `(d - ref_date) as f64 / 365.0` instead of `DayCounter::year_fraction`.
  - `NelsonSiegelSvensson<T>` (`crates/ql-fixed-income/src/nss.rs:3-36`) defines custom `zero_rate` and `discount_factor` methods without implementing `YieldTermStructure<T>`.
- **SOTA Comparison**:
  - OpenGamma **Strata** (`ResolvedTrade`, `RatesProvider`, `DiscountingSwapProductPricer`) and **QuantLib** (`Instrument::setupArguments` $\to$ `PricingEngine::calculate`) strictly separate static contractual terms (`*Spec`), calendar/convention resolution (`Resolved*`), immutable market snapshots (`MarketEnvironment<T>`), and pricing methodology (`PricingEngine<I, T>`).
- **Required Remediation**:
  1. Promote `MarketContext<T>` (`ql-curves/src/market_context.rs`) into the canonical multi-asset `MarketEnvironment<T>` carrying `valuation_date: Date`, discount/forward curves (`DiscountCurve<T>`), dividend/repo curves, spot quotes, volatility surfaces, hazard curves, and inflation curves, with `to_ad<U: Real>(&self)` for zero-friction lifting across `f64`, `Dual<N>`, `Dual2<N>`, `TapeVar<f64>`, and `SimdReal`.
  2. Define core traits `ContractSpec` (`type Resolved; fn resolve(&self, cal: &BitsetCalendar) -> Result<Self::Resolved, QlError>`) and `PricingEngine<I, T: Real>` (`fn evaluate(&self, instrument: &I, market: &MarketEnvironment<T>) -> Result<ValuationResult<T>, QlError>`).
  3. Implement `YieldTermStructure<T>` for `NelsonSiegelSvensson<T>` and wire `DayCounter` into `YieldTermStructure::forward_rate` and `bootstrap.rs`.

### ARCH-02 — **HIGH: Duplicate & Colliding Domain Types Across Crates**
- **Evidence**:
  - **4 incompatible Interest Rate Swap representations**:
    1. `ql_instruments::VanillaSwapSpec` / `ResolvedVanillaSwap` (`crates/ql-instruments/src/swap.rs:5-98`),
    2. `ql_curves::IrsSwapInstrument` (`crates/ql-curves/src/bootstrap.rs:77-155`),
    3. `ql_portfolio::PortfolioTrade::IrsSwap` (`crates/ql-portfolio/src/dag.rs:15-21`),
    4. `ql_xva::NettingSetSwap` (`crates/ql-xva/src/cube.rs:6-12`).
  - **2 colliding `BarrierType` enums**: `ql_instruments::BarrierType` (`crates/ql-instruments/src/barrier.rs:4-9`) vs. `ql_exotics::BarrierType` (`crates/ql-exotics/src/barrier.rs:5-10`).
  - **2 colliding `DupireLocalVol` structs**: `ql_vol::DupireLocalVol` (`crates/ql-vol/src/dupire.rs:4-37`) vs. `ql_models::DupireLocalVol` (`crates/ql-models/src/dupire.rs:4-41`).
  - **62 functions suppress `#[allow(clippy::too_many_arguments)]`** across the workspace, accepting 8–12 raw scalar `T` arguments (`spot, strike, vol, rate, div, expiry, ...`) where swapping two identical-type parameters compiles silently.
- **Required Remediation**:
  - Unify `BarrierType` in `ql-instruments` and re-export it in `ql-exotics`.
  - Unify `DupireLocalVol` in `ql-vol` and re-export it in `ql-models`.
  - Make `PortfolioTrade::IrsSwap`, `IrsSwapInstrument`, and `NettingSetSwap` wrap or convert directly to/from `ResolvedVanillaSwap` (or share a single `VanillaSwapLeg` schedule representation).

---

## 2. Error Handling, Domain Validation & Silent Solver Fallbacks

### ARCH-03 — **CRITICAL: Public Constructors & Solvers Panic (`assert!`/`unwrap`/`expect`) Instead of Returning Typed `Result`**
- **Evidence**:
  - `ql-ad`: `tape.rs:119,144,232,273,293,326,348` (`assert!`/`expect`), `ift.rs:21,22,45,50,54,58,65` (`assert!` + `expect("singular IFT Jacobian")`).
  - `ql-curves`: `curves.rs:64-68,79,144-147` (`assert!` on unsorted/non-positive discount factors), `bootstrap.rs:449` (`expect("curve calibration")`).
  - `ql-portfolio`: `dag.rs:202,223-233` (`assert!` on duplicate trade IDs and `expect("curve calibration")`), `allocation.rs:53-54,103-104` (`partial_cmp().unwrap()`, `try_inverse().unwrap()`), `risk_measures.rs:15-171` (`assert!` on empty slices and confidence bounds).
  - `ql-xva`: `cube.rs:57-65` (`assert!` bounds), `simm.rs:10-25` (`assert!` on threshold/correlation shapes).
  - `ql-credit`: `lib.rs:29-31,425` (`assert!`), `cdo.rs:130` (`expect("base correlation failed")`).
  - `ql-inflation`: `lib.rs:46,56,151` (`assert!` + `expect("cpi bootstrap")`).
  - `ql-time`: `schedule.rs:97-111` (`assert!`).
  - `ql-math`: `interpolation.rs:12-153` (`assert!`), `distributions.rs:148` (`assert!(p > 0.0 && p < 1.0)`), `tridiagonal.rs:13-63` (`assert!`).
  - `ql-vol`: `surface.rs:17-79` (`assert!` in `SviParams::new` and `SabrParams::new`).
  - `ql-dml`: `regression.rs:65,126` (`assert!` + `expect("regression solve")`), `pca.rs:10` (`assert!`).
  - `ql-fixed-income`: `nss.rs:66` (`expect("NSS calibration failed")`).
  - `ql-scripting`: `vm.rs:66-105` (`expect("invalid bytecode stack")`).
- **Required Remediation**:
  - Standardize on `thiserror` error enums across crates (and a top-level `QlError` or per-crate `try_*` / `Result`-returning constructors and solvers) so malformed market data, non-SPD matrices, or unbracketed roots never crash a server or TUI process.

### ARCH-04 — **CRITICAL: Silent `.unwrap_or(...)` Solver Fallbacks Corrupt Prices & Inject Bogus IFT Adjoints**
- **Evidence**:
  - `crates/ql-fixed-income/src/bond.rs:143`: `brent_solve(...).unwrap_or(guess)` feeds an unconverged guess into `ift_lift_1d`, where `r = residual(x)` has a huge non-zero primal residual and corrupts both the yield and all AD sensitivities.
  - `crates/ql-instruments/src/american.rs:32`: `brent_solve(...).unwrap_or_else(|_| seed.max(strike.val()))` feeds an unconverged critical boundary $S^*$ into `ift_lift_1d`.
  - `crates/ql-rates/src/g2.rs:279`: `brent_solve(...).unwrap_or(0.0)` substitutes $y^* = 0$ into the Jamshidian-style G2++ swaption integrand when root-finding fails.
  - `crates/ql-exotics/src/quanto_chooser.rs:87`: `brent_solve(...).unwrap_or(s.spot.val())` substitutes `spot` as the chooser critical strike/spot root on solver failure and lifts it via `ift_lift_1d`.
  - `crates/ql-exotics/src/compound.rs:64`: `brent_solve(...).unwrap_or_else(|_| k2.val().max(1e-10))` substitutes `k2` as the Geske critical spot $S^*$ on solver failure and lifts it via `ift_lift_1d`.
  - `crates/ql-credit/src/lib.rs:337`: `.unwrap_or_else(|| DMatrix::zeros(m, m))` silently zeroes the entire hazard-curve Jacobian when LU solve fails.
- **Required Remediation**:
  - Eliminate all silent `.unwrap_or(...)` solver fallbacks; widen/adapt root brackets analytically using monotone bounds (or Newton- safeguarded bisection) and propagate `Err(...)` if a root truly does not exist.

---

## 3. Incremental Portfolio DAG (`ql-portfolio`) & XVA / SIMM Engine (`ql-xva`)

### ARCH-05 — **HIGH: `PortfolioDag` Cache Leaks Memory, Conflates Telemetry Counters, & Fails to Canonicalize Float Hashes**
- **Evidence (`crates/ql-portfolio/src/dag.rs`)**:
  - **Unbounded cache (`lines 209-265`)**: `self.cache` is never pruned when trades mature or are removed from the portfolio slice; memory grows monotonically in long-running services. Add `self.cache.retain(|id, _| seen.contains(id))` (plus explicit `evict_trade` / `clear_cache`).
  - **Inconsistent telemetry (`lines 212, 241-265`)**: `self.trade_recomputes` is never reset across `evaluate_portfolio` calls (accumulating lifetime recomputes), whereas `self.curve_recomputes` is overwritten with `0` or `1` per call (`line 218`). Track both `last_trade_recomputes` (per-call) and `total_trade_recomputes` (lifetime) cleanly.
  - **Incomplete float hashing (`lines 10-17, 49-58, 94`)**: `hash_f64_slice` normalizes `-0.0` to `0.0` (though not `NaN` bit patterns), while `hash_program` (`line 54`) and `trade_market_hash` (`line 94`) call `.to_bits()` directly without normalizing `-0.0` or canonicalizing `NaN`. Use a single `canonical_f64_bits(x: f64) -> u64` helper across all hash functions.
  - **Toy 3-pillar curve & UI label mismatch**: `MarketSnapshot` (`lines 31-35`) hardcodes 3 OIS quotes and 3 forward quotes at `[1Y, 2Y, 5Y]`, while `ql-terminal/src/lib.rs:894-907` labels those 6 vector entries (3 OIS + 3 FWD) as `["1Y", "2Y", "5Y", "10Y", "20Y", "30Y"]`! Support arbitrary pillar vectors `(tenor, ois_quote, fwd_quote)` or resolved `VanillaSwapSpec` schedules in `PortfolioDag` and fix the TUI labels (`["OIS 1Y", "OIS 2Y", "OIS 5Y", "FWD 1Y", "FWD 2Y", "FWD 5Y"]`).

### ARCH-06 — **CRITICAL: `ql-xva` CSA Collateral Formula Returns Posted Collateral Instead of Remaining Exposure, & Violates Netting-Set Boundaries**
- **Evidence**:
  - **Inverted CSA collateralized exposure (`crates/ql-xva/src/metrics.rs:11-19`)**:
    ```rust
    pub fn collateralized_exposure<T: Real>(mtm: T, threshold_cpty: f64, mta: f64) -> T {
        let e = mtm.max_val(T::zero());
        let excess = e - threshold_cpty;
        if excess.val() > mta {
            excess
        } else {
            e
        }
    }
    ```
    When `mtm = 10.0`, `threshold_cpty = 2.0`, `mta = 0.5`, the collateral called is `excess = 8.0`, so the **remaining uncollateralized exposure** is `Threshold = 2.0` (`e - excess`)! Instead, `collateralized_exposure` returns `excess = 8.0` (the collateral posted!), meaning a lower threshold *increases* reported credit exposure!
    The exact continuous-CSA uncollateralized positive exposure is:
    $$E_{\text{coll}}(t) = \begin{cases} \max(V(t) - H_{\text{cpty}},\, 0) & \text{if } V(t) - H_{\text{cpty}} \le \text{MTA} \text{ (no margin delivery yet? No: if } V(t) - H_{\text{cpty}} > \text{MTA, collateral } C = V(t) - H_{\text{cpty}} \text{ is posted, leaving } H_{\text{cpty}}\text{!)} \end{cases}$$
    Specifically, under a standard CSA with threshold $H \ge 0$, independent amount $IA \ge 0$, and minimum transfer amount $\text{MTA} \ge 0$, when $V^+ - H > \text{MTA}$, collateral $C = V^+ - H$ is posted so remaining exposure is $\max(V^+ - C, 0) = \min(V^+, H)$; when $V^+ - H \le \text{MTA}$, no transfer occurs so remaining exposure is $V^+$. (And under smooth/continuous CSA with threshold $H$, remaining exposure is $\min(\max(V, 0), H + \text{MTA})$ or $\max(0, \min(V, H))$). Wait: let's check existing unit tests in `crates/ql-xva/tests/` so we preserve/upgrade any test expectations cleanly while providing both exact CSA formulas!
  - **Netting across distinct legal netting sets (`crates/ql-xva/src/metrics.rs:34-99`)**:
    `expected_positive_exposure`, `potential_future_exposure`, and `compute_xva` sum `cube.get(p, t, ns)` across all `ns in 0..cube.n_netting_sets` **before** applying `.max_val(T::zero())` / `collateralized_exposure`. Positive exposure in Netting Set A is illegally offset against negative MTM in Netting Set B! Each netting set `ns` must be floored/collateralized independently (`sum_ns collateralized_exposure(cube.get(p, t, ns), csa[ns])`) before summing across netting sets for a counterparty.
  - **Unweighted arithmetic mean for `epe` (`metrics.rs:107`)**:
    `epe` is computed as `ee.iter().sum() / T` instead of time-weighted trapezoidal integration $\frac{1}{T_{\max} - T_0} \sum_k \frac{EE(t_{k-1}) + EE(t_k)}{2} \Delta t_k$.
  - **ISDA SIMM v2.6 compliance (`crates/ql-xva/src/simm.rs:52-63` & `chebyshev.rs:306-319`)**:
    - `simm.rs:57` omits the pairwise concentration correlation factor $f_{kl} = \frac{\min(CR_k, CR_l)}{\max(CR_k, CR_l)}$ on off-diagonal terms ($k \ne l$).
    - `simm.rs` lacks the top-level 6-risk-class aggregation $\text{SIMM} = \sqrt{\sum_{r,s} \psi_{r,s} \text{IM}_r \text{IM}_s}$ using the official ISDA SIMM v2.6 cross-risk-class correlation matrix $\psi_{r,s}$ (`InterestRate`, `CreditQualifying`, `CreditNonQualifying`, `Equity`, `Commodity`, `FX`).
    - `simm_from_iter` (`crates/ql-xva/src/chebyshev.rs:306-319`) fails to net sensitivities across trades within the same risk factor before squaring, violating basic portfolio netting in SIMM.

---

## 4. Package Boundaries, Benchmark Harness & Terminal Architecture

### ARCH-07 — **HIGH: `ql-terminal` Depends on `ql-benches` in Production & Glob-Reexports 21 Crates**
- **Evidence**:
  - `crates/ql-terminal/Cargo.toml:26` lists `ql-benches.workspace = true` as a normal production dependency.
  - `crates/ql-terminal/src/lib.rs:3-28` suppresses `#[allow(ambiguous_glob_reexports)]` and glob-reexports 21 crates (`pub use ql_benches::*; pub use ql_exotics::*; ...`).
- **Required Remediation**:
  - Remove `pub use ql_benches::*;` from `ql_terminal::prelude` and provide explicit curated module/type exports instead of ambiguous workspace-wide globs.
  - Fix the TUI Tab 17 & 18 telemetry truncation bug in `crates/ql-terminal/src/lib.rs:1411-1438` where Phase 10 telemetry lines are pushed at indices 21–22 of `right_bottom` and then discarded because only `0..20` are rendered.

### ARCH-08 — **HIGH: `cargo test` Enforces Wall-Clock Speedup Thresholds & Benchmark Binaries Dirty Git-Tracked Files**
- **Evidence**:
  - `crates/ql-benches/tests/benches_suite.rs:3-70` runs all 20 full timing suites during `cargo test` and asserts wall-clock speedup ratios (`geomean >= 1.50`, `speedup >= 1.05`, etc.), making `cargo test` slow (>60s) and flaky under CPU load.
  - `crates/ql-benches/src/main.rs:10-20` and `crates/ql-benches/benches/quant_benchmarks.rs:9-14` unconditionally overwrite tracked repository files `BENCHMARK_RESULTS.md` and `knowledge_base/08_BENCHMARK_RESULTS.md`.
- **Required Remediation**:
  - In `crates/ql-benches/tests/benches_suite.rs`, run fast numerical parity & invariant verification (`verify_benchmark_suites_parity()`) during `cargo test` unless `QL_RUN_TIMING_GATES=1` is set.
  - Only write `BENCHMARK_RESULTS.md` and `knowledge_base/08_BENCHMARK_RESULTS.md` when `--write-report` or `QL_WRITE_BENCH_REPORT=1` is explicitly passed.
