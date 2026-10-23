# Phase 10 adversarial review — Rust engineering / architecture

**Scope.** Independent review of the workspace (23 crates plus the path repositories), with `cargo check --workspace --all-targets --release`, `cargo test --workspace --all-targets`, `cargo clippy --workspace --all-targets -- -D warnings`, and `cargo fmt --all -- --check`. Check, clippy, and fmt pass on the current nightly tree; the all-target test suite enters a timing benchmark that was still running after 60 seconds. This review intentionally treats numerical APIs as production/library APIs, not just demo code.

## Findings (prioritized)

### RUST-01 — **Critical: unstable nightly is a production ABI/toolchain requirement**

* **Evidence:** `rust-toolchain.toml:1-4` hard pins `channel = "nightly"`; `crates/ql-ad/src/lib.rs:1`, `crates/ql-math/src/lib.rs:1`, `crates/ql-mc/src/lib.rs:1`, and `crates/ql-instruments/src/lib.rs:1` use `#![feature(portable_simd)]`; `ql-ad/src/lib.rs:2` conditionally enables the still-unstable `autodiff` feature. There is no `rust-version`/MSRV declaration or stable fallback feature.
* **Impact:** downstream users, distro builds, rust-analyzer, and security scanners cannot build on stable Rust; nightly regressions can silently change SIMD/ABI behavior. A library that advertises a “modern Rust” SDK should not require an unspecified moving compiler.
* **Remediation:** isolate portable SIMD behind an opt-in `simd` feature and provide a stable scalar implementation (or use a stable SIMD abstraction); remove/contain `feature(autodiff)` experimental code; pin and document an exact tested nightly in CI and set `rust-version` for the default build.
* **Tests:** CI matrix on stable and pinned nightly; `cargo test --no-default-features` on stable; compile a minimal downstream crate against every published feature combination.

### RUST-02 — **High: `Real` is not a mathematically faithful generic scalar contract**

* **Evidence:** `crates/ql-ad/src/scalar.rs:8-89` requires `Add<f64>`, `Div<f64>`, etc., hard-wiring f64 and `Send + Sync`; comparison/branch methods use `val()` (`select_ge_zero`, `smooth_heaviside`, `min_val`, `max_val`). Many kernels branch on `.val()` (`ql-vol/src/black.rs:18-24, 58-79`, `ql-instruments/src/pde1d.rs:56-84`, `ql-fixed-income/src/callable_bond.rs:61-166`).
* **Impact:** non-smooth payoffs silently choose the primal branch and cannot represent derivatives through branch decisions; `smooth_heaviside` returns constant 0/1 outside the band, which is intentional but undocumented as a subgradient policy. `PartialEq/PartialOrd` on AD values is value-only. Third-party scalar types cannot implement the trait without f64 operator shims, and f32/decimal/complex semantics are excluded.
* **Remediation:** split `Real` into a scalar algebra (`Scalar`), differentiable operations, and a comparison/branch policy; document derivative conventions at kinks. Prefer `num_traits::Float`-style associated scalar constants or an explicit `Primal` type rather than f64-only operators. Mark branchy algorithms as primal-only or provide smooth alternatives.
* **Tests:** derivative tests at positive/negative/exact kink and transition boundaries for every branch helper; compile tests with a custom `Real` wrapper; property tests for algebra identities and finite/NaN behavior.

### RUST-03 — **High: generic AD is routinely discarded by `.val()`/f64 solver round-trips**

* **Evidence:** `crates/ql-vol/src/black.rs:55-79` and `ql-vol/src/bachelier.rs:28-52` solve implied vol using `forward.val()`, `strike.val()`, `expiry.val()`, and then lift once; `ql-portfolio/src/allocation.rs:98-133` converts covariance/views/omega to `DMatrix<f64>`, calls `try_inverse().unwrap()`, and converts results back with `T::from_f64`; `markowitz_*` and `solve_f64` do the same (`allocation.rs:42-69, 136-165`). Similar conversion is visible throughout `rg '.val()|from_f64' crates/*/src` (e.g. American critical-root and inflation calibration).
* **Impact:** sensitivities of the solved parameters with respect to all inputs are absent or only first-order ad hoc; second-order/forward/reverse AD callers receive constants. This violates the documented “one numerical kernel, many scalars” design and can produce plausible but wrong risk numbers.
* **Remediation:** make calibration/linear algebra generic over a primal scalar and use an explicit IFT sensitivity object for all quote parameters; expose separate `*_primal` APIs when scalar solving is intentional. Never claim generic AD on a function that projects to f64.
* **Tests:** seed each input independently with `Dual`/`Dual2` and compare gradients/Hessians with high-precision finite differences; add a test that fails if a nonzero seeded derivative is returned as zero.

### RUST-04 — **High: reverse tape identity and lifecycle are unsafe at the API level**

* **Evidence:** `crates/ql-ad/src/tape.rs:8` uses a global `AtomicU16` tape id (wraps after 65,535); `TapeVar` exposes public `idx`/`tape_id` (`lines 9-14`); `AdjointResult::gradient` only indexes by `idx` and never validates `tape_id` (`lines 28-34`); all implicit recording is a thread-local singleton (`lines 256-312`). `start_recording` overwrites an existing tape, there is no RAII scope/guard, and an early return/panic can leave recording active. `with_active` silently starts a tape and never stops it (`lines 268-274`).
* **Impact:** after id wrap, a stale variable can pass cross-tape assertions and read/write another tape's adjoint slot; callers can forge indices; nested recording clobbers outer tapes; `TapeVar: Send + Sync` (through `Real`) invites moving a value to another thread where the TLS tape is absent and operations panic. These are correctness and availability hazards, not merely ergonomics.
* **Remediation:** use a non-recycling generation token (`u64` plus owner identity), make fields private, bind `AdjointResult` to tape identity and reject mismatches, return `Result` instead of `assert!/expect`, and make recording an RAII `Tape<'a>`/`RecordingGuard` with explicit nesting. Either remove `Send + Sync` from `TapeVar`/`Real` or define a thread-safe tape handle and forbid cross-thread vars at the type level.
* **Tests:** create >65,535 tapes and verify stale variables are rejected; forged-index and cross-tape tests; nested guard restoration; panic/unwind cleanup; cross-thread compile/runtime tests.

### RUST-05 — **High: public numerical APIs panic on caller/data errors**

* **Evidence:** schedule construction asserts (`crates/ql-time/src/schedule.rs:96-110`); curve and interpolation constructors assert (`ql-curves/src/curves.rs:65-68,145-147`, `ql-math/src/interpolation.rs:27,115,185`); XVA cube indexing asserts (`ql-xva/src/cube.rs:56-65`); risk measures assert confidence/non-empty input (`ql-portfolio/src/risk_measures.rs:15,58,69,92-110,170-171`); tape/IFT and calibration use `expect("singular IFT Jacobian")` (`ql-ad/src/ift.rs:50-65`, `ql-curves/src/bootstrap.rs:449`); scripting VM uses `expect("invalid bytecode stack")` (`ql-scripting/src/vm.rs:66-105`).
* **Impact:** malformed market data can abort an interactive terminal, service process, or batch risk run. Singular/ill-conditioned market inputs are expected in production and should be reportable, not process-fatal.
* **Remediation:** define per-crate error enums with `thiserror`; return `Result` from constructors/solvers/indexing (or provide checked `get` plus explicitly panicking `get_unchecked`); validate dimensions, positivity, and finiteness before computation. Keep `assert!` only for internal invariants and use `debug_assert!` where appropriate.
* **Tests:** malformed dates/schedules, empty and mismatched matrices, singular Jacobians, invalid bytecode, out-of-bounds cube calls, NaN/Inf and confidence edge cases must return typed errors without panic (`catch_unwind` regression test).

### RUST-06 — **High: linear algebra silently masks failure and is numerically fragile**

* **Evidence:** `ql-portfolio/src/allocation.rs:42-69` implements hand-written Gaussian elimination, uses `partial_cmp(...).unwrap()`, clamps each pivot with `.max(1e-14)`, and has no residual/conditioning check. Black–Litterman calls `try_inverse().unwrap()` (`allocation.rs:102-104`); IFT factorization uses `expect` on `LU::solve` (`ql-ad/src/ift.rs:50-65`).
* **Impact:** negative/zero pivots are replaced by a positive value, producing silently wrong weights instead of detecting an indefinite/singular covariance. NaNs panic during ordering; large condition numbers produce unbounded risk sensitivities.
* **Remediation:** use pivoted/QR/Cholesky/SVD solves from `nalgebra`, validate symmetry/positive definiteness, return diagnostics (rank, residual, condition estimate), and use `total_cmp` only where NaN policy is explicit.
* **Tests:** singular, nearly singular, indefinite, asymmetric, NaN, and ill-conditioned matrices; compare residual norms and expected `Err` classifications.

### RUST-07 — **Medium: `Date` accepts arbitrary `i32` days but conversion arithmetic can overflow/panic**

* **Evidence:** `crates/ql-time/src/date.rs:93-111` exposes `Date::from_days` for every `i32`; `civil_from_days` adds `719468` (`lines 205-214`), and `days_from_civil` multiplies year/era values (`lines 218-228`) without checked arithmetic. `add_months` uses saturating intermediate arithmetic but then `from_ymd(...).expect("valid")` (`lines 156-166`). Only a lower `MIN` is documented; there is no upper bound.
* **Impact:** extreme or corrupted serialized dates can wrap in release builds, produce invalid Y/M/D, or panic in supposedly infallible date arithmetic. Calendar/year loops can become impractically large.
* **Remediation:** define a documented bounded epoch range with `Date::MAX`; make `from_days`/`add_days` checked or return `Result`; use checked integer arithmetic and reject out-of-range civil years.
* **Tests:** round trips at min/max, `i32::MIN/MAX`, add-month/day overflow, negative years, leap days, and fuzz arbitrary serialized values.

### RUST-08 — **Medium: Portfolio DAG cache is not a robust content-addressed graph**

* **Evidence:** `crates/ql-portfolio/src/dag.rs:229-241` keys `cache` only by user-supplied trade `id`; entries are never removed when trades disappear. `hash_f64_slice` (`lines 10-17`) normalizes `0.0` but has no canonical NaN policy; hashes are 64-bit with no collision verification. `evaluate_portfolio` asserts unique IDs (`lines 223-237`) and panics on calibration/cache assumptions.
* **Impact:** stale cache memory grows with every historical ID; accidental ID reuse across trade variants relies on hash collision resistance rather than type/key checking; NaN payloads and signed zero can make logically equal snapshots hash differently/equally unexpectedly. The mutable DAG is not `Sync`, so callers need external serialization.
* **Remediation:** key by `(TradeId, trade-kind, dependency digest)` and retain/evict only active IDs; use a canonical hash (or compare the full dependency object on hash hit), expose cache invalidation and metrics, return errors from calibration, and document single-writer/threading semantics.
* **Tests:** remove/re-add IDs, mutate trade kind/NaN/-0, forced hash collision test double, concurrent read/evaluate policy, and bounded-cache property test.

### RUST-09 — **Medium: benchmark claims are statistically and semantically brittle**

* **Evidence:** `ql-benches/src/harness.rs:23-48` uses `Instant` wall-clock timing, one warmup, 31 samples, no CPU/frequency control, and parity compares only the *last* iteration's return (`lines 63-78`). `ql-benches/tests/benches_suite.rs:3-69` runs all 20 timed suites as a normal integration test, hard-fails on hardware-dependent geometric/speedup thresholds, and therefore makes `cargo test --workspace --all-targets` take >60 s in this environment. `benchmark_risk_parity` (`ql-benches/src/suites.rs:1569-1580`) computes budget bump work but adds `(f(cov)-f(cov))*0.5`, i.e. identically zero, then discards it; the result/parity therefore does not validate budget sensitivities despite the “ninety covariance/budget bumps” description.
* **Impact:** flaky CI and misleading SOTA speedups; expensive timing is coupled to default tests; parity can pass while a sensitivity path is broken. Results in `BENCHMARK_RESULTS.md` are not portable evidence.
* **Remediation:** move timing to `criterion`/`iai-callgrind` benches excluded from `cargo test`; report confidence intervals and hardware/toolchain metadata; test all iteration outputs or deterministic checksums; remove absolute speedup gates from correctness tests; fix benchmark baseline to return the measured sensitivity.
* **Tests:** run with CPU affinity/frequency variation and assert only numerical parity in CI; a dedicated benchmark job records (does not gate) distributions; add a regression where budget sensitivities are nonzero and compared.

### RUST-10 — **Medium: terminal package has an architectural dependency on benchmarks**

* **Evidence:** `crates/ql-terminal/Cargo.toml` depends on `ql-benches.workspace = true` (lines 12-31), while `ql-benches` itself links almost every quantitative crate (`ql-benches/Cargo.toml:9-32`) and exposes benchmark runners as a library.
* **Impact:** a production workstation/SDK transitively compiles benchmark fixtures and all benchmark-only code, increasing build time, binary size, attack surface, and coupling delivery to research harnesses. It also makes feature/dependency cleanup difficult.
* **Remediation:** move shared scenario/pricing fixtures into a small `ql-examples`/`ql-fixtures` crate; keep `ql-benches` dev-only (`[[bench]]`/`dev-dependencies`), and remove it from terminal runtime dependencies.
* **Tests:** `cargo tree -e normal` for terminal must not include `ql-benches`; build-size and compile-time budgets in CI.

### RUST-11 — **Medium: allocation and hot-loop APIs create avoidable allocations and bounds checks**

* **Evidence:** `ArenaTape::gradient_from` allocates a full `Vec` every reverse pass (`ql-ad/src/tape.rs:236`); XVA `compute_xva` allocates `Vec<T>` and a second sorted `Vec<f64>` per time (`ql-xva/src/metrics.rs:89-111`); `potential_future_exposure` allocates/sorts per time (`metrics.rs:49-71`); `NpvCube::get/set` assert three bounds on every inner-loop access (`ql-xva/src/cube.rs:56-65`). `simulate_hw1f_exposure_cube` is scalar nested path/netting loops (`ql-xva/src/cube.rs:70-106`) despite SIMD/Rayon dependencies elsewhere.
* **Impact:** allocation and branch overhead dominate large XVA/AD workloads and increase latency variance; peak memory scales with paths × times × netting sets and can be attacker-controlled.
* **Remediation:** reuse scratch buffers/arena storage, expose checked indexing outside hot loops plus unsafe-free internal slices, stream quantiles/select instead of full sort where appropriate, and add bounded/parallel chunked simulation APIs.
* **Tests/benchmarks:** allocation-count tests, peak-memory tests at realistic cube sizes, and criterion benchmarks for get/set, reverse tape, PFE, and cube simulation.

### RUST-12 — **Medium: no explicit input finiteness/domain policy in pricing APIs**

* **Evidence:** Black/Bachelier and FX APIs check only expiry/vol sign (`ql-vol/src/black.rs:18-24`, `ql-vol/src/bachelier.rs:16-23`, `ql-fx/src/lib.rs:23-42`) and then evaluate `ln(forward/strike)`; strike, forward, discount, rates, and all inputs may be zero/negative/NaN. Public constructors in SVI/SABR use `assert!` for limited conditions (`ql-vol/src/surface.rs:17-27,72-78`).
* **Impact:** NaNs/Infs propagate to risk reports; zero/negative strike causes undefined logs; callers cannot distinguish invalid market data from valid zero-vol intrinsic cases.
* **Remediation:** typed validation (`finite`, positive forward/strike/expiry where required, bounded correlations), explicit policy for zero/negative rates, and `Result`-returning checked constructors. Keep fast unchecked kernels internal.
* **Tests:** boundary matrix for NaN, ±Inf, zero/negative strike/forward/discount/vol/expiry and put-call parity under zero vol.

### RUST-13 — **Low/Medium: API surface and semver/maintenance metadata are under-specified**

* **Evidence:** workspace package metadata (`Cargo.toml:30-40`) has edition/version/license but no `rust-version`, repository, documentation, authors, or publish strategy; dependencies are broad (`nalgebra = "0.34"`, `serde = "1"`, `rayon = "1.11"`). Public structs expose mutable fields widely (`TapeVar`, `ArenaTape`, `NpvCube`, market specs), and most public items have no rustdoc.
* **Impact:** semver compatibility is hard to maintain, users can violate invariants, dependency updates can alter numerical results, and docs.rs/API discoverability are weak. There is no CI workflow, cargo-deny/audit policy, or reproducible benchmark metadata in the repository.
* **Remediation:** add workspace metadata/MSRV, feature-gate heavy dependencies, document invariants and units, make fields private with validated constructors, add `cargo-deny`/`cargo-audit`, and publish API/benchmark compatibility policy.
* **Tests/tooling:** `cargo doc --workspace --no-deps`, `cargo semver-checks`, `cargo deny check`, `cargo audit`, and lockfile/update CI.

### RUST-14 — **Low/Medium: SIMD and numerical portability are not characterized**

* **Evidence:** `ql-ad/src/scalar.rs` and `ql-mc/src/rng.rs` implement hand-written `portable_simd` approximations (range reduction, erf, Box–Muller); `rust-toolchain.toml` gives no target-feature policy. SIMD approximations are tested against scalar in a narrow suite only (`ql-ad/tests/ad_suite.rs`), while benchmark output claims 1e-9 parity for MC.
* **Impact:** different CPU targets/compiler versions can produce materially different tails/Greeks; no documented accuracy envelope exists for extreme normal CDF/exp/erf. Reproducibility across x86/ARM/WASM is not guaranteed.
* **Remediation:** publish approximation error bounds and target policy, use shared tested math implementations, provide scalar fallback and deterministic mode, and test denormals/overflow/tails on x86_64, aarch64, and wasm where supported.
* **Tests:** exhaustive/random ulp tests against `libm`/MPFR across ±40 normal inputs and exp limits, plus cross-target golden vectors.

## Positive checks

* `cargo check --workspace --all-targets --release` completed successfully.
* `cargo clippy --workspace --all-targets -- -D warnings` and `cargo fmt --all -- --check` pass on the pinned nightly.
* Workspace lockfile is committed and dependencies are path-pinned for the in-tree `implied-vol`/`volsurf` repositories.
* Existing AD tests cover basic Dual/Dual2/tape finite differences and SIMD exp/CDF smoke behavior; these should be extended to the failure cases above.

## Recommended remediation order

1. Make tape identity/lifecycle and error handling safe (RUST-04/RUST-05/RUST-06).
2. Define an honest scalar/AD contract and eliminate `.val()` round-trips (RUST-02/RUST-03).
3. Establish stable/pinned toolchain and feature boundaries (RUST-01/RUST-10/RUST-13).
4. Separate correctness tests from timing benchmarks and repair benchmark parity (RUST-09).
5. Harden date/input domains and optimize/reuse hot-loop allocations (RUST-07/RUST-11/RUST-12/RUST-14).
