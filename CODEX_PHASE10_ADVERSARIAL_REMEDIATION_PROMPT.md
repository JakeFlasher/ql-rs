# Codex Remediation Prompt — ql-rs v0.2.0 Adversarial Review

You are remediating `/home/jakeshea/rust_quantlib` after an independent workspace-wide adversarial review. Read first:

- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_MATH.md`
- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_RUST_ARCH.md`
- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_ARCHITECTURE.md`
- `knowledge_base/PHASE10_ADVERSARIAL_REVIEW_SECOND_PASS.md`

Do not begin Phase 10 round-two feature work until every issue below is fixed,
tested, and independently reviewed. Do not use blanket lint suppressions,
weaken tests, or turn numerical failure into a plausible price.

## Quantitative correctness

1. Replace `ql-exotics/src/cliquet_extensible.rs`'s ATM-call/price clipping with true local-return cliquet pricing; honor floors/caps, global floor/cap and call/put orientation. Replace `WriterExtensibleOptionEngine`'s sum of unrelated vanillas with the Longstaff/Geske extension decomposition and bivariate-CDF boundary. Add zero-vol, equal-strike/maturity and Dual-Greek tests.
2. Implement Little-Heston-Trap/continuous complex-log handling in `ql-models/src/heston.rs`; enforce CF continuity and phi(0)=1. Compute analytic Heston/Bates c2/c4 cumulants for COS intervals (including vol-of-vol/rho and jump terms), with term-convergence and independent Fourier/QuantLib tests.
3. Implement Andersen QE/QE-M exactly in `ql-mc/src/paths.rs`, including sigma_v=0 exact integrated variance, martingale tests, Feller-edge tests and Black-limit convergence.
4. Remove ASB CDO ad-hoc first-moment renormalization in `ql-credit/src/cdo.rs`; improve quadrature/recursion accuracy and test conservation, bounds and independent LHP/QuantLib convergence.
5. Implement full ISDA SIMM bucket/inter-bucket aggregation, concentration factors, PSD validation and version metadata. Split XVA by legal netting set and support calibrated hazard/funding/discount curves, CSA schedules, MPOR, collateral currency/remuneration and time-integrated EPE, with independent golden vectors.
6. Make BAW, BS1993, BS2002 and Ju-Zhong genuinely distinct in `ql-instruments/src/american.rs`; remove shared delegation and silent root fallbacks. Validate published tables and Greeks.
7. Replace `ql-models/src/hybrid.rs`'s averaged Heston/BSM placeholder with the coupled Grzelak–Oosterlee H1-HW affine CF/COS engine; test eta→0, xi→0, rho signs and Dual6 Greeks.
8. Implement all Stulz min/max branches, Zhang correlation options, Li–Deng–Zhou J0+J1+J2 and Levy basket moment matching. Remove `is_call` no-op branches; test parity, K→0 Margrabe and QuantLib goldens.
9. Validate CMS spread conditional quadrature, correlation monotonicity, shifts and cap/floor/collar parity.
10. Implement real Philox/Sobol Monte Carlo for CAT-bond exponential/Pareto perils; `price_mc` must use `n_paths` and `seed`, return sampling error, and never alias `price_analytical`. Test Pareto convergence and lambda→0.
11. Implement canonical Black–Litterman posterior covariance, SPD validation, typed linear-solve errors and generic AD/IFT sensitivities. Add He–Litterman golden tests.
12. Harden inflation interpolation endpoint convention, true helper quote bootstrap, G2 root errors/monotonic brackets, double-barrier convergence/error bounds, BGK domain validation, inverse-normal tails, IFT conditioning and Dual/ Dual2 pow/sqrt domains.

## AD, tape, scalar and solver safety

- Redesign TapeVar identity/lifecycle: non-recycling generation token, private fields, tape-bound adjoints, RAII recording guards, nesting restoration, panic cleanup and explicit cross-thread policy. Test >65k tapes, forged indices, cross-tape use and unwind.
- Remove unjustified `.val()` round-trips from generic AD/IFT APIs. Use explicit primal APIs only when intentional and validated IFT objects for lifts. Test every seeded gradient/Hessian lane.
- Document/split Real branch and kink semantics; make sqrt(0), negative/zero pow, min/max and smoothing policy explicit.
- Replace hand Gaussian elimination, pivot clipping and unwrap inverses with validated pivoted/QR/Cholesky/SVD solves and typed rank/conditioning errors.
- Replace public numerical assert!/expect/unwrap_or fallback paths with typed errors or checked constructors. Invalid curves, schedules, dates, cubes, bytecode, covariance matrices and solver roots must not abort services.
- Harden Date range/arithmetic, discount/forward/strike/vol finiteness, correlations, PSD matrices and extrapolation policy.
- Publish SIMD erf/CDF/exp/RNG error bounds, scalar fallback and deterministic cross-target vectors; test tails ±6…±40.

## Portfolio/XVA/workflow architecture

Introduce a typed/versioned MarketData/PricingContext spanning calendars, currencies, fixings, conventions, curves, collateral and model IDs. Replace primitive-only composition with validated context and units.

Redesign PortfolioDag dependencies around typed factor IDs and full dependency digests; evict removed IDs; canonicalize NaN/-0; avoid hash-only collisions; expose provenance/trade-level Greeks/PnL explain; document concurrency safety; use structured calibration errors.

Separate terminal/UI from benchmark runtime dependencies. Add non-TTY-safe plain/JSON output, resize/input handling, deterministic frame scheduling and bounded allocation. Keep benchmark fixtures in a small examples/fixtures crate.

## Benchmarks and release integrity

Move timing out of normal cargo test; use a reproducible statistical runner with CPU/toolchain/hardware metadata, confidence intervals, allocation metrics and pinned independent baselines. Check complete deterministic vectors, not only last scalar. Suite 19/20 must checksum full vectors; budget bumps must be real nonzero measurements (not `(f(cov)-f(cov))*0.5`). Generated reports must not dirty source trees during tests.

Add stable-Rust scalar build or isolate nightly SIMD, pin/document MSRV/nightly, add CI matrices, docs, fuzz/property/differential tests, cargo-deny/audit, API docs and semver policy. Remove ql-benches from terminal runtime dependencies. Reuse tape/XVA scratch buffers and add memory/latency/concurrency tests.

## Verification gates

```bash
cargo +nightly fmt --all -- --check
cargo +nightly clippy --workspace --all-targets -- -D warnings
cargo +nightly test --workspace --release --no-fail-fast
cargo run -p ql-benches --release
cargo doc --workspace --no-deps
cargo tree -e normal -p ql-terminal
git diff --check
```

Keep the four review reports and this prompt as the audit trail until remediation is signed off.
