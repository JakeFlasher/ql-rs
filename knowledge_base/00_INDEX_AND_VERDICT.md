# ql-rs: Index and Release Verdict

## Verdict

`ql-rs` is a cohesive 23-crate, MIT/Apache-2.0 quantitative-finance workspace. It keeps the domain fidelity that made QuantLib useful while replacing hidden mutation with immutable values, generic scalar kernels, exact first- and second-order IFT risk, data-oriented storage, and explicit parallel execution. The implementation is accompanied by 30 numerical benchmark suites and a 66-paper research catalog.

## Three empirical pillars

### Pillar I — C++ QuantLib architectural debt

The source census in `01_QUANTLIB_CPP_ANATOMY_AND_DEBT.md` counts the C++ surface and identifies the costs that a clean Rust architecture should avoid:

- `ql/patterns/observable.hpp:63–159` stores observer relationships in heap-allocated `std::set` graphs; `ql/patterns/observable.hpp:290–406` adds recursive mutexes in the thread-safe variant.
- `ql/instrument.hpp:107–110`, `ql/patterns/lazyobject.hpp:131`, and `ql/pricingengine.hpp:72–73` show `mutable` cached results and mutable engine arguments behind `const` pricing calls, making shared engines non-reentrant.
- `ql/settings.hpp:37–118` and `IndexManager::instance()` expose process-wide evaluation dates and fixings, which makes scenario valuation and parallel horizon calculations implicit and fragile.
- `ql/types.hpp:50–82` fixes `Real` to a global `double` typedef, preventing one binary from composing `f64`, SIMD, and AD scalar instantiations.
- `ql/instrument.hpp:147–166` uses runtime `dynamic_cast` and `std::any` result maps; `Leg` remains a vector of heap-allocated polymorphic cashflows, blocking inlining and SIMD.

The Rust replacement retains the instrument/engine separation as pure functions while making scalar type and memory layout compile-time choices.

### Pillar II — RustQuant completeness and correctness audit

The source review in `02_RUST_QUANT_ECOSYSTEM_AUDIT.md` records concrete gaps in the educational RustQuant prototype: `RustQuant_pricing/src/lib.rs` and `RustQuant_models/src/lib.rs` are effectively empty exports; `RustQuant_data/src/curves.rs:392–765` is commented out; `src/surfaces.rs` contains conflicting definitions; schedule implementation in `RustQuant_time/src/schedule.rs:65–243` is commented; and several interpolation/model paths are `todo!()`. The simulation audit also records missing `dt.sqrt()` scaling and shared jump timing in `simulation.rs:40–100`, plus least-squares and portfolio valuation shortcuts.

These observations motivate ql-rs's explicit `Resolved*` values, mathematical integration suites, no blanket lint suppressions, and independent parity baselines rather than smoke-only finite checks.

### Pillar III — 2026 ecosystem fragmentation

The comparison in `03_CROSS_LANGUAGE_BENCHMARK_ARCHITECTURES.md` shows that the ecosystem has excellent fragments but no unified open standard: `rateslib` demonstrates manifold-shared Dual/IFT multi-curve risk; `quantsupport` demonstrates scalar traits, block-arena AAD, and scripting; `stochastic-rs` demonstrates SIMD/GPU stochastic engines; `volsurf` and `implied-vol` provide focused volatility components; and Strata/ORE provide immutable trade and exposure patterns. ql-rs integrates these ideas into one native workspace with a common `Real` trait, SoA instruments, AAD/IFT, XVA, exotics, commodities, structured credit, Lévy pricing, and a terminal SDK.

## Permanent knowledge base

| Document | Scope |
|---|---|
| `00_INDEX_AND_VERDICT.md` | Verdict, empirical pillars, and release index. |
| `01_QUANTLIB_CPP_ANATOMY_AND_DEBT.md` | C++ architecture census and modernization rationale. |
| `02_RUST_QUANT_ECOSYSTEM_AUDIT.md` | Rust ecosystem completeness and correctness comparison. |
| `03_CROSS_LANGUAGE_BENCHMARK_ARCHITECTURES.md` | Immutable, AAD, XVA, and scripting architecture patterns. |
| `04_SOTA_QUANT_COMPUTING_PARADIGMS.md` | Mathematical foundations for AD, IFT, PDE, MC, and XVA. |
| `05_BLUEPRINT_MODERN_RUST_QUANTLIB.md` | Detailed 23-crate software and mathematical blueprint. |
| `06_LITERATURE_AND_REFERENCES.md` | Annotated 66-paper and repository bibliography. |
| `07_BENCHMARK_METHODOLOGY_AND_STUDY.md` | Harness design, suite descriptions, and empirical interpretation. |
| `08_BENCHMARK_RESULTS.md` | Reproducible 30-suite benchmark output. |

## Release verification

```text
cargo +nightly fmt --all -- --check
cargo +nightly clippy --workspace --all-targets -- -D warnings
cargo +nightly test --workspace --release --no-fail-fast
cargo run -p ql-benches --release
```


## Phase 18 Institutional Expansion
The Phase 18 release adds MtM-XCCS and tenor-basis analytics, CTD-COLVA/WWR, ZABR/PDV volatility, KIKO/window barriers, Everest and capped timer options, weather degree-day and clean spark/dark spreads, minimum-torsion ENB and entropy pooling. The workstation now exposes Tab 28 and the benchmark suite contains 30 AD parity workloads.
