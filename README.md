# ql-rs — Modern Data-Oriented Quantitative Finance in Rust

`ql-rs` v0.2.0 is a 23-crate quantitative-finance workspace designed for transparent numerical methods, fast risk, and safe composition. It combines immutable value types with generic scalar kernels, forward and reverse automatic differentiation, Henrard-style IFT calibration, SIMD-friendly structure-of-arrays engines, and an incremental portfolio DAG.

## Design principles

- **Value semantics:** market snapshots, resolved trades, curves, surfaces, and valuation results are explicit values with no observer singletons or hidden mutable global state.
- **One numerical kernel, many scalars:** `Real`/`SimdReal` support ordinary `f64`, SIMD lanes, stack Dual/Dual2 numbers, and reverse tape variables.
- **IFT at solver boundaries:** calibration iterations stay in primal space; exact Jacobians and cached linear algebra map quote risk through converged roots.
- **Data-oriented execution:** SoA cashflows, grids, paths, and arenas maximize locality and vectorization while keeping hot loops allocation-free.
- **Incremental risk:** `PortfolioDag` uses content hashes and Rayon to recompute only trades affected by a market change.

## Workspace architecture

| Crate | Path | Layer | Key algorithms and types |
|---|---|---|---|
| `ql-time` | `crates/ql-time` | Foundation | Date(i32), BitsetCalendar POPCNT, schedules |
| `ql-ad` | `crates/ql-ad` | Foundation | Real/SimdReal, Dual, Dual2, ArenaTape, IFT |
| `ql-math` | `crates/ql-math` | Foundation | Roots, interpolation, distributions, quadrature |
| `ql-curves` | `crates/ql-curves` | Curves | OIS/IBOR bootstrap and IFT sensitivities |
| `ql-vol` | `crates/ql-vol` | Curves | Jäckel, SVI/eSSVI/SABR and Greeks |
| `ql-instruments` | `crates/ql-instruments` | Contracts | SoA legs, swaps, bonds, options |
| `ql-models` | `crates/ql-models` | Models | Heston, Bates, Hull–White, rough Heston COS |
| `ql-pde` | `crates/ql-pde` | Models | 2D Heston ADI PDE |
| `ql-mc` | `crates/ql-mc` | Models | SIMD GBM, QE Heston, AMC, RQMC |
| `ql-scripting` | `crates/ql-scripting` | Risk | Payoff AST, bytecode VM, Fuzzy AAD |
| `ql-xva` | `crates/ql-xva` | Risk | Exposure cube, CVA/DVA/FVA, SIMM/DIM/MVA |
| `ql-portfolio` | `crates/ql-portfolio` | Risk | Content-addressed incremental DAG |
| `ql-credit` | `crates/ql-credit` | Asset classes | CDS, hazard, LHP, ASB CDO, CDS options |
| `ql-inflation` | `crates/ql-inflation` | Asset classes | Seasonal CPI, ZCIIS, YYIIS, YoY caps/floors |
| `ql-fx` | `crates/ql-fx` | Asset classes | Garman–Kohlhagen, Vanna–Volga, FX forwards, XCCS |
| `ql-rates` | `crates/ql-rates` | Rates | Caplets, RFR/SABR, swaptions, G2++, Cheyette, CMS |
| `ql-dml` | `crates/ql-dml` | Risk | Differential PCA, surrogates, Cholesky AAD |
| `ql-exotics` | `crates/ql-exotics` | Asset classes | Barriers, double barriers, compound/chooser, Asians, lookbacks |
| `ql-commodities` | `crates/ql-commodities` | Asset classes | Futures curves, Schwartz–Smith, spreads |
| `ql-levy` | `crates/ql-levy` | Models | Merton, Kou, VG, CGMY COS |
| `ql-fixed-income` | `crates/ql-fixed-income` | Rates | NSS, callable/OAS, convertibles, CPI/MBS |
| `ql-benches` | `crates/ql-benches` | Delivery | 25-suite parity/timing harness |
| `ql-terminal` | `crates/ql-terminal` | Delivery | Unified SDK and Braille/ANSI workstation |

## Numerical capabilities

The workspace includes multi-curve OIS/IBOR bootstrapping, Black/Bachelier and arbitrage-aware volatility surfaces, Heston/rough-Heston and Lévy COS engines, ADI PDEs, SIMD Monte Carlo and American MC, scripted payoffs, XVA/SIMM/DIM/MVA, credit/inflation/FX/commodity analytics, barriers/Asians/lookbacks, callable and convertible bonds, and a pure-Rust ANSI/Braille risk workstation.

## Empirical performance

The offline benchmark harness checks numerical parity against independent baselines. Results below are generated from the current release tree; run the benchmark command on your hardware for fresh timings.





| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 191359.15 (IQR 825.95) | 2906.12 (IQR 4.22) | 65.85x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 33165.85 (IQR 66.85) | 2049.97 (IQR 9.35) | 16.18x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 240304.10 (IQR 254.90) | 4177.50 (IQR 79.22) | 57.52x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 92060.52 (IQR 212.15) | 30520.08 (IQR 1212.28) | 3.02x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1882636.20 (IQR 14183.57) | 726070.30 (IQR 1200.80) | 2.59x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 423469.50 (IQR 3806.00) | 149474.38 (IQR 3873.83) | 2.83x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 792706.43 (IQR 3418.25) | 170680.25 (IQR 1945.80) | 4.64x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 3031320.95 (IQR 34915.90) | 270267.03 (IQR 35252.75) | 11.22x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 361432.62 (IQR 926.53) | 15719.58 (IQR 194.90) | 22.99x | 1.842e-8 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 58082.75 (IQR 81.28) | 20574.50 (IQR 198.80) | 2.82x | 5.010e-5 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19518700.90 (IQR 49524.78) | 198823.65 (IQR 3407.58) | 98.17x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1757946.00 (IQR 1382.60) | 210830.80 (IQR 580.70) | 8.34x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37015.78 (IQR 55.25) | 11408.77 (IQR 25.92) | 3.24x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13553016.55 (IQR 31469.38) | 1202741.02 (IQR 853.02) | 11.27x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 545811.53 (IQR 449.57) | 128118.57 (IQR 104.03) | 4.26x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 207112.38 (IQR 150.73) | 57496.00 (IQR 215.60) | 3.60x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28978.08 (IQR 35.55) | 6003.07 (IQR 24.12) | 4.83x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1941696.55 (IQR 2432.95) | 11290.52 (IQR 179.38) | 171.98x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 528392.30 (IQR 273.97) | 56798.00 (IQR 42.00) | 9.30x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 6876535.22 (IQR 111791.90) | 148540.27 (IQR 1207.93) | 46.29x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 764350.45 (IQR 1579.50) | 293149.00 (IQR 146.50) | 2.61x | 5.596e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 504911.38 (IQR 976.08) | 7103.50 (IQR 64.03) | 71.08x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 803950.95 (IQR 6062.48) | 134035.23 (IQR 2523.83) | 6.00x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6676.00 (IQR 118.50) | 1312.60 (IQR 3.42) | 5.09x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_and_robust_markowitz_aad | VIX derivatives, Parisian barrier and robust Markowitz KKT sensitivities | 7330.05 (IQR 38.82) | 6496.85 (IQR 92.57) | 1.13x | 1.510e-8 |

## Quick start

```bash
cargo test --workspace
cargo run -p ql-benches --release
cargo run -p ql-terminal --release -- --frames 1
```

For strict CI verification:

```bash
cargo +nightly fmt --all -- --check
cargo +nightly clippy --workspace --all-targets -- -D warnings
cargo +nightly test --workspace --release --no-fail-fast
```

## Research and documentation

- `knowledge_base/00_INDEX_AND_VERDICT.md` through `08_BENCHMARK_RESULTS.md` contain the permanent architecture, numerical-method, literature, and benchmark references.
- `papers/01_*.pdf` through `papers/28_*.pdf` provide the academic and industry sources for the implemented formulas.
- `crates/ql-terminal` re-exports the stack through `ql_terminal::prelude` and can render a headless dashboard with `--plain` or `--dump-frame`.

## License

MIT OR Apache-2.0.
