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
| `ql-models` | `crates/ql-models` | Models | Heston, Bates, Double Heston COS, Bergomi 2F, Hull–White, rough Heston COS |
| `ql-pde` | `crates/ql-pde` | Models | 2D Heston ADI PDE |
| `ql-mc` | `crates/ql-mc` | Models | SIMD GBM, QE Heston, AMC, RQMC |
| `ql-scripting` | `crates/ql-scripting` | Risk | Payoff AST, bytecode VM, Fuzzy AAD |
| `ql-xva` | `crates/ql-xva` | Risk | Exposure cube, CVA/DVA/FVA, SIMM/DIM/MVA |
| `ql-portfolio` | `crates/ql-portfolio` | Risk | Content-addressed incremental DAG, HRP/NCO/Gerber/RMT, multi-asset Almgren–Chriss |
| `ql-credit` | `crates/ql-credit` | Asset classes | CDS, hazard, CDS/index swaptions, LHP, ASB CDO |
| `ql-inflation` | `crates/ql-inflation` | Asset classes | Seasonal CPI, Jarrow–Yildirim, inflation swaptions, ZCIIS, YYIIS, YoY caps/floors |
| `ql-fx` | `crates/ql-fx` | Asset classes | Garman–Kohlhagen, FX TARF, Vanna–Volga, FX forwards, XCCS |
| `ql-rates` | `crates/ql-rates` | Rates | Caplets, RFR/SABR, 3D SABR cube, mid-curve swaptions, G2++, Cheyette, CMS, range accruals |
| `ql-dml` | `crates/ql-dml` | Risk | Differential PCA, surrogates, Cholesky AAD |
| `ql-exotics` | `crates/ql-exotics` | Asset classes | Barriers, double barriers, compound/chooser, Asians, lookbacks, discrete dividends, timer and passport options |
| `ql-commodities` | `crates/ql-commodities` | Asset classes | Futures curves, Schwartz–Smith, 3-leg spreads |
| `ql-levy` | `crates/ql-levy` | Models | Merton, Kou, VG, CGMY COS, Lewis/Carr–Madan FFT |
| `ql-fixed-income` | `crates/ql-fixed-income` | Rates | NSS, callable/OAS, CoCo AT1 bonds, convertibles, CPI/MBS, specified pools and TBA rolls |
| `ql-benches` | `crates/ql-benches` | Delivery | 27-suite parity/timing harness |
| `ql-terminal` | `crates/ql-terminal` | Delivery | Unified SDK, 25-tab Braille/ANSI workstation |

## Numerical capabilities

The workspace includes multi-curve OIS/IBOR bootstrapping, Black/Bachelier and arbitrage-aware volatility surfaces, Heston/rough-Heston and Lévy COS engines, ADI PDEs, SIMD Monte Carlo and American MC, scripted payoffs, XVA/SIMM/DIM/MVA, credit/inflation/FX/commodity analytics, barriers/Asians/lookbacks, callable and convertible bonds, and a pure-Rust ANSI/Braille risk workstation.

## Empirical performance

The offline benchmark harness checks numerical parity against independent baselines. Results below are generated from the current release tree; run the benchmark command on your hardware for fresh timings.





| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 187501.90 (IQR 6997.60) | 2798.47 (IQR 7.18) | 67.00x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 33130.93 (IQR 114.10) | 1970.85 (IQR 25.97) | 16.81x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 233649.17 (IQR 6839.77) | 3911.40 (IQR 113.95) | 59.74x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 89099.05 (IQR 1659.72) | 29088.30 (IQR 101.48) | 3.06x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1856885.82 (IQR 31633.20) | 598525.10 (IQR 27845.75) | 3.10x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 417665.28 (IQR 10149.30) | 149259.67 (IQR 4846.65) | 2.80x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 770123.62 (IQR 13526.17) | 175216.60 (IQR 7034.50) | 4.40x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 3018741.42 (IQR 43784.97) | 298019.95 (IQR 62564.95) | 10.13x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 368935.65 (IQR 665.55) | 15745.60 (IQR 66.12) | 23.43x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 60582.25 (IQR 1161.00) | 20964.85 (IQR 72.90) | 2.89x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19485905.80 (IQR 185170.85) | 201202.73 (IQR 1886.95) | 96.85x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1759223.65 (IQR 9817.97) | 207128.98 (IQR 1598.25) | 8.49x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37005.10 (IQR 133.15) | 11446.50 (IQR 90.70) | 3.23x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13585619.95 (IQR 35366.33) | 1196401.35 (IQR 999.97) | 11.36x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 545102.95 (IQR 1854.93) | 128376.27 (IQR 180.73) | 4.25x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 206857.67 (IQR 406.77) | 53760.47 (IQR 119.85) | 3.85x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28793.28 (IQR 44.95) | 6007.95 (IQR 28.80) | 4.79x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1921057.90 (IQR 2953.65) | 12629.75 (IQR 206.67) | 152.11x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 530458.72 (IQR 597.12) | 57251.60 (IQR 1177.75) | 9.27x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7369731.33 (IQR 42367.58) | 160504.17 (IQR 912.42) | 45.92x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 546239.03 (IQR 986.35) | 264232.42 (IQR 3805.17) | 2.07x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 516778.08 (IQR 2347.92) | 8201.95 (IQR 105.15) | 63.01x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 1243635.57 (IQR 66424.98) | 141134.08 (IQR 292.15) | 8.81x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6596.73 (IQR 65.67) | 1399.62 (IQR 20.85) | 4.71x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 305248.35 (IQR 1204.30) | 46304.80 (IQR 278.58) | 6.59x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 326510.30 (IQR 526.00) | 69702.23 (IQR 537.98) | 4.68x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 108244.52 (IQR 279.07) | 67539.48 (IQR 202.73) | 1.60x | 3.449e-5 |

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
