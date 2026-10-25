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
| `ql-vol` | `crates/ql-vol` | Curves | Jäckel, SVI/eSSVI/SABR, free-boundary SABR density and Greeks |
| `ql-instruments` | `crates/ql-instruments` | Contracts | SoA legs, swaps, bonds, options |
| `ql-models` | `crates/ql-models` | Models | Heston, Bates, Double Heston COS, Bergomi 2F, Hull–White, rough Heston COS, QRH, Quintic OU SPX/VIX, ZABR/Feldman SABR, free-boundary SABR, Guyon PDV, Andreasen–Huge local-volatility grid |
| `ql-pde` | `crates/ql-pde` | Models | 2D Heston ADI PDE |
| `ql-mc` | `crates/ql-mc` | Models | SIMD GBM, QE Heston, AMC, RQMC |
| `ql-scripting` | `crates/ql-scripting` | Risk | Payoff AST, bytecode VM, Fuzzy AAD |
| `ql-xva` | `crates/ql-xva` | Risk | Exposure cube, CVA/DVA/FVA, dynamic APS SIMM/MVA and MPoR cashflow-spike CVA, BA-CVA capital, CTD-COLVA, Turlakov WWR CVA/FVA |
| `ql-portfolio` | `crates/ql-portfolio` | Risk | Content-addressed incremental DAG, HRP/NCO/Schur-HRP/OAS shrinkage, minimum-torsion ENB, entropy pooling, dynamic SIMM/MVA, multi-asset Almgren–Chriss/Avellaneda–Stoikov, EVaR, CDaR, factor risk attribution, Idzorek BL, fractional Kelly, dark/lit routing |
| `ql-credit` | `crates/ql-credit` | Asset classes | CDS, hazard, CDS/index swaptions, CMCDS market model, JLT rating migration, LHP, ASB CDO, synthetic CDO base correlation, Black–Cox structural credit |
| `ql-inflation` | `crates/ql-inflation` | Asset classes | Seasonal CPI, Jarrow–Yildirim ZCIIS cap/floors, inflation swaptions, ZCIIS, YYIIS, YoY caps/floors |
| `ql-fx` | `crates/ql-fx` | Asset classes | Garman–Kohlhagen, FX TARF, Vanna–Volga, FX forwards, MtM-resetting XCCS, KIKO/window barriers, PRDC hybrid notes, accumulator strips |
| `ql-rates` | `crates/ql-rates` | Rates | Caplets, RFR/SABR, tenor-basis swaps, 3D SABR cube, mid-curve swaptions, G2++, Cheyette Markovian swaptions, 1F Markov-functional swaptions, CMS, range accruals |
| `ql-dml` | `crates/ql-dml` | Risk | Differential PCA, surrogates, Cholesky AAD |
| `ql-exotics` | `crates/ql-exotics` | Asset classes | Barriers, double barriers, compound/chooser, Asians, lookbacks, discrete dividends, gamma/skewness swaps, Heston variance options, capped timer/passport, Everest/Himalaya and Napoleon options, Atlas and Altiplano cliquets |
| `ql-commodities` | `crates/ql-commodities` | Asset classes | Futures curves, Schwartz–Smith, Clewlow–Strickland strip swaptions, weather HDD/CDD, clean spark/dark spreads, 3-leg spreads, Average Price Options (APO), take-or-pay and BESS arbitrage DP |
| `ql-levy` | `crates/ql-levy` | Models | Merton, Kou, VG, CGMY COS, Lewis/Carr–Madan FFT |
| `ql-fixed-income` | `crates/ql-fixed-income` | Rates | NSS, callable/OAS, CoCo AT1 bonds, convertibles, AFV default-coupled PDE, CPI/MBS burnout/OAS, specified pools, TBA rolls, Treasury futures CTD, Agency CMO waterfalls, ASW and IFT z-spreads |
| `ql-benches` | `crates/ql-benches` | Delivery | 32-suite parity/timing harness |
| `ql-terminal` | `crates/ql-terminal` | Delivery | Unified SDK, 30-tab Braille/ANSI workstation |

## Numerical capabilities

The workspace includes multi-curve OIS/IBOR bootstrapping, Black/Bachelier and arbitrage-aware volatility surfaces, Heston/rough-Heston and Lévy COS engines, ADI PDEs, SIMD Monte Carlo and American MC, scripted payoffs, XVA/SIMM/DIM/MVA, credit/inflation/FX/commodity analytics, barriers/Asians/lookbacks, callable and convertible bonds, and a pure-Rust ANSI/Braille risk workstation.

## Empirical performance

The offline benchmark harness checks numerical parity against independent baselines. Results below are generated from the current release tree; run the benchmark command on your hardware for fresh timings.

**Benchmark summary:** 32 suites and 30 workstation tabs.






| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 199664.20 (IQR 5486.00) | 2982.80 (IQR 77.80) | 66.94x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 33805.60 (IQR 879.80) | 2111.40 (IQR 32.00) | 16.01x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 248352.00 (IQR 1836.20) | 4720.00 (IQR 896.20) | 52.62x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 94517.00 (IQR 715.60) | 30932.40 (IQR 237.80) | 3.06x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1942051.80 (IQR 39861.40) | 620029.40 (IQR 928.00) | 3.13x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 440064.00 (IQR 4456.00) | 152183.60 (IQR 367.80) | 2.89x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 834606.60 (IQR 5888.80) | 179247.40 (IQR 3617.80) | 4.66x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 3098761.80 (IQR 172667.00) | 296486.60 (IQR 46889.40) | 10.45x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 385609.00 (IQR 903.40) | 17768.20 (IQR 3784.40) | 21.70x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 62240.60 (IQR 339.80) | 21477.80 (IQR 21.60) | 2.90x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19840731.60 (IQR 380359.80) | 210310.20 (IQR 9105.00) | 94.34x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1808124.00 (IQR 12606.80) | 216657.00 (IQR 5483.80) | 8.35x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 38469.40 (IQR 3036.80) | 13015.60 (IQR 902.80) | 2.96x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13872917.00 (IQR 266866.20) | 1204883.60 (IQR 12942.00) | 11.51x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 546889.00 (IQR 248.20) | 129030.60 (IQR 1524.20) | 4.24x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 207456.00 (IQR 816.80) | 54849.40 (IQR 707.80) | 3.78x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28772.00 (IQR 357.60) | 6005.00 (IQR 178.20) | 4.79x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1934209.80 (IQR 3808.40) | 11436.00 (IQR 165.60) | 169.13x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 530347.40 (IQR 1456.60) | 57868.20 (IQR 400.60) | 9.16x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7059612.60 (IQR 30024.60) | 150565.00 (IQR 6688.20) | 46.89x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 554702.80 (IQR 1311.40) | 267300.40 (IQR 2345.80) | 2.08x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 516502.20 (IQR 14058.00) | 8365.00 (IQR 1113.20) | 61.75x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 1024173.20 (IQR 16237.20) | 133933.00 (IQR 5628.00) | 7.65x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6865.20 (IQR 211.40) | 1526.40 (IQR 544.80) | 4.50x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 307298.60 (IQR 6519.00) | 46001.40 (IQR 2119.60) | 6.68x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 324350.80 (IQR 1748.20) | 70127.40 (IQR 2116.80) | 4.63x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 136622.20 (IQR 2647.60) | 61931.80 (IQR 3897.60) | 2.21x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3867231.80 (IQR 57836.80) | 438093.80 (IQR 1782.80) | 8.83x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 84517.20 (IQR 4475.40) | 17526.20 (IQR 971.60) | 4.82x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 360839.40 (IQR 10866.20) | 44362.20 (IQR 1466.60) | 8.13x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 976201.60 (IQR 2800.60) | 140017.60 (IQR 4305.60) | 6.97x | 3.731e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 2017168.80 (IQR 7546.00) | 291824.80 (IQR 2789.00) | 6.91x | 1.598e-8 |

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


## Phase 18 Institutional Expansion
The Phase 18 release adds MtM-XCCS and tenor-basis analytics, CTD-COLVA/WWR, ZABR/PDV volatility, KIKO/window barriers, Everest and capped timer options, weather degree-day and clean spark/dark spreads, minimum-torsion ENB and entropy pooling. Phase 19 adds burnout-MBS/OAS, CMCDS, Quintic-OU VIX, AFV convertibles, Clewlow–Strickland swaptions, Schur-HRP and nonlinear shrinkage. Phase 20 adds Sinh-OU/CIR++ RFR, callable BK FRNs, no-Armageddon credit options, WWM-XVA, mixed rough Bergomi, Choi/Li–Deng–Zhou spreads, MRJD commodities, AFS LOB, EXVaR and PCRP. The workstation now exposes Tab 30 and the benchmark suite contains 32 AD parity workloads.
