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
| `ql-benches` | `crates/ql-benches` | Delivery | 33-suite parity/timing harness |
| `ql-terminal` | `crates/ql-terminal` | Delivery | Unified SDK, 31-tab Braille/ANSI workstation |

## Numerical capabilities

The workspace includes multi-curve OIS/IBOR bootstrapping, Black/Bachelier and arbitrage-aware volatility surfaces, Heston/rough-Heston and Lévy COS engines, ADI PDEs, SIMD Monte Carlo and American MC, scripted payoffs, XVA/SIMM/DIM/MVA, credit/inflation/FX/commodity analytics, barriers/Asians/lookbacks, callable and convertible bonds, and a pure-Rust ANSI/Braille risk workstation.

## Empirical performance

The offline benchmark harness checks numerical parity against independent baselines. Results below are generated from the current release tree; run the benchmark command on your hardware for fresh timings.

**Benchmark summary:** 33 suites and 31 workstation tabs.






| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 516141.40 (IQR 294397.80) | 4797.40 (IQR 151.00) | 107.59x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 52140.20 (IQR 6410.60) | 2989.40 (IQR 37.20) | 17.44x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 254857.20 (IQR 36866.80) | 4335.00 (IQR 377.80) | 58.79x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 95431.20 (IQR 734.20) | 30919.40 (IQR 1338.20) | 3.09x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1954890.20 (IQR 89546.40) | 635511.80 (IQR 40391.40) | 3.08x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 431306.40 (IQR 46358.20) | 149527.20 (IQR 2723.40) | 2.88x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 791744.80 (IQR 12134.80) | 174322.20 (IQR 12017.00) | 4.54x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 4052317.00 (IQR 1324040.80) | 380860.00 (IQR 120011.00) | 10.64x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 385905.80 (IQR 470.60) | 16632.80 (IQR 1990.20) | 23.20x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 62410.40 (IQR 818.00) | 21515.80 (IQR 95.00) | 2.90x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19448045.80 (IQR 325953.20) | 214642.40 (IQR 3389.60) | 90.61x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1866829.20 (IQR 62822.20) | 221106.40 (IQR 2406.20) | 8.44x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37851.40 (IQR 2946.20) | 11565.40 (IQR 2441.60) | 3.27x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13256473.20 (IQR 70684.60) | 1174278.80 (IQR 2735.40) | 11.29x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 534062.40 (IQR 1212.60) | 126678.80 (IQR 863.20) | 4.22x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 203410.80 (IQR 698.80) | 53508.80 (IQR 528.20) | 3.80x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28305.20 (IQR 330.00) | 5877.60 (IQR 202.60) | 4.82x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1879011.80 (IQR 5035.80) | 11018.60 (IQR 615.60) | 170.53x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 510498.80 (IQR 1573.80) | 53852.40 (IQR 2071.20) | 9.48x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 6767939.80 (IQR 8510.40) | 144871.20 (IQR 4263.20) | 46.72x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 534528.00 (IQR 6510.80) | 258350.60 (IQR 1118.80) | 2.07x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 499831.20 (IQR 5098.00) | 7895.40 (IQR 825.40) | 63.31x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 803582.80 (IQR 5424.60) | 134516.80 (IQR 5248.40) | 5.97x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6708.20 (IQR 6278.20) | 1417.60 (IQR 184.00) | 4.73x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 294678.40 (IQR 6408.40) | 48521.80 (IQR 39402.60) | 6.07x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 315927.80 (IQR 3547.60) | 69956.20 (IQR 1431.20) | 4.52x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 130847.20 (IQR 1314.20) | 61079.80 (IQR 1267.60) | 2.14x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3740509.60 (IQR 42855.60) | 423383.60 (IQR 801.20) | 8.83x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 84543.00 (IQR 2541.00) | 16875.80 (IQR 2748.20) | 5.01x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 283286.20 (IQR 4050.60) | 35813.80 (IQR 1369.60) | 7.91x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 951571.20 (IQR 2719.80) | 139885.20 (IQR 777.60) | 6.80x | 3.731e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1627688.60 (IQR 45339.60) | 220754.00 (IQR 2371.40) | 7.37x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Phase21 institutional 12-engine AAD vs FD | 478730.80 (IQR 4589.60) | 75338.80 (IQR 1749.20) | 6.35x | 6.899e-9 |

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
The Phase 18 release adds MtM-XCCS and tenor-basis analytics, CTD-COLVA/WWR, ZABR/PDV volatility, KIKO/window barriers, Everest and capped timer options, weather degree-day and clean spark/dark spreads, minimum-torsion ENB and entropy pooling. Phase 19 adds burnout-MBS/OAS, CMCDS, Quintic-OU VIX, AFV convertibles, Clewlow–Strickland swaptions, Schur-HRP and nonlinear shrinkage. Phase 20 adds Sinh-OU/CIR++ RFR, callable BK FRNs, no-Armageddon credit options, WWM-XVA, mixed rough Bergomi, Choi/Li–Deng–Zhou spreads, MRJD commodities, AFS LOB, EXVaR and PCRP. The workstation now exposes Tab 31 and the benchmark suite contains 33 AD parity workloads.
