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
| `ql-benches` | `crates/ql-benches` | Delivery | 41-suite parity/timing harness |
| `ql-terminal` | `crates/ql-terminal` | Delivery | Unified SDK, 39-tab Braille/ANSI workstation |

## Numerical capabilities

The workspace includes multi-curve OIS/IBOR bootstrapping, Black/Bachelier and arbitrage-aware volatility surfaces, Heston/rough-Heston and Lévy COS engines, ADI PDEs, SIMD Monte Carlo and American MC, scripted payoffs, XVA/SIMM/DIM/MVA, credit/inflation/FX/commodity analytics, barriers/Asians/lookbacks, callable and convertible bonds, and a pure-Rust ANSI/Braille risk workstation.

## Empirical performance

The offline benchmark harness checks numerical parity against independent baselines. Results below are generated from the current release tree; run the benchmark command on your hardware for fresh timings.

**Benchmark summary:** 41 suites and 39 workstation tabs.






| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 474938.20 (IQR 489367.00) | 18453.00 (IQR 5428.60) | 25.74x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 183074.20 (IQR 24027.00) | 11289.60 (IQR 134.60) | 16.22x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 758328.00 (IQR 676092.20) | 8681.40 (IQR 919.60) | 87.35x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 137953.40 (IQR 28893.80) | 35141.80 (IQR 3208.00) | 3.93x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1915313.20 (IQR 13676.00) | 648495.40 (IQR 1043.60) | 2.95x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 426783.80 (IQR 8887.00) | 149124.80 (IQR 244.40) | 2.86x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 806542.40 (IQR 6986.40) | 197525.60 (IQR 4885.40) | 4.08x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 2996492.60 (IQR 219969.20) | 295714.80 (IQR 78185.00) | 10.13x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 385082.60 (IQR 438.60) | 17223.00 (IQR 658.40) | 22.36x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 62201.40 (IQR 231.40) | 21733.80 (IQR 51.20) | 2.86x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19397238.60 (IQR 73700.40) | 209023.40 (IQR 13650.00) | 92.80x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1805604.20 (IQR 49643.00) | 217482.40 (IQR 3785.60) | 8.30x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 35617.40 (IQR 2250.60) | 10934.60 (IQR 478.20) | 3.26x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13289580.40 (IQR 240043.00) | 1177184.60 (IQR 12172.80) | 11.29x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 537094.00 (IQR 10178.80) | 126956.20 (IQR 1712.60) | 4.23x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 203348.20 (IQR 3629.80) | 53788.80 (IQR 3388.20) | 3.78x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28499.40 (IQR 667.00) | 5898.60 (IQR 184.60) | 4.83x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1900772.00 (IQR 16652.00) | 11479.80 (IQR 443.60) | 165.58x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 511181.80 (IQR 8724.40) | 54188.20 (IQR 432.60) | 9.43x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7523595.00 (IQR 306894.80) | 150494.20 (IQR 3211.40) | 49.99x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 548670.80 (IQR 1621.00) | 261143.20 (IQR 2463.60) | 2.10x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 520367.40 (IQR 7457.00) | 8183.80 (IQR 979.80) | 63.59x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 899665.60 (IQR 22742.20) | 147380.80 (IQR 3055.00) | 6.10x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6597.00 (IQR 371.60) | 1348.20 (IQR 121.40) | 4.89x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 306667.40 (IQR 740.20) | 46886.80 (IQR 1614.00) | 6.54x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 326667.40 (IQR 184.20) | 69407.60 (IQR 5721.20) | 4.71x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 131382.60 (IQR 830.40) | 62266.40 (IQR 1910.00) | 2.11x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3749798.60 (IQR 3925.60) | 430301.20 (IQR 3216.60) | 8.71x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 85815.40 (IQR 2432.40) | 16704.00 (IQR 1171.00) | 5.14x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 353395.00 (IQR 2814.20) | 42793.60 (IQR 438.20) | 8.26x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 973801.20 (IQR 11813.20) | 138255.00 (IQR 5197.80) | 7.04x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1657821.00 (IQR 26807.60) | 225507.20 (IQR 1673.40) | 7.35x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Phase21 institutional 12-engine AAD vs FD | 479141.00 (IQR 351.40) | 74864.80 (IQR 1657.00) | 6.40x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | Phase22 KRD/LPI/CoCo/Bergomi/CGMY/Guéant/Meucci/AAD | 2553658.40 (IQR 9134.80) | 666176.20 (IQR 5685.00) | 3.83x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | Phase23 RKHS/FMM/VV/LSV/JLT/CTD/QRH/Parisian/Aïd/Propagator/NCO/DML | 93734.00 (IQR 600.60) | 25611.60 (IQR 3599.80) | 3.66x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1346057.60 (IQR 8605.20) | 156660.80 (IQR 1533.40) | 8.59x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Phase25 Kalman/RTS/VECM/HJ-SDF/IFRS9/CPT/Hasbrouck/FV/Lando/LPPL/CB-Arb/Causal/Malliavin | 136169.00 (IQR 163.80) | 16530.20 (IQR 1709.60) | 8.24x | 8.967e-6 |
| 38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc | Phase26 MF-Rates/Wildcard/CTD/Pivot-TARF/Bergomi/Schönbucher/FRTB-CVA/Lookback/VTA/FLAM/ACD/HLZ/MLMC | 1105841.60 (IQR 3733.60) | 47897.40 (IQR 1076.00) | 23.09x | 4.996e-9 |
| 39_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt | Phase27 SOTA FTPL/G2/eSSVI/Offshore-CB/CoCo/MPOR/Double-Heston/Curran/Warin/Roll-Corwin/Busseti-Kelly/Diff-QP | 658880.20 (IQR 2968.40) | 36847.20 (IQR 2655.40) | 17.88x | 4.423e-10 |
| 40_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca | Phase28 BCBS/IRRBB/Turfus-RFR/Baran-MtM/Bourgey-BG/FEDS-PONV/Zeron-SIMM/Jacquier-rBergomi/Choi-Spread/BESS-SDP/Kratz-Dark/Föllmer/Diff-PCA | 269642.40 (IQR 641.00) | 29391.20 (IQR 3427.60) | 9.17x | 1.803e-9 |
| 41_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is | Phase29 Veronesi/Lando/Gatheral/Russian/Bates/Koziol/MFG/BCBS/CJP/ETF/Qian/Glasserman | 28421.40 (IQR 117.80) | 3942.00 (IQR 196.60) | 7.21x | 9.818e-11 |

Summary: 41 benchmark suites; speedups range from **2.10x** to **165.58x** (median **7.21x**, mean **18.21x**). All parity checks pass their domain tolerances; latest suite parity is **9.818e-11**.

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
The Phase 18 release adds MtM-XCCS and tenor-basis analytics, CTD-COLVA/WWR, ZABR/PDV volatility, KIKO/window barriers, Everest and capped timer options, weather degree-day and clean spark/dark spreads, minimum-torsion ENB and entropy pooling. Phase 19 adds burnout-MBS/OAS, CMCDS, Quintic-OU VIX, AFV convertibles, Clewlow–Strickland swaptions, Schur-HRP and nonlinear shrinkage. Phase 20 adds Sinh-OU/CIR++ RFR, callable BK FRNs, no-Armageddon credit options, WWM-XVA, mixed rough Bergomi, Choi/Li–Deng–Zhou spreads, MRJD commodities, AFS LOB, EXVaR and PCRP. The workstation now exposes Tab 39 and the benchmark suite contains 41 AD parity workloads, including book-backed Phase 25 SOTA engines.

Phase 27 adds twelve book-backed SOTA engines and Suite 39; workstation Tab 37 presents FTPL, G2++, eSSVI, Offshore-CB, CoCo, MPOR, Double-Heston, Curran, Warin, Roll-Corwin, Busseti-Kelly, and Diff-QP telemetry (`BOTTOM_ROWS = 34`).

Phase 28 adds twelve book-backed SOTA engines in Suite 40; workstation Tab 38 (`38:P28-SOTA`) presents IRRBB, Turfus-RFR, Baran-MtM, Bourgey-BG, FEDS-PONV, Zeron-SIMM, Jacquier-rBergomi, Choi-Spread, BESS-SDP, Kratz-Dark, Föllmer diagnostics, and Diff-PCA telemetry (`BOTTOM_ROWS = 35`).


Phase 29 adds twelve book-backed SOTA quantitative finance engines, Benchmark Suite 41, and workstation Tab 39 with `BOTTOM_ROWS = 36`.
