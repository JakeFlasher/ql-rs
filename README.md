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
| `ql-benches` | `crates/ql-benches` | Delivery | 46-suite parity/timing harness |
| `ql-terminal` | `crates/ql-terminal` | Delivery | Unified SDK, 44-tab Braille/ANSI workstation |

## Numerical capabilities

The workspace includes multi-curve OIS/IBOR bootstrapping, Black/Bachelier and arbitrage-aware volatility surfaces, Heston/rough-Heston and Lévy COS engines, ADI PDEs, SIMD Monte Carlo and American MC, scripted payoffs, XVA/SIMM/DIM/MVA, credit/inflation/FX/commodity analytics, barriers/Asians/lookbacks, callable and convertible bonds, and a pure-Rust ANSI/Braille risk workstation.

## Empirical performance

The offline benchmark harness checks numerical parity against independent baselines. Results below are generated from the current release tree; run the benchmark command on your hardware for fresh timings.

**Benchmark summary:** 46 suites and 44 workstation tabs.






| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 582243.20 (IQR 462122.20) | 5279.60 (IQR 195.60) | 110.28x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 56008.60 (IQR 6681.80) | 3294.00 (IQR 55.20) | 17.00x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 257112.80 (IQR 54797.40) | 4630.80 (IQR 140.20) | 55.52x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 95727.20 (IQR 280.80) | 31064.60 (IQR 297.60) | 3.08x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1963339.60 (IQR 62410.60) | 616445.40 (IQR 11655.00) | 3.18x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 440098.00 (IQR 6409.40) | 152344.60 (IQR 366.80) | 2.89x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 807873.00 (IQR 5966.00) | 187768.40 (IQR 3743.00) | 4.30x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 2986548.60 (IQR 136240.20) | 310032.60 (IQR 65022.00) | 9.63x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 385055.40 (IQR 715.60) | 19297.80 (IQR 3152.60) | 19.95x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 61831.60 (IQR 481.40) | 21548.40 (IQR 298.60) | 2.87x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19650491.80 (IQR 51583.40) | 209615.80 (IQR 9552.60) | 93.75x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1798602.80 (IQR 45035.80) | 214342.80 (IQR 6762.40) | 8.39x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37312.20 (IQR 2607.80) | 12345.20 (IQR 1306.40) | 3.02x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13502688.40 (IQR 359369.80) | 1176804.00 (IQR 6091.20) | 11.47x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 537164.00 (IQR 4190.00) | 128114.40 (IQR 1874.20) | 4.19x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 203033.40 (IQR 353.00) | 54284.20 (IQR 420.80) | 3.74x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28200.60 (IQR 236.80) | 6134.00 (IQR 304.20) | 4.60x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1912853.80 (IQR 17488.40) | 11464.20 (IQR 220.40) | 166.85x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 516968.20 (IQR 4513.80) | 55586.40 (IQR 843.40) | 9.30x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7058782.40 (IQR 223659.20) | 149949.00 (IQR 3448.60) | 47.07x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 540005.40 (IQR 7523.80) | 257200.40 (IQR 788.00) | 2.10x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 507776.80 (IQR 4469.60) | 8791.20 (IQR 851.20) | 57.76x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 799033.60 (IQR 5762.20) | 134059.40 (IQR 4977.00) | 5.96x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6606.40 (IQR 502.60) | 1287.60 (IQR 141.60) | 5.13x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 300635.00 (IQR 3761.60) | 46053.00 (IQR 2022.60) | 6.53x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 321482.00 (IQR 2016.80) | 68934.60 (IQR 5711.80) | 4.66x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 131279.80 (IQR 1317.60) | 61075.20 (IQR 910.80) | 2.15x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3812477.80 (IQR 73167.00) | 440019.60 (IQR 3594.80) | 8.66x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 85183.60 (IQR 1370.20) | 17435.40 (IQR 1346.80) | 4.89x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 361256.60 (IQR 15612.40) | 45613.20 (IQR 21231.60) | 7.92x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 974568.00 (IQR 1177.40) | 140149.00 (IQR 13075.80) | 6.95x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1680715.60 (IQR 8712.80) | 231608.60 (IQR 1678.80) | 7.26x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Monotone-Convex / Clark-FX / SA-CCR-KVA / rVol / Cartea / Paleologo AD | 494515.00 (IQR 17783.80) | 77267.80 (IQR 1343.80) | 6.40x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | KRD / LPI / CoCo / Bergomi / CGMY / Guéant / Meucci AD vs FD | 2631196.80 (IQR 54068.80) | 676246.80 (IQR 6257.00) | 3.89x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | RKHS / FMM / VV / LSV / JLT / CTD / QRH / Parisian / Aïd / Propagator / NCO / DML AD | 93230.60 (IQR 1787.60) | 27349.80 (IQR 10815.60) | 3.41x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1342894.20 (IQR 37386.60) | 155550.80 (IQR 2697.40) | 8.63x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Kalman-RTS / VECM / HJ-SDF / IFRS9 / CPT / Hasbrouck / FV / Lando / LPPL / CB-Arb / Causal / Malliavin AD | 137118.40 (IQR 10306.40) | 15212.20 (IQR 196.60) | 9.01x | 8.967e-6 |
| 38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc | MF-Rates / Wildcard / CTD / Pivot-TARF / Bergomi / Schönbucher / FRTB-CVA / Lookback / VTA / FLAM / ACD / HLZ / MLMC AD | 1169043.20 (IQR 10454.20) | 50913.40 (IQR 1043.00) | 22.96x | 4.996e-9 |
| 39_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt | FTPL / G2++ / eSSVI / Offshore-CB / CoCo / MPOR / Double-Heston / Curran / Warin / Roll-Corwin / Busseti-Kelly / Diff-QP AD | 694212.40 (IQR 9434.60) | 38285.40 (IQR 5471.80) | 18.13x | 4.423e-10 |
| 40_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca | BCBS-IRRBB / Turfus-RFR / Baran-MtM / Bourgey-BG / FEDS-PONV / Zeron-SIMM / Jacquier-rBergomi / Choi-Spread / BESS-SDP / Kratz-Dark / Föllmer / Diff-PCA AD | 278439.20 (IQR 604.60) | 25230.40 (IQR 3845.60) | 11.04x | 1.803e-9 |
| 41_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is | Veronesi / Lando / Gatheral-Corridor / Russian / Bates-SVJ / Koziol-CB / MFG / BCBS-Floor / CJP-LOB / ETF / Qian-IC / Glasserman-IS AD | 716226.40 (IQR 362.80) | 86492.40 (IQR 905.00) | 8.28x | 2.244e-9 |
| 42_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin | Romer-Taylor / Damodaran-CRP / Cumulant-VaR / LOB-LVaR / Prospect / Regret / Beneish / Campisi / Heckman / Elastic-Net / Reset-CB / FSB-BailIn AD | 21929.40 (IQR 56.60) | 2452.20 (IQR 146.60) | 8.94x | 8.962e-10 |
| 43_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi | KRD-LSC / CMO-PAC / Asquith-CB / ASU2024-CPDI / ASC470-Mod / ColVA-CTD / Hill-EVT / Bouchaud-TLF / Karnosky-FX / DGTW-Timing / Brock-Hommes / Penman-ReOI AD | 62041.40 (IQR 48.60) | 5580.20 (IQR 695.60) | 11.12x | 1.010e-9 |
| 44_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan | Cheyette-LV / ADB-NDIRS / Andreasen-ZABR / QRH-Zumbach / CMCDS-Morini / Turlakov-WWR / HHW-ADI / CappedCall-MW / Ahmadi-EVaR / Schur-HRP / BEKK-GARCH / Schilit-Sloan AD | 42730.20 (IQR 35.40) | 4533.40 (IQR 95.40) | 9.43x | 5.153e-8 |
| 45_rebonato_cms_leland_toft_calamos_ascot_spec_pool_roll_glosten_milgrom_potters_rie_egarch_hamilton_wooldridge_2sls_mcneil_euler_bacon_omega_damodaran_rd_mwg_arrow_debreu | Rebonato-CMS / Leland-Toft / Calamos-Ascot / SpecPool-Roll / Glosten-Milgrom / Potters-RIE / EGARCH-Hamilton / Wooldridge-2SLS / McNeil-Euler / Bacon-Omega / Damodaran-RD / Arrow-Debreu AD | 15015.60 (IQR 313.40) | 1747.20 (IQR 64.00) | 8.59x | 6.115e-9 |
| 46_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig | Tuckman-Repo / BK-FVA / DECS / MBS-IO-PO / Kyle-PIN / GP-Aim / GJR-DCC / LPPL / Spectral-Risk / Brinson-Carino / Prelec-CPT / Path-Sig AD | 1200.60 (IQR 51.60) | 2727.20 (IQR 53.20) | 0.44x | 0.000e0 |

Summary: 46 benchmark suites; speedups range from **0.44x** to **166.85x** (median **8.28x**, mean **18.07x**). All parity checks pass their domain tolerances; latest suite parity is **0.000e0**.

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


