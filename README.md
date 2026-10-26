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
| `ql-benches` | `crates/ql-benches` | Delivery | 48-suite parity/timing harness |
| `ql-terminal` | `crates/ql-terminal` | Delivery | Unified SDK, 46-tab Braille/ANSI workstation |

## Numerical capabilities

The workspace includes multi-curve OIS/IBOR bootstrapping, Black/Bachelier and arbitrage-aware volatility surfaces, Heston/rough-Heston and Lévy COS engines, ADI PDEs, SIMD Monte Carlo and American MC, scripted payoffs, XVA/SIMM/DIM/MVA, credit/inflation/FX/commodity analytics, barriers/Asians/lookbacks, callable and convertible bonds, and a pure-Rust ANSI/Braille risk workstation.

## Empirical performance

The offline benchmark harness checks numerical parity against independent baselines. Results below are generated from the current release tree; run the benchmark command on your hardware for fresh timings.

**Benchmark summary:** 48 suites and 46 workstation tabs.

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 306842.80 (IQR 42983.80) | 5207.00 (IQR 199.60) | 58.93x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 69678.60 (IQR 3643.60) | 4212.40 (IQR 133.20) | 16.54x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 492305.20 (IQR 183064.40) | 6914.80 (IQR 449.00) | 71.20x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 117290.00 (IQR 14471.60) | 32931.80 (IQR 1888.80) | 3.56x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1970034.60 (IQR 47265.00) | 612526.00 (IQR 6333.40) | 3.22x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 427237.20 (IQR 6678.60) | 149480.00 (IQR 382.60) | 2.86x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 808018.80 (IQR 16319.20) | 175336.00 (IQR 4139.00) | 4.61x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 3009438.40 (IQR 116022.60) | 376850.20 (IQR 77722.20) | 7.99x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 385018.20 (IQR 1879.00) | 16954.40 (IQR 2667.80) | 22.71x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 61914.20 (IQR 428.00) | 21651.60 (IQR 27.60) | 2.86x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19822078.80 (IQR 1350253.20) | 204269.40 (IQR 10298.60) | 97.04x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1774383.40 (IQR 10578.40) | 216599.80 (IQR 4168.00) | 8.19x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37295.40 (IQR 2507.20) | 13282.00 (IQR 1048.80) | 2.81x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13552864.40 (IQR 31958.60) | 1198719.60 (IQR 3553.00) | 11.31x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 562033.60 (IQR 10607.80) | 132451.00 (IQR 1400.60) | 4.24x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 210782.80 (IQR 282.40) | 55891.60 (IQR 788.40) | 3.77x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 29503.00 (IQR 346.40) | 6116.60 (IQR 209.60) | 4.82x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1927782.80 (IQR 43147.60) | 11532.80 (IQR 285.20) | 167.16x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 525350.20 (IQR 612.80) | 56663.60 (IQR 292.00) | 9.27x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7395442.20 (IQR 126192.80) | 157714.20 (IQR 5072.80) | 46.89x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 545854.80 (IQR 957.20) | 261798.60 (IQR 748.00) | 2.09x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 514043.60 (IQR 6796.80) | 8283.00 (IQR 792.20) | 62.06x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 1207612.80 (IQR 12186.00) | 134330.60 (IQR 3348.60) | 8.99x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6736.60 (IQR 1090.60) | 1306.40 (IQR 221.20) | 5.16x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 290917.60 (IQR 4818.20) | 42635.80 (IQR 1959.20) | 6.82x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 326115.80 (IQR 223.20) | 69748.60 (IQR 2503.20) | 4.68x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 133870.80 (IQR 267.40) | 65608.80 (IQR 428.20) | 2.04x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3813126.00 (IQR 16403.20) | 437064.00 (IQR 843.00) | 8.72x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 87386.20 (IQR 4491.40) | 17416.00 (IQR 1217.40) | 5.02x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 357670.20 (IQR 719.40) | 43664.60 (IQR 419.60) | 8.19x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 974783.40 (IQR 1627.60) | 141169.20 (IQR 8764.20) | 6.91x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1685686.20 (IQR 35505.20) | 231368.80 (IQR 1658.20) | 7.29x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Monotone-Convex / Clark-FX / SA-CCR-KVA / rVol / Cartea / Paleologo AD | 487709.80 (IQR 1487.60) | 76450.60 (IQR 1143.20) | 6.38x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | KRD / LPI / CoCo / Bergomi / CGMY / Guéant / Meucci AD vs FD | 2601388.20 (IQR 7408.20) | 674537.20 (IQR 6464.60) | 3.86x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | RKHS / FMM / VV / LSV / JLT / CTD / QRH / Parisian / Aïd / Propagator / NCO / DML AD | 96854.60 (IQR 2533.00) | 26392.20 (IQR 1694.80) | 3.67x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1356105.40 (IQR 16994.60) | 156120.60 (IQR 3168.60) | 8.69x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Kalman-RTS / VECM / HJ-SDF / IFRS9 / CPT / Hasbrouck / FV / Lando / LPPL / CB-Arb / Causal / Malliavin AD | 136234.80 (IQR 486.40) | 17331.00 (IQR 2207.60) | 7.86x | 8.967e-6 |
| 38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc | MF-Rates / Wildcard / CTD / Pivot-TARF / Bergomi / Schönbucher / FRTB-CVA / Lookback / VTA / FLAM / ACD / HLZ / MLMC AD | 1113467.40 (IQR 16078.20) | 48961.60 (IQR 1454.80) | 22.74x | 4.996e-9 |
| 39_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt | FTPL / G2++ / eSSVI / Offshore-CB / CoCo / MPOR / Double-Heston / Curran / Warin / Roll-Corwin / Busseti-Kelly / Diff-QP AD | 665655.60 (IQR 7259.40) | 37633.00 (IQR 4390.80) | 17.69x | 4.423e-10 |
| 40_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca | BCBS-IRRBB / Turfus-RFR / Baran-MtM / Bourgey-BG / FEDS-PONV / Zeron-SIMM / Jacquier-rBergomi / Choi-Spread / BESS-SDP / Kratz-Dark / Föllmer / Diff-PCA AD | 274228.80 (IQR 113.20) | 25832.20 (IQR 3815.00) | 10.62x | 1.803e-9 |
| 41_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is | Veronesi / Lando / Gatheral-Corridor / Russian / Bates-SVJ / Koziol-CB / MFG / BCBS-Floor / CJP-LOB / ETF / Qian-IC / Glasserman-IS AD | 697958.40 (IQR 1442.60) | 84786.80 (IQR 2657.80) | 8.23x | 2.244e-9 |
| 42_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin | Romer-Taylor / Damodaran-CRP / Cumulant-VaR / LOB-LVaR / Prospect / Regret / Beneish / Campisi / Heckman / Elastic-Net / Reset-CB / FSB-BailIn AD | 21568.40 (IQR 70.00) | 2398.60 (IQR 180.60) | 8.99x | 8.962e-10 |
| 43_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi | KRD-LSC / CMO-PAC / Asquith-CB / ASU2024-CPDI / ASC470-Mod / ColVA-CTD / Hill-EVT / Bouchaud-TLF / Karnosky-FX / DGTW-Timing / Brock-Hommes / Penman-ReOI AD | 60822.80 (IQR 415.80) | 5436.60 (IQR 362.60) | 11.19x | 1.010e-9 |
| 44_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan | Cheyette-LV / ADB-NDIRS / Andreasen-ZABR / QRH-Zumbach / CMCDS-Morini / Turlakov-WWR / HHW-ADI / CappedCall-MW / Ahmadi-EVaR / Schur-HRP / BEKK-GARCH / Schilit-Sloan AD | 41770.00 (IQR 306.00) | 4402.00 (IQR 98.60) | 9.49x | 5.153e-8 |
| 45_rebonato_cms_leland_toft_calamos_ascot_spec_pool_roll_glosten_milgrom_potters_rie_egarch_hamilton_wooldridge_2sls_mcneil_euler_bacon_omega_damodaran_rd_mwg_arrow_debreu | Rebonato-CMS / Leland-Toft / Calamos-Ascot / SpecPool-Roll / Glosten-Milgrom / Potters-RIE / EGARCH-Hamilton / Wooldridge-2SLS / McNeil-Euler / Bacon-Omega / Damodaran-RD / Arrow-Debreu AD | 14656.80 (IQR 57.80) | 1697.80 (IQR 79.60) | 8.63x | 6.115e-9 |
| 46_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig | Tuckman-Repo / BK-FVA / DECS / MBS-IO-PO / Kyle-PIN / GP-Aim / GJR-DCC / LPPL / Spectral-Risk / Brinson-Carino / Prelec-CPT / Path-Sig AD | 33970.00 (IQR 323.00) | 2707.80 (IQR 85.00) | 12.55x | 5.156e-9 |
| 47_hagan_west_mc_fabozzi_ctd_switch_afv_default_cb_lando_kth_default_gregory_sa_cva_gibson_schwartz_kirk_clark_fx_vv_bouzoubaa_napoleon_grinold_tc_flam_meucci_entropy_hastie_krr_loocv_koller_mckinsey_eva | Hagan-West-MC / CTD-Switch / AFV-CB / Kth-Default / SA-CVA-Euler / GS-Kirk-Spread / FX-Vanna-Volga / Napoleon-Cliquet / Transfer-Coef / Meucci-Entropy / KRR-LOOCV / McKinsey-EVA AD | 38419.80 (IQR 282.60) | 4285.00 (IQR 470.00) | 8.97x | 2.045e-8 |
| 48_brigo_csa_ctd_fabozzi_callable_oas_lando_jlt_migration_gatheral_corridor_varswap_glasserman_lsm_rogers_dual_cartea_lob_mm_lo_mackinlay_vr_lopez_dsr_purged_cv_roncalli_spinu_erc_ilmanen_vmc_penman_reoi_mccall_bellman | Multi-CSA-CTD / Callable-OAS / JLT-Migration / Corridor-VarSwap / LSM-Rogers-Dual / Cartea-LOB-MM / Lo-MacKinlay-VR / Lopez-DSR-PurgedCV / Roncalli-ERC / Ilmanen-VMC / Penman-ReOI / McCall-Bellman AD | 74023.80 (IQR 202.20) | 7684.20 (IQR 344.20) | 9.63x | 4.941e-8 |

Summary: 48 benchmark suites; speedups range from **2.04x** to **167.16x** (median **8.23x**, mean **17.44x**). All parity checks pass their domain tolerances; latest suite parity is **4.941e-8**.

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

