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
| `ql-benches` | `crates/ql-benches` | Delivery | 61-suite parity/timing harness |
| `ql-terminal` | `crates/ql-terminal` | Delivery | Unified SDK, 59-tab Braille/ANSI workstation |

## Numerical capabilities

The workspace includes multi-curve OIS/IBOR bootstrapping, Black/Bachelier and arbitrage-aware volatility surfaces, Heston/rough-Heston and Lévy COS engines, ADI PDEs, SIMD Monte Carlo and American MC, scripted payoffs, XVA/SIMM/DIM/MVA, credit/inflation/FX/commodity analytics, barriers/Asians/lookbacks, callable and convertible bonds, and a pure-Rust ANSI/Braille risk workstation.

## Empirical performance

The offline benchmark harness checks numerical parity against independent baselines. Results below are generated from the current release tree; run the benchmark command on your hardware for fresh timings.

**Benchmark summary:** 61 suites and 59 workstation tabs.

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 966334.80 (IQR 891727.40) | 5956.40 (IQR 208.40) | 162.23x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 76079.60 (IQR 8897.60) | 4843.00 (IQR 32.80) | 15.71x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 1126600.40 (IQR 661851.40) | 21265.40 (IQR 2058.60) | 52.98x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 317553.20 (IQR 65099.80) | 125965.60 (IQR 13000.20) | 2.52x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1915666.60 (IQR 814161.40) | 648402.00 (IQR 1087.40) | 2.95x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 428572.80 (IQR 7765.20) | 149475.40 (IQR 3099.60) | 2.87x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 787526.80 (IQR 5891.60) | 172510.20 (IQR 5323.20) | 4.57x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 3067289.00 (IQR 159755.60) | 318835.00 (IQR 114109.00) | 9.62x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 390470.40 (IQR 4610.40) | 17106.80 (IQR 482.80) | 22.83x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 62076.80 (IQR 971.40) | 21561.20 (IQR 41.40) | 2.88x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19695584.60 (IQR 666149.60) | 203768.20 (IQR 12932.60) | 96.66x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1766374.00 (IQR 15408.40) | 213995.00 (IQR 1378.80) | 8.25x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37317.60 (IQR 3225.40) | 12262.20 (IQR 830.60) | 3.04x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13528219.20 (IQR 22654.60) | 1193260.20 (IQR 3849.20) | 11.34x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 569719.40 (IQR 4388.20) | 130185.80 (IQR 463.00) | 4.38x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 206980.40 (IQR 225.40) | 54979.20 (IQR 885.00) | 3.76x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28791.60 (IQR 301.20) | 5929.80 (IQR 404.80) | 4.86x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1920838.40 (IQR 8992.80) | 11453.80 (IQR 433.20) | 167.70x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 520171.00 (IQR 815.80) | 55450.80 (IQR 199.40) | 9.38x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7156231.40 (IQR 62067.00) | 156280.60 (IQR 7790.60) | 45.79x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 572791.00 (IQR 11314.20) | 273346.20 (IQR 1252.60) | 2.10x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 535108.20 (IQR 4896.40) | 8761.80 (IQR 891.40) | 61.07x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 1052256.80 (IQR 32111.60) | 135978.20 (IQR 3829.40) | 7.74x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6993.80 (IQR 273.80) | 1476.00 (IQR 159.20) | 4.74x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 314285.60 (IQR 663.20) | 47596.60 (IQR 3026.40) | 6.60x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 333445.60 (IQR 1282.40) | 71258.00 (IQR 2584.80) | 4.68x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 134136.00 (IQR 3272.40) | 62525.20 (IQR 1250.60) | 2.15x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3812607.00 (IQR 1785.00) | 436981.60 (IQR 539.20) | 8.72x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 82305.60 (IQR 5558.20) | 17332.80 (IQR 1165.00) | 4.75x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 311520.80 (IQR 3311.60) | 38234.40 (IQR 404.40) | 8.15x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 969805.20 (IQR 771.80) | 138130.00 (IQR 4378.60) | 7.02x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1640143.60 (IQR 41505.80) | 242151.80 (IQR 13939.40) | 6.77x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Monotone-Convex / Clark-FX / SA-CCR-KVA / rVol / Cartea / Paleologo AD | 477533.00 (IQR 703.00) | 75524.40 (IQR 2010.60) | 6.32x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | KRD / LPI / CoCo / Bergomi / CGMY / Guéant / Meucci AD vs FD | 2633486.40 (IQR 83954.20) | 695046.60 (IQR 95689.40) | 3.79x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | RKHS / FMM / VV / LSV / JLT / CTD / QRH / Parisian / Aïd / Propagator / NCO / DML AD | 99691.60 (IQR 764.00) | 27237.40 (IQR 1520.80) | 3.66x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1394582.00 (IQR 16468.00) | 170405.20 (IQR 12217.20) | 8.18x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Kalman-RTS / VECM / HJ-SDF / IFRS9 / CPT / Hasbrouck / FV / Lando / LPPL / CB-Arb / Causal / Malliavin AD | 139576.60 (IQR 1035.60) | 17581.80 (IQR 2492.40) | 7.94x | 8.967e-6 |
| 38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc | MF-Rates / Wildcard / CTD / Pivot-TARF / Bergomi / Schönbucher / FRTB-CVA / Lookback / VTA / FLAM / ACD / HLZ / MLMC AD | 1123536.40 (IQR 15299.00) | 49194.60 (IQR 936.60) | 22.84x | 4.996e-9 |
| 39_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt | FTPL / G2++ / eSSVI / Offshore-CB / CoCo / MPOR / Double-Heston / Curran / Warin / Roll-Corwin / Busseti-Kelly / Diff-QP AD | 676104.80 (IQR 58756.80) | 48116.40 (IQR 4917.60) | 14.05x | 4.423e-10 |
| 40_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca | BCBS-IRRBB / Turfus-RFR / Baran-MtM / Bourgey-BG / FEDS-PONV / Zeron-SIMM / Jacquier-rBergomi / Choi-Spread / BESS-SDP / Kratz-Dark / Föllmer / Diff-PCA AD | 348214.00 (IQR 127329.80) | 26991.20 (IQR 3282.80) | 12.90x | 1.803e-9 |
| 41_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is | Veronesi / Lando / Gatheral-Corridor / Russian / Bates-SVJ / Koziol-CB / MFG / BCBS-Floor / CJP-LOB / ETF / Qian-IC / Glasserman-IS AD | 824708.80 (IQR 158672.20) | 85931.80 (IQR 1453.40) | 9.60x | 2.244e-9 |
| 42_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin | Romer-Taylor / Damodaran-CRP / Cumulant-VaR / LOB-LVaR / Prospect / Regret / Beneish / Campisi / Heckman / Elastic-Net / Reset-CB / FSB-BailIn AD | 21450.00 (IQR 34.60) | 2426.80 (IQR 105.80) | 8.84x | 8.962e-10 |
| 43_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi | KRD-LSC / CMO-PAC / Asquith-CB / ASU2024-CPDI / ASC470-Mod / ColVA-CTD / Hill-EVT / Bouchaud-TLF / Karnosky-FX / DGTW-Timing / Brock-Hommes / Penman-ReOI AD | 61211.60 (IQR 402.00) | 5519.20 (IQR 353.40) | 11.09x | 1.010e-9 |
| 44_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan | Cheyette-LV / ADB-NDIRS / Andreasen-ZABR / QRH-Zumbach / CMCDS-Morini / Turlakov-WWR / HHW-ADI / CappedCall-MW / Ahmadi-EVaR / Schur-HRP / BEKK-GARCH / Schilit-Sloan AD | 42137.60 (IQR 597.40) | 4412.20 (IQR 98.40) | 9.55x | 5.153e-8 |
| 45_rebonato_cms_leland_toft_calamos_ascot_spec_pool_roll_glosten_milgrom_potters_rie_egarch_hamilton_wooldridge_2sls_mcneil_euler_bacon_omega_damodaran_rd_mwg_arrow_debreu | Rebonato-CMS / Leland-Toft / Calamos-Ascot / SpecPool-Roll / Glosten-Milgrom / Potters-RIE / EGARCH-Hamilton / Wooldridge-2SLS / McNeil-Euler / Bacon-Omega / Damodaran-RD / Arrow-Debreu AD | 14675.00 (IQR 22.80) | 1697.00 (IQR 73.80) | 8.65x | 6.115e-9 |
| 46_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig | Tuckman-Repo / BK-FVA / DECS / MBS-IO-PO / Kyle-PIN / GP-Aim / GJR-DCC / LPPL / Spectral-Risk / Brinson-Carino / Prelec-CPT / Path-Sig AD | 34710.20 (IQR 372.60) | 2698.60 (IQR 108.20) | 12.86x | 5.156e-9 |
| 47_hagan_west_mc_fabozzi_ctd_switch_afv_default_cb_lando_kth_default_gregory_sa_cva_gibson_schwartz_kirk_clark_fx_vv_bouzoubaa_napoleon_grinold_tc_flam_meucci_entropy_hastie_krr_loocv_koller_mckinsey_eva | Hagan-West-MC / CTD-Switch / AFV-CB / Kth-Default / SA-CVA-Euler / GS-Kirk-Spread / FX-Vanna-Volga / Napoleon-Cliquet / Transfer-Coef / Meucci-Entropy / KRR-LOOCV / McKinsey-EVA AD | 38564.80 (IQR 321.80) | 4358.20 (IQR 439.40) | 8.85x | 2.045e-8 |
| 48_brigo_csa_ctd_fabozzi_callable_oas_lando_jlt_migration_gatheral_corridor_varswap_glasserman_lsm_rogers_dual_cartea_lob_mm_lo_mackinlay_vr_lopez_dsr_purged_cv_roncalli_spinu_erc_ilmanen_vmc_penman_reoi_mccall_bellman | Multi-CSA-CTD / Callable-OAS / JLT-Migration / Corridor-VarSwap / LSM-Rogers-Dual / Cartea-LOB-MM / Lo-MacKinlay-VR / Lopez-DSR-PurgedCV / Roncalli-ERC / Ilmanen-VMC / Penman-ReOI / McCall-Bellman AD | 73776.00 (IQR 265.00) | 7511.60 (IQR 291.60) | 9.82x | 4.941e-8 |
| 49_brigo_nss_krd_tuckman_tips_floor_schonbucher_rmv_rt_rfv_gatheral_sticky_lv_gregory_simm_ccp_mva_schwartz_smith_kalman_wystup_fx_dnt_bouzoubaa_timer_meucci_bl_woodbury_gueant_ac_sinh_zou_elastic_net_nk_ree | NSS-KRD / TIPS-Floor-BEI / RMV-RT-RFV / Sticky-LocalVol / SIMM-CCP-MVA / Schwartz-Smith-KF / FX-DNT-VV / Timer-VarClock / BL-Woodbury / AC-Sinh-Liq / ElasticNet-Group / NK-DSGE-REE AD | 40047.60 (IQR 345.40) | 3222.40 (IQR 127.20) | 12.43x | 1.279e-9 |
| 50_hagan_west_tension_duffie_lando_pac_cmo_fama_macbeth_shanken_gatheral_asvi_bgk_discrete_barrier_bouchaud_cross_impact_karnosky_singer_merton_jump_campbell_viceira_ohlson_lid_griewank_revolve | Hagan-West-Tension / Duffie-Lando / PAC-CMO / Fama-MacBeth-Shanken / Gatheral-aSVI / BGK-Barrier / Cross-Impact / Karnosky-Singer / Merton-Jump / Campbell-Viceira / Ohlson-LID / Revolve-BSB AD | 38658.20 (IQR 53.20) | 6888.00 (IQR 103.20) | 5.61x | 1.448e-8 |
| 51_brigo_xccy_mtm_fabozzi_spec_pool_io_po_lando_marshall_olkin_gatheral_vix_vvix_gregory_wwr_cva_geman_gibson_schwartz_clark_fx_fva_rubinstein_chooser_grinold_flam_tc_meucci_es_euler_james_stein_sure_campbell_cochrane_habit | XCCY-MtM-Basis / Spec-Pool-IO-PO / Marshall-Olkin / VIX-VVIX-Smile / WWR-CVA / Gibson-Schwartz / FX-FVA / Chooser / FLAM-TC / ES-Euler / James-Stein-SURE / Campbell-Cochrane AD | 7869.40 (IQR 286.80) | 1160.60 (IQR 78.40) | 6.78x | 2.010e-8 |
| 52_brigo_ois_fallback_fabozzi_collared_frn_schonbucher_cds_swaption_fep_bergomi_2f_skew_gregory_sa_ccr_ead_eydeland_spark_spread_wystup_prdc_geske_compound_qian_alpha_decay_roncalli_enb_johansen_vecm_bansal_yaron_lrr | OIS-Fallback / Collared-FRN / CDS-Swaption-FEP / Bergomi-2F-Skew / SA-CCR-EAD / Clean-Spark / FX-PRDC / Geske-Compound / Alpha-Decay-IR / Meucci-ENB / Johansen-VECM / Bansal-Yaron-LRR AD | 29544.80 (IQR 25.60) | 3174.00 (IQR 42.20) | 9.31x | 1.469e-9 |
| 53_andersen_cms_replication_tuckman_special_repo_okane_cds_upfront_gatheral_ssr_tanaka_gregory_saccr_kva_eydeland_crack_spread_wystup_prdc_cap_haug_eso_dilution_qian_spinu_erc_rockafellar_cvar_htf_pls_krylov_piazzesi_fomc | CMS-Replication / Special-Repo-CTD / ISDA-CDS-Upfront / SSR-Tanaka / SA-CCR-KVA / 3-Asset-Crack / PRDC-Cap / ESO-Dilution / Spinu-ERC / Rockafellar-CVaR / PLS-Krylov / FOMC-Hike AD | 19984.40 (IQR 10.80) | 1884.60 (IQR 257.60) | 10.60x | 2.693e-9 |
| 54_darbyshire_step_ois_turn_despiegeleer_cb_soft_call_okane_cds_index_skew_gatheral_varswap_truncation_gregory_fva_overlap_geman_working_storage_wystup_malz_fx_fly_goldman_sosin_lookback_ledoit_wolf_shrinkage_almgren_chriss_cara_tsay_estar_acm_cp | Step-OIS-Turn / CB-SoftCall-MW / CDX-Skew-Basis / VarSwap-Trunc / FVA-DVA-Overlap / Working-Storage / FX-Malz-Fly / GSG-Lookback / Ledoit-Wolf-CC / Almgren-CARA / ESTAR-Band / ACM-CP-TP AD | 13042.00 (IQR 6.80) | 1213.80 (IQR 38.80) | 10.74x | 1.054e-8 |
| 55_henrard_tenor_basis_fabozzi_msr_recapture_okane_cdo_base_corr_gatheral_fwd_start_smile_gregory_pykhtin_csa_mta_eydeland_gas_swing_wystup_fx_kiko_stulz_rainbow_jagannathan_ma_obizhaeva_wang_newey_west_hac_gabaix_disaster | Tenor-Basis-Mult / MSR-Recapture / CDO-Base-Corr / Fwd-Start-Smile / CSA-MTA-CVA / Gas-Swing / FX-Window-KIKO / Stulz-Rainbow / Jagannathan-Ma / Obizhaeva-Wang / Newey-West-HAC / Gabaix-Disaster AD | 106105.20 (IQR 300.60) | 9357.00 (IQR 397.40) | 11.34x | 3.026e-7 |
| 56_brigo_diff_swap_quanto_fabozzi_callable_muni_ytw_bielecki_step_up_bond_bkm_moment_strip_gregory_break_clause_cva_bessembinder_power_skew_wystup_fx_fader_longstaff_extendible_kan_zhou_3fund_sharpe_rbsa_diebold_mariano_hln_laubach_williams_rstar | Diff-Swap-Quanto / Callable-Muni-YTW / Step-Up-Bond / BKM-SKEW / Break-Clause-CVA / Power-Skew / FX-Fader / Longstaff-Extendible / Kan-Zhou-3Fund / Sharpe-RBSA / Diebold-Mariano-HLN / Laubach-Williams-r* AD | 60795.20 (IQR 303.20) | 7320.40 (IQR 51.00) | 8.30x | 8.489e-9 |
| 57_andersen_sofr_futures_fra_hayre_mbs_scurve_pedersen_cdx_fep_gatheral_dupire_total_var_gregory_simm_mva_routledge_stockout_baba_cip_basis_carr_chou_symmetry_pastor_stambaugh_treynor_black_engle_granger_ecm_epstein_zin_sdf | SOFR-Fut-FRA / Hayre-MBS-SCurve / CDX-Option-FEP / Dupire-TotalVar-LV / SIMM-MVA-Funding / RSS-Stockout / CIP-Basis-Wedge / Carr-Chou-Symmetry / Pastor-Stambaugh / Treynor-Black / Engle-Granger-ECM / Epstein-Zin-SDF AD | 6672.80 (IQR 18.80) | 1005.20 (IQR 34.00) | 6.64x | 8.576e-7 |
| 58_litterman_scheinkman_pca_fly_tf_split_cb_hull_white_double_t_cdo_roger_lee_moment_wings_brigo_capponi_cds_contagion_schwartz_ou_futures_wystup_dcd_yield_passport_option_grs_factor_ftest_bouchaud_sqrt_impact_bns_bipower_har_cj_lettau_ludvigson_cay | PCA-Fly-Carry / TF-Split-CB / Double-t-CDO / Roger-Lee-Wings / CDS-Contagion-CVA / Schwartz-OU-Fut / FX-DCD-Yield / Passport-Opt / GRS-Factor-FTest / Sqrt-Metaorder / BNS-HAR-CJ / Lettau-Ludvigson-cay AD | 9925.80 (IQR 11.20) | 1017.80 (IQR 34.40) | 9.75x | 3.312e-9 |
| 59_piterbarg_collateral_option_tuckman_rolldown_merton_kmv_dd_hagan_bartlett_sabr_pykhtin_mpor_spike_deng_locational_ftr_wystup_digital_overhedge_conze_fixed_lookback_tasche_hot_spots_es_frazzini_bab_pesaran_timmermann_lucas_mehra_prescott | Piterbarg-Coll-Opt / Tuckman-Rolldown / Merton-KMV-DD / Bartlett-SABR / Pykhtin-MPoR-Spike / Deng-FTR-Opt / Wystup-FX-Digital / Conze-Fixed-Lookback / Tasche-Hot-Spots-ES / Frazzini-BAB / Pesaran-Timmermann / Lucas-Mehra-Prescott AD | 10397.40 (IQR 7.80) | 1589.60 (IQR 32.20) | 6.54x | 6.269e-8 |
| 60_gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_ga_carr_wu_vol_drift_gregory_marginal_cva_working_brennan_storage_wystup_forward_extra_rubinstein_exchange_opt_hlz_factor_fdr_holden_subrahmanyam_kyle_bollerslev_garch_kurtosis_campbell_shiller_cape | GSW-Svensson-Curve / Fabozzi-MBS-ESF / Vasicek-Gordy-IRB-GA / Carr-Wu-Vol-Drift / Gregory-Marginal-CVA / Working-Brennan-Storage / Wystup-Forward-Extra / Rubinstein-Exchange-Opt / HLZ-Factor-FDR / Holden-Subrahmanyam-Kyle / Bollerslev-GARCH-Kurtosis / Campbell-Shiller-CAPE AD | 35897.00 (IQR 309.40) | 2669.60 (IQR 53.20) | 13.45x | 1.252e-8 |
| 61_kim_wright_shadow_rate_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_swap_burgard_kjaer_fva_carmona_durrleman_spark_wystup_fx_tarf_bgk_discrete_lookback_novy_marx_ff5_acerbi_spectral_es_corsi_har_rv_bll_news_noise | Kim-Wright-Shadow / Brennan-Schwartz-CB / Duffie-Singleton-RMV / Fukasawa-Gamma-Swap / Burgard-Kjaer-FVA / Carmona-Durrleman-Spark / Wystup-FX-TARF / BGK-Discrete-Lookback / Novy-Marx-FF5 / Acerbi-Spectral-ES / Corsi-HAR-RV / BLL-News-Noise AD | 27407.20 (IQR 5.60) | 1885.40 (IQR 54.80) | 14.54x | 2.428e-9 |

Summary: 61 benchmark suites; speedups range from **2.10x** to **167.70x** (median **8.65x**, mean **17.11x**). All parity checks pass their domain tolerances; latest suite parity is **2.428e-9**.

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

