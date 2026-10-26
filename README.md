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
| `ql-benches` | `crates/ql-benches` | Delivery | 63-suite parity/timing harness |
| `ql-terminal` | `crates/ql-terminal` | Delivery | Unified SDK, 61-tab Braille/ANSI workstation |

## Numerical capabilities

The workspace includes multi-curve OIS/IBOR bootstrapping, Black/Bachelier and arbitrage-aware volatility surfaces, Heston/rough-Heston and Lévy COS engines, ADI PDEs, SIMD Monte Carlo and American MC, scripted payoffs, XVA/SIMM/DIM/MVA, credit/inflation/FX/commodity analytics, barriers/Asians/lookbacks, callable and convertible bonds, and a pure-Rust ANSI/Braille risk workstation.

## Empirical performance

The offline benchmark harness checks numerical parity against independent baselines. Results below are generated from the current release tree; run the benchmark command on your hardware for fresh timings.

**Benchmark summary:** 63 suites and 61 workstation tabs.

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 200119.80 (IQR 419.20) | 3039.60 (IQR 164.40) | 65.84x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 36129.20 (IQR 4083.20) | 2169.60 (IQR 27.40) | 16.65x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 256666.60 (IQR 2815.00) | 5232.40 (IQR 1149.00) | 49.05x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 96126.40 (IQR 953.20) | 30717.20 (IQR 3871.80) | 3.13x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1997621.60 (IQR 2069.20) | 662866.60 (IQR 9468.60) | 3.01x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 437447.40 (IQR 3774.20) | 152327.20 (IQR 732.20) | 2.87x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 840522.60 (IQR 10623.40) | 174837.20 (IQR 4339.60) | 4.81x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 3198957.00 (IQR 11702.80) | 343566.60 (IQR 79415.40) | 9.31x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 384203.80 (IQR 497.40) | 19704.20 (IQR 2916.80) | 19.50x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 61763.20 (IQR 568.00) | 21544.20 (IQR 440.40) | 2.87x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19953039.20 (IQR 326693.80) | 212280.60 (IQR 9412.60) | 93.99x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1801067.00 (IQR 15661.00) | 218483.60 (IQR 4162.40) | 8.24x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37834.00 (IQR 3186.40) | 11761.40 (IQR 469.20) | 3.22x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13763854.80 (IQR 236692.80) | 1201931.40 (IQR 4231.00) | 11.45x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 547327.80 (IQR 2820.40) | 129250.80 (IQR 563.00) | 4.23x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 208928.80 (IQR 4085.80) | 54167.20 (IQR 1015.40) | 3.86x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28937.40 (IQR 304.60) | 6092.60 (IQR 381.80) | 4.75x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1918779.00 (IQR 2873.00) | 11603.00 (IQR 388.40) | 165.37x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 522954.00 (IQR 11336.40) | 56668.00 (IQR 378.80) | 9.23x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7005966.20 (IQR 78998.80) | 155662.00 (IQR 8862.40) | 45.01x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 547470.20 (IQR 7354.80) | 264457.60 (IQR 1297.00) | 2.07x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 510867.40 (IQR 3750.80) | 8109.60 (IQR 1209.80) | 63.00x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 805631.20 (IQR 7167.00) | 153309.00 (IQR 3282.60) | 5.25x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6812.00 (IQR 351.60) | 1598.40 (IQR 208.40) | 4.26x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 315772.00 (IQR 2576.60) | 47157.40 (IQR 2056.00) | 6.70x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 328511.80 (IQR 1482.00) | 71673.00 (IQR 2686.80) | 4.58x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 134280.00 (IQR 304.40) | 63109.60 (IQR 629.20) | 2.13x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3852189.60 (IQR 103320.40) | 449325.00 (IQR 10389.40) | 8.57x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 89740.80 (IQR 1039.20) | 18549.00 (IQR 986.40) | 4.84x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 375702.40 (IQR 14378.00) | 46410.40 (IQR 1809.80) | 8.10x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 1016005.40 (IQR 1081.40) | 146716.20 (IQR 1831.00) | 6.92x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1725815.00 (IQR 31694.00) | 233944.00 (IQR 7048.80) | 7.38x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Monotone-Convex / Clark-FX / SA-CCR-KVA / rVol / Cartea / Paleologo AD | 491828.20 (IQR 2194.20) | 77168.40 (IQR 903.00) | 6.37x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | KRD / LPI / CoCo / Bergomi / CGMY / Guéant / Meucci AD vs FD | 2588773.20 (IQR 11152.80) | 686241.80 (IQR 20441.00) | 3.77x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | RKHS / FMM / VV / LSV / JLT / CTD / QRH / Parisian / Aïd / Propagator / NCO / DML AD | 95247.60 (IQR 460.60) | 26763.60 (IQR 1053.40) | 3.56x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1364416.00 (IQR 2590.00) | 166323.80 (IQR 2541.20) | 8.20x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Kalman-RTS / VECM / HJ-SDF / IFRS9 / CPT / Hasbrouck / FV / Lando / LPPL / CB-Arb / Causal / Malliavin AD | 141743.40 (IQR 893.60) | 15978.60 (IQR 712.60) | 8.87x | 8.967e-6 |
| 38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc | MF-Rates / Wildcard / CTD / Pivot-TARF / Bergomi / Schönbucher / FRTB-CVA / Lookback / VTA / FLAM / ACD / HLZ / MLMC AD | 1123132.80 (IQR 3472.00) | 49695.20 (IQR 800.20) | 22.60x | 4.996e-9 |
| 39_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt | FTPL / G2++ / eSSVI / Offshore-CB / CoCo / MPOR / Double-Heston / Curran / Warin / Roll-Corwin / Busseti-Kelly / Diff-QP AD | 673877.80 (IQR 10483.80) | 38573.80 (IQR 974.40) | 17.47x | 4.423e-10 |
| 40_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca | BCBS-IRRBB / Turfus-RFR / Baran-MtM / Bourgey-BG / FEDS-PONV / Zeron-SIMM / Jacquier-rBergomi / Choi-Spread / BESS-SDP / Kratz-Dark / Föllmer / Diff-PCA AD | 273469.20 (IQR 1821.40) | 25199.80 (IQR 753.80) | 10.85x | 1.803e-9 |
| 41_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is | Veronesi / Lando / Gatheral-Corridor / Russian / Bates-SVJ / Koziol-CB / MFG / BCBS-Floor / CJP-LOB / ETF / Qian-IC / Glasserman-IS AD | 706259.20 (IQR 6884.60) | 90637.80 (IQR 1245.00) | 7.79x | 2.244e-9 |
| 42_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin | Romer-Taylor / Damodaran-CRP / Cumulant-VaR / LOB-LVaR / Prospect / Regret / Beneish / Campisi / Heckman / Elastic-Net / Reset-CB / FSB-BailIn AD | 22082.80 (IQR 238.00) | 2443.60 (IQR 190.40) | 9.04x | 8.962e-10 |
| 43_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi | KRD-LSC / CMO-PAC / Asquith-CB / ASU2024-CPDI / ASC470-Mod / ColVA-CTD / Hill-EVT / Bouchaud-TLF / Karnosky-FX / DGTW-Timing / Brock-Hommes / Penman-ReOI AD | 63192.20 (IQR 555.80) | 5694.60 (IQR 381.80) | 11.10x | 1.010e-9 |
| 44_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan | Cheyette-LV / ADB-NDIRS / Andreasen-ZABR / QRH-Zumbach / CMCDS-Morini / Turlakov-WWR / HHW-ADI / CappedCall-MW / Ahmadi-EVaR / Schur-HRP / BEKK-GARCH / Schilit-Sloan AD | 43151.20 (IQR 414.60) | 4536.20 (IQR 106.20) | 9.51x | 5.153e-8 |
| 45_rebonato_cms_leland_toft_calamos_ascot_spec_pool_roll_glosten_milgrom_potters_rie_egarch_hamilton_wooldridge_2sls_mcneil_euler_bacon_omega_damodaran_rd_mwg_arrow_debreu | Rebonato-CMS / Leland-Toft / Calamos-Ascot / SpecPool-Roll / Glosten-Milgrom / Potters-RIE / EGARCH-Hamilton / Wooldridge-2SLS / McNeil-Euler / Bacon-Omega / Damodaran-RD / Arrow-Debreu AD | 14973.60 (IQR 44.00) | 1741.40 (IQR 103.80) | 8.60x | 6.115e-9 |
| 46_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig | Tuckman-Repo / BK-FVA / DECS / MBS-IO-PO / Kyle-PIN / GP-Aim / GJR-DCC / LPPL / Spectral-Risk / Brinson-Carino / Prelec-CPT / Path-Sig AD | 35170.40 (IQR 315.40) | 2751.20 (IQR 76.20) | 12.78x | 5.156e-9 |
| 47_hagan_west_mc_fabozzi_ctd_switch_afv_default_cb_lando_kth_default_gregory_sa_cva_gibson_schwartz_kirk_clark_fx_vv_bouzoubaa_napoleon_grinold_tc_flam_meucci_entropy_hastie_krr_loocv_koller_mckinsey_eva | Hagan-West-MC / CTD-Switch / AFV-CB / Kth-Default / SA-CVA-Euler / GS-Kirk-Spread / FX-Vanna-Volga / Napoleon-Cliquet / Transfer-Coef / Meucci-Entropy / KRR-LOOCV / McKinsey-EVA AD | 38890.80 (IQR 62.20) | 4236.80 (IQR 130.00) | 9.18x | 2.045e-8 |
| 48_brigo_csa_ctd_fabozzi_callable_oas_lando_jlt_migration_gatheral_corridor_varswap_glasserman_lsm_rogers_dual_cartea_lob_mm_lo_mackinlay_vr_lopez_dsr_purged_cv_roncalli_spinu_erc_ilmanen_vmc_penman_reoi_mccall_bellman | Multi-CSA-CTD / Callable-OAS / JLT-Migration / Corridor-VarSwap / LSM-Rogers-Dual / Cartea-LOB-MM / Lo-MacKinlay-VR / Lopez-DSR-PurgedCV / Roncalli-ERC / Ilmanen-VMC / Penman-ReOI / McCall-Bellman AD | 75665.20 (IQR 554.60) | 8010.00 (IQR 308.80) | 9.45x | 4.941e-8 |
| 49_brigo_nss_krd_tuckman_tips_floor_schonbucher_rmv_rt_rfv_gatheral_sticky_lv_gregory_simm_ccp_mva_schwartz_smith_kalman_wystup_fx_dnt_bouzoubaa_timer_meucci_bl_woodbury_gueant_ac_sinh_zou_elastic_net_nk_ree | NSS-KRD / TIPS-Floor-BEI / RMV-RT-RFV / Sticky-LocalVol / SIMM-CCP-MVA / Schwartz-Smith-KF / FX-DNT-VV / Timer-VarClock / BL-Woodbury / AC-Sinh-Liq / ElasticNet-Group / NK-DSGE-REE AD | 40909.80 (IQR 331.80) | 3317.60 (IQR 131.40) | 12.33x | 1.279e-9 |
| 50_hagan_west_tension_duffie_lando_pac_cmo_fama_macbeth_shanken_gatheral_asvi_bgk_discrete_barrier_bouchaud_cross_impact_karnosky_singer_merton_jump_campbell_viceira_ohlson_lid_griewank_revolve | Hagan-West-Tension / Duffie-Lando / PAC-CMO / Fama-MacBeth-Shanken / Gatheral-aSVI / BGK-Barrier / Cross-Impact / Karnosky-Singer / Merton-Jump / Campbell-Viceira / Ohlson-LID / Revolve-BSB AD | 40292.00 (IQR 57.60) | 7298.00 (IQR 105.80) | 5.52x | 1.448e-8 |
| 51_brigo_xccy_mtm_fabozzi_spec_pool_io_po_lando_marshall_olkin_gatheral_vix_vvix_gregory_wwr_cva_geman_gibson_schwartz_clark_fx_fva_rubinstein_chooser_grinold_flam_tc_meucci_es_euler_james_stein_sure_campbell_cochrane_habit | XCCY-MtM-Basis / Spec-Pool-IO-PO / Marshall-Olkin / VIX-VVIX-Smile / WWR-CVA / Gibson-Schwartz / FX-FVA / Chooser / FLAM-TC / ES-Euler / James-Stein-SURE / Campbell-Cochrane AD | 8038.80 (IQR 14.60) | 1190.40 (IQR 29.20) | 6.75x | 2.010e-8 |
| 52_brigo_ois_fallback_fabozzi_collared_frn_schonbucher_cds_swaption_fep_bergomi_2f_skew_gregory_sa_ccr_ead_eydeland_spark_spread_wystup_prdc_geske_compound_qian_alpha_decay_roncalli_enb_johansen_vecm_bansal_yaron_lrr | OIS-Fallback / Collared-FRN / CDS-Swaption-FEP / Bergomi-2F-Skew / SA-CCR-EAD / Clean-Spark / FX-PRDC / Geske-Compound / Alpha-Decay-IR / Meucci-ENB / Johansen-VECM / Bansal-Yaron-LRR AD | 30535.60 (IQR 23.60) | 3270.60 (IQR 43.60) | 9.34x | 1.469e-9 |
| 53_andersen_cms_replication_tuckman_special_repo_okane_cds_upfront_gatheral_ssr_tanaka_gregory_saccr_kva_eydeland_crack_spread_wystup_prdc_cap_haug_eso_dilution_qian_spinu_erc_rockafellar_cvar_htf_pls_krylov_piazzesi_fomc | CMS-Replication / Special-Repo-CTD / ISDA-CDS-Upfront / SSR-Tanaka / SA-CCR-KVA / 3-Asset-Crack / PRDC-Cap / ESO-Dilution / Spinu-ERC / Rockafellar-CVaR / PLS-Krylov / FOMC-Hike AD | 20400.60 (IQR 38.00) | 1911.80 (IQR 41.40) | 10.67x | 2.693e-9 |
| 54_darbyshire_step_ois_turn_despiegeleer_cb_soft_call_okane_cds_index_skew_gatheral_varswap_truncation_gregory_fva_overlap_geman_working_storage_wystup_malz_fx_fly_goldman_sosin_lookback_ledoit_wolf_shrinkage_almgren_chriss_cara_tsay_estar_acm_cp | Step-OIS-Turn / CB-SoftCall-MW / CDX-Skew-Basis / VarSwap-Trunc / FVA-DVA-Overlap / Working-Storage / FX-Malz-Fly / GSG-Lookback / Ledoit-Wolf-CC / Almgren-CARA / ESTAR-Band / ACM-CP-TP AD | 13314.20 (IQR 381.20) | 1250.00 (IQR 45.80) | 10.65x | 1.054e-8 |
| 55_henrard_tenor_basis_fabozzi_msr_recapture_okane_cdo_base_corr_gatheral_fwd_start_smile_gregory_pykhtin_csa_mta_eydeland_gas_swing_wystup_fx_kiko_stulz_rainbow_jagannathan_ma_obizhaeva_wang_newey_west_hac_gabaix_disaster | Tenor-Basis-Mult / MSR-Recapture / CDO-Base-Corr / Fwd-Start-Smile / CSA-MTA-CVA / Gas-Swing / FX-Window-KIKO / Stulz-Rainbow / Jagannathan-Ma / Obizhaeva-Wang / Newey-West-HAC / Gabaix-Disaster AD | 108286.40 (IQR 497.20) | 9594.20 (IQR 326.00) | 11.29x | 3.026e-7 |
| 56_brigo_diff_swap_quanto_fabozzi_callable_muni_ytw_bielecki_step_up_bond_bkm_moment_strip_gregory_break_clause_cva_bessembinder_power_skew_wystup_fx_fader_longstaff_extendible_kan_zhou_3fund_sharpe_rbsa_diebold_mariano_hln_laubach_williams_rstar | Diff-Swap-Quanto / Callable-Muni-YTW / Step-Up-Bond / BKM-SKEW / Break-Clause-CVA / Power-Skew / FX-Fader / Longstaff-Extendible / Kan-Zhou-3Fund / Sharpe-RBSA / Diebold-Mariano-HLN / Laubach-Williams-r* AD | 62117.80 (IQR 379.00) | 7489.00 (IQR 51.80) | 8.29x | 8.489e-9 |
| 57_andersen_sofr_futures_fra_hayre_mbs_scurve_pedersen_cdx_fep_gatheral_dupire_total_var_gregory_simm_mva_routledge_stockout_baba_cip_basis_carr_chou_symmetry_pastor_stambaugh_treynor_black_engle_granger_ecm_epstein_zin_sdf | SOFR-Fut-FRA / Hayre-MBS-SCurve / CDX-Option-FEP / Dupire-TotalVar-LV / SIMM-MVA-Funding / RSS-Stockout / CIP-Basis-Wedge / Carr-Chou-Symmetry / Pastor-Stambaugh / Treynor-Black / Engle-Granger-ECM / Epstein-Zin-SDF AD | 6802.20 (IQR 7.80) | 1017.80 (IQR 40.80) | 6.68x | 8.576e-7 |
| 58_litterman_scheinkman_pca_fly_tf_split_cb_hull_white_double_t_cdo_roger_lee_moment_wings_brigo_capponi_cds_contagion_schwartz_ou_futures_wystup_dcd_yield_passport_option_grs_factor_ftest_bouchaud_sqrt_impact_bns_bipower_har_cj_lettau_ludvigson_cay | PCA-Fly-Carry / TF-Split-CB / Double-t-CDO / Roger-Lee-Wings / CDS-Contagion-CVA / Schwartz-OU-Fut / FX-DCD-Yield / Passport-Opt / GRS-Factor-FTest / Sqrt-Metaorder / BNS-HAR-CJ / Lettau-Ludvigson-cay AD | 10658.20 (IQR 230.60) | 1037.00 (IQR 36.80) | 10.28x | 3.312e-9 |
| 59_piterbarg_collateral_option_tuckman_rolldown_merton_kmv_dd_hagan_bartlett_sabr_pykhtin_mpor_spike_deng_locational_ftr_wystup_digital_overhedge_conze_fixed_lookback_tasche_hot_spots_es_frazzini_bab_pesaran_timmermann_lucas_mehra_prescott | Piterbarg-Coll-Opt / Tuckman-Rolldown / Merton-KMV-DD / Bartlett-SABR / Pykhtin-MPoR-Spike / Deng-FTR-Opt / Wystup-FX-Digital / Conze-Fixed-Lookback / Tasche-Hot-Spots-ES / Frazzini-BAB / Pesaran-Timmermann / Lucas-Mehra-Prescott AD | 10629.80 (IQR 235.20) | 1589.20 (IQR 34.60) | 6.69x | 6.269e-8 |
| 60_gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_ga_carr_wu_vol_drift_gregory_marginal_cva_working_brennan_storage_wystup_forward_extra_rubinstein_exchange_opt_hlz_factor_fdr_holden_subrahmanyam_kyle_bollerslev_garch_kurtosis_campbell_shiller_cape | GSW-Svensson-Curve / Fabozzi-MBS-ESF / Vasicek-Gordy-IRB-GA / Carr-Wu-Vol-Drift / Gregory-Marginal-CVA / Working-Brennan-Storage / Wystup-Forward-Extra / Rubinstein-Exchange-Opt / HLZ-Factor-FDR / Holden-Subrahmanyam-Kyle / Bollerslev-GARCH-Kurtosis / Campbell-Shiller-CAPE AD | 35988.60 (IQR 363.40) | 2658.60 (IQR 60.60) | 13.54x | 1.252e-8 |
| 61_kim_wright_shadow_rate_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_swap_burgard_kjaer_fva_carmona_durrleman_spark_wystup_fx_tarf_bgk_discrete_lookback_novy_marx_ff5_acerbi_spectral_es_corsi_har_rv_bll_news_noise | Kim-Wright-Shadow / Brennan-Schwartz-CB / Duffie-Singleton-RMV / Fukasawa-Gamma-Swap / Burgard-Kjaer-FVA / Carmona-Durrleman-Spark / Wystup-FX-TARF / BGK-Discrete-Lookback / Novy-Marx-FF5 / Acerbi-Spectral-ES / Corsi-HAR-RV / BLL-News-Noise AD | 27117.60 (IQR 34.60) | 1891.00 (IQR 52.20) | 14.34x | 2.428e-9 |
| 62_jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_regimes_hull_white_wwr_cva_gibson_schwartz_samuelson_reiner_quanto_fx_turnbull_wakeman_asian_amp_value_momentum_garleanu_pedersen_aim_bai_ng_icp2_bgg_accelerator | JPS-Unspanned-Macro / Kalotay-Refunding / Hull-White-NtD / Derman-Sticky-Regimes / Hull-White-WWR-CVA / Gibson-Schwartz-Samuelson / Reiner-Quanto-FX / Turnbull-Wakeman-Asian / AMP-Value-Momentum / Garleanu-Pedersen-Aim / Bai-Ng-ICp2 / BGG-Accelerator AD | 36423.40 (IQR 458.40) | 2666.60 (IQR 88.80) | 13.66x | 4.361e-7 |
| 63_lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_cds_fouque_multiscale_sv_gregory_bilateral_ftd_eydeland_vpp_heat_rate_wystup_fx_corridor_rubinstein_reiner_barrier_lns_gls_r2_bellini_expectile_bns_quarticity_kiyotaki_moore_cycle | Lyashenko-FMM-Caplet / Tuckman-CPI-Seasonality / Brigo-Alfonsi-CIR-CDS / Fouque-Multiscale-SV / Gregory-Bilateral-FtD / Eydeland-VPP-Heat-Rate / Wystup-FX-Corridor / Rubinstein-Reiner-Barrier / LNS-GLS-R2 / Bellini-Expectile / BNS-Quarticity / Kiyotaki-Moore-Cycle AD | 21819.60 (IQR 13.20) | 1789.20 (IQR 106.40) | 12.20x | 3.248e-9 |

Summary: 63 benchmark suites; speedups range from **2.07x** to **165.37x** (median **8.60x**, mean **15.26x**). All parity checks pass their domain tolerances; latest suite parity is **3.248e-9**.

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

