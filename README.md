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
| `ql-benches` | `crates/ql-benches` | Delivery | 60-suite parity/timing harness |
| `ql-terminal` | `crates/ql-terminal` | Delivery | Unified SDK, 58-tab Braille/ANSI workstation |

## Numerical capabilities

The workspace includes multi-curve OIS/IBOR bootstrapping, Black/Bachelier and arbitrage-aware volatility surfaces, Heston/rough-Heston and Lévy COS engines, ADI PDEs, SIMD Monte Carlo and American MC, scripted payoffs, XVA/SIMM/DIM/MVA, credit/inflation/FX/commodity analytics, barriers/Asians/lookbacks, callable and convertible bonds, and a pure-Rust ANSI/Braille risk workstation.

## Empirical performance

The offline benchmark harness checks numerical parity against independent baselines. Results below are generated from the current release tree; run the benchmark command on your hardware for fresh timings.

**Benchmark summary:** 60 suites and 58 workstation tabs.

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 357351.80 (IQR 192820.40) | 3972.40 (IQR 164.40) | 89.96x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 41249.60 (IQR 4331.60) | 2422.40 (IQR 7.60) | 17.03x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 247062.00 (IQR 6055.60) | 4581.00 (IQR 246.20) | 53.93x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 94230.40 (IQR 1926.60) | 30435.00 (IQR 367.20) | 3.10x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1949435.60 (IQR 41292.20) | 613741.60 (IQR 702.40) | 3.18x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 429696.60 (IQR 7660.20) | 149096.40 (IQR 296.00) | 2.88x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 797230.20 (IQR 9098.00) | 188630.40 (IQR 2708.60) | 4.23x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 3015653.20 (IQR 134939.80) | 366481.80 (IQR 86157.80) | 8.23x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 384207.00 (IQR 809.40) | 17163.80 (IQR 2106.20) | 22.38x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 61913.20 (IQR 270.00) | 21587.40 (IQR 49.40) | 2.87x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19569017.20 (IQR 111601.40) | 203997.80 (IQR 11174.60) | 95.93x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1769892.80 (IQR 10634.20) | 213335.80 (IQR 2975.00) | 8.30x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37045.60 (IQR 2814.00) | 12105.40 (IQR 1847.80) | 3.06x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13628533.80 (IQR 291808.80) | 1197017.00 (IQR 3521.00) | 11.39x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 567268.00 (IQR 539.40) | 128902.40 (IQR 1296.40) | 4.40x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 206871.00 (IQR 518.40) | 54937.80 (IQR 1164.60) | 3.77x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 29169.20 (IQR 349.20) | 5912.00 (IQR 212.80) | 4.93x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1930692.40 (IQR 6713.80) | 11528.60 (IQR 287.00) | 167.47x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 520397.00 (IQR 1216.60) | 55358.20 (IQR 340.60) | 9.40x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 6953813.20 (IQR 40939.20) | 151582.20 (IQR 3968.00) | 45.87x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 553531.20 (IQR 922.60) | 263131.20 (IQR 740.60) | 2.10x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 506562.80 (IQR 2356.00) | 8150.20 (IQR 1257.00) | 62.15x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 810113.40 (IQR 24622.20) | 142373.60 (IQR 2607.20) | 5.69x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6767.60 (IQR 187.40) | 1343.80 (IQR 169.40) | 5.04x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 303536.20 (IQR 721.20) | 45883.60 (IQR 1654.20) | 6.62x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 326208.20 (IQR 415.00) | 71468.80 (IQR 3151.40) | 4.56x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 133656.20 (IQR 1753.00) | 62036.60 (IQR 923.00) | 2.15x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3812511.00 (IQR 2034.20) | 437629.40 (IQR 407.40) | 8.71x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 87431.00 (IQR 587.00) | 17405.20 (IQR 952.80) | 5.02x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 360001.60 (IQR 11936.60) | 43889.60 (IQR 425.80) | 8.20x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 970886.80 (IQR 1356.40) | 138601.80 (IQR 1876.60) | 7.00x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1689856.60 (IQR 16722.80) | 231483.20 (IQR 2582.80) | 7.30x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Monotone-Convex / Clark-FX / SA-CCR-KVA / rVol / Cartea / Paleologo AD | 491562.20 (IQR 483.60) | 77179.80 (IQR 1175.40) | 6.37x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | KRD / LPI / CoCo / Bergomi / CGMY / Guéant / Meucci AD vs FD | 2585753.40 (IQR 2935.00) | 680670.00 (IQR 13344.20) | 3.80x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | RKHS / FMM / VV / LSV / JLT / CTD / QRH / Parisian / Aïd / Propagator / NCO / DML AD | 94208.40 (IQR 420.00) | 26671.40 (IQR 1337.20) | 3.53x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1362379.40 (IQR 1730.80) | 165274.20 (IQR 4195.20) | 8.24x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Kalman-RTS / VECM / HJ-SDF / IFRS9 / CPT / Hasbrouck / FV / Lando / LPPL / CB-Arb / Causal / Malliavin AD | 138838.40 (IQR 364.20) | 16424.40 (IQR 455.60) | 8.45x | 8.967e-6 |
| 38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc | MF-Rates / Wildcard / CTD / Pivot-TARF / Bergomi / Schönbucher / FRTB-CVA / Lookback / VTA / FLAM / ACD / HLZ / MLMC AD | 1124613.20 (IQR 137.60) | 48947.20 (IQR 1336.00) | 22.98x | 4.996e-9 |
| 39_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt | FTPL / G2++ / eSSVI / Offshore-CB / CoCo / MPOR / Double-Heston / Curran / Warin / Roll-Corwin / Busseti-Kelly / Diff-QP AD | 668970.80 (IQR 2409.40) | 39106.40 (IQR 731.20) | 17.11x | 4.423e-10 |
| 40_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca | BCBS-IRRBB / Turfus-RFR / Baran-MtM / Bourgey-BG / FEDS-PONV / Zeron-SIMM / Jacquier-rBergomi / Choi-Spread / BESS-SDP / Kratz-Dark / Föllmer / Diff-PCA AD | 278234.60 (IQR 2339.00) | 26017.00 (IQR 443.60) | 10.69x | 1.803e-9 |
| 41_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is | Veronesi / Lando / Gatheral-Corridor / Russian / Bates-SVJ / Koziol-CB / MFG / BCBS-Floor / CJP-LOB / ETF / Qian-IC / Glasserman-IS AD | 714907.60 (IQR 3260.20) | 90983.00 (IQR 1176.40) | 7.86x | 2.244e-9 |
| 42_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin | Romer-Taylor / Damodaran-CRP / Cumulant-VaR / LOB-LVaR / Prospect / Regret / Beneish / Campisi / Heckman / Elastic-Net / Reset-CB / FSB-BailIn AD | 22097.20 (IQR 52.40) | 2460.00 (IQR 207.60) | 8.98x | 8.962e-10 |
| 43_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi | KRD-LSC / CMO-PAC / Asquith-CB / ASU2024-CPDI / ASC470-Mod / ColVA-CTD / Hill-EVT / Bouchaud-TLF / Karnosky-FX / DGTW-Timing / Brock-Hommes / Penman-ReOI AD | 63316.80 (IQR 452.00) | 5693.80 (IQR 388.20) | 11.12x | 1.010e-9 |
| 44_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan | Cheyette-LV / ADB-NDIRS / Andreasen-ZABR / QRH-Zumbach / CMCDS-Morini / Turlakov-WWR / HHW-ADI / CappedCall-MW / Ahmadi-EVaR / Schur-HRP / BEKK-GARCH / Schilit-Sloan AD | 43266.00 (IQR 347.40) | 4510.40 (IQR 93.40) | 9.59x | 5.153e-8 |
| 45_rebonato_cms_leland_toft_calamos_ascot_spec_pool_roll_glosten_milgrom_potters_rie_egarch_hamilton_wooldridge_2sls_mcneil_euler_bacon_omega_damodaran_rd_mwg_arrow_debreu | Rebonato-CMS / Leland-Toft / Calamos-Ascot / SpecPool-Roll / Glosten-Milgrom / Potters-RIE / EGARCH-Hamilton / Wooldridge-2SLS / McNeil-Euler / Bacon-Omega / Damodaran-RD / Arrow-Debreu AD | 15266.40 (IQR 57.80) | 1724.00 (IQR 86.40) | 8.86x | 6.115e-9 |
| 46_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig | Tuckman-Repo / BK-FVA / DECS / MBS-IO-PO / Kyle-PIN / GP-Aim / GJR-DCC / LPPL / Spectral-Risk / Brinson-Carino / Prelec-CPT / Path-Sig AD | 35255.60 (IQR 340.60) | 2741.60 (IQR 116.00) | 12.86x | 5.156e-9 |
| 47_hagan_west_mc_fabozzi_ctd_switch_afv_default_cb_lando_kth_default_gregory_sa_cva_gibson_schwartz_kirk_clark_fx_vv_bouzoubaa_napoleon_grinold_tc_flam_meucci_entropy_hastie_krr_loocv_koller_mckinsey_eva | Hagan-West-MC / CTD-Switch / AFV-CB / Kth-Default / SA-CVA-Euler / GS-Kirk-Spread / FX-Vanna-Volga / Napoleon-Cliquet / Transfer-Coef / Meucci-Entropy / KRR-LOOCV / McKinsey-EVA AD | 39214.20 (IQR 466.80) | 4141.80 (IQR 128.80) | 9.47x | 2.045e-8 |
| 48_brigo_csa_ctd_fabozzi_callable_oas_lando_jlt_migration_gatheral_corridor_varswap_glasserman_lsm_rogers_dual_cartea_lob_mm_lo_mackinlay_vr_lopez_dsr_purged_cv_roncalli_spinu_erc_ilmanen_vmc_penman_reoi_mccall_bellman | Multi-CSA-CTD / Callable-OAS / JLT-Migration / Corridor-VarSwap / LSM-Rogers-Dual / Cartea-LOB-MM / Lo-MacKinlay-VR / Lopez-DSR-PurgedCV / Roncalli-ERC / Ilmanen-VMC / Penman-ReOI / McCall-Bellman AD | 73860.80 (IQR 445.00) | 7760.20 (IQR 288.00) | 9.52x | 4.941e-8 |
| 49_brigo_nss_krd_tuckman_tips_floor_schonbucher_rmv_rt_rfv_gatheral_sticky_lv_gregory_simm_ccp_mva_schwartz_smith_kalman_wystup_fx_dnt_bouzoubaa_timer_meucci_bl_woodbury_gueant_ac_sinh_zou_elastic_net_nk_ree | NSS-KRD / TIPS-Floor-BEI / RMV-RT-RFV / Sticky-LocalVol / SIMM-CCP-MVA / Schwartz-Smith-KF / FX-DNT-VV / Timer-VarClock / BL-Woodbury / AC-Sinh-Liq / ElasticNet-Group / NK-DSGE-REE AD | 40119.60 (IQR 265.00) | 3217.60 (IQR 129.60) | 12.47x | 1.279e-9 |
| 50_hagan_west_tension_duffie_lando_pac_cmo_fama_macbeth_shanken_gatheral_asvi_bgk_discrete_barrier_bouchaud_cross_impact_karnosky_singer_merton_jump_campbell_viceira_ohlson_lid_griewank_revolve | Hagan-West-Tension / Duffie-Lando / PAC-CMO / Fama-MacBeth-Shanken / Gatheral-aSVI / BGK-Barrier / Cross-Impact / Karnosky-Singer / Merton-Jump / Campbell-Viceira / Ohlson-LID / Revolve-BSB AD | 39394.60 (IQR 37.20) | 7115.40 (IQR 121.00) | 5.54x | 1.448e-8 |
| 51_brigo_xccy_mtm_fabozzi_spec_pool_io_po_lando_marshall_olkin_gatheral_vix_vvix_gregory_wwr_cva_geman_gibson_schwartz_clark_fx_fva_rubinstein_chooser_grinold_flam_tc_meucci_es_euler_james_stein_sure_campbell_cochrane_habit | XCCY-MtM-Basis / Spec-Pool-IO-PO / Marshall-Olkin / VIX-VVIX-Smile / WWR-CVA / Gibson-Schwartz / FX-FVA / Chooser / FLAM-TC / ES-Euler / James-Stein-SURE / Campbell-Cochrane AD | 7895.20 (IQR 221.80) | 1160.40 (IQR 25.60) | 6.80x | 2.010e-8 |
| 52_brigo_ois_fallback_fabozzi_collared_frn_schonbucher_cds_swaption_fep_bergomi_2f_skew_gregory_sa_ccr_ead_eydeland_spark_spread_wystup_prdc_geske_compound_qian_alpha_decay_roncalli_enb_johansen_vecm_bansal_yaron_lrr | OIS-Fallback / Collared-FRN / CDS-Swaption-FEP / Bergomi-2F-Skew / SA-CCR-EAD / Clean-Spark / FX-PRDC / Geske-Compound / Alpha-Decay-IR / Meucci-ENB / Johansen-VECM / Bansal-Yaron-LRR AD | 29934.60 (IQR 458.80) | 3171.20 (IQR 56.20) | 9.44x | 1.469e-9 |
| 53_andersen_cms_replication_tuckman_special_repo_okane_cds_upfront_gatheral_ssr_tanaka_gregory_saccr_kva_eydeland_crack_spread_wystup_prdc_cap_haug_eso_dilution_qian_spinu_erc_rockafellar_cvar_htf_pls_krylov_piazzesi_fomc | CMS-Replication / Special-Repo-CTD / ISDA-CDS-Upfront / SSR-Tanaka / SA-CCR-KVA / 3-Asset-Crack / PRDC-Cap / ESO-Dilution / Spinu-ERC / Rockafellar-CVaR / PLS-Krylov / FOMC-Hike AD | 19962.80 (IQR 19.00) | 1861.60 (IQR 38.40) | 10.72x | 2.693e-9 |
| 54_darbyshire_step_ois_turn_despiegeleer_cb_soft_call_okane_cds_index_skew_gatheral_varswap_truncation_gregory_fva_overlap_geman_working_storage_wystup_malz_fx_fly_goldman_sosin_lookback_ledoit_wolf_shrinkage_almgren_chriss_cara_tsay_estar_acm_cp | Step-OIS-Turn / CB-SoftCall-MW / CDX-Skew-Basis / VarSwap-Trunc / FVA-DVA-Overlap / Working-Storage / FX-Malz-Fly / GSG-Lookback / Ledoit-Wolf-CC / Almgren-CARA / ESTAR-Band / ACM-CP-TP AD | 13036.60 (IQR 352.00) | 1216.60 (IQR 26.80) | 10.72x | 1.054e-8 |
| 55_henrard_tenor_basis_fabozzi_msr_recapture_okane_cdo_base_corr_gatheral_fwd_start_smile_gregory_pykhtin_csa_mta_eydeland_gas_swing_wystup_fx_kiko_stulz_rainbow_jagannathan_ma_obizhaeva_wang_newey_west_hac_gabaix_disaster | Tenor-Basis-Mult / MSR-Recapture / CDO-Base-Corr / Fwd-Start-Smile / CSA-MTA-CVA / Gas-Swing / FX-Window-KIKO / Stulz-Rainbow / Jagannathan-Ma / Obizhaeva-Wang / Newey-West-HAC / Gabaix-Disaster AD | 106149.80 (IQR 401.60) | 9381.60 (IQR 284.00) | 11.31x | 3.026e-7 |
| 56_brigo_diff_swap_quanto_fabozzi_callable_muni_ytw_bielecki_step_up_bond_bkm_moment_strip_gregory_break_clause_cva_bessembinder_power_skew_wystup_fx_fader_longstaff_extendible_kan_zhou_3fund_sharpe_rbsa_diebold_mariano_hln_laubach_williams_rstar | Diff-Swap-Quanto / Callable-Muni-YTW / Step-Up-Bond / BKM-SKEW / Break-Clause-CVA / Power-Skew / FX-Fader / Longstaff-Extendible / Kan-Zhou-3Fund / Sharpe-RBSA / Diebold-Mariano-HLN / Laubach-Williams-r* AD | 60837.20 (IQR 294.60) | 7317.40 (IQR 56.80) | 8.31x | 8.489e-9 |
| 57_andersen_sofr_futures_fra_hayre_mbs_scurve_pedersen_cdx_fep_gatheral_dupire_total_var_gregory_simm_mva_routledge_stockout_baba_cip_basis_carr_chou_symmetry_pastor_stambaugh_treynor_black_engle_granger_ecm_epstein_zin_sdf | SOFR-Fut-FRA / Hayre-MBS-SCurve / CDX-Option-FEP / Dupire-TotalVar-LV / SIMM-MVA-Funding / RSS-Stockout / CIP-Basis-Wedge / Carr-Chou-Symmetry / Pastor-Stambaugh / Treynor-Black / Engle-Granger-ECM / Epstein-Zin-SDF AD | 6681.80 (IQR 6.80) | 999.80 (IQR 15.00) | 6.68x | 8.576e-7 |
| 58_litterman_scheinkman_pca_fly_tf_split_cb_hull_white_double_t_cdo_roger_lee_moment_wings_brigo_capponi_cds_contagion_schwartz_ou_futures_wystup_dcd_yield_passport_option_grs_factor_ftest_bouchaud_sqrt_impact_bns_bipower_har_cj_lettau_ludvigson_cay | PCA-Fly-Carry / TF-Split-CB / Double-t-CDO / Roger-Lee-Wings / CDS-Contagion-CVA / Schwartz-OU-Fut / FX-DCD-Yield / Passport-Opt / GRS-Factor-FTest / Sqrt-Metaorder / BNS-HAR-CJ / Lettau-Ludvigson-cay AD | 10079.40 (IQR 21.60) | 1002.60 (IQR 32.20) | 10.05x | 3.312e-9 |
| 59_piterbarg_collateral_option_tuckman_rolldown_merton_kmv_dd_hagan_bartlett_sabr_pykhtin_mpor_spike_deng_locational_ftr_wystup_digital_overhedge_conze_fixed_lookback_tasche_hot_spots_es_frazzini_bab_pesaran_timmermann_lucas_mehra_prescott | Piterbarg-Coll-Opt / Tuckman-Rolldown / Merton-KMV-DD / Bartlett-SABR / Pykhtin-MPoR-Spike / Deng-FTR-Opt / Wystup-FX-Digital / Conze-Fixed-Lookback / Tasche-Hot-Spots-ES / Frazzini-BAB / Pesaran-Timmermann / Lucas-Mehra-Prescott AD | 10395.40 (IQR 22.40) | 1592.80 (IQR 52.80) | 6.53x | 6.269e-8 |
| 60_gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_ga_carr_wu_vol_drift_gregory_marginal_cva_working_brennan_storage_wystup_forward_extra_rubinstein_exchange_opt_hlz_factor_fdr_holden_subrahmanyam_kyle_bollerslev_garch_kurtosis_campbell_shiller_cape | GSW-Svensson-Curve / Fabozzi-MBS-ESF / Vasicek-Gordy-IRB-GA / Carr-Wu-Vol-Drift / Gregory-Marginal-CVA / Working-Brennan-Storage / Wystup-Forward-Extra / Rubinstein-Exchange-Opt / HLZ-Factor-FDR / Holden-Subrahmanyam-Kyle / Bollerslev-GARCH-Kurtosis / Campbell-Shiller-CAPE AD | 35914.40 (IQR 305.00) | 2657.60 (IQR 55.40) | 13.51x | 1.252e-8 |

Summary: 60 benchmark suites; speedups range from **2.10x** to **167.47x** (median **8.31x**, mean **15.97x**). All parity checks pass their domain tolerances; latest suite parity is **1.252e-8**.

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

