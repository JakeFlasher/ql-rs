# Competitive Benchmark Synthesis

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 623404.60 (IQR 456824.00) | 5960.20 (IQR 265.80) | 104.59x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 63234.00 (IQR 6340.60) | 3769.00 (IQR 26.40) | 16.78x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 270071.40 (IQR 90407.80) | 4601.00 (IQR 229.20) | 58.70x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 95525.40 (IQR 767.60) | 30835.80 (IQR 702.60) | 3.10x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1925267.40 (IQR 32948.80) | 738728.80 (IQR 7178.20) | 2.61x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 429247.80 (IQR 8296.00) | 150665.40 (IQR 1375.80) | 2.85x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 792945.00 (IQR 8437.40) | 173048.20 (IQR 3055.80) | 4.58x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 2975972.00 (IQR 45642.40) | 262540.00 (IQR 95572.00) | 11.34x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 385407.80 (IQR 1151.60) | 17037.20 (IQR 2186.00) | 22.62x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 61935.60 (IQR 298.80) | 21630.80 (IQR 28.80) | 2.86x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19386521.80 (IQR 141681.20) | 213319.40 (IQR 15064.40) | 90.88x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1772647.40 (IQR 13796.20) | 216145.40 (IQR 2144.40) | 8.20x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37047.80 (IQR 2609.80) | 13225.40 (IQR 775.80) | 2.80x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13636450.40 (IQR 34782.40) | 1201616.00 (IQR 4959.60) | 11.35x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 544896.40 (IQR 502.80) | 129231.00 (IQR 1254.00) | 4.22x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 207009.20 (IQR 268.00) | 53475.00 (IQR 1641.40) | 3.87x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28207.00 (IQR 311.60) | 5810.20 (IQR 219.80) | 4.85x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1893872.60 (IQR 5812.00) | 11059.20 (IQR 212.40) | 171.25x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 511033.40 (IQR 812.20) | 54466.20 (IQR 2831.60) | 9.38x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7287315.40 (IQR 6979.40) | 152363.00 (IQR 4813.60) | 47.83x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 539725.80 (IQR 4061.40) | 259807.00 (IQR 3446.60) | 2.08x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 514605.80 (IQR 9118.80) | 7973.20 (IQR 926.60) | 64.54x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 794083.40 (IQR 5409.00) | 131210.00 (IQR 2961.80) | 6.05x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6473.80 (IQR 347.00) | 1431.40 (IQR 234.00) | 4.52x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 300609.40 (IQR 1405.00) | 46804.20 (IQR 2497.60) | 6.42x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 319740.40 (IQR 1384.80) | 76358.60 (IQR 7203.40) | 4.19x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 132236.40 (IQR 95.20) | 61382.40 (IQR 1079.00) | 2.15x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3751997.00 (IQR 10461.20) | 430921.20 (IQR 1136.60) | 8.71x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 82426.20 (IQR 588.40) | 17440.00 (IQR 1211.00) | 4.73x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 354776.40 (IQR 20366.60) | 44224.00 (IQR 811.40) | 8.02x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 968329.60 (IQR 13558.80) | 139494.20 (IQR 5507.80) | 6.94x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1655607.60 (IQR 2559.60) | 227793.40 (IQR 5257.60) | 7.27x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Phase21 institutional 12-engine AAD vs FD | 480798.80 (IQR 7102.80) | 75527.20 (IQR 2208.40) | 6.37x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | Phase22 KRD/LPI/CoCo/Bergomi/CGMY/Guéant/Meucci/AAD | 2599317.20 (IQR 14766.20) | 679671.80 (IQR 3858.40) | 3.82x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | Phase23 RKHS/FMM/VV/LSV/JLT/CTD/QRH/Parisian/Aïd/Propagator/NCO/DML | 96574.00 (IQR 478.00) | 26159.80 (IQR 1162.20) | 3.69x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1374747.60 (IQR 827.60) | 158218.00 (IQR 1066.60) | 8.69x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Phase25 Kalman/RTS/VECM/HJ-SDF/IFRS9/CPT/Hasbrouck/FV/Lando/LPPL/CB-Arb/Causal/Malliavin | 139273.40 (IQR 261.60) | 17739.60 (IQR 1985.20) | 7.85x | 8.967e-6 |
| 38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc | Phase26 MF-Rates/Wildcard/CTD/Pivot-TARF/Bergomi/Schönbucher/FRTB-CVA/Lookback/VTA/FLAM/ACD/HLZ/MLMC | 1130565.20 (IQR 1356.00) | 48664.60 (IQR 1275.60) | 23.23x | 4.996e-9 |
| 39_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt | Phase27 SOTA FTPL/G2/eSSVI/Offshore-CB/CoCo/MPOR/Double-Heston/Curran/Warin/Roll-Corwin/Busseti-Kelly/Diff-QP | 694906.80 (IQR 2440.40) | 37491.80 (IQR 2979.00) | 18.53x | 4.423e-10 |
| 40_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca | Phase28 BCBS/IRRBB/Turfus-RFR/Baran-MtM/Bourgey-BG/FEDS-PONV/Zeron-SIMM/Jacquier-rBergomi/Choi-Spread/BESS-SDP/Kratz-Dark/Föllmer/Diff-PCA | 289120.80 (IQR 134.60) | 25369.80 (IQR 3776.80) | 11.40x | 1.803e-9 |
| 41_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is | Phase29 Veronesi/Lando/Gatheral/Russian/Bates/Koziol/MFG/BCBS/CJP/ETF/Qian/Glasserman | 702521.80 (IQR 681.60) | 85061.20 (IQR 1059.60) | 8.26x | 2.244e-9 |
| 42_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin | Phase30 Romer/Taylor/Damodaran/Christoffersen/Foucault/Barberis/Hart/Beneish/Campisi/Heckman/Elastic/ResetCB/FSB | 24239.00 (IQR 92.40) | 2385.20 (IQR 185.00) | 10.16x | 8.962e-10 |
| 43_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi | Phase31 book-backed SOTA engines | 71249.80 (IQR 448.80) | 5540.40 (IQR 291.40) | 12.86x | 1.010e-9 |
| 44_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan | Phase32 book-backed SOTA engines | 51943.00 (IQR 442.40) | 4423.80 (IQR 93.40) | 11.74x | 5.153e-8 |
| 45_rebonato_cms_leland_toft_calamos_ascot_spec_pool_roll_glosten_milgrom_potters_rie_egarch_hamilton_wooldridge_2sls_mcneil_euler_bacon_omega_damodaran_rd_mwg_arrow_debreu | Phase33 book-backed SOTA engines | 22264.80 (IQR 264.80) | 1710.40 (IQR 74.00) | 13.02x | 1.461e-9 |

Summary: 45 benchmark suites; speedups range from **2.08x** to **171.25x** (median **8.02x**, mean **18.89x**). All parity checks pass their domain tolerances; latest suite parity is **1.461e-9**.
