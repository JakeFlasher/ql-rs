# Competitive Benchmark Synthesis

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 669940.60 (IQR 345552.60) | 6488.80 (IQR 319.80) | 103.25x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 66640.40 (IQR 12832.80) | 3915.60 (IQR 20.40) | 17.02x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 303605.60 (IQR 131438.80) | 5224.40 (IQR 620.80) | 58.11x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 95502.20 (IQR 137.80) | 30665.40 (IQR 337.00) | 3.11x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1990999.40 (IQR 51780.60) | 710170.20 (IQR 8733.60) | 2.80x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 428092.80 (IQR 9780.80) | 148743.20 (IQR 462.00) | 2.88x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 804884.60 (IQR 2800.20) | 197874.00 (IQR 436.20) | 4.07x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 2960121.20 (IQR 134008.80) | 292441.40 (IQR 43571.40) | 10.12x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 386849.40 (IQR 422182.60) | 16765.60 (IQR 695.60) | 23.07x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 62339.60 (IQR 316.40) | 21534.80 (IQR 33.20) | 2.89x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19467574.80 (IQR 285126.20) | 204053.40 (IQR 10174.00) | 95.40x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1767289.60 (IQR 11211.60) | 214320.40 (IQR 3218.20) | 8.25x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37749.80 (IQR 3454.60) | 13239.60 (IQR 632.00) | 2.85x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13588369.20 (IQR 277237.60) | 1205663.40 (IQR 2239.00) | 11.27x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 543640.20 (IQR 426.80) | 129428.40 (IQR 709.20) | 4.20x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 207499.80 (IQR 407.40) | 54039.00 (IQR 874.40) | 3.84x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28750.00 (IQR 349.40) | 6204.20 (IQR 210.20) | 4.63x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1934912.40 (IQR 1575.00) | 11566.60 (IQR 266.00) | 167.28x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 526825.60 (IQR 1329.20) | 56624.20 (IQR 284.20) | 9.30x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7071046.60 (IQR 89517.20) | 149084.00 (IQR 3558.20) | 47.43x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 538363.60 (IQR 285.00) | 257893.60 (IQR 754.20) | 2.09x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 501151.20 (IQR 2513.20) | 8010.80 (IQR 1077.40) | 62.56x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 791148.60 (IQR 8271.60) | 140135.80 (IQR 4725.80) | 5.65x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6575.80 (IQR 314.00) | 1470.80 (IQR 572.80) | 4.47x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 302768.40 (IQR 3423.00) | 44639.40 (IQR 3702.80) | 6.78x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 316652.80 (IQR 911.60) | 76410.00 (IQR 3661.60) | 4.14x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 130098.60 (IQR 1437.80) | 61912.40 (IQR 509.40) | 2.10x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3756181.80 (IQR 4494.20) | 430576.20 (IQR 2003.60) | 8.72x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 84135.00 (IQR 250.40) | 17570.00 (IQR 1009.00) | 4.79x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 303841.80 (IQR 8502.20) | 37858.80 (IQR 814.80) | 8.03x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 960115.40 (IQR 3052.60) | 137881.80 (IQR 2170.40) | 6.96x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1638212.40 (IQR 16498.00) | 227982.60 (IQR 2040.40) | 7.19x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Phase21 institutional 12-engine AAD vs FD | 485709.00 (IQR 367.20) | 76676.20 (IQR 1289.60) | 6.33x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | Phase22 KRD/LPI/CoCo/Bergomi/CGMY/Guéant/Meucci/AAD | 2582850.20 (IQR 4743.80) | 676091.40 (IQR 3418.20) | 3.82x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | Phase23 RKHS/FMM/VV/LSV/JLT/CTD/QRH/Parisian/Aïd/Propagator/NCO/DML | 95275.40 (IQR 1173.20) | 26789.20 (IQR 1225.60) | 3.56x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1359747.20 (IQR 1085.80) | 167247.20 (IQR 992.00) | 8.13x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Phase25 Kalman/RTS/VECM/HJ-SDF/IFRS9/CPT/Hasbrouck/FV/Lando/LPPL/CB-Arb/Causal/Malliavin | 138370.40 (IQR 789.40) | 15831.00 (IQR 351.60) | 8.74x | 8.967e-6 |
| 38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc | Phase26 MF-Rates/Wildcard/CTD/Pivot-TARF/Bergomi/Schönbucher/FRTB-CVA/Lookback/VTA/FLAM/ACD/HLZ/MLMC | 1101132.40 (IQR 3079.40) | 48199.20 (IQR 1173.80) | 22.85x | 4.996e-9 |
| 39_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt | Phase27 SOTA FTPL/G2/eSSVI/Offshore-CB/CoCo/MPOR/Double-Heston/Curran/Warin/Roll-Corwin/Busseti-Kelly/Diff-QP | 652693.20 (IQR 4712.40) | 37204.80 (IQR 973.60) | 17.54x | 4.423e-10 |
| 40_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca | Phase28 BCBS/IRRBB/Turfus-RFR/Baran-MtM/Bourgey-BG/FEDS-PONV/Zeron-SIMM/Jacquier-rBergomi/Choi-Spread/BESS-SDP/Kratz-Dark/Föllmer/Diff-PCA | 271196.80 (IQR 2582.20) | 24704.40 (IQR 789.60) | 10.98x | 1.803e-9 |
| 41_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is | Phase29 Veronesi/Lando/Gatheral/Russian/Bates/Koziol/MFG/BCBS/CJP/ETF/Qian/Glasserman | 684832.60 (IQR 716.80) | 87164.40 (IQR 954.20) | 7.86x | 2.244e-9 |

Summary: 41 benchmark suites; speedups range from **2.09x** to **167.28x** (median **7.19x**, mean **19.39x**). All parity checks pass their domain tolerances; latest suite parity is **2.244e-9**.
