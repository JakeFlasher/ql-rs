# Competitive Benchmark Synthesis

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 203887.40 (IQR 232.40) | 3033.60 (IQR 118.00) | 67.21x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 34491.00 (IQR 4289.20) | 2160.80 (IQR 42.40) | 15.96x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 254040.00 (IQR 2175.00) | 4821.00 (IQR 477.00) | 52.69x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 95960.00 (IQR 422.20) | 30909.20 (IQR 93.60) | 3.10x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 2000360.60 (IQR 30160.80) | 622555.60 (IQR 6381.20) | 3.21x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 440086.00 (IQR 4281.20) | 151196.00 (IQR 2760.40) | 2.91x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 804611.40 (IQR 12698.00) | 188336.20 (IQR 3770.40) | 4.27x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 2956727.80 (IQR 57410.40) | 347873.20 (IQR 83437.20) | 8.50x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 385022.40 (IQR 572.00) | 17910.60 (IQR 3079.60) | 21.50x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 61871.00 (IQR 326.60) | 21581.40 (IQR 34.80) | 2.87x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19548986.40 (IQR 133998.60) | 205986.20 (IQR 9659.20) | 94.90x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1771237.20 (IQR 11691.20) | 214688.20 (IQR 2068.80) | 8.25x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37191.60 (IQR 2278.20) | 13209.40 (IQR 871.80) | 2.82x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13577033.00 (IQR 19029.60) | 1196709.80 (IQR 2638.80) | 11.35x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 548092.40 (IQR 249.00) | 128603.60 (IQR 1109.60) | 4.26x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 206730.20 (IQR 419.40) | 54450.60 (IQR 782.40) | 3.80x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28772.00 (IQR 557.80) | 6013.60 (IQR 210.80) | 4.78x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1917364.20 (IQR 1396.60) | 11494.00 (IQR 188.60) | 166.81x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 524423.00 (IQR 805.60) | 56645.80 (IQR 237.60) | 9.26x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 6973607.40 (IQR 81265.60) | 151901.20 (IQR 4558.00) | 45.91x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 548642.80 (IQR 95.80) | 262068.60 (IQR 2689.80) | 2.09x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 511449.40 (IQR 6647.40) | 8068.80 (IQR 1105.80) | 63.39x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 822416.20 (IQR 2874.40) | 133896.00 (IQR 3070.20) | 6.14x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6628.80 (IQR 1032.00) | 1534.00 (IQR 270.20) | 4.32x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 311709.80 (IQR 6872.20) | 48142.80 (IQR 1985.20) | 6.47x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 327940.60 (IQR 724.60) | 69790.00 (IQR 7717.40) | 4.70x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 134674.60 (IQR 249.80) | 64007.40 (IQR 2540.80) | 2.10x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3821320.40 (IQR 27185.60) | 438505.60 (IQR 1111.60) | 8.71x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 83650.80 (IQR 3977.40) | 17422.60 (IQR 1089.00) | 4.80x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 328837.80 (IQR 4116.00) | 41200.20 (IQR 535.80) | 7.98x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 976274.60 (IQR 508.00) | 140448.00 (IQR 5646.20) | 6.95x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1655892.80 (IQR 1328.00) | 229454.00 (IQR 1430.00) | 7.22x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Phase21 institutional 12-engine AAD vs FD | 495141.00 (IQR 165.40) | 76739.40 (IQR 1029.20) | 6.45x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | Phase22 KRD/LPI/CoCo/Bergomi/CGMY/Guéant/Meucci/AAD | 2601296.60 (IQR 8922.60) | 677920.00 (IQR 4784.20) | 3.84x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | Phase23 RKHS/FMM/VV/LSV/JLT/CTD/QRH/Parisian/Aïd/Propagator/NCO/DML | 95309.20 (IQR 1298.00) | 25790.40 (IQR 1077.40) | 3.70x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1362364.60 (IQR 11972.40) | 161222.80 (IQR 1361.00) | 8.45x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Phase25 Kalman/RTS/VECM/HJ-SDF/IFRS9/CPT/Hasbrouck/FV/Lando/LPPL/CB-Arb/Causal/Malliavin | 141719.80 (IQR 420.80) | 16121.00 (IQR 1568.60) | 8.79x | 8.967e-6 |
| 38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc | Phase26 MF-Rates/Wildcard/CTD/Pivot-TARF/Bergomi/Schönbucher/FRTB-CVA/Lookback/VTA/FLAM/ACD/HLZ/MLMC | 1144013.20 (IQR 7207.20) | 50410.20 (IQR 1137.40) | 22.69x | 4.996e-9 |
| 39_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt | Phase27 SOTA FTPL/G2/eSSVI/Offshore-CB/CoCo/MPOR/Double-Heston/Curran/Warin/Roll-Corwin/Busseti-Kelly/Diff-QP | 678868.60 (IQR 13089.00) | 37516.00 (IQR 3599.20) | 18.10x | 4.423e-10 |
| 40_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca | Phase28 BCBS/IRRBB/Turfus-RFR/Baran-MtM/Bourgey-BG/FEDS-PONV/Zeron-SIMM/Jacquier-rBergomi/Choi-Spread/BESS-SDP/Kratz-Dark/Föllmer/Diff-PCA | 15080.40 (IQR 34.20) | 1933.80 (IQR 38.60) | 7.80x | 9.999e-8 |

Summary: 40 benchmark suites; speedups range from **2.09x** to **166.81x** (median **7.22x**, mean **18.48x**). All parity checks pass their domain tolerances; latest suite parity is **9.999e-8**.
