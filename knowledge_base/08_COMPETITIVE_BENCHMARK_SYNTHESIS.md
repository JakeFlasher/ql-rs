# Competitive Benchmark Synthesis

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 335632.20 (IQR 13311.60) | 3100.60 (IQR 67.60) | 108.25x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 34469.20 (IQR 4872.00) | 2158.20 (IQR 11.60) | 15.97x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 289480.00 (IQR 195438.20) | 8991.60 (IQR 652.20) | 32.19x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 218306.20 (IQR 38773.80) | 85275.20 (IQR 4850.00) | 2.56x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1991062.60 (IQR 297732.20) | 609861.80 (IQR 8641.00) | 3.26x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 430557.60 (IQR 5593.20) | 149081.60 (IQR 175.40) | 2.89x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 813515.20 (IQR 19477.60) | 174810.00 (IQR 5420.40) | 4.65x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 2937743.00 (IQR 128511.20) | 302097.60 (IQR 121501.60) | 9.72x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 631033.40 (IQR 448953.20) | 19166.20 (IQR 3414.40) | 32.92x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 61935.80 (IQR 503.20) | 21489.00 (IQR 39.80) | 2.88x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19490726.00 (IQR 235931.20) | 203563.60 (IQR 7749.40) | 95.75x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1756743.20 (IQR 24418.20) | 214420.80 (IQR 5354.60) | 8.19x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 36334.60 (IQR 2590.60) | 11607.20 (IQR 476.00) | 3.13x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13365095.00 (IQR 37892.60) | 1176616.00 (IQR 11977.00) | 11.36x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 541153.20 (IQR 6560.80) | 127978.00 (IQR 1574.60) | 4.23x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 209154.00 (IQR 7810.00) | 54562.80 (IQR 1021.40) | 3.83x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28210.00 (IQR 1275.40) | 5877.20 (IQR 179.20) | 4.80x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1908819.00 (IQR 12072.00) | 11012.60 (IQR 224.00) | 173.33x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 510236.20 (IQR 3247.60) | 54204.00 (IQR 248.40) | 9.41x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7291856.60 (IQR 10841.60) | 155216.20 (IQR 9636.40) | 46.98x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 543252.60 (IQR 688.60) | 261465.60 (IQR 730.20) | 2.08x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 516046.80 (IQR 18970.60) | 8318.20 (IQR 1508.60) | 62.04x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 1024429.40 (IQR 17224.20) | 134525.60 (IQR 4147.40) | 7.62x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6799.80 (IQR 573.60) | 1443.40 (IQR 305.80) | 4.71x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 307737.40 (IQR 5146.60) | 49825.00 (IQR 8736.80) | 6.18x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 309479.40 (IQR 677.80) | 74524.80 (IQR 3267.40) | 4.15x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 130928.60 (IQR 16940.40) | 61667.20 (IQR 3611.60) | 2.12x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3771129.60 (IQR 8771.60) | 441336.80 (IQR 8404.80) | 8.54x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 86358.40 (IQR 338.60) | 17577.40 (IQR 1029.00) | 4.91x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 307835.80 (IQR 2352.80) | 39534.00 (IQR 713.80) | 7.79x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 983291.80 (IQR 13712.60) | 139771.20 (IQR 2235.00) | 7.04x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1615914.00 (IQR 36066.80) | 222826.40 (IQR 3318.20) | 7.25x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Monotone-Convex / Clark-FX / SA-CCR-KVA / rVol / Cartea / Paleologo AD | 478295.60 (IQR 6512.20) | 76932.60 (IQR 2143.40) | 6.22x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | KRD / LPI / CoCo / Bergomi / CGMY / Guéant / Meucci AD vs FD | 2565844.40 (IQR 273873.40) | 675272.80 (IQR 27047.00) | 3.80x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | RKHS / FMM / VV / LSV / JLT / CTD / QRH / Parisian / Aïd / Propagator / NCO / DML AD | 92866.20 (IQR 1463.00) | 26213.60 (IQR 1444.00) | 3.54x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1344872.60 (IQR 13800.80) | 156780.00 (IQR 339.00) | 8.58x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Kalman-RTS / VECM / HJ-SDF / IFRS9 / CPT / Hasbrouck / FV / Lando / LPPL / CB-Arb / Causal / Malliavin AD | 136332.20 (IQR 1307.20) | 15656.00 (IQR 409.80) | 8.71x | 8.967e-6 |
| 38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc | MF-Rates / Wildcard / CTD / Pivot-TARF / Bergomi / Schönbucher / FRTB-CVA / Lookback / VTA / FLAM / ACD / HLZ / MLMC AD | 1102723.40 (IQR 11599.60) | 49627.20 (IQR 5973.40) | 22.22x | 4.996e-9 |
| 39_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt | FTPL / G2++ / eSSVI / Offshore-CB / CoCo / MPOR / Double-Heston / Curran / Warin / Roll-Corwin / Busseti-Kelly / Diff-QP AD | 661964.60 (IQR 3337.40) | 36853.40 (IQR 2220.80) | 17.96x | 4.423e-10 |
| 40_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca | BCBS-IRRBB / Turfus-RFR / Baran-MtM / Bourgey-BG / FEDS-PONV / Zeron-SIMM / Jacquier-rBergomi / Choi-Spread / BESS-SDP / Kratz-Dark / Föllmer / Diff-PCA AD | 267992.40 (IQR 4416.20) | 24900.80 (IQR 1898.40) | 10.76x | 1.803e-9 |
| 41_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is | Veronesi / Lando / Gatheral-Corridor / Russian / Bates-SVJ / Koziol-CB / MFG / BCBS-Floor / CJP-LOB / ETF / Qian-IC / Glasserman-IS AD | 676684.60 (IQR 16384.80) | 82632.60 (IQR 1730.40) | 8.19x | 2.244e-9 |
| 42_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin | Romer-Taylor / Damodaran-CRP / Cumulant-VaR / LOB-LVaR / Prospect / Regret / Beneish / Campisi / Heckman / Elastic-Net / Reset-CB / FSB-BailIn AD | 21250.00 (IQR 428.40) | 2352.20 (IQR 276.20) | 9.03x | 8.962e-10 |
| 43_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi | KRD-LSC / CMO-PAC / Asquith-CB / ASU2024-CPDI / ASC470-Mod / ColVA-CTD / Hill-EVT / Bouchaud-TLF / Karnosky-FX / DGTW-Timing / Brock-Hommes / Penman-ReOI AD | 60151.60 (IQR 1132.60) | 5488.60 (IQR 214.40) | 10.96x | 1.010e-9 |
| 44_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan | Cheyette-LV / ADB-NDIRS / Andreasen-ZABR / QRH-Zumbach / CMCDS-Morini / Turlakov-WWR / HHW-ADI / CappedCall-MW / Ahmadi-EVaR / Schur-HRP / BEKK-GARCH / Schilit-Sloan AD | 41827.80 (IQR 46.20) | 4400.60 (IQR 143.40) | 9.51x | 5.153e-8 |
| 45_rebonato_cms_leland_toft_calamos_ascot_spec_pool_roll_glosten_milgrom_potters_rie_egarch_hamilton_wooldridge_2sls_mcneil_euler_bacon_omega_damodaran_rd_mwg_arrow_debreu | Rebonato-CMS / Leland-Toft / Calamos-Ascot / SpecPool-Roll / Glosten-Milgrom / Potters-RIE / EGARCH-Hamilton / Wooldridge-2SLS / McNeil-Euler / Bacon-Omega / Damodaran-RD / Arrow-Debreu AD | 14439.20 (IQR 323.60) | 1666.20 (IQR 47.40) | 8.67x | 6.115e-9 |
| 46_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig | Tuckman-Repo / BK-FVA / DECS / MBS-IO-PO / Kyle-PIN / GP-Aim / GJR-DCC / LPPL / Spectral-Risk / Brinson-Carino / Prelec-CPT / Path-Sig AD | 3.60 (IQR 5.40) | 3.40 (IQR 4.40) | 1.06x | 0.000e0 |

Summary: 46 benchmark suites; speedups range from **1.06x** to **173.33x** (median **8.19x**, mean **18.04x**). All parity checks pass their domain tolerances; latest suite parity is **0.000e0**.
