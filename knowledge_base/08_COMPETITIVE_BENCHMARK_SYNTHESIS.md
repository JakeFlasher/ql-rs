# Competitive Benchmark Synthesis

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 599078.80 (IQR 464864.00) | 5671.20 (IQR 590.40) | 105.64x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 55966.20 (IQR 8162.00) | 3389.80 (IQR 581.60) | 16.51x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 258068.40 (IQR 70100.80) | 4909.20 (IQR 466.60) | 52.57x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 96213.40 (IQR 568.60) | 30653.40 (IQR 40.40) | 3.14x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1925010.60 (IQR 31568.00) | 719332.20 (IQR 22826.00) | 2.68x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 425916.20 (IQR 6435.40) | 149412.00 (IQR 406.80) | 2.85x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 790382.20 (IQR 6373.60) | 175694.20 (IQR 2991.20) | 4.50x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 3076423.00 (IQR 196590.00) | 385530.20 (IQR 10133.20) | 7.98x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 385446.20 (IQR 1792.80) | 19351.80 (IQR 3503.00) | 19.92x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 62125.60 (IQR 392.00) | 21525.40 (IQR 2.20) | 2.89x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19534930.00 (IQR 274488.20) | 203478.40 (IQR 11069.80) | 96.00x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1767838.20 (IQR 13294.60) | 214744.60 (IQR 2328.60) | 8.23x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37793.00 (IQR 2592.80) | 13252.40 (IQR 451.60) | 2.85x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13584031.20 (IQR 24571.40) | 1222073.40 (IQR 4123.80) | 11.12x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 556760.40 (IQR 618.60) | 131805.80 (IQR 739.00) | 4.22x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 206697.20 (IQR 1894.00) | 57993.40 (IQR 385.40) | 3.56x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28907.20 (IQR 336.40) | 5899.00 (IQR 175.60) | 4.90x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1918827.00 (IQR 5553.00) | 11807.40 (IQR 257.60) | 162.51x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 528255.20 (IQR 1322.00) | 56829.80 (IQR 421.80) | 9.30x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7073734.40 (IQR 34539.60) | 148396.80 (IQR 2730.80) | 47.67x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 546562.00 (IQR 967.40) | 261431.40 (IQR 1085.20) | 2.09x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 515453.60 (IQR 3598.40) | 8257.40 (IQR 861.40) | 62.42x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 803382.00 (IQR 4552.20) | 134004.00 (IQR 4732.80) | 6.00x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6715.20 (IQR 187.20) | 1443.00 (IQR 237.20) | 4.65x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 308159.40 (IQR 1469.80) | 47792.60 (IQR 945.80) | 6.45x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 325127.20 (IQR 671.00) | 77013.40 (IQR 5159.40) | 4.22x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 132911.60 (IQR 187.20) | 62197.60 (IQR 2144.60) | 2.14x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3762499.40 (IQR 2389.40) | 433113.40 (IQR 809.80) | 8.69x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 86086.20 (IQR 3903.20) | 17569.80 (IQR 877.00) | 4.90x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 295175.80 (IQR 4515.80) | 36866.80 (IQR 481.00) | 8.01x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 968071.80 (IQR 861.00) | 139046.60 (IQR 2886.80) | 6.96x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1622736.40 (IQR 4690.00) | 226494.20 (IQR 1702.80) | 7.16x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Phase21 institutional 12-engine AAD vs FD | 483923.60 (IQR 1568.20) | 76509.40 (IQR 1751.40) | 6.33x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | Phase22 KRD/LPI/CoCo/Bergomi/CGMY/Guéant/Meucci/AAD | 2591673.80 (IQR 6401.20) | 675151.60 (IQR 11078.40) | 3.84x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | Phase23 RKHS/FMM/VV/LSV/JLT/CTD/QRH/Parisian/Aïd/Propagator/NCO/DML | 93608.60 (IQR 618.40) | 26501.60 (IQR 2944.00) | 3.53x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1372567.40 (IQR 9626.80) | 157976.80 (IQR 4695.00) | 8.69x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Phase25 Kalman/RTS/VECM/HJ-SDF/IFRS9/CPT/Hasbrouck/FV/Lando/LPPL/CB-Arb/Causal/Malliavin | 139085.40 (IQR 349.60) | 17675.00 (IQR 1763.80) | 7.87x | 8.967e-6 |
| 38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc | Phase26 MF-Rates/Wildcard/CTD/Pivot-TARF/Bergomi/Schönbucher/FRTB-CVA/Lookback/VTA/FLAM/ACD/HLZ/MLMC | 1121599.40 (IQR 972.40) | 48922.20 (IQR 1077.40) | 22.93x | 4.996e-9 |

Summary: 38 benchmark suites; speedups range from **2.09x** to **162.51x** (median **6.96x**, mean **19.63x**). All parity checks pass their domain tolerances; Suite 38 parity is **4.996e-9**.
