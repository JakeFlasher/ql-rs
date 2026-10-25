# Competitive Benchmark Synthesis

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 321737.60 (IQR 116702.20) | 3745.80 (IQR 134.80) | 85.89x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 42682.00 (IQR 4728.40) | 2468.00 (IQR 26.60) | 17.29x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 252062.00 (IQR 11410.00) | 4570.40 (IQR 63.00) | 55.15x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 95885.80 (IQR 436.60) | 30812.80 (IQR 19.60) | 3.11x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1922537.20 (IQR 26284.00) | 619015.40 (IQR 365.80) | 3.11x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 439120.00 (IQR 5567.80) | 152305.00 (IQR 628.20) | 2.88x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 822766.00 (IQR 61022.40) | 175823.60 (IQR 4385.20) | 4.68x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 2973018.00 (IQR 223674.80) | 399627.60 (IQR 9614.20) | 7.44x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 385524.60 (IQR 433.60) | 17687.40 (IQR 3434.60) | 21.80x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 61814.60 (IQR 479.80) | 21481.80 (IQR 521.00) | 2.88x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19841582.00 (IQR 1826247.80) | 209804.60 (IQR 8381.60) | 94.57x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1771715.00 (IQR 48520.00) | 211607.20 (IQR 3105.60) | 8.37x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37368.40 (IQR 3082.40) | 11511.60 (IQR 858.80) | 3.25x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13648499.60 (IQR 471721.40) | 1201017.60 (IQR 9723.20) | 11.36x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 548818.80 (IQR 5015.00) | 129529.20 (IQR 1793.60) | 4.24x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 207355.00 (IQR 871.60) | 53867.40 (IQR 1052.60) | 3.85x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 29015.00 (IQR 303.20) | 6007.40 (IQR 214.20) | 4.83x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1946377.00 (IQR 6607.40) | 11342.40 (IQR 160.60) | 171.60x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 526594.60 (IQR 752.00) | 56589.20 (IQR 324.20) | 9.31x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7104226.40 (IQR 38680.40) | 149420.60 (IQR 3217.00) | 47.55x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 546031.00 (IQR 671.60) | 262094.40 (IQR 846.40) | 2.08x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 518720.60 (IQR 5199.80) | 8161.40 (IQR 1106.20) | 63.56x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 872834.40 (IQR 10554.20) | 135309.60 (IQR 3483.80) | 6.45x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6704.20 (IQR 865.80) | 1452.60 (IQR 193.20) | 4.62x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 312969.00 (IQR 6041.60) | 46894.60 (IQR 1892.20) | 6.67x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 324443.20 (IQR 548.00) | 70135.20 (IQR 2725.60) | 4.63x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 131879.00 (IQR 415.40) | 62494.20 (IQR 358.20) | 2.11x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3811919.60 (IQR 43123.40) | 432929.60 (IQR 819.80) | 8.80x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 81732.80 (IQR 1018.20) | 17385.60 (IQR 1330.20) | 4.70x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 285577.20 (IQR 6881.40) | 36151.00 (IQR 262.20) | 7.90x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 995302.80 (IQR 5311.60) | 142559.80 (IQR 955.80) | 6.98x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1600364.80 (IQR 57155.60) | 218985.80 (IQR 2017.00) | 7.31x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Phase21 institutional 12-engine AAD vs FD | 452014.80 (IQR 797.20) | 69531.20 (IQR 2036.40) | 6.50x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | Phase22 KRD/LPI/CoCo/Bergomi/CGMY/Guéant/Meucci/AAD | 2547111.00 (IQR 8574.00) | 669240.20 (IQR 7098.00) | 3.81x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | Phase23 RKHS/FMM/VV/LSV/JLT/CTD/QRH/Parisian/Aïd/Propagator/NCO/DML | 94881.00 (IQR 378.80) | 26021.60 (IQR 1498.00) | 3.65x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1344790.40 (IQR 25506.00) | 162776.40 (IQR 536.60) | 8.26x | 2.234e-7 |
| 37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin | Phase25 Kalman/RTS/VECM/HJ-SDF/IFRS9/CPT/Hasbrouck/FV/Lando/LPPL/CB-Arb/Causal/Malliavin | 9800.00 (IQR 42.00) | 3629.60 (IQR 304.20) | 2.70x | 7.588e-6 |

All 37 suites pass numerical parity within domain tolerances. The release exposes 35 workstation tabs; speedups range from 2.08x to 171.60x (6.50x median, 19.29x mean).
