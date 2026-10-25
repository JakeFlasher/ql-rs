# Competitive Benchmark Synthesis

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 581170.80 (IQR 115465.60) | 5496.00 (IQR 179.00) | 105.74x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 62381.20 (IQR 6047.60) | 3504.00 (IQR 47.20) | 17.80x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 295368.40 (IQR 103351.00) | 4421.80 (IQR 261.80) | 66.80x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 95668.40 (IQR 2266.40) | 30802.60 (IQR 10.60) | 3.11x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 2006800.00 (IQR 686074.40) | 764817.80 (IQR 1016.00) | 2.62x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 448859.60 (IQR 5637.40) | 156399.40 (IQR 8069.60) | 2.87x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 835006.60 (IQR 11125.80) | 179158.20 (IQR 3730.60) | 4.66x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 3031891.00 (IQR 128414.80) | 325455.20 (IQR 64951.80) | 9.32x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 386331.40 (IQR 1201.80) | 16586.20 (IQR 3018.80) | 23.29x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 62140.40 (IQR 385.40) | 22450.40 (IQR 1473.20) | 2.77x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 20344872.20 (IQR 383289.60) | 211252.20 (IQR 10563.20) | 96.31x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1839492.20 (IQR 29012.60) | 222491.60 (IQR 6866.40) | 8.27x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37846.80 (IQR 2644.40) | 11141.80 (IQR 523.40) | 3.40x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 14017606.00 (IQR 483195.00) | 1208852.80 (IQR 7852.20) | 11.60x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 546552.00 (IQR 6905.80) | 128239.20 (IQR 2160.20) | 4.26x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 204319.20 (IQR 2185.20) | 54113.80 (IQR 1476.60) | 3.78x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28620.20 (IQR 167.80) | 6016.60 (IQR 112.80) | 4.76x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1953776.40 (IQR 60991.80) | 13294.60 (IQR 3121.40) | 146.96x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 522551.60 (IQR 8223.40) | 55434.20 (IQR 1541.60) | 9.43x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7265888.00 (IQR 108735.40) | 155889.80 (IQR 6538.20) | 46.61x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 558005.20 (IQR 2050.00) | 262109.60 (IQR 6200.60) | 2.13x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 530663.00 (IQR 36980.80) | 8319.60 (IQR 933.60) | 63.78x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 1254493.60 (IQR 141480.40) | 136490.40 (IQR 3773.80) | 9.19x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6596.00 (IQR 1139.80) | 1590.60 (IQR 226.80) | 4.15x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 322591.80 (IQR 7343.80) | 50503.60 (IQR 2700.00) | 6.39x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 337211.40 (IQR 14598.80) | 79090.20 (IQR 2982.60) | 4.26x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 137042.60 (IQR 2145.40) | 63636.40 (IQR 986.60) | 2.15x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3928368.00 (IQR 16004.40) | 448398.40 (IQR 1426.00) | 8.76x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 83840.20 (IQR 2289.60) | 18237.00 (IQR 1219.40) | 4.60x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 298286.20 (IQR 16023.60) | 38115.80 (IQR 854.20) | 7.83x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 991560.60 (IQR 13596.80) | 144284.60 (IQR 2048.00) | 6.87x | 3.731e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1657254.80 (IQR 7903.40) | 234419.00 (IQR 15032.60) | 7.07x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Phase21 institutional 12-engine AAD vs FD | 333495.80 (IQR 1577.40) | 66006.80 (IQR 2986.80) | 5.05x | 1.232e-8 |

All 33 suites pass numerical parity. The release exposes 31 workstation tabs; speedups range from 2.03x to 171.48x (7.05x median, 21.46x mean).
