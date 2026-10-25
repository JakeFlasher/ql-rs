# Competitive Benchmark Synthesis

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 489434.20 (IQR 891382.40) | 7138.40 (IQR 2014.60) | 68.56x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 85306.40 (IQR 3027.40) | 5315.80 (IQR 80.40) | 16.05x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 655100.80 (IQR 303835.80) | 7993.80 (IQR 199.60) | 81.95x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 137648.80 (IQR 25920.60) | 38954.60 (IQR 6401.60) | 3.53x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1962367.80 (IQR 1538.40) | 749907.80 (IQR 801.80) | 2.62x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 446274.40 (IQR 16820.40) | 152585.00 (IQR 1735.60) | 2.92x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 841993.80 (IQR 6265.60) | 178456.20 (IQR 4262.00) | 4.72x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 3058976.40 (IQR 168572.20) | 283266.40 (IQR 100502.40) | 10.80x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 385353.40 (IQR 90073.60) | 16960.20 (IQR 628.60) | 22.72x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 61750.40 (IQR 535.00) | 21658.40 (IQR 438.00) | 2.85x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19688589.60 (IQR 61634.40) | 207529.40 (IQR 8269.80) | 94.87x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1776110.00 (IQR 30924.80) | 217008.40 (IQR 2160.40) | 8.18x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 37268.20 (IQR 2999.40) | 11534.60 (IQR 420.40) | 3.23x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13785416.00 (IQR 155144.80) | 1196479.40 (IQR 4898.20) | 11.52x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 546066.00 (IQR 461.20) | 127678.00 (IQR 1330.00) | 4.28x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 207360.60 (IQR 629.60) | 54054.00 (IQR 744.60) | 3.84x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 28770.40 (IQR 288.80) | 5917.80 (IQR 207.00) | 4.86x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1960101.20 (IQR 33437.60) | 11887.60 (IQR 303.00) | 164.89x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 527058.20 (IQR 10206.60) | 56583.40 (IQR 1675.00) | 9.31x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7577865.40 (IQR 85079.60) | 163455.80 (IQR 3971.20) | 46.36x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 544921.60 (IQR 1243.80) | 262186.00 (IQR 2270.00) | 2.08x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 514546.60 (IQR 2691.80) | 8225.20 (IQR 836.40) | 62.56x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 806658.40 (IQR 12280.80) | 142158.20 (IQR 3508.00) | 5.67x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 7155.60 (IQR 745.80) | 1548.00 (IQR 258.60) | 4.62x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 306709.20 (IQR 3729.60) | 46880.00 (IQR 1752.20) | 6.54x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 325022.40 (IQR 4988.20) | 78079.00 (IQR 3211.80) | 4.16x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 133864.80 (IQR 326.80) | 62167.20 (IQR 589.00) | 2.15x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3816840.00 (IQR 77254.80) | 445093.00 (IQR 5932.80) | 8.58x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 87292.80 (IQR 613.20) | 18276.00 (IQR 1081.20) | 4.78x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 365328.20 (IQR 778.20) | 45468.60 (IQR 270.20) | 8.03x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 989606.80 (IQR 625.60) | 142934.40 (IQR 2226.60) | 6.92x | 3.734e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1683578.20 (IQR 39020.40) | 232479.80 (IQR 1194.80) | 7.24x | 1.505e-8 |
| 33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad | Phase21 institutional 12-engine AAD vs FD | 487870.60 (IQR 277.80) | 76973.60 (IQR 1172.00) | 6.34x | 6.899e-9 |
| 34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad | Phase22 KRD/LPI/CoCo/Bergomi/CGMY/Guéant/Meucci/AAD | 2643024.00 (IQR 8449.00) | 692947.40 (IQR 5941.20) | 3.81x | 1.979e-9 |
| 35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml | Phase23 RKHS/FMM/VV/LSV/JLT/CTD/QRH/Parisian/Aïd/Propagator/NCO/DML | 99267.20 (IQR 331.80) | 27946.60 (IQR 1076.60) | 3.55x | 1.167e-9 |
| 36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp | MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD | 1387318.60 (IQR 16444.40) | 164746.20 (IQR 2255.40) | 8.42x | 2.234e-7 |

All 36 suites pass numerical parity. The release exposes 34 workstation tabs; speedups range from 2.08x to 164.89x (6.44x median, 19.82x mean).
