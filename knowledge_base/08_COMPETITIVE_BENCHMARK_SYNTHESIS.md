# Competitive Benchmark Synthesis

| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |
|---|---|---:|---:|---:|---:|
| 01_calendar_popcnt | day-by-day holidays vs. packed POPCNT calendar | 644775.60 (IQR 319673.80) | 16068.00 (IQR 6526.20) | 40.13x | 0.000e0 |
| 02_leg_simd_soa | boxed AoS virtual dispatch vs. contiguous SIMD SoA | 134694.40 (IQR 15612.00) | 7845.20 (IQR 323.20) | 17.17x | 1.490e-7 |
| 03_multicurve_ift_dv01 | 12 bump-and-rebootstrap valuations vs. cached IFT VJP | 472133.40 (IQR 260143.20) | 5406.00 (IQR 417.20) | 87.34x | 5.797e-5 |
| 04_implied_vol_ift_dual2 | nine Jäckel finite-difference solves vs. one Dual2 IFT solve | 96809.60 (IQR 7777.80) | 30660.00 (IQR 765.20) | 3.16x | 1.350e-5 |
| 05_heston_cos_dual5 | ten central COS valuations vs. one Dual5 COS valuation | 1972541.80 (IQR 56769.80) | 623163.20 (IQR 736.40) | 3.17x | 2.821e-4 |
| 06_gbm_mc_simd8 | scalar Philox path loop vs. 8-lane antithetic GBM | 439784.40 (IQR 31896.60) | 150089.40 (IQR 905.40) | 2.93x | 7.393e-11 |
| 07_scripted_autocall_fuzzy_aad | nine bumped scripted MC runs vs. one four-input fuzzy-AAD run | 811993.40 (IQR 46658.00) | 171243.20 (IQR 6147.00) | 4.74x | 2.914e-5 |
| 08_portfolio_dag_incremental | full bootstrap/revaluation vs. content-addressed incremental DAG | 2990560.00 (IQR 79210.20) | 256662.20 (IQR 70633.80) | 11.65x | 0.000e0 |
| 09_credit_cds_ift_cs01 | ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01 | 386122.80 (IQR 2875.00) | 16735.40 (IQR 1213.00) | 23.07x | 7.626e-10 |
| 10_fx_vanna_volga_dual3 | seven scalar FX smile quote bumps vs. one Dual3 smile pass | 62472.00 (IQR 3902.40) | 21709.60 (IQR 70.80) | 2.88x | 7.653e-7 |
| 11_capriotti_giles_cholesky_corr_aad | 66 correlation bumps vs. one adjoint Cholesky pass | 19346497.00 (IQR 908268.00) | 200257.40 (IQR 46819.00) | 96.61x | 8.326e-9 |
| 12_chebyshev_dynamic_simm_mva | brute-force dynamic margin vs. Chebyshev model-space margin | 1716185.20 (IQR 47297.20) | 209108.20 (IQR 47845.20) | 8.21x | 9.006e-10 |
| 13_reiner_rubinstein_barrier_dual5 | six finite-difference barrier bumps vs. one Dual5 analytical pass | 36303.40 (IQR 3911.40) | 11144.80 (IQR 478.80) | 3.26x | 1.854e-7 |
| 14_burley_owen_rqmc_variance_reduction | pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC | 13257669.20 (IQR 22952.00) | 1181286.80 (IQR 25941.60) | 11.22x | 2.719e-3 |
| 15_levy_cos_vg_dual5 | finite-difference VG COS strip Greeks vs. one Dual5 COS pass | 530881.60 (IQR 4789.60) | 124808.80 (IQR 2516.60) | 4.25x | 4.346e-7 |
| 16_asb_cdo_tranche_recursion_dual | eleven ASB loss recursions vs. one Dual5 tranche pass | 200003.40 (IQR 461.60) | 54314.60 (IQR 1532.60) | 3.68x | 2.850e-10 |
| 17_rfr_sabr_caplet_strip_dual4 | scalar RFR SABR strip bumps vs Dual4 effective-SABR strip | 27648.00 (IQR 311.80) | 5724.60 (IQR 550.80) | 4.83x | 8.401e-9 |
| 18_xccs_multicurrency_ift_dv01 | XCCS bump/bootstrap versus cached IFT basis risk | 1857689.80 (IQR 40791.80) | 10833.40 (IQR 199.60) | 171.48x | 4.775e-9 |
| 19_american_baw_ju_ift_dual5 | Ju/BAW scalar bumps versus one IFT Dual5 American pass | 515726.40 (IQR 35149.00) | 56066.60 (IQR 1780.00) | 9.20x | 5.363e-8 |
| 20_risk_parity_ift_allocation_aad | ninety covariance/budget bumps versus one IFT risk-parity solve | 7134205.60 (IQR 80318.40) | 145241.60 (IQR 2624.60) | 49.12x | 4.402e-9 |
| 21_bermudan_hw1f_tree_ift_dv01 | Bermudan Hull--White tree quote and parameter finite differences versus AD | 535520.00 (IQR 6125.80) | 263624.00 (IQR 43511.80) | 2.03x | 2.140e-9 |
| 22_sabr_svi_ift_and_qp_kkt_aad | SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD | 521116.80 (IQR 2739.20) | 8132.00 (IQR 1249.80) | 64.08x | 2.947e-7 |
| 23_slv_rbergomi_and_autocallable_cega_aad | SLV/rBergomi and worst-of autocallable sensitivities | 1195404.60 (IQR 51317.80) | 143973.80 (IQR 4323.40) | 8.30x | 2.027e-6 |
| 24_frtb_sbm_ima_and_cvar_execution_aad | FRTB/CVaR and optimal execution KKT sensitivities | 6446.40 (IQR 230.60) | 1413.00 (IQR 233.60) | 4.56x | 3.105e-7 |
| 25_monotone_convex_vix_parisian_robust_aad | monotone-convex curve plus VIX/Parisian and robust Markowitz AD | 309410.00 (IQR 41522.00) | 45926.20 (IQR 1618.00) | 6.74x | 1.502e-7 |
| 26_double_heston_tarf_coco_hrp_aad | Double Heston/TARF/CoCo/HRP single-pass AD | 316029.80 (IQR 4455.20) | 68277.40 (IQR 5014.60) | 4.63x | 1.206e-5 |
| 27_jy_inflation_tba_bergomi_nco_aad | JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD | 132924.60 (IQR 1577.20) | 57144.20 (IQR 1305.60) | 2.33x | 1.614e-7 |
| 28_cheyette_cmo_qrh_basecorr_cdar_aad | Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD | 3730820.80 (IQR 10318.20) | 423017.60 (IQR 1274.60) | 8.82x | 2.529e-9 |
| 29_asw_mf1f_ah_bess_ep_aad | ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD | 81174.40 (IQR 3148.80) | 16949.00 (IQR 1550.80) | 4.79x | 1.679e-9 |
| 30_mtmxccs_colva_zabr_pdv_weather_enb_aad | MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD | 311398.20 (IQR 43020.40) | 38527.20 (IQR 1299.80) | 8.08x | 2.986e-9 |
| 31_burnout_cmcds_quintic_afv_cs_schur_aad | Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD | 946515.20 (IQR 42575.40) | 136662.80 (IQR 2525.60) | 6.93x | 3.731e-9 |
| 32_sinhou_mb_rbergomi_choi_afs_expectile_aad | SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD | 1586399.60 (IQR 172533.40) | 221242.20 (IQR 12688.40) | 7.17x | 1.505e-8 |

All 32 suites pass numerical parity. The release exposes 30 workstation tabs; speedups range from 2.03x to 171.48x (7.05x median, 21.46x mean).
