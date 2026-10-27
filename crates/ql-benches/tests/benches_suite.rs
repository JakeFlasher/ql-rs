use ql_benches::harness::with_config_override;
use ql_benches::run_all_benchmark_suites_with_config;

#[test]
fn benchmarks_verify_parity_and_optional_sota_speedups() {
    // A short deterministic run is enough to exercise all 20 numerical
    // parity paths during normal `cargo test`; timing gates are inherently
    // host-dependent and opt in through QL_RUN_TIMING_GATES=1.
    let rows = run_all_benchmark_suites_with_config(1, 3);
    assert_eq!(rows.len(), 77);
    for row in &rows {
        assert!(!row.id.is_empty());
        assert!(!row.description.is_empty());
        assert!(!row.category.is_empty());
        assert!(row.baseline_median_ns.is_finite() && row.baseline_median_ns > 0.0);
        assert!(row.sota_median_ns.is_finite() && row.sota_median_ns > 0.0);
        assert!(row.speedup.is_finite() && row.speedup > 0.0);
        assert!(row.parity_ok, "{} parity flag is false", row.id);
        assert!(row.max_abs_diff.is_finite());
        assert!(
            row.max_abs_diff <= row.parity_tol,
            "{} parity diff={} tol={}",
            row.id,
            row.max_abs_diff,
            row.parity_tol
        );
    }

    if std::env::var("QL_RUN_TIMING_GATES").is_ok() {
        let log_mean = rows
            .iter()
            .map(|row| row.speedup.max(1e-12).ln())
            .sum::<f64>()
            / rows.len() as f64;
        assert!(log_mean.exp() >= 1.50, "geomean speedup={}", log_mean.exp());
        for row in &rows {
            assert!(row.speedup >= 1.05, "{} speedup={}", row.id, row.speedup);
        }
    }
}

#[test]
fn suite34_suite34_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("34_"))
        .expect("suite 34");
    assert!(row.parity_ok);
}
#[test]
fn suite35_suite35_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("35_"))
        .expect("suite 35");
    assert!(row.parity_ok);
}

#[test]
fn suite36_suite36_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("36_"))
        .expect("suite 36");
    assert!(row.parity_ok);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite37_suite37_workload_is_registered() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin(1)
    });
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite38_suite38_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("38_"))
        .expect("suite 38");
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite39_suite39_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("39_"))
        .expect("suite 39");
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite40_suite40_workload_is_registered() {
    let rows = run_all_benchmark_suites_with_config(1, 3);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("40_"))
        .expect("suite 40");
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite41_suite41_workload_is_registered() {
    let rows = ql_benches::run_all_benchmarks(1);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("41_"))
        .expect("suite 41");
    assert!(row.max_abs_diff < 1e-4);
}

#[test]
fn suite42_suite42_workload_is_registered() {
    let rows = ql_benches::run_all_benchmarks(1);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("42_"))
        .expect("suite 42");
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
}

#[test]
fn suite43_suite43_workload_is_registered() {
    let rows = ql_benches::run_all_benchmarks(1);
    let row = rows
        .iter()
        .find(|r| r.id.starts_with("43_"))
        .expect("suite 43");
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
}

#[test]
fn suite44_suite44_workload_is_registered() {
    let rows = ql_benches::run_all_benchmarks(1);
    assert!(rows.iter().any(|r| r.id.starts_with("44_")));
    let row = ql_benches::benchmark_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite45_suite45_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("45_")));
    let row = ql_benches::benchmark_rebonato_cms_lm_all_engines(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite46_suite46_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("46_")));
    let row = ql_benches::benchmark_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite47_suite47_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("47_")));
    let row =
        ql_benches::benchmark_hagan_west_ctd_afv_kth_sacva_gs_kirk_fx_vv_napoleon_tc_meucci_krr_eva(
            1,
        );
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite48_suite48_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("48_")));
    let row = ql_benches::benchmark_brigo_csa_ctd_fabozzi_callable_oas_lando_jlt_migration_gatheral_corridor_varswap_glasserman_lsm_rogers_dual_cartea_lob_mm_lo_mackinlay_vr_lopez_dsr_purged_cv_roncalli_spinu_erc_ilmanen_vmc_penman_reoi_mccall_bellman(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite49_suite49_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("49_")));
    let row = ql_benches::benchmark_brigo_nss_krd_tuckman_tips_floor_schonbucher_rmv_rt_rfv_gatheral_sticky_lv_gregory_simm_ccp_mva_schwartz_smith_kalman_wystup_fx_dnt_bouzoubaa_timer_meucci_bl_woodbury_gueant_ac_sinh_zou_elastic_net_nk_ree(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite50_suite50_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("50_")));
    let row = ql_benches::benchmark_hagan_west_tension_duffie_lando_pac_cmo_fama_macbeth_shanken_gatheral_asvi_bgk_discrete_barrier_bouchaud_cross_impact_karnosky_singer_merton_jump_campbell_viceira_ohlson_lid_griewank_revolve(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite51_suite51_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("51_")));
    let row = ql_benches::benchmark_brigo_xccy_mtm_fabozzi_spec_pool_io_po_lando_marshall_olkin_gatheral_vix_vvix_gregory_wwr_cva_geman_gibson_schwartz_clark_fx_fva_rubinstein_chooser_grinold_flam_tc_meucci_es_euler_james_stein_sure_campbell_cochrane_habit(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite52_suite52_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("52_")));
    let row = ql_benches::benchmark_brigo_ois_fallback_fabozzi_collared_frn_schonbucher_cds_swaption_fep_bergomi_2f_skew_gregory_sa_ccr_ead_eydeland_spark_spread_wystup_prdc_geske_compound_qian_alpha_decay_roncalli_enb_johansen_vecm_bansal_yaron_lrr(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite53_suite53_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("53_")));
    let row = ql_benches::benchmark_andersen_cms_replication_tuckman_special_repo_okane_cds_upfront_gatheral_ssr_tanaka_gregory_saccr_kva_eydeland_crack_spread_wystup_prdc_cap_haug_eso_dilution_qian_spinu_erc_rockafellar_cvar_htf_pls_krylov_piazzesi_fomc(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite54_suite54_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("54_")));
    let row = ql_benches::benchmark_darbyshire_step_ois_turn_despiegeleer_cb_soft_call_okane_cds_index_skew_gatheral_varswap_truncation_gregory_fva_overlap_geman_working_storage_wystup_malz_fx_fly_goldman_sosin_lookback_ledoit_wolf_shrinkage_almgren_chriss_cara_tsay_estar_acm_cp(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite55_suite55_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("55_")));
    let row = ql_benches::benchmark_henrard_tenor_basis_fabozzi_msr_recapture_okane_cdo_base_corr_gatheral_fwd_start_smile_gregory_pykhtin_csa_mta_eydeland_gas_swing_wystup_fx_kiko_stulz_rainbow_jagannathan_ma_obizhaeva_wang_newey_west_hac_gabaix_disaster(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite56_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("56_")));
    let row = ql_benches::benchmark_brigo_diff_swap_quanto_fabozzi_callable_muni_ytw_bielecki_step_up_bond_bkm_moment_strip_gregory_break_clause_cva_bessembinder_power_skew_wystup_fx_fader_longstaff_extendible_kan_zhou_3fund_sharpe_rbsa_diebold_mariano_hln_laubach_williams_rstar(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite57_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("57_")));
    let row = ql_benches::benchmark_andersen_sofr_futures_fra_hayre_mbs_scurve_pedersen_cdx_fep_gatheral_dupire_total_var_gregory_simm_mva_routledge_stockout_baba_cip_basis_carr_chou_symmetry_pastor_stambaugh_treynor_black_engle_granger_ecm_epstein_zin_sdf(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite58_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("58_")));
    let row = ql_benches::benchmark_litterman_scheinkman_pca_fly_tf_split_cb_hull_white_double_t_cdo_roger_lee_moment_wings_brigo_capponi_cds_contagion_schwartz_ou_futures_wystup_dcd_yield_passport_option_grs_factor_ftest_bouchaud_sqrt_impact_bns_bipower_har_cj_lettau_ludvigson_cay(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite59_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("59_")));
    let row = ql_benches::benchmark_piterbarg_collateral_option_tuckman_rolldown_merton_kmv_dd_hagan_bartlett_sabr_pykhtin_mpor_spike_deng_locational_ftr_wystup_digital_overhedge_conze_fixed_lookback_tasche_hot_spots_es_frazzini_bab_pesaran_timmermann_lucas_mehra_prescott(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite60_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("60_")));
    let row = ql_benches::benchmark_gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_ga_carr_wu_vol_drift_gregory_marginal_cva_working_brennan_storage_wystup_forward_extra_rubinstein_exchange_opt_hlz_factor_fdr_holden_subrahmanyam_kyle_bollerslev_garch_kurtosis_campbell_shiller_cape(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite61_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("61_")));
    let row = ql_benches::benchmark_kim_wright_shadow_rate_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_swap_burgard_kjaer_fva_carmona_durrleman_spark_wystup_fx_tarf_bgk_discrete_lookback_novy_marx_ff5_acerbi_spectral_es_corsi_har_rv_bll_news_noise(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite62_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("62_")));
    let row = ql_benches::benchmark_jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_regimes_hull_white_wwr_cva_gibson_schwartz_samuelson_reiner_quanto_fx_turnbull_wakeman_asian_amp_value_momentum_garleanu_pedersen_aim_bai_ng_icp2_bgg_accelerator(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite63_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("63_")));
    let row = ql_benches::benchmark_lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_cds_fouque_multiscale_sv_gregory_bilateral_ftd_eydeland_vpp_heat_rate_wystup_fx_corridor_rubinstein_reiner_barrier_lns_gls_r2_bellini_expectile_bns_quarticity_kiyotaki_moore_cycle(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite64_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("64_")));
    let row = ql_benches::benchmark_jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_albanese_sa_ccr_kva_rss_equilibrium_storage_du_tepper_verdelhan_cip_wilmott_cliquet_bsv_parametric_portfolio_ahmadi_javid_evar_engle_kroner_bekk_gabaix_disaster(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite65_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("65_")));
    let row = ql_benches::benchmark_bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_smile_brigo_pallavicini_ccp_xva_dmt_cointegrated_crack_piterbarg_prdc_tarn_bgk_discrete_barrier_connor_korajczyk_apt_afs_transient_impact_hamilton_markov_switching_bansal_yaron_lrr(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite66_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("66_")));
    let row = ql_benches::benchmark_andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_vol_swap_green_leverage_ratio_kva_cortazar_schwartz_3f_oil_fujii_collateralized_xccs_pechtl_haug_time_switch_boudt_cornish_fisher_mvar_rockafellar_quadrangle_hasbrouck_info_share_he_krishnamurthy_sdf(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite67_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("67_")));
    let row = ql_benches::benchmark_piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_green_kenyon_tva_eydeland_interconnector_spread_wystup_clark_fx_strangle_kunitomo_ikeda_curved_barrier_clmx_variance_decomp_kyle_multi_period_lambda_hansen_lunde_overnight_rv_goldstein_pauzner_bank_run(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite68_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("68_")));
    let row = ql_benches::benchmark_mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_burgard_kjaer_semi_replication_xva_benth_seasonal_hdd_cdd_swap_baba_packer_fx_swap_dollar_basis_conze_viswanathan_lookback_he_litterman_bl_posterior_ghs_importance_sampling_cvar_engle_granger_johansen_vecm_gabaix_rare_disasters_resilience(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite69_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("69_")));
    let row = ql_benches::benchmark_ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_gregory_simm_mva_carmona_clean_dark_spread_reiner_haug_quanto_call_hart_ross_soft_barrier_pastor_stambaugh_alpha_shrinkage_obizhaeva_wang_impulse_lob_diebold_mariano_west_hln_piazzesi_schneider_nominal_tp(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite70_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("70_")));
    let row = ql_benches::benchmark_smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_pykhtin_rosen_marginal_cva_hilliard_reis_jump_futures_gabaix_maggiori_fx_carry_cheuk_vorst_shout_option_grossman_zhou_drawdown_kelly_gls_tasche_cornish_fisher_raroc_andrews_ploberger_wald_break_jermann_quadrini_financial_shock(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite71_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("71_")));
    let row = ql_benches::benchmark_mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_kenyon_green_csa_lva_boogert_de_jong_gas_ratchet_garman_kohlhagen_baw_fx_hui_haug_binary_barrier_kan_robotti_hj_distance_cartea_jaimungal_order_flow_ait_sahalia_jacod_jump_test_wachter_disaster_volatility(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite72_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("72_")));
    let row = ql_benches::benchmark_hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_gregory_brigo_ftd_bilateral_cva_brennan_schwartz_mine_hysteresis_wystup_fx_fader_broadie_detemple_capped_call_meucci_entropy_pooling_cont_deguest_scandolo_robust_risk_pesaran_shin_smith_ardl_christiano_eichenbaum_evans_nk(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite73_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("73_")));
    let row = ql_benches::benchmark_mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_andersen_pykhtin_sokol_mpor_bjerksund_stensland_321_crack_carr_chou_symmetry_barrier_chesney_jeanblanc_yor_parisian_demiguel_brodie_norm_shrinkage_gueant_lehalle_tapia_mm_hansen_white_spa_test_brunnermeier_sannikov_paradox(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite74_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("74_")));
    let row = ql_benches::benchmark_filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_hull_white_burgard_kjaer_fva_windfall_casassus_collin_dufresne_convenience_yield_lustig_roussanov_verdelhan_fx_carry_carr_lee_corridor_varswap_garlappi_uppal_wang_multi_prior_mcneil_frey_evt_garch_patton_sjc_copula_pastor_veronesi_political_uncertainty(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite75_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("75_")));
    let row = ql_benches::benchmark_lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_brigo_capponi_pallavicini_contagion_wwr_routledge_seppi_spatt_storage_piterbarg_prdc_hybrid_broadie_glasserman_kou_discrete_barrier_brandt_santa_clara_valkanov_ppp_gatheral_schied_slynko_transient_impact_corsi_bpq_harq_lettau_ludvigson_cay(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite76_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("76_")));
    let row = ql_benches::benchmark_cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_crepey_bichuch_capponi_xva_bsde_trolle_schwartz_usv_commodity_ivashina_scharfstein_stein_cip_wilmott_windcliff_cliquet_skew_ledoit_wolf_nonlinear_shrinkage_adrian_brunnermeier_delta_covar_bollerslev_tauchen_zhou_vrp_campbell_vuolteenaho_bad_good_beta(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite77_workload_is_registered() {
    assert!(ql_benches::ALL_SUITES.iter().any(|s| s.starts_with("77_")));
    let row = ql_benches::benchmark_bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_andersen_duffie_song_debt_overhang_fva_pirrong_structural_inventory_spread_lipton_castagna_fx_dnt_vanna_volga_linetsky_step_barrier_koijen_yogo_demand_system_kyle_obizhaeva_invariance_gu_kelly_xiu_oos_r2_drechsler_savov_schnabl_deposit_channel(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}
