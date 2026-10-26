use ql_benches::harness::with_config_override;
use ql_benches::run_all_benchmark_suites_with_config;

#[test]
fn benchmarks_verify_parity_and_optional_sota_speedups() {
    // A short deterministic run is enough to exercise all 20 numerical
    // parity paths during normal `cargo test`; timing gates are inherently
    // host-dependent and opt in through QL_RUN_TIMING_GATES=1.
    let rows = run_all_benchmark_suites_with_config(1, 3);
    assert_eq!(rows.len(), 64);
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
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
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
