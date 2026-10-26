use ql_benches::harness::with_config_override;
use ql_benches::{
    benchmark_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin,
    benchmark_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad, run_all_benchmarks,
    suites::{
        benchmark_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp,
        benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad,
        benchmark_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml, suite35_workloads,
    },
};

#[test]
fn suite33_suite33_and_invariants() {
    let rows = with_config_override(1, 3, || run_all_benchmarks(1));
    assert_eq!(rows.len(), 55);
    let row = with_config_override(1, 3, || {
        benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad(1)
    });
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
    assert!(row.parity_ok);
}

#[test]
fn suite35_suite35_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml(1)
    });
    assert_eq!(
        row.id,
        "35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml"
    );
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
    let values = suite35_workloads([0.2_f64, -0.1, 0.3, -0.2]);
    assert_eq!(values.len(), 12);
    assert!(values.iter().all(|v| v.is_finite()));
}

#[test]
fn suite34_suite34_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad(1)
    });
    assert_eq!(row.id, "34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad");
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
}

#[test]
fn suite36_suite36_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp(1)
    });
    assert!(row.parity_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn suite37_suite37_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin(
            1,
        )
    });
    assert!(row.parity_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn suite38_suite38_and_invariants() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc(1)
    });
    assert!(row.parity_ok);
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
    assert!(row.speedup > 1.0);
    let x_dual = [
        ql_ad::Dual::<4>::variable(0.2, 0),
        ql_ad::Dual::<4>::variable(-0.1, 1),
        ql_ad::Dual::<4>::variable(0.3, 2),
        ql_ad::Dual::<4>::variable(-0.2, 3),
    ];
    let y_dual = ql_benches::suites::suite43_workloads(x_dual);
    for (idx, item) in y_dual.iter().enumerate() {
        let grad_l1: f64 = item.d.iter().map(|v| v.abs()).sum();
        assert!(grad_l1 > 1e-9, "entry {idx} has zero Dual<4> gradient");
    }
}

#[test]
fn suite39_suite39_and_invariants() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt(1)
    });
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite40_suite40_and_invariants() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca(1)
    });
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite41_suite41_and_invariants() {
    let row = ql_benches::benchmark_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is(1);
    assert!(row.max_abs_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn suite42_suite42_and_invariants() {
    let row = ql_benches::benchmark_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin(1);
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn suite43_suite43_and_invariants() {
    let row = ql_benches::benchmark_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi(1);
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn suite44_suite44_and_invariants() {
    use ql_ad::Dual;
    let x = [
        Dual::<4>::variable(0.2, 0),
        Dual::variable(-0.1, 1),
        Dual::variable(0.3, 2),
        Dual::variable(-0.2, 3),
    ];
    for y in ql_benches::suite44_workloads(x) {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    let row = ql_benches::benchmark_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite45_suite45_and_invariants() {
    use ql_ad::Dual;
    let x: [Dual<4>; 4] = [
        Dual::variable(0.2, 0),
        Dual::variable(-0.1, 1),
        Dual::variable(0.3, 2),
        Dual::variable(-0.2, 3),
    ];
    for y in ql_benches::suite45_workloads(x) {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    assert_eq!(ql_benches::SUITE_COUNT, 55);
    let row = ql_benches::benchmark_rebonato_cms_lm_all_engines(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite46_suite46_and_invariants() {
    use ql_ad::Dual;
    let x: [Dual<4>; 4] = [
        Dual::variable(0.2, 0),
        Dual::variable(-0.1, 1),
        Dual::variable(0.3, 2),
        Dual::variable(-0.2, 3),
    ];
    for y in ql_benches::suite46_workloads(x) {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    assert_eq!(ql_benches::SUITE_COUNT, 55);
    let row = ql_benches::benchmark_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite47_suite47_and_invariants() {
    use ql_ad::Dual;
    let x: [Dual<4>; 4] = [
        Dual::variable(0.2, 0),
        Dual::variable(-0.1, 1),
        Dual::variable(0.3, 2),
        Dual::variable(-0.2, 3),
    ];
    for y in ql_benches::suite47_workloads(x) {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    assert_eq!(ql_benches::SUITE_COUNT, 55);
    let row =
        ql_benches::benchmark_hagan_west_ctd_afv_kth_sacva_gs_kirk_fx_vv_napoleon_tc_meucci_krr_eva(
            1,
        );
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite48_suite48_and_invariants() {
    use ql_ad::Dual;
    let x: [Dual<4>; 4] = [
        Dual::variable(0.2, 0),
        Dual::variable(-0.1, 1),
        Dual::variable(0.3, 2),
        Dual::variable(-0.2, 3),
    ];
    for y in ql_benches::suite48_workloads(x) {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    assert_eq!(ql_benches::SUITE_COUNT, 55);
    let row = ql_benches::benchmark_brigo_csa_ctd_fabozzi_callable_oas_lando_jlt_migration_gatheral_corridor_varswap_glasserman_lsm_rogers_dual_cartea_lob_mm_lo_mackinlay_vr_lopez_dsr_purged_cv_roncalli_spinu_erc_ilmanen_vmc_penman_reoi_mccall_bellman(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite49_suite49_and_invariants() {
    use ql_ad::Dual;
    let x: [Dual<4>; 4] = [
        Dual::variable(0.2, 0),
        Dual::variable(-0.1, 1),
        Dual::variable(0.3, 2),
        Dual::variable(-0.2, 3),
    ];
    for y in ql_benches::suite49_workloads(x) {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    assert_eq!(ql_benches::SUITE_COUNT, 55);
    let row = ql_benches::benchmark_brigo_nss_krd_tuckman_tips_floor_schonbucher_rmv_rt_rfv_gatheral_sticky_lv_gregory_simm_ccp_mva_schwartz_smith_kalman_wystup_fx_dnt_bouzoubaa_timer_meucci_bl_woodbury_gueant_ac_sinh_zou_elastic_net_nk_ree(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite50_suite50_and_invariants() {
    use ql_ad::Dual;
    let x: [Dual<4>; 4] = [
        Dual::variable(0.2, 0),
        Dual::variable(-0.1, 1),
        Dual::variable(0.3, 2),
        Dual::variable(-0.2, 3),
    ];
    for y in ql_benches::suite50_workloads(x) {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    assert_eq!(ql_benches::SUITE_COUNT, 55);
    let row = ql_benches::benchmark_hagan_west_tension_duffie_lando_pac_cmo_fama_macbeth_shanken_gatheral_asvi_bgk_discrete_barrier_bouchaud_cross_impact_karnosky_singer_merton_jump_campbell_viceira_ohlson_lid_griewank_revolve(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite51_suite51_and_invariants() {
    use ql_ad::Dual;
    let x: [Dual<4>; 4] = [
        Dual::variable(0.2, 0),
        Dual::variable(-0.1, 1),
        Dual::variable(0.3, 2),
        Dual::variable(-0.2, 3),
    ];
    for y in ql_benches::suite51_workloads(x) {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    assert_eq!(ql_benches::SUITE_COUNT, 55);
    let row = ql_benches::benchmark_brigo_xccy_mtm_fabozzi_spec_pool_io_po_lando_marshall_olkin_gatheral_vix_vvix_gregory_wwr_cva_geman_gibson_schwartz_clark_fx_fva_rubinstein_chooser_grinold_flam_tc_meucci_es_euler_james_stein_sure_campbell_cochrane_habit(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite52_suite52_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite52_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite52_workloads(x2);
    for (a, b) in y_dual.iter().zip(y_dual2.iter()) {
        assert!((a.val - b.val).abs() < 1e-12);
        for k in 0..4 {
            assert!((a.d[k] - b.d1[k]).abs() < 1e-11);
        }
    }
    let guard = TapeScope::guard();
    let xt = TapeScope::with_active(|t| {
        [
            t.variable(x0[0]),
            t.variable(x0[1]),
            t.variable(x0[2]),
            t.variable(x0[3]),
        ]
    });
    let yt = ql_benches::suite52_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 55);
    let row = ql_benches::benchmark_brigo_ois_fallback_fabozzi_collared_frn_schonbucher_cds_swaption_fep_bergomi_2f_skew_gregory_sa_ccr_ead_eydeland_spark_spread_wystup_prdc_geske_compound_qian_alpha_decay_roncalli_enb_johansen_vecm_bansal_yaron_lrr(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite53_suite53_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite53_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite53_workloads(x2);
    for (a, b) in y_dual.iter().zip(y_dual2.iter()) {
        assert!((a.val - b.val).abs() < 1e-12);
        for k in 0..4 {
            assert!((a.d[k] - b.d1[k]).abs() < 1e-11);
        }
    }
    let guard = TapeScope::guard();
    let xt = TapeScope::with_active(|t| {
        [
            t.variable(x0[0]),
            t.variable(x0[1]),
            t.variable(x0[2]),
            t.variable(x0[3]),
        ]
    });
    let yt = ql_benches::suite53_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 55);
    let row = ql_benches::benchmark_andersen_cms_replication_tuckman_special_repo_okane_cds_upfront_gatheral_ssr_tanaka_gregory_saccr_kva_eydeland_crack_spread_wystup_prdc_cap_haug_eso_dilution_qian_spinu_erc_rockafellar_cvar_htf_pls_krylov_piazzesi_fomc(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite54_suite54_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite54_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite54_workloads(x2);
    for (a, b) in y_dual.iter().zip(y_dual2.iter()) {
        assert!((a.val - b.val).abs() < 1e-12);
        for k in 0..4 {
            assert!((a.d[k] - b.d1[k]).abs() < 1e-11);
        }
    }
    let guard = TapeScope::guard();
    let xt = TapeScope::with_active(|t| {
        [
            t.variable(x0[0]),
            t.variable(x0[1]),
            t.variable(x0[2]),
            t.variable(x0[3]),
        ]
    });
    let yt = ql_benches::suite54_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 55);
    let row = ql_benches::benchmark_darbyshire_step_ois_turn_despiegeleer_cb_soft_call_okane_cds_index_skew_gatheral_varswap_truncation_gregory_fva_overlap_geman_working_storage_wystup_malz_fx_fly_goldman_sosin_lookback_ledoit_wolf_shrinkage_almgren_chriss_cara_tsay_estar_acm_cp(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn suite55_suite55_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite55_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite55_workloads(x2);
    for (a, b) in y_dual.iter().zip(y_dual2.iter()) {
        assert!((a.val - b.val).abs() < 1e-12);
        for k in 0..4 {
            assert!((a.d[k] - b.d1[k]).abs() < 1e-11);
        }
    }
    let guard = TapeScope::guard();
    let xt = TapeScope::with_active(|t| {
        [
            t.variable(x0[0]),
            t.variable(x0[1]),
            t.variable(x0[2]),
            t.variable(x0[3]),
        ]
    });
    let yt = ql_benches::suite55_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 55);
    let row = ql_benches::benchmark_henrard_tenor_basis_fabozzi_msr_recapture_okane_cdo_base_corr_gatheral_fwd_start_smile_gregory_pykhtin_csa_mta_eydeland_gas_swing_wystup_fx_kiko_stulz_rainbow_jagannathan_ma_obizhaeva_wang_newey_west_hac_gabaix_disaster(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}
