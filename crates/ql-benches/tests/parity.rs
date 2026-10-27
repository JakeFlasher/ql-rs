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
    assert_eq!(rows.len(), 85);
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
}

#[test]
fn suite37_suite37_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin(
            1,
        )
    });
    assert!(row.parity_diff < 1e-4);
}

#[test]
fn suite38_suite38_and_invariants() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc(1)
    });
    assert!(row.parity_ok);
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
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
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite40_suite40_and_invariants() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca(1)
    });
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite41_suite41_and_invariants() {
    let row = ql_benches::benchmark_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is(1);
    assert!(row.max_abs_diff < 1e-4);
}

#[test]
fn suite42_suite42_and_invariants() {
    let row = ql_benches::benchmark_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin(1);
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
}

#[test]
fn suite43_suite43_and_invariants() {
    let row = ql_benches::benchmark_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi(1);
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
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
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
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
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_rebonato_cms_lm_all_engines(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
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
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
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
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row =
        ql_benches::benchmark_hagan_west_ctd_afv_kth_sacva_gs_kirk_fx_vv_napoleon_tc_meucci_krr_eva(
            1,
        );
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
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
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_brigo_csa_ctd_fabozzi_callable_oas_lando_jlt_migration_gatheral_corridor_varswap_glasserman_lsm_rogers_dual_cartea_lob_mm_lo_mackinlay_vr_lopez_dsr_purged_cv_roncalli_spinu_erc_ilmanen_vmc_penman_reoi_mccall_bellman(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
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
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_brigo_nss_krd_tuckman_tips_floor_schonbucher_rmv_rt_rfv_gatheral_sticky_lv_gregory_simm_ccp_mva_schwartz_smith_kalman_wystup_fx_dnt_bouzoubaa_timer_meucci_bl_woodbury_gueant_ac_sinh_zou_elastic_net_nk_ree(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
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
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_hagan_west_tension_duffie_lando_pac_cmo_fama_macbeth_shanken_gatheral_asvi_bgk_discrete_barrier_bouchaud_cross_impact_karnosky_singer_merton_jump_campbell_viceira_ohlson_lid_griewank_revolve(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
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
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_brigo_xccy_mtm_fabozzi_spec_pool_io_po_lando_marshall_olkin_gatheral_vix_vvix_gregory_wwr_cva_geman_gibson_schwartz_clark_fx_fva_rubinstein_chooser_grinold_flam_tc_meucci_es_euler_james_stein_sure_campbell_cochrane_habit(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
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
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_brigo_ois_fallback_fabozzi_collared_frn_schonbucher_cds_swaption_fep_bergomi_2f_skew_gregory_sa_ccr_ead_eydeland_spark_spread_wystup_prdc_geske_compound_qian_alpha_decay_roncalli_enb_johansen_vecm_bansal_yaron_lrr(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
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
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_andersen_cms_replication_tuckman_special_repo_okane_cds_upfront_gatheral_ssr_tanaka_gregory_saccr_kva_eydeland_crack_spread_wystup_prdc_cap_haug_eso_dilution_qian_spinu_erc_rockafellar_cvar_htf_pls_krylov_piazzesi_fomc(1);
    // Timing ratios are host-dependent; the numerical parity gate is the
    // deterministic invariant checked by this test.
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
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
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_darbyshire_step_ois_turn_despiegeleer_cb_soft_call_okane_cds_index_skew_gatheral_varswap_truncation_gregory_fva_overlap_geman_working_storage_wystup_malz_fx_fly_goldman_sosin_lookback_ledoit_wolf_shrinkage_almgren_chriss_cara_tsay_estar_acm_cp(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
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
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_henrard_tenor_basis_fabozzi_msr_recapture_okane_cdo_base_corr_gatheral_fwd_start_smile_gregory_pykhtin_csa_mta_eydeland_gas_swing_wystup_fx_kiko_stulz_rainbow_jagannathan_ma_obizhaeva_wang_newey_west_hac_gabaix_disaster(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite56_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite56_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite56_workloads(x2);
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
    let yt = ql_benches::suite56_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_brigo_diff_swap_quanto_fabozzi_callable_muni_ytw_bielecki_step_up_bond_bkm_moment_strip_gregory_break_clause_cva_bessembinder_power_skew_wystup_fx_fader_longstaff_extendible_kan_zhou_3fund_sharpe_rbsa_diebold_mariano_hln_laubach_williams_rstar(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
}

#[test]
fn suite57_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite57_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite57_workloads(x2);
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
    let yt = ql_benches::suite57_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_andersen_sofr_futures_fra_hayre_mbs_scurve_pedersen_cdx_fep_gatheral_dupire_total_var_gregory_simm_mva_routledge_stockout_baba_cip_basis_carr_chou_symmetry_pastor_stambaugh_treynor_black_engle_granger_ecm_epstein_zin_sdf(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite58_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite58_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite58_workloads(x2);
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
    let yt = ql_benches::suite58_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_litterman_scheinkman_pca_fly_tf_split_cb_hull_white_double_t_cdo_roger_lee_moment_wings_brigo_capponi_cds_contagion_schwartz_ou_futures_wystup_dcd_yield_passport_option_grs_factor_ftest_bouchaud_sqrt_impact_bns_bipower_har_cj_lettau_ludvigson_cay(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite59_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite59_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite59_workloads(x2);
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
    let yt = ql_benches::suite59_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_piterbarg_collateral_option_tuckman_rolldown_merton_kmv_dd_hagan_bartlett_sabr_pykhtin_mpor_spike_deng_locational_ftr_wystup_digital_overhedge_conze_fixed_lookback_tasche_hot_spots_es_frazzini_bab_pesaran_timmermann_lucas_mehra_prescott(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite60_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite60_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite60_workloads(x2);
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
    let yt = ql_benches::suite60_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_ga_carr_wu_vol_drift_gregory_marginal_cva_working_brennan_storage_wystup_forward_extra_rubinstein_exchange_opt_hlz_factor_fdr_holden_subrahmanyam_kyle_bollerslev_garch_kurtosis_campbell_shiller_cape(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite61_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite61_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite61_workloads(x2);
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
    let yt = ql_benches::suite61_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_kim_wright_shadow_rate_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_swap_burgard_kjaer_fva_carmona_durrleman_spark_wystup_fx_tarf_bgk_discrete_lookback_novy_marx_ff5_acerbi_spectral_es_corsi_har_rv_bll_news_noise(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite62_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite62_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite62_workloads(x2);
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
    let yt = ql_benches::suite62_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_regimes_hull_white_wwr_cva_gibson_schwartz_samuelson_reiner_quanto_fx_turnbull_wakeman_asian_amp_value_momentum_garleanu_pedersen_aim_bai_ng_icp2_bgg_accelerator(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite63_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite63_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite63_workloads(x2);
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
    let yt = ql_benches::suite63_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_cds_fouque_multiscale_sv_gregory_bilateral_ftd_eydeland_vpp_heat_rate_wystup_fx_corridor_rubinstein_reiner_barrier_lns_gls_r2_bellini_expectile_bns_quarticity_kiyotaki_moore_cycle(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite64_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite64_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite64_workloads(x2);
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
    let yt = ql_benches::suite64_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_albanese_sa_ccr_kva_rss_equilibrium_storage_du_tepper_verdelhan_cip_wilmott_cliquet_bsv_parametric_portfolio_ahmadi_javid_evar_engle_kroner_bekk_gabaix_disaster(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite65_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite65_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite65_workloads(x2);
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
    let yt = ql_benches::suite65_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_smile_brigo_pallavicini_ccp_xva_dmt_cointegrated_crack_piterbarg_prdc_tarn_bgk_discrete_barrier_connor_korajczyk_apt_afs_transient_impact_hamilton_markov_switching_bansal_yaron_lrr(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite66_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite66_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite66_workloads(x2);
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
    let yt = ql_benches::suite66_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_vol_swap_green_leverage_ratio_kva_cortazar_schwartz_3f_oil_fujii_collateralized_xccs_pechtl_haug_time_switch_boudt_cornish_fisher_mvar_rockafellar_quadrangle_hasbrouck_info_share_he_krishnamurthy_sdf(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite67_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite67_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite67_workloads(x2);
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
    let yt = ql_benches::suite67_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_green_kenyon_tva_eydeland_interconnector_spread_wystup_clark_fx_strangle_kunitomo_ikeda_curved_barrier_clmx_variance_decomp_kyle_multi_period_lambda_hansen_lunde_overnight_rv_goldstein_pauzner_bank_run(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite68_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite68_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite68_workloads(x2);
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
    let yt = ql_benches::suite68_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_burgard_kjaer_semi_replication_xva_benth_seasonal_hdd_cdd_swap_baba_packer_fx_swap_dollar_basis_conze_viswanathan_lookback_he_litterman_bl_posterior_ghs_importance_sampling_cvar_engle_granger_johansen_vecm_gabaix_rare_disasters_resilience(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite69_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite69_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite69_workloads(x2);
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
    let yt = ql_benches::suite69_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_gregory_simm_mva_carmona_clean_dark_spread_reiner_haug_quanto_call_hart_ross_soft_barrier_pastor_stambaugh_alpha_shrinkage_obizhaeva_wang_impulse_lob_diebold_mariano_west_hln_piazzesi_schneider_nominal_tp(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite70_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite70_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite70_workloads(x2);
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
    let yt = ql_benches::suite70_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_pykhtin_rosen_marginal_cva_hilliard_reis_jump_futures_gabaix_maggiori_fx_carry_cheuk_vorst_shout_option_grossman_zhou_drawdown_kelly_gls_tasche_cornish_fisher_raroc_andrews_ploberger_wald_break_jermann_quadrini_financial_shock(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite71_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite71_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite71_workloads(x2);
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
    let yt = ql_benches::suite71_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_kenyon_green_csa_lva_boogert_de_jong_gas_ratchet_garman_kohlhagen_baw_fx_hui_haug_binary_barrier_kan_robotti_hj_distance_cartea_jaimungal_order_flow_ait_sahalia_jacod_jump_test_wachter_disaster_volatility(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite72_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite72_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite72_workloads(x2);
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
    let yt = ql_benches::suite72_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_gregory_brigo_ftd_bilateral_cva_brennan_schwartz_mine_hysteresis_wystup_fx_fader_broadie_detemple_capped_call_meucci_entropy_pooling_cont_deguest_scandolo_robust_risk_pesaran_shin_smith_ardl_christiano_eichenbaum_evans_nk(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite73_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite73_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite73_workloads(x2);
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
    let yt = ql_benches::suite73_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_andersen_pykhtin_sokol_mpor_bjerksund_stensland_321_crack_carr_chou_symmetry_barrier_chesney_jeanblanc_yor_parisian_demiguel_brodie_norm_shrinkage_gueant_lehalle_tapia_mm_hansen_white_spa_test_brunnermeier_sannikov_paradox(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite74_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite74_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite74_workloads(x2);
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
    let yt = ql_benches::suite74_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_hull_white_burgard_kjaer_fva_windfall_casassus_collin_dufresne_convenience_yield_lustig_roussanov_verdelhan_fx_carry_carr_lee_corridor_varswap_garlappi_uppal_wang_multi_prior_mcneil_frey_evt_garch_patton_sjc_copula_pastor_veronesi_political_uncertainty(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite75_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite75_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite75_workloads(x2);
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
    let yt = ql_benches::suite75_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_brigo_capponi_pallavicini_contagion_wwr_routledge_seppi_spatt_storage_piterbarg_prdc_hybrid_broadie_glasserman_kou_discrete_barrier_brandt_santa_clara_valkanov_ppp_gatheral_schied_slynko_transient_impact_corsi_bpq_harq_lettau_ludvigson_cay(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite76_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite76_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite76_workloads(x2);
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
    let yt = ql_benches::suite76_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_crepey_bichuch_capponi_xva_bsde_trolle_schwartz_usv_commodity_ivashina_scharfstein_stein_cip_wilmott_windcliff_cliquet_skew_ledoit_wolf_nonlinear_shrinkage_adrian_brunnermeier_delta_covar_bollerslev_tauchen_zhou_vrp_campbell_vuolteenaho_bad_good_beta(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite77_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite77_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite77_workloads(x2);
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
    let yt = ql_benches::suite77_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_andersen_duffie_song_debt_overhang_fva_pirrong_structural_inventory_spread_lipton_castagna_fx_dnt_vanna_volga_linetsky_step_barrier_koijen_yogo_demand_system_kyle_obizhaeva_invariance_gu_kelly_xiu_oos_r2_drechsler_savov_schnabl_deposit_channel(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite78_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite78_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite78_workloads(x2);
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
    let yt = ql_benches::suite78_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_albanese_caenazzo_crepey_bilateral_kva_acharya_lochstoer_ramadorai_hedging_pressure_brunnermeier_nagel_pedersen_fx_crash_broadie_yamamoto_fgt_asian_kelly_pruitt_su_ipca_brownlees_engle_srisk_andersen_bollerslev_diebold_labys_arfima_farhi_gabaix_disaster_fx(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite79_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite79_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite79_workloads(x2);
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
    let yt = ql_benches::suite79_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_duffie_zhu_ccp_netting_gorton_hayashi_rouwenhorst_inventory_basis_della_corte_ramadorai_sarno_fx_vrp_pelsser_laplace_double_barrier_moreira_muir_volatility_managed_budish_cramton_shim_fba_ait_sahalia_mykland_zhang_tsrv_gertler_karadi_credit_policy(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite80_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite80_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite80_workloads(x2);
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
    let yt = ql_benches::suite80_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_bauer_rudebusch_qe_signaling_boyarchenko_fuster_lucca_mbs_oas_azizpour_giesecke_schwenkler_default_clustering_bollerslev_todorov_jump_fear_ghamami_glasserman_clearing_mva_szymanowska_de_roon_nijman_commodity_premia_menkhoff_sarno_schmeling_schrimpf_fx_vol_milevsky_posner_reciprocal_gamma_asian_kozak_nagel_santosh_sdf_shrinkage_acharya_pedersen_philippon_richardson_ses_bnhls_realized_kernel_brunnermeier_pedersen_margin_spiral(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite81_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite81_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite81_workloads(x2);
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
    let yt = ql_benches::suite81_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_gurkaynak_sack_swanson_fomc_path_bao_pan_wang_bond_illiquidity_collin_dufresne_goldstein_martin_credit_determinants_andersen_fusari_todorov_left_tail_capponi_cheng_ccp_margin_kang_rouwenhorst_tang_two_premiums_lustig_roussanov_verdelhan_dollar_carry_geman_yor_bessel_asian_tu_zhou_markowitz_talmud_biais_foucault_moinas_fast_trading_christensen_oomen_podolskij_uhf_jumps_gabaix_granular_residual(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite82_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite82_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite82_workloads(x2);
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
    let yt = ql_benches::suite82_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_nakamura_steinsson_fomc_info_dick_nielsen_feldhutter_lando_illiquidity_berndt_duffie_credit_premia_christoffersen_heston_jacobs_variance_kernel_brigo_francischello_pallavicini_closeout_xva_boons_prado_basis_momentum_colacito_croce_gavazzoni_ready_recursive_fx_fusai_meucci_erlang_asian_freyberger_neuhierl_weber_characteristics_kelly_jiang_tail_risk_jacod_li_mykland_podolskij_vetter_preaveraging_adrian_etula_muir_broker_dealer_sdf(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite83_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite83_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite83_workloads(x2);
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
    let yt = ql_benches::suite83_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_wright_inflation_uncertainty_damico_fan_kitsul_tips_liquidity_culp_nozawa_veronesi_pseudo_firm_bollerslev_todorov_xu_tail_predictability_piterbarg_lou_multicurrency_csa_ctd_basak_pavlova_commodity_financialization_mueller_tahbaz_salehi_vedolin_fomc_fx_kou_wang_jump_first_passage_barrier_ao_li_zheng_maxser_portfolio_menkveld_multi_venue_hft_patton_sheppard_realized_semivariance_bksy_macro_volatility_lrr(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite84_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite84_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite84_workloads(x2);
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
    let yt = ql_benches::suite84_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_duffie_singleton_swap_spread_chen_lesmond_wei_lot_illiquidity_giesecke_longstaff_schaefer_strebulaev_default_regime_drechsler_yaron_uncertainty_jump_vrp_glasserman_yang_entropy_wwr_cva_hamilton_wu_affine_commodity_premia_hassan_country_size_currency_safety_linetsky_whittaker_spectral_asian_barillas_shanken_factor_squared_sharpe_giglio_kelly_pruitt_systemic_pqr_mancini_threshold_realized_variance_constantinides_duffie_heterogeneous_sdf(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}

#[test]
fn suite85_and_invariants() {
    use ql_ad::{Dual, Dual2, TapeScope};
    let x0 = [0.2_f64, -0.1, 0.3, -0.2];
    let x: [Dual<4>; 4] = [
        Dual::variable(x0[0], 0),
        Dual::variable(x0[1], 1),
        Dual::variable(x0[2], 2),
        Dual::variable(x0[3], 3),
    ];
    let y_dual = ql_benches::suite85_workloads(x);
    for y in &y_dual {
        assert!(y.d.iter().all(|v| v.abs() > 1e-9));
    }
    let x2: [Dual2<4>; 4] = [
        Dual2::variable(x0[0], 0),
        Dual2::variable(x0[1], 1),
        Dual2::variable(x0[2], 2),
        Dual2::variable(x0[3], 3),
    ];
    let y_dual2 = ql_benches::suite85_workloads(x2);
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
    let yt = ql_benches::suite85_workloads(xt);
    for (idx, out) in yt.iter().enumerate() {
        let grad = TapeScope::gradient_from(*out);
        for k in 0..4 {
            assert!((grad.gradient(xt[k]) - y_dual[idx].d[k]).abs() < 1e-11);
        }
    }
    drop(guard);
    assert_eq!(ql_benches::SUITE_COUNT, 85);
    let row = ql_benches::benchmark_bikbov_chernov_macro_yield_krishnamurthy_vissing_jorgensen_treasury_convenience_bhamra_kuehn_strebulaev_credit_equity_pan_implicit_jump_fear_vrp_bo_capponi_contagion_bilateral_cva_singleton_commodity_investor_flows_verdelhan_habit_fx_premium_davydov_linetsky_cev_barrier_demiguel_garlappi_uppal_naive_diversification_engle_manganelli_asymmetric_caviar_barndorff_nielsen_multipower_variation_clt_kogan_papanikolaou_ist_shock_sdf(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-3);
}
