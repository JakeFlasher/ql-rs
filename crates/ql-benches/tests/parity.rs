use ql_benches::harness::with_config_override;
use ql_benches::{
    benchmark_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad,
    benchmark_phase25_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin,
    run_all_benchmarks,
    suites::{
        benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad,
        benchmark_phase24_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp,
        benchmark_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml, phase23_workloads,
    },
};

#[test]
fn phase21_suite33_and_invariants() {
    let rows = with_config_override(1, 3, || run_all_benchmarks(1));
    assert_eq!(rows.len(), 46);
    let row = with_config_override(1, 3, || {
        benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad(1)
    });
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
    assert!(row.parity_ok);
}

#[test]
fn phase23_suite35_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml(1)
    });
    assert_eq!(
        row.id,
        "35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml"
    );
    assert!(row.parity_ok && row.max_abs_diff < 1e-4);
    let values = phase23_workloads([0.2_f64, -0.1, 0.3, -0.2]);
    assert_eq!(values.len(), 12);
    assert!(values.iter().all(|v| v.is_finite()));
}

#[test]
fn phase22_suite34_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad(1)
    });
    assert_eq!(row.id, "34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad");
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
}

#[test]
fn phase24_suite36_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_phase24_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp(1)
    });
    assert!(row.parity_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn phase25_suite37_and_invariants() {
    let row = with_config_override(1, 3, || {
        benchmark_phase25_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin(1)
    });
    assert!(row.parity_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn phase26_suite38_and_invariants() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_phase26_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc(1)
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
    let y_dual = ql_benches::suites::phase31_workloads(x_dual);
    for (idx, item) in y_dual.iter().enumerate() {
        let grad_l1: f64 = item.d.iter().map(|v| v.abs()).sum();
        assert!(grad_l1 > 1e-9, "entry {idx} has zero Dual<4> gradient");
    }
}

#[test]
fn phase27_suite39_and_invariants() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_phase27_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt(1)
    });
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn phase28_suite40_and_invariants() {
    let row = with_config_override(1, 3, || {
        ql_benches::benchmark_phase28_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca(1)
    });
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn phase29_suite41_and_invariants() {
    let row = ql_benches::benchmark_phase29_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is(1);
    assert!(row.max_abs_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn phase30_suite42_and_invariants() {
    let row = ql_benches::benchmark_phase30_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin(1);
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn phase31_suite43_and_invariants() {
    let row = ql_benches::benchmark_phase31_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi(1);
    assert!(row.parity_ok);
    assert!(row.max_abs_diff < 1e-4);
    assert!(row.speedup > 1.0);
}

#[test]
fn phase32_suite44_and_invariants() {
    use ql_ad::Dual;
    let x = [
        Dual::<4>::variable(0.2, 0),
        Dual::variable(-0.1, 1),
        Dual::variable(0.3, 2),
        Dual::variable(-0.2, 3),
    ];
    for y in ql_benches::phase32_workloads(x) {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    let row = ql_benches::benchmark_phase32_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn phase33_suite45_and_invariants() {
    use ql_ad::Dual;
    let x: [Dual<4>; 4] = [
        Dual::variable(0.2, 0),
        Dual::variable(-0.1, 1),
        Dual::variable(0.3, 2),
        Dual::variable(-0.2, 3),
    ];
    for y in ql_benches::phase33_workloads(x) {
        assert!(y.d.iter().map(|v| v.abs()).sum::<f64>() > 1e-9);
    }
    assert_eq!(ql_benches::SUITE_COUNT, 46);
    let row = ql_benches::benchmark_phase33_rebonato_cms_lm_all_engines(1);
    assert!(row.parity_ok && row.max_abs_diff < 1e-4 && row.speedup > 1.0);
}

#[test]
fn phase34_suite46_and_invariants() {
    let x = [0.2_f64, -0.1, 0.3, -0.2];
    let y = ql_benches::phase34_workloads(x);
    assert_eq!(y.len(), 12);
    assert!(y.iter().all(|v| v.is_finite()));
}
