use ql_terminal::{
    BrailleCanvas, ChartMode, DeskWorkstation, ScenarioKind, TabId, render_diverging_bar,
    render_sparkline,
};

#[test]
fn braille_canvas_and_sparkline_rasterize_accurately() {
    let mut canvas = BrailleCanvas::new(8, 4);
    canvas.draw_line(
        0,
        0,
        canvas.width as isize - 1,
        canvas.height as isize - 1,
        [0, 220, 255],
    );
    let lines = canvas.render_lines(false);
    assert!(
        lines
            .iter()
            .any(|line| line.chars().any(|c| ('\u{2801}'..='\u{28ff}').contains(&c)))
    );
    let spark = render_sparkline(&[-1.0, 0.0, 1.0], 12);
    assert_eq!(spark.chars().count(), 12);
    assert!(spark.contains('▁') && spark.contains('█'));
    let bar = render_diverging_bar(-2.0, 4.0, 8, false);
    assert_eq!(bar.chars().count(), 17);
    assert_eq!(bar.chars().nth(8), Some('│'));
    let positive = render_diverging_bar(2.0, 4.0, 8, false);
    assert_eq!(positive.chars().count(), 17);
    assert_eq!(positive.chars().nth(8), Some('│'));
    let colored = render_diverging_bar(2.0, 4.0, 8, true);
    assert!(colored.contains("\x1b[38;2;0;230;140m"));
}

#[test]
fn workstation_multi_asset_ticks_and_stress_scenarios_are_consistent() {
    let mut ws = DeskWorkstation::new();
    let first = ws.step_tick();
    let micro = ws.step_tick();
    assert!(!micro.dag_curves_rebuilt);
    assert!(micro.dag_recomputed_trades < micro.dag_total_trades);
    let normal_plot = ws.render_dashboard(&first, false);
    ws.set_scenario(ScenarioKind::RatesPlus50Bp);
    let rates = ws.step_tick();
    assert!(rates.dag_curves_rebuilt);
    assert!(
        rates.total_npv != first.total_npv
            || rates.ir_dv01_buckets.iter().any(|(_, x)| x.abs() > 0.0)
    );
    let rates_plot = ws.render_dashboard(&rates, false);
    assert_ne!(normal_plot, rates_plot);
    ws.set_scenario(ScenarioKind::EquityCrashMinus15Pct);
    let equity = ws.step_tick();
    assert_ne!(equity.total_npv, first.total_npv);
    ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
    let credit = ws.step_tick();
    assert!(
        credit
            .credit_cs01_buckets
            .iter()
            .any(|(_, x)| x.abs() > 0.0)
    );
    ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
    let stag = ws.step_tick();
    assert!(stag.inflation_par_delta_sum.abs() > 0.0);
}

#[test]
fn dashboard_renders_all_quadrants_and_chart_modes() {
    let mut ws = DeskWorkstation::new();
    let frame = ws.step_tick();
    for mode in [
        ChartMode::YieldAndInflationCurves,
        ChartMode::RoughHestonAndFxSmiles,
        ChartMode::XvaExposureProfile,
    ] {
        ws.set_chart_mode(mode);
        let plain = ws.render_dashboard(&frame, false);
        assert!(plain.contains(mode.label()));
        assert!(plain.contains("1. BOOK SUMMARY"));
        assert!(plain.contains("2. BRAILLE CANVAS"));
        assert!(plain.contains("3. IFT BUCKETED"));
        assert!(plain.contains("4. CROSS-ASSET"));
        assert!(plain.contains("IRS 5Y PayFix"));
        assert!(plain.contains("Heston Call 1Y"));
        assert!(plain.contains("RoughHest H=.1"));
        assert!(plain.contains("Autocall Note"));
        assert!(plain.contains("CDS 5Y BuyProt"));
        assert!(plain.contains("ZCIIS 5Y CPI"));
        assert!(plain.contains("FX 25D VV Call"));
        assert!(plain.contains("LEVY-VG-COS-1Y"));
        assert!(plain.contains("CONV-TF-BOND-3Y"));
        assert!(
            plain
                .chars()
                .any(|c| ('\u{2801}'..='\u{28ff}').contains(&c))
        );
        let ansi = ws.render_dashboard(&frame, true);
        assert!(ansi.contains("\x1b[38;2;"));
    }
}

#[test]
fn suite23_tabs_and_live_telemetry_are_present() {
    assert_eq!(TabId::ALL.len(), 80);
    assert_eq!(TabId::ALL[18], TabId::BermudanVolIft);
    assert_eq!(TabId::ALL[19], TabId::SaCcrSwingQp);
    assert_eq!(TabId::ALL[22], TabId::MonotoneConvexVixParisianRobust);
    assert_eq!(TabId::ALL[24], TabId::JyInflationTbaBergomiNco);
    let mut ws = DeskWorkstation::new();
    let frame = ws.step_tick();
    assert!(frame.vix_parisian_robust_metric.is_finite());
    let dashboard = ws.render_dashboard(&frame, false);
    assert!(dashboard.contains("Bermudan HW") && dashboard.contains("SA-CCR EAD"));
    assert!(dashboard.contains("Commodity Swing") && dashboard.contains("Box Markowitz"));
    assert!(dashboard.contains("23 Monotone/VIX/Parisian/Robust"));
    assert!(dashboard.contains("25 JY-Inflation/TBA/Bergomi2F/NCO"));
    ws.set_scenario(ScenarioKind::RatesPlus50Bp);
    let stressed = ws.step_tick();
    assert_ne!(
        stressed.vix_parisian_robust_metric,
        frame.vix_parisian_robust_metric
    );
}

#[test]
fn suite29_tab_metric_and_rates_scenario_are_live() {
    assert_eq!(TabId::ALL.len(), 80);
    assert!(TabId::ALL.contains(&TabId::AswMf1fAndreasenHugeBessEp));
    let mut ws = DeskWorkstation::new();
    let frame = ws.step_tick();
    assert!(frame.asw_mf1f_ah_bess_ep_metric.is_finite());
    assert!(
        ws.render_dashboard(&frame, false)
            .contains("27 ASW/MF1F/AH-LV/BESS/Idzorek-DP")
    );
    let normal_next = ws.step_tick();
    assert_ne!(
        normal_next.asw_mf1f_ah_bess_ep_metric,
        frame.asw_mf1f_ah_bess_ep_metric
    );
    for scenario in [
        ScenarioKind::RatesPlus50Bp,
        ScenarioKind::EquityCrashMinus15Pct,
        ScenarioKind::CreditWideningPlus100Bp,
        ScenarioKind::StagflationPlus150Bp,
    ] {
        ws.set_scenario(scenario);
        let stressed = ws.step_tick();
        assert_ne!(
            stressed.asw_mf1f_ah_bess_ep_metric,
            normal_next.asw_mf1f_ah_bess_ep_metric
        );
    }
}

#[test]
fn suite30_tab28_terminal_renders_and_ticks_cleanly() {
    let mut ws = DeskWorkstation::new();
    let normal = ws.step_tick();
    assert!(
        ws.render_dashboard(&normal, false)
            .contains("BL / Risk Parity ERC")
    );
    assert!(normal.mtmxccs_colva_zabr_weather_enb_metric.is_finite());
    for scenario in [
        ScenarioKind::RatesPlus50Bp,
        ScenarioKind::EquityCrashMinus15Pct,
        ScenarioKind::CreditWideningPlus100Bp,
        ScenarioKind::StagflationPlus150Bp,
    ] {
        ws.set_scenario(scenario);
        let stressed = ws.step_tick();
        assert_ne!(
            stressed.mtmxccs_colva_zabr_weather_enb_metric,
            normal.mtmxccs_colva_zabr_weather_enb_metric
        );
    }
}

#[test]
fn suite31_tab29_terminal_renders_and_ticks_cleanly() {
    assert_eq!(TabId::ALL.len(), 80);
    assert_eq!(
        TabId::BurnoutCmcdsQuinticAfvCsSchur.label(),
        "29 Burnout-MBS / CMCDS / Quintic-OU / AFV-CB / CS-Swaption / Schur-HRP"
    );
    let mut ws = DeskWorkstation::new();
    let normal = ws.step_tick();
    assert!(normal.burnout_cmcds_quintic_afv_cs_schur_metric.is_finite());
    let dashboard = ws.render_dashboard(&normal, false);
    assert!(dashboard.contains("29 Burnout/CMCDS/Quintic/AFV/Schur"));
    assert!(dashboard.contains("BL / Risk Parity ERC"));
    let next = ws.step_tick();
    assert_ne!(
        next.burnout_cmcds_quintic_afv_cs_schur_metric,
        normal.burnout_cmcds_quintic_afv_cs_schur_metric
    );
    for scenario in [
        ScenarioKind::RatesPlus50Bp,
        ScenarioKind::EquityCrashMinus15Pct,
        ScenarioKind::CreditWideningPlus100Bp,
        ScenarioKind::StagflationPlus150Bp,
    ] {
        ws.set_scenario(scenario);
        let stressed = ws.step_tick();
        assert_ne!(
            stressed.burnout_cmcds_quintic_afv_cs_schur_metric,
            normal.burnout_cmcds_quintic_afv_cs_schur_metric
        );
    }
}

#[test]
fn suite32_tab30_terminal_renders_and_ticks_cleanly() {
    assert_eq!(TabId::ALL.len(), 80);
    assert_eq!(
        TabId::SinHouMbRBergomiChoiAfsExpectile.label(),
        "30 Sinh-OU / BK-FRN / Morini-Brigo / WWM-XVA / mrBergomi / Choi / AFS-LOB / EXVaR"
    );
    let mut ws = DeskWorkstation::new();
    let normal = ws.step_tick();
    assert!(
        ws.render_dashboard(&normal, false)
            .contains("30 SinhOU/MB/mrBergomi/Choi/EXVaR")
    );
    let next = ws.step_tick();
    assert_ne!(
        next.sinhou_mb_rbergomi_choi_afs_expectile_metric,
        normal.sinhou_mb_rbergomi_choi_afs_expectile_metric
    );
    for scenario in [
        ScenarioKind::RatesPlus50Bp,
        ScenarioKind::EquityCrashMinus15Pct,
        ScenarioKind::CreditWideningPlus100Bp,
        ScenarioKind::StagflationPlus150Bp,
    ] {
        ws.set_scenario(scenario);
        let stressed = ws.step_tick();
        assert_ne!(
            stressed.sinhou_mb_rbergomi_choi_afs_expectile_metric,
            normal.sinhou_mb_rbergomi_choi_afs_expectile_metric
        );
    }
}

#[test]
fn suite33_tab31_terminal_renders_and_ticks_cleanly() {
    assert_eq!(TabId::ALL.len(), 80);
    assert_eq!(
        TabId::MonotoneConvexClarkSaccrRvolCartea.label(),
        "31 HW-MC / Clark-FX / SA-CCR-KVA / BFG-rVol / Cartea-Alpha / Paleologo"
    );
    let mut ws = DeskWorkstation::new();
    let normal = ws.step_tick();
    assert!(
        normal
            .monotone_convex_clark_saccr_rvol_cartea_paleologo_metric
            .is_finite()
    );
    assert!(
        ws.render_dashboard(&normal, false)
            .contains("31 HW-MC/Clark/SACCR/rVol/Cartea/Paleo")
    );
    let next = ws.step_tick();
    assert_ne!(
        next.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric,
        normal.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric
    );
    for scenario in [
        ScenarioKind::RatesPlus50Bp,
        ScenarioKind::EquityCrashMinus15Pct,
        ScenarioKind::CreditWideningPlus100Bp,
        ScenarioKind::StagflationPlus150Bp,
    ] {
        ws.set_scenario(scenario);
        let stressed = ws.step_tick();
        assert!(
            stressed
                .monotone_convex_clark_saccr_rvol_cartea_paleologo_metric
                .is_finite()
        );
        assert_ne!(
            stressed.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric,
            normal.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric
        );
    }
}

#[test]
fn suite34_tab32_terminal_renders_and_ticks_cleanly() {
    assert_eq!(TabId::ALL.len(), 80);
    assert_eq!(
        TabId::KrdLpiCocoBergomi2fCgmyGueantMeucci.label(),
        "32 KRD / LPI / CoCo / Bergomi / CGMY / AAD"
    );
    let mut ws = DeskWorkstation::new();
    let normal = ws.step_tick();
    assert!(
        normal
            .krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
            .is_finite()
    );
    assert!(
        ws.render_dashboard(&normal, false)
            .contains("32 KRD/LPI/CoCo/Bergomi2F/CGMY/Guéant/Meucci")
    );
    let next = ws.step_tick();
    assert!(
        next.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
            .is_finite()
    );
    assert_ne!(
        next.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric,
        normal.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
    );
    for scenario in [
        ScenarioKind::RatesPlus50Bp,
        ScenarioKind::EquityCrashMinus15Pct,
        ScenarioKind::CreditWideningPlus100Bp,
        ScenarioKind::StagflationPlus150Bp,
    ] {
        ws.set_scenario(scenario);
        let stressed = ws.step_tick();
        assert_ne!(
            stressed.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric,
            normal.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
        );
    }
}

#[test]
fn suite35_tab33_terminal_renders_and_reacts() {
    assert_eq!(TabId::ALL.len(), 80);
    assert!(TabId::ALL.contains(&TabId::RkhsFmmVvLsvJltCtdQrhParisianAidNco));
    assert_eq!(
        TabId::RkhsFmmVvLsvJltCtdQrhParisianAidNco.label(),
        "33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML"
    );
    let mut ws = DeskWorkstation::new();
    let base = ws.step_tick();
    assert!(
        base.rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric
            .is_finite()
    );
    assert!(
        ws.render_dashboard(&base, false)
            .contains("33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML")
    );
    let tick = ws.step_tick();
    assert_ne!(
        base.rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric,
        tick.rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric
    );
    for scenario in [
        ScenarioKind::RatesPlus50Bp,
        ScenarioKind::EquityCrashMinus15Pct,
        ScenarioKind::CreditWideningPlus100Bp,
        ScenarioKind::StagflationPlus150Bp,
    ] {
        ws.set_scenario(scenario);
        let stressed = ws.step_tick();
        assert!(
            stressed
                .rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric
                .is_finite()
        );
        assert_ne!(
            base.rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric,
            stressed.rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric
        );
    }
}

#[test]
fn suite36_tab34_terminal_renders_and_reacts() {
    assert_eq!(
        TabId::MbsIoPoRmvMoXvaDilutionSigGarchPot.label(),
        "34 MBS-IO-PO / RMV / MO-XVA / Dilution / Sig / GARCH / POT / RIE / PIN / GP / Brinson / IERP"
    );
    let mut ws = DeskWorkstation::new();
    let mut values = Vec::new();
    for scenario in [
        ScenarioKind::RatesParallelPlus50Bp,
        ScenarioKind::VolShockPlus5Pts,
        ScenarioKind::CreditSpreadWidening100Bp,
        ScenarioKind::EquityCrashMinus15Pct,
        ScenarioKind::StagflationPlus150Bp,
    ] {
        ws.set_scenario(scenario);
        let f = ws.step_tick();
        values.push(f.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric);
        let frame = ws.render_dashboard(&f, true);
        assert!(frame.contains("34 MBS/RMV/MO/CB/Sig/HN/POT/RIE/PIN/GP/BC/IERP"));
    }
    assert!(values.iter().all(|v| v.is_finite()));
    assert!(values.windows(2).any(|w| (w[1] - w[0]).abs() > 1e-12));
}

#[test]
fn suite37_tab35_terminal_renders_and_reacts() {
    let mut ws = DeskWorkstation::new();
    assert_eq!(TabId::ALL.len(), 80);
    assert!(TabId::KalmanRtsVecmHjSdfIfrs9CptHasbrouck.label().starts_with("35 Kalman-RTS / VECM / HJ-SDF / IFRS9 / CPT / Hasbrouck / FV / Lando / LPPL / CB-Arb / Causal / Malliavin"));
    let mut values = Vec::new();
    for scenario in [
        ScenarioKind::RatesParallelPlus50Bp,
        ScenarioKind::VolShockPlus5Pts,
        ScenarioKind::CreditSpreadWidening100Bp,
        ScenarioKind::EquityCrashMinus15Pct,
        ScenarioKind::StagflationPlus150Bp,
    ] {
        ws.set_scenario(scenario);
        let f = ws.step_tick();
        values
            .push(f.kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric);
        assert!(
            f.kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric
                .is_finite()
        );
        assert!(ws.render_dashboard(&f, false).contains("35 Kalman/VECM"));
    }
    assert!(values.windows(2).any(|w| (w[1] - w[0]).abs() > 1e-12));
}

#[test]
fn suite38_tab36_terminal_renders_and_reacts() {
    let mut ws = DeskWorkstation::new();
    let normal = ws.step_tick();
    ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
    let stress = ws.step_tick();
    assert!(
        stress
            .mf_rates_wildcard_ctd_pivot_tarf_bergomi_metric
            .is_finite()
    );
    assert!(
        (stress.mf_rates_wildcard_ctd_pivot_tarf_bergomi_metric
            - normal.mf_rates_wildcard_ctd_pivot_tarf_bergomi_metric)
            .abs()
            > 1e-12
    );
    assert!(ws.render_dashboard(&stress, false).contains(
        "36 MF-Rates/Wildcard/TARF/Bergomi/Schönbucher/FRTB-CVA/Lookback/VTA/FLAM/ACD/HLZ/MLMC"
    ));
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(TabId::ALL.contains(&TabId::MfRatesWildcardCtdPivotTarfBergomi));
}

#[test]
fn suite39_tab37_terminal_renders_and_reacts() {
    let mut ws = DeskWorkstation::new();
    let normal = ws.step_tick();
    ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
    let stressed = ws.step_tick();
    assert!(
        stressed
            .cochrane_ftpl_g2_twist_essvi_offshore_cb_metric
            .is_finite()
    );
    assert!(
        (stressed.cochrane_ftpl_g2_twist_essvi_offshore_cb_metric
            - normal.cochrane_ftpl_g2_twist_essvi_offshore_cb_metric)
            .abs()
            > 1e-12
    );
    assert!(
        ws.render_dashboard(&stressed, false)
            .contains("37 FTPL/G2/eSSVI/CB/CoCo/MPOR/Double-Heston/Curran/Warin/Roll/Busseti/KKT")
    );
    assert_eq!(TabId::TAB_COUNT, 80);
}

#[test]
fn suite40_tab38_terminal_renders_and_reacts() {
    let mut ws = DeskWorkstation::default();
    let fa = ws.step_tick();
    let a = ws.render_dashboard(&fa, false);
    ws.scenario = ScenarioKind::StagflationPlus150Bp;
    let fb = ws.step_tick();
    let b = ws.render_dashboard(&fb, false);
    assert_eq!(TabId::ALL.len(), 80);
    assert!(
        a.contains("38 IRRBB/RFR/XCCY/BG/CoCo/SIMM/rBergomi/Spread/BESS/Dark/Risk/Diff-PCA")
            && b.contains("38 IRRBB/RFR/XCCY/BG/CoCo/SIMM/rBergomi/Spread/BESS/Dark/Risk/Diff-PCA")
    );
    assert_ne!(a, b);
}

#[test]
fn suite41_tab39_terminal_renders_and_reacts() {
    let mut ws = DeskWorkstation::default();
    let frame = ws.step_tick();
    let dashboard = ws.render_dashboard(&frame, false);
    assert_eq!(TabId::ALL.len(), 80);
    assert!(
        dashboard
            .contains("39 Veronesi/Lando/Corridor/Russian/Bates/Koziol/MFG/BCBS/CJP/ETF/IC/IS")
    );
    assert!(
        frame
            .veronesi_conv_lando_gen_gatheral_corridor_metric
            .is_finite()
    );
    let base = ws.veronesi_conv_lando_gen_gatheral_corridor_live_metric();
    ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
    let credit = ws.veronesi_conv_lando_gen_gatheral_corridor_live_metric();
    assert!((credit - base).abs() > 1e-12);
    ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
    let stag = ws.veronesi_conv_lando_gen_gatheral_corridor_live_metric();
    assert!((stag - base).abs() > 1e-12);
}

#[test]
fn suite42_tab40_terminal_renders_and_reacts() {
    let mut ws = DeskWorkstation::new();
    ws.set_scenario(ScenarioKind::Normal);
    let normal = ws.step_tick();
    ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
    let stressed = ws.step_tick();
    assert_eq!(TabId::TAB_COUNT, 80);
    assert_eq!(TabId::ALL.len(), 80);
    assert!(
        TabId::RomerTaylorDamodaranCrpChristoffersenVar
            .label()
            .starts_with("40 Taylor-Curve / Damodaran-CRP / Cumulant-VaR / LOB-LVaR / Prospect-CGO / Regret-Match / Beneish-HKAS32 / Campisi-FI / Heckman-GMM / Elastic-Net / Reset-CB / FSB-BailIn")
    );
    assert_ne!(
        normal.romer_taylor_damodaran_crp_christoffersen_var_metric,
        stressed.romer_taylor_damodaran_crp_christoffersen_var_metric
    );
    ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
    let stag = ws.step_tick();
    assert!(
        (stag.romer_taylor_damodaran_crp_christoffersen_var_metric
            - normal.romer_taylor_damodaran_crp_christoffersen_var_metric)
            .abs()
            > 1e-12
    );
    assert!(ws.render_dashboard(&stressed, false).contains(
        "40 Taylor/Damodaran/VaR/LVaR/CGO/Regret/Beneish/Campisi/Heckman/Elastic/ResetCB/FSB"
    ));
}

#[test]
fn suite43_tab41_terminal_renders_and_reacts() {
    let mut ws = DeskWorkstation::new();
    let normal = ws.step_tick();
    ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
    let stressed = ws.step_tick();
    assert_eq!(TabId::TAB_COUNT, 80);
    assert_eq!(TabId::ALL.len(), 80);
    assert!(TabId::ALL.contains(&TabId::FabozziKrdCmoPacAsquithCbAsc470Colva));
    assert!(
        TabId::FabozziKrdCmoPacAsquithCbAsc470Colva
            .label()
            .starts_with("41 KRD-LSC / CMO-PAC / Asquith-CB / ASU2024-CPDI / ASC470-Mod / ColVA-CTD / Hill-EVT / Bouchaud-TLF / Karnosky-FX / DGTW-Timing / Brock-Hommes / Penman-ReOI")
    );
    assert!(
        normal
            .fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric
            .is_finite()
    );
    assert!(
        (normal.fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric
            - stressed.fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric)
            .abs()
            > 1e-12
    );
    ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
    let stag = ws.step_tick();
    assert!(
        (stag.fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric
            - normal.fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric)
            .abs()
            > 1e-12
    );
    assert!(
        ws.render_dashboard(&stressed, false)
            .contains("41 KRD-LSC/CMO-PAC/Asquith-CB/ASU2024/ASC470/ColVA/Hill/Bouchaud/Karnosky/DGTW/Brock/Penman")
    );
}

#[test]
fn suite44_tab42_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert_eq!(TabId::ALL.len(), 80);
    assert!(
        TabId::CheyetteCapletAdbNdirsZabrQrhCmcdsAdi
            .label()
            .starts_with("42 Cheyette-LV / ADB-NDIRS / Andreasen-ZABR / QRH-Zumbach / CMCDS-Morini / Turlakov-WWR / HHW-ADI / CappedCall-MW / Ahmadi-EVaR / Schur-HRP / BEKK-GARCH / Schilit-Sloan")
    );
    let mut desk = DeskWorkstation::default();
    let base = desk.cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let credit = desk.cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let stag = desk.cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_live_metric();
    assert_ne!(base, credit);
    assert_ne!(base, stag);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("42 Cheyette/NDIRS/ZABR/QRH/CMCDS/WWR/HHW/CB/EVaR/Schur/BEKK/Quality")
    );
}

#[test]
fn suite45_tab43_terminal_renders_and_reacts() {
    assert_eq!(ql_terminal::TabId::TAB_COUNT, 80);
    assert_eq!(ql_terminal::TabId::ALL.len(), 80);
    assert!(
        ql_terminal::TabId::RebonatoCmsLelandToftAscotPoolRollRie
            .label()
            .starts_with("43 Rebonato-CMS / Leland-Toft / Calamos-Ascot / SpecPool-Roll / Glosten-Milgrom / Potters-RIE / EGARCH-Hamilton / Wooldridge-2SLS / McNeil-Euler / Bacon-Omega / Damodaran-RD / Arrow-Debreu")
    );
    let mut ws = ql_terminal::DeskWorkstation::new();
    let n = ws.rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric();
    let _ = ws.step_tick();
    ws.set_scenario(ql_terminal::ScenarioKind::CreditWideningPlus100Bp);
    let c = ws.rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric();
    ws.set_scenario(ql_terminal::ScenarioKind::StagflationPlus150Bp);
    let st = ws.rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric();
    assert!(n.is_finite() && (n - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
    let frame = ws.step_tick();
    assert!(
        ws.render_dashboard(&frame, false)
            .contains("43 Rebonato-CMS/Leland-Toft/Calamos/RIE/Arrow-Debreu")
    );
}

#[test]
fn suite46_tab44_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::TuckmanRepoBkFvaDecsIoPoKylePinGjrLppl
            .label()
            .starts_with("44 Tuckman-Repo")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("44 Tuckman-Repo")
    );
}

#[test]
fn suite47_tab45_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::HaganWestCtdSwitchAfvCbKthDefaultSaCva
            .label()
            .starts_with("45 Hagan-West-MC")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("45 Hagan-West-MC")
    );
}

#[test]
fn suite48_tab46_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::BrigoCsaCtdCallableOasJltCorridorLsm
            .label()
            .starts_with("46 Multi-CSA-CTD")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.brigo_csa_ctd_callable_oas_jlt_corridor_lsm_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.brigo_csa_ctd_callable_oas_jlt_corridor_lsm_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.brigo_csa_ctd_callable_oas_jlt_corridor_lsm_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("46 Multi-CSA-CTD")
    );
}

#[test]
fn suite49_tab47_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::BrigoNssKrdTipsFloorRmvTriangleStickyLv
            .label()
            .starts_with("47 NSS-KRD")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(desk.render_dashboard(&frame, false).contains("47 NSS-KRD"));
}

#[test]
fn suite50_tab48_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::HaganWestTensionPacCmoFamaMacBethAsviBgk
            .label()
            .starts_with("48 Tension-Curve")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("48 Tension-Curve")
    );
}

#[test]
fn suite51_tab49_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::BrigoXccyMtmSpecPoolMarshallOlkinVixVvix
            .label()
            .starts_with("49 XCCY-MtM")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(desk.render_dashboard(&frame, false).contains("49 XCCY-MtM"));
}

#[test]
fn suite52_tab50_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::BrigoOisFallbackCollaredFrnCdsFepSparkPrdc
            .label()
            .starts_with("50 OIS-Fallback")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("50 OIS-Fallback")
    );
}

#[test]
fn suite53_tab51_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::AndersenCmsSpecialRepoCdsUpfrontSsrKvaCrack
            .label()
            .starts_with("51 CMS-Replication")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.andersen_cms_special_repo_cds_upfront_ssr_kva_crack_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.andersen_cms_special_repo_cds_upfront_ssr_kva_crack_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.andersen_cms_special_repo_cds_upfront_ssr_kva_crack_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(desk.render_dashboard(&frame, false).contains("51 CMS-Rep"));
}

#[test]
fn suite54_tab52_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::DarbyshireStepOisCbSoftCallCdxSkewMalzFly
            .label()
            .starts_with("52 Step-OIS-Turn")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("52 Step-OIS-Turn")
    );
}

#[test]
fn suite55_tab53_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::HenrardTenorBasisMsrRecaptureCdoBaseCorrKiko
            .label()
            .starts_with("53 Tenor-Basis-Mult")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("53 Tenor-Basis/MSR-Recap")
    );
}

#[test]
fn suite56_tab54_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::BrigoDiffSwapCallableMuniBkmSkewBreakCva
            .label()
            .starts_with("54 Diff-Swap-Quanto")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.brigo_diff_swap_callable_muni_bkm_skew_break_cva_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.brigo_diff_swap_callable_muni_bkm_skew_break_cva_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.brigo_diff_swap_callable_muni_bkm_skew_break_cva_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("54 Diff-Swap/Muni-YTW")
    );
}

#[test]
fn suite57_tab55_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::AndersenSofrFutHayreMbsPedersenCdxDupireLv
            .label()
            .starts_with("55 SOFR-Fut-FRA")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("55 SOFR-Fut/Hayre-MBS")
    );
}

#[test]
fn suite58_tab56_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::LittermanScheinkmanPcaFlyTfSplitCbDoubleTCdoLeeWings
            .label()
            .starts_with("56 PCA-Fly-Carry")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("56 PCA-Fly/TF-CB")
    );
}

#[test]
fn suite59_tab57_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::PiterbargCollateralOptTuckmanRolldownMertonKmvBartlettSabr
            .label()
            .starts_with("57 Piterbarg-Coll-Opt")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("57 Piterbarg-Opt/Tuckman-Roll")
    );
}

#[test]
fn suite60_tab58_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::GswSvenssonCurveFabozziMbsEsfVasicekGordyIrbCarrWuVol
            .label()
            .starts_with("58 GSW-Svensson-Curve")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("58 GSW-Svensson/MBS-ESF")
    );
}

#[test]
fn suite61_tab59_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::KimWrightShadowBrennanSchwartzCbDuffieSingletonRmvFukasawaGamma
            .label()
            .starts_with("59 Kim-Wright-Shadow")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk
        .kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk
        .kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk
        .kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("59 Kim-Wright-Shadow/Brennan-Schwartz-CB")
    );
}

#[test]
fn suite62_tab60_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::JpsUnspannedMacroKalotayRefundingHullWhiteNtdDermanSticky
            .label()
            .starts_with("60 JPS-Unspanned-Macro")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("60 JPS-Unspanned-Macro/Kalotay-Refunding")
    );
}

#[test]
fn suite63_tab61_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::LyashenkoFmmCapletTuckmanCpiSeasonalityBrigoAlfonsiCirFouqueSv
            .label()
            .starts_with("61 Lyashenko-FMM-Caplet")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n =
        desk.lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c =
        desk.lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s =
        desk.lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("61 Lyashenko-FMM-Caplet/Tuckman-CPI-Seasonality")
    );
}

#[test]
fn suite64_tab62_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::JszCanonicalGdtsmTsiveriotisFernandesCbLaurentGregoryCdoBroadieJainVarswap
            .label()
            .starts_with("62 JSZ-Canonical-GDTSM")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("62 JSZ-Canonical-GDTSM/Tsiveriotis-Fernandes-CB")
    );
}

#[test]
fn suite65_tab63_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::BianchettiTwoCurveQuantoBdtCallableBondSchonbucherCdsSwaptionBatesSvj
            .label()
            .starts_with("63 Bianchetti-Two-Curve-Quanto")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("63 Bianchetti-Two-Curve-Quanto/BDT-Callable")
    );
}

#[test]
fn suite66_tab64_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::AndersenPiterbargCmsSpreadFabozziZBondCmoDuffieGarleanuAffineCdoCarrLee
            .label()
            .starts_with("64 Andersen-Piterbarg-CMS-Spread")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("64 Andersen-Piterbarg-CMS-Spread/Fabozzi-Z-Bond")
    );
}

#[test]
fn suite67_tab65_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::PiterbargCtdCollateralTuckmanMakeWholeBondBrigoMoriniCmcdsHaganShiftedSabr
            .label()
            .starts_with("65 Piterbarg-CTD-Collateral")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("65 Piterbarg-CTD-Collateral/Tuckman-Make-Whole")
    );
}

#[test]
fn suite68_tab66_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::MercurioFvaMultiCurveSwapIngersollContingentCallCbHullWhiteImpliedCopulaCdoGatheralJacquierSviJw
            .label()
            .starts_with("66 Mercurio-FVA-Swap")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("66 Mercurio-FVA-Swap/Ingersoll-CB")
    );
}

#[test]
fn suite69_tab67_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::AmetranoBianchettiTenorBasisFabozziHayreIoPoStripHullBrigoWwrCdsBossensVannaVolgaBarrier
            .label()
            .starts_with("67 Ametrano-Bianchetti-Tenor-Basis")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("67 Ametrano-Bianchetti-Tenor-Basis/Fabozzi-Hayre-IO-PO")
    );
}

#[test]
fn suite70_tab68_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::SmithWilsonUfrKalotayWilliamsFabozziCallableSchonbucherClaytonCdoBuehlerVarianceCurve
            .label()
            .starts_with("68 Smith-Wilson-UFR")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("68 Smith-Wilson-UFR/Kalotay-Williams-Fabozzi")
    );
}

#[test]
fn suite71_tab69_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::McCullochTaxCurveAndersenBuffumCbOkaneRecoveryLockCdsDermanKaniImpliedTree
            .label()
            .starts_with("69 McCulloch-Tax-Curve")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("69 McCulloch-Tax-Curve/Andersen-Buffum-CB")
    );
}

#[test]
fn suite72_tab70_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::HaganWestMonotoneConvexFabozziEffectiveDurationOasFryeAltmanDownturnLgdGatheralJacquierSviArbitrage
            .label()
            .starts_with("70 Hagan-West-Monotone-Convex")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("70 Hagan-West-Monotone-Convex/Fabozzi-Effective-Duration-OAS")
    );
}

#[test]
fn suite73_tab71_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::MercurioRussoIborOisBasisJarrowYildirimYoyInflationDuffieLandoHybridCdsBayerFrizGatheralRoughSkew
            .label()
            .starts_with("71 Mercurio-Russo-IBOR-OIS-Basis")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("71 Mercurio-Russo-IBOR-OIS-Basis/Jarrow-Yildirim-YoY-Inflation")
    );
}

#[test]
fn suite74_tab72_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::FilipovicLarssonTrolleLinearRationalDuffieKrishnamurthySpecialRepoDuffieSaitaWangStochasticDefaultCarrMadanNeubergerLogContract
            .label()
            .starts_with("72 Filipovic-Larsson-Trolle-LRTSM")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("72 Filipovic-Larsson-Trolle-LRTSM/Duffie-Krishnamurthy-Repo")
    );
}

#[test]
fn suite75_tab73_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::LyashenkoMercurioRfrInArrearsPennacchiSundaresanCocoAt1PedersenBrigoCdsIndexFepBroadieJainSeppSvjVolswap
            .label()
            .starts_with("73 Lyashenko-Mercurio-RFR")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("73 Lyashenko-Mercurio-RFR/Pennacchi-Sundaresan-CoCo")
    );
}

#[test]
fn suite76_tab74_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::CieslakPovalaTrendInflationCycleGabaixKrishnamurthyVigneronMbsOasErraisGieseckeGoldbergHawkesCreditMartinBondarenkoSvixLowerBound
            .label()
            .starts_with("74 Cieslak-Povala-Cycle")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("74 Cieslak-Povala-Cycle/Gabaix-Krishnamurthy-Vigneron-MBS")
    );
}

#[test]
fn suite77_tab75_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::BianchettiMercurioMultiCurveFraQuantoAcharyaAmihudBharathLiquidityStressBrigoAlfonsiCorrelatedCirCdsFukasawaAsymptoticSkewVarswap
            .label()
            .starts_with("75 Bianchetti-Mercurio-FRA-Quanto")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("75 Bianchetti-Mercurio-FRA-Quanto/Acharya-Amihud-Bharath-Liq")
    );
}

#[test]
fn suite78_tab76_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::HansonSteinGreenwoodVayanosDurationSupplyFleckensteinLongstaffLustigTipsPuzzleCollinDufresneGoldsteinYangSuperSeniorCdoTodorovOptionSpotVolatility
            .label()
            .starts_with("76 Hanson-Stein-Duration")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("76 Hanson-Stein-Duration/Fleckenstein-Longstaff-Lustig-TIPS")
    );
}

#[test]
fn suite79_tab77_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::KlinglerSyrstadSofrRepoSpikeLongstaffMithalNeisCdsBondBasisPanSingletonSovereignCdsMedvedevScailletShortMaturitySmile
            .label()
            .starts_with("77 Klingler-Syrstad-SOFR-Spike")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("77 Klingler-Syrstad-SOFR-Spike/Longstaff-Mithal-Neis-CDS-Bond")
    );
}

#[test]
fn suite80_tab78_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(TabId::BauerRudebuschQeSignalingBoyarchenkoFusterLuccaMbsOasAzizpourGieseckeSchwenklerDefaultClusteringBollerslevTodorovJumpFear.label().starts_with("78 Bauer-Rudebusch-QE"));
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.bauer_rudebusch_qe_signaling_boyarchenko_fuster_lucca_mbs_oas_azizpour_giesecke_schwenkler_default_clustering_bollerslev_todorov_jump_fear_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.bauer_rudebusch_qe_signaling_boyarchenko_fuster_lucca_mbs_oas_azizpour_giesecke_schwenkler_default_clustering_bollerslev_todorov_jump_fear_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let st = desk.bauer_rudebusch_qe_signaling_boyarchenko_fuster_lucca_mbs_oas_azizpour_giesecke_schwenkler_default_clustering_bollerslev_todorov_jump_fear_live_metric();
    assert!(n.is_finite() && (n - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("78 Bauer-Rudebusch-QE/Boyarchenko-Fuster-Lucca-MBS-OAS")
    );
}

#[test]
fn suite81_tab79_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::GurkaynakSackSwansonFomcPathBaoPanWangBondIlliquidityCollinDufresneGoldsteinMartinCreditDeterminantsAndersenFusariTodorovLeftTail
            .label()
            .starts_with("79 Gurkaynak-Sack-Swanson-FOMC")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.gurkaynak_sack_swanson_fomc_path_bao_pan_wang_bond_illiquidity_collin_dufresne_goldstein_martin_credit_determinants_andersen_fusari_todorov_left_tail_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.gurkaynak_sack_swanson_fomc_path_bao_pan_wang_bond_illiquidity_collin_dufresne_goldstein_martin_credit_determinants_andersen_fusari_todorov_left_tail_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.gurkaynak_sack_swanson_fomc_path_bao_pan_wang_bond_illiquidity_collin_dufresne_goldstein_martin_credit_determinants_andersen_fusari_todorov_left_tail_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("79 Gurkaynak-Sack-Swanson-FOMC/Bao-Pan-Wang-Bond-Illiquidity")
    );
}

#[test]
fn suite82_tab80_terminal_renders_and_reacts() {
    use ql_terminal::{DeskWorkstation, ScenarioKind, TabId};
    assert_eq!(TabId::TAB_COUNT, 80);
    assert!(
        TabId::NakamuraSteinssonFomcInfoDickNielsenFeldhutterLandoIlliquidityBerndtDuffieCreditPremiaChristoffersenHestonJacobsVarianceKernel
            .label()
            .starts_with("80 Nakamura-Steinsson-FOMC-Info")
    );
    let mut desk = DeskWorkstation::default();
    desk.scenario = ScenarioKind::Normal;
    let n = desk.nakamura_steinsson_fomc_info_dick_nielsen_feldhutter_lando_illiquidity_berndt_duffie_credit_premia_christoffersen_heston_jacobs_variance_kernel_live_metric();
    desk.scenario = ScenarioKind::CreditWideningPlus100Bp;
    let c = desk.nakamura_steinsson_fomc_info_dick_nielsen_feldhutter_lando_illiquidity_berndt_duffie_credit_premia_christoffersen_heston_jacobs_variance_kernel_live_metric();
    desk.scenario = ScenarioKind::StagflationPlus150Bp;
    let s = desk.nakamura_steinsson_fomc_info_dick_nielsen_feldhutter_lando_illiquidity_berndt_duffie_credit_premia_christoffersen_heston_jacobs_variance_kernel_live_metric();
    assert!((n - c).abs() > 1e-12 && (c - s).abs() > 1e-12);
    let frame = desk.step_tick();
    assert!(
        desk.render_dashboard(&frame, false)
            .contains("80 Nakamura-Steinsson-FOMC-Info/Dick-Nielsen-Feldhutter-Lando-Illiquidity")
    );
}
