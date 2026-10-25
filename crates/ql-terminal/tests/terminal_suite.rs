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
fn phase11_tabs_and_live_telemetry_are_present() {
    assert_eq!(TabId::ALL.len(), 33);
    assert_eq!(TabId::ALL[18], TabId::BermudanVolIft);
    assert_eq!(TabId::ALL[19], TabId::SaCcrSwingQp);
    assert_eq!(TabId::ALL[22], TabId::MonotoneConvexVixParisianRobust);
    assert_eq!(TabId::ALL[24], TabId::Phase15Sota);
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
fn phase17_tab_metric_and_rates_scenario_are_live() {
    assert_eq!(TabId::ALL.len(), 33);
    assert!(TabId::ALL.contains(&TabId::Phase17Sota));
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
fn phase18_tab28_terminal_renders_and_ticks_cleanly() {
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
fn phase19_tab29_terminal_renders_and_ticks_cleanly() {
    assert_eq!(TabId::ALL.len(), 33);
    assert_eq!(
        TabId::Phase19Sota.label(),
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
fn phase20_tab30_terminal_renders_and_ticks_cleanly() {
    assert_eq!(TabId::ALL.len(), 33);
    assert_eq!(
        TabId::Phase20Sota.label(),
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
fn phase21_tab31_terminal_renders_and_ticks_cleanly() {
    assert_eq!(TabId::ALL.len(), 33);
    assert_eq!(
        TabId::Phase21Sota.label(),
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
fn phase22_tab32_terminal_renders_and_ticks_cleanly() {
    assert_eq!(TabId::ALL.len(), 33);
    assert_eq!(
        TabId::Phase22Sota.label(),
        "32 Phase 22 SOTA / KRD / LPI / CoCo / Bergomi / CGMY / AAD"
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
fn phase23_tab33_terminal_renders_and_reacts() {
    assert_eq!(TabId::ALL.len(), 33);
    assert!(TabId::ALL.contains(&TabId::Phase23Sota));
    assert_eq!(
        TabId::Phase23Sota.label(),
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
