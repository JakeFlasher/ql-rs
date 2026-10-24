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
    assert_eq!(TabId::ALL.len(), 22);
    assert_eq!(TabId::ALL[18], TabId::BermudanVolIft);
    assert_eq!(TabId::ALL[19], TabId::SaCcrSwingQp);
    let mut ws = DeskWorkstation::new();
    let frame = ws.step_tick();
    let dashboard = ws.render_dashboard(&frame, false);
    assert!(dashboard.contains("Bermudan HW") && dashboard.contains("SA-CCR EAD"));
    assert!(dashboard.contains("Commodity Swing") && dashboard.contains("Box Markowitz"));
}
