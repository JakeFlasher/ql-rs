//! A dependency-free terminal dashboard and umbrella SDK for ql-rs.

#[allow(ambiguous_glob_reexports)]
pub mod prelude {
    pub use ql_ad::*;
    pub use ql_benches;
    pub use ql_commodities::*;
    pub use ql_credit::*;
    pub use ql_curves::*;
    pub use ql_dml::*;
    pub use ql_exotics::*;
    pub use ql_fixed_income::*;
    pub use ql_fx::*;
    pub use ql_inflation::*;
    pub use ql_instruments::*;
    pub use ql_levy::*;
    pub use ql_math::*;
    pub use ql_mc::*;
    pub use ql_models::*;
    pub use ql_pde::*;
    pub use ql_portfolio::*;
    pub use ql_rates::*;
    pub use ql_scripting::*;
    pub use ql_time::*;
    pub use ql_vol::*;
    pub use ql_xva::*;
}

use ql_ad::Real;
use ql_commodities::{SchwartzSmith2F, kirk_spread_option};
use ql_credit::{
    CalibratedCreditCurve, CdsQuote, CreditDefaultSwap, IsdaCdsEngine, ProtectionSide,
    asb_recursive_expected_tranche_loss,
};
use ql_curves::{DiscountCurve, FlatForwardCurve, YieldTermStructure};
use ql_exotics::{
    BarrierOptionSpec, BarrierType as ExoticBarrierType, MonitoringFrequency, price_barrier_option,
};
use ql_fixed_income::{ConvertibleBondSpec, price_convertible_bond_tf};
use ql_fx::ConstNotionalXccs;
use ql_fx::{FxSmileQuote, FxVannaVolgaSmile, garman_kohlhagen_spot_delta};
use ql_inflation::{CalibratedInflationCurve, MultiplicativeSeasonality, ZeroCouponInflationSwap};
use ql_levy::{VarianceGamma, cos_levy_price};
use ql_models::{HestonParams, HullWhite1F, LiftedHestonParams, lifted_heston_cos_price};
use ql_portfolio::gaussian_expected_shortfall;
use ql_portfolio::{MarketSnapshot, PortfolioDag, PortfolioRiskSummary, PortfolioTrade};
use ql_rates::{
    EuropeanSwaption, G2PlusPlus, RfrCapletSpec, RfrRateType, SwaptionSettlement,
    willems_backward_sabr_caplet_price,
};
use ql_time::Date;
use ql_vol::implied_black_vol;
use ql_xva::{
    NettingSetSwap, SimmBucket, compute_isda_simm, compute_xva, expected_positive_exposure,
    potential_future_exposure, simulate_hw1f_exposure_cube,
};
use std::fmt::Write as FmtWrite;
use std::sync::Arc;
use std::time::Instant;

const BRAILLE_DOTS: [[u32; 2]; 4] = [[0x01, 0x08], [0x02, 0x10], [0x04, 0x20], [0x40, 0x80]];

/// A 2x4 sub-pixel canvas represented by Unicode Braille cells.
#[derive(Clone, Debug)]
pub struct BrailleCanvas {
    pub cols: usize,
    pub rows: usize,
    pub width: usize,
    pub height: usize,
    pixels: Vec<bool>,
    colors: Vec<[u8; 3]>,
}
impl BrailleCanvas {
    pub fn new(cols: usize, rows: usize) -> Self {
        let width = cols.saturating_mul(2);
        let height = rows.saturating_mul(4);
        Self {
            cols,
            rows,
            width,
            height,
            pixels: vec![false; width * height],
            colors: vec![[0; 3]; width * height],
        }
    }
    pub fn clear(&mut self) {
        self.pixels.fill(false);
        self.colors.fill([0; 3]);
    }
    pub fn set_pixel(&mut self, x: usize, y: usize, color: [u8; 3]) {
        if x < self.width && y < self.height {
            let i = y * self.width + x;
            self.pixels[i] = true;
            self.colors[i] = color;
        }
    }
    pub fn draw_line(
        &mut self,
        mut x0: isize,
        mut y0: isize,
        x1: isize,
        y1: isize,
        color: [u8; 3],
    ) {
        let dx = (x1 - x0).abs();
        let sx = if x0 < x1 { 1 } else { -1 };
        let dy = -(y1 - y0).abs();
        let sy = if y0 < y1 { 1 } else { -1 };
        let mut err = dx + dy;
        loop {
            if x0 >= 0 && y0 >= 0 {
                self.set_pixel(x0 as usize, y0 as usize, color);
            }
            if x0 == x1 && y0 == y1 {
                break;
            }
            let twice = 2 * err;
            if twice >= dy {
                err += dy;
                x0 += sx;
            }
            if twice <= dx {
                err += dx;
                y0 += sy;
            }
        }
    }
    pub fn plot_series(&mut self, xs: &[f64], ys: &[f64], y_min: f64, y_max: f64, color: [u8; 3]) {
        if xs.is_empty()
            || ys.is_empty()
            || xs.len() != ys.len()
            || self.width == 0
            || self.height == 0
        {
            return;
        }
        let x_min = xs.iter().copied().fold(f64::INFINITY, f64::min);
        let x_max = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let dx = (x_max - x_min).abs().max(f64::EPSILON);
        let dy = (y_max - y_min).abs().max(f64::EPSILON);
        let canvas_width = self.width;
        let canvas_height = self.height;
        let point = |i: usize| {
            let x =
                (((xs[i] - x_min) / dx) * (canvas_width.saturating_sub(1) as f64)).round() as isize;
            let y_norm = ((ys[i] - y_min) / dy).clamp(0.0, 1.0);
            let y = ((1.0 - y_norm) * canvas_height.saturating_sub(1) as f64).round() as isize;
            (x, y)
        };
        let mut last = point(0);
        self.set_pixel(last.0.max(0) as usize, last.1.max(0) as usize, color);
        for i in 1..xs.len() {
            let next = point(i);
            self.draw_line(last.0, last.1, next.0, next.1, color);
            last = next;
        }
    }
    pub fn render_lines(&self, ansi_color: bool) -> Vec<String> {
        let mut lines = Vec::with_capacity(self.rows);
        for row in 0..self.rows {
            let mut line = String::new();
            for col in 0..self.cols {
                let mut mask = 0u32;
                let mut color = [0u8; 3];
                let mut active = false;
                for (py, dots) in BRAILLE_DOTS.iter().enumerate() {
                    for (px, dot) in dots.iter().enumerate() {
                        let x = col * 2 + px;
                        let y = row * 4 + py;
                        let idx = y * self.width + x;
                        if self.pixels[idx] {
                            mask |= *dot;
                            color = self.colors[idx];
                            active = true;
                        }
                    }
                }
                let ch = char::from_u32(0x2800 | mask).unwrap_or('\u{2800}');
                if ansi_color && active {
                    let _ = write!(
                        line,
                        "\x1b[38;2;{};{};{}m{}\x1b[0m",
                        color[0], color[1], color[2], ch
                    );
                } else {
                    line.push(ch);
                }
            }
            lines.push(line);
        }
        lines
    }
}

pub fn render_sparkline(values: &[f64], width: usize) -> String {
    if values.is_empty() || width == 0 {
        return String::new();
    }
    const LEVELS: [char; 8] = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];
    let lo = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::INFINITY, f64::min);
    let hi = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    if !lo.is_finite() || !hi.is_finite() {
        return LEVELS[0].to_string().repeat(width);
    }
    let span = (hi - lo).abs().max(f64::EPSILON);
    (0..width)
        .map(|i| {
            let j = i * values.len() / width;
            let v = values[j.min(values.len() - 1)];
            let v = if v.is_finite() { v } else { lo };
            let level = (((v - lo) / span) * 7.0).round().clamp(0.0, 7.0) as usize;
            LEVELS[level]
        })
        .collect()
}

pub fn render_diverging_bar(
    value: f64,
    max_abs: f64,
    half_width: usize,
    ansi_color: bool,
) -> String {
    let n = if max_abs > 0.0 {
        ((value.abs() / max_abs) * half_width as f64).round() as usize
    } else {
        0
    };
    let n = n.min(half_width);
    let mut out = String::new();
    let color = if value < 0.0 {
        [255, 90, 90]
    } else {
        [0, 230, 140]
    };
    if ansi_color {
        let _ = write!(out, "\x1b[38;2;{};{};{}m", color[0], color[1], color[2]);
    }
    if value < 0.0 {
        out.push_str(&" ".repeat(half_width - n));
        out.push_str(&"━".repeat(n));
    } else {
        out.push_str(&" ".repeat(half_width));
    }
    out.push('│');
    if value >= 0.0 {
        out.push_str(&"━".repeat(n));
        out.push_str(&" ".repeat(half_width - n));
    } else {
        out.push_str(&" ".repeat(half_width));
    }
    if ansi_color {
        out.push_str("\x1b[0m");
    }
    out
}

fn visible_width(s: &str) -> usize {
    let mut count = 0;
    let mut escape = false;
    for ch in s.chars() {
        if escape {
            if ch == 'm' {
                escape = false;
            }
        } else if ch == '\x1b' {
            escape = true;
        } else {
            count += 1;
        }
    }
    count
}

/// Pad a dashboard cell by terminal-visible characters, ignoring ANSI color
/// sequences.  This keeps the 2x2 layout aligned in both plain and color mode.
fn pad_visible(s: &str, width: usize) -> String {
    let visible = visible_width(s);
    if visible >= width {
        return s.to_string();
    }
    let mut out = s.to_string();
    out.push_str(&" ".repeat(width - visible));
    out
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ScenarioKind {
    Normal,
    RatesPlus50Bp,
    EquityCrashMinus15Pct,
    CreditWideningPlus100Bp,
    StagflationPlus150Bp,
}
impl ScenarioKind {
    pub fn next(self) -> Self {
        match self {
            Self::Normal => Self::RatesPlus50Bp,
            Self::RatesPlus50Bp => Self::EquityCrashMinus15Pct,
            Self::EquityCrashMinus15Pct => Self::CreditWideningPlus100Bp,
            Self::CreditWideningPlus100Bp => Self::StagflationPlus150Bp,
            Self::StagflationPlus150Bp => Self::Normal,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::Normal => "NORMAL",
            Self::RatesPlus50Bp => "RATES +50BP",
            Self::EquityCrashMinus15Pct => "EQUITY -15%",
            Self::CreditWideningPlus100Bp => "CREDIT +100BP",
            Self::StagflationPlus150Bp => "STAGFLATION +150BP",
        }
    }
}
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ChartMode {
    YieldAndInflationCurves,
    RoughHestonAndFxSmiles,
    XvaExposureProfile,
}
impl ChartMode {
    pub fn next(self) -> Self {
        match self {
            Self::YieldAndInflationCurves => Self::RoughHestonAndFxSmiles,
            Self::RoughHestonAndFxSmiles => Self::XvaExposureProfile,
            Self::XvaExposureProfile => Self::YieldAndInflationCurves,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::YieldAndInflationCurves => "Yield & Inflation Curves",
            Self::RoughHestonAndFxSmiles => "Rough Heston & FX Smiles",
            Self::XvaExposureProfile => "XVA Exposure Profile",
        }
    }
}

#[derive(Clone, Debug)]
pub struct WorkstationFrameData {
    pub tick: usize,
    pub scenario: ScenarioKind,
    pub chart_mode: ChartMode,
    pub eval_micros: f64,
    pub dag_recomputed_trades: usize,
    pub dag_total_trades: usize,
    pub dag_curves_rebuilt: bool,
    pub total_npv: f64,
    pub intraday_pnl: f64,
    pub eq_delta: f64,
    pub eq_vega: f64,
    pub fx_spot_delta: f64,
    pub inflation_par_delta_sum: f64,
    pub cva: f64,
    pub dva: f64,
    pub fva: f64,
    pub mva: f64,
    pub simm_initial_margin: f64,
    pub ir_dv01_buckets: Vec<(String, f64)>,
    pub credit_cs01_buckets: Vec<(String, f64)>,
    pub jump_to_default: f64,
    pub blotter_rows: Vec<(String, f64, String)>,
}

fn range2(a: &[f64], b: &[f64]) -> (f64, f64) {
    let lo = a.iter().chain(b).copied().fold(f64::INFINITY, f64::min);
    let hi = a.iter().chain(b).copied().fold(f64::NEG_INFINITY, f64::max);
    let pad = ((hi - lo).abs() * 0.08).max(1e-6);
    (lo - pad, hi + pad)
}

fn range3(a: &[f64], b: &[f64], c: &[f64]) -> (f64, f64) {
    let lo = a
        .iter()
        .chain(b)
        .chain(c)
        .copied()
        .fold(f64::INFINITY, f64::min);
    let hi = a
        .iter()
        .chain(b)
        .chain(c)
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let pad = ((hi - lo).abs() * 0.08).max(1e-6);
    (lo - pad, hi + pad)
}

pub struct DeskWorkstation {
    pub dag: PortfolioDag,
    pub base_market: MarketSnapshot,
    pub scenario: ScenarioKind,
    pub chart_mode: ChartMode,
    pub npv_history: Vec<f64>,
    trades: Vec<PortfolioTrade>,
    last_summary: PortfolioRiskSummary,
    tick: usize,
}
impl DeskWorkstation {
    pub fn new() -> Self {
        let ref_date = Date::from_ymd(2024, 1, 1).expect("valid workstation date");
        let market = MarketSnapshot {
            ref_date,
            ois_quotes: [0.02, 0.021, 0.022],
            fwd_quotes: [0.025, 0.026, 0.027],
            heston_params: HestonParams {
                v0: 0.04,
                kappa: 1.5,
                theta: 0.04,
                sigma_v: 0.3,
                rho: -0.7,
            },
            eq_spot: 100.0,
            eq_vol: 0.2,
        };
        let program =
            ql_scripting::BytecodeProgram::compile(&ql_scripting::ScriptedProduct::autocallable(
                vec![0.25, 0.5, 0.75, 1.0],
                105.0,
                95.0,
                0.04,
                100.0,
                85.0,
                2.0,
            ));
        let trades = vec![
            PortfolioTrade::IrsSwap {
                id: 0,
                notional: 1_000_000.0,
                fixed_rate: 0.02,
                maturity_years: 5,
            },
            PortfolioTrade::IrsSwap {
                id: 1,
                notional: 1_000_000.0,
                fixed_rate: 0.021,
                maturity_years: 10,
            },
            PortfolioTrade::IrsSwap {
                id: 2,
                notional: 1_000_000.0,
                fixed_rate: 0.022,
                maturity_years: 20,
            },
            PortfolioTrade::HestonEuropeanOption {
                id: 3,
                strike: 100.0,
                expiry: 1.0,
                is_call: true,
                notional: 1.0,
            },
            PortfolioTrade::HestonEuropeanOption {
                id: 4,
                strike: 100.0,
                expiry: 2.0,
                is_call: false,
                notional: 1.0,
            },
            PortfolioTrade::ScriptedNote {
                id: 5,
                program,
                expiry: 1.0,
            },
        ];
        let mut dag = PortfolioDag::new();
        let summary = dag.evaluate_portfolio(&market, &trades);
        Self {
            dag,
            base_market: market,
            scenario: ScenarioKind::Normal,
            chart_mode: ChartMode::YieldAndInflationCurves,
            npv_history: vec![summary.total_npv],
            trades,
            last_summary: summary,
            tick: 0,
        }
    }
    pub fn set_scenario(&mut self, scenario: ScenarioKind) {
        self.scenario = scenario;
    }
    pub fn set_chart_mode(&mut self, mode: ChartMode) {
        self.chart_mode = mode;
    }
    fn scenario_market(&self) -> MarketSnapshot {
        let mut s = self.base_market.clone();
        match self.scenario {
            ScenarioKind::Normal => {
                s.eq_spot += self.tick as f64 * 0.0001;
            }
            ScenarioKind::RatesPlus50Bp => {
                s.ois_quotes.iter_mut().for_each(|x| *x += 0.005);
                s.fwd_quotes.iter_mut().for_each(|x| *x += 0.005);
            }
            ScenarioKind::EquityCrashMinus15Pct => {
                s.eq_spot *= 0.85;
                s.heston_params.v0 *= 1.5;
            }
            ScenarioKind::CreditWideningPlus100Bp => {}
            ScenarioKind::StagflationPlus150Bp => {
                s.ois_quotes.iter_mut().for_each(|x| *x += 0.0025);
            }
        }
        s
    }
    fn credit_metrics(&self, scenario: ScenarioKind) -> (Vec<(String, f64)>, f64, f64, f64) {
        let d = self.base_market.ref_date;
        let mats = [1.0, 2.0, 3.0, 5.0, 10.0];
        let quotes: Vec<_> = mats
            .iter()
            .enumerate()
            .map(|(i, &maturity_years)| CdsQuote {
                maturity_years,
                par_spread: 0.01
                    + i as f64 * 0.001
                    + if scenario == ScenarioKind::CreditWideningPlus100Bp {
                        0.01
                    } else {
                        0.0
                    },
            })
            .collect();
        let discount = DiscountCurve::new(
            d,
            mats.to_vec(),
            mats.iter().map(|t| (-0.02 * t).exp()).collect(),
        );
        let calibrated = CalibratedCreditCurve::bootstrap(d, discount.clone(), 0.4, &quotes)
            .expect("credit calibration");
        let cds = CreditDefaultSwap {
            side: ProtectionSide::Buyer,
            notional: 1_000_000.0,
            running_spread: 0.005,
            recovery_rate: 0.4,
            payment_times: (1..=20).map(|i| i as f64 / 4.0).collect(),
            accrual_fractions: vec![0.25; 20],
        };
        let (_, cs01) = calibrated.value_and_par_cs01(|curve| {
            let dc = FlatForwardCurve::new(d, ql_ad::TapeVar::from_f64(0.02));
            IsdaCdsEngine::price(&cds, &dc, curve).npv
        });
        let labels = ["1Y", "2Y", "3Y", "5Y", "10Y"];
        let buckets = cs01
            .iter()
            .zip(labels)
            .map(|(x, label)| (label.to_string(), *x))
            .collect();
        let jtd = CalibratedCreditCurve::jump_to_default(&cds, 0.0);
        let cds_npv = IsdaCdsEngine::price(&cds, &discount, &calibrated.curve).npv;
        (
            buckets,
            jtd,
            calibrated.curve.hazard_rates.iter().sum(),
            cds_npv,
        )
    }
    fn inflation_metrics(&self, scenario: ScenarioKind) -> (f64, f64) {
        let d = self.base_market.ref_date;
        let discount = DiscountCurve::new(
            d,
            vec![1.0, 2.0, 5.0],
            vec![(-0.02f64).exp(), (-0.04f64).exp(), (-0.1f64).exp()],
        );
        let bump = if scenario == ScenarioKind::StagflationPlus150Bp {
            0.015
        } else {
            0.0
        };
        let quotes = [(1.0, 0.02 + bump), (2.0, 0.022 + bump), (5.0, 0.025 + bump)];
        let mut f = [1.0; 12];
        for (i, factor) in f.iter_mut().enumerate().take(11) {
            *factor = 1.0 + 0.006 * (i as f64 - 5.5) / 5.5;
        }
        f[11] = 1.0 / f[..11].iter().product::<f64>();
        let seasonality = MultiplicativeSeasonality::new(f).expect("seasonality");
        let curve = CalibratedInflationCurve::bootstrap_with_seasonality(
            d,
            100.0,
            discount.clone(),
            &quotes,
            Some(seasonality),
        )
        .expect("inflation calibration");
        let swap = ZeroCouponInflationSwap {
            notional: 1_000_000.0,
            fixed_rate: 0.021,
            maturity_years: 3.0,
            target_month: Some(7),
            pay_fixed: true,
        };
        let (_, delta) = curve.value_and_par_inflation_delta(|c| {
            let dc = discount.to_ad_curve(
                &discount
                    .discounts
                    .iter()
                    .map(|&x| ql_ad::TapeVar::from_f64(x))
                    .collect::<Vec<_>>(),
            );
            swap.npv(&dc, c)
        });
        let sum: f64 = delta.iter().sum::<f64>() * 1e-4;
        let npv = swap.npv(&discount, &curve.curve).val();
        (sum, npv)
    }
    fn inflation_calibration(
        &self,
        scenario: ScenarioKind,
    ) -> (DiscountCurve<f64>, CalibratedInflationCurve) {
        let d = self.base_market.ref_date;
        let discount = DiscountCurve::new(
            d,
            vec![1.0, 2.0, 5.0, 10.0],
            vec![
                (-0.02f64).exp(),
                (-0.04f64).exp(),
                (-0.1f64).exp(),
                (-0.2f64).exp(),
            ],
        );
        let bump = if scenario == ScenarioKind::StagflationPlus150Bp {
            0.015
        } else {
            0.0
        };
        let quotes = [
            (1.0, 0.02 + bump),
            (2.0, 0.022 + bump),
            (5.0, 0.025 + bump),
            (10.0, 0.027 + bump),
        ];
        let mut factors = [1.0; 12];
        for (i, factor) in factors.iter_mut().enumerate().take(11) {
            *factor = 1.0 + 0.006 * (i as f64 - 5.5) / 5.5;
        }
        factors[11] = 1.0 / factors[..11].iter().product::<f64>();
        let seasonality = MultiplicativeSeasonality::new(factors).expect("seasonality");
        let curve = CalibratedInflationCurve::bootstrap_with_seasonality(
            d,
            100.0,
            discount.clone(),
            &quotes,
            Some(seasonality),
        )
        .expect("inflation calibration");
        (discount, curve)
    }
    fn active_curve_system(
        &self,
        market: &MarketSnapshot,
    ) -> Arc<ql_curves::CalibratedCurveSystem> {
        ql_portfolio::build_curve_system(market)
    }
    fn exposure_cube(&self) -> ql_xva::NpvCube<f64> {
        let d = self.base_market.ref_date;
        let market = self.scenario_market();
        let curve = FlatForwardCurve::new(d, market.ois_quotes[0]);
        let model = HullWhite1F {
            a: 0.05,
            sigma: 0.01 + 0.25 * market.eq_vol,
        };
        let swap = NettingSetSwap {
            notional: 1_000_000.0,
            fixed_rate: market.ois_quotes[1],
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accruals: vec![1.0; 5],
            is_payer: true,
        };
        let times: Vec<f64> = (1..=16).map(|k| k as f64 * 5.0 / 16.0).collect();
        simulate_hw1f_exposure_cube(&model, &curve, &[swap], &times, 64, 17)
    }
    fn xva_metrics(&self) -> (f64, f64, f64) {
        let d = self.base_market.ref_date;
        let curve = FlatForwardCurve::new(d, 0.02);
        let model = HullWhite1F {
            a: 0.05,
            sigma: 0.01,
        };
        let swap = NettingSetSwap {
            notional: 1_000_000.0,
            fixed_rate: 0.02,
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accruals: vec![1.0; 5],
            is_payer: true,
        };
        let cube =
            simulate_hw1f_exposure_cube(&model, &curve, &[swap], &[1.0, 2.0, 3.0, 4.0], 32, 17);
        let report = compute_xva(&cube, &curve, None, 0.02, 0.6, 0.015, 0.6, 0.005);
        (report.cva.val(), report.dva.val(), report.fva.val())
    }
    fn mva_metric(&self) -> f64 {
        let d = self.base_market.ref_date;
        let market = self.scenario_market();
        let curve = FlatForwardCurve::new(d, market.ois_quotes[0]);
        let model = HullWhite1F {
            a: 0.05,
            sigma: 0.01 + 0.25 * market.eq_vol,
        };
        let swap = NettingSetSwap {
            notional: 1_000_000.0,
            fixed_rate: market.ois_quotes[1],
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accruals: vec![1.0; 5],
            is_payer: true,
        };
        let times: Vec<f64> = (1..=12).map(|k| k as f64 * 5.0 / 12.0).collect();
        ql_xva::compute_chebyshev_dim_and_mva(&model, &curve, &[swap], &times, 64, 6, 0.01, 29).mva
    }
    pub fn step_tick(&mut self) -> WorkstationFrameData {
        self.tick += 1;
        let market = self.scenario_market();
        let before_curve = self.last_summary.curve_recomputes;
        let before_trades = self.last_summary.trade_recomputes;
        let start = Instant::now();
        let summary = self.dag.evaluate_portfolio(&market, &self.trades);
        let elapsed = start.elapsed().as_secs_f64() * 1e6;
        let (credit_raw, jtd, _, cds_npv) = self.credit_metrics(self.scenario);
        let (cva, dva, fva) = self.xva_metrics();
        let mva = self.mva_metric();
        let (inflation, inflation_npv) = self.inflation_metrics(self.scenario);
        let active_curves = self.active_curve_system(&market);
        let fx = FxVannaVolgaSmile::new(
            market.eq_spot,
            0.03,
            0.01,
            1.0,
            FxSmileQuote {
                atm_vol: 0.1,
                rr_25d: 0.01,
                bf_25d: 0.005,
            },
        );
        let fx_delta =
            garman_kohlhagen_spot_delta(market.eq_spot, fx.k3, fx.sigma3, 1.0, 0.03, 0.01, true)
                .val();
        // Assemble the actual cross-asset SIMM buckets from the current
        // sensitivities, rather than applying one arbitrary scalar weight to
        // the whole IR ladder.  All values are expressed in the same dollar
        // sensitivity units before the risk weights are applied.
        let mut buckets = vec![SimmBucket {
            sensitivities: summary
                .total_par_dv01
                .iter()
                .map(|x| ql_ad::TapeVar::from_f64(*x))
                .collect(),
            risk_weights: vec![50.0; summary.total_par_dv01.len()],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0; summary.total_par_dv01.len()]; summary.total_par_dv01.len()],
        }];
        buckets.push(SimmBucket {
            sensitivities: credit_raw
                .iter()
                .map(|(_, value)| ql_ad::TapeVar::from_f64(*value * 1e-4))
                .collect(),
            risk_weights: vec![100.0; credit_raw.len()],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0; credit_raw.len()]; credit_raw.len()],
        });
        buckets.push(SimmBucket {
            sensitivities: vec![ql_ad::TapeVar::from_f64(
                summary.total_eq_delta * market.eq_spot,
            )],
            risk_weights: vec![0.15],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0]],
        });
        buckets.push(SimmBucket {
            sensitivities: vec![ql_ad::TapeVar::from_f64(fx_delta * market.eq_spot)],
            risk_weights: vec![0.08],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0]],
        });
        buckets.push(SimmBucket {
            sensitivities: vec![ql_ad::TapeVar::from_f64(inflation * 1e-4)],
            risk_weights: vec![50.0],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0]],
        });
        let inter_corr = vec![
            vec![1.0, 0.25, 0.2, 0.2, 0.2],
            vec![0.25, 1.0, 0.15, 0.2, 0.15],
            vec![0.2, 0.15, 1.0, 0.25, 0.15],
            vec![0.2, 0.2, 0.25, 1.0, 0.15],
            vec![0.2, 0.15, 0.15, 0.15, 1.0],
        ];
        let simm = compute_isda_simm(&buckets, &inter_corr).val();
        let adjusted_total_npv = summary.total_npv + cds_npv + inflation_npv;
        let prev = self
            .npv_history
            .last()
            .copied()
            .unwrap_or(adjusted_total_npv);
        self.npv_history.push(adjusted_total_npv);
        self.last_summary = summary.clone();
        let rate = active_curves.curves[0].zero_rate(1.0).val();
        let eq_vega = {
            let params = HestonParams {
                v0: ql_ad::Dual::<1>::variable(market.eq_vol * market.eq_vol, 0),
                kappa: ql_ad::Dual::from_value(market.heston_params.kappa),
                theta: ql_ad::Dual::from_value(market.heston_params.theta),
                sigma_v: ql_ad::Dual::from_value(market.heston_params.sigma_v),
                rho: ql_ad::Dual::from_value(market.heston_params.rho),
            };
            let v = ql_models::heston_cos_price(
                ql_ad::Dual::from_value(market.eq_spot),
                ql_ad::Dual::from_value(100.0),
                ql_ad::Dual::from_value(1.0),
                ql_ad::Dual::from_value(rate),
                ql_ad::Dual::from_value(0.0),
                &params,
                128,
                10.0,
                true,
            );
            v.d[0] * 2.0 * market.eq_vol
        };
        let ir_buckets: Vec<(String, f64)> = summary
            .total_par_dv01
            .iter()
            .enumerate()
            .map(|(i, x)| {
                (
                    format!(
                        "{}Y",
                        [1, 2, 5, 10, 20, 30].get(i).copied().unwrap_or(i + 1)
                    ),
                    *x * 1e-4,
                )
            })
            .collect();
        let credit_buckets: Vec<(String, f64)> = credit_raw
            .iter()
            .map(|(label, value)| (label.clone(), *value * 1e-4))
            .collect();
        let mut irs_npv = 0.0;
        if let (Some(ois), Some(fwd)) = (active_curves.curves.first(), active_curves.curves.get(1))
        {
            for m in 1..=5 {
                let t = m as f64;
                let forward = fwd.discount_factor(t - 1.0) / fwd.discount_factor(t) - 1.0;
                irs_npv += (forward - 0.02) * ois.discount_factor(t) * 1_000_000.0;
            }
        }
        let heston_params = HestonParams {
            v0: market.heston_params.v0,
            kappa: market.heston_params.kappa,
            theta: market.heston_params.theta,
            sigma_v: market.heston_params.sigma_v,
            rho: market.heston_params.rho,
        };
        let heston_npv = ql_models::heston_cos_price(
            market.eq_spot,
            100.0,
            1.0,
            rate,
            0.0,
            &heston_params,
            128,
            10.0,
            true,
        );
        let script_npv = match self.trades.get(5) {
            Some(PortfolioTrade::ScriptedNote { program, .. }) => ql_scripting::ScriptMcEngine {
                n_paths: 1024,
                seed: 42,
            }
            .price_gbm(program, market.eq_spot, market.eq_vol, rate, 0.0)
            .val(),
            _ => 0.0,
        };
        let (_, inflation_curve) = self.inflation_calibration(self.scenario);
        let zciis = ZeroCouponInflationSwap {
            notional: 1_000_000.0,
            fixed_rate: 0.021,
            maturity_years: 5.0,
            target_month: Some(7),
            pay_fixed: true,
        };
        let zciis_npv = zciis
            .npv(&inflation_curve.discount, &inflation_curve.curve)
            .val();
        let g2_swaption = EuropeanSwaption {
            notional: 1_000_000.0,
            strike: 0.025,
            expiry: 1.0,
            swap_times: (1..=6).map(|x| x as f64).collect(),
            is_payer: true,
            settlement: SwaptionSettlement::Physical,
        };
        let g2 = G2PlusPlus {
            a: 0.1,
            sigma: 0.01,
            b: 0.3,
            eta: 0.015,
            rho: -0.75,
        };
        let g2_npv = g2
            .swaption_price(&active_curves.curves[0], &g2_swaption, 16)
            .val();
        let vg = VarianceGamma {
            sigma: 0.18,
            nu: 0.20,
            theta: -0.10,
        };
        let levy_npv = cos_levy_price(
            &vg,
            market.eq_spot,
            100.0,
            rate,
            0.0,
            1.0,
            ql_instruments::OptionType::Call,
            128,
        )
        .val();
        let conv_spec = ConvertibleBondSpec {
            spot: market.eq_spot,
            conversion_ratio: 1.0,
            face_value: 100.0,
            coupon_rate: 0.03,
            coupon_times: (1..=3).map(|x| x as f64).collect(),
            maturity: 3.0,
            rate,
            div_yield: 0.0,
            credit_spread: 0.02,
            vol: market.eq_vol,
            call_price: None,
            put_price: None,
        };
        let (conv_npv, _) = price_convertible_bond_tf(&conv_spec, 48);
        let rfr_spec = RfrCapletSpec {
            notional: 1_000_000.0,
            forward_rate: 0.045,
            strike: 0.04,
            tau_0: 0.5,
            tau_1: 1.0,
            accrual_fraction: 0.5,
            discount_tau_1: active_curves.curves[0].discount_factor(1.0).val(),
            is_caplet: true,
            rate_type: RfrRateType::BackwardLooking,
        };
        let rfr_npv =
            willems_backward_sabr_caplet_price(&rfr_spec, 0.08, 0.5, -0.35, 0.45, 1.0).val();
        let xccs = ConstNotionalXccs {
            dom_notional: 1_100_000.0,
            for_notional: 1_000_000.0,
            spot_fx: 1.1,
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accrual_fractions: vec![1.0; 5],
            for_basis_spread: 0.001,
            dom_spread: 0.0,
            receive_foreign: true,
        };
        let xccs_npv = xccs
            .npv(
                &active_curves.curves[0],
                &active_curves.curves[1],
                &active_curves.curves[1],
            )
            .val();
        let es = gaussian_expected_shortfall(
            ql_ad::TapeVar::from_f64(summary.total_npv),
            ql_ad::TapeVar::from_f64(summary.total_eq_delta.abs().max(1.0)),
            0.975,
        )
        .val();
        let rough_params = LiftedHestonParams {
            v0: market.heston_params.v0,
            lambda: market.heston_params.kappa,
            theta: market.heston_params.theta,
            nu: market.heston_params.sigma_v,
            rho: market.heston_params.rho,
            hurst: 0.1,
            n_factors: 8,
            r_n: 2.5,
        };
        let rough_npv = lifted_heston_cos_price(
            market.eq_spot,
            100.0,
            1.0,
            rate,
            0.0,
            &rough_params,
            64,
            8.0,
            16,
            true,
        )
        .val();
        let barrier_spec = BarrierOptionSpec {
            option_type: ql_instruments::OptionType::Call,
            barrier_type: ExoticBarrierType::DownOut,
            spot: market.eq_spot,
            strike: 100.0,
            barrier: 90.0,
            rebate: 0.0,
            rate,
            div_yield: 0.0,
            vol: market.eq_vol,
            expiry: 1.0,
            monitoring: MonitoringFrequency::Discrete { observations: 12 },
        };
        let barrier_npv = price_barrier_option(&barrier_spec).val();
        let commodity = SchwartzSmith2F {
            chi0: 0.0,
            xi0: market.eq_spot.ln(),
            kappa: 1.0,
            sigma_chi: 0.3,
            mu_xi: 0.0,
            sigma_xi: 0.2,
            rho: -0.3,
        };
        let commodity_futures = commodity.futures_price(1.0);
        let commodity_option = commodity
            .price_futures_option(
                1.0,
                1.5,
                commodity_futures,
                rate,
                ql_instruments::OptionType::Call,
            )
            .val();
        let crack = kirk_spread_option(
            commodity_futures,
            3.0,
            0.0,
            0.3,
            0.25,
            0.2,
            rate,
            1.0,
            ql_instruments::OptionType::Call,
        )
        .val();
        let blotter_rows = vec![
            ("IRS 5Y PayFix".to_string(), irs_npv, "IFT".to_string()),
            ("Heston Call 1Y".to_string(), heston_npv, "COS".to_string()),
            (
                "RoughHest H=.1".to_string(),
                rough_npv,
                "Lift 8-Fac".to_string(),
            ),
            ("Autocall Note".to_string(), script_npv, "AAD".to_string()),
            ("CDS 5Y BuyProt".to_string(), cds_npv, "IFT".to_string()),
            (
                "ZCIIS 5Y CPI".to_string(),
                zciis_npv,
                "Seasonal".to_string(),
            ),
            (
                "FX 25D VV Call".to_string(),
                fx.price_vanna_volga(fx.k3, true).val(),
                "Dual3".to_string(),
            ),
            (
                "EXOT-BAR-DN-OUT-1Y".to_string(),
                barrier_npv,
                "Dual5 BGK".to_string(),
            ),
            (
                "CMDTY-WTI-SS2F-1Y".to_string(),
                commodity_option + crack,
                "SS2F/Kirk".to_string(),
            ),
            ("G2++ Swaption".to_string(), g2_npv, "G2 Quad".to_string()),
            ("LEVY-VG-COS-1Y".to_string(), levy_npv, "COS VG".to_string()),
            (
                "CONV-TF-BOND-3Y".to_string(),
                conv_npv.val(),
                "TF + ASB".to_string(),
            ),
            (
                "SOFR RFR SABR 1Y".to_string(),
                rfr_npv,
                "Willems".to_string(),
            ),
            (
                "EURUSD XCCS 5Y".to_string(),
                xccs_npv,
                "IFT Basis".to_string(),
            ),
            ("Portfolio ES 97.5%".to_string(), es, "Acerbi".to_string()),
        ];
        WorkstationFrameData {
            tick: self.tick,
            scenario: self.scenario,
            chart_mode: self.chart_mode,
            eval_micros: elapsed,
            dag_recomputed_trades: summary.trade_recomputes.saturating_sub(before_trades),
            dag_total_trades: self.trades.len(),
            dag_curves_rebuilt: summary.curve_recomputes != before_curve,
            total_npv: adjusted_total_npv,
            intraday_pnl: adjusted_total_npv - prev,
            eq_delta: summary.total_eq_delta,
            eq_vega,
            fx_spot_delta: fx_delta,
            inflation_par_delta_sum: inflation,
            cva,
            dva,
            fva,
            mva,
            simm_initial_margin: simm,
            ir_dv01_buckets: ir_buckets,
            credit_cs01_buckets: credit_buckets,
            jump_to_default: jtd,
            blotter_rows,
        }
    }
    fn chart(&self, frame: &WorkstationFrameData) -> BrailleCanvas {
        let mut canvas = BrailleCanvas::new(58, 10);
        let market = self.scenario_market();
        let curves = self.active_curve_system(&market);
        let x32: Vec<f64> = (0..32).map(|i| 0.25 + i as f64 * 9.75 / 31.0).collect();
        match self.chart_mode {
            ChartMode::YieldAndInflationCurves => {
                let (_, inflation) = self.inflation_calibration(self.scenario);
                let ois: Vec<f64> = x32
                    .iter()
                    .map(|&t| curves.curves[0].zero_rate(t).val())
                    .collect();
                let fwd: Vec<f64> = x32
                    .iter()
                    .map(|&t| curves.curves[1].zero_rate(t).val())
                    .collect();
                let cpi: Vec<f64> = x32
                    .iter()
                    .map(|&t| {
                        let month = (((t * 12.0).round() as i32 - 1).rem_euclid(12) + 1) as u8;
                        ZeroCouponInflationSwap {
                            notional: 1.0,
                            fixed_rate: 0.0,
                            maturity_years: t,
                            target_month: Some(month),
                            pay_fixed: true,
                        }
                        .fair_rate(&inflation.curve)
                        .val()
                    })
                    .collect();
                let (lo, hi) = range3(&ois, &fwd, &cpi);
                canvas.plot_series(&x32, &ois, lo, hi, [0, 220, 255]);
                canvas.plot_series(&x32, &fwd, lo, hi, [0, 255, 140]);
                canvas.plot_series(&x32, &cpi, lo, hi, [255, 200, 60]);
            }
            ChartMode::RoughHestonAndFxSmiles => {
                let strikes: Vec<f64> = (0..25).map(|i| 85.0 + i as f64 * 30.0 / 24.0).collect();
                let expiry = 1.0;
                let rate = curves.curves[0].zero_rate(expiry).val();
                let discount = (-rate * expiry).exp();
                let forward = market.eq_spot * (rate * expiry).exp();
                let rough = LiftedHestonParams {
                    v0: market.heston_params.v0,
                    lambda: market.heston_params.kappa,
                    theta: market.heston_params.theta,
                    nu: market.heston_params.sigma_v,
                    rho: market.heston_params.rho,
                    hurst: 0.1,
                    n_factors: 8,
                    r_n: 2.5,
                };
                let heston: Vec<f64> = strikes
                    .iter()
                    .map(|&k| {
                        let price = lifted_heston_cos_price(
                            market.eq_spot,
                            k,
                            expiry,
                            rate,
                            0.0,
                            &rough,
                            64,
                            8.0,
                            16,
                            true,
                        )
                        .val();
                        implied_black_vol(forward, k, expiry, price / discount, true)
                            .map(|v| v.val())
                            .unwrap_or(0.0)
                    })
                    .collect();
                let fx = FxVannaVolgaSmile::new(
                    100.0,
                    0.03,
                    0.01,
                    expiry,
                    FxSmileQuote {
                        atm_vol: 0.1 + 0.1 * market.eq_vol,
                        rr_25d: 0.01,
                        bf_25d: 0.005,
                    },
                );
                let fx_vol: Vec<f64> = strikes.iter().map(|&k| fx.implied_vol(k).val()).collect();
                let (lo, hi) = range2(&heston, &fx_vol);
                canvas.plot_series(&strikes, &heston, lo, hi, [0, 220, 255]);
                canvas.plot_series(&strikes, &fx_vol, lo, hi, [255, 200, 60]);
            }
            ChartMode::XvaExposureProfile => {
                let cube = self.exposure_cube();
                let epe = expected_positive_exposure(&cube);
                let pfe = potential_future_exposure(&cube, 0.975);
                let dim = {
                    let d = self.base_market.ref_date;
                    let active = self.scenario_market();
                    let curve = FlatForwardCurve::new(d, active.ois_quotes[0]);
                    let model = HullWhite1F {
                        a: 0.05,
                        sigma: 0.01 + 0.25 * active.eq_vol,
                    };
                    let swap = NettingSetSwap {
                        notional: 1_000_000.0,
                        fixed_rate: active.ois_quotes[1],
                        payment_times: (1..=5).map(|x| x as f64).collect(),
                        accruals: vec![1.0; 5],
                        is_payer: true,
                    };
                    let times = (1..=16).map(|k| k as f64 * 5.0 / 16.0).collect::<Vec<_>>();
                    ql_xva::compute_chebyshev_dim_and_mva(
                        &model,
                        &curve,
                        &[swap],
                        &times,
                        64,
                        6,
                        0.01,
                        29,
                    )
                    .eim_profile
                };
                let xs: Vec<f64> = cube.times.clone();
                let (lo, hi) = range3(&epe, &pfe, &dim);
                canvas.plot_series(&xs, &epe, lo, hi, [0, 255, 140]);
                canvas.plot_series(&xs, &pfe, lo, hi, [255, 90, 90]);
                canvas.plot_series(&xs, &dim, lo, hi, [180, 130, 255]);
            }
        }
        let _ = frame;
        canvas
    }
    pub fn render_dashboard(&self, frame: &WorkstationFrameData, ansi_color: bool) -> String {
        const LEFT: usize = 57;
        const RIGHT: usize = 58;
        let chart_lines = self.chart(frame).render_lines(ansi_color);
        let asb_el = asb_recursive_expected_tranche_loss(
            &[0.02_f64; 20],
            &[1_usize; 20],
            1.0 / 20.0,
            0.2,
            0.03,
            0.07,
        );
        let mut left_top = Vec::with_capacity(10);
        left_top.push("1. BOOK SUMMARY, XVA & ISDA SIMM".to_string());
        left_top.push(format!(
            "Total NPV: {:>14.2} {}",
            frame.total_npv,
            render_sparkline(&self.npv_history, 16)
        ));
        left_top.push(format!("Intraday PnL: {:>11.2}", frame.intraday_pnl));
        left_top.push(format!(
            "Eq Delta: {:>12.2}  Eq Vega: {:>10.2}",
            frame.eq_delta, frame.eq_vega
        ));
        left_top.push(format!(
            "FX Delta: {:>12.5}  CPI Δ01: {:>10.2}",
            frame.fx_spot_delta, frame.inflation_par_delta_sum
        ));
        left_top.push(format!(
            "XVA CVA/DVA/FVA: {:>7.2}/{:>7.2}/{:>7.2}  MVA:{:>7.2}",
            frame.cva, frame.dva, frame.fva, frame.mva
        ));
        left_top.push(format!("ISDA SIMM IM: {:>14.2}", frame.simm_initial_margin));
        left_top.push(format!("Jump-to-default: {:>12.2}", frame.jump_to_default));
        let cached = frame
            .dag_total_trades
            .saturating_sub(frame.dag_recomputed_trades);
        left_top.push(format!(
            "DAG cache hit: {cached}/{}",
            frame.dag_total_trades
        ));
        left_top.push(format!(
            "Eval: {:>8.0} µs  Mode: {}",
            frame.eval_micros,
            frame.scenario.label()
        ));
        let right_top = chart_lines;

        let max_ir = frame
            .ir_dv01_buckets
            .iter()
            .map(|x| x.1.abs())
            .fold(0.0, f64::max)
            .max(1.0);
        let max_cs = frame
            .credit_cs01_buckets
            .iter()
            .map(|x| x.1.abs())
            .fold(0.0, f64::max)
            .max(1.0);
        let mut left_bottom = Vec::with_capacity(11);
        for (label, value) in &frame.ir_dv01_buckets {
            left_bottom.push(format!(
                "IRS DV01 {label:>3}: {:>9.2} {}",
                value,
                render_diverging_bar(*value, max_ir, 10, ansi_color)
            ));
        }
        for (label, value) in &frame.credit_cs01_buckets {
            left_bottom.push(format!(
                "CDS CS01 {label:>3}: {:>9.2} {}",
                value,
                render_diverging_bar(*value, max_cs, 10, ansi_color)
            ));
        }
        while left_bottom.len() < 20 {
            left_bottom.push(String::new());
        }
        let mut right_bottom: Vec<String> = frame
            .blotter_rows
            .iter()
            .map(|(name, value, tag)| format!("{name:<20} {:>12.2} {tag:<10}", value))
            .collect();
        right_bottom.push("SOTA IFT / COS / AAD telemetry".to_string());
        right_bottom.push("Exotics / Commodities / Burley RQMC".to_string());
        right_bottom.push("Levy COS / Convertible TF / ASB CDO".to_string());
        right_bottom.push(format!("ASB tranche EL: {:.6}", asb_el));
        right_bottom.push(format!(
            "DAG recomputed: {}/{}",
            frame.dag_recomputed_trades, frame.dag_total_trades
        ));
        right_bottom.push(format!("Chart: {}", frame.chart_mode.label()));
        while right_bottom.len() < 20 {
            right_bottom.push(String::new());
        }

        let mut out = String::new();
        let top_left = pad_visible(
            &format!("─ ql-rs QUANT RISK WORKSTATION  Tick #{:04} ", frame.tick),
            LEFT,
        );
        let top_right = pad_visible(
            &format!("─ 2. BRAILLE CANVAS: {} ", frame.chart_mode.label()),
            RIGHT,
        );
        let _ = writeln!(out, "╭{}┬{}╮", top_left, top_right);
        for i in 0..10 {
            let l = left_top.get(i).map(String::as_str).unwrap_or("");
            let r = right_top.get(i).map(String::as_str).unwrap_or("");
            let _ = writeln!(out, "│{}│{}│", pad_visible(l, LEFT), pad_visible(r, RIGHT));
        }
        let _ = writeln!(
            out,
            "├{}┼{}┤",
            pad_visible(" 3. IFT BUCKETED DV01 & CREDIT CS01 LADDERS", LEFT),
            pad_visible(" 4. CROSS-ASSET LIVE BLOTTER & SOTA TELEMETRY", RIGHT)
        );
        for i in 0..20 {
            let _ = writeln!(
                out,
                "│{}│{}│",
                pad_visible(&left_bottom[i], LEFT),
                pad_visible(&right_bottom[i], RIGHT)
            );
        }
        let _ = writeln!(out, "╰{}┴{}╯", "─".repeat(LEFT), "─".repeat(RIGHT));
        out
    }
}
impl Default for DeskWorkstation {
    fn default() -> Self {
        Self::new()
    }
}
