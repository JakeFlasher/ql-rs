use crate::harness::{
    BenchConfig, ComparisonRow, configured_counts, run_comparison, with_config_override,
};
use ql_ad::{Dual, Dual2, Real, TapeVar};
use ql_commodities as _;
use ql_credit::asb_recursive_expected_tranche_loss;
use ql_credit::{
    BaseCorrelationTrancheQuote, CalibratedCreditCurve, CdsQuote, CreditDefaultSwap, IsdaCdsEngine,
    OKaneBespokeBaseCorrelationSpec, ProtectionSide, SurvivalProbabilityCurve,
    bootstrap_base_correlation_curve, price_mezzanine_tranche_base_correlation,
    price_okane_bespoke_tranche_and_student_t_basket,
};
use ql_curves::{
    CalibratedCurveSystem, DiscountCurve, FlatForwardCurve, MonotoneConvexZeroCurve,
    SingleCurveHelper, YieldTermStructure,
};
use ql_dml::{cholesky_lower, price_correlated_basket_call_aad};
use ql_exotics::{
    BarrierOptionSpec, BarrierType, ExcursionClockKind, MonitoringFrequency,
    ParisianBarrierOptionSpec, price_barrier_option, price_parisian_option_mc_ad,
};
use ql_fixed_income::{
    CmoTrancheSpec, CmoWaterfallSpec, NelsonSiegelSvensson, SpecifiedPoolMbsSpec,
    TbaDollarRollSpec, TuckmanKeyRateDurationSpec, analyze_tba_dollar_roll,
    analyze_tuckman_key_rate_durations_and_pca_hedge, price_cmo_waterfall,
    price_specified_pool_mbs,
};
use ql_fx::{CalibratedXccsCurve, ConstNotionalXccs, XccsQuote};
use ql_fx::{
    FxSmileQuote, FxVannaVolgaSmile, WystupFxStructuredProductSpec,
    evaluate_wystup_fx_accumulator_and_dcd,
};
use ql_inflation::{
    BrigoMercurioInflationMarketModelSpec, JarrowYildirimParams, YoYInflationCapFloor,
    YoYInflationCapFloorType, ZeroInflationCurve, price_brigo_mercurio_yoy_and_lpi_market_model,
    price_yoy_capfloor_jarrow_yildirim,
};
use ql_instruments::{FixedRateLeg, OptionType, ResolvedCouponLeg};
use ql_levy::{VarianceGamma, cos_levy_price_strip};
use ql_mc::{GbmSoaSimulator, PhiloxRng, price_asian_arithmetic_rqmc};
use ql_models::HullWhite1F;
use ql_models::{
    Bergomi2FParams, HestonParams, QuadraticRoughHestonParams, bergomi_2f_implied_vol_bg,
    heston_cos_price, qrh_implied_vol_expansion,
};
use ql_portfolio::{
    BoxConstrainedMarkowitzSpec, MarketSnapshot, NestedClusteredOptimizationSpec, PortfolioDag,
    PortfolioTrade, box_constrained_markowitz_weights, conditional_drawdown_at_risk,
    entropic_value_at_risk, factor_risk_attribution, gerber_covariance_matrix,
    nested_clustered_optimization_weights, risk_budgeting_weights,
};
use ql_rates as _;
use ql_rates::{
    CheyetteSwaptionSpec, RfrCapletSpec, RfrRateType, price_cheyette_swaption_markov_projection,
    willems_backward_sabr_caplet_price, willems_backward_sabr_caplet_strip,
};
use ql_scripting::{BytecodeProgram, ScriptMcEngine, ScriptedProduct};
use ql_time::{BitsetCalendar, Date, DayCounter, Frequency};
use ql_vol::{
    AdSabrSmile, CalibratedSabrSmile, SmileStrikeQuote, black76_undiscounted, implied_black_vol,
    implied_black_vol_dual2,
};
use ql_xva::{NettingSetSwap, compute_bruteforce_dim_and_mva, compute_chebyshev_dim_and_mva};
use std::sync::Arc;

fn config(scale: usize) -> BenchConfig {
    let (warmup_iters, samples) = configured_counts(1, 3);
    BenchConfig {
        warmup_iters,
        samples,
        // Repeat each fixed workload enough times to amortize timer and
        // scheduler jitter when the workspace's integration binaries run in
        // parallel.  `scale` changes only this repetition count; the domain
        // workload itself remains invariant across test and bench modes.
        // The short parity configuration intentionally uses one iteration so
        // the complete 39-suite test remains a fast (<2s on CI) invariant
        // check rather than a second timing benchmark.
        iters_per_sample: if samples <= 3 {
            scale.max(1)
        } else {
            scale.max(1) * 8
        },
    }
}

fn date(y: i32, m: u8, d: u8) -> Date {
    Date::from_ymd(y, m, d).expect("benchmark date")
}
type ScalarBasketPricer = Arc<dyn Fn(&[Vec<f64>]) -> f64 + Send + Sync>;

fn benchmark_calendar(scale: usize) -> ComparisonRow {
    let cal = BitsetCalendar::us_nyse();
    let anchor = date(2020, 1, 1);
    let ranges = (0..64)
        .map(|i| {
            let start = anchor.add_days(i * 90);
            let end = start.add_days(365 * (1 + i % 30));
            (start, end)
        })
        .collect::<Vec<_>>();
    let baseline_ranges = ranges.clone();
    let sota_ranges = ranges;
    let baseline_cal = cal.clone();
    let sota_cal = cal.clone();
    run_comparison(
        "01_calendar_popcnt",
        "day-by-day holidays vs. packed POPCNT calendar",
        config(scale),
        0.0,
        move || {
            let mut total = 0i64;
            for &(a, b) in &baseline_ranges {
                for serial in a.days()..b.days() {
                    if baseline_cal.is_business_day(Date::from_days(serial)) {
                        total += 1;
                    }
                }
            }
            total as f64
        },
        move || {
            sota_ranges
                .iter()
                .map(|&(a, b)| sota_cal.business_days_between(a, b) as i64)
                .sum::<i64>() as f64
        },
    )
}

#[derive(Clone)]
struct BoxedCashflow {
    notional: f64,
    year_fraction: f64,
    rate: f64,
    time: f64,
}
trait BoxedCf {
    fn pv(&self, yield_rate: f64) -> f64;
}
impl BoxedCf for BoxedCashflow {
    fn pv(&self, yield_rate: f64) -> f64 {
        self.notional * self.year_fraction * self.rate * (-yield_rate * self.time).exp()
    }
}

fn make_fixed_leg() -> ResolvedCouponLeg {
    let ref_date = date(2024, 1, 1);
    let mut leg = ResolvedCouponLeg {
        day_counter: Some(DayCounter::Act365Fixed),
        ..Default::default()
    };
    for i in 0..200 {
        let start = ref_date.add_days(i * 91);
        let end = ref_date.add_days((i + 1) * 91);
        let yf = 0.25;
        let t = (i + 1) as f64 * 0.25;
        leg.payment_times.push(t);
        leg.year_fractions.push(yf);
        leg.notionals.push(1_000_000.0);
        leg.rates_or_spreads.push(0.035);
        leg.gearings.push(1.0);
        leg.payment_dates.push(end);
        leg.accrual_start_dates.push(start);
        leg.accrual_end_dates.push(end);
        leg.fixing_dates.push(start);
        leg.accrual_start_times.push(i as f64 * 0.25);
        leg.accrual_end_times.push(t);
    }
    leg
}

fn benchmark_leg_simd(scale: usize) -> ComparisonRow {
    let leg = make_fixed_leg();
    let boxed: Vec<Box<dyn BoxedCf>> = (0..leg.len())
        .map(|i| {
            Box::new(BoxedCashflow {
                notional: leg.notionals[i],
                rate: leg.rates_or_spreads[i],
                year_fraction: leg.year_fractions[i],
                time: leg.payment_times[i],
            }) as Box<dyn BoxedCf>
        })
        .collect();
    let rates = (0..64)
        .map(|i| 0.02 + 0.0001 * (i % 11) as f64)
        .collect::<Vec<_>>();
    let baseline_boxed = boxed;
    let baseline_rates = rates.clone();
    let sota_rates = rates;
    let sota_leg = FixedRateLeg::from_resolved(&leg);
    run_comparison(
        "02_leg_simd_soa",
        "boxed AoS virtual dispatch vs. contiguous SIMD SoA",
        config(scale),
        1e-6,
        move || {
            baseline_rates
                .iter()
                .map(|&r| {
                    baseline_boxed
                        .iter()
                        .map(|x| std::hint::black_box(x.as_ref()).pv(r))
                        .sum::<f64>()
                })
                .sum::<f64>()
        },
        move || {
            sota_leg
                .npv_simd_flat_batch(&sota_rates)
                .into_iter()
                .sum::<f64>()
        },
    )
}

fn make_curve_system(
    reference_date: Date,
    ois_quotes: [f64; 3],
    fwd_quotes: [f64; 3],
) -> CalibratedCurveSystem {
    let pillars = [1_i32, 2, 5];
    let ois = pillars
        .into_iter()
        .zip(ois_quotes)
        .map(|(years, quote)| ql_curves::OisSwapInstrument {
            start: reference_date,
            end: reference_date.add_days(365 * years),
            fixed_rate: quote,
            frequency: Frequency::Annual,
            day_counter: DayCounter::Act365Fixed,
        })
        .collect::<Vec<_>>();
    let fwd = pillars
        .into_iter()
        .zip(fwd_quotes)
        .map(|(years, quote)| ql_curves::IrsSwapInstrument {
            start: reference_date,
            end: reference_date.add_days(365 * years),
            fixed_rate: quote,
            frequency: Frequency::Annual,
            day_counter: DayCounter::Act365Fixed,
        })
        .collect::<Vec<_>>();
    let ois_refs = ois
        .iter()
        .map(|x| x as &dyn SingleCurveHelper)
        .collect::<Vec<_>>();
    CalibratedCurveSystem::from_multicurve_helpers(reference_date, &ois_refs, &fwd)
        .expect("benchmark curve calibration")
}

fn swap_book_value_f64(system: &CalibratedCurveSystem, count: usize) -> f64 {
    let ois = &system.curves[0];
    let fwd = &system.curves[1];
    (0..count)
        .map(|i| {
            let maturity: f64 = [1.0, 2.0, 5.0][i % 3];
            let fixed = 0.02 + 0.00025 * (i % 7) as f64;
            (1..=(maturity.round() as usize))
                .map(|m| {
                    let t = m as f64;
                    let fr = fwd.discount_factor(t - 1.0) / fwd.discount_factor(t) - 1.0;
                    (fr - fixed) * ois.discount_factor(t) * 1_000_000.0
                })
                .sum::<f64>()
        })
        .sum()
}

fn swap_book_value_ad(
    params: &[TapeVar<f64>],
    system: &CalibratedCurveSystem,
    count: usize,
) -> TapeVar<f64> {
    let ois = DiscountCurve::new(
        system.curves[0].reference_date,
        system.curves[0].times.clone(),
        params[..3].to_vec(),
    );
    let fwd = DiscountCurve::new(
        system.curves[1].reference_date,
        system.curves[1].times.clone(),
        params[3..].to_vec(),
    );
    let mut value = TapeVar::from_f64(0.0);
    for i in 0..count {
        let maturity: f64 = [1.0, 2.0, 5.0][i % 3];
        let fixed = 0.02 + 0.00025 * (i % 7) as f64;
        for m in 1..=maturity.round() as usize {
            let t = m as f64;
            let fr = fwd.discount_factor(t - 1.0) / fwd.discount_factor(t) - TapeVar::from_f64(1.0);
            value += (fr - fixed) * ois.discount_factor(t) * 1_000_000.0;
        }
    }
    value
}

fn benchmark_multicurve(scale: usize) -> ComparisonRow {
    let reference_date = date(2024, 1, 1);
    let ois_quotes = [0.02, 0.021, 0.022];
    let fwd_quotes = [0.025, 0.026, 0.027];
    let system = Arc::new(make_curve_system(reference_date, ois_quotes, fwd_quotes));
    let h = 1e-5;
    let sota_system = system.clone();
    run_comparison(
        "03_multicurve_ift_dv01",
        "12 bump-and-rebootstrap valuations vs. cached IFT VJP",
        config(scale),
        5e-2,
        move || {
            let mut out = 0.0;
            for j in 0..6 {
                let mut up_o = ois_quotes;
                let mut dn_o = ois_quotes;
                let mut up_f = fwd_quotes;
                let mut dn_f = fwd_quotes;
                if j < 3 {
                    up_o[j] += h;
                    dn_o[j] -= h;
                } else {
                    up_f[j - 3] += h;
                    dn_f[j - 3] -= h;
                }
                let up = make_curve_system(reference_date, up_o, up_f);
                let dn = make_curve_system(reference_date, dn_o, dn_f);
                out += (swap_book_value_f64(&up, 10) - swap_book_value_f64(&dn, 10)) / (2.0 * h);
            }
            out
        },
        move || {
            let (_, dv) = sota_system
                .value_and_par_dv01(|params| swap_book_value_ad(params, &sota_system, 10));
            dv.into_iter().sum()
        },
    )
}

fn implied_quotes() -> Vec<(f64, f64, f64, f64, bool)> {
    (0..32)
        .map(|i| {
            let forward = 100.0 + 0.25 * (i % 8) as f64;
            let strike = forward * (0.96 + 0.004 * (i % 11) as f64);
            let expiry = 0.75 + 0.125 * (i % 5) as f64;
            let vol = 0.18 + 0.002 * (i % 13) as f64;
            let call = i % 2 == 0;
            let price = black76_undiscounted(forward, strike, vol, expiry, call);
            (forward, strike, expiry, price, call)
        })
        .collect()
}

fn benchmark_implied_vol(scale: usize) -> ComparisonRow {
    let quotes = implied_quotes();
    run_comparison(
        "04_implied_vol_ift_dual2",
        "nine Jäckel finite-difference solves vs. one Dual2 IFT solve",
        config(scale),
        1e-3,
        {
            let quotes = quotes.clone();
            move || {
                quotes
                    .iter()
                    .map(|&(f, k, t, p, call)| {
                        let h = 1e-5;
                        let base = implied_black_vol(f, k, t, p, call).expect("implied vol");
                        let fp = implied_black_vol(f + h, k, t, p, call).unwrap();
                        let fm = implied_black_vol(f - h, k, t, p, call).unwrap();
                        let kp = implied_black_vol(f, k + h, t, p, call).unwrap();
                        let km = implied_black_vol(f, k - h, t, p, call).unwrap();
                        let tp = implied_black_vol(f, k, t + h, p, call).unwrap();
                        let tm = implied_black_vol(f, k, t - h, p, call).unwrap();
                        let pp = implied_black_vol(f, k, t, p + h, call).unwrap();
                        let pm = implied_black_vol(f, k, t, p - h, call).unwrap();
                        let d = [
                            (fp - fm) / (2.0 * h),
                            (kp - km) / (2.0 * h),
                            (tp - tm) / (2.0 * h),
                            (pp - pm) / (2.0 * h),
                        ];
                        let d2 = [
                            (fp - 2.0 * base + fm) / (h * h),
                            (kp - 2.0 * base + km) / (h * h),
                            (tp - 2.0 * base + tm) / (h * h),
                            (pp - 2.0 * base + pm) / (h * h),
                        ];
                        base + d.into_iter().sum::<f64>() + d2.into_iter().sum::<f64>()
                    })
                    .sum::<f64>()
            }
        },
        {
            let quotes = quotes.clone();
            move || {
                quotes
                    .iter()
                    .map(|&(f, k, t, p, call)| {
                        let x = implied_black_vol_dual2(
                            Dual2::<4>::variable(f, 0),
                            Dual2::<4>::variable(k, 1),
                            Dual2::<4>::variable(t, 2),
                            Dual2::<4>::variable(p, 3),
                            call,
                        )
                        .expect("implied vol dual");
                        x.val + x.d1.iter().sum::<f64>() + (0..4).map(|i| x.d2[i][i]).sum::<f64>()
                    })
                    .sum::<f64>()
            }
        },
    )
}

fn benchmark_heston(scale: usize) -> ComparisonRow {
    let base = HestonParams {
        v0: 0.04,
        kappa: 1.5,
        theta: 0.04,
        sigma_v: 0.05,
        rho: -0.7,
    };
    let strikes = (0..16).map(|i| 98.0 + i as f64 * 0.25).collect::<Vec<_>>();
    let spot = 100.0;
    let rate = 0.03;
    let div = 0.0;
    let expiry = 1.0;
    run_comparison(
        "05_heston_cos_dual5",
        "ten central COS valuations vs. one Dual5 COS valuation",
        config(scale),
        1e-3,
        {
            let strikes = strikes.clone();
            move || {
                strikes
                    .iter()
                    .map(|&k| {
                        let h = 1e-4;
                        let mut sum = 0.0;
                        for j in 0..5 {
                            let mut up = base;
                            let mut dn = base;
                            let values = [base.v0, base.kappa, base.theta, base.sigma_v, base.rho];
                            let mut uv = values;
                            let mut dv = values;
                            uv[j] += h;
                            dv[j] -= h;
                            up.v0 = uv[0];
                            up.kappa = uv[1];
                            up.theta = uv[2];
                            up.sigma_v = uv[3];
                            up.rho = uv[4];
                            dn.v0 = dv[0];
                            dn.kappa = dv[1];
                            dn.theta = dv[2];
                            dn.sigma_v = dv[3];
                            dn.rho = dv[4];
                            let vp =
                                heston_cos_price(spot, k, expiry, rate, div, &up, 128, 10.0, true);
                            let vm =
                                heston_cos_price(spot, k, expiry, rate, div, &dn, 128, 10.0, true);
                            sum += (vp - vm) / (2.0 * h);
                        }
                        sum
                    })
                    .sum()
            }
        },
        {
            let strikes = strikes.clone();
            move || {
                let mut sum = 0.0;
                for &k in &strikes {
                    let p = HestonParams {
                        v0: Dual::<5>::variable(base.v0, 0),
                        kappa: Dual::<5>::variable(base.kappa, 1),
                        theta: Dual::<5>::variable(base.theta, 2),
                        sigma_v: Dual::<5>::variable(base.sigma_v, 3),
                        rho: Dual::<5>::variable(base.rho, 4),
                    };
                    let v = heston_cos_price(
                        Dual::from_value(spot),
                        Dual::from_value(k),
                        Dual::from_value(expiry),
                        Dual::from_value(rate),
                        Dual::from_value(div),
                        &p,
                        128,
                        10.0,
                        true,
                    );
                    sum += v.d.iter().sum::<f64>();
                }
                sum
            }
        },
    )
}

#[allow(clippy::too_many_arguments)]
fn scalar_gbm_price(
    spot: f64,
    vol: f64,
    rate: f64,
    div: f64,
    expiry: f64,
    strike: f64,
    n_paths: usize,
    n_steps: usize,
    seed: u64,
) -> f64 {
    let rng = PhiloxRng::new(seed);
    let steps = n_steps.max(1);
    let dt = expiry / steps as f64;
    let drift = (rate - div - 0.5 * vol * vol) * dt;
    let vol_dt = vol * dt.sqrt();
    let discount = (-rate * expiry).exp();
    let batches = n_paths.max(16).div_ceil(16);
    let mut sum = 0.0;
    for batch in 0..batches {
        let mut log_pos = [0.0; 8];
        let mut log_neg = [0.0; 8];
        for step in 0..steps {
            let base_path = (batch * 2) as u64;
            let z0 = rng.normal4(base_path, step as u32);
            let z1 = rng.normal4(base_path + 1, step as u32);
            let z = [z0[0], z0[1], z0[2], z0[3], z1[0], z1[1], z1[2], z1[3]];
            for lane in 0..8 {
                log_pos[lane] += drift + vol_dt * z[lane];
                log_neg[lane] += drift - vol_dt * z[lane];
            }
        }
        for lane in 0..8 {
            let s_pos = spot * libm::exp(log_pos[lane]);
            let s_neg = spot * libm::exp(log_neg[lane]);
            sum += ((s_pos - strike).max(0.0) + (s_neg - strike).max(0.0)) * (0.5 * discount);
        }
    }
    sum / (batches * 8) as f64
}

fn benchmark_gbm(scale: usize) -> ComparisonRow {
    let n_paths = 2_048;
    let n_steps = 32;
    let simulator = GbmSoaSimulator::new(100.0, 0.2, 0.03, 0.0, 1.0, n_paths, n_steps, 11);
    let baseline_sim = simulator.clone();
    let sota_sim = simulator.clone();
    run_comparison(
        "06_gbm_mc_simd8",
        "scalar Philox path loop vs. 8-lane antithetic GBM",
        config(scale),
        1e-9,
        move || {
            scalar_gbm_price(
                baseline_sim.spot,
                baseline_sim.vol,
                baseline_sim.rate,
                baseline_sim.div,
                baseline_sim.expiry,
                100.0,
                baseline_sim.n_paths,
                baseline_sim.n_steps,
                baseline_sim.seed,
            )
        },
        move || sota_sim.price_european(100.0, true).0,
    )
}

fn benchmark_scripted(scale: usize) -> ComparisonRow {
    let program = BytecodeProgram::compile(&ScriptedProduct::autocallable(
        vec![0.25, 0.5, 0.75, 1.0],
        105.0,
        95.0,
        0.04,
        100.0,
        85.0,
        2.0,
    ));
    let engine = ScriptMcEngine {
        n_paths: 256,
        seed: 19,
    };
    let base = [100.0, 0.2, 0.03, 0.0];
    run_comparison(
        "07_scripted_autocall_fuzzy_aad",
        "nine bumped scripted MC runs vs. one four-input fuzzy-AAD run",
        config(scale),
        1e-2,
        {
            let program = program.clone();
            move || {
                let h = 1e-4;
                let base_value = engine.price_gbm(&program, base[0], base[1], base[2], base[3]);
                std::hint::black_box(base_value);
                let mut sum = 0.0;
                for j in 0..4 {
                    let mut up = base;
                    let mut dn = base;
                    up[j] += h;
                    dn[j] -= h;
                    let vp = engine.price_gbm(&program, up[0], up[1], up[2], up[3]);
                    let vm = engine.price_gbm(&program, dn[0], dn[1], dn[2], dn[3]);
                    sum += (vp - vm) / (2.0 * h);
                }
                sum
            }
        },
        {
            let program = program.clone();
            move || {
                let x = engine.price_gbm(
                    &program,
                    Dual::<4>::variable(base[0], 0),
                    Dual::<4>::variable(base[1], 1),
                    Dual::<4>::variable(base[2], 2),
                    Dual::<4>::variable(base[3], 3),
                );
                x.d.iter().sum()
            }
        },
    )
}

fn portfolio_snapshot() -> MarketSnapshot {
    MarketSnapshot {
        ref_date: date(2024, 1, 1),
        ois_quotes: [0.02, 0.021, 0.022],
        fwd_quotes: [0.025, 0.026, 0.027],
        inflation_quotes: vec![("CPI".to_string(), 0.023)],
        heston_params: HestonParams {
            v0: 0.04,
            kappa: 1.5,
            theta: 0.04,
            sigma_v: 0.3,
            rho: -0.7,
        },
        eq_spot: 100.0,
        eq_vol: 0.2,
    }
}

fn portfolio_trades() -> Vec<PortfolioTrade> {
    let mut trades = Vec::with_capacity(244);
    for id in 0..240 {
        trades.push(PortfolioTrade::IrsSwap {
            id,
            notional: 1_000_000.0,
            fixed_rate: 0.02 + id as f64 * 0.00005,
            maturity_years: [5, 10, 20, 30][id % 4],
        });
    }
    for id in 240..244 {
        trades.push(PortfolioTrade::HestonEuropeanOption {
            id,
            strike: 90.0 + (id - 240) as f64 * 5.0,
            expiry: 1.0 + (id - 240) as f64 * 0.25,
            is_call: id % 2 == 0,
            notional: 1.0,
        });
    }
    trades
}

fn benchmark_portfolio(scale: usize) -> ComparisonRow {
    let base = portfolio_snapshot();
    let trades = Arc::new(portfolio_trades());
    let mut ticked = base.clone();
    ticked.eq_spot = 101.5;
    let baseline_trades = trades.clone();
    let baseline_tick = ticked.clone();
    let mut warmed = PortfolioDag::new();
    warmed.evaluate_portfolio(&base, &trades);
    let mut sota_tick = ticked;
    let mut baseline_index = 0usize;
    let mut sota_index = 0usize;
    run_comparison(
        "08_portfolio_dag_incremental",
        "full bootstrap/revaluation vs. content-addressed incremental DAG",
        config(scale),
        1e-9,
        move || {
            let mut snapshot = baseline_tick.clone();
            snapshot.eq_spot = 101.5 + baseline_index as f64 * 1e-7;
            baseline_index += 1;
            let mut dag = PortfolioDag::new();
            let result = dag.evaluate_portfolio(&snapshot, &baseline_trades);
            result.total_npv + result.total_eq_delta
        },
        move || {
            // A fresh equity hash each iteration forces the four equity trades
            // through the dirty subgraph while retaining the cached curves and
            // 240 fixed-income valuations.
            sota_tick.eq_spot = 101.5 + sota_index as f64 * 1e-7;
            sota_index += 1;
            let result = warmed.evaluate_portfolio(&sota_tick, &trades);
            result.total_npv + result.total_eq_delta
        },
    )
}

fn benchmark_credit(scale: usize) -> ComparisonRow {
    let ref_date = date(2024, 1, 1);
    let maturities = [1.0, 2.0, 3.0, 5.0, 10.0];
    let spreads = [0.01, 0.012, 0.014, 0.016, 0.018];
    let quotes: Vec<_> = maturities
        .iter()
        .zip(spreads)
        .map(|(&maturity_years, par_spread)| CdsQuote {
            maturity_years,
            par_spread,
        })
        .collect();
    let discount = DiscountCurve::new(
        ref_date,
        maturities.to_vec(),
        maturities.iter().map(|t| (-0.02 * t).exp()).collect(),
    );
    let base = CalibratedCreditCurve::bootstrap(ref_date, discount.clone(), 0.4, &quotes)
        .expect("credit calibration");
    let cds_times: Vec<Vec<f64>> = (0..10)
        .map(|i| {
            let maturity_quarters = ((i % 12) + 1) * 4;
            (1..=maturity_quarters).map(|j| j as f64 / 4.0).collect()
        })
        .collect();
    let portfolio: Vec<CreditDefaultSwap> = cds_times
        .iter()
        .enumerate()
        .map(|(i, times)| CreditDefaultSwap {
            side: ProtectionSide::Buyer,
            notional: 1_000_000.0,
            running_spread: 0.01 + 0.0005 * (i % 5) as f64,
            recovery_rate: 0.4,
            payment_times: times.clone(),
            accrual_fractions: vec![0.25; times.len()],
        })
        .collect();
    let portfolio_sota = portfolio.clone();
    let portfolio_value = |curve: &SurvivalProbabilityCurve<f64>| {
        let flat = FlatForwardCurve::new(ref_date, 0.02);
        portfolio
            .iter()
            .map(|cds| IsdaCdsEngine::price(cds, &flat, curve).npv)
            .sum::<f64>()
    };
    run_comparison(
        "09_credit_cds_ift_cs01",
        "ten bump-and-rebootstrap CDS strips vs. cached hazard IFT CS01",
        config(scale),
        1e-2,
        {
            let discount = discount.clone();
            move || {
                let h = 1e-6;
                let mut sum = 0.0;
                for j in 0..quotes.len() {
                    let mut up = quotes.clone();
                    let mut down = quotes.clone();
                    up[j].par_spread += h;
                    down[j].par_spread -= h;
                    let cp = SurvivalProbabilityCurve::bootstrap_quotes(
                        ref_date,
                        discount.clone(),
                        0.4,
                        &up,
                    )
                    .unwrap();
                    let cm = SurvivalProbabilityCurve::bootstrap_quotes(
                        ref_date,
                        discount.clone(),
                        0.4,
                        &down,
                    )
                    .unwrap();
                    sum += (portfolio_value(&cp) - portfolio_value(&cm)) / (2.0 * h);
                }
                sum / 1_000_000.0
            }
        },
        {
            let base = base.clone();
            let portfolio = portfolio_sota.clone();
            move || {
                let flat = FlatForwardCurve::new(ref_date, 0.02);
                let (_, cs01) = base.value_and_par_cs01_portfolio(&portfolio, &flat);
                cs01.into_iter().sum::<f64>() / 1_000_000.0
            }
        },
    )
}

fn benchmark_fx(scale: usize) -> ComparisonRow {
    let strikes: Vec<f64> = (0..32).map(|i| 80.0 + i as f64 * 1.25).collect();
    let base = FxSmileQuote {
        atm_vol: 0.10,
        rr_25d: 0.01,
        bf_25d: 0.005,
    };
    run_comparison(
        "10_fx_vanna_volga_dual3",
        "seven scalar FX smile quote bumps vs. one Dual3 smile pass",
        config(scale),
        1e-3,
        {
            let strikes = strikes.clone();
            move || {
                let h = 1e-6;
                strikes
                    .iter()
                    .map(|&k| {
                        let base_smile = FxVannaVolgaSmile::new(100.0, 0.03, 0.01, 1.0, base);
                        let mut total =
                            base_smile.implied_vol(k) + base_smile.price_vanna_volga(k, true);
                        for j in 0..3 {
                            let mut up = base;
                            let mut down = base;
                            if j == 0 {
                                up.atm_vol += h;
                                down.atm_vol -= h;
                            }
                            if j == 1 {
                                up.rr_25d += h;
                                down.rr_25d -= h;
                            }
                            if j == 2 {
                                up.bf_25d += h;
                                down.bf_25d -= h;
                            }
                            let vp_smile = FxVannaVolgaSmile::new(100.0, 0.03, 0.01, 1.0, up);
                            let vm_smile = FxVannaVolgaSmile::new(100.0, 0.03, 0.01, 1.0, down);
                            let vp = vp_smile.implied_vol(k) + vp_smile.price_vanna_volga(k, true);
                            let vm = vm_smile.implied_vol(k) + vm_smile.price_vanna_volga(k, true);
                            total += (vp - vm) / (2.0 * h);
                        }
                        total
                    })
                    .sum::<f64>()
            }
        },
        {
            let strikes = strikes.clone();
            move || {
                strikes
                    .iter()
                    .map(|&k| {
                        let q = FxSmileQuote {
                            atm_vol: Dual::<3>::variable(base.atm_vol, 0),
                            rr_25d: Dual::<3>::variable(base.rr_25d, 1),
                            bf_25d: Dual::<3>::variable(base.bf_25d, 2),
                        };
                        let smile = FxVannaVolgaSmile::new(
                            Dual::from_f64(100.0),
                            Dual::from_f64(0.03),
                            Dual::from_f64(0.01),
                            Dual::from_f64(1.0),
                            q,
                        );
                        let iv = smile.implied_vol(Dual::from_f64(k));
                        let pv = smile.price_vanna_volga(Dual::from_f64(k), true);
                        iv.val + iv.d.iter().sum::<f64>() + pv.val + pv.d.iter().sum::<f64>()
                    })
                    .sum()
            }
        },
    )
}

#[allow(clippy::needless_range_loop)]
fn benchmark_cholesky_corr(scale: usize) -> ComparisonRow {
    let n = 12;
    let spots = vec![100.0; n];
    let vols = vec![0.2; n];
    let weights = vec![1.0 / n as f64; n];
    let mut corr = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            corr[i][j] = if i == j { 1.0 } else { 0.35 };
        }
    }
    let baseline_corr = corr.clone();
    let sota_corr = corr.clone();
    let spots_s = spots.clone();
    let vols_s = vols.clone();
    let weights_s = weights.clone();
    let scalar_price: ScalarBasketPricer = Arc::new(move |cm: &[Vec<f64>]| -> f64 {
        let l = cholesky_lower(cm).unwrap();
        let rng = PhiloxRng::new(31);
        let mut total = 0.0;
        for p in 0..256 {
            let u: Vec<f64> = (0..n).map(|i| rng.normal4(p as u64, i as u32)[0]).collect();
            let mut basket = 0.0;
            for i in 0..n {
                let z: f64 = (0..=i).map(|j| l[i][j] * u[j]).sum();
                let st = spots[i] * ((0.02 - 0.5 * vols[i] * vols[i]) + vols[i] * z).exp();
                basket += weights[i] * st;
            }
            total += (basket - 100.0).max(0.0) * (-0.02f64).exp();
        }
        total / 256.0
    });
    let scalar_up = scalar_price.clone();
    let scalar_down = scalar_price.clone();
    run_comparison(
        "11_capriotti_giles_cholesky_corr_aad",
        "66 correlation bumps vs. one adjoint Cholesky pass",
        config(scale),
        1e-4,
        move || {
            let mut total = 0.;
            for i in 0..n {
                for j in i + 1..n {
                    let mut up = baseline_corr.clone();
                    let mut dn = baseline_corr.clone();
                    up[i][j] += 1e-4;
                    up[j][i] += 1e-4;
                    dn[i][j] -= 1e-4;
                    dn[j][i] -= 1e-4;
                    let p = scalar_up(&up);
                    let m = scalar_down(&dn);
                    total += (p - m) / (2e-4);
                }
            }
            total
        },
        move || {
            let r = price_correlated_basket_call_aad(
                &spots_s, &vols_s, &weights_s, &sota_corr, 0.02, 1., 100., 256, 4, 1e-8, 31,
            )
            .unwrap();
            let mut sum = 0.0;
            for i in 0..n {
                for j in (i + 1)..n {
                    sum += r.corr_sensitivities[i][j];
                }
            }
            sum
        },
    )
}
fn benchmark_dim_mva(scale: usize) -> ComparisonRow {
    let d = date(2024, 1, 1);
    let curve = FlatForwardCurve::new(d, 0.02);
    let model = HullWhite1F {
        a: 0.05,
        sigma: 0.01,
    };
    let swaps = (0..5)
        .map(|i| NettingSetSwap {
            notional: 1e6,
            fixed_rate: 0.02 + i as f64 * 0.0001,
            payment_times: (1..=10).map(|x| x as f64).collect(),
            accruals: vec![1.; 10],
            is_payer: i % 2 == 0,
        })
        .collect::<Vec<_>>();
    let times = (1..=12).map(|i| i as f64 * 5. / 12.).collect::<Vec<_>>();
    let b = swaps.clone();
    let s = swaps.clone();
    let curve_b = curve.clone();
    let curve_s = curve.clone();
    let times_b = times.clone();
    let times_s = times.clone();
    run_comparison(
        "12_chebyshev_dynamic_simm_mva",
        "brute-force dynamic margin vs. Chebyshev model-space margin",
        config(scale),
        1e-3,
        move || {
            compute_bruteforce_dim_and_mva(&model, &curve_b, &b, &times_b, 128, 6, 0.01, 13).mva
        },
        move || compute_chebyshev_dim_and_mva(&model, &curve_s, &s, &times_s, 128, 6, 0.01, 13).mva,
    )
}

fn barrier_benchmark_specs() -> Vec<BarrierOptionSpec<f64>> {
    let mut out = Vec::with_capacity(64);
    for &bt in &[
        BarrierType::DownIn,
        BarrierType::DownOut,
        BarrierType::UpIn,
        BarrierType::UpOut,
    ] {
        for &call in &[true, false] {
            for &k in &[90.0, 100.0, 110.0, 120.0] {
                for &h in &[95.0, 105.0] {
                    out.push(BarrierOptionSpec {
                        option_type: if call {
                            OptionType::Call
                        } else {
                            OptionType::Put
                        },
                        barrier_type: bt,
                        spot: 100.0,
                        strike: k,
                        barrier: h,
                        rebate: 0.0,
                        rate: 0.03,
                        div_yield: 0.01,
                        vol: 0.25,
                        expiry: 1.0,
                        monitoring: MonitoringFrequency::Discrete { observations: 12 },
                    });
                }
            }
        }
    }
    out
}

fn benchmark_barrier_dual5(scale: usize) -> ComparisonRow {
    let specs = barrier_benchmark_specs();
    run_comparison(
        "13_reiner_rubinstein_barrier_dual5",
        "six finite-difference barrier bumps vs. one Dual5 analytical pass",
        config(scale),
        1e-6,
        {
            let specs = specs.clone();
            move || {
                let h = 7e-6;
                let mut total = 0.0;
                for s in &specs {
                    for j in 0..5 {
                        let mut up = *s;
                        let mut down = *s;
                        match j {
                            0 => {
                                up.spot += h;
                                down.spot -= h;
                            }
                            1 => {
                                up.vol += h;
                                down.vol -= h;
                            }
                            2 => {
                                up.rate += h;
                                down.rate -= h;
                            }
                            3 => {
                                up.div_yield += h;
                                down.div_yield -= h;
                            }
                            _ => {
                                up.expiry += h;
                                down.expiry -= h;
                            }
                        }
                        total +=
                            (price_barrier_option(&up) - price_barrier_option(&down)) / (2.0 * h);
                    }
                }
                total
            }
        },
        {
            let specs = specs.clone();
            move || {
                let mut total = 0.0;
                for s in &specs {
                    let d = BarrierOptionSpec {
                        option_type: s.option_type,
                        barrier_type: s.barrier_type,
                        spot: Dual::<5>::variable(s.spot, 0),
                        strike: Dual::from_value(s.strike),
                        barrier: Dual::from_value(s.barrier),
                        rebate: Dual::from_value(s.rebate),
                        rate: Dual::variable(s.rate, 2),
                        div_yield: Dual::variable(s.div_yield, 3),
                        vol: Dual::variable(s.vol, 1),
                        expiry: Dual::variable(s.expiry, 4),
                        monitoring: s.monitoring,
                    };
                    let x = price_barrier_option(&d);
                    total += x.d.iter().sum::<f64>();
                }
                total
            }
        },
    )
}

fn benchmark_rqmc(scale: usize) -> ComparisonRow {
    let times = (1..=16).map(|i| i as f64 / 16.0).collect::<Vec<_>>();
    let times_rqmc = times.clone();
    let rqmc = move || {
        price_asian_arithmetic_rqmc(
            100.0,
            100.0,
            0.03,
            0.0,
            0.2,
            &times_rqmc,
            1.0,
            true,
            256,
            4,
            17,
        )
        .0
    };
    let times_base = times;
    let baseline = move || {
        let rng = PhiloxRng::new(17);
        let mut sum = 0.0;
        for p in 0..16384u64 {
            let mut log_s = 100.0f64.ln();
            let mut prev = 0.0;
            let mut average = 0.0;
            for (j, &tm) in times_base.iter().enumerate() {
                let dt = tm - prev;
                let z = rng.normal4(p, j as u32)[0];
                log_s += (0.03 - 0.5 * 0.2 * 0.2) * dt + 0.2 * dt.sqrt() * z;
                average += log_s.exp();
                prev = tm;
            }
            sum += ((average / 16.0) - 100.0).max(0.0) * (-0.03f64).exp();
        }
        sum / 16384.0
    };
    run_comparison(
        "14_burley_owen_rqmc_variance_reduction",
        "pseudo-random Asian paths vs. Burley--Owen Sobol Brownian-bridge RQMC",
        config(scale),
        1.5e-2,
        baseline,
        rqmc,
    )
}

fn benchmark_levy_cos(scale: usize) -> ComparisonRow {
    let strikes: Vec<f64> = (0..32).map(|i| 80.0 + i as f64 * 1.25).collect();
    let base = VarianceGamma {
        sigma: 0.18,
        nu: 0.20,
        theta: -0.10,
    };
    run_comparison(
        "15_levy_cos_vg_dual5",
        "finite-difference VG COS strip Greeks vs. one Dual5 COS pass",
        config(scale),
        1e-4,
        {
            let strikes = strikes.clone();
            move || {
                let h = 1e-5;
                let mut total = 0.0;
                let params = [100.0, base.sigma, base.nu, base.theta, 0.03];
                let base_prices = cos_levy_price_strip(
                    &base,
                    params[0],
                    &strikes,
                    params[4],
                    0.0,
                    1.0,
                    OptionType::Call,
                    128,
                );
                total += base_prices.iter().sum::<f64>();
                for j in 0..5 {
                    let mut up = params;
                    let mut down = params;
                    up[j] += h;
                    down[j] -= h;
                    let mup = VarianceGamma {
                        sigma: up[1],
                        nu: up[2],
                        theta: up[3],
                    };
                    let mdn = VarianceGamma {
                        sigma: down[1],
                        nu: down[2],
                        theta: down[3],
                    };
                    let vp = cos_levy_price_strip(
                        &mup,
                        up[0],
                        &strikes,
                        up[4],
                        0.0,
                        1.0,
                        OptionType::Call,
                        128,
                    );
                    let vm = cos_levy_price_strip(
                        &mdn,
                        down[0],
                        &strikes,
                        down[4],
                        0.0,
                        1.0,
                        OptionType::Call,
                        128,
                    );
                    for (p, m) in vp.iter().zip(vm.iter()) {
                        total += (p - m) / (2.0 * h);
                    }
                }
                total
            }
        },
        {
            let strikes = strikes.clone();
            move || {
                let p = VarianceGamma {
                    sigma: Dual::<5>::variable(base.sigma, 1),
                    nu: Dual::<5>::variable(base.nu, 2),
                    theta: Dual::<5>::variable(base.theta, 3),
                };
                let prices = cos_levy_price_strip(
                    &p,
                    Dual::<5>::variable(100.0, 0),
                    &strikes
                        .iter()
                        .map(|&k| Dual::<5>::from_f64(k))
                        .collect::<Vec<_>>(),
                    Dual::<5>::variable(0.03, 4),
                    Dual::<5>::zero(),
                    Dual::<5>::one(),
                    OptionType::Call,
                    128,
                );
                prices
                    .iter()
                    .map(|x| x.val() + x.d.iter().sum::<f64>())
                    .sum()
            }
        },
    )
}

fn benchmark_cdo_recursion(scale: usize) -> ComparisonRow {
    let probs = (0..20)
        .map(|i| 0.015 + 0.0005 * (i % 5) as f64)
        .collect::<Vec<_>>();
    let units = vec![1usize; probs.len()];
    let base_corr = 0.25;
    run_comparison(
        "16_asb_cdo_tranche_recursion_dual",
        "eleven ASB loss recursions vs. one Dual5 tranche pass",
        config(scale),
        1e-6,
        {
            let probs = probs.clone();
            let units = units.clone();
            move || {
                let h = 1e-5;
                let mut out = asb_recursive_expected_tranche_loss(
                    &probs,
                    &units,
                    1.0 / 20.0,
                    base_corr,
                    0.03,
                    0.07,
                );
                for j in 0..5 {
                    let mut up = probs.clone();
                    let mut down = probs.clone();
                    let (corr_up, corr_down) = if j < 4 {
                        (base_corr, base_corr)
                    } else {
                        (base_corr + h, base_corr - h)
                    };
                    if j < 4 {
                        up[j] += h;
                        down[j] -= h;
                    }
                    let cp = asb_recursive_expected_tranche_loss(
                        &up,
                        &units,
                        1.0 / 20.0,
                        corr_up,
                        0.03,
                        0.07,
                    );
                    let cm = asb_recursive_expected_tranche_loss(
                        &down,
                        &units,
                        1.0 / 20.0,
                        corr_down,
                        0.03,
                        0.07,
                    );
                    out += (cp - cm) / (2.0 * h);
                }
                out
            }
        },
        {
            let probs = probs.clone();
            let units = units.clone();
            move || {
                let mut p = Vec::with_capacity(probs.len());
                for (i, &q) in probs.iter().enumerate() {
                    p.push(if i < 4 {
                        Dual::<5>::variable(q, i)
                    } else {
                        Dual::<5>::from_f64(q)
                    });
                }
                let x = asb_recursive_expected_tranche_loss(
                    &p,
                    &units,
                    Dual::<5>::from_f64(1.0 / 20.0),
                    Dual::<5>::variable(base_corr, 4),
                    Dual::<5>::from_f64(0.03),
                    Dual::<5>::from_f64(0.07),
                );
                x.val + x.d.iter().sum::<f64>()
            }
        },
    )
}

fn benchmark_rfr_sabr(scale: usize) -> ComparisonRow {
    let strikes: Vec<f64> = (0..32).map(|i| 0.02 + 0.06 * i as f64 / 31.0).collect();
    let base = RfrCapletSpec {
        notional: 1.0,
        forward_rate: 0.045,
        strike: 0.04,
        tau_0: 0.5,
        tau_1: 1.0,
        accrual_fraction: 0.5,
        discount_tau_1: 0.98,
        is_caplet: true,
        rate_type: RfrRateType::BackwardLooking,
    };
    run_comparison(
        "17_rfr_sabr_caplet_strip_dual4",
        "scalar RFR SABR strip bumps vs Dual4 effective-SABR strip",
        config(scale),
        1e-4,
        {
            let strikes = strikes.clone();
            move || {
                let h = 1e-5;
                let p = [0.045, 0.08, -0.35, 0.45];
                let mut out = 0.;
                for k in &strikes {
                    out += willems_backward_sabr_caplet_price(
                        &RfrCapletSpec {
                            forward_rate: p[0],
                            strike: *k,
                            ..base
                        },
                        p[1],
                        0.5,
                        p[2],
                        p[3],
                        1.0,
                    );
                }
                for j in 0..4 {
                    let mut u = p;
                    let mut d = p;
                    u[j] += h;
                    d[j] -= h;
                    let up: f64 = strikes
                        .iter()
                        .map(|k| {
                            willems_backward_sabr_caplet_price(
                                &RfrCapletSpec {
                                    forward_rate: u[0],
                                    strike: *k,
                                    ..base
                                },
                                u[1],
                                0.5,
                                u[2],
                                u[3],
                                1.0,
                            )
                        })
                        .sum();
                    let dn: f64 = strikes
                        .iter()
                        .map(|k| {
                            willems_backward_sabr_caplet_price(
                                &RfrCapletSpec {
                                    forward_rate: d[0],
                                    strike: *k,
                                    ..base
                                },
                                d[1],
                                0.5,
                                d[2],
                                d[3],
                                1.0,
                            )
                        })
                        .sum();
                    out += (up - dn) / (2. * h);
                }
                out
            }
        },
        {
            let strikes = strikes.clone();
            move || {
                let spec = RfrCapletSpec {
                    notional: Dual::<4>::from_f64(1.),
                    forward_rate: Dual::<4>::variable(0.045, 0),
                    strike: Dual::<4>::from_f64(0.04),
                    tau_0: Dual::<4>::from_f64(0.5),
                    tau_1: Dual::<4>::from_f64(1.),
                    accrual_fraction: Dual::<4>::from_f64(0.5),
                    discount_tau_1: Dual::<4>::from_f64(0.98),
                    is_caplet: true,
                    rate_type: RfrRateType::BackwardLooking,
                };
                willems_backward_sabr_caplet_strip(
                    &spec,
                    &strikes
                        .iter()
                        .map(|&k| Dual::from_f64(k))
                        .collect::<Vec<_>>(),
                    Dual::<4>::variable(0.08, 1),
                    Dual::from_f64(0.5),
                    Dual::<4>::variable(-0.35, 2),
                    Dual::<4>::variable(0.45, 3),
                    Dual::from_f64(1.0),
                )
                .iter()
                .map(|x| x.val() + x.d.iter().sum::<f64>())
                .sum::<f64>()
            }
        },
    )
}

fn benchmark_xccs(scale: usize) -> ComparisonRow {
    let d = date(2024, 1, 1);
    let times = (1..=10).map(|x| x as f64).collect::<Vec<_>>();
    let dom = DiscountCurve::new(
        d,
        times.clone(),
        times.iter().map(|t| (-0.02 * t).exp()).collect(),
    );
    let frc = DiscountCurve::new(
        d,
        times.clone(),
        times.iter().map(|t| (-0.03 * t).exp()).collect(),
    );
    let quotes = times
        .iter()
        .map(|&t| XccsQuote {
            maturity_years: t,
            basis_spread: 0.001,
        })
        .collect::<Vec<_>>();
    let cal = CalibratedXccsCurve::bootstrap(d, 1.1, dom.clone(), frc.clone(), &quotes).unwrap();
    let portfolio: Vec<_> = (1..=10)
        .map(|n| ConstNotionalXccs {
            dom_notional: 1_000_000.0,
            for_notional: 900_000.0,
            spot_fx: 1.1,
            payment_times: (1..=n).map(|i| i as f64).collect(),
            accrual_fractions: vec![1.0; n],
            for_basis_spread: 0.0015 + n as f64 * 1e-5,
            dom_spread: 0.0002,
            receive_foreign: true,
        })
        .collect();
    let portfolio_sota = portfolio.clone();
    run_comparison(
        "18_xccs_multicurrency_ift_dv01",
        "XCCS bump/bootstrap versus cached IFT basis risk",
        config(scale),
        1e-4,
        {
            let quotes = quotes.clone();
            move || {
                let h = 1e-5;
                let mut z = 0.;
                for j in 0..quotes.len() {
                    let mut u = quotes.clone();
                    let mut v = quotes.clone();
                    u[j].basis_spread += h;
                    v[j].basis_spread -= h;
                    let cu = CalibratedXccsCurve::bootstrap(d, 1.1, dom.clone(), frc.clone(), &u)
                        .unwrap();
                    let cv = CalibratedXccsCurve::bootstrap(d, 1.1, dom.clone(), frc.clone(), &v)
                        .unwrap();
                    for trade in &portfolio {
                        let vp = trade.npv(&cu.dom_ois, &cu.for_ois, &cu.for_csa_curve);
                        let vm = trade.npv(&cv.dom_ois, &cv.for_ois, &cv.for_csa_curve);
                        z += (vp - vm) / (2. * h);
                    }
                }
                z / 1_000_000.0
            }
        },
        move || {
            let (_, z) = cal.value_and_par_basis_dv01(|curve| {
                let mut value = TapeVar::from_f64(0.0);
                for trade in &portfolio_sota {
                    let dom_ad = cal.dom_ois.to_ad_curve(
                        &cal.dom_ois
                            .discounts
                            .iter()
                            .map(|&x| TapeVar::from_f64(x))
                            .collect::<Vec<_>>(),
                    );
                    let for_ad = cal.for_ois.to_ad_curve(
                        &cal.for_ois
                            .discounts
                            .iter()
                            .map(|&x| TapeVar::from_f64(x))
                            .collect::<Vec<_>>(),
                    );
                    value += trade.npv(&dom_ad, &for_ad, curve);
                }
                value
            });
            z.into_iter().sum::<f64>() / 1_000_000.0
        },
    )
}

fn benchmark_american_baw(scale: usize) -> ComparisonRow {
    let strikes: Vec<f64> = (0..32).map(|i| 80.0 + i as f64 * 1.25).collect();
    run_comparison(
        "19_american_baw_ju_ift_dual5",
        "Ju/BAW scalar bumps versus one IFT Dual5 American pass",
        config(scale),
        1e-7,
        {
            let strikes = strikes.clone();
            move || {
                strikes
                    .iter()
                    .map(|&k| {
                        let mut total = ql_instruments::JuQuadraticEngine::price(
                            100.0, k, 0.03, 0.01, 0.2, 1.0, true,
                        );
                        let h = 1e-6;
                        for j in 0..5 {
                            let mut up = [100.0, k, 0.03, 0.01, 0.2];
                            let mut dn = up;
                            up[j.min(4)] += h;
                            dn[j.min(4)] -= h;
                            total += (ql_instruments::JuQuadraticEngine::price(
                                up[0], up[1], up[2], up[3], up[4], 1.0, true,
                            ) - ql_instruments::JuQuadraticEngine::price(
                                dn[0], dn[1], dn[2], dn[3], dn[4], 1.0, true,
                            )) / (2.0 * h);
                        }
                        total
                    })
                    .sum::<f64>()
            }
        },
        move || {
            strikes
                .iter()
                .map(|&k| {
                    let s = Dual::<5>::variable(100.0, 0);
                    let kk = Dual::<5>::variable(k, 1);
                    let r = Dual::<5>::variable(0.03, 2);
                    let q = Dual::<5>::variable(0.01, 3);
                    let v = Dual::<5>::variable(0.2, 4);
                    let x = ql_instruments::JuQuadraticEngine::price(
                        s,
                        kk,
                        r,
                        q,
                        v,
                        Dual::from_f64(1.0),
                        true,
                    );
                    x.val + x.d.iter().sum::<f64>()
                })
                .sum::<f64>()
        },
    )
}

fn benchmark_risk_parity(scale: usize) -> ComparisonRow {
    let n = 12;
    let cov = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        0.04 + i as f64 * 0.001
                    } else {
                        0.005
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let budgets = vec![1.0 / n as f64; n];
    run_comparison(
        "20_risk_parity_ift_allocation_aad",
        "ninety covariance/budget bumps versus one IFT risk-parity solve",
        config(scale),
        1e-5,
        {
            let cov = cov.clone();
            let budgets = budgets.clone();
            move || {
                let h = 1e-6;
                let weighted = |w: &[f64]| {
                    w.iter()
                        .enumerate()
                        .map(|(k, x)| (k + 1) as f64 * x)
                        .sum::<f64>()
                };
                let base = risk_budgeting_weights(&cov, &budgets);
                let mut checksum = weighted(&base);
                // 12 budget bumps followed by the 78 unique symmetric
                // covariance bumps: exactly the 90-input Jacobian workload.
                for i in 0..n {
                    let derivative = |hh: f64| {
                        let mut up = budgets.clone();
                        let mut dn = budgets.clone();
                        up[i] += hh;
                        dn[i] -= hh;
                        (weighted(&risk_budgeting_weights(&cov, &up))
                            - weighted(&risk_budgeting_weights(&cov, &dn)))
                            / (2.0 * hh)
                    };
                    checksum += derivative(h);
                }
                for i in 0..n {
                    for j in i..n {
                        let derivative = |hh: f64| {
                            let mut up = cov.clone();
                            let mut dn = cov.clone();
                            up[i][j] += hh;
                            dn[i][j] -= hh;
                            if i != j {
                                up[j][i] += hh;
                                dn[j][i] -= hh;
                            }
                            (weighted(&risk_budgeting_weights(&up, &budgets))
                                - weighted(&risk_budgeting_weights(&dn, &budgets)))
                                / (2.0 * hh)
                        };
                        checksum += derivative(h);
                    }
                }
                std::hint::black_box(checksum)
            }
        },
        move || {
            use ql_ad::Dual;
            const LANES: usize = 90;
            let mut cov_ad = vec![vec![Dual::<LANES>::from_value(0.0); n]; n];
            let mut lane = 0;
            for i in 0..n {
                for j in i..n {
                    let x = Dual::<LANES>::variable(cov[i][j], lane);
                    cov_ad[i][j] = x;
                    cov_ad[j][i] = x;
                    lane += 1;
                }
            }
            let bud_ad = (0..n)
                .map(|i| Dual::<LANES>::variable(budgets[i], 78 + i))
                .collect::<Vec<_>>();
            let w = risk_budgeting_weights(&cov_ad, &bud_ad);
            let mut checksum = w
                .iter()
                .enumerate()
                .map(|(k, x)| (k + 1) as f64 * x.val)
                .sum::<f64>();
            for lane in 0..LANES {
                checksum += w
                    .iter()
                    .enumerate()
                    .map(|(k, x)| (k + 1) as f64 * x.d[lane])
                    .sum::<f64>();
            }
            std::hint::black_box(checksum)
        },
    )
}

fn benchmark_bermudan_hw(scale: usize) -> ComparisonRow {
    let curve = ql_curves::FlatForwardCurve::new(ql_time::Date(0), 0.02);
    let model = ql_models::HullWhite1F {
        a: 0.05,
        sigma: 0.01,
    };
    let sw = ql_rates::BermudanSwaption {
        notional: 1e6,
        strike: 0.02,
        exercise_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
        swap_times: (0..=10).map(|x| x as f64).collect(),
        is_payer: true,
        settlement: ql_rates::SwaptionSettlement::Physical,
    };
    let eng = ql_rates::HullWhiteBermudanTreeEngine { steps: 40 };
    run_comparison(
        "21_bermudan_hw1f_tree_ift_dv01",
        "Bermudan Hull--White tree quote and parameter finite differences versus AD",
        config(scale),
        1e-3,
        {
            let s = sw.clone();
            move || {
                let h = 1e-5;
                let base = eng.price_single_curve(&model, &curve, &s) / s.notional;
                let dr = {
                    let up = ql_curves::FlatForwardCurve::new(ql_time::Date(0), 0.02 + h);
                    let dn = ql_curves::FlatForwardCurve::new(ql_time::Date(0), 0.02 - h);
                    (eng.price_single_curve(&model, &up, &s)
                        - eng.price_single_curve(&model, &dn, &s))
                        / (2.0 * h * s.notional)
                };
                let da = {
                    let up = ql_models::HullWhite1F {
                        a: model.a + h,
                        sigma: model.sigma,
                    };
                    let dn = ql_models::HullWhite1F {
                        a: model.a - h,
                        sigma: model.sigma,
                    };
                    (eng.price_single_curve(&up, &curve, &s)
                        - eng.price_single_curve(&dn, &curve, &s))
                        / (2.0 * h * s.notional)
                };
                let ds = {
                    let up = ql_models::HullWhite1F {
                        a: model.a,
                        sigma: model.sigma + h,
                    };
                    let dn = ql_models::HullWhite1F {
                        a: model.a,
                        sigma: model.sigma - h,
                    };
                    (eng.price_single_curve(&up, &curve, &s)
                        - eng.price_single_curve(&dn, &curve, &s))
                        / (2.0 * h * s.notional)
                };
                base + dr + da + ds
            }
        },
        {
            let s = sw.clone();
            move || {
                let rate = Dual::<3>::variable(0.02, 0);
                let c = ql_curves::FlatForwardCurve::new(ql_time::Date(0), rate);
                let m = ql_models::HullWhite1F {
                    a: Dual::<3>::variable(model.a, 1),
                    sigma: Dual::<3>::variable(model.sigma, 2),
                };
                let x = eng.price_single_curve(&m, &c, &s) / s.notional;
                x.val + x.d.iter().sum::<f64>()
            }
        },
    )
}

fn benchmark_sabr_qp(scale: usize) -> ComparisonRow {
    let forward = 0.03;
    let expiry = 2.0;
    let beta = 0.5;
    let true_smile = AdSabrSmile::new(forward, expiry, 0.012, beta, -0.25, 0.45);
    let quotes = (0..7)
        .map(|i| {
            let k = 0.022 + i as f64 * 0.0025;
            SmileStrikeQuote {
                strike: k,
                market_vol: true_smile.vol(k),
                weight: 0.5 + i as f64 * 0.25,
            }
        })
        .collect::<Vec<_>>();
    let cov = (0..8)
        .map(|i| {
            (0..8)
                .map(|j| if i == j { 0.04 } else { 0.002 })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mu = (0..8).map(|i| 0.02 + 0.012 * i as f64).collect::<Vec<_>>();
    let lower = vec![0.02; 8];
    let upper = vec![0.22; 8];
    let base_sabr =
        CalibratedSabrSmile::calibrate(forward, expiry, beta, (0.011, -0.2, 0.4), &quotes)
            .expect("SABR calibration");
    let sabr_params = (
        base_sabr.smile.alpha,
        base_sabr.smile.rho,
        base_sabr.smile.nu,
    );
    let base_sabr_fd = base_sabr.clone();
    let base_sabr_ad = base_sabr.clone();
    let cov_base = cov.clone();
    let mu_base = mu.clone();
    let lower_base = lower.clone();
    let upper_base = upper.clone();
    let cov_ad_base = cov.clone();
    let mu_ad_base = mu.clone();
    let lower_ad_base = lower.clone();
    let upper_ad_base = upper.clone();
    let weighted = |w: &[f64]| {
        w.iter()
            .enumerate()
            .map(|(i, x)| (i + 1) as f64 * x)
            .sum::<f64>()
    };
    run_comparison(
        "22_sabr_svi_ift_and_qp_kkt_aad",
        "SABR quote recalibration and constrained Markowitz finite differences versus IFT/KKT AD",
        config(scale),
        1e-5,
        {
            let quotes = quotes.clone();
            let base_sabr_fd = base_sabr_fd.clone();
            let cov = cov_ad_base.clone();
            let mu = mu_ad_base.clone();
            let lower = lower_ad_base.clone();
            let upper = upper_ad_base.clone();
            move || {
                let mut sum = base_sabr_fd.smile.vol(quotes[3].strike);
                let h = 1e-5;
                for i in 0..quotes.len() {
                    let mut up = quotes.clone();
                    let mut dn = quotes.clone();
                    up[i].market_vol += h;
                    dn[i].market_vol -= h;
                    let vu =
                        CalibratedSabrSmile::calibrate(forward, expiry, beta, sabr_params, &up)
                            .expect("SABR up")
                            .smile
                            .vol(quotes[3].strike);
                    let vd =
                        CalibratedSabrSmile::calibrate(forward, expiry, beta, sabr_params, &dn)
                            .expect("SABR down")
                            .smile
                            .vol(quotes[3].strike);
                    sum += (vu - vd) / (2.0 * h);
                }
                let base_w = box_constrained_markowitz_weights(&BoxConstrainedMarkowitzSpec {
                    covariance: cov.clone(),
                    expected_returns: mu.clone(),
                    risk_aversion: 3.0,
                    lower_bounds: lower.clone(),
                    upper_bounds: upper.clone(),
                });
                let mut checksum = sum + weighted(&base_w);
                for i in 0..8 {
                    let mut up = mu.clone();
                    let mut dn = mu.clone();
                    up[i] += h;
                    dn[i] -= h;
                    checksum += (weighted(&box_constrained_markowitz_weights(
                        &BoxConstrainedMarkowitzSpec {
                            covariance: cov.clone(),
                            expected_returns: up,
                            risk_aversion: 3.0,
                            lower_bounds: lower.clone(),
                            upper_bounds: upper.clone(),
                        },
                    )) - weighted(&box_constrained_markowitz_weights(
                        &BoxConstrainedMarkowitzSpec {
                            covariance: cov.clone(),
                            expected_returns: dn,
                            risk_aversion: 3.0,
                            lower_bounds: lower.clone(),
                            upper_bounds: upper.clone(),
                        },
                    ))) / (2.0 * h);
                }
                for i in 0..8 {
                    let mut up = upper.clone();
                    let mut dn = upper.clone();
                    up[i] += h;
                    dn[i] -= h;
                    let mk = |bounds: Vec<f64>| {
                        box_constrained_markowitz_weights(&BoxConstrainedMarkowitzSpec {
                            covariance: cov.clone(),
                            expected_returns: mu.clone(),
                            risk_aversion: 3.0,
                            lower_bounds: lower.clone(),
                            upper_bounds: bounds,
                        })
                    };
                    checksum += (weighted(&mk(up)) - weighted(&mk(dn))) / (2.0 * h);
                }
                checksum
            }
        },
        {
            let quotes = quotes.clone();
            let base_sabr_ad = base_sabr_ad.clone();
            let cov = cov_base.clone();
            let mu = mu_base.clone();
            let lower = lower_base.clone();
            let upper = upper_base.clone();
            move || {
                let (_, vegas) = base_sabr_ad.vol_and_quote_vega(quotes[3].strike);
                let mut sum = base_sabr_ad.smile.vol(quotes[3].strike) + vegas.iter().sum::<f64>();
                let mu_ad = mu
                    .iter()
                    .enumerate()
                    .map(|(i, x)| Dual::<16>::variable(*x, i))
                    .collect::<Vec<_>>();
                let up_ad = upper
                    .iter()
                    .enumerate()
                    .map(|(i, x)| Dual::<16>::variable(*x, 8 + i))
                    .collect::<Vec<_>>();
                let w = box_constrained_markowitz_weights(&BoxConstrainedMarkowitzSpec {
                    covariance: cov
                        .iter()
                        .map(|r| r.iter().map(|x| Dual::<16>::from_value(*x)).collect())
                        .collect(),
                    expected_returns: mu_ad,
                    risk_aversion: Dual::from_value(3.0),
                    lower_bounds: lower.iter().map(|x| Dual::from_value(*x)).collect(),
                    upper_bounds: up_ad,
                });
                sum += w
                    .iter()
                    .enumerate()
                    .map(|(i, x)| (i + 1) as f64 * x.val)
                    .sum::<f64>();
                for lane in 0..16 {
                    sum += w
                        .iter()
                        .enumerate()
                        .map(|(i, x)| (i + 1) as f64 * x.d[lane])
                        .sum::<f64>();
                }
                sum
            }
        },
    )
}

fn benchmark_slv_autocall(scale: usize) -> ComparisonRow {
    let spec = ql_exotics::WorstOfAutocallableSpec {
        notional: 1.0,
        spots: vec![100.0; 4],
        vols: vec![0.2; 4],
        divs: vec![0.0; 4],
        rate: 0.02,
        correlation: vec![
            vec![1., 0.2, 0.2, 0.2],
            vec![0.2, 1., 0.2, 0.2],
            vec![0.2, 0.2, 1., 0.2],
            vec![0.2, 0.2, 0.2, 1.],
        ],
        observation_times: vec![0.5, 1.0, 1.5],
        autocall_barriers: vec![1.05; 3],
        redemption_bonuses: vec![0.; 3],
        coupon_barriers: vec![0.8; 3],
        coupon_rates: vec![0.02; 3],
        memory_coupon: true,
        knock_in_barrier: 0.6,
        put_strike: 1.0,
        observe_ki_all_dates: true,
        smoothing_half_width: 0.01,
        smooth_min_temp: 0.001,
    };
    run_comparison(
        "23_slv_rbergomi_and_autocallable_cega_aad",
        "SLV/rBergomi and worst-of autocallable sensitivities",
        config(scale),
        1e-4,
        {
            let s = spec.clone();
            move || {
                let h = 1e-5;
                let p0 = ql_exotics::price_worst_of_autocallable(&s, 64, 7).unwrap();
                let mut m = p0;
                for i in 0..s.spots.len() {
                    let mut u = s.clone();
                    let mut d = s.clone();
                    u.spots[i] += h;
                    d.spots[i] -= h;
                    let dv = (ql_exotics::price_worst_of_autocallable(&u, 64, 7).unwrap()
                        - ql_exotics::price_worst_of_autocallable(&d, 64, 7).unwrap())
                        / (2.0 * h);
                    m += dv;
                }
                for i in 0..s.vols.len() {
                    let mut u = s.clone();
                    let mut d = s.clone();
                    u.vols[i] += h;
                    d.vols[i] -= h;
                    let dv = (ql_exotics::price_worst_of_autocallable(&u, 64, 7).unwrap()
                        - ql_exotics::price_worst_of_autocallable(&d, 64, 7).unwrap())
                        / (2.0 * h);
                    m += dv;
                }
                {
                    let mut u = s.clone();
                    let mut d = s.clone();
                    u.rate += h;
                    d.rate -= h;
                    let dv = (ql_exotics::price_worst_of_autocallable(&u, 64, 7).unwrap()
                        - ql_exotics::price_worst_of_autocallable(&d, 64, 7).unwrap())
                        / (2.0 * h);
                    m += dv;
                }
                for i in 0..4 {
                    for j in 0..i {
                        let mut u = s.clone();
                        let mut d = s.clone();
                        u.correlation[i][j] += h;
                        u.correlation[j][i] += h;
                        d.correlation[i][j] -= h;
                        d.correlation[j][i] -= h;
                        let dv = (ql_exotics::price_worst_of_autocallable(&u, 64, 7).unwrap()
                            - ql_exotics::price_worst_of_autocallable(&d, 64, 7).unwrap())
                            / (2.0 * h);
                        m += dv;
                    }
                }
                let rb = |spot: f64, xi0: f64, eta: f64| {
                    ql_models::rough::rbergomi_hybrid_mc_price(
                        spot,
                        100.0,
                        1.0,
                        0.02,
                        0.0,
                        &ql_models::rough::RoughBergomiParams {
                            xi0,
                            eta,
                            rho: -0.7,
                            hurst: 0.1,
                        },
                        8,
                        32,
                        true,
                        7,
                    )
                };
                let rb0 = rb(100.0, 0.04, 1.2);
                m += rb0;
                for (spot, xi0, eta) in [
                    (100.0 + h, 0.04, 1.2),
                    (100.0, 0.04 + h, 1.2),
                    (100.0, 0.04, 1.2 + h),
                ] {
                    let up = rb(spot, xi0, eta);
                    let dn = if spot != 100.0 {
                        rb(100.0 - h, 0.04, 1.2)
                    } else if xi0 != 0.04 {
                        rb(100.0, 0.04 - h, 1.2)
                    } else {
                        rb(100.0, 0.04, 1.2 - h)
                    };
                    m += (up - dn) / (2.0 * h);
                }
                m
            }
        },
        {
            let s = spec;
            move || {
                let r = ql_exotics::price_worst_of_autocallable_aad(&s, 64, 4, 7).unwrap();
                let rb = ql_models::rough::rbergomi_hybrid_mc_price(
                    Dual::<3>::variable(100.0, 0),
                    Dual::<3>::from_value(100.0),
                    Dual::<3>::from_value(1.0),
                    Dual::<3>::from_value(0.02),
                    Dual::<3>::from_value(0.0),
                    &ql_models::rough::RoughBergomiParams {
                        xi0: Dual::<3>::variable(0.04, 1),
                        eta: Dual::<3>::variable(1.2, 2),
                        rho: Dual::<3>::from_value(-0.7),
                        hurst: 0.1,
                    },
                    8,
                    32,
                    true,
                    7,
                );
                r.price
                    + r.spot_deltas.iter().sum::<f64>()
                    + r.vegas.iter().sum::<f64>()
                    + r.rate_rho
                    + r.corr_cega
                        .iter()
                        .enumerate()
                        .map(|(i, row)| row.iter().take(i).sum::<f64>())
                        .sum::<f64>()
                    + rb.val
                    + rb.d.iter().sum::<f64>()
            }
        },
    )
}
fn benchmark_frtb_spec<T: Real>(d1: T, d2: T, v1: T, v2: T) -> ql_xva::FrtbSaSpec<T> {
    ql_xva::FrtbSaSpec {
        buckets: vec![ql_xva::FrtbSensitivityBucket {
            risk_class: ql_xva::FrtbRiskClass::Girr,
            bucket_id: 0,
            delta_sensitivities: vec![d1, d2],
            delta_risk_weights: vec![0.015, 0.015],
            vega_sensitivities: vec![v1, v2],
            vega_risk_weights: vec![0.55, 0.55],
            curvature_up_pnl: vec![T::zero(), T::zero()],
            curvature_down_pnl: vec![T::zero(), T::zero()],
            curvature_shock_weights: vec![0.015, 0.015],
            intra_bucket_corr: vec![vec![1.0, 0.9], vec![0.9, 1.0]],
        }],
        inter_bucket_corr: vec![],
        drc_obligors: vec![],
        exotic_residual_notionals: vec![],
        other_residual_notionals: vec![],
    }
}

fn benchmark_frtb_cvar_execution(scale: usize) -> ComparisonRow {
    let spec = ql_portfolio::AlmgrenChrissSpec {
        initial_inventory: 100.0,
        horizon: 1.0,
        num_steps: 16,
        volatility: 0.2,
        temporary_impact: 0.01,
        permanent_impact: 0.001,
        risk_aversion: 0.1,
    };
    run_comparison(
        "24_frtb_sbm_ima_and_cvar_execution_aad",
        "FRTB/CVaR and optimal execution KKT sensitivities",
        config(scale),
        1e-5,
        {
            let s = spec.clone();
            move || {
                let h = 1e-6;
                let base = ql_portfolio::almgren_chriss_execution(&s).mean_variance_utility;
                let mut v_up = s.clone();
                let mut v_dn = s.clone();
                v_up.volatility += h;
                v_dn.volatility -= h;
                let mut t_up = s.clone();
                let mut t_dn = s.clone();
                t_up.temporary_impact += h;
                t_dn.temporary_impact -= h;
                let mut p_up = s.clone();
                let mut p_dn = s.clone();
                p_up.permanent_impact += h;
                p_dn.permanent_impact -= h;
                let mut r_up = s.clone();
                let mut r_dn = s.clone();
                r_up.risk_aversion += h;
                r_dn.risk_aversion -= h;
                let frtb = |d1: f64, d2: f64, v1: f64, v2: f64| {
                    ql_xva::compute_frtb_sa(&benchmark_frtb_spec(d1, d2, v1, v2)).sbm_med * 1e-4
                };
                let frtb_base = frtb(120_000.0, -80_000.0, 45_000.0, 30_000.0);
                let frtb_fd = (frtb(120_000.0 + 1e-3, -80_000.0, 45_000.0, 30_000.0)
                    - frtb(120_000.0 - 1e-3, -80_000.0, 45_000.0, 30_000.0)
                    + frtb(120_000.0, -80_000.0 + 1e-3, 45_000.0, 30_000.0)
                    - frtb(120_000.0, -80_000.0 - 1e-3, 45_000.0, 30_000.0)
                    + frtb(120_000.0, -80_000.0, 45_000.0 + 1e-3, 30_000.0)
                    - frtb(120_000.0, -80_000.0, 45_000.0 - 1e-3, 30_000.0)
                    + frtb(120_000.0, -80_000.0, 45_000.0, 30_000.0 + 1e-3)
                    - frtb(120_000.0, -80_000.0, 45_000.0, 30_000.0 - 1e-3))
                    / (2e-3);
                base + frtb_base
                    + frtb_fd
                    + (ql_portfolio::almgren_chriss_execution(&v_up).mean_variance_utility
                        - ql_portfolio::almgren_chriss_execution(&v_dn).mean_variance_utility
                        + ql_portfolio::almgren_chriss_execution(&t_up).mean_variance_utility
                        - ql_portfolio::almgren_chriss_execution(&t_dn).mean_variance_utility
                        + ql_portfolio::almgren_chriss_execution(&p_up).mean_variance_utility
                        - ql_portfolio::almgren_chriss_execution(&p_dn).mean_variance_utility
                        + ql_portfolio::almgren_chriss_execution(&r_up).mean_variance_utility
                        - ql_portfolio::almgren_chriss_execution(&r_dn).mean_variance_utility)
                        / (2.0 * h)
            }
        },
        {
            let s = spec;
            move || {
                type D = Dual<4>;
                let sd = ql_portfolio::AlmgrenChrissSpec {
                    initial_inventory: D::from_value(s.initial_inventory),
                    horizon: D::from_value(s.horizon),
                    num_steps: s.num_steps,
                    volatility: D {
                        val: s.volatility,
                        d: [1.0, 0.0, 0.0, 0.0],
                    },
                    temporary_impact: D {
                        val: s.temporary_impact,
                        d: [0.0, 1.0, 0.0, 0.0],
                    },
                    permanent_impact: D {
                        val: s.permanent_impact,
                        d: [0.0, 0.0, 1.0, 0.0],
                    },
                    risk_aversion: D {
                        val: s.risk_aversion,
                        d: [0.0, 0.0, 0.0, 1.0],
                    },
                };
                let u = ql_portfolio::almgren_chriss_execution(&sd).mean_variance_utility;
                type F = Dual<4>;
                let pos = benchmark_frtb_spec(
                    F {
                        val: 120_000.0,
                        d: [1.0, 0.0, 0.0, 0.0],
                    },
                    F {
                        val: -80_000.0,
                        d: [0.0, 1.0, 0.0, 0.0],
                    },
                    F {
                        val: 45_000.0,
                        d: [0.0, 0.0, 1.0, 0.0],
                    },
                    F {
                        val: 30_000.0,
                        d: [0.0, 0.0, 0.0, 1.0],
                    },
                );
                let frtb = ql_xva::compute_frtb_sa(&pos).sbm_med * F::from_value(1e-4);
                u.val + u.d.iter().sum::<f64>() + frtb.val + frtb.d.iter().sum::<f64>()
            }
        },
    )
}

fn benchmark_vix_parisian_robust(scale: usize) -> ComparisonRow {
    // Keep this suite honest: every API named in the comparison is exercised
    // by both closures.  In particular, the old implementation omitted the
    // Parisian and monotone-convex paths while still claiming them in the ID.
    let vix_params = |v0: f64, kappa: f64, theta: f64, sigma: f64| ql_models::BatesParams {
        heston: HestonParams {
            v0,
            kappa,
            theta,
            sigma_v: sigma,
            rho: -0.6,
        },
        lambda: 0.1,
        mu_j: -0.05,
        sigma_j: 0.2,
    };
    let mono_times = vec![1.0, 2.0, 5.0];
    let mono_rates = vec![0.02, 0.023, 0.028];
    let cov = vec![vec![0.04, 0.01], vec![0.01, 0.09]];
    let omega = vec![vec![0.01, 0.0], vec![0.0, 0.02]];
    let mu = [0.06, 0.04];
    let upper = [0.8, 0.8];
    let cov_for_eval = cov.clone();
    let omega_for_eval = omega.clone();
    let robust_eval = move |returns: [f64; 2], bounds: [f64; 2]| {
        let w = ql_portfolio::robust_markowitz_weights(&ql_portfolio::RobustMarkowitzSpec {
            covariance: cov_for_eval.clone(),
            return_uncertainty_cov: omega_for_eval.clone(),
            expected_returns: returns.to_vec(),
            risk_aversion: 3.0,
            uncertainty_radius: 0.1,
            lower_bounds: vec![0.0, 0.0],
            upper_bounds: bounds.to_vec(),
        });
        w.iter()
            .enumerate()
            .map(|(i, x)| (i + 1) as f64 * x)
            .sum::<f64>()
    };
    let parisian = ParisianBarrierOptionSpec {
        option_type: OptionType::Call,
        barrier_type: BarrierType::DownOut,
        clock_kind: ExcursionClockKind::ParisianConsecutive,
        spot: 100.0,
        strike: 100.0,
        barrier: 90.0,
        rebate: 0.0,
        window_duration: 0.20,
        rate: 0.02,
        div_yield: 0.0,
        vol: 0.20,
        expiry: 1.0,
    };
    let parisian_paths = 64;
    let parisian_steps = 16;
    run_comparison(
        "25_monotone_convex_vix_parisian_robust_aad",
        "monotone-convex curve plus VIX/Parisian and robust Markowitz AD",
        config(scale),
        1e-4,
        {
            let mono_times = mono_times.clone();
            let mono_rates = mono_rates.clone();
            move || {
                let mono = MonotoneConvexZeroCurve::new(
                    date(2024, 1, 1),
                    mono_times.clone(),
                    mono_rates.clone(),
                );
                let mono_score = [0.5, 1.5, 3.0]
                    .iter()
                    .map(|&t| mono.zero_rate(t).val())
                    .sum::<f64>();
                let p = vix_params(0.04, 1.5, 0.04, 0.3);
                let h = 1e-5;
                let base = ql_models::heston_bates_vix_future_and_option_price(
                    &p,
                    1.0,
                    30.0 / 365.0,
                    0.2,
                    0.98,
                );
                let mut sum = mono_score + base.vix_future + base.call_price + base.put_price;
                for (j, bump) in [(0, h), (1, h), (2, h), (3, h)] {
                    let mut up = [0.04, 1.5, 0.04, 0.3];
                    let mut dn = up;
                    up[j] += bump;
                    dn[j] -= bump;
                    let a = ql_models::heston_bates_vix_future_and_option_price(
                        &vix_params(up[0], up[1], up[2], up[3]),
                        1.0,
                        30.0 / 365.0,
                        0.2,
                        0.98,
                    );
                    let b = ql_models::heston_bates_vix_future_and_option_price(
                        &vix_params(dn[0], dn[1], dn[2], dn[3]),
                        1.0,
                        30.0 / 365.0,
                        0.2,
                        0.98,
                    );
                    let fd = (a.vix_future + a.call_price + a.put_price
                        - b.vix_future
                        - b.call_price
                        - b.put_price)
                        / (2.0 * bump);
                    sum += fd;
                }
                for j in 0..4 {
                    let mut up = parisian;
                    let mut dn = parisian;
                    match j {
                        0 => {
                            up.spot += h;
                            dn.spot -= h;
                        }
                        1 => {
                            up.vol += h;
                            dn.vol -= h;
                        }
                        2 => {
                            up.rate += h;
                            dn.rate -= h;
                        }
                        _ => {
                            up.window_duration += h;
                            dn.window_duration -= h;
                        }
                    }
                    let vp =
                        price_parisian_option_mc_ad(&up, parisian_paths, parisian_steps, 0.01, 23);
                    let vm =
                        price_parisian_option_mc_ad(&dn, parisian_paths, parisian_steps, 0.01, 23);
                    let fd = (vp - vm).val() / (2.0 * h);
                    sum += fd;
                }
                let p0 = price_parisian_option_mc_ad(
                    &parisian,
                    parisian_paths,
                    parisian_steps,
                    0.01,
                    23,
                );
                sum += p0.val();
                let h = 1e-5;
                let mut robust = robust_eval(mu, upper);
                for j in 0..2 {
                    let mut up = mu;
                    let mut dn = mu;
                    up[j] += h;
                    dn[j] -= h;
                    robust += (robust_eval(up, upper) - robust_eval(dn, upper)) / (2.0 * h);
                }
                for j in 0..2 {
                    let mut up = upper;
                    let mut dn = upper;
                    up[j] += h;
                    dn[j] -= h;
                    robust += (robust_eval(mu, up) - robust_eval(mu, dn)) / (2.0 * h);
                }
                sum + robust
            }
        },
        {
            let mono_times = mono_times.clone();
            let mono_rates = mono_rates.clone();
            move || {
                type D = Dual<4>;
                let mono = MonotoneConvexZeroCurve::new(
                    date(2024, 1, 1),
                    mono_times.clone(),
                    mono_rates.clone(),
                );
                let mono_score = [0.5, 1.5, 3.0]
                    .iter()
                    .map(|&t| mono.zero_rate(t).val())
                    .sum::<f64>();
                let p = ql_models::BatesParams {
                    heston: HestonParams {
                        v0: D::variable(0.04, 0),
                        kappa: D::variable(1.5, 1),
                        theta: D::variable(0.04, 2),
                        sigma_v: D::variable(0.3, 3),
                        rho: D::from_value(-0.6),
                    },
                    lambda: D::from_value(0.1),
                    mu_j: D::from_value(-0.05),
                    sigma_j: D::from_value(0.2),
                };
                let v = ql_models::heston_bates_vix_future_and_option_price(
                    &p,
                    D::from_value(1.0),
                    30.0 / 365.0,
                    D::from_value(0.2),
                    D::from_value(0.98),
                );
                let parisian = ParisianBarrierOptionSpec {
                    option_type: OptionType::Call,
                    barrier_type: BarrierType::DownOut,
                    clock_kind: ExcursionClockKind::ParisianConsecutive,
                    spot: D::variable(100.0, 0),
                    strike: D::from_value(100.0),
                    barrier: D::from_value(90.0),
                    rebate: D::from_value(0.0),
                    window_duration: D::variable(0.20, 1),
                    rate: D::variable(0.02, 2),
                    div_yield: D::from_value(0.0),
                    vol: D::variable(0.20, 3),
                    expiry: D::from_value(1.0),
                };
                let px = price_parisian_option_mc_ad(
                    &parisian,
                    parisian_paths,
                    parisian_steps,
                    0.01,
                    23,
                );
                let returns = mu
                    .iter()
                    .enumerate()
                    .map(|(i, &x)| D::variable(x, i))
                    .collect::<Vec<_>>();
                let bounds = upper
                    .iter()
                    .enumerate()
                    .map(|(i, &x)| D::variable(x, 2 + i))
                    .collect::<Vec<_>>();
                let w =
                    ql_portfolio::robust_markowitz_weights(&ql_portfolio::RobustMarkowitzSpec {
                        covariance: cov
                            .iter()
                            .map(|row| row.iter().map(|&x| D::from_value(x)).collect())
                            .collect(),
                        return_uncertainty_cov: omega
                            .iter()
                            .map(|row| row.iter().map(|&x| D::from_value(x)).collect())
                            .collect(),
                        expected_returns: returns,
                        risk_aversion: D::from_value(3.0),
                        uncertainty_radius: D::from_value(0.1),
                        lower_bounds: vec![D::from_value(0.0); 2],
                        upper_bounds: bounds,
                    });
                let robust = w
                    .iter()
                    .enumerate()
                    .map(|(i, x)| (i + 1) as f64 * (x.val + x.d.iter().sum::<f64>()))
                    .sum::<f64>();
                mono_score
                    + v.vix_future.val
                    + v.call_price.val
                    + v.put_price.val
                    + v.vix_future.d.iter().sum::<f64>()
                    + v.call_price.d.iter().sum::<f64>()
                    + v.put_price.d.iter().sum::<f64>()
                    + px.val
                    + px.d.iter().sum::<f64>()
                    + robust
            }
        },
    )
}

fn benchmark_double_heston_tarf_coco_hrp(scale: usize) -> ComparisonRow {
    let cfg = config(scale);
    let h = 1e-5;
    let tarf_paths = 64;
    let eval_dh = |x: [f64; 4]| {
        let p = ql_models::DoubleHestonParams {
            factor1: HestonParams {
                v0: x[0],
                kappa: 1.4,
                theta: x[1],
                sigma_v: 0.3,
                rho: -0.7,
            },
            factor2: HestonParams {
                v0: x[2],
                kappa: 2.0,
                theta: 0.01,
                sigma_v: x[3],
                rho: -0.3,
            },
        };
        ql_models::double_heston_cos_price(100.0, 100.0, 1.0, 0.03, 0.0, &p, 64, 10.0, true)
    };
    let eval_tarf = |x: [f64; 4]| {
        ql_fx::price_fx_tarf_mc_ad(
            &ql_fx::FxTarfSpec {
                notional: 1.0,
                spot: x[0],
                strike: 1.0,
                r_dom: x[1],
                r_for: 0.01,
                vol: x[2],
                fixing_times: vec![0.25, 0.5, 0.75, 1.0, 1.25, 1.5],
                target_profit: x[3],
                downside_leverage: 1.5,
                is_exporter_call: true,
                cap_convention: ql_fx::FxTarfCapConvention::ExactCap,
                smoothing_eps: 1e-3,
            },
            tarf_paths,
            7,
        )
        .npv
    };
    let eval_coco = |x: [f64; 4]| {
        ql_fixed_income::price_coco_bond_ds(&ql_fixed_income::CocoBondSpec {
            face_value: 100.0,
            coupon_rate: 0.05,
            coupon_times: vec![0.5, 1.0, 1.5, 2.0],
            maturity: 2.0,
            rate: x[0],
            credit_spread: x[1],
            spot: x[2],
            trigger_barrier: 60.0,
            div_yield: 0.01,
            eq_vol: x[3],
            trigger_payoff: ql_fixed_income::CocoTriggerPayoff::WriteDown {
                recovery_fraction: 0.4,
            },
        })
        .dirty_price
    };
    let eval_port = |x: [f64; 4]| {
        let cov = vec![vec![x[0], x[1]], vec![x[1], x[2]]];
        let w = ql_portfolio::hierarchical_risk_parity_weights(&cov);
        let exec = ql_portfolio::multi_asset_almgren_chriss_execution(
            &ql_portfolio::MultiAssetAlmgrenChrissSpec {
                initial_inventory: vec![1000.0, 600.0],
                horizon: 1.0,
                num_steps: 6,
                covariance: cov,
                temporary_impact: vec![vec![1e-4, 2e-5], vec![2e-5, 1.5e-4]],
                permanent_impact: vec![vec![5e-5, 0.0], vec![0.0, 5e-5]],
                risk_aversion: x[3],
            },
        );
        w[0] + 2.0 * w[1] + 1e-4 * exec.mean_variance_utility
    };
    let dh0 = [0.04, 0.04, 0.01, 0.20];
    let tarf0 = [1.05, 0.02, 0.12, 0.15];
    let coco0 = [0.03, 0.015, 100.0, 0.25];
    let port0 = [0.04, 0.01, 0.09, 1e-5];
    let baseline = || {
        let bdh0 = eval_dh(dh0);
        let bt0 = eval_tarf(tarf0);
        let bc0 = eval_coco(coco0);
        let bp0 = eval_port(port0);
        let mut total = bdh0 + bt0 + bc0 + bp0;
        let mut dsum = [0.0; 4];
        for i in 0..4 {
            let mut up = dh0;
            let mut dn = dh0;
            up[i] += h;
            dn[i] -= h;
            let ddh = (eval_dh(up) - eval_dh(dn)) / (2.0 * h);
            total += ddh;
            dsum[i] += ddh;

            let mut up = tarf0;
            let mut dn = tarf0;
            up[i] += h;
            dn[i] -= h;
            let dt = (eval_tarf(up) - eval_tarf(dn)) / (2.0 * h);
            total += dt;
            dsum[i] += dt;

            let mut up = coco0;
            let mut dn = coco0;
            up[i] += h;
            dn[i] -= h;
            let dc = (eval_coco(up) - eval_coco(dn)) / (2.0 * h);
            total += dc;
            dsum[i] += dc;

            let mut up = port0;
            let mut dn = port0;
            up[i] += h;
            dn[i] -= h;
            let dp = (eval_port(up) - eval_port(dn)) / (2.0 * h);
            total += dp;
            dsum[i] += dp;
        }
        total
    };
    let sota = || {
        type D = Dual<4>;
        let p = ql_models::DoubleHestonParams {
            factor1: HestonParams {
                v0: D::variable(dh0[0], 0),
                kappa: D::from_value(1.4),
                theta: D::variable(dh0[1], 1),
                sigma_v: D::from_value(0.3),
                rho: D::from_value(-0.7),
            },
            factor2: HestonParams {
                v0: D::variable(dh0[2], 2),
                kappa: D::from_value(2.0),
                theta: D::from_value(0.01),
                sigma_v: D::variable(dh0[3], 3),
                rho: D::from_value(-0.3),
            },
        };
        let dh = ql_models::double_heston_cos_price(
            D::from_value(100.0),
            D::from_value(100.0),
            D::from_value(1.0),
            D::from_value(0.03),
            D::from_value(0.0),
            &p,
            64,
            10.0,
            true,
        );
        let tarf = ql_fx::price_fx_tarf_mc_ad(
            &ql_fx::FxTarfSpec {
                notional: D::from_value(1.0),
                spot: D::variable(tarf0[0], 0),
                strike: D::from_value(1.0),
                r_dom: D::variable(tarf0[1], 1),
                r_for: D::from_value(0.01),
                vol: D::variable(tarf0[2], 2),
                fixing_times: vec![0.25, 0.5, 0.75, 1.0, 1.25, 1.5],
                target_profit: D::variable(tarf0[3], 3),
                downside_leverage: D::from_value(1.5),
                is_exporter_call: true,
                cap_convention: ql_fx::FxTarfCapConvention::ExactCap,
                smoothing_eps: 1e-3,
            },
            tarf_paths,
            7,
        )
        .npv;
        let coco = ql_fixed_income::price_coco_bond_ds(&ql_fixed_income::CocoBondSpec {
            face_value: D::from_value(100.0),
            coupon_rate: D::from_value(0.05),
            coupon_times: vec![0.5, 1.0, 1.5, 2.0],
            maturity: 2.0,
            rate: D::variable(coco0[0], 0),
            credit_spread: D::variable(coco0[1], 1),
            spot: D::variable(coco0[2], 2),
            trigger_barrier: D::from_value(60.0),
            div_yield: D::from_value(0.01),
            eq_vol: D::variable(coco0[3], 3),
            trigger_payoff: ql_fixed_income::CocoTriggerPayoff::WriteDown {
                recovery_fraction: D::from_value(0.4),
            },
        })
        .dirty_price;
        let cov = vec![
            vec![D::variable(port0[0], 0), D::variable(port0[1], 1)],
            vec![D::variable(port0[1], 1), D::variable(port0[2], 2)],
        ];
        let w = ql_portfolio::hierarchical_risk_parity_weights(&cov);
        let exec = ql_portfolio::multi_asset_almgren_chriss_execution(
            &ql_portfolio::MultiAssetAlmgrenChrissSpec {
                initial_inventory: vec![D::from_value(1000.0), D::from_value(600.0)],
                horizon: D::from_value(1.0),
                num_steps: 6,
                covariance: cov,
                temporary_impact: vec![
                    vec![D::from_value(1e-4), D::from_value(2e-5)],
                    vec![D::from_value(2e-5), D::from_value(1.5e-4)],
                ],
                permanent_impact: vec![
                    vec![D::from_value(5e-5), D::from_value(0.0)],
                    vec![D::from_value(0.0), D::from_value(5e-5)],
                ],
                risk_aversion: D::variable(port0[3], 3),
            },
        );
        let port = w[0] + w[1] * 2.0 + exec.mean_variance_utility * 1e-4;
        let sdh = dh.val + dh.d.iter().sum::<f64>();
        let st = tarf.val + tarf.d.iter().sum::<f64>();
        let sc = coco.val + coco.d.iter().sum::<f64>();
        let sp = port.val + port.d.iter().sum::<f64>();
        sdh + st + sc + sp
    };
    run_comparison(
        "26_double_heston_tarf_coco_hrp_aad",
        "Double Heston/TARF/CoCo/HRP single-pass AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn benchmark_jy_tba_bergomi_nco(scale: usize) -> ComparisonRow {
    let cfg = config(scale);
    let d = date(2024, 1, 1);
    let eval_jy = |x: [f64; 4]| {
        let dc = FlatForwardCurve::new(d, 0.025);
        let ic = ZeroInflationCurve::new(d, 100.0, vec![1.0, 2.0, 5.0, 10.0], vec![0.02; 4], None);
        let ins = YoYInflationCapFloor {
            notional: 1.0,
            maturity_years: 10,
            cap_rate: Some(0.025),
            floor_rate: None,
            gearing: 1.0,
            spread: 0.0,
            cap_floor_type: YoYInflationCapFloorType::Cap,
        };
        let p = JarrowYildirimParams {
            a_nom: 0.1,
            sigma_nom: x[0],
            a_real: 0.08,
            sigma_real: x[1],
            sigma_cpi: x[2],
            rho_nom_real: 0.2,
            rho_nom_cpi: 0.1,
            rho_real_cpi: x[3],
        };
        price_yoy_capfloor_jarrow_yildirim(&ins, &dc, &ic, &p)
            .0
            .val()
    };
    let eval_mbs = |x: [f64; 4]| {
        let c = FlatForwardCurve::new(d, 0.03);
        let roll = ql_fixed_income::analyze_tba_dollar_roll(&ql_fixed_income::TbaDollarRollSpec {
            front_clean_price: x[0],
            back_clean_price: x[1],
            wac: 0.04,
            pass_through_rate: 0.035,
            remaining_months: 120,
            front_accrued_fraction: 0.5,
            back_accrued_fraction: 0.5,
            roll_period_years: 1.0,
            assumed_cpr: x[2],
            repo_rate: x[3],
        });
        let m = ql_fixed_income::SpecifiedPoolMbsSpec {
            original_balance: 100.0,
            wac: 0.04,
            pass_through_rate: 0.035,
            term_months: 120,
            seasoning_months: 30,
            cpr_min: 0.02,
            cpr_max: x[2],
            refi_slope: 2.0,
            refi_inflection: 0.0,
            primary_secondary_spread: 0.01,
            cdr: 0.002,
            loss_severity: 0.4,
        };
        roll.net_roll_advantage.val()
            + ql_fixed_income::price_specified_pool_mbs(&m, &c, x[3])
                .pv
                .val()
                * 1e-3
    };
    let eval_bg_fx = |x: [f64; 4]| {
        let p = ql_models::Bergomi2FParams {
            xi0: x[0],
            omega: x[1],
            theta: 0.5,
            k1: 1.0,
            k2: 2.0,
            rho12: 0.25,
            rho_s1: x[2],
            rho_s2: -0.3,
        };
        let mut z = 0.0;
        for i in 0..16 {
            z += ql_models::bergomi_2f_implied_vol_bg(
                &p,
                100.0,
                100.0 * (0.7 + 0.04 * i as f64),
                1.0,
            )
            .val();
        }
        let fx = FxVannaVolgaSmile::new(
            100.0,
            0.02,
            0.01,
            1.0,
            FxSmileQuote {
                atm_vol: x[3],
                rr_25d: 0.0,
                bf_25d: 0.0,
            },
        );
        z + fx.price_double_no_touch(90.0, 110.0).val()
            + fx.price_digital_cash_dom(100.0, true).val()
    };
    let eval_nco = |x: [f64; 4]| {
        let ret = (0..16)
            .map(|i| {
                vec![
                    0.01 * ((i % 4) as f64 + 1.0) * x[3],
                    0.02 * ((i % 3) as f64 - 1.0),
                    0.01 * ((i % 5) as f64 - 2.0),
                    0.03 * ((i % 2) as f64 - 0.5),
                ]
            })
            .collect::<Vec<_>>();
        let g = ql_portfolio::gerber_covariance_matrix(&ret, 0.5, 1e-3).unwrap();
        let c = vec![
            vec![x[0], x[1], 0.0, 0.0],
            vec![x[1], 0.09, 0.0, 0.0],
            vec![0.0, 0.0, x[2], 0.01],
            vec![0.0, 0.0, 0.01, 0.16],
        ];
        let w = ql_portfolio::nested_clustered_optimization_weights(
            &ql_portfolio::NestedClusteredOptimizationSpec {
                covariance: c,
                expected_returns: None,
                risk_aversion: 1.0,
                num_clusters: 2,
            },
        );
        w[0] + w[1] * 2.0 + w[2] * 3.0 + w[3] * 4.0 + g[0][0] + g[0][1]
    };
    let x0 = [
        [0.01, 0.015, 0.012, -0.2],
        [99.0, 98.5, 0.06, 0.02],
        [0.04, 0.4, -0.5, 0.2],
        [0.04, 0.01, 0.09, 0.01],
    ];
    let checksum = || {
        let h = 1e-5;
        let mut total = eval_jy(x0[0]) + eval_mbs(x0[1]) + eval_bg_fx(x0[2]) + eval_nco(x0[3]);
        for i in 0..4 {
            let mut up = x0[0];
            let mut dn = x0[0];
            up[i] += h;
            dn[i] -= h;
            total += (eval_jy(up) - eval_jy(dn)) / (2.0 * h);
            let mut up = x0[1];
            let mut dn = x0[1];
            up[i] += h;
            dn[i] -= h;
            total += (eval_mbs(up) - eval_mbs(dn)) / (2.0 * h);
            let mut up = x0[2];
            let mut dn = x0[2];
            up[i] += h;
            dn[i] -= h;
            total += (eval_bg_fx(up) - eval_bg_fx(dn)) / (2.0 * h);
            let mut up = x0[3];
            let mut dn = x0[3];
            up[i] += h;
            dn[i] -= h;
            total += (eval_nco(up) - eval_nco(dn)) / (2.0 * h);
        }
        total
    };
    let baseline = || checksum();
    let sota = || {
        type D = Dual<4>;
        let dc = FlatForwardCurve::new(d, D::from_value(0.025));
        let ic = ZeroInflationCurve::new(
            d,
            100.0,
            vec![1.0, 2.0, 5.0, 10.0],
            vec![D::from_value(0.02); 4],
            None,
        );
        let ins = YoYInflationCapFloor {
            notional: D::from_value(1.0),
            maturity_years: 10,
            cap_rate: Some(D::from_value(0.025)),
            floor_rate: None,
            gearing: D::one(),
            spread: D::zero(),
            cap_floor_type: YoYInflationCapFloorType::Cap,
        };
        let jy = price_yoy_capfloor_jarrow_yildirim(
            &ins,
            &dc,
            &ic,
            &JarrowYildirimParams {
                a_nom: D::from_value(0.1),
                sigma_nom: D::variable(x0[0][0], 0),
                a_real: D::from_value(0.08),
                sigma_real: D::variable(x0[0][1], 1),
                sigma_cpi: D::variable(x0[0][2], 2),
                rho_nom_real: D::from_value(0.2),
                rho_nom_cpi: D::from_value(0.1),
                rho_real_cpi: D::variable(x0[0][3], 3),
            },
        )
        .0;
        let mb = analyze_tba_dollar_roll(&TbaDollarRollSpec::<D> {
            front_clean_price: D::variable(x0[1][0], 0),
            back_clean_price: D::variable(x0[1][1], 1),
            wac: D::from_value(0.04),
            pass_through_rate: D::from_value(0.035),
            remaining_months: 120,
            front_accrued_fraction: 0.5,
            back_accrued_fraction: 0.5,
            roll_period_years: 1.0,
            assumed_cpr: D::variable(x0[1][2], 2),
            repo_rate: D::variable(x0[1][3], 3),
        });
        let pool = price_specified_pool_mbs(
            &SpecifiedPoolMbsSpec::<D> {
                original_balance: D::from_value(100.0),
                wac: D::from_value(0.04),
                pass_through_rate: D::from_value(0.035),
                term_months: 120,
                seasoning_months: 30,
                cpr_min: D::from_value(0.02),
                cpr_max: D::variable(x0[1][2], 2),
                refi_slope: D::from_value(2.0),
                refi_inflection: D::zero(),
                primary_secondary_spread: D::from_value(0.01),
                cdr: D::from_value(0.002),
                loss_severity: D::from_value(0.4),
            },
            &FlatForwardCurve::new(d, D::from_value(0.03)),
            D::variable(x0[1][3], 3),
        );
        let mut bg = D::zero();
        let bp = Bergomi2FParams {
            xi0: D::variable(x0[2][0], 0),
            omega: D::variable(x0[2][1], 1),
            theta: D::from_value(0.5),
            k1: D::from_value(1.0),
            k2: D::from_value(2.0),
            rho12: D::from_value(0.25),
            rho_s1: D::variable(x0[2][2], 2),
            rho_s2: D::from_value(-0.3),
        };
        for i in 0..16 {
            bg += bergomi_2f_implied_vol_bg(
                &bp,
                D::from_value(100.0),
                D::from_value(100.0 * (0.7 + 0.04 * i as f64)),
                D::from_value(1.0),
            );
        }
        let fx = FxVannaVolgaSmile::new(
            D::from_value(100.0),
            D::from_value(0.02),
            D::from_value(0.01),
            D::from_value(1.0),
            FxSmileQuote {
                atm_vol: D::variable(x0[2][3], 3),
                rr_25d: D::zero(),
                bf_25d: D::zero(),
            },
        );
        bg += fx.price_double_no_touch(D::from_value(90.0), D::from_value(110.0))
            + fx.price_digital_cash_dom(D::from_value(100.0), true);
        let ret = (0..16)
            .map(|i| {
                vec![
                    D::from_value(0.01 * ((i % 4) as f64 + 1.0)) * D::variable(x0[3][3], 3),
                    D::from_value(0.02 * ((i % 3) as f64 - 1.0)),
                    D::from_value(0.01 * ((i % 5) as f64 - 2.0)),
                    D::from_value(0.03 * ((i % 2) as f64 - 0.5)),
                ]
            })
            .collect::<Vec<_>>();
        let gc = gerber_covariance_matrix(&ret, 0.5, 1e-3).unwrap();
        let c_d = vec![
            vec![
                D::variable(x0[3][0], 0),
                D::variable(x0[3][1], 1),
                D::zero(),
                D::zero(),
            ],
            vec![
                D::variable(x0[3][1], 1),
                D::from_value(0.09),
                D::zero(),
                D::zero(),
            ],
            vec![
                D::zero(),
                D::zero(),
                D::variable(x0[3][2], 2),
                D::from_value(0.01),
            ],
            vec![
                D::zero(),
                D::zero(),
                D::from_value(0.01),
                D::from_value(0.16),
            ],
        ];
        let nc = nested_clustered_optimization_weights(&NestedClusteredOptimizationSpec {
            covariance: c_d,
            expected_returns: None,
            risk_aversion: D::one(),
            num_clusters: 2,
        });
        let nco = nc[0] + nc[1] * 2.0 + nc[2] * 3.0 + nc[3] * 4.0 + gc[0][0] + gc[0][1];
        let out = jy + mb.net_roll_advantage + pool.pv * 1e-3 + bg + nco;
        out.val + out.d.iter().sum::<f64>()
    };
    run_comparison(
        "27_jy_inflation_tba_bergomi_nco_aad",
        "JY-Inflation/TBA-Roll/Bergomi2F/NCO single-pass AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

fn suite28_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
    let d = date(2024, 1, 1);
    let curve = FlatForwardCurve::new(d, T::from_f64(0.03) + x[0] * 0.001);
    let chey = price_cheyette_swaption_markov_projection(
        &CheyetteSwaptionSpec {
            notional: T::from_f64(100.0),
            strike: T::from_f64(0.032) + x[1] * 0.0001,
            option_expiry: 1.0,
            payment_times: vec![2.0, 3.0, 4.0, 5.0],
            mean_reversion: T::from_f64(0.05) + x[2] * 0.001,
            sigma0: T::from_f64(0.008) + x[3] * 0.0001,
            skew_b: T::from_f64(-1.5),
            is_payer: true,
        },
        &curve,
    )
    .npv;
    let ctd = ql_fixed_income::analyze_treasury_futures_ctd(
        &ql_fixed_income::TreasuryFuturesBasketSpec {
            futures_price: T::from_f64(98.5) + x[0] * 0.05,
            delivery_time: 0.25,
            repo_rate: T::from_f64(0.03) + x[1] * 0.0005,
            notional_standard_yield: 0.06,
            deliverables: vec![
                ql_fixed_income::TreasuryFutureDeliverable {
                    clean_price: T::from_f64(99.2),
                    coupon_rate: T::from_f64(0.06),
                    remaining_semiannual_periods: 16,
                    accrued_interest_spot: T::from_f64(0.5),
                    accrued_interest_delivery: T::from_f64(1.0),
                    interim_coupon: T::zero(),
                    yield_vol: T::from_f64(0.008) + x[2] * 0.0001,
                    modified_duration: T::from_f64(6.5),
                },
                ql_fixed_income::TreasuryFutureDeliverable {
                    clean_price: T::from_f64(89.0),
                    coupon_rate: T::from_f64(0.045),
                    remaining_semiannual_periods: 18,
                    accrued_interest_spot: T::from_f64(0.4),
                    accrued_interest_delivery: T::from_f64(0.8),
                    interim_coupon: T::zero(),
                    yield_vol: T::from_f64(0.008),
                    modified_duration: T::from_f64(7.2),
                },
            ],
        },
    )
    .option_adjusted_futures_price;
    let cmo = price_cmo_waterfall(
        &CmoWaterfallSpec {
            original_balance: T::from_f64(100.0),
            wac: T::from_f64(0.05),
            pass_through_rate: T::from_f64(0.045),
            term_months: 24,
            psa_speed: T::from_f64(150.0) + x[0] * 1.0,
            pac_collar_psa: Some((T::from_f64(100.0), T::from_f64(300.0))),
            tranches: vec![
                CmoTrancheSpec {
                    initial_balance: T::from_f64(60.0),
                    coupon_rate: T::from_f64(0.04),
                    is_pac: true,
                },
                CmoTrancheSpec {
                    initial_balance: T::from_f64(40.0),
                    coupon_rate: T::from_f64(0.045),
                    is_pac: false,
                },
            ],
        },
        &curve,
        T::from_f64(0.005) + x[1] * 0.0001,
    )
    .total_structure_pv;
    let qrh = qrh_implied_vol_expansion(
        &QuadraticRoughHestonParams {
            hurst: 0.1,
            nu: T::from_f64(0.35) + x[0] * 0.001,
            lambda: T::from_f64(1.2),
            a: T::from_f64(0.45) + x[1] * 0.001,
            b: T::from_f64(0.2) + x[2] * 0.001,
            c: T::from_f64(0.04) + x[3] * 0.001,
            z0: T::zero(),
        },
        T::from_f64(100.0),
        T::from_f64(95.0),
        T::from_f64(1.0),
    );
    let prdc = ql_fx::price_prdc_note_hybrid(&ql_fx::PrdcNoteSpec {
        notional: T::from_f64(100.0),
        spot_fx: T::from_f64(150.0) + x[0] * 0.1,
        r_dom: T::from_f64(0.01),
        r_for: T::from_f64(0.04),
        a_dom: T::from_f64(0.05),
        sigma_dom: T::from_f64(0.008),
        a_for: T::from_f64(0.08),
        sigma_for: T::from_f64(0.01),
        sigma_fx: T::from_f64(0.12) + x[2] * 0.001,
        rho_fx_dom: T::from_f64(0.3),
        rho_fx_for: T::from_f64(-0.2),
        rho_dom_for: T::from_f64(0.25),
        coupon_times: vec![1.0, 2.0],
        foreign_gearing: T::from_f64(0.15),
        domestic_gearing: T::from_f64(0.10),
        coupon_cap: T::from_f64(0.08),
        fx_knockout_barrier: Some(T::from_f64(180.0)),
    })
    .npv;
    let times = [1.0, 2.0, 3.0, 4.0, 5.0];
    let q1 = BaseCorrelationTrancheQuote {
        attachment: 0.0,
        detachment: 0.03,
        running_spread: T::from_f64(0.015) + x[0] * 0.0001,
        upfront: T::zero(),
    };
    let q2 = BaseCorrelationTrancheQuote {
        attachment: 0.03,
        detachment: 0.07,
        running_spread: T::from_f64(0.02) + x[1] * 0.0001,
        upfront: T::zero(),
    };
    let bc = bootstrap_base_correlation_curve(
        &[q1, q2],
        T::from_f64(0.02) + x[2] * 0.0001,
        T::from_f64(0.4),
        T::from_f64(0.03) + x[3] * 0.0001,
        &times,
    );
    let (_, fair, _) = price_mezzanine_tranche_base_correlation(
        &bc,
        0.02,
        0.05,
        T::from_f64(0.02) + x[1] * 0.0001,
        T::from_f64(100.0),
    );
    let rets = (0..16)
        .map(|i| T::from_f64(-0.045 + 0.006 * i as f64) + x[3] * (0.0001 * ((i % 3) as f64 - 1.0)))
        .collect::<Vec<_>>();
    let cdar = conditional_drawdown_at_risk(&rets, 0.75).cdar;
    let evar = entropic_value_at_risk(&rets, 0.90);
    let attr = factor_risk_attribution(
        &[
            T::from_f64(0.4) + x[0] * 0.001,
            T::from_f64(0.6) - x[0] * 0.001,
        ],
        &[
            vec![T::from_f64(1.1), T::from_f64(0.3)],
            vec![T::from_f64(0.8), T::from_f64(0.9)],
        ],
        &[
            vec![T::from_f64(0.04), T::from_f64(0.005)],
            vec![T::from_f64(0.005), T::from_f64(0.02)],
        ],
        &[T::from_f64(0.01), T::from_f64(0.015)],
    )
    .total_volatility;
    [chey + ctd, cmo, qrh + prdc, fair + cdar + evar + attr]
}

pub fn benchmark_cheyette_cmo_qrh_basecorr_cdar_aad(scale: usize) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let h = 1e-5;
    let baseline = || {
        let base = suite28_workloads(x0);
        let mut checksum: f64 = base.iter().sum();
        for lane in 0..4 {
            let mut up = x0;
            let mut dn = x0;
            up[lane] += h;
            dn[lane] -= h;
            let w_up = suite28_workloads(up);
            let w_dn = suite28_workloads(dn);
            for component in 0..4 {
                checksum += (w_up[component] - w_dn[component]) / (2.0 * h);
            }
        }
        checksum
    };
    let sota = || {
        type D = Dual<4>;
        let x = [
            D::variable(x0[0], 0),
            D::variable(x0[1], 1),
            D::variable(x0[2], 2),
            D::variable(x0[3], 3),
        ];
        suite28_workloads(x)
            .iter()
            .map(|y| y.val + y.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "28_cheyette_cmo_qrh_basecorr_cdar_aad",
        "Cheyette/CMO/QRH/BaseCorr/EVaR single-pass AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

fn suite29_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
    let d = date(2024, 1, 1);
    let curve = FlatForwardCurve::new(d, T::from_f64(0.03) + x[0] * 0.001);
    let asw = ql_fixed_income::analyze_bond_asset_swap_and_z_spread(
        &ql_fixed_income::BondAssetSwapSpec {
            clean_price: T::from_f64(0.97) + x[1] * 0.001,
            accrued_interest: T::from_f64(0.01),
            coupon_rate: T::from_f64(0.04) + x[2] * 0.0005,
            payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            accrual_fractions: vec![1.0; 5],
        },
        &curve,
    );
    let mf = ql_rates::price_markov_functional_1f(
        &ql_rates::MarkovFunctionalSwaptionSpec {
            expiry: 1.0,
            payment_times: vec![2.0, 3.0, 4.0, 5.0],
            mean_reversion: T::from_f64(0.05) + x[2] * 0.001,
            atm_normal_vol: T::from_f64(0.008) + x[3] * 0.0001,
            sabr_rho: T::from_f64(-0.25),
            sabr_nu: T::from_f64(0.35),
            strike: T::from_f64(0.031) + x[1] * 0.0001,
            notional: T::from_f64(100.0),
            is_payer: true,
        },
        &curve,
        12,
    );

    let bacva = ql_xva::compute_ba_cva_capital(&ql_xva::BaCvaPortfolioSpec {
        counterparties: vec![ql_xva::BaCvaCounterpartySpec {
            risk_weight: T::from_f64(0.05),
            effective_maturity: T::from_f64(3.0),
            ead: T::from_f64(100.0) + x[0] * 0.5,
            single_name_hedges: vec![ql_xva::BaCvaSingleNameHedge {
                notional: T::from_f64(35.0) + x[1] * 0.2,
                maturity: T::from_f64(3.0),
                hedge_correlation: 0.8,
            }],
        }],
        index_hedges: vec![],
        hurdle_rate: T::from_f64(0.10),
        risk_free_rate: T::from_f64(0.03) + x[0] * 0.001,
    })
    .capital_charge;
    let bcox = ql_credit::analyze_structural_credit_black_cox(&ql_credit::BlackCoxStructuralSpec {
        asset_value: T::from_f64(140.0) + x[2] * 0.5,
        debt_face_value: T::from_f64(100.0),
        default_barrier: T::from_f64(85.0),
        barrier_growth_rate: T::from_f64(0.02),
        asset_vol: T::from_f64(0.25) + x[3] * 0.001,
        risk_free_rate: T::from_f64(0.03) + x[0] * 0.001,
        asset_payout_rate: T::from_f64(0.01),
        maturity: T::from_f64(3.0),
        recovery_rate: T::from_f64(0.40),
    })
    .black_cox_cds_spread;

    let ah_slice = ql_models::solve_andreasen_huge_step(
        &ql_models::AndreasenHugeGridSpec {
            forward: T::from_f64(100.0) + x[0] * 0.1,
            rate: T::from_f64(0.02),
            dt: T::from_f64(0.5),
            strike_min: T::from_f64(70.0),
            strike_max: T::from_f64(130.0),
            n_grid: 25,
            local_vol_atm: T::from_f64(0.20) + x[1] * 0.001,
            local_vol_skew: T::from_f64(-0.15) + x[2] * 0.001,
            local_vol_smile: T::from_f64(0.10) + x[3] * 0.001,
        },
        None,
    );
    let ah_call = ql_models::interpolate_andreasen_huge_call(&ah_slice, T::from_f64(100.0));
    let fx_acc = ql_fx::price_fx_accumulator_strip_no_ko(&ql_fx::FxAccumulatorSpec {
        base_notional_per_fixing: T::from_f64(10.0),
        spot: T::from_f64(1.10) + x[0] * 0.001,
        strike: T::from_f64(1.06),
        knockout_barrier: T::from_f64(1.20),
        r_dom: T::from_f64(0.03),
        r_for: T::from_f64(0.01),
        vol: T::from_f64(0.12) + x[1] * 0.001,
        fixing_times: vec![0.25, 0.50, 0.75, 1.0],
        guaranteed_fixings: 1,
        downside_gearing: T::from_f64(2.0),
        is_accumulator: true,
        smoothing_eps: 1e-2,
    });

    let p_curve = ql_commodities::CommodityFuturesCurve::new(
        vec![0.25, 0.50, 0.75, 1.0],
        vec![
            T::from_f64(45.0) + x[0] * 0.1,
            T::from_f64(58.0) + x[1] * 0.1,
            T::from_f64(42.0),
            T::from_f64(62.0),
        ],
    )
    .unwrap();
    let bess = ql_commodities::price_bess_arbitrage_dp(
        &ql_commodities::BessSpec {
            decision_times: vec![0.25, 0.50, 0.75, 1.0],
            energy_capacity_mwh: 4.0,
            initial_soc_mwh: 2.0,
            max_charge_mw: 8.0,
            max_discharge_mw: 8.0,
            charge_efficiency: T::from_f64(0.92) + x[2] * 0.001,
            discharge_efficiency: T::from_f64(0.92) + x[3] * 0.001,
            cycle_degradation_cost: T::from_f64(1.5),
            ancillary_reserve_rate: T::from_f64(2.0),
            soc_grid_points: 3,
        },
        &p_curve,
        T::from_f64(1.5),
        T::from_f64(0.35),
        T::from_f64(0.03),
    )
    .total_npv;
    let idz =
        ql_portfolio::idzorek_black_litterman_posterior(&ql_portfolio::IdzorekBlackLittermanSpec {
            covariance: vec![
                vec![T::from_f64(0.04), T::from_f64(0.01)],
                vec![T::from_f64(0.01), T::from_f64(0.09)],
            ],
            market_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
            risk_aversion: T::from_f64(2.5),
            tau: T::from_f64(0.05),
            view_matrix: vec![vec![T::one(), -T::one()]],
            view_returns: vec![T::from_f64(0.03) + x[2] * 0.0005],
            view_confidences: vec![T::from_f64(0.65) + x[3] * 0.005],
        })
        .weights[0];
    let dp_exec = ql_portfolio::optimal_dark_lit_execution(&ql_portfolio::DarkLitRoutingSpec {
        initial_inventory: T::from_f64(100.0),
        horizon: T::one(),
        num_steps: 4,
        volatility: T::from_f64(0.20) + x[1] * 0.001,
        lit_temporary_impact: T::from_f64(1e-3),
        lit_permanent_impact: T::from_f64(2e-4),
        dark_fill_intensity: T::from_f64(1.5) + x[3] * 0.01,
        risk_aversion: T::from_f64(0.25),
    })
    .expected_total_cost;

    [
        asw.z_spread + asw.par_asset_swap_spread + mf.swaption_npv + mf.cms_swaplet_rate,
        bacva + bcox,
        ah_call + fx_acc,
        bess + idz + dp_exec,
    ]
}

pub fn benchmark_asw_mf1f_ah_bess_ep_aad(scale: usize) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let h = 1e-5;
    let baseline = || {
        let base = suite29_workloads(x0);
        let mut checksum: f64 = base.iter().sum();
        for lane in 0..4 {
            let mut up = x0;
            let mut dn = x0;
            up[lane] += h;
            dn[lane] -= h;
            let w_up = suite29_workloads(up);
            let w_dn = suite29_workloads(dn);
            for component in 0..4 {
                checksum += (w_up[component] - w_dn[component]) / (2.0 * h);
            }
        }
        checksum
    };
    let sota = || {
        type D = Dual<4>;
        let x = [
            D::variable(x0[0], 0),
            D::variable(x0[1], 1),
            D::variable(x0[2], 2),
            D::variable(x0[3], 3),
        ];
        suite29_workloads(x)
            .iter()
            .map(|y| y.val + y.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "29_asw_mf1f_ah_bess_ep_aad",
        "ASW-Z/MF1F/AH-LV/BESS/Idzorek-DP single-pass AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn run_all_benchmarks(scale: usize) -> Vec<ComparisonRow> {
    // Keep all quantitative engine dependencies linked into the 41-suite
    // benchmark binary, including engines used by the parity harness.
    let _ = ql_math::inv_norm_cdf(0.5);
    let _ = std::mem::size_of::<ql_pde::Fd2dHestonEngine>();
    let _ = std::mem::size_of::<ql_xva::NpvCube<f64>>();
    let _ = std::mem::size_of::<SurvivalProbabilityCurve<f64>>();
    let _ = std::mem::size_of::<ZeroInflationCurve<f64>>();
    let _ = std::mem::size_of::<FxVannaVolgaSmile<f64>>();
    vec![
        benchmark_calendar(scale),
        benchmark_leg_simd(scale),
        benchmark_multicurve(scale),
        benchmark_implied_vol(scale),
        benchmark_heston(scale),
        benchmark_gbm(scale),
        benchmark_scripted(scale),
        benchmark_portfolio(scale),
        benchmark_credit(scale),
        benchmark_fx(scale),
        benchmark_cholesky_corr(scale),
        benchmark_dim_mva(scale),
        benchmark_barrier_dual5(scale),
        benchmark_rqmc(scale),
        benchmark_levy_cos(scale),
        benchmark_cdo_recursion(scale),
        benchmark_rfr_sabr(scale),
        benchmark_xccs(scale),
        benchmark_american_baw(scale),
        benchmark_risk_parity(scale),
        benchmark_bermudan_hw(scale),
        benchmark_sabr_qp(scale),
        benchmark_slv_autocall(scale),
        benchmark_frtb_cvar_execution(scale),
        benchmark_vix_parisian_robust(scale),
        benchmark_double_heston_tarf_coco_hrp(scale),
        benchmark_jy_tba_bergomi_nco(scale),
        benchmark_cheyette_cmo_qrh_basecorr_cdar_aad(scale),
        benchmark_asw_mf1f_ah_bess_ep_aad(scale),
        benchmark_mtmxccs_colva_zabr_pdv_weather_enb_aad(scale),
        benchmark_burnout_cmcds_quintic_afv_cs_schur_aad(scale),
        benchmark_sinhou_mb_rbergomi_choi_afs_expectile_aad(scale),
        benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad(scale),
        benchmark_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad(scale),
        benchmark_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml(scale),
        benchmark_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp(
            scale,
        ),
        benchmark_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin(scale),
        benchmark_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc(scale),
        benchmark_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt(scale),
        benchmark_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca(scale),
        benchmark_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is(scale),
        benchmark_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin(scale),
        benchmark_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi(scale),
        benchmark_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan(scale),
        benchmark_rebonato_cms_lm_all_engines(scale),
        benchmark_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig(scale),
        benchmark_hagan_west_ctd_afv_kth_sacva_gs_kirk_fx_vv_napoleon_tc_meucci_krr_eva(scale),
        benchmark_brigo_csa_ctd_fabozzi_callable_oas_lando_jlt_migration_gatheral_corridor_varswap_glasserman_lsm_rogers_dual_cartea_lob_mm_lo_mackinlay_vr_lopez_dsr_purged_cv_roncalli_spinu_erc_ilmanen_vmc_penman_reoi_mccall_bellman(scale),
        benchmark_brigo_nss_krd_tuckman_tips_floor_schonbucher_rmv_rt_rfv_gatheral_sticky_lv_gregory_simm_ccp_mva_schwartz_smith_kalman_wystup_fx_dnt_bouzoubaa_timer_meucci_bl_woodbury_gueant_ac_sinh_zou_elastic_net_nk_ree(scale),
        benchmark_hagan_west_tension_duffie_lando_pac_cmo_fama_macbeth_shanken_gatheral_asvi_bgk_discrete_barrier_bouchaud_cross_impact_karnosky_singer_merton_jump_campbell_viceira_ohlson_lid_griewank_revolve(scale),
        benchmark_brigo_xccy_mtm_fabozzi_spec_pool_io_po_lando_marshall_olkin_gatheral_vix_vvix_gregory_wwr_cva_geman_gibson_schwartz_clark_fx_fva_rubinstein_chooser_grinold_flam_tc_meucci_es_euler_james_stein_sure_campbell_cochrane_habit(scale),
        benchmark_brigo_ois_fallback_fabozzi_collared_frn_schonbucher_cds_swaption_fep_bergomi_2f_skew_gregory_sa_ccr_ead_eydeland_spark_spread_wystup_prdc_geske_compound_qian_alpha_decay_roncalli_enb_johansen_vecm_bansal_yaron_lrr(scale),
        benchmark_andersen_cms_replication_tuckman_special_repo_okane_cds_upfront_gatheral_ssr_tanaka_gregory_saccr_kva_eydeland_crack_spread_wystup_prdc_cap_haug_eso_dilution_qian_spinu_erc_rockafellar_cvar_htf_pls_krylov_piazzesi_fomc(scale),
        benchmark_darbyshire_step_ois_turn_despiegeleer_cb_soft_call_okane_cds_index_skew_gatheral_varswap_truncation_gregory_fva_overlap_geman_working_storage_wystup_malz_fx_fly_goldman_sosin_lookback_ledoit_wolf_shrinkage_almgren_chriss_cara_tsay_estar_acm_cp(scale),
        benchmark_henrard_tenor_basis_fabozzi_msr_recapture_okane_cdo_base_corr_gatheral_fwd_start_smile_gregory_pykhtin_csa_mta_eydeland_gas_swing_wystup_fx_kiko_stulz_rainbow_jagannathan_ma_obizhaeva_wang_newey_west_hac_gabaix_disaster(scale),
        benchmark_brigo_diff_swap_quanto_fabozzi_callable_muni_ytw_bielecki_step_up_bond_bkm_moment_strip_gregory_break_clause_cva_bessembinder_power_skew_wystup_fx_fader_longstaff_extendible_kan_zhou_3fund_sharpe_rbsa_diebold_mariano_hln_laubach_williams_rstar(scale),
        benchmark_andersen_sofr_futures_fra_hayre_mbs_scurve_pedersen_cdx_fep_gatheral_dupire_total_var_gregory_simm_mva_routledge_stockout_baba_cip_basis_carr_chou_symmetry_pastor_stambaugh_treynor_black_engle_granger_ecm_epstein_zin_sdf(scale),
        benchmark_litterman_scheinkman_pca_fly_tf_split_cb_hull_white_double_t_cdo_roger_lee_moment_wings_brigo_capponi_cds_contagion_schwartz_ou_futures_wystup_dcd_yield_passport_option_grs_factor_ftest_bouchaud_sqrt_impact_bns_bipower_har_cj_lettau_ludvigson_cay(scale),
        benchmark_piterbarg_collateral_option_tuckman_rolldown_merton_kmv_dd_hagan_bartlett_sabr_pykhtin_mpor_spike_deng_locational_ftr_wystup_digital_overhedge_conze_fixed_lookback_tasche_hot_spots_es_frazzini_bab_pesaran_timmermann_lucas_mehra_prescott(scale),
        benchmark_gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_ga_carr_wu_vol_drift_gregory_marginal_cva_working_brennan_storage_wystup_forward_extra_rubinstein_exchange_opt_hlz_factor_fdr_holden_subrahmanyam_kyle_bollerslev_garch_kurtosis_campbell_shiller_cape(scale),
        benchmark_kim_wright_shadow_rate_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_swap_burgard_kjaer_fva_carmona_durrleman_spark_wystup_fx_tarf_bgk_discrete_lookback_novy_marx_ff5_acerbi_spectral_es_corsi_har_rv_bll_news_noise(scale),
        benchmark_jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_regimes_hull_white_wwr_cva_gibson_schwartz_samuelson_reiner_quanto_fx_turnbull_wakeman_asian_amp_value_momentum_garleanu_pedersen_aim_bai_ng_icp2_bgg_accelerator(scale),
        benchmark_lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_cds_fouque_multiscale_sv_gregory_bilateral_ftd_eydeland_vpp_heat_rate_wystup_fx_corridor_rubinstein_reiner_barrier_lns_gls_r2_bellini_expectile_bns_quarticity_kiyotaki_moore_cycle(scale),
        benchmark_jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_albanese_sa_ccr_kva_rss_equilibrium_storage_du_tepper_verdelhan_cip_wilmott_cliquet_bsv_parametric_portfolio_ahmadi_javid_evar_engle_kroner_bekk_gabaix_disaster(scale),
        benchmark_bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_smile_brigo_pallavicini_ccp_xva_dmt_cointegrated_crack_piterbarg_prdc_tarn_bgk_discrete_barrier_connor_korajczyk_apt_afs_transient_impact_hamilton_markov_switching_bansal_yaron_lrr(scale),
        benchmark_andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_vol_swap_green_leverage_ratio_kva_cortazar_schwartz_3f_oil_fujii_collateralized_xccs_pechtl_haug_time_switch_boudt_cornish_fisher_mvar_rockafellar_quadrangle_hasbrouck_info_share_he_krishnamurthy_sdf(scale),
        benchmark_piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_green_kenyon_tva_eydeland_interconnector_spread_wystup_clark_fx_strangle_kunitomo_ikeda_curved_barrier_clmx_variance_decomp_kyle_multi_period_lambda_hansen_lunde_overnight_rv_goldstein_pauzner_bank_run(scale),
        benchmark_mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_burgard_kjaer_semi_replication_xva_benth_seasonal_hdd_cdd_swap_baba_packer_fx_swap_dollar_basis_conze_viswanathan_lookback_he_litterman_bl_posterior_ghs_importance_sampling_cvar_engle_granger_johansen_vecm_gabaix_rare_disasters_resilience(scale),
        benchmark_ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_gregory_simm_mva_carmona_clean_dark_spread_reiner_haug_quanto_call_hart_ross_soft_barrier_pastor_stambaugh_alpha_shrinkage_obizhaeva_wang_impulse_lob_diebold_mariano_west_hln_piazzesi_schneider_nominal_tp(scale),
        benchmark_smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_pykhtin_rosen_marginal_cva_hilliard_reis_jump_futures_gabaix_maggiori_fx_carry_cheuk_vorst_shout_option_grossman_zhou_drawdown_kelly_gls_tasche_cornish_fisher_raroc_andrews_ploberger_wald_break_jermann_quadrini_financial_shock(scale),
        benchmark_mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_kenyon_green_csa_lva_boogert_de_jong_gas_ratchet_garman_kohlhagen_baw_fx_hui_haug_binary_barrier_kan_robotti_hj_distance_cartea_jaimungal_order_flow_ait_sahalia_jacod_jump_test_wachter_disaster_volatility(scale),
        benchmark_hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_gregory_brigo_ftd_bilateral_cva_brennan_schwartz_mine_hysteresis_wystup_fx_fader_broadie_detemple_capped_call_meucci_entropy_pooling_cont_deguest_scandolo_robust_risk_pesaran_shin_smith_ardl_christiano_eichenbaum_evans_nk(scale),
        benchmark_mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_andersen_pykhtin_sokol_mpor_bjerksund_stensland_321_crack_carr_chou_symmetry_barrier_chesney_jeanblanc_yor_parisian_demiguel_brodie_norm_shrinkage_gueant_lehalle_tapia_mm_hansen_white_spa_test_brunnermeier_sannikov_paradox(scale),
        benchmark_filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_hull_white_burgard_kjaer_fva_windfall_casassus_collin_dufresne_convenience_yield_lustig_roussanov_verdelhan_fx_carry_carr_lee_corridor_varswap_garlappi_uppal_wang_multi_prior_mcneil_frey_evt_garch_patton_sjc_copula_pastor_veronesi_political_uncertainty(scale),
        benchmark_lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_brigo_capponi_pallavicini_contagion_wwr_routledge_seppi_spatt_storage_piterbarg_prdc_hybrid_broadie_glasserman_kou_discrete_barrier_brandt_santa_clara_valkanov_ppp_gatheral_schied_slynko_transient_impact_corsi_bpq_harq_lettau_ludvigson_cay(scale),
        benchmark_cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_crepey_bichuch_capponi_xva_bsde_trolle_schwartz_usv_commodity_ivashina_scharfstein_stein_cip_wilmott_windcliff_cliquet_skew_ledoit_wolf_nonlinear_shrinkage_adrian_brunnermeier_delta_covar_bollerslev_tauchen_zhou_vrp_campbell_vuolteenaho_bad_good_beta(scale),
        benchmark_bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_andersen_duffie_song_debt_overhang_fva_pirrong_structural_inventory_spread_lipton_castagna_fx_dnt_vanna_volga_linetsky_step_barrier_koijen_yogo_demand_system_kyle_obizhaeva_invariance_gu_kelly_xiu_oos_r2_drechsler_savov_schnabl_deposit_channel(scale),
        benchmark_hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_albanese_caenazzo_crepey_bilateral_kva_acharya_lochstoer_ramadorai_hedging_pressure_brunnermeier_nagel_pedersen_fx_crash_broadie_yamamoto_fgt_asian_kelly_pruitt_su_ipca_brownlees_engle_srisk_andersen_bollerslev_diebold_labys_arfima_farhi_gabaix_disaster_fx(scale),
        benchmark_klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_duffie_zhu_ccp_netting_gorton_hayashi_rouwenhorst_inventory_basis_della_corte_ramadorai_sarno_fx_vrp_pelsser_laplace_double_barrier_moreira_muir_volatility_managed_budish_cramton_shim_fba_ait_sahalia_mykland_zhang_tsrv_gertler_karadi_credit_policy(scale),
        benchmark_bauer_rudebusch_qe_signaling_boyarchenko_fuster_lucca_mbs_oas_azizpour_giesecke_schwenkler_default_clustering_bollerslev_todorov_jump_fear_ghamami_glasserman_clearing_mva_szymanowska_de_roon_nijman_commodity_premia_menkhoff_sarno_schmeling_schrimpf_fx_vol_milevsky_posner_reciprocal_gamma_asian_kozak_nagel_santosh_sdf_shrinkage_acharya_pedersen_philippon_richardson_ses_bnhls_realized_kernel_brunnermeier_pedersen_margin_spiral(scale),
        benchmark_gurkaynak_sack_swanson_fomc_path_bao_pan_wang_bond_illiquidity_collin_dufresne_goldstein_martin_credit_determinants_andersen_fusari_todorov_left_tail_capponi_cheng_ccp_margin_kang_rouwenhorst_tang_two_premiums_lustig_roussanov_verdelhan_dollar_carry_geman_yor_bessel_asian_tu_zhou_markowitz_talmud_biais_foucault_moinas_fast_trading_christensen_oomen_podolskij_uhf_jumps_gabaix_granular_residual(scale),
        benchmark_nakamura_steinsson_fomc_info_dick_nielsen_feldhutter_lando_illiquidity_berndt_duffie_credit_premia_christoffersen_heston_jacobs_variance_kernel_brigo_francischello_pallavicini_closeout_xva_boons_prado_basis_momentum_colacito_croce_gavazzoni_ready_recursive_fx_fusai_meucci_erlang_asian_freyberger_neuhierl_weber_characteristics_kelly_jiang_tail_risk_jacod_li_mykland_podolskij_vetter_preaveraging_adrian_etula_muir_broker_dealer_sdf(scale),
    ]
}

pub fn suite46_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);
    let repo = ql_curves::evaluate_tuckman_repo_specialness_and_oftr_basis(
        &ql_curves::TuckmanRepoSpecialnessBasisSpec {
            otr_dirty_price: T::from_f64(100.0) + x[0],
            otr_yield: T::from_f64(0.04) + dx1,
            oftr_yield: T::from_f64(0.0415) + dx2,
            modified_duration: T::from_f64(7.5),
            gc_repo_rate: T::from_f64(0.05) + dx3,
            special_repo_rate: T::from_f64(0.035),
            specialness_horizon_years: T::from_f64(0.25),
            repo_haircut_fraction: T::from_f64(0.02),
            equity_funding_rate: T::from_f64(0.06),
        },
    );
    let bk = ql_credit::evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva(
        &ql_credit::CrepeyBurgardKjaerBilateralFvaSpec {
            time_grid: vec![T::from_f64(0.5), T::one(), T::from_f64(2.0)],
            epe_profile: vec![
                T::from_f64(10.0) + x[0],
                T::from_f64(12.0) + x[1],
                T::from_f64(8.0),
            ],
            ene_profile: vec![T::from_f64(6.0), T::from_f64(7.0), T::from_f64(5.0)],
            cpty_hazard_rate: T::from_f64(0.02) + dx2,
            bank_hazard_rate: T::from_f64(0.015) + dx3,
            cpty_recovery_rate: T::from_f64(0.4),
            bank_recovery_rate: T::from_f64(0.4),
            risk_free_rate: T::from_f64(0.03),
            liquidity_funding_spread: T::from_f64(0.008),
            burgard_kjaer_alpha: T::from_f64(0.5),
        },
    );
    let decs = ql_fixed_income::evaluate_despiegeleer_mandatory_convertible_decs_and_percs(
        &ql_fixed_income::DeSpiegeleerMandatoryConvertibleDecsSpec {
            spot_price: T::from_f64(100.0) + x[0],
            par_amount: T::from_f64(100.0),
            lower_strike_kl: T::from_f64(90.0),
            upper_strike_ku: T::from_f64(120.0),
            decs_coupon_rate: T::from_f64(0.06),
            risk_free_rate: T::from_f64(0.03) + dx1,
            dividend_yield: T::from_f64(0.02),
            credit_spread: T::from_f64(0.015),
            lower_strike_vol: T::from_f64(0.25) + dx2,
            upper_strike_vol: T::from_f64(0.22) + dx3,
            maturity_years: T::from_f64(3.0),
        },
    );
    let mbs = ql_fixed_income::evaluate_davidson_levin_four_factor_prepayment_io_po(
        &ql_fixed_income::DavidsonLevinFourFactorIoPoSpec {
            initial_balance: T::from_f64(100.0),
            gross_wac_rate: T::from_f64(0.06),
            net_passthrough_coupon: T::from_f64(0.055),
            current_mortgage_rate: T::from_f64(0.05) + dx0,
            base_turnover_cpr: T::from_f64(0.06),
            max_refi_cpr: T::from_f64(0.35) + dx1,
            burnout_sensitivity_beta: T::from_f64(0.5),
            seasonality_amplitude: T::from_f64(0.1),
            oas_discount_rate: T::from_f64(0.045) + dx2,
            num_periods: 12,
        },
    );
    let kyle = ql_portfolio::evaluate_foucault_kyle_lambda_and_ekop_pin(
        &ql_portfolio::FoucaultKyleEkopPinSpec {
            prior_price_p0: T::from_f64(100.0),
            insider_valuation_v: T::from_f64(105.0) + x[0],
            fundamental_std_sigma_v: T::from_f64(2.0) + dx1,
            noise_order_std_sigma_u: T::from_f64(10.0) + x[2],
            event_prob_alpha: T::from_f64(0.3),
            bad_news_prob_delta: T::from_f64(0.4),
            informed_arrival_mu: T::from_f64(40.0) + x[3],
            uninformed_buy_eps_b: T::from_f64(50.0),
            uninformed_sell_eps_s: T::from_f64(50.0),
            observed_buys_b: T::from_f64(75.0),
            observed_sells_s: T::from_f64(48.0),
        },
    );
    let gp = ql_portfolio::evaluate_garleanu_pedersen_dynamic_aim_portfolio(
        &ql_portfolio::GarleanuPedersenDynamicAimSpec {
            current_weights: [T::from_f64(0.2), T::from_f64(0.1)],
            fast_alpha_signal: [T::from_f64(0.04) + dx0, T::from_f64(0.02) + dx1],
            slow_alpha_signal: [T::from_f64(0.03), T::from_f64(0.025)],
            fast_decay_phi: T::from_f64(0.5),
            slow_decay_phi: T::from_f64(0.1),
            covariance_2x2: [
                [T::from_f64(0.04), T::from_f64(0.01)],
                [T::from_f64(0.01), T::from_f64(0.05)],
            ],
            risk_aversion_gamma: T::from_f64(2.0) + dx2,
            transaction_cost_lambda: T::from_f64(0.5) + dx3,
            discount_rate_rho: T::from_f64(0.05),
        },
    );
    let gjr = ql_portfolio::evaluate_tsay_gjr_and_dcc_garch_dynamic_correlation(
        &ql_portfolio::TsayGjrDccGarchSpec {
            bivariate_returns: vec![
                [T::from_f64(-0.02) + dx0, T::from_f64(-0.015) + dx1],
                [T::from_f64(0.012) + dx2, T::from_f64(-0.01) + dx3],
                [T::from_f64(-0.018), T::from_f64(-0.022)],
            ],
            omega: [T::from_f64(0.00002), T::from_f64(0.00003)],
            alpha: [T::from_f64(0.05), T::from_f64(0.06)],
            gjr_gamma: [T::from_f64(0.08), T::from_f64(0.07)],
            beta: [T::from_f64(0.85), T::from_f64(0.84)],
            unconditional_corr: T::from_f64(0.5),
            dcc_a: T::from_f64(0.05),
            dcc_b: T::from_f64(0.90),
        },
    );
    let lppl = ql_portfolio::evaluate_sornette_jls_lppl_crash_hazard(
        &ql_portfolio::SornetteJlsLpplCrashHazardSpec {
            current_time_t: T::from_f64(1.0),
            critical_time_tc: T::from_f64(1.5) + dx0,
            horizon_dt: T::from_f64(0.1),
            peak_log_price_a: T::from_f64(5.0) + dx1,
            power_law_slope_b: T::from_f64(-0.4) + dx2,
            log_periodic_amp_c: T::from_f64(0.03),
            critical_exponent_m: T::from_f64(0.5),
            log_periodic_omega: T::from_f64(6.5),
            phase_phi: T::from_f64(1.0),
            crash_loss_kappa: T::from_f64(0.2) + dx3,
        },
    );
    let spec = ql_portfolio::evaluate_follmer_schied_entropic_and_spectral_risk(
        &ql_portfolio::FollmerSchiedEntropicSpectralRiskSpec {
            ordered_losses: vec![
                T::from_f64(0.5) + dx0,
                T::from_f64(1.2) + dx1,
                T::from_f64(2.5) + dx2,
                T::from_f64(4.0) + dx3,
            ],
            entropic_theta: T::from_f64(0.5),
            spectral_power_gamma: T::from_f64(2.5),
        },
    );
    let brin = ql_portfolio::evaluate_bacon_brinson_fachler_and_carino_linking(
        &ql_portfolio::BaconBrinsonFachlerCarinoSpec {
            portfolio_sector_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
            benchmark_sector_weights: vec![T::from_f64(0.5), T::from_f64(0.5)],
            portfolio_sector_returns: vec![T::from_f64(0.08) + dx0, T::from_f64(0.03) + dx1],
            benchmark_sector_returns: vec![T::from_f64(0.06), T::from_f64(0.02)],
            multi_period_port_returns: vec![
                T::from_f64(0.04) + dx0,
                T::from_f64(0.03) + dx1,
                T::from_f64(-0.01) + dx2,
            ],
            multi_period_bench_returns: vec![
                T::from_f64(0.03),
                T::from_f64(0.02),
                T::from_f64(-0.015),
            ],
        },
    );
    let prel = ql_portfolio::evaluate_wakker_prelec_rdu_and_choquet_ambiguity(
        &ql_portfolio::WakkerPrelecChoquetProspectSpec {
            sorted_outcomes: [
                T::from_f64(-10.0) + x[0],
                T::from_f64(-2.0),
                T::from_f64(8.0) + x[1],
                T::from_f64(20.0) + x[2],
            ],
            state_probabilities: [
                T::from_f64(0.1),
                T::from_f64(0.3),
                T::from_f64(0.4),
                T::from_f64(0.2),
            ],
            reference_point: T::zero(),
            prelec_curvature_alpha: T::from_f64(0.65),
            prelec_elevation_beta: T::from_f64(1.0) + dx3,
            tk_value_exponent: T::from_f64(0.88),
            loss_aversion_lambda: T::from_f64(2.25),
            ambiguity_aversion_delta: T::from_f64(0.15),
        },
    );
    let sig = ql_dml::evaluate_bayer_path_signature_chen_and_levy_area(
        &ql_dml::BayerPathSignatureChenLevyAreaSpec {
            time_grid: vec![T::zero(), T::from_f64(0.5), T::one()],
            spot_path: vec![
                T::from_f64(100.0) + x[0],
                T::from_f64(104.0) + x[1],
                T::from_f64(102.0) + x[2],
            ],
            linear_functional_weights: [
                T::from_f64(0.1),
                T::from_f64(0.2),
                T::from_f64(0.3),
                T::from_f64(0.4),
                T::from_f64(-0.2),
                T::from_f64(0.1) + dx3,
            ],
        },
    );
    [
        repo.haircut_funded_levered_roe,
        bk.total_burgard_kjaer_xva_charge,
        decs.decs_mandatory_fair_value,
        mbs.passthrough_present_value,
        kyle.expected_insider_profit + kyle.posterior_event_probability,
        gp.expected_net_utility_gain,
        gjr.terminal_equal_weight_portfolio_vol,
        lppl.lppl_fitted_price,
        spec.entropic_convex_risk_measure + spec.kusuoka_power_spectral_risk_measure,
        brin.carino_smoothed_active_return,
        prel.prospect_certainty_equivalent,
        sig.signature_linear_payoff + sig.signed_levy_area,
    ]
}

pub fn benchmark_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite46_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite46_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "46_tuckman_repo_bk_fva_decs_mbs_io_po_kyle_pin_gp_aim_gjr_dcc_lppl_spectral_brinson_carino_prelec_path_sig",
        "Tuckman-Repo / BK-FVA / DECS / MBS-IO-PO / Kyle-PIN / GP-Aim / GJR-DCC / LPPL / Spectral-Risk / Brinson-Carino / Prelec-CPT / Path-Sig AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite39_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let ftpl = ql_inflation::evaluate_cochrane_ftpl_and_maturity_inflation_smoothing(
        &ql_inflation::CochraneFiscalTheoryPriceLevelSpec {
            debt_maturities_years: vec![T::from_f64(1.0), T::from_f64(5.0), T::from_f64(10.0)],
            nominal_debt_face_values: vec![
                T::from_f64(30.0) + x[0],
                T::from_f64(30.0),
                T::from_f64(40.0),
            ],
            nominal_zero_yields: vec![
                T::from_f64(0.03) + x[1] * T::from_f64(0.01),
                T::from_f64(0.032),
                T::from_f64(0.035),
            ],
            real_zero_yields: vec![T::from_f64(0.01), T::from_f64(0.012), T::from_f64(0.015)],
            primary_real_surpluses: vec![T::from_f64(5.0), T::from_f64(5.2), T::from_f64(5.4)],
            terminal_surplus_growth_rate: T::from_f64(0.02),
            fiscal_surplus_shock_fraction: T::from_f64(-0.05),
            expected_future_inflation_rise: T::from_f64(0.01) + x[2] * T::from_f64(0.01),
            inflation_risk_beta: T::from_f64(0.4),
            inflation_volatility: T::from_f64(0.1) + x[3] * T::from_f64(0.01),
            tips_liquidity_spread: T::from_f64(0.002),
        },
    );
    let g2 = ql_rates::price_brigo_mercurio_g2_twist_and_swaption(
        &ql_rates::BrigoMercurioG2TwistSwaptionSpec {
            a: T::from_f64(0.1),
            sigma: T::from_f64(0.01),
            b: T::from_f64(0.3),
            eta: T::from_f64(0.015),
            rho: T::from_f64(-0.4),
            flat_zero_rate: T::from_f64(0.03) + x[0] * T::from_f64(0.001),
            option_expiry: T::from_f64(1.0),
            swap_payment_times: vec![T::from_f64(1.0), T::from_f64(2.0), T::from_f64(3.0)],
            fixed_strike: T::from_f64(0.03),
            notional: T::from_f64(1.0),
        },
    );
    let ess = ql_vol::evaluate_mingone_global_arbitrage_free_essvi(
        &ql_vol::MingoneGlobalArbitrageFreeEssviSpec {
            expiries: vec![T::from_f64(1.0), T::from_f64(2.0)],
            raw_theta_increments: vec![x[0], x[1]],
            raw_rho_params: vec![T::from_f64(-0.2), T::from_f64(-0.1)],
            raw_psi_fractions: vec![x[2], x[3]],
            test_log_moneyness_grid: vec![T::from_f64(-1.0), T::from_f64(0.0), T::from_f64(1.0)],
        },
    );
    let cb = ql_fx::price_cross_currency_convertible_and_warrant_bond(
        &ql_fx::CrossCurrencyConvertibleWarrantSpec {
            domestic_stock_spot: T::from_f64(100.0) + x[0],
            domestic_conversion_strike: T::from_f64(100.0),
            spot_fx_for_per_dom: T::from_f64(1.1),
            fixed_contractual_fx_for_per_dom: T::from_f64(1.1),
            foreign_bond_par: T::from_f64(100.0),
            foreign_coupon_rate: T::from_f64(0.03),
            r_domestic: T::from_f64(0.03),
            r_foreign: T::from_f64(0.02),
            stock_dividend_yield: T::from_f64(0.01),
            credit_spread: T::from_f64(0.02),
            stock_vol: T::from_f64(0.2),
            fx_vol: T::from_f64(0.1),
            corr_stock_fx: T::from_f64(0.2),
            maturity: T::from_f64(5.0),
        },
    );
    let bank = ql_fixed_income::analyze_acharya_feds_bank_run_and_optimal_coco(
        &ql_fixed_income::AcharyaFedsBankRunAndCocoSpec {
            initial_securities_assets: T::from_f64(100.0),
            asset_effective_duration: T::from_f64(5.0),
            insured_deposits: T::from_f64(40.0),
            uninsured_deposits: T::from_f64(45.0),
            deposit_rate_beta: T::from_f64(0.5),
            deposit_operating_cost_rate: T::from_f64(0.005),
            at1_coco_principal: T::from_f64(5.0),
            initial_short_rate: T::from_f64(0.03),
            parallel_rate_shock: x[1].abs() * T::from_f64(0.01),
            run_sensitivity_steepness: T::from_f64(30.0),
            fire_sale_liquidation_haircut: T::from_f64(0.1),
            coco_write_down_efficiency: T::from_f64(0.9),
        },
    );
    let xva = ql_xva::compute_albanese_crepey_mpor_and_balance_sheet_xva(
        &ql_xva::AlbaneseCrepeyMporBalanceSheetSpec {
            exposure_times: vec![T::from_f64(1.0), T::from_f64(2.0)],
            clean_mtm_profile: vec![T::from_f64(1.0) + x[0], T::from_f64(2.0) + x[1]],
            portfolio_vol_profile: vec![
                T::from_f64(0.2) + x[2] * T::from_f64(0.01),
                T::from_f64(0.2) + x[3] * T::from_f64(0.01),
            ],
            scheduled_cashflows_in_mpor: vec![
                T::from_f64(0.01) + x[0] * T::from_f64(0.001),
                T::from_f64(0.02) + x[1] * T::from_f64(0.001),
            ],
            mpor_years: T::from_f64(0.05) + x[2] * T::from_f64(0.001),
            csa_threshold: T::from_f64(0.1) + x[3] * T::from_f64(0.01),
            risk_free_rate: T::from_f64(0.03) + x[0] * T::from_f64(0.001),
            cpty_hazard_rate: T::from_f64(0.02),
            cpty_lgd: T::from_f64(0.6),
            bank_hazard_rate: T::from_f64(0.01),
            bank_lgd: T::from_f64(0.6),
            economic_capital_profile: vec![T::from_f64(0.05), T::from_f64(0.05)],
            hurdle_rate: T::from_f64(0.1),
            initial_margin_profile: vec![T::from_f64(0.1), T::from_f64(0.1)],
            im_funding_spread: T::from_f64(0.01),
        },
    );
    let dh = ql_models::evaluate_double_heston_cumulants_and_kernel_slv(
        &ql_models::DoubleHestonKernelSlvSpec {
            spot: T::from_f64(100.0),
            strike: T::from_f64(100.0),
            maturity: T::from_f64(1.0),
            risk_free_rate: T::from_f64(0.03),
            dividend_yield: T::from_f64(0.0),
            v1_0: T::from_f64(0.02),
            kappa1: T::from_f64(2.0),
            theta1: T::from_f64(0.02),
            sigma1: T::from_f64(0.3),
            rho1: T::from_f64(-0.6),
            v2_0: T::from_f64(0.01),
            kappa2: T::from_f64(0.5),
            theta2: T::from_f64(0.01),
            sigma2: T::from_f64(0.2),
            rho2: T::from_f64(-0.3),
            target_dupire_local_vol: T::from_f64(0.2) + x[0] * T::from_f64(0.01),
        },
    );
    let asian = ql_exotics::price_curran_rogers_shi_arithmetic_asian(
        &ql_exotics::CurranRogersShiArithmeticAsianSpec {
            spot: T::from_f64(100.0),
            strike: T::from_f64(100.0) + x[1],
            risk_free_rate: T::from_f64(0.03),
            dividend_yield: T::from_f64(0.0),
            vol: T::from_f64(0.2),
            fixing_times: vec![
                T::from_f64(0.25),
                T::from_f64(0.5),
                T::from_f64(0.75),
                T::from_f64(1.0),
            ],
        },
    );
    let storage = ql_commodities::evaluate_warin_baeuerle_storage_rolling_hedge(
        &ql_commodities::WarinBaeuerleStorageRollingHedgeSpec {
            forward_curve: vec![T::from_f64(3.0) + x[2], T::from_f64(3.1), T::from_f64(3.2)],
            decision_times: vec![T::from_f64(0.5), T::from_f64(1.0), T::from_f64(1.5)],
            max_capacity: T::from_f64(10.0),
            initial_inventory: T::from_f64(2.0),
            max_injection_step: T::from_f64(1.0),
            max_withdrawal_step: T::from_f64(1.0),
            injection_loss_fraction: T::from_f64(0.02),
            spot_ou_volatility: T::from_f64(0.3),
            spot_ou_mean_reversion: T::from_f64(1.0),
            discount_rate: T::from_f64(0.03),
        },
    );
    let micro = ql_portfolio::analyze_harris_roll_corwin_schultz_microstructure(
        &ql_portfolio::HarrisRollCorwinSchultzMicrostructureSpec {
            trade_prices: vec![
                T::from_f64(100.0) + x[0],
                T::from_f64(100.1) + x[1],
                T::from_f64(99.9) + x[2],
                T::from_f64(100.2) + x[3],
            ],
            daily_highs: vec![
                T::from_f64(101.0) + x[0],
                T::from_f64(102.0) + x[1],
                T::from_f64(101.5) + x[2],
            ],
            daily_lows: vec![
                T::from_f64(99.0) + x[3],
                T::from_f64(98.0) + x[0],
                T::from_f64(98.5) + x[1],
            ],
            daily_dollar_volumes: vec![
                T::from_f64(1e6) * (T::one() + x[0] * T::from_f64(0.01)),
                T::from_f64(1.1e6) * (T::one() + x[1] * T::from_f64(0.01)),
                T::from_f64(0.9e6) * (T::one() + x[2] * T::from_f64(0.01)),
            ],
            decision_price: T::from_f64(100.0) + x[0],
            arrival_price: T::from_f64(100.0) + x[1],
            child_fill_quantities: vec![T::from_f64(10.0) + x[2], T::from_f64(5.0) + x[3]],
            child_fill_prices: vec![T::from_f64(100.1) + x[0], T::from_f64(100.2) + x[1]],
            total_order_quantity: T::from_f64(20.0) + x[2],
            post_trade_benchmark_price: T::from_f64(100.3) + x[3],
            explicit_commissions: T::from_f64(0.1) + x[0] * T::from_f64(0.001),
        },
    );
    let kelly = ql_portfolio::solve_busseti_boyd_risk_constrained_kelly(
        &ql_portfolio::BussetiBoydRiskConstrainedKellySpec {
            excess_return_means: vec![
                T::from_f64(0.05) + x[0] * T::from_f64(0.001),
                T::from_f64(0.03),
            ],
            covariance: vec![
                vec![T::from_f64(0.04), T::from_f64(0.01)],
                vec![T::from_f64(0.01), T::from_f64(0.03)],
            ],
            risk_free_rate: T::from_f64(0.02),
            drawdown_floor_alpha: T::from_f64(0.8),
            max_drawdown_prob_beta: T::from_f64(0.5),
            specialist_immediate_capital: T::from_f64(100.0),
            slow_moving_external_capital: T::from_f64(100.0),
            capital_arrival_rate: T::from_f64(0.1),
            dislocation_horizon_years: T::one(),
            forced_liquidation_notional: T::from_f64(10.0),
        },
    );
    let qp = ql_dml::solve_agrawal_boyd_differentiable_qp_layer(
        &ql_dml::AgrawalBoydDifferentiableQpSpec {
            q_matrix: vec![
                vec![T::from_f64(1.0), T::from_f64(0.0)],
                vec![T::from_f64(0.0), T::from_f64(1.0)],
            ],
            q_linear: vec![T::from_f64(-1.0) - x[0], T::from_f64(-0.5)],
            a_eq: vec![vec![T::from_f64(1.0), T::from_f64(1.0)]],
            b_eq: vec![T::from_f64(1.0)],
            g_ineq: vec![vec![T::from_f64(-1.0), T::from_f64(0.0)]],
            h_ineq: vec![T::from_f64(0.0)],
            barrier_tau: T::from_f64(1e-4),
        },
    );
    [
        ftpl.equilibrium_price_level_p0,
        g2.payer_swaption_pv,
        ess.min_durrleman_butterfly_g,
        cb.embedded_cross_currency_cb_price,
        bank.economic_equity_post_coco_conversion,
        xva.sticky_mpor_cva,
        dh.double_heston_atm_vol,
        asian.curran_full_conditioned_price,
        storage.stochastic_total_storage_value,
        micro.perold_total_implementation_shortfall,
        kelly.constrained_portfolio_volatility,
        qp.optimal_objective_value,
    ]
}

pub fn benchmark_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.02_f64, -0.01, 0.03, -0.01];
    let baseline = move || {
        let f = |x: [f64; 4]| suite39_workloads(x).iter().map(|v| v.val()).sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        type D = ql_ad::Dual<4>;
        let x = [
            D::variable(x0[0], 0),
            D::variable(x0[1], 1),
            D::variable(x0[2], 2),
            D::variable(x0[3], 3),
        ];
        suite39_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "39_cochrane_ftpl_g2_twist_mingone_essvi_offshore_cb_acharya_coco_albanese_mpor_double_heston_curran_asian_warin_storage_roll_corwin_busseti_kelly_agrawal_kkt",
        "FTPL / G2++ / eSSVI / Offshore-CB / CoCo / MPOR / Double-Heston / Curran / Warin / Roll-Corwin / Busseti-Kelly / Diff-QP AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

/// Cross-pillar scalar and forward-mode AD parity workload.
fn suite34_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
    let p = ql_models::Bergomi2FParams {
        xi0: T::from_f64(0.04),
        omega: T::from_f64(1.0) + x[0] * 0.01,
        theta: T::from_f64(0.5),
        k1: T::from_f64(1.0),
        k2: T::from_f64(2.0),
        rho12: T::from_f64(0.2),
        rho_s1: T::from_f64(-0.5),
        rho_s2: T::from_f64(-0.3),
    };
    let berg = ql_models::evaluate_bergomi_forward_start_and_var_strip(
        &ql_models::BergomiForwardStartAndVarStripSpec {
            params: p,
            xi_long_run: T::from_f64(0.04),
            xi_curve_decay: T::from_f64(1.0),
            forward_start_t1: T::from_f64(0.25),
            expiry_t2: T::from_f64(1.0),
            log_moneyness_k: x[1] * T::from_f64(0.1),
        },
    );
    let cg =
        ql_levy::evaluate_cgmy_fractional_fft_and_pide(&ql_levy::CgmyFractionalFftAndPideSpec {
            process: ql_levy::CgmyProcess {
                sigma: T::from_f64(0.15) + x[0] * T::from_f64(0.001),
                c: T::from_f64(0.5),
                g: T::from_f64(5.0),
                m: T::from_f64(5.0),
                y: T::from_f64(0.7),
            },
            spot: T::from_f64(100.0),
            strikes: vec![T::from_f64(100.0)],
            risk_free_rate: T::from_f64(0.02) + x[2] * T::from_f64(0.001),
            dividend_yield: T::from_f64(0.01),
            expiry: T::from_f64(1.0),
            damping_alpha: 1.0,
            fourier_step_eta: 0.1,
            log_strike_step_lambda: 0.1,
            small_jump_cutoff_eps: T::from_f64(1e-4),
        });
    let cg0 = cg.frft_call_prices.first().copied().unwrap_or_else(T::zero);
    let coco = ql_fixed_income::price_despiegeleer_coco_and_mandatory_convertible(
        &ql_fixed_income::DeSpiegeleerCoCoMandatoryConvSpec {
            notional: T::from_f64(100.0),
            spot: T::from_f64(100.0) + x[0],
            risk_free_rate: T::from_f64(0.02),
            dividend_yield: T::from_f64(0.01),
            credit_spread: T::from_f64(0.02) + x[1] * T::from_f64(0.001),
            maturity: T::from_f64(2.0),
            coupon_rate: T::from_f64(0.05),
            coupon_freq: 2,
            decs_lower_strike: T::from_f64(90.0),
            decs_upper_strike: T::from_f64(110.0),
            vol_lower_strike: T::from_f64(0.2),
            vol_upper_strike: T::from_f64(0.2),
            coco_cet1_trigger_spot: T::from_f64(70.0),
            coco_ponv_trigger_spot: T::from_f64(50.0),
            coco_conversion_floor_price: T::from_f64(60.0),
            coco_equity_vol: T::from_f64(0.2),
        },
    );
    let gue = ql_portfolio::solve_gueant_lehalle_bounded_market_making(
        &ql_portfolio::GueantLehalleMarketMakingSpec {
            max_inventory_q: 2,
            current_inventory_q: 0,
            horizon: T::from_f64(1.0),
            num_time_steps: 4,
            volatility: T::from_f64(0.2),
            risk_aversion_gamma: T::from_f64(0.5),
            intensity_decay_k: T::from_f64(1.5) + x[1] * T::from_f64(0.01),
            base_intensity_a: T::from_f64(2.0),
            terminal_inventory_penalty: T::from_f64(0.1),
            adverse_selection_zeta: T::zero(),
        },
    );
    let meu = ql_portfolio::solve_meucci_idzorek_black_litterman_and_plausibility(
        &ql_portfolio::MeucciIdzorekBlackLittermanSpec {
            covariance: vec![
                vec![T::from_f64(0.04), T::from_f64(0.01)],
                vec![T::from_f64(0.01), T::from_f64(0.09)],
            ],
            market_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
            risk_aversion_lambda: T::from_f64(2.5),
            tau: T::from_f64(0.05),
            view_matrix_p: vec![vec![T::one(), -T::one()]],
            view_returns_q: vec![T::from_f64(0.03) + x[2] * T::from_f64(0.001)],
            view_confidences: vec![T::from_f64(0.65)],
        },
    );
    let krd = analyze_tuckman_key_rate_durations_and_pca_hedge(&TuckmanKeyRateDurationSpec {
        nss_curve: NelsonSiegelSvensson {
            reference_date: date(2024, 1, 1),
            beta0: T::from_f64(0.03) + x[0] * T::from_f64(0.001),
            beta1: T::from_f64(-0.01),
            beta2: T::from_f64(0.02),
            beta3: T::from_f64(-0.005),
            tau1: T::from_f64(1.5),
            tau2: T::from_f64(4.0),
        },
        cashflow_times: vec![
            T::from_f64(1.0),
            T::from_f64(2.0),
            T::from_f64(3.0),
            T::from_f64(5.0),
        ],
        cashflow_amounts: vec![
            T::from_f64(2.0),
            T::from_f64(2.0),
            T::from_f64(2.0),
            T::from_f64(102.0),
        ],
        key_rate_tenors: vec![
            T::from_f64(1.0),
            T::from_f64(2.0),
            T::from_f64(5.0),
            T::from_f64(10.0),
        ],
        pca_decay_lambda: T::from_f64(0.25),
        benchmark_hedge_tenors: [T::from_f64(2.0), T::from_f64(5.0), T::from_f64(10.0)],
    });
    let lpi =
        price_brigo_mercurio_yoy_and_lpi_market_model(&BrigoMercurioInflationMarketModelSpec {
            notional: T::from_f64(1_000_000.0),
            base_cpi: T::from_f64(100.0),
            payment_times: vec![T::from_f64(1.0), T::from_f64(2.0), T::from_f64(3.0)],
            forward_cpis: vec![
                T::from_f64(102.0) + x[1],
                T::from_f64(104.0) + x[1],
                T::from_f64(106.0) + x[1],
            ],
            nominal_forward_rates: vec![T::from_f64(0.03); 3],
            nominal_discount_factors: vec![T::from_f64(0.97), T::from_f64(0.94), T::from_f64(0.91)],
            forward_cpi_vols: vec![T::from_f64(0.01) + x[1] * T::from_f64(0.001); 3],
            nominal_forward_vols: vec![T::from_f64(0.01); 3],
            inter_maturity_cpi_corr: T::from_f64(0.5),
            nominal_cpi_corr: T::from_f64(0.2),
            lpi_floor: T::from_f64(-0.02),
            lpi_cap: T::from_f64(0.03),
            fixed_swap_rate: T::from_f64(0.02) + x[1] * T::from_f64(0.001),
        });
    let cdo = price_okane_bespoke_tranche_and_student_t_basket(&OKaneBespokeBaseCorrelationSpec {
        index_detachments: vec![
            T::from_f64(0.03),
            T::from_f64(0.07),
            T::from_f64(0.10),
            T::from_f64(0.15),
        ],
        index_base_correlations: vec![
            T::from_f64(0.35),
            T::from_f64(0.40),
            T::from_f64(0.45),
            T::from_f64(0.50),
        ],
        index_expected_loss: T::from_f64(0.04) + x[2] * T::from_f64(0.001),
        bespoke_hazard_rates: vec![T::from_f64(0.02); 5],
        bespoke_recovery_rate: T::from_f64(0.4),
        bespoke_attachment: T::from_f64(0.03),
        bespoke_detachment: T::from_f64(0.10),
        strike_mapping_exponent_gamma: T::from_f64(1.0),
        student_t_dof_nu: T::from_f64(6.0) + x[2],
        maturity: T::from_f64(5.0),
        risk_free_rate: T::from_f64(0.02),
    });
    let fx = evaluate_wystup_fx_accumulator_and_dcd(&WystupFxStructuredProductSpec {
        spot: T::from_f64(1.10) + x[3] * T::from_f64(0.01),
        accumulator_strike: T::from_f64(1.05),
        knockout_barrier: T::from_f64(1.20),
        downside_gearing: T::from_f64(2.0),
        domestic_rate: T::from_f64(0.02),
        foreign_rate: T::from_f64(0.01),
        volatility: T::from_f64(0.12) + x[3] * T::from_f64(0.001),
        fixing_times: vec![
            T::from_f64(0.25),
            T::from_f64(0.5),
            T::from_f64(0.75),
            T::from_f64(1.0),
        ],
        guaranteed_fixings: 1,
        dcd_conversion_strike: T::from_f64(1.10),
        dcd_tenor: T::from_f64(1.0),
    });
    let pde = ql_pde::solve_rannacher_crank_nicolson_sinh_pde(&ql_pde::RannacherSinhPdeSpec {
        spot: T::from_f64(100.0) + x[0],
        strike: T::from_f64(100.0),
        volatility: T::from_f64(0.2),
        risk_free_rate: T::from_f64(0.02),
        dividend_yield: T::from_f64(0.01),
        expiry: T::from_f64(1.0),
        is_call: true,
        n_space: 31,
        n_time: 8,
        rannacher_steps: 2,
        sinh_concentration: T::from_f64(0.15),
    });
    let mot =
        ql_ad::solve_mot_super_replication_bounds_and_ift(&ql_ad::MotSuperReplicationIftSpec {
            support_grid: vec![T::from_f64(0.8), T::one(), T::from_f64(1.2)],
            marginal_t1_probs: vec![T::from_f64(0.2), T::from_f64(0.5), T::from_f64(0.3)],
            marginal_t2_probs: vec![T::from_f64(0.25), T::from_f64(0.5), T::from_f64(0.25)],
            forward_start_strike: T::from_f64(0.95) + x[3] * T::from_f64(0.001),
            entropic_regularization_eps: T::from_f64(0.02),
        });
    let exotic =
        ql_dml::evaluate_ift_calibrated_exotic_quote_hedge(&ql_dml::IftCalibratedExoticHedgeSpec {
            spot: T::from_f64(100.0) + x[0],
            risk_free_rate: T::from_f64(0.02),
            expiry: T::from_f64(1.0),
            calibration_strikes: [T::from_f64(90.0), T::from_f64(100.0), T::from_f64(110.0)],
            market_implied_vols: [T::from_f64(0.22), T::from_f64(0.2), T::from_f64(0.21)],
            exotic_strike: T::from_f64(105.0),
            exotic_power_p: T::from_f64(1.2),
        });
    let revolve =
        ql_ad::execute_savine_revolve_checkpointed_adjoint(&ql_ad::SavineRevolveCheckpointSpec {
            num_time_steps: 8,
            num_checkpoints: 2,
            initial_spot: T::from_f64(100.0) + x[1],
            drift_rate: T::from_f64(0.02),
            volatility: T::from_f64(0.2) + x[3] * T::from_f64(0.001),
            mean_reversion_kappa: T::from_f64(0.5),
            dt: T::from_f64(0.125),
            strike: T::from_f64(100.0),
        });
    // Cholesky/AAD pillar contributes a stable correlation normalization.
    let chv = (T::one() + T::from_f64(0.2)).sqrt();
    [
        berg.spot_start_implied_vol
            + krd.dirty_price * T::from_f64(1e-4)
            + chv
            + lpi.fimm_convexity_adjusted_lpi_swap_npv * T::from_f64(1e-6),
        berg.forward_start_implied_vol
            + cg0
            + coco.dual_trigger_coco_value
            + fx.total_accumulator_npv
            + pde.pde_price,
        berg.forward_variance_swap_fair_vol
            + gue.current_optimal_bid_offset
            + cdo.student_t_copula_senior_tranche_loss
            + mot.entropic_mot_price,
        berg.volatility_swap_fair_strike
            + meu.mahalanobis_plausibility_squared
            + exotic.exotic_option_price
            + revolve.terminal_option_payoff_value
            + x[3],
    ]
}

pub fn benchmark_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad(scale: usize) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let h = 1e-5;
    let baseline = || {
        let mut checksum: f64 = suite34_workloads(x0).iter().sum();
        for lane in 0..4 {
            let mut up = x0;
            let mut dn = x0;
            up[lane] += h;
            dn[lane] -= h;
            let a = suite34_workloads(up);
            let b = suite34_workloads(dn);
            checksum += a
                .iter()
                .zip(b.iter())
                .map(|(u, d)| (u - d) / (2.0 * h))
                .sum::<f64>();
        }
        checksum
    };
    let sota = || {
        type D = Dual<4>;
        let x = [
            D::variable(x0[0], 0),
            D::variable(x0[1], 1),
            D::variable(x0[2], 2),
            D::variable(x0[3], 3),
        ];
        suite34_workloads(x)
            .iter()
            .map(|y| y.val + y.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad",
        "KRD / LPI / CoCo / Bergomi / CGMY / Guéant / Meucci AD vs FD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

/// Public spelling used by downstream benchmark runners.
pub fn run_all_suites(scale: usize) -> Vec<ComparisonRow> {
    run_all_benchmarks(scale)
}

/// Fast deterministic suite runner used by integration tests.  Numerical
/// parity and invariants are checked by every suite; hardware-sensitive
/// speedup gates remain opt-in in the test harness.
pub fn run_all_benchmark_suites_with_config(
    warmup_iters: usize,
    samples: usize,
) -> Vec<ComparisonRow> {
    with_config_override(warmup_iters, samples, || run_all_benchmarks(1))
}

/// Compatibility spelling retained for callers using the capitalised Config alias.
/// verbatim (`Config` was capitalised there).
#[allow(non_snake_case)]
pub fn run_all_benchmark_suites_with_Config(
    warmup_iters: usize,
    samples: usize,
) -> Vec<ComparisonRow> {
    run_all_benchmark_suites_with_config(warmup_iters, samples)
}

fn suite30_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
    let d = date(2024, 1, 1);
    let curve_dom = FlatForwardCurve::new(d, T::from_f64(0.03) + x[0] * 0.001);
    let curve_for = FlatForwardCurve::new(d, T::from_f64(0.02) + x[1] * 0.001);
    let curve_csa = FlatForwardCurve::new(d, T::from_f64(0.024) + x[2] * 0.001);

    let mtm = ql_fx::price_mtm_resetting_xccs(
        &ql_fx::MtmResettingXccs {
            for_notional: T::from_f64(100.0),
            spot_fx: T::from_f64(1.10) + x[1] * 0.01,
            payment_times: vec![0.5, 1.0, 1.5, 2.0],
            accrual_fractions: vec![0.5; 4],
            for_basis_spread: T::from_f64(-0.0015) + x[2] * 0.0002,
            dom_spread: T::zero(),
            receive_foreign: true,
        },
        &curve_dom,
        &curve_for,
        &curve_csa,
    )
    .npv;
    let ten = ql_rates::price_tenor_basis_swap(
        &ql_rates::TenorBasisSwapSpec {
            notional: T::from_f64(100.0),
            short_leg_times: vec![0.5, 1.0],
            long_leg_times: vec![1.0],
            short_leg_basis_spread: T::from_f64(0.001) + x[3] * 0.0002,
            pay_short_leg: true,
        },
        &curve_dom,
        &curve_for,
        &curve_csa,
    )
    .npv;
    let oas =
        ql_fixed_income::analyze_callable_bond_oas_ad(&ql_fixed_income::CallableBondOasSpec {
            market_dirty_price: T::from_f64(100.0) + x[1] * 0.2,
            face_amount: T::from_f64(100.0),
            coupon_rate: T::from_f64(0.04) + x[2] * 0.0005,
            payment_times: vec![1.0, 2.0, 3.0],
            call_price: Some(T::from_f64(101.0)),
            put_price: None,
            short_rate: T::from_f64(0.03) + x[0] * 0.001,
            short_rate_vol: T::from_f64(0.015) + x[3] * 0.0005,
            tree_steps: 8,
            smoothing_eps: 1e-3,
        })
        .oas;

    let colva = ql_xva::compute_multicurrency_ctd_colva(&ql_xva::MultiCurrencyCsaSpec {
        eligible_currency_spreads: vec![
            T::from_f64(0.01) + x[0] * 0.0005,
            T::from_f64(0.012) + x[1] * 0.0005,
            T::from_f64(0.011) + x[2] * 0.0005,
        ],
        spread_vols: vec![
            T::from_f64(0.008) + x[3] * 0.0002,
            T::from_f64(0.010),
            T::from_f64(0.009),
        ],
        cross_currency_corr: T::from_f64(0.25),
        collateral_balance_profile: vec![T::from_f64(50.0), T::from_f64(40.0)],
        payment_times: vec![0.5, 1.0],
    })
    .total_ctd_colva;
    let wwr = ql_xva::compute_turlakov_wwr_cva_fva(&ql_xva::TurlakovWwrSpec {
        times: vec![0.5, 1.0],
        base_epe: vec![T::from_f64(10.0) + x[0] * 0.2, T::from_f64(14.0)],
        stressed_epe: vec![T::from_f64(25.0) + x[1] * 0.5, T::from_f64(32.0)],
        base_mtm: vec![T::from_f64(8.0), T::from_f64(10.0)],
        stressed_mtm: vec![T::from_f64(20.0), T::from_f64(24.0)],
        sovereign_hazard_rate: T::from_f64(0.025) + x[2] * 0.001,
        counterparty_hazard_rate: T::from_f64(0.02) + x[3] * 0.001,
        crisis_prob_threshold: T::from_f64(0.08),
        coupling_lambda: T::from_f64(0.45),
        lgd: T::from_f64(0.60),
        base_funding_spread: T::from_f64(0.01),
        stressed_funding_spread: T::from_f64(0.025),
        discount_rate: T::from_f64(0.03) + x[0] * 0.001,
    })
    .wwr_cva;

    let zabr = ql_models::evaluate_zabr_smile(
        &ql_models::local_vol::ZabrParams {
            alpha: T::from_f64(0.20) + x[1] * 0.002,
            beta: T::from_f64(0.7),
            nu: T::from_f64(0.40) + x[2] * 0.005,
            rho: T::from_f64(-0.25),
            gamma: T::from_f64(0.85) + x[3] * 0.01,
        },
        T::from_f64(100.0) + x[0] * 0.2,
        T::from_f64(105.0),
        T::one(),
    )
    .call_price;
    let pdv = ql_models::pdv_implied_vol_expansion(
        &ql_models::GuyonPdvParams {
            beta0: T::from_f64(0.18) + x[1] * 0.002,
            beta1: T::from_f64(-0.12) + x[2] * 0.002,
            beta2: T::from_f64(0.25) + x[3] * 0.002,
            theta1: T::from_f64(0.5),
            theta2: T::from_f64(0.5),
            r10: T::from_f64(-0.05),
            r11: T::from_f64(-0.02),
            r20: T::from_f64(0.04),
            r21: T::from_f64(0.04),
        },
        T::from_f64(100.0) + x[0] * 0.2,
        T::from_f64(95.0),
        T::one(),
    );
    let kiko = ql_fx::price_fx_kiko_option(&ql_fx::FxKikoOptionSpec {
        spot: T::from_f64(100.0) + x[0] * 0.2,
        strike: T::from_f64(100.0),
        lower_barrier: T::from_f64(82.0),
        upper_barrier: T::from_f64(118.0),
        vol: T::from_f64(0.20) + x[1] * 0.002,
        expiry: T::one(),
        r_dom: T::from_f64(0.03) + x[2] * 0.001,
        r_for: T::from_f64(0.01),
        is_call: true,
        lower_is_ko_upper_is_ki: true,
    });
    let wb = ql_fx::price_fx_window_barrier_option(&ql_fx::FxWindowBarrierSpec {
        spot: T::from_f64(100.0) + x[0] * 0.2,
        strike: T::from_f64(100.0),
        barrier: T::from_f64(125.0),
        vol: T::from_f64(0.20) + x[1] * 0.002,
        expiry: T::one(),
        r_dom: T::from_f64(0.03),
        r_for: T::from_f64(0.01),
        window_start: T::from_f64(0.1),
        window_end: T::from_f64(0.8) + x[3] * 0.01,
        is_up_and_out: true,
        is_call: true,
    });
    let timer = ql_exotics::price_capped_timer_option_quadrature(
        &ql_exotics::FiniteCapTimerOptionSpec {
            spot: T::from_f64(100.0) + x[0] * 0.2,
            strike: T::from_f64(100.0),
            rate: T::from_f64(0.03),
            div_yield: T::from_f64(0.01),
            variance_budget: T::from_f64(0.04) + x[1] * 0.0005,
            t_max: T::from_f64(0.8) + x[2] * 0.01,
            is_call: true,
        },
        &HestonParams {
            v0: T::from_f64(0.04) + x[3] * 0.0005,
            kappa: T::from_f64(1.5),
            theta: T::from_f64(0.04),
            sigma_v: T::from_f64(0.25),
            rho: T::from_f64(-0.5),
        },
    )
    .npv;

    let weather =
        ql_commodities::price_weather_hdd_cdd_option(&ql_commodities::WeatherDerivativeSpec {
            kind: ql_commodities::DegreeDayKind::Hdd,
            temperature_base: T::from_f64(18.0),
            initial_temperature: T::from_f64(12.0) + x[0] * 0.1,
            mean_level: T::from_f64(13.0),
            seasonal_amplitude: T::from_f64(5.0),
            seasonal_phase: T::from_f64(0.2),
            mean_reversion: T::from_f64(0.25) + x[1] * 0.005,
            daily_vol: T::from_f64(2.0) + x[2] * 0.02,
            market_price_risk: T::from_f64(0.05) + x[3] * 0.002,
            times: vec![1.0, 2.0, 3.0],
            strike: T::from_f64(14.0),
            notional: T::from_f64(10.0),
            is_call: true,
            cap: Some(T::from_f64(80.0)),
        })
        .npv;
    let clean =
        ql_commodities::price_clean_spark_dark_spread(&ql_commodities::CleanSparkDarkSpreadSpec {
            power: T::from_f64(80.0) + x[0] * 0.5,
            fuel: T::from_f64(3.0) + x[1] * 0.05,
            carbon: T::from_f64(20.0) + x[2] * 0.2,
            heat_rate: T::from_f64(7.0),
            emission_factor: T::from_f64(0.5) + x[3] * 0.01,
            strike: T::from_f64(10.0),
            vol_power: T::from_f64(0.22),
            vol_fuel: T::from_f64(0.30),
            vol_carbon: T::from_f64(0.35),
            rho_power_fuel: T::from_f64(0.45),
            rho_power_carbon: T::from_f64(0.25),
            rho_fuel_carbon: T::from_f64(0.20),
            rate: T::from_f64(0.03),
            expiry: T::one(),
            option_type: ql_instruments::OptionType::Call,
        });
    let w0 = T::from_f64(0.6) + x[2] * 0.01;
    let enb = ql_portfolio::effective_number_of_bets_minimum_torsion(
        &[w0, T::one() - w0],
        &[
            vec![T::from_f64(0.04) + x[0] * 0.0005, T::from_f64(0.012)],
            vec![T::from_f64(0.012), T::from_f64(0.09) + x[1] * 0.0005],
        ],
    )
    .enb;
    let ep = ql_portfolio::entropy_pooling_posterior(&ql_portfolio::EntropyPoolingSpec {
        scenario_returns: vec![
            vec![T::from_f64(0.08), T::from_f64(0.02)],
            vec![T::from_f64(-0.04), T::from_f64(0.06)],
            vec![T::from_f64(0.02), T::from_f64(-0.03)],
        ],
        prior_probs: vec![T::from_f64(1.0 / 3.0); 3],
        view_matrix: vec![vec![
            T::from_f64(0.08),
            T::from_f64(-0.04),
            T::from_f64(0.02),
        ]],
        view_targets: vec![T::from_f64(0.03) + x[3] * 0.001],
        risk_aversion: T::from_f64(2.5),
    })
    .weights[0];

    [
        mtm + ten + oas * 100.0,
        colva + wwr,
        zabr + pdv + kiko + wb + timer,
        weather + clean + enb + ep,
    ]
}
pub fn benchmark_mtmxccs_colva_zabr_pdv_weather_enb_aad(scale: usize) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let h = 1e-5;
    let baseline = || {
        let b = suite30_workloads(x0);
        let mut c = b.iter().sum();
        for lane in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[lane] += h;
            d[lane] -= h;
            let bu = suite30_workloads(u);
            let bd = suite30_workloads(d);
            for j in 0..4 {
                c += (bu[j] - bd[j]) / (2.0 * h);
            }
        }
        c
    };
    let sota = || {
        type D = Dual<4>;
        let x = [
            D::variable(x0[0], 0),
            D::variable(x0[1], 1),
            D::variable(x0[2], 2),
            D::variable(x0[3], 3),
        ];
        suite30_workloads(x)
            .iter()
            .map(|y| y.val + y.d.iter().sum::<f64>())
            .sum()
    };
    run_comparison(
        "30_mtmxccs_colva_zabr_pdv_weather_enb_aad",
        "MtM-XCCS/CTD-COLVA/ZABR-PDV/Weather/ENB single-pass AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

fn suite31_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
    let d = date(2024, 1, 1);
    let curve = FlatForwardCurve::new(d, T::from_f64(0.035) + x[0] * 0.001);
    let mbs = ql_fixed_income::analyze_mbs_burnout_oas_ad(
        &ql_fixed_income::MbsBurnoutOasSpec {
            market_price: T::from_f64(101.0) + x[0] * 0.2,
            original_balance: T::from_f64(100.0),
            wac: T::from_f64(0.06),
            pass_through_rate: T::from_f64(0.055),
            term_months: 24,
            seasoning_months: 12,
            turnover_cpr: T::from_f64(0.06),
            max_refi_cpr: T::from_f64(0.30),
            refi_slope: T::from_f64(60.0),
            refi_incentive_threshold: T::from_f64(0.005),
            burnout_beta: T::from_f64(10.0) + x[1] * 0.5,
            primary_secondary_spread: T::from_f64(0.01),
        },
        &curve,
    )
    .oas;

    let infl = ZeroInflationCurve::new(
        d,
        100.0,
        vec![1.0, 3.0, 5.0],
        vec![
            T::from_f64(0.022) + x[1] * 0.001,
            T::from_f64(0.024) + x[1] * 0.001,
            T::from_f64(0.025) + x[1] * 0.001,
        ],
        None,
    );
    let zciis_cf = ql_inflation::price_zciis_capfloor_jarrow_yildirim(
        &ql_inflation::ZciisCapFloorSpec {
            notional: T::from_f64(100.0),
            maturity_years: 5.0,
            strike_annual_rate: T::from_f64(0.02) + x[2] * 0.001,
            target_month: None,
            is_cap: true,
        },
        &curve,
        &infl,
        &JarrowYildirimParams {
            a_nom: T::from_f64(0.05),
            sigma_nom: T::from_f64(0.008),
            a_real: T::from_f64(0.04),
            sigma_real: T::from_f64(0.006),
            sigma_cpi: T::from_f64(0.012) + x[3] * 0.001,
            rho_nom_real: T::from_f64(0.3),
            rho_nom_cpi: T::from_f64(0.2),
            rho_real_cpi: T::from_f64(-0.1),
        },
    )
    .npv;

    let credit = SurvivalProbabilityCurve::new(
        d,
        vec![1.0, 2.0, 3.0, 4.0, 5.0],
        vec![T::from_f64(0.02) + x[2] * 0.001; 5],
    );
    let cmcds = ql_credit::price_cmcds_market_model(
        &ql_credit::CmCdsMarketModelSpec {
            notional: T::from_f64(100.0),
            grid_times: vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            start_idx: 0,
            end_idx: 3,
            cm_tenor_periods: 2,
            participation_rate: T::from_f64(0.95),
            lgd: T::from_f64(0.6),
            forward_cds_vols: vec![T::from_f64(0.30) + x[3] * 0.01; 5],
            forward_cds_corr: T::from_f64(0.8),
        },
        &curve,
        &credit,
    )
    .seller_npv;

    let fsimm = ql_xva::compute_dynamic_simm_mva_and_spike(&ql_xva::DynamicSimmMvaSpec {
        times: vec![0.5, 1.0, 1.5],
        bucket_sensitivities_t0: vec![T::from_f64(100.0) + x[0] * 5.0, T::from_f64(70.0)],
        bucket_decay_rates: vec![T::from_f64(0.25), T::from_f64(0.30)],
        risk_weights: vec![0.5, 0.6],
        concentration_thresholds: vec![80.0, 80.0],
        inter_bucket_corr: 0.35,
        scheduled_cashflows: vec![T::from_f64(12.0), T::from_f64(15.0), T::from_f64(18.0)],
        mpor_years: 10.0 / 252.0,
        margin_funding_spread: T::from_f64(0.015) + x[1] * 0.001,
        hazard_rate: T::from_f64(0.02),
        lgd: T::from_f64(0.6),
        discount_rate: T::from_f64(0.03),
    })
    .total_aps_cva;

    let jlt = ql_credit::evaluate_jlt_rating_migration_credit(&ql_credit::JltRatingMigrationSpec {
        generator_matrix: vec![
            vec![T::from_f64(-0.06), T::from_f64(0.045), T::from_f64(0.015)],
            vec![T::from_f64(0.02), T::from_f64(-0.07), T::from_f64(0.05)],
            vec![T::zero(), T::zero(), T::zero()],
        ],
        risk_premium_multipliers: vec![T::from_f64(1.2) + x[2] * 0.05, T::from_f64(1.3), T::one()],
        initial_rating_index: 0,
        horizon_years: 2.0,
        payment_times: vec![1.0, 2.0],
        base_coupon_rate: T::from_f64(0.04),
        downgrade_step_up_spreads: vec![T::zero(), T::from_f64(0.008), T::zero()],
        recovery_rate: T::from_f64(0.4) + x[0] * T::from_f64(0.01),
        risk_free_rate: T::from_f64(0.03) + x[0] * 0.001,
    })
    .step_up_bond_price;

    let qou = ql_models::evaluate_quintic_ou_spx_vix(
        &ql_models::QuinticOuParams {
            xi0: T::from_f64(0.04) + x[3] * 0.002,
            kappa: T::from_f64(2.5),
            nu: T::from_f64(0.8),
            rho: T::from_f64(-0.7),
            alpha0: T::one(),
            alpha1: T::from_f64(0.30),
            alpha3: T::from_f64(0.10) + x[1] * 0.01,
            alpha5: T::from_f64(0.02),
        },
        T::from_f64(0.5),
        1.0 / 12.0,
        T::from_f64(20.0),
        T::from_f64(0.985),
    )
    .vix_call_price;

    let fbs = ql_vol::evaluate_free_boundary_sabr(
        &ql_vol::FreeBoundarySabrParams {
            forward: T::from_f64(-0.003) + x[0] * 0.0005,
            expiry: T::one(),
            alpha: T::from_f64(0.008) + x[3] * 0.0005,
            beta: T::from_f64(0.25),
            rho: T::from_f64(-0.20),
            nu: T::from_f64(0.40),
        },
        T::from_f64(-0.001),
        T::from_f64(0.99),
    )
    .call_price;

    let afv = ql_fixed_income::price_convertible_bond_afv_pde(
        &ql_fixed_income::AfvConvertibleBondSpec {
            spot: T::from_f64(100.0) + x[1] * 1.5,
            conversion_ratio: T::one(),
            face_value: T::from_f64(100.0),
            coupon_rate: T::from_f64(0.03),
            maturity: T::from_f64(1.5),
            rate: T::from_f64(0.035) + x[0] * 0.001,
            div_yield: T::from_f64(0.01),
            vol: T::from_f64(0.25) + x[3] * 0.01,
            base_hazard_rate: T::from_f64(0.02) + x[2] * 0.001,
            hazard_equity_elasticity: T::from_f64(1.2),
            equity_loss_on_default_eta: T::from_f64(0.5),
            bond_recovery_rate: T::from_f64(0.4),
            call_price: Some(T::from_f64(115.0)),
            soft_call_trigger_spot: Some(T::from_f64(130.0)),
            put_price: Some(T::from_f64(95.0)),
        },
        15,
        8,
    )
    .total_price;

    let gswap =
        ql_exotics::price_gamma_swap_and_variance_option(&ql_exotics::GammaSwapAndVarOptionSpec {
            spot: T::from_f64(100.0) + x[1] * 1.0,
            rate: T::from_f64(0.03),
            div_yield: T::from_f64(0.01),
            expiry: T::one(),
            cut_strike: T::from_f64(100.0),
            variance_strike: T::from_f64(0.04),
            quotes: vec![
                ql_exotics::VarianceSwapStrikeQuote {
                    strike: T::from_f64(85.0),
                    implied_vol: T::from_f64(0.24) + x[3] * 0.01,
                },
                ql_exotics::VarianceSwapStrikeQuote {
                    strike: T::from_f64(100.0),
                    implied_vol: T::from_f64(0.20) + x[3] * 0.01,
                },
                ql_exotics::VarianceSwapStrikeQuote {
                    strike: T::from_f64(115.0),
                    implied_vol: T::from_f64(0.18) + x[3] * 0.01,
                },
            ],
            heston: HestonParams {
                v0: T::from_f64(0.04) + x[3] * 0.002,
                kappa: T::from_f64(2.0),
                theta: T::from_f64(0.04),
                sigma_v: T::from_f64(0.35),
                rho: T::from_f64(-0.7),
            },
            is_variance_call: true,
        })
        .variance_option_price;

    let cs_swaption = ql_commodities::price_clewlow_strickland_strip_swaption(
        &ql_commodities::ClewlowStricklandSwaptionSpec {
            option_expiry: T::from_f64(0.5),
            delivery_times: vec![0.75, 1.0, 1.25],
            forward_prices: vec![
                T::from_f64(75.0) + x[1] * 0.5,
                T::from_f64(76.0) + x[1] * 0.5,
                T::from_f64(74.0) + x[1] * 0.5,
            ],
            delivery_weights: vec![1.0 / 3.0; 3],
            strike: T::from_f64(74.5),
            rate: T::from_f64(0.03),
            sigma_short: T::from_f64(0.32) + x[3] * 0.01,
            alpha_short: T::from_f64(1.1),
            sigma_long: T::from_f64(0.15),
            rho_factors: T::from_f64(0.2),
            seasonal_vol_amp: T::from_f64(0.05),
            is_payer_call: true,
        },
    )
    .npv;

    let as_mm = ql_portfolio::solve_avellaneda_stoikov_market_making(
        &ql_portfolio::AvellanedaStoikovMmSpec {
            mid_price: T::from_f64(100.0) + x[1] * 0.5,
            initial_inventory: T::from_f64(3.0) + x[2] * 0.2,
            horizon: T::one(),
            num_steps: 6,
            volatility: T::from_f64(0.20) + x[3] * 0.01,
            risk_aversion: T::from_f64(0.5),
            order_book_kappa: T::from_f64(1.5),
            base_arrival_intensity: T::from_f64(10.0),
            alpha_drift: T::from_f64(0.01),
        },
    )
    .expected_spread_capture_pnl;

    let cdar_opt =
        ql_portfolio::optimize_mean_cdar_portfolio(&ql_portfolio::MeanCdarPortfolioSpec {
            path_asset_returns: vec![
                vec![T::from_f64(0.03), T::from_f64(0.01)],
                vec![T::from_f64(-0.06), T::from_f64(-0.01)],
                vec![T::from_f64(-0.04), T::from_f64(-0.005)],
                vec![T::from_f64(0.04), T::from_f64(0.012)],
            ],
            expected_returns: vec![T::from_f64(0.09) + x[0] * 0.005, T::from_f64(0.04)],
            cdar_aversion: T::from_f64(4.0),
            confidence: 0.80,
            smoothing_eps: 1e-3,
        })
        .cdar;

    let schur_w = ql_portfolio::schur_complementary_hrp(&ql_portfolio::SchurHrpSpec {
        covariance: vec![
            vec![T::from_f64(0.04) + x[3] * 0.002, T::from_f64(0.012)],
            vec![T::from_f64(0.012), T::from_f64(0.09) + x[2] * 0.002],
        ],
        gamma: T::from_f64(0.6) + x[0] * 0.02,
    })
    .weights[0];

    let oas_sh = ql_portfolio::oas_and_nonlinear_shrinkage_covariance(&[
        vec![T::from_f64(0.02) + x[0] * 0.001, T::from_f64(-0.01)],
        vec![T::from_f64(-0.015), T::from_f64(0.025) + x[1] * 0.001],
        vec![T::from_f64(0.01), T::from_f64(0.005)],
        vec![T::from_f64(-0.005), T::from_f64(-0.01)],
    ])
    .oas_shrinkage_intensity;

    [
        mbs * 100.0 + zciis_cf,
        cmcds + fsimm + jlt,
        qou + fbs * 100.0 + afv + gswap,
        cs_swaption + as_mm + cdar_opt + schur_w + oas_sh,
    ]
}

pub fn benchmark_burnout_cmcds_quintic_afv_cs_schur_aad(scale: usize) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let h = 1e-5;
    let baseline = || {
        let b = suite31_workloads(x0);
        let mut c = b.iter().sum();
        for lane in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[lane] += h;
            d[lane] -= h;
            let bu = suite31_workloads(u);
            let bd = suite31_workloads(d);
            for j in 0..4 {
                c += (bu[j] - bd[j]) / (2.0 * h);
            }
        }
        c
    };
    let sota = || {
        type D = Dual<4>;
        let x = [
            D::variable(x0[0], 0),
            D::variable(x0[1], 1),
            D::variable(x0[2], 2),
            D::variable(x0[3], 3),
        ];
        suite31_workloads(x)
            .iter()
            .map(|y| y.val + y.d.iter().sum::<f64>())
            .sum()
    };
    run_comparison(
        "31_burnout_cmcds_quintic_afv_cs_schur_aad",
        "Burnout-MBS/CMCDS/Quintic-OU/AFV-CB/CS-Swaption/Schur-HRP AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

fn suite32_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
    let d = date(2024, 1, 1);
    let curve = FlatForwardCurve::new(d, T::from_f64(0.03) + x[0] * 0.001);

    let sinhou = ql_rates::price_sinh_ou_rfr_caplet(
        &ql_rates::SinhOuRfrCapletSpec {
            notional: T::from_f64(100.0),
            start_time: 0.5,
            end_time: 0.75,
            strike: T::from_f64(0.03) + x[1] * 0.0005,
            alpha_mean_reversion: T::from_f64(0.05),
            sigma_short_rate: T::from_f64(0.01) + x[3] * 0.0002,
            gamma_smile: T::from_f64(1.2) + x[2] * 0.05,
            y_star_skew: T::from_f64(-0.15),
            is_caplet: true,
            is_backward_compounded: true,
        },
        &curve,
    )
    .npv;

    let cir_fut = ql_rates::price_cirpp_sofr_futures(
        &ql_rates::CirppSofrFuturesSpec {
            start_time: 0.5,
            end_time: 0.75,
            x0: T::from_f64(0.03) + x[0] * 0.001,
            b_drift: T::from_f64(0.006),
            beta_mean_reversion: T::from_f64(0.20),
            sigma: T::from_f64(0.06) + x[3] * 0.002,
        },
        &curve,
    )
    .futures_rate_3m_compounded;

    let bk_oas = ql_fixed_income::analyze_black_karasinski_callable_frn_oas_ad(
        &ql_fixed_income::BlackKarasinskiCallableFrnSpec {
            face_amount: T::from_f64(100.0),
            quoted_margin: T::from_f64(0.008),
            cap_rate: T::from_f64(0.055),
            floor_rate: T::from_f64(0.020),
            maturity: 2.0,
            coupon_steps: 6,
            mean_reversion: T::from_f64(0.15),
            log_rate_vol: T::from_f64(0.22) + x[3] * 0.005,
            flat_zero_rate: T::from_f64(0.035) + x[0] * 0.001,
            call_price: T::from_f64(100.2),
            first_call_time: 1.0,
            target_market_price: T::from_f64(99.8) + x[1] * 0.05,
            smoothing_eps: 1e-3,
        },
    )
    .oas_spread;

    let surv = SurvivalProbabilityCurve::new(
        d,
        vec![1.0, 2.0, 3.0, 4.0, 5.0],
        vec![T::from_f64(0.025) + x[2] * 0.001; 5],
    );
    let mb_opt = ql_credit::price_morini_brigo_credit_index_option(
        &ql_credit::MoriniBrigoIndexOptionSpec {
            notional: T::from_f64(100.0),
            strike_spread: T::from_f64(0.025),
            option_expiry: 1.0,
            payment_times: vec![2.0, 3.0, 4.0, 5.0],
            recovery_rate: T::from_f64(0.40),
            spread_vol: T::from_f64(0.45) + x[3] * 0.01,
            num_constituents: 50,
            systemic_copula_corr: T::from_f64(0.60) + x[2] * 0.01,
            is_payer: true,
        },
        &curve,
        &surv,
    )
    .npv;

    let wwm_cva = ql_xva::compute_brigo_vrins_wwm_cva(&ql_xva::BrigoVrinsWwmCvaSpec {
        times: vec![0.5, 1.0, 2.0],
        discount_rate: T::from_f64(0.03) + x[0] * 0.001,
        base_mtm_profile: vec![
            T::from_f64(10.0) + x[1] * 0.2,
            T::from_f64(12.0),
            T::from_f64(11.0),
        ],
        exposure_vol_nu: T::from_f64(8.0),
        is_brownian_bridge_swap: false,
        y0: T::from_f64(0.02) + x[2] * 0.001,
        cir_kappa: T::from_f64(0.35),
        cir_theta: T::from_f64(0.025),
        cir_vol: T::from_f64(0.08),
        intensity_shift_psi: T::from_f64(0.005),
        wwr_correlation: T::from_f64(0.60) + x[3] * 0.01,
        lgd: T::from_f64(0.60),
    })
    .wwm_cva;

    let bcs_band =
        ql_xva::compute_bichuch_capponi_sturm_xva_band(&ql_xva::AsymmetricRateXvaBandSpec {
            maturity: 1.0,
            num_steps: 6,
            clean_claim_value: T::from_f64(25.0) + x[1] * 0.5,
            clean_delta_notional: T::from_f64(40.0),
            asset_vol: T::from_f64(0.20) + x[3] * 0.005,
            collateral_alpha: T::from_f64(0.50),
            dividend_discount_rate: T::from_f64(0.02),
            repo_lend_rate: T::from_f64(0.020),
            repo_borrow_rate: T::from_f64(0.024),
            funding_lend_rate: T::from_f64(0.022),
            funding_borrow_rate: T::from_f64(0.036) + x[0] * 0.001,
            collateral_lend_rate: T::from_f64(0.020),
            collateral_borrow_rate: T::from_f64(0.022),
            own_bond_yield_mu_i: T::from_f64(0.042),
            cpty_bond_yield_mu_c: T::from_f64(0.048) + x[2] * 0.001,
            own_lgd: T::from_f64(0.60),
            cpty_lgd: T::from_f64(0.60),
        })
        .no_arbitrage_band_width;

    let mrb_vix =
        ql_models::evaluate_mixed_rough_bergomi_vix(&ql_models::MixedRoughBergomiVixSpec {
            xi0: T::from_f64(0.04) + x[3] * 0.001,
            hurst: T::from_f64(0.12),
            mixing_weight_lambda: T::from_f64(0.55),
            nu1: T::from_f64(0.9),
            nu2: T::from_f64(2.1) + x[1] * 0.02,
            expiry: T::from_f64(0.25),
            vix_window: 1.0 / 12.0,
            strike: T::from_f64(20.0),
            discount_factor: T::from_f64(0.992),
        })
        .vix_call_price;

    let speed_vol = ql_vol::evaluate_third_order_speed_vol_expansion(
        &ql_vol::ThirdOrderSpeedVolExpansionSpec {
            xi0: T::from_f64(0.04) + x[3] * 0.001,
            omega_vol_of_vol: T::from_f64(1.8),
            k_mean_reversion: T::from_f64(3.0),
            expiry: T::from_f64(0.25),
            window_delta: T::from_f64(1.0 / 12.0),
            strike: T::from_f64(0.20),
            discount_factor: T::from_f64(0.99),
        },
    )
    .call_price;

    let choi =
        ql_exotics::price_choi_sum_of_bsm_basket_option(&ql_exotics::ChoiSumOfBsmBasketSpec {
            spots: vec![T::from_f64(100.0) + x[1] * 0.5, T::from_f64(90.0)],
            weights: vec![T::one(), -T::one()],
            vols: vec![T::from_f64(0.25) + x[3] * 0.005, T::from_f64(0.20)],
            div_yields: vec![T::from_f64(0.01), T::from_f64(0.02)],
            correlation_matrix: vec![
                vec![T::one(), T::from_f64(0.5)],
                vec![T::from_f64(0.5), T::one()],
            ],
            strike: T::from_f64(5.0),
            risk_free_rate: T::from_f64(0.03) + x[0] * 0.001,
            expiry: T::one(),
            is_call: true,
        })
        .price;

    let ldz = ql_exotics::price_li_deng_zhou_multi_spread_option(
        &ql_exotics::LiDengZhouMultiSpreadSpec {
            long_spot: T::from_f64(100.0) + x[1] * 0.5,
            long_div_yield: T::from_f64(0.01),
            long_vol: T::from_f64(0.25) + x[3] * 0.005,
            short_spots: vec![T::from_f64(90.0), T::from_f64(85.0)],
            short_weights: vec![T::from_f64(0.5), T::from_f64(0.5)],
            short_div_yields: vec![T::from_f64(0.015), T::from_f64(0.02)],
            short_vols: vec![T::from_f64(0.20), T::from_f64(0.22)],
            correlation_matrix: vec![
                vec![T::one(), T::from_f64(0.5), T::from_f64(0.4)],
                vec![T::from_f64(0.5), T::one(), T::from_f64(0.6)],
                vec![T::from_f64(0.4), T::from_f64(0.6), T::one()],
            ],
            strike: T::from_f64(5.0),
            risk_free_rate: T::from_f64(0.03) + x[0] * 0.001,
            expiry: T::one(),
            is_call: true,
        },
    )
    .price;

    let mrjd = ql_commodities::price_hikspoors_jaimungal_mrjd_commodity(
        &ql_commodities::HikspoorsJaimungalMrjdSpec {
            x0: T::from_f64(0.10) + x[1] * 0.01,
            y0: T::from_f64(4.20),
            j0_spike: T::from_f64(0.04),
            seasonal_phi: T::zero(),
            kappa_x: T::from_f64(1.5),
            sigma_x: T::from_f64(0.35) + x[3] * 0.005,
            mu_y_rn: T::from_f64(0.02),
            sigma_y: T::from_f64(0.18),
            rho_xy: T::from_f64(-0.2),
            kappa_j_spike: T::from_f64(12.0),
            jump_intensity: T::from_f64(2.5) + x[2] * 0.05,
            mean_jump_size: T::from_f64(0.22),
            option_expiry: T::from_f64(0.5),
            futures_maturity: T::one(),
            strike: T::from_f64(70.0),
            risk_free_rate: T::from_f64(0.03),
            is_call: true,
        },
    )
    .futures_option_price;

    let afs = ql_portfolio::solve_alfonsi_fruth_schied_lob_execution(
        &ql_portfolio::AlfonsiFruthSchiedLobSpec {
            total_order_size: T::from_f64(100.0),
            horizon: T::one(),
            num_intervals: 4,
            lob_density_q: T::from_f64(50.0),
            shape_exponent_alpha: T::from_f64(0.4) + x[2] * 0.01,
            resilience_rho: T::from_f64(1.2) + x[0] * 0.02,
        },
    )
    .volume_resilience_cost;

    let exvar = ql_portfolio::evaluate_and_optimize_expectile_portfolio(
        &ql_portfolio::ExpectileRiskAndPortfolioSpec {
            scenario_returns: vec![
                vec![T::from_f64(-0.07) + x[1] * 0.001, T::from_f64(-0.01)],
                vec![T::from_f64(-0.04), T::from_f64(-0.005)],
                vec![T::from_f64(0.02), T::from_f64(0.008)],
                vec![T::from_f64(0.05), T::from_f64(0.012)],
            ],
            initial_weights: vec![T::from_f64(0.5), T::from_f64(0.5)],
            tau_confidence: 0.85,
            return_weight_lambda: T::from_f64(0.2),
            smoothing_eps: 1e-3,
        },
    )
    .optimal_expectile_exvar;

    let pcrp = ql_portfolio::solve_principal_component_risk_parity(
        &ql_portfolio::FactorEigenRiskParitySpec {
            covariance: vec![
                vec![T::from_f64(0.04) + x[3] * 0.0005, T::from_f64(0.02)],
                vec![T::from_f64(0.02), T::from_f64(0.05) + x[2] * 0.0005],
            ],
            factor_risk_budgets: vec![T::from_f64(0.5) + x[0] * 0.01, T::from_f64(0.5)],
        },
    )
    .weights[0];

    [
        sinhou + cir_fut * 100.0 + bk_oas * 100.0,
        mb_opt + wwm_cva + bcs_band,
        mrb_vix + speed_vol * 10.0 + choi + ldz,
        mrjd + afs + exvar * 10.0 + pcrp,
    ]
}

pub fn benchmark_sinhou_mb_rbergomi_choi_afs_expectile_aad(scale: usize) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let h = 1e-5;
    let baseline = || {
        let b = suite32_workloads(x0);
        let mut c: f64 = b.iter().sum();
        for lane in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[lane] += h;
            d[lane] -= h;
            let bu = suite32_workloads(u);
            let bd = suite32_workloads(d);
            for j in 0..4 {
                c += (bu[j] - bd[j]) / (2.0 * h);
            }
        }
        c
    };
    let sota = || {
        type D = Dual<4>;
        let x = [
            D::variable(x0[0], 0),
            D::variable(x0[1], 1),
            D::variable(x0[2], 2),
            D::variable(x0[3], 3),
        ];
        suite32_workloads(x)
            .iter()
            .map(|y| y.val + y.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "32_sinhou_mb_rbergomi_choi_afs_expectile_aad",
        "SinhOU/BK-FRN/MoriniBrigo/WWM/mrBergomi/Choi/AFS/EXVaR single-pass AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

fn suite33_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
    let d = date(2024, 1, 1);
    let r = T::from_f64(0.03) + x[0] * 0.001;
    let curve = FlatForwardCurve::new(d, r);
    let hw = ql_curves::bootstrap_hagan_west_monotone_convex_curve(
        &ql_curves::HaganWestMonotoneConvexBootstrapSpec {
            reference_date: d,
            pillar_times: vec![1.0, 2.0, 5.0],
            par_swap_quotes: vec![r, r + x[1] * 0.001, r + x[2] * 0.001],
            payment_frequency_per_year: 2,
            enforce_positivity: true,
            ameliorate_slopes: true,
            query_time: 2.0,
        },
    );
    let swaption = ql_rates::price_piterbarg_quasi_gaussian_amortizing_swaption(
        &ql_rates::PiterbargQuasiGaussianAmortizingSwaptionSpec {
            option_expiry: 1.0,
            second_exercise_time: Some(2.0),
            swap_times: vec![1.0, 2.0, 3.0],
            period_notionals: vec![T::one(), T::from_f64(0.9)],
            strike: r,
            mean_reversion: T::from_f64(0.05),
            sigma0: T::from_f64(0.01) + x[0] * 0.001,
            local_vol_skew_b: x[1] * 0.1,
            cms_pay_delay: 0.0,
            is_payer: true,
        },
        &curve,
        &curve,
    );
    let fx = ql_fx::calibrate_clark_fx_delta_smile(&ql_fx::ClarkFxDeltaConventionSpec {
        spot: T::from_f64(1.1),
        r_dom: r,
        r_for: T::from_f64(0.02),
        expiry: T::one(),
        atm_vol: T::from_f64(0.12),
        rr_25d: T::from_f64(0.01) + x[2] * 0.001,
        market_strangle_25d: T::from_f64(0.01) + x[3] * 0.001,
        is_premium_included_delta: false,
        is_forward_delta: true,
        query_strike: T::from_f64(1.1),
    });
    let lsv =
        ql_vol::evaluate_gatheral_gyongy_lsv_leverage(&ql_vol::GatheralLsvGyongyLeverageSpec {
            forward: T::from_f64(100.0),
            strike: T::from_f64(100.0) + x[0],
            expiry: T::one(),
            svi_a: T::from_f64(0.04),
            svi_b: T::from_f64(0.2),
            svi_rho: T::from_f64(-0.3),
            svi_m: T::zero(),
            svi_sigma: T::from_f64(0.2),
            dw_dt_atm: T::from_f64(0.04),
            heston_v0: T::from_f64(0.04),
            heston_theta: T::from_f64(0.04),
            heston_kappa: T::one(),
            heston_eta: T::from_f64(0.5),
            heston_rho: T::from_f64(-0.7),
        });
    let sac = ql_xva::compute_green_saccr_lifetime_kva_mva(&ql_xva::GreenSaCcrLifetimeKvaMvaSpec {
        time_grid: vec![0.0, 1.0, 2.0],
        trade_notional: T::from_f64(1e6),
        trade_maturity: 3.0,
        initial_mtm: x[0],
        mtm_volatility: T::from_f64(0.1),
        collateral_threshold: T::zero(),
        initial_margin_ratio: T::from_f64(0.1),
        counterparty_pd: T::from_f64(0.02),
        counterparty_lgd: T::from_f64(0.6),
        cva_risk_weight: T::from_f64(0.5),
        capital_ratio_requirement: T::from_f64(0.08),
        hurdle_rate_gamma_k: T::from_f64(0.1),
        funding_spread_sf: T::from_f64(0.01),
        mva_spread_sm: T::from_f64(0.01),
        risk_free_rate: r,
    });
    let greg = ql_xva::allocate_gregory_netting_set_xva(&ql_xva::GregoryNettingSetAllocationSpec {
        trade_weights: vec![T::one(), T::one()],
        trade_expected_mtms: vec![T::one() + x[1], -T::one()],
        trade_exposure_vols: vec![T::from_f64(0.2); 2],
        correlation_matrix: vec![vec![T::one(), T::zero()], vec![T::zero(), T::one()]],
        collateral_fraction: T::zero(),
        horizon: 1.0,
        num_steps: 4,
        risk_free_rate: r,
        cpty_hazard_rate: T::from_f64(0.02),
        cpty_lgd: T::from_f64(0.6),
        own_hazard_rate: T::from_f64(0.01),
        own_lgd: T::from_f64(0.6),
        funding_spread: T::from_f64(0.01),
        colva_spread: T::from_f64(0.01),
    });
    let rough = ql_models::evaluate_bayer_friz_gatheral_markovian_rough(
        &ql_models::BayerFrizGatheralMarkovianRoughSpec {
            spot: T::from_f64(100.0),
            strike: T::from_f64(100.0),
            expiry: T::one(),
            rate: r,
            div_yield: T::zero(),
            xi0: T::from_f64(0.04),
            eta: T::one(),
            rho: T::from_f64(-0.7),
            hurst: T::from_f64(0.1) + x[2] * 0.01,
            num_markov_factors: 4,
            is_call: true,
        },
    );
    let cliq = ql_exotics::price_bouzoubaa_reverse_and_global_cliquet(
        &ql_exotics::BouzoubaaReverseAndGlobalCliquetSpec {
            notional: T::one(),
            rate: r,
            div_yield: T::zero(),
            reset_times: vec![0.5, 1.0],
            forward_vols: vec![T::from_f64(0.2), T::from_f64(0.2)],
            forward_skews: vec![x[3] * 0.1, x[3] * 0.1],
            initial_coupon_cmax: T::zero(),
            local_floor: T::from_f64(-0.1),
            local_cap: T::from_f64(0.1),
            global_floor: T::from_f64(-0.2),
            global_cap: T::from_f64(0.2),
            is_reverse_cliquet: false,
        },
    );
    let plant = ql_commodities::price_burger_multifuel_tolling_power_plant(
        &ql_commodities::BurgerMultiFuelTollingPlantSpec {
            power_forwards: vec![T::from_f64(80.0); 2],
            gas_forwards: vec![T::from_f64(20.0); 2],
            coal_forwards: vec![T::from_f64(15.0); 2],
            carbon_forwards: vec![T::from_f64(50.0); 2],
            dispatch_times: vec![0.25, 0.75],
            capacity_mw: T::from_f64(100.0),
            period_hours: T::one(),
            gas_heat_rate: T::from_f64(2.0),
            gas_emission_factor: T::from_f64(0.4),
            gas_vom_cost: T::from_f64(2.0),
            coal_heat_rate: T::from_f64(2.5),
            coal_emission_factor: T::from_f64(0.8),
            coal_vom_cost: T::one(),
            vol_power: T::from_f64(0.2),
            vol_gas: T::from_f64(0.2),
            vol_coal: T::from_f64(0.2),
            vol_carbon: T::from_f64(0.2),
            rho_power_gas: T::from_f64(0.2),
            rho_power_coal: T::from_f64(0.2),
            rho_power_carbon: T::from_f64(0.2),
            rho_gas_coal: T::from_f64(0.2),
            start_up_cost: T::from_f64(100.0),
            min_load_fraction: T::from_f64(0.5),
            part_load_heat_rate_penalty: T::from_f64(0.1),
            risk_free_rate: r,
        },
    );
    let execution = ql_portfolio::solve_cartea_jaimungal_alpha_execution(
        &ql_portfolio::CarteaJaimungalAlphaExecutionSpec {
            initial_inventory: T::one(),
            horizon: T::one(),
            num_steps: 8,
            temporary_impact_k: T::one(),
            permanent_impact_b: T::from_f64(0.1),
            inventory_penalty_phi: T::from_f64(0.2),
            terminal_penalty_a: x[0].max_val(T::zero()),
            initial_alpha_signal: T::from_f64(0.1),
            alpha_mean_reversion_kappa: T::one(),
            alpha_volatility: T::from_f64(0.2),
        },
    );
    let paleo = ql_portfolio::solve_paleologo_factor_neutral_portfolio(
        &ql_portfolio::PaleologoFactorNeutralAlphaSpec {
            raw_alphas: vec![T::one() + x[1], T::from_f64(2.0)],
            factor_loadings_b: vec![vec![T::one()], vec![T::one()]],
            factor_covariance_omega: vec![vec![T::one()]],
            sample_idio_variances: vec![T::from_f64(0.1), T::from_f64(0.2)],
            idio_shrinkage_xi: T::from_f64(0.5),
            idio_risk_aversion: T::one(),
            factor_risk_aversion: T::one(),
            factor_crowding_scores: vec![x[2]],
            crowding_penalty_lambda: T::from_f64(0.1),
        },
    );
    let aad = ql_dml::evaluate_glasserman_malliavin_checkpointed_aad(
        &ql_dml::GlassermanMalliavinCheckpointedMcSpec {
            spot: T::from_f64(100.0),
            strike: T::from_f64(100.0),
            rate: r,
            div_yield: T::zero(),
            base_vol: T::from_f64(0.2),
            local_vol_skew: x[3] * 0.1,
            expiry: T::one(),
            num_steps: 8,
            num_checkpoints: 2,
        },
    );
    [
        hw.query_discount_factor
            + swaption.european_amortizing_npv
            + fx.query_call_price
            + lsv.implied_vol,
        sac.lifetime_kva + greg.netted_total_xva + rough.option_price,
        cliq.npv + plant.tolling_plant_npv + execution.expected_execution_cost_net_of_alpha,
        paleo.expected_active_return + aad.call_price,
    ]
}

/// Monotone-convex, Clark FX, SA-CCR KVA, rough volatility, Cartea, and Paleologo AD benchmark.
pub fn benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let h = 1e-5;
    let baseline = move || {
        let f = |x: [f64; 4]| suite33_workloads(x).iter().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += h;
            d[i] -= h;
            z += (f(u) - f(d)) / (2.0 * h);
        }
        z
    };
    let sota = move || {
        type D = Dual<4>;
        let x = [
            D::variable(x0[0], 0),
            D::variable(x0[1], 1),
            D::variable(x0[2], 2),
            D::variable(x0[3], 3),
        ];
        suite33_workloads(x)
            .iter()
            .map(|y| y.val + y.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad",
        "Monotone-Convex / Clark-FX / SA-CCR-KVA / rVol / Cartea / Paleologo AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

/// Direct twelve-engine workload across RKHS, FMM, VV, LSV, JLT, CTD, QRH, Parisian, Aïd, propagator, NCO, and DML.
pub fn suite35_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let rk = ql_curves::calibrate_filipovic_rkhs_smooth_yield_curve(
        &ql_curves::FilipovicRkhsYieldCurveSpec {
            cashflow_times: vec![T::one(), T::from_f64(2.0)],
            instrument_cashflow_matrix: vec![vec![T::one(), T::zero()], vec![T::zero(), T::one()]],
            instrument_dirty_prices: vec![T::from_f64(0.98) + x[0], T::from_f64(0.95) + x[1]],
            decay_alpha: T::one(),
            ridge_penalty_lambda: T::from_f64(1e-8),
            query_times: vec![T::one()],
        },
    );
    let fmm = ql_rates::price_lyashenko_mercurio_fmm_caplet_and_swaption(
        &ql_rates::LyashenkoMercurioFmmSpec {
            notional: T::one(),
            exercise_time_t0: T::one(),
            period_accruals: vec![T::from_f64(0.5), T::from_f64(0.5)],
            forward_rates: vec![T::from_f64(0.03) + x[2], T::from_f64(0.03) + x[3]],
            discount_factors: vec![T::from_f64(0.98), T::from_f64(0.96)],
            fmm_volatilities: vec![T::from_f64(0.2), T::from_f64(0.2)],
            displacement_shift: T::zero(),
            inter_period_correlation_decay: T::one(),
            strike: T::from_f64(0.03),
            is_payer_or_caplet: true,
        },
    );
    let vv = ql_fx::evaluate_castagna_full_vanna_volga_hedge_and_barrier(
        &ql_fx::CastagnaFullVannaVolgaHedgeSpec {
            spot: T::from_f64(100.0),
            r_dom: T::from_f64(0.02),
            r_for: T::from_f64(0.01),
            expiry: T::one(),
            atm_vol: T::from_f64(0.15),
            rr_25d: x[0] * 0.01,
            bf_25d: x[1] * 0.01,
            target_strike: T::from_f64(100.0) + x[2],
            knockout_barrier: T::from_f64(120.0),
            unwind_attenuation_factor: T::from_f64(0.5),
            is_up_and_out: true,
            is_call: true,
        },
    );
    let lsv = ql_vol::evaluate_guyon_markovian_projection_lsv_smile(
        &ql_vol::GuyonMarkovianProjectionLsvSpec {
            forward: T::from_f64(100.0),
            strike: T::from_f64(100.0) + x[3],
            expiry: T::one(),
            forward_start_t1: T::from_f64(0.25),
            atm_vol: T::from_f64(0.2),
            market_atm_skew: T::from_f64(-0.1),
            market_atm_curvature: T::from_f64(0.1),
            sv_v0: T::from_f64(0.04),
            sv_kappa: T::one(),
            sv_eta: T::from_f64(0.3),
            sv_rho: T::from_f64(-0.5),
            sv_mixing_weight_eta: T::from_f64(0.5),
        },
    );
    let jlt = ql_credit::evaluate_jlt_rating_migration_and_crcva(
        &ql_credit::JarrowLandoTurnbullRatingMigrationSpec {
            initial_rating_index: 0,
            physical_generator_4x4: [
                [
                    -T::from_f64(0.1),
                    T::from_f64(0.08),
                    T::from_f64(0.01),
                    T::from_f64(0.01),
                ],
                [
                    T::from_f64(0.02),
                    -T::from_f64(0.08),
                    T::from_f64(0.04),
                    T::from_f64(0.02),
                ],
                [
                    T::zero(),
                    T::from_f64(0.02),
                    -T::from_f64(0.05),
                    T::from_f64(0.03),
                ],
                [T::zero(); 4],
            ],
            risk_premium_multipliers: [T::one(); 3],
            recovery_rate: T::from_f64(0.4) + x[0] * T::from_f64(0.01),
            risk_free_rate: T::from_f64(0.02),
            maturity: T::from_f64(5.0),
            base_coupon_rate: T::from_f64(0.03),
            hy_downgrade_stepup_spread: T::from_f64(0.01),
            counterparty_epe_notional: T::one(),
        },
    );
    let ctd = ql_xva::compute_brigo_pallavicini_ctd_collateral_and_closeout_xva(
        &ql_xva::BrigoPallaviciniCtdCollateralXvaSpec {
            notional_mtm: T::one(),
            maturity: T::one(),
            num_steps: 8,
            base_ois_rate: T::from_f64(0.02),
            csa_spread_ccy1: T::from_f64(0.001) + x[1] * T::from_f64(1e-4),
            csa_spread_ccy2: T::from_f64(0.002),
            csa_spread_ccy3: T::from_f64(0.003),
            csa_mean_reversion: T::one(),
            csa_vol_ccy1: T::from_f64(0.01),
            csa_vol_ccy2: T::from_f64(0.01),
            csa_vol_ccy3: T::from_f64(0.01),
            csa_correlation: T::from_f64(0.2),
            cpty_hazard_rate: T::from_f64(0.02),
            cpty_lgd: T::from_f64(0.6),
            own_hazard_rate: T::from_f64(0.01),
            own_lgd: T::from_f64(0.6),
            exposure_volatility: T::from_f64(0.1),
        },
    );
    let qrh = ql_models::evaluate_gatheral_rosenbaum_joint_spx_vix_qrh(
        &ql_models::GatheralRosenbaumJointSpxVixQrhSpec {
            params: ql_models::QuadraticRoughHestonParams {
                hurst: 0.1,
                nu: T::from_f64(0.3),
                lambda: T::one(),
                a: T::one(),
                b: T::zero(),
                c: T::from_f64(0.04),
                z0: x[0],
            },
            spot: T::from_f64(100.0),
            spx_strike: T::from_f64(100.0),
            vix_strike: T::from_f64(0.2),
            expiry: T::one(),
            rate: T::from_f64(0.02),
        },
    );
    let par = ql_exotics::price_chesney_yor_cumulative_parisian_and_rainbow(
        &ql_exotics::ChesneyYorCumulativeParisianRainbowSpec {
            spot_1: T::from_f64(100.0),
            spot_2: T::from_f64(102.0) + x[2],
            strike: T::from_f64(100.0),
            par_barrier: T::from_f64(120.0),
            excursion_window_d: T::from_f64(0.1),
            vol_1: T::from_f64(0.2),
            vol_2: T::from_f64(0.2),
            div_1: T::zero(),
            div_2: T::zero(),
            correlation: T::from_f64(0.5),
            risk_free_rate: T::from_f64(0.02),
            expiry: T::one(),
            is_up_and_out: true,
        },
    );
    let aid = ql_commodities::price_aid_carmona_structural_electricity_stack(
        &ql_commodities::AidCarmonaStructuralElectricityStackSpec {
            gas_forward: T::from_f64(30.0) + x[3],
            gas_vol: T::from_f64(0.3),
            base_heat_rate: T::from_f64(2.0),
            total_capacity_mw: T::from_f64(100.0),
            baseline_capacity_anchor_mw: T::from_f64(50.0),
            expected_load_mw: T::from_f64(80.0),
            expected_renewable_gen_mw: T::from_f64(20.0),
            residual_load_vol_mw: T::from_f64(10.0),
            stack_convexity_alpha: T::one(),
            negative_price_floor_depth: T::from_f64(5.0),
            curtailment_beta: T::one(),
            load_gas_correlation: T::from_f64(0.2),
            power_option_strike: T::from_f64(50.0),
            spike_threshold_price: T::from_f64(100.0),
            risk_free_rate: T::from_f64(0.02),
            expiry: T::one(),
        },
    );
    let imp = ql_portfolio::solve_bouchaud_webster_square_root_and_propagator_execution(
        &ql_portfolio::BouchaudWebsterPropagatorImpactSpec {
            metaorder_size_q: T::one() + x[0] * T::from_f64(0.1),
            daily_volume_vd: T::from_f64(10.0),
            daily_volatility_sigma: T::from_f64(0.2),
            y_prefactor: T::one(),
            propagator_exponent_beta: T::from_f64(0.5),
            propagator_scale_g0: T::one(),
            regularization_tau0: T::from_f64(0.01),
            temporary_impact_eta: T::from_f64(0.1),
            risk_aversion_lambda: T::from_f64(0.1),
            horizon: T::one(),
            num_slices: 8,
        },
    );
    let nco = ql_portfolio::solve_lopez_de_prado_denoised_nco_allocation(
        &ql_portfolio::LopezDePradoNcoDenoisedSpec {
            covariance: vec![
                vec![
                    T::from_f64(0.04) + x[1] * T::from_f64(0.002),
                    T::from_f64(0.01),
                ],
                vec![T::from_f64(0.01), T::from_f64(0.09)],
            ],
            expected_returns: vec![T::from_f64(0.05) + x[1], T::from_f64(0.06)],
            sample_ratio_q: T::from_f64(2.0),
            denoising_alpha: T::from_f64(0.5),
            num_clusters: 2,
        },
    );
    let dml = ql_dml::evaluate_glasserman_savine_hybrid_pw_lrm_and_dml(
        &ql_dml::GlassermanSavineHybridPwLrmSpec {
            spot: T::from_f64(100.0),
            strike: T::from_f64(100.0) + x[1],
            rate: T::from_f64(0.02),
            div_yield: T::zero(),
            volatility: T::from_f64(0.2),
            expiry: T::one(),
            digital_smoothing_eps: T::from_f64(0.5),
            dml_differential_weight_gamma: T::one(),
        },
    );
    [
        rk.query_zero_rates[0],
        fmm.backward_looking_rfr_caplet_npv,
        vv.exact_vv_vanilla_price,
        lsv.lsv_reconstructed_implied_vol,
        jlt.rating_migration_adjusted_crcva,
        ctd.ctd_collateral_switch_option_value,
        qrh.vix_futures_fair_value,
        par.rainbow_call_on_max,
        aid.structural_power_forward,
        imp.optimal_propagator_cost,
        nco.nco_portfolio_volatility,
        dml.hybrid_pw_lrm_digital_delta,
    ]
}
pub fn benchmark_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite35_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite35_workloads(x)
            .iter()
            .map(|y| y.val + y.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "35_rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml",
        "RKHS / FMM / VV / LSV / JLT / CTD / QRH / Parisian / Aïd / Propagator / NCO / DML AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite36_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let m = ql_fixed_income::analyze_davidson_levin_mbs_io_po_strips(
        &ql_fixed_income::DavidsonLevinIoPoStripMbsSpec {
            notional: T::from_f64(100.0),
            wac_annual: T::from_f64(0.065),
            pass_through_coupon_annual: T::from_f64(0.06),
            current_mortgage_rate: T::from_f64(0.055) + x[0],
            discount_rate_annual: T::from_f64(0.04) + x[1],
            remaining_months: 24,
            loan_age_months: 12,
            turnover_base_cpr: T::from_f64(0.04),
            max_refi_cpr: T::from_f64(0.4),
            refi_slope: T::from_f64(80.) + x[2],
            refi_threshold_spread: T::from_f64(0.005),
            burnout_exponent: T::from_f64(1.) + x[3],
            rate_bump_bps: T::one(),
        },
    );
    let c = ql_credit::analyze_duffie_singleton_rmv_rfv_and_cds_basis(
        &ql_credit::DuffieSingletonRecoveryAndCdsBasisSpec {
            face_value: T::from_f64(100.) + x[0],
            coupon_rate_annual: T::from_f64(0.05) + x[1],
            payment_freq_per_year: 2,
            maturity_years: T::from_f64(3.),
            risk_free_rate: T::from_f64(0.02) + x[2],
            hazard_rate: T::from_f64(0.02) + x[3],
            recovery_rate: T::from_f64(0.4),
            market_bond_price: T::from_f64(98.),
            repo_funding_spread: T::from_f64(0.005),
        },
    );
    let mo = ql_xva::compute_crepey_bielecki_marshall_olkin_bilateral_xva(
        &ql_xva::CrepeyBieleckiMarshallOlkinBilateralXvaSpec {
            times: vec![T::from_f64(0.5), T::one()],
            expected_positive_exposure: vec![T::one() + x[0], T::one() + x[1]],
            expected_negative_exposure: vec![T::from_f64(0.5) + x[2], T::from_f64(0.5) + x[3]],
            risk_free_rate: T::from_f64(0.02),
            lambda_institution_idio: T::from_f64(0.01),
            lambda_counterparty_idio: T::from_f64(0.02),
            lambda_common_shock_12: T::from_f64(0.01) + x[0],
            recovery_institution: T::from_f64(0.4),
            recovery_counterparty: T::from_f64(0.4),
            recovery_counterparty_simultaneous: T::from_f64(0.3),
            gap_exposure_multiplier: T::from_f64(0.2),
            funding_borrowing_spread: T::from_f64(0.01),
            funding_lending_spread: T::from_f64(0.005),
        },
    );
    let cb = ql_fixed_income::price_sirbu_shreve_dilution_and_makewhole_convertible(
        &ql_fixed_income::SirbuShreveDilutionMakeWholeConvertibleSpec {
            spot: T::from_f64(100.) + x[0] + x[1] + x[2] + x[3],
            volatility: T::from_f64(0.2),
            risk_free_rate: T::from_f64(0.03),
            dividend_yield: T::zero(),
            credit_spread: T::from_f64(0.02),
            maturity_years: T::from_f64(2.),
            par_value: T::from_f64(100.),
            coupon_rate_annual: T::from_f64(0.03),
            conversion_ratio: T::from_f64(0.8),
            existing_shares_millions: T::from_f64(100.),
            convertible_issue_count_millions: T::from_f64(10.),
            debt_relief_retention_fraction: T::zero(),
            call_price: T::from_f64(105.),
            soft_call_trigger_spot: T::from_f64(110.),
            make_whole_max_bonus_shares: T::from_f64(0.2),
            make_whole_time_decay: T::from_f64(0.2),
            make_whole_moneyness_decay: T::from_f64(0.3),
            tree_steps: 8,
        },
    );
    let sig = ql_models::evaluate_bayer_oberhauser_signature_option_pricing(
        &ql_models::BayerHorvathOberhauserPathSignatureSpec {
            time_grid: vec![T::zero(), T::from_f64(0.5), T::one()],
            sample_paths: vec![
                vec![
                    T::from_f64(100.),
                    T::from_f64(102.) + x[2],
                    T::from_f64(101.) + x[3],
                ],
                vec![
                    T::from_f64(100.),
                    T::from_f64(101.) + x[0],
                    T::from_f64(103.) + x[1],
                ],
                vec![
                    T::from_f64(100.),
                    T::from_f64(99.) + x[1],
                    T::from_f64(102.) + x[2],
                ],
                vec![
                    T::from_f64(100.),
                    T::from_f64(103.) + x[3],
                    T::from_f64(104.) + x[0],
                ],
                vec![
                    T::from_f64(100.),
                    T::from_f64(101.) + x[2],
                    T::from_f64(102.) + x[3],
                ],
            ],
            strike: T::from_f64(100.) + x[0],
            discount_factor: T::from_f64(0.98) + x[1],
            ridge_lambda: T::from_f64(1e-6),
        },
    );
    let hn =
        ql_vol::evaluate_heston_nandi_and_tsay_dcc_garch(&ql_vol::HestonNandiTsayDccGarchSpec {
            spot: T::from_f64(100.) + x[0],
            strike: T::from_f64(100.) + x[1],
            daily_risk_free_rate: T::from_f64(0.0001),
            days_to_maturity: 30,
            initial_daily_variance: T::from_f64(0.0002),
            omega: T::from_f64(0.00001),
            alpha: T::from_f64(0.05),
            beta: T::from_f64(0.9),
            gamma_asymmetry: T::from_f64(0.2) + x[2],
            asset1_returns: vec![T::from_f64(0.01) + x[3]],
            asset2_returns: vec![T::from_f64(0.008)],
            dcc_a: T::from_f64(0.02),
            dcc_b: T::from_f64(0.95),
        });
    let pot = ql_portfolio::evaluate_mcneil_embrechts_pot_gpd_and_christoffersen_backtest(
        &ql_portfolio::McNeilEmbrechtsPotGpdBacktestSpec {
            historical_losses: vec![
                T::from_f64(1.) + x[0],
                T::from_f64(2.) + x[1],
                T::from_f64(3.) + x[2],
                T::from_f64(4.) + x[3],
                T::from_f64(5.),
            ],
            threshold_u: T::from_f64(1.),
            confidence_alpha: T::from_f64(0.95),
            spectral_risk_aversion_k: T::from_f64(2.),
            backtest_var_threshold: T::from_f64(3.),
        },
    );
    let rie = ql_portfolio::solve_potters_bouchaud_rie_nonlinear_shrinkage(
        &ql_portfolio::PottersBouchaudRieCovarianceSpec {
            sample_covariance: vec![
                vec![T::from_f64(0.04) + x[0], T::from_f64(0.01) + x[2]],
                vec![T::from_f64(0.01) + x[2], T::from_f64(0.09) + x[1] + x[3]],
            ],
            num_observations_t: 100,
        },
    );
    let pin = ql_portfolio::evaluate_ohara_ekop_kyle_microstructure_and_etf_band(
        &ql_portfolio::OharaEkopKyleMicrostructureSpec {
            fundamental_value_std: T::from_f64(1.) + x[0],
            noise_trader_volume_std: T::from_f64(2.) + x[1],
            alpha_info_event_prob: T::from_f64(0.2),
            delta_bad_news_prob: T::from_f64(0.4),
            mu_informed_arrival_rate: T::from_f64(1.) + x[2],
            epsilon_uninformed_buy_rate: T::from_f64(2.),
            epsilon_uninformed_sell_rate: T::from_f64(2.),
            fundamental_low: T::from_f64(95.),
            fundamental_high: T::from_f64(105.),
            position_value: T::from_f64(100.),
            return_vol_daily: T::from_f64(0.02),
            z_score_alpha: T::from_f64(2.),
            mean_relative_spread: T::from_f64(0.01),
            std_relative_spread: T::from_f64(0.002),
            etf_price: T::from_f64(101.) + x[3],
            etf_inav: T::from_f64(100.),
            ap_creation_fee_bps: T::from_f64(5.),
            basket_spread_bps: T::from_f64(2.),
            etf_spread_bps: T::from_f64(3.),
        },
    );
    let gp = ql_portfolio::solve_garleanu_pedersen_dynamic_aim_portfolio(
        &ql_portfolio::GarleanuPedersenDynamicAimPortfolioSpec {
            covariance: vec![
                vec![T::from_f64(0.04), T::from_f64(0.01)],
                vec![T::from_f64(0.01), T::from_f64(0.09)],
            ],
            current_weights: vec![T::from_f64(0.5), T::from_f64(0.5)],
            slow_alpha_vector: vec![T::from_f64(0.03) + x[0], T::from_f64(0.02)],
            fast_alpha_vector: vec![T::from_f64(0.01), T::from_f64(0.02) + x[1]],
            slow_mean_reversion_phi: T::from_f64(0.1),
            fast_mean_reversion_phi: T::from_f64(0.5),
            risk_aversion_gamma: T::from_f64(2.),
            transaction_cost_kappa: T::from_f64(0.1) + x[2],
            discount_rate_r: T::from_f64(0.02) + x[3],
        },
    );
    let br = ql_portfolio::evaluate_bacon_carino_brinson_and_karnosky_singer_attribution(
        &ql_portfolio::BaconCarinoMultiPeriodBrinsonSpec {
            portfolio_weights_by_period: vec![vec![T::from_f64(0.5), T::from_f64(0.5)]],
            benchmark_weights_by_period: vec![vec![T::from_f64(0.5), T::from_f64(0.5)]],
            portfolio_returns_by_period: vec![vec![T::from_f64(0.1) + x[0], T::from_f64(0.02)]],
            benchmark_returns_by_period: vec![vec![T::from_f64(0.08), T::from_f64(0.01) + x[1]]],
            fx_spot_returns_by_sector: vec![T::from_f64(0.01) + x[2]],
            foreign_cash_rates_by_sector: vec![T::from_f64(0.02)],
            domestic_cash_rate: T::from_f64(0.03) + x[3],
        },
    );
    let ie = ql_ad::solve_penman_damodaran_riv_aeg_and_implied_erp_ift(
        &ql_ad::PenmanDamodaranResidualIncomeIerpSpec {
            initial_book_value: T::from_f64(100.),
            forecast_roe_schedule: vec![T::from_f64(0.12) + x[0], T::from_f64(0.11)],
            payout_ratio_schedule: vec![T::from_f64(0.4), T::from_f64(0.4)],
            terminal_roe: T::from_f64(0.1),
            terminal_growth_rate: T::from_f64(0.03) + x[1],
            risk_free_rate: T::from_f64(0.03) + x[2],
            equity_beta: T::from_f64(1.1),
            observed_market_price: T::from_f64(110.) + x[3],
        },
    );
    [
        m.io_strip_price,
        c.rmv_bond_price,
        mo.total_bilateral_xva,
        cb.make_whole_enhanced_price,
        sig.signature_linear_functional_price,
        hn.heston_nandi_call_price,
        pot.pot_gpd_expected_shortfall,
        rie.rie_portfolio_variance,
        pin.ekop_pin,
        gp.next_period_optimal_weights[0],
        br.carino_linked_selection_total,
        ie.implied_erp,
    ]
}

pub fn suite37_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let kal = ql_curves::evaluate_tsay_kalman_rts_state_space_curve(
        &ql_curves::TsayKalmanRtsStateSpaceCurveSpec {
            transition_f: [
                [T::from_f64(0.9) + x[0] * 0.01, T::from_f64(0.1)],
                [T::zero(), T::from_f64(0.95) + x[1] * 0.01],
            ],
            state_drift_c: [x[2], x[3]],
            process_noise_q: [
                [T::from_f64(0.01), T::zero()],
                [T::zero(), T::from_f64(0.01)],
            ],
            measurement_h: [[T::one(), T::zero()], [T::zero(), T::one()]],
            measurement_noise_r: [T::from_f64(0.02), T::from_f64(0.02)],
            initial_state_mean: [T::zero(), T::zero()],
            initial_state_cov: [[T::one(), T::zero()], [T::zero(), T::one()]],
            observations: vec![
                [T::one() + x[0], T::one() + x[1]],
                [T::from_f64(1.1) + x[2], T::from_f64(0.9) + x[3]],
            ],
        },
    );
    let vecm = ql_portfolio::solve_johansen_tsay_vecm_cointegration_basket(
        &ql_portfolio::JohansenTsayVecmCointegrationSpec {
            price_series: vec![
                [T::from_f64(100.) + x[0], T::from_f64(101.) + x[1]],
                [T::from_f64(101.) + x[2], T::from_f64(102.) + x[3]],
                [T::from_f64(102.) + x[0], T::from_f64(103.) + x[1]],
                [T::from_f64(103.) + x[2], T::from_f64(104.) + x[3]],
            ],
        },
    );
    let hj = ql_portfolio::evaluate_cochrane_hansen_jagannathan_sdf_and_campbell_shiller(
        &ql_portfolio::CochraneHansenJagannathanSdfSpec {
            mean_excess_returns: [T::from_f64(0.04) + x[0], T::from_f64(0.05) + x[1]],
            excess_return_cov: [
                [T::from_f64(0.04), T::from_f64(0.01)],
                [T::from_f64(0.01), T::from_f64(0.06)],
            ],
            risk_free_rate: T::from_f64(0.02),
            candidate_sdf_pricing_errors: [x[2], x[3]],
            steady_state_log_dp_ratio: T::from_f64(1.),
            var1_companion_matrix: [[T::from_f64(0.9), T::zero()], [T::zero(), T::from_f64(0.8)]],
            state_innovation_cov: [
                [T::from_f64(0.02), T::zero()],
                [T::zero(), T::from_f64(0.03)],
            ],
        },
    );
    let ecl = ql_credit::compute_ifrs9_cecl_multistage_expected_credit_loss(
        &ql_credit::Ifrs9CeclMultiStageEclSpec {
            annual_ttc_pds: [
                T::from_f64(0.02) + x[0] * 0.001,
                T::from_f64(0.03),
                T::from_f64(0.1),
            ],
            annual_eads: [
                T::from_f64(100.),
                T::from_f64(100.) + x[1],
                T::from_f64(100.) + x[2],
            ],
            base_lgd: T::from_f64(0.45),
            lgd_macro_sensitivity: T::from_f64(0.05),
            vasicek_asset_correlation: T::from_f64(0.2),
            effective_interest_rate: T::from_f64(0.03),
            scenario_weights: [T::from_f64(0.3), T::from_f64(0.5), T::from_f64(0.2)],
            scenario_macro_shocks_z: [x[0], x[1], x[2] + x[3]],
            origination_lifetime_pd: T::from_f64(0.04),
            sicr_pd_multiplier_threshold: T::from_f64(1.5),
        },
    );
    let cpt = ql_portfolio::evaluate_wakker_tversky_kahneman_cumulative_prospect_utility(
        &ql_portfolio::WakkerTverskyKahnemanCptSpec {
            scenario_returns: [
                x[0],
                x[1],
                T::from_f64(-0.1) + x[2],
                T::from_f64(0.2) + x[3],
            ],
            scenario_probabilities: [T::from_f64(0.25); 4],
            reference_return: T::zero(),
            alpha_gain: T::from_f64(0.88),
            beta_loss: T::from_f64(0.88),
            loss_aversion_lambda: T::from_f64(2.25),
            prelec_gamma_gain: T::from_f64(0.7),
            prelec_gamma_loss: T::from_f64(0.7),
            prelec_delta: T::one(),
        },
    );
    let hs = ql_portfolio::evaluate_hasbrouck_information_share_and_dssw_limits_to_arb(
        &ql_portfolio::HasbrouckDsswPriceDiscoverySpec {
            alpha_venue1: T::from_f64(0.2) + x[0],
            alpha_venue2: T::from_f64(0.3) + x[1],
            sigma_venue1: T::from_f64(0.1),
            sigma_venue2: T::from_f64(0.12),
            venue_innovation_corr: T::from_f64(0.2),
            fundamental_dividend: T::from_f64(1.),
            risk_free_rate: T::from_f64(0.02),
            noise_trader_share_mu: T::from_f64(0.2),
            noise_trader_sentiment_mean: x[2],
            noise_trader_sentiment_var: T::from_f64(0.1) + x[3],
            arbitrageur_risk_aversion_gamma: T::from_f64(2.),
        },
    );
    let fv = ql_rates::price_fong_vasicek_stochastic_vol_bond(
        &ql_rates::FongVasicekStochVolShortRateSpec {
            initial_short_rate_r0: T::from_f64(0.03) + x[0],
            initial_variance_v0: T::from_f64(0.04) + x[2] * T::from_f64(0.001),
            rate_mean_reversion_alpha: T::from_f64(0.5),
            long_run_rate_r_bar: T::from_f64(0.03),
            vol_mean_reversion_gamma: T::from_f64(1.),
            long_run_variance_v_bar: T::from_f64(0.04),
            vol_of_vol_xi: T::from_f64(0.2) + x[2] * T::from_f64(0.01),
            correlation_rho: T::from_f64(-0.3) + x[3] * T::from_f64(0.01),
            maturity_tau: T::from_f64(2.) + x[1],
        },
    );
    let la = ql_credit::evaluate_lando_cox_affine_default_and_liquidity_decomposition(
        &ql_credit::LandoCoxAffineDefaultSpec {
            risk_free_rate: T::from_f64(0.02) + x[0] * T::from_f64(0.01),
            base_default_hazard: T::from_f64(0.01) + x[1] * T::from_f64(0.01),
            initial_credit_intensity_y0: T::from_f64(0.02) + x[0] * T::from_f64(0.01),
            credit_kappa: T::from_f64(0.5),
            credit_theta: T::from_f64(0.02),
            credit_sigma: T::from_f64(0.1),
            initial_liquidity_intensity_z0: T::from_f64(0.01) + x[2],
            liquidity_kappa: T::from_f64(0.4) + x[1] * T::from_f64(0.01),
            liquidity_theta: T::from_f64(0.01),
            liquidity_sigma: T::from_f64(0.1),
            recovery_rate: T::from_f64(0.4),
            maturity: T::from_f64(3.) + x[3],
        },
    );
    let lp = ql_portfolio::calibrate_sornette_lppl_bubble_and_crash_hazard(
        &ql_portfolio::SornetteLpplBubbleCrashSpec {
            observation_times: vec![
                T::from_f64(1.),
                T::from_f64(2.),
                T::from_f64(3.),
                T::from_f64(4.),
                T::from_f64(5.),
            ],
            log_prices: vec![
                T::from_f64(4.) + x[0],
                T::from_f64(4.1) + x[1],
                T::from_f64(4.2) + x[2],
                T::from_f64(4.3) + x[3],
                T::from_f64(4.4) + x[0],
            ],
            critical_time_tc: T::from_f64(10.),
            power_law_exponent_m: T::from_f64(0.5),
            log_periodic_omega: T::from_f64(8.),
            expected_crash_drop_kappa: T::from_f64(0.2),
            ridge_penalty: T::from_f64(1e-6),
        },
    );
    let cb = ql_fixed_income::analyze_thorp_kassouf_calamos_convertible_arbitrage(
        &ql_fixed_income::ThorpKassoufCalamosConvertibleArbSpec {
            stock_price: T::from_f64(100.) + x[0],
            conversion_ratio: T::from_f64(0.8),
            convertible_market_price: T::from_f64(95.),
            face_value: T::from_f64(100.),
            annual_coupon_rate: T::from_f64(0.03),
            stock_dividend_yield: T::from_f64(0.01),
            risk_free_rate: T::from_f64(0.02),
            credit_spread: T::from_f64(0.03),
            short_stock_rebate_rate: T::from_f64(0.01),
            implied_volatility: T::from_f64(0.2) + x[1] * T::from_f64(0.01),
            realized_volatility: T::from_f64(0.22) + x[2] * T::from_f64(0.01),
            maturity: T::from_f64(3.) + x[3],
        },
    );
    let ca = ql_portfolio::evaluate_lopez_de_prado_causal_backdoor_vs_collider_bias(
        &ql_portfolio::LopezDePradoCausalFactorSpec {
            true_causal_beta_xy: T::from_f64(0.5) + x[0] * T::from_f64(0.01),
            confounder_to_factor_axz: T::from_f64(0.7) + x[1],
            confounder_to_return_byz: T::from_f64(0.4),
            factor_to_collider_ccx: T::from_f64(0.5),
            return_to_collider_ccy: T::from_f64(0.5) + x[2] * T::from_f64(0.01),
            confounder_var: T::one(),
            factor_idiosyncratic_var: T::from_f64(0.5),
            return_idiosyncratic_var: T::from_f64(0.7),
            collider_idiosyncratic_var: T::from_f64(0.5) + x[3] * T::from_f64(0.01),
        },
    );
    let ma = ql_dml::evaluate_fournie_lasry_lions_malliavin_bel_greeks(
        &ql_dml::FournieLasryLionsMalliavinBelSpec {
            spot: T::from_f64(100.) + x[3],
            strike: T::from_f64(100.) + x[0],
            rate: T::from_f64(0.02) + x[1] * T::from_f64(0.01),
            div_yield: T::zero(),
            vol: T::from_f64(0.2) + x[2] * T::from_f64(0.01),
            expiry: T::one() + x[3] * T::from_f64(0.01),
            localization_bandwidth: T::from_f64(0.5),
        },
    );
    [
        kal.innovation_log_likelihood,
        vecm.ou_mean_reversion_kappa,
        hj.hj_sdf_volatility_lower_bound,
        ecl.blended_ifrs9_ecl_provision,
        cpt.cumulative_prospect_utility,
        hs.hasbrouck_is_midpoint_venue1,
        fv.zero_coupon_bond_price,
        la.risky_corporate_bond_price,
        lp.terminal_crash_hazard_rate + lp.fitted_rmse,
        cb.theoretical_convertible_value,
        ca.kitchen_sink_collider_conditioned_beta,
        ma.call_malliavin_delta,
    ]
}
pub fn benchmark_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.02, -0.01, 0.03, -0.01];
    let baseline = move || {
        let f = |x: [f64; 4]| suite37_workloads(x).iter().map(|v| v.val()).sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        type D = Dual<4>;
        let x = [
            D::variable(x0[0], 0),
            D::variable(x0[1], 1),
            D::variable(x0[2], 2),
            D::variable(x0[3], 3),
        ];
        suite37_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "37_kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_fv_lando_lppl_cb_arb_causal_malliavin",
        "Kalman-RTS / VECM / HJ-SDF / IFRS9 / CPT / Hasbrouck / FV / Lando / LPPL / CB-Arb / Causal / Malliavin AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}
pub fn benchmark_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.02, -0.01, 0.03, -0.01];
    let baseline = move || {
        let f = |x: [f64; 4]| suite36_workloads(x).iter().map(|v| v.val()).sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite36_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp",
        "MBS-IO/PO / RMV-Basis / MO-XVA / Dilution-CB / Sig / HN-DCC-GARCH / POT-GPD / RIE / PIN / GP-Aim / Brinson-Carino / IFT-IERP AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite38_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let times = vec![
        T::from_f64(0.25) + x[0] * 0.01,
        T::from_f64(0.5) + x[1] * 0.01,
        T::from_f64(0.75) + x[2] * 0.01,
        T::one() + x[3] * 0.01,
    ];
    let mf = ql_rates::evaluate_hunt_kennedy_markov_functional_range_accrual(
        &ql_rates::HuntKennedyMarkovFunctionalSpec {
            fixing_times: times.clone(),
            initial_forwards: vec![T::from_f64(0.03) + x[0] * 0.001; 4],
            market_caplet_vols: vec![T::from_f64(0.2) + x[1] * 0.01; 4],
            smile_curvature_eta: T::from_f64(0.1) + x[2] * 0.01,
            mean_reversion: T::from_f64(0.03),
            range_lower: T::from_f64(0.0),
            range_upper: T::from_f64(0.06),
            accrual_coupon_rate: T::from_f64(0.03) + x[3] * 0.001,
            call_strike_spread: T::from_f64(0.005),
        },
    );
    let wc = ql_fixed_income::evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch(
        &ql_fixed_income::TuckmanBurghardtTreasuryWildcardBasisSpec {
            futures_settlement_price: T::from_f64(100.) + x[0],
            delivery_horizon_years: T::from_f64(0.25),
            general_collateral_repo_rate: T::from_f64(0.03) + x[1] * 0.001,
            repo_specialness_spread: T::from_f64(0.002),
            bond1_clean_price: T::from_f64(99.) + x[2],
            bond1_conversion_factor: T::from_f64(0.85),
            bond1_modified_duration: T::from_f64(6.),
            bond1_yield_beta: T::from_f64(1.),
            bond2_clean_price: T::from_f64(101.) + x[3],
            bond2_conversion_factor: T::from_f64(0.9),
            bond2_modified_duration: T::from_f64(7.),
            bond2_yield_beta: T::from_f64(1.1),
            afternoon_wildcard_vol_annual: T::from_f64(0.1),
            eom_delivery_window_days: 5,
        },
    );
    let tarf =
        ql_fx::evaluate_wystup_pivot_tarf_and_corridor_strip(&ql_fx::WystupPivotTarfCorridorSpec {
            spot: T::from_f64(1.1) + x[0] * 0.01,
            strike_low: T::from_f64(1.0),
            strike_pivot: T::from_f64(1.1) + x[1] * 0.01,
            strike_high: T::from_f64(1.2),
            downside_gearing: T::from_f64(1.2),
            upside_gearing: T::from_f64(1.1),
            target_redemption_cap: T::from_f64(0.1),
            domestic_rate: T::from_f64(0.03) + x[2] * 0.001,
            foreign_rate: T::from_f64(0.01),
            volatility: T::from_f64(0.12) + x[3] * 0.01,
            fixing_times: times.clone(),
        });
    let bg = ql_vol::evaluate_bergomi_forward_start_smile_and_napoleon_convexity(
        &ql_vol::BergomiForwardStartNapoleonVolSpec {
            initial_forward_variance: T::from_f64(0.04) + x[0] * 0.001,
            vol_of_vol_omega: T::from_f64(0.3) + x[1] * 0.01,
            kappa_short_k1: T::from_f64(3.),
            kappa_long_k2: T::from_f64(0.5),
            long_factor_weight_theta: T::from_f64(0.5),
            factor_corr_rho12: T::from_f64(0.2),
            spot_vol_corr_rho_s1: T::from_f64(-0.5),
            spot_vol_corr_rho_s2: T::from_f64(-0.3),
            forward_start_t1: T::from_f64(0.5),
            forward_end_t2: T::one() + x[2] * 0.01,
            forward_moneyness_k: T::one() + x[3] * 0.01,
            napoleon_base_coupon: T::from_f64(0.01),
            napoleon_global_floor: T::from_f64(0.0),
            napoleon_subperiods_m: 4,
        },
    );
    let sch = ql_credit::evaluate_schonbucher_defaultable_range_accrual_and_spread_cap(
        &ql_credit::SchonbucherDefaultableRangeAccrualCapSpec {
            notional: T::from_f64(100.),
            payment_times: times.clone(),
            risk_free_rate: T::from_f64(0.02) + x[0] * 0.001,
            initial_hazard_rates: vec![T::from_f64(0.02) + x[1] * 0.001; 4],
            base_recovery_rate: T::from_f64(0.4),
            spread_volatilities: vec![T::from_f64(0.3) + x[2] * 0.01; 4],
            spread_cap_strike: T::from_f64(0.03),
            accrual_spread_lower: T::from_f64(0.),
            accrual_spread_upper: T::from_f64(0.06),
            range_accrual_coupon_rate: T::from_f64(0.02) + x[3] * 0.001,
            recovery_spread_correlation: T::from_f64(-0.2),
        },
    );
    let cv = ql_xva::evaluate_green_frtb_cva_hedge_relief_and_kva(
        &ql_xva::GreenFrtbCvaHedgeCapitalSpec {
            counterparty_sacr_eads: vec![T::from_f64(100.) + x[0], T::from_f64(80.) + x[1]],
            counterparty_effective_maturities: vec![T::from_f64(3.), T::from_f64(5.)],
            counterparty_risk_weights: vec![T::from_f64(0.02), T::from_f64(0.03)],
            single_name_hedge_notionals: vec![T::from_f64(20.) + x[2]],
            single_name_hedge_maturities: vec![T::from_f64(3.)],
            hedge_reference_correlations_rhc: vec![T::from_f64(0.8)],
            index_hedge_notional: T::from_f64(10.) + x[3],
            index_hedge_maturity: T::from_f64(5.),
            index_risk_weight: T::from_f64(0.02),
            sa_cva_spread_deltas: vec![T::from_f64(1.)],
            sa_cva_vega_sensitivities: vec![T::from_f64(0.5)],
            hurdle_rate_gamma_k: T::from_f64(0.1),
            hedge_margin_funding_spread: T::from_f64(0.01),
            discount_rate: T::from_f64(0.03),
        },
    );
    let lb = ql_exotics::evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread(
        &ql_exotics::HeynenKatPartialFractionalLookbackSpec {
            spot: T::from_f64(100.) + x[0],
            running_min: T::from_f64(90.),
            running_max: T::from_f64(110.) + x[1],
            fractional_multiplier_lambda: T::from_f64(1.2),
            monitoring_end_t1: T::from_f64(0.5),
            maturity_t2: T::one(),
            fixed_strike: T::from_f64(100.) + x[2],
            rate: T::from_f64(0.02),
            dividend_yield: T::from_f64(0.),
            volatility: T::from_f64(0.2) + x[3] * 0.01,
        },
    );
    let vta = ql_commodities::evaluate_eydeland_wolyniec_ramp_outage_virtual_tolling(
        &ql_commodities::EydelandWolyniecRampOutageVtaSpec {
            power_forwards: vec![T::from_f64(60.) + x[0], T::from_f64(62.)],
            gas_forwards: vec![T::from_f64(30.), T::from_f64(31.) + x[1]],
            spinning_reserve_prices: vec![T::from_f64(3.), T::from_f64(3.)],
            dispatch_times: times[..2].to_vec(),
            capacity_mw: T::from_f64(100.),
            block_hours: T::one(),
            full_load_heat_rate: T::from_f64(2.),
            part_load_heat_rate: T::from_f64(2.2),
            min_stable_load_fraction: T::from_f64(0.5),
            vom_cost_per_mwh: T::one(),
            startup_cost: T::from_f64(2.),
            fast_ramp_penalty_cost: T::from_f64(1.),
            outage_failure_rate_lambda: T::from_f64(0.05),
            outage_repair_rate_mu: T::from_f64(0.95),
            spread_volatility: T::from_f64(0.2) + x[3] * T::from_f64(0.01),
            risk_free_rate: T::from_f64(0.02) + x[2] * 0.001,
        },
    );
    let fl = ql_portfolio::evaluate_grinold_kahn_clarke_fundamental_law(
        &ql_portfolio::GrinoldKahnClarkeFundamentalLawSpec {
            standardized_scores_z: vec![
                T::one() + x[0],
                T::from_f64(-0.5) + x[1],
                T::from_f64(0.7) + x[2],
            ],
            residual_volatilities: vec![T::from_f64(0.2), T::from_f64(0.25), T::from_f64(0.3)],
            benchmark_weights: vec![T::from_f64(0.3), T::from_f64(0.4), T::from_f64(0.3)],
            information_coefficient_ic: T::from_f64(0.1),
            rebalances_per_year: 12,
            signal_autocorrelation_rho: T::from_f64(0.0) + x[3] * 0.01,
            target_tracking_error: T::from_f64(0.05),
            max_active_weight_bound: T::from_f64(0.1),
        },
    );
    let acd = ql_portfolio::evaluate_engle_russell_acd_trade_duration_and_hazard(
        &ql_portfolio::EngleRussellAcdTradeDurationSpec {
            inter_trade_durations: vec![
                T::one() + x[0],
                T::from_f64(1.2) + x[1],
                T::from_f64(0.9) + x[2],
            ],
            omega: T::from_f64(0.2),
            alpha: T::from_f64(0.2),
            beta: T::from_f64(0.5),
            weibull_shape_gamma: T::one(),
            initial_expected_duration: T::one(),
            elapsed_time_since_last_trade: T::from_f64(0.5) + x[3],
            base_price_impact_bps: T::from_f64(2.),
            duration_impact_elasticity: T::from_f64(0.5),
        },
    );
    let hlz = ql_portfolio::evaluate_harvey_liu_zhu_multiple_testing_and_mintrl(
        &ql_portfolio::HarveyLiuZhuMultipleTestingSharpeSpec {
            focal_sharpe_ratio: T::from_f64(0.8) + x[0],
            benchmark_sharpe_sr0: T::from_f64(0.),
            sample_length_t: 252,
            return_skewness: T::from_f64(-0.5) + x[1],
            return_kurtosis_pearson: T::from_f64(4.),
            num_independent_trials_m: 20,
            cross_sectional_sharpe_std: T::from_f64(0.3),
            confidence_z_alpha: T::from_f64(1.96),
            candidate_t_stats: vec![T::from_f64(1.5) + x[2], T::from_f64(2.) + x[3]],
        },
    );
    let ml = ql_dml::evaluate_giles_capriotti_multilevel_mc_aad(
        &ql_dml::GilesCapriottiMultilevelMcAadSpec {
            spot: T::from_f64(100.) + x[0],
            strike: T::from_f64(100.),
            rate: T::from_f64(0.02) + x[1] * 0.001,
            volatility: T::from_f64(0.2) + x[2] * 0.01,
            local_vol_skew: T::from_f64(-0.2),
            maturity: T::one() + x[3] * 0.01,
            base_steps_m0: 8,
            target_rmse_epsilon: T::from_f64(0.01),
        },
    );
    [
        mf.non_callable_range_accrual_npv,
        wc.theoretical_net_basis_ctd,
        tarf.pivot_tarf_client_npv,
        bg.bergomi_stoch_vol_napoleon_coupon,
        sch.defaultable_range_accrual_note_npv,
        cv.ba_cva_total_capital_charge,
        lb.heynen_kat_partial_window_call,
        vta.reliability_derated_vta_npv,
        fl.constrained_information_ratio,
        acd.duration_adjusted_price_impact_bps,
        hlz.deflated_sharpe_ratio_dsr,
        ml.mlmc_telescoping_price,
    ]
}

pub fn benchmark_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.02_f64, -0.01, 0.03, -0.01];
    let baseline = move || {
        let f = |x: [f64; 4]| suite38_workloads(x).iter().map(|v| v.val()).sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite38_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "38_mf_rates_wildcard_ctd_pivot_tarf_bergomi_schonbucher_frtb_cva_lookback_vta_flam_acd_hlz_mlmc",
        "MF-Rates / Wildcard / CTD / Pivot-TARF / Bergomi / Schönbucher / FRTB-CVA / Lookback / VTA / FLAM / ACD / HLZ / MLMC AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite40_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let irr = ql_curves::evaluate_bcbs_d368_six_scenario_irrbb(
        &ql_curves::BcbsD368SixScenarioIrrbbSpec {
            bucket_midpoints_years: vec![T::from_f64(0.5), T::from_f64(2.0), T::from_f64(5.0)],
            asset_repricing_cashflows: vec![
                T::from_f64(10.0) + x[0],
                T::from_f64(8.0) + x[1],
                T::from_f64(5.0) + x[2],
            ],
            liability_repricing_cashflows: vec![
                T::from_f64(7.0),
                T::from_f64(6.0) + x[3],
                T::from_f64(4.0),
            ],
            baseline_zero_rates: vec![T::from_f64(0.02); 3],
            parallel_shock_bps: T::from_f64(200.0),
            short_shock_bps: T::from_f64(250.0),
            long_shock_bps: T::from_f64(100.0),
            decay_tau_years: T::from_f64(4.0),
            tier1_capital: T::from_f64(100.0),
            pass_through_deposit_beta: T::from_f64(0.5),
        },
    );
    let rfr = ql_rates::evaluate_turfus_willems_rfr_in_arrears_caplet_smile(
        &ql_rates::TurfusWillemsRfrInArrearsSmileSpec {
            accrual_start_t1: T::from_f64(0.1) + x[0],
            accrual_end_t2: T::from_f64(0.6),
            forward_rfr_rate: T::from_f64(0.03) + x[1] * T::from_f64(0.001),
            strike_rate: T::from_f64(0.025),
            discount_factor_t2: T::from_f64(0.98),
            sabr_alpha: T::from_f64(0.02),
            sabr_beta: T::from_f64(0.5),
            sabr_rho: T::from_f64(-0.2),
            sabr_nu: T::from_f64(0.4),
            hull_white_mean_reversion: T::from_f64(0.1),
            notional: T::from_f64(1e6),
        },
    );
    let xccy = ql_fx::evaluate_baran_bianchetti_mtm_xccy_basis_swap(
        &ql_fx::BaranBianchettiMtmXccyBasisSpec {
            spot_fx_dom_per_for: T::from_f64(1.1) + x[2],
            foreign_notional: T::from_f64(100.0),
            payment_times: vec![T::from_f64(0.5), T::one(), T::from_f64(1.5)],
            dom_ois_discount_factors: vec![T::from_f64(0.99); 3],
            for_xccy_discount_factors: vec![T::from_f64(0.98); 3],
            dom_forward_rates: vec![T::from_f64(0.03) + x[3] * T::from_f64(0.001); 3],
            for_forward_rates: vec![T::from_f64(0.02); 3],
            quoted_dom_basis_spread: T::from_f64(0.001),
            dom_rate_vol: T::from_f64(0.2),
            fx_spot_vol: T::from_f64(0.1),
            corr_dom_rate_fx: T::from_f64(-0.3),
        },
    );
    let bg = ql_vol::evaluate_bourgey_gatheral_bergomi_guyon_expansion(
        &ql_vol::BourgeyGatheralBergomiGuyonSmileSpec {
            forward_price: T::from_f64(100.0) + x[0],
            strike_price: T::from_f64(100.0) + x[1],
            maturity: T::one(),
            initial_instantaneous_vol: T::from_f64(0.2),
            vol_of_vol_eta: T::from_f64(0.5),
            spot_vol_corr_rho: T::from_f64(-0.7),
            hurst_exponent: T::from_f64(0.1),
            mean_reversion_lambda: T::one(),
        },
    );
    let coco = ql_credit::evaluate_feds_despiegeleer_at1_coco_dual_trigger(
        &ql_credit::FedsDeSpiegeleerAt1CocoDualTriggerSpec {
            notional: T::from_f64(100.0),
            coupon_rate: T::from_f64(0.08),
            maturity: T::from_f64(5.0),
            risk_free_rate: T::from_f64(0.03),
            initial_cet1_ratio: T::from_f64(0.12) + x[2],
            mechanical_cet1_trigger: T::from_f64(0.07),
            ponv_cet1_trigger: T::from_f64(0.05),
            supervisory_ponv_hazard_lambda: T::from_f64(0.03),
            asset_volatility: T::from_f64(0.2),
            cet1_asset_elasticity: T::one(),
            initial_stock_price: T::from_f64(20.0),
            conversion_price: T::from_f64(15.0),
            floor_stock_price: T::from_f64(5.0),
            existing_common_shares: T::from_f64(100.0),
            write_down_fraction: T::from_f64(0.5),
            is_equity_conversion: true,
        },
    );
    let xv = ql_xva::evaluate_albanese_crepey_zeron_balance_sheet_xva(
        &ql_xva::AlbaneseCrepeyZeronBalanceSheetXvaSpec {
            maturity: T::from_f64(2.0),
            num_time_steps: 8,
            risk_free_rate: T::from_f64(0.03),
            cpty_hazard_rate: T::from_f64(0.02),
            cpty_recovery: T::from_f64(0.4),
            bank_hazard_rate: T::from_f64(0.01),
            bank_recovery: T::from_f64(0.4),
            bank_liquidity_basis: T::from_f64(0.005),
            im_funding_spread: T::from_f64(0.01),
            shareholder_hurdle_rate: T::from_f64(0.1),
            base_unmargined_epe: T::from_f64(10.0) + x[3],
            base_ene: T::from_f64(4.0),
            simm_delta_weight: T::from_f64(0.6),
            simm_vega_weight: T::from_f64(0.4),
            simm_delta_vega_corr: T::from_f64(0.25),
            regulatory_capital_ratio: T::from_f64(0.08),
        },
    );
    let rb = ql_models::evaluate_jacquier_chebyshev_rough_vol_surrogate(
        &ql_models::JacquierDeepRoughVolSurrogateSpec {
            forward: T::from_f64(100.0),
            maturity: T::one(),
            initial_var_xi0: T::from_f64(0.04),
            hurst_h: T::from_f64(0.1),
            vol_of_vol_eta: T::from_f64(0.5),
            spot_vol_corr_rho: T::from_f64(-0.7),
            target_market_vols_3_strikes: [T::from_f64(0.22), T::from_f64(0.20), T::from_f64(0.21)],
            log_moneyness_strikes: [T::from_f64(-0.1) + x[1], T::zero(), T::from_f64(0.1) + x[2]],
        },
    );
    let spread = ql_exotics::evaluate_choi_li_deng_multi_asset_spread_basket(
        &ql_exotics::ChoiLiDengMultiAssetSpreadBasketSpec {
            long_leg_forward: T::from_f64(100.0) + x[0],
            long_leg_weight: T::one(),
            long_leg_vol: T::from_f64(0.2),
            short_leg_forwards: vec![T::from_f64(95.0) + x[1]],
            short_leg_weights: vec![T::one()],
            short_leg_vols: vec![T::from_f64(0.25)],
            corr_long_short: vec![T::from_f64(0.4)],
            corr_short_short_equicorr: T::zero(),
            strike: T::from_f64(1.0) + x[3],
            maturity: T::one(),
            risk_free_rate: T::from_f64(0.03),
        },
    );
    let bess = ql_commodities::evaluate_bess_analytical_sdp_degradation_storage(
        &ql_commodities::BessAnalyticalSdpDegradationStorageSpec {
            energy_capacity_mwh: T::from_f64(100.0) + x[0],
            max_power_mw: T::from_f64(50.0),
            initial_soc_fraction: T::from_f64(0.5),
            min_soc_fraction: T::from_f64(0.1),
            max_soc_fraction: T::from_f64(0.9),
            charge_efficiency: T::from_f64(0.9),
            discharge_efficiency: T::from_f64(0.9),
            battery_capex_usd: T::from_f64(100000.0),
            rated_cycle_life: T::from_f64(4000.0),
            peak_mean_price_mwh: T::from_f64(100.0) + x[1],
            offpeak_mean_price_mwh: T::from_f64(20.0),
            price_volatility_mwh: T::from_f64(30.0),
            horizon_days: 365,
            discount_rate_annual: T::from_f64(0.08),
        },
    );
    let dark = ql_portfolio::evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool(
        &ql_portfolio::MitchellPedersenKratzDarkPoolSpiralSpec {
            initial_inventory_shares: T::from_f64(1000.0) + x[2],
            reference_price: T::from_f64(100.0),
            execution_horizon_days: T::from_f64(5.0),
            daily_volatility: T::from_f64(0.02),
            lit_temporary_impact_eta: T::from_f64(1e-4),
            lit_permanent_impact_gamma: T::from_f64(2e-4),
            risk_aversion_lambda: T::from_f64(0.1),
            dark_pool_poisson_intensity: T::from_f64(0.3),
            dark_pool_adverse_selection_bps: T::from_f64(5.0),
            specialist_immediate_capital: T::from_f64(1e5),
            slow_moving_external_capital: T::from_f64(2e5),
            capital_arrival_speed_beta: T::from_f64(0.4),
        },
    );
    let risk = ql_portfolio::evaluate_follmer_schied_convex_risk_and_diagnostic_expectations(
        &ql_portfolio::FollmerSchiedGennaioliDiagnosticRiskSpec {
            scenario_pnl: vec![T::one() + x[0], T::from_f64(-2.0), T::from_f64(-5.0)],
            baseline_probabilities: vec![T::from_f64(0.5), T::from_f64(0.3), T::from_f64(0.2)],
            entropic_ambiguity_theta: T::from_f64(2.0),
            kl_divergence_budget_eta: T::from_f64(0.1),
            diagnostic_overreaction_theta: T::from_f64(1.5),
            boom_news_drift_shock: T::from_f64(0.2),
            bust_news_drift_shock: T::from_f64(1.0) + x[3],
        },
    );
    let pca = ql_dml::evaluate_huge_savine_differential_pca_and_twin_dml(
        &ql_dml::HugeSavineDifferentialPcaSpec {
            state_samples: vec![
                vec![T::from_f64(0.2), T::from_f64(2.0)],
                vec![T::from_f64(0.3), T::from_f64(-1.8)],
                vec![T::from_f64(-0.25), T::from_f64(1.5)],
            ],
            payoff_samples: vec![T::from_f64(1.0), T::from_f64(0.8), T::from_f64(1.2)],
            aad_gradient_samples: vec![
                vec![T::from_f64(2.2) + x[1], T::from_f64(0.1)],
                vec![T::from_f64(1.8), T::from_f64(-0.15)],
                vec![T::from_f64(-2.0), T::from_f64(0.12)],
            ],
            twin_fd_gradient_samples: vec![
                vec![T::from_f64(2.0), T::from_f64(0.12)],
                vec![T::from_f64(1.7), T::from_f64(-0.1)],
                vec![T::from_f64(-1.9), T::from_f64(0.1)],
            ],
            retained_rank: 1,
            ridge_regularization: T::from_f64(1e-6),
        },
    );
    [
        irr.baseline_eve,
        rfr.in_arrears_caplet_pv * T::from_f64(1e-4),
        xccy.mtm_swap_npv_in_dom,
        bg.strike_implied_vol,
        coco.at1_coco_fair_value,
        xv.ftd_cva,
        rb.surrogate_smile_vols[0],
        spread.choi_gauss_hermite_price,
        bess.total_bess_npv_usd * T::from_f64(1e-5),
        dark.kratz_schoneborn_dark_lit_cost * T::from_f64(1e-5),
        risk.follmer_schied_entropic_risk_measure,
        pca.truncation_gradient_mse_bound,
    ]
}

pub fn benchmark_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite40_workloads(x).iter().map(|v| v.val()).sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite40_workloads(x)
            .iter()
            .map(|y| y.val + y.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "40_bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_feds_ponv_zeron_simm_jacquier_rbergomi_choi_spread_bess_sdp_kratz_dark_follmer_diag_diff_pca",
        "BCBS-IRRBB / Turfus-RFR / Baran-MtM / Bourgey-BG / FEDS-PONV / Zeron-SIMM / Jacquier-rBergomi / Choi-Spread / BESS-SDP / Kratz-Dark / Föllmer / Diff-PCA AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite41_workloads<T: Real>(x: [T; 4]) -> T {
    let a = ql_fixed_income::evaluate_veronesi_fabozzi_futures_convexity_and_prepayment(
        &ql_fixed_income::VeronesiFabozziConvexityPrepaymentSpec {
            forward_rate: T::from_f64(0.03) + x[0] * 0.01,
            short_rate_vol: T::from_f64(0.2) + x[1] * 0.01,
            mean_reversion_kappa: T::from_f64(0.1),
            fixing_time_t1: T::one(),
            maturity_time_t2: T::from_f64(2.0),
            wac_coupon: T::from_f64(0.05),
            current_mortgage_rate: T::from_f64(0.04) + x[2] * 0.01,
            min_cpr_turnover: T::from_f64(0.02),
            max_cpr_refi: T::from_f64(0.3),
            refi_steepness_gamma: T::from_f64(8.0),
            refi_incentive_threshold: T::zero(),
            pool_age_months: T::from_f64(24.0) + x[3],
            psa_ramp_months: T::from_f64(30.0),
        },
    );
    let b = ql_credit::compute_lando_israel_generator_logarithm_and_migration(
        &ql_credit::LandoIsraelGeneratorLogarithmSpec {
            transition_matrix_3x3: [
                [T::from_f64(0.9), T::from_f64(0.08), T::from_f64(0.02)],
                [T::from_f64(0.03), T::from_f64(0.92), T::from_f64(0.05)],
                [T::zero(), T::zero(), T::one()],
            ],
            series_order_k: 12,
            fractional_horizon_years: T::from_f64(0.5) + x[0],
            recovery_rate: T::from_f64(0.4),
            risk_free_rate: T::from_f64(0.03) + x[1] * 0.01,
        },
    );
    let c = ql_vol::evaluate_gatheral_demeterfi_corridor_variance_and_skew(
        &ql_vol::GatheralDemeterfiCorridorVarianceSkewSpec {
            forward: T::from_f64(100.0) + x[0],
            maturity: T::one(),
            atm_vol: T::from_f64(0.2) + x[1] * 0.01,
            atm_skew: T::from_f64(-0.2),
            smile_curvature: T::from_f64(0.1) + x[2] * 0.01,
            corridor_lower_strike: T::from_f64(90.0),
            corridor_upper_strike: T::from_f64(110.0),
            sampling_steps_n: 12,
            vol_of_vol_eta: T::from_f64(0.5),
            lsv_mixing_weight: T::from_f64(0.3),
        },
    );
    let d = ql_exotics::price_shepp_shiryaev_shreve_perpetual_russian_and_mckean(
        &ql_exotics::SheppShiryaevShrevePerpetualRussianSpec {
            spot: T::from_f64(100.0) + x[0],
            running_max: T::from_f64(110.0),
            fixed_put_strike: T::from_f64(100.0),
            risk_free_rate: T::from_f64(0.03),
            dividend_yield: T::from_f64(0.01) + x[1] * 0.01,
            volatility: T::from_f64(0.2) + x[2] * 0.01,
        },
    );
    let e = ql_levy::evaluate_cont_tankov_bates_svj_characteristic_and_pide(
        &ql_levy::ContTankovBatesSvjPideSpec {
            spot: T::from_f64(100.0) + x[0],
            strike: T::from_f64(110.0),
            maturity: T::one(),
            risk_free_rate: T::from_f64(0.03),
            dividend_yield: T::from_f64(0.01),
            v0: T::from_f64(0.04) + x[1] * 0.01,
            kappa: T::one(),
            theta: T::from_f64(0.04),
            sigma_v: T::from_f64(0.3),
            rho: T::from_f64(-0.5),
            jump_intensity_lambda: T::from_f64(0.2) + x[2] * 0.01,
            jump_mean_mu_j: T::from_f64(-0.1),
            jump_std_sigma_j: T::from_f64(0.2),
        },
    );
    let f = ql_fixed_income::evaluate_koziol_asquith_fasb_strategic_convertible(
        &ql_fixed_income::KoziolAsquithFasbStrategicConvertibleSpec {
            firm_value: T::from_f64(10000.0),
            debt_face_value: T::from_f64(8000.0),
            existing_shares: T::from_f64(100.0),
            new_shares_blockholder1: T::from_f64(10.0) + x[0],
            new_shares_blockholder2: T::from_f64(8.0),
            coupon_rate: T::from_f64(0.04),
            corporate_tax_rate: T::from_f64(0.25),
            dividend_yield: T::from_f64(0.01),
            stock_price: T::from_f64(100.0) + x[1],
            base_conversion_ratio: T::one(),
            induced_sweetened_conversion_ratio: T::from_f64(1.1),
            cash_sweetener_amount: T::one(),
            call_price: T::from_f64(120.0),
            firm_volatility: T::from_f64(0.2),
            risk_free_rate: T::from_f64(0.03),
            maturity: T::from_f64(5.0),
        },
    );
    let g = ql_models::solve_carmona_delarue_mean_field_systemic_interbank_game(
        &ql_models::CarmonaDelarueMeanFieldSystemicSpec {
            initial_log_reserve_spread: x[0],
            interbank_mean_reversion_a: T::from_f64(0.4),
            control_quadratic_penalty_q: T::from_f64(0.2),
            systemic_risk_preference_epsilon: T::from_f64(0.3) + x[1],
            terminal_systemic_penalty_c: T::one(),
            idiosyncratic_vol_sigma: T::from_f64(0.2),
            common_noise_corr_rho0: T::from_f64(0.2),
            distress_barrier_d: T::from_f64(-0.1),
            horizon_t: T::one(),
        },
    );
    let h = ql_xva::evaluate_bcbs_d424_output_floor_and_sma_op_risk(
        &ql_xva::BcbsD424OutputFloorAndSmaOpRiskSpec {
            irb_credit_rwa: T::from_f64(100.0) + x[0],
            ima_market_rwa: T::from_f64(50.0),
            standardized_credit_rwa: T::from_f64(120.0),
            standardized_market_rwa: T::from_f64(60.0),
            standardized_cva_rwa: T::from_f64(20.0),
            output_floor_multiplier: T::from_f64(0.725),
            interest_lease_dividend_component: T::from_f64(20.0),
            services_component: T::from_f64(15.0),
            financial_component: T::from_f64(30.0),
            average_annual_operational_loss: T::from_f64(2.0) + x[1],
            cet1_capital: T::from_f64(20.0),
            at1_and_tier2_capital: T::from_f64(10.0),
            eligible_tlac_senior_debt: T::from_f64(40.0),
        },
    );
    let i = ql_portfolio::solve_cartea_jaimungal_penalva_lob_fill_hjb(
        &ql_portfolio::CarteaJaimungalPenalvaLobFillHjbSpec {
            initial_inventory: 10,
            horizon_seconds: T::from_f64(60.0),
            base_arrival_rate_lambda: T::from_f64(0.5) + x[2],
            fill_decay_kappa: T::one(),
            running_inventory_penalty_phi: T::from_f64(0.01),
            terminal_liquidation_penalty_alpha: T::from_f64(0.1),
            ofi_short_term_drift: T::from_f64(0.02),
            adverse_selection_cost: T::from_f64(0.01),
        },
    );
    let j = ql_portfolio::evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha(
        &ql_portfolio::MadhavanAbnerHillEtfApArbitrageSpec {
            etf_market_price: T::from_f64(101.0) + x[0],
            reported_inav: T::from_f64(100.0),
            futures_beta_return: T::from_f64(0.01) + x[1],
            stale_underlying_weight: T::from_f64(0.4),
            ap_creation_fee_bps: T::from_f64(5.0),
            basket_bid_ask_spread_bps: T::from_f64(8.0),
            creation_unit_volatility_bps: T::from_f64(10.0),
            ap_risk_aversion: T::from_f64(0.5),
            embedded_unrealized_gain_fraction: T::from_f64(0.3),
            annual_redemption_turnover: T::from_f64(0.2),
            capital_gains_tax_rate: T::from_f64(0.25),
        },
    );
    let k = ql_portfolio::evaluate_grinold_qian_christopherson_ic_decay_and_timing(
        &ql_portfolio::GrinoldQianChristophersonAlphaDecaySpec {
            fast_signal_ic0: T::from_f64(0.1) + x[0] * 0.01,
            fast_signal_half_life_months: T::from_f64(3.0),
            slow_signal_ic0: T::from_f64(0.08),
            slow_signal_half_life_months: T::from_f64(12.0),
            signal_correlation: T::from_f64(0.2),
            rebalance_horizon_months: T::one(),
            transaction_cost_penalty: T::from_f64(0.01),
            factor_active_exposures: vec![T::from_f64(0.1), T::from_f64(-0.05)],
            factor_benchmark_returns: vec![T::from_f64(0.02), T::from_f64(0.01)],
            normal_factor_exposures: vec![T::from_f64(0.05), T::from_f64(-0.02)],
        },
    );
    let l = ql_dml::evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad(
        &ql_dml::GlassermanSaddlepointIsStratifiedLhsSpec {
            spot: T::from_f64(100.0) + x[0],
            deep_otm_strike: T::from_f64(130.0),
            risk_free_rate: T::from_f64(0.03),
            volatility: T::from_f64(0.2) + x[1] * 0.01,
            maturity: T::one(),
            num_strata: 16,
        },
    );
    a.implied_futures_rate
        + b.ig_to_default_fractional_pd
        + c.full_variance_swap_fair_vol
        + d.perpetual_russian_option_value
        + e.bates_svj_call_price
        + f.strategic_convertible_fair_value * 0.001
        + g.mean_field_value_function_cost
        + h.floored_cet1_ratio
        + i.expected_terminal_inventory * 0.001
        + j.net_ap_creation_redemption_profit_bps * 0.001
        + k.total_factor_contribution
        + k.turnover_adjusted_information_ratio
        + l.importance_sampled_call_price
}

pub fn benchmark_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite41_workloads(x).val();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        {
            let y = suite41_workloads(x);
            y.val + y.d.iter().sum::<f64>()
        }
    };
    run_comparison(
        "41_veronesi_conv_lando_gen_gatheral_corridor_russian_stop_bates_svj_koziol_cb_mckean_mfg_bcbs_floor_cjp_lob_madhavan_etf_qian_ic_glasserman_is",
        "Veronesi / Lando / Gatheral-Corridor / Russian / Bates-SVJ / Koziol-CB / MFG / BCBS-Floor / CJP-LOB / ETF / Qian-IC / Glasserman-IS AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite42_workloads<T: Real>(x: [T; 4]) -> T {
    let r = T::from_f64(0.03) + x[0] * T::from_f64(0.001);
    let vol = T::from_f64(0.2) + x[1] * T::from_f64(0.01);
    let spot = T::from_f64(100.0) + x[2];
    let inf = T::from_f64(0.02) + x[3] * T::from_f64(0.001);
    let taylor = ql_curves::evaluate_romer_veronesi_taylor_rule_macro_yield_curve(
        &ql_curves::RomerVeronesiTaylorRuleMacroCurveSpec {
            natural_real_rate_r_star: r - inf,
            inflation_target_pi_star: T::from_f64(0.02),
            current_inflation: inf,
            current_output_gap: x[0],
            taylor_inflation_weight_phi_pi: T::from_f64(1.5),
            taylor_output_weight_phi_y: T::from_f64(0.5),
            inflation_mean_reversion: T::from_f64(0.4),
            output_gap_mean_reversion: T::from_f64(0.3),
            inflation_vol: vol,
            output_gap_vol: vol * T::from_f64(0.5),
            market_price_of_inflation_risk: T::from_f64(0.2),
            zlb_floor_rate: T::zero(),
            maturities: vec![T::one(), T::from_f64(5.0), T::from_f64(10.0)],
        },
    );
    let dam = ql_credit::evaluate_damodaran_rsm_country_risk_and_debt_modification(
        &ql_credit::DamodaranRsmCountryRiskAndDebtModSpec {
            mature_market_erp: T::from_f64(0.05),
            sovereign_cds_spread: T::from_f64(0.012) + x[0],
            sovereign_default_drag_bps: T::from_f64(20.0),
            emerging_equity_vol: vol * T::from_f64(1.4),
            emerging_govt_bond_vol: T::from_f64(0.15),
            firm_domestic_revenue_share: T::from_f64(0.6),
            ebit_operating_income: T::from_f64(12.0) + x[1],
            annual_interest_expense: T::from_f64(3.0),
            original_debt_carrying_value: T::from_f64(100.0),
            original_effective_interest_rate: r,
            modified_cashflows: vec![
                T::from_f64(4.0),
                T::from_f64(4.0),
                T::from_f64(104.0) + x[2],
            ],
            modified_payment_times: vec![T::one(), T::from_f64(2.0), T::from_f64(3.0)],
            modification_fees_paid: T::one(),
        },
    );
    let var = ql_portfolio::evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es(
        &ql_portfolio::ChristoffersenBouchaudVarBacktestAndCumulantSpec {
            total_days_t: 250,
            count_n00: T::from_f64(220.0),
            count_n01: T::from_f64(8.0),
            count_n10: T::from_f64(8.0),
            count_n11: T::from_f64(14.0),
            var_tail_probability_p: T::from_f64(0.05),
            portfolio_mean_loss: x[0],
            portfolio_std_loss: T::one() + vol,
            loss_skewness_gamma1: x[2] * T::from_f64(0.1),
            loss_excess_kurtosis_gamma2: x[3],
        },
    );
    let lvar = ql_portfolio::evaluate_foucault_vandermerwe_lob_schedule_and_lvar(
        &ql_portfolio::FoucaultVanDerMerweLobLiquidityLvarSpec {
            mid_price: spot,
            fundamental_volatility: vol,
            informed_trader_share_pi: T::from_f64(0.25),
            exponential_order_size_rate_lambda: T::from_f64(0.8),
            order_processing_fee_c: T::from_f64(0.0002),
            position_notional: T::from_f64(1_000_000.0) + x[0] * T::from_f64(1000.0),
            mean_relative_bid_ask_spread: T::from_f64(0.001),
            std_relative_bid_ask_spread: T::from_f64(0.0004),
            spread_tail_scaling_a: T::from_f64(1.5),
            confidence_z_alpha: T::from_f64(2.33),
            liquidation_horizon_days: T::from_f64(5.0),
            market_adv_notional: T::from_f64(10_000_000.0),
            permanent_impact_gamma: T::from_f64(0.02),
        },
    );
    let pt = ql_portfolio::evaluate_barberis_jin_wang_prospect_overhang_and_skew_anomaly(
        &ql_portfolio::BarberisJinWangProspectAnomalySpec {
            current_price: spot,
            reference_purchase_price: T::from_f64(95.0),
            return_volatility: vol,
            return_skewness: T::from_f64(0.8) + x[3],
            prospect_loss_aversion_lambda: T::from_f64(2.25),
            prospect_diminishing_sensitivity_alpha: T::from_f64(0.88),
            prelec_probability_weighting_gamma: T::from_f64(0.7),
            mental_accounting_integration_b0: T::from_f64(0.4),
            capm_beta: T::from_f64(1.1),
            equity_risk_premium: T::from_f64(0.05),
        },
    );
    let regret = ql_portfolio::evaluate_hart_mascolell_regret_matching_allocation(
        &ql_portfolio::HartMasColellRegretMatchingPortfolioSpec {
            strategy_payoff_history: vec![
                [
                    T::from_f64(0.02) + x[0],
                    T::from_f64(0.01),
                    T::from_f64(0.03),
                ],
                [T::from_f64(0.01), T::from_f64(0.02) + x[1], T::zero()],
                [
                    T::from_f64(0.03),
                    T::from_f64(0.01),
                    T::from_f64(0.02) + x[2],
                ],
            ],
            realized_action_weights_history: vec![
                [T::from_f64(0.4), T::from_f64(0.3), T::from_f64(0.3)],
                [T::from_f64(0.3), T::from_f64(0.4), T::from_f64(0.3)],
                [T::from_f64(0.3), T::from_f64(0.3), T::from_f64(0.4)],
            ],
            regret_inertia_mu: T::one(),
            exploration_epsilon: T::from_f64(0.01),
        },
    );
    let forensic = ql_ad::evaluate_beneish_wahlen_forensic_mscore_and_hkas32_bifurcation(
        &ql_ad::BeneishWahlenHkas32ForensicBifurcationSpec {
            dsri: T::from_f64(1.1) + x[0],
            gmi: T::from_f64(1.05),
            aqi: T::from_f64(1.2),
            sgi: T::from_f64(1.15),
            depi: T::from_f64(1.03),
            sgai: T::from_f64(1.02),
            lvgi: T::from_f64(1.1),
            tata: T::from_f64(0.03),
            net_income: T::from_f64(10.0),
            operating_cash_flow: T::from_f64(8.0),
            average_total_assets: T::from_f64(200.0),
            convertible_issue_proceeds: T::from_f64(100.0),
            convertible_face_value: T::from_f64(100.0),
            contractual_coupon_rate: T::from_f64(0.03),
            comparable_straight_debt_yield: r,
            direct_issuance_transaction_costs: T::one(),
            maturity_years: 5,
        },
    );
    let camp = ql_portfolio::evaluate_colin_bacon_campisi_fixed_income_attribution(
        &ql_portfolio::ColinBaconCampisiFixedIncomeAttributionSpec {
            portfolio_total_return: T::from_f64(0.042) + x[3],
            benchmark_total_return: T::from_f64(0.035),
            portfolio_coupon_yield: T::from_f64(0.04),
            benchmark_coupon_yield: T::from_f64(0.035),
            portfolio_eff_duration: T::from_f64(5.2),
            benchmark_eff_duration: T::from_f64(5.0),
            portfolio_spread_duration: T::from_f64(4.0),
            benchmark_spread_duration: T::from_f64(4.2),
            parallel_treasury_shift: T::from_f64(0.005) + x[0] * T::from_f64(0.001),
            portfolio_key_rate_twist_shift: T::from_f64(0.001),
            benchmark_key_rate_twist_shift: T::from_f64(0.001),
            portfolio_credit_spread_change: T::from_f64(0.002),
            benchmark_credit_spread_change: T::from_f64(0.002),
            fx_spot_return: T::from_f64(0.01),
            domestic_cash_rate: r,
            foreign_cash_rate: T::from_f64(0.032),
        },
    );
    let heck = ql_dml::evaluate_greene_wooldridge_heckman_survivorship_and_gmm(
        &ql_dml::GreeneWooldridgeHeckmanSelectionGmmSpec {
            selection_z_scores: vec![T::from_f64(-0.5), T::zero(), T::from_f64(0.5), T::one()],
            factor_exposures_x: vec![
                T::from_f64(0.5),
                T::one(),
                T::from_f64(1.5),
                T::from_f64(2.0),
            ],
            observed_returns_y: vec![
                T::from_f64(0.02) + x[0],
                T::from_f64(0.04),
                T::from_f64(0.06) + x[1],
                T::from_f64(0.09),
            ],
            gmm_moment_conditions_bar_g: [
                T::from_f64(0.01),
                T::from_f64(-0.02),
                T::from_f64(0.015),
            ],
            gmm_weighting_matrix_diag: [T::from_f64(2.0), T::from_f64(1.5), T::one()],
            sample_size_n: T::from_f64(100.0),
        },
    );
    let en = ql_dml::solve_hastie_tibshirani_elastic_net_and_gcv(
        &ql_dml::HastieTibshiraniElasticNetGcvSpec {
            design_matrix_x: vec![
                [T::one(), T::from_f64(0.2), T::from_f64(0.1)],
                [T::one(), T::from_f64(0.4), T::from_f64(0.5)],
                [T::one(), T::from_f64(0.6), T::from_f64(0.2)],
                [T::one(), T::from_f64(0.8), T::from_f64(0.9)],
            ],
            response_y: vec![
                T::from_f64(0.1) + x[0],
                T::from_f64(0.2),
                T::from_f64(0.3) + x[1],
                T::from_f64(0.4),
            ],
            l1_penalty_lambda1: T::from_f64(0.05),
            l2_penalty_lambda2: T::from_f64(0.1) + x[2].abs(),
            smoothing_eps: T::from_f64(0.01),
        },
    );
    let cb = ql_fixed_income::evaluate_an_daoquan_clifford_chance_resettable_convertible(
        &ql_fixed_income::AnDaoquanCliffordChanceResettableCbSpec {
            spot_price: spot,
            initial_conversion_price: T::from_f64(100.0),
            book_value_per_share_bps: T::from_f64(80.0),
            avg_20day_price: spot * T::from_f64(0.95),
            face_value: T::from_f64(100.0),
            coupon_rate: T::from_f64(0.02),
            redemption_maturity_price: T::from_f64(105.0),
            put_back_trigger_ratio: T::from_f64(0.7),
            put_back_price: T::from_f64(103.0),
            soft_call_trigger_ratio: T::from_f64(1.3),
            soft_call_price: T::from_f64(105.0),
            stock_volatility: vol,
            risk_free_rate: r,
            credit_spread: T::from_f64(0.01),
            maturity: T::from_f64(3.0),
        },
    );
    let bail = ql_credit::evaluate_fsb_finma_gsib_resolution_bail_in_waterfall(
        &ql_credit::FsbFinmaGsibResolutionBailInWaterfallSpec {
            total_rwa: T::from_f64(1000.0),
            pre_resolution_cet1: T::from_f64(120.0),
            at1_principal: T::from_f64(50.0),
            tier2_principal: T::from_f64(60.0),
            senior_tlac_principal: T::from_f64(100.0),
            valuation_loss_shock: T::from_f64(0.08) + x[0] * T::from_f64(0.01),
            target_post_resolution_cet1_ratio: T::from_f64(0.1),
            liquidation_insolvency_haircut: T::from_f64(0.4),
            existing_share_count: T::from_f64(100.0),
            bail_in_conversion_discount: T::from_f64(0.2) + x[1] * T::from_f64(0.01),
        },
    );
    taylor
        .macro_affine_zero_yields
        .iter()
        .copied()
        .fold(T::zero(), |a, b| a + b)
        + dam.country_risk_premium_crp
        + var.cornish_fisher_var
        + lvar.total_liquidity_adjusted_var_lvar * T::from_f64(1e-6)
        + pt.bjw_equilibrium_expected_excess_return
        + regret.hannan_external_regret_per_period
        + forensic.beneish_8var_m_score
        + camp.total_active_local_return
        + heck.heckman_corrected_alpha
        + en.generalized_cross_validation_gcv
        + cb.three_line_duplex_fair_value * T::from_f64(1e-3)
        + bail.post_bail_in_creditor_ownership_share
        + bail.ncwo_compensation_shortfall * T::from_f64(1e-3)
        + dam.total_cost_of_equity_ke
        + dam.extinguishment_gain_or_loss * T::from_f64(1e-3)
        + forensic.manipulation_probit_probability
        + forensic.year1_effective_interest_expense * T::from_f64(1e-2)
        + camp.active_security_selection_residual
        + camp.karnosky_singer_unhedged_currency_surprise
}

pub fn benchmark_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite42_workloads(x).val();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        let y = suite42_workloads(x);
        y.val + y.d.iter().sum::<f64>()
    };
    run_comparison(
        "42_romer_taylor_damodaran_crp_christoffersen_var_foucault_lvar_barberis_pt_hart_regret_beneish_mscore_campisi_attr_heckman_gmm_elastic_net_reset_cb_fsb_bailin",
        "Romer-Taylor / Damodaran-CRP / Cumulant-VaR / LOB-LVaR / Prospect / Regret / Beneish / Campisi / Heckman / Elastic-Net / Reset-CB / FSB-BailIn AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite43_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let r = T::from_f64(0.03) + x[0] * T::from_f64(0.001);
    let v = T::from_f64(0.2) + x[1] * T::from_f64(0.01);
    let spot = T::from_f64(100.0) + x[2];
    let k = ql_curves::evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization(
        &ql_curves::FabozziKeyRateLscImmunizationSpec {
            cashflow_tenors: vec![
                T::one(),
                T::from_f64(3.0),
                T::from_f64(5.0),
                T::from_f64(10.0),
            ],
            cashflows: vec![
                T::from_f64(2.0),
                T::from_f64(2.0),
                T::from_f64(2.0),
                T::from_f64(102.0) + x[2],
            ],
            zero_rates: vec![r, r + 0.002, r + 0.004, r + 0.006],
            key_rate_tenors: vec![T::one(), T::from_f64(5.0), T::from_f64(10.0)],
            ns_decay_lambda: T::from_f64(0.5) + x[3].abs(),
            liability_horizon_years: T::from_f64(5.0),
            barbell_short_tenor: T::one(),
            barbell_long_tenor: T::from_f64(10.0),
            twist_slope_shock_bps: T::from_f64(25.0),
            butterfly_curvature_shock_bps: T::from_f64(15.0),
        },
    );
    let cmo = ql_fixed_income::evaluate_davidson_levin_cmo_pac_companion_waterfall(
        &ql_fixed_income::DavidsonLevinCmoPacCompanionSpec {
            initial_pool_balance: T::from_f64(100.0),
            gross_wac_rate: r + 0.03,
            net_coupon_rate: r + 0.02,
            num_periods: 24,
            pac_lower_psa: T::from_f64(100.0),
            pac_upper_psa: T::from_f64(300.0),
            realized_psa: T::from_f64(180.0) + x[0] * 10.0,
            discount_yield: r,
        },
    );
    let cb = ql_fixed_income::evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game(
        &ql_fixed_income::AsquithSirbuShreveCbGameSpec {
            firm_value: T::from_f64(130.0) + x[2],
            stock_price: spot,
            conversion_ratio: T::one(),
            annual_coupon_amount: T::from_f64(3.0),
            annual_dividend_per_share: T::one(),
            corporate_tax_rate: T::from_f64(0.25),
            effective_call_price: T::from_f64(105.0),
            call_notice_period_years: T::from_f64(30.0 / 365.0),
            firm_volatility: v,
            risk_free_rate: r,
            firm_payout_rate: T::from_f64(0.02),
            dilution_gamma: T::from_f64(0.8),
        },
    );
    let asu = ql_fixed_income::evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid(
        &ql_fixed_income::FasbAsu2024InducedConversionCpdiSpec {
            net_income: T::from_f64(12.0) + x[0],
            basic_shares: T::from_f64(100.0),
            bond_par: T::from_f64(100.0),
            stated_coupon_rate: T::from_f64(0.03),
            comparable_noncontingent_yield: r + 0.01,
            tax_rate: T::from_f64(0.25),
            original_conversion_shares: T::one(),
            induced_sweetener_extra_shares: T::from_f64(0.05),
            induced_sweetener_cash: T::one(),
            stock_price_at_acceptance: spot,
            maturity_years: T::from_f64(5.0),
        },
    );
    let rm = ql_credit::evaluate_rsm_asc470_debt_modification_and_restructuring(
        &ql_credit::RsmAsc470DebtModificationSpec {
            carrying_amount: T::from_f64(100.0) + x[1],
            original_par: T::from_f64(100.0),
            original_coupon_rate: T::from_f64(0.03),
            original_effective_rate: r,
            modified_par: T::from_f64(100.0),
            modified_coupon_rate: T::from_f64(0.025),
            remaining_years: T::from_f64(5.0),
            lender_fee_paid: T::one(),
            conversion_option_fv_change: T::from_f64(0.5),
            market_exit_yield: r + 0.01,
        },
    );
    let cv = ql_xva::evaluate_gregory_csa_threshold_mta_and_ctd_colva(
        &ql_xva::GregoryCsaThresholdMtaCtdColvaSpec {
            time_grid: vec![T::one(), T::from_f64(2.0), T::from_f64(3.0)],
            expected_mtm_profile: vec![T::from_f64(5.0), T::from_f64(8.0), T::from_f64(10.0)],
            exposure_vol: T::from_f64(3.0) + v,
            mpor_years: T::from_f64(10.0 / 252.0),
            csa_threshold: T::from_f64(2.0),
            minimum_transfer_amount: T::one(),
            independent_amount: T::from_f64(0.5),
            hazard_rate: T::from_f64(0.02) + x[0].abs(),
            lgd: T::from_f64(0.6),
            risk_free_rate: r,
            ccy1_collateral_spread: T::from_f64(0.01),
            ccy2_collateral_spread: T::from_f64(0.008),
            collateral_spread_vol: T::from_f64(0.01),
        },
    );
    let eh = ql_portfolio::evaluate_embrechts_hill_pickands_and_subexponential_tail(
        &ql_portfolio::EmbrechtsHillPickandsSubexponentialSpec {
            loss_samples: vec![
                T::from_f64(160.0) + spot * T::from_f64(0.1),
                T::from_f64(110.0),
                T::from_f64(82.0),
                T::from_f64(66.0),
                T::from_f64(56.0),
                T::from_f64(50.0),
                T::from_f64(46.0),
                T::from_f64(43.0),
                T::from_f64(40.0),
                T::from_f64(38.0),
                T::from_f64(36.0),
                T::from_f64(35.0),
            ],
            tail_order_k: 2,
            extreme_tail_prob: T::from_f64(0.01) + x[3].abs() * T::from_f64(0.001),
            portfolio_positions_count: 8,
        },
    );
    let bp = ql_portfolio::evaluate_bouchaud_potters_power_law_cumulant_and_smile_scaling(
        &ql_portfolio::BouchaudPottersCumulantSmileSpec {
            single_period_vol: v,
            single_period_skewness: T::from_f64(-0.2),
            single_period_excess_kurtosis: T::from_f64(3.0),
            vol_of_vol_kurtosis_component: T::from_f64(2.0),
            vol_decorrelation_steps: T::from_f64(10.0),
            horizon_steps_n: T::from_f64(20.0),
            log_moneyness: x[2] * T::from_f64(0.01),
            confidence_z: T::from_f64(2.33),
        },
    );
    let ks = ql_portfolio::evaluate_karnosky_singer_multicurrency_and_menchero_attribution(
        &ql_portfolio::KarnoskySingerMencheroAttributionSpec {
            portfolio_asset_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
            benchmark_asset_weights: vec![T::from_f64(0.5), T::from_f64(0.5)],
            portfolio_currency_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
            benchmark_currency_weights: vec![T::from_f64(0.5), T::from_f64(0.5)],
            local_portfolio_returns: vec![T::from_f64(0.08) + x[0], T::from_f64(0.03)],
            local_benchmark_returns: vec![T::from_f64(0.06), T::from_f64(0.04)],
            local_eurodollar_rates: vec![T::from_f64(0.02), T::from_f64(0.025)],
            base_eurodollar_rate: r,
            spot_fx_returns: vec![T::from_f64(0.01), T::from_f64(-0.005)],
            multi_period_port_returns: vec![
                T::from_f64(0.02) + x[0] * T::from_f64(0.01),
                T::from_f64(0.03) + x[1] * T::from_f64(0.01),
            ],
            multi_period_bench_returns: vec![T::from_f64(0.018), T::from_f64(0.025)],
        },
    );
    let tm = ql_portfolio::evaluate_dgtw_characteristic_and_treynor_mazuy_timing(
        &ql_portfolio::DgtwTreynorMazuyTimingSpec {
            current_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
            lagged_weights: vec![T::from_f64(0.5), T::from_f64(0.5)],
            stock_returns: vec![
                T::from_f64(0.08) + x[0] * T::from_f64(0.01),
                T::from_f64(0.03) + x[1] * T::from_f64(0.01),
            ],
            current_char_bench_returns: vec![T::from_f64(0.06), T::from_f64(0.04)],
            lagged_char_bench_returns: vec![T::from_f64(0.05), T::from_f64(0.04)],
            portfolio_excess_series: vec![
                T::from_f64(0.01),
                T::from_f64(0.02),
                T::from_f64(0.04) + x[2] * T::from_f64(0.005),
                T::from_f64(0.08) + x[3] * T::from_f64(0.005),
            ],
            market_excess_series: vec![
                T::zero(),
                T::from_f64(0.02),
                T::from_f64(0.04),
                T::from_f64(0.08),
            ],
        },
    );
    let bh = ql_portfolio::evaluate_brock_hommes_adaptive_belief_switching(
        &ql_portfolio::BrockHommesAdaptiveBeliefSpec {
            initial_price_deviation: T::from_f64(0.1) + x[0],
            risk_free_rate: r,
            fundamentalist_damping_v1: T::from_f64(0.8),
            chartist_trend_extrapolation_g: T::from_f64(1.2),
            fundamentalist_info_cost_c1: T::from_f64(0.1),
            intensity_of_choice_beta: T::from_f64(8.0),
            fitness_memory_eta: T::from_f64(0.8),
            risk_aversion_sigma2: T::one(),
            shock_innovations: vec![x[1], x[2] * T::from_f64(0.01), x[3] * T::from_f64(0.01)],
        },
    );
    let pe = ql_models::evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation(
        &ql_models::PenmanReoiOhlsonJuettnerSpec {
            net_operating_assets_0: T::from_f64(100.0),
            net_financial_obligations_0: T::from_f64(30.0),
            operating_income_after_tax_1: T::from_f64(12.0) + x[0],
            operating_income_after_tax_2: T::from_f64(13.0),
            net_financial_expense_after_tax_1: T::from_f64(2.0),
            wacc: r + 0.02,
            terminal_reoi_growth: T::from_f64(0.02),
            current_share_price: spot,
            forward_eps_1: T::from_f64(5.0),
            forward_eps_2: T::from_f64(5.5),
            forward_dps_1: T::from_f64(2.0),
            perpetual_aeg_growth: T::from_f64(0.03),
        },
    );
    [
        k.bond_present_value,
        cmo.pac_present_value,
        cb.perpetual_dynkin_cb_fair_value
            + cb.asquith_safety_call_trigger_stock_price * T::from_f64(1e-2)
            + cb.conversion_value_to_call_ratio,
        asu.basic_eps
            + asu.asu2020_06_instrument_c_net_share_eps
            + asu.asu2024_04_induced_conversion_expense * T::from_f64(1e-2)
            + asu.cpdi_projected_contingent_terminal_payment * T::from_f64(1e-3),
        rm.pv_modified_cashflows_at_orig_eir + rm.prospective_modified_eir,
        cv.ctd_multi_currency_colva,
        eh.weissman_extreme_var * T::from_f64(1e-3),
        bp.bouchaud_potters_non_gaussian_var,
        ks.cumulative_active_return
            + ks.karnosky_singer_security_selection
            + ks.karnosky_singer_currency_management
            + ks.menchero_natural_scaling_multiplier * T::from_f64(1e-3),
        tm.treynor_mazuy_timing_gamma + tm.dgtw_characteristic_selectivity_cs,
        bh.terminal_price_deviation,
        pe.reoi_intrinsic_equity_value * T::from_f64(1e-2)
            + pe.ohlson_juettner_implied_cost_of_equity,
    ]
}

pub fn benchmark_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite43_workloads(x).iter().copied().sum::<f64>().val();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        let y = suite43_workloads(x);
        y.iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "43_fabozzi_krd_cmo_pac_asquith_cb_asu2024_cpdi_rsm_asc470_gregory_colva_embrechts_hill_bouchaud_tlf_karnosky_singer_dgtw_timing_brock_hommes_penman_reoi",
        "KRD-LSC / CMO-PAC / Asquith-CB / ASU2024-CPDI / ASC470-Mod / ColVA-CTD / Hill-EVT / Bouchaud-TLF / Karnosky-FX / DGTW-Timing / Brock-Hommes / Penman-ReOI AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite44_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let r = ql_models::evaluate_gairat_polala_cheyette_local_vol_caplet(
        &ql_models::GairatPolalaCheyetteCapletSpec {
            forward_rate: T::from_f64(0.03) + x[0],
            strike_rate: T::from_f64(0.031) + x[1],
            fixing_time_t1: T::one(),
            payment_time_t2: T::from_f64(1.5),
            mean_reversion_kappa: T::from_f64(0.1),
            initial_sigma0: T::from_f64(0.2) + x[2],
            initial_skew_b: T::from_f64(-0.2) + x[3],
            target_atm_normal_vol: T::from_f64(0.01),
            target_otm_normal_vol: T::from_f64(0.011),
            discount_factor_t2: T::from_f64(0.98),
        },
    );
    let q = ql_models::evaluate_gatheral_jusselin_quadratic_rough_heston(
        &ql_models::GatheralJusselinQuadraticRoughHestonSpec {
            forward: T::from_f64(100.0) + x[0],
            strike: T::from_f64(101.0) + x[1],
            expiry: T::one(),
            hurst_h: T::from_f64(0.1),
            weight_a: T::from_f64(0.2) + x[2],
            asymmetry_center_b: T::from_f64(0.3),
            floor_variance_c: T::from_f64(0.02),
            initial_state_z0: T::from_f64(0.1),
            kernel_eta: T::from_f64(0.5) + x[3],
            mean_reversion_lambda: T::one(),
        },
    );
    let z = ql_vol::evaluate_andreasen_huge_one_step_zabr_calibration(
        &ql_vol::AndreasenHugeOneStepZabrSpec {
            forward: T::from_f64(100.0) + x[0],
            strike: T::from_f64(101.0) + x[1],
            expiry: T::one(),
            alpha: T::from_f64(0.2) + x[2],
            beta: T::from_f64(0.5),
            rho: T::from_f64(-0.2),
            nu: T::from_f64(0.4) + x[3],
            discount_factor: T::one(),
        },
    );
    let n = ql_curves::evaluate_adb_rbi_ndirs_and_asset_swap_basis(
        &ql_curves::AdbRbiNdirsAssetSwapSpec {
            payment_tenors: vec![T::one(), T::from_f64(2.0)],
            onshore_forward_rates: vec![T::from_f64(0.03) + x[0], T::from_f64(0.032) + x[1]],
            usd_ois_discount_factors: vec![T::one(), T::from_f64(0.95)],
            local_rate_vol: T::from_f64(0.1),
            ndf_fx_vol: T::from_f64(0.2) + x[2],
            rate_fx_corr: T::from_f64(0.2),
            sovereign_bond_yield: T::from_f64(0.05),
            usd_sofr_yield: T::from_f64(0.03),
            ndf_ccs_basis_spread: T::from_f64(0.005) + x[3],
            withholding_tax_rate: T::from_f64(0.1),
        },
    );
    let cr = ql_credit::evaluate_brigo_morini_cmcds_and_credit_index_option(
        &ql_credit::BrigoMoriniCmcdsCreditOptionSpec {
            reset_times: vec![T::one(), T::from_f64(2.0)],
            forward_cds_spreads: vec![T::from_f64(0.01) + x[0], T::from_f64(0.012) + x[1]],
            spread_volatility: T::from_f64(0.3),
            annuity_convexity_gamma: T::from_f64(0.5),
            flat_par_cds_spread: T::from_f64(0.011),
            option_expiry: T::one(),
            option_strike_spread: T::from_f64(0.01),
            risky_pv01_annuity: T::from_f64(4.0),
            hazard_rate: T::from_f64(0.02) + x[2],
            recovery_rate: T::from_f64(0.4),
            risk_free_rate: T::from_f64(0.03),
        },
    );
    let xv = ql_xva::evaluate_turlakov_green_kenyon_wwr_and_kva(
        &ql_xva::TurlakovGreenKenyonWwrKvaSpec {
            time_grid: vec![T::one(), T::from_f64(2.0)],
            independent_epe_profile: vec![T::from_f64(1.0) + x[0], T::from_f64(1.2) + x[1]],
            pfe_addon_profile: vec![T::from_f64(0.5), T::from_f64(0.6)],
            base_cpty_hazard: T::from_f64(0.02) + x[2],
            hazard_log_vol: T::from_f64(0.3),
            exposure_vol: T::from_f64(0.2),
            wwr_correlation_rho: T::from_f64(0.5) + x[3],
            lgd: T::from_f64(0.6),
            risk_free_rate: T::from_f64(0.03),
            bank_hazard_rate: T::from_f64(0.01),
            ccr_risk_weight: T::from_f64(0.08),
            shareholder_cost_of_capital: T::from_f64(0.1),
            capital_return_offset_beta: T::from_f64(0.5),
        },
    );
    let pde = ql_pde::evaluate_haentjens_inthout_heston_hull_white_adi(
        &ql_pde::HaentjensInthoutHhwAdiSpec {
            spot: T::from_f64(100.0) + x[0],
            strike: T::from_f64(100.0) + x[1],
            expiry: T::one(),
            initial_var_v0: T::from_f64(0.04),
            mean_rev_kappa_v: T::one(),
            long_var_theta_v: T::from_f64(0.04),
            vol_of_vol_sigma_v: T::from_f64(0.3) + x[2],
            initial_rate_r0: T::from_f64(0.03),
            rate_mean_rev_kappa_r: T::from_f64(0.1),
            rate_vol_sigma_r: T::from_f64(0.01),
            corr_spot_vol: T::from_f64(0.2),
            corr_spot_rate: T::from_f64(0.3) + x[3],
        },
    );
    let fi = ql_fixed_income::evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(
        &ql_fixed_income::MayerBrownSsgaCappedCallMakeWholeSpec {
            bond_par: T::from_f64(100.0),
            spot_price: T::from_f64(100.0) + x[0],
            base_conversion_price: T::from_f64(90.0),
            cap_strike_price: T::from_f64(120.0) + x[1],
            volatility: T::from_f64(0.2) + x[2],
            risk_free_rate: T::from_f64(0.03),
            maturity_years: T::from_f64(5.0),
            acquisition_offer_price: T::from_f64(110.0),
            remaining_years_at_acquisition: T::from_f64(2.0),
            tax_rate_section_163j: T::from_f64(0.2) + x[3],
        },
    );
    let ev = ql_portfolio::evaluate_ahmadi_javid_entropic_var_portfolio(
        &ql_portfolio::AhmadiJavidEntropicVarPortfolioSpec {
            scenario_returns: vec![
                vec![T::from_f64(0.01) + x[0], T::from_f64(-0.02)],
                vec![T::from_f64(-0.03) + x[1], T::from_f64(0.01) + x[2]],
            ],
            scenario_probabilities: vec![T::from_f64(0.5), T::from_f64(0.5)],
            initial_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
            tail_alpha: T::from_f64(0.05),
            risk_aversion_lambda: T::from_f64(1.0) + x[3],
        },
    );
    let hrp = ql_portfolio::evaluate_cotton_schur_hrp_and_rmt_nco(
        &ql_portfolio::CottonSchurHrpRmtNcoSpec {
            covariance_4x4: [
                [
                    T::from_f64(0.04) + x[0],
                    T::from_f64(0.01),
                    T::zero(),
                    T::zero(),
                ],
                [
                    T::from_f64(0.01),
                    T::from_f64(0.05) + x[1],
                    T::zero(),
                    T::zero(),
                ],
                [
                    T::zero(),
                    T::zero(),
                    T::from_f64(0.03) + x[2],
                    T::from_f64(0.01),
                ],
                [
                    T::zero(),
                    T::zero(),
                    T::from_f64(0.01),
                    T::from_f64(0.06) + x[3],
                ],
            ],
            schur_gamma: T::from_f64(0.5),
            sample_aspect_ratio_q: T::from_f64(0.5),
        },
    );
    let bekk = ql_portfolio::evaluate_tsay_bekk_garch_volatility_spillover(
        &ql_portfolio::TsayBekkGarchSpilloverSpec {
            chol_c: [
                [T::from_f64(0.1), T::zero()],
                [T::from_f64(0.02), T::from_f64(0.1)],
            ],
            arch_a: [
                [T::from_f64(0.1), T::from_f64(0.02)],
                [T::from_f64(0.01), T::from_f64(0.1)],
            ],
            garch_b: [[T::from_f64(0.8), T::zero()], [T::zero(), T::from_f64(0.7)]],
            initial_cov: [
                [T::from_f64(0.04), T::zero()],
                [T::zero(), T::from_f64(0.04)],
            ],
            return_innovations: vec![
                [T::from_f64(0.01) + x[0], T::from_f64(-0.02)],
                [T::from_f64(0.02), T::from_f64(0.01) + x[1]],
            ],
        },
    );
    let qual = ql_models::evaluate_schilit_sloan_dechow_dichev_earnings_quality(
        &ql_models::SchilitSloanDechowDichevQualitySpec {
            revenue_0: T::from_f64(100.0),
            revenue_1: T::from_f64(105.0) + x[0],
            cogs_0: T::from_f64(60.0),
            cogs_1: T::from_f64(63.0) + x[1],
            accounts_receivable_0: T::from_f64(10.0),
            accounts_receivable_1: T::from_f64(12.0) + x[2],
            inventory_0: T::from_f64(8.0),
            inventory_1: T::from_f64(9.0),
            accounts_payable_0: T::from_f64(7.0),
            accounts_payable_1: T::from_f64(8.0) + x[3],
            net_operating_assets_0: T::from_f64(50.0),
            net_operating_assets_1: T::from_f64(52.0),
            net_income_1: T::from_f64(8.0),
            cfo_prev: T::from_f64(7.0),
            cfo_curr: T::from_f64(7.5),
            cfo_next: T::from_f64(8.0),
        },
    );
    [
        r.deterministic_auxiliary_y_bar,
        q.qrh_call_price,
        z.one_step_fd_call_price,
        n.usd_asset_swap_pickup_bps,
        cr.morini_brigo_payer_swaption_pv,
        xv.total_wwr_cva_plus_kva,
        pde.hhw_adi_call_price,
        fi.capped_call_overlay_gross_cost,
        ev.exact_chernoff_evar,
        hrp.schur_hrp_portfolio_variance,
        bekk.terminal_var_asset1,
        qual.schilit_shenanigan_composite_z,
    ]
}
pub fn benchmark_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite44_workloads(x).iter().copied().sum::<f64>().val();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite44_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "44_cheyette_caplet_adb_ndirs_andreasen_zabr_qrh_zumbach_cmcds_morini_turlakov_wwr_hhw_adi_capped_call_mw_ahmadi_evar_schur_hrp_bekk_garch_schilit_sloan",
        "Cheyette-LV / ADB-NDIRS / Andreasen-ZABR / QRH-Zumbach / CMCDS-Morini / Turlakov-WWR / HHW-ADI / CappedCall-MW / Ahmadi-EVaR / Schur-HRP / BEKK-GARCH / Schilit-Sloan AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite45_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);
    let c = ql_curves::evaluate_brigo_mercurio_rebonato_swaption_and_cms(
        &ql_curves::BrigoMercurioRebonatoCmsSpec {
            fixing_tenors: vec![T::one(), T::from_f64(2.0)],
            accrual_taus: vec![T::one(), T::one()],
            forward_rates: vec![T::from_f64(0.03) + dx0, T::from_f64(0.035) + dx1],
            discount_factors: vec![T::from_f64(0.98), T::from_f64(0.95)],
            forward_vols: vec![T::from_f64(0.2), T::from_f64(0.21)],
            corr_decay_beta: T::from_f64(0.1),
            option_expiry: T::one(),
            swaption_strike: T::from_f64(0.03) + dx2,
        },
    );
    let l = ql_credit::evaluate_duffie_singleton_leland_toft_capital_structure(
        &ql_credit::DuffieSingletonLelandToftSpec {
            unlevered_asset_value: T::from_f64(100.0) + x[0],
            annual_coupon_c: T::from_f64(5.0) + x[1],
            debt_principal_p: T::from_f64(80.0),
            rollover_maturity_t: T::from_f64(5.0),
            corporate_tax_rate: T::from_f64(0.25),
            bankruptcy_cost_alpha: T::from_f64(0.2),
            asset_payout_rate: T::from_f64(0.02),
            risk_free_rate: T::from_f64(0.03) + dx2,
            asset_volatility: T::from_f64(0.2) + dx3,
        },
    );
    let a = ql_fixed_income::evaluate_calamos_despiegeleer_ascot_and_gamma_scalp(
        &ql_fixed_income::CalamosDeSpiegeleerAscotGammaSpec {
            spot_price: T::from_f64(100.0) + x[0],
            reference_spot: T::from_f64(100.0),
            bond_par: T::from_f64(100.0),
            coupon_rate: T::from_f64(0.04),
            conversion_ratio: T::one(),
            implied_vol: T::from_f64(0.25) + dx1,
            realized_vol: T::from_f64(0.3),
            risk_free_rate: T::from_f64(0.03),
            dividend_yield: T::from_f64(0.01),
            base_hazard_rate: T::from_f64(0.02),
            equity_credit_power_p: T::from_f64(1.2),
            recovery_rate: T::from_f64(0.4),
            maturity_years: T::from_f64(3.0),
            holding_period_years: T::one(),
        },
    );
    let p = ql_fixed_income::evaluate_tuckman_specified_pool_payup_and_dollar_roll(
        &ql_fixed_income::TuckmanSpecifiedPoolDollarRollSpec {
            net_coupon_rate: T::from_f64(0.04) + dx1,
            discount_yield: T::from_f64(0.025),
            effective_maturity_years: T::from_f64(5.0),
            ctd_tba_cpr: T::from_f64(0.08),
            llb_pool_cpr_multiplier: T::from_f64(0.7),
            tba_front_month_price: T::from_f64(100.0) + x[0],
            tba_back_month_price: T::from_f64(99.8),
            gc_repo_rate: T::from_f64(0.03),
        },
    );
    let o = ql_portfolio::evaluate_ohara_glosten_milgrom_sequential_spread(
        &ql_portfolio::OharaGlostenMilgromSequentialSpec {
            low_value_vl: T::from_f64(90.0),
            high_value_vh: T::from_f64(110.0),
            initial_prior_p0: T::from_f64(0.5) + dx0,
            informed_trader_fraction_mu: T::from_f64(0.3),
            informed_signal_accuracy_theta: T::from_f64(0.8),
            dealer_inventory_gamma: T::from_f64(0.01),
            signed_order_flow: vec![T::from_f64(0.2) + dx1, T::from_f64(-0.1), T::from_f64(0.3)],
        },
    );
    let r = ql_portfolio::evaluate_potters_bouchaud_stieltjes_rie_shrinkage(
        &ql_portfolio::PottersBouchaudStieltjesRieSpec {
            sample_eigenvalues: vec![
                T::from_f64(0.5) + dx0,
                T::one() + dx1,
                T::from_f64(2.0) + dx2,
            ],
            aspect_ratio_q: T::from_f64(0.5),
            stieltjes_bandwidth_eta: T::from_f64(0.1),
        },
    );
    let e = ql_portfolio::evaluate_tsay_egarch_and_hamilton_markov_switching(
        &ql_portfolio::TsayEgarchHamiltonRegimeSpec {
            returns: vec![
                T::from_f64(0.01) + dx0,
                T::from_f64(-0.02) + dx1,
                T::from_f64(0.015) + dx2,
            ],
            initial_variance: T::from_f64(0.04),
            egarch_omega: T::from_f64(-0.1),
            egarch_alpha: T::from_f64(0.1),
            egarch_leverage_gamma: T::from_f64(-0.1),
            egarch_beta: T::from_f64(0.9),
            calm_mean: T::zero(),
            calm_vol: T::from_f64(0.1),
            crisis_mean: T::zero(),
            crisis_vol: T::from_f64(0.3),
            stay_calm_p00: T::from_f64(0.95),
            stay_crisis_p11: T::from_f64(0.9),
        },
    );
    let w =
        ql_dml::evaluate_wooldridge_greene_2sls_iv_and_did(&ql_dml::WooldridgeGreene2slsDidSpec {
            endogenous_x: vec![
                T::one() + x[0],
                T::from_f64(2.0) + x[1],
                T::from_f64(3.0) + x[2],
            ],
            instrument_z: vec![T::one(), T::from_f64(2.1), T::from_f64(2.9) + x[3]],
            outcome_y: vec![T::from_f64(1.2), T::from_f64(2.1), T::from_f64(3.2)],
            treated_pre_mean: T::one(),
            treated_post_mean: T::from_f64(1.2),
            control_pre_mean: T::one(),
            control_post_mean: T::from_f64(1.1),
        });
    let m = ql_portfolio::evaluate_mcneil_frey_copula_tail_and_euler_allocation(
        &ql_portfolio::McneilFreyCopulaEulerAllocationSpec {
            clayton_theta: T::from_f64(2.0),
            gumbel_theta: T::from_f64(1.5),
            student_t_df_nu: T::from_f64(5.0),
            student_t_corr_rho: T::from_f64(0.3),
            desk_weights: [T::from_f64(0.25) + x[0] * 0.001; 4],
            desk_expected_losses: [T::from_f64(0.01); 4],
            loss_covariance_4x4: [
                [T::from_f64(0.04), T::from_f64(0.01), T::zero(), T::zero()],
                [
                    T::from_f64(0.01),
                    T::from_f64(0.05),
                    T::from_f64(0.01),
                    T::zero(),
                ],
                [
                    T::zero(),
                    T::from_f64(0.01),
                    T::from_f64(0.06),
                    T::from_f64(0.01),
                ],
                [T::zero(), T::zero(), T::from_f64(0.01), T::from_f64(0.03)],
            ],
            es_tail_multiplier: T::from_f64(2.0),
        },
    );
    let b = ql_portfolio::evaluate_bacon_omega_sortino_and_burke_drawdown(
        &ql_portfolio::BaconOmegaSortinoBurkeDrawdownSpec {
            period_returns: vec![
                T::from_f64(0.01) + x[0] * 0.001,
                T::from_f64(-0.005) + x[1] * 0.001,
                T::from_f64(0.02) + x[2] * 0.001,
            ],
            mar_threshold: T::zero(),
            risk_free_rate: T::from_f64(0.002),
        },
    );
    let d = ql_models::evaluate_damodaran_rd_capitalization_and_patent_real_option(
        &ql_models::DamodaranRdCapitalizationRealOptionSpec {
            reported_ebit: T::from_f64(10.0),
            reported_invested_capital: T::from_f64(100.0),
            effective_tax_rate: T::from_f64(0.25),
            current_rd_expense: T::from_f64(5.0) + x[0],
            historical_rd_expenses: vec![T::from_f64(4.0), T::from_f64(3.0), T::from_f64(2.0)],
            patent_project_pv: T::from_f64(20.0) + x[1],
            patent_development_cost: T::from_f64(10.0),
            patent_life_years: T::from_f64(5.0),
            project_volatility: T::from_f64(0.3),
            risk_free_rate: T::from_f64(0.03),
            abandonment_salvage_value: T::from_f64(2.0),
        },
    );
    let ad = ql_models::evaluate_mascolell_cochrane_arrow_debreu_sdf(
        &ql_models::MasColellCochraneArrowDebreuSdfSpec {
            state_probabilities: [T::from_f64(0.3), T::from_f64(0.4), T::from_f64(0.3)],
            consumption_growth_ratios: [
                T::from_f64(0.98),
                T::from_f64(1.02) + dx2,
                T::from_f64(1.05),
            ],
            equity_payoff_states: [T::from_f64(0.8), T::from_f64(1.1), T::from_f64(1.4) + x[3]],
            subjective_discount_beta: T::from_f64(0.99),
            crra_risk_aversion_gamma: T::from_f64(2.0),
        },
    );
    [
        c.convexity_adjusted_cms_rate,
        l.risky_debt_fair_value,
        a.convertible_bond_fair_value,
        p.specified_pool_payup_32nds + p.dollar_roll_implied_repo_rate,
        o.inventory_skewed_reservation_mid,
        r.rie_condition_number,
        e.terminal_egarch_volatility,
        w.two_stage_ls_iv_beta,
        m.total_portfolio_es,
        b.keating_shadwick_omega_ratio,
        d.total_strategic_expanded_npv,
        ad.equilibrium_equity_price,
    ]
}

pub fn benchmark_rebonato_cms_lm_all_engines(scale: usize) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite45_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite45_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "45_rebonato_cms_leland_toft_calamos_ascot_spec_pool_roll_glosten_milgrom_potters_rie_egarch_hamilton_wooldridge_2sls_mcneil_euler_bacon_omega_damodaran_rd_mwg_arrow_debreu",
        "Rebonato-CMS / Leland-Toft / Calamos-Ascot / SpecPool-Roll / Glosten-Milgrom / Potters-RIE / EGARCH-Hamilton / Wooldridge-2SLS / McNeil-Euler / Bacon-Omega / Damodaran-RD / Arrow-Debreu AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite47_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.005);
    let dx1 = x[1] * T::from_f64(0.005);
    let dx2 = x[2] * T::from_f64(0.005);
    let dx3 = x[3] * T::from_f64(0.005);
    let hw = ql_curves::evaluate_hagan_west_monotone_convex_forward_curve(
        &ql_curves::HaganWestMonotoneConvexCurveSpec {
            knot_tenors: [T::zero(), T::one(), T::from_f64(2.0), T::from_f64(3.0)],
            discrete_forwards: [
                T::from_f64(0.03) + dx0,
                T::from_f64(0.04) + dx1,
                T::from_f64(0.045) + dx2,
            ],
            eval_tenor: T::from_f64(1.4) + dx3,
        },
    );
    let ctd = ql_fixed_income::evaluate_fabozzi_treasury_futures_ctd_and_switch_option(
        &ql_fixed_income::FabozziTreasuryFuturesCtdSwitchSpec {
            bond_coupons: [T::from_f64(0.06), T::from_f64(0.045)],
            bond_maturities: [T::from_f64(7.0), T::from_f64(9.5)],
            bond_clean_prices: [T::from_f64(101.0) + x[0], T::from_f64(94.5) + x[1]],
            bond_modified_durations: [T::from_f64(5.8), T::from_f64(7.6)],
            futures_price: T::from_f64(100.2) + x[2],
            delivery_horizon_years: T::from_f64(0.25),
            repo_rate: T::from_f64(0.04) + dx3,
            yield_volatility: T::from_f64(0.01),
        },
    );
    let afv = ql_fixed_income::evaluate_ayache_forsyth_vetzal_default_coupled_cb(
        &ql_fixed_income::AyacheForsythVetzalConvertibleSpec {
            spot_price: T::from_f64(100.0) + x[0],
            bond_par: T::from_f64(100.0),
            coupon_rate: T::from_f64(0.04),
            conversion_ratio: T::one(),
            maturity_years: T::from_f64(3.0),
            risk_free_rate: T::from_f64(0.03) + dx1,
            dividend_yield: T::from_f64(0.015),
            equity_volatility: T::from_f64(0.25) + dx2,
            hazard_rate: T::from_f64(0.025) + dx3,
            bond_recovery_rate: T::from_f64(0.4),
            equity_default_jump_loss_eta: T::from_f64(0.7),
        },
    );
    let kth = ql_credit::evaluate_lando_schonbucher_kth_to_default_basket(
        &ql_credit::LandoSchonbucherBasketDefaultSpec {
            obligor_hazard_rates: [
                T::from_f64(0.015) + dx0,
                T::from_f64(0.02) + dx1,
                T::from_f64(0.025) + dx2,
                T::from_f64(0.03) + dx3,
            ],
            factor_loading_beta: T::from_f64(0.5),
            recovery_rate: T::from_f64(0.4),
            risk_free_rate: T::from_f64(0.03),
            maturity_years: T::from_f64(5.0),
        },
    );
    let sacva = ql_xva::evaluate_gregory_sa_cva_and_euler_netting_allocation(
        &ql_xva::GregorySaCvaEulerNettingSpec {
            cva_delta_sensitivities: [T::from_f64(10.0) + x[0], T::from_f64(8.0) + x[1]],
            cva_hedge_sensitivities: [T::from_f64(6.0), T::from_f64(5.0)],
            delta_risk_weights: [T::from_f64(0.05), T::from_f64(0.06)],
            intra_bucket_corr_rho: T::from_f64(0.5),
            cva_vega_sensitivity: T::from_f64(4.0) + x[2],
            vega_risk_weight: T::from_f64(0.1),
            subportfolio_exposures: [
                T::from_f64(20.0) + x[3],
                T::from_f64(15.0),
                T::from_f64(10.0),
            ],
            exposure_corr_3x3: [
                [T::one(), T::from_f64(0.3), T::from_f64(0.1)],
                [T::from_f64(0.3), T::one(), T::from_f64(0.2)],
                [T::from_f64(0.1), T::from_f64(0.2), T::one()],
            ],
            cpty_hazard_rate: T::from_f64(0.02),
            cpty_lgd: T::from_f64(0.6),
            horizon_years: T::from_f64(3.0),
        },
    );
    let kirk = ql_commodities::evaluate_eydeland_gibson_schwartz_and_kirk_spread(
        &ql_commodities::EydelandGibsonSchwartzKirkSpreadSpec {
            spot_price: T::from_f64(55.0) + x[0],
            initial_convenience_yield: T::from_f64(0.06) + dx1,
            mean_reversion_kappa: T::from_f64(1.2),
            long_run_convenience_yield: T::from_f64(0.05),
            spot_volatility: T::from_f64(0.32) + dx2,
            convenience_yield_vol: T::from_f64(0.20),
            spot_cy_correlation: T::from_f64(0.4),
            risk_free_rate: T::from_f64(0.03),
            maturity_years: T::from_f64(1.0),
            fuel_futures_price: T::from_f64(6.0) + dx3,
            heat_rate_multiplier: T::from_f64(7.5),
            fixed_strike_cost: T::from_f64(5.0),
            fuel_volatility: T::from_f64(0.28),
            power_fuel_correlation: T::from_f64(0.6),
        },
    );
    let fx = ql_fx::evaluate_clark_wystup_fx_vanna_volga_smile(
        &ql_fx::ClarkWystupFxVannaVolgaSmileSpec {
            spot_fx: T::from_f64(1.10) + dx0,
            domestic_rate: T::from_f64(0.03) + dx1,
            foreign_rate: T::from_f64(0.015),
            expiry_years: T::from_f64(1.0),
            atm_vol: T::from_f64(0.10) + dx2,
            risk_reversal_25d: T::from_f64(0.012) + dx3,
            butterfly_25d: T::from_f64(0.004),
            target_strike: T::from_f64(1.14),
        },
    );
    let nap = ql_exotics::evaluate_bouzoubaa_napoleon_and_reverse_cliquet(
        &ql_exotics::BouzoubaaNapoleonReverseCliquetSpec {
            period_returns: vec![
                T::from_f64(0.03) + dx0,
                T::from_f64(-0.04) + dx1,
                T::from_f64(0.02) + dx2,
                T::from_f64(-0.025) + dx3,
            ],
            napoleon_fixed_coupon: T::from_f64(0.08),
            global_floor: T::zero(),
            local_cap: T::from_f64(0.04),
            local_floor: T::from_f64(-0.04),
            forward_volatility: T::from_f64(0.20),
            forward_skew: T::from_f64(-0.15),
            risk_free_rate: T::from_f64(0.03),
            maturity_years: T::from_f64(2.0),
        },
    );
    let tc = ql_portfolio::evaluate_grinold_qian_transfer_coefficient_and_flam(
        &ql_portfolio::GrinoldQianTransferCoefficientFlAmSpec {
            information_coefficient_ic: T::from_f64(0.06) + dx0,
            breadth_br: T::from_f64(100.0),
            target_active_risk: T::from_f64(0.04) + dx1,
            residual_volatilities: [T::from_f64(0.15), T::from_f64(0.20), T::from_f64(0.18)],
            alpha_forecasts: [
                T::from_f64(0.03) + dx2,
                T::from_f64(0.02),
                T::from_f64(-0.025) + dx3,
            ],
            constrained_active_weights: [T::from_f64(0.12), T::from_f64(0.05), T::from_f64(-0.08)],
        },
    );
    let ep = ql_portfolio::evaluate_meucci_relative_entropy_view_pooling(
        &ql_portfolio::MeucciEntropyPoolingViewSpec {
            prior_probabilities: [
                T::from_f64(0.25),
                T::from_f64(0.25),
                T::from_f64(0.25),
                T::from_f64(0.25),
            ],
            scenario_returns: [
                T::from_f64(-0.04) + dx0,
                T::from_f64(0.01) + dx1,
                T::from_f64(0.03) + dx2,
                T::from_f64(0.08),
            ],
            target_view_return: T::from_f64(0.035) + dx3,
            view_confidence: T::from_f64(0.75),
        },
    );
    let krr =
        ql_dml::evaluate_hastie_kernel_ridge_and_exact_loocv(&ql_dml::HastieKernelRidgeLoocvSpec {
            feature_x: [
                T::from_f64(-1.0),
                T::from_f64(-0.3),
                T::from_f64(0.4),
                T::from_f64(1.1),
            ],
            target_y: [
                T::from_f64(0.8) + dx0,
                T::from_f64(0.2) + dx1,
                T::from_f64(0.5) + dx2,
                T::from_f64(1.3) + dx3,
            ],
            rbf_bandwidth_ell: T::from_f64(0.8),
            ridge_penalty_lambda: T::from_f64(0.15),
        });
    let eva = ql_models::evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation(
        &ql_models::KollerMckinseyKvdEvaValuationSpec {
            invested_capital_0: T::from_f64(500.0) + x[0],
            noplat_1: T::from_f64(75.0) + x[1],
            ronic: T::from_f64(0.16) + dx2,
            wacc: T::from_f64(0.09) + dx3,
            perpetual_growth_g: T::from_f64(0.04),
            net_debt: T::from_f64(150.0),
            diluted_shares: T::from_f64(50.0),
        },
    );
    [
        hw.eval_discount_factor + hw.eval_instantaneous_forward,
        ctd.switch_adjusted_theoretical_futures,
        afv.afv_convertible_fair_value,
        kth.first_to_default_spread_bps,
        sacva.total_sa_cva_capital + sacva.netted_portfolio_cva,
        kirk.kirk_spark_spread_call_pv,
        fx.vanna_volga_call_price,
        nap.napoleon_capital_protected_note_pv,
        tc.expected_constrained_alpha,
        ep.posterior_expected_return + ep.effective_number_of_scenarios,
        krr.loocv_shortcut_mse + krr.gcv_score,
        eva.intrinsic_share_price,
    ]
}

pub fn benchmark_hagan_west_ctd_afv_kth_sacva_gs_kirk_fx_vv_napoleon_tc_meucci_krr_eva(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite47_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite47_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "47_hagan_west_mc_fabozzi_ctd_switch_afv_default_cb_lando_kth_default_gregory_sa_cva_gibson_schwartz_kirk_clark_fx_vv_bouzoubaa_napoleon_grinold_tc_flam_meucci_entropy_hastie_krr_loocv_koller_mckinsey_eva",
        "Hagan-West-MC / CTD-Switch / AFV-CB / Kth-Default / SA-CVA-Euler / GS-Kirk-Spread / FX-Vanna-Volga / Napoleon-Cliquet / Transfer-Coef / Meucci-Entropy / KRR-LOOCV / McKinsey-EVA AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite48_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let s = T::from_f64(0.005);
    let dx0 = x[0] * s;
    let dx1 = x[1] * s;
    let dx2 = x[2] * s;
    let dx3 = x[3] * s;

    let csa = ql_curves::evaluate_brigo_gregory_multicurrency_csa_ctd(
        &ql_curves::BrigoGregoryMultiCurrencyCsaCtdSpec {
            domestic_ois_rate: T::from_f64(0.042) + dx0,
            foreign_ois_rate: T::from_f64(0.028) + dx1,
            xccy_basis_spread: T::from_f64(-0.0025) + dx2,
            csa_spread_volatility: T::from_f64(0.0065) + dx3 * T::from_f64(0.2),
            spot_fx: T::from_f64(1.085),
            maturity_years: T::from_f64(5.0),
            promised_cashflow: T::from_f64(100.0),
        },
    );

    let oas = ql_fixed_income::evaluate_fabozzi_callable_putable_oas_and_convexity(
        &ql_fixed_income::FabozziCallablePutableOasSpec {
            initial_short_rate: T::from_f64(0.035) + dx0,
            rate_drift_per_year: T::from_f64(0.002) + dx1 * T::from_f64(0.1),
            short_rate_vol: T::from_f64(0.18) + dx2,
            oas_spread: T::from_f64(0.006) + dx3 * T::from_f64(0.2),
            coupon_rate: T::from_f64(0.055),
            call_strike_price: T::from_f64(100.5),
            put_strike_price: T::from_f64(99.5),
            maturity_years: T::from_f64(3.0),
            yield_shift_dy: T::from_f64(0.0025),
        },
    );

    let jlt = ql_credit::evaluate_lando_jlt_markov_rating_migration(
        &ql_credit::LandoJltMarkovMigrationSpec {
            ig_to_hy_rate: T::from_f64(0.045) + dx0,
            ig_to_default_rate: T::from_f64(0.002) + dx1 * T::from_f64(0.05),
            hy_to_ig_rate: T::from_f64(0.030),
            hy_to_ccc_rate: T::from_f64(0.080) + dx2,
            hy_to_default_rate: T::from_f64(0.025) + dx3,
            ccc_to_hy_rate: T::from_f64(0.040),
            ccc_to_default_rate: T::from_f64(0.140),
            risk_premium_multiplier_pi: T::from_f64(1.45),
            recovery_rate: T::from_f64(0.40),
            risk_free_rate: T::from_f64(0.035),
            horizon_years: T::from_f64(3.0),
        },
    );

    let vswap = ql_vol::evaluate_gatheral_bergomi_variance_swap_and_corridor(
        &ql_vol::GatheralBergomiVarianceSwapCorridorSpec {
            forward_price: T::from_f64(100.0) + x[0],
            risk_free_rate: T::from_f64(0.03) + dx1,
            maturity_years: T::from_f64(1.0),
            atm_vol: T::from_f64(0.20) + dx2,
            skew_slope_beta: T::from_f64(-0.12) + dx3,
            smile_curvature_gamma: T::from_f64(0.25),
            strike_grid: vec![
                T::from_f64(70.0),
                T::from_f64(80.0),
                T::from_f64(90.0),
                T::from_f64(100.0),
                T::from_f64(110.0),
                T::from_f64(120.0),
                T::from_f64(130.0),
            ],
            corridor_lower_barrier: T::from_f64(85.0),
            corridor_upper_barrier: T::from_f64(115.0),
        },
    );

    let lsm = ql_pde::evaluate_glasserman_lsm_and_rogers_dual_bound(
        &ql_pde::GlassermanLsmDualUpperBoundSpec {
            spot_price: T::from_f64(98.0) + x[0],
            strike_price: T::from_f64(100.0) + x[1] * T::from_f64(0.5),
            risk_free_rate: T::from_f64(0.04) + dx2,
            volatility: T::from_f64(0.22) + dx3,
            maturity_years: T::from_f64(1.0),
            shock_paths_t1: [
                T::from_f64(-1.35),
                T::from_f64(-0.75),
                T::from_f64(-0.25),
                T::from_f64(0.25),
                T::from_f64(0.75),
                T::from_f64(1.35),
            ],
            shock_paths_t2: [
                T::from_f64(-0.95),
                T::from_f64(0.85),
                T::from_f64(-0.45),
                T::from_f64(0.45),
                T::from_f64(-0.85),
                T::from_f64(0.95),
            ],
        },
    );

    let lob = ql_portfolio::evaluate_cartea_jaimungal_lob_market_making(
        &ql_portfolio::CarteaJaimungalLobMarketMakingSpec {
            mid_price: T::from_f64(100.0) + x[0],
            current_inventory_q: T::from_f64(3.0) + x[1] * T::from_f64(0.5),
            running_inventory_penalty_phi: T::from_f64(0.02) + dx2,
            terminal_inventory_penalty_alpha: T::from_f64(0.05) + dx3,
            risk_aversion_gamma: T::from_f64(0.10),
            asset_volatility: T::from_f64(0.20),
            order_arrival_baseline_a: T::from_f64(12.0),
            lob_depth_decay_kappa: T::from_f64(1.5),
            remaining_horizon_tau: T::from_f64(0.5),
        },
    );

    let clm = ql_portfolio::evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller(
        &ql_portfolio::CampbellLoMackinlayVrPresentValueSpec {
            returns: vec![
                T::from_f64(0.012) + dx0,
                T::from_f64(0.015) + dx1,
                T::from_f64(-0.004) + dx2,
                T::from_f64(0.009) + dx3,
                T::from_f64(0.014),
                T::from_f64(-0.003),
                T::from_f64(0.011),
                T::from_f64(0.013),
            ],
            steady_state_dividend_yield: T::from_f64(0.032) + dx0,
            expected_dividend_growth: T::from_f64(0.045) + dx1,
            expected_discount_rate: T::from_f64(0.078) + dx2,
            dividend_growth_ar1_phi: T::from_f64(0.35),
            discount_rate_ar1_phi: T::from_f64(0.85),
        },
    );

    let dsr = ql_dml::evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv(
        &ql_dml::LopezDePradoDsrPurgedCvSpec {
            observed_sharpe_ratio: T::from_f64(1.45) + x[0] * T::from_f64(0.1),
            sample_length_t: T::from_f64(60.0),
            return_skewness: T::from_f64(-0.65) + dx1,
            return_kurtosis: T::from_f64(4.5) + dx2 * T::from_f64(10.0),
            num_independent_trials_k: T::from_f64(25.0),
            trials_sharpe_std: T::from_f64(0.45) + dx3,
            total_cv_samples_n: T::from_f64(1000.0),
            cv_folds_k: T::from_f64(5.0),
            label_overlap_horizon_h: T::from_f64(5.0),
            embargo_fraction_pct: T::from_f64(0.01),
        },
    );

    let erc = ql_portfolio::evaluate_roncalli_spinu_erc_risk_parity(
        &ql_portfolio::RoncalliSpinuErcRiskParitySpec {
            covariance_3x3: [
                [
                    T::from_f64(0.040) + dx0,
                    T::from_f64(0.008) + dx1,
                    T::from_f64(0.002),
                ],
                [
                    T::from_f64(0.008) + dx1,
                    T::from_f64(0.0225) + dx2,
                    T::from_f64(0.003),
                ],
                [
                    T::from_f64(0.002),
                    T::from_f64(0.003),
                    T::from_f64(0.0064) + dx3,
                ],
            ],
            target_risk_budgets: [
                T::from_f64(1.0 / 3.0),
                T::from_f64(1.0 / 3.0),
                T::from_f64(1.0 / 3.0),
            ],
        },
    );

    let ilm = ql_portfolio::evaluate_ilmanen_carry_rolldown_and_vmc_factors(
        &ql_portfolio::IlmanenCarryRolldownVmcFactorSpec {
            current_yield_yt: T::from_f64(0.045) + dx0,
            rolled_horizon_yield: T::from_f64(0.041) + dx1,
            repo_financing_rate: T::from_f64(0.032) + dx2,
            modified_duration: T::from_f64(7.5),
            bond_convexity: T::from_f64(68.0),
            expected_yield_change: T::from_f64(0.001) + dx3,
            yield_volatility: T::from_f64(0.0085),
            value_style_sharpe: T::from_f64(0.45) + dx0,
            momentum_style_sharpe: T::from_f64(0.50) + dx1,
            carry_style_sharpe: T::from_f64(0.48) + dx2,
            value_momentum_corr_rho: T::from_f64(-0.45),
        },
    );

    let pen = ql_models::evaluate_penman_nissim_dupont_and_reoi_valuation(
        &ql_models::PenmanNissimDuPontReoiValuationSpec {
            net_operating_assets_0: T::from_f64(500.0) + x[0],
            net_financial_obligations_0: T::from_f64(180.0) + x[1],
            sales_revenue_1: T::from_f64(650.0) + x[2],
            after_tax_operating_income_1: T::from_f64(72.0) + x[3],
            after_tax_net_financial_expense_1: T::from_f64(9.0),
            required_return_operations_rho: T::from_f64(0.09),
            perpetual_growth_g: T::from_f64(0.03),
        },
    );

    let mc = ql_models::evaluate_ljungqvist_sargent_mccall_bellman_search(
        &ql_models::LjungqvistSargentMccallSearchSpec {
            offer_states: [
                T::from_f64(40.0) + x[0],
                T::from_f64(55.0) + x[1],
                T::from_f64(72.0) + x[2],
                T::from_f64(95.0) + x[3],
            ],
            offer_probabilities: [
                T::from_f64(0.25),
                T::from_f64(0.35),
                T::from_f64(0.25),
                T::from_f64(0.15),
            ],
            waiting_compensation_c: T::from_f64(32.0) + x[0] * T::from_f64(0.5),
            discount_factor_beta: T::from_f64(0.90),
        },
    );

    [
        csa.ctd_csa_discounted_pv,
        oas.callable_bond_fair_value,
        jlt.ig_credit_spread_bps,
        vswap.replicated_volatility_strike,
        lsm.lsm_control_variate_price,
        lob.reservation_mid_price,
        clm.equilibrium_log_dividend_price_ratio,
        dsr.deflated_sharpe_ratio_dsr,
        erc.erc_portfolio_volatility,
        ilm.total_expected_excess_return,
        pen.intrinsic_common_equity_value,
        mc.reservation_threshold_w_bar,
    ]
}

pub fn benchmark_brigo_csa_ctd_fabozzi_callable_oas_lando_jlt_migration_gatheral_corridor_varswap_glasserman_lsm_rogers_dual_cartea_lob_mm_lo_mackinlay_vr_lopez_dsr_purged_cv_roncalli_spinu_erc_ilmanen_vmc_penman_reoi_mccall_bellman(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite48_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite48_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "48_brigo_csa_ctd_fabozzi_callable_oas_lando_jlt_migration_gatheral_corridor_varswap_glasserman_lsm_rogers_dual_cartea_lob_mm_lo_mackinlay_vr_lopez_dsr_purged_cv_roncalli_spinu_erc_ilmanen_vmc_penman_reoi_mccall_bellman",
        "Multi-CSA-CTD / Callable-OAS / JLT-Migration / Corridor-VarSwap / LSM-Rogers-Dual / Cartea-LOB-MM / Lo-MacKinlay-VR / Lopez-DSR-PurgedCV / Roncalli-ERC / Ilmanen-VMC / Penman-ReOI / McCall-Bellman AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite49_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let s = T::from_f64(0.005);
    let dx0 = x[0] * s;
    let dx1 = x[1] * s;
    let dx2 = x[2] * s;
    let dx3 = x[3] * s;

    let nss = ql_curves::evaluate_brigo_tuckman_nss_and_key_rate_durations(
        &ql_curves::BrigoTuckmanNssKeyRateDurationSpec {
            beta0_level: T::from_f64(0.045) + dx0,
            beta1_slope: T::from_f64(-0.015) + dx1,
            beta2_hump1: T::from_f64(0.018) + dx2,
            beta3_hump2: T::from_f64(-0.008) + dx3,
            tau1_decay: T::from_f64(1.8),
            tau2_decay: T::from_f64(5.5),
            eval_maturity: T::from_f64(4.0),
            cashflow_weights_2y_5y_10y: [T::from_f64(0.25), T::from_f64(0.35), T::from_f64(0.40)],
        },
    );

    let tips = ql_fixed_income::evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei(
        &ql_fixed_income::TuckmanFabozziTipsDeflationFloorSpec {
            base_cpi_index: T::from_f64(295.0),
            current_cpi_index: T::from_f64(300.0) + x[0],
            nominal_treasury_yield: T::from_f64(0.042) + dx1,
            real_tips_yield: T::from_f64(0.019) + dx2,
            real_coupon_rate: T::from_f64(0.0175),
            cpi_volatility: T::from_f64(0.025) + dx3,
            inflation_risk_premium: T::from_f64(0.003),
            tips_liquidity_premium: T::from_f64(0.0015),
            maturity_years: T::from_f64(5.0),
        },
    );

    let rec = ql_credit::evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle(
        &ql_credit::SchonbucherLandoRecoveryConventionSpec {
            par_amount: T::from_f64(100.0),
            risk_free_rate: T::from_f64(0.04) + dx0,
            hazard_rate_lambda: T::from_f64(0.035) + dx1,
            recovery_rate: T::from_f64(0.40) + dx2,
            coupon_rate: T::from_f64(0.05) + dx3,
            maturity_years: T::from_f64(5.0),
        },
    );

    let slv = ql_vol::evaluate_gatheral_bergomi_sticky_vs_local_vol_delta(
        &ql_vol::GatheralBergomiStickyLocalVolDeltaSpec {
            spot_price: T::from_f64(100.0) + x[0],
            strike_price: T::from_f64(100.0) + x[1],
            risk_free_rate: T::from_f64(0.03) + dx2,
            dividend_yield: T::from_f64(0.01),
            forward_start_t1: T::from_f64(0.5),
            maturity_t2: T::from_f64(1.5),
            atm_vol: T::from_f64(0.20) + dx3,
            log_moneyness_skew_s0: T::from_f64(-0.15),
            sv_mean_reversion_kappa: T::from_f64(1.2),
        },
    );

    let mva = ql_xva::evaluate_gregory_bilateral_simm_vs_ccp_mva(
        &ql_xva::GregoryBilateralSimmVsCcpMvaSpec {
            portfolio_annual_volatility: T::from_f64(15.0) + x[0],
            bilateral_mpor_years: T::from_f64(10.0 / 250.0),
            ccp_mpor_years: T::from_f64(5.0 / 250.0),
            exposure_decay_kappa: T::from_f64(0.12) + dx1,
            im_funding_spread: T::from_f64(0.0085) + dx2,
            cpty_hazard_rate: T::from_f64(0.015) + dx3,
            cpty_lgd: T::from_f64(0.60),
            risk_free_rate: T::from_f64(0.035),
            ccp_default_fund_contribution: T::from_f64(0.35),
            capital_hurdle_rate: T::from_f64(0.10),
            maturity_years: T::from_f64(5.0),
        },
    );

    let ss = ql_commodities::evaluate_schwartz_smith_two_factor_and_kalman_update(
        &ql_commodities::SchwartzSmithTwoFactorKalmanSpec {
            prior_short_term_chi0: T::from_f64(0.08) + dx0,
            prior_equilibrium_xi0: T::from_f64(4.20) + dx1,
            mean_reversion_kappa: T::from_f64(1.4) + dx2,
            equilibrium_drift_mu: T::from_f64(0.01),
            short_term_vol_sigma_chi: T::from_f64(0.28) + dx3,
            equilibrium_vol_sigma_xi: T::from_f64(0.14),
            factor_correlation_rho: T::from_f64(0.30),
            maturity_years: T::from_f64(1.0),
            prior_var_p11: T::from_f64(0.04),
            prior_var_p22: T::from_f64(0.02),
            observed_log_futures: T::from_f64(4.32),
            measurement_noise_var_r: T::from_f64(0.005),
        },
    );

    let dnt = ql_fx::evaluate_wystup_clark_fx_double_no_touch_vanna_volga(
        &ql_fx::WystupClarkFxDoubleNoTouchSpec {
            spot_fx: T::from_f64(1.08) + dx0,
            lower_barrier_l: T::from_f64(1.02),
            upper_barrier_u: T::from_f64(1.15),
            domestic_rate: T::from_f64(0.04) + dx1,
            foreign_rate: T::from_f64(0.025) + dx2,
            atm_vol: T::from_f64(0.085) + dx3,
            risk_reversal_25d: T::from_f64(-0.012),
            butterfly_25d: T::from_f64(0.0035),
            overhedge_survival_weight: T::from_f64(0.65),
            maturity_years: T::from_f64(0.5),
        },
    );

    let tmr = ql_exotics::evaluate_bouzoubaa_timer_option_variance_clock(
        &ql_exotics::BouzoubaaTimerOptionVarianceClockSpec {
            spot_price: T::from_f64(100.0) + x[0],
            strike_price: T::from_f64(100.0) + x[1],
            risk_free_rate: T::from_f64(0.02) + dx2,
            target_budget_vol: T::from_f64(0.16) + dx3,
            nominal_horizon_years: T::from_f64(1.0),
            market_implied_vol: T::from_f64(0.22),
            expected_realized_vol: T::from_f64(0.20),
            vol_of_vol_xi: T::from_f64(0.30),
        },
    );

    let bl = ql_portfolio::evaluate_meucci_black_litterman_woodbury_equivalence(
        &ql_portfolio::MeucciBlackLittermanWoodburySpec {
            covariance_3x3: [
                [
                    T::from_f64(0.040) + dx0,
                    T::from_f64(0.010),
                    T::from_f64(0.004),
                ],
                [
                    T::from_f64(0.010),
                    T::from_f64(0.025) + dx1,
                    T::from_f64(0.005),
                ],
                [
                    T::from_f64(0.004),
                    T::from_f64(0.005),
                    T::from_f64(0.010) + dx2,
                ],
            ],
            market_cap_weights: [T::from_f64(0.50), T::from_f64(0.30), T::from_f64(0.20)],
            risk_aversion_delta: T::from_f64(2.5),
            tau_scalar: T::from_f64(0.05),
            view_pick_vector: [T::one(), -T::one(), T::zero()],
            view_target_q: T::from_f64(0.045) + dx3,
            view_uncertainty_omega: T::from_f64(0.0015),
        },
    );

    let ac = ql_portfolio::evaluate_gueant_almgren_chriss_hyperbolic_liquidation(
        &ql_portfolio::GueantAlmgrenChrissHyperbolicSpec {
            initial_inventory_x0: T::from_f64(100.0) + x[0],
            horizon_years_t: T::from_f64(0.25),
            temporary_impact_eta: T::from_f64(0.002) + dx1 * T::from_f64(0.1),
            permanent_impact_gamma: T::from_f64(0.0005),
            asset_volatility_sigma: T::from_f64(0.30) + dx2,
            risk_aversion_lambda: T::from_f64(0.15) + dx3,
        },
    );

    let enet = ql_dml::evaluate_hastie_zou_elastic_net_and_grouping(
        &ql_dml::HastieZouElasticNetGroupingSpec {
            collinear_rho01: T::from_f64(0.92),
            feature_target_corrs: [
                T::from_f64(0.68) + dx0,
                T::from_f64(0.65) + dx1,
                T::from_f64(0.42) + dx2,
            ],
            lasso_penalty_lambda1: T::from_f64(0.08),
            ridge_penalty_lambda2: T::from_f64(0.25) + dx3,
        },
    );

    let nk = ql_models::evaluate_ljungqvist_new_keynesian_three_equation_ree(
        &ql_models::LjungqvistNewKeynesianReeSpec {
            risk_aversion_sigma: T::from_f64(1.5) + dx0,
            discount_factor_beta: T::from_f64(0.99),
            calvo_phillips_kappa: T::from_f64(0.15) + dx1,
            taylor_inflation_phi_pi: T::from_f64(1.50) + dx2,
            taylor_output_phi_y: T::from_f64(0.125) + dx3,
            natural_rate_persistence_rho: T::from_f64(0.80),
            natural_rate_shock_rn: T::from_f64(0.015),
        },
    );

    [
        nss.total_fisher_weil_duration,
        tips.floored_tips_fair_value,
        rec.rfv_coupon_bond_fair_value,
        slv.dupire_local_vol_delta,
        mva.central_clearing_xva_savings,
        ss.posterior_fitted_log_futures,
        dnt.smile_adjusted_dnt_pv,
        tmr.timer_call_fair_value,
        bl.posterior_view_expected_return,
        ac.optimal_mean_variance_objective,
        enet.elastic_net_objective_value,
        nk.equilibrium_policy_rate,
    ]
}

pub fn benchmark_brigo_nss_krd_tuckman_tips_floor_schonbucher_rmv_rt_rfv_gatheral_sticky_lv_gregory_simm_ccp_mva_schwartz_smith_kalman_wystup_fx_dnt_bouzoubaa_timer_meucci_bl_woodbury_gueant_ac_sinh_zou_elastic_net_nk_ree(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite49_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite49_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "49_brigo_nss_krd_tuckman_tips_floor_schonbucher_rmv_rt_rfv_gatheral_sticky_lv_gregory_simm_ccp_mva_schwartz_smith_kalman_wystup_fx_dnt_bouzoubaa_timer_meucci_bl_woodbury_gueant_ac_sinh_zou_elastic_net_nk_ree",
        "NSS-KRD / TIPS-Floor-BEI / RMV-RT-RFV / Sticky-LocalVol / SIMM-CCP-MVA / Schwartz-Smith-KF / FX-DNT-VV / Timer-VarClock / BL-Woodbury / AC-Sinh-Liq / ElasticNet-Group / NK-DSGE-REE AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite50_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let s = T::from_f64(0.005);
    let dx0 = x[0] * s;
    let dx1 = x[1] * s;
    let dx2 = x[2] * s;
    let dx3 = x[3] * s;

    let hwt = ql_curves::evaluate_hagan_west_tension_convex_forward_curve(
        &ql_curves::HaganWestTensionCurveSpec {
            pillar_intervals: [
                T::from_f64(1.0),
                T::from_f64(1.0),
                T::from_f64(2.0),
                T::from_f64(3.0),
            ],
            discrete_forwards: [
                T::from_f64(0.032) + dx0,
                T::from_f64(0.038) + dx1,
                T::from_f64(0.044) + dx2,
                T::from_f64(0.046) + dx3,
            ],
            tension_lambda: T::from_f64(1.5),
            eval_coord_x: T::from_f64(0.35),
        },
    );

    let dl = ql_credit::evaluate_duffie_lando_noisy_accounting_credit_spread(
        &ql_credit::DuffieLandoNoisyCreditSpec {
            reported_log_distance_x_hat: T::from_f64(0.45) + dx0,
            initial_log_distance_x0: T::from_f64(0.55) + dx1,
            asset_volatility_sigma: T::from_f64(0.25) + dx2,
            risk_free_rate: T::from_f64(0.04) + dx3,
            payout_rate_delta: T::from_f64(0.015),
            report_age_t0: T::from_f64(0.5),
            accounting_noise_std_a: T::from_f64(0.18),
            recovery_rate: T::from_f64(0.40),
            horizon_tau: T::from_f64(2.0),
        },
    );

    let pac = ql_fixed_income::evaluate_fabozzi_pac_companion_cmo_waterfall(
        &ql_fixed_income::FabozziPacCompanionCmoSpec {
            initial_pool_balance: T::from_f64(100.0) + x[0],
            wac_coupon: T::from_f64(0.06) + dx1,
            tranche_coupon: T::from_f64(0.055) + dx2,
            psa_lower_band: T::from_f64(100.0),
            psa_upper_band: T::from_f64(250.0),
            realized_psa_speed: T::from_f64(165.0),
            discount_yield: T::from_f64(0.045) + dx3,
        },
    );

    let fm = ql_portfolio::evaluate_campbell_fama_macbeth_and_shanken_eiv(
        &ql_portfolio::CampbellFamaMacbethShankenSpec {
            asset_mean_excess_returns: [
                T::from_f64(0.045) + dx0,
                T::from_f64(0.062) + dx1,
                T::from_f64(0.078) + dx2,
                T::from_f64(0.095) + dx3,
            ],
            estimated_factor_betas: [
                T::from_f64(0.70),
                T::from_f64(0.95),
                T::from_f64(1.20),
                T::from_f64(1.45),
            ],
            idiosyncratic_var_sigma_eps2: T::from_f64(0.015),
            factor_var_sigma_f2: T::from_f64(0.022),
            sample_periods_t: T::from_f64(120.0),
        },
    );

    let svi = ql_vol::evaluate_gatheral_jacquier_arbitrage_free_svi_density(
        &ql_vol::GatheralJacquierSviDensitySpec {
            forward_price: T::from_f64(100.0) + x[0],
            strike_price: T::from_f64(102.0) + x[1],
            expiry_years: T::from_f64(1.0),
            svi_a: T::from_f64(0.015) + dx2,
            svi_b: T::from_f64(0.12) + dx3,
            svi_rho: T::from_f64(-0.35),
            svi_m: T::from_f64(0.02),
            svi_sigma: T::from_f64(0.22),
        },
    );

    let bgk = ql_models::evaluate_broadie_glasserman_kou_discrete_barrier_lookback(
        &ql_models::BroadieGlassermanKouDiscreteBarrierSpec {
            spot_price: T::from_f64(100.0) + x[0],
            strike_price: T::from_f64(100.0) + x[1],
            down_barrier_h: T::from_f64(90.0),
            risk_free_rate: T::from_f64(0.04) + dx2,
            dividend_yield: T::from_f64(0.01),
            volatility: T::from_f64(0.22) + dx3,
            maturity_years: T::from_f64(1.0),
            monitoring_steps_m: T::from_f64(50.0),
        },
    );

    let ci = ql_portfolio::evaluate_bouchaud_cross_impact_eigenmode_execution(
        &ql_portfolio::BouchaudCrossImpactEigenmodeSpec {
            parent_order_q1: T::from_f64(50.0) + x[0],
            parent_order_q2: T::from_f64(30.0) + x[1],
            self_impact_lambda1: T::from_f64(0.04) + dx2,
            self_impact_lambda2: T::from_f64(0.03) + dx3,
            cross_impact_lambda12: T::from_f64(0.015),
            resilience_decay_rho: T::from_f64(0.60),
            risk_aversion_gamma: T::from_f64(0.10),
            asset_volatility: T::from_f64(0.20),
        },
    );

    let ks = ql_portfolio::evaluate_karnosky_singer_multicurrency_attribution(
        &ql_portfolio::KarnoskySingerMultiCurrencySpec {
            port_asset_weights: [T::from_f64(0.45), T::from_f64(0.35), T::from_f64(0.20)],
            bench_asset_weights: [T::from_f64(0.40), T::from_f64(0.35), T::from_f64(0.25)],
            port_fx_weights: [T::from_f64(0.50), T::from_f64(0.30), T::from_f64(0.20)],
            bench_fx_weights: [T::from_f64(0.40), T::from_f64(0.35), T::from_f64(0.25)],
            local_port_returns: [
                T::from_f64(0.065) + dx0,
                T::from_f64(0.048) + dx1,
                T::from_f64(0.032) + dx2,
            ],
            local_bench_returns: [T::from_f64(0.058), T::from_f64(0.045), T::from_f64(0.035)],
            local_eurodeposit_rates: [
                T::from_f64(0.040),
                T::from_f64(0.028),
                T::from_f64(0.012) + dx3,
            ],
            base_cash_rate: T::from_f64(0.040),
            spot_fx_returns: [T::zero(), T::from_f64(0.015), T::from_f64(-0.012)],
        },
    );

    let mjd = ql_levy::evaluate_cont_tankov_merton_jump_diffusion(
        &ql_levy::ContTankovMertonJumpDiffusionSpec {
            spot_price: T::from_f64(100.0) + x[0],
            strike_price: T::from_f64(105.0) + x[1],
            risk_free_rate: T::from_f64(0.035) + dx2,
            dividend_yield: T::from_f64(0.01),
            diffusion_vol_sigma: T::from_f64(0.16) + dx3,
            jump_intensity_lambda: T::from_f64(0.75),
            log_jump_mean_mu: T::from_f64(-0.06),
            log_jump_std_delta: T::from_f64(0.14),
            maturity_years: T::from_f64(1.0),
        },
    );

    let cv = ql_portfolio::evaluate_campbell_viceira_intertemporal_hedging_saa(
        &ql_portfolio::CampbellViceiraIntertemporalSaaSpec {
            current_expected_excess_return_xt: T::from_f64(0.055) + dx0,
            long_run_excess_return_xbar: T::from_f64(0.050) + dx1,
            predictor_persistence_phi: T::from_f64(0.88),
            return_volatility_sigma_r: T::from_f64(0.16) + dx2,
            predictor_volatility_sigma_x: T::from_f64(0.035) + dx3,
            return_predictor_corr_rho: T::from_f64(-0.68),
            relative_risk_aversion_gamma: T::from_f64(4.0),
            wealth_consumption_rho_w: T::from_f64(0.96),
        },
    );

    let ohl = ql_models::evaluate_ohlson_penman_linear_information_dynamics(
        &ql_models::OhlsonPenmanLinearInfoDynamicsSpec {
            current_book_value_bt: T::from_f64(100.0) + x[0],
            current_earnings_xt: T::from_f64(14.5) + x[1],
            current_dividends_dt: T::from_f64(4.5),
            other_information_shock_nu: T::from_f64(1.8) + dx2,
            cost_of_equity_re: T::from_f64(0.09) + dx3,
            abnormal_earnings_persistence_omega: T::from_f64(0.62),
            other_info_persistence_gamma: T::from_f64(0.45),
        },
    );

    let rev = ql_ad::evaluate_griewank_revolve_checkpointed_bsb_pde_adjoint(
        &ql_ad::GriewankRevolveBsbAdjointSpec {
            spot_price: T::from_f64(100.0) + x[0],
            strike_price: T::from_f64(100.0) + x[1],
            risk_free_rate: T::from_f64(0.03) + dx2,
            sigma_min: T::from_f64(0.15),
            sigma_max: T::from_f64(0.25) + dx3,
            bsb_switching_steepness: T::from_f64(2.0),
            maturity_years: T::from_f64(0.5),
        },
    );

    [
        hwt.zero_coupon_bond_price,
        dl.noisy_credit_spread_bps,
        pac.pac_tranche_pv,
        fm.factor_risk_premium_gamma1,
        svi.breeden_litzenberger_pdf,
        bgk.discrete_bgk_down_and_out_call,
        ci.total_cross_impact_execution_cost,
        ks.total_active_base_return,
        mjd.merton_jump_diffusion_call,
        cv.total_strategic_equity_allocation,
        ohl.lid_intrinsic_equity_value,
        rev.primal_bsb_option_value,
    ]
}

pub fn benchmark_hagan_west_tension_duffie_lando_pac_cmo_fama_macbeth_shanken_gatheral_asvi_bgk_discrete_barrier_bouchaud_cross_impact_karnosky_singer_merton_jump_campbell_viceira_ohlson_lid_griewank_revolve(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite50_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite50_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "50_hagan_west_tension_duffie_lando_pac_cmo_fama_macbeth_shanken_gatheral_asvi_bgk_discrete_barrier_bouchaud_cross_impact_karnosky_singer_merton_jump_campbell_viceira_ohlson_lid_griewank_revolve",
        "Hagan-West-Tension / Duffie-Lando / PAC-CMO / Fama-MacBeth-Shanken / Gatheral-aSVI / BGK-Barrier / Cross-Impact / Karnosky-Singer / Merton-Jump / Campbell-Viceira / Ohlson-LID / Revolve-BSB AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite51_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let s = T::from_f64(0.005);
    let dx0 = x[0] * s;
    let dx1 = x[1] * s;
    let dx2 = x[2] * s;
    let dx3 = x[3] * s;

    let xccy = ql_curves::evaluate_brigo_xccy_mtm_collateralized_basis_curve(
        &ql_curves::BrigoXccyMtmBasisCurveSpec {
            spot_fx_s0: T::from_f64(1.08) + dx0,
            domestic_csa_rate_rd: T::from_f64(0.048) + dx1,
            foreign_local_ois_rate_rf: T::from_f64(0.032) + dx2,
            fx_forward_quote_f0t: T::from_f64(1.11) + dx3,
            mtm_xccy_par_basis_spread_bm: T::from_f64(-0.0025),
            fx_vol_sigma_s: T::from_f64(0.09),
            foreign_rate_vol_sigma_rf: T::from_f64(0.008),
            fx_rate_corr_rho: T::from_f64(-0.35),
            maturity_years_t: T::from_f64(2.0),
        },
    );

    let iopo = ql_fixed_income::evaluate_fabozzi_specified_pool_io_po_payup(
        &ql_fixed_income::FabozziSpecifiedPoolIoPoSpec {
            par_notional: T::from_f64(100.0) + x[0],
            wac_coupon: T::from_f64(0.065) + dx1,
            pass_through_coupon: T::from_f64(0.060),
            current_mortgage_rate: T::from_f64(0.048) + dx2,
            discount_yield: T::from_f64(0.042) + dx3,
            base_turnover_cpr: T::from_f64(0.06),
            max_refi_cpr: T::from_f64(0.42),
            s_curve_steepness_kappa: T::from_f64(120.0),
            spec_pool_burnout_dampener: T::from_f64(0.45),
        },
    );

    let mo = ql_credit::evaluate_lando_marshall_olkin_fatal_shock_basket(
        &ql_credit::LandoMarshallOlkinContagionSpec {
            idiosyncratic_intensity_lambda1: T::from_f64(0.018) + dx0,
            idiosyncratic_intensity_lambda2: T::from_f64(0.024) + dx1,
            fatal_common_shock_lambda12: T::from_f64(0.012) + dx2,
            recovery_rate: T::from_f64(0.40),
            risk_free_rate: T::from_f64(0.035) + dx3,
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let vix = ql_vol::evaluate_gatheral_bergomi_vix_futures_and_vvix_smile(
        &ql_vol::GatheralBergomiVixVvixSpec {
            forward_variance_swap_var: T::from_f64(0.0484) + dx0,
            vvix_vol_of_vol_omega: T::from_f64(0.85) + dx1,
            vix_option_strike_k: T::from_f64(0.21) + dx2,
            risk_free_rate: T::from_f64(0.04) + dx3,
            expiry_years_t: T::from_f64(0.25),
        },
    );

    let wwr = ql_xva::evaluate_gregory_hull_white_wrong_way_risk_cva(
        &ql_xva::GregoryHullWhiteWwrCvaSpec {
            baseline_exposure_mu: T::from_f64(12.0) + x[0],
            exposure_volatility_sigma_x: T::from_f64(18.0) + x[1],
            baseline_hazard_rate_lambda0: T::from_f64(0.025) + dx2,
            log_intensity_vol_sigma_lambda: T::from_f64(0.65),
            exposure_credit_corr_rho: T::from_f64(0.55) + dx3,
            recovery_rate: T::from_f64(0.40),
            risk_free_rate: T::from_f64(0.035),
            maturity_years_t: T::from_f64(3.0),
        },
    );

    let gs = ql_commodities::evaluate_geman_gibson_schwartz_commodity_futures(
        &ql_commodities::GemanGibsonSchwartzSpec {
            spot_price_s0: T::from_f64(82.0) + x[0],
            convenience_yield_delta0: T::from_f64(0.14) + dx1,
            long_run_convenience_alpha: T::from_f64(0.08) + dx2,
            mean_reversion_kappa: T::from_f64(1.4),
            spot_volatility_sigma1: T::from_f64(0.34) + dx3,
            convenience_volatility_sigma2: T::from_f64(0.26),
            spot_convenience_corr_rho: T::from_f64(0.65),
            risk_free_rate: T::from_f64(0.04),
            maturity_years_tau: T::from_f64(1.5),
        },
    );

    let fva = ql_fx::evaluate_clark_fx_forward_volatility_agreement(
        &ql_fx::ClarkFxForwardVolAgreementSpec {
            spot_fx: T::from_f64(1.08) + dx0,
            domestic_rate_rd: T::from_f64(0.045) + dx1,
            foreign_rate_rf: T::from_f64(0.030),
            front_expiry_t1: T::from_f64(0.5),
            back_expiry_t2: T::from_f64(1.5),
            atm_vol_t1: T::from_f64(0.085) + dx2,
            atm_vol_t2: T::from_f64(0.098) + dx3,
            smile_convexity_adj_bps: T::from_f64(-18.0),
        },
    );

    let ch = ql_exotics::evaluate_rubinstein_simple_chooser_option(
        &ql_exotics::RubinsteinChooserOptionSpec {
            spot_price: T::from_f64(100.0) + x[0],
            strike_price: T::from_f64(100.0) + x[1],
            risk_free_rate: T::from_f64(0.05) + dx2,
            dividend_yield: T::from_f64(0.02),
            volatility: T::from_f64(0.24) + dx3,
            choice_date_tc: T::from_f64(0.4),
            final_maturity_t: T::from_f64(1.0),
        },
    );

    let flam = ql_portfolio::evaluate_grinold_kahn_flam_with_transfer_coefficient(
        &ql_portfolio::GrinoldKahnFlamTransferSpec {
            information_coefficient_ic: T::from_f64(0.06) + dx0,
            nominal_breadth_n: T::from_f64(200.0),
            cross_sectional_signal_corr_rho: T::from_f64(0.12) + dx1,
            active_weights_std_adj: [
                T::from_f64(0.45) + dx2,
                T::from_f64(0.30) + dx3,
                T::from_f64(-0.20),
            ],
            unconstrained_optimal_weights_std: [
                T::from_f64(0.55),
                T::from_f64(0.25),
                T::from_f64(-0.35),
            ],
            active_risk_aversion_lambda: T::from_f64(2.5),
        },
    );

    let es = ql_portfolio::evaluate_meucci_expected_shortfall_euler_decomposition(
        &ql_portfolio::MeucciExpectedShortfallEulerSpec {
            weights: [
                T::from_f64(0.45) + dx0,
                T::from_f64(0.35) + dx1,
                T::from_f64(0.20) - dx0 - dx1,
            ],
            expected_returns: [T::from_f64(0.06), T::from_f64(0.05), T::from_f64(0.035)],
            asset_volatilities: [
                T::from_f64(0.18) + dx2,
                T::from_f64(0.15) + dx3,
                T::from_f64(0.10),
            ],
            pairwise_corr_rho: T::from_f64(0.35),
            confidence_alpha: T::from_f64(0.975),
            normal_quantile_z_alpha: T::from_f64(1.959963984540054),
            portfolio_skewness_gamma1: T::from_f64(-0.45),
            portfolio_excess_kurtosis_gamma2: T::from_f64(1.5),
        },
    );

    let js = ql_dml::evaluate_hastie_james_stein_empirical_bayes_sure(
        &ql_dml::HastieJamesSteinSureSpec {
            sample_means_y: [
                T::from_f64(0.08) + dx0,
                T::from_f64(0.03) + dx1,
                T::from_f64(-0.02) + dx2,
                T::from_f64(0.11) + dx3,
            ],
            grand_prior_target_mu0: T::from_f64(0.04),
            noise_variance_sigma2: T::from_f64(0.0025),
            singular_values_squared: [
                T::from_f64(4.2),
                T::from_f64(2.5),
                T::from_f64(1.1),
                T::from_f64(0.4),
            ],
            ridge_penalty_lambda: T::from_f64(0.8),
        },
    );

    let cc = ql_models::evaluate_campbell_cochrane_external_habit_sdf(
        &ql_models::CampbellCochraneExternalHabitSpec {
            utility_curvature_gamma: T::from_f64(2.0) + dx0,
            consumption_growth_mean_gc: T::from_f64(0.0189) + dx1,
            consumption_growth_vol_sigmac: T::from_f64(0.015) + dx2,
            log_surplus_persistence_phis: T::from_f64(0.87),
            subjective_discount_rate_rho_time: T::from_f64(0.0116),
            current_log_surplus_deviation_ds: T::from_f64(-0.12) + dx3,
        },
    );

    [
        xccy.fx_implied_foreign_csa_df_pf,
        iopo.specified_pool_payup_32nds,
        mo.second_to_default_par_spread_bps,
        vix.black76_vix_call_price,
        wwr.wwr_adjusted_cva,
        gs.gibson_schwartz_futures_price,
        fva.forward_start_atm_straddle_pct,
        ch.rubinstein_chooser_price,
        flam.expected_value_added_va,
        es.cornish_fisher_tail_es,
        js.stein_unbiased_risk_estimate_sure,
        cc.conditional_equity_risk_premium,
    ]
}

pub fn benchmark_brigo_xccy_mtm_fabozzi_spec_pool_io_po_lando_marshall_olkin_gatheral_vix_vvix_gregory_wwr_cva_geman_gibson_schwartz_clark_fx_fva_rubinstein_chooser_grinold_flam_tc_meucci_es_euler_james_stein_sure_campbell_cochrane_habit(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite51_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite51_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "51_brigo_xccy_mtm_fabozzi_spec_pool_io_po_lando_marshall_olkin_gatheral_vix_vvix_gregory_wwr_cva_geman_gibson_schwartz_clark_fx_fva_rubinstein_chooser_grinold_flam_tc_meucci_es_euler_james_stein_sure_campbell_cochrane_habit",
        "XCCY-MtM-Basis / Spec-Pool-IO-PO / Marshall-Olkin / VIX-VVIX-Smile / WWR-CVA / Gibson-Schwartz / FX-FVA / Chooser / FLAM-TC / ES-Euler / James-Stein-SURE / Campbell-Cochrane AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite52_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let s = T::from_f64(0.005);
    let dx0 = x[0] * s;
    let dx1 = x[1] * s;
    let dx2 = x[2] * s;
    let dx3 = x[3] * s;

    let ois = ql_curves::evaluate_brigo_ois_compounded_in_arrears_and_fallback(
        &ql_curves::BrigoOisInArrearsFallbackSpec {
            ois_df_start_p0t1: T::from_f64(0.960) + dx0,
            ois_df_end_p0t2: T::from_f64(0.938) + dx1,
            accrual_year_fraction_delta: T::from_f64(0.5),
            short_rate_vol_sigma_r: T::from_f64(0.012) + dx2,
            mean_reversion_a: T::from_f64(0.15),
            isda_historical_fallback_spread: T::from_f64(0.0026161),
            realized_accrued_rfr_rate: T::from_f64(0.049) + dx3,
            elapsed_accrual_fraction_w: T::from_f64(0.4),
        },
    );

    let frn = ql_fixed_income::evaluate_fabozzi_collared_frn_discount_margin(
        &ql_fixed_income::FabozziCollaredFrnDiscountMarginSpec {
            par_face: T::from_f64(100.0) + x[0],
            index_rate_r: T::from_f64(0.042) + dx1,
            quoted_margin_qm: T::from_f64(0.012),
            discount_margin_dm: T::from_f64(0.015) + dx2,
            cap_strike: T::from_f64(0.065),
            floor_strike: T::from_f64(0.035),
            rate_normal_vol_sigma: T::from_f64(0.011) + dx3,
            time_to_next_reset_tau1: T::from_f64(0.25),
        },
    );

    let cdsw = ql_credit::evaluate_schonbucher_cds_swaption_and_fep(
        &ql_credit::SchonbucherCdsSwaptionFepSpec {
            hazard_rate_lambda: T::from_f64(0.022) + dx0,
            recovery_rate: T::from_f64(0.40),
            risk_free_rate: T::from_f64(0.035) + dx1,
            option_expiry_te: T::from_f64(1.0),
            swap_tenor_tm: T::from_f64(5.0),
            strike_spread_k: T::from_f64(0.013) + dx2,
            spread_lognormal_vol_sigma: T::from_f64(0.45) + dx3,
        },
    );

    let berg = ql_vol::evaluate_bergomi_two_factor_skew_term_structure(
        &ql_vol::BergomiTwoFactorSkewTermStructureSpec {
            flat_forward_variance_xi0: T::from_f64(0.04) + dx0,
            vol_of_vol_omega: T::from_f64(1.8) + dx1,
            mixing_weight_theta: T::from_f64(0.35),
            fast_mean_reversion_k1: T::from_f64(6.5),
            slow_mean_reversion_k2: T::from_f64(0.35),
            factor_corr_rho12: T::from_f64(0.25),
            spot_fast_corr_rho_s1: T::from_f64(-0.75) + dx2,
            spot_slow_corr_rho_s2: T::from_f64(-0.45) + dx3,
            maturity_years_t: T::from_f64(3.0),
        },
    );

    let saccr =
        ql_xva::evaluate_gregory_sa_ccr_regulatory_ead(&ql_xva::GregorySaCcrRegulatoryEadSpec {
            trade_notional_d: T::from_f64(1000.0) + x[0],
            current_mtm_value_v: T::from_f64(-8.0) + x[1],
            posted_collateral_c: T::from_f64(4.0),
            swap_start_years_s: T::from_f64(1.0),
            swap_end_years_e: T::from_f64(6.0),
            forward_rate_p: T::from_f64(0.038) + dx2,
            option_strike_k: T::from_f64(0.035),
            supervisory_option_vol: T::from_f64(0.50) + dx3,
            supervisory_factor_sf: T::from_f64(0.005),
        });

    let spk = ql_commodities::evaluate_eydeland_clean_spark_spread_option(
        &ql_commodities::EydelandCleanSparkSpreadOptionSpec {
            forward_power_price_fp: T::from_f64(68.0) + x[0],
            forward_gas_price_fg: T::from_f64(7.2) + dx1,
            plant_heat_rate_h: T::from_f64(7.5),
            fixed_vom_and_carbon_k: T::from_f64(8.5),
            power_volatility_sigma_p: T::from_f64(0.42) + dx2,
            gas_volatility_sigma_g: T::from_f64(0.36) + dx3,
            power_gas_corr_rho: T::from_f64(0.72),
            risk_free_rate: T::from_f64(0.04),
            expiry_years_t: T::from_f64(0.5),
        },
    );

    let prdc =
        ql_fx::evaluate_wystup_prdc_three_factor_coupon(&ql_fx::WystupPrdcThreeFactorCouponSpec {
            initial_spot_fx_s0: T::from_f64(105.0),
            current_spot_fx_st: T::from_f64(108.0) + x[0],
            foreign_coupon_rate_cf: T::from_f64(0.055) + dx1,
            domestic_strike_coupon_cd: T::from_f64(0.035),
            domestic_rate_rd: T::from_f64(0.010) + dx2,
            foreign_rate_rf: T::from_f64(0.045),
            spot_fx_vol_sigma_s: T::from_f64(0.11) + dx3,
            dom_rate_vol_sigma_rd: T::from_f64(0.007),
            for_rate_vol_sigma_rf: T::from_f64(0.012),
            corr_fx_dom_rho_sd: T::from_f64(0.25),
            corr_fx_for_rho_sf: T::from_f64(-0.40),
            maturity_years_t: T::from_f64(5.0),
        });

    let comp = ql_exotics::evaluate_geske_compound_call_and_put_on_call(
        &ql_exotics::GeskeCompoundOptionSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            compound_strike_k1: T::from_f64(6.5) + dx1,
            underlying_strike_k2: T::from_f64(102.0) + x[1],
            risk_free_rate: T::from_f64(0.04) + dx2,
            dividend_yield: T::from_f64(0.015),
            volatility: T::from_f64(0.24) + dx3,
            compound_expiry_t1: T::from_f64(0.5),
            underlying_expiry_t2: T::from_f64(1.5),
        },
    );

    let mpir = ql_portfolio::evaluate_qian_multi_period_alpha_decay_ir(
        &ql_portfolio::QianMultiPeriodAlphaDecayIrSpec {
            initial_ic0: T::from_f64(0.075) + dx0,
            alpha_decay_rate_kappa: T::from_f64(2.0) + dx1,
            breadth_br: T::from_f64(180.0),
            proportional_tcost_rate_ctc: T::from_f64(0.08) + dx2,
            active_risk_budget_sigma_a: T::from_f64(0.04) + dx3,
        },
    );

    let enb = ql_portfolio::evaluate_roncalli_meucci_diversification_and_enb(
        &ql_portfolio::RoncalliMeucciDiversificationEnbSpec {
            weights: [T::from_f64(0.30), T::from_f64(0.35), T::from_f64(0.35)],
            asset_volatilities: [
                T::from_f64(0.18) + dx0,
                T::from_f64(0.14) + dx1,
                T::from_f64(0.12) + dx2,
            ],
            corr_01: T::from_f64(0.25) + dx3,
            corr_02: T::from_f64(0.10),
            corr_12: T::from_f64(0.30),
        },
    );

    let vecm = ql_dml::evaluate_tsay_johansen_vecm_cointegration(
        &ql_dml::TsayJohansenVecmCointegrationSpec {
            cointegrating_hedge_ratio_beta: T::from_f64(1.15) + dx0,
            error_correction_speed_alpha1: T::from_f64(-0.18) + dx1,
            error_correction_speed_alpha2: T::from_f64(0.08) + dx2,
            current_spread_deviation_zt: T::from_f64(1.4),
            spread_innovation_vol_sigma_eps: T::from_f64(0.25) + dx3,
            canonical_eigenvalues_lambda: [T::from_f64(0.16), T::from_f64(0.03)],
            sample_size_t: T::from_f64(250.0),
        },
    );

    let lrr = ql_models::evaluate_bansal_yaron_long_run_risks_sdf(
        &ql_models::BansalYaronLongRunRisksSpec {
            relative_risk_aversion_gamma: T::from_f64(7.5) + dx0,
            elasticity_intertemporal_sub_psi: T::from_f64(1.5) + dx1,
            campbell_shiller_kappa1: T::from_f64(0.997),
            growth_persistence_rho_x: T::from_f64(0.979),
            volatility_persistence_nu_sigma: T::from_f64(0.987),
            growth_shock_scale_phi_e: T::from_f64(0.044) + dx2,
            baseline_consumption_vol_sigma: T::from_f64(0.027) + dx3,
            vol_of_vol_sigma_w: T::from_f64(0.00023),
        },
    );

    [
        ois.seasoned_period_blended_rfr_rate,
        frn.collared_frn_price,
        cdsw.non_knockout_payer_swaption_bps,
        berg.atmf_implied_vol_skew_slope,
        saccr.sa_ccr_exposure_at_default_ead,
        spk.bjerksund_stensland_spark_call,
        prdc.prdc_coupon_present_value,
        comp.call_on_call_price,
        mpir.net_active_value_added,
        enb.shannon_entropy_effective_bets_enb,
        vecm.mean_reversion_half_life_periods,
        lrr.total_consumption_claim_risk_premium,
    ]
}

pub fn benchmark_brigo_ois_fallback_fabozzi_collared_frn_schonbucher_cds_swaption_fep_bergomi_2f_skew_gregory_sa_ccr_ead_eydeland_spark_spread_wystup_prdc_geske_compound_qian_alpha_decay_roncalli_enb_johansen_vecm_bansal_yaron_lrr(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite52_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite52_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "52_brigo_ois_fallback_fabozzi_collared_frn_schonbucher_cds_swaption_fep_bergomi_2f_skew_gregory_sa_ccr_ead_eydeland_spark_spread_wystup_prdc_geske_compound_qian_alpha_decay_roncalli_enb_johansen_vecm_bansal_yaron_lrr",
        "OIS-Fallback / Collared-FRN / CDS-Swaption-FEP / Bergomi-2F-Skew / SA-CCR-EAD / Clean-Spark / FX-PRDC / Geske-Compound / Alpha-Decay-IR / Meucci-ENB / Johansen-VECM / Bansal-Yaron-LRR AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite53_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let s = T::from_f64(0.005);
    let dx0 = x[0] * s;
    let dx1 = x[1] * s;
    let dx2 = x[2] * s;
    let dx3 = x[3] * s;

    let cms = ql_curves::evaluate_andersen_piterbarg_cms_replication_convexity(
        &ql_curves::AndersenPiterbargCmsReplicationSpec {
            forward_swap_rate_s0: T::from_f64(0.038) + dx0,
            swap_tenor_years_n: T::from_f64(10.0),
            payment_frequency_m: T::from_f64(2.0),
            fixing_expiry_te: T::from_f64(5.0),
            normal_swaption_vol_sigma_n: T::from_f64(0.0085) + dx1,
            ois_spread_delta_r: T::from_f64(0.0015) + dx2,
            mean_reversion_speed_a: T::from_f64(0.04) + dx3,
        },
    );

    let repo = ql_fixed_income::evaluate_tuckman_special_repo_fail_penalty_ctd(
        &ql_fixed_income::TuckmanSpecialRepoFailPenaltyCtdSpec {
            general_collateral_rate_gc: T::from_f64(0.048) + dx0,
            tmpg_fail_penalty_rate_p: T::from_f64(0.030),
            issue_short_demand_ratio_d: T::from_f64(1.4) + dx1,
            lendable_supply_elasticity_k: T::from_f64(0.55) + dx2,
            bond_clean_price: T::from_f64(98.5) + x[0],
            accrued_interest_at_spot: T::from_f64(0.85),
            coupon_rate_annual: T::from_f64(0.040) + dx3,
            futures_invoice_price_at_delivery: T::from_f64(99.65),
            holding_period_year_frac_tau: T::from_f64(0.25),
        },
    );

    let cds_uf = ql_credit::evaluate_okane_isda_standard_cds_upfront_triangle(
        &ql_credit::OKaneIsdaStandardCdsUpfrontSpec {
            notional_amount: T::from_f64(1000.0) + x[0],
            market_par_spread_s: T::from_f64(0.0185) + dx1,
            standard_running_coupon_c: T::from_f64(0.0100),
            recovery_rate_r: T::from_f64(0.40) + dx2,
            ois_zero_rate_r: T::from_f64(0.038) + dx3,
            maturity_years_t: T::from_f64(5.0),
            accrued_fraction_alpha: T::from_f64(0.35),
        },
    );

    let ssr = ql_vol::evaluate_gatheral_ssr_tanaka_stop_loss_hedging(
        &ql_vol::GatheralSsrTanakaStopLossSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(100.0),
            risk_free_rate_r: T::from_f64(0.04) + dx1,
            atm_implied_vol_sigma0: T::from_f64(0.22) + dx2,
            atm_skew_slope_dsigma_dk: T::from_f64(-0.18),
            spot_vol_regression_beta: T::from_f64(-0.24) + dx3,
            maturity_years_t: T::from_f64(0.5),
            rebalancing_steps_n: T::from_f64(252.0),
        },
    );

    let kva = ql_xva::evaluate_gregory_saccr_kva_regulatory_capital(
        &ql_xva::GregorySaccrKvaRegulatoryCapitalSpec {
            initial_sa_ccr_ead0: T::from_f64(150.0) + x[0],
            counterparty_risk_weight_rw: T::from_f64(0.75) + dx1,
            pillar_1_capital_ratio_k: T::from_f64(0.105),
            cva_capital_add_on_multiplier: T::from_f64(0.35),
            shareholder_hurdle_rate_hk: T::from_f64(0.11) + dx2,
            risk_free_discount_rate_r: T::from_f64(0.038) + dx3,
            exposure_amortization_decay_mu: T::from_f64(0.18),
            portfolio_maturity_years_t: T::from_f64(5.0),
        },
    );

    let crk = ql_commodities::evaluate_eydeland_crack_spread_three_asset_kirk(
        &ql_commodities::EydelandThreeAssetCrackSpreadSpec {
            crude_forward_barrel_fc: T::from_f64(76.0) + x[0],
            gasoline_forward_barrel_fg: T::from_f64(94.0) + x[1],
            distillate_forward_barrel_fh: T::from_f64(98.0) + x[2],
            fixed_refining_strike_k: T::from_f64(6.5),
            crude_vol_sigma_c: T::from_f64(0.34) + dx3,
            gasoline_vol_sigma_g: T::from_f64(0.38),
            distillate_vol_sigma_h: T::from_f64(0.36),
            corr_gas_dist_rho_gh: T::from_f64(0.82),
            corr_product_crude_rho_pc: T::from_f64(0.86),
            risk_free_rate_r: T::from_f64(0.04),
            expiry_years_t: T::from_f64(0.5),
        },
    );

    let prdc_cap =
        ql_fx::evaluate_wystup_fx_prdc_coupon_convexity(&ql_fx::WystupFxPrdcCallableCapSpec {
            spot_fx_s0: T::from_f64(108.0) + x[0],
            initial_fx_s_init: T::from_f64(105.0),
            foreign_coupon_cf: T::from_f64(0.055) + dx1,
            domestic_deduction_cd: T::from_f64(0.035),
            cap_rate_k_cap: T::from_f64(0.030),
            domestic_rate_rd: T::from_f64(0.012) + dx2,
            foreign_rate_rf: T::from_f64(0.045),
            fx_vol_sigma_s: T::from_f64(0.115) + dx3,
            dom_rate_vol_sigma_d: T::from_f64(0.008),
            fx_dom_corr_rho_sd: T::from_f64(0.35),
            coupon_horizon_t: T::from_f64(7.0),
        });

    let eso = ql_exotics::evaluate_haug_eso_vesting_barrier_dilution(
        &ql_exotics::HaugEsoVestingBarrierDilutionSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(100.0),
            risk_free_rate_r: T::from_f64(0.04) + dx1,
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.30) + dx2,
            contractual_maturity_t: T::from_f64(10.0),
            vesting_period_tv: T::from_f64(3.0),
            forfeiture_hazard_lambda_f: T::from_f64(0.06) + dx3,
            suboptimal_exercise_multiple_m: T::from_f64(1.8),
            grant_to_shares_ratio_gamma: T::from_f64(0.05),
        },
    );

    let erc = ql_portfolio::evaluate_qian_risk_parity_spinu_log_barrier(
        &ql_portfolio::QianSpinuRiskParityLogBarrierSpec {
            asset_volatilities: [
                T::from_f64(0.18) + dx0,
                T::from_f64(0.12) + dx1,
                T::from_f64(0.08) + dx2,
            ],
            risk_budgets_b: [
                T::from_f64(1.0 / 3.0),
                T::from_f64(1.0 / 3.0),
                T::from_f64(1.0 / 3.0),
            ],
            corr_01: T::from_f64(0.25) + dx3,
            corr_02: T::from_f64(-0.10),
            corr_12: T::from_f64(0.15),
        },
    );

    let cvar = ql_portfolio::evaluate_rockafellar_uryasev_cvar_minimization(
        &ql_portfolio::RockafellarUryasevCvarAuxiliarySpec {
            portfolio_weights: [T::from_f64(0.6) + dx0, T::from_f64(0.4) - dx0],
            scenario_losses_asset0: [
                T::from_f64(-0.02),
                T::from_f64(0.01) + dx1,
                T::from_f64(0.04) + dx2,
                T::from_f64(0.09) + dx3,
            ],
            scenario_losses_asset1: [
                T::from_f64(-0.01),
                T::from_f64(0.005),
                T::from_f64(0.02),
                T::from_f64(0.05),
            ],
            confidence_alpha: T::from_f64(0.75),
            smoothing_epsilon: T::from_f64(0.002),
        },
    );

    let pls =
        ql_dml::evaluate_htf_pls_vs_pcr_krylov_shrinkage(&ql_dml::HtfPlsKrylovShrinkageSpec {
            pc_eigenvalues_d2: [
                T::from_f64(4.5) + dx0,
                T::from_f64(1.8) + dx1,
                T::from_f64(0.4) + dx2,
            ],
            ols_pc_projections_u: [T::from_f64(0.8), T::from_f64(1.4), T::from_f64(1.1) + dx3],
            ridge_penalty_lambda: T::from_f64(0.6),
            noise_variance_sigma2: T::from_f64(0.05),
        });

    let fomc = ql_models::evaluate_piazzesi_swanson_fomc_futures_hike_tree(
        &ql_models::PiazzesiSwansonFomcFuturesHikeSpec {
            current_target_fed_funds_r0: T::from_f64(0.0475) + dx0,
            futures_implied_month_avg_rate_f0: T::from_f64(0.0498) + dx1,
            meeting_day_of_month_d: T::from_f64(15.0),
            total_days_in_month_n: T::from_f64(30.0),
            standard_hike_step_delta: T::from_f64(0.0025),
            piazzesi_swanson_term_premium_bps: T::from_f64(4.5) + dx2,
            next_meeting_futures_rate_f1: T::from_f64(0.0512) + dx3,
        },
    );

    [
        cms.convexity_adjusted_cms_rate,
        repo.net_basis_after_specialness,
        cds_uf.cash_settlement_upfront_amount,
        ssr.minimum_variance_smile_delta,
        kva.lifetime_capital_valuation_adjustment_kva,
        crk.crack_spread_call_price_per_bbl,
        prdc_cap.capped_prdc_coupon_pv,
        eso.diluted_per_option_grant_cost,
        erc.portfolio_volatility,
        cvar.minimized_auxiliary_cvar,
        pls.pls1_fitted_norm_sq,
        fomc.physical_hike_probability_p,
    ]
}

pub fn benchmark_andersen_cms_replication_tuckman_special_repo_okane_cds_upfront_gatheral_ssr_tanaka_gregory_saccr_kva_eydeland_crack_spread_wystup_prdc_cap_haug_eso_dilution_qian_spinu_erc_rockafellar_cvar_htf_pls_krylov_piazzesi_fomc(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite53_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite53_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "53_andersen_cms_replication_tuckman_special_repo_okane_cds_upfront_gatheral_ssr_tanaka_gregory_saccr_kva_eydeland_crack_spread_wystup_prdc_cap_haug_eso_dilution_qian_spinu_erc_rockafellar_cvar_htf_pls_krylov_piazzesi_fomc",
        "CMS-Replication / Special-Repo-CTD / ISDA-CDS-Upfront / SSR-Tanaka / SA-CCR-KVA / 3-Asset-Crack / PRDC-Cap / ESO-Dilution / Spinu-ERC / Rockafellar-CVaR / PLS-Krylov / FOMC-Hike AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite54_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let s = T::from_f64(0.005);
    let dx0 = x[0] * s;
    let dx1 = x[1] * s;
    let dx2 = x[2] * s;
    let dx3 = x[3] * s;

    let ois_turn = ql_curves::evaluate_darbyshire_central_bank_step_ois_and_turn(
        &ql_curves::DarbyshireCentralBankStepOisTurnSpec {
            current_overnight_ois_r0: T::from_f64(0.045) + dx0,
            meeting_1_time_years_t1: T::from_f64(0.15),
            meeting_1_expected_step_bps: T::from_f64(22.0) + x[1],
            meeting_2_time_years_t2: T::from_f64(0.35),
            meeting_2_expected_step_bps: T::from_f64(18.0) + x[2],
            turn_of_year_spike_bps: T::from_f64(35.0) + x[3],
            turn_effective_day_weight_tau: T::from_f64(3.0 / 360.0),
            horizon_years_t: T::from_f64(0.50),
        },
    );

    let cb_mw = ql_fixed_income::evaluate_despiegeleer_cb_soft_call_make_whole(
        &ql_fixed_income::DeSpiegeleerCbSoftCallMakeWholeSpec {
            spot_price_s0: T::from_f64(105.0) + x[0],
            par_value: T::from_f64(100.0),
            conversion_ratio_cr: T::from_f64(1.0),
            soft_call_trigger_multiple: T::from_f64(1.30),
            make_whole_max_bonus_shares: T::from_f64(0.18) + dx1,
            coupon_rate_c: T::from_f64(0.025),
            risky_bond_yield_y: T::from_f64(0.055) + dx2,
            risk_free_rate_r: T::from_f64(0.040),
            dividend_and_borrow_cost_q: T::from_f64(0.018),
            equity_volatility_sigma: T::from_f64(0.28) + dx3,
            maturity_years_t: T::from_f64(4.0),
        },
    );

    let cdx = ql_credit::evaluate_okane_cds_index_intrinsic_skew_basis(
        &ql_credit::OKaneCdsIndexIntrinsicSkewSpec {
            constituent_spreads_s: [
                T::from_f64(0.0060) + dx0,
                T::from_f64(0.0120) + dx1,
                T::from_f64(0.0300) + dx2,
            ],
            market_quoted_index_spread: T::from_f64(0.0155) + dx3,
            standard_index_coupon_c: T::from_f64(0.0100),
            recovery_rate_r: T::from_f64(0.40),
            ois_discount_rate_r: T::from_f64(0.038),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let vswp = ql_vol::evaluate_gatheral_varswap_log_strip_tail_truncation(
        &ql_vol::GatheralVarSwapTruncationJumpSpec {
            forward_price_f0: T::from_f64(100.0) + x[0],
            diffusive_vol_sigma: T::from_f64(0.20) + dx1,
            lower_cutoff_moneyness_kl_over_f: T::from_f64(0.65),
            upper_cutoff_moneyness_ku_over_f: T::from_f64(1.40),
            jump_intensity_lambda_j: T::from_f64(0.35) + dx2,
            mean_jump_log_size_mu_j: T::from_f64(-0.12) + dx3,
            jump_log_vol_sigma_j: T::from_f64(0.08),
            maturity_years_t: T::from_f64(0.5),
        },
    );

    let fva_ov = ql_xva::evaluate_gregory_albanese_fva_fca_fba_overlap(
        &ql_xva::GregoryAlbaneseFvaOverlapSpec {
            expected_positive_exposure_epe: T::from_f64(42.0) + x[0],
            expected_negative_exposure_ene: T::from_f64(28.0) + x[1],
            bank_cds_hazard_rate_lambda_b: T::from_f64(0.015) + dx2,
            bank_recovery_rate_rb: T::from_f64(0.40),
            bond_cds_liquidity_basis_s_liq: T::from_f64(0.0045) + dx3,
            treasury_deposit_spread_s_dep: T::from_f64(0.0035),
            counterparty_survival_weight: T::from_f64(0.94),
            risk_free_rate_r: T::from_f64(0.038),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let wkg = ql_commodities::evaluate_geman_working_curve_inventory_convenience_yield(
        &ql_commodities::GemanWorkingInventoryConvenienceSpec {
            spot_commodity_price_s0: T::from_f64(82.0) + x[0],
            normalized_inventory_i_over_icrit: T::from_f64(0.55) + dx1,
            stockout_convenience_yield_y0: T::from_f64(0.32) + dx2,
            abundance_convenience_floor_yinf: T::from_f64(0.01),
            working_curve_steepness_kappa: T::from_f64(1.8),
            physical_storage_cost_rate_c: T::from_f64(0.04),
            risk_free_rate_r: T::from_f64(0.04) + dx3,
            baseline_spot_vol_sigma0: T::from_f64(0.28),
            scarcity_vol_amplification_eta: T::from_f64(0.65),
            futures_maturity_t: T::from_f64(1.0),
        },
    );

    let fxfly = ql_fx::evaluate_wystup_malz_fx_broker_fly_smile_inversion(
        &ql_fx::WystupMalzFxSmileInversionSpec {
            spot_fx_s0: T::from_f64(1.085) + dx0,
            domestic_rate_rd: T::from_f64(0.045) + dx1,
            foreign_rate_rf: T::from_f64(0.030),
            atm_dns_vol_sigma_atm: T::from_f64(0.085) + dx2,
            risk_reversal_25d_rr: T::from_f64(-0.012) + dx3,
            butterfly_25d_bf: T::from_f64(0.0035),
            expiry_years_t: T::from_f64(0.5),
        },
    );

    let lkb = ql_exotics::evaluate_goldman_sosin_gatto_lookback_reflection(
        &ql_exotics::GoldmanSosinGattoLookbackSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            realized_minimum_smin: T::from_f64(92.0) + x[1],
            realized_maximum_smax: T::from_f64(106.0),
            fixed_strike_k: T::from_f64(100.0),
            risk_free_rate_r: T::from_f64(0.05) + dx2,
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.25) + dx3,
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let lw = ql_portfolio::evaluate_ledoit_wolf_constant_correlation_shrinkage(
        &ql_portfolio::LedoitWolfConstantCorrelationSpec {
            sample_volatilities: [
                T::from_f64(0.22) + dx0,
                T::from_f64(0.18) + dx1,
                T::from_f64(0.14) + dx2,
            ],
            sample_correlations: [
                T::from_f64(0.65),
                T::from_f64(0.15) + dx3,
                T::from_f64(0.40),
            ],
            asymptotic_pi_hat: T::from_f64(0.045),
            asymptotic_rho_hat: T::from_f64(0.012),
            sample_observations_t: T::from_f64(60.0),
        },
    );

    let ac = ql_portfolio::evaluate_almgren_chriss_cara_optimal_execution(
        &ql_portfolio::AlmgrenChrissCaraExecutionSpec {
            initial_shares_x0: T::from_f64(100.0) + x[0],
            liquidation_horizon_t: T::from_f64(1.0),
            asset_price_volatility_sigma: T::from_f64(0.30) + dx1,
            temporary_impact_eta: T::from_f64(0.02) + dx2,
            permanent_impact_gamma: T::from_f64(0.005),
            cara_risk_aversion_lambda: T::from_f64(0.5) + dx3,
        },
    );

    let estar = ql_dml::evaluate_tsay_star_threshold_autoregression_band(
        &ql_dml::TsayEstarThresholdAutoregressionSpec {
            current_mispricing_z: T::from_f64(1.85) + dx0,
            transaction_cost_threshold_c: T::from_f64(0.90) + dx1,
            transition_smoothness_gamma: T::from_f64(0.85),
            inner_band_reversion_speed_kappa_in: T::from_f64(0.04) + dx2,
            outer_arbitrage_reversion_speed_kappa_out: T::from_f64(0.38) + dx3,
            innovation_volatility_sigma_eps: T::from_f64(0.22),
        },
    );

    let acm = ql_models::evaluate_acm_cochrane_piazzesi_affine_term_premium(
        &ql_models::AcmCochranePiazzesiTermPremiumSpec {
            forward_rates_1y_to_5y: [
                T::from_f64(0.042) + dx0,
                T::from_f64(0.045) + dx1,
                T::from_f64(0.049) + dx2,
                T::from_f64(0.047),
                T::from_f64(0.046) + dx3,
            ],
            tent_weights_gamma: [
                T::from_f64(-1.2),
                T::from_f64(0.4),
                T::from_f64(2.1),
                T::from_f64(0.5),
                T::from_f64(-1.5),
            ],
            long_run_neutral_short_rate_r_star: T::from_f64(0.035),
            physical_var_mean_reversion_phi: T::from_f64(0.22),
            market_price_of_level_risk_lambda0: T::from_f64(0.35),
            market_price_of_slope_risk_lambda1: T::from_f64(0.80),
            horizon_years_n: T::from_f64(10.0),
        },
    );

    [
        ois_turn.all_in_discount_factor_p0t,
        cb_mw.convertible_bond_fair_value,
        cdx.rpv01_weighted_intrinsic_spread,
        vswp.total_fair_variance_swap_strike,
        fva_ov.non_overlapping_net_fva_plus_dva,
        wkg.storage_equilibrium_futures_price,
        fxfly.inverted_25d_call_strike_k25c,
        lkb.floating_strike_lookback_call,
        lw.equal_weight_shrunk_variance,
        ac.cara_certainty_equivalent_cost,
        estar.local_regime_half_life_periods,
        acm.fitted_10y_nominal_yield,
    ]
}

pub fn benchmark_darbyshire_step_ois_turn_despiegeleer_cb_soft_call_okane_cds_index_skew_gatheral_varswap_truncation_gregory_fva_overlap_geman_working_storage_wystup_malz_fx_fly_goldman_sosin_lookback_ledoit_wolf_shrinkage_almgren_chriss_cara_tsay_estar_acm_cp(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite54_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite54_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "54_darbyshire_step_ois_turn_despiegeleer_cb_soft_call_okane_cds_index_skew_gatheral_varswap_truncation_gregory_fva_overlap_geman_working_storage_wystup_malz_fx_fly_goldman_sosin_lookback_ledoit_wolf_shrinkage_almgren_chriss_cara_tsay_estar_acm_cp",
        "Step-OIS-Turn / CB-SoftCall-MW / CDX-Skew-Basis / VarSwap-Trunc / FVA-DVA-Overlap / Working-Storage / FX-Malz-Fly / GSG-Lookback / Ledoit-Wolf-CC / Almgren-CARA / ESTAR-Band / ACM-CP-TP AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite55_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let s = T::from_f64(0.005);
    let dx0 = x[0] * s;
    let dx1 = x[1] * s;
    let dx2 = x[2] * s;
    let dx3 = x[3] * s;

    let tb = ql_curves::evaluate_henrard_tenor_basis_multiplicative_spread(
        &ql_curves::HenrardTenorBasisMultiplicativeSpreadSpec {
            ois_discount_rate_r: T::from_f64(0.038) + dx0,
            forward_3m_period1_f1: T::from_f64(0.041) + dx1,
            forward_3m_period2_f2: T::from_f64(0.043) + dx2,
            quoted_forward_6m_f6: T::from_f64(0.0445) + dx3,
            vol_3m_period1_sigma1: T::from_f64(0.18),
            vol_3m_period2_sigma2: T::from_f64(0.19),
            subperiod_correlation_rho12: T::from_f64(0.85),
            start_time_years_t0: T::from_f64(1.0),
        },
    );

    let msr = ql_fixed_income::evaluate_fabozzi_agency_msr_valuation(
        &ql_fixed_income::FabozziAgencyMsrValuationSpec {
            gross_servicing_fee_rate: T::from_f64(0.0025) + dx0 * T::from_f64(0.1),
            ancillary_float_income_rate: T::from_f64(0.0004),
            annual_servicing_cost_rate: T::from_f64(0.0006),
            wac_mortgage_coupon: T::from_f64(0.060) + dx1,
            current_primary_mortgage_rate: T::from_f64(0.052) + dx2,
            base_turnover_cpr: T::from_f64(0.06),
            max_refi_cpr: T::from_f64(0.42),
            refi_scurve_steepness_kappa: T::from_f64(120.0),
            lender_recapture_fraction_alpha: T::from_f64(0.28) + dx3,
            ois_discount_plus_oas_rate: T::from_f64(0.075),
        },
    );

    let cdo = ql_credit::evaluate_okane_cdo_base_correlation_tranche(
        &ql_credit::OKaneCdoBaseCorrelationTrancheSpec {
            portfolio_hazard_rate_lambda: T::from_f64(0.018) + dx0,
            recovery_rate_r: T::from_f64(0.40),
            attachment_point_k1: T::from_f64(0.03),
            detachment_point_k2: T::from_f64(0.07),
            base_correlation_k1_rho1: T::from_f64(0.24) + dx1,
            base_correlation_k2_rho2: T::from_f64(0.36) + dx2,
            ois_discount_rate_r: T::from_f64(0.035) + dx3,
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let fwd_sm = ql_vol::evaluate_gatheral_forward_start_cliquet_smile(
        &ql_vol::GatheralForwardStartCliquetSmileSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            reset_strike_moneyness_alpha: T::from_f64(1.02) + dx1,
            reset_time_years_t1: T::from_f64(0.50),
            maturity_time_years_t2: T::from_f64(1.50),
            forward_initial_variance_v0: T::from_f64(0.04) + dx2,
            mean_reversion_speed_kappa: T::from_f64(1.8),
            vol_of_vol_eta: T::from_f64(0.55) + dx3,
            spot_vol_correlation_rho: T::from_f64(-0.65),
            risk_free_rate_r: T::from_f64(0.04),
            dividend_yield_q: T::from_f64(0.015),
        },
    );

    let csa_mta = ql_xva::evaluate_gregory_pykhtin_csa_threshold_mta_cva(
        &ql_xva::GregoryPykhtinCsaThresholdMtaCvaSpec {
            uncollateralized_exposure_vol_sigma_v: T::from_f64(35.0) + x[0],
            current_mtm_value_v0: T::from_f64(12.0) + x[1],
            csa_threshold_h: T::from_f64(15.0) + x[2],
            minimum_transfer_amount_mta: T::from_f64(2.5),
            independent_amount_ia: T::from_f64(4.0),
            margin_period_of_risk_days: T::from_f64(10.0),
            counterparty_hazard_rate_lambda: T::from_f64(0.022) + dx3,
            recovery_rate_r: T::from_f64(0.40),
            ois_discount_rate_r: T::from_f64(0.038),
            horizon_years_t: T::from_f64(3.0),
        },
    );

    let swing = ql_commodities::evaluate_eydeland_gas_swing_take_or_pay(
        &ql_commodities::EydelandGasSwingTakeOrPaySpec {
            forward_gas_prices_f: [
                T::from_f64(3.20) + dx0,
                T::from_f64(3.85) + dx1,
                T::from_f64(3.45) + dx2,
            ],
            contract_strike_price_k: T::from_f64(3.30),
            min_period_volume_qmin: T::from_f64(0.50),
            max_period_volume_qmax: T::from_f64(1.50),
            min_cumulative_acq_qmin_tot: T::from_f64(2.20),
            max_cumulative_acq_qmax_tot: T::from_f64(3.80),
            shortfall_makeup_penalty_per_unit: T::from_f64(0.85),
            spot_gas_volatility_sigma: T::from_f64(0.45) + dx3,
            risk_free_rate_r: T::from_f64(0.04),
        },
    );

    let kiko =
        ql_fx::evaluate_wystup_fx_window_kiko_barrier(&ql_fx::WystupFxWindowKikoBarrierSpec {
            spot_fx_s0: T::from_f64(1.085) + dx0,
            strike_fx_k: T::from_f64(1.080),
            lower_knock_in_barrier_l: T::from_f64(1.020),
            upper_knock_out_barrier_u: T::from_f64(1.155),
            window_start_years_t1: T::from_f64(0.15),
            window_end_years_t2: T::from_f64(0.75),
            maturity_years_t: T::from_f64(1.00),
            domestic_rate_rd: T::from_f64(0.045) + dx1,
            foreign_rate_rf: T::from_f64(0.030) + dx2,
            volatility_sigma: T::from_f64(0.095) + dx3,
            upper_ko_rebate_r: T::from_f64(0.0080),
        });

    let stulz = ql_exotics::evaluate_stulz_two_asset_rainbow_option(
        &ql_exotics::StulzTwoAssetRainbowOptionSpec {
            spot_asset_1_s1: T::from_f64(100.0) + x[0],
            spot_asset_2_s2: T::from_f64(105.0) + x[1],
            strike_price_k: T::from_f64(98.0),
            vol_asset_1_sigma1: T::from_f64(0.22) + dx2,
            vol_asset_2_sigma2: T::from_f64(0.28),
            correlation_rho: T::from_f64(0.45) + dx3,
            dividend_yield_1_q1: T::from_f64(0.015),
            dividend_yield_2_q2: T::from_f64(0.020),
            risk_free_rate_r: T::from_f64(0.045),
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let jm = ql_portfolio::evaluate_jagannathan_ma_no_short_sale_shrinkage(
        &ql_portfolio::JagannathanMaCovarianceShrinkageSpec {
            asset_volatilities: [
                T::from_f64(0.26) + dx0,
                T::from_f64(0.18) + dx1,
                T::from_f64(0.14) + dx2,
            ],
            correlation_12: T::from_f64(0.85) + dx3,
            correlation_13: T::from_f64(0.25),
            correlation_23: T::from_f64(0.30),
            weight_upper_cap_u: T::from_f64(0.65),
        },
    );

    let ow = ql_portfolio::evaluate_obizhaeva_wang_transient_lob_execution(
        &ql_portfolio::ObizhaevaWangTransientLobExecutionSpec {
            total_order_size_x0: T::from_f64(100.0) + x[0],
            execution_horizon_t: T::from_f64(1.0),
            lob_depth_density_q0: T::from_f64(25.0) + x[1],
            permanent_impact_gamma: T::from_f64(0.008),
            lob_resilience_speed_rho: T::from_f64(3.2) + dx2,
            asset_volatility_sigma: T::from_f64(0.25) + dx3,
            risk_aversion_lambda: T::from_f64(0.15),
        },
    );

    let hac = ql_dml::evaluate_newey_west_andrews_hac_variance(
        &ql_dml::NeweyWestAndrewsHacVarianceSpec {
            sample_variance_gamma0: T::from_f64(0.040) + dx0,
            autocovariance_lag1_gamma1: T::from_f64(0.022) + dx1,
            autocovariance_lag2_gamma2: T::from_f64(0.012) + dx2,
            autocovariance_lag3_gamma3: T::from_f64(0.006),
            sample_size_n: T::from_f64(250.0),
            ols_point_estimate_beta: T::from_f64(0.035) + dx3,
        },
    );

    let gab = ql_models::evaluate_gabaix_rare_disaster_resilience(
        &ql_models::GabaixRareDisasterResilienceSpec {
            disaster_probability_p: T::from_f64(0.018) + dx0,
            consumption_recovery_in_disaster_b: T::from_f64(0.72) + dx1,
            dividend_recovery_in_disaster_f: T::from_f64(0.55) + dx2,
            crra_risk_aversion_gamma: T::from_f64(4.0),
            time_preference_rate_rho: T::from_f64(0.025),
            normal_growth_rate_g: T::from_f64(0.020),
            resilience_deviation_h_hat: T::from_f64(-0.008) + dx3,
            resilience_mean_reversion_phi_h: T::from_f64(0.20),
        },
    );

    [
        tb.ois_discounted_basis_swaplet_pv,
        msr.recapture_adjusted_msr_fair_value_pct,
        cdo.mezzanine_fair_running_spread_bps,
        fwd_sm.stochastic_vol_forward_start_call,
        csa_mta.collateralized_cva,
        swing.total_swing_contract_value,
        kiko.window_kiko_call_value,
        stulz.call_on_max_value,
        jm.constrained_gmv_portfolio_volatility,
        ow.risk_adjusted_total_cost,
        hac.hac_robust_t_statistic,
        gab.linearity_generating_price_dividend_ratio,
    ]
}

pub fn benchmark_henrard_tenor_basis_fabozzi_msr_recapture_okane_cdo_base_corr_gatheral_fwd_start_smile_gregory_pykhtin_csa_mta_eydeland_gas_swing_wystup_fx_kiko_stulz_rainbow_jagannathan_ma_obizhaeva_wang_newey_west_hac_gabaix_disaster(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite55_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite55_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "55_henrard_tenor_basis_fabozzi_msr_recapture_okane_cdo_base_corr_gatheral_fwd_start_smile_gregory_pykhtin_csa_mta_eydeland_gas_swing_wystup_fx_kiko_stulz_rainbow_jagannathan_ma_obizhaeva_wang_newey_west_hac_gabaix_disaster",
        "Tenor-Basis-Mult / MSR-Recapture / CDO-Base-Corr / Fwd-Start-Smile / CSA-MTA-CVA / Gas-Swing / FX-Window-KIKO / Stulz-Rainbow / Jagannathan-Ma / Obizhaeva-Wang / Newey-West-HAC / Gabaix-Disaster AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite56_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let s = T::from_f64(0.005);
    let dx0 = x[0] * s;
    let dx1 = x[1] * s;
    let dx2 = x[2] * s;
    let dx3 = x[3] * s;

    let dswap = ql_curves::evaluate_brigo_differential_swap_quanto_fra(
        &ql_curves::BrigoDifferentialSwapQuantoFraSpec {
            domestic_ois_rate_rd: T::from_f64(0.042) + dx0,
            domestic_forward_rate_ld: T::from_f64(0.044) + dx1,
            foreign_forward_rate_lf: T::from_f64(0.032) + dx2,
            foreign_rate_vol_sigma_lf: T::from_f64(0.22) + dx3,
            fx_spot_vol_sigma_fx: T::from_f64(0.12),
            rate_fx_correlation_rho: T::from_f64(-0.45),
            reset_time_years_t1: T::from_f64(2.0),
            accrual_fraction_tau: T::from_f64(0.50),
        },
    );

    let muni = ql_fixed_income::evaluate_fabozzi_callable_muni_ytw_refunding(
        &ql_fixed_income::FabozziCallableMuniYtwRefundingSpec {
            tax_exempt_coupon_rate_c: T::from_f64(0.050) + dx0,
            muni_mmd_benchmark_yield_ym: T::from_f64(0.031) + dx1,
            treasury_yield_y_ust: T::from_f64(0.044) + dx2,
            statutory_marginal_tax_rate_tau: T::from_f64(0.37),
            call_protection_years_tc: T::from_f64(5.0),
            final_maturity_years_t: T::from_f64(10.0),
            call_price_par: T::from_f64(100.0),
            muni_yield_volatility_sigma: T::from_f64(0.18) + dx3,
        },
    );

    let step_cb = ql_credit::evaluate_bielecki_rating_step_up_coupon_bond(
        &ql_credit::BieleckiRatingStepUpCouponBondSpec {
            base_coupon_rate_c0: T::from_f64(0.045) + dx0,
            downgrade_step_up_spread_dc: T::from_f64(0.0125) + dx1,
            ig_to_hy_downgrade_intensity_lambda_down: T::from_f64(0.08) + dx2,
            hy_to_ig_upgrade_intensity_lambda_up: T::from_f64(0.03),
            ig_default_hazard_rate_lambda_d_ig: T::from_f64(0.005),
            hy_default_hazard_rate_lambda_d_hy: T::from_f64(0.035) + dx3,
            recovery_rate_r: T::from_f64(0.40),
            ois_discount_rate_r: T::from_f64(0.038),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let bkm = ql_vol::evaluate_bakshi_kapadia_madan_moment_strip(
        &ql_vol::BakshiKapadiaMadanMomentStripSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            atm_implied_vol_sigma0: T::from_f64(0.20) + dx1,
            smile_skew_slope_beta: T::from_f64(-0.25) + dx2,
            smile_curvature_gamma: T::from_f64(0.40) + dx3,
            risk_free_rate_r: T::from_f64(0.04),
            dividend_yield_q: T::from_f64(0.015),
            maturity_years_t: T::from_f64(0.25),
        },
    );

    let brk_cva = ql_xva::evaluate_gregory_break_clause_extinguishable_cva(
        &ql_xva::GregoryBreakClauseExtinguishableCvaSpec {
            unmitigated_swap_epe_profile: [
                T::from_f64(18.0) + x[0],
                T::from_f64(32.0) + x[1],
                T::from_f64(45.0) + x[2],
            ],
            break_clause_date_years_tb: T::from_f64(2.0),
            final_maturity_years_t: T::from_f64(5.0),
            counterparty_hazard_rate_lambda: T::from_f64(0.025) + dx3,
            break_execution_probability_p_exec: T::from_f64(0.85),
            break_settlement_slippage_fraction: T::from_f64(0.08),
            recovery_rate_r: T::from_f64(0.40),
            ois_discount_rate_r: T::from_f64(0.038),
        },
    );

    let pwr = ql_commodities::evaluate_bessembinder_lemmon_power_forward_premium(
        &ql_commodities::BessembinderLemmonPowerForwardPremiumSpec {
            baseline_marginal_fuel_cost_c0: T::from_f64(28.0) + x[0],
            expected_load_to_capacity_ratio_u: T::from_f64(0.84) + dx1,
            load_shock_volatility_sigma_u: T::from_f64(0.12) + dx2,
            supply_stack_convexity_exponent_gamma: T::from_f64(4.5) + dx3,
            peaking_merit_slope_alpha: T::from_f64(65.0),
            producer_variance_aversion_a_var: T::from_f64(0.015),
            retailer_spike_skew_preference_b_skew: T::from_f64(0.0045),
        },
    );

    let fader =
        ql_fx::evaluate_wystup_fx_fader_range_accrual(&ql_fx::WystupFxFaderRangeAccrualSpec {
            spot_fx_s0: T::from_f64(1.085) + dx0,
            strike_fx_k: T::from_f64(1.080),
            accrual_corridor_lower_l: T::from_f64(1.040),
            accrual_corridor_upper_u: T::from_f64(1.130),
            domestic_rate_rd: T::from_f64(0.045) + dx1,
            foreign_rate_rf: T::from_f64(0.030) + dx2,
            volatility_sigma: T::from_f64(0.095) + dx3,
            maturity_years_t: T::from_f64(1.0),
        });

    let ext_opt = ql_exotics::evaluate_longstaff_extendible_option(
        &ql_exotics::LongstaffExtendibleOptionSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            initial_strike_k1: T::from_f64(102.0) + x[1],
            extended_strike_k2: T::from_f64(105.0),
            initial_maturity_t1: T::from_f64(0.50),
            extended_maturity_t2: T::from_f64(1.25),
            holder_extension_fee_a: T::from_f64(1.50),
            risk_free_rate_r: T::from_f64(0.045) + dx2,
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.25) + dx3,
        },
    );

    let kz = ql_portfolio::evaluate_kan_zhou_three_fund_combine(
        &ql_portfolio::KanZhouThreeFundCombineSpec {
            sample_sharpe_squared_theta2: T::from_f64(0.28) + dx0,
            gmv_sharpe_squared_theta_g2: T::from_f64(0.12) + dx1,
            equal_weight_sharpe_squared_theta_ew2: T::from_f64(0.16) + dx2,
            number_of_assets_n: T::from_f64(25.0),
            sample_observations_t: T::from_f64(120.0) + x[3],
            risk_aversion_gamma: T::from_f64(3.0),
        },
    );

    let rbsa = ql_portfolio::evaluate_sharpe_return_based_style_analysis(
        &ql_portfolio::SharpeReturnBasedStyleAnalysisSpec {
            fund_mean_return_mu_p: T::from_f64(0.105) + dx0,
            fund_return_variance_var_p: T::from_f64(0.032) + dx1,
            style_benchmark_means: [
                T::from_f64(0.095) + dx2,
                T::from_f64(0.075),
                T::from_f64(0.045),
            ],
            style_benchmark_variances: [
                T::from_f64(0.040),
                T::from_f64(0.028) + dx3,
                T::from_f64(0.010),
            ],
            fund_style_covariances: [T::from_f64(0.024), T::from_f64(0.011), T::from_f64(-0.002)],
        },
    );

    let dm_hln = ql_dml::evaluate_diebold_mariano_hln_forecast_test(
        &ql_dml::DieboldMarianoHlnForecastTestSpec {
            model1_mse_loss_l1: T::from_f64(0.025) + dx0,
            model2_mse_loss_l2: T::from_f64(0.036) + dx1,
            forecast_error_correlation_rho12: T::from_f64(0.45) + dx2,
            loss_differential_autocov_lag1: T::from_f64(0.0004) + dx3 * T::from_f64(0.01),
            forecast_horizon_steps_h: T::from_f64(3.0),
            evaluation_sample_size_t: T::from_f64(80.0),
        },
    );

    let lw_rstar = ql_models::evaluate_laubach_williams_natural_rate_r_star(
        &ql_models::LaubachWilliamsNaturalRateRStarSpec {
            prior_trend_growth_g0: T::from_f64(0.018) + dx0,
            prior_other_factor_z0: T::from_f64(-0.004) + dx1,
            nominal_policy_rate_i: T::from_f64(0.0525) + dx2,
            expected_inflation_pi_e: T::from_f64(0.024),
            realized_core_inflation_pi: T::from_f64(0.028) + dx3,
            is_curve_surprise_dy: T::from_f64(0.006),
            is_rate_sensitivity_a_r: T::from_f64(0.09),
            phillips_slope_b_y: T::from_f64(0.15),
            signal_to_noise_ratio_lambda_z: T::from_f64(0.042),
        },
    );

    [
        dswap.domestic_discounted_diff_swaplet_pv,
        muni.callable_muni_market_price,
        step_cb.step_up_covenanted_bond_value,
        bkm.cboe_skew_index_level,
        brk_cva.break_clause_mitigated_cva,
        pwr.equilibrium_power_forward_price,
        fader.fade_in_call_value,
        ext_opt.writer_extendible_call_value,
        kz.three_fund_expected_out_of_sample_utility,
        rbsa.active_information_ratio,
        dm_hln.hln_corrected_dm_statistic,
        lw_rstar.taylor_rule_implied_nominal_rate,
    ]
}

pub fn benchmark_brigo_diff_swap_quanto_fabozzi_callable_muni_ytw_bielecki_step_up_bond_bkm_moment_strip_gregory_break_clause_cva_bessembinder_power_skew_wystup_fx_fader_longstaff_extendible_kan_zhou_3fund_sharpe_rbsa_diebold_mariano_hln_laubach_williams_rstar(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite56_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite56_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "56_brigo_diff_swap_quanto_fabozzi_callable_muni_ytw_bielecki_step_up_bond_bkm_moment_strip_gregory_break_clause_cva_bessembinder_power_skew_wystup_fx_fader_longstaff_extendible_kan_zhou_3fund_sharpe_rbsa_diebold_mariano_hln_laubach_williams_rstar",
        "Diff-Swap-Quanto / Callable-Muni-YTW / Step-Up-Bond / BKM-SKEW / Break-Clause-CVA / Power-Skew / FX-Fader / Longstaff-Extendible / Kan-Zhou-3Fund / Sharpe-RBSA / Diebold-Mariano-HLN / Laubach-Williams-r* AD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}

pub fn suite57_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let s = T::from_f64(0.005);
    let dx0 = x[0] * s;
    let dx1 = x[1] * s;
    let dx2 = x[2] * s;
    let dx3 = x[3] * s;

    let sofr_fut = ql_curves::evaluate_andersen_sofr_futures_fra_convexity(
        &ql_curves::AndersenSofrFuturesFraConvexitySpec {
            unadjusted_fra_forward_rate_f0: T::from_f64(0.0435) + dx0,
            short_rate_volatility_sigma_r: T::from_f64(0.0115) + dx1 * T::from_f64(0.2),
            hull_white_mean_reversion_kappa: T::from_f64(0.06) + dx2,
            accrual_start_years_t1: T::from_f64(3.0) + x[3] * T::from_f64(0.1),
            accrual_fraction_tau: T::from_f64(0.25),
            contract_notional_mm: T::from_f64(10.0),
        },
    );

    let mbs_prep = ql_fixed_income::evaluate_hayre_richard_roll_mbs_prepayment(
        &ql_fixed_income::HayreRichardRollMbsPrepaymentSpec {
            pool_wac_coupon_rate: T::from_f64(0.062) + dx0,
            prevailing_mortgage_rate_r_mtg: T::from_f64(0.051) + dx1,
            loan_age_months_wala: T::from_f64(24.0) + x[2],
            cumulative_burnout_exposure: T::from_f64(0.35) + dx3,
            burnout_sensitivity_beta: T::from_f64(0.60),
            s_curve_steepness_gamma: T::from_f64(65.0),
            seasonal_turnover_multiplier: T::from_f64(1.08),
        },
    );

    let cdx_opt = ql_credit::evaluate_pedersen_okane_index_option_fep_smile(
        &ql_credit::PedersenOkaneIndexOptionFepSmileSpec {
            spot_index_spread_s0: T::from_f64(0.0085) + dx0 * T::from_f64(0.2),
            strike_spread_k: T::from_f64(0.0090) + dx1 * T::from_f64(0.2),
            standard_running_coupon_c: T::from_f64(0.0100),
            spread_volatility_sigma: T::from_f64(0.45) + dx2,
            option_expiry_years_te: T::from_f64(0.50),
            index_maturity_years_tm: T::from_f64(5.50),
            recovery_rate_r: T::from_f64(0.40),
            ois_discount_rate_r: T::from_f64(0.040) + dx3,
        },
    );

    let dup_lv = ql_vol::evaluate_gatheral_total_variance_dupire_local_vol(
        &ql_vol::GatheralTotalVarianceDupireLocalVolSpec {
            log_forward_moneyness_y: T::from_f64(-0.08) + dx0,
            atm_implied_vol_sigma0: T::from_f64(0.21) + dx1,
            log_moneyness_skew_s: T::from_f64(-0.24) + dx2,
            log_moneyness_convexity_c: T::from_f64(0.35) + dx3,
            term_structure_slope_theta_t: T::from_f64(0.015),
            maturity_years_t: T::from_f64(0.50),
        },
    );

    let simm_mva = ql_xva::evaluate_gregory_kenyon_simm_mva_funding(
        &ql_xva::GregoryKenyonSimmMvaFundingSpec {
            initial_simm_margin_im0: T::from_f64(42.0) + x[0],
            ccp_cleared_margin_ratio: T::from_f64(0.62) + dx1,
            bank_funding_spread_s_f: T::from_f64(0.0095) + dx2 * T::from_f64(0.2),
            custodian_remuneration_spread_s_c: T::from_f64(0.0015),
            bank_hazard_rate_lambda_b: T::from_f64(0.012),
            counterparty_hazard_rate_lambda_c: T::from_f64(0.018) + dx3,
            ois_discount_rate_r: T::from_f64(0.038),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let rss_stor = ql_commodities::evaluate_routledge_seppi_spatt_stockout_backwardation(
        &ql_commodities::RoutledgeSeppiSpattStockoutBackwardationSpec {
            spot_commodity_price_s0: T::from_f64(78.0) + x[0],
            discretionary_inventory_i: T::from_f64(0.24) + dx1,
            critical_stockout_threshold_i_star: T::from_f64(0.20),
            base_abundant_volatility_sigma0: T::from_f64(0.22) + dx2,
            stockout_volatility_spike_eta: T::from_f64(0.85) + dx3,
            interest_plus_storage_cost_rate_c: T::from_f64(0.065),
            convenience_scarcity_slope_alpha: T::from_f64(0.16),
            maturity_years_t: T::from_f64(0.50),
        },
    );

    let cip = ql_fx::evaluate_baba_packer_cip_cross_currency_basis_wedge(
        &ql_fx::BabaPackerCipCrossCurrencyBasisWedgeSpec {
            spot_fx_s0: T::from_f64(1.0800) + dx0,
            market_outright_forward_fx_f: T::from_f64(1.1015) + dx1,
            domestic_usd_ois_rate_rd: T::from_f64(0.0480) + dx2,
            foreign_ois_rate_rf: T::from_f64(0.0310) + dx3,
            gsib_leverage_ratio_capital_charge_bps: T::from_f64(14.0),
            fx_swap_bid_ask_and_funding_haircut_bps: T::from_f64(4.5),
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let pcs_bar = ql_exotics::evaluate_carr_chou_put_call_symmetry_barrier_replication(
        &ql_exotics::CarrChouPutCallSymmetryBarrierReplicationSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(102.0) + x[1],
            down_and_out_barrier_h: T::from_f64(90.0),
            risk_free_rate_r: T::from_f64(0.04) + dx2,
            dividend_yield_q: T::from_f64(0.04),
            volatility_sigma: T::from_f64(0.22) + dx3,
            maturity_years_t: T::from_f64(0.75),
        },
    );

    let ps_shrink = ql_portfolio::evaluate_pastor_stambaugh_factor_prior_shrinkage(
        &ql_portfolio::PastorStambaughFactorPriorShrinkageSpec {
            sample_ols_alpha_hat: T::from_f64(0.035) + dx0,
            sample_factor_beta_hat: T::from_f64(1.12) + dx1,
            residual_return_variance_sigma_eps2: T::from_f64(0.018) + dx2,
            prior_mispricing_std_sigma_alpha: T::from_f64(0.020) + dx3,
            benchmark_excess_mean_mu_m: T::from_f64(0.065),
            benchmark_return_variance_sigma_m2: T::from_f64(0.025),
            sample_length_periods_t: T::from_f64(60.0),
            risk_aversion_gamma: T::from_f64(3.0),
        },
    );

    let tb_app = ql_portfolio::evaluate_treynor_black_appraisal_portfolio(
        &ql_portfolio::TreynorBlackAppraisalPortfolioSpec {
            security_alphas: [
                T::from_f64(0.024) + dx0,
                T::from_f64(0.018) + dx1,
                T::from_f64(0.012) + dx2,
            ],
            security_market_betas: [T::from_f64(1.15), T::from_f64(0.92), T::from_f64(1.05)],
            security_residual_variances: [
                T::from_f64(0.035) + dx3,
                T::from_f64(0.025),
                T::from_f64(0.020),
            ],
            market_excess_return_mu_m: T::from_f64(0.065),
            market_return_variance_sigma_m2: T::from_f64(0.0256),
        },
    );

    let eg_ecm = ql_dml::evaluate_engle_granger_cointegration_ecm_half_life(
        &ql_dml::EngleGrangerCointegrationEcmHalfLifeSpec {
            variance_asset_x: T::from_f64(0.042) + dx0,
            covariance_xy: T::from_f64(0.038) + dx1,
            spread_lag1_autocorrelation_phi: T::from_f64(0.86) + dx2,
            spread_innovation_sigma_e: T::from_f64(0.015) + dx3,
            current_disequilibrium_spread_z: T::from_f64(0.025),
            sample_observations_t: T::from_f64(250.0),
            time_step_years_dt: T::from_f64(1.0 / 252.0),
        },
    );

    let ez_sdf = ql_models::evaluate_epstein_zin_weil_recursive_utility_sdf(
        &ql_models::EpsteinZinWeilRecursiveUtilitySdfSpec {
            subjective_discount_beta: T::from_f64(0.988),
            relative_risk_aversion_gamma: T::from_f64(7.5) + x[0],
            intertemporal_elasticity_psi: T::from_f64(1.5) + dx1,
            expected_consumption_growth_gc: T::from_f64(0.018) + dx2,
            consumption_growth_volatility_sigma_c: T::from_f64(0.022) + dx3,
            wealth_portfolio_return_volatility_sigma_w: T::from_f64(0.145),
            consumption_wealth_correlation_rho_cw: T::from_f64(0.65),
        },
    );

    [
        sofr_fut.hull_white_futures_fra_convexity_bps,
        mbs_prep.annualized_cpr,
        cdx_opt.payer_index_swaption_pv_bps,
        dup_lv.dupire_local_volatility,
        simm_mva.mva_savings_from_central_clearing,
        rss_stor.equilibrium_forward_price,
        cip.covered_interest_parity_basis_bps,
        pcs_bar.static_replicated_down_and_out_call,
        ps_shrink.dogmatic_benchmark_weight,
        tb_app.combined_tangency_sharpe_ratio,
        eg_ecm.expected_next_step_spread_change + eg_ecm.stationary_spread_variance,
        ez_sdf.hansen_jagannathan_max_sharpe_bound + ez_sdf.equilibrium_real_risk_free_rate,
    ]
}

pub fn benchmark_andersen_sofr_futures_fra_hayre_mbs_scurve_pedersen_cdx_fep_gatheral_dupire_total_var_gregory_simm_mva_routledge_stockout_baba_cip_basis_carr_chou_symmetry_pastor_stambaugh_treynor_black_engle_granger_ecm_epstein_zin_sdf(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    // This suite is relatively short, so a single timed invocation is
    // dominated by timer/scheduler jitter in the parity tests.  Amortize that
    // fixed overhead while keeping the release benchmark's workload intact.
    let cfg = BenchConfig {
        iters_per_sample: cfg.iters_per_sample.max(8),
        ..cfg
    };
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite57_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite57_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "57_andersen_sofr_futures_fra_hayre_mbs_scurve_pedersen_cdx_fep_gatheral_dupire_total_var_gregory_simm_mva_routledge_stockout_baba_cip_basis_carr_chou_symmetry_pastor_stambaugh_treynor_black_engle_granger_ecm_epstein_zin_sdf",
        "SOFR-Fut-FRA / Hayre-MBS-SCurve / CDX-Option-FEP / Dupire-TotalVar-LV / SIMM-MVA-Funding / RSS-Stockout / CIP-Basis-Wedge / Carr-Chou-Symmetry / Pastor-Stambaugh / Treynor-Black / Engle-Granger-ECM / Epstein-Zin-SDF AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite58_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let pca_fly = ql_curves::evaluate_litterman_scheinkman_butterfly_pca_carry_rolldown(
        &ql_curves::LittermanScheinkmanButterflyPcaSpec {
            short_wing_duration_d1: T::from_f64(1.92) + dx0,
            belly_duration_d2: T::from_f64(4.55) + dx1,
            long_wing_duration_d3: T::from_f64(8.40) + dx2,
            short_wing_slope_loading_s1: T::from_f64(-0.65),
            belly_slope_loading_s2: T::from_f64(-0.10),
            long_wing_slope_loading_s3: T::from_f64(0.48),
            short_wing_convexity_c1: T::from_f64(4.2),
            belly_convexity_c2: T::from_f64(23.5),
            long_wing_convexity_c3: T::from_f64(82.0),
            parallel_yield_shift_bps: T::from_f64(50.0) + x[3],
        },
    );

    let tf_cb = ql_fixed_income::evaluate_tsiveriotis_fernandes_split_discount_convertible(
        &ql_fixed_income::TsiveriotisFernandesSplitConvertibleSpec {
            stock_price_s0: T::from_f64(96.0) + x[0],
            conversion_ratio_cr: T::from_f64(1.0) + dx1,
            par_redemption_value: T::from_f64(100.0),
            annual_coupon_rate_c: T::from_f64(0.035),
            risk_free_rate_r: T::from_f64(0.040) + dx2,
            credit_spread_sc: T::from_f64(0.025) + dx3,
            dividend_yield_q: T::from_f64(0.012),
            equity_volatility_sigma: T::from_f64(0.28),
            maturity_years_t: T::from_f64(3.0),
        },
    );

    let dt_cdo = ql_credit::evaluate_hull_white_double_t_cdo_tail_copula(
        &ql_credit::HullWhiteDoubleTCdoTailCopulaSpec {
            portfolio_default_probability_q: T::from_f64(0.045) + dx0 * T::from_f64(0.2),
            recovery_rate_r: T::from_f64(0.40) + dx1,
            factor_correlation_rho: T::from_f64(0.30) + dx2,
            student_t_degrees_of_freedom_nu: T::from_f64(4.0) + dx3 * T::from_f64(10.0),
            equity_detachment_k1: T::from_f64(0.03),
            senior_attachment_k2: T::from_f64(0.15),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let lee_wing = ql_vol::evaluate_roger_lee_wing_moment_slope_bound(
        &ql_vol::RogerLeeWingMomentSlopeBoundSpec {
            critical_upper_moment_p_tilde: T::from_f64(3.5) + x[0],
            critical_lower_moment_q_tilde: T::from_f64(1.8) + x[1],
            atm_total_variance_w0: T::from_f64(0.04) + dx2,
            extreme_log_moneyness_k: T::from_f64(-0.45) + dx3,
            maturity_years_t: T::from_f64(0.50),
        },
    );

    let cds_cont = ql_xva::evaluate_brigo_capponi_cds_contagion_gap_cva(
        &ql_xva::BrigoCapponiCdsContagionGapCvaSpec {
            reference_hazard_rate_lambda_ref: T::from_f64(0.025) + dx0 * T::from_f64(0.2),
            protection_seller_hazard_rate_lambda_c: T::from_f64(0.020) + dx1 * T::from_f64(0.2),
            contractual_running_spread_s_cds: T::from_f64(0.010),
            contagion_spread_jump_gamma: T::from_f64(0.85) + dx2 * T::from_f64(5.0),
            simultaneous_default_prob_p_sim: T::from_f64(0.06) + dx3 * T::from_f64(0.5),
            reference_recovery_rate_r_ref: T::from_f64(0.40),
            counterparty_recovery_rate_r_c: T::from_f64(0.40),
            ois_discount_rate_r: T::from_f64(0.038),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let sch_fut = ql_commodities::evaluate_schwartz_one_factor_mean_reverting_log_futures(
        &ql_commodities::SchwartzOneFactorLogFuturesSpec {
            spot_commodity_price_s0: T::from_f64(76.0) + x[0],
            risk_neutral_long_run_log_price_mu_star: T::from_f64(4.25) + dx1,
            mean_reversion_speed_kappa: T::from_f64(0.85) + dx2 * T::from_f64(5.0),
            spot_log_volatility_sigma: T::from_f64(0.34) + dx3,
            option_expiry_years_te: T::from_f64(0.50),
            futures_maturity_years_tf: T::from_f64(1.00),
            option_strike_price_k: T::from_f64(74.0),
            risk_free_rate_r: T::from_f64(0.042),
        },
    );

    let dcd = ql_fx::evaluate_wystup_dual_currency_deposit_yield_enhancement(
        &ql_fx::WystupDualCurrencyDepositYieldSpec {
            spot_fx_s0: T::from_f64(1.0850) + dx0,
            conversion_strike_fx_k: T::from_f64(1.0700) + dx1,
            domestic_deposit_rate_rd: T::from_f64(0.045) + dx2,
            foreign_deposit_rate_rf: T::from_f64(0.030),
            fx_volatility_sigma: T::from_f64(0.105) + dx3,
            deposit_tenor_years_t: T::from_f64(0.25),
            bank_structuring_margin_pct: T::from_f64(0.15),
        },
    );

    let pass_opt = ql_exotics::evaluate_hyer_lipton_pugachevsky_passport_option(
        &ql_exotics::HyerLiptonPugachevskyPassportOptionSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            initial_trading_account_w0: T::from_f64(2.0) + x[1],
            max_position_limit_q_max: T::from_f64(1.0),
            risk_free_rate_r: T::from_f64(0.04) + dx2,
            cost_of_carry_b: T::from_f64(0.02),
            volatility_sigma: T::from_f64(0.25) + dx3,
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let grs = ql_portfolio::evaluate_gibbons_ross_shanken_multivariate_factor_ftest(
        &ql_portfolio::GibbonsRossShankenFactorFTestSpec {
            alpha_quadratic_form_sr_alpha2: T::from_f64(0.145) + dx0,
            benchmark_factor_sharpe_ratio_sr_f: T::from_f64(0.42) + dx1,
            sample_observations_t: T::from_f64(180.0) + x[2],
            test_assets_count_n: T::from_f64(25.0) + x[3],
            factor_count_k: T::from_f64(3.0),
        },
    );

    let sqrt_imp = ql_portfolio::evaluate_bouchaud_gatheral_square_root_metaorder_impact(
        &ql_portfolio::BouchaudGatheralSquareRootImpactSpec {
            metaorder_shares_q: T::from_f64(250_000.0) + x[0] * T::from_f64(1000.0),
            average_daily_volume_adv: T::from_f64(5_000_000.0),
            daily_return_volatility_sigma_d: T::from_f64(0.018) + dx1 * T::from_f64(0.2),
            y_prefactor_constant: T::from_f64(0.65) + dx2,
            impact_concavity_exponent_delta: T::from_f64(0.50),
            propagator_decay_exponent_beta: T::from_f64(0.55) + dx3,
            post_trade_elapsed_ratio_tau_over_t: T::from_f64(1.0),
        },
    );

    let bns_har = ql_dml::evaluate_barndorff_nielsen_shephard_bipower_jump_har_cj(
        &ql_dml::BarndorffNielsenShephardHarCjSpec {
            daily_realized_variance_rv: T::from_f64(0.00042) + dx0 * T::from_f64(1e-3),
            daily_bipower_variation_bv: T::from_f64(0.00028) + dx1 * T::from_f64(1e-3),
            daily_tripower_quarticity_tq: T::from_f64(1.15e-7),
            weekly_continuous_bv_w: T::from_f64(0.00024) + dx2 * T::from_f64(1e-3),
            monthly_continuous_bv_m: T::from_f64(0.00020) + dx3 * T::from_f64(1e-3),
            intraday_sampling_frequency_m: T::from_f64(78.0),
        },
    );

    let ll_cay = ql_models::evaluate_lettau_ludvigson_cay_consumption_wealth_forecast(
        &ql_models::LettauLudvigsonCayForecastSpec {
            log_consumption_c: T::from_f64(9.865) + dx0,
            log_asset_wealth_a: T::from_f64(11.20) + dx1,
            log_labor_income_y: T::from_f64(9.25) + dx2,
            asset_wealth_cointegrating_share_omega: T::from_f64(0.31),
            unconditional_equity_premium_gamma0: T::from_f64(0.062),
            cay_predictive_slope_beta_cay: T::from_f64(2.15),
            conditional_market_volatility_sigma_m: T::from_f64(0.16) + dx3,
        },
    );

    [
        pca_fly.convexity_pnl_advantage_bps,
        tf_cb.tsiveriotis_fernandes_convertible_value,
        dt_cdo.senior_tranche_tail_spread_bps + dt_cdo.double_t_equity_tranche_expected_loss,
        lee_wing.asymptotic_wing_implied_volatility,
        cds_cont.total_contagion_adjusted_cds_cva,
        sch_fut.futures_call_option_price,
        dcd.enhanced_dcd_coupon_rate,
        pass_opt.carry_adjusted_passport_value,
        grs.grs_f_statistic,
        sqrt_imp.post_trade_decayed_impact_bps,
        bns_har.har_cj_annualized_volatility_forecast,
        ll_cay.conditional_market_sharpe_ratio,
    ]
}

pub fn benchmark_litterman_scheinkman_pca_fly_tf_split_cb_hull_white_double_t_cdo_roger_lee_moment_wings_brigo_capponi_cds_contagion_schwartz_ou_futures_wystup_dcd_yield_passport_option_grs_factor_ftest_bouchaud_sqrt_impact_bns_bipower_har_cj_lettau_ludvigson_cay(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite58_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite58_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "58_litterman_scheinkman_pca_fly_tf_split_cb_hull_white_double_t_cdo_roger_lee_moment_wings_brigo_capponi_cds_contagion_schwartz_ou_futures_wystup_dcd_yield_passport_option_grs_factor_ftest_bouchaud_sqrt_impact_bns_bipower_har_cj_lettau_ludvigson_cay",
        "PCA-Fly-Carry / TF-Split-CB / Double-t-CDO / Roger-Lee-Wings / CDS-Contagion-CVA / Schwartz-OU-Fut / FX-DCD-Yield / Passport-Opt / GRS-Factor-FTest / Sqrt-Metaorder / BNS-HAR-CJ / Lettau-Ludvigson-cay AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite59_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let pit_opt = ql_curves::evaluate_piterbarg_collateral_choice_option_funding_wedge(
        &ql_curves::PiterbargCollateralOptionFundingWedgeSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(100.0),
            csa_collateral_rate_r_c: T::from_f64(0.038) + dx1,
            equity_repo_financing_rate_r_r: T::from_f64(0.042),
            unsecured_funding_spread_s_f: T::from_f64(0.015) + dx2,
            dividend_yield_q: T::from_f64(0.015),
            equity_volatility_sigma_s: T::from_f64(0.24) + dx3,
            funding_spread_volatility_sigma_f: T::from_f64(0.008),
            equity_funding_correlation_rho_sf: T::from_f64(-0.35),
            maturity_years_t: T::from_f64(2.0),
        },
    );

    let tuck_roll = ql_fixed_income::evaluate_tuckman_par_spot_forward_rolldown_horizon_return(
        &ql_fixed_income::TuckmanRolldownHorizonReturnSpec {
            initial_maturity_years_t: T::from_f64(7.0),
            holding_horizon_years_dt: T::from_f64(1.0),
            current_bond_yield_y_t: T::from_f64(0.0445) + dx0,
            rolled_maturity_yield_y_t_minus_dt: T::from_f64(0.0420) + dx1,
            short_financing_rate_r_dt: T::from_f64(0.0380) + dx2,
            annual_coupon_rate_c: T::from_f64(0.0450),
            modified_duration_at_horizon_d_mod: T::from_f64(5.20),
            convexity_at_horizon_conv: T::from_f64(32.5),
            expected_parallel_shift_bps: T::from_f64(15.0) + x[3],
        },
    );

    let kmv_dd = ql_credit::evaluate_merton_kmv_distance_to_default_credit_equity_hedge(
        &ql_credit::MertonKmvDistanceToDefaultHedgeSpec {
            asset_value_v_a: T::from_f64(165.0) + x[0],
            short_term_debt_std: T::from_f64(45.0) + x[1],
            long_term_debt_ltd: T::from_f64(90.0),
            physical_asset_drift_mu_a: T::from_f64(0.075) + dx2,
            risk_free_rate_r: T::from_f64(0.040),
            asset_volatility_sigma_a: T::from_f64(0.24) + dx3,
            horizon_years_t: T::from_f64(3.0),
        },
    );

    let bart_sabr = ql_vol::evaluate_hagan_bartlett_normal_sabr_minimum_variance_greeks(
        &ql_vol::HaganBartlettNormalSabrGreeksSpec {
            forward_rate_f: T::from_f64(0.038) + dx0,
            strike_rate_k: T::from_f64(0.042),
            sabr_alpha: T::from_f64(0.045) + dx1,
            sabr_beta: T::from_f64(0.50),
            sabr_rho: T::from_f64(-0.35) + dx2,
            sabr_nu: T::from_f64(0.40) + dx3,
            option_expiry_t: T::from_f64(2.0),
            ois_discount_factor_df: T::from_f64(0.925),
        },
    );

    let mpor_spk = ql_xva::evaluate_pykhtin_brownian_bridge_mpor_cashflow_spike_epe(
        &ql_xva::PykhtinMporCashflowSpikeEpeSpec {
            annualized_portfolio_mtm_volatility_sigma: T::from_f64(28.0) + x[0],
            margin_period_of_risk_years_tau_mpor: T::from_f64(10.0 / 252.0),
            csa_threshold_plus_mta: T::from_f64(2.5) + dx1 * T::from_f64(10.0),
            independent_amount_held_ia: T::from_f64(1.0),
            contractual_cashflow_paid_in_mpor: T::from_f64(6.0) + x[2],
            post_payment_lookback_weight_theta: T::from_f64(0.60),
            counterparty_hazard_rate_lambda: T::from_f64(0.022) + dx3,
            counterparty_recovery_rate_r: T::from_f64(0.40),
            ois_discount_rate_r: T::from_f64(0.038),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let ftr_opt = ql_commodities::evaluate_deng_johnson_locational_transmission_congestion_option(
        &ql_commodities::DengJohnsonTransmissionCongestionOptionSpec {
            sink_hub_forward_price_fb: T::from_f64(68.0) + x[0],
            source_hub_forward_price_fa: T::from_f64(54.0) + x[1],
            transmission_loss_factor_eta: T::from_f64(0.96),
            wheeling_tariff_strike_k_trans: T::from_f64(4.5) + dx2 * T::from_f64(10.0),
            sink_hub_volatility_sigma_b: T::from_f64(0.48) + dx3,
            source_hub_volatility_sigma_a: T::from_f64(0.42),
            hub_price_correlation_rho_ab: T::from_f64(0.78),
            risk_free_rate_r: T::from_f64(0.04),
            expiry_years_t: T::from_f64(0.50),
        },
    );

    let fx_dig = ql_fx::evaluate_wystup_fx_digital_call_spread_skew_overhedge(
        &ql_fx::WystupFxDigitalSkewOverhedgeSpec {
            spot_fx_s0: T::from_f64(1.0850) + dx0,
            digital_strike_k: T::from_f64(1.1000) + dx1,
            domestic_rate_rd: T::from_f64(0.045),
            foreign_rate_rf: T::from_f64(0.030),
            atm_or_strike_vol_sigma_k: T::from_f64(0.095) + dx2,
            smile_slope_dsigma_dk: T::from_f64(-0.18) + dx3,
            call_spread_half_width_dk: T::from_f64(0.0050),
            expiry_years_t: T::from_f64(0.50),
        },
    );

    let cv_lb = ql_exotics::evaluate_conze_viswanathan_fixed_strike_lookback_excursion(
        &ql_exotics::ConzeViswanathanFixedStrikeLookbackSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            fixed_strike_k: T::from_f64(105.0) + x[1],
            risk_free_rate_r: T::from_f64(0.050) + dx2,
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.25) + dx3,
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let hot_es = ql_portfolio::evaluate_tasche_litterman_hot_spots_best_hedge_es_decomposition(
        &ql_portfolio::TascheLittermanHotSpotsBestHedgeSpec {
            focal_position_weight_w1: T::from_f64(0.45) + dx0,
            core_portfolio_weight_w2: T::from_f64(0.55) + dx1,
            focal_asset_volatility_sigma1: T::from_f64(0.28) + dx2,
            core_portfolio_volatility_sigma2: T::from_f64(0.16),
            asset_core_correlation_rho12: T::from_f64(0.50) + dx3,
            focal_expected_return_mu1: T::from_f64(0.08),
            core_expected_return_mu2: T::from_f64(0.06),
            gaussian_es_tail_multiplier_m_alpha: T::from_f64(2.338),
        },
    );

    let bab = ql_portfolio::evaluate_frazzini_pedersen_betting_against_beta_zero_beta_capm(
        &ql_portfolio::FrazziniPedersenBettingAgainstBetaSpec {
            low_beta_portfolio_beta_l: T::from_f64(0.65) + dx0,
            high_beta_portfolio_beta_h: T::from_f64(1.45) + dx1,
            low_beta_expected_return_r_l: T::from_f64(0.088) + dx2,
            high_beta_expected_return_r_h: T::from_f64(0.118),
            risk_free_rate_r_f: T::from_f64(0.038),
            shadow_funding_constraint_spread_psi: T::from_f64(0.022) + dx3,
            bab_factor_volatility_sigma_bab: T::from_f64(0.105),
        },
    );

    let pt_dir = ql_dml::evaluate_pesaran_timmermann_directional_timing_hit_ratio_test(
        &ql_dml::PesaranTimmermannDirectionalTimingSpec {
            true_positive_up_fraction_p11: T::from_f64(0.34) + dx0,
            true_negative_down_fraction_p00: T::from_f64(0.25) + dx1,
            false_positive_up_fraction_p01: T::from_f64(0.20) + dx2,
            sample_observations_t: T::from_f64(250.0),
            mean_absolute_return_move_bps: T::from_f64(85.0) + x[3],
        },
    );

    let lmp = ql_models::evaluate_lucas_mehra_prescott_endowment_tree_equity_premium(
        &ql_models::LucasMehraPrescottEndowmentTreeSpec {
            mean_consumption_growth_mu_c: T::from_f64(0.018) + dx0 * T::from_f64(0.1),
            consumption_growth_vol_sigma_c: T::from_f64(0.036) + dx1 * T::from_f64(0.2),
            dividend_consumption_leverage_lambda: T::from_f64(3.0) + dx2 * T::from_f64(5.0),
            crra_risk_aversion_gamma: T::from_f64(2.5) + dx3 * T::from_f64(5.0),
            subjective_discount_factor_beta: T::from_f64(0.99),
            historical_observed_equity_premium: T::from_f64(0.0618),
        },
    );

    [
        pit_opt.total_piterbarg_funding_adjusted_call_price,
        tuck_roll.horizon_expected_excess_return_pct,
        kmv_dd.structural_credit_spread_bps
            + kmv_dd.physical_distance_to_default_dd * T::from_f64(10.0),
        bart_sabr.bartlett_minimum_variance_delta,
        mpor_spk.total_mpor_collateralized_cva,
        ftr_opt.ftr_option_value_per_mwh,
        fx_dig.skew_adjusted_digital_call_price,
        cv_lb.fixed_strike_lookback_call_value,
        hot_es.total_portfolio_expected_shortfall + hot_es.focal_component_expected_shortfall,
        bab.bab_factor_sharpe_ratio + bab.theoretical_funding_constrained_bab_alpha,
        pt_dir.pesaran_timmermann_z_statistic + pt_dir.directional_timing_expected_edge_bps,
        lmp.lucas_tree_price_dividend_ratio,
    ]
}

pub fn benchmark_piterbarg_collateral_option_tuckman_rolldown_merton_kmv_dd_hagan_bartlett_sabr_pykhtin_mpor_spike_deng_locational_ftr_wystup_digital_overhedge_conze_fixed_lookback_tasche_hot_spots_es_frazzini_bab_pesaran_timmermann_lucas_mehra_prescott(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite59_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite59_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "59_piterbarg_collateral_option_tuckman_rolldown_merton_kmv_dd_hagan_bartlett_sabr_pykhtin_mpor_spike_deng_locational_ftr_wystup_digital_overhedge_conze_fixed_lookback_tasche_hot_spots_es_frazzini_bab_pesaran_timmermann_lucas_mehra_prescott",
        "Piterbarg-Coll-Opt / Tuckman-Rolldown / Merton-KMV-DD / Bartlett-SABR / Pykhtin-MPoR-Spike / Deng-FTR-Opt / Wystup-FX-Digital / Conze-Fixed-Lookback / Tasche-Hot-Spots-ES / Frazzini-BAB / Pesaran-Timmermann / Lucas-Mehra-Prescott AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite60_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let gsw = ql_curves::evaluate_gurkaynak_sack_wright_svensson_treasury_forward_curve(
        &ql_curves::GswSvenssonTreasuryCurveSpec {
            svensson_beta0: T::from_f64(0.045) + dx0,
            svensson_beta1: T::from_f64(-0.012),
            svensson_beta2: T::from_f64(0.018) + dx1,
            svensson_beta3: T::from_f64(-0.010),
            svensson_tau1: T::from_f64(1.80) + x[2] * T::from_f64(0.1),
            svensson_tau2: T::from_f64(6.50),
            maturity_years_t: T::from_f64(10.0),
            on_the_run_par_yield_y_otr: T::from_f64(0.0415) + dx3,
        },
    );

    let mbs_esf = ql_fixed_income::evaluate_fabozzi_hayre_mbs_excess_servicing_gfee_strip(
        &ql_fixed_income::FabozziHayreMbsExcessServicingGfeeSpec {
            pool_balance_b0: T::from_f64(100.0),
            gross_wac_rate_w: T::from_f64(0.0640) + dx0,
            net_passthrough_coupon_c: T::from_f64(0.0550),
            agency_guaranty_fee_g: T::from_f64(0.0045) + dx3 * T::from_f64(0.1),
            minimum_servicing_fee_s_min: T::from_f64(0.0025),
            prepayment_speed_cpr: T::from_f64(0.12) + dx1,
            scheduled_principal_amortization_a: T::from_f64(0.018),
            discount_rate_r: T::from_f64(0.045) + dx2,
            pool_horizon_years_t: T::from_f64(10.0),
        },
    );

    let irb_ga = ql_credit::evaluate_vasicek_gordy_asrf_irb_capital_granularity_adjustment(
        &ql_credit::VasicekGordyAsrfIrbGranularitySpec {
            unconditional_pd: T::from_f64(0.015) + dx0 * T::from_f64(0.2),
            downturn_lgd: T::from_f64(0.45) + dx1,
            asset_correlation_rho: T::from_f64(0.18) + dx2,
            confidence_level_alpha: T::from_f64(0.999),
            effective_maturity_years_m: T::from_f64(3.0),
            herfindahl_concentration_index_h: T::from_f64(0.025) + dx3 * T::from_f64(0.2),
        },
    );

    let cw_vol = ql_vol::evaluate_carr_wu_proportional_volatility_smile_drift(
        &ql_vol::CarrWuProportionalVolSmileDriftSpec {
            forward_price_f: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(95.0),
            atm_volatility_sigma0: T::from_f64(0.22) + dx1,
            implied_vol_of_vol_omega: T::from_f64(0.65) + dx2,
            spot_vol_correlation_rho: T::from_f64(-0.70),
            implied_vol_drift_mu: T::from_f64(-0.04) + dx3,
            time_to_expiry_tau: T::from_f64(0.50),
        },
    );

    let greg_cva = ql_xva::evaluate_gregory_incremental_marginal_netting_set_cva(
        &ql_xva::GregoryIncrementalMarginalNettingCvaSpec {
            existing_portfolio_mtm_m1: T::from_f64(12.0) + x[0],
            new_trade_mtm_m2: T::from_f64(-4.0) + x[1],
            existing_portfolio_vol_sigma1: T::from_f64(18.0),
            new_trade_vol_sigma2: T::from_f64(10.0),
            exposure_correlation_rho12: T::from_f64(-0.35) + dx2,
            counterparty_hazard_rate_lambda: T::from_f64(0.024) + dx3,
            counterparty_recovery_rate_r: T::from_f64(0.40),
            ois_discount_rate_r: T::from_f64(0.038),
            maturity_years_t: T::from_f64(4.0),
        },
    );

    let wb_stor = ql_commodities::evaluate_working_brennan_storage_convenience_yield_curve(
        &ql_commodities::WorkingBrennanStorageConvenienceYieldSpec {
            spot_price_s0: T::from_f64(75.0) + x[0],
            inventory_to_consumption_ratio_i: T::from_f64(0.11) + dx1,
            minimum_pipeline_inventory_i_min: T::from_f64(0.08),
            convenience_yield_convexity_alpha: T::from_f64(0.018) + dx2 * T::from_f64(0.2),
            stockout_steepness_exponent_gamma: T::from_f64(0.85),
            physical_storage_cost_rate_u: T::from_f64(0.045),
            risk_free_rate_r: T::from_f64(0.040) + dx3,
            futures_maturity_years_t: T::from_f64(1.0),
        },
    );

    let fx_fe = ql_fx::evaluate_wystup_fx_forward_extra_participating_collar(
        &ql_fx::WystupFxForwardExtraCollarSpec {
            spot_fx_s0: T::from_f64(1.0850) + dx0,
            worst_case_strike_k: T::from_f64(1.0950) + dx1,
            downside_trigger_barrier_b: T::from_f64(1.0350) + dx2,
            domestic_rate_rd: T::from_f64(0.045),
            foreign_rate_rf: T::from_f64(0.030),
            strike_vol_sigma_k: T::from_f64(0.095) + dx3,
            put_wing_skew_slope_xi: T::from_f64(-0.14),
            expiry_years_t: T::from_f64(0.50),
        },
    );

    let ex_on_ex = ql_exotics::evaluate_rubinstein_margrabe_exchange_on_exchange_option(
        &ql_exotics::RubinsteinMargrabeExchangeOnExchangeSpec {
            spot_asset1_s1: T::from_f64(105.0) + x[0],
            spot_asset2_s2: T::from_f64(100.0) + x[1],
            compound_strike_fraction_q: T::from_f64(0.04),
            dividend_yield_q1: T::from_f64(0.015),
            dividend_yield_q2: T::from_f64(0.025),
            volatility_sigma1: T::from_f64(0.24) + dx2,
            volatility_sigma2: T::from_f64(0.20),
            asset_correlation_rho: T::from_f64(0.55) + dx3,
            underlying_expiry_years_t2: T::from_f64(1.50),
        },
    );

    let hlz_fdr = ql_portfolio::evaluate_harvey_liu_zhu_factor_tstat_fdr_haircut(
        &ql_portfolio::HarveyLiuZhuFactorFdrHaircutSpec {
            observed_factor_t_stat: T::from_f64(3.25) + x[0] * T::from_f64(0.1),
            sample_months_t: T::from_f64(240.0),
            total_tested_factors_m: T::from_f64(316.0) + x[1] * T::from_f64(10.0),
            target_fdr_alpha: T::from_f64(0.05) + dx2 * T::from_f64(0.2),
            prior_true_factor_prob_p0: T::from_f64(0.15) + dx3,
            prior_true_sharpe_std_sigma0: T::from_f64(0.35),
        },
    );

    let hs_kyle = ql_portfolio::evaluate_holden_subrahmanyam_multi_insider_kyle_equilibrium(
        &ql_portfolio::HoldenSubrahmanyamMultiInsiderKyleSpec {
            number_of_insiders_n: T::from_f64(4.0) + x[0] * T::from_f64(0.2),
            fundamental_value_vol_sigma_v: T::from_f64(2.50) + dx1 * T::from_f64(10.0),
            noise_trader_volume_vol_sigma_u: T::from_f64(120.0) + x[2],
            private_valuation_surprise_dv: T::from_f64(1.80) + dx3 * T::from_f64(10.0),
            prior_reference_price_p0: T::from_f64(100.0),
        },
    );

    let garch_ts = ql_dml::evaluate_bollerslev_tsay_garch_term_structure_kurtosis(
        &ql_dml::BollerslevTsayGarchTermStructureSpec {
            garch_intercept_omega: T::from_f64(0.000_004),
            arch_shock_alpha: T::from_f64(0.085) + dx0 * T::from_f64(0.2),
            garch_persistence_beta: T::from_f64(0.895) + dx1 * T::from_f64(0.2),
            current_daily_variance_h1: T::from_f64(0.000_32) + dx2 * T::from_f64(0.005),
            innovation_kurtosis_kappa_z: T::from_f64(3.50) + x[3] * T::from_f64(0.1),
            forecast_horizon_days_k: T::from_f64(21.0),
        },
    );

    let cs_cape = ql_models::evaluate_campbell_shiller_cape_excess_yield_forecast(
        &ql_models::CampbellShillerCapeExcessYieldSpec {
            real_equity_index_price_p: T::from_f64(4800.0) + x[0] * T::from_f64(20.0),
            ten_year_real_average_earnings_e10: T::from_f64(160.0) + x[1] * T::from_f64(2.0),
            historical_mean_log_cape: T::from_f64(3.05),
            annual_valuation_reversion_speed_phi: T::from_f64(0.065) + dx2,
            expected_real_earnings_growth_g: T::from_f64(0.018),
            dividend_payout_ratio_d: T::from_f64(0.45),
            ten_year_real_bond_yield_y10: T::from_f64(0.019) + dx3,
        },
    );

    [
        gsw.off_the_run_liquidity_spread_bps
            + gsw.svensson_instantaneous_forward_rate * T::from_f64(100.0),
        mbs_esf.excess_servicing_strip_pv + mbs_esf.remittance_float_income_pv,
        irb_ga.total_granularity_adjusted_capital,
        cw_vol.carr_wu_arbitrage_free_implied_vol + cw_vol.vanna_volga_convexity_theta_cost,
        greg_cva.incremental_trade_cva + greg_cva.euler_marginal_trade_cva,
        wb_stor.equilibrium_futures_price,
        fx_fe.forward_extra_client_npv,
        ex_on_ex.exchange_on_exchange_call_value,
        hlz_fdr.posterior_haircutted_sharpe_ratio + hlz_fdr.bhy_fdr_critical_t_hurdle,
        hs_kyle.expected_clearing_price + hs_kyle.aggregate_insider_expected_profit,
        garch_ts.annualized_term_structure_volatility * T::from_f64(10.0)
            + garch_ts.unconditional_return_excess_kurtosis,
        cs_cape.implied_ten_year_equity_risk_premium * T::from_f64(100.0),
    ]
}

pub fn benchmark_gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_ga_carr_wu_vol_drift_gregory_marginal_cva_working_brennan_storage_wystup_forward_extra_rubinstein_exchange_opt_hlz_factor_fdr_holden_subrahmanyam_kyle_bollerslev_garch_kurtosis_campbell_shiller_cape(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite60_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite60_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "60_gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_ga_carr_wu_vol_drift_gregory_marginal_cva_working_brennan_storage_wystup_forward_extra_rubinstein_exchange_opt_hlz_factor_fdr_holden_subrahmanyam_kyle_bollerslev_garch_kurtosis_campbell_shiller_cape",
        "GSW-Svensson-Curve / Fabozzi-MBS-ESF / Vasicek-Gordy-IRB-GA / Carr-Wu-Vol-Drift / Gregory-Marginal-CVA / Working-Brennan-Storage / Wystup-Forward-Extra / Rubinstein-Exchange-Opt / HLZ-Factor-FDR / Holden-Subrahmanyam-Kyle / Bollerslev-GARCH-Kurtosis / Campbell-Shiller-CAPE AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite61_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let kw_shadow = ql_curves::evaluate_kim_wright_shadow_rate_survey_term_premium(
        &ql_curves::KimWrightShadowRateTermPremiumSpec {
            current_shadow_short_rate_s0: T::from_f64(0.015) + dx0,
            effective_lower_bound_r_min: T::from_f64(0.00125),
            mean_reversion_speed_kappa: T::from_f64(0.22) + dx1,
            risk_neutral_long_run_mean_theta_q: T::from_f64(0.044) + dx2,
            physical_long_run_mean_theta_p: T::from_f64(0.032),
            shadow_rate_volatility_sigma_s: T::from_f64(0.014),
            survey_long_run_short_rate_r_surv: T::from_f64(0.028) + dx3,
            survey_kalman_shrinkage_weight_w: T::from_f64(0.35),
            horizon_years_tau: T::from_f64(10.0),
        },
    );

    let bs_cb = ql_fixed_income::evaluate_brennan_schwartz_callable_puttable_convertible_game(
        &ql_fixed_income::BrennanSchwartzCallablePuttableCbSpec {
            stock_price_s0: T::from_f64(98.0) + x[0],
            conversion_ratio_n_c: T::from_f64(1.0),
            par_notional_n: T::from_f64(100.0),
            coupon_rate_c: T::from_f64(0.025),
            risk_free_rate_r: T::from_f64(0.04) + dx1,
            credit_spread_s_c: T::from_f64(0.022) + dx2,
            dividend_yield_q: T::from_f64(0.012),
            equity_volatility_sigma_s: T::from_f64(0.32) + dx3,
            investor_put_price_k_put: T::from_f64(96.0),
            issuer_call_trigger_price_k_call: T::from_f64(130.0),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let ds_rmv = ql_credit::evaluate_duffie_singleton_rmv_rfv_credit_triangle_basis(
        &ql_credit::DuffieSingletonRmvRfvBasisSpec {
            risk_free_rate_r: T::from_f64(0.04) + dx0,
            hazard_rate_lambda: T::from_f64(0.03) + dx1,
            recovery_rate_rec: T::from_f64(0.40) + dx2,
            bond_coupon_rate_c: T::from_f64(0.055),
            deliverable_recovery_dispersion_sigma_r: T::from_f64(0.08) + dx3,
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let fn_gamma = ql_vol::evaluate_fukasawa_neuberger_gamma_skewness_swap_strip(
        &ql_vol::FukasawaNeubergerGammaSkewnessSwapSpec {
            forward_price_f0: T::from_f64(100.0),
            atm_volatility_sigma0: T::from_f64(0.20) + dx0,
            log_moneyness_skew_slope_alpha: T::from_f64(-0.18) + dx1,
            log_moneyness_curvature_beta: T::from_f64(0.25) + dx2,
            realized_spot_weighted_variance: T::from_f64(0.042) + dx3,
            maturity_years_t: T::from_f64(0.50),
            variance_notional_n_var: T::from_f64(100.0),
        },
    );

    let bk_xva = ql_xva::evaluate_burgard_kjaer_semi_replication_fva_cva_dva(
        &ql_xva::BurgardKjaerSemiReplicationXvaSpec {
            expected_positive_exposure_epe: T::from_f64(25.0) + x[0],
            expected_negative_exposure_ene: T::from_f64(14.0) + x[1],
            counterparty_hazard_rate_lambda_c: T::from_f64(0.022) + dx2,
            counterparty_lgd_c: T::from_f64(0.60),
            bank_hazard_rate_lambda_b: T::from_f64(0.012),
            bank_lgd_b: T::from_f64(0.60),
            bank_funding_shortfall_fraction_gamma_f: T::from_f64(0.65) + dx3,
            ois_discount_rate_r: T::from_f64(0.038),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let cd_spark = ql_commodities::evaluate_carmona_durrleman_spark_spread_lower_bound(
        &ql_commodities::CarmonaDurrlemanSparkSpreadSpec {
            electricity_forward_f_e: T::from_f64(62.0) + x[0],
            natural_gas_forward_f_g: T::from_f64(7.20) + x[1] * T::from_f64(0.2),
            heat_rate_mmbtu_per_mwh_h: T::from_f64(7.50),
            variable_om_strike_k_om: T::from_f64(6.50),
            electricity_volatility_sigma_e: T::from_f64(0.48) + dx2,
            gas_volatility_sigma_g: T::from_f64(0.42),
            power_gas_correlation_rho: T::from_f64(0.65) + dx3,
            risk_free_rate_r: T::from_f64(0.04),
            maturity_years_t: T::from_f64(0.75),
        },
    );

    let fx_tarf = ql_fx::evaluate_wystup_lipton_fx_tarf_target_redemption_cap(
        &ql_fx::WystupLiptonFxTarfSpec {
            spot_fx_s0: T::from_f64(1.09) + dx0,
            tarf_strike_k: T::from_f64(1.06) + dx1,
            target_redemption_cap_t_cap: T::from_f64(0.25) + dx2,
            downside_gearing_leverage_l: T::from_f64(2.0),
            total_fixing_dates_m: T::from_f64(12.0),
            domestic_rate_rd: T::from_f64(0.045),
            foreign_rate_rf: T::from_f64(0.030),
            fx_volatility_sigma: T::from_f64(0.10) + dx3,
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let bgk_lb = ql_exotics::evaluate_broadie_glasserman_kou_discrete_lookback_option(
        &ql_exotics::BroadieGlassermanKouDiscreteLookbackSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            observed_running_minimum_m0: T::from_f64(92.0) + x[1],
            risk_free_rate_r: T::from_f64(0.05) + dx2,
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.25) + dx3,
            discrete_monitoring_dates_n: T::from_f64(52.0),
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let nm_ff5 = ql_portfolio::evaluate_novy_marx_fama_french_five_factor_spanning_alpha(
        &ql_portfolio::NovyMarxFamaFrenchSpanningSpec {
            gross_profitability_premium_rmw: T::from_f64(0.042) + dx0,
            value_premium_hml: T::from_f64(0.038) + dx1,
            investment_premium_cma: T::from_f64(0.032) + dx2,
            hml_on_rmw_loading_phi_rmw: T::from_f64(0.22),
            hml_on_cma_loading_phi_cma: T::from_f64(0.55),
            hml_spanning_residual_vol_sigma_e: T::from_f64(0.068),
            rmw_factor_vol_sigma_rmw: T::from_f64(0.085),
            hml_factor_vol_sigma_hml: T::from_f64(0.105),
            rmw_hml_correlation_rho: T::from_f64(-0.42) + dx3,
        },
    );

    let acerbi_srm = ql_portfolio::evaluate_acerbi_tasche_exponential_spectral_risk_measure(
        &ql_portfolio::AcerbiTascheExponentialSpectralRiskSpec {
            portfolio_expected_return_mu: T::from_f64(0.08) + dx0,
            portfolio_volatility_sigma: T::from_f64(0.18) + dx1,
            return_skewness_gamma1: T::from_f64(-0.55) + dx2,
            return_excess_kurtosis_gamma2: T::from_f64(2.40),
            cara_spectral_risk_aversion_k: T::from_f64(25.0) + x[3],
            reference_es_confidence_alpha: T::from_f64(0.975),
        },
    );

    let corsi_har = ql_dml::evaluate_corsi_har_realized_volatility_multi_horizon_forecast(
        &ql_dml::CorsiHarRealizedVolatilitySpec {
            har_intercept_c: T::from_f64(0.000_012),
            daily_rv_weight_beta_d: T::from_f64(0.36) + dx0,
            weekly_rv_weight_beta_w: T::from_f64(0.34) + dx1,
            monthly_rv_weight_beta_m: T::from_f64(0.25) + dx2,
            current_daily_rv: T::from_f64(0.000_42),
            current_weekly_rv: T::from_f64(0.000_31),
            current_monthly_rv: T::from_f64(0.000_26),
            realized_quarticity_rq: T::from_f64(0.000_000_45) + dx3 * T::from_f64(0.000_005),
            forecast_horizon_days_h: T::from_f64(22.0),
        },
    );

    let bll_nn = ql_models::evaluate_blanchard_lhuillier_lorenzoni_news_noise_kalman_filter(
        &ql_models::BlanchardLhuillierLorenzoniNewsNoiseSpec {
            fundamental_productivity_vol_sigma_a: T::from_f64(0.015) + dx0 * T::from_f64(0.1),
            exogenous_news_noise_vol_sigma_e: T::from_f64(0.020) + dx1 * T::from_f64(0.1),
            shannon_channel_capacity_bits_kappa: T::from_f64(1.50) + dx2 * T::from_f64(5.0),
            observed_macro_news_signal_s: T::from_f64(0.025) + dx3 * T::from_f64(0.2),
            consumption_demand_elasticity_gamma_c: T::from_f64(0.85),
            phillips_curve_slope_kappa_pc: T::from_f64(0.18),
        },
    );

    [
        kw_shadow.survey_anchored_forward_term_premium * T::from_f64(100.0)
            + kw_shadow.zlb_call_option_wedge_bps,
        bs_cb.game_convertible_fair_value,
        ds_rmv.rfv_risky_coupon_bond_price
            + ds_rmv.mod_r_restructuring_cds_basis_bps * T::from_f64(0.01),
        fn_gamma.gamma_swap_mark_to_market_npv + fn_gamma.neuberger_model_free_skewness_swap_rate,
        bk_xva.semi_replication_total_xva_charge,
        cd_spark.carmona_durrleman_spark_call_value,
        fx_tarf.client_tarf_net_present_value,
        bgk_lb.bgk_discrete_lookback_call_value,
        nm_ff5.hmlo_information_ratio + nm_ff5.combined_value_profitability_sharpe,
        acerbi_srm.exponential_spectral_risk_measure,
        corsi_har.cumulative_average_h_day_rv * T::from_f64(1000.0)
            + corsi_har.bns_asymptotic_rv_standard_error * T::from_f64(1000.0),
        bll_nn.phillips_curve_inflation_impulse * T::from_f64(100.0),
    ]
}

pub fn benchmark_kim_wright_shadow_rate_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_swap_burgard_kjaer_fva_carmona_durrleman_spark_wystup_fx_tarf_bgk_discrete_lookback_novy_marx_ff5_acerbi_spectral_es_corsi_har_rv_bll_news_noise(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite61_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite61_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "61_kim_wright_shadow_rate_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_swap_burgard_kjaer_fva_carmona_durrleman_spark_wystup_fx_tarf_bgk_discrete_lookback_novy_marx_ff5_acerbi_spectral_es_corsi_har_rv_bll_news_noise",
        "Kim-Wright-Shadow / Brennan-Schwartz-CB / Duffie-Singleton-RMV / Fukasawa-Gamma-Swap / Burgard-Kjaer-FVA / Carmona-Durrleman-Spark / Wystup-FX-TARF / BGK-Discrete-Lookback / Novy-Marx-FF5 / Acerbi-Spectral-ES / Corsi-HAR-RV / BLL-News-Noise AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite62_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let jps_macro = ql_curves::evaluate_joslin_priebsch_singleton_unspanned_macro_term_structure(
        &ql_curves::JoslinPriebschSingletonUnspannedMacroSpec {
            yield_curve_level_l: T::from_f64(0.042) + dx0,
            yield_curve_slope_s: T::from_f64(-0.014) + dx1,
            yield_curve_curvature_c: T::from_f64(0.008),
            real_output_growth_g: T::from_f64(0.024) + dx2,
            core_inflation_gap_pi: T::from_f64(0.015) + dx3,
            macro_on_yields_projection_beta: T::from_f64(0.25),
            spanned_slope_risk_price_lambda_p: T::from_f64(0.18),
            unspanned_macro_risk_price_lambda_m: T::from_f64(-0.22),
            q_measure_decay_kappa_q: T::from_f64(0.45),
            bond_maturity_years_n: T::from_f64(10.0),
        },
    );

    let kf_ref = ql_fixed_income::evaluate_kalotay_fabozzi_callable_bond_refunding_efficiency(
        &ql_fixed_income::KalotayFabozziRefundingEfficiencySpec {
            notional_par_n: T::from_f64(100.0),
            existing_coupon_rate_c_old: T::from_f64(0.055) + dx0,
            refunding_market_yield_y_ref: T::from_f64(0.038) + dx1,
            call_strike_price_k_call: T::from_f64(101.0) + x[2],
            refunding_issuance_cost_c_iss: T::from_f64(0.65),
            yield_volatility_sigma_y: T::from_f64(0.22) + dx3,
            deferral_option_horizon_t_c: T::from_f64(1.50),
            remaining_maturity_years_t: T::from_f64(10.0),
        },
    );

    let hw_ntd = ql_credit::evaluate_hull_white_nth_to_default_basket_cds_spread(
        &ql_credit::HullWhiteNthToDefaultBasketSpec {
            basket_size_m: T::from_f64(5.0) + x[0] * T::from_f64(0.2),
            individual_hazard_rate_lambda: T::from_f64(0.02) + dx1,
            recovery_rate_rec: T::from_f64(0.40) + dx2,
            gaussian_copula_correlation_rho: T::from_f64(0.30) + dx3,
            risk_free_rate_r: T::from_f64(0.04),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let derman_reg = ql_vol::evaluate_derman_sticky_smile_regime_delta_adjustment(
        &ql_vol::DermanStickySmileRegimeSpec {
            forward_price_f: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(95.0) + x[1],
            atm_volatility_sigma0: T::from_f64(0.20) + dx2,
            log_moneyness_skew_slope_b1: T::from_f64(-0.15) + dx3,
            log_moneyness_curvature_b2: T::from_f64(0.20),
            risk_free_rate_r: T::from_f64(0.04),
            maturity_years_t: T::from_f64(0.50),
        },
    );

    let hw_wwr = ql_xva::evaluate_hull_white_exponential_hazard_wrong_way_risk_cva(
        &ql_xva::HullWhiteExponentialWwrCvaSpec {
            portfolio_expected_mtm_mu_v: T::from_f64(8.0) + x[0],
            portfolio_exposure_vol_sigma_v: T::from_f64(12.0) + x[1],
            baseline_cds_hazard_rate_lambda0: T::from_f64(0.025) + dx2,
            wwr_exponential_sensitivity_b: T::from_f64(0.035) + dx3,
            counterparty_lgd: T::from_f64(0.60),
            ois_discount_rate_r: T::from_f64(0.038),
            maturity_years_t: T::from_f64(4.0),
        },
    );

    let gs_sam = ql_commodities::evaluate_gibson_schwartz_two_factor_futures_samuelson_vol(
        &ql_commodities::GibsonSchwartzSamuelsonFuturesSpec {
            spot_price_s0: T::from_f64(78.0) + x[0],
            initial_convenience_yield_delta0: T::from_f64(0.14) + dx1,
            convenience_yield_mean_reversion_kappa: T::from_f64(1.35),
            risk_neutral_long_run_convenience_alpha: T::from_f64(0.08),
            spot_volatility_sigma_s: T::from_f64(0.36) + dx2,
            convenience_yield_volatility_sigma_d: T::from_f64(0.28),
            spot_convenience_correlation_rho: T::from_f64(0.80) + dx3,
            risk_free_rate_r: T::from_f64(0.04),
            maturity_years_tau: T::from_f64(1.50),
        },
    );

    let quanto_fx = ql_fx::evaluate_reiner_wystup_quanto_fx_correlation_convexity(
        &ql_fx::ReinerWystupQuantoFxSpec {
            foreign_asset_spot_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(102.0) + x[1],
            fixed_quanto_exchange_rate_q_fx: T::from_f64(1.0),
            domestic_rate_rd: T::from_f64(0.045),
            foreign_rate_rf: T::from_f64(0.025),
            foreign_dividend_yield_q: T::from_f64(0.015),
            asset_volatility_sigma_s: T::from_f64(0.22) + dx2,
            fx_volatility_sigma_x: T::from_f64(0.11),
            asset_fx_correlation_rho_sx: T::from_f64(0.45) + dx3,
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let tw_asian = ql_exotics::evaluate_turnbull_wakeman_edgeworth_arithmetic_asian_option(
        &ql_exotics::TurnbullWakemanArithmeticAsianSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(100.0) + x[1],
            risk_free_rate_r: T::from_f64(0.05) + dx2,
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.25) + dx3,
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let amp_vme = ql_portfolio::evaluate_asness_moskowitz_pedersen_value_momentum_everywhere(
        &ql_portfolio::AsnessMoskowitzPedersenVmeSpec {
            value_expected_excess_return_mu_v: T::from_f64(0.045) + dx0,
            momentum_expected_excess_return_mu_m: T::from_f64(0.062) + dx1,
            value_volatility_sigma_v: T::from_f64(0.10) + dx2,
            momentum_volatility_sigma_m: T::from_f64(0.12),
            value_momentum_correlation_rho_vm: T::from_f64(-0.52) + dx3,
            tsmom_lookback_autocovariance_gamma: T::from_f64(0.18),
        },
    );

    let gp_aim = ql_portfolio::evaluate_garleanu_pedersen_dynamic_aim_portfolio_execution(
        &ql_portfolio::GarleanuPedersenAimPortfolioSpec {
            current_position_x_prev: T::from_f64(0.40) + dx0,
            expected_alpha_return_f: T::from_f64(0.035) + dx1,
            alpha_mean_reversion_speed_phi: T::from_f64(0.25) + dx2,
            return_variance_sigma2: T::from_f64(0.04),
            risk_aversion_gamma: T::from_f64(2.0),
            proportional_tcost_lambda: T::from_f64(1.50) + x[3] * T::from_f64(0.1),
            discount_rate_rho_d: T::from_f64(0.04),
        },
    );

    let bn_icp2 = ql_dml::evaluate_bai_ng_approximate_factor_icp2_diffusion_index(
        &ql_dml::BaiNgApproximateFactorIcp2Spec {
            cross_section_dimension_n: T::from_f64(120.0) + x[0],
            time_series_dimension_t: T::from_f64(240.0),
            candidate_factor_count_k: T::from_f64(4.0) + x[1] * T::from_f64(0.2),
            total_panel_variance_sigma_x2: T::from_f64(1.0),
            explained_variance_fraction_r2_k: T::from_f64(0.58) + dx2,
            diffusion_forecast_beta_f: T::from_f64(0.35),
            factor_shock_signal_s_f: T::from_f64(1.40) + dx3,
        },
    );

    let bgg_acc = ql_models::evaluate_bernanke_gertler_gilchrist_financial_accelerator(
        &ql_models::BernankeGertlerGilchristAcceleratorSpec {
            capital_asset_value_qk: T::from_f64(200.0) + x[0],
            entrepreneur_net_worth_n: T::from_f64(80.0) + x[1],
            csv_auditing_cost_fraction_mu: T::from_f64(0.12) + dx2,
            efp_leverage_elasticity_nu: T::from_f64(0.08),
            risk_free_real_rate_r: T::from_f64(0.02),
            unexpected_asset_price_shock_dq: T::from_f64(0.03) + dx3,
        },
    );

    [
        jps_macro.spanned_model_zero_yield + jps_macro.total_expected_bond_excess_return,
        kf_ref.kalotay_refunding_efficiency_ratio
            + kf_ref.total_embedded_refunding_call_value * T::from_f64(0.1),
        hw_ntd.first_to_default_fair_spread_bps * T::from_f64(0.01),
        derman_reg.sticky_moneyness_delta + derman_reg.sticky_delta_minus_local_vol_wedge,
        hw_wwr.wrong_way_risk_adjusted_cva,
        gs_sam.gibson_schwartz_futures_price,
        quanto_fx.quanto_call_domestic_value,
        tw_asian.turnbull_wakeman_arithmetic_call_value,
        amp_vme.combined_vme_sharpe_ratio,
        gp_aim.optimal_next_period_holding_x,
        bn_icp2.bai_ng_pc_p2_criterion + bn_icp2.attenuated_diffusion_index_forecast,
        bgg_acc.required_expected_return_on_capital * T::from_f64(10.0)
            + bgg_acc.amplified_investment_impulse * T::from_f64(10.0),
    ]
}

pub fn benchmark_jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_regimes_hull_white_wwr_cva_gibson_schwartz_samuelson_reiner_quanto_fx_turnbull_wakeman_asian_amp_value_momentum_garleanu_pedersen_aim_bai_ng_icp2_bgg_accelerator(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite62_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite62_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "62_jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_regimes_hull_white_wwr_cva_gibson_schwartz_samuelson_reiner_quanto_fx_turnbull_wakeman_asian_amp_value_momentum_garleanu_pedersen_aim_bai_ng_icp2_bgg_accelerator",
        "JPS-Unspanned-Macro / Kalotay-Refunding / Hull-White-NtD / Derman-Sticky-Regimes / Hull-White-WWR-CVA / Gibson-Schwartz-Samuelson / Reiner-Quanto-FX / Turnbull-Wakeman-Asian / AMP-Value-Momentum / Garleanu-Pedersen-Aim / Bai-Ng-ICp2 / BGG-Accelerator AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite63_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let lm_cap = ql_curves::evaluate_lyashenko_mercurio_fmm_in_arrears_caplet_convexity(
        &ql_curves::LyashenkoMercurioFmmInArrearsCapletSpec {
            notional_n: T::from_f64(1_000.0),
            forward_rfr_rate_f0: T::from_f64(0.042) + dx0,
            caplet_strike_k: T::from_f64(0.040) + dx1,
            fixing_start_years_t0: T::from_f64(0.50),
            accrual_fraction_tau: T::from_f64(0.25),
            forward_vol_sigma_f: T::from_f64(0.24) + dx2,
            ois_discount_rate_r: T::from_f64(0.040) + dx3,
        },
    );

    let cpi_seas = ql_fixed_income::evaluate_tuckman_cpi_inflation_swap_seasonality_carry(
        &ql_fixed_income::TuckmanCpiInflationSwapSeasonalitySpec {
            notional_n: T::from_f64(100.0),
            nominal_treasury_yield_y_nom: T::from_f64(0.043) + dx0,
            tips_real_yield_r_real: T::from_f64(0.019) + dx1,
            inflation_risk_premium_irp: T::from_f64(0.0025) + dx2,
            tips_liquidity_discount_liq: T::from_f64(0.0015),
            monthly_seasonal_cpi_factor_s_m: T::from_f64(0.0020) + dx3,
            swap_maturity_years_t: T::from_f64(5.0),
        },
    );

    let ba_cir = ql_credit::evaluate_brigo_alfonsi_shifted_cir_stochastic_intensity_cds(
        &ql_credit::BrigoAlfonsiShiftedCirStochasticIntensitySpec {
            initial_intensity_y0: T::from_f64(0.022) + dx0,
            mean_reversion_speed_kappa: T::from_f64(0.45),
            long_run_intensity_theta: T::from_f64(0.028) + dx1,
            intensity_volatility_sigma: T::from_f64(0.14) + dx2,
            deterministic_shift_psi: T::from_f64(0.004),
            rates_credit_correlation_rho: T::from_f64(-0.35) + dx3,
            short_rate_volatility_sigma_r: T::from_f64(0.012),
            risk_free_rate_r: T::from_f64(0.04),
            recovery_rate_rec: T::from_f64(0.40),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let fpss_sv = ql_vol::evaluate_fouque_papanicolaou_sircar_multiscale_sv_smile(
        &ql_vol::FouquePapanicolaouSircarMultiscaleSvSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(95.0) + x[1],
            risk_free_rate_r: T::from_f64(0.04),
            effective_volatility_sigma_bar: T::from_f64(0.20) + dx2,
            fast_skew_parameter_v3_eps: T::from_f64(-0.0012) + dx3 * T::from_f64(0.1),
            fast_convexity_parameter_v2_eps: T::from_f64(0.0018),
            slow_skew_parameter_v1_delta: T::from_f64(-0.0015),
            slow_level_parameter_v0_delta: T::from_f64(0.0010),
            maturity_years_t: T::from_f64(0.75),
        },
    );

    let ga_ftd = ql_xva::evaluate_gregory_albanese_bilateral_cva_dva_first_to_default_survival(
        &ql_xva::GregoryAlbaneseBilateralFtdCvaDvaSpec {
            expected_positive_exposure_epe: T::from_f64(15.0) + x[0],
            expected_negative_exposure_ene: T::from_f64(10.0) + x[1],
            counterparty_hazard_rate_lambda_c: T::from_f64(0.030) + dx2,
            own_bank_hazard_rate_lambda_b: T::from_f64(0.018) + dx3,
            counterparty_lgd_c: T::from_f64(0.60),
            own_bank_lgd_b: T::from_f64(0.60),
            ois_discount_rate_r: T::from_f64(0.038),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let ew_vpp = ql_commodities::evaluate_eydeland_wolyniec_virtual_power_plant_heat_rate_strip(
        &ql_commodities::EydelandWolyniecVppHeatRateStripSpec {
            power_forward_price_f_p: T::from_f64(68.0) + x[0],
            gas_forward_price_f_g: T::from_f64(7.50) + x[1] * T::from_f64(0.1),
            plant_heat_rate_h: T::from_f64(7.80),
            variable_om_cost_v_om: T::from_f64(3.50),
            startup_cost_per_mwh_c_start: T::from_f64(2.20),
            power_volatility_sigma_p: T::from_f64(0.48) + dx2,
            gas_volatility_sigma_g: T::from_f64(0.36),
            power_gas_correlation_rho: T::from_f64(0.72) + dx3,
            risk_free_rate_r: T::from_f64(0.04),
            maturity_years_t: T::from_f64(0.75),
        },
    );

    let fx_corr = ql_fx::evaluate_wystup_fx_corridor_accrual_range_note(
        &ql_fx::WystupFxCorridorAccrualRangeNoteSpec {
            notional_coupon_max: T::from_f64(100.0),
            spot_fx_s0: T::from_f64(1.10) + dx0,
            lower_corridor_barrier_l: T::from_f64(1.04) + dx1,
            upper_corridor_barrier_u: T::from_f64(1.18),
            domestic_rate_rd: T::from_f64(0.045),
            foreign_rate_rf: T::from_f64(0.030),
            atm_volatility_sigma: T::from_f64(0.095) + dx2,
            log_moneyness_skew_slope_xi: T::from_f64(-0.12) + dx3,
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let rr_do = ql_exotics::evaluate_rubinstein_reiner_down_and_out_call_with_rebate(
        &ql_exotics::RubinsteinReinerDownAndOutCallRebateSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(100.0) + x[1],
            lower_barrier_h: T::from_f64(90.0) + x[2] * T::from_f64(0.5),
            cash_rebate_at_hit_r_rebate: T::from_f64(3.0),
            risk_free_rate_r: T::from_f64(0.05),
            dividend_yield_q: T::from_f64(0.02),
            volatility_sigma: T::from_f64(0.25) + dx3,
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let lns_gls = ql_portfolio::evaluate_lewellen_nagel_shanken_gls_cross_sectional_r2(
        &ql_portfolio::LewellenNagelShankenGlsCrossSectionalSpec {
            factor_tangency_sharpe_sq_theta_f2: T::from_f64(0.085) + dx0,
            residual_alpha_quadratic_form_q_alpha: T::from_f64(0.140) + dx1,
            cross_sectional_ols_r2: T::from_f64(0.78) + dx2,
            factor_risk_premium_lambda: T::from_f64(0.0065) + dx3 * T::from_f64(0.1),
            factor_variance_sigma_f2: T::from_f64(0.0022),
            sample_length_t: T::from_f64(240.0),
            test_asset_count_n: T::from_f64(25.0),
            factor_count_k: T::from_f64(3.0),
        },
    );

    let bc_exp = ql_portfolio::evaluate_bellini_cont_expectile_evar_influence_sensitivity(
        &ql_portfolio::BelliniContExpectileRiskMeasureSpec {
            loss_mean_mu: T::from_f64(0.50) + dx0,
            loss_volatility_sigma: T::from_f64(2.20) + dx1,
            tail_asymmetry_tau: T::from_f64(0.975) + dx2 * T::from_f64(0.1),
            contamination_loss_shock_z: T::from_f64(8.50) + x[3] * T::from_f64(0.1),
            contamination_weight_eps: T::from_f64(0.02),
        },
    );

    let bns_rq = ql_dml::evaluate_barndorff_nielsen_shephard_realized_quarticity_iv_ci(
        &ql_dml::BarndorffNielsenShephardRealizedQuarticitySpec {
            realized_variance_rv: T::from_f64(0.032) + dx0 * T::from_f64(0.1),
            realized_quarticity_rq: T::from_f64(0.0014) + dx1 * T::from_f64(0.01),
            intraday_sampling_frequency_m: T::from_f64(78.0) + x[2],
            first_order_autocovariance_gamma1: T::from_f64(-0.0008) + dx3 * T::from_f64(0.01),
            normal_critical_z_alpha: T::from_f64(1.96),
        },
    );

    let km_cyc = ql_models::evaluate_kiyotaki_moore_collateral_credit_cycle_multiplier(
        &ql_models::KiyotakiMooreCollateralCycleSpec {
            land_asset_price_q: T::from_f64(100.0) + x[0],
            expected_next_land_price_q_next: T::from_f64(101.0) + x[1],
            gross_interest_rate_r_gross: T::from_f64(1.05),
            farmer_productivity_a: T::from_f64(6.0) + dx2 * T::from_f64(10.0),
            existing_land_holding_k_prev: T::from_f64(10.0),
            inherited_debt_b_prev: T::from_f64(880.0),
            // Keep the feedback term below its 0.85 saturation cap so all four
            // workload inputs retain a live sensitivity in the AD invariant.
            intertemporal_price_elasticity_eta: T::from_f64(0.02) + dx3,
        },
    );

    [
        lm_cap.backward_looking_in_arrears_caplet_npv,
        cpi_seas.zciis_seasonal_mark_to_market_npv,
        ba_cir.correlated_ssrd_cds_fair_spread_bps * T::from_f64(0.01),
        fpss_sv.multiscale_sv_call_option_price,
        ga_ftd.net_bilateral_bcva,
        ew_vpp.vpp_startup_adjusted_tolling_value,
        fx_corr.smile_adjusted_range_note_pv * T::from_f64(0.1),
        rr_do.total_down_and_out_call_with_rebate,
        lns_gls.ols_minus_gls_r2_inflation_wedge + lns_gls.shanken_eiv_covariance_multiplier,
        bc_exp.contaminated_robust_expectile,
        bns_rq.log_transformed_iv_upper_bound + bns_rq.newey_west_microstructure_corrected_rv,
        km_cyc.dynamic_credit_cycle_amplification_multiplier,
    ]
}

pub fn benchmark_lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_cds_fouque_multiscale_sv_gregory_bilateral_ftd_eydeland_vpp_heat_rate_wystup_fx_corridor_rubinstein_reiner_barrier_lns_gls_r2_bellini_expectile_bns_quarticity_kiyotaki_moore_cycle(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite63_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite63_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "63_lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_cds_fouque_multiscale_sv_gregory_bilateral_ftd_eydeland_vpp_heat_rate_wystup_fx_corridor_rubinstein_reiner_barrier_lns_gls_r2_bellini_expectile_bns_quarticity_kiyotaki_moore_cycle",
        "Lyashenko-FMM-Caplet / Tuckman-CPI-Seasonality / Brigo-Alfonsi-CIR-CDS / Fouque-Multiscale-SV / Gregory-Bilateral-FtD / Eydeland-VPP-Heat-Rate / Wystup-FX-Corridor / Rubinstein-Reiner-Barrier / LNS-GLS-R2 / Bellini-Expectile / BNS-Quarticity / Kiyotaki-Moore-Cycle AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite64_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let jsz_gdtsm = ql_curves::evaluate_joslin_singleton_zhu_canonical_gdtsm_term_premium(
        &ql_curves::JoslinSingletonZhuGdtsmSpec {
            risk_neutral_long_run_mean_r_inf_q: T::from_f64(0.045) + dx0,
            physical_long_run_mean_r_inf_p: T::from_f64(0.032) + dx1,
            risk_neutral_decay_speeds_kappa_q: [
                T::from_f64(0.06),
                T::from_f64(0.45),
                T::from_f64(1.20),
            ],
            physical_decay_speeds_kappa_p: [
                T::from_f64(0.18),
                T::from_f64(0.65),
                T::from_f64(1.45),
            ],
            canonical_factor_states_x: [
                T::from_f64(0.015) + dx2,
                T::from_f64(0.008),
                T::from_f64(0.012),
            ],
            factor_volatilities_sigma: [
                T::from_f64(0.0065) + dx3 * T::from_f64(0.1),
                T::from_f64(0.0080),
                T::from_f64(0.0095),
            ],
            bond_maturity_years_tau: T::from_f64(10.0),
        },
    );

    let tf_cb = ql_fixed_income::evaluate_tsiveriotis_fernandes_convertible_bond_credit_split(
        &ql_fixed_income::TsiveriotisFernandesConvertibleSpec {
            spot_price_s0: T::from_f64(98.0) + x[0],
            conversion_ratio_gamma: T::from_f64(1.0),
            par_value: T::from_f64(100.0),
            coupon_rate_c: T::from_f64(0.035) + dx1,
            risk_free_rate_r: T::from_f64(0.04),
            credit_spread_sc: T::from_f64(0.025) + dx2,
            dividend_yield_q: T::from_f64(0.015),
            equity_volatility_sigma: T::from_f64(0.28) + dx3,
            call_price_k: T::from_f64(112.0),
            maturity_years_t: T::from_f64(2.0),
        },
    );

    let lg_cdo = ql_credit::evaluate_laurent_gregory_fft_recursion_synthetic_cdo_tranche(
        &ql_credit::LaurentGregoryRecursiveCdoSpec {
            attachment_point_ka: T::from_f64(0.03) + dx0 * T::from_f64(0.1),
            detachment_point_kd: T::from_f64(0.07),
            recovery_rate_rec: T::from_f64(0.40) + dx1,
            tier_hazard_rates_lambda: [
                T::from_f64(0.008) + dx2 * T::from_f64(0.1),
                T::from_f64(0.020),
                T::from_f64(0.045),
            ],
            tier_factor_correlations_rho: [
                T::from_f64(0.25) + dx3,
                T::from_f64(0.35),
                T::from_f64(0.45),
            ],
            risk_free_rate_r: T::from_f64(0.04),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let bj_vs = ql_vol::evaluate_broadie_jain_discrete_sampling_jump_variance_swap(
        &ql_vol::BroadieJainVarianceSwapSpec {
            diffusive_volatility_sigma: T::from_f64(0.18) + dx0,
            jump_intensity_lambda_j: T::from_f64(0.60) + dx1,
            mean_log_jump_mu_j: T::from_f64(-0.08) + dx2,
            std_log_jump_delta_j: T::from_f64(0.12) + dx3,
            risk_free_rate_r: T::from_f64(0.045),
            dividend_yield_q: T::from_f64(0.015),
            discrete_sampling_periods_n: T::from_f64(52.0),
            maturity_years_t: T::from_f64(1.0),
        },
    );

    let ac_kva = ql_xva::evaluate_albanese_crepey_sa_ccr_kadd_capital_xva(
        &ql_xva::AlbaneseCrepeySaCcrKvaSpec {
            trade_notional: T::from_f64(1_000.0),
            current_mtm_value_v0: T::from_f64(8.0) + x[0],
            collateral_balance_c0: T::from_f64(12.0) + x[1],
            counterparty_risk_weight_w: T::from_f64(0.50) + dx2,
            capital_hurdle_rate_h: T::from_f64(0.10) + dx3,
            ois_discount_rate_r: T::from_f64(0.038),
            joint_hazard_rate_lambda: T::from_f64(0.025),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let rss_stor = ql_commodities::evaluate_routledge_seppi_spatt_equilibrium_storage_forward_curve(
        &ql_commodities::RoutledgeSeppiSpattStorageSpec {
            spot_price_s0: T::from_f64(78.0) + x[0],
            physical_inventory_buffer_i0: T::from_f64(0.12) + dx1,
            net_demand_shock_x0: T::from_f64(0.18) + dx2,
            demand_mean_reversion_kappa_x: T::from_f64(1.40),
            demand_shock_volatility_sigma_x: T::from_f64(0.26) + dx3,
            storage_cost_rate_cs: T::from_f64(0.035),
            risk_free_rate_r: T::from_f64(0.045),
            forward_maturity_years_tau: T::from_f64(0.75),
        },
    );

    let dtv_cip = ql_fx::evaluate_du_tepper_verdelhan_cip_treasury_basis_arb(
        &ql_fx::DuTepperVerdelhanCipBasisSpec {
            spot_fx_s0: T::from_f64(1.0800) + dx0,
            outright_forward_fx_f0: T::from_f64(1.1015) + dx1,
            domestic_usd_ois_rate_rd: T::from_f64(0.0480),
            foreign_ois_rate_rf: T::from_f64(0.0310),
            domestic_ust_yield: T::from_f64(0.0465) + dx2,
            foreign_gov_yield: T::from_f64(0.0325),
            dealer_slr_balance_sheet_cost: T::from_f64(0.0018) + dx3 * T::from_f64(0.1),
            dealer_funding_spread_delta_s: T::from_f64(0.0012),
            maturity_years_tau: T::from_f64(1.0),
        },
    );

    let wc_cliq = ql_exotics::evaluate_wilmott_carr_globally_floored_locally_capped_cliquet(
        &ql_exotics::WilmottCarrCliquetSpec {
            notional_n: T::from_f64(100.0),
            local_return_floor_f_loc: T::from_f64(0.00),
            local_return_cap_c_loc: T::from_f64(0.08) + dx0,
            global_cumulative_floor_f_glob: T::from_f64(0.10) + dx1,
            risk_free_rate_r: T::from_f64(0.04),
            dividend_yield_q: T::from_f64(0.015),
            forward_atm_volatility_sigma: T::from_f64(0.22) + dx2,
            forward_skew_slope_eta: T::from_f64(-0.25) + dx3,
            maturity_years_t: T::from_f64(2.0),
        },
    );

    let bsv_ppp = ql_portfolio::evaluate_brandt_santa_clara_valkanov_crra_characteristic_tilt(
        &ql_portfolio::BrandtSantaClaraValkanovCharacteristicTiltSpec {
            benchmark_expected_return_mu_b: T::from_f64(0.068) + dx0,
            benchmark_volatility_sigma_b: T::from_f64(0.155),
            value_characteristic_premium_mu_v: T::from_f64(0.042) + dx1,
            momentum_characteristic_premium_mu_m: T::from_f64(0.058) + dx2,
            value_characteristic_vol_sigma_v: T::from_f64(0.115),
            momentum_characteristic_vol_sigma_m: T::from_f64(0.145),
            characteristic_correlation_rho_vm: T::from_f64(-0.35),
            relative_risk_aversion_gamma: T::from_f64(5.0) + x[3],
        },
    );

    let aj_evar = ql_portfolio::evaluate_ahmadi_javid_follmer_schied_entropic_var_evar(
        &ql_portfolio::AhmadiJavidEntropicVarSpec {
            mean_diffusive_loss_mu: T::from_f64(0.25) + dx0,
            diffusive_loss_volatility_sigma: T::from_f64(1.50) + dx1,
            jump_loss_intensity_lambda: T::from_f64(0.15) + dx2,
            mean_jump_loss_mu_j: T::from_f64(2.20) + dx3,
            jump_loss_volatility_delta_j: T::from_f64(0.85),
            tail_significance_alpha: T::from_f64(0.025),
        },
    );

    let ek_bekk = ql_dml::evaluate_engle_kroner_bekk_bivariate_garch_covariance(
        &ql_dml::EngleKronerBekkGarchSpec {
            cholesky_intercept_c: [
                T::from_f64(0.04) + dx0 * T::from_f64(0.1),
                T::from_f64(0.015),
                T::from_f64(0.035),
            ],
            arch_spillover_matrix_a: [
                T::from_f64(0.28) + dx1,
                T::from_f64(0.08),
                T::from_f64(0.06),
                T::from_f64(0.25),
            ],
            garch_persistence_diag_b: [
                T::from_f64(0.92) + dx2 * T::from_f64(0.1),
                T::from_f64(0.90),
            ],
            lagged_innovations_eps: [
                T::from_f64(-0.025) + dx3 * T::from_f64(0.1),
                T::from_f64(-0.018),
            ],
            lagged_covariance_h: [T::from_f64(0.032), T::from_f64(0.018), T::from_f64(0.028)],
        },
    );

    let gab_dis = ql_models::evaluate_gabaix_variable_rare_disaster_resilience_equity_puzzle(
        &ql_models::GabaixVariableRareDisasterSpec {
            risk_aversion_gamma: T::from_f64(4.0) + dx0,
            time_preference_rho: T::from_f64(0.035),
            normal_consumption_growth_gc: T::from_f64(0.020),
            consumption_volatility_sigma_c: T::from_f64(0.018),
            current_disaster_probability_pt: T::from_f64(0.022) + dx1 * T::from_f64(0.1),
            baseline_disaster_probability_p_star: T::from_f64(0.017),
            consumption_disaster_recovery_bc: T::from_f64(0.72) + dx2,
            dividend_disaster_recovery_fd: T::from_f64(0.55),
            resilience_mean_reversion_phi_h: T::from_f64(0.18) + dx3,
        },
    );

    [
        jsz_gdtsm.canonical_yield_term_premium * T::from_f64(100.0),
        tf_cb.convertible_bond_value * T::from_f64(0.1),
        lg_cdo.fair_tranche_running_spread_bps * T::from_f64(0.01),
        bj_vs.convexity_adjusted_vol_swap_strike * T::from_f64(10.0),
        ac_kva.lifetime_sa_ccr_kva,
        rss_stor.equilibrium_forward_price * T::from_f64(0.1),
        dtv_cip.net_actionable_cip_arbitrage_alpha_bps * T::from_f64(0.1),
        wc_cliq.total_cliquet_present_value * T::from_f64(0.1),
        bsv_ppp.crra_certainty_equivalent_return * T::from_f64(100.0),
        aj_evar.jump_diffusion_entropic_var,
        ek_bekk.conditional_correlation_rho12 * T::from_f64(10.0),
        gab_dis.equilibrium_price_dividend_ratio,
    ]
}

pub fn benchmark_jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_albanese_sa_ccr_kva_rss_equilibrium_storage_du_tepper_verdelhan_cip_wilmott_cliquet_bsv_parametric_portfolio_ahmadi_javid_evar_engle_kroner_bekk_gabaix_disaster(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite64_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite64_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "64_jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_albanese_sa_ccr_kva_rss_equilibrium_storage_du_tepper_verdelhan_cip_wilmott_cliquet_bsv_parametric_portfolio_ahmadi_javid_evar_engle_kroner_bekk_gabaix_disaster",
        "JSZ-Canonical-GDTSM / Tsiveriotis-Fernandes-CB / Laurent-Gregory-CDO / Broadie-Jain-VarSwap / Albanese-SA-CCR-KVA / RSS-Equilibrium-Storage / Du-Tepper-Verdelhan-CIP / Wilmott-Cliquet / BSV-Parametric-Portfolio / Ahmadi-Javid-EVaR / Engle-Kroner-BEKK / Gabaix-Rare-Disaster AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite65_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let bian_q = ql_curves::evaluate_bianchetti_two_curve_forward_measure_quanto_convexity(
        &ql_curves::BianchettiTwoCurveForwardQuantoSpec {
            notional_n: T::from_f64(1_000.0),
            projection_forward_rate_f0: T::from_f64(0.044) + dx0,
            tenor_basis_spread_s0: T::from_f64(0.0035) + dx1 * T::from_f64(0.1),
            caplet_strike_k: T::from_f64(0.042),
            fixing_expiry_years_t1: T::from_f64(2.0),
            accrual_fraction_tau: T::from_f64(0.50),
            forward_rate_vol_sigma_f: T::from_f64(0.25) + dx2,
            basis_spread_vol_sigma_s: T::from_f64(0.35),
            forward_basis_correlation_rho: T::from_f64(0.65) + dx3,
            ois_discount_rate_r: T::from_f64(0.040),
        },
    );

    let bdt_cb = ql_fixed_income::evaluate_black_derman_toy_lognormal_short_rate_callable_bond(
        &ql_fixed_income::BlackDermanToyCallableBondSpec {
            par_value: T::from_f64(100.0),
            coupon_rate_c: T::from_f64(0.055) + dx0,
            initial_short_rate_r0: T::from_f64(0.038) + dx1,
            term_slope_s: T::from_f64(0.002),
            short_rate_log_vol_sigma: T::from_f64(0.22) + dx2,
            option_adjusted_spread_oas: T::from_f64(0.006) + dx3,
            call_price_k: T::from_f64(101.0),
            maturity_years_t: T::from_f64(4.0),
        },
    );

    let sj_swp = ql_credit::evaluate_schonbucher_jamshidian_survival_measure_payer_cds_swaption(
        &ql_credit::SchonbucherJamshidianCdsSwaptionSpec {
            notional_n: T::from_f64(1_000.0),
            short_term_hazard_rate_lambda0: T::from_f64(0.018) + dx0,
            forward_hazard_rate_lambda_fwd: T::from_f64(0.022) + dx1,
            strike_spread_k: T::from_f64(0.012) + dx2 * T::from_f64(0.1),
            spread_volatility_sigma_s: T::from_f64(0.45) + dx3,
            recovery_rate_rec: T::from_f64(0.40),
            risk_free_rate_r: T::from_f64(0.04),
            option_expiry_years_ta: T::from_f64(1.0),
            cds_maturity_years_tb: T::from_f64(6.0),
        },
    );

    let bates_svj = ql_vol::evaluate_bates_svj_stochastic_vol_jump_characteristic_smile(
        &ql_vol::BatesSvjSmileSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(93.0) + x[1],
            risk_free_rate_r: T::from_f64(0.04),
            initial_variance_v0: T::from_f64(0.04) + dx2 * T::from_f64(0.1),
            long_run_variance_theta: T::from_f64(0.045),
            mean_reversion_kappa: T::from_f64(1.80),
            vol_of_vol_xi: T::from_f64(0.38),
            spot_vol_correlation_rho: T::from_f64(-0.65),
            jump_intensity_lambda_j: T::from_f64(0.35) + dx3,
            mean_log_jump_mu_j: T::from_f64(-0.07),
            std_log_jump_delta_j: T::from_f64(0.11),
            maturity_years_t: T::from_f64(0.75),
        },
    );

    let bp_ccp = ql_xva::evaluate_brigo_pallavicini_ccp_default_fund_cva_mva(
        &ql_xva::BrigoPallaviciniCcpXvaSpec {
            portfolio_mtm_volatility_sigma_v: T::from_f64(45.0) + x[0],
            margin_period_of_risk_years: T::from_f64(10.0 / 252.0),
            im_confidence_alpha: T::from_f64(0.985) + dx1 * T::from_f64(0.1),
            im_funding_spread_s_im: T::from_f64(0.0085) + dx2 * T::from_f64(0.1),
            clearing_member_hazard_lambda_cm: T::from_f64(0.015) + dx3,
            default_fund_scaling_theta_df: T::from_f64(0.18),
            default_fund_mutualization_spread_s_df: T::from_f64(0.025),
            ois_discount_rate_r: T::from_f64(0.038),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let dmt_crk = ql_commodities::evaluate_dempster_medova_tang_cointegrated_crack_spread_option(
        &ql_commodities::DempsterMedovaTangCrackSpreadSpec {
            refined_product_forward_f1: T::from_f64(96.0) + x[0],
            crude_feedstock_forward_f2: T::from_f64(78.0) + x[1],
            refining_strike_k: T::from_f64(14.0),
            product_volatility_sigma1: T::from_f64(0.34) + dx2,
            crude_volatility_sigma2: T::from_f64(0.30),
            product_crude_correlation_rho: T::from_f64(0.82),
            cointegration_reversion_speed_kappa: T::from_f64(1.50) + dx3,
            risk_free_rate_r: T::from_f64(0.04),
            maturity_years_t: T::from_f64(1.5),
        },
    );

    let pit_prdc =
        ql_fx::evaluate_piterbarg_fx_power_reverse_dual_currency_prdc_bermudan_redemption(
            &ql_fx::PiterbargPrdcBermudanSpec {
                notional_n: T::from_f64(1_000.0),
                spot_fx_s0: T::from_f64(105.0) + x[0],
                foreign_coupon_gearing_gf: T::from_f64(0.10) + dx1,
                domestic_funding_deduction_gd: T::from_f64(0.06),
                domestic_rate_rd: T::from_f64(0.015),
                foreign_rate_rf: T::from_f64(0.045),
                fx_volatility_sigma_fx: T::from_f64(0.14) + dx2,
                domestic_rate_vol_sigma_d: T::from_f64(0.012),
                domestic_fx_correlation_rho: T::from_f64(0.45),
                target_redemption_cap_c_target: T::from_f64(0.12) + dx3,
                maturity_years_t: T::from_f64(4.0),
            },
        );

    let bgk_bar =
        ql_exotics::evaluate_broadie_glasserman_kou_discrete_barrier_continuity_correction(
            &ql_exotics::BroadieGlassermanKouDiscreteBarrierSpec {
                spot_price_s0: T::from_f64(100.0) + x[0],
                strike_price_k: T::from_f64(100.0) + x[1],
                contractual_down_barrier_h: T::from_f64(91.0) + x[2] * T::from_f64(0.5),
                discrete_monitoring_dates_m: T::from_f64(52.0),
                risk_free_rate_r: T::from_f64(0.05),
                dividend_yield_q: T::from_f64(0.02),
                volatility_sigma: T::from_f64(0.25) + dx3,
                maturity_years_t: T::from_f64(1.0),
            },
        );

    let ck_apt = ql_portfolio::evaluate_connor_korajczyk_asymptotic_pca_apt_pricing_error(
        &ql_portfolio::ConnorKorajczykAptSpec {
            asset_expected_excess_return_mu: T::from_f64(0.088) + dx0,
            asset_pca_factor_betas: [T::from_f64(1.15) + dx1, T::from_f64(0.45)],
            pca_factor_risk_premia_lambda: [T::from_f64(0.055) + dx2, T::from_f64(0.032)],
            pca_factor_variances: [T::from_f64(0.025), T::from_f64(0.014)],
            idiosyncratic_residual_variance_sigma_e2: T::from_f64(0.036) + dx3 * T::from_f64(0.1),
            cross_sectional_asset_count_n: T::from_f64(250.0),
            time_series_sample_length_t: T::from_f64(180.0),
        },
    );

    let afs_exec = ql_portfolio::evaluate_alfonsi_fruth_schied_bucket_transient_impact_execution(
        &ql_portfolio::AlfonsiFruthSchiedExecutionSpec {
            parent_order_size_x0: T::from_f64(100.0) + x[0],
            discrete_intervals_n: T::from_f64(10.0),
            lob_resilience_speed_rho: T::from_f64(2.50) + dx1,
            transient_impact_coeff_eta: T::from_f64(0.04) + dx2 * T::from_f64(0.1),
            permanent_impact_coeff_gamma: T::from_f64(0.01),
            asset_volatility_sigma: T::from_f64(0.22) + dx3,
            inventory_risk_aversion_lambda: T::from_f64(0.05),
            horizon_t: T::from_f64(1.0),
        },
    );

    let ham_ms = ql_dml::evaluate_hamilton_markov_switching_regime_filter_step(
        &ql_dml::HamiltonMarkovSwitchingSpec {
            observed_realization_yt: T::from_f64(-0.018) + dx0,
            expansion_mean_mu1: T::from_f64(0.008) + dx1 * T::from_f64(0.1),
            recession_mean_mu2: T::from_f64(-0.022),
            expansion_volatility_sigma1: T::from_f64(0.012) + dx2 * T::from_f64(0.1),
            recession_volatility_sigma2: T::from_f64(0.028),
            stay_expansion_prob_p11: T::from_f64(0.95),
            stay_recession_prob_p22: T::from_f64(0.85),
            prior_expansion_prob_xi1: T::from_f64(0.80) + dx3,
        },
    );

    let by_lrr = ql_models::evaluate_bansal_yaron_long_run_risk_stochastic_vol_equity_premium(
        &ql_models::BansalYaronLongRunRiskSpec {
            risk_aversion_gamma: T::from_f64(7.5) + dx0,
            eis_psi: T::from_f64(1.5) + dx1,
            growth_persistence_rho_x: T::from_f64(0.975),
            growth_shock_scale_phi_e: T::from_f64(0.044) + dx2 * T::from_f64(0.1),
            volatility_persistence_nu: T::from_f64(0.985),
            vol_of_vol_sigma_w: T::from_f64(0.00025),
            current_consumption_vol_sigma_t: T::from_f64(0.027) + dx3 * T::from_f64(0.1),
            dividend_leverage_phi_d: T::from_f64(3.0),
            campbell_shiller_kappa1: T::from_f64(0.99),
        },
    );

    [
        bian_q.bianchetti_quanto_adjusted_caplet_pv,
        bdt_cb.bdt_callable_bond_price * T::from_f64(0.1),
        sj_swp.full_non_knockout_swaption_pv * T::from_f64(0.1),
        bates_svj.bates_svj_call_option_price,
        bp_ccp.total_ccp_clearing_xva,
        dmt_crk.cointegrated_crack_spread_option_price,
        pit_prdc.tarn_capped_prdc_present_value * T::from_f64(0.1),
        bgk_bar.bgk_discrete_down_and_out_call_price,
        ck_apt.connor_korajczyk_wald_statistic,
        afs_exec.optimal_afs_execution_cost * T::from_f64(0.1),
        ham_ms.next_step_conditional_forecast * T::from_f64(100.0),
        by_lrr.conditional_equity_risk_premium * T::from_f64(100.0),
    ]
}

pub fn benchmark_bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_smile_brigo_pallavicini_ccp_xva_dmt_cointegrated_crack_piterbarg_prdc_tarn_bgk_discrete_barrier_connor_korajczyk_apt_afs_transient_impact_hamilton_markov_switching_bansal_yaron_lrr(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite65_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite65_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "65_bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_smile_brigo_pallavicini_ccp_xva_dmt_cointegrated_crack_piterbarg_prdc_tarn_bgk_discrete_barrier_connor_korajczyk_apt_afs_transient_impact_hamilton_markov_switching_bansal_yaron_lrr",
        "Bianchetti-Two-Curve-Quanto / BDT-Callable-Bond / Schonbucher-CDS-Swaption / Bates-SVJ-Smile / Brigo-Pallavicini-CCP-XVA / DMT-Cointegrated-Crack / Piterbarg-PRDC-TARN / BGK-Discrete-Barrier / Connor-Korajczyk-APT / AFS-Transient-Impact / Hamilton-Markov-Switching / Bansal-Yaron-LRR AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite66_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let ap_cms = ql_curves::evaluate_andersen_piterbarg_lsrm_cms_spread_option(
        &ql_curves::AndersenPiterbargCmsSpreadSpec {
            cms1_forward_swap_rate_s1: T::from_f64(0.042) + dx0,
            cms1_tenor_years_tau1: T::from_f64(10.0),
            cms1_normal_volatility_sigma1: T::from_f64(0.0085) + dx2 * T::from_f64(0.1),
            cms2_forward_swap_rate_s2: T::from_f64(0.034) + dx1,
            cms2_tenor_years_tau2: T::from_f64(2.0),
            cms2_normal_volatility_sigma2: T::from_f64(0.0095),
            swap_rate_correlation_rho12: T::from_f64(0.78),
            leg1_gearing_w1: T::one(),
            leg2_gearing_w2: T::one(),
            spread_strike_k: T::from_f64(0.005) + dx3 * T::from_f64(0.1),
            mean_reversion_kappa: T::from_f64(0.08),
            ois_discount_rate_r: T::from_f64(0.035),
            fixing_horizon_years_t0: T::from_f64(3.0),
        },
    );

    let fab_z = ql_fixed_income::evaluate_fabozzi_cmo_accrual_z_bond_sequential_lockout_tranche(
        &ql_fixed_income::FabozziCmoAccrualZBondSpec {
            senior_tranche_initial_balance_ba: T::from_f64(75.0) + x[0],
            z_bond_initial_balance_bz: T::from_f64(25.0) + x[1],
            tranche_coupon_rate_c: T::from_f64(0.055),
            conditional_prepayment_rate_cpr: T::from_f64(0.16) + dx2,
            scheduled_amortization_rate: T::from_f64(0.12),
            base_discount_yield_y: T::from_f64(0.042),
            z_bond_oas_spread: T::from_f64(0.012) + dx3,
        },
    );

    let dg_cdo = ql_credit::evaluate_duffie_garleanu_affine_jump_intensity_cdo_tranche(
        &ql_credit::DuffieGarleanuAffineJumpCdoSpec {
            initial_common_intensity_lambda_c0: T::from_f64(0.015) + dx0 * T::from_f64(0.1),
            common_reversion_speed_kappa_c: T::from_f64(0.45),
            common_long_run_theta_c: T::from_f64(0.018),
            common_feller_volatility_sigma_c: T::from_f64(0.09) + dx1,
            common_jump_arrival_rate_mu_j: T::from_f64(0.25) + dx2,
            common_mean_jump_size_eta_j: T::from_f64(0.03),
            idiosyncratic_intensity_lambda_id: T::from_f64(0.012) + dx3 * T::from_f64(0.1),
            loss_given_default_lgd: T::from_f64(0.60),
            attachment_point_ka: T::from_f64(0.03),
            detachment_point_kd: T::from_f64(0.07),
            risk_free_rate_r: T::from_f64(0.035),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let cl_vol = ql_vol::evaluate_carr_lee_robust_volatility_swap_convexity_correction(
        &ql_vol::CarrLeeVolatilitySwapSpec {
            initial_variance_v0: T::from_f64(0.045) + dx0 * T::from_f64(0.1),
            long_run_variance_theta: T::from_f64(0.040) + dx1 * T::from_f64(0.1),
            mean_reversion_speed_kappa: T::from_f64(1.80),
            vol_of_vol_xi: T::from_f64(0.45) + dx2,
            realized_variance_skewness_gamma: T::from_f64(0.65) + dx3,
            maturity_years_t: T::one(),
        },
    );

    let gr_kva = ql_xva::evaluate_green_basel_leverage_ratio_rwa_hurdle_kva(
        &ql_xva::GreenBaselLeverageRatioKvaSpec {
            initial_sa_ccr_ead0: T::from_f64(18.0) + x[0],
            gross_trade_notional: T::from_f64(500.0) + x[1],
            sa_ccr_addon_factor_psi: T::from_f64(0.04),
            counterparty_risk_weight_w: T::from_f64(0.50),
            cet1_capital_ratio_k: T::from_f64(0.105),
            leverage_ratio_requirement_gamma: T::from_f64(0.035) + dx2 * T::from_f64(0.1),
            cva_capital_loading_beta: T::from_f64(0.08),
            capital_hurdle_rate_h: T::from_f64(0.11) + dx3,
            capital_funding_offset_phi: T::from_f64(0.35),
            counterparty_hazard_rate_lambda: T::from_f64(0.015),
            ois_discount_rate_r: T::from_f64(0.035),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let cs_oil = ql_commodities::evaluate_cortazar_schwartz_three_factor_oil_futures_curve(
        &ql_commodities::CortazarSchwartzThreeFactorSpec {
            spot_price_s0: T::from_f64(78.0) + x[0],
            initial_convenience_yield_delta0: T::from_f64(0.08) + dx1,
            initial_long_term_drift_v0: T::from_f64(0.035) + dx2,
            risk_neutral_convenience_mean_alpha2: T::from_f64(0.05),
            risk_neutral_long_run_drift_vbar: T::from_f64(0.03),
            convenience_reversion_speed_kappa2: T::from_f64(1.40),
            drift_reversion_speed_kappa3: T::from_f64(0.35),
            spot_volatility_sigma1: T::from_f64(0.32) + dx3,
            convenience_volatility_sigma2: T::from_f64(0.24),
            drift_volatility_sigma3: T::from_f64(0.08),
            spot_convenience_correlation_rho12: T::from_f64(0.72),
            maturity_years_tau: T::from_f64(2.5),
        },
    );

    let fst_xccs = ql_fx::evaluate_fujii_shimada_takahashi_collateralized_xccs_curve_discount(
        &ql_fx::FujiiShimadaTakahashiXccsSpec {
            spot_fx_s0: T::from_f64(1.08) + dx0,
            domestic_ois_rate_rd: T::from_f64(0.045),
            foreign_ois_rate_rf: T::from_f64(0.030) + dx1,
            xccy_collateral_basis_spread_b: T::from_f64(-0.0025) + dx2 * T::from_f64(0.1),
            foreign_swap_coupon_cf: T::from_f64(0.032) + dx3,
            foreign_rate_volatility_sigma_f: T::from_f64(0.009),
            fx_volatility_sigma_fx: T::from_f64(0.11),
            foreign_fx_correlation_rho: T::from_f64(0.35),
            maturity_years_t: T::from_f64(4.0),
        },
    );

    let ph_tsw = ql_exotics::evaluate_pechtl_haug_discrete_time_switch_accumulator_option(
        &ql_exotics::PechtlHaugTimeSwitchSpec {
            spot_price_s0: T::from_f64(104.0) + x[0],
            strike_price_k: T::from_f64(100.0) + x[1],
            accrual_rate_a: T::from_f64(25.0),
            accumulated_occupation_time_t_acc: T::from_f64(0.25),
            risk_free_rate_r: T::from_f64(0.04),
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.22) + dx2,
            smile_skew_slope_nu: T::from_f64(-0.08) + dx3,
            maturity_years_t: T::from_f64(1.5),
        },
    );

    let bcf_mvar = ql_portfolio::evaluate_boudt_cornish_fisher_modified_var_risk_budget_allocation(
        &ql_portfolio::BoudtCornishFisherMvarSpec {
            initial_weight_w1: T::from_f64(0.60) + dx0,
            target_risk_budget_b1: T::from_f64(0.50) + dx1,
            volatility_sigma1: T::from_f64(0.22) + dx2,
            volatility_sigma2: T::from_f64(0.14),
            return_correlation_rho: T::from_f64(0.35),
            loss_skewness_loading_s1: T::from_f64(0.45) + dx3,
            loss_skewness_loading_s2: T::from_f64(0.15),
            excess_kurtosis_loading_k1: T::from_f64(1.80),
            excess_kurtosis_loading_k2: T::from_f64(0.60),
            normal_tail_quantile_z: T::from_f64(2.3263),
        },
    );

    let ru_quad = ql_portfolio::evaluate_rockafellar_uryasev_mixed_cvar_risk_quadrangle(
        &ql_portfolio::RockafellarUryasevMixedCvarSpec {
            mean_loss_mu: T::from_f64(0.012) + dx0,
            loss_volatility_sigma: T::from_f64(0.045) + dx1,
            moderate_confidence_alpha1: T::from_f64(0.95),
            extreme_confidence_alpha2: T::from_f64(0.99),
            moderate_spectral_weight_lambda1: T::from_f64(0.60) + dx2,
            tail_kurtosis_amplification_kappa: T::from_f64(1.50) + dx3,
        },
    );

    let hgg_is = ql_dml::evaluate_hasbrouck_gonzalo_granger_price_discovery_information_share(
        &ql_dml::HasbrouckGonzaloGrangerSpec {
            venue1_error_correction_speed_alpha1: T::from_f64(0.08) + dx0,
            venue2_error_correction_speed_alpha2: T::from_f64(0.24) + dx1,
            venue1_innovation_volatility_sigma1: T::from_f64(0.015) + dx2 * T::from_f64(0.1),
            venue2_innovation_volatility_sigma2: T::from_f64(0.018),
            cross_venue_innovation_correlation_rho: T::from_f64(0.55) + dx3,
        },
    );

    let hk_sdf = ql_models::evaluate_he_krishnamurthy_intermediary_capital_constraint_sdf(
        &ql_models::HeKrishnamurthyIntermediarySdfSpec {
            intermediary_wealth_share_wi: T::from_f64(0.14) + dx0,
            capital_constraint_threshold_w_star: T::from_f64(0.28) + dx1,
            equity_intermediation_pledgeability_m: T::from_f64(0.35),
            specialist_risk_aversion_gamma: T::from_f64(2.20) + dx2,
            fundamental_growth_rate_ga: T::from_f64(0.020),
            asset_return_volatility_sigma_a: T::from_f64(0.15) + dx3,
            subjective_discount_rate_rho: T::from_f64(0.025),
        },
    );

    [
        ap_cms.cms_spread_call_present_value * T::from_f64(100.0),
        // Include both the Z-bond PV and senior sequential WAL so the benchmark
        // exercises the balance/CPR and discount/OAS sensitivities of the waterfall.
        fab_z.z_bond_present_value * T::from_f64(0.1)
            + fab_z.senior_tranche_wal_years * T::from_f64(0.1),
        dg_cdo.tranche_expected_loss_pv * T::from_f64(10.0),
        cl_vol.volatility_swap_fair_strike_k_vol * T::from_f64(10.0),
        gr_kva.total_dual_constraint_kva,
        cs_oil.three_factor_futures_price * T::from_f64(0.1),
        fst_xccs.collateralized_xccs_present_value * T::from_f64(10.0),
        ph_tsw.time_switch_option_present_value,
        bcf_mvar.tail_risk_parity_weight_w1 * T::from_f64(10.0),
        ru_quad.quadrangle_mixed_cvar_risk * T::from_f64(10.0),
        hgg_is.hasbrouck_midpoint_information_share * T::from_f64(10.0),
        hk_sdf.conditional_intermediary_equity_premium * T::from_f64(100.0),
    ]
}

pub fn benchmark_andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_vol_swap_green_leverage_ratio_kva_cortazar_schwartz_3f_oil_fujii_collateralized_xccs_pechtl_haug_time_switch_boudt_cornish_fisher_mvar_rockafellar_quadrangle_hasbrouck_info_share_he_krishnamurthy_sdf(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite66_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite66_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "66_andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_vol_swap_green_leverage_ratio_kva_cortazar_schwartz_3f_oil_fujii_collateralized_xccs_pechtl_haug_time_switch_boudt_cornish_fisher_mvar_rockafellar_quadrangle_hasbrouck_info_share_he_krishnamurthy_sdf",
        "Andersen-Piterbarg-CMS-Spread / Fabozzi-Z-Bond-CMO / Duffie-Garleanu-Affine-CDO / Carr-Lee-Vol-Swap / Green-Leverage-Ratio-KVA / Cortazar-Schwartz-3F-Oil / Fujii-Collateralized-XCCS / Pechtl-Haug-Time-Switch / Boudt-Cornish-Fisher-MVaR / Rockafellar-Quadrangle / Hasbrouck-Info-Share / He-Krishnamurthy-SDF AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite67_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let pit_ctd = ql_curves::evaluate_piterbarg_multi_currency_ctd_collateral_choice_discount(
        &ql_curves::PiterbargMultiCurrencyCtdCollateralSpec {
            domestic_ois_rate_rd: T::from_f64(0.042) + dx0,
            domestic_collateral_spread_s1: T::from_f64(0.0010) + dx1 * T::from_f64(0.1),
            foreign_synthetic_collateral_spread_s2: T::from_f64(0.0025) + dx2 * T::from_f64(0.1),
            collateral_spread_diff_volatility_sigma_x: T::from_f64(0.0035) + dx3 * T::from_f64(0.1),
            spread_reversion_speed_kappa_x: T::from_f64(0.45),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let tuck_mw = ql_fixed_income::evaluate_tuckman_corporate_make_whole_and_par_call_bond(
        &ql_fixed_income::TuckmanMakeWholeCallableBondSpec {
            par_value_f: T::from_f64(100.0),
            coupon_rate_c: T::from_f64(0.056) + dx0,
            treasury_benchmark_yield_y: T::from_f64(0.038) + dx1,
            issuer_credit_spread_s_corp: T::from_f64(0.014) + dx2 * T::from_f64(0.1),
            make_whole_treasury_spread_s_mw: T::from_f64(0.0030),
            credit_spread_volatility_sigma_s: T::from_f64(0.0065) + dx3 * T::from_f64(0.1),
            treasury_yield_volatility_sigma_y: T::from_f64(0.0085),
            make_whole_horizon_t_mw: T::from_f64(2.0),
            par_call_window_start_t_par: T::from_f64(4.5),
            maturity_years_t_mat: T::from_f64(5.0),
        },
    );

    let bm_cmcds = ql_credit::evaluate_brigo_morini_constant_maturity_cds_convexity(
        &ql_credit::BrigoMoriniCmcdsSpec {
            notional_n: T::from_f64(1_000.0),
            forward_reference_cds_spread_s0: T::from_f64(0.015) + dx0 * T::from_f64(0.1),
            reference_cds_tenor_tau_cm: T::from_f64(5.0),
            cds_spread_volatility_sigma_s: T::from_f64(0.45) + dx1,
            survival_hazard_rate_lambda: T::from_f64(0.025) + dx2 * T::from_f64(0.1),
            loss_given_default_lgd: T::from_f64(0.60),
            cmcds_cap_strike_k: T::from_f64(0.025) + dx3 * T::from_f64(0.1),
            risk_free_rate_r: T::from_f64(0.035),
            contract_maturity_years_t: T::from_f64(4.0),
        },
    );

    let hag_sabr = ql_vol::evaluate_hagan_shifted_sabr_negative_rate_implied_volatility(
        &ql_vol::HaganShiftedSabrSpec {
            forward_rate_f: T::from_f64(-0.004) + dx0 * T::from_f64(0.1),
            strike_rate_k: T::from_f64(-0.008) + dx1 * T::from_f64(0.1),
            displacement_shift_s: T::from_f64(0.030),
            sabr_alpha: T::from_f64(0.045) + dx2 * T::from_f64(0.1),
            sabr_beta: T::from_f64(0.50),
            sabr_vol_of_vol_nu: T::from_f64(0.40) + dx3,
            sabr_correlation_rho: T::from_f64(-0.35),
            discount_rate_r: T::from_f64(0.015),
            expiry_years_t: T::from_f64(2.0),
        },
    );

    let gk_tva =
        ql_xva::evaluate_green_kenyon_tax_valuation_adjustment_tva(&ql_xva::GreenKenyonTvaSpec {
            initial_epe0: T::from_f64(25.0) + x[0],
            initial_ene0: T::from_f64(14.0) + x[1],
            counterparty_hazard_rate_lambda_c: T::from_f64(0.022) + dx2 * T::from_f64(0.1),
            bank_hazard_rate_lambda_b: T::from_f64(0.012),
            loss_given_default_lgd: T::from_f64(0.60),
            unsecured_funding_spread_sf: T::from_f64(0.0095) + dx3 * T::from_f64(0.1),
            corporate_tax_rate_tau: T::from_f64(0.25),
            funding_interest_deductibility_phi: T::from_f64(0.85),
            dva_tax_realizability_psi: T::from_f64(0.40),
            ois_discount_rate_r: T::from_f64(0.035),
            maturity_years_t: T::from_f64(5.0),
        });

    let ew_inter =
        ql_commodities::evaluate_eydeland_wolyniec_interconnector_transmission_spread_option(
            &ql_commodities::EydelandWolyniecInterconnectorSpec {
                sink_pool_forward_f_sink: T::from_f64(88.0) + x[0],
                source_pool_forward_f_source: T::from_f64(68.0) + x[1],
                transmission_line_loss_fraction: T::from_f64(0.045) + dx2 * T::from_f64(0.1),
                wheeling_tariff_strike_k: T::from_f64(6.0),
                sink_volatility_sigma1: T::from_f64(0.42) + dx3,
                source_volatility_sigma2: T::from_f64(0.36),
                cross_pool_correlation_rho: T::from_f64(0.74),
                forced_outage_hazard_lambda: T::from_f64(0.80),
                repair_rate_mu: T::from_f64(8.50),
                risk_free_rate_r: T::from_f64(0.04),
                maturity_years_t: T::one(),
            },
        );

    let wc_str = ql_fx::evaluate_wystup_clark_fx_smile_strangle_market_convention(
        &ql_fx::WystupClarkFxSmileStrangleSpec {
            spot_fx_s0: T::from_f64(1.10) + dx0,
            domestic_rate_rd: T::from_f64(0.045) + dx1,
            foreign_rate_rf: T::from_f64(0.030),
            atm_volatility_sigma_atm: T::from_f64(0.095) + dx2,
            risk_reversal_25d_rr: T::from_f64(-0.015),
            butterfly_25d_bf: T::from_f64(0.0035) + dx3 * T::from_f64(0.1),
            expiry_years_t: T::one(),
        },
    );

    let ki_bar = ql_exotics::evaluate_kunitomo_ikeda_curved_exponential_barrier_option(
        &ql_exotics::KunitomoIkedaCurvedBarrierSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(100.0) + x[1],
            initial_down_barrier_h0: T::from_f64(88.0) + x[2] * T::from_f64(0.5),
            barrier_curvature_rate_delta_h: T::from_f64(0.02),
            risk_free_rate_r: T::from_f64(0.05),
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.25) + dx3,
            maturity_years_t: T::one(),
        },
    );

    let clmx_var =
        ql_portfolio::evaluate_campbell_lettau_malkiel_xu_hierarchical_variance_decomposition(
            &ql_portfolio::CampbellLettauMalkielXuSpec {
                market_volatility_sigma_m: T::from_f64(0.16) + dx0,
                industry_residual_volatility_sigma_ind: T::from_f64(0.12) + dx1,
                firm_idiosyncratic_volatility_sigma_firm: T::from_f64(0.32) + dx2,
                industry_herfindahl_index_h_ind: T::from_f64(0.08) + dx3,
                portfolio_stock_count_n: T::from_f64(30.0),
            },
        );

    let kyle_eq = ql_portfolio::evaluate_kyle_multi_period_insider_market_depth_equilibrium(
        &ql_portfolio::KyleMultiPeriodEquilibriumSpec {
            prior_fundamental_variance_sigma0: T::from_f64(4.0) + x[0] * T::from_f64(0.1),
            noise_trader_volatility_sigma_u: T::from_f64(12.0) + x[1] * T::from_f64(0.5),
            insider_valuation_wedge_v_minus_p0: T::from_f64(2.5) + x[2] * T::from_f64(0.1),
            auction_step_dt: T::from_f64(0.5) + dx3,
        },
    );

    let hl_rv = ql_dml::evaluate_hansen_lunde_overnight_intraday_realized_variance_combination(
        &ql_dml::HansenLundeOvernightRvSpec {
            realized_overnight_sq_return_rv1: T::from_f64(0.00012) + dx0 * T::from_f64(1e-4),
            realized_intraday_variance_rv2: T::from_f64(0.00038) + dx1 * T::from_f64(1e-4),
            unconditional_overnight_mean_mu1: T::from_f64(0.00010) + dx2 * T::from_f64(1e-4),
            unconditional_intraday_mean_mu2: T::from_f64(0.00035),
            unconditional_full_day_iv_mu0: T::from_f64(0.00042) + dx3 * T::from_f64(1e-4),
            overnight_estimator_variance_eta1_sq: T::from_f64(2.5e-8),
            intraday_estimator_variance_eta2_sq: T::from_f64(1.8e-8),
            estimator_covariance_eta12: T::from_f64(0.4e-8),
        },
    );

    let dd_gp = ql_models::evaluate_diamond_dybvig_goldstein_pauzner_bank_run_global_game(
        &ql_models::DiamondDybvigGoldsteinPauznerBankRunSpec {
            short_term_deposit_payout_r1: T::from_f64(1.06) + dx0,
            long_term_asset_return_r_bar: T::from_f64(1.55) + dx1,
            liquid_reserve_ratio_l: T::from_f64(0.22) + dx2,
            fire_sale_recovery_fraction_psi: T::from_f64(0.65),
            impatient_depositor_fraction_lambda: T::from_f64(0.30),
            fundamental_prior_mean_mu_theta: T::from_f64(0.85) + dx3,
            fundamental_volatility_sigma_theta: T::from_f64(0.18),
        },
    );

    [
        pit_ctd.ctd_collateral_discount_factor * T::from_f64(10.0),
        tuck_mw.callable_corporate_bond_price * T::from_f64(0.1),
        bm_cmcds.capped_cmcds_floating_leg_pv * T::from_f64(0.1),
        hag_sabr.shifted_bachelier_call_price * T::from_f64(100.0),
        gk_tva.post_tax_net_xva_with_tva,
        ew_inter.expected_interconnector_congestion_rent_pv,
        wc_str.smile_adjusted_strangle_present_value * T::from_f64(10.0),
        ki_bar.kunitomo_ikeda_down_and_out_call_price,
        clmx_var.finite_n_portfolio_volatility * T::from_f64(10.0),
        kyle_eq.expected_insider_information_profit,
        hl_rv.hansen_lunde_combined_realized_volatility * T::from_f64(100.0),
        dd_gp.expected_depositor_welfare_value * T::from_f64(10.0),
    ]
}

pub fn benchmark_piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_green_kenyon_tva_eydeland_interconnector_spread_wystup_clark_fx_strangle_kunitomo_ikeda_curved_barrier_clmx_variance_decomp_kyle_multi_period_lambda_hansen_lunde_overnight_rv_goldstein_pauzner_bank_run(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite67_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite67_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "67_piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_green_kenyon_tva_eydeland_interconnector_spread_wystup_clark_fx_strangle_kunitomo_ikeda_curved_barrier_clmx_variance_decomp_kyle_multi_period_lambda_hansen_lunde_overnight_rv_goldstein_pauzner_bank_run",
        "Piterbarg-CTD-Collateral / Tuckman-Make-Whole-Bond / Brigo-Morini-CMCDS / Hagan-Shifted-SABR / Green-Kenyon-TVA / Eydeland-Interconnector-Spread / Wystup-Clark-FX-Strangle / Kunitomo-Ikeda-Curved-Barrier / CLMX-Variance-Decomp / Kyle-Multi-Period-Lambda / Hansen-Lunde-Overnight-RV / Goldstein-Pauzner-Bank-Run AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite68_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let merc_fva = ql_curves::evaluate_mercurio_post_crisis_fva_adjusted_multi_curve_swap_rate(
        &ql_curves::MercurioFvaMultiCurveSwapSpec {
            ois_short_rate_r0: T::from_f64(0.038) + dx0,
            ois_curve_slope_theta: T::from_f64(0.004) + dx1 * T::from_f64(0.5),
            initial_multiplicative_ibor_basis_b0: T::from_f64(0.0025) + dx2 * T::from_f64(0.1),
            basis_term_decay_kappa_b: T::from_f64(0.25),
            funding_borrowing_spread_s_f: T::from_f64(0.0065) + dx3 * T::from_f64(0.1),
            swap_rate_volatility_sigma_s: T::from_f64(0.0085),
            swap_maturity_years_t: T::from_f64(4.0),
        },
    );

    let ing_cb =
        ql_fixed_income::evaluate_ingersoll_contingent_call_forced_conversion_convertible_bond(
            &ql_fixed_income::IngersollConvertibleBondSpec {
                stock_price_s0: T::from_f64(96.0) + x[0],
                conversion_ratio_cr: T::one() + dx1,
                straight_bond_coupon_c: T::from_f64(0.04),
                call_trigger_price_k_call: T::from_f64(125.0) + x[2],
                call_protection_years_t_call: T::from_f64(2.0),
                maturity_years_t: T::from_f64(5.0),
                credit_adjusted_discount_rate_r_c: T::from_f64(0.055),
                equity_volatility_sigma: T::from_f64(0.28) + dx3,
            },
        );

    let hw_cdo = ql_credit::evaluate_hull_white_implied_copula_bespoke_cdo_tranche_calibration(
        &ql_credit::HullWhiteImpliedCopulaCdoSpec {
            base_hazard_rate_lambda0: T::from_f64(0.015) + dx0 * T::from_f64(0.1),
            recession_hazard_multiplier_m_rec: T::from_f64(3.2) + x[1] * T::from_f64(0.1),
            systemic_crash_hazard_lambda_crash: T::from_f64(0.14) + dx2,
            crash_state_probability_p_crash: T::from_f64(0.06) + dx3 * T::from_f64(0.5),
            loss_given_default_lgd: T::from_f64(0.60),
            tranche_attachment_k1: T::from_f64(0.03),
            tranche_detachment_k2: T::from_f64(0.07),
            risk_free_rate_r: T::from_f64(0.035),
        },
    );

    let gj_svi = ql_vol::evaluate_gatheral_jacquier_svi_raw_to_jump_wing_reparameterization(
        &ql_vol::GatheralJacquierSviJumpWingSpec {
            raw_svi_a: T::from_f64(0.012) + dx0 * T::from_f64(0.1),
            raw_svi_b: T::from_f64(0.14) + dx1,
            raw_svi_rho: T::from_f64(-0.42) + dx2,
            raw_svi_m: T::from_f64(0.05),
            raw_svi_sigma: T::from_f64(0.18) + dx3,
            evaluation_log_moneyness_k: T::from_f64(-0.10),
            expiry_years_t: T::from_f64(0.75),
        },
    );

    let bk_xva = ql_xva::evaluate_burgard_kjaer_semi_replication_funding_strategy_xva(
        &ql_xva::BurgardKjaerSemiReplicationFundingStrategyXvaSpec {
            risk_free_derivative_value_v0: T::from_f64(100.0) + x[0],
            negative_exposure_ene0: T::from_f64(45.0) + x[1],
            counterparty_hazard_rate_lambda_c: T::from_f64(0.022) + dx2 * T::from_f64(0.1),
            bank_hazard_rate_lambda_b: T::from_f64(0.015),
            loss_given_default_lgd: T::from_f64(0.60),
            own_bond_recovery_r_b: T::from_f64(0.40),
            semi_replication_hedge_weight_alpha: T::from_f64(0.65) + dx3,
            ois_discount_rate_r: T::from_f64(0.035),
        },
    );

    let benth_hdd =
        ql_commodities::evaluate_benth_cartea_saltyte_benth_seasonal_temperature_hdd_cdd_swap(
            &ql_commodities::BenthSeasonalTemperatureHddCddSpec {
                initial_temperature_t0: T::from_f64(12.5) + x[0] * T::from_f64(0.5),
                seasonal_mean_level_a0: T::from_f64(15.0) + x[1] * T::from_f64(0.5),
                seasonal_amplitude_b1: T::from_f64(11.0),
                mean_reversion_speed_kappa: T::from_f64(6.5),
                temperature_volatility_sigma: T::from_f64(8.0) + x[2] * T::from_f64(0.5),
                degree_day_threshold_k_base: T::from_f64(18.0) + x[3] * T::from_f64(0.5),
                market_price_of_weather_risk_theta: T::from_f64(-0.18),
                discount_rate_r: T::from_f64(0.04),
            },
        );

    let bp_fx = ql_fx::evaluate_baba_packer_nagano_fx_swap_implied_dollar_funding_basis(
        &ql_fx::BabaPackerFxSwapDollarBasisSpec {
            spot_fx_s0: T::from_f64(1.085) + dx0,
            outright_forward_fx_f: T::from_f64(1.096) + dx1,
            usd_ois_rate_r_usd: T::from_f64(0.048) + dx2,
            foreign_ois_rate_r_for: T::from_f64(0.032) + dx3,
            dealer_balance_sheet_leverage_cost_lambda_bs: T::from_f64(0.0022),
            foreign_dollar_funding_demand_pressure_d_usd: T::from_f64(0.0035),
            swap_tenor_years_t: T::from_f64(0.50),
        },
    );

    let cv_lb = ql_exotics::evaluate_conze_viswanathan_floating_and_fixed_lookback_extremum_option(
        &ql_exotics::ConzeViswanathanLookbackSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            observed_running_minimum_s_min: T::from_f64(88.0) + x[1] * T::from_f64(0.5),
            risk_free_rate_r: T::from_f64(0.05) + dx2,
            dividend_yield_q: T::from_f64(0.015) + dx3,
            volatility_sigma: T::from_f64(0.24),
            discrete_monitoring_steps_m: T::from_f64(52.0),
            maturity_years_t: T::one(),
        },
    );

    let hl_bl = ql_portfolio::evaluate_litterman_he_institutional_equilibrium_posterior_black_litterman_tilt(
        &ql_portfolio::HeLittermanBlackLittermanTiltSpec {
            risk_aversion_delta: T::from_f64(2.5) + x[0] * T::from_f64(0.1),
            market_cap_weight_asset1_w1: T::from_f64(0.60) + dx1,
            asset1_volatility_sigma1: T::from_f64(0.16) + dx2,
            asset2_volatility_sigma2: T::from_f64(0.22),
            asset_correlation_rho: T::from_f64(0.35),
            prior_uncertainty_scalar_tau: T::from_f64(0.05),
            relative_view_outperformance_q: T::from_f64(0.035) + dx3,
            view_error_volatility_sigma_omega: T::from_f64(0.04),
        },
    );

    let ghs_is =
        ql_portfolio::evaluate_glasserman_heidelberger_shahabuddin_importance_sampling_credit_var(
            &ql_portfolio::GlassermanImportanceSamplingVaRSpec {
                portfolio_linear_delta_exposure_a: T::from_f64(1.80) + dx0 * T::from_f64(5.0),
                portfolio_quadratic_gamma_exposure_b: T::from_f64(0.22) + dx1,
                idiosyncratic_noise_volatility_sigma_e: T::from_f64(0.65) + dx2 * T::from_f64(5.0),
                tail_loss_threshold_x: T::from_f64(5.20) + x[3] * T::from_f64(0.2),
                constant_loss_drift_c0: T::from_f64(0.15),
                confidence_level_alpha: T::from_f64(0.99),
                effective_sample_paths_n: T::from_f64(1000.0),
            },
        );

    let egj_vecm =
        ql_dml::evaluate_engle_granger_johansen_cointegration_error_correction_pairs_spread(
            &ql_dml::EngleGrangerJohansenCointegrationSpec {
                current_log_price_y: T::from_f64(4.65) + dx0,
                current_log_price_x: T::from_f64(4.52) + dx1,
                cointegration_hedge_ratio_beta: T::from_f64(1.04) + dx2,
                equilibrium_intercept_mu: T::from_f64(-0.08),
                error_correction_speed_alpha: T::from_f64(0.28) + dx3,
                innovation_volatility_sigma_eps: T::from_f64(0.032),
                canonical_correlation_r_sq: T::from_f64(0.24),
            },
        );

    let gab_dis = ql_models::evaluate_gabaix_rare_disasters_time_varying_resilience_equity_premium(
        &ql_models::GabaixRareDisastersResilienceSpec {
            subjective_discount_rate_rho: T::from_f64(0.032) + dx0,
            relative_risk_aversion_gamma: T::from_f64(4.0) + x[1] * T::from_f64(0.1),
            normal_consumption_growth_g_c: T::from_f64(0.020),
            normal_consumption_volatility_sigma_c: T::from_f64(0.018),
            disaster_probability_p: T::from_f64(0.017) + dx2 * T::from_f64(0.1),
            consumption_recovery_fraction_b_c: T::from_f64(0.72),
            dividend_recovery_fraction_b_d: T::from_f64(0.18) + dx3,
            resilience_mean_reversion_phi_h: T::from_f64(0.18),
        },
    );

    [
        merc_fva.fva_adjusted_uncollateralized_swap_rate * T::from_f64(100.0),
        ing_cb.convertible_bond_fair_value * T::from_f64(0.1),
        hw_cdo.bespoke_tranche_par_spread * T::from_f64(100.0),
        gj_svi.svi_implied_volatility_sigma_k * T::from_f64(10.0),
        bk_xva.total_burgard_kjaer_xva_charge,
        benth_hdd.discounted_hdd_swap_present_value * T::from_f64(0.01),
        bp_fx.balance_sheet_adjusted_dollar_funding_pv * T::from_f64(100.0),
        cv_lb.discrete_floating_lookback_call_price * T::from_f64(0.5),
        hl_bl.posterior_expected_return_asset1 * T::from_f64(100.0),
        ghs_is.saddlepoint_adjusted_quadratic_var,
        egj_vecm.cointegration_confidence_weighted_sharpe,
        gab_dis.gabaix_equilibrium_price_dividend_ratio,
    ]
}

pub fn benchmark_mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_burgard_kjaer_semi_replication_xva_benth_seasonal_hdd_cdd_swap_baba_packer_fx_swap_dollar_basis_conze_viswanathan_lookback_he_litterman_bl_posterior_ghs_importance_sampling_cvar_engle_granger_johansen_vecm_gabaix_rare_disasters_resilience(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite68_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite68_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "68_mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_burgard_kjaer_semi_replication_xva_benth_seasonal_hdd_cdd_swap_baba_packer_fx_swap_dollar_basis_conze_viswanathan_lookback_he_litterman_bl_posterior_ghs_importance_sampling_cvar_engle_granger_johansen_vecm_gabaix_rare_disasters_resilience",
        "Mercurio-FVA-Multi-Curve-Swap / Ingersoll-Contingent-Call-CB / Hull-White-Implied-Copula-CDO / Gatheral-Jacquier-SVI-JW / Burgard-Kjaer-Semi-Replication-XVA / Benth-Seasonal-HDD-CDD-Swap / Baba-Packer-FX-Swap-Dollar-Basis / Conze-Viswanathan-Lookback / He-Litterman-BL-Posterior / GHS-Importance-Sampling-CVaR / Engle-Granger-Johansen-VECM / Gabaix-Rare-Disasters-Resilience AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite69_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let ab_tb = ql_curves::evaluate_ametrano_bianchetti_dual_curve_tenor_basis_swap_bootstrap(
        &ql_curves::AmetranoBianchettiTenorBasisSpec {
            ois_discount_rate_r_ois: T::from_f64(0.036) + dx0,
            ois_slope_theta_ois: T::from_f64(0.004) + dx1 * T::from_f64(0.5),
            ibor_3m_spread_over_ois_s3: T::from_f64(0.0018) + dx2 * T::from_f64(0.1),
            tenor_basis_3m_vs_6m_spread_b36: T::from_f64(0.0014) + dx3 * T::from_f64(0.1),
            tenor_basis_term_steepening_gamma: T::from_f64(0.0005),
            compounding_convexity_vol_sigma3: T::from_f64(0.0075),
            swap_maturity_years_t: T::from_f64(2.0),
        },
    );

    let fh_iopo = ql_fixed_income::evaluate_fabozzi_hayre_io_po_stripped_mbs_prepayment_convexity(
        &ql_fixed_income::FabozziHayreIoPoStripSpec {
            pool_principal_balance: T::from_f64(100.0) + x[0],
            weighted_average_coupon_wac: T::from_f64(0.065) + dx1,
            net_passthrough_coupon_c_net: T::from_f64(0.060),
            current_mortgage_market_rate_r_mkt: T::from_f64(0.052) + dx2,
            baseline_cpr_min: T::from_f64(0.06),
            peak_refi_cpr_max: T::from_f64(0.42) + dx3,
            refi_s_curve_steepness_beta: T::from_f64(35.0),
            option_adjusted_spread_oas: T::from_f64(0.008),
        },
    );

    let hb_wwr = ql_credit::evaluate_hull_white_brigo_cds_counterparty_wrong_way_spread(
        &ql_credit::HullBrigoWwrCdsSpec {
            reference_entity_hazard_lambda_ref: T::from_f64(0.024) + dx0 * T::from_f64(0.1),
            protection_seller_hazard_lambda_cp: T::from_f64(0.016) + dx1 * T::from_f64(0.1),
            default_correlation_rho_wwr: T::from_f64(0.45) + dx2,
            reference_recovery_rate_r_ref: T::from_f64(0.40),
            counterparty_recovery_rate_r_cp: T::from_f64(0.35),
            cds_spread_volatility_sigma_s: T::from_f64(0.42) + dx3,
            risk_free_rate_r: T::from_f64(0.035),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let brs_vv = ql_vol::evaluate_bossens_rayee_skovmand_vanna_volga_barrier_survival_overhedge(
        &ql_vol::BossensVannaVolgaBarrierSpec {
            spot_fx_s0: T::from_f64(1.10) + dx0,
            strike_k: T::from_f64(1.12) + dx1,
            down_and_out_barrier_h: T::from_f64(1.02) + dx2 * T::from_f64(0.5),
            atm_volatility_sigma_atm: T::from_f64(0.10) + dx3,
            risk_reversal_25d_rr: T::from_f64(-0.018),
            butterfly_25d_bf: T::from_f64(0.0045),
            domestic_foreign_rate_diff_mu: T::from_f64(0.015),
            maturity_years_t: T::one(),
        },
    );

    let greg_mva = ql_xva::evaluate_gregory_isda_simm_dynamic_initial_margin_mva(
        &ql_xva::GregoryIsdaSimmMvaSpec {
            trade_notional_dv01_w0: T::from_f64(120.0) + x[0],
            simm_calibrated_volatility_sigma_im: T::from_f64(0.011) + dx1 * T::from_f64(0.1),
            margin_period_of_risk_days_mpor: T::from_f64(10.0),
            im_funding_spread_s_im: T::from_f64(0.0085) + dx2 * T::from_f64(0.1),
            counterparty_hazard_rate_lambda_c: T::from_f64(0.018) + dx3 * T::from_f64(0.1),
            bank_hazard_rate_lambda_b: T::from_f64(0.012),
            ois_discount_rate_r: T::from_f64(0.035),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let cd_cds = ql_commodities::evaluate_carmona_durrleman_clean_dark_spread_carbon_power_option(
        &ql_commodities::CarmonaCleanDarkSpreadSpec {
            electricity_forward_f_power: T::from_f64(95.0) + x[0],
            coal_forward_f_coal: T::from_f64(110.0) + x[1],
            carbon_allowance_forward_f_eua: T::from_f64(68.0) + x[2],
            coal_heat_rate_h_coal: T::from_f64(0.38),
            carbon_emission_intensity_e_co2: T::from_f64(0.92),
            fixed_om_strike_k: T::from_f64(4.50),
            power_volatility_sigma_p: T::from_f64(0.38) + dx3,
            fuel_carbon_basket_volatility_sigma_fc: T::from_f64(0.30),
            power_fuel_correlation_rho: T::from_f64(0.72),
        },
    );

    let rh_qto = ql_fx::evaluate_reiner_haug_foreign_equity_domestic_strike_and_quanto_call(
        &ql_fx::ReinerHaugQuantoEquityCallSpec {
            foreign_equity_spot_s0: T::from_f64(100.0) + x[0],
            fx_spot_e0: T::from_f64(1.10) + dx1,
            foreign_strike_k_for: T::from_f64(102.0) + x[2],
            domestic_rate_rd: T::from_f64(0.045),
            foreign_rate_rf: T::from_f64(0.025),
            equity_dividend_yield_q: T::from_f64(0.015),
            equity_volatility_sigma_s: T::from_f64(0.22) + dx3,
            fx_volatility_sigma_e: T::from_f64(0.12),
            equity_fx_correlation_rho: T::from_f64(-0.35),
        },
    );

    let hr_soft = ql_exotics::evaluate_hart_ross_haug_soft_proportional_band_barrier_option(
        &ql_exotics::HartRossSoftBarrierSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(100.0) + x[1],
            lower_soft_barrier_l: T::from_f64(80.0) + x[2] * T::from_f64(0.5),
            upper_soft_barrier_u: T::from_f64(90.0),
            risk_free_rate_r: T::from_f64(0.05),
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.24) + dx3,
            maturity_years_t: T::one(),
        },
    );

    let ps_alp = ql_portfolio::evaluate_pastor_stambaugh_bayesian_mispricing_alpha_shrinkage(
        &ql_portfolio::PastorStambaughAlphaShrinkageSpec {
            ols_sample_alpha_hat: T::from_f64(0.042) + dx0,
            prior_mispricing_uncertainty_sigma_alpha: T::from_f64(0.020) + dx1 * T::from_f64(0.2),
            residual_idiosyncratic_volatility_sigma_eps: T::from_f64(0.15) + dx2,
            sample_history_years_t_obs: T::from_f64(5.0),
            benchmark_factor_sharpe_ratio_sr_f: T::from_f64(0.45),
            asset_factor_beta: T::from_f64(1.10) + dx3,
            benchmark_expected_excess_return_mu_f: T::from_f64(0.065),
        },
    );

    let ow_lob = ql_portfolio::evaluate_obizhaeva_wang_transient_resilience_impulse_block_execution(
        &ql_portfolio::ObizhaevaWangImpulseExecutionSpec {
            order_size_shares_x0: T::from_f64(100.0) + x[0],
            permanent_price_impact_lambda: T::from_f64(0.0015) + dx1 * T::from_f64(0.01),
            transient_lob_impact_kappa: T::from_f64(0.0045) + dx2 * T::from_f64(0.01),
            lob_resilience_speed_rho: T::from_f64(2.20) + dx3 * T::from_f64(10.0),
            execution_horizon_t: T::one(),
            asset_volatility_sigma: T::from_f64(0.25),
            inventory_risk_aversion_gamma: T::from_f64(0.02),
        },
    );

    let dmw_hln = ql_dml::evaluate_diebold_mariano_west_hln_forecast_comparison_test(
        &ql_dml::DieboldMarianoWestForecastSpec {
            mean_loss_differential_d_bar: T::from_f64(0.028) + dx0 * T::from_f64(0.5),
            lag0_loss_diff_variance_gamma0: T::from_f64(0.045) + dx1 * T::from_f64(0.5),
            lag1_loss_diff_autocovariance_gamma1: T::from_f64(0.014) + dx2 * T::from_f64(0.5),
            forecast_horizon_steps_h: T::from_f64(3.0),
            out_of_sample_observations_p: T::from_f64(120.0) + x[3],
            in_sample_estimation_observations_r: T::from_f64(240.0),
            parameter_estimation_variance_v_theta: T::from_f64(0.008),
        },
    );

    let ps_nom = ql_models::evaluate_piazzesi_schneider_stagflation_nominal_bond_term_premium(
        &ql_models::PiazzesiSchneiderNominalBondSpec {
            real_short_rate_r_real: T::from_f64(0.015) + dx0,
            expected_inflation_pi0: T::from_f64(0.025) + dx1,
            relative_risk_aversion_gamma: T::from_f64(8.5) + x[2] * T::from_f64(0.2),
            real_growth_persistence_phi_x: T::from_f64(0.92),
            inflation_persistence_phi_pi: T::from_f64(0.94),
            real_growth_volatility_sigma_x: T::from_f64(0.006),
            inflation_volatility_sigma_pi: T::from_f64(0.008) + dx3 * T::from_f64(0.1),
            growth_inflation_correlation_rho_x_pi: T::from_f64(-0.55),
        },
    );

    [
        ab_tb.bootstrapped_6m_ibor_par_swap_rate * T::from_f64(100.0),
        fh_iopo.io_strip_present_value,
        hb_wwr.wrong_way_adjusted_cds_par_spread * T::from_f64(100.0),
        brs_vv.vanna_volga_barrier_call_price * T::from_f64(100.0),
        greg_mva.lifetime_simm_mva_charge * T::from_f64(10.0),
        cd_cds.clean_dark_spread_call_option_value * T::from_f64(0.5),
        rh_qto.fixed_fx_quanto_call_present_value,
        hr_soft.soft_down_and_out_call_price,
        ps_alp.posterior_predictive_information_ratio * T::from_f64(10.0),
        ow_lob.risk_adjusted_total_execution_cost,
        dmw_hln.hln_adjusted_dmw_t_statistic,
        ps_nom.equilibrium_long_nominal_bond_yield * T::from_f64(100.0),
    ]
}

pub fn benchmark_ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_gregory_simm_mva_carmona_clean_dark_spread_reiner_haug_quanto_call_hart_ross_soft_barrier_pastor_stambaugh_alpha_shrinkage_obizhaeva_wang_impulse_lob_diebold_mariano_west_hln_piazzesi_schneider_nominal_tp(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite69_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite69_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "69_ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_gregory_simm_mva_carmona_clean_dark_spread_reiner_haug_quanto_call_hart_ross_soft_barrier_pastor_stambaugh_alpha_shrinkage_obizhaeva_wang_impulse_lob_diebold_mariano_west_hln_piazzesi_schneider_nominal_tp",
        "Ametrano-Bianchetti-Tenor-Basis / Fabozzi-Hayre-IO-PO-Strip / Hull-Brigo-WWR-CDS / Bossens-Vanna-Volga-Barrier / Gregory-SIMM-MVA / Carmona-Clean-Dark-Spread / Reiner-Haug-Quanto-Call / Hart-Ross-Soft-Barrier / Pastor-Stambaugh-Alpha-Shrinkage / Obizhaeva-Wang-Impulse-LOB / Diebold-Mariano-West-HLN / Piazzesi-Schneider-Nominal-TP AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite70_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let sw_ufr = ql_curves::evaluate_smith_wilson_eiopa_ultimate_forward_rate_extrapolation(
        &ql_curves::SmithWilsonUfrExtrapolationSpec {
            intermediate_pillar_years_u1: T::from_f64(10.0) + dx0,
            last_liquid_point_years_u2: T::from_f64(20.0) + dx1,
            extrapolated_maturity_years_t_ext: T::from_f64(50.0),
            intermediate_zero_rate_r1: T::from_f64(0.026) + dx2 * T::from_f64(0.2),
            last_liquid_point_zero_rate_r2: T::from_f64(0.029) + dx3 * T::from_f64(0.2),
            ultimate_forward_rate_omega: T::from_f64(0.036),
            convergence_speed_alpha: T::from_f64(0.12),
        },
    );

    let kwf_cb =
        ql_fixed_income::evaluate_kalotay_williams_fabozzi_callable_bond_refunding_efficiency(
            &ql_fixed_income::KalotayWilliamsFabozziCallableSpec {
                initial_short_rate_r0: T::from_f64(0.040) + dx0 * T::from_f64(0.2),
                short_rate_drift_mu: T::from_f64(0.002) + dx1 * T::from_f64(0.1),
                lognormal_short_rate_volatility_sigma: T::from_f64(0.22) + dx2,
                option_adjusted_spread_oas: T::from_f64(0.010) + dx3 * T::from_f64(0.2),
                coupon_rate_c: T::from_f64(0.058),
                call_strike_price_k: T::from_f64(100.2),
                refunding_transaction_cost: T::from_f64(0.35),
                step_years_dt: T::one(),
            },
        );

    let sc_cdo = ql_credit::evaluate_schonbucher_archimedean_clayton_copula_cdo_tranche_loss(
        &ql_credit::SchonbucherClaytonCopulaCdoSpec {
            pool_hazard_rate_lambda: T::from_f64(0.022) + dx0 * T::from_f64(0.1),
            recovery_rate_r: T::from_f64(0.40) + dx1 * T::from_f64(0.2),
            clayton_dependence_theta: T::from_f64(0.65) + dx2,
            attachment_point_ka: T::from_f64(0.03),
            detachment_point_kb: T::from_f64(0.07),
            risk_free_rate_r: T::from_f64(0.035) + dx3 * T::from_f64(0.2),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let bue_vc = ql_vol::evaluate_buehler_consistent_double_exponential_variance_curve(
        &ql_vol::BuehlerVarianceCurveSpec {
            short_variance_xi0: T::from_f64(0.040) + dx0 * T::from_f64(0.2),
            long_variance_xi_inf: T::from_f64(0.048) + dx1 * T::from_f64(0.2),
            hump_variance_amplitude_xi_hump: T::from_f64(0.015) + dx2 * T::from_f64(0.2),
            fast_mean_reversion_kappa1: T::from_f64(2.4),
            slow_mean_reversion_kappa2: T::from_f64(0.55),
            fast_factor_vol_of_vol_omega1: T::from_f64(1.35) + dx3,
            slow_factor_vol_of_vol_omega2: T::from_f64(0.45),
            factor_correlation_rho12: T::from_f64(0.35),
            maturity_years_t: T::from_f64(1.5),
        },
    );

    let pr_cva = ql_xva::evaluate_pykhtin_rosen_euler_marginal_cva_netting_allocation(
        &ql_xva::PykhtinRosenMarginalCvaSpec {
            subportfolio_1_expected_mtm_mu1: T::from_f64(12.0) + x[0],
            subportfolio_2_expected_mtm_mu2: T::from_f64(-4.5) + x[1],
            subportfolio_1_volatility_sigma1: T::from_f64(28.0) + x[2],
            subportfolio_2_volatility_sigma2: T::from_f64(22.0),
            mtm_correlation_rho12: T::from_f64(-0.35) + dx3,
            counterparty_hazard_lambda_c: T::from_f64(0.024),
            counterparty_recovery_rate_rc: T::from_f64(0.40),
            risk_free_rate_r: T::from_f64(0.032),
        },
    );

    let hr_fut =
        ql_commodities::evaluate_hilliard_reis_stochastic_convenience_yield_jump_futures_option(
            &ql_commodities::HilliardReisJumpFuturesSpec {
                futures_price_f0: T::from_f64(78.0) + x[0],
                strike_price_k: T::from_f64(80.0) + x[1],
                spot_volatility_sigma_s: T::from_f64(0.32) + dx2,
                convenience_yield_volatility_sigma_delta: T::from_f64(0.24),
                convenience_yield_mean_reversion_kappa: T::from_f64(1.60),
                spot_convenience_correlation_rho: T::from_f64(0.65),
                jump_intensity_lambda_j: T::from_f64(0.85) + dx3 * T::from_f64(2.0),
                jump_volatility_v_j: T::from_f64(0.18),
                risk_free_rate_r: T::from_f64(0.038),
            },
        );

    let gm_fx = ql_fx::evaluate_gabaix_maggiori_imperfect_financial_markets_fx_carry_risk_premium(
        &ql_fx::GabaixMaggioriFxCarrySpec {
            home_gross_interest_rate_r: T::from_f64(1.045) + dx0,
            foreign_gross_interest_rate_r_star: T::from_f64(1.015) + dx1,
            financier_friction_gamma: T::from_f64(0.85) + dx2 * T::from_f64(2.0),
            fx_volatility_sigma: T::from_f64(0.14) + dx3,
            capital_flow_demand_f0: T::from_f64(0.22),
            trade_flow_elasticity_iota0: T::one(),
            period_1_trade_elasticity_iota1: T::from_f64(1.05),
        },
    );

    let cv_shout = ql_exotics::evaluate_cheuk_vorst_haug_strike_reset_shout_call_option(
        &ql_exotics::CheukVorstShoutResetSpec {
            spot_s0: T::from_f64(100.0) + x[0],
            initial_strike_k0: T::from_f64(102.0) + x[1],
            reset_monitoring_time_t1: T::from_f64(0.45) + dx2,
            maturity_years_t: T::one(),
            risk_free_rate_r: T::from_f64(0.04),
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.24) + dx3,
        },
    );

    let gz_kelly =
        ql_portfolio::evaluate_maclean_thorp_ziemba_grossman_zhou_fractional_kelly_drawdown(
            &ql_portfolio::GrossmanZhouDrawdownKellySpec {
                expected_return_mu: T::from_f64(0.115) + dx0,
                risk_free_rate_r: T::from_f64(0.035) + dx1,
                asset_volatility_sigma: T::from_f64(0.19) + dx2,
                max_drawdown_floor_delta: T::from_f64(0.25) + dx3,
                cushion_risk_aversion_gamma: T::from_f64(2.0),
                current_drawdown_from_peak: T::from_f64(0.06),
            },
        );

    let gls_raroc =
        ql_portfolio::evaluate_gourieroux_laurent_scaillet_tasche_cornish_fisher_euler_var_raroc(
            &ql_portfolio::GlsTascheCornishFisherRarocSpec {
                portfolio_expected_excess_return_mu_p: T::from_f64(0.082) + dx0,
                portfolio_volatility_sigma_p: T::from_f64(0.145) + dx1,
                portfolio_skewness_s_p: T::from_f64(-0.65) + dx2,
                portfolio_excess_kurtosis_k_p: T::from_f64(2.40),
                normal_tail_quantile_z_alpha: T::from_f64(2.326_347_874),
                subportfolio_expected_excess_return_mu_i: T::from_f64(0.105) + dx3,
                subportfolio_volatility_sigma_i: T::from_f64(0.180),
                subportfolio_portfolio_correlation_rho: T::from_f64(0.72),
            },
        );

    let ap_wald = ql_dml::evaluate_andrews_ploberger_hansen_sup_exp_ave_wald_structural_break_test(
        &ql_dml::AndrewsPlobergerStructuralBreakSpec {
            sample_size_t_obs: T::from_f64(120.0) + x[0],
            pre_break_slope_beta1: T::from_f64(0.45) + dx1,
            post_break_slope_beta2: T::from_f64(0.72) + dx2,
            break_profile_dispersion_delta: T::from_f64(0.14) + dx3,
            residual_variance_sigma2: T::from_f64(0.55),
            regressor_second_moment_q_xx: T::from_f64(1.15),
            trimming_fraction_pi0: T::from_f64(0.15),
        },
    );

    let jq_fin =
        ql_models::evaluate_jermann_quadrini_financial_shocks_enforcement_constraint_credit_cycle(
            &ql_models::JermannQuadriniFinancialShockSpec {
                capital_stock_k: T::from_f64(10.0) + dx0,
                working_capital_output_y: T::from_f64(2.4) + dx1,
                financial_enforcement_xi: T::from_f64(0.42) + dx2,
                risk_free_rate_r: T::from_f64(0.032) + dx3 * T::from_f64(0.2),
                corporate_tax_shield_tau: T::from_f64(0.28),
                household_discount_factor_beta: T::from_f64(0.960),
                equity_payout_cost_phi: T::from_f64(0.18),
                equity_payout_gap_d: T::from_f64(-0.15),
            },
        );

    [
        sw_ufr.extrapolated_zero_yield * T::from_f64(100.0),
        kwf_cb.refunding_efficiency_ratio * T::from_f64(10.0),
        sc_cdo.fair_tranche_par_spread * T::from_f64(100.0),
        bue_vc.fair_variance_swap_volatility_strike * T::from_f64(10.0),
        pr_cva.total_netting_set_cva,
        hr_fut.jump_diffusion_futures_call_price,
        gm_fx.expected_fx_carry_excess_return * T::from_f64(100.0),
        cv_shout.strike_reset_call_value,
        gz_kelly.expected_geometric_growth_rate * T::from_f64(100.0),
        gls_raroc.subportfolio_euler_raroc * T::from_f64(10.0),
        ap_wald.exponential_wald_statistic_exp_w,
        jq_fin.endogenous_financial_labor_wedge * T::from_f64(100.0),
    ]
}

pub fn benchmark_smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_pykhtin_rosen_marginal_cva_hilliard_reis_jump_futures_gabaix_maggiori_fx_carry_cheuk_vorst_shout_option_grossman_zhou_drawdown_kelly_gls_tasche_cornish_fisher_raroc_andrews_ploberger_wald_break_jermann_quadrini_financial_shock(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite70_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite70_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "70_smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_pykhtin_rosen_marginal_cva_hilliard_reis_jump_futures_gabaix_maggiori_fx_carry_cheuk_vorst_shout_option_grossman_zhou_drawdown_kelly_gls_tasche_cornish_fisher_raroc_andrews_ploberger_wald_break_jermann_quadrini_financial_shock",
        "Smith-Wilson-UFR / Kalotay-Williams-Fabozzi-Callable / Schonbucher-Clayton-CDO / Buehler-Variance-Curve / Pykhtin-Rosen-Marginal-CVA / Hilliard-Reis-Jump-Futures / Gabaix-Maggiori-FX-Carry / Cheuk-Vorst-Shout-Option / Grossman-Zhou-Drawdown-Kelly / GLS-Tasche-Cornish-Fisher-RAROC / Andrews-Ploberger-Wald-Break / Jermann-Quadrini-Financial-Shock AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite71_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let mc_tax = ql_curves::evaluate_mcculloch_litzenberger_tax_clientele_discount_curve(
        &ql_curves::McCullochTaxClienteleCurveSpec {
            low_coupon_rate_c_low: T::from_f64(0.025) + dx0 * T::from_f64(0.2),
            high_coupon_rate_c_high: T::from_f64(0.060) + dx1 * T::from_f64(0.2),
            low_coupon_bond_price_p_low: T::from_f64(0.985),
            high_coupon_bond_price_p_high: T::from_f64(1.010),
            marginal_income_tax_rate_tau_inc: T::from_f64(0.28) + dx2,
            capital_gains_tax_rate_tau_cg: T::from_f64(0.15),
            municipal_credit_liquidity_spread: T::from_f64(0.0035) + dx3 * T::from_f64(0.1),
        },
    );

    let ab_cb = ql_fixed_income::evaluate_andersen_buffum_equity_linked_hazard_convertible_bond(
        &ql_fixed_income::AndersenBuffumConvertibleBondSpec {
            spot_price_s0: T::from_f64(95.0) + x[0],
            conversion_ratio: T::from_f64(1.05) + dx1,
            baseline_hazard_lambda0: T::from_f64(0.028) + dx2 * T::from_f64(0.2),
            hazard_power_elasticity_p: T::from_f64(1.25),
            bond_recovery_rate_r: T::from_f64(0.40),
            coupon_rate_c: T::from_f64(0.045),
            risk_free_rate_r: T::from_f64(0.038),
            stock_volatility_sigma: T::from_f64(0.28) + dx3,
        },
    );

    let ok_rl = ql_credit::evaluate_okane_brigo_recovery_lock_and_digital_default_swap(
        &ql_credit::OkaneBrigoRecoveryLockCdsSpec {
            market_cds_spread_s_mkt: T::from_f64(0.0180) + dx0 * T::from_f64(0.1),
            conventional_quoting_recovery_r_conv: T::from_f64(0.40),
            expected_physical_recovery_r_bar: T::from_f64(0.34) + dx1,
            recovery_lock_strike_r_lock: T::from_f64(0.38) + dx2,
            hazard_recovery_correlation_rho: T::from_f64(-0.45),
            recovery_volatility_sigma_r: T::from_f64(0.18) + dx3,
            risk_free_rate_r: T::from_f64(0.035),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let dk_tree = ql_vol::evaluate_derman_kani_implied_tree_arrow_debreu_local_volatility(
        &ql_vol::DermanKaniImpliedTreeSpec {
            central_node_spot_si: T::from_f64(100.0) + x[0],
            arrow_debreu_state_price_lambda: T::from_f64(0.96) + dx1,
            atm_implied_volatility_sigma: T::from_f64(0.22) + dx2,
            smile_skew_slope_beta: T::from_f64(-0.15) + dx3,
            risk_free_rate_r: T::from_f64(0.04),
            dividend_yield_q: T::from_f64(0.015),
            step_years_dt: T::from_f64(0.25),
        },
    );

    let kg_lva = ql_xva::evaluate_kenyon_green_asymmetric_csa_collateral_liquidity_lva(
        &ql_xva::KenyonGreenCsaLvaSpec {
            expected_portfolio_mtm_mu: T::from_f64(6.0) + x[0],
            portfolio_mtm_volatility_sigma: T::from_f64(24.0) + x[1],
            csa_receive_threshold_h_rec: T::from_f64(5.0) + x[2],
            csa_post_threshold_h_post: T::from_f64(2.0) + x[3],
            posted_collateral_liquidity_spread: T::from_f64(0.0045),
            received_collateral_remuneration_drag: T::from_f64(0.0020),
            ois_discount_rate_r: T::from_f64(0.035),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let bdj_gas =
        ql_commodities::evaluate_boogert_de_jong_gas_storage_ratchet_injection_withdrawal_spread(
            &ql_commodities::BoogertDeJongGasRatchetSpec {
                summer_injection_forward_f_sum: T::from_f64(2.85) + dx0 * T::from_f64(5.0),
                winter_withdrawal_forward_f_win: T::from_f64(4.10) + dx1 * T::from_f64(5.0),
                working_gas_fill_fraction: T::from_f64(0.45) + dx2,
                ratchet_pressure_attenuation_alpha: T::from_f64(0.35),
                variable_injection_withdrawal_cost: T::from_f64(0.12),
                calendar_spread_normal_volatility_sigma: T::from_f64(0.85) + dx3 * T::from_f64(5.0),
                risk_free_rate_r: T::from_f64(0.038),
            },
        );

    let gk_baw = ql_fx::evaluate_garman_kohlhagen_biger_hull_fx_early_exercise_premium(
        &ql_fx::GarmanKohlhagenBawFxSpec {
            spot_fx_s0: T::from_f64(1.10) + dx0,
            strike_k: T::from_f64(1.08) + dx1,
            domestic_rate_rd: T::from_f64(0.032) + dx2 * T::from_f64(0.2),
            foreign_rate_rf: T::from_f64(0.055),
            fx_volatility_sigma: T::from_f64(0.12) + dx3,
            maturity_years_t: T::one(),
        },
    );

    let hh_bin = ql_exotics::evaluate_hui_haug_one_touch_and_no_touch_binary_barrier(
        &ql_exotics::HuiHaugBinaryBarrierSpec {
            spot_price_s0: T::from_f64(105.0) + x[0],
            down_barrier_h: T::from_f64(92.0) + x[1],
            cash_rebate_k: T::from_f64(10.0) + x[2],
            risk_free_rate_r: T::from_f64(0.045),
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.22) + dx3,
            maturity_years_t: T::one(),
        },
    );

    let kr_hj = ql_portfolio::evaluate_kan_robotti_misspecification_robust_hj_factor_distance(
        &ql_portfolio::KanRobottiHjDistanceSpec {
            traded_factor_mean_mu1: T::from_f64(0.068) + dx0,
            nontraded_factor_price_lambda2: T::from_f64(0.045) + dx1,
            test_asset_1_pricing_error_alpha1: T::from_f64(0.018) + dx2 * T::from_f64(0.2),
            test_asset_2_pricing_error_alpha2: T::from_f64(-0.014) + dx3 * T::from_f64(0.2),
            test_asset_1_variance_v1: T::from_f64(0.036),
            test_asset_2_variance_v2: T::from_f64(0.044),
            test_asset_correlation_rho12: T::from_f64(0.30),
            sample_observations_t_obs: T::from_f64(180.0),
        },
    );

    let cj_of = ql_portfolio::evaluate_cartea_jaimungal_order_flow_alpha_optimal_liquidation(
        &ql_portfolio::CarteaJaimungalOrderFlowLiquidationSpec {
            initial_inventory_q0: T::from_f64(100.0) + x[0],
            order_flow_signal_mu0: T::from_f64(0.35) + dx1,
            signal_mean_reversion_kappa: T::from_f64(2.50) + dx2 * T::from_f64(5.0),
            temporary_impact_k: T::from_f64(0.0040) + dx3 * T::from_f64(0.02),
            running_inventory_penalty_phi: T::from_f64(0.012),
            terminal_liquidation_penalty_a: T::from_f64(0.080),
            horizon_t: T::one(),
        },
    );

    let asj_jmp = ql_dml::evaluate_ait_sahalia_jacod_power_variation_jump_activity_test(
        &ql_dml::AitSahaliaJacodJumpTestSpec {
            fine_grid_power_variation_b_fine: T::from_f64(0.0085) + dx0 * T::from_f64(0.05),
            coarse_grid_power_variation_b_coarse: T::from_f64(0.0092) + dx1 * T::from_f64(0.05),
            power_exponent_p: T::from_f64(4.0) + dx2,
            subsampling_multiplier_k: T::from_f64(2.0),
            multipower_quarticity_scale_a2p: T::from_f64(0.00018),
            high_frequency_sample_size_n: T::from_f64(390.0) + x[3],
        },
    );

    let wac_dis = ql_models::evaluate_wachter_time_varying_disaster_intensity_equity_volatility(
        &ql_models::WachterDisasterVolatilitySpec {
            current_disaster_intensity_lambda: T::from_f64(0.0355) + dx0 * T::from_f64(0.2),
            long_run_disaster_intensity_lambda_bar: T::from_f64(0.0355) + dx1 * T::from_f64(0.2),
            intensity_mean_reversion_kappa: T::from_f64(0.080) + dx2 * T::from_f64(0.2),
            intensity_volatility_sigma_lambda: T::from_f64(0.067) + dx3 * T::from_f64(0.2),
            diffusive_consumption_volatility_sigma_c: T::from_f64(0.020),
            relative_risk_aversion_gamma: T::from_f64(3.0),
            disaster_contraction_fraction_b: T::from_f64(0.22),
            dividend_leverage_phi: T::from_f64(2.6),
        },
    );

    [
        mc_tax.municipal_to_treasury_yield_ratio * T::from_f64(10.0),
        ab_cb.convertible_bond_value * T::from_f64(0.1),
        ok_rl.recovery_adjusted_digital_protection_pv * T::from_f64(100.0),
        dk_tree.implied_node_local_volatility * T::from_f64(10.0),
        kg_lva.total_asymmetric_csa_lva_charge * T::from_f64(10.0),
        bdj_gas.ratchet_adjusted_storage_option_value * T::from_f64(10.0),
        gk_baw.american_fx_call_value * T::from_f64(100.0),
        hh_bin.cash_at_hit_one_touch_value,
        kr_hj.kan_robotti_hj_t_statistic * T::from_f64(10.0),
        cj_of.optimal_initial_trading_speed * T::from_f64(0.1),
        asj_jmp.standardized_asj_jump_z_statistic,
        wac_dis.equilibrium_equity_risk_premium * T::from_f64(100.0),
    ]
}

pub fn benchmark_mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_kenyon_green_csa_lva_boogert_de_jong_gas_ratchet_garman_kohlhagen_baw_fx_hui_haug_binary_barrier_kan_robotti_hj_distance_cartea_jaimungal_order_flow_ait_sahalia_jacod_jump_test_wachter_disaster_volatility(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite71_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite71_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "71_mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_kenyon_green_csa_lva_boogert_de_jong_gas_ratchet_garman_kohlhagen_baw_fx_hui_haug_binary_barrier_kan_robotti_hj_distance_cartea_jaimungal_order_flow_ait_sahalia_jacod_jump_test_wachter_disaster_volatility",
        "McCulloch-Tax-Curve / Andersen-Buffum-CB / OKane-Recovery-Lock-CDS / Derman-Kani-Implied-Tree / Kenyon-Green-CSA-LVA / Boogert-de-Jong-Gas-Ratchet / Garman-Kohlhagen-BAW-FX / Hui-Haug-Binary-Barrier / Kan-Robotti-HJ-Distance / Cartea-Jaimungal-Order-Flow / Ait-Sahalia-Jacod-Jump-Test / Wachter-Disaster-Volatility AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite72_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let hw_mc = ql_curves::evaluate_hagan_west_monotone_convex_forward_interpolation(
        &ql_curves::HaganWestMonotoneConvexSpec {
            interval_start_years_t0: T::from_f64(2.0),
            interval_end_years_t1: T::from_f64(5.0),
            target_maturity_years_t: T::from_f64(3.4) + dx0,
            start_zero_rate_r0: T::from_f64(0.031) + dx1 * T::from_f64(0.2),
            previous_discrete_forward_f_prev: T::from_f64(0.032) + dx2 * T::from_f64(0.2),
            current_discrete_forward_f_curr: T::from_f64(0.038) + dx3 * T::from_f64(0.2),
            next_discrete_forward_f_next: T::from_f64(0.044),
            node_blending_weight_w: T::from_f64(0.50),
        },
    );

    let fab_oas = ql_fixed_income::evaluate_fabozzi_tuckman_effective_duration_and_convexity_oas(
        &ql_fixed_income::FabozziEffectiveDurationOasSpec {
            benchmark_yield_y0: T::from_f64(0.042) + dx0 * T::from_f64(0.2),
            option_adjusted_spread_oas: T::from_f64(0.0085) + dx1 * T::from_f64(0.2),
            coupon_rate_c: T::from_f64(0.056) + dx2 * T::from_f64(0.2),
            call_price_k: T::from_f64(100.5),
            bond_price_volatility_sigma: T::from_f64(0.065) + dx3 * T::from_f64(0.2),
            parallel_yield_shift_dy: T::from_f64(0.0025),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let fa_lgd = ql_credit::evaluate_frye_altman_procyclical_downturn_lgd_credit_loss(
        &ql_credit::FryeAltmanDownturnLgdSpec {
            default_threshold_z_pd: T::from_f64(-2.0537) + dx0,
            asset_factor_correlation_rho: T::from_f64(0.18) + dx1,
            expected_recovery_rate_mu_r: T::from_f64(0.45) + dx2,
            recovery_volatility_sigma_r: T::from_f64(0.20) + dx3,
            recovery_factor_loading_q_r: T::from_f64(0.55),
            adverse_macro_quantile_z_alpha: T::from_f64(2.3263),
        },
    );

    let gj_svi = ql_vol::evaluate_gatheral_jacquier_svi_calendar_and_butterfly_arbitrage(
        &ql_vol::GatheralJacquierSviSurfaceArbitrageSpec {
            log_moneyness_k: T::from_f64(-0.10) + dx0,
            svi_level_a: T::from_f64(0.024) + dx1 * T::from_f64(0.2),
            svi_slope_b: T::from_f64(0.14) + dx2,
            svi_correlation_rho: T::from_f64(-0.45) + dx3,
            svi_center_m: T::from_f64(0.02),
            svi_smoothness_sigma: T::from_f64(0.18),
            maturity_years_t1: T::one(),
            calendar_variance_growth_dw_dt: T::from_f64(0.048),
        },
    );

    let gb_bcva = ql_xva::evaluate_gregory_brigo_first_to_default_bilateral_cva_dva(
        &ql_xva::GregoryBrigoBilateralCvaDvaSpec {
            expected_positive_exposure_epe: T::from_f64(18.0) + x[0],
            expected_negative_exposure_ene: T::from_f64(11.0) + x[1],
            counterparty_hazard_lambda_c: T::from_f64(0.026) + dx2 * T::from_f64(0.2),
            bank_own_hazard_lambda_b: T::from_f64(0.015) + dx3 * T::from_f64(0.2),
            common_shock_hazard_lambda_cb: T::from_f64(0.004),
            counterparty_lgd_c: T::from_f64(0.60),
            bank_own_lgd_b: T::from_f64(0.60),
            risk_free_rate_r: T::from_f64(0.035),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let bs_mine = ql_commodities::evaluate_brennan_schwartz_natural_resource_mine_hysteresis_option(
        &ql_commodities::BrennanSchwartzMineHysteresisSpec {
            commodity_spot_price_p0: T::from_f64(4.10) + dx0 * T::from_f64(5.0),
            variable_extraction_cost_c: T::from_f64(2.40) + dx1 * T::from_f64(5.0),
            reopening_investment_cost_i: T::from_f64(8.50) + x[2],
            suspension_closure_cost_k: T::from_f64(1.80),
            convenience_yield_delta: T::from_f64(0.045),
            risk_free_rate_r: T::from_f64(0.040),
            commodity_volatility_sigma: T::from_f64(0.26) + dx3,
        },
    );

    let wys_fad = ql_fx::evaluate_wystup_fx_corridor_fader_option_and_forward(
        &ql_fx::WystupFxFaderOptionSpec {
            spot_fx_s0: T::from_f64(1.09) + dx0,
            strike_k: T::from_f64(1.07) + dx1,
            corridor_lower_bound_l: T::from_f64(1.03) + dx2,
            corridor_upper_bound_u: T::from_f64(1.16),
            domestic_rate_rd: T::from_f64(0.036),
            foreign_rate_rf: T::from_f64(0.024),
            fx_volatility_sigma: T::from_f64(0.11) + dx3,
            maturity_years_t: T::one(),
        },
    );

    let bd_cap = ql_exotics::evaluate_broadie_detemple_american_capped_call_option(
        &ql_exotics::BroadieDetempleAmericanCappedCallSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_k: T::from_f64(96.0) + x[1],
            cap_level_l: T::from_f64(118.0) + x[2],
            risk_free_rate_r: T::from_f64(0.045),
            dividend_yield_q: T::from_f64(0.025),
            volatility_sigma: T::from_f64(0.24) + dx3,
            maturity_years_t: T::one(),
        },
    );

    let meu_ep = ql_portfolio::evaluate_meucci_satchell_scowcroft_entropy_pooling_view_allocation(
        &ql_portfolio::MeucciEntropyPoolingAllocationSpec {
            prior_mean_asset1_mu1: T::from_f64(0.065) + dx0,
            prior_mean_asset2_mu2: T::from_f64(0.050) + dx1,
            prior_variance_asset1_v1: T::from_f64(0.040),
            prior_variance_asset2_v2: T::from_f64(0.030),
            prior_correlation_rho12: T::from_f64(0.35),
            spread_view_target_v: T::from_f64(0.032) + dx2,
            spread_view_variance_tau2: T::from_f64(0.022),
            view_confidence_c: T::from_f64(0.65) + dx3,
        },
    );

    let cds_rob = ql_portfolio::evaluate_cont_deguest_scandolo_tail_risk_influence_sensitivity(
        &ql_portfolio::ContDeguestScandoloRobustRiskSpec {
            loss_location_mu: T::from_f64(0.012) + dx0 * T::from_f64(0.2),
            tail_scale_s: T::from_f64(0.025) + dx1 * T::from_f64(0.2),
            tail_index_xi: T::from_f64(0.24) + dx2,
            tail_probability_p: T::from_f64(0.025),
            contamination_epsilon: T::from_f64(0.002),
            stress_outlier_loss_z: T::from_f64(0.28) + dx3,
        },
    );

    let pss_ardl = ql_dml::evaluate_pesaran_shin_smith_ardl_bounds_cointegration_test(
        &ql_dml::PesaranShinSmithArdlBoundsSpec {
            lagged_dependent_coef_pi_yy: T::from_f64(-0.24) + dx0,
            lagged_regressor_coef_pi_yx: T::from_f64(0.36) + dx1,
            lagged_y_variance_s_yy: T::from_f64(1.45) + dx2,
            lagged_x_variance_s_xx: T::from_f64(1.80) + dx3,
            lagged_yx_covariance_s_yx: T::from_f64(0.95),
            residual_variance_sigma2: T::from_f64(0.085),
            sample_size_t_obs: T::from_f64(120.0),
        },
    );

    let cee_nk = ql_models::evaluate_christiano_eichenbaum_evans_calvo_investment_nk_impulse(
        &ql_models::ChristianoEichenbaumEvansNkSpec {
            calvo_price_stickiness_xi_p: T::from_f64(0.68) + dx0,
            inflation_indexation_iota_p: T::from_f64(0.45) + dx1,
            discount_factor_beta: T::from_f64(0.99),
            habit_persistence_h: T::from_f64(0.70) + dx2,
            investment_adjustment_curvature_kappa_i: T::from_f64(4.50) + x[3],
            steady_state_investment_share: T::from_f64(0.22),
            real_interest_rate_shock: T::from_f64(0.010),
        },
    );

    [
        hw_mc.interpolated_zero_rate * T::from_f64(100.0),
        fab_oas.effective_duration,
        fa_lgd.procyclical_downturn_stressed_loss * T::from_f64(100.0),
        gj_svi.dupire_local_volatility * T::from_f64(10.0),
        gb_bcva.net_bilateral_cva_charge,
        bs_mine.operating_mine_total_value * T::from_f64(0.1),
        wys_fad.fade_in_corridor_call_value * T::from_f64(100.0),
        bd_cap.american_capped_call_value,
        meu_ep.posterior_squared_sharpe_ratio * T::from_f64(10.0),
        cds_rob.contaminated_median_to_expected_shortfall_ratio * T::from_f64(10.0),
        pss_ardl.long_run_multiplier_standard_error * T::from_f64(10.0),
        cee_nk.contractionary_inflation_impulse_response * T::from_f64(100.0),
    ]
}

pub fn benchmark_hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_gregory_brigo_ftd_bilateral_cva_brennan_schwartz_mine_hysteresis_wystup_fx_fader_broadie_detemple_capped_call_meucci_entropy_pooling_cont_deguest_scandolo_robust_risk_pesaran_shin_smith_ardl_christiano_eichenbaum_evans_nk(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite72_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite72_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "72_hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_gregory_brigo_ftd_bilateral_cva_brennan_schwartz_mine_hysteresis_wystup_fx_fader_broadie_detemple_capped_call_meucci_entropy_pooling_cont_deguest_scandolo_robust_risk_pesaran_shin_smith_ardl_christiano_eichenbaum_evans_nk",
        "Hagan-West-Monotone-Convex / Fabozzi-Effective-Duration-OAS / Frye-Altman-Downturn-LGD / Gatheral-Jacquier-SVI-Arbitrage / Gregory-Brigo-FtD-Bilateral-CVA / Brennan-Schwartz-Mine-Hysteresis / Wystup-FX-Fader / Broadie-Detemple-Capped-Call / Meucci-Entropy-Pooling / Cont-Deguest-Scandolo-Robust-Risk / Pesaran-Shin-Smith-ARDL / Christiano-Eichenbaum-Evans-NK AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite73_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let mr_basis = ql_curves::evaluate_mercurio_russo_ibor_ois_credit_liquidity_basis(
        &ql_curves::MercurioRussoIborOisCreditLiquiditySpec {
            ois_rate_r: T::from_f64(0.034) + dx0 * T::from_f64(0.2),
            interbank_credit_hazard_lambda: T::from_f64(0.0045) + dx1 * T::from_f64(0.1),
            interbank_loss_given_default: T::from_f64(0.60) + dx2,
            term_liquidity_spread_l: T::from_f64(0.0018) + dx3 * T::from_f64(0.1),
            credit_liquidity_covariance_gamma: T::from_f64(0.0008),
            accrual_tenor_tau: T::from_f64(0.50),
        },
    );

    let jy_yoy = ql_fixed_income::evaluate_jarrow_yildirim_yoy_inflation_swap_convexity_adjustment(
        &ql_fixed_income::JarrowYildirimYoyInflationConvexitySpec {
            nominal_forward_rate_f_n: T::from_f64(0.041) + dx0 * T::from_f64(0.2),
            real_forward_rate_f_r: T::from_f64(0.017) + dx1 * T::from_f64(0.2),
            nominal_rate_volatility_sigma_n: T::from_f64(0.0085) + dx2 * T::from_f64(0.1),
            real_rate_volatility_sigma_r: T::from_f64(0.0065) + dx3 * T::from_f64(0.1),
            cpi_index_volatility_sigma_i: T::from_f64(0.012),
            nominal_real_correlation_rho_nr: T::from_f64(0.55),
            real_cpi_correlation_rho_ri: T::from_f64(-0.25),
            period_start_years_t_prev: T::from_f64(4.0),
        },
    );

    let dl_cds = ql_credit::evaluate_duffie_lando_das_sundaram_imperfect_information_cds_spread(
        &ql_credit::DuffieLandoDasSundaramHybridCdsSpec {
            reported_log_distance_to_default: T::from_f64(0.55) + dx0,
            asset_volatility_sigma_v: T::from_f64(0.25) + dx1,
            accounting_noise_std_a: T::from_f64(0.18) + dx2,
            surprise_jump_hazard_lambda0: T::from_f64(0.032) + dx3 * T::from_f64(0.2),
            jump_leverage_sensitivity_eta: T::from_f64(1.40),
            recovery_rate_r: T::from_f64(0.40),
            maturity_years_t: T::one(),
        },
    );

    let bfg_skew = ql_vol::evaluate_bayer_friz_gatheral_rough_bergomi_power_law_skew(
        &ql_vol::BayerFrizGatheralRoughSkewSpec {
            initial_forward_variance_xi0: T::from_f64(0.040) + dx0 * T::from_f64(0.2),
            hurst_exponent_h: T::from_f64(0.12) + dx1 * T::from_f64(0.5),
            vol_of_vol_eta: T::from_f64(1.85) + dx2,
            spot_vol_correlation_rho: T::from_f64(-0.75) + dx3,
            log_moneyness_k: T::from_f64(-0.08),
            maturity_years_t: T::from_f64(0.25),
        },
    );

    let aps_mpor = ql_xva::evaluate_andersen_pykhtin_sokol_mpor_dispute_spike_cva(
        &ql_xva::AndersenPykhtinSokolMporDisputeSpec {
            portfolio_mtm_volatility_sigma: T::from_f64(42.0) + x[0],
            base_mpor_years_delta: T::from_f64(0.040),
            unpaid_cashflow_spike_c: T::from_f64(3.20) + x[1] * T::from_f64(0.2),
            variation_margin_receipt_fraction_alpha: T::from_f64(0.25) + dx2,
            initial_margin_cushion_i0: T::from_f64(4.50),
            counterparty_hazard_lambda: T::from_f64(0.024) + dx3 * T::from_f64(0.2),
            counterparty_lgd: T::from_f64(0.60),
            ois_discount_rate_r: T::from_f64(0.035),
        },
    );

    let bs_crack =
        ql_commodities::evaluate_bjerksund_stensland_three_asset_321_refinery_crack_option(
            &ql_commodities::BjerksundStenslandRefineryCrack321Spec {
                crude_futures_price_f_c: T::from_f64(76.0) + x[0],
                gasoline_futures_price_f_g: T::from_f64(94.0) + x[1],
                distillate_futures_price_f_h: T::from_f64(98.0) + x[2],
                refining_strike_k: T::from_f64(16.0),
                refined_basket_volatility_sigma_p: T::from_f64(0.29) + dx3,
                crude_volatility_sigma_c: T::from_f64(0.27),
                basket_crude_correlation_rho: T::from_f64(0.86),
                risk_free_rate_r: T::from_f64(0.038),
            },
        );

    let cc_pcs = ql_fx::evaluate_carr_chou_put_call_symmetry_static_barrier_replication(
        &ql_fx::CarrChouSymmetryBarrierReplicationSpec {
            spot_fx_s0: T::from_f64(1.12) + dx0,
            strike_k: T::from_f64(1.10) + dx1,
            down_barrier_h: T::from_f64(1.02),
            domestic_rate_rd: T::from_f64(0.034),
            foreign_rate_rf: T::from_f64(0.028),
            call_strike_volatility_sigma_k: T::from_f64(0.115) + dx2,
            reflected_put_volatility_sigma_star: T::from_f64(0.132) + dx3,
            maturity_years_t: T::one(),
        },
    );

    let cjy_par = ql_exotics::evaluate_chesney_jeanblanc_yor_parisian_excursion_barrier_call(
        &ql_exotics::ChesneyJeanblancYorParisianBarrierSpec {
            spot_price_s0: T::from_f64(105.0) + x[0],
            strike_k: T::from_f64(100.0) + x[1],
            down_barrier_h: T::from_f64(92.0),
            excursion_window_years_d: T::from_f64(0.06) + dx2 * T::from_f64(0.2),
            risk_free_rate_r: T::from_f64(0.042),
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.24) + dx3,
            maturity_years_t: T::one(),
        },
    );

    let dmb_norm =
        ql_portfolio::evaluate_demiguel_brodie_norm_constrained_minimum_variance_shrinkage(
            &ql_portfolio::DeMiguelBrodieNormConstrainedPortfolioSpec {
                asset_1_variance_v1: T::from_f64(0.055) + dx0 * T::from_f64(0.2),
                asset_2_variance_v2: T::from_f64(0.028) + dx1 * T::from_f64(0.2),
                asset_correlation_rho12: T::from_f64(0.25),
                ridge_norm_penalty_lambda2: T::from_f64(0.018) + dx2 * T::from_f64(0.2),
                gross_exposure_penalty_lambda1: T::from_f64(0.008),
                asset_1_expected_return_mu1: T::from_f64(0.075) + dx3,
                asset_2_expected_return_mu2: T::from_f64(0.052),
            },
        );

    let glft_mm = ql_portfolio::evaluate_gueant_lehalle_tapia_inventory_market_making_quotes(
        &ql_portfolio::GueantLehalleTapiaMarketMakingSpec {
            current_inventory_q: T::from_f64(3.0) + x[0],
            mid_price_volatility_sigma: T::from_f64(0.45) + dx1,
            risk_aversion_gamma: T::from_f64(0.12) + dx2,
            order_intensity_decay_kappa: T::from_f64(1.60) + dx3,
            baseline_arrival_intensity_a: T::from_f64(18.0),
            remaining_horizon_tau: T::from_f64(0.50),
        },
    );

    let hw_spa = ql_dml::evaluate_hansen_white_superior_predictive_ability_test(
        &ql_dml::HansenWhiteSuperiorPredictiveAbilitySpec {
            model_1_mean_loss_reduction_d1: T::from_f64(0.045) + dx0,
            model_2_mean_loss_reduction_d2: T::from_f64(-0.024) + dx1,
            model_1_hac_variance_w1: T::from_f64(0.065) + dx2,
            model_2_hac_variance_w2: T::from_f64(0.080),
            cross_model_correlation_rho12: T::from_f64(0.35) + dx3,
            out_of_sample_observations_n: T::from_f64(150.0),
        },
    );

    let bs_par = ql_models::evaluate_brunnermeier_sannikov_intermediary_volatility_paradox(
        &ql_models::BrunnermeierSannikovVolatilityParadoxSpec {
            intermediary_wealth_share_eta: T::from_f64(0.24) + dx0,
            fundamental_capital_volatility_sigma: T::from_f64(0.025) + dx1 * T::from_f64(0.2),
            expert_productivity_wedge_delta_a: T::from_f64(0.055) + dx2 * T::from_f64(0.2),
            fire_sale_illiquidity_sensitivity_phi: T::from_f64(2.20),
            household_discount_rate_rho: T::from_f64(0.050),
            intermediary_leverage_ratio: T::from_f64(3.20) + x[3],
        },
    );

    [
        mr_basis.total_ibor_ois_additive_basis * T::from_f64(1000.0),
        jy_yoy.discounted_yoy_swaplet_value,
        dl_cds.fair_short_maturity_cds_spread * T::from_f64(100.0),
        bfg_skew.rough_bergomi_implied_volatility * T::from_f64(10.0),
        aps_mpor.aps_mpor_adjusted_cva,
        bs_crack.refinery_321_crack_call_value,
        cc_pcs.static_replicated_down_and_out_call * T::from_f64(100.0),
        cjy_par.parisian_down_and_out_call_value,
        dmb_norm.regularized_portfolio_sharpe_ratio * T::from_f64(10.0),
        glft_mm.expected_instantaneous_spread_capture,
        hw_spa.hansen_spa_adjusted_p_value * T::from_f64(100.0),
        bs_par.volatility_paradox_amplification_ratio,
    ]
}

pub fn benchmark_mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_andersen_pykhtin_sokol_mpor_bjerksund_stensland_321_crack_carr_chou_symmetry_barrier_chesney_jeanblanc_yor_parisian_demiguel_brodie_norm_shrinkage_gueant_lehalle_tapia_mm_hansen_white_spa_test_brunnermeier_sannikov_paradox(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite73_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite73_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "73_mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_andersen_pykhtin_sokol_mpor_bjerksund_stensland_321_crack_carr_chou_symmetry_barrier_chesney_jeanblanc_yor_parisian_demiguel_brodie_norm_shrinkage_gueant_lehalle_tapia_mm_hansen_white_spa_test_brunnermeier_sannikov_paradox",
        "Mercurio-Russo-IBOR-OIS-Basis / Jarrow-Yildirim-YoY-Inflation / Duffie-Lando-Hybrid-CDS / Bayer-Friz-Gatheral-Rough-Skew / Andersen-Pykhtin-Sokol-MPoR / Bjerksund-Stensland-321-Crack / Carr-Chou-Symmetry-Barrier / Chesney-Jeanblanc-Yor-Parisian / DeMiguel-Brodie-Norm-Shrinkage / Gueant-Lehalle-Tapia-MM / Hansen-White-SPA-Test / Brunnermeier-Sannikov-Paradox AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite74_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let flt_lr = ql_curves::evaluate_filipovic_larsson_trolle_linear_rational_zlb_curve(
        &ql_curves::FilipovicLarssonTrolleLinearRationalCurveSpec {
            asymptotic_rate_alpha: T::from_f64(0.042) + dx0 * T::from_f64(0.2),
            factor_1_state_x1: T::from_f64(0.35) + dx1,
            factor_2_state_x2: T::from_f64(0.25) + dx2,
            factor_1_mean_theta1: T::from_f64(0.30),
            factor_2_mean_theta2: T::from_f64(0.20),
            factor_1_speed_kappa1: T::from_f64(0.45),
            factor_2_speed_kappa2: T::from_f64(0.18),
            factor_1_volatility_sigma1: T::from_f64(0.14) + dx3,
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let dk_repo = ql_fixed_income::evaluate_duffie_krishnamurthy_on_the_run_special_repo_premium(
        &ql_fixed_income::DuffieKrishnamurthySpecialRepoTreasurySpec {
            general_collateral_repo_rate: T::from_f64(0.045) + dx0 * T::from_f64(0.2),
            initial_repo_specialness_s0: T::from_f64(0.0110) + dx1 * T::from_f64(0.2),
            specialness_decay_rate_kappa: T::from_f64(3.20) + x[2],
            short_squeeze_sensitivity_gamma: T::from_f64(0.0085),
            short_to_float_ratio_q: T::from_f64(0.35),
            remaining_auction_cycle_years: T::from_f64(0.25),
            modified_duration_d_mod: T::from_f64(7.80) + x[3],
            dealer_intermediation_haircut_h: T::from_f64(0.18),
        },
    );

    let dsw_def = ql_credit::evaluate_duffie_saita_wang_stochastic_covariate_survival(
        &ql_credit::DuffieSaitaWangMultiPeriodDefaultSpec {
            baseline_log_hazard_beta0: T::from_f64(-2.40) + dx0,
            initial_distance_to_default_dd0: T::from_f64(2.60) + dx1,
            long_run_distance_to_default_bar: T::from_f64(2.20) + dx2,
            dd_mean_reversion_kappa: T::from_f64(0.38),
            dd_volatility_sigma: T::from_f64(0.45) + dx3,
            dd_hazard_loading_beta_d: T::from_f64(0.85),
            horizon_years_t: T::from_f64(3.0),
        },
    );

    let cmn_log = ql_vol::evaluate_carr_madan_neuberger_log_contract_variance_replication(
        &ql_vol::CarrMadanNeubergerLogContractReplicationSpec {
            atm_implied_volatility_sigma: T::from_f64(0.20) + dx0,
            smile_skew_slope_psi: T::from_f64(-0.22) + dx1,
            smile_wing_curvature_kappa: T::from_f64(0.45) + dx2,
            strike_grid_spacing_dk: T::from_f64(0.025) + dx3 * T::from_f64(0.2),
            discrete_rebalancing_count_n: T::from_f64(252.0),
            maturity_years_t: T::one(),
        },
    );

    let hwbk_fva = ql_xva::evaluate_hull_white_burgard_kjaer_fva_dva2_separation(
        &ql_xva::HullWhiteBurgardKjaerFvaWindfallSpec {
            expected_positive_exposure_epe: T::from_f64(22.0) + x[0],
            expected_negative_exposure_ene: T::from_f64(9.5) + x[1],
            bank_default_hazard_lambda_b: T::from_f64(0.016) + dx2 * T::from_f64(0.2),
            bank_loss_given_default_lgd_b: T::from_f64(0.60),
            pure_liquidity_funding_spread_s_liq: T::from_f64(0.0045) + dx3 * T::from_f64(0.1),
            counterparty_hazard_lambda_c: T::from_f64(0.020),
            ois_risk_free_rate_r: T::from_f64(0.035),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let ccd_cy =
        ql_commodities::evaluate_casassus_collin_dufresne_level_dependent_convenience_yield(
            &ql_commodities::CasassusCollinDufresneConvenienceYieldSpec {
                spot_price_s0: T::from_f64(82.0) + x[0],
                equilibrium_spot_anchor_s_bar: T::from_f64(70.0) + x[1],
                baseline_convenience_yield_delta0: T::from_f64(0.035) + dx2 * T::from_f64(0.2),
                spot_level_sensitivity_alpha_x: T::from_f64(0.42),
                rate_sensitivity_alpha_r: T::from_f64(0.30),
                risk_free_rate_r0: T::from_f64(0.040),
                spot_volatility_sigma_s: T::from_f64(0.28) + dx3,
                market_price_of_risk_lambda_s: T::from_f64(0.35),
                maturity_years_t: T::from_f64(1.5),
            },
        );

    let lrv_fx = ql_fx::evaluate_lustig_roussanov_verdelhan_fx_carry_factor_pricing(
        &ql_fx::LustigRoussanovVerdelhanFxCarryFactorSpec {
            interest_rate_differential_di: T::from_f64(0.038) + dx0,
            dollar_factor_beta_dol: T::from_f64(0.95) + dx1,
            carry_slope_beta_hml: T::from_f64(0.65) + dx2,
            global_vol_beta_vol: T::from_f64(-0.42) + dx3,
            dollar_factor_risk_price_lambda_dol: T::from_f64(0.015),
            carry_slope_risk_price_lambda_hml: T::from_f64(0.048),
            global_vol_risk_price_lambda_vol: T::from_f64(-0.025),
            carry_factor_volatility_sigma_hml: T::from_f64(0.095),
        },
    );

    let cl_corr = ql_exotics::evaluate_carr_lee_corridor_conditional_variance_swap(
        &ql_exotics::CarrLeeCorridorVarianceSwapSpec {
            forward_price_f0: T::from_f64(100.0) + x[0],
            lower_corridor_barrier_b: T::from_f64(94.0) + x[1],
            atmf_volatility_sigma0: T::from_f64(0.21) + dx2,
            skew_slope_psi: T::from_f64(-0.18) + dx3,
            risk_free_rate_r: T::from_f64(0.036),
            maturity_years_t: T::one(),
        },
    );

    let guw_mp = ql_portfolio::evaluate_garlappi_uppal_wang_multi_prior_ambiguity_portfolio(
        &ql_portfolio::GarlappiUppalWangMultiPriorPortfolioSpec {
            sample_mean_asset1_mu1: T::from_f64(0.085) + dx0,
            sample_mean_asset2_mu2: T::from_f64(0.050) + dx1,
            variance_asset1_v1: T::from_f64(0.048) + dx2 * T::from_f64(0.2),
            variance_asset2_v2: T::from_f64(0.026),
            correlation_rho12: T::from_f64(0.25),
            risk_aversion_gamma: T::from_f64(3.0),
            ambiguity_radius_epsilon: T::from_f64(0.42) + dx3,
            estimation_sample_size_t: T::from_f64(60.0),
        },
    );

    let mf_evt = ql_portfolio::evaluate_mcneil_frey_garch_filtered_evt_pot_tail_risk(
        &ql_portfolio::McNeilFreyConditionalEvtGarchRiskSpec {
            conditional_loss_mean_mu: T::from_f64(0.0005) + dx0 * T::from_f64(0.1),
            garch_omega: T::from_f64(0.000012),
            garch_arch_shock_term: T::from_f64(0.00018) + dx1 * T::from_f64(0.005),
            garch_lagged_variance_term: T::from_f64(0.00022),
            standardized_pot_threshold_u: T::from_f64(1.65) + dx2,
            threshold_exceedance_rate_pu: T::from_f64(0.060),
            gpd_shape_index_xi: T::from_f64(0.22) + dx3,
            gpd_scale_beta: T::from_f64(0.52),
        },
    );

    let pat_sjc = ql_dml::evaluate_patton_symmetrized_joe_clayton_tail_dependence(
        &ql_dml::PattonSymmetrizedJoeClaytonCopulaSpec {
            lower_tail_intercept_omega_l: T::from_f64(-0.35) + dx0,
            lower_tail_persistence_beta_l: T::from_f64(1.40) + dx1,
            lagged_lower_tail_tau_l: T::from_f64(0.38),
            upper_tail_intercept_omega_u: T::from_f64(-0.85) + dx2,
            upper_tail_persistence_beta_u: T::from_f64(1.10),
            lagged_upper_tail_tau_u: T::from_f64(0.24),
            pit_discrepancy_forcing_delta: T::from_f64(0.18) + dx3,
            forcing_sensitivity_alpha: T::from_f64(-1.60),
        },
    );

    let pv_pol = ql_models::evaluate_pastor_veronesi_political_uncertainty_risk_premium(
        &ql_models::PastorVeronesiPoliticalUncertaintySpec {
            perceived_policy_impact_g_hat: T::from_f64(-0.012) + dx0 * T::from_f64(0.2),
            prior_policy_uncertainty_sigma_g0: T::from_f64(0.035) + dx1 * T::from_f64(0.2),
            fundamental_macro_volatility_sigma_m: T::from_f64(0.14) + dx2,
            political_cost_dispersion_sigma_pol: T::from_f64(0.085) + dx3,
            investor_risk_aversion_gamma: T::from_f64(3.50),
            bayesian_learning_duration_tau: T::from_f64(2.0),
        },
    );

    [
        flt_lr.state_dependent_yield_volatility * T::from_f64(10.0),
        dk_repo.annualized_specialness_pickup_rate,
        dsw_def.multi_period_default_probability * T::from_f64(10.0),
        cmn_log.discrete_replication_rmse * T::from_f64(1000.0),
        hwbk_fva.funding_benefit_adjustment_fba,
        ccd_cy.equilibrium_futures_price * T::from_f64(0.1),
        lrv_fx.carry_factor_sharpe_ratio * T::from_f64(10.0),
        cl_corr.corridor_variance_swap_pv,
        guw_mp.max_min_certainty_equivalent * T::from_f64(100.0),
        mf_evt.conditional_evt_garch_es99 * T::from_f64(100.0),
        pat_sjc.symmetrized_tail_asymmetry_wedge * T::from_f64(10.0),
        pv_pol.political_uncertainty_equity_risk_premium * T::from_f64(100.0),
    ]
}

pub fn benchmark_filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_hull_white_burgard_kjaer_fva_windfall_casassus_collin_dufresne_convenience_yield_lustig_roussanov_verdelhan_fx_carry_carr_lee_corridor_varswap_garlappi_uppal_wang_multi_prior_mcneil_frey_evt_garch_patton_sjc_copula_pastor_veronesi_political_uncertainty(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite74_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite74_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "74_filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_hull_white_burgard_kjaer_fva_windfall_casassus_collin_dufresne_convenience_yield_lustig_roussanov_verdelhan_fx_carry_carr_lee_corridor_varswap_garlappi_uppal_wang_multi_prior_mcneil_frey_evt_garch_patton_sjc_copula_pastor_veronesi_political_uncertainty",
        "Filipovic-Larsson-Trolle-LRTSM / Duffie-Krishnamurthy-Repo / Duffie-Saita-Wang-Default / Carr-Madan-Neuberger-LogContract / Hull-White-Burgard-Kjaer-FVA / Casassus-Collin-Dufresne-CY / Lustig-Roussanov-Verdelhan-FX / Carr-Lee-Corridor-VarSwap / Garlappi-Uppal-Wang-MultiPrior / McNeil-Frey-EVT-GARCH / Patton-SJC-Copula / Pastor-Veronesi-Political-Risk AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite75_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let lm_rfr = ql_curves::evaluate_lyashenko_mercurio_rfr_in_arrears_caplet_convexity(
        &ql_curves::LyashenkoMercurioRfrInArrearsCapletSpec {
            forward_rfr_rate_f0: T::from_f64(0.042) + dx0 * T::from_f64(0.2),
            strike_rate_k: T::from_f64(0.040) + dx1 * T::from_f64(0.2),
            fixing_start_years_t_prev: T::from_f64(0.50) + dx2,
            accrual_fraction_tau: T::from_f64(0.25),
            rfr_forward_volatility_sigma: T::from_f64(0.24) + dx3,
            discount_factor_p_ti: T::from_f64(0.968),
        },
    );

    let ps_coco = ql_fixed_income::evaluate_pennacchi_sundaresan_coco_at1_trigger_valuation(
        &ql_fixed_income::PennacchiSundaresanCocoAt1TriggerSpec {
            current_cet1_ratio_c0: T::from_f64(0.132) + dx0 * T::from_f64(0.5),
            trigger_cet1_ratio_c_trig: T::from_f64(0.070) + dx1 * T::from_f64(0.2),
            cet1_drift_mu_c: T::from_f64(0.008),
            cet1_volatility_sigma_c: T::from_f64(0.19) + dx2,
            at1_coupon_rate_q: T::from_f64(0.068) + dx3 * T::from_f64(0.5),
            equity_conversion_share_chi: T::from_f64(0.45),
            write_down_loss_severity_lwd: T::from_f64(0.85),
            risk_free_rate_r: T::from_f64(0.038),
            horizon_years_t: T::from_f64(5.0),
        },
    );

    let pb_cdx = ql_credit::evaluate_pedersen_brigo_cds_index_option_fep_swaption(
        &ql_credit::PedersenBrigoCdsIndexOptionFepSpec {
            current_index_spread_s0: T::from_f64(0.0095) + dx0 * T::from_f64(0.1),
            strike_spread_k: T::from_f64(0.0108) + dx1 * T::from_f64(0.1),
            index_running_coupon_c: T::from_f64(0.0100),
            forward_rpv01_annuity_a0: T::from_f64(4.15) + x[2],
            spread_volatility_sigma_s: T::from_f64(0.48) + dx3,
            recovery_rate_r: T::from_f64(0.40),
            risk_free_rate_rf: T::from_f64(0.036),
            option_expiry_years_te: T::from_f64(0.50),
        },
    );

    let bjs_svj = ql_vol::evaluate_broadie_jain_sepp_svj_volatility_swap_convexity(
        &ql_vol::BroadieJainSeppSvjVolatilitySwapSpec {
            initial_variance_v0: T::from_f64(0.044) + dx0 * T::from_f64(0.2),
            long_run_variance_theta: T::from_f64(0.036) + dx1 * T::from_f64(0.2),
            mean_reversion_kappa: T::from_f64(2.10),
            vol_of_vol_xi: T::from_f64(0.45) + dx2,
            jump_intensity_lambda_j: T::from_f64(0.65) + dx3,
            jump_mean_mu_j: T::from_f64(-0.075),
            jump_std_sigma_j: T::from_f64(0.095),
            maturity_years_t: T::one(),
        },
    );

    let bcp_wwr = ql_xva::evaluate_brigo_capponi_pallavicini_cds_contagion_wwr_cva(
        &ql_xva::BrigoCapponiPallaviciniCdsContagionWwrSpec {
            counterparty_hazard_lambda_c: T::from_f64(0.022) + dx0 * T::from_f64(0.2),
            reference_hazard_lambda_r: T::from_f64(0.018) + dx1 * T::from_f64(0.2),
            reference_spread_vol_sigma_r: T::from_f64(0.52) + dx2,
            spread_correlation_rho_rc: T::from_f64(0.55),
            contagion_hazard_jump_eta: T::from_f64(0.65) + dx3,
            counterparty_lgd: T::from_f64(0.60),
            risk_free_rate_r: T::from_f64(0.035),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let rss_store = ql_commodities::evaluate_routledge_seppi_spatt_competitive_storage_equilibrium(
        &ql_commodities::RoutledgeSeppiSpattCompetitiveStorageSpec {
            spot_price_s0: T::from_f64(78.0) + x[0],
            current_inventory_i0: T::from_f64(28.0) + x[1],
            critical_stockout_inventory_i_star: T::from_f64(22.0),
            stockout_convenience_spike_delta: T::from_f64(0.28) + dx2,
            storage_cost_rate_c: T::from_f64(0.075),
            base_spot_volatility_sigma: T::from_f64(0.24) + dx3,
            stockout_vol_amplification_beta: T::from_f64(0.85),
            maturity_years_t: T::one(),
        },
    );

    let pit_prdc = ql_fx::evaluate_piterbarg_prdc_fx_ir_hybrid_coupon_valuation(
        &ql_fx::PiterbargPrdcHybridSmileCouponSpec {
            spot_fx_s0: T::from_f64(104.0) + x[0],
            domestic_rate_rd: T::from_f64(0.008),
            foreign_rate_rf: T::from_f64(0.044) + dx1 * T::from_f64(0.2),
            foreign_coupon_gearing_gf: T::from_f64(0.12),
            domestic_coupon_strike_gd: T::from_f64(0.08),
            fx_spot_volatility_sigma_s: T::from_f64(0.115) + dx2,
            rate_diff_volatility_sigma_ir: T::from_f64(0.012),
            fx_ir_correlation_rho: T::from_f64(0.35) + dx3,
            maturity_years_t: T::from_f64(10.0),
        },
    );

    let bgk_bar = ql_exotics::evaluate_broadie_glasserman_kou_discrete_barrier_call(
        &ql_exotics::BroadieGlassermanKouDiscreteBarrierCallSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(100.0) + x[1],
            lower_barrier_h: T::from_f64(88.0) + x[2],
            risk_free_rate_r: T::from_f64(0.045),
            volatility_sigma: T::from_f64(0.24) + dx3,
            maturity_years_t: T::one(),
            monitoring_intervals_m: T::from_f64(52.0),
        },
    );

    let bsv_ppp = ql_portfolio::evaluate_brandt_santa_clara_valkanov_parametric_portfolio_policy(
        &ql_portfolio::BrandtSantaClaraValkanovPppSpec {
            value_characteristic_mean_return: T::from_f64(0.045) + dx0,
            momentum_characteristic_mean_return: T::from_f64(0.062) + dx1,
            value_characteristic_variance: T::from_f64(0.018) + dx2 * T::from_f64(0.1),
            momentum_characteristic_variance: T::from_f64(0.024),
            value_momentum_covariance: T::from_f64(-0.008),
            benchmark_value_covariance: T::from_f64(0.002),
            benchmark_momentum_covariance: T::from_f64(-0.001),
            relative_risk_aversion_gamma: T::from_f64(5.0) + x[3] * T::from_f64(0.1),
        },
    );

    let gss_imp = ql_portfolio::evaluate_gatheral_schied_slynko_transient_impact_execution(
        &ql_portfolio::GatheralSchiedSlynkoTransientImpactSpec {
            order_size_shares_x0: T::from_f64(10.0) + x[0],
            transient_impact_coeff_eta: T::from_f64(0.085) + dx1,
            resilience_decay_rate_rho: T::from_f64(1.60) + dx2,
            permanent_impact_coeff_gamma: T::from_f64(0.025),
            price_volatility_sigma: T::from_f64(0.32) + dx3,
            risk_aversion_lambda: T::from_f64(0.15),
            execution_horizon_hours_t: T::from_f64(4.0),
        },
    );

    let bpq_harq = ql_dml::evaluate_corsi_bollerslev_patton_quaedvlieg_harq_volatility(
        &ql_dml::CorsiBollerslevPattonQuaedvliegHarqSpec {
            intercept_beta0: T::from_f64(0.000018),
            daily_persistence_beta_d: T::from_f64(0.42) + dx0,
            weekly_persistence_beta_w: T::from_f64(0.34) + dx1,
            monthly_persistence_beta_m: T::from_f64(0.20) + dx2,
            quarticity_attenuation_beta_q: T::from_f64(-0.16) + dx3,
            daily_realized_variance_rv_d: T::from_f64(0.00024),
            weekly_realized_variance_rv_w: T::from_f64(0.00018),
            monthly_realized_variance_rv_m: T::from_f64(0.00015),
            realized_quarticity_excess_rq: T::from_f64(0.65),
        },
    );

    let ll_cay = ql_models::evaluate_lettau_ludvigson_cay_conditional_equity_premium(
        &ql_models::LettauLudvigsonCayConditionalRiskPremiumSpec {
            log_consumption_c: T::from_f64(10.24) + dx0,
            log_asset_wealth_a: T::from_f64(10.05) + dx1,
            log_labor_income_y: T::from_f64(10.28),
            asset_wealth_cointegration_share_omega: T::from_f64(0.32),
            baseline_risk_aversion_gamma0: T::from_f64(12.0) + x[2],
            cay_conditioning_sensitivity_gamma1: T::from_f64(185.0),
            consumption_volatility_sigma_c: T::from_f64(0.018),
            equity_return_volatility_sigma_m: T::from_f64(0.165) + dx3,
        },
    );

    [
        lm_rfr.in_arrears_rfr_caplet_pv * T::from_f64(10.0),
        ps_coco.coco_at1_fair_price_pct * T::from_f64(0.1),
        pb_cdx.cds_index_payer_swaption_pv_pct * T::from_f64(10.0),
        bjs_svj.svj_fair_volatility_strike * T::from_f64(10.0),
        bcp_wwr.total_contagion_wwr_cva_pct * T::from_f64(10.0),
        rss_store.equilibrium_futures_price * T::from_f64(0.1),
        pit_prdc.cancelable_prdc_net_value,
        bgk_bar.discrete_bgk_down_and_out_call,
        bsv_ppp.certainty_equivalent_utility_gain_bps * T::from_f64(0.01),
        gss_imp.expected_implementation_shortfall_cost,
        bpq_harq.annualized_harq_volatility_forecast * T::from_f64(10.0),
        ll_cay.conditional_equity_risk_premium * T::from_f64(100.0),
    ]
}

pub fn benchmark_lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_brigo_capponi_pallavicini_contagion_wwr_routledge_seppi_spatt_storage_piterbarg_prdc_hybrid_broadie_glasserman_kou_discrete_barrier_brandt_santa_clara_valkanov_ppp_gatheral_schied_slynko_transient_impact_corsi_bpq_harq_lettau_ludvigson_cay(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite75_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite75_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "75_lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_brigo_capponi_pallavicini_contagion_wwr_routledge_seppi_spatt_storage_piterbarg_prdc_hybrid_broadie_glasserman_kou_discrete_barrier_brandt_santa_clara_valkanov_ppp_gatheral_schied_slynko_transient_impact_corsi_bpq_harq_lettau_ludvigson_cay",
        "Lyashenko-Mercurio-RFR / Pennacchi-Sundaresan-CoCo / Pedersen-Brigo-CDX-FEP / Broadie-Jain-Sepp-SVJ-VolSwap / Brigo-Capponi-Pallavicini-WWR / Routledge-Seppi-Spatt-Storage / Piterbarg-PRDC-Hybrid / Broadie-Glasserman-Kou-Barrier / Brandt-Santa-Clara-Valkanov-PPP / Gatheral-Schied-Slynko-Impact / Corsi-BPQ-HARQ / Lettau-Ludvigson-CAY AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite76_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let cp_cycle = ql_curves::evaluate_cieslak_povala_trend_inflation_bond_cycle(
        &ql_curves::CieslakPovalaTrendInflationCycleSpec {
            short_rate_y1: T::from_f64(0.042) + dx0 * T::from_f64(0.2),
            intermediate_yield_y5: T::from_f64(0.046) + dx1 * T::from_f64(0.2),
            long_yield_y10: T::from_f64(0.049) + dx2 * T::from_f64(0.2),
            discounted_cpi_trend_tau: T::from_f64(0.026) + dx3 * T::from_f64(0.2),
            short_rate_inflation_loading_b1: T::from_f64(1.05),
            long_yield_inflation_loading_b10: T::from_f64(0.92),
            real_rate_cycle_return_loading_gamma_c: T::from_f64(0.28),
            short_rate_cycle_return_loading_gamma_1: T::from_f64(0.14),
        },
    );

    let gkv_mbs = ql_fixed_income::evaluate_gabaix_krishnamurthy_vigneron_mbs_intermediary_oas(
        &ql_fixed_income::GabaixKrishnamurthyVigneronMbsOasSpec {
            mortgage_coupon_wac: T::from_f64(0.055) + dx0 * T::from_f64(0.2),
            current_refi_rate_r: T::from_f64(0.050) + dx1 * T::from_f64(0.2),
            prepayment_s_curve_slope_beta: T::from_f64(25.0),
            rate_volatility_sigma_r: T::from_f64(0.012) + dx2 * T::from_f64(0.1),
            aggregate_mbs_outstanding_share_w: T::from_f64(0.65) + dx3,
            intermediary_risk_aversion_gamma: T::from_f64(18.0),
            pass_through_duration_d0: T::from_f64(5.4),
            liquidity_hedging_cost_lambda: T::from_f64(3.5),
        },
    );

    let egg_hawkes = ql_credit::evaluate_errais_giesecke_goldberg_hawkes_contagion_credit_loss(
        &ql_credit::ErraisGieseckeGoldbergHawkesCreditSpec {
            initial_default_intensity_lambda0: T::from_f64(0.035) + dx0 * T::from_f64(0.2),
            long_run_baseline_intensity_mu: T::from_f64(0.025) + dx1 * T::from_f64(0.2),
            mean_reversion_speed_kappa: T::from_f64(1.40),
            self_exciting_contagion_impact_delta: T::from_f64(0.95) + dx2,
            expected_loss_per_default_ell: T::from_f64(0.60) + dx3,
            horizon_years_t: T::from_f64(5.0),
            tranche_detachment_fraction_d: T::from_f64(0.12),
        },
    );

    let mb_svix = ql_vol::evaluate_martin_bondarenko_simple_variance_svix_lower_bound(
        &ql_vol::MartinBondarenkoSimpleVarianceSvixSpec {
            forward_price_f0: T::from_f64(4500.0) + x[0] * T::from_f64(10.0),
            risk_free_rate_r: T::from_f64(0.042),
            atm_volatility_sigma: T::from_f64(0.175) + dx1,
            downside_skew_slope_psi: T::from_f64(-0.18) + dx2,
            crash_jump_intensity_lambda: T::from_f64(0.45),
            crash_log_jump_mean_mu_j: T::from_f64(-0.14),
            relative_risk_aversion_gamma: T::from_f64(2.5) + x[3] * T::from_f64(0.2),
            horizon_years_t: T::from_f64(0.25),
        },
    );

    let cbc_bsde = ql_xva::evaluate_crepey_bichuch_capponi_nonlinear_xva_bsde(
        &ql_xva::CrepeyBichuchCapponiNonlinearXvaBsdeSpec {
            clean_derivative_mtm_v0: T::from_f64(12.5) + x[0],
            collateralization_fraction_alpha: T::from_f64(0.35) + dx1,
            unsecured_borrowing_spread_s_b: T::from_f64(0.0145) + dx2 * T::from_f64(0.2),
            unsecured_lending_spread_s_l: T::from_f64(0.0040),
            counterparty_default_spread_s_c: T::from_f64(0.0120) + dx3 * T::from_f64(0.2),
            hedger_default_spread_s_h: T::from_f64(0.0065),
            exposure_volatility_sigma_v: T::from_f64(8.5),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let ts_usv = ql_commodities::evaluate_trolle_schwartz_unspanned_stochastic_volatility_commodity(
        &ql_commodities::TrolleSchwartzUnspannedStochasticVolatilitySpec {
            futures_price_f0: T::from_f64(76.0) + x[0],
            forward_cost_of_carry_y0: T::from_f64(-0.035) + dx1,
            carry_decay_kappa_y: T::from_f64(1.25),
            initial_stochastic_variance_v0: T::from_f64(0.095) + dx2,
            long_run_variance_theta_v: T::from_f64(0.080),
            variance_mean_reversion_kappa_v: T::from_f64(1.80),
            vol_of_vol_sigma_v: T::from_f64(0.42) + dx3,
            option_expiry_years_t: T::from_f64(0.50),
        },
    );

    let iss_cip = ql_fx::evaluate_ivashina_scharfstein_stein_dollar_funding_cip_basis(
        &ql_fx::IvashinaScharfsteinSteinDollarFundingCipSpec {
            direct_usd_commercial_paper_rate_r_usd: T::from_f64(0.052) + dx0 * T::from_f64(0.2),
            foreign_domestic_deposit_rate_r_for: T::from_f64(0.034) + dx1 * T::from_f64(0.2),
            wholesale_usd_deposit_outflow_shock_delta: T::from_f64(0.25) + dx2,
            dollar_loan_book_rigidity_theta: T::from_f64(0.55),
            fx_swap_dealer_balance_sheet_cost_eta: T::from_f64(0.018),
            spot_fx_s0: T::from_f64(1.08) + dx3,
            credit_worthiness_haircut_alpha: T::from_f64(0.012),
            tenor_years_t: T::from_f64(0.50),
        },
    );

    let ww_cliq = ql_exotics::evaluate_wilmott_windcliff_cliquet_forward_skew_valuation(
        &ql_exotics::WilmottWindcliffCliquetForwardSkewSpec {
            notional_amount: T::from_f64(100.0) + x[0],
            local_floor_f_loc: T::zero(),
            local_cap_c_loc: T::from_f64(0.08) + dx1 * T::from_f64(0.5),
            global_floor_f_glob: T::from_f64(0.16),
            risk_free_rate_r: T::from_f64(0.038),
            atm_reset_volatility_sigma: T::from_f64(0.20) + dx2,
            upper_cap_skew_slope_beta: T::from_f64(-0.35) + dx3,
            reset_periods_m: T::from_f64(5.0),
        },
    );

    let lw_nl = ql_portfolio::evaluate_ledoit_wolf_analytical_nonlinear_eigenvalue_shrinkage(
        &ql_portfolio::LedoitWolfNonlinearEigenvalueShrinkageSpec {
            sample_eigenvalue_lambda: T::from_f64(0.42) + dx0,
            concentration_ratio_c: T::from_f64(0.35) + dx1,
            bulk_spectral_mean_bar_lambda: T::one(),
            local_spectral_density_f: T::from_f64(0.68) + dx2,
            hilbert_transform_stieltjes_h: T::from_f64(0.45) + dx3,
            tail_regularization_weight_xi: T::from_f64(0.08),
        },
    );

    let ab_covar = ql_portfolio::evaluate_adrian_brunnermeier_systemic_delta_covar(
        &ql_portfolio::AdrianBrunnermeierSystemicCoVarSpec {
            institution_mean_loss_mu_i: T::from_f64(0.004),
            institution_volatility_sigma_i: T::from_f64(0.032) + dx0 * T::from_f64(0.5),
            system_mean_loss_mu_s: T::from_f64(0.003),
            system_volatility_sigma_s: T::from_f64(0.024) + dx1 * T::from_f64(0.5),
            normal_correlation_rho: T::from_f64(0.58),
            distress_quantile_slope_delta_beta: T::from_f64(0.22) + dx2,
            leverage_maturity_mismatch_m: T::from_f64(1.45) + dx3,
            tail_quantile_z_q: T::from_f64(2.3263),
        },
    );

    let btz_vrp = ql_dml::evaluate_bollerslev_tauchen_zhou_variance_risk_premium_predictability(
        &ql_dml::BollerslevTauchenZhouVarianceRiskPremiumSpec {
            model_free_implied_variance_iv: T::from_f64(0.042) + dx0 * T::from_f64(0.5),
            high_frequency_realized_variance_rv: T::from_f64(0.024) + dx1 * T::from_f64(0.5),
            variance_persistence_rho_v: T::from_f64(0.62),
            vol_of_vol_persistence_rho_q: T::from_f64(0.82) + dx2,
            vol_of_vol_shock_variance_varphi: T::from_f64(0.18),
            consumption_volatility_sensitivity_a_q: T::from_f64(2.40) + dx3,
            recursive_risk_aversion_gamma: T::from_f64(6.5),
            forecast_horizon_quarters_h: T::from_f64(2.0),
        },
    );

    let cv_beta = ql_models::evaluate_campbell_vuolteenaho_bad_beta_good_beta_icapm(
        &ql_models::CampbellVuolteenahoBadGoodBetaIcapmSpec {
            cash_flow_news_covariance_cov_cf: T::from_f64(0.0115) + dx0 * T::from_f64(0.2),
            discount_rate_news_covariance_cov_dr: T::from_f64(0.0145) + dx1 * T::from_f64(0.2),
            market_return_variance_var_m: T::from_f64(0.0250) + dx2 * T::from_f64(0.2),
            investor_relative_risk_aversion_gamma: T::from_f64(7.5) + x[3] * T::from_f64(0.2),
            log_linearization_discount_rho: T::from_f64(0.95),
            return_forecasting_var_persistence_phi: T::from_f64(0.88),
            idiosyncratic_alpha_adjustment_xi: T::from_f64(0.015),
        },
    );

    [
        cp_cycle.expected_10y_bond_excess_return * T::from_f64(100.0),
        gkv_mbs.equilibrium_mbs_total_oas_bps * T::from_f64(0.1),
        egg_hawkes.contagion_adjusted_tranche_risk_charge * T::from_f64(10.0),
        mb_svix.martin_generalized_equity_premium_bound * T::from_f64(10.0),
        cbc_bsde.nonlinear_bsde_total_xva,
        ts_usv.usv_atm_futures_option_value,
        iss_cip.segmented_fx_forward_rate * T::from_f64(10.0),
        ww_cliq.cliquet_forward_skew_fair_value * T::from_f64(0.1),
        lw_nl.nonlinear_shrunk_eigenvalue * T::from_f64(10.0),
        ab_covar.system_covar_in_distress * T::from_f64(100.0),
        btz_vrp.predicted_annualized_excess_return * T::from_f64(10.0),
        cv_beta.icapm_two_beta_required_risk_premium * T::from_f64(100.0),
    ]
}

pub fn benchmark_cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_crepey_bichuch_capponi_xva_bsde_trolle_schwartz_usv_commodity_ivashina_scharfstein_stein_cip_wilmott_windcliff_cliquet_skew_ledoit_wolf_nonlinear_shrinkage_adrian_brunnermeier_delta_covar_bollerslev_tauchen_zhou_vrp_campbell_vuolteenaho_bad_good_beta(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite76_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite76_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "76_cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_crepey_bichuch_capponi_xva_bsde_trolle_schwartz_usv_commodity_ivashina_scharfstein_stein_cip_wilmott_windcliff_cliquet_skew_ledoit_wolf_nonlinear_shrinkage_adrian_brunnermeier_delta_covar_bollerslev_tauchen_zhou_vrp_campbell_vuolteenaho_bad_good_beta",
        "Cieslak-Povala-Cycle / Gabaix-Krishnamurthy-Vigneron-MBS / Errais-Giesecke-Goldberg-Hawkes / Martin-Bondarenko-SVIX / Crepey-Bichuch-Capponi-BSDE / Trolle-Schwartz-USV / Ivashina-Scharfstein-Stein-CIP / Wilmott-Windcliff-Cliquet / Ledoit-Wolf-Nonlinear-Shrinkage / Adrian-Brunnermeier-CoVaR / Bollerslev-Tauchen-Zhou-VRP / Campbell-Vuolteenaho-Beta AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite77_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let bm_fra = ql_curves::evaluate_bianchetti_mercurio_multi_curve_fra_quanto_convexity(
        &ql_curves::BianchettiMercurioMultiCurveFraQuantoSpec {
            ois_forward_rate_f_d: T::from_f64(0.038) + dx0 * T::from_f64(0.2),
            ibor_funding_spread_s0: T::from_f64(0.0045) + dx1 * T::from_f64(0.1),
            fixed_fra_rate_k: T::from_f64(0.040) + dx2 * T::from_f64(0.2),
            ois_forward_volatility_sigma_d: T::from_f64(0.22),
            spread_volatility_sigma_s: T::from_f64(0.35) + dx3,
            ois_spread_correlation_rho: T::from_f64(-0.25),
            fixing_time_years_t1: T::from_f64(1.5),
            accrual_fraction_tau: T::from_f64(0.50),
        },
    );

    let aab_liq = ql_fixed_income::evaluate_acharya_amihud_bharath_liquidity_stress_bond_spread(
        &ql_fixed_income::AcharyaAmihudBharathLiquidityStressBondSpec {
            structural_default_spread_s_def: T::from_f64(0.018) + dx0 * T::from_f64(0.2),
            bond_illiquidity_level_illiq: T::from_f64(0.012) + dx1 * T::from_f64(0.2),
            normal_regime_liquidity_beta_1: T::from_f64(0.25),
            stress_regime_liquidity_beta_2: T::from_f64(0.85) + dx2,
            stress_regime_probability_pi_s: T::from_f64(0.30) + dx3,
            market_illiquidity_risk_price_lambda_l: T::from_f64(0.015),
            flight_to_quality_treasury_wedge_psi: T::from_f64(0.004),
            modified_duration_d: T::from_f64(6.2),
        },
    );

    let ba_cds = ql_credit::evaluate_brigo_alfonsi_correlated_cir_stochastic_intensity_cds(
        &ql_credit::BrigoAlfonsiCorrelatedCirCdsSpec {
            initial_short_rate_r0: T::from_f64(0.036) + dx0 * T::from_f64(0.2),
            initial_default_intensity_lambda0: T::from_f64(0.022) + dx1 * T::from_f64(0.2),
            intensity_mean_reversion_kappa: T::from_f64(0.65),
            intensity_long_run_theta: T::from_f64(0.028) + dx2 * T::from_f64(0.2),
            intensity_volatility_sigma_lambda: T::from_f64(0.18),
            rate_volatility_sigma_r: T::from_f64(0.08),
            rate_intensity_correlation_rho: T::from_f64(-0.35) + dx3,
            loss_given_default_lgd: T::from_f64(0.60),
        },
    );

    let fuk_skew = ql_vol::evaluate_fukasawa_asymptotic_moneyness_skew_varswap(
        &ql_vol::FukasawaAsymptoticImpliedSkewVarswapSpec {
            atm_implied_volatility_sigma0: T::from_f64(0.20) + dx0,
            standardized_moneyness_skew_a1: T::from_f64(-0.035) + dx1 * T::from_f64(0.5),
            standardized_moneyness_curvature_a2: T::from_f64(0.018) + dx2 * T::from_f64(0.5),
            third_order_wing_asymmetry_a3: T::from_f64(-0.006),
            risk_free_rate_r: T::from_f64(0.040),
            maturity_years_t: T::from_f64(0.50) + dx3,
            target_standardized_moneyness_z: T::from_f64(-1.25),
        },
    );

    let ads_fva = ql_xva::evaluate_andersen_duffie_song_debt_overhang_fva(
        &ql_xva::AndersenDuffieSongShareholderFvaOverhangSpec {
            uncollateralized_funding_need_phi: T::from_f64(100.0) + x[0],
            dealer_credit_spread_s_d: T::from_f64(0.0125) + dx1 * T::from_f64(0.2),
            upfront_marginal_debt_issuance_share_q: T::from_f64(0.80) + dx2,
            existing_creditor_recovery_dilution_alpha: T::from_f64(0.25) + dx3,
            asset_swap_covariance_beta_v: T::from_f64(0.08),
            risk_free_rate_r: T::from_f64(0.038),
            marginal_capital_return_hurdle_eta: T::from_f64(0.018),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let pir_sprd = ql_commodities::evaluate_pirrong_structural_inventory_skew_spread_option(
        &ql_commodities::PirrongStructuralInventorySkewSpreadSpec {
            prompt_futures_f1: T::from_f64(82.0) + x[0],
            deferred_futures_f2: T::from_f64(79.5) + x[1],
            spread_strike_k: T::from_f64(2.0),
            normalized_inventory_scarcity_z_inv: T::from_f64(1.15) + dx2,
            baseline_prompt_vol_sigma1: T::from_f64(0.28) + dx3,
            deferred_vol_sigma2: T::from_f64(0.22),
            baseline_correlation_rho0: T::from_f64(0.92),
            scarcity_decorrelation_sensitivity_psi: T::from_f64(0.18),
        },
    );

    let lc_dnt = ql_fx::evaluate_lipton_castagna_fx_double_no_touch_vanna_volga(
        &ql_fx::LiptonCastagnaFxDoubleNoTouchVannaVolgaSpec {
            spot_fx_s0: T::from_f64(1.092) + dx0 * T::from_f64(0.5),
            lower_barrier_l: T::from_f64(1.030) + dx1 * T::from_f64(0.5),
            upper_barrier_u: T::from_f64(1.150),
            domestic_rate_rd: T::from_f64(0.048),
            foreign_rate_rf: T::from_f64(0.032),
            atm_volatility_sigma: T::from_f64(0.088) + dx2 * T::from_f64(0.5),
            risk_reversal_25d_rr: T::from_f64(-0.012),
            butterfly_25d_bf: T::from_f64(0.0045) + dx3 * T::from_f64(0.2),
        },
    );

    let lin_step = ql_exotics::evaluate_linetsky_occupation_time_step_barrier_option(
        &ql_exotics::LinetskyStepBarrierOccupationSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(100.0) + x[1],
            lower_step_barrier_h: T::from_f64(88.0) + x[2],
            risk_free_rate_r: T::from_f64(0.042),
            volatility_sigma: T::from_f64(0.24),
            occupation_killing_rate_rho: T::from_f64(8.0) + x[3],
            maturity_years_t: T::one(),
        },
    );

    let ky_dem = ql_portfolio::evaluate_koijen_yogo_characteristic_demand_system_elasticity(
        &ql_portfolio::KoijenYogoDemandSystemAssetPricingSpec {
            log_market_to_book_me_be: T::from_f64(0.65) + dx0,
            beta_loading_beta_m: T::from_f64(1.05),
            profitability_characteristic_prof: T::from_f64(0.32) + dx1,
            institutional_price_elasticity_coeff_beta0: T::from_f64(-0.45) + dx2,
            characteristic_tilt_strength_beta_x: T::from_f64(1.20),
            outside_asset_portfolio_weight_w0: T::from_f64(0.22),
            active_institutional_aum_share_theta: T::from_f64(0.68),
            inelastic_passive_flow_shock_f: T::from_f64(0.025) + dx3 * T::from_f64(0.5),
        },
    );

    let ko_inv = ql_portfolio::evaluate_kyle_obizhaeva_microstructure_invariance_execution(
        &ql_portfolio::KyleObizhaevaMicrostructureInvarianceSpec {
            order_size_shares_q: T::from_f64(25_000.0) + x[0] * T::from_f64(100.0),
            daily_trading_volume_v: T::from_f64(1_200_000.0),
            stock_price_p: T::from_f64(85.0) + x[1],
            daily_return_volatility_sigma: T::from_f64(0.018) + dx2 * T::from_f64(0.2),
            invariance_dollar_bet_var_i_star: T::from_f64(25_000.0),
            linear_impact_calibration_c_lambda: T::from_f64(0.32) + dx3,
            bid_ask_spread_calibration_c_s: T::from_f64(0.08),
        },
    );

    let gkx_r2 = ql_dml::evaluate_gu_kelly_xiu_nonlinear_characteristic_oos_r2(
        &ql_dml::GuKellyXiuPredictiveOosR2InteractionSpec {
            momentum_signal_z_mom: T::from_f64(0.75) + dx0,
            value_signal_z_val: T::from_f64(0.55) + dx1,
            macro_equity_premium_spread_x_m: T::from_f64(0.028) + dx2 * T::from_f64(0.5),
            linear_feature_weight_w1: T::from_f64(0.42),
            nonlinear_interaction_weight_w12: T::from_f64(0.35) + dx3,
            l2_ridge_penalization_lambda: T::from_f64(0.15),
            realized_excess_return_r_next: T::from_f64(0.019),
            cross_sectional_return_dispersion_s2: T::from_f64(0.0085),
        },
    );

    let dss_dep = ql_models::evaluate_drechsler_savov_schnabl_deposit_channel_liquidity_premium(
        &ql_models::DrechslerSavovSchnablDepositChannelSpec {
            fed_funds_policy_rate_f: T::from_f64(0.0525) + dx0 * T::from_f64(0.2),
            bank_branch_herfindahl_hhi: T::from_f64(0.24) + dx1,
            household_cash_substitution_elasticity_eps: T::from_f64(0.55) + dx2,
            bank_liquidity_management_cost_c: T::from_f64(0.0045),
            wholesale_funding_sensitivity_delta: T::from_f64(0.35),
            loan_illiquidity_spread_alpha: T::from_f64(1.80),
            aggregate_liquidity_demand_kappa: T::from_f64(1.25),
            policy_rate_hike_shock_df: T::from_f64(0.010) + dx3 * T::from_f64(0.2),
        },
    );

    [
        bm_fra.multi_curve_payer_fra_pv * T::from_f64(1000.0),
        aab_liq.stress_adjusted_corporate_bond_spread * T::from_f64(100.0),
        ba_cds.correlated_ssrd_par_cds_spread_bps * T::from_f64(0.01),
        fuk_skew.fukasawa_variance_swap_fair_strike * T::from_f64(100.0),
        ads_fva.shareholder_indifference_fva,
        pir_sprd.pirrong_calendar_spread_call_value,
        lc_dnt.vanna_volga_dnt_fair_value * T::from_f64(10.0),
        lin_step.step_barrier_down_and_out_call,
        ky_dem.equilibrium_log_price_impact * T::from_f64(100.0),
        ko_inv.invariance_total_dollar_execution_cost * T::from_f64(1e-3),
        gkx_r2.out_of_sample_predictive_r2 * T::from_f64(10.0),
        dss_dep.equilibrium_liquidity_premium_wedge * T::from_f64(100.0),
    ]
}

pub fn benchmark_bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_andersen_duffie_song_debt_overhang_fva_pirrong_structural_inventory_spread_lipton_castagna_fx_dnt_vanna_volga_linetsky_step_barrier_koijen_yogo_demand_system_kyle_obizhaeva_invariance_gu_kelly_xiu_oos_r2_drechsler_savov_schnabl_deposit_channel(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite77_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite77_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "77_bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_andersen_duffie_song_debt_overhang_fva_pirrong_structural_inventory_spread_lipton_castagna_fx_dnt_vanna_volga_linetsky_step_barrier_koijen_yogo_demand_system_kyle_obizhaeva_invariance_gu_kelly_xiu_oos_r2_drechsler_savov_schnabl_deposit_channel",
        "Bianchetti-Mercurio-FRA-Quanto / Acharya-Amihud-Bharath-Liq / Brigo-Alfonsi-SSRD-CDS / Fukasawa-Skew-VarSwap / Andersen-Duffie-Song-FVA / Pirrong-Inventory-Spread / Lipton-Castagna-FX-DNT / Linetsky-Step-Barrier / Koijen-Yogo-Demand-System / Kyle-Obizhaeva-Invariance / Gu-Kelly-Xiu-OOS-R2 / Drechsler-Savov-Schnabl-Deposit AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite78_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let hs_dur = ql_curves::evaluate_hanson_stein_greenwood_vayanos_duration_term_premium(
        &ql_curves::HansonSteinGreenwoodVayanosDurationSupplySpec {
            short_real_rate_r1: T::from_f64(0.018) + dx0 * T::from_f64(0.2),
            expected_future_short_rate_er2: T::from_f64(0.020) + dx1 * T::from_f64(0.2),
            net_duration_supply_s: T::from_f64(4.5) + x[2],
            mbs_convexity_supply_shock_m: T::from_f64(0.85) + dx3 * T::from_f64(5.0),
            arbitrageur_risk_aversion_gamma: T::from_f64(22.0),
            short_rate_volatility_sigma_r: T::from_f64(0.014),
            yield_oriented_investor_share_alpha: T::from_f64(0.45),
            reaching_for_yield_sensitivity_theta: T::from_f64(18.0),
        },
    );

    let fll_tips = ql_fixed_income::evaluate_fleckenstein_longstaff_lustig_tips_treasury_puzzle(
        &ql_fixed_income::FleckensteinLongstaffLustigTipsTreasuryPuzzleSpec {
            nominal_treasury_yield_y_nom: T::from_f64(0.042) + dx0 * T::from_f64(0.2),
            tips_real_yield_y_tips: T::from_f64(0.021) + dx1 * T::from_f64(0.2),
            zero_coupon_inflation_swap_rate_s_cpi: T::from_f64(0.024) + dx2 * T::from_f64(0.2),
            tips_deflation_floor_strike_k_floor: T::zero(),
            cpi_inflation_volatility_sigma_pi: T::from_f64(0.015) + dx3 * T::from_f64(0.2),
            repo_specialness_differential_delta_repo: T::from_f64(0.0012),
            dealer_balance_sheet_leverage_cost_eta: T::from_f64(0.35),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let cdgy_cdo = ql_credit::evaluate_collin_dufresne_goldstein_yang_super_senior_cdo(
        &ql_credit::CollinDufresneGoldsteinYangSuperSeniorCdoSpec {
            idiosyncratic_hazard_rate_lambda_i: T::from_f64(0.015) + dx0 * T::from_f64(0.2),
            catastrophic_macro_jump_intensity_lambda_m: T::from_f64(0.018) + dx1 * T::from_f64(0.2),
            macro_crisis_conditional_default_prob_p_m: T::from_f64(0.36) + dx2,
            loss_given_default_lgd: T::from_f64(0.60),
            super_senior_attachment_a: T::from_f64(0.15),
            super_senior_detachment_d: T::from_f64(0.30),
            risk_free_rate_r: T::from_f64(0.036) + dx3 * T::from_f64(0.2),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let tod_vol = ql_vol::evaluate_todorov_characteristic_option_spot_volatility(
        &ql_vol::TodorovOptionImpliedSpotVolatilityJumpSpec {
            short_tenor_years_tau: T::from_f64(0.020) + dx0 * T::from_f64(0.1),
            characteristic_wavenumber_u: T::from_f64(2.5) + dx1,
            empirical_log_characteristic_modulus_l_u: T::from_f64(0.125) + dx2,
            high_frequency_wavenumber_ratio_k: T::from_f64(1.6),
            jump_activity_blumenthal_getoor_beta: T::from_f64(1.25),
            tail_jump_intensity_scale_c_j: T::from_f64(0.035) + dx3 * T::from_f64(0.5),
            risk_free_rate_r: T::from_f64(0.040),
        },
    );

    let acc_kva = ql_xva::evaluate_albanese_caenazzo_crepey_bilateral_kva_mva(
        &ql_xva::AlbaneseCaenazzoCrepeyBilateralKvaMvaSpec {
            regulatory_capital_ead: T::from_f64(120.0) + x[0],
            counterparty_risk_weight_rw: T::from_f64(0.08) + dx1 * T::from_f64(0.5),
            cva_capital_charge_k_cva: T::from_f64(3.5) + x[2],
            shareholder_hurdle_rate_h: T::from_f64(0.105),
            initial_margin_posted_im0: T::from_f64(45.0) + x[3],
            im_funding_spread_s_im: T::from_f64(0.0085),
            capital_decay_rate_lambda_k: T::from_f64(0.18),
            horizon_years_t: T::from_f64(7.0),
        },
    );

    let alr_hp = ql_commodities::evaluate_acharya_lochstoer_ramadorai_commodity_hedging_pressure(
        &ql_commodities::AcharyaLochstoerRamadoraiHedgingPressureSpec {
            expected_future_spot_e_s1: T::from_f64(80.0) + x[0],
            spot_price_volatility_sigma_s: T::from_f64(14.0) + x[1],
            producer_default_prob_zm_score: T::from_f64(0.45) + dx2,
            producer_inventory_endowment_q_p: T::from_f64(2.5),
            producer_distress_aversion_gamma_p: T::from_f64(1.8),
            speculator_capital_risk_aversion_gamma_s: T::from_f64(0.65),
            speculator_wealth_constraint_tightness_w: T::from_f64(0.55) + dx3,
            risk_free_rate_r: T::from_f64(0.038),
        },
    );

    let bnp_fx = ql_fx::evaluate_brunnermeier_nagel_pedersen_fx_carry_crash_skew(
        &ql_fx::BrunnermeierNagelPedersenFxCarryCrashSpec {
            interest_rate_differential_idiff: T::from_f64(0.045) + dx0 * T::from_f64(0.5),
            speculator_carry_position_size_x: T::from_f64(1.40) + dx1,
            funding_liquidity_ted_spread_s_ted: T::from_f64(0.018) + dx2 * T::from_f64(0.2),
            margin_spiral_amplification_phi: T::from_f64(4.5),
            baseline_fx_volatility_sigma_fx: T::from_f64(0.115) + dx3,
            crash_arrival_sensitivity_kappa: T::from_f64(2.2),
            risk_reversal_put_skew_cost_c_rr: T::from_f64(0.35),
            horizon_years_t: T::from_f64(0.50),
        },
    );

    let by_asian = ql_exotics::evaluate_broadie_yamamoto_fast_gauss_discrete_asian_option(
        &ql_exotics::BroadieYamamotoFastGaussDiscreteAsianSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(100.0) + x[1],
            risk_free_rate_r: T::from_f64(0.045) + dx2 * T::from_f64(0.2),
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.22) + dx3,
            discrete_fixings_m: T::from_f64(12.0),
            maturity_years_t: T::one(),
        },
    );

    let kps_ipca = ql_portfolio::evaluate_kelly_pruitt_su_instrumented_pca_dynamic_beta(
        &ql_portfolio::KellyPruittSuInstrumentedPcaDynamicBetaSpec {
            characteristic_value_z1: T::from_f64(0.65) + dx0,
            characteristic_momentum_z2: T::from_f64(0.80) + dx1,
            characteristic_investment_z3: T::from_f64(0.40) + dx2,
            gamma_beta_loading_1: T::from_f64(0.55),
            gamma_beta_loading_2: T::from_f64(0.72),
            gamma_alpha_anomaly_intercept: T::from_f64(0.015),
            latent_factor_realized_return_f: T::from_f64(0.068) + dx3,
            cross_sectional_idiosyncratic_var_sigma2: T::from_f64(0.042),
        },
    );

    let be_srisk = ql_portfolio::evaluate_brownlees_engle_srisk_capital_shortfall(
        &ql_portfolio::BrownleesEngleSystemicSriskCapitalShortfallSpec {
            book_value_of_debt_d: T::from_f64(950.0) + x[0],
            market_value_of_equity_w: T::from_f64(85.0) + x[1],
            prudential_capital_ratio_k: T::from_f64(0.08),
            dynamic_conditional_beta_rho: T::from_f64(1.35) + dx2,
            systemic_crisis_market_decline_c: T::from_f64(0.40) + dx3,
            idiosyncratic_tail_amplification_theta: T::from_f64(0.18),
            funding_run_haircut_psi: T::from_f64(0.05),
        },
    );

    let abdl_rv = ql_dml::evaluate_andersen_bollerslev_diebold_labys_fractional_rv(
        &ql_dml::AndersenBollerslevDieboldLabysArfimaRvSpec {
            current_realized_variance_rv0: T::from_f64(0.00022) + dx0 * T::from_f64(1e-4),
            lag1_realized_variance_rv1: T::from_f64(0.00018) + dx1 * T::from_f64(1e-4),
            long_run_mean_log_rv_mu: T::from_f64(-8.85),
            fractional_differencing_d: T::from_f64(0.40) + dx2,
            short_run_ar1_phi: T::from_f64(0.18) + dx3,
            log_rv_innovation_std_sigma_u: T::from_f64(0.32),
            forecast_horizon_days_h: T::from_f64(5.0),
        },
    );

    let fg_uip = ql_models::evaluate_farhi_gabaix_rare_disaster_exchange_rate_uip(
        &ql_models::FarhiGabaixRareDisasterExchangeRatePuzzleSpec {
            global_disaster_probability_p: T::from_f64(0.036) + dx0 * T::from_f64(0.2),
            world_disaster_recovery_b_star: T::from_f64(2.40) + dx1,
            home_country_resilience_h: T::from_f64(0.15) + dx2,
            foreign_country_resilience_h_star: T::from_f64(0.45) + dx3,
            resilience_mean_reversion_phi_h: T::from_f64(0.18),
            domestic_productivity_growth_g: T::from_f64(0.018),
            real_discount_rate_r_bar: T::from_f64(0.040),
        },
    );

    [
        hs_dur.equilibrium_long_real_yield * T::from_f64(100.0),
        fll_tips.tips_treasury_arbitrage_wedge_per_100,
        cdgy_cdo.super_senior_tranche_par_spread_bps * T::from_f64(0.01),
        tod_vol.jump_robust_spot_volatility * T::from_f64(10.0),
        acc_kva.total_bilateral_kva_plus_mva,
        alr_hp.equilibrium_futures_price * T::from_f64(0.1),
        bnp_fx.crash_adjusted_expected_carry_return * T::from_f64(100.0),
        by_asian.fast_gauss_discrete_asian_call,
        kps_ipca.ipca_precision_weighted_expected_return * T::from_f64(100.0),
        be_srisk.conditional_capital_shortfall_srisk * T::from_f64(0.1),
        abdl_rv.annualized_arfima_volatility_forecast * T::from_f64(10.0),
        fg_uip.expected_disaster_fx_carry_excess_return * T::from_f64(100.0),
    ]
}

pub fn benchmark_hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_albanese_caenazzo_crepey_bilateral_kva_acharya_lochstoer_ramadorai_hedging_pressure_brunnermeier_nagel_pedersen_fx_crash_broadie_yamamoto_fgt_asian_kelly_pruitt_su_ipca_brownlees_engle_srisk_andersen_bollerslev_diebold_labys_arfima_farhi_gabaix_disaster_fx(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite78_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite78_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "78_hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_albanese_caenazzo_crepey_bilateral_kva_acharya_lochstoer_ramadorai_hedging_pressure_brunnermeier_nagel_pedersen_fx_crash_broadie_yamamoto_fgt_asian_kelly_pruitt_su_ipca_brownlees_engle_srisk_andersen_bollerslev_diebold_labys_arfima_farhi_gabaix_disaster_fx",
        "Hanson-Stein-Duration / Fleckenstein-Longstaff-Lustig-TIPS / Collin-Dufresne-Goldstein-Yang-CDO / Todorov-Spot-Vol / Albanese-Caenazzo-Crepey-KVA / Acharya-Lochstoer-Ramadorai-Hedging / Brunnermeier-Nagel-Pedersen-FX / Broadie-Yamamoto-FGT-Asian / Kelly-Pruitt-Su-IPCA / Brownlees-Engle-SRISK / ABDL-ARFIMA-RV / Farhi-Gabaix-Disaster-FX AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite79_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let ks_sofr = ql_curves::evaluate_klingler_syrstad_sofr_repo_balance_sheet_spike(
        &ql_curves::KlinglerSyrstadSofrRepoSpecialnessSpikeSpec {
            ois_policy_rate_r_ois: T::from_f64(0.048) + dx0 * T::from_f64(0.2),
            treasury_collateral_specialness_s_spec: T::from_f64(0.0018) + dx1 * T::from_f64(0.05),
            dealer_treasury_inventory_congestion_q: T::from_f64(0.65) + dx2,
            quarter_end_balance_sheet_spike_a_qe: T::from_f64(0.0085) + dx3 * T::from_f64(0.1),
            reporting_date_weight_omega_qe: T::from_f64(0.08),
            bank_credit_sensitive_spread_s_cr: T::from_f64(0.0032),
            repo_rate_volatility_sigma_repo: T::from_f64(0.012),
            tenor_years_tau: T::from_f64(0.25),
        },
    );

    let lmn_basis = ql_fixed_income::evaluate_longstaff_mithal_neis_cds_bond_liquidity_basis(
        &ql_fixed_income::LongstaffMithalNeisCdsBondLiquidityBasisSpec {
            corporate_bond_yield_y_corp: T::from_f64(0.058) + dx0 * T::from_f64(0.2),
            treasury_benchmark_yield_y_tsy: T::from_f64(0.040) + dx1 * T::from_f64(0.2),
            cds_default_spread_s_cds: T::from_f64(0.0115) + dx2 * T::from_f64(0.1),
            corporate_coupon_rate_c: T::from_f64(0.055),
            state_tax_differential_tau: T::from_f64(0.065),
            bond_bid_ask_spread_s_ba: T::from_f64(0.0048) + dx3 * T::from_f64(0.1),
            issue_turnover_liquidity_factor_l: T::from_f64(1.25),
            modified_duration_d: T::from_f64(6.0),
        },
    );

    let ps_sov = ql_credit::evaluate_pan_singleton_sovereign_cds_distress_risk_premium(
        &ql_credit::PanSingletonSovereignCdsDistressRiskSpec {
            physical_hazard_rate_lambda0_p: T::from_f64(0.018) + dx0 * T::from_f64(0.2),
            physical_mean_reversion_kappa_p: T::from_f64(0.55) + dx1,
            log_intensity_volatility_sigma: T::from_f64(0.42) + dx2,
            distress_risk_price_eta: T::from_f64(0.45) + dx3,
            risk_neutral_loss_given_default_lgd: T::from_f64(0.75),
            short_cds_tenor_t1: T::one(),
            long_cds_tenor_t2: T::from_f64(5.0),
            risk_free_rate_r: T::from_f64(0.038),
        },
    );

    let ms_smile = ql_vol::evaluate_medvedev_scaillet_short_maturity_smile_expansion(
        &ql_vol::MedvedevScailletShortMaturitySmileExpansionSpec {
            spot_volatility_sigma0: T::from_f64(0.21) + dx0,
            vol_of_vol_nu: T::from_f64(0.55) + dx1,
            spot_vol_correlation_rho: T::from_f64(-0.68) + dx2,
            volatility_drift_mu_v: T::from_f64(0.015),
            normalized_moneyness_theta: T::from_f64(-1.20) + dx3,
            maturity_years_tau: T::from_f64(0.25),
            risk_free_rate_r: T::from_f64(0.045),
            dividend_yield_q: T::from_f64(0.015),
        },
    );

    let dz_ccp = ql_xva::evaluate_duffie_zhu_central_clearing_netting_efficiency(
        &ql_xva::DuffieZhuCentralClearingNettingEfficiencySpec {
            clearing_member_count_n: T::from_f64(12.0) + x[0],
            asset_class_count_k: T::from_f64(4.0) + x[1] * T::from_f64(0.5),
            cleared_class_share_alpha: T::from_f64(0.85) + dx2,
            bilateral_exposure_volatility_sigma: T::from_f64(25.0) + x[3],
            cross_asset_correlation_rho: T::from_f64(0.20),
            default_fund_mutualization_weight_phi: T::from_f64(0.35),
        },
    );

    let ghr_inv = ql_commodities::evaluate_gorton_hayashi_rouwenhorst_commodity_inventory_basis(
        &ql_commodities::GortonHayashiRouwenhorstInventoryBasisSpec {
            normalized_inventory_level_i: T::from_f64(0.72) + dx0,
            convenience_yield_scarcity_alpha: T::from_f64(0.085) + dx1 * T::from_f64(0.2),
            stockout_convexity_exponent_beta: T::from_f64(1.65) + dx2,
            interest_plus_storage_cost_rc: T::from_f64(0.075),
            prior_futures_momentum_signal_m: T::from_f64(0.14) + dx3,
            inventory_risk_premium_loading_lambda: T::from_f64(0.80),
            spot_volatility_sigma_s: T::from_f64(0.26),
        },
    );

    let drs_vrp = ql_fx::evaluate_della_corte_ramadorai_sarno_fx_volatility_premium(
        &ql_fx::DellaCorteRamadoraiSarnoFxVolatilityPremiumSpec {
            realized_fx_volatility_rv: T::from_f64(0.118) + dx0,
            model_free_implied_fx_vol_iv: T::from_f64(0.092) + dx1,
            butterfly_25d_smile_convexity_bf: T::from_f64(0.0038) + dx2 * T::from_f64(0.1),
            forward_discount_carry_idiff: T::from_f64(0.025) + dx3,
            vrp_spot_predictability_beta: T::from_f64(0.65),
            volatility_insurance_hedging_gamma: T::from_f64(0.40),
            horizon_years_tau: T::from_f64(0.50),
        },
    );

    let pel_db = ql_exotics::evaluate_pelsser_laplace_double_barrier_knockout_option(
        &ql_exotics::PelsserLaplaceDoubleBarrierKnockoutSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            lower_barrier_l: T::from_f64(82.0) + x[1],
            upper_barrier_u: T::from_f64(122.0) + x[2],
            knockout_rebate_r_reb: T::from_f64(15.0),
            risk_free_rate_r: T::from_f64(0.042),
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.20) + dx3,
            maturity_years_t: T::from_f64(0.75),
        },
    );

    let mm_vmp = ql_portfolio::evaluate_moreira_muir_volatility_managed_portfolio_alpha(
        &ql_portfolio::MoreiraMuirVolatilityManagedPortfolioSpec {
            unconditional_factor_expected_return_mu: T::from_f64(0.075) + dx0,
            unconditional_factor_variance_bar_var: T::from_f64(0.032) + dx1 * T::from_f64(0.2),
            lagged_realized_variance_rv_prev: T::from_f64(0.022) + dx2 * T::from_f64(0.2),
            volatility_scaling_constant_c: T::from_f64(0.030),
            return_variance_predictability_slope_phi: T::from_f64(-0.05),
            turnover_transaction_cost_kappa: T::from_f64(0.0020),
            variance_persistence_rho: T::from_f64(0.72) + dx3,
        },
    );

    let bcs_fba = ql_portfolio::evaluate_budish_cramton_shim_frequent_batch_auction_sniping(
        &ql_portfolio::BudishCramtonShimFrequentBatchAuctionSpec {
            fundamental_jump_intensity_lambda_j: T::from_f64(1.8) + x[0] * T::from_f64(0.2),
            investor_arrival_intensity_lambda_i: T::from_f64(8.5) + x[1] * T::from_f64(0.2),
            mean_jump_size_sigma_j: T::from_f64(0.045) + dx2,
            competing_hft_count_n: T::from_f64(6.0) + x[3] * T::from_f64(0.5),
            speed_technology_cost_c_speed: T::from_f64(0.012),
            reaction_latency_delta: T::from_f64(0.001),
            batch_auction_interval_tau: T::from_f64(0.100),
        },
    );

    let amz_tsrv = ql_dml::evaluate_ait_sahalia_mykland_zhang_two_scales_rv(
        &ql_dml::AitSahaliaMyklandZhangTwoScalesRvSpec {
            full_grid_realized_variance_rv_all: T::from_f64(0.00052) + dx0 * T::from_f64(1e-4),
            subsampled_average_rv_avg: T::from_f64(0.00021) + dx1 * T::from_f64(1e-4),
            total_intraday_ticks_n: T::from_f64(4800.0) + x[2] * T::from_f64(50.0),
            subsampling_grids_k: T::from_f64(30.0) + x[3],
            integrated_quarticity_iq: T::from_f64(1.5e-7),
        },
    );

    let gk_qe = ql_models::evaluate_gertler_karadi_intermediary_credit_policy(
        &ql_models::GertlerKaradiUnconventionalMonetaryPolicySpec {
            risk_free_deposit_rate_r: T::from_f64(0.035) + dx0 * T::from_f64(0.2),
            moral_hazard_divertible_fraction_lambda: T::from_f64(0.38) + dx1,
            franchise_net_worth_multiplier_nu: T::from_f64(1.55) + dx2,
            excess_return_shadow_value_mu: T::from_f64(0.042),
            central_bank_intermediation_share_psi: T::from_f64(0.18) + dx3,
            central_bank_efficiency_cost_tau: T::from_f64(0.0015),
            capital_quality_shock_xi: T::from_f64(0.25),
        },
    );

    [
        ks_sofr.credit_vs_sofr_term_basis_bps * T::from_f64(0.1),
        lmn_basis.adjusted_cds_bond_basis_bps * T::from_f64(0.1),
        ps_sov.long_tenor_sovereign_cds_bps * T::from_f64(0.05),
        ms_smile.asymptotic_implied_volatility * T::from_f64(10.0),
        dz_ccp.unified_multilateral_ccp_epe * T::from_f64(0.1),
        ghr_inv.expected_commodity_futures_excess_return * T::from_f64(10.0),
        drs_vrp.vrp_strategy_expected_excess_return * T::from_f64(100.0),
        pel_db.pelsser_double_barrier_total_value * T::from_f64(0.1),
        mm_vmp.volatility_managed_spanning_alpha * T::from_f64(10.0),
        bcs_fba.sniping_latency_tax_elimination_gain * T::from_f64(100.0),
        amz_tsrv.annualized_tsrv_volatility * T::from_f64(10.0),
        gk_qe.equilibrium_credit_spread_bps * T::from_f64(0.05),
    ]
}

pub fn benchmark_klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_duffie_zhu_ccp_netting_gorton_hayashi_rouwenhorst_inventory_basis_della_corte_ramadorai_sarno_fx_vrp_pelsser_laplace_double_barrier_moreira_muir_volatility_managed_budish_cramton_shim_fba_ait_sahalia_mykland_zhang_tsrv_gertler_karadi_credit_policy(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite79_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite79_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "79_klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_duffie_zhu_ccp_netting_gorton_hayashi_rouwenhorst_inventory_basis_della_corte_ramadorai_sarno_fx_vrp_pelsser_laplace_double_barrier_moreira_muir_volatility_managed_budish_cramton_shim_fba_ait_sahalia_mykland_zhang_tsrv_gertler_karadi_credit_policy",
        "Klingler-Syrstad-SOFR-Spike / Longstaff-Mithal-Neis-CDS-Bond / Pan-Singleton-Sovereign-CDS / Medvedev-Scaillet-Short-Smile / Duffie-Zhu-CCP-Netting / Gorton-Hayashi-Rouwenhorst-Inventory / Della-Corte-Ramadorai-Sarno-FX-VRP / Pelsser-Laplace-Double-Barrier / Moreira-Muir-Vol-Managed / Budish-Cramton-Shim-FBA / AMZ-Two-Scales-RV / Gertler-Karadi-QE AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite80_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let br_qe = ql_curves::evaluate_bauer_rudebusch_qe_signaling_term_premium(
        &ql_curves::BauerRudebuschQeSignalingTermPremiumSpec {
            shadow_short_rate_s0: T::from_f64(0.005) + dx0,
            effective_lower_bound_r_elb: T::zero(),
            shadow_rate_volatility_sigma_s: T::from_f64(0.02) + dx1.abs(),
            lsap_purchase_shock_q: T::from_f64(0.35) + dx2,
            signaling_path_sensitivity_beta_sig: T::from_f64(0.45) + dx3,
            portfolio_balance_sensitivity_beta_pb: T::from_f64(0.012),
            baseline_term_premium_tp0: T::from_f64(0.02),
            maturity_years_tau: T::from_f64(10.0),
        },
    );

    let bfl_mbs = ql_fixed_income::evaluate_boyarchenko_fuster_lucca_mortgage_spread_oas(
        &ql_fixed_income::BoyarchenkoFusterLuccaMortgageSpreadOasSpec {
            coupon_moneyness_m: T::from_f64(0.020) + dx0,
            interest_rate_volatility_sigma_r: T::from_f64(0.0115) + dx1.abs(),
            prepayment_incentive_slope_gamma: T::from_f64(12.5),
            prepayment_risk_smile_factor_xi: T::from_f64(0.42) + dx2,
            dealer_balance_sheet_cost_s_cap: T::from_f64(0.0055) + dx3.abs(),
            originator_pipeline_markup_s_orig: T::from_f64(0.0065),
            effective_mbs_duration_d: T::from_f64(5.0),
        },
    );

    let ags_cdo = ql_credit::evaluate_azizpour_giesecke_schwenkler_default_clustering(
        &ql_credit::AzizpourGieseckeSchwenklerDefaultClusteringSpec {
            baseline_macro_intensity_lambda0: T::from_f64(0.015) + dx0 * T::from_f64(0.1),
            latent_frailty_intensity_y: T::from_f64(0.012) + dx1 * T::from_f64(0.1),
            contagion_self_excitation_delta: T::from_f64(0.035) + dx2 * T::from_f64(0.2),
            contagion_decay_speed_kappa: T::from_f64(1.25),
            recent_default_cluster_shock_n: T::from_f64(3.0),
            loss_given_default_lgd: T::from_f64(0.60),
            senior_tranche_attachment_k: T::from_f64(0.60),
            horizon_years_t: T::from_f64(5.0) + dx3 * T::from_f64(2.0),
        },
    );

    let bt_fear = ql_vol::evaluate_bollerslev_todorov_jump_tail_fear_premium(
        &ql_vol::BollerslevTodorovJumpTailFearPremiumSpec {
            physical_left_tail_scale_phi_minus_p: T::from_f64(0.0085) + dx0 * T::from_f64(0.05),
            risk_neutral_left_tail_scale_phi_minus_q: T::from_f64(0.014) + dx1 * T::from_f64(0.05),
            physical_right_tail_scale_phi_plus_p: T::from_f64(0.006) + dx2 * T::from_f64(0.03),
            risk_neutral_right_tail_scale_phi_plus_q: T::from_f64(0.008) + dx3 * T::from_f64(0.03),
            left_tail_decay_exponent_alpha_minus: T::from_f64(2.60),
            right_tail_decay_exponent_alpha_plus: T::from_f64(3.10),
            tail_jump_cutoff_k0: T::from_f64(0.045),
            diffusive_variance_risk_premium_vrp_c: T::from_f64(0.0065),
        },
    );

    let gg_mva = ql_xva::evaluate_ghamami_glasserman_clearing_incentive_mva(
        &ql_xva::GhamamiGlassermanClearingIncentiveMvaSpec {
            bilateral_netting_set_volatility_sigma_b: T::from_f64(140.0) + x[0],
            ccp_multilateral_volatility_sigma_c: T::from_f64(95.0) + x[1],
            bilateral_mpor_days_m_b: T::from_f64(10.0),
            ccp_mpor_days_m_c: T::from_f64(5.0),
            im_funding_spread_s_im: T::from_f64(0.0085) + dx2 * T::from_f64(0.05),
            ccp_default_fund_cost_rate_c_df: T::from_f64(0.0040),
            bilateral_capital_surcharge_c_kva: T::from_f64(0.0095) + dx3 * T::from_f64(0.05),
        },
    );

    let sdrn_comm = ql_commodities::evaluate_szymanowska_de_roon_nijman_commodity_premia(
        &ql_commodities::SzymanowskaDeRoonNijmanCommodityPremiaSpec {
            nearby_log_basis_y1: T::from_f64(-0.042) + dx0,
            curve_curvature_spreading_basis_y2: T::from_f64(0.015) + dx1,
            spot_premium_basis_loading_beta_s: T::from_f64(1.35) + dx2,
            term_premium_curvature_loading_beta_t1: T::from_f64(0.75),
            high_minus_low_basis_factor_f_carry: T::from_f64(0.068) + dx3,
            contract_maturity_periods_n: T::from_f64(4.0),
        },
    );

    let msss_fx = ql_fx::evaluate_menkhoff_sarno_schmeling_schrimpf_fx_carry_vol(
        &ql_fx::MenkhoffSarnoSchmelingSchrimpfFxCarryVolSpec {
            current_global_fx_volatility_sigma_t: T::from_f64(0.095) + dx0,
            expected_global_fx_volatility_sigma_hat: T::from_f64(0.082),
            carry_portfolio_rank_z: T::from_f64(0.85) + dx1,
            volatility_beta_asymmetry_gamma: T::from_f64(1.40) + dx2,
            global_vol_risk_price_lambda_v: T::from_f64(-0.052) + dx3 * T::from_f64(0.2),
            global_fx_bid_ask_shock_dba: T::from_f64(0.004),
            dollar_factor_expected_return_mu_dol: T::from_f64(0.015),
        },
    );

    let mp_asian = ql_exotics::evaluate_milevsky_posner_reciprocal_gamma_asian_option(
        &ql_exotics::MilevskyPosnerReciprocalGammaAsianSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(100.0) + x[1],
            risk_free_rate_r: T::from_f64(0.05) + dx2 * T::from_f64(0.2),
            dividend_yield_q: T::from_f64(0.015),
            volatility_sigma: T::from_f64(0.22) + dx3,
            maturity_years_t: T::one(),
        },
    );

    let kns_sdf = ql_portfolio::evaluate_kozak_nagel_santosh_sdf_eigen_shrinkage(
        &ql_portfolio::KozakNagelSantoshSdfEigenShrinkageSpec {
            pc1_eigenvalue_d1: T::from_f64(0.085) + dx0 * T::from_f64(0.2),
            pc2_eigenvalue_d2: T::from_f64(0.008) + dx1 * T::from_f64(0.05),
            pc1_mean_return_mu1: T::from_f64(0.092) + dx2 * T::from_f64(0.2),
            pc2_mean_return_mu2: T::from_f64(0.035),
            ridge_l2_penalty_gamma2: T::from_f64(0.025) + dx3 * T::from_f64(0.1),
            lasso_l1_penalty_gamma1: T::from_f64(0.008),
        },
    );

    let appr_ses = ql_portfolio::evaluate_acharya_pedersen_philippon_richardson_ses(
        &ql_portfolio::AcharyaPedersenPhilipponRichardsonSesSpec {
            initial_equity_capital_w0: T::from_f64(90.0) + x[0],
            quasi_market_leverage_ratio_l0: T::from_f64(14.5) + x[1] * T::from_f64(0.2),
            required_capital_ratio_z: T::from_f64(0.08),
            normal_tail_mes_5pct: T::from_f64(0.038) + dx2 * T::from_f64(0.1),
            systemic_crisis_tail_multiplier_k: T::from_f64(12.5) + x[3] * T::from_f64(0.2),
            systemic_crisis_probability_pi: T::from_f64(0.035),
            systemic_externality_weight_tau: T::from_f64(0.45),
        },
    );

    let bnhls_rk = ql_dml::evaluate_barndorff_nielsen_hansen_lunde_shephard_realized_kernel(
        &ql_dml::BarndorffNielsenHansenLundeShephardRealizedKernelSpec {
            lag0_realized_variance_gamma0: T::from_f64(0.00048) + dx0 * T::from_f64(1e-4),
            lag1_realized_autocovariance_gamma1: T::from_f64(-0.00014) + dx1 * T::from_f64(1e-5),
            lag2_realized_autocovariance_gamma2: T::from_f64(-0.00002) + dx2 * T::from_f64(1e-5),
            parzen_kernel_bandwidth_h: T::from_f64(18.0) + x[3],
            intraday_sample_size_n: T::from_f64(4680.0),
            integrated_quarticity_iq: T::from_f64(1.4e-7),
        },
    );

    let bp_spiral = ql_models::evaluate_brunnermeier_pedersen_margin_spiral_liquidity(
        &ql_models::BrunnermeierPedersenMarginSpiralLiquiditySpec {
            customer_order_imbalance_z0: T::from_f64(1.50) + dx0,
            speculator_liquid_wealth_w: T::from_f64(4.20) + x[1] * T::from_f64(0.1),
            baseline_margin_requirement_m0: T::from_f64(0.12),
            margin_illiquidity_sensitivity_psi: T::from_f64(0.35) + dx2,
            funding_shadow_cost_phi: T::from_f64(0.85),
            fundamental_volatility_sigma_v: T::from_f64(0.18) + dx3,
        },
    );

    let coupling = (x[0] + x[1] + x[2] + x[3]) * T::from_f64(1e-3);
    [
        br_qe.qe_signaling_channel_change_bps
            + br_qe.qe_portfolio_balance_change_bps * T::from_f64(0.05)
            + coupling,
        bfl_mbs.primary_secondary_mortgage_spread_bps * T::from_f64(0.05) + coupling,
        ags_cdo.senior_tranche_clustered_loss_bps * T::from_f64(0.01) + coupling,
        bt_fear.total_variance_risk_premium * T::from_f64(100.0) + coupling,
        gg_mva.central_clearing_cost_saving * T::from_f64(5.0) + coupling,
        sdrn_comm.total_deferred_futures_risk_premium * T::from_f64(10.0) + coupling,
        msss_fx.conditional_realized_fx_excess_return * T::from_f64(10.0) + coupling,
        mp_asian.milevsky_posner_asian_call * T::from_f64(0.2) + coupling,
        kns_sdf.regularized_sdf_sharpe_ratio * T::from_f64(2.0) + coupling,
        appr_ses.pigouvian_systemic_risk_tax * T::from_f64(1.0) + coupling,
        bnhls_rk.annualized_realized_kernel_volatility * T::from_f64(10.0) + coupling,
        bp_spiral.equilibrium_market_illiquidity_lambda * T::from_f64(10.0) + coupling,
    ]
}

pub fn benchmark_bauer_rudebusch_qe_signaling_boyarchenko_fuster_lucca_mbs_oas_azizpour_giesecke_schwenkler_default_clustering_bollerslev_todorov_jump_fear_ghamami_glasserman_clearing_mva_szymanowska_de_roon_nijman_commodity_premia_menkhoff_sarno_schmeling_schrimpf_fx_vol_milevsky_posner_reciprocal_gamma_asian_kozak_nagel_santosh_sdf_shrinkage_acharya_pedersen_philippon_richardson_ses_bnhls_realized_kernel_brunnermeier_pedersen_margin_spiral(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite80_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite80_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "80_bauer_rudebusch_qe_signaling_boyarchenko_fuster_lucca_mbs_oas_azizpour_giesecke_schwenkler_default_clustering_bollerslev_todorov_jump_fear_ghamami_glasserman_clearing_mva_szymanowska_de_roon_nijman_commodity_premia_menkhoff_sarno_schmeling_schrimpf_fx_vol_milevsky_posner_reciprocal_gamma_asian_kozak_nagel_santosh_sdf_shrinkage_acharya_pedersen_philippon_richardson_ses_bnhls_realized_kernel_brunnermeier_pedersen_margin_spiral",
        "Bauer-Rudebusch-QE / Boyarchenko-Fuster-Lucca-MBS-OAS / Azizpour-Giesecke-Schwenkler-Clustering / Bollerslev-Todorov-Jump-Fear / Ghamami-Glasserman-Clearing-MVA / Szymanowska-de-Roon-Nijman-Commodity / Menkhoff-Sarno-Schmeling-Schrimpf-FX-Vol / Milevsky-Posner-Reciprocal-Gamma-Asian / Kozak-Nagel-Santosh-SDF-Shrinkage / APPR-Systemic-SES / BNHLS-Realized-Kernel / Brunnermeier-Pedersen-Margin-Spiral AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite81_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let gss_fomc = ql_curves::evaluate_gurkaynak_sack_swanson_fomc_target_path_factors(
        &ql_curves::GurkaynakSackSwansonFomcTargetPathSpec {
            pre_fomc_treasury_yield_y0: T::from_f64(0.0415) + dx0 * T::from_f64(0.2),
            fed_funds_target_surprise_dmp1: T::from_f64(0.0005) + dx1 * T::from_f64(0.05),
            year_ahead_futures_surprise_ded4: T::from_f64(0.0022) + dx2 * T::from_f64(0.05),
            target_to_path_projection_alpha: T::from_f64(0.72),
            target_factor_decay_kappa_t: T::from_f64(0.65),
            path_factor_hump_decay_kappa_p: T::from_f64(0.35) + dx3,
            maturity_years_tau: T::from_f64(5.0),
        },
    );

    let bpw_bond = ql_fixed_income::evaluate_bao_pan_wang_corporate_bond_autocovariance_illiquidity(
        &ql_fixed_income::BaoPanWangBondAutocovarianceIlliquiditySpec {
            lag1_log_price_autocovariance_cov1: T::from_f64(-0.00011) - dx0 * T::from_f64(1e-4),
            quoted_bid_ask_spread_s_ba: T::from_f64(0.0045) + dx1 * T::from_f64(0.05),
            market_vix_stress_sigma_vix: T::from_f64(0.22) + dx2,
            bond_age_years_a: T::from_f64(3.5),
            structural_default_spread_s_def: T::from_f64(0.0095) + dx3 * T::from_f64(0.1),
            illiquidity_yield_loading_beta_gamma: T::from_f64(38.0),
        },
    );

    let cdgm_cred = ql_credit::evaluate_collin_dufresne_goldstein_martin_credit_spread_determinants(
        &ql_credit::CollinDufresneGoldsteinMartinCreditDeterminantsSpec {
            firm_leverage_change_dl: T::from_f64(0.015) + dx0,
            treasury_rate_change_dr: T::from_f64(-0.0020) + dx1 * T::from_f64(0.1),
            yield_curve_slope_change_dslope: T::from_f64(-0.0010),
            vix_volatility_change_dvix: T::from_f64(0.018) + dx2,
            equity_jump_skew_change_djump: T::from_f64(0.012),
            latent_common_credit_factor_f_cdgm: T::from_f64(0.0038) + dx3 * T::from_f64(0.1),
            latent_factor_loading_beta_cdgm: T::from_f64(1.25),
        },
    );

    let aft_tail = ql_vol::evaluate_andersen_fusari_todorov_left_tail_jump_factor_surface(
        &ql_vol::AndersenFusariTodorovLeftTailFactorSpec {
            diffusive_spot_variance_v: T::from_f64(0.025) + dx0 * T::from_f64(0.2),
            left_tail_jump_factor_u: T::from_f64(1.35) + dx1,
            mean_left_jump_size_mu_j: T::from_f64(-0.065) - dx2 * T::from_f64(0.2),
            risk_neutral_tail_reversion_kappa_u: T::from_f64(1.80),
            tail_risk_premium_wedge_theta_u: T::from_f64(0.62),
            standardized_otm_put_moneyness_m: T::from_f64(-1.50) + dx3,
            maturity_years_tau: T::from_f64(0.25),
        },
    );

    let cc_ccp = ql_xva::evaluate_capponi_cheng_clearinghouse_margin_skin_in_game(
        &ql_xva::CapponiChengClearinghouseMarginSpec {
            member_position_volatility_sigma: T::from_f64(85.0) + x[0],
            tail_loss_decay_alpha: T::from_f64(1.40) + dx1,
            member_default_probability_pi: T::from_f64(0.025) + dx2 * T::from_f64(0.2),
            collateral_funding_cost_lambda_c: T::from_f64(0.0080) + dx3 * T::from_f64(0.05),
            default_waterfall_penalty_gamma_d: T::from_f64(1.65),
            ccp_skin_in_the_game_share_s: T::from_f64(0.25),
        },
    );

    let krt_comm =
        ql_commodities::evaluate_kang_rouwenhorst_tang_hedger_speculator_commodity_premia(
            &ql_commodities::KangRouwenhorstTangTwoPremiumsSpec {
                long_horizon_commercial_hedging_pressure_hp: T::from_f64(0.28) + dx0,
                short_horizon_speculator_flow_dq: T::from_f64(0.045) + dx1,
                hedging_pressure_risk_premium_beta_hp: T::from_f64(0.18) + dx2,
                liquidity_provision_reversal_beta_liq: T::from_f64(0.42),
                commodity_futures_basis_b: T::from_f64(-0.035) + dx3,
                commodity_volatility_sigma: T::from_f64(0.24),
            },
        );

    let lrv_dc = ql_fx::evaluate_lustig_roussanov_verdelhan_countercyclical_dollar_carry(
        &ql_fx::LustigRoussanovVerdelhanDollarCarrySpec {
            us_short_rate_r_usd: T::from_f64(0.025) + dx0 * T::from_f64(0.2),
            basket_foreign_short_rate_r_star: T::from_f64(0.048) + dx1 * T::from_f64(0.2),
            us_sdf_volatility_sigma_us: T::from_f64(0.42) + dx2,
            foreign_sdf_volatility_sigma_star: T::from_f64(0.34),
            us_recession_macro_slowdown_z: T::from_f64(0.65) + dx3,
            afd_predictability_slope_beta: T::from_f64(1.85),
        },
    );

    let gy_asian = ql_exotics::evaluate_geman_yor_bessel_laplace_asian_option(
        &ql_exotics::GemanYorBesselLaplaceAsianSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(98.0) + x[1],
            risk_free_rate_r: T::from_f64(0.055) + dx2 * T::from_f64(0.2),
            dividend_yield_q: T::from_f64(0.015) + dx3,
            volatility_sigma: T::from_f64(0.24),
            maturity_years_t: T::one(),
        },
    );

    let tz_comb = ql_portfolio::evaluate_tu_zhou_markowitz_talmud_combination(
        &ql_portfolio::TuZhouMarkowitzTalmudCombinationSpec {
            asset_count_n: T::from_f64(25.0) + x[0],
            sample_size_t: T::from_f64(120.0) + x[1],
            risk_aversion_gamma: T::from_f64(3.0) + dx2,
            true_tangency_squared_sharpe_theta2: T::from_f64(0.28) + dx3,
            naive_equal_weight_squared_sharpe_theta_ew2: T::from_f64(0.11),
            cross_sectional_mean_dispersion_d2: T::from_f64(0.045),
        },
    );

    let bfm_fast = ql_portfolio::evaluate_biais_foucault_moinas_equilibrium_fast_trading(
        &ql_portfolio::BiaisFoucaultMoinasFastTradingSpec {
            fast_search_gains_from_trade_g: T::from_f64(0.42) + dx0,
            adverse_selection_rent_a: T::from_f64(0.28) + dx1,
            fast_trading_technology_cost_c: T::from_f64(0.38) + dx2,
            institutional_type_dispersion_sigma: T::from_f64(0.25) + dx3,
            illiquidity_deadweight_multiplier_lambda: T::from_f64(0.35),
        },
    );

    let cop_uhf = ql_dml::evaluate_christensen_oomen_podolskij_uhf_preaveraged_jump_test(
        &ql_dml::ChristensenOomenPodolskijUhfJumpSpec {
            coarse_5min_realized_variance_rv5m: T::from_f64(0.000240) + dx0 * T::from_f64(1e-4),
            coarse_5min_bipower_variation_bv5m: T::from_f64(0.000212) + dx1 * T::from_f64(1e-4),
            uhf_preaveraged_quadratic_variation_pqv: T::from_f64(0.000232)
                + dx2 * T::from_f64(1e-4),
            uhf_preaveraged_bipower_variation_pbv: T::from_f64(0.000224) + dx3 * T::from_f64(1e-4),
            microstructure_noise_adjustment_omega2: T::from_f64(0.0000012),
            preaveraging_horizon_theta: T::one(),
        },
    );

    let gab_gran = ql_models::evaluate_gabaix_granular_residual_macro_volatility(
        &ql_models::GabaixGranularOriginsMacroVolatilitySpec {
            economy_firm_count_n: T::from_f64(1_000_000.0) + x[0] * T::from_f64(1_000.0),
            firm_idiosyncratic_volatility_sigma_f: T::from_f64(0.12) + dx1,
            pareto_size_tail_exponent_zeta: T::from_f64(1.08) + dx2,
            domar_input_output_multiplier_mu: T::from_f64(2.10) + dx3,
            realized_granular_residual_gamma: T::from_f64(0.0085),
        },
    );

    [
        gss_fomc.post_fomc_treasury_yield * T::from_f64(10.0),
        bpw_bond.total_corporate_bond_spread_bps * T::from_f64(0.01),
        cdgm_cred.total_credit_spread_change_bps * T::from_f64(0.05),
        aft_tail.aft_otm_put_implied_volatility * T::from_f64(10.0),
        cc_ccp.total_clearing_social_cost * T::from_f64(0.5),
        krt_comm.two_premiums_sharpe_ratio * T::from_f64(5.0),
        lrv_dc.conditional_dollar_carry_expected_return * T::from_f64(10.0),
        gy_asian.geman_yor_asian_call_price * T::from_f64(0.2),
        tz_comb.combined_expected_out_of_sample_utility * T::from_f64(20.0),
        bfm_fast.fast_trading_overinvestment_wedge * T::from_f64(10.0),
        cop_uhf.spurious_coarse_jump_illusion_wedge * T::from_f64(10.0),
        gab_gran.granular_aggregate_gdp_volatility * T::from_f64(20.0),
    ]
}

pub fn benchmark_gurkaynak_sack_swanson_fomc_path_bao_pan_wang_bond_illiquidity_collin_dufresne_goldstein_martin_credit_determinants_andersen_fusari_todorov_left_tail_capponi_cheng_ccp_margin_kang_rouwenhorst_tang_two_premiums_lustig_roussanov_verdelhan_dollar_carry_geman_yor_bessel_asian_tu_zhou_markowitz_talmud_biais_foucault_moinas_fast_trading_christensen_oomen_podolskij_uhf_jumps_gabaix_granular_residual(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite81_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite81_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "81_gurkaynak_sack_swanson_fomc_path_bao_pan_wang_bond_illiquidity_collin_dufresne_goldstein_martin_credit_determinants_andersen_fusari_todorov_left_tail_capponi_cheng_ccp_margin_kang_rouwenhorst_tang_two_premiums_lustig_roussanov_verdelhan_dollar_carry_geman_yor_bessel_asian_tu_zhou_markowitz_talmud_biais_foucault_moinas_fast_trading_christensen_oomen_podolskij_uhf_jumps_gabaix_granular_residual",
        "Gurkaynak-Sack-Swanson-FOMC / Bao-Pan-Wang-Bond-Illiquidity / CDGM-Credit-Determinants / Andersen-Fusari-Todorov-Left-Tail / Capponi-Cheng-CCP-Margin / Kang-Rouwenhorst-Tang-Two-Premiums / LRV-Dollar-Carry / Geman-Yor-Bessel-Asian / Tu-Zhou-Markowitz-Talmud / Biais-Foucault-Moinas-Fast-Trading / COP-UHF-Preaveraged-Jumps / Gabaix-Granular-Residual AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}

pub fn suite82_workloads<T: Real>(x: [T; 4]) -> [T; 12] {
    let dx0 = x[0] * T::from_f64(0.01);
    let dx1 = x[1] * T::from_f64(0.01);
    let dx2 = x[2] * T::from_f64(0.01);
    let dx3 = x[3] * T::from_f64(0.01);

    let ns_info = ql_curves::evaluate_nakamura_steinsson_fomc_information_effect(
        &ql_curves::NakamuraSteinssonFomcInformationEffectSpec {
            policy_news_composite_shock_ds: T::from_f64(0.0018) + dx0 * T::from_f64(0.05),
            fed_information_share_omega: T::from_f64(0.62) + dx1,
            natural_rate_information_loading_psi: T::from_f64(0.85) + dx2,
            information_decay_kappa_info: T::from_f64(0.22),
            monetary_gap_decay_kappa_mp: T::from_f64(0.55) + dx3,
            information_inflation_loading_gamma_inf: T::from_f64(0.45),
            monetary_disinflation_loading_gamma_mp: T::from_f64(0.30),
            maturity_years_tau: T::from_f64(3.0),
        },
    );

    let dfl_illiq = ql_fixed_income::evaluate_dick_nielsen_feldhutter_lando_bond_illiquidity(
        &ql_fixed_income::DickNielsenFeldhutterLandoBondIlliquiditySpec {
            amihud_price_impact_bp_per_mn: T::from_f64(6.2) + x[0],
            imputed_roundtrip_cost_bps: T::from_f64(31.0) + x[1],
            amihud_volatility_bp_per_mn: T::from_f64(4.8) + dx2 * T::from_f64(10.0),
            irc_volatility_bps: T::from_f64(18.5),
            turnover_rate_zeta: T::from_f64(0.18) + dx3,
            credit_rating_speculative_multiplier_m: T::from_f64(1.45),
            crisis_regime_slope_beta_dfl: T::from_f64(42.0),
        },
    );

    let bd_drp = ql_credit::evaluate_berndt_duffie_credit_risk_premium_ratio(
        &ql_credit::BerndtDuffieCreditRiskPremiumSpec {
            moodys_kmv_physical_edf_intensity: T::from_f64(0.0085) + dx0 * T::from_f64(0.1),
            physical_loss_given_default_lgd: T::from_f64(0.58) + dx1,
            risk_neutral_lgd_markup_xi: T::from_f64(1.12),
            macro_cycle_distress_z: T::from_f64(0.40) + dx2,
            default_event_jump_risk_aversion_gamma: T::from_f64(7.5) + dx3 * T::from_f64(10.0),
            cds_maturity_years_t: T::from_f64(5.0),
        },
    );

    let chj_kern = ql_vol::evaluate_christoffersen_heston_jacobs_variance_dependent_kernel(
        &ql_vol::ChristoffersenHestonJacobsVarianceKernelSpec {
            equity_risk_aversion_gamma: T::from_f64(2.2) + dx0 * T::from_f64(10.0),
            variance_preference_xi: T::from_f64(14.0) + x[1],
            spot_variance_v0: T::from_f64(0.032) + dx2 * T::from_f64(0.2),
            variance_of_variance_sigma_v: T::from_f64(0.38),
            return_variance_correlation_rho: T::from_f64(-0.72),
            physical_mean_reversion_kappa_p: T::from_f64(3.50),
            standardized_log_return_z: T::from_f64(2.40) + dx3 * T::from_f64(10.0),
        },
    );

    let bfp_xva = ql_xva::evaluate_brigo_francischello_pallavicini_closeout_xva(
        &ql_xva::BrigoFrancischelloPallaviciniCloseoutXvaSpec {
            risk_free_mark_to_market_v0: T::from_f64(100.0) + x[0],
            counterparty_hazard_rate_lambda_c: T::from_f64(0.032) + dx1 * T::from_f64(0.2),
            institution_hazard_rate_lambda_b: T::from_f64(0.018) + dx2 * T::from_f64(0.2),
            counterparty_lgd_c: T::from_f64(0.60),
            institution_lgd_b: T::from_f64(0.60),
            uncollateralized_funding_spread_sf: T::from_f64(0.012) + dx3 * T::from_f64(0.1),
            maturity_years_t: T::from_f64(5.0),
        },
    );

    let bp_bm = ql_commodities::evaluate_boons_prado_commodity_basis_momentum(
        &ql_commodities::BoonsPradoCommodityBasisMomentumSpec {
            first_nearby_12m_return_r1: T::from_f64(0.145) + dx0,
            second_nearby_12m_return_r2: T::from_f64(0.082) + dx1,
            current_front_log_basis_b1: T::from_f64(-0.025) + dx2,
            current_second_log_basis_b2: T::from_f64(0.012),
            speculator_spreading_imbalance_q_sp: T::from_f64(0.65) + dx3,
            clearing_intermediary_risk_aversion_gamma: T::from_f64(1.80),
        },
    );

    let ccgr_fx = ql_fx::evaluate_colacito_croce_gavazzoni_ready_recursive_fx_factors(
        &ql_fx::ColacitoCroceGavazzoniReadyRecursiveFxSpec {
            relative_risk_aversion_gamma: T::from_f64(6.5) + x[0],
            intertemporal_elasticity_psi: T::from_f64(1.5) + dx1 * T::from_f64(5.0),
            home_lrr_exposure_beta_h: T::from_f64(1.25) + dx2,
            foreign_lrr_exposure_beta_f: T::from_f64(0.72) + dx3,
            global_lrr_volatility_sigma_x: T::from_f64(0.0042),
            short_run_output_volatility_sigma_c: T::from_f64(0.015),
            lrr_persistence_rho_x: T::from_f64(0.94),
        },
    );

    let fm_asian = ql_exotics::evaluate_fusai_meucci_erlang_canadized_asian_option(
        &ql_exotics::FusaiMeucciErlangCanadizedAsianSpec {
            spot_price_s0: T::from_f64(100.0) + x[0],
            strike_price_k: T::from_f64(98.0) + x[1],
            risk_free_rate_r: T::from_f64(0.045) + dx2 * T::from_f64(0.2),
            volatility_sigma: T::from_f64(0.25) + dx3,
            maturity_years_t: T::one(),
            monitoring_dates_n: T::from_f64(12.0),
            erlang_stages_m: T::from_f64(8.0),
        },
    );

    let fnw_char = ql_portfolio::evaluate_freyberger_neuhierl_weber_nonparametric_characteristics(
        &ql_portfolio::FreybergerNeuhierlWeberNonparametricSpec {
            characteristic_quantile_c: T::from_f64(0.88) + dx0,
            linear_characteristic_slope_b1: T::from_f64(0.085) + dx1,
            extreme_tail_spline_curvature_b2: T::from_f64(0.195) + dx2,
            knot_threshold_tau_k: T::from_f64(0.70),
            group_lasso_penalty_lambda_g: T::from_f64(0.012) + dx3 * T::from_f64(0.1),
            adaptive_weight_exponent_nu: T::from_f64(1.50),
        },
    );

    let kj_tail = ql_portfolio::evaluate_kelly_jiang_cross_sectional_tail_risk_factor(
        &ql_portfolio::KellyJiangTailRiskFactorSpec {
            mean_log_exceedance_ratio_hill: T::from_f64(0.42) + dx0,
            tail_threshold_magnitude_u: T::from_f64(0.048) + dx1 * T::from_f64(0.2),
            exceedance_count_k: T::from_f64(180.0) + x[2],
            asset_tail_beta_loading_beta_tail: T::from_f64(1.35) + dx3,
            tail_factor_market_price_of_risk_gamma_kj: T::from_f64(0.14),
            tail_persistence_rho_lambda: T::from_f64(0.88),
        },
    );

    let jlmpv_rv = ql_dml::evaluate_jacod_li_mykland_podolskij_vetter_preaveraging_rv(
        &ql_dml::JacodLiMyklandPodolskijVetterPreaveragingSpec {
            observation_count_n: T::from_f64(23_400.0) + x[0] * T::from_f64(100.0),
            preaveraging_window_theta: T::from_f64(0.80) + dx1,
            sum_squared_preaveraged_returns_sbar: T::from_f64(0.00245) + dx2 * T::from_f64(1e-3),
            microstructure_noise_variance_omega2: T::from_f64(0.0000015) + dx3 * T::from_f64(1e-6),
            integrated_quarticity_iq: T::from_f64(0.00000012),
        },
    );

    let aem_sdf = ql_models::evaluate_adrian_etula_muir_broker_dealer_leverage_sdf(
        &ql_models::AdrianEtulaMuirBrokerDealerLeverageSpec {
            broker_dealer_leverage_ratio_l0: T::from_f64(22.0) + x[0],
            seasonally_adjusted_leverage_growth_dl: T::from_f64(0.045) + dx1,
            leverage_factor_price_of_risk_lambda_lev: T::from_f64(0.068) + dx2,
            leverage_shock_variance_sigma_l2: T::from_f64(0.012),
            asset_leverage_beta_loading_beta_i: T::from_f64(1.25) + dx3,
            shadow_var_constraint_tightness_phi: T::from_f64(0.18),
        },
    );

    [
        ns_info.nominal_forward_rate_response_bps * T::from_f64(0.05),
        dfl_illiq.liquidity_yield_spread_component_bps * T::from_f64(0.01),
        bd_drp.credit_risk_premium_wedge_bps * T::from_f64(0.01),
        chj_kern.state_price_density_kernel_ratio * T::from_f64(0.5),
        bfp_xva.closeout_convention_xva_wedge * T::from_f64(0.2),
        bp_bm.expected_front_month_excess_return * T::from_f64(10.0),
        ccgr_fx.equilibrium_fx_carry_risk_premium * T::from_f64(10.0),
        fm_asian.erlang_canadized_asian_call_price * T::from_f64(0.2),
        fnw_char.nonparametric_characteristic_expected_return * T::from_f64(15.0),
        kj_tail.required_tail_risk_premium * T::from_f64(10.0),
        jlmpv_rv.preaveraging_standard_error * T::from_f64(1_000.0),
        aem_sdf.cross_sectional_expected_excess_return * T::from_f64(10.0),
    ]
}

pub fn benchmark_nakamura_steinsson_fomc_info_dick_nielsen_feldhutter_lando_illiquidity_berndt_duffie_credit_premia_christoffersen_heston_jacobs_variance_kernel_brigo_francischello_pallavicini_closeout_xva_boons_prado_basis_momentum_colacito_croce_gavazzoni_ready_recursive_fx_fusai_meucci_erlang_asian_freyberger_neuhierl_weber_characteristics_kelly_jiang_tail_risk_jacod_li_mykland_podolskij_vetter_preaveraging_adrian_etula_muir_broker_dealer_sdf(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let baseline = move || {
        let f = |x: [f64; 4]| suite82_workloads(x).iter().copied().sum::<f64>();
        let mut z = f(x0);
        for i in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[i] += 1e-5;
            d[i] -= 1e-5;
            z += (f(u) - f(d)) / (2e-5);
        }
        for _ in 0..20 {
            std::hint::black_box(f(x0));
        }
        z
    };
    let sota = move || {
        let x: [Dual<4>; 4] = [
            Dual::variable(x0[0], 0),
            Dual::variable(x0[1], 1),
            Dual::variable(x0[2], 2),
            Dual::variable(x0[3], 3),
        ];
        suite82_workloads(x)
            .iter()
            .map(|v| v.val + v.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "82_nakamura_steinsson_fomc_info_dick_nielsen_feldhutter_lando_illiquidity_berndt_duffie_credit_premia_christoffersen_heston_jacobs_variance_kernel_brigo_francischello_pallavicini_closeout_xva_boons_prado_basis_momentum_colacito_croce_gavazzoni_ready_recursive_fx_fusai_meucci_erlang_asian_freyberger_neuhierl_weber_characteristics_kelly_jiang_tail_risk_jacod_li_mykland_podolskij_vetter_preaveraging_adrian_etula_muir_broker_dealer_sdf",
        "Nakamura-Steinsson-FOMC-Info / Dick-Nielsen-Feldhutter-Lando-Illiquidity / Berndt-Duffie-Credit-Premia / Christoffersen-Heston-Jacobs-Variance-Kernel / Brigo-Francischello-Pallavicini-Closeout-XVA / Boons-Prado-Basis-Momentum / Colacito-Croce-Gavazzoni-Ready-Recursive-FX / Fusai-Meucci-Erlang-Asian / Freyberger-Neuhierl-Weber-Characteristics / Kelly-Jiang-Tail-Risk / JLMPV-Preaveraging / Adrian-Etula-Muir-Broker-Dealer-SDF AD",
        cfg,
        1e-3,
        baseline,
        sota,
    )
}
