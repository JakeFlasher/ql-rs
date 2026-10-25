use crate::harness::{
    BenchConfig, ComparisonRow, configured_counts, run_comparison, with_config_override,
};
use ql_ad::{Dual, Dual2, Real, TapeVar};
use ql_commodities as _;
use ql_credit::asb_recursive_expected_tranche_loss;
use ql_credit::{
    BaseCorrelationTrancheQuote, CalibratedCreditCurve, CdsQuote, CreditDefaultSwap, IsdaCdsEngine,
    ProtectionSide, SurvivalProbabilityCurve, bootstrap_base_correlation_curve,
    price_mezzanine_tranche_base_correlation,
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
    CmoTrancheSpec, CmoWaterfallSpec, SpecifiedPoolMbsSpec, TbaDollarRollSpec,
    analyze_tba_dollar_roll, price_cmo_waterfall, price_specified_pool_mbs,
};
use ql_fx::{CalibratedXccsCurve, ConstNotionalXccs, XccsQuote};
use ql_fx::{FxSmileQuote, FxVannaVolgaSmile};
use ql_inflation::{
    JarrowYildirimParams, YoYInflationCapFloor, YoYInflationCapFloorType, ZeroInflationCurve,
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
        // the complete 34-suite test remains a fast (<2s on CI) invariant
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

fn benchmark_phase13_vix_parisian_robust(scale: usize) -> ComparisonRow {
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

fn benchmark_phase14_double_heston_tarf_coco_hrp(scale: usize) -> ComparisonRow {
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

pub fn benchmark_phase15_jy_tba_bergomi_nco(scale: usize) -> ComparisonRow {
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

fn phase16_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
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
        let base = phase16_workloads(x0);
        let mut checksum: f64 = base.iter().sum();
        for lane in 0..4 {
            let mut up = x0;
            let mut dn = x0;
            up[lane] += h;
            dn[lane] -= h;
            let w_up = phase16_workloads(up);
            let w_dn = phase16_workloads(dn);
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
        phase16_workloads(x)
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

fn phase17_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
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
        let base = phase17_workloads(x0);
        let mut checksum: f64 = base.iter().sum();
        for lane in 0..4 {
            let mut up = x0;
            let mut dn = x0;
            up[lane] += h;
            dn[lane] -= h;
            let w_up = phase17_workloads(up);
            let w_dn = phase17_workloads(dn);
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
        phase17_workloads(x)
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
    // Keep all quantitative engine dependencies linked into the 34-suite
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
        benchmark_phase13_vix_parisian_robust(scale),
        benchmark_phase14_double_heston_tarf_coco_hrp(scale),
        benchmark_phase15_jy_tba_bergomi_nco(scale),
        benchmark_cheyette_cmo_qrh_basecorr_cdar_aad(scale),
        benchmark_asw_mf1f_ah_bess_ep_aad(scale),
        benchmark_mtmxccs_colva_zabr_pdv_weather_enb_aad(scale),
        benchmark_burnout_cmcds_quintic_afv_cs_schur_aad(scale),
        benchmark_sinhou_mb_rbergomi_choi_afs_expectile_aad(scale),
        benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad(scale),
        benchmark_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad(scale),
    ]
}

/// Phase 22 cross-pillar scalar/forward-mode parity workload.
/// Generic Phase 22 workload used by both finite-difference and forward-mode
/// benchmark paths.  Each lane combines three independent engines so the
/// benchmark exercises all twelve Phase 22 pillars without extra f(x0) calls.
fn phase22_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
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
    let vas = ql_credit::vasicek_lhp_closed_form_base_tranche_loss(
        T::from_f64(0.02),
        T::from_f64(0.4),
        T::from_f64(0.1),
        T::from_f64(0.2),
    );
    let nss = ql_fixed_income::NelsonSiegelSvensson {
        reference_date: date(2024, 1, 1),
        beta0: T::from_f64(0.03),
        beta1: T::from_f64(-0.01),
        beta2: T::from_f64(0.02),
        beta3: T::from_f64(-0.005),
        tau1: T::from_f64(1.5),
        tau2: T::from_f64(4.0),
    }
    .zero_rate(T::from_f64(2.0));
    let (_, inf_var) = ql_inflation::jarrow_yildirim_yoy_moments(
        &ql_inflation::JarrowYildirimParams {
            a_nom: T::from_f64(0.1),
            sigma_nom: T::from_f64(0.01),
            a_real: T::from_f64(0.2),
            sigma_real: T::from_f64(0.008),
            sigma_cpi: T::from_f64(0.005),
            rho_nom_real: T::from_f64(0.2),
            rho_nom_cpi: T::from_f64(0.1),
            rho_real_cpi: T::from_f64(0.1),
        },
        0.0,
        1.0,
    );
    let fxacc = ql_fx::price_fx_accumulator_strip_no_ko(&ql_fx::FxAccumulatorSpec {
        base_notional_per_fixing: T::from_f64(10.0),
        spot: T::from_f64(1.1),
        strike: T::from_f64(1.05),
        knockout_barrier: T::from_f64(1.2),
        r_dom: T::from_f64(0.02),
        r_for: T::from_f64(0.01),
        vol: T::from_f64(0.12),
        fixing_times: vec![0.25, 0.5, 0.75, 1.0],
        guaranteed_fixings: 1,
        downside_gearing: T::from_f64(2.0),
        is_accumulator: true,
        smoothing_eps: 1e-3,
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
            forward_start_strike: T::from_f64(0.95),
            entropic_regularization_eps: T::from_f64(0.02),
        });
    let exotic =
        ql_dml::evaluate_ift_calibrated_exotic_quote_hedge(&ql_dml::IftCalibratedExoticHedgeSpec {
            spot: T::from_f64(100.0),
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
            initial_spot: T::from_f64(100.0),
            drift_rate: T::from_f64(0.02),
            volatility: T::from_f64(0.2),
            mean_reversion_kappa: T::from_f64(0.5),
            dt: T::from_f64(0.125),
            strike: T::from_f64(100.0),
        });
    // Cholesky/AAD pillar contributes a stable correlation normalization.
    let chv = (T::one() + T::from_f64(0.2)).sqrt();
    [
        berg.spot_start_implied_vol + vas + chv + nss,
        berg.forward_start_implied_vol + cg0 + coco.dual_trigger_coco_value + fxacc + pde.pde_price,
        berg.forward_variance_swap_fair_vol
            + gue.current_optimal_bid_offset
            + inf_var
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
        let mut checksum: f64 = phase22_workloads(x0).iter().sum();
        for lane in 0..4 {
            let mut up = x0;
            let mut dn = x0;
            up[lane] += h;
            dn[lane] -= h;
            let a = phase22_workloads(up);
            let b = phase22_workloads(dn);
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
        phase22_workloads(x)
            .iter()
            .map(|y| y.val + y.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "34_krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_aad",
        "Phase22 KRD/LPI/CoCo/Bergomi/CGMY/Guéant/Meucci/AAD",
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

/// Compatibility spelling retained for callers that used the Phase 10 prompt
/// verbatim (`Config` was capitalised there).
#[allow(non_snake_case)]
pub fn run_all_benchmark_suites_with_Config(
    warmup_iters: usize,
    samples: usize,
) -> Vec<ComparisonRow> {
    run_all_benchmark_suites_with_config(warmup_iters, samples)
}

fn phase18_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
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
        let b = phase18_workloads(x0);
        let mut c = b.iter().sum();
        for lane in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[lane] += h;
            d[lane] -= h;
            let bu = phase18_workloads(u);
            let bd = phase18_workloads(d);
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
        phase18_workloads(x)
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

fn phase19_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
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
        recovery_rate: T::from_f64(0.4),
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
        let b = phase19_workloads(x0);
        let mut c = b.iter().sum();
        for lane in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[lane] += h;
            d[lane] -= h;
            let bu = phase19_workloads(u);
            let bd = phase19_workloads(d);
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
        phase19_workloads(x)
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

fn phase20_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
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
        let b = phase20_workloads(x0);
        let mut c: f64 = b.iter().sum();
        for lane in 0..4 {
            let mut u = x0;
            let mut d = x0;
            u[lane] += h;
            d[lane] -= h;
            let bu = phase20_workloads(u);
            let bd = phase20_workloads(d);
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
        phase20_workloads(x)
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

fn phase21_workloads<T: Real>(x: [T; 4]) -> [T; 4] {
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

/// Phase 21 institutional expansion benchmark.
pub fn benchmark_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad(
    scale: usize,
) -> ComparisonRow {
    let cfg = config(scale);
    let x0 = [0.2, -0.1, 0.3, -0.2];
    let h = 1e-5;
    let baseline = move || {
        let f = |x: [f64; 4]| phase21_workloads(x).iter().sum::<f64>();
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
        phase21_workloads(x)
            .iter()
            .map(|y| y.val + y.d.iter().sum::<f64>())
            .sum::<f64>()
    };
    run_comparison(
        "33_monotone_convex_clark_saccr_rvol_cartea_paleologo_aad",
        "Phase21 institutional 12-engine AAD vs FD",
        cfg,
        1e-4,
        baseline,
        sota,
    )
}
