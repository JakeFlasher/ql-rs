use crate::harness::{
    BenchConfig, ComparisonRow, configured_counts, run_comparison, with_config_override,
};
use ql_ad::{Dual, Dual2, Real, TapeVar};
use ql_commodities as _;
use ql_credit::asb_recursive_expected_tranche_loss;
use ql_credit::{
    CalibratedCreditCurve, CdsQuote, CreditDefaultSwap, IsdaCdsEngine, ProtectionSide,
    SurvivalProbabilityCurve,
};
use ql_curves::{
    CalibratedCurveSystem, DiscountCurve, FlatForwardCurve, SingleCurveHelper, YieldTermStructure,
};
use ql_dml::{cholesky_lower, price_correlated_basket_call_aad};
use ql_exotics::{BarrierOptionSpec, BarrierType, MonitoringFrequency, price_barrier_option};
use ql_fx::{CalibratedXccsCurve, ConstNotionalXccs, XccsQuote};
use ql_fx::{FxSmileQuote, FxVannaVolgaSmile};
use ql_inflation::ZeroInflationCurve;
use ql_instruments::{FixedRateLeg, OptionType, ResolvedCouponLeg};
use ql_levy::{VarianceGamma, cos_levy_price_strip};
use ql_mc::{GbmSoaSimulator, PhiloxRng, price_asian_arithmetic_rqmc};
use ql_models::HullWhite1F;
use ql_models::{HestonParams, heston_cos_price};
use ql_portfolio::{MarketSnapshot, PortfolioDag, PortfolioTrade, risk_budgeting_weights};
use ql_rates as _;
use ql_rates::{
    RfrCapletSpec, RfrRateType, willems_backward_sabr_caplet_price,
    willems_backward_sabr_caplet_strip,
};
use ql_scripting::{BytecodeProgram, ScriptMcEngine, ScriptedProduct};
use ql_time::{BitsetCalendar, Date, DayCounter, Frequency};
use ql_vol::{black76_undiscounted, implied_black_vol, implied_black_vol_dual2};
use ql_xva::{NettingSetSwap, compute_bruteforce_dim_and_mva, compute_chebyshev_dim_and_mva};
use std::sync::Arc;

fn config(scale: usize) -> BenchConfig {
    let (warmup_iters, samples) = configured_counts(1, 31);
    BenchConfig {
        warmup_iters,
        samples,
        // Repeat each fixed workload enough times to amortize timer and
        // scheduler jitter when the workspace's integration binaries run in
        // parallel.  `scale` changes only this repetition count; the domain
        // workload itself remains invariant across test and bench modes.
        // The short parity configuration intentionally uses one iteration so
        // the complete 20-suite test remains a fast (<2s on CI) invariant
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
        let mut s_pos = [spot; 8];
        let mut s_neg = [spot; 8];
        for step in 0..steps {
            let base_path = (batch * 2) as u64;
            let z0 = rng.normal4(base_path, step as u32);
            let z1 = rng.normal4(base_path + 1, step as u32);
            let z = [z0[0], z0[1], z0[2], z0[3], z1[0], z1[1], z1[2], z1[3]];
            for lane in 0..8 {
                s_pos[lane] *= libm::exp(drift + vol_dt * z[lane]);
                s_neg[lane] *= libm::exp(drift - vol_dt * z[lane]);
            }
        }
        for lane in 0..8 {
            sum += ((s_pos[lane] - strike).max(0.0) + (s_neg[lane] - strike).max(0.0))
                * (0.5 * discount);
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
                let h = 1e-5;
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
                let h = 1e-5;
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
                for &k in &strikes {
                    total += ql_levy::cos_levy_price(
                        &base,
                        params[0],
                        k,
                        params[4],
                        0.0,
                        1.0,
                        OptionType::Call,
                        128,
                    );
                }
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
                    for &k in &strikes {
                        let vp = ql_levy::cos_levy_price(
                            &mup,
                            up[0],
                            k,
                            up[4],
                            0.0,
                            1.0,
                            OptionType::Call,
                            128,
                        );
                        let vm = ql_levy::cos_levy_price(
                            &mdn,
                            down[0],
                            k,
                            down[4],
                            0.0,
                            1.0,
                            OptionType::Call,
                            128,
                        );
                        total += (vp - vm) / (2.0 * h);
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
        1e-8,
        {
            let cov = cov.clone();
            let budgets = budgets.clone();
            move || {
                let h = 1e-4;
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
                    let d1 = derivative(h);
                    let d2 = derivative(h * 0.5);
                    checksum += (4.0 * d2 - d1) / 3.0;
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
                        let d1 = derivative(h);
                        let d2 = derivative(h * 0.5);
                        checksum += (4.0 * d2 - d1) / 3.0;
                    }
                }
                std::hint::black_box(checksum)
            }
        },
        move || {
            use ql_ad::Dual;
            const LANES: usize = 90;
            let mut lane = 0;
            let cov_ad = (0..n)
                .map(|i| {
                    (0..n)
                        .map(|j| {
                            let x = Dual::<LANES>::variable(cov[i][j], lane);
                            if j >= i {
                                lane += 1;
                            }
                            x
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            // The upper-triangle lane assignment above is shared by the
            // symmetric inputs by explicitly seeding a second pass below.
            let mut cov_ad = cov_ad;
            lane = 0;
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

pub fn run_all_benchmarks(scale: usize) -> Vec<ComparisonRow> {
    // Keep all quantitative engine dependencies linked into the 20-suite
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
    ]
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
