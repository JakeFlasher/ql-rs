use ql_ad::Real;
use ql_curves::{
    DiscountCurve, DiscountCurveId, ForwardIndexId, MarketContext, YieldTermStructure,
};
use ql_instruments::{ResolvedVanillaSwap, SwapType, VanillaSwapSpec};
use ql_models::{HestonParams, heston_cos_price};
use ql_scripting::{BytecodeProgram, OpCode};
use ql_time::{Date, DayCounter, Frequency};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

/// Return a stable representation of an IEEE-754 value for content hashes.
///
/// Rust's `f64::to_bits` preserves the sign bit of zero and all NaN payloads,
/// which is useful for bit-level work but not for cache keys: `-0.0` and
/// `0.0` have the same numerical meaning and two NaNs should not create two
/// independent market states.  Keep this helper as the single canonicalisation
/// point used by every portfolio hash.
#[inline]
pub fn canonical_f64_bits(x: f64) -> u64 {
    if x == 0.0 {
        0.0_f64.to_bits()
    } else if x.is_nan() {
        0x7ff8_0000_0000_0000_u64
    } else {
        x.to_bits()
    }
}

pub fn hash_f64_slice(v: &[f64]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &x in v {
        h ^= canonical_f64_bits(x);
        h = h.wrapping_mul(0x0100000001b3);
    }
    h
}

#[inline]
fn hash_word(mut h: u64, word: u64) -> u64 {
    h ^= word;
    h.wrapping_mul(0x0100000001b3)
}

fn hash_program(program: &BytecodeProgram) -> u64 {
    let mut h = hash_f64_slice(&program.times);
    h = hash_word(h, program.num_vars as u64);
    for op in &program.ops {
        let tag = match op {
            OpCode::BeginStep { .. } => 0,
            OpCode::LoadConst(_) => 1,
            OpCode::LoadSpot => 2,
            OpCode::LoadVar(_) => 3,
            OpCode::StoreVar(_) => 4,
            OpCode::Add => 5,
            OpCode::Sub => 6,
            OpCode::Mul => 7,
            OpCode::Div => 8,
            OpCode::Max => 9,
            OpCode::Min => 10,
            OpCode::SmoothStep(_) => 11,
            OpCode::Pay { .. } => 12,
            OpCode::Autocall { .. } => 13,
        };
        h = hash_word(h, tag);
        match op {
            OpCode::BeginStep { step_idx }
            | OpCode::LoadVar(step_idx)
            | OpCode::StoreVar(step_idx) => h = hash_word(h, *step_idx as u64),
            OpCode::LoadConst(x) | OpCode::SmoothStep(x) => {
                h = hash_word(h, canonical_f64_bits(*x))
            }
            OpCode::Pay { discount_time } => h = hash_word(h, canonical_f64_bits(*discount_time)),
            OpCode::Autocall {
                barrier,
                half_width,
                discount_time,
            } => {
                h = hash_word(h, canonical_f64_bits(*barrier));
                h = hash_word(h, canonical_f64_bits(*half_width));
                h = hash_word(h, canonical_f64_bits(*discount_time));
            }
            OpCode::LoadSpot
            | OpCode::Add
            | OpCode::Sub
            | OpCode::Mul
            | OpCode::Div
            | OpCode::Max
            | OpCode::Min => {}
        }
    }
    h
}

fn trade_dependency_hash(trade: &PortfolioTrade, curve_hash: u64, equity_hash: u64) -> u64 {
    let h = match trade {
        PortfolioTrade::IrsSwap {
            notional,
            fixed_rate,
            maturity_years,
            ..
        } => hash_f64_slice(&[*notional, *fixed_rate, *maturity_years as f64]),
        PortfolioTrade::HestonEuropeanOption {
            strike,
            expiry,
            is_call,
            notional,
            ..
        } => hash_f64_slice(&[
            *strike,
            *expiry,
            *notional,
            if *is_call { 1.0 } else { 0.0 },
        ]),
        PortfolioTrade::ScriptedNote {
            program, expiry, ..
        } => hash_word(hash_program(program), canonical_f64_bits(*expiry)),
    };
    match trade {
        PortfolioTrade::IrsSwap { .. } => hash_word(curve_hash, h),
        PortfolioTrade::HestonEuropeanOption { .. } | PortfolioTrade::ScriptedNote { .. } => {
            hash_word(curve_hash ^ equity_hash.rotate_left(13), h)
        }
    }
}

/// Public content hash for a trade under a market snapshot.  Keeping this
/// helper next to the dependency hash makes it possible for callers to build
/// external cache keys using the same canonical float treatment as the DAG.
pub fn trade_market_hash(trade: &PortfolioTrade, market: &MarketSnapshot) -> u64 {
    trade_dependency_hash(trade, market.curve_hash(), market.equity_hash())
}
#[derive(Clone, Debug)]
pub struct MarketSnapshot {
    pub ref_date: Date,
    pub ois_quotes: [f64; 3],
    pub fwd_quotes: [f64; 3],
    pub heston_params: HestonParams<f64>,
    pub eq_spot: f64,
    pub eq_vol: f64,
}
impl MarketSnapshot {
    pub fn curve_hash(&self) -> u64 {
        hash_word(
            hash_f64_slice(&[self.ois_quotes.as_slice(), self.fwd_quotes.as_slice()].concat()),
            self.ref_date.0 as u32 as u64,
        )
    }
    pub fn equity_hash(&self) -> u64 {
        hash_f64_slice(&[
            self.heston_params.v0,
            self.heston_params.kappa,
            self.heston_params.theta,
            self.heston_params.sigma_v,
            self.heston_params.rho,
            self.eq_spot,
            self.eq_vol,
        ])
    }

    /// Materialize the snapshot in the shared market-context representation
    /// used by portfolio and pricing services.
    pub fn to_market_context(&self) -> MarketContext<f64> {
        let times = vec![1.0, 2.0, 5.0];
        let ois = self
            .ois_quotes
            .iter()
            .zip(times.iter())
            .map(|(&r, &t)| (-r * t).exp())
            .collect();
        let fwd = self
            .fwd_quotes
            .iter()
            .zip(times.iter())
            .map(|(&r, &t)| (-r * t).exp())
            .collect();
        MarketContext::new()
            .with_reference_date(self.ref_date)
            .with_discount_curve(
                DiscountCurveId::USD_SOFR,
                DiscountCurve::new(self.ref_date, times.clone(), ois),
            )
            .with_forward_curve(
                ForwardIndexId::USD_SOFR,
                DiscountCurve::new(self.ref_date, times, fwd),
            )
            .with_spot("EQ", self.eq_spot)
            .with_vol("EQ", self.eq_vol)
    }
}
#[derive(Clone, Debug)]
pub enum PortfolioTrade {
    IrsSwap {
        id: usize,
        notional: f64,
        fixed_rate: f64,
        maturity_years: usize,
    },
    HestonEuropeanOption {
        id: usize,
        strike: f64,
        expiry: f64,
        is_call: bool,
        notional: f64,
    },
    ScriptedNote {
        id: usize,
        program: BytecodeProgram,
        expiry: f64,
    },
}
#[derive(Clone, Debug)]
pub struct TradeValuation {
    pub id: usize,
    pub npv: f64,
    pub par_dv01: Vec<f64>,
    pub eq_delta: f64,
}
#[derive(Clone, Debug, Default)]
pub struct PortfolioRiskSummary {
    pub total_npv: f64,
    pub total_par_dv01: Vec<f64>,
    pub total_eq_delta: f64,
    pub curve_recomputes: usize,
    pub trade_recomputes: usize,
    /// Per-call counters.  The legacy fields above remain cumulative for
    /// dashboard/backward compatibility; new telemetry should consume these
    /// explicit last-call values (or the corresponding DAG fields).
    pub last_curve_recomputes: usize,
    pub last_trade_recomputes: usize,
    pub total_curve_recomputes: usize,
    pub total_trade_recomputes: usize,
}

/// Bootstrap the OIS and forwarding curves used by the portfolio dependency
/// graph.  The three pillars are deliberately the contractual 1Y/2Y/5Y
/// maturities represented by the snapshot quote arrays (rather than array
/// indices interpreted as years).
pub fn build_curve_system(s: &MarketSnapshot) -> Arc<ql_curves::CalibratedCurveSystem> {
    try_build_curve_system(s).expect("curve calibration")
}

/// Fallible form of [`build_curve_system`].  Portfolio services should use
/// this form so malformed quotes are reported to the caller instead of being
/// converted into a panic (or, worse, an unconverged valuation).
pub fn try_build_curve_system(
    s: &MarketSnapshot,
) -> Result<Arc<ql_curves::CalibratedCurveSystem>, String> {
    let reference_date = s.ref_date;
    let mut ois = Vec::with_capacity(3);
    let mut fwd = Vec::with_capacity(3);
    for (years, quote) in [1_i32, 2, 5].into_iter().zip(s.ois_quotes) {
        let end = reference_date.add_days(365 * years);
        ois.push(ql_curves::OisSwapInstrument {
            start: reference_date,
            end,
            fixed_rate: quote,
            frequency: Frequency::Annual,
            day_counter: DayCounter::Act365Fixed,
        });
    }
    for (years, quote) in [1_i32, 2, 5].into_iter().zip(s.fwd_quotes) {
        let end = reference_date.add_days(365 * years);
        fwd.push(ql_curves::IrsSwapInstrument {
            start: reference_date,
            end,
            fixed_rate: quote,
            frequency: Frequency::Annual,
            day_counter: DayCounter::Act365Fixed,
        });
    }
    let ois_refs: Vec<&dyn ql_curves::SingleCurveHelper> = ois
        .iter()
        .map(|helper| helper as &dyn ql_curves::SingleCurveHelper)
        .collect();
    ql_curves::CalibratedCurveSystem::from_multicurve_helpers(reference_date, &ois_refs, &fwd)
        .map(Arc::new)
        .map_err(|e| format!("curve calibration failed: {e:?}"))
}

#[derive(Default)]
pub struct PortfolioDag {
    cached_curves: Option<(u64, Arc<ql_curves::CalibratedCurveSystem>)>,
    curve_recomputes: usize,
    cache: HashMap<usize, (u64, TradeValuation)>,
    trade_recomputes: usize,
    /// Lifetime number of curve recalibrations.  Unlike `curve_recomputes`,
    /// this value is never reset and is suitable for service telemetry.
    pub total_curve_recomputes: usize,
    /// Lifetime number of dirty trade valuations.
    pub total_trade_recomputes: usize,
}
impl PortfolioDag {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn evaluate_portfolio(
        &mut self,
        s: &MarketSnapshot,
        trades: &[PortfolioTrade],
    ) -> PortfolioRiskSummary {
        // Keep the original infallible API for callers that have already
        // validated their market snapshot.  New services should call
        // `try_evaluate_portfolio`, which returns a typed error instead.
        let mut report = self
            .try_evaluate_portfolio(trades, s)
            .unwrap_or_else(|e| panic!("portfolio evaluation failed: {e}"));
        // Both APIs expose per-call counters.  Lifetime totals are available
        // through the explicitly named total_* fields.
        report.curve_recomputes = report.last_curve_recomputes;
        report.trade_recomputes = report.last_trade_recomputes;
        report.total_curve_recomputes = self.total_curve_recomputes;
        report.total_trade_recomputes = self.total_trade_recomputes;
        report
    }

    /// Evaluate a portfolio and report duplicate IDs or curve calibration
    /// failures without panicking.  The argument order follows the canonical
    /// workflow (`trades`, then immutable market snapshot); the legacy
    /// `evaluate_portfolio(market, trades)` wrapper remains available.
    pub fn try_evaluate_portfolio(
        &mut self,
        trades: &[PortfolioTrade],
        s: &MarketSnapshot,
    ) -> Result<PortfolioRiskSummary, String> {
        let mut ids = std::collections::HashSet::with_capacity(trades.len());
        for trade in trades {
            let id = match trade {
                PortfolioTrade::IrsSwap { id, .. }
                | PortfolioTrade::HestonEuropeanOption { id, .. }
                | PortfolioTrade::ScriptedNote { id, .. } => *id,
            };
            if !ids.insert(id) {
                return Err(format!(
                    "portfolio trade ids must be unique (duplicate id {id})"
                ));
            }
        }
        // Per-call counters are deliberately separate from the cumulative
        // counters exposed above.  This makes cache hit/miss telemetry
        // unambiguous even when a DAG is reused for thousands of ticks.
        self.trade_recomputes = 0;
        self.curve_recomputes = 0;
        let ch = s.curve_hash();
        if self.cached_curves.as_ref().map(|entry| entry.0) != Some(ch) {
            let curves = try_build_curve_system(s)?;
            self.cached_curves = Some((ch, curves));
            self.curve_recomputes += 1;
            self.total_curve_recomputes += 1;
        }
        let curves = self
            .cached_curves
            .as_ref()
            .ok_or_else(|| "curve cache unavailable after calibration".to_string())?
            .1
            .clone();
        let dirty: Vec<_> = trades
            .iter()
            .filter_map(|t| {
                let id = match t {
                    PortfolioTrade::IrsSwap { id, .. }
                    | PortfolioTrade::HestonEuropeanOption { id, .. }
                    | PortfolioTrade::ScriptedNote { id, .. } => *id,
                };
                let h = trade_market_hash(t, s);
                (self.cache.get(&id).map(|x| x.0) != Some(h)).then_some((t, h))
            })
            .collect();
        self.trade_recomputes = dirty.len();
        self.total_trade_recomputes += dirty.len();
        // Evict removed/matured trades before inserting newly valued entries.
        // This is essential for long-lived services where portfolios are
        // replaced incrementally and IDs are not reused.
        self.cache.retain(|id, _| ids.contains(id));
        let vals: Vec<_> = dirty
            .par_iter()
            .map(|(t, h)| (*h, self.value_trade(t, s, &curves)))
            .collect();
        for (h, v) in vals {
            self.cache.insert(v.id, (h, v));
        }
        let mut out = PortfolioRiskSummary {
            // The fallible API reports per-call counters.  The compatibility
            // wrapper above converts its legacy fields to cumulative values.
            curve_recomputes: self.curve_recomputes,
            trade_recomputes: self.trade_recomputes,
            last_curve_recomputes: self.curve_recomputes,
            last_trade_recomputes: self.trade_recomputes,
            total_curve_recomputes: self.total_curve_recomputes,
            total_trade_recomputes: self.total_trade_recomputes,
            ..Default::default()
        };
        for t in trades {
            let id = match t {
                PortfolioTrade::IrsSwap { id, .. }
                | PortfolioTrade::HestonEuropeanOption { id, .. }
                | PortfolioTrade::ScriptedNote { id, .. } => *id,
            };
            if let Some((_, v)) = self.cache.get(&id) {
                out.total_npv += v.npv;
                out.total_eq_delta += v.eq_delta;
                if out.total_par_dv01.len() < v.par_dv01.len() {
                    out.total_par_dv01.resize(v.par_dv01.len(), 0.);
                }
                for (i, x) in v.par_dv01.iter().enumerate() {
                    out.total_par_dv01[i] += x;
                }
            }
        }
        Ok(out)
    }

    /// Convenience alias for integrations that mirror the legacy
    /// `evaluate_portfolio(market, trades)` argument order.
    pub fn try_evaluate_portfolio_market_first(
        &mut self,
        s: &MarketSnapshot,
        trades: &[PortfolioTrade],
    ) -> Result<PortfolioRiskSummary, String> {
        self.try_evaluate_portfolio(trades, s)
    }

    /// Remove all cached trade valuations while retaining calibrated curves.
    pub fn clear_cache(&mut self) {
        self.cache.clear();
    }

    /// Number of currently cached (live) trade valuations.
    pub fn cached_trade_count(&self) -> usize {
        self.cache.len()
    }
    fn value_trade(
        &self,
        t: &PortfolioTrade,
        s: &MarketSnapshot,
        curves: &ql_curves::CalibratedCurveSystem,
    ) -> TradeValuation {
        match t {
            PortfolioTrade::IrsSwap {
                id,
                notional,
                fixed_rate,
                maturity_years,
            } => {
                let n = *maturity_years;
                let (npv, dv) = curves.value_and_par_dv01(|params| {
                    let ois_n = curves.ois_pillars;
                    let ois = ql_curves::DiscountCurve::new(
                        s.ref_date,
                        curves.curves[0].times.clone(),
                        params[..ois_n].to_vec(),
                    );
                    let fwd = ql_curves::DiscountCurve::new(
                        s.ref_date,
                        curves.curves[1].times.clone(),
                        params[ois_n..].to_vec(),
                    );
                    let spec = VanillaSwapSpec {
                        effective_date: s.ref_date,
                        termination_date: s
                            .ref_date
                            .add_months((n as i32).saturating_mul(12), false),
                        notional: *notional,
                        fixed_rate_or_coupon: *fixed_rate,
                        fixed_tenor: ql_time::Period::Years(1),
                        float_tenor: ql_time::Period::Years(1),
                        calendar: Arc::new(ql_time::BitsetCalendar::weekend_only()),
                        bdc: ql_time::BusinessDayConvention::Unadjusted,
                        fixed_dc: DayCounter::Act365Fixed,
                        float_dc: DayCounter::Act365Fixed,
                        swap_type: SwapType::Payer,
                    };
                    let resolved: ResolvedVanillaSwap = spec.resolve(s.ref_date);
                    resolved.npv(&ois, &fwd)
                });
                TradeValuation {
                    id: *id,
                    npv,
                    par_dv01: dv,
                    eq_delta: 0.,
                }
            }
            PortfolioTrade::HestonEuropeanOption {
                id,
                strike,
                expiry,
                is_call,
                notional,
            } => {
                let (npv, par_dv01) = curves.value_and_par_dv01(|params| {
                    let ois = ql_curves::DiscountCurve::new(
                        s.ref_date,
                        curves.curves[0].times.clone(),
                        params[..3].to_vec(),
                    );
                    let r = -ois.discount_factor(*expiry).ln() / expiry.max(1e-6);
                    let hp = HestonParams {
                        v0: ql_ad::TapeVar::from_f64(s.heston_params.v0),
                        kappa: ql_ad::TapeVar::from_f64(s.heston_params.kappa),
                        theta: ql_ad::TapeVar::from_f64(s.heston_params.theta),
                        sigma_v: ql_ad::TapeVar::from_f64(s.heston_params.sigma_v),
                        rho: ql_ad::TapeVar::from_f64(s.heston_params.rho),
                    };
                    heston_cos_price(
                        ql_ad::TapeVar::from_f64(s.eq_spot),
                        ql_ad::TapeVar::from_f64(*strike),
                        ql_ad::TapeVar::from_f64(*expiry),
                        r,
                        ql_ad::TapeVar::zero(),
                        &hp,
                        128,
                        10.0,
                        *is_call,
                    ) * *notional
                });
                let spot = ql_ad::Dual::<1>::variable(s.eq_spot, 0);
                let params = HestonParams {
                    v0: ql_ad::Dual::from_value(s.heston_params.v0),
                    kappa: ql_ad::Dual::from_value(s.heston_params.kappa),
                    theta: ql_ad::Dual::from_value(s.heston_params.theta),
                    sigma_v: ql_ad::Dual::from_value(s.heston_params.sigma_v),
                    rho: ql_ad::Dual::from_value(s.heston_params.rho),
                };
                let r = -curves.curves[0].discount_factor(*expiry).ln() / expiry.max(1e-6);
                let price = heston_cos_price(
                    spot,
                    ql_ad::Dual::from_value(*strike),
                    ql_ad::Dual::from_value(*expiry),
                    ql_ad::Dual::from_value(r),
                    ql_ad::Dual::zero(),
                    &params,
                    128,
                    10.0,
                    *is_call,
                );
                TradeValuation {
                    id: *id,
                    npv,
                    par_dv01,
                    eq_delta: *notional * price.d[0],
                }
            }
            PortfolioTrade::ScriptedNote {
                id,
                program,
                expiry,
            } => {
                let (npv, par_dv01) = curves.value_and_par_dv01(|params| {
                    let ois = ql_curves::DiscountCurve::new(
                        s.ref_date,
                        curves.curves[0].times.clone(),
                        params[..3].to_vec(),
                    );
                    let r = -ois.discount_factor(*expiry).ln() / expiry.max(1e-6);
                    ql_scripting::ScriptMcEngine {
                        n_paths: 1024,
                        seed: 42,
                    }
                    .price_gbm(
                        program,
                        ql_ad::TapeVar::from_f64(s.eq_spot),
                        ql_ad::TapeVar::from_f64(s.eq_vol),
                        r,
                        ql_ad::TapeVar::zero(),
                    )
                });
                let r = -curves.curves[0].discount_factor(*expiry).ln() / expiry.max(1e-6);
                let spot = ql_ad::Dual::<1>::variable(s.eq_spot, 0);
                let price = ql_scripting::ScriptMcEngine {
                    n_paths: 1024,
                    seed: 42,
                }
                .price_gbm(
                    program,
                    spot,
                    ql_ad::Dual::from_value(s.eq_vol),
                    ql_ad::Dual::from_value(r),
                    ql_ad::Dual::zero(),
                );
                TradeValuation {
                    id: *id,
                    npv,
                    par_dv01,
                    eq_delta: price.d[0],
                }
            }
        }
    }
}
