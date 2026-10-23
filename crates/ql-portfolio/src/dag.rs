use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_models::{HestonParams, heston_cos_price};
use ql_scripting::{BytecodeProgram, OpCode};
use ql_time::{Date, DayCounter, Frequency};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;
pub fn hash_f64_slice(v: &[f64]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &x in v {
        h ^= if x == 0. { 0 } else { x.to_bits() };
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
            OpCode::LoadConst(x) | OpCode::SmoothStep(x) => h = hash_word(h, x.to_bits()),
            OpCode::Pay { discount_time } => h = hash_word(h, discount_time.to_bits()),
            OpCode::Autocall {
                barrier,
                half_width,
                discount_time,
            } => {
                h = hash_word(h, barrier.to_bits());
                h = hash_word(h, half_width.to_bits());
                h = hash_word(h, discount_time.to_bits());
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
        } => hash_word(hash_program(program), expiry.to_bits()),
    };
    match trade {
        PortfolioTrade::IrsSwap { .. } => hash_word(curve_hash, h),
        PortfolioTrade::HestonEuropeanOption { .. } | PortfolioTrade::ScriptedNote { .. } => {
            hash_word(curve_hash ^ equity_hash.rotate_left(13), h)
        }
    }
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
}

/// Bootstrap the OIS and forwarding curves used by the portfolio dependency
/// graph.  The three pillars are deliberately the contractual 1Y/2Y/5Y
/// maturities represented by the snapshot quote arrays (rather than array
/// indices interpreted as years).
pub fn build_curve_system(s: &MarketSnapshot) -> Arc<ql_curves::CalibratedCurveSystem> {
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
    Arc::new(
        ql_curves::CalibratedCurveSystem::from_multicurve_helpers(reference_date, &ois_refs, &fwd)
            .expect("curve calibration"),
    )
}

#[derive(Default)]
pub struct PortfolioDag {
    cached_curves: Option<(u64, Arc<ql_curves::CalibratedCurveSystem>)>,
    curve_recomputes: usize,
    cache: HashMap<usize, (u64, TradeValuation)>,
    trade_recomputes: usize,
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
        let mut ids = std::collections::HashSet::with_capacity(trades.len());
        assert!(
            trades.iter().all(|trade| {
                let id = match trade {
                    PortfolioTrade::IrsSwap { id, .. }
                    | PortfolioTrade::HestonEuropeanOption { id, .. }
                    | PortfolioTrade::ScriptedNote { id, .. } => *id,
                };
                ids.insert(id)
            }),
            "portfolio trade ids must be unique"
        );
        let ch = s.curve_hash();
        if self.cached_curves.as_ref().map(|entry| entry.0) != Some(ch) {
            let curves = build_curve_system(s);
            self.cached_curves = Some((ch, curves));
            self.curve_recomputes += 1;
        }
        let curves = self.cached_curves.as_ref().expect("curve cache").1.clone();
        let eh = s.equity_hash();
        let dirty: Vec<_> = trades
            .iter()
            .filter_map(|t| {
                let id = match t {
                    PortfolioTrade::IrsSwap { id, .. }
                    | PortfolioTrade::HestonEuropeanOption { id, .. }
                    | PortfolioTrade::ScriptedNote { id, .. } => *id,
                };
                let h = trade_dependency_hash(t, ch, eh);
                (self.cache.get(&id).map(|x| x.0) != Some(h)).then_some((t, h))
            })
            .collect();
        let vals: Vec<_> = dirty
            .par_iter()
            .map(|(t, h)| (*h, self.value_trade(t, s, &curves)))
            .collect();
        for (h, v) in vals {
            self.cache.insert(v.id, (h, v));
            self.trade_recomputes += 1;
        }
        let mut out = PortfolioRiskSummary {
            curve_recomputes: self.curve_recomputes,
            trade_recomputes: self.trade_recomputes,
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
        out
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
                    let ois = ql_curves::DiscountCurve::new(
                        s.ref_date,
                        curves.curves[0].times.clone(),
                        params[..3].to_vec(),
                    );
                    let fwd = ql_curves::DiscountCurve::new(
                        s.ref_date,
                        curves.curves[1].times.clone(),
                        params[3..].to_vec(),
                    );
                    let mut v = ql_ad::TapeVar::from_f64(0.);
                    for m in 1..=n {
                        let t = m as f64;
                        let fr = fwd.discount_factor(t - 1.) / fwd.discount_factor(t)
                            - ql_ad::TapeVar::from_f64(1.);
                        v += (fr - ql_ad::TapeVar::from_f64(*fixed_rate))
                            * ois.discount_factor(t)
                            * ql_ad::TapeVar::from_f64(*notional);
                    }
                    v
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
