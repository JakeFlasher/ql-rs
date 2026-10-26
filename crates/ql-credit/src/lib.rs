#![allow(clippy::needless_range_loop)]
pub mod phase34;
pub use phase34::*;
mod phase33;
pub use phase33::*;
pub mod phase32;
use nalgebra::{DMatrix, DVector};
pub use phase32::*;
use ql_ad::{Real, TapeScope, TapeVar};
use ql_curves::YieldTermStructure;
use ql_math::brent_solve;
use ql_time::Date;
use thiserror::Error;

mod phase30;
#[allow(unused_imports)]
pub use phase30::*;
mod phase31;
pub use phase31::*;

pub mod cdo;
pub mod cds_option;
pub use cdo::*;
pub use cds_option::*;

#[derive(Clone, Debug)]
pub struct JarrowLandoTurnbullRatingMigrationSpec<T: Real = f64> {
    pub initial_rating_index: usize,
    pub physical_generator_4x4: [[T; 4]; 4],
    pub risk_premium_multipliers: [T; 3],
    pub recovery_rate: T,
    pub risk_free_rate: T,
    pub maturity: T,
    pub base_coupon_rate: T,
    pub hy_downgrade_stepup_spread: T,
    pub counterparty_epe_notional: T,
}
#[derive(Clone, Debug)]
pub struct JarrowLandoTurnbullRatingMigrationReport<T: Real = f64> {
    pub transition_matrix_at_maturity: [[T; 4]; 4],
    pub max_row_sum_stochastic_error: T,
    pub cumulative_default_probabilities_by_rating: [T; 3],
    pub credit_spreads_bps_by_rating: [T; 3],
    pub downgrade_probability_at_maturity: T,
    pub vanilla_risky_bond_price: T,
    pub stepup_rating_sensitive_bond_price: T,
    pub stepup_clause_option_value: T,
    pub default_only_cva: T,
    pub rating_migration_adjusted_crcva: T,
}
fn mm4<T: Real>(a: &[[T; 4]; 4], b: &[[T; 4]; 4]) -> [[T; 4]; 4] {
    let mut c = [[T::zero(); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
fn solve4m<T: Real>(mut a: [[T; 4]; 4], mut b: [[T; 4]; 4]) -> [[T; 4]; 4] {
    for i in 0..4 {
        let mut p = i;
        for r in i + 1..4 {
            if a[r][i].abs().val() > a[p][i].abs().val() {
                p = r;
            }
        }
        a.swap(i, p);
        b.swap(i, p);
        let d = if a[i][i].abs().val() < 1e-14 {
            if a[i][i].val() < 0.0 {
                T::from_f64(-1e-14)
            } else {
                T::from_f64(1e-14)
            }
        } else {
            a[i][i]
        };
        for j in i..4 {
            a[i][j] /= d;
        }
        for j in 0..4 {
            b[i][j] /= d;
        }
        for r in 0..4 {
            if r != i {
                let f = a[r][i];
                for j in i..4 {
                    a[r][j] -= f * a[i][j];
                }
                for j in 0..4 {
                    b[r][j] -= f * b[i][j];
                }
            }
        }
    }
    b
}
pub fn evaluate_jlt_rating_migration_and_crcva<T: Real>(
    s: &JarrowLandoTurnbullRatingMigrationSpec<T>,
) -> JarrowLandoTurnbullRatingMigrationReport<T> {
    let mut q = [[T::zero(); 4]; 4];
    for i in 0..3 {
        let mut sm = T::zero();
        for j in 0..4 {
            if j != i {
                q[i][j] = s.physical_generator_4x4[i][j] * s.risk_premium_multipliers[i.min(2)];
                sm += q[i][j];
            }
        }
        q[i][i] = -sm;
    }
    let mut a = [[T::zero(); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            a[i][j] = q[i][j] * s.maturity / T::from_f64(16.0);
        }
    }
    let id = {
        let mut z = [[T::zero(); 4]; 4];
        for i in 0..4 {
            z[i][i] = T::one();
        }
        z
    };
    let a2 = mm4(&a, &a);
    let a4 = mm4(&a2, &a2);
    let a6 = mm4(&a4, &a2);
    let c0 = T::one();
    let c1 = T::from_f64(0.5);
    let c2 = T::from_f64(5.0 / 44.0);
    let c3 = T::from_f64(1.0 / 66.0);
    let c4 = T::from_f64(1.0 / 792.0);
    let c5 = T::from_f64(1.0 / 15840.0);
    let c6 = T::from_f64(1.0 / 665280.0);
    let mut poly = [[T::zero(); 4]; 4];
    let mut v = [[T::zero(); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            poly[i][j] = c1 * id[i][j] + c3 * a2[i][j] + c5 * a4[i][j];
            v[i][j] = c0 * id[i][j] + c2 * a2[i][j] + c4 * a4[i][j] + c6 * a6[i][j];
        }
    }
    let u = mm4(&a, &poly);
    let mut vm = [[T::zero(); 4]; 4];
    let mut vp = [[T::zero(); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            vm[i][j] = v[i][j] - u[i][j];
            vp[i][j] = v[i][j] + u[i][j];
        }
    }
    let mut p = solve4m(vm, vp);
    for _ in 0..4 {
        p = mm4(&p, &p);
    }
    let mut err = T::zero();
    for i in 0..4 {
        let mut rs = T::zero();
        for j in 0..4 {
            rs += p[i][j];
        }
        err = err.max_val((rs - T::one()).abs());
    }
    let mut def = [T::zero(); 3];
    let mut spr = [T::zero(); 3];
    for i in 0..3 {
        def[i] = p[i][3].max_val(T::zero());
        spr[i] = -T::from_f64(1e4) / s.maturity.max_val(T::from_f64(1e-12))
            * (T::one() - (T::one() - s.recovery_rate) * def[i])
                .max_val(T::from_f64(1e-12))
                .ln();
    }
    let i0 = s.initial_rating_index.min(2);
    let mut down = T::zero();
    for j in i0 + 1..3 {
        down += p[i0][j];
    }
    let disc = (-s.risk_free_rate * s.maturity).exp();
    let van = disc
        * ((T::one() - def[i0]) * (T::one() + s.base_coupon_rate * s.maturity)
            + s.recovery_rate * def[i0]);
    let step = van + disc * p[i0][2] * s.hy_downgrade_stepup_spread * s.maturity;
    let cva = s.counterparty_epe_notional * (T::one() - s.recovery_rate) * def[i0] * disc;
    let mig = cva + s.counterparty_epe_notional * (T::one() - s.recovery_rate) * down * disc;
    JarrowLandoTurnbullRatingMigrationReport {
        transition_matrix_at_maturity: p,
        max_row_sum_stochastic_error: err,
        cumulative_default_probabilities_by_rating: def,
        credit_spreads_bps_by_rating: spr,
        downgrade_probability_at_maturity: down,
        vanilla_risky_bond_price: van,
        stepup_rating_sensitive_bond_price: step,
        stepup_clause_option_value: (step - van).max_val(T::zero()),
        default_only_cva: cva,
        rating_migration_adjusted_crcva: mig,
    }
}

#[derive(Debug, Error, Clone)]
pub enum CreditError {
    #[error("invalid credit curve input")]
    InvalidInput,
    #[error("credit bootstrap solver failed: {0}")]
    Solver(String),
}

#[derive(Clone, Debug)]
pub struct SurvivalProbabilityCurve<T: Real = f64> {
    pub reference_date: Date,
    pub pillar_times: Vec<f64>,
    pub hazard_rates: Vec<T>,
}
impl<T: Real> SurvivalProbabilityCurve<T> {
    pub fn new(reference_date: Date, pillar_times: Vec<f64>, hazard_rates: Vec<T>) -> Self {
        assert_eq!(pillar_times.len(), hazard_rates.len());
        assert!(!pillar_times.is_empty());
        assert!(pillar_times.windows(2).all(|w| w[1] > w[0]));
        Self {
            reference_date,
            pillar_times,
            hazard_rates,
        }
    }
    pub fn hazard_rate(&self, t: f64) -> T {
        let i = self
            .pillar_times
            .partition_point(|x| *x < t)
            .min(self.hazard_rates.len() - 1);
        self.hazard_rates[i]
    }
    pub fn survival_probability(&self, t: f64) -> T {
        if t <= 0.0 {
            return T::one();
        }
        let mut acc = T::zero();
        let mut left = 0.0;
        for (i, &right) in self.pillar_times.iter().enumerate() {
            if t <= right {
                acc += self.hazard_rates[i] * (t - left);
                return (-acc).exp();
            }
            acc += self.hazard_rates[i] * (right - left);
            left = right;
        }
        acc += *self.hazard_rates.last().unwrap() * (t - left);
        (-acc).exp()
    }
    pub fn default_probability(&self, t1: f64, t2: f64) -> T {
        self.survival_probability(t1) - self.survival_probability(t2)
    }
    pub fn to_ad<S: Real>(&self, nodes: &[S]) -> SurvivalProbabilityCurve<S> {
        SurvivalProbabilityCurve::new(
            self.reference_date,
            self.pillar_times.clone(),
            nodes.to_vec(),
        )
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ProtectionSide {
    Buyer,
    Seller,
}
#[derive(Clone, Debug)]
pub struct CreditDefaultSwap {
    pub side: ProtectionSide,
    pub notional: f64,
    pub running_spread: f64,
    pub recovery_rate: f64,
    pub payment_times: Vec<f64>,
    pub accrual_fractions: Vec<f64>,
}
#[derive(Clone, Debug)]
pub struct CdsValuation<T: Real> {
    pub protection_leg_npv: T,
    pub premium_leg_npv: T,
    pub accrual_on_default_npv: T,
    pub risky_annuity: T,
    pub npv: T,
    pub fair_spread: T,
}

#[allow(clippy::too_many_arguments)]
fn interval_terms_from_values<T: Real>(
    p0: T,
    p1: T,
    q0: T,
    q1: T,
    hh: T,
    a: f64,
    b: f64,
    coupon_start: f64,
) -> (T, T) {
    let fh = (p0 / p1).ln();
    let x = fh + hh;
    let h = if x.val().abs() < 1e-4 {
        hh * p0 * q0 * (T::one() - x * 0.5 + x * x / 6.0 - x * x * x / 24.0 + x.powf(4.0) / 120.0)
    } else {
        hh / x * (p0 * q0 - p1 * q1)
    };
    let du = b - a;
    let tau0 = a - coupon_start;
    let aod = if x.val().abs() < 1e-4 {
        hh * p0
            * q0
            * (T::from_f64(tau0) * (T::one() - x * 0.5 + x * x / 6.0 - x * x * x / 24.0)
                + T::from_f64(du) * (T::from_f64(0.5) - x / 3.0 + x * x / 8.0 - x * x * x / 30.0))
    } else {
        hh / x
            * (T::from_f64(du) * ((p0 * q0 - p1 * q1) / x - p1 * q1)
                + T::from_f64(tau0) * (p0 * q0 - p1 * q1))
    };
    (h, aod)
}

pub struct IsdaCdsEngine;
impl IsdaCdsEngine {
    pub fn price<T: Real>(
        cds: &CreditDefaultSwap,
        discount: &impl YieldTermStructure<T>,
        credit: &SurvivalProbabilityCurve<T>,
    ) -> CdsValuation<T> {
        // Accumulate normalized legs first and apply notional exactly once at
        // the end.  This avoids magnifying cancellation round-off by a large
        // notional during every coupon operation.
        let mut protection_unit = T::zero();
        let mut aod_unit = T::zero();
        let mut premium_unit = T::zero();
        let mut prev = 0.0;
        let mut p0 = discount.discount_factor(0.0);
        let mut q0 = T::one();
        let mut pillar_idx = 0usize;
        let n_pillars = credit.pillar_times.len();
        for (j, &end) in cds.payment_times.iter().enumerate() {
            let alpha = cds.accrual_fractions.get(j).copied().unwrap_or(end - prev);
            let period_len = (end - prev).max(1e-12);
            let accrual_scale = alpha / period_len;
            while pillar_idx < n_pillars && credit.pillar_times[pillar_idx] <= prev {
                pillar_idx += 1;
            }
            let mut a = prev;
            while a < end {
                let segment = pillar_idx.min(n_pillars - 1);
                let next_pillar = credit.pillar_times.get(pillar_idx).copied();
                let b = next_pillar.filter(|&x| x < end).unwrap_or(end);
                let hh = credit.hazard_rates[segment] * (b - a);
                let q1 = q0 * (-hh).exp();
                let p1 = discount.discount_factor(b);
                let (pr, ao) = interval_terms_from_values(p0, p1, q0, q1, hh, a, b, prev);
                protection_unit += T::from_f64(1.0 - cds.recovery_rate) * pr;
                aod_unit += T::from_f64(cds.running_spread * accrual_scale) * ao;
                p0 = p1;
                q0 = q1;
                a = b;
                if next_pillar == Some(b) {
                    pillar_idx += 1;
                }
            }
            premium_unit += T::from_f64(cds.running_spread * alpha) * p0 * q0;
            prev = end;
        }
        let notional = T::from_f64(cds.notional);
        let protection = notional * protection_unit;
        let aod = notional * aod_unit;
        let premium = notional * premium_unit;
        let risky_annuity =
            (premium + aod) / T::from_f64(cds.notional * cds.running_spread.max(1e-16));
        let fair_spread = protection / T::from_f64(cds.notional) / risky_annuity;
        let net = protection - premium - aod;
        let npv = if cds.side == ProtectionSide::Buyer {
            net
        } else {
            -net
        };
        CdsValuation {
            protection_leg_npv: protection,
            premium_leg_npv: premium,
            accrual_on_default_npv: aod,
            risky_annuity,
            npv,
            fair_spread,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CdsQuote {
    pub maturity_years: f64,
    pub par_spread: f64,
}

fn benchmark_cds(quote: CdsQuote, recovery_rate: f64) -> CreditDefaultSwap {
    let n = (quote.maturity_years * 4.0).round() as usize;
    CreditDefaultSwap {
        side: ProtectionSide::Buyer,
        notional: 1.0,
        running_spread: quote.par_spread,
        recovery_rate,
        payment_times: (1..=n).map(|i| i as f64 / 4.0).collect(),
        accrual_fractions: vec![0.25; n],
    }
}

impl SurvivalProbabilityCurve<f64> {
    /// Sequentially bootstrap piecewise-constant hazard rates from par CDS
    /// quotes.  This deliberately only performs the one-dimensional solves;
    /// the calibrated wrapper below adds the cached IFT Jacobian.
    pub fn bootstrap_quotes(
        reference_date: Date,
        discount: ql_curves::DiscountCurve<f64>,
        recovery_rate: f64,
        quotes: &[CdsQuote],
    ) -> Result<Self, CreditError> {
        if quotes.is_empty()
            || quotes
                .iter()
                .any(|q| q.maturity_years <= 0.0 || q.par_spread <= 0.0)
            || quotes
                .windows(2)
                .any(|w| w[1].maturity_years <= w[0].maturity_years)
        {
            return Err(CreditError::InvalidInput);
        }
        let times: Vec<f64> = quotes.iter().map(|q| q.maturity_years).collect();
        let mut hazards = vec![0.01; times.len()];
        for k in 0..times.len() {
            let prior = hazards.clone();
            let quote = quotes[k];
            let f = |lambda: f64| {
                let mut trial = prior.clone();
                trial[k] = lambda.max(1e-12);
                let curve = Self::new(reference_date, times.clone(), trial);
                IsdaCdsEngine::price(&benchmark_cds(quote, recovery_rate), &discount, &curve).npv
            };
            let root = brent_solve(f, 1e-10, 5.0, 1e-15, 300)
                .map_err(|e| CreditError::Solver(e.to_string()))?;
            let mut best = root;
            let mut best_abs = f(root).abs();
            for bits in [root.to_bits().saturating_sub(1), root.to_bits() + 1] {
                let candidate = f64::from_bits(bits);
                if candidate > 1e-10 && candidate < 5.0 {
                    let abs = f(candidate).abs();
                    if abs < best_abs {
                        best = candidate;
                        best_abs = abs;
                    }
                }
            }
            hazards[k] = best;
        }
        Ok(Self::new(reference_date, times, hazards))
    }
}

#[derive(Clone, Debug)]
pub struct CalibratedCreditCurve {
    pub curve: SurvivalProbabilityCurve<f64>,
    pub quotes: Vec<CdsQuote>,
    pub sensitivity: DMatrix<f64>,
    pub discount: ql_curves::DiscountCurve<f64>,
    pub recovery_rate: f64,
}
impl CalibratedCreditCurve {
    pub fn bootstrap(
        reference_date: Date,
        discount: ql_curves::DiscountCurve<f64>,
        recovery_rate: f64,
        quotes: &[CdsQuote],
    ) -> Result<Self, CreditError> {
        let curve = SurvivalProbabilityCurve::bootstrap_quotes(
            reference_date,
            discount.clone(),
            recovery_rate,
            quotes,
        )?;
        let times = curve.pillar_times.clone();
        let hazards = curve.hazard_rates.clone();
        let m = quotes.len();
        let mut jl = DMatrix::zeros(m, m);
        let mut js = DMatrix::zeros(m, m);
        for j in 0..m {
            for i in 0..m {
                let hcurve = SurvivalProbabilityCurve::new(
                    reference_date,
                    times.clone(),
                    hazards
                        .iter()
                        .enumerate()
                        .map(|(k, &x)| {
                            if k == j {
                                ql_ad::Dual::<1>::variable(x, 0)
                            } else {
                                ql_ad::Dual::<1>::from_value(x)
                            }
                        })
                        .collect(),
                );
                let dcurve = discount.to_ad_curve(
                    &discount
                        .discounts
                        .iter()
                        .map(|&x| ql_ad::Dual::<1>::from_value(x))
                        .collect::<Vec<_>>(),
                );
                let v = IsdaCdsEngine::price(
                    &benchmark_cds(quotes[i], recovery_rate),
                    &dcurve,
                    &hcurve,
                );
                jl[(i, j)] = v.npv.d[0];
            }
        }
        for i in 0..m {
            let v =
                IsdaCdsEngine::price(&benchmark_cds(quotes[i], recovery_rate), &discount, &curve);
            js[(i, i)] = -v.risky_annuity.val();
        }
        let sens = jl
            .clone()
            .lu()
            .solve(&(-js))
            .ok_or_else(|| CreditError::Solver("singular credit IFT Jacobian".into()))?;
        Ok(Self {
            curve,
            quotes: quotes.to_vec(),
            sensitivity: sens,
            discount,
            recovery_rate,
        })
    }
    pub fn value_and_par_cs01<F>(&self, f: F) -> (f64, Vec<f64>)
    where
        F: FnOnce(&SurvivalProbabilityCurve<TapeVar<f64>>) -> TapeVar<f64>,
    {
        let guard = TapeScope::guard_with_capacity(16_384);
        let vars: Vec<_> = self
            .curve
            .hazard_rates
            .iter()
            .map(|&x| TapeScope::with_active(|t| t.input(x)))
            .collect();
        let c = self.curve.to_ad(&vars);
        let out = f(&c);
        let adj = TapeScope::gradient_from(out);
        let g: Vec<f64> = vars.iter().map(|v| adj.gradient(*v)).collect();
        let _ = guard.stop();
        let q = self.sensitivity.transpose() * DVector::from_vec(g);
        (out.val, q.iter().copied().collect())
    }

    /// Price a portfolio of CDSs and map its hazard-rate adjoint to par
    /// quote CS01s.  Unlike the fully generic closure API above, this path
    /// shares every survival-curve node across the portfolio.  It is still a
    /// single reverse tape pass, but avoids rebuilding the same piecewise
    /// hazard integration once for every trade.
    pub fn value_and_par_cs01_portfolio<D>(
        &self,
        portfolio: &[CreditDefaultSwap],
        discount: &D,
    ) -> (f64, Vec<f64>)
    where
        D: YieldTermStructure<f64>,
    {
        let mut times = vec![0.0];
        for cds in portfolio {
            let mut prev = 0.0;
            for &end in &cds.payment_times {
                times.push(end);
                times.extend(
                    self.curve
                        .pillar_times
                        .iter()
                        .copied()
                        .filter(|&x| x > prev && x < end),
                );
                prev = end;
            }
        }
        times.sort_by(f64::total_cmp);
        times.dedup_by(|a, b| (*a - *b).abs() < 1e-14);

        let guard = TapeScope::guard_with_capacity(16_384);
        let vars: Vec<_> = self
            .curve
            .hazard_rates
            .iter()
            .map(|&x| TapeScope::with_active(|t| t.input(x)))
            .collect();
        let mut q_cache = Vec::with_capacity(times.len());
        for &t in &times {
            let mut acc = TapeVar::from_f64(0.0);
            let mut left = 0.0;
            for (i, &right) in self.curve.pillar_times.iter().enumerate() {
                if t <= right {
                    acc += vars[i] * (t - left);
                    break;
                }
                acc += vars[i] * (right - left);
                left = right;
            }
            if t > *self.curve.pillar_times.last().unwrap() {
                let i = vars.len() - 1;
                acc += vars[i] * (t - left);
            }
            q_cache.push((-acc).exp());
        }
        let lookup = |t: f64| -> TapeVar<f64> {
            let i = times
                .binary_search_by(|x| x.total_cmp(&t))
                .expect("cached CDS time");
            q_cache[i]
        };

        // The common institutional CDS strip is a regular schedule starting
        // at time zero.  Build cumulative unit protection, accrual-on-default
        // and premium sums once for that grid; every trade then becomes a
        // handful of scalar tape operations instead of replaying all coupon
        // intervals.  Irregular schedules use the general path below.
        let regular_dt = portfolio.first().and_then(|cds| {
            let dt = cds.payment_times.first().copied()?;
            let regular = portfolio.iter().all(|trade| {
                trade.payment_times.iter().enumerate().all(|(i, &t)| {
                    (t - (i as f64 + 1.0) * dt).abs() < 1e-12
                        && (trade.accrual_fractions.get(i).copied().unwrap_or(dt) - dt).abs()
                            < 1e-12
                })
            });
            regular.then_some(dt)
        });
        if let Some(dt) = regular_dt {
            let max_n = portfolio
                .iter()
                .map(|trade| trade.payment_times.len())
                .max()
                .unwrap_or(0);
            let mut prot = vec![TapeVar::from_f64(0.0); max_n + 1];
            let mut aod = vec![TapeVar::from_f64(0.0); max_n + 1];
            let mut prem = vec![TapeVar::from_f64(0.0); max_n + 1];
            let mut q0 = lookup(0.0);
            let mut p0 = discount.discount_factor(0.0);
            let mut pillar_idx = 0usize;
            for n in 0..max_n {
                let prev = n as f64 * dt;
                let end = (n as f64 + 1.0) * dt;
                while pillar_idx < self.curve.pillar_times.len()
                    && self.curve.pillar_times[pillar_idx] <= prev
                {
                    pillar_idx += 1;
                }
                let mut a = prev;
                let mut protection = TapeVar::from_f64(0.0);
                let mut accrual = TapeVar::from_f64(0.0);
                while a < end {
                    let segment = pillar_idx.min(self.curve.hazard_rates.len() - 1);
                    let next = self.curve.pillar_times.get(pillar_idx).copied();
                    let b = next.filter(|&x| x < end).unwrap_or(end);
                    let q1 = lookup(b);
                    let p1 = discount.discount_factor(b);
                    let hh = vars[segment] * (b - a);
                    let (pr, ao) = interval_terms_from_values(
                        TapeVar::from_f64(p0),
                        TapeVar::from_f64(p1),
                        q0,
                        q1,
                        hh,
                        a,
                        b,
                        prev,
                    );
                    protection += pr;
                    accrual += ao;
                    p0 = p1;
                    q0 = q1;
                    a = b;
                    if next == Some(b) {
                        pillar_idx += 1;
                    }
                }
                let p_end = TapeVar::from_f64(p0);
                prem[n + 1] = prem[n] + dt * p_end * q0;
                prot[n + 1] = prot[n] + protection;
                aod[n + 1] = aod[n] + accrual;
            }
            let mut total = TapeVar::from_f64(0.0);
            for cds in portfolio {
                let n = cds.payment_times.len();
                let net = TapeVar::from_f64(cds.notional)
                    * (TapeVar::from_f64(1.0 - cds.recovery_rate) * prot[n]
                        - TapeVar::from_f64(cds.running_spread) * (prem[n] + aod[n]));
                total += if cds.side == ProtectionSide::Buyer {
                    net
                } else {
                    -net
                };
            }
            let adj = TapeScope::gradient_from(total);
            let gradient: Vec<f64> = vars.iter().map(|v| adj.gradient(*v)).collect();
            let _ = guard.stop();
            let q = self.sensitivity.transpose() * DVector::from_vec(gradient);
            return (total.val, q.iter().copied().collect());
        }

        let mut total = TapeVar::from_f64(0.0);
        for cds in portfolio {
            let mut protection_unit = TapeVar::from_f64(0.0);
            let mut aod_unit = TapeVar::from_f64(0.0);
            let mut premium_unit = TapeVar::from_f64(0.0);
            let mut prev = 0.0;
            let mut q0 = lookup(0.0);
            let mut p0 = discount.discount_factor(0.0);
            let mut pillar_idx = 0usize;
            for (j, &end) in cds.payment_times.iter().enumerate() {
                let alpha = cds.accrual_fractions.get(j).copied().unwrap_or(end - prev);
                let scale = alpha / (end - prev).max(1e-12);
                while pillar_idx < self.curve.pillar_times.len()
                    && self.curve.pillar_times[pillar_idx] <= prev
                {
                    pillar_idx += 1;
                }
                let mut a = prev;
                while a < end {
                    let segment = pillar_idx.min(self.curve.hazard_rates.len() - 1);
                    let next = self.curve.pillar_times.get(pillar_idx).copied();
                    let b = next.filter(|&x| x < end).unwrap_or(end);
                    let q1 = lookup(b);
                    let p1 = discount.discount_factor(b);
                    let hh = vars[segment] * (b - a);
                    let (pr, ao) = interval_terms_from_values(
                        TapeVar::from_f64(p0),
                        TapeVar::from_f64(p1),
                        q0,
                        q1,
                        hh,
                        a,
                        b,
                        prev,
                    );
                    protection_unit += (1.0 - cds.recovery_rate) * pr;
                    aod_unit += cds.running_spread * scale * ao;
                    p0 = p1;
                    q0 = q1;
                    a = b;
                    if next == Some(b) {
                        pillar_idx += 1;
                    }
                }
                premium_unit += cds.running_spread * alpha * p0 * q0;
                prev = end;
            }
            let net = TapeVar::from_f64(cds.notional) * (protection_unit - premium_unit - aod_unit);
            total += if cds.side == ProtectionSide::Buyer {
                net
            } else {
                -net
            };
        }
        let adj = TapeScope::gradient_from(total);
        let gradient: Vec<f64> = vars.iter().map(|v| adj.gradient(*v)).collect();
        let _ = guard.stop();
        let q = self.sensitivity.transpose() * DVector::from_vec(gradient);
        (total.val, q.iter().copied().collect())
    }
    pub fn jump_to_default(cds: &CreditDefaultSwap, current_npv: f64) -> f64 {
        let sign = if cds.side == ProtectionSide::Buyer {
            1.
        } else {
            -1.
        };
        sign * cds.notional * (1. - cds.recovery_rate) - current_npv
    }
}

#[derive(Copy, Clone, Debug)]
pub struct BlackCoxStructuralSpec<T: Real = f64> {
    pub asset_value: T,
    pub debt_face_value: T,
    pub default_barrier: T,
    pub barrier_growth_rate: T,
    pub asset_vol: T,
    pub risk_free_rate: T,
    pub asset_payout_rate: T,
    pub maturity: T,
    pub recovery_rate: T,
}

#[derive(Copy, Clone, Debug)]
pub struct StructuralCreditReport<T: Real> {
    pub equity_value: T,
    pub equity_vol: T,
    pub merton_distance_to_default: T,
    pub merton_survival_prob: T,
    pub merton_credit_spread: T,
    pub black_cox_survival_prob: T,
    pub black_cox_cds_spread: T,
}

pub fn analyze_structural_credit_black_cox<T: Real>(
    spec: &BlackCoxStructuralSpec<T>,
) -> StructuralCreditReport<T> {
    let v0 = spec.asset_value.max_val(T::from_f64(1e-8));
    let d_face = spec.debt_face_value.max_val(T::from_f64(1e-8));
    let t = spec.maturity.max_val(T::from_f64(1e-6));
    let sig = spec.asset_vol.max_val(T::from_f64(1e-6));
    let st = sig * t.sqrt();
    let r = spec.risk_free_rate;
    let q = spec.asset_payout_rate;

    let d1 = ((v0 / d_face).ln() + (r - q + sig * sig * 0.5) * t) / st;
    let d2 = d1 - st;
    let v_disc = v0 * (-q * t).exp();
    let d_disc = d_face * (-r * t).exp();
    let eq_val = (v_disc * d1.norm_cdf() - d_disc * d2.norm_cdf()).max_val(T::from_f64(1e-8));
    let eq_vol = sig * (v_disc * d1.norm_cdf()) / eq_val;
    let debt_val = (v_disc - eq_val).max_val(T::from_f64(1e-8));
    let merton_surv = d2.norm_cdf();
    let merton_spread = -(debt_val / d_disc).min_val(T::one()).ln() / t;

    let b0 = spec.default_barrier.max_val(T::from_f64(1e-8));
    let bc_surv = if b0.val() >= v0.val() {
        T::from_f64(1e-8)
    } else {
        let log_b = (b0 / v0).ln();
        let mu_tilde = r - q - spec.barrier_growth_rate - sig * sig * 0.5;
        let k_t = (d_face / (b0 * (spec.barrier_growth_rate * t).exp()))
            .ln()
            .max_val(T::zero());
        let z1 = (-log_b - k_t + mu_tilde * t) / st;
        let z2 = (log_b - k_t + mu_tilde * t) / st;
        let exp_arg = (mu_tilde * log_b * 2.0 / (sig * sig)).min_val(T::from_f64(50.0));
        (z1.norm_cdf() - exp_arg.exp() * z2.norm_cdf())
            .max_val(T::from_f64(1e-8))
            .min_val(T::one())
    };
    let lgd = (T::one() - spec.recovery_rate).max_val(T::zero());
    let bc_spread = -lgd * bc_surv.ln() / t;

    StructuralCreditReport {
        equity_value: eq_val,
        equity_vol: eq_vol,
        merton_distance_to_default: d2,
        merton_survival_prob: merton_surv,
        merton_credit_spread: merton_spread,
        black_cox_survival_prob: bc_surv,
        black_cox_cds_spread: bc_spread,
    }
}

#[derive(Clone, Debug)]
pub struct JltRatingMigrationSpec<T: Real = f64> {
    pub generator_matrix: Vec<Vec<T>>,
    pub risk_premium_multipliers: Vec<T>,
    pub initial_rating_index: usize,
    pub horizon_years: f64,
    pub payment_times: Vec<f64>,
    pub base_coupon_rate: T,
    pub downgrade_step_up_spreads: Vec<T>,
    pub recovery_rate: T,
    pub risk_free_rate: T,
}

#[derive(Clone, Debug)]
pub struct JltMigrationValuation<T: Real> {
    pub transition_probabilities: Vec<T>,
    pub survival_probability: T,
    pub default_probability: T,
    pub downgrade_probability: T,
    pub step_up_bond_price: T,
    pub step_up_coupon_pv_uplift: T,
    pub implied_credit_spread: T,
}

fn mat_mul_generic<T: Real>(a: &[Vec<T>], b: &[Vec<T>]) -> Vec<Vec<T>> {
    let m = a.len();
    let mut out = vec![vec![T::zero(); m]; m];
    for i in 0..m {
        for k in 0..m {
            let aik = a[i][k];
            for j in 0..m {
                out[i][j] += aik * b[k][j];
            }
        }
    }
    out
}

fn jlt_matrix_exp<T: Real>(q_tilde: &[Vec<T>], t: f64) -> Vec<Vec<T>> {
    let m = q_tilde.len();
    let scale = 16.0_f64;
    let dt = t.max(0.0) / scale;
    let mut a = vec![vec![T::zero(); m]; m];
    let mut sum = vec![vec![T::zero(); m]; m];
    for i in 0..m {
        sum[i][i] = T::one();
        for j in 0..m {
            a[i][j] = q_tilde[i][j] * dt;
        }
    }
    let mut term = sum.clone();
    for deg in 1..=8 {
        term = mat_mul_generic(&term, &a);
        let inv_k = 1.0 / deg as f64;
        for i in 0..m {
            for j in 0..m {
                term[i][j] *= T::from_f64(inv_k);
                sum[i][j] += term[i][j];
            }
        }
    }
    for _ in 0..4 {
        sum = mat_mul_generic(&sum, &sum);
    }
    sum
}

pub fn evaluate_jlt_rating_migration_credit<T: Real>(
    spec: &JltRatingMigrationSpec<T>,
) -> JltMigrationValuation<T> {
    let m = spec.generator_matrix.len();
    if m < 2
        || spec.generator_matrix.iter().any(|r| r.len() != m)
        || spec.risk_premium_multipliers.len() < m
        || spec.downgrade_step_up_spreads.len() < m
        || spec.initial_rating_index >= m - 1
    {
        return JltMigrationValuation {
            transition_probabilities: Vec::new(),
            survival_probability: T::one(),
            default_probability: T::zero(),
            downgrade_probability: T::zero(),
            step_up_bond_price: T::one(),
            step_up_coupon_pv_uplift: T::zero(),
            implied_credit_spread: T::zero(),
        };
    }
    let mut q_tilde = vec![vec![T::zero(); m]; m];
    for i in 0..(m - 1) {
        let pi_i = spec.risk_premium_multipliers[i].max_val(T::zero());
        let mut row_sum = T::zero();
        for j in 0..m {
            if i != j {
                let q_ij = spec.generator_matrix[i][j].max_val(T::zero()) * pi_i;
                q_tilde[i][j] = q_ij;
                row_sum += q_ij;
            }
        }
        q_tilde[i][i] = -row_sum;
    }

    let i0 = spec.initial_rating_index;
    let p_horizon = jlt_matrix_exp(&q_tilde, spec.horizon_years.max(1e-6));
    let row_h = p_horizon[i0].clone();
    let p_def_h = row_h[m - 1]
        .max_val(T::zero())
        .min_val(T::from_f64(0.999999));
    let p_surv_h = (T::one() - p_def_h).max_val(T::from_f64(1e-8));
    let mut p_down = T::zero();
    for p in row_h.iter().take(m - 1).skip(i0 + 1) {
        p_down += *p;
    }

    let mut base_pv = T::zero();
    let mut step_uplift = T::zero();
    let mut prev_t = 0.0;
    let mut prev_def = T::zero();
    for (idx, &t) in spec.payment_times.iter().enumerate() {
        let dt = (t - prev_t).max(1e-6);
        let p_t = jlt_matrix_exp(&q_tilde, t);
        let row_t = &p_t[i0];
        let def_t = row_t[m - 1];
        let surv_t = T::one() - def_t;
        let df = (-spec.risk_free_rate * t).exp();
        base_pv += spec.base_coupon_rate * dt * surv_t * df;
        for j in 0..(m - 1) {
            step_uplift += row_t[j] * spec.downgrade_step_up_spreads[j] * dt * df;
        }
        base_pv += spec.recovery_rate * (def_t - prev_def).max_val(T::zero()) * df;
        if idx + 1 == spec.payment_times.len() {
            base_pv += surv_t * df;
        }
        prev_t = t;
        prev_def = def_t;
    }
    let imp_spread =
        -p_surv_h.ln() / spec.horizon_years.max(1e-6) * (T::one() - spec.recovery_rate);

    JltMigrationValuation {
        transition_probabilities: row_h,
        survival_probability: p_surv_h,
        default_probability: p_def_h,
        downgrade_probability: p_down,
        step_up_bond_price: base_pv + step_uplift,
        step_up_coupon_pv_uplift: step_uplift,
        implied_credit_spread: imp_spread,
    }
}

#[derive(Clone, Debug)]
pub struct DuffieSingletonRecoveryAndCdsBasisSpec<T: Real = f64> {
    pub face_value: T,
    pub coupon_rate_annual: T,
    pub payment_freq_per_year: usize,
    pub maturity_years: T,
    pub risk_free_rate: T,
    pub hazard_rate: T,
    pub recovery_rate: T,
    pub market_bond_price: T,
    pub repo_funding_spread: T,
}
#[derive(Clone, Debug)]
pub struct DuffieSingletonRecoveryAndCdsBasisReport<T: Real = f64> {
    pub rmv_bond_price: T,
    pub rfv_bond_price: T,
    pub rt_bond_price: T,
    pub risk_free_bond_price: T,
    pub duffie_singleton_spread: T,
    pub par_cds_spread: T,
    pub bond_z_spread: T,
    pub raw_cds_bond_basis: T,
    pub funding_adjusted_net_basis: T,
    pub rmv_minus_rfv_premium_bias: T,
}
pub type DuffieSingletonLarixRecoveryAndBasisSpec<T = f64> =
    DuffieSingletonRecoveryAndCdsBasisSpec<T>;
pub type DuffieSingletonLarixRecoveryAndBasisReport<T = f64> =
    DuffieSingletonRecoveryAndCdsBasisReport<T>;
pub fn analyze_duffie_singleton_rmv_rfv_and_cds_basis<T: Real>(
    s: &DuffieSingletonRecoveryAndCdsBasisSpec<T>,
) -> DuffieSingletonRecoveryAndCdsBasisReport<T> {
    let f = s.payment_freq_per_year.max(1);
    let m = (s.maturity_years.val() * f as f64).round().max(1.0) as usize;
    let dt = s.maturity_years / T::from_f64(m as f64);
    let c = s.face_value * s.coupon_rate_annual * dt;
    let l = T::one() - s.recovery_rate;
    let ds = s.hazard_rate * l;
    let mut rmv = T::zero();
    let mut rfv = T::zero();
    let mut rt = T::zero();
    let mut ust = T::zero();
    let mut rec = T::zero();
    let mut ann = T::zero();
    for j in 1..=m {
        let t = dt * T::from_f64(j as f64);
        let t0 = dt * T::from_f64((j - 1) as f64);
        let d = (-s.risk_free_rate * t).exp();
        let surv = (-s.hazard_rate * t).exp();
        let surv0 = (-s.hazard_rate * t0).exp();
        let cdf = surv0 - surv;
        rmv += c * (-(s.risk_free_rate + ds) * t).exp();
        rfv += c * (-(s.risk_free_rate + s.hazard_rate) * t).exp();
        rfv += s.recovery_rate
            * s.face_value
            * (-s.risk_free_rate * (t - dt / T::from_f64(2.0))).exp()
            * cdf;
        let q = surv + s.recovery_rate * (T::one() - surv);
        rt += c * d * q;
        ust += c * d;
        rec += (-s.risk_free_rate * (t - dt / T::from_f64(2.0))).exp() * cdf;
        ann += dt * d * (surv + T::from_f64(0.5) * cdf);
    }
    let t = s.maturity_years;
    let surv = (-s.hazard_rate * t).exp();
    rmv += s.face_value * (-(s.risk_free_rate + ds) * t).exp();
    rfv += s.face_value * (-(s.risk_free_rate + s.hazard_rate) * t).exp();
    rt +=
        s.face_value * (-s.risk_free_rate * t).exp() * (surv + s.recovery_rate * (T::one() - surv));
    ust += s.face_value * (-s.risk_free_rate * t).exp();
    let mut z = ds;
    for _ in 0..12 {
        let mut pv = T::zero();
        let mut dp = T::zero();
        for j in 1..=m {
            let tt = dt * T::from_f64(j as f64);
            let cf = if j == m { c + s.face_value } else { c };
            let e = (-(s.risk_free_rate + z) * tt).exp();
            pv += cf * e;
            dp += -tt * cf * e;
        }
        let den = if dp.val().abs() < 1e-12 {
            T::from_f64(-1e-12)
        } else {
            dp
        };
        z -= (pv - s.market_bond_price) / den;
    }
    let cds = l * rec / ann.max_val(T::from_f64(1e-12));
    let basis = cds - z;
    DuffieSingletonRecoveryAndCdsBasisReport {
        rmv_bond_price: rmv,
        rfv_bond_price: rfv,
        rt_bond_price: rt,
        risk_free_bond_price: ust,
        duffie_singleton_spread: ds,
        par_cds_spread: cds,
        bond_z_spread: z,
        raw_cds_bond_basis: basis,
        funding_adjusted_net_basis: basis - s.repo_funding_spread,
        rmv_minus_rfv_premium_bias: rmv - rfv,
    }
}
/// Compatibility entry point for Duffie-Singleton RMV/RFV recovery and CDS-bond basis analysis.
pub fn evaluate_duffie_singleton_larix_recovery_and_basis<T: Real>(
    s: &DuffieSingletonRecoveryAndCdsBasisSpec<T>,
) -> DuffieSingletonRecoveryAndCdsBasisReport<T> {
    analyze_duffie_singleton_rmv_rfv_and_cds_basis(s)
}

#[derive(Clone, Debug)]
pub struct Ifrs9CeclMultiStageEclSpec<T: Real = f64> {
    pub annual_ttc_pds: [T; 3],
    pub annual_eads: [T; 3],
    pub base_lgd: T,
    pub lgd_macro_sensitivity: T,
    pub vasicek_asset_correlation: T,
    pub effective_interest_rate: T,
    pub scenario_weights: [T; 3],
    pub scenario_macro_shocks_z: [T; 3],
    pub origination_lifetime_pd: T,
    pub sicr_pd_multiplier_threshold: T,
}
#[derive(Clone, Debug)]
pub struct Ifrs9CeclMultiStageEclReport<T: Real = f64> {
    pub stage1_12m_ecl: T,
    pub stage2_lifetime_ecl: T,
    pub stage3_credit_impaired_ecl: T,
    pub current_weighted_lifetime_pd: T,
    pub sicr_stage2_allocation_weight: T,
    pub blended_ifrs9_ecl_provision: T,
    pub stage1_gross_interest_income: T,
    pub stage3_net_eir_interest_income: T,
}
fn norm_inv<T: Real>(p: T) -> T {
    // Acklam's rational approximation followed by Newton refinement in T.
    let pp = p
        .max_val(T::from_f64(1e-15))
        .min_val(T::one() - T::from_f64(1e-15));
    let a = [
        -39.6968302866538,
        220.946098424521,
        -275.928510446969,
        138.357751867269,
        -30.6647980661472,
        2.50662827745924,
    ];
    let b = [
        -54.4760987982241,
        161.585836858041,
        -155.698979859887,
        66.8013118877197,
        -13.2806815528857,
    ];
    let c = [
        -0.00778489400243029,
        -0.322396458041136,
        -2.40075827716184,
        -2.54973253934373,
        4.37466414146497,
        2.93816398269878,
    ];
    let d = [
        0.00778469570904146,
        0.32246712907004,
        2.445134137143,
        3.75440866190742,
    ];
    let mut x;
    if pp.val() < 0.02425 {
        let q = (-T::from_f64(2.0) * pp.ln()).sqrt();
        x = (((((T::from_f64(c[0]) * q + T::from_f64(c[1])) * q + T::from_f64(c[2])) * q
            + T::from_f64(c[3]))
            * q
            + T::from_f64(c[4]))
            * q
            + T::from_f64(c[5]))
            / ((((T::from_f64(d[0]) * q + T::from_f64(d[1])) * q + T::from_f64(d[2])) * q
                + T::from_f64(d[3]))
                * q
                + T::one());
    } else if pp.val() > 0.97575 {
        let q = (-T::from_f64(2.0) * (T::one() - pp).ln()).sqrt();
        x = -(((((T::from_f64(c[0]) * q + T::from_f64(c[1])) * q + T::from_f64(c[2])) * q
            + T::from_f64(c[3]))
            * q
            + T::from_f64(c[4]))
            * q
            + T::from_f64(c[5]))
            / ((((T::from_f64(d[0]) * q + T::from_f64(d[1])) * q + T::from_f64(d[2])) * q
                + T::from_f64(d[3]))
                * q
                + T::one());
    } else {
        let q = pp - T::from_f64(0.5);
        let r = q * q;
        x = (((((T::from_f64(a[0]) * r + T::from_f64(a[1])) * r + T::from_f64(a[2])) * r
            + T::from_f64(a[3]))
            * r
            + T::from_f64(a[4]))
            * r
            + T::from_f64(a[5]))
            * q
            / (((((T::from_f64(b[0]) * r + T::from_f64(b[1])) * r + T::from_f64(b[2])) * r
                + T::from_f64(b[3]))
                * r
                + T::from_f64(b[4]))
                * r
                + T::one());
    }
    for _ in 0..2 {
        x -= (x.norm_cdf() - pp) / x.norm_pdf().max_val(T::from_f64(1e-15));
    }
    x
}
pub fn compute_ifrs9_cecl_multistage_expected_credit_loss<T: Real>(
    s: &Ifrs9CeclMultiStageEclSpec<T>,
) -> Ifrs9CeclMultiStageEclReport<T> {
    let rho = s
        .vasicek_asset_correlation
        .max_val(T::from_f64(1e-8))
        .min_val(T::from_f64(0.999999));
    let sr = rho.sqrt();
    let id = (T::one() - rho).sqrt();
    let wsum = (s.scenario_weights[0] + s.scenario_weights[1] + s.scenario_weights[2])
        .max_val(T::from_f64(1e-14));
    let mut stage1 = T::zero();
    let mut stage2 = T::zero();
    let mut stage3 = T::zero();
    let mut life_pd = T::zero();
    for k in 0..3 {
        let w = s.scenario_weights[k] / wsum;
        let z = s.scenario_macro_shocks_z[k];
        let mut surv = T::one();
        let mut s1k = T::zero();
        let mut s2k = T::zero();
        for h in 0..3 {
            let pttc = s.annual_ttc_pds[h]
                .max_val(T::from_f64(1e-12))
                .min_val(T::one() - T::from_f64(1e-12));
            let pit = ((norm_inv(pttc) - sr * z) / id).norm_cdf();
            let d = surv * pit;
            surv *= T::one() - pit;
            life_pd += w * d;
            let lgd = (s.base_lgd + s.lgd_macro_sensitivity * (pit - pttc).max_val(T::zero()))
                .max_val(T::from_f64(0.01))
                .min_val(T::from_f64(0.99));
            let disc = (T::one() + s.effective_interest_rate).powf(-(h as f64 + 1.0));
            let loss = s.annual_eads[h] * lgd * d * disc;
            s2k += loss;
            if h == 0 {
                s1k = loss;
            }
        }
        stage1 += w * s1k;
        stage2 += w * s2k;
        let pttc = s.annual_ttc_pds[0]
            .max_val(T::from_f64(1e-12))
            .min_val(T::one() - T::from_f64(1e-12));
        let pit = ((norm_inv(pttc) - sr * z) / id).norm_cdf();
        let lgd = (s.base_lgd + s.lgd_macro_sensitivity * (pit - pttc).max_val(T::zero()))
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.99));
        stage3 += w * s.annual_eads[0] * lgd;
    }
    let mult = life_pd / (s.origination_lifetime_pd + T::from_f64(1e-14));
    let omega = T::one()
        / (T::one() + (-(mult - s.sicr_pd_multiplier_threshold) * T::from_f64(20.0)).exp());
    let blend = (T::one() - omega) * stage1 + omega * stage2;
    Ifrs9CeclMultiStageEclReport {
        stage1_12m_ecl: stage1,
        stage2_lifetime_ecl: stage2,
        stage3_credit_impaired_ecl: stage3,
        current_weighted_lifetime_pd: life_pd,
        sicr_stage2_allocation_weight: omega,
        blended_ifrs9_ecl_provision: blend,
        stage1_gross_interest_income: s.annual_eads[0] * s.effective_interest_rate,
        stage3_net_eir_interest_income: (s.annual_eads[0] - stage3).max_val(T::zero())
            * s.effective_interest_rate,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct LandoCoxAffineDefaultSpec<T: Real = f64> {
    pub risk_free_rate: T,
    pub base_default_hazard: T,
    pub initial_credit_intensity_y0: T,
    pub credit_kappa: T,
    pub credit_theta: T,
    pub credit_sigma: T,
    pub initial_liquidity_intensity_z0: T,
    pub liquidity_kappa: T,
    pub liquidity_theta: T,
    pub liquidity_sigma: T,
    pub recovery_rate: T,
    pub maturity: T,
}
#[derive(Clone, Copy, Debug)]
pub struct LandoCoxAffineDefaultReport<T: Real = f64> {
    pub survival_probability: T,
    pub liquidity_discount_factor: T,
    pub defaultable_zero_coupon_price: T,
    pub recovery_leg_pv: T,
    pub risky_corporate_bond_price: T,
    pub total_yield_spread_bps: T,
    pub pure_credit_spread_bps: T,
    pub pure_liquidity_spread_bps: T,
    pub credit_spread_fraction: T,
}
fn cir_surv<T: Real>(k: T, th: T, sig: T, y: T, t: T) -> T {
    if k.val().abs() + th.val().abs() + sig.val().abs() + y.val().abs() < 1e-14 {
        return T::one();
    }
    let e = (-k * t).exp();
    let g = (k * k + T::from_f64(2.0) * sig * sig).sqrt();
    let b = T::from_f64(2.0) * (e - T::one()) / ((k + g) * (e - T::one()) + T::from_f64(2.0) * g);
    let a = ((T::from_f64(2.0) * g * ((k + g) * t / T::from_f64(2.0)).exp())
        / ((k + g) * (e - T::one()) + T::from_f64(2.0) * g))
        .pow(T::from_f64(2.0) * k * th / (sig * sig + T::from_f64(1e-14)));
    a * (-b * y).exp()
}
pub fn evaluate_lando_cox_affine_default_and_liquidity_decomposition<T: Real>(
    s: &LandoCoxAffineDefaultSpec<T>,
) -> LandoCoxAffineDefaultReport<T> {
    let t = s.maturity;
    let eps = T::from_f64(1e-14);
    let surv = (-s.base_default_hazard * t).exp()
        * cir_surv(
            s.credit_kappa,
            s.credit_theta,
            s.credit_sigma,
            s.initial_credit_intensity_y0,
            t,
        );
    let liq = cir_surv(
        s.liquidity_kappa,
        s.liquidity_theta,
        s.liquidity_sigma,
        s.initial_liquidity_intensity_z0,
        t,
    );
    let df = (-s.risk_free_rate * t).exp();
    let defp = df * surv;
    const X: [f64; 16] = [
        -0.98940093499165,
        -0.944575023073233,
        -0.865631202387832,
        -0.755404408355003,
        -0.617876244402644,
        -0.458016777657227,
        -0.281603550779259,
        -0.095012509837637,
        0.095012509837637,
        0.281603550779259,
        0.458016777657227,
        0.617876244402644,
        0.755404408355003,
        0.865631202387832,
        0.944575023073233,
        0.98940093499165,
    ];
    const W: [f64; 16] = [
        0.027152459411754,
        0.062253523938648,
        0.095158511682493,
        0.124628971255534,
        0.149595988816577,
        0.169156519395003,
        0.182603415044924,
        0.189450610455069,
        0.189450610455069,
        0.182603415044924,
        0.169156519395003,
        0.149595988816577,
        0.124628971255534,
        0.095158511682493,
        0.062253523938648,
        0.027152459411754,
    ];
    let cir_b = |k: T, sig: T, u: T| -> T {
        let g = (k * k + T::from_f64(2.0) * sig * sig).sqrt();
        let e = (-k * u).exp();
        T::from_f64(2.0) * (e - T::one()) / ((k + g) * (e - T::one()) + T::from_f64(2.0) * g)
    };
    let mut rec = T::zero();
    for i in 0..16 {
        let u = t * (T::from_f64(X[i]) + T::one()) / T::from_f64(2.0);
        let b = cir_b(s.credit_kappa, s.credit_sigma, u);
        let g = (s.credit_kappa * s.credit_kappa
            + T::from_f64(2.0) * s.credit_sigma * s.credit_sigma)
            .sqrt();
        let e = (g * u).exp();
        let den = (s.credit_kappa + g) * (-s.credit_kappa * u).exp() - (s.credit_kappa + g)
            + T::from_f64(2.0) * g;
        let bp = T::from_f64(4.0) * g * g * e / (den * den).max_val(eps);
        let sd = (-s.base_default_hazard * u).exp()
            * cir_surv(
                s.credit_kappa,
                s.credit_theta,
                s.credit_sigma,
                s.initial_credit_intensity_y0,
                u,
            );
        let sl = cir_surv(
            s.liquidity_kappa,
            s.liquidity_theta,
            s.liquidity_sigma,
            s.initial_liquidity_intensity_z0,
            u,
        );
        let hz = s.base_default_hazard - s.credit_kappa * s.credit_theta * b
            + bp * s.initial_credit_intensity_y0;
        let q = s.recovery_rate * (-s.risk_free_rate * u).exp() * sd * sl * hz;
        rec += q * t * T::from_f64(0.5 * W[i]);
    }
    let credit = -(surv.max_val(eps)).ln() / t.max_val(eps);
    let lspread = -(liq.max_val(eps)).ln() / t.max_val(eps);
    let total = credit + lspread;
    let risky = df * (-total * t).exp() + rec;
    let frac = if total.val().abs() <= 1e-14 {
        T::one()
    } else {
        credit / total
    };
    LandoCoxAffineDefaultReport {
        survival_probability: surv,
        liquidity_discount_factor: liq,
        defaultable_zero_coupon_price: defp,
        recovery_leg_pv: rec,
        risky_corporate_bond_price: risky,
        total_yield_spread_bps: total * T::from_f64(1e4),
        pure_credit_spread_bps: credit * T::from_f64(1e4),
        pure_liquidity_spread_bps: lspread * T::from_f64(1e4),
        credit_spread_fraction: frac,
    }
}

#[derive(Clone, Debug)]
pub struct FedsDeSpiegeleerAt1CocoDualTriggerSpec<T: Real = f64> {
    pub notional: T,
    pub coupon_rate: T,
    pub maturity: T,
    pub risk_free_rate: T,
    pub initial_cet1_ratio: T,
    pub mechanical_cet1_trigger: T,
    pub ponv_cet1_trigger: T,
    pub supervisory_ponv_hazard_lambda: T,
    pub asset_volatility: T,
    pub cet1_asset_elasticity: T,
    pub initial_stock_price: T,
    pub conversion_price: T,
    pub floor_stock_price: T,
    pub existing_common_shares: T,
    pub write_down_fraction: T,
    pub is_equity_conversion: bool,
}

#[derive(Clone, Debug)]
pub struct FedsDeSpiegeleerAt1CocoDualTriggerReport<T: Real = f64> {
    pub mechanical_cet1_trigger_prob: T,
    pub total_dual_trigger_prob: T,
    pub trigger_implied_stock_price: T,
    pub conversion_share_count: T,
    pub equity_dilution_fraction: T,
    pub effective_trigger_recovery_rate: T,
    pub coupon_leg_pv: T,
    pub principal_survival_pv: T,
    pub trigger_recovery_pv: T,
    pub at1_coco_fair_value: T,
    pub at1_coco_credit_spread_bps: T,
    pub write_down_vs_conversion_pv_wedge: T,
}

pub fn evaluate_feds_despiegeleer_at1_coco_dual_trigger<T: Real>(
    s: &FedsDeSpiegeleerAt1CocoDualTriggerSpec<T>,
) -> FedsDeSpiegeleerAt1CocoDualTriggerReport<T> {
    let t = s.maturity.max_val(T::from_f64(1e-8));
    let sigma_a = s.asset_volatility.abs().max_val(T::from_f64(1e-8));
    let b_mech = (s.initial_cet1_ratio - s.mechanical_cet1_trigger).max_val(T::from_f64(1e-6))
        / s.cet1_asset_elasticity.abs().max_val(T::from_f64(1e-6));
    let mu_a = s.risk_free_rate - T::from_f64(0.5) * sigma_a * sigma_a;
    let st = sigma_a * t.sqrt();
    let mech = ((-b_mech - mu_a * t) / st).norm_cdf()
        + (-T::from_f64(2.0) * mu_a * b_mech / (sigma_a * sigma_a)).exp()
            * ((-b_mech + mu_a * t) / st).norm_cdf();
    let ponv_gap = (s.mechanical_cet1_trigger - s.ponv_cet1_trigger).abs();
    let sup_hazard = s.supervisory_ponv_hazard_lambda.max_val(T::zero()) * (T::one() + ponv_gap);
    let sup = T::one() - (-sup_hazard * t).exp();
    let total = (mech + (T::one() - mech) * sup).min_val(T::one());
    let trigger_stock = s.initial_stock_price * (-b_mech).exp();
    let eff_conv_price = trigger_stock
        .max_val(s.floor_stock_price)
        .min_val(s.conversion_price.max_val(s.floor_stock_price));
    let shares = s.notional / eff_conv_price.max_val(T::from_f64(1e-8));
    let dilution = shares / (s.existing_common_shares + shares).max_val(T::from_f64(1e-12));
    let diluted_stock = (s.existing_common_shares * trigger_stock + s.notional)
        / (s.existing_common_shares + shares).max_val(T::from_f64(1e-12));
    let recovery_conv = (shares * diluted_stock / s.notional.max_val(T::from_f64(1e-12)))
        .min_val(T::one() - T::from_f64(1e-10));
    let recovery_wd = (T::one() - s.write_down_fraction).max_val(T::zero());
    let recovery = if s.is_equity_conversion {
        recovery_conv
    } else {
        recovery_wd
    };
    let df = (-s.risk_free_rate * t).exp();
    let h_trig = -(T::one() - total).max_val(T::from_f64(1e-8)).ln() / t;
    let rt = s.risk_free_rate + h_trig;
    let ann_surv = (T::one() - (-rt * t).exp()) / rt.max_val(T::from_f64(1e-10));
    let coupon = s.notional * s.coupon_rate * ann_surv;
    let principal = s.notional * (T::one() - total) * df;
    let recpv = s.notional * total * recovery * df;
    let fair = coupon + principal + recpv;
    let spread = h_trig * (T::one() - recovery) * T::from_f64(1e4);
    let wedge = s.notional * total * df * (recovery_conv - recovery_wd);
    FedsDeSpiegeleerAt1CocoDualTriggerReport {
        mechanical_cet1_trigger_prob: mech,
        total_dual_trigger_prob: total,
        trigger_implied_stock_price: trigger_stock,
        conversion_share_count: shares,
        equity_dilution_fraction: dilution,
        effective_trigger_recovery_rate: recovery,
        coupon_leg_pv: coupon,
        principal_survival_pv: principal,
        trigger_recovery_pv: recpv,
        at1_coco_fair_value: fair,
        at1_coco_credit_spread_bps: spread,
        write_down_vs_conversion_pv_wedge: wedge,
    }
}

#[derive(Clone, Debug)]
pub struct LandoIsraelGeneratorLogarithmSpec<T: Real = f64> {
    pub transition_matrix_3x3: [[T; 3]; 3],
    pub series_order_k: usize,
    pub fractional_horizon_years: T,
    pub recovery_rate: T,
    pub risk_free_rate: T,
}
#[derive(Clone, Debug)]
pub struct LandoIsraelGeneratorLogarithmReport<T: Real = f64> {
    pub raw_generator_matrix: [[T; 3]; 3],
    pub regularized_generator_matrix: [[T; 3]; 3],
    pub generator_row_sum_max_error: T,
    pub off_diagonal_negativity_repaired: T,
    pub fractional_transition_matrix: [[T; 3]; 3],
    pub ig_to_default_fractional_pd: T,
    pub hy_to_default_fractional_pd: T,
    pub migration_adjusted_ig_spread_bps: T,
}
fn lando_mm3<T: Real>(a: &[[T; 3]; 3], b: &[[T; 3]; 3]) -> [[T; 3]; 3] {
    let mut c = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
fn lando_exp3<T: Real>(a: &[[T; 3]; 3], tau: T) -> [[T; 3]; 3] {
    let mut out = [[T::zero(); 3]; 3];
    for i in 0..3 {
        out[i][i] = T::one();
    }
    let mut term = out;
    let at = {
        let mut z = [[T::zero(); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                z[i][j] = a[i][j] * tau;
            }
        }
        z
    };
    let mut fact = 1.0;
    for n in 1..=12 {
        term = lando_mm3(&term, &at);
        fact *= n as f64;
        for i in 0..3 {
            for j in 0..3 {
                out[i][j] += term[i][j] / fact;
            }
        }
    }
    out
}
pub fn compute_lando_israel_generator_logarithm_and_migration<T: Real>(
    s: &LandoIsraelGeneratorLogarithmSpec<T>,
) -> LandoIsraelGeneratorLogarithmReport<T> {
    let mut x = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            x[i][j] = s.transition_matrix_3x3[i][j] - if i == j { T::one() } else { T::zero() };
        }
    }
    let mut raw = [[T::zero(); 3]; 3];
    let mut power = x;
    let order = s.series_order_k.max(1);
    for k in 1..=order {
        let sign = if k % 2 == 1 { 1.0 } else { -1.0 };
        for i in 0..3 {
            for j in 0..3 {
                raw[i][j] += power[i][j] * (sign / k as f64);
            }
        }
        power = lando_mm3(&power, &x);
    }
    let mut reg = [[T::zero(); 3]; 3];
    let mut repaired = T::zero();
    let smooth_max = |a: T, b: T, eps: f64| {
        (a + b + ((a - b) * (a - b) + T::from_f64(eps * eps)).sqrt()) / T::from_f64(2.0)
    };
    for i in 0..3 {
        let mut row = T::zero();
        for j in 0..3 {
            if i != j {
                let v = smooth_max(raw[i][j], T::zero(), 1e-6);
                repaired += smooth_max(T::zero() - raw[i][j], T::zero(), 1e-6);
                reg[i][j] = v;
                row += v;
            }
        }
        reg[i][i] = -row;
    }
    let mut err = T::zero();
    for i in 0..3 {
        let mut sum = T::zero();
        for j in 0..3 {
            sum += reg[i][j];
        }
        err = err.max_val(sum.abs());
    }
    let frac = lando_exp3(&reg, s.fractional_horizon_years);
    let igpd = frac[0][2].max_val(T::zero());
    let hypd = frac[1][2].max_val(T::zero());
    let tau = s.fractional_horizon_years.max_val(T::from_f64(1e-6));
    let cum_pd = frac[0][2]
        .min_val(T::from_f64(0.999999))
        .max_val(T::from_f64(1e-12));
    let lgd = T::one() - s.recovery_rate;
    let risky_annuity_adj = T::one() + s.risk_free_rate * tau * T::from_f64(0.5);
    let spread = -(T::one() - lgd * cum_pd).ln() / tau * risky_annuity_adj * T::from_f64(10000.0);
    LandoIsraelGeneratorLogarithmReport {
        raw_generator_matrix: raw,
        regularized_generator_matrix: reg,
        generator_row_sum_max_error: err,
        off_diagonal_negativity_repaired: repaired,
        fractional_transition_matrix: frac,
        ig_to_default_fractional_pd: igpd,
        hy_to_default_fractional_pd: hypd,
        migration_adjusted_ig_spread_bps: spread,
    }
}

#[cfg(test)]
mod phase23_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn jlt_rows_stochastic() {
        let s = JarrowLandoTurnbullRatingMigrationSpec {
            initial_rating_index: 0,
            physical_generator_4x4: [
                [-0.1, 0.08, 0.01, 0.01],
                [0.02, -0.08, 0.04, 0.02],
                [0.0, 0.02, -0.05, 0.03],
                [0.0; 4],
            ],
            risk_premium_multipliers: [1.0; 3],
            recovery_rate: 0.4,
            risk_free_rate: 0.02,
            maturity: 5.0,
            base_coupon_rate: 0.03,
            hy_downgrade_stepup_spread: 0.01,
            counterparty_epe_notional: 1.0,
        };
        let r = evaluate_jlt_rating_migration_and_crcva(&s);
        for row in r.transition_matrix_at_maturity {
            assert!((row.iter().sum::<f64>() - 1.0).abs() < 1e-8);
        }
    }

    #[test]
    fn dual_fd_parity() {
        let mk = |x: f64| {
            let mut s = JarrowLandoTurnbullRatingMigrationSpec {
                initial_rating_index: 0,
                physical_generator_4x4: [
                    [-0.1, 0.08, 0.01, 0.01],
                    [0.02, -0.08, 0.04, 0.02],
                    [0.0, 0.02, -0.05, 0.03],
                    [0.0; 4],
                ],
                risk_premium_multipliers: [1.0, 1.0, 1.0],
                recovery_rate: 0.4,
                risk_free_rate: 0.02,
                maturity: 5.0,
                base_coupon_rate: 0.03,
                hy_downgrade_stepup_spread: 0.01,
                counterparty_epe_notional: 1.0,
            };
            s.risk_premium_multipliers[0] = x;
            s
        };
        let h = 1e-5;
        let fd = (evaluate_jlt_rating_migration_and_crcva(&mk(1.0 + h))
            .credit_spreads_bps_by_rating[0]
            - evaluate_jlt_rating_migration_and_crcva(&mk(1.0 - h)).credit_spreads_bps_by_rating
                [0])
            / (2.0 * h);
        let ad = evaluate_jlt_rating_migration_and_crcva(&JarrowLandoTurnbullRatingMigrationSpec {
            initial_rating_index: 0,
            physical_generator_4x4: [
                [
                    -Dual::<1>::from_value(0.1),
                    Dual::<1>::from_value(0.08),
                    Dual::<1>::from_value(0.01),
                    Dual::<1>::from_value(0.01),
                ],
                [
                    Dual::<1>::from_value(0.02),
                    -Dual::<1>::from_value(0.08),
                    Dual::<1>::from_value(0.04),
                    Dual::<1>::from_value(0.02),
                ],
                [
                    Dual::<1>::zero(),
                    Dual::<1>::from_value(0.02),
                    -Dual::<1>::from_value(0.05),
                    Dual::<1>::from_value(0.03),
                ],
                [Dual::<1>::zero(); 4],
            ],
            risk_premium_multipliers: [
                Dual::<1>::variable(1.0, 0),
                Dual::<1>::one(),
                Dual::<1>::one(),
            ],
            recovery_rate: Dual::<1>::from_value(0.4),
            risk_free_rate: Dual::<1>::from_value(0.02),
            maturity: Dual::<1>::from_value(5.0),
            base_coupon_rate: Dual::<1>::from_value(0.03),
            hy_downgrade_stepup_spread: Dual::<1>::from_value(0.01),
            counterparty_epe_notional: Dual::<1>::one(),
        })
        .credit_spreads_bps_by_rating[0];
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase17_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn black_cox_first_passage_dominance_and_ad() {
        let s = BlackCoxStructuralSpec {
            asset_value: 140.0,
            debt_face_value: 100.0,
            default_barrier: 100.0 * (-0.02 * 3.0f64).exp(),
            barrier_growth_rate: 0.02,
            asset_vol: 0.25,
            risk_free_rate: 0.03,
            asset_payout_rate: 0.01,
            maturity: 3.0,
            recovery_rate: 0.40,
        };
        let rep = analyze_structural_credit_black_cox(&s);
        assert!(rep.black_cox_survival_prob <= rep.merton_survival_prob + 1e-12);
        assert!(rep.black_cox_cds_spread > 0.0);

        let h = 1e-5;
        let eval_v = |v: f64| {
            let mut sp = s;
            sp.asset_value = v;
            analyze_structural_credit_black_cox(&sp).black_cox_cds_spread
        };
        let fd = (eval_v(140.0 + h) - eval_v(140.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sd = BlackCoxStructuralSpec {
            asset_value: D::variable(140.0, 0),
            debt_face_value: D::from_value(100.0),
            default_barrier: D::from_value(100.0 * (-0.02 * 3.0f64).exp()),
            barrier_growth_rate: D::from_value(0.02),
            asset_vol: D::from_value(0.25),
            risk_free_rate: D::from_value(0.03),
            asset_payout_rate: D::from_value(0.01),
            maturity: D::from_value(3.0),
            recovery_rate: D::from_value(0.40),
        };
        let ad = analyze_structural_credit_black_cox(&sd).black_cox_cds_spread;
        assert!(ad.d[0] < 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase19_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn jlt_markov_generator_conservation_step_up_uplift_and_ad() {
        // 4 states: 0 = IG-A, 1 = IG-BBB, 2 = HY-BB, 3 = Default (absorbing)
        let q = vec![
            vec![-0.06, 0.045, 0.010, 0.005],
            vec![0.020, -0.09, 0.055, 0.015],
            vec![0.005, 0.025, -0.11, 0.080],
            vec![0.0, 0.0, 0.0, 0.0],
        ];
        let spec = JltRatingMigrationSpec {
            generator_matrix: q.clone(),
            risk_premium_multipliers: vec![1.2, 1.4, 1.6, 1.0],
            initial_rating_index: 0,
            horizon_years: 3.0,
            payment_times: vec![1.0, 2.0, 3.0],
            base_coupon_rate: 0.04,
            downgrade_step_up_spreads: vec![0.0, 0.005, 0.015, 0.0],
            recovery_rate: 0.40,
            risk_free_rate: 0.03,
        };
        let val = evaluate_jlt_rating_migration_credit(&spec);
        let row_sum: f64 = val.transition_probabilities.iter().sum();
        assert!((row_sum - 1.0).abs() < 1e-10);
        assert!(val.transition_probabilities.iter().all(|&p| p >= 0.0));
        assert!(val.step_up_coupon_pv_uplift > 0.0);

        let h = 1e-5;
        let eval_pi0 = |pi0: f64| {
            let mut s = spec.clone();
            s.risk_premium_multipliers[0] = pi0;
            evaluate_jlt_rating_migration_credit(&s).implied_credit_spread
        };
        let fd = (eval_pi0(1.2 + h) - eval_pi0(1.2 - h)) / (2.0 * h);
        type D = Dual<1>;
        let q_d: Vec<Vec<D>> = q
            .iter()
            .map(|r| r.iter().map(|&x| D::from_value(x)).collect())
            .collect();
        let spec_d = JltRatingMigrationSpec {
            generator_matrix: q_d,
            risk_premium_multipliers: vec![
                D::variable(1.2, 0),
                D::from_value(1.4),
                D::from_value(1.6),
                D::one(),
            ],
            initial_rating_index: 0,
            horizon_years: 3.0,
            payment_times: vec![1.0, 2.0, 3.0],
            base_coupon_rate: D::from_value(0.04),
            downgrade_step_up_spreads: vec![
                D::zero(),
                D::from_value(0.005),
                D::from_value(0.015),
                D::zero(),
            ],
            recovery_rate: D::from_value(0.40),
            risk_free_rate: D::from_value(0.03),
        };
        let ad = evaluate_jlt_rating_migration_credit(&spec_d).implied_credit_spread;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase24_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(h: T, c: T, r: T) -> DuffieSingletonRecoveryAndCdsBasisSpec<T> {
        DuffieSingletonRecoveryAndCdsBasisSpec {
            face_value: T::from_f64(100.0),
            coupon_rate_annual: c,
            payment_freq_per_year: 2,
            maturity_years: T::from_f64(3.0),
            risk_free_rate: r,
            hazard_rate: h,
            recovery_rate: T::from_f64(0.4),
            market_bond_price: T::from_f64(98.0),
            repo_funding_spread: T::from_f64(0.005),
        }
    }
    #[test]
    fn phase24() {
        let s = mk(0.02, 0.1, 0.02);
        let r = analyze_duffie_singleton_rmv_rfv_and_cds_basis(&s);
        assert!(r.rmv_bond_price > r.rfv_bond_price);
        let z = mk(0.02, 0.0, 0.02);
        let rz = analyze_duffie_singleton_rmv_rfv_and_cds_basis(&z);
        assert!(rz.rmv_bond_price < rz.rfv_bond_price);
        let h = 1e-5;
        let fd = (analyze_duffie_singleton_rmv_rfv_and_cds_basis(&mk(0.02 + h, 0.03, 0.02))
            .rmv_bond_price
            - analyze_duffie_singleton_rmv_rfv_and_cds_basis(&mk(0.02 - h, 0.03, 0.02))
                .rmv_bond_price)
            / (2.0 * h);
        let ad = analyze_duffie_singleton_rmv_rfv_and_cds_basis(&mk(
            Dual::<1>::variable(0.02, 0),
            Dual::from_value(0.03),
            Dual::from_value(0.02),
        ))
        .rmv_bond_price;
        assert!((ad.d[0] - fd).abs() < 1e-5);
        let mut r0s = mk(0.02, 0.03, 0.02);
        r0s.recovery_rate = 0.0;
        let r0 = analyze_duffie_singleton_rmv_rfv_and_cds_basis(&r0s);
        assert!((r0.rmv_bond_price - r0.rfv_bond_price).abs() < 1e-12);
        assert!((r0.rmv_bond_price - r0.rt_bond_price).abs() < 1e-12);
        let fd_basis = (analyze_duffie_singleton_rmv_rfv_and_cds_basis(&mk(0.02 + h, 0.03, 0.02))
            .raw_cds_bond_basis
            - analyze_duffie_singleton_rmv_rfv_and_cds_basis(&mk(0.02 - h, 0.03, 0.02))
                .raw_cds_bond_basis)
            / (2.0 * h);
        let ad_basis = analyze_duffie_singleton_rmv_rfv_and_cds_basis(&mk(
            Dual::<1>::variable(0.02, 0),
            Dual::from_value(0.03),
            Dual::from_value(0.02),
        ))
        .raw_cds_bond_basis;
        assert!((ad_basis.d[0] - fd_basis).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase25_tests {
    use super::*;
    use ql_ad::Dual;
    fn ecl<T: Real>(x: T) -> Ifrs9CeclMultiStageEclSpec<T> {
        Ifrs9CeclMultiStageEclSpec {
            annual_ttc_pds: [x, T::from_f64(0.04), T::from_f64(0.1)],
            annual_eads: [T::from_f64(100.0), T::from_f64(100.0), T::from_f64(100.0)],
            base_lgd: T::from_f64(0.45),
            lgd_macro_sensitivity: T::from_f64(0.05),
            vasicek_asset_correlation: T::from_f64(0.2),
            effective_interest_rate: T::from_f64(0.03),
            scenario_weights: [T::from_f64(0.3), T::from_f64(0.5), T::from_f64(0.2)],
            scenario_macro_shocks_z: [T::from_f64(-1.0), T::zero(), T::from_f64(1.0)],
            origination_lifetime_pd: T::from_f64(0.04),
            sicr_pd_multiplier_threshold: T::from_f64(1.5),
        }
    }
    #[test]
    fn ifrs9_dual_fd() {
        let h = 1e-5;
        let fd = (compute_ifrs9_cecl_multistage_expected_credit_loss(&ecl(0.02 + h))
            .stage1_12m_ecl
            .val()
            - compute_ifrs9_cecl_multistage_expected_credit_loss(&ecl(0.02 - h))
                .stage1_12m_ecl
                .val())
            / (2.0 * h);
        let ad =
            compute_ifrs9_cecl_multistage_expected_credit_loss(&ecl(Dual::<1>::variable(0.02, 0)))
                .stage1_12m_ecl;
        assert!((ad.d[0] - fd).abs() < 1e-5);
        let r = compute_ifrs9_cecl_multistage_expected_credit_loss(&ecl(0.02));
        assert!(r.stage1_12m_ecl <= r.stage2_lifetime_ecl);
        assert!(r.stage1_12m_ecl <= r.blended_ifrs9_ecl_provision);
        assert!(r.blended_ifrs9_ecl_provision <= r.stage2_lifetime_ecl);
        assert!(r.stage2_lifetime_ecl <= r.stage3_credit_impaired_ecl);
        assert!(r.stage3_net_eir_interest_income < r.stage1_gross_interest_income);
    }
    fn lando<T: Real>(x: T) -> LandoCoxAffineDefaultSpec<T> {
        LandoCoxAffineDefaultSpec {
            risk_free_rate: T::from_f64(0.02),
            base_default_hazard: T::from_f64(0.01),
            initial_credit_intensity_y0: x,
            credit_kappa: T::from_f64(0.5),
            credit_theta: T::from_f64(0.02),
            credit_sigma: T::from_f64(0.1),
            initial_liquidity_intensity_z0: T::from_f64(0.01),
            liquidity_kappa: T::from_f64(0.4),
            liquidity_theta: T::from_f64(0.01),
            liquidity_sigma: T::from_f64(0.1),
            recovery_rate: T::from_f64(0.4),
            maturity: T::from_f64(3.0),
        }
    }
    #[test]
    fn lando_dual_fd_and_spread_decomposition() {
        let mut zl = lando(0.02);
        zl.initial_liquidity_intensity_z0 = 0.0;
        zl.liquidity_theta = 0.0;
        let zero_liq = evaluate_lando_cox_affine_default_and_liquidity_decomposition(&zl);
        assert!(zero_liq.pure_liquidity_spread_bps.abs() < 1e-12);
        assert!(
            (zero_liq.total_yield_spread_bps
                - zero_liq.pure_credit_spread_bps
                - zero_liq.pure_liquidity_spread_bps)
                .abs()
                < 1e-12
        );
        let h = 1e-5;
        let fd = (evaluate_lando_cox_affine_default_and_liquidity_decomposition(&lando(0.02 + h))
            .survival_probability
            .val()
            - evaluate_lando_cox_affine_default_and_liquidity_decomposition(&lando(0.02 - h))
                .survival_probability
                .val())
            / (2.0 * h);
        let ad = evaluate_lando_cox_affine_default_and_liquidity_decomposition(&lando(
            Dual::<1>::variable(0.02, 0),
        ))
        .survival_probability;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase28_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(x: T) -> FedsDeSpiegeleerAt1CocoDualTriggerSpec<T> {
        FedsDeSpiegeleerAt1CocoDualTriggerSpec {
            notional: T::from_f64(100.0),
            coupon_rate: T::from_f64(0.08),
            maturity: T::from_f64(5.0),
            risk_free_rate: T::from_f64(0.03),
            initial_cet1_ratio: T::from_f64(0.12) + x,
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
        }
    }
    #[test]
    fn coco_invariants_and_dual() {
        let r = evaluate_feds_despiegeleer_at1_coco_dual_trigger(&mk(0.0));
        assert!(r.total_dual_trigger_prob >= r.mechanical_cet1_trigger_prob);
        assert!(r.at1_coco_credit_spread_bps > 0.0);
        assert!(r.equity_dilution_fraction > 0.0 && r.equity_dilution_fraction < 1.0);
        let h = 1e-5;
        let fd = (evaluate_feds_despiegeleer_at1_coco_dual_trigger(&mk(h)).at1_coco_fair_value
            - evaluate_feds_despiegeleer_at1_coco_dual_trigger(&mk(-h)).at1_coco_fair_value)
            / (2.0 * h);
        let ad = evaluate_feds_despiegeleer_at1_coco_dual_trigger(&mk(Dual::<1>::variable(0.0, 0)))
            .at1_coco_fair_value;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase29_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn lando_invariants_and_dual() {
        fn mk<T: Real>(x: T) -> LandoIsraelGeneratorLogarithmSpec<T> {
            LandoIsraelGeneratorLogarithmSpec {
                transition_matrix_3x3: [
                    [T::from_f64(0.9), T::from_f64(0.08), T::from_f64(0.02)],
                    [T::from_f64(0.03), T::from_f64(0.92), T::from_f64(0.05)],
                    [T::zero(), T::zero(), T::one()],
                ],
                series_order_k: 12,
                fractional_horizon_years: x,
                recovery_rate: T::from_f64(0.4),
                risk_free_rate: T::from_f64(0.03),
            }
        }
        let r = compute_lando_israel_generator_logarithm_and_migration(&mk(0.5));
        for row in r.fractional_transition_matrix {
            assert!((row.iter().sum::<f64>() - 1.0).abs() < 1e-10);
        }
        assert!(r.generator_row_sum_max_error < 1e-12);
        let h = 1e-5;
        let fd = (compute_lando_israel_generator_logarithm_and_migration(&mk(0.5 + h))
            .migration_adjusted_ig_spread_bps
            - compute_lando_israel_generator_logarithm_and_migration(&mk(0.5 - h))
                .migration_adjusted_ig_spread_bps)
            / (2.0 * h);
        let ad = compute_lando_israel_generator_logarithm_and_migration(&mk(Dual::<1>::variable(
            0.5, 0,
        )))
        .migration_adjusted_ig_spread_bps;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
