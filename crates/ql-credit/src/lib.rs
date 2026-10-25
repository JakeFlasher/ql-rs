use nalgebra::{DMatrix, DVector};
use ql_ad::{Real, TapeScope, TapeVar};
use ql_curves::YieldTermStructure;
use ql_math::brent_solve;
use ql_time::Date;
use thiserror::Error;

pub mod cdo;
pub mod cds_option;
pub use cdo::*;
pub use cds_option::*;

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
    let p_def_h = row_h[m - 1].max_val(T::zero()).min_val(T::from_f64(0.999999));
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
    let imp_spread = -p_surv_h.ln() / spec.horizon_years.max(1e-6) * (T::one() - spec.recovery_rate);

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
