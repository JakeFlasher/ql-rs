pub use crate::phase31::{
    DavidsonLevinCmoPacCompanionReport, DavidsonLevinCmoPacCompanionSpec,
    FasbAsu2024InducedConversionCpdiReport, FasbAsu2024InducedConversionCpdiSpec,
    evaluate_davidson_levin_cmo_pac_companion_waterfall,
    evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid,
};
pub use crate::phase33::{
    TuckmanSpecifiedPoolDollarRollReport, TuckmanSpecifiedPoolDollarRollSpec,
    evaluate_tuckman_specified_pool_payup_and_dollar_roll,
};
use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_inflation::ZeroInflationCurve;
use ql_math::brent_solve;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixedIncomeError {
    InvalidInput,
    SolverFailed,
}

pub fn sinking_notionals<T: Real>(
    n_periods: usize,
    coupon_rate: T,
    frequency: usize,
    initial_notional: T,
) -> Vec<T> {
    if n_periods == 0 {
        return vec![initial_notional];
    }
    let r = coupon_rate / frequency.max(1) as f64;
    let payment = if r.val().abs() < 1e-14 {
        initial_notional / n_periods as f64
    } else {
        initial_notional * r / (T::one() - (T::one() + r).powf(-(n_periods as f64)))
    };
    let mut out = Vec::with_capacity(n_periods + 1);
    out.push(initial_notional);
    for k in 0..n_periods {
        let next = if k + 1 == n_periods {
            T::zero()
        } else {
            out[k] * (T::one() + r) - payment
        };
        out.push(next.max_val(T::zero()));
    }
    out
}

#[derive(Clone, Debug)]
pub struct AmortizingFixedRateBond<T: Real> {
    pub initial_notional: T,
    pub coupon_rate: T,
    pub frequency: usize,
    pub payment_times: Vec<f64>,
}

pub fn price_amortizing_bond<T: Real>(
    bond: &AmortizingFixedRateBond<T>,
    curve: &impl YieldTermStructure<T>,
) -> T {
    let notionals = sinking_notionals(
        bond.payment_times.len(),
        bond.coupon_rate,
        bond.frequency,
        bond.initial_notional,
    );
    let mut pv = T::zero();
    for (i, &t) in bond.payment_times.iter().enumerate() {
        let coupon = notionals[i] * bond.coupon_rate / bond.frequency.max(1) as f64;
        let principal = notionals[i] - notionals[i + 1];
        pv += (coupon + principal) * curve.discount_factor(t);
    }
    pv
}

#[derive(Clone, Debug)]
pub struct CpiBondSpec<T: Real> {
    pub face_amount: T,
    pub base_cpi: T,
    pub real_coupon_rate: T,
    pub payment_times: Vec<f64>,
    pub accrual_fractions: Vec<f64>,
    pub deflation_floor: bool,
}

pub fn price_cpi_bond<T: Real>(
    spec: &CpiBondSpec<T>,
    nominal_discount: &impl YieldTermStructure<T>,
    inflation_curve: &ZeroInflationCurve<T>,
) -> T {
    let mut pv = T::zero();
    for (i, &t) in spec.payment_times.iter().enumerate() {
        // Keep month arithmetic in a wide integer domain.  Casting directly
        // to `u8` wrapped for maturities beyond ~21 years (e.g. month 256
        // became month 0), selecting the wrong seasonal factor.  The
        // inflation helper expects the lagged month `(index - 1) mod 12`;
        // month zero is retained as January for a zero-time payment.
        let idx = (t * 12.0).round();
        let month = if idx <= 0.0 {
            1u8
        } else {
            (((idx - 1.0).rem_euclid(12.0)) as u8) + 1
        };
        let ratio =
            inflation_curve.forward_cpi(t, Some(month)) / spec.base_cpi.max_val(T::from_f64(1e-12));
        let accrual = spec.accrual_fractions.get(i).copied().unwrap_or_else(|| {
            if i == 0 {
                t
            } else {
                t - spec.payment_times[i - 1]
            }
        });
        let mut cf = spec.face_amount * spec.real_coupon_rate * accrual * ratio;
        if i + 1 == spec.payment_times.len() {
            cf += spec.face_amount
                * if spec.deflation_floor {
                    ratio.max_val(T::one())
                } else {
                    ratio
                };
        }
        pv += cf * nominal_discount.discount_factor(t);
    }
    pv
}

pub fn breakeven_inflation_rate<T: Real>(
    spec: &CpiBondSpec<T>,
    nominal_price: T,
    nominal_discount: &impl YieldTermStructure<T>,
) -> T {
    let t_last = *spec.payment_times.last().unwrap_or(&1.0);
    if !t_last.is_finite() || t_last <= 0.0 || spec.face_amount.val().abs() <= 1e-14 {
        return T::zero();
    }
    let target = nominal_price
        / (spec.face_amount * nominal_discount.discount_factor(t_last)).max_val(T::from_f64(1e-12));
    let mut pi = target.powf(1.0 / t_last) - T::one();
    if spec.real_coupon_rate.val().abs() <= 1e-14 {
        return pi;
    }
    for _ in 0..4 {
        let mut pv = T::zero();
        let mut dpv = T::zero();
        for (i, &t) in spec.payment_times.iter().enumerate() {
            let acc = spec.accrual_fractions.get(i).copied().unwrap_or_else(|| {
                if i == 0 {
                    t
                } else {
                    t - spec.payment_times[i - 1]
                }
            });
            let growth = (T::one() + pi).powf(t);
            let dgrowth = (T::one() + pi).powf(t - 1.0) * t;
            let mut real_cf = spec.face_amount * spec.real_coupon_rate * acc;
            if i + 1 == spec.payment_times.len() {
                real_cf += spec.face_amount;
            }
            let df = nominal_discount.discount_factor(t);
            pv += real_cf * growth * df;
            dpv += real_cf * dgrowth * df;
        }
        pi -= (pv - nominal_price) / dpv.max_val(T::from_f64(1e-12));
    }
    pi
}

#[derive(Copy, Clone, Debug)]
pub enum PrepaymentModel<T: Real> {
    ConstantCpr(T),
    Psa { psa_multiplier: T },
}

#[derive(Clone, Debug)]
pub struct AgencyMbsSpec<T: Real> {
    pub original_balance: T,
    pub wac: T,
    pub pass_through_rate: T,
    pub term_months: usize,
    pub seasoning_months: usize,
    pub prepayment: PrepaymentModel<T>,
}

fn smm<T: Real>(model: PrepaymentModel<T>, age: usize) -> T {
    let cpr = match model {
        PrepaymentModel::ConstantCpr(x) => x,
        PrepaymentModel::Psa { psa_multiplier } => {
            psa_multiplier * 0.06 * (age as f64 / 30.0).min(1.0)
        }
    }
    .max_val(T::zero())
    .min_val(T::one());
    T::one() - (T::one() - cpr).powf(1.0 / 12.0)
}

pub fn price_agency_mbs<T: Real>(
    spec: &AgencyMbsSpec<T>,
    curve: &impl YieldTermStructure<T>,
    oas: T,
) -> (T, T) {
    if !spec.original_balance.val().is_finite()
        || spec.original_balance.val() <= 1e-14
        || spec.term_months == 0
    {
        return (T::zero(), T::zero());
    }
    let mut balance = spec.original_balance;
    let mut pv = T::zero();
    let mut wal = T::zero();
    let rw = spec.wac / 12.0;
    for m in 1..=spec.term_months {
        if balance.val() <= 1e-14 {
            break;
        }
        let rem = (spec.term_months - m + 1) as f64;
        let payment = if rw.val().abs() < 1e-14 {
            balance / rem
        } else {
            balance * rw / (T::one() - (T::one() + rw).powf(-rem))
        };
        let scheduled_principal = (payment - balance * rw).min_val(balance);
        let prepay = smm(spec.prepayment, spec.seasoning_months + m)
            * (balance - scheduled_principal).max_val(T::zero());
        let coupon = balance * spec.pass_through_rate / 12.0;
        let principal = scheduled_principal + prepay;
        let t = m as f64 / 12.0;
        let df = curve.discount_factor(t) * (-oas * t).exp();
        pv += (coupon + principal) * df;
        wal += T::from_f64(t) * principal;
        balance = (balance - principal).max_val(T::zero());
    }
    (pv, wal / spec.original_balance)
}

/// Davidson--Levin four-factor agency MBS strip analytics.
#[derive(Clone, Debug)]
pub struct DavidsonLevinIoPoStripMbsSpec<T: Real = f64> {
    pub notional: T,
    pub wac_annual: T,
    pub pass_through_coupon_annual: T,
    pub current_mortgage_rate: T,
    pub discount_rate_annual: T,
    pub remaining_months: usize,
    pub loan_age_months: usize,
    pub turnover_base_cpr: T,
    pub max_refi_cpr: T,
    pub refi_slope: T,
    pub refi_threshold_spread: T,
    pub burnout_exponent: T,
    pub rate_bump_bps: T,
}
#[derive(Clone, Debug)]
pub struct DavidsonLevinIoPoStripMbsReport<T: Real = f64> {
    pub io_strip_price: T,
    pub po_strip_price: T,
    pub pass_through_price: T,
    pub strip_replication_error: T,
    pub io_effective_duration: T,
    pub po_effective_duration: T,
    pub pass_through_effective_duration: T,
    pub io_effective_convexity: T,
    pub po_effective_convexity: T,
    pub initial_cpr: T,
    pub weighted_average_life_years: T,
}
fn price_davidson_levin_at<T: Real>(
    s: &DavidsonLevinIoPoStripMbsSpec<T>,
    dr: T,
) -> (T, T, T, T, T) {
    let n = s.remaining_months.max(1);
    let cg = s.wac_annual / 12.0;
    let cn = s.pass_through_coupon_annual / 12.0;
    let rm = s.current_mortgage_rate + dr;
    let yd = (s.discount_rate_annual + dr) / 12.0;
    let seas = [
        0.94, 0.92, 0.96, 1.0, 1.05, 1.10, 1.12, 1.10, 1.04, 0.98, 0.91, 0.88,
    ];
    let mut b = T::one();
    let mut p = s.notional;
    let mut io = T::zero();
    let mut po = T::zero();
    let mut pt = T::zero();
    let mut wal = T::zero();
    let mut cpr0 = T::zero();
    for m in 1..=n {
        let age = (s.loan_age_months + m) as f64;
        let a = T::from_f64((age / 30.0).min(1.0));
        let dc = s.wac_annual - rm - s.refi_threshold_spread;
        let sig = (T::one() + (-s.refi_slope * dc).exp()).powf(-1.0);
        let w = if m == 1 {
            T::one()
        } else {
            (s.burnout_exponent.max_val(T::zero()) * b.max_val(T::from_f64(1e-12)).ln()).exp()
        };
        let cpr = ((s.turnover_base_cpr + s.max_refi_cpr * sig * w)
            * a
            * T::from_f64(seas[(m - 1) % 12]))
        .min_val(T::from_f64(0.95));
        if m == 1 {
            cpr0 = cpr;
        }
        let smm = T::one() - (T::one() - cpr).powf(1.0 / 12.0);
        let rem = (n - m + 1) as f64;
        let payment = if cg.val().abs() < 1e-14 {
            p / T::from_f64(rem)
        } else {
            let growth = (T::one() + cg).powf(rem);
            p * cg * growth / (growth - T::one())
        };
        let interest = p * cg;
        let sched = (payment - interest).max_val(T::zero());
        let prepay = (p - sched).max_val(T::zero()) * smm;
        let po_cf = sched + prepay;
        let io_cf = p * cn;
        let df = (T::one() + yd).powf(-(m as f64));
        io += io_cf * df;
        po += po_cf * df;
        pt += (io_cf + po_cf) * df;
        wal += T::from_f64(m as f64 / 12.0) * po_cf / s.notional.max_val(T::from_f64(1e-12));
        p = (p - po_cf).max_val(T::zero());
        b *= T::one() - smm;
    }
    (io, po, pt, cpr0, wal)
}
pub fn analyze_davidson_levin_mbs_io_po_strips<T: Real>(
    s: &DavidsonLevinIoPoStripMbsSpec<T>,
) -> DavidsonLevinIoPoStripMbsReport<T> {
    let h = s.rate_bump_bps.max_val(T::one()) * 1e-4;
    let z = price_davidson_levin_at(s, T::zero());
    let up = price_davidson_levin_at(s, h);
    let dn = price_davidson_levin_at(s, -h);
    let dur = |u: T, d: T, v: T| (d - u) / (T::from_f64(2.0) * h * v.max_val(T::from_f64(1e-12)));
    let conv =
        |u: T, d: T, v: T| (u - T::from_f64(2.0) * v + d) / (h * h * v.max_val(T::from_f64(1e-12)));
    let io_d = dur(up.0, dn.0, z.0);
    let po_d = dur(up.1, dn.1, z.1);
    let pt_d = dur(up.2, dn.2, z.2);
    DavidsonLevinIoPoStripMbsReport {
        io_strip_price: z.0,
        po_strip_price: z.1,
        pass_through_price: z.2,
        strip_replication_error: z.0 + z.1 - z.2,
        io_effective_duration: io_d,
        po_effective_duration: po_d,
        pass_through_effective_duration: pt_d,
        io_effective_convexity: conv(up.0, dn.0, z.0),
        po_effective_convexity: conv(up.1, dn.1, z.1),
        initial_cpr: z.3,
        weighted_average_life_years: z.4,
    }
}

pub fn solve_mbs_oas(
    spec: &AgencyMbsSpec<f64>,
    curve: &impl YieldTermStructure<f64>,
    market_price: f64,
) -> Result<f64, FixedIncomeError> {
    let f = |x: f64| price_agency_mbs(spec, curve, x).0 - market_price;
    brent_solve(f, -1.0, 1.0, 1e-13, 200).map_err(|_| FixedIncomeError::SolverFailed)
}

#[derive(Clone, Debug)]
pub struct SpecifiedPoolMbsSpec<T: Real = f64> {
    pub original_balance: T,
    pub wac: T,
    pub pass_through_rate: T,
    pub term_months: usize,
    pub seasoning_months: usize,
    pub cpr_min: T,
    pub cpr_max: T,
    pub refi_slope: T,
    pub refi_inflection: T,
    pub primary_secondary_spread: T,
    pub cdr: T,
    pub loss_severity: T,
}
#[derive(Copy, Clone, Debug)]
pub struct SpecifiedPoolValuation<T: Real> {
    pub pv: T,
    pub wal: T,
    pub total_prepaid_principal: T,
    pub total_credit_loss: T,
}

pub fn price_specified_pool_mbs<T: Real>(
    spec: &SpecifiedPoolMbsSpec<T>,
    curve: &impl YieldTermStructure<T>,
    oas: T,
) -> SpecifiedPoolValuation<T> {
    if spec.term_months == 0 || spec.original_balance.val() <= 0.0 {
        return SpecifiedPoolValuation {
            pv: T::zero(),
            wal: T::zero(),
            total_prepaid_principal: T::zero(),
            total_credit_loss: T::zero(),
        };
    }
    let mut bal = spec.original_balance;
    let mut pv = T::zero();
    let mut wal_num = T::zero();
    let mut tot_prin = T::zero();
    let mut prep = T::zero();
    let mut loss = T::zero();
    let rw = spec.wac / 12.0;
    for m in 1..=spec.term_months {
        if bal.val() <= 1e-14 {
            break;
        }
        let t = m as f64 / 12.0;
        let age = (spec.seasoning_months + m) as f64;
        let r_mtg = curve.zero_rate(t) + spec.primary_secondary_spread;
        let inc = spec.wac - r_mtg;
        let mut cpr = spec.cpr_min
            + (spec.cpr_max - spec.cpr_min)
                / (T::one() + (-spec.refi_slope * (inc - spec.refi_inflection)).exp());
        cpr *= T::from_f64((age / 30.0).min(1.0));
        cpr = cpr.max_val(T::zero()).min_val(T::from_f64(0.99));
        let smm = T::one() - (T::one() - cpr).powf(1.0 / 12.0);
        let cdr = spec.cdr.max_val(T::zero()).min_val(T::from_f64(0.99));
        let mdr = T::one() - (T::one() - cdr).powf(1.0 / 12.0);
        let def = bal * mdr;
        let credit_loss = def * spec.loss_severity;
        let perf = bal - def;
        let rem = (spec.term_months - m + 1) as f64;
        let pmt = if rw.val().abs() < 1e-14 {
            perf / rem
        } else {
            perf * rw / (T::one() - (T::one() + rw).powf(-rem))
        };
        let sched = (pmt - perf * rw).max_val(T::zero()).min_val(perf);
        let pre = smm * (perf - sched).max_val(T::zero());
        let coupon = perf * spec.pass_through_rate / 12.0;
        let prin_cf = sched + pre + (T::one() - spec.loss_severity) * def;
        let cf = coupon + prin_cf;
        pv += cf * curve.discount_factor(t) * (-oas * t).exp();
        wal_num += T::from_f64(t) * prin_cf;
        tot_prin += prin_cf;
        prep += pre;
        loss += credit_loss;
        bal = (perf - sched - pre).max_val(T::zero());
    }
    SpecifiedPoolValuation {
        pv,
        wal: wal_num / tot_prin.max_val(T::from_f64(1e-12)),
        total_prepaid_principal: prep,
        total_credit_loss: loss,
    }
}

#[derive(Copy, Clone, Debug)]
pub struct TbaDollarRollSpec<T: Real = f64> {
    pub front_clean_price: T,
    pub back_clean_price: T,
    pub wac: T,
    pub pass_through_rate: T,
    pub remaining_months: usize,
    pub front_accrued_fraction: f64,
    pub back_accrued_fraction: f64,
    pub roll_period_years: f64,
    pub assumed_cpr: T,
    pub repo_rate: T,
}
#[derive(Copy, Clone, Debug)]
pub struct TbaDollarRollReport<T: Real> {
    pub price_drop: T,
    pub breakeven_drop: T,
    pub implied_repo_rate: T,
    pub net_roll_advantage: T,
    pub surviving_factor: T,
}
pub fn analyze_tba_dollar_roll<T: Real>(s: &TbaDollarRollSpec<T>) -> TbaDollarRollReport<T> {
    let rem = T::from_f64(s.remaining_months.max(1) as f64);
    let rw = s.wac / 12.0;
    let pmt = if rw.val().abs() < 1e-14 {
        T::one() / rem
    } else {
        rw / (T::one() - (T::one() + rw).powf(-rem.val()))
    };
    let a_sched = (pmt - rw).max_val(T::zero()).min_val(T::one());
    let smm = T::one()
        - (T::one() - s.assumed_cpr.max_val(T::zero()).min_val(T::from_f64(0.99))).powf(1.0 / 12.0);
    let surv = (T::one() - a_sched) * (T::one() - smm);
    let monthly_coupon = s.pass_through_rate / 12.0 * 100.0;
    let i1 = s.front_clean_price + monthly_coupon * s.front_accrued_fraction;
    let i2 = surv * (s.back_clean_price + monthly_coupon * s.back_accrued_fraction);
    let cf = monthly_coupon + (T::one() - surv) * 100.0;
    let dt = s.roll_period_years.max(1e-6);
    let imp = (cf + i2 - i1) / (i1 * dt);
    let net = i1 * (T::one() + s.repo_rate * dt) - i2 - cf;
    let pstar = (i1 * (T::one() + s.repo_rate * dt) - cf) / surv.max_val(T::from_f64(1e-12))
        - monthly_coupon * s.back_accrued_fraction;
    TbaDollarRollReport {
        price_drop: s.front_clean_price - s.back_clean_price,
        breakeven_drop: s.front_clean_price - pstar,
        implied_repo_rate: imp,
        net_roll_advantage: net,
        surviving_factor: surv,
    }
}

#[derive(Clone, Debug)]
pub struct MbsBurnoutOasSpec<T: Real = f64> {
    pub market_price: T,
    pub original_balance: T,
    pub wac: T,
    pub pass_through_rate: T,
    pub term_months: usize,
    pub seasoning_months: usize,
    pub turnover_cpr: T,
    pub max_refi_cpr: T,
    pub refi_slope: T,
    pub refi_incentive_threshold: T,
    pub burnout_beta: T,
    pub primary_secondary_spread: T,
}

#[derive(Copy, Clone, Debug)]
pub struct MbsBurnoutOasReport<T: Real> {
    pub oas: T,
    pub model_price_zero_oas: T,
    pub wal: T,
    pub terminal_burnout_factor: T,
    pub option_adjusted_duration: T,
    pub option_adjusted_convexity: T,
}

pub fn price_mbs_with_burnout<T: Real>(
    spec: &MbsBurnoutOasSpec<T>,
    curve: &impl YieldTermStructure<T>,
    oas: T,
) -> (T, T, T) {
    let n = spec.term_months.max(1);
    let mut bal = spec.original_balance;
    let mut pv = T::zero();
    let mut prin_time = T::zero();
    let mut prin_total = T::zero();
    let mut cum_burnout = T::zero();
    let mut burnout_factor = T::one();
    let eps = T::from_f64(1e-4);

    for m in 1..=n {
        let t = m as f64 / 12.0;
        let rem = (n - m + 1) as f64;
        let r_m = spec.wac / 12.0;
        let pmt_factor = if r_m.val().abs() < 1e-12 {
            T::one() / rem
        } else {
            let one_plus = (T::one() + r_m).max_val(T::from_f64(1e-8));
            r_m / (T::one() - (-one_plus.ln() * rem).exp())
        };
        let sched_rate = (pmt_factor - r_m).max_val(T::zero()).min_val(T::one());
        let df_base = curve.discount_factor(t).max_val(T::from_f64(1e-12));
        let z_t = -df_base.ln() / t.max(1e-8);
        let incentive = spec.wac - (z_t + spec.primary_secondary_spread);
        burnout_factor = (-spec.burnout_beta * cum_burnout).exp();
        let refi_arg = (-spec.refi_slope * (incentive - spec.refi_incentive_threshold))
            .max_val(T::from_f64(-40.0))
            .min_val(T::from_f64(40.0));
        let refi_cpr = burnout_factor * spec.max_refi_cpr / (T::one() + refi_arg.exp());
        let psa_ramp = ((spec.seasoning_months + m) as f64 / 30.0).min(1.0);
        let cpr = ((spec.turnover_cpr + refi_cpr) * psa_ramp)
            .max_val(T::zero())
            .min_val(T::from_f64(0.98));
        let smm = T::one() - ((T::one() - cpr).max_val(T::from_f64(1e-8)).ln() / 12.0).exp();

        let interest_pass = bal * (spec.pass_through_rate / 12.0);
        let sched_prin = bal * sched_rate;
        let prepay_prin = (bal - sched_prin) * smm;
        let total_prin = sched_prin + prepay_prin;
        let cf = interest_pass + total_prin;
        let df = df_base * (-oas * t).exp();
        pv += cf * df;
        prin_time += total_prin * t;
        prin_total += total_prin;
        bal = (bal - total_prin).max_val(T::zero());

        let pos_inc = (incentive + (incentive * incentive + eps * eps).sqrt()) * 0.5;
        cum_burnout += pos_inc / 12.0;
    }
    let wal = prin_time / prin_total.max_val(T::from_f64(1e-12));
    (pv, wal, burnout_factor)
}

pub fn analyze_mbs_burnout_oas_ad<T: Real>(
    spec: &MbsBurnoutOasSpec<T>,
    curve: &impl YieldTermStructure<T>,
) -> MbsBurnoutOasReport<T> {
    struct ScalarDiscount<'a> {
        dfs: &'a [f64],
    }
    impl YieldTermStructure<f64> for ScalarDiscount<'_> {
        fn reference_date(&self) -> ql_time::Date {
            ql_time::Date(0)
        }
        fn discount_factor(&self, t: f64) -> f64 {
            let idx = ((t * 12.0).round() as usize).saturating_sub(1);
            self.dfs
                .get(idx)
                .copied()
                .unwrap_or_else(|| self.dfs.last().copied().unwrap_or(1.0))
        }
        fn zero_rate(&self, t: f64) -> f64 {
            let tt = t.max(1e-8);
            -self.discount_factor(t).ln() / tt
        }
        fn instantaneous_forward(&self, t: f64) -> f64 {
            self.zero_rate(t)
        }
    }

    let n = spec.term_months.max(1);
    let dfs_f64: Vec<f64> = (1..=n)
        .map(|m| curve.discount_factor(m as f64 / 12.0).val())
        .collect();
    let sc_curve = ScalarDiscount { dfs: &dfs_f64 };
    let spec_f64 = MbsBurnoutOasSpec {
        market_price: spec.market_price.val(),
        original_balance: spec.original_balance.val(),
        wac: spec.wac.val(),
        pass_through_rate: spec.pass_through_rate.val(),
        term_months: n,
        seasoning_months: spec.seasoning_months,
        turnover_cpr: spec.turnover_cpr.val(),
        max_refi_cpr: spec.max_refi_cpr.val(),
        refi_slope: spec.refi_slope.val(),
        refi_incentive_threshold: spec.refi_incentive_threshold.val(),
        burnout_beta: spec.burnout_beta.val(),
        primary_secondary_spread: spec.primary_secondary_spread.val(),
    };

    let mut lo = -0.10_f64;
    let mut hi = 0.25_f64;
    for _ in 0..48 {
        let mid = 0.5 * (lo + hi);
        let p_mid = price_mbs_with_burnout(&spec_f64, &sc_curve, mid).0;
        if p_mid > spec_f64.market_price {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mut oas_root = 0.5 * (lo + hi);
    let h_fd = 1e-5;
    for _ in 0..5 {
        let p0 = price_mbs_with_burnout(&spec_f64, &sc_curve, oas_root).0 - spec_f64.market_price;
        let pp = price_mbs_with_burnout(&spec_f64, &sc_curve, oas_root + h_fd).0;
        let pm = price_mbs_with_burnout(&spec_f64, &sc_curve, oas_root - h_fd).0;
        let dp = (pp - pm) / (2.0 * h_fd);
        if dp.abs() > 1e-12 {
            oas_root = (oas_root - p0 / dp).clamp(-0.15, 0.35);
        }
    }
    let pp_f64 = price_mbs_with_burnout(&spec_f64, &sc_curve, oas_root + h_fd).0;
    let p0_f64 = price_mbs_with_burnout(&spec_f64, &sc_curve, oas_root).0;
    let pm_f64 = price_mbs_with_burnout(&spec_f64, &sc_curve, oas_root - h_fd).0;
    let f_xx = (pp_f64 - 2.0 * p0_f64 + pm_f64) / (h_fd * h_fd);

    let oas_t = T::from_f64(oas_root);
    let h_t = T::from_f64(h_fd);
    let (p_up_oas, _, _) = price_mbs_with_burnout(spec, curve, oas_t + h_t);
    let (p_dn_oas, _, _) = price_mbs_with_burnout(spec, curve, oas_t - h_t);
    let f_x = (p_up_oas - p_dn_oas) / (2.0 * h_fd);

    let oas = ql_ad::ift_lift_1d(
        oas_root,
        |s| price_mbs_with_burnout(spec, curve, s).0 - spec.market_price,
        f_x,
        f_xx,
    );
    let (model_price_zero_oas, _, _) = price_mbs_with_burnout(spec, curve, T::zero());
    let (p_mid_t, wal, terminal_burnout_factor) = price_mbs_with_burnout(spec, curve, oas);

    let dy = 1e-4;
    struct ShiftedCurve<'a, C> {
        base: &'a C,
        shift: f64,
    }
    impl<T: Real, C: YieldTermStructure<T>> YieldTermStructure<T> for ShiftedCurve<'_, C> {
        fn reference_date(&self) -> ql_time::Date {
            self.base.reference_date()
        }
        fn discount_factor(&self, t: f64) -> T {
            self.base.discount_factor(t) * (-self.shift * t).exp()
        }
        fn zero_rate(&self, t: f64) -> T {
            self.base.zero_rate(t) + T::from_f64(self.shift)
        }
        fn instantaneous_forward(&self, t: f64) -> T {
            self.base.instantaneous_forward(t) + T::from_f64(self.shift)
        }
    }
    let c_up = ShiftedCurve {
        base: curve,
        shift: dy,
    };
    let c_dn = ShiftedCurve {
        base: curve,
        shift: -dy,
    };
    let (p_rate_up, _, _) = price_mbs_with_burnout(spec, &c_up, oas);
    let (p_rate_dn, _, _) = price_mbs_with_burnout(spec, &c_dn, oas);
    let p_denom = p_mid_t.max_val(T::from_f64(1e-8));
    let oad = (p_rate_dn - p_rate_up) / (p_denom * (2.0 * dy));
    let oac = (p_rate_up - p_mid_t * 2.0 + p_rate_dn) / (p_denom * (dy * dy));

    MbsBurnoutOasReport {
        oas,
        model_price_zero_oas,
        wal,
        terminal_burnout_factor,
        option_adjusted_duration: oad,
        option_adjusted_convexity: oac,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct TuckmanBurghardtTreasuryWildcardBasisSpec<T: Real = f64> {
    pub futures_settlement_price: T,
    pub delivery_horizon_years: T,
    pub general_collateral_repo_rate: T,
    pub repo_specialness_spread: T,
    pub bond1_clean_price: T,
    pub bond1_conversion_factor: T,
    pub bond1_modified_duration: T,
    pub bond1_yield_beta: T,
    pub bond2_clean_price: T,
    pub bond2_conversion_factor: T,
    pub bond2_modified_duration: T,
    pub bond2_yield_beta: T,
    pub afternoon_wildcard_vol_annual: T,
    pub eom_delivery_window_days: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct TuckmanBurghardtTreasuryWildcardBasisReport<T: Real = f64> {
    pub bond1_cf_adjusted_forward: T,
    pub bond2_cf_adjusted_forward: T,
    pub special_repo_carry_pickup: T,
    pub critical_yield_shift_ctd_switch_bps: T,
    pub quality_switch_option_value: T,
    pub afternoon_wildcard_option_value: T,
    pub eom_timing_carry_option_value: T,
    pub theoretical_net_basis_ctd: T,
}
pub fn evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch<T: Real>(
    s: &TuckmanBurghardtTreasuryWildcardBasisSpec<T>,
) -> TuckmanBurghardtTreasuryWildcardBasisReport<T> {
    const Z_GH9: [f64; 9] = [
        -4.51274586339978,
        -3.20542900285647,
        -2.07684798313409,
        -1.02325566378913,
        0.0,
        1.02325566378913,
        2.07684798313409,
        3.20542900285647,
        4.51274586339978,
    ];
    const W_GH9: [f64; 9] = [
        0.00002234584404,
        0.00278914132131,
        0.04991640676364,
        0.24409745015871,
        0.40634931182460,
        0.24409745015871,
        0.04991640676364,
        0.00278914132131,
        0.00002234584404,
    ];
    let t = s.delivery_horizon_years.max_val(T::from_f64(1e-8));
    let carry = (s.general_collateral_repo_rate - s.repo_specialness_spread) * t;
    let cf1 = s.bond1_conversion_factor.max_val(T::from_f64(1e-8));
    let cf2 = s.bond2_conversion_factor.max_val(T::from_f64(1e-8));
    let f1 = s.bond1_clean_price * (T::one() + carry) / cf1;
    let f2 = s.bond2_clean_price * (T::one() + carry) / cf2;
    let b1 = f1 - s.futures_settlement_price;
    let b2 = f2 - s.futures_settlement_price;
    let dv1 = f1 * s.bond1_modified_duration * s.bond1_yield_beta;
    let dv2 = f2 * s.bond2_modified_duration * s.bond2_yield_beta;
    let den = (dv1 - dv2).abs().max_val(T::from_f64(1e-6));
    let shift_bps = ((f2 - f1).abs() / den) * T::from_f64(1e4);
    let sig_y = s.afternoon_wildcard_vol_annual.abs() * T::from_f64(0.1) * t.sqrt();
    let base_min = f1.min_val(f2);
    let mut qopt = T::zero();
    for k in 0..9 {
        let dy = sig_y * T::from_f64(Z_GH9[k]);
        let f1_k = f1 - dv1 * dy;
        let f2_k = f2 - dv2 * dy;
        let diff = base_min - f1_k.min_val(f2_k);
        let pos = T::from_f64(0.5) * (diff + (diff * diff + T::from_f64(1e-12)).sqrt());
        qopt += T::from_f64(W_GH9[k]) * pos;
    }
    let cf_ctd = cf1.min_val(cf2);
    let tail = (T::one() - cf_ctd).abs() / cf_ctd;
    let wild = s.afternoon_wildcard_vol_annual.abs()
        * (t / T::from_f64(252.0)).sqrt()
        * s.futures_settlement_price
        * (T::from_f64(0.1) + tail);
    let eom = s.futures_settlement_price
        * s.afternoon_wildcard_vol_annual.abs()
        * T::from_f64(s.eom_delivery_window_days.max(1) as f64 / 252.0).sqrt()
        * T::from_f64(0.01);
    TuckmanBurghardtTreasuryWildcardBasisReport {
        bond1_cf_adjusted_forward: f1,
        bond2_cf_adjusted_forward: f2,
        special_repo_carry_pickup: s.repo_specialness_spread * t * s.futures_settlement_price,
        critical_yield_shift_ctd_switch_bps: shift_bps,
        quality_switch_option_value: qopt,
        afternoon_wildcard_option_value: wild,
        eom_timing_carry_option_value: eom,
        theoretical_net_basis_ctd: b1.min_val(b2) - qopt - wild - eom,
    }
}

#[derive(Clone, Debug)]
pub struct VeronesiFabozziConvexityPrepaymentSpec<T: Real = f64> {
    pub forward_rate: T,
    pub short_rate_vol: T,
    pub mean_reversion_kappa: T,
    pub fixing_time_t1: T,
    pub maturity_time_t2: T,
    pub wac_coupon: T,
    pub current_mortgage_rate: T,
    pub min_cpr_turnover: T,
    pub max_cpr_refi: T,
    pub refi_steepness_gamma: T,
    pub refi_incentive_threshold: T,
    pub pool_age_months: T,
    pub psa_ramp_months: T,
}

#[derive(Clone, Debug)]
pub struct VeronesiFabozziConvexityPrepaymentReport<T: Real = f64> {
    pub vasicek_b_t1_t2: T,
    pub futures_forward_convexity_bias_bps: T,
    pub implied_futures_rate: T,
    pub psa_seasoning_multiplier: T,
    pub refi_incentive_spread: T,
    pub scurve_cpr: T,
    pub single_monthly_mortality_smm: T,
    pub effective_prepayment_duration: T,
}

pub fn evaluate_veronesi_fabozzi_futures_convexity_and_prepayment<T: Real>(
    s: &VeronesiFabozziConvexityPrepaymentSpec<T>,
) -> VeronesiFabozziConvexityPrepaymentReport<T> {
    let k = s.mean_reversion_kappa.max_val(T::from_f64(1e-10));
    let dt = (s.maturity_time_t2 - s.fixing_time_t1).max_val(T::from_f64(1e-10));
    let e1 = (-k * dt).exp();
    let e0 = (-k * s.fixing_time_t1).exp();
    let b12 = (T::one() - e1) / k;
    let b01 = (T::one() - e0) / k;
    let e2 = (-T::from_f64(2.0) * k * s.fixing_time_t1).exp();
    let bracket = b01 * b01 + b12 / (T::from_f64(2.0) * k) * (T::one() - e2);
    let bias = s.short_rate_vol * s.short_rate_vol / 2.0 * b12 / dt * bracket;
    let age = s.pool_age_months.max_val(T::zero());
    let ramp = s.psa_ramp_months.max_val(T::from_f64(1e-10));
    let u = age / ramp;
    let smooth_min =
        |a: T, b: T| (a + b - ((a - b) * (a - b) + T::from_f64(1e-8)).sqrt()) / T::from_f64(2.0);
    let smooth_max =
        |a: T, b: T| (a + b + ((a - b) * (a - b) + T::from_f64(1e-8)).sqrt()) / T::from_f64(2.0);
    let psa = smooth_max(T::zero(), smooth_min(u, T::one()));
    let incentive = s.wac_coupon - s.current_mortgage_rate;
    let logistic = T::one()
        / (T::one() + (-s.refi_steepness_gamma * (incentive - s.refi_incentive_threshold)).exp());
    let cpr = psa
        * (s.min_cpr_turnover + (s.max_cpr_refi - s.min_cpr_turnover) * logistic)
            .max_val(T::zero())
            .min_val(T::from_f64(0.999999));
    let smm = T::one() - (T::one() - cpr).powf(1.0 / 12.0);
    VeronesiFabozziConvexityPrepaymentReport {
        vasicek_b_t1_t2: b12,
        futures_forward_convexity_bias_bps: bias * 10000.0,
        implied_futures_rate: s.forward_rate + bias,
        psa_seasoning_multiplier: psa,
        refi_incentive_spread: incentive,
        scurve_cpr: cpr,
        single_monthly_mortality_smm: smm,
        effective_prepayment_duration: {
            let rem_years = (T::from_f64(30.0) - age / T::from_f64(12.0)).max_val(T::one());
            let eff_hazard = (s.current_mortgage_rate + cpr).max_val(T::from_f64(1e-6));
            (T::one() - (-eff_hazard * rem_years).exp()) / eff_hazard
        },
    }
}

#[derive(Clone, Debug)]
pub struct TreasuryFutureDeliverable<T: Real = f64> {
    pub clean_price: T,
    pub coupon_rate: T,
    pub remaining_semiannual_periods: usize,
    pub accrued_interest_spot: T,
    pub accrued_interest_delivery: T,
    pub interim_coupon: T,
    pub yield_vol: T,
    pub modified_duration: T,
}

#[derive(Clone, Debug)]
pub struct TreasuryFuturesBasketSpec<T: Real = f64> {
    pub futures_price: T,
    pub delivery_time: f64,
    pub repo_rate: T,
    pub notional_standard_yield: f64,
    pub deliverables: Vec<TreasuryFutureDeliverable<T>>,
}

#[derive(Copy, Clone, Debug)]
pub struct DeliverableBondMetrics<T: Real> {
    pub conversion_factor: T,
    pub invoice_price: T,
    pub gross_basis: T,
    pub carry: T,
    pub net_basis: T,
    pub implied_repo_rate: T,
    pub forward_cf_adjusted_price: T,
}

#[derive(Clone, Debug)]
pub struct TreasuryFuturesCtdAnalysis<T: Real> {
    pub ctd_index: usize,
    pub ctd_conversion_factor: T,
    pub ctd_invoice_price: T,
    pub ctd_net_basis: T,
    pub ctd_implied_repo_rate: T,
    pub switch_option_value: T,
    pub option_adjusted_futures_price: T,
    pub bond_metrics: Vec<DeliverableBondMetrics<T>>,
}

pub fn analyze_treasury_futures_ctd<T: Real>(
    spec: &TreasuryFuturesBasketSpec<T>,
) -> TreasuryFuturesCtdAnalysis<T> {
    if spec.deliverables.is_empty() {
        return TreasuryFuturesCtdAnalysis {
            ctd_index: 0,
            ctd_conversion_factor: T::one(),
            ctd_invoice_price: T::zero(),
            ctd_net_basis: T::zero(),
            ctd_implied_repo_rate: T::zero(),
            switch_option_value: T::zero(),
            option_adjusted_futures_price: T::zero(),
            bond_metrics: Vec::new(),
        };
    }
    let dt = spec.delivery_time.max(1e-6);
    let y_semi = (spec.notional_standard_yield.max(1e-6)) * 0.5;
    let disc_semi = 1.0 / (1.0 + y_semi);
    let mut metrics = Vec::with_capacity(spec.deliverables.len());
    let mut best_idx = 0usize;
    let mut best_irr = f64::NEG_INFINITY;
    for (idx, b) in spec.deliverables.iter().enumerate() {
        let m_per = b.remaining_semiannual_periods.max(1);
        let df_n = disc_semi.powi(m_per as i32);
        let ann_n = (1.0 - df_n) / y_semi;
        let cf = b.coupon_rate * 0.5 * ann_n + df_n;
        let invoice = spec.futures_price * cf + b.accrued_interest_delivery;
        let dirty = b.clean_price + b.accrued_interest_spot;
        let gross_basis = b.clean_price - spec.futures_price * cf;
        let carry = (b.accrued_interest_delivery - b.accrued_interest_spot + b.interim_coupon)
            - dirty * spec.repo_rate * dt;
        let net_basis = gross_basis - carry;
        let irr = (invoice + b.interim_coupon - dirty) / (dirty.max_val(T::from_f64(1e-8)) * dt);
        let fwd_clean = b.clean_price - carry;
        let fwd_adj = fwd_clean / cf.max_val(T::from_f64(1e-8));
        if irr.val() > best_irr {
            best_irr = irr.val();
            best_idx = idx;
        }
        metrics.push(DeliverableBondMetrics {
            conversion_factor: cf,
            invoice_price: invoice,
            gross_basis,
            carry,
            net_basis,
            implied_repo_rate: irr,
            forward_cf_adjusted_price: fwd_adj,
        });
    }
    let ctd = metrics[best_idx];
    let switch_val = if metrics.len() >= 2 {
        let mut runner_idx = if best_idx == 0 { 1 } else { 0 };
        for i in 0..metrics.len() {
            if i != best_idx
                && metrics[i].forward_cf_adjusted_price.val()
                    < metrics[runner_idx].forward_cf_adjusted_price.val()
            {
                runner_idx = i;
            }
        }
        let b1 = &spec.deliverables[best_idx];
        let b2 = &spec.deliverables[runner_idx];
        let m1 = ctd.forward_cf_adjusted_price;
        let m2 = metrics[runner_idx].forward_cf_adjusted_price;
        let s1 = m1 * b1.modified_duration * b1.yield_vol;
        let s2 = m2 * b2.modified_duration * b2.yield_vol;
        let sig_sp =
            ((s1 * s1 + s2 * s2 - s1 * s2 * 1.8).max_val(T::from_f64(1e-12))).sqrt() * dt.sqrt();
        let diff = m1 - m2;
        let spread = -((diff * diff + 1e-12).sqrt());
        let d = spread / sig_sp;
        (spread * d.norm_cdf() + sig_sp * d.norm_pdf()).max_val(T::zero())
    } else {
        T::zero()
    };
    TreasuryFuturesCtdAnalysis {
        ctd_index: best_idx,
        ctd_conversion_factor: ctd.conversion_factor,
        ctd_invoice_price: ctd.invoice_price,
        ctd_net_basis: ctd.net_basis,
        ctd_implied_repo_rate: ctd.implied_repo_rate,
        switch_option_value: switch_val,
        option_adjusted_futures_price: ctd.forward_cf_adjusted_price - switch_val,
        bond_metrics: metrics,
    }
}

#[derive(Clone, Debug)]
pub struct CmoTrancheSpec<T: Real = f64> {
    pub initial_balance: T,
    pub coupon_rate: T,
    pub is_pac: bool,
}

#[derive(Clone, Debug)]
pub struct CmoWaterfallSpec<T: Real = f64> {
    pub original_balance: T,
    pub wac: T,
    pub pass_through_rate: T,
    pub term_months: usize,
    pub psa_speed: T,
    pub pac_collar_psa: Option<(T, T)>,
    pub tranches: Vec<CmoTrancheSpec<T>>,
}

#[derive(Copy, Clone, Debug)]
pub struct CmoTrancheValuation<T: Real> {
    pub pv: T,
    pub wal: T,
    pub total_principal: T,
    pub total_interest: T,
}

#[derive(Clone, Debug)]
pub struct CmoWaterfallValuation<T: Real> {
    pub pool_pv: T,
    pub pool_wal: T,
    pub tranches: Vec<CmoTrancheValuation<T>>,
    pub io_strip_pv: T,
    pub total_structure_pv: T,
}

pub fn price_cmo_waterfall<T: Real>(
    s: &CmoWaterfallSpec<T>,
    curve: &impl YieldTermStructure<T>,
    oas: T,
) -> CmoWaterfallValuation<T> {
    let n_m = s.term_months.max(1);
    let rw = s.wac / 12.0;
    let sim_pool_prin = |psa_mult: T| -> Vec<T> {
        let mut bal = s.original_balance;
        let mut flows = Vec::with_capacity(n_m);
        for m in 1..=n_m {
            let rem = (n_m - m + 1) as f64;
            let pmt = if rw.val().abs() < 1e-14 {
                bal / rem
            } else {
                bal * rw / (T::one() - (T::one() + rw).powf(-rem))
            };
            let sched = (pmt - bal * rw).max_val(T::zero()).min_val(bal);
            let base_cpr = (0.002 * m as f64).min(0.06);
            let psa_scale = if psa_mult.val() > 10.0 { 0.01 } else { 1.0 };
            let cpr = (psa_mult * (base_cpr * psa_scale))
                .max_val(T::zero())
                .min_val(T::from_f64(0.99));
            let smm = T::one() - (T::one() - cpr).powf(1.0 / 12.0);
            let pre = smm * (bal - sched).max_val(T::zero());
            let prin = (sched + pre).min_val(bal);
            flows.push(prin);
            bal = (bal - prin).max_val(T::zero());
        }
        flows
    };
    let pac_sched = if let Some((psa_lo, psa_hi)) = s.pac_collar_psa {
        let f_lo = sim_pool_prin(psa_lo);
        let f_hi = sim_pool_prin(psa_hi);
        f_lo.into_iter()
            .zip(f_hi)
            .map(|(a, b)| a.min_val(b))
            .collect::<Vec<_>>()
    } else {
        vec![T::zero(); n_m]
    };
    let pool_prin = sim_pool_prin(s.psa_speed);
    let k_tr = s.tranches.len();
    let mut tr_bal: Vec<T> = s.tranches.iter().map(|t| t.initial_balance).collect();
    let mut tr_pv = vec![T::zero(); k_tr];
    let mut tr_wal_num = vec![T::zero(); k_tr];
    let mut tr_prin_tot = vec![T::zero(); k_tr];
    let mut tr_int_tot = vec![T::zero(); k_tr];
    let mut io_pv = T::zero();
    let mut pool_pv = T::zero();
    let mut pool_wal_num = T::zero();
    let mut pool_bal = s.original_balance;
    for m in 1..=n_m {
        let t = m as f64 / 12.0;
        let df = curve.discount_factor(t) * (-oas * t).exp();
        let p_tot = pool_prin[m - 1];
        let pool_int = pool_bal * s.pass_through_rate / 12.0;
        pool_pv += (p_tot + pool_int) * df;
        pool_wal_num += p_tot * t;
        let mut tr_int_raw = vec![T::zero(); k_tr];
        let mut tr_int_sum = T::zero();
        for i in 0..k_tr {
            let icf = tr_bal[i] * s.tranches[i].coupon_rate / 12.0;
            tr_int_raw[i] = icf;
            tr_int_sum += icf;
        }
        let int_scale = if tr_int_sum.val() > pool_int.val() && tr_int_sum.val() > 1e-14 {
            pool_int / tr_int_sum
        } else {
            T::one()
        };
        let mut paid_int_sum = T::zero();
        for i in 0..k_tr {
            let icf = tr_int_raw[i] * int_scale;
            tr_int_tot[i] += icf;
            tr_pv[i] += icf * df;
            paid_int_sum += icf;
        }
        io_pv += (pool_int - paid_int_sum).max_val(T::zero()) * df;
        let mut rem_p = p_tot;
        if s.pac_collar_psa.is_some() {
            let mut rem_target = pac_sched[m - 1];
            for i in 0..k_tr {
                if s.tranches[i].is_pac {
                    let alloc = rem_p.min_val(rem_target).min_val(tr_bal[i]);
                    tr_bal[i] = (tr_bal[i] - alloc).max_val(T::zero());
                    rem_p = (rem_p - alloc).max_val(T::zero());
                    rem_target = (rem_target - alloc).max_val(T::zero());
                    tr_pv[i] += alloc * df;
                    tr_wal_num[i] += alloc * t;
                    tr_prin_tot[i] += alloc;
                }
            }
            for i in 0..k_tr {
                if !s.tranches[i].is_pac {
                    let alloc = rem_p.min_val(tr_bal[i]);
                    tr_bal[i] = (tr_bal[i] - alloc).max_val(T::zero());
                    rem_p = (rem_p - alloc).max_val(T::zero());
                    tr_pv[i] += alloc * df;
                    tr_wal_num[i] += alloc * t;
                    tr_prin_tot[i] += alloc;
                }
            }
        }
        for i in 0..k_tr {
            if rem_p.val() <= 1e-14 {
                break;
            }
            let alloc = rem_p.min_val(tr_bal[i]);
            tr_bal[i] = (tr_bal[i] - alloc).max_val(T::zero());
            rem_p = (rem_p - alloc).max_val(T::zero());
            tr_pv[i] += alloc * df;
            tr_wal_num[i] += alloc * t;
            tr_prin_tot[i] += alloc;
        }
        pool_bal = (pool_bal - p_tot).max_val(T::zero());
    }
    let mut tranches = Vec::with_capacity(k_tr);
    let mut sum_tr_pv = io_pv;
    for i in 0..k_tr {
        sum_tr_pv += tr_pv[i];
        tranches.push(CmoTrancheValuation {
            pv: tr_pv[i],
            wal: tr_wal_num[i] / tr_prin_tot[i].max_val(T::from_f64(1e-12)),
            total_principal: tr_prin_tot[i],
            total_interest: tr_int_tot[i],
        });
    }
    CmoWaterfallValuation {
        pool_pv,
        pool_wal: pool_wal_num / s.original_balance.max_val(T::from_f64(1e-12)),
        tranches,
        io_strip_pv: io_pv,
        total_structure_pv: sum_tr_pv,
    }
}

#[derive(Clone, Debug)]
pub struct BondAssetSwapSpec<T: Real = f64> {
    pub clean_price: T,
    pub accrued_interest: T,
    pub coupon_rate: T,
    pub payment_times: Vec<f64>,
    pub accrual_fractions: Vec<f64>,
}

#[derive(Copy, Clone, Debug)]
pub struct BondSpreadReport<T: Real> {
    pub dirty_price: T,
    pub theoretical_dirty_price: T,
    pub float_annuity: T,
    pub par_asset_swap_spread: T,
    pub market_value_asset_swap_spread: T,
    pub z_spread: T,
    pub spread_dv01: T,
}

pub fn analyze_bond_asset_swap_and_z_spread<T: Real>(
    spec: &BondAssetSwapSpec<T>,
    ois_curve: &impl YieldTermStructure<T>,
) -> BondSpreadReport<T> {
    let n = spec.payment_times.len().min(spec.accrual_fractions.len());
    if n == 0 {
        return BondSpreadReport {
            dirty_price: T::zero(),
            theoretical_dirty_price: T::zero(),
            float_annuity: T::zero(),
            par_asset_swap_spread: T::zero(),
            market_value_asset_swap_spread: T::zero(),
            z_spread: T::zero(),
            spread_dv01: T::zero(),
        };
    }
    let par_scale = if spec.clean_price.val() > 2.0 {
        100.0
    } else {
        1.0
    };
    let dirty = (spec.clean_price + spec.accrued_interest) / par_scale;
    let mut ann = T::zero();
    let mut theo_dirty = T::zero();
    let mut dfs = Vec::with_capacity(n);
    let mut cfs = Vec::with_capacity(n);
    let mut dfs_f64 = Vec::with_capacity(n);
    let mut cfs_f64 = Vec::with_capacity(n);
    for i in 0..n {
        let t = spec.payment_times[i].max(1e-8);
        let tau = spec.accrual_fractions[i].max(0.0);
        let df = ois_curve.discount_factor(t);
        ann += df * tau;
        let cf = if i + 1 == n {
            spec.coupon_rate * tau + T::one()
        } else {
            spec.coupon_rate * tau
        };
        theo_dirty += cf * df;
        dfs.push(df);
        cfs.push(cf);
        dfs_f64.push(df.val());
        cfs_f64.push(cf.val());
    }
    let ann_safe = ann.max_val(T::from_f64(1e-12));
    let dirty_safe = dirty.max_val(T::from_f64(1e-12));
    let par_asw = (theo_dirty - dirty) / ann_safe;
    let mv_asw = (theo_dirty - dirty) / (dirty_safe * ann_safe);

    let dirty_f64 = dirty.val();
    let eval_f64 = |z: f64| -> (f64, f64, f64) {
        let mut r = -dirty_f64;
        let mut dr = 0.0;
        let mut d2r = 0.0;
        for i in 0..n {
            let t = spec.payment_times[i].max(1e-8);
            let term = cfs_f64[i] * dfs_f64[i] * (-z * t).exp();
            r += term;
            dr -= t * term;
            d2r += t * t * term;
        }
        (r, dr, d2r)
    };
    let mut lo = -0.10f64;
    let mut hi = 0.50f64;
    let f_lo = eval_f64(lo).0;
    for _ in 0..52 {
        let mid = 0.5 * (lo + hi);
        let f_mid = eval_f64(mid).0;
        if (f_mid >= 0.0) == (f_lo >= 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mut z_root = 0.5 * (lo + hi);
    for _ in 0..4 {
        let (r, dr, _) = eval_f64(z_root);
        if dr.abs() > 1e-14 {
            z_root = (z_root - r / dr).clamp(-0.20, 1.0);
        }
    }
    let (_, _, f_zz_f64) = eval_f64(z_root);
    let mut f_z = T::zero();
    for i in 0..n {
        let t = spec.payment_times[i].max(1e-8);
        f_z -= cfs[i] * dfs[i] * (-z_root * t).exp() * t;
    }
    let residual_at = |z_t: T| -> T {
        let mut pv = T::zero();
        for i in 0..n {
            let t = spec.payment_times[i].max(1e-8);
            pv += cfs[i] * dfs[i] * (-z_t * t).exp();
        }
        pv - dirty
    };
    let z_spread = ql_ad::ift_lift_1d(z_root, residual_at, f_z, f_zz_f64);
    let mut dv01 = T::zero();
    for i in 0..n {
        let t = spec.payment_times[i].max(1e-8);
        dv01 += cfs[i] * dfs[i] * (-z_spread * t).exp() * t;
    }
    BondSpreadReport {
        dirty_price: dirty,
        theoretical_dirty_price: theo_dirty,
        float_annuity: ann,
        par_asset_swap_spread: par_asw,
        market_value_asset_swap_spread: mv_asw,
        z_spread,
        spread_dv01: dv01,
    }
}

#[cfg(test)]
mod phase17_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;

    #[test]
    fn bond_asw_and_z_spread_identities_and_ift_ad() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let ois = FlatForwardCurve::new(d, 0.03);
        let times = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let taus = vec![1.0; 5];
        let mut theo = 0.0;
        for &t in &times {
            theo += 0.04 * ois.discount_factor(t);
        }
        theo += ois.discount_factor(5.0);
        let at_theo = analyze_bond_asset_swap_and_z_spread(
            &BondAssetSwapSpec {
                clean_price: theo,
                accrued_interest: 0.0,
                coupon_rate: 0.04,
                payment_times: times.clone(),
                accrual_fractions: taus.clone(),
            },
            &ois,
        );
        assert!(at_theo.par_asset_swap_spread.abs() < 1e-12);
        assert!(at_theo.market_value_asset_swap_spread.abs() < 1e-12);
        assert!(at_theo.z_spread.abs() < 1e-12);

        let at_par = analyze_bond_asset_swap_and_z_spread(
            &BondAssetSwapSpec {
                clean_price: 1.0,
                accrued_interest: 0.0,
                coupon_rate: 0.04,
                payment_times: times.clone(),
                accrual_fractions: taus.clone(),
            },
            &ois,
        );
        assert!(
            (at_par.par_asset_swap_spread - at_par.market_value_asset_swap_spread).abs() < 1e-12
        );

        let h = 1e-5;
        let eval_z = |p: f64| {
            analyze_bond_asset_swap_and_z_spread(
                &BondAssetSwapSpec {
                    clean_price: p,
                    accrued_interest: 0.01,
                    coupon_rate: 0.04,
                    payment_times: times.clone(),
                    accrual_fractions: taus.clone(),
                },
                &ois,
            )
            .z_spread
        };
        let fd = (eval_z(0.96 + h) - eval_z(0.96 - h)) / (2.0 * h);
        type D = Dual<1>;
        let ois_d = FlatForwardCurve::new(d, D::from_value(0.03));
        let rep_d = analyze_bond_asset_swap_and_z_spread(
            &BondAssetSwapSpec {
                clean_price: D::variable(0.96, 0),
                accrued_interest: D::from_value(0.01),
                coupon_rate: D::from_value(0.04),
                payment_times: times,
                accrual_fractions: taus,
            },
            &ois_d,
        );
        assert!((rep_d.z_spread.d[0] - fd).abs() < 1e-6);
        assert!((rep_d.z_spread.d[0] + 1.0 / rep_d.spread_dv01.val).abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase19_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;

    #[test]
    fn mbs_burnout_oas_zero_recovery_burnout_monotonicity_and_ift_ad() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let curve = FlatForwardCurve::new(d, 0.035);
        let mut spec = MbsBurnoutOasSpec {
            market_price: 100.0,
            original_balance: 100.0,
            wac: 0.060,
            pass_through_rate: 0.055,
            term_months: 60,
            seasoning_months: 12,
            turnover_cpr: 0.06,
            max_refi_cpr: 0.35,
            refi_slope: 80.0,
            refi_incentive_threshold: 0.005,
            burnout_beta: 12.0,
            primary_secondary_spread: 0.010,
        };
        let (p0, wal_with_burnout, b_term) = price_mbs_with_burnout(&spec, &curve, 0.0);
        spec.market_price = p0;
        let rep0 = analyze_mbs_burnout_oas_ad(&spec, &curve);
        assert!(rep0.oas.abs() < 1e-9);
        assert!(b_term < 1.0 && b_term > 0.0);
        assert!(rep0.option_adjusted_duration > 0.0);

        let mut no_burn = spec.clone();
        no_burn.burnout_beta = 0.0;
        let (_, wal_no_burnout, b_no) = price_mbs_with_burnout(&no_burn, &curve, 0.0);
        assert!((b_no - 1.0).abs() < 1e-12);
        // Burnout slows in-the-money refinancing over time, extending WAL
        assert!(wal_with_burnout > wal_no_burnout);

        let h = 1e-5;
        let eval_oas = |mp: f64| {
            let mut s = spec.clone();
            s.market_price = mp;
            analyze_mbs_burnout_oas_ad(&s, &curve).oas
        };
        let fd = (eval_oas(101.5 + h) - eval_oas(101.5 - h)) / (2.0 * h);
        type D = Dual<1>;
        let curve_d = FlatForwardCurve::new(d, D::from_value(0.035));
        let spec_d = MbsBurnoutOasSpec {
            market_price: D::variable(101.5, 0),
            original_balance: D::from_value(100.0),
            wac: D::from_value(0.060),
            pass_through_rate: D::from_value(0.055),
            term_months: 60,
            seasoning_months: 12,
            turnover_cpr: D::from_value(0.06),
            max_refi_cpr: D::from_value(0.35),
            refi_slope: D::from_value(80.0),
            refi_incentive_threshold: D::from_value(0.005),
            burnout_beta: D::from_value(12.0),
            primary_secondary_spread: D::from_value(0.010),
        };
        let ad = analyze_mbs_burnout_oas_ad(&spec_d, &curve_d).oas;
        assert!(ad.d[0] < 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase24_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(r: T) -> DavidsonLevinIoPoStripMbsSpec<T> {
        DavidsonLevinIoPoStripMbsSpec {
            notional: T::from_f64(100.0),
            wac_annual: T::from_f64(0.065),
            pass_through_coupon_annual: T::from_f64(0.06),
            current_mortgage_rate: r,
            discount_rate_annual: T::from_f64(0.04),
            remaining_months: 24,
            loan_age_months: 12,
            turnover_base_cpr: T::from_f64(0.04),
            max_refi_cpr: T::from_f64(0.45),
            refi_slope: T::from_f64(80.0),
            refi_threshold_spread: T::from_f64(0.005),
            burnout_exponent: T::from_f64(1.0),
            rate_bump_bps: T::from_f64(1.0),
        }
    }
    #[test]
    fn phase24_mbs_invariants_and_ad() {
        let s = mk(0.055);
        let r = analyze_davidson_levin_mbs_io_po_strips(&s);
        assert!(r.strip_replication_error.abs() < 1e-11);
        assert!(r.io_effective_duration < 0.);
        assert!(
            r.po_effective_duration > r.pass_through_effective_duration
                && r.pass_through_effective_duration > 0.
        );
        let h = 1e-5;
        let fd = (analyze_davidson_levin_mbs_io_po_strips(&mk(0.055 + h)).io_strip_price
            - analyze_davidson_levin_mbs_io_po_strips(&mk(0.055 - h)).io_strip_price)
            / (2.0 * h);
        let ad = analyze_davidson_levin_mbs_io_po_strips(&mk(Dual::<1>::variable(0.055, 0)))
            .io_strip_price;
        assert!((ad.d[0] - fd).abs() < 1e-5);
        assert!((r.io_strip_price + r.po_strip_price - r.pass_through_price).abs() < 1e-10);
        let fd_po = (analyze_davidson_levin_mbs_io_po_strips(&mk(0.055 + h)).po_strip_price
            - analyze_davidson_levin_mbs_io_po_strips(&mk(0.055 - h)).po_strip_price)
            / (2.0 * h);
        let ad_po = analyze_davidson_levin_mbs_io_po_strips(&mk(Dual::<1>::variable(0.055, 0)))
            .po_strip_price;
        assert!((ad_po.d[0] - fd_po).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase26_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(p: T) -> TuckmanBurghardtTreasuryWildcardBasisSpec<T> {
        TuckmanBurghardtTreasuryWildcardBasisSpec {
            futures_settlement_price: p,
            delivery_horizon_years: T::from_f64(0.25),
            general_collateral_repo_rate: T::from_f64(0.03),
            repo_specialness_spread: T::from_f64(0.002),
            bond1_clean_price: T::from_f64(95.0),
            bond1_conversion_factor: T::from_f64(0.85),
            bond1_modified_duration: T::from_f64(6.0),
            bond1_yield_beta: T::one(),
            bond2_clean_price: T::from_f64(101.0),
            bond2_conversion_factor: T::from_f64(0.9),
            bond2_modified_duration: T::from_f64(7.0),
            bond2_yield_beta: T::from_f64(1.1),
            afternoon_wildcard_vol_annual: T::from_f64(0.1),
            eom_delivery_window_days: 5,
        }
    }
    #[test]
    fn wildcard_invariants_and_dual_fd() {
        let r = evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch(&spec(100.0));
        assert!(
            r.quality_switch_option_value > 0.0
                && r.afternoon_wildcard_option_value > 0.0
                && r.eom_timing_carry_option_value > 0.0
        );
        assert!(r.special_repo_carry_pickup > 0.0);
        let h = 1e-5;
        let fd = (evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch(&spec(100.0 + h))
            .theoretical_net_basis_ctd
            - evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch(&spec(100.0 - h))
                .theoretical_net_basis_ctd)
            / (2.0 * h);
        let ad = evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch(&spec(
            Dual::<1>::variable(100.0, 0),
        ))
        .theoretical_net_basis_ctd;
        assert!(ad.d[0].abs() > 1e-4 && (ad.d[0] - fd).abs() < 1e-5);
        let mut up = spec(100.0);
        let mut dn = spec(100.0);
        up.bond1_clean_price += h;
        dn.bond1_clean_price -= h;
        let fd_q = (evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch(&up)
            .quality_switch_option_value
            - evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch(&dn)
                .quality_switch_option_value)
            / (2.0 * h);
        let mut sd = spec(Dual::<1>::from_f64(100.0));
        sd.bond1_clean_price = Dual::<1>::variable(95.0, 0);
        let ad_q = evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch(&sd)
            .quality_switch_option_value;
        println!(
            "q fd={} ad={} base={}",
            fd_q, ad_q.d[0], r.quality_switch_option_value
        );
        assert!(ad_q.d[0].abs() > 1e-4 && (ad_q.d[0] - fd_q).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase29_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn veronesi_invariants_and_dual() {
        fn mk<T: Real>(x: T) -> VeronesiFabozziConvexityPrepaymentSpec<T> {
            VeronesiFabozziConvexityPrepaymentSpec {
                forward_rate: x,
                short_rate_vol: T::from_f64(0.2),
                mean_reversion_kappa: T::from_f64(0.1),
                fixing_time_t1: T::one(),
                maturity_time_t2: T::from_f64(2.0),
                wac_coupon: T::from_f64(0.05),
                current_mortgage_rate: T::from_f64(0.04),
                min_cpr_turnover: T::from_f64(0.02),
                max_cpr_refi: T::from_f64(0.3),
                refi_steepness_gamma: T::from_f64(8.0),
                refi_incentive_threshold: T::zero(),
                pool_age_months: T::from_f64(24.0),
                psa_ramp_months: T::from_f64(30.0),
            }
        }
        let r = evaluate_veronesi_fabozzi_futures_convexity_and_prepayment(&mk(0.03));
        assert!(r.futures_forward_convexity_bias_bps > 0.0);
        let h = 1e-5;
        let fd = (evaluate_veronesi_fabozzi_futures_convexity_and_prepayment(&mk(0.03 + h))
            .implied_futures_rate
            - evaluate_veronesi_fabozzi_futures_convexity_and_prepayment(&mk(0.03 - h))
                .implied_futures_rate)
            / (2.0 * h);
        let ad = evaluate_veronesi_fabozzi_futures_convexity_and_prepayment(&mk(
            Dual::<1>::variable(0.03, 0),
        ))
        .implied_futures_rate;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        let mut zero = mk(0.03);
        zero.short_rate_vol = 0.0;
        let rz = evaluate_veronesi_fabozzi_futures_convexity_and_prepayment(&zero);
        assert!(rz.futures_forward_convexity_bias_bps.abs() < 1e-12);
        assert!(
            (1.0 - (1.0 - rz.single_monthly_mortality_smm).powi(12) - rz.scurve_cpr).abs() < 1e-12
        );
    }
}

#[cfg(test)]
mod phase15_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;
    #[test]
    fn specified_pool_and_tba_invariants() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let c = FlatForwardCurve::new(d, 0.03);
        let s = TbaDollarRollSpec {
            front_clean_price: 99.0,
            back_clean_price: 98.5,
            wac: 0.04,
            pass_through_rate: 0.035,
            remaining_months: 120,
            front_accrued_fraction: 0.5,
            back_accrued_fraction: 0.5,
            roll_period_years: 1.0 / 12.0,
            assumed_cpr: 0.06,
            repo_rate: 0.02,
        };
        let r = analyze_tba_dollar_roll(&s);
        let mut s2 = s;
        s2.back_clean_price = 99.0 - r.breakeven_drop;
        let r2 = analyze_tba_dollar_roll(&s2);
        assert!(r2.net_roll_advantage.abs() < 1e-10);
        assert!((r2.implied_repo_rate - s.repo_rate).abs() < 1e-10);
        let h = 1e-5;
        let mut s_up = s;
        let mut s_dn = s;
        s_up.front_clean_price += h;
        s_dn.front_clean_price -= h;
        let fd = (analyze_tba_dollar_roll(&s_up).implied_repo_rate
            - analyze_tba_dollar_roll(&s_dn).implied_repo_rate)
            / (2.0 * h);
        let sd = TbaDollarRollSpec::<Dual<1>> {
            front_clean_price: Dual::variable(99.0, 0),
            back_clean_price: Dual::from_value(98.5),
            wac: Dual::from_value(0.04),
            pass_through_rate: Dual::from_value(0.035),
            remaining_months: 120,
            front_accrued_fraction: 0.5,
            back_accrued_fraction: 0.5,
            roll_period_years: 1.0 / 12.0,
            assumed_cpr: Dual::from_value(0.06),
            repo_rate: Dual::from_value(0.02),
        };
        let ad = analyze_tba_dollar_roll(&sd).implied_repo_rate;
        assert!((ad.d[0] - fd).abs() < 1e-6);
        let m = SpecifiedPoolMbsSpec {
            original_balance: 100.0,
            wac: 0.04,
            pass_through_rate: 0.035,
            term_months: 12,
            seasoning_months: 30,
            cpr_min: 0.06,
            cpr_max: 0.06,
            refi_slope: 0.0,
            refi_inflection: 0.0,
            primary_secondary_spread: 0.0,
            cdr: 0.01,
            loss_severity: 0.4,
        };
        let val = price_specified_pool_mbs(&m, &c, 0.0);
        assert!(val.pv > 90.0 && val.wal > 0.4 && val.wal <= 1.0);
    }
}

#[cfg(test)]
mod phase16_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;
    #[test]
    fn treasury_ctd_and_cmo_waterfall_conservation() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let c = FlatForwardCurve::new(d, 0.035);
        let basket = TreasuryFuturesBasketSpec {
            futures_price: 98.5,
            delivery_time: 0.25,
            repo_rate: 0.03,
            notional_standard_yield: 0.06,
            deliverables: vec![
                TreasuryFutureDeliverable {
                    clean_price: 99.2,
                    coupon_rate: 0.06,
                    remaining_semiannual_periods: 16,
                    accrued_interest_spot: 0.5,
                    accrued_interest_delivery: 1.0,
                    interim_coupon: 0.0,
                    yield_vol: 0.008,
                    modified_duration: 6.5,
                },
                TreasuryFutureDeliverable {
                    clean_price: 89.0,
                    coupon_rate: 0.045,
                    remaining_semiannual_periods: 18,
                    accrued_interest_spot: 0.4,
                    accrued_interest_delivery: 0.8,
                    interim_coupon: 0.0,
                    yield_vol: 0.008,
                    modified_duration: 7.2,
                },
            ],
        };
        let ctd = analyze_treasury_futures_ctd(&basket);
        assert!((ctd.bond_metrics[0].conversion_factor - 1.0).abs() < 1e-12);
        let dirty0 =
            basket.deliverables[0].clean_price + basket.deliverables[0].accrued_interest_spot;
        let bnoc0 = dirty0
            * basket.delivery_time
            * (basket.repo_rate - ctd.bond_metrics[0].implied_repo_rate);
        assert!((ctd.bond_metrics[0].net_basis - bnoc0).abs() < 1e-12);
        assert!(ctd.switch_option_value > 0.0);
        let cmo = CmoWaterfallSpec {
            original_balance: 100.0,
            wac: 0.055,
            pass_through_rate: 0.05,
            term_months: 60,
            psa_speed: 150.0,
            pac_collar_psa: Some((100.0, 300.0)),
            tranches: vec![
                CmoTrancheSpec {
                    initial_balance: 60.0,
                    coupon_rate: 0.045,
                    is_pac: true,
                },
                CmoTrancheSpec {
                    initial_balance: 40.0,
                    coupon_rate: 0.045,
                    is_pac: false,
                },
            ],
        };
        let w1 = price_cmo_waterfall(&cmo, &c, 0.0);
        let mut cmo2 = cmo.clone();
        cmo2.psa_speed = 250.0;
        let w2 = price_cmo_waterfall(&cmo2, &c, 0.0);
        assert!((w1.total_structure_pv - w1.pool_pv).abs() < 1e-10);
        assert!((w1.tranches[0].wal - w2.tranches[0].wal).abs() < 1e-10);
        let h = 1e-5;
        let mut b_up = basket.clone();
        let mut b_dn = basket.clone();
        b_up.futures_price += h;
        b_dn.futures_price -= h;
        let fd = (analyze_treasury_futures_ctd(&b_up).ctd_implied_repo_rate
            - analyze_treasury_futures_ctd(&b_dn).ctd_implied_repo_rate)
            / (2.0 * h);
        type D = Dual<1>;
        let bd = TreasuryFuturesBasketSpec {
            futures_price: D::variable(98.5, 0),
            delivery_time: 0.25,
            repo_rate: D::from_value(0.03),
            notional_standard_yield: 0.06,
            deliverables: basket
                .deliverables
                .iter()
                .map(|x| TreasuryFutureDeliverable {
                    clean_price: D::from_value(x.clean_price),
                    coupon_rate: D::from_value(x.coupon_rate),
                    remaining_semiannual_periods: x.remaining_semiannual_periods,
                    accrued_interest_spot: D::from_value(x.accrued_interest_spot),
                    accrued_interest_delivery: D::from_value(x.accrued_interest_delivery),
                    interim_coupon: D::from_value(x.interim_coupon),
                    yield_vol: D::from_value(x.yield_vol),
                    modified_duration: D::from_value(x.modified_duration),
                })
                .collect(),
        };
        let ad = analyze_treasury_futures_ctd(&bd).ctd_implied_repo_rate;
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}
pub use crate::phase34::*;
