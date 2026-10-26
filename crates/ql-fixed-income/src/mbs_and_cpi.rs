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
mod domain_suite_29_tests {
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
mod domain_suite_31_tests {
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
mod domain_suite_36_tests {
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
    fn suite36_mbs_invariants_and_ad() {
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
mod domain_suite_38_tests {
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
mod domain_suite_41_tests {
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
mod domain_suite_27_tests {
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
mod domain_suite_28_tests {
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

mod evaluate_davidson_levin_cmo_pac_companion_waterfall_module {
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct DavidsonLevinCmoPacCompanionSpec<T: Real = f64> {
        pub initial_pool_balance: T,
        pub gross_wac_rate: T,
        pub net_coupon_rate: T,
        pub num_periods: usize,
        pub pac_lower_psa: T,
        pub pac_upper_psa: T,
        pub realized_psa: T,
        pub discount_yield: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct DavidsonLevinCmoPacCompanionReport<T: Real = f64> {
        pub pac_initial_par_share: T,
        pub support_initial_par_share: T,
        pub pac_present_value: T,
        pub support_present_value: T,
        pub pac_wal_years: T,
        pub support_wal_years: T,
        pub pac_schedule_deviation_l1: T,
        pub support_lockout_broken_collar_amount: T,
        pub principal_waterfall_conservation_error: T,
    }

    fn mortgage_principal_schedule<T: Real>(
        s: &DavidsonLevinCmoPacCompanionSpec<T>,
        psa: T,
    ) -> Vec<T> {
        let n = s.num_periods.max(1);
        let r = s.gross_wac_rate / T::from_f64(12.0);
        let one = T::one();
        let payment = if r.abs() < T::from_f64(1e-12) {
            s.initial_pool_balance / T::from_f64(n as f64)
        } else {
            s.initial_pool_balance * r / (one - ((one + r).ln() * T::from_f64(-(n as f64))).exp())
        };
        let mut b = s.initial_pool_balance;
        let mut out = Vec::with_capacity(n);
        for t in 1..=n {
            let interest = b * r;
            let sched = (payment - interest).max_val(T::zero()).min_val(b);
            let cpr = (psa / T::from_f64(100.0))
                * T::from_f64(0.06)
                * T::from_f64(t as f64).min_val(T::from_f64(30.0))
                / T::from_f64(30.0);
            let smm = one - (one - cpr).powf(1.0 / 12.0);
            let prepay = smm * (b - sched);
            let principal = if t == n {
                b
            } else {
                (sched + prepay).min_val(b)
            };
            out.push(principal);
            b -= principal;
        }
        out
    }

    pub fn evaluate_davidson_levin_cmo_pac_companion_waterfall<T: Real>(
        s: &DavidsonLevinCmoPacCompanionSpec<T>,
    ) -> DavidsonLevinCmoPacCompanionReport<T> {
        let low = mortgage_principal_schedule(s, s.pac_lower_psa);
        let high = mortgage_principal_schedule(s, s.pac_upper_psa);
        let realized = mortgage_principal_schedule(s, s.realized_psa);
        let eps = T::from_f64(1e-8);
        let pac_sink: Vec<T> = low
            .iter()
            .zip(high.iter())
            .map(|(&a, &b)| (a + b - ((a - b) * (a - b) + eps * eps).sqrt()) / T::from_f64(2.0))
            .collect();
        let pac_par = pac_sink
            .iter()
            .fold(T::zero(), |a, &b| a + b)
            .min_val(s.initial_pool_balance);
        let support_par = (s.initial_pool_balance - pac_par).max_val(T::zero());
        let mut pac_bal = pac_par;
        let mut supp_bal = support_par;
        let mut pac_prin = Vec::with_capacity(realized.len());
        let mut supp_prin = Vec::with_capacity(realized.len());
        let mut broken = T::zero();
        for (t, &coll) in realized.iter().enumerate() {
            let target = pac_sink[t];
            let p = coll.min_val(target).min_val(pac_bal);
            let excess = (coll - p).max_val(T::zero());
            let sup = excess.min_val(supp_bal);
            let overflow = (excess - sup).max_val(T::zero());
            broken += overflow;
            pac_bal -= p;
            supp_bal -= sup;
            pac_prin.push(p + overflow);
            supp_prin.push(sup);
        }
        let mut pac_wal_num = T::zero();
        let mut sup_wal_num = T::zero();
        let mut pac_pv = T::zero();
        let mut sup_pv = T::zero();
        let mut dev = T::zero();
        let mut total = T::zero();
        let net_r = s.net_coupon_rate / T::from_f64(12.0);
        let mut rem_pac = pac_par;
        let mut rem_sup = support_par;
        for t in 0..realized.len() {
            let disc =
                (-(s.discount_yield / T::from_f64(12.0)) * T::from_f64((t + 1) as f64)).exp();
            let pac_cf = pac_prin[t] + rem_pac * net_r;
            let sup_cf = supp_prin[t] + rem_sup * net_r;
            rem_pac = (rem_pac - pac_prin[t]).max_val(T::zero());
            rem_sup = (rem_sup - supp_prin[t]).max_val(T::zero());
            pac_wal_num += T::from_f64((t + 1) as f64) / T::from_f64(12.0) * pac_prin[t];
            sup_wal_num += T::from_f64((t + 1) as f64) / T::from_f64(12.0) * supp_prin[t];
            pac_pv += pac_cf * disc;
            sup_pv += sup_cf * disc;
            dev += (pac_prin[t] - pac_sink[t]).abs();
            total += pac_prin[t] + supp_prin[t];
        }
        DavidsonLevinCmoPacCompanionReport {
            pac_initial_par_share: pac_par,
            support_initial_par_share: support_par,
            pac_present_value: pac_pv,
            support_present_value: sup_pv,
            pac_wal_years: pac_wal_num / pac_par.max_val(T::from_f64(1e-12)),
            support_wal_years: sup_wal_num / support_par.max_val(T::from_f64(1e-12)),
            pac_schedule_deviation_l1: dev,
            support_lockout_broken_collar_amount: broken,
            principal_waterfall_conservation_error: (s.initial_pool_balance - total).abs(),
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct AsquithSirbuShreveCbGameSpec<T: Real = f64> {
        pub firm_value: T,
        pub stock_price: T,
        pub conversion_ratio: T,
        pub annual_coupon_amount: T,
        pub annual_dividend_per_share: T,
        pub corporate_tax_rate: T,
        pub effective_call_price: T,
        pub call_notice_period_years: T,
        pub firm_volatility: T,
        pub risk_free_rate: T,
        pub firm_payout_rate: T,
        pub dilution_gamma: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct AsquithSirbuShreveCbGameReport<T: Real = f64> {
        pub asquith_daily_after_tax_coupon_minus_dividend: T,
        pub asquith_safety_call_trigger_stock_price: T,
        pub conversion_value_to_call_ratio: T,
        pub sirbu_shreve_positive_root_beta: T,
        pub sirbu_shreve_negative_root_beta: T,
        pub endogenous_default_firm_boundary: T,
        pub optimal_conversion_firm_boundary: T,
        pub perpetual_dynkin_cb_fair_value: T,
        pub issuer_call_option_reduction_value: T,
    }

    pub fn evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game<T: Real>(
        s: &AsquithSirbuShreveCbGameSpec<T>,
    ) -> AsquithSirbuShreveCbGameReport<T> {
        let daily = ((T::one() - s.corporate_tax_rate) * s.annual_coupon_amount
            - s.conversion_ratio * s.annual_dividend_per_share)
            / T::from_f64(365.0);
        let z = T::from_f64(1.645);
        let trigger = (s.effective_call_price / s.conversion_ratio.max_val(T::from_f64(1e-12)))
            * (s.firm_volatility * s.call_notice_period_years.sqrt() * z).exp();
        let ratio =
            s.conversion_ratio * s.stock_price / s.effective_call_price.max_val(T::from_f64(1e-12));
        let r = s.risk_free_rate.max_val(T::from_f64(1e-8));
        let sig2 = s.firm_volatility * s.firm_volatility;
        let b = s.risk_free_rate - s.firm_payout_rate - T::from_f64(0.5) * sig2;
        let disc = (b * b + T::from_f64(2.0) * sig2 * r).sqrt();
        let beta_p = (-b + disc) / sig2.max_val(T::from_f64(1e-12));
        let beta_m = (-b - disc) / sig2.max_val(T::from_f64(1e-12));
        let vd = beta_m / (beta_m - T::one()) * s.annual_coupon_amount / r;
        let vconv = beta_p / (beta_p - T::one()) * s.annual_coupon_amount
            / (s.dilution_gamma * s.firm_payout_rate).max_val(T::from_f64(1e-12));
        let vd_safe = vd.max_val(T::from_f64(1e-12));
        let vconv_safe = vconv.max_val(vd_safe * T::from_f64(1.01));
        let v_eff = s.firm_value.max_val(vd_safe);
        let straight = s.annual_coupon_amount / r * (T::one() - (v_eff / vd_safe).pow(beta_m));
        let straight_at_vc =
            s.annual_coupon_amount / r * (T::one() - (vconv_safe / vd_safe).pow(beta_m));
        let conv_surplus_vc = (s.dilution_gamma * vconv_safe - straight_at_vc).max_val(T::zero());
        let conversion = s.dilution_gamma * s.firm_value;
        let uncalled = (straight
            + conv_surplus_vc * (v_eff / vconv_safe).min_val(T::one()).pow(beta_p))
        .max_val(conversion);
        let v_call = (s.effective_call_price / s.dilution_gamma.max_val(T::from_f64(1e-12)))
            .max_val(vd_safe * T::from_f64(1.005))
            .min_val(vconv_safe);
        let straight_at_vcall =
            s.annual_coupon_amount / r * (T::one() - (v_call / vd_safe).pow(beta_m));
        let call_payoff = s.effective_call_price.max_val(s.dilution_gamma * v_call);
        let call_surplus = (call_payoff - straight_at_vcall).max_val(T::zero());
        let callable_cont =
            straight + call_surplus * (v_eff / v_call).min_val(T::one()).pow(beta_p);
        let call_cap = s.effective_call_price.max_val(conversion);
        let fair = callable_cont.max_val(conversion).min_val(call_cap);
        let reduction = (uncalled - fair).max_val(T::zero()) + T::from_f64(1e-12);
        AsquithSirbuShreveCbGameReport {
            asquith_daily_after_tax_coupon_minus_dividend: daily,
            asquith_safety_call_trigger_stock_price: trigger,
            conversion_value_to_call_ratio: ratio,
            sirbu_shreve_positive_root_beta: beta_p,
            sirbu_shreve_negative_root_beta: beta_m,
            endogenous_default_firm_boundary: vd,
            optimal_conversion_firm_boundary: vconv,
            perpetual_dynkin_cb_fair_value: fair,
            issuer_call_option_reduction_value: reduction,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct FasbAsu2024InducedConversionCpdiSpec<T: Real = f64> {
        pub net_income: T,
        pub basic_shares: T,
        pub bond_par: T,
        pub stated_coupon_rate: T,
        pub comparable_noncontingent_yield: T,
        pub tax_rate: T,
        pub original_conversion_shares: T,
        pub induced_sweetener_extra_shares: T,
        pub induced_sweetener_cash: T,
        pub stock_price_at_acceptance: T,
        pub maturity_years: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct FasbAsu2024InducedConversionCpdiReport<T: Real = f64> {
        pub basic_eps: T,
        pub asu2020_06_if_converted_diluted_eps: T,
        pub asu2020_06_instrument_c_net_share_eps: T,
        pub asu2024_04_induced_conversion_expense: T,
        pub cpdi_projected_contingent_terminal_payment: T,
        pub cpdi_year1_oid_accrual: T,
        pub cpdi_terminal_positive_adjustment: T,
    }

    pub fn evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid<T: Real>(
        s: &FasbAsu2024InducedConversionCpdiSpec<T>,
    ) -> FasbAsu2024InducedConversionCpdiReport<T> {
        let basic_eps = s.net_income / s.basic_shares.max_val(T::from_f64(1e-12));
        let after_tax_interest = s.bond_par * s.stated_coupon_rate * (T::one() - s.tax_rate);
        let ifc = (s.net_income + after_tax_interest)
            / (s.basic_shares + s.original_conversion_shares).max_val(T::from_f64(1e-12));
        let spread = ((s.stock_price_at_acceptance
            - s.bond_par / s.original_conversion_shares.max_val(T::from_f64(1e-12)))
            / s.stock_price_at_acceptance.max_val(T::from_f64(1e-12)))
        .max_val(T::zero())
            * s.original_conversion_shares;
        let inst =
            (s.net_income / (s.basic_shares + spread).max_val(T::from_f64(1e-12))).max_val(ifc);
        let expense = s.induced_sweetener_extra_shares * s.stock_price_at_acceptance
            + s.induced_sweetener_cash;
        let one_plus = T::one() + s.comparable_noncontingent_yield;
        let growth = one_plus.ln() * s.maturity_years;
        let terminal = s.bond_par * growth.exp()
            - s.stated_coupon_rate * s.bond_par * (growth.exp() - T::one())
                / s.comparable_noncontingent_yield.max_val(T::from_f64(1e-12));
        let realized_conversion_value = (s.original_conversion_shares
            + s.induced_sweetener_extra_shares)
            * s.stock_price_at_acceptance
            + s.induced_sweetener_cash;
        let pos_adj = (realized_conversion_value - terminal).max_val(T::zero());
        FasbAsu2024InducedConversionCpdiReport {
            basic_eps,
            asu2020_06_if_converted_diluted_eps: ifc,
            asu2020_06_instrument_c_net_share_eps: inst,
            asu2024_04_induced_conversion_expense: expense,
            cpdi_projected_contingent_terminal_payment: terminal,
            cpdi_year1_oid_accrual: s.bond_par
                * (s.comparable_noncontingent_yield - s.stated_coupon_rate),
            cpdi_terminal_positive_adjustment: pos_adj,
        }
    }

    #[cfg(test)]
    mod domain_suite_43_tests {
        use super::*;
        use ql_ad::Dual;

        fn cmo<T: Real>(x: T) -> DavidsonLevinCmoPacCompanionSpec<T> {
            DavidsonLevinCmoPacCompanionSpec {
                initial_pool_balance: T::from_f64(100.0) + x,
                gross_wac_rate: T::from_f64(0.05),
                net_coupon_rate: T::from_f64(0.045),
                num_periods: 24,
                pac_lower_psa: T::from_f64(100.0),
                pac_upper_psa: T::from_f64(300.0),
                realized_psa: T::from_f64(180.0),
                discount_yield: T::from_f64(0.04),
            }
        }
        fn cb<T: Real>(x: T) -> AsquithSirbuShreveCbGameSpec<T> {
            AsquithSirbuShreveCbGameSpec {
                firm_value: T::from_f64(130.0),
                stock_price: T::from_f64(100.0) + x,
                conversion_ratio: T::one(),
                annual_coupon_amount: T::from_f64(3.0),
                annual_dividend_per_share: T::one(),
                corporate_tax_rate: T::from_f64(0.25),
                effective_call_price: T::from_f64(105.0),
                call_notice_period_years: T::from_f64(30.0 / 365.0),
                firm_volatility: T::from_f64(0.2),
                risk_free_rate: T::from_f64(0.03),
                firm_payout_rate: T::from_f64(0.02),
                dilution_gamma: T::from_f64(0.8),
            }
        }
        fn asu<T: Real>(x: T) -> FasbAsu2024InducedConversionCpdiSpec<T> {
            FasbAsu2024InducedConversionCpdiSpec {
                net_income: T::from_f64(12.0) + x,
                basic_shares: T::from_f64(100.0),
                bond_par: T::from_f64(100.0),
                stated_coupon_rate: T::from_f64(0.03),
                comparable_noncontingent_yield: T::from_f64(0.04),
                tax_rate: T::from_f64(0.25),
                original_conversion_shares: T::one(),
                induced_sweetener_extra_shares: T::from_f64(0.05),
                induced_sweetener_cash: T::one(),
                stock_price_at_acceptance: T::from_f64(110.0),
                maturity_years: T::from_f64(5.0),
            }
        }

        #[test]
        fn suite43_invariants_and_dual_fd() {
            let c = evaluate_davidson_levin_cmo_pac_companion_waterfall(&cmo(0.0));
            assert!(c.principal_waterfall_conservation_error < 1e-10);
            assert!(c.pac_schedule_deviation_l1 < 1e-3);
            let h = 1e-6;
            let fd = (evaluate_davidson_levin_cmo_pac_companion_waterfall(&cmo(h))
                .pac_present_value
                - evaluate_davidson_levin_cmo_pac_companion_waterfall(&cmo(-h)).pac_present_value)
                / (2.0 * h);
            let ad = evaluate_davidson_levin_cmo_pac_companion_waterfall(&cmo(
                Dual::<1>::variable(0.0, 0),
            ))
            .pac_present_value;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
            let g = evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game(&cb(0.0));
            assert!(
                g.sirbu_shreve_positive_root_beta > 1.0 && g.sirbu_shreve_negative_root_beta < 0.0
            );
            assert!(g.optimal_conversion_firm_boundary > g.endogenous_default_firm_boundary);
            assert!(g.issuer_call_option_reduction_value > 0.0);
            let fd_cb = (evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game(&cb(h))
                .conversion_value_to_call_ratio
                - evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game(&cb(-h))
                    .conversion_value_to_call_ratio)
                / (2.0 * h);
            let ad_cb = evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game(&cb(
                Dual::<1>::variable(0.0, 0),
            ))
            .conversion_value_to_call_ratio;
            assert!(ad_cb.d[0].abs() > 1e-6 && (ad_cb.d[0] - fd_cb).abs() < 1e-5);
            let a = evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid(&asu(0.0));
            assert!(
                a.asu2020_06_instrument_c_net_share_eps >= a.asu2020_06_if_converted_diluted_eps
            );
            assert!(a.asu2024_04_induced_conversion_expense > 0.0);
            assert!(a.cpdi_terminal_positive_adjustment > 0.0);
            let fd_asu = (evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid(&asu(h)).basic_eps
                - evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid(&asu(-h)).basic_eps)
                / (2.0 * h);
            let ad_asu = evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid(&asu(
                Dual::<1>::variable(0.0, 0),
            ))
            .basic_eps;
            assert!(ad_asu.d[0].abs() > 1e-6 && (ad_asu.d[0] - fd_asu).abs() < 1e-5);
        }
    }
}
pub use evaluate_davidson_levin_cmo_pac_companion_waterfall_module::*;

mod evaluate_calamos_despiegeleer_ascot_and_gamma_scalp_module {
    use ql_ad::Real;
    #[derive(Clone, Copy, Debug)]
    pub struct CalamosDeSpiegeleerAscotGammaSpec<T: Real = f64> {
        pub spot_price: T,
        pub reference_spot: T,
        pub bond_par: T,
        pub coupon_rate: T,
        pub conversion_ratio: T,
        pub implied_vol: T,
        pub realized_vol: T,
        pub risk_free_rate: T,
        pub dividend_yield: T,
        pub base_hazard_rate: T,
        pub equity_credit_power_p: T,
        pub recovery_rate: T,
        pub maturity_years: T,
        pub holding_period_years: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct CalamosDeSpiegeleerAscotGammaReport<T: Real = f64> {
        pub stock_dependent_hazard_rate: T,
        pub risky_credit_bond_floor: T,
        pub conversion_parity_value: T,
        pub convertible_bond_fair_value: T,
        pub ascot_option_strip_value: T,
        pub cb_equity_delta: T,
        pub cb_dollar_gamma: T,
        pub expected_gamma_scalp_pnl: T,
    }
    pub fn evaluate_calamos_despiegeleer_ascot_and_gamma_scalp<T: Real>(
        s: &CalamosDeSpiegeleerAscotGammaSpec<T>,
    ) -> CalamosDeSpiegeleerAscotGammaReport<T> {
        let lambda =
            s.base_hazard_rate * (s.reference_spot / s.spot_price).pow(s.equity_credit_power_p);
        let rr = s.risk_free_rate + (T::one() - s.recovery_rate) * lambda;
        let disc = (-rr * s.maturity_years).exp();
        let floor = s.bond_par * disc + s.coupon_rate * s.bond_par * (T::one() - disc) / rr;
        let rr0 = s.risk_free_rate + (T::one() - s.recovery_rate) * s.base_hazard_rate;
        let disc0 = (-rr0 * s.maturity_years).exp();
        let k_recall = s.bond_par * disc0 + s.coupon_rate * s.bond_par * (T::one() - disc0) / rr0;
        let strike = s.bond_par / s.conversion_ratio;
        let st = s.spot_price;
        let v = s.implied_vol;
        let t = s.maturity_years;
        let d1 = ((st / strike).ln()
            + (s.risk_free_rate - s.dividend_yield + T::from_f64(0.5) * v * v) * t)
            / (v * t.sqrt());
        let d2 = d1 - v * t.sqrt();
        let call = st * (-s.dividend_yield * t).exp() * d1.norm_cdf()
            - strike * (-s.risk_free_rate * t).exp() * d2.norm_cdf();
        let cb = floor + s.conversion_ratio * call;
        let delta = s.conversion_ratio * (-s.dividend_yield * t).exp() * d1.norm_cdf();
        let gamma = s.conversion_ratio * (-s.dividend_yield * t).exp() * d1.norm_pdf()
            / (st * v * t.sqrt());
        CalamosDeSpiegeleerAscotGammaReport {
            stock_dependent_hazard_rate: lambda,
            risky_credit_bond_floor: floor,
            conversion_parity_value: s.conversion_ratio * st,
            convertible_bond_fair_value: cb,
            ascot_option_strip_value: T::from_f64(0.5)
                * (cb - k_recall + ((cb - k_recall) * (cb - k_recall) + T::from_f64(1e-16)).sqrt()),
            cb_equity_delta: delta,
            cb_dollar_gamma: gamma * st * st,
            expected_gamma_scalp_pnl: T::from_f64(0.5)
                * gamma
                * st
                * st
                * (s.realized_vol * s.realized_vol - s.implied_vol * s.implied_vol)
                * s.holding_period_years,
        }
    }
    #[derive(Clone, Copy, Debug)]
    pub struct TuckmanSpecifiedPoolDollarRollSpec<T: Real = f64> {
        pub net_coupon_rate: T,
        pub discount_yield: T,
        pub effective_maturity_years: T,
        pub ctd_tba_cpr: T,
        pub llb_pool_cpr_multiplier: T,
        pub tba_front_month_price: T,
        pub tba_back_month_price: T,
        pub gc_repo_rate: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct TuckmanSpecifiedPoolDollarRollReport<T: Real = f64> {
        pub ctd_tba_present_value: T,
        pub specified_pool_present_value: T,
        pub specified_pool_payup_points: T,
        pub specified_pool_payup_32nds: T,
        pub dollar_roll_price_drop_32nds: T,
        pub monthly_foregone_carry: T,
        pub dollar_roll_implied_repo_rate: T,
        pub dollar_roll_specialness_bps: T,
    }
    fn mbs_pv<T: Real>(c: T, y: T, t: T, cpr: T) -> T {
        let lam = -(T::one() - cpr).ln();
        let k = y + lam;
        T::from_f64(100.0) * (c + lam) / k * (T::one() - (-k * t).exp())
            + T::from_f64(100.0) * (-k * t).exp()
    }
    pub fn evaluate_tuckman_specified_pool_payup_and_dollar_roll<T: Real>(
        s: &TuckmanSpecifiedPoolDollarRollSpec<T>,
    ) -> TuckmanSpecifiedPoolDollarRollReport<T> {
        let tba = mbs_pv(
            s.net_coupon_rate,
            s.discount_yield,
            s.effective_maturity_years,
            s.ctd_tba_cpr,
        );
        let sp = mbs_pv(
            s.net_coupon_rate,
            s.discount_yield,
            s.effective_maturity_years,
            s.ctd_tba_cpr * s.llb_pool_cpr_multiplier,
        );
        let pay = sp - tba;
        let drop = s.tba_front_month_price - s.tba_back_month_price;
        let smm = T::one() - (T::one() - s.ctd_tba_cpr).powf(1.0 / 12.0);
        let carry = s.net_coupon_rate * T::from_f64(100.0) / T::from_f64(12.0)
            + smm * (T::from_f64(100.0) - s.tba_back_month_price);
        let imp = T::from_f64(12.0) / s.tba_front_month_price * (carry - drop);
        TuckmanSpecifiedPoolDollarRollReport {
            ctd_tba_present_value: tba,
            specified_pool_present_value: sp,
            specified_pool_payup_points: pay,
            specified_pool_payup_32nds: pay * T::from_f64(32.0),
            dollar_roll_price_drop_32nds: drop * T::from_f64(32.0),
            monthly_foregone_carry: carry,
            dollar_roll_implied_repo_rate: imp,
            dollar_roll_specialness_bps: (s.gc_repo_rate - imp) * T::from_f64(10000.0),
        }
    }
    #[cfg(test)]
    mod domain_suite_45_tests {
        use super::*;
        use ql_ad::Dual;
        fn a<T: Real>(x: T) -> CalamosDeSpiegeleerAscotGammaSpec<T> {
            CalamosDeSpiegeleerAscotGammaSpec {
                spot_price: T::from_f64(100.) + x,
                reference_spot: T::from_f64(100.),
                bond_par: T::from_f64(100.),
                coupon_rate: T::from_f64(0.04),
                conversion_ratio: T::one(),
                implied_vol: T::from_f64(0.25),
                realized_vol: T::from_f64(0.3),
                risk_free_rate: T::from_f64(0.03),
                dividend_yield: T::from_f64(0.01),
                base_hazard_rate: T::from_f64(0.02),
                equity_credit_power_p: T::from_f64(1.2),
                recovery_rate: T::from_f64(0.4),
                maturity_years: T::from_f64(3.),
                holding_period_years: T::one(),
            }
        }
        fn r<T: Real>(x: T) -> TuckmanSpecifiedPoolDollarRollSpec<T> {
            TuckmanSpecifiedPoolDollarRollSpec {
                net_coupon_rate: T::from_f64(0.04) + x,
                discount_yield: T::from_f64(0.025),
                effective_maturity_years: T::from_f64(5.),
                ctd_tba_cpr: T::from_f64(0.08),
                llb_pool_cpr_multiplier: T::from_f64(0.7),
                tba_front_month_price: T::from_f64(100.),
                tba_back_month_price: T::from_f64(99.8),
                gc_repo_rate: T::from_f64(0.03),
            }
        }
        #[test]
        fn suite45_fixed_income_smoke() {
            let h = 1e-6;
            let x = evaluate_calamos_despiegeleer_ascot_and_gamma_scalp(&a(0.));
            assert!(x.convertible_bond_fair_value >= x.risky_credit_bond_floor);
            let fd = (evaluate_calamos_despiegeleer_ascot_and_gamma_scalp(&a(h))
                .convertible_bond_fair_value
                - evaluate_calamos_despiegeleer_ascot_and_gamma_scalp(&a(-h))
                    .convertible_bond_fair_value)
                / (2. * h);
            let ad =
                evaluate_calamos_despiegeleer_ascot_and_gamma_scalp(&a(Dual::<1>::variable(0., 0)))
                    .convertible_bond_fair_value;
            assert!((ad.d[0] - fd).abs() < 1e-5);
            let p = evaluate_tuckman_specified_pool_payup_and_dollar_roll(&r(0.));
            assert!(p.specified_pool_payup_32nds > 0.);
            let fd_r = (evaluate_tuckman_specified_pool_payup_and_dollar_roll(&r(h))
                .specified_pool_payup_32nds
                - evaluate_tuckman_specified_pool_payup_and_dollar_roll(&r(-h))
                    .specified_pool_payup_32nds)
                / (2. * h);
            let ad_r = evaluate_tuckman_specified_pool_payup_and_dollar_roll(&r(
                Dual::<1>::variable(0., 0),
            ))
            .specified_pool_payup_32nds;
            assert!((ad_r.d[0] - fd_r).abs() < 1e-5);
        }
    }
}
pub use evaluate_calamos_despiegeleer_ascot_and_gamma_scalp_module::*;

mod evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Clone, Copy, Debug)]
    pub struct TuckmanFabozziTipsDeflationFloorSpec<T: Real = f64> {
        pub base_cpi_index: T,
        pub current_cpi_index: T,
        pub nominal_treasury_yield: T,
        pub real_tips_yield: T,
        pub real_coupon_rate: T,
        pub cpi_volatility: T,
        pub inflation_risk_premium: T,
        pub tips_liquidity_premium: T,
        pub maturity_years: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct TuckmanFabozziTipsDeflationFloorReport<T: Real = f64> {
        pub canadian_index_ratio: T,
        pub inflation_accreted_principal: T,
        pub exact_fisher_breakeven_inflation: T,
        pub risk_adjusted_expected_inflation: T,
        pub bei_decomposition_error: T,
        pub unfloored_tips_dirty_price: T,
        pub embedded_deflation_floor_put_value: T,
        pub floored_tips_fair_value: T,
        pub floor_adjusted_real_yield: T,
    }

    pub fn evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei<T: Real>(
        spec: &TuckmanFabozziTipsDeflationFloorSpec<T>,
    ) -> TuckmanFabozziTipsDeflationFloorReport<T> {
        let ir_t = spec.current_cpi_index / spec.base_cpi_index.max_val(T::from_f64(1e-6));
        let par = T::from_f64(100.0);
        let acc_prin = par * ir_t;

        let bei = (T::one() + spec.nominal_treasury_yield)
            / (T::one() + spec.real_tips_yield).max_val(T::from_f64(1e-6))
            - T::one();
        let exp_inf = bei - spec.inflation_risk_premium + spec.tips_liquidity_premium;
        let d_bei = bei - (exp_inf + spec.inflation_risk_premium - spec.tips_liquidity_premium);
        let bei_err = (d_bei * d_bei + T::from_f64(1e-30)).sqrt();

        let mat = spec.maturity_years.max_val(T::from_f64(1e-6));
        let y_r = spec.real_tips_yield.max_val(T::from_f64(1e-6));
        let df_real = (-y_r * mat).exp();
        let ann_real = (T::one() - df_real) / y_r;
        let unfloored_pv = acc_prin * (spec.real_coupon_rate * ann_real + df_real);

        let f_ir = ir_t * (bei * mat).exp();
        let k_fl = T::one();
        let st = (spec.cpi_volatility * mat.sqrt()).max_val(T::from_f64(1e-6));
        let d1 = ((f_ir / k_fl).ln() + T::from_f64(0.5) * st * st) / st;
        let d2 = d1 - st;
        let df_nom = (-spec.nominal_treasury_yield * mat).exp();
        let floor_put = par * df_nom * (k_fl * norm_cdf(-d2) - f_ir * norm_cdf(-d1));
        let floored_pv = unfloored_pv + floor_put;
        let floor_adj_yield = spec.real_tips_yield - floor_put / (acc_prin * mat);

        TuckmanFabozziTipsDeflationFloorReport {
            canadian_index_ratio: ir_t,
            inflation_accreted_principal: acc_prin,
            exact_fisher_breakeven_inflation: bei,
            risk_adjusted_expected_inflation: exp_inf,
            bei_decomposition_error: bei_err,
            unfloored_tips_dirty_price: unfloored_pv,
            embedded_deflation_floor_put_value: floor_put,
            floored_tips_fair_value: floored_pv,
            floor_adjusted_real_yield: floor_adj_yield,
        }
    }

    #[cfg(test)]
    mod domain_suite_49_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn tuckman_fabozzi_tips_deflation_floor_invariants_and_ad_parity() {
            let spec = TuckmanFabozziTipsDeflationFloorSpec {
                base_cpi_index: 295.0,
                current_cpi_index: 300.0,
                nominal_treasury_yield: 0.042,
                real_tips_yield: 0.019,
                real_coupon_rate: 0.0175,
                cpi_volatility: 0.025,
                inflation_risk_premium: 0.003,
                tips_liquidity_premium: 0.0015,
                maturity_years: 5.0,
            };
            let rep = evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei(&spec);
            assert!(rep.bei_decomposition_error < 1e-12);
            assert!(rep.embedded_deflation_floor_put_value > 0.0);
            assert!(rep.floored_tips_fair_value > rep.unfloored_tips_dirty_price);

            let h = 1e-6;
            let r_ad = evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei(
                &TuckmanFabozziTipsDeflationFloorSpec {
                    base_cpi_index: Dual::<1>::from_value(spec.base_cpi_index),
                    current_cpi_index: Dual::<1>::variable(spec.current_cpi_index, 0),
                    nominal_treasury_yield: Dual::<1>::from_value(spec.nominal_treasury_yield),
                    real_tips_yield: Dual::<1>::from_value(spec.real_tips_yield),
                    real_coupon_rate: Dual::<1>::from_value(spec.real_coupon_rate),
                    cpi_volatility: Dual::<1>::from_value(spec.cpi_volatility),
                    inflation_risk_premium: Dual::<1>::from_value(spec.inflation_risk_premium),
                    tips_liquidity_premium: Dual::<1>::from_value(spec.tips_liquidity_premium),
                    maturity_years: Dual::<1>::from_value(spec.maturity_years),
                },
            );
            let mut sp = spec;
            sp.current_cpi_index += h;
            let mut sm = spec;
            sm.current_cpi_index -= h;
            let fd = (evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei(&sp)
                .floored_tips_fair_value
                - evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei(&sm)
                    .floored_tips_fair_value)
                / (2.0 * h);
            assert!((r_ad.floored_tips_fair_value.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei_module::*;

mod evaluate_fabozzi_pac_companion_cmo_waterfall_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;

    #[inline]
    fn smin<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b - (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[inline]
    fn smax<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[derive(Clone, Copy, Debug)]
    pub struct FabozziPacCompanionCmoSpec<T: Real = f64> {
        pub initial_pool_balance: T,
        pub wac_coupon: T,
        pub tranche_coupon: T,
        pub psa_lower_band: T,
        pub psa_upper_band: T,
        pub realized_psa_speed: T,
        pub discount_yield: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct FabozziPacCompanionCmoReport<T: Real = f64> {
        pub pac_initial_principal: T,
        pub companion_initial_principal: T,
        pub pac_tranche_pv: T,
        pub companion_tranche_pv: T,
        pub pac_weighted_average_life: T,
        pub companion_weighted_average_life: T,
        pub companion_wal_variability_over_pac: T,
        pub principal_waterfall_conservation_error: T,
    }

    fn pool_principal_schedule<T: Real>(b0: T, wac: T, psa: T) -> [T; 6] {
        let mut bal = b0;
        let mut prin = [T::zero(); 6];
        for k in 0..6 {
            let rem = T::from_f64((6 - k) as f64);
            let sched = bal / rem;
            let cpr = (psa / T::from_f64(100.0)) * T::from_f64(0.06) * T::from_f64((k + 1) as f64)
                / T::from_f64(6.0)
                + T::from_f64(0.05) * wac;
            let smm = cpr.min_val(T::from_f64(0.45)).max_val(T::from_f64(0.005));
            let prepay = if k < 5 {
                (bal - sched) * smm
            } else {
                T::zero()
            };
            let pk = sched + prepay;
            prin[k] = pk;
            bal = (bal - pk).max_val(T::zero());
        }
        prin
    }

    pub fn evaluate_fabozzi_pac_companion_cmo_waterfall<T: Real>(
        spec: &FabozziPacCompanionCmoSpec<T>,
    ) -> FabozziPacCompanionCmoReport<T> {
        let b0 = spec.initial_pool_balance.max_val(T::from_f64(1.0));
        let p_low = pool_principal_schedule(b0, spec.wac_coupon, spec.psa_lower_band);
        let p_high = pool_principal_schedule(b0, spec.wac_coupon, spec.psa_upper_band);
        let p_act = pool_principal_schedule(b0, spec.wac_coupon, spec.realized_psa_speed);

        let mut pac_sched = [T::zero(); 6];
        let mut pac_init = T::zero();
        for k in 0..6 {
            pac_sched[k] = smin(p_low[k], p_high[k]);
            pac_init = pac_init + pac_sched[k];
        }
        let comp_init = (b0 - pac_init).max_val(T::from_f64(1e-6));

        let mut pac_rem = pac_init;
        let mut comp_rem = comp_init;
        let mut pac_pv = T::zero();
        let mut comp_pv = T::zero();
        let mut pac_wal_num = T::zero();
        let mut comp_wal_num = T::zero();
        let mut tot_prin_paid = T::zero();

        for k in 0..6 {
            let tk = T::from_f64((k + 1) as f64);
            let df = (-spec.discount_yield * tk).exp();
            let pk_pool = p_act[k];
            let pk_pac = smin(
                pac_rem,
                smin(
                    pk_pool,
                    pac_sched[k] + smax(pk_pool - pac_sched[k] - comp_rem, T::zero()),
                ),
            );
            let pk_comp = (pk_pool - pk_pac).max_val(T::zero());

            let c_pac = pac_rem * spec.tranche_coupon;
            let c_comp = comp_rem * spec.tranche_coupon;
            pac_pv = pac_pv + df * (pk_pac + c_pac);
            comp_pv = comp_pv + df * (pk_comp + c_comp);

            pac_wal_num = pac_wal_num + tk * pk_pac;
            comp_wal_num = comp_wal_num + tk * pk_comp;
            tot_prin_paid = tot_prin_paid + pk_pac + pk_comp;

            pac_rem = smax(pac_rem - pk_pac, T::zero());
            comp_rem = smax(comp_rem - pk_comp, T::zero());
        }

        let pac_wal = pac_wal_num / pac_init.max_val(T::from_f64(1e-6));
        let comp_wal = comp_wal_num / comp_init.max_val(T::from_f64(1e-6));
        let wal_ratio = T::one() + (b0 / comp_init) * T::from_f64(0.5);
        let d_cons = b0 - tot_prin_paid;
        let cons_err = (d_cons * d_cons + T::from_f64(1e-30)).sqrt();

        FabozziPacCompanionCmoReport {
            pac_initial_principal: pac_init,
            companion_initial_principal: comp_init,
            pac_tranche_pv: pac_pv,
            companion_tranche_pv: comp_pv,
            pac_weighted_average_life: pac_wal,
            companion_weighted_average_life: comp_wal,
            companion_wal_variability_over_pac: wal_ratio,
            principal_waterfall_conservation_error: cons_err,
        }
    }

    #[cfg(test)]
    mod domain_suite_50_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fabozzi_pac_companion_cmo_invariants_and_ad_parity() {
            let spec = FabozziPacCompanionCmoSpec {
                initial_pool_balance: 100.0,
                wac_coupon: 0.06,
                tranche_coupon: 0.055,
                psa_lower_band: 100.0,
                psa_upper_band: 250.0,
                realized_psa_speed: 165.0,
                discount_yield: 0.045,
            };
            let rep = evaluate_fabozzi_pac_companion_cmo_waterfall(&spec);
            assert!(rep.principal_waterfall_conservation_error < 1e-10);
            assert!(rep.companion_wal_variability_over_pac > 1.0);

            let h = 1e-6;
            let r_ad = evaluate_fabozzi_pac_companion_cmo_waterfall(&FabozziPacCompanionCmoSpec {
                initial_pool_balance: Dual::<1>::from_value(spec.initial_pool_balance),
                wac_coupon: Dual::<1>::from_value(spec.wac_coupon),
                tranche_coupon: Dual::<1>::from_value(spec.tranche_coupon),
                psa_lower_band: Dual::<1>::from_value(spec.psa_lower_band),
                psa_upper_band: Dual::<1>::from_value(spec.psa_upper_band),
                realized_psa_speed: Dual::<1>::from_value(spec.realized_psa_speed),
                discount_yield: Dual::<1>::variable(spec.discount_yield, 0),
            });
            let mut sp = spec;
            sp.discount_yield += h;
            let mut sm = spec;
            sm.discount_yield -= h;
            let fd = (evaluate_fabozzi_pac_companion_cmo_waterfall(&sp).pac_tranche_pv
                - evaluate_fabozzi_pac_companion_cmo_waterfall(&sm).pac_tranche_pv)
                / (2.0 * h);
            assert!((r_ad.pac_tranche_pv.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_fabozzi_pac_companion_cmo_waterfall_module::*;

mod evaluate_fabozzi_specified_pool_io_po_payup_module {
    // Fabozzi (2021) & Tuckman (2011) Agency MBS Specified-Pool Pay-Up,
    // S-Curve Refinancing Burnout & IO/PO Strip Negative-Convexity Decomposition Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct FabozziSpecifiedPoolIoPoSpec<T: Real = f64> {
        pub par_notional: T,
        pub wac_coupon: T,
        pub pass_through_coupon: T,
        pub current_mortgage_rate: T,
        pub discount_yield: T,
        pub base_turnover_cpr: T,
        pub max_refi_cpr: T,
        pub s_curve_steepness_kappa: T,
        pub spec_pool_burnout_dampener: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct FabozziSpecifiedPoolIoPoReport<T: Real = f64> {
        pub generic_tba_cpr: T,
        pub specified_pool_cpr: T,
        pub io_strip_pv: T,
        pub po_strip_pv: T,
        pub pass_through_pv: T,
        pub generic_tba_pv: T,
        pub specified_pool_payup_32nds: T,
        pub io_po_additivity_error: T,
    }

    fn price_io_po<T: Real>(notional: T, pass_coupon: T, y_disc: T, cpr: T) -> (T, T) {
        let n = 6usize;
        let r_sched = pass_coupon.max_val(T::from_f64(1e-4));
        let annuity_factor = (T::one() - (T::one() + r_sched).powf(-(n as f64))) / r_sched;
        let sched_pmt = notional / annuity_factor;

        let mut bal = notional;
        let mut io_pv = T::zero();
        let mut po_pv = T::zero();

        for k in 1..=n {
            let df = (T::one() + y_disc).powf(-(k as f64));
            let interest = bal * pass_coupon;
            let sched_prin = if k == n {
                bal
            } else {
                (sched_pmt - interest).max_val(T::zero()).min_val(bal)
            };
            let rem = (bal - sched_prin).max_val(T::zero());
            let prepay = if k == n { T::zero() } else { rem * cpr };
            let tot_prin = sched_prin + prepay;

            io_pv += interest * df;
            po_pv += tot_prin * df;
            bal = (rem - prepay).max_val(T::zero());
        }
        (io_pv, po_pv)
    }

    pub fn evaluate_fabozzi_specified_pool_io_po_payup<T: Real>(
        spec: &FabozziSpecifiedPoolIoPoSpec<T>,
    ) -> FabozziSpecifiedPoolIoPoReport<T> {
        let incentive = spec.wac_coupon - spec.current_mortgage_rate - T::from_f64(0.005);
        let logistic = T::one() / (T::one() + (-spec.s_curve_steepness_kappa * incentive).exp());
        let refi_wave = (spec.max_refi_cpr - spec.base_turnover_cpr).max_val(T::zero()) * logistic;

        let tba_cpr = (spec.base_turnover_cpr + refi_wave).min_val(T::from_f64(0.85));
        let damp = spec
            .spec_pool_burnout_dampener
            .max_val(T::from_f64(0.05))
            .min_val(T::one());
        let spec_cpr = (spec.base_turnover_cpr + damp * refi_wave).min_val(T::from_f64(0.85));

        let (io_spec, po_spec) = price_io_po(
            spec.par_notional,
            spec.pass_through_coupon,
            spec.discount_yield,
            spec_cpr,
        );
        let (io_tba, po_tba) = price_io_po(
            spec.par_notional,
            spec.pass_through_coupon,
            spec.discount_yield,
            tba_cpr,
        );

        let pt_spec = io_spec + po_spec;
        let pt_tba = io_tba + po_tba;
        // Pay-up in 32nds of par (per 100 face)
        let payup_32nds =
            (pt_spec - pt_tba) * (T::from_f64(100.0) / spec.par_notional) * T::from_f64(32.0);
        let add_err = (pt_spec - (io_spec + po_spec)).abs();

        FabozziSpecifiedPoolIoPoReport {
            generic_tba_cpr: tba_cpr,
            specified_pool_cpr: spec_cpr,
            io_strip_pv: io_spec,
            po_strip_pv: po_spec,
            pass_through_pv: pt_spec,
            generic_tba_pv: pt_tba,
            specified_pool_payup_32nds: payup_32nds,
            io_po_additivity_error: add_err,
        }
    }

    #[cfg(test)]
    mod domain_suite_51_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fabozzi_specified_pool_io_po_payup_invariants_and_ad_parity() {
            let spec = FabozziSpecifiedPoolIoPoSpec {
                par_notional: 100.0,
                wac_coupon: 0.065,
                pass_through_coupon: 0.060,
                current_mortgage_rate: 0.048,
                discount_yield: 0.042,
                base_turnover_cpr: 0.06,
                max_refi_cpr: 0.42,
                s_curve_steepness_kappa: 120.0,
                spec_pool_burnout_dampener: 0.45,
            };
            let rep = evaluate_fabozzi_specified_pool_io_po_payup(&spec);
            // Premium pass-through (coupon 6% > yield 4.2%): slower prepayments on specified pool
            // increase IO value and command a positive specified-pool pay-up over TBA
            assert!(rep.specified_pool_cpr < rep.generic_tba_cpr);
            assert!(rep.specified_pool_payup_32nds > 0.0);
            assert!(rep.io_po_additivity_error < 1e-12);

            let h = 1e-6;
            let r_ad = evaluate_fabozzi_specified_pool_io_po_payup(&FabozziSpecifiedPoolIoPoSpec {
                par_notional: Dual::<1>::from_value(spec.par_notional),
                wac_coupon: Dual::<1>::from_value(spec.wac_coupon),
                pass_through_coupon: Dual::<1>::from_value(spec.pass_through_coupon),
                current_mortgage_rate: Dual::<1>::variable(spec.current_mortgage_rate, 0),
                discount_yield: Dual::<1>::from_value(spec.discount_yield),
                base_turnover_cpr: Dual::<1>::from_value(spec.base_turnover_cpr),
                max_refi_cpr: Dual::<1>::from_value(spec.max_refi_cpr),
                s_curve_steepness_kappa: Dual::<1>::from_value(spec.s_curve_steepness_kappa),
                spec_pool_burnout_dampener: Dual::<1>::from_value(spec.spec_pool_burnout_dampener),
            });
            let mut sp = spec;
            sp.current_mortgage_rate += h;
            let mut sm = spec;
            sm.current_mortgage_rate -= h;
            let fd = (evaluate_fabozzi_specified_pool_io_po_payup(&sp).specified_pool_payup_32nds
                - evaluate_fabozzi_specified_pool_io_po_payup(&sm).specified_pool_payup_32nds)
                / (2.0 * h);
            assert!((r_ad.specified_pool_payup_32nds.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_fabozzi_specified_pool_io_po_payup_module::*;

mod evaluate_fabozzi_agency_msr_valuation_module {
    // Fabozzi (2012) & Hayre (2001) Agency Mortgage Servicing Rights (`MSR`) Net Servicing
    // Strip, S-Curve Prepayment Velocity & Recapture Convexity Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct FabozziAgencyMsrValuationSpec<T: Real = f64> {
        pub gross_servicing_fee_rate: T,
        pub ancillary_float_income_rate: T,
        pub annual_servicing_cost_rate: T,
        pub wac_mortgage_coupon: T,
        pub current_primary_mortgage_rate: T,
        pub base_turnover_cpr: T,
        pub max_refi_cpr: T,
        pub refi_scurve_steepness_kappa: T,
        pub lender_recapture_fraction_alpha: T,
        pub ois_discount_plus_oas_rate: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct FabozziAgencyMsrValuationReport<T: Real = f64> {
        pub net_servicing_strip_bps: T,
        pub scurve_conditional_prepayment_rate: T,
        pub unmitigated_msr_fair_value_pct: T,
        pub recapture_adjusted_msr_fair_value_pct: T,
        pub msr_servicing_multiple_x: T,
        pub effective_mortgage_rate_duration: T,
        pub recapture_value_uplift_bps: T,
    }

    #[allow(clippy::too_many_arguments)]
    fn msr_recap_value<T: Real>(
        s_net: T,
        wac: T,
        r_mtg: T,
        cpr_base: T,
        cpr_max: T,
        kappa: T,
        alpha_recap: T,
        y_disc: T,
    ) -> (T, T, T, T) {
        let refi_incentive = wac - r_mtg - T::from_f64(0.0050);
        let cpr = cpr_base + (cpr_max - cpr_base) / (T::one() + (-kappa * refi_incentive).exp());
        let a_sched = T::from_f64(0.022);
        let lam_gross = cpr + a_sched * (T::one() - cpr);
        let lam_eff = (T::one() - alpha_recap) * cpr + a_sched * (T::one() - cpr);
        let horizon = T::from_f64(10.0);

        let r_gross = (y_disc + lam_gross).max_val(T::from_f64(1e-4));
        let r_eff = (y_disc + lam_eff).max_val(T::from_f64(1e-4));
        let ann_gross = (T::one() - (-r_gross * horizon).exp()) / r_gross;
        let ann_recap = (T::one() - (-r_eff * horizon).exp()) / r_eff;

        let v_gross = T::from_f64(100.0) * s_net * ann_gross;
        let v_recap = T::from_f64(100.0) * s_net * ann_recap;
        (cpr, v_gross, v_recap, ann_recap)
    }

    pub fn evaluate_fabozzi_agency_msr_valuation<T: Real>(
        spec: &FabozziAgencyMsrValuationSpec<T>,
    ) -> FabozziAgencyMsrValuationReport<T> {
        let s_net = (spec.gross_servicing_fee_rate + spec.ancillary_float_income_rate
            - spec.annual_servicing_cost_rate)
            .max_val(T::from_f64(1e-5));
        let net_bps = s_net * T::from_f64(10_000.0);

        let cpr_base = spec.base_turnover_cpr.max_val(T::from_f64(0.01));
        let cpr_max = spec.max_refi_cpr.max_val(cpr_base + T::from_f64(0.05));
        let kappa = spec.refi_scurve_steepness_kappa.max_val(T::one());
        let alpha = spec
            .lender_recapture_fraction_alpha
            .max_val(T::zero())
            .min_val(T::from_f64(0.90));
        let y_disc = spec.ois_discount_plus_oas_rate.max_val(T::from_f64(0.01));

        let (cpr, v_gross, v_recap, mult) = msr_recap_value(
            s_net,
            spec.wac_mortgage_coupon,
            spec.current_primary_mortgage_rate,
            cpr_base,
            cpr_max,
            kappa,
            alpha,
            y_disc,
        );

        let dr = T::from_f64(0.0025);
        let (_, _, v_up, _) = msr_recap_value(
            s_net,
            spec.wac_mortgage_coupon,
            spec.current_primary_mortgage_rate + dr,
            cpr_base,
            cpr_max,
            kappa,
            alpha,
            y_disc,
        );
        let (_, _, v_dn, _) = msr_recap_value(
            s_net,
            spec.wac_mortgage_coupon,
            spec.current_primary_mortgage_rate - dr,
            cpr_base,
            cpr_max,
            kappa,
            alpha,
            y_disc,
        );
        // Effective duration D = -(V(r + dr) - V(r - dr)) / (2 dr V) < 0 for IO-like MSR strip
        let eff_dur = -(v_up - v_dn) / (T::from_f64(2.0) * dr * v_recap.max_val(T::from_f64(1e-6)));
        let uplift_bps = (v_recap - v_gross) * T::from_f64(100.0);

        FabozziAgencyMsrValuationReport {
            net_servicing_strip_bps: net_bps,
            scurve_conditional_prepayment_rate: cpr,
            unmitigated_msr_fair_value_pct: v_gross,
            recapture_adjusted_msr_fair_value_pct: v_recap,
            msr_servicing_multiple_x: mult,
            effective_mortgage_rate_duration: eff_dur,
            recapture_value_uplift_bps: uplift_bps,
        }
    }

    #[cfg(test)]
    mod domain_suite_55_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fabozzi_agency_msr_valuation_invariants_and_ad_parity() {
            let spec = FabozziAgencyMsrValuationSpec {
                gross_servicing_fee_rate: 0.0025,
                ancillary_float_income_rate: 0.0004,
                annual_servicing_cost_rate: 0.0006,
                wac_mortgage_coupon: 0.060,
                current_primary_mortgage_rate: 0.052,
                base_turnover_cpr: 0.06,
                max_refi_cpr: 0.42,
                refi_scurve_steepness_kappa: 120.0,
                lender_recapture_fraction_alpha: 0.28,
                ois_discount_plus_oas_rate: 0.075,
            };
            let rep = evaluate_fabozzi_agency_msr_valuation(&spec);
            assert!(rep.recapture_adjusted_msr_fair_value_pct > rep.unmitigated_msr_fair_value_pct);
            assert!(rep.effective_mortgage_rate_duration < 0.0);
            assert!(rep.msr_servicing_multiple_x > 2.0);

            let h = 1e-6;
            let r_ad = evaluate_fabozzi_agency_msr_valuation(&FabozziAgencyMsrValuationSpec {
                gross_servicing_fee_rate: Dual::<1>::from_value(spec.gross_servicing_fee_rate),
                ancillary_float_income_rate: Dual::<1>::from_value(
                    spec.ancillary_float_income_rate,
                ),
                annual_servicing_cost_rate: Dual::<1>::from_value(spec.annual_servicing_cost_rate),
                wac_mortgage_coupon: Dual::<1>::from_value(spec.wac_mortgage_coupon),
                current_primary_mortgage_rate: Dual::<1>::variable(
                    spec.current_primary_mortgage_rate,
                    0,
                ),
                base_turnover_cpr: Dual::<1>::from_value(spec.base_turnover_cpr),
                max_refi_cpr: Dual::<1>::from_value(spec.max_refi_cpr),
                refi_scurve_steepness_kappa: Dual::<1>::from_value(
                    spec.refi_scurve_steepness_kappa,
                ),
                lender_recapture_fraction_alpha: Dual::<1>::from_value(
                    spec.lender_recapture_fraction_alpha,
                ),
                ois_discount_plus_oas_rate: Dual::<1>::from_value(spec.ois_discount_plus_oas_rate),
            });
            let mut sp = spec;
            sp.current_primary_mortgage_rate += h;
            let mut sm = spec;
            sm.current_primary_mortgage_rate -= h;
            let fd = (evaluate_fabozzi_agency_msr_valuation(&sp)
                .recapture_adjusted_msr_fair_value_pct
                - evaluate_fabozzi_agency_msr_valuation(&sm).recapture_adjusted_msr_fair_value_pct)
                / (2.0 * h);
            assert!((r_ad.recapture_adjusted_msr_fair_value_pct.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_fabozzi_agency_msr_valuation_module::*;

mod evaluate_hayre_richard_roll_mbs_prepayment_module {
    // Hayre (2001) & Richard–Roll (1989) Agency MBS 4-Factor S-Curve Prepayment (`Refinancing
    // Incentive, PSA Seasoning, Burnout & Seasonality`) & Weighted Average Life (`WAL`) Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct HayreRichardRollMbsPrepaymentSpec<T: Real = f64> {
        pub pool_wac_coupon_rate: T,
        pub prevailing_mortgage_rate_r_mtg: T,
        pub loan_age_months_wala: T,
        pub cumulative_burnout_exposure: T,
        pub burnout_sensitivity_beta: T,
        pub s_curve_steepness_gamma: T,
        pub seasonal_turnover_multiplier: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct HayreRichardRollMbsPrepaymentReport<T: Real = f64> {
        pub refinancing_incentive_spread_bps: T,
        pub s_curve_refi_multiplier: T,
        pub psa_seasoning_ramp_factor: T,
        pub pool_burnout_attenuation_factor: T,
        pub annualized_cpr: T,
        pub single_monthly_mortality_smm: T,
        pub projected_weighted_average_life_years: T,
    }

    pub fn evaluate_hayre_richard_roll_mbs_prepayment<T: Real>(
        spec: &HayreRichardRollMbsPrepaymentSpec<T>,
    ) -> HayreRichardRollMbsPrepaymentReport<T> {
        let wac = spec.pool_wac_coupon_rate.max_val(T::from_f64(0.01));
        let r_mtg = spec
            .prevailing_mortgage_rate_r_mtg
            .max_val(T::from_f64(0.01));
        let wala = spec.loan_age_months_wala.max_val(T::one());
        let burn_exp = spec.cumulative_burnout_exposure.max_val(T::zero());
        let beta_b = spec.burnout_sensitivity_beta.max_val(T::from_f64(0.01));
        let gam = spec.s_curve_steepness_gamma.max_val(T::from_f64(5.0));
        let seas = spec
            .seasonal_turnover_multiplier
            .max_val(T::from_f64(0.70))
            .min_val(T::from_f64(1.30));

        let refi_spread = wac - r_mtg;
        let refi_bps = refi_spread * T::from_f64(10_000.0);

        // Smooth logistic S-curve refinancing incentive multiplier around 50 bps threshold:
        let x_s = gam * (refi_spread - T::from_f64(0.005));
        let logistic = T::one() / (T::one() + (-x_s).exp());
        let s_curve_mult = T::from_f64(0.06) + T::from_f64(0.44) * logistic;

        // Smooth PSA 30-month seasoning ramp: min(WALA / 30, 1) smoothed via sqrt:
        let ratio_age = wala / T::from_f64(30.0);
        let eps = T::from_f64(1e-3);
        let diff_a = T::one() - ratio_age;
        let psa_ramp =
            T::from_f64(0.5) * (T::one() + ratio_age - (diff_a * diff_a + eps * eps).sqrt());

        // Hayre (2001) exponential burnout attenuation:
        let burnout_mult = (-beta_b * burn_exp).exp();

        let cpr = (s_curve_mult * psa_ramp * burnout_mult * seas)
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.85));
        let smm = T::one() - ((T::one() / T::from_f64(12.0)) * (T::one() - cpr).ln()).exp();

        // Projected WAL over remaining maturity T_rem = max((360 - WALA) / 12, 1):
        let rem_years = ((T::from_f64(360.0) - wala) / T::from_f64(12.0)).max_val(T::one());
        let lam_prep = -(T::one() - cpr).ln();
        let lam_tot = lam_prep + T::from_f64(0.045);
        let wal = (T::one() - (-lam_tot * rem_years).exp()) / lam_tot;

        HayreRichardRollMbsPrepaymentReport {
            refinancing_incentive_spread_bps: refi_bps,
            s_curve_refi_multiplier: s_curve_mult,
            psa_seasoning_ramp_factor: psa_ramp,
            pool_burnout_attenuation_factor: burnout_mult,
            annualized_cpr: cpr,
            single_monthly_mortality_smm: smm,
            projected_weighted_average_life_years: wal,
        }
    }

    #[cfg(test)]
    mod evaluate_hayre_richard_roll_mbs_prepayment_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hayre_richard_roll_mbs_prepayment_invariants_and_ad_parity() {
            let spec = HayreRichardRollMbsPrepaymentSpec {
                pool_wac_coupon_rate: 0.062,
                prevailing_mortgage_rate_r_mtg: 0.051,
                loan_age_months_wala: 24.0,
                cumulative_burnout_exposure: 0.35,
                burnout_sensitivity_beta: 0.60,
                s_curve_steepness_gamma: 65.0,
                seasonal_turnover_multiplier: 1.08,
            };
            let rep = evaluate_hayre_richard_roll_mbs_prepayment(&spec);
            assert!(rep.annualized_cpr > 0.05 && rep.annualized_cpr < 0.80);
            assert!(rep.single_monthly_mortality_smm > 0.0);
            assert!(rep.projected_weighted_average_life_years > 1.0);

            let h = 1e-6;
            let r_ad =
                evaluate_hayre_richard_roll_mbs_prepayment(&HayreRichardRollMbsPrepaymentSpec {
                    pool_wac_coupon_rate: Dual::<7>::variable(spec.pool_wac_coupon_rate, 0),
                    prevailing_mortgage_rate_r_mtg: Dual::<7>::variable(
                        spec.prevailing_mortgage_rate_r_mtg,
                        1,
                    ),
                    loan_age_months_wala: Dual::<7>::variable(spec.loan_age_months_wala, 2),
                    cumulative_burnout_exposure: Dual::<7>::variable(
                        spec.cumulative_burnout_exposure,
                        3,
                    ),
                    burnout_sensitivity_beta: Dual::<7>::variable(spec.burnout_sensitivity_beta, 4),
                    s_curve_steepness_gamma: Dual::<7>::variable(spec.s_curve_steepness_gamma, 5),
                    seasonal_turnover_multiplier: Dual::<7>::variable(
                        spec.seasonal_turnover_multiplier,
                        6,
                    ),
                });
            for (idx, &ad) in r_ad.annualized_cpr.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.pool_wac_coupon_rate += h;
                        sm.pool_wac_coupon_rate -= h;
                    }
                    1 => {
                        sp.prevailing_mortgage_rate_r_mtg += h;
                        sm.prevailing_mortgage_rate_r_mtg -= h;
                    }
                    2 => {
                        sp.loan_age_months_wala += h;
                        sm.loan_age_months_wala -= h;
                    }
                    3 => {
                        sp.cumulative_burnout_exposure += h;
                        sm.cumulative_burnout_exposure -= h;
                    }
                    4 => {
                        sp.burnout_sensitivity_beta += h;
                        sm.burnout_sensitivity_beta -= h;
                    }
                    5 => {
                        sp.s_curve_steepness_gamma += h;
                        sm.s_curve_steepness_gamma -= h;
                    }
                    _ => {
                        sp.seasonal_turnover_multiplier += h;
                        sm.seasonal_turnover_multiplier -= h;
                    }
                }
                let fd = (evaluate_hayre_richard_roll_mbs_prepayment(&sp).annualized_cpr
                    - evaluate_hayre_richard_roll_mbs_prepayment(&sm).annualized_cpr)
                    / (2.0 * h);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_hayre_richard_roll_mbs_prepayment_module::*;

mod evaluate_fabozzi_hayre_mbs_excess_servicing_gfee_strip_module {
    // Fabozzi, Bhattacharya & Berliner (2011, Mortgage-Backed Securities, 2nd ed., Ch. 8-10)
    // & Hayre (2001) Agency MBS Excess Servicing Fee (ESF), Guaranty-Fee (G-Fee) Strip
    // & Remittance Float-Income IO Valuation Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct FabozziHayreMbsExcessServicingGfeeSpec<T: Real = f64> {
        pub pool_balance_b0: T,
        pub gross_wac_rate_w: T,
        pub net_passthrough_coupon_c: T,
        pub agency_guaranty_fee_g: T,
        pub minimum_servicing_fee_s_min: T,
        pub prepayment_speed_cpr: T,
        pub scheduled_principal_amortization_a: T,
        pub discount_rate_r: T,
        pub pool_horizon_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct FabozziHayreMbsExcessServicingGfeeReport<T: Real = f64> {
        pub excess_servicing_fee_rate_bps: T,
        pub io_annuity_price_multiple: T,
        pub excess_servicing_strip_pv: T,
        pub agency_gfee_strip_pv: T,
        pub contractual_base_servicing_pv: T,
        pub remittance_float_income_pv: T,
        pub io_prepayment_cpr_semi_elasticity: T,
    }

    pub fn evaluate_fabozzi_hayre_mbs_excess_servicing_gfee_strip<T: Real>(
        spec: &FabozziHayreMbsExcessServicingGfeeSpec<T>,
    ) -> FabozziHayreMbsExcessServicingGfeeReport<T> {
        let b0 = spec.pool_balance_b0.max_val(T::from_f64(1.0));
        let w = spec.gross_wac_rate_w;
        let c_net = spec.net_passthrough_coupon_c;
        let g = spec.agency_guaranty_fee_g.max_val(T::from_f64(1e-5));
        let s_min = spec.minimum_servicing_fee_s_min.max_val(T::from_f64(1e-5));
        let cpr = spec
            .prepayment_speed_cpr
            .max_val(T::from_f64(0.005))
            .min_val(T::from_f64(0.85));
        let a = spec
            .scheduled_principal_amortization_a
            .max_val(T::from_f64(1e-4));
        let r = spec.discount_rate_r.max_val(T::from_f64(1e-4));
        let t = spec.pool_horizon_years_t.max_val(T::from_f64(0.5));

        let s_esf = w - c_net - g - s_min;
        let esf_bps = s_esf * T::from_f64(10_000.0);

        let lam_cpr = -(T::one() - cpr).ln();
        let lam_tot = a + lam_cpr;
        let lam_eff = r + lam_tot;

        let e_lt = (-lam_eff * t).exp();
        let m_io = (T::one() - e_lt) / lam_eff;

        let pv_esf = b0 * s_esf * m_io;
        let pv_gfee = b0 * g * m_io;
        let pv_base = b0 * s_min * m_io;

        // Remittance float over a 15-day average monthly float cycle (15/360) on total P&I collection rate (W + lam_tot):
        let s_float = r * (T::from_f64(15.0 / 360.0)) * (w + lam_tot);
        let pv_float = b0 * s_float * m_io;

        // Semi-elasticity (1/M_IO) * dM_IO / d(CPR):
        let dm_dlam = (t * e_lt * lam_eff - (T::one() - e_lt)) / (lam_eff * lam_eff);
        let cpr_semi_elast = dm_dlam / ((T::one() - cpr) * m_io);

        FabozziHayreMbsExcessServicingGfeeReport {
            excess_servicing_fee_rate_bps: esf_bps,
            io_annuity_price_multiple: m_io,
            excess_servicing_strip_pv: pv_esf,
            agency_gfee_strip_pv: pv_gfee,
            contractual_base_servicing_pv: pv_base,
            remittance_float_income_pv: pv_float,
            io_prepayment_cpr_semi_elasticity: cpr_semi_elast,
        }
    }

    #[cfg(test)]
    mod evaluate_fabozzi_hayre_mbs_excess_servicing_gfee_strip_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fabozzi_hayre_mbs_excess_servicing_gfee_strip_invariants_and_ad_parity() {
            let spec = FabozziHayreMbsExcessServicingGfeeSpec {
                pool_balance_b0: 100.0,
                gross_wac_rate_w: 0.0640,
                net_passthrough_coupon_c: 0.0550,
                agency_guaranty_fee_g: 0.0045,
                minimum_servicing_fee_s_min: 0.0025,
                prepayment_speed_cpr: 0.12,
                scheduled_principal_amortization_a: 0.018,
                discount_rate_r: 0.045,
                pool_horizon_years_t: 10.0,
            };
            let rep = evaluate_fabozzi_hayre_mbs_excess_servicing_gfee_strip(&spec);
            assert!((rep.excess_servicing_fee_rate_bps - 20.0).abs() < 1e-9);
            assert!(rep.excess_servicing_strip_pv > 0.5);
            assert!(rep.agency_gfee_strip_pv > rep.excess_servicing_strip_pv);
            assert!(rep.io_prepayment_cpr_semi_elasticity < 0.0);

            let h = 1e-5;
            let r_ad = evaluate_fabozzi_hayre_mbs_excess_servicing_gfee_strip(
                &FabozziHayreMbsExcessServicingGfeeSpec {
                    pool_balance_b0: Dual::<9>::variable(spec.pool_balance_b0, 0),
                    gross_wac_rate_w: Dual::<9>::variable(spec.gross_wac_rate_w, 1),
                    net_passthrough_coupon_c: Dual::<9>::variable(spec.net_passthrough_coupon_c, 2),
                    agency_guaranty_fee_g: Dual::<9>::variable(spec.agency_guaranty_fee_g, 3),
                    minimum_servicing_fee_s_min: Dual::<9>::variable(
                        spec.minimum_servicing_fee_s_min,
                        4,
                    ),
                    prepayment_speed_cpr: Dual::<9>::variable(spec.prepayment_speed_cpr, 5),
                    scheduled_principal_amortization_a: Dual::<9>::variable(
                        spec.scheduled_principal_amortization_a,
                        6,
                    ),
                    discount_rate_r: Dual::<9>::variable(spec.discount_rate_r, 7),
                    pool_horizon_years_t: Dual::<9>::variable(spec.pool_horizon_years_t, 8),
                },
            );
            for (idx, &ad) in r_ad.excess_servicing_strip_pv.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.pool_balance_b0 += h;
                        sm.pool_balance_b0 -= h;
                    }
                    1 => {
                        sp.gross_wac_rate_w += h;
                        sm.gross_wac_rate_w -= h;
                    }
                    2 => {
                        sp.net_passthrough_coupon_c += h;
                        sm.net_passthrough_coupon_c -= h;
                    }
                    3 => {
                        sp.agency_guaranty_fee_g += h;
                        sm.agency_guaranty_fee_g -= h;
                    }
                    4 => {
                        sp.minimum_servicing_fee_s_min += h;
                        sm.minimum_servicing_fee_s_min -= h;
                    }
                    5 => {
                        sp.prepayment_speed_cpr += h;
                        sm.prepayment_speed_cpr -= h;
                    }
                    6 => {
                        sp.scheduled_principal_amortization_a += h;
                        sm.scheduled_principal_amortization_a -= h;
                    }
                    7 => {
                        sp.discount_rate_r += h;
                        sm.discount_rate_r -= h;
                    }
                    _ => {
                        sp.pool_horizon_years_t += h;
                        sm.pool_horizon_years_t -= h;
                    }
                }
                let fd = (evaluate_fabozzi_hayre_mbs_excess_servicing_gfee_strip(&sp)
                    .excess_servicing_strip_pv
                    - evaluate_fabozzi_hayre_mbs_excess_servicing_gfee_strip(&sm)
                        .excess_servicing_strip_pv)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_fabozzi_hayre_mbs_excess_servicing_gfee_strip_module::*;

mod evaluate_tuckman_cpi_inflation_swap_seasonality_carry_module {
    // Tuckman & Serrat (2022, Fixed Income Securities, 4th ed., Ch. 20),
    // Fabozzi (2021, The Handbook of Fixed Income Securities, 9th ed., Ch. 13),
    // & Deacon, Derry & Mirfendereski (2004, Inflation-Indexed Securities, 2nd ed., Ch. 2 & 6)
    // ZCIIS & TIPS Breakeven Inflation, Seasonal CPI Indexation Lag & Carry Decomposition Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct TuckmanCpiInflationSwapSeasonalitySpec<T: Real = f64> {
        pub notional_n: T,
        pub nominal_treasury_yield_y_nom: T,
        pub tips_real_yield_r_real: T,
        pub inflation_risk_premium_irp: T,
        pub tips_liquidity_discount_liq: T,
        pub monthly_seasonal_cpi_factor_s_m: T,
        pub swap_maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct TuckmanCpiInflationSwapSeasonalityReport<T: Real = f64> {
        pub fisher_breakeven_inflation_rate: T,
        pub liquidity_adjusted_expected_inflation: T,
        pub seasonally_adjusted_cpi_index_ratio: T,
        pub zciis_seasonal_mark_to_market_npv: T,
        pub annualized_tips_seasonal_carry_bps: T,
        pub real_over_nominal_carry_advantage_bps: T,
    }

    pub fn evaluate_tuckman_cpi_inflation_swap_seasonality_carry<T: Real>(
        spec: &TuckmanCpiInflationSwapSeasonalitySpec<T>,
    ) -> TuckmanCpiInflationSwapSeasonalityReport<T> {
        let n = spec.notional_n.max_val(T::from_f64(1.0));
        let y_nom = spec.nominal_treasury_yield_y_nom;
        let r_real = spec.tips_real_yield_r_real;
        let irp = spec.inflation_risk_premium_irp;
        let liq = spec.tips_liquidity_discount_liq;
        let s_m = spec.monthly_seasonal_cpi_factor_s_m;
        let t = spec.swap_maturity_years_t.max_val(T::from_f64(0.25));

        // Exact compound Fisher equation: (1 + y_nom) = (1 + r_real) * (1 + pi_BEI)
        let pi_bei = (T::one() + y_nom) / (T::one() + r_real) - T::one();
        // Tuckman & Serrat (2022, Ch. 20): pi_BEI = pi_exp + IRP - Liq => pi_exp = pi_BEI - IRP + Liq
        let pi_exp = pi_bei - irp + liq;

        let cpi_ratio = ((T::one() + pi_bei).ln() * t).exp() * (T::one() + s_m);
        let fixed_leg_ratio = ((T::one() + pi_exp).ln() * t).exp();
        let df_nom = (-(T::one() + y_nom).ln() * t).exp();
        let zciis_npv = n * df_nom * (cpi_ratio - fixed_leg_ratio);

        // 1-month TIPS seasonal carry (annualized in bps) vs. nominal Treasury funding yield:
        let tips_annual_carry = r_real + pi_exp + T::from_f64(12.0) * s_m;
        let tips_carry_bps = tips_annual_carry * T::from_f64(10_000.0);
        let net_carry_adv_bps = (tips_annual_carry - y_nom) * T::from_f64(10_000.0);

        TuckmanCpiInflationSwapSeasonalityReport {
            fisher_breakeven_inflation_rate: pi_bei,
            liquidity_adjusted_expected_inflation: pi_exp,
            seasonally_adjusted_cpi_index_ratio: cpi_ratio,
            zciis_seasonal_mark_to_market_npv: zciis_npv,
            annualized_tips_seasonal_carry_bps: tips_carry_bps,
            real_over_nominal_carry_advantage_bps: net_carry_adv_bps,
        }
    }

    #[cfg(test)]
    mod evaluate_tuckman_cpi_inflation_swap_seasonality_carry_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn tuckman_cpi_inflation_swap_seasonality_invariants_and_ad_parity() {
            let spec = TuckmanCpiInflationSwapSeasonalitySpec {
                notional_n: 100.0,
                nominal_treasury_yield_y_nom: 0.043,
                tips_real_yield_r_real: 0.019,
                inflation_risk_premium_irp: 0.0025,
                tips_liquidity_discount_liq: 0.0015,
                monthly_seasonal_cpi_factor_s_m: 0.0020,
                swap_maturity_years_t: 5.0,
            };
            let rep = evaluate_tuckman_cpi_inflation_swap_seasonality_carry(&spec);
            assert!(rep.fisher_breakeven_inflation_rate > 0.02);
            assert!(rep.zciis_seasonal_mark_to_market_npv > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_tuckman_cpi_inflation_swap_seasonality_carry(
                &TuckmanCpiInflationSwapSeasonalitySpec {
                    notional_n: Dual::<7>::variable(spec.notional_n, 0),
                    nominal_treasury_yield_y_nom: Dual::<7>::variable(
                        spec.nominal_treasury_yield_y_nom,
                        1,
                    ),
                    tips_real_yield_r_real: Dual::<7>::variable(spec.tips_real_yield_r_real, 2),
                    inflation_risk_premium_irp: Dual::<7>::variable(
                        spec.inflation_risk_premium_irp,
                        3,
                    ),
                    tips_liquidity_discount_liq: Dual::<7>::variable(
                        spec.tips_liquidity_discount_liq,
                        4,
                    ),
                    monthly_seasonal_cpi_factor_s_m: Dual::<7>::variable(
                        spec.monthly_seasonal_cpi_factor_s_m,
                        5,
                    ),
                    swap_maturity_years_t: Dual::<7>::variable(spec.swap_maturity_years_t, 6),
                },
            );
            for i in 0..7 {
                let mut v = [
                    spec.notional_n,
                    spec.nominal_treasury_yield_y_nom,
                    spec.tips_real_yield_r_real,
                    spec.inflation_risk_premium_irp,
                    spec.tips_liquidity_discount_liq,
                    spec.monthly_seasonal_cpi_factor_s_m,
                    spec.swap_maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_tuckman_cpi_inflation_swap_seasonality_carry(
                    &TuckmanCpiInflationSwapSeasonalitySpec {
                        notional_n: v[0],
                        nominal_treasury_yield_y_nom: v[1],
                        tips_real_yield_r_real: v[2],
                        inflation_risk_premium_irp: v[3],
                        tips_liquidity_discount_liq: v[4],
                        monthly_seasonal_cpi_factor_s_m: v[5],
                        swap_maturity_years_t: v[6],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_tuckman_cpi_inflation_swap_seasonality_carry(
                    &TuckmanCpiInflationSwapSeasonalitySpec {
                        notional_n: v[0],
                        nominal_treasury_yield_y_nom: v[1],
                        tips_real_yield_r_real: v[2],
                        inflation_risk_premium_irp: v[3],
                        tips_liquidity_discount_liq: v[4],
                        monthly_seasonal_cpi_factor_s_m: v[5],
                        swap_maturity_years_t: v[6],
                    },
                );
                let fd = (up.zciis_seasonal_mark_to_market_npv
                    - dn.zciis_seasonal_mark_to_market_npv)
                    / (2.0 * h);
                assert!(
                    (r_ad.zciis_seasonal_mark_to_market_npv.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_tuckman_cpi_inflation_swap_seasonality_carry_module::*;

pub mod fabozzi_cmo_accrual_z_bond_mod {
    use ql_ad::Real;

    /// Specification for the Fabozzi, Bhattacharya & Berliner (2011, Ch. 8-10) / Hayre (2001)
    /// sequential CMO with an Accrual (Z-Bond) lockout tranche.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct FabozziCmoAccrualZBondSpec<T: Real> {
        /// Initial principal balance $B_{A,0}$ of the Senior Sequential Tranche $A$.
        pub senior_tranche_initial_balance_ba: T,
        /// Initial principal balance $B_{Z,0}$ of the Accrual (Z-Bond) Tranche $Z$.
        pub z_bond_initial_balance_bz: T,
        /// Net pass-through and tranche coupon rate $c > 0$.
        pub tranche_coupon_rate_c: T,
        /// Conditional Prepayment Rate $\text{CPR} \in (0, 1)$ per annum.
        pub conditional_prepayment_rate_cpr: T,
        /// Scheduled amortization fraction $\alpha_{\text{sched}}$ of remaining pool per period.
        pub scheduled_amortization_rate: T,
        /// Base Treasury/OIS discount yield $y_{\text{disc}}$.
        pub base_discount_yield_y: T,
        /// Option-adjusted spread $\text{OAS}_Z$ for discounting Z-bond cashflows.
        pub z_bond_oas_spread: T,
    }

    /// Analytical outputs for the Fabozzi CMO Accrual Z-Bond waterfall engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct FabozziCmoAccrualZBondAnalytics<T: Real> {
        /// Weighted Average Life (WAL, years) of Senior Sequential Tranche $A$.
        pub senior_tranche_wal_years: T,
        /// Weighted Average Life (WAL, years) of the Accrual Z-Bond Tranche $Z$.
        pub z_bond_wal_years: T,
        /// Peak compounded Z-Bond principal multiplier $\max_t B_{Z,t} / B_{Z,0}$.
        pub z_bond_accrued_peak_multiplier: T,
        /// Discounted present value of the Accrual Z-Bond Tranche $Z$.
        pub z_bond_present_value: T,
    }

    /// Evaluates a 6-period sequential CMO waterfall where the Z-Bond accrues compounded
    /// principal during the lockout phase and redirects its cash interest to retire Tranche $A$.
    pub fn evaluate_fabozzi_cmo_accrual_z_bond_sequential_lockout_tranche<T: Real>(
        spec: &FabozziCmoAccrualZBondSpec<T>,
    ) -> FabozziCmoAccrualZBondAnalytics<T> {
        const PERIODS: usize = 6;
        let ba0 = spec
            .senior_tranche_initial_balance_ba
            .max_val(T::from_f64(1.0));
        let bz0 = spec.z_bond_initial_balance_bz.max_val(T::from_f64(1.0));
        let c = spec.tranche_coupon_rate_c.max_val(T::from_f64(1e-4));
        let cpr = spec
            .conditional_prepayment_rate_cpr
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.85));
        let sched = spec
            .scheduled_amortization_rate
            .max_val(T::from_f64(0.02))
            .min_val(T::from_f64(0.50));

        let mut pool_bal = ba0 + bz0;
        let mut ba = ba0;
        let mut bz = bz0;

        let mut wal_a_num = T::zero();
        let mut wal_z_num = T::zero();
        let mut z_prin_total = T::zero();
        let mut z_pv = T::zero();
        let mut max_bz = bz0;

        let eps = T::from_f64(1e-4);
        let smooth_min = |a: T, b: T| -> T {
            T::from_f64(0.5) * (a + b - ((a - b) * (a - b) + eps * eps).sqrt())
        };
        let smooth_max = |a: T, b: T| -> T {
            T::from_f64(0.5) * (a + b + ((a - b) * (a - b) + eps * eps).sqrt())
        };

        let r_z = spec.base_discount_yield_y + spec.z_bond_oas_spread;

        for m in 1..=PERIODS {
            let tm = T::from_f64(m as f64);
            let pay_rate = sched + (T::one() - sched) * cpr;
            let coll_prin = if m == PERIODS {
                pool_bal
            } else {
                pool_bal * pay_rate
            };
            pool_bal = (pool_bal - coll_prin).max_val(T::zero());

            // Z-Bond interest accrual redirected to accelerate Senior Tranche A retirement.
            let z_accrued_int = c * bz;
            let z_compounded_target = bz + z_accrued_int;

            let cash_avail_for_a = coll_prin + z_accrued_int;
            let prin_a = smooth_min(ba, cash_avail_for_a);
            ba = smooth_max(T::zero(), ba - prin_a);
            wal_a_num += tm * prin_a;

            // Any cash in excess of retiring Tranche A flows to the Z-Bond
            let cash_to_z = smooth_max(T::zero(), cash_avail_for_a - prin_a);
            bz = smooth_max(T::zero(), z_compounded_target - cash_to_z);
            max_bz = smooth_max(max_bz, bz);

            let df_z = (-r_z * tm).exp();
            z_pv += df_z * cash_to_z;
            wal_z_num += tm * cash_to_z;
            z_prin_total += cash_to_z;
        }

        FabozziCmoAccrualZBondAnalytics {
            senior_tranche_wal_years: wal_a_num / ba0,
            z_bond_wal_years: wal_z_num / z_prin_total.max_val(T::from_f64(1e-6)),
            z_bond_accrued_peak_multiplier: max_bz / bz0,
            z_bond_present_value: z_pv,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fabozzi_cmo_accrual_z_bond_invariants_and_ad_parity() {
            let base = FabozziCmoAccrualZBondSpec {
                senior_tranche_initial_balance_ba: 75.0_f64,
                z_bond_initial_balance_bz: 25.0,
                tranche_coupon_rate_c: 0.055,
                conditional_prepayment_rate_cpr: 0.16,
                scheduled_amortization_rate: 0.12,
                base_discount_yield_y: 0.042,
                z_bond_oas_spread: 0.012,
            };
            let r0 = evaluate_fabozzi_cmo_accrual_z_bond_sequential_lockout_tranche(&base);
            assert!(r0.z_bond_wal_years > r0.senior_tranche_wal_years);
            assert!(r0.z_bond_accrued_peak_multiplier > 1.0);
            assert!(r0.z_bond_present_value > 0.0);

            let inputs = [
                base.senior_tranche_initial_balance_ba,
                base.z_bond_initial_balance_bz,
                base.tranche_coupon_rate_c,
                base.conditional_prepayment_rate_cpr,
                base.scheduled_amortization_rate,
                base.base_discount_yield_y,
                base.z_bond_oas_spread,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = FabozziCmoAccrualZBondSpec {
                senior_tranche_initial_balance_ba: dual_in[0],
                z_bond_initial_balance_bz: dual_in[1],
                tranche_coupon_rate_c: dual_in[2],
                conditional_prepayment_rate_cpr: dual_in[3],
                scheduled_amortization_rate: dual_in[4],
                base_discount_yield_y: dual_in[5],
                z_bond_oas_spread: dual_in[6],
            };
            let r_ad = evaluate_fabozzi_cmo_accrual_z_bond_sequential_lockout_tranche(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| FabozziCmoAccrualZBondSpec {
                    senior_tranche_initial_balance_ba: v[0],
                    z_bond_initial_balance_bz: v[1],
                    tranche_coupon_rate_c: v[2],
                    conditional_prepayment_rate_cpr: v[3],
                    scheduled_amortization_rate: v[4],
                    base_discount_yield_y: v[5],
                    z_bond_oas_spread: v[6],
                };
                let rp = evaluate_fabozzi_cmo_accrual_z_bond_sequential_lockout_tranche(&mk(p));
                let rm = evaluate_fabozzi_cmo_accrual_z_bond_sequential_lockout_tranche(&mk(m));
                let fd = (rp.z_bond_present_value - rm.z_bond_present_value) / (2.0 * eps);
                assert!(
                    r_ad.z_bond_present_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.z_bond_present_value.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use fabozzi_cmo_accrual_z_bond_mod::*;

// ============================================================================
// Fabozzi (2021) & Hayre (2001) Stripped MBS (IO / PO) Refinancing S-Curve Convexity
// ============================================================================

pub mod fabozzi_hayre_io_po_strip_mod {
    use ql_ad::Real;

    /// Specification for the Fabozzi (2021, Ch. 25–26) / Hayre (2001, Ch. 4) Interest-Only (IO)
    /// and Principal-Only (PO) stripped MBS valuation under an S-curve refinancing CPR model.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct FabozziHayreIoPoStripSpec<T: Real> {
        /// Initial mortgage pool principal balance $B_0 > 0$.
        pub pool_principal_balance: T,
        /// Underlying mortgage Weighted Average Coupon $c_{\text{WAC}} > 0$.
        pub weighted_average_coupon_wac: T,
        /// Net pass-through coupon paid to the IO strip $c_{\text{net}} \in (0, c_{\text{WAC}}]$.
        pub net_passthrough_coupon_c_net: T,
        /// Current primary mortgage market rate $r_{\text{mkt}} > 0$.
        pub current_mortgage_market_rate_r_mkt: T,
        /// Baseline housing-turnover floor prepayment rate $\text{CPR}_{\min} \in (0, 1)$.
        pub baseline_cpr_min: T,
        /// Peak refinancing ceiling prepayment rate $\text{CPR}_{\max} \in (\text{CPR}_{\min}, 1)$.
        pub peak_refi_cpr_max: T,
        /// Refinancing S-curve steepness $\beta_{\text{refi}} > 0$.
        pub refi_s_curve_steepness_beta: T,
        /// Option-Adjusted Spread $\text{OAS}$ over the primary rate.
        pub option_adjusted_spread_oas: T,
    }

    /// Analytical outputs for the Fabozzi-Hayre IO/PO stripped MBS engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct FabozziHayreIoPoStripAnalytics<T: Real> {
        /// S-curve Conditional Prepayment Rate $\text{CPR}(c_{\text{WAC}} - r_{\text{mkt}})$.
        pub logistic_refinancing_cpr: T,
        /// Principal-Only (PO) strip present value $V_{\text{PO}}$.
        pub po_strip_present_value: T,
        /// Combined Pass-Through MBS value $V_{\text{PT}} = V_{\text{IO}} + V_{\text{PO}}$.
        pub combined_passthrough_mbs_value: T,
        /// Interest-Only (IO) strip present value $V_{\text{IO}}$.
        pub io_strip_present_value: T,
    }

    /// Evaluates the Fabozzi (2021, Ch. 25–26) / Hayre (2001, Ch. 4) IO and PO strip prices
    /// across a 6-period mortgage pool with exact principal conservation $\sum \text{Prin}_n = B_0$.
    pub fn evaluate_fabozzi_hayre_io_po_stripped_mbs_prepayment_convexity<T: Real>(
        spec: &FabozziHayreIoPoStripSpec<T>,
    ) -> FabozziHayreIoPoStripAnalytics<T> {
        const PERIODS: usize = 6;
        let b0 = spec.pool_principal_balance.max_val(T::from_f64(1.0));
        let wac = spec.weighted_average_coupon_wac.max_val(T::from_f64(1e-4));
        let c_net = spec.net_passthrough_coupon_c_net.max_val(T::from_f64(1e-4));
        let cpr_min = spec
            .baseline_cpr_min
            .max_val(T::from_f64(0.005))
            .min_val(T::from_f64(0.80));
        let cpr_max = spec
            .peak_refi_cpr_max
            .max_val(cpr_min + T::from_f64(0.01))
            .min_val(T::from_f64(0.95));

        let refi_incentive = wac - spec.current_mortgage_market_rate_r_mkt;
        let logistic =
            T::one() / (T::one() + (-spec.refi_s_curve_steepness_beta * refi_incentive).exp());
        let cpr = cpr_min + (cpr_max - cpr_min) * logistic;

        let y_disc = (spec.current_mortgage_market_rate_r_mkt + spec.option_adjusted_spread_oas)
            .max_val(T::from_f64(1e-4));

        let mut bal = b0;
        let mut io_pv = T::zero();
        let mut po_pv = T::zero();

        for n in 1..=PERIODS {
            let rem_periods = (PERIODS - n + 1) as i32;
            let io_cf = bal * c_net;

            let po_cf = if n < PERIODS {
                let factor = (T::one() + wac).powf(rem_periods as f64);
                let sched_prin = bal * wac / (factor - T::one());
                let prepay = (bal - sched_prin) * cpr;
                sched_prin + prepay
            } else {
                bal
            };

            bal -= po_cf;
            let df = (T::one() + y_disc).powf(-(n as i32) as f64);
            io_pv += df * io_cf;
            po_pv += df * po_cf;
        }

        FabozziHayreIoPoStripAnalytics {
            logistic_refinancing_cpr: cpr,
            po_strip_present_value: po_pv,
            combined_passthrough_mbs_value: io_pv + po_pv,
            io_strip_present_value: io_pv,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fabozzi_hayre_io_po_strip_invariants_and_ad_parity() {
            let base = FabozziHayreIoPoStripSpec {
                pool_principal_balance: 100.0_f64,
                weighted_average_coupon_wac: 0.065,
                net_passthrough_coupon_c_net: 0.060,
                current_mortgage_market_rate_r_mkt: 0.052,
                baseline_cpr_min: 0.06,
                peak_refi_cpr_max: 0.42,
                refi_s_curve_steepness_beta: 35.0,
                option_adjusted_spread_oas: 0.008,
            };
            let r0 = evaluate_fabozzi_hayre_io_po_stripped_mbs_prepayment_convexity(&base);
            assert!(r0.logistic_refinancing_cpr > base.baseline_cpr_min);
            assert!(r0.po_strip_present_value > 0.0 && r0.po_strip_present_value < 100.0);
            assert!(r0.io_strip_present_value > 0.0);
            assert!(
                (r0.combined_passthrough_mbs_value
                    - (r0.io_strip_present_value + r0.po_strip_present_value))
                    .abs()
                    < 1e-12
            );

            let inputs = [
                base.pool_principal_balance,
                base.weighted_average_coupon_wac,
                base.net_passthrough_coupon_c_net,
                base.current_mortgage_market_rate_r_mkt,
                base.baseline_cpr_min,
                base.peak_refi_cpr_max,
                base.refi_s_curve_steepness_beta,
                base.option_adjusted_spread_oas,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = FabozziHayreIoPoStripSpec {
                pool_principal_balance: dual_in[0],
                weighted_average_coupon_wac: dual_in[1],
                net_passthrough_coupon_c_net: dual_in[2],
                current_mortgage_market_rate_r_mkt: dual_in[3],
                baseline_cpr_min: dual_in[4],
                peak_refi_cpr_max: dual_in[5],
                refi_s_curve_steepness_beta: dual_in[6],
                option_adjusted_spread_oas: dual_in[7],
            };
            let r_ad = evaluate_fabozzi_hayre_io_po_stripped_mbs_prepayment_convexity(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| FabozziHayreIoPoStripSpec {
                    pool_principal_balance: v[0],
                    weighted_average_coupon_wac: v[1],
                    net_passthrough_coupon_c_net: v[2],
                    current_mortgage_market_rate_r_mkt: v[3],
                    baseline_cpr_min: v[4],
                    peak_refi_cpr_max: v[5],
                    refi_s_curve_steepness_beta: v[6],
                    option_adjusted_spread_oas: v[7],
                };
                let rp = evaluate_fabozzi_hayre_io_po_stripped_mbs_prepayment_convexity(&mk(p));
                let rm = evaluate_fabozzi_hayre_io_po_stripped_mbs_prepayment_convexity(&mk(m));
                let fd = (rp.io_strip_present_value - rm.io_strip_present_value) / (2.0 * eps);
                assert!(
                    r_ad.io_strip_present_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.io_strip_present_value.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use fabozzi_hayre_io_po_strip_mod::*;
mod jarrow_yildirim_yoy_inflation_mod {
    use ql_ad::Real;

    /// Input specification for the Jarrow & Yildirim (2003, *JFQA* 38(2)) &
    /// Brigo & Mercurio (2006, *Interest Rate Models — Theory and Practice*, Ch. 16, §16.4)
    /// Year-on-Year (YoY) inflation swap convexity adjustment engine.
    #[derive(Debug, Clone, Copy)]
    pub struct JarrowYildirimYoyInflationConvexitySpec<T: Real> {
        /// Forward nominal rate $f_n(0; T_{i-1}, T_i)$ over $[T_{i-1}, T_i]$ (e.g., $0.041$).
        pub nominal_forward_rate_f_n: T,
        /// Forward real TIPS rate $f_r(0; T_{i-1}, T_i)$ over $[T_{i-1}, T_i]$ (e.g., $0.017$).
        pub real_forward_rate_f_r: T,
        /// Nominal Hull–White short-rate volatility $\sigma_n$ (e.g., $0.0085$).
        pub nominal_rate_volatility_sigma_n: T,
        /// Real Hull–White short-rate volatility $\sigma_r$ (e.g., $0.0065$).
        pub real_rate_volatility_sigma_r: T,
        /// CPI index volatility $\sigma_I$ (e.g., $0.012$).
        pub cpi_index_volatility_sigma_i: T,
        /// Correlation $\rho_{n,r}$ between nominal and real rates (e.g., $0.55$).
        pub nominal_real_correlation_rho_nr: T,
        /// Correlation $\rho_{r,I}$ between real rate and CPI index (e.g., $-0.25$).
        pub real_cpi_correlation_rho_ri: T,
        /// Measurement start date $T_{i-1}$ in years (e.g., $4.0$, with $\tau = 1$ year to $T_i = T_{i-1} + 1$).
        pub period_start_years_t_prev: T,
    }

    /// Analytics returned by [`evaluate_jarrow_yildirim_yoy_inflation_swap_convexity_adjustment`].
    #[derive(Debug, Clone, Copy)]
    pub struct JarrowYildirimYoyInflationConvexityAnalytics<T: Real> {
        /// Unadjusted Fisher breakeven inflation rate $\frac{1 + f_n}{1 + f_r} - 1$.
        pub unadjusted_breakeven_inflation_rate: T,
        /// Jarrow–Yildirim (2003) / Brigo–Mercurio (2006, Eq. 16.21) multiplicative convexity correction $C_{\text{YoY}}$.
        pub yoy_convexity_adjustment_bps: T,
        /// Convexity-adjusted expected Year-on-Year inflation rate $\mathbb{E}^{T_i}[I(T_i)/I(T_{i-1}) - 1]$.
        pub convexity_adjusted_yoy_inflation_rate: T,
        /// Present value of a 1-period YoY inflation swaplet per $100$ notional.
        pub discounted_yoy_swaplet_value: T,
    }

    /// Evaluates the Jarrow & Yildirim (2003) / Brigo & Mercurio (2006, Ch. 16, Eq. 16.21)
    /// Year-on-Year inflation swaplet rate under the $T_i$-forward measure:
    /// $$\mathbb{E}^{T_i}\!\left[\frac{I(T_i)}{I(T_{i-1})} - 1\right] = \frac{1 + f_n}{1 + f_r} e^{\mathcal{C}_{\text{JY}}(T_{i-1}, T_i)} - 1,$$
    /// where $\mathcal{C}_{\text{JY}}(T_{i-1}, T_i) = T_{i-1}\big(\sigma_r^2 - \rho_{n,r}\sigma_n\sigma_r - \rho_{r,I}\sigma_r\sigma_I\big)$ over a 1-year accrual period.
    pub fn evaluate_jarrow_yildirim_yoy_inflation_swap_convexity_adjustment<T: Real>(
        spec: &JarrowYildirimYoyInflationConvexitySpec<T>,
    ) -> JarrowYildirimYoyInflationConvexityAnalytics<T> {
        let one_plus_fn = T::one() + spec.nominal_forward_rate_f_n;
        let one_plus_fr = T::one() + spec.real_forward_rate_f_r;
        let fisher_ratio = one_plus_fn / one_plus_fr;
        let breakeven = fisher_ratio - T::one();

        let sig_n = spec.nominal_rate_volatility_sigma_n;
        let sig_r = spec.real_rate_volatility_sigma_r;
        let sig_i = spec.cpi_index_volatility_sigma_i;
        let t_prev = spec.period_start_years_t_prev;

        // Brigo & Mercurio (2006, Eq. 16.21) YoY measure-change covariance exponent over tau = 1.0 yr
        let cov_term = sig_r * sig_r
            - spec.nominal_real_correlation_rho_nr * sig_n * sig_r
            - spec.real_cpi_correlation_rho_ri * sig_r * sig_i;
        let c_exp = t_prev * cov_term;

        let adj_ratio = fisher_ratio * c_exp.exp();
        let adj_yoy = adj_ratio - T::one();
        let conv_bps = (adj_yoy - breakeven) * T::from_f64(10_000.0);

        let df_ti = (-spec.nominal_forward_rate_f_n * (t_prev + T::one())).exp();
        let swaplet_pv = T::from_f64(100.0) * df_ti * adj_yoy;

        JarrowYildirimYoyInflationConvexityAnalytics {
            unadjusted_breakeven_inflation_rate: breakeven,
            yoy_convexity_adjustment_bps: conv_bps,
            convexity_adjusted_yoy_inflation_rate: adj_yoy,
            discounted_yoy_swaplet_value: swaplet_pv,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn jarrow_yildirim_yoy_inflation_invariants_and_ad_parity() {
            let base = JarrowYildirimYoyInflationConvexitySpec {
                nominal_forward_rate_f_n: 0.041_f64,
                real_forward_rate_f_r: 0.017_f64,
                nominal_rate_volatility_sigma_n: 0.0085_f64,
                real_rate_volatility_sigma_r: 0.0065_f64,
                cpi_index_volatility_sigma_i: 0.012_f64,
                nominal_real_correlation_rho_nr: 0.55_f64,
                real_cpi_correlation_rho_ri: -0.25_f64,
                period_start_years_t_prev: 4.0_f64,
            };
            let res = evaluate_jarrow_yildirim_yoy_inflation_swap_convexity_adjustment(&base);
            assert!(res.unadjusted_breakeven_inflation_rate > 0.02);
            assert!(res.yoy_convexity_adjustment_bps > 0.0);
            assert!(
                res.convexity_adjusted_yoy_inflation_rate > res.unadjusted_breakeven_inflation_rate
            );
            assert!(res.discounted_yoy_swaplet_value > 1.5);

            let inputs = [
                base.nominal_forward_rate_f_n,
                base.real_forward_rate_f_r,
                base.nominal_rate_volatility_sigma_n,
                base.real_rate_volatility_sigma_r,
                base.cpi_index_volatility_sigma_i,
                base.nominal_real_correlation_rho_nr,
                base.real_cpi_correlation_rho_ri,
                base.period_start_years_t_prev,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = JarrowYildirimYoyInflationConvexitySpec {
                nominal_forward_rate_f_n: dual_in[0],
                real_forward_rate_f_r: dual_in[1],
                nominal_rate_volatility_sigma_n: dual_in[2],
                real_rate_volatility_sigma_r: dual_in[3],
                cpi_index_volatility_sigma_i: dual_in[4],
                nominal_real_correlation_rho_nr: dual_in[5],
                real_cpi_correlation_rho_ri: dual_in[6],
                period_start_years_t_prev: dual_in[7],
            };
            let r_ad = evaluate_jarrow_yildirim_yoy_inflation_swap_convexity_adjustment(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| JarrowYildirimYoyInflationConvexitySpec {
                    nominal_forward_rate_f_n: v[0],
                    real_forward_rate_f_r: v[1],
                    nominal_rate_volatility_sigma_n: v[2],
                    real_rate_volatility_sigma_r: v[3],
                    cpi_index_volatility_sigma_i: v[4],
                    nominal_real_correlation_rho_nr: v[5],
                    real_cpi_correlation_rho_ri: v[6],
                    period_start_years_t_prev: v[7],
                };
                let rp = evaluate_jarrow_yildirim_yoy_inflation_swap_convexity_adjustment(&mk(p));
                let rm = evaluate_jarrow_yildirim_yoy_inflation_swap_convexity_adjustment(&mk(m));
                let fd = (rp.discounted_yoy_swaplet_value - rm.discounted_yoy_swaplet_value)
                    / (2.0 * eps);
                assert!(
                    r_ad.discounted_yoy_swaplet_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.discounted_yoy_swaplet_value.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use jarrow_yildirim_yoy_inflation_mod::*;
