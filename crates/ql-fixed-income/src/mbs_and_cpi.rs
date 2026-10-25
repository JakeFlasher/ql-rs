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
        let cf = coupon + sched + pre + (T::one() - spec.loss_severity) * def;
        pv += cf * curve.discount_factor(t) * (-oas * t).exp();
        wal_num += T::from_f64(t) * (sched + pre);
        prep += pre;
        loss += credit_loss;
        bal = (perf - sched - pre).max_val(T::zero());
    }
    SpecifiedPoolValuation {
        pv,
        wal: wal_num / spec.original_balance.max_val(T::from_f64(1e-12)),
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
    let i1 = s.front_clean_price + s.pass_through_rate * (100.0 * s.front_accrued_fraction);
    let i2 = surv * (s.back_clean_price + s.pass_through_rate * (100.0 * s.back_accrued_fraction));
    let cf = (s.pass_through_rate / 12.0 + T::one() - surv) * 100.0;
    let dt = s.roll_period_years.max(1e-6);
    let imp = (cf + i2 - i1) / (i1 * dt);
    let net = i1 * (T::one() + s.repo_rate * dt) - i2 - cf;
    let pstar = (i1 * (T::one() + s.repo_rate * dt) - cf) / surv.max_val(T::from_f64(1e-12))
        - s.pass_through_rate * (100.0 * s.back_accrued_fraction);
    TbaDollarRollReport {
        price_drop: s.front_clean_price - s.back_clean_price,
        breakeven_drop: s.front_clean_price - pstar,
        implied_repo_rate: imp,
        net_roll_advantage: net,
        surviving_factor: surv,
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
            back_clean_price: 98.0,
            wac: 0.04,
            pass_through_rate: 0.035,
            remaining_months: 120,
            front_accrued_fraction: 0.5,
            back_accrued_fraction: 0.5,
            roll_period_years: 1.0,
            assumed_cpr: 0.06,
            repo_rate: 0.02,
        };
        let r = analyze_tba_dollar_roll(&s);
        let mut s2 = s;
        s2.back_clean_price = 99.0 - r.breakeven_drop;
        let r2 = analyze_tba_dollar_roll(&s2);
        assert!(r2.net_roll_advantage.abs() < 1e-10);
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
            cdr: 0.0,
            loss_severity: 0.4,
        };
        assert!(price_specified_pool_mbs(&m, &c, 0.0).pv.is_finite());
        let sd = TbaDollarRollSpec::<Dual<1>> {
            front_clean_price: Dual::variable(99.0, 0),
            back_clean_price: Dual::from_value(98.0),
            wac: Dual::from_value(0.04),
            pass_through_rate: Dual::from_value(0.035),
            remaining_months: 120,
            front_accrued_fraction: 0.5,
            back_accrued_fraction: 0.5,
            roll_period_years: 1.0,
            assumed_cpr: Dual::from_value(0.06),
            repo_rate: Dual::from_value(0.02),
        };
        assert!(analyze_tba_dollar_roll(&sd).implied_repo_rate.d[0].is_finite());
    }
}
