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
        let idx = (t * 12.0).round() as u8;
        let ratio =
            inflation_curve.forward_cpi(t, Some((idx.saturating_sub(1) % 12) + 1)) / spec.base_cpi;
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
    let t = *spec.payment_times.last().unwrap_or(&1.0);
    let target = nominal_price / (spec.face_amount * nominal_discount.discount_factor(t));
    target.powf(1.0 / t) - T::one()
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
    };
    T::one() - (T::one() - cpr).powf(1.0 / 12.0)
}

pub fn price_agency_mbs<T: Real>(
    spec: &AgencyMbsSpec<T>,
    curve: &impl YieldTermStructure<T>,
    oas: T,
) -> (T, T) {
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
