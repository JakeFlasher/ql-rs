//! Black's option on a forward CDS spread, including front-end protection.

use ql_ad::Real;

/// Black (1976) CDS option.  `risky_annuity` is the PV of one unit of
/// running spread, `default_prob_to_expiry` is the cumulative default
/// probability to option expiry, and `discount_to_expiry` is its discount
/// factor.  A payer option is a call on the spread; a receiver is a put.
#[allow(clippy::too_many_arguments)]
pub fn black_cds_option_price<T: Real>(
    forward_spread: T,
    strike_spread: T,
    risky_annuity: T,
    vol: T,
    expiry: T,
    default_prob_to_expiry: T,
    discount_to_expiry: T,
    recovery: T,
    notional: T,
    is_payer: bool,
    knocks_out: bool,
) -> T {
    let sqrt_t = expiry.max_val(T::zero()).sqrt();
    let std = vol.abs() * sqrt_t;
    let intrinsic = if is_payer {
        (forward_spread - strike_spread).max_val(T::zero())
    } else {
        (strike_spread - forward_spread).max_val(T::zero())
    };
    let option = if std.val() < 1e-14 || forward_spread.val() <= 0.0 || strike_spread.val() <= 0.0 {
        risky_annuity * intrinsic
    } else {
        let d1 = ((forward_spread / strike_spread).ln() + std * std * 0.5) / std;
        let d2 = d1 - std;
        if is_payer {
            risky_annuity * (forward_spread * d1.norm_cdf() - strike_spread * d2.norm_cdf())
        } else {
            risky_annuity * (strike_spread * (-d2).norm_cdf() - forward_spread * (-d1).norm_cdf())
        }
    } * notional;
    let fep = if !knocks_out && is_payer {
        notional * (T::one() - recovery) * default_prob_to_expiry * discount_to_expiry
    } else {
        T::zero()
    };
    option + fep
}

#[derive(Clone, Debug)]
pub struct CdsSwaptionSpec<T: Real = f64> {
    pub notional: T,
    pub strike_spread: T,
    pub option_expiry: f64,
    pub maturity: f64,
    pub recovery: T,
    pub running_spread: T,
    pub vol: T,
    pub payment_times: Vec<f64>,
    pub knocks_out_on_default: bool,
    pub is_payer: bool,
}
#[derive(Copy, Clone, Debug)]
pub struct CdsSwaptionValuation<T: Real> {
    pub npv: T,
    pub forward_spread: T,
    pub risky_annuity: T,
    pub protection_pv: T,
    pub fep: T,
}
pub fn price_cds_swaption_from_curves<T: Real>(
    s: &CdsSwaptionSpec<T>,
    discount: &impl ql_curves::YieldTermStructure<T>,
    credit: &crate::SurvivalProbabilityCurve<T>,
) -> CdsSwaptionValuation<T> {
    let te = s.option_expiry;
    let mut prot = T::zero();
    let mut ann = T::zero();
    let mut prev = credit.survival_probability(te);
    for &t in &s.payment_times {
        if t <= te || t > s.maturity {
            continue;
        }
        let surv = credit.survival_probability(t);
        let df = discount.discount_factor(t);
        prot += (T::one() - s.recovery) * df * (prev - surv);
        ann += df * (t - te) * surv;
        prev = surv;
    }
    let fep = if s.knocks_out_on_default {
        T::zero()
    } else {
        (T::one() - s.recovery)
            * discount.discount_factor(te)
            * (T::one() - credit.survival_probability(te))
    };
    let fwd = (prot + fep) / ann.max_val(T::from_f64(1e-14));
    let npv = black_cds_option_price(
        fwd,
        s.strike_spread,
        ann,
        s.vol,
        T::from_f64(te),
        T::one() - credit.survival_probability(te),
        discount.discount_factor(te),
        s.recovery,
        s.notional,
        s.is_payer,
        true,
    );
    CdsSwaptionValuation {
        npv,
        forward_spread: fwd,
        risky_annuity: ann,
        protection_pv: prot,
        fep,
    }
}
