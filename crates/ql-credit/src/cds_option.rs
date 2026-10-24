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
    let mut prev_surv = credit.survival_probability(te);
    let mut prev_t = te;
    for &t in &s.payment_times {
        if t <= te || t > s.maturity + 1e-12 {
            continue;
        }
        let dt = (t - prev_t).max(0.0);
        let surv = credit.survival_probability(t);
        let df = discount.discount_factor(t);
        prot += (T::one() - s.recovery) * df * (prev_surv - surv);
        ann += df * T::from_f64(dt) * (prev_surv + surv) * 0.5;
        prev_surv = surv;
        prev_t = t;
    }
    let fep = if s.knocks_out_on_default {
        T::zero()
    } else {
        (T::one() - s.recovery)
            * discount.discount_factor(te)
            * (T::one() - credit.survival_probability(te))
    };

    let ann_safe = ann.max_val(T::from_f64(1e-14));
    let eff_fwd = (prot + fep) / ann_safe;
    let npv = black_cds_option_price(
        eff_fwd,
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
        forward_spread: eff_fwd,
        risky_annuity: ann,
        protection_pv: prot,
        fep,
    }
}

#[cfg(test)]
mod phase14_tests {
    use super::*;
    use crate::SurvivalProbabilityCurve;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;
    #[test]
    fn phase14_cds_swaption_curve_backed() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let disc = FlatForwardCurve::new(d, 0.02);
        let credit = SurvivalProbabilityCurve::new(d, vec![1.0, 2.0, 5.0], vec![0.02, 0.02, 0.02]);
        let mut s = CdsSwaptionSpec {
            notional: 100.0,
            strike_spread: 0.012,
            option_expiry: 1.0,
            maturity: 5.0,
            recovery: 0.4,
            running_spread: 0.012,
            vol: 0.25,
            payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            knocks_out_on_default: false,
            is_payer: true,
        };
        let pay = price_cds_swaption_from_curves(&s, &disc, &credit);
        s.is_payer = false;
        let rec = price_cds_swaption_from_curves(&s, &disc, &credit);
        let expected_diff =
            s.notional * (pay.protection_pv + pay.fep - s.strike_spread * pay.risky_annuity);
        assert!((pay.npv - rec.npv - expected_diff).abs() < 1e-10);

        let eval_vol = |v: f64| {
            let mut spec = s.clone();
            spec.is_payer = true;
            spec.vol = v;
            price_cds_swaption_from_curves(&spec, &disc, &credit).npv
        };
        let h = 1e-5;
        let fd_vega = (eval_vol(0.25 + h) - eval_vol(0.25 - h)) / (2.0 * h);
        type D = Dual<1>;
        let disc_d = FlatForwardCurve::new(d, D::from_value(0.02));
        let credit_d = SurvivalProbabilityCurve::new(
            d,
            vec![1.0, 2.0, 5.0],
            vec![
                D::from_value(0.02),
                D::from_value(0.02),
                D::from_value(0.02),
            ],
        );
        let sd = CdsSwaptionSpec {
            notional: D::from_value(100.0),
            strike_spread: D::from_value(0.012),
            option_expiry: 1.0,
            maturity: 5.0,
            recovery: D::from_value(0.4),
            running_spread: D::from_value(0.012),
            vol: D::variable(0.25, 0),
            payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            knocks_out_on_default: false,
            is_payer: true,
        };
        let ad = price_cds_swaption_from_curves(&sd, &disc_d, &credit_d).npv;
        assert!((ad.d[0] - fd_vega).abs() < 1e-6);
    }
}
