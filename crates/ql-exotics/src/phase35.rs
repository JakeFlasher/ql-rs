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

#[derive(Clone, Debug)]
pub struct BouzoubaaNapoleonReverseCliquetSpec<T: Real = f64> {
    pub period_returns: Vec<T>,
    pub napoleon_fixed_coupon: T,
    pub global_floor: T,
    pub local_cap: T,
    pub local_floor: T,
    pub forward_volatility: T,
    pub forward_skew: T,
    pub risk_free_rate: T,
    pub maturity_years: T,
}

#[derive(Clone, Copy, Debug)]
pub struct BouzoubaaNapoleonReverseCliquetReport<T: Real = f64> {
    pub worst_negative_period_return: T,
    pub cumulative_negative_return_sum: T,
    pub napoleon_coupon_payoff: T,
    pub reverse_cliquet_coupon_payoff: T,
    pub standard_capped_cliquet_payoff: T,
    pub forward_vol_cliquet_convexity_cost: T,
    pub napoleon_capital_protected_note_pv: T,
}

pub fn evaluate_bouzoubaa_napoleon_and_reverse_cliquet<T: Real>(
    s: &BouzoubaaNapoleonReverseCliquetSpec<T>,
) -> BouzoubaaNapoleonReverseCliquetReport<T> {
    let mut w_min = T::zero();
    let mut neg_sum = T::zero();
    let mut std_sum = T::zero();
    for &rm in &s.period_returns {
        let neg_m = smin(rm, T::zero());
        w_min = smin(w_min, neg_m);
        neg_sum += neg_m;
        std_sum += smax(s.local_floor, smin(s.local_cap, rm));
    }
    let c_nap = smax(s.global_floor, s.napoleon_fixed_coupon + w_min);
    let c_rev = smax(s.global_floor, s.napoleon_fixed_coupon + neg_sum);
    let c_std = smax(s.global_floor, std_sum);
    let m = (s.period_returns.len().max(1)) as f64;
    let cost = T::from_f64(0.5)
        * s.forward_volatility
        * (s.maturity_years / T::from_f64(m)).sqrt()
        * T::from_f64((m + 1.0).ln())
        * (T::one() - s.forward_skew);
    let pv = (-s.risk_free_rate * s.maturity_years).exp() * (T::one() + c_nap - cost);
    BouzoubaaNapoleonReverseCliquetReport {
        worst_negative_period_return: w_min,
        cumulative_negative_return_sum: neg_sum,
        napoleon_coupon_payoff: c_nap,
        reverse_cliquet_coupon_payoff: c_rev,
        standard_capped_cliquet_payoff: c_std,
        forward_vol_cliquet_convexity_cost: cost,
        napoleon_capital_protected_note_pv: pv,
    }
}

#[cfg(test)]
mod phase35_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> BouzoubaaNapoleonReverseCliquetSpec<T> {
        BouzoubaaNapoleonReverseCliquetSpec {
            period_returns: vec![
                T::from_f64(0.03),
                T::from_f64(-0.04) + x,
                T::from_f64(0.02),
                T::from_f64(-0.025),
            ],
            napoleon_fixed_coupon: T::from_f64(0.08),
            global_floor: T::zero(),
            local_cap: T::from_f64(0.04),
            local_floor: T::from_f64(-0.04),
            forward_volatility: T::from_f64(0.20),
            forward_skew: T::from_f64(-0.15),
            risk_free_rate: T::from_f64(0.03),
            maturity_years: T::from_f64(2.0),
        }
    }
    #[test]
    fn phase35_exotics_invariants_and_dual_fd() {
        let r = evaluate_bouzoubaa_napoleon_and_reverse_cliquet(&spec(0.0));
        assert!(r.napoleon_coupon_payoff >= 0.0);
        assert!(r.napoleon_coupon_payoff >= r.reverse_cliquet_coupon_payoff - 1e-8);
        let h = 1e-6;
        let fd = (evaluate_bouzoubaa_napoleon_and_reverse_cliquet(&spec(h))
            .napoleon_capital_protected_note_pv
            - evaluate_bouzoubaa_napoleon_and_reverse_cliquet(&spec(-h))
                .napoleon_capital_protected_note_pv)
            / (2.0 * h);
        let ad =
            evaluate_bouzoubaa_napoleon_and_reverse_cliquet(&spec(Dual::<1>::variable(0.0, 0)))
                .napoleon_capital_protected_note_pv;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
