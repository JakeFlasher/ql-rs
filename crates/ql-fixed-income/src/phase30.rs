use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct AnDaoquanCliffordChanceResettableCbSpec<T: Real = f64> {
    pub spot_price: T,
    pub initial_conversion_price: T,
    pub book_value_per_share_bps: T,
    pub avg_20day_price: T,
    pub face_value: T,
    pub coupon_rate: T,
    pub redemption_maturity_price: T,
    pub put_back_trigger_ratio: T,
    pub put_back_price: T,
    pub soft_call_trigger_ratio: T,
    pub soft_call_price: T,
    pub stock_volatility: T,
    pub risk_free_rate: T,
    pub credit_spread: T,
    pub maturity: T,
}
#[derive(Clone, Copy, Debug)]
pub struct AnDaoquanCliffordChanceResettableCbReport<T: Real = f64> {
    pub straight_bond_floor: T,
    pub initial_conversion_parity: T,
    pub legally_floored_reset_conversion_price: T,
    pub post_reset_conversion_parity: T,
    pub downward_reset_option_uplift: T,
    pub conditional_put_back_protection_value: T,
    pub issuer_soft_call_cap_giveup: T,
    pub three_line_duplex_fair_value: T,
}
pub fn evaluate_an_daoquan_clifford_chance_resettable_convertible<T: Real>(
    s: &AnDaoquanCliffordChanceResettableCbSpec<T>,
) -> AnDaoquanCliffordChanceResettableCbReport<T> {
    let floor = s.face_value
        * (-s.credit_spread * s.maturity).exp()
        * (T::one() + s.coupon_rate * s.maturity * T::from_f64(0.5));
    let init_price = s.initial_conversion_price.max_val(T::from_f64(1e-6));
    let initial = s.face_value / init_price * s.spot_price;
    let candidate = s.avg_20day_price.max_val(s.book_value_per_share_bps);
    let post_price = init_price.min_val(candidate).max_val(T::from_f64(1e-6));
    let post = s.face_value / post_price * s.spot_price;
    let trigger_gap = s.initial_conversion_price * s.put_back_trigger_ratio - s.avg_20day_price;
    let reset_weight = T::one() / (T::one() + (-trigger_gap * T::from_f64(0.2)).exp());
    let uplift = (post - initial).max_val(T::zero())
        * (T::from_f64(0.5) + T::from_f64(0.5) * reset_weight)
        + s.stock_volatility * s.maturity.sqrt() * T::from_f64(4.0);
    let put = (s.put_back_price - initial.min_val(floor)).max_val(T::zero())
        * (-s.risk_free_rate * s.maturity).exp()
        * (T::from_f64(0.25) + T::from_f64(0.5) * reset_weight);
    let cap = (post - s.initial_conversion_price * s.soft_call_trigger_ratio).max_val(T::zero())
        * T::from_f64(0.2);
    let fair = floor.max_val(initial) + uplift + put - cap;
    AnDaoquanCliffordChanceResettableCbReport {
        straight_bond_floor: floor,
        initial_conversion_parity: initial,
        legally_floored_reset_conversion_price: post_price,
        post_reset_conversion_parity: post,
        downward_reset_option_uplift: uplift,
        conditional_put_back_protection_value: put,
        issuer_soft_call_cap_giveup: cap,
        three_line_duplex_fair_value: fair,
    }
}

#[cfg(test)]
mod phase30_tests {
    use super::*;
    use ql_ad::Dual;
    fn make<T: Real>(x: T) -> AnDaoquanCliffordChanceResettableCbSpec<T> {
        AnDaoquanCliffordChanceResettableCbSpec {
            spot_price: T::from_f64(90.0) + x,
            initial_conversion_price: T::from_f64(100.0),
            book_value_per_share_bps: T::from_f64(80.0),
            avg_20day_price: T::from_f64(85.0),
            face_value: T::from_f64(100.0),
            coupon_rate: T::from_f64(0.02),
            redemption_maturity_price: T::from_f64(105.0),
            put_back_trigger_ratio: T::from_f64(0.7),
            put_back_price: T::from_f64(103.0),
            soft_call_trigger_ratio: T::from_f64(1.3),
            soft_call_price: T::from_f64(105.0),
            stock_volatility: T::from_f64(0.2),
            risk_free_rate: T::from_f64(0.03),
            credit_spread: T::from_f64(0.01),
            maturity: T::from_f64(3.0),
        }
    }
    #[test]
    fn phase30_reset_cb_dual_fd() {
        let r = evaluate_an_daoquan_clifford_chance_resettable_convertible(&make(0.0));
        assert!(r.legally_floored_reset_conversion_price >= 80.0);
        assert!(r.post_reset_conversion_parity > r.initial_conversion_parity);
        let h = 1e-6;
        let fd = (evaluate_an_daoquan_clifford_chance_resettable_convertible(&make(h))
            .post_reset_conversion_parity
            - evaluate_an_daoquan_clifford_chance_resettable_convertible(&make(-h))
                .post_reset_conversion_parity)
            / (2.0 * h);
        let ad = evaluate_an_daoquan_clifford_chance_resettable_convertible(&make(
            Dual::<1>::variable(0.0, 0),
        ))
        .post_reset_conversion_parity;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
