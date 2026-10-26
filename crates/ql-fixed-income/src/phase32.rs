use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct MayerBrownSsgaCappedCallMakeWholeSpec<T: Real = f64> {
    pub bond_par: T,
    pub spot_price: T,
    pub base_conversion_price: T,
    pub cap_strike_price: T,
    pub volatility: T,
    pub risk_free_rate: T,
    pub maturity_years: T,
    pub acquisition_offer_price: T,
    pub remaining_years_at_acquisition: T,
    pub tax_rate_section_163j: T,
}
#[derive(Clone, Copy, Debug)]
pub struct MayerBrownSsgaCappedCallMakeWholeReport<T: Real = f64> {
    pub base_conversion_ratio: T,
    pub lower_call_option_value: T,
    pub upper_cap_call_option_value: T,
    pub capped_call_overlay_gross_cost: T,
    pub after_tax_integrated_capped_call_cost: T,
    pub unhedged_net_dilution_shares_at_offer: T,
    pub capped_call_hedged_dilution_shares_at_offer: T,
    pub make_whole_additional_shares: T,
    pub total_fundamental_change_conversion_value: T,
}
fn call<T: Real>(s: T, k: T, v: T, r: T, t: T) -> T {
    let vol = (v * t.sqrt()).max_val(T::from_f64(1e-8));
    let d1 = ((s / k).ln() + (r + T::from_f64(0.5) * v * v) * t) / vol;
    let d2 = d1 - vol;
    s * d1.norm_cdf() - k * (-r * t).exp() * d2.norm_cdf()
}
pub fn evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole<T: Real>(
    s: &MayerBrownSsgaCappedCallMakeWholeSpec<T>,
) -> MayerBrownSsgaCappedCallMakeWholeReport<T> {
    let ratio = s.bond_par / s.base_conversion_price.max_val(T::from_f64(1e-8));
    let lower = ratio
        * call(
            s.spot_price,
            s.base_conversion_price,
            s.volatility,
            s.risk_free_rate,
            s.maturity_years,
        );
    let upper = ratio
        * call(
            s.spot_price,
            s.cap_strike_price,
            s.volatility,
            s.risk_free_rate,
            s.maturity_years,
        );
    let gross = (lower - upper).max_val(T::zero());
    let after = gross * (T::one() - s.tax_rate_section_163j);
    let s_offer_safe = s.acquisition_offer_price.max_val(T::from_f64(1e-8));
    let unhedged = ratio
        * ((s.acquisition_offer_price - s.base_conversion_price).max_val(T::zero()) / s_offer_safe);
    let hedged = ratio
        * ((s.acquisition_offer_price - s.cap_strike_price).max_val(T::zero()) / s_offer_safe);
    let intrinsic =
        (s.acquisition_offer_price - s.base_conversion_price).max_val(T::zero()) * ratio;
    let call_at_offer = ratio
        * call(
            s.acquisition_offer_price,
            s.base_conversion_price,
            s.volatility,
            s.risk_free_rate,
            s.remaining_years_at_acquisition,
        );
    let mw = (call_at_offer - intrinsic).max_val(T::zero()) / s_offer_safe;
    MayerBrownSsgaCappedCallMakeWholeReport {
        base_conversion_ratio: ratio,
        lower_call_option_value: lower,
        upper_cap_call_option_value: upper,
        capped_call_overlay_gross_cost: gross,
        after_tax_integrated_capped_call_cost: after,
        unhedged_net_dilution_shares_at_offer: unhedged,
        capped_call_hedged_dilution_shares_at_offer: hedged,
        make_whole_additional_shares: mw,
        total_fundamental_change_conversion_value: (ratio + mw) * s.acquisition_offer_price,
    }
}
#[cfg(test)]
mod phase32_tests {
    use super::*;
    use ql_ad::Dual;

    fn spec<T: Real>(x: T) -> MayerBrownSsgaCappedCallMakeWholeSpec<T> {
        MayerBrownSsgaCappedCallMakeWholeSpec {
            bond_par: T::from_f64(100.0),
            spot_price: T::from_f64(100.0) + x,
            base_conversion_price: T::from_f64(90.0),
            cap_strike_price: T::from_f64(120.0),
            volatility: T::from_f64(0.2),
            risk_free_rate: T::from_f64(0.03),
            maturity_years: T::from_f64(5.0),
            acquisition_offer_price: T::from_f64(110.0),
            remaining_years_at_acquisition: T::from_f64(2.0),
            tax_rate_section_163j: T::from_f64(0.2),
        }
    }

    #[test]
    fn overlay_positive() {
        let r = evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(&spec(0.0));
        assert!(r.capped_call_overlay_gross_cost > 0.0);
        assert!(
            r.capped_call_hedged_dilution_shares_at_offer
                <= r.unhedged_net_dilution_shares_at_offer
        );
        assert!(r.make_whole_additional_shares > 0.0);
        let h = 1e-6;
        let fd = (evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(&spec(h))
            .capped_call_overlay_gross_cost
            - evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(&spec(-h))
                .capped_call_overlay_gross_cost)
            / (2.0 * h);
        let ad = evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(&spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .capped_call_overlay_gross_cost;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
