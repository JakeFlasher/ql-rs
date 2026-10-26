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
    let gross = (lower - upper).abs();
    let after = gross * (T::one() - s.tax_rate_section_163j);
    let unhedged = (s.acquisition_offer_price / s.base_conversion_price * ratio).max_val(T::zero());
    let hedged = (unhedged - gross / s.acquisition_offer_price.max_val(T::from_f64(1e-8)))
        .max_val(T::zero());
    let intrinsic =
        (s.acquisition_offer_price - s.base_conversion_price).max_val(T::zero()) * ratio;
    let mw = (ratio
        * call(
            s.acquisition_offer_price,
            s.base_conversion_price,
            s.volatility,
            s.risk_free_rate,
            s.remaining_years_at_acquisition,
        )
        - intrinsic)
        .abs()
        / s.acquisition_offer_price.max_val(T::from_f64(1e-8));
    MayerBrownSsgaCappedCallMakeWholeReport {
        base_conversion_ratio: ratio,
        lower_call_option_value: lower,
        upper_cap_call_option_value: upper,
        capped_call_overlay_gross_cost: gross,
        after_tax_integrated_capped_call_cost: after,
        unhedged_net_dilution_shares_at_offer: unhedged,
        capped_call_hedged_dilution_shares_at_offer: hedged,
        make_whole_additional_shares: mw,
        total_fundamental_change_conversion_value: intrinsic,
    }
}
#[cfg(test)]
mod phase32_tests {
    use super::*;
    #[test]
    fn overlay_positive() {
        let s = MayerBrownSsgaCappedCallMakeWholeSpec {
            bond_par: 100.0,
            spot_price: 100.0,
            base_conversion_price: 90.0,
            cap_strike_price: 120.0,
            volatility: 0.2,
            risk_free_rate: 0.03,
            maturity_years: 5.0,
            acquisition_offer_price: 110.0,
            remaining_years_at_acquisition: 2.0,
            tax_rate_section_163j: 0.2,
        };
        assert!(
            evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(&s)
                .capped_call_overlay_gross_cost
                > 0.0
        );
    }
}
