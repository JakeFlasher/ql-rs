use ql_ad::Real;
#[derive(Clone, Debug)]
pub struct BrigoMoriniCmcdsCreditOptionSpec<T: Real = f64> {
    pub reset_times: Vec<T>,
    pub forward_cds_spreads: Vec<T>,
    pub spread_volatility: T,
    pub annuity_convexity_gamma: T,
    pub flat_par_cds_spread: T,
    pub option_expiry: T,
    pub option_strike_spread: T,
    pub risky_pv01_annuity: T,
    pub hazard_rate: T,
    pub recovery_rate: T,
    pub risk_free_rate: T,
}
#[derive(Clone, Copy, Debug)]
pub struct BrigoMoriniCmcdsCreditOptionReport<T: Real = f64> {
    pub cmcds_convexity_uplift_bps: T,
    pub cmcds_equilibrium_participation_rate: T,
    pub front_end_protection_pv: T,
    pub fep_adjusted_forward_spread: T,
    pub fep_adjusted_strike_spread: T,
    pub morini_brigo_payer_swaption_pv: T,
    pub naive_unadjusted_swaption_pv: T,
}
pub fn evaluate_brigo_morini_cmcds_and_credit_index_option<T: Real>(
    s: &BrigoMoriniCmcdsCreditOptionSpec<T>,
) -> BrigoMoriniCmcdsCreditOptionReport<T> {
    let mut den = T::zero();
    let mut num = T::zero();
    let mut up = T::zero();
    for i in 0..s.reset_times.len().min(s.forward_cds_spreads.len()) {
        let t = s.reset_times[i];
        let u = (s.annuity_convexity_gamma * s.spread_volatility * s.spread_volatility * t).exp();
        den += s.forward_cds_spreads[i] * u;
        num += s.flat_par_cds_spread;
        up += s.forward_cds_spreads[i] * (u - T::one());
    }
    let part = (num / den.max_val(T::from_f64(1e-8))).min_val(T::one());
    let df = (-s.risk_free_rate * s.option_expiry).exp();
    let surv = (-s.hazard_rate * s.option_expiry).exp();
    let fep = (T::one() - s.recovery_rate) * (T::one() - surv) * df;
    let adjf = s.flat_par_cds_spread + fep / s.risky_pv01_annuity.max_val(T::from_f64(1e-8));
    let adjk = s.option_strike_spread + fep / s.risky_pv01_annuity.max_val(T::from_f64(1e-8));
    let vol = s.spread_volatility * s.option_expiry.sqrt();
    let d1 = ((adjf / adjk.max_val(T::from_f64(1e-8))).ln() + T::from_f64(0.5) * vol * vol)
        / vol.max_val(T::from_f64(1e-8));
    let d2 = d1 - vol;
    let pv = s.risky_pv01_annuity * (adjf * d1.norm_cdf() - adjk * d2.norm_cdf()) * df;
    let naive = s.risky_pv01_annuity
        * (s.flat_par_cds_spread * d1.norm_cdf() - s.option_strike_spread * d2.norm_cdf())
        * df;
    BrigoMoriniCmcdsCreditOptionReport {
        cmcds_convexity_uplift_bps: up / den.max_val(T::from_f64(1e-8)) * T::from_f64(1e4),
        cmcds_equilibrium_participation_rate: part,
        front_end_protection_pv: fep,
        fep_adjusted_forward_spread: adjf,
        fep_adjusted_strike_spread: adjk,
        morini_brigo_payer_swaption_pv: pv,
        naive_unadjusted_swaption_pv: naive,
    }
}
#[cfg(test)]
mod phase32_tests {
    use super::*;
    #[test]
    fn fep_positive() {
        let s = BrigoMoriniCmcdsCreditOptionSpec {
            reset_times: vec![1.0],
            forward_cds_spreads: vec![0.01],
            spread_volatility: 0.2,
            annuity_convexity_gamma: 0.5,
            flat_par_cds_spread: 0.01,
            option_expiry: 1.0,
            option_strike_spread: 0.01,
            risky_pv01_annuity: 4.0,
            hazard_rate: 0.02,
            recovery_rate: 0.4,
            risk_free_rate: 0.02,
        };
        assert!(
            evaluate_brigo_morini_cmcds_and_credit_index_option(&s).fep_adjusted_forward_spread
                > 0.01
        );
    }
}
