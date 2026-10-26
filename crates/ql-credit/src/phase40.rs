//! Schönbucher (2003) & Lando (2004) Credit-Default Swaption (`Payer`/`Receiver`)
//! Risky-Annuity Numéraire & Front-End Protection (`FEP`) Black-76 Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct SchonbucherCdsSwaptionFepSpec<T: Real = f64> {
    pub hazard_rate_lambda: T,
    pub recovery_rate: T,
    pub risk_free_rate: T,
    pub option_expiry_te: T,
    pub swap_tenor_tm: T,
    pub strike_spread_k: T,
    pub spread_lognormal_vol_sigma: T,
}

#[derive(Debug, Clone, Copy)]
pub struct SchonbucherCdsSwaptionFepReport<T: Real = f64> {
    pub forward_risky_annuity_rpv01: T,
    pub knockout_forward_cds_spread: T,
    pub front_end_protection_pv: T,
    pub fep_adjusted_forward_spread: T,
    pub knockout_payer_swaption_bps: T,
    pub non_knockout_payer_swaption_bps: T,
    pub receiver_swaption_bps: T,
    pub risky_put_call_parity_error: T,
}

pub fn evaluate_schonbucher_cds_swaption_and_fep<T: Real>(
    spec: &SchonbucherCdsSwaptionFepSpec<T>,
) -> SchonbucherCdsSwaptionFepReport<T> {
    let lam = spec.hazard_rate_lambda.max_val(T::from_f64(1e-5));
    let r = spec.risk_free_rate.max_val(T::from_f64(1e-4));
    let lgd = T::one() - spec.recovery_rate;
    let te = spec.option_expiry_te.max_val(T::from_f64(0.05));
    let tm = spec.swap_tenor_tm.max_val(T::from_f64(0.5));
    let k = spec.strike_spread_k.max_val(T::from_f64(1e-5));
    let sig = spec.spread_lognormal_vol_sigma.max_val(T::from_f64(1e-4));

    // Forward survival-contingent risky annuity A(0; T_e, T_e + T_M):
    let u = r + lam;
    let df_surv_te = (-u * te).exp();
    let fwd_rpv01 = df_surv_te * (T::one() - (-u * tm).exp()) / u;

    // Knockout forward CDS par spread S_fwd = (1 - R) * lambda:
    let s_fwd = lgd * lam;

    // Front-End Protection (FEP) PV for default prior to option expiry T_e:
    let fep_pv = lgd * lam * (T::one() - df_surv_te) / u;
    let s_fep = s_fwd + fep_pv / fwd_rpv01.max_val(T::from_f64(1e-8));

    // Black-76 under risky annuity numéraire for knockout payer/receiver swaption:
    let std = sig * te.sqrt();
    let d1_ko = ((s_fwd / k).ln() + T::from_f64(0.5) * std * std) / std;
    let d2_ko = d1_ko - std;
    let payer_ko = fwd_rpv01 * (s_fwd * norm_cdf(d1_ko) - k * norm_cdf(d2_ko));
    let receiver_ko = fwd_rpv01 * (k * norm_cdf(-d2_ko) - s_fwd * norm_cdf(-d1_ko));

    // Non-knockout payer swaption using FEP-adjusted forward spread S_FEP:
    let d1_fep = ((s_fep / k).ln() + T::from_f64(0.5) * std * std) / std;
    let d2_fep = d1_fep - std;
    let payer_fep = fwd_rpv01 * (s_fep * norm_cdf(d1_fep) - k * norm_cdf(d2_fep));

    // Risky-annuity Put-Call Parity: Payer_KO - Receiver_KO == A(0; T_e, T_M) * (S_fwd - K)
    let pcp_err = ((payer_ko - receiver_ko) - fwd_rpv01 * (s_fwd - k)).abs();

    SchonbucherCdsSwaptionFepReport {
        forward_risky_annuity_rpv01: fwd_rpv01,
        knockout_forward_cds_spread: s_fwd,
        front_end_protection_pv: fep_pv,
        fep_adjusted_forward_spread: s_fep,
        knockout_payer_swaption_bps: payer_ko * T::from_f64(1e4),
        non_knockout_payer_swaption_bps: payer_fep * T::from_f64(1e4),
        receiver_swaption_bps: receiver_ko * T::from_f64(1e4),
        risky_put_call_parity_error: pcp_err,
    }
}

#[cfg(test)]
mod phase40_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn schonbucher_cds_swaption_fep_invariants_and_ad_parity() {
        let spec = SchonbucherCdsSwaptionFepSpec {
            hazard_rate_lambda: 0.022,
            recovery_rate: 0.40,
            risk_free_rate: 0.035,
            option_expiry_te: 1.0,
            swap_tenor_tm: 5.0,
            strike_spread_k: 0.013,
            spread_lognormal_vol_sigma: 0.45,
        };
        let rep = evaluate_schonbucher_cds_swaption_and_fep(&spec);
        assert!(rep.risky_put_call_parity_error < 1e-12);
        assert!(rep.non_knockout_payer_swaption_bps > rep.knockout_payer_swaption_bps);

        let h = 1e-7;
        let r_ad = evaluate_schonbucher_cds_swaption_and_fep(&SchonbucherCdsSwaptionFepSpec {
            hazard_rate_lambda: Dual::<1>::variable(spec.hazard_rate_lambda, 0),
            recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
            risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
            option_expiry_te: Dual::<1>::from_value(spec.option_expiry_te),
            swap_tenor_tm: Dual::<1>::from_value(spec.swap_tenor_tm),
            strike_spread_k: Dual::<1>::from_value(spec.strike_spread_k),
            spread_lognormal_vol_sigma: Dual::<1>::from_value(spec.spread_lognormal_vol_sigma),
        });
        let mut sp = spec;
        sp.hazard_rate_lambda += h;
        let mut sm = spec;
        sm.hazard_rate_lambda -= h;
        let fd = (evaluate_schonbucher_cds_swaption_and_fep(&sp).non_knockout_payer_swaption_bps
            - evaluate_schonbucher_cds_swaption_and_fep(&sm).non_knockout_payer_swaption_bps)
            / (2.0 * h);
        assert!((r_ad.non_knockout_payer_swaption_bps.d[0] - fd).abs() < 1e-5);
    }
}
