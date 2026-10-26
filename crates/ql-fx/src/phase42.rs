//! Wystup (2017) & Malz (1997) FX 25-Delta Risk-Reversal (`RR25`) & Butterfly (`BF25`)
//! Broker-Fly vs. Smile-Strangle Strike Inversion Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct WystupMalzFxSmileInversionSpec<T: Real = f64> {
    pub spot_fx_s0: T,
    pub domestic_rate_rd: T,
    pub foreign_rate_rf: T,
    pub atm_dns_vol_sigma_atm: T,
    pub risk_reversal_25d_rr: T,
    pub butterfly_25d_bf: T,
    pub expiry_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct WystupMalzFxSmileInversionReport<T: Real = f64> {
    pub forward_fx_rate: T,
    pub atm_delta_neutral_strike_k_dns: T,
    pub implied_25d_call_vol: T,
    pub implied_25d_put_vol: T,
    pub inverted_25d_call_strike_k25c: T,
    pub inverted_25d_put_strike_k25p: T,
    pub smile_strangle_minus_broker_fly_bps: T,
    pub call_25d_target_delta_error: T,
}

pub fn evaluate_wystup_malz_fx_broker_fly_smile_inversion<T: Real>(
    spec: &WystupMalzFxSmileInversionSpec<T>,
) -> WystupMalzFxSmileInversionReport<T> {
    let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-4));
    let rd = spec.domestic_rate_rd;
    let rf = spec.foreign_rate_rf;
    let t = spec.expiry_years_t.max_val(T::from_f64(0.02));
    let sqrt_t = t.sqrt();
    let df_f = (-rf * t).exp();

    let fwd = s0 * ((rd - rf) * t).exp();
    let sig_atm = spec.atm_dns_vol_sigma_atm.max_val(T::from_f64(1e-4));
    let rr = spec.risk_reversal_25d_rr;
    let bf = spec.butterfly_25d_bf;

    // Unadjusted delta-neutral straddle (DNS) ATM strike:
    let k_dns = fwd * (T::from_f64(0.5) * sig_atm * sig_atm * t).exp();

    // 25-delta call and put implied volatilities from ATM, 25D RR, and 25D BF:
    let sig_25c = (sig_atm + bf + T::from_f64(0.5) * rr).max_val(T::from_f64(1e-4));
    let sig_25p = (sig_atm + bf - T::from_f64(0.5) * rr).max_val(T::from_f64(1e-4));

    // Exact spot-delta inversion for 25D call and 25D put strikes:
    // Spot Call Delta = e^{-r_f T} Phi(d_1) = 0.25 => d_1 = Phi^{-1}(0.25 * e^{r_f T})
    // Solve via 6 Newton steps in d_1 starting from -0.67448975:
    let target_cdf = (T::from_f64(0.25) / df_f)
        .max_val(T::from_f64(0.05))
        .min_val(T::from_f64(0.95));
    let mut alpha_c = T::from_f64(0.6744897501960817);
    for _ in 0..6 {
        let cdf_val = norm_cdf(-alpha_c);
        let pdf_val =
            (-T::from_f64(0.5) * alpha_c * alpha_c).exp() * T::from_f64(0.3989422804014327);
        alpha_c += (cdf_val - target_cdf) / pdf_val.max_val(T::from_f64(1e-6));
    }

    let k_25c = fwd * (alpha_c * sig_25c * sqrt_t + T::from_f64(0.5) * sig_25c * sig_25c * t).exp();
    let k_25p =
        fwd * (-alpha_c * sig_25p * sqrt_t + T::from_f64(0.5) * sig_25p * sig_25p * t).exp();

    let d1_c = ((fwd / k_25c).ln() + T::from_f64(0.5) * sig_25c * sig_25c * t) / (sig_25c * sqrt_t);
    let delta_err = (df_f * norm_cdf(d1_c) - T::from_f64(0.25)).abs();

    // Wystup (2017) Smile-Strangle vs. Broker-Fly convexity adjustment in bps:
    let fly_conv_bps =
        T::from_f64(0.125) * (alpha_c * alpha_c) * (rr * rr) * sqrt_t / sig_atm * T::from_f64(1e4);

    WystupMalzFxSmileInversionReport {
        forward_fx_rate: fwd,
        atm_delta_neutral_strike_k_dns: k_dns,
        implied_25d_call_vol: sig_25c,
        implied_25d_put_vol: sig_25p,
        inverted_25d_call_strike_k25c: k_25c,
        inverted_25d_put_strike_k25p: k_25p,
        smile_strangle_minus_broker_fly_bps: fly_conv_bps,
        call_25d_target_delta_error: delta_err,
    }
}

#[cfg(test)]
mod phase42_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn wystup_malz_fx_broker_fly_smile_inversion_invariants_and_ad_parity() {
        let spec = WystupMalzFxSmileInversionSpec {
            spot_fx_s0: 1.085,
            domestic_rate_rd: 0.045,
            foreign_rate_rf: 0.030,
            atm_dns_vol_sigma_atm: 0.085,
            risk_reversal_25d_rr: -0.012,
            butterfly_25d_bf: 0.0035,
            expiry_years_t: 0.5,
        };
        let rep = evaluate_wystup_malz_fx_broker_fly_smile_inversion(&spec);
        assert!(rep.inverted_25d_call_strike_k25c > rep.atm_delta_neutral_strike_k_dns);
        assert!(rep.inverted_25d_put_strike_k25p < rep.atm_delta_neutral_strike_k_dns);
        assert!(rep.call_25d_target_delta_error < 1e-10);

        let h = 1e-6;
        let r_ad =
            evaluate_wystup_malz_fx_broker_fly_smile_inversion(&WystupMalzFxSmileInversionSpec {
                spot_fx_s0: Dual::<1>::from_value(spec.spot_fx_s0),
                domestic_rate_rd: Dual::<1>::from_value(spec.domestic_rate_rd),
                foreign_rate_rf: Dual::<1>::from_value(spec.foreign_rate_rf),
                atm_dns_vol_sigma_atm: Dual::<1>::variable(spec.atm_dns_vol_sigma_atm, 0),
                risk_reversal_25d_rr: Dual::<1>::from_value(spec.risk_reversal_25d_rr),
                butterfly_25d_bf: Dual::<1>::from_value(spec.butterfly_25d_bf),
                expiry_years_t: Dual::<1>::from_value(spec.expiry_years_t),
            });
        let mut sp = spec;
        sp.atm_dns_vol_sigma_atm += h;
        let mut sm = spec;
        sm.atm_dns_vol_sigma_atm -= h;
        let fd = (evaluate_wystup_malz_fx_broker_fly_smile_inversion(&sp)
            .inverted_25d_call_strike_k25c
            - evaluate_wystup_malz_fx_broker_fly_smile_inversion(&sm)
                .inverted_25d_call_strike_k25c)
            / (2.0 * h);
        assert!((r_ad.inverted_25d_call_strike_k25c.d[0] - fd).abs() < 1e-6);
    }
}
