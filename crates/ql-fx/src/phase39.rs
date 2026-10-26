//! Clark (2011) & Wystup (2017) FX Forward-Volatility Agreement (`FVA`),
//! Forward-Starting Smile Roll-Down & Calendar-Spread No-Arbitrage Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct ClarkFxForwardVolAgreementSpec<T: Real = f64> {
    pub spot_fx: T,
    pub domestic_rate_rd: T,
    pub foreign_rate_rf: T,
    pub front_expiry_t1: T,
    pub back_expiry_t2: T,
    pub atm_vol_t1: T,
    pub atm_vol_t2: T,
    pub smile_convexity_adj_bps: T,
}

#[derive(Debug, Clone, Copy)]
pub struct ClarkFxForwardVolAgreementReport<T: Real = f64> {
    pub total_variance_t1: T,
    pub total_variance_t2: T,
    pub calendar_variance_margin: T,
    pub forward_atm_implied_vol: T,
    pub fva_fair_vol_strike: T,
    pub forward_start_atm_straddle_pct: T,
}

pub fn evaluate_clark_fx_forward_volatility_agreement<T: Real>(
    spec: &ClarkFxForwardVolAgreementSpec<T>,
) -> ClarkFxForwardVolAgreementReport<T> {
    let t1 = spec.front_expiry_t1.max_val(T::from_f64(0.05));
    let t2 = spec.back_expiry_t2.max_val(t1 + T::from_f64(0.05));
    let dt = t2 - t1;

    let w1 = spec.atm_vol_t1 * spec.atm_vol_t1 * t1;
    let w2 = spec.atm_vol_t2 * spec.atm_vol_t2 * t2;
    let cal_margin = w2 - w1;

    // Forward implied variance sigma_fwd^2(T1, T2) = (sigma(T2)^2 T2 - sigma(T1)^2 T1) / (T2 - T1)
    let fwd_var = (cal_margin / dt).max_val(T::from_f64(1e-8));
    let fwd_vol = fwd_var.sqrt();

    // Smile-convexity-adjusted FVA strike (Clark Ch. 8):
    let fva_strike = fwd_vol + spec.smile_convexity_adj_bps * T::from_f64(1e-4);

    // Forward-starting ATM-forward straddle setting strike at T1 equal to F(T1, T2):
    // Discounted to t=0 in % of initial spot S0:
    let std = fva_strike * dt.sqrt();
    let d1 = T::from_f64(0.5) * std;
    let d2 = -d1;
    let df_d_t1 = (-spec.domestic_rate_rd * t1).exp();
    let df_f_dt = (-spec.foreign_rate_rf * dt).exp();
    let straddle_pct =
        T::from_f64(2.0) * df_d_t1 * df_f_dt * (norm_cdf(d1) - norm_cdf(d2)) * spec.spot_fx;

    ClarkFxForwardVolAgreementReport {
        total_variance_t1: w1,
        total_variance_t2: w2,
        calendar_variance_margin: cal_margin,
        forward_atm_implied_vol: fwd_vol,
        fva_fair_vol_strike: fva_strike,
        forward_start_atm_straddle_pct: straddle_pct,
    }
}

#[cfg(test)]
mod phase39_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn clark_fx_forward_vol_agreement_invariants_and_ad_parity() {
        let spec = ClarkFxForwardVolAgreementSpec {
            spot_fx: 1.08,
            domestic_rate_rd: 0.045,
            foreign_rate_rf: 0.030,
            front_expiry_t1: 0.5,
            back_expiry_t2: 1.5,
            atm_vol_t1: 0.085,
            atm_vol_t2: 0.098,
            smile_convexity_adj_bps: -18.0,
        };
        let rep = evaluate_clark_fx_forward_volatility_agreement(&spec);
        assert!(rep.calendar_variance_margin > 0.0);
        assert!(rep.forward_atm_implied_vol > spec.atm_vol_t2);
        assert!(rep.forward_start_atm_straddle_pct > 0.0);

        let h = 1e-6;
        let r_ad =
            evaluate_clark_fx_forward_volatility_agreement(&ClarkFxForwardVolAgreementSpec {
                spot_fx: Dual::<1>::from_value(spec.spot_fx),
                domestic_rate_rd: Dual::<1>::from_value(spec.domestic_rate_rd),
                foreign_rate_rf: Dual::<1>::from_value(spec.foreign_rate_rf),
                front_expiry_t1: Dual::<1>::from_value(spec.front_expiry_t1),
                back_expiry_t2: Dual::<1>::from_value(spec.back_expiry_t2),
                atm_vol_t1: Dual::<1>::from_value(spec.atm_vol_t1),
                atm_vol_t2: Dual::<1>::variable(spec.atm_vol_t2, 0),
                smile_convexity_adj_bps: Dual::<1>::from_value(spec.smile_convexity_adj_bps),
            });
        let mut sp = spec;
        sp.atm_vol_t2 += h;
        let mut sm = spec;
        sm.atm_vol_t2 -= h;
        let fd = (evaluate_clark_fx_forward_volatility_agreement(&sp)
            .forward_start_atm_straddle_pct
            - evaluate_clark_fx_forward_volatility_agreement(&sm).forward_start_atm_straddle_pct)
            / (2.0 * h);
        assert!((r_ad.forward_start_atm_straddle_pct.d[0] - fd).abs() < 1e-6);
    }
}
