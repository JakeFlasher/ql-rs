//! Bouzoubaa & Osseiran (2010) & Haug (2007) Rubinstein (1991) Chooser Option
//! Put–Call Parity Decomposition & Decision-Flexibility Premium Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct RubinsteinChooserOptionSpec<T: Real = f64> {
    pub spot_price: T,
    pub strike_price: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub volatility: T,
    pub choice_date_tc: T,
    pub final_maturity_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct RubinsteinChooserOptionReport<T: Real = f64> {
    pub underlying_call_t: T,
    pub underlying_put_t: T,
    pub choice_date_effective_put_tc: T,
    pub rubinstein_chooser_price: T,
    pub straddle_upper_bound: T,
    pub decision_flexibility_premium: T,
}

pub fn evaluate_rubinstein_simple_chooser_option<T: Real>(
    spec: &RubinsteinChooserOptionSpec<T>,
) -> RubinsteinChooserOptionReport<T> {
    let s = spec.spot_price.max_val(T::from_f64(1e-4));
    let k = spec.strike_price.max_val(T::from_f64(1e-4));
    let r = spec.risk_free_rate;
    let q = spec.dividend_yield;
    let v = spec.volatility.max_val(T::from_f64(1e-4));
    let tc = spec.choice_date_tc.max_val(T::from_f64(0.01));
    let t = spec.final_maturity_t.max_val(tc + T::from_f64(0.01));

    // Standard Black-Scholes call and put for maturity T:
    let std_t = v * t.sqrt();
    let d1 = ((s / k).ln() + (r - q + T::from_f64(0.5) * v * v) * t) / std_t;
    let d2 = d1 - std_t;
    let df_q_t = (-q * t).exp();
    let df_r_t = (-r * t).exp();
    let call_t = s * df_q_t * norm_cdf(d1) - k * df_r_t * norm_cdf(d2);
    let put_t = k * df_r_t * norm_cdf(-d2) - s * df_q_t * norm_cdf(-d1);

    // By Rubinstein (1991) put-call parity at choice date T_c:
    // max(C(S_Tc, K, T-Tc), P(S_Tc, K, T-Tc)) = C(S_Tc, K, T-Tc) + e^{-q(T-Tc)} max(0, K e^{-(r-q)(T-Tc)} - S_Tc)
    let std_tc = v * tc.sqrt();
    let y1 = ((s / k).ln() + (r - q) * t + T::from_f64(0.5) * v * v * tc) / std_tc;
    let y2 = y1 - std_tc;
    let put_tc = k * df_r_t * norm_cdf(-y2) - s * df_q_t * norm_cdf(-y1);

    let chooser = call_t + put_tc;
    let straddle = call_t + put_t;
    let max_cp = if call_t > put_t { call_t } else { put_t };
    let flex_prem = chooser - max_cp;

    RubinsteinChooserOptionReport {
        underlying_call_t: call_t,
        underlying_put_t: put_t,
        choice_date_effective_put_tc: put_tc,
        rubinstein_chooser_price: chooser,
        straddle_upper_bound: straddle,
        decision_flexibility_premium: flex_prem,
    }
}

#[cfg(test)]
mod phase39_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn rubinstein_chooser_option_invariants_and_ad_parity() {
        let spec = RubinsteinChooserOptionSpec {
            spot_price: 100.0,
            strike_price: 100.0,
            risk_free_rate: 0.05,
            dividend_yield: 0.02,
            volatility: 0.24,
            choice_date_tc: 0.4,
            final_maturity_t: 1.0,
        };
        let rep = evaluate_rubinstein_simple_chooser_option(&spec);
        // No-arbitrage Sandwich: max(C, P) < Chooser < C + P (Straddle)
        assert!(rep.decision_flexibility_premium > 0.0);
        assert!(rep.rubinstein_chooser_price < rep.straddle_upper_bound);

        let h = 1e-6;
        let r_ad = evaluate_rubinstein_simple_chooser_option(&RubinsteinChooserOptionSpec {
            spot_price: Dual::<1>::variable(spec.spot_price, 0),
            strike_price: Dual::<1>::from_value(spec.strike_price),
            risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
            dividend_yield: Dual::<1>::from_value(spec.dividend_yield),
            volatility: Dual::<1>::from_value(spec.volatility),
            choice_date_tc: Dual::<1>::from_value(spec.choice_date_tc),
            final_maturity_t: Dual::<1>::from_value(spec.final_maturity_t),
        });
        let mut sp = spec;
        sp.spot_price += h;
        let mut sm = spec;
        sm.spot_price -= h;
        let fd = (evaluate_rubinstein_simple_chooser_option(&sp).rubinstein_chooser_price
            - evaluate_rubinstein_simple_chooser_option(&sm).rubinstein_chooser_price)
            / (2.0 * h);
        assert!((r_ad.rubinstein_chooser_price.d[0] - fd).abs() < 1e-6);
    }
}
