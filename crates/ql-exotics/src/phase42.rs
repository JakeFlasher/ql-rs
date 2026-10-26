//! Haug (2007) & Goldman–Sosin–Gatto (1979) Floating-Strike & Fixed-Strike
//! Lookback Option Reflection-Principle Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct GoldmanSosinGattoLookbackSpec<T: Real = f64> {
    pub spot_price_s0: T,
    pub realized_minimum_smin: T,
    pub realized_maximum_smax: T,
    pub fixed_strike_k: T,
    pub risk_free_rate_r: T,
    pub dividend_yield_q: T,
    pub volatility_sigma: T,
    pub maturity_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct GoldmanSosinGattoLookbackReport<T: Real = f64> {
    pub floating_strike_lookback_call: T,
    pub floating_strike_lookback_put: T,
    pub fixed_strike_lookback_call: T,
    pub reflection_principle_bonus_value: T,
    pub lookback_over_vanilla_ratio: T,
    pub minimum_reset_sensitivity: T,
}

pub fn evaluate_goldman_sosin_gatto_lookback_reflection<T: Real>(
    spec: &GoldmanSosinGattoLookbackSpec<T>,
) -> GoldmanSosinGattoLookbackReport<T> {
    let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
    let s_min = spec
        .realized_minimum_smin
        .max_val(T::from_f64(1e-4))
        .min_val(s0);
    let s_max = spec.realized_maximum_smax.max_val(s0);
    let k = spec.fixed_strike_k.max_val(T::from_f64(1e-4));
    let r = spec.risk_free_rate_r;
    let q = spec.dividend_yield_q;
    let b = (r - q).max_val(T::from_f64(0.005));
    let v = spec.volatility_sigma.max_val(T::from_f64(1e-4));
    let t = spec.maturity_years_t.max_val(T::from_f64(0.05));
    let std = v * t.sqrt();

    let df_r = (-r * t).exp();
    let df_q = (-q * t).exp();
    let lam = T::from_f64(2.0) * b / (v * v);

    // Goldman-Sosin-Gatto (1979) / Haug (2007, Ch. 2.1) Floating-Strike Lookback Call:
    let a1 = ((s0 / s_min).ln() + (b + T::from_f64(0.5) * v * v) * t) / std;
    let a2 = a1 - std;
    let a3 = ((s0 / s_min).ln() + (-b + T::from_f64(0.5) * v * v) * t) / std;

    let vanilla_smin = s0 * df_q * norm_cdf(a1) - s_min * df_r * norm_cdf(a2);
    let refl_call =
        (s0 / lam) * ((s0 / s_min).pow(-lam) * df_r * norm_cdf(-a3) - df_q * norm_cdf(-a1));
    let c_float = vanilla_smin + refl_call;

    // Floating-Strike Lookback Put (struck at realized maximum S_max):
    let b1 = ((s0 / s_max).ln() + (b + T::from_f64(0.5) * v * v) * t) / std;
    let b2 = b1 - std;
    let b3 = ((s0 / s_max).ln() + (-b + T::from_f64(0.5) * v * v) * t) / std;
    let p_float = s_max * df_r * norm_cdf(-b2) - s0 * df_q * norm_cdf(-b1)
        + (s0 / lam) * (df_q * norm_cdf(b1) - (s0 / s_max).pow(-lam) * df_r * norm_cdf(b3));

    // Conze-Viswanathan (1991) Fixed-Strike Lookback Call:
    let m_eff = s_max.max_val(k);
    let e1 = ((s0 / m_eff).ln() + (b + T::from_f64(0.5) * v * v) * t) / std;
    let e2 = e1 - std;
    let e3 = ((s0 / m_eff).ln() + (-b + T::from_f64(0.5) * v * v) * t) / std;
    let c_fixed = (m_eff - k) * df_r + s0 * df_q * norm_cdf(e1) - m_eff * df_r * norm_cdf(e2)
        + (s0 / lam) * (df_q * norm_cdf(e1) - (s0 / m_eff).pow(-lam) * df_r * norm_cdf(e3));

    let ratio = c_float / vanilla_smin.max_val(T::from_f64(1e-4));
    let min_sens = -df_r * norm_cdf(a2);

    GoldmanSosinGattoLookbackReport {
        floating_strike_lookback_call: c_float,
        floating_strike_lookback_put: p_float,
        fixed_strike_lookback_call: c_fixed,
        reflection_principle_bonus_value: refl_call,
        lookback_over_vanilla_ratio: ratio,
        minimum_reset_sensitivity: min_sens,
    }
}

#[cfg(test)]
mod phase42_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn goldman_sosin_gatto_lookback_reflection_invariants_and_ad_parity() {
        let spec = GoldmanSosinGattoLookbackSpec {
            spot_price_s0: 100.0,
            realized_minimum_smin: 92.0,
            realized_maximum_smax: 106.0,
            fixed_strike_k: 100.0,
            risk_free_rate_r: 0.05,
            dividend_yield_q: 0.015,
            volatility_sigma: 0.25,
            maturity_years_t: 1.0,
        };
        let rep = evaluate_goldman_sosin_gatto_lookback_reflection(&spec);
        assert!(rep.reflection_principle_bonus_value > 0.0);
        assert!(rep.lookback_over_vanilla_ratio > 1.0);
        assert!(rep.fixed_strike_lookback_call > 0.0);

        let h = 1e-6;
        let r_ad =
            evaluate_goldman_sosin_gatto_lookback_reflection(&GoldmanSosinGattoLookbackSpec {
                spot_price_s0: Dual::<1>::variable(spec.spot_price_s0, 0),
                realized_minimum_smin: Dual::<1>::from_value(spec.realized_minimum_smin),
                realized_maximum_smax: Dual::<1>::from_value(spec.realized_maximum_smax),
                fixed_strike_k: Dual::<1>::from_value(spec.fixed_strike_k),
                risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
                dividend_yield_q: Dual::<1>::from_value(spec.dividend_yield_q),
                volatility_sigma: Dual::<1>::from_value(spec.volatility_sigma),
                maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
            });
        let mut sp = spec;
        sp.spot_price_s0 += h;
        let mut sm = spec;
        sm.spot_price_s0 -= h;
        let fd = (evaluate_goldman_sosin_gatto_lookback_reflection(&sp)
            .floating_strike_lookback_call
            - evaluate_goldman_sosin_gatto_lookback_reflection(&sm).floating_strike_lookback_call)
            / (2.0 * h);
        assert!((r_ad.floating_strike_lookback_call.d[0] - fd).abs() < 1e-5);
    }
}
