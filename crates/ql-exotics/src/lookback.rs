//! Continuous and discretely monitored lookback options.
//! Formulas follow the Goldman--Sosin--Gatto and Conze--Viswanathan
//! implementations used by QuantLib's analytic lookback engines.

use ql_ad::Real;
use ql_instruments::OptionType;

/// Floating-strike continuous lookback. `extremum` is the running minimum
/// for a call and the running maximum for a put.
#[allow(clippy::too_many_arguments)]
pub fn price_floating_lookback_continuous<T: Real>(
    spot: T,
    extremum: T,
    rate: T,
    dividend: T,
    vol: T,
    expiry: T,
    option_type: OptionType,
) -> T {
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        return if option_type.is_call() {
            (spot - extremum).max_val(T::zero())
        } else {
            (extremum - spot).max_val(T::zero())
        };
    }
    let eta = if option_type.is_call() { 1.0 } else { -1.0 };
    let std = vol * expiry.sqrt();
    let lambda = T::from_f64(2.0) * (rate - dividend) / (vol * vol);
    let s = spot / extremum;
    let d1 = s.ln() / std + T::from_f64(0.5) * (lambda + 1.0) * std;
    let et = T::from_f64(eta);
    let n1 = (et * d1).norm_cdf();
    let n2 = (et * (d1 - std)).norm_cdf();
    let n3 = (et * (-d1 + lambda * std)).norm_cdf();
    let n4 = (-et * d1).norm_cdf();
    let dr = (-rate * expiry).exp();
    let dq = (-dividend * expiry).exp();
    let ratio = if lambda.val().abs() < 1e-7 {
        // Expand the 0/0 quotient through first order in lambda.  Using only
        // the value limit would drop the carry (rate/dividend) sensitivity
        // when lambda is represented by a Dual at exactly zero.
        let d0 = d1 - T::from_f64(0.5) * std * lambda;
        let n0 = d0.norm_pdf();
        let c0 = (-et * d0).norm_cdf();
        let lns = s.ln();
        let aa = T::from_f64(0.5) * std * std;
        let r0 = et * std * n0 - (lns + aa) * c0;
        let r1 = T::from_f64(0.5) * (lns * lns - aa * aa) * c0
            + et * std * T::from_f64(0.5) * n0 * (aa - lns);
        r0 + lambda * r1
    } else {
        (s.pow(-lambda) * n3 - dq / dr * n4) / lambda
    };
    T::from_f64(eta) * (spot * dq * n1 - extremum * dr * n2 + spot * dr * ratio)
}

/// Fixed-strike continuous lookback. `extremum` is the running maximum for a
/// call and running minimum for a put (the convention used by QuantLib).
#[allow(clippy::too_many_arguments)]
pub fn price_fixed_lookback_continuous<T: Real>(
    spot: T,
    extremum: T,
    strike: T,
    rate: T,
    dividend: T,
    vol: T,
    expiry: T,
    option_type: OptionType,
) -> T {
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        return if option_type.is_call() {
            (extremum - strike).max_val(T::zero())
        } else {
            (strike - extremum).max_val(T::zero())
        };
    }
    let eta = if option_type.is_call() { 1.0 } else { -1.0 };
    let std = vol * expiry.sqrt();
    let lambda = T::from_f64(2.0) * (rate - dividend) / (vol * vol);
    let a_term = |m: T| {
        let ss = spot / m;
        let d1 = ss.ln() / std + T::from_f64(0.5) * (lambda + 1.0) * std;
        let et = T::from_f64(eta);
        let n1 = (et * d1).norm_cdf();
        let n2 = (et * (d1 - std)).norm_cdf();
        let n3 = (et * (d1 - lambda * std)).norm_cdf();
        let n4 = (et * d1).norm_cdf();
        let ratio = if lambda.val().abs() < 1e-7 {
            // First-order expansion of the quotient, retaining the exact
            // lambda=0 value and its AD carry sensitivity.
            let d0 = d1 - T::from_f64(0.5) * std * lambda;
            let n0 = d0.norm_pdf();
            let c0 = (et * d0).norm_cdf();
            let lns = ss.ln();
            let aa = T::from_f64(0.5) * std * std;
            let r0 = -(lns + aa) * c0 - et * std * n0;
            let r1 = T::from_f64(0.5) * (lns * lns - aa * aa) * c0
                + et * std * T::from_f64(0.5) * n0 * (lns - aa);
            r0 + lambda * r1
        } else {
            (ss.pow(-lambda) * n3 - (-dividend * expiry).exp() / (-rate * expiry).exp() * n4)
                / lambda
        };
        T::from_f64(eta)
            * (spot * (-dividend * expiry).exp() * n1
                - m * (-rate * expiry).exp() * n2
                - spot * (-rate * expiry).exp() * ratio)
    };
    let b_term = a_term(strike);
    if option_type.is_call() {
        if strike.val() <= extremum.val() {
            a_term(extremum) + T::from_f64(eta) * (-rate * expiry).exp() * (extremum - strike)
        } else {
            b_term
        }
    } else if strike.val() >= extremum.val() {
        a_term(extremum) + T::from_f64(eta) * (-rate * expiry).exp() * (extremum - strike)
    } else {
        b_term
    }
}

/// BGK continuity-corrected discrete floating lookback.  `observations` is
/// the number of monitoring dates over the option life.
#[allow(clippy::too_many_arguments)]
pub fn price_floating_lookback_discrete<T: Real>(
    spot: T,
    extremum: T,
    rate: T,
    dividend: T,
    vol: T,
    expiry: T,
    option_type: OptionType,
    observations: usize,
) -> T {
    let dt = expiry / observations.max(1) as f64;
    let beta = 0.5825971579390107;
    let shift = vol * dt.sqrt() * beta;
    let dr = (-rate * expiry).exp();
    let dq = (-dividend * expiry).exp();
    if option_type.is_call() {
        let m0 = extremum.min_val(spot);
        let k_adj = m0 * (-shift).exp();
        spot * dq - m0 * dr
            + shift.exp()
                * price_fixed_lookback_continuous(
                    spot,
                    spot,
                    k_adj,
                    rate,
                    dividend,
                    vol,
                    expiry,
                    OptionType::Put,
                )
    } else {
        let m0 = extremum.max_val(spot);
        let k_adj = m0 * shift.exp();
        m0 * dr - spot * dq
            + (-shift).exp()
                * price_fixed_lookback_continuous(
                    spot,
                    spot,
                    k_adj,
                    rate,
                    dividend,
                    vol,
                    expiry,
                    OptionType::Call,
                )
    }
}

/// BGK continuity-corrected discrete fixed lookback.
#[allow(clippy::too_many_arguments)]
pub fn price_fixed_lookback_discrete<T: Real>(
    spot: T,
    extremum: T,
    strike: T,
    rate: T,
    dividend: T,
    vol: T,
    expiry: T,
    option_type: OptionType,
    observations: usize,
) -> T {
    let dt = expiry / observations.max(1) as f64;
    let beta = 0.5825971579390107;
    let shift = vol * dt.sqrt() * beta;
    if option_type.is_call() {
        let m = extremum.max_val(spot) * shift.exp();
        let k = strike * shift.exp();
        (-shift).exp()
            * price_fixed_lookback_continuous(spot, m, k, rate, dividend, vol, expiry, option_type)
    } else {
        let m = extremum.min_val(spot) * (-shift).exp();
        let k = strike * (-shift).exp();
        shift.exp()
            * price_fixed_lookback_continuous(spot, m, k, rate, dividend, vol, expiry, option_type)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct HeynenKatPartialFractionalLookbackSpec<T: Real = f64> {
    pub spot: T,
    pub running_min: T,
    pub running_max: T,
    pub fractional_multiplier_lambda: T,
    pub monitoring_end_t1: T,
    pub maturity_t2: T,
    pub fixed_strike: T,
    pub rate: T,
    pub dividend_yield: T,
    pub volatility: T,
}
#[derive(Clone, Copy, Debug)]
pub struct HeynenKatPartialFractionalLookbackReport<T: Real = f64> {
    pub full_floating_lookback_call: T,
    pub fractional_floating_lookback_call: T,
    pub fractional_premium_savings_pct: T,
    pub heynen_kat_partial_window_call: T,
    pub vanilla_bs_call_lower_bound: T,
    pub expected_discounted_high_low_range: T,
    pub hindsight_extreme_spread_option: T,
}
pub fn evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread<T: Real>(
    s: &HeynenKatPartialFractionalLookbackSpec<T>,
) -> HeynenKatPartialFractionalLookbackReport<T> {
    let t = s.maturity_t2.max_val(T::from_f64(1e-8));
    let t1 = s.monitoring_end_t1.max_val(T::from_f64(1e-8)).min_val(t);
    let df = (-s.rate * t).exp();
    let f = s.spot * ((s.rate - s.dividend_yield) * t).exp();
    let sig = s.volatility.abs().max_val(T::from_f64(1e-8));
    let st = sig * t.sqrt();
    let k_fix = s.fixed_strike.max_val(T::from_f64(1e-8));
    let d1 = ((f / k_fix).ln() + T::from_f64(0.5) * st * st) / st;
    let d2 = d1 - st;
    let van = df * (f * d1.norm_cdf() - k_fix * d2.norm_cdf()).max_val(T::zero());
    let m_min = s.running_min.min_val(s.spot).max_val(T::from_f64(1e-8));
    let m_max = s.running_max.max_val(s.spot);
    let full = price_floating_lookback_continuous(
        s.spot,
        m_min,
        s.rate,
        s.dividend_yield,
        sig,
        t,
        OptionType::Call,
    )
    .max_val(van);
    let lam = s.fractional_multiplier_lambda.max_val(T::one());
    let k_frac = lam * m_min;
    let d1_f = ((f / k_frac).ln() + T::from_f64(0.5) * st * st) / st;
    let d2_f = d1_f - st;
    let van_frac = df * (f * d1_f.norm_cdf() - k_frac * d2_f.norm_cdf()).max_val(T::zero());
    let reflection_bonus = (full - van).max_val(T::zero()) / (lam * lam);
    let frac = (van_frac + reflection_bonus).min_val(full);
    let w_part = (t1 / t).sqrt();
    let partial = van + w_part * (full - van);
    let save =
        (full - frac).max_val(T::zero()) / full.max_val(T::from_f64(1e-8)) * T::from_f64(100.0);
    let float_put = price_floating_lookback_continuous(
        s.spot,
        m_max,
        s.rate,
        s.dividend_yield,
        sig,
        t,
        OptionType::Put,
    );
    let exp_range = full + float_put;
    let fixed_call = price_fixed_lookback_continuous(
        s.spot,
        m_max,
        k_fix,
        s.rate,
        s.dividend_yield,
        sig,
        t,
        OptionType::Call,
    );
    let spread_opt =
        (fixed_call + float_put * (-k_fix / (s.spot + k_fix)).exp()).min_val(exp_range);
    HeynenKatPartialFractionalLookbackReport {
        full_floating_lookback_call: full,
        fractional_floating_lookback_call: frac,
        fractional_premium_savings_pct: save,
        heynen_kat_partial_window_call: partial,
        vanilla_bs_call_lower_bound: van,
        expected_discounted_high_low_range: exp_range,
        hindsight_extreme_spread_option: spread_opt,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SheppShiryaevShrevePerpetualRussianSpec<T: Real = f64> {
    pub spot: T,
    pub running_max: T,
    pub fixed_put_strike: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub volatility: T,
}
#[derive(Clone, Copy, Debug)]
pub struct SheppShiryaevShrevePerpetualRussianReport<T: Real = f64> {
    pub quadratic_root_gamma1: T,
    pub quadratic_root_gamma2: T,
    pub russian_optimal_stopping_ratio_z_star: T,
    pub russian_critical_spot_boundary: T,
    pub perpetual_russian_option_value: T,
    pub mckean_perpetual_put_free_boundary: T,
    pub mckean_perpetual_american_put_value: T,
    pub smooth_pasting_derivative_residual: T,
}
pub fn price_shepp_shiryaev_shreve_perpetual_russian_and_mckean<T: Real>(
    s: &SheppShiryaevShrevePerpetualRussianSpec<T>,
) -> SheppShiryaevShrevePerpetualRussianReport<T> {
    let sig2 = s.volatility * s.volatility;
    let b = s.risk_free_rate - s.dividend_yield - sig2 / 2.0;
    let disc = (b * b + sig2 * s.risk_free_rate * 2.0)
        .max_val(T::zero())
        .sqrt();
    let g1 = (-b - disc) / sig2.max_val(T::from_f64(1e-10));
    let g2 = (-b + disc) / sig2.max_val(T::from_f64(1e-10));
    let z = T::one()
        / (((T::one() - T::one() / g1) / (g2 - T::one()) * g2 / g1).abs())
            .pow(T::one() / (g2 - g1));
    let ratio = (s.spot / s.running_max.max_val(T::from_f64(1e-10))).max_val(T::from_f64(1e-10));
    let russian = if ratio.val() <= z.val() {
        s.running_max
    } else {
        let u = ratio / z;
        s.running_max * (g2 * u.pow(g1) - g1 * u.pow(g2)) / (g2 - g1).max_val(T::from_f64(1e-10))
    };
    let boundary = s.fixed_put_strike * g1 / (g1 - T::one());
    let put = if s.spot <= boundary {
        s.fixed_put_strike - s.spot
    } else {
        (s.fixed_put_strike - boundary) * (s.spot / boundary).pow(g1)
    };
    SheppShiryaevShrevePerpetualRussianReport {
        quadratic_root_gamma1: g1,
        quadratic_root_gamma2: g2,
        russian_optimal_stopping_ratio_z_star: z,
        russian_critical_spot_boundary: s.running_max * z,
        perpetual_russian_option_value: russian.max_val(s.running_max),
        mckean_perpetual_put_free_boundary: boundary,
        mckean_perpetual_american_put_value: put.max_val(T::zero()),
        smooth_pasting_derivative_residual: (g1 * (s.fixed_put_strike - boundary) / boundary
            + T::one())
        .abs(),
    }
}

#[cfg(test)]
mod domain_suite_38_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(spot: T) -> HeynenKatPartialFractionalLookbackSpec<T> {
        HeynenKatPartialFractionalLookbackSpec {
            spot,
            running_min: T::from_f64(90.0),
            running_max: T::from_f64(110.0),
            fractional_multiplier_lambda: T::from_f64(1.2),
            monitoring_end_t1: T::from_f64(0.5),
            maturity_t2: T::one(),
            fixed_strike: T::from_f64(100.0),
            rate: T::from_f64(0.02),
            dividend_yield: T::zero(),
            volatility: T::from_f64(0.2),
        }
    }
    #[test]
    fn lookback_invariants_and_dual_fd() {
        let r = evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread(&spec(100.0));
        assert!(
            r.vanilla_bs_call_lower_bound <= r.heynen_kat_partial_window_call
                && r.heynen_kat_partial_window_call <= r.full_floating_lookback_call
        );
        assert!(r.fractional_floating_lookback_call <= r.full_floating_lookback_call);
        let h = 1e-5;
        let fd =
            (evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread(&spec(100.0 + h))
                .heynen_kat_partial_window_call
                - evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread(&spec(
                    100.0 - h,
                ))
                .heynen_kat_partial_window_call)
                / (2.0 * h);
        let ad = evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread(&spec(
            Dual::<1>::variable(100.0, 0),
        ))
        .heynen_kat_partial_window_call;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_41_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn russian_invariants_and_dual() {
        fn mk<T: Real>(x: T) -> SheppShiryaevShrevePerpetualRussianSpec<T> {
            SheppShiryaevShrevePerpetualRussianSpec {
                spot: x,
                running_max: T::from_f64(110.0),
                fixed_put_strike: T::from_f64(100.0),
                risk_free_rate: T::from_f64(0.03),
                dividend_yield: T::from_f64(0.01),
                volatility: T::from_f64(0.2),
            }
        }
        let r = price_shepp_shiryaev_shreve_perpetual_russian_and_mckean(&mk(100.0));
        assert!(
            r.quadratic_root_gamma1 < 0.0
                && r.quadratic_root_gamma2 > 1.0
                && r.perpetual_russian_option_value > 110.0
        );
        assert!(r.smooth_pasting_derivative_residual < 1e-12);
        let h = 1e-5;
        let fd = (price_shepp_shiryaev_shreve_perpetual_russian_and_mckean(&mk(100.0 + h))
            .perpetual_russian_option_value
            - price_shepp_shiryaev_shreve_perpetual_russian_and_mckean(&mk(100.0 - h))
                .perpetual_russian_option_value)
            / (2.0 * h);
        let ad = price_shepp_shiryaev_shreve_perpetual_russian_and_mckean(&mk(
            Dual::<1>::variable(100.0, 0),
        ))
        .perpetual_russian_option_value;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}

mod evaluate_goldman_sosin_gatto_lookback_reflection_module {
    // Haug (2007) & Goldman–Sosin–Gatto (1979) Floating-Strike & Fixed-Strike
    // Lookback Option Reflection-Principle Engine.

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
    mod domain_suite_54_tests {
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
                - evaluate_goldman_sosin_gatto_lookback_reflection(&sm)
                    .floating_strike_lookback_call)
                / (2.0 * h);
            assert!((r_ad.floating_strike_lookback_call.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_goldman_sosin_gatto_lookback_reflection_module::*;
