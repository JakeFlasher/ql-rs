//! Haug (2007) & Bouzoubaa–Osseiran (2010) Geske (1979) Compound Option
//! (`Call-on-Call` & `Put-on-Call`) Critical-Spot Inversion & Parity Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct GeskeCompoundOptionSpec<T: Real = f64> {
    pub spot_price_s0: T,
    pub compound_strike_k1: T,
    pub underlying_strike_k2: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub volatility: T,
    pub compound_expiry_t1: T,
    pub underlying_expiry_t2: T,
}

#[derive(Debug, Clone, Copy)]
pub struct GeskeCompoundOptionReport<T: Real = f64> {
    pub critical_spot_s_star: T,
    pub critical_spot_inversion_error: T,
    pub temporal_correlation_rho12: T,
    pub underlying_bs_call_price: T,
    pub call_on_call_price: T,
    pub put_on_call_price: T,
    pub geske_compound_parity_error: T,
}

fn bs_call<T: Real>(s: T, k: T, r: T, q: T, v: T, tau: T) -> (T, T) {
    let std = v * tau.sqrt();
    let d1 = ((s / k).ln() + (r - q + T::from_f64(0.5) * v * v) * tau) / std;
    let d2 = d1 - std;
    let df_q = (-q * tau).exp();
    let df_r = (-r * tau).exp();
    let price = s * df_q * norm_cdf(d1) - k * df_r * norm_cdf(d2);
    let delta = df_q * norm_cdf(d1);
    (price, delta)
}

fn bvn_cdf<T: Real>(a: T, b: T, rho: T) -> T {
    const X: [f64; 5] = [
        0.1488743389816312,
        0.4333953941292472,
        0.6794095682990244,
        0.8650633666889845,
        0.9739065285171717,
    ];
    const W: [f64; 5] = [
        0.2955242247147529,
        0.2692667193099963,
        0.219086362515982,
        0.1494513491505806,
        0.0666713443086881,
    ];
    let pa = norm_cdf(a);
    let pb = norm_cdf(b);
    let angle = (rho / (T::one() - rho * rho).max_val(T::from_f64(1e-8)).sqrt()).atan();
    let mut sum = T::zero();
    for i in 0..5 {
        for &sign in &[-1.0_f64, 1.0_f64] {
            let x = T::from_f64(sign * X[i]);
            let theta = angle * (T::one() - x) * T::from_f64(0.5);
            let s = theta.sin();
            let denom = (T::one() - s * s).max_val(T::from_f64(1e-8));
            let exponent = (s * a * b - T::from_f64(0.5) * (a * a + b * b)) / denom;
            sum += exponent.exp() * T::from_f64(W[i]);
        }
    }
    pa * pb + angle * sum / T::from_f64(4.0 * std::f64::consts::PI)
}

pub fn evaluate_geske_compound_call_and_put_on_call<T: Real>(
    spec: &GeskeCompoundOptionSpec<T>,
) -> GeskeCompoundOptionReport<T> {
    let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
    let k1 = spec.compound_strike_k1.max_val(T::from_f64(1e-4));
    let k2 = spec.underlying_strike_k2.max_val(T::from_f64(1e-4));
    let r = spec.risk_free_rate;
    let q = spec.dividend_yield;
    let v = spec.volatility.max_val(T::from_f64(1e-4));
    let t1 = spec.compound_expiry_t1.max_val(T::from_f64(0.05));
    let t2 = spec.underlying_expiry_t2.max_val(t1 + T::from_f64(0.05));
    let tau21 = t2 - t1;

    // 6-step Newton inversion for critical underlying price S^* at T_1 where C_BS(S^*, K_2, T_2-T_1) = K_1:
    let mut s_star = k2 + k1;
    for _ in 0..6 {
        let (c_val, c_del) = bs_call(s_star, k2, r, q, v, tau21);
        s_star =
            (s_star - (c_val - k1) / c_del.max_val(T::from_f64(1e-6))).max_val(T::from_f64(1e-4));
    }
    let (c_check, _) = bs_call(s_star, k2, r, q, v, tau21);
    let inv_err = (c_check - k1).abs();

    let rho = (t1 / t2).sqrt();
    let (c_under, _) = bs_call(s0, k2, r, q, v, t2);

    // Geske (1979) exact closed-form bivariate normal pricing for Call-on-Call & Put-on-Call:
    let std1 = v * t1.sqrt();
    let std2 = v * t2.sqrt();
    let a1 = ((s0 / s_star).ln() + (r - q + T::from_f64(0.5) * v * v) * t1) / std1;
    let a2 = a1 - std1;
    let b1 = ((s0 / k2).ln() + (r - q + T::from_f64(0.5) * v * v) * t2) / std2;
    let b2 = b1 - std2;
    let df_q2 = (-q * t2).exp();
    let df_r2 = (-r * t2).exp();
    let df_r1 = (-r * t1).exp();

    let cc = s0 * df_q2 * bvn_cdf(a1, b1, rho)
        - k2 * df_r2 * bvn_cdf(a2, b2, rho)
        - k1 * df_r1 * norm_cdf(a2);
    let pc = k2 * df_r2 * bvn_cdf(-a2, b2, -rho) - s0 * df_q2 * bvn_cdf(-a1, b1, -rho)
        + k1 * df_r1 * norm_cdf(-a2);
    let parity_err = ((cc - pc) - (c_under - k1 * df_r1)).abs();

    GeskeCompoundOptionReport {
        critical_spot_s_star: s_star,
        critical_spot_inversion_error: inv_err,
        temporal_correlation_rho12: rho,
        underlying_bs_call_price: c_under,
        call_on_call_price: cc,
        put_on_call_price: pc,
        geske_compound_parity_error: parity_err,
    }
}

#[cfg(test)]
mod phase40_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn geske_compound_option_invariants_and_ad_parity() {
        let spec = GeskeCompoundOptionSpec {
            spot_price_s0: 100.0,
            compound_strike_k1: 6.5,
            underlying_strike_k2: 102.0,
            risk_free_rate: 0.04,
            dividend_yield: 0.015,
            volatility: 0.24,
            compound_expiry_t1: 0.5,
            underlying_expiry_t2: 1.5,
        };
        let rep = evaluate_geske_compound_call_and_put_on_call(&spec);
        assert!(rep.critical_spot_inversion_error < 1e-10);
        assert!(rep.geske_compound_parity_error < 1e-12);
        assert!(rep.call_on_call_price < rep.underlying_bs_call_price);

        let h = 1e-6;
        let r_ad = evaluate_geske_compound_call_and_put_on_call(&GeskeCompoundOptionSpec {
            spot_price_s0: Dual::<1>::variable(spec.spot_price_s0, 0),
            compound_strike_k1: Dual::<1>::from_value(spec.compound_strike_k1),
            underlying_strike_k2: Dual::<1>::from_value(spec.underlying_strike_k2),
            risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
            dividend_yield: Dual::<1>::from_value(spec.dividend_yield),
            volatility: Dual::<1>::from_value(spec.volatility),
            compound_expiry_t1: Dual::<1>::from_value(spec.compound_expiry_t1),
            underlying_expiry_t2: Dual::<1>::from_value(spec.underlying_expiry_t2),
        });
        let mut sp = spec;
        sp.spot_price_s0 += h;
        let mut sm = spec;
        sm.spot_price_s0 -= h;
        let fd = (evaluate_geske_compound_call_and_put_on_call(&sp).call_on_call_price
            - evaluate_geske_compound_call_and_put_on_call(&sm).call_on_call_price)
            / (2.0 * h);
        assert!((r_ad.call_on_call_price.d[0] - fd).abs() < 1e-5);
    }
}
