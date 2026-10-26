use ql_ad::{Real, ift_lift_1d};
use ql_instruments::OptionType;
use ql_math::{bivariate_normal_cdf, brent_solve};
use ql_vol::black_scholes_price;

#[derive(Copy, Clone, Debug)]
pub struct CompoundOptionSpec<T: Real> {
    pub mother_type: OptionType,
    pub daughter_type: OptionType,
    pub spot: T,
    pub strike_mother: T,
    pub strike_daughter: T,
    pub rate: T,
    pub div_yield: T,
    pub vol: T,
    pub expiry_mother: T,
    pub expiry_daughter: T,
}

fn daughter_delta<T: Real>(spot: T, spec: &CompoundOptionSpec<T>, tau: T) -> T {
    if tau.val() <= 0.0 {
        let sign = if spec.daughter_type.is_call() {
            1.0
        } else {
            -1.0
        };
        return T::from_f64(sign)
            * (-spec.div_yield * tau).exp()
            * if (spot - spec.strike_daughter).val() * sign > 0.0 {
                T::one()
            } else {
                T::zero()
            };
    }
    let std = spec.vol * tau.sqrt();
    let d1 = ((spot / spec.strike_daughter).ln()
        + (spec.rate - spec.div_yield + spec.vol * spec.vol * 0.5) * tau)
        / std;
    let sign = if spec.daughter_type.is_call() {
        1.0
    } else {
        -1.0
    };
    T::from_f64(sign) * (-spec.div_yield * tau).exp() * (T::from_f64(sign) * d1).norm_cdf()
}

pub fn price_compound_option<T: Real>(spec: &CompoundOptionSpec<T>) -> T {
    let tau12 = spec.expiry_daughter - spec.expiry_mother;
    assert!(tau12.val() >= 0.0);
    let k1 = spec.strike_mother;
    let k2 = spec.strike_daughter;
    let f = |s: f64| {
        black_scholes_price(
            s,
            k2.val(),
            spec.vol.val(),
            tau12.val(),
            spec.rate.val(),
            spec.div_yield.val(),
            spec.daughter_type.is_call(),
        ) - k1.val()
    };
    let root = match brent_solve(f, 1e-10, 1e5_f64.max(k2.val() * 100.0), 1e-12, 200) {
        Ok(x) => x,
        Err(_) => {
            let mut x = k2.val().max(1e-8);
            for _ in 0..64 {
                let fx = f(x);
                let h = (x.abs() * 1e-5).max(1e-6);
                let dx = (f(x + h) - f(x - h)) / (2.0 * h);
                if !fx.is_finite() || !dx.is_finite() || dx.abs() < 1e-14 {
                    break;
                }
                let xn = (x - fx / dx).max(1e-10);
                if (xn - x).abs() < 1e-12 {
                    x = xn;
                    break;
                }
                x = xn;
            }
            x
        }
    };
    let d1_root = ((root / k2.val()).ln()
        + (spec.rate.val() - spec.div_yield.val() + 0.5 * spec.vol.val() * spec.vol.val())
            * tau12.val())
        / (spec.vol.val() * tau12.val().sqrt());
    let fxx = (-spec.div_yield.val() * tau12.val()).exp() * (-0.5 * d1_root * d1_root).exp()
        / (2.0 * std::f64::consts::PI).sqrt()
        / (root * spec.vol.val() * tau12.val().sqrt());
    let root_ad = ift_lift_1d(
        root,
        |x| {
            black_scholes_price(
                x,
                k2,
                spec.vol,
                tau12,
                spec.rate,
                spec.div_yield,
                spec.daughter_type.is_call(),
            ) - k1
        },
        daughter_delta(T::from_f64(root), spec, tau12),
        fxx,
    );
    let t1 = spec.expiry_mother;
    let t2 = spec.expiry_daughter;
    let sqrt1 = t1.sqrt();
    let sqrt2 = t2.sqrt();
    let rho = (t1 / t2).sqrt();
    // `e1`/`e2` are the two critical-boundary d-values at the mother's
    // exercise date.  Writing them in the usual S/S* form is important: it
    // keeps the Geske put-call identities exact and gives a smooth lifted
    // derivative when the critical spot is obtained by IFT.
    let e1 = ((spec.spot / root_ad).ln()
        + (spec.rate - spec.div_yield + spec.vol * spec.vol * 0.5) * t1)
        / (spec.vol * sqrt1);
    let e2 = e1 - spec.vol * sqrt1;
    let dp = ((spec.spot / k2).ln()
        + (spec.rate - spec.div_yield + spec.vol * spec.vol * 0.5) * t2)
        / (spec.vol * sqrt2);
    let dm = dp - spec.vol * sqrt2;
    let (call_on_daughter, daughter_vanilla) = if spec.daughter_type.is_call() {
        (
            spec.spot * (-spec.div_yield * t2).exp() * bivariate_normal_cdf(dp, e1, rho)
                - k2 * (-spec.rate * t2).exp() * bivariate_normal_cdf(dm, e2, rho)
                - k1 * (-spec.rate * t1).exp() * e2.norm_cdf(),
            black_scholes_price(spec.spot, k2, spec.vol, t2, spec.rate, spec.div_yield, true),
        )
    } else {
        (
            k2 * (-spec.rate * t2).exp() * bivariate_normal_cdf(-dm, -e2, rho)
                - spec.spot * (-spec.div_yield * t2).exp() * bivariate_normal_cdf(-dp, -e1, rho)
                - k1 * (-spec.rate * t1).exp() * (-e2).norm_cdf(),
            black_scholes_price(
                spec.spot,
                k2,
                spec.vol,
                t2,
                spec.rate,
                spec.div_yield,
                false,
            ),
        )
    };
    if spec.mother_type.is_call() {
        call_on_daughter
    } else {
        // Geske's compound put-call parity is exact for either daughter:
        // C_on_X - P_on_X = X - K1 exp(-r T1).
        call_on_daughter - daughter_vanilla + k1 * (-spec.rate * t1).exp()
    }
}

mod evaluate_geske_compound_call_and_put_on_call_module {
    // Haug (2007) & Bouzoubaa–Osseiran (2010) Geske (1979) Compound Option
    // (`Call-on-Call` & `Put-on-Call`) Critical-Spot Inversion & Parity Engine.

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
            s_star = (s_star - (c_val - k1) / c_del.max_val(T::from_f64(1e-6)))
                .max_val(T::from_f64(1e-4));
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
    mod domain_suite_52_tests {
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
}
pub use evaluate_geske_compound_call_and_put_on_call_module::*;
