use ql_ad::{Dual, Real};
pub use ql_instruments::{BarrierType, OptionType};
use ql_vol::black_scholes_price;

pub const BGK_BETA: f64 = 0.5825971579390107;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MonitoringFrequency {
    Continuous,
    Discrete { observations: usize },
}

#[derive(Clone, Copy, Debug)]
pub struct BarrierOptionSpec<T: Real> {
    pub option_type: OptionType,
    pub barrier_type: BarrierType,
    pub spot: T,
    pub strike: T,
    pub barrier: T,
    pub rebate: T,
    pub rate: T,
    pub div_yield: T,
    pub vol: T,
    pub expiry: T,
    pub monitoring: MonitoringFrequency,
}

#[inline]
pub fn bgk_adjusted_barrier<T: Real>(
    spot: T,
    barrier: T,
    vol: T,
    expiry: T,
    observations: usize,
) -> T {
    if observations == 0 {
        return barrier;
    }
    let dt = expiry / observations as f64;
    let shift = vol * dt.sqrt() * BGK_BETA;
    if barrier.val() > spot.val() {
        barrier * shift.exp()
    } else if barrier.val() < spot.val() {
        barrier * (-shift).exp()
    } else {
        barrier
    }
}

fn price_continuous<T: Real>(spec: &BarrierOptionSpec<T>, h: T) -> T {
    let s = spec.spot;
    let k = spec.strike;
    let r = spec.rate;
    let q = spec.div_yield;
    let v = spec.vol;
    let t = spec.expiry;
    if t.val() <= 0.0 || v.val() <= 0.0 {
        return black_scholes_price(s, k, v, t, r, q, spec.option_type.is_call());
    }
    let b = r - q;
    let vv = v * v;
    let st = v * t.sqrt();
    let mu = (b - T::from_f64(0.5) * vv) / vv;
    let lambda = (mu * mu + T::from_f64(2.0) * r / vv).sqrt();
    let x1 = (s / k).ln() / st + (T::one() + mu) * st;
    let x2 = (s / h).ln() / st + (T::one() + mu) * st;
    let y1 = (h * h / (s * k)).ln() / st + (T::one() + mu) * st;
    let y2 = (h / s).ln() / st + (T::one() + mu) * st;
    let z = (h / s).ln() / st + lambda * st;
    let dr = (-r * t).exp();
    let dq = ((b - r) * t).exp();
    let ratio = h / s;
    let phi = if spec.option_type.is_call() {
        1.0
    } else {
        -1.0
    };
    let eta = if matches!(
        spec.barrier_type,
        BarrierType::DownIn | BarrierType::DownOut
    ) {
        1.0
    } else {
        -1.0
    };
    let a = T::from_f64(phi)
        * (s * dq * (T::from_f64(phi) * x1).norm_cdf()
            - k * dr * (T::from_f64(phi) * (x1 - st)).norm_cdf());
    let bb = T::from_f64(phi)
        * (s * dq * (T::from_f64(phi) * x2).norm_cdf()
            - k * dr * (T::from_f64(phi) * (x2 - st)).norm_cdf());
    let c = T::from_f64(phi)
        * (s * dq * ratio.pow((mu + 1.0) * 2.0) * (T::from_f64(eta) * y1).norm_cdf()
            - k * dr * ratio.pow(mu * 2.0) * (T::from_f64(eta) * (y1 - st)).norm_cdf());
    let d = T::from_f64(phi)
        * (s * dq * ratio.pow((mu + 1.0) * 2.0) * (T::from_f64(eta) * y2).norm_cdf()
            - k * dr * ratio.pow(mu * 2.0) * (T::from_f64(eta) * (y2 - st)).norm_cdf());
    let e = spec.rebate
        * dr
        * ((T::from_f64(eta) * (x2 - st)).norm_cdf()
            - ratio.pow(mu * 2.0) * (T::from_f64(eta) * (y2 - st)).norm_cdf());
    let f = spec.rebate
        * (ratio.pow(mu + lambda) * (T::from_f64(eta) * z).norm_cdf()
            + ratio.pow(mu - lambda)
                * (T::from_f64(eta) * (z - T::from_f64(2.0) * lambda * st)).norm_cdf());
    let down = matches!(
        spec.barrier_type,
        BarrierType::DownIn | BarrierType::DownOut
    );
    let in_barrier = matches!(spec.barrier_type, BarrierType::DownIn | BarrierType::UpIn);
    let ge = k.val() >= h.val();
    match (down, in_barrier, spec.option_type.is_call(), ge) {
        (true, true, true, true) => c + e,
        (true, true, true, false) => a - bb + d + e,
        (true, false, true, true) => a - c + f,
        (true, false, true, false) => bb - d + f,
        (false, true, true, true) => a + e,
        (false, true, true, false) => bb - c + d + e,
        (false, false, true, true) => f,
        (false, false, true, false) => a - bb + c - d + f,
        (true, true, false, true) => bb - c + d + e,
        (true, true, false, false) => a + e,
        (true, false, false, true) => a - bb + c - d + f,
        (true, false, false, false) => f,
        (false, true, false, true) => a - bb + d + e,
        (false, true, false, false) => c + e,
        (false, false, false, true) => bb - d + f,
        (false, false, false, false) => a - c + f,
    }
}

pub fn price_barrier_option<T: Real>(spec: &BarrierOptionSpec<T>) -> T {
    let adjusted = match spec.monitoring {
        MonitoringFrequency::Continuous => spec.barrier,
        MonitoringFrequency::Discrete { observations } => {
            bgk_adjusted_barrier(spec.spot, spec.barrier, spec.vol, spec.expiry, observations)
        }
    };
    let down = matches!(
        spec.barrier_type,
        BarrierType::DownIn | BarrierType::DownOut
    );
    // Contractual breach is tested against the contractual barrier.  The BGK
    // continuity correction only changes the analytic barrier used after the
    // initial-state check; using it here incorrectly knocks in/out a trade
    // that starts between the true and adjusted levels.
    let breached = if down {
        spec.spot.val() <= spec.barrier.val()
    } else {
        spec.spot.val() >= spec.barrier.val()
    };
    let knocked_in = matches!(spec.barrier_type, BarrierType::DownIn | BarrierType::UpIn);
    if breached {
        return if knocked_in {
            black_scholes_price(
                spec.spot,
                spec.strike,
                spec.vol,
                spec.expiry,
                spec.rate,
                spec.div_yield,
                spec.option_type.is_call(),
            )
        } else {
            spec.rebate
        };
    }
    let adjusted = if down {
        adjusted.min_val(spec.spot)
    } else {
        adjusted.max_val(spec.spot)
    };
    price_continuous(spec, adjusted)
}

pub fn barrier_option_greeks_dual5(spec: &BarrierOptionSpec<f64>) -> (f64, [f64; 5]) {
    let d = BarrierOptionSpec {
        option_type: spec.option_type,
        barrier_type: spec.barrier_type,
        spot: Dual::<5>::variable(spec.spot, 0),
        strike: Dual::from_value(spec.strike),
        barrier: Dual::from_value(spec.barrier),
        rebate: Dual::from_value(spec.rebate),
        rate: Dual::variable(spec.rate, 2),
        div_yield: Dual::variable(spec.div_yield, 3),
        vol: Dual::variable(spec.vol, 1),
        expiry: Dual::variable(spec.expiry, 4),
        monitoring: spec.monitoring,
    };
    let p = price_barrier_option(&d);
    (p.val, [p.d[0], p.d[1], p.d[2], p.d[3], -p.d[4]])
}

mod evaluate_carr_chou_put_call_symmetry_barrier_replication_module {
    // Carr & Chou (1997, JF) & Haug (2007, Ch. 4.17) Static Options Replication of
    // Down-and-Out Call via Put-Call Symmetry (`Geometric Reflection Strike H^2 / K`) Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct CarrChouPutCallSymmetryBarrierReplicationSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub strike_price_k: T,
        pub down_and_out_barrier_h: T,
        pub risk_free_rate_r: T,
        pub dividend_yield_q: T,
        pub volatility_sigma: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct CarrChouPutCallSymmetryBarrierReplicationReport<T: Real = f64> {
        pub vanilla_call_at_k: T,
        pub reflected_put_strike_h2_over_k: T,
        pub symmetry_reflection_ratio_weight: T,
        pub static_reflected_put_hedge_cost: T,
        pub static_replicated_down_and_out_call: T,
        pub boundary_residual_at_barrier: T,
    }

    fn bs_call<T: Real>(s: T, k: T, r: T, q: T, t: T, sig: T) -> T {
        let std = sig * t.sqrt();
        let d1 = ((s / k).ln() + (r - q + T::from_f64(0.5) * sig * sig) * t) / std;
        let d2 = d1 - std;
        s * (-q * t).exp() * norm_cdf(d1) - k * (-r * t).exp() * norm_cdf(d2)
    }

    fn bs_put<T: Real>(s: T, k: T, r: T, q: T, t: T, sig: T) -> T {
        let std = sig * t.sqrt();
        let d1 = ((s / k).ln() + (r - q + T::from_f64(0.5) * sig * sig) * t) / std;
        let d2 = d1 - std;
        k * (-r * t).exp() * norm_cdf(-d2) - s * (-q * t).exp() * norm_cdf(-d1)
    }

    pub fn evaluate_carr_chou_put_call_symmetry_barrier_replication<T: Real>(
        spec: &CarrChouPutCallSymmetryBarrierReplicationSpec<T>,
    ) -> CarrChouPutCallSymmetryBarrierReplicationReport<T> {
        let h_bar = spec.down_and_out_barrier_h.max_val(T::from_f64(1.0));
        let s0 = spec.spot_price_s0.max_val(h_bar + T::from_f64(0.5));
        let k = spec.strike_price_k.max_val(h_bar + T::from_f64(0.1));
        let r = spec.risk_free_rate_r;
        let q = spec.dividend_yield_q;
        let sig = spec.volatility_sigma.max_val(T::from_f64(0.02));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let c_van = bs_call(s0, k, r, q, t, sig);

        // Reiner-Rubinstein (1991) / Carr-Chou (1997) / Haug (2007, Ch. 4.17.1) reflection strike K_ref = H^2 / K
        // and general reflection weight w_sym(S0) = (K / H) * (H / S0)^{2(r - q) / sigma^2} with swapped carry (q, r):
        let k_ref = (h_bar * h_bar) / k;
        let drift_exp = (T::from_f64(2.0) * (r - q)) / (sig * sig);
        let sym_weight = (k / h_bar) * (drift_exp * (h_bar / s0).ln()).exp();

        let put_ref_s0 = bs_put(s0, k_ref, q, r, t, sig);
        let hedge_cost = sym_weight * put_ref_s0;
        let doc_val = (c_van - hedge_cost).max_val(T::zero());

        // Verify exact Put-Call Symmetry / Duality boundary matching at S = H (where (H / S)^{2(r - q)/sigma^2} = 1):
        let c_at_h = bs_call(h_bar, k, r, q, t, sig);
        let p_at_h = bs_put(h_bar, k_ref, q, r, t, sig);
        let bound_res = c_at_h - (k / h_bar) * p_at_h;

        CarrChouPutCallSymmetryBarrierReplicationReport {
            vanilla_call_at_k: c_van,
            reflected_put_strike_h2_over_k: k_ref,
            symmetry_reflection_ratio_weight: sym_weight,
            static_reflected_put_hedge_cost: hedge_cost,
            static_replicated_down_and_out_call: doc_val,
            boundary_residual_at_barrier: bound_res,
        }
    }

    #[cfg(test)]
    mod evaluate_carr_chou_put_call_symmetry_barrier_replication_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn carr_chou_put_call_symmetry_barrier_replication_invariants_and_ad_parity() {
            let spec = CarrChouPutCallSymmetryBarrierReplicationSpec {
                spot_price_s0: 100.0,
                strike_price_k: 102.0,
                down_and_out_barrier_h: 90.0,
                risk_free_rate_r: 0.045,
                dividend_yield_q: 0.020,
                volatility_sigma: 0.22,
                maturity_years_t: 0.75,
            };
            let rep = evaluate_carr_chou_put_call_symmetry_barrier_replication(&spec);
            assert!(rep.static_replicated_down_and_out_call > 0.0);
            assert!(rep.static_replicated_down_and_out_call < rep.vanilla_call_at_k);
            assert!(rep.boundary_residual_at_barrier.abs() < 1e-12);

            let h = 1e-6;
            let r_ad = evaluate_carr_chou_put_call_symmetry_barrier_replication(
                &CarrChouPutCallSymmetryBarrierReplicationSpec {
                    spot_price_s0: Dual::<7>::variable(spec.spot_price_s0, 0),
                    strike_price_k: Dual::<7>::variable(spec.strike_price_k, 1),
                    down_and_out_barrier_h: Dual::<7>::variable(spec.down_and_out_barrier_h, 2),
                    risk_free_rate_r: Dual::<7>::variable(spec.risk_free_rate_r, 3),
                    dividend_yield_q: Dual::<7>::variable(spec.dividend_yield_q, 4),
                    volatility_sigma: Dual::<7>::variable(spec.volatility_sigma, 5),
                    maturity_years_t: Dual::<7>::variable(spec.maturity_years_t, 6),
                },
            );
            for (idx, &ad) in r_ad
                .static_replicated_down_and_out_call
                .d
                .iter()
                .enumerate()
            {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_price_s0 += h;
                        sm.spot_price_s0 -= h;
                    }
                    1 => {
                        sp.strike_price_k += h;
                        sm.strike_price_k -= h;
                    }
                    2 => {
                        sp.down_and_out_barrier_h += h;
                        sm.down_and_out_barrier_h -= h;
                    }
                    3 => {
                        sp.risk_free_rate_r += h;
                        sm.risk_free_rate_r -= h;
                    }
                    4 => {
                        sp.dividend_yield_q += h;
                        sm.dividend_yield_q -= h;
                    }
                    5 => {
                        sp.volatility_sigma += h;
                        sm.volatility_sigma -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_carr_chou_put_call_symmetry_barrier_replication(&sp)
                    .static_replicated_down_and_out_call
                    - evaluate_carr_chou_put_call_symmetry_barrier_replication(&sm)
                        .static_replicated_down_and_out_call)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_carr_chou_put_call_symmetry_barrier_replication_module::*;

mod evaluate_conze_viswanathan_fixed_strike_lookback_excursion_module {
    // Conze & Viswanathan (1991, J. Finance, "Path Dependent Options: The Case of Lookback Options")
    // & Haug (2007, Ch. 4.15.1) Fixed-Strike Lookback Call Option & Running-Maximum Excursion Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct ConzeViswanathanFixedStrikeLookbackSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub fixed_strike_k: T,
        pub risk_free_rate_r: T,
        pub dividend_yield_q: T,
        pub volatility_sigma: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct ConzeViswanathanFixedStrikeLookbackReport<T: Real = f64> {
        pub standard_european_call_value: T,
        pub maximum_excursion_reflection_premium: T,
        pub fixed_strike_lookback_call_value: T,
        pub lookback_to_european_call_multiplier: T,
        pub reflection_scale_exponent_lambda: T,
    }

    pub fn evaluate_conze_viswanathan_fixed_strike_lookback_excursion<T: Real>(
        spec: &ConzeViswanathanFixedStrikeLookbackSpec<T>,
    ) -> ConzeViswanathanFixedStrikeLookbackReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1.0));
        let k = spec.fixed_strike_k.max_val(s0);
        let r = spec.risk_free_rate_r;
        let q = spec.dividend_yield_q;
        let sig = spec.volatility_sigma.max_val(T::from_f64(0.05));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let b = (r - q).max_val(T::from_f64(0.005));
        let sig2 = sig * sig;
        let std = sig * t.sqrt();

        let d1 = ((s0 / k).ln() + (b + T::from_f64(0.5) * sig2) * t) / std;
        let d2 = d1 - std;

        let df_r = (-r * t).exp();
        let df_q = (-q * t).exp();

        let euro_call = s0 * df_q * norm_cdf(d1) - k * df_r * norm_cdf(d2);

        // Conze & Viswanathan (1991, Eq. 9) / Haug (2007, Ch. 4.15.1) for K >= M_0:
        let lambda = (T::from_f64(2.0) * b) / sig2;
        let d3 = d1 - lambda * std;
        let excursion_prem = (s0 / lambda)
            * df_r
            * ((b * t).exp() * norm_cdf(d1) - (s0 / k).pow(-lambda) * norm_cdf(d3));
        let lookback_call = euro_call + excursion_prem;
        let mult = lookback_call / euro_call.max_val(T::from_f64(1e-4));

        ConzeViswanathanFixedStrikeLookbackReport {
            standard_european_call_value: euro_call,
            maximum_excursion_reflection_premium: excursion_prem,
            fixed_strike_lookback_call_value: lookback_call,
            lookback_to_european_call_multiplier: mult,
            reflection_scale_exponent_lambda: lambda,
        }
    }

    #[cfg(test)]
    mod evaluate_conze_viswanathan_fixed_strike_lookback_excursion_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn conze_viswanathan_fixed_strike_lookback_invariants_and_ad_parity() {
            let spec = ConzeViswanathanFixedStrikeLookbackSpec {
                spot_price_s0: 100.0,
                fixed_strike_k: 105.0,
                risk_free_rate_r: 0.050,
                dividend_yield_q: 0.015,
                volatility_sigma: 0.25,
                maturity_years_t: 1.0,
            };
            let rep = evaluate_conze_viswanathan_fixed_strike_lookback_excursion(&spec);
            assert!(rep.maximum_excursion_reflection_premium > 0.0);
            assert!(rep.fixed_strike_lookback_call_value > rep.standard_european_call_value);
            assert!(rep.lookback_to_european_call_multiplier > 1.5);

            let h = 1e-5;
            let r_ad = evaluate_conze_viswanathan_fixed_strike_lookback_excursion(
                &ConzeViswanathanFixedStrikeLookbackSpec {
                    spot_price_s0: Dual::<6>::variable(spec.spot_price_s0, 0),
                    fixed_strike_k: Dual::<6>::variable(spec.fixed_strike_k, 1),
                    risk_free_rate_r: Dual::<6>::variable(spec.risk_free_rate_r, 2),
                    dividend_yield_q: Dual::<6>::variable(spec.dividend_yield_q, 3),
                    volatility_sigma: Dual::<6>::variable(spec.volatility_sigma, 4),
                    maturity_years_t: Dual::<6>::variable(spec.maturity_years_t, 5),
                },
            );
            for (idx, &ad) in r_ad.fixed_strike_lookback_call_value.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_price_s0 += h;
                        sm.spot_price_s0 -= h;
                    }
                    1 => {
                        sp.fixed_strike_k += h;
                        sm.fixed_strike_k -= h;
                    }
                    2 => {
                        sp.risk_free_rate_r += h;
                        sm.risk_free_rate_r -= h;
                    }
                    3 => {
                        sp.dividend_yield_q += h;
                        sm.dividend_yield_q -= h;
                    }
                    4 => {
                        sp.volatility_sigma += h;
                        sm.volatility_sigma -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_conze_viswanathan_fixed_strike_lookback_excursion(&sp)
                    .fixed_strike_lookback_call_value
                    - evaluate_conze_viswanathan_fixed_strike_lookback_excursion(&sm)
                        .fixed_strike_lookback_call_value)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_conze_viswanathan_fixed_strike_lookback_excursion_module::*;

mod evaluate_broadie_glasserman_kou_discrete_lookback_option_module {
    // Broadie, Glasserman & Kou (1999, Finance and Stochastics 3(1), "Connecting Discrete and
    // Continuous Path-Dependent Options", Theorem 2), Goldman, Sosin & Gatto (1979, JF 34(5)),
    // & Haug (2007, The Complete Guide to Option Pricing Formulas, Ch. 4.16.4)
    // Continuity-Corrected Discrete Floating-Strike Lookback Call Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct BroadieGlassermanKouDiscreteLookbackSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub observed_running_minimum_m0: T,
        pub risk_free_rate_r: T,
        pub dividend_yield_q: T,
        pub volatility_sigma: T,
        pub discrete_monitoring_dates_n: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BroadieGlassermanKouDiscreteLookbackReport<T: Real = f64> {
        pub continuous_gsg_lookback_call_value: T,
        pub bgk_discrete_lookback_call_value: T,
        pub discrete_monitoring_price_discount: T,
        pub riemann_zeta_shift_factor: T,
        pub effective_shifted_minimum_m_eff: T,
    }

    fn gsg_floating_lookback_call_inner<T: Real>(s0: T, m: T, r: T, q: T, sig: T, t: T) -> T {
        let b = (r - q).max_val(T::from_f64(0.002));
        let var = sig * sig;
        let sqrt_t = t.sqrt();
        let std = sig * sqrt_t;
        let df_r = (-r * t).exp();
        let df_q = (-q * t).exp();

        let a1 = ((s0 / m).ln() + (b + T::from_f64(0.5) * var) * t) / std;
        let a2 = a1 - std;
        let a3 = -a1 + (T::from_f64(2.0) * b * sqrt_t) / sig;
        let lambda = (T::from_f64(2.0) * b) / var;
        let ratio_pow = (-lambda * (s0 / m).ln()).exp();

        s0 * df_q * norm_cdf(a1) - m * df_r * norm_cdf(a2)
            + (s0 / lambda) * (df_r * ratio_pow * norm_cdf(a3) - df_q * norm_cdf(-a1))
    }

    pub fn evaluate_broadie_glasserman_kou_discrete_lookback_option<T: Real>(
        spec: &BroadieGlassermanKouDiscreteLookbackSpec<T>,
    ) -> BroadieGlassermanKouDiscreteLookbackReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1.0));
        let m0 = spec
            .observed_running_minimum_m0
            .max_val(T::from_f64(0.5))
            .min_val(s0);
        let r = spec.risk_free_rate_r;
        let q = spec.dividend_yield_q;
        let sig = spec.volatility_sigma.max_val(T::from_f64(0.02));
        let n = spec.discrete_monitoring_dates_n.max_val(T::from_f64(2.0));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let c_cont = gsg_floating_lookback_call_inner(s0, m0, r, q, sig, t);

        // Broadie-Glasserman-Kou (1999) constant beta_1 = -zeta(1/2) / sqrt(2*pi) = 0.5825971579390107:
        let beta1 = T::from_f64(0.582_597_157_939_010_7);
        let dt = t / n;
        let shift = (-beta1 * sig * dt.sqrt()).exp();

        // Smooth minimum min(S0, m0 * shift) to keep C^infinity differentiability:
        // The continuity correction maps the observed discrete minimum to a
        // lower effective continuous minimum (BGK's beta_1 shift).
        let m_shifted = m0 * shift;
        let diff = s0 - m_shifted;
        let m_eff = T::from_f64(0.5) * (s0 + m_shifted - (diff * diff + T::from_f64(1e-6)).sqrt());

        let df_q = (-q * t).exp();
        let c_cont_shifted = gsg_floating_lookback_call_inner(s0, m_eff, r, q, sig, t);
        let c_disc = c_cont_shifted / shift + s0 * df_q * (T::one() - T::one() / shift);
        let discount = c_cont - c_disc;

        BroadieGlassermanKouDiscreteLookbackReport {
            continuous_gsg_lookback_call_value: c_cont,
            bgk_discrete_lookback_call_value: c_disc,
            discrete_monitoring_price_discount: discount,
            riemann_zeta_shift_factor: shift,
            effective_shifted_minimum_m_eff: m_eff,
        }
    }

    #[cfg(test)]
    mod evaluate_broadie_glasserman_kou_discrete_lookback_option_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bgk_discrete_lookback_invariants_and_ad_parity() {
            let spec = BroadieGlassermanKouDiscreteLookbackSpec {
                spot_price_s0: 100.0,
                observed_running_minimum_m0: 92.0,
                risk_free_rate_r: 0.05,
                dividend_yield_q: 0.015,
                volatility_sigma: 0.25,
                discrete_monitoring_dates_n: 52.0,
                maturity_years_t: 1.0,
            };
            let rep = evaluate_broadie_glasserman_kou_discrete_lookback_option(&spec);
            assert!(rep.bgk_discrete_lookback_call_value < rep.continuous_gsg_lookback_call_value);
            assert!(rep.discrete_monitoring_price_discount > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_broadie_glasserman_kou_discrete_lookback_option(
                &BroadieGlassermanKouDiscreteLookbackSpec {
                    spot_price_s0: Dual::<6>::variable(spec.spot_price_s0, 0),
                    observed_running_minimum_m0: Dual::<6>::variable(
                        spec.observed_running_minimum_m0,
                        1,
                    ),
                    risk_free_rate_r: Dual::<6>::variable(spec.risk_free_rate_r, 2),
                    dividend_yield_q: Dual::from_f64(spec.dividend_yield_q),
                    volatility_sigma: Dual::<6>::variable(spec.volatility_sigma, 3),
                    discrete_monitoring_dates_n: Dual::<6>::variable(
                        spec.discrete_monitoring_dates_n,
                        4,
                    ),
                    maturity_years_t: Dual::<6>::variable(spec.maturity_years_t, 5),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.spot_price_s0,
                    spec.observed_running_minimum_m0,
                    spec.risk_free_rate_r,
                    spec.volatility_sigma,
                    spec.discrete_monitoring_dates_n,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_broadie_glasserman_kou_discrete_lookback_option(
                    &BroadieGlassermanKouDiscreteLookbackSpec {
                        spot_price_s0: v[0],
                        observed_running_minimum_m0: v[1],
                        risk_free_rate_r: v[2],
                        volatility_sigma: v[3],
                        discrete_monitoring_dates_n: v[4],
                        maturity_years_t: v[5],
                        ..spec
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_broadie_glasserman_kou_discrete_lookback_option(
                    &BroadieGlassermanKouDiscreteLookbackSpec {
                        spot_price_s0: v[0],
                        observed_running_minimum_m0: v[1],
                        risk_free_rate_r: v[2],
                        volatility_sigma: v[3],
                        discrete_monitoring_dates_n: v[4],
                        maturity_years_t: v[5],
                        ..spec
                    },
                );
                let fd = (up.bgk_discrete_lookback_call_value
                    - dn.bgk_discrete_lookback_call_value)
                    / (2.0 * h);
                assert!(
                    (r_ad.bgk_discrete_lookback_call_value.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_broadie_glasserman_kou_discrete_lookback_option_module::*;

mod evaluate_rubinstein_reiner_down_and_out_call_with_rebate_module {
    // Rubinstein & Reiner (1991, Risk 4(8), "Breaking Down the Barriers"),
    // Rich (1994, Journal of Derivatives 2(1)), & Haug (2007, The Complete Guide to Option
    // Pricing Formulas, 2nd ed., Section 4.17.1)
    // Down-and-Out Call with First-Passage Hit Rebate & In-Out Reflection Parity Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct RubinsteinReinerDownAndOutCallRebateSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub strike_price_k: T,
        pub lower_barrier_h: T,
        pub cash_rebate_at_hit_r_rebate: T,
        pub risk_free_rate_r: T,
        pub dividend_yield_q: T,
        pub volatility_sigma: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct RubinsteinReinerDownAndOutCallRebateReport<T: Real = f64> {
        pub vanilla_black_scholes_call_a: T,
        pub down_and_in_reflection_call_c: T,
        pub pure_down_and_out_call_no_rebate: T,
        pub first_passage_hit_rebate_pv_f: T,
        pub total_down_and_out_call_with_rebate: T,
        pub barrier_to_spot_distance_ratio: T,
    }

    pub fn evaluate_rubinstein_reiner_down_and_out_call_with_rebate<T: Real>(
        spec: &RubinsteinReinerDownAndOutCallRebateSpec<T>,
    ) -> RubinsteinReinerDownAndOutCallRebateReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(10.0));
        let h_bar = spec
            .lower_barrier_h
            .max_val(T::from_f64(1.0))
            .min_val(s0 * T::from_f64(0.98));
        let k = spec.strike_price_k.max_val(h_bar + T::from_f64(0.5));
        let rebate = spec.cash_rebate_at_hit_r_rebate.max_val(T::zero());
        let r = spec.risk_free_rate_r.max_val(T::from_f64(0.001));
        let q = spec.dividend_yield_q;
        let sig = spec.volatility_sigma.max_val(T::from_f64(0.02));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let b = r - q;
        let var = sig * sig;
        let sqrt_t = t.sqrt();
        let std_t = sig * sqrt_t;

        // Rubinstein & Reiner (1991) / Haug (2007, Section 4.17.1) parameters:
        let mu = (b - T::from_f64(0.5) * var) / var;
        let lam = (mu * mu + (T::from_f64(2.0) * r) / var).sqrt();

        let df_r = (-r * t).exp();
        let df_q = ((b - r) * t).exp();
        let ln_h_over_s = (h_bar / s0).ln();

        // Term A: Vanilla Black-Scholes Call:
        let x1 = (s0 / k).ln() / std_t + (T::one() + mu) * std_t;
        let term_a = s0 * df_q * norm_cdf(x1) - k * df_r * norm_cdf(x1 - std_t);

        // Term C (for K > H): Down-and-In Reflection Call:
        let y1 = ((h_bar * h_bar) / (s0 * k)).ln() / std_t + (T::one() + mu) * std_t;
        let hs_2mu = (T::from_f64(2.0) * mu * ln_h_over_s).exp();
        let hs_2mup2 = (T::from_f64(2.0) * (T::one() + mu) * ln_h_over_s).exp();
        let term_c = s0 * df_q * hs_2mup2 * norm_cdf(y1) - k * df_r * hs_2mu * norm_cdf(y1 - std_t);

        let do_no_rebate = term_a - term_c;

        // Term F: First-Passage Laplace Transform Cash Rebate Paid at Hit (eta = +1 for down barrier):
        let z = ln_h_over_s / std_t + lam * std_t;
        let hs_mu_plus_lam = ((mu + lam) * ln_h_over_s).exp();
        let hs_mu_minus_lam = ((mu - lam) * ln_h_over_s).exp();
        let term_f = rebate
            * (hs_mu_plus_lam * norm_cdf(z)
                + hs_mu_minus_lam * norm_cdf(z - T::from_f64(2.0) * lam * std_t));

        let total_do = do_no_rebate + term_f;
        let dist_ratio = (s0 - h_bar) / s0;

        RubinsteinReinerDownAndOutCallRebateReport {
            vanilla_black_scholes_call_a: term_a,
            down_and_in_reflection_call_c: term_c,
            pure_down_and_out_call_no_rebate: do_no_rebate,
            first_passage_hit_rebate_pv_f: term_f,
            total_down_and_out_call_with_rebate: total_do,
            barrier_to_spot_distance_ratio: dist_ratio,
        }
    }

    #[cfg(test)]
    mod evaluate_rubinstein_reiner_down_and_out_call_with_rebate_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn rubinstein_reiner_down_and_out_call_rebate_invariants_and_ad_parity() {
            let spec = RubinsteinReinerDownAndOutCallRebateSpec {
                spot_price_s0: 100.0,
                strike_price_k: 100.0,
                lower_barrier_h: 90.0,
                cash_rebate_at_hit_r_rebate: 3.0,
                risk_free_rate_r: 0.05,
                dividend_yield_q: 0.02,
                volatility_sigma: 0.25,
                maturity_years_t: 1.0,
            };
            let rep = evaluate_rubinstein_reiner_down_and_out_call_with_rebate(&spec);
            assert!(
                (rep.pure_down_and_out_call_no_rebate + rep.down_and_in_reflection_call_c
                    - rep.vanilla_black_scholes_call_a)
                    .abs()
                    < 1e-12
            );
            assert!(rep.first_passage_hit_rebate_pv_f > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_rubinstein_reiner_down_and_out_call_with_rebate(
                &RubinsteinReinerDownAndOutCallRebateSpec {
                    spot_price_s0: Dual::<8>::variable(spec.spot_price_s0, 0),
                    strike_price_k: Dual::<8>::variable(spec.strike_price_k, 1),
                    lower_barrier_h: Dual::<8>::variable(spec.lower_barrier_h, 2),
                    cash_rebate_at_hit_r_rebate: Dual::<8>::variable(
                        spec.cash_rebate_at_hit_r_rebate,
                        3,
                    ),
                    risk_free_rate_r: Dual::<8>::variable(spec.risk_free_rate_r, 4),
                    dividend_yield_q: Dual::<8>::variable(spec.dividend_yield_q, 5),
                    volatility_sigma: Dual::<8>::variable(spec.volatility_sigma, 6),
                    maturity_years_t: Dual::<8>::variable(spec.maturity_years_t, 7),
                },
            );
            for i in 0..8 {
                let mut v = [
                    spec.spot_price_s0,
                    spec.strike_price_k,
                    spec.lower_barrier_h,
                    spec.cash_rebate_at_hit_r_rebate,
                    spec.risk_free_rate_r,
                    spec.dividend_yield_q,
                    spec.volatility_sigma,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_rubinstein_reiner_down_and_out_call_with_rebate(
                    &RubinsteinReinerDownAndOutCallRebateSpec {
                        spot_price_s0: v[0],
                        strike_price_k: v[1],
                        lower_barrier_h: v[2],
                        cash_rebate_at_hit_r_rebate: v[3],
                        risk_free_rate_r: v[4],
                        dividend_yield_q: v[5],
                        volatility_sigma: v[6],
                        maturity_years_t: v[7],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_rubinstein_reiner_down_and_out_call_with_rebate(
                    &RubinsteinReinerDownAndOutCallRebateSpec {
                        spot_price_s0: v[0],
                        strike_price_k: v[1],
                        lower_barrier_h: v[2],
                        cash_rebate_at_hit_r_rebate: v[3],
                        risk_free_rate_r: v[4],
                        dividend_yield_q: v[5],
                        volatility_sigma: v[6],
                        maturity_years_t: v[7],
                    },
                );
                let fd = (up.total_down_and_out_call_with_rebate
                    - dn.total_down_and_out_call_with_rebate)
                    / (2.0 * h);
                assert!(
                    (r_ad.total_down_and_out_call_with_rebate.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_rubinstein_reiner_down_and_out_call_with_rebate_module::*;

pub mod broadie_glasserman_kou_barrier_suite65_mod {
    use ql_ad::Real;
    use ql_math::norm_cdf;

    /// Specification for the Broadie–Glasserman–Kou (1997) discrete-monitoring barrier continuity
    /// correction engine.
    ///
    /// Reference: Broadie, Glasserman & Kou (1997, *Mathematical Finance* 7(4), "A Continuity
    /// Correction for Discrete Barrier Options"), Haug (2007, Ch. 4.17.6).
    #[derive(Debug, Clone, Copy)]
    pub struct BroadieGlassermanKouDiscreteBarrierSpec<T: Real> {
        /// Underlying spot price $S_0$.
        pub spot_price_s0: T,
        /// Strike price $K$.
        pub strike_price_k: T,
        /// Contractual lower knock-out barrier $H < S_0$.
        pub contractual_down_barrier_h: T,
        /// Number of discrete monitoring dates $m$ over $[0, T]$.
        pub discrete_monitoring_dates_m: T,
        /// Risk-free interest rate $r$.
        pub risk_free_rate_r: T,
        /// Continuous dividend yield $q$.
        pub dividend_yield_q: T,
        /// Volatility $\sigma$.
        pub volatility_sigma: T,
        /// Maturity $T$ (in years).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_broadie_glasserman_kou_discrete_barrier_continuity_correction`].
    #[derive(Debug, Clone, Copy)]
    pub struct BroadieGlassermanKouDiscreteBarrierAnalytics<T: Real> {
        /// BGK-shifted effective continuous barrier $H_m = H \exp(-\beta_{\text{BGK}} \sigma \sqrt{T/m})$.
        pub bgk_shifted_effective_barrier: T,
        /// Continuous-monitoring down-and-out call value $V_{\text{DO}}(S_0, K, H)$.
        pub continuous_down_and_out_call_price: T,
        /// Broadie–Glasserman–Kou discretely monitored down-and-out call value $V_{\text{DO}}(S_0, K, H_m)$.
        pub bgk_discrete_down_and_out_call_price: T,
        /// Discrete-monitoring survival valuation premium $V_{\text{DO}}(H_m) - V_{\text{DO}}(H) > 0$.
        pub discrete_monitoring_survival_premium: T,
    }

    /// Evaluates the Broadie–Glasserman–Kou (1997, Theorem 1) Riemann-zeta continuity correction
    /// $\beta_{\text{BGK}} = -\zeta(1/2)/\sqrt{2\pi} \approx 0.582597157939$ for a discretely monitored
    /// down-and-out call option ($K \ge H$).
    pub fn evaluate_broadie_glasserman_kou_discrete_barrier_continuity_correction<T: Real>(
        spec: &BroadieGlassermanKouDiscreteBarrierSpec<T>,
    ) -> BroadieGlassermanKouDiscreteBarrierAnalytics<T> {
        let beta_bgk = T::from_f64(0.5825971579390107);
        let t = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let m = spec.discrete_monitoring_dates_m.max_val(T::from_f64(1.0));
        let sig = spec.volatility_sigma.max_val(T::from_f64(1e-4));
        let sig2 = sig * sig;
        let sqrt_t = t.sqrt();
        let v_sqrt = sig * sqrt_t;

        let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
        let k = spec.strike_price_k.max_val(T::from_f64(1e-4));
        let h = spec.contractual_down_barrier_h.max_val(T::from_f64(1e-4));

        // BGK downward shift for a down-and-out barrier: H_m = H * exp(-beta_BGK * sigma * sqrt(T / m))
        let h_m = h * (-beta_bgk * sig * (t / m).sqrt()).exp();

        let r = spec.risk_free_rate_r;
        let q = spec.dividend_yield_q;
        let df_r = (-r * t).exp();
        let df_q = (-q * t).exp();

        // Standard Black-Scholes call C_BS(S_0, K)
        let d1 = ((s0 / k).ln() + (r - q + T::from_f64(0.5) * sig2) * t) / v_sqrt;
        let d2 = d1 - v_sqrt;
        let c_bs = s0 * df_q * norm_cdf(d1) - k * df_r * norm_cdf(d2);

        // Reflection exponent 2 * mu_drift / sigma^2 where mu_drift = r - q - 0.5 * sigma^2
        let two_mu = (T::from_f64(2.0) * (r - q)) / sig2 - T::one();

        let eval_do_at_barrier = |b: T| -> T {
            let s_refl = (b * b) / s0;
            let y1 = ((s_refl / k).ln() + (r - q + T::from_f64(0.5) * sig2) * t) / v_sqrt;
            let y2 = y1 - v_sqrt;
            let refl_weight = (two_mu * (b / s0).ln()).exp();
            let c_refl = refl_weight * (s_refl * df_q * norm_cdf(y1) - k * df_r * norm_cdf(y2));
            c_bs - c_refl
        };

        let v_cont = eval_do_at_barrier(h);
        let v_disc = eval_do_at_barrier(h_m);

        BroadieGlassermanKouDiscreteBarrierAnalytics {
            bgk_shifted_effective_barrier: h_m,
            continuous_down_and_out_call_price: v_cont,
            bgk_discrete_down_and_out_call_price: v_disc,
            discrete_monitoring_survival_premium: v_disc - v_cont,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bgk_discrete_barrier_invariants_and_ad_parity() {
            let base = BroadieGlassermanKouDiscreteBarrierSpec {
                spot_price_s0: 100.0_f64,
                strike_price_k: 100.0,
                contractual_down_barrier_h: 91.0,
                discrete_monitoring_dates_m: 52.0,
                risk_free_rate_r: 0.05,
                dividend_yield_q: 0.02,
                volatility_sigma: 0.25,
                maturity_years_t: 1.0,
            };
            let r0 = evaluate_broadie_glasserman_kou_discrete_barrier_continuity_correction(&base);
            assert!(r0.bgk_shifted_effective_barrier < base.contractual_down_barrier_h);
            assert!(
                r0.bgk_discrete_down_and_out_call_price > r0.continuous_down_and_out_call_price
            );
            assert!(r0.discrete_monitoring_survival_premium > 0.0);

            let inputs = [
                base.spot_price_s0,
                base.strike_price_k,
                base.contractual_down_barrier_h,
                base.discrete_monitoring_dates_m,
                base.risk_free_rate_r,
                base.dividend_yield_q,
                base.volatility_sigma,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BroadieGlassermanKouDiscreteBarrierSpec {
                spot_price_s0: dual_in[0],
                strike_price_k: dual_in[1],
                contractual_down_barrier_h: dual_in[2],
                discrete_monitoring_dates_m: dual_in[3],
                risk_free_rate_r: dual_in[4],
                dividend_yield_q: dual_in[5],
                volatility_sigma: dual_in[6],
                maturity_years_t: dual_in[7],
            };
            let r_ad =
                evaluate_broadie_glasserman_kou_discrete_barrier_continuity_correction(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| BroadieGlassermanKouDiscreteBarrierSpec {
                    spot_price_s0: v[0],
                    strike_price_k: v[1],
                    contractual_down_barrier_h: v[2],
                    discrete_monitoring_dates_m: v[3],
                    risk_free_rate_r: v[4],
                    dividend_yield_q: v[5],
                    volatility_sigma: v[6],
                    maturity_years_t: v[7],
                };
                let rp =
                    evaluate_broadie_glasserman_kou_discrete_barrier_continuity_correction(&mk(p));
                let rm =
                    evaluate_broadie_glasserman_kou_discrete_barrier_continuity_correction(&mk(m));
                let fd = (rp.bgk_discrete_down_and_out_call_price
                    - rm.bgk_discrete_down_and_out_call_price)
                    / (2.0 * eps);
                assert!(
                    r_ad.bgk_discrete_down_and_out_call_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.bgk_discrete_down_and_out_call_price.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use broadie_glasserman_kou_barrier_suite65_mod::*;

pub mod kunitomo_ikeda_curved_barrier_mod {
    use ql_ad::Real;
    use ql_math::norm_cdf;

    /// Specification for the Kunitomo & Ikeda (1992, Math. Finance 2(4)) / Haug (2007, Ch. 4.17.4)
    /// down-and-out call option with an exponentially curved barrier $H(t) = H_0 e^{\delta_H t}$.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct KunitomoIkedaCurvedBarrierSpec<T: Real> {
        /// Spot price $S_0 > H_0$.
        pub spot_price_s0: T,
        /// Option strike price $K > H_0 e^{\delta_H T}$.
        pub strike_price_k: T,
        /// Initial down barrier level $H_0 < S_0$ at $t=0$.
        pub initial_down_barrier_h0: T,
        /// Exponential barrier curvature growth rate $\delta_H$ ($H(t) = H_0 e^{\delta_H t}$).
        pub barrier_curvature_rate_delta_h: T,
        /// Continuous risk-free rate $r$.
        pub risk_free_rate_r: T,
        /// Continuous dividend yield $q$.
        pub dividend_yield_q: T,
        /// Black-Scholes volatility $\sigma > 0$.
        pub volatility_sigma: T,
        /// Option maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Kunitomo-Ikeda curved exponential barrier option engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct KunitomoIkedaCurvedBarrierAnalytics<T: Real> {
        /// Terminal curved barrier level $H(T) = H_0 e^{\delta_H T}$.
        pub terminal_curved_barrier_level: T,
        /// Curvature-adjusted reflection exponent $\lambda_H = (r - q - \delta_H)/\sigma^2 - 1/2$.
        pub barrier_curvature_reflection_exponent: T,
        /// Knock-out reflection discount wedge relative to vanilla Black-Scholes call.
        pub curved_barrier_knockout_discount: T,
        /// Kunitomo-Ikeda down-and-out call option price.
        pub kunitomo_ikeda_down_and_out_call_price: T,
    }

    /// Evaluates the Kunitomo & Ikeda (1992) closed-form price of a down-and-out call option
    /// with an exponentially curved boundary $H(t) = H_0 e^{\delta_H t}$.
    pub fn evaluate_kunitomo_ikeda_curved_exponential_barrier_option<T: Real>(
        spec: &KunitomoIkedaCurvedBarrierSpec<T>,
    ) -> KunitomoIkedaCurvedBarrierAnalytics<T> {
        let t = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
        let k = spec.strike_price_k.max_val(T::from_f64(1e-4));
        let h0 = spec
            .initial_down_barrier_h0
            .max_val(T::from_f64(1e-4))
            .min_val(s0 * T::from_f64(0.99));
        let sig = spec.volatility_sigma.max_val(T::from_f64(1e-4));
        let sig2 = sig * sig;
        let v_sqrt = sig * t.sqrt();

        let ht = h0 * (spec.barrier_curvature_rate_delta_h * t).exp();
        let b = spec.risk_free_rate_r - spec.dividend_yield_q;
        let lam_h = (b - spec.barrier_curvature_rate_delta_h) / sig2 - T::from_f64(0.5);

        let d1 = ((s0 / k).ln() + (b + T::from_f64(0.5) * sig2) * t) / v_sqrt;
        let d2 = d1 - v_sqrt;
        let df_r = (-spec.risk_free_rate_r * t).exp();
        let df_q = (-spec.dividend_yield_q * t).exp();
        let vanilla = s0 * df_q * norm_cdf(d1) - k * df_r * norm_cdf(d2);

        // Reflection image arguments under curved barrier X_t = ln(S_t / H(t)):
        // In ln S_T = X_T + ln H_0 + delta_H * T, the +delta_H * T boundary shift cancels
        // the -delta_H * T drift shift in mu_X * T, ensuring y1 = d1 and V(H_0, 0) = 0 at S_0 = H_0.
        let y1 = ((h0 * h0 / (s0 * k)).ln() + (b + T::from_f64(0.5) * sig2) * t) / v_sqrt;
        let y2 = y1 - v_sqrt;

        let ratio = h0 / s0;
        let pow_s = ratio.pow(T::from_f64(2.0) * (lam_h + T::one()));
        let pow_k = ratio.pow(T::from_f64(2.0) * lam_h);

        let knock_out_reflection =
            s0 * df_q * pow_s * norm_cdf(y1) - k * df_r * pow_k * norm_cdf(y2);
        let raw_do = vanilla - knock_out_reflection;
        let do_call = T::from_f64(0.5) * (raw_do + (raw_do * raw_do + T::from_f64(1e-12)).sqrt());

        KunitomoIkedaCurvedBarrierAnalytics {
            terminal_curved_barrier_level: ht,
            barrier_curvature_reflection_exponent: lam_h,
            curved_barrier_knockout_discount: knock_out_reflection,
            kunitomo_ikeda_down_and_out_call_price: do_call,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn kunitomo_ikeda_curved_barrier_invariants_and_ad_parity() {
            let base = KunitomoIkedaCurvedBarrierSpec {
                spot_price_s0: 100.0_f64,
                strike_price_k: 100.0,
                initial_down_barrier_h0: 88.0,
                barrier_curvature_rate_delta_h: 0.02,
                risk_free_rate_r: 0.05,
                dividend_yield_q: 0.015,
                volatility_sigma: 0.25,
                maturity_years_t: 1.0,
            };
            let r0 = evaluate_kunitomo_ikeda_curved_exponential_barrier_option(&base);
            assert!(r0.terminal_curved_barrier_level > base.initial_down_barrier_h0);
            assert!(r0.curved_barrier_knockout_discount > 0.0);
            assert!(r0.kunitomo_ikeda_down_and_out_call_price > 0.0);

            let inputs = [
                base.spot_price_s0,
                base.strike_price_k,
                base.initial_down_barrier_h0,
                base.barrier_curvature_rate_delta_h,
                base.risk_free_rate_r,
                base.dividend_yield_q,
                base.volatility_sigma,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = KunitomoIkedaCurvedBarrierSpec {
                spot_price_s0: dual_in[0],
                strike_price_k: dual_in[1],
                initial_down_barrier_h0: dual_in[2],
                barrier_curvature_rate_delta_h: dual_in[3],
                risk_free_rate_r: dual_in[4],
                dividend_yield_q: dual_in[5],
                volatility_sigma: dual_in[6],
                maturity_years_t: dual_in[7],
            };
            let r_ad = evaluate_kunitomo_ikeda_curved_exponential_barrier_option(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| KunitomoIkedaCurvedBarrierSpec {
                    spot_price_s0: v[0],
                    strike_price_k: v[1],
                    initial_down_barrier_h0: v[2],
                    barrier_curvature_rate_delta_h: v[3],
                    risk_free_rate_r: v[4],
                    dividend_yield_q: v[5],
                    volatility_sigma: v[6],
                    maturity_years_t: v[7],
                };
                let rp = evaluate_kunitomo_ikeda_curved_exponential_barrier_option(&mk(p));
                let rm = evaluate_kunitomo_ikeda_curved_exponential_barrier_option(&mk(m));
                let fd = (rp.kunitomo_ikeda_down_and_out_call_price
                    - rm.kunitomo_ikeda_down_and_out_call_price)
                    / (2.0 * eps);
                assert!(
                    r_ad.kunitomo_ikeda_down_and_out_call_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.kunitomo_ikeda_down_and_out_call_price.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use kunitomo_ikeda_curved_barrier_mod::*;

// ============================================================================
// Hart & Ross (1994) & Haug (2007, Ch. 4.17.3) Soft Proportional Band Barrier Option
// ============================================================================

pub mod hart_ross_soft_barrier_mod {
    use ql_ad::Real;
    use ql_math::distributions::norm_cdf;

    /// Specification for the Hart & Ross (1994, *Risk* 7(5)) / Haug (2007, Section 4.17.3)
    /// soft down-and-out call option knocking out proportionally across $[L, U]$.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HartRossSoftBarrierSpec<T: Real> {
        /// Current spot price $S_0 > 0$.
        pub spot_price_s0: T,
        /// Option strike price $K > 0$.
        pub strike_price_k: T,
        /// Lower soft-barrier edge $L > 0$ (100% knocked out below $L$).
        pub lower_soft_barrier_l: T,
        /// Upper soft-barrier edge $U \in (L, S_0)$ (knock-out begins below $U$).
        pub upper_soft_barrier_u: T,
        /// Risk-free interest rate $r > 0$.
        pub risk_free_rate_r: T,
        /// Continuous dividend yield $q \ge 0$.
        pub dividend_yield_q: T,
        /// Volatility $\sigma > 0$.
        pub volatility_sigma: T,
        /// Option maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Hart-Ross soft proportional band barrier option engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HartRossSoftBarrierAnalytics<T: Real> {
        /// European Black-Scholes vanilla call price $C_{\text{BS}}$.
        pub european_vanilla_call_price: T,
        /// Soft down-and-in call value $C_{\text{soft-DI}}(L, U)$.
        pub soft_down_and_in_call_value: T,
        /// Effective midpoint barrier level $H_{\text{mid}} = \frac{1}{2}(L + U)$.
        pub effective_band_midpoint_barrier: T,
        /// Soft down-and-out call price $C_{\text{soft-DO}}(L, U) = C_{\text{BS}} - C_{\text{soft-DI}}(L, U)$.
        pub soft_down_and_out_call_price: T,
    }

    /// Evaluates the Hart & Ross (1994) / Haug (2007, Ch. 4.17.3) soft down-and-out call price
    /// by integrating the down-and-in reflection formula across the proportional knock-out band $[L, U]$.
    pub fn evaluate_hart_ross_haug_soft_proportional_band_barrier_option<T: Real>(
        spec: &HartRossSoftBarrierSpec<T>,
    ) -> HartRossSoftBarrierAnalytics<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
        let k = spec.strike_price_k.max_val(T::from_f64(1e-4));
        let l_bar = spec
            .lower_soft_barrier_l
            .max_val(T::from_f64(1e-4))
            .min_val(s0 * T::from_f64(0.95));
        let u_bar = spec
            .upper_soft_barrier_u
            .max_val(l_bar + T::from_f64(0.5))
            .min_val(s0 * T::from_f64(0.99));

        let sig = spec.volatility_sigma.max_val(T::from_f64(1e-4));
        let sig2 = sig * sig;
        let t = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let v_sqrt = sig * t.sqrt();
        let b = spec.risk_free_rate_r - spec.dividend_yield_q;
        let mu = (b + T::from_f64(0.5) * sig2) / sig2;

        let d1 = ((s0 / k).ln() + (b + T::from_f64(0.5) * sig2) * t) / v_sqrt;
        let d2 = d1 - v_sqrt;
        let df_r = (-spec.risk_free_rate_r * t).exp();
        let df_q = (-spec.dividend_yield_q * t).exp();
        let vanilla = s0 * df_q * norm_cdf(d1) - k * df_r * norm_cdf(d2);

        // 4-point Gauss-Legendre integration of C_DI(H) over H in [L, U]
        let gl_nodes = [
            0.069_431_844_202_973_71,
            0.330_009_478_207_571_87,
            0.669_990_521_792_428_1,
            0.930_568_155_797_026_2,
        ];
        let gl_weights = [
            0.173_927_422_568_726_92,
            0.326_072_577_431_273_1,
            0.326_072_577_431_273_1,
            0.173_927_422_568_726_92,
        ];

        let mut soft_di = T::zero();
        for (&node_f64, &w_f64) in gl_nodes.iter().zip(gl_weights.iter()) {
            let s = T::from_f64(node_f64);
            let w = T::from_f64(w_f64);
            let h = l_bar + s * (u_bar - l_bar);

            let y1 = ((h * h / (s0 * k)).ln() + (b + T::from_f64(0.5) * sig2) * t) / v_sqrt;
            let y2 = y1 - v_sqrt;
            let ratio = h / s0;
            let pow_s = ratio.pow(T::from_f64(2.0) * mu);
            let pow_k = ratio.pow(T::from_f64(2.0) * (mu - T::one()));
            let di_h = s0 * df_q * pow_s * norm_cdf(y1) - k * df_r * pow_k * norm_cdf(y2);
            soft_di += w * di_h;
        }

        let raw_do = vanilla - soft_di;
        let soft_do = T::from_f64(0.5) * (raw_do + (raw_do * raw_do + T::from_f64(1e-12)).sqrt());

        HartRossSoftBarrierAnalytics {
            european_vanilla_call_price: vanilla,
            soft_down_and_in_call_value: soft_di,
            effective_band_midpoint_barrier: T::from_f64(0.5) * (l_bar + u_bar),
            soft_down_and_out_call_price: soft_do,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hart_ross_soft_barrier_invariants_and_ad_parity() {
            let base = HartRossSoftBarrierSpec {
                spot_price_s0: 100.0_f64,
                strike_price_k: 100.0,
                lower_soft_barrier_l: 80.0,
                upper_soft_barrier_u: 90.0,
                risk_free_rate_r: 0.05,
                dividend_yield_q: 0.015,
                volatility_sigma: 0.24,
                maturity_years_t: 1.0,
            };
            let r0 = evaluate_hart_ross_haug_soft_proportional_band_barrier_option(&base);
            assert!(r0.soft_down_and_in_call_value > 0.0);
            assert!(r0.soft_down_and_out_call_price > 0.0);
            assert!(r0.soft_down_and_out_call_price < r0.european_vanilla_call_price);

            let inputs = [
                base.spot_price_s0,
                base.strike_price_k,
                base.lower_soft_barrier_l,
                base.upper_soft_barrier_u,
                base.risk_free_rate_r,
                base.dividend_yield_q,
                base.volatility_sigma,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HartRossSoftBarrierSpec {
                spot_price_s0: dual_in[0],
                strike_price_k: dual_in[1],
                lower_soft_barrier_l: dual_in[2],
                upper_soft_barrier_u: dual_in[3],
                risk_free_rate_r: dual_in[4],
                dividend_yield_q: dual_in[5],
                volatility_sigma: dual_in[6],
                maturity_years_t: dual_in[7],
            };
            let r_ad = evaluate_hart_ross_haug_soft_proportional_band_barrier_option(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| HartRossSoftBarrierSpec {
                    spot_price_s0: v[0],
                    strike_price_k: v[1],
                    lower_soft_barrier_l: v[2],
                    upper_soft_barrier_u: v[3],
                    risk_free_rate_r: v[4],
                    dividend_yield_q: v[5],
                    volatility_sigma: v[6],
                    maturity_years_t: v[7],
                };
                let rp = evaluate_hart_ross_haug_soft_proportional_band_barrier_option(&mk(p));
                let rm = evaluate_hart_ross_haug_soft_proportional_band_barrier_option(&mk(m));
                let fd = (rp.soft_down_and_out_call_price - rm.soft_down_and_out_call_price)
                    / (2.0 * eps);
                assert!(
                    r_ad.soft_down_and_out_call_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.soft_down_and_out_call_price.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use hart_ross_soft_barrier_mod::*;
mod hui_haug_binary_barrier_mod {
    use ql_ad::Real;

    /// Input specification for the Hui (1996, *Applied Financial Economics* 6(4)) /
    /// Reiner & Rubinstein (1991, *Risk* 4(9)) / Haug (2007, §4.19.1–§4.19.4) binary barrier
    /// options (down-and-in cash-at-hit One-Touch vs. down-and-out cash-at-expiry No-Touch).
    #[derive(Debug, Clone, Copy)]
    pub struct HuiHaugBinaryBarrierSpec<T: Real> {
        /// Current spot price $S_0$ (e.g., $105.0$).
        pub spot_price_s0: T,
        /// Lower barrier $H < S_0$ (e.g., $92.0$).
        pub down_barrier_h: T,
        /// Cash rebate $K_{\text{cash}} > 0$ paid upon hit / survival (e.g., $10.0$).
        pub cash_rebate_k: T,
        /// Risk-free interest rate $r$ (e.g., $0.045$).
        pub risk_free_rate_r: T,
        /// Continuous dividend / foreign yield $q$ (e.g., $0.015$).
        pub dividend_yield_q: T,
        /// Volatility $\sigma$ (e.g., $0.22$).
        pub volatility_sigma: T,
        /// Option maturity $T$ in years (e.g., $1.0$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_hui_haug_one_touch_and_no_touch_binary_barrier`].
    #[derive(Debug, Clone, Copy)]
    pub struct HuiHaugBinaryBarrierAnalytics<T: Real> {
        /// Risk-neutral survival (no-touch) probability $P_{\text{NT}}(0, T)$ over $[0, T]$.
        pub risk_neutral_no_touch_probability: T,
        /// Down-and-out cash-at-expiry No-Touch value $K_{\text{cash}} e^{-r T} P_{\text{NT}}(0, T)$.
        pub cash_at_expiry_no_touch_value: T,
        /// Down-and-in cash-at-expiry One-Touch value $K_{\text{cash}} e^{-r T} (1 - P_{\text{NT}}(0, T))$.
        pub cash_at_expiry_one_touch_value: T,
        /// Down-and-in cash-at-hit One-Touch value $K_{\text{cash}} \mathbb{E}^{\mathbb{Q}}[e^{-r\tau_H}\mathbf{1}_{\{\tau_H \le T\}}]$.
        pub cash_at_hit_one_touch_value: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        T::from_f64(0.5) * (T::one() + (x * T::from_f64(0.707_106_781_186_547_5)).erf())
    }

    /// Evaluates the exact Reiner & Rubinstein (1991) / Hui (1996) / Haug (2007, §4.19.1–§4.19.4)
    /// binary down-and-out No-Touch (cash-at-expiry), down-and-in One-Touch (cash-at-expiry), and
    /// down-and-in One-Touch (cash-at-hit) options, satisfying exact discounted rebate conservation
    /// $V_{\text{NT,exp}} + V_{\text{OT,exp}} \equiv K_{\text{cash}} e^{-r T}$ and early-hit discounting gain
    /// $V_{\text{OT,hit}} > V_{\text{OT,exp}}$ for $r > 0$.
    pub fn evaluate_hui_haug_one_touch_and_no_touch_binary_barrier<T: Real>(
        spec: &HuiHaugBinaryBarrierSpec<T>,
    ) -> HuiHaugBinaryBarrierAnalytics<T> {
        let s = spec.spot_price_s0;
        let h = spec.down_barrier_h;
        let k_cash = spec.cash_rebate_k;
        let r = spec.risk_free_rate_r;
        let q = spec.dividend_yield_q;
        let b = r - q;
        let sig = spec.volatility_sigma;
        let sig2 = sig * sig;
        let t = spec.maturity_years_t;
        let std_dev = sig * t.sqrt();

        // Drift parameter mu = (b - 0.5*sig^2) / sig^2 (Haug 2007 §4.19 notation)
        let mu = (b - T::from_f64(0.5) * sig2) / sig2;
        let log_s_h = (s / h).ln();
        let x2 = log_s_h / std_dev + (T::one() + mu) * std_dev;
        let y2 = -log_s_h / std_dev + (T::one() + mu) * std_dev;

        // Risk-neutral survival probability P(min_{0<=u<=T} S_u > H) = Phi(x2 - std_dev) - (H/S)^(2*mu) * Phi(y2 - std_dev)
        let d2_x = x2 - std_dev;
        let d2_y = y2 - std_dev;
        let h_over_s = h / s;
        let p_nt = std_norm_cdf(d2_x) - h_over_s.pow(T::from_f64(2.0) * mu) * std_norm_cdf(d2_y);

        let df = (-r * t).exp();
        let v_nt_exp = k_cash * df * p_nt;
        let v_ot_exp = k_cash * df * (T::one() - p_nt);

        // Cash-at-hit One-Touch (Haug §4.19.1, formula 1):
        // lam_hit = sqrt(mu^2 + 2*r / sig^2)
        // z = ln(H/S) / (sig*sqrt(T)) + lam_hit * sig * sqrt(T)
        // V_OT_hit = K * [ (H/S)^(mu + lam_hit) * Phi(z) + (H/S)^(mu - lam_hit) * Phi(z - 2*lam_hit*sig*sqrt(T)) ]
        let lam_hit = (mu * mu + T::from_f64(2.0) * r / sig2).sqrt();
        let z_hit = -log_s_h / std_dev + lam_hit * std_dev;
        let z_hit_minus = z_hit - T::from_f64(2.0) * lam_hit * std_dev;

        let v_ot_hit = k_cash
            * (h_over_s.pow(mu + lam_hit) * std_norm_cdf(z_hit)
                + h_over_s.pow(mu - lam_hit) * std_norm_cdf(z_hit_minus));

        HuiHaugBinaryBarrierAnalytics {
            risk_neutral_no_touch_probability: p_nt,
            cash_at_expiry_no_touch_value: v_nt_exp,
            cash_at_expiry_one_touch_value: v_ot_exp,
            cash_at_hit_one_touch_value: v_ot_hit,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hui_haug_binary_barrier_invariants_and_ad_parity() {
            let base = HuiHaugBinaryBarrierSpec {
                spot_price_s0: 105.0_f64,
                down_barrier_h: 92.0_f64,
                cash_rebate_k: 10.0_f64,
                risk_free_rate_r: 0.045_f64,
                dividend_yield_q: 0.015_f64,
                volatility_sigma: 0.22_f64,
                maturity_years_t: 1.0_f64,
            };
            let res = evaluate_hui_haug_one_touch_and_no_touch_binary_barrier(&base);
            assert!(
                res.risk_neutral_no_touch_probability > 0.0
                    && res.risk_neutral_no_touch_probability < 1.0
            );
            let df = (-base.risk_free_rate_r * base.maturity_years_t).exp();
            assert!(
                (res.cash_at_expiry_no_touch_value + res.cash_at_expiry_one_touch_value
                    - base.cash_rebate_k * df)
                    .abs()
                    < 1e-12
            );
            assert!(res.cash_at_hit_one_touch_value > res.cash_at_expiry_one_touch_value);

            let inputs = [
                base.spot_price_s0,
                base.down_barrier_h,
                base.cash_rebate_k,
                base.risk_free_rate_r,
                base.dividend_yield_q,
                base.volatility_sigma,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HuiHaugBinaryBarrierSpec {
                spot_price_s0: dual_in[0],
                down_barrier_h: dual_in[1],
                cash_rebate_k: dual_in[2],
                risk_free_rate_r: dual_in[3],
                dividend_yield_q: dual_in[4],
                volatility_sigma: dual_in[5],
                maturity_years_t: dual_in[6],
            };
            let r_ad = evaluate_hui_haug_one_touch_and_no_touch_binary_barrier(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| HuiHaugBinaryBarrierSpec {
                    spot_price_s0: v[0],
                    down_barrier_h: v[1],
                    cash_rebate_k: v[2],
                    risk_free_rate_r: v[3],
                    dividend_yield_q: v[4],
                    volatility_sigma: v[5],
                    maturity_years_t: v[6],
                };
                let rp = evaluate_hui_haug_one_touch_and_no_touch_binary_barrier(&mk(p));
                let rm = evaluate_hui_haug_one_touch_and_no_touch_binary_barrier(&mk(m));
                let fd =
                    (rp.cash_at_hit_one_touch_value - rm.cash_at_hit_one_touch_value) / (2.0 * eps);
                assert!(
                    r_ad.cash_at_hit_one_touch_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.cash_at_hit_one_touch_value.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use hui_haug_binary_barrier_mod::*;
mod chesney_jeanblanc_yor_parisian_mod {
    use ql_ad::Real;

    /// Input specification for the Chesney, Jeanblanc-Picqué & Yor (1997, *Advances in
    /// Applied Probability* 29(1)) & Hugonnier (1999) Parisian down-and-out call option
    /// where knock-out requires a continuous excursion below $H$ of duration $D_{\text{excur}} > 0$.
    #[derive(Debug, Clone, Copy)]
    pub struct ChesneyJeanblancYorParisianBarrierSpec<T: Real> {
        /// Current spot price $S_0 > H$ (e.g., $105.0$).
        pub spot_price_s0: T,
        /// Option strike $K$ (e.g., $100.0$).
        pub strike_k: T,
        /// Lower Parisian barrier $H < S_0$ (e.g., $92.0$).
        pub down_barrier_h: T,
        /// Required continuous excursion window $D_{\text{excur}} > 0$ in years (e.g., $0.06$).
        pub excursion_window_years_d: T,
        /// Risk-free interest rate $r$ (e.g., $0.042$).
        pub risk_free_rate_r: T,
        /// Continuous dividend yield $q$ (e.g., $0.015$).
        pub dividend_yield_q: T,
        /// Volatility $\sigma$ (e.g., $0.24$).
        pub volatility_sigma: T,
        /// Option maturity $T > D_{\text{excur}}$ in years (e.g., $1.0$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_chesney_jeanblanc_yor_parisian_excursion_barrier_call`].
    #[derive(Debug, Clone, Copy)]
    pub struct ChesneyJeanblancYorParisianBarrierAnalytics<T: Real> {
        /// Standard instantaneous-knockout down-and-out call price $C_{\text{DOC}}(S_0, K, H, T)$ ($D = 0$).
        pub instantaneous_down_and_out_call: T,
        /// Parisian excursion-shifted effective barrier $H_{\text{Par}}(D) = H \exp(-\sigma \sqrt{\frac{\pi D}{2}}) < H$.
        pub effective_parisian_shifted_barrier: T,
        /// Parisian excursion survival gain $C_{\text{Parisian}} - C_{\text{DOC}} > 0$.
        pub parisian_excursion_survival_premium: T,
        /// Total Chesney–Jeanblanc–Yor Parisian down-and-out call value.
        pub parisian_down_and_out_call_value: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        T::from_f64(0.5) * (T::one() + (x * T::from_f64(0.707_106_781_186_547_5)).erf())
    }

    #[inline]
    fn doc_price<T: Real>(s0: T, k: T, h: T, r: T, q: T, sig: T, t: T) -> T {
        let sig2 = sig * sig;
        let std_dev = sig * t.sqrt();
        let half = T::from_f64(0.5);
        let d1 = ((s0 / k).ln() + (r - q + half * sig2) * t) / std_dev;
        let d2 = d1 - std_dev;
        let c_van = s0 * (-q * t).exp() * std_norm_cdf(d1) - k * (-r * t).exp() * std_norm_cdf(d2);

        let lam = (r - q + half * sig2) / sig2;
        let y = ((h * h) / (s0 * k)).ln() / std_dev + lam * std_dev;
        let h_over_s = h / s0;
        let c_di = s0 * (-q * t).exp() * h_over_s.pow(T::from_f64(2.0) * lam) * std_norm_cdf(y)
            - k * (-r * t).exp()
                * h_over_s.pow(T::from_f64(2.0) * lam - T::from_f64(2.0))
                * std_norm_cdf(y - std_dev);
        c_van - c_di
    }

    /// Evaluates the Chesney, Jeanblanc-Picqué & Yor (1997, Proposition 4) & Anderluh–van Weeren
    /// Parisian down-and-out call option via the Laplace-excursion shifted barrier
    /// $H_{\text{Par}}(D) = H \exp\!\left(-\sigma\sqrt{\frac{\pi D}{2}}\right)$, which satisfies
    /// $C_{\text{DOC}}(H) < C_{\text{Parisian}}(H, D) < C_{\text{BS}}$ for all $D > 0$.
    pub fn evaluate_chesney_jeanblanc_yor_parisian_excursion_barrier_call<T: Real>(
        spec: &ChesneyJeanblancYorParisianBarrierSpec<T>,
    ) -> ChesneyJeanblancYorParisianBarrierAnalytics<T> {
        let s0 = spec.spot_price_s0;
        let k = spec.strike_k;
        let h = spec.down_barrier_h;
        let d_win = spec.excursion_window_years_d;
        let r = spec.risk_free_rate_r;
        let q = spec.dividend_yield_q;
        let sig = spec.volatility_sigma;
        let t = spec.maturity_years_t;

        let c_doc = doc_price(s0, k, h, r, q, sig, t);

        // Chesney-Jeanblanc-Yor Brownian meander excursion shift sqrt(pi / 2) * sigma * sqrt(D)
        let sqrt_half_pi = T::from_f64(1.253_314_137_315_500_3);
        let h_par = h * (-sqrt_half_pi * sig * d_win.sqrt()).exp();
        let c_par = doc_price(s0, k, h_par, r, q, sig, t);
        let surv_prem = c_par - c_doc;

        ChesneyJeanblancYorParisianBarrierAnalytics {
            instantaneous_down_and_out_call: c_doc,
            effective_parisian_shifted_barrier: h_par,
            parisian_excursion_survival_premium: surv_prem,
            parisian_down_and_out_call_value: c_par,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn chesney_jeanblanc_yor_parisian_invariants_and_ad_parity() {
            let base = ChesneyJeanblancYorParisianBarrierSpec {
                spot_price_s0: 105.0_f64,
                strike_k: 100.0_f64,
                down_barrier_h: 92.0_f64,
                excursion_window_years_d: 0.06_f64,
                risk_free_rate_r: 0.042_f64,
                dividend_yield_q: 0.015_f64,
                volatility_sigma: 0.24_f64,
                maturity_years_t: 1.0_f64,
            };
            let res = evaluate_chesney_jeanblanc_yor_parisian_excursion_barrier_call(&base);
            assert!(res.effective_parisian_shifted_barrier < base.down_barrier_h);
            assert!(res.parisian_excursion_survival_premium > 0.0);
            assert!(res.parisian_down_and_out_call_value > res.instantaneous_down_and_out_call);

            let inputs = [
                base.spot_price_s0,
                base.strike_k,
                base.down_barrier_h,
                base.excursion_window_years_d,
                base.risk_free_rate_r,
                base.dividend_yield_q,
                base.volatility_sigma,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = ChesneyJeanblancYorParisianBarrierSpec {
                spot_price_s0: dual_in[0],
                strike_k: dual_in[1],
                down_barrier_h: dual_in[2],
                excursion_window_years_d: dual_in[3],
                risk_free_rate_r: dual_in[4],
                dividend_yield_q: dual_in[5],
                volatility_sigma: dual_in[6],
                maturity_years_t: dual_in[7],
            };
            let r_ad = evaluate_chesney_jeanblanc_yor_parisian_excursion_barrier_call(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| ChesneyJeanblancYorParisianBarrierSpec {
                    spot_price_s0: v[0],
                    strike_k: v[1],
                    down_barrier_h: v[2],
                    excursion_window_years_d: v[3],
                    risk_free_rate_r: v[4],
                    dividend_yield_q: v[5],
                    volatility_sigma: v[6],
                    maturity_years_t: v[7],
                };
                let rp = evaluate_chesney_jeanblanc_yor_parisian_excursion_barrier_call(&mk(p));
                let rm = evaluate_chesney_jeanblanc_yor_parisian_excursion_barrier_call(&mk(m));
                let fd = (rp.parisian_down_and_out_call_value
                    - rm.parisian_down_and_out_call_value)
                    / (2.0 * eps);
                assert!(
                    r_ad.parisian_down_and_out_call_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.parisian_down_and_out_call_value.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use chesney_jeanblanc_yor_parisian_mod::*;

mod broadie_glasserman_kou_discrete_barrier_mod {
    use ql_ad::Real;

    /// Input specification for the Broadie, Glasserman & Kou (1997, *Mathematical Finance* 7(4),
    /// *"A Continuity Correction for Discrete Barrier Options"*) discrete down-and-out call
    /// continuity correction ($\beta = -\zeta(1/2)/\sqrt{2\pi} \approx 0.5826$) engine.
    #[derive(Debug, Clone, Copy)]
    pub struct BroadieGlassermanKouDiscreteBarrierCallSpec<T: Real> {
        /// Underlying spot price $S_0 > H$ (e.g., $100.0$).
        pub spot_price_s0: T,
        /// Option strike price $K \ge H$ (e.g., $100.0$).
        pub strike_price_k: T,
        /// Contractual lower knock-out barrier $H < S_0$ (e.g., $88.0$).
        pub lower_barrier_h: T,
        /// Risk-free continuously compounded rate $r > 0$ (e.g., $0.045$).
        pub risk_free_rate_r: T,
        /// Underlying return volatility $\sigma > 0$ (e.g., $0.24$).
        pub volatility_sigma: T,
        /// Option maturity $T > 0$ in years (e.g., $1.0$).
        pub maturity_years_t: T,
        /// Number of discrete monitoring dates $M > 0$ over $[0, T]$ (e.g., $52.0$).
        pub monitoring_intervals_m: T,
    }

    /// Analytics returned by [`evaluate_broadie_glasserman_kou_discrete_barrier_call`].
    #[derive(Debug, Clone, Copy)]
    pub struct BroadieGlassermanKouDiscreteBarrierCallAnalytics<T: Real> {
        /// Continuity-corrected discrete barrier $H_{\text{adj}} = H \exp(-\beta_{\text{BGK}}\sigma\sqrt{T/M}) < H$.
        pub shifted_discrete_barrier_h_adj: T,
        /// Continuous-monitoring down-and-out call price $C_{\text{DO}}^{\text{cont}}(H)$.
        pub continuous_down_and_out_call: T,
        /// Discrete-monitoring Broadie–Glasserman–Kou down-and-out call price $C_{\text{DO}}^{\text{disc}}(H_{\text{adj}})$.
        pub discrete_bgk_down_and_out_call: T,
        /// Discrete vs. continuous monitoring price pickup $C_{\text{DO}}^{\text{disc}} - C_{\text{DO}}^{\text{cont}} > 0$.
        pub discrete_monitoring_price_pickup: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        let z = x * T::from_f64(0.707_106_781_186_547_5);
        T::from_f64(0.5) * (T::one() + z.erf())
    }

    #[inline]
    fn bs_call<T: Real>(s: T, k: T, r: T, sig: T, t: T) -> T {
        let std_t = sig * t.sqrt();
        let d1 = ((s / k).ln() + (r + T::from_f64(0.5) * sig * sig) * t) / std_t;
        let d2 = d1 - std_t;
        s * std_norm_cdf(d1) - k * (-r * t).exp() * std_norm_cdf(d2)
    }

    /// Evaluates Theorem 1 of Broadie, Glasserman & Kou (1997): shifting the lower knock-out
    /// barrier $H \mapsto H_{\text{adj}} = H \exp(-\beta_{\text{BGK}}\sigma\sqrt{T/M})$ with
    /// $\beta_{\text{BGK}} = -\zeta(1/2)/\sqrt{2\pi} \approx 0.5825971579$ prices the $M$-date
    /// discretely monitored down-and-out call to $O(M^{-1})$ accuracy.
    pub fn evaluate_broadie_glasserman_kou_discrete_barrier_call<T: Real>(
        spec: &BroadieGlassermanKouDiscreteBarrierCallSpec<T>,
    ) -> BroadieGlassermanKouDiscreteBarrierCallAnalytics<T> {
        let s0 = spec.spot_price_s0;
        let k = spec.strike_price_k;
        let h = spec.lower_barrier_h;
        let r = spec.risk_free_rate_r;
        let sig = spec.volatility_sigma;
        let t = spec.maturity_years_t;
        let m = spec.monitoring_intervals_m;

        let beta_bgk = T::from_f64(0.582_597_157_939_010_7);
        let dt_sqrt = (t / m).sqrt();
        let h_adj = h * (-beta_bgk * sig * dt_sqrt).exp();

        let vanilla = bs_call(s0, k, r, sig, t);
        let refl_exp = T::from_f64(2.0) * r / (sig * sig) - T::one();

        let h_ratio_cont = (h / s0).pow(refl_exp);
        let c_refl_cont = bs_call(h * h / s0, k, r, sig, t);
        let c_do_cont = vanilla - h_ratio_cont * c_refl_cont;

        let h_ratio_disc = (h_adj / s0).pow(refl_exp);
        let c_refl_disc = bs_call(h_adj * h_adj / s0, k, r, sig, t);
        let c_do_disc = vanilla - h_ratio_disc * c_refl_disc;

        BroadieGlassermanKouDiscreteBarrierCallAnalytics {
            shifted_discrete_barrier_h_adj: h_adj,
            continuous_down_and_out_call: c_do_cont,
            discrete_bgk_down_and_out_call: c_do_disc,
            discrete_monitoring_price_pickup: c_do_disc - c_do_cont,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn broadie_glasserman_kou_discrete_barrier_invariants_and_ad_parity() {
            let base = BroadieGlassermanKouDiscreteBarrierCallSpec {
                spot_price_s0: 100.0_f64,
                strike_price_k: 100.0_f64,
                lower_barrier_h: 88.0_f64,
                risk_free_rate_r: 0.045_f64,
                volatility_sigma: 0.24_f64,
                maturity_years_t: 1.0_f64,
                monitoring_intervals_m: 52.0_f64,
            };
            let res = evaluate_broadie_glasserman_kou_discrete_barrier_call(&base);
            assert!(res.shifted_discrete_barrier_h_adj < base.lower_barrier_h);
            assert!(res.continuous_down_and_out_call > 0.0);
            assert!(res.discrete_bgk_down_and_out_call > res.continuous_down_and_out_call);
            assert!(res.discrete_monitoring_price_pickup > 0.0);

            let inputs = [
                base.spot_price_s0,
                base.strike_price_k,
                base.lower_barrier_h,
                base.risk_free_rate_r,
                base.volatility_sigma,
                base.maturity_years_t,
                base.monitoring_intervals_m,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BroadieGlassermanKouDiscreteBarrierCallSpec {
                spot_price_s0: dual_in[0],
                strike_price_k: dual_in[1],
                lower_barrier_h: dual_in[2],
                risk_free_rate_r: dual_in[3],
                volatility_sigma: dual_in[4],
                maturity_years_t: dual_in[5],
                monitoring_intervals_m: dual_in[6],
            };
            let ad = evaluate_broadie_glasserman_kou_discrete_barrier_call(&spec_ad)
                .discrete_bgk_down_and_out_call;
            let eps = 1e-7;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let f_p = evaluate_broadie_glasserman_kou_discrete_barrier_call(
                    &BroadieGlassermanKouDiscreteBarrierCallSpec {
                        spot_price_s0: p[0],
                        strike_price_k: p[1],
                        lower_barrier_h: p[2],
                        risk_free_rate_r: p[3],
                        volatility_sigma: p[4],
                        maturity_years_t: p[5],
                        monitoring_intervals_m: p[6],
                    },
                )
                .discrete_bgk_down_and_out_call;
                let f_m = evaluate_broadie_glasserman_kou_discrete_barrier_call(
                    &BroadieGlassermanKouDiscreteBarrierCallSpec {
                        spot_price_s0: m[0],
                        strike_price_k: m[1],
                        lower_barrier_h: m[2],
                        risk_free_rate_r: m[3],
                        volatility_sigma: m[4],
                        maturity_years_t: m[5],
                        monitoring_intervals_m: m[6],
                    },
                )
                .discrete_bgk_down_and_out_call;
                let fd = (f_p - f_m) / (2.0 * eps);
                assert!(ad.d[i].abs() > 1e-6, "zero derivative at index {i}");
                assert!(
                    (ad.d[i] - fd).abs() < 1e-4,
                    "mismatch at index {i}: ad={} fd={}",
                    ad.d[i],
                    fd
                );
            }
        }
    }
}
pub use broadie_glasserman_kou_discrete_barrier_mod::{
    BroadieGlassermanKouDiscreteBarrierCallAnalytics, BroadieGlassermanKouDiscreteBarrierCallSpec,
    evaluate_broadie_glasserman_kou_discrete_barrier_call,
};
