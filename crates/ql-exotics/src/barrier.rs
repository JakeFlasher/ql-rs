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
