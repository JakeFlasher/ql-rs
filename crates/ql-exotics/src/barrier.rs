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
