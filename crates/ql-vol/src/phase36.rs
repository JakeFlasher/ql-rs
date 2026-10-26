#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
use ql_ad::Real;
use ql_math::norm_cdf;

#[inline]
fn smooth_ind<T: Real>(z: T) -> T {
    T::one() / (T::one() + (-z).exp())
}

#[derive(Clone, Debug)]
pub struct GatheralBergomiVarianceSwapCorridorSpec<T: Real = f64> {
    pub forward_price: T,
    pub risk_free_rate: T,
    pub maturity_years: T,
    pub atm_vol: T,
    pub skew_slope_beta: T,
    pub smile_curvature_gamma: T,
    pub strike_grid: Vec<T>,
    pub corridor_lower_barrier: T,
    pub corridor_upper_barrier: T,
}

#[derive(Clone, Copy, Debug)]
pub struct GatheralBergomiVarianceSwapCorridorReport<T: Real = f64> {
    pub replicated_variance_strike: T,
    pub replicated_volatility_strike: T,
    pub atm_variance_baseline: T,
    pub smile_convexity_vol_uplift_bps: T,
    pub corridor_variance_strike: T,
    pub corridor_volatility_strike: T,
    pub corridor_variance_share: T,
    pub log_contract_static_hedge_cost: T,
}

pub fn evaluate_gatheral_bergomi_variance_swap_and_corridor<T: Real>(
    spec: &GatheralBergomiVarianceSwapCorridorSpec<T>,
) -> GatheralBergomiVarianceSwapCorridorReport<T> {
    let f0 = spec.forward_price.max_val(T::from_f64(1e-6));
    let mat = spec.maturity_years.max_val(T::from_f64(1e-6));
    let df = (-spec.risk_free_rate * mat).exp();
    let sqrt_t = mat.sqrt();
    let m = spec.strike_grid.len();
    let mut sum_full = T::zero();
    let mut sum_corr = T::zero();

    for j in 0..m {
        let kj = spec.strike_grid[j].max_val(T::from_f64(1e-6));
        let dk = if m <= 1 {
            kj * T::from_f64(0.1)
        } else if j == 0 {
            spec.strike_grid[1] - spec.strike_grid[0]
        } else if j + 1 == m {
            spec.strike_grid[m - 1] - spec.strike_grid[m - 2]
        } else {
            T::from_f64(0.5) * (spec.strike_grid[j + 1] - spec.strike_grid[j - 1])
        };
        let k_log = (kj / f0).ln();
        let sig_k = (spec.atm_vol
            + spec.skew_slope_beta * k_log
            + spec.smile_curvature_gamma * k_log * k_log)
            .max_val(T::from_f64(0.02));
        let st = sig_k * sqrt_t;
        let d1 = (-k_log + T::from_f64(0.5) * st * st) / st;
        let d2 = d1 - st;
        let call = df * (f0 * norm_cdf(d1) - kj * norm_cdf(d2));
        let put = df * (kj * norm_cdf(-d2) - f0 * norm_cdf(-d1));
        let w_call = smooth_ind(T::from_f64(40.0) * (kj - f0) / f0);
        let q_otm = (T::one() - w_call) * put + w_call * call;
        let term = (dk / (kj * kj)) * q_otm;
        sum_full += term;

        let ind_lo = smooth_ind(T::from_f64(20.0) * (kj - spec.corridor_lower_barrier) / f0);
        let ind_hi = smooth_ind(T::from_f64(20.0) * (spec.corridor_upper_barrier - kj) / f0);
        sum_corr += ind_lo * ind_hi * term;
    }

    let scale = T::from_f64(2.0) * (spec.risk_free_rate * mat).exp() / mat;
    let k_var = (scale * sum_full).max_val(T::from_f64(1e-10));
    let sig_var = k_var.sqrt();
    let atm_var = spec.atm_vol * spec.atm_vol;
    let uplift_bps = (sig_var - spec.atm_vol) * T::from_f64(1e4);
    let k_var_corr = (scale * sum_corr)
        .max_val(T::from_f64(1e-12))
        .min_val(k_var);
    let sig_var_corr = k_var_corr.sqrt();
    let share = k_var_corr / k_var;

    GatheralBergomiVarianceSwapCorridorReport {
        replicated_variance_strike: k_var,
        replicated_volatility_strike: sig_var,
        atm_variance_baseline: atm_var,
        smile_convexity_vol_uplift_bps: uplift_bps,
        corridor_variance_strike: k_var_corr,
        corridor_volatility_strike: sig_var_corr,
        corridor_variance_share: share,
        log_contract_static_hedge_cost: T::from_f64(2.0) * sum_full,
    }
}

#[cfg(test)]
mod phase36_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn gatheral_bergomi_variance_swap_corridor_invariants_and_ad_parity() {
        let spec = GatheralBergomiVarianceSwapCorridorSpec {
            forward_price: 100.0,
            risk_free_rate: 0.03,
            maturity_years: 1.0,
            atm_vol: 0.20,
            skew_slope_beta: -0.12,
            smile_curvature_gamma: 0.25,
            strike_grid: vec![70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0],
            corridor_lower_barrier: 85.0,
            corridor_upper_barrier: 115.0,
        };
        let rep = evaluate_gatheral_bergomi_variance_swap_and_corridor(&spec);
        assert!(rep.replicated_volatility_strike > 0.0);
        assert!(rep.corridor_variance_strike <= rep.replicated_variance_strike);
        assert!(rep.corridor_variance_share > 0.0 && rep.corridor_variance_share < 1.0);

        let h = 1e-6;
        let r_ad = evaluate_gatheral_bergomi_variance_swap_and_corridor(
            &GatheralBergomiVarianceSwapCorridorSpec {
                forward_price: Dual::<1>::from_value(spec.forward_price),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                maturity_years: Dual::<1>::from_value(spec.maturity_years),
                atm_vol: Dual::<1>::variable(spec.atm_vol, 0),
                skew_slope_beta: Dual::<1>::from_value(spec.skew_slope_beta),
                smile_curvature_gamma: Dual::<1>::from_value(spec.smile_curvature_gamma),
                strike_grid: spec
                    .strike_grid
                    .iter()
                    .copied()
                    .map(Dual::<1>::from_value)
                    .collect(),
                corridor_lower_barrier: Dual::<1>::from_value(spec.corridor_lower_barrier),
                corridor_upper_barrier: Dual::<1>::from_value(spec.corridor_upper_barrier),
            },
        );
        let mut sp = spec.clone();
        sp.atm_vol += h;
        let mut sm = spec;
        sm.atm_vol -= h;
        let _fd = (evaluate_gatheral_bergomi_variance_swap_and_corridor(&sp)
            .replicated_volatility_strike
            - evaluate_gatheral_bergomi_variance_swap_and_corridor(&sm)
                .replicated_volatility_strike)
            / (2.0 * h);
        assert!((r_ad.replicated_volatility_strike.d[0] - _fd).abs() < 1e-5);
    }
}
