//! Gatheral (2006) & Bergomi–Guyon (2012) Skew-Stickiness Ratio (`SSR`),
//! Smile-Adjusted Minimum-Variance Delta & Tanaka Local-Time Stop-Loss Slippage Engine.

use ql_ad::Real;
use ql_math::{norm_cdf, norm_pdf};

#[derive(Debug, Clone, Copy)]
pub struct GatheralSsrTanakaStopLossSpec<T: Real = f64> {
    pub spot_price_s0: T,
    pub strike_price_k: T,
    pub risk_free_rate_r: T,
    pub atm_implied_vol_sigma0: T,
    pub atm_skew_slope_dsigma_dk: T,
    pub spot_vol_regression_beta: T,
    pub maturity_years_t: T,
    pub rebalancing_steps_n: T,
}

#[derive(Debug, Clone, Copy)]
pub struct GatheralSsrTanakaStopLossReport<T: Real = f64> {
    pub skew_stickiness_ratio_ssr: T,
    pub black_scholes_delta: T,
    pub black_scholes_vega: T,
    pub sticky_strike_delta_adjustment: T,
    pub minimum_variance_smile_delta: T,
    pub tanaka_local_time_expected_cost: T,
    pub discrete_stop_loss_slippage_std: T,
}

pub fn evaluate_gatheral_ssr_tanaka_stop_loss_hedging<T: Real>(
    spec: &GatheralSsrTanakaStopLossSpec<T>,
) -> GatheralSsrTanakaStopLossReport<T> {
    let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
    let k = spec.strike_price_k.max_val(T::from_f64(1e-4));
    let r = spec.risk_free_rate_r;
    let sig = spec.atm_implied_vol_sigma0.max_val(T::from_f64(1e-4));
    let t = spec.maturity_years_t.max_val(T::from_f64(0.02));
    let sqrt_t = t.sqrt();
    let std = sig * sqrt_t;

    let d1 = ((s0 / k).ln() + (r + T::from_f64(0.5) * sig * sig) * t) / std;
    let bs_delta = norm_cdf(d1);
    let bs_vega = s0 * norm_pdf(d1) * sqrt_t;

    // Skew-Stickiness Ratio (Bergomi 2016, Ch. 2): SSR = beta_{ATM} / (d sigma / d k)|_{k=0}:
    let skew_mag = spec
        .atm_skew_slope_dsigma_dk
        .abs()
        .max_val(T::from_f64(1e-4));
    let ssr = (spec.spot_vol_regression_beta / (-skew_mag)).max_val(T::from_f64(0.1));

    // Sticky-strike delta adjustment (SSR = 0) vs. Minimum-Variance delta (SSR = R_skew):
    let dsigma_ds = spec.atm_skew_slope_dsigma_dk / s0;
    let sticky_strike_adj = -bs_vega * dsigma_ds;
    let mv_delta = bs_delta + bs_vega * ssr * dsigma_ds;

    // Tanaka-Meyer local-time expected cost of continuous stop-loss hedging = 0.5 * sigma^2 * K^2 * G_BS:
    let df = (-r * t).exp();
    let d2 = d1 - std;
    let bs_call = s0 * bs_delta - k * df * norm_cdf(d2);
    let intrinsic = (s0 - k * df).max_val(T::zero());
    let tanaka_cost = (bs_call - intrinsic).max_val(T::zero());
    let n_steps = spec.rebalancing_steps_n.max_val(T::one());
    let slippage_std = tanaka_cost / n_steps.pow(T::from_f64(0.25));

    GatheralSsrTanakaStopLossReport {
        skew_stickiness_ratio_ssr: ssr,
        black_scholes_delta: bs_delta,
        black_scholes_vega: bs_vega,
        sticky_strike_delta_adjustment: sticky_strike_adj,
        minimum_variance_smile_delta: mv_delta,
        tanaka_local_time_expected_cost: tanaka_cost,
        discrete_stop_loss_slippage_std: slippage_std,
    }
}

#[cfg(test)]
mod phase41_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn gatheral_ssr_tanaka_stop_loss_invariants_and_ad_parity() {
        let spec = GatheralSsrTanakaStopLossSpec {
            spot_price_s0: 100.0,
            strike_price_k: 100.0,
            risk_free_rate_r: 0.04,
            atm_implied_vol_sigma0: 0.22,
            atm_skew_slope_dsigma_dk: -0.18,
            spot_vol_regression_beta: -0.24,
            maturity_years_t: 0.5,
            rebalancing_steps_n: 252.0,
        };
        let rep = evaluate_gatheral_ssr_tanaka_stop_loss_hedging(&spec);
        assert!(rep.skew_stickiness_ratio_ssr > 1.0);
        assert!(rep.minimum_variance_smile_delta < rep.black_scholes_delta);
        assert!(rep.tanaka_local_time_expected_cost > 0.0);

        let h = 1e-6;
        let r_ad = evaluate_gatheral_ssr_tanaka_stop_loss_hedging(&GatheralSsrTanakaStopLossSpec {
            spot_price_s0: Dual::<1>::variable(spec.spot_price_s0, 0),
            strike_price_k: Dual::<1>::from_value(spec.strike_price_k),
            risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
            atm_implied_vol_sigma0: Dual::<1>::from_value(spec.atm_implied_vol_sigma0),
            atm_skew_slope_dsigma_dk: Dual::<1>::from_value(spec.atm_skew_slope_dsigma_dk),
            spot_vol_regression_beta: Dual::<1>::from_value(spec.spot_vol_regression_beta),
            maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
            rebalancing_steps_n: Dual::<1>::from_value(spec.rebalancing_steps_n),
        });
        let mut sp = spec;
        sp.spot_price_s0 += h;
        let mut sm = spec;
        sm.spot_price_s0 -= h;
        let fd = (evaluate_gatheral_ssr_tanaka_stop_loss_hedging(&sp).minimum_variance_smile_delta
            - evaluate_gatheral_ssr_tanaka_stop_loss_hedging(&sm).minimum_variance_smile_delta)
            / (2.0 * h);
        assert!((r_ad.minimum_variance_smile_delta.d[0] - fd).abs() < 1e-5);
    }
}
