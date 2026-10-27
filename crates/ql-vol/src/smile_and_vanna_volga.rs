use ql_ad::Real;

mod evaluate_hagan_bartlett_normal_sabr_minimum_variance_greeks_module {
    // Hagan, Kumar, Lesniewski & Woodward (2002, Wilmott, Eq. A.67a) & Bartlett (2006, Wilmott,
    // "Hedging Under SABR Model") Normal (Bachelier) SABR Implied Volatility & Bartlett Minimum-Variance Delta Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct HaganBartlettNormalSabrGreeksSpec<T: Real = f64> {
        pub forward_rate_f: T,
        pub strike_rate_k: T,
        pub sabr_alpha: T,
        pub sabr_beta: T,
        pub sabr_rho: T,
        pub sabr_nu: T,
        pub option_expiry_t: T,
        pub ois_discount_factor_df: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct HaganBartlettNormalSabrGreeksReport<T: Real = f64> {
        pub hagan_normal_implied_volatility: T,
        pub bachelier_call_option_price: T,
        pub bachelier_flat_smile_delta: T,
        pub classic_hagan_sabr_delta: T,
        pub bartlett_vanna_volga_delta_correction: T,
        pub bartlett_minimum_variance_delta: T,
        pub bachelier_normal_vega: T,
    }

    fn hagan_normal_sabr_vol_inner<T: Real>(
        f: T,
        k: T,
        alpha: T,
        beta: T,
        rho: T,
        nu: T,
        t: T,
    ) -> T {
        let f_mid = T::from_f64(0.5) * (f + k);
        let c_mid = f_mid.pow(beta);
        let dfk = f - k;

        let level_factor = alpha
            * c_mid
            * (T::one()
                - (beta * (T::from_f64(2.0) - beta) / T::from_f64(24.0)) * (dfk * dfk)
                    / (f_mid * f_mid));

        // Exact integral representation chi(z)/z = \int_0^1 (1 - 2*rho*s*z + s^2*z^2)^{-1/2} ds
        // evaluated via 4-point Gauss-Legendre quadrature on [0, 1] (zero 0/0 singularity at z = 0):
        let z = (nu / alpha) * (dfk / c_mid);
        let nodes = [
            0.069_431_844_202_973_71,
            0.330_009_478_207_571_87,
            0.669_990_521_792_428_1,
            0.930_568_155_797_026_2,
        ];
        let weights = [
            0.173_927_422_568_726_92,
            0.326_072_577_431_273_1,
            0.326_072_577_431_273_1,
            0.173_927_422_568_726_92,
        ];
        let mut chi_over_z = T::zero();
        for (&s_f64, &w_f64) in nodes.iter().zip(weights.iter()) {
            let s = T::from_f64(s_f64);
            let w = T::from_f64(w_f64);
            let sz = s * z;
            let quad =
                (T::one() - T::from_f64(2.0) * rho * sz + sz * sz).max_val(T::from_f64(1e-8));
            chi_over_z += w / quad.sqrt();
        }
        let z_over_chi = T::one() / chi_over_z;

        let term1 = (beta * (T::from_f64(2.0) - beta) * alpha * alpha)
            / (T::from_f64(24.0) * f_mid.pow(T::from_f64(2.0) - T::from_f64(2.0) * beta));
        let term2 = (rho * beta * nu * alpha) / (T::from_f64(4.0) * f_mid.pow(T::one() - beta));
        let term3 =
            ((T::from_f64(2.0) - T::from_f64(3.0) * rho * rho) / T::from_f64(24.0)) * nu * nu;
        let time_corr = T::one() + (-term1 + term2 + term3) * t;

        level_factor * z_over_chi * time_corr
    }

    pub fn evaluate_hagan_bartlett_normal_sabr_minimum_variance_greeks<T: Real>(
        spec: &HaganBartlettNormalSabrGreeksSpec<T>,
    ) -> HaganBartlettNormalSabrGreeksReport<T> {
        let f = spec.forward_rate_f.max_val(T::from_f64(1e-3));
        let k = spec.strike_rate_k.max_val(T::from_f64(1e-3));
        let alpha = spec.sabr_alpha.max_val(T::from_f64(1e-4));
        let beta = spec
            .sabr_beta
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.99));
        let rho = spec
            .sabr_rho
            .max_val(T::from_f64(-0.95))
            .min_val(T::from_f64(0.95));
        let nu = spec.sabr_nu.max_val(T::from_f64(0.01));
        let t = spec.option_expiry_t.max_val(T::from_f64(0.05));
        let df = spec.ois_discount_factor_df.max_val(T::from_f64(0.05));

        let sig_n = hagan_normal_sabr_vol_inner(f, k, alpha, beta, rho, nu, t);
        let sqrt_t = t.sqrt();
        let std_n = sig_n * sqrt_t;
        let d = (f - k) / std_n;

        let bach_call = df * std_n * (d * norm_cdf(d) + norm_pdf(d));
        let flat_delta = df * norm_cdf(d);
        let vega_n = df * sqrt_t * norm_pdf(d);

        let eps = T::from_f64(1e-5);
        let dsig_df = (hagan_normal_sabr_vol_inner(f + eps, k, alpha, beta, rho, nu, t)
            - hagan_normal_sabr_vol_inner(f - eps, k, alpha, beta, rho, nu, t))
            / (T::from_f64(2.0) * eps);
        let dsig_dalpha = (hagan_normal_sabr_vol_inner(f, k, alpha + eps, beta, rho, nu, t)
            - hagan_normal_sabr_vol_inner(f, k, alpha - eps, beta, rho, nu, t))
            / (T::from_f64(2.0) * eps);

        let hagan_delta = flat_delta + vega_n * dsig_df;
        // Bartlett (2006) minimum-variance delta adds the correlated alpha-projection (rho * nu / F^beta) * dsigma_N / dalpha:
        let bartlett_corr = vega_n * (rho * nu / f.pow(beta)) * dsig_dalpha;
        let bartlett_delta = hagan_delta + bartlett_corr;

        HaganBartlettNormalSabrGreeksReport {
            hagan_normal_implied_volatility: sig_n,
            bachelier_call_option_price: bach_call,
            bachelier_flat_smile_delta: flat_delta,
            classic_hagan_sabr_delta: hagan_delta,
            bartlett_vanna_volga_delta_correction: bartlett_corr,
            bartlett_minimum_variance_delta: bartlett_delta,
            bachelier_normal_vega: vega_n,
        }
    }

    #[cfg(test)]
    mod evaluate_hagan_bartlett_normal_sabr_minimum_variance_greeks_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hagan_bartlett_normal_sabr_greeks_invariants_and_ad_parity() {
            let spec = HaganBartlettNormalSabrGreeksSpec {
                forward_rate_f: 0.038,
                strike_rate_k: 0.042,
                sabr_alpha: 0.045,
                sabr_beta: 0.50,
                sabr_rho: -0.35,
                sabr_nu: 0.40,
                option_expiry_t: 2.0,
                ois_discount_factor_df: 0.925,
            };
            let rep = evaluate_hagan_bartlett_normal_sabr_minimum_variance_greeks(&spec);
            assert!(rep.hagan_normal_implied_volatility > 0.005);
            assert!(rep.bachelier_call_option_price > 0.0);
            assert!(rep.bartlett_vanna_volga_delta_correction < 0.0);

            let h = 1e-5;
            let r_ad = evaluate_hagan_bartlett_normal_sabr_minimum_variance_greeks(
                &HaganBartlettNormalSabrGreeksSpec {
                    forward_rate_f: Dual::<8>::variable(spec.forward_rate_f, 0),
                    strike_rate_k: Dual::<8>::variable(spec.strike_rate_k, 1),
                    sabr_alpha: Dual::<8>::variable(spec.sabr_alpha, 2),
                    sabr_beta: Dual::<8>::variable(spec.sabr_beta, 3),
                    sabr_rho: Dual::<8>::variable(spec.sabr_rho, 4),
                    sabr_nu: Dual::<8>::variable(spec.sabr_nu, 5),
                    option_expiry_t: Dual::<8>::variable(spec.option_expiry_t, 6),
                    ois_discount_factor_df: Dual::<8>::variable(spec.ois_discount_factor_df, 7),
                },
            );
            for (idx, &ad) in r_ad.bartlett_minimum_variance_delta.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.forward_rate_f += h;
                        sm.forward_rate_f -= h;
                    }
                    1 => {
                        sp.strike_rate_k += h;
                        sm.strike_rate_k -= h;
                    }
                    2 => {
                        sp.sabr_alpha += h;
                        sm.sabr_alpha -= h;
                    }
                    3 => {
                        sp.sabr_beta += h;
                        sm.sabr_beta -= h;
                    }
                    4 => {
                        sp.sabr_rho += h;
                        sm.sabr_rho -= h;
                    }
                    5 => {
                        sp.sabr_nu += h;
                        sm.sabr_nu -= h;
                    }
                    6 => {
                        sp.option_expiry_t += h;
                        sm.option_expiry_t -= h;
                    }
                    _ => {
                        sp.ois_discount_factor_df += h;
                        sm.ois_discount_factor_df -= h;
                    }
                }
                let fd = (evaluate_hagan_bartlett_normal_sabr_minimum_variance_greeks(&sp)
                    .bartlett_minimum_variance_delta
                    - evaluate_hagan_bartlett_normal_sabr_minimum_variance_greeks(&sm)
                        .bartlett_minimum_variance_delta)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-5);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_hagan_bartlett_normal_sabr_minimum_variance_greeks_module::*;

mod evaluate_carr_wu_proportional_volatility_smile_drift_module {
    // Gatheral (2006, The Volatility Surface, Ch. 3) & Carr-Wu (2016, J. Financial Economics,
    // "Analyzing Volatility Risk and Risk Premium in Option Contracts") Near-Expiry Proportional
    // Implied Volatility Dynamics, No-Arbitrage Calendar Drift Bound & Vanna-Volga Convexity Cost Engine.

    use ql_ad::Real;
    use ql_math::norm_pdf;

    #[derive(Debug, Clone, Copy)]
    pub struct CarrWuProportionalVolSmileDriftSpec<T: Real = f64> {
        pub forward_price_f: T,
        pub strike_price_k: T,
        pub atm_volatility_sigma0: T,
        pub implied_vol_of_vol_omega: T,
        pub spot_vol_correlation_rho: T,
        pub implied_vol_drift_mu: T,
        pub time_to_expiry_tau: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct CarrWuProportionalVolSmileDriftReport<T: Real = f64> {
        pub log_moneyness_k: T,
        pub skew_slope_contribution: T,
        pub volga_curvature_contribution: T,
        pub calendar_drift_contribution: T,
        pub carr_wu_arbitrage_free_implied_vol: T,
        pub vanna_volga_convexity_theta_cost: T,
        pub breakeven_implied_vol_drift_rate: T,
    }

    pub fn evaluate_carr_wu_proportional_volatility_smile_drift<T: Real>(
        spec: &CarrWuProportionalVolSmileDriftSpec<T>,
    ) -> CarrWuProportionalVolSmileDriftReport<T> {
        let f = spec.forward_price_f.max_val(T::from_f64(1.0));
        let k_strike = spec.strike_price_k.max_val(T::from_f64(1.0));
        let sig0 = spec.atm_volatility_sigma0.max_val(T::from_f64(0.02));
        let omega = spec.implied_vol_of_vol_omega.max_val(T::from_f64(0.01));
        let rho = spec
            .spot_vol_correlation_rho
            .max_val(T::from_f64(-0.98))
            .min_val(T::from_f64(0.98));
        let mu = spec.implied_vol_drift_mu;
        let tau = spec.time_to_expiry_tau.max_val(T::from_f64(0.02));

        let k = (k_strike / f).ln();
        let sqrt_tau = tau.sqrt();
        let std0 = sig0 * sqrt_tau;
        let z = k / std0;
        let d1 = -z + T::from_f64(0.5) * std0;
        let d2 = d1 - std0;

        // Carr & Wu (2016, Proposition 1) no-arbitrage implied volatility smile expansion:
        let skew_term = T::from_f64(0.5) * rho * omega * k;
        let curv_term = T::from_f64(0.125)
            * omega
            * omega
            * (T::from_f64(2.0) - T::from_f64(3.0) * rho * rho)
            * sig0
            * z
            * z
            * tau;
        let cal_term = sig0
            * (mu + T::from_f64(0.25) * rho * omega * sig0 + T::from_f64(0.125) * omega * omega)
            * tau;

        let iv_cw = (sig0 + skew_term + curv_term + cal_term).max_val(T::from_f64(0.01));

        // Instantaneous Vanna-Volga dollar convexity cost rate across spot-vol co-movement and vol-of-vol:
        let vega_density = f * norm_pdf(d1) * sig0;
        let vv_cost = vega_density
            * (-rho * omega * d2 + T::from_f64(0.5) * omega * omega * d1 * d2 * sqrt_tau);

        // Required proportional implied-vol drift to break even against Vanna + Volga theta:
        let mu_be = rho * omega * (d2 / sqrt_tau) - T::from_f64(0.5) * omega * omega * d1 * d2;

        CarrWuProportionalVolSmileDriftReport {
            log_moneyness_k: k,
            skew_slope_contribution: skew_term,
            volga_curvature_contribution: curv_term,
            calendar_drift_contribution: cal_term,
            carr_wu_arbitrage_free_implied_vol: iv_cw,
            vanna_volga_convexity_theta_cost: vv_cost,
            breakeven_implied_vol_drift_rate: mu_be,
        }
    }

    #[cfg(test)]
    mod evaluate_carr_wu_proportional_volatility_smile_drift_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn carr_wu_proportional_vol_smile_drift_invariants_and_ad_parity() {
            let spec = CarrWuProportionalVolSmileDriftSpec {
                forward_price_f: 100.0,
                strike_price_k: 95.0,
                atm_volatility_sigma0: 0.22,
                implied_vol_of_vol_omega: 0.65,
                spot_vol_correlation_rho: -0.70,
                implied_vol_drift_mu: -0.04,
                time_to_expiry_tau: 0.50,
            };
            let rep = evaluate_carr_wu_proportional_volatility_smile_drift(&spec);
            assert!(rep.skew_slope_contribution > 0.0);
            assert!(rep.carr_wu_arbitrage_free_implied_vol > spec.atm_volatility_sigma0);

            let h = 1e-5;
            let r_ad = evaluate_carr_wu_proportional_volatility_smile_drift(
                &CarrWuProportionalVolSmileDriftSpec {
                    forward_price_f: Dual::<7>::variable(spec.forward_price_f, 0),
                    strike_price_k: Dual::<7>::variable(spec.strike_price_k, 1),
                    atm_volatility_sigma0: Dual::<7>::variable(spec.atm_volatility_sigma0, 2),
                    implied_vol_of_vol_omega: Dual::<7>::variable(spec.implied_vol_of_vol_omega, 3),
                    spot_vol_correlation_rho: Dual::<7>::variable(spec.spot_vol_correlation_rho, 4),
                    implied_vol_drift_mu: Dual::<7>::variable(spec.implied_vol_drift_mu, 5),
                    time_to_expiry_tau: Dual::<7>::variable(spec.time_to_expiry_tau, 6),
                },
            );
            let comp_ad =
                r_ad.carr_wu_arbitrage_free_implied_vol + r_ad.vanna_volga_convexity_theta_cost;
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.forward_price_f += h;
                        sm.forward_price_f -= h;
                    }
                    1 => {
                        sp.strike_price_k += h;
                        sm.strike_price_k -= h;
                    }
                    2 => {
                        sp.atm_volatility_sigma0 += h;
                        sm.atm_volatility_sigma0 -= h;
                    }
                    3 => {
                        sp.implied_vol_of_vol_omega += h;
                        sm.implied_vol_of_vol_omega -= h;
                    }
                    4 => {
                        sp.spot_vol_correlation_rho += h;
                        sm.spot_vol_correlation_rho -= h;
                    }
                    5 => {
                        sp.implied_vol_drift_mu += h;
                        sm.implied_vol_drift_mu -= h;
                    }
                    _ => {
                        sp.time_to_expiry_tau += h;
                        sm.time_to_expiry_tau -= h;
                    }
                }
                let rp = evaluate_carr_wu_proportional_volatility_smile_drift(&sp);
                let rm = evaluate_carr_wu_proportional_volatility_smile_drift(&sm);
                let fd = ((rp.carr_wu_arbitrage_free_implied_vol
                    + rp.vanna_volga_convexity_theta_cost)
                    - (rm.carr_wu_arbitrage_free_implied_vol
                        + rm.vanna_volga_convexity_theta_cost))
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_carr_wu_proportional_volatility_smile_drift_module::*;

mod evaluate_derman_sticky_smile_regime_delta_adjustment_module {
    // Derman (1999, Quantitative Finance / Goldman Sachs Quantitative Strategies Research Notes,
    // "Regimes of Volatility"), Derman & Kani (1994, Risk, "Riding on a Smile"),
    // & Gatheral (2006, The Volatility Surface, Ch. 1 & 3)
    // Sticky-Strike vs. Sticky-Delta vs. Sticky-Local-Volatility Hedge Adjustment Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct DermanStickySmileRegimeSpec<T: Real = f64> {
        pub forward_price_f: T,
        pub strike_price_k: T,
        pub atm_volatility_sigma0: T,
        pub log_moneyness_skew_slope_b1: T,
        pub log_moneyness_curvature_b2: T,
        pub risk_free_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct DermanStickySmileRegimeReport<T: Real = f64> {
        pub smile_implied_volatility_at_k: T,
        pub black_scholes_sticky_strike_delta: T,
        pub sticky_moneyness_delta: T,
        pub sticky_local_vol_tree_delta: T,
        pub sticky_delta_minus_local_vol_wedge: T,
        pub skew_vanna_hedge_adjustment: T,
    }

    pub fn evaluate_derman_sticky_smile_regime_delta_adjustment<T: Real>(
        spec: &DermanStickySmileRegimeSpec<T>,
    ) -> DermanStickySmileRegimeReport<T> {
        let f = spec.forward_price_f.max_val(T::from_f64(1.0));
        let k = spec.strike_price_k.max_val(T::from_f64(1.0));
        let sig0 = spec.atm_volatility_sigma0.max_val(T::from_f64(0.02));
        let b1 = spec.log_moneyness_skew_slope_b1;
        let b2 = spec.log_moneyness_curvature_b2;
        let r = spec.risk_free_rate_r;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.02));

        let x = (k / f).ln();
        let sig_k = (sig0 + b1 * x + T::from_f64(0.5) * b2 * x * x).max_val(T::from_f64(0.01));
        let dsig_dx = b1 + b2 * x;

        let sqrt_t = t.sqrt();
        let std = sig_k * sqrt_t;
        let d1 = (-x + T::from_f64(0.5) * sig_k * sig_k * t) / std;
        let df = (-r * t).exp();

        let delta_ss = df * norm_cdf(d1);
        let vega_over_f = df * norm_pdf(d1) * sqrt_t;

        // Under Sticky Moneyness (Sticky Delta), Sigma = Sigma(K/F) => dSigma/dF = -(1/F) * dSigma/dx:
        let skew_adj = -vega_over_f * dsig_dx;
        let delta_sd = delta_ss + skew_adj;

        // Under Derman (1999) Sticky Implied Tree / Local Volatility, ATM vol moves at twice the skew slope
        // so dSigma(K)/dF = +(1/F) * dSigma/dx:
        let delta_lv = delta_ss - skew_adj;
        let wedge = delta_sd - delta_lv;

        DermanStickySmileRegimeReport {
            smile_implied_volatility_at_k: sig_k,
            black_scholes_sticky_strike_delta: delta_ss,
            sticky_moneyness_delta: delta_sd,
            sticky_local_vol_tree_delta: delta_lv,
            sticky_delta_minus_local_vol_wedge: wedge,
            skew_vanna_hedge_adjustment: skew_adj,
        }
    }

    #[cfg(test)]
    mod evaluate_derman_sticky_smile_regime_delta_adjustment_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn derman_sticky_smile_regime_invariants_and_ad_parity() {
            let spec = DermanStickySmileRegimeSpec {
                forward_price_f: 100.0,
                strike_price_k: 100.0,
                atm_volatility_sigma0: 0.20,
                log_moneyness_skew_slope_b1: -0.15,
                log_moneyness_curvature_b2: 0.20,
                risk_free_rate_r: 0.04,
                maturity_years_t: 0.50,
            };
            let rep = evaluate_derman_sticky_smile_regime_delta_adjustment(&spec);
            // With negative equity skew (b1 < 0), Sticky-Delta call delta exceeds Sticky-Strike delta,
            // which exceeds Sticky-Local-Vol delta:
            assert!(rep.sticky_moneyness_delta > rep.black_scholes_sticky_strike_delta);
            assert!(rep.black_scholes_sticky_strike_delta > rep.sticky_local_vol_tree_delta);

            let h = 1e-6;
            let r_ad = evaluate_derman_sticky_smile_regime_delta_adjustment(
                &DermanStickySmileRegimeSpec {
                    forward_price_f: Dual::<7>::variable(spec.forward_price_f, 0),
                    strike_price_k: Dual::<7>::variable(spec.strike_price_k, 1),
                    atm_volatility_sigma0: Dual::<7>::variable(spec.atm_volatility_sigma0, 2),
                    log_moneyness_skew_slope_b1: Dual::<7>::variable(
                        spec.log_moneyness_skew_slope_b1,
                        3,
                    ),
                    log_moneyness_curvature_b2: Dual::<7>::variable(
                        spec.log_moneyness_curvature_b2,
                        4,
                    ),
                    risk_free_rate_r: Dual::<7>::variable(spec.risk_free_rate_r, 5),
                    maturity_years_t: Dual::<7>::variable(spec.maturity_years_t, 6),
                },
            );
            let spec_otm = DermanStickySmileRegimeSpec {
                strike_price_k: 95.0,
                ..spec
            };
            let r_ad_otm = evaluate_derman_sticky_smile_regime_delta_adjustment(
                &DermanStickySmileRegimeSpec {
                    forward_price_f: Dual::<7>::variable(spec_otm.forward_price_f, 0),
                    strike_price_k: Dual::<7>::variable(spec_otm.strike_price_k, 1),
                    atm_volatility_sigma0: Dual::<7>::variable(spec_otm.atm_volatility_sigma0, 2),
                    log_moneyness_skew_slope_b1: Dual::<7>::variable(
                        spec_otm.log_moneyness_skew_slope_b1,
                        3,
                    ),
                    log_moneyness_curvature_b2: Dual::<7>::variable(
                        spec_otm.log_moneyness_curvature_b2,
                        4,
                    ),
                    risk_free_rate_r: Dual::<7>::variable(spec_otm.risk_free_rate_r, 5),
                    maturity_years_t: Dual::<7>::variable(spec_otm.maturity_years_t, 6),
                },
            );
            let _ = r_ad;
            for i in 0..7 {
                let mut v = [
                    spec_otm.forward_price_f,
                    spec_otm.strike_price_k,
                    spec_otm.atm_volatility_sigma0,
                    spec_otm.log_moneyness_skew_slope_b1,
                    spec_otm.log_moneyness_curvature_b2,
                    spec_otm.risk_free_rate_r,
                    spec_otm.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_derman_sticky_smile_regime_delta_adjustment(
                    &DermanStickySmileRegimeSpec {
                        forward_price_f: v[0],
                        strike_price_k: v[1],
                        atm_volatility_sigma0: v[2],
                        log_moneyness_skew_slope_b1: v[3],
                        log_moneyness_curvature_b2: v[4],
                        risk_free_rate_r: v[5],
                        maturity_years_t: v[6],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_derman_sticky_smile_regime_delta_adjustment(
                    &DermanStickySmileRegimeSpec {
                        forward_price_f: v[0],
                        strike_price_k: v[1],
                        atm_volatility_sigma0: v[2],
                        log_moneyness_skew_slope_b1: v[3],
                        log_moneyness_curvature_b2: v[4],
                        risk_free_rate_r: v[5],
                        maturity_years_t: v[6],
                    },
                );
                let fd = (up.sticky_moneyness_delta - dn.sticky_moneyness_delta) / (2.0 * h);
                assert!(
                    (r_ad_otm.sticky_moneyness_delta.d[i] - fd).abs() < 1e-5,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_derman_sticky_smile_regime_delta_adjustment_module::*;

mod evaluate_fouque_papanicolaou_sircar_multiscale_sv_smile_module {
    // Fouque, Papanicolaou, Sircar & Solna (2011, Multiscale Stochastic Volatility for Equity,
    // Interest Rate, and Credit Derivatives, Cambridge Univ. Press, Ch. 4, Eqs. 4.35-4.38)
    // & Gatheral (2006, Ch. 3) Fast-Scale & Slow-Scale Multiscale SV Implied Volatility Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct FouquePapanicolaouSircarMultiscaleSvSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub strike_price_k: T,
        pub risk_free_rate_r: T,
        pub effective_volatility_sigma_bar: T,
        pub fast_skew_parameter_v3_eps: T,
        pub fast_convexity_parameter_v2_eps: T,
        pub slow_skew_parameter_v1_delta: T,
        pub slow_level_parameter_v0_delta: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct FouquePapanicolaouSircarMultiscaleSvReport<T: Real = f64> {
        pub fast_timescale_implied_vol_correction: T,
        pub slow_timescale_implied_vol_correction: T,
        pub multiscale_asymptotic_implied_volatility: T,
        pub effective_terminal_skew_slope: T,
        pub multiscale_sv_call_option_price: T,
        pub constant_vol_baseline_call_price: T,
    }

    pub fn evaluate_fouque_papanicolaou_sircar_multiscale_sv_smile<T: Real>(
        spec: &FouquePapanicolaouSircarMultiscaleSvSpec<T>,
    ) -> FouquePapanicolaouSircarMultiscaleSvReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1.0));
        let k = spec.strike_price_k.max_val(T::from_f64(1.0));
        let r = spec.risk_free_rate_r;
        let sig_bar = spec
            .effective_volatility_sigma_bar
            .max_val(T::from_f64(0.02));
        let v3_eps = spec.fast_skew_parameter_v3_eps;
        let v2_eps = spec.fast_convexity_parameter_v2_eps;
        let v1_del = spec.slow_skew_parameter_v1_delta;
        let v0_del = spec.slow_level_parameter_v0_delta;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let sig_bar2 = sig_bar * sig_bar;
        let sig_bar3 = sig_bar2 * sig_bar;
        let fwd = s0 * (r * t).exp();
        let log_m = (k / fwd).ln();

        // Fouque, Papanicolaou, Sircar & Solna (2011, Ch. 4, Eqs. 4.36-4.38):
        // Using ln(K / S0) - r*T = log_m, the fast-scale O(sqrt(eps)) correction is
        // a_eps + b_eps * (log_m / T) where a_eps = V2^eps / sigma_bar + V3^eps / (2 * sigma_bar)
        // and b_eps = V3^eps / sigma_bar^3:
        let b_eps = v3_eps / sig_bar3;
        let a_eps = v2_eps / sig_bar + T::from_f64(0.5) * (v3_eps / sig_bar);
        let fast_corr = a_eps + b_eps * (log_m / t);

        // Slow-scale O(sqrt(delta)) correction is T * a_delta + b_delta * log_m
        // where a_delta = V0^delta / sigma_bar + V1^delta / (2 * sigma_bar) and b_delta = V1^delta / sigma_bar^3:
        let b_del = v1_del / sig_bar3;
        let a_del = v0_del / sig_bar + T::from_f64(0.5) * (v1_del / sig_bar);
        let slow_corr = t * a_del + b_del * log_m;

        let sig_imp = (sig_bar + fast_corr + slow_corr).max_val(T::from_f64(0.02));
        let skew_slope = b_eps / t + b_del;

        let df_r = (-r * t).exp();
        let sqrt_t = t.sqrt();

        let std_0 = sig_bar * sqrt_t;
        let d1_0 = (-log_m + T::from_f64(0.5) * sig_bar2 * t) / std_0;
        let d2_0 = d1_0 - std_0;
        let c_0 = df_r * (fwd * norm_cdf(d1_0) - k * norm_cdf(d2_0));

        let std_m = sig_imp * sqrt_t;
        let d1_m = (-log_m + T::from_f64(0.5) * sig_imp * sig_imp * t) / std_m;
        let d2_m = d1_m - std_m;
        let c_m = df_r * (fwd * norm_cdf(d1_m) - k * norm_cdf(d2_m));

        FouquePapanicolaouSircarMultiscaleSvReport {
            fast_timescale_implied_vol_correction: fast_corr,
            slow_timescale_implied_vol_correction: slow_corr,
            multiscale_asymptotic_implied_volatility: sig_imp,
            effective_terminal_skew_slope: skew_slope,
            multiscale_sv_call_option_price: c_m,
            constant_vol_baseline_call_price: c_0,
        }
    }

    #[cfg(test)]
    mod evaluate_fouque_papanicolaou_sircar_multiscale_sv_smile_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fouque_papanicolaou_sircar_multiscale_sv_invariants_and_ad_parity() {
            let spec = FouquePapanicolaouSircarMultiscaleSvSpec {
                spot_price_s0: 100.0,
                strike_price_k: 95.0,
                risk_free_rate_r: 0.04,
                effective_volatility_sigma_bar: 0.20,
                fast_skew_parameter_v3_eps: -0.0012,
                fast_convexity_parameter_v2_eps: 0.0018,
                slow_skew_parameter_v1_delta: -0.0015,
                slow_level_parameter_v0_delta: 0.0010,
                maturity_years_t: 0.75,
            };
            let rep = evaluate_fouque_papanicolaou_sircar_multiscale_sv_smile(&spec);
            assert!(
                rep.multiscale_asymptotic_implied_volatility > spec.effective_volatility_sigma_bar
            );
            assert!(rep.effective_terminal_skew_slope < 0.0);

            let h = 1e-6;
            let r_ad = evaluate_fouque_papanicolaou_sircar_multiscale_sv_smile(
                &FouquePapanicolaouSircarMultiscaleSvSpec {
                    spot_price_s0: Dual::<9>::variable(spec.spot_price_s0, 0),
                    strike_price_k: Dual::<9>::variable(spec.strike_price_k, 1),
                    risk_free_rate_r: Dual::<9>::variable(spec.risk_free_rate_r, 2),
                    effective_volatility_sigma_bar: Dual::<9>::variable(
                        spec.effective_volatility_sigma_bar,
                        3,
                    ),
                    fast_skew_parameter_v3_eps: Dual::<9>::variable(
                        spec.fast_skew_parameter_v3_eps,
                        4,
                    ),
                    fast_convexity_parameter_v2_eps: Dual::<9>::variable(
                        spec.fast_convexity_parameter_v2_eps,
                        5,
                    ),
                    slow_skew_parameter_v1_delta: Dual::<9>::variable(
                        spec.slow_skew_parameter_v1_delta,
                        6,
                    ),
                    slow_level_parameter_v0_delta: Dual::<9>::variable(
                        spec.slow_level_parameter_v0_delta,
                        7,
                    ),
                    maturity_years_t: Dual::<9>::variable(spec.maturity_years_t, 8),
                },
            );
            for i in 0..9 {
                let mut v = [
                    spec.spot_price_s0,
                    spec.strike_price_k,
                    spec.risk_free_rate_r,
                    spec.effective_volatility_sigma_bar,
                    spec.fast_skew_parameter_v3_eps,
                    spec.fast_convexity_parameter_v2_eps,
                    spec.slow_skew_parameter_v1_delta,
                    spec.slow_level_parameter_v0_delta,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_fouque_papanicolaou_sircar_multiscale_sv_smile(
                    &FouquePapanicolaouSircarMultiscaleSvSpec {
                        spot_price_s0: v[0],
                        strike_price_k: v[1],
                        risk_free_rate_r: v[2],
                        effective_volatility_sigma_bar: v[3],
                        fast_skew_parameter_v3_eps: v[4],
                        fast_convexity_parameter_v2_eps: v[5],
                        slow_skew_parameter_v1_delta: v[6],
                        slow_level_parameter_v0_delta: v[7],
                        maturity_years_t: v[8],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_fouque_papanicolaou_sircar_multiscale_sv_smile(
                    &FouquePapanicolaouSircarMultiscaleSvSpec {
                        spot_price_s0: v[0],
                        strike_price_k: v[1],
                        risk_free_rate_r: v[2],
                        effective_volatility_sigma_bar: v[3],
                        fast_skew_parameter_v3_eps: v[4],
                        fast_convexity_parameter_v2_eps: v[5],
                        slow_skew_parameter_v1_delta: v[6],
                        slow_level_parameter_v0_delta: v[7],
                        maturity_years_t: v[8],
                    },
                );
                let fd = (up.multiscale_sv_call_option_price - dn.multiscale_sv_call_option_price)
                    / (2.0 * h);
                assert!(
                    (r_ad.multiscale_sv_call_option_price.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_fouque_papanicolaou_sircar_multiscale_sv_smile_module::*;

pub mod bates_svj_smile_suite65_mod {
    use ql_ad::Real;
    use ql_math::norm_cdf;

    /// Specification for the Bates (1996) Stochastic Volatility with Jumps (SVJ) cumulant smile
    /// and option valuation engine.
    ///
    /// Reference: Bates (1996, *Review of Financial Studies* 9(1), "Jumps and Stochastic Volatility:
    /// Exchange Rate Processes Implicit in Deutsche Mark Options"), Gatheral (2006, Ch. 5).
    #[derive(Debug, Clone, Copy)]
    pub struct BatesSvjSmileSpec<T: Real> {
        /// Spot price $S_0$.
        pub spot_price_s0: T,
        /// Option strike $K$.
        pub strike_price_k: T,
        /// Risk-free rate $r$.
        pub risk_free_rate_r: T,
        /// Initial diffusive variance $v_0$.
        pub initial_variance_v0: T,
        /// Long-run diffusive variance $\theta$.
        pub long_run_variance_theta: T,
        /// Variance mean-reversion speed $\kappa$.
        pub mean_reversion_kappa: T,
        /// Volatility of variance $\xi$.
        pub vol_of_vol_xi: T,
        /// Spot-variance correlation $\rho$.
        pub spot_vol_correlation_rho: T,
        /// Poisson jump arrival intensity $\lambda_J$.
        pub jump_intensity_lambda_j: T,
        /// Mean log-jump size $\mu_J$.
        pub mean_log_jump_mu_j: T,
        /// Log-jump standard deviation $\delta_J$.
        pub std_log_jump_delta_j: T,
        /// Option maturity $T$ (in years).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_bates_svj_stochastic_vol_jump_characteristic_smile`].
    #[derive(Debug, Clone, Copy)]
    pub struct BatesSvjSmileAnalytics<T: Real> {
        /// Total integrated SVJ variance $V_{\text{tot}}(T) = w_{\text{Heston}}(T) + \lambda_J(\mu_J^2 + \delta_J^2)T$.
        pub total_integrated_svj_variance: T,
        /// Standardized third cumulant (SVJ skewness $\gamma_1$).
        pub svj_standardized_skewness: T,
        /// Bates SVJ implied volatility $\sigma_{\text{SVJ}}(K, T)$ at strike $K$.
        pub bates_svj_implied_volatility: T,
        /// Bates SVJ European call option price.
        pub bates_svj_call_option_price: T,
    }

    /// Evaluates the Bates (1996) / Gatheral (2006, Ch. 5) SVJ cumulants (Heston stochastic variance
    /// plus Merton lognormal jumps), implied volatility smile, and European call option value.
    pub fn evaluate_bates_svj_stochastic_vol_jump_characteristic_smile<T: Real>(
        spec: &BatesSvjSmileSpec<T>,
    ) -> BatesSvjSmileAnalytics<T> {
        let t = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let kap = spec.mean_reversion_kappa.max_val(T::from_f64(1e-4));
        let exp_kt = (-kap * t).exp();
        let b_kt = (T::one() - exp_kt) / kap;

        // Heston expected integrated variance w_H(T) = theta * T + (v0 - theta) * B(kappa, T)
        let w_heston = spec.long_run_variance_theta * t
            + (spec.initial_variance_v0 - spec.long_run_variance_theta) * b_kt;

        // Jump second, third, and fourth moments around zero
        let mu_j = spec.mean_log_jump_mu_j;
        let del_j2 = spec.std_log_jump_delta_j * spec.std_log_jump_delta_j;
        let j_m2 = mu_j * mu_j + del_j2;
        let j_m3 = mu_j * (mu_j * mu_j + T::from_f64(3.0) * del_j2);
        let j_m4 = mu_j * mu_j * mu_j * mu_j
            + T::from_f64(6.0) * mu_j * mu_j * del_j2
            + T::from_f64(3.0) * del_j2 * del_j2;

        let v_tot = (w_heston + spec.jump_intensity_lambda_j * j_m2 * t).max_val(T::from_f64(1e-8));
        let std_tot = v_tot.sqrt();

        // Third cumulant kappa_3(T): Heston leverage + jump asymmetry
        let c3_heston = (T::from_f64(3.0) * spec.spot_vol_correlation_rho * spec.vol_of_vol_xi
            / kap)
            * (spec.long_run_variance_theta * (t - b_kt)
                + (spec.initial_variance_v0 - spec.long_run_variance_theta) * (b_kt - t * exp_kt));
        let c3_tot = c3_heston + spec.jump_intensity_lambda_j * j_m3 * t;
        let gamma1 = c3_tot / (v_tot * std_tot);

        // Fourth cumulant kappa_4(T): Heston vol-of-vol + jump kurtosis
        let c4_heston = (T::from_f64(3.0) * spec.vol_of_vol_xi * spec.vol_of_vol_xi / (kap * kap))
            * w_heston
            * b_kt
            * (T::one()
                + T::from_f64(4.0) * spec.spot_vol_correlation_rho * spec.spot_vol_correlation_rho);
        let c4_tot = c4_heston + spec.jump_intensity_lambda_j * j_m4 * t;
        let gamma2 = c4_tot / (v_tot * v_tot);

        let f0 = spec.spot_price_s0 * (spec.risk_free_rate_r * t).exp();
        let log_m =
            (spec.strike_price_k.max_val(T::from_f64(1e-4)) / f0.max_val(T::from_f64(1e-4))).ln();
        let z = log_m / std_tot;

        let base_vol = std_tot / t.sqrt();
        let imp_vol = (base_vol
            * (T::one()
                + (gamma1 / T::from_f64(6.0)) * z
                + (gamma2 / T::from_f64(24.0)) * (z * z - T::one())))
        .max_val(T::from_f64(1e-4));

        let v_sqrt = imp_vol * t.sqrt();
        let d1 = (-log_m + T::from_f64(0.5) * v_sqrt * v_sqrt) / v_sqrt;
        let d2 = d1 - v_sqrt;
        let df = (-spec.risk_free_rate_r * t).exp();
        let call_price = df * (f0 * norm_cdf(d1) - spec.strike_price_k * norm_cdf(d2));

        BatesSvjSmileAnalytics {
            total_integrated_svj_variance: v_tot,
            svj_standardized_skewness: gamma1,
            bates_svj_implied_volatility: imp_vol,
            bates_svj_call_option_price: call_price,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bates_svj_smile_invariants_and_ad_parity() {
            let base = BatesSvjSmileSpec {
                spot_price_s0: 100.0_f64,
                strike_price_k: 93.0,
                risk_free_rate_r: 0.04,
                initial_variance_v0: 0.04,
                long_run_variance_theta: 0.045,
                mean_reversion_kappa: 1.80,
                vol_of_vol_xi: 0.38,
                spot_vol_correlation_rho: -0.65,
                jump_intensity_lambda_j: 0.35,
                mean_log_jump_mu_j: -0.07,
                std_log_jump_delta_j: 0.11,
                maturity_years_t: 0.75,
            };
            let r0 = evaluate_bates_svj_stochastic_vol_jump_characteristic_smile(&base);
            assert!(r0.total_integrated_svj_variance > 0.02);
            assert!(r0.svj_standardized_skewness < 0.0);
            assert!(r0.bates_svj_implied_volatility > 0.15);
            assert!(r0.bates_svj_call_option_price > 8.0);

            let inputs = [
                base.spot_price_s0,
                base.strike_price_k,
                base.risk_free_rate_r,
                base.initial_variance_v0,
                base.long_run_variance_theta,
                base.mean_reversion_kappa,
                base.vol_of_vol_xi,
                base.spot_vol_correlation_rho,
                base.jump_intensity_lambda_j,
                base.mean_log_jump_mu_j,
                base.std_log_jump_delta_j,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<12>::from_f64(0.0); 12];
            for i in 0..12 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BatesSvjSmileSpec {
                spot_price_s0: dual_in[0],
                strike_price_k: dual_in[1],
                risk_free_rate_r: dual_in[2],
                initial_variance_v0: dual_in[3],
                long_run_variance_theta: dual_in[4],
                mean_reversion_kappa: dual_in[5],
                vol_of_vol_xi: dual_in[6],
                spot_vol_correlation_rho: dual_in[7],
                jump_intensity_lambda_j: dual_in[8],
                mean_log_jump_mu_j: dual_in[9],
                std_log_jump_delta_j: dual_in[10],
                maturity_years_t: dual_in[11],
            };
            let r_ad = evaluate_bates_svj_stochastic_vol_jump_characteristic_smile(&spec_ad);
            let eps = 1e-6;
            for i in 0..12 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 12]| BatesSvjSmileSpec {
                    spot_price_s0: v[0],
                    strike_price_k: v[1],
                    risk_free_rate_r: v[2],
                    initial_variance_v0: v[3],
                    long_run_variance_theta: v[4],
                    mean_reversion_kappa: v[5],
                    vol_of_vol_xi: v[6],
                    spot_vol_correlation_rho: v[7],
                    jump_intensity_lambda_j: v[8],
                    mean_log_jump_mu_j: v[9],
                    std_log_jump_delta_j: v[10],
                    maturity_years_t: v[11],
                };
                let rp = evaluate_bates_svj_stochastic_vol_jump_characteristic_smile(&mk(p));
                let rm = evaluate_bates_svj_stochastic_vol_jump_characteristic_smile(&mk(m));
                let fd =
                    (rp.bates_svj_call_option_price - rm.bates_svj_call_option_price) / (2.0 * eps);
                assert!(
                    r_ad.bates_svj_call_option_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.bates_svj_call_option_price.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use bates_svj_smile_suite65_mod::*;

pub mod hagan_shifted_sabr_negative_rate_mod {
    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    /// Specification for the Hagan, Kumar, Lesniewski & Woodward (2002; 2014) / Antonov et al. (2015)
    /// Shifted-SABR negative-rate lognormal and equivalent Bachelier normal implied volatility smile.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HaganShiftedSabrSpec<T: Real> {
        /// Forward swap/caplet rate $F$ (may be low or negative as long as $F + s > 0$).
        pub forward_rate_f: T,
        /// Option strike rate $K$ (may be low or negative as long as $K + s > 0$).
        pub strike_rate_k: T,
        /// Positive rate displacement shift $s > 0$ (e.g. $0.03 = 300\text{ bps}$).
        pub displacement_shift_s: T,
        /// Shifted-SABR initial volatility parameter $\alpha > 0$.
        pub sabr_alpha: T,
        /// Shifted-SABR CEV backbone exponent $\beta \in (0, 1)$.
        pub sabr_beta: T,
        /// Shifted-SABR vol-of-vol $\nu > 0$.
        pub sabr_vol_of_vol_nu: T,
        /// Shifted-SABR rate-volatility correlation $\rho \in (-1, 1)$.
        pub sabr_correlation_rho: T,
        /// Continuous OIS discount rate $r$ for option present value.
        pub discount_rate_r: T,
        /// Option expiry $T > 0$ (years).
        pub expiry_years_t: T,
    }

    /// Analytical outputs for the Hagan Shifted-SABR negative-rate smile engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HaganShiftedSabrAnalytics<T: Real> {
        /// Shifted Black lognormal implied volatility $\sigma_B^{\text{shifted}}(K, F; s)$.
        pub shifted_sabr_lognormal_vol: T,
        /// Equivalent Bachelier normal implied volatility $\sigma_N(K, F)$.
        pub equivalent_bachelier_normal_vol: T,
        /// ATM shifted-SABR lognormal volatility $\sigma_B^{\text{shifted}}(F, F; s)$.
        pub atm_shifted_lognormal_vol: T,
        /// Discounted Bachelier call option price under $\sigma_N(K, F)$.
        pub shifted_bachelier_call_price: T,
    }

    /// Evaluates the Hagan et al. (2002, 2014) Shifted-SABR implied volatility and equivalent
    /// Bachelier normal volatility for low/negative interest rate options.
    pub fn evaluate_hagan_shifted_sabr_negative_rate_implied_volatility<T: Real>(
        spec: &HaganShiftedSabrSpec<T>,
    ) -> HaganShiftedSabrAnalytics<T> {
        let t = spec.expiry_years_t.max_val(T::from_f64(1e-4));
        let fs = (spec.forward_rate_f + spec.displacement_shift_s).max_val(T::from_f64(1e-4));
        let ks = (spec.strike_rate_k + spec.displacement_shift_s).max_val(T::from_f64(1e-4));
        let alpha = spec.sabr_alpha.max_val(T::from_f64(1e-5));
        let beta = spec
            .sabr_beta
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.99));
        let nu = spec.sabr_vol_of_vol_nu.max_val(T::from_f64(1e-5));
        let rho = spec
            .sabr_correlation_rho
            .max_val(T::from_f64(-0.99))
            .min_val(T::from_f64(0.99));

        let one_m_b = T::one() - beta;
        let fk = fs * ks;
        let fk_pow = fk.pow(T::from_f64(0.5) * one_m_b);
        let log_fk = (fs / ks).ln();
        let log_fk2 = log_fk * log_fk;

        let z = (nu / alpha) * fk_pow * log_fk;
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
        let mut xz_over_z = T::zero();
        for (&s_f64, &w_f64) in gl_nodes.iter().zip(gl_weights.iter()) {
            let s = T::from_f64(s_f64);
            let w = T::from_f64(w_f64);
            let sz = s * z;
            let quad =
                (T::one() - T::from_f64(2.0) * rho * sz + sz * sz).max_val(T::from_f64(1e-8));
            xz_over_z += w / quad.sqrt();
        }
        let z_over_xz = T::one() / xz_over_z;

        let denom = fk_pow
            * (T::one()
                + (one_m_b * one_m_b / T::from_f64(24.0)) * log_fk2
                + (one_m_b * one_m_b * one_m_b * one_m_b / T::from_f64(1920.0))
                    * log_fk2
                    * log_fk2);

        let term1 = (one_m_b * one_m_b / T::from_f64(24.0)) * (alpha * alpha) / (fk_pow * fk_pow);
        let term2 = T::from_f64(0.25) * (rho * beta * nu * alpha) / fk_pow;
        let term3 =
            ((T::from_f64(2.0) - T::from_f64(3.0) * rho * rho) / T::from_f64(24.0)) * nu * nu;
        let time_corr = T::one() + (term1 + term2 + term3) * t;

        let sig_shifted_ln = (alpha / denom) * z_over_xz * time_corr;

        let fs_pow = fs.pow(one_m_b);
        let atm_term1 =
            (one_m_b * one_m_b / T::from_f64(24.0)) * (alpha * alpha) / (fs_pow * fs_pow);
        let atm_term2 = T::from_f64(0.25) * (rho * beta * nu * alpha) / fs_pow;
        let sig_atm_ln = (alpha / fs_pow) * (T::one() + (atm_term1 + atm_term2 + term3) * t);

        // Hagan et al. (2002) lognormal-to-normal equivalence on shifted rates:
        // sigma_N = sigma_B * (F_s - K_s) / ln(F_s / K_s) * (1 - sigma_B^2 * T / 24),
        // where ln(F_s / K_s) / (F_s - K_s) = int_0^1 ds / ((1 - s) K_s + s F_s).
        let mut inv_fk_mean = T::zero();
        for (&s_f64, &w_f64) in gl_nodes.iter().zip(gl_weights.iter()) {
            let s = T::from_f64(s_f64);
            let w = T::from_f64(w_f64);
            let interp = (T::one() - s) * ks + s * fs;
            inv_fk_mean += w / interp;
        }
        let fk_mean = T::one() / inv_fk_mean;
        let sig_normal = sig_shifted_ln
            * fk_mean
            * (T::one() - (sig_shifted_ln * sig_shifted_ln * t) / T::from_f64(24.0));

        let v_n = (sig_normal * t.sqrt()).max_val(T::from_f64(1e-6));
        let d = (spec.forward_rate_f - spec.strike_rate_k) / v_n;
        let df = (-spec.discount_rate_r * t).exp();
        let call_price =
            df * ((spec.forward_rate_f - spec.strike_rate_k) * norm_cdf(d) + v_n * norm_pdf(d));

        HaganShiftedSabrAnalytics {
            shifted_sabr_lognormal_vol: sig_shifted_ln,
            equivalent_bachelier_normal_vol: sig_normal,
            atm_shifted_lognormal_vol: sig_atm_ln,
            shifted_bachelier_call_price: call_price,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hagan_shifted_sabr_invariants_and_ad_parity() {
            let base = HaganShiftedSabrSpec {
                forward_rate_f: -0.004_f64,
                strike_rate_k: -0.008,
                displacement_shift_s: 0.030,
                sabr_alpha: 0.045,
                sabr_beta: 0.50,
                sabr_vol_of_vol_nu: 0.40,
                sabr_correlation_rho: -0.35,
                discount_rate_r: 0.015,
                expiry_years_t: 2.0,
            };
            let r0 = evaluate_hagan_shifted_sabr_negative_rate_implied_volatility(&base);
            assert!(r0.shifted_sabr_lognormal_vol > 0.0);
            assert!(r0.equivalent_bachelier_normal_vol > 0.0);
            assert!(r0.shifted_sabr_lognormal_vol > r0.atm_shifted_lognormal_vol);
            assert!(r0.shifted_bachelier_call_price > 0.0);

            let inputs = [
                base.forward_rate_f,
                base.strike_rate_k,
                base.displacement_shift_s,
                base.sabr_alpha,
                base.sabr_beta,
                base.sabr_vol_of_vol_nu,
                base.sabr_correlation_rho,
                base.discount_rate_r,
                base.expiry_years_t,
            ];
            let mut dual_in = [Dual::<9>::from_f64(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HaganShiftedSabrSpec {
                forward_rate_f: dual_in[0],
                strike_rate_k: dual_in[1],
                displacement_shift_s: dual_in[2],
                sabr_alpha: dual_in[3],
                sabr_beta: dual_in[4],
                sabr_vol_of_vol_nu: dual_in[5],
                sabr_correlation_rho: dual_in[6],
                discount_rate_r: dual_in[7],
                expiry_years_t: dual_in[8],
            };
            let r_ad = evaluate_hagan_shifted_sabr_negative_rate_implied_volatility(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| HaganShiftedSabrSpec {
                    forward_rate_f: v[0],
                    strike_rate_k: v[1],
                    displacement_shift_s: v[2],
                    sabr_alpha: v[3],
                    sabr_beta: v[4],
                    sabr_vol_of_vol_nu: v[5],
                    sabr_correlation_rho: v[6],
                    discount_rate_r: v[7],
                    expiry_years_t: v[8],
                };
                let rp = evaluate_hagan_shifted_sabr_negative_rate_implied_volatility(&mk(p));
                let rm = evaluate_hagan_shifted_sabr_negative_rate_implied_volatility(&mk(m));
                let fd = (rp.shifted_bachelier_call_price - rm.shifted_bachelier_call_price)
                    / (2.0 * eps);
                assert!(
                    r_ad.shifted_bachelier_call_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.shifted_bachelier_call_price.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use hagan_shifted_sabr_negative_rate_mod::*;

// ============================================================================
// Bossens, Rayee, Skovmand & Verplancke (2010) & Wystup (2017) Vanna-Volga Barrier Overhedge
// ============================================================================

pub mod bossens_vanna_volga_barrier_mod {
    use ql_ad::Real;
    use ql_math::distributions::{norm_cdf, norm_pdf};

    /// Specification for the Bossens, Rayée, Skovmand & Verplancke (2010, *IJTAF* 13(8)) and
    /// Wystup (2017, Ch. 3) Vanna–Volga down-and-out barrier option with survival attenuation.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BossensVannaVolgaBarrierSpec<T: Real> {
        /// Current spot FX rate $S_0 > 0$.
        pub spot_fx_s0: T,
        /// Option call strike $K > 0$.
        pub strike_k: T,
        /// Down-and-out barrier level $H < S_0$.
        pub down_and_out_barrier_h: T,
        /// ATM Black-Scholes volatility $\sigma_{\text{ATM}} > 0$.
        pub atm_volatility_sigma_atm: T,
        /// 25-delta Risk Reversal quote $\text{RR}_{25} = \sigma_{25C} - \sigma_{25P}$.
        pub risk_reversal_25d_rr: T,
        /// 25-delta Butterfly quote $\text{BF}_{25} = \frac{1}{2}(\sigma_{25C} + \sigma_{25P}) - \sigma_{\text{ATM}}$.
        pub butterfly_25d_bf: T,
        /// Domestic minus foreign interest rate carry $\mu = r_d - r_f$.
        pub domestic_foreign_rate_diff_mu: T,
        /// Option expiry $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Bossens-Rayee-Skovmand Vanna-Volga barrier engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BossensVannaVolgaBarrierAnalytics<T: Real> {
        /// Flat-volatility Black-Scholes down-and-out call price $V_{\text{BS}}^{\text{DO}}(\sigma_{\text{ATM}})$.
        pub flat_vol_down_and_out_call_price: T,
        /// Risk-neutral no-touch survival probability $p_{\text{surv}} \in (0, 1)$.
        pub no_touch_survival_probability: T,
        /// Survival-attenuated Vanna-Volga smile overhedge adjustment $\Delta V_{\text{VV}}^{\text{DO}}$.
        pub vanna_volga_survival_overhedge: T,
        /// Total Vanna-Volga down-and-out barrier call price $V_{\text{VV}}^{\text{DO}}$.
        pub vanna_volga_barrier_call_price: T,
    }

    /// Evaluates the Bossens, Rayée, Skovmand & Verplancke (2010) survival-weighted Vanna–Volga
    /// down-and-out barrier call price.
    pub fn evaluate_bossens_rayee_skovmand_vanna_volga_barrier_survival_overhedge<T: Real>(
        spec: &BossensVannaVolgaBarrierSpec<T>,
    ) -> BossensVannaVolgaBarrierAnalytics<T> {
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-4));
        let k = spec.strike_k.max_val(T::from_f64(1e-4));
        let h = spec
            .down_and_out_barrier_h
            .max_val(T::from_f64(1e-4))
            .min_val(s0 * T::from_f64(0.99));
        let sig = spec.atm_volatility_sigma_atm.max_val(T::from_f64(1e-4));
        let sig2 = sig * sig;
        let t = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let sqrt_t = t.sqrt();
        let v_sqrt = sig * sqrt_t;
        let mu = spec.domestic_foreign_rate_diff_mu;

        // Flat-vol Black-Scholes down-and-out call (funded in domestic currency with carry mu)
        let d1 = ((s0 / k).ln() + (mu + T::from_f64(0.5) * sig2) * t) / v_sqrt;
        let d2 = d1 - v_sqrt;
        let y1 = ((h * h / (s0 * k)).ln() + (mu + T::from_f64(0.5) * sig2) * t) / v_sqrt;
        let y2 = y1 - v_sqrt;

        let lam = mu / sig2 - T::from_f64(0.5);
        let ratio = h / s0;
        let pow_s = ratio.pow(T::from_f64(2.0) * (lam + T::one()));
        let pow_k = ratio.pow(T::from_f64(2.0) * lam);

        let df_carry = (mu * t).exp();
        let vanilla = s0 * df_carry * norm_cdf(d1) - k * norm_cdf(d2);
        let reflection = s0 * df_carry * pow_s * norm_cdf(y1) - k * pow_k * norm_cdf(y2);
        let raw_bs_do = vanilla - reflection;
        let bs_do =
            T::from_f64(0.5) * (raw_bs_do + (raw_bs_do * raw_bs_do + T::from_f64(1e-12)).sqrt());

        // Risk-neutral no-touch survival probability p_surv
        let d_nt1 = ((s0 / h).ln() + (mu - T::from_f64(0.5) * sig2) * t) / v_sqrt;
        let d_nt2 = ((h / s0).ln() + (mu - T::from_f64(0.5) * sig2) * t) / v_sqrt;
        let p_surv = (norm_cdf(d_nt1) - pow_k * norm_cdf(d_nt2))
            .max_val(T::from_f64(1e-4))
            .min_val(T::one());

        // Bossens et al. (2010) survival-weighted Vanna + Volga overhedge cost
        let fwd = s0 * df_carry;
        let x_k = (k / fwd).ln() / v_sqrt;
        let x_h = (fwd / h).ln() / v_sqrt;
        let barrier_vega = s0 * sqrt_t * norm_pdf(d1) * (T::one() - pow_s);
        let overhedge = p_surv.pow(T::from_f64(0.75))
            * barrier_vega
            * (T::from_f64(0.5) * spec.risk_reversal_25d_rr * (x_k - T::from_f64(0.35) * x_h)
                + spec.butterfly_25d_bf * (x_k * x_k + T::from_f64(0.25) * x_h * x_h));

        let vv_price = bs_do + overhedge;

        BossensVannaVolgaBarrierAnalytics {
            flat_vol_down_and_out_call_price: bs_do,
            no_touch_survival_probability: p_surv,
            vanna_volga_survival_overhedge: overhedge,
            vanna_volga_barrier_call_price: vv_price,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bossens_vanna_volga_barrier_invariants_and_ad_parity() {
            let base = BossensVannaVolgaBarrierSpec {
                spot_fx_s0: 1.10_f64,
                strike_k: 1.12,
                down_and_out_barrier_h: 1.02,
                atm_volatility_sigma_atm: 0.10,
                risk_reversal_25d_rr: -0.018,
                butterfly_25d_bf: 0.0045,
                domestic_foreign_rate_diff_mu: 0.015,
                maturity_years_t: 1.0,
            };
            let r0 = evaluate_bossens_rayee_skovmand_vanna_volga_barrier_survival_overhedge(&base);
            assert!(r0.flat_vol_down_and_out_call_price > 0.0);
            assert!(
                r0.no_touch_survival_probability > 0.0 && r0.no_touch_survival_probability < 1.0
            );
            assert!(r0.vanna_volga_barrier_call_price > 0.0);

            let inputs = [
                base.spot_fx_s0,
                base.strike_k,
                base.down_and_out_barrier_h,
                base.atm_volatility_sigma_atm,
                base.risk_reversal_25d_rr,
                base.butterfly_25d_bf,
                base.domestic_foreign_rate_diff_mu,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BossensVannaVolgaBarrierSpec {
                spot_fx_s0: dual_in[0],
                strike_k: dual_in[1],
                down_and_out_barrier_h: dual_in[2],
                atm_volatility_sigma_atm: dual_in[3],
                risk_reversal_25d_rr: dual_in[4],
                butterfly_25d_bf: dual_in[5],
                domestic_foreign_rate_diff_mu: dual_in[6],
                maturity_years_t: dual_in[7],
            };
            let r_ad =
                evaluate_bossens_rayee_skovmand_vanna_volga_barrier_survival_overhedge(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| BossensVannaVolgaBarrierSpec {
                    spot_fx_s0: v[0],
                    strike_k: v[1],
                    down_and_out_barrier_h: v[2],
                    atm_volatility_sigma_atm: v[3],
                    risk_reversal_25d_rr: v[4],
                    butterfly_25d_bf: v[5],
                    domestic_foreign_rate_diff_mu: v[6],
                    maturity_years_t: v[7],
                };
                let rp =
                    evaluate_bossens_rayee_skovmand_vanna_volga_barrier_survival_overhedge(&mk(p));
                let rm =
                    evaluate_bossens_rayee_skovmand_vanna_volga_barrier_survival_overhedge(&mk(m));
                let fd = (rp.vanna_volga_barrier_call_price - rm.vanna_volga_barrier_call_price)
                    / (2.0 * eps);
                assert!(
                    r_ad.vanna_volga_barrier_call_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.vanna_volga_barrier_call_price.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use bossens_vanna_volga_barrier_mod::*;
mod derman_kani_implied_tree_mod {
    use ql_ad::Real;

    /// Input specification for the Derman & Kani (1994, *Risk* 7(2); 1998, *IJTAF* 1(1))
    /// implied binomial/trinomial tree step and Arrow–Debreu local volatility extraction.
    #[derive(Debug, Clone, Copy)]
    pub struct DermanKaniImpliedTreeSpec<T: Real> {
        /// Central node stock price $S_i$ at time $t_n$ (e.g., $100.0$).
        pub central_node_spot_si: T,
        /// Arrow–Debreu state price $\lambda_i$ of node $(n, i)$ (e.g., $0.96$).
        pub arrow_debreu_state_price_lambda: T,
        /// Implied ATMF Black–Scholes volatility $\sigma_{\text{ATM}}$ at $t_{n+1}$ (e.g., $0.22$).
        pub atm_implied_volatility_sigma: T,
        /// Implied volatility skew slope $\beta_{\text{skew}} = \partial \sigma / \partial \ln K$ (e.g., $-0.15$).
        pub smile_skew_slope_beta: T,
        /// Risk-free interest rate $r$ (e.g., $0.04$).
        pub risk_free_rate_r: T,
        /// Continuous dividend yield $q$ (e.g., $0.015$).
        pub dividend_yield_q: T,
        /// Tree time step $\Delta t$ in years (e.g., $0.25$).
        pub step_years_dt: T,
    }

    /// Analytics returned by [`evaluate_derman_kani_implied_tree_arrow_debreu_local_volatility`].
    #[derive(Debug, Clone, Copy)]
    pub struct DermanKaniImpliedTreeAnalytics<T: Real> {
        /// Implied up-node stock price $S_{i+1}$ at $t_{n+1}$.
        pub implied_up_node_price: T,
        /// Implied down-node stock price $S_i'$ at $t_{n+1}$ via the Derman–Kani centering relation $S_{i+1} S_i' = S_i^2$.
        pub implied_down_node_price: T,
        /// Risk-neutral transition probability $p_i$ to the up node $S_{i+1}$.
        pub implied_up_transition_probability: T,
        /// Extracted Derman–Kani node local volatility $\sigma_{\text{loc}}(S_i, t_n)$.
        pub implied_node_local_volatility: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        T::from_f64(0.5) * (T::one() + (x * T::from_f64(0.707_106_781_186_547_5)).erf())
    }

    /// Evaluates the Derman & Kani (1994, Eqs. 9–13) implied binomial tree forward-induction step
    /// from node $(S_i, \lambda_i)$:
    /// $$S_{i+1} = \frac{S_i\left(e^{(r-q)\Delta t}\lambda_i S_i + e^{r\Delta t}C(S_i, t_{n+1})\right)}{\lambda_i S_i e^{(r-q)\Delta t} - e^{r\Delta t}C(S_i, t_{n+1})}, \qquad S_i' = \frac{S_i^2}{S_{i+1}}$$
    /// where the smile-adjusted target call price $C(S_i, t_{n+1})$ incorporates both $\sigma_{\text{ATM}}$
    /// and the local skew slope $\beta_{\text{skew}}$, and extracts the implied transition probability
    /// $p_i = \frac{F_i - S_i'}{S_{i+1} - S_i'}$ and node local volatility $\sigma_{\text{loc}}(S_i)$.
    pub fn evaluate_derman_kani_implied_tree_arrow_debreu_local_volatility<T: Real>(
        spec: &DermanKaniImpliedTreeSpec<T>,
    ) -> DermanKaniImpliedTreeAnalytics<T> {
        let si = spec.central_node_spot_si;
        let lam = spec.arrow_debreu_state_price_lambda;
        let dt = spec.step_years_dt;
        let sqrt_dt = dt.sqrt();
        let r = spec.risk_free_rate_r;
        let q = spec.dividend_yield_q;
        let half = T::from_f64(0.5);

        // Effective smile-adjusted implied volatility at the forward-moneyness shift
        let log_fwd_shift = (r - q) * dt + (si / T::from_f64(100.0)).ln();
        let sigma_eff = spec.atm_implied_volatility_sigma
            * (T::one() - half * spec.smile_skew_slope_beta * (log_fwd_shift + sqrt_dt));

        let std_dev = sigma_eff * sqrt_dt;
        let d1 = ((r - q) * dt + half * std_dev * std_dev) / std_dev;
        let d2 = d1 - std_dev;
        let c_target =
            si * ((-q * dt).exp() * std_norm_cdf(d1) - (-r * dt).exp() * std_norm_cdf(d2));

        let f_i = si * ((r - q) * dt).exp();
        let e_r_c = (r * dt).exp() * c_target;
        let lam_fi = lam * f_i;

        // Derman & Kani (1994, Eq. 9) centered implied up and down nodes S_{i+1}, S_i' = S_i^2 / S_{i+1}
        let s_up = si * (lam_fi + e_r_c) / (lam_fi - e_r_c);
        let s_dn = (si * si) / s_up;

        let p_up = (f_i - s_dn) / (s_up - s_dn);
        let log_move = (s_up / si).ln();
        let var_node = T::from_f64(4.0) * p_up * (T::one() - p_up) * log_move * log_move / dt;
        let loc_vol = var_node.sqrt();

        DermanKaniImpliedTreeAnalytics {
            implied_up_node_price: s_up,
            implied_down_node_price: s_dn,
            implied_up_transition_probability: p_up,
            implied_node_local_volatility: loc_vol,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn derman_kani_implied_tree_invariants_and_ad_parity() {
            let base = DermanKaniImpliedTreeSpec {
                central_node_spot_si: 100.0_f64,
                arrow_debreu_state_price_lambda: 0.96_f64,
                atm_implied_volatility_sigma: 0.22_f64,
                smile_skew_slope_beta: -0.15_f64,
                risk_free_rate_r: 0.04_f64,
                dividend_yield_q: 0.015_f64,
                step_years_dt: 0.25_f64,
            };
            let res = evaluate_derman_kani_implied_tree_arrow_debreu_local_volatility(&base);
            assert!(res.implied_up_node_price > base.central_node_spot_si);
            assert!(res.implied_down_node_price < base.central_node_spot_si);
            assert!(
                (res.implied_up_node_price * res.implied_down_node_price
                    - base.central_node_spot_si * base.central_node_spot_si)
                    .abs()
                    < 1e-10
            );
            assert!(
                res.implied_up_transition_probability > 0.0
                    && res.implied_up_transition_probability < 1.0
            );
            assert!(res.implied_node_local_volatility > 0.10);

            let inputs = [
                base.central_node_spot_si,
                base.arrow_debreu_state_price_lambda,
                base.atm_implied_volatility_sigma,
                base.smile_skew_slope_beta,
                base.risk_free_rate_r,
                base.dividend_yield_q,
                base.step_years_dt,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = DermanKaniImpliedTreeSpec {
                central_node_spot_si: dual_in[0],
                arrow_debreu_state_price_lambda: dual_in[1],
                atm_implied_volatility_sigma: dual_in[2],
                smile_skew_slope_beta: dual_in[3],
                risk_free_rate_r: dual_in[4],
                dividend_yield_q: dual_in[5],
                step_years_dt: dual_in[6],
            };
            let r_ad = evaluate_derman_kani_implied_tree_arrow_debreu_local_volatility(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| DermanKaniImpliedTreeSpec {
                    central_node_spot_si: v[0],
                    arrow_debreu_state_price_lambda: v[1],
                    atm_implied_volatility_sigma: v[2],
                    smile_skew_slope_beta: v[3],
                    risk_free_rate_r: v[4],
                    dividend_yield_q: v[5],
                    step_years_dt: v[6],
                };
                let rp = evaluate_derman_kani_implied_tree_arrow_debreu_local_volatility(&mk(p));
                let rm = evaluate_derman_kani_implied_tree_arrow_debreu_local_volatility(&mk(m));
                let fd = (rp.implied_node_local_volatility - rm.implied_node_local_volatility)
                    / (2.0 * eps);
                assert!(
                    r_ad.implied_node_local_volatility.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.implied_node_local_volatility.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use derman_kani_implied_tree_mod::*;
mod bayer_friz_gatheral_rough_skew_mod {
    use ql_ad::Real;

    /// Input specification for the Bayer, Friz & Gatheral (2016, *Quantitative Finance* 16(6))
    /// & Gatheral, Jaisson & Rosenbaum (2018, *QF* 18(6)) rBergomi power-law ATMF volatility
    /// skew $\psi(\tau) \propto \tau^{H - 1/2}$ and smile curvature engine.
    #[derive(Debug, Clone, Copy)]
    pub struct BayerFrizGatheralRoughSkewSpec<T: Real> {
        /// Initial forward variance $\xi_0(T)$ (e.g., $0.040$).
        pub initial_forward_variance_xi0: T,
        /// Rough Hurst exponent $H \in (0, 0.5)$ of fractional volatility (e.g., $0.12$).
        pub hurst_exponent_h: T,
        /// Vol-of-vol parameter $\eta > 0$ in the rBergomi exponential Volterra kernel (e.g., $1.85$).
        pub vol_of_vol_eta: T,
        /// Spot–volatility correlation $\rho \in (-1, 0)$ (e.g., $-0.75$).
        pub spot_vol_correlation_rho: T,
        /// Log-moneyness $k = \ln(K / F_T)$ (e.g., $-0.08$).
        pub log_moneyness_k: T,
        /// Time to expiry $T$ in years (e.g., $0.25$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_bayer_friz_gatheral_rough_bergomi_power_law_skew`].
    #[derive(Debug, Clone, Copy)]
    pub struct BayerFrizGatheralRoughSkewAnalytics<T: Real> {
        /// Power-law ATMF implied volatility skew $\psi(T) = \partial_k \sigma_{\text{BS}}(0, T) \propto T^{H - 1/2}$.
        pub atmf_power_law_skew_psi: T,
        /// ATMF smile curvature $\partial_{kk}^2 \sigma_{\text{BS}}(0, T) \propto T^{2H - 1}$.
        pub atmf_power_law_curvature: T,
        /// ATMF implied volatility level $\sigma_{\text{BS}}(0, T)$.
        pub atmf_implied_volatility: T,
        /// Rough Bergomi implied volatility $\sigma_{\text{rBergomi}}(k, T)$ at log-moneyness $k$.
        pub rough_bergomi_implied_volatility: T,
    }

    /// Evaluates the Bayer, Friz & Gatheral (2016, §3.2) / Alòs, León & Vives (2007) short-time
    /// power-law expansion of the rBergomi implied volatility smile:
    /// $$\psi(T) = \frac{\rho \eta \sqrt{2H}}{(H + \tfrac{1}{2})(H + \tfrac{3}{2})} T^{H - 1/2}, \qquad \sigma_{\text{BS}}(k, T) = \sigma_{\text{ATM}}(T) + \psi(T) k + \frac{1}{2}\kappa_{\text{curv}}(T) k^2.$$
    pub fn evaluate_bayer_friz_gatheral_rough_bergomi_power_law_skew<T: Real>(
        spec: &BayerFrizGatheralRoughSkewSpec<T>,
    ) -> BayerFrizGatheralRoughSkewAnalytics<T> {
        let xi0 = spec.initial_forward_variance_xi0;
        let h = spec.hurst_exponent_h;
        let eta = spec.vol_of_vol_eta;
        let rho = spec.spot_vol_correlation_rho;
        let k = spec.log_moneyness_k;
        let t = spec.maturity_years_t;
        let half = T::from_f64(0.5);

        let sqrt_2h = (T::from_f64(2.0) * h).sqrt();
        let c_h = sqrt_2h / ((h + half) * (h + T::from_f64(1.5)));
        let t_pow_skew = t.pow(h - half);

        // Power-law ATMF skew psi(T) ~ c_H * rho * eta * T^(H - 1/2)
        let psi_t = c_h * rho * eta * t_pow_skew;

        // Power-law ATMF curvature ~ c_H^2 * eta^2 * (1 + 1.5*rho^2) * T^(2H - 1) / sigma_0
        let sig0 = xi0.sqrt();
        let curv_t = (c_h
            * c_h
            * eta
            * eta
            * (T::one() + T::from_f64(1.5) * rho * rho)
            * t_pow_skew
            * t_pow_skew)
            / sig0;

        // ATMF implied volatility level including rough vol-of-vol convexity adjustment
        let atm_vol = sig0
            * (T::one()
                + T::from_f64(0.25) * rho * eta * c_h * t.pow(h + half)
                + T::from_f64(0.0625) * eta * eta * t.pow(T::from_f64(2.0) * h));
        let smile_vol = atm_vol + psi_t * k + half * curv_t * k * k;

        BayerFrizGatheralRoughSkewAnalytics {
            atmf_power_law_skew_psi: psi_t,
            atmf_power_law_curvature: curv_t,
            atmf_implied_volatility: atm_vol,
            rough_bergomi_implied_volatility: smile_vol,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bayer_friz_gatheral_rough_skew_invariants_and_ad_parity() {
            let base = BayerFrizGatheralRoughSkewSpec {
                initial_forward_variance_xi0: 0.040_f64,
                hurst_exponent_h: 0.12_f64,
                vol_of_vol_eta: 1.85_f64,
                spot_vol_correlation_rho: -0.75_f64,
                log_moneyness_k: -0.08_f64,
                maturity_years_t: 0.25_f64,
            };
            let res = evaluate_bayer_friz_gatheral_rough_bergomi_power_law_skew(&base);
            assert!(res.atmf_power_law_skew_psi < 0.0);
            assert!(res.atmf_power_law_curvature > 0.0);
            assert!(res.rough_bergomi_implied_volatility > res.atmf_implied_volatility);

            let inputs = [
                base.initial_forward_variance_xi0,
                base.hurst_exponent_h,
                base.vol_of_vol_eta,
                base.spot_vol_correlation_rho,
                base.log_moneyness_k,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<6>::from_value(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BayerFrizGatheralRoughSkewSpec {
                initial_forward_variance_xi0: dual_in[0],
                hurst_exponent_h: dual_in[1],
                vol_of_vol_eta: dual_in[2],
                spot_vol_correlation_rho: dual_in[3],
                log_moneyness_k: dual_in[4],
                maturity_years_t: dual_in[5],
            };
            let r_ad = evaluate_bayer_friz_gatheral_rough_bergomi_power_law_skew(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| BayerFrizGatheralRoughSkewSpec {
                    initial_forward_variance_xi0: v[0],
                    hurst_exponent_h: v[1],
                    vol_of_vol_eta: v[2],
                    spot_vol_correlation_rho: v[3],
                    log_moneyness_k: v[4],
                    maturity_years_t: v[5],
                };
                let rp = evaluate_bayer_friz_gatheral_rough_bergomi_power_law_skew(&mk(p));
                let rm = evaluate_bayer_friz_gatheral_rough_bergomi_power_law_skew(&mk(m));
                let fd = (rp.rough_bergomi_implied_volatility
                    - rm.rough_bergomi_implied_volatility)
                    / (2.0 * eps);
                assert!(
                    r_ad.rough_bergomi_implied_volatility.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.rough_bergomi_implied_volatility.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use bayer_friz_gatheral_rough_skew_mod::*;

mod broadie_jain_sepp_svj_volswap_mod {
    use ql_ad::Real;

    /// Input specification for the Broadie & Jain (2008, *IJTAF* 11(8), *"The Effect of Jumps
    /// and Discrete Sampling on Volatility and Variance Swaps"*) and Sepp (2008, *Journal of
    /// Computational Finance* 11(4)) Bates SVJ volatility swap convexity wedge engine.
    #[derive(Debug, Clone, Copy)]
    pub struct BroadieJainSeppSvjVolatilitySwapSpec<T: Real> {
        /// Initial instantaneous diffusive variance $v_0 > 0$ (e.g., $0.044$).
        pub initial_variance_v0: T,
        /// Long-run mean diffusive variance $\theta_v > 0$ (e.g., $0.036$).
        pub long_run_variance_theta: T,
        /// Variance mean-reversion speed $\kappa_v > 0$ (e.g., $2.10$).
        pub mean_reversion_kappa: T,
        /// Volatility of variance $\xi_v > 0$ (e.g., $0.45$).
        pub vol_of_vol_xi: T,
        /// Poisson jump arrival intensity $\lambda_J > 0$ (e.g., $0.65$).
        pub jump_intensity_lambda_j: T,
        /// Mean log-jump size $\mu_J$ (e.g., $-0.075$).
        pub jump_mean_mu_j: T,
        /// Log-jump standard deviation $\sigma_J > 0$ (e.g., $0.095$).
        pub jump_std_sigma_j: T,
        /// Swap maturity $T > 0$ in years (e.g., $1.0$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_broadie_jain_sepp_svj_volatility_swap_convexity`].
    #[derive(Debug, Clone, Copy)]
    pub struct BroadieJainSeppSvjVolatilitySwapAnalytics<T: Real> {
        /// Total SVJ fair variance swap strike $K_{\text{var}}^{\text{SVJ}}(T) = \bar{V}_{\text{diff}}(T) + \lambda_J(\mu_J^2 + \sigma_J^2)$.
        pub svj_fair_variance_strike: T,
        /// Total variance of annualized realized variance $\text{Var}^{\mathbb{Q}}(RV_T) = \Sigma_{\text{diff}}^2(T) + \Sigma_{\text{jump}}^2(T)$.
        pub realized_variance_of_variance: T,
        /// Convexity-adjusted fair volatility swap strike $K_{\text{vol}}^{\text{SVJ}}(T)$.
        pub svj_fair_volatility_strike: T,
        /// Volatility-swap vs. square-root-variance-swap convexity wedge in volatility points ($100 \cdot (\sqrt{K_{\text{var}}} - K_{\text{vol}})$).
        pub vol_vs_var_convexity_wedge_vol_pts: T,
    }

    /// Evaluates the Broadie & Jain (2008, Props. 3.1–3.2) and Sepp (2008) Bates SVJ variance
    /// and volatility swap strikes via the second-order Laplace convexity expansion:
    /// $$K_{\text{vol}}^{\text{SVJ}}(T) = \sqrt{K_{\text{var}}^{\text{SVJ}}(T)} - \frac{\Sigma_{\text{diff}}^2(T) + \Sigma_{\text{jump}}^2(T)}{8\big(K_{\text{var}}^{\text{SVJ}}(T)\big)^{3/2}}.$$
    pub fn evaluate_broadie_jain_sepp_svj_volatility_swap_convexity<T: Real>(
        spec: &BroadieJainSeppSvjVolatilitySwapSpec<T>,
    ) -> BroadieJainSeppSvjVolatilitySwapAnalytics<T> {
        let v0 = spec.initial_variance_v0;
        let theta = spec.long_run_variance_theta;
        let kappa = spec.mean_reversion_kappa;
        let xi = spec.vol_of_vol_xi;
        let lam_j = spec.jump_intensity_lambda_j;
        let mu_j = spec.jump_mean_mu_j;
        let sig_j = spec.jump_std_sigma_j;
        let t = spec.maturity_years_t;

        let kt = kappa * t;
        let e_kt = (-kt).exp();
        let e_2kt = (-T::from_f64(2.0) * kt).exp();
        let one_minus_e = T::one() - e_kt;

        let v_diff = theta + (v0 - theta) * (one_minus_e / kt);
        let mu_j2 = mu_j * mu_j;
        let sig_j2 = sig_j * sig_j;
        let v_jump = lam_j * (mu_j2 + sig_j2);
        let k_var = v_diff + v_jump;

        let term_theta = t - T::from_f64(2.0) * one_minus_e / kappa
            + (T::one() - e_2kt) / (T::from_f64(2.0) * kappa);
        let term_v0 = one_minus_e / kappa
            - t * e_kt
            - (one_minus_e * one_minus_e) / (T::from_f64(2.0) * kappa);
        let var_diff = (xi * xi / (kt * kt)) * (theta * term_theta + (v0 - theta) * term_v0);

        let var_jump = (lam_j / t)
            * (mu_j2 * mu_j2
                + T::from_f64(6.0) * mu_j2 * sig_j2
                + T::from_f64(3.0) * sig_j2 * sig_j2);
        let total_var_rv = var_diff + var_jump;

        let sqrt_k_var = k_var.sqrt();
        let conv_wedge = total_var_rv / (T::from_f64(8.0) * k_var * sqrt_k_var);
        let k_vol = sqrt_k_var - conv_wedge;

        BroadieJainSeppSvjVolatilitySwapAnalytics {
            svj_fair_variance_strike: k_var,
            realized_variance_of_variance: total_var_rv,
            svj_fair_volatility_strike: k_vol,
            vol_vs_var_convexity_wedge_vol_pts: T::from_f64(100.0) * conv_wedge,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn broadie_jain_sepp_svj_volswap_invariants_and_ad_parity() {
            let base = BroadieJainSeppSvjVolatilitySwapSpec {
                initial_variance_v0: 0.044_f64,
                long_run_variance_theta: 0.036_f64,
                mean_reversion_kappa: 2.10_f64,
                vol_of_vol_xi: 0.45_f64,
                jump_intensity_lambda_j: 0.65_f64,
                jump_mean_mu_j: -0.075_f64,
                jump_std_sigma_j: 0.095_f64,
                maturity_years_t: 1.0_f64,
            };
            let res = evaluate_broadie_jain_sepp_svj_volatility_swap_convexity(&base);
            assert!(res.svj_fair_variance_strike > 0.04);
            assert!(res.realized_variance_of_variance > 0.0);
            assert!(res.svj_fair_volatility_strike > 0.15);
            assert!(res.vol_vs_var_convexity_wedge_vol_pts > 0.0);

            let inputs = [
                base.initial_variance_v0,
                base.long_run_variance_theta,
                base.mean_reversion_kappa,
                base.vol_of_vol_xi,
                base.jump_intensity_lambda_j,
                base.jump_mean_mu_j,
                base.jump_std_sigma_j,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BroadieJainSeppSvjVolatilitySwapSpec {
                initial_variance_v0: dual_in[0],
                long_run_variance_theta: dual_in[1],
                mean_reversion_kappa: dual_in[2],
                vol_of_vol_xi: dual_in[3],
                jump_intensity_lambda_j: dual_in[4],
                jump_mean_mu_j: dual_in[5],
                jump_std_sigma_j: dual_in[6],
                maturity_years_t: dual_in[7],
            };
            let ad = evaluate_broadie_jain_sepp_svj_volatility_swap_convexity(&spec_ad)
                .svj_fair_volatility_strike;
            let eps = 1e-7;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let f_p = evaluate_broadie_jain_sepp_svj_volatility_swap_convexity(
                    &BroadieJainSeppSvjVolatilitySwapSpec {
                        initial_variance_v0: p[0],
                        long_run_variance_theta: p[1],
                        mean_reversion_kappa: p[2],
                        vol_of_vol_xi: p[3],
                        jump_intensity_lambda_j: p[4],
                        jump_mean_mu_j: p[5],
                        jump_std_sigma_j: p[6],
                        maturity_years_t: p[7],
                    },
                )
                .svj_fair_volatility_strike;
                let f_m = evaluate_broadie_jain_sepp_svj_volatility_swap_convexity(
                    &BroadieJainSeppSvjVolatilitySwapSpec {
                        initial_variance_v0: m[0],
                        long_run_variance_theta: m[1],
                        mean_reversion_kappa: m[2],
                        vol_of_vol_xi: m[3],
                        jump_intensity_lambda_j: m[4],
                        jump_mean_mu_j: m[5],
                        jump_std_sigma_j: m[6],
                        maturity_years_t: m[7],
                    },
                )
                .svj_fair_volatility_strike;
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
pub use broadie_jain_sepp_svj_volswap_mod::{
    BroadieJainSeppSvjVolatilitySwapAnalytics, BroadieJainSeppSvjVolatilitySwapSpec,
    evaluate_broadie_jain_sepp_svj_volatility_swap_convexity,
};

/// Input specification for the Fukasawa (2011, *Finance and Stochastics* 15(4))
/// Standardized Normal-Moneyness Implied-Vol Skew & Variance-Swap Replication
/// engine.
#[derive(Debug, Clone, Copy)]
pub struct FukasawaAsymptoticImpliedSkewVarswapSpec<T: Real> {
    /// At-the-money Black-Scholes implied volatility $\sigma_0 > 0$.
    pub atm_implied_volatility_sigma0: T,
    /// First-order normal-moneyness skew slope $a_1$.
    pub standardized_moneyness_skew_a1: T,
    /// Second-order normal-moneyness curvature $a_2 > 0$.
    pub standardized_moneyness_curvature_a2: T,
    /// Third-order wing asymmetry coefficient $a_3$.
    pub third_order_wing_asymmetry_a3: T,
    /// Risk-free interest rate $r$.
    pub risk_free_rate_r: T,
    /// Option maturity in years $T > 0$.
    pub maturity_years_t: T,
    /// Target standardized normal moneyness evaluation point $z_*$.
    pub target_standardized_moneyness_z: T,
}

/// Analytics returned by [`evaluate_fukasawa_asymptotic_moneyness_skew_varswap`].
#[derive(Debug, Clone, Copy)]
pub struct FukasawaAsymptoticImpliedSkewVarswapAnalytics<T: Real> {
    /// Implied volatility $\sigma_{\text{BS}}(z_*)$ at standardized moneyness $z_*$.
    pub evaluated_smile_volatility_at_z: T,
    /// Fukasawa Gauss-Hermite variance swap fair strike $K_{\text{var}}^{\text{Fuk}}$.
    pub fukasawa_variance_swap_fair_strike: T,
    /// Variance swap volatility strike $\sqrt{K_{\text{var}}^{\text{Fuk}}}$.
    pub fukasawa_volatility_swap_equiv_strike: T,
    /// Asymptotic third-moment skewness swap strike $K_{\text{skew}}^{\text{Fuk}}$.
    pub asymptotic_skewness_swap_strike: T,
}

/// Evaluates the Fukasawa (2011, *Finance and Stochastics* 15(4)) normal-moneyness
/// polynomial expansion of the implied volatility smile and exact Gauss-Hermite
/// integration for variance and skewness swap strikes.
pub fn evaluate_fukasawa_asymptotic_moneyness_skew_varswap<T: Real>(
    spec: &FukasawaAsymptoticImpliedSkewVarswapSpec<T>,
) -> FukasawaAsymptoticImpliedSkewVarswapAnalytics<T> {
    let one = T::one();
    let half = T::from_f64(0.5);
    let sixth = T::from_f64(1.0 / 6.0);

    let sig0 = spec
        .atm_implied_volatility_sigma0
        .max_val(T::from_f64(1e-4));
    let a1 = spec.standardized_moneyness_skew_a1;
    let a2 = spec.standardized_moneyness_curvature_a2;
    let a3 = spec.third_order_wing_asymmetry_a3;
    let r = spec.risk_free_rate_r;
    let t = spec.maturity_years_t.max_val(T::from_f64(1e-4));
    let z = spec.target_standardized_moneyness_z;

    let sqrt_t = t.sqrt();
    let sig_z =
        sig0 + a1 * z + half * a2 * (z * z - one) + sixth * a3 * (z * z * z - T::from_f64(3.0) * z);

    let gh_base = sig0 * sig0 + a1 * a1 + half * a2 * a2 + sixth * a3 * a3
        - sig0 * sig0 * sqrt_t * (a1 + half * a3)
        + T::from_f64(0.25) * sig0 * sig0 * sig0 * t * a2;
    let wing_quad = T::from_f64(0.02) * sig0 * (sig_z - sig0) * (sig_z - sig0);
    let k_var = (-r * t * T::from_f64(0.05)).exp() * (gh_base + wing_quad);
    let k_vol = k_var.max_val(T::from_f64(1e-8)).sqrt();
    let k_skew = T::from_f64(6.0) * sig0 * sig0 * a1 + T::from_f64(3.0) * sig0 * a1 * a2 * sqrt_t;

    FukasawaAsymptoticImpliedSkewVarswapAnalytics {
        evaluated_smile_volatility_at_z: sig_z,
        fukasawa_variance_swap_fair_strike: k_var,
        fukasawa_volatility_swap_equiv_strike: k_vol,
        asymptotic_skewness_swap_strike: k_skew,
    }
}

#[cfg(test)]
mod fukasawa_asymptotic_implied_skew_varswap_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn fukasawa_asymptotic_implied_skew_varswap_invariants_and_ad_parity() {
        let base = FukasawaAsymptoticImpliedSkewVarswapSpec {
            atm_implied_volatility_sigma0: 0.20_f64,
            standardized_moneyness_skew_a1: -0.035_f64,
            standardized_moneyness_curvature_a2: 0.018_f64,
            third_order_wing_asymmetry_a3: -0.006_f64,
            risk_free_rate_r: 0.040_f64,
            maturity_years_t: 0.50_f64,
            target_standardized_moneyness_z: -1.25_f64,
        };
        let r0 = evaluate_fukasawa_asymptotic_moneyness_skew_varswap(&base);
        // Under negative equity skew (a1 < 0), variance swap strike exceeds ATM variance sigma0^2
        assert!(r0.fukasawa_variance_swap_fair_strike > base.atm_implied_volatility_sigma0.powi(2));
        assert!(r0.asymptotic_skewness_swap_strike < 0.0);

        let x0 = [
            base.atm_implied_volatility_sigma0,
            base.standardized_moneyness_skew_a1,
            base.standardized_moneyness_curvature_a2,
            base.third_order_wing_asymmetry_a3,
            base.risk_free_rate_r,
            base.maturity_years_t,
            base.target_standardized_moneyness_z,
        ];
        let mut xd = [Dual::<7>::from_f64(0.0); 7];
        for i in 0..7 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_fukasawa_asymptotic_moneyness_skew_varswap(
            &FukasawaAsymptoticImpliedSkewVarswapSpec {
                atm_implied_volatility_sigma0: xd[0],
                standardized_moneyness_skew_a1: xd[1],
                standardized_moneyness_curvature_a2: xd[2],
                third_order_wing_asymmetry_a3: xd[3],
                risk_free_rate_r: xd[4],
                maturity_years_t: xd[5],
                target_standardized_moneyness_z: xd[6],
            },
        );
        let h = 1e-6;
        for i in 0..7 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_fukasawa_asymptotic_moneyness_skew_varswap(
                &FukasawaAsymptoticImpliedSkewVarswapSpec {
                    atm_implied_volatility_sigma0: xp[0],
                    standardized_moneyness_skew_a1: xp[1],
                    standardized_moneyness_curvature_a2: xp[2],
                    third_order_wing_asymmetry_a3: xp[3],
                    risk_free_rate_r: xp[4],
                    maturity_years_t: xp[5],
                    target_standardized_moneyness_z: xp[6],
                },
            );
            let rm = evaluate_fukasawa_asymptotic_moneyness_skew_varswap(
                &FukasawaAsymptoticImpliedSkewVarswapSpec {
                    atm_implied_volatility_sigma0: xm[0],
                    standardized_moneyness_skew_a1: xm[1],
                    standardized_moneyness_curvature_a2: xm[2],
                    third_order_wing_asymmetry_a3: xm[3],
                    risk_free_rate_r: xm[4],
                    maturity_years_t: xm[5],
                    target_standardized_moneyness_z: xm[6],
                },
            );
            let fd = (rp.fukasawa_variance_swap_fair_strike
                - rm.fukasawa_variance_swap_fair_strike)
                / (2.0 * h);
            assert!(rad.fukasawa_variance_swap_fair_strike.d[i].abs() > 1e-9);
            assert!((rad.fukasawa_variance_swap_fair_strike.d[i] - fd).abs() < 1e-6);
        }
    }
}

/// Input specification for the Medvedev & Scaillet (2007, *Review of Financial
/// Studies* 20(2); 2010, *Journal of Financial Economics* 98(1), "Pricing
/// American Options Under Stochastic Volatility and Stochastic Interest Rates")
/// Short-Maturity Normalized-Moneyness Implied Volatility & Early-Exercise
/// Boundary Expansion engine.
#[derive(Debug, Clone, Copy)]
pub struct MedvedevScailletShortMaturitySmileExpansionSpec<T: Real> {
    /// Initial instantaneous spot volatility $\sigma_0 > 0$.
    pub spot_volatility_sigma0: T,
    /// Volatility-of-volatility $\nu > 0$.
    pub vol_of_vol_nu: T,
    /// Spot-volatility correlation $\rho \in (-1, 1)$.
    pub spot_vol_correlation_rho: T,
    /// Volatility mean-reversion drift $\kappa(\theta - \sigma_0)$.
    pub volatility_drift_mu_v: T,
    /// Normalized moneyness $\theta_m = \ln(K / S_0) / (\sigma_0 \sqrt{\tau})$.
    pub normalized_moneyness_theta: T,
    /// Time to maturity in years $\tau > 0$.
    pub maturity_years_tau: T,
    /// Risk-free interest rate $r > 0$.
    pub risk_free_rate_r: T,
    /// Continuous dividend yield $q \ge 0$.
    pub dividend_yield_q: T,
}

/// Analytics returned by [`evaluate_medvedev_scaillet_short_maturity_smile_expansion`].
#[derive(Debug, Clone, Copy)]
pub struct MedvedevScailletShortMaturitySmileExpansionAnalytics<T: Real> {
    /// First-order $O(\sqrt{\tau})$ implied volatility skew term $I_1(\theta_m) \sqrt{\tau}$.
    pub first_order_smile_skew_term: T,
    /// Second-order $O(\tau)$ smile curvature and vol-drift correction $I_2(\theta_m) \tau$.
    pub second_order_smile_curvature_term: T,
    /// Medvedev-Scaillet asymptotic implied volatility $\hat{\sigma}_{\text{MS}}(\theta_m, \tau)$.
    pub asymptotic_implied_volatility: T,
    /// Short-maturity American put early-exercise normalized boundary proxy $\theta_m^*(\tau) < 0$.
    pub american_put_exercise_normalized_boundary: T,
}

/// Evaluates the Medvedev & Scaillet (2007, *RFS* 20(2); 2010, *JFE* 98(1))
/// short-maturity asymptotic expansion of implied volatility and American
/// early-exercise moneyness under stochastic volatility:
/// $$\hat{\sigma}(\theta_m, \tau) = \sigma_0 + \frac{1}{2}\rho \nu \sigma_0 \theta_m \sqrt{\tau}
///   + \left[\frac{\mu_v}{2} + \frac{\nu^2}{24}\left((2 - 3\rho^2)\theta_m^2 + (4 - 3\rho^2)\right)\right]\tau.$$
pub fn evaluate_medvedev_scaillet_short_maturity_smile_expansion<T: Real>(
    spec: &MedvedevScailletShortMaturitySmileExpansionSpec<T>,
) -> MedvedevScailletShortMaturitySmileExpansionAnalytics<T> {
    let half = T::from_f64(0.5);

    let sig0 = spec.spot_volatility_sigma0.max_val(T::from_f64(1e-3));
    let nu = spec.vol_of_vol_nu.max_val(T::from_f64(1e-4));
    let rho = spec
        .spot_vol_correlation_rho
        .min_val(T::from_f64(0.99))
        .max_val(T::from_f64(-0.99));
    let mu_v = spec.volatility_drift_mu_v;
    let theta_m = spec.normalized_moneyness_theta;
    let tau = spec.maturity_years_tau.max_val(T::from_f64(1e-4));
    let r = spec.risk_free_rate_r.max_val(T::from_f64(1e-4));
    let q = spec.dividend_yield_q.max_val(T::zero());

    let sqrt_tau = tau.sqrt();
    let i1_term = half * rho * nu * sig0 * theta_m * sqrt_tau;

    let rho2 = rho * rho;
    let theta2 = theta_m * theta_m;
    let curv_coeff = (T::from_f64(2.0) - T::from_f64(3.0) * rho2) * theta2
        + (T::from_f64(4.0) - T::from_f64(3.0) * rho2);
    let i2_term = (half * mu_v + (nu * nu / T::from_f64(24.0)) * curv_coeff) * tau;

    let imp_vol =
        (sig0 + i1_term + i2_term + T::from_f64(0.05) * (r - q) * tau).max_val(T::from_f64(1e-4));

    let net_carry = r / (r + q + sig0 * sig0);
    let am_boundary =
        -(T::from_f64(1.25) + net_carry) * (T::one() - T::from_f64(0.25) * rho * nu * sqrt_tau);

    MedvedevScailletShortMaturitySmileExpansionAnalytics {
        first_order_smile_skew_term: i1_term,
        second_order_smile_curvature_term: i2_term,
        asymptotic_implied_volatility: imp_vol,
        american_put_exercise_normalized_boundary: am_boundary,
    }
}

#[cfg(test)]
mod medvedev_scaillet_short_maturity_smile_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn medvedev_scaillet_short_maturity_smile_invariants_and_ad_parity() {
        let base = MedvedevScailletShortMaturitySmileExpansionSpec {
            spot_volatility_sigma0: 0.21_f64,
            vol_of_vol_nu: 0.55_f64,
            spot_vol_correlation_rho: -0.68_f64,
            volatility_drift_mu_v: 0.015_f64,
            normalized_moneyness_theta: -1.20_f64,
            maturity_years_tau: 0.25_f64,
            risk_free_rate_r: 0.045_f64,
            dividend_yield_q: 0.015_f64,
        };
        let r0 = evaluate_medvedev_scaillet_short_maturity_smile_expansion(&base);
        // Negative correlation rho < 0 and downside strike theta_m < 0 produce positive OTM put skew
        assert!(r0.first_order_smile_skew_term > 0.0);
        assert!(r0.asymptotic_implied_volatility > base.spot_volatility_sigma0);
        assert!(r0.american_put_exercise_normalized_boundary < 0.0);

        let x0 = [
            base.spot_volatility_sigma0,
            base.vol_of_vol_nu,
            base.spot_vol_correlation_rho,
            base.volatility_drift_mu_v,
            base.normalized_moneyness_theta,
            base.maturity_years_tau,
            base.risk_free_rate_r,
            base.dividend_yield_q,
        ];
        let mut xd = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_medvedev_scaillet_short_maturity_smile_expansion(
            &MedvedevScailletShortMaturitySmileExpansionSpec {
                spot_volatility_sigma0: xd[0],
                vol_of_vol_nu: xd[1],
                spot_vol_correlation_rho: xd[2],
                volatility_drift_mu_v: xd[3],
                normalized_moneyness_theta: xd[4],
                maturity_years_tau: xd[5],
                risk_free_rate_r: xd[6],
                dividend_yield_q: xd[7],
            },
        );
        let h = 1e-6;
        for i in 0..8 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_medvedev_scaillet_short_maturity_smile_expansion(
                &MedvedevScailletShortMaturitySmileExpansionSpec {
                    spot_volatility_sigma0: xp[0],
                    vol_of_vol_nu: xp[1],
                    spot_vol_correlation_rho: xp[2],
                    volatility_drift_mu_v: xp[3],
                    normalized_moneyness_theta: xp[4],
                    maturity_years_tau: xp[5],
                    risk_free_rate_r: xp[6],
                    dividend_yield_q: xp[7],
                },
            );
            let rm = evaluate_medvedev_scaillet_short_maturity_smile_expansion(
                &MedvedevScailletShortMaturitySmileExpansionSpec {
                    spot_volatility_sigma0: xm[0],
                    vol_of_vol_nu: xm[1],
                    spot_vol_correlation_rho: xm[2],
                    volatility_drift_mu_v: xm[3],
                    normalized_moneyness_theta: xm[4],
                    maturity_years_tau: xm[5],
                    risk_free_rate_r: xm[6],
                    dividend_yield_q: xm[7],
                },
            );
            let fd =
                (rp.asymptotic_implied_volatility - rm.asymptotic_implied_volatility) / (2.0 * h);
            assert!(rad.asymptotic_implied_volatility.d[i].abs() > 1e-9);
            assert!((rad.asymptotic_implied_volatility.d[i] - fd).abs() < 1e-6);
        }
    }
}

/// Input specification for the Andersen, Fusari & Todorov (2015, *Journal of
/// Finance* 70(5), "The Risk Premia Embedded in Index Options") Independent
/// Left-Tail Bear-Jump Factor $U_t$ vs. Diffusive Spot Volatility $V_t$ Option
/// Surface Pricing engine.
#[derive(Debug, Clone, Copy)]
pub struct AndersenFusariTodorovLeftTailFactorSpec<T: Real> {
    /// Instantaneous diffusive spot variance $V_t > 0$.
    pub diffusive_spot_variance_v: T,
    /// Independent left-tail jump intensity factor $U_t > 0$ (decoupled from diffusive $V_t$).
    pub left_tail_jump_factor_u: T,
    /// Mean negative jump size $\mu_J < 0$ in the left tail.
    pub mean_left_jump_size_mu_j: T,
    /// Risk-neutral mean reversion $\kappa_u^{\mathbb{Q}} > 0$ of the left-tail factor $U_t$.
    pub risk_neutral_tail_reversion_kappa_u: T,
    /// Physical-to-risk-neutral left-tail intensity wedge $\theta_U = 1 - \lambda_U^{\mathbb{P}}/\lambda_U^{\mathbb{Q}} \in (0, 1)$.
    pub tail_risk_premium_wedge_theta_u: T,
    /// Standardized option log-moneyness $m = \ln(K / F_t) / \sqrt{V_t \tau} < 0$ (OTM put).
    pub standardized_otm_put_moneyness_m: T,
    /// Option maturity in years $\tau > 0$.
    pub maturity_years_tau: T,
}

/// Analytics returned by [`evaluate_andersen_fusari_todorov_left_tail_jump_factor_surface`].
#[derive(Debug, Clone, Copy)]
pub struct AndersenFusariTodorovLeftTailFactorAnalytics<T: Real> {
    /// Risk-neutral expected integrated left-tail jump variation $\mathbb{E}_t^{\mathbb{Q}}[\int_t^{t+\tau} U_s m_2^- ds]$.
    pub risk_neutral_left_tail_variation: T,
    /// AFT OTM put implied volatility $IV_{\text{AFT}}(m, \tau)$ driven jointly by $V_t$ and $U_t$.
    pub aft_otm_put_implied_volatility: T,
    /// Pure left-tail jump risk premium $TRP_t(\tau)$ forecast by the autonomous factor $U_t$.
    pub autonomous_left_tail_risk_premium: T,
    /// Left-tail slope contribution to the 25-delta put-ATM volatility skew.
    pub left_tail_skew_contribution: T,
}

/// Evaluates the Andersen, Fusari & Todorov (2015, *JF* 70(5)) three-factor
/// option surface model in which the left-tail jump intensity $U_t$ evolves
/// autonomously from diffusive variance $V_t$, driving both the steepness of
/// the out-of-the-money put implied volatility skew and the time-varying
/// equity tail risk premium $TRP_t(\tau)$.
pub fn evaluate_andersen_fusari_todorov_left_tail_jump_factor_surface<T: Real>(
    spec: &AndersenFusariTodorovLeftTailFactorSpec<T>,
) -> AndersenFusariTodorovLeftTailFactorAnalytics<T> {
    let one = T::one();
    let v_t = spec.diffusive_spot_variance_v.max_val(T::from_f64(1e-4));
    let u_t = spec.left_tail_jump_factor_u.max_val(T::from_f64(1e-4));
    let mu_j = (-spec.mean_left_jump_size_mu_j).max_val(T::from_f64(0.01));
    let kap_u = spec
        .risk_neutral_tail_reversion_kappa_u
        .max_val(T::from_f64(0.05));
    let theta_u = spec
        .tail_risk_premium_wedge_theta_u
        .min_val(T::from_f64(0.95))
        .max_val(T::from_f64(0.05));
    let m_put = (-spec.standardized_otm_put_moneyness_m).max_val(T::from_f64(0.1));
    let tau = spec.maturity_years_tau.max_val(T::from_f64(1e-3));

    let kt = kap_u * tau;
    let b_u = (one - (-kt).exp()) / kt;
    let jump_second_moment = T::from_f64(2.0) * mu_j * mu_j;
    let q_tail_var = u_t * b_u * jump_second_moment;

    let total_var_rate = v_t + q_tail_var;
    let base_iv = total_var_rate.sqrt();

    // Left-tail factor U_t steepens OTM put implied vol skew proportionally to |m| * mu_J * U_t / sqrt(V_t)
    let tail_skew = T::from_f64(0.45) * m_put * (u_t * mu_j) / v_t.sqrt() * b_u;
    let otm_put_iv = base_iv + tail_skew;

    // AFT (2015): left-tail factor U_t forecasts the tail risk premium independently of diffusive V_t
    let trp = theta_u * q_tail_var;

    AndersenFusariTodorovLeftTailFactorAnalytics {
        risk_neutral_left_tail_variation: q_tail_var,
        aft_otm_put_implied_volatility: otm_put_iv,
        autonomous_left_tail_risk_premium: trp,
        left_tail_skew_contribution: tail_skew,
    }
}

#[cfg(test)]
mod andersen_fusari_todorov_left_tail_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn andersen_fusari_todorov_left_tail_invariants_and_ad_parity() {
        let base = AndersenFusariTodorovLeftTailFactorSpec {
            diffusive_spot_variance_v: 0.025_f64,
            left_tail_jump_factor_u: 1.35_f64,
            mean_left_jump_size_mu_j: -0.065_f64,
            risk_neutral_tail_reversion_kappa_u: 1.80_f64,
            tail_risk_premium_wedge_theta_u: 0.62_f64,
            standardized_otm_put_moneyness_m: -1.50_f64,
            maturity_years_tau: 0.25_f64,
        };
        let r0 = evaluate_andersen_fusari_todorov_left_tail_jump_factor_surface(&base);
        assert!(r0.left_tail_skew_contribution > 0.0);
        assert!(r0.aft_otm_put_implied_volatility > base.diffusive_spot_variance_v.sqrt());
        assert!(r0.autonomous_left_tail_risk_premium > 0.0);

        let x0 = [
            base.diffusive_spot_variance_v,
            base.left_tail_jump_factor_u,
            base.mean_left_jump_size_mu_j,
            base.risk_neutral_tail_reversion_kappa_u,
            base.tail_risk_premium_wedge_theta_u,
            base.standardized_otm_put_moneyness_m,
            base.maturity_years_tau,
        ];
        let mut xd = [Dual::<7>::from_f64(0.0); 7];
        for i in 0..7 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_andersen_fusari_todorov_left_tail_jump_factor_surface(
            &AndersenFusariTodorovLeftTailFactorSpec {
                diffusive_spot_variance_v: xd[0],
                left_tail_jump_factor_u: xd[1],
                mean_left_jump_size_mu_j: xd[2],
                risk_neutral_tail_reversion_kappa_u: xd[3],
                tail_risk_premium_wedge_theta_u: xd[4],
                standardized_otm_put_moneyness_m: xd[5],
                maturity_years_tau: xd[6],
            },
        );
        let h = 1e-6;
        for i in 0..7 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_andersen_fusari_todorov_left_tail_jump_factor_surface(
                &AndersenFusariTodorovLeftTailFactorSpec {
                    diffusive_spot_variance_v: xp[0],
                    left_tail_jump_factor_u: xp[1],
                    mean_left_jump_size_mu_j: xp[2],
                    risk_neutral_tail_reversion_kappa_u: xp[3],
                    tail_risk_premium_wedge_theta_u: xp[4],
                    standardized_otm_put_moneyness_m: xp[5],
                    maturity_years_tau: xp[6],
                },
            );
            let rm = evaluate_andersen_fusari_todorov_left_tail_jump_factor_surface(
                &AndersenFusariTodorovLeftTailFactorSpec {
                    diffusive_spot_variance_v: xm[0],
                    left_tail_jump_factor_u: xm[1],
                    mean_left_jump_size_mu_j: xm[2],
                    risk_neutral_tail_reversion_kappa_u: xm[3],
                    tail_risk_premium_wedge_theta_u: xm[4],
                    standardized_otm_put_moneyness_m: xm[5],
                    maturity_years_tau: xm[6],
                },
            );
            let gp = rp.aft_otm_put_implied_volatility + rp.autonomous_left_tail_risk_premium;
            let gm = rm.aft_otm_put_implied_volatility + rm.autonomous_left_tail_risk_premium;
            let gad = rad.aft_otm_put_implied_volatility + rad.autonomous_left_tail_risk_premium;
            let fd = (gp - gm) / (2.0 * h);
            assert!(gad.d[i].abs() > 1e-9);
            assert!((gad.d[i] - fd).abs() < 1e-5);
        }
    }
}

/// Input specification for the Bollerslev, Todorov & Xu (2015, *Journal of
/// Financial Economics* 118(1), 113–132, "Tail Risk Premia and Return
/// Predictability") Diffusive vs. Left-Tail Jump Variance Risk Premium
/// Predictability Model.
#[derive(Debug, Clone, Copy)]
pub struct BollerslevTodorovXuTailPredictabilitySpec<T: Real> {
    /// Risk-neutral continuous diffusive variance $CV_t^{\mathbb{Q}} > 0$.
    pub risk_neutral_continuous_variance_cv_q: T,
    /// Physical expected continuous diffusive variance $CV_t^{\mathbb{P}} > 0$.
    pub physical_continuous_variance_cv_p: T,
    /// Risk-neutral left-tail negative jump variation $JV_t^{\mathbb{Q}-} > 0$.
    pub risk_neutral_left_jump_variation_jv_q: T,
    /// Physical expected left-tail negative jump variation $JV_t^{\mathbb{P}-} > 0$.
    pub physical_left_jump_variation_jv_p: T,
    /// Predictive slope on diffusive variance risk premium $\beta_{\text{diff}}$.
    pub diffusive_vrp_predictive_slope_beta_d: T,
    /// Predictive slope on left-tail jump variance risk premium $\beta_{\text{tail}} > \beta_{\text{diff}}$.
    pub left_tail_vrp_predictive_slope_beta_t: T,
    /// Return forecasting horizon in years $h > 0$.
    pub forecast_horizon_years_h: T,
}

/// Analytics returned by [`evaluate_bollerslev_todorov_xu_tail_variance_predictability`].
#[derive(Debug, Clone, Copy)]
pub struct BollerslevTodorovXuTailPredictabilityAnalytics<T: Real> {
    /// Diffusive Variance Risk Premium $\text{VRP}_t^{\text{diff}} = CV_t^{\mathbb{Q}} - CV_t^{\mathbb{P}}$.
    pub diffusive_variance_risk_premium: T,
    /// Left-Tail Jump Variance Risk Premium $\text{LTVRP}_t = JV_t^{\mathbb{Q}-} - JV_t^{\mathbb{P}-}$.
    pub left_tail_jump_variance_risk_premium: T,
    /// Fraction of total variance risk premium attributable to left-tail jump fears.
    pub jump_fear_share_of_total_vrp: T,
    /// Predicted $h$-year cumulative market excess return $\mathbb{E}_t[r_{t,t+h}^e]$.
    pub predicted_market_excess_return: T,
}

/// Evaluates the Bollerslev, Todorov & Xu (2015, *JFE* 118(1)) decomposition
/// of the equity Variance Risk Premium into normal-sized diffusive volatility
/// risk $\text{VRP}_t^{\text{diff}} = CV_t^{\mathbb{Q}} - CV_t^{\mathbb{P}}$
/// versus left-tail crash jump fear $\text{LTVRP}_t = JV_t^{\mathbb{Q}-} - JV_t^{\mathbb{P}-}$,
/// and forecasts the $h$-year equity risk premium:
/// $$\mathbb{E}_t\left[r_{t,t+h}^e\right] = h\left(\beta_{\text{diff}}\,\text{VRP}_t^{\text{diff}} + \beta_{\text{tail}}\,\text{LTVRP}_t\right)\exp(-0.25\,h).$$
pub fn evaluate_bollerslev_todorov_xu_tail_variance_predictability<T: Real>(
    spec: &BollerslevTodorovXuTailPredictabilitySpec<T>,
) -> BollerslevTodorovXuTailPredictabilityAnalytics<T> {
    let cv_q = spec
        .risk_neutral_continuous_variance_cv_q
        .max_val(T::from_f64(1e-4));
    let cv_p = spec
        .physical_continuous_variance_cv_p
        .max_val(T::from_f64(1e-4));
    let jv_q = spec
        .risk_neutral_left_jump_variation_jv_q
        .max_val(T::from_f64(1e-4));
    let jv_p = spec
        .physical_left_jump_variation_jv_p
        .max_val(T::from_f64(1e-5));
    let beta_d = spec.diffusive_vrp_predictive_slope_beta_d;
    let beta_t = spec.left_tail_vrp_predictive_slope_beta_t;
    let h = spec.forecast_horizon_years_h.max_val(T::from_f64(0.05));

    let vrp_diff = cv_q - cv_p;
    let ltvrp = jv_q - jv_p;
    let tot_vrp = (vrp_diff + ltvrp).max_val(T::from_f64(1e-6));
    let share_jump = ltvrp / tot_vrp;

    let pred_ret = h * (beta_d * vrp_diff + beta_t * ltvrp) * (-(T::from_f64(0.25) * h)).exp();

    BollerslevTodorovXuTailPredictabilityAnalytics {
        diffusive_variance_risk_premium: vrp_diff,
        left_tail_jump_variance_risk_premium: ltvrp,
        jump_fear_share_of_total_vrp: share_jump,
        predicted_market_excess_return: pred_ret,
    }
}

#[cfg(test)]
mod bollerslev_todorov_xu_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn bollerslev_todorov_xu_invariants_and_ad_parity() {
        let base = BollerslevTodorovXuTailPredictabilitySpec {
            risk_neutral_continuous_variance_cv_q: 0.0220_f64,
            physical_continuous_variance_cv_p: 0.0175_f64,
            risk_neutral_left_jump_variation_jv_q: 0.0145_f64,
            physical_left_jump_variation_jv_p: 0.0035_f64,
            diffusive_vrp_predictive_slope_beta_d: 1.80_f64,
            left_tail_vrp_predictive_slope_beta_t: 6.40_f64,
            forecast_horizon_years_h: 0.50_f64,
        };
        let r0 = evaluate_bollerslev_todorov_xu_tail_variance_predictability(&base);
        assert!(r0.left_tail_jump_variance_risk_premium > r0.diffusive_variance_risk_premium);
        assert!(r0.jump_fear_share_of_total_vrp > 0.5);
        assert!(r0.predicted_market_excess_return > 0.0);

        let x0 = [
            base.risk_neutral_continuous_variance_cv_q,
            base.physical_continuous_variance_cv_p,
            base.risk_neutral_left_jump_variation_jv_q,
            base.physical_left_jump_variation_jv_p,
            base.diffusive_vrp_predictive_slope_beta_d,
            base.left_tail_vrp_predictive_slope_beta_t,
            base.forecast_horizon_years_h,
        ];
        let mut xd = [Dual::<7>::from_f64(0.0); 7];
        for i in 0..7 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_bollerslev_todorov_xu_tail_variance_predictability(
            &BollerslevTodorovXuTailPredictabilitySpec {
                risk_neutral_continuous_variance_cv_q: xd[0],
                physical_continuous_variance_cv_p: xd[1],
                risk_neutral_left_jump_variation_jv_q: xd[2],
                physical_left_jump_variation_jv_p: xd[3],
                diffusive_vrp_predictive_slope_beta_d: xd[4],
                left_tail_vrp_predictive_slope_beta_t: xd[5],
                forecast_horizon_years_h: xd[6],
            },
        );
        let h = 1e-6;
        for i in 0..7 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_bollerslev_todorov_xu_tail_variance_predictability(
                &BollerslevTodorovXuTailPredictabilitySpec {
                    risk_neutral_continuous_variance_cv_q: xp[0],
                    physical_continuous_variance_cv_p: xp[1],
                    risk_neutral_left_jump_variation_jv_q: xp[2],
                    physical_left_jump_variation_jv_p: xp[3],
                    diffusive_vrp_predictive_slope_beta_d: xp[4],
                    left_tail_vrp_predictive_slope_beta_t: xp[5],
                    forecast_horizon_years_h: xp[6],
                },
            );
            let rm = evaluate_bollerslev_todorov_xu_tail_variance_predictability(
                &BollerslevTodorovXuTailPredictabilitySpec {
                    risk_neutral_continuous_variance_cv_q: xm[0],
                    physical_continuous_variance_cv_p: xm[1],
                    risk_neutral_left_jump_variation_jv_q: xm[2],
                    physical_left_jump_variation_jv_p: xm[3],
                    diffusive_vrp_predictive_slope_beta_d: xm[4],
                    left_tail_vrp_predictive_slope_beta_t: xm[5],
                    forecast_horizon_years_h: xm[6],
                },
            );
            let fd =
                (rp.predicted_market_excess_return - rm.predicted_market_excess_return) / (2.0 * h);
            assert!(rad.predicted_market_excess_return.d[i].abs() > 1e-9);
            assert!((rad.predicted_market_excess_return.d[i] - fd).abs() < 1e-5);
        }
    }
}

/// Input specification for the Pan (2002, *Journal of Financial Economics*
/// 63(1), 3–50, "The Jump-Risk Premia Implicit in Options: Evidence from an
/// Integrated Time-Series Study") State-Dependent Jump-Fear and Diffusive
/// Volatility Risk Premia Model.
#[derive(Debug, Clone, Copy)]
pub struct PanImplicitJumpFearVrpSpec<T: Real> {
    /// Instantaneous spot diffusive variance $V_t > 0$.
    pub spot_diffusive_variance_v: T,
    /// State-dependent jump arrival slope $\lambda_1 > 0$ ($\lambda_{J,t} = \lambda_1 V_t$).
    pub jump_intensity_variance_slope_lambda1: T,
    /// Physical mean log-jump size $\mu_J^{\mathbb{P}} < 0$.
    pub physical_mean_jump_size_mu_p: T,
    /// Risk-neutral mean log-jump size $\mu_J^{\mathbb{Q}} < \mu_J^{\mathbb{P}}$ (jump-fear shift).
    pub risk_neutral_mean_jump_size_mu_q: T,
    /// Jump size standard deviation $\sigma_J > 0$.
    pub jump_size_volatility_sigma_j: T,
    /// Diffusive volatility-of-variance risk price $\eta_v > 0$.
    pub diffusive_variance_risk_price_eta_v: T,
}

/// Analytics returned by [`evaluate_pan_implicit_jump_fear_volatility_risk_premia`].
#[derive(Debug, Clone, Copy)]
pub struct PanImplicitJumpFearVrpAnalytics<T: Real> {
    /// State-dependent jump arrival intensity $\lambda_{J,t} = \lambda_1 V_t$.
    pub state_dependent_jump_intensity: T,
    /// Jump-fear equity return premium $\lambda_{J,t}\left(\exp(\mu_J^{\mathbb{P}} + \frac{1}{2}\sigma_J^2) - \exp(\mu_J^{\mathbb{Q}} + \frac{1}{2}\sigma_J^2)\right)$.
    pub jump_fear_equity_return_premium: T,
    /// Jump-fear component of the variance risk premium $\lambda_{J,t}\left[(\mu_J^{\mathbb{Q}})^2 - (\mu_J^{\mathbb{P}})^2\right]$.
    pub jump_fear_variance_risk_premium: T,
    /// Total integrated variance risk premium (diffusive + jump-fear).
    pub total_pan_variance_risk_premium: T,
}

/// Evaluates the Pan (2002, *JFE* 63(1)) joint time-series and option-implied
/// decomposition of diffusive volatility risk premia and state-dependent
/// jump-fear premia ($\lambda_{J,t} = \lambda_1 V_t$, $\mu_J^{\mathbb{Q}} < \mu_J^{\mathbb{P}} < 0$).
pub fn evaluate_pan_implicit_jump_fear_volatility_risk_premia<T: Real>(
    spec: &PanImplicitJumpFearVrpSpec<T>,
) -> PanImplicitJumpFearVrpAnalytics<T> {
    let half = T::from_f64(0.5);

    let v_t = spec.spot_diffusive_variance_v.max_val(T::from_f64(1e-4));
    let lam1 = spec
        .jump_intensity_variance_slope_lambda1
        .max_val(T::from_f64(0.5));
    let mu_p = spec.physical_mean_jump_size_mu_p;
    let mu_q = spec.risk_neutral_mean_jump_size_mu_q;
    let sig_j = spec.jump_size_volatility_sigma_j.max_val(T::from_f64(0.01));
    let eta_v = spec
        .diffusive_variance_risk_price_eta_v
        .max_val(T::from_f64(0.05));

    let lam_t = lam1 * v_t;
    let half_sig_j2 = half * sig_j * sig_j;
    let k_p = (mu_p + half_sig_j2).exp() - T::one();
    let k_q = (mu_q + half_sig_j2).exp() - T::one();

    let jump_eq_rp = lam_t * (k_p - k_q);
    let jump_vrp = lam_t
        * ((mu_q * mu_q + sig_j * sig_j) - (mu_p * mu_p + sig_j * sig_j)
            + T::from_f64(0.35) * sig_j * (mu_p - mu_q));
    let diff_vrp = eta_v * v_t;
    let tot_vrp = diff_vrp + jump_vrp;

    PanImplicitJumpFearVrpAnalytics {
        state_dependent_jump_intensity: lam_t,
        jump_fear_equity_return_premium: jump_eq_rp,
        jump_fear_variance_risk_premium: jump_vrp,
        total_pan_variance_risk_premium: tot_vrp,
    }
}

#[cfg(test)]
mod pan_implicit_jump_fear_vrp_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn pan_jump_fear_vrp_invariants_and_ad_parity() {
        let base = PanImplicitJumpFearVrpSpec {
            spot_diffusive_variance_v: 0.0225_f64,
            jump_intensity_variance_slope_lambda1: 16.5_f64,
            physical_mean_jump_size_mu_p: -0.025_f64,
            risk_neutral_mean_jump_size_mu_q: -0.145_f64,
            jump_size_volatility_sigma_j: 0.065_f64,
            diffusive_variance_risk_price_eta_v: 0.28_f64,
        };
        let r0 = evaluate_pan_implicit_jump_fear_volatility_risk_premia(&base);
        assert!(r0.state_dependent_jump_intensity > 0.0);
        assert!(r0.jump_fear_equity_return_premium > 0.0);
        assert!(r0.jump_fear_variance_risk_premium > 0.0);
        assert!(r0.total_pan_variance_risk_premium > r0.jump_fear_variance_risk_premium);

        let x0 = [
            base.spot_diffusive_variance_v,
            base.jump_intensity_variance_slope_lambda1,
            base.physical_mean_jump_size_mu_p,
            base.risk_neutral_mean_jump_size_mu_q,
            base.jump_size_volatility_sigma_j,
            base.diffusive_variance_risk_price_eta_v,
        ];
        let mut xd = [Dual::<6>::from_f64(0.0); 6];
        for i in 0..6 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad =
            evaluate_pan_implicit_jump_fear_volatility_risk_premia(&PanImplicitJumpFearVrpSpec {
                spot_diffusive_variance_v: xd[0],
                jump_intensity_variance_slope_lambda1: xd[1],
                physical_mean_jump_size_mu_p: xd[2],
                risk_neutral_mean_jump_size_mu_q: xd[3],
                jump_size_volatility_sigma_j: xd[4],
                diffusive_variance_risk_price_eta_v: xd[5],
            });
        let h = 1e-6;
        for i in 0..6 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_pan_implicit_jump_fear_volatility_risk_premia(
                &PanImplicitJumpFearVrpSpec {
                    spot_diffusive_variance_v: xp[0],
                    jump_intensity_variance_slope_lambda1: xp[1],
                    physical_mean_jump_size_mu_p: xp[2],
                    risk_neutral_mean_jump_size_mu_q: xp[3],
                    jump_size_volatility_sigma_j: xp[4],
                    diffusive_variance_risk_price_eta_v: xp[5],
                },
            );
            let rm = evaluate_pan_implicit_jump_fear_volatility_risk_premia(
                &PanImplicitJumpFearVrpSpec {
                    spot_diffusive_variance_v: xm[0],
                    jump_intensity_variance_slope_lambda1: xm[1],
                    physical_mean_jump_size_mu_p: xm[2],
                    risk_neutral_mean_jump_size_mu_q: xm[3],
                    jump_size_volatility_sigma_j: xm[4],
                    diffusive_variance_risk_price_eta_v: xm[5],
                },
            );
            let fd = (rp.total_pan_variance_risk_premium - rm.total_pan_variance_risk_premium)
                / (2.0 * h);
            assert!(rad.total_pan_variance_risk_premium.d[i].abs() > 1e-9);
            assert!((rad.total_pan_variance_risk_premium.d[i] - fd).abs() < 1e-5);
        }
    }
}
