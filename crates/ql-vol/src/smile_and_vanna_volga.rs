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
