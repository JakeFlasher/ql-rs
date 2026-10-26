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
