mod evaluate_clark_wystup_fx_vanna_volga_smile_module {
    use ql_ad::Real;
    use ql_math::{inv_norm_cdf, norm_cdf};

    #[derive(Clone, Copy, Debug)]
    pub struct ClarkWystupFxVannaVolgaSmileSpec<T: Real = f64> {
        pub spot_fx: T,
        pub domestic_rate: T,
        pub foreign_rate: T,
        pub expiry_years: T,
        pub atm_vol: T,
        pub risk_reversal_25d: T,
        pub butterfly_25d: T,
        pub target_strike: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct ClarkWystupFxVannaVolgaSmileReport<T: Real = f64> {
        pub vol_25d_put: T,
        pub vol_25d_call: T,
        pub pillar_reconstruction_error: T,
        pub strike_25d_put: T,
        pub strike_atm_dns: T,
        pub strike_25d_call: T,
        pub vanna_volga_weights: [T; 3],
        pub weights_partition_of_unity_error: T,
        pub vanna_volga_implied_vol: T,
        pub vanna_volga_call_price: T,
        pub flat_atm_bs_call_price: T,
        pub smile_overhedge_adjustment: T,
    }

    fn gk_call<T: Real>(s: T, k: T, rd: T, rf: T, sigma: T, t: T) -> T {
        let st = sigma * t.sqrt();
        let d1 = ((s / k).ln() + (rd - rf + T::from_f64(0.5) * sigma * sigma) * t) / st;
        let d2 = d1 - st;
        s * (-rf * t).exp() * norm_cdf(d1) - k * (-rd * t).exp() * norm_cdf(d2)
    }

    pub fn evaluate_clark_wystup_fx_vanna_volga_smile<T: Real>(
        s: &ClarkWystupFxVannaVolgaSmileSpec<T>,
    ) -> ClarkWystupFxVannaVolgaSmileReport<T> {
        let t = s.expiry_years.max_val(T::from_f64(1e-6));
        let half = T::from_f64(0.5);
        let sig_c = s.atm_vol + s.butterfly_25d + half * s.risk_reversal_25d;
        let sig_p = s.atm_vol + s.butterfly_25d - half * s.risk_reversal_25d;
        let rec1 = (sig_c - sig_p) - s.risk_reversal_25d;
        let rec2 = half * (sig_c + sig_p) - s.atm_vol - s.butterfly_25d;
        let rec_err = (rec1 * rec1 + rec2 * rec2 + T::from_f64(1e-30)).sqrt();
        let fwd = s.spot_fx * ((s.domestic_rate - s.foreign_rate) * t).exp();
        let alpha_d = inv_norm_cdf(
            (T::from_f64(0.25) * (s.foreign_rate * t).exp())
                .max_val(T::from_f64(0.01))
                .min_val(T::from_f64(0.49)),
        );
        let k_atm = fwd * (half * s.atm_vol * s.atm_vol * t).exp();
        let k1 = fwd * (alpha_d * sig_p * t.sqrt() + half * sig_p * sig_p * t).exp();
        let k3 = fwd * (-alpha_d * sig_c * t.sqrt() + half * sig_c * sig_c * t).exp();
        let x1 = k1.ln();
        let x2 = k_atm.ln();
        let x3 = k3.ln();
        let x = s.target_strike.ln();
        let w1 = (x - x2) * (x - x3) / ((x1 - x2) * (x1 - x3));
        let w2 = (x - x1) * (x - x3) / ((x2 - x1) * (x2 - x3));
        let w3 = (x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2));
        let w_diff = T::one() - (w1 + w2 + w3);
        let w_err = (w_diff * w_diff + T::from_f64(1e-30)).sqrt();
        let sig_vv = w1 * sig_p + w2 * s.atm_vol + w3 * sig_c;
        let flat_c = gk_call(
            s.spot_fx,
            s.target_strike,
            s.domestic_rate,
            s.foreign_rate,
            s.atm_vol,
            t,
        );
        let vv_c = gk_call(
            s.spot_fx,
            s.target_strike,
            s.domestic_rate,
            s.foreign_rate,
            sig_vv,
            t,
        );
        ClarkWystupFxVannaVolgaSmileReport {
            vol_25d_put: sig_p,
            vol_25d_call: sig_c,
            pillar_reconstruction_error: rec_err,
            strike_25d_put: k1,
            strike_atm_dns: k_atm,
            strike_25d_call: k3,
            vanna_volga_weights: [w1, w2, w3],
            weights_partition_of_unity_error: w_err,
            vanna_volga_implied_vol: sig_vv,
            vanna_volga_call_price: vv_c,
            flat_atm_bs_call_price: flat_c,
            smile_overhedge_adjustment: vv_c - flat_c,
        }
    }

    #[cfg(test)]
    mod domain_suite_47_tests {
        use super::*;
        use ql_ad::Dual;
        fn spec<T: Real>(x: T) -> ClarkWystupFxVannaVolgaSmileSpec<T> {
            ClarkWystupFxVannaVolgaSmileSpec {
                spot_fx: T::from_f64(1.10) + x,
                domestic_rate: T::from_f64(0.03),
                foreign_rate: T::from_f64(0.015),
                expiry_years: T::from_f64(1.0),
                atm_vol: T::from_f64(0.10),
                risk_reversal_25d: T::from_f64(0.012),
                butterfly_25d: T::from_f64(0.004),
                target_strike: T::from_f64(1.14),
            }
        }
        #[test]
        fn suite47_fx_invariants_and_dual_fd() {
            let r = evaluate_clark_wystup_fx_vanna_volga_smile(&spec(0.0));
            assert!(r.pillar_reconstruction_error < 1e-12);
            assert!(r.weights_partition_of_unity_error < 1e-12);
            assert!(r.strike_25d_put < r.strike_atm_dns && r.strike_atm_dns < r.strike_25d_call);
            let h = 1e-6;
            let fd = (evaluate_clark_wystup_fx_vanna_volga_smile(&spec(h)).vanna_volga_call_price
                - evaluate_clark_wystup_fx_vanna_volga_smile(&spec(-h)).vanna_volga_call_price)
                / (2.0 * h);
            let ad = evaluate_clark_wystup_fx_vanna_volga_smile(&spec(Dual::<1>::variable(0.0, 0)))
                .vanna_volga_call_price;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_clark_wystup_fx_vanna_volga_smile_module::*;

mod evaluate_clark_fx_forward_volatility_agreement_module {
    // Clark (2011) & Wystup (2017) FX Forward-Volatility Agreement (`FVA`),
    // Forward-Starting Smile Roll-Down & Calendar-Spread No-Arbitrage Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct ClarkFxForwardVolAgreementSpec<T: Real = f64> {
        pub spot_fx: T,
        pub domestic_rate_rd: T,
        pub foreign_rate_rf: T,
        pub front_expiry_t1: T,
        pub back_expiry_t2: T,
        pub atm_vol_t1: T,
        pub atm_vol_t2: T,
        pub smile_convexity_adj_bps: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct ClarkFxForwardVolAgreementReport<T: Real = f64> {
        pub total_variance_t1: T,
        pub total_variance_t2: T,
        pub calendar_variance_margin: T,
        pub forward_atm_implied_vol: T,
        pub fva_fair_vol_strike: T,
        pub forward_start_atm_straddle_pct: T,
    }

    pub fn evaluate_clark_fx_forward_volatility_agreement<T: Real>(
        spec: &ClarkFxForwardVolAgreementSpec<T>,
    ) -> ClarkFxForwardVolAgreementReport<T> {
        let t1 = spec.front_expiry_t1.max_val(T::from_f64(0.05));
        let t2 = spec.back_expiry_t2.max_val(t1 + T::from_f64(0.05));
        let dt = t2 - t1;

        let w1 = spec.atm_vol_t1 * spec.atm_vol_t1 * t1;
        let w2 = spec.atm_vol_t2 * spec.atm_vol_t2 * t2;
        let cal_margin = w2 - w1;

        // Forward implied variance sigma_fwd^2(T1, T2) = (sigma(T2)^2 T2 - sigma(T1)^2 T1) / (T2 - T1)
        let fwd_var = (cal_margin / dt).max_val(T::from_f64(1e-8));
        let fwd_vol = fwd_var.sqrt();

        // Smile-convexity-adjusted FVA strike (Clark Ch. 8):
        let fva_strike = fwd_vol + spec.smile_convexity_adj_bps * T::from_f64(1e-4);

        // Forward-starting ATM-forward straddle setting strike at T1 equal to F(T1, T2):
        // Discounted to t=0 in % of initial spot S0:
        let std = fva_strike * dt.sqrt();
        let d1 = T::from_f64(0.5) * std;
        let d2 = -d1;
        let df_d_t1 = (-spec.domestic_rate_rd * t1).exp();
        let df_f_dt = (-spec.foreign_rate_rf * dt).exp();
        let straddle_pct =
            T::from_f64(2.0) * df_d_t1 * df_f_dt * (norm_cdf(d1) - norm_cdf(d2)) * spec.spot_fx;

        ClarkFxForwardVolAgreementReport {
            total_variance_t1: w1,
            total_variance_t2: w2,
            calendar_variance_margin: cal_margin,
            forward_atm_implied_vol: fwd_vol,
            fva_fair_vol_strike: fva_strike,
            forward_start_atm_straddle_pct: straddle_pct,
        }
    }

    #[cfg(test)]
    mod domain_suite_51_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn clark_fx_forward_vol_agreement_invariants_and_ad_parity() {
            let spec = ClarkFxForwardVolAgreementSpec {
                spot_fx: 1.08,
                domestic_rate_rd: 0.045,
                foreign_rate_rf: 0.030,
                front_expiry_t1: 0.5,
                back_expiry_t2: 1.5,
                atm_vol_t1: 0.085,
                atm_vol_t2: 0.098,
                smile_convexity_adj_bps: -18.0,
            };
            let rep = evaluate_clark_fx_forward_volatility_agreement(&spec);
            assert!(rep.calendar_variance_margin > 0.0);
            assert!(rep.forward_atm_implied_vol > spec.atm_vol_t2);
            assert!(rep.forward_start_atm_straddle_pct > 0.0);

            let h = 1e-6;
            let r_ad =
                evaluate_clark_fx_forward_volatility_agreement(&ClarkFxForwardVolAgreementSpec {
                    spot_fx: Dual::<1>::from_value(spec.spot_fx),
                    domestic_rate_rd: Dual::<1>::from_value(spec.domestic_rate_rd),
                    foreign_rate_rf: Dual::<1>::from_value(spec.foreign_rate_rf),
                    front_expiry_t1: Dual::<1>::from_value(spec.front_expiry_t1),
                    back_expiry_t2: Dual::<1>::from_value(spec.back_expiry_t2),
                    atm_vol_t1: Dual::<1>::from_value(spec.atm_vol_t1),
                    atm_vol_t2: Dual::<1>::variable(spec.atm_vol_t2, 0),
                    smile_convexity_adj_bps: Dual::<1>::from_value(spec.smile_convexity_adj_bps),
                });
            let mut sp = spec;
            sp.atm_vol_t2 += h;
            let mut sm = spec;
            sm.atm_vol_t2 -= h;
            let fd = (evaluate_clark_fx_forward_volatility_agreement(&sp)
                .forward_start_atm_straddle_pct
                - evaluate_clark_fx_forward_volatility_agreement(&sm)
                    .forward_start_atm_straddle_pct)
                / (2.0 * h);
            assert!((r_ad.forward_start_atm_straddle_pct.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_clark_fx_forward_volatility_agreement_module::*;

mod evaluate_wystup_malz_fx_broker_fly_smile_inversion_module {
    // Wystup (2017) & Malz (1997) FX 25-Delta Risk-Reversal (`RR25`) & Butterfly (`BF25`)
    // Broker-Fly vs. Smile-Strangle Strike Inversion Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct WystupMalzFxSmileInversionSpec<T: Real = f64> {
        pub spot_fx_s0: T,
        pub domestic_rate_rd: T,
        pub foreign_rate_rf: T,
        pub atm_dns_vol_sigma_atm: T,
        pub risk_reversal_25d_rr: T,
        pub butterfly_25d_bf: T,
        pub expiry_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct WystupMalzFxSmileInversionReport<T: Real = f64> {
        pub forward_fx_rate: T,
        pub atm_delta_neutral_strike_k_dns: T,
        pub implied_25d_call_vol: T,
        pub implied_25d_put_vol: T,
        pub inverted_25d_call_strike_k25c: T,
        pub inverted_25d_put_strike_k25p: T,
        pub smile_strangle_minus_broker_fly_bps: T,
        pub call_25d_target_delta_error: T,
    }

    pub fn evaluate_wystup_malz_fx_broker_fly_smile_inversion<T: Real>(
        spec: &WystupMalzFxSmileInversionSpec<T>,
    ) -> WystupMalzFxSmileInversionReport<T> {
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-4));
        let rd = spec.domestic_rate_rd;
        let rf = spec.foreign_rate_rf;
        let t = spec.expiry_years_t.max_val(T::from_f64(0.02));
        let sqrt_t = t.sqrt();
        let df_f = (-rf * t).exp();

        let fwd = s0 * ((rd - rf) * t).exp();
        let sig_atm = spec.atm_dns_vol_sigma_atm.max_val(T::from_f64(1e-4));
        let rr = spec.risk_reversal_25d_rr;
        let bf = spec.butterfly_25d_bf;

        // Unadjusted delta-neutral straddle (DNS) ATM strike:
        let k_dns = fwd * (T::from_f64(0.5) * sig_atm * sig_atm * t).exp();

        // 25-delta call and put implied volatilities from ATM, 25D RR, and 25D BF:
        let sig_25c = (sig_atm + bf + T::from_f64(0.5) * rr).max_val(T::from_f64(1e-4));
        let sig_25p = (sig_atm + bf - T::from_f64(0.5) * rr).max_val(T::from_f64(1e-4));

        // Exact spot-delta inversion for 25D call and 25D put strikes:
        // Spot Call Delta = e^{-r_f T} Phi(d_1) = 0.25 => d_1 = Phi^{-1}(0.25 * e^{r_f T})
        // Solve via 6 Newton steps in d_1 starting from -0.67448975:
        let target_cdf = (T::from_f64(0.25) / df_f)
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.95));
        let mut alpha_c = T::from_f64(0.6744897501960817);
        for _ in 0..6 {
            let cdf_val = norm_cdf(-alpha_c);
            let pdf_val =
                (-T::from_f64(0.5) * alpha_c * alpha_c).exp() * T::from_f64(0.3989422804014327);
            alpha_c += (cdf_val - target_cdf) / pdf_val.max_val(T::from_f64(1e-6));
        }

        let k_25c =
            fwd * (alpha_c * sig_25c * sqrt_t + T::from_f64(0.5) * sig_25c * sig_25c * t).exp();
        let k_25p =
            fwd * (-alpha_c * sig_25p * sqrt_t + T::from_f64(0.5) * sig_25p * sig_25p * t).exp();

        let d1_c =
            ((fwd / k_25c).ln() + T::from_f64(0.5) * sig_25c * sig_25c * t) / (sig_25c * sqrt_t);
        let delta_err = (df_f * norm_cdf(d1_c) - T::from_f64(0.25)).abs();

        // Wystup (2017) Smile-Strangle vs. Broker-Fly convexity adjustment in bps:
        let fly_conv_bps = T::from_f64(0.125) * (alpha_c * alpha_c) * (rr * rr) * sqrt_t / sig_atm
            * T::from_f64(1e4);

        WystupMalzFxSmileInversionReport {
            forward_fx_rate: fwd,
            atm_delta_neutral_strike_k_dns: k_dns,
            implied_25d_call_vol: sig_25c,
            implied_25d_put_vol: sig_25p,
            inverted_25d_call_strike_k25c: k_25c,
            inverted_25d_put_strike_k25p: k_25p,
            smile_strangle_minus_broker_fly_bps: fly_conv_bps,
            call_25d_target_delta_error: delta_err,
        }
    }

    #[cfg(test)]
    mod domain_suite_54_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn wystup_malz_fx_broker_fly_smile_inversion_invariants_and_ad_parity() {
            let spec = WystupMalzFxSmileInversionSpec {
                spot_fx_s0: 1.085,
                domestic_rate_rd: 0.045,
                foreign_rate_rf: 0.030,
                atm_dns_vol_sigma_atm: 0.085,
                risk_reversal_25d_rr: -0.012,
                butterfly_25d_bf: 0.0035,
                expiry_years_t: 0.5,
            };
            let rep = evaluate_wystup_malz_fx_broker_fly_smile_inversion(&spec);
            assert!(rep.inverted_25d_call_strike_k25c > rep.atm_delta_neutral_strike_k_dns);
            assert!(rep.inverted_25d_put_strike_k25p < rep.atm_delta_neutral_strike_k_dns);
            assert!(rep.call_25d_target_delta_error < 1e-10);

            let h = 1e-6;
            let r_ad = evaluate_wystup_malz_fx_broker_fly_smile_inversion(
                &WystupMalzFxSmileInversionSpec {
                    spot_fx_s0: Dual::<1>::from_value(spec.spot_fx_s0),
                    domestic_rate_rd: Dual::<1>::from_value(spec.domestic_rate_rd),
                    foreign_rate_rf: Dual::<1>::from_value(spec.foreign_rate_rf),
                    atm_dns_vol_sigma_atm: Dual::<1>::variable(spec.atm_dns_vol_sigma_atm, 0),
                    risk_reversal_25d_rr: Dual::<1>::from_value(spec.risk_reversal_25d_rr),
                    butterfly_25d_bf: Dual::<1>::from_value(spec.butterfly_25d_bf),
                    expiry_years_t: Dual::<1>::from_value(spec.expiry_years_t),
                },
            );
            let mut sp = spec;
            sp.atm_dns_vol_sigma_atm += h;
            let mut sm = spec;
            sm.atm_dns_vol_sigma_atm -= h;
            let fd = (evaluate_wystup_malz_fx_broker_fly_smile_inversion(&sp)
                .inverted_25d_call_strike_k25c
                - evaluate_wystup_malz_fx_broker_fly_smile_inversion(&sm)
                    .inverted_25d_call_strike_k25c)
                / (2.0 * h);
            assert!((r_ad.inverted_25d_call_strike_k25c.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_wystup_malz_fx_broker_fly_smile_inversion_module::*;
