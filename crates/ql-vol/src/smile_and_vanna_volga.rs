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
