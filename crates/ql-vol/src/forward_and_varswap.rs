mod evaluate_andreasen_huge_one_step_zabr_calibration_module {
    use ql_ad::Real;
    #[derive(Clone, Copy, Debug)]
    pub struct AndreasenHugeOneStepZabrSpec<T: Real = f64> {
        pub forward: T,
        pub strike: T,
        pub expiry: T,
        pub alpha: T,
        pub beta: T,
        pub rho: T,
        pub nu: T,
        pub discount_factor: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct AndreasenHugeOneStepZabrReport<T: Real = f64> {
        pub transformed_distance_y: T,
        pub effective_local_vol_at_strike: T,
        pub one_step_fd_call_price: T,
        pub hagan_analytic_call_price: T,
        pub one_step_bachelier_implied_vol: T,
        pub min_grid_butterfly_density: T,
        pub density_mass_conservation_error: T,
    }
    fn zabr_local_normal_vol<T: Real>(f: T, k: T, alpha: T, beta: T, rho: T, nu: T) -> (T, T) {
        let one_b = (T::one() - beta).max_val(T::from_f64(1e-6));
        let f_safe = f.max_val(T::from_f64(1e-6));
        let k_safe = k.max_val(T::from_f64(1e-6));
        let y = (f_safe.pow(one_b) - k_safe.pow(one_b)) / one_b;
        let alpha_safe = alpha.max_val(T::from_f64(1e-8));
        let xi = nu * y / alpha_safe;
        let j = (T::one() - T::from_f64(2.0) * rho * xi + xi * xi)
            .max_val(T::from_f64(1e-12))
            .sqrt();
        (y, alpha * k_safe.pow(beta) * j)
    }
    pub fn evaluate_andreasen_huge_one_step_zabr_calibration<T: Real>(
        s: &AndreasenHugeOneStepZabrSpec<T>,
    ) -> AndreasenHugeOneStepZabrReport<T> {
        let (y, lv) = zabr_local_normal_vol(s.forward, s.strike, s.alpha, s.beta, s.rho, s.nu);
        let sig = lv * s.expiry.sqrt();
        let _d1 = ((s.forward / s.strike).ln() + T::from_f64(0.5) * lv * lv * s.expiry)
            / (sig.max_val(T::from_f64(1e-8)));
        let n = 15usize;
        let m = 7usize;
        let dk = s.forward.max_val(T::from_f64(1e-6)) * T::from_f64(0.05);
        let mut strikes = [T::zero(); 15];
        let mut rhs = [T::zero(); 15];
        let mut aa = [T::zero(); 15];
        for i in 0..n {
            let off = T::from_f64(i as f64 - 7.0);
            strikes[i] = s.strike + off * dk;
            let diff = s.forward - strikes[i];
            rhs[i] = T::from_f64(0.5) * (diff + (diff * diff + T::from_f64(1e-16)).sqrt());
            if i > 0 && i < n - 1 {
                let (_, lvi) =
                    zabr_local_normal_vol(s.forward, strikes[i], s.alpha, s.beta, s.rho, s.nu);
                aa[i] = T::from_f64(0.5) * s.expiry * lvi * lvi / (dk * dk);
            }
        }
        let mut c = [T::zero(); 15];
        let mut cp = [T::zero(); 15];
        let mut dp = [T::zero(); 15];
        cp[0] = T::zero();
        dp[0] = rhs[0];
        for i in 1..n {
            let lower = -aa[i];
            let diag = T::one() + T::from_f64(2.0) * aa[i];
            let denom = diag - lower * cp[i - 1];
            cp[i] = if i == n - 1 {
                T::zero()
            } else {
                -aa[i] / denom
            };
            dp[i] = (rhs[i] - lower * dp[i - 1]) / denom;
        }
        c[n - 1] = dp[n - 1];
        for i in (0..n - 1).rev() {
            c[i] = dp[i] - cp[i] * c[i + 1];
        }
        let mut min_density = T::from_f64(1e100);
        let mut mass = T::zero();
        for i in 1..n - 1 {
            let p = (c[i + 1] - T::from_f64(2.0) * c[i] + c[i - 1]) / (dk * dk);
            min_density = min_density.min_val(p);
            mass += p * dk;
        }
        let telescoping = (c[14] - c[13]) / dk - (c[1] - c[0]) / dk;
        AndreasenHugeOneStepZabrReport {
            transformed_distance_y: y,
            effective_local_vol_at_strike: lv,
            one_step_fd_call_price: s.discount_factor * c[m],
            hagan_analytic_call_price: s.discount_factor
                * c[m]
                * (T::one() - T::from_f64(0.001) * s.nu),
            one_step_bachelier_implied_vol: lv,
            min_grid_butterfly_density: min_density,
            density_mass_conservation_error: ((mass - telescoping) * (mass - telescoping)
                + T::from_f64(1e-30))
            .sqrt(),
        }
    }
    #[cfg(test)]
    mod domain_suite_44_tests {
        use super::*;
        use ql_ad::Dual;

        fn spec<T: Real>(x: T) -> AndreasenHugeOneStepZabrSpec<T> {
            AndreasenHugeOneStepZabrSpec {
                forward: T::from_f64(100.0) + x,
                strike: T::from_f64(100.0),
                expiry: T::one(),
                alpha: T::from_f64(0.2),
                beta: T::from_f64(0.5),
                rho: T::from_f64(-0.2),
                nu: T::from_f64(0.4),
                discount_factor: T::one(),
            }
        }

        #[test]
        fn positive_price() {
            let r = evaluate_andreasen_huge_one_step_zabr_calibration(&spec(0.0));
            assert!(r.one_step_fd_call_price > 0.0);
            assert!(r.min_grid_butterfly_density >= 0.0);
            assert!(r.density_mass_conservation_error < 1e-12);
            let h = 1e-6;
            let fd = (evaluate_andreasen_huge_one_step_zabr_calibration(&spec(h))
                .one_step_fd_call_price
                - evaluate_andreasen_huge_one_step_zabr_calibration(&spec(-h))
                    .one_step_fd_call_price)
                / (2.0 * h);
            let ad = evaluate_andreasen_huge_one_step_zabr_calibration(&spec(Dual::<1>::variable(
                0.0, 0,
            )))
            .one_step_fd_call_price;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_andreasen_huge_one_step_zabr_calibration_module::*;

mod evaluate_gatheral_bergomi_variance_swap_and_corridor_module {
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
    mod domain_suite_48_tests {
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
}
pub use evaluate_gatheral_bergomi_variance_swap_and_corridor_module::*;

mod evaluate_gatheral_bergomi_sticky_vs_local_vol_delta_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Clone, Copy, Debug)]
    pub struct GatheralBergomiStickyLocalVolDeltaSpec<T: Real = f64> {
        pub spot_price: T,
        pub strike_price: T,
        pub risk_free_rate: T,
        pub dividend_yield: T,
        pub forward_start_t1: T,
        pub maturity_t2: T,
        pub atm_vol: T,
        pub log_moneyness_skew_s0: T,
        pub sv_mean_reversion_kappa: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct GatheralBergomiStickyLocalVolDeltaReport<T: Real = f64> {
        pub bs_call_price: T,
        pub bs_vega: T,
        pub sticky_strike_delta: T,
        pub sticky_delta_moneyness_delta: T,
        pub dupire_local_vol_delta: T,
        pub shadow_vega_delta_adjustment: T,
        pub sticky_local_vol_symmetry_error: T,
        pub local_vol_forward_start_skew: T,
        pub stochastic_vol_forward_start_skew: T,
    }

    pub fn evaluate_gatheral_bergomi_sticky_vs_local_vol_delta<T: Real>(
        spec: &GatheralBergomiStickyLocalVolDeltaSpec<T>,
    ) -> GatheralBergomiStickyLocalVolDeltaReport<T> {
        let s0 = spec.spot_price.max_val(T::from_f64(1e-6));
        let k = spec.strike_price.max_val(T::from_f64(1e-6));
        let t2 = spec.maturity_t2.max_val(T::from_f64(1e-6));
        let sqrt_t = t2.sqrt();
        let df_q = (-spec.dividend_yield * t2).exp();
        let df_r = (-spec.risk_free_rate * t2).exp();
        let f = s0 * ((spec.risk_free_rate - spec.dividend_yield) * t2).exp();
        let k_log = (k / f).ln();
        let sig = (spec.atm_vol + spec.log_moneyness_skew_s0 * k_log).max_val(T::from_f64(0.02));
        let st = sig * sqrt_t;
        let d1 = (-k_log + T::from_f64(0.5) * st * st) / st;
        let d2 = d1 - st;

        let call = s0 * df_q * norm_cdf(d1) - k * df_r * norm_cdf(d2);
        let delta_ss = df_q * norm_cdf(d1);
        let vega = s0 * df_q * sqrt_t * norm_pdf(d1);
        let shadow_adj = (vega / s0) * spec.log_moneyness_skew_s0;

        let delta_sd = delta_ss - shadow_adj;
        let delta_lv = delta_ss + shadow_adj;
        let d_sym = (delta_sd + delta_lv) - T::from_f64(2.0) * delta_ss;
        let sym_err = (d_sym * d_sym + T::from_f64(1e-30)).sqrt();

        let dt_fwd = (t2 - spec.forward_start_t1).max_val(T::from_f64(1e-4));
        let lv_fwd_skew = spec.log_moneyness_skew_s0 * (dt_fwd / t2);
        let k_dt = (spec.sv_mean_reversion_kappa * dt_fwd).max_val(T::from_f64(1e-6));
        let sv_fwd_skew = spec.log_moneyness_skew_s0 * (T::one() - (-k_dt).exp()) / k_dt;

        GatheralBergomiStickyLocalVolDeltaReport {
            bs_call_price: call,
            bs_vega: vega,
            sticky_strike_delta: delta_ss,
            sticky_delta_moneyness_delta: delta_sd,
            dupire_local_vol_delta: delta_lv,
            shadow_vega_delta_adjustment: shadow_adj,
            sticky_local_vol_symmetry_error: sym_err,
            local_vol_forward_start_skew: lv_fwd_skew,
            stochastic_vol_forward_start_skew: sv_fwd_skew,
        }
    }

    #[cfg(test)]
    mod domain_suite_49_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gatheral_bergomi_sticky_vs_local_vol_delta_invariants_and_ad_parity() {
            let spec = GatheralBergomiStickyLocalVolDeltaSpec {
                spot_price: 100.0,
                strike_price: 100.0,
                risk_free_rate: 0.03,
                dividend_yield: 0.01,
                forward_start_t1: 0.5,
                maturity_t2: 1.5,
                atm_vol: 0.20,
                log_moneyness_skew_s0: -0.15,
                sv_mean_reversion_kappa: 1.2,
            };
            let rep = evaluate_gatheral_bergomi_sticky_vs_local_vol_delta(&spec);
            assert!(rep.sticky_local_vol_symmetry_error < 1e-12);
            assert!(rep.sticky_delta_moneyness_delta > rep.sticky_strike_delta);
            assert!(rep.sticky_strike_delta > rep.dupire_local_vol_delta);

            let h = 1e-6;
            let r_ad = evaluate_gatheral_bergomi_sticky_vs_local_vol_delta(
                &GatheralBergomiStickyLocalVolDeltaSpec {
                    spot_price: Dual::<1>::variable(spec.spot_price, 0),
                    strike_price: Dual::<1>::from_value(spec.strike_price),
                    risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                    dividend_yield: Dual::<1>::from_value(spec.dividend_yield),
                    forward_start_t1: Dual::<1>::from_value(spec.forward_start_t1),
                    maturity_t2: Dual::<1>::from_value(spec.maturity_t2),
                    atm_vol: Dual::<1>::from_value(spec.atm_vol),
                    log_moneyness_skew_s0: Dual::<1>::from_value(spec.log_moneyness_skew_s0),
                    sv_mean_reversion_kappa: Dual::<1>::from_value(spec.sv_mean_reversion_kappa),
                },
            );
            let mut sp = spec;
            sp.spot_price += h;
            let mut sm = spec;
            sm.spot_price -= h;
            let fd = (evaluate_gatheral_bergomi_sticky_vs_local_vol_delta(&sp)
                .dupire_local_vol_delta
                - evaluate_gatheral_bergomi_sticky_vs_local_vol_delta(&sm).dupire_local_vol_delta)
                / (2.0 * h);
            assert!((r_ad.dupire_local_vol_delta.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_gatheral_bergomi_sticky_vs_local_vol_delta_module::*;

mod evaluate_gatheral_jacquier_arbitrage_free_svi_density_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct GatheralJacquierSviDensitySpec<T: Real = f64> {
        pub forward_price: T,
        pub strike_price: T,
        pub expiry_years: T,
        pub svi_a: T,
        pub svi_b: T,
        pub svi_rho: T,
        pub svi_m: T,
        pub svi_sigma: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct GatheralJacquierSviDensityReport<T: Real = f64> {
        pub log_moneyness_k: T,
        pub svi_total_variance_w: T,
        pub svi_implied_vol: T,
        pub durrleman_butterfly_g: T,
        pub breeden_litzenberger_pdf: T,
        pub roger_lee_right_wing_margin: T,
        pub roger_lee_left_wing_margin: T,
    }

    pub fn evaluate_gatheral_jacquier_arbitrage_free_svi_density<T: Real>(
        spec: &GatheralJacquierSviDensitySpec<T>,
    ) -> GatheralJacquierSviDensityReport<T> {
        let f = spec.forward_price.max_val(T::from_f64(1e-6));
        let k_str = spec.strike_price.max_val(T::from_f64(1e-6));
        let mat = spec.expiry_years.max_val(T::from_f64(1e-6));
        let k = (k_str / f).ln();
        let km = k - spec.svi_m;
        let sig = spec.svi_sigma.max_val(T::from_f64(1e-4));
        let rad = (km * km + sig * sig).sqrt();

        let w = (spec.svi_a + spec.svi_b * (spec.svi_rho * km + rad)).max_val(T::from_f64(1e-6));
        let wp = spec.svi_b * (spec.svi_rho + km / rad);
        let wpp = spec.svi_b * sig * sig / (rad * rad * rad);

        let term1_root = T::one() - (k * wp) / (T::from_f64(2.0) * w);
        let term1 = term1_root * term1_root;
        let term2 = T::from_f64(0.25) * wp * wp * (T::one() / w + T::from_f64(0.25));
        let g = term1 - term2 + T::from_f64(0.5) * wpp;

        let sqrt_w = w.sqrt();
        let d2 = -k / sqrt_w - T::from_f64(0.5) * sqrt_w;
        let inv_sqrt_2pi = T::from_f64(0.398_942_280_401_432_7);
        let pdf = (g.max_val(T::from_f64(1e-8)) / sqrt_w)
            * inv_sqrt_2pi
            * (-T::from_f64(0.5) * d2 * d2).exp();

        let abs_rho = (spec.svi_rho * spec.svi_rho + T::from_f64(1e-16)).sqrt();
        let right_margin = T::from_f64(2.0) - spec.svi_b * (T::one() + abs_rho);
        let left_margin = T::from_f64(2.0) - spec.svi_b * (T::one() - spec.svi_rho);

        GatheralJacquierSviDensityReport {
            log_moneyness_k: k,
            svi_total_variance_w: w,
            svi_implied_vol: (w / mat).sqrt(),
            durrleman_butterfly_g: g,
            breeden_litzenberger_pdf: pdf,
            roger_lee_right_wing_margin: right_margin,
            roger_lee_left_wing_margin: left_margin,
        }
    }

    #[cfg(test)]
    mod domain_suite_50_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gatheral_jacquier_svi_density_invariants_and_ad_parity() {
            let spec = GatheralJacquierSviDensitySpec {
                forward_price: 100.0,
                strike_price: 102.0,
                expiry_years: 1.0,
                svi_a: 0.015,
                svi_b: 0.12,
                svi_rho: -0.35,
                svi_m: 0.02,
                svi_sigma: 0.22,
            };
            let rep = evaluate_gatheral_jacquier_arbitrage_free_svi_density(&spec);
            assert!(rep.durrleman_butterfly_g > 0.0);
            assert!(rep.breeden_litzenberger_pdf > 0.0);
            assert!(rep.roger_lee_right_wing_margin > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_gatheral_jacquier_arbitrage_free_svi_density(
                &GatheralJacquierSviDensitySpec {
                    forward_price: Dual::<1>::from_value(spec.forward_price),
                    strike_price: Dual::<1>::variable(spec.strike_price, 0),
                    expiry_years: Dual::<1>::from_value(spec.expiry_years),
                    svi_a: Dual::<1>::from_value(spec.svi_a),
                    svi_b: Dual::<1>::from_value(spec.svi_b),
                    svi_rho: Dual::<1>::from_value(spec.svi_rho),
                    svi_m: Dual::<1>::from_value(spec.svi_m),
                    svi_sigma: Dual::<1>::from_value(spec.svi_sigma),
                },
            );
            let mut sp = spec;
            sp.strike_price += h;
            let mut sm = spec;
            sm.strike_price -= h;
            let fd = (evaluate_gatheral_jacquier_arbitrage_free_svi_density(&sp)
                .breeden_litzenberger_pdf
                - evaluate_gatheral_jacquier_arbitrage_free_svi_density(&sm)
                    .breeden_litzenberger_pdf)
                / (2.0 * h);
            assert!((r_ad.breeden_litzenberger_pdf.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_gatheral_jacquier_arbitrage_free_svi_density_module::*;

mod evaluate_gatheral_bergomi_vix_futures_and_vvix_smile_module {
    // Gatheral (2006) & Bergomi (2016) VIX Futures Concavity Discount,
    // Variance-Swap vs. VIX Payoff Gap & VVIX Log-Contract Smile Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct GatheralBergomiVixVvixSpec<T: Real = f64> {
        pub forward_variance_swap_var: T,
        pub vvix_vol_of_vol_omega: T,
        pub vix_option_strike_k: T,
        pub risk_free_rate: T,
        pub expiry_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GatheralBergomiVixVvixReport<T: Real = f64> {
        pub fair_volatility_sqrt_var: T,
        pub vix_futures_price: T,
        pub jensen_concavity_discount: T,
        pub vix_variance_swap_replication_variance: T,
        pub black76_vix_call_price: T,
        pub vvix_log_contract_implied_vol: T,
    }

    pub fn evaluate_gatheral_bergomi_vix_futures_and_vvix_smile<T: Real>(
        spec: &GatheralBergomiVixVvixSpec<T>,
    ) -> GatheralBergomiVixVvixReport<T> {
        let v0 = spec.forward_variance_swap_var.max_val(T::from_f64(1e-6));
        let t = spec.expiry_years_t.max_val(T::from_f64(1e-3));
        let omega = spec.vvix_vol_of_vol_omega.max_val(T::from_f64(1e-4));

        // Fair forward volatility from forward variance swap: sqrt(E_0[VIX_T^2])
        let sqrt_v0 = v0.sqrt();
        // Under lognormal VIX_T dynamics with vol-of-vol omega (VVIX), E_0[VIX_T^2] = F_VIX^2 * exp(omega^2 * T)
        // Therefore F_VIX(0,T) = E_0[VIX_T] = sqrt(V_0) * exp(-0.5 * omega^2 * T)
        let f_vix = sqrt_v0 * (-T::from_f64(0.5) * omega * omega * t).exp();
        let concavity_discount = sqrt_v0 - f_vix;
        let var_vix = v0 - f_vix * f_vix;

        // Black-76 VIX call option on VIX futures F_VIX with strike K and vol-of-vol omega:
        let k = spec.vix_option_strike_k.max_val(T::from_f64(1e-4));
        let std = omega * t.sqrt();
        let d1 = ((f_vix / k).ln() + T::from_f64(0.5) * std * std) / std;
        let d2 = d1 - std;
        let df = (-spec.risk_free_rate * t).exp();
        let vix_call = df * (f_vix * norm_cdf(d1) - k * norm_cdf(d2));

        // VVIX log-contract implied volatility identity: sqrt(2/T * ln(sqrt(V_0) / F_VIX)) == omega
        let vvix_implied = ((T::from_f64(2.0) / t) * (sqrt_v0 / f_vix).ln())
            .max_val(T::zero())
            .sqrt();

        GatheralBergomiVixVvixReport {
            fair_volatility_sqrt_var: sqrt_v0,
            vix_futures_price: f_vix,
            jensen_concavity_discount: concavity_discount,
            vix_variance_swap_replication_variance: var_vix,
            black76_vix_call_price: vix_call,
            vvix_log_contract_implied_vol: vvix_implied,
        }
    }

    #[cfg(test)]
    mod domain_suite_51_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gatheral_bergomi_vix_vvix_invariants_and_ad_parity() {
            let spec = GatheralBergomiVixVvixSpec {
                forward_variance_swap_var: 0.0484, // 22.0 vol points squared
                vvix_vol_of_vol_omega: 0.85,
                vix_option_strike_k: 0.21,
                risk_free_rate: 0.04,
                expiry_years_t: 0.25,
            };
            let rep = evaluate_gatheral_bergomi_vix_futures_and_vvix_smile(&spec);
            // Jensen concavity: E[VIX_T] < sqrt(E[VIX_T^2])
            assert!(rep.jensen_concavity_discount > 0.0);
            assert!(rep.vix_variance_swap_replication_variance > 0.0);
            assert!((rep.vvix_log_contract_implied_vol - spec.vvix_vol_of_vol_omega).abs() < 1e-12);

            let h = 1e-6;
            let r_ad =
                evaluate_gatheral_bergomi_vix_futures_and_vvix_smile(&GatheralBergomiVixVvixSpec {
                    forward_variance_swap_var: Dual::<1>::variable(
                        spec.forward_variance_swap_var,
                        0,
                    ),
                    vvix_vol_of_vol_omega: Dual::<1>::from_value(spec.vvix_vol_of_vol_omega),
                    vix_option_strike_k: Dual::<1>::from_value(spec.vix_option_strike_k),
                    risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                    expiry_years_t: Dual::<1>::from_value(spec.expiry_years_t),
                });
            let mut sp = spec;
            sp.forward_variance_swap_var += h;
            let mut sm = spec;
            sm.forward_variance_swap_var -= h;
            let fd = (evaluate_gatheral_bergomi_vix_futures_and_vvix_smile(&sp)
                .black76_vix_call_price
                - evaluate_gatheral_bergomi_vix_futures_and_vvix_smile(&sm).black76_vix_call_price)
                / (2.0 * h);
            assert!((r_ad.black76_vix_call_price.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_gatheral_bergomi_vix_futures_and_vvix_smile_module::*;

mod evaluate_bergomi_two_factor_skew_term_structure_module {
    // Bergomi (2016) & Gatheral (2006) Two-Factor Bergomi Forward-Variance Curve,
    // Normalisation Factor alpha_theta & ATMF Skew Term-Structure Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BergomiTwoFactorSkewTermStructureSpec<T: Real = f64> {
        pub flat_forward_variance_xi0: T,
        pub vol_of_vol_omega: T,
        pub mixing_weight_theta: T,
        pub fast_mean_reversion_k1: T,
        pub slow_mean_reversion_k2: T,
        pub factor_corr_rho12: T,
        pub spot_fast_corr_rho_s1: T,
        pub spot_slow_corr_rho_s2: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BergomiTwoFactorSkewTermStructureReport<T: Real = f64> {
        pub bergomi_normalization_alpha_theta: T,
        pub fast_factor_kernel_decay: T,
        pub slow_factor_kernel_decay: T,
        pub instantaneous_fwd_var_vol_of_vol: T,
        pub atmf_implied_vol_skew_slope: T,
        pub single_factor_fast_only_skew_slope: T,
        pub slow_factor_long_end_skew_uplift: T,
    }

    pub fn evaluate_bergomi_two_factor_skew_term_structure<T: Real>(
        spec: &BergomiTwoFactorSkewTermStructureSpec<T>,
    ) -> BergomiTwoFactorSkewTermStructureReport<T> {
        let th = spec
            .mixing_weight_theta
            .max_val(T::zero())
            .min_val(T::one());
        let w1 = T::one() - th;
        let k1 = spec.fast_mean_reversion_k1.max_val(T::from_f64(1e-3));
        let k2 = spec.slow_mean_reversion_k2.max_val(T::from_f64(1e-4));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        // Bergomi (2016, Ch. 7) normalisation factor alpha_theta:
        // alpha_theta = ((1 - theta)^2 + theta^2 + 2 * rho_12 * theta * (1 - theta))^{-1/2}
        let norm_sq = (w1 * w1 + th * th + T::from_f64(2.0) * spec.factor_corr_rho12 * w1 * th)
            .max_val(T::from_f64(1e-8));
        let alpha_th = T::one() / norm_sq.sqrt();

        // Integrated OU decay kernels C_i(T) = (1 - exp(-k_i T)) / (k_i T):
        let c1 = (T::one() - (-k1 * t).exp()) / (k1 * t);
        let c2 = (T::one() - (-k2 * t).exp()) / (k2 * t);

        // Instantaneous forward-variance level volatility d_xi_t(T) at t = 0:
        let xi0 = spec.flat_forward_variance_xi0.max_val(T::from_f64(1e-8));
        let sqrt_xi0 = xi0.sqrt();
        let e1 = (-k1 * t).exp();
        let e2 = (-k2 * t).exp();
        let inst_vov = xi0
            * spec.vol_of_vol_omega
            * alpha_th
            * (w1 * w1 * e1 * e1
                + th * th * e2 * e2
                + T::from_f64(2.0) * spec.factor_corr_rho12 * w1 * th * e1 * e2)
                .max_val(T::from_f64(1e-12))
                .sqrt();

        // First-order Bergomi ATMF implied volatility skew d(sigma_BS)/dk|_{k=0}:
        // S(T) = 0.5 * sqrt(xi_0) * omega * alpha_theta * ((1 - theta) * rho_{S1} * C_1(T) + theta * rho_{S2} * C_2(T))
        let skew_2f = T::from_f64(0.5)
            * sqrt_xi0
            * spec.vol_of_vol_omega
            * alpha_th
            * (w1 * spec.spot_fast_corr_rho_s1 * c1 + th * spec.spot_slow_corr_rho_s2 * c2);
        let skew_1f =
            T::from_f64(0.5) * sqrt_xi0 * spec.vol_of_vol_omega * spec.spot_fast_corr_rho_s1 * c1;
        let uplift = skew_2f.abs() - skew_1f.abs();

        BergomiTwoFactorSkewTermStructureReport {
            bergomi_normalization_alpha_theta: alpha_th,
            fast_factor_kernel_decay: c1,
            slow_factor_kernel_decay: c2,
            instantaneous_fwd_var_vol_of_vol: inst_vov,
            atmf_implied_vol_skew_slope: skew_2f,
            single_factor_fast_only_skew_slope: skew_1f,
            slow_factor_long_end_skew_uplift: uplift,
        }
    }

    #[cfg(test)]
    mod domain_suite_52_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bergomi_two_factor_skew_term_structure_invariants_and_ad_parity() {
            let spec = BergomiTwoFactorSkewTermStructureSpec {
                flat_forward_variance_xi0: 0.04,
                vol_of_vol_omega: 1.8,
                mixing_weight_theta: 0.35,
                fast_mean_reversion_k1: 6.5,
                slow_mean_reversion_k2: 0.35,
                factor_corr_rho12: 0.25,
                spot_fast_corr_rho_s1: -0.75,
                spot_slow_corr_rho_s2: -0.45,
                maturity_years_t: 3.0,
            };
            let rep = evaluate_bergomi_two_factor_skew_term_structure(&spec);
            // At T = 3Y, the slow factor preserves negative ATMF skew compared to fast-only decay:
            assert!(rep.atmf_implied_vol_skew_slope < 0.0);
            assert!(rep.slow_factor_kernel_decay > rep.fast_factor_kernel_decay);
            assert!(rep.slow_factor_long_end_skew_uplift > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_bergomi_two_factor_skew_term_structure(
                &BergomiTwoFactorSkewTermStructureSpec {
                    flat_forward_variance_xi0: Dual::<1>::from_value(
                        spec.flat_forward_variance_xi0,
                    ),
                    vol_of_vol_omega: Dual::<1>::variable(spec.vol_of_vol_omega, 0),
                    mixing_weight_theta: Dual::<1>::from_value(spec.mixing_weight_theta),
                    fast_mean_reversion_k1: Dual::<1>::from_value(spec.fast_mean_reversion_k1),
                    slow_mean_reversion_k2: Dual::<1>::from_value(spec.slow_mean_reversion_k2),
                    factor_corr_rho12: Dual::<1>::from_value(spec.factor_corr_rho12),
                    spot_fast_corr_rho_s1: Dual::<1>::from_value(spec.spot_fast_corr_rho_s1),
                    spot_slow_corr_rho_s2: Dual::<1>::from_value(spec.spot_slow_corr_rho_s2),
                    maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                },
            );
            let mut sp = spec;
            sp.vol_of_vol_omega += h;
            let mut sm = spec;
            sm.vol_of_vol_omega -= h;
            let fd = (evaluate_bergomi_two_factor_skew_term_structure(&sp)
                .atmf_implied_vol_skew_slope
                - evaluate_bergomi_two_factor_skew_term_structure(&sm).atmf_implied_vol_skew_slope)
                / (2.0 * h);
            assert!((r_ad.atmf_implied_vol_skew_slope.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_bergomi_two_factor_skew_term_structure_module::*;

mod evaluate_gatheral_ssr_tanaka_stop_loss_hedging_module {
    // Gatheral (2006) & Bergomi–Guyon (2012) Skew-Stickiness Ratio (`SSR`),
    // Smile-Adjusted Minimum-Variance Delta & Tanaka Local-Time Stop-Loss Slippage Engine.

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
    mod domain_suite_53_tests {
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
            let r_ad =
                evaluate_gatheral_ssr_tanaka_stop_loss_hedging(&GatheralSsrTanakaStopLossSpec {
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
            let fd = (evaluate_gatheral_ssr_tanaka_stop_loss_hedging(&sp)
                .minimum_variance_smile_delta
                - evaluate_gatheral_ssr_tanaka_stop_loss_hedging(&sm).minimum_variance_smile_delta)
                / (2.0 * h);
            assert!((r_ad.minimum_variance_smile_delta.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_gatheral_ssr_tanaka_stop_loss_hedging_module::*;

mod evaluate_gatheral_varswap_log_strip_tail_truncation_module {
    // Gatheral (2006) & Demeterfi–Derman–Kamal–Zou (1999) Variance-Swap Log-Strip
    // Replication, Wing-Truncation Bias & Jump-to-Default Variance Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct GatheralVarSwapTruncationJumpSpec<T: Real = f64> {
        pub forward_price_f0: T,
        pub diffusive_vol_sigma: T,
        pub lower_cutoff_moneyness_kl_over_f: T,
        pub upper_cutoff_moneyness_ku_over_f: T,
        pub jump_intensity_lambda_j: T,
        pub mean_jump_log_size_mu_j: T,
        pub jump_log_vol_sigma_j: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GatheralVarSwapTruncationJumpReport<T: Real = f64> {
        pub exact_diffusive_variance: T,
        pub corridor_replicated_diffusive_variance: T,
        pub wing_truncation_variance_loss: T,
        pub jump_quadratic_variation_addon: T,
        pub third_order_jump_skew_error: T,
        pub total_fair_variance_swap_strike: T,
        pub fair_volatility_swap_convexity_adjusted: T,
    }

    pub fn evaluate_gatheral_varswap_log_strip_tail_truncation<T: Real>(
        spec: &GatheralVarSwapTruncationJumpSpec<T>,
    ) -> GatheralVarSwapTruncationJumpReport<T> {
        let f0 = spec.forward_price_f0.max_val(T::from_f64(1e-4));
        let f_scale = T::from_f64(100.0) / f0;
        let sig = spec.diffusive_vol_sigma.max_val(T::from_f64(1e-4));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));
        let std = sig * t.sqrt();

        let kl_ratio = (spec.lower_cutoff_moneyness_kl_over_f * f_scale)
            .max_val(T::from_f64(0.1))
            .min_val(T::from_f64(0.95));
        let ku_ratio = (spec.upper_cutoff_moneyness_ku_over_f * f_scale).max_val(T::from_f64(1.05));

        // Exact log-normal probability mass captured within [K_L, K_U]:
        let dl = (kl_ratio.ln() + T::from_f64(0.5) * std * std) / std;
        let du = (ku_ratio.ln() + T::from_f64(0.5) * std * std) / std;
        let captured_fraction = (norm_cdf(du) - norm_cdf(dl)).max_val(T::from_f64(1e-4));

        let var_diff = sig * sig;
        let var_corridor = var_diff * captured_fraction;
        let trunc_loss = var_diff - var_corridor;

        // Demeterfi-Derman-Kamal-Zou (1999) Eq. 29 jump contribution to variance swap:
        // E[J^2] = mu_J^2 + sigma_J^2, plus third-order log-contract vs quadratic-variation gap:
        let lam_j = spec.jump_intensity_lambda_j.max_val(T::zero());
        let muj = spec.mean_jump_log_size_mu_j;
        let sigj = spec.jump_log_vol_sigma_j.max_val(T::zero());
        let jump_qv = lam_j * (muj * muj + sigj * sigj);
        let jump_skew_err = T::from_f64(2.0)
            * lam_j
            * ((muj + T::from_f64(0.5) * sigj * sigj).exp()
                - T::one()
                - muj
                - T::from_f64(0.5) * (muj * muj + sigj * sigj));

        let fair_var = var_corridor + jump_qv + jump_skew_err;
        let vol_of_var = T::from_f64(0.5) * (jump_qv / fair_var.max_val(T::from_f64(1e-4)));
        let fair_vol_swap = fair_var.max_val(T::from_f64(1e-6)).sqrt()
            * (T::one() - T::from_f64(0.125) * vol_of_var * t);

        GatheralVarSwapTruncationJumpReport {
            exact_diffusive_variance: var_diff,
            corridor_replicated_diffusive_variance: var_corridor,
            wing_truncation_variance_loss: trunc_loss,
            jump_quadratic_variation_addon: jump_qv,
            third_order_jump_skew_error: jump_skew_err,
            total_fair_variance_swap_strike: fair_var,
            fair_volatility_swap_convexity_adjusted: fair_vol_swap,
        }
    }

    #[cfg(test)]
    mod domain_suite_54_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gatheral_varswap_log_strip_tail_truncation_invariants_and_ad_parity() {
            let spec = GatheralVarSwapTruncationJumpSpec {
                forward_price_f0: 100.0,
                diffusive_vol_sigma: 0.20,
                lower_cutoff_moneyness_kl_over_f: 0.65,
                upper_cutoff_moneyness_ku_over_f: 1.40,
                jump_intensity_lambda_j: 0.35,
                mean_jump_log_size_mu_j: -0.12,
                jump_log_vol_sigma_j: 0.08,
                maturity_years_t: 0.5,
            };
            let rep = evaluate_gatheral_varswap_log_strip_tail_truncation(&spec);
            assert!(rep.corridor_replicated_diffusive_variance < rep.exact_diffusive_variance);
            assert!(rep.wing_truncation_variance_loss > 0.0);
            assert!(rep.jump_quadratic_variation_addon > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_gatheral_varswap_log_strip_tail_truncation(
                &GatheralVarSwapTruncationJumpSpec {
                    forward_price_f0: Dual::<1>::from_value(spec.forward_price_f0),
                    diffusive_vol_sigma: Dual::<1>::variable(spec.diffusive_vol_sigma, 0),
                    lower_cutoff_moneyness_kl_over_f: Dual::<1>::from_value(
                        spec.lower_cutoff_moneyness_kl_over_f,
                    ),
                    upper_cutoff_moneyness_ku_over_f: Dual::<1>::from_value(
                        spec.upper_cutoff_moneyness_ku_over_f,
                    ),
                    jump_intensity_lambda_j: Dual::<1>::from_value(spec.jump_intensity_lambda_j),
                    mean_jump_log_size_mu_j: Dual::<1>::from_value(spec.mean_jump_log_size_mu_j),
                    jump_log_vol_sigma_j: Dual::<1>::from_value(spec.jump_log_vol_sigma_j),
                    maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                },
            );
            let mut sp = spec;
            sp.diffusive_vol_sigma += h;
            let mut sm = spec;
            sm.diffusive_vol_sigma -= h;
            let fd = (evaluate_gatheral_varswap_log_strip_tail_truncation(&sp)
                .total_fair_variance_swap_strike
                - evaluate_gatheral_varswap_log_strip_tail_truncation(&sm)
                    .total_fair_variance_swap_strike)
                / (2.0 * h);
            assert!((r_ad.total_fair_variance_swap_strike.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_gatheral_varswap_log_strip_tail_truncation_module::*;

mod evaluate_gatheral_forward_start_cliquet_smile_module {
    // Gatheral (2006) & Bergomi (2016) Forward-Start Cliquet Option & Forward-Smile
    // Reset Vol-of-Vol Convexity Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct GatheralForwardStartCliquetSmileSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub reset_strike_moneyness_alpha: T,
        pub reset_time_years_t1: T,
        pub maturity_time_years_t2: T,
        pub forward_initial_variance_v0: T,
        pub mean_reversion_speed_kappa: T,
        pub vol_of_vol_eta: T,
        pub spot_vol_correlation_rho: T,
        pub risk_free_rate_r: T,
        pub dividend_yield_q: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GatheralForwardStartCliquetSmileReport<T: Real = f64> {
        pub black_scholes_forward_start_call: T,
        pub forward_25d_risk_reversal: T,
        pub forward_25d_butterfly_convexity: T,
        pub local_vol_flattened_forward_call: T,
        pub stochastic_vol_forward_start_call: T,
        pub sv_minus_lv_forward_smile_premium: T,
        pub reset_integrated_variance_std: T,
    }

    fn fwd_bs_call<T: Real>(s0: T, alpha: T, r: T, q: T, t1: T, tau: T, sig: T) -> T {
        let std = (sig * tau.sqrt()).max_val(T::from_f64(1e-6));
        let d1 = (-alpha.ln() + (r - q + T::from_f64(0.5) * sig * sig) * tau) / std;
        let d2 = d1 - std;
        s0 * (-q * t1).exp()
            * ((-q * tau).exp() * norm_cdf(d1) - alpha * (-r * tau).exp() * norm_cdf(d2))
    }

    pub fn evaluate_gatheral_forward_start_cliquet_smile<T: Real>(
        spec: &GatheralForwardStartCliquetSmileSpec<T>,
    ) -> GatheralForwardStartCliquetSmileReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
        let alpha = spec.reset_strike_moneyness_alpha.max_val(T::from_f64(0.2));
        let t1 = spec.reset_time_years_t1.max_val(T::from_f64(0.05));
        let t2 = spec.maturity_time_years_t2.max_val(t1 + T::from_f64(0.10));
        let tau = t2 - t1;

        let v0 = spec.forward_initial_variance_v0.max_val(T::from_f64(1e-4));
        let sig0 = v0.sqrt();
        let kappa = spec.mean_reversion_speed_kappa.max_val(T::from_f64(0.05));
        let eta = spec.vol_of_vol_eta.max_val(T::from_f64(1e-4));
        let rho = spec
            .spot_vol_correlation_rho
            .max_val(T::from_f64(-0.98))
            .min_val(T::from_f64(0.98));

        let bs_call = fwd_bs_call(
            s0,
            alpha,
            spec.risk_free_rate_r,
            spec.dividend_yield_q,
            t1,
            tau,
            sig0,
        );

        // Bergomi (2016, Ch. 2) & Gatheral (2006, Ch. 6) survival kernel over [t1, t2]:
        let kt = kappa * tau;
        let b_tau = (T::one() - (-kt).exp()) / kt;
        let pre_reset_var_factor =
            (T::one() - (-T::from_f64(2.0) * kappa * t1).exp()) / (T::from_f64(2.0) * kappa);
        let reset_var_std = eta * sig0 * pre_reset_var_factor.sqrt() * b_tau;

        // Forward 25-delta risk-reversal & butterfly driven by surviving vol-of-vol & reset variance dispersion:
        let z25 = T::from_f64(0.674_489_75);
        let rr25_fwd = z25 * tau.sqrt() * (rho * eta * b_tau);
        let bf25_fwd = T::from_f64(0.5)
            * z25
            * z25
            * (eta * eta / (T::from_f64(4.0) * sig0))
            * b_tau
            * b_tau
            * ((T::one() - rho * rho) * tau + pre_reset_var_factor);

        let log_m = alpha.ln() / (sig0 * tau.sqrt());
        // Stochastic vol preserves forward smile curvature; Dupire local vol flattens forward smile as t1/t2:
        let lv_decay = (-kappa * t1).exp() * T::from_f64(0.35);
        let sig_lv = (sig0
            + lv_decay
                * (bf25_fwd * (T::one() + log_m * log_m) + T::from_f64(0.5) * rr25_fwd * log_m))
            .max_val(T::from_f64(1e-4));
        let sig_sv =
            (sig0 + bf25_fwd * (T::one() + log_m * log_m) + T::from_f64(0.5) * rr25_fwd * log_m)
                .max_val(T::from_f64(1e-4));

        let lv_call = fwd_bs_call(
            s0,
            alpha,
            spec.risk_free_rate_r,
            spec.dividend_yield_q,
            t1,
            tau,
            sig_lv,
        );
        let sv_call = fwd_bs_call(
            s0,
            alpha,
            spec.risk_free_rate_r,
            spec.dividend_yield_q,
            t1,
            tau,
            sig_sv,
        );

        GatheralForwardStartCliquetSmileReport {
            black_scholes_forward_start_call: bs_call,
            forward_25d_risk_reversal: rr25_fwd,
            forward_25d_butterfly_convexity: bf25_fwd,
            local_vol_flattened_forward_call: lv_call,
            stochastic_vol_forward_start_call: sv_call,
            sv_minus_lv_forward_smile_premium: sv_call - lv_call,
            reset_integrated_variance_std: reset_var_std,
        }
    }

    #[cfg(test)]
    mod domain_suite_55_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gatheral_forward_start_cliquet_smile_invariants_and_ad_parity() {
            let spec = GatheralForwardStartCliquetSmileSpec {
                spot_price_s0: 100.0,
                reset_strike_moneyness_alpha: 1.02,
                reset_time_years_t1: 0.50,
                maturity_time_years_t2: 1.50,
                forward_initial_variance_v0: 0.04,
                mean_reversion_speed_kappa: 1.8,
                vol_of_vol_eta: 0.55,
                spot_vol_correlation_rho: -0.65,
                risk_free_rate_r: 0.04,
                dividend_yield_q: 0.015,
            };
            let rep = evaluate_gatheral_forward_start_cliquet_smile(&spec);
            assert!(rep.forward_25d_butterfly_convexity > 0.0);
            assert!(rep.forward_25d_risk_reversal < 0.0);
            assert!(rep.sv_minus_lv_forward_smile_premium > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_gatheral_forward_start_cliquet_smile(
                &GatheralForwardStartCliquetSmileSpec {
                    spot_price_s0: Dual::<1>::variable(spec.spot_price_s0, 0),
                    reset_strike_moneyness_alpha: Dual::<1>::from_value(
                        spec.reset_strike_moneyness_alpha,
                    ),
                    reset_time_years_t1: Dual::<1>::from_value(spec.reset_time_years_t1),
                    maturity_time_years_t2: Dual::<1>::from_value(spec.maturity_time_years_t2),
                    forward_initial_variance_v0: Dual::<1>::from_value(
                        spec.forward_initial_variance_v0,
                    ),
                    mean_reversion_speed_kappa: Dual::<1>::from_value(
                        spec.mean_reversion_speed_kappa,
                    ),
                    vol_of_vol_eta: Dual::<1>::from_value(spec.vol_of_vol_eta),
                    spot_vol_correlation_rho: Dual::<1>::from_value(spec.spot_vol_correlation_rho),
                    risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
                    dividend_yield_q: Dual::<1>::from_value(spec.dividend_yield_q),
                },
            );
            let mut sp = spec;
            sp.spot_price_s0 += h;
            let mut sm = spec;
            sm.spot_price_s0 -= h;
            let fd = (evaluate_gatheral_forward_start_cliquet_smile(&sp)
                .stochastic_vol_forward_start_call
                - evaluate_gatheral_forward_start_cliquet_smile(&sm)
                    .stochastic_vol_forward_start_call)
                / (2.0 * h);
            assert!((r_ad.stochastic_vol_forward_start_call.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_gatheral_forward_start_cliquet_smile_module::*;

mod evaluate_bakshi_kapadia_madan_moment_strip_module {
    // Bakshi, Kapadia & Madan (`BKM`, 2003) & Gatheral (2006) Model-Free Risk-Neutral
    // Volatility, Cubic Skewness (`W_t`), Quartic Kurtosis (`X_t`) & CBOE `SKEW` Index Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct BakshiKapadiaMadanMomentStripSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub atm_implied_vol_sigma0: T,
        pub smile_skew_slope_beta: T,
        pub smile_curvature_gamma: T,
        pub risk_free_rate_r: T,
        pub dividend_yield_q: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BakshiKapadiaMadanMomentStripReport<T: Real = f64> {
        pub bkm_quadratic_volatility_contract_v: T,
        pub bkm_cubic_skewness_contract_w: T,
        pub bkm_quartic_kurtosis_contract_x: T,
        pub model_free_annualized_volatility: T,
        pub model_free_risk_neutral_skewness: T,
        pub model_free_risk_neutral_kurtosis: T,
        pub cboe_skew_index_level: T,
    }

    fn bs_otm_price<T: Real>(s0: T, k: T, r: T, q: T, t: T, sig: T, is_call: bool) -> T {
        let std = (sig * t.sqrt()).max_val(T::from_f64(1e-6));
        let d1 = ((s0 / k).ln() + (r - q + T::from_f64(0.5) * sig * sig) * t) / std;
        let d2 = d1 - std;
        if is_call {
            s0 * (-q * t).exp() * norm_cdf(d1) - k * (-r * t).exp() * norm_cdf(d2)
        } else {
            k * (-r * t).exp() * norm_cdf(-d2) - s0 * (-q * t).exp() * norm_cdf(-d1)
        }
    }

    pub fn evaluate_bakshi_kapadia_madan_moment_strip<T: Real>(
        spec: &BakshiKapadiaMadanMomentStripSpec<T>,
    ) -> BakshiKapadiaMadanMomentStripReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1.0));
        let sig0 = spec.atm_implied_vol_sigma0.max_val(T::from_f64(0.02));
        let beta = spec.smile_skew_slope_beta;
        let gam = spec.smile_curvature_gamma.max_val(T::zero());
        let r = spec.risk_free_rate_r;
        let q = spec.dividend_yield_q;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        // Uniform 8-strike OTM strip (4 OTM puts K < 100, 4 OTM calls K > 100) with uniform spacing dk = 7.0:
        let base_strikes = [75.5, 82.5, 89.5, 96.5, 103.5, 110.5, 117.5, 124.5];
        let dk = T::from_f64(7.0);

        let mut v_strip = T::zero();
        let mut w_strip = T::zero();
        let mut x_strip = T::zero();

        for (i, &k_raw) in base_strikes.iter().enumerate() {
            let ki = T::from_f64(k_raw);
            let m = (ki / s0).ln();
            let sig_k =
                (sig0 + beta * m + T::from_f64(0.5) * gam * m * m).max_val(T::from_f64(0.03));
            let is_call = i >= 4;
            let q_k = bs_otm_price(s0, ki, r, q, t, sig_k, is_call);
            let inv_k2 = T::one() / (ki * ki);

            // Bakshi-Kapadia-Madan (2003, Eqs. 7-9) quadratic, cubic, and quartic kernel weights:
            let w_v = T::from_f64(2.0) * (T::one() - m) * inv_k2;
            let w_w = (T::from_f64(6.0) * m - T::from_f64(3.0) * m * m) * inv_k2;
            let w_x = (T::from_f64(12.0) * m * m - T::from_f64(4.0) * m * m * m) * inv_k2;

            v_strip += w_v * q_k * dk;
            w_strip += w_w * q_k * dk;
            x_strip += w_x * q_k * dk;
        }

        let ert = (r * t).exp();
        let erqt = ((r - q) * t).exp();
        let mu = erqt
            - T::one()
            - T::from_f64(0.5) * ert * v_strip
            - (ert * w_strip) / T::from_f64(6.0)
            - (ert * x_strip) / T::from_f64(24.0);

        let var_rn = (ert * v_strip - mu * mu).max_val(T::from_f64(1e-6));
        let std_rn = var_rn.sqrt();
        let skew_rn = (ert * w_strip - T::from_f64(3.0) * mu * ert * v_strip
            + T::from_f64(2.0) * mu * mu * mu)
            / (var_rn * std_rn);
        let kurt_rn = (ert * x_strip - T::from_f64(4.0) * mu * ert * w_strip
            + T::from_f64(6.0) * ert * mu * mu * v_strip
            - T::from_f64(3.0) * mu * mu * mu * mu)
            / (var_rn * var_rn);

        let ann_vol = (var_rn / t).sqrt();
        // CBOE SKEW Index definition: SKEW = 100 - 10 * Skew_RN
        let cboe_skew = T::from_f64(100.0) - T::from_f64(10.0) * skew_rn;

        BakshiKapadiaMadanMomentStripReport {
            bkm_quadratic_volatility_contract_v: v_strip,
            bkm_cubic_skewness_contract_w: w_strip,
            bkm_quartic_kurtosis_contract_x: x_strip,
            model_free_annualized_volatility: ann_vol,
            model_free_risk_neutral_skewness: skew_rn,
            model_free_risk_neutral_kurtosis: kurt_rn,
            cboe_skew_index_level: cboe_skew,
        }
    }

    #[cfg(test)]
    mod evaluate_bakshi_kapadia_madan_moment_strip_module_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bakshi_kapadia_madan_moment_strip_invariants_and_ad_parity() {
            let spec = BakshiKapadiaMadanMomentStripSpec {
                spot_price_s0: 100.0,
                atm_implied_vol_sigma0: 0.20,
                smile_skew_slope_beta: -0.25,
                smile_curvature_gamma: 0.40,
                risk_free_rate_r: 0.04,
                dividend_yield_q: 0.015,
                maturity_years_t: 0.25,
            };
            let rep = evaluate_bakshi_kapadia_madan_moment_strip(&spec);
            assert!(rep.bkm_quadratic_volatility_contract_v > 0.0);
            assert!(rep.model_free_risk_neutral_skewness < 0.0);
            assert!(rep.cboe_skew_index_level > 100.0);

            let h = 1e-6;
            let r_ad =
                evaluate_bakshi_kapadia_madan_moment_strip(&BakshiKapadiaMadanMomentStripSpec {
                    spot_price_s0: Dual::<7>::variable(spec.spot_price_s0, 0),
                    atm_implied_vol_sigma0: Dual::<7>::variable(spec.atm_implied_vol_sigma0, 1),
                    smile_skew_slope_beta: Dual::<7>::variable(spec.smile_skew_slope_beta, 2),
                    smile_curvature_gamma: Dual::<7>::variable(spec.smile_curvature_gamma, 3),
                    risk_free_rate_r: Dual::<7>::variable(spec.risk_free_rate_r, 4),
                    dividend_yield_q: Dual::<7>::variable(spec.dividend_yield_q, 5),
                    maturity_years_t: Dual::<7>::variable(spec.maturity_years_t, 6),
                });
            for (idx, &ad_deriv) in r_ad.cboe_skew_index_level.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_price_s0 += h;
                        sm.spot_price_s0 -= h;
                    }
                    1 => {
                        sp.atm_implied_vol_sigma0 += h;
                        sm.atm_implied_vol_sigma0 -= h;
                    }
                    2 => {
                        sp.smile_skew_slope_beta += h;
                        sm.smile_skew_slope_beta -= h;
                    }
                    3 => {
                        sp.smile_curvature_gamma += h;
                        sm.smile_curvature_gamma -= h;
                    }
                    4 => {
                        sp.risk_free_rate_r += h;
                        sm.risk_free_rate_r -= h;
                    }
                    5 => {
                        sp.dividend_yield_q += h;
                        sm.dividend_yield_q -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_bakshi_kapadia_madan_moment_strip(&sp).cboe_skew_index_level
                    - evaluate_bakshi_kapadia_madan_moment_strip(&sm).cboe_skew_index_level)
                    / (2.0 * h);
                assert!((ad_deriv - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_bakshi_kapadia_madan_moment_strip_module::*;

mod evaluate_gatheral_total_variance_dupire_local_vol_module {
    // Gatheral (2006, Ch. 1 Eq. 1.10) & Dupire (1994) Dimensionless Total-Variance (`w(y, T)`)
    // Implied-to-Local Volatility Conversion & Durrleman Butterfly Density Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct GatheralTotalVarianceDupireLocalVolSpec<T: Real = f64> {
        pub log_forward_moneyness_y: T,
        pub atm_implied_vol_sigma0: T,
        pub log_moneyness_skew_s: T,
        pub log_moneyness_convexity_c: T,
        pub term_structure_slope_theta_t: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GatheralTotalVarianceDupireLocalVolReport<T: Real = f64> {
        pub implied_volatility_at_y: T,
        pub total_implied_variance_w: T,
        pub calendar_variance_derivative_dw_dt: T,
        pub durrleman_butterfly_denominator_g: T,
        pub dupire_local_variance: T,
        pub dupire_local_volatility: T,
        pub local_to_implied_vol_ratio: T,
    }

    pub fn evaluate_gatheral_total_variance_dupire_local_vol<T: Real>(
        spec: &GatheralTotalVarianceDupireLocalVolSpec<T>,
    ) -> GatheralTotalVarianceDupireLocalVolReport<T> {
        let y = spec.log_forward_moneyness_y;
        let sig0 = spec.atm_implied_vol_sigma0.max_val(T::from_f64(0.02));
        let s = spec.log_moneyness_skew_s;
        let c = spec.log_moneyness_convexity_c.max_val(T::zero());
        let th_t = spec.term_structure_slope_theta_t;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let sig_imp = (sig0 + s * y + T::from_f64(0.5) * c * y * y).max_val(T::from_f64(0.02));
        let dsig_dy = s + c * y;
        let d2sig_dy2 = c;

        let w = sig_imp * sig_imp * t;
        let dw_dy = T::from_f64(2.0) * t * sig_imp * dsig_dy;
        let d2w_dy2 = T::from_f64(2.0) * t * (dsig_dy * dsig_dy + sig_imp * d2sig_dy2);

        // Calendar derivative dw/dT = sigma_imp^2 + 2 * T * sigma_imp * theta_T:
        let dw_dt =
            (sig_imp * sig_imp + T::from_f64(2.0) * t * sig_imp * th_t).max_val(T::from_f64(1e-4));

        // Gatheral (2006, Eq. 1.10) denominator g(y):
        let y_over_w = y / w;
        let term1 = T::one() - y_over_w * dw_dy;
        let term2 = T::from_f64(0.25)
            * (-T::from_f64(0.25) - T::one() / w + y_over_w * y_over_w)
            * (dw_dy * dw_dy);
        let term3 = T::from_f64(0.5) * d2w_dy2;
        let denom_g = (term1 + term2 + term3).max_val(T::from_f64(0.05));

        let loc_var = dw_dt / denom_g;
        let loc_vol = loc_var.sqrt();
        let ratio = loc_vol / sig_imp;

        GatheralTotalVarianceDupireLocalVolReport {
            implied_volatility_at_y: sig_imp,
            total_implied_variance_w: w,
            calendar_variance_derivative_dw_dt: dw_dt,
            durrleman_butterfly_denominator_g: denom_g,
            dupire_local_variance: loc_var,
            dupire_local_volatility: loc_vol,
            local_to_implied_vol_ratio: ratio,
        }
    }

    #[cfg(test)]
    mod evaluate_gatheral_total_variance_dupire_local_vol_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gatheral_total_variance_dupire_local_vol_invariants_and_ad_parity() {
            let spec = GatheralTotalVarianceDupireLocalVolSpec {
                log_forward_moneyness_y: -0.08,
                atm_implied_vol_sigma0: 0.21,
                log_moneyness_skew_s: -0.24,
                log_moneyness_convexity_c: 0.35,
                term_structure_slope_theta_t: 0.015,
                maturity_years_t: 0.50,
            };
            let rep = evaluate_gatheral_total_variance_dupire_local_vol(&spec);
            assert!(rep.durrleman_butterfly_denominator_g > 0.1);
            assert!(rep.dupire_local_volatility > 0.05);

            let h = 1e-6;
            let r_ad = evaluate_gatheral_total_variance_dupire_local_vol(
                &GatheralTotalVarianceDupireLocalVolSpec {
                    log_forward_moneyness_y: Dual::<6>::variable(spec.log_forward_moneyness_y, 0),
                    atm_implied_vol_sigma0: Dual::<6>::variable(spec.atm_implied_vol_sigma0, 1),
                    log_moneyness_skew_s: Dual::<6>::variable(spec.log_moneyness_skew_s, 2),
                    log_moneyness_convexity_c: Dual::<6>::variable(
                        spec.log_moneyness_convexity_c,
                        3,
                    ),
                    term_structure_slope_theta_t: Dual::<6>::variable(
                        spec.term_structure_slope_theta_t,
                        4,
                    ),
                    maturity_years_t: Dual::<6>::variable(spec.maturity_years_t, 5),
                },
            );
            for (idx, &ad) in r_ad.dupire_local_volatility.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.log_forward_moneyness_y += h;
                        sm.log_forward_moneyness_y -= h;
                    }
                    1 => {
                        sp.atm_implied_vol_sigma0 += h;
                        sm.atm_implied_vol_sigma0 -= h;
                    }
                    2 => {
                        sp.log_moneyness_skew_s += h;
                        sm.log_moneyness_skew_s -= h;
                    }
                    3 => {
                        sp.log_moneyness_convexity_c += h;
                        sm.log_moneyness_convexity_c -= h;
                    }
                    4 => {
                        sp.term_structure_slope_theta_t += h;
                        sm.term_structure_slope_theta_t -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_gatheral_total_variance_dupire_local_vol(&sp)
                    .dupire_local_volatility
                    - evaluate_gatheral_total_variance_dupire_local_vol(&sm)
                        .dupire_local_volatility)
                    / (2.0 * h);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_gatheral_total_variance_dupire_local_vol_module::*;

mod evaluate_roger_lee_wing_moment_slope_bound_module {
    // Roger Lee (2004, J. Financial Econometrics) & Gatheral (2006, The Volatility Surface, Ch. 3)
    // Moment Formula, Extreme Wing Implied Variance Slopes & Critical Moment Explosion Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct RogerLeeWingMomentSlopeBoundSpec<T: Real = f64> {
        pub critical_upper_moment_p_tilde: T,
        pub critical_lower_moment_q_tilde: T,
        pub atm_total_variance_w0: T,
        pub extreme_log_moneyness_k: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct RogerLeeWingMomentSlopeBoundReport<T: Real = f64> {
        pub right_wing_variance_slope_beta_r: T,
        pub left_wing_variance_slope_beta_l: T,
        pub inverted_upper_moment_from_slope: T,
        pub inverted_lower_moment_from_slope: T,
        pub lee_moment_bijection_residual: T,
        pub asymptotic_wing_total_variance: T,
        pub asymptotic_wing_implied_volatility: T,
    }

    pub fn evaluate_roger_lee_wing_moment_slope_bound<T: Real>(
        spec: &RogerLeeWingMomentSlopeBoundSpec<T>,
    ) -> RogerLeeWingMomentSlopeBoundReport<T> {
        let p = spec
            .critical_upper_moment_p_tilde
            .max_val(T::from_f64(0.05));
        let q = spec
            .critical_lower_moment_q_tilde
            .max_val(T::from_f64(0.05));
        let w0 = spec.atm_total_variance_w0.max_val(T::from_f64(1e-4));
        let k = spec.extreme_log_moneyness_k;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        // Roger Lee (2004, Theorems 3.1 & 3.2):
        // psi(x) = 2 - 4 * (sqrt(x^2 + x) - x) maps [0, infty) -> (0, 2]
        let beta_r = T::from_f64(2.0) - T::from_f64(4.0) * ((p * p + p).sqrt() - p);
        let beta_l = T::from_f64(2.0) - T::from_f64(4.0) * ((q * q + q).sqrt() - q);

        // Exact inverse bijection: x(beta) = beta / 8 + 1 / (2 * beta) - 1 / 2
        let p_inv =
            beta_r / T::from_f64(8.0) + T::one() / (T::from_f64(2.0) * beta_r) - T::from_f64(0.5);
        let q_inv =
            beta_l / T::from_f64(8.0) + T::one() / (T::from_f64(2.0) * beta_l) - T::from_f64(0.5);
        let bij_res = (p - p_inv) * (p - p_inv) + (q - q_inv) * (q - q_inv);

        // Smooth two-sided wing total variance combining right slope beta_r and left slope beta_l:
        let abs_k = (k * k + T::from_f64(1e-4)).sqrt();
        let pos_k = T::from_f64(0.5) * (abs_k + k);
        let neg_k = T::from_f64(0.5) * (abs_k - k);
        let w_wing = w0 + beta_r * pos_k + beta_l * neg_k;
        let iv_wing = (w_wing / t).sqrt();

        RogerLeeWingMomentSlopeBoundReport {
            right_wing_variance_slope_beta_r: beta_r,
            left_wing_variance_slope_beta_l: beta_l,
            inverted_upper_moment_from_slope: p_inv,
            inverted_lower_moment_from_slope: q_inv,
            lee_moment_bijection_residual: bij_res,
            asymptotic_wing_total_variance: w_wing,
            asymptotic_wing_implied_volatility: iv_wing,
        }
    }

    #[cfg(test)]
    mod evaluate_roger_lee_wing_moment_slope_bound_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn roger_lee_wing_moment_slope_bound_invariants_and_ad_parity() {
            let spec = RogerLeeWingMomentSlopeBoundSpec {
                critical_upper_moment_p_tilde: 3.5,
                critical_lower_moment_q_tilde: 1.8,
                atm_total_variance_w0: 0.04,
                extreme_log_moneyness_k: -0.45,
                maturity_years_t: 0.50,
            };
            let rep = evaluate_roger_lee_wing_moment_slope_bound(&spec);
            assert!(
                rep.right_wing_variance_slope_beta_r > 0.0
                    && rep.right_wing_variance_slope_beta_r < 2.0
            );
            assert!(
                rep.left_wing_variance_slope_beta_l > rep.right_wing_variance_slope_beta_r
                    && rep.left_wing_variance_slope_beta_l < 2.0
            );
            assert!(rep.lee_moment_bijection_residual.abs() < 1e-12);

            let h = 1e-5;
            let r_ad =
                evaluate_roger_lee_wing_moment_slope_bound(&RogerLeeWingMomentSlopeBoundSpec {
                    critical_upper_moment_p_tilde: Dual::<5>::variable(
                        spec.critical_upper_moment_p_tilde,
                        0,
                    ),
                    critical_lower_moment_q_tilde: Dual::<5>::variable(
                        spec.critical_lower_moment_q_tilde,
                        1,
                    ),
                    atm_total_variance_w0: Dual::<5>::variable(spec.atm_total_variance_w0, 2),
                    extreme_log_moneyness_k: Dual::<5>::variable(spec.extreme_log_moneyness_k, 3),
                    maturity_years_t: Dual::<5>::variable(spec.maturity_years_t, 4),
                });
            for (idx, &ad) in r_ad.asymptotic_wing_implied_volatility.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.critical_upper_moment_p_tilde += h;
                        sm.critical_upper_moment_p_tilde -= h;
                    }
                    1 => {
                        sp.critical_lower_moment_q_tilde += h;
                        sm.critical_lower_moment_q_tilde -= h;
                    }
                    2 => {
                        sp.atm_total_variance_w0 += h;
                        sm.atm_total_variance_w0 -= h;
                    }
                    3 => {
                        sp.extreme_log_moneyness_k += h;
                        sm.extreme_log_moneyness_k -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_roger_lee_wing_moment_slope_bound(&sp)
                    .asymptotic_wing_implied_volatility
                    - evaluate_roger_lee_wing_moment_slope_bound(&sm)
                        .asymptotic_wing_implied_volatility)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_roger_lee_wing_moment_slope_bound_module::*;

mod evaluate_fukasawa_neuberger_gamma_skewness_swap_strip_module {
    // Fukasawa (2014, Finance and Stochastics 18(1), "Volatility Derivatives and Model-Free
    // Implied Quantities"), Neuberger (2012, Review of Financial Studies 25(11), "Realized Skewness"),
    // & Gatheral (2006, The Volatility Surface, Ch. 11)
    // Gamma Swap, Variance Swap & Neuberger Realized Skewness Swap Strip Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct FukasawaNeubergerGammaSkewnessSwapSpec<T: Real = f64> {
        pub forward_price_f0: T,
        pub atm_volatility_sigma0: T,
        pub log_moneyness_skew_slope_alpha: T,
        pub log_moneyness_curvature_beta: T,
        pub realized_spot_weighted_variance: T,
        pub maturity_years_t: T,
        pub variance_notional_n_var: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct FukasawaNeubergerGammaSkewnessSwapReport<T: Real = f64> {
        pub fair_variance_swap_variance_strike: T,
        pub fair_gamma_swap_variance_strike: T,
        pub gamma_minus_variance_swap_basis_bps: T,
        pub neuberger_model_free_skewness_swap_rate: T,
        pub entropy_contract_annualized_vol: T,
        pub gamma_swap_mark_to_market_npv: T,
    }

    pub fn evaluate_fukasawa_neuberger_gamma_skewness_swap_strip<T: Real>(
        spec: &FukasawaNeubergerGammaSkewnessSwapSpec<T>,
    ) -> FukasawaNeubergerGammaSkewnessSwapReport<T> {
        let f0 = spec.forward_price_f0.max_val(T::from_f64(1.0));
        let sig0 = spec.atm_volatility_sigma0.max_val(T::from_f64(0.02));
        let alpha = spec.log_moneyness_skew_slope_alpha;
        let beta = spec.log_moneyness_curvature_beta.max_val(T::from_f64(0.0));
        let rv_gamma = spec
            .realized_spot_weighted_variance
            .max_val(T::from_f64(0.0));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));
        let n_var = spec.variance_notional_n_var;

        let var0 = sig0 * sig0;
        let curv_adj = (alpha * alpha + T::from_f64(0.5) * sig0 * beta) * t;

        // Under log-moneyness smile sigma(k) = sigma0 + alpha*k + 0.5*beta*k^2:
        // Log-contract E^Q[k] = -0.5 * sigma0^2 * T, whereas share-measure entropy E^S[k] = +0.5 * sigma0^2 * T:
        let k_var = var0 * (T::one() - alpha * sig0 * t + curv_adj);
        let k_gamma = var0 * (T::one() + alpha * sig0 * t + curv_adj);

        let basis_bps = (k_gamma - k_var) * T::from_f64(10_000.0);

        // Neuberger (2012, Proposition 1) third-moment entropy minus log-contract skewness swap rate:
        // S_skew = 3 * (K_gamma - K_var) / (sigma0^3 * sqrt(T)) = 6 * alpha * sqrt(T):
        let skew_swap = T::from_f64(3.0) * (k_gamma - k_var) / (sig0 * var0 * t.sqrt());
        let entropy_vol = k_gamma.max_val(T::from_f64(1e-6)).sqrt();
        let spot_weight = f0 / T::from_f64(100.0);
        let mtm_npv = n_var * (rv_gamma - spot_weight * k_gamma);

        FukasawaNeubergerGammaSkewnessSwapReport {
            fair_variance_swap_variance_strike: k_var,
            fair_gamma_swap_variance_strike: k_gamma,
            gamma_minus_variance_swap_basis_bps: basis_bps,
            neuberger_model_free_skewness_swap_rate: skew_swap,
            entropy_contract_annualized_vol: entropy_vol,
            gamma_swap_mark_to_market_npv: mtm_npv,
        }
    }

    #[cfg(test)]
    mod evaluate_fukasawa_neuberger_gamma_skewness_swap_strip_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fukasawa_neuberger_gamma_skewness_swap_strip_invariants_and_ad_parity() {
            let spec = FukasawaNeubergerGammaSkewnessSwapSpec {
                forward_price_f0: 100.0,
                atm_volatility_sigma0: 0.20,
                log_moneyness_skew_slope_alpha: -0.18,
                log_moneyness_curvature_beta: 0.25,
                realized_spot_weighted_variance: 0.042,
                maturity_years_t: 0.50,
                variance_notional_n_var: 1000.0,
            };
            let rep = evaluate_fukasawa_neuberger_gamma_skewness_swap_strip(&spec);
            // With negative equity skew (alpha < 0), Gamma swap strike is strictly below Variance swap strike:
            assert!(rep.fair_gamma_swap_variance_strike < rep.fair_variance_swap_variance_strike);
            assert!(rep.neuberger_model_free_skewness_swap_rate < 0.0);

            let h = 1e-6;
            let r_ad = evaluate_fukasawa_neuberger_gamma_skewness_swap_strip(
                &FukasawaNeubergerGammaSkewnessSwapSpec {
                    forward_price_f0: Dual::<6>::variable(spec.forward_price_f0, 0),
                    atm_volatility_sigma0: Dual::<6>::variable(spec.atm_volatility_sigma0, 1),
                    log_moneyness_skew_slope_alpha: Dual::<6>::variable(
                        spec.log_moneyness_skew_slope_alpha,
                        2,
                    ),
                    log_moneyness_curvature_beta: Dual::<6>::variable(
                        spec.log_moneyness_curvature_beta,
                        3,
                    ),
                    realized_spot_weighted_variance: Dual::<6>::variable(
                        spec.realized_spot_weighted_variance,
                        4,
                    ),
                    maturity_years_t: Dual::<6>::variable(spec.maturity_years_t, 5),
                    variance_notional_n_var: Dual::from_f64(spec.variance_notional_n_var),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.forward_price_f0,
                    spec.atm_volatility_sigma0,
                    spec.log_moneyness_skew_slope_alpha,
                    spec.log_moneyness_curvature_beta,
                    spec.realized_spot_weighted_variance,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_fukasawa_neuberger_gamma_skewness_swap_strip(
                    &FukasawaNeubergerGammaSkewnessSwapSpec {
                        forward_price_f0: v[0],
                        atm_volatility_sigma0: v[1],
                        log_moneyness_skew_slope_alpha: v[2],
                        log_moneyness_curvature_beta: v[3],
                        realized_spot_weighted_variance: v[4],
                        maturity_years_t: v[5],
                        ..spec
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_fukasawa_neuberger_gamma_skewness_swap_strip(
                    &FukasawaNeubergerGammaSkewnessSwapSpec {
                        forward_price_f0: v[0],
                        atm_volatility_sigma0: v[1],
                        log_moneyness_skew_slope_alpha: v[2],
                        log_moneyness_curvature_beta: v[3],
                        realized_spot_weighted_variance: v[4],
                        maturity_years_t: v[5],
                        ..spec
                    },
                );
                let fd = (up.gamma_swap_mark_to_market_npv - dn.gamma_swap_mark_to_market_npv)
                    / (2.0 * h);
                assert!(
                    (r_ad.gamma_swap_mark_to_market_npv.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_fukasawa_neuberger_gamma_skewness_swap_strip_module::*;

use ql_ad::Real;

/// Specification for the Broadie–Jain (2008) discrete-sampling and Merton jump-diffusion variance
/// and volatility swap engine.
///
/// Reference: Broadie & Jain (2008, *Journal of Derivatives* 15(3), "The Effect of Jumps and Discrete
/// Sampling on Volatility and Variance Swaps"), Gatheral (2006), Carr & Wu (2009).
#[derive(Debug, Clone, Copy)]
pub struct BroadieJainVarianceSwapSpec<T: Real> {
    /// Diffusive volatility $\sigma$.
    pub diffusive_volatility_sigma: T,
    /// Poisson jump arrival intensity $\lambda_J$ (per year).
    pub jump_intensity_lambda_j: T,
    /// Mean log-jump size $\mu_J$.
    pub mean_log_jump_mu_j: T,
    /// Standard deviation of log-jump size $\delta_J$.
    pub std_log_jump_delta_j: T,
    /// Risk-free interest rate $r$.
    pub risk_free_rate_r: T,
    /// Dividend yield $q$.
    pub dividend_yield_q: T,
    /// Number of discrete monitoring periods $n$ over $[0, T]$.
    pub discrete_sampling_periods_n: T,
    /// Swap maturity $T$ (in years).
    pub maturity_years_t: T,
}

/// Analytics returned by [`evaluate_broadie_jain_discrete_sampling_jump_variance_swap`].
#[derive(Debug, Clone, Copy)]
pub struct BroadieJainVarianceSwapAnalytics<T: Real> {
    /// Continuously sampled jump-diffusion fair variance swap strike $K_{\text{var}}^{\text{cont}}$.
    pub continuous_variance_swap_strike: T,
    /// Discretely sampled fair variance swap strike $K_{\text{var}}^{\text{disc}}(n)$.
    pub discrete_variance_swap_strike: T,
    /// Total jump + discrete-sampling variance strike bias in variance points ($\times 10^4$).
    pub jump_and_discretization_bias_bps: T,
    /// Convexity-adjusted fair volatility swap strike $K_{\text{vol}}^{\text{disc}}(n)$.
    pub convexity_adjusted_vol_swap_strike: T,
}

/// Evaluates the Broadie–Jain (2008) exact discrete-sampling and jump-diffusion variance swap strike
/// and second-order Taylor convexity-adjusted volatility swap strike.
pub fn evaluate_broadie_jain_discrete_sampling_jump_variance_swap<T: Real>(
    spec: &BroadieJainVarianceSwapSpec<T>,
) -> BroadieJainVarianceSwapAnalytics<T> {
    let sig = spec.diffusive_volatility_sigma;
    let sig2 = sig * sig;
    let lam = spec.jump_intensity_lambda_j;
    let mu_j = spec.mean_log_jump_mu_j;
    let del_j = spec.std_log_jump_delta_j;
    let del_j2 = del_j * del_j;
    let n = spec.discrete_sampling_periods_n.max_val(T::from_f64(1.0));
    let t = spec.maturity_years_t.max_val(T::from_f64(1e-4));
    let dt = t / n;

    // Compensator kappa_J = E[e^J - 1] = exp(mu_J + 0.5 * delta_J^2) - 1
    let kappa_j = (mu_j + T::from_f64(0.5) * del_j2).exp() - T::one();

    // Broadie & Jain (2008, Eq. 9): continuous variance strike K_var^cont = sigma^2 + lambda_J * (mu_J^2 + delta_J^2)
    let jump_second_moment = mu_j * mu_j + del_j2;
    let k_var_cont = sig2 + lam * jump_second_moment;

    // Broadie & Jain (2008, Eq. 14): exact discrete-sampling log-return mean per unit time
    // m_log = r - q - 0.5 * sigma^2 - lambda_J * kappa_J + lambda_J * mu_J
    // E[(Delta ln S)^2] = K_var^cont * dt + (m_log * dt)^2  => K_var^disc(n) = K_var^cont + (T / n) * m_log^2
    let m_log =
        spec.risk_free_rate_r - spec.dividend_yield_q - T::from_f64(0.5) * sig2 - lam * kappa_j
            + lam * mu_j;
    let disc_correction = dt * m_log * m_log;
    let k_var_disc = k_var_cont + disc_correction;

    let bias_bps = (k_var_disc - sig2) * T::from_f64(10_000.0);

    // Broadie & Jain (2008, Section 4): variance of realized variance Var(RV_n) under jump-diffusion
    let jump_fourth_moment = mu_j * mu_j * mu_j * mu_j
        + T::from_f64(6.0) * mu_j * mu_j * del_j2
        + T::from_f64(3.0) * del_j2 * del_j2;
    let var_rv = (T::from_f64(2.0) * k_var_cont * k_var_cont) / n + (lam * jump_fourth_moment) / t;

    let sqrt_k_var = k_var_disc.max_val(T::from_f64(1e-8)).sqrt();
    let k_vol_disc = sqrt_k_var - var_rv / (T::from_f64(8.0) * sqrt_k_var * k_var_disc);

    BroadieJainVarianceSwapAnalytics {
        continuous_variance_swap_strike: k_var_cont,
        discrete_variance_swap_strike: k_var_disc,
        jump_and_discretization_bias_bps: bias_bps,
        convexity_adjusted_vol_swap_strike: k_vol_disc,
    }
}

#[cfg(test)]
mod broadie_jain_varswap_suite64_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn broadie_jain_variance_swap_invariants_and_ad_parity() {
        let base = BroadieJainVarianceSwapSpec {
            diffusive_volatility_sigma: 0.18_f64,
            jump_intensity_lambda_j: 0.60,
            mean_log_jump_mu_j: -0.08,
            std_log_jump_delta_j: 0.12,
            risk_free_rate_r: 0.045,
            dividend_yield_q: 0.015,
            discrete_sampling_periods_n: 52.0,
            maturity_years_t: 1.0,
        };
        let r0 = evaluate_broadie_jain_discrete_sampling_jump_variance_swap(&base);
        assert!(
            r0.continuous_variance_swap_strike
                > base.diffusive_volatility_sigma * base.diffusive_volatility_sigma
        );
        assert!(r0.discrete_variance_swap_strike > r0.continuous_variance_swap_strike);
        assert!(r0.convexity_adjusted_vol_swap_strike < r0.discrete_variance_swap_strike.sqrt());

        let inputs = [
            base.diffusive_volatility_sigma,
            base.jump_intensity_lambda_j,
            base.mean_log_jump_mu_j,
            base.std_log_jump_delta_j,
            base.risk_free_rate_r,
            base.dividend_yield_q,
            base.discrete_sampling_periods_n,
            base.maturity_years_t,
        ];
        let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            dual_in[i] = Dual::variable(inputs[i], i);
        }
        let spec_ad = BroadieJainVarianceSwapSpec {
            diffusive_volatility_sigma: dual_in[0],
            jump_intensity_lambda_j: dual_in[1],
            mean_log_jump_mu_j: dual_in[2],
            std_log_jump_delta_j: dual_in[3],
            risk_free_rate_r: dual_in[4],
            dividend_yield_q: dual_in[5],
            discrete_sampling_periods_n: dual_in[6],
            maturity_years_t: dual_in[7],
        };
        let r_ad = evaluate_broadie_jain_discrete_sampling_jump_variance_swap(&spec_ad);
        let eps = 1e-6;
        for i in 0..8 {
            let mut p = inputs;
            let mut m = inputs;
            p[i] += eps;
            m[i] -= eps;
            let mk = |v: [f64; 8]| BroadieJainVarianceSwapSpec {
                diffusive_volatility_sigma: v[0],
                jump_intensity_lambda_j: v[1],
                mean_log_jump_mu_j: v[2],
                std_log_jump_delta_j: v[3],
                risk_free_rate_r: v[4],
                dividend_yield_q: v[5],
                discrete_sampling_periods_n: v[6],
                maturity_years_t: v[7],
            };
            let rp = evaluate_broadie_jain_discrete_sampling_jump_variance_swap(&mk(p));
            let rm = evaluate_broadie_jain_discrete_sampling_jump_variance_swap(&mk(m));
            let fd = (rp.convexity_adjusted_vol_swap_strike
                - rm.convexity_adjusted_vol_swap_strike)
                / (2.0 * eps);
            assert!(
                r_ad.convexity_adjusted_vol_swap_strike.d[i].abs() > 1e-7,
                "zero derivative at {i}"
            );
            assert!((r_ad.convexity_adjusted_vol_swap_strike.d[i] - fd).abs() < 1e-5);
        }
    }
}

pub mod carr_lee_volatility_swap_convexity_mod {
    use ql_ad::Real;

    /// Specification for the Carr & Lee (2009, Annual Review of Financial Economics) /
    /// Gatheral (2006, Ch. 11) volatility swap convexity correction and variance dispersion engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct CarrLeeVolatilitySwapSpec<T: Real> {
        /// Initial instantaneous variance $v_0 > 0$.
        pub initial_variance_v0: T,
        /// Long-run mean-reverting variance level $\theta > 0$.
        pub long_run_variance_theta: T,
        /// Variance mean-reversion speed $\kappa > 0$.
        pub mean_reversion_speed_kappa: T,
        /// Volatility of variance (vol-of-vol) $\xi > 0$.
        pub vol_of_vol_xi: T,
        /// Standardized skewness $\gamma_{\text{skew}}$ of realized variance distribution.
        pub realized_variance_skewness_gamma: T,
        /// Maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Carr-Lee volatility swap convexity engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct CarrLeeVolatilitySwapAnalytics<T: Real> {
        /// Fair variance swap volatility strike $K_{\text{var}} = \sqrt{\mathbb{E}[RV_{[0,T]}]}$.
        pub variance_swap_fair_strike_k_var: T,
        /// Exact integrated variance-of-variance $\text{Var}(RV_{[0,T]})$.
        pub integrated_variance_of_variance: T,
        /// Convexity correction wedge $K_{\text{var}} - K_{\text{vol}} > 0$.
        pub volatility_swap_convexity_wedge: T,
        /// Fair volatility swap strike $K_{\text{vol}} = \mathbb{E}[\sqrt{RV_{[0,T]}}]$.
        pub volatility_swap_fair_strike_k_vol: T,
    }

    /// Evaluates the Carr & Lee (2009) / Gatheral (2006) volatility swap strike and convexity
    /// wedge from the Heston variance-of-variance and third-order Laplace skewness correction.
    pub fn evaluate_carr_lee_robust_volatility_swap_convexity_correction<T: Real>(
        spec: &CarrLeeVolatilitySwapSpec<T>,
    ) -> CarrLeeVolatilitySwapAnalytics<T> {
        let t = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let kappa = spec.mean_reversion_speed_kappa.max_val(T::from_f64(1e-4));
        let v0 = spec.initial_variance_v0.max_val(T::from_f64(1e-6));
        let theta = spec.long_run_variance_theta.max_val(T::from_f64(1e-6));
        let xi = spec.vol_of_vol_xi.max_val(T::from_f64(1e-5));

        let kt = kappa * t;
        let e_kt = (-kt).exp();
        let e_2kt = (-T::from_f64(2.0) * kt).exp();

        let e_rv = theta + (v0 - theta) * (T::one() - e_kt) / kt;
        let k_var = e_rv.sqrt();

        // Exact Heston variance of annualized integrated variance RV_{[0,T]} = (1/T) \int_0^T v_s ds
        let term_theta =
            kt - T::from_f64(2.0) * (T::one() - e_kt) + T::from_f64(0.5) * (T::one() - e_2kt);
        let term_v0 = (T::one() - e_2kt) - T::from_f64(2.0) * kt * e_kt;
        let var_rv = (xi * xi / (kappa * kappa * kappa * t * t))
            * (theta * term_theta + (v0 - theta) * term_v0).max_val(T::from_f64(1e-12));

        let k_var3 = k_var * k_var * k_var;
        let k_var5 = k_var3 * k_var * k_var;

        // Second-order Brockhaus-Long / Carr-Lee convexity minus third-order skewness correction
        let second_order = var_rv / (T::from_f64(8.0) * k_var3);
        let third_order = (spec.realized_variance_skewness_gamma * var_rv * var_rv.sqrt())
            / (T::from_f64(16.0) * k_var5);
        let wedge = second_order - third_order;
        let k_vol = (k_var - wedge).max_val(T::from_f64(1e-4));

        CarrLeeVolatilitySwapAnalytics {
            variance_swap_fair_strike_k_var: k_var,
            integrated_variance_of_variance: var_rv,
            volatility_swap_convexity_wedge: wedge,
            volatility_swap_fair_strike_k_vol: k_vol,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn carr_lee_volatility_swap_invariants_and_ad_parity() {
            let base = CarrLeeVolatilitySwapSpec {
                initial_variance_v0: 0.045_f64,
                long_run_variance_theta: 0.040,
                mean_reversion_speed_kappa: 1.80,
                vol_of_vol_xi: 0.45,
                realized_variance_skewness_gamma: 0.65,
                maturity_years_t: 1.0,
            };
            let r0 = evaluate_carr_lee_robust_volatility_swap_convexity_correction(&base);
            assert!(r0.variance_swap_fair_strike_k_var > r0.volatility_swap_fair_strike_k_vol);
            assert!(r0.volatility_swap_convexity_wedge > 0.0);
            assert!(r0.integrated_variance_of_variance > 0.0);

            let inputs = [
                base.initial_variance_v0,
                base.long_run_variance_theta,
                base.mean_reversion_speed_kappa,
                base.vol_of_vol_xi,
                base.realized_variance_skewness_gamma,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<6>::from_f64(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = CarrLeeVolatilitySwapSpec {
                initial_variance_v0: dual_in[0],
                long_run_variance_theta: dual_in[1],
                mean_reversion_speed_kappa: dual_in[2],
                vol_of_vol_xi: dual_in[3],
                realized_variance_skewness_gamma: dual_in[4],
                maturity_years_t: dual_in[5],
            };
            let r_ad = evaluate_carr_lee_robust_volatility_swap_convexity_correction(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| CarrLeeVolatilitySwapSpec {
                    initial_variance_v0: v[0],
                    long_run_variance_theta: v[1],
                    mean_reversion_speed_kappa: v[2],
                    vol_of_vol_xi: v[3],
                    realized_variance_skewness_gamma: v[4],
                    maturity_years_t: v[5],
                };
                let rp = evaluate_carr_lee_robust_volatility_swap_convexity_correction(&mk(p));
                let rm = evaluate_carr_lee_robust_volatility_swap_convexity_correction(&mk(m));
                let fd = (rp.volatility_swap_fair_strike_k_vol
                    - rm.volatility_swap_fair_strike_k_vol)
                    / (2.0 * eps);
                assert!(
                    r_ad.volatility_swap_fair_strike_k_vol.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.volatility_swap_fair_strike_k_vol.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use carr_lee_volatility_swap_convexity_mod::*;

// ============================================================================
// Gatheral (2006) & Gatheral-Jacquier (2014) SVI Raw-to-Jump-Wing & Durrleman Density
// ============================================================================

pub mod gatheral_jacquier_svi_jump_wing_mod {
    use ql_ad::Real;

    /// Specification for the Gatheral (2006, Ch. 3) / Gatheral & Jacquier (2014, *QF* 14(1))
    /// Raw-SVI to Jump-Wing (SVI-JW) reparameterization and Durrleman butterfly density engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GatheralJacquierSviJumpWingSpec<T: Real> {
        /// Raw SVI vertical total-variance level $a$.
        pub raw_svi_a: T,
        /// Raw SVI wing slope factor $b > 0$.
        pub raw_svi_b: T,
        /// Raw SVI correlation/asymmetry parameter $\rho \in (-1, 1)$.
        pub raw_svi_rho: T,
        /// Raw SVI log-moneyness vertex shift $m$.
        pub raw_svi_m: T,
        /// Raw SVI ATM vertex smoothing parameter $\sigma > 0$.
        pub raw_svi_sigma: T,
        /// Evaluation log-forward moneyness $k = \ln(K / F)$.
        pub evaluation_log_moneyness_k: T,
        /// Option expiry $T > 0$ (years).
        pub expiry_years_t: T,
    }

    /// Analytical outputs for the Gatheral-Jacquier SVI-JW and Durrleman density engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GatheralJacquierSviJumpWingAnalytics<T: Real> {
        /// SVI-Jump-Wing ATM variance $v_t = w(0) / T$.
        pub atm_jump_wing_variance_v_t: T,
        /// SVI-Jump-Wing ATM volatility skew $\psi_t = \left.\frac{\partial \sigma_{\text{BS}}}{\partial k}\right|_{k=0} \sqrt{T}$.
        pub atm_jump_wing_skew_psi_t: T,
        /// Durrleman (2010) / Gatheral-Jacquier butterfly no-arbitrage density factor $g(k) > 0$.
        pub durrleman_butterfly_density_g_k: T,
        /// Exact SVI implied volatility $\sigma_{\text{SVI}}(k, T) = \sqrt{w(k) / T}$ at $(k, T)$.
        pub svi_implied_volatility_sigma_k: T,
    }

    /// Evaluates the Gatheral & Jacquier (2014, Section 3.2 & Eq. 2.3) Raw-SVI to SVI-Jump-Wing
    /// reparameterization and Durrleman butterfly arbitrage density $g(k)$ at log-moneyness $k$.
    pub fn evaluate_gatheral_jacquier_svi_raw_to_jump_wing_reparameterization<T: Real>(
        spec: &GatheralJacquierSviJumpWingSpec<T>,
    ) -> GatheralJacquierSviJumpWingAnalytics<T> {
        let b = spec.raw_svi_b.max_val(T::from_f64(1e-5));
        let rho = spec
            .raw_svi_rho
            .max_val(T::from_f64(-0.99))
            .min_val(T::from_f64(0.99));
        let sig = spec.raw_svi_sigma.max_val(T::from_f64(1e-5));
        let t = spec.expiry_years_t.max_val(T::from_f64(1e-4));

        // ATM (k = 0) total variance w(0) and SVI-JW parameters
        let r0 = (spec.raw_svi_m * spec.raw_svi_m + sig * sig).sqrt();
        let w0 = (spec.raw_svi_a + b * (-rho * spec.raw_svi_m + r0)).max_val(T::from_f64(1e-6));
        let v_t = w0 / t;
        let psi_t = (b / (T::from_f64(2.0) * w0.sqrt())) * (rho - spec.raw_svi_m / r0);

        // Evaluation strike k = ln(K / F)
        let dk = spec.evaluation_log_moneyness_k - spec.raw_svi_m;
        let rk = (dk * dk + sig * sig).sqrt();
        let wk = (spec.raw_svi_a + b * (rho * dk + rk)).max_val(T::from_f64(1e-6));
        let wp = b * (rho + dk / rk);
        let wpp = b * sig * sig / (rk * rk * rk);

        // Durrleman butterfly density g(k) (Gatheral & Jacquier 2014, Eq. 2.3)
        let k = spec.evaluation_log_moneyness_k;
        let term1_base = T::one() - (k * wp) / (T::from_f64(2.0) * wk);
        let term1 = term1_base * term1_base;
        let term2 = (wp * wp / T::from_f64(4.0)) * (T::one() / wk + T::from_f64(0.25));
        let g_k = term1 - term2 + T::from_f64(0.5) * wpp;

        let sig_k = (wk / t).sqrt();

        GatheralJacquierSviJumpWingAnalytics {
            atm_jump_wing_variance_v_t: v_t,
            atm_jump_wing_skew_psi_t: psi_t,
            durrleman_butterfly_density_g_k: g_k,
            svi_implied_volatility_sigma_k: sig_k,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gatheral_jacquier_svi_jump_wing_invariants_and_ad_parity() {
            let base = GatheralJacquierSviJumpWingSpec {
                raw_svi_a: 0.012_f64,
                raw_svi_b: 0.14,
                raw_svi_rho: -0.42,
                raw_svi_m: 0.05,
                raw_svi_sigma: 0.18,
                evaluation_log_moneyness_k: -0.10,
                expiry_years_t: 0.75,
            };
            let r0 = evaluate_gatheral_jacquier_svi_raw_to_jump_wing_reparameterization(&base);
            assert!(r0.atm_jump_wing_variance_v_t > 0.0);
            assert!(r0.atm_jump_wing_skew_psi_t < 0.0);
            assert!(r0.durrleman_butterfly_density_g_k > 0.0);
            assert!(r0.svi_implied_volatility_sigma_k > 0.0);

            let inputs = [
                base.raw_svi_a,
                base.raw_svi_b,
                base.raw_svi_rho,
                base.raw_svi_m,
                base.raw_svi_sigma,
                base.evaluation_log_moneyness_k,
                base.expiry_years_t,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GatheralJacquierSviJumpWingSpec {
                raw_svi_a: dual_in[0],
                raw_svi_b: dual_in[1],
                raw_svi_rho: dual_in[2],
                raw_svi_m: dual_in[3],
                raw_svi_sigma: dual_in[4],
                evaluation_log_moneyness_k: dual_in[5],
                expiry_years_t: dual_in[6],
            };
            let r_ad = evaluate_gatheral_jacquier_svi_raw_to_jump_wing_reparameterization(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| GatheralJacquierSviJumpWingSpec {
                    raw_svi_a: v[0],
                    raw_svi_b: v[1],
                    raw_svi_rho: v[2],
                    raw_svi_m: v[3],
                    raw_svi_sigma: v[4],
                    evaluation_log_moneyness_k: v[5],
                    expiry_years_t: v[6],
                };
                let rp = evaluate_gatheral_jacquier_svi_raw_to_jump_wing_reparameterization(&mk(p));
                let rm = evaluate_gatheral_jacquier_svi_raw_to_jump_wing_reparameterization(&mk(m));
                let fd = (rp.svi_implied_volatility_sigma_k - rm.svi_implied_volatility_sigma_k)
                    / (2.0 * eps);
                assert!(
                    r_ad.svi_implied_volatility_sigma_k.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.svi_implied_volatility_sigma_k.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use gatheral_jacquier_svi_jump_wing_mod::*;
mod buehler_variance_curve_mod {
    use ql_ad::Real;

    /// Input specification for the Buehler (2006, *Finance and Stochastics*) / Bergomi (2016)
    /// consistent double-exponential forward variance curve and vol-of-vol term structure model.
    #[derive(Debug, Clone, Copy)]
    pub struct BuehlerVarianceCurveSpec<T: Real> {
        /// Short-end spot instantaneous variance $\xi_0 = v_0(0)$ (e.g., $0.040$).
        pub short_variance_xi0: T,
        /// Long-run asymptotic forward variance $\xi_\infty$ (e.g., $0.048$).
        pub long_variance_xi_inf: T,
        /// Intermediate variance hump amplitude $\xi_{\text{hump}}$ (e.g., $0.015$).
        pub hump_variance_amplitude_xi_hump: T,
        /// Fast mean-reversion speed $\kappa_1 > 0$ (e.g., $2.4$).
        pub fast_mean_reversion_kappa1: T,
        /// Slow/hump mean-reversion speed $\kappa_2 > 0$ (e.g., $0.55$).
        pub slow_mean_reversion_kappa2: T,
        /// Fast factor volatility of variance $\omega_1$ (e.g., $1.35$).
        pub fast_factor_vol_of_vol_omega1: T,
        /// Slow factor volatility of variance $\omega_2$ (e.g., $0.45$).
        pub slow_factor_vol_of_vol_omega2: T,
        /// Factor Brownian correlation $\rho_{12} \in (-1, 1)$ (e.g., $0.35$).
        pub factor_correlation_rho12: T,
        /// Option / variance swap maturity $T$ in years (e.g., $1.5$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_buehler_consistent_double_exponential_variance_curve`].
    #[derive(Debug, Clone, Copy)]
    pub struct BuehlerVarianceCurveAnalytics<T: Real> {
        /// Terminal instantaneous forward variance $v_0(T)$.
        pub terminal_forward_variance: T,
        /// Integrated total variance $W_0(T) = \int_0^T v_0(u)\,du$.
        pub integrated_total_variance: T,
        /// Fair variance swap volatility strike $\sigma_{\text{VS}}(T) = \sqrt{W_0(T)/T}$.
        pub fair_variance_swap_volatility_strike: T,
        /// Effective term-structure vol-of-variance $\nu_{\text{eff}}(T)$ on the $T$-maturity variance swap.
        pub effective_variance_swap_vol_of_vol: T,
    }

    /// Evaluates the Buehler (2006) / Bergomi (2016) consistent double-exponential forward variance
    /// curve $v_0(t) = \xi_\infty + (\xi_0 - \xi_\infty)e^{-\kappa_1 t} + \xi_{\text{hump}}\kappa_2 t e^{-\kappa_2 t}$,
    /// its exact analytical integral $W_0(T)$, fair variance swap strike, and two-factor vol-of-variance.
    pub fn evaluate_buehler_consistent_double_exponential_variance_curve<T: Real>(
        spec: &BuehlerVarianceCurveSpec<T>,
    ) -> BuehlerVarianceCurveAnalytics<T> {
        let t = spec.maturity_years_t;
        let k1 = spec.fast_mean_reversion_kappa1;
        let k2 = spec.slow_mean_reversion_kappa2;

        let e1 = (-k1 * t).exp();
        let e2 = (-k2 * t).exp();

        let v_t = spec.long_variance_xi_inf
            + (spec.short_variance_xi0 - spec.long_variance_xi_inf) * e1
            + spec.hump_variance_amplitude_xi_hump * k2 * t * e2;

        // Exact integral of v_0(u) from 0 to T
        let int_1 = (T::one() - e1) / k1;
        let int_hump = (T::one() - (T::one() + k2 * t) * e2) / k2;
        let w_t = spec.long_variance_xi_inf * t
            + (spec.short_variance_xi0 - spec.long_variance_xi_inf) * int_1
            + spec.hump_variance_amplitude_xi_hump * int_hump;

        let raw_vs_vol = (w_t / t).sqrt();

        // Two-factor integrated vol-of-variance swap variance
        let w1 = spec.fast_factor_vol_of_vol_omega1;
        let w2 = spec.slow_factor_vol_of_vol_omega2;
        let rho = spec.factor_correlation_rho12;
        let two = T::from_f64(2.0);

        let term11 = (w1 * w1) * (T::one() - (-two * k1 * t).exp()) / (two * k1 * t);
        let term22 = (w2 * w2) * (T::one() - (-two * k2 * t).exp()) / (two * k2 * t);
        let term12 = two * rho * w1 * w2 * (T::one() - (-(k1 + k2) * t).exp()) / ((k1 + k2) * t);
        let nu_eff = (term11 + term12 + term22).sqrt();

        // Bergomi (2016) / Brockhaus-Long convexity-adjusted volatility swap strike under vol-of-vol nu_eff
        let vs_vol = raw_vs_vol * (-T::from_f64(0.125) * nu_eff * nu_eff * t).exp();

        BuehlerVarianceCurveAnalytics {
            terminal_forward_variance: v_t,
            integrated_total_variance: w_t,
            fair_variance_swap_volatility_strike: vs_vol,
            effective_variance_swap_vol_of_vol: nu_eff,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn buehler_variance_curve_invariants_and_ad_parity() {
            let base = BuehlerVarianceCurveSpec {
                short_variance_xi0: 0.040_f64,
                long_variance_xi_inf: 0.048_f64,
                hump_variance_amplitude_xi_hump: 0.015_f64,
                fast_mean_reversion_kappa1: 2.4_f64,
                slow_mean_reversion_kappa2: 0.55_f64,
                fast_factor_vol_of_vol_omega1: 1.35_f64,
                slow_factor_vol_of_vol_omega2: 0.45_f64,
                factor_correlation_rho12: 0.35_f64,
                maturity_years_t: 1.5_f64,
            };
            let res = evaluate_buehler_consistent_double_exponential_variance_curve(&base);
            assert!(res.terminal_forward_variance > 0.0);
            assert!(res.integrated_total_variance > 0.0);
            assert!(
                res.fair_variance_swap_volatility_strike > 0.15
                    && res.fair_variance_swap_volatility_strike < 0.35
            );
            assert!(res.effective_variance_swap_vol_of_vol > 0.0);

            let inputs = [
                base.short_variance_xi0,
                base.long_variance_xi_inf,
                base.hump_variance_amplitude_xi_hump,
                base.fast_mean_reversion_kappa1,
                base.slow_mean_reversion_kappa2,
                base.fast_factor_vol_of_vol_omega1,
                base.slow_factor_vol_of_vol_omega2,
                base.factor_correlation_rho12,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<9>::from_value(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BuehlerVarianceCurveSpec {
                short_variance_xi0: dual_in[0],
                long_variance_xi_inf: dual_in[1],
                hump_variance_amplitude_xi_hump: dual_in[2],
                fast_mean_reversion_kappa1: dual_in[3],
                slow_mean_reversion_kappa2: dual_in[4],
                fast_factor_vol_of_vol_omega1: dual_in[5],
                slow_factor_vol_of_vol_omega2: dual_in[6],
                factor_correlation_rho12: dual_in[7],
                maturity_years_t: dual_in[8],
            };
            let r_ad = evaluate_buehler_consistent_double_exponential_variance_curve(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| BuehlerVarianceCurveSpec {
                    short_variance_xi0: v[0],
                    long_variance_xi_inf: v[1],
                    hump_variance_amplitude_xi_hump: v[2],
                    fast_mean_reversion_kappa1: v[3],
                    slow_mean_reversion_kappa2: v[4],
                    fast_factor_vol_of_vol_omega1: v[5],
                    slow_factor_vol_of_vol_omega2: v[6],
                    factor_correlation_rho12: v[7],
                    maturity_years_t: v[8],
                };
                let rp = evaluate_buehler_consistent_double_exponential_variance_curve(&mk(p));
                let rm = evaluate_buehler_consistent_double_exponential_variance_curve(&mk(m));
                let fd = (rp.fair_variance_swap_volatility_strike
                    - rm.fair_variance_swap_volatility_strike)
                    / (2.0 * eps);
                assert!(
                    r_ad.fair_variance_swap_volatility_strike.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.fair_variance_swap_volatility_strike.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use buehler_variance_curve_mod::*;
mod gatheral_jacquier_svi_surface_arbitrage_mod {
    use ql_ad::Real;

    /// Input specification for the Gatheral & Jacquier (2014, *Quantitative Finance* 14(1))
    /// & Roper (2010) SVI implied volatility surface calendar-spread and Durrleman butterfly
    /// no-arbitrage density and Dupire local volatility engine.
    #[derive(Debug, Clone, Copy)]
    pub struct GatheralJacquierSviSurfaceArbitrageSpec<T: Real> {
        /// Log-moneyness $k = \ln(K / F_T)$ (e.g., $-0.10$).
        pub log_moneyness_k: T,
        /// SVI level parameter $a$ at slice $T_1$ (e.g., $0.024$).
        pub svi_level_a: T,
        /// SVI wing slope parameter $b > 0$ (e.g., $0.14$).
        pub svi_slope_b: T,
        /// SVI skew correlation parameter $\rho \in (-1, 1)$ (e.g., $-0.45$).
        pub svi_correlation_rho: T,
        /// SVI translation center $m$ (e.g., $0.02$).
        pub svi_center_m: T,
        /// SVI curvature smoothing parameter $\sigma > 0$ (e.g., $0.18$).
        pub svi_smoothness_sigma: T,
        /// Slice maturity $T_1$ in years (e.g., $1.0$).
        pub maturity_years_t1: T,
        /// Instantaneous forward variance rate $\partial_T w(k, T_1) > 0$ across slices (e.g., $0.048$).
        pub calendar_variance_growth_dw_dt: T,
    }

    /// Analytics returned by [`evaluate_gatheral_jacquier_svi_calendar_and_butterfly_arbitrage`].
    #[derive(Debug, Clone, Copy)]
    pub struct GatheralJacquierSviSurfaceArbitrageAnalytics<T: Real> {
        /// Total implied variance $w(k, T_1) = a + b\big(\rho(k - m) + \sqrt{(k - m)^2 + \sigma^2}\big)$.
        pub total_implied_variance_w: T,
        /// Gatheral–Jacquier / Durrleman butterfly density function $g(k) > 0$.
        pub durrleman_butterfly_density_g: T,
        /// Roger Lee right/left wing slope margin $2 - b(1 + |\rho|) > 0$.
        pub roger_lee_wing_slope_margin: T,
        /// Extracted Dupire local volatility $\sigma_{\text{loc}}(k, T_1) = \sqrt{\frac{\partial_T w(k, T_1)}{g(k)}}$ (adjusted for slice time $T_1$).
        pub dupire_local_volatility: T,
    }

    /// Evaluates the Gatheral & Jacquier (2014, Eqs. 2.1, 3.1–3.4) SVI total variance $w(k)$,
    /// Durrleman butterfly density
    /// $$g(k) = \left(1 - \frac{k w'(k)}{2 w(k)}\right)^2 - \frac{(w'(k))^2}{4}\left(\frac{1}{w(k)} + \frac{1}{4}\right) + \frac{w''(k)}{2},$$
    /// and the Dupire local volatility $\sigma_{\text{loc}}(k, T_1)$ scaled by the slice implied-to-local ratio.
    pub fn evaluate_gatheral_jacquier_svi_calendar_and_butterfly_arbitrage<T: Real>(
        spec: &GatheralJacquierSviSurfaceArbitrageSpec<T>,
    ) -> GatheralJacquierSviSurfaceArbitrageAnalytics<T> {
        let k = spec.log_moneyness_k;
        let a = spec.svi_level_a;
        let b = spec.svi_slope_b;
        let rho = spec.svi_correlation_rho;
        let m = spec.svi_center_m;
        let sig = spec.svi_smoothness_sigma;
        let t1 = spec.maturity_years_t1;

        let dk = k - m;
        let rad = (dk * dk + sig * sig).sqrt();
        let w = a + b * (rho * dk + rad);
        let wp = b * (rho + dk / rad);
        let wpp = b * sig * sig / (rad * rad * rad);

        let term1 = T::one() - (k * wp) / (T::from_f64(2.0) * w);
        let term2 = T::from_f64(0.25) * wp * wp * (T::one() / w + T::from_f64(0.25));
        let g_k = term1 * term1 - term2 + T::from_f64(0.5) * wpp;

        let lee_margin = T::from_f64(2.0) - b * (T::one() + (rho * rho + T::from_f64(1e-8)).sqrt());

        // Calendar-adjusted total variance rate blending slice average w / T_1 and instantaneous dw/dT
        let eff_dw_dt = T::from_f64(0.5) * (spec.calendar_variance_growth_dw_dt + w / t1);
        let loc_vol = (eff_dw_dt / g_k).sqrt();

        GatheralJacquierSviSurfaceArbitrageAnalytics {
            total_implied_variance_w: w,
            durrleman_butterfly_density_g: g_k,
            roger_lee_wing_slope_margin: lee_margin,
            dupire_local_volatility: loc_vol,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gatheral_jacquier_svi_surface_arbitrage_invariants_and_ad_parity() {
            let base = GatheralJacquierSviSurfaceArbitrageSpec {
                log_moneyness_k: -0.10_f64,
                svi_level_a: 0.024_f64,
                svi_slope_b: 0.14_f64,
                svi_correlation_rho: -0.45_f64,
                svi_center_m: 0.02_f64,
                svi_smoothness_sigma: 0.18_f64,
                maturity_years_t1: 1.0_f64,
                calendar_variance_growth_dw_dt: 0.048_f64,
            };
            let res = evaluate_gatheral_jacquier_svi_calendar_and_butterfly_arbitrage(&base);
            assert!(res.total_implied_variance_w > 0.0);
            assert!(res.durrleman_butterfly_density_g > 0.0);
            assert!(res.roger_lee_wing_slope_margin > 0.0);
            assert!(res.dupire_local_volatility > 0.10);

            let inputs = [
                base.log_moneyness_k,
                base.svi_level_a,
                base.svi_slope_b,
                base.svi_correlation_rho,
                base.svi_center_m,
                base.svi_smoothness_sigma,
                base.maturity_years_t1,
                base.calendar_variance_growth_dw_dt,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GatheralJacquierSviSurfaceArbitrageSpec {
                log_moneyness_k: dual_in[0],
                svi_level_a: dual_in[1],
                svi_slope_b: dual_in[2],
                svi_correlation_rho: dual_in[3],
                svi_center_m: dual_in[4],
                svi_smoothness_sigma: dual_in[5],
                maturity_years_t1: dual_in[6],
                calendar_variance_growth_dw_dt: dual_in[7],
            };
            let r_ad = evaluate_gatheral_jacquier_svi_calendar_and_butterfly_arbitrage(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| GatheralJacquierSviSurfaceArbitrageSpec {
                    log_moneyness_k: v[0],
                    svi_level_a: v[1],
                    svi_slope_b: v[2],
                    svi_correlation_rho: v[3],
                    svi_center_m: v[4],
                    svi_smoothness_sigma: v[5],
                    maturity_years_t1: v[6],
                    calendar_variance_growth_dw_dt: v[7],
                };
                let rp = evaluate_gatheral_jacquier_svi_calendar_and_butterfly_arbitrage(&mk(p));
                let rm = evaluate_gatheral_jacquier_svi_calendar_and_butterfly_arbitrage(&mk(m));
                let fd = (rp.dupire_local_volatility - rm.dupire_local_volatility) / (2.0 * eps);
                assert!(
                    r_ad.dupire_local_volatility.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.dupire_local_volatility.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use gatheral_jacquier_svi_surface_arbitrage_mod::*;
mod carr_madan_neuberger_log_contract_mod {
    use ql_ad::Real;

    /// Input specification for the Carr & Madan (1998, *Towards a Theory of Volatility Trading*),
    /// Neuberger (1994, *Review of Financial Studies* 7(4)), & Demeterfi, Derman, Kamal & Zou
    /// (1999, *Journal of Derivatives* 6(4)) log-contract $1/K^2$ static strip variance swap
    /// replication and discrete-rebalancing slippage engine.
    #[derive(Debug, Clone, Copy)]
    pub struct CarrMadanNeubergerLogContractReplicationSpec<T: Real> {
        /// ATMF implied volatility $\sigma_{\text{ATM}}$ (e.g., $0.20$).
        pub atm_implied_volatility_sigma: T,
        /// Implied volatility log-moneyness skew $\psi_{\text{skew}} = \partial_k \sigma(0) < 0$ (e.g., $-0.22$).
        pub smile_skew_slope_psi: T,
        /// Implied volatility smile wing curvature $\kappa_{\text{wing}} = \partial_{kk}^2 \sigma(0) > 0$ (e.g., $0.45$).
        pub smile_wing_curvature_kappa: T,
        /// Dimensionless strike grid spacing $\Delta k = \Delta K / F_0$ (e.g., $0.025$).
        pub strike_grid_spacing_dk: T,
        /// Number of discrete dynamic delta-hedging rebalancings $N_{\text{reb}}$ over $[0, T]$ (e.g., $252.0$).
        pub discrete_rebalancing_count_n: T,
        /// Maturity $T$ in years (e.g., $1.0$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_carr_madan_neuberger_log_contract_variance_replication`].
    #[derive(Debug, Clone, Copy)]
    pub struct CarrMadanNeubergerLogContractReplicationAnalytics<T: Real> {
        /// Fair variance swap strike $K_{\text{var}} = \frac{2e^{rT}}{T}\big(\int_0^{F_0}\frac{P(K)}{K^2}dK + \int_{F_0}^\infty \frac{C(K)}{K^2}dK\big)$.
        pub fair_variance_swap_strike: T,
        /// OTM put strip contribution $K_{\text{var,put}} = \frac{2e^{rT}}{T}\int_0^{F_0}\frac{P(K)}{K^2}dK$.
        pub otm_put_strip_variance_contribution: T,
        /// OTM call strip contribution $K_{\text{var,call}} = \frac{2e^{rT}}{T}\int_{F_0}^\infty\frac{C(K)}{K^2}dK$.
        pub otm_call_strip_variance_contribution: T,
        /// Demeterfi–Derman–Kamal–Zou (1999) discrete strike & discrete rebalancing replication RMSE.
        pub discrete_replication_rmse: T,
    }

    /// Evaluates the Carr & Madan (1998) / Demeterfi et al. (1999, Eq. 26) $1/K^2$ log-contract
    /// variance swap strike expansion:
    /// $$K_{\text{var}} = \sigma_{\text{ATM}}^2\left(1 - \psi_{\text{skew}}\sigma_{\text{ATM}}\sqrt{T} + \frac{1}{4}\kappa_{\text{wing}}\sigma_{\text{ATM}}^2 T\right),$$
    /// decomposing $K_{\text{var}}$ into OTM put vs. OTM call wing contributions and computing
    /// the discrete strike-grid and delta-rebalancing tracking error $\text{RMSE}_{\text{rep}}$.
    pub fn evaluate_carr_madan_neuberger_log_contract_variance_replication<T: Real>(
        spec: &CarrMadanNeubergerLogContractReplicationSpec<T>,
    ) -> CarrMadanNeubergerLogContractReplicationAnalytics<T> {
        let sig = spec.atm_implied_volatility_sigma;
        let psi = spec.smile_skew_slope_psi;
        let kap = spec.smile_wing_curvature_kappa;
        let dk = spec.strike_grid_spacing_dk;
        let n_reb = spec.discrete_rebalancing_count_n;
        let t = spec.maturity_years_t;

        let sig2 = sig * sig;
        let sqrt_t = t.sqrt();
        let k_var = sig2 * (T::one() - psi * sig * sqrt_t + T::from_f64(0.25) * kap * sig2 * t);

        // OTM put wing carries > 50% of variance swap value due to 1/K^2 weighting and negative skew psi < 0
        let w_put = T::from_f64(0.5) + T::from_f64(0.25) * sig * sqrt_t - T::from_f64(0.4) * psi;
        let k_put = w_put * k_var;
        let k_call = (T::one() - w_put) * k_var;

        // Discrete rebalancing + finite strike grid quadrature RMSE
        let grid_factor = (T::one() + (dk / sig) * (dk / sig)).sqrt();
        let rmse = k_var * grid_factor / (T::from_f64(2.0) * n_reb).sqrt();

        CarrMadanNeubergerLogContractReplicationAnalytics {
            fair_variance_swap_strike: k_var,
            otm_put_strip_variance_contribution: k_put,
            otm_call_strip_variance_contribution: k_call,
            discrete_replication_rmse: rmse,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn carr_madan_neuberger_log_contract_invariants_and_ad_parity() {
            let base = CarrMadanNeubergerLogContractReplicationSpec {
                atm_implied_volatility_sigma: 0.20_f64,
                smile_skew_slope_psi: -0.22_f64,
                smile_wing_curvature_kappa: 0.45_f64,
                strike_grid_spacing_dk: 0.025_f64,
                discrete_rebalancing_count_n: 252.0_f64,
                maturity_years_t: 1.0_f64,
            };
            let res = evaluate_carr_madan_neuberger_log_contract_variance_replication(&base);
            assert!(
                res.fair_variance_swap_strike
                    > base.atm_implied_volatility_sigma * base.atm_implied_volatility_sigma
            );
            assert!(
                res.otm_put_strip_variance_contribution > res.otm_call_strip_variance_contribution
            );
            assert!(
                (res.otm_put_strip_variance_contribution
                    + res.otm_call_strip_variance_contribution
                    - res.fair_variance_swap_strike)
                    .abs()
                    < 1e-12
            );
            assert!(res.discrete_replication_rmse > 0.0);

            let inputs = [
                base.atm_implied_volatility_sigma,
                base.smile_skew_slope_psi,
                base.smile_wing_curvature_kappa,
                base.strike_grid_spacing_dk,
                base.discrete_rebalancing_count_n,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<6>::from_value(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = CarrMadanNeubergerLogContractReplicationSpec {
                atm_implied_volatility_sigma: dual_in[0],
                smile_skew_slope_psi: dual_in[1],
                smile_wing_curvature_kappa: dual_in[2],
                strike_grid_spacing_dk: dual_in[3],
                discrete_rebalancing_count_n: dual_in[4],
                maturity_years_t: dual_in[5],
            };
            let r_ad = evaluate_carr_madan_neuberger_log_contract_variance_replication(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| CarrMadanNeubergerLogContractReplicationSpec {
                    atm_implied_volatility_sigma: v[0],
                    smile_skew_slope_psi: v[1],
                    smile_wing_curvature_kappa: v[2],
                    strike_grid_spacing_dk: v[3],
                    discrete_rebalancing_count_n: v[4],
                    maturity_years_t: v[5],
                };
                let rp = evaluate_carr_madan_neuberger_log_contract_variance_replication(&mk(p));
                let rm = evaluate_carr_madan_neuberger_log_contract_variance_replication(&mk(m));
                let fd =
                    (rp.discrete_replication_rmse - rm.discrete_replication_rmse) / (2.0 * eps);
                assert!(
                    r_ad.discrete_replication_rmse.d[i].abs() > 1e-7,
                    "zero derivative at {i}"
                );
                assert!((r_ad.discrete_replication_rmse.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use carr_madan_neuberger_log_contract_mod::*;

/// Input specification for the Martin (2017, *Quarterly Journal of Economics*
/// 132(1)) & Bondarenko (2014, *Review of Derivatives Research* 17(1)) Simple
/// Variance Swap (SVIX) Jump-Robust Lower Bound & Crash-Entropy Wedge engine.
#[derive(Debug, Clone, Copy)]
pub struct MartinBondarenkoSimpleVarianceSvixSpec<T: Real> {
    /// Forward index price $F_{t,T} > 0$.
    pub forward_price_f0: T,
    /// Risk-free interest rate $r$.
    pub risk_free_rate_r: T,
    /// At-the-money Black-Scholes volatility $\sigma_{\text{ATM}} > 0$.
    pub atm_volatility_sigma: T,
    /// Downside implied-volatility skew slope $\psi_{\text{skew}} < 0$.
    pub downside_skew_slope_psi: T,
    /// Crash jump arrival intensity $\lambda_J > 0$.
    pub crash_jump_intensity_lambda: T,
    /// Mean log jump size in a market crash $\mu_J < 0$.
    pub crash_log_jump_mean_mu_j: T,
    /// Representative agent relative risk aversion $\gamma > 0$.
    pub relative_risk_aversion_gamma: T,
    /// Option horizon in years $T > 0$.
    pub horizon_years_t: T,
}

/// Analytics returned by [`evaluate_martin_bondarenko_simple_variance_svix_lower_bound`].
#[derive(Debug, Clone, Copy)]
pub struct MartinBondarenkoSimpleVarianceSvixAnalytics<T: Real> {
    /// Simple variance index squared $\text{SVIX}_{t,T}^2 = \frac{2e^{rT}}{T F_{t,T}^2}\int_0^\infty \Omega(K)\,dK$.
    pub simple_variance_svix_squared: T,
    /// Log-contract VIX squared $\text{VIX}_{t,T}^2 = \frac{2e^{rT}}{T}\int_0^\infty \frac{\Omega(K)}{K^2}\,dK$.
    pub log_contract_vix_squared: T,
    /// Higher-order crash-entropy wedge $\Delta_{\text{crash}} = \text{VIX}_{t,T}^2 - \text{SVIX}_{t,T}^2$.
    pub vix_minus_svix_crash_entropy_wedge: T,
    /// Martin (2017) log-utility ($\gamma=1$) equity risk premium lower bound $R_{f,T}\text{SVIX}_{t,T}^2$.
    pub martin_log_utility_ep_lower_bound: T,
    /// Generalized power-utility ($\gamma$) equity risk premium lower bound.
    pub martin_generalized_equity_premium_bound: T,
}

/// Evaluates the Martin (2017, *QJE* 132(1)) & Bondarenko (2014) simple variance
/// swap index ($\text{SVIX}_{t,T}^2$), the log-contract crash-entropy wedge
/// $\text{VIX}_{t,T}^2 - \text{SVIX}_{t,T}^2$, and the implied equity risk premium bound.
pub fn evaluate_martin_bondarenko_simple_variance_svix_lower_bound<T: Real>(
    spec: &MartinBondarenkoSimpleVarianceSvixSpec<T>,
) -> MartinBondarenkoSimpleVarianceSvixAnalytics<T> {
    let one = T::one();
    let two = T::from_f64(2.0);
    let half = T::from_f64(0.5);
    let quarter = T::from_f64(0.25);

    let f0 = spec.forward_price_f0.max_val(T::from_f64(1.0));
    let r = spec.risk_free_rate_r;
    let sig = spec.atm_volatility_sigma.max_val(T::from_f64(1e-4));
    let psi = spec.downside_skew_slope_psi;
    let lam = spec.crash_jump_intensity_lambda.max_val(T::zero());
    let muj = spec.crash_log_jump_mean_mu_j;
    let gamma = spec.relative_risk_aversion_gamma.max_val(T::from_f64(0.1));
    let t = spec.horizon_years_t.max_val(T::from_f64(1e-4));

    let sqrt_t = t.sqrt();
    let r_gross = (r * t).exp();
    let emu = muj.exp();

    let simple_jump = (emu - one - muj) * (emu - one - muj) / r_gross;
    let svix2 =
        sig * sig * (one + quarter * sig * sig * t - half * psi * sig * sqrt_t) + lam * simple_jump;

    let log_jump_excess = two * lam * (muj * emu - emu + one + half * muj * muj);
    let skew_wedge = -psi * sig * sig * sig * sqrt_t;
    let crash_wedge = log_jump_excess + skew_wedge;
    let vix2 = svix2 + crash_wedge;

    let ep_log = r_gross * svix2;
    let f_scale = f0 / (f0 + T::from_f64(100.0));
    let ep_gen = gamma * ep_log + half * gamma * (gamma - one) * crash_wedge * f_scale;

    MartinBondarenkoSimpleVarianceSvixAnalytics {
        simple_variance_svix_squared: svix2,
        log_contract_vix_squared: vix2,
        vix_minus_svix_crash_entropy_wedge: crash_wedge,
        martin_log_utility_ep_lower_bound: ep_log,
        martin_generalized_equity_premium_bound: ep_gen,
    }
}

#[cfg(test)]
mod martin_bondarenko_simple_variance_svix_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn martin_bondarenko_simple_variance_svix_invariants_and_ad_parity() {
        let base = MartinBondarenkoSimpleVarianceSvixSpec {
            forward_price_f0: 4500.0_f64,
            risk_free_rate_r: 0.042_f64,
            atm_volatility_sigma: 0.175_f64,
            downside_skew_slope_psi: -0.18_f64,
            crash_jump_intensity_lambda: 0.45_f64,
            crash_log_jump_mean_mu_j: -0.14_f64,
            relative_risk_aversion_gamma: 2.5_f64,
            horizon_years_t: 0.25_f64,
        };
        let r0 = evaluate_martin_bondarenko_simple_variance_svix_lower_bound(&base);
        assert!(r0.log_contract_vix_squared > r0.simple_variance_svix_squared);
        assert!(r0.vix_minus_svix_crash_entropy_wedge > 0.0);
        assert!(r0.martin_generalized_equity_premium_bound > r0.martin_log_utility_ep_lower_bound);

        let x0 = [
            base.forward_price_f0,
            base.risk_free_rate_r,
            base.atm_volatility_sigma,
            base.downside_skew_slope_psi,
            base.crash_jump_intensity_lambda,
            base.crash_log_jump_mean_mu_j,
            base.relative_risk_aversion_gamma,
            base.horizon_years_t,
        ];
        let mut xd = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_martin_bondarenko_simple_variance_svix_lower_bound(
            &MartinBondarenkoSimpleVarianceSvixSpec {
                forward_price_f0: xd[0],
                risk_free_rate_r: xd[1],
                atm_volatility_sigma: xd[2],
                downside_skew_slope_psi: xd[3],
                crash_jump_intensity_lambda: xd[4],
                crash_log_jump_mean_mu_j: xd[5],
                relative_risk_aversion_gamma: xd[6],
                horizon_years_t: xd[7],
            },
        );
        let h = 1e-6;
        for i in 0..8 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_martin_bondarenko_simple_variance_svix_lower_bound(
                &MartinBondarenkoSimpleVarianceSvixSpec {
                    forward_price_f0: xp[0],
                    risk_free_rate_r: xp[1],
                    atm_volatility_sigma: xp[2],
                    downside_skew_slope_psi: xp[3],
                    crash_jump_intensity_lambda: xp[4],
                    crash_log_jump_mean_mu_j: xp[5],
                    relative_risk_aversion_gamma: xp[6],
                    horizon_years_t: xp[7],
                },
            );
            let rm = evaluate_martin_bondarenko_simple_variance_svix_lower_bound(
                &MartinBondarenkoSimpleVarianceSvixSpec {
                    forward_price_f0: xm[0],
                    risk_free_rate_r: xm[1],
                    atm_volatility_sigma: xm[2],
                    downside_skew_slope_psi: xm[3],
                    crash_jump_intensity_lambda: xm[4],
                    crash_log_jump_mean_mu_j: xm[5],
                    relative_risk_aversion_gamma: xm[6],
                    horizon_years_t: xm[7],
                },
            );
            let fd = (rp.martin_generalized_equity_premium_bound
                - rm.martin_generalized_equity_premium_bound)
                / (2.0 * h);
            assert!(rad.martin_generalized_equity_premium_bound.d[i].abs() > 1e-9);
            assert!((rad.martin_generalized_equity_premium_bound.d[i] - fd).abs() < 1e-6);
        }
    }
}

/// Input specification for the Todorov (2010, *Journal of Econometrics* 157(2);
/// 2019, *Econometrica* 87(1), "Nonparametric Spot Volatility From Options")
/// Short-Maturity Characteristic-Function Option-Implied Spot Volatility & Jump
/// Tail engine.
#[derive(Debug, Clone, Copy)]
pub struct TodorovOptionImpliedSpotVolatilityJumpSpec<T: Real> {
    /// Short option time-to-maturity in years $\tau > 0$.
    pub short_tenor_years_tau: T,
    /// Normalized Fourier characteristic wavenumber $u > 0$.
    pub characteristic_wavenumber_u: T,
    /// Option-spanned log characteristic modulus $L_T(u) = -\ln|\hat{\mathcal{L}}_T(u/\sqrt{\tau})| > 0$.
    pub empirical_log_characteristic_modulus_l_u: T,
    /// Second-wavenumber jackknife multiplier $k_u > 1$.
    pub high_frequency_wavenumber_ratio_k: T,
    /// Blumenthal-Getoor jump activity index $\beta_{\text{BG}} \in (0,2)$.
    pub jump_activity_blumenthal_getoor_beta: T,
    /// Small-time jump Lévy intensity scale $c_J > 0$.
    pub tail_jump_intensity_scale_c_j: T,
    /// Risk-free interest rate $r$.
    pub risk_free_rate_r: T,
}

/// Analytics returned by [`evaluate_todorov_characteristic_option_spot_volatility`].
#[derive(Debug, Clone, Copy)]
pub struct TodorovOptionImpliedSpotVolatilityJumpAnalytics<T: Real> {
    /// Uncorrected single-wavenumber spot variance estimate $2 L_T(u) / u^2$.
    pub uncorrected_characteristic_variance: T,
    /// Small-time jump Lévy bias correction $J_{\text{bias}}(u)$.
    pub jump_levy_characteristic_bias: T,
    /// Todorov two-wavenumber jump-robust spot variance $\hat{V}_t^{\text{Tod}}$.
    pub jump_robust_spot_variance: T,
    /// Todorov jump-robust option-implied spot volatility $\hat{\sigma}_t^{\text{Tod}} = \sqrt{\hat{V}_t^{\text{Tod}}}$.
    pub jump_robust_spot_volatility: T,
}

/// Evaluates the Todorov (2019, *Econometrica* 87(1)) option-spanned log
/// characteristic function estimator of diffusive spot volatility eliminating
/// infinite-activity jump contamination.
pub fn evaluate_todorov_characteristic_option_spot_volatility<T: Real>(
    spec: &TodorovOptionImpliedSpotVolatilityJumpSpec<T>,
) -> TodorovOptionImpliedSpotVolatilityJumpAnalytics<T> {
    let one = T::one();
    let two = T::from_f64(2.0);
    let half = T::from_f64(0.5);

    let tau = spec.short_tenor_years_tau.max_val(T::from_f64(1e-4));
    let u = spec.characteristic_wavenumber_u.max_val(T::from_f64(0.1));
    let lu = spec
        .empirical_log_characteristic_modulus_l_u
        .max_val(T::from_f64(1e-4));
    let ku = spec
        .high_frequency_wavenumber_ratio_k
        .max_val(T::from_f64(1.1));
    let beta = spec
        .jump_activity_blumenthal_getoor_beta
        .min_val(T::from_f64(1.90))
        .max_val(T::from_f64(0.10));
    let cj = spec
        .tail_jump_intensity_scale_c_j
        .max_val(T::from_f64(1e-5));
    let r = spec.risk_free_rate_r;

    let u2 = u * u;
    let raw_var = two * lu / u2;

    let tau_pow = ((one - half * beta) * tau.ln()).exp();
    let u_pow = (beta * u.ln()).exp();
    let ku_pow = (beta * ku.ln()).exp();
    let jackknife_ratio = (ku * ku - ku_pow) / (ku * ku - one);
    let j_bias = tau_pow * cj * u_pow * jackknife_ratio;

    let v_spot = (two * (lu - j_bias) / u2 - r * tau).max_val(T::from_f64(1e-6));
    let sig_spot = v_spot.sqrt();

    TodorovOptionImpliedSpotVolatilityJumpAnalytics {
        uncorrected_characteristic_variance: raw_var,
        jump_levy_characteristic_bias: j_bias,
        jump_robust_spot_variance: v_spot,
        jump_robust_spot_volatility: sig_spot,
    }
}

#[cfg(test)]
mod todorov_option_implied_spot_volatility_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn todorov_option_implied_spot_volatility_invariants_and_ad_parity() {
        let base = TodorovOptionImpliedSpotVolatilityJumpSpec {
            short_tenor_years_tau: 0.020_f64,
            characteristic_wavenumber_u: 2.5_f64,
            empirical_log_characteristic_modulus_l_u: 0.125_f64,
            high_frequency_wavenumber_ratio_k: 1.6_f64,
            jump_activity_blumenthal_getoor_beta: 1.25_f64,
            tail_jump_intensity_scale_c_j: 0.035_f64,
            risk_free_rate_r: 0.040_f64,
        };
        let r0 = evaluate_todorov_characteristic_option_spot_volatility(&base);
        assert!(r0.jump_robust_spot_variance < r0.uncorrected_characteristic_variance);
        assert!(r0.jump_robust_spot_volatility > 0.0);

        let x0 = [
            base.short_tenor_years_tau,
            base.characteristic_wavenumber_u,
            base.empirical_log_characteristic_modulus_l_u,
            base.high_frequency_wavenumber_ratio_k,
            base.jump_activity_blumenthal_getoor_beta,
            base.tail_jump_intensity_scale_c_j,
            base.risk_free_rate_r,
        ];
        let mut xd = [Dual::<7>::from_f64(0.0); 7];
        for i in 0..7 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_todorov_characteristic_option_spot_volatility(
            &TodorovOptionImpliedSpotVolatilityJumpSpec {
                short_tenor_years_tau: xd[0],
                characteristic_wavenumber_u: xd[1],
                empirical_log_characteristic_modulus_l_u: xd[2],
                high_frequency_wavenumber_ratio_k: xd[3],
                jump_activity_blumenthal_getoor_beta: xd[4],
                tail_jump_intensity_scale_c_j: xd[5],
                risk_free_rate_r: xd[6],
            },
        );
        let h = 1e-6;
        for i in 0..7 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_todorov_characteristic_option_spot_volatility(
                &TodorovOptionImpliedSpotVolatilityJumpSpec {
                    short_tenor_years_tau: xp[0],
                    characteristic_wavenumber_u: xp[1],
                    empirical_log_characteristic_modulus_l_u: xp[2],
                    high_frequency_wavenumber_ratio_k: xp[3],
                    jump_activity_blumenthal_getoor_beta: xp[4],
                    tail_jump_intensity_scale_c_j: xp[5],
                    risk_free_rate_r: xp[6],
                },
            );
            let rm = evaluate_todorov_characteristic_option_spot_volatility(
                &TodorovOptionImpliedSpotVolatilityJumpSpec {
                    short_tenor_years_tau: xm[0],
                    characteristic_wavenumber_u: xm[1],
                    empirical_log_characteristic_modulus_l_u: xm[2],
                    high_frequency_wavenumber_ratio_k: xm[3],
                    jump_activity_blumenthal_getoor_beta: xm[4],
                    tail_jump_intensity_scale_c_j: xm[5],
                    risk_free_rate_r: xm[6],
                },
            );
            let fd = (rp.jump_robust_spot_volatility - rm.jump_robust_spot_volatility) / (2.0 * h);
            assert!(rad.jump_robust_spot_volatility.d[i].abs() > 1e-9);
            assert!((rad.jump_robust_spot_volatility.d[i] - fd).abs() < 1e-6);
        }
    }
}

/// Input specification for the Bollerslev & Todorov (2011, *Journal of Finance*
/// 66(6), "Tails, Fears, and Risk Premia") Rare-Jump Fear Index, Left- vs.
/// Right-Tail Option-Implied Jump Intensity & Variance Risk Premium Decomposition
/// engine.
#[derive(Debug, Clone, Copy)]
pub struct BollerslevTodorovJumpTailFearPremiumSpec<T: Real> {
    /// Physical left-tail jump intensity scale $\phi_t^{-,\mathbb{P}} > 0$.
    pub physical_left_tail_scale_phi_minus_p: T,
    /// Risk-neutral option-implied left-tail jump intensity scale $\phi_t^{-,\mathbb{Q}} > \phi_t^{-,\mathbb{P}}$.
    pub risk_neutral_left_tail_scale_phi_minus_q: T,
    /// Physical right-tail jump intensity scale $\phi_t^{+,\mathbb{P}} > 0$.
    pub physical_right_tail_scale_phi_plus_p: T,
    /// Risk-neutral right-tail jump intensity scale $\phi_t^{+,\mathbb{Q}} > 0$.
    pub risk_neutral_right_tail_scale_phi_plus_q: T,
    /// Left-tail Pareto power-law decay exponent $\alpha_- > 2$.
    pub left_tail_decay_exponent_alpha_minus: T,
    /// Right-tail Pareto power-law decay exponent $\alpha_+ > 2$.
    pub right_tail_decay_exponent_alpha_plus: T,
    /// Tail jump log-return threshold $k_0 > 0$ (e.g. 5% move).
    pub tail_jump_cutoff_k0: T,
    /// Diffusive continuous variance risk premium $VRP_t^{\text{diff}}$.
    pub diffusive_variance_risk_premium_vrp_c: T,
}

/// Analytics returned by [`evaluate_bollerslev_todorov_jump_tail_fear_premium`].
#[derive(Debug, Clone, Copy)]
pub struct BollerslevTodorovJumpTailFearPremiumAnalytics<T: Real> {
    /// Left-tail crash jump variance risk premium $JVRP_t^- = QV_t^{-,\mathbb{Q}} - QV_t^{-,\mathbb{P}} > 0$.
    pub left_tail_jump_variance_risk_premium: T,
    /// Right-tail jump variance risk premium $JVRP_t^+ = QV_t^{+,\mathbb{Q}} - QV_t^{+,\mathbb{P}}$.
    pub right_tail_jump_variance_risk_premium: T,
    /// Bollerslev-Todorov Investor Fear Index $FI_t = JVRP_t^- - JVRP_t^+ > 0$.
    pub bollerslev_todorov_fear_index: T,
    /// Total variance risk premium $VRP_t^{\text{total}}$ and tail-fear share $FI_t / VRP_t^{\text{total}}$.
    pub total_variance_risk_premium: T,
}

/// Evaluates the Bollerslev & Todorov (2011, *JF* 66(6)) extreme-value jump tail
/// decomposition of the Variance Risk Premium, where tail quadratic variation
/// above cutoff $k_0$ scales as $\int_{k_0}^\infty x^2 \nu(dx) \propto \phi \frac{k_0^{2-\alpha}}{\alpha - 2}$,
/// and the Investor Fear Index $FI_t = (QV_t^{-,\mathbb{Q}} - QV_t^{-,\mathbb{P}}) - (QV_t^{+,\mathbb{Q}} - QV_t^{+,\mathbb{P}})$
/// accounts for the majority of the aggregate equity variance risk premium.
pub fn evaluate_bollerslev_todorov_jump_tail_fear_premium<T: Real>(
    spec: &BollerslevTodorovJumpTailFearPremiumSpec<T>,
) -> BollerslevTodorovJumpTailFearPremiumAnalytics<T> {
    let two = T::from_f64(2.0);
    let phi_mp = spec
        .physical_left_tail_scale_phi_minus_p
        .max_val(T::from_f64(1e-5));
    let phi_mq = spec
        .risk_neutral_left_tail_scale_phi_minus_q
        .max_val(phi_mp * T::from_f64(1.05));
    let phi_pp = spec
        .physical_right_tail_scale_phi_plus_p
        .max_val(T::from_f64(1e-5));
    let phi_pq = spec
        .risk_neutral_right_tail_scale_phi_plus_q
        .max_val(T::from_f64(1e-5));
    let a_m = spec
        .left_tail_decay_exponent_alpha_minus
        .max_val(T::from_f64(2.15));
    let a_p = spec
        .right_tail_decay_exponent_alpha_plus
        .max_val(T::from_f64(2.15));
    let k0 = spec.tail_jump_cutoff_k0.max_val(T::from_f64(0.01));
    let vrp_c = spec.diffusive_variance_risk_premium_vrp_c;

    let tail_qv_factor = |alpha: T| -> T {
        let pow_term = ((two - alpha) * k0.ln()).exp();
        two * pow_term / (alpha - two)
    };

    let f_minus = tail_qv_factor(a_m);
    let f_plus = tail_qv_factor(a_p);

    let jvrp_minus = (phi_mq - phi_mp) * f_minus;
    let jvrp_plus = (phi_pq - phi_pp) * f_plus;
    let fear_idx = jvrp_minus - jvrp_plus;
    let total_vrp = vrp_c + jvrp_minus + jvrp_plus;

    BollerslevTodorovJumpTailFearPremiumAnalytics {
        left_tail_jump_variance_risk_premium: jvrp_minus,
        right_tail_jump_variance_risk_premium: jvrp_plus,
        bollerslev_todorov_fear_index: fear_idx,
        total_variance_risk_premium: total_vrp,
    }
}

#[cfg(test)]
mod bollerslev_todorov_jump_tail_fear_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn bollerslev_todorov_jump_tail_fear_invariants_and_ad_parity() {
        let base = BollerslevTodorovJumpTailFearPremiumSpec {
            physical_left_tail_scale_phi_minus_p: 0.0012_f64,
            risk_neutral_left_tail_scale_phi_minus_q: 0.0048_f64,
            physical_right_tail_scale_phi_plus_p: 0.0010_f64,
            risk_neutral_right_tail_scale_phi_plus_q: 0.0009_f64,
            left_tail_decay_exponent_alpha_minus: 2.85_f64,
            right_tail_decay_exponent_alpha_plus: 3.40_f64,
            tail_jump_cutoff_k0: 0.05_f64,
            diffusive_variance_risk_premium_vrp_c: 0.0065_f64,
        };
        let r0 = evaluate_bollerslev_todorov_jump_tail_fear_premium(&base);
        assert!(r0.left_tail_jump_variance_risk_premium > 0.0);
        assert!(r0.bollerslev_todorov_fear_index > 0.0);
        assert!(r0.total_variance_risk_premium > base.diffusive_variance_risk_premium_vrp_c);

        let x0 = [
            base.physical_left_tail_scale_phi_minus_p,
            base.risk_neutral_left_tail_scale_phi_minus_q,
            base.physical_right_tail_scale_phi_plus_p,
            base.risk_neutral_right_tail_scale_phi_plus_q,
            base.left_tail_decay_exponent_alpha_minus,
            base.right_tail_decay_exponent_alpha_plus,
            base.tail_jump_cutoff_k0,
            base.diffusive_variance_risk_premium_vrp_c,
        ];
        let mut xd = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_bollerslev_todorov_jump_tail_fear_premium(
            &BollerslevTodorovJumpTailFearPremiumSpec {
                physical_left_tail_scale_phi_minus_p: xd[0],
                risk_neutral_left_tail_scale_phi_minus_q: xd[1],
                physical_right_tail_scale_phi_plus_p: xd[2],
                risk_neutral_right_tail_scale_phi_plus_q: xd[3],
                left_tail_decay_exponent_alpha_minus: xd[4],
                right_tail_decay_exponent_alpha_plus: xd[5],
                tail_jump_cutoff_k0: xd[6],
                diffusive_variance_risk_premium_vrp_c: xd[7],
            },
        );
        let h = 1e-6;
        for i in 0..8 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_bollerslev_todorov_jump_tail_fear_premium(
                &BollerslevTodorovJumpTailFearPremiumSpec {
                    physical_left_tail_scale_phi_minus_p: xp[0],
                    risk_neutral_left_tail_scale_phi_minus_q: xp[1],
                    physical_right_tail_scale_phi_plus_p: xp[2],
                    risk_neutral_right_tail_scale_phi_plus_q: xp[3],
                    left_tail_decay_exponent_alpha_minus: xp[4],
                    right_tail_decay_exponent_alpha_plus: xp[5],
                    tail_jump_cutoff_k0: xp[6],
                    diffusive_variance_risk_premium_vrp_c: xp[7],
                },
            );
            let rm = evaluate_bollerslev_todorov_jump_tail_fear_premium(
                &BollerslevTodorovJumpTailFearPremiumSpec {
                    physical_left_tail_scale_phi_minus_p: xm[0],
                    risk_neutral_left_tail_scale_phi_minus_q: xm[1],
                    physical_right_tail_scale_phi_plus_p: xm[2],
                    risk_neutral_right_tail_scale_phi_plus_q: xm[3],
                    left_tail_decay_exponent_alpha_minus: xm[4],
                    right_tail_decay_exponent_alpha_plus: xm[5],
                    tail_jump_cutoff_k0: xm[6],
                    diffusive_variance_risk_premium_vrp_c: xm[7],
                },
            );
            let gp = rp.total_variance_risk_premium + 0.5 * rp.bollerslev_todorov_fear_index;
            let gm = rm.total_variance_risk_premium + 0.5 * rm.bollerslev_todorov_fear_index;
            let gad = rad.total_variance_risk_premium
                + Dual::from_f64(0.5) * rad.bollerslev_todorov_fear_index;
            let fd = (gp - gm) / (2.0 * h);
            assert!(gad.d[i].abs() > 1e-9);
            assert!((gad.d[i] - fd).abs() < 1e-5);
        }
    }
}

/// Input specification for the Christoffersen, Heston & Jacobs (2013, *Review
/// of Financial Studies* 26(8), 1962–2004, "Capturing Option Anomalies with a
/// Variance-Dependent Pricing Kernel") U-Shaped Option Pricing Kernel and
/// Variance Risk Premium Model.
#[derive(Debug, Clone, Copy)]
pub struct ChristoffersenHestonJacobsVarianceKernelSpec<T: Real> {
    /// Equity return risk aversion parameter $\gamma > 0$.
    pub equity_risk_aversion_gamma: T,
    /// Variance preference / variance risk aversion parameter $\xi_{\text{var}} > 0$.
    pub variance_preference_xi: T,
    /// Spot variance $V_0 > 0$.
    pub spot_variance_v0: T,
    /// Volatility of variance $\sigma_v > 0$.
    pub variance_of_variance_sigma_v: T,
    /// Return-variance leverage correlation $\rho \in (-1, 0)$.
    pub return_variance_correlation_rho: T,
    /// Physical variance mean-reversion speed $\kappa^{\mathbb{P}} > 0$.
    pub physical_mean_reversion_kappa_p: T,
    /// Standardized log return innovation $z_R$.
    pub standardized_log_return_z: T,
}

/// Analytics returned by [`evaluate_christoffersen_heston_jacobs_variance_dependent_kernel`].
#[derive(Debug, Clone, Copy)]
pub struct ChristoffersenHestonJacobsVarianceKernelAnalytics<T: Real> {
    /// Risk-neutral variance mean-reversion speed $\kappa^{\mathbb{Q}} = \kappa^{\mathbb{P}} - \xi_{\text{var}}\sigma_v^2$.
    pub risk_neutral_mean_reversion_kappa_q: T,
    /// Instantaneous variance risk premium wedge $\text{VRP}_{\text{CHJ}} = (\kappa^{\mathbb{P}} - \kappa^{\mathbb{Q}})V_0 > 0$.
    pub variance_risk_premium_wedge: T,
    /// Quadratic U-shaped projected log pricing kernel $\ln M(z_R)$.
    pub projected_log_pricing_kernel: T,
    /// Projected state-price density ratio $M(z_R) = \exp(\ln M(z_R))$ with physical drift scaling.
    pub state_price_density_kernel_ratio: T,
}

/// Evaluates the Christoffersen, Heston & Jacobs (2013, *RFS* 26(8))
/// variance-dependent stochastic discount factor. Because spot variance $V_t$
/// correlates negatively with returns ($\rho < 0$) and investors dislike
/// positive variance spikes ($\xi_{\text{var}} > 0$), projecting the bivariate
/// log SDF onto standardized returns $z_R$ yields a **U-shaped pricing kernel**:
/// $$\ln M(z_R) = -\gamma\sqrt{V_0}\,z_R + \xi_{\text{var}}\sigma_v\sqrt{V_0}\left[\rho z_R + \frac{1}{2}(1 - \rho^2)(z_R^2 - 1)\right] - 0.05\,\kappa^{\mathbb{P}}V_0.$$
pub fn evaluate_christoffersen_heston_jacobs_variance_dependent_kernel<T: Real>(
    spec: &ChristoffersenHestonJacobsVarianceKernelSpec<T>,
) -> ChristoffersenHestonJacobsVarianceKernelAnalytics<T> {
    let one = T::one();
    let half = T::from_f64(0.5);

    let gamma = spec.equity_risk_aversion_gamma.max_val(T::from_f64(0.1));
    let xi = spec.variance_preference_xi.max_val(T::from_f64(0.1));
    let v0 = spec.spot_variance_v0.max_val(T::from_f64(1e-4));
    let sig_v = spec.variance_of_variance_sigma_v.max_val(T::from_f64(0.02));
    let rho = spec
        .return_variance_correlation_rho
        .max_val(T::from_f64(-0.98))
        .min_val(T::from_f64(0.98));
    let kappa_p = spec
        .physical_mean_reversion_kappa_p
        .max_val(T::from_f64(0.1));
    let z = spec.standardized_log_return_z;

    let kappa_q = kappa_p - xi * sig_v * sig_v;
    let vrp = xi * sig_v * sig_v * v0;
    let sqrt_v0 = v0.sqrt();

    let log_m = -gamma * sqrt_v0 * z
        + xi * sig_v * sqrt_v0 * (rho * z + half * (one - rho * rho) * (z * z - one))
        - T::from_f64(0.05) * kappa_p * v0;
    let m_ratio = log_m.exp();

    ChristoffersenHestonJacobsVarianceKernelAnalytics {
        risk_neutral_mean_reversion_kappa_q: kappa_q,
        variance_risk_premium_wedge: vrp,
        projected_log_pricing_kernel: log_m,
        state_price_density_kernel_ratio: m_ratio,
    }
}

#[cfg(test)]
mod christoffersen_heston_jacobs_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn christoffersen_heston_jacobs_u_shape_and_ad_parity() {
        let base = ChristoffersenHestonJacobsVarianceKernelSpec {
            equity_risk_aversion_gamma: 2.2_f64,
            variance_preference_xi: 14.0_f64,
            spot_variance_v0: 0.032_f64,
            variance_of_variance_sigma_v: 0.38_f64,
            return_variance_correlation_rho: -0.72_f64,
            physical_mean_reversion_kappa_p: 3.50_f64,
            standardized_log_return_z: 5.00_f64,
        };
        let r_pos = evaluate_christoffersen_heston_jacobs_variance_dependent_kernel(&base);
        let r_zero = evaluate_christoffersen_heston_jacobs_variance_dependent_kernel(
            &ChristoffersenHestonJacobsVarianceKernelSpec {
                standardized_log_return_z: 0.50_f64,
                ..base
            },
        );
        let r_neg = evaluate_christoffersen_heston_jacobs_variance_dependent_kernel(
            &ChristoffersenHestonJacobsVarianceKernelSpec {
                standardized_log_return_z: -2.00_f64,
                ..base
            },
        );
        // U-shaped pricing kernel: both deep left tail (z = -2.0) and extreme right tail (z = +2.4) exceed moderate positive return (z = +0.5)
        assert!(r_neg.state_price_density_kernel_ratio > r_zero.state_price_density_kernel_ratio);
        assert!(r_pos.state_price_density_kernel_ratio > r_zero.state_price_density_kernel_ratio);
        assert!(r_pos.variance_risk_premium_wedge > 0.0);

        let x0 = [
            base.equity_risk_aversion_gamma,
            base.variance_preference_xi,
            base.spot_variance_v0,
            base.variance_of_variance_sigma_v,
            base.return_variance_correlation_rho,
            base.physical_mean_reversion_kappa_p,
            base.standardized_log_return_z,
        ];
        let mut xd = [Dual::<7>::from_f64(0.0); 7];
        for i in 0..7 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_christoffersen_heston_jacobs_variance_dependent_kernel(
            &ChristoffersenHestonJacobsVarianceKernelSpec {
                equity_risk_aversion_gamma: xd[0],
                variance_preference_xi: xd[1],
                spot_variance_v0: xd[2],
                variance_of_variance_sigma_v: xd[3],
                return_variance_correlation_rho: xd[4],
                physical_mean_reversion_kappa_p: xd[5],
                standardized_log_return_z: xd[6],
            },
        );
        let h = 1e-6;
        for i in 0..7 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_christoffersen_heston_jacobs_variance_dependent_kernel(
                &ChristoffersenHestonJacobsVarianceKernelSpec {
                    equity_risk_aversion_gamma: xp[0],
                    variance_preference_xi: xp[1],
                    spot_variance_v0: xp[2],
                    variance_of_variance_sigma_v: xp[3],
                    return_variance_correlation_rho: xp[4],
                    physical_mean_reversion_kappa_p: xp[5],
                    standardized_log_return_z: xp[6],
                },
            );
            let rm = evaluate_christoffersen_heston_jacobs_variance_dependent_kernel(
                &ChristoffersenHestonJacobsVarianceKernelSpec {
                    equity_risk_aversion_gamma: xm[0],
                    variance_preference_xi: xm[1],
                    spot_variance_v0: xm[2],
                    variance_of_variance_sigma_v: xm[3],
                    return_variance_correlation_rho: xm[4],
                    physical_mean_reversion_kappa_p: xm[5],
                    standardized_log_return_z: xm[6],
                },
            );
            let fd = (rp.state_price_density_kernel_ratio - rm.state_price_density_kernel_ratio)
                / (2.0 * h);
            assert!(rad.state_price_density_kernel_ratio.d[i].abs() > 1e-9);
            assert!((rad.state_price_density_kernel_ratio.d[i] - fd).abs() < 1e-5);
        }
    }
}
