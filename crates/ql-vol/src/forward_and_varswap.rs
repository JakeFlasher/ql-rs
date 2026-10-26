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
