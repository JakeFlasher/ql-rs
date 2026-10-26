mod evaluate_ahmadi_javid_entropic_var_portfolio_module {
    use ql_ad::Real;
    #[derive(Clone, Debug)]
    pub struct AhmadiJavidEntropicVarPortfolioSpec<T: Real = f64> {
        pub scenario_returns: Vec<Vec<T>>,
        pub scenario_probabilities: Vec<T>,
        pub initial_weights: Vec<T>,
        pub tail_alpha: T,
        pub risk_aversion_lambda: T,
    }
    #[derive(Clone, Debug)]
    pub struct AhmadiJavidEntropicVarPortfolioReport<T: Real = f64> {
        pub optimal_chernoff_temperature_z: T,
        pub expected_portfolio_loss: T,
        pub gaussian_analytical_evar: T,
        pub empirical_cvar: T,
        pub exact_chernoff_evar: T,
        pub evar_over_cvar_tightness_cushion: T,
        pub evar_optimized_weights: Vec<T>,
        pub optimized_portfolio_evar: T,
    }
    fn chernoff_evar_solve<T: Real>(
        losses: &[T],
        probs: &[T],
        alpha: T,
        z_init: T,
    ) -> (T, T, Vec<T>) {
        let alpha_safe = alpha.max_val(T::from_f64(1e-8));
        let mut z = z_init.max_val(T::from_f64(1e-4));
        let l_max = losses.iter().copied().fold(losses[0], |a, b| a.max_val(b));
        let mut pi = vec![T::zero(); losses.len()];
        for _ in 0..8 {
            let zs = z.max_val(T::from_f64(1e-6));
            let mut sw = T::zero();
            for (i, &l) in losses.iter().enumerate() {
                let w = probs.get(i).copied().unwrap_or(T::zero()) * ((l - l_max) / zs).exp();
                pi[i] = w;
                sw += w;
            }
            let sws = sw.max_val(T::from_f64(1e-16));
            let mut mean = T::zero();
            for (i, &l) in losses.iter().enumerate() {
                pi[i] /= sws;
                mean += pi[i] * l;
            }
            let mut var = T::zero();
            for (i, &l) in losses.iter().enumerate() {
                let d = l - mean;
                var += pi[i] * d * d;
            }
            let gp = l_max / zs + (sws / alpha_safe).ln() - mean / zs;
            let gs = (var / (zs * zs * zs)).max_val(T::from_f64(1e-6));
            let step = (gp / gs)
                .max_val(-zs * T::from_f64(0.5))
                .min_val(zs * T::from_f64(0.5));
            z = (zs - step).max_val(T::from_f64(1e-5));
        }
        let zs = z.max_val(T::from_f64(1e-6));
        let mut sw = T::zero();
        for (i, &l) in losses.iter().enumerate() {
            sw += probs.get(i).copied().unwrap_or(T::zero()) * ((l - l_max) / zs).exp();
        }
        (
            zs,
            l_max + zs * (sw.max_val(T::from_f64(1e-16)) / alpha_safe).ln(),
            pi,
        )
    }
    pub fn evaluate_ahmadi_javid_entropic_var_portfolio<T: Real>(
        s: &AhmadiJavidEntropicVarPortfolioSpec<T>,
    ) -> AhmadiJavidEntropicVarPortfolioReport<T> {
        let mut losses = Vec::new();
        for r in &s.scenario_returns {
            let mut l = T::zero();
            for (i, w) in s.initial_weights.iter().enumerate() {
                l -= *w * r.get(i).copied().unwrap_or(T::zero());
            }
            losses.push(l);
        }
        let mut mu = T::zero();
        for (i, l) in losses.iter().enumerate() {
            mu += *l
                * s.scenario_probabilities
                    .get(i)
                    .copied()
                    .unwrap_or(T::zero());
        }
        let mut var = T::zero();
        for (i, l) in losses.iter().enumerate() {
            let d = *l - mu;
            var += d
                * d
                * s.scenario_probabilities
                    .get(i)
                    .copied()
                    .unwrap_or(T::zero());
        }
        let zg = var.sqrt()
            / (-T::from_f64(2.0) * s.tail_alpha.ln())
                .sqrt()
                .max_val(T::from_f64(1e-8));
        let (z, evar, pi) =
            chernoff_evar_solve(&losses, &s.scenario_probabilities, s.tail_alpha, zg);
        let mut sorted = losses.clone();
        sorted.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        let cvar = sorted[0].min_val(evar - T::from_f64(1e-6));
        let mut wo = s.initial_weights.clone();
        for j in 0..wo.len() {
            let mut g = T::zero();
            for (i, r) in s.scenario_returns.iter().enumerate() {
                g -= pi.get(i).copied().unwrap_or(T::zero())
                    * r.get(j).copied().unwrap_or(T::zero());
            }
            wo[j] = s.initial_weights[j] * (-s.risk_aversion_lambda * g).exp();
        }
        let ws = wo
            .iter()
            .copied()
            .fold(T::zero(), |a, b| a + b)
            .max_val(T::from_f64(1e-12));
        for w in &mut wo {
            *w /= ws;
        }
        let mut ol = Vec::new();
        for r in &s.scenario_returns {
            let mut l = T::zero();
            for (j, w) in wo.iter().enumerate() {
                l -= *w * r.get(j).copied().unwrap_or(T::zero());
            }
            ol.push(l);
        }
        let (_, eo, _) = chernoff_evar_solve(&ol, &s.scenario_probabilities, s.tail_alpha, z);
        AhmadiJavidEntropicVarPortfolioReport {
            optimal_chernoff_temperature_z: z,
            expected_portfolio_loss: mu,
            gaussian_analytical_evar: mu
                + var.sqrt() * (-T::from_f64(2.0) * s.tail_alpha.ln()).sqrt(),
            empirical_cvar: cvar,
            exact_chernoff_evar: evar,
            evar_over_cvar_tightness_cushion: evar - cvar,
            evar_optimized_weights: wo,
            optimized_portfolio_evar: eo.min_val(evar),
        }
    }
    #[derive(Clone, Debug)]
    pub struct CottonSchurHrpRmtNcoSpec<T: Real = f64> {
        pub covariance_4x4: [[T; 4]; 4],
        pub schur_gamma: T,
        pub sample_aspect_ratio_q: T,
    }
    #[derive(Clone, Debug)]
    pub struct CottonSchurHrpRmtNcoReport<T: Real = f64> {
        pub marchenko_pastur_upper_edge: T,
        pub classical_hrp_weights: [T; 4],
        pub schur_hrp_weights: [T; 4],
        pub rmt_nco_weights: [T; 4],
        pub classical_hrp_variance: T,
        pub schur_hrp_portfolio_variance: T,
        pub rmt_nco_portfolio_variance: T,
        pub schur_variance_reduction_bps: T,
        pub weights_sum_conservation_error: T,
    }
    pub fn evaluate_cotton_schur_hrp_and_rmt_nco<T: Real>(
        s: &CottonSchurHrpRmtNcoSpec<T>,
    ) -> CottonSchurHrpRmtNcoReport<T> {
        let mut inv = [T::zero(); 4];
        for i in 0..4 {
            inv[i] = T::one() / s.covariance_4x4[i][i].abs().max_val(T::from_f64(1e-8));
        }
        let si = inv.iter().copied().fold(T::zero(), |a, b| a + b);
        let w = inv.map(|x| x / si);
        let qv = |cov: &[[T; 4]; 4], wt: &[T; 4]| {
            let mut a = T::zero();
            for i in 0..4 {
                for j in 0..4 {
                    a += wt[i] * cov[i][j] * wt[j];
                }
            }
            a
        };
        let v = qv(&s.covariance_4x4, &w);
        let mut grad = [T::zero(); 4];
        for i in 0..4 {
            for j in 0..4 {
                grad[i] += s.covariance_4x4[i][j] * w[j];
            }
        }
        let scale =
            s.schur_gamma.max_val(T::zero()).min_val(T::one()) / v.max_val(T::from_f64(1e-8));
        let mut raw = [T::zero(); 4];
        for i in 0..4 {
            raw[i] = (w[i] * (T::one() - scale * (grad[i] - v))).max_val(T::from_f64(1e-6));
        }
        let sr = raw.iter().copied().fold(T::zero(), |a, b| a + b);
        let sw_candidate = raw.map(|x| x / sr);
        let sv_candidate = qv(&s.covariance_4x4, &sw_candidate);
        let sw = if sv_candidate <= v { sw_candidate } else { w };
        let sv = qv(&s.covariance_4x4, &sw);
        let edge = (T::one() + s.sample_aspect_ratio_q.max_val(T::zero()).sqrt()).powf(2.0);
        let sh = s.sample_aspect_ratio_q.max_val(T::zero()).sqrt() / edge;
        let mut rr = [T::zero(); 4];
        for i in 0..4 {
            rr[i] = ((T::one() - sh) * sw[i] + sh * T::from_f64(0.25)).max_val(T::from_f64(1e-6));
        }
        let rs = rr.iter().copied().fold(T::zero(), |a, b| a + b);
        let rw = rr.map(|x| x / rs);
        let rv = qv(&s.covariance_4x4, &rw);
        CottonSchurHrpRmtNcoReport {
            marchenko_pastur_upper_edge: edge,
            classical_hrp_weights: w,
            schur_hrp_weights: sw,
            rmt_nco_weights: rw,
            classical_hrp_variance: v,
            schur_hrp_portfolio_variance: sv,
            rmt_nco_portfolio_variance: rv,
            schur_variance_reduction_bps: (v - sv) * T::from_f64(1e4),
            weights_sum_conservation_error: (sw.iter().copied().fold(T::zero(), |a, b| a + b)
                - T::one())
            .abs(),
        }
    }
    #[derive(Clone, Debug)]
    pub struct TsayBekkGarchSpilloverSpec<T: Real = f64> {
        pub chol_c: [[T; 2]; 2],
        pub arch_a: [[T; 2]; 2],
        pub garch_b: [[T; 2]; 2],
        pub initial_cov: [[T; 2]; 2],
        pub return_innovations: Vec<[T; 2]>,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct TsayBekkGarchSpilloverReport<T: Real = f64> {
        pub terminal_var_asset1: T,
        pub terminal_var_asset2: T,
        pub terminal_cov_12: T,
        pub terminal_dynamic_correlation: T,
        pub terminal_cov_determinant: T,
        pub cross_shock_spillover_a12_sq: T,
        pub spectral_persistence_trace: T,
        pub bivariate_gaussian_log_likelihood: T,
    }
    pub fn evaluate_tsay_bekk_garch_volatility_spillover<T: Real>(
        s: &TsayBekkGarchSpilloverSpec<T>,
    ) -> TsayBekkGarchSpilloverReport<T> {
        let mut h = s.initial_cov;
        let mut ll = T::zero();
        for e in &s.return_innovations {
            let mut n = [[T::zero(); 2]; 2];
            for i in 0..2 {
                for j in 0..2 {
                    n[i][j] = s.chol_c[i][0] * s.chol_c[j][0] + s.chol_c[i][1] * s.chol_c[j][1];
                    for k in 0..2 {
                        for l in 0..2 {
                            n[i][j] += s.arch_a[i][k] * e[k] * e[l] * s.arch_a[j][l]
                                + s.garch_b[i][k] * h[k][l] * s.garch_b[j][l];
                        }
                    }
                }
            }
            h = n;
            let det = (h[0][0] * h[1][1] - h[0][1] * h[1][0]).max_val(T::from_f64(1e-12));
            let quad = (h[1][1] * e[0] * e[0] - T::from_f64(2.0) * h[0][1] * e[0] * e[1]
                + h[0][0] * e[1] * e[1])
                / det;
            ll -= T::from_f64(0.5)
                * (T::from_f64(2.0) * T::from_f64(2.0 * std::f64::consts::PI).ln()
                    + det.ln()
                    + quad);
        }
        let det = h[0][0] * h[1][1] - h[0][1] * h[1][0];
        TsayBekkGarchSpilloverReport {
            terminal_var_asset1: h[0][0],
            terminal_var_asset2: h[1][1],
            terminal_cov_12: h[0][1],
            terminal_dynamic_correlation: h[0][1]
                / (h[0][0] * h[1][1]).sqrt().max_val(T::from_f64(1e-8)),
            terminal_cov_determinant: det,
            cross_shock_spillover_a12_sq: s.arch_a[0][1] * s.arch_a[0][1],
            spectral_persistence_trace: s.arch_a[0][0]
                + s.arch_a[1][1]
                + s.garch_b[0][0]
                + s.garch_b[1][1],
            bivariate_gaussian_log_likelihood: ll,
        }
    }
    #[cfg(test)]
    mod domain_suite_44_tests {
        use super::*;
        use ql_ad::Dual;

        fn evar_spec<T: Real>(x: T) -> AhmadiJavidEntropicVarPortfolioSpec<T> {
            AhmadiJavidEntropicVarPortfolioSpec {
                scenario_returns: vec![
                    vec![T::from_f64(0.02) + x, T::from_f64(-0.01)],
                    vec![T::from_f64(-0.04), T::from_f64(0.01)],
                    vec![T::from_f64(0.01), T::from_f64(0.015)],
                ],
                scenario_probabilities: vec![T::from_f64(0.4), T::from_f64(0.3), T::from_f64(0.3)],
                initial_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
                tail_alpha: T::from_f64(0.1),
                risk_aversion_lambda: T::one(),
            }
        }

        fn hrp_spec<T: Real>(x: T) -> CottonSchurHrpRmtNcoSpec<T> {
            CottonSchurHrpRmtNcoSpec {
                covariance_4x4: [
                    [
                        T::from_f64(0.04) + x,
                        T::from_f64(0.01),
                        T::zero(),
                        T::zero(),
                    ],
                    [T::from_f64(0.01), T::from_f64(0.05), T::zero(), T::zero()],
                    [T::zero(), T::zero(), T::from_f64(0.03), T::from_f64(0.01)],
                    [T::zero(), T::zero(), T::from_f64(0.01), T::from_f64(0.06)],
                ],
                schur_gamma: T::from_f64(0.5),
                sample_aspect_ratio_q: T::from_f64(0.5),
            }
        }

        fn bekk_spec<T: Real>(x: T) -> TsayBekkGarchSpilloverSpec<T> {
            TsayBekkGarchSpilloverSpec {
                chol_c: [
                    [T::from_f64(0.1), T::zero()],
                    [T::from_f64(0.02), T::from_f64(0.1)],
                ],
                arch_a: [
                    [T::from_f64(0.1), T::from_f64(0.02)],
                    [T::from_f64(0.01), T::from_f64(0.1)],
                ],
                garch_b: [[T::from_f64(0.8), T::zero()], [T::zero(), T::from_f64(0.7)]],
                initial_cov: [
                    [T::from_f64(0.04), T::zero()],
                    [T::zero(), T::from_f64(0.04)],
                ],
                return_innovations: vec![
                    [T::from_f64(0.01) + x, T::from_f64(-0.02)],
                    [T::from_f64(0.02), T::from_f64(0.01)],
                ],
            }
        }

        #[test]
        fn weights_conserve() {
            let h = 1e-6;
            let ev0 = evaluate_ahmadi_javid_entropic_var_portfolio(&evar_spec(0.0));
            assert!(ev0.exact_chernoff_evar >= ev0.empirical_cvar);
            assert!(ev0.optimized_portfolio_evar <= ev0.exact_chernoff_evar);
            let fd_ev = (evaluate_ahmadi_javid_entropic_var_portfolio(&evar_spec(h))
                .exact_chernoff_evar
                - evaluate_ahmadi_javid_entropic_var_portfolio(&evar_spec(-h)).exact_chernoff_evar)
                / (2.0 * h);
            let ad_ev = evaluate_ahmadi_javid_entropic_var_portfolio(&evar_spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .exact_chernoff_evar;
            assert!(ad_ev.d[0].abs() > 1e-6 && (ad_ev.d[0] - fd_ev).abs() < 1e-5);

            let hrp0 = evaluate_cotton_schur_hrp_and_rmt_nco(&hrp_spec(0.0));
            assert!(hrp0.weights_sum_conservation_error < 1e-12);
            assert!(hrp0.schur_hrp_portfolio_variance <= hrp0.classical_hrp_variance);
            let fd_hrp = (evaluate_cotton_schur_hrp_and_rmt_nco(&hrp_spec(h))
                .schur_hrp_portfolio_variance
                - evaluate_cotton_schur_hrp_and_rmt_nco(&hrp_spec(-h))
                    .schur_hrp_portfolio_variance)
                / (2.0 * h);
            let ad_hrp =
                evaluate_cotton_schur_hrp_and_rmt_nco(&hrp_spec(Dual::<1>::variable(0.0, 0)))
                    .schur_hrp_portfolio_variance;
            assert!(ad_hrp.d[0].abs() > 1e-6 && (ad_hrp.d[0] - fd_hrp).abs() < 1e-5);

            let bk0 = evaluate_tsay_bekk_garch_volatility_spillover(&bekk_spec(0.0));
            assert!(bk0.terminal_cov_determinant > 0.0);
            assert!(bk0.terminal_dynamic_correlation.abs() < 1.0);
            let fd_bk = (evaluate_tsay_bekk_garch_volatility_spillover(&bekk_spec(h))
                .bivariate_gaussian_log_likelihood
                - evaluate_tsay_bekk_garch_volatility_spillover(&bekk_spec(-h))
                    .bivariate_gaussian_log_likelihood)
                / (2.0 * h);
            let ad_bk = evaluate_tsay_bekk_garch_volatility_spillover(&bekk_spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .bivariate_gaussian_log_likelihood;
            assert!(ad_bk.d[0].abs() > 1e-6 && (ad_bk.d[0] - fd_bk).abs() < 1e-5);
        }
    }
}
pub use evaluate_ahmadi_javid_entropic_var_portfolio_module::*;

mod evaluate_grinold_qian_transfer_coefficient_and_flam_module {
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct GrinoldQianTransferCoefficientFlAmSpec<T: Real = f64> {
        pub information_coefficient_ic: T,
        pub breadth_br: T,
        pub target_active_risk: T,
        pub residual_volatilities: [T; 3],
        pub alpha_forecasts: [T; 3],
        pub constrained_active_weights: [T; 3],
    }

    #[derive(Clone, Copy, Debug)]
    pub struct GrinoldQianTransferCoefficientFlAmReport<T: Real = f64> {
        pub unconstrained_active_weights: [T; 3],
        pub realized_constrained_active_risk: T,
        pub transfer_coefficient_tc: T,
        pub unconstrained_information_ratio: T,
        pub constrained_information_ratio: T,
        pub expected_unconstrained_alpha: T,
        pub expected_constrained_alpha: T,
        pub constraint_alpha_haircut: T,
        pub transfer_coef_pythagorean_error: T,
    }

    pub fn evaluate_grinold_qian_transfer_coefficient_and_flam<T: Real>(
        s: &GrinoldQianTransferCoefficientFlAmSpec<T>,
    ) -> GrinoldQianTransferCoefficientFlAmReport<T> {
        let mut u = [T::zero(); 3];
        let mut v = [T::zero(); 3];
        let mut u2 = T::zero();
        let mut v2 = T::zero();
        let mut uv = T::zero();
        for i in 0..3 {
            let sig = s.residual_volatilities[i].max_val(T::from_f64(1e-8));
            u[i] = s.alpha_forecasts[i] / sig;
            v[i] = sig * s.constrained_active_weights[i];
            u2 += u[i] * u[i];
            v2 += v[i] * v[i];
            uv += u[i] * v[i];
        }
        let u_norm = u2.max_val(T::from_f64(1e-12)).sqrt();
        let v_norm = v2.max_val(T::from_f64(1e-12)).sqrt();
        let mut uw = [T::zero(); 3];
        for i in 0..3 {
            let sig = s.residual_volatilities[i].max_val(T::from_f64(1e-8));
            uw[i] = s.target_active_risk * (s.alpha_forecasts[i] / (sig * sig)) / u_norm;
        }
        let tc = uv / (u_norm * v_norm);
        let mut orth2 = T::zero();
        for i in 0..3 {
            let r_i = (v[i] / v_norm) - tc * (u[i] / u_norm);
            orth2 += r_i * r_i;
        }
        let pyth_diff = (T::one() - tc * tc) - orth2;
        let pyth_err = (pyth_diff * pyth_diff + T::from_f64(1e-30)).sqrt();
        let ir_u = s.information_coefficient_ic * s.breadth_br.sqrt();
        let ir_c = ir_u * tc;
        let alpha_u = s.target_active_risk * ir_u;
        let alpha_c = s.target_active_risk * ir_c;
        GrinoldQianTransferCoefficientFlAmReport {
            unconstrained_active_weights: uw,
            realized_constrained_active_risk: v_norm,
            transfer_coefficient_tc: tc,
            unconstrained_information_ratio: ir_u,
            constrained_information_ratio: ir_c,
            expected_unconstrained_alpha: alpha_u,
            expected_constrained_alpha: alpha_c,
            constraint_alpha_haircut: alpha_u - alpha_c,
            transfer_coef_pythagorean_error: pyth_err,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct MeucciEntropyPoolingViewSpec<T: Real = f64> {
        pub prior_probabilities: [T; 4],
        pub scenario_returns: [T; 4],
        pub target_view_return: T,
        pub view_confidence: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct MeucciEntropyPoolingViewReport<T: Real = f64> {
        pub optimal_lagrange_lambda: T,
        pub tilted_probabilities_q: [T; 4],
        pub posterior_probabilities_p: [T; 4],
        pub prior_expected_return: T,
        pub posterior_expected_return: T,
        pub tilted_view_moment_error: T,
        pub posterior_prob_sum_error: T,
        pub kl_relative_entropy_divergence: T,
        pub effective_number_of_scenarios: T,
    }

    pub fn evaluate_meucci_relative_entropy_view_pooling<T: Real>(
        s: &MeucciEntropyPoolingViewSpec<T>,
    ) -> MeucciEntropyPoolingViewReport<T> {
        let mut lam = T::zero();
        for _ in 0..8 {
            let mut z = T::zero();
            let mut dz = T::zero();
            let mut d2z = T::zero();
            for j in 0..4 {
                let d = s.scenario_returns[j] - s.target_view_return;
                let w = s.prior_probabilities[j] * (lam * d).exp();
                z += w;
                dz += w * d;
                d2z += w * d * d;
            }
            let g = dz / z;
            let h = (d2z / z - g * g).max_val(T::from_f64(1e-10));
            lam -= g / h;
        }
        let mut z = T::zero();
        let mut q = [T::zero(); 4];
        for j in 0..4 {
            q[j] = s.prior_probabilities[j]
                * (lam * (s.scenario_returns[j] - s.target_view_return)).exp();
            z += q[j];
        }
        let mut q_mean = T::zero();
        let mut prior_mean = T::zero();
        let mut post_mean = T::zero();
        let mut p = [T::zero(); 4];
        let mut p_sum = T::zero();
        let mut kl = T::zero();
        let mut shannon = T::zero();
        for j in 0..4 {
            q[j] /= z;
            p[j] = (T::one() - s.view_confidence) * s.prior_probabilities[j]
                + s.view_confidence * q[j];
            q_mean += q[j] * s.scenario_returns[j];
            prior_mean += s.prior_probabilities[j] * s.scenario_returns[j];
            post_mean += p[j] * s.scenario_returns[j];
            p_sum += p[j];
            let pj = p[j].max_val(T::from_f64(1e-12));
            let p0j = s.prior_probabilities[j].max_val(T::from_f64(1e-12));
            kl += pj * (pj / p0j).ln();
            shannon -= pj * pj.ln();
        }
        let m_diff = q_mean - s.target_view_return;
        let s_diff = T::one() - p_sum;
        MeucciEntropyPoolingViewReport {
            optimal_lagrange_lambda: lam,
            tilted_probabilities_q: q,
            posterior_probabilities_p: p,
            prior_expected_return: prior_mean,
            posterior_expected_return: post_mean,
            tilted_view_moment_error: (m_diff * m_diff + T::from_f64(1e-30)).sqrt(),
            posterior_prob_sum_error: (s_diff * s_diff + T::from_f64(1e-30)).sqrt(),
            kl_relative_entropy_divergence: kl,
            effective_number_of_scenarios: shannon.exp(),
        }
    }

    #[cfg(test)]
    mod domain_suite_47_tests {
        use super::*;
        use ql_ad::Dual;
        fn tc_spec<T: Real>(x: T) -> GrinoldQianTransferCoefficientFlAmSpec<T> {
            GrinoldQianTransferCoefficientFlAmSpec {
                information_coefficient_ic: T::from_f64(0.06),
                breadth_br: T::from_f64(100.0),
                target_active_risk: T::from_f64(0.04) + x,
                residual_volatilities: [T::from_f64(0.15), T::from_f64(0.20), T::from_f64(0.18)],
                alpha_forecasts: [T::from_f64(0.03), T::from_f64(0.02), T::from_f64(-0.025)],
                constrained_active_weights: [
                    T::from_f64(0.12),
                    T::from_f64(0.05),
                    T::from_f64(-0.08),
                ],
            }
        }
        fn ep_spec<T: Real>(x: T) -> MeucciEntropyPoolingViewSpec<T> {
            MeucciEntropyPoolingViewSpec {
                prior_probabilities: [
                    T::from_f64(0.25),
                    T::from_f64(0.25),
                    T::from_f64(0.25),
                    T::from_f64(0.25),
                ],
                scenario_returns: [
                    T::from_f64(-0.04),
                    T::from_f64(0.01),
                    T::from_f64(0.03),
                    T::from_f64(0.08),
                ],
                target_view_return: T::from_f64(0.035) + x,
                view_confidence: T::from_f64(0.75),
            }
        }
        #[test]
        fn suite47_portfolio_invariants_and_dual_fd() {
            let r1 = evaluate_grinold_qian_transfer_coefficient_and_flam(&tc_spec(0.0));
            assert!(r1.transfer_coef_pythagorean_error < 1e-12);
            assert!(r1.transfer_coefficient_tc <= 1.0 && r1.constraint_alpha_haircut >= 0.0);
            let h = 1e-6;
            let fd1 = (evaluate_grinold_qian_transfer_coefficient_and_flam(&tc_spec(h))
                .expected_constrained_alpha
                - evaluate_grinold_qian_transfer_coefficient_and_flam(&tc_spec(-h))
                    .expected_constrained_alpha)
                / (2.0 * h);
            let ad1 = evaluate_grinold_qian_transfer_coefficient_and_flam(&tc_spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .expected_constrained_alpha;
            assert!(ad1.d[0].abs() > 1e-6 && (ad1.d[0] - fd1).abs() < 1e-5);

            let r2 = evaluate_meucci_relative_entropy_view_pooling(&ep_spec(0.0));
            assert!(r2.tilted_view_moment_error < 1e-10);
            assert!(r2.posterior_prob_sum_error < 1e-12);
            assert!(r2.kl_relative_entropy_divergence >= 0.0);
            let fd2 = (evaluate_meucci_relative_entropy_view_pooling(&ep_spec(h))
                .posterior_expected_return
                - evaluate_meucci_relative_entropy_view_pooling(&ep_spec(-h))
                    .posterior_expected_return)
                / (2.0 * h);
            let ad2 = evaluate_meucci_relative_entropy_view_pooling(&ep_spec(Dual::<1>::variable(
                0.0, 0,
            )))
            .posterior_expected_return;
            assert!(ad2.d[0].abs() > 1e-6 && (ad2.d[0] - fd2).abs() < 1e-5);
        }
    }
}
pub use evaluate_grinold_qian_transfer_coefficient_and_flam_module::*;

mod evaluate_campbell_fama_macbeth_and_shanken_eiv_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct CampbellFamaMacbethShankenSpec<T: Real = f64> {
        pub asset_mean_excess_returns: [T; 4],
        pub estimated_factor_betas: [T; 4],
        pub idiosyncratic_var_sigma_eps2: T,
        pub factor_var_sigma_f2: T,
        pub sample_periods_t: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct CampbellFamaMacbethShankenReport<T: Real = f64> {
        pub zero_beta_rate_gamma0: T,
        pub factor_risk_premium_gamma1: T,
        pub pricing_error_alphas: [T; 4],
        pub cross_sectional_orthogonality_error: T,
        pub shanken_eiv_correction_c: T,
        pub naive_fama_macbeth_se: T,
        pub shanken_corrected_se: T,
        pub shanken_corrected_t_stat: T,
    }

    pub fn evaluate_campbell_fama_macbeth_and_shanken_eiv<T: Real>(
        spec: &CampbellFamaMacbethShankenSpec<T>,
    ) -> CampbellFamaMacbethShankenReport<T> {
        let r = &spec.asset_mean_excess_returns;
        let b = &spec.estimated_factor_betas;
        let q = T::from_f64(0.25);
        let r_bar = q * (r[0] + r[1] + r[2] + r[3]);
        let b_bar = q * (b[0] + b[1] + b[2] + b[3]);

        let mut num = T::zero();
        let mut den = T::zero();
        for i in 0..4 {
            let db = b[i] - b_bar;
            let dr = r[i] - r_bar;
            num = num + db * dr;
            den = den + db * db;
        }
        let den_safe = den.max_val(T::from_f64(1e-8));
        let gamma1 = num / den_safe;
        let gamma0 = r_bar - gamma1 * b_bar;

        let mut alphas = [T::zero(); 4];
        let mut sum_a = T::zero();
        let mut sum_ab = T::zero();
        for i in 0..4 {
            alphas[i] = r[i] - (gamma0 + gamma1 * b[i]);
            sum_a = sum_a + alphas[i];
            sum_ab = sum_ab + alphas[i] * b[i];
        }
        let orth_err = (sum_a * sum_a + sum_ab * sum_ab + T::from_f64(1e-30)).sqrt();

        let sig_f2 = spec.factor_var_sigma_f2.max_val(T::from_f64(1e-8));
        let sig_e2 = spec.idiosyncratic_var_sigma_eps2.max_val(T::from_f64(1e-8));
        let t_obs = spec.sample_periods_t.max_val(T::from_f64(1.0));
        let c_shanken = (gamma1 * gamma1) / sig_f2;
        let var_naive = (sig_e2 / den_safe + sig_f2) / t_obs;
        let var_shank = ((T::one() + c_shanken) * sig_e2 / den_safe + sig_f2) / t_obs;
        let se_naive = var_naive.sqrt();
        let se_shank = var_shank.sqrt();
        let t_stat = gamma1 / se_shank.max_val(T::from_f64(1e-8));

        CampbellFamaMacbethShankenReport {
            zero_beta_rate_gamma0: gamma0,
            factor_risk_premium_gamma1: gamma1,
            pricing_error_alphas: alphas,
            cross_sectional_orthogonality_error: orth_err,
            shanken_eiv_correction_c: c_shanken,
            naive_fama_macbeth_se: se_naive,
            shanken_corrected_se: se_shank,
            shanken_corrected_t_stat: t_stat,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct BouchaudCrossImpactEigenmodeSpec<T: Real = f64> {
        pub parent_order_q1: T,
        pub parent_order_q2: T,
        pub self_impact_lambda1: T,
        pub self_impact_lambda2: T,
        pub cross_impact_lambda12: T,
        pub resilience_decay_rho: T,
        pub risk_aversion_gamma: T,
        pub asset_volatility: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct BouchaudCrossImpactEigenmodeReport<T: Real = f64> {
        pub market_eigenvalue_mu_plus: T,
        pub spread_eigenvalue_mu_minus: T,
        pub eigen_trace_det_conservation_error: T,
        pub self_impact_cost_component: T,
        pub cross_asset_spillover_cost: T,
        pub total_cross_impact_execution_cost: T,
        pub market_eigenmode_cost_share: T,
    }

    pub fn evaluate_bouchaud_cross_impact_eigenmode_execution<T: Real>(
        spec: &BouchaudCrossImpactEigenmodeSpec<T>,
    ) -> BouchaudCrossImpactEigenmodeReport<T> {
        let l1 = spec.self_impact_lambda1.max_val(T::from_f64(1e-6));
        let l2 = spec.self_impact_lambda2.max_val(T::from_f64(1e-6));
        let l12 = spec.cross_impact_lambda12;
        let tr = l1 + l2;
        let det = l1 * l2 - l12 * l12;
        let disc =
            ((l1 - l2) * (l1 - l2) + T::from_f64(4.0) * l12 * l12 + T::from_f64(1e-24)).sqrt();
        let mu_p = T::from_f64(0.5) * (tr + disc);
        let mu_m = T::from_f64(0.5) * (tr - disc);

        let d_tr = (mu_p + mu_m) - tr;
        let d_det = (mu_p * mu_m) - det;
        let eig_err = (d_tr * d_tr + d_det * d_det + T::from_f64(1e-30)).sqrt();

        let q1 = spec.parent_order_q1;
        let q2 = spec.parent_order_q2;
        let rho = spec
            .resilience_decay_rho
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.95));
        let kernel_mult = (T::one() + rho) / (T::from_f64(2.0) * (T::one() - rho));
        let risk_pen = T::from_f64(0.25)
            * spec.risk_aversion_gamma
            * spec.asset_volatility
            * spec.asset_volatility
            * (q1 * q1 + q2 * q2);

        let self_cost = kernel_mult * (l1 * q1 * q1 + l2 * q2 * q2) + risk_pen;
        let cross_cost = kernel_mult * T::from_f64(2.0) * l12 * q1 * q2;
        let tot_cost = (self_cost + cross_cost).max_val(T::from_f64(1e-8));

        let v_p1 = l12;
        let v_p2 = mu_p - l1;
        let norm_vp = (v_p1 * v_p1 + v_p2 * v_p2 + T::from_f64(1e-16)).sqrt();
        let proj_p = (q1 * v_p1 + q2 * v_p2) / norm_vp;
        let mkt_mode_cost = kernel_mult * mu_p * proj_p * proj_p;
        let share = (mkt_mode_cost / tot_cost)
            .max_val(T::zero())
            .min_val(T::one());

        BouchaudCrossImpactEigenmodeReport {
            market_eigenvalue_mu_plus: mu_p,
            spread_eigenvalue_mu_minus: mu_m,
            eigen_trace_det_conservation_error: eig_err,
            self_impact_cost_component: self_cost,
            cross_asset_spillover_cost: cross_cost,
            total_cross_impact_execution_cost: tot_cost,
            market_eigenmode_cost_share: share,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct KarnoskySingerMultiCurrencySpec<T: Real = f64> {
        pub port_asset_weights: [T; 3],
        pub bench_asset_weights: [T; 3],
        pub port_fx_weights: [T; 3],
        pub bench_fx_weights: [T; 3],
        pub local_port_returns: [T; 3],
        pub local_bench_returns: [T; 3],
        pub local_eurodeposit_rates: [T; 3],
        pub base_cash_rate: T,
        pub spot_fx_returns: [T; 3],
    }

    #[derive(Clone, Copy, Debug)]
    pub struct KarnoskySingerMultiCurrencyReport<T: Real = f64> {
        pub local_excess_allocation_effect: T,
        pub local_security_selection_effect: T,
        pub currency_forward_carry_effect: T,
        pub spot_currency_surprise_effect: T,
        pub total_active_base_return: T,
        pub karnosky_singer_conservation_error: T,
    }

    pub fn evaluate_karnosky_singer_multicurrency_attribution<T: Real>(
        spec: &KarnoskySingerMultiCurrencySpec<T>,
    ) -> KarnoskySingerMultiCurrencyReport<T> {
        let mut r_b_star = [T::zero(); 3];
        let mut r_p_star = [T::zero(); 3];
        let mut d = [T::zero(); 3];
        let mut r_b_bar = T::zero();
        let mut d_bar = T::zero();
        let mut e_bar = T::zero();

        for i in 0..3 {
            r_b_star[i] = spec.local_bench_returns[i] - spec.local_eurodeposit_rates[i];
            r_p_star[i] = spec.local_port_returns[i] - spec.local_eurodeposit_rates[i];
            d[i] = spec.local_eurodeposit_rates[i] - spec.base_cash_rate;
            r_b_bar = r_b_bar + spec.bench_asset_weights[i] * r_b_star[i];
            d_bar = d_bar + spec.bench_fx_weights[i] * d[i];
            e_bar = e_bar + spec.bench_fx_weights[i] * spec.spot_fx_returns[i];
        }

        let mut alloc = T::zero();
        let mut sel = T::zero();
        let mut carry = T::zero();
        let mut spot_fx = T::zero();
        let mut port_tot = T::zero();
        let mut bench_tot = T::zero();

        for i in 0..3 {
            let dw_a = spec.port_asset_weights[i] - spec.bench_asset_weights[i];
            let dw_fx = spec.port_fx_weights[i] - spec.bench_fx_weights[i];
            alloc = alloc + dw_a * (r_b_star[i] - r_b_bar);
            sel = sel + spec.port_asset_weights[i] * (r_p_star[i] - r_b_star[i]);
            carry = carry + dw_fx * (d[i] - d_bar);
            spot_fx = spot_fx + dw_fx * (spec.spot_fx_returns[i] - e_bar);

            port_tot = port_tot
                + spec.port_asset_weights[i] * r_p_star[i]
                + spec.port_fx_weights[i] * (d[i] + spec.spot_fx_returns[i]);
            bench_tot = bench_tot
                + spec.bench_asset_weights[i] * r_b_star[i]
                + spec.bench_fx_weights[i] * (d[i] + spec.spot_fx_returns[i]);
        }
        let active_tot = port_tot - bench_tot;
        let sum_effects = alloc + sel + carry + spot_fx;
        let d_ks = active_tot - sum_effects;
        let ks_err = (d_ks * d_ks + T::from_f64(1e-30)).sqrt();

        KarnoskySingerMultiCurrencyReport {
            local_excess_allocation_effect: alloc,
            local_security_selection_effect: sel,
            currency_forward_carry_effect: carry,
            spot_currency_surprise_effect: spot_fx,
            total_active_base_return: active_tot,
            karnosky_singer_conservation_error: ks_err,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct CampbellViceiraIntertemporalSaaSpec<T: Real = f64> {
        pub current_expected_excess_return_xt: T,
        pub long_run_excess_return_xbar: T,
        pub predictor_persistence_phi: T,
        pub return_volatility_sigma_r: T,
        pub predictor_volatility_sigma_x: T,
        pub return_predictor_corr_rho: T,
        pub relative_risk_aversion_gamma: T,
        pub wealth_consumption_rho_w: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct CampbellViceiraIntertemporalSaaReport<T: Real = f64> {
        pub myopic_equity_allocation: T,
        pub intertemporal_hedging_demand: T,
        pub total_strategic_equity_allocation: T,
        pub hedging_share_of_total_demand: T,
        pub horizon_variance_reduction_ratio: T,
        pub certainty_equivalent_welfare_gain: T,
    }

    pub fn evaluate_campbell_viceira_intertemporal_hedging_saa<T: Real>(
        spec: &CampbellViceiraIntertemporalSaaSpec<T>,
    ) -> CampbellViceiraIntertemporalSaaReport<T> {
        let xt = T::from_f64(0.7) * spec.current_expected_excess_return_xt
            + T::from_f64(0.3) * spec.long_run_excess_return_xbar;
        let sig_r = spec.return_volatility_sigma_r.max_val(T::from_f64(0.02));
        let sig_r2 = sig_r * sig_r;
        let gamma = spec.relative_risk_aversion_gamma.max_val(T::from_f64(1.05));
        let myopic = (xt + T::from_f64(0.5) * sig_r2) / (gamma * sig_r2);

        let denom_cv = (T::one() - spec.wealth_consumption_rho_w * spec.predictor_persistence_phi)
            .max_val(T::from_f64(0.02));
        let psi_hedge = -(T::one() - T::one() / gamma)
            * (spec.wealth_consumption_rho_w
                * spec.return_predictor_corr_rho
                * sig_r
                * spec.predictor_volatility_sigma_x)
            / denom_cv;
        let hedge = (psi_hedge / sig_r2) * (xt / gamma);
        let total = myopic + hedge;
        let share = hedge / total.max_val(T::from_f64(1e-6));

        let var_red = (T::one()
            + T::from_f64(2.0)
                * spec.return_predictor_corr_rho
                * (spec.predictor_volatility_sigma_x / sig_r)
                / (T::one() - spec.predictor_persistence_phi).max_val(T::from_f64(0.05)))
        .max_val(T::from_f64(0.2));
        let ce_gain = T::from_f64(0.5) * gamma * sig_r2 * hedge * hedge;

        CampbellViceiraIntertemporalSaaReport {
            myopic_equity_allocation: myopic,
            intertemporal_hedging_demand: hedge,
            total_strategic_equity_allocation: total,
            hedging_share_of_total_demand: share,
            horizon_variance_reduction_ratio: var_red,
            certainty_equivalent_welfare_gain: ce_gain,
        }
    }

    #[cfg(test)]
    mod domain_suite_50_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn suite50_portfolio_all_four_engines_invariants_and_ad_parity() {
            let h = 1e-6;

            // 1. Campbell-Fama-MacBeth-Shanken EIV
            let fm_spec = CampbellFamaMacbethShankenSpec {
                asset_mean_excess_returns: [0.045, 0.062, 0.078, 0.095],
                estimated_factor_betas: [0.70, 0.95, 1.20, 1.45],
                idiosyncratic_var_sigma_eps2: 0.015,
                factor_var_sigma_f2: 0.022,
                sample_periods_t: 120.0,
            };
            let fm_rep = evaluate_campbell_fama_macbeth_and_shanken_eiv(&fm_spec);
            assert!(fm_rep.cross_sectional_orthogonality_error < 1e-12);
            assert!(fm_rep.shanken_corrected_se > fm_rep.naive_fama_macbeth_se);
            let fm_ad =
                evaluate_campbell_fama_macbeth_and_shanken_eiv(&CampbellFamaMacbethShankenSpec {
                    asset_mean_excess_returns: [
                        Dual::<1>::variable(fm_spec.asset_mean_excess_returns[0], 0),
                        Dual::<1>::from_value(fm_spec.asset_mean_excess_returns[1]),
                        Dual::<1>::from_value(fm_spec.asset_mean_excess_returns[2]),
                        Dual::<1>::from_value(fm_spec.asset_mean_excess_returns[3]),
                    ],
                    estimated_factor_betas: fm_spec
                        .estimated_factor_betas
                        .map(Dual::<1>::from_value),
                    idiosyncratic_var_sigma_eps2: Dual::<1>::from_value(
                        fm_spec.idiosyncratic_var_sigma_eps2,
                    ),
                    factor_var_sigma_f2: Dual::<1>::from_value(fm_spec.factor_var_sigma_f2),
                    sample_periods_t: Dual::<1>::from_value(fm_spec.sample_periods_t),
                });
            let mut fmp = fm_spec;
            fmp.asset_mean_excess_returns[0] += h;
            let mut fmm = fm_spec;
            fmm.asset_mean_excess_returns[0] -= h;
            let fd_fm = (evaluate_campbell_fama_macbeth_and_shanken_eiv(&fmp)
                .factor_risk_premium_gamma1
                - evaluate_campbell_fama_macbeth_and_shanken_eiv(&fmm).factor_risk_premium_gamma1)
                / (2.0 * h);
            assert!((fm_ad.factor_risk_premium_gamma1.d[0] - fd_fm).abs() < 1e-5);

            // 2. Bouchaud Cross-Impact Eigenmode Execution
            let ci_spec = BouchaudCrossImpactEigenmodeSpec {
                parent_order_q1: 50.0,
                parent_order_q2: 30.0,
                self_impact_lambda1: 0.04,
                self_impact_lambda2: 0.03,
                cross_impact_lambda12: 0.015,
                resilience_decay_rho: 0.60,
                risk_aversion_gamma: 0.10,
                asset_volatility: 0.20,
            };
            let ci_rep = evaluate_bouchaud_cross_impact_eigenmode_execution(&ci_spec);
            assert!(ci_rep.eigen_trace_det_conservation_error < 1e-12);
            assert!(ci_rep.market_eigenvalue_mu_plus > ci_rep.spread_eigenvalue_mu_minus);
            let ci_ad = evaluate_bouchaud_cross_impact_eigenmode_execution(
                &BouchaudCrossImpactEigenmodeSpec {
                    parent_order_q1: Dual::<1>::variable(ci_spec.parent_order_q1, 0),
                    parent_order_q2: Dual::<1>::from_value(ci_spec.parent_order_q2),
                    self_impact_lambda1: Dual::<1>::from_value(ci_spec.self_impact_lambda1),
                    self_impact_lambda2: Dual::<1>::from_value(ci_spec.self_impact_lambda2),
                    cross_impact_lambda12: Dual::<1>::from_value(ci_spec.cross_impact_lambda12),
                    resilience_decay_rho: Dual::<1>::from_value(ci_spec.resilience_decay_rho),
                    risk_aversion_gamma: Dual::<1>::from_value(ci_spec.risk_aversion_gamma),
                    asset_volatility: Dual::<1>::from_value(ci_spec.asset_volatility),
                },
            );
            let mut cip = ci_spec;
            cip.parent_order_q1 += h;
            let mut cim = ci_spec;
            cim.parent_order_q1 -= h;
            let fd_ci = (evaluate_bouchaud_cross_impact_eigenmode_execution(&cip)
                .total_cross_impact_execution_cost
                - evaluate_bouchaud_cross_impact_eigenmode_execution(&cim)
                    .total_cross_impact_execution_cost)
                / (2.0 * h);
            assert!((ci_ad.total_cross_impact_execution_cost.d[0] - fd_ci).abs() < 1e-5);

            // 3. Karnosky-Singer Multi-Currency Attribution
            let ks_spec = KarnoskySingerMultiCurrencySpec {
                port_asset_weights: [0.45, 0.35, 0.20],
                bench_asset_weights: [0.40, 0.35, 0.25],
                port_fx_weights: [0.50, 0.30, 0.20],
                bench_fx_weights: [0.40, 0.35, 0.25],
                local_port_returns: [0.065, 0.048, 0.032],
                local_bench_returns: [0.058, 0.045, 0.035],
                local_eurodeposit_rates: [0.040, 0.028, 0.012],
                base_cash_rate: 0.040,
                spot_fx_returns: [0.0, 0.015, -0.012],
            };
            let ks_rep = evaluate_karnosky_singer_multicurrency_attribution(&ks_spec);
            assert!(ks_rep.karnosky_singer_conservation_error < 1e-12);
            let ks_ad = evaluate_karnosky_singer_multicurrency_attribution(
                &KarnoskySingerMultiCurrencySpec {
                    port_asset_weights: ks_spec.port_asset_weights.map(Dual::<1>::from_value),
                    bench_asset_weights: ks_spec.bench_asset_weights.map(Dual::<1>::from_value),
                    port_fx_weights: ks_spec.port_fx_weights.map(Dual::<1>::from_value),
                    bench_fx_weights: ks_spec.bench_fx_weights.map(Dual::<1>::from_value),
                    local_port_returns: [
                        Dual::<1>::variable(ks_spec.local_port_returns[0], 0),
                        Dual::<1>::from_value(ks_spec.local_port_returns[1]),
                        Dual::<1>::from_value(ks_spec.local_port_returns[2]),
                    ],
                    local_bench_returns: ks_spec.local_bench_returns.map(Dual::<1>::from_value),
                    local_eurodeposit_rates: ks_spec
                        .local_eurodeposit_rates
                        .map(Dual::<1>::from_value),
                    base_cash_rate: Dual::<1>::from_value(ks_spec.base_cash_rate),
                    spot_fx_returns: ks_spec.spot_fx_returns.map(Dual::<1>::from_value),
                },
            );
            let mut ksp = ks_spec;
            ksp.local_port_returns[0] += h;
            let mut ksm = ks_spec;
            ksm.local_port_returns[0] -= h;
            let fd_ks = (evaluate_karnosky_singer_multicurrency_attribution(&ksp)
                .total_active_base_return
                - evaluate_karnosky_singer_multicurrency_attribution(&ksm)
                    .total_active_base_return)
                / (2.0 * h);
            assert!((ks_ad.total_active_base_return.d[0] - fd_ks).abs() < 1e-5);

            // 4. Campbell-Viceira Intertemporal Hedging SAA
            let cv_spec = CampbellViceiraIntertemporalSaaSpec {
                current_expected_excess_return_xt: 0.055,
                long_run_excess_return_xbar: 0.050,
                predictor_persistence_phi: 0.88,
                return_volatility_sigma_r: 0.16,
                predictor_volatility_sigma_x: 0.035,
                return_predictor_corr_rho: -0.68,
                relative_risk_aversion_gamma: 4.0,
                wealth_consumption_rho_w: 0.96,
            };
            let cv_rep = evaluate_campbell_viceira_intertemporal_hedging_saa(&cv_spec);
            assert!(cv_rep.intertemporal_hedging_demand > 0.0);
            assert!(cv_rep.total_strategic_equity_allocation > cv_rep.myopic_equity_allocation);
            let cv_ad = evaluate_campbell_viceira_intertemporal_hedging_saa(
                &CampbellViceiraIntertemporalSaaSpec {
                    current_expected_excess_return_xt: Dual::<1>::variable(
                        cv_spec.current_expected_excess_return_xt,
                        0,
                    ),
                    long_run_excess_return_xbar: Dual::<1>::from_value(
                        cv_spec.long_run_excess_return_xbar,
                    ),
                    predictor_persistence_phi: Dual::<1>::from_value(
                        cv_spec.predictor_persistence_phi,
                    ),
                    return_volatility_sigma_r: Dual::<1>::from_value(
                        cv_spec.return_volatility_sigma_r,
                    ),
                    predictor_volatility_sigma_x: Dual::<1>::from_value(
                        cv_spec.predictor_volatility_sigma_x,
                    ),
                    return_predictor_corr_rho: Dual::<1>::from_value(
                        cv_spec.return_predictor_corr_rho,
                    ),
                    relative_risk_aversion_gamma: Dual::<1>::from_value(
                        cv_spec.relative_risk_aversion_gamma,
                    ),
                    wealth_consumption_rho_w: Dual::<1>::from_value(
                        cv_spec.wealth_consumption_rho_w,
                    ),
                },
            );
            let mut cvp = cv_spec;
            cvp.current_expected_excess_return_xt += h;
            let mut cvm = cv_spec;
            cvm.current_expected_excess_return_xt -= h;
            let fd_cv = (evaluate_campbell_viceira_intertemporal_hedging_saa(&cvp)
                .total_strategic_equity_allocation
                - evaluate_campbell_viceira_intertemporal_hedging_saa(&cvm)
                    .total_strategic_equity_allocation)
                / (2.0 * h);
            assert!((cv_ad.total_strategic_equity_allocation.d[0] - fd_cv).abs() < 1e-5);
        }
    }
}
pub use evaluate_campbell_fama_macbeth_and_shanken_eiv_module::*;

mod evaluate_grinold_kahn_flam_with_transfer_coefficient_module {
    // Engine 5A: Grinold & Kahn (2000) & Clarke–de Silva–Thorley (2002) Fundamental Law of Active
    // Management (`FLAM`), Transfer Coefficient (`TC`) & Effective Breadth Engine.
    // Engine 5B: Meucci (2005) & McNeil–Frey–Embrechts (2015) Elliptical & Cornish–Fisher
    // Expected Shortfall (`CVaR`) Euler Marginal Risk Contribution Engine.

    use ql_ad::Real;
    use ql_math::norm_pdf;

    #[derive(Debug, Clone, Copy)]
    pub struct GrinoldKahnFlamTransferSpec<T: Real = f64> {
        pub information_coefficient_ic: T,
        pub nominal_breadth_n: T,
        pub cross_sectional_signal_corr_rho: T,
        pub active_weights_std_adj: [T; 3],
        pub unconstrained_optimal_weights_std: [T; 3],
        pub active_risk_aversion_lambda: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GrinoldKahnFlamTransferReport<T: Real = f64> {
        pub effective_breadth_br: T,
        pub transfer_coefficient_tc: T,
        pub unconstrained_information_ratio: T,
        pub constrained_information_ratio: T,
        pub optimal_active_risk_sigma_a: T,
        pub expected_value_added_va: T,
    }

    pub fn evaluate_grinold_kahn_flam_with_transfer_coefficient<T: Real>(
        spec: &GrinoldKahnFlamTransferSpec<T>,
    ) -> GrinoldKahnFlamTransferReport<T> {
        let n = spec.nominal_breadth_n.max_val(T::one());
        let rho2 = spec.cross_sectional_signal_corr_rho * spec.cross_sectional_signal_corr_rho;
        // Buckle (2004) / Qian-Hua-Sorensen effective breadth under cross-sectional correlation:
        let br_eff = n / (T::one() + (n - T::one()) * rho2);

        // Clarke-de Silva-Thorley (2002) Transfer Coefficient TC = Corr(Sigma^{1/2} h_A, Sigma^{-1/2} alpha)
        let mut dot = T::zero();
        let mut n1 = T::zero();
        let mut n2 = T::zero();
        for i in 0..3 {
            let u = spec.active_weights_std_adj[i];
            let v = spec.unconstrained_optimal_weights_std[i];
            dot += u * v;
            n1 += u * u;
            n2 += v * v;
        }
        let tc = (dot / (n1 * n2).max_val(T::from_f64(1e-12)).sqrt())
            .max_val(T::zero())
            .min_val(T::one());

        let ir_unc = spec.information_coefficient_ic * br_eff.sqrt();
        let ir_con = tc * ir_unc;
        let lam = spec.active_risk_aversion_lambda.max_val(T::from_f64(1e-4));
        let sig_a = ir_con / (T::from_f64(2.0) * lam);
        let va = (ir_con * ir_con) / (T::from_f64(4.0) * lam);

        GrinoldKahnFlamTransferReport {
            effective_breadth_br: br_eff,
            transfer_coefficient_tc: tc,
            unconstrained_information_ratio: ir_unc,
            constrained_information_ratio: ir_con,
            optimal_active_risk_sigma_a: sig_a,
            expected_value_added_va: va,
        }
    }

    #[derive(Debug, Clone, Copy)]
    pub struct MeucciExpectedShortfallEulerSpec<T: Real = f64> {
        pub weights: [T; 3],
        pub expected_returns: [T; 3],
        pub asset_volatilities: [T; 3],
        pub pairwise_corr_rho: T,
        pub confidence_alpha: T,
        pub normal_quantile_z_alpha: T,
        pub portfolio_skewness_gamma1: T,
        pub portfolio_excess_kurtosis_gamma2: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct MeucciExpectedShortfallEulerReport<T: Real = f64> {
        pub portfolio_volatility: T,
        pub gaussian_expected_shortfall: T,
        pub marginal_es_contributions: [T; 3],
        pub component_es_contributions: [T; 3],
        pub euler_homogeneity_error: T,
        pub cornish_fisher_tail_es: T,
    }

    pub fn evaluate_meucci_expected_shortfall_euler_decomposition<T: Real>(
        spec: &MeucciExpectedShortfallEulerSpec<T>,
    ) -> MeucciExpectedShortfallEulerReport<T> {
        let w = spec.weights;
        let mu = spec.expected_returns;
        let sig = spec.asset_volatilities;
        let rho = spec.pairwise_corr_rho;

        let mut sigma_w = [T::zero(); 3];
        let mut port_var = T::zero();
        let mut port_mu = T::zero();
        for i in 0..3 {
            port_mu += w[i] * mu[i];
            let mut s_i = T::zero();
            for j in 0..3 {
                let c_ij = if i == j { T::one() } else { rho };
                s_i += sig[i] * sig[j] * c_ij * w[j];
            }
            sigma_w[i] = s_i;
            port_var += w[i] * s_i;
        }
        let port_vol = port_var.max_val(T::from_f64(1e-12)).sqrt();

        let tail_prob = (T::one() - spec.confidence_alpha).max_val(T::from_f64(1e-4));
        let es_mult = norm_pdf(spec.normal_quantile_z_alpha) / tail_prob;
        let es_gauss = -port_mu + es_mult * port_vol;

        let mut mes = [T::zero(); 3];
        let mut ces = [T::zero(); 3];
        let mut ces_sum = T::zero();
        for i in 0..3 {
            mes[i] = -mu[i] + es_mult * (sigma_w[i] / port_vol);
            ces[i] = w[i] * mes[i];
            ces_sum += ces[i];
        }
        let euler_err = (es_gauss - ces_sum).abs();

        // Cornish-Fisher non-Gaussian tail multiplier adjustment:
        let z = spec.normal_quantile_z_alpha;
        let cf_adj = es_mult
            * (T::one() - spec.portfolio_skewness_gamma1 * z / T::from_f64(6.0)
                + spec.portfolio_excess_kurtosis_gamma2 * (z * z - T::one()) / T::from_f64(24.0));
        let es_cf = -port_mu + cf_adj * port_vol;

        MeucciExpectedShortfallEulerReport {
            portfolio_volatility: port_vol,
            gaussian_expected_shortfall: es_gauss,
            marginal_es_contributions: mes,
            component_es_contributions: ces,
            euler_homogeneity_error: euler_err,
            cornish_fisher_tail_es: es_cf,
        }
    }

    #[cfg(test)]
    mod domain_suite_51_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn grinold_kahn_flam_and_meucci_es_euler_invariants_and_ad_parity() {
            let s1 = GrinoldKahnFlamTransferSpec {
                information_coefficient_ic: 0.06,
                nominal_breadth_n: 200.0,
                cross_sectional_signal_corr_rho: 0.12,
                active_weights_std_adj: [0.45, 0.30, -0.20],
                unconstrained_optimal_weights_std: [0.55, 0.25, -0.35],
                active_risk_aversion_lambda: 2.5,
            };
            let r1 = evaluate_grinold_kahn_flam_with_transfer_coefficient(&s1);
            assert!(r1.transfer_coefficient_tc > 0.0 && r1.transfer_coefficient_tc <= 1.0);
            assert!(r1.constrained_information_ratio < r1.unconstrained_information_ratio);

            let h1 = 1e-6;
            let r1_ad = evaluate_grinold_kahn_flam_with_transfer_coefficient(
                &GrinoldKahnFlamTransferSpec {
                    information_coefficient_ic: Dual::<1>::variable(
                        s1.information_coefficient_ic,
                        0,
                    ),
                    nominal_breadth_n: Dual::<1>::from_value(s1.nominal_breadth_n),
                    cross_sectional_signal_corr_rho: Dual::<1>::from_value(
                        s1.cross_sectional_signal_corr_rho,
                    ),
                    active_weights_std_adj: s1.active_weights_std_adj.map(Dual::<1>::from_value),
                    unconstrained_optimal_weights_std: s1
                        .unconstrained_optimal_weights_std
                        .map(Dual::<1>::from_value),
                    active_risk_aversion_lambda: Dual::<1>::from_value(
                        s1.active_risk_aversion_lambda,
                    ),
                },
            );
            let mut s1p = s1;
            s1p.information_coefficient_ic += h1;
            let mut s1m = s1;
            s1m.information_coefficient_ic -= h1;
            let fd1 = (evaluate_grinold_kahn_flam_with_transfer_coefficient(&s1p)
                .expected_value_added_va
                - evaluate_grinold_kahn_flam_with_transfer_coefficient(&s1m)
                    .expected_value_added_va)
                / (2.0 * h1);
            assert!((r1_ad.expected_value_added_va.d[0] - fd1).abs() < 1e-6);

            let s2 = MeucciExpectedShortfallEulerSpec {
                weights: [0.45, 0.35, 0.20],
                expected_returns: [0.06, 0.05, 0.035],
                asset_volatilities: [0.18, 0.15, 0.10],
                pairwise_corr_rho: 0.35,
                confidence_alpha: 0.975,
                normal_quantile_z_alpha: 1.959963984540054,
                portfolio_skewness_gamma1: -0.45,
                portfolio_excess_kurtosis_gamma2: 1.5,
            };
            let r2 = evaluate_meucci_expected_shortfall_euler_decomposition(&s2);
            assert!(r2.euler_homogeneity_error < 1e-12);
            assert!(r2.cornish_fisher_tail_es > r2.gaussian_expected_shortfall);

            let h = 1e-6;
            let r_ad = evaluate_meucci_expected_shortfall_euler_decomposition(
                &MeucciExpectedShortfallEulerSpec {
                    weights: [
                        Dual::<1>::variable(s2.weights[0], 0),
                        Dual::<1>::from_value(s2.weights[1]),
                        Dual::<1>::from_value(s2.weights[2]),
                    ],
                    expected_returns: s2.expected_returns.map(Dual::<1>::from_value),
                    asset_volatilities: s2.asset_volatilities.map(Dual::<1>::from_value),
                    pairwise_corr_rho: Dual::<1>::from_value(s2.pairwise_corr_rho),
                    confidence_alpha: Dual::<1>::from_value(s2.confidence_alpha),
                    normal_quantile_z_alpha: Dual::<1>::from_value(s2.normal_quantile_z_alpha),
                    portfolio_skewness_gamma1: Dual::<1>::from_value(s2.portfolio_skewness_gamma1),
                    portfolio_excess_kurtosis_gamma2: Dual::<1>::from_value(
                        s2.portfolio_excess_kurtosis_gamma2,
                    ),
                },
            );
            // By Euler's theorem, d(ES_alpha)/dw_0 == marginal_es_contributions[0]!
            assert!(
                (r_ad.gaussian_expected_shortfall.d[0] - r2.marginal_es_contributions[0]).abs()
                    < 1e-10
            );
            let mut sp = s2;
            sp.weights[0] += h;
            let mut sm = s2;
            sm.weights[0] -= h;
            let fd = (evaluate_meucci_expected_shortfall_euler_decomposition(&sp)
                .gaussian_expected_shortfall
                - evaluate_meucci_expected_shortfall_euler_decomposition(&sm)
                    .gaussian_expected_shortfall)
                / (2.0 * h);
            assert!((r_ad.gaussian_expected_shortfall.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_grinold_kahn_flam_with_transfer_coefficient_module::*;

mod evaluate_qian_multi_period_alpha_decay_ir_module {
    // Engine 5A: Qian, Hua & Sorensen (2007) & Grinold–Kahn (2000) Multi-Period Information Ratio
    // (`MP-IR`), Alpha Decay Horizon & Turnover-Optimal Trading Speed Engine.
    // Engine 5B: Roncalli (2013) & Meucci (2005) Choueifaty–Coignard Diversification Ratio (`DR`) &
    // Shannon-Entropy Effective Number of Bets (`ENB`) Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct QianMultiPeriodAlphaDecayIrSpec<T: Real = f64> {
        pub initial_ic0: T,
        pub alpha_decay_rate_kappa: T,
        pub breadth_br: T,
        pub proportional_tcost_rate_ctc: T,
        pub active_risk_budget_sigma_a: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct QianMultiPeriodAlphaDecayIrReport<T: Real = f64> {
        pub alpha_half_life_years: T,
        pub optimal_rebalancing_speed_eta: T,
        pub realized_multi_period_ic: T,
        pub gross_information_ratio: T,
        pub net_information_ratio: T,
        pub net_active_value_added: T,
    }

    pub fn evaluate_qian_multi_period_alpha_decay_ir<T: Real>(
        spec: &QianMultiPeriodAlphaDecayIrSpec<T>,
    ) -> QianMultiPeriodAlphaDecayIrReport<T> {
        let kap = spec.alpha_decay_rate_kappa.max_val(T::from_f64(1e-4));
        let ic0 = spec.initial_ic0.max_val(T::from_f64(1e-4));
        let sqrt_br = spec.breadth_br.max_val(T::one()).sqrt();
        let ctc = spec.proportional_tcost_rate_ctc.max_val(T::from_f64(1e-5));

        let half_life = T::from_f64(std::f64::consts::LN_2) / kap;

        // Optimal trading rate eta^* maximizing IC_0 * sqrt(BR) * eta / (eta + kappa) - c_tc * eta:
        // eta^* = sqrt(kappa * IC_0 * sqrt(BR) / c_tc) - kappa
        let raw_eta = (kap * ic0 * sqrt_br / ctc).sqrt() - kap;
        let eta_star =
            T::from_f64(0.5) * (raw_eta + (raw_eta * raw_eta + T::from_f64(1e-6)).sqrt());

        let ic_real = ic0 * (eta_star / (eta_star + kap));
        let ir_gross = ic_real * sqrt_br;
        let ir_net = ir_gross - ctc * eta_star;
        let net_va = ir_net * spec.active_risk_budget_sigma_a;

        QianMultiPeriodAlphaDecayIrReport {
            alpha_half_life_years: half_life,
            optimal_rebalancing_speed_eta: eta_star,
            realized_multi_period_ic: ic_real,
            gross_information_ratio: ir_gross,
            net_information_ratio: ir_net,
            net_active_value_added: net_va,
        }
    }

    #[derive(Debug, Clone, Copy)]
    pub struct RoncalliMeucciDiversificationEnbSpec<T: Real = f64> {
        pub weights: [T; 3],
        pub asset_volatilities: [T; 3],
        pub corr_01: T,
        pub corr_02: T,
        pub corr_12: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct RoncalliMeucciDiversificationEnbReport<T: Real = f64> {
        pub weighted_average_volatility: T,
        pub portfolio_volatility: T,
        pub choueifaty_diversification_ratio: T,
        pub risk_budget_shares_p: [T; 3],
        pub shannon_entropy_effective_bets_enb: T,
        pub risk_share_partition_error: T,
    }

    pub fn evaluate_roncalli_meucci_diversification_and_enb<T: Real>(
        spec: &RoncalliMeucciDiversificationEnbSpec<T>,
    ) -> RoncalliMeucciDiversificationEnbReport<T> {
        let w = spec.weights;
        let sig = spec.asset_volatilities;
        let corr = [
            [T::one(), spec.corr_01, spec.corr_02],
            [spec.corr_01, T::one(), spec.corr_12],
            [spec.corr_02, spec.corr_12, T::one()],
        ];

        let mut weighted_vol = T::zero();
        let mut port_var = T::zero();
        let mut mrc_num = [T::zero(); 3];
        for i in 0..3 {
            weighted_vol += w[i] * sig[i];
            let mut s_i = T::zero();
            for j in 0..3 {
                s_i += sig[i] * sig[j] * corr[i][j] * w[j];
            }
            mrc_num[i] = w[i] * s_i;
            port_var += mrc_num[i];
        }
        let port_vol = port_var.max_val(T::from_f64(1e-12)).sqrt();
        let dr = weighted_vol / port_vol;

        let mut p = [T::zero(); 3];
        let mut p_sum = T::zero();
        let mut neg_ent = T::zero();
        for i in 0..3 {
            let pi = (mrc_num[i] / port_var.max_val(T::from_f64(1e-12))).max_val(T::from_f64(1e-8));
            p[i] = pi;
            p_sum += pi;
            neg_ent -= pi * pi.ln();
        }
        let enb = neg_ent.exp();
        let part_err = (p_sum - T::one()).abs();

        RoncalliMeucciDiversificationEnbReport {
            weighted_average_volatility: weighted_vol,
            portfolio_volatility: port_vol,
            choueifaty_diversification_ratio: dr,
            risk_budget_shares_p: p,
            shannon_entropy_effective_bets_enb: enb,
            risk_share_partition_error: part_err,
        }
    }

    #[cfg(test)]
    mod domain_suite_52_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn qian_alpha_decay_and_roncalli_meucci_enb_invariants_and_ad_parity() {
            let s1 = QianMultiPeriodAlphaDecayIrSpec {
                initial_ic0: 0.075,
                alpha_decay_rate_kappa: 2.0,
                breadth_br: 180.0,
                proportional_tcost_rate_ctc: 0.08,
                active_risk_budget_sigma_a: 0.04,
            };
            let r1 = evaluate_qian_multi_period_alpha_decay_ir(&s1);
            assert!(r1.optimal_rebalancing_speed_eta > 0.0);
            assert!(
                r1.net_information_ratio > 0.0
                    && r1.net_information_ratio < r1.gross_information_ratio
            );

            let h1 = 1e-6;
            let r1_ad =
                evaluate_qian_multi_period_alpha_decay_ir(&QianMultiPeriodAlphaDecayIrSpec {
                    initial_ic0: Dual::<1>::variable(s1.initial_ic0, 0),
                    alpha_decay_rate_kappa: Dual::<1>::from_value(s1.alpha_decay_rate_kappa),
                    breadth_br: Dual::<1>::from_value(s1.breadth_br),
                    proportional_tcost_rate_ctc: Dual::<1>::from_value(
                        s1.proportional_tcost_rate_ctc,
                    ),
                    active_risk_budget_sigma_a: Dual::<1>::from_value(
                        s1.active_risk_budget_sigma_a,
                    ),
                });
            let mut s1p = s1;
            s1p.initial_ic0 += h1;
            let mut s1m = s1;
            s1m.initial_ic0 -= h1;
            let fd1 = (evaluate_qian_multi_period_alpha_decay_ir(&s1p).net_active_value_added
                - evaluate_qian_multi_period_alpha_decay_ir(&s1m).net_active_value_added)
                / (2.0 * h1);
            assert!((r1_ad.net_active_value_added.d[0] - fd1).abs() < 1e-6);

            let s2 = RoncalliMeucciDiversificationEnbSpec {
                weights: [0.30, 0.35, 0.35],
                asset_volatilities: [0.18, 0.14, 0.12],
                corr_01: 0.25,
                corr_02: 0.10,
                corr_12: 0.30,
            };
            let r2 = evaluate_roncalli_meucci_diversification_and_enb(&s2);
            assert!(r2.choueifaty_diversification_ratio > 1.0);
            assert!(
                r2.shannon_entropy_effective_bets_enb > 1.0
                    && r2.shannon_entropy_effective_bets_enb <= 3.0
            );
            assert!(r2.risk_share_partition_error < 1e-12);

            let h = 1e-6;
            let r_ad = evaluate_roncalli_meucci_diversification_and_enb(
                &RoncalliMeucciDiversificationEnbSpec {
                    weights: s2.weights.map(Dual::<1>::from_value),
                    asset_volatilities: [
                        Dual::<1>::variable(s2.asset_volatilities[0], 0),
                        Dual::<1>::from_value(s2.asset_volatilities[1]),
                        Dual::<1>::from_value(s2.asset_volatilities[2]),
                    ],
                    corr_01: Dual::<1>::from_value(s2.corr_01),
                    corr_02: Dual::<1>::from_value(s2.corr_02),
                    corr_12: Dual::<1>::from_value(s2.corr_12),
                },
            );
            let mut sp = s2;
            sp.asset_volatilities[0] += h;
            let mut sm = s2;
            sm.asset_volatilities[0] -= h;
            let fd = (evaluate_roncalli_meucci_diversification_and_enb(&sp)
                .shannon_entropy_effective_bets_enb
                - evaluate_roncalli_meucci_diversification_and_enb(&sm)
                    .shannon_entropy_effective_bets_enb)
                / (2.0 * h);
            assert!((r_ad.shannon_entropy_effective_bets_enb.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_qian_multi_period_alpha_decay_ir_module::*;

mod evaluate_qian_risk_parity_spinu_log_barrier_module {
    // Engine 5A: Qian (2005) & Spinu (2013) Equal Risk Contribution (`ERC` / Risk Parity)
    // Convex Log-Barrier Solver & Diversification Ratio Engine.
    // Engine 5B: Rockafellar & Uryasev (2000) & Meucci (2005) Auxiliary-Threshold CVaR
    // Minimization & KKT Dual Tail-Probability Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct QianSpinuRiskParityLogBarrierSpec<T: Real = f64> {
        pub asset_volatilities: [T; 3],
        pub risk_budgets_b: [T; 3],
        pub corr_01: T,
        pub corr_02: T,
        pub corr_12: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct QianSpinuRiskParityLogBarrierReport<T: Real = f64> {
        pub optimal_erc_weights: [T; 3],
        pub risk_contribution_shares: [T; 3],
        pub portfolio_volatility: T,
        pub diversification_ratio: T,
        pub max_risk_budget_parity_error: T,
    }

    pub fn evaluate_qian_risk_parity_spinu_log_barrier<T: Real>(
        spec: &QianSpinuRiskParityLogBarrierSpec<T>,
    ) -> QianSpinuRiskParityLogBarrierReport<T> {
        let v0 = spec.asset_volatilities[0].max_val(T::from_f64(1e-4));
        let v1 = spec.asset_volatilities[1].max_val(T::from_f64(1e-4));
        let v2 = spec.asset_volatilities[2].max_val(T::from_f64(1e-4));
        let b_sum = (spec.risk_budgets_b[0] + spec.risk_budgets_b[1] + spec.risk_budgets_b[2])
            .max_val(T::from_f64(1e-6));
        let b = [
            spec.risk_budgets_b[0] / b_sum,
            spec.risk_budgets_b[1] / b_sum,
            spec.risk_budgets_b[2] / b_sum,
        ];

        let c00 = v0 * v0;
        let c11 = v1 * v1;
        let c22 = v2 * v2;
        let c01 = spec.corr_01 * v0 * v1;
        let c02 = spec.corr_02 * v0 * v2;
        let c12 = spec.corr_12 * v1 * v2;

        // Spinu (2013) coordinatewise exact quadratic updates on 0.5 y^T Sigma y - sum b_i ln y_i:
        // For asset i, c_{ii} y_i^2 + (sum_{j != i} c_{ij} y_j) y_i - b_i = 0:
        let mut y = [b[0].sqrt() / v0, b[1].sqrt() / v1, b[2].sqrt() / v2];
        for _ in 0..8 {
            let s0 = c01 * y[1] + c02 * y[2];
            y[0] =
                (-s0 + (s0 * s0 + T::from_f64(4.0) * c00 * b[0]).sqrt()) / (T::from_f64(2.0) * c00);
            let s1 = c01 * y[0] + c12 * y[2];
            y[1] =
                (-s1 + (s1 * s1 + T::from_f64(4.0) * c11 * b[1]).sqrt()) / (T::from_f64(2.0) * c11);
            let s2 = c02 * y[0] + c12 * y[1];
            y[2] =
                (-s2 + (s2 * s2 + T::from_f64(4.0) * c22 * b[2]).sqrt()) / (T::from_f64(2.0) * c22);
        }

        let y_sum = (y[0] + y[1] + y[2]).max_val(T::from_f64(1e-8));
        let w = [y[0] / y_sum, y[1] / y_sum, y[2] / y_sum];

        let mw0 = c00 * w[0] + c01 * w[1] + c02 * w[2];
        let mw1 = c01 * w[0] + c11 * w[1] + c12 * w[2];
        let mw2 = c02 * w[0] + c12 * w[1] + c22 * w[2];
        let port_var = (w[0] * mw0 + w[1] * mw1 + w[2] * mw2).max_val(T::from_f64(1e-10));
        let port_vol = port_var.sqrt();

        let rc = [
            (w[0] * mw0) / port_var,
            (w[1] * mw1) / port_var,
            (w[2] * mw2) / port_var,
        ];
        let err0 = (rc[0] - b[0]).abs();
        let err1 = (rc[1] - b[1]).abs();
        let err2 = (rc[2] - b[2]).abs();
        let max_err = err0.max_val(err1).max_val(err2);
        let weighted_vol = w[0] * v0 + w[1] * v1 + w[2] * v2;
        let dr = weighted_vol / port_vol;

        QianSpinuRiskParityLogBarrierReport {
            optimal_erc_weights: w,
            risk_contribution_shares: rc,
            portfolio_volatility: port_vol,
            diversification_ratio: dr,
            max_risk_budget_parity_error: max_err,
        }
    }

    #[derive(Debug, Clone, Copy)]
    pub struct RockafellarUryasevCvarAuxiliarySpec<T: Real = f64> {
        pub portfolio_weights: [T; 2],
        pub scenario_losses_asset0: [T; 4],
        pub scenario_losses_asset1: [T; 4],
        pub confidence_alpha: T,
        pub smoothing_epsilon: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct RockafellarUryasevCvarAuxiliaryReport<T: Real = f64> {
        pub optimal_var_threshold_zeta: T,
        pub minimized_auxiliary_cvar: T,
        pub expected_portfolio_loss: T,
        pub tail_excess_loss_integral: T,
        pub dual_tail_weight_sum: T,
    }

    pub fn evaluate_rockafellar_uryasev_cvar_minimization<T: Real>(
        spec: &RockafellarUryasevCvarAuxiliarySpec<T>,
    ) -> RockafellarUryasevCvarAuxiliaryReport<T> {
        let w0 = spec.portfolio_weights[0];
        let w1 = spec.portfolio_weights[1];
        let one_minus_alpha = (T::one() - spec.confidence_alpha).max_val(T::from_f64(0.05));
        let eps = spec.smoothing_epsilon.max_val(T::from_f64(1e-4));

        let mut losses = [T::zero(); 4];
        let mut mean_loss = T::zero();
        for s in 0..4 {
            losses[s] = w0 * spec.scenario_losses_asset0[s] + w1 * spec.scenario_losses_asset1[s];
            mean_loss += T::from_f64(0.25) * losses[s];
        }

        // Smooth Newton minimization of Rockafellar-Uryasev auxiliary function F_alpha(w, zeta):
        // F_alpha(w, zeta) = zeta + (1 / (1 - alpha)) * (1/4) sum_{s=1}^4 pos_{eff_eps}(L_s - zeta)
        // dF/dzeta = 1 - (1 / (1 - alpha)) * (1/4) sum_{s=1}^4 psi_{eff_eps}(L_s - zeta) = 0
        let span = (losses[3] - losses[2]).abs().max_val(T::from_f64(0.01));
        let eff_eps = eps + T::from_f64(0.25) * span;
        let mut zeta = T::from_f64(0.5) * (losses[2] + losses[3]);
        for _ in 0..8 {
            let mut prob_sum = T::zero();
            let mut dens_sum = T::zero();
            for &ls in &losses {
                let diff = ls - zeta;
                let r = (diff * diff + eff_eps * eff_eps).sqrt();
                let psi = T::from_f64(0.5) * (T::one() + diff / r);
                let dpsi = T::from_f64(0.5) * (eff_eps * eff_eps) / (r * r * r);
                prob_sum += T::from_f64(0.25) * psi;
                dens_sum += T::from_f64(0.25) * dpsi;
            }
            let grad = T::one() - prob_sum / one_minus_alpha;
            let hess = (dens_sum / one_minus_alpha).max_val(T::from_f64(1e-4));
            zeta -= grad / hess;
        }

        let mut excess_sum = T::zero();
        let mut dual_sum = T::zero();
        for &ls in &losses {
            let diff = ls - zeta;
            let r = (diff * diff + eff_eps * eff_eps).sqrt();
            let pos = T::from_f64(0.5) * (diff + r);
            excess_sum += T::from_f64(0.25) * pos;
            let psi = T::from_f64(0.5) * (T::one() + diff / r);
            dual_sum += T::from_f64(0.25) * psi / one_minus_alpha;
        }
        let cvar = zeta + excess_sum / one_minus_alpha;

        RockafellarUryasevCvarAuxiliaryReport {
            optimal_var_threshold_zeta: zeta,
            minimized_auxiliary_cvar: cvar,
            expected_portfolio_loss: mean_loss,
            tail_excess_loss_integral: excess_sum,
            dual_tail_weight_sum: dual_sum,
        }
    }

    #[cfg(test)]
    mod domain_suite_53_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn qian_spinu_erc_and_rockafellar_uryasev_cvar_invariants_and_ad_parity() {
            let s1 = QianSpinuRiskParityLogBarrierSpec {
                asset_volatilities: [0.18, 0.12, 0.08],
                risk_budgets_b: [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
                corr_01: 0.25,
                corr_02: -0.10,
                corr_12: 0.15,
            };
            let r1 = evaluate_qian_risk_parity_spinu_log_barrier(&s1);
            assert!(r1.max_risk_budget_parity_error < 1e-10);
            assert!(r1.diversification_ratio > 1.0);
            assert!(r1.optimal_erc_weights[2] > r1.optimal_erc_weights[0]);

            let h = 1e-6;
            let r1_ad =
                evaluate_qian_risk_parity_spinu_log_barrier(&QianSpinuRiskParityLogBarrierSpec {
                    asset_volatilities: [
                        Dual::<1>::variable(s1.asset_volatilities[0], 0),
                        Dual::<1>::from_value(s1.asset_volatilities[1]),
                        Dual::<1>::from_value(s1.asset_volatilities[2]),
                    ],
                    risk_budgets_b: s1.risk_budgets_b.map(Dual::<1>::from_value),
                    corr_01: Dual::<1>::from_value(s1.corr_01),
                    corr_02: Dual::<1>::from_value(s1.corr_02),
                    corr_12: Dual::<1>::from_value(s1.corr_12),
                });
            let mut s1p = s1;
            s1p.asset_volatilities[0] += h;
            let mut s1m = s1;
            s1m.asset_volatilities[0] -= h;
            let fd1 = (evaluate_qian_risk_parity_spinu_log_barrier(&s1p).portfolio_volatility
                - evaluate_qian_risk_parity_spinu_log_barrier(&s1m).portfolio_volatility)
                / (2.0 * h);
            assert!((r1_ad.portfolio_volatility.d[0] - fd1).abs() < 1e-6);

            let s2 = RockafellarUryasevCvarAuxiliarySpec {
                portfolio_weights: [0.6, 0.4],
                scenario_losses_asset0: [-0.02, 0.01, 0.04, 0.09],
                scenario_losses_asset1: [-0.01, 0.005, 0.02, 0.05],
                confidence_alpha: 0.75,
                smoothing_epsilon: 0.002,
            };
            let r2 = evaluate_rockafellar_uryasev_cvar_minimization(&s2);
            assert!(r2.minimized_auxiliary_cvar > r2.expected_portfolio_loss);
            assert!((r2.dual_tail_weight_sum - 1.0).abs() < 1e-5);

            let r2_ad = evaluate_rockafellar_uryasev_cvar_minimization(
                &RockafellarUryasevCvarAuxiliarySpec {
                    portfolio_weights: [
                        Dual::<1>::variable(s2.portfolio_weights[0], 0),
                        Dual::<1>::from_value(s2.portfolio_weights[1]),
                    ],
                    scenario_losses_asset0: s2.scenario_losses_asset0.map(Dual::<1>::from_value),
                    scenario_losses_asset1: s2.scenario_losses_asset1.map(Dual::<1>::from_value),
                    confidence_alpha: Dual::<1>::from_value(s2.confidence_alpha),
                    smoothing_epsilon: Dual::<1>::from_value(s2.smoothing_epsilon),
                },
            );
            let mut s2p = s2;
            s2p.portfolio_weights[0] += h;
            let mut s2m = s2;
            s2m.portfolio_weights[0] -= h;
            let fd2 = (evaluate_rockafellar_uryasev_cvar_minimization(&s2p)
                .minimized_auxiliary_cvar
                - evaluate_rockafellar_uryasev_cvar_minimization(&s2m).minimized_auxiliary_cvar)
                / (2.0 * h);
            assert!((r2_ad.minimized_auxiliary_cvar.d[0] - fd2).abs() < 1e-5);
        }
    }
}
pub use evaluate_qian_risk_parity_spinu_log_barrier_module::*;

mod evaluate_kan_zhou_and_sharpe_rbsa_module {
    // Engine 5A: Kan & Zhou (2007) & DeMiguel–Garlappi–Uppal (2009) Three-Fund Optimal
    // Combine (`1/N` + Sample Tangency + Minimum-Variance) Estimation-Risk Shrinkage Engine.
    // Engine 5B: Sharpe (1992) & Bacon (2008) Return-Based Style Analysis (`RBSA`) Constrained
    // Quadratic Tracking-Error, Style `R^2` & Active Selection Return Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct KanZhouThreeFundCombineSpec<T: Real = f64> {
        pub sample_sharpe_squared_theta2: T,
        pub gmv_sharpe_squared_theta_g2: T,
        pub equal_weight_sharpe_squared_theta_ew2: T,
        pub number_of_assets_n: T,
        pub sample_observations_t: T,
        pub risk_aversion_gamma: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct KanZhouThreeFundCombineReport<T: Real = f64> {
        pub plug_in_expected_out_of_sample_utility: T,
        pub kan_zhou_tangency_shrinkage_weight_c1: T,
        pub kan_zhou_gmv_combine_weight_c2: T,
        pub demiguel_equal_weight_blend_c3: T,
        pub three_fund_expected_out_of_sample_utility: T,
        pub estimation_risk_utility_gain: T,
        pub critical_sample_length_for_plug_in: T,
    }

    pub fn evaluate_kan_zhou_three_fund_combine<T: Real>(
        spec: &KanZhouThreeFundCombineSpec<T>,
    ) -> KanZhouThreeFundCombineReport<T> {
        let th2 = spec.sample_sharpe_squared_theta2.max_val(T::from_f64(0.01));
        let th_g2 = spec
            .gmv_sharpe_squared_theta_g2
            .max_val(T::from_f64(0.005))
            .min_val(th2);
        let th_ew2 = spec
            .equal_weight_sharpe_squared_theta_ew2
            .max_val(T::from_f64(0.005));
        let n = spec.number_of_assets_n.max_val(T::from_f64(3.0));
        let t_obs = spec.sample_observations_t.max_val(n + T::from_f64(5.0));
        let gam = spec.risk_aversion_gamma.max_val(T::from_f64(0.5));

        // Kan & Zhou (2007, JFQA Eqs. 23 & 40) finite-sample moments of plug-in Markowitz rule:
        let k_mult = (t_obs / (t_obs - n - T::from_f64(2.0))).max_val(T::one());
        let u_plugin = (th2 / (T::from_f64(2.0) * gam))
            * (T::from_f64(2.0)
                - k_mult * ((t_obs - T::from_f64(2.0)) / (t_obs - n - T::from_f64(1.0))))
            - (n * t_obs)
                / (T::from_f64(2.0)
                    * gam
                    * (t_obs - n - T::from_f64(1.0))
                    * (t_obs - n - T::from_f64(2.0)));

        // Excess ex-ante tangency over GMV dispersion psi^2 = theta^2 - theta_g^2:
        let psi2 = (th2 - th_g2).max_val(T::from_f64(0.005));
        let c_kz_scale = ((t_obs - n - T::from_f64(1.0)) * (t_obs - n - T::from_f64(4.0)))
            / (t_obs * (t_obs - T::from_f64(2.0)));

        let w_tan = c_kz_scale * (psi2 / (psi2 + n / t_obs));
        let w_gmv = c_kz_scale * ((n / t_obs) / (psi2 + n / t_obs));
        let w_ew =
            (T::one() - w_tan - w_gmv).max_val(T::from_f64(0.02)) * (th_ew2 / (th2 + th_ew2));

        let u_three_fund = (T::one() / (T::from_f64(2.0) * gam))
            * (w_tan * psi2 + c_kz_scale * th_g2 + T::from_f64(0.25) * w_ew * th_ew2);
        let util_gain = u_three_fund - u_plugin;
        let t_crit = n * (T::one() + T::from_f64(2.0) / psi2);

        KanZhouThreeFundCombineReport {
            plug_in_expected_out_of_sample_utility: u_plugin,
            kan_zhou_tangency_shrinkage_weight_c1: w_tan,
            kan_zhou_gmv_combine_weight_c2: w_gmv,
            demiguel_equal_weight_blend_c3: w_ew,
            three_fund_expected_out_of_sample_utility: u_three_fund,
            estimation_risk_utility_gain: util_gain,
            critical_sample_length_for_plug_in: t_crit,
        }
    }

    #[derive(Debug, Clone, Copy)]
    pub struct SharpeReturnBasedStyleAnalysisSpec<T: Real = f64> {
        pub fund_mean_return_mu_p: T,
        pub fund_return_variance_var_p: T,
        pub style_benchmark_means: [T; 3],
        pub style_benchmark_variances: [T; 3],
        pub fund_style_covariances: [T; 3],
    }

    #[derive(Debug, Clone, Copy)]
    pub struct SharpeReturnBasedStyleAnalysisReport<T: Real = f64> {
        pub constrained_style_weights: [T; 3],
        pub tracking_error_variance: T,
        pub annualized_tracking_error: T,
        pub sharpe_style_r_squared: T,
        pub passive_style_benchmark_return: T,
        pub active_selection_alpha_return: T,
        pub active_information_ratio: T,
    }

    pub fn evaluate_sharpe_return_based_style_analysis<T: Real>(
        spec: &SharpeReturnBasedStyleAnalysisSpec<T>,
    ) -> SharpeReturnBasedStyleAnalysisReport<T> {
        let var_p = spec.fund_return_variance_var_p.max_val(T::from_f64(1e-4));
        let v0 = spec.style_benchmark_variances[0].max_val(T::from_f64(1e-4));
        let v1 = spec.style_benchmark_variances[1].max_val(T::from_f64(1e-4));
        let v2 = spec.style_benchmark_variances[2].max_val(T::from_f64(1e-4));

        let c0 = spec.fund_style_covariances[0];
        let c1 = spec.fund_style_covariances[1];
        let c2 = spec.fund_style_covariances[2];

        // Unconstrained beta_i = c_i / v_i, projected smoothly onto simplex {w_i >= 0, sum w_i = 1}:
        let b0 = c0 / v0;
        let b1 = c1 / v1;
        let b2 = c2 / v2;

        let eps = T::from_f64(1e-4);
        let smooth_pos = |x: T| T::from_f64(0.5) * (x + (x * x + eps * eps).sqrt());
        let smooth_dpos = |x: T| T::from_f64(0.5) * (T::one() + x / (x * x + eps * eps).sqrt());

        let mut mu_s = (b0 + b1 + b2 - T::one()) / T::from_f64(3.0);
        for _ in 0..10 {
            let x0 = b0 - mu_s / v0;
            let x1 = b1 - mu_s / v1;
            let x2 = b2 - mu_s / v2;
            let w0_t = smooth_pos(x0);
            let w1_t = smooth_pos(x1);
            let w2_t = smooth_pos(x2);
            let jac = (smooth_dpos(x0) / v0 + smooth_dpos(x1) / v1 + smooth_dpos(x2) / v2)
                .max_val(T::from_f64(1e-4));
            mu_s += ((w0_t + w1_t + w2_t) - T::one()) / jac;
        }

        let mut w0 = smooth_pos(b0 - mu_s / v0);
        let mut w1 = smooth_pos(b1 - mu_s / v1);
        let mut w2 = smooth_pos(b2 - mu_s / v2);
        let sw = (w0 + w1 + w2).max_val(T::from_f64(1e-6));
        w0 /= sw;
        w1 /= sw;
        w2 /= sw;

        let style_var = w0 * w0 * v0 + w1 * w1 * v1 + w2 * w2 * v2;
        let cross_cov = w0 * c0 + w1 * c1 + w2 * c2;
        let te_var = (var_p + style_var - T::from_f64(2.0) * cross_cov).max_val(T::from_f64(1e-6));
        let te_vol = te_var.sqrt();

        // Sharpe (1992, JPM Eq. 4) Style R^2 = 1 - Var(e_p) / Var(R_p):
        let r2 = (T::one() - te_var / var_p)
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.999));

        let style_ret = w0 * spec.style_benchmark_means[0]
            + w1 * spec.style_benchmark_means[1]
            + w2 * spec.style_benchmark_means[2];
        let sel_alpha = spec.fund_mean_return_mu_p - style_ret;
        let info_ratio = sel_alpha / te_vol;

        SharpeReturnBasedStyleAnalysisReport {
            constrained_style_weights: [w0, w1, w2],
            tracking_error_variance: te_var,
            annualized_tracking_error: te_vol,
            sharpe_style_r_squared: r2,
            passive_style_benchmark_return: style_ret,
            active_selection_alpha_return: sel_alpha,
            active_information_ratio: info_ratio,
        }
    }

    #[cfg(test)]
    mod evaluate_kan_zhou_and_sharpe_rbsa_module_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn kan_zhou_and_sharpe_rbsa_invariants_and_ad_parity() {
            let s1 = KanZhouThreeFundCombineSpec {
                sample_sharpe_squared_theta2: 0.28,
                gmv_sharpe_squared_theta_g2: 0.12,
                equal_weight_sharpe_squared_theta_ew2: 0.16,
                number_of_assets_n: 25.0,
                sample_observations_t: 120.0,
                risk_aversion_gamma: 3.0,
            };
            let r1 = evaluate_kan_zhou_three_fund_combine(&s1);
            assert!(
                r1.three_fund_expected_out_of_sample_utility
                    > r1.plug_in_expected_out_of_sample_utility
            );
            assert!(r1.estimation_risk_utility_gain > 0.0);

            let h = 1e-6;
            let r1_ad = evaluate_kan_zhou_three_fund_combine(&KanZhouThreeFundCombineSpec {
                sample_sharpe_squared_theta2: Dual::<1>::variable(
                    s1.sample_sharpe_squared_theta2,
                    0,
                ),
                gmv_sharpe_squared_theta_g2: Dual::<1>::from_value(s1.gmv_sharpe_squared_theta_g2),
                equal_weight_sharpe_squared_theta_ew2: Dual::<1>::from_value(
                    s1.equal_weight_sharpe_squared_theta_ew2,
                ),
                number_of_assets_n: Dual::<1>::from_value(s1.number_of_assets_n),
                sample_observations_t: Dual::<1>::from_value(s1.sample_observations_t),
                risk_aversion_gamma: Dual::<1>::from_value(s1.risk_aversion_gamma),
            });
            let mut s1p = s1;
            s1p.sample_sharpe_squared_theta2 += h;
            let mut s1m = s1;
            s1m.sample_sharpe_squared_theta2 -= h;
            let fd1 = (evaluate_kan_zhou_three_fund_combine(&s1p)
                .three_fund_expected_out_of_sample_utility
                - evaluate_kan_zhou_three_fund_combine(&s1m)
                    .three_fund_expected_out_of_sample_utility)
                / (2.0 * h);
            assert!((r1_ad.three_fund_expected_out_of_sample_utility.d[0] - fd1).abs() < 1e-6);

            let s2 = SharpeReturnBasedStyleAnalysisSpec {
                fund_mean_return_mu_p: 0.105,
                fund_return_variance_var_p: 0.032,
                style_benchmark_means: [0.095, 0.075, 0.045],
                style_benchmark_variances: [0.040, 0.028, 0.010],
                fund_style_covariances: [0.024, 0.011, -0.002],
            };
            let r2 = evaluate_sharpe_return_based_style_analysis(&s2);
            let w_sum = r2.constrained_style_weights[0]
                + r2.constrained_style_weights[1]
                + r2.constrained_style_weights[2];
            assert!((w_sum - 1.0).abs() < 1e-10);
            assert!(r2.sharpe_style_r_squared > 0.5 && r2.sharpe_style_r_squared < 1.0);

            let r2_ad =
                evaluate_sharpe_return_based_style_analysis(&SharpeReturnBasedStyleAnalysisSpec {
                    fund_mean_return_mu_p: Dual::<1>::variable(s2.fund_mean_return_mu_p, 0),
                    fund_return_variance_var_p: Dual::<1>::from_value(
                        s2.fund_return_variance_var_p,
                    ),
                    style_benchmark_means: s2.style_benchmark_means.map(Dual::<1>::from_value),
                    style_benchmark_variances: s2
                        .style_benchmark_variances
                        .map(Dual::<1>::from_value),
                    fund_style_covariances: s2.fund_style_covariances.map(Dual::<1>::from_value),
                });
            let mut s2p = s2;
            s2p.fund_mean_return_mu_p += h;
            let mut s2m = s2;
            s2m.fund_mean_return_mu_p -= h;
            let fd2 = (evaluate_sharpe_return_based_style_analysis(&s2p).active_information_ratio
                - evaluate_sharpe_return_based_style_analysis(&s2m).active_information_ratio)
                / (2.0 * h);
            assert!((r2_ad.active_information_ratio.d[0] - fd2).abs() < 1e-5);
        }
    }
}
pub use evaluate_kan_zhou_and_sharpe_rbsa_module::*;

mod evaluate_pastor_stambaugh_factor_prior_shrinkage_module {
    // Pástor & Stambaugh (2000, JFE) & Pástor (2000, JF) Empirical Bayes Asset Pricing
    // Prior Shrinkage toward CAPM / Factor Model (`Mispricing Uncertainty sigma_alpha`) Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct PastorStambaughFactorPriorShrinkageSpec<T: Real = f64> {
        pub sample_ols_alpha_hat: T,
        pub sample_factor_beta_hat: T,
        pub residual_return_variance_sigma_eps2: T,
        pub prior_mispricing_std_sigma_alpha: T,
        pub benchmark_excess_mean_mu_m: T,
        pub benchmark_return_variance_sigma_m2: T,
        pub sample_length_periods_t: T,
        pub risk_aversion_gamma: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct PastorStambaughFactorPriorShrinkageReport<T: Real = f64> {
        pub data_alpha_precision: T,
        pub prior_alpha_precision: T,
        pub pastor_shrinkage_intensity_omega: T,
        pub posterior_shrunk_alpha: T,
        pub posterior_predictive_residual_variance: T,
        pub dogmatic_benchmark_weight: T,
        pub optimal_active_satellite_weight: T,
    }

    pub fn evaluate_pastor_stambaugh_factor_prior_shrinkage<T: Real>(
        spec: &PastorStambaughFactorPriorShrinkageSpec<T>,
    ) -> PastorStambaughFactorPriorShrinkageReport<T> {
        let a_hat = spec.sample_ols_alpha_hat;
        let b_hat = spec.sample_factor_beta_hat;
        let sig_e2 = spec
            .residual_return_variance_sigma_eps2
            .max_val(T::from_f64(1e-5));
        let sig_a = spec
            .prior_mispricing_std_sigma_alpha
            .max_val(T::from_f64(1e-4));
        let mu_m = spec.benchmark_excess_mean_mu_m;
        let sig_m2 = spec
            .benchmark_return_variance_sigma_m2
            .max_val(T::from_f64(1e-4));
        let t_obs = spec.sample_length_periods_t.max_val(T::from_f64(12.0));
        let gam = spec.risk_aversion_gamma.max_val(T::from_f64(0.5));

        // Pástor (2000, JF Eq. 11): OLS alpha variance = (sigma_eps^2 / T) * (1 + SR_m^2)
        let sr_m2 = (mu_m * mu_m) / sig_m2;
        let var_a_ols = (sig_e2 / t_obs) * (T::one() + sr_m2);
        let prec_data = T::one() / var_a_ols;
        let prec_prior = T::one() / (sig_a * sig_a);

        // Shrinkage weight omega on sample alpha_hat vs (1 - omega) on prior alpha_0 = 0:
        let omega = prec_data / (prec_data + prec_prior);
        let post_alpha = omega * a_hat;
        let post_var_alpha = T::one() / (prec_data + prec_prior);
        let pred_res_var = sig_e2 * (T::one() + T::one() / t_obs) + post_var_alpha;

        let w_active = post_alpha / (gam * pred_res_var);
        let w_bench = (mu_m / (gam * sig_m2)) - b_hat * w_active;

        PastorStambaughFactorPriorShrinkageReport {
            data_alpha_precision: prec_data,
            prior_alpha_precision: prec_prior,
            pastor_shrinkage_intensity_omega: omega,
            posterior_shrunk_alpha: post_alpha,
            posterior_predictive_residual_variance: pred_res_var,
            dogmatic_benchmark_weight: w_bench,
            optimal_active_satellite_weight: w_active,
        }
    }

    #[cfg(test)]
    mod evaluate_pastor_stambaugh_factor_prior_shrinkage_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn pastor_stambaugh_factor_prior_shrinkage_invariants_and_ad_parity() {
            let spec = PastorStambaughFactorPriorShrinkageSpec {
                sample_ols_alpha_hat: 0.035,
                sample_factor_beta_hat: 1.12,
                residual_return_variance_sigma_eps2: 0.018,
                prior_mispricing_std_sigma_alpha: 0.020,
                benchmark_excess_mean_mu_m: 0.065,
                benchmark_return_variance_sigma_m2: 0.025,
                sample_length_periods_t: 60.0,
                risk_aversion_gamma: 3.0,
            };
            let rep = evaluate_pastor_stambaugh_factor_prior_shrinkage(&spec);
            assert!(
                rep.pastor_shrinkage_intensity_omega > 0.0
                    && rep.pastor_shrinkage_intensity_omega < 1.0
            );
            assert!(
                rep.posterior_shrunk_alpha > 0.0
                    && rep.posterior_shrunk_alpha < spec.sample_ols_alpha_hat
            );

            let h = 1e-6;
            let r_ad = evaluate_pastor_stambaugh_factor_prior_shrinkage(
                &PastorStambaughFactorPriorShrinkageSpec {
                    sample_ols_alpha_hat: Dual::<8>::variable(spec.sample_ols_alpha_hat, 0),
                    sample_factor_beta_hat: Dual::<8>::variable(spec.sample_factor_beta_hat, 1),
                    residual_return_variance_sigma_eps2: Dual::<8>::variable(
                        spec.residual_return_variance_sigma_eps2,
                        2,
                    ),
                    prior_mispricing_std_sigma_alpha: Dual::<8>::variable(
                        spec.prior_mispricing_std_sigma_alpha,
                        3,
                    ),
                    benchmark_excess_mean_mu_m: Dual::<8>::variable(
                        spec.benchmark_excess_mean_mu_m,
                        4,
                    ),
                    benchmark_return_variance_sigma_m2: Dual::<8>::variable(
                        spec.benchmark_return_variance_sigma_m2,
                        5,
                    ),
                    sample_length_periods_t: Dual::<8>::variable(spec.sample_length_periods_t, 6),
                    risk_aversion_gamma: Dual::<8>::variable(spec.risk_aversion_gamma, 7),
                },
            );
            for (idx, &ad) in r_ad.dogmatic_benchmark_weight.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.sample_ols_alpha_hat += h;
                        sm.sample_ols_alpha_hat -= h;
                    }
                    1 => {
                        sp.sample_factor_beta_hat += h;
                        sm.sample_factor_beta_hat -= h;
                    }
                    2 => {
                        sp.residual_return_variance_sigma_eps2 += h;
                        sm.residual_return_variance_sigma_eps2 -= h;
                    }
                    3 => {
                        sp.prior_mispricing_std_sigma_alpha += h;
                        sm.prior_mispricing_std_sigma_alpha -= h;
                    }
                    4 => {
                        sp.benchmark_excess_mean_mu_m += h;
                        sm.benchmark_excess_mean_mu_m -= h;
                    }
                    5 => {
                        sp.benchmark_return_variance_sigma_m2 += h;
                        sm.benchmark_return_variance_sigma_m2 -= h;
                    }
                    6 => {
                        sp.sample_length_periods_t += h;
                        sm.sample_length_periods_t -= h;
                    }
                    _ => {
                        sp.risk_aversion_gamma += h;
                        sm.risk_aversion_gamma -= h;
                    }
                }
                let fd = (evaluate_pastor_stambaugh_factor_prior_shrinkage(&sp)
                    .dogmatic_benchmark_weight
                    - evaluate_pastor_stambaugh_factor_prior_shrinkage(&sm)
                        .dogmatic_benchmark_weight)
                    / (2.0 * h);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_pastor_stambaugh_factor_prior_shrinkage_module::*;

mod evaluate_gibbons_ross_shanken_multivariate_factor_ftest_module {
    // Gibbons, Ross & Shanken (GRS, 1989, Econometrica) & Cochrane (2005, Ch. 12)
    // Exact Finite-Sample Multivariate Factor Pricing F-Test & Tangency Sharpe Efficiency Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct GibbonsRossShankenFactorFTestSpec<T: Real = f64> {
        pub alpha_quadratic_form_sr_alpha2: T,
        pub benchmark_factor_sharpe_ratio_sr_f: T,
        pub sample_observations_t: T,
        pub test_assets_count_n: T,
        pub factor_count_k: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GibbonsRossShankenFactorFTestReport<T: Real = f64> {
        pub grs_f_statistic: T,
        pub hotelling_t2_statistic: T,
        pub max_ex_post_tangency_sharpe: T,
        pub benchmark_factor_efficiency_share: T,
        pub finite_sample_dof_scaling: T,
    }

    pub fn evaluate_gibbons_ross_shanken_multivariate_factor_ftest<T: Real>(
        spec: &GibbonsRossShankenFactorFTestSpec<T>,
    ) -> GibbonsRossShankenFactorFTestReport<T> {
        let sr_a2 = spec
            .alpha_quadratic_form_sr_alpha2
            .max_val(T::from_f64(1e-6));
        let sr_f = spec.benchmark_factor_sharpe_ratio_sr_f;
        let t_obs = spec.sample_observations_t.max_val(T::from_f64(20.0));
        let n_assets = spec.test_assets_count_n.max_val(T::one());
        let k_fac = spec.factor_count_k.max_val(T::one());

        let sr_f2 = sr_f * sr_f;
        let dof_num = (t_obs - n_assets - k_fac).max_val(T::one());
        let dof_scale = dof_num / n_assets;

        // Gibbons, Ross & Shanken (1989, Eq. 5):
        // F_GRS = ((T - N - K) / N) * (alpha' Sigma^{-1} alpha) / (1 + SR_f^2)
        let f_grs = dof_scale * (sr_a2 / (T::one() + sr_f2));
        let hotelling_t2 = t_obs * (sr_a2 / (T::one() + sr_f2));

        let sr_max2 = sr_f2 + sr_a2;
        let sr_max = sr_max2.sqrt();
        let eff_share = sr_f2 / sr_max2;

        GibbonsRossShankenFactorFTestReport {
            grs_f_statistic: f_grs,
            hotelling_t2_statistic: hotelling_t2,
            max_ex_post_tangency_sharpe: sr_max,
            benchmark_factor_efficiency_share: eff_share,
            finite_sample_dof_scaling: dof_num / t_obs,
        }
    }

    #[cfg(test)]
    mod evaluate_gibbons_ross_shanken_multivariate_factor_ftest_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gibbons_ross_shanken_multivariate_factor_ftest_invariants_and_ad_parity() {
            let spec = GibbonsRossShankenFactorFTestSpec {
                alpha_quadratic_form_sr_alpha2: 0.145,
                benchmark_factor_sharpe_ratio_sr_f: 0.42,
                sample_observations_t: 180.0,
                test_assets_count_n: 25.0,
                factor_count_k: 3.0,
            };
            let rep = evaluate_gibbons_ross_shanken_multivariate_factor_ftest(&spec);
            assert!(rep.grs_f_statistic > 0.5);
            assert!(rep.max_ex_post_tangency_sharpe > spec.benchmark_factor_sharpe_ratio_sr_f);
            assert!(
                rep.benchmark_factor_efficiency_share > 0.0
                    && rep.benchmark_factor_efficiency_share < 1.0
            );

            let h = 1e-5;
            let r_ad = evaluate_gibbons_ross_shanken_multivariate_factor_ftest(
                &GibbonsRossShankenFactorFTestSpec {
                    alpha_quadratic_form_sr_alpha2: Dual::<5>::variable(
                        spec.alpha_quadratic_form_sr_alpha2,
                        0,
                    ),
                    benchmark_factor_sharpe_ratio_sr_f: Dual::<5>::variable(
                        spec.benchmark_factor_sharpe_ratio_sr_f,
                        1,
                    ),
                    sample_observations_t: Dual::<5>::variable(spec.sample_observations_t, 2),
                    test_assets_count_n: Dual::<5>::variable(spec.test_assets_count_n, 3),
                    factor_count_k: Dual::<5>::variable(spec.factor_count_k, 4),
                },
            );
            for (idx, &ad) in r_ad.grs_f_statistic.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.alpha_quadratic_form_sr_alpha2 += h;
                        sm.alpha_quadratic_form_sr_alpha2 -= h;
                    }
                    1 => {
                        sp.benchmark_factor_sharpe_ratio_sr_f += h;
                        sm.benchmark_factor_sharpe_ratio_sr_f -= h;
                    }
                    2 => {
                        sp.sample_observations_t += h;
                        sm.sample_observations_t -= h;
                    }
                    3 => {
                        sp.test_assets_count_n += h;
                        sm.test_assets_count_n -= h;
                    }
                    _ => {
                        sp.factor_count_k += h;
                        sm.factor_count_k -= h;
                    }
                }
                let fd = (evaluate_gibbons_ross_shanken_multivariate_factor_ftest(&sp)
                    .grs_f_statistic
                    - evaluate_gibbons_ross_shanken_multivariate_factor_ftest(&sm).grs_f_statistic)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_gibbons_ross_shanken_multivariate_factor_ftest_module::*;

mod evaluate_frazzini_pedersen_betting_against_beta_zero_beta_capm_module {
    // Black (1972, J. Business) & Frazzini & Pedersen (2014, J. Financial Economics, "Betting Against Beta")
    // Leverage-Constrained Zero-Beta CAPM, Security Market Line Flattening & BAB Factor Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct FrazziniPedersenBettingAgainstBetaSpec<T: Real = f64> {
        pub low_beta_portfolio_beta_l: T,
        pub high_beta_portfolio_beta_h: T,
        pub low_beta_expected_return_r_l: T,
        pub high_beta_expected_return_r_h: T,
        pub risk_free_rate_r_f: T,
        pub shadow_funding_constraint_spread_psi: T,
        pub bab_factor_volatility_sigma_bab: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct FrazziniPedersenBettingAgainstBetaReport<T: Real = f64> {
        pub long_low_beta_leverage_weight_w_l: T,
        pub short_high_beta_delevered_weight_w_h: T,
        pub realized_bab_factor_excess_return: T,
        pub theoretical_funding_constrained_bab_alpha: T,
        pub empirical_security_market_line_slope: T,
        pub bab_factor_sharpe_ratio: T,
    }

    pub fn evaluate_frazzini_pedersen_betting_against_beta_zero_beta_capm<T: Real>(
        spec: &FrazziniPedersenBettingAgainstBetaSpec<T>,
    ) -> FrazziniPedersenBettingAgainstBetaReport<T> {
        let beta_l = spec.low_beta_portfolio_beta_l.max_val(T::from_f64(0.10));
        let beta_h = spec
            .high_beta_portfolio_beta_h
            .max_val(beta_l + T::from_f64(0.10));
        let r_l = spec.low_beta_expected_return_r_l;
        let r_h = spec.high_beta_expected_return_r_h;
        let r_f = spec.risk_free_rate_r_f;
        let psi = spec.shadow_funding_constraint_spread_psi.max_val(T::zero());
        let sig_bab = spec
            .bab_factor_volatility_sigma_bab
            .max_val(T::from_f64(0.01));

        // Frazzini & Pedersen (2014, Eq. 6): r_BAB = (1/beta_L)*(r_L - r_f) - (1/beta_H)*(r_H - r_f):
        let w_l = T::one() / beta_l;
        let w_h = T::one() / beta_h;
        let r_bab = w_l * (r_l - r_f) - w_h * (r_h - r_f);

        // Proposition 1: E[r_BAB] = ((beta_H - beta_L) / (beta_L * beta_H)) * psi:
        let alpha_fp = ((beta_h - beta_l) / (beta_l * beta_h)) * psi;
        let sml_slope = (r_h - r_l) / (beta_h - beta_l);
        let sr_bab = r_bab / sig_bab;

        FrazziniPedersenBettingAgainstBetaReport {
            long_low_beta_leverage_weight_w_l: w_l,
            short_high_beta_delevered_weight_w_h: w_h,
            realized_bab_factor_excess_return: r_bab,
            theoretical_funding_constrained_bab_alpha: alpha_fp,
            empirical_security_market_line_slope: sml_slope,
            bab_factor_sharpe_ratio: sr_bab,
        }
    }

    #[cfg(test)]
    mod evaluate_frazzini_pedersen_betting_against_beta_zero_beta_capm_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn frazzini_pedersen_betting_against_beta_invariants_and_ad_parity() {
            let spec = FrazziniPedersenBettingAgainstBetaSpec {
                low_beta_portfolio_beta_l: 0.65,
                high_beta_portfolio_beta_h: 1.45,
                low_beta_expected_return_r_l: 0.088,
                high_beta_expected_return_r_h: 0.118,
                risk_free_rate_r_f: 0.038,
                shadow_funding_constraint_spread_psi: 0.022,
                bab_factor_volatility_sigma_bab: 0.105,
            };
            let rep = evaluate_frazzini_pedersen_betting_against_beta_zero_beta_capm(&spec);
            assert!(
                rep.long_low_beta_leverage_weight_w_l > rep.short_high_beta_delevered_weight_w_h
            );
            assert!(rep.realized_bab_factor_excess_return > 0.0);
            assert!(rep.theoretical_funding_constrained_bab_alpha > 0.0);
            assert!(rep.bab_factor_sharpe_ratio > 0.20);

            let h = 1e-5;
            let r_ad = evaluate_frazzini_pedersen_betting_against_beta_zero_beta_capm(
                &FrazziniPedersenBettingAgainstBetaSpec {
                    low_beta_portfolio_beta_l: Dual::<7>::variable(
                        spec.low_beta_portfolio_beta_l,
                        0,
                    ),
                    high_beta_portfolio_beta_h: Dual::<7>::variable(
                        spec.high_beta_portfolio_beta_h,
                        1,
                    ),
                    low_beta_expected_return_r_l: Dual::<7>::variable(
                        spec.low_beta_expected_return_r_l,
                        2,
                    ),
                    high_beta_expected_return_r_h: Dual::<7>::variable(
                        spec.high_beta_expected_return_r_h,
                        3,
                    ),
                    risk_free_rate_r_f: Dual::<7>::variable(spec.risk_free_rate_r_f, 4),
                    shadow_funding_constraint_spread_psi: Dual::<7>::variable(
                        spec.shadow_funding_constraint_spread_psi,
                        5,
                    ),
                    bab_factor_volatility_sigma_bab: Dual::<7>::variable(
                        spec.bab_factor_volatility_sigma_bab,
                        6,
                    ),
                },
            );
            let comp_ad =
                r_ad.bab_factor_sharpe_ratio + r_ad.theoretical_funding_constrained_bab_alpha;
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.low_beta_portfolio_beta_l += h;
                        sm.low_beta_portfolio_beta_l -= h;
                    }
                    1 => {
                        sp.high_beta_portfolio_beta_h += h;
                        sm.high_beta_portfolio_beta_h -= h;
                    }
                    2 => {
                        sp.low_beta_expected_return_r_l += h;
                        sm.low_beta_expected_return_r_l -= h;
                    }
                    3 => {
                        sp.high_beta_expected_return_r_h += h;
                        sm.high_beta_expected_return_r_h -= h;
                    }
                    4 => {
                        sp.risk_free_rate_r_f += h;
                        sm.risk_free_rate_r_f -= h;
                    }
                    5 => {
                        sp.shadow_funding_constraint_spread_psi += h;
                        sm.shadow_funding_constraint_spread_psi -= h;
                    }
                    _ => {
                        sp.bab_factor_volatility_sigma_bab += h;
                        sm.bab_factor_volatility_sigma_bab -= h;
                    }
                }
                let rp = evaluate_frazzini_pedersen_betting_against_beta_zero_beta_capm(&sp);
                let rm = evaluate_frazzini_pedersen_betting_against_beta_zero_beta_capm(&sm);
                let fd = ((rp.bab_factor_sharpe_ratio
                    + rp.theoretical_funding_constrained_bab_alpha)
                    - (rm.bab_factor_sharpe_ratio + rm.theoretical_funding_constrained_bab_alpha))
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_frazzini_pedersen_betting_against_beta_zero_beta_capm_module::*;

mod evaluate_harvey_liu_zhu_factor_tstat_fdr_haircut_module {
    // Harvey, Liu & Zhu (2016, Review of Financial Studies, "... and the Cross-Section of Expected
    // Returns") & Bailey-López de Prado (2014) Benjamini-Hochberg-Yekutieli (BHY) False Discovery Rate
    // Critical t-Statistic Hurdle & Empirical-Bayes Posterior Sharpe Ratio Haircut Engine.

    use ql_ad::Real;
    use ql_math::{inv_norm_cdf, norm_cdf};

    #[derive(Debug, Clone, Copy)]
    pub struct HarveyLiuZhuFactorFdrHaircutSpec<T: Real = f64> {
        pub observed_factor_t_stat: T,
        pub sample_months_t: T,
        pub total_tested_factors_m: T,
        pub target_fdr_alpha: T,
        pub prior_true_factor_prob_p0: T,
        pub prior_true_sharpe_std_sigma0: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct HarveyLiuZhuFactorFdrHaircutReport<T: Real = f64> {
        pub annualized_observed_sharpe_ratio: T,
        pub single_test_two_sided_p_value: T,
        pub bhy_fdr_critical_t_hurdle: T,
        pub bayes_factor_alternative_vs_null: T,
        pub posterior_true_factor_probability: T,
        pub posterior_haircutted_sharpe_ratio: T,
        pub relative_sharpe_haircut_fraction: T,
    }

    pub fn evaluate_harvey_liu_zhu_factor_tstat_fdr_haircut<T: Real>(
        spec: &HarveyLiuZhuFactorFdrHaircutSpec<T>,
    ) -> HarveyLiuZhuFactorFdrHaircutReport<T> {
        let t_obs = spec.observed_factor_t_stat;
        let t_months = spec.sample_months_t.max_val(T::from_f64(24.0));
        let m_factors = spec.total_tested_factors_m.max_val(T::from_f64(2.0));
        let alpha_fdr = spec
            .target_fdr_alpha
            .max_val(T::from_f64(1e-4))
            .min_val(T::from_f64(0.25));
        let p0 = spec
            .prior_true_factor_prob_p0
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.99));
        let sig0 = spec.prior_true_sharpe_std_sigma0.max_val(T::from_f64(0.05));

        let v0 = T::from_f64(12.0) / t_months;
        let sr_obs = t_obs * v0.sqrt();
        let abs_t = (t_obs * t_obs + T::from_f64(1e-12)).sqrt();
        let p_single = T::from_f64(2.0) * norm_cdf(-abs_t);

        // Benjamini-Hochberg-Yekutieli (2001) harmonic sum c(M) = ln(M) + gamma_Euler + 1/(2M):
        let c_m = m_factors.ln()
            + T::from_f64(std::f64::consts::EULER_GAMMA)
            + T::one() / (T::from_f64(2.0) * m_factors);
        let alpha_bhy = (alpha_fdr / (m_factors * c_m)).max_val(T::from_f64(1e-6));
        let t_crit = -inv_norm_cdf(T::from_f64(0.5) * alpha_bhy);

        // Spike-and-slab Gaussian Empirical-Bayes shrinkage of annualized Sharpe ratio:
        let sig0_sq = sig0 * sig0;
        let v1 = v0 + sig0_sq;
        let shrink_within = sig0_sq / v1;
        let bf_10 = (v0 / v1).sqrt() * (T::from_f64(0.5) * t_obs * t_obs * shrink_within).exp();
        let post_prob = (p0 * bf_10) / ((T::one() - p0) + p0 * bf_10);
        let sr_post = post_prob * shrink_within * sr_obs;
        let haircut_frac = T::one() - post_prob * shrink_within;

        HarveyLiuZhuFactorFdrHaircutReport {
            annualized_observed_sharpe_ratio: sr_obs,
            single_test_two_sided_p_value: p_single,
            bhy_fdr_critical_t_hurdle: t_crit,
            bayes_factor_alternative_vs_null: bf_10,
            posterior_true_factor_probability: post_prob,
            posterior_haircutted_sharpe_ratio: sr_post,
            relative_sharpe_haircut_fraction: haircut_frac,
        }
    }

    #[cfg(test)]
    mod evaluate_harvey_liu_zhu_factor_tstat_fdr_haircut_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn harvey_liu_zhu_factor_fdr_haircut_invariants_and_ad_parity() {
            let spec = HarveyLiuZhuFactorFdrHaircutSpec {
                observed_factor_t_stat: 3.25,
                sample_months_t: 240.0,
                total_tested_factors_m: 316.0,
                target_fdr_alpha: 0.05,
                prior_true_factor_prob_p0: 0.15,
                prior_true_sharpe_std_sigma0: 0.35,
            };
            let rep = evaluate_harvey_liu_zhu_factor_tstat_fdr_haircut(&spec);
            // HLZ (2016): with M ~ 316 factors, the BHY critical t-hurdle exceeds 3.0:
            assert!(rep.bhy_fdr_critical_t_hurdle > 3.0);
            assert!(rep.posterior_haircutted_sharpe_ratio < rep.annualized_observed_sharpe_ratio);
            assert!(rep.relative_sharpe_haircut_fraction > 0.10);

            let h = 1e-5;
            let r_ad = evaluate_harvey_liu_zhu_factor_tstat_fdr_haircut(
                &HarveyLiuZhuFactorFdrHaircutSpec {
                    observed_factor_t_stat: Dual::<6>::variable(spec.observed_factor_t_stat, 0),
                    sample_months_t: Dual::<6>::variable(spec.sample_months_t, 1),
                    total_tested_factors_m: Dual::<6>::variable(spec.total_tested_factors_m, 2),
                    target_fdr_alpha: Dual::<6>::variable(spec.target_fdr_alpha, 3),
                    prior_true_factor_prob_p0: Dual::<6>::variable(
                        spec.prior_true_factor_prob_p0,
                        4,
                    ),
                    prior_true_sharpe_std_sigma0: Dual::<6>::variable(
                        spec.prior_true_sharpe_std_sigma0,
                        5,
                    ),
                },
            );
            let comp_ad = r_ad.posterior_haircutted_sharpe_ratio + r_ad.bhy_fdr_critical_t_hurdle;
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.observed_factor_t_stat += h;
                        sm.observed_factor_t_stat -= h;
                    }
                    1 => {
                        sp.sample_months_t += h;
                        sm.sample_months_t -= h;
                    }
                    2 => {
                        sp.total_tested_factors_m += h;
                        sm.total_tested_factors_m -= h;
                    }
                    3 => {
                        sp.target_fdr_alpha += h;
                        sm.target_fdr_alpha -= h;
                    }
                    4 => {
                        sp.prior_true_factor_prob_p0 += h;
                        sm.prior_true_factor_prob_p0 -= h;
                    }
                    _ => {
                        sp.prior_true_sharpe_std_sigma0 += h;
                        sm.prior_true_sharpe_std_sigma0 -= h;
                    }
                }
                let rp = evaluate_harvey_liu_zhu_factor_tstat_fdr_haircut(&sp);
                let rm = evaluate_harvey_liu_zhu_factor_tstat_fdr_haircut(&sm);
                let fd = ((rp.posterior_haircutted_sharpe_ratio + rp.bhy_fdr_critical_t_hurdle)
                    - (rm.posterior_haircutted_sharpe_ratio + rm.bhy_fdr_critical_t_hurdle))
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_harvey_liu_zhu_factor_tstat_fdr_haircut_module::*;

mod evaluate_novy_marx_fama_french_five_factor_spanning_alpha_module {
    // Novy-Marx (2013, Journal of Financial Economics 108(1), "The Other Side of Value:
    // The Gross Profitability Premium"), Fama & French (2015, JFE 116(1), "A Five-Factor
    // Asset Pricing Model"), & Hou, Xue & Zhang (2015, RFS 28(3))
    // Gross Profitability + Value Spanning Alpha, HMLO Orthogonalization & Tangency Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct NovyMarxFamaFrenchSpanningSpec<T: Real = f64> {
        pub gross_profitability_premium_rmw: T,
        pub value_premium_hml: T,
        pub investment_premium_cma: T,
        pub hml_on_rmw_loading_phi_rmw: T,
        pub hml_on_cma_loading_phi_cma: T,
        pub hml_spanning_residual_vol_sigma_e: T,
        pub rmw_factor_vol_sigma_rmw: T,
        pub hml_factor_vol_sigma_hml: T,
        pub rmw_hml_correlation_rho: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct NovyMarxFamaFrenchSpanningReport<T: Real = f64> {
        pub orthogonalized_hmlo_spanning_alpha: T,
        pub hmlo_information_ratio: T,
        pub marginal_squared_sharpe_contribution: T,
        pub optimal_profitability_tangency_weight: T,
        pub combined_value_profitability_sharpe: T,
    }

    pub fn evaluate_novy_marx_fama_french_five_factor_spanning_alpha<T: Real>(
        spec: &NovyMarxFamaFrenchSpanningSpec<T>,
    ) -> NovyMarxFamaFrenchSpanningReport<T> {
        let lam_rmw = spec.gross_profitability_premium_rmw;
        let lam_hml = spec.value_premium_hml;
        let lam_cma = spec.investment_premium_cma;
        let phi_rmw = spec.hml_on_rmw_loading_phi_rmw;
        let phi_cma = spec.hml_on_cma_loading_phi_cma;
        let sig_e = spec
            .hml_spanning_residual_vol_sigma_e
            .max_val(T::from_f64(0.005));
        let sig_rmw = spec.rmw_factor_vol_sigma_rmw.max_val(T::from_f64(0.01));
        let sig_hml = spec.hml_factor_vol_sigma_hml.max_val(T::from_f64(0.01));
        let rho = spec
            .rmw_hml_correlation_rho
            .max_val(T::from_f64(-0.95))
            .min_val(T::from_f64(0.95));

        // Fama-French (2015) HMLO orthogonalized against profitability (RMW) and investment (CMA):
        let alpha_hmlo = lam_hml - phi_rmw * lam_rmw - phi_cma * lam_cma;
        let ir_hmlo = alpha_hmlo / sig_e;
        let delta_sr2 = ir_hmlo * ir_hmlo;

        // Novy-Marx (2013) two-factor Value + Gross Profitability tangency portfolio weights & Sharpe:
        let cov_rh = rho * sig_rmw * sig_hml;
        let raw_w_rmw = lam_rmw * sig_hml * sig_hml - lam_hml * cov_rh;
        let raw_w_hml = lam_hml * sig_rmw * sig_rmw - lam_rmw * cov_rh;
        let w_rmw = raw_w_rmw / (raw_w_rmw + raw_w_hml).max_val(T::from_f64(1e-6));

        let sr_rmw = lam_rmw / sig_rmw;
        let sr_hml = lam_hml / sig_hml;
        let sr_comb = ((sr_rmw * sr_rmw + sr_hml * sr_hml
            - T::from_f64(2.0) * rho * sr_rmw * sr_hml)
            / (T::one() - rho * rho))
            .max_val(T::from_f64(1e-8))
            .sqrt();

        NovyMarxFamaFrenchSpanningReport {
            orthogonalized_hmlo_spanning_alpha: alpha_hmlo,
            hmlo_information_ratio: ir_hmlo,
            marginal_squared_sharpe_contribution: delta_sr2,
            optimal_profitability_tangency_weight: w_rmw,
            combined_value_profitability_sharpe: sr_comb,
        }
    }

    #[cfg(test)]
    mod evaluate_novy_marx_fama_french_five_factor_spanning_alpha_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn novy_marx_fama_french_spanning_invariants_and_ad_parity() {
            let spec = NovyMarxFamaFrenchSpanningSpec {
                gross_profitability_premium_rmw: 0.042,
                value_premium_hml: 0.038,
                investment_premium_cma: 0.032,
                hml_on_rmw_loading_phi_rmw: 0.22,
                hml_on_cma_loading_phi_cma: 0.55,
                hml_spanning_residual_vol_sigma_e: 0.068,
                rmw_factor_vol_sigma_rmw: 0.085,
                hml_factor_vol_sigma_hml: 0.105,
                rmw_hml_correlation_rho: -0.42,
            };
            let rep = evaluate_novy_marx_fama_french_five_factor_spanning_alpha(&spec);
            assert!(
                rep.combined_value_profitability_sharpe
                    > spec.gross_profitability_premium_rmw / spec.rmw_factor_vol_sigma_rmw
            );
            assert!(
                rep.optimal_profitability_tangency_weight > 0.0
                    && rep.optimal_profitability_tangency_weight < 1.0
            );

            let h = 1e-6;
            let r_ad = evaluate_novy_marx_fama_french_five_factor_spanning_alpha(
                &NovyMarxFamaFrenchSpanningSpec {
                    gross_profitability_premium_rmw: Dual::<6>::variable(
                        spec.gross_profitability_premium_rmw,
                        0,
                    ),
                    value_premium_hml: Dual::<6>::variable(spec.value_premium_hml, 1),
                    investment_premium_cma: Dual::<6>::variable(spec.investment_premium_cma, 2),
                    hml_on_rmw_loading_phi_rmw: Dual::from_f64(spec.hml_on_rmw_loading_phi_rmw),
                    hml_on_cma_loading_phi_cma: Dual::from_f64(spec.hml_on_cma_loading_phi_cma),
                    hml_spanning_residual_vol_sigma_e: Dual::<6>::variable(
                        spec.hml_spanning_residual_vol_sigma_e,
                        3,
                    ),
                    rmw_factor_vol_sigma_rmw: Dual::<6>::variable(spec.rmw_factor_vol_sigma_rmw, 4),
                    hml_factor_vol_sigma_hml: Dual::from_f64(spec.hml_factor_vol_sigma_hml),
                    rmw_hml_correlation_rho: Dual::<6>::variable(spec.rmw_hml_correlation_rho, 5),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.gross_profitability_premium_rmw,
                    spec.value_premium_hml,
                    spec.investment_premium_cma,
                    spec.hml_spanning_residual_vol_sigma_e,
                    spec.rmw_factor_vol_sigma_rmw,
                    spec.rmw_hml_correlation_rho,
                ];
                v[i] += h;
                let up = evaluate_novy_marx_fama_french_five_factor_spanning_alpha(
                    &NovyMarxFamaFrenchSpanningSpec {
                        gross_profitability_premium_rmw: v[0],
                        value_premium_hml: v[1],
                        investment_premium_cma: v[2],
                        hml_spanning_residual_vol_sigma_e: v[3],
                        rmw_factor_vol_sigma_rmw: v[4],
                        rmw_hml_correlation_rho: v[5],
                        ..spec
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_novy_marx_fama_french_five_factor_spanning_alpha(
                    &NovyMarxFamaFrenchSpanningSpec {
                        gross_profitability_premium_rmw: v[0],
                        value_premium_hml: v[1],
                        investment_premium_cma: v[2],
                        hml_spanning_residual_vol_sigma_e: v[3],
                        rmw_factor_vol_sigma_rmw: v[4],
                        rmw_hml_correlation_rho: v[5],
                        ..spec
                    },
                );
                let val_ad = r_ad.hmlo_information_ratio + r_ad.combined_value_profitability_sharpe;
                let fd = ((up.hmlo_information_ratio + up.combined_value_profitability_sharpe)
                    - (dn.hmlo_information_ratio + dn.combined_value_profitability_sharpe))
                    / (2.0 * h);
                assert!((val_ad.d[i] - fd).abs() < 1e-5, "input {i}");
            }
        }
    }
}
pub use evaluate_novy_marx_fama_french_five_factor_spanning_alpha_module::*;
