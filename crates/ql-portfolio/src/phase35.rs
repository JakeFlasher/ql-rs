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
        q[j] =
            s.prior_probabilities[j] * (lam * (s.scenario_returns[j] - s.target_view_return)).exp();
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
        p[j] = (T::one() - s.view_confidence) * s.prior_probabilities[j] + s.view_confidence * q[j];
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
mod phase35_tests {
    use super::*;
    use ql_ad::Dual;
    fn tc_spec<T: Real>(x: T) -> GrinoldQianTransferCoefficientFlAmSpec<T> {
        GrinoldQianTransferCoefficientFlAmSpec {
            information_coefficient_ic: T::from_f64(0.06),
            breadth_br: T::from_f64(100.0),
            target_active_risk: T::from_f64(0.04) + x,
            residual_volatilities: [T::from_f64(0.15), T::from_f64(0.20), T::from_f64(0.18)],
            alpha_forecasts: [T::from_f64(0.03), T::from_f64(0.02), T::from_f64(-0.025)],
            constrained_active_weights: [T::from_f64(0.12), T::from_f64(0.05), T::from_f64(-0.08)],
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
    fn phase35_portfolio_invariants_and_dual_fd() {
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
        let ad2 =
            evaluate_meucci_relative_entropy_view_pooling(&ep_spec(Dual::<1>::variable(0.0, 0)))
                .posterior_expected_return;
        assert!(ad2.d[0].abs() > 1e-6 && (ad2.d[0] - fd2).abs() < 1e-5);
    }
}
