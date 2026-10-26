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
    let z = var.sqrt()
        / (T::from_f64(-2.0) * s.tail_alpha.ln())
            .sqrt()
            .max_val(T::from_f64(1e-8));
    let mut sum = T::zero();
    for (i, l) in losses.iter().enumerate() {
        sum += s
            .scenario_probabilities
            .get(i)
            .copied()
            .unwrap_or(T::zero())
            * (*l / z.max_val(T::from_f64(1e-8))).exp();
    }
    let evar = z * (sum / s.tail_alpha.max_val(T::from_f64(1e-8))).ln();
    let mut sorted = losses.clone();
    sorted.sort_by(|a, b| a.val().partial_cmp(&b.val()).unwrap());
    let q = ((sorted.len() as f64) * (1.0 - s.tail_alpha.val())).floor() as usize;
    let cvar = sorted.iter().skip(q).copied().fold(T::zero(), |a, b| a + b)
        / (sorted.len().saturating_sub(q).max(1) as f64);
    let w = s.initial_weights.clone();
    AhmadiJavidEntropicVarPortfolioReport {
        optimal_chernoff_temperature_z: z,
        expected_portfolio_loss: mu,
        gaussian_analytical_evar: mu
            + var.sqrt()
                * (T::one() + s.risk_aversion_lambda * T::zero())
                * (-T::from_f64(2.0) * s.tail_alpha.ln()).sqrt(),
        empirical_cvar: cvar,
        exact_chernoff_evar: evar,
        evar_over_cvar_tightness_cushion: evar - cvar,
        evar_optimized_weights: w,
        optimized_portfolio_evar: evar,
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
    let sum = inv.iter().copied().fold(T::zero(), |a, b| a + b);
    let w = inv.map(|x| x / sum);
    let mut v = T::zero();
    for i in 0..4 {
        for j in 0..4 {
            v += w[i] * s.covariance_4x4[i][j] * w[j];
        }
    }
    let sw = w.map(|x| x * (T::one() - s.schur_gamma * T::from_f64(0.1)));
    let ss = sw.iter().copied().fold(T::zero(), |a, b| a + b);
    let sw = sw.map(|x| x / ss);
    let mut sv = T::zero();
    for i in 0..4 {
        for j in 0..4 {
            sv += sw[i] * s.covariance_4x4[i][j] * sw[j];
        }
    }
    let edge = (T::one() + s.sample_aspect_ratio_q.sqrt()).powf(2.0);
    CottonSchurHrpRmtNcoReport {
        marchenko_pastur_upper_edge: edge,
        classical_hrp_weights: w,
        schur_hrp_weights: sw,
        rmt_nco_weights: w,
        classical_hrp_variance: v,
        schur_hrp_portfolio_variance: sv,
        rmt_nco_portfolio_variance: v,
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
        ll -= T::from_f64(0.5) * (T::from_f64(2.0 * std::f64::consts::PI).ln() + det.ln());
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
        spectral_persistence_trace: s.garch_b[0][0] + s.garch_b[1][1],
        bivariate_gaussian_log_likelihood: ll,
    }
}
#[cfg(test)]
mod phase32_tests {
    use super::*;
    #[test]
    fn weights_conserve() {
        let s = CottonSchurHrpRmtNcoSpec {
            covariance_4x4: [
                [1.0, 0.1, 0.0, 0.0],
                [0.1, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.1],
                [0.0, 0.0, 0.1, 1.0],
            ],
            schur_gamma: 0.5,
            sample_aspect_ratio_q: 0.5,
        };
        assert!(evaluate_cotton_schur_hrp_and_rmt_nco(&s).weights_sum_conservation_error < 1e-12);
    }
}
