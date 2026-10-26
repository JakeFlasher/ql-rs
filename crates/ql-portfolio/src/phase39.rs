//! Engine 5A: Grinold & Kahn (2000) & Clarke–de Silva–Thorley (2002) Fundamental Law of Active
//! Management (`FLAM`), Transfer Coefficient (`TC`) & Effective Breadth Engine.
//! Engine 5B: Meucci (2005) & McNeil–Frey–Embrechts (2015) Elliptical & Cornish–Fisher
//! Expected Shortfall (`CVaR`) Euler Marginal Risk Contribution Engine.

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
mod phase39_tests {
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
        let r1_ad =
            evaluate_grinold_kahn_flam_with_transfer_coefficient(&GrinoldKahnFlamTransferSpec {
                information_coefficient_ic: Dual::<1>::variable(s1.information_coefficient_ic, 0),
                nominal_breadth_n: Dual::<1>::from_value(s1.nominal_breadth_n),
                cross_sectional_signal_corr_rho: Dual::<1>::from_value(
                    s1.cross_sectional_signal_corr_rho,
                ),
                active_weights_std_adj: s1.active_weights_std_adj.map(Dual::<1>::from_value),
                unconstrained_optimal_weights_std: s1
                    .unconstrained_optimal_weights_std
                    .map(Dual::<1>::from_value),
                active_risk_aversion_lambda: Dual::<1>::from_value(s1.active_risk_aversion_lambda),
            });
        let mut s1p = s1;
        s1p.information_coefficient_ic += h1;
        let mut s1m = s1;
        s1m.information_coefficient_ic -= h1;
        let fd1 = (evaluate_grinold_kahn_flam_with_transfer_coefficient(&s1p)
            .expected_value_added_va
            - evaluate_grinold_kahn_flam_with_transfer_coefficient(&s1m).expected_value_added_va)
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
            (r_ad.gaussian_expected_shortfall.d[0] - r2.marginal_es_contributions[0]).abs() < 1e-10
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
