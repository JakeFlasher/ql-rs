//! Bergomi (2016) & Gatheral (2006) Two-Factor Bergomi Forward-Variance Curve,
//! Normalisation Factor alpha_theta & ATMF Skew Term-Structure Engine.

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
mod phase40_tests {
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
                flat_forward_variance_xi0: Dual::<1>::from_value(spec.flat_forward_variance_xi0),
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
        let fd = (evaluate_bergomi_two_factor_skew_term_structure(&sp).atmf_implied_vol_skew_slope
            - evaluate_bergomi_two_factor_skew_term_structure(&sm).atmf_implied_vol_skew_slope)
            / (2.0 * h);
        assert!((r_ad.atmf_implied_vol_skew_slope.d[0] - fd).abs() < 1e-6);
    }
}
