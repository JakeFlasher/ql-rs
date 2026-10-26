//! Brigo & Mercurio (2006) & Gregory (2020) Collateralized Cross-Currency Basis Swap (`XCCY`)
//! Mark-to-Market (`MtM`) Principal-Reset & FX-Forward Implied Basis-Spread Bootstrapper.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct BrigoXccyMtmBasisCurveSpec<T: Real = f64> {
    pub spot_fx_s0: T,
    pub domestic_csa_rate_rd: T,
    pub foreign_local_ois_rate_rf: T,
    pub fx_forward_quote_f0t: T,
    pub mtm_xccy_par_basis_spread_bm: T,
    pub fx_vol_sigma_s: T,
    pub foreign_rate_vol_sigma_rf: T,
    pub fx_rate_corr_rho: T,
    pub maturity_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct BrigoXccyMtmBasisCurveReport<T: Real = f64> {
    pub domestic_csa_df_pd: T,
    pub fx_implied_foreign_csa_df_pf: T,
    pub fx_implied_ccy_basis_spread: T,
    pub mtm_stripped_foreign_csa_zero_rate: T,
    pub constant_notional_vs_mtm_convexity_adj_bps: T,
    pub constant_notional_par_basis_spread_bc: T,
    pub covered_interest_parity_residual: T,
}

pub fn evaluate_brigo_xccy_mtm_collateralized_basis_curve<T: Real>(
    spec: &BrigoXccyMtmBasisCurveSpec<T>,
) -> BrigoXccyMtmBasisCurveReport<T> {
    let t = spec.maturity_years_t.max_val(T::from_f64(0.1));
    let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-6));
    let f0t = spec.fx_forward_quote_f0t.max_val(T::from_f64(1e-6));

    // Domestic CSA discount factor P_{d,CSA}(0,T) = exp(-r_d T)
    let pd = (-spec.domestic_csa_rate_rd * t).exp();
    // By Covered Interest Parity under domestic CSA collateral:
    // F_{FX}(0,T) = S_0 * P_{f,CSA}(0,T) / P_{d,CSA}(0,T)
    let pf_fx = pd * (f0t / s0);
    let r_f_csa_fx = -pf_fx.ln() / t;
    let fx_implied_basis = spec.foreign_local_ois_rate_rf - r_f_csa_fx;

    // MtM XCCY swap par basis spread b_M shifts the foreign collateralized curve:
    let r_f_csa_mtm = spec.foreign_local_ois_rate_rf - spec.mtm_xccy_par_basis_spread_bm;
    let pf_mtm = (-r_f_csa_mtm * t).exp();
    let f_mtm = s0 * (pf_mtm / pd);
    let cip_residual =
        (f0t * (-fx_implied_basis * t + spec.mtm_xccy_par_basis_spread_bm * t).exp() - f_mtm).abs();

    // Constant-notional vs. Mark-to-Market XCCY convexity adjustment:
    // Under correlated FX and foreign short rate, constant-notional principal lacks the
    // periodic FX reset N_{d,i} = N_f S(t_{i-1})/S(0), inducing covariance adjustment
    // Delta_b_conv = -rho_{S,r_f} * sigma_S * sigma_{r_f} * T / 2
    let conv_adj = -T::from_f64(0.5)
        * spec.fx_rate_corr_rho
        * spec.fx_vol_sigma_s
        * spec.foreign_rate_vol_sigma_rf
        * t;
    let conv_adj_bps = conv_adj * T::from_f64(1e4);
    let bc = spec.mtm_xccy_par_basis_spread_bm + conv_adj;

    BrigoXccyMtmBasisCurveReport {
        domestic_csa_df_pd: pd,
        fx_implied_foreign_csa_df_pf: pf_fx,
        fx_implied_ccy_basis_spread: fx_implied_basis,
        mtm_stripped_foreign_csa_zero_rate: r_f_csa_mtm,
        constant_notional_vs_mtm_convexity_adj_bps: conv_adj_bps,
        constant_notional_par_basis_spread_bc: bc,
        covered_interest_parity_residual: cip_residual,
    }
}

#[cfg(test)]
mod phase39_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn brigo_xccy_mtm_basis_curve_invariants_and_ad_parity() {
        let spec = BrigoXccyMtmBasisCurveSpec {
            spot_fx_s0: 1.08,
            domestic_csa_rate_rd: 0.048,
            foreign_local_ois_rate_rf: 0.032,
            fx_forward_quote_f0t: 1.08 * ((0.048 - (0.032 + 0.0025)) * 2.0_f64).exp(),
            mtm_xccy_par_basis_spread_bm: -0.0025,
            fx_vol_sigma_s: 0.09,
            foreign_rate_vol_sigma_rf: 0.008,
            fx_rate_corr_rho: -0.35,
            maturity_years_t: 2.0,
        };
        let rep = evaluate_brigo_xccy_mtm_collateralized_basis_curve(&spec);
        assert!(rep.covered_interest_parity_residual < 1e-12);
        assert!(
            (rep.fx_implied_ccy_basis_spread - spec.mtm_xccy_par_basis_spread_bm).abs() < 1e-12
        );
        assert!(rep.constant_notional_vs_mtm_convexity_adj_bps > 0.0);

        let h = 1e-6;
        let r_ad =
            evaluate_brigo_xccy_mtm_collateralized_basis_curve(&BrigoXccyMtmBasisCurveSpec {
                spot_fx_s0: Dual::<1>::from_value(spec.spot_fx_s0),
                domestic_csa_rate_rd: Dual::<1>::variable(spec.domestic_csa_rate_rd, 0),
                foreign_local_ois_rate_rf: Dual::<1>::from_value(spec.foreign_local_ois_rate_rf),
                fx_forward_quote_f0t: Dual::<1>::from_value(spec.fx_forward_quote_f0t),
                mtm_xccy_par_basis_spread_bm: Dual::<1>::from_value(
                    spec.mtm_xccy_par_basis_spread_bm,
                ),
                fx_vol_sigma_s: Dual::<1>::from_value(spec.fx_vol_sigma_s),
                foreign_rate_vol_sigma_rf: Dual::<1>::from_value(spec.foreign_rate_vol_sigma_rf),
                fx_rate_corr_rho: Dual::<1>::from_value(spec.fx_rate_corr_rho),
                maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
            });
        let mut sp = spec;
        sp.domestic_csa_rate_rd += h;
        let mut sm = spec;
        sm.domestic_csa_rate_rd -= h;
        let fd = (evaluate_brigo_xccy_mtm_collateralized_basis_curve(&sp)
            .fx_implied_foreign_csa_df_pf
            - evaluate_brigo_xccy_mtm_collateralized_basis_curve(&sm).fx_implied_foreign_csa_df_pf)
            / (2.0 * h);
        assert!((r_ad.fx_implied_foreign_csa_df_pf.d[0] - fd).abs() < 1e-7);
    }
}
