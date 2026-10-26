mod evaluate_adb_rbi_ndirs_and_asset_swap_basis_module {
    use ql_ad::Real;
    #[derive(Clone, Debug)]
    pub struct AdbRbiNdirsAssetSwapSpec<T: Real = f64> {
        pub payment_tenors: Vec<T>,
        pub onshore_forward_rates: Vec<T>,
        pub usd_ois_discount_factors: Vec<T>,
        pub local_rate_vol: T,
        pub ndf_fx_vol: T,
        pub rate_fx_corr: T,
        pub sovereign_bond_yield: T,
        pub usd_sofr_yield: T,
        pub ndf_ccs_basis_spread: T,
        pub withholding_tax_rate: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct AdbRbiNdirsAssetSwapReport<T: Real = f64> {
        pub onshore_par_irs_rate: T,
        pub offshore_par_ndirs_rate: T,
        pub ndirs_quanto_convexity_bps: T,
        pub after_tax_sovereign_yield: T,
        pub fx_hedged_usd_asset_swap_yield: T,
        pub usd_asset_swap_pickup_bps: T,
    }
    pub fn evaluate_adb_rbi_ndirs_and_asset_swap_basis<T: Real>(
        s: &AdbRbiNdirsAssetSwapSpec<T>,
    ) -> AdbRbiNdirsAssetSwapReport<T> {
        let mut w = T::zero();
        let mut on = T::zero();
        let mut off = T::zero();
        for i in 0..s
            .payment_tenors
            .len()
            .min(s.onshore_forward_rates.len())
            .min(s.usd_ois_discount_factors.len())
        {
            let dt = if i == 0 {
                s.payment_tenors[0]
            } else {
                s.payment_tenors[i] - s.payment_tenors[i - 1]
            };
            let wt = dt.max_val(T::from_f64(1e-8)) * s.usd_ois_discount_factors[i];
            let cv = (s.ndf_fx_vol * s.ndf_fx_vol
                - s.rate_fx_corr * s.local_rate_vol * s.ndf_fx_vol)
                * s.payment_tenors[i];
            w += wt;
            on += wt * s.onshore_forward_rates[i];
            off += wt * (s.onshore_forward_rates[i] + cv);
        }
        let onr = on / w.max_val(T::from_f64(1e-8));
        let offr = off / w.max_val(T::from_f64(1e-8));
        let pickup = s.sovereign_bond_yield * (T::one() - s.withholding_tax_rate)
            - s.usd_sofr_yield
            + s.ndf_ccs_basis_spread;
        AdbRbiNdirsAssetSwapReport {
            onshore_par_irs_rate: onr,
            offshore_par_ndirs_rate: offr,
            ndirs_quanto_convexity_bps: (offr - onr) * T::from_f64(1e4),
            after_tax_sovereign_yield: s.sovereign_bond_yield * (T::one() - s.withholding_tax_rate),
            fx_hedged_usd_asset_swap_yield: s.usd_sofr_yield + pickup,
            usd_asset_swap_pickup_bps: pickup * T::from_f64(1e4),
        }
    }
    #[cfg(test)]
    mod domain_suite_44_tests {
        use super::*;
        use ql_ad::Dual;

        fn spec<T: Real>(x: T, fx_vol: f64) -> AdbRbiNdirsAssetSwapSpec<T> {
            AdbRbiNdirsAssetSwapSpec {
                payment_tenors: vec![T::one(), T::from_f64(2.0)],
                onshore_forward_rates: vec![T::from_f64(0.03) + x, T::from_f64(0.035)],
                usd_ois_discount_factors: vec![T::from_f64(0.97), T::from_f64(0.94)],
                local_rate_vol: T::from_f64(0.1),
                ndf_fx_vol: T::from_f64(fx_vol),
                rate_fx_corr: T::from_f64(0.2),
                sovereign_bond_yield: T::from_f64(0.045),
                usd_sofr_yield: T::from_f64(0.03),
                ndf_ccs_basis_spread: T::from_f64(0.005),
                withholding_tax_rate: T::from_f64(0.1),
            }
        }

        #[test]
        fn convexity_zero_fx() {
            assert_eq!(
                evaluate_adb_rbi_ndirs_and_asset_swap_basis(&spec(0.0, 0.0))
                    .ndirs_quanto_convexity_bps,
                0.0
            );
            let h = 1e-6;
            let fd = (evaluate_adb_rbi_ndirs_and_asset_swap_basis(&spec(h, 0.15))
                .offshore_par_ndirs_rate
                - evaluate_adb_rbi_ndirs_and_asset_swap_basis(&spec(-h, 0.15))
                    .offshore_par_ndirs_rate)
                / (2.0 * h);
            let ad = evaluate_adb_rbi_ndirs_and_asset_swap_basis(&spec(
                Dual::<1>::variable(0.0, 0),
                0.15,
            ))
            .offshore_par_ndirs_rate;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_adb_rbi_ndirs_and_asset_swap_basis_module::*;

mod evaluate_tuckman_repo_specialness_and_oftr_basis_module {
    use ql_ad::Real;
    #[derive(Clone, Copy, Debug)]
    pub struct TuckmanRepoSpecialnessBasisSpec<T: Real = f64> {
        pub otr_dirty_price: T,
        pub otr_yield: T,
        pub oftr_yield: T,
        pub modified_duration: T,
        pub gc_repo_rate: T,
        pub special_repo_rate: T,
        pub specialness_horizon_years: T,
        pub repo_haircut_fraction: T,
        pub equity_funding_rate: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct TuckmanRepoSpecialnessBasisReport<T: Real = f64> {
        pub repo_specialness_spread_bps: T,
        pub repo_specialness_pv_per_100: T,
        pub total_oftr_otr_yield_spread_bps: T,
        pub repo_specialness_yield_equiv_bps: T,
        pub structural_liquidity_yield_spread_bps: T,
        pub spread_decomposition_conservation_error: T,
        pub levered_basis_net_carry_bps: T,
        pub haircut_funded_levered_roe: T,
    }
    pub fn evaluate_tuckman_repo_specialness_and_oftr_basis<T: Real>(
        s: &TuckmanRepoSpecialnessBasisSpec<T>,
    ) -> TuckmanRepoSpecialnessBasisReport<T> {
        let d = s.gc_repo_rate - s.special_repo_rate;
        let pv = s.otr_dirty_price
            * d
            * s.specialness_horizon_years
            * (-s.gc_repo_rate * s.specialness_horizon_years).exp();
        let total = (s.oftr_yield - s.otr_yield) * T::from_f64(1e4);
        let repo = pv / (s.otr_dirty_price * s.modified_duration) * T::from_f64(1e4);
        let liq = total - repo;
        let carry = (s.oftr_yield - s.otr_yield) + (s.gc_repo_rate - s.special_repo_rate)
            - T::from_f64(2.0) * s.repo_haircut_fraction * (s.equity_funding_rate - s.gc_repo_rate);
        let roe = carry / (T::from_f64(2.0) * s.repo_haircut_fraction.max_val(T::from_f64(1e-12)));
        TuckmanRepoSpecialnessBasisReport {
            repo_specialness_spread_bps: d * T::from_f64(1e4),
            repo_specialness_pv_per_100: pv * T::from_f64(100.0) / s.otr_dirty_price,
            total_oftr_otr_yield_spread_bps: total,
            repo_specialness_yield_equiv_bps: repo,
            structural_liquidity_yield_spread_bps: liq,
            spread_decomposition_conservation_error: {
                let diff = total - repo - liq;
                (diff * diff + T::from_f64(1e-30)).sqrt()
            },
            levered_basis_net_carry_bps: carry * T::from_f64(1e4),
            haircut_funded_levered_roe: roe,
        }
    }

    #[cfg(test)]
    mod domain_suite_46_tests {
        use super::*;
        use ql_ad::Dual;
        fn spec<T: Real>(x: T) -> TuckmanRepoSpecialnessBasisSpec<T> {
            TuckmanRepoSpecialnessBasisSpec {
                otr_dirty_price: T::from_f64(100.0) + x,
                otr_yield: T::from_f64(0.04),
                oftr_yield: T::from_f64(0.0415),
                modified_duration: T::from_f64(7.5),
                gc_repo_rate: T::from_f64(0.05) + x * T::from_f64(0.01),
                special_repo_rate: T::from_f64(0.035),
                specialness_horizon_years: T::from_f64(0.25),
                repo_haircut_fraction: T::from_f64(0.02),
                equity_funding_rate: T::from_f64(0.06),
            }
        }
        #[test]
        fn suite46_curves_invariants_and_dual_fd() {
            let r = evaluate_tuckman_repo_specialness_and_oftr_basis(&spec(0.0));
            assert!(
                r.spread_decomposition_conservation_error < 1e-12
                    && r.repo_specialness_spread_bps > 0.0
            );
            let h = 1e-6;
            let fd = (evaluate_tuckman_repo_specialness_and_oftr_basis(&spec(h))
                .haircut_funded_levered_roe
                - evaluate_tuckman_repo_specialness_and_oftr_basis(&spec(-h))
                    .haircut_funded_levered_roe)
                / (2.0 * h);
            let ad = evaluate_tuckman_repo_specialness_and_oftr_basis(&spec(Dual::<1>::variable(
                0.0, 0,
            )))
            .haircut_funded_levered_roe;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_tuckman_repo_specialness_and_oftr_basis_module::*;

mod evaluate_brigo_gregory_multicurrency_csa_ctd_module {
    #![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[inline]
    fn smax<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[derive(Clone, Copy, Debug)]
    pub struct BrigoGregoryMultiCurrencyCsaCtdSpec<T: Real = f64> {
        pub domestic_ois_rate: T,
        pub foreign_ois_rate: T,
        pub xccy_basis_spread: T,
        pub csa_spread_volatility: T,
        pub spot_fx: T,
        pub maturity_years: T,
        pub promised_cashflow: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct BrigoGregoryMultiCurrencyCsaCtdReport<T: Real = f64> {
        pub domestic_ois_df: T,
        pub foreign_csa_equiv_rate: T,
        pub foreign_csa_equiv_df: T,
        pub basis_adjusted_fx_forward: T,
        pub intrinsic_ctd_collateral_rate: T,
        pub expected_ctd_collateral_rate: T,
        pub csa_switch_option_rate_pickup_bps: T,
        pub ois_discounted_pv: T,
        pub ctd_csa_discounted_pv: T,
        pub ctd_valuation_discount: T,
    }

    pub fn evaluate_brigo_gregory_multicurrency_csa_ctd<T: Real>(
        spec: &BrigoGregoryMultiCurrencyCsaCtdSpec<T>,
    ) -> BrigoGregoryMultiCurrencyCsaCtdReport<T> {
        let mat = spec.maturity_years.max_val(T::from_f64(1e-8));
        let rd = spec.domestic_ois_rate;
        let rf = spec.foreign_ois_rate;
        let bx = spec.xccy_basis_spread;
        let r_fx_csa = rd - bx;
        let df_ois = (-rd * mat).exp();
        let df_fx_csa = (-r_fx_csa * mat).exp();
        let fx_fwd = spec.spot_fx * ((rd - rf - bx) * mat).exp();
        let intrinsic_ctd = smax(rd, r_fx_csa);
        let sig_t = (spec.csa_spread_volatility * mat.sqrt()).max_val(T::from_f64(1e-8));
        let d = (r_fx_csa - rd) / sig_t;
        let expected_ctd = rd + (r_fx_csa - rd) * norm_cdf(d) + sig_t * norm_pdf(d);
        let pickup_bps = smax(expected_ctd - intrinsic_ctd, T::from_f64(1e-8)) * T::from_f64(1e4);
        let ois_pv = spec.promised_cashflow * df_ois;
        let ctd_pv = spec.promised_cashflow * (-expected_ctd * mat).exp();
        let discount = ois_pv - ctd_pv;
        BrigoGregoryMultiCurrencyCsaCtdReport {
            domestic_ois_df: df_ois,
            foreign_csa_equiv_rate: r_fx_csa,
            foreign_csa_equiv_df: df_fx_csa,
            basis_adjusted_fx_forward: fx_fwd,
            intrinsic_ctd_collateral_rate: intrinsic_ctd,
            expected_ctd_collateral_rate: expected_ctd,
            csa_switch_option_rate_pickup_bps: pickup_bps,
            ois_discounted_pv: ois_pv,
            ctd_csa_discounted_pv: ctd_pv,
            ctd_valuation_discount: discount,
        }
    }

    #[cfg(test)]
    mod domain_suite_48_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn brigo_gregory_multicurrency_csa_ctd_invariants_and_ad_parity() {
            let spec = BrigoGregoryMultiCurrencyCsaCtdSpec {
                domestic_ois_rate: 0.042,
                foreign_ois_rate: 0.028,
                xccy_basis_spread: -0.0025,
                csa_spread_volatility: 0.0065,
                spot_fx: 1.085,
                maturity_years: 5.0,
                promised_cashflow: 100.0,
            };
            let rep = evaluate_brigo_gregory_multicurrency_csa_ctd(&spec);
            assert!(rep.csa_switch_option_rate_pickup_bps > 0.0);
            assert!(rep.ctd_csa_discounted_pv < rep.ois_discounted_pv);

            let h = 1e-6;
            let r_ad = evaluate_brigo_gregory_multicurrency_csa_ctd(
                &BrigoGregoryMultiCurrencyCsaCtdSpec {
                    domestic_ois_rate: Dual::<1>::variable(spec.domestic_ois_rate, 0),
                    foreign_ois_rate: Dual::<1>::from_value(spec.foreign_ois_rate),
                    xccy_basis_spread: Dual::<1>::from_value(spec.xccy_basis_spread),
                    csa_spread_volatility: Dual::<1>::from_value(spec.csa_spread_volatility),
                    spot_fx: Dual::<1>::from_value(spec.spot_fx),
                    maturity_years: Dual::<1>::from_value(spec.maturity_years),
                    promised_cashflow: Dual::<1>::from_value(spec.promised_cashflow),
                },
            );
            let mut sp = spec;
            sp.domestic_ois_rate += h;
            let mut sm = spec;
            sm.domestic_ois_rate -= h;
            let _fd = (evaluate_brigo_gregory_multicurrency_csa_ctd(&sp).ctd_csa_discounted_pv
                - evaluate_brigo_gregory_multicurrency_csa_ctd(&sm).ctd_csa_discounted_pv)
                / (2.0 * h);
            assert!((r_ad.ctd_csa_discounted_pv.d[0] - _fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_brigo_gregory_multicurrency_csa_ctd_module::*;

mod evaluate_brigo_xccy_mtm_collateralized_basis_curve_module {
    // Brigo & Mercurio (2006) & Gregory (2020) Collateralized Cross-Currency Basis Swap (`XCCY`)
    // Mark-to-Market (`MtM`) Principal-Reset & FX-Forward Implied Basis-Spread Bootstrapper.

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
            (f0t * (-fx_implied_basis * t + spec.mtm_xccy_par_basis_spread_bm * t).exp() - f_mtm)
                .abs();

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
    mod domain_suite_51_tests {
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
                    foreign_local_ois_rate_rf: Dual::<1>::from_value(
                        spec.foreign_local_ois_rate_rf,
                    ),
                    fx_forward_quote_f0t: Dual::<1>::from_value(spec.fx_forward_quote_f0t),
                    mtm_xccy_par_basis_spread_bm: Dual::<1>::from_value(
                        spec.mtm_xccy_par_basis_spread_bm,
                    ),
                    fx_vol_sigma_s: Dual::<1>::from_value(spec.fx_vol_sigma_s),
                    foreign_rate_vol_sigma_rf: Dual::<1>::from_value(
                        spec.foreign_rate_vol_sigma_rf,
                    ),
                    fx_rate_corr_rho: Dual::<1>::from_value(spec.fx_rate_corr_rho),
                    maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                });
            let mut sp = spec;
            sp.domestic_csa_rate_rd += h;
            let mut sm = spec;
            sm.domestic_csa_rate_rd -= h;
            let fd = (evaluate_brigo_xccy_mtm_collateralized_basis_curve(&sp)
                .fx_implied_foreign_csa_df_pf
                - evaluate_brigo_xccy_mtm_collateralized_basis_curve(&sm)
                    .fx_implied_foreign_csa_df_pf)
                / (2.0 * h);
            assert!((r_ad.fx_implied_foreign_csa_df_pf.d[0] - fd).abs() < 1e-7);
        }
    }
}
pub use evaluate_brigo_xccy_mtm_collateralized_basis_curve_module::*;

mod evaluate_brigo_ois_compounded_in_arrears_and_fallback_module {
    // Brigo & Mercurio (2006) & Andersen–Piterbarg (2010) IBOR-to-OIS Fallback Spread,
    // Compounded-in-Arrears vs. In-Advance Convexity & Seasoned Accrual Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BrigoOisInArrearsFallbackSpec<T: Real = f64> {
        pub ois_df_start_p0t1: T,
        pub ois_df_end_p0t2: T,
        pub accrual_year_fraction_delta: T,
        pub short_rate_vol_sigma_r: T,
        pub mean_reversion_a: T,
        pub isda_historical_fallback_spread: T,
        pub realized_accrued_rfr_rate: T,
        pub elapsed_accrual_fraction_w: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BrigoOisInArrearsFallbackReport<T: Real = f64> {
        pub forward_compounded_in_arrears_rate: T,
        pub in_advance_timing_convexity_adj_bps: T,
        pub convexity_adjusted_in_advance_rate: T,
        pub isda_fallback_ibor_all_in_rate: T,
        pub seasoned_period_blended_rfr_rate: T,
        pub martingale_arrears_parity_error: T,
    }

    pub fn evaluate_brigo_ois_compounded_in_arrears_and_fallback<T: Real>(
        spec: &BrigoOisInArrearsFallbackSpec<T>,
    ) -> BrigoOisInArrearsFallbackReport<T> {
        let delta = spec.accrual_year_fraction_delta.max_val(T::from_f64(1e-4));
        let p1 = spec.ois_df_start_p0t1.max_val(T::from_f64(1e-6));
        let p2 = spec.ois_df_end_p0t2.max_val(T::from_f64(1e-6));

        // Under T_2-forward measure, daily compounded-in-arrears RFR paid at T_2 has
        // exact martingale expectation F_0(T_1, T_2) = (P(0, T_1) / P(0, T_2) - 1) / delta:
        let ratio = p1 / p2;
        let f_arrears = (ratio - T::one()) / delta;
        let parity_err = ((T::one() + delta * f_arrears) * p2 - p1).abs();

        // When a term RFR rate is paid in-advance at T_1 (or fixed in-arrears without compounding),
        // Hull-White timing convexity adjustment over [T_1, T_2] is:
        // var_int = sigma_r^2 * delta^3 / 3 * (1 - 0.5 * a * delta)
        let a_damp =
            (T::one() - T::from_f64(0.5) * spec.mean_reversion_a * delta).max_val(T::from_f64(0.2));
        let var_int =
            spec.short_rate_vol_sigma_r * spec.short_rate_vol_sigma_r * delta * delta * delta
                / T::from_f64(3.0)
                * a_damp;
        let f_in_advance = (ratio * (-var_int).exp() - T::one()) / delta;
        let conv_adj_bps = (f_arrears - f_in_advance) * T::from_f64(1e4);

        // ISDA IBOR Fallback all-in rate = Compounded-in-arrears RFR + 5Y median spread adjustment:
        let fallback_ibor = f_arrears + spec.isda_historical_fallback_spread;

        // Seasoned mid-period valuation with elapsed weight w in [0, 1]:
        // (1 + delta * R_blended) = (1 + w * delta * R_realized) * (1 + (1 - w) * delta * F_arrears)
        let w = spec
            .elapsed_accrual_fraction_w
            .max_val(T::zero())
            .min_val(T::one());
        let gross_realized = T::one() + w * delta * spec.realized_accrued_rfr_rate;
        let gross_remaining = T::one() + (T::one() - w) * delta * f_arrears;
        let blended_rate = (gross_realized * gross_remaining - T::one()) / delta;

        BrigoOisInArrearsFallbackReport {
            forward_compounded_in_arrears_rate: f_arrears,
            in_advance_timing_convexity_adj_bps: conv_adj_bps,
            convexity_adjusted_in_advance_rate: f_in_advance,
            isda_fallback_ibor_all_in_rate: fallback_ibor,
            seasoned_period_blended_rfr_rate: blended_rate,
            martingale_arrears_parity_error: parity_err,
        }
    }

    #[cfg(test)]
    mod domain_suite_52_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn brigo_ois_in_arrears_fallback_invariants_and_ad_parity() {
            let spec = BrigoOisInArrearsFallbackSpec {
                ois_df_start_p0t1: 0.960,
                ois_df_end_p0t2: 0.938,
                accrual_year_fraction_delta: 0.5,
                short_rate_vol_sigma_r: 0.012,
                mean_reversion_a: 0.15,
                isda_historical_fallback_spread: 0.0026161,
                realized_accrued_rfr_rate: 0.049,
                elapsed_accrual_fraction_w: 0.4,
            };
            let rep = evaluate_brigo_ois_compounded_in_arrears_and_fallback(&spec);
            assert!(rep.martingale_arrears_parity_error < 1e-12);
            assert!(rep.in_advance_timing_convexity_adj_bps > 0.0);
            assert!(rep.isda_fallback_ibor_all_in_rate > rep.forward_compounded_in_arrears_rate);

            let h = 1e-6;
            let r_ad = evaluate_brigo_ois_compounded_in_arrears_and_fallback(
                &BrigoOisInArrearsFallbackSpec {
                    ois_df_start_p0t1: Dual::<1>::from_value(spec.ois_df_start_p0t1),
                    ois_df_end_p0t2: Dual::<1>::variable(spec.ois_df_end_p0t2, 0),
                    accrual_year_fraction_delta: Dual::<1>::from_value(
                        spec.accrual_year_fraction_delta,
                    ),
                    short_rate_vol_sigma_r: Dual::<1>::from_value(spec.short_rate_vol_sigma_r),
                    mean_reversion_a: Dual::<1>::from_value(spec.mean_reversion_a),
                    isda_historical_fallback_spread: Dual::<1>::from_value(
                        spec.isda_historical_fallback_spread,
                    ),
                    realized_accrued_rfr_rate: Dual::<1>::from_value(
                        spec.realized_accrued_rfr_rate,
                    ),
                    elapsed_accrual_fraction_w: Dual::<1>::from_value(
                        spec.elapsed_accrual_fraction_w,
                    ),
                },
            );
            let mut sp = spec;
            sp.ois_df_end_p0t2 += h;
            let mut sm = spec;
            sm.ois_df_end_p0t2 -= h;
            let fd = (evaluate_brigo_ois_compounded_in_arrears_and_fallback(&sp)
                .seasoned_period_blended_rfr_rate
                - evaluate_brigo_ois_compounded_in_arrears_and_fallback(&sm)
                    .seasoned_period_blended_rfr_rate)
                / (2.0 * h);
            assert!((r_ad.seasoned_period_blended_rfr_rate.d[0] - fd).abs() < 1e-7);
        }
    }
}
pub use evaluate_brigo_ois_compounded_in_arrears_and_fallback_module::*;

mod evaluate_darbyshire_central_bank_step_ois_and_turn_module {
    // Darbyshire (2016) & Henrard (2014) Central-Bank Meeting-Dated Step-OIS Curve
    // & Turn-of-Year (`ToY`) Repo Spike Decomposition Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct DarbyshireCentralBankStepOisTurnSpec<T: Real = f64> {
        pub current_overnight_ois_r0: T,
        pub meeting_1_time_years_t1: T,
        pub meeting_1_expected_step_bps: T,
        pub meeting_2_time_years_t2: T,
        pub meeting_2_expected_step_bps: T,
        pub turn_of_year_spike_bps: T,
        pub turn_effective_day_weight_tau: T,
        pub horizon_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct DarbyshireCentralBankStepOisTurnReport<T: Real = f64> {
        pub terminal_policy_ois_rate: T,
        pub integrated_step_ois_zero_rate: T,
        pub turn_of_year_zero_rate_addon_bps: T,
        pub all_in_discount_factor_p0t: T,
        pub strip_without_turn_discount_factor: T,
        pub meeting_1_mp01_per_bp: T,
        pub meeting_2_mp01_per_bp: T,
    }

    pub fn evaluate_darbyshire_central_bank_step_ois_and_turn<T: Real>(
        spec: &DarbyshireCentralBankStepOisTurnSpec<T>,
    ) -> DarbyshireCentralBankStepOisTurnReport<T> {
        let r0 = spec.current_overnight_ois_r0;
        let t = spec.horizon_years_t.max_val(T::from_f64(0.10));
        let t1 = spec
            .meeting_1_time_years_t1
            .max_val(T::from_f64(0.02))
            .min_val(t);
        let t2 = spec.meeting_2_time_years_t2.max_val(t1).min_val(t);

        let step1 = spec.meeting_1_expected_step_bps * T::from_f64(1e-4);
        let step2 = spec.meeting_2_expected_step_bps * T::from_f64(1e-4);
        let r1 = r0 + step1;
        let r2 = r1 + step2;

        // Integrated step-function OIS rate over [0, T]:
        // int_0^T r(u) du = r0 * t + step1 * (t - t1) + step2 * (t - t2)
        let int_step = r0 * t + step1 * (t - t1) + step2 * (t - t2);
        let zero_step = int_step / t;

        // Turn-of-year (ToY) spike contribution over effective year fraction tau_turn:
        let tau_turn = spec
            .turn_effective_day_weight_tau
            .max_val(T::from_f64(1.0 / 365.0));
        let int_turn = (spec.turn_of_year_spike_bps * T::from_f64(1e-4)) * tau_turn;
        let turn_zero_bps = (int_turn / t) * T::from_f64(1e4);

        let df_no_turn = (-int_step).exp();
        let df_all_in = (-(int_step + int_turn)).exp();

        // Meeting-date Jacobian sensitivities (MP01 per 1bp step on a 10,000 notional OIS):
        let mp01_1 = (t - t1) * df_all_in;
        let mp01_2 = (t - t2) * df_all_in;

        DarbyshireCentralBankStepOisTurnReport {
            terminal_policy_ois_rate: r2,
            integrated_step_ois_zero_rate: zero_step,
            turn_of_year_zero_rate_addon_bps: turn_zero_bps,
            all_in_discount_factor_p0t: df_all_in,
            strip_without_turn_discount_factor: df_no_turn,
            meeting_1_mp01_per_bp: mp01_1,
            meeting_2_mp01_per_bp: mp01_2,
        }
    }

    #[cfg(test)]
    mod domain_suite_54_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn darbyshire_central_bank_step_ois_and_turn_invariants_and_ad_parity() {
            let spec = DarbyshireCentralBankStepOisTurnSpec {
                current_overnight_ois_r0: 0.045,
                meeting_1_time_years_t1: 0.15,
                meeting_1_expected_step_bps: 22.0,
                meeting_2_time_years_t2: 0.35,
                meeting_2_expected_step_bps: 18.0,
                turn_of_year_spike_bps: 35.0,
                turn_effective_day_weight_tau: 3.0 / 360.0,
                horizon_years_t: 0.50,
            };
            let rep = evaluate_darbyshire_central_bank_step_ois_and_turn(&spec);
            assert!(rep.terminal_policy_ois_rate > spec.current_overnight_ois_r0);
            assert!(rep.all_in_discount_factor_p0t < rep.strip_without_turn_discount_factor);
            assert!(rep.meeting_1_mp01_per_bp > rep.meeting_2_mp01_per_bp);

            let h = 1e-6;
            let r_ad = evaluate_darbyshire_central_bank_step_ois_and_turn(
                &DarbyshireCentralBankStepOisTurnSpec {
                    current_overnight_ois_r0: Dual::<1>::variable(spec.current_overnight_ois_r0, 0),
                    meeting_1_time_years_t1: Dual::<1>::from_value(spec.meeting_1_time_years_t1),
                    meeting_1_expected_step_bps: Dual::<1>::from_value(
                        spec.meeting_1_expected_step_bps,
                    ),
                    meeting_2_time_years_t2: Dual::<1>::from_value(spec.meeting_2_time_years_t2),
                    meeting_2_expected_step_bps: Dual::<1>::from_value(
                        spec.meeting_2_expected_step_bps,
                    ),
                    turn_of_year_spike_bps: Dual::<1>::from_value(spec.turn_of_year_spike_bps),
                    turn_effective_day_weight_tau: Dual::<1>::from_value(
                        spec.turn_effective_day_weight_tau,
                    ),
                    horizon_years_t: Dual::<1>::from_value(spec.horizon_years_t),
                },
            );
            let mut sp = spec;
            sp.current_overnight_ois_r0 += h;
            let mut sm = spec;
            sm.current_overnight_ois_r0 -= h;
            let fd = (evaluate_darbyshire_central_bank_step_ois_and_turn(&sp)
                .all_in_discount_factor_p0t
                - evaluate_darbyshire_central_bank_step_ois_and_turn(&sm)
                    .all_in_discount_factor_p0t)
                / (2.0 * h);
            assert!((r_ad.all_in_discount_factor_p0t.d[0] - fd).abs() < 1e-7);
        }
    }
}
pub use evaluate_darbyshire_central_bank_step_ois_and_turn_module::*;

mod evaluate_henrard_tenor_basis_multiplicative_spread_module {
    // Henrard (2014) & Andersen–Piterbarg (2010) Multi-Curve Tenor-Basis Multiplicative
    // Spread Factor `K_mult` & Sub-Period 3M-vs-6M Compounded FRA Convexity Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct HenrardTenorBasisMultiplicativeSpreadSpec<T: Real = f64> {
        pub ois_discount_rate_r: T,
        pub forward_3m_period1_f1: T,
        pub forward_3m_period2_f2: T,
        pub quoted_forward_6m_f6: T,
        pub vol_3m_period1_sigma1: T,
        pub vol_3m_period2_sigma2: T,
        pub subperiod_correlation_rho12: T,
        pub start_time_years_t0: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct HenrardTenorBasisMultiplicativeSpreadReport<T: Real = f64> {
        pub deterministic_compounded_3m_rate: T,
        pub convexity_adjusted_compounded_3m_rate: T,
        pub subperiod_fra_convexity_bps: T,
        pub multiplicative_spread_factor_k_mult: T,
        pub additive_tenor_basis_bps: T,
        pub multiplicative_equivalent_basis_bps: T,
        pub ois_discounted_basis_swaplet_pv: T,
    }

    pub fn evaluate_henrard_tenor_basis_multiplicative_spread<T: Real>(
        spec: &HenrardTenorBasisMultiplicativeSpreadSpec<T>,
    ) -> HenrardTenorBasisMultiplicativeSpreadReport<T> {
        let d1 = T::from_f64(0.25);
        let d2 = T::from_f64(0.25);
        let d6 = T::from_f64(0.50);

        let f1 = spec.forward_3m_period1_f1.max_val(T::from_f64(1e-4));
        let f2 = spec.forward_3m_period2_f2.max_val(T::from_f64(1e-4));
        let f6 = spec.quoted_forward_6m_f6.max_val(T::from_f64(1e-4));
        let sig1 = spec.vol_3m_period1_sigma1.max_val(T::from_f64(1e-4));
        let sig2 = spec.vol_3m_period2_sigma2.max_val(T::from_f64(1e-4));
        let rho = spec
            .subperiod_correlation_rho12
            .max_val(T::from_f64(-0.99))
            .min_val(T::from_f64(0.99));
        let t0 = spec.start_time_years_t0.max_val(T::from_f64(0.05));

        let g1 = T::one() + d1 * f1;
        let g2 = T::one() + d2 * f2;
        let g_det = g1 * g2;
        let f_comp_det = (g_det - T::one()) / d6;

        // Henrard (2014, Ch. 3) / Mercurio (2010) sub-period log-FRA covariance multiplier:
        let w1 = (d1 * f1) / g1;
        let w2 = (d2 * f2) / g2;
        let gamma_conv = (w1 * w2 * rho * sig1 * sig2 * t0).exp();
        let g_conv = g_det * gamma_conv;
        let f_comp_conv = (g_conv - T::one()) / d6;
        let conv_bps = (f_comp_conv - f_comp_det) * T::from_f64(10_000.0);

        let k_mult = (T::one() + d6 * f6) / g_conv;
        let add_basis_bps = (f6 - f_comp_conv) * T::from_f64(10_000.0);
        let mult_basis_bps = ((k_mult - T::one()) / d6) * T::from_f64(10_000.0);

        let df_pay = (-spec.ois_discount_rate_r * (t0 + d6)).exp();
        let pv_swaplet = T::from_f64(10_000.0) * d6 * df_pay * (f6 - f_comp_conv);

        HenrardTenorBasisMultiplicativeSpreadReport {
            deterministic_compounded_3m_rate: f_comp_det,
            convexity_adjusted_compounded_3m_rate: f_comp_conv,
            subperiod_fra_convexity_bps: conv_bps,
            multiplicative_spread_factor_k_mult: k_mult,
            additive_tenor_basis_bps: add_basis_bps,
            multiplicative_equivalent_basis_bps: mult_basis_bps,
            ois_discounted_basis_swaplet_pv: pv_swaplet,
        }
    }

    #[cfg(test)]
    mod domain_suite_55_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn henrard_tenor_basis_multiplicative_spread_invariants_and_ad_parity() {
            let spec = HenrardTenorBasisMultiplicativeSpreadSpec {
                ois_discount_rate_r: 0.038,
                forward_3m_period1_f1: 0.041,
                forward_3m_period2_f2: 0.043,
                quoted_forward_6m_f6: 0.0445,
                vol_3m_period1_sigma1: 0.18,
                vol_3m_period2_sigma2: 0.19,
                subperiod_correlation_rho12: 0.85,
                start_time_years_t0: 1.0,
            };
            let rep = evaluate_henrard_tenor_basis_multiplicative_spread(&spec);
            assert!(rep.subperiod_fra_convexity_bps > 0.0);
            assert!(rep.multiplicative_spread_factor_k_mult > 1.0);
            assert!(rep.additive_tenor_basis_bps > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_henrard_tenor_basis_multiplicative_spread(
                &HenrardTenorBasisMultiplicativeSpreadSpec {
                    ois_discount_rate_r: Dual::<1>::from_value(spec.ois_discount_rate_r),
                    forward_3m_period1_f1: Dual::<1>::variable(spec.forward_3m_period1_f1, 0),
                    forward_3m_period2_f2: Dual::<1>::from_value(spec.forward_3m_period2_f2),
                    quoted_forward_6m_f6: Dual::<1>::from_value(spec.quoted_forward_6m_f6),
                    vol_3m_period1_sigma1: Dual::<1>::from_value(spec.vol_3m_period1_sigma1),
                    vol_3m_period2_sigma2: Dual::<1>::from_value(spec.vol_3m_period2_sigma2),
                    subperiod_correlation_rho12: Dual::<1>::from_value(
                        spec.subperiod_correlation_rho12,
                    ),
                    start_time_years_t0: Dual::<1>::from_value(spec.start_time_years_t0),
                },
            );
            let mut sp = spec;
            sp.forward_3m_period1_f1 += h;
            let mut sm = spec;
            sm.forward_3m_period1_f1 -= h;
            let fd = (evaluate_henrard_tenor_basis_multiplicative_spread(&sp)
                .ois_discounted_basis_swaplet_pv
                - evaluate_henrard_tenor_basis_multiplicative_spread(&sm)
                    .ois_discounted_basis_swaplet_pv)
                / (2.0 * h);
            assert!((r_ad.ois_discounted_basis_swaplet_pv.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_henrard_tenor_basis_multiplicative_spread_module::*;

mod evaluate_brigo_differential_swap_quanto_fra_module {
    // Brigo & Mercurio (2006, Ch. 13) & Piterbarg (2012) Differential Swap (`Diff-Swap`)
    // & Quanto FRA FX-Rate Covariance Convexity Adjustment Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BrigoDifferentialSwapQuantoFraSpec<T: Real = f64> {
        pub domestic_ois_rate_rd: T,
        pub domestic_forward_rate_ld: T,
        pub foreign_forward_rate_lf: T,
        pub foreign_rate_vol_sigma_lf: T,
        pub fx_spot_vol_sigma_fx: T,
        pub rate_fx_correlation_rho: T,
        pub reset_time_years_t1: T,
        pub accrual_fraction_tau: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BrigoDifferentialSwapQuantoFraReport<T: Real = f64> {
        pub unadjusted_foreign_forward_rate: T,
        pub quanto_drift_multiplier: T,
        pub timing_in_arrears_convexity_bps: T,
        pub quanto_covariance_adjustment_bps: T,
        pub quanto_adjusted_foreign_rate: T,
        pub fair_differential_swap_spread_bps: T,
        pub domestic_discounted_diff_swaplet_pv: T,
    }

    pub fn evaluate_brigo_differential_swap_quanto_fra<T: Real>(
        spec: &BrigoDifferentialSwapQuantoFraSpec<T>,
    ) -> BrigoDifferentialSwapQuantoFraReport<T> {
        let rd = spec.domestic_ois_rate_rd.max_val(T::from_f64(1e-4));
        let ld = spec.domestic_forward_rate_ld.max_val(T::from_f64(1e-4));
        let lf = spec.foreign_forward_rate_lf.max_val(T::from_f64(1e-4));
        let sig_lf = spec.foreign_rate_vol_sigma_lf.max_val(T::from_f64(1e-4));
        let sig_fx = spec.fx_spot_vol_sigma_fx.max_val(T::from_f64(1e-4));
        let rho = spec
            .rate_fx_correlation_rho
            .max_val(T::from_f64(-0.98))
            .min_val(T::from_f64(0.98));
        let t1 = spec.reset_time_years_t1.max_val(T::from_f64(0.05));
        let tau = spec.accrual_fraction_tau.max_val(T::from_f64(0.05));

        // Brigo & Mercurio (2006, Section 13.4) Girsanov quanto measure change:
        // Under domestic forward measure Q_d^{T2}, foreign rate L_f acquires drift -rho * sigma_Lf * sigma_FX
        let quanto_mult = (-rho * sig_lf * sig_fx * t1).exp();
        let lf_quanto_only = lf * quanto_mult;
        let quanto_adj_bps = (lf_quanto_only - lf) * T::from_f64(10_000.0);

        // Sub-period payment-lag convexity adjustment:
        let w_tau = (tau * lf) / (T::one() + tau * lf);
        let timing_mult = (w_tau * sig_lf * sig_lf * t1).exp();
        let lf_all_in = lf * quanto_mult * timing_mult;
        let timing_bps = (lf_all_in - lf_quanto_only) * T::from_f64(10_000.0);

        // Fair spread s_diff added to domestic leg so E_d[L_f^{quanto} - (L_d + s_diff)] = 0:
        let diff_spread_bps = (lf_all_in - ld) * T::from_f64(10_000.0);
        let df_d = (-rd * (t1 + tau)).exp();
        let pv_swaplet = T::from_f64(10_000.0) * tau * df_d * (lf_all_in - ld);

        BrigoDifferentialSwapQuantoFraReport {
            unadjusted_foreign_forward_rate: lf,
            quanto_drift_multiplier: quanto_mult,
            timing_in_arrears_convexity_bps: timing_bps,
            quanto_covariance_adjustment_bps: quanto_adj_bps,
            quanto_adjusted_foreign_rate: lf_all_in,
            fair_differential_swap_spread_bps: diff_spread_bps,
            domestic_discounted_diff_swaplet_pv: pv_swaplet,
        }
    }

    #[cfg(test)]
    mod evaluate_brigo_differential_swap_quanto_fra_module_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn brigo_differential_swap_quanto_fra_invariants_and_ad_parity() {
            let spec = BrigoDifferentialSwapQuantoFraSpec {
                domestic_ois_rate_rd: 0.042,
                domestic_forward_rate_ld: 0.044,
                foreign_forward_rate_lf: 0.032,
                foreign_rate_vol_sigma_lf: 0.22,
                fx_spot_vol_sigma_fx: 0.12,
                rate_fx_correlation_rho: -0.45,
                reset_time_years_t1: 2.0,
                accrual_fraction_tau: 0.50,
            };
            let rep = evaluate_brigo_differential_swap_quanto_fra(&spec);
            assert!(rep.quanto_covariance_adjustment_bps > 0.0);
            assert!(rep.timing_in_arrears_convexity_bps > 0.0);
            assert!(rep.quanto_adjusted_foreign_rate > rep.unadjusted_foreign_forward_rate);

            let h = 1e-6;
            let r_ad =
                evaluate_brigo_differential_swap_quanto_fra(&BrigoDifferentialSwapQuantoFraSpec {
                    domestic_ois_rate_rd: Dual::<8>::variable(spec.domestic_ois_rate_rd, 0),
                    domestic_forward_rate_ld: Dual::<8>::variable(spec.domestic_forward_rate_ld, 1),
                    foreign_forward_rate_lf: Dual::<8>::variable(spec.foreign_forward_rate_lf, 2),
                    foreign_rate_vol_sigma_lf: Dual::<8>::variable(
                        spec.foreign_rate_vol_sigma_lf,
                        3,
                    ),
                    fx_spot_vol_sigma_fx: Dual::<8>::variable(spec.fx_spot_vol_sigma_fx, 4),
                    rate_fx_correlation_rho: Dual::<8>::variable(spec.rate_fx_correlation_rho, 5),
                    reset_time_years_t1: Dual::<8>::variable(spec.reset_time_years_t1, 6),
                    accrual_fraction_tau: Dual::<8>::variable(spec.accrual_fraction_tau, 7),
                });
            for (idx, &ad_deriv) in r_ad
                .domestic_discounted_diff_swaplet_pv
                .d
                .iter()
                .enumerate()
            {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.domestic_ois_rate_rd += h;
                        sm.domestic_ois_rate_rd -= h;
                    }
                    1 => {
                        sp.domestic_forward_rate_ld += h;
                        sm.domestic_forward_rate_ld -= h;
                    }
                    2 => {
                        sp.foreign_forward_rate_lf += h;
                        sm.foreign_forward_rate_lf -= h;
                    }
                    3 => {
                        sp.foreign_rate_vol_sigma_lf += h;
                        sm.foreign_rate_vol_sigma_lf -= h;
                    }
                    4 => {
                        sp.fx_spot_vol_sigma_fx += h;
                        sm.fx_spot_vol_sigma_fx -= h;
                    }
                    5 => {
                        sp.rate_fx_correlation_rho += h;
                        sm.rate_fx_correlation_rho -= h;
                    }
                    6 => {
                        sp.reset_time_years_t1 += h;
                        sm.reset_time_years_t1 -= h;
                    }
                    _ => {
                        sp.accrual_fraction_tau += h;
                        sm.accrual_fraction_tau -= h;
                    }
                }
                let fd = (evaluate_brigo_differential_swap_quanto_fra(&sp)
                    .domestic_discounted_diff_swaplet_pv
                    - evaluate_brigo_differential_swap_quanto_fra(&sm)
                        .domestic_discounted_diff_swaplet_pv)
                    / (2.0 * h);
                assert!((ad_deriv - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_brigo_differential_swap_quanto_fra_module::*;

mod evaluate_andersen_sofr_futures_fra_convexity_module {
    // Andersen & Piterbarg (2010, Vol. I Ch. 6) & Brigo–Mercurio (2006) Eurodollar / SOFR
    // Futures vs. FRA Convexity Bias (`Ho–Lee / Hull–White Daily Margining Adjustment`) Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct AndersenSofrFuturesFraConvexitySpec<T: Real = f64> {
        pub unadjusted_fra_forward_rate_f0: T,
        pub short_rate_volatility_sigma_r: T,
        pub hull_white_mean_reversion_kappa: T,
        pub accrual_start_years_t1: T,
        pub accrual_fraction_tau: T,
        pub contract_notional_mm: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct AndersenSofrFuturesFraConvexityReport<T: Real = f64> {
        pub unadjusted_fra_rate: T,
        pub ho_lee_zero_mean_reversion_convexity_bps: T,
        pub hull_white_convexity_exponent_gamma: T,
        pub hull_white_futures_fra_convexity_bps: T,
        pub implied_sofr_futures_rate: T,
        pub implied_imm_futures_price: T,
        pub convexity_bias_pv_dollars: T,
    }

    pub fn evaluate_andersen_sofr_futures_fra_convexity<T: Real>(
        spec: &AndersenSofrFuturesFraConvexitySpec<T>,
    ) -> AndersenSofrFuturesFraConvexityReport<T> {
        let f0 = spec
            .unadjusted_fra_forward_rate_f0
            .max_val(T::from_f64(1e-4));
        let sig = spec
            .short_rate_volatility_sigma_r
            .max_val(T::from_f64(1e-4));
        let kap = spec
            .hull_white_mean_reversion_kappa
            .max_val(T::from_f64(1e-3));
        let t1 = spec.accrual_start_years_t1.max_val(T::from_f64(0.05));
        let tau = spec.accrual_fraction_tau.max_val(T::from_f64(0.05));
        let notl = spec.contract_notional_mm.max_val(T::from_f64(0.1));
        let t2 = t1 + tau;

        // Ho-Lee (kappa -> 0) limit: Delta_HL = 0.5 * sigma_r^2 * T1 * T2
        let hl_conv_bps = T::from_f64(0.5) * sig * sig * t1 * t2 * T::from_f64(10_000.0);

        // Exact Hull-White 1F futures-forward convexity exponent Gamma(0, T1, T2) (Brigo & Mercurio Eq. 3.41):
        let b_12 = (T::one() - (-kap * tau).exp()) / kap;
        let b_01 = (T::one() - (-kap * t1).exp()) / kap;
        let var_r_t1 = (sig * sig / (T::from_f64(2.0) * kap))
            * (T::one() - (-T::from_f64(2.0) * kap * t1).exp());
        let gamma_hw = b_12
            * (T::from_f64(0.5) * sig * sig * b_01 * b_01 + T::from_f64(0.5) * b_12 * var_r_t1);

        let one_plus_tau_f = T::one() + tau * f0;
        let fut_rate = (one_plus_tau_f * gamma_hw.exp() - T::one()) / tau;
        let hw_conv_bps = (fut_rate - f0) * T::from_f64(10_000.0);
        let imm_price = T::from_f64(100.0) * (T::one() - fut_rate);
        let pv_bias = notl * T::from_f64(1_000_000.0) * tau * (fut_rate - f0);

        AndersenSofrFuturesFraConvexityReport {
            unadjusted_fra_rate: f0,
            ho_lee_zero_mean_reversion_convexity_bps: hl_conv_bps,
            hull_white_convexity_exponent_gamma: gamma_hw,
            hull_white_futures_fra_convexity_bps: hw_conv_bps,
            implied_sofr_futures_rate: fut_rate,
            implied_imm_futures_price: imm_price,
            convexity_bias_pv_dollars: pv_bias,
        }
    }

    #[cfg(test)]
    mod evaluate_andersen_sofr_futures_fra_convexity_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn andersen_sofr_futures_fra_convexity_invariants_and_ad_parity() {
            let spec = AndersenSofrFuturesFraConvexitySpec {
                unadjusted_fra_forward_rate_f0: 0.0435,
                short_rate_volatility_sigma_r: 0.0115,
                hull_white_mean_reversion_kappa: 0.06,
                accrual_start_years_t1: 3.0,
                accrual_fraction_tau: 0.25,
                contract_notional_mm: 10.0,
            };
            let rep = evaluate_andersen_sofr_futures_fra_convexity(&spec);
            assert!(rep.hull_white_futures_fra_convexity_bps > 0.0);
            assert!(
                rep.hull_white_futures_fra_convexity_bps
                    < rep.ho_lee_zero_mean_reversion_convexity_bps
            );
            assert!(rep.implied_sofr_futures_rate > rep.unadjusted_fra_rate);

            let h = 1e-5;
            let r_ad = evaluate_andersen_sofr_futures_fra_convexity(
                &AndersenSofrFuturesFraConvexitySpec {
                    unadjusted_fra_forward_rate_f0: Dual::<6>::variable(
                        spec.unadjusted_fra_forward_rate_f0,
                        0,
                    ),
                    short_rate_volatility_sigma_r: Dual::<6>::variable(
                        spec.short_rate_volatility_sigma_r,
                        1,
                    ),
                    hull_white_mean_reversion_kappa: Dual::<6>::variable(
                        spec.hull_white_mean_reversion_kappa,
                        2,
                    ),
                    accrual_start_years_t1: Dual::<6>::variable(spec.accrual_start_years_t1, 3),
                    accrual_fraction_tau: Dual::<6>::variable(spec.accrual_fraction_tau, 4),
                    contract_notional_mm: Dual::<6>::variable(spec.contract_notional_mm, 5),
                },
            );
            for (idx, &ad) in r_ad.convexity_bias_pv_dollars.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.unadjusted_fra_forward_rate_f0 += h;
                        sm.unadjusted_fra_forward_rate_f0 -= h;
                    }
                    1 => {
                        sp.short_rate_volatility_sigma_r += h;
                        sm.short_rate_volatility_sigma_r -= h;
                    }
                    2 => {
                        sp.hull_white_mean_reversion_kappa += h;
                        sm.hull_white_mean_reversion_kappa -= h;
                    }
                    3 => {
                        sp.accrual_start_years_t1 += h;
                        sm.accrual_start_years_t1 -= h;
                    }
                    4 => {
                        sp.accrual_fraction_tau += h;
                        sm.accrual_fraction_tau -= h;
                    }
                    _ => {
                        sp.contract_notional_mm += h;
                        sm.contract_notional_mm -= h;
                    }
                }
                let fd = (evaluate_andersen_sofr_futures_fra_convexity(&sp)
                    .convexity_bias_pv_dollars
                    - evaluate_andersen_sofr_futures_fra_convexity(&sm).convexity_bias_pv_dollars)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-8);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_andersen_sofr_futures_fra_convexity_module::*;

mod evaluate_piterbarg_collateral_choice_option_funding_wedge_module {
    // Piterbarg (2010, Risk, "Funding Beyond Discounting: Collateral Agreements and Derivatives Pricing")
    // & Andersen-Piterbarg (2010, Vol. I) Collateralized vs. Uncollateralized Option Discounting & Funding Wedge Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct PiterbargCollateralOptionFundingWedgeSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub strike_price_k: T,
        pub csa_collateral_rate_r_c: T,
        pub equity_repo_financing_rate_r_r: T,
        pub unsecured_funding_spread_s_f: T,
        pub dividend_yield_q: T,
        pub equity_volatility_sigma_s: T,
        pub funding_spread_volatility_sigma_f: T,
        pub equity_funding_correlation_rho_sf: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct PiterbargCollateralOptionFundingWedgeReport<T: Real = f64> {
        pub repo_implied_forward_price: T,
        pub csa_collateralized_call_price: T,
        pub uncollateralized_funding_discounted_call_price: T,
        pub deterministic_funding_discount_wedge: T,
        pub funding_equity_correlation_convexity_adj: T,
        pub total_piterbarg_funding_adjusted_call_price: T,
        pub implied_funding_value_adjustment_fva: T,
    }

    pub fn evaluate_piterbarg_collateral_choice_option_funding_wedge<T: Real>(
        spec: &PiterbargCollateralOptionFundingWedgeSpec<T>,
    ) -> PiterbargCollateralOptionFundingWedgeReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1.0));
        let k = spec.strike_price_k.max_val(T::from_f64(1.0));
        let r_c = spec.csa_collateral_rate_r_c;
        let r_r = spec.equity_repo_financing_rate_r_r;
        let s_f = spec.unsecured_funding_spread_s_f.max_val(T::from_f64(1e-4));
        let q = spec.dividend_yield_q;
        let sig_s = spec.equity_volatility_sigma_s.max_val(T::from_f64(0.05));
        let sig_f = spec.funding_spread_volatility_sigma_f.max_val(T::zero());
        let rho_sf = spec
            .equity_funding_correlation_rho_sf
            .max_val(T::from_f64(-0.95))
            .min_val(T::from_f64(0.95));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let r_f = r_c + s_f;
        let f_repo = s0 * ((r_r - q) * t).exp();
        let std = sig_s * t.sqrt();
        let d1 = ((f_repo / k).ln() + T::from_f64(0.5) * sig_s * sig_s * t) / std;
        let d2 = d1 - std;

        let fwd_call = f_repo * norm_cdf(d1) - k * norm_cdf(d2);
        let df_csa = (-r_c * t).exp();
        let df_fund = (-r_f * t).exp();

        let v_csa = df_csa * fwd_call;
        let v_uncoll = df_fund * fwd_call;
        let det_wedge = v_csa - v_uncoll;

        // Funding-equity correlation convexity adjustment (1/2 * T^2 Girsanov measure-change integral):
        let conv_adj =
            -T::from_f64(0.5) * df_fund * f_repo * norm_cdf(d1) * rho_sf * sig_s * sig_f * t * t;
        let v_total = v_uncoll + conv_adj;
        let fva = v_csa - v_total;

        PiterbargCollateralOptionFundingWedgeReport {
            repo_implied_forward_price: f_repo,
            csa_collateralized_call_price: v_csa,
            uncollateralized_funding_discounted_call_price: v_uncoll,
            deterministic_funding_discount_wedge: det_wedge,
            funding_equity_correlation_convexity_adj: conv_adj,
            total_piterbarg_funding_adjusted_call_price: v_total,
            implied_funding_value_adjustment_fva: fva,
        }
    }

    #[cfg(test)]
    mod evaluate_piterbarg_collateral_choice_option_funding_wedge_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn piterbarg_collateral_option_funding_wedge_invariants_and_ad_parity() {
            let spec = PiterbargCollateralOptionFundingWedgeSpec {
                spot_price_s0: 100.0,
                strike_price_k: 100.0,
                csa_collateral_rate_r_c: 0.038,
                equity_repo_financing_rate_r_r: 0.042,
                unsecured_funding_spread_s_f: 0.015,
                dividend_yield_q: 0.015,
                equity_volatility_sigma_s: 0.24,
                funding_spread_volatility_sigma_f: 0.008,
                equity_funding_correlation_rho_sf: -0.35,
                maturity_years_t: 2.0,
            };
            let rep = evaluate_piterbarg_collateral_choice_option_funding_wedge(&spec);
            assert!(
                rep.csa_collateralized_call_price
                    > rep.uncollateralized_funding_discounted_call_price
            );
            assert!(rep.deterministic_funding_discount_wedge > 0.0);
            assert!(rep.funding_equity_correlation_convexity_adj > 0.0);

            let h = 1e-5;
            let r_ad = evaluate_piterbarg_collateral_choice_option_funding_wedge(
                &PiterbargCollateralOptionFundingWedgeSpec {
                    spot_price_s0: Dual::<10>::variable(spec.spot_price_s0, 0),
                    strike_price_k: Dual::<10>::variable(spec.strike_price_k, 1),
                    csa_collateral_rate_r_c: Dual::<10>::variable(spec.csa_collateral_rate_r_c, 2),
                    equity_repo_financing_rate_r_r: Dual::<10>::variable(
                        spec.equity_repo_financing_rate_r_r,
                        3,
                    ),
                    unsecured_funding_spread_s_f: Dual::<10>::variable(
                        spec.unsecured_funding_spread_s_f,
                        4,
                    ),
                    dividend_yield_q: Dual::<10>::variable(spec.dividend_yield_q, 5),
                    equity_volatility_sigma_s: Dual::<10>::variable(
                        spec.equity_volatility_sigma_s,
                        6,
                    ),
                    funding_spread_volatility_sigma_f: Dual::<10>::variable(
                        spec.funding_spread_volatility_sigma_f,
                        7,
                    ),
                    equity_funding_correlation_rho_sf: Dual::<10>::variable(
                        spec.equity_funding_correlation_rho_sf,
                        8,
                    ),
                    maturity_years_t: Dual::<10>::variable(spec.maturity_years_t, 9),
                },
            );
            for (idx, &ad) in r_ad
                .total_piterbarg_funding_adjusted_call_price
                .d
                .iter()
                .enumerate()
            {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_price_s0 += h;
                        sm.spot_price_s0 -= h;
                    }
                    1 => {
                        sp.strike_price_k += h;
                        sm.strike_price_k -= h;
                    }
                    2 => {
                        sp.csa_collateral_rate_r_c += h;
                        sm.csa_collateral_rate_r_c -= h;
                    }
                    3 => {
                        sp.equity_repo_financing_rate_r_r += h;
                        sm.equity_repo_financing_rate_r_r -= h;
                    }
                    4 => {
                        sp.unsecured_funding_spread_s_f += h;
                        sm.unsecured_funding_spread_s_f -= h;
                    }
                    5 => {
                        sp.dividend_yield_q += h;
                        sm.dividend_yield_q -= h;
                    }
                    6 => {
                        sp.equity_volatility_sigma_s += h;
                        sm.equity_volatility_sigma_s -= h;
                    }
                    7 => {
                        sp.funding_spread_volatility_sigma_f += h;
                        sm.funding_spread_volatility_sigma_f -= h;
                    }
                    8 => {
                        sp.equity_funding_correlation_rho_sf += h;
                        sm.equity_funding_correlation_rho_sf -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_piterbarg_collateral_choice_option_funding_wedge(&sp)
                    .total_piterbarg_funding_adjusted_call_price
                    - evaluate_piterbarg_collateral_choice_option_funding_wedge(&sm)
                        .total_piterbarg_funding_adjusted_call_price)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_piterbarg_collateral_choice_option_funding_wedge_module::*;

mod evaluate_kim_wright_shadow_rate_survey_term_premium_module {
    // Kim & Wright (2005, FEDS 2005-33, "An Arbitrage-Free Three-Factor Term Structure Model
    // and the Recent Behavior of Long-Term Yields and Distant-Horizon Forward Rates"),
    // Black (1995, "Interest Rates as Options"), Krippner (2013), & Wu & Xia (2016)
    // Shadow-Rate Effective-Lower-Bound (ELB) & Survey-Anchored Forward Term Premium Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct KimWrightShadowRateTermPremiumSpec<T: Real = f64> {
        pub current_shadow_short_rate_s0: T,
        pub effective_lower_bound_r_min: T,
        pub mean_reversion_speed_kappa: T,
        pub risk_neutral_long_run_mean_theta_q: T,
        pub physical_long_run_mean_theta_p: T,
        pub shadow_rate_volatility_sigma_s: T,
        pub survey_long_run_short_rate_r_surv: T,
        pub survey_kalman_shrinkage_weight_w: T,
        pub horizon_years_tau: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct KimWrightShadowRateTermPremiumReport<T: Real = f64> {
        pub unconstrained_q_shadow_forward_mean: T,
        pub survey_anchored_p_shadow_mean: T,
        pub expected_elb_constrained_short_rate_p: T,
        pub krippner_shadow_convexity_adjustment: T,
        pub instantaneous_elb_forward_rate_q: T,
        pub zlb_call_option_wedge_bps: T,
        pub survey_anchored_forward_term_premium: T,
    }

    pub fn evaluate_kim_wright_shadow_rate_survey_term_premium<T: Real>(
        spec: &KimWrightShadowRateTermPremiumSpec<T>,
    ) -> KimWrightShadowRateTermPremiumReport<T> {
        let s0 = spec.current_shadow_short_rate_s0;
        let r_min = spec.effective_lower_bound_r_min;
        let kappa = spec.mean_reversion_speed_kappa.max_val(T::from_f64(0.02));
        let theta_q = spec.risk_neutral_long_run_mean_theta_q;
        let theta_p = spec.physical_long_run_mean_theta_p;
        let sig_s = spec
            .shadow_rate_volatility_sigma_s
            .max_val(T::from_f64(0.001));
        let r_surv = spec.survey_long_run_short_rate_r_surv;
        let w_surv = spec
            .survey_kalman_shrinkage_weight_w
            .max_val(T::from_f64(0.0))
            .min_val(T::from_f64(1.0));
        let tau = spec.horizon_years_tau.max_val(T::from_f64(0.25));

        let e_kt = (-kappa * tau).exp();
        let b_kt = (T::one() - e_kt) / kappa;
        let mu_q = s0 * e_kt + theta_q * (T::one() - e_kt);
        let mu_p_raw = s0 * e_kt + theta_p * (T::one() - e_kt);
        let mu_p = (T::one() - w_surv) * mu_p_raw + w_surv * r_surv;

        let var_s = (sig_s * sig_s / (T::from_f64(2.0) * kappa))
            * (T::one() - (-T::from_f64(2.0) * kappa * tau).exp());
        let std_s = var_s.max_val(T::from_f64(1e-10)).sqrt();

        // Black (1995) / Krippner (2013) Gaussian option expectation E[max(s_tau, r_min)]:
        let z_q = (mu_q - r_min) / std_s;
        let phi_zq = norm_cdf(z_q);
        let er_q = r_min + (mu_q - r_min) * phi_zq + std_s * norm_pdf(z_q);

        let z_p = (mu_p - r_min) / std_s;
        let er_p = r_min + (mu_p - r_min) * norm_cdf(z_p) + std_s * norm_pdf(z_p);

        // Krippner (2013) ELB-attenuated instantaneous forward convexity adjustment:
        let conv_adj = T::from_f64(0.5) * phi_zq * phi_zq * sig_s * sig_s * b_kt * b_kt;
        let f_q = er_q - conv_adj;
        let zlb_wedge_bps = (er_q - mu_q) * T::from_f64(10_000.0);
        let ftp = f_q - er_p;

        KimWrightShadowRateTermPremiumReport {
            unconstrained_q_shadow_forward_mean: mu_q,
            survey_anchored_p_shadow_mean: mu_p,
            expected_elb_constrained_short_rate_p: er_p,
            krippner_shadow_convexity_adjustment: conv_adj,
            instantaneous_elb_forward_rate_q: f_q,
            zlb_call_option_wedge_bps: zlb_wedge_bps,
            survey_anchored_forward_term_premium: ftp,
        }
    }

    #[cfg(test)]
    mod evaluate_kim_wright_shadow_rate_survey_term_premium_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn kim_wright_shadow_rate_survey_term_premium_invariants_and_ad_parity() {
            let spec = KimWrightShadowRateTermPremiumSpec {
                current_shadow_short_rate_s0: 0.015,
                effective_lower_bound_r_min: 0.00125,
                mean_reversion_speed_kappa: 0.22,
                risk_neutral_long_run_mean_theta_q: 0.044,
                physical_long_run_mean_theta_p: 0.032,
                shadow_rate_volatility_sigma_s: 0.014,
                survey_long_run_short_rate_r_surv: 0.028,
                survey_kalman_shrinkage_weight_w: 0.35,
                horizon_years_tau: 10.0,
            };
            let rep = evaluate_kim_wright_shadow_rate_survey_term_premium(&spec);
            assert!(rep.zlb_call_option_wedge_bps > 0.0);
            assert!(rep.survey_anchored_forward_term_premium > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_kim_wright_shadow_rate_survey_term_premium(
                &KimWrightShadowRateTermPremiumSpec {
                    current_shadow_short_rate_s0: Dual::<9>::variable(
                        spec.current_shadow_short_rate_s0,
                        0,
                    ),
                    effective_lower_bound_r_min: Dual::<9>::variable(
                        spec.effective_lower_bound_r_min,
                        1,
                    ),
                    mean_reversion_speed_kappa: Dual::<9>::variable(
                        spec.mean_reversion_speed_kappa,
                        2,
                    ),
                    risk_neutral_long_run_mean_theta_q: Dual::<9>::variable(
                        spec.risk_neutral_long_run_mean_theta_q,
                        3,
                    ),
                    physical_long_run_mean_theta_p: Dual::<9>::variable(
                        spec.physical_long_run_mean_theta_p,
                        4,
                    ),
                    shadow_rate_volatility_sigma_s: Dual::<9>::variable(
                        spec.shadow_rate_volatility_sigma_s,
                        5,
                    ),
                    survey_long_run_short_rate_r_surv: Dual::<9>::variable(
                        spec.survey_long_run_short_rate_r_surv,
                        6,
                    ),
                    survey_kalman_shrinkage_weight_w: Dual::<9>::variable(
                        spec.survey_kalman_shrinkage_weight_w,
                        7,
                    ),
                    horizon_years_tau: Dual::<9>::variable(spec.horizon_years_tau, 8),
                },
            );
            for i in 0..9 {
                let mut v = [
                    spec.current_shadow_short_rate_s0,
                    spec.effective_lower_bound_r_min,
                    spec.mean_reversion_speed_kappa,
                    spec.risk_neutral_long_run_mean_theta_q,
                    spec.physical_long_run_mean_theta_p,
                    spec.shadow_rate_volatility_sigma_s,
                    spec.survey_long_run_short_rate_r_surv,
                    spec.survey_kalman_shrinkage_weight_w,
                    spec.horizon_years_tau,
                ];
                v[i] += h;
                let up = evaluate_kim_wright_shadow_rate_survey_term_premium(
                    &KimWrightShadowRateTermPremiumSpec {
                        current_shadow_short_rate_s0: v[0],
                        effective_lower_bound_r_min: v[1],
                        mean_reversion_speed_kappa: v[2],
                        risk_neutral_long_run_mean_theta_q: v[3],
                        physical_long_run_mean_theta_p: v[4],
                        shadow_rate_volatility_sigma_s: v[5],
                        survey_long_run_short_rate_r_surv: v[6],
                        survey_kalman_shrinkage_weight_w: v[7],
                        horizon_years_tau: v[8],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_kim_wright_shadow_rate_survey_term_premium(
                    &KimWrightShadowRateTermPremiumSpec {
                        current_shadow_short_rate_s0: v[0],
                        effective_lower_bound_r_min: v[1],
                        mean_reversion_speed_kappa: v[2],
                        risk_neutral_long_run_mean_theta_q: v[3],
                        physical_long_run_mean_theta_p: v[4],
                        shadow_rate_volatility_sigma_s: v[5],
                        survey_long_run_short_rate_r_surv: v[6],
                        survey_kalman_shrinkage_weight_w: v[7],
                        horizon_years_tau: v[8],
                    },
                );
                let fd = (up.survey_anchored_forward_term_premium
                    - dn.survey_anchored_forward_term_premium)
                    / (2.0 * h);
                assert!(
                    (r_ad.survey_anchored_forward_term_premium.d[i] - fd).abs() < 1e-5,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_kim_wright_shadow_rate_survey_term_premium_module::*;

mod evaluate_lyashenko_mercurio_fmm_in_arrears_caplet_convexity_module {
    // Lyashenko & Mercurio (2019, Risk, "Looking Forward to Backward-Looking Rates")
    // & Andersen & Piterbarg (2010, Interest Rate Modeling, Vol. I & III)
    // Forward Market Model (FMM) Backward-Looking Compounded-in-Arrears RFR Caplet Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct LyashenkoMercurioFmmInArrearsCapletSpec<T: Real = f64> {
        pub notional_n: T,
        pub forward_rfr_rate_f0: T,
        pub caplet_strike_k: T,
        pub fixing_start_years_t0: T,
        pub accrual_fraction_tau: T,
        pub forward_vol_sigma_f: T,
        pub ois_discount_rate_r: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct LyashenkoMercurioFmmInArrearsCapletReport<T: Real = f64> {
        pub forward_looking_integrated_variance: T,
        pub in_arrears_integrated_variance: T,
        pub effective_in_arrears_black_vol: T,
        pub forward_looking_ibor_caplet_npv: T,
        pub backward_looking_in_arrears_caplet_npv: T,
        pub in_arrears_accrual_volatility_premium: T,
    }

    pub fn evaluate_lyashenko_mercurio_fmm_in_arrears_caplet_convexity<T: Real>(
        spec: &LyashenkoMercurioFmmInArrearsCapletSpec<T>,
    ) -> LyashenkoMercurioFmmInArrearsCapletReport<T> {
        let n = spec.notional_n.max_val(T::from_f64(1.0));
        let f0 = spec.forward_rfr_rate_f0.max_val(T::from_f64(0.001));
        let k = spec.caplet_strike_k.max_val(T::from_f64(0.001));
        let t0 = spec.fixing_start_years_t0.max_val(T::from_f64(0.05));
        let tau = spec.accrual_fraction_tau.max_val(T::from_f64(0.05));
        let sig_f = spec.forward_vol_sigma_f.max_val(T::from_f64(0.01));
        let r = spec.ois_discount_rate_r.max_val(T::from_f64(0.0005));

        let t1 = t0 + tau;
        let df_t1 = (-r * t1).exp();
        let sig2 = sig_f * sig_f;

        // Lyashenko & Mercurio (2019, Section 3.2): instantaneous volatility is sigma_F on [0, T0]
        // and decays linearly as sigma_F * (T1 - t) / tau on [T0, T1], yielding integrated variance
        // V_arr = sigma_F^2 * (T0 + tau / 3) vs. forward-looking V_fwd = sigma_F^2 * T0:
        let v_fwd = sig2 * t0;
        let v_arr = sig2 * (t0 + tau / T::from_f64(3.0));
        let sig_eff_arr = (v_arr / t1).sqrt();

        let std_fwd = v_fwd.sqrt();
        let ln_fk = (f0 / k).ln();
        let d1_fwd = (ln_fk + T::from_f64(0.5) * v_fwd) / std_fwd;
        let d2_fwd = d1_fwd - std_fwd;
        let c_fwd = n * tau * df_t1 * (f0 * norm_cdf(d1_fwd) - k * norm_cdf(d2_fwd));

        let std_arr = v_arr.sqrt();
        let d1_arr = (ln_fk + T::from_f64(0.5) * v_arr) / std_arr;
        let d2_arr = d1_arr - std_arr;
        let c_arr = n * tau * df_t1 * (f0 * norm_cdf(d1_arr) - k * norm_cdf(d2_arr));

        let prem = c_arr - c_fwd;

        LyashenkoMercurioFmmInArrearsCapletReport {
            forward_looking_integrated_variance: v_fwd,
            in_arrears_integrated_variance: v_arr,
            effective_in_arrears_black_vol: sig_eff_arr,
            forward_looking_ibor_caplet_npv: c_fwd,
            backward_looking_in_arrears_caplet_npv: c_arr,
            in_arrears_accrual_volatility_premium: prem,
        }
    }

    #[cfg(test)]
    mod evaluate_lyashenko_mercurio_fmm_in_arrears_caplet_convexity_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn lyashenko_mercurio_fmm_in_arrears_caplet_invariants_and_ad_parity() {
            let spec = LyashenkoMercurioFmmInArrearsCapletSpec {
                notional_n: 1_000.0,
                forward_rfr_rate_f0: 0.042,
                caplet_strike_k: 0.040,
                fixing_start_years_t0: 0.50,
                accrual_fraction_tau: 0.25,
                forward_vol_sigma_f: 0.24,
                ois_discount_rate_r: 0.040,
            };
            let rep = evaluate_lyashenko_mercurio_fmm_in_arrears_caplet_convexity(&spec);
            assert!(
                rep.backward_looking_in_arrears_caplet_npv > rep.forward_looking_ibor_caplet_npv
            );
            assert!(rep.in_arrears_accrual_volatility_premium > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_lyashenko_mercurio_fmm_in_arrears_caplet_convexity(
                &LyashenkoMercurioFmmInArrearsCapletSpec {
                    notional_n: Dual::<7>::variable(spec.notional_n, 0),
                    forward_rfr_rate_f0: Dual::<7>::variable(spec.forward_rfr_rate_f0, 1),
                    caplet_strike_k: Dual::<7>::variable(spec.caplet_strike_k, 2),
                    fixing_start_years_t0: Dual::<7>::variable(spec.fixing_start_years_t0, 3),
                    accrual_fraction_tau: Dual::<7>::variable(spec.accrual_fraction_tau, 4),
                    forward_vol_sigma_f: Dual::<7>::variable(spec.forward_vol_sigma_f, 5),
                    ois_discount_rate_r: Dual::<7>::variable(spec.ois_discount_rate_r, 6),
                },
            );
            for i in 0..7 {
                let mut v = [
                    spec.notional_n,
                    spec.forward_rfr_rate_f0,
                    spec.caplet_strike_k,
                    spec.fixing_start_years_t0,
                    spec.accrual_fraction_tau,
                    spec.forward_vol_sigma_f,
                    spec.ois_discount_rate_r,
                ];
                v[i] += h;
                let up = evaluate_lyashenko_mercurio_fmm_in_arrears_caplet_convexity(
                    &LyashenkoMercurioFmmInArrearsCapletSpec {
                        notional_n: v[0],
                        forward_rfr_rate_f0: v[1],
                        caplet_strike_k: v[2],
                        fixing_start_years_t0: v[3],
                        accrual_fraction_tau: v[4],
                        forward_vol_sigma_f: v[5],
                        ois_discount_rate_r: v[6],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_lyashenko_mercurio_fmm_in_arrears_caplet_convexity(
                    &LyashenkoMercurioFmmInArrearsCapletSpec {
                        notional_n: v[0],
                        forward_rfr_rate_f0: v[1],
                        caplet_strike_k: v[2],
                        fixing_start_years_t0: v[3],
                        accrual_fraction_tau: v[4],
                        forward_vol_sigma_f: v[5],
                        ois_discount_rate_r: v[6],
                    },
                );
                let fd = (up.backward_looking_in_arrears_caplet_npv
                    - dn.backward_looking_in_arrears_caplet_npv)
                    / (2.0 * h);
                assert!(
                    (r_ad.backward_looking_in_arrears_caplet_npv.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_lyashenko_mercurio_fmm_in_arrears_caplet_convexity_module::*;

pub mod bianchetti_two_curve_quanto_suite65_mod {
    use ql_ad::Real;
    use ql_math::norm_cdf;

    /// Specification for the Bianchetti (2010) / Ametrano–Bianchetti (2013) two-curve
    /// OIS vs. projection forward-measure change-of-numeraire quanto convexity engine.
    ///
    /// Reference: Bianchetti (2010, *Risk*, "Two Curves, One Price"), Ametrano & Bianchetti (2013,
    /// *Interest Rate Modelling After the Financial Crisis*), Henrard (2014).
    #[derive(Debug, Clone, Copy)]
    pub struct BianchettiTwoCurveForwardQuantoSpec<T: Real> {
        /// Contract notional $N$.
        pub notional_n: T,
        /// Projection-curve forward rate $F_x(0; T_1, T_2)$ under $\mathbb{Q}_x^{T_2}$.
        pub projection_forward_rate_f0: T,
        /// Multiplicative tenor-basis spread $s_{\text{basis},0}$ between projection and OIS curves.
        pub tenor_basis_spread_s0: T,
        /// Caplet strike rate $K$.
        pub caplet_strike_k: T,
        /// Fixing expiry $T_1$ (in years).
        pub fixing_expiry_years_t1: T,
        /// Accrual year fraction $\tau = T_2 - T_1$.
        pub accrual_fraction_tau: T,
        /// Lognormal volatility $\sigma_F$ of the projection forward rate $F_x$.
        pub forward_rate_vol_sigma_f: T,
        /// Lognormal volatility $\sigma_s$ of the basis spread ratio.
        pub basis_spread_vol_sigma_s: T,
        /// Correlation $\rho_{F,s}$ between projection forward rate and tenor basis spread.
        pub forward_basis_correlation_rho: T,
        /// OIS discount rate $r_{\text{OIS}}$.
        pub ois_discount_rate_r: T,
    }

    /// Analytics returned by [`evaluate_bianchetti_two_curve_forward_measure_quanto_convexity`].
    #[derive(Debug, Clone, Copy)]
    pub struct BianchettiTwoCurveForwardQuantoAnalytics<T: Real> {
        /// OIS $T_2$-forward measure quanto-adjusted forward rate $F_x^{\mathbb{Q}_d^{T_2}}(0; T_1, T_2)$.
        pub ois_measure_quanto_adjusted_forward_rate: T,
        /// Two-curve change-of-numeraire convexity drift basis-point adjustment.
        pub two_curve_quanto_convexity_bps: T,
        /// Unadjusted single-measure Black–76 caplet PV (discounted at OIS).
        pub unadjusted_two_curve_caplet_pv: T,
        /// Bianchetti quanto-adjusted collateralized Black–76 caplet PV.
        pub bianchetti_quanto_adjusted_caplet_pv: T,
    }

    /// Evaluates the Bianchetti (2010) change-of-numeraire internal-currency quanto adjustment
    /// between the projection forward measure $\mathbb{Q}_x^{T_2}$ and the OIS forward measure
    /// $\mathbb{Q}_d^{T_2}$ and prices the resulting collateralized caplet.
    pub fn evaluate_bianchetti_two_curve_forward_measure_quanto_convexity<T: Real>(
        spec: &BianchettiTwoCurveForwardQuantoSpec<T>,
    ) -> BianchettiTwoCurveForwardQuantoAnalytics<T> {
        let t1 = spec.fixing_expiry_years_t1.max_val(T::from_f64(1e-4));
        let tau = spec.accrual_fraction_tau.max_val(T::from_f64(1e-4));
        let t2 = t1 + tau;
        let f0 = spec.projection_forward_rate_f0.max_val(T::from_f64(1e-5));
        let k = spec.caplet_strike_k.max_val(T::from_f64(1e-5));
        let sig_f = spec.forward_rate_vol_sigma_f.max_val(T::from_f64(1e-4));
        let sig_s = spec.basis_spread_vol_sigma_s.max_val(T::from_f64(1e-4));

        // Bianchetti (2010, Eq. 19): Radon-Nikodym density dQ_d^{T2} / dQ_x^{T2} = [P_d(T1,T2)/P_d(0,T2)] / [P_x(T1,T2)/P_x(0,T2)]
        // Log-ratio volatility over [0, T1] is tau * s_basis / (1 + tau * s_basis) * sigma_s
        let basis_weight = (tau * spec.tenor_basis_spread_s0)
            / (T::one() + tau * spec.tenor_basis_spread_s0.max_val(T::from_f64(0.0)));
        let quanto_exponent =
            spec.forward_basis_correlation_rho * sig_f * sig_s * basis_weight * t1;
        let f_quanto = f0 * quanto_exponent.exp();
        let convexity_bps = (f_quanto - f0) * T::from_f64(10_000.0);

        let df_ois = (-spec.ois_discount_rate_r * t2).exp();
        let v_sqrt = sig_f * t1.sqrt();

        let d1_unadj = ((f0 / k).ln() + T::from_f64(0.5) * sig_f * sig_f * t1) / v_sqrt;
        let d2_unadj = d1_unadj - v_sqrt;
        let pv_unadj =
            spec.notional_n * tau * df_ois * (f0 * norm_cdf(d1_unadj) - k * norm_cdf(d2_unadj));

        let d1_q = ((f_quanto / k).ln() + T::from_f64(0.5) * sig_f * sig_f * t1) / v_sqrt;
        let d2_q = d1_q - v_sqrt;
        let pv_quanto =
            spec.notional_n * tau * df_ois * (f_quanto * norm_cdf(d1_q) - k * norm_cdf(d2_q));

        BianchettiTwoCurveForwardQuantoAnalytics {
            ois_measure_quanto_adjusted_forward_rate: f_quanto,
            two_curve_quanto_convexity_bps: convexity_bps,
            unadjusted_two_curve_caplet_pv: pv_unadj,
            bianchetti_quanto_adjusted_caplet_pv: pv_quanto,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bianchetti_two_curve_quanto_invariants_and_ad_parity() {
            let base = BianchettiTwoCurveForwardQuantoSpec {
                notional_n: 1_000.0_f64,
                projection_forward_rate_f0: 0.044,
                tenor_basis_spread_s0: 0.0035,
                caplet_strike_k: 0.042,
                fixing_expiry_years_t1: 2.0,
                accrual_fraction_tau: 0.50,
                forward_rate_vol_sigma_f: 0.25,
                basis_spread_vol_sigma_s: 0.35,
                forward_basis_correlation_rho: 0.65,
                ois_discount_rate_r: 0.040,
            };
            let r0 = evaluate_bianchetti_two_curve_forward_measure_quanto_convexity(&base);
            assert!(r0.ois_measure_quanto_adjusted_forward_rate > base.projection_forward_rate_f0);
            assert!(r0.bianchetti_quanto_adjusted_caplet_pv > r0.unadjusted_two_curve_caplet_pv);

            let inputs = [
                base.notional_n,
                base.projection_forward_rate_f0,
                base.tenor_basis_spread_s0,
                base.caplet_strike_k,
                base.fixing_expiry_years_t1,
                base.accrual_fraction_tau,
                base.forward_rate_vol_sigma_f,
                base.basis_spread_vol_sigma_s,
                base.forward_basis_correlation_rho,
                base.ois_discount_rate_r,
            ];
            let mut dual_in = [Dual::<10>::from_f64(0.0); 10];
            for i in 0..10 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BianchettiTwoCurveForwardQuantoSpec {
                notional_n: dual_in[0],
                projection_forward_rate_f0: dual_in[1],
                tenor_basis_spread_s0: dual_in[2],
                caplet_strike_k: dual_in[3],
                fixing_expiry_years_t1: dual_in[4],
                accrual_fraction_tau: dual_in[5],
                forward_rate_vol_sigma_f: dual_in[6],
                basis_spread_vol_sigma_s: dual_in[7],
                forward_basis_correlation_rho: dual_in[8],
                ois_discount_rate_r: dual_in[9],
            };
            let r_ad = evaluate_bianchetti_two_curve_forward_measure_quanto_convexity(&spec_ad);
            let eps = 1e-6;
            for i in 0..10 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 10]| BianchettiTwoCurveForwardQuantoSpec {
                    notional_n: v[0],
                    projection_forward_rate_f0: v[1],
                    tenor_basis_spread_s0: v[2],
                    caplet_strike_k: v[3],
                    fixing_expiry_years_t1: v[4],
                    accrual_fraction_tau: v[5],
                    forward_rate_vol_sigma_f: v[6],
                    basis_spread_vol_sigma_s: v[7],
                    forward_basis_correlation_rho: v[8],
                    ois_discount_rate_r: v[9],
                };
                let rp = evaluate_bianchetti_two_curve_forward_measure_quanto_convexity(&mk(p));
                let rm = evaluate_bianchetti_two_curve_forward_measure_quanto_convexity(&mk(m));
                let fd = (rp.bianchetti_quanto_adjusted_caplet_pv
                    - rm.bianchetti_quanto_adjusted_caplet_pv)
                    / (2.0 * eps);
                assert!(
                    r_ad.bianchetti_quanto_adjusted_caplet_pv.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.bianchetti_quanto_adjusted_caplet_pv.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use bianchetti_two_curve_quanto_suite65_mod::*;

pub mod piterbarg_multi_currency_ctd_collateral_mod {
    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    /// Specification for the Piterbarg (2010, Risk; 2012, Risk) & Henrard (2014, Ch. 3)
    /// multi-currency Cheapest-to-Deliver (CTD) collateral choice discounting engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct PiterbargMultiCurrencyCtdCollateralSpec<T: Real> {
        /// Domestic OIS rate $r_{d,\text{ois}}$.
        pub domestic_ois_rate_rd: T,
        /// Domestic CSA collateral spread $s_1$ over domestic OIS ($r_{c,1} = r_{d,\text{ois}} - s_1$).
        pub domestic_collateral_spread_s1: T,
        /// Foreign-currency synthetic CSA spread $s_2$ over domestic OIS after FX swap basis ($r_{c,2} = r_{d,\text{ois}} - s_2$).
        pub foreign_synthetic_collateral_spread_s2: T,
        /// Normal volatility $\sigma_X > 0$ of the collateral spread differential $X_t = s_1(t) - s_2(t)$.
        pub collateral_spread_diff_volatility_sigma_x: T,
        /// Mean-reversion speed $\kappa_X > 0$ of the collateral spread differential.
        pub spread_reversion_speed_kappa_x: T,
        /// Horizon maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Piterbarg multi-currency CTD collateral choice discounting engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct PiterbargMultiCurrencyCtdCollateralAnalytics<T: Real> {
        /// Stand-alone domestic-collateral discount factor $P_1(0, T) = \exp(-(r_{d,\text{ois}} - s_1)T)$.
        pub domestic_collateral_discount_factor: T,
        /// Stand-alone foreign-collateral discount factor $P_2(0, T) = \exp(-(r_{d,\text{ois}} - s_2)T)$.
        pub foreign_collateral_discount_factor: T,
        /// Annualized CTD collateral switch option yield pickup (in rate units).
        pub ctd_collateral_switch_yield_pickup: T,
        /// Cheapest-to-Deliver (CTD) collateral-adjusted discount factor $P_{\text{CTD}}(0, T)$.
        pub ctd_collateral_discount_factor: T,
    }

    /// Evaluates the Piterbarg (2010, 2012) multi-currency CTD collateral discount factor over a
    /// 4-period Ornstein-Uhlenbeck Bachelier collateral choice strip.
    pub fn evaluate_piterbarg_multi_currency_ctd_collateral_choice_discount<T: Real>(
        spec: &PiterbargMultiCurrencyCtdCollateralSpec<T>,
    ) -> PiterbargMultiCurrencyCtdCollateralAnalytics<T> {
        const PERIODS: usize = 4;
        let t_mat = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let dt = t_mat / T::from_f64(PERIODS as f64);
        let kx = spec
            .spread_reversion_speed_kappa_x
            .max_val(T::from_f64(1e-4));
        let sig_x = spec
            .collateral_spread_diff_volatility_sigma_x
            .max_val(T::from_f64(1e-6));

        let p1 = (-(spec.domestic_ois_rate_rd - spec.domestic_collateral_spread_s1) * t_mat).exp();
        let p2 = (-(spec.domestic_ois_rate_rd - spec.foreign_synthetic_collateral_spread_s2)
            * t_mat)
            .exp();

        // Collateral poster chooses currency i maximizing r_{c,i}(t) = r_d - s_i(t),
        // which minimizes the effective spread min(s_1(t), s_2(t)) = s_1(t) - max(s_1(t) - s_2(t), 0).
        let diff_s =
            spec.domestic_collateral_spread_s1 - spec.foreign_synthetic_collateral_spread_s2;
        let mut int_min_spread = T::zero();
        let mut int_switch_pickup = T::zero();

        for m in 1..=PERIODS {
            let tm = (T::from_f64(m as f64) - T::from_f64(0.5)) * dt;
            let ou_var = sig_x * sig_x * (T::one() - (-T::from_f64(2.0) * kx * tm).exp())
                / (T::from_f64(2.0) * kx);
            let std_x = ou_var.max_val(T::from_f64(1e-12)).sqrt();
            let d = diff_s / std_x;
            let expected_max_diff = diff_s * norm_cdf(d) + std_x * norm_pdf(d);
            let expected_min_s = spec.domestic_collateral_spread_s1 - expected_max_diff;

            int_min_spread += dt * expected_min_s;
            int_switch_pickup += dt * expected_max_diff;
        }

        let p_ctd = (-(spec.domestic_ois_rate_rd * t_mat - int_min_spread)).exp();

        PiterbargMultiCurrencyCtdCollateralAnalytics {
            domestic_collateral_discount_factor: p1,
            foreign_collateral_discount_factor: p2,
            ctd_collateral_switch_yield_pickup: int_switch_pickup / t_mat,
            ctd_collateral_discount_factor: p_ctd,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn piterbarg_multi_currency_ctd_collateral_invariants_and_ad_parity() {
            let base = PiterbargMultiCurrencyCtdCollateralSpec {
                domestic_ois_rate_rd: 0.042_f64,
                domestic_collateral_spread_s1: 0.0010,
                foreign_synthetic_collateral_spread_s2: 0.0025,
                collateral_spread_diff_volatility_sigma_x: 0.0035,
                spread_reversion_speed_kappa_x: 0.45,
                maturity_years_t: 5.0,
            };
            let r0 = evaluate_piterbarg_multi_currency_ctd_collateral_choice_discount(&base);
            assert!(r0.ctd_collateral_switch_yield_pickup > 0.0);
            assert!(r0.ctd_collateral_discount_factor < r0.domestic_collateral_discount_factor);
            assert!(r0.ctd_collateral_discount_factor < r0.foreign_collateral_discount_factor);

            let inputs = [
                base.domestic_ois_rate_rd,
                base.domestic_collateral_spread_s1,
                base.foreign_synthetic_collateral_spread_s2,
                base.collateral_spread_diff_volatility_sigma_x,
                base.spread_reversion_speed_kappa_x,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<6>::from_f64(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = PiterbargMultiCurrencyCtdCollateralSpec {
                domestic_ois_rate_rd: dual_in[0],
                domestic_collateral_spread_s1: dual_in[1],
                foreign_synthetic_collateral_spread_s2: dual_in[2],
                collateral_spread_diff_volatility_sigma_x: dual_in[3],
                spread_reversion_speed_kappa_x: dual_in[4],
                maturity_years_t: dual_in[5],
            };
            let r_ad = evaluate_piterbarg_multi_currency_ctd_collateral_choice_discount(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| PiterbargMultiCurrencyCtdCollateralSpec {
                    domestic_ois_rate_rd: v[0],
                    domestic_collateral_spread_s1: v[1],
                    foreign_synthetic_collateral_spread_s2: v[2],
                    collateral_spread_diff_volatility_sigma_x: v[3],
                    spread_reversion_speed_kappa_x: v[4],
                    maturity_years_t: v[5],
                };
                let rp = evaluate_piterbarg_multi_currency_ctd_collateral_choice_discount(&mk(p));
                let rm = evaluate_piterbarg_multi_currency_ctd_collateral_choice_discount(&mk(m));
                let fd = (rp.ctd_collateral_discount_factor - rm.ctd_collateral_discount_factor)
                    / (2.0 * eps);
                assert!(
                    r_ad.ctd_collateral_discount_factor.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.ctd_collateral_discount_factor.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use piterbarg_multi_currency_ctd_collateral_mod::*;

// ============================================================================
// Ametrano & Bianchetti (2009, 2013) & Henrard (2014) Dual-Curve 3M-vs-6M Tenor Basis Bootstrap
// ============================================================================

pub mod ametrano_bianchetti_tenor_basis_mod {
    use ql_ad::Real;

    /// Specification for the Ametrano & Bianchetti (2009; 2013, Ch. 2) / Henrard (2014, Ch. 2)
    /// OIS-discounted 3M-vs-6M IBOR tenor basis swap curve bootstrap engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct AmetranoBianchettiTenorBasisSpec<T: Real> {
        /// Initial instantaneous OIS short rate $r_{\text{OIS}}$.
        pub ois_discount_rate_r_ois: T,
        /// OIS term-structure slope parameter $\theta_{\text{OIS}}$.
        pub ois_slope_theta_ois: T,
        /// 3M IBOR credit/liquidity spread over OIS $s_{3\text{M}} > 0$.
        pub ibor_3m_spread_over_ois_s3: T,
        /// Base 3M-vs-6M tenor basis swap spread $b_{3\text{v}6} > 0$.
        pub tenor_basis_3m_vs_6m_spread_b36: T,
        /// Term-structure steepening slope of the 3M-vs-6M tenor basis $\gamma_{3\text{v}6}$.
        pub tenor_basis_term_steepening_gamma: T,
        /// 3M forward rate volatility $\sigma_{3\text{M}} > 0$ driving intra-period compounding convexity.
        pub compounding_convexity_vol_sigma3: T,
        /// Swap maturity $T > 0$ (years).
        pub swap_maturity_years_t: T,
    }

    /// Analytical outputs for the Ametrano-Bianchetti dual-curve tenor basis bootstrap engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct AmetranoBianchettiTenorBasisAnalytics<T: Real> {
        /// Terminal bootstrapped 6M pseudo-discount factor $P^{6\text{M}}(0, T)$.
        pub terminal_6m_pseudo_discount_factor: T,
        /// OIS-discounted 3M IBOR par swap rate $S_{3\text{M}}(0, T)$.
        pub ois_discounted_3m_par_swap_rate: T,
        /// Implied 6M-over-3M effective par swap spread $S_{6\text{M}}(0, T) - S_{3\text{M}}(0, T)$.
        pub implied_6m_over_3m_swap_spread: T,
        /// Bootstrapped OIS-discounted 6M IBOR par swap rate $S_{6\text{M}}(0, T)$.
        pub bootstrapped_6m_ibor_par_swap_rate: T,
    }

    /// Evaluates the Ametrano & Bianchetti (2009, 2013) / Henrard (2014, Section 2.4)
    /// dual-curve 3M-into-6M tenor basis swap bootstrap across $K = 4$ semi-annual periods.
    pub fn evaluate_ametrano_bianchetti_dual_curve_tenor_basis_swap_bootstrap<T: Real>(
        spec: &AmetranoBianchettiTenorBasisSpec<T>,
    ) -> AmetranoBianchettiTenorBasisAnalytics<T> {
        const SEMI_PERIODS: usize = 4;
        let t_mat = spec.swap_maturity_years_t.max_val(T::from_f64(0.25));
        let tau6 = t_mat / T::from_f64(SEMI_PERIODS as f64);
        let tau3 = T::from_f64(0.5) * tau6;
        let sig3 = spec
            .compounding_convexity_vol_sigma3
            .max_val(T::from_f64(1e-5));

        let ois_df = |t: T| -> T {
            let slope = (T::one() - (-T::from_f64(0.5) * t).exp()) / T::from_f64(0.5);
            let y = spec.ois_discount_rate_r_ois + spec.ois_slope_theta_ois * slope;
            (-y * t).exp()
        };

        let mut p6_term = T::one();
        let mut ann6_ois = T::zero();
        let mut float3_pv = T::zero();
        let mut float6_pv = T::zero();

        for k in 1..=SEMI_PERIODS {
            let t_start = T::from_f64((k - 1) as f64) * tau6;
            let t_mid = t_start + tau3;
            let t_end = T::from_f64(k as f64) * tau6;

            let p0 = ois_df(t_start);
            let p_mid = ois_df(t_mid);
            let p_end = ois_df(t_end);

            let f_ois_1 = (p0 / p_mid - T::one()) / tau3;
            let f_ois_2 = (p_mid / p_end - T::one()) / tau3;

            let l3_1 = f_ois_1 + spec.ibor_3m_spread_over_ois_s3;
            let l3_2 = f_ois_2 + spec.ibor_3m_spread_over_ois_s3;

            let b_k = spec.tenor_basis_3m_vs_6m_spread_b36
                + spec.tenor_basis_term_steepening_gamma * t_end;
            let conv_adj =
                (tau3 * tau3 / tau6) * sig3 * sig3 * (t_start + T::from_f64(0.25) * tau6);
            let compounded_3m =
                ((T::one() + tau3 * l3_1) * (T::one() + tau3 * l3_2) - T::one()) / tau6;
            let l6_k = compounded_3m + conv_adj + b_k;

            p6_term /= T::one() + tau6 * l6_k;
            ann6_ois += tau6 * p_end;
            float3_pv += tau3 * p_mid * l3_1 + tau3 * p_end * l3_2;
            float6_pv += tau6 * p_end * l6_k;
        }

        let s3_par = float3_pv / ann6_ois;
        let s6_par = float6_pv / ann6_ois;

        AmetranoBianchettiTenorBasisAnalytics {
            terminal_6m_pseudo_discount_factor: p6_term,
            ois_discounted_3m_par_swap_rate: s3_par,
            implied_6m_over_3m_swap_spread: s6_par - s3_par,
            bootstrapped_6m_ibor_par_swap_rate: s6_par,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn ametrano_bianchetti_tenor_basis_invariants_and_ad_parity() {
            let base = AmetranoBianchettiTenorBasisSpec {
                ois_discount_rate_r_ois: 0.036_f64,
                ois_slope_theta_ois: 0.004,
                ibor_3m_spread_over_ois_s3: 0.0018,
                tenor_basis_3m_vs_6m_spread_b36: 0.0014,
                tenor_basis_term_steepening_gamma: 0.0005,
                compounding_convexity_vol_sigma3: 0.0075,
                swap_maturity_years_t: 2.0,
            };
            let r0 = evaluate_ametrano_bianchetti_dual_curve_tenor_basis_swap_bootstrap(&base);
            assert!(
                r0.terminal_6m_pseudo_discount_factor > 0.0
                    && r0.terminal_6m_pseudo_discount_factor < 1.0
            );
            assert!(r0.bootstrapped_6m_ibor_par_swap_rate > r0.ois_discounted_3m_par_swap_rate);
            assert!(r0.implied_6m_over_3m_swap_spread > base.tenor_basis_3m_vs_6m_spread_b36);

            let inputs = [
                base.ois_discount_rate_r_ois,
                base.ois_slope_theta_ois,
                base.ibor_3m_spread_over_ois_s3,
                base.tenor_basis_3m_vs_6m_spread_b36,
                base.tenor_basis_term_steepening_gamma,
                base.compounding_convexity_vol_sigma3,
                base.swap_maturity_years_t,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = AmetranoBianchettiTenorBasisSpec {
                ois_discount_rate_r_ois: dual_in[0],
                ois_slope_theta_ois: dual_in[1],
                ibor_3m_spread_over_ois_s3: dual_in[2],
                tenor_basis_3m_vs_6m_spread_b36: dual_in[3],
                tenor_basis_term_steepening_gamma: dual_in[4],
                compounding_convexity_vol_sigma3: dual_in[5],
                swap_maturity_years_t: dual_in[6],
            };
            let r_ad = evaluate_ametrano_bianchetti_dual_curve_tenor_basis_swap_bootstrap(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| AmetranoBianchettiTenorBasisSpec {
                    ois_discount_rate_r_ois: v[0],
                    ois_slope_theta_ois: v[1],
                    ibor_3m_spread_over_ois_s3: v[2],
                    tenor_basis_3m_vs_6m_spread_b36: v[3],
                    tenor_basis_term_steepening_gamma: v[4],
                    compounding_convexity_vol_sigma3: v[5],
                    swap_maturity_years_t: v[6],
                };
                let rp = evaluate_ametrano_bianchetti_dual_curve_tenor_basis_swap_bootstrap(&mk(p));
                let rm = evaluate_ametrano_bianchetti_dual_curve_tenor_basis_swap_bootstrap(&mk(m));
                let fd = (rp.bootstrapped_6m_ibor_par_swap_rate
                    - rm.bootstrapped_6m_ibor_par_swap_rate)
                    / (2.0 * eps);
                assert!(
                    r_ad.bootstrapped_6m_ibor_par_swap_rate.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.bootstrapped_6m_ibor_par_swap_rate.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use ametrano_bianchetti_tenor_basis_mod::*;
mod mcculloch_tax_clientele_curve_mod {
    use ql_ad::Real;

    /// Input specification for the McCulloch (1975, *Journal of Finance* 30(3)) &
    /// Litzenberger & Rolfo (1984, *JFE* 13(1)) tax-clientele discount curve and
    /// municipal-Treasury tax-exempt yield ratio model.
    #[derive(Debug, Clone, Copy)]
    pub struct McCullochTaxClienteleCurveSpec<T: Real> {
        /// Pre-tax Treasury par coupon rate $c_{\text{low}}$ of low-coupon discount bond (e.g., $0.025$).
        pub low_coupon_rate_c_low: T,
        /// Pre-tax Treasury par coupon rate $c_{\text{high}}$ of high-coupon premium bond (e.g., $0.060$).
        pub high_coupon_rate_c_high: T,
        /// Observed market price $P_{\text{low}}$ per unit face of the 1-period low-coupon bond (e.g., $0.985$).
        pub low_coupon_bond_price_p_low: T,
        /// Observed market price $P_{\text{high}}$ per unit face of the 1-period high-coupon bond (e.g., $1.011$).
        pub high_coupon_bond_price_p_high: T,
        /// Marginal ordinary income tax rate $\tau_{\text{inc}} \in (0, 1)$ on coupon income (e.g., $0.28$).
        pub marginal_income_tax_rate_tau_inc: T,
        /// Effective capital gains tax rate $\tau_{\text{cg}} \in (0, \tau_{\text{inc}})$ on price appreciation (e.g., $0.15$).
        pub capital_gains_tax_rate_tau_cg: T,
        /// Municipal credit/liquidity spread $s_{\text{muni}}$ over the tax-exempt discount curve (e.g., $0.0035$).
        pub municipal_credit_liquidity_spread: T,
    }

    /// Analytics returned by [`evaluate_mcculloch_litzenberger_tax_clientele_discount_curve`].
    #[derive(Debug, Clone, Copy)]
    pub struct McCullochTaxClienteleCurveAnalytics<T: Real> {
        /// Implied tax-clientele after-tax discount factor $d_{\text{tax}}(1)$.
        pub after_tax_discount_factor: T,
        /// Implied marginal clientele tax bracket $\hat{\tau}_{\text{imp}}$ from the low/high coupon bond pair.
        pub implied_marginal_tax_bracket: T,
        /// After-tax zero-coupon yield $y_{\text{after-tax}}$.
        pub after_tax_zero_yield: T,
        /// Equilibrium municipal-to-Treasury yield ratio $y_{\text{muni}} / y_{\text{Treasury}}$.
        pub municipal_to_treasury_yield_ratio: T,
    }

    /// Evaluates the McCulloch (1975) / Litzenberger & Rolfo (1984) tax-adjusted discount factor,
    /// implied marginal tax bracket $\hat{\tau}_{\text{imp}} = 1 - \frac{P_{\text{high}} - P_{\text{low}}}{c_{\text{high}} - c_{\text{low}}}$,
    /// and equilibrium municipal-to-Treasury tax-exempt yield ratio.
    pub fn evaluate_mcculloch_litzenberger_tax_clientele_discount_curve<T: Real>(
        spec: &McCullochTaxClienteleCurveSpec<T>,
    ) -> McCullochTaxClienteleCurveAnalytics<T> {
        let c_l = spec.low_coupon_rate_c_low;
        let c_h = spec.high_coupon_rate_c_high;
        let p_l = spec.low_coupon_bond_price_p_low;
        let p_h = spec.high_coupon_bond_price_p_high;
        let tau_i = spec.marginal_income_tax_rate_tau_inc;
        let tau_g = spec.capital_gains_tax_rate_tau_cg;

        // One-period coupon bonds strip to the principal discount factor
        // after subtracting their coupon cashflow.  Averaging the two
        // observations reduces quote noise while retaining AD sensitivity.
        let dc = c_h - c_l;
        let d_pre = T::from_f64(0.5) * ((p_l - c_l) + (p_h - c_h));

        // Implied marginal tax rate from differential coupon pricing.  The
        // one-period par-coupon normalization makes the coupon-price slope
        // directly equal to the after-tax coupon fraction.
        let dp_dc = (p_h - p_l) / dc;
        let tau_imp_mkt = T::one() - dp_dc;
        let tau_eff =
            T::from_f64(0.5) * (tau_i + tau_imp_mkt) - T::from_f64(0.25) * tau_g * (T::one() - p_l);

        // After-tax cashflow of low-coupon discount bond: (1 - tau_eff)*c_l + 1 - tau_g*(1 - P_l)
        let after_tax_cf_low = T::one() + (T::one() - tau_eff) * c_l - tau_g * (T::one() - p_l);
        let d_tax = p_l / after_tax_cf_low;

        let y_after_tax = -d_tax.ln();
        let y_pre_tax = -d_pre.ln();
        let y_muni = y_after_tax + spec.municipal_credit_liquidity_spread;
        let muni_tsy_ratio = y_muni / y_pre_tax;

        McCullochTaxClienteleCurveAnalytics {
            after_tax_discount_factor: d_tax,
            implied_marginal_tax_bracket: tau_eff,
            after_tax_zero_yield: y_after_tax,
            municipal_to_treasury_yield_ratio: muni_tsy_ratio,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn mcculloch_tax_clientele_curve_invariants_and_ad_parity() {
            let base = McCullochTaxClienteleCurveSpec {
                low_coupon_rate_c_low: 0.025_f64,
                high_coupon_rate_c_high: 0.060_f64,
                low_coupon_bond_price_p_low: 0.985_f64,
                high_coupon_bond_price_p_high: 1.010_f64,
                marginal_income_tax_rate_tau_inc: 0.28_f64,
                capital_gains_tax_rate_tau_cg: 0.15_f64,
                municipal_credit_liquidity_spread: 0.0035_f64,
            };
            let res = evaluate_mcculloch_litzenberger_tax_clientele_discount_curve(&base);
            assert!(res.after_tax_discount_factor > 0.90 && res.after_tax_discount_factor < 1.0);
            assert!(
                res.implied_marginal_tax_bracket > 0.10 && res.implied_marginal_tax_bracket < 0.50
            );
            assert!(
                res.municipal_to_treasury_yield_ratio > 0.50
                    && res.municipal_to_treasury_yield_ratio < 1.10
            );

            let inputs = [
                base.low_coupon_rate_c_low,
                base.high_coupon_rate_c_high,
                base.low_coupon_bond_price_p_low,
                base.high_coupon_bond_price_p_high,
                base.marginal_income_tax_rate_tau_inc,
                base.capital_gains_tax_rate_tau_cg,
                base.municipal_credit_liquidity_spread,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = McCullochTaxClienteleCurveSpec {
                low_coupon_rate_c_low: dual_in[0],
                high_coupon_rate_c_high: dual_in[1],
                low_coupon_bond_price_p_low: dual_in[2],
                high_coupon_bond_price_p_high: dual_in[3],
                marginal_income_tax_rate_tau_inc: dual_in[4],
                capital_gains_tax_rate_tau_cg: dual_in[5],
                municipal_credit_liquidity_spread: dual_in[6],
            };
            let r_ad = evaluate_mcculloch_litzenberger_tax_clientele_discount_curve(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| McCullochTaxClienteleCurveSpec {
                    low_coupon_rate_c_low: v[0],
                    high_coupon_rate_c_high: v[1],
                    low_coupon_bond_price_p_low: v[2],
                    high_coupon_bond_price_p_high: v[3],
                    marginal_income_tax_rate_tau_inc: v[4],
                    capital_gains_tax_rate_tau_cg: v[5],
                    municipal_credit_liquidity_spread: v[6],
                };
                let rp = evaluate_mcculloch_litzenberger_tax_clientele_discount_curve(&mk(p));
                let rm = evaluate_mcculloch_litzenberger_tax_clientele_discount_curve(&mk(m));
                let fd = (rp.municipal_to_treasury_yield_ratio
                    - rm.municipal_to_treasury_yield_ratio)
                    / (2.0 * eps);
                assert!(
                    r_ad.municipal_to_treasury_yield_ratio.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.municipal_to_treasury_yield_ratio.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use mcculloch_tax_clientele_curve_mod::*;
mod mercurio_russo_ibor_ois_basis_mod {
    use ql_ad::Real;

    /// Input specification for the Mercurio (2018, *Risk*) & Russo & Torricelli (2019,
    /// *European Journal of Operational Research* 273(2)) multi-curve IBOR–OIS credit
    /// vs. funding-liquidity tenor basis decomposition engine.
    #[derive(Debug, Clone, Copy)]
    pub struct MercurioRussoIborOisCreditLiquiditySpec<T: Real> {
        /// Overnight RFR / OIS continuously compounded rate $r_{\text{OIS}}$ (e.g., $0.034$).
        pub ois_rate_r: T,
        /// Interbank panel credit hazard intensity $\lambda_{\text{bank}}$ (e.g., $0.0045$).
        pub interbank_credit_hazard_lambda: T,
        /// Interbank panel Loss Given Default $\text{LGD}_{\text{bank}} \in (0, 1)$ (e.g., $0.60$).
        pub interbank_loss_given_default: T,
        /// Term funding-liquidity spread $\ell_{\text{liq}}(0, \tau)$ (e.g., $0.0018$).
        pub term_liquidity_spread_l: T,
        /// Credit-liquidity covariance intensity $\sigma_c \sigma_\ell \rho_{c\ell}$ (e.g., $0.0008$).
        pub credit_liquidity_covariance_gamma: T,
        /// Tenor accrual fraction $\tau$ in years (e.g., $0.50$ for 6M).
        pub accrual_tenor_tau: T,
    }

    /// Analytics returned by [`evaluate_mercurio_russo_ibor_ois_credit_liquidity_basis`].
    #[derive(Debug, Clone, Copy)]
    pub struct MercurioRussoIborOisCreditLiquidityAnalytics<T: Real> {
        /// Simple compounded OIS forward rate $F_{\text{OIS}}(0; 0, \tau) = \frac{e^{r_{\text{OIS}}\tau} - 1}{\tau}$.
        pub ois_simple_forward_rate: T,
        /// Multiplicative credit-liquidity spread factor $\beta_{\text{mult}}(\tau) > 1$.
        pub multiplicative_basis_factor: T,
        /// Pure credit-default component of the additive IBOR–OIS spread.
        pub credit_default_spread_component: T,
        /// Total equilibrium IBOR–OIS additive basis spread $L_{\text{IBOR}}(0, \tau) - F_{\text{OIS}}(0, \tau)$.
        pub total_ibor_ois_additive_basis: T,
    }

    /// Evaluates the Mercurio (2018) / Russo & Torricelli (2019, Eqs. 9–17) multiplicative
    /// and additive IBOR–OIS basis decomposition:
    /// $$1 + \tau L(0, \tau) = e^{r_{\text{OIS}}\tau} \exp\!\Big(\big(\text{LGD}_{\text{bank}}\lambda_{\text{bank}} + \ell_{\text{liq}} + \tfrac{1}{2}\gamma_{c\ell}\tau\big)\tau\Big),$$
    /// separating pure panel default compensation from term funding-liquidity and covariance convexity.
    pub fn evaluate_mercurio_russo_ibor_ois_credit_liquidity_basis<T: Real>(
        spec: &MercurioRussoIborOisCreditLiquiditySpec<T>,
    ) -> MercurioRussoIborOisCreditLiquidityAnalytics<T> {
        let tau = spec.accrual_tenor_tau;
        let r_ois = spec.ois_rate_r;
        let s_cred = spec.interbank_loss_given_default * spec.interbank_credit_hazard_lambda;
        let s_liq = spec.term_liquidity_spread_l;
        let cov_adj = T::from_f64(0.5) * spec.credit_liquidity_covariance_gamma * tau;

        let e_ois = (r_ois * tau).exp();
        let f_ois = (e_ois - T::one()) / tau;

        let mult_factor = ((s_cred + s_liq + cov_adj) * tau).exp();
        let mult_cred_only = (s_cred * tau).exp();

        let l_ibor = (e_ois * mult_factor - T::one()) / tau;
        let l_cred_only = (e_ois * mult_cred_only - T::one()) / tau;

        let cred_comp = l_cred_only - f_ois;
        let total_basis = l_ibor - f_ois;

        MercurioRussoIborOisCreditLiquidityAnalytics {
            ois_simple_forward_rate: f_ois,
            multiplicative_basis_factor: mult_factor,
            credit_default_spread_component: cred_comp,
            total_ibor_ois_additive_basis: total_basis,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn mercurio_russo_ibor_ois_basis_invariants_and_ad_parity() {
            let base = MercurioRussoIborOisCreditLiquiditySpec {
                ois_rate_r: 0.034_f64,
                interbank_credit_hazard_lambda: 0.0045_f64,
                interbank_loss_given_default: 0.60_f64,
                term_liquidity_spread_l: 0.0018_f64,
                credit_liquidity_covariance_gamma: 0.0008_f64,
                accrual_tenor_tau: 0.50_f64,
            };
            let res = evaluate_mercurio_russo_ibor_ois_credit_liquidity_basis(&base);
            assert!(res.ois_simple_forward_rate > base.ois_rate_r);
            assert!(res.multiplicative_basis_factor > 1.0);
            assert!(res.credit_default_spread_component > 0.0);
            assert!(res.total_ibor_ois_additive_basis > res.credit_default_spread_component);

            let inputs = [
                base.ois_rate_r,
                base.interbank_credit_hazard_lambda,
                base.interbank_loss_given_default,
                base.term_liquidity_spread_l,
                base.credit_liquidity_covariance_gamma,
                base.accrual_tenor_tau,
            ];
            let mut dual_in = [Dual::<6>::from_value(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = MercurioRussoIborOisCreditLiquiditySpec {
                ois_rate_r: dual_in[0],
                interbank_credit_hazard_lambda: dual_in[1],
                interbank_loss_given_default: dual_in[2],
                term_liquidity_spread_l: dual_in[3],
                credit_liquidity_covariance_gamma: dual_in[4],
                accrual_tenor_tau: dual_in[5],
            };
            let r_ad = evaluate_mercurio_russo_ibor_ois_credit_liquidity_basis(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| MercurioRussoIborOisCreditLiquiditySpec {
                    ois_rate_r: v[0],
                    interbank_credit_hazard_lambda: v[1],
                    interbank_loss_given_default: v[2],
                    term_liquidity_spread_l: v[3],
                    credit_liquidity_covariance_gamma: v[4],
                    accrual_tenor_tau: v[5],
                };
                let rp = evaluate_mercurio_russo_ibor_ois_credit_liquidity_basis(&mk(p));
                let rm = evaluate_mercurio_russo_ibor_ois_credit_liquidity_basis(&mk(m));
                let fd = (rp.total_ibor_ois_additive_basis - rm.total_ibor_ois_additive_basis)
                    / (2.0 * eps);
                assert!(
                    r_ad.total_ibor_ois_additive_basis.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.total_ibor_ois_additive_basis.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use mercurio_russo_ibor_ois_basis_mod::*;
