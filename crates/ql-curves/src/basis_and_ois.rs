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
