mod evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation_module {
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct PenmanReoiOhlsonJuettnerSpec<T: Real = f64> {
        pub net_operating_assets_0: T,
        pub net_financial_obligations_0: T,
        pub operating_income_after_tax_1: T,
        pub operating_income_after_tax_2: T,
        pub net_financial_expense_after_tax_1: T,
        pub wacc: T,
        pub terminal_reoi_growth: T,
        pub current_share_price: T,
        pub forward_eps_1: T,
        pub forward_eps_2: T,
        pub forward_dps_1: T,
        pub perpetual_aeg_growth: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct PenmanReoiOhlsonJuettnerReport<T: Real = f64> {
        pub return_on_net_operating_assets_rnoa: T,
        pub net_borrowing_cost_nbc: T,
        pub financial_leverage_flev: T,
        pub operating_spread: T,
        pub penman_roe: T,
        pub dupont_identity_conservation_error: T,
        pub residual_operating_income_1: T,
        pub reoi_enterprise_value: T,
        pub reoi_intrinsic_equity_value: T,
        pub ohlson_juettner_implied_cost_of_equity: T,
    }

    pub fn evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation<T: Real>(
        s: &PenmanReoiOhlsonJuettnerSpec<T>,
    ) -> PenmanReoiOhlsonJuettnerReport<T> {
        let noa = s.net_operating_assets_0.max_val(T::from_f64(1e-12));
        let nfo = s.net_financial_obligations_0;
        let cse = (noa - nfo).max_val(T::from_f64(1e-12));
        let rnoa = s.operating_income_after_tax_1 / noa;
        let nbc = s.net_financial_expense_after_tax_1 / nfo.abs().max_val(T::from_f64(1e-12));
        let flev = nfo / cse;
        let spread = rnoa - nbc;
        let roe = rnoa + flev * spread;
        let roe_direct =
            (s.operating_income_after_tax_1 - s.net_financial_expense_after_tax_1) / cse;
        let reoi1 = s.operating_income_after_tax_1 - s.wacc * noa;
        let noa1 = noa * (T::one() + s.terminal_reoi_growth);
        let reoi2 = s.operating_income_after_tax_2 - s.wacc * noa1;
        let rw = s.wacc.max_val(T::from_f64(1e-8));
        let ev = noa
            + reoi1 / (T::one() + rw)
            + reoi2 / ((T::one() + rw) * (rw - s.terminal_reoi_growth).max_val(T::from_f64(1e-8)));
        let g = s.perpetual_aeg_growth;
        let a = (g + s.forward_dps_1 / s.current_share_price.max_val(T::from_f64(1e-12)))
            / T::from_f64(2.0);
        let growth_component =
            (s.forward_eps_2 - s.forward_eps_1) / s.forward_eps_1.max_val(T::from_f64(1e-12)) - g;
        let disc = (a * a
            + s.forward_eps_1 / s.current_share_price.max_val(T::from_f64(1e-12))
                * growth_component)
            .max_val(T::zero())
            .sqrt();
        PenmanReoiOhlsonJuettnerReport {
            return_on_net_operating_assets_rnoa: rnoa,
            net_borrowing_cost_nbc: nbc,
            financial_leverage_flev: flev,
            operating_spread: spread,
            penman_roe: roe,
            dupont_identity_conservation_error: (roe_direct - roe).abs(),
            residual_operating_income_1: reoi1,
            reoi_enterprise_value: ev,
            reoi_intrinsic_equity_value: ev - nfo,
            ohlson_juettner_implied_cost_of_equity: a + disc,
        }
    }

    #[cfg(test)]
    mod domain_suite_43_tests {
        use super::*;
        use ql_ad::Dual;

        fn spec<T: Real>(x: T) -> PenmanReoiOhlsonJuettnerSpec<T> {
            PenmanReoiOhlsonJuettnerSpec {
                net_operating_assets_0: T::from_f64(100.0) + x,
                net_financial_obligations_0: T::from_f64(30.0),
                operating_income_after_tax_1: T::from_f64(12.0),
                operating_income_after_tax_2: T::from_f64(13.0),
                net_financial_expense_after_tax_1: T::from_f64(2.0),
                wacc: T::from_f64(0.06),
                terminal_reoi_growth: T::from_f64(0.02),
                current_share_price: T::from_f64(100.0),
                forward_eps_1: T::from_f64(5.0),
                forward_eps_2: T::from_f64(5.5),
                forward_dps_1: T::from_f64(2.0),
                perpetual_aeg_growth: T::from_f64(0.03),
            }
        }

        #[test]
        fn suite43_invariants_and_dual_fd() {
            let r = evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation(&spec(0.0));
            assert!(r.dupont_identity_conservation_error < 1e-12);
            assert!(r.ohlson_juettner_implied_cost_of_equity > 0.03);
            let h = 1e-6;
            let fd = (evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation(&spec(h))
                .reoi_enterprise_value
                - evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation(&spec(-h))
                    .reoi_enterprise_value)
                / (2.0 * h);
            let ad = evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation(&spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .reoi_enterprise_value;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation_module::*;

mod evaluate_damodaran_rd_capitalization_and_patent_real_option_module {
    use ql_ad::Real;
    #[derive(Clone, Debug)]
    pub struct DamodaranRdCapitalizationRealOptionSpec<T: Real = f64> {
        pub reported_ebit: T,
        pub reported_invested_capital: T,
        pub effective_tax_rate: T,
        pub current_rd_expense: T,
        pub historical_rd_expenses: Vec<T>,
        pub patent_project_pv: T,
        pub patent_development_cost: T,
        pub patent_life_years: T,
        pub project_volatility: T,
        pub risk_free_rate: T,
        pub abandonment_salvage_value: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct DamodaranRdCapitalizationRealOptionReport<T: Real = f64> {
        pub unamortized_research_asset: T,
        pub current_rd_amortization: T,
        pub adjusted_operating_ebit: T,
        pub adjusted_invested_capital: T,
        pub reported_after_tax_roic: T,
        pub rd_adjusted_after_tax_roic: T,
        pub delay_cost_dividend_yield: T,
        pub patent_commercialization_option_value: T,
        pub abandonment_put_option_value: T,
        pub total_strategic_expanded_npv: T,
    }
    pub fn evaluate_damodaran_rd_capitalization_and_patent_real_option<T: Real>(
        s: &DamodaranRdCapitalizationRealOptionSpec<T>,
    ) -> DamodaranRdCapitalizationRealOptionReport<T> {
        let l = s.historical_rd_expenses.len();
        let mut a = s.current_rd_expense;
        for (k, v) in s.historical_rd_expenses.iter().enumerate() {
            a += *v * T::from_f64((l - k - 1) as f64) / (l as f64);
        }
        let amort = s
            .historical_rd_expenses
            .iter()
            .copied()
            .fold(T::zero(), |x, y| x + y)
            / (l as f64);
        let e = s.reported_ebit + s.current_rd_expense - amort;
        let ic = s.reported_invested_capital + a;
        let rat = T::one() - s.effective_tax_rate;
        let q = T::one() / s.patent_life_years;
        let v = s.patent_project_pv;
        let k = s.patent_development_cost;
        let t = s.patent_life_years;
        let vol = s.project_volatility;
        let d1 = ((v / k).ln() + (s.risk_free_rate - q + T::from_f64(0.5) * vol * vol) * t)
            / (vol * t.sqrt());
        let d2 = d1 - vol * t.sqrt();
        let call =
            v * (-q * t).exp() * d1.norm_cdf() - k * (-s.risk_free_rate * t).exp() * d2.norm_cdf();
        let lq = s.abandonment_salvage_value.max_val(T::from_f64(1e-8));
        let pd1 = ((v / lq).ln() + (s.risk_free_rate - q + T::from_f64(0.5) * vol * vol) * t)
            / (vol * t.sqrt());
        let pd2 = pd1 - vol * t.sqrt();
        let put = lq * (-s.risk_free_rate * t).exp() * (-pd2).norm_cdf()
            - v * (-q * t).exp() * (-pd1).norm_cdf();
        DamodaranRdCapitalizationRealOptionReport {
            unamortized_research_asset: a,
            current_rd_amortization: amort,
            adjusted_operating_ebit: e,
            adjusted_invested_capital: ic,
            reported_after_tax_roic: rat * s.reported_ebit / s.reported_invested_capital,
            rd_adjusted_after_tax_roic: rat * e / ic,
            delay_cost_dividend_yield: q,
            patent_commercialization_option_value: call,
            abandonment_put_option_value: put,
            total_strategic_expanded_npv: call + put,
        }
    }
    #[derive(Clone, Copy, Debug)]
    pub struct MasColellCochraneArrowDebreuSdfSpec<T: Real = f64> {
        pub state_probabilities: [T; 3],
        pub consumption_growth_ratios: [T; 3],
        pub equity_payoff_states: [T; 3],
        pub subjective_discount_beta: T,
        pub crra_risk_aversion_gamma: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct MasColellCochraneArrowDebreuSdfReport<T: Real = f64> {
        pub arrow_debreu_state_prices: [T; 3],
        pub risk_neutral_probabilities_q: [T; 3],
        pub risk_free_gross_rate: T,
        pub equilibrium_equity_price: T,
        pub expected_physical_equity_return: T,
        pub ccapm_equity_risk_premium: T,
        pub sdf_volatility_to_mean_ratio: T,
        pub equity_sharpe_ratio: T,
        pub hansen_jagannathan_bound_slack: T,
        pub q_measure_conservation_error: T,
    }
    pub fn evaluate_mascolell_cochrane_arrow_debreu_sdf<T: Real>(
        s: &MasColellCochraneArrowDebreuSdfSpec<T>,
    ) -> MasColellCochraneArrowDebreuSdfReport<T> {
        let mut m = [T::zero(); 3];
        let mut psi = [T::zero(); 3];
        for i in 0..3 {
            m[i] = s.subjective_discount_beta
                * s.consumption_growth_ratios[i].pow(-s.crra_risk_aversion_gamma);
            psi[i] = s.state_probabilities[i] * m[i];
        }
        let p = psi.iter().copied().fold(T::zero(), |a, x| a + x);
        let rf = T::one() / p;
        let mut q = [T::zero(); 3];
        for i in 0..3 {
            q[i] = psi[i] / p;
        }
        let ep_raw = s
            .state_probabilities
            .iter()
            .enumerate()
            .map(|(i, pr)| *pr * s.equity_payoff_states[i])
            .fold(T::zero(), |a, x| a + x);
        let ep = s
            .state_probabilities
            .iter()
            .enumerate()
            .map(|(i, pr)| *pr * m[i] * s.equity_payoff_states[i])
            .fold(T::zero(), |a, x| a + x);
        let ep_safe = ep.max_val(T::from_f64(1e-12));
        let er = ep_raw / ep_safe;
        let cov = T::one() - p * er;
        let prem = -rf * cov;
        let vm = s
            .state_probabilities
            .iter()
            .enumerate()
            .map(|(i, pr)| *pr * (m[i] - p) * (m[i] - p))
            .fold(T::zero(), |a, x| a + x)
            .sqrt();
        let vr = s
            .state_probabilities
            .iter()
            .enumerate()
            .map(|(i, pr)| {
                let ri = s.equity_payoff_states[i] / ep_safe;
                *pr * (ri - er) * (ri - er)
            })
            .fold(T::zero(), |a, x| a + x)
            .sqrt()
            .max_val(T::from_f64(1e-12));
        let sharpe = ((er - rf) * (er - rf) + T::from_f64(1e-30)).sqrt() / vr;
        MasColellCochraneArrowDebreuSdfReport {
            arrow_debreu_state_prices: psi,
            risk_neutral_probabilities_q: q,
            risk_free_gross_rate: rf,
            equilibrium_equity_price: ep,
            expected_physical_equity_return: er,
            ccapm_equity_risk_premium: prem,
            sdf_volatility_to_mean_ratio: vm / p,
            equity_sharpe_ratio: sharpe,
            hansen_jagannathan_bound_slack: vm / p - sharpe,
            q_measure_conservation_error: {
                let q_err = q.iter().copied().fold(T::zero(), |a, x| a + x) - T::one();
                (q_err * q_err + T::from_f64(1e-30)).sqrt()
            },
        }
    }
    #[cfg(test)]
    mod domain_suite_45_tests {
        use super::*;
        use ql_ad::Dual;
        fn s<T: Real>(x: T) -> DamodaranRdCapitalizationRealOptionSpec<T> {
            DamodaranRdCapitalizationRealOptionSpec {
                reported_ebit: T::from_f64(10.),
                reported_invested_capital: T::from_f64(100.),
                effective_tax_rate: T::from_f64(0.25),
                current_rd_expense: T::from_f64(5.) + x,
                historical_rd_expenses: vec![T::from_f64(4.), T::from_f64(3.), T::from_f64(2.)],
                patent_project_pv: T::from_f64(20.),
                patent_development_cost: T::from_f64(10.),
                patent_life_years: T::from_f64(5.),
                project_volatility: T::from_f64(0.3),
                risk_free_rate: T::from_f64(0.03),
                abandonment_salvage_value: T::from_f64(8.),
            }
        }
        fn sdf_spec<T: Real>(x: T) -> MasColellCochraneArrowDebreuSdfSpec<T> {
            MasColellCochraneArrowDebreuSdfSpec {
                state_probabilities: [T::from_f64(0.3), T::from_f64(0.4), T::from_f64(0.3)],
                consumption_growth_ratios: [
                    T::from_f64(0.98),
                    T::from_f64(1.02) + x,
                    T::from_f64(1.05),
                ],
                equity_payoff_states: [T::from_f64(0.8), T::from_f64(1.1), T::from_f64(1.4)],
                subjective_discount_beta: T::from_f64(0.99),
                crra_risk_aversion_gamma: T::from_f64(2.0),
            }
        }
        #[test]
        fn suite45_models_smoke_and_dual_fd() {
            let s0 = evaluate_mascolell_cochrane_arrow_debreu_sdf(&sdf_spec(0.0));
            assert!(s0.q_measure_conservation_error < 1e-12);
            assert!(
                s0.hansen_jagannathan_bound_slack >= -1e-12 && s0.ccapm_equity_risk_premium > 0.0
            );
            let h = 1e-6;
            let fd = (evaluate_damodaran_rd_capitalization_and_patent_real_option(&s(h))
                .rd_adjusted_after_tax_roic
                - evaluate_damodaran_rd_capitalization_and_patent_real_option(&s(-h))
                    .rd_adjusted_after_tax_roic)
                / (2. * h);
            let ad = evaluate_damodaran_rd_capitalization_and_patent_real_option(&s(
                Dual::<1>::variable(0., 0),
            ))
            .rd_adjusted_after_tax_roic;
            assert!((ad.d[0] - fd).abs() < 1e-5);

            let fd_sdf = (evaluate_mascolell_cochrane_arrow_debreu_sdf(&sdf_spec(h))
                .equilibrium_equity_price
                - evaluate_mascolell_cochrane_arrow_debreu_sdf(&sdf_spec(-h))
                    .equilibrium_equity_price)
                / (2. * h);
            let ad_sdf =
                evaluate_mascolell_cochrane_arrow_debreu_sdf(&sdf_spec(Dual::<1>::variable(0., 0)))
                    .equilibrium_equity_price;
            assert!((ad_sdf.d[0] - fd_sdf).abs() < 1e-5);
        }
    }
}
pub use evaluate_damodaran_rd_capitalization_and_patent_real_option_module::*;

mod evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation_module {
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct KollerMckinseyKvdEvaValuationSpec<T: Real = f64> {
        pub invested_capital_0: T,
        pub noplat_1: T,
        pub ronic: T,
        pub wacc: T,
        pub perpetual_growth_g: T,
        pub net_debt: T,
        pub diluted_shares: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct KollerMckinseyKvdEvaValuationReport<T: Real = f64> {
        pub existing_roic_1: T,
        pub reinvestment_rate_ir: T,
        pub free_cash_flow_1: T,
        pub economic_profit_eva_1: T,
        pub kvd_dcf_enterprise_value: T,
        pub economic_profit_eva_enterprise_value: T,
        pub dcf_eva_equivalence_error: T,
        pub ev_to_invested_capital_multiple: T,
        pub intrinsic_share_price: T,
    }

    pub fn evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation<T: Real>(
        s: &KollerMckinseyKvdEvaValuationSpec<T>,
    ) -> KollerMckinseyKvdEvaValuationReport<T> {
        let roic1 = s.noplat_1 / s.invested_capital_0;
        let ir = s.perpetual_growth_g / s.ronic;
        let fcf1 = s.noplat_1 * (T::one() - ir);
        let spread = (s.wacc - s.perpetual_growth_g).max_val(T::from_f64(1e-8));
        let v_dcf = fcf1 / spread;
        let eva1 = s.invested_capital_0 * (roic1 - s.wacc);
        let v_eva = s.invested_capital_0
            + eva1 / s.wacc
            + (s.noplat_1 * ir * (s.ronic - s.wacc)) / (s.wacc * spread);
        let diff = v_dcf - v_eva;
        KollerMckinseyKvdEvaValuationReport {
            existing_roic_1: roic1,
            reinvestment_rate_ir: ir,
            free_cash_flow_1: fcf1,
            economic_profit_eva_1: eva1,
            kvd_dcf_enterprise_value: v_dcf,
            economic_profit_eva_enterprise_value: v_eva,
            dcf_eva_equivalence_error: (diff * diff + T::from_f64(1e-30)).sqrt(),
            ev_to_invested_capital_multiple: v_dcf / s.invested_capital_0,
            intrinsic_share_price: (v_dcf - s.net_debt) / s.diluted_shares,
        }
    }

    #[cfg(test)]
    mod domain_suite_47_tests {
        use super::*;
        use ql_ad::Dual;
        fn spec<T: Real>(x: T) -> KollerMckinseyKvdEvaValuationSpec<T> {
            KollerMckinseyKvdEvaValuationSpec {
                invested_capital_0: T::from_f64(500.0),
                noplat_1: T::from_f64(75.0) + x,
                ronic: T::from_f64(0.16),
                wacc: T::from_f64(0.09),
                perpetual_growth_g: T::from_f64(0.04),
                net_debt: T::from_f64(150.0),
                diluted_shares: T::from_f64(50.0),
            }
        }
        #[test]
        fn suite47_models_invariants_and_dual_fd() {
            let r = evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation(&spec(0.0));
            assert!(r.dcf_eva_equivalence_error < 1e-10);
            assert!(r.ev_to_invested_capital_multiple > 1.0);
            let h = 1e-6;
            let fd = (evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation(&spec(h))
                .intrinsic_share_price
                - evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation(&spec(-h))
                    .intrinsic_share_price)
                / (2.0 * h);
            let ad = evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation(&spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .intrinsic_share_price;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation_module::*;

mod evaluate_penman_nissim_dupont_and_reoi_valuation_module {
    #![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
    use ql_ad::Real;

    #[inline]
    fn smax<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[inline]
    fn smooth_ind<T: Real>(z: T) -> T {
        T::one() / (T::one() + (-z).exp())
    }

    #[derive(Clone, Copy, Debug)]
    pub struct PenmanNissimDuPontReoiValuationSpec<T: Real = f64> {
        pub net_operating_assets_0: T,
        pub net_financial_obligations_0: T,
        pub sales_revenue_1: T,
        pub after_tax_operating_income_1: T,
        pub after_tax_net_financial_expense_1: T,
        pub required_return_operations_rho: T,
        pub perpetual_growth_g: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct PenmanNissimDuPontReoiValuationReport<T: Real = f64> {
        pub operating_profit_margin_pm: T,
        pub asset_turnover_ato: T,
        pub return_on_net_operating_assets_rnoa: T,
        pub net_borrowing_cost_nbc: T,
        pub financial_leverage_flev: T,
        pub operating_spread: T,
        pub return_on_common_equity_roe: T,
        pub dupont_roe_conservation_error: T,
        pub residual_operating_income_1: T,
        pub reoi_enterprise_value: T,
        pub fcf_enterprise_value: T,
        pub reoi_fcf_equivalence_error: T,
        pub intrinsic_common_equity_value: T,
    }

    pub fn evaluate_penman_nissim_dupont_and_reoi_valuation<T: Real>(
        spec: &PenmanNissimDuPontReoiValuationSpec<T>,
    ) -> PenmanNissimDuPontReoiValuationReport<T> {
        let noa0 = spec.net_operating_assets_0.max_val(T::from_f64(1e-6));
        let nfo0 = spec.net_financial_obligations_0.max_val(T::from_f64(1e-6));
        let cse0 = (noa0 - nfo0).max_val(T::from_f64(1e-6));
        let sales1 = spec.sales_revenue_1.max_val(T::from_f64(1e-6));
        let oi1 = spec.after_tax_operating_income_1;
        let nfe1 = spec.after_tax_net_financial_expense_1;

        let pm = oi1 / sales1;
        let ato = sales1 / noa0;
        let rnoa = pm * ato;
        let nbc = nfe1 / nfo0;
        let flev = nfo0 / cse0;
        let spread = rnoa - nbc;
        let roe_direct = (oi1 - nfe1) / cse0;
        let roe_np = rnoa + flev * spread;
        let d_roe = roe_direct - roe_np;
        let roe_err = (d_roe * d_roe + T::from_f64(1e-30)).sqrt();

        let denom = (spec.required_return_operations_rho - spec.perpetual_growth_g)
            .max_val(T::from_f64(1e-6));
        let reoi1 = oi1 - spec.required_return_operations_rho * noa0;
        let v_reoi = noa0 + reoi1 / denom;
        let fcf1 = oi1 - spec.perpetual_growth_g * noa0;
        let v_fcf = fcf1 / denom;
        let d_val = v_reoi - v_fcf;
        let val_err = (d_val * d_val + T::from_f64(1e-30)).sqrt();

        PenmanNissimDuPontReoiValuationReport {
            operating_profit_margin_pm: pm,
            asset_turnover_ato: ato,
            return_on_net_operating_assets_rnoa: rnoa,
            net_borrowing_cost_nbc: nbc,
            financial_leverage_flev: flev,
            operating_spread: spread,
            return_on_common_equity_roe: roe_direct,
            dupont_roe_conservation_error: roe_err,
            residual_operating_income_1: reoi1,
            reoi_enterprise_value: v_reoi,
            fcf_enterprise_value: v_fcf,
            reoi_fcf_equivalence_error: val_err,
            intrinsic_common_equity_value: v_reoi - nfo0,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct LjungqvistSargentMccallSearchSpec<T: Real = f64> {
        pub offer_states: [T; 4],
        pub offer_probabilities: [T; 4],
        pub waiting_compensation_c: T,
        pub discount_factor_beta: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct LjungqvistSargentMccallSearchReport<T: Real = f64> {
        pub reservation_threshold_w_bar: T,
        pub bellman_lifetime_value: T,
        pub option_value_of_waiting: T,
        pub acceptance_probability_per_period: T,
        pub expected_stopping_duration_periods: T,
        pub bellman_reservation_residual: T,
    }

    pub fn evaluate_ljungqvist_sargent_mccall_bellman_search<T: Real>(
        spec: &LjungqvistSargentMccallSearchSpec<T>,
    ) -> LjungqvistSargentMccallSearchReport<T> {
        let p_sum = (spec.offer_probabilities[0]
            + spec.offer_probabilities[1]
            + spec.offer_probabilities[2]
            + spec.offer_probabilities[3])
            .max_val(T::from_f64(1e-8));
        let p = [
            spec.offer_probabilities[0] / p_sum,
            spec.offer_probabilities[1] / p_sum,
            spec.offer_probabilities[2] / p_sum,
            spec.offer_probabilities[3] / p_sum,
        ];
        let beta = spec.discount_factor_beta.min_val(T::from_f64(0.999));
        let one_m_b = (T::one() - beta).max_val(T::from_f64(1e-6));
        let mut w_bar = p[0] * spec.offer_states[0]
            + p[1] * spec.offer_states[1]
            + p[2] * spec.offer_states[2]
            + p[3] * spec.offer_states[3];

        for _ in 0..80 {
            let mut exp_max = T::zero();
            for i in 0..4 {
                exp_max += p[i] * smax(spec.offer_states[i], w_bar);
            }
            w_bar = one_m_b * spec.waiting_compensation_c + beta * exp_max;
        }

        let mut exp_max_final = T::zero();
        let mut exp_excess = T::zero();
        let mut acc_prob = T::zero();
        for i in 0..4 {
            exp_max_final += p[i] * smax(spec.offer_states[i], w_bar);
            exp_excess += p[i] * smax(spec.offer_states[i] - w_bar, T::zero());
            acc_prob += p[i] * smooth_ind(T::from_f64(30.0) * (spec.offer_states[i] - w_bar));
        }
        let rhs = one_m_b * spec.waiting_compensation_c + beta * exp_max_final;
        let d_res = w_bar - rhs;
        let res_err = (d_res * d_res + T::from_f64(1e-30)).sqrt();

        LjungqvistSargentMccallSearchReport {
            reservation_threshold_w_bar: w_bar,
            bellman_lifetime_value: w_bar / one_m_b,
            option_value_of_waiting: (beta / one_m_b) * exp_excess,
            acceptance_probability_per_period: acc_prob,
            expected_stopping_duration_periods: T::one() / acc_prob.max_val(T::from_f64(1e-6)),
            bellman_reservation_residual: res_err,
        }
    }

    #[cfg(test)]
    mod domain_suite_48_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn penman_and_mccall_models_invariants_and_ad_parity() {
            let h = 1e-6;

            let pen_spec = PenmanNissimDuPontReoiValuationSpec {
                net_operating_assets_0: 500.0,
                net_financial_obligations_0: 180.0,
                sales_revenue_1: 650.0,
                after_tax_operating_income_1: 72.0,
                after_tax_net_financial_expense_1: 9.0,
                required_return_operations_rho: 0.09,
                perpetual_growth_g: 0.03,
            };
            let pen_rep = evaluate_penman_nissim_dupont_and_reoi_valuation(&pen_spec);
            assert!(pen_rep.dupont_roe_conservation_error < 1e-12);
            assert!(pen_rep.reoi_fcf_equivalence_error < 1e-10);
            assert!(
                pen_rep.return_on_common_equity_roe > pen_rep.return_on_net_operating_assets_rnoa
            );

            let pen_ad = evaluate_penman_nissim_dupont_and_reoi_valuation(
                &PenmanNissimDuPontReoiValuationSpec {
                    net_operating_assets_0: Dual::<1>::from_value(pen_spec.net_operating_assets_0),
                    net_financial_obligations_0: Dual::<1>::from_value(
                        pen_spec.net_financial_obligations_0,
                    ),
                    sales_revenue_1: Dual::<1>::from_value(pen_spec.sales_revenue_1),
                    after_tax_operating_income_1: Dual::<1>::variable(
                        pen_spec.after_tax_operating_income_1,
                        0,
                    ),
                    after_tax_net_financial_expense_1: Dual::<1>::from_value(
                        pen_spec.after_tax_net_financial_expense_1,
                    ),
                    required_return_operations_rho: Dual::<1>::from_value(
                        pen_spec.required_return_operations_rho,
                    ),
                    perpetual_growth_g: Dual::<1>::from_value(pen_spec.perpetual_growth_g),
                },
            );
            let mut pp = pen_spec;
            pp.after_tax_operating_income_1 += h;
            let mut pm = pen_spec;
            pm.after_tax_operating_income_1 -= h;
            let _fd_pen = (evaluate_penman_nissim_dupont_and_reoi_valuation(&pp)
                .intrinsic_common_equity_value
                - evaluate_penman_nissim_dupont_and_reoi_valuation(&pm)
                    .intrinsic_common_equity_value)
                / (2.0 * h);
            assert!((pen_ad.intrinsic_common_equity_value.d[0] - _fd_pen).abs() < 1e-5);

            let mc_spec = LjungqvistSargentMccallSearchSpec {
                offer_states: [40.0, 55.0, 72.0, 95.0],
                offer_probabilities: [0.25, 0.35, 0.25, 0.15],
                waiting_compensation_c: 32.0,
                discount_factor_beta: 0.90,
            };
            let mc_rep = evaluate_ljungqvist_sargent_mccall_bellman_search(&mc_spec);
            assert!(mc_rep.bellman_reservation_residual < 1e-9);
            assert!(mc_rep.reservation_threshold_w_bar > mc_spec.waiting_compensation_c);
            assert!(mc_rep.expected_stopping_duration_periods >= 1.0);

            let mc_ad = evaluate_ljungqvist_sargent_mccall_bellman_search(
                &LjungqvistSargentMccallSearchSpec {
                    offer_states: mc_spec.offer_states.map(Dual::<1>::from_value),
                    offer_probabilities: mc_spec.offer_probabilities.map(Dual::<1>::from_value),
                    waiting_compensation_c: Dual::<1>::variable(mc_spec.waiting_compensation_c, 0),
                    discount_factor_beta: Dual::<1>::from_value(mc_spec.discount_factor_beta),
                },
            );
            let mut mp = mc_spec;
            mp.waiting_compensation_c += h;
            let mut mm = mc_spec;
            mm.waiting_compensation_c -= h;
            let _fd_mc = (evaluate_ljungqvist_sargent_mccall_bellman_search(&mp)
                .reservation_threshold_w_bar
                - evaluate_ljungqvist_sargent_mccall_bellman_search(&mm)
                    .reservation_threshold_w_bar)
                / (2.0 * h);
            assert!((mc_ad.reservation_threshold_w_bar.d[0] - _fd_mc).abs() < 1e-5);
        }
    }
}
pub use evaluate_penman_nissim_dupont_and_reoi_valuation_module::*;
