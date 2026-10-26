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
mod phase33_tests {
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
    fn phase33_models_smoke_and_dual_fd() {
        let s0 = evaluate_mascolell_cochrane_arrow_debreu_sdf(&sdf_spec(0.0));
        assert!(s0.q_measure_conservation_error < 1e-12);
        assert!(s0.hansen_jagannathan_bound_slack >= -1e-12 && s0.ccapm_equity_risk_premium > 0.0);
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
            - evaluate_mascolell_cochrane_arrow_debreu_sdf(&sdf_spec(-h)).equilibrium_equity_price)
            / (2. * h);
        let ad_sdf =
            evaluate_mascolell_cochrane_arrow_debreu_sdf(&sdf_spec(Dual::<1>::variable(0., 0)))
                .equilibrium_equity_price;
        assert!((ad_sdf.d[0] - fd_sdf).abs() < 1e-5);
    }
}
