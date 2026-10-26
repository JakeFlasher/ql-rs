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
    let lq = s.abandonment_salvage_value;
    let pd1 = ((lq / v).ln() + (s.risk_free_rate - q + T::from_f64(0.5) * vol * vol) * t)
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
    let er = s
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
    let cov = ep - p * er;
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
        .map(|(i, pr)| *pr * (s.equity_payoff_states[i] - er) * (s.equity_payoff_states[i] - er))
        .fold(T::zero(), |a, x| a + x)
        .sqrt();
    let sharpe = (er - rf).abs() / vr;
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
        q_measure_conservation_error: (q.iter().copied().fold(T::zero(), |a, x| a + x) - T::one())
            .abs(),
    }
}
#[cfg(test)]
mod phase33_tests {
    use super::*;
    #[test]
    fn phase33_models_smoke() {
        let s = MasColellCochraneArrowDebreuSdfSpec {
            state_probabilities: [0.3, 0.4, 0.3],
            consumption_growth_ratios: [0.98, 1.02, 1.05],
            equity_payoff_states: [0.8, 1.1, 1.4],
            subjective_discount_beta: 0.99,
            crra_risk_aversion_gamma: 2.0,
        };
        assert!(
            evaluate_mascolell_cochrane_arrow_debreu_sdf(&s).q_measure_conservation_error < 1e-12
        );
    }
}
