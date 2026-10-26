use ql_ad::Real;

trait InvNormCdf: Real {
    fn inv_norm_cdf(self) -> Self {
        let target = self;
        let mut q = Self::zero();
        for _ in 0..8 {
            q -= (q.norm_cdf() - target) / q.norm_pdf().max_val(Self::from_f64(1e-12));
        }
        q
    }
}
impl<T: Real> InvNormCdf for T {}

#[derive(Clone, Debug)]
pub struct ChristoffersenBouchaudVarBacktestAndCumulantSpec<T: Real = f64> {
    pub total_days_t: usize,
    pub count_n00: T,
    pub count_n01: T,
    pub count_n10: T,
    pub count_n11: T,
    pub var_tail_probability_p: T,
    pub portfolio_mean_loss: T,
    pub portfolio_std_loss: T,
    pub loss_skewness_gamma1: T,
    pub loss_excess_kurtosis_gamma2: T,
}
#[derive(Clone, Debug)]
pub struct ChristoffersenBouchaudVarBacktestAndCumulantReport<T: Real = f64> {
    pub empirical_hit_rate: T,
    pub transition_prob_pi01: T,
    pub transition_prob_pi11: T,
    pub kupiec_unconditional_lr_uc: T,
    pub christoffersen_independence_lr_ind: T,
    pub christoffersen_conditional_coverage_lr_cc: T,
    pub lr_decomposition_identity_error: T,
    pub gaussian_var: T,
    pub cornish_fisher_var: T,
    pub gram_charlier_expected_shortfall: T,
}
pub fn evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es<T: Real>(
    s: &ChristoffersenBouchaudVarBacktestAndCumulantSpec<T>,
) -> ChristoffersenBouchaudVarBacktestAndCumulantReport<T> {
    let n = s.count_n01 + s.count_n11;
    let t = T::from_f64(s.total_days_t.max(1) as f64);
    let hit = n / t;
    let pi01 = s.count_n01 / (s.count_n00 + s.count_n01).max_val(T::from_f64(1e-12));
    let pi11 = s.count_n11 / (s.count_n10 + s.count_n11).max_val(T::from_f64(1e-12));
    let p = s
        .var_tail_probability_p
        .max_val(T::from_f64(1e-12))
        .min_val(T::one() - T::from_f64(1e-12));
    let loglik = |q: T| {
        s.count_n01 * q.max_val(T::from_f64(1e-12)).ln()
            + s.count_n00 * (T::one() - q).max_val(T::from_f64(1e-12)).ln()
            + s.count_n11 * q.max_val(T::from_f64(1e-12)).ln()
            + s.count_n10 * (T::one() - q).max_val(T::from_f64(1e-12)).ln()
    };
    let lr_uc = T::from_f64(2.0)
        * (loglik(hit) - (n * p.ln() + (t - n) * (T::one() - p).max_val(T::from_f64(1e-12)).ln()));
    let lr_ind = T::from_f64(2.0)
        * (s.count_n01 * pi01.max_val(T::from_f64(1e-12)).ln()
            + s.count_n00 * (T::one() - pi01).max_val(T::from_f64(1e-12)).ln()
            + s.count_n11 * pi11.max_val(T::from_f64(1e-12)).ln()
            + s.count_n10 * (T::one() - pi11).max_val(T::from_f64(1e-12)).ln()
            - loglik(hit));
    let cc = lr_uc + lr_ind;
    let z = (-p.inv_norm_cdf()).max_val(T::from_f64(1e-6));
    let gaussian = s.portfolio_mean_loss + s.portfolio_std_loss * z;
    let cf = gaussian
        + s.portfolio_std_loss
            * (T::from_f64(1.0 / 6.0) * (z * z - T::one()) * s.loss_skewness_gamma1
                + T::from_f64(1.0 / 24.0)
                    * (z * z * z - T::from_f64(3.0) * z)
                    * s.loss_excess_kurtosis_gamma2
                - T::from_f64(1.0 / 36.0)
                    * (T::from_f64(2.0) * z * z * z - T::from_f64(5.0) * z)
                    * s.loss_skewness_gamma1
                    * s.loss_skewness_gamma1);
    let es = s.portfolio_mean_loss
        + s.portfolio_std_loss
            * (z.norm_pdf() / p)
            * (T::one()
                + T::from_f64(1.0 / 6.0) * (z * z - T::one()) * s.loss_skewness_gamma1
                + T::from_f64(1.0 / 24.0)
                    * (z * z * z - T::from_f64(3.0) * z)
                    * s.loss_excess_kurtosis_gamma2);
    ChristoffersenBouchaudVarBacktestAndCumulantReport {
        empirical_hit_rate: hit,
        transition_prob_pi01: pi01,
        transition_prob_pi11: pi11,
        kupiec_unconditional_lr_uc: lr_uc,
        christoffersen_independence_lr_ind: lr_ind,
        christoffersen_conditional_coverage_lr_cc: cc,
        lr_decomposition_identity_error: cc - lr_uc - lr_ind,
        gaussian_var: gaussian,
        cornish_fisher_var: cf,
        gram_charlier_expected_shortfall: es,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FoucaultVanDerMerweLobLiquidityLvarSpec<T: Real = f64> {
    pub mid_price: T,
    pub fundamental_volatility: T,
    pub informed_trader_share_pi: T,
    pub exponential_order_size_rate_lambda: T,
    pub order_processing_fee_c: T,
    pub position_notional: T,
    pub mean_relative_bid_ask_spread: T,
    pub std_relative_bid_ask_spread: T,
    pub spread_tail_scaling_a: T,
    pub confidence_z_alpha: T,
    pub liquidation_horizon_days: T,
    pub market_adv_notional: T,
    pub permanent_impact_gamma: T,
}
#[derive(Clone, Copy, Debug)]
pub struct FoucaultVanDerMerweLobLiquidityLvarReport<T: Real = f64> {
    pub glosten_sandas_small_trade_half_spread: T,
    pub glosten_sandas_marginal_slope: T,
    pub equilibrium_quoted_ask_depth_q: T,
    pub pure_market_risk_var: T,
    pub exogenous_spread_liquidity_cost: T,
    pub endogenous_market_impact_cost: T,
    pub total_liquidity_adjusted_var_lvar: T,
    pub liquidity_share_of_total_lvar: T,
}
pub fn evaluate_foucault_vandermerwe_lob_schedule_and_lvar<T: Real>(
    s: &FoucaultVanDerMerweLobLiquidityLvarSpec<T>,
) -> FoucaultVanDerMerweLobLiquidityLvarReport<T> {
    let pi_ratio = s.informed_trader_share_pi
        / (T::one() - s.informed_trader_share_pi).max_val(T::from_f64(1e-6));
    let lam = s
        .exponential_order_size_rate_lambda
        .abs()
        .max_val(T::from_f64(1e-12));
    let mid = s.mid_price.abs().max_val(T::from_f64(1e-12));
    let half = s.mean_relative_bid_ask_spread * mid * T::from_f64(0.5)
        + s.order_processing_fee_c * mid
        + pi_ratio * s.fundamental_volatility * mid / lam * T::from_f64(0.1);
    let slope = pi_ratio * s.fundamental_volatility * mid / lam;
    let depth = (T::one() - s.informed_trader_share_pi).max_val(T::from_f64(1e-6)) / lam;
    let pure = s.position_notional
        * s.fundamental_volatility
        * s.confidence_z_alpha
        * (s.liquidation_horizon_days / T::from_f64(252.0)).sqrt();
    let exo = s.position_notional
        * (T::from_f64(0.5)
            * (s.mean_relative_bid_ask_spread
                + s.spread_tail_scaling_a * s.std_relative_bid_ask_spread)
            + s.order_processing_fee_c);
    let impact = s.permanent_impact_gamma
        * s.position_notional.abs()
        * (s.position_notional.abs()
            / (s.market_adv_notional.abs() * s.liquidation_horizon_days.max_val(T::one()))
                .max_val(T::from_f64(1e-12)))
        .sqrt()
        + s.position_notional.abs() * (slope * depth / mid) * T::from_f64(0.01);
    let total = pure + exo + impact;
    FoucaultVanDerMerweLobLiquidityLvarReport {
        glosten_sandas_small_trade_half_spread: half,
        glosten_sandas_marginal_slope: slope,
        equilibrium_quoted_ask_depth_q: depth,
        pure_market_risk_var: pure,
        exogenous_spread_liquidity_cost: exo,
        endogenous_market_impact_cost: impact,
        total_liquidity_adjusted_var_lvar: total,
        liquidity_share_of_total_lvar: (exo + impact) / total.abs().max_val(T::from_f64(1e-12)),
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BarberisJinWangProspectAnomalySpec<T: Real = f64> {
    pub current_price: T,
    pub reference_purchase_price: T,
    pub return_volatility: T,
    pub return_skewness: T,
    pub prospect_loss_aversion_lambda: T,
    pub prospect_diminishing_sensitivity_alpha: T,
    pub prelec_probability_weighting_gamma: T,
    pub mental_accounting_integration_b0: T,
    pub capm_beta: T,
    pub equity_risk_premium: T,
}
#[derive(Clone, Copy, Debug)]
pub struct BarberisJinWangProspectAnomalyReport<T: Real = f64> {
    pub capital_gains_overhang_cgo: T,
    pub prior_gain_cushion_utility: T,
    pub prior_loss_pain_amplification: T,
    pub lottery_skewness_overvaluation_drag: T,
    pub disposition_effect_underreaction_alpha: T,
    pub bjw_equilibrium_expected_excess_return: T,
    pub bjw_anomaly_alpha_vs_capm: T,
}
pub fn evaluate_barberis_jin_wang_prospect_overhang_and_skew_anomaly<T: Real>(
    s: &BarberisJinWangProspectAnomalySpec<T>,
) -> BarberisJinWangProspectAnomalyReport<T> {
    let cgo = (s.current_price - s.reference_purchase_price)
        / s.reference_purchase_price.abs().max_val(T::from_f64(1e-12));
    let alpha = s
        .prospect_diminishing_sensitivity_alpha
        .max_val(T::from_f64(0.01));
    let eps = T::from_f64(1e-6);
    let pos_cgo = (cgo + (cgo * cgo + eps * eps).sqrt()) * T::from_f64(0.5);
    let neg_cgo = (-cgo + (cgo * cgo + eps * eps).sqrt()) * T::from_f64(0.5);
    let gain = (pos_cgo + eps).pow(alpha) - eps.pow(alpha);
    let pain = s.prospect_loss_aversion_lambda * ((neg_cgo + eps).pow(alpha) - eps.pow(alpha));
    let drag = -(s.return_skewness.max_val(T::zero())
        * (T::one() - s.prelec_probability_weighting_gamma).max_val(T::zero())
        * s.return_volatility);
    let under = s.mental_accounting_integration_b0 * (cgo + T::from_f64(0.25) * (gain - pain))
        / (T::one() + s.return_volatility.abs());
    let expected = s.capm_beta * s.equity_risk_premium + under + drag;
    BarberisJinWangProspectAnomalyReport {
        capital_gains_overhang_cgo: cgo,
        prior_gain_cushion_utility: gain,
        prior_loss_pain_amplification: pain,
        lottery_skewness_overvaluation_drag: drag,
        disposition_effect_underreaction_alpha: under,
        bjw_equilibrium_expected_excess_return: expected,
        bjw_anomaly_alpha_vs_capm: expected - s.capm_beta * s.equity_risk_premium,
    }
}

#[derive(Clone, Debug)]
pub struct HartMasColellRegretMatchingPortfolioSpec<T: Real = f64> {
    pub strategy_payoff_history: Vec<[T; 3]>,
    pub realized_action_weights_history: Vec<[T; 3]>,
    pub regret_inertia_mu: T,
    pub exploration_epsilon: T,
}
#[derive(Clone, Debug)]
pub struct HartMasColellRegretMatchingPortfolioReport<T: Real = f64> {
    pub cumulative_counterfactual_regrets: [T; 3],
    pub positive_regret_potentials: [T; 3],
    pub next_period_regret_matching_weights: [T; 3],
    pub weight_simplex_sum_error: T,
    pub blackwell_approachability_distance: T,
    pub hannan_external_regret_per_period: T,
    pub theoretical_no_regret_upper_bound: T,
}
pub fn evaluate_hart_mascolell_regret_matching_allocation<T: Real>(
    s: &HartMasColellRegretMatchingPortfolioSpec<T>,
) -> HartMasColellRegretMatchingPortfolioReport<T> {
    let mut reg = [T::zero(); 3];
    for (pay, w) in s
        .strategy_payoff_history
        .iter()
        .zip(s.realized_action_weights_history.iter())
    {
        let realized = pay[0] * w[0] + pay[1] * w[1] + pay[2] * w[2];
        for k in 0..3 {
            reg[k] += s.regret_inertia_mu * (pay[k] - realized);
        }
    }
    let pos = [
        reg[0].max_val(T::zero()),
        reg[1].max_val(T::zero()),
        reg[2].max_val(T::zero()),
    ];
    let denom = pos[0] + pos[1] + pos[2] + s.exploration_epsilon * T::from_f64(3.0);
    let w = [
        (pos[0] + s.exploration_epsilon) / denom,
        (pos[1] + s.exploration_epsilon) / denom,
        (pos[2] + s.exploration_epsilon) / denom,
    ];
    let best = pos.iter().copied().fold(T::zero(), |a, b| a.max_val(b));
    let mut payoff_min = T::from_f64(f64::INFINITY);
    let mut payoff_max = T::from_f64(f64::NEG_INFINITY);
    for pay in &s.strategy_payoff_history {
        for value in pay {
            payoff_min = payoff_min.min_val(*value);
            payoff_max = payoff_max.max_val(*value);
        }
    }
    let n_periods = T::from_f64(s.strategy_payoff_history.len().max(1) as f64);
    let bound =
        (payoff_max - payoff_min) * (T::from_f64(2.0) * T::from_f64(3.0).ln() / n_periods).sqrt();
    HartMasColellRegretMatchingPortfolioReport {
        cumulative_counterfactual_regrets: reg,
        positive_regret_potentials: pos,
        next_period_regret_matching_weights: w,
        weight_simplex_sum_error: w[0] + w[1] + w[2] - T::one(),
        blackwell_approachability_distance: (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2])
            .sqrt()
            / n_periods,
        hannan_external_regret_per_period: best
            / T::from_f64(s.strategy_payoff_history.len().max(1) as f64),
        theoretical_no_regret_upper_bound: bound,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ColinBaconCampisiFixedIncomeAttributionSpec<T: Real = f64> {
    pub portfolio_total_return: T,
    pub benchmark_total_return: T,
    pub portfolio_coupon_yield: T,
    pub benchmark_coupon_yield: T,
    pub portfolio_eff_duration: T,
    pub benchmark_eff_duration: T,
    pub portfolio_spread_duration: T,
    pub benchmark_spread_duration: T,
    pub parallel_treasury_shift: T,
    pub portfolio_key_rate_twist_shift: T,
    pub benchmark_key_rate_twist_shift: T,
    pub portfolio_credit_spread_change: T,
    pub benchmark_credit_spread_change: T,
    pub fx_spot_return: T,
    pub domestic_cash_rate: T,
    pub foreign_cash_rate: T,
}
#[derive(Clone, Copy, Debug)]
pub struct ColinBaconCampisiFixedIncomeAttributionReport<T: Real = f64> {
    pub active_income_effect: T,
    pub active_treasury_shift_effect: T,
    pub active_treasury_twist_effect: T,
    pub active_spread_allocation_effect: T,
    pub active_security_selection_residual: T,
    pub total_active_local_return: T,
    pub campisi_identity_conservation_error: T,
    pub karnosky_singer_forward_hedge_carry: T,
    pub karnosky_singer_unhedged_currency_surprise: T,
}
pub fn evaluate_colin_bacon_campisi_fixed_income_attribution<T: Real>(
    s: &ColinBaconCampisiFixedIncomeAttributionSpec<T>,
) -> ColinBaconCampisiFixedIncomeAttributionReport<T> {
    let active = s.portfolio_total_return - s.benchmark_total_return;
    let income = s.portfolio_coupon_yield - s.benchmark_coupon_yield;
    let treasury =
        -(s.portfolio_eff_duration - s.benchmark_eff_duration) * s.parallel_treasury_shift;
    let twist = -(s.portfolio_eff_duration * s.portfolio_key_rate_twist_shift
        - s.benchmark_eff_duration * s.benchmark_key_rate_twist_shift);
    let spread = -(s.portfolio_spread_duration * s.portfolio_credit_spread_change
        - s.benchmark_spread_duration * s.benchmark_credit_spread_change);
    let residual = active - income - treasury - twist - spread;
    let carry = s.domestic_cash_rate - s.foreign_cash_rate;
    let currency = s.fx_spot_return - carry;
    ColinBaconCampisiFixedIncomeAttributionReport {
        active_income_effect: income,
        active_treasury_shift_effect: treasury,
        active_treasury_twist_effect: twist,
        active_spread_allocation_effect: spread,
        active_security_selection_residual: residual,
        total_active_local_return: income + treasury + twist + spread + residual,
        campisi_identity_conservation_error: income + treasury + twist + spread + residual - active,
        karnosky_singer_forward_hedge_carry: carry,
        karnosky_singer_unhedged_currency_surprise: currency,
    }
}

#[cfg(test)]
mod phase30_tests {
    use super::*;
    use ql_ad::Dual;
    fn var_spec<T: Real>(x: T) -> ChristoffersenBouchaudVarBacktestAndCumulantSpec<T> {
        ChristoffersenBouchaudVarBacktestAndCumulantSpec {
            total_days_t: 100,
            count_n00: T::from_f64(80.0),
            count_n01: T::from_f64(8.0),
            count_n10: T::from_f64(8.0),
            count_n11: T::from_f64(4.0),
            var_tail_probability_p: T::from_f64(0.05),
            portfolio_mean_loss: x,
            portfolio_std_loss: T::one(),
            loss_skewness_gamma1: T::zero(),
            loss_excess_kurtosis_gamma2: T::zero(),
        }
    }
    #[test]
    fn phase30_portfolio_invariants_and_dual_fd() {
        let r = evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es(&var_spec(0.1));
        assert!(r.lr_decomposition_identity_error.abs() < 1e-12);
        assert!((r.cornish_fisher_var - r.gaussian_var).abs() < 1e-12);
        let h = 1e-6;
        let fd =
            (evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es(&var_spec(0.1 + h))
                .gaussian_var
                - evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es(&var_spec(
                    0.1 - h,
                ))
                .gaussian_var)
                / (2.0 * h);
        let ad = evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es(&var_spec(
            Dual::<1>::variable(0.1, 0),
        ))
        .gaussian_var;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        let mut vup = var_spec(0.1);
        let mut vdn = var_spec(0.1);
        vup.portfolio_mean_loss += h;
        vdn.portfolio_mean_loss -= h;
        let fd_var = (evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es(&vup)
            .cornish_fisher_var
            - evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es(&vdn)
                .cornish_fisher_var)
            / (2.0 * h);
        let ad_var = evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es(&var_spec(
            Dual::<1>::variable(0.1, 0),
        ))
        .cornish_fisher_var;
        assert!(ad_var.d[0].abs() > 1e-6 && (ad_var.d[0] - fd_var).abs() < 1e-5);
        let l = evaluate_foucault_vandermerwe_lob_schedule_and_lvar(
            &FoucaultVanDerMerweLobLiquidityLvarSpec {
                mid_price: 100.0,
                fundamental_volatility: 0.2,
                informed_trader_share_pi: 0.2,
                exponential_order_size_rate_lambda: 1.0,
                order_processing_fee_c: 0.01,
                position_notional: 100.0,
                mean_relative_bid_ask_spread: 0.01,
                std_relative_bid_ask_spread: 0.002,
                spread_tail_scaling_a: 1.0,
                confidence_z_alpha: 2.0,
                liquidation_horizon_days: 1.0,
                market_adv_notional: 1000.0,
                permanent_impact_gamma: 0.1,
            },
        );
        assert!(
            (l.total_liquidity_adjusted_var_lvar
                - l.pure_market_risk_var
                - l.exogenous_spread_liquidity_cost
                - l.endogenous_market_impact_cost)
                .abs()
                < 1e-12
        );
        let mut lup = FoucaultVanDerMerweLobLiquidityLvarSpec {
            mid_price: 100.0,
            fundamental_volatility: 0.2,
            informed_trader_share_pi: 0.2,
            exponential_order_size_rate_lambda: 1.0,
            order_processing_fee_c: 0.01,
            position_notional: 100.0,
            mean_relative_bid_ask_spread: 0.01,
            std_relative_bid_ask_spread: 0.002,
            spread_tail_scaling_a: 1.0,
            confidence_z_alpha: 2.0,
            liquidation_horizon_days: 1.0,
            market_adv_notional: 1000.0,
            permanent_impact_gamma: 0.1,
        };
        let mut ldn = lup;
        lup.order_processing_fee_c += h;
        ldn.order_processing_fee_c -= h;
        let fd_l = (evaluate_foucault_vandermerwe_lob_schedule_and_lvar(&lup)
            .total_liquidity_adjusted_var_lvar
            - evaluate_foucault_vandermerwe_lob_schedule_and_lvar(&ldn)
                .total_liquidity_adjusted_var_lvar)
            / (2.0 * h);
        let lad: FoucaultVanDerMerweLobLiquidityLvarSpec<Dual<1>> =
            FoucaultVanDerMerweLobLiquidityLvarSpec {
                mid_price: Dual::from_f64(100.0),
                fundamental_volatility: Dual::from_f64(0.2),
                informed_trader_share_pi: Dual::from_f64(0.2),
                exponential_order_size_rate_lambda: Dual::from_f64(1.0),
                order_processing_fee_c: Dual::<1>::variable(0.01, 0),
                position_notional: Dual::from_f64(100.0),
                mean_relative_bid_ask_spread: Dual::from_f64(0.01),
                std_relative_bid_ask_spread: Dual::from_f64(0.002),
                spread_tail_scaling_a: Dual::from_f64(1.0),
                confidence_z_alpha: Dual::from_f64(2.0),
                liquidation_horizon_days: Dual::from_f64(1.0),
                market_adv_notional: Dual::from_f64(1000.0),
                permanent_impact_gamma: Dual::from_f64(0.1),
            };
        let ad_l = evaluate_foucault_vandermerwe_lob_schedule_and_lvar(&lad)
            .total_liquidity_adjusted_var_lvar;
        assert!(ad_l.d[0].abs() > 1e-6 && (ad_l.d[0] - fd_l).abs() < 1e-5);
        let p = evaluate_barberis_jin_wang_prospect_overhang_and_skew_anomaly(
            &BarberisJinWangProspectAnomalySpec {
                current_price: 110.0,
                reference_purchase_price: 100.0,
                return_volatility: 0.2,
                return_skewness: 0.8,
                prospect_loss_aversion_lambda: 2.0,
                prospect_diminishing_sensitivity_alpha: 0.8,
                prelec_probability_weighting_gamma: 0.7,
                mental_accounting_integration_b0: 0.4,
                capm_beta: 1.0,
                equity_risk_premium: 0.05,
            },
        );
        assert!(p.capital_gains_overhang_cgo > 0.0 && p.lottery_skewness_overvaluation_drag < 0.0);
        let mut pup = BarberisJinWangProspectAnomalySpec {
            current_price: 110.0,
            reference_purchase_price: 100.0,
            return_volatility: 0.2,
            return_skewness: 0.8,
            prospect_loss_aversion_lambda: 2.0,
            prospect_diminishing_sensitivity_alpha: 0.8,
            prelec_probability_weighting_gamma: 0.7,
            mental_accounting_integration_b0: 0.4,
            capm_beta: 1.0,
            equity_risk_premium: 0.05,
        };
        let mut pdn = pup;
        pup.current_price += h;
        pdn.current_price -= h;
        let fd_p = (evaluate_barberis_jin_wang_prospect_overhang_and_skew_anomaly(&pup)
            .bjw_equilibrium_expected_excess_return
            - evaluate_barberis_jin_wang_prospect_overhang_and_skew_anomaly(&pdn)
                .bjw_equilibrium_expected_excess_return)
            / (2.0 * h);
        let pad: BarberisJinWangProspectAnomalySpec<Dual<1>> = BarberisJinWangProspectAnomalySpec {
            current_price: Dual::<1>::variable(110.0, 0),
            reference_purchase_price: Dual::from_f64(100.0),
            return_volatility: Dual::from_f64(0.2),
            return_skewness: Dual::from_f64(0.8),
            prospect_loss_aversion_lambda: Dual::from_f64(2.0),
            prospect_diminishing_sensitivity_alpha: Dual::from_f64(0.8),
            prelec_probability_weighting_gamma: Dual::from_f64(0.7),
            mental_accounting_integration_b0: Dual::from_f64(0.4),
            capm_beta: Dual::from_f64(1.0),
            equity_risk_premium: Dual::from_f64(0.05),
        };
        let ad_p = evaluate_barberis_jin_wang_prospect_overhang_and_skew_anomaly(&pad)
            .bjw_equilibrium_expected_excess_return;
        assert!(ad_p.d[0].abs() > 1e-6 && (ad_p.d[0] - fd_p).abs() < 1e-5);
        let hrt = evaluate_hart_mascolell_regret_matching_allocation(
            &HartMasColellRegretMatchingPortfolioSpec {
                strategy_payoff_history: vec![[0.1, 0.2, 0.0], [0.4, 0.1, 0.3]],
                realized_action_weights_history: vec![[0.3, 0.4, 0.3], [0.4, 0.3, 0.3]],
                regret_inertia_mu: 1.0,
                exploration_epsilon: 0.01,
            },
        );
        assert!(
            hrt.next_period_regret_matching_weights
                .iter()
                .all(|x| *x >= 0.0)
                && hrt.weight_simplex_sum_error.abs() < 1e-12
        );
        let mut rup = HartMasColellRegretMatchingPortfolioSpec {
            strategy_payoff_history: vec![[0.1, 0.2, 0.0], [0.4, 0.1, 0.3]],
            realized_action_weights_history: vec![[0.3, 0.4, 0.3], [0.4, 0.3, 0.3]],
            regret_inertia_mu: 1.0,
            exploration_epsilon: 0.01,
        };
        let mut rdn = rup.clone();
        rup.regret_inertia_mu += h;
        rdn.regret_inertia_mu -= h;
        let fd_r = (evaluate_hart_mascolell_regret_matching_allocation(&rup)
            .blackwell_approachability_distance
            - evaluate_hart_mascolell_regret_matching_allocation(&rdn)
                .blackwell_approachability_distance)
            / (2.0 * h);
        let rad: HartMasColellRegretMatchingPortfolioSpec<Dual<1>> =
            HartMasColellRegretMatchingPortfolioSpec {
                strategy_payoff_history: vec![
                    [
                        Dual::from_f64(0.1),
                        Dual::from_f64(0.2),
                        Dual::from_f64(0.0),
                    ],
                    [
                        Dual::from_f64(0.4),
                        Dual::from_f64(0.1),
                        Dual::from_f64(0.3),
                    ],
                ],
                realized_action_weights_history: vec![
                    [
                        Dual::from_f64(0.3),
                        Dual::from_f64(0.4),
                        Dual::from_f64(0.3),
                    ],
                    [
                        Dual::from_f64(0.4),
                        Dual::from_f64(0.3),
                        Dual::from_f64(0.3),
                    ],
                ],
                regret_inertia_mu: Dual::<1>::variable(1.0, 0),
                exploration_epsilon: Dual::from_f64(0.01),
            };
        let ad_r = evaluate_hart_mascolell_regret_matching_allocation(&rad)
            .blackwell_approachability_distance;
        assert!(ad_r.d[0].abs() > 1e-6 && (ad_r.d[0] - fd_r).abs() < 1e-5);
        let c = evaluate_colin_bacon_campisi_fixed_income_attribution(
            &ColinBaconCampisiFixedIncomeAttributionSpec {
                portfolio_total_return: 0.04,
                benchmark_total_return: 0.03,
                portfolio_coupon_yield: 0.04,
                benchmark_coupon_yield: 0.03,
                portfolio_eff_duration: 5.0,
                benchmark_eff_duration: 4.0,
                portfolio_spread_duration: 3.0,
                benchmark_spread_duration: 3.0,
                parallel_treasury_shift: 0.01,
                portfolio_key_rate_twist_shift: 0.001,
                benchmark_key_rate_twist_shift: 0.002,
                portfolio_credit_spread_change: 0.002,
                benchmark_credit_spread_change: 0.001,
                fx_spot_return: 0.01,
                domestic_cash_rate: 0.02,
                foreign_cash_rate: 0.03,
            },
        );
        assert!(c.campisi_identity_conservation_error.abs() < 1e-12);
        let mut cup = ColinBaconCampisiFixedIncomeAttributionSpec {
            portfolio_total_return: 0.04,
            benchmark_total_return: 0.03,
            portfolio_coupon_yield: 0.04,
            benchmark_coupon_yield: 0.03,
            portfolio_eff_duration: 5.0,
            benchmark_eff_duration: 4.0,
            portfolio_spread_duration: 3.0,
            benchmark_spread_duration: 3.0,
            parallel_treasury_shift: 0.01,
            portfolio_key_rate_twist_shift: 0.001,
            benchmark_key_rate_twist_shift: 0.002,
            portfolio_credit_spread_change: 0.002,
            benchmark_credit_spread_change: 0.001,
            fx_spot_return: 0.01,
            domestic_cash_rate: 0.02,
            foreign_cash_rate: 0.03,
        };
        let mut cdn = cup;
        cup.parallel_treasury_shift += h;
        cdn.parallel_treasury_shift -= h;
        let fd_c = (evaluate_colin_bacon_campisi_fixed_income_attribution(&cup)
            .active_security_selection_residual
            - evaluate_colin_bacon_campisi_fixed_income_attribution(&cdn)
                .active_security_selection_residual)
            / (2.0 * h);
        let cad: ColinBaconCampisiFixedIncomeAttributionSpec<Dual<1>> =
            ColinBaconCampisiFixedIncomeAttributionSpec {
                portfolio_total_return: Dual::from_f64(0.04),
                benchmark_total_return: Dual::from_f64(0.03),
                portfolio_coupon_yield: Dual::from_f64(0.04),
                benchmark_coupon_yield: Dual::from_f64(0.03),
                portfolio_eff_duration: Dual::from_f64(5.0),
                benchmark_eff_duration: Dual::from_f64(4.0),
                portfolio_spread_duration: Dual::from_f64(3.0),
                benchmark_spread_duration: Dual::from_f64(3.0),
                parallel_treasury_shift: Dual::<1>::variable(0.01, 0),
                portfolio_key_rate_twist_shift: Dual::from_f64(0.001),
                benchmark_key_rate_twist_shift: Dual::from_f64(0.002),
                portfolio_credit_spread_change: Dual::from_f64(0.002),
                benchmark_credit_spread_change: Dual::from_f64(0.001),
                fx_spot_return: Dual::from_f64(0.01),
                domestic_cash_rate: Dual::from_f64(0.02),
                foreign_cash_rate: Dual::from_f64(0.03),
            };
        let ad_c = evaluate_colin_bacon_campisi_fixed_income_attribution(&cad)
            .active_security_selection_residual;
        assert!(ad_c.d[0].abs() > 1e-6 && (ad_c.d[0] - fd_c).abs() < 1e-5);
    }
}
