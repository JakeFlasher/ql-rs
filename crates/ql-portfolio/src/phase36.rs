#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct CarteaJaimungalLobMarketMakingSpec<T: Real = f64> {
    pub mid_price: T,
    pub current_inventory_q: T,
    pub running_inventory_penalty_phi: T,
    pub terminal_inventory_penalty_alpha: T,
    pub risk_aversion_gamma: T,
    pub asset_volatility: T,
    pub order_arrival_baseline_a: T,
    pub lob_depth_decay_kappa: T,
    pub remaining_horizon_tau: T,
}

#[derive(Clone, Copy, Debug)]
pub struct CarteaJaimungalLobMarketMakingReport<T: Real = f64> {
    pub riccati_inventory_coeff_theta: T,
    pub optimal_bid_depth_delta_minus: T,
    pub optimal_ask_depth_delta_plus: T,
    pub optimal_bid_quote: T,
    pub optimal_ask_quote: T,
    pub reservation_mid_price: T,
    pub inventory_quote_skew: T,
    pub total_quoted_spread: T,
    pub expected_instantaneous_spread_capture: T,
}

pub fn evaluate_cartea_jaimungal_lob_market_making<T: Real>(
    spec: &CarteaJaimungalLobMarketMakingSpec<T>,
) -> CarteaJaimungalLobMarketMakingReport<T> {
    let kappa = spec.lob_depth_decay_kappa.max_val(T::from_f64(1e-6));
    let tau = spec.remaining_horizon_tau.max_val(T::from_f64(1e-6));
    let inv_e = (-T::one()).exp();
    let theta = T::from_f64(0.5)
        * (spec.running_inventory_penalty_phi
            + spec.risk_aversion_gamma * spec.asset_volatility * spec.asset_volatility)
        * tau
        + spec.terminal_inventory_penalty_alpha
            / (T::one()
                + kappa
                    * spec.terminal_inventory_penalty_alpha
                    * spec.order_arrival_baseline_a
                    * inv_e
                    * tau);
    let two = T::from_f64(2.0);
    let q = spec.current_inventory_q;
    let delta_plus = T::one() / kappa + (two * q - T::one()) * theta;
    let delta_minus = T::one() / kappa - (two * q + T::one()) * theta;
    let ask_quote = spec.mid_price + delta_plus;
    let bid_quote = spec.mid_price - delta_minus;
    let res_mid = T::from_f64(0.5) * (ask_quote + bid_quote);
    let skew = -two * q * theta;
    let spread = ask_quote - bid_quote;
    let lam_plus = spec.order_arrival_baseline_a * (-kappa * delta_plus).exp();
    let lam_minus = spec.order_arrival_baseline_a * (-kappa * delta_minus).exp();
    let edge = lam_plus * delta_plus + lam_minus * delta_minus;

    CarteaJaimungalLobMarketMakingReport {
        riccati_inventory_coeff_theta: theta,
        optimal_bid_depth_delta_minus: delta_minus,
        optimal_ask_depth_delta_plus: delta_plus,
        optimal_bid_quote: bid_quote,
        optimal_ask_quote: ask_quote,
        reservation_mid_price: res_mid,
        inventory_quote_skew: skew,
        total_quoted_spread: spread,
        expected_instantaneous_spread_capture: edge,
    }
}

#[derive(Clone, Debug)]
pub struct CampbellLoMackinlayVrPresentValueSpec<T: Real = f64> {
    pub returns: Vec<T>,
    pub steady_state_dividend_yield: T,
    pub expected_dividend_growth: T,
    pub expected_discount_rate: T,
    pub dividend_growth_ar1_phi: T,
    pub discount_rate_ar1_phi: T,
}

#[derive(Clone, Copy, Debug)]
pub struct CampbellLoMackinlayVrPresentValueReport<T: Real = f64> {
    pub one_period_variance: T,
    pub overlapping_two_period_variance: T,
    pub lo_mackinlay_variance_ratio_vr2: T,
    pub homoskedastic_z_statistic: T,
    pub heteroskedastic_theta_variance: T,
    pub heteroskedastic_z_star_statistic: T,
    pub campbell_shiller_rho: T,
    pub equilibrium_log_dividend_price_ratio: T,
    pub discount_rate_news_multiplier: T,
    pub cashflow_news_multiplier: T,
}

pub fn evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller<T: Real>(
    spec: &CampbellLoMackinlayVrPresentValueSpec<T>,
) -> CampbellLoMackinlayVrPresentValueReport<T> {
    let n = spec.returns.len().max(3);
    let n_t = T::from_f64(n as f64);
    let mut sum_r = T::zero();
    for &r in &spec.returns {
        sum_r = sum_r + r;
    }
    let mu = sum_r / n_t;
    let mut sum_sq1 = T::zero();
    for &r in &spec.returns {
        let d = r - mu;
        sum_sq1 = sum_sq1 + d * d;
    }
    let var1 = (sum_sq1 / T::from_f64((n - 1) as f64)).max_val(T::from_f64(1e-12));

    let mut sum_sq2 = T::zero();
    let mut sum_theta = T::zero();
    let two = T::from_f64(2.0);
    for t in 1..spec.returns.len() {
        let d2 = spec.returns[t] + spec.returns[t - 1] - two * mu;
        sum_sq2 = sum_sq2 + d2 * d2;
        let u_t = (spec.returns[t] - mu) * (spec.returns[t] - mu) / var1;
        let u_tm1 = (spec.returns[t - 1] - mu) * (spec.returns[t - 1] - mu) / var1;
        sum_theta = sum_theta + u_t * u_tm1;
    }
    let var2 = sum_sq2 / T::from_f64((2 * (n - 2)) as f64);
    let vr2 = var2 / var1;
    let z_homo = n_t.sqrt() * (vr2 - T::one());
    let theta_var = (T::from_f64(4.0) / n_t * sum_theta).max_val(T::from_f64(1e-12));
    let z_star = z_homo / theta_var.sqrt();

    let dy = spec.steady_state_dividend_yield.max_val(T::from_f64(1e-6));
    let rho = T::one() / (T::one() + dy);
    let k_cs = -rho.ln() - (T::one() - rho) * (T::one() / rho - T::one()).ln();
    let cf_mult =
        T::one() / (T::one() - rho * spec.dividend_growth_ar1_phi).max_val(T::from_f64(1e-6));
    let dr_mult =
        T::one() / (T::one() - rho * spec.discount_rate_ar1_phi).max_val(T::from_f64(1e-6));
    let eq_dp = -k_cs / (T::one() - rho) - spec.expected_dividend_growth * cf_mult
        + spec.expected_discount_rate * dr_mult;

    CampbellLoMackinlayVrPresentValueReport {
        one_period_variance: var1,
        overlapping_two_period_variance: var2,
        lo_mackinlay_variance_ratio_vr2: vr2,
        homoskedastic_z_statistic: z_homo,
        heteroskedastic_theta_variance: theta_var,
        heteroskedastic_z_star_statistic: z_star,
        campbell_shiller_rho: rho,
        equilibrium_log_dividend_price_ratio: eq_dp,
        discount_rate_news_multiplier: dr_mult,
        cashflow_news_multiplier: cf_mult,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RoncalliSpinuErcRiskParitySpec<T: Real = f64> {
    pub covariance_3x3: [[T; 3]; 3],
    pub target_risk_budgets: [T; 3],
}

#[derive(Clone, Copy, Debug)]
pub struct RoncalliSpinuErcRiskParityReport<T: Real = f64> {
    pub erc_weights: [T; 3],
    pub inverse_vol_naive_weights: [T; 3],
    pub relative_risk_contributions: [T; 3],
    pub erc_portfolio_volatility: T,
    pub risk_budget_parity_l1_error: T,
    pub weights_sum_conservation_error: T,
    pub choueifaty_diversification_ratio: T,
    pub spinu_convex_objective_val: T,
}

pub fn evaluate_roncalli_spinu_erc_risk_parity<T: Real>(
    spec: &RoncalliSpinuErcRiskParitySpec<T>,
) -> RoncalliSpinuErcRiskParityReport<T> {
    let cov = &spec.covariance_3x3;
    let b_sum =
        (spec.target_risk_budgets[0] + spec.target_risk_budgets[1] + spec.target_risk_budgets[2])
            .max_val(T::from_f64(1e-8));
    let b = [
        spec.target_risk_budgets[0] / b_sum,
        spec.target_risk_budgets[1] / b_sum,
        spec.target_risk_budgets[2] / b_sum,
    ];
    let vols = [
        cov[0][0].max_val(T::from_f64(1e-10)).sqrt(),
        cov[1][1].max_val(T::from_f64(1e-10)).sqrt(),
        cov[2][2].max_val(T::from_f64(1e-10)).sqrt(),
    ];
    let mut x = [b[0] / vols[0], b[1] / vols[1], b[2] / vols[2]];
    let inv_sum = x[0] + x[1] + x[2];
    let inv_vol_w = [x[0] / inv_sum, x[1] / inv_sum, x[2] / inv_sum];

    let two = T::from_f64(2.0);
    let four = T::from_f64(4.0);
    for _ in 0..16 {
        for i in 0..3 {
            let mut c_i = T::zero();
            for j in 0..3 {
                if j != i {
                    c_i = c_i + cov[i][j] * x[j];
                }
            }
            let sig_ii = cov[i][i].max_val(T::from_f64(1e-10));
            x[i] = (-c_i + (c_i * c_i + four * sig_ii * b[i]).sqrt()) / (two * sig_ii);
        }
    }

    let mut x_sig_x = T::zero();
    for i in 0..3 {
        for j in 0..3 {
            x_sig_x = x_sig_x + x[i] * cov[i][j] * x[j];
        }
    }
    let spinu_obj =
        T::from_f64(0.5) * x_sig_x - (b[0] * x[0].ln() + b[1] * x[1].ln() + b[2] * x[2].ln());

    let x_tot = x[0] + x[1] + x[2];
    let w = [x[0] / x_tot, x[1] / x_tot, x[2] / x_tot];
    let mut sig_w = [T::zero(); 3];
    for i in 0..3 {
        sig_w[i] = cov[i][0] * w[0] + cov[i][1] * w[1] + cov[i][2] * w[2];
    }
    let port_var =
        (w[0] * sig_w[0] + w[1] * sig_w[1] + w[2] * sig_w[2]).max_val(T::from_f64(1e-12));
    let port_vol = port_var.sqrt();
    let rrc = [
        w[0] * sig_w[0] / port_var,
        w[1] * sig_w[1] / port_var,
        w[2] * sig_w[2] / port_var,
    ];
    let mut l1_err = T::zero();
    for i in 0..3 {
        let d = rrc[i] - b[i];
        l1_err = l1_err + (d * d + T::from_f64(1e-30)).sqrt();
    }
    let w_diff = (w[0] + w[1] + w[2]) - T::one();
    let w_err = (w_diff * w_diff + T::from_f64(1e-30)).sqrt();
    let dr = (w[0] * vols[0] + w[1] * vols[1] + w[2] * vols[2]) / port_vol;

    RoncalliSpinuErcRiskParityReport {
        erc_weights: w,
        inverse_vol_naive_weights: inv_vol_w,
        relative_risk_contributions: rrc,
        erc_portfolio_volatility: port_vol,
        risk_budget_parity_l1_error: l1_err,
        weights_sum_conservation_error: w_err,
        choueifaty_diversification_ratio: dr,
        spinu_convex_objective_val: spinu_obj,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct IlmanenCarryRolldownVmcFactorSpec<T: Real = f64> {
    pub current_yield_yt: T,
    pub rolled_horizon_yield: T,
    pub repo_financing_rate: T,
    pub modified_duration: T,
    pub bond_convexity: T,
    pub expected_yield_change: T,
    pub yield_volatility: T,
    pub value_style_sharpe: T,
    pub momentum_style_sharpe: T,
    pub carry_style_sharpe: T,
    pub value_momentum_corr_rho: T,
}

#[derive(Clone, Copy, Debug)]
pub struct IlmanenCarryRolldownVmcFactorReport<T: Real = f64> {
    pub pure_carry_return: T,
    pub rolldown_return: T,
    pub expected_valuation_return: T,
    pub convexity_bias_return: T,
    pub total_expected_excess_return: T,
    pub carry_plus_rolldown_breakeven_bps: T,
    pub standalone_average_style_sharpe: T,
    pub composite_vmc_sharpe_ratio: T,
    pub style_diversification_sharpe_gain: T,
}

pub fn evaluate_ilmanen_carry_rolldown_and_vmc_factors<T: Real>(
    spec: &IlmanenCarryRolldownVmcFactorSpec<T>,
) -> IlmanenCarryRolldownVmcFactorReport<T> {
    let carry = spec.current_yield_yt - spec.repo_financing_rate;
    let rolldown = -spec.modified_duration * (spec.rolled_horizon_yield - spec.current_yield_yt);
    let val_ret = -spec.modified_duration * spec.expected_yield_change;
    let conv_bias = T::from_f64(0.5)
        * spec.bond_convexity
        * (spec.expected_yield_change * spec.expected_yield_change
            + spec.yield_volatility * spec.yield_volatility);
    let total_ret = carry + rolldown + val_ret + conv_bias;
    let breakeven_bps =
        (carry + rolldown) / spec.modified_duration.max_val(T::from_f64(1e-6)) * T::from_f64(1e4);

    let sum_sr = spec.value_style_sharpe + spec.momentum_style_sharpe + spec.carry_style_sharpe;
    let avg_sr = sum_sr / T::from_f64(3.0);
    let denom = (T::from_f64(3.0) + T::from_f64(2.0) * spec.value_momentum_corr_rho)
        .max_val(T::from_f64(0.2))
        .sqrt();
    let vmc_sr = sum_sr / denom;
    let sr_gain = vmc_sr - avg_sr;

    IlmanenCarryRolldownVmcFactorReport {
        pure_carry_return: carry,
        rolldown_return: rolldown,
        expected_valuation_return: val_ret,
        convexity_bias_return: conv_bias,
        total_expected_excess_return: total_ret,
        carry_plus_rolldown_breakeven_bps: breakeven_bps,
        standalone_average_style_sharpe: avg_sr,
        composite_vmc_sharpe_ratio: vmc_sr,
        style_diversification_sharpe_gain: sr_gain,
    }
}

#[cfg(test)]
mod phase36_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn phase36_portfolio_all_four_engines_invariants_and_ad_parity() {
        let h = 1e-6;

        // 1. Cartea-Jaimungal LOB Market Making
        let lob_spec = CarteaJaimungalLobMarketMakingSpec {
            mid_price: 100.0,
            current_inventory_q: 3.0,
            running_inventory_penalty_phi: 0.02,
            terminal_inventory_penalty_alpha: 0.05,
            risk_aversion_gamma: 0.10,
            asset_volatility: 0.20,
            order_arrival_baseline_a: 12.0,
            lob_depth_decay_kappa: 1.5,
            remaining_horizon_tau: 0.5,
        };
        let lob_rep = evaluate_cartea_jaimungal_lob_market_making(&lob_spec);
        assert!(lob_rep.inventory_quote_skew < 0.0);
        assert!(
            (lob_rep.optimal_ask_quote - lob_rep.optimal_bid_quote - lob_rep.total_quoted_spread)
                .abs()
                < 1e-12
        );
        let lob_ad =
            evaluate_cartea_jaimungal_lob_market_making(&CarteaJaimungalLobMarketMakingSpec {
                mid_price: Dual::<1>::from_value(lob_spec.mid_price),
                current_inventory_q: Dual::<1>::variable(lob_spec.current_inventory_q, 0),
                running_inventory_penalty_phi: Dual::<1>::from_value(
                    lob_spec.running_inventory_penalty_phi,
                ),
                terminal_inventory_penalty_alpha: Dual::<1>::from_value(
                    lob_spec.terminal_inventory_penalty_alpha,
                ),
                risk_aversion_gamma: Dual::<1>::from_value(lob_spec.risk_aversion_gamma),
                asset_volatility: Dual::<1>::from_value(lob_spec.asset_volatility),
                order_arrival_baseline_a: Dual::<1>::from_value(lob_spec.order_arrival_baseline_a),
                lob_depth_decay_kappa: Dual::<1>::from_value(lob_spec.lob_depth_decay_kappa),
                remaining_horizon_tau: Dual::<1>::from_value(lob_spec.remaining_horizon_tau),
            });
        let mut lp = lob_spec;
        lp.current_inventory_q += h;
        let mut lm = lob_spec;
        lm.current_inventory_q -= h;
        let _fd_lob = (evaluate_cartea_jaimungal_lob_market_making(&lp).reservation_mid_price
            - evaluate_cartea_jaimungal_lob_market_making(&lm).reservation_mid_price)
            / (2.0 * h);
        assert!((lob_ad.reservation_mid_price.d[0] - _fd_lob).abs() < 1e-5);

        // 2. Campbell-Lo-MacKinlay VR & Campbell-Shiller
        let clm_spec = CampbellLoMackinlayVrPresentValueSpec {
            returns: vec![0.012, 0.015, -0.004, 0.009, 0.014, -0.003, 0.011, 0.013],
            steady_state_dividend_yield: 0.032,
            expected_dividend_growth: 0.045,
            expected_discount_rate: 0.078,
            dividend_growth_ar1_phi: 0.35,
            discount_rate_ar1_phi: 0.85,
        };
        let clm_rep = evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller(&clm_spec);
        assert!(clm_rep.campbell_shiller_rho > 0.9 && clm_rep.campbell_shiller_rho < 1.0);
        assert!(clm_rep.lo_mackinlay_variance_ratio_vr2 > 0.0);
        assert!(clm_rep.discount_rate_news_multiplier > 1.0);
        let clm_ad = evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller(
            &CampbellLoMackinlayVrPresentValueSpec {
                returns: clm_spec
                    .returns
                    .iter()
                    .copied()
                    .map(Dual::<1>::from_value)
                    .collect(),
                steady_state_dividend_yield: Dual::<1>::variable(
                    clm_spec.steady_state_dividend_yield,
                    0,
                ),
                expected_dividend_growth: Dual::<1>::from_value(clm_spec.expected_dividend_growth),
                expected_discount_rate: Dual::<1>::from_value(clm_spec.expected_discount_rate),
                dividend_growth_ar1_phi: Dual::<1>::from_value(clm_spec.dividend_growth_ar1_phi),
                discount_rate_ar1_phi: Dual::<1>::from_value(clm_spec.discount_rate_ar1_phi),
            },
        );
        let mut cp = clm_spec.clone();
        cp.steady_state_dividend_yield += h;
        let mut cm = clm_spec;
        cm.steady_state_dividend_yield -= h;
        let _fd_clm = (evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller(&cp)
            .equilibrium_log_dividend_price_ratio
            - evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller(&cm)
                .equilibrium_log_dividend_price_ratio)
            / (2.0 * h);
        assert!((clm_ad.equilibrium_log_dividend_price_ratio.d[0] - _fd_clm).abs() < 1e-5);

        // 3. Roncalli Spinu ERC Risk Parity
        let erc_spec = RoncalliSpinuErcRiskParitySpec {
            covariance_3x3: [
                [0.040, 0.008, 0.002],
                [0.008, 0.0225, 0.003],
                [0.002, 0.003, 0.0064],
            ],
            target_risk_budgets: [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        };
        let erc_rep = evaluate_roncalli_spinu_erc_risk_parity(&erc_spec);
        assert!(erc_rep.risk_budget_parity_l1_error < 1e-8);
        assert!(erc_rep.weights_sum_conservation_error < 1e-12);
        assert!(erc_rep.choueifaty_diversification_ratio > 1.0);
        let mut cov_ad = erc_spec
            .covariance_3x3
            .map(|r| r.map(Dual::<1>::from_value));
        cov_ad[0][0] = Dual::<1>::variable(erc_spec.covariance_3x3[0][0], 0);
        let erc_ad = evaluate_roncalli_spinu_erc_risk_parity(&RoncalliSpinuErcRiskParitySpec {
            covariance_3x3: cov_ad,
            target_risk_budgets: erc_spec.target_risk_budgets.map(Dual::<1>::from_value),
        });
        let mut ep = erc_spec;
        ep.covariance_3x3[0][0] += h;
        let mut em = erc_spec;
        em.covariance_3x3[0][0] -= h;
        let _fd_erc = (evaluate_roncalli_spinu_erc_risk_parity(&ep).erc_portfolio_volatility
            - evaluate_roncalli_spinu_erc_risk_parity(&em).erc_portfolio_volatility)
            / (2.0 * h);
        assert!((erc_ad.erc_portfolio_volatility.d[0] - _fd_erc).abs() < 1e-5);

        // 4. Ilmanen Carry, Roll-Down & VMC Factors
        let ilm_spec = IlmanenCarryRolldownVmcFactorSpec {
            current_yield_yt: 0.045,
            rolled_horizon_yield: 0.041,
            repo_financing_rate: 0.032,
            modified_duration: 7.5,
            bond_convexity: 68.0,
            expected_yield_change: 0.001,
            yield_volatility: 0.0085,
            value_style_sharpe: 0.45,
            momentum_style_sharpe: 0.50,
            carry_style_sharpe: 0.48,
            value_momentum_corr_rho: -0.45,
        };
        let ilm_rep = evaluate_ilmanen_carry_rolldown_and_vmc_factors(&ilm_spec);
        assert!(ilm_rep.rolldown_return > 0.0);
        assert!(ilm_rep.composite_vmc_sharpe_ratio > ilm_rep.standalone_average_style_sharpe);
        let ilm_ad =
            evaluate_ilmanen_carry_rolldown_and_vmc_factors(&IlmanenCarryRolldownVmcFactorSpec {
                current_yield_yt: Dual::<1>::variable(ilm_spec.current_yield_yt, 0),
                rolled_horizon_yield: Dual::<1>::from_value(ilm_spec.rolled_horizon_yield),
                repo_financing_rate: Dual::<1>::from_value(ilm_spec.repo_financing_rate),
                modified_duration: Dual::<1>::from_value(ilm_spec.modified_duration),
                bond_convexity: Dual::<1>::from_value(ilm_spec.bond_convexity),
                expected_yield_change: Dual::<1>::from_value(ilm_spec.expected_yield_change),
                yield_volatility: Dual::<1>::from_value(ilm_spec.yield_volatility),
                value_style_sharpe: Dual::<1>::from_value(ilm_spec.value_style_sharpe),
                momentum_style_sharpe: Dual::<1>::from_value(ilm_spec.momentum_style_sharpe),
                carry_style_sharpe: Dual::<1>::from_value(ilm_spec.carry_style_sharpe),
                value_momentum_corr_rho: Dual::<1>::from_value(ilm_spec.value_momentum_corr_rho),
            });
        let mut ip = ilm_spec;
        ip.current_yield_yt += h;
        let mut im = ilm_spec;
        im.current_yield_yt -= h;
        let _fd_ilm = (evaluate_ilmanen_carry_rolldown_and_vmc_factors(&ip)
            .total_expected_excess_return
            - evaluate_ilmanen_carry_rolldown_and_vmc_factors(&im).total_expected_excess_return)
            / (2.0 * h);
        assert!((ilm_ad.total_expected_excess_return.d[0] - _fd_ilm).abs() < 1e-5);
    }
}
