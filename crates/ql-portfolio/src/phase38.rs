#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct CampbellFamaMacbethShankenSpec<T: Real = f64> {
    pub asset_mean_excess_returns: [T; 4],
    pub estimated_factor_betas: [T; 4],
    pub idiosyncratic_var_sigma_eps2: T,
    pub factor_var_sigma_f2: T,
    pub sample_periods_t: T,
}

#[derive(Clone, Copy, Debug)]
pub struct CampbellFamaMacbethShankenReport<T: Real = f64> {
    pub zero_beta_rate_gamma0: T,
    pub factor_risk_premium_gamma1: T,
    pub pricing_error_alphas: [T; 4],
    pub cross_sectional_orthogonality_error: T,
    pub shanken_eiv_correction_c: T,
    pub naive_fama_macbeth_se: T,
    pub shanken_corrected_se: T,
    pub shanken_corrected_t_stat: T,
}

pub fn evaluate_campbell_fama_macbeth_and_shanken_eiv<T: Real>(
    spec: &CampbellFamaMacbethShankenSpec<T>,
) -> CampbellFamaMacbethShankenReport<T> {
    let r = &spec.asset_mean_excess_returns;
    let b = &spec.estimated_factor_betas;
    let q = T::from_f64(0.25);
    let r_bar = q * (r[0] + r[1] + r[2] + r[3]);
    let b_bar = q * (b[0] + b[1] + b[2] + b[3]);

    let mut num = T::zero();
    let mut den = T::zero();
    for i in 0..4 {
        let db = b[i] - b_bar;
        let dr = r[i] - r_bar;
        num = num + db * dr;
        den = den + db * db;
    }
    let den_safe = den.max_val(T::from_f64(1e-8));
    let gamma1 = num / den_safe;
    let gamma0 = r_bar - gamma1 * b_bar;

    let mut alphas = [T::zero(); 4];
    let mut sum_a = T::zero();
    let mut sum_ab = T::zero();
    for i in 0..4 {
        alphas[i] = r[i] - (gamma0 + gamma1 * b[i]);
        sum_a = sum_a + alphas[i];
        sum_ab = sum_ab + alphas[i] * b[i];
    }
    let orth_err = (sum_a * sum_a + sum_ab * sum_ab + T::from_f64(1e-30)).sqrt();

    let sig_f2 = spec.factor_var_sigma_f2.max_val(T::from_f64(1e-8));
    let sig_e2 = spec.idiosyncratic_var_sigma_eps2.max_val(T::from_f64(1e-8));
    let t_obs = spec.sample_periods_t.max_val(T::from_f64(1.0));
    let c_shanken = (gamma1 * gamma1) / sig_f2;
    let var_naive = (sig_e2 / den_safe + sig_f2) / t_obs;
    let var_shank = ((T::one() + c_shanken) * sig_e2 / den_safe + sig_f2) / t_obs;
    let se_naive = var_naive.sqrt();
    let se_shank = var_shank.sqrt();
    let t_stat = gamma1 / se_shank.max_val(T::from_f64(1e-8));

    CampbellFamaMacbethShankenReport {
        zero_beta_rate_gamma0: gamma0,
        factor_risk_premium_gamma1: gamma1,
        pricing_error_alphas: alphas,
        cross_sectional_orthogonality_error: orth_err,
        shanken_eiv_correction_c: c_shanken,
        naive_fama_macbeth_se: se_naive,
        shanken_corrected_se: se_shank,
        shanken_corrected_t_stat: t_stat,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BouchaudCrossImpactEigenmodeSpec<T: Real = f64> {
    pub parent_order_q1: T,
    pub parent_order_q2: T,
    pub self_impact_lambda1: T,
    pub self_impact_lambda2: T,
    pub cross_impact_lambda12: T,
    pub resilience_decay_rho: T,
    pub risk_aversion_gamma: T,
    pub asset_volatility: T,
}

#[derive(Clone, Copy, Debug)]
pub struct BouchaudCrossImpactEigenmodeReport<T: Real = f64> {
    pub market_eigenvalue_mu_plus: T,
    pub spread_eigenvalue_mu_minus: T,
    pub eigen_trace_det_conservation_error: T,
    pub self_impact_cost_component: T,
    pub cross_asset_spillover_cost: T,
    pub total_cross_impact_execution_cost: T,
    pub market_eigenmode_cost_share: T,
}

pub fn evaluate_bouchaud_cross_impact_eigenmode_execution<T: Real>(
    spec: &BouchaudCrossImpactEigenmodeSpec<T>,
) -> BouchaudCrossImpactEigenmodeReport<T> {
    let l1 = spec.self_impact_lambda1.max_val(T::from_f64(1e-6));
    let l2 = spec.self_impact_lambda2.max_val(T::from_f64(1e-6));
    let l12 = spec.cross_impact_lambda12;
    let tr = l1 + l2;
    let det = l1 * l2 - l12 * l12;
    let disc = ((l1 - l2) * (l1 - l2) + T::from_f64(4.0) * l12 * l12 + T::from_f64(1e-24)).sqrt();
    let mu_p = T::from_f64(0.5) * (tr + disc);
    let mu_m = T::from_f64(0.5) * (tr - disc);

    let d_tr = (mu_p + mu_m) - tr;
    let d_det = (mu_p * mu_m) - det;
    let eig_err = (d_tr * d_tr + d_det * d_det + T::from_f64(1e-30)).sqrt();

    let q1 = spec.parent_order_q1;
    let q2 = spec.parent_order_q2;
    let rho = spec
        .resilience_decay_rho
        .max_val(T::from_f64(0.05))
        .min_val(T::from_f64(0.95));
    let kernel_mult = (T::one() + rho) / (T::from_f64(2.0) * (T::one() - rho));
    let risk_pen = T::from_f64(0.25)
        * spec.risk_aversion_gamma
        * spec.asset_volatility
        * spec.asset_volatility
        * (q1 * q1 + q2 * q2);

    let self_cost = kernel_mult * (l1 * q1 * q1 + l2 * q2 * q2) + risk_pen;
    let cross_cost = kernel_mult * T::from_f64(2.0) * l12 * q1 * q2;
    let tot_cost = (self_cost + cross_cost).max_val(T::from_f64(1e-8));

    let v_p1 = l12;
    let v_p2 = mu_p - l1;
    let norm_vp = (v_p1 * v_p1 + v_p2 * v_p2 + T::from_f64(1e-16)).sqrt();
    let proj_p = (q1 * v_p1 + q2 * v_p2) / norm_vp;
    let mkt_mode_cost = kernel_mult * mu_p * proj_p * proj_p;
    let share = (mkt_mode_cost / tot_cost)
        .max_val(T::zero())
        .min_val(T::one());

    BouchaudCrossImpactEigenmodeReport {
        market_eigenvalue_mu_plus: mu_p,
        spread_eigenvalue_mu_minus: mu_m,
        eigen_trace_det_conservation_error: eig_err,
        self_impact_cost_component: self_cost,
        cross_asset_spillover_cost: cross_cost,
        total_cross_impact_execution_cost: tot_cost,
        market_eigenmode_cost_share: share,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct KarnoskySingerMultiCurrencySpec<T: Real = f64> {
    pub port_asset_weights: [T; 3],
    pub bench_asset_weights: [T; 3],
    pub port_fx_weights: [T; 3],
    pub bench_fx_weights: [T; 3],
    pub local_port_returns: [T; 3],
    pub local_bench_returns: [T; 3],
    pub local_eurodeposit_rates: [T; 3],
    pub base_cash_rate: T,
    pub spot_fx_returns: [T; 3],
}

#[derive(Clone, Copy, Debug)]
pub struct KarnoskySingerMultiCurrencyReport<T: Real = f64> {
    pub local_excess_allocation_effect: T,
    pub local_security_selection_effect: T,
    pub currency_forward_carry_effect: T,
    pub spot_currency_surprise_effect: T,
    pub total_active_base_return: T,
    pub karnosky_singer_conservation_error: T,
}

pub fn evaluate_karnosky_singer_multicurrency_attribution<T: Real>(
    spec: &KarnoskySingerMultiCurrencySpec<T>,
) -> KarnoskySingerMultiCurrencyReport<T> {
    let mut r_b_star = [T::zero(); 3];
    let mut r_p_star = [T::zero(); 3];
    let mut d = [T::zero(); 3];
    let mut r_b_bar = T::zero();
    let mut d_bar = T::zero();
    let mut e_bar = T::zero();

    for i in 0..3 {
        r_b_star[i] = spec.local_bench_returns[i] - spec.local_eurodeposit_rates[i];
        r_p_star[i] = spec.local_port_returns[i] - spec.local_eurodeposit_rates[i];
        d[i] = spec.local_eurodeposit_rates[i] - spec.base_cash_rate;
        r_b_bar = r_b_bar + spec.bench_asset_weights[i] * r_b_star[i];
        d_bar = d_bar + spec.bench_fx_weights[i] * d[i];
        e_bar = e_bar + spec.bench_fx_weights[i] * spec.spot_fx_returns[i];
    }

    let mut alloc = T::zero();
    let mut sel = T::zero();
    let mut carry = T::zero();
    let mut spot_fx = T::zero();
    let mut port_tot = T::zero();
    let mut bench_tot = T::zero();

    for i in 0..3 {
        let dw_a = spec.port_asset_weights[i] - spec.bench_asset_weights[i];
        let dw_fx = spec.port_fx_weights[i] - spec.bench_fx_weights[i];
        alloc = alloc + dw_a * (r_b_star[i] - r_b_bar);
        sel = sel + spec.port_asset_weights[i] * (r_p_star[i] - r_b_star[i]);
        carry = carry + dw_fx * (d[i] - d_bar);
        spot_fx = spot_fx + dw_fx * (spec.spot_fx_returns[i] - e_bar);

        port_tot = port_tot
            + spec.port_asset_weights[i] * r_p_star[i]
            + spec.port_fx_weights[i] * (d[i] + spec.spot_fx_returns[i]);
        bench_tot = bench_tot
            + spec.bench_asset_weights[i] * r_b_star[i]
            + spec.bench_fx_weights[i] * (d[i] + spec.spot_fx_returns[i]);
    }
    let active_tot = port_tot - bench_tot;
    let sum_effects = alloc + sel + carry + spot_fx;
    let d_ks = active_tot - sum_effects;
    let ks_err = (d_ks * d_ks + T::from_f64(1e-30)).sqrt();

    KarnoskySingerMultiCurrencyReport {
        local_excess_allocation_effect: alloc,
        local_security_selection_effect: sel,
        currency_forward_carry_effect: carry,
        spot_currency_surprise_effect: spot_fx,
        total_active_base_return: active_tot,
        karnosky_singer_conservation_error: ks_err,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CampbellViceiraIntertemporalSaaSpec<T: Real = f64> {
    pub current_expected_excess_return_xt: T,
    pub long_run_excess_return_xbar: T,
    pub predictor_persistence_phi: T,
    pub return_volatility_sigma_r: T,
    pub predictor_volatility_sigma_x: T,
    pub return_predictor_corr_rho: T,
    pub relative_risk_aversion_gamma: T,
    pub wealth_consumption_rho_w: T,
}

#[derive(Clone, Copy, Debug)]
pub struct CampbellViceiraIntertemporalSaaReport<T: Real = f64> {
    pub myopic_equity_allocation: T,
    pub intertemporal_hedging_demand: T,
    pub total_strategic_equity_allocation: T,
    pub hedging_share_of_total_demand: T,
    pub horizon_variance_reduction_ratio: T,
    pub certainty_equivalent_welfare_gain: T,
}

pub fn evaluate_campbell_viceira_intertemporal_hedging_saa<T: Real>(
    spec: &CampbellViceiraIntertemporalSaaSpec<T>,
) -> CampbellViceiraIntertemporalSaaReport<T> {
    let xt = T::from_f64(0.7) * spec.current_expected_excess_return_xt
        + T::from_f64(0.3) * spec.long_run_excess_return_xbar;
    let sig_r = spec.return_volatility_sigma_r.max_val(T::from_f64(0.02));
    let sig_r2 = sig_r * sig_r;
    let gamma = spec.relative_risk_aversion_gamma.max_val(T::from_f64(1.05));
    let myopic = (xt + T::from_f64(0.5) * sig_r2) / (gamma * sig_r2);

    let denom_cv = (T::one() - spec.wealth_consumption_rho_w * spec.predictor_persistence_phi)
        .max_val(T::from_f64(0.02));
    let psi_hedge = -(T::one() - T::one() / gamma)
        * (spec.wealth_consumption_rho_w
            * spec.return_predictor_corr_rho
            * sig_r
            * spec.predictor_volatility_sigma_x)
        / denom_cv;
    let hedge = (psi_hedge / sig_r2) * (xt / gamma);
    let total = myopic + hedge;
    let share = hedge / total.max_val(T::from_f64(1e-6));

    let var_red = (T::one()
        + T::from_f64(2.0)
            * spec.return_predictor_corr_rho
            * (spec.predictor_volatility_sigma_x / sig_r)
            / (T::one() - spec.predictor_persistence_phi).max_val(T::from_f64(0.05)))
    .max_val(T::from_f64(0.2));
    let ce_gain = T::from_f64(0.5) * gamma * sig_r2 * hedge * hedge;

    CampbellViceiraIntertemporalSaaReport {
        myopic_equity_allocation: myopic,
        intertemporal_hedging_demand: hedge,
        total_strategic_equity_allocation: total,
        hedging_share_of_total_demand: share,
        horizon_variance_reduction_ratio: var_red,
        certainty_equivalent_welfare_gain: ce_gain,
    }
}

#[cfg(test)]
mod phase38_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn phase38_portfolio_all_four_engines_invariants_and_ad_parity() {
        let h = 1e-6;

        // 1. Campbell-Fama-MacBeth-Shanken EIV
        let fm_spec = CampbellFamaMacbethShankenSpec {
            asset_mean_excess_returns: [0.045, 0.062, 0.078, 0.095],
            estimated_factor_betas: [0.70, 0.95, 1.20, 1.45],
            idiosyncratic_var_sigma_eps2: 0.015,
            factor_var_sigma_f2: 0.022,
            sample_periods_t: 120.0,
        };
        let fm_rep = evaluate_campbell_fama_macbeth_and_shanken_eiv(&fm_spec);
        assert!(fm_rep.cross_sectional_orthogonality_error < 1e-12);
        assert!(fm_rep.shanken_corrected_se > fm_rep.naive_fama_macbeth_se);
        let fm_ad =
            evaluate_campbell_fama_macbeth_and_shanken_eiv(&CampbellFamaMacbethShankenSpec {
                asset_mean_excess_returns: [
                    Dual::<1>::variable(fm_spec.asset_mean_excess_returns[0], 0),
                    Dual::<1>::from_value(fm_spec.asset_mean_excess_returns[1]),
                    Dual::<1>::from_value(fm_spec.asset_mean_excess_returns[2]),
                    Dual::<1>::from_value(fm_spec.asset_mean_excess_returns[3]),
                ],
                estimated_factor_betas: fm_spec.estimated_factor_betas.map(Dual::<1>::from_value),
                idiosyncratic_var_sigma_eps2: Dual::<1>::from_value(
                    fm_spec.idiosyncratic_var_sigma_eps2,
                ),
                factor_var_sigma_f2: Dual::<1>::from_value(fm_spec.factor_var_sigma_f2),
                sample_periods_t: Dual::<1>::from_value(fm_spec.sample_periods_t),
            });
        let mut fmp = fm_spec;
        fmp.asset_mean_excess_returns[0] += h;
        let mut fmm = fm_spec;
        fmm.asset_mean_excess_returns[0] -= h;
        let fd_fm = (evaluate_campbell_fama_macbeth_and_shanken_eiv(&fmp)
            .factor_risk_premium_gamma1
            - evaluate_campbell_fama_macbeth_and_shanken_eiv(&fmm).factor_risk_premium_gamma1)
            / (2.0 * h);
        assert!((fm_ad.factor_risk_premium_gamma1.d[0] - fd_fm).abs() < 1e-5);

        // 2. Bouchaud Cross-Impact Eigenmode Execution
        let ci_spec = BouchaudCrossImpactEigenmodeSpec {
            parent_order_q1: 50.0,
            parent_order_q2: 30.0,
            self_impact_lambda1: 0.04,
            self_impact_lambda2: 0.03,
            cross_impact_lambda12: 0.015,
            resilience_decay_rho: 0.60,
            risk_aversion_gamma: 0.10,
            asset_volatility: 0.20,
        };
        let ci_rep = evaluate_bouchaud_cross_impact_eigenmode_execution(&ci_spec);
        assert!(ci_rep.eigen_trace_det_conservation_error < 1e-12);
        assert!(ci_rep.market_eigenvalue_mu_plus > ci_rep.spread_eigenvalue_mu_minus);
        let ci_ad =
            evaluate_bouchaud_cross_impact_eigenmode_execution(&BouchaudCrossImpactEigenmodeSpec {
                parent_order_q1: Dual::<1>::variable(ci_spec.parent_order_q1, 0),
                parent_order_q2: Dual::<1>::from_value(ci_spec.parent_order_q2),
                self_impact_lambda1: Dual::<1>::from_value(ci_spec.self_impact_lambda1),
                self_impact_lambda2: Dual::<1>::from_value(ci_spec.self_impact_lambda2),
                cross_impact_lambda12: Dual::<1>::from_value(ci_spec.cross_impact_lambda12),
                resilience_decay_rho: Dual::<1>::from_value(ci_spec.resilience_decay_rho),
                risk_aversion_gamma: Dual::<1>::from_value(ci_spec.risk_aversion_gamma),
                asset_volatility: Dual::<1>::from_value(ci_spec.asset_volatility),
            });
        let mut cip = ci_spec;
        cip.parent_order_q1 += h;
        let mut cim = ci_spec;
        cim.parent_order_q1 -= h;
        let fd_ci = (evaluate_bouchaud_cross_impact_eigenmode_execution(&cip)
            .total_cross_impact_execution_cost
            - evaluate_bouchaud_cross_impact_eigenmode_execution(&cim)
                .total_cross_impact_execution_cost)
            / (2.0 * h);
        assert!((ci_ad.total_cross_impact_execution_cost.d[0] - fd_ci).abs() < 1e-5);

        // 3. Karnosky-Singer Multi-Currency Attribution
        let ks_spec = KarnoskySingerMultiCurrencySpec {
            port_asset_weights: [0.45, 0.35, 0.20],
            bench_asset_weights: [0.40, 0.35, 0.25],
            port_fx_weights: [0.50, 0.30, 0.20],
            bench_fx_weights: [0.40, 0.35, 0.25],
            local_port_returns: [0.065, 0.048, 0.032],
            local_bench_returns: [0.058, 0.045, 0.035],
            local_eurodeposit_rates: [0.040, 0.028, 0.012],
            base_cash_rate: 0.040,
            spot_fx_returns: [0.0, 0.015, -0.012],
        };
        let ks_rep = evaluate_karnosky_singer_multicurrency_attribution(&ks_spec);
        assert!(ks_rep.karnosky_singer_conservation_error < 1e-12);
        let ks_ad =
            evaluate_karnosky_singer_multicurrency_attribution(&KarnoskySingerMultiCurrencySpec {
                port_asset_weights: ks_spec.port_asset_weights.map(Dual::<1>::from_value),
                bench_asset_weights: ks_spec.bench_asset_weights.map(Dual::<1>::from_value),
                port_fx_weights: ks_spec.port_fx_weights.map(Dual::<1>::from_value),
                bench_fx_weights: ks_spec.bench_fx_weights.map(Dual::<1>::from_value),
                local_port_returns: [
                    Dual::<1>::variable(ks_spec.local_port_returns[0], 0),
                    Dual::<1>::from_value(ks_spec.local_port_returns[1]),
                    Dual::<1>::from_value(ks_spec.local_port_returns[2]),
                ],
                local_bench_returns: ks_spec.local_bench_returns.map(Dual::<1>::from_value),
                local_eurodeposit_rates: ks_spec.local_eurodeposit_rates.map(Dual::<1>::from_value),
                base_cash_rate: Dual::<1>::from_value(ks_spec.base_cash_rate),
                spot_fx_returns: ks_spec.spot_fx_returns.map(Dual::<1>::from_value),
            });
        let mut ksp = ks_spec;
        ksp.local_port_returns[0] += h;
        let mut ksm = ks_spec;
        ksm.local_port_returns[0] -= h;
        let fd_ks = (evaluate_karnosky_singer_multicurrency_attribution(&ksp)
            .total_active_base_return
            - evaluate_karnosky_singer_multicurrency_attribution(&ksm).total_active_base_return)
            / (2.0 * h);
        assert!((ks_ad.total_active_base_return.d[0] - fd_ks).abs() < 1e-5);

        // 4. Campbell-Viceira Intertemporal Hedging SAA
        let cv_spec = CampbellViceiraIntertemporalSaaSpec {
            current_expected_excess_return_xt: 0.055,
            long_run_excess_return_xbar: 0.050,
            predictor_persistence_phi: 0.88,
            return_volatility_sigma_r: 0.16,
            predictor_volatility_sigma_x: 0.035,
            return_predictor_corr_rho: -0.68,
            relative_risk_aversion_gamma: 4.0,
            wealth_consumption_rho_w: 0.96,
        };
        let cv_rep = evaluate_campbell_viceira_intertemporal_hedging_saa(&cv_spec);
        assert!(cv_rep.intertemporal_hedging_demand > 0.0);
        assert!(cv_rep.total_strategic_equity_allocation > cv_rep.myopic_equity_allocation);
        let cv_ad = evaluate_campbell_viceira_intertemporal_hedging_saa(
            &CampbellViceiraIntertemporalSaaSpec {
                current_expected_excess_return_xt: Dual::<1>::variable(
                    cv_spec.current_expected_excess_return_xt,
                    0,
                ),
                long_run_excess_return_xbar: Dual::<1>::from_value(
                    cv_spec.long_run_excess_return_xbar,
                ),
                predictor_persistence_phi: Dual::<1>::from_value(cv_spec.predictor_persistence_phi),
                return_volatility_sigma_r: Dual::<1>::from_value(cv_spec.return_volatility_sigma_r),
                predictor_volatility_sigma_x: Dual::<1>::from_value(
                    cv_spec.predictor_volatility_sigma_x,
                ),
                return_predictor_corr_rho: Dual::<1>::from_value(cv_spec.return_predictor_corr_rho),
                relative_risk_aversion_gamma: Dual::<1>::from_value(
                    cv_spec.relative_risk_aversion_gamma,
                ),
                wealth_consumption_rho_w: Dual::<1>::from_value(cv_spec.wealth_consumption_rho_w),
            },
        );
        let mut cvp = cv_spec;
        cvp.current_expected_excess_return_xt += h;
        let mut cvm = cv_spec;
        cvm.current_expected_excess_return_xt -= h;
        let fd_cv = (evaluate_campbell_viceira_intertemporal_hedging_saa(&cvp)
            .total_strategic_equity_allocation
            - evaluate_campbell_viceira_intertemporal_hedging_saa(&cvm)
                .total_strategic_equity_allocation)
            / (2.0 * h);
        assert!((cv_ad.total_strategic_equity_allocation.d[0] - fd_cv).abs() < 1e-5);
    }
}
