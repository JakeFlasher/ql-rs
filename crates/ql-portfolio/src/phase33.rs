use ql_ad::Real;
#[derive(Clone, Debug)]
pub struct OharaGlostenMilgromSequentialSpec<T: Real = f64> {
    pub low_value_vl: T,
    pub high_value_vh: T,
    pub initial_prior_p0: T,
    pub informed_trader_fraction_mu: T,
    pub informed_signal_accuracy_theta: T,
    pub dealer_inventory_gamma: T,
    pub signed_order_flow: Vec<T>,
}
#[derive(Clone, Copy, Debug)]
pub struct OharaGlostenMilgromSequentialReport<T: Real = f64> {
    pub initial_glosten_milgrom_ask: T,
    pub initial_glosten_milgrom_bid: T,
    pub initial_adverse_selection_spread: T,
    pub terminal_posterior_belief_p: T,
    pub terminal_adverse_selection_spread: T,
    pub terminal_dealer_inventory: T,
    pub inventory_skewed_reservation_mid: T,
    pub information_learning_entropy_reduction: T,
}
pub fn evaluate_ohara_glosten_milgrom_sequential_spread<T: Real>(
    s: &OharaGlostenMilgromSequentialSpec<T>,
) -> OharaGlostenMilgromSequentialReport<T> {
    let d = s.high_value_vh - s.low_value_vl;
    let b = s.informed_trader_fraction_mu;
    let th = s.informed_signal_accuracy_theta;
    let mut p = s.initial_prior_p0;
    let mut q = T::zero();
    let buy_h = b * th + T::from_f64(0.5) * (T::one() - b);
    let buy_l = b * (T::one() - th) + T::from_f64(0.5) * (T::one() - b);
    let pb0 = (p * buy_h + (T::one() - p) * buy_l).max_val(T::from_f64(1e-12));
    let ps0 = (T::one() - pb0).max_val(T::from_f64(1e-12));
    let ask0 = s.low_value_vl + (p * buy_h / pb0) * d;
    let bid0 = s.low_value_vl + (p * buy_l / ps0) * d;
    let h0 = -(p * p.ln() + (T::one() - p) * (T::one() - p).ln());
    for &flow in &s.signed_order_flow {
        let pb = (p * buy_h + (T::one() - p) * buy_l).max_val(T::from_f64(1e-12));
        let ps = (T::one() - pb).max_val(T::from_f64(1e-12));
        let pp = p * buy_h / pb;
        let pm = p * buy_l / ps;
        p = (T::one() + flow) / T::from_f64(2.0) * pp + (T::one() - flow) / T::from_f64(2.0) * pm;
        q -= flow;
    }
    let pb_t = (p * buy_h + (T::one() - p) * buy_l).max_val(T::from_f64(1e-12));
    let ps_t = (T::one() - pb_t).max_val(T::from_f64(1e-12));
    let ask = s.low_value_vl + (p * buy_h / pb_t) * d;
    let bid = s.low_value_vl + (p * buy_l / ps_t) * d;
    let mid = s.low_value_vl + p * d - s.dealer_inventory_gamma * d * d * p * (T::one() - p) * q;
    let h = -(p * p.ln() + (T::one() - p) * (T::one() - p).ln());
    OharaGlostenMilgromSequentialReport {
        initial_glosten_milgrom_ask: ask0,
        initial_glosten_milgrom_bid: bid0,
        initial_adverse_selection_spread: ask0 - bid0,
        terminal_posterior_belief_p: p,
        terminal_adverse_selection_spread: ask - bid,
        terminal_dealer_inventory: q,
        inventory_skewed_reservation_mid: mid,
        information_learning_entropy_reduction: h0 - h,
    }
}
#[derive(Clone, Debug)]
pub struct PottersBouchaudStieltjesRieSpec<T: Real = f64> {
    pub sample_eigenvalues: Vec<T>,
    pub aspect_ratio_q: T,
    pub stieltjes_bandwidth_eta: T,
}
#[derive(Clone, Debug)]
pub struct PottersBouchaudStieltjesRieReport<T: Real = f64> {
    pub rie_shrunk_eigenvalues: Vec<T>,
    pub sample_min_eigenvalue: T,
    pub rie_min_eigenvalue: T,
    pub sample_condition_number: T,
    pub rie_condition_number: T,
    pub trace_conservation_error: T,
    pub markowitz_oos_variance_inflation_factor: T,
    pub inverse_variance_portfolio_risk_reduction: T,
}
pub fn evaluate_potters_bouchaud_stieltjes_rie_shrinkage<T: Real>(
    s: &PottersBouchaudStieltjesRieSpec<T>,
) -> PottersBouchaudStieltjesRieReport<T> {
    let n = s.sample_eigenvalues.len();
    let mut shr = Vec::with_capacity(n);
    for i in 0..n {
        let l = s.sample_eigenvalues[i];
        let eta = s.stieltjes_bandwidth_eta * l / (n as f64).sqrt();
        let mut re = T::zero();
        let mut im = T::zero();
        for &lj in &s.sample_eigenvalues {
            let d = lj - l;
            let den = d * d + eta * eta;
            re += d / den;
            im += eta / den;
        }
        re /= T::from_f64(n as f64);
        im /= T::from_f64(n as f64);
        let a = T::one() - s.aspect_ratio_q + s.aspect_ratio_q * (l * re + eta * im);
        let b = s.aspect_ratio_q * (l * im - eta * re);
        shr.push(l / (a * a + b * b));
    }
    let suml = s
        .sample_eigenvalues
        .iter()
        .copied()
        .fold(T::zero(), |a, x| a + x);
    let sums = shr.iter().copied().fold(T::zero(), |a, x| a + x);
    for x in &mut shr {
        *x *= suml / sums;
    }
    let min_s = s
        .sample_eigenvalues
        .iter()
        .copied()
        .reduce(|a, b| a.min_val(b))
        .unwrap();
    let max_s = s
        .sample_eigenvalues
        .iter()
        .copied()
        .reduce(|a, b| a.max_val(b))
        .unwrap();
    let min_r = shr.iter().copied().reduce(|a, b| a.min_val(b)).unwrap();
    let max_r = shr.iter().copied().reduce(|a, b| a.max_val(b)).unwrap();
    PottersBouchaudStieltjesRieReport {
        rie_shrunk_eigenvalues: shr.clone(),
        sample_min_eigenvalue: min_s,
        rie_min_eigenvalue: min_r,
        sample_condition_number: max_s / min_s,
        rie_condition_number: max_r / min_r,
        trace_conservation_error: (shr.iter().copied().fold(T::zero(), |a, x| a + x) - suml).abs(),
        markowitz_oos_variance_inflation_factor: T::one() / (T::one() - s.aspect_ratio_q),
        inverse_variance_portfolio_risk_reduction: T::one()
            - min_s / min_r.max_val(T::from_f64(1e-12)),
    }
}
#[derive(Clone, Debug)]
pub struct TsayEgarchHamiltonRegimeSpec<T: Real = f64> {
    pub returns: Vec<T>,
    pub initial_variance: T,
    pub egarch_omega: T,
    pub egarch_alpha: T,
    pub egarch_leverage_gamma: T,
    pub egarch_beta: T,
    pub calm_mean: T,
    pub calm_vol: T,
    pub crisis_mean: T,
    pub crisis_vol: T,
    pub stay_calm_p00: T,
    pub stay_crisis_p11: T,
}
#[derive(Clone, Copy, Debug)]
pub struct TsayEgarchHamiltonRegimeReport<T: Real = f64> {
    pub terminal_egarch_log_variance: T,
    pub terminal_egarch_volatility: T,
    pub egarch_negative_shock_asymmetry_ratio: T,
    pub ergodic_crisis_probability: T,
    pub terminal_filtered_crisis_probability: T,
    pub regime_blended_conditional_vol: T,
    pub hamilton_filter_log_likelihood: T,
}
pub fn evaluate_tsay_egarch_and_hamilton_markov_switching<T: Real>(
    s: &TsayEgarchHamiltonRegimeSpec<T>,
) -> TsayEgarchHamiltonRegimeReport<T> {
    let mut lh = s.initial_variance.ln();
    let mut ll = T::zero();
    let mut pc =
        (T::one() - s.stay_calm_p00) / (T::from_f64(2.0) - s.stay_calm_p00 - s.stay_crisis_p11);
    for &r in &s.returns {
        let h = lh.exp();
        let z = r / h.sqrt();
        lh = s.egarch_omega
            + s.egarch_beta * lh
            + s.egarch_alpha
                * ((z * z + T::from_f64(1e-12)).sqrt()
                    - T::from_f64((2.0 / std::f64::consts::PI).sqrt()))
            + s.egarch_leverage_gamma * z;
        let p0 = s.stay_calm_p00 * (T::one() - pc) + (T::one() - s.stay_crisis_p11) * pc;
        let p1 = T::one() - p0;
        let f0 = ((-T::from_f64(0.5) * (r - s.calm_mean) * (r - s.calm_mean)
            / (s.calm_vol * s.calm_vol))
            .exp())
            / (T::from_f64(2.0 * std::f64::consts::PI).sqrt() * s.calm_vol);
        let f1 = ((-T::from_f64(0.5) * (r - s.crisis_mean) * (r - s.crisis_mean)
            / (s.crisis_vol * s.crisis_vol))
            .exp())
            / (T::from_f64(2.0 * std::f64::consts::PI).sqrt() * s.crisis_vol);
        let f = p0 * f0 + p1 * f1;
        pc = p1 * f1 / f;
        ll += f.ln();
    }
    TsayEgarchHamiltonRegimeReport {
        terminal_egarch_log_variance: lh,
        terminal_egarch_volatility: (lh.exp()).sqrt(),
        egarch_negative_shock_asymmetry_ratio: (s.egarch_alpha - s.egarch_leverage_gamma)
            / (s.egarch_alpha + s.egarch_leverage_gamma).max_val(T::from_f64(1e-6)),
        ergodic_crisis_probability: (T::one() - s.stay_calm_p00)
            / (T::from_f64(2.0) - s.stay_calm_p00 - s.stay_crisis_p11),
        terminal_filtered_crisis_probability: pc,
        regime_blended_conditional_vol: ((T::one() - pc) * s.calm_vol * s.calm_vol
            + pc * s.crisis_vol * s.crisis_vol)
            .sqrt(),
        hamilton_filter_log_likelihood: ll,
    }
}
#[derive(Clone, Copy, Debug)]
pub struct McneilFreyCopulaEulerAllocationSpec<T: Real = f64> {
    pub clayton_theta: T,
    pub gumbel_theta: T,
    pub student_t_df_nu: T,
    pub student_t_corr_rho: T,
    pub desk_weights: [T; 4],
    pub desk_expected_losses: [T; 4],
    pub loss_covariance_4x4: [[T; 4]; 4],
    pub es_tail_multiplier: T,
}
#[derive(Clone, Copy, Debug)]
pub struct McneilFreyCopulaEulerAllocationReport<T: Real = f64> {
    pub clayton_lower_tail_lambda: T,
    pub gumbel_upper_tail_lambda: T,
    pub student_t_tail_lambda: T,
    pub total_portfolio_es: T,
    pub standalone_es_sum: T,
    pub euler_desk_allocations: [T; 4],
    pub portfolio_diversification_benefit: T,
    pub euler_additivity_error: T,
}
pub fn evaluate_mcneil_frey_copula_tail_and_euler_allocation<T: Real>(
    s: &McneilFreyCopulaEulerAllocationSpec<T>,
) -> McneilFreyCopulaEulerAllocationReport<T> {
    let cl = T::from_f64(2.0).pow(-T::one() / s.clayton_theta);
    let gu = T::from_f64(2.0) - T::from_f64(2.0).pow(T::one() / s.gumbel_theta);
    let arg = ((s.student_t_df_nu + T::one()) * (T::one() - s.student_t_corr_rho)
        / (T::one() + s.student_t_corr_rho))
        .sqrt();
    let st = T::from_f64(2.0) * (-arg).norm_cdf();
    let mut sw = [T::zero(); 4];
    for i in 0..4 {
        for j in 0..4 {
            sw[i] += s.loss_covariance_4x4[i][j] * s.desk_weights[j];
        }
    }
    let var = s
        .desk_weights
        .iter()
        .enumerate()
        .map(|(i, w)| *w * sw[i])
        .fold(T::zero(), |a, x| a + x);
    let vol = var.sqrt();
    let mean = s
        .desk_weights
        .iter()
        .zip(s.desk_expected_losses.iter())
        .map(|(w, m)| *w * *m)
        .fold(T::zero(), |a, x| a + x);
    let es = mean + s.es_tail_multiplier * vol;
    let mut rc = [T::zero(); 4];
    for i in 0..4 {
        rc[i] = s.desk_weights[i] * s.desk_expected_losses[i]
            + s.es_tail_multiplier * s.desk_weights[i] * sw[i] / vol;
    }
    let stand = s
        .desk_weights
        .iter()
        .enumerate()
        .map(|(i, w)| {
            *w * s.desk_expected_losses[i]
                + s.es_tail_multiplier * w.abs() * s.loss_covariance_4x4[i][i].sqrt()
        })
        .fold(T::zero(), |a, x| a + x);
    let sum_rc = rc.iter().copied().fold(T::zero(), |a, x| a + x);
    McneilFreyCopulaEulerAllocationReport {
        clayton_lower_tail_lambda: cl,
        gumbel_upper_tail_lambda: gu,
        student_t_tail_lambda: st,
        total_portfolio_es: es,
        standalone_es_sum: stand,
        euler_desk_allocations: rc,
        portfolio_diversification_benefit: stand - es,
        euler_additivity_error: (es - sum_rc).abs(),
    }
}
#[derive(Clone, Debug)]
pub struct BaconOmegaSortinoBurkeDrawdownSpec<T: Real = f64> {
    pub period_returns: Vec<T>,
    pub mar_threshold: T,
    pub risk_free_rate: T,
}
#[derive(Clone, Copy, Debug)]
pub struct BaconOmegaSortinoBurkeDrawdownReport<T: Real = f64> {
    pub mean_return: T,
    pub upside_potential_lpm1: T,
    pub downside_shortfall_lpm1: T,
    pub keating_shadwick_omega_ratio: T,
    pub omega_put_call_parity_error: T,
    pub downside_deviation: T,
    pub sortino_ratio: T,
    pub kappa_three_ratio: T,
    pub ulcer_drawdown_index: T,
    pub martin_ulcer_ratio: T,
    pub burke_drawdown_ratio: T,
}
pub fn evaluate_bacon_omega_sortino_and_burke_drawdown<T: Real>(
    s: &BaconOmegaSortinoBurkeDrawdownSpec<T>,
) -> BaconOmegaSortinoBurkeDrawdownReport<T> {
    let n = s.period_returns.len();
    let mut mean = T::zero();
    let mut up = T::zero();
    let mut dn = T::zero();
    let mut l2 = T::zero();
    let mut l3 = T::zero();
    let mut w = T::one();
    let mut peak = w;
    let mut dds = T::zero();
    for &r in &s.period_returns {
        mean += r;
        let x = r - s.mar_threshold;
        let g = (x + (x * x + T::from_f64(1e-16)).sqrt()) / T::from_f64(2.0);
        let d = (-x + (x * x + T::from_f64(1e-16)).sqrt()) / T::from_f64(2.0);
        up += g;
        dn += d;
        l2 += d * d;
        l3 += d * d * d;
        w *= T::one() + r;
        peak = peak.max_val(w);
        let dd = (peak - w) / peak;
        dds += dd * dd;
    }
    let nf = T::from_f64(n as f64);
    mean /= nf;
    up /= nf;
    dn /= nf;
    let ui = (dds / nf).sqrt();
    let ex = mean - s.risk_free_rate;
    BaconOmegaSortinoBurkeDrawdownReport {
        mean_return: mean,
        upside_potential_lpm1: up,
        downside_shortfall_lpm1: dn,
        keating_shadwick_omega_ratio: up / dn,
        omega_put_call_parity_error: (up - dn - (mean - s.mar_threshold)).abs(),
        downside_deviation: (l2 / nf).sqrt(),
        sortino_ratio: (mean - s.mar_threshold) / (l2 / nf).sqrt(),
        kappa_three_ratio: (mean - s.mar_threshold) / (l3 / nf).powf(1.0 / 3.0),
        ulcer_drawdown_index: ui,
        martin_ulcer_ratio: ex / ui,
        burke_drawdown_ratio: ex / dds.sqrt(),
    }
}
#[cfg(test)]
mod phase33_tests {
    use super::*;
    #[test]
    fn phase33_portfolio_smoke() {
        let s = BaconOmegaSortinoBurkeDrawdownSpec {
            period_returns: vec![0.01, -0.005, 0.02],
            mar_threshold: 0.0,
            risk_free_rate: 0.001,
        };
        assert!(
            evaluate_bacon_omega_sortino_and_burke_drawdown(&s).keating_shadwick_omega_ratio > 1.0
        );
    }
}
#[cfg(test)]
mod phase33_tests_ad {
    use super::*;
    use ql_ad::Dual;
    fn s<T: Real>(x: T) -> OharaGlostenMilgromSequentialSpec<T> {
        OharaGlostenMilgromSequentialSpec {
            low_value_vl: T::from_f64(90.),
            high_value_vh: T::from_f64(110.),
            initial_prior_p0: T::from_f64(0.5) + x,
            informed_trader_fraction_mu: T::from_f64(0.3),
            informed_signal_accuracy_theta: T::from_f64(0.8),
            dealer_inventory_gamma: T::from_f64(0.01),
            signed_order_flow: vec![T::from_f64(0.4), T::from_f64(0.3)],
        }
    }
    #[test]
    fn dual_fd() {
        let h = 1e-6;
        let fd = (evaluate_ohara_glosten_milgrom_sequential_spread(&s(h))
            .inventory_skewed_reservation_mid
            - evaluate_ohara_glosten_milgrom_sequential_spread(&s(-h))
                .inventory_skewed_reservation_mid)
            / (2. * h);
        let ad = evaluate_ohara_glosten_milgrom_sequential_spread(&s(Dual::<1>::variable(0., 0)))
            .inventory_skewed_reservation_mid;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
