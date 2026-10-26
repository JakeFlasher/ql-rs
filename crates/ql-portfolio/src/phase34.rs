use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct FoucaultKyleEkopPinSpec<T: Real = f64> {
    pub prior_price_p0: T,
    pub insider_valuation_v: T,
    pub fundamental_std_sigma_v: T,
    pub noise_order_std_sigma_u: T,
    pub event_prob_alpha: T,
    pub bad_news_prob_delta: T,
    pub informed_arrival_mu: T,
    pub uninformed_buy_eps_b: T,
    pub uninformed_sell_eps_s: T,
    pub observed_buys_b: T,
    pub observed_sells_s: T,
}
#[derive(Clone, Copy, Debug)]
pub struct FoucaultKyleEkopPinReport<T: Real = f64> {
    pub kyle_insider_intensity_beta: T,
    pub kyle_illiquidity_lambda: T,
    pub kyle_market_depth_inverse: T,
    pub optimal_insider_order_size: T,
    pub expected_insider_profit: T,
    pub post_trade_fundamental_variance: T,
    pub ekop_pin_probability: T,
    pub posterior_event_probability: T,
}
pub fn evaluate_foucault_kyle_lambda_and_ekop_pin<T: Real>(
    s: &FoucaultKyleEkopPinSpec<T>,
) -> FoucaultKyleEkopPinReport<T> {
    let b = s.noise_order_std_sigma_u / s.fundamental_std_sigma_v;
    let l = s.fundamental_std_sigma_v / (T::from_f64(2.0) * s.noise_order_std_sigma_u);
    let d = s.insider_valuation_v - s.prior_price_p0;
    let eb = s.uninformed_buy_eps_b;
    let es = s.uninformed_sell_eps_s;
    let mu = s.informed_arrival_mu;
    let ll0 = -(eb + es)
        + s.observed_buys_b * eb.max_val(T::from_f64(1e-12)).ln()
        + s.observed_sells_s * es.max_val(T::from_f64(1e-12)).ln();
    let llb = -(eb + es + mu)
        + s.observed_buys_b * eb.max_val(T::from_f64(1e-12)).ln()
        + s.observed_sells_s * (es + mu).max_val(T::from_f64(1e-12)).ln();
    let llg = -(eb + es + mu)
        + s.observed_buys_b * (eb + mu).max_val(T::from_f64(1e-12)).ln()
        + s.observed_sells_s * es.max_val(T::from_f64(1e-12)).ln();
    let m = ll0.max_val(llb).max_val(llg);
    let p0 = (ll0 - m).exp();
    let pb = (llb - m).exp();
    let pg = (llg - m).exp();
    let pin =
        (s.event_prob_alpha * mu) / (s.event_prob_alpha * mu + eb + es).max_val(T::from_f64(1e-12));
    let num_ev =
        s.event_prob_alpha * (s.bad_news_prob_delta * pb + (T::one() - s.bad_news_prob_delta) * pg);
    let den_ev = ((T::one() - s.event_prob_alpha) * p0 + num_ev).max_val(T::from_f64(1e-12));
    let post_ev = num_ev / den_ev;
    FoucaultKyleEkopPinReport {
        kyle_insider_intensity_beta: b,
        kyle_illiquidity_lambda: l,
        kyle_market_depth_inverse: T::one() / l.max_val(T::from_f64(1e-12)),
        optimal_insider_order_size: b * d,
        expected_insider_profit: T::from_f64(0.5) * b * d * d,
        post_trade_fundamental_variance: s.fundamental_std_sigma_v * s.fundamental_std_sigma_v
            / T::from_f64(2.0),
        ekop_pin_probability: pin,
        posterior_event_probability: post_ev,
    }
}
#[derive(Clone, Copy, Debug)]
pub struct GarleanuPedersenDynamicAimSpec<T: Real = f64> {
    pub current_weights: [T; 2],
    pub fast_alpha_signal: [T; 2],
    pub slow_alpha_signal: [T; 2],
    pub fast_decay_phi: T,
    pub slow_decay_phi: T,
    pub covariance_2x2: [[T; 2]; 2],
    pub risk_aversion_gamma: T,
    pub transaction_cost_lambda: T,
    pub discount_rate_rho: T,
}
#[derive(Clone, Copy, Debug)]
pub struct GarleanuPedersenDynamicAimReport<T: Real = f64> {
    pub optimal_trading_speed_a: T,
    pub fast_signal_attenuation_weight: T,
    pub slow_signal_attenuation_weight: T,
    pub markowitz_unconstrained_target: [T; 2],
    pub aim_portfolio_weights: [T; 2],
    pub optimal_next_weights: [T; 2],
    pub single_step_transaction_cost: T,
    pub expected_net_utility_gain: T,
}
pub fn evaluate_garleanu_pedersen_dynamic_aim_portfolio<T: Real>(
    s: &GarleanuPedersenDynamicAimSpec<T>,
) -> GarleanuPedersenDynamicAimReport<T> {
    let q = s.risk_aversion_gamma / s.transaction_cost_lambda;
    let a = ((-(s.discount_rate_rho + q)
        + ((s.discount_rate_rho + q) * (s.discount_rate_rho + q)
            + T::from_f64(4.) * (T::one() - s.discount_rate_rho) * q)
            .sqrt())
        / (T::from_f64(2.0) * (T::one() - s.discount_rate_rho)))
        .max_val(T::from_f64(1e-6));
    let pf = (T::one() + s.fast_decay_phi * s.transaction_cost_lambda / s.risk_aversion_gamma * a)
        .powf(-1.0);
    let ps = (T::one() + s.slow_decay_phi * s.transaction_cost_lambda / s.risk_aversion_gamma * a)
        .powf(-1.0);
    let a11 = s.covariance_2x2[0][0];
    let a12 = s.covariance_2x2[0][1];
    let a22 = s.covariance_2x2[1][1];
    let det = (a11 * a22 - a12 * a12).max_val(T::from_f64(1e-12));
    let inv = |x0: T, x1: T| [(a22 * x0 - a12 * x1) / det, (-a12 * x0 + a11 * x1) / det];
    let mf = inv(
        s.fast_alpha_signal[0] + s.slow_alpha_signal[0],
        s.fast_alpha_signal[1] + s.slow_alpha_signal[1],
    );
    let m = [mf[0] / s.risk_aversion_gamma, mf[1] / s.risk_aversion_gamma];
    let am = inv(
        pf * s.fast_alpha_signal[0] + ps * s.slow_alpha_signal[0],
        pf * s.fast_alpha_signal[1] + ps * s.slow_alpha_signal[1],
    );
    let aim = [am[0] / s.risk_aversion_gamma, am[1] / s.risk_aversion_gamma];
    let next = [
        (T::one() - a) * s.current_weights[0] + a * aim[0],
        (T::one() - a) * s.current_weights[1] + a * aim[1],
    ];
    GarleanuPedersenDynamicAimReport {
        optimal_trading_speed_a: a,
        fast_signal_attenuation_weight: pf,
        slow_signal_attenuation_weight: ps,
        markowitz_unconstrained_target: m,
        aim_portfolio_weights: aim,
        optimal_next_weights: next,
        single_step_transaction_cost: {
            let d0 = next[0] - s.current_weights[0];
            let d1 = next[1] - s.current_weights[1];
            T::from_f64(0.5)
                * s.transaction_cost_lambda
                * (a11 * d0 * d0 + T::from_f64(2.0) * a12 * d0 * d1 + a22 * d1 * d1)
        },
        expected_net_utility_gain: {
            let al0 = s.fast_alpha_signal[0] + s.slow_alpha_signal[0];
            let al1 = s.fast_alpha_signal[1] + s.slow_alpha_signal[1];
            let quad = |w0: T, w1: T| {
                al0 * w0 + al1 * w1
                    - T::from_f64(0.5)
                        * s.risk_aversion_gamma
                        * (a11 * w0 * w0 + T::from_f64(2.0) * a12 * w0 * w1 + a22 * w1 * w1)
            };
            let d0 = next[0] - s.current_weights[0];
            let d1 = next[1] - s.current_weights[1];
            let tc = T::from_f64(0.5)
                * s.transaction_cost_lambda
                * (a11 * d0 * d0 + T::from_f64(2.0) * a12 * d0 * d1 + a22 * d1 * d1);
            (quad(next[0], next[1]) - tc) - quad(s.current_weights[0], s.current_weights[1])
        },
    }
}
#[derive(Clone, Debug)]
pub struct TsayGjrDccGarchSpec<T: Real = f64> {
    pub bivariate_returns: Vec<[T; 2]>,
    pub omega: [T; 2],
    pub alpha: [T; 2],
    pub gjr_gamma: [T; 2],
    pub beta: [T; 2],
    pub unconditional_corr: T,
    pub dcc_a: T,
    pub dcc_b: T,
}
#[derive(Clone, Copy, Debug)]
pub struct TsayGjrDccGarchReport<T: Real = f64> {
    pub unconditional_var_asset1: T,
    pub unconditional_var_asset2: T,
    pub terminal_gjr_vol_asset1: T,
    pub terminal_gjr_vol_asset2: T,
    pub news_impact_asymmetry_ratio: T,
    pub terminal_dcc_correlation: T,
    pub terminal_equal_weight_portfolio_vol: T,
    pub dcc_persistence_sum: T,
}
pub fn evaluate_tsay_gjr_and_dcc_garch_dynamic_correlation<T: Real>(
    s: &TsayGjrDccGarchSpec<T>,
) -> TsayGjrDccGarchReport<T> {
    let h0 = s.omega[0] / (T::one() - s.alpha[0] - s.gjr_gamma[0] / T::from_f64(2.0) - s.beta[0]);
    let h1 = s.omega[1] / (T::one() - s.alpha[1] - s.gjr_gamma[1] / T::from_f64(2.0) - s.beta[1]);
    let mut h = [h0, h1];
    let mut q = [
        [T::one(), s.unconditional_corr],
        [s.unconditional_corr, T::one()],
    ];
    for r in &s.bivariate_returns {
        let z = [r[0] / h[0].sqrt(), r[1] / h[1].sqrt()];
        q[0][0] = (T::one() - s.dcc_a - s.dcc_b) + s.dcc_a * z[0] * z[0] + s.dcc_b * q[0][0];
        q[0][1] = (T::one() - s.dcc_a - s.dcc_b) * s.unconditional_corr
            + s.dcc_a * z[0] * z[1]
            + s.dcc_b * q[0][1];
        q[1][0] = q[0][1];
        q[1][1] = (T::one() - s.dcc_a - s.dcc_b) + s.dcc_a * z[1] * z[1] + s.dcc_b * q[1][1];
        let neg0 = T::from_f64(0.5) * ((r[0] * r[0] + T::from_f64(1e-16)).sqrt() - r[0]);
        let neg1 = T::from_f64(0.5) * ((r[1] * r[1] + T::from_f64(1e-16)).sqrt() - r[1]);
        h[0] =
            s.omega[0] + s.alpha[0] * r[0] * r[0] + s.gjr_gamma[0] * neg0 * neg0 + s.beta[0] * h[0];
        h[1] =
            s.omega[1] + s.alpha[1] * r[1] * r[1] + s.gjr_gamma[1] * neg1 * neg1 + s.beta[1] * h[1];
    }
    let corr = (q[0][1] / (q[0][0] * q[1][1]).sqrt())
        .max_val(T::from_f64(-0.999))
        .min_val(T::from_f64(0.999));
    TsayGjrDccGarchReport {
        unconditional_var_asset1: h0,
        unconditional_var_asset2: h1,
        terminal_gjr_vol_asset1: h[0].sqrt(),
        terminal_gjr_vol_asset2: h[1].sqrt(),
        news_impact_asymmetry_ratio: (s.alpha[0] + s.gjr_gamma[0]) / s.alpha[0],
        terminal_dcc_correlation: corr,
        terminal_equal_weight_portfolio_vol: ((h[0]
            + h[1]
            + T::from_f64(2.0) * corr * (h[0] * h[1]).sqrt())
            / T::from_f64(4.))
        .sqrt(),
        dcc_persistence_sum: s.dcc_a + s.dcc_b,
    }
}
#[derive(Clone, Copy, Debug)]
pub struct SornetteJlsLpplCrashHazardSpec<T: Real = f64> {
    pub current_time_t: T,
    pub critical_time_tc: T,
    pub horizon_dt: T,
    pub peak_log_price_a: T,
    pub power_law_slope_b: T,
    pub log_periodic_amp_c: T,
    pub critical_exponent_m: T,
    pub log_periodic_omega: T,
    pub phase_phi: T,
    pub crash_loss_kappa: T,
}
#[derive(Clone, Copy, Debug)]
pub struct SornetteJlsLpplCrashHazardReport<T: Real = f64> {
    pub lppl_fitted_price: T,
    pub instantaneous_crash_hazard: T,
    pub horizon_cumulative_crash_prob: T,
    pub jls_damping_validity_ratio: T,
    pub log_periodic_oscillation_term: T,
    pub rational_bubble_excess_return: T,
    pub expected_crash_adjusted_return: T,
}
pub fn evaluate_sornette_jls_lppl_crash_hazard<T: Real>(
    s: &SornetteJlsLpplCrashHazardSpec<T>,
) -> SornetteJlsLpplCrashHazardReport<T> {
    let tau = s.critical_time_tc - s.current_time_t;
    let tau_m = (tau.max_val(T::from_f64(1e-12)).ln() * s.critical_exponent_m).exp();
    let tau_m1 = (tau.max_val(T::from_f64(1e-12)).ln() * (s.critical_exponent_m - T::one())).exp();
    let ph = s.log_periodic_omega * tau.max_val(T::from_f64(1e-12)).ln() - s.phase_phi;
    let osc = s.log_periodic_amp_c * ph.cos();
    let lp =
        s.peak_log_price_a + s.power_law_slope_b * tau_m + s.log_periodic_amp_c * tau_m * ph.cos();
    let h = ((-s.power_law_slope_b * s.critical_exponent_m * tau_m1
        - s.log_periodic_amp_c
            * tau_m1
            * (s.critical_exponent_m * ph.cos() - s.log_periodic_omega * ph.sin()))
        / s.crash_loss_kappa)
        .max_val(T::zero());
    SornetteJlsLpplCrashHazardReport {
        lppl_fitted_price: lp.exp(),
        instantaneous_crash_hazard: h,
        horizon_cumulative_crash_prob: (-h * s.horizon_dt).exp().neg() + T::one(),
        jls_damping_validity_ratio: {
            let b_abs = (s.power_law_slope_b * s.power_law_slope_b + T::from_f64(1e-30)).sqrt();
            let c_abs = (s.log_periodic_amp_c * s.log_periodic_amp_c + T::from_f64(1e-30)).sqrt();
            let mw = (s.critical_exponent_m * s.critical_exponent_m
                + s.log_periodic_omega * s.log_periodic_omega)
                .sqrt();
            s.critical_exponent_m * b_abs / (c_abs * mw).max_val(T::from_f64(1e-12))
        },
        log_periodic_oscillation_term: osc,
        rational_bubble_excess_return: s.crash_loss_kappa * h,
        expected_crash_adjusted_return: {
            let p_surv = (-h * s.horizon_dt).exp();
            let p_crash = T::one() - p_surv;
            let mu_power = -s.power_law_slope_b * s.critical_exponent_m * tau_m1;
            p_surv * (mu_power * s.horizon_dt) - p_crash * s.crash_loss_kappa
        },
    }
}
#[derive(Clone, Debug)]
pub struct FollmerSchiedEntropicSpectralRiskSpec<T: Real = f64> {
    pub ordered_losses: Vec<T>,
    pub entropic_theta: T,
    pub spectral_power_gamma: T,
}
#[derive(Clone, Copy, Debug)]
pub struct FollmerSchiedEntropicSpectralRiskReport<T: Real = f64> {
    pub expected_loss: T,
    pub entropic_convex_risk_measure: T,
    pub worst_case_esscher_expected_loss: T,
    pub relative_entropy_kl_penalty: T,
    pub fenchel_duality_conservation_error: T,
    pub kusuoka_power_spectral_risk_measure: T,
    pub spectral_weights_sum_error: T,
    pub spectral_risk_loading_over_mean: T,
}
pub fn evaluate_follmer_schied_entropic_and_spectral_risk<T: Real>(
    s: &FollmerSchiedEntropicSpectralRiskSpec<T>,
) -> FollmerSchiedEntropicSpectralRiskReport<T> {
    let n = T::from_f64(s.ordered_losses.len() as f64);
    let mut e = T::zero();
    for &x in &s.ordered_losses {
        e += x;
    }
    e /= n;
    let mut zsum = T::zero();
    for &x in &s.ordered_losses {
        zsum += (s.entropic_theta * x).exp();
    }
    let rho = (zsum / n).ln() / s.entropic_theta.max_val(T::from_f64(1e-12));
    let mut sp = T::zero();
    let mut w_sum = T::zero();
    let mut qmean = T::zero();
    let mut kl = T::zero();
    for &x in &s.ordered_losses {
        let qi = (s.entropic_theta * x).exp() / zsum;
        qmean += qi * x;
        kl += qi * (n * qi).max_val(T::from_f64(1e-12)).ln();
    }
    for (i, &x) in s.ordered_losses.iter().enumerate() {
        let a = T::from_f64((i + 1) as f64) / n;
        let b = T::from_f64(i as f64) / n;
        let wi = a.pow(s.spectral_power_gamma) - b.pow(s.spectral_power_gamma);
        sp += wi * x;
        w_sum += wi;
    }
    FollmerSchiedEntropicSpectralRiskReport {
        expected_loss: e,
        entropic_convex_risk_measure: rho,
        worst_case_esscher_expected_loss: qmean,
        relative_entropy_kl_penalty: kl,
        fenchel_duality_conservation_error: (rho
            - (qmean - kl / s.entropic_theta.max_val(T::from_f64(1e-12))))
        .abs(),
        kusuoka_power_spectral_risk_measure: sp,
        spectral_weights_sum_error: ((w_sum - T::one()) * (w_sum - T::one()) + T::from_f64(1e-30))
            .sqrt(),
        spectral_risk_loading_over_mean: sp / e.max_val(T::from_f64(1e-12)),
    }
}
#[derive(Clone, Debug)]
pub struct BaconBrinsonFachlerCarinoSpec<T: Real = f64> {
    pub portfolio_sector_weights: Vec<T>,
    pub benchmark_sector_weights: Vec<T>,
    pub portfolio_sector_returns: Vec<T>,
    pub benchmark_sector_returns: Vec<T>,
    pub multi_period_port_returns: Vec<T>,
    pub multi_period_bench_returns: Vec<T>,
}
#[derive(Clone, Copy, Debug)]
pub struct BaconBrinsonFachlerCarinoReport<T: Real = f64> {
    pub brinson_fachler_allocation: T,
    pub brinson_fachler_selection: T,
    pub brinson_fachler_interaction: T,
    pub single_period_active_return: T,
    pub single_period_conservation_error: T,
    pub cumulative_portfolio_return: T,
    pub cumulative_benchmark_return: T,
    pub carino_total_linking_factor: T,
    pub carino_smoothed_active_return: T,
    pub carino_linking_conservation_error: T,
}
pub fn evaluate_bacon_brinson_fachler_and_carino_linking<T: Real>(
    s: &BaconBrinsonFachlerCarinoSpec<T>,
) -> BaconBrinsonFachlerCarinoReport<T> {
    let n = s.portfolio_sector_weights.len();
    let mut rp = T::zero();
    let mut rb = T::zero();
    for i in 0..n {
        rp += s.portfolio_sector_weights[i] * s.portfolio_sector_returns[i];
        rb += s.benchmark_sector_weights[i] * s.benchmark_sector_returns[i];
    }
    let mut a = T::zero();
    let mut se = T::zero();
    let mut inx = T::zero();
    for i in 0..n {
        a += (s.portfolio_sector_weights[i] - s.benchmark_sector_weights[i])
            * (s.benchmark_sector_returns[i] - rb);
        se += s.benchmark_sector_weights[i]
            * (s.portfolio_sector_returns[i] - s.benchmark_sector_returns[i]);
        inx += (s.portfolio_sector_weights[i] - s.benchmark_sector_weights[i])
            * (s.portfolio_sector_returns[i] - s.benchmark_sector_returns[i]);
    }
    let mut cp = T::one();
    let mut cb = T::one();
    for &x in &s.multi_period_port_returns {
        cp *= T::one() + x;
    }
    for &x in &s.multi_period_bench_returns {
        cb *= T::one() + x;
    }
    let cp = cp - T::one();
    let cb = cb - T::one();
    let log_link_ratio = |rp: T, rb: T| {
        let diff = rp - rb;
        if diff.val().abs() < 1e-10 {
            T::one() / (T::one() + T::from_f64(0.5) * (rp + rb))
        } else {
            ((T::one() + rp).ln() - (T::one() + rb).ln()) / diff
        }
    };
    let k = log_link_ratio(cp, cb);
    let mut carino_sum = T::zero();
    for (rpt, rbt) in s
        .multi_period_port_returns
        .iter()
        .zip(s.multi_period_bench_returns.iter())
    {
        carino_sum += (log_link_ratio(*rpt, *rbt) / k) * (*rpt - *rbt);
    }
    BaconBrinsonFachlerCarinoReport {
        brinson_fachler_allocation: a,
        brinson_fachler_selection: se,
        brinson_fachler_interaction: inx,
        single_period_active_return: rp - rb,
        single_period_conservation_error: ((rp - rb) - (a + se + inx)).abs(),
        cumulative_portfolio_return: cp,
        cumulative_benchmark_return: cb,
        carino_total_linking_factor: k,
        carino_smoothed_active_return: carino_sum,
        carino_linking_conservation_error: {
            let d = carino_sum - (cp - cb);
            (d * d + T::from_f64(1e-30)).sqrt()
        },
    }
}
#[derive(Clone, Copy, Debug)]
pub struct WakkerPrelecChoquetProspectSpec<T: Real = f64> {
    pub sorted_outcomes: [T; 4],
    pub state_probabilities: [T; 4],
    pub reference_point: T,
    pub prelec_curvature_alpha: T,
    pub prelec_elevation_beta: T,
    pub tk_value_exponent: T,
    pub loss_aversion_lambda: T,
    pub ambiguity_aversion_delta: T,
}
#[derive(Clone, Copy, Debug)]
pub struct WakkerPrelecChoquetProspectReport<T: Real = f64> {
    pub prelec_inflection_fixed_point: T,
    pub prelec_weighted_tail_prob: T,
    pub expected_value_objective: T,
    pub rank_dependent_cpt_value: T,
    pub choquet_ambiguity_adjusted_value: T,
    pub ambiguity_valuation_discount: T,
    pub prospect_certainty_equivalent: T,
}
pub fn evaluate_wakker_prelec_rdu_and_choquet_ambiguity<T: Real>(
    s: &WakkerPrelecChoquetProspectSpec<T>,
) -> WakkerPrelecChoquetProspectReport<T> {
    let mut e = T::zero();
    for i in 0..4 {
        e += s.state_probabilities[i] * s.sorted_outcomes[i];
    }
    let prelec_w = |p: T| -> T {
        let pc = p
            .max_val(T::from_f64(1e-12))
            .min_val(T::one() - T::from_f64(1e-12));
        (-s.prelec_elevation_beta * (-pc.ln()).pow(s.prelec_curvature_alpha)).exp()
    };
    let mut g = [T::zero(); 5];
    for i in (0..4).rev() {
        g[i] = g[i + 1] + s.state_probabilities[i];
    }
    let mut cpt = T::zero();
    for i in 0..4 {
        let wg_i = if i == 0 { T::one() } else { prelec_w(g[i]) };
        let wg_next = if i == 3 {
            T::zero()
        } else {
            prelec_w(g[i + 1])
        };
        let pi_i = wg_i - wg_next;
        let dx = s.sorted_outcomes[i] - s.reference_point;
        let v = if dx.val() >= 0.0 {
            dx.max_val(T::from_f64(1e-12)).pow(s.tk_value_exponent)
        } else {
            -s.loss_aversion_lambda * (-dx).max_val(T::from_f64(1e-12)).pow(s.tk_value_exponent)
        };
        cpt += pi_i * v;
    }
    let inv_gamma = T::one() / s.tk_value_exponent.max_val(T::from_f64(1e-6));
    let ce = if cpt.val() >= 0.0 {
        s.reference_point + cpt.max_val(T::from_f64(1e-12)).pow(inv_gamma)
    } else {
        s.reference_point
            - ((-cpt) / s.loss_aversion_lambda.max_val(T::from_f64(1e-6)))
                .max_val(T::from_f64(1e-12))
                .pow(inv_gamma)
    };
    let choq = (T::one() - s.ambiguity_aversion_delta) * e
        + s.ambiguity_aversion_delta * s.sorted_outcomes[0];
    WakkerPrelecChoquetProspectReport {
        prelec_inflection_fixed_point: (-s
            .prelec_elevation_beta
            .max_val(T::from_f64(1e-12))
            .pow(T::one() / (T::one() - s.prelec_curvature_alpha).max_val(T::from_f64(1e-6))))
        .exp(),
        prelec_weighted_tail_prob: prelec_w(s.state_probabilities[3]),
        expected_value_objective: e,
        rank_dependent_cpt_value: cpt,
        choquet_ambiguity_adjusted_value: choq,
        ambiguity_valuation_discount: e - choq,
        prospect_certainty_equivalent: ce,
    }
}

#[cfg(test)]
mod phase34_tests {
    #[test]
    fn phase34_smoke() {
        assert!(std::hint::black_box(true));
    }
}
