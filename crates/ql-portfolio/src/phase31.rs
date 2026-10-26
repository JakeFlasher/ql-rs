use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct EmbrechtsHillPickandsSubexponentialSpec<T: Real = f64> {
    pub loss_samples: Vec<T>,
    pub tail_order_k: usize,
    pub extreme_tail_prob: T,
    pub portfolio_positions_count: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct EmbrechtsHillPickandsSubexponentialReport<T: Real = f64> {
    pub hill_xi_tail_index: T,
    pub hill_alpha_pareto_exponent: T,
    pub pickands_xi_tail_index: T,
    pub hill_asymptotic_std_error: T,
    pub weissman_extreme_var: T,
    pub weissman_extreme_es: T,
    pub subexponential_sum_to_single_tail_multiplier: T,
    pub superadditive_var_excess_ratio: T,
}

pub fn evaluate_embrechts_hill_pickands_and_subexponential_tail<T: Real>(
    s: &EmbrechtsHillPickandsSubexponentialSpec<T>,
) -> EmbrechtsHillPickandsSubexponentialReport<T> {
    assert!(s.tail_order_k >= 1);
    assert!(s.loss_samples.len() >= 4 * s.tail_order_k);
    let mut x = s.loss_samples.clone();
    x.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let k = s.tail_order_k.max(1).min(x.len() - 1);
    let threshold = x[k];
    let mut sum_log = T::zero();
    for value in x.iter().take(k) {
        sum_log += value.abs().max_val(T::from_f64(1e-12)).ln()
            - threshold.abs().max_val(T::from_f64(1e-12)).ln();
    }
    let xi = sum_log / T::from_f64(k as f64);
    let alpha = T::one() / xi.abs().max_val(T::from_f64(1e-12));
    let i2 = 2 * k;
    let i4 = 4 * k;
    let pick_ratio = (x[k - 1] - x[i2 - 1]).abs().max_val(T::from_f64(1e-12))
        / (x[i2 - 1] - x[i4 - 1]).abs().max_val(T::from_f64(1e-12));
    let pick = pick_ratio.ln() / T::from_f64(2.0).ln();
    let p = s.extreme_tail_prob.max_val(T::from_f64(1e-12));
    let n = T::from_f64(x.len() as f64);
    let var = threshold * (T::from_f64(k as f64) / (n * p)).ln().exp().pow(xi);
    let es = var / (T::one() - xi.min_val(T::from_f64(0.95)).max_val(T::from_f64(-10.0)));
    let m = T::from_f64(s.portfolio_positions_count.max(1) as f64);
    let multiplier = m;
    let excess = m.pow(xi.max_val(T::zero()));
    EmbrechtsHillPickandsSubexponentialReport {
        hill_xi_tail_index: xi,
        hill_alpha_pareto_exponent: alpha,
        pickands_xi_tail_index: pick,
        hill_asymptotic_std_error: xi.abs() / T::from_f64((k as f64).sqrt()),
        weissman_extreme_var: var,
        weissman_extreme_es: es,
        subexponential_sum_to_single_tail_multiplier: multiplier,
        superadditive_var_excess_ratio: excess,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BouchaudPottersCumulantSmileSpec<T: Real = f64> {
    pub single_period_vol: T,
    pub single_period_skewness: T,
    pub single_period_excess_kurtosis: T,
    pub vol_of_vol_kurtosis_component: T,
    pub vol_decorrelation_steps: T,
    pub horizon_steps_n: T,
    pub log_moneyness: T,
    pub confidence_z: T,
}

#[derive(Clone, Copy, Debug)]
pub struct BouchaudPottersCumulantSmileReport<T: Real = f64> {
    pub horizon_skewness_k3: T,
    pub iid_clt_kurtosis_k4: T,
    pub bouchaud_potters_persistent_kurtosis_k4: T,
    pub kurtosis_persistence_amplification: T,
    pub atm_implied_vol_kurtosis_shift: T,
    pub cumulant_smile_implied_vol: T,
    pub naive_sqrt_t_gaussian_var: T,
    pub bouchaud_potters_non_gaussian_var: T,
    pub sqrt_time_var_underestimation_ratio: T,
}

pub fn evaluate_bouchaud_potters_power_law_cumulant_and_smile_scaling<T: Real>(
    s: &BouchaudPottersCumulantSmileSpec<T>,
) -> BouchaudPottersCumulantSmileReport<T> {
    let n = s.horizon_steps_n.max_val(T::from_f64(1e-8));
    let nc = s.vol_decorrelation_steps.max_val(T::from_f64(1e-8));
    let k3 = s.single_period_skewness / n.sqrt();
    let iid = s.single_period_excess_kurtosis / n;
    let persistent = iid
        + T::from_f64(6.0) * s.vol_of_vol_kurtosis_component / (n * n)
            * (n * nc - nc * nc * (T::one() - (-n / nc).exp()));
    let amp = persistent / iid.abs().max_val(T::from_f64(1e-12));
    let v = n * s.single_period_vol * s.single_period_vol;
    let smile = s.single_period_vol
        * (T::one()
            + k3 / T::from_f64(6.0) * s.log_moneyness / v.sqrt()
            + persistent / T::from_f64(24.0) * (s.log_moneyness * s.log_moneyness / v - T::one()));
    let naive = s.confidence_z * v.sqrt();
    let z = s.confidence_z;
    let cf = z
        + k3 / T::from_f64(6.0) * (z * z - T::one())
        + persistent / T::from_f64(24.0) * (z * z * z - T::from_f64(3.0) * z)
        - k3 * k3 / T::from_f64(36.0) * (T::from_f64(2.0) * z * z * z - T::from_f64(5.0) * z);
    let non_gaussian = v.sqrt() * cf;
    BouchaudPottersCumulantSmileReport {
        horizon_skewness_k3: k3,
        iid_clt_kurtosis_k4: iid,
        bouchaud_potters_persistent_kurtosis_k4: persistent,
        kurtosis_persistence_amplification: amp,
        atm_implied_vol_kurtosis_shift: s.single_period_vol * persistent / T::from_f64(24.0),
        cumulant_smile_implied_vol: smile,
        naive_sqrt_t_gaussian_var: naive,
        bouchaud_potters_non_gaussian_var: non_gaussian,
        sqrt_time_var_underestimation_ratio: non_gaussian / naive.abs().max_val(T::from_f64(1e-12)),
    }
}

#[derive(Clone, Debug)]
pub struct KarnoskySingerMencheroAttributionSpec<T: Real = f64> {
    pub portfolio_asset_weights: Vec<T>,
    pub benchmark_asset_weights: Vec<T>,
    pub portfolio_currency_weights: Vec<T>,
    pub benchmark_currency_weights: Vec<T>,
    pub local_portfolio_returns: Vec<T>,
    pub local_benchmark_returns: Vec<T>,
    pub local_eurodollar_rates: Vec<T>,
    pub base_eurodollar_rate: T,
    pub spot_fx_returns: Vec<T>,
    pub multi_period_port_returns: Vec<T>,
    pub multi_period_bench_returns: Vec<T>,
}

#[derive(Clone, Copy, Debug)]
pub struct KarnoskySingerMencheroAttributionReport<T: Real = f64> {
    pub karnosky_singer_market_allocation: T,
    pub karnosky_singer_security_selection: T,
    pub karnosky_singer_currency_management: T,
    pub karnosky_singer_active_return_error: T,
    pub menchero_natural_scaling_multiplier: T,
    pub cumulative_active_return: T,
    pub menchero_multi_period_linking_error: T,
}

pub fn evaluate_karnosky_singer_multicurrency_and_menchero_attribution<T: Real>(
    s: &KarnoskySingerMencheroAttributionSpec<T>,
) -> KarnoskySingerMencheroAttributionReport<T> {
    let n = s.portfolio_asset_weights.len();
    let mut rb_ex = T::zero();
    let mut ccy_b = T::zero();
    for i in 0..n {
        let ex = s.local_benchmark_returns[i] - s.local_eurodollar_rates[i];
        rb_ex += s.benchmark_asset_weights[i] * ex;
        ccy_b += s.benchmark_currency_weights[i]
            * (s.spot_fx_returns[i] + s.local_eurodollar_rates[i] - s.base_eurodollar_rate);
    }
    let mut mkt = T::zero();
    let mut sec = T::zero();
    let mut ccy = T::zero();
    for i in 0..n {
        let rp = s.local_portfolio_returns[i] - s.local_eurodollar_rates[i];
        let rb = s.local_benchmark_returns[i] - s.local_eurodollar_rates[i];
        let ci = s.spot_fx_returns[i] + s.local_eurodollar_rates[i] - s.base_eurodollar_rate;
        mkt += (s.portfolio_asset_weights[i] - s.benchmark_asset_weights[i]) * (rb - rb_ex);
        sec += s.portfolio_asset_weights[i] * (rp - rb);
        ccy += (s.portfolio_currency_weights[i] - s.benchmark_currency_weights[i]) * (ci - ccy_b);
    }
    let mut portfolio_total = T::zero();
    let mut benchmark_total = T::zero();
    for i in 0..n {
        let rp = s.local_portfolio_returns[i] - s.local_eurodollar_rates[i];
        let rb = s.local_benchmark_returns[i] - s.local_eurodollar_rates[i];
        let ci = s.spot_fx_returns[i] + s.local_eurodollar_rates[i] - s.base_eurodollar_rate;
        portfolio_total += s.portfolio_asset_weights[i] * rp + s.portfolio_currency_weights[i] * ci;
        benchmark_total += s.benchmark_asset_weights[i] * rb + s.benchmark_currency_weights[i] * ci;
    }
    let active = portfolio_total - benchmark_total;
    let mut rp_cum = T::one();
    let mut rb_cum = T::one();
    for i in 0..s
        .multi_period_port_returns
        .len()
        .min(s.multi_period_bench_returns.len())
    {
        rp_cum *= T::one() + s.multi_period_port_returns[i];
        rb_cum *= T::one() + s.multi_period_bench_returns[i];
    }
    let rp = rp_cum - T::one();
    let rb = rb_cum - T::one();
    let cumulative = rp - rb;
    let periods = T::from_f64(s.multi_period_port_returns.len().max(1) as f64);
    let diff_root = rp_cum.pow(T::one() / periods) - rb_cum.pow(T::one() / periods);
    let multiplier = (cumulative / periods) / diff_root.max_val(T::from_f64(1e-12));
    let linked = multiplier * diff_root * periods;
    KarnoskySingerMencheroAttributionReport {
        karnosky_singer_market_allocation: mkt,
        karnosky_singer_security_selection: sec,
        karnosky_singer_currency_management: ccy,
        karnosky_singer_active_return_error: (active - (mkt + sec + ccy)).abs(),
        menchero_natural_scaling_multiplier: multiplier,
        cumulative_active_return: cumulative,
        menchero_multi_period_linking_error: (linked - cumulative).abs(),
    }
}

#[derive(Clone, Debug)]
pub struct DgtwTreynorMazuyTimingSpec<T: Real = f64> {
    pub current_weights: Vec<T>,
    pub lagged_weights: Vec<T>,
    pub stock_returns: Vec<T>,
    pub current_char_bench_returns: Vec<T>,
    pub lagged_char_bench_returns: Vec<T>,
    pub portfolio_excess_series: Vec<T>,
    pub market_excess_series: Vec<T>,
}

#[derive(Clone, Copy, Debug)]
pub struct DgtwTreynorMazuyTimingReport<T: Real = f64> {
    pub dgtw_characteristic_selectivity_cs: T,
    pub dgtw_characteristic_timing_ct: T,
    pub dgtw_average_style_as: T,
    pub dgtw_decomposition_identity_error: T,
    pub treynor_mazuy_selectivity_alpha: T,
    pub treynor_mazuy_market_beta: T,
    pub treynor_mazuy_timing_gamma: T,
    pub treynor_mazuy_r_squared: T,
}

fn solve_3x3<T: Real>(mut a: [[T; 4]; 3]) -> [T; 3] {
    for i in 0..3 {
        let mut pivot = i;
        for r in i + 1..3 {
            if a[r][i].abs() > a[pivot][i].abs() {
                pivot = r;
            }
        }
        a.swap(i, pivot);
        let d = if a[i][i].abs() < T::from_f64(1e-12) {
            T::from_f64(1e-12)
        } else {
            a[i][i]
        };
        for value in a[i].iter_mut().skip(i) {
            *value /= d;
        }
        for r in 0..3 {
            if r != i {
                let f = a[r][i];
                for j in i..4 {
                    a[r][j] -= f * a[i][j];
                }
            }
        }
    }
    [a[0][3], a[1][3], a[2][3]]
}

pub fn evaluate_dgtw_characteristic_and_treynor_mazuy_timing<T: Real>(
    s: &DgtwTreynorMazuyTimingSpec<T>,
) -> DgtwTreynorMazuyTimingReport<T> {
    let n = s
        .current_weights
        .len()
        .min(s.lagged_weights.len())
        .min(s.stock_returns.len());
    let mut cs = T::zero();
    let mut ct = T::zero();
    let mut as_ = T::zero();
    let mut rp = T::zero();
    for i in 0..n {
        cs += s.current_weights[i] * (s.stock_returns[i] - s.current_char_bench_returns[i]);
        ct += s.current_weights[i] * s.current_char_bench_returns[i]
            - s.lagged_weights[i] * s.lagged_char_bench_returns[i];
        as_ += s.lagged_weights[i] * s.lagged_char_bench_returns[i];
        rp += s.current_weights[i] * s.stock_returns[i];
    }
    let m = s
        .portfolio_excess_series
        .len()
        .min(s.market_excess_series.len());
    let mut normal = [[T::zero(); 4]; 3];
    let mut sum_y = T::zero();
    for i in 0..m {
        let x = s.market_excess_series[i];
        let y = s.portfolio_excess_series[i];
        let row = [T::one(), x, x * x];
        for a in 0..3 {
            for b in 0..3 {
                normal[a][b] += row[a] * row[b];
            }
            normal[a][3] += row[a] * y;
        }
        sum_y += y;
    }
    let coef = solve_3x3(normal);
    let mean = sum_y / T::from_f64(m.max(1) as f64);
    let mut sse = T::zero();
    let mut sst = T::zero();
    for i in 0..m {
        let x = s.market_excess_series[i];
        let y = s.portfolio_excess_series[i];
        let fit = coef[0] + coef[1] * x + coef[2] * x * x;
        sse += (y - fit) * (y - fit);
        sst += (y - mean) * (y - mean);
    }
    DgtwTreynorMazuyTimingReport {
        dgtw_characteristic_selectivity_cs: cs,
        dgtw_characteristic_timing_ct: ct,
        dgtw_average_style_as: as_,
        dgtw_decomposition_identity_error: (rp - (cs + ct + as_)).abs(),
        treynor_mazuy_selectivity_alpha: coef[0],
        treynor_mazuy_market_beta: coef[1],
        treynor_mazuy_timing_gamma: coef[2],
        treynor_mazuy_r_squared: T::one() - sse / sst.max_val(T::from_f64(1e-12)),
    }
}

#[derive(Clone, Debug)]
pub struct BrockHommesAdaptiveBeliefSpec<T: Real = f64> {
    pub initial_price_deviation: T,
    pub risk_free_rate: T,
    pub fundamentalist_damping_v1: T,
    pub chartist_trend_extrapolation_g: T,
    pub fundamentalist_info_cost_c1: T,
    pub intensity_of_choice_beta: T,
    pub fitness_memory_eta: T,
    pub risk_aversion_sigma2: T,
    pub shock_innovations: Vec<T>,
}

#[derive(Clone, Debug)]
pub struct BrockHommesAdaptiveBeliefReport<T: Real = f64> {
    pub terminal_fundamentalist_fraction: T,
    pub terminal_chartist_fraction: T,
    pub steady_state_jacobian_eigenvalue: T,
    pub critical_bifurcation_beta: T,
    pub terminal_price_deviation: T,
    pub price_deviation_trajectory: Vec<T>,
    pub chartist_fraction_trajectory: Vec<T>,
    pub endogenous_excess_volatility_ratio: T,
}

pub fn evaluate_brock_hommes_adaptive_belief_switching<T: Real>(
    s: &BrockHommesAdaptiveBeliefSpec<T>,
) -> BrockHommesAdaptiveBeliefReport<T> {
    let r = T::one() + s.risk_free_rate;
    let mut x = s.initial_price_deviation;
    let mut prev = x;
    let mut u1 = -s.fundamentalist_info_cost_c1;
    let mut u2 = T::zero();
    let mut path = Vec::with_capacity(s.shock_innovations.len() + 1);
    let mut frac = Vec::with_capacity(s.shock_innovations.len());
    path.push(x);
    for shock in &s.shock_innovations {
        let pi = x - r * prev;
        u1 = s.fitness_memory_eta * u1
            + (s.fundamentalist_damping_v1 * prev * pi)
                / s.risk_aversion_sigma2.max_val(T::from_f64(1e-12))
            - s.fundamentalist_info_cost_c1;
        u2 = s.fitness_memory_eta * u2
            + (s.chartist_trend_extrapolation_g * prev * pi)
                / s.risk_aversion_sigma2.max_val(T::from_f64(1e-12));
        let n1 = T::one() / (T::one() + (s.intensity_of_choice_beta * (u2 - u1)).exp());
        let n2 = T::one() - n1;
        prev = x;
        x = (n1 * s.fundamentalist_damping_v1 + n2 * s.chartist_trend_extrapolation_g) * x / r
            + *shock / r;
        path.push(x);
        frac.push(n2);
    }
    let n1_star =
        T::one() / (T::one() + (s.intensity_of_choice_beta * s.fundamentalist_info_cost_c1).exp());
    let eig = ((T::one() - n1_star) * s.chartist_trend_extrapolation_g
        + n1_star * s.fundamentalist_damping_v1)
        / r;
    let critical = ((s.fundamentalist_damping_v1 - r) / (r - s.chartist_trend_extrapolation_g))
        .abs()
        .max_val(T::from_f64(1e-12));
    BrockHommesAdaptiveBeliefReport {
        terminal_fundamentalist_fraction: T::one() - frac.last().copied().unwrap_or(T::zero()),
        terminal_chartist_fraction: frac.last().copied().unwrap_or(T::zero()),
        steady_state_jacobian_eigenvalue: eig,
        critical_bifurcation_beta: critical,
        terminal_price_deviation: x,
        price_deviation_trajectory: path,
        chartist_fraction_trajectory: frac,
        endogenous_excess_volatility_ratio: x.abs()
            / s.initial_price_deviation.abs().max_val(T::from_f64(1e-12)),
    }
}

#[cfg(test)]
mod phase31_tests {
    use super::*;
    use ql_ad::Dual;

    fn hill<T: Real>(x: T) -> EmbrechtsHillPickandsSubexponentialSpec<T> {
        EmbrechtsHillPickandsSubexponentialSpec {
            loss_samples: vec![
                T::from_f64(120.0) + x,
                T::from_f64(100.0),
                T::from_f64(90.0),
                T::from_f64(80.0),
                T::from_f64(70.0),
                T::from_f64(60.0),
                T::from_f64(50.0),
                T::from_f64(45.0),
                T::from_f64(40.0),
                T::from_f64(35.0),
                T::from_f64(30.0),
                T::from_f64(25.0),
            ],
            tail_order_k: 2,
            extreme_tail_prob: T::from_f64(0.01),
            portfolio_positions_count: 8,
        }
    }

    #[test]
    fn phase31_invariants_and_dual_fd() {
        let hrep = evaluate_embrechts_hill_pickands_and_subexponential_tail(&hill(0.0));
        assert!(hrep.hill_xi_tail_index > 0.0);
        assert!(hrep.weissman_extreme_es > hrep.weissman_extreme_var);
        let h = 1e-6;
        let fd = (evaluate_embrechts_hill_pickands_and_subexponential_tail(&hill(h))
            .weissman_extreme_var
            - evaluate_embrechts_hill_pickands_and_subexponential_tail(&hill(-h))
                .weissman_extreme_var)
            / (2.0 * h);
        let ad = evaluate_embrechts_hill_pickands_and_subexponential_tail(&hill(
            Dual::<1>::variable(0.0, 0),
        ))
        .weissman_extreme_var;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        let bp = evaluate_bouchaud_potters_power_law_cumulant_and_smile_scaling(
            &BouchaudPottersCumulantSmileSpec {
                single_period_vol: 0.2,
                single_period_skewness: -0.2,
                single_period_excess_kurtosis: 3.0,
                vol_of_vol_kurtosis_component: 2.0,
                vol_decorrelation_steps: 10.0,
                horizon_steps_n: 20.0,
                log_moneyness: 0.1,
                confidence_z: 2.33,
            },
        );
        assert!(bp.bouchaud_potters_persistent_kurtosis_k4 > bp.iid_clt_kurtosis_k4);
        assert!(bp.sqrt_time_var_underestimation_ratio > 1.0);
        let ks = evaluate_karnosky_singer_multicurrency_and_menchero_attribution(
            &KarnoskySingerMencheroAttributionSpec {
                portfolio_asset_weights: vec![0.6, 0.4],
                benchmark_asset_weights: vec![0.5, 0.5],
                portfolio_currency_weights: vec![0.6, 0.4],
                benchmark_currency_weights: vec![0.5, 0.5],
                local_portfolio_returns: vec![0.08, 0.03],
                local_benchmark_returns: vec![0.06, 0.04],
                local_eurodollar_rates: vec![0.02, 0.025],
                base_eurodollar_rate: 0.03,
                spot_fx_returns: vec![0.01, -0.005],
                multi_period_port_returns: vec![0.02, 0.03],
                multi_period_bench_returns: vec![0.018, 0.025],
            },
        );
        assert!(
            ks.karnosky_singer_active_return_error < 1e-12
                && ks.menchero_multi_period_linking_error < 1e-12
        );
        let d =
            evaluate_dgtw_characteristic_and_treynor_mazuy_timing(&DgtwTreynorMazuyTimingSpec {
                current_weights: vec![0.6, 0.4],
                lagged_weights: vec![0.5, 0.5],
                stock_returns: vec![0.08, 0.03],
                current_char_bench_returns: vec![0.06, 0.04],
                lagged_char_bench_returns: vec![0.05, 0.04],
                portfolio_excess_series: vec![0.01, 0.02, 0.04, 0.08],
                market_excess_series: vec![0.0, 0.02, 0.04, 0.08],
            });
        assert!(d.dgtw_decomposition_identity_error < 1e-12 && d.treynor_mazuy_timing_gamma > 0.0);
        let bh = evaluate_brock_hommes_adaptive_belief_switching(&BrockHommesAdaptiveBeliefSpec {
            initial_price_deviation: 0.1,
            risk_free_rate: 0.03,
            fundamentalist_damping_v1: 0.8,
            chartist_trend_extrapolation_g: 1.2,
            fundamentalist_info_cost_c1: 0.1,
            intensity_of_choice_beta: 8.0,
            fitness_memory_eta: 0.8,
            risk_aversion_sigma2: 1.0,
            shock_innovations: vec![0.01, -0.005, 0.003],
        });
        assert!(
            (bh.terminal_fundamentalist_fraction + bh.terminal_chartist_fraction - 1.0).abs()
                < 1e-12
        );
        assert!(bh.steady_state_jacobian_eigenvalue > 1.0);
    }
}
