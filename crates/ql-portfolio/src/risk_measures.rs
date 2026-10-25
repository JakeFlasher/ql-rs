//! Coherent empirical and parametric risk measures plus Euler allocations.

use ql_ad::{Real, TapeScope, TapeVar};

#[derive(Clone, Debug, PartialEq)]
pub struct EulerRiskAllocation {
    pub es: f64,
    pub marginal_es: Vec<f64>,
    pub component_es: Vec<f64>,
}

pub type EulerEsAllocation = EulerRiskAllocation;

fn validate_confidence(confidence: f64) -> f64 {
    assert!(confidence > 0.5 && confidence < 1.0);
    1.0 - confidence
}

/// Population semi-variance conditional on observations below their mean,
/// with the same N/(N-1) correction used by QuantLib risk statistics.
pub fn semi_variance(returns: &[f64]) -> f64 {
    if returns.is_empty() {
        return 0.0;
    }
    let mean = returns.iter().sum::<f64>() / returns.len() as f64;
    regret(returns, mean)
}

pub fn semi_deviation(returns: &[f64]) -> f64 {
    semi_variance(returns).sqrt()
}

pub fn downside_variance(returns: &[f64]) -> f64 {
    regret(returns, 0.0)
}

pub fn downside_deviation(returns: &[f64]) -> f64 {
    downside_variance(returns).sqrt()
}

pub fn regret(returns: &[f64], target: f64) -> f64 {
    let below: Vec<f64> = returns
        .iter()
        .copied()
        .filter(|x| *x < target)
        .map(|x| (x - target).powi(2))
        .collect();
    if below.len() <= 1 {
        return 0.0;
    }
    below.iter().sum::<f64>() / (below.len() as f64 - 1.0)
}

/// Acerbi--Tasche fractional-tail empirical VaR, with returns interpreted as
/// PnL (negative observations are losses).
pub fn empirical_var(returns: &[f64], confidence: f64) -> f64 {
    let alpha = validate_confidence(confidence);
    assert!(!returns.is_empty());
    let mut sorted = returns.to_vec();
    sorted.sort_by(f64::total_cmp);
    let index = ((sorted.len() as f64 * alpha).floor() as usize).min(sorted.len() - 1);
    -sorted[index]
}

/// Coherent finite-sample expected shortfall using the exact fractional-tail
/// estimator from Acerbi--Tasche (2002), Eq. (13).
pub fn empirical_expected_shortfall(returns: &[f64], confidence: f64) -> f64 {
    let alpha = validate_confidence(confidence);
    assert!(!returns.is_empty());
    let mut sorted = returns.to_vec();
    sorted.sort_by(f64::total_cmp);
    let n = sorted.len();
    let mass = n as f64 * alpha;
    let whole = mass.floor() as usize;
    let theta = mass - whole as f64;
    let mut numerator = sorted[..whole.min(n)].iter().sum::<f64>();
    if whole < n {
        numerator += theta * sorted[whole];
    }
    -numerator / mass
}

pub fn value_at_risk(returns: &[f64], confidence: f64) -> f64 {
    empirical_var(returns, confidence)
}

pub fn expected_shortfall(returns: &[f64], confidence: f64) -> f64 {
    empirical_expected_shortfall(returns, confidence)
}

pub fn gaussian_var<T: Real>(mean_pnl: T, std_pnl: T, confidence: f64) -> T {
    assert!(confidence > 0.5 && confidence < 1.0);
    mean_pnl.neg() + std_pnl * ql_math::inv_norm_cdf(confidence)
}

pub fn gaussian_expected_shortfall<T: Real>(mean_pnl: T, std_pnl: T, confidence: f64) -> T {
    assert!(confidence > 0.5 && confidence < 1.0);
    let z = ql_math::inv_norm_cdf(confidence);
    let phi = std::f64::consts::TAU.sqrt().recip() * (-0.5 * z * z).exp();
    mean_pnl.neg() + std_pnl * (phi / (1.0 - confidence))
}

pub fn cornish_fisher_var<T: Real>(
    mean_pnl: T,
    std_pnl: T,
    skew: T,
    excess_kurtosis: T,
    confidence: f64,
) -> T {
    assert!(confidence > 0.5 && confidence < 1.0);
    let z = T::from_f64(ql_math::inv_norm_cdf(confidence));
    let z2 = z * z;
    let z3 = z2 * z;
    let adjusted = z - (z2 - 1.0) * skew / 6.0 + (z3 - z * 3.0) * excess_kurtosis / 24.0
        - (z3 * 2.0 - z * 5.0) * skew * skew / 36.0;
    mean_pnl.neg() + std_pnl * adjusted
}

/// Euler allocation of Gaussian ES through one reverse tape pass.
pub fn euler_gaussian_es_allocation_aad(
    means: &[f64],
    covariance: &[Vec<f64>],
    weights: &[f64],
    confidence: f64,
) -> EulerRiskAllocation {
    assert_eq!(means.len(), weights.len());
    assert_eq!(covariance.len(), weights.len());
    assert!(covariance.iter().all(|row| row.len() == weights.len()));
    let guard = TapeScope::guard_with_capacity(weights.len() * weights.len() * 4 + 32);
    let vars: Vec<TapeVar<f64>> = weights
        .iter()
        .map(|&w| TapeScope::with_active(|tape| tape.input(w)))
        .collect();
    let mut mean = TapeVar::zero();
    for (x, &m) in vars.iter().zip(means) {
        mean += *x * m;
    }
    let mut variance = TapeVar::zero();
    for i in 0..vars.len() {
        for j in 0..vars.len() {
            variance += vars[i] * vars[j] * covariance[i][j];
        }
    }
    let es =
        gaussian_expected_shortfall(mean, variance.max_val(TapeVar::zero()).sqrt(), confidence);
    let adjoints = TapeScope::gradient_from(es);
    let marginal: Vec<f64> = vars.iter().map(|v| adjoints.gradient(*v)).collect();
    let component: Vec<f64> = weights
        .iter()
        .zip(&marginal)
        .map(|(&w, &m)| w * m)
        .collect();
    let result = EulerRiskAllocation {
        es: es.val,
        marginal_es: marginal,
        component_es: component,
    };
    let _ = guard.stop();
    result
}

/// Exact scenario Euler allocation using the fractional Acerbi--Tasche tail
/// weights. Rows are asset returns and `weights` are portfolio exposures.
pub fn euler_scenario_es_allocation(
    scenarios: &[Vec<f64>],
    weights: &[f64],
    confidence: f64,
) -> EulerRiskAllocation {
    let alpha = validate_confidence(confidence);
    assert!(!scenarios.is_empty());
    assert!(scenarios.iter().all(|row| row.len() == weights.len()));
    let mut pnl: Vec<(f64, usize)> = scenarios
        .iter()
        .enumerate()
        .map(|(i, row)| (row.iter().zip(weights).map(|(r, w)| r * w).sum(), i))
        .collect();
    pnl.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mass = scenarios.len() as f64 * alpha;
    let whole = mass.floor() as usize;
    let theta = mass - whole as f64;
    let mut marginal = vec![0.0; weights.len()];
    for (_, index) in pnl.iter().take(whole) {
        for (j, value) in scenarios[*index].iter().enumerate() {
            marginal[j] -= *value / mass;
        }
    }
    if whole < pnl.len() {
        let (_, index) = pnl[whole];
        for (j, value) in scenarios[index].iter().enumerate() {
            marginal[j] -= theta * *value / mass;
        }
    }
    let es = -pnl.iter().take(whole).map(|(x, _)| x).sum::<f64>() / mass
        - if whole < pnl.len() {
            theta * pnl[whole].0 / mass
        } else {
            0.0
        };
    let component: Vec<f64> = weights
        .iter()
        .zip(&marginal)
        .map(|(&w, &m)| w * m)
        .collect();
    EulerRiskAllocation {
        es,
        marginal_es: marginal,
        component_es: component,
    }
}

pub struct EmpiricalRiskStatistics;
impl EmpiricalRiskStatistics {
    pub fn semi_variance(x: &[f64]) -> f64 {
        semi_variance(x)
    }
    pub fn semi_deviation(x: &[f64]) -> f64 {
        semi_deviation(x)
    }
    pub fn downside_variance(x: &[f64]) -> f64 {
        downside_variance(x)
    }
    pub fn downside_deviation(x: &[f64]) -> f64 {
        downside_deviation(x)
    }
    pub fn regret(x: &[f64], target: f64) -> f64 {
        regret(x, target)
    }
    pub fn empirical_var(x: &[f64], confidence: f64) -> f64 {
        empirical_var(x, confidence)
    }
    pub fn empirical_expected_shortfall(x: &[f64], confidence: f64) -> f64 {
        empirical_expected_shortfall(x, confidence)
    }
    pub fn value_at_risk(x: &[f64], confidence: f64) -> f64 {
        empirical_var(x, confidence)
    }
    pub fn expected_shortfall(x: &[f64], confidence: f64) -> f64 {
        empirical_expected_shortfall(x, confidence)
    }
}

pub fn parametric_gaussian_evar<T: Real>(mean: T, std_dev: T, confidence: f64) -> T {
    let alpha = validate_confidence(confidence);
    -mean + std_dev * (-2.0 * alpha.ln()).sqrt()
}

pub fn entropic_value_at_risk<T: Real>(returns: &[T], confidence: f64) -> T {
    assert!(!returns.is_empty());
    let alpha = validate_confidence(confidence);
    let losses: Vec<T> = returns.iter().map(|&r| -r).collect();
    let max_idx = losses
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.val().total_cmp(&b.val()))
        .map(|(i, _)| i)
        .unwrap_or(0);
    let l_max_t = losses[max_idx];
    let l_max = l_max_t.val();
    let min_idx = losses
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.val().total_cmp(&b.val()))
        .map(|(i, _)| i)
        .unwrap_or(0);
    let l_min_t = losses[min_idx];
    let l_min = l_min_t.val();
    if (l_max - l_min).abs() < 1e-14 {
        return losses[0];
    }
    let n_f = losses.len() as f64;
    let m_max = losses
        .iter()
        .filter(|l| (l.val() - l_max).abs() <= 1e-14)
        .count() as f64;
    if alpha <= m_max / n_f {
        return l_max_t;
    }
    let obj_f64 = |z: f64| -> f64 {
        let sum_e: f64 = losses.iter().map(|l| (z * (l.val() - l_max)).exp()).sum();
        l_max + ((sum_e / n_f).ln() - alpha.ln()) / z
    };
    let mut a = 1e-4f64;
    let mut b = 500.0 / (l_max - l_min).max(1e-4);
    let gr = 0.5 * (5.0f64.sqrt() - 1.0);
    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);
    let mut fc = obj_f64(c);
    let mut fd = obj_f64(d);
    for _ in 0..64 {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = obj_f64(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = obj_f64(d);
        }
    }
    let z_star = T::from_f64(0.5 * (a + b));
    let mut sum_exp = T::zero();
    for &l in &losses {
        sum_exp += (z_star * (l - l_max_t)).exp();
    }
    l_max_t + ((sum_exp / n_f).ln() - alpha.ln()) / z_star
}

#[derive(Clone, Debug)]
pub struct DrawdownRiskMetrics<T: Real> {
    pub cdar: T,
    pub max_drawdown: T,
    pub average_drawdown: T,
    pub drawdowns: Vec<T>,
}

pub fn conditional_drawdown_at_risk<T: Real>(
    period_returns: &[T],
    confidence: f64,
) -> DrawdownRiskMetrics<T> {
    assert!(!period_returns.is_empty());
    let alpha = validate_confidence(confidence);
    let n = period_returns.len();
    let mut cum = Vec::with_capacity(n + 1);
    cum.push(T::zero());
    let mut acc = T::zero();
    for &r in period_returns {
        acc += r;
        cum.push(acc);
    }
    let mut drawdowns = Vec::with_capacity(n);
    let mut avg_dd = T::zero();
    let mut max_idx = 0usize;
    for t in 1..=n {
        let mut peak_idx = 0usize;
        let mut peak_val = f64::NEG_INFINITY;
        for (s, w_s) in cum.iter().enumerate().take(t + 1) {
            if w_s.val() > peak_val {
                peak_val = w_s.val();
                peak_idx = s;
            }
        }
        let dd = (cum[peak_idx] - cum[t]).max_val(T::zero());
        if dd.val() > drawdowns.get(max_idx).map(|x: &T| x.val()).unwrap_or(-1.0) {
            max_idx = t - 1;
        }
        avg_dd += dd;
        drawdowns.push(dd);
    }
    let max_dd = drawdowns[max_idx];
    avg_dd /= T::from_f64(n as f64);
    let mut indexed: Vec<(f64, usize)> = drawdowns
        .iter()
        .enumerate()
        .map(|(i, d)| (d.val(), i))
        .collect();
    indexed.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mass = alpha * n as f64;
    let whole = mass.floor() as usize;
    let rem = mass - whole as f64;
    let mut tail_sum = T::zero();
    for &(_, idx) in indexed.iter().take(whole) {
        tail_sum += drawdowns[idx];
    }
    if whole < indexed.len() {
        tail_sum += drawdowns[indexed[whole].1] * rem;
    }
    DrawdownRiskMetrics {
        cdar: tail_sum / mass,
        max_drawdown: max_dd,
        average_drawdown: avg_dd,
        drawdowns,
    }
}

#[derive(Clone, Debug)]
pub struct FactorRiskAttribution<T: Real> {
    pub total_volatility: T,
    pub systematic_volatility: T,
    pub idiosyncratic_volatility: T,
    pub portfolio_betas: Vec<T>,
    pub factor_risk_contributions: Vec<T>,
    pub idiosyncratic_risk_contribution: T,
    pub r_squared: T,
}

#[allow(clippy::needless_range_loop)]
pub fn factor_risk_attribution<T: Real>(
    weights: &[T],
    factor_loadings: &[Vec<T>],
    factor_covariance: &[Vec<T>],
    idiosyncratic_vars: &[T],
) -> FactorRiskAttribution<T> {
    let n = weights
        .len()
        .min(factor_loadings.len())
        .min(idiosyncratic_vars.len());
    let k = factor_covariance.len();
    let mut betas = vec![T::zero(); k];
    for i in 0..n {
        for f in 0..k.min(factor_loadings[i].len()) {
            betas[f] += weights[i] * factor_loadings[i][f];
        }
    }
    let mut sigma_f_beta = vec![T::zero(); k];
    let mut sys_var = T::zero();
    for a in 0..k {
        for b in 0..k.min(factor_covariance[a].len()) {
            sigma_f_beta[a] += factor_covariance[a][b] * betas[b];
        }
        sys_var += betas[a] * sigma_f_beta[a];
    }
    let mut idio_var = T::zero();
    for i in 0..n {
        idio_var += weights[i] * weights[i] * idiosyncratic_vars[i];
    }
    let tot_var = (sys_var + idio_var).max_val(T::from_f64(1e-16));
    let tot_vol = tot_var.sqrt();
    let frc: Vec<T> = (0..k)
        .map(|a| betas[a] * sigma_f_beta[a] / tot_vol)
        .collect();
    let irc = idio_var / tot_vol;
    FactorRiskAttribution {
        total_volatility: tot_vol,
        systematic_volatility: sys_var.max_val(T::zero()).sqrt(),
        idiosyncratic_volatility: idio_var.max_val(T::zero()).sqrt(),
        portfolio_betas: betas,
        factor_risk_contributions: frc,
        idiosyncratic_risk_contribution: irc,
        r_squared: sys_var / tot_var,
    }
}

#[derive(Clone, Debug)]
pub struct MinimumTorsionEnbReport<T: Real> {
    pub enb: T,
    pub risk_contributions: Vec<T>,
    pub probabilities: Vec<T>,
    pub portfolio_variance: T,
}

fn invert_square_matrix<T: Real>(mat: &[Vec<T>]) -> Result<Vec<Vec<T>>, String> {
    let n = mat.len();
    let mut a = mat.to_vec();
    let mut inv = vec![vec![T::zero(); n]; n];
    for (i, row) in inv.iter_mut().enumerate().take(n) {
        row[i] = T::one();
    }
    for col in 0..n {
        let mut pivot = col;
        let mut max_abs = a[col][col].val().abs();
        for (r, row) in a.iter().enumerate().take(n).skip(col + 1) {
            let v = row[col].val().abs();
            if v > max_abs {
                max_abs = v;
                pivot = r;
            }
        }
        if max_abs < 1e-14 {
            return Err("singular matrix in Denman-Beavers".into());
        }
        if pivot != col {
            a.swap(pivot, col);
            inv.swap(pivot, col);
        }
        let diag = a[col][col];
        for j in 0..n {
            a[col][j] /= diag;
            inv[col][j] /= diag;
        }
        for r in 0..n {
            if r == col {
                continue;
            }
            let factor = a[r][col];
            for j in 0..n {
                let a_cj = a[col][j];
                let inv_cj = inv[col][j];
                a[r][j] -= factor * a_cj;
                inv[r][j] -= factor * inv_cj;
            }
        }
    }
    Ok(inv)
}

pub fn try_effective_number_of_bets_minimum_torsion<T: Real>(
    weights: &[T],
    cov: &[Vec<T>],
) -> Result<MinimumTorsionEnbReport<T>, String> {
    let n = weights.len();
    if n == 0
        || cov.len() != n
        || cov.iter().any(|r| r.len() != n)
        || weights.iter().any(|x| !x.val().is_finite())
        || cov.iter().flatten().any(|x| !x.val().is_finite())
    {
        return Err("invalid dimensions".into());
    }
    let mut sig = Vec::with_capacity(n);
    for (i, row) in cov.iter().enumerate().take(n) {
        if row[i].val() <= 0.0 {
            return Err("non-positive diagonal variance".into());
        }
        sig.push(row[i].sqrt());
    }
    let mut y_mat = vec![vec![T::zero(); n]; n];
    let mut z_mat = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        z_mat[i][i] = T::one();
        for j in 0..n {
            y_mat[i][j] = if i == j {
                T::one()
            } else {
                (cov[i][j] / (sig[i] * sig[j]))
                    .max_val(T::from_f64(-0.999))
                    .min_val(T::from_f64(0.999))
            };
        }
    }
    // 8-step Denman-Beavers quadratic iteration for C^{1/2} (y_mat) and C^{-1/2} (z_mat)
    for _ in 0..8 {
        let y_inv = invert_square_matrix(&y_mat)?;
        let z_inv = invert_square_matrix(&z_mat)?;
        for i in 0..n {
            for j in 0..n {
                y_mat[i][j] = (y_mat[i][j] + z_inv[i][j]) * 0.5;
                z_mat[i][j] = (z_mat[i][j] + y_inv[i][j]) * 0.5;
            }
        }
    }
    // Uncorrelated Riccati / minimum-torsion bet exposures b_i = \sum_j C^{1/2}_{ij} \sigma_j w_j
    let mut v = Vec::with_capacity(n);
    let mut total = T::zero();
    for row in y_mat.iter().take(n) {
        let mut b_i = T::zero();
        for j in 0..n {
            b_i += row[j] * sig[j] * weights[j];
        }
        let v_i = b_i * b_i;
        v.push(v_i);
        total += v_i;
    }
    let total_safe = total.max_val(T::from_f64(1e-30));
    let p: Vec<T> = v.iter().map(|&x| x / total_safe).collect();
    let ent = p
        .iter()
        .map(|&x| {
            let q = x.max_val(T::from_f64(1e-16));
            -q * q.ln()
        })
        .fold(T::zero(), |a, b| a + b);
    Ok(MinimumTorsionEnbReport {
        enb: ent.exp(),
        risk_contributions: v,
        probabilities: p,
        portfolio_variance: total,
    })
}

pub fn effective_number_of_bets_minimum_torsion<T: Real>(
    weights: &[T],
    cov: &[Vec<T>],
) -> MinimumTorsionEnbReport<T> {
    try_effective_number_of_bets_minimum_torsion(weights, cov).expect("invalid covariance")
}

#[derive(Clone, Debug)]
pub struct MeanCdarPortfolioSpec<T: Real = f64> {
    pub path_asset_returns: Vec<Vec<T>>,
    pub expected_returns: Vec<T>,
    pub cdar_aversion: T,
    pub confidence: f64,
    pub smoothing_eps: f64,
}

#[derive(Clone, Debug)]
pub struct MeanCdarPortfolioResult<T: Real> {
    pub weights: Vec<T>,
    pub cdar: T,
    pub dar_threshold: T,
    pub max_drawdown: T,
    pub expected_return: T,
    pub drawdown_path: Vec<T>,
}

pub fn optimize_mean_cdar_portfolio<T: Real>(
    spec: &MeanCdarPortfolioSpec<T>,
) -> MeanCdarPortfolioResult<T> {
    let t_len = spec.path_asset_returns.len();
    let n = spec.expected_returns.len();
    if t_len == 0 || n == 0 || spec.path_asset_returns.iter().any(|r| r.len() < n) {
        return MeanCdarPortfolioResult {
            weights: Vec::new(),
            cdar: T::zero(),
            dar_threshold: T::zero(),
            max_drawdown: T::zero(),
            expected_return: T::zero(),
            drawdown_path: Vec::new(),
        };
    }

    let eps_f = spec.smoothing_eps.clamp(1e-5, 2e-4);
    let eps = T::from_f64(eps_f);
    let tail = (1.0 - spec.confidence.clamp(0.5, 0.99)).max(1.0 / t_len as f64);
    let mut logits = vec![T::zero(); n];

    let eval_drawdowns_and_grad = |w: &[T]| -> (Vec<T>, T, T, Vec<T>) {
        let mut cum_w = T::zero();
        let mut peak = T::zero();
        let mut d_cum_dw = vec![T::zero(); n];
        let mut d_peak_dw = vec![T::zero(); n];
        let mut dd = Vec::with_capacity(t_len);
        let mut d_dd_dw = vec![vec![T::zero(); n]; t_len];

        for (t, row) in spec.path_asset_returns.iter().enumerate() {
            for i in 0..n {
                cum_w += row[i] * w[i];
                d_cum_dw[i] += row[i];
            }
            let diff = peak - cum_w;
            let rad = (diff * diff + eps * eps).sqrt();
            let w_peak = (T::one() + diff / rad) * 0.5;
            let w_cum = T::one() - w_peak;
            peak = (peak + cum_w + rad - eps) * 0.5;
            for i in 0..n {
                d_peak_dw[i] = w_peak * d_peak_dw[i] + w_cum * d_cum_dw[i];
                d_dd_dw[t][i] = d_peak_dw[i] - d_cum_dw[i];
            }
            dd.push((peak - cum_w).max_val(T::zero()));
        }

        // Solve Rockafellar-Uryasev first-order quantile optimality condition:
        // sum_t 0.5 * (1 + (dd_t - zeta) / sqrt((dd_t - zeta)^2 + eps^2)) = tail * t_len
        let max_dd_f = dd.iter().map(|x| x.val()).fold(0.0_f64, f64::max);
        let target_cnt = tail * t_len as f64;
        let mut lo = -eps_f;
        let mut hi = max_dd_f + eps_f;
        for _ in 0..32 {
            let mid = 0.5 * (lo + hi);
            let cnt: f64 = dd
                .iter()
                .map(|d_t| {
                    let u = d_t.val() - mid;
                    0.5 * (1.0 + u / (u * u + eps_f * eps_f).sqrt())
                })
                .sum();
            if cnt > target_cnt {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let mut zeta = T::from_f64(0.5 * (lo + hi));
        for _ in 0..6 {
            let mut sum_prob = T::zero();
            let mut sum_dens = T::zero();
            for &d_t in &dd {
                let u = d_t - zeta;
                let r = (u * u + eps * eps).sqrt();
                sum_prob += (T::one() + u / r) * 0.5;
                sum_dens += (eps * eps) / (r * r * r * 2.0);
            }
            let g_z = T::one() - sum_prob / target_cnt;
            let h_z = (sum_dens / target_cnt).max_val(T::from_f64(1e-6));
            zeta = (zeta - g_z / h_z)
                .max_val(T::zero())
                .min_val(T::from_f64(max_dd_f));
        }

        let mut cdar_val = zeta;
        let mut grad_w = vec![T::zero(); n];
        for i in 0..n {
            grad_w[i] = -spec.expected_returns[i];
        }
        let scale = T::one() / target_cnt;
        for t in 0..t_len {
            let u = dd[t] - zeta;
            let r = (u * u + eps * eps).sqrt();
            // Zero-anchored smooth hinge \tilde{h}_eps(u) = (u + sqrt(u^2 + eps^2) - eps) / 2
            let hinge = (u + r - eps) * 0.5;
            let d_hinge = (T::one() + u / r) * 0.5;
            cdar_val += hinge * scale;
            for i in 0..n {
                grad_w[i] += spec.cdar_aversion * d_hinge * scale * d_dd_dw[t][i];
            }
        }
        (dd, zeta, cdar_val, grad_w)
    };

    let softmax = |lg: &[T]| -> Vec<T> {
        let mx = lg.iter().map(|x| x.val()).fold(f64::NEG_INFINITY, f64::max);
        let mut ex = vec![T::zero(); n];
        let mut s = T::zero();
        for i in 0..n {
            ex[i] = (lg[i] - mx).exp();
            s += ex[i];
        }
        ex.into_iter().map(|v| v / s).collect()
    };

    for step in 0..18 {
        let w = softmax(&logits);
        let (_, _, _, grad_w) = eval_drawdowns_and_grad(&w);
        let mut dot = T::zero();
        for i in 0..n {
            dot += w[i] * grad_w[i];
        }
        let lr = 2.5 / (1.0 + 0.15 * step as f64);
        for i in 0..n {
            logits[i] -= w[i] * (grad_w[i] - dot) * lr;
        }
    }

    let weights = softmax(&logits);
    let (drawdown_path, dar_threshold, cdar, _) = eval_drawdowns_and_grad(&weights);
    let max_dd = drawdown_path
        .iter()
        .copied()
        .fold(T::zero(), |a, b| a.max_val(b));
    let mut exp_ret = T::zero();
    for i in 0..n {
        exp_ret += weights[i] * spec.expected_returns[i];
    }

    MeanCdarPortfolioResult {
        weights,
        cdar,
        dar_threshold,
        max_drawdown: max_dd,
        expected_return: exp_ret,
        drawdown_path,
    }
}

#[cfg(test)]
mod phase16_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn evar_cdar_and_factor_euler_invariants() {
        let rets = vec![
            -0.05, -0.035, -0.025, -0.02, -0.015, -0.01, -0.005, 0.0, 0.005, 0.01, 0.015, 0.02,
            0.025, 0.03, 0.035, 0.04,
        ];
        let var90 = empirical_var(&rets, 0.90);
        let es90 = empirical_expected_shortfall(&rets, 0.90);
        let evar90 = entropic_value_at_risk(&rets, 0.90);
        assert!(var90 <= es90 + 1e-12 && es90 <= evar90 + 1e-12 && evar90 <= 0.05 + 1e-12);
        let dd = conditional_drawdown_at_risk(&rets, 0.75);
        assert!(dd.average_drawdown <= dd.cdar + 1e-12 && dd.cdar <= dd.max_drawdown + 1e-12);
        let w = vec![0.4, 0.6];
        let b = vec![vec![1.1, 0.3], vec![0.8, 0.9]];
        let fc = vec![vec![0.04, 0.005], vec![0.005, 0.02]];
        let idio = vec![0.01, 0.015];
        let attr = factor_risk_attribution(&w, &b, &fc, &idio);
        let sum_rc = attr.factor_risk_contributions.iter().sum::<f64>()
            + attr.idiosyncratic_risk_contribution;
        assert!((sum_rc - attr.total_volatility).abs() < 1e-12);
        let h = 1e-5;
        let eval_r0 = |r0: f64| {
            let mut r = rets.clone();
            r[0] = r0;
            entropic_value_at_risk(&r, 0.90)
        };
        let fd = (eval_r0(-0.05 + h) - eval_r0(-0.05 - h)) / (2.0 * h);
        type D = Dual<1>;
        let mut rd: Vec<D> = rets.iter().map(|&x| D::from_value(x)).collect();
        rd[0] = D::variable(-0.05, 0);
        let ad = entropic_value_at_risk(&rd, 0.90);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn enb_riccati_root_identity_concentration_and_ad() {
        let c_id = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        assert!(
            (effective_number_of_bets_minimum_torsion(&[0.5, 0.5], &c_id).enb - 2.0).abs() < 1e-12
        );
        assert!(
            (effective_number_of_bets_minimum_torsion(&[1.0, 0.0], &c_id).enb - 1.0).abs() < 1e-12
        );

        let c_corr = vec![vec![0.04, 0.015], vec![0.015, 0.09]];
        let rep = effective_number_of_bets_minimum_torsion(&[0.6, 0.4], &c_corr);
        assert!(rep.enb > 1.0 && rep.enb <= 2.0);
        assert!((rep.probabilities[0] + rep.probabilities[1] - 1.0).abs() < 1e-12);

        let h = 1e-5;
        let eval_w = |w0: f64| {
            effective_number_of_bets_minimum_torsion(&[w0, 1.0 - w0], &c_corr).enb
        };
        let fd = (eval_w(0.6 + h) - eval_w(0.6 - h)) / (2.0 * h);
        type D = Dual<1>;
        let w_d = [D::variable(0.6, 0), D::from_value(1.0) - D::variable(0.6, 0)];
        let c_d = vec![
            vec![D::from_value(0.04), D::from_value(0.015)],
            vec![D::from_value(0.015), D::from_value(0.09)],
        ];
        let ad = effective_number_of_bets_minimum_torsion(&w_d, &c_d).enb;
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase19_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn mean_cdar_optimization_simplex_aversion_monotonicity_and_ad() {
        let path = vec![
            vec![0.04, 0.01],
            vec![-0.08, -0.01],
            vec![-0.06, -0.005],
            vec![0.05, 0.012],
            vec![-0.05, -0.008],
            vec![0.03, 0.010],
        ];
        let low_av = MeanCdarPortfolioSpec {
            path_asset_returns: path.clone(),
            expected_returns: vec![0.10, 0.04],
            cdar_aversion: 0.2,
            confidence: 0.80,
            smoothing_eps: 1e-3,
        };
        let high_av = MeanCdarPortfolioSpec {
            cdar_aversion: 8.0,
            ..low_av.clone()
        };
        let r_low = optimize_mean_cdar_portfolio(&low_av);
        let r_high = optimize_mean_cdar_portfolio(&high_av);
        assert!((r_high.weights.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!(r_high.weights.iter().all(|&w| w >= 0.0));
        // Higher CDaR aversion shifts weight toward the low-drawdown asset 1 and lowers portfolio CDaR
        assert!(r_high.cdar <= r_low.cdar);
        assert!(r_high.weights[1] > r_low.weights[1]);
        assert!(r_high.dar_threshold <= r_high.cdar + 1e-6);
        assert!(r_high.cdar <= r_high.max_drawdown + 1e-6);

        let h = 1e-5;
        let eval_mu0 = |mu0: f64| {
            let mut s = high_av.clone();
            s.expected_returns[0] = mu0;
            optimize_mean_cdar_portfolio(&s).weights[0]
        };
        let fd = (eval_mu0(0.10 + h) - eval_mu0(0.10 - h)) / (2.0 * h);
        type D = Dual<1>;
        let path_d: Vec<Vec<D>> = path
            .iter()
            .map(|r| r.iter().map(|&x| D::from_value(x)).collect())
            .collect();
        let spec_d = MeanCdarPortfolioSpec {
            path_asset_returns: path_d,
            expected_returns: vec![D::variable(0.10, 0), D::from_value(0.04)],
            cdar_aversion: D::from_value(8.0),
            confidence: 0.80,
            smoothing_eps: 1e-3,
        };
        let ad = optimize_mean_cdar_portfolio(&spec_d).weights[0];
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
