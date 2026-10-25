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
        let eval_w =
            |w0: f64| effective_number_of_bets_minimum_torsion(&[w0, 1.0 - w0], &c_corr).enb;
        let fd = (eval_w(0.6 + h) - eval_w(0.6 - h)) / (2.0 * h);
        type D = Dual<1>;
        let w_d = [
            D::variable(0.6, 0),
            D::from_value(1.0) - D::variable(0.6, 0),
        ];
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

#[derive(Clone, Debug)]
pub struct ExpectileRiskAndPortfolioSpec<T: Real = f64> {
    pub scenario_returns: Vec<Vec<T>>,
    pub initial_weights: Vec<T>,
    pub tau_confidence: f64,
    pub return_weight_lambda: T,
    pub smoothing_eps: f64,
}

#[derive(Clone, Debug)]
pub struct ExpectileRiskAndPortfolioReport<T: Real = f64> {
    pub initial_expectile_exvar: T,
    pub initial_component_exvar: Vec<T>,
    pub optimal_weights: Vec<T>,
    pub optimal_expectile_exvar: T,
    pub optimal_component_exvar: Vec<T>,
    pub optimal_expected_return: T,
    pub omega_gain_loss_ratio: T,
}

fn solve_expectile_and_grad<T: Real>(
    scenario_returns: &[Vec<T>],
    weights: &[T],
    tau: f64,
    eps_f64: f64,
) -> (T, Vec<T>, T) {
    let s_cnt = scenario_returns.len();
    let n = weights.len();
    if s_cnt == 0 || n == 0 {
        return (T::zero(), vec![T::zero(); n], T::one());
    }
    let tau_c = tau.clamp(0.5, 0.995);
    let eps = T::from_f64(eps_f64.max(1e-5));
    let losses: Vec<T> = scenario_returns
        .iter()
        .map(|row| {
            -row.iter()
                .zip(weights)
                .fold(T::zero(), |acc, (&r, &w)| acc + r * w)
        })
        .collect();

    // 32-step f64 bisection + 6-step T: Real Newton root solve on the asymmetric expectile balance
    let mut lo = losses.iter().map(|l| l.val()).fold(f64::INFINITY, f64::min) - 1.0;
    let mut hi = losses
        .iter()
        .map(|l| l.val())
        .fold(f64::NEG_INFINITY, f64::max)
        + 1.0;
    let eval_f64 = |e: f64| -> f64 {
        let mut s = 0.0_f64;
        for l in &losses {
            let u = l.val() - e;
            let ind = 0.5 * (1.0 + u / (u * u + eps_f64 * eps_f64).sqrt());
            let w_asym = (1.0 - tau_c) + (2.0 * tau_c - 1.0) * ind;
            s += w_asym * u;
        }
        s
    };
    for _ in 0..32 {
        let mid = 0.5 * (lo + hi);
        if eval_f64(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mut e_star = T::from_f64(0.5 * (lo + hi));
    for _ in 0..6 {
        let mut f = T::zero();
        let mut df_de = T::zero();
        for &l in &losses {
            let u = l - e_star;
            let rad = (u * u + eps * eps).sqrt();
            let ind = (T::one() + u / rad) * 0.5;
            let dind_du = eps * eps / (rad * rad * rad) * 0.5;
            let w_asym = T::from_f64(1.0 - tau_c) + ind * (2.0 * tau_c - 1.0);
            let dw_du = dind_du * (2.0 * tau_c - 1.0);
            f += w_asym * u;
            df_de -= w_asym + dw_du * u;
        }
        e_star -= f / df_de.min_val(T::from_f64(-1e-12));
    }

    // Exact Euler representation at the converged expectile root:
    // sum_s w_asym_s * (L_s - e_star) = 0  ==>  e_star = (sum_s w_asym_s L_s) / (sum_s w_asym_s)
    // Evaluating w_asym_s at the converged primal root gives exact 1-homogeneity
    // sum_j w_j * grad_j == e_star to machine precision (< 1e-14).
    let mut denom = T::zero();
    let mut asym_weights = Vec::with_capacity(s_cnt);
    let mut weighted_loss = T::zero();
    let mut pos_sum = T::zero();
    let mut neg_sum = T::zero();
    for &l in &losses {
        let u = l - e_star;
        let rad = (u * u + eps * eps).sqrt();
        let ind = (T::one() + u / rad) * 0.5;
        let w_asym = T::from_f64(1.0 - tau_c) + ind * (2.0 * tau_c - 1.0);
        asym_weights.push(w_asym);
        denom += w_asym;
        weighted_loss += w_asym * l;
        pos_sum += (u + rad) * 0.5;
        neg_sum += (-u + rad) * 0.5;
    }
    let denom_safe = denom.max_val(T::from_f64(1e-12));
    let e_fixed = weighted_loss / denom_safe;
    let mut grad = vec![T::zero(); n];
    for j in 0..n {
        let mut num = T::zero();
        for s in 0..s_cnt {
            num -= asym_weights[s] * scenario_returns[s][j];
        }
        grad[j] = num / denom_safe;
    }
    let omega = neg_sum / pos_sum.max_val(T::from_f64(1e-12));
    (e_fixed, grad, omega)
}

/// Bellini & Bignozzi (2015) / Ziegel (2016) Coherent & Elicitable Expectile Value-at-Risk
/// ($\mathrm{EXVaR}_\tau$), Euler Risk Decomposition, and Mean-Expectile Portfolio Optimization.
pub fn evaluate_and_optimize_expectile_portfolio<T: Real>(
    spec: &ExpectileRiskAndPortfolioSpec<T>,
) -> ExpectileRiskAndPortfolioReport<T> {
    let s_cnt = spec.scenario_returns.len();
    let n = spec.initial_weights.len();
    if s_cnt == 0 || n == 0 {
        return ExpectileRiskAndPortfolioReport {
            initial_expectile_exvar: T::zero(),
            initial_component_exvar: Vec::new(),
            optimal_weights: Vec::new(),
            optimal_expectile_exvar: T::zero(),
            optimal_component_exvar: Vec::new(),
            optimal_expected_return: T::zero(),
            omega_gain_loss_ratio: T::one(),
        };
    }
    let eps_f64 = spec.smoothing_eps.max(1e-4);
    let (init_exvar, init_grad, _) = solve_expectile_and_grad(
        &spec.scenario_returns,
        &spec.initial_weights,
        spec.tau_confidence,
        eps_f64,
    );
    let init_comp: Vec<T> = (0..n)
        .map(|j| spec.initial_weights[j] * init_grad[j])
        .collect();

    let mean_rets: Vec<T> = (0..n)
        .map(|j| {
            spec.scenario_returns
                .iter()
                .fold(T::zero(), |acc, r| acc + r[j])
                / (s_cnt as f64)
        })
        .collect();

    // Smooth entropic mirror descent (softmax over logits) in T: Real so all Dual<N> tangents
    // propagate smoothly through the convex simplex optimization without hard clipping.
    let mut logits: Vec<T> = spec
        .initial_weights
        .iter()
        .map(|&wj| wj.max_val(T::from_f64(1e-6)).ln())
        .collect();
    let mut w = vec![T::from_f64(1.0 / n as f64); n];
    for iter in 0..24 {
        let max_l = logits.iter().fold(logits[0], |a, &b| a.max_val(b));
        let exps: Vec<T> = logits.iter().map(|&l| (l - max_l).exp()).collect();
        let sum_e = exps
            .iter()
            .fold(T::zero(), |a, &b| a + b)
            .max_val(T::from_f64(1e-14));
        for j in 0..n {
            w[j] = exps[j] / sum_e;
        }
        let (_, g_exv, _) =
            solve_expectile_and_grad(&spec.scenario_returns, &w, spec.tau_confidence, eps_f64);
        let step = 1.2 / (1.0 + 0.15 * iter as f64);
        for j in 0..n {
            let g_j = g_exv[j] - spec.return_weight_lambda * mean_rets[j];
            logits[j] -= g_j * step;
        }
    }
    let max_l = logits.iter().fold(logits[0], |a, &b| a.max_val(b));
    let exps: Vec<T> = logits.iter().map(|&l| (l - max_l).exp()).collect();
    let sum_e = exps
        .iter()
        .fold(T::zero(), |a, &b| a + b)
        .max_val(T::from_f64(1e-14));
    for j in 0..n {
        w[j] = exps[j] / sum_e;
    }

    let (opt_exvar, opt_grad, omega) =
        solve_expectile_and_grad(&spec.scenario_returns, &w, spec.tau_confidence, eps_f64);
    let opt_comp: Vec<T> = (0..n).map(|j| w[j] * opt_grad[j]).collect();
    let opt_ret = w
        .iter()
        .zip(&mean_rets)
        .fold(T::zero(), |acc, (&wj, &mj)| acc + wj * mj);

    ExpectileRiskAndPortfolioReport {
        initial_expectile_exvar: init_exvar,
        initial_component_exvar: init_comp,
        optimal_weights: w,
        optimal_expectile_exvar: opt_exvar,
        optimal_component_exvar: opt_comp,
        optimal_expected_return: opt_ret,
        omega_gain_loss_ratio: omega,
    }
}

#[cfg(test)]
mod phase20_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn expectile_exvar_euler_additivity_tau_monotonicity_and_ad() {
        let scen = vec![
            vec![-0.07, -0.01],
            vec![-0.05, -0.008],
            vec![-0.02, 0.005],
            vec![0.01, 0.009],
            vec![0.04, 0.012],
            vec![0.06, 0.015],
        ];
        let spec50 = ExpectileRiskAndPortfolioSpec {
            scenario_returns: scen.clone(),
            initial_weights: vec![0.5, 0.5],
            tau_confidence: 0.50,
            return_weight_lambda: 0.2,
            smoothing_eps: 1e-4,
        };
        let spec90 = ExpectileRiskAndPortfolioSpec {
            tau_confidence: 0.90,
            ..spec50.clone()
        };
        let r50 = evaluate_and_optimize_expectile_portfolio(&spec50);
        let r90 = evaluate_and_optimize_expectile_portfolio(&spec90);
        // Higher confidence tau_confidence = 0.90 increases Expectile EXVaR
        assert!(r90.initial_expectile_exvar > r50.initial_expectile_exvar);
        // Exact Euler homogeneity: sum of component EXVaR == total EXVaR
        let sum_init: f64 = r90.initial_component_exvar.iter().sum();
        let sum_opt: f64 = r90.optimal_component_exvar.iter().sum();
        assert!((sum_init - r90.initial_expectile_exvar).abs() < 1e-12);
        assert!((sum_opt - r90.optimal_expectile_exvar).abs() < 1e-12);

        let h = 1e-5;
        let eval_w0 = |w0: f64| {
            let mut s = spec90.clone();
            s.initial_weights = vec![w0, 1.0 - w0];
            evaluate_and_optimize_expectile_portfolio(&s).initial_expectile_exvar
        };
        let fd = (eval_w0(0.5 + h) - eval_w0(0.5 - h)) / (2.0 * h);
        type D = Dual<1>;
        let scen_d: Vec<Vec<D>> = scen
            .iter()
            .map(|r| r.iter().map(|&x| D::from_value(x)).collect())
            .collect();
        let w0_d = D::variable(0.5, 0);
        let spec_d = ExpectileRiskAndPortfolioSpec {
            scenario_returns: scen_d,
            initial_weights: vec![w0_d, D::one() - w0_d],
            tau_confidence: 0.90,
            return_weight_lambda: D::from_value(0.2),
            smoothing_eps: 1e-4,
        };
        let ad = evaluate_and_optimize_expectile_portfolio(&spec_d).initial_expectile_exvar;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[derive(Clone, Debug)]
pub struct McNeilEmbrechtsPotGpdBacktestSpec<T: Real = f64> {
    pub historical_losses: Vec<T>,
    pub threshold_u: T,
    pub confidence_alpha: T,
    pub spectral_risk_aversion_k: T,
    pub backtest_var_threshold: T,
}
#[derive(Clone, Debug)]
pub struct McNeilEmbrechtsPotGpdBacktestReport<T: Real = f64> {
    pub exceedance_count: usize,
    pub fitted_gpd_xi: T,
    pub fitted_gpd_beta: T,
    pub pot_gpd_var: T,
    pub pot_gpd_expected_shortfall: T,
    pub spectral_risk_measure: T,
    pub kupiec_lr_uc: T,
    pub christoffersen_lr_ind: T,
    pub christoffersen_lr_cc: T,
    pub transition_counts_n00_n01_n10_n11: [usize; 4],
}
pub fn evaluate_mcneil_embrechts_pot_gpd_and_christoffersen_backtest<T: Real>(
    s: &McNeilEmbrechtsPotGpdBacktestSpec<T>,
) -> McNeilEmbrechtsPotGpdBacktestReport<T> {
    let mut ex = s
        .historical_losses
        .iter()
        .copied()
        .filter(|x| x.val() > s.threshold_u.val())
        .map(|x| x - s.threshold_u)
        .collect::<Vec<_>>();
    ex.sort_by(|a, b| a.val().partial_cmp(&b.val()).unwrap());
    let n = ex.len();
    let mut m0 = T::zero();
    let mut m1 = T::zero();
    for (j, y) in ex.iter().enumerate() {
        m0 += *y / T::from_f64(n.max(1) as f64);
        if n > 1 {
            m1 += T::from_f64((n - j - 1) as f64 / (n - 1) as f64) * *y / T::from_f64(n as f64);
        }
    }
    let den = m0 - T::from_f64(2.0) * m1;
    let xi = (T::from_f64(2.0) - m0 / den.max_val(T::from_f64(1e-12)))
        .max_val(T::from_f64(-0.45))
        .min_val(T::from_f64(0.85));
    let beta =
        (T::from_f64(2.0) * m0 * m1 / den.max_val(T::from_f64(1e-12))).max_val(T::from_f64(1e-12));
    let nn = s.historical_losses.len().max(1) as f64;
    let tail = (T::from_f64(nn / n.max(1) as f64) * (T::one() - s.confidence_alpha))
        .max_val(T::from_f64(1e-12));
    let var = s.threshold_u
        + if xi.val().abs() < 1e-10 {
            beta * (-tail.ln())
        } else {
            beta / xi * ((-xi * tail.ln()).exp() - T::one())
        };
    let es = (var + beta - xi * s.threshold_u) / (T::one() - xi).max_val(T::from_f64(1e-12));
    let mut sorted = s.historical_losses.clone();
    sorted.sort_by(|a, b| a.val().partial_cmp(&b.val()).unwrap());
    let k = s.spectral_risk_aversion_k;
    let mut spec = T::zero();
    for (i, x) in sorted.iter().enumerate() {
        let i0 = i as f64 / nn;
        let i1 = (i + 1) as f64 / nn;
        let w = if k.val().abs() < 1e-12 {
            T::from_f64(1.0 / nn)
        } else {
            ((-k * (T::one() - T::from_f64(i1))).exp() - (-k * (T::one() - T::from_f64(i0))).exp())
                / (T::one() - (-k).exp())
        };
        spec += w * *x;
    }
    let mut c = [0usize; 4];
    for w in s.historical_losses.windows(2) {
        let prev = w[0].val() > s.backtest_var_threshold.val();
        let cur = w[1].val() > s.backtest_var_threshold.val();
        c[(prev as usize) * 2 + (cur as usize)] += 1;
    }
    let nt = c.iter().sum::<usize>().max(1);
    let n_ex = c[1] + c[3];
    let p0 = (T::one() - s.confidence_alpha)
        .max_val(T::from_f64(1e-12))
        .min_val(T::one() - T::from_f64(1e-12));
    let ph = T::from_f64(n_ex as f64 / nt as f64)
        .max_val(T::from_f64(1e-12))
        .min_val(T::one() - T::from_f64(1e-12));
    let ll = |n0: usize, n1: usize, p: T| {
        T::from_f64(n0 as f64) * (T::one() - p).ln() + T::from_f64(n1 as f64) * p.ln()
    };
    let kupiec = T::from_f64(2.0) * (ll(nt - n_ex, n_ex, ph) - ll(nt - n_ex, n_ex, p0));
    let pi0 = T::from_f64(c[1] as f64 / (c[0] + c[1]).max(1) as f64)
        .max_val(T::from_f64(1e-12))
        .min_val(T::one() - T::from_f64(1e-12));
    let pi1 = T::from_f64(c[3] as f64 / (c[2] + c[3]).max(1) as f64)
        .max_val(T::from_f64(1e-12))
        .min_val(T::one() - T::from_f64(1e-12));
    let pall = T::from_f64((c[1] + c[3]) as f64 / nt as f64)
        .max_val(T::from_f64(1e-12))
        .min_val(T::one() - T::from_f64(1e-12));
    let ind = T::from_f64(2.0)
        * (ll(c[0], c[1], pi0) + ll(c[2], c[3], pi1) - ll(c[0] + c[2], c[1] + c[3], pall));
    McNeilEmbrechtsPotGpdBacktestReport {
        exceedance_count: n,
        fitted_gpd_xi: xi,
        fitted_gpd_beta: beta,
        pot_gpd_var: var,
        pot_gpd_expected_shortfall: es,
        spectral_risk_measure: spec,
        kupiec_lr_uc: kupiec.max_val(T::zero()),
        christoffersen_lr_ind: ind.max_val(T::zero()),
        christoffersen_lr_cc: (kupiec + ind).max_val(T::zero()),
        transition_counts_n00_n01_n10_n11: c,
    }
}

#[derive(Clone, Debug)]
pub struct BaconCarinoMultiPeriodBrinsonSpec<T: Real = f64> {
    pub portfolio_weights_by_period: Vec<Vec<T>>,
    pub benchmark_weights_by_period: Vec<Vec<T>>,
    pub portfolio_returns_by_period: Vec<Vec<T>>,
    pub benchmark_returns_by_period: Vec<Vec<T>>,
    pub fx_spot_returns_by_sector: Vec<T>,
    pub foreign_cash_rates_by_sector: Vec<T>,
    pub domestic_cash_rate: T,
}
#[derive(Clone, Debug)]
pub struct BaconCarinoMultiPeriodBrinsonReport<T: Real = f64> {
    pub cumulative_portfolio_return: T,
    pub cumulative_benchmark_return: T,
    pub total_geometric_active_return: T,
    pub carino_linked_allocation_total: T,
    pub carino_linked_selection_total: T,
    pub carino_linked_interaction_total: T,
    pub multi_period_linking_residual: T,
    pub karnosky_singer_currency_allocation: T,
    pub karnosky_singer_local_excess_allocation: T,
}
fn psi_bc<T: Real>(x: T, y: T) -> T {
    let d = x - y;
    let m = (T::from_f64(2.0) + x + y) / T::from_f64(2.0);
    if d.val().abs() > 1e-7 {
        ((T::one() + x).ln() - (T::one() + y).ln()) / d
    } else {
        (T::one() + d * d / (T::from_f64(12.0) * m * m)) / m
    }
}
pub fn evaluate_bacon_carino_brinson_and_karnosky_singer_attribution<T: Real>(
    s: &BaconCarinoMultiPeriodBrinsonSpec<T>,
) -> BaconCarinoMultiPeriodBrinsonReport<T> {
    let m = s.portfolio_weights_by_period.len();
    let mut rp = T::one();
    let mut rb = T::one();
    let mut aa = T::zero();
    let mut ss = T::zero();
    let mut ii = T::zero();
    for t in 0..m {
        let wp = &s.portfolio_weights_by_period[t];
        let wb = &s.benchmark_weights_by_period[t];
        let rp_t = wp
            .iter()
            .zip(&s.portfolio_returns_by_period[t])
            .map(|(w, r)| *w * *r)
            .fold(T::zero(), |a, x| a + x);
        let rb_t = wb
            .iter()
            .zip(&s.benchmark_returns_by_period[t])
            .map(|(w, r)| *w * *r)
            .fold(T::zero(), |a, x| a + x);
        rp *= T::one() + rp_t;
        rb *= T::one() + rb_t;
        let k = psi_bc(rp_t, rb_t);
        for i in 0..wp.len().min(wb.len()) {
            let rb_i = s.benchmark_returns_by_period[t][i];
            let active_i = s.portfolio_returns_by_period[t][i] - rb_i;
            aa += (wp[i] - wb[i]) * (rb_i - rb_t) * k;
            ss += wb[i] * active_i * k;
            ii += (wp[i] - wb[i]) * active_i * k;
        }
    }
    let a = rp - T::one();
    let b = rb - T::one();
    let active = a - b;
    let k_total = psi_bc(a, b).max_val(T::from_f64(1e-12));
    let alloc = aa / k_total;
    let sel = ss / k_total;
    let inter = ii / k_total;
    let residual = active - (alloc + sel + inter);
    let mut rc = Vec::new();
    let mut rl = Vec::new();
    for i in 0..s.fx_spot_returns_by_sector.len() {
        let foreign = s
            .foreign_cash_rates_by_sector
            .get(i)
            .copied()
            .unwrap_or(T::zero());
        rc.push(s.fx_spot_returns_by_sector[i] + foreign - s.domestic_cash_rate);
        let br = s
            .benchmark_returns_by_period
            .first()
            .and_then(|v| v.get(i))
            .copied()
            .unwrap_or(T::zero());
        rl.push(br - foreign);
    }
    let wb0 = s
        .benchmark_weights_by_period
        .first()
        .cloned()
        .unwrap_or_default();
    let wp0 = s
        .portfolio_weights_by_period
        .first()
        .cloned()
        .unwrap_or_default();
    let mut rcb = T::zero();
    let mut rlb = T::zero();
    for i in 0..rc.len().min(wb0.len()) {
        rcb += wb0[i] * rc[i];
        rlb += wb0[i] * rl[i];
    }
    let mut kcur = T::zero();
    let mut kloc = T::zero();
    for i in 0..rc.len().min(wp0.len()).min(wb0.len()) {
        let dw = wp0[i] - wb0[i];
        kcur += dw * (rc[i] - rcb);
        kloc += dw * (rl[i] - rlb);
    }
    BaconCarinoMultiPeriodBrinsonReport {
        cumulative_portfolio_return: a,
        cumulative_benchmark_return: b,
        total_geometric_active_return: active,
        carino_linked_allocation_total: alloc,
        carino_linked_selection_total: sel,
        carino_linked_interaction_total: inter,
        multi_period_linking_residual: residual,
        karnosky_singer_currency_allocation: kcur,
        karnosky_singer_local_excess_allocation: kloc,
    }
}

#[cfg(test)]
mod phase24_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(a: T) -> McNeilEmbrechtsPotGpdBacktestSpec<T> {
        McNeilEmbrechtsPotGpdBacktestSpec {
            historical_losses: vec![
                T::from_f64(1.0),
                T::from_f64(2.0),
                T::from_f64(3.0) + a,
                T::from_f64(4.0),
                T::from_f64(5.0),
            ],
            threshold_u: T::from_f64(1.0),
            confidence_alpha: T::from_f64(0.95),
            spectral_risk_aversion_k: T::from_f64(2.0),
            backtest_var_threshold: T::from_f64(3.0),
        }
    }
    #[test]
    fn pot_gpd_phase24() {
        let r = evaluate_mcneil_embrechts_pot_gpd_and_christoffersen_backtest(&mk(0.0));
        assert!(r.pot_gpd_expected_shortfall > r.pot_gpd_var && r.pot_gpd_var > 1.0);
        assert!(
            (r.christoffersen_lr_cc - (r.kupiec_lr_uc + r.christoffersen_lr_ind)).abs() < 1e-12
        );
        let h = 1e-5;
        let fd = (evaluate_mcneil_embrechts_pot_gpd_and_christoffersen_backtest(
            &McNeilEmbrechtsPotGpdBacktestSpec {
                confidence_alpha: 0.95 + h,
                ..mk(0.0)
            },
        )
        .pot_gpd_expected_shortfall
            - evaluate_mcneil_embrechts_pot_gpd_and_christoffersen_backtest(
                &McNeilEmbrechtsPotGpdBacktestSpec {
                    confidence_alpha: 0.95 - h,
                    ..mk(0.0)
                },
            )
            .pot_gpd_expected_shortfall)
            / (2.0 * h);
        let ad = evaluate_mcneil_embrechts_pot_gpd_and_christoffersen_backtest(
            &McNeilEmbrechtsPotGpdBacktestSpec {
                confidence_alpha: Dual::<1>::variable(0.95, 0),
                ..mk(Dual::zero())
            },
        )
        .pot_gpd_expected_shortfall;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
    fn br_spec<T: Real>(x: T) -> BaconCarinoMultiPeriodBrinsonSpec<T> {
        BaconCarinoMultiPeriodBrinsonSpec {
            portfolio_weights_by_period: vec![vec![T::from_f64(0.6), T::from_f64(0.4)]; 2],
            benchmark_weights_by_period: vec![vec![T::from_f64(0.5), T::from_f64(0.5)]; 2],
            portfolio_returns_by_period: vec![
                vec![T::from_f64(0.10) + x, T::from_f64(0.02)],
                vec![T::from_f64(0.04), T::from_f64(0.01)],
            ],
            benchmark_returns_by_period: vec![
                vec![T::from_f64(0.08), T::from_f64(0.01)],
                vec![T::from_f64(0.03), T::from_f64(0.015)],
            ],
            fx_spot_returns_by_sector: vec![T::from_f64(0.01), T::from_f64(-0.005)],
            foreign_cash_rates_by_sector: vec![T::from_f64(0.02), T::from_f64(0.01)],
            domestic_cash_rate: T::from_f64(0.03),
        }
    }
    #[test]
    fn brinson_carino_phase24() {
        let r = evaluate_bacon_carino_brinson_and_karnosky_singer_attribution(&br_spec(0.0));
        assert!(r.multi_period_linking_residual.abs() < 1e-12);
        let h = 1e-5;
        let fd = (evaluate_bacon_carino_brinson_and_karnosky_singer_attribution(&br_spec(h))
            .carino_linked_selection_total
            - evaluate_bacon_carino_brinson_and_karnosky_singer_attribution(&br_spec(-h))
                .carino_linked_selection_total)
            / (2.0 * h);
        let ad = evaluate_bacon_carino_brinson_and_karnosky_singer_attribution(&br_spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .carino_linked_selection_total;
        assert!((ad.d[0] - fd).abs() < 1e-9);
    }
}
