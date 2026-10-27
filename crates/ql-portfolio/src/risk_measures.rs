//! Coherent empirical and parametric risk measures plus Euler allocations.
#![allow(clippy::needless_range_loop)]
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

#[derive(Clone, Debug)]
pub struct HarveyLiuZhuMultipleTestingSharpeSpec<T: Real = f64> {
    pub focal_sharpe_ratio: T,
    pub benchmark_sharpe_sr0: T,
    pub sample_length_t: usize,
    pub return_skewness: T,
    pub return_kurtosis_pearson: T,
    pub num_independent_trials_m: usize,
    pub cross_sectional_sharpe_std: T,
    pub confidence_z_alpha: T,
    pub candidate_t_stats: Vec<T>,
}
#[derive(Clone, Debug)]
pub struct HarveyLiuZhuMultipleTestingSharpeReport<T: Real = f64> {
    pub sharpe_standard_error: T,
    pub probabilistic_sharpe_ratio_psr: T,
    pub minimum_track_record_length_mintrl: T,
    pub expected_max_false_sharpe_sr0_star: T,
    pub deflated_sharpe_ratio_dsr: T,
    pub raw_two_sided_p_values: Vec<T>,
    pub bonferroni_adjusted_p_values: Vec<T>,
    pub holm_adjusted_p_values: Vec<T>,
    pub bhy_fdr_adjusted_p_values: Vec<T>,
    pub bhy_haircut_sharpe_ratios: Vec<T>,
}
fn hlz_p<T: Real>(z: T) -> T {
    T::from_f64(2.) * (T::one() - z.abs().norm_cdf())
}
fn norm_inv<T: Real>(p: T) -> T {
    // Acklam's rational approximation followed by Newton refinement in T.
    let pp = p
        .max_val(T::from_f64(1e-15))
        .min_val(T::one() - T::from_f64(1e-15));
    let a = [
        -39.6968302866538,
        220.946098424521,
        -275.928510446969,
        138.357751867269,
        -30.6647980661472,
        2.50662827745924,
    ];
    let b = [
        -54.4760987982241,
        161.585836858041,
        -155.698979859887,
        66.8013118877197,
        -13.2806815528857,
    ];
    let c = [
        -0.00778489400243029,
        -0.322396458041136,
        -2.40075827716184,
        -2.54973253934373,
        4.37466414146497,
        2.93816398269878,
    ];
    let d = [
        0.00778469570904146,
        0.32246712907004,
        2.445134137143,
        3.75440866190742,
    ];
    let mut x;
    if pp.val() < 0.02425 {
        let q = (-T::from_f64(2.0) * pp.ln()).sqrt();
        x = (((((T::from_f64(c[0]) * q + T::from_f64(c[1])) * q + T::from_f64(c[2])) * q
            + T::from_f64(c[3]))
            * q
            + T::from_f64(c[4]))
            * q
            + T::from_f64(c[5]))
            / ((((T::from_f64(d[0]) * q + T::from_f64(d[1])) * q + T::from_f64(d[2])) * q
                + T::from_f64(d[3]))
                * q
                + T::one());
    } else if pp.val() > 0.97575 {
        let q = (-T::from_f64(2.0) * (T::one() - pp).ln()).sqrt();
        x = -(((((T::from_f64(c[0]) * q + T::from_f64(c[1])) * q + T::from_f64(c[2])) * q
            + T::from_f64(c[3]))
            * q
            + T::from_f64(c[4]))
            * q
            + T::from_f64(c[5]))
            / ((((T::from_f64(d[0]) * q + T::from_f64(d[1])) * q + T::from_f64(d[2])) * q
                + T::from_f64(d[3]))
                * q
                + T::one());
    } else {
        let q = pp - T::from_f64(0.5);
        let r = q * q;
        x = (((((T::from_f64(a[0]) * r + T::from_f64(a[1])) * r + T::from_f64(a[2])) * r
            + T::from_f64(a[3]))
            * r
            + T::from_f64(a[4]))
            * r
            + T::from_f64(a[5]))
            * q
            / (((((T::from_f64(b[0]) * r + T::from_f64(b[1])) * r + T::from_f64(b[2])) * r
                + T::from_f64(b[3]))
                * r
                + T::from_f64(b[4]))
                * r
                + T::one());
    }
    for _ in 0..2 {
        x -= (x.norm_cdf() - pp) / x.norm_pdf().max_val(T::from_f64(1e-15));
    }
    x
}
pub fn evaluate_harvey_liu_zhu_multiple_testing_and_mintrl<T: Real>(
    s: &HarveyLiuZhuMultipleTestingSharpeSpec<T>,
) -> HarveyLiuZhuMultipleTestingSharpeReport<T> {
    let t = T::from_f64(s.sample_length_t.max(1) as f64);
    let se = ((T::one() - s.return_skewness * s.focal_sharpe_ratio
        + ((s.return_kurtosis_pearson - T::from_f64(3.0)) * T::from_f64(0.25))
            * s.focal_sharpe_ratio
            * s.focal_sharpe_ratio)
        / t)
        .abs()
        .sqrt();
    let psr = ((s.focal_sharpe_ratio - s.benchmark_sharpe_sr0) / se.max_val(T::from_f64(1e-12)))
        .norm_cdf();
    let mintrl = ((s.confidence_z_alpha * se.max_val(T::from_f64(1e-12))
        / (s.focal_sharpe_ratio - s.benchmark_sharpe_sr0).max_val(T::from_f64(1e-8)))
    .powf(2.0)
        * t)
        .max_val(T::one());
    let k_len = s.candidate_t_stats.len().max(1);
    let m = s.num_independent_trials_m.max(k_len);
    let maxfalse = if m <= 1 {
        s.benchmark_sharpe_sr0
    } else {
        let gamma_em = std::f64::consts::EULER_GAMMA;
        let m_f = m as f64;
        let z1 = norm_inv(T::from_f64(1.0 - 1.0 / m_f));
        let z2 = norm_inv(T::from_f64(1.0 - 1.0 / (m_f * std::f64::consts::E)));
        s.benchmark_sharpe_sr0
            + s.cross_sectional_sharpe_std
                * (T::from_f64(1.0 - gamma_em) * z1 + T::from_f64(gamma_em) * z2)
    };
    let dsr = ((s.focal_sharpe_ratio - maxfalse) / se.max_val(T::from_f64(1e-12))).norm_cdf();
    let mut raw = Vec::with_capacity(s.candidate_t_stats.len());
    for z in &s.candidate_t_stats {
        raw.push(hlz_p(*z).min_val(T::one()).max_val(T::from_f64(1e-15)));
    }
    let n = raw.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        raw[a]
            .val()
            .partial_cmp(&raw[b].val())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let c_m: f64 = (1..=m).map(|j| 1.0 / j as f64).sum();
    let mut bon = vec![T::zero(); n];
    let mut holm = vec![T::zero(); n];
    let mut bhy = vec![T::zero(); n];
    let mut hair = vec![T::zero(); n];
    let mut running_holm = T::zero();
    for (rank, &idx) in order.iter().enumerate() {
        bon[idx] = (raw[idx] * T::from_f64(m as f64)).min_val(T::one());
        let cand = (raw[idx] * T::from_f64((m - rank) as f64)).min_val(T::one());
        running_holm = running_holm.max_val(cand);
        holm[idx] = running_holm;
    }
    let mut running_bhy = T::one();
    for (rank, &idx) in order.iter().enumerate().rev() {
        let denom = (rank + 1 + m - n) as f64;
        let cand = (raw[idx] * T::from_f64((m as f64) * c_m / denom)).min_val(T::one());
        running_bhy = running_bhy.min_val(cand);
        bhy[idx] = running_bhy;
        let u_prob = (T::one() - T::from_f64(0.5) * bhy[idx])
            .min_val(T::from_f64(1.0 - 1e-12))
            .max_val(T::from_f64(0.5 + 1e-12));
        let t_adj = norm_inv(u_prob);
        let t_orig = s.candidate_t_stats[idx].abs().max_val(T::from_f64(1e-8));
        hair[idx] = s.focal_sharpe_ratio * (t_adj / t_orig).min_val(T::one());
    }
    HarveyLiuZhuMultipleTestingSharpeReport {
        sharpe_standard_error: se,
        probabilistic_sharpe_ratio_psr: psr,
        minimum_track_record_length_mintrl: mintrl,
        expected_max_false_sharpe_sr0_star: maxfalse,
        deflated_sharpe_ratio_dsr: dsr,
        raw_two_sided_p_values: raw,
        bonferroni_adjusted_p_values: bon,
        holm_adjusted_p_values: holm,
        bhy_fdr_adjusted_p_values: bhy,
        bhy_haircut_sharpe_ratios: hair,
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

#[derive(Clone, Debug)]
pub struct CochraneHansenJagannathanSdfSpec<T: Real = f64> {
    pub mean_excess_returns: [T; 2],
    pub excess_return_cov: [[T; 2]; 2],
    pub risk_free_rate: T,
    pub candidate_sdf_pricing_errors: [T; 2],
    pub steady_state_log_dp_ratio: T,
    pub var1_companion_matrix: [[T; 2]; 2],
    pub state_innovation_cov: [[T; 2]; 2],
}
#[derive(Clone, Debug)]
pub struct CochraneHansenJagannathanSdfReport<T: Real = f64> {
    pub mean_sdf_m_bar: T,
    pub max_sharpe_ratio: T,
    pub hj_sdf_volatility_lower_bound: T,
    pub linear_sdf_loadings_b: [T; 2],
    pub hj_specification_distance: T,
    pub campbell_shiller_rho: T,
    pub discount_rate_news_variance_share: T,
    pub cash_flow_news_variance_share: T,
}
fn inv2_rm<T: Real>(a: [[T; 2]; 2]) -> [[T; 2]; 2] {
    let d = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let z = T::one() / (d + T::from_f64(1e-14));
    [[a[1][1] * z, -a[0][1] * z], [-a[1][0] * z, a[0][0] * z]]
}
pub fn evaluate_cochrane_hansen_jagannathan_sdf_and_campbell_shiller<T: Real>(
    s: &CochraneHansenJagannathanSdfSpec<T>,
) -> CochraneHansenJagannathanSdfReport<T> {
    let m = T::one() / (T::one() + s.risk_free_rate);
    let iv = inv2_rm(s.excess_return_cov);
    let q0 = iv[0][0] * s.mean_excess_returns[0] + iv[0][1] * s.mean_excess_returns[1];
    let q1 = iv[1][0] * s.mean_excess_returns[0] + iv[1][1] * s.mean_excess_returns[1];
    let sr = (s.mean_excess_returns[0] * q0 + s.mean_excess_returns[1] * q1).sqrt();
    let b = [m * q0, m * q1];
    let av = [
        [
            s.excess_return_cov[0][0] + s.mean_excess_returns[0] * s.mean_excess_returns[0],
            s.excess_return_cov[0][1] + s.mean_excess_returns[0] * s.mean_excess_returns[1],
        ],
        [
            s.excess_return_cov[1][0] + s.mean_excess_returns[1] * s.mean_excess_returns[0],
            s.excess_return_cov[1][1] + s.mean_excess_returns[1] * s.mean_excess_returns[1],
        ],
    ];
    let ai = inv2_rm(av);
    let d0 =
        ai[0][0] * s.candidate_sdf_pricing_errors[0] + ai[0][1] * s.candidate_sdf_pricing_errors[1];
    let d1 =
        ai[1][0] * s.candidate_sdf_pricing_errors[0] + ai[1][1] * s.candidate_sdf_pricing_errors[1];
    let dist =
        (s.candidate_sdf_pricing_errors[0] * d0 + s.candidate_sdf_pricing_errors[1] * d1).sqrt();
    let rho = T::one() / (T::one() + s.steady_state_log_dp_ratio.exp());
    let a = s.var1_companion_matrix;
    let mmat = [
        [T::one() - rho * a[0][0], -rho * a[0][1]],
        [-rho * a[1][0], T::one() - rho * a[1][1]],
    ];
    let mi = inv2_rm(mmat);
    let lam = [
        rho * (a[0][0] * mi[0][0] + a[0][1] * mi[1][0]),
        rho * (a[0][0] * mi[0][1] + a[0][1] * mi[1][1]),
    ];
    let cf = [T::one() + lam[0], lam[1]];
    let qform = |v: [T; 2]| -> T {
        v[0] * (s.state_innovation_cov[0][0] * v[0] + s.state_innovation_cov[0][1] * v[1])
            + v[1] * (s.state_innovation_cov[1][0] * v[0] + s.state_innovation_cov[1][1] * v[1])
    };
    let vdr = qform(lam).max_val(T::zero());
    let vcf = qform(cf).max_val(T::zero());
    let dr = vdr / (vdr + vcf).max_val(T::from_f64(1e-14));
    CochraneHansenJagannathanSdfReport {
        mean_sdf_m_bar: m,
        max_sharpe_ratio: sr,
        hj_sdf_volatility_lower_bound: m * sr,
        linear_sdf_loadings_b: b,
        hj_specification_distance: dist,
        campbell_shiller_rho: rho,
        discount_rate_news_variance_share: dr,
        cash_flow_news_variance_share: T::one() - dr,
    }
}

#[derive(Clone, Debug)]
pub struct WakkerTverskyKahnemanCptSpec<T: Real = f64> {
    pub scenario_returns: [T; 4],
    pub scenario_probabilities: [T; 4],
    pub reference_return: T,
    pub alpha_gain: T,
    pub beta_loss: T,
    pub loss_aversion_lambda: T,
    pub prelec_gamma_gain: T,
    pub prelec_gamma_loss: T,
    pub prelec_delta: T,
}
#[derive(Clone, Debug)]
pub struct WakkerTverskyKahnemanCptReport<T: Real = f64> {
    pub decision_weights: [T; 4],
    pub prospect_values: [T; 4],
    pub cumulative_prospect_utility: T,
    pub cpt_certainty_equivalent: T,
    pub loss_tail_overweighting_ratio: T,
    pub gain_tail_overweighting_ratio: T,
}
pub fn evaluate_wakker_tversky_kahneman_cumulative_prospect_utility<T: Real>(
    s: &WakkerTverskyKahnemanCptSpec<T>,
) -> WakkerTverskyKahnemanCptReport<T> {
    let mut idx = [0usize, 1, 2, 3];
    idx.sort_by(|&i, &j| {
        s.scenario_returns[i]
            .val()
            .partial_cmp(&s.scenario_returns[j].val())
            .unwrap()
    });
    let mut w = [T::zero(); 4];
    let mut v = [T::zero(); 4];
    let prelec = |p: T, gamma: T| {
        if p.val() <= 1e-15 {
            T::zero()
        } else if p.val() >= 1.0 - 1e-15 {
            T::one()
        } else {
            (-s.prelec_delta * (-p.ln()).pow(gamma)).exp()
        }
    };
    for (k, &i) in idx.iter().enumerate() {
        let gain = s.scenario_returns[i].val() >= s.reference_return.val();
        let gamma = if gain {
            s.prelec_gamma_gain
        } else {
            s.prelec_gamma_loss
        };
        let mut cum_prob = T::zero();
        let mut next_prob = T::zero();
        if gain {
            for &j in idx.iter().skip(k) {
                cum_prob += s.scenario_probabilities[j].max_val(T::from_f64(1e-14));
            }
            for &j in idx.iter().skip(k + 1) {
                next_prob += s.scenario_probabilities[j].max_val(T::from_f64(1e-14));
            }
        } else {
            for &j in idx.iter().take(k + 1) {
                cum_prob += s.scenario_probabilities[j].max_val(T::from_f64(1e-14));
            }
            for &j in idx.iter().take(k) {
                next_prob += s.scenario_probabilities[j].max_val(T::from_f64(1e-14));
            }
        }
        w[i] = prelec(cum_prob, gamma) - prelec(next_prob, gamma);
        if gain {
            let x = s.scenario_returns[i] - s.reference_return;
            v[i] = x.max_val(T::zero()).pow(s.alpha_gain);
        } else {
            let x = s.reference_return - s.scenario_returns[i];
            v[i] = -s.loss_aversion_lambda * x.max_val(T::zero()).pow(s.beta_loss);
        }
    }
    let cum = (0..4).fold(T::zero(), |z, i| z + w[i] * v[i]);
    let ce = if cum.val() >= 0.0 {
        cum.max_val(T::zero()).pow(T::one() / s.alpha_gain)
    } else {
        -((-cum) / s.loss_aversion_lambda).pow(T::one() / s.beta_loss)
    };
    let loss_ratio = w[idx[0]] / (s.scenario_probabilities[idx[0]] + T::from_f64(1e-14));
    let gain_ratio = w[idx[3]] / (s.scenario_probabilities[idx[3]] + T::from_f64(1e-14));
    WakkerTverskyKahnemanCptReport {
        decision_weights: w,
        prospect_values: v,
        cumulative_prospect_utility: cum,
        cpt_certainty_equivalent: s.reference_return + ce,
        loss_tail_overweighting_ratio: loss_ratio,
        gain_tail_overweighting_ratio: gain_ratio,
    }
}

#[derive(Clone, Debug)]
pub struct SornetteLpplBubbleCrashSpec<T: Real = f64> {
    pub observation_times: Vec<T>,
    pub log_prices: Vec<T>,
    pub critical_time_tc: T,
    pub power_law_exponent_m: T,
    pub log_periodic_omega: T,
    pub expected_crash_drop_kappa: T,
    pub ridge_penalty: T,
}
#[derive(Clone, Debug)]
pub struct SornetteLpplBubbleCrashReport<T: Real = f64> {
    pub linear_params_a_b_c1_c2: [T; 4],
    pub oscillation_amplitude_c: T,
    pub fitted_rmse: T,
    pub sornette_damping_parameter: T,
    pub terminal_crash_hazard_rate: T,
    pub cumulative_crash_probability: T,
}
fn solve4_lp<T: Real>(mut a: [[T; 5]; 4]) -> [T; 4] {
    for i in 0..4 {
        let mut p = i;
        for r in i + 1..4 {
            if a[r][i].abs().val() > a[p][i].abs().val() {
                p = r;
            }
        }
        a.swap(i, p);
        let mut d = a[i][i];
        if d.val().abs() < 1e-14 {
            d = T::from_f64(if d.val().is_sign_negative() {
                -1e-14
            } else {
                1e-14
            });
        }
        for value in a[i].iter_mut().skip(i) {
            *value /= d;
        }
        let pivot = a[i];
        for (r, row) in a.iter_mut().enumerate() {
            if r != i {
                let f = row[i];
                for (j, value) in row.iter_mut().enumerate().skip(i) {
                    *value -= f * pivot[j];
                }
            }
        }
    }
    [a[0][4], a[1][4], a[2][4], a[3][4]]
}
pub fn calibrate_sornette_lppl_bubble_and_crash_hazard<T: Real>(
    s: &SornetteLpplBubbleCrashSpec<T>,
) -> SornetteLpplBubbleCrashReport<T> {
    let n = s.observation_times.len().max(1);
    let mut g = [[T::zero(); 5]; 4];
    for (&t, &y) in s.observation_times.iter().zip(s.log_prices.iter()) {
        let dt = (s.critical_time_tc - t).max_val(T::from_f64(1e-12));
        let x = dt.pow(s.power_law_exponent_m);
        let th = s.log_periodic_omega * dt.ln();
        let z = [T::one(), x, x * th.cos(), x * th.sin()];
        for i in 0..4 {
            for j in 0..4 {
                g[i][j] += z[i] * z[j];
            }
            g[i][4] += z[i] * y;
        }
    }
    for (i, row) in g.iter_mut().enumerate() {
        row[i] += s.ridge_penalty;
    }
    let p = solve4_lp(g);
    let mut e = T::zero();
    for (&t, &y) in s.observation_times.iter().zip(s.log_prices.iter()) {
        let dt = (s.critical_time_tc - t).max_val(T::from_f64(1e-12));
        let x = dt.pow(s.power_law_exponent_m);
        let th = s.log_periodic_omega * dt.ln();
        let fit = p[0] + p[1] * x + p[2] * x * th.cos() + p[3] * x * th.sin();
        let d = fit - y;
        e += d * d;
    }
    let rm = (e / n as f64).sqrt();
    let c = (p[2] * p[2] + p[3] * p[3] + T::from_f64(1e-24)).sqrt();
    let damp = s.power_law_exponent_m.abs() * p[1].abs()
        / (s.log_periodic_omega.abs() * c + T::from_f64(1e-14));
    let tn = *s.observation_times.last().unwrap_or(&T::zero());
    let tau = (s.critical_time_tc - tn).max_val(T::from_f64(1e-8));
    let th = s.log_periodic_omega * tau.ln();
    let hp = -tau.pow(s.power_law_exponent_m - T::one())
        * (s.power_law_exponent_m * p[1]
            + (s.power_law_exponent_m * p[2] - s.log_periodic_omega * p[3]) * th.cos()
            + (s.power_law_exponent_m * p[3] + s.log_periodic_omega * p[2]) * th.sin());
    let terminal =
        (hp / s.expected_crash_drop_kappa.max_val(T::from_f64(1e-12))).max_val(T::zero());
    let eval_fit = |t_obs: T| -> T {
        let dt = (s.critical_time_tc - t_obs).max_val(T::from_f64(1e-8));
        let x = dt.pow(s.power_law_exponent_m);
        let th_obs = s.log_periodic_omega * dt.ln();
        p[0] + p[1] * x + p[2] * x * th_obs.cos() + p[3] * x * th_obs.sin()
    };
    let t0 = s.observation_times.first().copied().unwrap_or(T::zero());
    let drop = (eval_fit(tn) - eval_fit(t0)).max_val(T::zero());
    let prob = T::one() - (-drop / s.expected_crash_drop_kappa.max_val(T::from_f64(1e-12))).exp();
    SornetteLpplBubbleCrashReport {
        linear_params_a_b_c1_c2: p,
        oscillation_amplitude_c: c,
        fitted_rmse: rm,
        sornette_damping_parameter: damp,
        terminal_crash_hazard_rate: terminal,
        cumulative_crash_probability: prob,
    }
}

#[derive(Clone, Debug)]
pub struct FollmerSchiedGennaioliDiagnosticRiskSpec<T: Real = f64> {
    pub scenario_pnl: Vec<T>,
    pub baseline_probabilities: Vec<T>,
    pub entropic_ambiguity_theta: T,
    pub kl_divergence_budget_eta: T,
    pub diagnostic_overreaction_theta: T,
    pub boom_news_drift_shock: T,
    pub bust_news_drift_shock: T,
}

#[derive(Clone, Debug)]
pub struct FollmerSchiedGennaioliDiagnosticRiskReport<T: Real = f64> {
    pub rational_expected_loss: T,
    pub follmer_schied_entropic_risk_measure: T,
    pub worst_case_kl_divergence: T,
    pub ambiguity_penalty_premium: T,
    pub kl_ball_robust_superhedging_bound: T,
    pub diagnostic_boom_perceived_risk: T,
    pub diagnostic_bust_perceived_risk: T,
    pub gennaioli_shleifer_crisis_belief_reversal_spike: T,
    pub boom_tail_neglect_bias: T,
}

pub fn evaluate_follmer_schied_convex_risk_and_diagnostic_expectations<T: Real>(
    s: &FollmerSchiedGennaioliDiagnosticRiskSpec<T>,
) -> FollmerSchiedGennaioliDiagnosticRiskReport<T> {
    let n = s.scenario_pnl.len().min(s.baseline_probabilities.len());
    assert!(n > 0);
    let mut ps = T::zero();
    for i in 0..n {
        ps += s.baseline_probabilities[i];
    }
    let ps = ps.max_val(T::from_f64(1e-12));
    let mut losses = Vec::with_capacity(n);
    let mut exp_loss = T::zero();
    for i in 0..n {
        let l = -s.scenario_pnl[i];
        losses.push(l);
        exp_loss += s.baseline_probabilities[i] / ps * l;
    }
    let theta = s.entropic_ambiguity_theta.abs().max_val(T::from_f64(1e-8));
    let lmax = losses
        .iter()
        .copied()
        .fold(T::from_f64(-1e300), |a, x| a.max_val(x));
    let mut weights = Vec::with_capacity(n);
    let mut z = T::zero();
    for i in 0..n {
        let w = s.baseline_probabilities[i] / ps * ((losses[i] - lmax) / theta).exp();
        weights.push(w);
        z += w;
    }
    let entropic = lmax + theta * z.ln();
    let penalty = entropic - exp_loss;
    let mut kl = T::zero();
    let mut var = T::zero();
    for i in 0..n {
        let q = weights[i] / z.max_val(T::from_f64(1e-12));
        let p = s.baseline_probabilities[i] / ps;
        kl += q * (q / p.max_val(T::from_f64(1e-12))).ln();
        let d = losses[i] - exp_loss;
        var += p * d * d;
    }
    let sigma = (var + T::from_f64(1e-8)).sqrt();
    let mut boom_w = Vec::with_capacity(n);
    let mut bust_w = Vec::with_capacity(n);
    let mut zb = T::zero();
    let mut zs = T::zero();
    for i in 0..n {
        let p = s.baseline_probabilities[i] / ps;
        let d = (losses[i] - exp_loss) / sigma;
        let wb = p * (-s.diagnostic_overreaction_theta * s.boom_news_drift_shock.abs() * d).exp();
        let ws = p * (s.diagnostic_overreaction_theta * s.bust_news_drift_shock.abs() * d).exp();
        boom_w.push(wb);
        bust_w.push(ws);
        zb += wb;
        zs += ws;
    }
    let mut boom = T::zero();
    let mut bust = T::zero();
    for i in 0..n {
        boom += boom_w[i] / zb.max_val(T::from_f64(1e-12)) * losses[i];
        bust += bust_w[i] / zs.max_val(T::from_f64(1e-12)) * losses[i];
    }
    let robust = entropic + theta * s.kl_divergence_budget_eta.max_val(T::zero());
    FollmerSchiedGennaioliDiagnosticRiskReport {
        rational_expected_loss: exp_loss,
        follmer_schied_entropic_risk_measure: entropic,
        worst_case_kl_divergence: kl,
        ambiguity_penalty_premium: penalty,
        kl_ball_robust_superhedging_bound: robust,
        diagnostic_boom_perceived_risk: boom,
        diagnostic_bust_perceived_risk: bust,
        gennaioli_shleifer_crisis_belief_reversal_spike: bust - boom,
        boom_tail_neglect_bias: exp_loss - boom,
    }
}

#[cfg(test)]
mod domain_suite_28_tests {
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
mod domain_suite_30_tests {
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
mod domain_suite_31_tests {
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

#[cfg(test)]
mod domain_suite_32_tests {
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

#[cfg(test)]
mod domain_suite_36_tests {
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
    fn pot_gpd_suite36() {
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
    fn brinson_carino_suite36() {
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

#[cfg(test)]
mod domain_suite_37_tests {
    use super::*;
    use ql_ad::Dual;
    fn hj<T: Real>(x: T) -> CochraneHansenJagannathanSdfSpec<T> {
        CochraneHansenJagannathanSdfSpec {
            mean_excess_returns: [x, T::from_f64(0.05)],
            excess_return_cov: [
                [T::from_f64(0.04), T::from_f64(0.01)],
                [T::from_f64(0.01), T::from_f64(0.06)],
            ],
            risk_free_rate: T::from_f64(0.02),
            candidate_sdf_pricing_errors: [T::from_f64(0.01), T::from_f64(-0.01)],
            steady_state_log_dp_ratio: T::from_f64(1.0),
            var1_companion_matrix: [[T::from_f64(0.9), T::zero()], [T::zero(), T::from_f64(0.8)]],
            state_innovation_cov: [
                [T::from_f64(0.02), T::zero()],
                [T::zero(), T::from_f64(0.03)],
            ],
        }
    }
    #[test]
    fn hj_dual_fd_and_share_partition() {
        let r = evaluate_cochrane_hansen_jagannathan_sdf_and_campbell_shiller(&hj(0.04));
        assert!(
            (r.discount_rate_news_variance_share + r.cash_flow_news_variance_share - 1.0).abs()
                < 1e-12
        );
        let h = 1e-5;
        let fd = (evaluate_cochrane_hansen_jagannathan_sdf_and_campbell_shiller(&hj(0.04 + h))
            .max_sharpe_ratio
            .val()
            - evaluate_cochrane_hansen_jagannathan_sdf_and_campbell_shiller(&hj(0.04 - h))
                .max_sharpe_ratio
                .val())
            / (2.0 * h);
        let ad = evaluate_cochrane_hansen_jagannathan_sdf_and_campbell_shiller(&hj(
            Dual::<1>::variable(0.04, 0),
        ))
        .max_sharpe_ratio;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
    fn cpt_spec<T: Real>(l: T, mid: T) -> WakkerTverskyKahnemanCptSpec<T> {
        WakkerTverskyKahnemanCptSpec {
            scenario_returns: [
                T::from_f64(-0.2),
                T::from_f64(-0.05),
                mid,
                T::from_f64(0.15),
            ],
            scenario_probabilities: [
                T::from_f64(0.1),
                T::from_f64(0.4),
                T::from_f64(0.4),
                T::from_f64(0.1),
            ],
            reference_return: T::zero(),
            alpha_gain: T::from_f64(0.88),
            beta_loss: T::from_f64(0.88),
            loss_aversion_lambda: l,
            prelec_gamma_gain: T::from_f64(0.7),
            prelec_gamma_loss: T::from_f64(0.7),
            prelec_delta: T::one(),
        }
    }
    #[test]
    fn cpt_lambda_dual_fd() {
        let r = evaluate_wakker_tversky_kahneman_cumulative_prospect_utility(&cpt_spec(2.25, 0.03));
        assert!(r.loss_tail_overweighting_ratio > 1.0);
        assert!(r.gain_tail_overweighting_ratio > 1.0);
        let h = 1e-5;
        let fd = (evaluate_wakker_tversky_kahneman_cumulative_prospect_utility(&cpt_spec(
            2.25,
            0.03 + h,
        ))
        .cumulative_prospect_utility
            - evaluate_wakker_tversky_kahneman_cumulative_prospect_utility(&cpt_spec(
                2.25,
                0.03 - h,
            ))
            .cumulative_prospect_utility)
            / (2.0 * h);
        let ad = evaluate_wakker_tversky_kahneman_cumulative_prospect_utility(&cpt_spec(
            Dual::<1>::from_value(2.25),
            Dual::<1>::variable(0.03, 0),
        ))
        .cumulative_prospect_utility;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
    fn lppl_noiseless<T: Real>(tc: T) -> SornetteLpplBubbleCrashSpec<T> {
        let times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let log_prices = times
            .iter()
            .map(|&t| {
                let dt = T::from_f64(10.0 - t);
                let x = dt.pow(T::from_f64(0.5));
                let th = T::from_f64(6.0) * dt.ln();
                T::from_f64(5.0) - T::from_f64(0.4) * x + T::from_f64(0.02) * x * th.cos()
                    - T::from_f64(0.01) * x * th.sin()
            })
            .collect();
        SornetteLpplBubbleCrashSpec {
            observation_times: times.iter().map(|&t| T::from_f64(t)).collect(),
            log_prices,
            critical_time_tc: tc,
            power_law_exponent_m: T::from_f64(0.5),
            log_periodic_omega: T::from_f64(6.0),
            expected_crash_drop_kappa: T::from_f64(0.2),
            ridge_penalty: T::zero(),
        }
    }
    #[test]
    fn lppl_dual_fd() {
        let h = 1e-5;
        let fd = (calibrate_sornette_lppl_bubble_and_crash_hazard(&lppl_noiseless(10.0 + h))
            .terminal_crash_hazard_rate
            - calibrate_sornette_lppl_bubble_and_crash_hazard(&lppl_noiseless(10.0 - h))
                .terminal_crash_hazard_rate)
            / (2.0 * h);
        let r = calibrate_sornette_lppl_bubble_and_crash_hazard(&lppl_noiseless(10.0));
        assert!(r.fitted_rmse < 1e-10);
        assert!(r.linear_params_a_b_c1_c2[1] < 0.0);
        assert!(r.sornette_damping_parameter >= 1.0);
        let ad = calibrate_sornette_lppl_bubble_and_crash_hazard(&lppl_noiseless(
            Dual::<1>::variable(10.0, 0),
        ))
        .terminal_crash_hazard_rate;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_38_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec(sr: f64) -> HarveyLiuZhuMultipleTestingSharpeSpec<f64> {
        HarveyLiuZhuMultipleTestingSharpeSpec {
            focal_sharpe_ratio: sr,
            benchmark_sharpe_sr0: 0.0,
            sample_length_t: 252,
            return_skewness: -0.5,
            return_kurtosis_pearson: 4.0,
            num_independent_trials_m: 1,
            cross_sectional_sharpe_std: 0.3,
            confidence_z_alpha: 1.96,
            candidate_t_stats: vec![1.5],
        }
    }
    #[test]
    fn hlz_invariants_and_dual_fd() {
        let r = evaluate_harvey_liu_zhu_multiple_testing_and_mintrl(&spec(0.8));
        assert!(r.raw_two_sided_p_values[0] <= r.bhy_fdr_adjusted_p_values[0]);
        assert!((r.bhy_haircut_sharpe_ratios[0] - 0.8).abs() < 1e-12);
        let h = 1e-5;
        let fd = (evaluate_harvey_liu_zhu_multiple_testing_and_mintrl(&spec(0.8 + h))
            .probabilistic_sharpe_ratio_psr
            - evaluate_harvey_liu_zhu_multiple_testing_and_mintrl(&spec(0.8 - h))
                .probabilistic_sharpe_ratio_psr)
            / (2.0 * h);
        let sd = spec(0.8);
        let ad = evaluate_harvey_liu_zhu_multiple_testing_and_mintrl(
            &HarveyLiuZhuMultipleTestingSharpeSpec {
                focal_sharpe_ratio: Dual::<1>::variable(0.8, 0),
                benchmark_sharpe_sr0: Dual::from_value(0.0),
                sample_length_t: sd.sample_length_t,
                return_skewness: Dual::from_value(sd.return_skewness),
                return_kurtosis_pearson: Dual::from_value(sd.return_kurtosis_pearson),
                num_independent_trials_m: sd.num_independent_trials_m,
                cross_sectional_sharpe_std: Dual::from_value(sd.cross_sectional_sharpe_std),
                confidence_z_alpha: Dual::from_value(sd.confidence_z_alpha),
                candidate_t_stats: vec![Dual::from_value(1.5)],
            },
        )
        .probabilistic_sharpe_ratio_psr;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_40_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(x: T) -> FollmerSchiedGennaioliDiagnosticRiskSpec<T> {
        FollmerSchiedGennaioliDiagnosticRiskSpec {
            scenario_pnl: vec![T::from_f64(1.0) + x, T::from_f64(-2.0), T::from_f64(-5.0)],
            baseline_probabilities: vec![T::from_f64(0.5), T::from_f64(0.3), T::from_f64(0.2)],
            entropic_ambiguity_theta: T::from_f64(2.0),
            kl_divergence_budget_eta: T::from_f64(0.1),
            diagnostic_overreaction_theta: T::from_f64(1.5),
            boom_news_drift_shock: T::from_f64(0.2),
            bust_news_drift_shock: T::from_f64(1.0),
        }
    }
    #[test]
    fn risk_invariants_and_dual() {
        let r = evaluate_follmer_schied_convex_risk_and_diagnostic_expectations(&mk(0.0));
        assert!(r.follmer_schied_entropic_risk_measure >= r.rational_expected_loss);
        assert!(r.diagnostic_bust_perceived_risk > r.diagnostic_boom_perceived_risk);
        assert!(r.gennaioli_shleifer_crisis_belief_reversal_spike > 0.0);
        let h = 1e-5;
        let fd = (evaluate_follmer_schied_convex_risk_and_diagnostic_expectations(&mk(h))
            .follmer_schied_entropic_risk_measure
            - evaluate_follmer_schied_convex_risk_and_diagnostic_expectations(&mk(-h))
                .follmer_schied_entropic_risk_measure)
            / (2.0 * h);
        let ad = evaluate_follmer_schied_convex_risk_and_diagnostic_expectations(&mk(
            Dual::<1>::variable(0.0, 0),
        ))
        .follmer_schied_entropic_risk_measure;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

mod evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es_module {
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
            * (loglik(hit)
                - (n * p.ln() + (t - n) * (T::one() - p).max_val(T::from_f64(1e-12)).ln()));
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
        let bound = (payoff_max - payoff_min)
            * (T::from_f64(2.0) * T::from_f64(3.0).ln() / n_periods).sqrt();
        HartMasColellRegretMatchingPortfolioReport {
            cumulative_counterfactual_regrets: reg,
            positive_regret_potentials: pos,
            next_period_regret_matching_weights: w,
            weight_simplex_sum_error: w[0] + w[1] + w[2] - T::one(),
            blackwell_approachability_distance: (pos[0] * pos[0]
                + pos[1] * pos[1]
                + pos[2] * pos[2])
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
            campisi_identity_conservation_error: income + treasury + twist + spread + residual
                - active,
            karnosky_singer_forward_hedge_carry: carry,
            karnosky_singer_unhedged_currency_surprise: currency,
        }
    }

    #[cfg(test)]
    mod domain_suite_42_tests {
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
        fn suite42_portfolio_invariants_and_dual_fd() {
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
            assert!(
                p.capital_gains_overhang_cgo > 0.0 && p.lottery_skewness_overvaluation_drag < 0.0
            );
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
            let pad: BarberisJinWangProspectAnomalySpec<Dual<1>> =
                BarberisJinWangProspectAnomalySpec {
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
}
pub use evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es_module::*;

mod evaluate_embrechts_hill_pickands_and_subexponential_tail_module {
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
        let multiplier = m.pow(xi.max_val(T::zero()));
        let excess = multiplier / m.sqrt();
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
                + persistent / T::from_f64(24.0)
                    * (s.log_moneyness * s.log_moneyness / v - T::one()));
        let naive = s.confidence_z * v.sqrt();
        let z = s.confidence_z;
        let loss_k3 = -k3;
        let cf = z
            + loss_k3 / T::from_f64(6.0) * (z * z - T::one())
            + persistent / T::from_f64(24.0) * (z * z * z - T::from_f64(3.0) * z)
            - loss_k3 * loss_k3 / T::from_f64(36.0)
                * (T::from_f64(2.0) * z * z * z - T::from_f64(5.0) * z);
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
            sqrt_time_var_underestimation_ratio: non_gaussian
                / naive.abs().max_val(T::from_f64(1e-12)),
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
            ccy +=
                (s.portfolio_currency_weights[i] - s.benchmark_currency_weights[i]) * (ci - ccy_b);
        }
        let mut portfolio_total = T::zero();
        let mut benchmark_total = T::zero();
        for i in 0..n {
            let rp = s.local_portfolio_returns[i] - s.local_eurodollar_rates[i];
            let rb = s.local_benchmark_returns[i] - s.local_eurodollar_rates[i];
            let ci = s.spot_fx_returns[i] + s.local_eurodollar_rates[i] - s.base_eurodollar_rate;
            portfolio_total +=
                s.portfolio_asset_weights[i] * rp + s.portfolio_currency_weights[i] * ci;
            benchmark_total +=
                s.benchmark_asset_weights[i] * rb + s.benchmark_currency_weights[i] * ci;
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
        let n_mp = s
            .multi_period_port_returns
            .len()
            .min(s.multi_period_bench_returns.len());
        let periods = T::from_f64(n_mp.max(1) as f64);
        let diff_root = rp_cum.pow(T::one() / periods) - rb_cum.pow(T::one() / periods);
        let diff_root_safe = if diff_root >= T::zero() {
            diff_root.max_val(T::from_f64(1e-12))
        } else {
            diff_root.min_val(T::from_f64(-1e-12))
        };
        let multiplier = (cumulative / periods) / diff_root_safe;
        let mut sum_a = T::zero();
        let mut sum_a2 = T::zero();
        for i in 0..n_mp {
            let a_t = s.multi_period_port_returns[i] - s.multi_period_bench_returns[i];
            sum_a += a_t;
            sum_a2 += a_t * a_t;
        }
        let c_scale = (cumulative - multiplier * sum_a) / sum_a2.max_val(T::from_f64(1e-12));
        let mut linked = T::zero();
        for i in 0..n_mp {
            let a_t = s.multi_period_port_returns[i] - s.multi_period_bench_returns[i];
            let alpha_t = c_scale * a_t;
            linked += (multiplier + alpha_t) * a_t;
        }
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
        let n1_star = T::one()
            / (T::one() + (s.intensity_of_choice_beta * s.fundamentalist_info_cost_c1).exp());
        let eig = ((T::one() - n1_star) * s.chartist_trend_extrapolation_g
            + n1_star * s.fundamentalist_damping_v1)
            / r;
        let ratio = ((r - s.fundamentalist_damping_v1)
            / (s.chartist_trend_extrapolation_g - r)
                .abs()
                .max_val(T::from_f64(1e-12)))
        .abs()
        .max_val(T::from_f64(1e-12));
        let critical = ratio.ln().abs()
            / s.fundamentalist_info_cost_c1
                .abs()
                .max_val(T::from_f64(1e-12));
        let n_shocks = T::from_f64(s.shock_innovations.len().max(1) as f64);
        let traj_rms = (path.iter().skip(1).fold(T::zero(), |acc, &v| acc + v * v) / n_shocks)
            .max_val(T::from_f64(1e-24))
            .sqrt();
        let shock_rms = (s
            .shock_innovations
            .iter()
            .fold(T::zero(), |acc, &v| acc + v * v)
            / n_shocks)
            .max_val(T::from_f64(1e-24))
            .sqrt();
        BrockHommesAdaptiveBeliefReport {
            terminal_fundamentalist_fraction: T::one() - frac.last().copied().unwrap_or(T::zero()),
            terminal_chartist_fraction: frac.last().copied().unwrap_or(T::zero()),
            steady_state_jacobian_eigenvalue: eig,
            critical_bifurcation_beta: critical,
            terminal_price_deviation: x,
            price_deviation_trajectory: path,
            chartist_fraction_trajectory: frac,
            endogenous_excess_volatility_ratio: traj_rms / shock_rms.max_val(T::from_f64(1e-12)),
        }
    }

    #[cfg(test)]
    mod domain_suite_43_tests {
        use super::*;
        use ql_ad::Dual;

        fn hill<T: Real>(x: T) -> EmbrechtsHillPickandsSubexponentialSpec<T> {
            EmbrechtsHillPickandsSubexponentialSpec {
                loss_samples: vec![
                    T::from_f64(160.0) + x,
                    T::from_f64(110.0),
                    T::from_f64(82.0),
                    T::from_f64(66.0),
                    T::from_f64(56.0),
                    T::from_f64(50.0),
                    T::from_f64(46.0),
                    T::from_f64(43.0),
                    T::from_f64(40.0),
                    T::from_f64(38.0),
                    T::from_f64(36.0),
                    T::from_f64(35.0),
                ],
                tail_order_k: 2,
                extreme_tail_prob: T::from_f64(0.01),
                portfolio_positions_count: 8,
            }
        }

        #[test]
        fn suite43_invariants_and_dual_fd() {
            let hrep = evaluate_embrechts_hill_pickands_and_subexponential_tail(&hill(0.0));
            assert!(hrep.hill_xi_tail_index > 0.0 && hrep.pickands_xi_tail_index > 0.0);
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
            fn bp_spec<T: Real>(x: T) -> BouchaudPottersCumulantSmileSpec<T> {
                BouchaudPottersCumulantSmileSpec {
                    single_period_vol: T::from_f64(0.2) + x,
                    single_period_skewness: T::from_f64(-0.2),
                    single_period_excess_kurtosis: T::from_f64(3.0),
                    vol_of_vol_kurtosis_component: T::from_f64(2.0),
                    vol_decorrelation_steps: T::from_f64(10.0),
                    horizon_steps_n: T::from_f64(20.0),
                    log_moneyness: T::from_f64(0.1),
                    confidence_z: T::from_f64(2.33),
                }
            }
            let fd_bp =
                (evaluate_bouchaud_potters_power_law_cumulant_and_smile_scaling(&bp_spec(h))
                    .bouchaud_potters_non_gaussian_var
                    - evaluate_bouchaud_potters_power_law_cumulant_and_smile_scaling(&bp_spec(-h))
                        .bouchaud_potters_non_gaussian_var)
                    / (2.0 * h);
            let ad_bp = evaluate_bouchaud_potters_power_law_cumulant_and_smile_scaling(&bp_spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .bouchaud_potters_non_gaussian_var;
            assert!(ad_bp.d[0].abs() > 1e-6 && (ad_bp.d[0] - fd_bp).abs() < 1e-5);
            fn ks_spec<T: Real>(x: T) -> KarnoskySingerMencheroAttributionSpec<T> {
                KarnoskySingerMencheroAttributionSpec {
                    portfolio_asset_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
                    benchmark_asset_weights: vec![T::from_f64(0.5), T::from_f64(0.5)],
                    portfolio_currency_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
                    benchmark_currency_weights: vec![T::from_f64(0.5), T::from_f64(0.5)],
                    local_portfolio_returns: vec![T::from_f64(0.08) + x, T::from_f64(0.03)],
                    local_benchmark_returns: vec![T::from_f64(0.06), T::from_f64(0.04)],
                    local_eurodollar_rates: vec![T::from_f64(0.02), T::from_f64(0.025)],
                    base_eurodollar_rate: T::from_f64(0.03),
                    spot_fx_returns: vec![T::from_f64(0.01), T::from_f64(-0.005)],
                    multi_period_port_returns: vec![T::from_f64(0.02), T::from_f64(0.03)],
                    multi_period_bench_returns: vec![T::from_f64(0.018), T::from_f64(0.025)],
                }
            }
            let fd_ks =
                (evaluate_karnosky_singer_multicurrency_and_menchero_attribution(&ks_spec(h))
                    .karnosky_singer_security_selection
                    - evaluate_karnosky_singer_multicurrency_and_menchero_attribution(&ks_spec(
                        -h,
                    ))
                    .karnosky_singer_security_selection)
                    / (2.0 * h);
            let ad_ks = evaluate_karnosky_singer_multicurrency_and_menchero_attribution(&ks_spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .karnosky_singer_security_selection;
            assert!(ad_ks.d[0].abs() > 1e-6 && (ad_ks.d[0] - fd_ks).abs() < 1e-5);
            let d = evaluate_dgtw_characteristic_and_treynor_mazuy_timing(
                &DgtwTreynorMazuyTimingSpec {
                    current_weights: vec![0.6, 0.4],
                    lagged_weights: vec![0.5, 0.5],
                    stock_returns: vec![0.08, 0.03],
                    current_char_bench_returns: vec![0.06, 0.04],
                    lagged_char_bench_returns: vec![0.05, 0.04],
                    portfolio_excess_series: vec![0.01, 0.02, 0.04, 0.08],
                    market_excess_series: vec![0.0, 0.02, 0.04, 0.08],
                },
            );
            assert!(
                d.dgtw_decomposition_identity_error < 1e-12 && d.treynor_mazuy_timing_gamma > 0.0
            );
            fn tm_spec<T: Real>(x: T) -> DgtwTreynorMazuyTimingSpec<T> {
                DgtwTreynorMazuyTimingSpec {
                    current_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
                    lagged_weights: vec![T::from_f64(0.5), T::from_f64(0.5)],
                    stock_returns: vec![T::from_f64(0.08), T::from_f64(0.03)],
                    current_char_bench_returns: vec![T::from_f64(0.06), T::from_f64(0.04)],
                    lagged_char_bench_returns: vec![T::from_f64(0.05), T::from_f64(0.04)],
                    portfolio_excess_series: vec![
                        T::from_f64(0.01),
                        T::from_f64(0.02),
                        T::from_f64(0.04),
                        T::from_f64(0.08) + x,
                    ],
                    market_excess_series: vec![
                        T::zero(),
                        T::from_f64(0.02),
                        T::from_f64(0.04),
                        T::from_f64(0.08),
                    ],
                }
            }
            let fd_tm = (evaluate_dgtw_characteristic_and_treynor_mazuy_timing(&tm_spec(h))
                .treynor_mazuy_timing_gamma
                - evaluate_dgtw_characteristic_and_treynor_mazuy_timing(&tm_spec(-h))
                    .treynor_mazuy_timing_gamma)
                / (2.0 * h);
            let ad_tm = evaluate_dgtw_characteristic_and_treynor_mazuy_timing(&tm_spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .treynor_mazuy_timing_gamma;
            assert!(ad_tm.d[0].abs() > 1e-6 && (ad_tm.d[0] - fd_tm).abs() < 1e-5);
            let bh =
                evaluate_brock_hommes_adaptive_belief_switching(&BrockHommesAdaptiveBeliefSpec {
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
            fn bh_spec<T: Real>(x: T) -> BrockHommesAdaptiveBeliefSpec<T> {
                BrockHommesAdaptiveBeliefSpec {
                    initial_price_deviation: T::from_f64(0.1) + x,
                    risk_free_rate: T::from_f64(0.03),
                    fundamentalist_damping_v1: T::from_f64(0.8),
                    chartist_trend_extrapolation_g: T::from_f64(1.2),
                    fundamentalist_info_cost_c1: T::from_f64(0.1),
                    intensity_of_choice_beta: T::from_f64(8.0),
                    fitness_memory_eta: T::from_f64(0.8),
                    risk_aversion_sigma2: T::one(),
                    shock_innovations: vec![
                        T::from_f64(0.01),
                        T::from_f64(-0.005),
                        T::from_f64(0.003),
                    ],
                }
            }
            let fd_bh = (evaluate_brock_hommes_adaptive_belief_switching(&bh_spec(h))
                .terminal_price_deviation
                - evaluate_brock_hommes_adaptive_belief_switching(&bh_spec(-h))
                    .terminal_price_deviation)
                / (2.0 * h);
            let ad_bh = evaluate_brock_hommes_adaptive_belief_switching(&bh_spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .terminal_price_deviation;
            assert!(ad_bh.d[0].abs() > 1e-6 && (ad_bh.d[0] - fd_bh).abs() < 1e-5);
        }
    }
}
pub use evaluate_embrechts_hill_pickands_and_subexponential_tail_module::*;

mod evaluate_ohara_glosten_milgrom_sequential_spread_module {
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
            p = (T::one() + flow) / T::from_f64(2.0) * pp
                + (T::one() - flow) / T::from_f64(2.0) * pm;
            q -= flow;
        }
        let pb_t = (p * buy_h + (T::one() - p) * buy_l).max_val(T::from_f64(1e-12));
        let ps_t = (T::one() - pb_t).max_val(T::from_f64(1e-12));
        let ask = s.low_value_vl + (p * buy_h / pb_t) * d;
        let bid = s.low_value_vl + (p * buy_l / ps_t) * d;
        let mid =
            s.low_value_vl + p * d - s.dealer_inventory_gamma * d * d * p * (T::one() - p) * q;
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
            *x = (*x + suml / T::from_f64(n as f64)) * T::from_f64(0.5);
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
            trace_conservation_error: {
                let t_err = shr.iter().copied().fold(T::zero(), |a, x| a + x) - suml;
                (t_err * t_err + T::from_f64(1e-30)).sqrt()
            },
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
            euler_additivity_error: {
                let e_err = es - sum_rc;
                (e_err * e_err + T::from_f64(1e-30)).sqrt()
            },
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
            omega_put_call_parity_error: {
                let o_err = up - dn - (mean - s.mar_threshold);
                (o_err * o_err + T::from_f64(1e-30)).sqrt()
            },
            downside_deviation: (l2 / nf).sqrt(),
            sortino_ratio: (mean - s.mar_threshold) / (l2 / nf).sqrt(),
            kappa_three_ratio: (mean - s.mar_threshold)
                / (l3 / nf).max_val(T::from_f64(1e-18)).powf(1.0 / 3.0),
            ulcer_drawdown_index: ui,
            martin_ulcer_ratio: ex / ui,
            burke_drawdown_ratio: ex / dds.sqrt(),
        }
    }
    #[cfg(test)]
    mod domain_suite_45_tests {
        use super::*;
        use ql_ad::Dual;
        fn gm_spec<T: Real>(x: T) -> OharaGlostenMilgromSequentialSpec<T> {
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
        fn rie_spec<T: Real>(x: T) -> PottersBouchaudStieltjesRieSpec<T> {
            PottersBouchaudStieltjesRieSpec {
                sample_eigenvalues: vec![
                    T::from_f64(0.4) + x,
                    T::from_f64(0.9),
                    T::from_f64(1.5),
                    T::from_f64(2.4),
                ],
                aspect_ratio_q: T::from_f64(0.5),
                stieltjes_bandwidth_eta: T::from_f64(0.25),
            }
        }
        fn egarch_spec<T: Real>(x: T) -> TsayEgarchHamiltonRegimeSpec<T> {
            TsayEgarchHamiltonRegimeSpec {
                returns: vec![
                    T::from_f64(0.01) + x,
                    T::from_f64(-0.025),
                    T::from_f64(0.015),
                ],
                initial_variance: T::from_f64(0.04),
                egarch_omega: T::from_f64(-0.1),
                egarch_alpha: T::from_f64(0.15),
                egarch_leverage_gamma: T::from_f64(-0.08),
                egarch_beta: T::from_f64(0.92),
                calm_mean: T::zero(),
                calm_vol: T::from_f64(0.12),
                crisis_mean: T::from_f64(-0.01),
                crisis_vol: T::from_f64(0.30),
                stay_calm_p00: T::from_f64(0.95),
                stay_crisis_p11: T::from_f64(0.90),
            }
        }
        fn copula_spec<T: Real>(x: T) -> McneilFreyCopulaEulerAllocationSpec<T> {
            McneilFreyCopulaEulerAllocationSpec {
                clayton_theta: T::from_f64(2.0) + x,
                gumbel_theta: T::from_f64(1.5),
                student_t_df_nu: T::from_f64(5.0),
                student_t_corr_rho: T::from_f64(0.3),
                desk_weights: [
                    T::from_f64(0.25) + x,
                    T::from_f64(0.25),
                    T::from_f64(0.25),
                    T::from_f64(0.25),
                ],
                desk_expected_losses: [T::from_f64(0.01); 4],
                loss_covariance_4x4: [
                    [T::from_f64(0.04), T::from_f64(0.01), T::zero(), T::zero()],
                    [
                        T::from_f64(0.01),
                        T::from_f64(0.05),
                        T::from_f64(0.01),
                        T::zero(),
                    ],
                    [
                        T::zero(),
                        T::from_f64(0.01),
                        T::from_f64(0.06),
                        T::from_f64(0.01),
                    ],
                    [T::zero(), T::zero(), T::from_f64(0.01), T::from_f64(0.03)],
                ],
                es_tail_multiplier: T::from_f64(2.0),
            }
        }
        fn bacon_spec<T: Real>(x: T) -> BaconOmegaSortinoBurkeDrawdownSpec<T> {
            BaconOmegaSortinoBurkeDrawdownSpec {
                period_returns: vec![
                    T::from_f64(0.015) + x,
                    T::from_f64(-0.008),
                    T::from_f64(0.022),
                    T::from_f64(-0.004),
                ],
                mar_threshold: T::zero(),
                risk_free_rate: T::from_f64(0.001),
            }
        }
        #[test]
        fn suite45_portfolio_all_five_engines_invariants_and_dual_fd() {
            let h = 1e-6;
            let gm0 = evaluate_ohara_glosten_milgrom_sequential_spread(&gm_spec(0.));
            assert!(
                gm0.initial_adverse_selection_spread > 0.0 && gm0.terminal_posterior_belief_p > 0.5
            );
            let fd_gm = (evaluate_ohara_glosten_milgrom_sequential_spread(&gm_spec(h))
                .inventory_skewed_reservation_mid
                - evaluate_ohara_glosten_milgrom_sequential_spread(&gm_spec(-h))
                    .inventory_skewed_reservation_mid)
                / (2. * h);
            let ad_gm = evaluate_ohara_glosten_milgrom_sequential_spread(&gm_spec(
                Dual::<1>::variable(0., 0),
            ))
            .inventory_skewed_reservation_mid;
            assert!((ad_gm.d[0] - fd_gm).abs() < 1e-5);

            let rie0 = evaluate_potters_bouchaud_stieltjes_rie_shrinkage(&rie_spec(0.));
            assert!(
                rie0.trace_conservation_error < 1e-12
                    && rie0.rie_condition_number < rie0.sample_condition_number
                    && rie0.inverse_variance_portfolio_risk_reduction > 0.0
            );
            let fd_rie = (evaluate_potters_bouchaud_stieltjes_rie_shrinkage(&rie_spec(h))
                .rie_condition_number
                - evaluate_potters_bouchaud_stieltjes_rie_shrinkage(&rie_spec(-h))
                    .rie_condition_number)
                / (2. * h);
            let ad_rie = evaluate_potters_bouchaud_stieltjes_rie_shrinkage(&rie_spec(
                Dual::<1>::variable(0., 0),
            ))
            .rie_condition_number;
            assert!((ad_rie.d[0] - fd_rie).abs() < 1e-5);

            let eg0 = evaluate_tsay_egarch_and_hamilton_markov_switching(&egarch_spec(0.));
            assert!(
                eg0.egarch_negative_shock_asymmetry_ratio > 1.0
                    && eg0.terminal_filtered_crisis_probability > 0.0
                    && eg0.terminal_filtered_crisis_probability < 1.0
            );
            let fd_eg = (evaluate_tsay_egarch_and_hamilton_markov_switching(&egarch_spec(h))
                .terminal_egarch_volatility
                - evaluate_tsay_egarch_and_hamilton_markov_switching(&egarch_spec(-h))
                    .terminal_egarch_volatility)
                / (2. * h);
            let ad_eg = evaluate_tsay_egarch_and_hamilton_markov_switching(&egarch_spec(
                Dual::<1>::variable(0., 0),
            ))
            .terminal_egarch_volatility;
            assert!((ad_eg.d[0] - fd_eg).abs() < 1e-5);

            let cop0 = evaluate_mcneil_frey_copula_tail_and_euler_allocation(&copula_spec(0.));
            assert!(
                cop0.euler_additivity_error < 1e-12
                    && cop0.portfolio_diversification_benefit > 0.0
                    && cop0.clayton_lower_tail_lambda > 0.0
                    && cop0.clayton_lower_tail_lambda < 1.0
            );
            let fd_cop = (evaluate_mcneil_frey_copula_tail_and_euler_allocation(&copula_spec(h))
                .total_portfolio_es
                - evaluate_mcneil_frey_copula_tail_and_euler_allocation(&copula_spec(-h))
                    .total_portfolio_es)
                / (2. * h);
            let ad_cop = evaluate_mcneil_frey_copula_tail_and_euler_allocation(&copula_spec(
                Dual::<1>::variable(0., 0),
            ))
            .total_portfolio_es;
            assert!((ad_cop.d[0] - fd_cop).abs() < 1e-5);

            let bac0 = evaluate_bacon_omega_sortino_and_burke_drawdown(&bacon_spec(0.));
            assert!(
                bac0.keating_shadwick_omega_ratio > 1.0 && bac0.omega_put_call_parity_error < 1e-12
            );
            let fd_bac = (evaluate_bacon_omega_sortino_and_burke_drawdown(&bacon_spec(h))
                .keating_shadwick_omega_ratio
                - evaluate_bacon_omega_sortino_and_burke_drawdown(&bacon_spec(-h))
                    .keating_shadwick_omega_ratio)
                / (2. * h);
            let ad_bac = evaluate_bacon_omega_sortino_and_burke_drawdown(&bacon_spec(
                Dual::<1>::variable(0., 0),
            ))
            .keating_shadwick_omega_ratio;
            assert!((ad_bac.d[0] - fd_bac).abs() < 1e-5);
        }
    }
}
pub use evaluate_ohara_glosten_milgrom_sequential_spread_module::*;

mod evaluate_foucault_kyle_lambda_and_ekop_pin_module {
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
        let pin = (s.event_prob_alpha * mu)
            / (s.event_prob_alpha * mu + eb + es).max_val(T::from_f64(1e-12));
        let num_ev = s.event_prob_alpha
            * (s.bad_news_prob_delta * pb + (T::one() - s.bad_news_prob_delta) * pg);
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
        let pf = (T::one()
            + s.fast_decay_phi * s.transaction_cost_lambda / s.risk_aversion_gamma * a)
            .powf(-1.0);
        let ps = (T::one()
            + s.slow_decay_phi * s.transaction_cost_lambda / s.risk_aversion_gamma * a)
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
        let h0 =
            s.omega[0] / (T::one() - s.alpha[0] - s.gjr_gamma[0] / T::from_f64(2.0) - s.beta[0]);
        let h1 =
            s.omega[1] / (T::one() - s.alpha[1] - s.gjr_gamma[1] / T::from_f64(2.0) - s.beta[1]);
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
            h[0] = s.omega[0]
                + s.alpha[0] * r[0] * r[0]
                + s.gjr_gamma[0] * neg0 * neg0
                + s.beta[0] * h[0];
            h[1] = s.omega[1]
                + s.alpha[1] * r[1] * r[1]
                + s.gjr_gamma[1] * neg1 * neg1
                + s.beta[1] * h[1];
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
        let tau_m1 =
            (tau.max_val(T::from_f64(1e-12)).ln() * (s.critical_exponent_m - T::one())).exp();
        let ph = s.log_periodic_omega * tau.max_val(T::from_f64(1e-12)).ln() - s.phase_phi;
        let osc = s.log_periodic_amp_c * ph.cos();
        let lp = s.peak_log_price_a
            + s.power_law_slope_b * tau_m
            + s.log_periodic_amp_c * tau_m * ph.cos();
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
                let c_abs =
                    (s.log_periodic_amp_c * s.log_periodic_amp_c + T::from_f64(1e-30)).sqrt();
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
            spectral_weights_sum_error: ((w_sum - T::one()) * (w_sum - T::one())
                + T::from_f64(1e-30))
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
            if diff.abs() < T::from_f64(1e-10) {
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
            let v = if dx >= T::zero() {
                dx.max_val(T::from_f64(1e-12)).pow(s.tk_value_exponent)
            } else {
                -s.loss_aversion_lambda * (-dx).max_val(T::from_f64(1e-12)).pow(s.tk_value_exponent)
            };
            cpt += pi_i * v;
        }
        let inv_gamma = T::one() / s.tk_value_exponent.max_val(T::from_f64(1e-6));
        let ce = if cpt >= T::zero() {
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
    mod domain_suite_46_tests {
        use super::*;
        use ql_ad::Dual;
        fn eval_all<T: Real>(x: T) -> [T; 7] {
            let dx = x * T::from_f64(0.01);
            let kyle = evaluate_foucault_kyle_lambda_and_ekop_pin(&FoucaultKyleEkopPinSpec {
                prior_price_p0: T::from_f64(100.0),
                insider_valuation_v: T::from_f64(105.0) + x,
                fundamental_std_sigma_v: T::from_f64(2.0),
                noise_order_std_sigma_u: T::from_f64(10.0),
                event_prob_alpha: T::from_f64(0.3),
                bad_news_prob_delta: T::from_f64(0.4),
                informed_arrival_mu: T::from_f64(40.0),
                uninformed_buy_eps_b: T::from_f64(50.0),
                uninformed_sell_eps_s: T::from_f64(50.0),
                observed_buys_b: T::from_f64(75.0),
                observed_sells_s: T::from_f64(48.0),
            });
            let gp =
                evaluate_garleanu_pedersen_dynamic_aim_portfolio(&GarleanuPedersenDynamicAimSpec {
                    current_weights: [T::from_f64(0.2), T::from_f64(0.1)],
                    fast_alpha_signal: [T::from_f64(0.04) + dx, T::from_f64(0.02)],
                    slow_alpha_signal: [T::from_f64(0.03), T::from_f64(0.025)],
                    fast_decay_phi: T::from_f64(0.5),
                    slow_decay_phi: T::from_f64(0.1),
                    covariance_2x2: [
                        [T::from_f64(0.04), T::from_f64(0.01)],
                        [T::from_f64(0.01), T::from_f64(0.05)],
                    ],
                    risk_aversion_gamma: T::from_f64(2.0),
                    transaction_cost_lambda: T::from_f64(0.5),
                    discount_rate_rho: T::from_f64(0.05),
                });
            let gjr = evaluate_tsay_gjr_and_dcc_garch_dynamic_correlation(&TsayGjrDccGarchSpec {
                bivariate_returns: vec![
                    [T::from_f64(-0.02) + dx, T::from_f64(-0.015)],
                    [T::from_f64(0.012), T::from_f64(-0.01)],
                    [T::from_f64(-0.018), T::from_f64(-0.022)],
                ],
                omega: [T::from_f64(0.00002), T::from_f64(0.00003)],
                alpha: [T::from_f64(0.05), T::from_f64(0.06)],
                gjr_gamma: [T::from_f64(0.08), T::from_f64(0.07)],
                beta: [T::from_f64(0.85), T::from_f64(0.84)],
                unconditional_corr: T::from_f64(0.5),
                dcc_a: T::from_f64(0.05),
                dcc_b: T::from_f64(0.90),
            });
            let lppl = evaluate_sornette_jls_lppl_crash_hazard(&SornetteJlsLpplCrashHazardSpec {
                current_time_t: T::from_f64(1.0),
                critical_time_tc: T::from_f64(1.5) + dx,
                horizon_dt: T::from_f64(0.1),
                peak_log_price_a: T::from_f64(5.0),
                power_law_slope_b: T::from_f64(-0.4),
                log_periodic_amp_c: T::from_f64(0.03),
                critical_exponent_m: T::from_f64(0.5),
                log_periodic_omega: T::from_f64(6.5),
                phase_phi: T::from_f64(1.0),
                crash_loss_kappa: T::from_f64(0.2),
            });
            let spec = evaluate_follmer_schied_entropic_and_spectral_risk(
                &FollmerSchiedEntropicSpectralRiskSpec {
                    ordered_losses: vec![
                        T::from_f64(0.5),
                        T::from_f64(1.2),
                        T::from_f64(2.5),
                        T::from_f64(4.0) + dx,
                    ],
                    entropic_theta: T::from_f64(0.5),
                    spectral_power_gamma: T::from_f64(2.5),
                },
            );
            let brin =
                evaluate_bacon_brinson_fachler_and_carino_linking(&BaconBrinsonFachlerCarinoSpec {
                    portfolio_sector_weights: vec![T::from_f64(0.6), T::from_f64(0.4)],
                    benchmark_sector_weights: vec![T::from_f64(0.5), T::from_f64(0.5)],
                    portfolio_sector_returns: vec![T::from_f64(0.08) + dx, T::from_f64(0.03)],
                    benchmark_sector_returns: vec![T::from_f64(0.06), T::from_f64(0.02)],
                    multi_period_port_returns: vec![
                        T::from_f64(0.04) + dx,
                        T::from_f64(0.03),
                        T::from_f64(-0.01),
                    ],
                    multi_period_bench_returns: vec![
                        T::from_f64(0.03),
                        T::from_f64(0.02),
                        T::from_f64(-0.015),
                    ],
                });
            let prel = evaluate_wakker_prelec_rdu_and_choquet_ambiguity(
                &WakkerPrelecChoquetProspectSpec {
                    sorted_outcomes: [
                        T::from_f64(-10.0),
                        T::from_f64(-2.0),
                        T::from_f64(8.0),
                        T::from_f64(20.0) + x,
                    ],
                    state_probabilities: [
                        T::from_f64(0.1),
                        T::from_f64(0.3),
                        T::from_f64(0.4),
                        T::from_f64(0.2),
                    ],
                    reference_point: T::zero(),
                    prelec_curvature_alpha: T::from_f64(0.65),
                    prelec_elevation_beta: T::from_f64(1.0),
                    tk_value_exponent: T::from_f64(0.88),
                    loss_aversion_lambda: T::from_f64(2.25),
                    ambiguity_aversion_delta: T::from_f64(0.15),
                },
            );
            [
                kyle.expected_insider_profit,
                gp.expected_net_utility_gain,
                gjr.terminal_equal_weight_portfolio_vol,
                lppl.lppl_fitted_price,
                spec.entropic_convex_risk_measure + spec.kusuoka_power_spectral_risk_measure,
                brin.carino_smoothed_active_return,
                prel.prospect_certainty_equivalent,
            ]
        }
        #[test]
        fn suite46_portfolio_invariants_and_dual_fd() {
            let h = 1e-6;
            let up = eval_all(h);
            let dn = eval_all(-h);
            let ad = eval_all(Dual::<1>::variable(0.0, 0));
            for i in 0..7 {
                let fd = (up[i] - dn[i]) / (2.0 * h);
                assert!(ad[i].d[0].abs() > 1e-7 && (ad[i].d[0] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_foucault_kyle_lambda_and_ekop_pin_module::*;

mod evaluate_treynor_black_appraisal_portfolio_module {
    // Treynor & Black (1973, JB) & Grinold–Kahn (2000) Active-Passive Portfolio Separation,
    // Security Appraisal Ratio (`alpha_i / sigma_eps_i`) & Squared-Sharpe Theorem Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct TreynorBlackAppraisalPortfolioSpec<T: Real = f64> {
        pub security_alphas: [T; 3],
        pub security_market_betas: [T; 3],
        pub security_residual_variances: [T; 3],
        pub market_excess_return_mu_m: T,
        pub market_return_variance_sigma_m2: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct TreynorBlackAppraisalPortfolioReport<T: Real = f64> {
        pub active_portfolio_relative_weights: [T; 3],
        pub aggregate_active_alpha: T,
        pub aggregate_active_residual_variance: T,
        pub active_portfolio_appraisal_ratio: T,
        pub beta_adjusted_active_allocation_wa: T,
        pub combined_tangency_sharpe_ratio: T,
        pub sharpe_squared_orthogonality_residual: T,
    }

    pub fn evaluate_treynor_black_appraisal_portfolio<T: Real>(
        spec: &TreynorBlackAppraisalPortfolioSpec<T>,
    ) -> TreynorBlackAppraisalPortfolioReport<T> {
        let a0 = spec.security_alphas[0];
        let a1 = spec.security_alphas[1];
        let a2 = spec.security_alphas[2];

        let b0 = spec.security_market_betas[0];
        let b1 = spec.security_market_betas[1];
        let b2 = spec.security_market_betas[2];

        let v0 = spec.security_residual_variances[0].max_val(T::from_f64(1e-5));
        let v1 = spec.security_residual_variances[1].max_val(T::from_f64(1e-5));
        let v2 = spec.security_residual_variances[2].max_val(T::from_f64(1e-5));

        let mu_m = spec.market_excess_return_mu_m.max_val(T::from_f64(1e-4));
        let sig_m2 = spec
            .market_return_variance_sigma_m2
            .max_val(T::from_f64(1e-4));

        // Treynor-Black relative weights within the active portfolio: z_i = alpha_i / sigma_{eps, i}^2
        let z0 = a0 / v0;
        let z1 = a1 / v1;
        let z2 = a2 / v2;
        let z_sum = z0 + z1 + z2;
        let denom_z = (z_sum * z_sum + T::from_f64(1e-12)).sqrt();

        let w0 = z0 / denom_z;
        let w1 = z1 / denom_z;
        let w2 = z2 / denom_z;

        let alpha_a = w0 * a0 + w1 * a1 + w2 * a2;
        let beta_a = w0 * b0 + w1 * b1 + w2 * b2;
        let var_eps_a = (w0 * w0 * v0 + w1 * w1 * v1 + w2 * w2 * v2).max_val(T::from_f64(1e-6));

        // Aggregate Appraisal Ratio squared AR_A^2 = sum_i (alpha_i^2 / sigma_{eps, i}^2):
        let ar2_sum = (a0 * a0) / v0 + (a1 * a1) / v1 + (a2 * a2) / v2;
        let ar_a = ar2_sum.sqrt();

        // Unadjusted active weight w_A^0 and beta-adjusted weight w_A* (Treynor & Black 1973, Eqs. 16-18):
        let w_a0 = (alpha_a / var_eps_a) / (mu_m / sig_m2);
        let w_a_star = w_a0 / (T::one() + (T::one() - beta_a) * w_a0);

        let sr_m2 = (mu_m * mu_m) / sig_m2;
        let sr_p = (sr_m2 + ar2_sum).sqrt();
        let ar2_from_agg = (alpha_a * alpha_a) / var_eps_a;
        let ortho_res = (sr_p * sr_p) - (sr_m2 + ar2_from_agg);

        TreynorBlackAppraisalPortfolioReport {
            active_portfolio_relative_weights: [w0, w1, w2],
            aggregate_active_alpha: alpha_a,
            aggregate_active_residual_variance: var_eps_a,
            active_portfolio_appraisal_ratio: ar_a,
            beta_adjusted_active_allocation_wa: w_a_star,
            combined_tangency_sharpe_ratio: sr_p,
            sharpe_squared_orthogonality_residual: ortho_res,
        }
    }

    #[cfg(test)]
    mod evaluate_treynor_black_appraisal_portfolio_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn treynor_black_appraisal_portfolio_invariants_and_ad_parity() {
            let spec = TreynorBlackAppraisalPortfolioSpec {
                security_alphas: [0.024, 0.018, 0.012],
                security_market_betas: [1.15, 0.92, 1.05],
                security_residual_variances: [0.035, 0.025, 0.020],
                market_excess_return_mu_m: 0.065,
                market_return_variance_sigma_m2: 0.0256,
            };
            let rep = evaluate_treynor_black_appraisal_portfolio(&spec);
            assert!(rep.active_portfolio_appraisal_ratio > 0.0);
            assert!(rep.sharpe_squared_orthogonality_residual.abs() < 1e-8);

            let h = 1e-6;
            let r_ad =
                evaluate_treynor_black_appraisal_portfolio(&TreynorBlackAppraisalPortfolioSpec {
                    security_alphas: [
                        Dual::<11>::variable(spec.security_alphas[0], 0),
                        Dual::<11>::variable(spec.security_alphas[1], 1),
                        Dual::<11>::variable(spec.security_alphas[2], 2),
                    ],
                    security_market_betas: [
                        Dual::<11>::variable(spec.security_market_betas[0], 3),
                        Dual::<11>::variable(spec.security_market_betas[1], 4),
                        Dual::<11>::variable(spec.security_market_betas[2], 5),
                    ],
                    security_residual_variances: [
                        Dual::<11>::variable(spec.security_residual_variances[0], 6),
                        Dual::<11>::variable(spec.security_residual_variances[1], 7),
                        Dual::<11>::variable(spec.security_residual_variances[2], 8),
                    ],
                    market_excess_return_mu_m: Dual::<11>::variable(
                        spec.market_excess_return_mu_m,
                        9,
                    ),
                    market_return_variance_sigma_m2: Dual::<11>::variable(
                        spec.market_return_variance_sigma_m2,
                        10,
                    ),
                });
            for (idx, &ad) in r_ad.beta_adjusted_active_allocation_wa.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0..=2 => {
                        sp.security_alphas[idx] += h;
                        sm.security_alphas[idx] -= h;
                    }
                    3..=5 => {
                        sp.security_market_betas[idx - 3] += h;
                        sm.security_market_betas[idx - 3] -= h;
                    }
                    6..=8 => {
                        sp.security_residual_variances[idx - 6] += h;
                        sm.security_residual_variances[idx - 6] -= h;
                    }
                    9 => {
                        sp.market_excess_return_mu_m += h;
                        sm.market_excess_return_mu_m -= h;
                    }
                    _ => {
                        sp.market_return_variance_sigma_m2 += h;
                        sm.market_return_variance_sigma_m2 -= h;
                    }
                }
                let fd = (evaluate_treynor_black_appraisal_portfolio(&sp)
                    .beta_adjusted_active_allocation_wa
                    - evaluate_treynor_black_appraisal_portfolio(&sm)
                        .beta_adjusted_active_allocation_wa)
                    / (2.0 * h);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_treynor_black_appraisal_portfolio_module::*;

mod evaluate_tasche_litterman_hot_spots_best_hedge_es_decomposition_module {
    // Litterman (1996, Risk, "Hot Spots and Hedges") & Tasche (2002, J. Banking & Finance)
    // Component Expected Shortfall (Hot Spots), Best-Hedge Trade Sizing & Post-Hedge Risk Reduction Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct TascheLittermanHotSpotsBestHedgeSpec<T: Real = f64> {
        pub focal_position_weight_w1: T,
        pub core_portfolio_weight_w2: T,
        pub focal_asset_volatility_sigma1: T,
        pub core_portfolio_volatility_sigma2: T,
        pub asset_core_correlation_rho12: T,
        pub focal_expected_return_mu1: T,
        pub core_expected_return_mu2: T,
        pub gaussian_es_tail_multiplier_m_alpha: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct TascheLittermanHotSpotsBestHedgeReport<T: Real = f64> {
        pub total_portfolio_expected_shortfall: T,
        pub focal_marginal_expected_shortfall: T,
        pub focal_component_expected_shortfall: T,
        pub focal_hot_spot_percentage_share: T,
        pub litterman_best_hedge_trade_delta_w1: T,
        pub post_best_hedge_expected_shortfall: T,
        pub best_hedge_es_reduction: T,
    }

    pub fn evaluate_tasche_litterman_hot_spots_best_hedge_es_decomposition<T: Real>(
        spec: &TascheLittermanHotSpotsBestHedgeSpec<T>,
    ) -> TascheLittermanHotSpotsBestHedgeReport<T> {
        let w1 = spec.focal_position_weight_w1;
        let w2 = spec.core_portfolio_weight_w2;
        let sig1 = spec
            .focal_asset_volatility_sigma1
            .max_val(T::from_f64(0.01));
        let sig2 = spec
            .core_portfolio_volatility_sigma2
            .max_val(T::from_f64(0.01));
        let rho = spec
            .asset_core_correlation_rho12
            .max_val(T::from_f64(-0.95))
            .min_val(T::from_f64(0.95));
        let mu1 = spec.focal_expected_return_mu1;
        let mu2 = spec.core_expected_return_mu2;
        let m_alpha = spec
            .gaussian_es_tail_multiplier_m_alpha
            .max_val(T::from_f64(1.5));

        let cov1p = w1 * sig1 * sig1 + w2 * rho * sig1 * sig2;
        let var_p = (w1 * w1 * sig1 * sig1
            + T::from_f64(2.0) * w1 * w2 * rho * sig1 * sig2
            + w2 * w2 * sig2 * sig2)
            .max_val(T::from_f64(1e-6));
        let sig_p = var_p.sqrt();

        let es_total = m_alpha * sig_p - (w1 * mu1 + w2 * mu2);
        let mes1 = m_alpha * (cov1p / sig_p) - mu1;
        let ces1 = w1 * mes1;
        let hot_spot_pct = (ces1 / es_total.max_val(T::from_f64(1e-4))) * T::from_f64(100.0);

        // Litterman (1996) Best-Hedge trade size in asset 1 that orthogonalizes portfolio against asset 1:
        let delta_w1_best = -cov1p / (sig1 * sig1);
        let w1_hedged = w1 + delta_w1_best;
        let sig_p_hedged =
            (w2.abs() * sig2 * (T::one() - rho * rho).sqrt()).max_val(T::from_f64(1e-4));
        let es_hedged = m_alpha * sig_p_hedged - (w1_hedged * mu1 + w2 * mu2);
        let es_reduction = es_total - es_hedged;

        TascheLittermanHotSpotsBestHedgeReport {
            total_portfolio_expected_shortfall: es_total,
            focal_marginal_expected_shortfall: mes1,
            focal_component_expected_shortfall: ces1,
            focal_hot_spot_percentage_share: hot_spot_pct,
            litterman_best_hedge_trade_delta_w1: delta_w1_best,
            post_best_hedge_expected_shortfall: es_hedged,
            best_hedge_es_reduction: es_reduction,
        }
    }

    #[cfg(test)]
    mod evaluate_tasche_litterman_hot_spots_best_hedge_es_decomposition_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn tasche_litterman_hot_spots_best_hedge_invariants_and_ad_parity() {
            let spec = TascheLittermanHotSpotsBestHedgeSpec {
                focal_position_weight_w1: 0.45,
                core_portfolio_weight_w2: 0.55,
                focal_asset_volatility_sigma1: 0.28,
                core_portfolio_volatility_sigma2: 0.16,
                asset_core_correlation_rho12: 0.50,
                focal_expected_return_mu1: 0.08,
                core_expected_return_mu2: 0.06,
                gaussian_es_tail_multiplier_m_alpha: 2.338,
            };
            let rep = evaluate_tasche_litterman_hot_spots_best_hedge_es_decomposition(&spec);
            assert!(rep.focal_hot_spot_percentage_share > 50.0);
            assert!(rep.litterman_best_hedge_trade_delta_w1 < -0.45);
            assert!(rep.best_hedge_es_reduction > 0.0);

            let h = 1e-5;
            let r_ad = evaluate_tasche_litterman_hot_spots_best_hedge_es_decomposition(
                &TascheLittermanHotSpotsBestHedgeSpec {
                    focal_position_weight_w1: Dual::<8>::variable(spec.focal_position_weight_w1, 0),
                    core_portfolio_weight_w2: Dual::<8>::variable(spec.core_portfolio_weight_w2, 1),
                    focal_asset_volatility_sigma1: Dual::<8>::variable(
                        spec.focal_asset_volatility_sigma1,
                        2,
                    ),
                    core_portfolio_volatility_sigma2: Dual::<8>::variable(
                        spec.core_portfolio_volatility_sigma2,
                        3,
                    ),
                    asset_core_correlation_rho12: Dual::<8>::variable(
                        spec.asset_core_correlation_rho12,
                        4,
                    ),
                    focal_expected_return_mu1: Dual::<8>::variable(
                        spec.focal_expected_return_mu1,
                        5,
                    ),
                    core_expected_return_mu2: Dual::<8>::variable(spec.core_expected_return_mu2, 6),
                    gaussian_es_tail_multiplier_m_alpha: Dual::<8>::variable(
                        spec.gaussian_es_tail_multiplier_m_alpha,
                        7,
                    ),
                },
            );
            let comp_ad =
                r_ad.total_portfolio_expected_shortfall + r_ad.focal_component_expected_shortfall;
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.focal_position_weight_w1 += h;
                        sm.focal_position_weight_w1 -= h;
                    }
                    1 => {
                        sp.core_portfolio_weight_w2 += h;
                        sm.core_portfolio_weight_w2 -= h;
                    }
                    2 => {
                        sp.focal_asset_volatility_sigma1 += h;
                        sm.focal_asset_volatility_sigma1 -= h;
                    }
                    3 => {
                        sp.core_portfolio_volatility_sigma2 += h;
                        sm.core_portfolio_volatility_sigma2 -= h;
                    }
                    4 => {
                        sp.asset_core_correlation_rho12 += h;
                        sm.asset_core_correlation_rho12 -= h;
                    }
                    5 => {
                        sp.focal_expected_return_mu1 += h;
                        sm.focal_expected_return_mu1 -= h;
                    }
                    6 => {
                        sp.core_expected_return_mu2 += h;
                        sm.core_expected_return_mu2 -= h;
                    }
                    _ => {
                        sp.gaussian_es_tail_multiplier_m_alpha += h;
                        sm.gaussian_es_tail_multiplier_m_alpha -= h;
                    }
                }
                let rp = evaluate_tasche_litterman_hot_spots_best_hedge_es_decomposition(&sp);
                let rm = evaluate_tasche_litterman_hot_spots_best_hedge_es_decomposition(&sm);
                let fd = ((rp.total_portfolio_expected_shortfall
                    + rp.focal_component_expected_shortfall)
                    - (rm.total_portfolio_expected_shortfall
                        + rm.focal_component_expected_shortfall))
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_tasche_litterman_hot_spots_best_hedge_es_decomposition_module::*;

mod evaluate_acerbi_tasche_exponential_spectral_risk_measure_module {
    // Acerbi (2002, Journal of Banking & Finance 26(7), "Spectral Measures of Risk:
    // A Coherent Representation of Subjective Risk Aversion"), Acerbi & Tasche (2002, JBF 26(7)),
    // & McNeil, Frey & Embrechts (2015, Quantitative Risk Management, Ch. 8)
    // Coherent Exponential Spectral Risk Measure (SRM) & Cornish-Fisher Tail Distortion Engine.

    use ql_ad::Real;
    use ql_math::{inv_norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct AcerbiTascheExponentialSpectralRiskSpec<T: Real = f64> {
        pub portfolio_expected_return_mu: T,
        pub portfolio_volatility_sigma: T,
        pub return_skewness_gamma1: T,
        pub return_excess_kurtosis_gamma2: T,
        pub cara_spectral_risk_aversion_k: T,
        pub reference_es_confidence_alpha: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct AcerbiTascheExponentialSpectralRiskReport<T: Real = f64> {
        pub gaussian_expected_shortfall_es: T,
        pub cornish_fisher_tail_var_alpha: T,
        pub exponential_spectral_risk_measure: T,
        pub spectral_tail_aversion_add_on: T,
        pub terminal_spectral_weight_at_one: T,
    }

    fn cornish_fisher_loss_quantile<T: Real>(z: T, gamma1: T, gamma2: T) -> T {
        let z2 = z * z;
        let z3 = z2 * z;
        // For loss L = -R, skewness of L is -gamma1:
        z - (gamma1 / T::from_f64(6.0)) * (z2 - T::one())
            + (gamma2 / T::from_f64(24.0)) * (z3 - T::from_f64(3.0) * z)
            - (gamma1 * gamma1 / T::from_f64(36.0)) * (T::from_f64(2.0) * z3 - T::from_f64(5.0) * z)
    }

    pub fn evaluate_acerbi_tasche_exponential_spectral_risk_measure<T: Real>(
        spec: &AcerbiTascheExponentialSpectralRiskSpec<T>,
    ) -> AcerbiTascheExponentialSpectralRiskReport<T> {
        let mu = spec.portfolio_expected_return_mu;
        let sig = spec.portfolio_volatility_sigma.max_val(T::from_f64(0.005));
        let g1 = spec.return_skewness_gamma1;
        let g2 = spec.return_excess_kurtosis_gamma2.max_val(T::from_f64(0.0));
        let k = spec.cara_spectral_risk_aversion_k.max_val(T::from_f64(1.0));
        let alpha = spec
            .reference_es_confidence_alpha
            .max_val(T::from_f64(0.90))
            .min_val(T::from_f64(0.999));

        let z_a = inv_norm_cdf(alpha);
        let es_gauss = -mu + sig * (norm_pdf(z_a) / (T::one() - alpha));
        let var_cf = -mu + sig * cornish_fisher_loss_quantile(z_a, g1, g2);

        // Exact importance-sampling inverse distortion u in (0, 1) -> p(u) = 1 + ln(e^{-k} + u*(1 - e^{-k})) / k
        // under Acerbi (2002) exponential spectrum phi_k(p) = k * e^{-k*(1 - p)} / (1 - e^{-k}):
        const GL_NODES: [f64; 4] = [
            -0.861_136_311_594_052_6,
            -0.339_981_043_584_856_26,
            0.339_981_043_584_856_26,
            0.861_136_311_594_052_6,
        ];
        const GL_WEIGHTS: [f64; 4] = [
            0.347_854_845_137_453_85,
            0.652_145_154_862_546_1,
            0.652_145_154_862_546_1,
            0.347_854_845_137_453_85,
        ];
        let ek = (-k).exp();
        let one_minus_ek = (T::one() - ek).max_val(T::from_f64(1e-8));
        let mut q_spec = T::zero();
        for (&xi, &wi) in GL_NODES.iter().zip(GL_WEIGHTS.iter()) {
            let u = T::from_f64(0.5 * (1.0 + xi));
            let p_u = (T::one() + (ek + u * one_minus_ek).ln() / k)
                .max_val(T::from_f64(1e-4))
                .min_val(T::from_f64(0.9995));
            let z_u = inv_norm_cdf(p_u);
            q_spec += T::from_f64(0.5 * wi) * cornish_fisher_loss_quantile(z_u, g1, g2);
        }

        let srm = -mu + sig * q_spec;
        let add_on = srm - es_gauss;
        let phi_1 = k / one_minus_ek;

        AcerbiTascheExponentialSpectralRiskReport {
            gaussian_expected_shortfall_es: es_gauss,
            cornish_fisher_tail_var_alpha: var_cf,
            exponential_spectral_risk_measure: srm,
            spectral_tail_aversion_add_on: add_on,
            terminal_spectral_weight_at_one: phi_1,
        }
    }

    #[cfg(test)]
    mod evaluate_acerbi_tasche_exponential_spectral_risk_measure_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn acerbi_tasche_exponential_spectral_risk_invariants_and_ad_parity() {
            let spec = AcerbiTascheExponentialSpectralRiskSpec {
                portfolio_expected_return_mu: 0.08,
                portfolio_volatility_sigma: 0.18,
                return_skewness_gamma1: -0.55,
                return_excess_kurtosis_gamma2: 2.40,
                cara_spectral_risk_aversion_k: 25.0,
                reference_es_confidence_alpha: 0.975,
            };
            let rep = evaluate_acerbi_tasche_exponential_spectral_risk_measure(&spec);
            assert!(rep.exponential_spectral_risk_measure > 0.20);
            assert!(rep.terminal_spectral_weight_at_one > 24.9);

            let h = 1e-6;
            let r_ad = evaluate_acerbi_tasche_exponential_spectral_risk_measure(
                &AcerbiTascheExponentialSpectralRiskSpec {
                    portfolio_expected_return_mu: Dual::<6>::variable(
                        spec.portfolio_expected_return_mu,
                        0,
                    ),
                    portfolio_volatility_sigma: Dual::<6>::variable(
                        spec.portfolio_volatility_sigma,
                        1,
                    ),
                    return_skewness_gamma1: Dual::<6>::variable(spec.return_skewness_gamma1, 2),
                    return_excess_kurtosis_gamma2: Dual::<6>::variable(
                        spec.return_excess_kurtosis_gamma2,
                        3,
                    ),
                    cara_spectral_risk_aversion_k: Dual::<6>::variable(
                        spec.cara_spectral_risk_aversion_k,
                        4,
                    ),
                    reference_es_confidence_alpha: Dual::<6>::variable(
                        spec.reference_es_confidence_alpha,
                        5,
                    ),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.portfolio_expected_return_mu,
                    spec.portfolio_volatility_sigma,
                    spec.return_skewness_gamma1,
                    spec.return_excess_kurtosis_gamma2,
                    spec.cara_spectral_risk_aversion_k,
                    spec.reference_es_confidence_alpha,
                ];
                v[i] += h;
                let up = evaluate_acerbi_tasche_exponential_spectral_risk_measure(
                    &AcerbiTascheExponentialSpectralRiskSpec {
                        portfolio_expected_return_mu: v[0],
                        portfolio_volatility_sigma: v[1],
                        return_skewness_gamma1: v[2],
                        return_excess_kurtosis_gamma2: v[3],
                        cara_spectral_risk_aversion_k: v[4],
                        reference_es_confidence_alpha: v[5],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_acerbi_tasche_exponential_spectral_risk_measure(
                    &AcerbiTascheExponentialSpectralRiskSpec {
                        portfolio_expected_return_mu: v[0],
                        portfolio_volatility_sigma: v[1],
                        return_skewness_gamma1: v[2],
                        return_excess_kurtosis_gamma2: v[3],
                        cara_spectral_risk_aversion_k: v[4],
                        reference_es_confidence_alpha: v[5],
                    },
                );
                let val_ad =
                    r_ad.exponential_spectral_risk_measure + r_ad.gaussian_expected_shortfall_es;
                let fd = ((up.exponential_spectral_risk_measure
                    + up.gaussian_expected_shortfall_es)
                    - (dn.exponential_spectral_risk_measure + dn.gaussian_expected_shortfall_es))
                    / (2.0 * h);
                assert!((val_ad.d[i] - fd).abs() < 1e-4, "input {i}");
            }
        }
    }
}
pub use evaluate_acerbi_tasche_exponential_spectral_risk_measure_module::*;

mod evaluate_bellini_cont_expectile_evar_influence_sensitivity_module {
    // Bellini, Klar, Müller & Rosazza Gianin (2014, Finance and Stochastics 18(4),
    // "Generalized Quantiles as Risk Measures"), Newey & Powell (1987, Econometrica 55(4)),
    // & Cont, Deguest & Scandolo (2010, Mathematical Finance 20(4))
    // Elicitable Coherent Expectile Risk Measure (EVaR) & Hampel Influence Function Engine.

    use ql_ad::Real;
    use ql_math::{inv_norm_cdf, norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct BelliniContExpectileRiskMeasureSpec<T: Real = f64> {
        pub loss_mean_mu: T,
        pub loss_volatility_sigma: T,
        pub tail_asymmetry_tau: T,
        pub contamination_loss_shock_z: T,
        pub contamination_weight_eps: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BelliniContExpectileRiskMeasureReport<T: Real = f64> {
        pub standardized_gaussian_expectile: T,
        pub coherent_expectile_risk_measure: T,
        pub gaussian_value_at_risk_benchmark: T,
        pub hampel_expectile_influence_function: T,
        pub contaminated_robust_expectile: T,
        pub expectile_to_var_tail_ratio: T,
    }

    pub fn evaluate_bellini_cont_expectile_evar_influence_sensitivity<T: Real>(
        spec: &BelliniContExpectileRiskMeasureSpec<T>,
    ) -> BelliniContExpectileRiskMeasureReport<T> {
        let mu = spec.loss_mean_mu;
        let sig = spec.loss_volatility_sigma.max_val(T::from_f64(1e-4));
        let tau = spec
            .tail_asymmetry_tau
            .max_val(T::from_f64(0.501))
            .min_val(T::from_f64(0.999));
        let z_shock = spec.contamination_loss_shock_z;
        let eps = spec
            .contamination_weight_eps
            .max_val(T::from_f64(1e-4))
            .min_val(T::from_f64(0.20));

        // Solve Newey & Powell (1987) / Bellini et al. (2014) standardized Gaussian expectile equation:
        // g(e) = (1 - 2*tau)*phi(e) + ((1 - tau)*Phi(e) + tau*(1 - Phi(e)))*e = 0 via 4 Newton steps:
        let z_q = inv_norm_cdf(tau);
        let mut e_std = z_q;
        let one_minus_tau = T::one() - tau;
        let one_minus_two_tau = T::one() - T::from_f64(2.0) * tau;

        for _ in 0..4 {
            let phi_e = norm_pdf(e_std);
            let cdf_e = norm_cdf(e_std);
            let gp = one_minus_tau * cdf_e + tau * (T::one() - cdf_e);
            let g = one_minus_two_tau * phi_e + gp * e_std;
            e_std -= g / gp;
        }

        let evar = mu + sig * e_std;
        let var_q = mu + sig * z_q;

        // Cont, Deguest & Scandolo (2010): Hampel influence function of expectile at tail shock z_shock > evar:
        let cdf_final = norm_cdf(e_std);
        let denom_if = one_minus_tau * cdf_final + tau * (T::one() - cdf_final);
        let diff_z = z_shock - evar;
        let pos_diff = T::from_f64(0.5) * (diff_z + (diff_z * diff_z + T::from_f64(1e-4)).sqrt());
        let hampel_if = (tau * pos_diff) / denom_if;
        let evar_contam = evar + eps * hampel_if;
        let ratio = e_std / z_q.max_val(T::from_f64(1e-4));

        BelliniContExpectileRiskMeasureReport {
            standardized_gaussian_expectile: e_std,
            coherent_expectile_risk_measure: evar,
            gaussian_value_at_risk_benchmark: var_q,
            hampel_expectile_influence_function: hampel_if,
            contaminated_robust_expectile: evar_contam,
            expectile_to_var_tail_ratio: ratio,
        }
    }

    #[cfg(test)]
    mod evaluate_bellini_cont_expectile_evar_influence_sensitivity_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bellini_cont_expectile_risk_measure_invariants_and_ad_parity() {
            let spec = BelliniContExpectileRiskMeasureSpec {
                loss_mean_mu: 0.50,
                loss_volatility_sigma: 2.20,
                tail_asymmetry_tau: 0.975,
                contamination_loss_shock_z: 8.50,
                contamination_weight_eps: 0.02,
            };
            let rep = evaluate_bellini_cont_expectile_evar_influence_sensitivity(&spec);
            assert!(rep.standardized_gaussian_expectile > 1.0);
            assert!(rep.contaminated_robust_expectile > rep.coherent_expectile_risk_measure);

            let h = 1e-6;
            let r_ad = evaluate_bellini_cont_expectile_evar_influence_sensitivity(
                &BelliniContExpectileRiskMeasureSpec {
                    loss_mean_mu: Dual::<5>::variable(spec.loss_mean_mu, 0),
                    loss_volatility_sigma: Dual::<5>::variable(spec.loss_volatility_sigma, 1),
                    tail_asymmetry_tau: Dual::<5>::variable(spec.tail_asymmetry_tau, 2),
                    contamination_loss_shock_z: Dual::<5>::variable(
                        spec.contamination_loss_shock_z,
                        3,
                    ),
                    contamination_weight_eps: Dual::<5>::variable(spec.contamination_weight_eps, 4),
                },
            );
            for i in 0..5 {
                let mut v = [
                    spec.loss_mean_mu,
                    spec.loss_volatility_sigma,
                    spec.tail_asymmetry_tau,
                    spec.contamination_loss_shock_z,
                    spec.contamination_weight_eps,
                ];
                v[i] += h;
                let up = evaluate_bellini_cont_expectile_evar_influence_sensitivity(
                    &BelliniContExpectileRiskMeasureSpec {
                        loss_mean_mu: v[0],
                        loss_volatility_sigma: v[1],
                        tail_asymmetry_tau: v[2],
                        contamination_loss_shock_z: v[3],
                        contamination_weight_eps: v[4],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_bellini_cont_expectile_evar_influence_sensitivity(
                    &BelliniContExpectileRiskMeasureSpec {
                        loss_mean_mu: v[0],
                        loss_volatility_sigma: v[1],
                        tail_asymmetry_tau: v[2],
                        contamination_loss_shock_z: v[3],
                        contamination_weight_eps: v[4],
                    },
                );
                let fd = (up.contaminated_robust_expectile - dn.contaminated_robust_expectile)
                    / (2.0 * h);
                assert!(
                    (r_ad.contaminated_robust_expectile.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_bellini_cont_expectile_evar_influence_sensitivity_module::*;

/// Specification for the Ahmadi-Javid (2012) / Föllmer–Schied (2002) Entropic Value-at-Risk ($\text{EVaR}_{1-\alpha}$)
/// coherent risk measure under compound Poisson-Gaussian jump losses.
///
/// Reference: Ahmadi-Javid (2012, *Journal of Optimization Theory and Applications* 155(3),
/// "Entropic Value-at-Risk: A New Coherent Risk Measure"), Föllmer & Schied (2002, *Finance & Stochastics* 6(4)).
#[derive(Debug, Clone, Copy)]
pub struct AhmadiJavidEntropicVarSpec<T: Real> {
    /// Expected diffusive loss $\mu_L$.
    pub mean_diffusive_loss_mu: T,
    /// Diffusive loss standard deviation $\sigma_L$.
    pub diffusive_loss_volatility_sigma: T,
    /// Poisson tail jump intensity $\lambda_J$.
    pub jump_loss_intensity_lambda: T,
    /// Mean jump loss severity $\mu_J$.
    pub mean_jump_loss_mu_j: T,
    /// Jump loss severity standard deviation $\delta_J$.
    pub jump_loss_volatility_delta_j: T,
    /// Tail significance level $\alpha \in (0, 1)$ (e.g. $\alpha = 0.01$ for $99\%$ EVaR).
    pub tail_significance_alpha: T,
}

/// Analytics returned by [`evaluate_ahmadi_javid_follmer_schied_entropic_var_evar`].
#[derive(Debug, Clone, Copy)]
pub struct AhmadiJavidEntropicVarAnalytics<T: Real> {
    /// Pure Gaussian Entropic Value-at-Risk benchmark $\mu_L + \sigma_L \sqrt{-2 \ln \alpha}$.
    pub gaussian_benchmark_evar: T,
    /// Optimal Chernoff / Esscher risk-aversion saddlepoint parameter $z^* > 0$.
    pub optimal_chernoff_saddlepoint_z: T,
    /// Exact Jump-Diffusion Entropic Value-at-Risk $\text{EVaR}_{1-\alpha}(L) = \inf_{z>0} \frac{\psi_L(z) - \ln \alpha}{z}$.
    pub jump_diffusion_entropic_var: T,
    /// Föllmer–Schied worst-case Esscher-tilted jump arrival intensity $\lambda_J^* = \lambda_J e^{\mu_J z^* + \frac{1}{2}\delta_J^2 (z^*)^2}$.
    pub worst_case_tilted_jump_intensity: T,
}

/// Evaluates the Ahmadi-Javid (2012) / Föllmer–Schied (2002) Entropic Value-at-Risk ($\text{EVaR}_{1-\alpha}$)
/// via 5 Newton–Raphson steps on the cumulant-generating saddlepoint equation $z \psi_L'(z) - \psi_L(z) + \ln \alpha = 0$.
pub fn evaluate_ahmadi_javid_follmer_schied_entropic_var_evar<T: Real>(
    spec: &AhmadiJavidEntropicVarSpec<T>,
) -> AhmadiJavidEntropicVarAnalytics<T> {
    let alpha = spec
        .tail_significance_alpha
        .max_val(T::from_f64(1e-5))
        .min_val(T::from_f64(0.5));
    let neg_ln_alpha = -alpha.ln();
    let sig = spec
        .diffusive_loss_volatility_sigma
        .max_val(T::from_f64(1e-4));
    let sig2 = sig * sig;
    let lam = spec.jump_loss_intensity_lambda;
    let mu_j = spec.mean_jump_loss_mu_j;
    let del_j = spec.jump_loss_volatility_delta_j;
    let del_j2 = del_j * del_j;

    let gauss_evar = spec.mean_diffusive_loss_mu + sig * (T::from_f64(2.0) * neg_ln_alpha).sqrt();

    // Initial guess from total second cumulant var_total = sig^2 + lam * (mu_j^2 + del_j^2)
    let var_total = sig2 + lam * (mu_j * mu_j + del_j2);
    let mut z = ((T::from_f64(2.0) * neg_ln_alpha) / var_total.max_val(T::from_f64(1e-6))).sqrt();

    // 5 exact Newton steps on g(z) = z * psi'(z) - psi(z) - (-ln alpha) = 0, where g'(z) = z * psi''(z)
    for _ in 0..5 {
        let j_exp_arg = (mu_j * z + T::from_f64(0.5) * del_j2 * z * z).min_val(T::from_f64(8.0));
        let j_mgf = j_exp_arg.exp();
        let j_d1 = mu_j + del_j2 * z;
        let psi = spec.mean_diffusive_loss_mu * z
            + T::from_f64(0.5) * sig2 * z * z
            + lam * (j_mgf - T::one());
        let psi_p = spec.mean_diffusive_loss_mu + sig2 * z + lam * j_mgf * j_d1;
        let psi_pp = sig2 + lam * j_mgf * (j_d1 * j_d1 + del_j2);

        let g = z * psi_p - psi - neg_ln_alpha;
        let gp = (z * psi_pp).max_val(T::from_f64(1e-6));
        z = (z - g / gp).max_val(T::from_f64(1e-4));
    }

    let j_exp_arg = (mu_j * z + T::from_f64(0.5) * del_j2 * z * z).min_val(T::from_f64(8.0));
    let j_mgf = j_exp_arg.exp();
    let psi_star = spec.mean_diffusive_loss_mu * z
        + T::from_f64(0.5) * sig2 * z * z
        + lam * (j_mgf - T::one());
    let evar = (psi_star + neg_ln_alpha) / z;
    let tilted_lambda = lam * j_mgf;

    AhmadiJavidEntropicVarAnalytics {
        gaussian_benchmark_evar: gauss_evar,
        optimal_chernoff_saddlepoint_z: z,
        jump_diffusion_entropic_var: evar,
        worst_case_tilted_jump_intensity: tilted_lambda,
    }
}

#[cfg(test)]
mod ahmadi_javid_evar_suite64_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn ahmadi_javid_entropic_var_invariants_and_ad_parity() {
        let base = AhmadiJavidEntropicVarSpec {
            mean_diffusive_loss_mu: 0.25_f64,
            diffusive_loss_volatility_sigma: 1.50,
            jump_loss_intensity_lambda: 0.15,
            mean_jump_loss_mu_j: 2.20,
            jump_loss_volatility_delta_j: 0.85,
            tail_significance_alpha: 0.025,
        };
        let r0 = evaluate_ahmadi_javid_follmer_schied_entropic_var_evar(&base);
        assert!(r0.gaussian_benchmark_evar > base.mean_diffusive_loss_mu);
        assert!(r0.jump_diffusion_entropic_var > r0.gaussian_benchmark_evar);
        assert!(r0.optimal_chernoff_saddlepoint_z > 0.0);
        assert!(r0.worst_case_tilted_jump_intensity > base.jump_loss_intensity_lambda);

        let inputs = [
            base.mean_diffusive_loss_mu,
            base.diffusive_loss_volatility_sigma,
            base.jump_loss_intensity_lambda,
            base.mean_jump_loss_mu_j,
            base.jump_loss_volatility_delta_j,
            base.tail_significance_alpha,
        ];
        let mut dual_in = [Dual::<6>::from_f64(0.0); 6];
        for i in 0..6 {
            dual_in[i] = Dual::variable(inputs[i], i);
        }
        let spec_ad = AhmadiJavidEntropicVarSpec {
            mean_diffusive_loss_mu: dual_in[0],
            diffusive_loss_volatility_sigma: dual_in[1],
            jump_loss_intensity_lambda: dual_in[2],
            mean_jump_loss_mu_j: dual_in[3],
            jump_loss_volatility_delta_j: dual_in[4],
            tail_significance_alpha: dual_in[5],
        };
        let r_ad = evaluate_ahmadi_javid_follmer_schied_entropic_var_evar(&spec_ad);
        let eps = 1e-6;
        for i in 0..6 {
            let mut p = inputs;
            let mut m = inputs;
            p[i] += eps;
            m[i] -= eps;
            let mk = |v: [f64; 6]| AhmadiJavidEntropicVarSpec {
                mean_diffusive_loss_mu: v[0],
                diffusive_loss_volatility_sigma: v[1],
                jump_loss_intensity_lambda: v[2],
                mean_jump_loss_mu_j: v[3],
                jump_loss_volatility_delta_j: v[4],
                tail_significance_alpha: v[5],
            };
            let rp = evaluate_ahmadi_javid_follmer_schied_entropic_var_evar(&mk(p));
            let rm = evaluate_ahmadi_javid_follmer_schied_entropic_var_evar(&mk(m));
            let fd =
                (rp.jump_diffusion_entropic_var - rm.jump_diffusion_entropic_var) / (2.0 * eps);
            assert!(
                r_ad.jump_diffusion_entropic_var.d[i].abs() > 1e-6,
                "zero derivative at {i}"
            );
            assert!((r_ad.jump_diffusion_entropic_var.d[i] - fd).abs() < 1e-4);
        }
    }
}

pub mod rockafellar_uryasev_mixed_cvar_quadrangle_mod {
    use ql_ad::Real;
    use ql_math::norm_pdf;

    /// Specification for the Rockafellar & Uryasev (2013, SORMS) / Föllmer & Schied (2016, Ch. 4)
    /// Fundamental Risk Quadrangle Mixed-CVaR spectral risk, deviation, and regret engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct RockafellarUryasevMixedCvarSpec<T: Real> {
        /// Expected portfolio loss mean $\mu_L$.
        pub mean_loss_mu: T,
        /// Portfolio loss scale / standard deviation $\sigma_L > 0$.
        pub loss_volatility_sigma: T,
        /// Moderate tail confidence level $\alpha_1 \in (0.5, 0.99)$ (e.g. 0.95).
        pub moderate_confidence_alpha1: T,
        /// Extreme tail confidence level $\alpha_2 \in (\alpha_1, 0.999)$ (e.g. 0.99).
        pub extreme_confidence_alpha2: T,
        /// Spectral mixture weight $\lambda_1 \in (0, 1)$ on $\text{CVaR}_{\alpha_1}$ (with $\lambda_2 = 1 - \lambda_1$).
        pub moderate_spectral_weight_lambda1: T,
        /// Excess tail kurtosis amplification $\kappa_{\text{tail}} \ge 0$.
        pub tail_kurtosis_amplification_kappa: T,
    }

    /// Analytical outputs for the Rockafellar-Uryasev Mixed-CVaR Risk Quadrangle engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct RockafellarUryasevMixedCvarAnalytics<T: Real> {
        /// Moderate tail Conditional Value-at-Risk $\text{CVaR}_{\alpha_1}(L)$.
        pub moderate_tail_cvar_alpha1: T,
        /// Extreme tail Conditional Value-at-Risk $\text{CVaR}_{\alpha_2}(L)$.
        pub extreme_tail_cvar_alpha2: T,
        /// Quadrangle Mixed-CVaR Deviation $\mathcal{D}(L) = \mathcal{R}(L - \mathbb{E}[L])$.
        pub quadrangle_mixed_deviation: T,
        /// Coherent Mixed-CVaR Spectral Risk $\mathcal{R}(L) = \lambda_1 \text{CVaR}_{\alpha_1} + (1-\lambda_1)\text{CVaR}_{\alpha_2}$.
        pub quadrangle_mixed_cvar_risk: T,
    }

    #[allow(clippy::excessive_precision)]
    fn inv_norm_cdf_tail<T: Real>(p: T) -> T {
        let p_clamped = p
            .max_val(T::from_f64(1e-9))
            .min_val(T::from_f64(1.0 - 1e-9));
        let q = p_clamped - T::from_f64(0.5);
        if q.abs() <= T::from_f64(0.425) {
            let r = T::from_f64(0.180625) - q * q;
            let num = ((((((T::from_f64(2.5090809287301227e3) * r
                + T::from_f64(3.3430575583588128e4))
                * r
                + T::from_f64(6.7265770927008700e4))
                * r
                + T::from_f64(4.5921953931549871e4))
                * r
                + T::from_f64(1.3731693765509461e4))
                * r
                + T::from_f64(1.9715909503065514e3))
                * r
                + T::from_f64(1.3314166789178437e2))
                * r
                + T::from_f64(3.3871328727963666e0);
            let den = ((((((T::from_f64(5.2264952788528545e3) * r
                + T::from_f64(2.8729085735721942e4))
                * r
                + T::from_f64(3.9307895800092710e4))
                * r
                + T::from_f64(2.1213794301586595e4))
                * r
                + T::from_f64(5.3941960214247511e3))
                * r
                + T::from_f64(6.8718700749205790e2))
                * r
                + T::from_f64(4.2313330701600911e1))
                * r
                + T::one();
            q * num / den
        } else {
            let r_tail = if q < T::zero() {
                p_clamped
            } else {
                T::one() - p_clamped
            };
            let s = (-r_tail.ln()).sqrt() - T::from_f64(1.6);
            let num = ((((((T::from_f64(7.7454501427834140e-4) * s
                + T::from_f64(2.2723844989269184e-2))
                * s
                + T::from_f64(2.4178072517745061e-1))
                * s
                + T::from_f64(1.2704582524523683e0))
                * s
                + T::from_f64(3.6478483247632046e0))
                * s
                + T::from_f64(5.7694972214606914e0))
                * s
                + T::from_f64(4.6303378461565452e0))
                * s
                + T::from_f64(1.4234371107496835e0);
            let den = ((((((T::from_f64(1.0507500716444168e-9) * s
                + T::from_f64(5.4759380849953449e-4))
                * s
                + T::from_f64(1.5198666563616457e-2))
                * s
                + T::from_f64(1.4810397642748007e-1))
                * s
                + T::from_f64(6.8976733498510000e-1))
                * s
                + T::from_f64(1.6763848301838038e0))
                * s
                + T::from_f64(2.0531916266377588e0))
                * s
                + T::one();
            let z = num / den;
            if q < T::zero() { -z } else { z }
        }
    }

    /// Evaluates the Rockafellar & Uryasev (2013) Mixed-CVaR Risk Quadrangle functional.
    pub fn evaluate_rockafellar_uryasev_mixed_cvar_risk_quadrangle<T: Real>(
        spec: &RockafellarUryasevMixedCvarSpec<T>,
    ) -> RockafellarUryasevMixedCvarAnalytics<T> {
        let sig = spec.loss_volatility_sigma.max_val(T::from_f64(1e-5));
        let a1 = spec
            .moderate_confidence_alpha1
            .max_val(T::from_f64(0.55))
            .min_val(T::from_f64(0.985));
        let a2 = spec
            .extreme_confidence_alpha2
            .max_val(T::from_f64(0.90))
            .min_val(T::from_f64(0.999));
        let lam1 = spec
            .moderate_spectral_weight_lambda1
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.99));
        let kap = spec.tail_kurtosis_amplification_kappa.max_val(T::zero());

        let z1 = inv_norm_cdf_tail(a1);
        let z2 = inv_norm_cdf_tail(a2);

        let es1_std = (norm_pdf(z1) / (T::one() - a1))
            * (T::one() + (kap / T::from_f64(24.0)) * (z1 * z1 * z1 - T::from_f64(3.0) * z1));
        let es2_std = (norm_pdf(z2) / (T::one() - a2))
            * (T::one() + (kap / T::from_f64(24.0)) * (z2 * z2 * z2 - T::from_f64(3.0) * z2));

        let cvar1 = spec.mean_loss_mu + sig * es1_std;
        let cvar2 = spec.mean_loss_mu + sig * es2_std;

        let mixed_dev = sig * (lam1 * es1_std + (T::one() - lam1) * es2_std);
        let mixed_risk = spec.mean_loss_mu + mixed_dev;

        RockafellarUryasevMixedCvarAnalytics {
            moderate_tail_cvar_alpha1: cvar1,
            extreme_tail_cvar_alpha2: cvar2,
            quadrangle_mixed_deviation: mixed_dev,
            quadrangle_mixed_cvar_risk: mixed_risk,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn rockafellar_uryasev_mixed_cvar_invariants_and_ad_parity() {
            let base = RockafellarUryasevMixedCvarSpec {
                mean_loss_mu: 0.012_f64,
                loss_volatility_sigma: 0.045,
                moderate_confidence_alpha1: 0.95,
                extreme_confidence_alpha2: 0.99,
                moderate_spectral_weight_lambda1: 0.60,
                tail_kurtosis_amplification_kappa: 1.50,
            };
            let r0 = evaluate_rockafellar_uryasev_mixed_cvar_risk_quadrangle(&base);
            assert!(r0.extreme_tail_cvar_alpha2 > r0.moderate_tail_cvar_alpha1);
            assert!(
                r0.quadrangle_mixed_cvar_risk > r0.moderate_tail_cvar_alpha1
                    && r0.quadrangle_mixed_cvar_risk < r0.extreme_tail_cvar_alpha2
            );
            assert!(
                (r0.quadrangle_mixed_cvar_risk
                    - (base.mean_loss_mu + r0.quadrangle_mixed_deviation))
                    .abs()
                    < 1e-12
            );

            let inputs = [
                base.mean_loss_mu,
                base.loss_volatility_sigma,
                base.moderate_confidence_alpha1,
                base.extreme_confidence_alpha2,
                base.moderate_spectral_weight_lambda1,
                base.tail_kurtosis_amplification_kappa,
            ];
            let mut dual_in = [Dual::<6>::from_f64(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = RockafellarUryasevMixedCvarSpec {
                mean_loss_mu: dual_in[0],
                loss_volatility_sigma: dual_in[1],
                moderate_confidence_alpha1: dual_in[2],
                extreme_confidence_alpha2: dual_in[3],
                moderate_spectral_weight_lambda1: dual_in[4],
                tail_kurtosis_amplification_kappa: dual_in[5],
            };
            let r_ad = evaluate_rockafellar_uryasev_mixed_cvar_risk_quadrangle(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| RockafellarUryasevMixedCvarSpec {
                    mean_loss_mu: v[0],
                    loss_volatility_sigma: v[1],
                    moderate_confidence_alpha1: v[2],
                    extreme_confidence_alpha2: v[3],
                    moderate_spectral_weight_lambda1: v[4],
                    tail_kurtosis_amplification_kappa: v[5],
                };
                let rp = evaluate_rockafellar_uryasev_mixed_cvar_risk_quadrangle(&mk(p));
                let rm = evaluate_rockafellar_uryasev_mixed_cvar_risk_quadrangle(&mk(m));
                let fd =
                    (rp.quadrangle_mixed_cvar_risk - rm.quadrangle_mixed_cvar_risk) / (2.0 * eps);
                assert!(
                    r_ad.quadrangle_mixed_cvar_risk.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.quadrangle_mixed_cvar_risk.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use rockafellar_uryasev_mixed_cvar_quadrangle_mod::*;

// ============================================================================
// Glasserman, Heidelberger & Shahabuddin (2000, 2002) Importance Sampling Credit/Delta-Gamma VaR
// ============================================================================

pub mod glasserman_importance_sampling_var_mod {
    use ql_ad::Real;
    use ql_math::distributions::{norm_cdf, norm_pdf};

    /// Specification for the Glasserman, Heidelberger & Shahabuddin (2000, *Management Science* 46(10);
    /// 2002, *Math. Finance* 12(3)) and Glasserman (2004, Ch. 9) exponential-twisting saddlepoint
    /// importance-sampling tail risk and Value-at-Risk engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GlassermanImportanceSamplingVaRSpec<T: Real> {
        /// First-order systemic factor sensitivity $a_{\text{delta}} > 0$.
        pub portfolio_linear_delta_exposure_a: T,
        /// Second-order quadratic convexity / concentration exposure $b_{\text{gamma}} \in (0, 0.45)$.
        pub portfolio_quadratic_gamma_exposure_b: T,
        /// Conditional idiosyncratic loss volatility $\sigma_\varepsilon > 0$.
        pub idiosyncratic_noise_volatility_sigma_e: T,
        /// Extreme tail loss threshold $x > 0$ for $P(L > x)$ saddlepoint twisting.
        pub tail_loss_threshold_x: T,
        /// Constant portfolio loss drift offset $c_0$.
        pub constant_loss_drift_c0: T,
        /// Value-at-Risk confidence level $\alpha \in (0.90, 0.999)$.
        pub confidence_level_alpha: T,
        /// Importance-sampling Monte Carlo path budget $N_{\text{IS}} > 0$.
        pub effective_sample_paths_n: T,
    }

    /// Analytical outputs for the Glasserman-Heidelberger-Shahabuddin importance-sampling VaR engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct GlassermanImportanceSamplingVaRAnalytics<T: Real> {
        /// Optimal exponential-twisting saddlepoint parameter $\theta^* > 0$ solving $\psi'(\theta^*) = x$.
        pub optimal_exponential_twisting_theta: T,
        /// Lugannani-Rice / saddlepoint tail exceedance probability $P(L > x)$.
        pub saddlepoint_tail_exceedance_probability: T,
        /// Importance-sampling estimator standard error $\text{SE}_{\text{IS}}(x)$ for $N_{\text{IS}}$ paths.
        pub importance_sampling_standard_error: T,
        /// Saddlepoint-adjusted quadratic Value-at-Risk $\text{VaR}_\alpha^{\text{IS}}$.
        pub saddlepoint_adjusted_quadratic_var: T,
    }

    /// Evaluates the Glasserman, Heidelberger & Shahabuddin (2000, Section 3; Glasserman 2004, Ch. 9.2)
    /// cumulant generating function $\psi(\theta)$, saddlepoint twisting parameter $\theta^*$,
    /// Lugannani-Rice tail probability, and importance-sampling-adjusted quadratic VaR.
    pub fn evaluate_glasserman_heidelberger_shahabuddin_importance_sampling_credit_var<T: Real>(
        spec: &GlassermanImportanceSamplingVaRSpec<T>,
    ) -> GlassermanImportanceSamplingVaRAnalytics<T> {
        let a = spec
            .portfolio_linear_delta_exposure_a
            .max_val(T::from_f64(0.1));
        let b = spec
            .portfolio_quadratic_gamma_exposure_b
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.45));
        let sig_e = spec
            .idiosyncratic_noise_volatility_sigma_e
            .max_val(T::from_f64(0.05));
        let x_excess = (spec.tail_loss_threshold_x - spec.constant_loss_drift_c0)
            .max_val(b + T::from_f64(0.5));

        // Exact root of E_theta[a*Z + b*Z^2 + sig_e*eps] = x_excess via 4 Newton steps on psi'(theta) = x
        let theta_max = T::from_f64(0.44) / b;
        let mut theta = ((x_excess - b)
            / (a * a + sig_e * sig_e + T::from_f64(2.0) * b * x_excess))
            .max_val(T::from_f64(0.02))
            .min_val(theta_max);

        for _ in 0..4 {
            let denom = (T::one() - T::from_f64(2.0) * theta * b).max_val(T::from_f64(0.08));
            let psi_p = b / denom
                + (theta * a * a * (T::one() - theta * b)) / (denom * denom)
                + theta * sig_e * sig_e;
            let psi_pp = (T::from_f64(2.0) * b * b) / (denom * denom)
                + (a * a) / (denom * denom * denom)
                + sig_e * sig_e;
            theta = (theta - (psi_p - x_excess) / psi_pp)
                .max_val(T::from_f64(0.01))
                .min_val(theta_max);
        }

        let denom = (T::one() - T::from_f64(2.0) * theta * b).max_val(T::from_f64(0.08));
        let psi = -T::from_f64(0.5) * denom.ln()
            + (theta * theta * a * a) / (T::from_f64(2.0) * denom)
            + T::from_f64(0.5) * theta * theta * sig_e * sig_e;
        let psi_pp = (T::from_f64(2.0) * b * b) / (denom * denom)
            + (a * a) / (denom * denom * denom)
            + sig_e * sig_e;

        let rate_fn = (theta * x_excess - psi).max_val(T::from_f64(1e-6));
        let w_lr = (T::from_f64(2.0) * rate_fn).sqrt();
        let u_lr = theta * psi_pp.sqrt();
        let tail_prob = (T::one() - norm_cdf(w_lr)
            + norm_pdf(w_lr) * (T::one() / u_lr - T::one() / w_lr))
            .max_val(T::from_f64(1e-9));

        // Asymptotic second moment M_2(theta*) of the importance-sampling likelihood ratio estimator
        let n_paths = spec.effective_sample_paths_n.max_val(T::one());
        let m2_is = (-T::from_f64(2.0) * rate_fn).exp()
            / (T::from_f64(2.0) * T::from_f64(std::f64::consts::PI).sqrt() * u_lr);
        let se_is = (m2_is / n_paths).sqrt();

        // Smooth normal quantile approximation for alpha in (0.90, 0.999) + saddlepoint tail coupling
        let p_tail = (T::one() - spec.confidence_level_alpha)
            .max_val(T::from_f64(1e-4))
            .min_val(T::from_f64(0.15));
        let t_q = (-T::from_f64(2.0) * p_tail.ln()).sqrt();
        let z_alpha = t_q
            - (T::from_f64(2.515517)
                + T::from_f64(0.802853) * t_q
                + T::from_f64(0.010328) * t_q * t_q)
                / (T::one()
                    + T::from_f64(1.432788) * t_q
                    + T::from_f64(0.189269) * t_q * t_q
                    + T::from_f64(0.001308) * t_q * t_q * t_q);

        let var_quad = spec.constant_loss_drift_c0
            + a * z_alpha
            + b * z_alpha * z_alpha
            + (sig_e * sig_e * z_alpha) / (T::from_f64(2.0) * (a * a + sig_e * sig_e).sqrt())
            + T::from_f64(10.0) * se_is;

        GlassermanImportanceSamplingVaRAnalytics {
            optimal_exponential_twisting_theta: theta,
            saddlepoint_tail_exceedance_probability: tail_prob,
            importance_sampling_standard_error: se_is,
            saddlepoint_adjusted_quadratic_var: var_quad,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn glasserman_importance_sampling_var_invariants_and_ad_parity() {
            let base = GlassermanImportanceSamplingVaRSpec {
                portfolio_linear_delta_exposure_a: 1.80_f64,
                portfolio_quadratic_gamma_exposure_b: 0.22,
                idiosyncratic_noise_volatility_sigma_e: 0.65,
                tail_loss_threshold_x: 5.20,
                constant_loss_drift_c0: 0.15,
                confidence_level_alpha: 0.99,
                effective_sample_paths_n: 1000.0,
            };
            let r0 =
                evaluate_glasserman_heidelberger_shahabuddin_importance_sampling_credit_var(&base);
            assert!(r0.optimal_exponential_twisting_theta > 0.0);
            assert!(
                r0.saddlepoint_tail_exceedance_probability > 0.0
                    && r0.saddlepoint_tail_exceedance_probability < 0.2
            );
            assert!(r0.importance_sampling_standard_error > 0.0);
            assert!(r0.saddlepoint_adjusted_quadratic_var > base.tail_loss_threshold_x);

            let inputs = [
                base.portfolio_linear_delta_exposure_a,
                base.portfolio_quadratic_gamma_exposure_b,
                base.idiosyncratic_noise_volatility_sigma_e,
                base.tail_loss_threshold_x,
                base.constant_loss_drift_c0,
                base.confidence_level_alpha,
                base.effective_sample_paths_n,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GlassermanImportanceSamplingVaRSpec {
                portfolio_linear_delta_exposure_a: dual_in[0],
                portfolio_quadratic_gamma_exposure_b: dual_in[1],
                idiosyncratic_noise_volatility_sigma_e: dual_in[2],
                tail_loss_threshold_x: dual_in[3],
                constant_loss_drift_c0: dual_in[4],
                confidence_level_alpha: dual_in[5],
                effective_sample_paths_n: dual_in[6],
            };
            let r_ad = evaluate_glasserman_heidelberger_shahabuddin_importance_sampling_credit_var(
                &spec_ad,
            );
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| GlassermanImportanceSamplingVaRSpec {
                    portfolio_linear_delta_exposure_a: v[0],
                    portfolio_quadratic_gamma_exposure_b: v[1],
                    idiosyncratic_noise_volatility_sigma_e: v[2],
                    tail_loss_threshold_x: v[3],
                    constant_loss_drift_c0: v[4],
                    confidence_level_alpha: v[5],
                    effective_sample_paths_n: v[6],
                };
                let rp =
                    evaluate_glasserman_heidelberger_shahabuddin_importance_sampling_credit_var(
                        &mk(p),
                    );
                let rm =
                    evaluate_glasserman_heidelberger_shahabuddin_importance_sampling_credit_var(
                        &mk(m),
                    );
                let fd = (rp.saddlepoint_adjusted_quadratic_var
                    - rm.saddlepoint_adjusted_quadratic_var)
                    / (2.0 * eps);
                assert!(
                    r_ad.saddlepoint_adjusted_quadratic_var.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.saddlepoint_adjusted_quadratic_var.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use glasserman_importance_sampling_var_mod::*;
mod gls_tasche_cornish_fisher_raroc_mod {
    use ql_ad::Real;

    /// Input specification for the Gouriéroux, Laurent & Scaillet (2000, *J. Empirical Finance*) /
    /// Tasche (1999; 2008) Cornish–Fisher Modified VaR Euler decomposition and RAROC engine.
    #[derive(Debug, Clone, Copy)]
    pub struct GlsTascheCornishFisherRarocSpec<T: Real> {
        /// Portfolio expected excess return $\mu_P$ (e.g., $0.082$).
        pub portfolio_expected_excess_return_mu_p: T,
        /// Portfolio volatility $\sigma_P$ (e.g., $0.145$).
        pub portfolio_volatility_sigma_p: T,
        /// Portfolio return skewness $S_P$ (negative for crash risk, e.g., $-0.65$).
        pub portfolio_skewness_s_p: T,
        /// Portfolio excess kurtosis $K_P$ (e.g., $2.40$).
        pub portfolio_excess_kurtosis_k_p: T,
        /// Standard normal tail quantile $z_\alpha > 0$ (e.g., $2.326347874$ for 99% VaR).
        pub normal_tail_quantile_z_alpha: T,
        /// Subportfolio expected excess return $\mu_i$ (e.g., $0.105$).
        pub subportfolio_expected_excess_return_mu_i: T,
        /// Subportfolio return volatility $\sigma_i$ (e.g., $0.180$).
        pub subportfolio_volatility_sigma_i: T,
        /// Correlation $\rho_{i,P}$ between subportfolio $i$ and total portfolio $P$ (e.g., $0.72$).
        pub subportfolio_portfolio_correlation_rho: T,
    }

    /// Analytics returned by [`evaluate_gourieroux_laurent_scaillet_tasche_cornish_fisher_euler_var_raroc`].
    #[derive(Debug, Clone, Copy)]
    pub struct GlsTascheCornishFisherRarocAnalytics<T: Real> {
        /// Cornish–Fisher four-moment adjusted tail quantile $z_{\text{CF}}(\alpha)$.
        pub cornish_fisher_adjusted_quantile: T,
        /// Portfolio Modified Value-at-Risk $\text{mVaR}_\alpha(P) = -\mu_P + z_{\text{CF}}(\alpha)\sigma_P$.
        pub portfolio_modified_var: T,
        /// Marginal Modified VaR $\partial \text{mVaR}_\alpha / \partial w_i$ of subportfolio $i$.
        pub subportfolio_marginal_modified_var: T,
        /// Tasche (1999) Euler Risk-Adjusted Return on Capital $\text{RAROC}_i = \mu_i / (\partial \text{mVaR}_\alpha / \partial w_i)$.
        pub subportfolio_euler_raroc: T,
    }

    /// Evaluates the Favre & Galeano (2002) / Gouriéroux, Laurent & Scaillet (2000) / Tasche (1999)
    /// Cornish–Fisher Modified Value-at-Risk, marginal risk contribution $\partial \text{mVaR}_\alpha / \partial w_i$,
    /// and Euler-consistent RAROC for a skewed, leptokurtic subportfolio.
    pub fn evaluate_gourieroux_laurent_scaillet_tasche_cornish_fisher_euler_var_raroc<T: Real>(
        spec: &GlsTascheCornishFisherRarocSpec<T>,
    ) -> GlsTascheCornishFisherRarocAnalytics<T> {
        let z = spec.normal_tail_quantile_z_alpha;
        let z2 = z * z;
        let z3 = z2 * z;
        // For left-tail loss (-R_P), loss skewness is -S_P > 0 when return skewness S_P < 0
        let skew_loss = -spec.portfolio_skewness_s_p;
        let kurt = spec.portfolio_excess_kurtosis_k_p;

        let z_cf = z
            + (z2 - T::one()) / T::from_f64(6.0) * skew_loss
            + (z3 - T::from_f64(3.0) * z) / T::from_f64(24.0) * kurt
            - (T::from_f64(2.0) * z3 - T::from_f64(5.0) * z) / T::from_f64(36.0)
                * skew_loss
                * skew_loss;

        let mvar_p =
            -spec.portfolio_expected_excess_return_mu_p + z_cf * spec.portfolio_volatility_sigma_p;

        // Portfolio Cornish-Fisher capital intensity factor mVaR_P / (z_alpha * sigma_P)
        let port_capital_factor = mvar_p / (z * spec.portfolio_volatility_sigma_p);
        // Marginal volatility contribution d(sigma_P)/dw_i = rho_{i,P} * sigma_i
        let marg_sigma_i =
            spec.subportfolio_portfolio_correlation_rho * spec.subportfolio_volatility_sigma_i;
        let marg_mvar_i = (-spec.subportfolio_expected_excess_return_mu_i + z_cf * marg_sigma_i)
            * port_capital_factor;

        let raroc_i = spec.subportfolio_expected_excess_return_mu_i / marg_mvar_i;

        GlsTascheCornishFisherRarocAnalytics {
            cornish_fisher_adjusted_quantile: z_cf,
            portfolio_modified_var: mvar_p,
            subportfolio_marginal_modified_var: marg_mvar_i,
            subportfolio_euler_raroc: raroc_i,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gls_tasche_cornish_fisher_raroc_invariants_and_ad_parity() {
            let base = GlsTascheCornishFisherRarocSpec {
                portfolio_expected_excess_return_mu_p: 0.082_f64,
                portfolio_volatility_sigma_p: 0.145_f64,
                portfolio_skewness_s_p: -0.65_f64,
                portfolio_excess_kurtosis_k_p: 2.40_f64,
                normal_tail_quantile_z_alpha: 2.326_347_874_f64,
                subportfolio_expected_excess_return_mu_i: 0.105_f64,
                subportfolio_volatility_sigma_i: 0.180_f64,
                subportfolio_portfolio_correlation_rho: 0.72_f64,
            };
            let res =
                evaluate_gourieroux_laurent_scaillet_tasche_cornish_fisher_euler_var_raroc(&base);
            assert!(res.cornish_fisher_adjusted_quantile > base.normal_tail_quantile_z_alpha);
            assert!(res.portfolio_modified_var > 0.20);
            assert!(res.subportfolio_marginal_modified_var > 0.10);
            assert!(res.subportfolio_euler_raroc > 0.10);

            let inputs = [
                base.portfolio_expected_excess_return_mu_p,
                base.portfolio_volatility_sigma_p,
                base.portfolio_skewness_s_p,
                base.portfolio_excess_kurtosis_k_p,
                base.normal_tail_quantile_z_alpha,
                base.subportfolio_expected_excess_return_mu_i,
                base.subportfolio_volatility_sigma_i,
                base.subportfolio_portfolio_correlation_rho,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GlsTascheCornishFisherRarocSpec {
                portfolio_expected_excess_return_mu_p: dual_in[0],
                portfolio_volatility_sigma_p: dual_in[1],
                portfolio_skewness_s_p: dual_in[2],
                portfolio_excess_kurtosis_k_p: dual_in[3],
                normal_tail_quantile_z_alpha: dual_in[4],
                subportfolio_expected_excess_return_mu_i: dual_in[5],
                subportfolio_volatility_sigma_i: dual_in[6],
                subportfolio_portfolio_correlation_rho: dual_in[7],
            };
            let r_ad = evaluate_gourieroux_laurent_scaillet_tasche_cornish_fisher_euler_var_raroc(
                &spec_ad,
            );
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| GlsTascheCornishFisherRarocSpec {
                    portfolio_expected_excess_return_mu_p: v[0],
                    portfolio_volatility_sigma_p: v[1],
                    portfolio_skewness_s_p: v[2],
                    portfolio_excess_kurtosis_k_p: v[3],
                    normal_tail_quantile_z_alpha: v[4],
                    subportfolio_expected_excess_return_mu_i: v[5],
                    subportfolio_volatility_sigma_i: v[6],
                    subportfolio_portfolio_correlation_rho: v[7],
                };
                let rp = evaluate_gourieroux_laurent_scaillet_tasche_cornish_fisher_euler_var_raroc(
                    &mk(p),
                );
                let rm = evaluate_gourieroux_laurent_scaillet_tasche_cornish_fisher_euler_var_raroc(
                    &mk(m),
                );
                let fd = (rp.subportfolio_euler_raroc - rm.subportfolio_euler_raroc) / (2.0 * eps);
                assert!(
                    r_ad.subportfolio_euler_raroc.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.subportfolio_euler_raroc.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use gls_tasche_cornish_fisher_raroc_mod::*;
mod cont_deguest_scandolo_robust_risk_mod {
    use ql_ad::Real;

    /// Input specification for the Cont, Deguest & Scandolo (2010, *Mathematical Finance* 20(4))
    /// & Kou, Peng & Heyde (2013, *Mathematics of Operations Research* 38(3)) Hampel influence
    /// function and tail robustness comparison of VaR, Expected Shortfall, and Median Shortfall ($\text{MS}_\alpha$).
    #[derive(Debug, Clone, Copy)]
    pub struct ContDeguestScandoloRobustRiskSpec<T: Real> {
        /// Portfolio location loss parameter $\mu_L$ (e.g., $0.012$).
        pub loss_location_mu: T,
        /// Generalized Pareto / heavy-tail scale parameter $s_L > 0$ (e.g., $0.025$).
        pub tail_scale_s: T,
        /// Extreme-value tail shape index $\xi \in (0, 0.5)$ (e.g., $0.24$).
        pub tail_index_xi: T,
        /// Tail tail-probability level $p = 1 - \alpha \in (0, 0.10)$ (e.g., $0.025$ for $97.5\%$ ES/MS).
        pub tail_probability_p: T,
        /// Dirac point-mass contamination weight $\varepsilon_{\text{cont}} > 0$ (e.g., $0.002$).
        pub contamination_epsilon: T,
        /// Extreme stress outlier loss level $z_{\text{outlier}} > \text{VaR}_{1-p}$ (e.g., $0.28$).
        pub stress_outlier_loss_z: T,
    }

    /// Analytics returned by [`evaluate_cont_deguest_scandolo_tail_risk_influence_sensitivity`].
    #[derive(Debug, Clone, Copy)]
    pub struct ContDeguestScandoloRobustRiskAnalytics<T: Real> {
        /// Baseline Generalized Pareto Value-at-Risk $\text{VaR}_{1-p}$.
        pub baseline_value_at_risk: T,
        /// Baseline Kou–Peng–Heyde Median Shortfall $\text{MS}_{1-p} = \text{VaR}_{1-p/2}$.
        pub baseline_median_shortfall: T,
        /// Baseline Expected Shortfall $\text{ES}_{1-p}$.
        pub baseline_expected_shortfall: T,
        /// Contamination-perturbed Median Shortfall vs. ES robustness ratio $\text{MS}_{1-p}^{(\varepsilon)} / \text{ES}_{1-p}^{(\varepsilon)}$.
        pub contaminated_median_to_expected_shortfall_ratio: T,
    }

    /// Evaluates the Cont, Deguest & Scandolo (2010, Proposition 3.2) Hampel influence functions
    /// for $\text{VaR}_{1-p}$, Median Shortfall $\text{MS}_{1-p} = \text{VaR}_{1-p/2}$ (Kou, Peng & Heyde 2013),
    /// and Expected Shortfall $\text{ES}_{1-p}$ under a heavy-tailed Generalized Pareto loss law:
    /// - $\text{IF}(z; \text{MS}_{1-p}) = \frac{1 - p/2}{f(\text{VaR}_{1-p/2})}$ is **bounded** in $z$, whereas
    /// - $\text{IF}(z; \text{ES}_{1-p}) = \text{VaR}_{1-p} - \text{ES}_{1-p} + \frac{z - \text{VaR}_{1-p}}{p}$ grows **linearly** with the outlier $z$.
    pub fn evaluate_cont_deguest_scandolo_tail_risk_influence_sensitivity<T: Real>(
        spec: &ContDeguestScandoloRobustRiskSpec<T>,
    ) -> ContDeguestScandoloRobustRiskAnalytics<T> {
        let mu = spec.loss_location_mu;
        let s = spec.tail_scale_s;
        let xi = spec.tail_index_xi;
        let p = spec.tail_probability_p;
        let eps = spec.contamination_epsilon;
        let z_out = spec.stress_outlier_loss_z;

        // Generalized Pareto quantile VaR_{1-p} = mu + (s / xi) * (p^(-xi) - 1)
        let p_neg_xi = p.pow(-xi);
        let var_p = mu + (s / xi) * (p_neg_xi - T::one());

        // Median Shortfall MS_{1-p} = VaR_{1 - p/2} (Kou, Peng & Heyde 2013)
        let p_half = T::from_f64(0.5) * p;
        let p_half_neg_xi = p_half.pow(-xi);
        let ms_p = mu + (s / xi) * (p_half_neg_xi - T::one());

        // Expected Shortfall ES_{1-p} = (VaR_{1-p} + s - xi*mu) / (1 - xi)
        let es_p = (var_p + s - xi * mu) / (T::one() - xi);

        // Hampel influence function at outlier z_out > MS_{1-p}:
        // Density at MS_{1-p} is f(MS_{1-p}) = (1 / s) * (p / 2)^(1 + xi)
        // IF(z; MS_{1-p}) = (1 - p/2) / f(MS_{1-p})
        let inv_f_ms = s * p_half.pow(-(T::one() + xi));
        let if_ms = (T::one() - p_half) * inv_f_ms;

        // IF(z_out; ES_{1-p}) = VaR_{1-p} - ES_{1-p} + (z_out - VaR_{1-p}) / p (Cont, Deguest & Scandolo 2010)
        let if_es = var_p - es_p + (z_out - var_p) / p;

        let ms_cont = ms_p + eps * if_ms;
        let es_cont = es_p + eps * if_es;
        let rob_ratio = ms_cont / es_cont;

        ContDeguestScandoloRobustRiskAnalytics {
            baseline_value_at_risk: var_p,
            baseline_median_shortfall: ms_p,
            baseline_expected_shortfall: es_p,
            contaminated_median_to_expected_shortfall_ratio: rob_ratio,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn cont_deguest_scandolo_robust_risk_invariants_and_ad_parity() {
            let base = ContDeguestScandoloRobustRiskSpec {
                loss_location_mu: 0.012_f64,
                tail_scale_s: 0.025_f64,
                tail_index_xi: 0.24_f64,
                tail_probability_p: 0.025_f64,
                contamination_epsilon: 0.002_f64,
                stress_outlier_loss_z: 0.28_f64,
            };
            let res = evaluate_cont_deguest_scandolo_tail_risk_influence_sensitivity(&base);
            assert!(res.baseline_median_shortfall > res.baseline_value_at_risk);
            assert!(res.baseline_expected_shortfall > res.baseline_median_shortfall);
            assert!(
                res.contaminated_median_to_expected_shortfall_ratio > 0.5
                    && res.contaminated_median_to_expected_shortfall_ratio < 1.0
            );

            let inputs = [
                base.loss_location_mu,
                base.tail_scale_s,
                base.tail_index_xi,
                base.tail_probability_p,
                base.contamination_epsilon,
                base.stress_outlier_loss_z,
            ];
            let mut dual_in = [Dual::<6>::from_value(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = ContDeguestScandoloRobustRiskSpec {
                loss_location_mu: dual_in[0],
                tail_scale_s: dual_in[1],
                tail_index_xi: dual_in[2],
                tail_probability_p: dual_in[3],
                contamination_epsilon: dual_in[4],
                stress_outlier_loss_z: dual_in[5],
            };
            let r_ad = evaluate_cont_deguest_scandolo_tail_risk_influence_sensitivity(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| ContDeguestScandoloRobustRiskSpec {
                    loss_location_mu: v[0],
                    tail_scale_s: v[1],
                    tail_index_xi: v[2],
                    tail_probability_p: v[3],
                    contamination_epsilon: v[4],
                    stress_outlier_loss_z: v[5],
                };
                let rp = evaluate_cont_deguest_scandolo_tail_risk_influence_sensitivity(&mk(p));
                let rm = evaluate_cont_deguest_scandolo_tail_risk_influence_sensitivity(&mk(m));
                let fd = (rp.contaminated_median_to_expected_shortfall_ratio
                    - rm.contaminated_median_to_expected_shortfall_ratio)
                    / (2.0 * eps);
                assert!(
                    r_ad.contaminated_median_to_expected_shortfall_ratio.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!(
                    (r_ad.contaminated_median_to_expected_shortfall_ratio.d[i] - fd).abs() < 1e-5
                );
            }
        }
    }
}
pub use cont_deguest_scandolo_robust_risk_mod::*;
mod mcneil_frey_evt_garch_mod {
    use ql_ad::Real;

    /// Input specification for the McNeil & Frey (2000, *Journal of Empirical Finance* 7(3),
    /// *"Estimation of Tail-Related Risk Measures for Heteroscedastic Financial Time Series"*)
    /// two-stage GARCH-filtered Peaks-Over-Threshold (POT) Generalized Pareto conditional
    /// $\text{VaR}$ and Expected Shortfall engine.
    #[derive(Debug, Clone, Copy)]
    pub struct McNeilFreyConditionalEvtGarchRiskSpec<T: Real> {
        /// Conditional 1-step loss drift forecast $\mu_{t+1}$ (e.g., $0.0005$).
        pub conditional_loss_mean_mu: T,
        /// GARCH$(1,1)$ unconditional variance intercept $\omega > 0$ (e.g., $0.000012$).
        pub garch_omega: T,
        /// GARCH$(1,1)$ squared-shock reaction coefficient $\alpha_{\text{ARCH}} \times \epsilon_t^2$ (e.g., $0.00018$).
        pub garch_arch_shock_term: T,
        /// GARCH$(1,1)$ lagged conditional variance term $\beta_{\text{GARCH}} \times \sigma_t^2$ (e.g., $0.00022$).
        pub garch_lagged_variance_term: T,
        /// Standardized residual POT threshold $u_z > 0$ (e.g., $1.65$).
        pub standardized_pot_threshold_u: T,
        /// Tail exceedance fraction $p_u = N_u / n \in (1 - q, 1)$ above $u_z$ (e.g., $0.060$).
        pub threshold_exceedance_rate_pu: T,
        /// GPD tail shape index $\xi \in (0, 1)$ on standardized residuals (e.g., $0.22$).
        pub gpd_shape_index_xi: T,
        /// GPD scale parameter $\beta_{\text{GPD}} > 0$ on standardized residuals (e.g., $0.52$).
        pub gpd_scale_beta: T,
    }

    /// Analytics returned by [`evaluate_mcneil_frey_garch_filtered_evt_pot_tail_risk`].
    #[derive(Debug, Clone, Copy)]
    pub struct McNeilFreyConditionalEvtGarchRiskAnalytics<T: Real> {
        /// 1-step-ahead GARCH conditional volatility forecast $\sigma_{t+1}$.
        pub conditional_garch_volatility: T,
        /// Standardized 99% EVT-POT innovation quantile $z_{0.99}$.
        pub standardized_evt_quantile_z: T,
        /// Dynamic conditional 99% Value-at-Risk $\text{VaR}_{t+1}^{0.99} = \mu_{t+1} + \sigma_{t+1} z_{0.99}$.
        pub conditional_evt_garch_var99: T,
        /// Dynamic conditional 99% Expected Shortfall $\text{ES}_{t+1}^{0.99} = \mu_{t+1} + \sigma_{t+1} s_{0.99}$.
        pub conditional_evt_garch_es99: T,
    }

    /// Evaluates the McNeil & Frey (2000, Eqs. 11–13) two-stage GARCH-EVT conditional tail
    /// risk estimator at confidence level $q = 0.99$:
    /// $$z_q = u_z + \frac{\beta_{\text{GPD}}}{\xi}\left(\left(\frac{1 - q}{p_u}\right)^{-\xi} - 1\right), \qquad s_q = \frac{z_q}{1 - \xi} + \frac{\beta_{\text{GPD}} - \xi u_z}{1 - \xi},$$
    /// scaling standardized tail risk by the 1-step GARCH volatility $\sigma_{t+1} = \sqrt{\omega + \alpha\epsilon_t^2 + \beta\sigma_t^2}$.
    pub fn evaluate_mcneil_frey_garch_filtered_evt_pot_tail_risk<T: Real>(
        spec: &McNeilFreyConditionalEvtGarchRiskSpec<T>,
    ) -> McNeilFreyConditionalEvtGarchRiskAnalytics<T> {
        let sig_next =
            (spec.garch_omega + spec.garch_arch_shock_term + spec.garch_lagged_variance_term)
                .sqrt();
        let u_z = spec.standardized_pot_threshold_u;
        let p_u = spec.threshold_exceedance_rate_pu;
        let xi = spec.gpd_shape_index_xi;
        let beta = spec.gpd_scale_beta;

        // 99% tail probability 1 - q = 0.01
        let tail_ratio = p_u / T::from_f64(0.01);
        let z_q = u_z + (beta / xi) * (tail_ratio.pow(xi) - T::one());
        let s_q = (z_q + beta - xi * u_z) / (T::one() - xi);

        let var99 = spec.conditional_loss_mean_mu + sig_next * z_q;
        let es99 = spec.conditional_loss_mean_mu + sig_next * s_q;

        McNeilFreyConditionalEvtGarchRiskAnalytics {
            conditional_garch_volatility: sig_next,
            standardized_evt_quantile_z: z_q,
            conditional_evt_garch_var99: var99,
            conditional_evt_garch_es99: es99,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn mcneil_frey_evt_garch_invariants_and_ad_parity() {
            let base = McNeilFreyConditionalEvtGarchRiskSpec {
                conditional_loss_mean_mu: 0.0005_f64,
                garch_omega: 0.000012_f64,
                garch_arch_shock_term: 0.00018_f64,
                garch_lagged_variance_term: 0.00022_f64,
                standardized_pot_threshold_u: 1.65_f64,
                threshold_exceedance_rate_pu: 0.060_f64,
                gpd_shape_index_xi: 0.22_f64,
                gpd_scale_beta: 0.52_f64,
            };
            let res = evaluate_mcneil_frey_garch_filtered_evt_pot_tail_risk(&base);
            assert!(res.conditional_garch_volatility > 0.015);
            assert!(res.standardized_evt_quantile_z > base.standardized_pot_threshold_u);
            assert!(res.conditional_evt_garch_es99 > res.conditional_evt_garch_var99);

            let inputs = [
                base.conditional_loss_mean_mu,
                base.garch_omega,
                base.garch_arch_shock_term,
                base.garch_lagged_variance_term,
                base.standardized_pot_threshold_u,
                base.threshold_exceedance_rate_pu,
                base.gpd_shape_index_xi,
                base.gpd_scale_beta,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = McNeilFreyConditionalEvtGarchRiskSpec {
                conditional_loss_mean_mu: dual_in[0],
                garch_omega: dual_in[1],
                garch_arch_shock_term: dual_in[2],
                garch_lagged_variance_term: dual_in[3],
                standardized_pot_threshold_u: dual_in[4],
                threshold_exceedance_rate_pu: dual_in[5],
                gpd_shape_index_xi: dual_in[6],
                gpd_scale_beta: dual_in[7],
            };
            let r_ad = evaluate_mcneil_frey_garch_filtered_evt_pot_tail_risk(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| McNeilFreyConditionalEvtGarchRiskSpec {
                    conditional_loss_mean_mu: v[0],
                    garch_omega: v[1],
                    garch_arch_shock_term: v[2],
                    garch_lagged_variance_term: v[3],
                    standardized_pot_threshold_u: v[4],
                    threshold_exceedance_rate_pu: v[5],
                    gpd_shape_index_xi: v[6],
                    gpd_scale_beta: v[7],
                };
                let rp = evaluate_mcneil_frey_garch_filtered_evt_pot_tail_risk(&mk(p));
                let rm = evaluate_mcneil_frey_garch_filtered_evt_pot_tail_risk(&mk(m));
                let fd =
                    (rp.conditional_evt_garch_es99 - rm.conditional_evt_garch_es99) / (2.0 * eps);
                assert!(
                    r_ad.conditional_evt_garch_es99.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.conditional_evt_garch_es99.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use mcneil_frey_evt_garch_mod::*;
