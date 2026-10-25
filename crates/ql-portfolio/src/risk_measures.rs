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
