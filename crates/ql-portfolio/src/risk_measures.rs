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
