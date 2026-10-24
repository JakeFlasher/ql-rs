//! Structured-credit loss models.
//!
//! The implementations in this module deliberately keep the quadrature nodes
//! in `f64` (the integration measure is deterministic) while evaluating all
//! model parameters and loss amounts in `T`.  This makes the routines useful
//! for both ordinary pricing and forward/reverse automatic differentiation.

use ql_ad::Real;
use std::sync::LazyLock;

/// Standard-normal inverse evaluated in the AD domain using two Newton steps.
fn inv_norm_cdf_ad<T: Real>(p: T) -> T {
    let mut x = T::from_f64(ql_math::inv_norm_cdf(p.val().clamp(1e-14, 1.0 - 1e-14)));
    for _ in 0..2 {
        x -= (x.norm_cdf() - p) / x.norm_pdf().max_val(T::from_f64(1e-15));
    }
    x
}

/// 32-point Gauss--Hermite nodes and weights for a standard normal integral.
/// The values are generated once per call using the stable Newton recurrence
/// for the physicists' Hermite polynomial.  Mapping `x -> sqrt(2)x` and
/// dividing the weights by `sqrt(pi)` gives `E[f(Z)]`, `Z~N(0,1)`.
static NORMAL_HERMITE_32: LazyLock<([f64; 32], [f64; 32])> = LazyLock::new(|| {
    // Golub--Welsch for the Hermite Jacobi matrix.  Its eigenvalues are the
    // roots of the physicists' H_32 and the squared first eigenvector entries
    // are already normalized weights for a standard-normal expectation.
    let n = 32usize;
    let mut jacobi = nalgebra::DMatrix::<f64>::zeros(n, n);
    for i in 0..n - 1 {
        let off = ((i + 1) as f64 / 2.0).sqrt();
        jacobi[(i, i + 1)] = off;
        jacobi[(i + 1, i)] = off;
    }
    let eig = jacobi.symmetric_eigen();
    let mut pairs: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            (
                eig.eigenvalues[i] * 2.0_f64.sqrt(),
                eig.eigenvectors[(0, i)].powi(2),
            )
        })
        .collect();
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut nodes = [0.0; 32];
    let mut weights = [0.0; 32];
    for (i, (node, weight)) in pairs.into_iter().enumerate() {
        nodes[i] = node;
        weights[i] = weight;
    }
    (nodes, weights)
});

fn normal_hermite_32() -> &'static ([f64; 32], [f64; 32]) {
    &NORMAL_HERMITE_32
}

fn high_corr_normal_rule() -> (Vec<f64>, Vec<f64>) {
    // A wider, deterministic Gauss--Legendre rule resolves the narrow
    // conditional-loss transition that appears as rho approaches one.  It is
    // used only in that regime; ordinary calibrations retain the cached
    // 32-node Hermite rule.
    let n = 192usize;
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..n.div_ceil(2) {
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..32 {
            let mut p0 = 1.0;
            let mut p1 = z;
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            let dp = n as f64 * (z * p1 - p0) / (z * z - 1.0);
            z -= p1 / dp;
        }
        let mut p0 = 1.0;
        let mut p1 = z;
        for k in 2..=n {
            let p2 = ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64;
            p0 = p1;
            p1 = p2;
        }
        let dp = n as f64 * (z * p1 - p0) / (z * z - 1.0);
        let ww = 2.0 / ((1.0 - z * z) * dp * dp);
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = ww;
        w[n - 1 - i] = ww;
    }
    let inv_sqrt_2pi = 1.0 / std::f64::consts::TAU.sqrt();
    for i in 0..n {
        let z = 8.0 * x[i];
        x[i] = z;
        w[i] *= 8.0 * (-0.5 * z * z).exp() * inv_sqrt_2pi;
    }
    (x, w)
}

/// Vasicek large-homogeneous-portfolio expected tranche loss.  The return is
/// a loss fraction of the total portfolio notional (not divided by tranche
/// thickness), which is the convention used by QuantLib's LHP model.
pub fn gaussian_lhp_expected_tranche_loss<T: Real>(
    default_prob: T,
    recovery: T,
    correlation: T,
    attachment: T,
    detachment: T,
) -> T {
    let (nodes, weights) = normal_hermite_32();
    let p_inv = inv_norm_cdf_ad(default_prob);
    let rho = correlation
        .max_val(T::zero())
        .min_val(T::from_f64(1.0 - 1e-14));
    let sr = rho.sqrt();
    let den = (T::one() - rho).sqrt();
    let loss_unit = T::one() - recovery;
    let width = (detachment - attachment).max_val(T::zero());
    let mut out = T::zero();
    for (&z, &w) in nodes.iter().zip(weights.iter()) {
        let cond_p = ((p_inv - sr * z) / den).norm_cdf();
        let pool_loss = loss_unit * cond_p;
        let tranche = (pool_loss - attachment).max_val(T::zero()).min_val(width);
        out += tranche * w;
    }
    out
}

/// LHP loss as a fraction of tranche notional (rather than pool notional).
pub fn gaussian_lhp_expected_tranche_loss_fraction<T: Real>(
    default_prob: T,
    recovery: T,
    correlation: T,
    attachment: T,
    detachment: T,
) -> T {
    let width = (detachment - attachment).max_val(T::from_f64(1e-15));
    gaussian_lhp_expected_tranche_loss(default_prob, recovery, correlation, attachment, detachment)
        / width
}

/// LHP expected loss in currency units for a pool with `total_notional`.
pub fn gaussian_lhp_expected_tranche_loss_notional<T: Real>(
    default_prob: T,
    recovery: T,
    correlation: T,
    attachment: T,
    detachment: T,
    total_notional: T,
) -> T {
    total_notional
        * gaussian_lhp_expected_tranche_loss(
            default_prob,
            recovery,
            correlation,
            attachment,
            detachment,
        )
}

/// Exact Andersen--Sidenius--Basu recursion for a one-factor Gaussian loss
/// distribution. `loss_units[k]` is the integer loss size of issuer `k` in
/// units of `unit_notional`; attachment and detachment are monetary losses.
pub fn asb_recursive_expected_tranche_loss<T: Real>(
    default_probs: &[T],
    loss_units: &[usize],
    unit_notional: T,
    correlation: T,
    attachment: T,
    detachment: T,
) -> T {
    assert_eq!(default_probs.len(), loss_units.len());
    let max_loss: usize = loss_units.iter().sum();
    if max_loss == 0 {
        return T::zero();
    }
    let (nodes32, weights32) = normal_hermite_32();
    let high_rule = (correlation.val() > 0.8).then(high_corr_normal_rule);
    let (nodes, weights): (&[f64], &[f64]) = if let Some((ref n, ref w)) = high_rule {
        (n.as_slice(), w.as_slice())
    } else {
        (&nodes32[..], &weights32[..])
    };
    let q_invs: Vec<T> = default_probs.iter().map(|&q| inv_norm_cdf_ad(q)).collect();
    let rho = correlation
        .max_val(T::zero())
        .min_val(T::from_f64(1.0 - 1e-14));
    let sr = rho.sqrt();
    let den = (T::one() - rho).sqrt();
    let mut out = T::zero();
    let mut dist = vec![T::zero(); max_loss + 1];
    for (&z, &quad_w) in nodes.iter().zip(weights.iter()) {
        dist.fill(T::zero());
        dist[0] = T::one();
        let mut reached = 0usize;
        for (&q_inv, &w_loss) in q_invs.iter().zip(loss_units.iter()) {
            if w_loss == 0 {
                continue;
            }
            let pk = ((q_inv - sr * z) / den).norm_cdf();
            let one_minus_pk = T::one() - pk;
            for l in (0..=reached).rev() {
                let cur = dist[l];
                dist[l + w_loss] += cur * pk;
                dist[l] = cur * one_minus_pk;
            }
            reached += w_loss;
        }
        let width = (detachment - attachment).max_val(T::zero());
        let mut conditional = T::zero();
        for (l, &prob) in dist[..=reached].iter().enumerate() {
            let loss = unit_notional * l as f64;
            conditional += prob * (loss - attachment).max_val(T::zero()).min_val(width);
        }
        out += conditional * quad_w;
    }
    out.max_val(T::zero())
        .min_val((detachment - attachment).max_val(T::zero()))
}

/// Backwards-compatible short spelling used by the structured-credit API.
pub fn asb_recursive_tranche_loss<T: Real>(
    default_probs: &[T],
    loss_units: &[usize],
    unit_notional: T,
    correlation: T,
    attachment: T,
    detachment: T,
) -> T {
    asb_recursive_expected_tranche_loss(
        default_probs,
        loss_units,
        unit_notional,
        correlation,
        attachment,
        detachment,
    )
}

#[derive(Clone, Debug)]
pub struct SyntheticCdoTranche<T: Real> {
    pub attachment: T,
    pub detachment: T,
    pub running_spread: T,
    pub payment_times: Vec<T>,
    pub accruals: Vec<T>,
}

impl<T: Real> SyntheticCdoTranche<T> {
    pub fn price(
        &self,
        default_probs: &[T],
        loss_units: &[usize],
        unit_notional: T,
        correlation: T,
        discount_factors: &[T],
    ) -> (T, T, T) {
        price_cdo_tranche(
            self,
            default_probs,
            loss_units,
            unit_notional,
            correlation,
            discount_factors,
        )
    }
}

/// Price a synthetic tranche.  `discount_factors` and `expected_loss`'s
/// horizon are aligned with `payment_times`; midpoint default accrual is
/// represented by the loss increment between consecutive dates.
pub fn price_cdo_tranche<T: Real>(
    tranche: &SyntheticCdoTranche<T>,
    default_probs: &[T],
    loss_units: &[usize],
    unit_notional: T,
    correlation: T,
    discount_factors: &[T],
) -> (T, T, T) {
    let mut protection = T::zero();
    let mut premium = T::zero();
    let mut prev_loss = T::zero();
    let n_issuers = loss_units.len();
    let has_term_structure =
        n_issuers > 0 && default_probs.len() == n_issuers * tranche.payment_times.len();
    for (i, &t) in tranche.payment_times.iter().enumerate() {
        let q = if has_term_structure {
            default_probs[i * n_issuers..(i + 1) * n_issuers].to_vec()
        } else {
            default_probs.to_vec()
        };
        let loss = asb_recursive_expected_tranche_loss(
            &q,
            loss_units,
            unit_notional,
            correlation,
            tranche.attachment,
            tranche.detachment,
        );
        let df = discount_factors.get(i).copied().unwrap_or(T::one());
        let accrual = tranche.accruals.get(i).copied().unwrap_or_else(|| {
            if i == 0 {
                t
            } else {
                t - tranche.payment_times[i - 1]
            }
        });
        protection += df * (loss - prev_loss);
        premium += df
            * tranche.running_spread
            * accrual
            * ((tranche.detachment - tranche.attachment).max_val(T::zero())
                - (loss + prev_loss) * 0.5);
        prev_loss = loss;
    }
    let fair = if premium.val().abs() > 1e-15 {
        protection / premium * tranche.running_spread
    } else {
        T::zero()
    };
    (protection, premium, fair)
}

#[derive(Clone, Debug)]
pub struct KthToDefaultBasketSpec<T: Real> {
    pub k: usize,
    pub notional: T,
    pub running_spread: T,
    pub correlation: T,
    pub recovery_rates: Vec<T>,
    pub payment_times: Vec<f64>,
    pub accruals: Vec<f64>,
}
#[derive(Copy, Clone, Debug)]
pub struct KthToDefaultValuation<T: Real> {
    pub protection_leg_npv: T,
    pub premium_leg_npv: T,
    pub risky_annuity: T,
    pub fair_spread: T,
    pub npv_buyer: T,
}
pub fn kth_to_default_survival_probability<T: Real>(
    k: usize,
    default_probs: &[T],
    correlation: T,
) -> T {
    // The survival leg of a k-th-to-default contract is the probability that
    // fewer than k names have defaulted.  Under the one-factor Gaussian
    // copula, condition on the common factor and evaluate the resulting
    // Poisson-binomial distribution with a short dynamic program.  This keeps
    // the function differentiable in all probabilities and in correlation.
    if k == 0 || default_probs.is_empty() || k > default_probs.len() {
        return if k == 0 { T::zero() } else { T::one() };
    }
    let (nodes32, weights32) = normal_hermite_32();
    let high_rule = (correlation.val() > 0.8).then(high_corr_normal_rule);
    let (nodes, weights): (&[f64], &[f64]) = if let Some((ref n, ref w)) = high_rule {
        (n.as_slice(), w.as_slice())
    } else {
        (&nodes32[..], &weights32[..])
    };
    let rho = correlation
        .max_val(T::zero())
        .min_val(T::from_f64(1.0 - 1e-14));
    let sr = rho.sqrt();
    let den = (T::one() - rho).sqrt();
    let q_invs: Vec<T> = default_probs.iter().map(|&q| inv_norm_cdf_ad(q)).collect();
    let mut out = T::zero();
    // `dist[j]` stores the conditional probability of exactly j defaults;
    // only j < k is needed for the survival event.
    let mut dist = vec![T::zero(); k];
    for (&z, &weight) in nodes.iter().zip(weights.iter()) {
        dist.fill(T::zero());
        dist[0] = T::one();
        for q_inv in &q_invs {
            let p = ((*q_inv - sr * z) / den)
                .norm_cdf()
                .max_val(T::zero())
                .min_val(T::one());
            // Reverse iteration prevents an issuer from being counted twice.
            for j in (1..k).rev() {
                dist[j] = dist[j] * (T::one() - p) + dist[j - 1] * p;
            }
            dist[0] *= T::one() - p;
        }
        out += dist.iter().copied().fold(T::zero(), |a, x| a + x) * weight;
    }
    out.max_val(T::zero()).min_val(T::one())
}

pub fn price_kth_to_default_basket<T: Real>(
    s: &KthToDefaultBasketSpec<T>,
    issuer_curves: &[crate::SurvivalProbabilityCurve<T>],
    discount: &impl ql_curves::YieldTermStructure<T>,
) -> KthToDefaultValuation<T> {
    if issuer_curves.is_empty() || s.k == 0 || s.k > issuer_curves.len() {
        return KthToDefaultValuation {
            protection_leg_npv: T::zero(),
            premium_leg_npv: T::zero(),
            risky_annuity: T::zero(),
            fair_spread: T::zero(),
            npv_buyer: T::zero(),
        };
    }
    let mut protection = T::zero();
    let mut risky_annuity = T::zero();
    let mut prev_survival = T::one();
    let mut prev_ps = vec![T::zero(); issuer_curves.len()];
    for (i, &t) in s.payment_times.iter().enumerate() {
        if !t.is_finite() || t < 0.0 {
            continue;
        }
        let ps: Vec<T> = issuer_curves
            .iter()
            .map(|c| T::one() - c.survival_probability(t))
            .collect();
        let survival = kth_to_default_survival_probability(s.k, &ps, s.correlation);
        let df = discount.discount_factor(t);
        let mut w_lgd = T::zero();
        let mut w_sum = T::zero();
        for m in 0..ps.len() {
            let dp = (ps[m]
                - if m < prev_ps.len() {
                    prev_ps[m]
                } else {
                    T::zero()
                })
            .max_val(T::zero());
            let sub: Vec<T> = ps
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != m)
                .map(|(_, p)| *p)
                .collect();
            let cav = if s.k == 1 {
                T::one()
            } else {
                kth_to_default_survival_probability(s.k - 1, &sub, s.correlation)
            };
            let wm = dp * cav;
            w_sum += wm;
            w_lgd += wm * (T::one() - s.recovery_rates.get(m).copied().unwrap_or(T::from_f64(0.4)));
        }
        let event = (prev_survival - survival).max_val(T::zero());
        let eff_lgd = if w_sum.val() > 1e-14 {
            w_lgd / w_sum
        } else {
            T::from_f64(0.4)
        };
        protection += s.notional * df * eff_lgd * event;
        let accrual = T::from_f64(*s.accruals.get(i).unwrap_or(&1.0));
        // Premium is paid while the basket is alive.  Midpoint default
        // accrual can be enabled by passing the usual half-period survival;
        // the endpoint survival is the conservative standard convention.
        risky_annuity += s.notional * df * accrual * (prev_survival + survival) * 0.5;
        prev_survival = survival;
        prev_ps = ps;
    }
    let premium = risky_annuity * s.running_spread;
    let fair = if risky_annuity.val().abs() > 1e-14 {
        protection / risky_annuity
    } else {
        T::zero()
    };
    KthToDefaultValuation {
        protection_leg_npv: protection,
        premium_leg_npv: premium,
        risky_annuity,
        fair_spread: fair,
        npv_buyer: protection - premium,
    }
}

#[cfg(test)]
mod phase13_tests {
    use super::*;
    #[test]
    fn kth_finite() {
        assert!(
            kth_to_default_survival_probability(1, &[0.1, 0.2], 0.3)
                .val()
                .is_finite()
        );
    }
}
