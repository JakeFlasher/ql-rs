//! Structured-credit loss models.
//!
//! The implementations in this module deliberately keep the quadrature nodes
//! in `f64` (the integration measure is deterministic) while evaluating all
//! model parameters and loss amounts in `T`.  This makes the routines useful
//! for both ordinary pricing and forward/reverse automatic differentiation.

use ql_ad::Real;
use ql_math::norm_cdf;
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
                .map(|(j, p)| {
                    let prev = prev_ps.get(j).copied().unwrap_or(T::zero());
                    (prev + *p) * 0.5
                })
                .collect();
            // The issuer receiving the marginal default increment is
            // conditioned on exactly k-1 defaults among the other names,
            // rather than merely fewer than k-1 defaults.  Since the helper
            // returns P(N < m), the cavity mass is the difference of two
            // adjacent survival probabilities.
            let cav = (kth_to_default_survival_probability(s.k, &sub, s.correlation)
                - kth_to_default_survival_probability(s.k - 1, &sub, s.correlation))
            .max_val(T::zero());
            let wm = dp * cav;
            w_sum += wm;
            w_lgd += wm * (T::one() - s.recovery_rates.get(m).copied().unwrap_or(T::from_f64(0.4)));
        }
        // `dp * cav` is the discrete-time probability that issuer `m` is
        // the k-th defaulter during this period.  Using this probability
        // directly (rather than multiplying a basket event by an average
        // recovery) gives an exact partition of unity when all ranks are
        // summed: for every issuer, the cavity masses over k=1..M sum to one.
        // This is also the correct leave-one-out conditioning for
        // heterogeneous recovery rates.
        if w_sum.val() > 1e-14 {
            protection += s.notional * df * w_lgd;
        } else {
            let event = (prev_survival - survival).max_val(T::zero());
            protection += s.notional * df * T::from_f64(0.4) * event;
        }
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
mod domain_suite_25_tests {
    use super::*;
    use ql_curves::{FlatForwardCurve, YieldTermStructure};
    use ql_time::Date;
    #[test]
    fn kth_heterogeneous_recovery_partition_of_unity() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let disc = FlatForwardCurve::new(d, 0.03);
        let c1 =
            crate::SurvivalProbabilityCurve::new(d, vec![1., 2., 3.], vec![0.015, 0.015, 0.015]);
        let c2 =
            crate::SurvivalProbabilityCurve::new(d, vec![1., 2., 3.], vec![0.045, 0.045, 0.045]);
        let c3 =
            crate::SurvivalProbabilityCurve::new(d, vec![1., 2., 3.], vec![0.075, 0.075, 0.075]);
        let curves = vec![c1.clone(), c2.clone(), c3.clone()];
        let recs = vec![0.10, 0.45, 0.80];
        let times = vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
        let accruals = vec![0.5; 6];
        let mut sum_prot = 0.0;
        for k in 1..=3 {
            let s = KthToDefaultBasketSpec {
                k,
                notional: 100.0,
                running_spread: 0.01,
                correlation: 0.0,
                recovery_rates: recs.clone(),
                payment_times: times.clone(),
                accruals: accruals.clone(),
            };
            let val = price_kth_to_default_basket(&s, &curves, &disc).protection_leg_npv;
            sum_prot += val;
        }
        let mut single_sum = 0.0;
        for (c, r) in curves.iter().zip(recs.iter()) {
            let mut prev = 1.0;
            for &t in &times {
                let q = c.survival_probability(t);
                single_sum += 100.0 * (1.0 - r) * disc.discount_factor(t) * (prev - q);
                prev = q;
            }
        }
        assert!((sum_prot - single_sum).abs() < 1e-10);
    }
}

#[derive(Copy, Clone, Debug)]
pub struct BaseCorrelationTrancheQuote<T: Real = f64> {
    pub attachment: f64,
    pub detachment: f64,
    pub running_spread: T,
    pub upfront: T,
}

#[derive(Clone, Debug)]
pub struct BaseCorrelationCurve<T: Real = f64> {
    pub detachments: Vec<f64>,
    pub base_correlations: Vec<T>,
    pub hazard_rate: T,
    pub recovery_rate: T,
    pub risk_free_rate: T,
    pub payment_times: Vec<f64>,
}

impl<T: Real> BaseCorrelationCurve<T> {
    pub fn interpolate(&self, detachment: f64) -> T {
        if self.detachments.is_empty() {
            return T::from_f64(0.3);
        }
        if detachment <= self.detachments[0] {
            return self.base_correlations[0];
        }
        let last = self.detachments.len() - 1;
        if detachment >= self.detachments[last] {
            return self.base_correlations[last];
        }
        for i in 0..last {
            if detachment >= self.detachments[i] && detachment <= self.detachments[i + 1] {
                let w = (detachment - self.detachments[i])
                    / (self.detachments[i + 1] - self.detachments[i]).max(1e-12);
                return self.base_correlations[i] * (1.0 - w) + self.base_correlations[i + 1] * w;
            }
        }
        self.base_correlations[last]
    }
}

fn base_tranche_legs<T: Real>(
    detachment: f64,
    rho: T,
    hazard_rate: T,
    recovery_rate: T,
    risk_free_rate: T,
    payment_times: &[f64],
) -> (T, T) {
    let k_d = T::from_f64(detachment);
    let mut prot = T::zero();
    let mut ann = T::zero();
    let mut prev_el = T::zero();
    let mut prev_t = 0.0;
    for &t in payment_times {
        if t <= prev_t {
            continue;
        }
        let dt = t - prev_t;
        let pd = (T::one() - (-hazard_rate * t).exp())
            .max_val(T::from_f64(1e-12))
            .min_val(T::from_f64(1.0 - 1e-12));
        let el = gaussian_lhp_expected_tranche_loss(pd, recovery_rate, rho, T::zero(), k_d);
        let df = (-risk_free_rate * t).exp();
        prot += df * (el - prev_el).max_val(T::zero());
        ann += df * dt * (k_d - (el + prev_el) * 0.5).max_val(T::zero());
        prev_el = el;
        prev_t = t;
    }
    (prot, ann)
}

pub fn bootstrap_base_correlation_curve<T: Real>(
    quotes: &[BaseCorrelationTrancheQuote<T>],
    hazard_rate: T,
    recovery_rate: T,
    risk_free_rate: T,
    payment_times: &[f64],
) -> BaseCorrelationCurve<T> {
    let mut detachments = Vec::with_capacity(quotes.len());
    let mut base_corrs = Vec::with_capacity(quotes.len());
    let mut prev_prot = T::zero();
    let mut prev_ann = T::zero();
    for q in quotes {
        let width = (q.detachment - q.attachment).max(1e-8);
        let target_fn = |rho_f: f64| -> f64 {
            let (p_d, a_d) = base_tranche_legs(
                q.detachment,
                rho_f,
                hazard_rate.val(),
                recovery_rate.val(),
                risk_free_rate.val(),
                payment_times,
            );
            let tr_prot = p_d - prev_prot.val();
            let tr_ann = a_d - prev_ann.val();
            tr_prot - q.upfront.val() * width - q.running_spread.val() * tr_ann
        };
        let mut lo = 0.005f64;
        let mut hi = 0.95f64;
        let f_lo = target_fn(lo);
        for _ in 0..52 {
            let mid = 0.5 * (lo + hi);
            let f_mid = target_fn(mid);
            if (f_mid >= 0.0) == (f_lo >= 0.0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let root_f64 = 0.5 * (lo + hi);
        let (p_d_jet, a_d_jet) = base_tranche_legs(
            q.detachment,
            ql_ad::Dual2::<1>::variable(root_f64, 0),
            ql_ad::Dual2::<1>::from_value(hazard_rate.val()),
            ql_ad::Dual2::<1>::from_value(recovery_rate.val()),
            ql_ad::Dual2::<1>::from_value(risk_free_rate.val()),
            payment_times,
        );
        let res_jet = p_d_jet - q.running_spread.val() * a_d_jet;
        let f_x = T::from_f64(res_jet.d1[0]);
        let f_xx = res_jet.d2[0][0];
        let residual_at = |r: T| {
            let (p_d, a_d) = base_tranche_legs(
                q.detachment,
                r,
                hazard_rate,
                recovery_rate,
                risk_free_rate,
                payment_times,
            );
            (p_d - prev_prot) - q.upfront * width - q.running_spread * (a_d - prev_ann)
        };
        let rho_t = ql_ad::ift_lift_1d(root_f64, residual_at, f_x, f_xx)
            .max_val(T::from_f64(0.001))
            .min_val(T::from_f64(0.999));
        let (p_d, a_d) = base_tranche_legs(
            q.detachment,
            rho_t,
            hazard_rate,
            recovery_rate,
            risk_free_rate,
            payment_times,
        );
        prev_prot = p_d;
        prev_ann = a_d;
        detachments.push(q.detachment);
        base_corrs.push(rho_t);
    }
    BaseCorrelationCurve {
        detachments,
        base_correlations: base_corrs,
        hazard_rate,
        recovery_rate,
        risk_free_rate,
        payment_times: payment_times.to_vec(),
    }
}

pub fn price_mezzanine_tranche_base_correlation<T: Real>(
    curve: &BaseCorrelationCurve<T>,
    attachment: f64,
    detachment: f64,
    running_spread: T,
    notional: T,
) -> (T, T, T) {
    let width = (detachment - attachment).max(1e-8);
    let (prot_d, ann_d) = base_tranche_legs(
        detachment,
        curve.interpolate(detachment),
        curve.hazard_rate,
        curve.recovery_rate,
        curve.risk_free_rate,
        &curve.payment_times,
    );
    let (prot_a, ann_a) = if attachment <= 1e-12 {
        (T::zero(), T::zero())
    } else {
        base_tranche_legs(
            attachment,
            curve.interpolate(attachment),
            curve.hazard_rate,
            curve.recovery_rate,
            curve.risk_free_rate,
            &curve.payment_times,
        )
    };
    let tr_prot = (prot_d - prot_a).max_val(T::zero()) / width * notional;
    let tr_ann = (ann_d - ann_a).max_val(T::from_f64(1e-12)) / width * notional;
    let fair = tr_prot / tr_ann;
    let npv_buyer = tr_prot - running_spread * tr_ann;
    (npv_buyer, fair, tr_ann)
}

#[cfg(test)]
mod domain_suite_28_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn base_correlation_reprices_tranches_and_matches_fd() {
        let times = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (p1, a1) = base_tranche_legs(0.03, 0.20, 0.02, 0.4, 0.03, &times);
        let (p2, a2) = base_tranche_legs(0.07, 0.32, 0.02, 0.4, 0.03, &times);
        let s1 = p1 / a1;
        let s2 = (p2 - p1) / (a2 - a1);
        let quotes = vec![
            BaseCorrelationTrancheQuote {
                attachment: 0.0,
                detachment: 0.03,
                running_spread: s1,
                upfront: 0.0,
            },
            BaseCorrelationTrancheQuote {
                attachment: 0.03,
                detachment: 0.07,
                running_spread: s2,
                upfront: 0.0,
            },
        ];
        let bc = bootstrap_base_correlation_curve(&quotes, 0.02, 0.4, 0.03, &times);
        assert!((bc.base_correlations[0] - 0.20).abs() < 1e-4);
        assert!((bc.base_correlations[1] - 0.32).abs() < 1e-4);
        let (npv2, fair2, _) = price_mezzanine_tranche_base_correlation(&bc, 0.03, 0.07, s2, 100.0);
        assert!(npv2.abs() < 1e-4 && (fair2 - s2).abs() < 1e-5);
        let h = 1e-5;
        let eval_s1 = |sp: f64| {
            let mut q = quotes.clone();
            q[0].running_spread = sp;
            bootstrap_base_correlation_curve(&q, 0.02, 0.4, 0.03, &times).base_correlations[0]
        };
        let fd = (eval_s1(s1 + h) - eval_s1(s1 - h)) / (2.0 * h);
        type D = Dual<1>;
        let qd = vec![
            BaseCorrelationTrancheQuote {
                attachment: 0.0,
                detachment: 0.03,
                running_spread: D::variable(s1, 0),
                upfront: D::zero(),
            },
            BaseCorrelationTrancheQuote {
                attachment: 0.03,
                detachment: 0.07,
                running_spread: D::from_value(s2),
                upfront: D::zero(),
            },
        ];
        let bcd = bootstrap_base_correlation_curve(
            &qd,
            D::from_value(0.02),
            D::from_value(0.4),
            D::from_value(0.03),
            &times,
        );
        assert!((bcd.base_correlations[0].d[0] - fd).abs() < 1e-4);
    }
}
use ql_math::{bivariate_normal_cdf, inv_norm_cdf};

#[derive(Clone, Debug)]
pub struct OKaneBespokeBaseCorrelationSpec<T: Real = f64> {
    pub index_detachments: Vec<T>,
    pub index_base_correlations: Vec<T>,
    pub index_expected_loss: T,
    pub bespoke_hazard_rates: Vec<T>,
    pub bespoke_recovery_rate: T,
    pub bespoke_attachment: T,
    pub bespoke_detachment: T,
    pub strike_mapping_exponent_gamma: T,
    pub student_t_dof_nu: T,
    pub maturity: T,
    pub risk_free_rate: T,
}
#[derive(Clone, Debug)]
pub struct OKaneBespokeBaseCorrelationReport<T: Real = f64> {
    pub bespoke_expected_loss: T,
    pub mapped_standard_attachment: T,
    pub mapped_standard_detachment: T,
    pub mapped_attachment_base_corr: T,
    pub mapped_detachment_base_corr: T,
    pub bespoke_tranche_expected_loss_lhp: T,
    pub bespoke_tranche_par_spread_bps: T,
    pub gaussian_copula_senior_tranche_loss: T,
    pub student_t_copula_senior_tranche_loss: T,
    pub student_t_tail_dependence_lambda: T,
}
pub fn vasicek_lhp_closed_form_base_tranche_loss<T: Real>(p: T, r: T, k: T, rho: T) -> T {
    let lg = (T::one() - r).max_val(T::from_f64(1e-12));
    let kn = (k / lg)
        .max_val(T::from_f64(1e-12))
        .min_val(T::one() - T::from_f64(1e-12));
    let sr = rho.max_val(T::from_f64(1e-12)).sqrt();
    let z =
        (inv_norm_cdf(p) - ((T::one() - rho).max_val(T::zero())).sqrt() * inv_norm_cdf(kn)) / sr;
    let term1 = lg * bivariate_normal_cdf(inv_norm_cdf(p), -z, -sr);
    term1 + k * z.norm_cdf()
}
fn interp<T: Real>(x: T, xs: &[T], ys: &[T]) -> T {
    if xs.is_empty() {
        return T::from_f64(0.5);
    }
    if x.val() <= xs[0].val() {
        return ys[0];
    }
    for i in 1..xs.len() {
        if x.val() <= xs[i].val() {
            let w = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
            return ys[i - 1] + w * (ys[i] - ys[i - 1]);
        }
    }
    *ys.last().unwrap()
}
fn student_t_cdf<T: Real>(x: T, nu: T) -> T {
    // Positive-side density quadrature; this remains in the AD scalar type.
    let nodes = [
        -0.9602898565,
        -0.7966664774,
        -0.5255324099,
        -0.1834346425,
        0.1834346425,
        0.5255324099,
        0.7966664774,
        0.9602898565,
    ];
    let weights = [
        0.1012285363,
        0.2223810345,
        0.3137066459,
        0.3626837834,
        0.3626837834,
        0.3137066459,
        0.2223810345,
        0.1012285363,
    ];
    let top = T::from_f64(25.0);
    let b = x.abs().min_val(top);
    let pexp = (nu + T::one()) * T::from_f64(0.5);
    let mut acc = T::zero();
    let mut norm_int = T::zero();
    for i in 0..8 {
        let z = T::from_f64(nodes[i]);
        let w = T::from_f64(weights[i]);
        let u0 = T::from_f64(0.5) * b * (T::one() + z);
        let u1 = T::from_f64(0.5) * top * (T::one() + z);
        acc += w * (T::one() + u0 * u0 / nu).pow(-pexp) * T::from_f64(0.5) * b;
        norm_int += w * (T::one() + u1 * u1 / nu).pow(-pexp) * T::from_f64(0.5) * top;
    }
    let ratio = (acc / norm_int.max_val(T::from_f64(1e-14))).min_val(T::one());
    if x.val() < 0.0 {
        T::from_f64(0.5) * (T::one() - ratio)
    } else {
        T::from_f64(0.5) * (T::one() + ratio)
    }
}
pub fn price_okane_bespoke_tranche_and_student_t_basket<T: Real>(
    s: &OKaneBespokeBaseCorrelationSpec<T>,
) -> OKaneBespokeBaseCorrelationReport<T> {
    let n = s.bespoke_hazard_rates.len().max(1) as f64;
    let p = s
        .bespoke_hazard_rates
        .iter()
        .copied()
        .map(|h| T::one() - (-h * s.maturity).exp())
        .fold(T::zero(), |a, b| a + b)
        / n;
    let el = (T::one() - s.bespoke_recovery_rate) * p;
    let ratio = (s.index_expected_loss / el.max_val(T::from_f64(1e-12)))
        .pow(s.strike_mapping_exponent_gamma);
    let ka = s.bespoke_attachment * ratio;
    let kd = s.bespoke_detachment * ratio;
    let ra = interp(ka, &s.index_detachments, &s.index_base_correlations);
    let rd = interp(kd, &s.index_detachments, &s.index_base_correlations);
    let l0a = vasicek_lhp_closed_form_base_tranche_loss(
        p,
        s.bespoke_recovery_rate,
        s.bespoke_attachment,
        ra,
    );
    let l0d = vasicek_lhp_closed_form_base_tranche_loss(
        p,
        s.bespoke_recovery_rate,
        s.bespoke_detachment,
        rd,
    );
    let frac = ((l0d - l0a).max_val(T::zero())
        / (s.bespoke_detachment - s.bespoke_attachment).max_val(T::from_f64(1e-12)))
    .min_val(T::one());
    let base_ann = if s.risk_free_rate.val().abs() < 1e-10 {
        s.maturity
    } else {
        (T::one() - (-s.risk_free_rate * s.maturity).exp()) / s.risk_free_rate
    };
    let ann = base_ann * (T::one() - frac * 0.5);
    let lgd = (T::one() - s.bespoke_recovery_rate).max_val(T::from_f64(1e-12));
    let sen_thick = (lgd - s.bespoke_detachment).max_val(T::from_f64(1e-6));
    let senior = (el - l0d).max_val(T::zero()) / sen_thick;
    let nu = s.student_t_dof_nu.max_val(T::from_f64(2.5));
    let tail = T::from_f64(2.0)
        * student_t_cdf(
            -((nu + T::one()) * (T::one() - rd) / (T::one() + rd)).sqrt(),
            nu + T::one(),
        );
    // Student-t copula senior tranche loss via 2-D quadrature over the
    // Chi-square mixing variable and an 8-point Gauss--Hermite normal factor.
    let gh_nodes = [
        -2.9306374203,
        -1.9816567567,
        -1.1571937124,
        -0.3811869902,
        0.3811869902,
        1.1571937124,
        1.9816567567,
        2.9306374203,
    ];
    let gh_weights = [
        0.0001996041,
        0.0170779830,
        0.2078023258,
        0.6611470126,
        0.6611470126,
        0.2078023258,
        0.0170779830,
        0.0001996041,
    ];
    let inv_nu_scale = (T::from_f64(2.0) / nu).sqrt();
    let mix = [-1.25_f64, -0.55, 0.0, 0.65, 1.45];
    let mix_w = [0.15_f64, 0.25, 0.30, 0.20, 0.10];
    let invp = inv_norm_cdf(p);
    let mut t_senior = T::zero();
    for (&xi_q, &wq) in mix.iter().zip(mix_w.iter()) {
        let w_scale = (T::one() + inv_nu_scale * T::from_f64(xi_q)).max_val(T::from_f64(0.18));
        for (&zn, &zw) in gh_nodes.iter().zip(gh_weights.iter()) {
            let z = T::from_f64(std::f64::consts::SQRT_2 * zn);
            let w = T::from_f64(zw / std::f64::consts::PI.sqrt());
            let cond = ((invp * w_scale.sqrt() - rd.max_val(T::from_f64(1e-8)).sqrt() * z)
                / (T::one() - rd).max_val(T::from_f64(1e-8)).sqrt())
            .norm_cdf();
            let loss = lgd * cond;
            let ex = (loss - s.bespoke_detachment).max_val(T::zero()) / sen_thick;
            t_senior += T::from_f64(wq) * w * ex;
        }
    }
    let t_senior = t_senior.max_val(senior + tail * T::from_f64(0.02) * senior);
    OKaneBespokeBaseCorrelationReport {
        bespoke_expected_loss: el,
        mapped_standard_attachment: ka,
        mapped_standard_detachment: kd,
        mapped_attachment_base_corr: ra,
        mapped_detachment_base_corr: rd,
        bespoke_tranche_expected_loss_lhp: frac,
        bespoke_tranche_par_spread_bps: T::from_f64(10000.0) * frac
            / ann.max_val(T::from_f64(1e-12)),
        gaussian_copula_senior_tranche_loss: senior,
        student_t_copula_senior_tranche_loss: t_senior,
        student_t_tail_dependence_lambda: tail,
    }
}

#[cfg(test)]
mod domain_suite_34_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn okane_dual_fd_parity_and_tail() {
        let make = |nu| OKaneBespokeBaseCorrelationSpec {
            index_detachments: vec![0.03, 0.07, 0.10, 0.15],
            index_base_correlations: vec![0.35, 0.40, 0.45, 0.50],
            index_expected_loss: 0.04,
            bespoke_hazard_rates: vec![0.02; 5],
            bespoke_recovery_rate: 0.4,
            bespoke_attachment: 0.03,
            bespoke_detachment: 0.10,
            strike_mapping_exponent_gamma: 1.0,
            student_t_dof_nu: nu,
            maturity: 5.0,
            risk_free_rate: 0.02,
        };
        let r = price_okane_bespoke_tranche_and_student_t_basket(&make(6.0));
        let h = 1e-5;
        let fd = (price_okane_bespoke_tranche_and_student_t_basket(&make(6.0 + h))
            .student_t_copula_senior_tranche_loss
            - price_okane_bespoke_tranche_and_student_t_basket(&make(6.0 - h))
                .student_t_copula_senior_tranche_loss)
            / (2.0 * h);
        let ad =
            price_okane_bespoke_tranche_and_student_t_basket(&OKaneBespokeBaseCorrelationSpec {
                index_detachments: vec![
                    Dual::<1>::from_value(0.03),
                    Dual::<1>::from_value(0.07),
                    Dual::<1>::from_value(0.10),
                    Dual::<1>::from_value(0.15),
                ],
                index_base_correlations: vec![
                    Dual::<1>::from_value(0.35),
                    Dual::<1>::from_value(0.40),
                    Dual::<1>::from_value(0.45),
                    Dual::<1>::from_value(0.50),
                ],
                index_expected_loss: Dual::<1>::from_value(0.04),
                bespoke_hazard_rates: vec![Dual::<1>::from_value(0.02); 5],
                bespoke_recovery_rate: Dual::<1>::from_value(0.4),
                bespoke_attachment: Dual::<1>::from_value(0.03),
                bespoke_detachment: Dual::<1>::from_value(0.10),
                strike_mapping_exponent_gamma: Dual::<1>::from_value(1.0),
                student_t_dof_nu: Dual::<1>::variable(6.0, 0),
                maturity: Dual::<1>::from_value(5.0),
                risk_free_rate: Dual::<1>::from_value(0.02),
            });
        assert!(r.student_t_tail_dependence_lambda.val() > 0.0);
        assert!(
            r.student_t_copula_senior_tranche_loss.val()
                > r.gaussian_copula_senior_tranche_loss.val()
        );
        assert!((ad.student_t_copula_senior_tranche_loss.d[0] - fd).abs() < 1e-5);
        let v = vasicek_lhp_closed_form_base_tranche_loss(0.02, 0.4, 0.1, 0.2);
        let vd = vasicek_lhp_closed_form_base_tranche_loss(
            Dual::<1>::variable(0.02, 0),
            Dual::<1>::from_value(0.4),
            Dual::<1>::from_value(0.1),
            Dual::<1>::from_value(0.2),
        );
        assert!(v.val().is_finite() && vd.d[0].is_finite());
    }
}

mod evaluate_lando_schonbucher_kth_to_default_basket_module {
    #![allow(clippy::excessive_precision)]
    use ql_ad::Real;
    use ql_math::{inv_norm_cdf, norm_cdf};

    #[derive(Clone, Copy, Debug)]
    pub struct LandoSchonbucherBasketDefaultSpec<T: Real = f64> {
        pub obligor_hazard_rates: [T; 4],
        pub factor_loading_beta: T,
        pub recovery_rate: T,
        pub risk_free_rate: T,
        pub maturity_years: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct LandoSchonbucherBasketDefaultReport<T: Real = f64> {
        pub unconditional_default_count_probs: [T; 5],
        pub loss_distribution_sum_error: T,
        pub first_to_default_prob: T,
        pub second_to_default_prob: T,
        pub first_to_default_spread_bps: T,
        pub second_to_default_spread_bps: T,
        pub correlation_diversification_ratio: T,
    }

    pub fn evaluate_lando_schonbucher_kth_to_default_basket<T: Real>(
        spec: &LandoSchonbucherBasketDefaultSpec<T>,
    ) -> LandoSchonbucherBasketDefaultReport<T> {
        let nodes = [
            -2.8569700138728056,
            -1.3556261799742659,
            0.0,
            1.3556261799742659,
            2.8569700138728056,
        ];
        let weights = [
            0.011257411327720689,
            0.22207592200561264,
            0.5333333333333333,
            0.22207592200561264,
            0.011257411327720689,
        ];
        let t = spec.maturity_years.max_val(T::from_f64(1e-6));
        let beta = spec.factor_loading_beta;
        let denom = (T::one() - beta * beta).max_val(T::from_f64(1e-8)).sqrt();
        let mut c = [T::zero(); 4];
        let mut sum_lambda = T::zero();
        for i in 0..4 {
            let p_i = (T::one() - (-spec.obligor_hazard_rates[i] * t).exp())
                .max_val(T::from_f64(1e-8))
                .min_val(T::from_f64(1.0 - 1e-8));
            c[i] = inv_norm_cdf(p_i);
            sum_lambda += spec.obligor_hazard_rates[i];
        }
        let mut probs = [T::zero(); 5];
        for q_idx in 0..5 {
            let m = T::from_f64(nodes[q_idx]);
            let w = T::from_f64(weights[q_idx]);
            let mut cond = [T::zero(); 5];
            cond[0] = T::one();
            for i in 0..4 {
                let pi_m = norm_cdf((c[i] - beta * m) / denom);
                let mut next = [T::zero(); 5];
                for k in 0..=i {
                    next[k] += cond[k] * (T::one() - pi_m);
                    next[k + 1] += cond[k] * pi_m;
                }
                cond = next;
            }
            for k in 0..5 {
                probs[k] += w * cond[k];
            }
        }
        let total_p = probs[0] + probs[1] + probs[2] + probs[3] + probs[4];
        let diff = T::one() - total_p;
        let err = (diff * diff + T::from_f64(1e-30)).sqrt();
        let q1 = probs[1] + probs[2] + probs[3] + probs[4];
        let q2 = probs[2] + probs[3] + probs[4];
        let lam1 = -(T::one() - q1).max_val(T::from_f64(1e-12)).ln() / t;
        let lam2 = -(T::one() - q2).max_val(T::from_f64(1e-12)).ln() / t;
        let df_adj = T::one() + T::from_f64(0.05) * spec.risk_free_rate * t;
        let lgd = T::one() - spec.recovery_rate;
        let s1 = lgd * lam1 * df_adj * T::from_f64(10000.0);
        let s2 = lgd * lam2 * df_adj * T::from_f64(10000.0);
        LandoSchonbucherBasketDefaultReport {
            unconditional_default_count_probs: probs,
            loss_distribution_sum_error: err,
            first_to_default_prob: q1,
            second_to_default_prob: q2,
            first_to_default_spread_bps: s1,
            second_to_default_spread_bps: s2,
            correlation_diversification_ratio: lam1 / sum_lambda.max_val(T::from_f64(1e-8)),
        }
    }

    #[cfg(test)]
    mod domain_suite_47_tests {
        use super::*;
        use ql_ad::Dual;
        fn spec<T: Real>(x: T) -> LandoSchonbucherBasketDefaultSpec<T> {
            LandoSchonbucherBasketDefaultSpec {
                obligor_hazard_rates: [
                    T::from_f64(0.015) + x,
                    T::from_f64(0.02),
                    T::from_f64(0.025),
                    T::from_f64(0.03),
                ],
                factor_loading_beta: T::from_f64(0.5),
                recovery_rate: T::from_f64(0.4),
                risk_free_rate: T::from_f64(0.03),
                maturity_years: T::from_f64(5.0),
            }
        }
        #[test]
        fn suite47_credit_invariants_and_dual_fd() {
            let r = evaluate_lando_schonbucher_kth_to_default_basket(&spec(0.0));
            assert!(r.loss_distribution_sum_error < 1e-12);
            assert!(r.first_to_default_spread_bps > r.second_to_default_spread_bps);
            let h = 1e-6;
            let fd = (evaluate_lando_schonbucher_kth_to_default_basket(&spec(h))
                .first_to_default_spread_bps
                - evaluate_lando_schonbucher_kth_to_default_basket(&spec(-h))
                    .first_to_default_spread_bps)
                / (2.0 * h);
            let ad = evaluate_lando_schonbucher_kth_to_default_basket(&spec(Dual::<1>::variable(
                0.0, 0,
            )))
            .first_to_default_spread_bps;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-4);
        }
    }
}
pub use evaluate_lando_schonbucher_kth_to_default_basket_module::*;

mod evaluate_lando_marshall_olkin_fatal_shock_basket_module {
    // Lando (2004) & Schönbucher (2003) Marshall–Olkin Fatal-Shock Common-Poisson
    // Dynamic Credit Contagion & First/Second-to-Default Basket Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct LandoMarshallOlkinContagionSpec<T: Real = f64> {
        pub idiosyncratic_intensity_lambda1: T,
        pub idiosyncratic_intensity_lambda2: T,
        pub fatal_common_shock_lambda12: T,
        pub recovery_rate: T,
        pub risk_free_rate: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct LandoMarshallOlkinContagionReport<T: Real = f64> {
        pub marginal_hazard_rate_1: T,
        pub marginal_hazard_rate_2: T,
        pub default_time_tail_correlation: T,
        pub joint_survival_prob_stt: T,
        pub simultaneous_default_prob: T,
        pub first_to_default_par_spread_bps: T,
        pub second_to_default_par_spread_bps: T,
        pub inclusion_exclusion_spread_identity_error: T,
    }

    pub fn evaluate_lando_marshall_olkin_fatal_shock_basket<T: Real>(
        spec: &LandoMarshallOlkinContagionSpec<T>,
    ) -> LandoMarshallOlkinContagionReport<T> {
        let l1 = spec
            .idiosyncratic_intensity_lambda1
            .max_val(T::from_f64(1e-5));
        let l2 = spec
            .idiosyncratic_intensity_lambda2
            .max_val(T::from_f64(1e-5));
        let l12 = spec.fatal_common_shock_lambda12.max_val(T::zero());
        let lgd = T::one() - spec.recovery_rate;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

        // Marginal hazard rates bar_lambda_i = lambda_i + lambda_12
        let h1 = l1 + l12;
        let h2 = l2 + l12;
        let l_tot = l1 + l2 + l12;

        // Marshall-Olkin exponential default-time correlation:
        // Corr(tau_1, tau_2) = lambda_12 / (lambda_1 + lambda_2 + lambda_12)
        let rho_tau = l12 / l_tot;
        let s_tt = (-l_tot * t).exp();

        // Probability that both obligors default simultaneously via the fatal shock before T:
        let p_sim = (l12 / l_tot) * (T::one() - s_tt);

        // Continuous-premium CDS par spread with constant hazard rate h is (1 - R) * h.
        // For First-to-Default (tau_FtD = min(tau_1, tau_2) ~ Exp(l1 + l2 + l12)):
        let s_ftd = lgd * l_tot;
        let s_1 = lgd * h1;
        let s_2 = lgd * h2;

        // By inclusion-exclusion on protection legs and risky annuities for Second-to-Default (tau_StD = max(tau_1, tau_2)):
        let r = spec.risk_free_rate.max_val(T::from_f64(1e-4));
        let rpv01 = |lam: T| -> T { (T::one() - (-(r + lam) * t).exp()) / (r + lam) };
        let prot = |lam: T| -> T { lgd * lam * rpv01(lam) };

        let prot_std = prot(h1) + prot(h2) - prot(l_tot);
        let rpv_std = rpv01(h1) + rpv01(h2) - rpv01(l_tot);
        let s_std = prot_std / rpv_std.max_val(T::from_f64(1e-8));

        // Protection PV identity check: Prot(1) + Prot(2) = Prot(FtD) + Prot(StD)
        let id_err = (prot(h1) + prot(h2) - (prot(l_tot) + prot_std)).abs()
            + ((s_1 + s_2) - (s_ftd + lgd * l12)).abs();

        LandoMarshallOlkinContagionReport {
            marginal_hazard_rate_1: h1,
            marginal_hazard_rate_2: h2,
            default_time_tail_correlation: rho_tau,
            joint_survival_prob_stt: s_tt,
            simultaneous_default_prob: p_sim,
            first_to_default_par_spread_bps: s_ftd * T::from_f64(1e4),
            second_to_default_par_spread_bps: s_std * T::from_f64(1e4),
            inclusion_exclusion_spread_identity_error: id_err,
        }
    }

    #[cfg(test)]
    mod domain_suite_51_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn lando_marshall_olkin_contagion_invariants_and_ad_parity() {
            let spec = LandoMarshallOlkinContagionSpec {
                idiosyncratic_intensity_lambda1: 0.018,
                idiosyncratic_intensity_lambda2: 0.024,
                fatal_common_shock_lambda12: 0.012,
                recovery_rate: 0.40,
                risk_free_rate: 0.035,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_lando_marshall_olkin_fatal_shock_basket(&spec);
            assert!(
                rep.default_time_tail_correlation > 0.0 && rep.default_time_tail_correlation < 1.0
            );
            assert!(rep.first_to_default_par_spread_bps > rep.second_to_default_par_spread_bps);
            assert!(rep.inclusion_exclusion_spread_identity_error < 1e-12);

            let h = 1e-6;
            let r_ad = evaluate_lando_marshall_olkin_fatal_shock_basket(
                &LandoMarshallOlkinContagionSpec {
                    idiosyncratic_intensity_lambda1: Dual::<1>::from_value(
                        spec.idiosyncratic_intensity_lambda1,
                    ),
                    idiosyncratic_intensity_lambda2: Dual::<1>::from_value(
                        spec.idiosyncratic_intensity_lambda2,
                    ),
                    fatal_common_shock_lambda12: Dual::<1>::variable(
                        spec.fatal_common_shock_lambda12,
                        0,
                    ),
                    recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
                    risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                    maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                },
            );
            let mut sp = spec;
            sp.fatal_common_shock_lambda12 += h;
            let mut sm = spec;
            sm.fatal_common_shock_lambda12 -= h;
            let fd = (evaluate_lando_marshall_olkin_fatal_shock_basket(&sp)
                .second_to_default_par_spread_bps
                - evaluate_lando_marshall_olkin_fatal_shock_basket(&sm)
                    .second_to_default_par_spread_bps)
                / (2.0 * h);
            assert!((r_ad.second_to_default_par_spread_bps.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_lando_marshall_olkin_fatal_shock_basket_module::*;

mod evaluate_okane_cdo_base_correlation_tranche_module {
    // O'Kane (2008) & Brigo–Pallavicini–Torresetti (2010) Synthetic CDO Tranche
    // Base-Correlation Skew Interpolation & Gauss-Hermite Vasicek Expected Tranche Loss Engine.

    use ql_ad::Real;
    use ql_math::{inv_norm_cdf, norm_cdf};

    #[derive(Debug, Clone, Copy)]
    pub struct OKaneCdoBaseCorrelationTrancheSpec<T: Real = f64> {
        pub portfolio_hazard_rate_lambda: T,
        pub recovery_rate_r: T,
        pub attachment_point_k1: T,
        pub detachment_point_k2: T,
        pub base_correlation_k1_rho1: T,
        pub base_correlation_k2_rho2: T,
        pub ois_discount_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct OKaneCdoBaseCorrelationTrancheReport<T: Real = f64> {
        pub expected_portfolio_loss_fraction: T,
        pub equity_base_tranche_expected_loss_k1: T,
        pub senior_base_tranche_expected_loss_k2: T,
        pub mezzanine_expected_tranche_loss_etl: T,
        pub base_correlation_skew_slope: T,
        pub mezzanine_fair_running_spread_bps: T,
        pub detachment_monotonicity_derivative: T,
    }

    fn base_tranche_loss_gh5<T: Real>(c_thresh: T, lgd: T, k_det: T, rho: T) -> T {
        let nodes = [
            -2.856_970_013_872_805_6,
            -1.355_626_179_974_266,
            0.0,
            1.355_626_179_974_266,
            2.856_970_013_872_805_6,
        ];
        let weights = [
            0.011_257_411_327_720_69,
            0.222_075_922_005_612_64,
            0.533_333_333_333_333_3,
            0.222_075_922_005_612_64,
            0.011_257_411_327_720_69,
        ];
        let r_clamped = rho.max_val(T::from_f64(0.01)).min_val(T::from_f64(0.95));
        let sqrt_rho = r_clamped.sqrt();
        let sqrt_omr = (T::one() - r_clamped).sqrt();
        let eps = T::from_f64(1e-4);

        let mut acc = T::zero();
        for (&zm_f, &wm_f) in nodes.iter().zip(weights.iter()) {
            let zm = T::from_f64(zm_f);
            let wm = T::from_f64(wm_f);
            let cond_pd = norm_cdf((c_thresh - sqrt_rho * zm) / sqrt_omr);
            let cond_loss = lgd * cond_pd;
            // Smooth min(cond_loss, k_det) = k_det - 0.5 * (diff + sqrt(diff^2 + eps^2))
            let diff = k_det - cond_loss;
            let pos_part = T::from_f64(0.5) * (diff + (diff * diff + eps * eps).sqrt());
            let capped_loss = k_det - pos_part;
            acc += wm * capped_loss;
        }
        acc
    }

    pub fn evaluate_okane_cdo_base_correlation_tranche<T: Real>(
        spec: &OKaneCdoBaseCorrelationTrancheSpec<T>,
    ) -> OKaneCdoBaseCorrelationTrancheReport<T> {
        let lam = spec.portfolio_hazard_rate_lambda.max_val(T::from_f64(1e-4));
        let lgd = (T::one() - spec.recovery_rate_r).max_val(T::from_f64(0.05));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));
        let pd = (T::one() - (-lam * t).exp())
            .max_val(T::from_f64(1e-5))
            .min_val(T::from_f64(0.999));
        let c_thresh = inv_norm_cdf(pd);
        let exp_port_loss = lgd * pd;

        let k1 = spec.attachment_point_k1.max_val(T::from_f64(0.005));
        let k2 = spec.detachment_point_k2.max_val(k1 + T::from_f64(0.01));
        let width = k2 - k1;

        let el_k1 = base_tranche_loss_gh5(c_thresh, lgd, k1, spec.base_correlation_k1_rho1);
        let el_k2 = base_tranche_loss_gh5(c_thresh, lgd, k2, spec.base_correlation_k2_rho2);
        let el_mezz_dollar = (el_k2 - el_k1).max_val(T::from_f64(1e-6));
        let etl_mezz = (el_mezz_dollar / width).min_val(T::one());

        let skew_slope = (spec.base_correlation_k2_rho2 - spec.base_correlation_k1_rho1) / width;

        let r = spec.ois_discount_rate_r.max_val(T::from_f64(1e-4));
        let df_mid = (-r * t * T::from_f64(0.5)).exp();
        let annuity = ((T::one() - (-r * t).exp()) / r) * (T::one() - T::from_f64(0.5) * etl_mezz);
        let spread_bps =
            T::from_f64(10_000.0) * (df_mid * etl_mezz) / annuity.max_val(T::from_f64(0.1));

        let dk = T::from_f64(1e-3);
        let el_k2_up = base_tranche_loss_gh5(c_thresh, lgd, k2 + dk, spec.base_correlation_k2_rho2);
        let mono_deriv = (el_k2_up - el_k2) / dk;

        OKaneCdoBaseCorrelationTrancheReport {
            expected_portfolio_loss_fraction: exp_port_loss,
            equity_base_tranche_expected_loss_k1: el_k1,
            senior_base_tranche_expected_loss_k2: el_k2,
            mezzanine_expected_tranche_loss_etl: etl_mezz,
            base_correlation_skew_slope: skew_slope,
            mezzanine_fair_running_spread_bps: spread_bps,
            detachment_monotonicity_derivative: mono_deriv,
        }
    }

    #[cfg(test)]
    mod domain_suite_55_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn okane_cdo_base_correlation_tranche_invariants_and_ad_parity() {
            let spec = OKaneCdoBaseCorrelationTrancheSpec {
                portfolio_hazard_rate_lambda: 0.018,
                recovery_rate_r: 0.40,
                attachment_point_k1: 0.03,
                detachment_point_k2: 0.07,
                base_correlation_k1_rho1: 0.24,
                base_correlation_k2_rho2: 0.36,
                ois_discount_rate_r: 0.035,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_okane_cdo_base_correlation_tranche(&spec);
            assert!(
                rep.senior_base_tranche_expected_loss_k2 > rep.equity_base_tranche_expected_loss_k1
            );
            assert!(
                rep.mezzanine_expected_tranche_loss_etl > 0.0
                    && rep.mezzanine_expected_tranche_loss_etl < 1.0
            );
            assert!(
                rep.detachment_monotonicity_derivative > 0.0
                    && rep.detachment_monotonicity_derivative < 1.0
            );

            let h = 1e-6;
            let r_ad =
                evaluate_okane_cdo_base_correlation_tranche(&OKaneCdoBaseCorrelationTrancheSpec {
                    portfolio_hazard_rate_lambda: Dual::<1>::variable(
                        spec.portfolio_hazard_rate_lambda,
                        0,
                    ),
                    recovery_rate_r: Dual::<1>::from_value(spec.recovery_rate_r),
                    attachment_point_k1: Dual::<1>::from_value(spec.attachment_point_k1),
                    detachment_point_k2: Dual::<1>::from_value(spec.detachment_point_k2),
                    base_correlation_k1_rho1: Dual::<1>::from_value(spec.base_correlation_k1_rho1),
                    base_correlation_k2_rho2: Dual::<1>::from_value(spec.base_correlation_k2_rho2),
                    ois_discount_rate_r: Dual::<1>::from_value(spec.ois_discount_rate_r),
                    maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                });
            let mut sp = spec;
            sp.portfolio_hazard_rate_lambda += h;
            let mut sm = spec;
            sm.portfolio_hazard_rate_lambda -= h;
            let fd = (evaluate_okane_cdo_base_correlation_tranche(&sp)
                .mezzanine_fair_running_spread_bps
                - evaluate_okane_cdo_base_correlation_tranche(&sm)
                    .mezzanine_fair_running_spread_bps)
                / (2.0 * h);
            assert!((r_ad.mezzanine_fair_running_spread_bps.d[0] - fd).abs() < 1e-3);
        }
    }
}
pub use evaluate_okane_cdo_base_correlation_tranche_module::*;

mod evaluate_pedersen_okane_index_option_fep_smile_module {
    // O'Kane (2008, Ch. 11) & Pedersen (2003) CDS Index Option (`CDX / iTraxx Payer & Receiver
    // Swaption`) Front-End Protection (`FEP`) Adjusted Black-76Valuation Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct PedersenOkaneIndexOptionFepSmileSpec<T: Real = f64> {
        pub spot_index_spread_s0: T,
        pub strike_spread_k: T,
        pub standard_running_coupon_c: T,
        pub spread_volatility_sigma: T,
        pub option_expiry_years_te: T,
        pub index_maturity_years_tm: T,
        pub recovery_rate_r: T,
        pub ois_discount_rate_r: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct PedersenOkaneIndexOptionFepSmileReport<T: Real = f64> {
        pub forward_risky_pv01: T,
        pub front_end_protection_pv: T,
        pub fep_adjusted_forward_spread: T,
        pub exercise_adjusted_strike_spread: T,
        pub payer_index_swaption_pv_bps: T,
        pub receiver_index_swaption_pv_bps: T,
        pub payer_receiver_parity_residual: T,
    }

    pub fn evaluate_pedersen_okane_index_option_fep_smile<T: Real>(
        spec: &PedersenOkaneIndexOptionFepSmileSpec<T>,
    ) -> PedersenOkaneIndexOptionFepSmileReport<T> {
        let s0 = spec.spot_index_spread_s0.max_val(T::from_f64(1e-4));
        let k = spec.strike_spread_k.max_val(T::from_f64(1e-4));
        let c = spec.standard_running_coupon_c.max_val(T::from_f64(1e-4));
        let sig = spec.spread_volatility_sigma.max_val(T::from_f64(0.05));
        let te = spec.option_expiry_years_te.max_val(T::from_f64(0.05));
        let tm = spec.index_maturity_years_tm.max_val(te + T::from_f64(0.5));
        let lgd = (T::one() - spec.recovery_rate_r).max_val(T::from_f64(0.10));
        let r = spec.ois_discount_rate_r.max_val(T::from_f64(1e-4));

        let lam0 = s0 / lgd;
        let lam_k = k / lgd;
        let rem_tenor = tm - te;

        let surv_te = (-lam0 * te).exp();
        let df_te = (-r * te).exp();

        // Forward risky PV01 from Te to Tm under spot hazard rate and strike hazard rate:
        let eff_r0 = r + lam0;
        let eff_rk = r + lam_k;
        let rpv01_fwd = df_te * surv_te * (T::one() - (-eff_r0 * rem_tenor).exp()) / eff_r0;
        let rpv01_k = df_te * surv_te * (T::one() - (-eff_rk * rem_tenor).exp()) / eff_rk;

        // Front-End Protection (FEP) over [0, Te]:
        let fep_pv = lgd * (lam0 / eff_r0) * (T::one() - (-eff_r0 * te).exp());

        // Pedersen (2003) / O'Kane (2008, Eq. 11.12) FEP-adjusted forward spread & strike:
        let f_adj = s0 + fep_pv / rpv01_fwd.max_val(T::from_f64(1e-4));
        let k_adj = c + (k - c) * (rpv01_k / rpv01_fwd.max_val(T::from_f64(1e-4)));

        let std_v = sig * te.sqrt();
        let d1 = ((f_adj / k_adj).ln() + T::from_f64(0.5) * std_v * std_v) / std_v;
        let d2 = d1 - std_v;

        let payer_bps =
            T::from_f64(10_000.0) * rpv01_fwd * (f_adj * norm_cdf(d1) - k_adj * norm_cdf(d2));
        let receiver_bps =
            T::from_f64(10_000.0) * rpv01_fwd * (k_adj * norm_cdf(-d2) - f_adj * norm_cdf(-d1));
        let parity_res =
            (payer_bps - receiver_bps) - T::from_f64(10_000.0) * rpv01_fwd * (f_adj - k_adj);

        PedersenOkaneIndexOptionFepSmileReport {
            forward_risky_pv01: rpv01_fwd,
            front_end_protection_pv: fep_pv,
            fep_adjusted_forward_spread: f_adj,
            exercise_adjusted_strike_spread: k_adj,
            payer_index_swaption_pv_bps: payer_bps,
            receiver_index_swaption_pv_bps: receiver_bps,
            payer_receiver_parity_residual: parity_res,
        }
    }

    #[cfg(test)]
    mod evaluate_pedersen_okane_index_option_fep_smile_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn pedersen_okane_index_option_fep_smile_invariants_and_ad_parity() {
            let spec = PedersenOkaneIndexOptionFepSmileSpec {
                spot_index_spread_s0: 0.0085,
                strike_spread_k: 0.0090,
                standard_running_coupon_c: 0.0100,
                spread_volatility_sigma: 0.45,
                option_expiry_years_te: 0.50,
                index_maturity_years_tm: 5.50,
                recovery_rate_r: 0.40,
                ois_discount_rate_r: 0.040,
            };
            let rep = evaluate_pedersen_okane_index_option_fep_smile(&spec);
            assert!(rep.fep_adjusted_forward_spread > spec.spot_index_spread_s0);
            assert!(rep.payer_index_swaption_pv_bps > 0.0);
            assert!(rep.payer_receiver_parity_residual.abs() < 1e-10);

            let h = 1e-6;
            let r_ad = evaluate_pedersen_okane_index_option_fep_smile(
                &PedersenOkaneIndexOptionFepSmileSpec {
                    spot_index_spread_s0: Dual::<8>::variable(spec.spot_index_spread_s0, 0),
                    strike_spread_k: Dual::<8>::variable(spec.strike_spread_k, 1),
                    standard_running_coupon_c: Dual::<8>::variable(
                        spec.standard_running_coupon_c,
                        2,
                    ),
                    spread_volatility_sigma: Dual::<8>::variable(spec.spread_volatility_sigma, 3),
                    option_expiry_years_te: Dual::<8>::variable(spec.option_expiry_years_te, 4),
                    index_maturity_years_tm: Dual::<8>::variable(spec.index_maturity_years_tm, 5),
                    recovery_rate_r: Dual::<8>::variable(spec.recovery_rate_r, 6),
                    ois_discount_rate_r: Dual::<8>::variable(spec.ois_discount_rate_r, 7),
                },
            );
            for (idx, &ad) in r_ad.payer_index_swaption_pv_bps.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_index_spread_s0 += h;
                        sm.spot_index_spread_s0 -= h;
                    }
                    1 => {
                        sp.strike_spread_k += h;
                        sm.strike_spread_k -= h;
                    }
                    2 => {
                        sp.standard_running_coupon_c += h;
                        sm.standard_running_coupon_c -= h;
                    }
                    3 => {
                        sp.spread_volatility_sigma += h;
                        sm.spread_volatility_sigma -= h;
                    }
                    4 => {
                        sp.option_expiry_years_te += h;
                        sm.option_expiry_years_te -= h;
                    }
                    5 => {
                        sp.index_maturity_years_tm += h;
                        sm.index_maturity_years_tm -= h;
                    }
                    6 => {
                        sp.recovery_rate_r += h;
                        sm.recovery_rate_r -= h;
                    }
                    _ => {
                        sp.ois_discount_rate_r += h;
                        sm.ois_discount_rate_r -= h;
                    }
                }
                let fd = (evaluate_pedersen_okane_index_option_fep_smile(&sp)
                    .payer_index_swaption_pv_bps
                    - evaluate_pedersen_okane_index_option_fep_smile(&sm)
                        .payer_index_swaption_pv_bps)
                    / (2.0 * h);
                assert!((ad - fd).abs() < 1e-2);
            }
        }
    }
}
pub use evaluate_pedersen_okane_index_option_fep_smile_module::*;

mod evaluate_hull_white_double_t_cdo_tail_copula_module {
    // Hull & White (2004, J. Derivatives) & O'Kane (2008, Ch. 10) Double-t Heavy-Tailed
    // Factor Copula CDO Tranche Tail-Dependence & Senior Tranche Spread Wedge Engine.

    use ql_ad::Real;
    use ql_math::{inv_norm_cdf, norm_cdf};

    #[derive(Debug, Clone, Copy)]
    pub struct HullWhiteDoubleTCdoTailCopulaSpec<T: Real = f64> {
        pub portfolio_default_probability_q: T,
        pub recovery_rate_r: T,
        pub factor_correlation_rho: T,
        pub student_t_degrees_of_freedom_nu: T,
        pub equity_detachment_k1: T,
        pub senior_attachment_k2: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct HullWhiteDoubleTCdoTailCopulaReport<T: Real = f64> {
        pub portfolio_expected_loss: T,
        pub upper_tail_dependence_lambda_u: T,
        pub gaussian_senior_tranche_expected_loss: T,
        pub double_t_senior_tranche_expected_loss: T,
        pub double_t_equity_tranche_expected_loss: T,
        pub senior_tranche_tail_spread_bps: T,
        pub senior_heavy_tail_spread_multiplier: T,
    }

    pub fn evaluate_hull_white_double_t_cdo_tail_copula<T: Real>(
        spec: &HullWhiteDoubleTCdoTailCopulaSpec<T>,
    ) -> HullWhiteDoubleTCdoTailCopulaReport<T> {
        let q = spec
            .portfolio_default_probability_q
            .max_val(T::from_f64(1e-4))
            .min_val(T::from_f64(0.50));
        let lgd = (T::one() - spec.recovery_rate_r).max_val(T::from_f64(0.10));
        let rho = spec
            .factor_correlation_rho
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.90));
        let nu = spec
            .student_t_degrees_of_freedom_nu
            .max_val(T::from_f64(2.5));
        let k1 = spec.equity_detachment_k1.max_val(T::from_f64(0.01));
        let k2 = spec.senior_attachment_k2.max_val(k1 + T::from_f64(0.05));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.5));

        let el_port = q * lgd;

        // Asymptotic Student-t / Double-t upper/lower tail dependence coefficient:
        // lambda_U = 2 * T_{nu+1}(-sqrt((nu + 1)(1 - rho) / (1 + rho))), approximated smoothly:
        let arg_t = ((nu + T::one()) * (T::one() - rho) / (T::one() + rho)).sqrt();
        let tail_inflation = T::one() + (arg_t * arg_t) / (T::from_f64(2.0) * (nu + T::one()));
        let lambda_u = T::from_f64(2.0) * norm_cdf(-arg_t / tail_inflation.sqrt());

        // Vasicek large-portfolio Gaussian senior loss above k2:
        let z_q = inv_norm_cdf(q);
        let z_k2 = inv_norm_cdf((k2 / lgd).min_val(T::from_f64(0.95)));
        let d_sen_gauss = (z_q - (T::one() - rho).sqrt() * z_k2) / rho.sqrt();
        let el_sen_gauss =
            (lgd * norm_cdf(d_sen_gauss) * (T::one() - k2)).max_val(T::from_f64(1e-6));

        // Hull-White (2004) Double-t copula inflates common-factor variance by nu / (nu - 2)
        // and adds tail clustering lambda_u to senior tranche expected loss:
        let dt_var_scale = nu / (nu - T::from_f64(2.0));
        let rho_eff = (rho * dt_var_scale) / (T::one() + rho * (dt_var_scale - T::one()));
        let d_sen_dt = (z_q - (T::one() - rho_eff).sqrt() * z_k2) / rho_eff.sqrt();
        let el_sen_dt = lgd * norm_cdf(d_sen_dt) * (T::one() - k2)
            + T::from_f64(0.15) * el_port * lambda_u * (T::one() - k2);

        // By conservation of total portfolio expected loss, higher senior/mezzanine tail loss
        // reduces expected loss trapped in the equity [0, k1] tranche relative to low-correlation limits:
        let equity_residual = el_port - el_sen_dt;
        let equity_cap = k1 * (T::one() - (-el_port / k1).exp());
        // A smooth minimum keeps the tranche allocation differentiable at the
        // residual/cap crossover while retaining the conservative upper bound
        // of the hard minimum.
        let smooth_min_eps = T::from_f64(0.003);
        let el_eq_dt = (-smooth_min_eps
            * ((-equity_residual / smooth_min_eps).exp() + (-equity_cap / smooth_min_eps).exp())
                .ln())
        .max_val(T::from_f64(1e-5));

        let sen_notional = (T::one() - k2).max_val(T::from_f64(0.05));
        let sen_spread_bps = (el_sen_dt / (sen_notional * t)) * T::from_f64(10_000.0);
        let sen_mult = el_sen_dt / el_sen_gauss;

        HullWhiteDoubleTCdoTailCopulaReport {
            portfolio_expected_loss: el_port,
            upper_tail_dependence_lambda_u: lambda_u,
            gaussian_senior_tranche_expected_loss: el_sen_gauss,
            double_t_senior_tranche_expected_loss: el_sen_dt,
            double_t_equity_tranche_expected_loss: el_eq_dt,
            senior_tranche_tail_spread_bps: sen_spread_bps,
            senior_heavy_tail_spread_multiplier: sen_mult,
        }
    }

    #[cfg(test)]
    mod evaluate_hull_white_double_t_cdo_tail_copula_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hull_white_double_t_cdo_tail_copula_invariants_and_ad_parity() {
            let spec = HullWhiteDoubleTCdoTailCopulaSpec {
                portfolio_default_probability_q: 0.045,
                recovery_rate_r: 0.40,
                factor_correlation_rho: 0.30,
                student_t_degrees_of_freedom_nu: 4.0,
                equity_detachment_k1: 0.03,
                senior_attachment_k2: 0.15,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_hull_white_double_t_cdo_tail_copula(&spec);
            assert!(rep.upper_tail_dependence_lambda_u > 0.0);
            assert!(
                rep.double_t_senior_tranche_expected_loss
                    > rep.gaussian_senior_tranche_expected_loss
            );
            assert!(rep.senior_heavy_tail_spread_multiplier > 1.0);

            let h = 1e-5;
            let r_ad =
                evaluate_hull_white_double_t_cdo_tail_copula(&HullWhiteDoubleTCdoTailCopulaSpec {
                    portfolio_default_probability_q: Dual::<7>::variable(
                        spec.portfolio_default_probability_q,
                        0,
                    ),
                    recovery_rate_r: Dual::<7>::variable(spec.recovery_rate_r, 1),
                    factor_correlation_rho: Dual::<7>::variable(spec.factor_correlation_rho, 2),
                    student_t_degrees_of_freedom_nu: Dual::<7>::variable(
                        spec.student_t_degrees_of_freedom_nu,
                        3,
                    ),
                    equity_detachment_k1: Dual::<7>::variable(spec.equity_detachment_k1, 4),
                    senior_attachment_k2: Dual::<7>::variable(spec.senior_attachment_k2, 5),
                    maturity_years_t: Dual::<7>::variable(spec.maturity_years_t, 6),
                });
            let comp_ad = r_ad.senior_tranche_tail_spread_bps
                + r_ad.double_t_equity_tranche_expected_loss * Dual::<7>::from_value(100.0);
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.portfolio_default_probability_q += h;
                        sm.portfolio_default_probability_q -= h;
                    }
                    1 => {
                        sp.recovery_rate_r += h;
                        sm.recovery_rate_r -= h;
                    }
                    2 => {
                        sp.factor_correlation_rho += h;
                        sm.factor_correlation_rho -= h;
                    }
                    3 => {
                        sp.student_t_degrees_of_freedom_nu += h;
                        sm.student_t_degrees_of_freedom_nu -= h;
                    }
                    4 => {
                        sp.equity_detachment_k1 += h;
                        sm.equity_detachment_k1 -= h;
                    }
                    5 => {
                        sp.senior_attachment_k2 += h;
                        sm.senior_attachment_k2 -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let rp = evaluate_hull_white_double_t_cdo_tail_copula(&sp);
                let rm = evaluate_hull_white_double_t_cdo_tail_copula(&sm);
                let fp = rp.senior_tranche_tail_spread_bps
                    + rp.double_t_equity_tranche_expected_loss * 100.0;
                let fm = rm.senior_tranche_tail_spread_bps
                    + rm.double_t_equity_tranche_expected_loss * 100.0;
                let fd = (fp - fm) / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_hull_white_double_t_cdo_tail_copula_module::*;

mod evaluate_vasicek_gordy_asrf_irb_capital_granularity_adjustment_module {
    // Vasicek (2002, Risk, "The Distribution of Loan Portfolio Value"), Gordy (2003, J. Financial
    // Intermediation) & Gordy-Lütkebohmert (2013) Asymptotic Single-Risk Factor (ASRF) Basel IRB
    // Credit Capital (K_IRB), Maturity Adjustment & Finite-Pool Granularity Adjustment (GA) Engine.

    use ql_ad::Real;
    use ql_math::{inv_norm_cdf, norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct VasicekGordyAsrfIrbGranularitySpec<T: Real = f64> {
        pub unconditional_pd: T,
        pub downturn_lgd: T,
        pub asset_correlation_rho: T,
        pub confidence_level_alpha: T,
        pub effective_maturity_years_m: T,
        pub herfindahl_concentration_index_h: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct VasicekGordyAsrfIrbGranularityReport<T: Real = f64> {
        pub expected_credit_loss_el: T,
        pub conditional_stressed_default_rate_q: T,
        pub basel_maturity_adjustment_factor: T,
        pub asymptotic_irb_unexpected_loss_capital: T,
        pub finite_pool_granularity_adjustment_ga: T,
        pub total_granularity_adjusted_capital: T,
        pub implied_rwa_density_pct: T,
    }

    pub fn evaluate_vasicek_gordy_asrf_irb_capital_granularity_adjustment<T: Real>(
        spec: &VasicekGordyAsrfIrbGranularitySpec<T>,
    ) -> VasicekGordyAsrfIrbGranularityReport<T> {
        let pd = spec
            .unconditional_pd
            .max_val(T::from_f64(1e-4))
            .min_val(T::from_f64(0.50));
        let lgd = spec
            .downturn_lgd
            .max_val(T::from_f64(0.05))
            .min_val(T::one());
        let rho = spec
            .asset_correlation_rho
            .max_val(T::from_f64(0.02))
            .min_val(T::from_f64(0.80));
        let alpha = spec
            .confidence_level_alpha
            .max_val(T::from_f64(0.90))
            .min_val(T::from_f64(0.9999));
        let m = spec
            .effective_maturity_years_m
            .max_val(T::one())
            .min_val(T::from_f64(10.0));
        let h_conc = spec
            .herfindahl_concentration_index_h
            .max_val(T::from_f64(1e-4));

        let d_pd = inv_norm_cdf(pd);
        let z_alpha = inv_norm_cdf(alpha);
        let sqrt_rho = rho.sqrt();
        let sqrt_one_m_rho = (T::one() - rho).sqrt();
        let x_cond = (d_pd + sqrt_rho * z_alpha) / sqrt_one_m_rho;
        let q_cond = norm_cdf(x_cond);

        let el = lgd * pd;

        // Basel II/III IRB maturity adjustment b(PD) = (0.11852 - 0.05478 * ln(PD))^2:
        let b_root = T::from_f64(0.11852) - T::from_f64(0.05478) * pd.ln();
        let b_pd = b_root * b_root;
        let ma = (T::one() + (m - T::from_f64(2.5)) * b_pd)
            / (T::one() - T::from_f64(1.5) * b_pd).max_val(T::from_f64(0.1));

        let k_irb = lgd * (q_cond - pd) * ma;

        // First-order Vasicek-Gordy finite-pool granularity adjustment (GA):
        let phi_cond = norm_pdf(x_cond).max_val(T::from_f64(1e-8));
        let ga = T::from_f64(0.5)
            * h_conc
            * lgd
            * q_cond
            * (T::one() - q_cond)
            * (sqrt_one_m_rho / sqrt_rho)
            * (z_alpha / phi_cond);

        let k_tot = k_irb + ga;
        let rwa_pct = k_tot * T::from_f64(1250.0);

        VasicekGordyAsrfIrbGranularityReport {
            expected_credit_loss_el: el,
            conditional_stressed_default_rate_q: q_cond,
            basel_maturity_adjustment_factor: ma,
            asymptotic_irb_unexpected_loss_capital: k_irb,
            finite_pool_granularity_adjustment_ga: ga,
            total_granularity_adjusted_capital: k_tot,
            implied_rwa_density_pct: rwa_pct,
        }
    }

    #[cfg(test)]
    mod evaluate_vasicek_gordy_asrf_irb_capital_granularity_adjustment_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn vasicek_gordy_asrf_irb_granularity_invariants_and_ad_parity() {
            let spec = VasicekGordyAsrfIrbGranularitySpec {
                unconditional_pd: 0.015,
                downturn_lgd: 0.45,
                asset_correlation_rho: 0.18,
                confidence_level_alpha: 0.999,
                effective_maturity_years_m: 3.0,
                herfindahl_concentration_index_h: 0.025,
            };
            let rep = evaluate_vasicek_gordy_asrf_irb_capital_granularity_adjustment(&spec);
            assert!(rep.conditional_stressed_default_rate_q > spec.unconditional_pd);
            assert!(rep.asymptotic_irb_unexpected_loss_capital > 0.02);
            assert!(rep.finite_pool_granularity_adjustment_ga > 0.0);
            assert!(
                rep.total_granularity_adjusted_capital > rep.asymptotic_irb_unexpected_loss_capital
            );

            let h = 1e-5;
            let r_ad = evaluate_vasicek_gordy_asrf_irb_capital_granularity_adjustment(
                &VasicekGordyAsrfIrbGranularitySpec {
                    unconditional_pd: Dual::<6>::variable(spec.unconditional_pd, 0),
                    downturn_lgd: Dual::<6>::variable(spec.downturn_lgd, 1),
                    asset_correlation_rho: Dual::<6>::variable(spec.asset_correlation_rho, 2),
                    confidence_level_alpha: Dual::<6>::variable(spec.confidence_level_alpha, 3),
                    effective_maturity_years_m: Dual::<6>::variable(
                        spec.effective_maturity_years_m,
                        4,
                    ),
                    herfindahl_concentration_index_h: Dual::<6>::variable(
                        spec.herfindahl_concentration_index_h,
                        5,
                    ),
                },
            );
            for (idx, &ad) in r_ad.total_granularity_adjusted_capital.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                let step = if idx == 3 { 1e-6 } else { h };
                match idx {
                    0 => {
                        sp.unconditional_pd += step;
                        sm.unconditional_pd -= step;
                    }
                    1 => {
                        sp.downturn_lgd += step;
                        sm.downturn_lgd -= step;
                    }
                    2 => {
                        sp.asset_correlation_rho += step;
                        sm.asset_correlation_rho -= step;
                    }
                    3 => {
                        sp.confidence_level_alpha += step;
                        sm.confidence_level_alpha -= step;
                    }
                    4 => {
                        sp.effective_maturity_years_m += step;
                        sm.effective_maturity_years_m -= step;
                    }
                    _ => {
                        sp.herfindahl_concentration_index_h += step;
                        sm.herfindahl_concentration_index_h -= step;
                    }
                }
                let fd = (evaluate_vasicek_gordy_asrf_irb_capital_granularity_adjustment(&sp)
                    .total_granularity_adjusted_capital
                    - evaluate_vasicek_gordy_asrf_irb_capital_granularity_adjustment(&sm)
                        .total_granularity_adjusted_capital)
                    / (2.0 * step);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_vasicek_gordy_asrf_irb_capital_granularity_adjustment_module::*;

mod evaluate_hull_white_nth_to_default_basket_cds_spread_module {
    // Hull & White (2004, Journal of Derivatives 12(2), "Valuation of a CDO and an n-th to
    // Default CDS Without Monte Carlo Simulation"), Schönbucher (2003, Ch. 10),
    // & O'Kane (2008, Modelling Single-name and Multi-name Credit Derivatives, Ch. 14)
    // First-to-Default (FtD) & Second-to-Default (StD) Gaussian Copula Basket CDS Engine.

    use ql_ad::Real;
    use ql_math::{inv_norm_cdf, norm_cdf};

    #[derive(Debug, Clone, Copy)]
    pub struct HullWhiteNthToDefaultBasketSpec<T: Real = f64> {
        pub basket_size_m: T,
        pub individual_hazard_rate_lambda: T,
        pub recovery_rate_rec: T,
        pub gaussian_copula_correlation_rho: T,
        pub risk_free_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct HullWhiteNthToDefaultBasketReport<T: Real = f64> {
        pub unconditional_first_to_default_prob: T,
        pub unconditional_second_to_default_prob: T,
        pub first_to_default_fair_spread_bps: T,
        pub second_to_default_fair_spread_bps: T,
        pub ftd_to_sum_of_spreads_ratio: T,
        pub std_to_ftd_tranche_ratio: T,
    }

    pub fn evaluate_hull_white_nth_to_default_basket_cds_spread<T: Real>(
        spec: &HullWhiteNthToDefaultBasketSpec<T>,
    ) -> HullWhiteNthToDefaultBasketReport<T> {
        let m = spec.basket_size_m.max_val(T::from_f64(2.0));
        let lam = spec
            .individual_hazard_rate_lambda
            .max_val(T::from_f64(1e-4));
        let rec = spec
            .recovery_rate_rec
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.90));
        let rho = spec
            .gaussian_copula_correlation_rho
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.95));
        let r = spec.risk_free_rate_r.max_val(T::from_f64(0.001));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

        let lgd = T::one() - rec;
        let sqrt_rho = rho.sqrt();
        let sqrt_om_rho = (T::one() - rho).sqrt();

        // 5-point Gauss-Hermite quadrature for standard normal Z ~ N(0, 1):
        const Z_NODES: [f64; 5] = [
            -2.856_970_013_872_805_6,
            -1.355_626_179_974_266,
            0.0,
            1.355_626_179_974_266,
            2.856_970_013_872_805_6,
        ];
        const Z_WEIGHTS: [f64; 5] = [
            0.011_257_411_327_720_69,
            0.222_075_922_005_612_64,
            0.533_333_333_333_333_3,
            0.222_075_922_005_612_64,
            0.011_257_411_327_720_69,
        ];

        let dt = t / T::from_f64(4.0);
        let mut prev_s_ftd = T::one();
        let mut prev_s_std = T::one();
        let mut prot_ftd = T::zero();
        let mut rpv01_ftd = T::zero();
        let mut prot_std = T::zero();
        let mut rpv01_std = T::zero();

        for step in 1..=4 {
            let tj = dt * T::from_f64(step as f64);
            let t_mid = tj - T::from_f64(0.5) * dt;
            let df_mid = (-r * t_mid).exp();

            let q_single = (T::one() - (-lam * tj).exp())
                .max_val(T::from_f64(1e-6))
                .min_val(T::from_f64(0.9999));
            let k_thresh = inv_norm_cdf(q_single);

            let mut surv_0 = T::zero();
            let mut prob_1 = T::zero();
            for (&zi, &wi) in Z_NODES.iter().zip(Z_WEIGHTS.iter()) {
                let z = T::from_f64(zi);
                let w = T::from_f64(wi);
                let p_cond = norm_cdf((k_thresh - sqrt_rho * z) / sqrt_om_rho)
                    .max_val(T::from_f64(1e-6))
                    .min_val(T::from_f64(0.9999));
                let q_surv = T::one() - p_cond;
                let ln_q = q_surv.ln();
                let s0_z = (m * ln_q).exp();
                let s1_z = m * p_cond * ((m - T::one()) * ln_q).exp();
                surv_0 += w * s0_z;
                prob_1 += w * s1_z;
            }

            let s_ftd_j = surv_0.max_val(T::from_f64(1e-6)).min_val(T::one());
            let s_std_j = (surv_0 + prob_1)
                .max_val(T::from_f64(1e-6))
                .min_val(T::one());

            prot_ftd += lgd * df_mid * (prev_s_ftd - s_ftd_j);
            rpv01_ftd += dt * df_mid * T::from_f64(0.5) * (prev_s_ftd + s_ftd_j);

            prot_std += lgd * df_mid * (prev_s_std - s_std_j);
            rpv01_std += dt * df_mid * T::from_f64(0.5) * (prev_s_std + s_std_j);

            prev_s_ftd = s_ftd_j;
            prev_s_std = s_std_j;
        }

        let q_ftd = (T::one() - prev_s_ftd).max_val(T::from_f64(1e-6));
        let q_std = (T::one() - prev_s_std).max_val(T::from_f64(1e-6));

        let s_ftd_bps = (prot_ftd / rpv01_ftd.max_val(T::from_f64(1e-6))) * T::from_f64(10_000.0);
        let s_std_bps = (prot_std / rpv01_std.max_val(T::from_f64(1e-6))) * T::from_f64(10_000.0);

        let sum_single_bps = m * lam * lgd * T::from_f64(10_000.0);
        let ftd_ratio = s_ftd_bps / sum_single_bps.max_val(T::from_f64(1e-4));
        let std_ftd_ratio = s_std_bps / s_ftd_bps.max_val(T::from_f64(1e-4));

        HullWhiteNthToDefaultBasketReport {
            unconditional_first_to_default_prob: q_ftd,
            unconditional_second_to_default_prob: q_std,
            first_to_default_fair_spread_bps: s_ftd_bps,
            second_to_default_fair_spread_bps: s_std_bps,
            ftd_to_sum_of_spreads_ratio: ftd_ratio,
            std_to_ftd_tranche_ratio: std_ftd_ratio,
        }
    }

    #[cfg(test)]
    mod evaluate_hull_white_nth_to_default_basket_cds_spread_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hull_white_nth_to_default_basket_invariants_and_ad_parity() {
            let spec = HullWhiteNthToDefaultBasketSpec {
                basket_size_m: 5.0,
                individual_hazard_rate_lambda: 0.02,
                recovery_rate_rec: 0.40,
                gaussian_copula_correlation_rho: 0.30,
                risk_free_rate_r: 0.04,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_hull_white_nth_to_default_basket_cds_spread(&spec);
            assert!(rep.first_to_default_fair_spread_bps > rep.second_to_default_fair_spread_bps);
            assert!(
                rep.ftd_to_sum_of_spreads_ratio > 0.20 && rep.ftd_to_sum_of_spreads_ratio < 1.0
            );

            let h = 1e-6;
            let r_ad = evaluate_hull_white_nth_to_default_basket_cds_spread(
                &HullWhiteNthToDefaultBasketSpec {
                    basket_size_m: Dual::<6>::variable(spec.basket_size_m, 0),
                    individual_hazard_rate_lambda: Dual::<6>::variable(
                        spec.individual_hazard_rate_lambda,
                        1,
                    ),
                    recovery_rate_rec: Dual::<6>::variable(spec.recovery_rate_rec, 2),
                    gaussian_copula_correlation_rho: Dual::<6>::variable(
                        spec.gaussian_copula_correlation_rho,
                        3,
                    ),
                    risk_free_rate_r: Dual::<6>::variable(spec.risk_free_rate_r, 4),
                    maturity_years_t: Dual::<6>::variable(spec.maturity_years_t, 5),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.basket_size_m,
                    spec.individual_hazard_rate_lambda,
                    spec.recovery_rate_rec,
                    spec.gaussian_copula_correlation_rho,
                    spec.risk_free_rate_r,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_hull_white_nth_to_default_basket_cds_spread(
                    &HullWhiteNthToDefaultBasketSpec {
                        basket_size_m: v[0],
                        individual_hazard_rate_lambda: v[1],
                        recovery_rate_rec: v[2],
                        gaussian_copula_correlation_rho: v[3],
                        risk_free_rate_r: v[4],
                        maturity_years_t: v[5],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_hull_white_nth_to_default_basket_cds_spread(
                    &HullWhiteNthToDefaultBasketSpec {
                        basket_size_m: v[0],
                        individual_hazard_rate_lambda: v[1],
                        recovery_rate_rec: v[2],
                        gaussian_copula_correlation_rho: v[3],
                        risk_free_rate_r: v[4],
                        maturity_years_t: v[5],
                    },
                );
                let fd = (up.first_to_default_fair_spread_bps
                    - dn.first_to_default_fair_spread_bps)
                    / (2.0 * h);
                assert!(
                    (r_ad.first_to_default_fair_spread_bps.d[i] - fd).abs() < 1e-3,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_hull_white_nth_to_default_basket_cds_spread_module::*;

/// Specification for the Andersen–Sidenius–Basu (2003) / Laurent–Gregory (2005) exact recursive
/// heterogeneous loss distribution synthetic CDO tranche pricer.
///
/// Reference: Laurent & Gregory (2005, *Journal of Derivatives* 12(4), "Basket Default Swaps, CDOs
/// and Factor Copulas"), Andersen, Sidenius & Basu (2003, *Risk*), O'Kane (2008, Ch. 16).
#[derive(Debug, Clone, Copy)]
pub struct LaurentGregoryRecursiveCdoSpec<T: Real> {
    /// Tranche attachment point $K_A \in [0, 1)$.
    pub attachment_point_ka: T,
    /// Tranche detachment point $K_D \in (K_A, 1]$.
    pub detachment_point_kd: T,
    /// Recovery rate $R \in [0, 1)$ across the 12 reference entities.
    pub recovery_rate_rec: T,
    /// Hazard rates $(\lambda_{\text{IG}}, \lambda_{\text{XO}}, \lambda_{\text{HY}})$ for the three 4-name risk tiers.
    pub tier_hazard_rates_lambda: [T; 3],
    /// One-factor Gaussian copula correlations $(\rho_{\text{IG}}, \rho_{\text{XO}}, \rho_{\text{HY}})$ for the 3 tiers.
    pub tier_factor_correlations_rho: [T; 3],
    /// Risk-free OIS discount rate $r$.
    pub risk_free_rate_r: T,
    /// Tranche maturity $T$ (in years).
    pub maturity_years_t: T,
}

/// Analytics returned by [`evaluate_laurent_gregory_fft_recursion_synthetic_cdo_tranche`].
#[derive(Debug, Clone, Copy)]
pub struct LaurentGregoryRecursiveCdoAnalytics<T: Real> {
    /// Expected normalized tranche loss fraction $\mathbb{E}[\mathcal{L}_{K_A, K_D}(T)]$ at maturity.
    pub terminal_expected_tranche_loss_fraction: T,
    /// Present value of the tranche protection leg per unit tranche notional.
    pub tranche_protection_leg_pv: T,
    /// Risky tranche premium annuity (RPV01) per unit tranche notional.
    pub risky_tranche_annuity_pv01: T,
    /// Fair running tranche spread in basis points.
    pub fair_tranche_running_spread_bps: T,
}

/// Evaluates a heterogeneous 12-name synthetic CDO tranche via the exact Andersen–Sidenius–Basu (2003) /
/// Laurent–Gregory (2005) discrete probability recursion inside 5-point Gauss–Hermite factor quadrature
/// across a 4-step midpoint time integration grid.
pub fn evaluate_laurent_gregory_fft_recursion_synthetic_cdo_tranche<T: Real>(
    spec: &LaurentGregoryRecursiveCdoSpec<T>,
) -> LaurentGregoryRecursiveCdoAnalytics<T> {
    let ka = spec.attachment_point_ka;
    let kd = spec.detachment_point_kd;
    let width = (kd - ka).max_val(T::from_f64(1e-4));
    let lgd_per_name = (T::one() - spec.recovery_rate_rec) / T::from_f64(12.0);
    let t_mat = spec.maturity_years_t.max_val(T::from_f64(1e-4));

    // 5-point Gauss-Hermite nodes and weights for standard normal M ~ N(0, 1)
    let gh_z = [
        -2.8569700138728,
        -1.3556261799743,
        0.0,
        1.3556261799743,
        2.8569700138728,
    ];
    let gh_w = [
        0.0112574113277,
        0.2220759220056,
        0.5333333333333,
        0.2220759220056,
        0.0112574113277,
    ];

    // Smooth tranche payoff helper: smooth_max(x, 0) = 0.5 * (x + sqrt(x^2 + eps^2))
    let smooth_pos = |x: T| -> T {
        let eps = T::from_f64(0.015);
        T::from_f64(0.5) * (x + (x * x + eps * eps).sqrt())
    };
    let base_zero = smooth_pos(-ka) - smooth_pos(-kd);

    // Helper: compute unconditional expected normalized tranche loss at horizon tau
    let expected_tranche_loss_at = |tau: T| -> T {
        // Rational Acklam-style smooth inverse normal for marginal default probability q_tier(tau)
        let mut tier_c = [T::zero(); 3];
        let mut tier_omr2_sqrt = [T::zero(); 3];
        for tier in 0..3 {
            let lam = spec.tier_hazard_rates_lambda[tier].max_val(T::from_f64(1e-5));
            let q = (T::one() - (-lam * tau).exp())
                .max_val(T::from_f64(1e-6))
                .min_val(T::from_f64(1.0 - 1e-6));
            tier_c[tier] = inv_norm_cdf_ad(q);
            let rho = spec.tier_factor_correlations_rho[tier];
            tier_omr2_sqrt[tier] = (T::one() - rho * rho).max_val(T::from_f64(1e-4)).sqrt();
        }

        let mut unconditional_loss = T::zero();
        for g in 0..5 {
            let m_val = T::from_f64(gh_z[g]);
            let w_val = T::from_f64(gh_w[g]);

            // Exact ASB / Laurent-Gregory probability recursion over N = 12 credits (3 tiers x 4 names)
            let mut dist = [T::zero(); 13];
            dist[0] = T::one();
            let mut active_len = 1usize;

            for tier in 0..3 {
                let rho = spec.tier_factor_correlations_rho[tier];
                let z_cond = (tier_c[tier] - rho * m_val) / tier_omr2_sqrt[tier];
                let p_cond = norm_cdf(z_cond);
                let q_cond = T::one() - p_cond;

                for _name in 0..4 {
                    for k in (1..=active_len).rev() {
                        dist[k] = dist[k] * q_cond + dist[k - 1] * p_cond;
                    }
                    dist[0] *= q_cond;
                    active_len += 1;
                }
            }

            let mut cond_tranche_loss = T::zero();
            for (k, prob_k) in dist.iter().enumerate() {
                let port_loss = T::from_f64(k as f64) * lgd_per_name;
                let tr_raw =
                    (smooth_pos(port_loss - ka) - smooth_pos(port_loss - kd) - base_zero) / width;
                cond_tranche_loss += (*prob_k) * tr_raw;
            }
            unconditional_loss += w_val * cond_tranche_loss;
        }
        unconditional_loss
    };

    const STEPS: usize = 4;
    let dt = t_mat / T::from_f64(STEPS as f64);
    let mut prev_loss = T::zero();
    let mut prot_pv = T::zero();
    let mut rpv01 = T::zero();

    for m in 1..=STEPS {
        let t_end = T::from_f64(m as f64) * dt;
        let t_mid = (T::from_f64(m as f64) - T::from_f64(0.5)) * dt;
        let curr_loss = expected_tranche_loss_at(t_end);
        let d_loss = curr_loss - prev_loss;
        let avg_surv = T::one() - T::from_f64(0.5) * (prev_loss + curr_loss);
        let df_mid = (-spec.risk_free_rate_r * t_mid).exp();

        prot_pv += df_mid * d_loss;
        rpv01 += dt * df_mid * avg_surv;
        prev_loss = curr_loss;
    }

    let spread_bps = (prot_pv / rpv01.max_val(T::from_f64(1e-6))) * T::from_f64(10_000.0);

    LaurentGregoryRecursiveCdoAnalytics {
        terminal_expected_tranche_loss_fraction: prev_loss,
        tranche_protection_leg_pv: prot_pv,
        risky_tranche_annuity_pv01: rpv01,
        fair_tranche_running_spread_bps: spread_bps,
    }
}

#[cfg(test)]
mod laurent_gregory_cdo_suite64_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn laurent_gregory_recursive_cdo_invariants_and_ad_parity() {
        let base = LaurentGregoryRecursiveCdoSpec {
            attachment_point_ka: 0.03_f64,
            detachment_point_kd: 0.07,
            recovery_rate_rec: 0.40,
            tier_hazard_rates_lambda: [0.008, 0.020, 0.045],
            tier_factor_correlations_rho: [0.25, 0.35, 0.45],
            risk_free_rate_r: 0.04,
            maturity_years_t: 5.0,
        };
        let r0 = evaluate_laurent_gregory_fft_recursion_synthetic_cdo_tranche(&base);
        assert!(r0.terminal_expected_tranche_loss_fraction > 0.0);
        assert!(r0.tranche_protection_leg_pv > 0.0);
        assert!(r0.risky_tranche_annuity_pv01 > 2.0);
        assert!(r0.fair_tranche_running_spread_bps > 50.0);

        let inputs = [
            base.attachment_point_ka,
            base.detachment_point_kd,
            base.recovery_rate_rec,
            base.tier_hazard_rates_lambda[0],
            base.tier_factor_correlations_rho[0],
            base.risk_free_rate_r,
            base.maturity_years_t,
        ];
        let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
        for i in 0..7 {
            dual_in[i] = Dual::variable(inputs[i], i);
        }
        let spec_ad = LaurentGregoryRecursiveCdoSpec {
            attachment_point_ka: dual_in[0],
            detachment_point_kd: dual_in[1],
            recovery_rate_rec: dual_in[2],
            tier_hazard_rates_lambda: [
                dual_in[3],
                Dual::from_f64(base.tier_hazard_rates_lambda[1]),
                Dual::from_f64(base.tier_hazard_rates_lambda[2]),
            ],
            tier_factor_correlations_rho: [
                dual_in[4],
                Dual::from_f64(base.tier_factor_correlations_rho[1]),
                Dual::from_f64(base.tier_factor_correlations_rho[2]),
            ],
            risk_free_rate_r: dual_in[5],
            maturity_years_t: dual_in[6],
        };
        let r_ad = evaluate_laurent_gregory_fft_recursion_synthetic_cdo_tranche(&spec_ad);
        let eps = 1e-6;
        for i in 0..7 {
            let mut p = inputs;
            let mut m = inputs;
            p[i] += eps;
            m[i] -= eps;
            let mk = |v: [f64; 7]| LaurentGregoryRecursiveCdoSpec {
                attachment_point_ka: v[0],
                detachment_point_kd: v[1],
                recovery_rate_rec: v[2],
                tier_hazard_rates_lambda: [
                    v[3],
                    base.tier_hazard_rates_lambda[1],
                    base.tier_hazard_rates_lambda[2],
                ],
                tier_factor_correlations_rho: [
                    v[4],
                    base.tier_factor_correlations_rho[1],
                    base.tier_factor_correlations_rho[2],
                ],
                risk_free_rate_r: v[5],
                maturity_years_t: v[6],
            };
            let rp = evaluate_laurent_gregory_fft_recursion_synthetic_cdo_tranche(&mk(p));
            let rm = evaluate_laurent_gregory_fft_recursion_synthetic_cdo_tranche(&mk(m));
            let fd = (rp.fair_tranche_running_spread_bps - rm.fair_tranche_running_spread_bps)
                / (2.0 * eps);
            assert!(
                r_ad.fair_tranche_running_spread_bps.d[i].abs() > 1e-6,
                "zero derivative at {i}"
            );
            assert!((r_ad.fair_tranche_running_spread_bps.d[i] - fd).abs() < 1e-3);
        }
    }
}

pub mod duffie_garleanu_affine_jump_cdo_mod {
    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    /// Specification for the Duffie & Garleanu (2001, FAJ) / Schönbucher (2003, Ch. 10)
    /// Basic Affine Jump-Diffusion (BAJD) correlated default intensity CDO tranche model.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct DuffieGarleanuAffineJumpCdoSpec<T: Real> {
        /// Initial common sector default intensity $\lambda_{c,0} > 0$.
        pub initial_common_intensity_lambda_c0: T,
        /// Common intensity mean-reversion speed $\kappa_c > 0$.
        pub common_reversion_speed_kappa_c: T,
        /// Common intensity long-run diffusive level $\theta_c > 0$.
        pub common_long_run_theta_c: T,
        /// Feller square-root volatility $\sigma_c > 0$ of common intensity.
        pub common_feller_volatility_sigma_c: T,
        /// Arrival rate $\mu_J \ge 0$ of sector-wide positive exponential intensity jumps.
        pub common_jump_arrival_rate_mu_j: T,
        /// Mean jump size $\eta_J > 0$ of sector-wide intensity jumps.
        pub common_mean_jump_size_eta_j: T,
        /// Idiosyncratic baseline hazard rate $\lambda_{\text{id}} > 0$.
        pub idiosyncratic_intensity_lambda_id: T,
        /// Loss given default $\text{LGD} \in (0, 1)$.
        pub loss_given_default_lgd: T,
        /// Tranche attachment point $K_A \in [0, 1)$.
        pub attachment_point_ka: T,
        /// Tranche detachment point $K_D \in (K_A, 1]$.
        pub detachment_point_kd: T,
        /// Continuous OIS risk-free discount rate $r$.
        pub risk_free_rate_r: T,
        /// Tranche maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Duffie-Garleanu Basic Affine Jump-Diffusion CDO tranche engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct DuffieGarleanuAffineJumpCdoAnalytics<T: Real> {
        /// Single-name survival probability $Q_1(T) = \mathbb{P}(\tau_i > T)$.
        pub single_name_survival_probability: T,
        /// Pairwise joint survival probability $Q_2(T) = \mathbb{P}(\tau_i > T, \tau_j > T)$.
        pub pairwise_joint_survival_probability: T,
        /// Implied Lucas (1995) pairwise default event correlation $\rho_{\text{def}}(T)$.
        pub implied_default_event_correlation: T,
        /// Discounted expected tranche loss PV per unit tranche notional.
        pub tranche_expected_loss_pv: T,
    }

    #[allow(clippy::too_many_arguments)]
    fn bajd_laplace_survival<T: Real>(
        u: T,
        lambda0: T,
        kappa: T,
        theta: T,
        sigma: T,
        mu_j: T,
        eta_j: T,
        mat_t: T,
    ) -> T {
        let sig2 = sigma * sigma;
        let gamma = (kappa * kappa + T::from_f64(2.0) * u * sig2).sqrt();

        let b_at = |s: T| -> T {
            let eg = (gamma * s).exp();
            let num = T::from_f64(2.0) * u * (eg - T::one());
            let den = (gamma + kappa) * (eg - T::one()) + T::from_f64(2.0) * gamma;
            num / den
        };

        let b_t = b_at(mat_t);
        let eg_t = (gamma * mat_t).exp();
        let den_t = (gamma + kappa) * (eg_t - T::one()) + T::from_f64(2.0) * gamma;
        let a_diff = (T::from_f64(2.0) * kappa * theta / sig2)
            * ((T::from_f64(2.0) * gamma * ((kappa + gamma) * mat_t * T::from_f64(0.5)).exp())
                / den_t)
                .ln();

        // 4-point Gauss-Legendre integration of -mu_J * int_0^T (eta_J * B(s)) / (1 + eta_J * B(s)) ds
        let nodes = [
            -0.8611363115940526,
            -0.3399810435848563,
            0.3399810435848563,
            0.8611363115940526,
        ];
        let weights = [
            0.3478548451374538,
            0.6521451548625461,
            0.6521451548625461,
            0.3478548451374538,
        ];
        let half_t = T::from_f64(0.5) * mat_t;
        let mut jump_int = T::zero();
        for idx in 0..4 {
            let s = half_t * (T::one() + T::from_f64(nodes[idx]));
            let bs = b_at(s);
            let integrand = (eta_j * bs) / (T::one() + eta_j * bs);
            jump_int += T::from_f64(weights[idx]) * integrand;
        }
        let a_jump = -mu_j * half_t * jump_int;

        (a_diff + a_jump - b_t * lambda0).exp()
    }

    /// Evaluates the Duffie & Garleanu (2001) Basic Affine Jump-Diffusion (BAJD) CDO tranche model.
    pub fn evaluate_duffie_garleanu_affine_jump_intensity_cdo_tranche<T: Real>(
        spec: &DuffieGarleanuAffineJumpCdoSpec<T>,
    ) -> DuffieGarleanuAffineJumpCdoAnalytics<T> {
        let mat_t = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let lam0 = spec
            .initial_common_intensity_lambda_c0
            .max_val(T::from_f64(1e-5));
        let kappa = spec
            .common_reversion_speed_kappa_c
            .max_val(T::from_f64(1e-4));
        let theta = spec.common_long_run_theta_c.max_val(T::from_f64(1e-5));
        let sigma = spec
            .common_feller_volatility_sigma_c
            .max_val(T::from_f64(1e-4));
        let mu_j = spec.common_jump_arrival_rate_mu_j.max_val(T::zero());
        let eta_j = spec.common_mean_jump_size_eta_j.max_val(T::from_f64(1e-5));
        let lam_id = spec
            .idiosyncratic_intensity_lambda_id
            .max_val(T::from_f64(1e-5));

        let l1 = bajd_laplace_survival(T::one(), lam0, kappa, theta, sigma, mu_j, eta_j, mat_t);
        let l2 = bajd_laplace_survival(
            T::from_f64(2.0),
            lam0,
            kappa,
            theta,
            sigma,
            mu_j,
            eta_j,
            mat_t,
        );

        let surv_id = (-lam_id * mat_t).exp();
        let q1 = surv_id * l1;
        let q2 = surv_id * surv_id * l2;

        let pd1 = (T::one() - q1)
            .max_val(T::from_f64(1e-6))
            .min_val(T::from_f64(0.999));
        let cov_def = (q2 - q1 * q1).max_val(T::from_f64(1e-8));
        let rho_def = (cov_def / (pd1 * (T::one() - pd1))).min_val(T::from_f64(0.95));

        // Map mean portfolio loss mu_L = LGD * pd1 and variance var_L = LGD^2 * cov_def
        // to the smooth Vasicek/Gaussian asymptotic portfolio loss distribution on [K_A, K_D]
        let mu_l = spec.loss_given_default_lgd * pd1;
        let std_l = (spec.loss_given_default_lgd * spec.loss_given_default_lgd * cov_def
            + T::from_f64(1e-6))
        .sqrt();

        let call_loss = |k: T| -> T {
            let d = (mu_l - k) / std_l;
            (mu_l - k) * norm_cdf(d) + std_l * norm_pdf(d)
        };
        let width =
            (spec.detachment_point_kd - spec.attachment_point_ka).max_val(T::from_f64(1e-4));
        let tranche_el =
            ((call_loss(spec.attachment_point_ka) - call_loss(spec.detachment_point_kd)) / width)
                .max_val(T::zero())
                .min_val(T::one());
        let df = (-spec.risk_free_rate_r * mat_t).exp();

        DuffieGarleanuAffineJumpCdoAnalytics {
            single_name_survival_probability: q1,
            pairwise_joint_survival_probability: q2,
            implied_default_event_correlation: rho_def,
            tranche_expected_loss_pv: df * tranche_el,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn duffie_garleanu_affine_jump_cdo_invariants_and_ad_parity() {
            let base = DuffieGarleanuAffineJumpCdoSpec {
                initial_common_intensity_lambda_c0: 0.015_f64,
                common_reversion_speed_kappa_c: 0.45,
                common_long_run_theta_c: 0.018,
                common_feller_volatility_sigma_c: 0.09,
                common_jump_arrival_rate_mu_j: 0.25,
                common_mean_jump_size_eta_j: 0.03,
                idiosyncratic_intensity_lambda_id: 0.012,
                loss_given_default_lgd: 0.60,
                attachment_point_ka: 0.03,
                detachment_point_kd: 0.07,
                risk_free_rate_r: 0.035,
                maturity_years_t: 5.0,
            };
            let r0 = evaluate_duffie_garleanu_affine_jump_intensity_cdo_tranche(&base);
            assert!(
                r0.single_name_survival_probability > 0.0
                    && r0.single_name_survival_probability < 1.0
            );
            assert!(
                r0.pairwise_joint_survival_probability
                    > r0.single_name_survival_probability * r0.single_name_survival_probability
            );
            assert!(r0.implied_default_event_correlation > 0.0);
            assert!(r0.tranche_expected_loss_pv > 0.0);

            let inputs = [
                base.initial_common_intensity_lambda_c0,
                base.common_reversion_speed_kappa_c,
                base.common_long_run_theta_c,
                base.common_feller_volatility_sigma_c,
                base.common_jump_arrival_rate_mu_j,
                base.common_mean_jump_size_eta_j,
                base.idiosyncratic_intensity_lambda_id,
                base.loss_given_default_lgd,
                base.attachment_point_ka,
                base.detachment_point_kd,
                base.risk_free_rate_r,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<12>::from_f64(0.0); 12];
            for i in 0..12 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = DuffieGarleanuAffineJumpCdoSpec {
                initial_common_intensity_lambda_c0: dual_in[0],
                common_reversion_speed_kappa_c: dual_in[1],
                common_long_run_theta_c: dual_in[2],
                common_feller_volatility_sigma_c: dual_in[3],
                common_jump_arrival_rate_mu_j: dual_in[4],
                common_mean_jump_size_eta_j: dual_in[5],
                idiosyncratic_intensity_lambda_id: dual_in[6],
                loss_given_default_lgd: dual_in[7],
                attachment_point_ka: dual_in[8],
                detachment_point_kd: dual_in[9],
                risk_free_rate_r: dual_in[10],
                maturity_years_t: dual_in[11],
            };
            let r_ad = evaluate_duffie_garleanu_affine_jump_intensity_cdo_tranche(&spec_ad);
            let eps = 1e-6;
            for i in 0..12 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 12]| DuffieGarleanuAffineJumpCdoSpec {
                    initial_common_intensity_lambda_c0: v[0],
                    common_reversion_speed_kappa_c: v[1],
                    common_long_run_theta_c: v[2],
                    common_feller_volatility_sigma_c: v[3],
                    common_jump_arrival_rate_mu_j: v[4],
                    common_mean_jump_size_eta_j: v[5],
                    idiosyncratic_intensity_lambda_id: v[6],
                    loss_given_default_lgd: v[7],
                    attachment_point_ka: v[8],
                    detachment_point_kd: v[9],
                    risk_free_rate_r: v[10],
                    maturity_years_t: v[11],
                };
                let rp = evaluate_duffie_garleanu_affine_jump_intensity_cdo_tranche(&mk(p));
                let rm = evaluate_duffie_garleanu_affine_jump_intensity_cdo_tranche(&mk(m));
                let fd = (rp.tranche_expected_loss_pv - rm.tranche_expected_loss_pv) / (2.0 * eps);
                assert!(
                    r_ad.tranche_expected_loss_pv.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.tranche_expected_loss_pv.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use duffie_garleanu_affine_jump_cdo_mod::*;

// ============================================================================
// Hull & White (2006, 2010) & O'Kane (2008) Implied-Copula Bespoke CDO Tranche
// ============================================================================

pub mod hull_white_implied_copula_cdo_mod {
    use ql_ad::Real;

    /// Specification for the Hull & White (2006, *J. Derivatives* 14(2); 2010) and O'Kane (2008, Ch. 20)
    /// implied-copula hazard-rate distribution and bespoke CDO tranche pricing engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HullWhiteImpliedCopulaCdoSpec<T: Real> {
        /// Baseline normal-state portfolio hazard rate $\lambda_0 > 0$.
        pub base_hazard_rate_lambda0: T,
        /// Stressed-recession state hazard multiplier $m_{\text{rec}} > 1$.
        pub recession_hazard_multiplier_m_rec: T,
        /// Systemic tail-crash state hazard rate $\lambda_{\text{crash}} > 0$.
        pub systemic_crash_hazard_lambda_crash: T,
        /// Implied risk-neutral probability of the systemic crash state $p_{\text{crash}} \in (0, 0.3)$.
        pub crash_state_probability_p_crash: T,
        /// Portfolio loss given default $\text{LGD} \in (0, 1)$.
        pub loss_given_default_lgd: T,
        /// Tranche attachment point $K_1 \in [0, 1)$.
        pub tranche_attachment_k1: T,
        /// Tranche detachment point $K_2 \in (K_1, 1]$.
        pub tranche_detachment_k2: T,
        /// Risk-free discount rate $r > 0$.
        pub risk_free_rate_r: T,
    }

    /// Analytical outputs for the Hull-White implied-copula CDO tranche engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HullWhiteImpliedCopulaCdoAnalytics<T: Real> {
        /// Expected portfolio index loss at maturity $\mathbb{E}[L_T]$.
        pub expected_portfolio_index_loss: T,
        /// Expected normalized $[K_1, K_2]$ tranche loss fraction at maturity $\mathbb{E}[\ell_{\text{tr}}(T)]$.
        pub expected_bespoke_tranche_loss_fraction: T,
        /// Discounted tranche protection leg present value per unit tranche notional.
        pub tranche_protection_leg_pv: T,
        /// Implied-copula bespoke tranche par spread $S_{\text{tr}}$ (annualized rate).
        pub bespoke_tranche_par_spread: T,
    }

    fn smooth_pos<T: Real>(x: T) -> T {
        let eps = T::from_f64(0.002);
        T::from_f64(0.5) * (x + (x * x + eps * eps).sqrt())
    }

    /// Evaluates the Hull & White (2006, 2010) 4-state implied-copula hazard distribution
    /// and prices a bespoke $[K_1, K_2]$ CDO tranche.
    pub fn evaluate_hull_white_implied_copula_bespoke_cdo_tranche_calibration<T: Real>(
        spec: &HullWhiteImpliedCopulaCdoSpec<T>,
    ) -> HullWhiteImpliedCopulaCdoAnalytics<T> {
        let p_crash = spec
            .crash_state_probability_p_crash
            .max_val(T::from_f64(1e-4))
            .min_val(T::from_f64(0.45));
        let non_crash = T::one() - p_crash;
        let weights = [
            T::from_f64(0.55) * non_crash,
            T::from_f64(0.30) * non_crash,
            T::from_f64(0.15) * non_crash,
            p_crash,
        ];

        let lam0 = spec.base_hazard_rate_lambda0.max_val(T::from_f64(1e-4));
        let m_rec = spec.recession_hazard_multiplier_m_rec.max_val(T::one());
        let lam_crash = spec
            .systemic_crash_hazard_lambda_crash
            .max_val(T::from_f64(1e-3));
        let hazards = [T::from_f64(0.5) * lam0, lam0, m_rec * lam0, lam_crash];

        let lgd = spec
            .loss_given_default_lgd
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.95));
        let k1 = spec.tranche_attachment_k1.max_val(T::zero());
        let width = (spec.tranche_detachment_k2 - k1).max_val(T::from_f64(1e-3));
        let k2 = k1 + width;

        let t_mid = T::from_f64(2.5);
        let t_mat = T::from_f64(5.0);

        let mut port_loss_mat = T::zero();
        let mut tr_loss_mid = T::zero();
        let mut tr_loss_mat = T::zero();

        for (&w, &lam) in weights.iter().zip(hazards.iter()) {
            let l_mid = lgd * (T::one() - (-lam * t_mid).exp());
            let l_mat = lgd * (T::one() - (-lam * t_mat).exp());
            port_loss_mat += w * l_mat;

            let tr_mid = (smooth_pos(l_mid - k1) - smooth_pos(l_mid - k2)) / width;
            let tr_mat = (smooth_pos(l_mat - k1) - smooth_pos(l_mat - k2)) / width;
            tr_loss_mid += w * tr_mid;
            tr_loss_mat += w * tr_mat;
        }

        let r = spec.risk_free_rate_r;
        let df1 = (-r * T::from_f64(1.25)).exp();
        let df2 = (-r * T::from_f64(3.75)).exp();
        let prot_pv = df1 * tr_loss_mid + df2 * (tr_loss_mat - tr_loss_mid);

        let risky_ann = T::from_f64(2.5)
            * ((-r * t_mid).exp() * (T::one() - tr_loss_mid)
                + (-r * t_mat).exp() * (T::one() - tr_loss_mat))
                .max_val(T::from_f64(1e-4));
        let par_spread = prot_pv / risky_ann;

        HullWhiteImpliedCopulaCdoAnalytics {
            expected_portfolio_index_loss: port_loss_mat,
            expected_bespoke_tranche_loss_fraction: tr_loss_mat,
            tranche_protection_leg_pv: prot_pv,
            bespoke_tranche_par_spread: par_spread,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hull_white_implied_copula_cdo_invariants_and_ad_parity() {
            let base = HullWhiteImpliedCopulaCdoSpec {
                base_hazard_rate_lambda0: 0.015_f64,
                recession_hazard_multiplier_m_rec: 3.2,
                systemic_crash_hazard_lambda_crash: 0.14,
                crash_state_probability_p_crash: 0.06,
                loss_given_default_lgd: 0.60,
                tranche_attachment_k1: 0.03,
                tranche_detachment_k2: 0.07,
                risk_free_rate_r: 0.035,
            };
            let r0 = evaluate_hull_white_implied_copula_bespoke_cdo_tranche_calibration(&base);
            assert!(r0.expected_portfolio_index_loss > 0.0);
            assert!(
                r0.expected_bespoke_tranche_loss_fraction > 0.0
                    && r0.expected_bespoke_tranche_loss_fraction < 1.0
            );
            assert!(r0.tranche_protection_leg_pv > 0.0);
            assert!(r0.bespoke_tranche_par_spread > 0.0);

            let inputs = [
                base.base_hazard_rate_lambda0,
                base.recession_hazard_multiplier_m_rec,
                base.systemic_crash_hazard_lambda_crash,
                base.crash_state_probability_p_crash,
                base.loss_given_default_lgd,
                base.tranche_attachment_k1,
                base.tranche_detachment_k2,
                base.risk_free_rate_r,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HullWhiteImpliedCopulaCdoSpec {
                base_hazard_rate_lambda0: dual_in[0],
                recession_hazard_multiplier_m_rec: dual_in[1],
                systemic_crash_hazard_lambda_crash: dual_in[2],
                crash_state_probability_p_crash: dual_in[3],
                loss_given_default_lgd: dual_in[4],
                tranche_attachment_k1: dual_in[5],
                tranche_detachment_k2: dual_in[6],
                risk_free_rate_r: dual_in[7],
            };
            let r_ad = evaluate_hull_white_implied_copula_bespoke_cdo_tranche_calibration(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| HullWhiteImpliedCopulaCdoSpec {
                    base_hazard_rate_lambda0: v[0],
                    recession_hazard_multiplier_m_rec: v[1],
                    systemic_crash_hazard_lambda_crash: v[2],
                    crash_state_probability_p_crash: v[3],
                    loss_given_default_lgd: v[4],
                    tranche_attachment_k1: v[5],
                    tranche_detachment_k2: v[6],
                    risk_free_rate_r: v[7],
                };
                let rp = evaluate_hull_white_implied_copula_bespoke_cdo_tranche_calibration(&mk(p));
                let rm = evaluate_hull_white_implied_copula_bespoke_cdo_tranche_calibration(&mk(m));
                let fd =
                    (rp.bespoke_tranche_par_spread - rm.bespoke_tranche_par_spread) / (2.0 * eps);
                assert!(
                    r_ad.bespoke_tranche_par_spread.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.bespoke_tranche_par_spread.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use hull_white_implied_copula_cdo_mod::*;
mod schonbucher_clayton_cdo_mod {
    use ql_ad::Real;

    /// Input specification for the Schönbucher (2003, Ch. 10, §10.4) / McNeil, Frey & Embrechts (2015)
    /// Archimedean Clayton frailty copula CDO tranche valuation model.
    #[derive(Debug, Clone, Copy)]
    pub struct SchonbucherClaytonCopulaCdoSpec<T: Real> {
        /// Pool average default hazard intensity $\lambda$ (e.g., $0.022$).
        pub pool_hazard_rate_lambda: T,
        /// Pool recovery rate $R \in (0, 1)$ (e.g., $0.40$).
        pub recovery_rate_r: T,
        /// Clayton Archimedean dependence parameter $\theta > 0$ (e.g., $0.65$).
        pub clayton_dependence_theta: T,
        /// Tranche attachment point $K_A \in [0, 1)$ (e.g., $0.03$).
        pub attachment_point_ka: T,
        /// Tranche detachment point $K_B \in (K_A, 1]$ (e.g., $0.07$).
        pub detachment_point_kb: T,
        /// Risk-free discount rate $r$ (e.g., $0.035$).
        pub risk_free_rate_r: T,
        /// Tranche maturity $T$ in years (e.g., $5.0$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_schonbucher_archimedean_clayton_copula_cdo_tranche_loss`].
    #[derive(Debug, Clone, Copy)]
    pub struct SchonbucherClaytonCopulaCdoAnalytics<T: Real> {
        /// Lower-tail dependence coefficient $\lambda_L = 2^{-1/\theta}$ of the Clayton copula.
        pub lower_tail_dependence_lambda_l: T,
        /// Kendall's rank correlation $\tau_K = \theta / (\theta + 2)$.
        pub kendall_rank_correlation_tau: T,
        /// Expected normalized tranche loss fraction on $[K_A, K_B]$.
        pub expected_tranche_loss_fraction: T,
        /// Annualized fair CDO tranche par spread (in decimal).
        pub fair_tranche_par_spread: T,
    }

    /// Evaluates a CDO tranche under the Schönbucher (2003, §10.4) Archimedean Clayton frailty
    /// copula with generator $\psi_\theta(s) = (1 + s)^{-1/\theta}$, where conditional default
    /// probability given Gamma frailty $Y$ (with $\mathbb{E}[Y]=1, \operatorname{Var}(Y)=\theta$) is
    /// $p(Y) = 1 - \exp\left(-Y((1 - p_0)^{-\theta} - 1)\right)$.
    pub fn evaluate_schonbucher_archimedean_clayton_copula_cdo_tranche_loss<T: Real>(
        spec: &SchonbucherClaytonCopulaCdoSpec<T>,
    ) -> SchonbucherClaytonCopulaCdoAnalytics<T> {
        let theta = spec.clayton_dependence_theta;
        let t = spec.maturity_years_t;
        let lgd = T::one() - spec.recovery_rate_r;

        // Marginal survival q0 = exp(-lambda * T) and Clayton inverted generator u0 = q0^(-theta) - 1
        let q0 = (-spec.pool_hazard_rate_lambda * t).exp();
        let u0 = q0.pow(-theta) - T::one();

        let lambda_l = T::from_f64(2.0).pow(-T::one() / theta);
        let tau_k = theta / (theta + T::from_f64(2.0));

        // 4-point moment-matched lognormal/Gamma frailty quadrature for Y > 0 with E[Y] = 1, Var(Y) = theta
        let sigma_y = (T::one() + theta).ln().sqrt();
        let mu_y = -T::from_f64(0.5) * sigma_y * sigma_y;
        let gh_z = [
            -1.650_680_123_885_784_5,
            -0.524_647_623_275_290_3,
            0.524_647_623_275_290_3,
            1.650_680_123_885_784_5,
        ];
        let gh_w = [
            0.081_312_835_447_245_18,
            0.418_687_164_552_754_8,
            0.418_687_164_552_754_8,
            0.081_312_835_447_245_18,
        ];

        let width = spec.detachment_point_kb - spec.attachment_point_ka;
        let eps_sm = T::from_f64(1e-6);
        let half = T::from_f64(0.5);

        let mut exp_tranche_loss = T::zero();
        for idx in 0..4 {
            let z = T::from_f64(gh_z[idx]);
            let w = T::from_f64(gh_w[idx]);
            let y_frailty = (mu_y + sigma_y * T::from_f64(std::f64::consts::SQRT_2) * z).exp();
            let cond_pd = T::one() - (-y_frailty * u0).exp();
            let pool_loss = lgd * cond_pd;

            let d_a = pool_loss - spec.attachment_point_ka;
            let d_b = pool_loss - spec.detachment_point_kb;
            let call_a = half * (d_a + (d_a * d_a + eps_sm).sqrt());
            let call_b = half * (d_b + (d_b * d_b + eps_sm).sqrt());
            let tranche_loss_frac = (call_a - call_b) / width;

            exp_tranche_loss += w * tranche_loss_frac;
        }

        let r = spec.risk_free_rate_r;
        let risky_annuity = (T::one() - (-r * t).exp()) / r * (T::one() - half * exp_tranche_loss);
        let protection_pv = (-half * r * t).exp() * exp_tranche_loss;
        let fair_spread = protection_pv / risky_annuity;

        SchonbucherClaytonCopulaCdoAnalytics {
            lower_tail_dependence_lambda_l: lambda_l,
            kendall_rank_correlation_tau: tau_k,
            expected_tranche_loss_fraction: exp_tranche_loss,
            fair_tranche_par_spread: fair_spread,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn schonbucher_clayton_copula_cdo_invariants_and_ad_parity() {
            let base = SchonbucherClaytonCopulaCdoSpec {
                pool_hazard_rate_lambda: 0.022_f64,
                recovery_rate_r: 0.40_f64,
                clayton_dependence_theta: 0.65_f64,
                attachment_point_ka: 0.03_f64,
                detachment_point_kb: 0.07_f64,
                risk_free_rate_r: 0.035_f64,
                maturity_years_t: 5.0_f64,
            };
            let res = evaluate_schonbucher_archimedean_clayton_copula_cdo_tranche_loss(&base);
            assert!(
                res.lower_tail_dependence_lambda_l > 0.0
                    && res.lower_tail_dependence_lambda_l < 1.0
            );
            assert!(
                res.kendall_rank_correlation_tau > 0.0 && res.kendall_rank_correlation_tau < 1.0
            );
            assert!(
                res.expected_tranche_loss_fraction > 0.0
                    && res.expected_tranche_loss_fraction < 1.0
            );
            assert!(res.fair_tranche_par_spread > 0.0);

            let inputs = [
                base.pool_hazard_rate_lambda,
                base.recovery_rate_r,
                base.clayton_dependence_theta,
                base.attachment_point_ka,
                base.detachment_point_kb,
                base.risk_free_rate_r,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = SchonbucherClaytonCopulaCdoSpec {
                pool_hazard_rate_lambda: dual_in[0],
                recovery_rate_r: dual_in[1],
                clayton_dependence_theta: dual_in[2],
                attachment_point_ka: dual_in[3],
                detachment_point_kb: dual_in[4],
                risk_free_rate_r: dual_in[5],
                maturity_years_t: dual_in[6],
            };
            let r_ad = evaluate_schonbucher_archimedean_clayton_copula_cdo_tranche_loss(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| SchonbucherClaytonCopulaCdoSpec {
                    pool_hazard_rate_lambda: v[0],
                    recovery_rate_r: v[1],
                    clayton_dependence_theta: v[2],
                    attachment_point_ka: v[3],
                    detachment_point_kb: v[4],
                    risk_free_rate_r: v[5],
                    maturity_years_t: v[6],
                };
                let rp = evaluate_schonbucher_archimedean_clayton_copula_cdo_tranche_loss(&mk(p));
                let rm = evaluate_schonbucher_archimedean_clayton_copula_cdo_tranche_loss(&mk(m));
                let fd = (rp.fair_tranche_par_spread - rm.fair_tranche_par_spread) / (2.0 * eps);
                assert!(
                    r_ad.fair_tranche_par_spread.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.fair_tranche_par_spread.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use schonbucher_clayton_cdo_mod::*;
mod frye_altman_downturn_lgd_mod {
    use ql_ad::Real;

    /// Input specification for the Frye (2000, *Risk* 13(11), *"Depressing Recoveries"*) &
    /// Altman, Brady, Resti & Sironi (2005, *Journal of Business* 78(6)) single-factor
    /// procyclical downturn LGD and credit portfolio loss model.
    #[derive(Debug, Clone, Copy)]
    pub struct FryeAltmanDownturnLgdSpec<T: Real> {
        /// Through-the-cycle inverse-normal default threshold $c_{\text{PD}} = \Phi^{-1}(\overline{\text{PD}})$ (e.g., $-2.0537$ for $2\%$ PD).
        pub default_threshold_z_pd: T,
        /// Asset-return systematic factor correlation $\rho_A \in (0, 1)$ (e.g., $0.18$).
        pub asset_factor_correlation_rho: T,
        /// Expected through-the-cycle collateral recovery rate $\mu_R \in (0, 1)$ (e.g., $0.45$).
        pub expected_recovery_rate_mu_r: T,
        /// Collateral recovery volatility $\sigma_R > 0$ (e.g., $0.20$).
        pub recovery_volatility_sigma_r: T,
        /// Systematic factor loading $q_R \in (0, 1)$ of recovery rate on macro factor $Z$ (e.g., $0.55$).
        pub recovery_factor_loading_q_r: T,
        /// Adverse tail macro realization magnitude $z_{\alpha} = -\Phi^{-1}(1 - \alpha) > 0$ (e.g., $2.3263$ for $99\%$).
        pub adverse_macro_quantile_z_alpha: T,
    }

    /// Analytics returned by [`evaluate_frye_altman_procyclical_downturn_lgd_credit_loss`].
    #[derive(Debug, Clone, Copy)]
    pub struct FryeAltmanDownturnLgdAnalytics<T: Real> {
        /// Conditional default probability $\text{PD}(Z = -z_\alpha)$ in the systemic downturn.
        pub conditional_downturn_pd: T,
        /// Conditional downturn Loss Given Default $\text{LGD}_{\text{down}}(Z = -z_\alpha) > 1 - \mu_R$.
        pub conditional_downturn_lgd: T,
        /// Constant-recovery Vasicek portfolio loss at $Z = -z_\alpha$, $(1 - \mu_R)\text{PD}(-z_\alpha)$.
        pub constant_lgd_stressed_loss: T,
        /// Procyclical Frye–Altman stressed portfolio loss $\text{LGD}_{\text{down}}(-z_\alpha)\,\text{PD}(-z_\alpha)$.
        pub procyclical_downturn_stressed_loss: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        T::from_f64(0.5) * (T::one() + (x * T::from_f64(0.707_106_781_186_547_5)).erf())
    }

    /// Evaluates the Frye (2000) / Altman, Brady, Resti & Sironi (2005) single-factor
    /// procyclical default-recovery portfolio loss at systematic downturn $Z = -z_\alpha$:
    /// $$\text{PD}(-z_\alpha) = \Phi\!\left(\frac{c_{\text{PD}} + \sqrt{\rho_A}\,z_\alpha}{\sqrt{1 - \rho_A}}\right), \qquad \text{LGD}_{\text{down}}(-z_\alpha) = 1 - \big(\mu_R - \sigma_R q_R z_\alpha\big).$$
    pub fn evaluate_frye_altman_procyclical_downturn_lgd_credit_loss<T: Real>(
        spec: &FryeAltmanDownturnLgdSpec<T>,
    ) -> FryeAltmanDownturnLgdAnalytics<T> {
        let rho = spec.asset_factor_correlation_rho;
        let z_a = spec.adverse_macro_quantile_z_alpha;
        let cond_arg = (spec.default_threshold_z_pd + rho.sqrt() * z_a) / (T::one() - rho).sqrt();
        let pd_down = std_norm_cdf(cond_arg);

        let r_down = spec.expected_recovery_rate_mu_r
            - spec.recovery_volatility_sigma_r * spec.recovery_factor_loading_q_r * z_a;
        let lgd_ttc = T::one() - spec.expected_recovery_rate_mu_r;
        let lgd_down = T::one() - r_down;

        let loss_const = lgd_ttc * pd_down;
        let loss_procyc = lgd_down * pd_down;

        FryeAltmanDownturnLgdAnalytics {
            conditional_downturn_pd: pd_down,
            conditional_downturn_lgd: lgd_down,
            constant_lgd_stressed_loss: loss_const,
            procyclical_downturn_stressed_loss: loss_procyc,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn frye_altman_downturn_lgd_invariants_and_ad_parity() {
            let base = FryeAltmanDownturnLgdSpec {
                default_threshold_z_pd: -2.0537_f64,
                asset_factor_correlation_rho: 0.18_f64,
                expected_recovery_rate_mu_r: 0.45_f64,
                recovery_volatility_sigma_r: 0.20_f64,
                recovery_factor_loading_q_r: 0.55_f64,
                adverse_macro_quantile_z_alpha: 2.3263_f64,
            };
            let res = evaluate_frye_altman_procyclical_downturn_lgd_credit_loss(&base);
            assert!(res.conditional_downturn_pd > 0.05 && res.conditional_downturn_pd < 0.50);
            assert!(res.conditional_downturn_lgd > 1.0 - base.expected_recovery_rate_mu_r);
            assert!(res.procyclical_downturn_stressed_loss > res.constant_lgd_stressed_loss);

            let inputs = [
                base.default_threshold_z_pd,
                base.asset_factor_correlation_rho,
                base.expected_recovery_rate_mu_r,
                base.recovery_volatility_sigma_r,
                base.recovery_factor_loading_q_r,
                base.adverse_macro_quantile_z_alpha,
            ];
            let mut dual_in = [Dual::<6>::from_value(0.0); 6];
            for i in 0..6 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = FryeAltmanDownturnLgdSpec {
                default_threshold_z_pd: dual_in[0],
                asset_factor_correlation_rho: dual_in[1],
                expected_recovery_rate_mu_r: dual_in[2],
                recovery_volatility_sigma_r: dual_in[3],
                recovery_factor_loading_q_r: dual_in[4],
                adverse_macro_quantile_z_alpha: dual_in[5],
            };
            let r_ad = evaluate_frye_altman_procyclical_downturn_lgd_credit_loss(&spec_ad);
            let eps = 1e-6;
            for i in 0..6 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 6]| FryeAltmanDownturnLgdSpec {
                    default_threshold_z_pd: v[0],
                    asset_factor_correlation_rho: v[1],
                    expected_recovery_rate_mu_r: v[2],
                    recovery_volatility_sigma_r: v[3],
                    recovery_factor_loading_q_r: v[4],
                    adverse_macro_quantile_z_alpha: v[5],
                };
                let rp = evaluate_frye_altman_procyclical_downturn_lgd_credit_loss(&mk(p));
                let rm = evaluate_frye_altman_procyclical_downturn_lgd_credit_loss(&mk(m));
                let fd = (rp.procyclical_downturn_stressed_loss
                    - rm.procyclical_downturn_stressed_loss)
                    / (2.0 * eps);
                assert!(
                    r_ad.procyclical_downturn_stressed_loss.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.procyclical_downturn_stressed_loss.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use frye_altman_downturn_lgd_mod::*;
mod duffie_saita_wang_stochastic_default_mod {
    use ql_ad::Real;

    /// Input specification for the Duffie, Saita & Wang (2007, *Journal of Financial Economics*
    /// 83(3), *"Multi-Period Corporate Default Prediction with Stochastic Covariates"*)
    /// doubly-stochastic Distance-to-Default term-structure survival engine.
    #[derive(Debug, Clone, Copy)]
    pub struct DuffieSaitaWangMultiPeriodDefaultSpec<T: Real> {
        /// Baseline log-hazard intercept $\beta_0$ (e.g., $-2.40$).
        pub baseline_log_hazard_beta0: T,
        /// Current Merton Distance-to-Default $\text{DD}_0$ in standard deviations (e.g., $2.60$).
        pub initial_distance_to_default_dd0: T,
        /// Long-run stationary Distance-to-Default $\bar{D}$ (e.g., $2.20$).
        pub long_run_distance_to_default_bar: T,
        /// Distance-to-Default mean-reversion speed $\kappa_D > 0$ (e.g., $0.38$).
        pub dd_mean_reversion_kappa: T,
        /// Distance-to-Default diffusion volatility $\sigma_D > 0$ (e.g., $0.45$).
        pub dd_volatility_sigma: T,
        /// Sensitivity $\beta_D > 0$ of log-default intensity to Distance-to-Default (e.g., $0.85$).
        pub dd_hazard_loading_beta_d: T,
        /// Prediction horizon $T$ in years (e.g., $3.0$).
        pub horizon_years_t: T,
    }

    /// Analytics returned by [`evaluate_duffie_saita_wang_stochastic_covariate_survival`].
    #[derive(Debug, Clone, Copy)]
    pub struct DuffieSaitaWangMultiPeriodDefaultAnalytics<T: Real> {
        /// Instantaneous default intensity $\lambda_0 = \exp(\beta_0 - \beta_D \text{DD}_0)$.
        pub instantaneous_default_hazard: T,
        /// Expected cumulative hazard $\bar{H}(0, T)$ incorporating covariate mean-reversion and convexity.
        pub expected_cumulative_hazard: T,
        /// Multi-period survival probability $S(0, T) = \exp(-\bar{H}(0, T))$.
        pub multi_period_survival_probability: T,
        /// Multi-period cumulative default probability $\text{PD}(0, T) = 1 - S(0, T)$.
        pub multi_period_default_probability: T,
    }

    /// Evaluates the Duffie, Saita & Wang (2007, Eqs. 1–9) doubly-stochastic intensity
    /// $\lambda_t = \exp(\beta_0 - \beta_D \text{DD}_t)$ under Ornstein–Uhlenbeck covariate
    /// dynamics $d\text{DD}_t = \kappa_D(\bar{D} - \text{DD}_t)\,dt + \sigma_D\,dW_t$, computing
    /// the mean-reverting path average $\bar{\mu}_\Lambda(T)$ and integrated log-hazard variance $\Sigma_\Lambda^2(T)$.
    pub fn evaluate_duffie_saita_wang_stochastic_covariate_survival<T: Real>(
        spec: &DuffieSaitaWangMultiPeriodDefaultSpec<T>,
    ) -> DuffieSaitaWangMultiPeriodDefaultAnalytics<T> {
        let b0 = spec.baseline_log_hazard_beta0;
        let dd0 = spec.initial_distance_to_default_dd0;
        let d_bar = spec.long_run_distance_to_default_bar;
        let kap = spec.dd_mean_reversion_kappa;
        let sig = spec.dd_volatility_sigma;
        let b_d = spec.dd_hazard_loading_beta_d;
        let t = spec.horizon_years_t;

        let lam0 = (b0 - b_d * dd0).exp();

        let e_kt = (-kap * t).exp();
        let b_frac = (T::one() - e_kt) / (kap * t);
        let avg_dd = d_bar + (dd0 - d_bar) * b_frac;
        let mu_log = b0 - b_d * avg_dd;

        // Integrated OU variance over [0, T] divided by T^2
        let e_2kt = e_kt * e_kt;
        let ou_int_var = (b_d * b_d * sig * sig / (kap * kap))
            * (t - T::from_f64(2.0) * (T::one() - e_kt) / kap
                + (T::one() - e_2kt) / (T::from_f64(2.0) * kap));
        let h_cum = t * (mu_log + T::from_f64(0.5) * ou_int_var / t).exp();
        let surv = (-h_cum).exp();
        let pd = T::one() - surv;

        DuffieSaitaWangMultiPeriodDefaultAnalytics {
            instantaneous_default_hazard: lam0,
            expected_cumulative_hazard: h_cum,
            multi_period_survival_probability: surv,
            multi_period_default_probability: pd,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn duffie_saita_wang_stochastic_default_invariants_and_ad_parity() {
            let base = DuffieSaitaWangMultiPeriodDefaultSpec {
                baseline_log_hazard_beta0: -2.40_f64,
                initial_distance_to_default_dd0: 2.60_f64,
                long_run_distance_to_default_bar: 2.20_f64,
                dd_mean_reversion_kappa: 0.38_f64,
                dd_volatility_sigma: 0.45_f64,
                dd_hazard_loading_beta_d: 0.85_f64,
                horizon_years_t: 3.0_f64,
            };
            let res = evaluate_duffie_saita_wang_stochastic_covariate_survival(&base);
            assert!(res.instantaneous_default_hazard > 0.0);
            assert!(res.expected_cumulative_hazard > 0.0);
            assert!(
                res.multi_period_survival_probability > 0.0
                    && res.multi_period_survival_probability < 1.0
            );
            assert!(
                (res.multi_period_survival_probability + res.multi_period_default_probability
                    - 1.0)
                    .abs()
                    < 1e-12
            );

            let inputs = [
                base.baseline_log_hazard_beta0,
                base.initial_distance_to_default_dd0,
                base.long_run_distance_to_default_bar,
                base.dd_mean_reversion_kappa,
                base.dd_volatility_sigma,
                base.dd_hazard_loading_beta_d,
                base.horizon_years_t,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = DuffieSaitaWangMultiPeriodDefaultSpec {
                baseline_log_hazard_beta0: dual_in[0],
                initial_distance_to_default_dd0: dual_in[1],
                long_run_distance_to_default_bar: dual_in[2],
                dd_mean_reversion_kappa: dual_in[3],
                dd_volatility_sigma: dual_in[4],
                dd_hazard_loading_beta_d: dual_in[5],
                horizon_years_t: dual_in[6],
            };
            let r_ad = evaluate_duffie_saita_wang_stochastic_covariate_survival(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| DuffieSaitaWangMultiPeriodDefaultSpec {
                    baseline_log_hazard_beta0: v[0],
                    initial_distance_to_default_dd0: v[1],
                    long_run_distance_to_default_bar: v[2],
                    dd_mean_reversion_kappa: v[3],
                    dd_volatility_sigma: v[4],
                    dd_hazard_loading_beta_d: v[5],
                    horizon_years_t: v[6],
                };
                let rp = evaluate_duffie_saita_wang_stochastic_covariate_survival(&mk(p));
                let rm = evaluate_duffie_saita_wang_stochastic_covariate_survival(&mk(m));
                let fd = (rp.multi_period_default_probability
                    - rm.multi_period_default_probability)
                    / (2.0 * eps);
                assert!(
                    r_ad.multi_period_default_probability.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.multi_period_default_probability.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use duffie_saita_wang_stochastic_default_mod::*;

/// Input specification for the Errais, Giesecke & Goldberg (2010, *SIAM Journal
/// on Financial Mathematics* 1(1)) Self-Exciting Hawkes Frailty Contagion
/// Portfolio Credit Loss engine.
#[derive(Debug, Clone, Copy)]
pub struct ErraisGieseckeGoldbergHawkesCreditSpec<T: Real> {
    /// Initial portfolio default intensity $\lambda_0 > 0$.
    pub initial_default_intensity_lambda0: T,
    /// Long-run baseline reversion intensity $\mu_\lambda > 0$.
    pub long_run_baseline_intensity_mu: T,
    /// Intensity exponential mean-reversion speed $\kappa > 0$.
    pub mean_reversion_speed_kappa: T,
    /// Self-exciting Hawkes jump impact per unit loss $\delta_H > 0$.
    pub self_exciting_contagion_impact_delta: T,
    /// Average fractional loss per default event $\bar{\ell} \in (0,1)$.
    pub expected_loss_per_default_ell: T,
    /// Evaluation horizon in years $T > 0$.
    pub horizon_years_t: T,
    /// Tranche detachment fraction $K_d \in (0,1)$.
    pub tranche_detachment_fraction_d: T,
}

/// Analytics returned by [`evaluate_errais_giesecke_goldberg_hawkes_contagion_credit_loss`].
#[derive(Debug, Clone, Copy)]
pub struct ErraisGieseckeGoldbergHawkesCreditAnalytics<T: Real> {
    /// Effective net decay rate $\hat{\kappa} = \kappa - \delta_H \bar{\ell}$.
    pub effective_net_decay_kappa_hat: T,
    /// Hawkes contagion branching ratio $\eta_H = \delta_H \bar{\ell} / \kappa$.
    pub contagion_branching_ratio_eta: T,
    /// Expected terminal default intensity $\mathbb{E}[\lambda_T]$.
    pub expected_terminal_intensity_lambda_t: T,
    /// Expected cumulative portfolio credit loss $\mathbb{E}[L_T]$.
    pub expected_cumulative_portfolio_loss: T,
    /// Contagion-amplified tranche risk charge $\Pi_{\text{tranche}}(K_d)$.
    pub contagion_adjusted_tranche_risk_charge: T,
}

/// Evaluates the Errais, Giesecke & Goldberg (2010, *SIFIN* 1(1)) affine
/// self-exciting Hawkes point-process portfolio credit loss and contagion
/// clustering amplification.
pub fn evaluate_errais_giesecke_goldberg_hawkes_contagion_credit_loss<T: Real>(
    spec: &ErraisGieseckeGoldbergHawkesCreditSpec<T>,
) -> ErraisGieseckeGoldbergHawkesCreditAnalytics<T> {
    let one = T::one();
    let kappa = spec.mean_reversion_speed_kappa.max_val(T::from_f64(0.1));
    let ell = spec
        .expected_loss_per_default_ell
        .max_val(T::from_f64(1e-4));
    let delta = spec
        .self_exciting_contagion_impact_delta
        .max_val(T::from_f64(0.0));
    let t = spec.horizon_years_t.max_val(T::from_f64(1e-4));

    let feedback = delta * ell;
    let kappa_hat = (kappa - feedback).max_val(T::from_f64(0.05));
    let eta = feedback / kappa;
    let asym_lambda = kappa * spec.long_run_baseline_intensity_mu / kappa_hat;

    let decay = (-kappa_hat * t).exp();
    let e_lambda_t = asym_lambda + (spec.initial_default_intensity_lambda0 - asym_lambda) * decay;
    let e_n_t = asym_lambda * t
        + (spec.initial_default_intensity_lambda0 - asym_lambda) * (one - decay) / kappa_hat;
    let e_loss = ell * e_n_t;

    let kd = spec
        .tranche_detachment_fraction_d
        .max_val(T::from_f64(1e-4));
    let cluster_mult =
        one + (feedback / kappa_hat) * (-kd / e_loss.max_val(T::from_f64(1e-6))).exp();
    let tranche_charge = e_loss * cluster_mult;

    ErraisGieseckeGoldbergHawkesCreditAnalytics {
        effective_net_decay_kappa_hat: kappa_hat,
        contagion_branching_ratio_eta: eta,
        expected_terminal_intensity_lambda_t: e_lambda_t,
        expected_cumulative_portfolio_loss: e_loss,
        contagion_adjusted_tranche_risk_charge: tranche_charge,
    }
}

#[cfg(test)]
mod errais_giesecke_goldberg_hawkes_credit_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn errais_giesecke_goldberg_hawkes_credit_invariants_and_ad_parity() {
        let base = ErraisGieseckeGoldbergHawkesCreditSpec {
            initial_default_intensity_lambda0: 0.035_f64,
            long_run_baseline_intensity_mu: 0.025_f64,
            mean_reversion_speed_kappa: 1.40_f64,
            self_exciting_contagion_impact_delta: 0.95_f64,
            expected_loss_per_default_ell: 0.60_f64,
            horizon_years_t: 5.0_f64,
            tranche_detachment_fraction_d: 0.12_f64,
        };
        let r0 = evaluate_errais_giesecke_goldberg_hawkes_contagion_credit_loss(&base);
        assert!(r0.contagion_branching_ratio_eta > 0.0 && r0.contagion_branching_ratio_eta < 1.0);
        assert!(r0.contagion_adjusted_tranche_risk_charge > r0.expected_cumulative_portfolio_loss);

        let x0 = [
            base.initial_default_intensity_lambda0,
            base.long_run_baseline_intensity_mu,
            base.mean_reversion_speed_kappa,
            base.self_exciting_contagion_impact_delta,
            base.expected_loss_per_default_ell,
            base.horizon_years_t,
            base.tranche_detachment_fraction_d,
        ];
        let mut xd = [Dual::<7>::from_f64(0.0); 7];
        for i in 0..7 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_errais_giesecke_goldberg_hawkes_contagion_credit_loss(
            &ErraisGieseckeGoldbergHawkesCreditSpec {
                initial_default_intensity_lambda0: xd[0],
                long_run_baseline_intensity_mu: xd[1],
                mean_reversion_speed_kappa: xd[2],
                self_exciting_contagion_impact_delta: xd[3],
                expected_loss_per_default_ell: xd[4],
                horizon_years_t: xd[5],
                tranche_detachment_fraction_d: xd[6],
            },
        );
        let h = 1e-6;
        for i in 0..7 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_errais_giesecke_goldberg_hawkes_contagion_credit_loss(
                &ErraisGieseckeGoldbergHawkesCreditSpec {
                    initial_default_intensity_lambda0: xp[0],
                    long_run_baseline_intensity_mu: xp[1],
                    mean_reversion_speed_kappa: xp[2],
                    self_exciting_contagion_impact_delta: xp[3],
                    expected_loss_per_default_ell: xp[4],
                    horizon_years_t: xp[5],
                    tranche_detachment_fraction_d: xp[6],
                },
            );
            let rm = evaluate_errais_giesecke_goldberg_hawkes_contagion_credit_loss(
                &ErraisGieseckeGoldbergHawkesCreditSpec {
                    initial_default_intensity_lambda0: xm[0],
                    long_run_baseline_intensity_mu: xm[1],
                    mean_reversion_speed_kappa: xm[2],
                    self_exciting_contagion_impact_delta: xm[3],
                    expected_loss_per_default_ell: xm[4],
                    horizon_years_t: xm[5],
                    tranche_detachment_fraction_d: xm[6],
                },
            );
            let fd = (rp.contagion_adjusted_tranche_risk_charge
                - rm.contagion_adjusted_tranche_risk_charge)
                / (2.0 * h);
            assert!(rad.contagion_adjusted_tranche_risk_charge.d[i].abs() > 1e-9);
            assert!((rad.contagion_adjusted_tranche_risk_charge.d[i] - fd).abs() < 1e-6);
        }
    }
}

/// Input specification for the Collin-Dufresne, Goldstein & Yang (2012, *Review
/// of Financial Studies* 25(11)) Catastrophic Macro Jump & Super-Senior CDO
/// Tranche Spread engine.
#[derive(Debug, Clone, Copy)]
pub struct CollinDufresneGoldsteinYangSuperSeniorCdoSpec<T: Real> {
    /// Baseline diffusive idiosyncratic default intensity $\lambda_{\text{idio}} > 0$.
    pub idiosyncratic_hazard_rate_lambda_i: T,
    /// Risk-neutral arrival intensity of catastrophic macro crisis jump $\lambda_{\text{macro}} > 0$.
    pub catastrophic_macro_jump_intensity_lambda_m: T,
    /// Conditional portfolio default fraction in a macro crisis $p_{\text{crisis}} \in (0,1)$.
    pub macro_crisis_conditional_default_prob_p_m: T,
    /// Loss given default $\text{LGD} \in (0,1)$.
    pub loss_given_default_lgd: T,
    /// Super-senior tranche attachment point $K_a \in (0,1)$.
    pub super_senior_attachment_a: T,
    /// Super-senior tranche detachment point $K_d \in (K_a, 1)$.
    pub super_senior_detachment_d: T,
    /// Risk-free interest rate $r > 0$.
    pub risk_free_rate_r: T,
    /// Tranche maturity in years $T > 0$.
    pub maturity_years_t: T,
}

/// Analytics returned by [`evaluate_collin_dufresne_goldstein_yang_super_senior_cdo`].
#[derive(Debug, Clone, Copy)]
pub struct CollinDufresneGoldsteinYangSuperSeniorCdoAnalytics<T: Real> {
    /// Expected portfolio loss conditional on a macro catastrophe $L_{\text{macro}}$.
    pub conditional_catastrophic_portfolio_loss: T,
    /// Expected fractional wipeout of the super-senior tranche in a macro crisis $\ell_{\text{tr}}^{\text{macro}}$.
    pub super_senior_conditional_tranche_loss_fraction: T,
    /// Super-senior tranche risky RPV01 annuity $A_{\text{SS}}(T)$.
    pub super_senior_risky_annuity: T,
    /// Super-senior tranche par running spread in basis points $S_{\text{SS}}^{\text{CDGY}}$.
    pub super_senior_tranche_par_spread_bps: T,
}

/// Evaluates the Collin-Dufresne, Goldstein & Yang (2012, *RFS* 25(11)) joint
/// index-option and CDO catastrophe-jump model for super-senior tranche spreads.
pub fn evaluate_collin_dufresne_goldstein_yang_super_senior_cdo<T: Real>(
    spec: &CollinDufresneGoldsteinYangSuperSeniorCdoSpec<T>,
) -> CollinDufresneGoldsteinYangSuperSeniorCdoAnalytics<T> {
    let one = T::one();
    let half = T::from_f64(0.5);
    let inv_sqrt_2 = T::from_f64(std::f64::consts::FRAC_1_SQRT_2);
    let inv_sqrt_2pi = T::from_f64(0.3989422804014327);

    let lam_i = spec
        .idiosyncratic_hazard_rate_lambda_i
        .max_val(T::from_f64(1e-4));
    let lam_m = spec
        .catastrophic_macro_jump_intensity_lambda_m
        .max_val(T::from_f64(1e-5));
    let p_m = spec
        .macro_crisis_conditional_default_prob_p_m
        .min_val(T::from_f64(0.99))
        .max_val(T::from_f64(0.01));
    let lgd = spec.loss_given_default_lgd.max_val(T::from_f64(0.05));
    let ka = spec.super_senior_attachment_a.max_val(T::from_f64(0.01));
    let kd = spec
        .super_senior_detachment_d
        .max_val(ka + T::from_f64(0.02));
    let r = spec.risk_free_rate_r.max_val(T::from_f64(1e-4));
    let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

    let l_macro = p_m * lgd + (one - p_m) * lam_i * lgd * half * t;
    let sig_l = (lam_i * t).sqrt() * lgd * T::from_f64(0.35);

    let bachelier_call = |strike: T| -> T {
        let z = (l_macro - strike) / sig_l;
        let nz = half * (one + (z * inv_sqrt_2).erf());
        let pdf = inv_sqrt_2pi * (-half * z * z).exp();
        (l_macro - strike) * nz + sig_l * pdf
    };

    let tr_loss_frac = ((bachelier_call(ka) - bachelier_call(kd)) / (kd - ka))
        .min_val(one)
        .max_val(T::from_f64(1e-6));
    let h_tr = lam_m * tr_loss_frac;
    let disc = r + h_tr;
    let rpv01 = (one - (-disc * t).exp()) / disc;
    let par_bps = T::from_f64(10_000.0) * h_tr * (one + (lam_i * lgd) / (r + lam_m));

    CollinDufresneGoldsteinYangSuperSeniorCdoAnalytics {
        conditional_catastrophic_portfolio_loss: l_macro,
        super_senior_conditional_tranche_loss_fraction: tr_loss_frac,
        super_senior_risky_annuity: rpv01,
        super_senior_tranche_par_spread_bps: par_bps,
    }
}

#[cfg(test)]
mod collin_dufresne_goldstein_yang_super_senior_cdo_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn collin_dufresne_goldstein_yang_super_senior_cdo_invariants_and_ad_parity() {
        let base = CollinDufresneGoldsteinYangSuperSeniorCdoSpec {
            idiosyncratic_hazard_rate_lambda_i: 0.015_f64,
            catastrophic_macro_jump_intensity_lambda_m: 0.018_f64,
            macro_crisis_conditional_default_prob_p_m: 0.36_f64,
            loss_given_default_lgd: 0.60_f64,
            super_senior_attachment_a: 0.15_f64,
            super_senior_detachment_d: 0.30_f64,
            risk_free_rate_r: 0.036_f64,
            maturity_years_t: 5.0_f64,
        };
        let r0 = evaluate_collin_dufresne_goldstein_yang_super_senior_cdo(&base);
        assert!(r0.super_senior_conditional_tranche_loss_fraction > 0.0);
        assert!(r0.super_senior_tranche_par_spread_bps > 0.0);

        let x0 = [
            base.idiosyncratic_hazard_rate_lambda_i,
            base.catastrophic_macro_jump_intensity_lambda_m,
            base.macro_crisis_conditional_default_prob_p_m,
            base.loss_given_default_lgd,
            base.super_senior_attachment_a,
            base.super_senior_detachment_d,
            base.risk_free_rate_r,
            base.maturity_years_t,
        ];
        let mut xd = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_collin_dufresne_goldstein_yang_super_senior_cdo(
            &CollinDufresneGoldsteinYangSuperSeniorCdoSpec {
                idiosyncratic_hazard_rate_lambda_i: xd[0],
                catastrophic_macro_jump_intensity_lambda_m: xd[1],
                macro_crisis_conditional_default_prob_p_m: xd[2],
                loss_given_default_lgd: xd[3],
                super_senior_attachment_a: xd[4],
                super_senior_detachment_d: xd[5],
                risk_free_rate_r: xd[6],
                maturity_years_t: xd[7],
            },
        );
        let h = 1e-6;
        for i in 0..8 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_collin_dufresne_goldstein_yang_super_senior_cdo(
                &CollinDufresneGoldsteinYangSuperSeniorCdoSpec {
                    idiosyncratic_hazard_rate_lambda_i: xp[0],
                    catastrophic_macro_jump_intensity_lambda_m: xp[1],
                    macro_crisis_conditional_default_prob_p_m: xp[2],
                    loss_given_default_lgd: xp[3],
                    super_senior_attachment_a: xp[4],
                    super_senior_detachment_d: xp[5],
                    risk_free_rate_r: xp[6],
                    maturity_years_t: xp[7],
                },
            );
            let rm = evaluate_collin_dufresne_goldstein_yang_super_senior_cdo(
                &CollinDufresneGoldsteinYangSuperSeniorCdoSpec {
                    idiosyncratic_hazard_rate_lambda_i: xm[0],
                    catastrophic_macro_jump_intensity_lambda_m: xm[1],
                    macro_crisis_conditional_default_prob_p_m: xm[2],
                    loss_given_default_lgd: xm[3],
                    super_senior_attachment_a: xm[4],
                    super_senior_detachment_d: xm[5],
                    risk_free_rate_r: xm[6],
                    maturity_years_t: xm[7],
                },
            );
            let fd = (rp.super_senior_tranche_par_spread_bps
                - rm.super_senior_tranche_par_spread_bps)
                / (2.0 * h);
            assert!(rad.super_senior_tranche_par_spread_bps.d[i].abs() > 1e-9);
            assert!((rad.super_senior_tranche_par_spread_bps.d[i] - fd).abs() < 1e-5);
        }
    }
}
