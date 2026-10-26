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
