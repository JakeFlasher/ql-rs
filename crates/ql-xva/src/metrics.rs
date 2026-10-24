use crate::NpvCube;
use ql_ad::Real;
use ql_curves::YieldTermStructure;
#[derive(Copy, Clone, Debug)]
pub struct CsaSpec {
    pub threshold_cpty: f64,
    pub threshold_self: f64,
    pub mta: f64,
}
impl CsaSpec {
    /// Remaining signed uncollateralized exposure for a netting set.
    ///
    /// Variation margin is called only when the positive/negative exposure
    /// exceeds the threshold plus the minimum transfer amount.  Once a call
    /// is made the residual is the threshold (not the collateral call itself).
    pub fn collateralized_exposure<T: Real>(&self, x: T) -> T {
        if x.val() >= 0.0 {
            collateralized_exposure(x, self.threshold_cpty, self.mta)
        } else {
            -collateralized_exposure(-x, self.threshold_self, self.mta)
        }
    }
}

/// Positive remaining exposure after a one-way CSA margin call.
///
/// The branch is selected using the primal value while all arithmetic in the
/// selected branch remains in `T`, preserving Dual/Dual2/TapeVar tangents.
/// For a zero-threshold, zero-MTA CSA this is identically zero for positive
/// mark-to-market values, as expected for continuous variation margin.
pub fn collateralized_exposure<T: Real>(mtm: T, threshold_cpty: f64, mta: f64) -> T {
    let e = mtm.max_val(T::zero());
    let threshold = T::from_f64(threshold_cpty.max(0.0));
    let excess = e - threshold;
    if excess.val() > mta.max(0.0) {
        e.min_val(threshold)
    } else {
        e
    }
}
#[derive(Clone, Debug)]
pub struct XvaReport<T: Real> {
    pub ee: Vec<T>,
    pub ene: Vec<T>,
    pub pfe_95: Vec<f64>,
    pub epe: T,
    pub cva: T,
    pub dva: T,
    pub fva: T,
}

/// Expected positive exposure by simulation time, aggregated across all
/// netting sets.  The result is intentionally `f64` for charting and report
/// consumers; the AD-preserving XVA valuation remains available through
/// `compute_xva`.
pub fn expected_positive_exposure<T: Real>(cube: &NpvCube<T>) -> Vec<f64> {
    (0..cube.times.len())
        .map(|k| {
            let mut sum = 0.0;
            for path in 0..cube.n_paths {
                let mut value = 0.0;
                for ns in 0..cube.n_netting_sets {
                    // Legal netting sets may not offset one another.  Floor
                    // each set before summing to the counterparty profile.
                    value += cube.get(ns, k, path).val().max(0.0);
                }
                sum += value;
            }
            sum / cube.n_paths.max(1) as f64
        })
        .collect()
}

/// Empirical potential future exposure by simulation time.
pub fn potential_future_exposure<T: Real>(cube: &NpvCube<T>, quantile: f64) -> Vec<f64> {
    let q = quantile.clamp(0.0, 1.0);
    (0..cube.times.len())
        .map(|k| {
            let mut values = Vec::with_capacity(cube.n_paths);
            for path in 0..cube.n_paths {
                let mut value = 0.0;
                for ns in 0..cube.n_netting_sets {
                    value += cube.get(ns, k, path).val().max(0.0);
                }
                values.push(value);
            }
            values.sort_by(f64::total_cmp);
            let index = (q * values.len().saturating_sub(1) as f64).round() as usize;
            values
                .get(index.min(values.len().saturating_sub(1)))
                .copied()
                .unwrap_or(0.0)
        })
        .collect()
}
#[allow(clippy::too_many_arguments)]
pub fn compute_xva<T: Real>(
    cube: &NpvCube<T>,
    curve: &impl YieldTermStructure<T>,
    csa: Option<&CsaSpec>,
    hc: T,
    lc: T,
    hs: T,
    ls: T,
    fs: T,
) -> XvaReport<T> {
    let mut ee = vec![];
    let mut ene = vec![];
    let mut pfe = vec![];
    let mut cva = T::zero();
    let mut dva = T::zero();
    let mut fva = T::zero();
    for k in 0..cube.times.len() {
        // Keep positive and negative exposures separate all the way through
        // aggregation.  Summing signed MTMs first would illegally net two
        // distinct legal netting sets (e.g. +10 and -8 -> +2 instead of +10).
        let mut pos_vals = Vec::with_capacity(cube.n_paths);
        let mut neg_vals = Vec::with_capacity(cube.n_paths);
        for p in 0..cube.n_paths {
            let mut pos = T::zero();
            let mut neg = T::zero();
            for ns in 0..cube.n_netting_sets {
                let v = cube.get(ns, k, p);
                match csa {
                    Some(c) => {
                        pos += collateralized_exposure(v, c.threshold_cpty, c.mta);
                        neg -= collateralized_exposure(-v, c.threshold_self, c.mta);
                    }
                    None => {
                        pos += v.max_val(T::zero());
                        neg += v.min_val(T::zero());
                    }
                }
            }
            pos_vals.push(pos);
            neg_vals.push(neg);
        }
        let e = pos_vals.iter().copied().fold(T::zero(), |a, b| a + b)
            / T::from_f64(cube.n_paths.max(1) as f64);
        let n = neg_vals.iter().copied().fold(T::zero(), |a, b| a + b)
            / T::from_f64(cube.n_paths.max(1) as f64);
        let mut raw = pos_vals.iter().map(|x| x.val()).collect::<Vec<_>>();
        raw.sort_by(|a, b| a.total_cmp(b));
        let idx = ((0.95 * (raw.len().saturating_sub(1) as f64)).round() as usize)
            .min(raw.len().saturating_sub(1));
        ee.push(e);
        ene.push(n);
        pfe.push(raw.get(idx).copied().unwrap_or(0.));
        let prev = if k == 0 { 0. } else { cube.times[k - 1] };
        let t = cube.times[k];
        let df = curve.discount_factor(t);
        let sc0 = (-hc * prev).exp();
        let sc = (-hc * t).exp();
        let sb0 = (-hs * prev).exp();
        let sb = (-hs * t).exp();
        cva += lc * e * df * (sc0 - sc) * sb0;
        dva -= ls * n * df * (sb0 - sb) * sc0;
        fva += fs * (e + n) * df * sc * sb * T::from_f64(t - prev);
    }
    // EPE is a time integral, not an unweighted arithmetic mean.  Include
    // EE(0)=0 and use trapezoidal integration on the supplied irregular grid.
    let epe = if let Some(&t_last) = cube.times.last() {
        if t_last > 0.0 && ee.len() == cube.times.len() {
            let mut integral = T::zero();
            let mut prev_t = 0.0;
            let mut prev_ee = T::zero();
            for (&t, &e_t) in cube.times.iter().zip(&ee) {
                let dt = (t - prev_t).max(0.0);
                integral += (prev_ee + e_t) * (0.5 * dt);
                prev_t = t;
                prev_ee = e_t;
            }
            integral / T::from_f64(t_last)
        } else {
            T::zero()
        }
    } else {
        T::zero()
    };
    XvaReport {
        ee,
        ene,
        pfe_95: pfe,
        epe,
        cva,
        dva,
        fva,
    }
}

#[allow(clippy::too_many_arguments)]
pub fn compute_kva<T: Real>(
    cube: &NpvCube<T>,
    curve: &impl YieldTermStructure<T>,
    csa: Option<&CsaSpec>,
    hc: T,
    hs: T,
    cpty_risk_weight: T,
    capital_ratio: T,
    hurdle_rate: T,
    add_on_profile: Option<&[T]>,
) -> T {
    let mut kva = T::zero();
    for k in 0..cube.times.len() {
        let t = cube.times[k];
        let prev = if k == 0 { 0.0 } else { cube.times[k - 1] };
        let dt = (t - prev).max(0.0);
        let mut ee = T::zero();
        for p in 0..cube.n_paths {
            let mut pos = T::zero();
            for ns in 0..cube.n_netting_sets {
                let x = cube.get(ns, k, p);
                pos += if let Some(c) = csa {
                    collateralized_exposure(x, c.threshold_cpty, c.mta)
                } else {
                    x.max_val(T::zero())
                };
            }
            ee += pos;
        }
        ee /= T::from_f64(cube.n_paths.max(1) as f64);
        let addon = add_on_profile
            .and_then(|x| x.get(k).copied())
            .unwrap_or(T::zero());
        let ead = T::from_f64(1.4) * (ee + addon);
        let capital = T::from_f64(12.5) * cpty_risk_weight * capital_ratio * ead;
        kva += hurdle_rate * curve.discount_factor(t) * (-(hc + hs) * t).exp() * capital * dt;
    }
    kva
}

#[allow(clippy::too_many_arguments)]
pub fn compute_xva_with_wwr<T: Real>(
    cube: &NpvCube<T>,
    curve: &impl YieldTermStructure<T>,
    csa: Option<&CsaSpec>,
    hc: T,
    lc: T,
    hs: T,
    ls: T,
    fs: T,
    wwr_sensitivity: T,
) -> XvaReport<T> {
    let mut ee = Vec::new();
    let mut ene = Vec::new();
    let mut pfe = Vec::new();
    let mut sc_prev = vec![T::one(); cube.n_paths];
    let mut cva = T::zero();
    let mut dva = T::zero();
    let mut fva = T::zero();
    for k in 0..cube.times.len() {
        let t = cube.times[k];
        let prev = if k == 0 { 0.0 } else { cube.times[k - 1] };
        let dt = (t - prev).max(0.0);
        let mut pos = Vec::with_capacity(cube.n_paths);
        let mut neg = Vec::with_capacity(cube.n_paths);
        for p in 0..cube.n_paths {
            let mut po = T::zero();
            let mut ne = T::zero();
            for ns in 0..cube.n_netting_sets {
                let x = cube.get(ns, k, p);
                if let Some(c) = csa {
                    po += collateralized_exposure(x, c.threshold_cpty, c.mta);
                    ne -= collateralized_exposure(-x, c.threshold_self, c.mta);
                } else {
                    po += x.max_val(T::zero());
                    ne += x.min_val(T::zero());
                }
            }
            pos.push(po);
            neg.push(ne);
        }
        let mean = pos.iter().copied().fold(T::zero(), |a, x| a + x)
            / T::from_f64(cube.n_paths.max(1) as f64);
        let mut beta = Vec::with_capacity(cube.n_paths);
        for x in &pos {
            let z = (wwr_sensitivity * (*x - mean))
                .max_val(T::from_f64(-20.0))
                .min_val(T::from_f64(20.0));
            beta.push(z.exp() * dt);
        }
        let den = T::from_f64(cube.n_paths.max(1) as f64);
        let mut sc = Vec::with_capacity(cube.n_paths);
        if dt <= 0.0 {
            sc.extend_from_slice(&sc_prev);
        } else {
            let target = (-hc * t).exp();
            let mean_sc = sc_prev
                .iter()
                .zip(&beta)
                .fold(T::zero(), |a, (s, b)| a + *s * (*b / T::from_f64(dt)))
                / den;
            let mut u = (hc * (-hc * prev).exp() / mean_sc.max_val(T::from_f64(1e-12)))
                .max_val(T::from_f64(1e-12));
            for _ in 0..8 {
                let mut g = T::zero();
                let mut gp = T::zero();
                for p in 0..cube.n_paths {
                    let e = (-u * beta[p]).exp();
                    g += sc_prev[p] * e;
                    gp -= sc_prev[p] * beta[p] * e;
                }
                g = g / den - target;
                gp /= den;
                if gp.val().abs() > 1e-14 {
                    u = (u - g / gp).max_val(T::from_f64(1e-12));
                }
            }
            for p in 0..cube.n_paths {
                sc.push(sc_prev[p] * (-u * beta[p]).exp());
            }
        }
        let e = pos.iter().copied().fold(T::zero(), |a, x| a + x)
            / T::from_f64(cube.n_paths.max(1) as f64);
        let n = neg.iter().copied().fold(T::zero(), |a, x| a + x)
            / T::from_f64(cube.n_paths.max(1) as f64);
        let mut raw = pos.iter().map(|x| x.val()).collect::<Vec<_>>();
        raw.sort_by(f64::total_cmp);
        let idx = (0.95 * (raw.len().saturating_sub(1) as f64)).round() as usize;
        ee.push(e);
        ene.push(n);
        pfe.push(raw.get(idx).copied().unwrap_or(0.));
        let df = curve.discount_factor(t);
        let sb0 = (-hs * prev).exp();
        let sb = (-hs * t).exp();
        let mut cva_step = T::zero();
        let mut dva_step = T::zero();
        let mut fva_step = T::zero();
        for p in 0..cube.n_paths {
            cva_step += pos[p] * (sc_prev[p] - sc[p]);
            dva_step += neg[p] * (sb0 - sb) * sc_prev[p];
            fva_step += (pos[p] + neg[p]) * sc[p];
        }
        cva += lc * df * sb0 * cva_step / den;
        dva -= ls * df * dva_step / den;
        fva += fs * df * sb * fva_step / den * dt;
        sc_prev = sc;
    }
    let epe = if let Some(&last) = cube.times.last() {
        if last > 0.0 {
            let mut integral = T::zero();
            let mut prev_t = 0.0;
            let mut prev_e = T::zero();
            for (t, e_t) in cube.times.iter().copied().zip(ee.iter().copied()) {
                let dt_t = (t - prev_t).max(0.0);
                integral += (prev_e + e_t) * T::from_f64(0.5 * dt_t);
                prev_t = t;
                prev_e = e_t;
            }
            integral / T::from_f64(last)
        } else {
            T::zero()
        }
    } else {
        T::zero()
    };
    XvaReport {
        ee,
        ene,
        pfe_95: pfe,
        epe,
        cva,
        dva,
        fva,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CsaTerms<T: Real> {
    pub collateral: T,
    pub threshold: T,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum XvaCloseOutConvention {
    RiskFree,
    Replacement,
}
#[derive(Clone, Debug)]
pub struct BilateralFundingSpec<T: Real> {
    pub lgd_counterparty: T,
    pub lgd_own: T,
    pub funding_borrow_spread: T,
    pub funding_lend_spread: T,
    pub collateral_pay_spread: T,
    pub collateral_rec_spread: T,
    pub close_out: XvaCloseOutConvention,
}
#[derive(Clone, Debug)]
pub struct BilateralXvaReport<T: Real> {
    pub cva: T,
    pub dva: T,
    pub fca: T,
    pub fba: T,
    pub colva: T,
    pub total_xva: T,
    pub fixed_point_iterations: usize,
}
pub fn compute_bilateral_xva<T: Real>(
    cube: &NpvCube<T>,
    csa: Option<&CsaTerms<T>>,
    ois_dfs: &[T],
    cpty_survival: &[T],
    own_survival: &[T],
    spec: &BilateralFundingSpec<T>,
) -> BilateralXvaReport<T> {
    let mut cva = T::zero();
    let mut dva = T::zero();
    let mut fca = T::zero();
    let mut fba = T::zero();
    let mut colva = T::zero();
    let mut u = T::zero();
    let mut iterations = 0;
    for k in 0..cube.times.len() {
        let t = cube.times[k];
        let prev = if k == 0 { 0. } else { cube.times[k - 1] };
        let dt = t - prev;
        if dt <= 0. {
            continue;
        }
        let df = ois_dfs.get(k).copied().unwrap_or(T::one());
        let sc0 = cpty_survival
            .get(k.saturating_sub(1))
            .copied()
            .unwrap_or(T::one());
        let sc = cpty_survival.get(k).copied().unwrap_or(sc0);
        let sb0 = own_survival
            .get(k.saturating_sub(1))
            .copied()
            .unwrap_or(T::one());
        let sb = own_survival.get(k).copied().unwrap_or(sb0);
        let den = T::from_f64(cube.n_paths.max(1) as f64);
        let mut pos = T::zero();
        let mut neg = T::zero();
        let mut xp = T::zero();
        let mut xn = T::zero();
        for p in 0..cube.n_paths {
            let v = (0..cube.n_netting_sets)
                .map(|ns| cube.get(ns, k, p))
                .fold(T::zero(), |a, z| a + z);
            let x = if let Some(c) = csa {
                c.collateral + (v - c.collateral - c.threshold).max_val(T::zero())
                    - (-v + c.collateral - c.threshold).max_val(T::zero())
            } else {
                T::zero()
            };
            let e = v - u - x;
            pos += e.max_val(T::zero());
            neg += (-e).max_val(T::zero());
            xp += x.max_val(T::zero());
            xn += (-x).max_val(T::zero());
        }
        pos /= den;
        neg /= den;
        xp /= den;
        xn /= den;
        let cv = spec.lgd_counterparty * df * sb0 * (sc0 - sc) * pos;
        let dv = spec.lgd_own * df * sc0 * (sb0 - sb) * neg;
        let fc = spec.funding_borrow_spread * df * sc * sb * pos * T::from_f64(dt);
        let fb = spec.funding_lend_spread * df * sc * sb * neg * T::from_f64(dt);
        let co = (spec.collateral_pay_spread * xp + spec.collateral_rec_spread * xn)
            * df
            * sc
            * sb
            * T::from_f64(dt);
        if matches!(spec.close_out, XvaCloseOutConvention::Replacement) {
            let next = -(cv - dv + fc - fb + co);
            iterations += 1;
            if (next - u).val().abs() < 1e-13 {
                u = next;
                break;
            }
            u = next;
        }
        cva += cv;
        dva += dv;
        fca += fc;
        fba += fb;
        colva += co;
    }
    if matches!(spec.close_out, XvaCloseOutConvention::Replacement) {
        cva = -u;
    }
    let total = cva - dva + fca - fba + colva;
    BilateralXvaReport {
        cva,
        dva,
        fca,
        fba,
        colva,
        total_xva: total,
        fixed_point_iterations: iterations,
    }
}

// Phase 12 Round 4 audit coverage.

// Phase 12 Round 4 surgical audit coverage.
