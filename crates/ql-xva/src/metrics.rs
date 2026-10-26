#![allow(clippy::too_many_arguments)]
use crate::NpvCube;
pub use crate::phase35::*;
pub use crate::phase37::*;
use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct AlbaneseCrepeyMporBalanceSheetSpec<T: Real = f64> {
    pub exposure_times: Vec<T>,
    pub clean_mtm_profile: Vec<T>,
    pub portfolio_vol_profile: Vec<T>,
    pub scheduled_cashflows_in_mpor: Vec<T>,
    pub mpor_years: T,
    pub csa_threshold: T,
    pub risk_free_rate: T,
    pub cpty_hazard_rate: T,
    pub cpty_lgd: T,
    pub bank_hazard_rate: T,
    pub bank_lgd: T,
    pub economic_capital_profile: Vec<T>,
    pub hurdle_rate: T,
    pub initial_margin_profile: Vec<T>,
    pub im_funding_spread: T,
}
#[derive(Clone, Debug)]
pub struct AlbaneseCrepeyMporBalanceSheetReport<T: Real = f64> {
    pub classical_mpor_epe_profile: Vec<T>,
    pub sticky_mpor_epe_profile: Vec<T>,
    pub classical_mpor_cva: T,
    pub sticky_mpor_cva: T,
    pub mpor_coupon_spike_cva_amplification: T,
    pub unadjusted_fca_without_capital_offset: T,
    pub albanese_crepey_shareholder_fca: T,
    pub creditor_fda_wealth_transfer: T,
    pub mva: T,
    pub kva: T,
    pub shareholder_funds_transfer_price_ftp: T,
}
pub fn compute_albanese_crepey_mpor_and_balance_sheet_xva<T: Real>(
    s: &AlbaneseCrepeyMporBalanceSheetSpec<T>,
) -> AlbaneseCrepeyMporBalanceSheetReport<T> {
    let n = s.exposure_times.len();
    let mut cls = Vec::with_capacity(n);
    let mut sty = Vec::with_capacity(n);
    let smooth_pos = |x: T| T::from_f64(0.5) * (x + (x * x + T::from_f64(1e-12)).sqrt());
    for i in 0..n {
        let sig = (s
            .portfolio_vol_profile
            .get(i)
            .copied()
            .unwrap_or(T::zero())
            .abs()
            * s.mpor_years.max_val(T::from_f64(1e-8)).sqrt())
        .max_val(T::from_f64(1e-8));
        let mtm = s.clean_mtm_profile.get(i).copied().unwrap_or_else(T::zero);
        let mu_c = s.csa_threshold + T::from_f64(0.1) * (mtm - s.csa_threshold);
        let zc = mu_c / sig;
        let ec = mu_c * zc.norm_cdf() + sig * zc.norm_pdf();
        cls.push(ec);
        let cf = s
            .scheduled_cashflows_in_mpor
            .get(i)
            .copied()
            .unwrap_or_else(T::zero);
        let mu_s = mu_c + smooth_pos(cf);
        let zs = mu_s / sig;
        let es = mu_s * zs.norm_cdf() + sig * zs.norm_pdf();
        sty.push(es);
    }
    let mut cva = T::zero();
    let mut scva = T::zero();
    let mut fca_unadj = T::zero();
    let mut fca_sh = T::zero();
    let mut kva = T::zero();
    let mut mva = T::zero();
    let sf = s.bank_lgd * s.bank_hazard_rate;
    let mut prev_t = T::zero();
    for i in 0..n {
        let ti = s.exposure_times[i];
        let dt = (ti - prev_t).max_val(T::from_f64(1e-6));
        prev_t = ti;
        let df = (-s.risk_free_rate * ti).exp() * dt;
        cva += s.cpty_lgd * s.cpty_hazard_rate * df * cls[i];
        scva += s.cpty_lgd * s.cpty_hazard_rate * df * sty[i];
        let cap = s
            .economic_capital_profile
            .get(i)
            .copied()
            .unwrap_or_else(T::zero);
        let im = s
            .initial_margin_profile
            .get(i)
            .copied()
            .unwrap_or_else(T::zero);
        fca_unadj += sf * df * sty[i];
        fca_sh += sf * df * smooth_pos(sty[i] - cap);
        kva += s.hurdle_rate * df * cap;
        mva += s.im_funding_spread * df * im;
    }
    let fda = (fca_unadj - fca_sh).max_val(T::zero());
    AlbaneseCrepeyMporBalanceSheetReport {
        classical_mpor_epe_profile: cls,
        sticky_mpor_epe_profile: sty,
        classical_mpor_cva: cva,
        sticky_mpor_cva: scva,
        mpor_coupon_spike_cva_amplification: scva / cva.max_val(T::from_f64(1e-12)) - T::one(),
        unadjusted_fca_without_capital_offset: fca_unadj,
        albanese_crepey_shareholder_fca: fca_sh,
        creditor_fda_wealth_transfer: fda,
        mva,
        kva,
        shareholder_funds_transfer_price_ftp: scva + fca_sh + mva + kva,
    }
}

use ql_curves::YieldTermStructure;

#[derive(Clone, Debug)]
pub struct BrigoPallaviciniCtdCollateralXvaSpec<T: Real = f64> {
    pub notional_mtm: T,
    pub maturity: T,
    pub num_steps: usize,
    pub base_ois_rate: T,
    pub csa_spread_ccy1: T,
    pub csa_spread_ccy2: T,
    pub csa_spread_ccy3: T,
    pub csa_mean_reversion: T,
    pub csa_vol_ccy1: T,
    pub csa_vol_ccy2: T,
    pub csa_vol_ccy3: T,
    pub csa_correlation: T,
    pub cpty_hazard_rate: T,
    pub cpty_lgd: T,
    pub own_hazard_rate: T,
    pub own_lgd: T,
    pub exposure_volatility: T,
}
#[derive(Clone, Debug)]
pub struct BrigoPallaviciniCtdCollateralXvaReport<T: Real = f64> {
    pub single_ccy1_discount_factor: T,
    pub intrinsic_max_csa_discount_factor: T,
    pub ctd_stochastic_csa_discount_factor: T,
    pub ctd_collateral_switch_option_value: T,
    pub effective_ctd_spread_pickup_bps: T,
    pub risk_free_closeout_bcva: T,
    pub replacement_cost_closeout_bcva: T,
    pub closeout_contagion_spread_impact: T,
}
fn clark_max2<T: Real>(mu1: T, v1_sq: T, mu2: T, v2_sq: T, rho12: T) -> (T, T, T) {
    let v1 = v1_sq.max_val(T::from_f64(1e-16)).sqrt();
    let v2 = v2_sq.max_val(T::from_f64(1e-16)).sqrt();
    let a_sq = (v1_sq + v2_sq - T::from_f64(2.0) * rho12 * v1 * v2).max_val(T::from_f64(1e-14));
    let a = a_sq.sqrt();
    let d = (mu1 - mu2) / a;
    let p1 = d.norm_cdf();
    let p2 = (-d).norm_cdf();
    let pdf = d.norm_pdf();
    let mu = mu1 * p1 + mu2 * p2 + a * pdf;
    let second = (mu1 * mu1 + v1_sq) * p1 + (mu2 * mu2 + v2_sq) * p2 + (mu1 + mu2) * a * pdf;
    ((mu), (second - mu * mu).max_val(T::from_f64(1e-14)), p1)
}
pub fn compute_brigo_pallavicini_ctd_collateral_and_closeout_xva<T: Real>(
    s: &BrigoPallaviciniCtdCollateralXvaSpec<T>,
) -> BrigoPallaviciniCtdCollateralXvaReport<T> {
    let n = s.num_steps.max(1);
    let dt = s.maturity / T::from_f64(n as f64);
    let mut integ = T::zero();
    let maxs = s
        .csa_spread_ccy1
        .max_val(s.csa_spread_ccy2)
        .max_val(s.csa_spread_ccy3);
    let mut mus = Vec::with_capacity(n);
    for k in 0..n {
        let t = T::from_f64(k as f64 + 0.5) * dt;
        let g = (T::one() - (-T::from_f64(2.0) * s.csa_mean_reversion * t).exp())
            / (T::from_f64(2.0) * s.csa_mean_reversion.max_val(T::from_f64(1e-12)));
        let v1 = s.csa_vol_ccy1 * s.csa_vol_ccy1 * g;
        let v2 = s.csa_vol_ccy2 * s.csa_vol_ccy2 * g;
        let v3 = s.csa_vol_ccy3 * s.csa_vol_ccy3 * g;
        let (mu12, v12, p1) = clark_max2(
            s.csa_spread_ccy1,
            v1,
            s.csa_spread_ccy2,
            v2,
            s.csa_correlation,
        );
        let v12s = v12.max_val(T::from_f64(1e-16)).sqrt();
        let rho123 = ((s.csa_vol_ccy1 * g.sqrt() * s.csa_correlation * p1
            + s.csa_vol_ccy2 * g.sqrt() * s.csa_correlation * (T::one() - p1))
            / v12s)
            .max_val(T::from_f64(-0.999))
            .min_val(T::from_f64(0.999));
        let (mu123_raw, _, _) = clark_max2(mu12, v12, s.csa_spread_ccy3, v3, rho123);
        let mu123 = if s.csa_correlation.val() > 0.999
            && (s.csa_spread_ccy1 - s.csa_spread_ccy2).abs().val() < 1e-12
            && (s.csa_spread_ccy2 - s.csa_spread_ccy3).abs().val() < 1e-12
            && (s.csa_vol_ccy1 - s.csa_vol_ccy2).abs().val() < 1e-12
            && (s.csa_vol_ccy2 - s.csa_vol_ccy3).abs().val() < 1e-12
        {
            maxs
        } else {
            mu123_raw
        };
        mus.push(mu123);
        integ += mu123 * dt;
    }
    let d1 = (-s.base_ois_rate * s.maturity).exp();
    let intrinsic = (-(s.base_ois_rate + maxs) * s.maturity).exp();
    let stochastic = (-s.base_ois_rate * s.maturity - integ).exp();
    let ctd = (intrinsic - stochastic).max_val(T::zero());
    let pickup = T::from_f64(1e4) * (intrinsic.ln() - stochastic.ln())
        / s.maturity.max_val(T::from_f64(1e-12));
    let mut cva = T::zero();
    let mut theta_next = T::zero();
    let mut theta = T::zero();
    for kk in (0..n).rev() {
        let t = T::from_f64(kk as f64 + 0.5) * dt;
        let mu123 = mus[kk];
        let disc = (-(s.base_ois_rate + mu123) * t).exp();
        let surv = (-(s.cpty_hazard_rate + s.own_hazard_rate) * t).exp();
        let sig = (s.exposure_volatility * t.sqrt()).max_val(T::from_f64(1e-6));
        let m = s.notional_mtm;
        let z = m / sig;
        let ep = m * z.norm_cdf() + sig * z.norm_pdf();
        let en = ep - m;
        cva += disc
            * surv
            * (s.cpty_hazard_rate * s.cpty_lgd * ep - s.own_hazard_rate * s.own_lgd * en)
            * dt;
        let mut th = theta_next;
        for _ in 0..6 {
            let mm = s.notional_mtm - th;
            let zz = mm / sig;
            let epp = mm * zz.norm_cdf() + sig * zz.norm_pdf();
            let enn = epp - mm;
            th = theta_next
                * (-(s.base_ois_rate + s.cpty_hazard_rate + s.own_hazard_rate) * dt).exp()
                + (s.cpty_hazard_rate * s.cpty_lgd * epp - s.own_hazard_rate * s.own_lgd * enn)
                    * dt;
        }
        theta = th;
        theta_next = th;
    }
    let rep = theta;
    BrigoPallaviciniCtdCollateralXvaReport {
        single_ccy1_discount_factor: d1,
        intrinsic_max_csa_discount_factor: intrinsic,
        ctd_stochastic_csa_discount_factor: stochastic,
        ctd_collateral_switch_option_value: ctd,
        effective_ctd_spread_pickup_bps: pickup,
        risk_free_closeout_bcva: cva,
        replacement_cost_closeout_bcva: rep,
        closeout_contagion_spread_impact: (rep - cva).abs(),
    }
}

#[derive(Clone, Debug)]
pub struct GregoryNettingSetAllocationSpec<T: Real = f64> {
    pub trade_weights: Vec<T>,
    pub trade_expected_mtms: Vec<T>,
    pub trade_exposure_vols: Vec<T>,
    pub correlation_matrix: Vec<Vec<T>>,
    pub collateral_fraction: T,
    pub horizon: f64,
    pub num_steps: usize,
    pub risk_free_rate: T,
    pub cpty_hazard_rate: T,
    pub cpty_lgd: T,
    pub own_hazard_rate: T,
    pub own_lgd: T,
    pub funding_spread: T,
    pub colva_spread: T,
}
#[derive(Clone, Debug)]
pub struct GregoryNettingSetAllocationReport<T: Real = f64> {
    pub netted_total_xva: T,
    pub standalone_total_xva: T,
    pub netting_benefit: T,
    pub netting_factor: T,
    pub euler_allocated_cva: Vec<T>,
    pub euler_allocated_dva: Vec<T>,
    pub euler_allocated_fva: Vec<T>,
    pub euler_allocated_total_xva: Vec<T>,
    pub incremental_xva_by_trade: Vec<T>,
}
pub fn allocate_gregory_netting_set_xva<T: Real>(
    s: &GregoryNettingSetAllocationSpec<T>,
) -> GregoryNettingSetAllocationReport<T> {
    let n = s.trade_weights.len();
    let a = (T::one() - s.collateral_fraction).max_val(T::from_f64(1e-12));
    let mut mu = T::zero();
    let mut var = T::zero();
    for i in 0..n {
        mu += s.trade_weights[i] * s.trade_expected_mtms[i];
        for j in 0..n {
            var += s.trade_weights[i]
                * s.trade_weights[j]
                * s.correlation_matrix
                    .get(i)
                    .and_then(|r| r.get(j))
                    .copied()
                    .unwrap_or(T::zero())
                * s.trade_exposure_vols[i]
                * s.trade_exposure_vols[j];
        }
    }
    let mut cva = T::zero();
    let mut dva = T::zero();
    let mut fva = T::zero();
    let mut ac = vec![T::zero(); n];
    let mut ad = vec![T::zero(); n];
    let mut af = vec![T::zero(); n];
    let steps = s.num_steps.max(1);
    for k in 1..=steps {
        let t = s.horizon * k as f64 / steps as f64;
        let sk = (t * (1.0 - 0.5 * t / s.horizon.max(1e-8))).max(0.0).sqrt();
        let m = a * mu * sk;
        let sig = (a * a * var).max_val(T::from_f64(1e-14)).sqrt() * T::from_f64(sk);
        let d = m / sig.max_val(T::from_f64(1e-10));
        let epe = m * d.norm_cdf() + sig * d.norm_pdf();
        let ene = -m * (-d).norm_cdf() + sig * d.norm_pdf();
        let wt = T::from_f64(s.horizon / steps as f64)
            * (-(s.risk_free_rate + s.cpty_hazard_rate + s.own_hazard_rate) * t).exp();
        let wc = wt * s.cpty_hazard_rate * s.cpty_lgd;
        let wd = -wt * s.own_hazard_rate * s.own_lgd;
        let wf = wt * s.funding_spread + wt * s.colva_spread * s.collateral_fraction / a;
        cva += wc * epe;
        dva += wd * ene;
        fva += wf * (epe + ene);
        for i in 0..n {
            let cross = (0..n)
                .map(|j| {
                    s.trade_weights[j]
                        * s.correlation_matrix
                            .get(i)
                            .and_then(|r| r.get(j))
                            .copied()
                            .unwrap_or(T::zero())
                        * s.trade_exposure_vols[i]
                        * s.trade_exposure_vols[j]
                })
                .fold(T::zero(), |x, y| x + y);
            let common = s.trade_weights[i] * a * s.trade_expected_mtms[i] * sk * d.norm_cdf();
            let risk = s.trade_weights[i] * a * a * cross / sig.max_val(T::from_f64(1e-10))
                * T::from_f64(sk * sk)
                * d.norm_pdf();
            ac[i] += wc * (common + risk);
            ad[i] += wd
                * (-s.trade_weights[i] * a * s.trade_expected_mtms[i] * sk * (-d).norm_cdf()
                    + risk);
            af[i] += wf
                * (common
                    + risk
                    + (-s.trade_weights[i] * a * s.trade_expected_mtms[i] * sk * (-d).norm_cdf()
                        + risk));
        }
    }
    let total = cva + dva + fva;
    let stand = s
        .trade_weights
        .iter()
        .enumerate()
        .map(|(i, _w)| {
            let mut one = s.clone();
            for j in 0..n {
                one.trade_weights[j] = if i == j {
                    s.trade_weights[j]
                } else {
                    T::zero()
                };
            }
            allocate_gregory_netting_set_xva_single(&one)
        })
        .fold(T::zero(), |x, y| x + y);
    let alloc_total: Vec<T> = (0..n).map(|i| ac[i] + ad[i] + af[i]).collect();
    let mut inc = Vec::with_capacity(n);
    for i in 0..n {
        let mut minus = s.clone();
        minus.trade_weights[i] = T::zero();
        inc.push(total - allocate_gregory_netting_set_xva_core(&minus));
    }
    GregoryNettingSetAllocationReport {
        netted_total_xva: total,
        standalone_total_xva: stand,
        netting_benefit: (stand - total).max_val(T::zero()),
        netting_factor: total / stand.max_val(T::from_f64(1e-12)),
        euler_allocated_cva: ac,
        euler_allocated_dva: ad,
        euler_allocated_fva: af,
        euler_allocated_total_xva: alloc_total,
        incremental_xva_by_trade: inc,
    }
}
fn allocate_gregory_netting_set_xva_single<T: Real>(s: &GregoryNettingSetAllocationSpec<T>) -> T {
    let x = s.clone();
    let n = x.trade_weights.len();
    let mut total = T::zero();
    for i in 0..n {
        if x.trade_weights[i].val() != 0.0 {
            let mut y = x.clone();
            for (j, weight) in y.trade_weights.iter_mut().enumerate().take(n) {
                *weight = if i == j {
                    x.trade_weights[i]
                } else {
                    T::zero()
                };
            }
            total += allocate_gregory_netting_set_xva_core(&y);
        }
    }
    total
}
fn allocate_gregory_netting_set_xva_core<T: Real>(s: &GregoryNettingSetAllocationSpec<T>) -> T {
    let n = s.trade_weights.len();
    let a = (T::one() - s.collateral_fraction).max_val(T::from_f64(1e-12));
    let mut mu = T::zero();
    let mut var = T::zero();
    for i in 0..n {
        mu += s.trade_weights[i] * s.trade_expected_mtms[i];
        for j in 0..n {
            var += s.trade_weights[i]
                * s.trade_weights[j]
                * s.correlation_matrix[i][j]
                * s.trade_exposure_vols[i]
                * s.trade_exposure_vols[j];
        }
    }
    let mut x = T::zero();
    for k in 1..=s.num_steps.max(1) {
        let t = s.horizon * k as f64 / s.num_steps.max(1) as f64;
        let sk = (t * (1. - 0.5 * t / s.horizon.max(1e-8))).sqrt();
        let m = a * mu * sk;
        let sig = (a * a * var).sqrt() * T::from_f64(sk);
        let d = m / sig.max_val(T::from_f64(1e-10));
        let e = m * d.norm_cdf() + sig * d.norm_pdf();
        let en = -m * (-d).norm_cdf() + sig * d.norm_pdf();
        let wt = T::from_f64(s.horizon / s.num_steps.max(1) as f64)
            * (-(s.risk_free_rate + s.cpty_hazard_rate + s.own_hazard_rate) * t).exp();
        x += wt
            * (s.cpty_hazard_rate * s.cpty_lgd * e - s.own_hazard_rate * s.own_lgd * en
                + (s.funding_spread + s.colva_spread * s.collateral_fraction / a) * (e + en));
    }
    x
}
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
pub fn compute_correlated_first_to_default_bilateral_xva<T: Real>(
    cube: &NpvCube<T>,
    curve: &impl YieldTermStructure<T>,
    csa: Option<&CsaSpec>,
    hc: T,
    lc: T,
    hs: T,
    ls: T,
    fs: T,
    rho_cb: T,
) -> XvaReport<T> {
    let mut ee = Vec::new();
    let mut ene = Vec::new();
    let mut pfe = Vec::new();
    let mut cva = T::zero();
    let mut dva = T::zero();
    let mut fva = T::zero();
    let inv_paths = T::from_f64(1.0 / cube.n_paths.max(1) as f64);
    let joint = |tc: T, tb: T| {
        if tc.val() <= 0.0 || tb.val() <= 0.0 {
            return T::zero();
        }
        if tc.val() >= 1.0 - 1e-14 {
            return tb;
        }
        if tb.val() >= 1.0 - 1e-14 {
            return tc;
        }
        let ic = T::from_f64(ql_math::inv_norm_cdf(tc.val().clamp(1e-14, 1.0 - 1e-14)));
        let ib = T::from_f64(ql_math::inv_norm_cdf(tb.val().clamp(1e-14, 1.0 - 1e-14)));
        let mut xc = ic;
        let mut xb = ib;
        for _ in 0..2 {
            xc -= (xc.norm_cdf() - tc) / xc.norm_pdf().max_val(T::from_f64(1e-14));
            xb -= (xb.norm_cdf() - tb) / xb.norm_pdf().max_val(T::from_f64(1e-14));
        }
        ql_math::bivariate_normal_cdf(
            xc,
            xb,
            rho_cb
                .max_val(T::from_f64(-0.999))
                .min_val(T::from_f64(0.999)),
        )
    };
    for k in 0..cube.times.len() {
        let t = cube.times[k];
        let prev = if k == 0 { 0.0 } else { cube.times[k - 1] };
        let dt = (t - prev).max(0.0);
        let mut pos = T::zero();
        let mut neg = T::zero();
        let mut vals = Vec::with_capacity(cube.n_paths);
        for p in 0..cube.n_paths {
            let mut ep = T::zero();
            let mut en = T::zero();
            for ns in 0..cube.n_netting_sets {
                let x = cube.get(ns, k, p);
                ep += if let Some(c) = csa {
                    collateralized_exposure(x, c.threshold_cpty, c.mta)
                } else {
                    x.max_val(T::zero())
                };
                en += if let Some(c) = csa {
                    -collateralized_exposure(-x, c.threshold_self, c.mta)
                } else {
                    x.min_val(T::zero())
                };
            }
            pos += ep;
            neg += en;
            vals.push(ep.val());
        }
        let ep = pos * inv_paths;
        let en = neg * inv_paths;
        vals.sort_by(f64::total_cmp);
        let ix = (0.95 * vals.len().saturating_sub(1) as f64).round() as usize;
        ee.push(ep);
        ene.push(en);
        pfe.push(vals.get(ix).copied().unwrap_or(0.0));
        let sc0 = (-hc * T::from_f64(prev)).exp();
        let sc1 = (-hc * T::from_f64(t)).exp();
        let sb0 = (-hs * T::from_f64(prev)).exp();
        let sb1 = (-hs * T::from_f64(t)).exp();
        let j00 = joint(sc0, sb0);
        let j10 = joint(sc1, sb0);
        let j01 = joint(sc0, sb1);
        let j11 = joint(sc1, sb1);
        let c_inc = ((j00 - j10) + (j01 - j11)).max_val(T::zero()) * 0.5;
        let b_inc = ((j00 - j01) + (j10 - j11)).max_val(T::zero()) * 0.5;
        let df = curve.discount_factor(t);
        cva += lc * ep * df * c_inc;
        dva -= ls * en * df * b_inc;
        fva += fs * (ep + en) * df * j11 * T::from_f64(dt);
    }
    let epe = if let Some(&last) = cube.times.last() {
        let mut acc = T::zero();
        let mut pt = 0.0;
        let mut pe = T::zero();
        for (&t, &x) in cube.times.iter().zip(&ee) {
            acc += (pe + x) * (t - pt).max(0.0) * 0.5;
            pt = t;
            pe = x;
        }
        if last > 0.0 {
            acc / T::from_f64(last)
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

#[derive(Clone, Debug)]
pub struct MultiCurrencyCsaSpec<T: Real = f64> {
    pub eligible_currency_spreads: Vec<T>,
    pub spread_vols: Vec<T>,
    pub cross_currency_corr: T,
    pub collateral_balance_profile: Vec<T>,
    pub payment_times: Vec<f64>,
}

#[derive(Copy, Clone, Debug)]
pub struct MultiCurrencyColvaReport<T: Real> {
    pub intrinsic_ctd_colva: T,
    pub switch_option_colva: T,
    pub total_ctd_colva: T,
    pub effective_ctd_spread: T,
}

pub fn compute_multicurrency_ctd_colva<T: Real>(
    s: &MultiCurrencyCsaSpec<T>,
) -> MultiCurrencyColvaReport<T> {
    let n = s.payment_times.len();
    let m_ccy = s.eligible_currency_spreads.len();
    if m_ccy == 0 || n == 0 {
        return MultiCurrencyColvaReport {
            intrinsic_ctd_colva: T::zero(),
            switch_option_colva: T::zero(),
            total_ctd_colva: T::zero(),
            effective_ctd_spread: T::zero(),
        };
    }
    let rho = s
        .cross_currency_corr
        .max_val(T::from_f64(-1.0 / (m_ccy.max(2) as f64 - 1.0) + 1e-4))
        .min_val(T::from_f64(0.999));
    let mut mx = s.eligible_currency_spreads[0];
    for &sp in &s.eligible_currency_spreads[1..] {
        let diff = sp - mx;
        mx = (mx + sp + (diff * diff + T::from_f64(1e-16)).sqrt()) * 0.5;
    }

    let mut intr = T::zero();
    let mut sw = T::zero();
    let mut bal = T::zero();
    let mut prev_t = 0.0;

    for i in 0..n {
        let t_f64 = s.payment_times[i].max(prev_t);
        let dt = (t_f64 - prev_t).max(0.0);
        prev_t = t_f64;
        let t = T::from_f64(t_f64);

        let vols: Vec<T> = (0..m_ccy)
            .map(|j| {
                let v = s.spread_vols.get(j).copied().unwrap_or(T::zero());
                (v * v + T::from_f64(1e-16)).sqrt()
            })
            .collect();
        let mut mu_y = s.eligible_currency_spreads[0];
        let mut var_y = vols[0] * vols[0] * t;
        let mut cov_y: Vec<T> = (0..m_ccy)
            .map(|k| {
                if k == 0 {
                    var_y
                } else {
                    rho * vols[0] * vols[k] * t
                }
            })
            .collect();

        for j in 1..m_ccy {
            let mu_x = s.eligible_currency_spreads[j];
            let var_x = vols[j] * vols[j] * t;
            let cov_yx = cov_y[j];
            let a2 = (var_y + var_x - cov_yx * 2.0).max_val(T::zero());
            let a = (a2 + T::from_f64(1e-16)).sqrt();
            let alpha = (mu_y - mu_x) / a;
            let p = alpha.norm_cdf();
            let q = T::one() - p;
            let phi = alpha.norm_pdf();

            let mu_next = mu_y * p + mu_x * q + a * phi;
            let ey2 =
                (mu_y * mu_y + var_y) * p + (mu_x * mu_x + var_x) * q + (mu_y + mu_x) * a * phi;
            let var_next = (ey2 - mu_next * mu_next).max_val(T::zero());

            for k in (j + 1)..m_ccy {
                let cov_xk = rho * vols[j] * vols[k] * t;
                cov_y[k] = cov_y[k] * p + cov_xk * q;
            }
            mu_y = mu_next;
            var_y = var_next;
        }

        let em = mu_y.max_val(mx);
        let b = s
            .collateral_balance_profile
            .get(i)
            .copied()
            .unwrap_or(T::one());
        intr += mx * b * dt;
        sw += (em - mx) * b * dt;
        bal += b * dt;
    }
    let total = intr + sw;
    let eff = total / bal.max_val(T::from_f64(1e-14));
    MultiCurrencyColvaReport {
        intrinsic_ctd_colva: intr,
        switch_option_colva: sw,
        total_ctd_colva: total,
        effective_ctd_spread: eff,
    }
}

#[derive(Clone, Debug)]
pub struct TurlakovWwrSpec<T: Real = f64> {
    pub times: Vec<f64>,
    pub base_epe: Vec<T>,
    pub stressed_epe: Vec<T>,
    pub base_mtm: Vec<T>,
    pub stressed_mtm: Vec<T>,
    pub sovereign_hazard_rate: T,
    pub counterparty_hazard_rate: T,
    pub crisis_prob_threshold: T,
    pub coupling_lambda: T,
    pub lgd: T,
    pub base_funding_spread: T,
    pub stressed_funding_spread: T,
    pub discount_rate: T,
}

#[derive(Copy, Clone, Debug)]
pub struct TurlakovWwrReport<T: Real> {
    pub wwr_cva: T,
    pub base_cva: T,
    pub wwr_fva: T,
    pub base_fva: T,
}

pub fn compute_turlakov_wwr_cva_fva<T: Real>(s: &TurlakovWwrSpec<T>) -> TurlakovWwrReport<T> {
    let n = s.times.len().min(s.base_epe.len());
    let mut wc = T::zero();
    let mut bc = T::zero();
    let mut wf = T::zero();
    let mut bf = T::zero();
    let mut prev_t = 0.0;
    let p_thresh = s.crisis_prob_threshold.max_val(T::from_f64(1e-12));
    let lam_c = s.coupling_lambda.max_val(T::zero()).min_val(T::one());

    for i in 0..n {
        let t = s.times[i].max(prev_t);
        let dt = T::from_f64((t - prev_t).max(0.0));
        let ps = T::one() - (-s.sovereign_hazard_rate * t).exp();
        let w_normal = T::one() - (ps / p_thresh).tanh();
        let w_crisis = T::one() - w_normal;
        let w_stress = lam_c * w_crisis;

        let base_e = s.base_epe[i];
        let stress_e = s.stressed_epe.get(i).copied().unwrap_or(base_e);
        let eff_epe = base_e * (T::one() - w_stress) + stress_e * w_stress;

        let df = (-s.discount_rate * t).exp();
        let surv_prev = (-s.counterparty_hazard_rate * prev_t).exp();
        let surv_curr = (-s.counterparty_hazard_rate * t).exp();
        let dpd = (surv_prev - surv_curr).max_val(T::zero());

        bc += s.lgd * dpd * df * base_e;
        wc += s.lgd * dpd * df * eff_epe;

        let mt = s.base_mtm.get(i).copied().unwrap_or(T::zero());
        let sm = s.stressed_mtm.get(i).copied().unwrap_or(mt);
        let base_fund = mt * s.base_funding_spread;
        let stress_fund = sm * s.stressed_funding_spread;
        let eff_fund = base_fund * (T::one() - w_stress) + stress_fund * w_stress;

        bf += base_fund * df * dt;
        wf += eff_fund * df * dt;
        prev_t = t;
    }
    TurlakovWwrReport {
        wwr_cva: wc,
        base_cva: bc,
        wwr_fva: wf,
        base_fva: bf,
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
    let mut u_cont = T::zero();
    let mut iterations = 0;
    let den = T::from_f64(cube.n_paths.max(1) as f64);
    let ns_count = T::from_f64(cube.n_netting_sets.max(1) as f64);
    for k in (0..cube.times.len()).rev() {
        let t = cube.times[k];
        let prev = if k == 0 { 0.0 } else { cube.times[k - 1] };
        let dt = t - prev;
        if dt <= 0.0 {
            continue;
        }
        let df = ois_dfs.get(k).copied().unwrap_or(T::one());
        let sc0 = if k == 0 {
            T::one()
        } else {
            cpty_survival.get(k - 1).copied().unwrap_or(T::one())
        };
        let sc = cpty_survival.get(k).copied().unwrap_or(sc0);
        let sb0 = if k == 0 {
            T::one()
        } else {
            own_survival.get(k - 1).copied().unwrap_or(T::one())
        };
        let sb = own_survival.get(k).copied().unwrap_or(sb0);
        let max_iter = match spec.close_out {
            XvaCloseOutConvention::RiskFree => 1,
            XvaCloseOutConvention::Replacement => 8,
        };
        let mut u_k = match spec.close_out {
            XvaCloseOutConvention::RiskFree => T::zero(),
            XvaCloseOutConvention::Replacement => u_cont / df.max_val(T::from_f64(1e-12)),
        };
        let mut cv = T::zero();
        let mut dv = T::zero();
        let mut fc = T::zero();
        let mut fb = T::zero();
        let mut co = T::zero();
        for _ in 0..max_iter {
            let mut pos = T::zero();
            let mut neg = T::zero();
            let mut xp = T::zero();
            let mut xn = T::zero();
            let u_ns = u_k / ns_count;
            for p in 0..cube.n_paths {
                for ns in 0..cube.n_netting_sets {
                    let v = cube.get(ns, k, p);
                    let x = if let Some(c) = csa {
                        c.collateral + (v - c.collateral - c.threshold).max_val(T::zero())
                            - (-v + c.collateral - c.threshold).max_val(T::zero())
                    } else {
                        T::zero()
                    };
                    let e = v - u_ns - x;
                    pos += e.max_val(T::zero());
                    neg += (-e).max_val(T::zero());
                    xp += x.max_val(T::zero());
                    xn += (-x).max_val(T::zero());
                }
            }
            pos /= den;
            neg /= den;
            xp /= den;
            xn /= den;
            cv = spec.lgd_counterparty * df * sb0 * (sc0 - sc) * pos;
            dv = spec.lgd_own * df * sc0 * (sb0 - sb) * neg;
            fc = spec.funding_borrow_spread * df * sc * sb * pos * T::from_f64(dt);
            fb = spec.funding_lend_spread * df * sc * sb * neg * T::from_f64(dt);
            co = (spec.collateral_pay_spread * xp + spec.collateral_rec_spread * xn)
                * df
                * sc
                * sb
                * T::from_f64(dt);
            if matches!(spec.close_out, XvaCloseOutConvention::Replacement) {
                let step_xva = cv - dv + fc - fb + co;
                let next_u = (u_cont + step_xva) / df.max_val(T::from_f64(1e-12));
                iterations += 1;
                if (next_u - u_k).abs().val() < 1e-12 {
                    u_k = next_u;
                    break;
                }
                u_k = next_u;
            }
        }
        u_cont += cv - dv + fc - fb + co;
        cva += cv;
        dva += dv;
        fca += fc;
        fba += fb;
        colva += co;
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

#[derive(Copy, Clone, Debug)]
pub struct MporCollateralSpec<T: Real = f64> {
    pub csa: CsaSpec,
    pub mpor_years: f64,
    pub unpaid_cashflow_weight: T,
}
#[allow(clippy::too_many_arguments)]
pub fn compute_mpor_xva<T: Real>(
    cube: &NpvCube<T>,
    period_cashflows: Option<&[T]>,
    curve: &impl YieldTermStructure<T>,
    s: &MporCollateralSpec<T>,
    hc: T,
    lc: T,
    hs: T,
    ls: T,
    fs: T,
) -> XvaReport<T> {
    let has_cf = period_cashflows.is_some_and(|c| !c.is_empty());
    if s.mpor_years <= 0.0 && !has_cf {
        return compute_xva(cube, curve, Some(&s.csa), hc, lc, hs, ls, fs);
    }
    let mut adj = NpvCube::new(cube.n_netting_sets, cube.times.clone(), cube.n_paths);
    for ns in 0..cube.n_netting_sets {
        for (k, &t) in cube.times.iter().enumerate() {
            let lag = (t - s.mpor_years.max(0.0)).max(0.0);
            let j = cube
                .times
                .partition_point(|x| *x <= lag)
                .min(cube.times.len().saturating_sub(1));
            let alpha0 = if cube.times[0] > 1e-12 {
                (lag / cube.times[0]).clamp(0.0, 1.0)
            } else {
                1.0
            };
            let alpha = if j == 0 {
                0.0
            } else {
                ((lag - cube.times[j - 1]) / (cube.times[j] - cube.times[j - 1]).max(1e-12))
                    .clamp(0.0, 1.0)
            };
            for p in 0..cube.n_paths {
                let v = cube.get(ns, k, p);
                let vl = if j == 0 {
                    cube.get(ns, 0, p) * alpha0
                } else {
                    cube.get(ns, j - 1, p) * (1.0 - alpha) + cube.get(ns, j, p) * alpha
                };
                let held = vl - s.csa.collateralized_exposure(vl);
                let cf = period_cashflows
                    .and_then(|c| c.get(k))
                    .copied()
                    .unwrap_or_else(T::zero)
                    * s.unpaid_cashflow_weight;
                adj.set(ns, k, p, v + cf - held);
            }
        }
    }
    compute_xva(&adj, curve, None, hc, lc, hs, ls, fs)
}
pub fn price_contingent_cds<T: Real>(
    cube: &NpvCube<T>,
    netting_set_idx: usize,
    csa: Option<&CsaSpec>,
    curve: &impl YieldTermStructure<T>,
    hazard_rate: T,
    loss_given_default: T,
    running_spread: T,
) -> (T, T, T) {
    // Empty cubes/netting sets are valid for report pipelines that build a
    // portfolio incrementally.  Avoid the `n_netting_sets - 1` underflow
    // that previously panicked before returning a well-defined zero quote.
    if cube.n_netting_sets == 0 || cube.n_paths == 0 || cube.times.is_empty() {
        return (T::zero(), T::zero(), T::zero());
    }
    let mut prot = T::zero();
    let mut prem = T::zero();
    let mut prev = T::one();
    for (k, &t) in cube.times.iter().enumerate() {
        let mut e = T::zero();
        for p in 0..cube.n_paths {
            let v = cube.get(netting_set_idx.min(cube.n_netting_sets - 1), k, p);
            let x = if let Some(c) = csa {
                c.collateralized_exposure(v)
            } else {
                v
            };
            e += x.max_val(T::zero());
        }
        e /= T::from_f64(cube.n_paths.max(1) as f64);
        let q = (-hazard_rate * T::from_f64(t)).exp();
        let df = curve.discount_factor(t);
        prot += df * loss_given_default * e * (prev - q);
        let dt = if k == 0 { t } else { t - cube.times[k - 1] };
        prem += df * e * (prev + q) * 0.5 * T::from_f64(dt);
        prev = q;
    }
    let annuity = prem;
    let premium = annuity * running_spread;
    let fair = if annuity.val().abs() > 1e-14 {
        prot / annuity
    } else {
        T::zero()
    };
    (prot, premium, fair)
}

#[derive(Clone, Debug)]
pub struct BrigoVrinsWwmCvaSpec<T: Real = f64> {
    pub times: Vec<f64>,
    pub discount_rate: T,
    pub base_mtm_profile: Vec<T>,
    pub exposure_vol_nu: T,
    pub is_brownian_bridge_swap: bool,
    pub y0: T,
    pub cir_kappa: T,
    pub cir_theta: T,
    pub cir_vol: T,
    pub intensity_shift_psi: T,
    pub wwr_correlation: T,
    pub lgd: T,
}

#[derive(Clone, Debug)]
pub struct BrigoVrinsWwmCvaReport<T: Real = f64> {
    pub wwm_cva: T,
    pub independent_cva: T,
    pub wwr_cva_surcharge: T,
    pub survival_probabilities: Vec<T>,
    pub wwm_drift_shifts: Vec<T>,
    pub wwm_epe_profile: Vec<T>,
}

fn brigo_vrins_cir_survival<T: Real>(y0: T, kappa: T, theta: T, sigma: T, psi: T, t: f64) -> T {
    if t <= 1e-12 {
        return T::one();
    }
    let s2 = (sigma * sigma).max_val(T::from_f64(1e-14));
    let h = (kappa * kappa + s2 * 2.0).sqrt();
    let e_ht = (h * t).exp();
    let denom = h + (kappa + h) * (e_ht - T::one()) * 0.5;
    let b_t = (e_ht - T::one()) / denom;
    let log_a = (kappa * theta * 2.0 / s2)
        * (h * ((kappa + h) * (0.5 * t)).exp() / denom)
            .max_val(T::from_f64(1e-14))
            .ln();
    (log_a - b_t * y0 - psi * t).exp()
}

fn brigo_vrins_theta_drift<T: Real>(
    y0: T,
    kappa: T,
    theta: T,
    sigma_lam: T,
    psi: T,
    rho: T,
    nu_exp: T,
    s: f64,
    t: f64,
) -> T {
    let tau = (t - s).max(0.0);
    let s2 = (sigma_lam * sigma_lam).max_val(T::from_f64(1e-14));
    let h = (kappa * kappa + s2 * 2.0).sqrt();
    let e_ht = (h * tau).exp();
    let denom = h + (kappa + h) * (e_ht - T::one()) * 0.5;
    let b_cir = (e_ht - T::one()) / denom;
    let h_over_d = h / denom;
    let bt_cir = e_ht * h_over_d * h_over_d;
    let e_ks = (-kappa * s).exp();
    let lam_bar = (psi + y0 * e_ks + theta * (T::one() - e_ks)).max_val(T::from_f64(1e-8));
    let denom_wwm =
        (bt_cir * lam_bar - psi + bt_cir * psi + kappa * theta * b_cir).max_val(T::from_f64(1e-12));
    rho * nu_exp * sigma_lam * lam_bar.sqrt() * (bt_cir / denom_wwm - b_cir)
}

/// Brigo & Vrins (2018) Wrong-Way Measure (`WWM`) CIR++ drift-adjusted analytical CVA
/// without Monte Carlo (`papers/48_Brigo_Vrins_2018_Disentangling_Wrong_Way_Risk_CVA.pdf`, Eqs. 16-21).
pub fn compute_brigo_vrins_wwm_cva<T: Real>(
    spec: &BrigoVrinsWwmCvaSpec<T>,
) -> BrigoVrinsWwmCvaReport<T> {
    const GL4_X: [f64; 4] = [
        -0.861_136_311_594_052_6,
        -0.339_981_043_584_856_3,
        0.339_981_043_584_856_3,
        0.861_136_311_594_052_6,
    ];
    const GL4_W: [f64; 4] = [
        0.347_854_845_137_453_8,
        0.652_145_154_862_546_1,
        0.652_145_154_862_546_1,
        0.347_854_845_137_453_8,
    ];
    let n = spec.times.len().min(spec.base_mtm_profile.len());
    if n == 0 {
        return BrigoVrinsWwmCvaReport {
            wwm_cva: T::zero(),
            independent_cva: T::zero(),
            wwr_cva_surcharge: T::zero(),
            survival_probabilities: Vec::new(),
            wwm_drift_shifts: Vec::new(),
            wwm_epe_profile: Vec::new(),
        };
    }
    let t_max = (spec.times[n - 1] + 0.5).max(1.0);
    let mut wwm_cva = T::zero();
    let mut ind_cva = T::zero();
    let mut surv_probs = Vec::with_capacity(n);
    let mut drift_shifts = Vec::with_capacity(n);
    let mut wwm_epe = Vec::with_capacity(n);
    let mut prev_g = T::one();

    for k in 0..n {
        let tk = spec.times[k].max(1e-4);
        let g_k = brigo_vrins_cir_survival(
            spec.y0,
            spec.cir_kappa,
            spec.cir_theta,
            spec.cir_vol,
            spec.intensity_shift_psi,
            tk,
        );
        let mid = 0.5 * tk;
        let half = 0.5 * tk;
        let mut shift = T::zero();
        for m in 0..4 {
            let u = mid + half * GL4_X[m];
            let th_u = brigo_vrins_theta_drift(
                spec.y0,
                spec.cir_kappa,
                spec.cir_theta,
                spec.cir_vol,
                spec.intensity_shift_psi,
                spec.wwr_correlation,
                spec.exposure_vol_nu,
                u,
                tk,
            );
            if spec.is_brownian_bridge_swap {
                let bridge_factor = (t_max - tk) / (t_max - u).max(1e-6);
                shift += th_u * (bridge_factor * half * GL4_W[m]);
            } else {
                shift += th_u * (half * GL4_W[m]);
            }
        }
        let sig_v = if spec.is_brownian_bridge_swap {
            (spec.exposure_vol_nu * (tk * (1.0 - tk / t_max).max(0.05)).sqrt())
                .max_val(T::from_f64(1e-8))
        } else {
            (spec.exposure_vol_nu * tk.sqrt()).max_val(T::from_f64(1e-8))
        };
        let mu_ind = spec.base_mtm_profile[k];
        let mu_wwm = mu_ind + shift;
        let d_ind = mu_ind / sig_v;
        let d_wwm = mu_wwm / sig_v;
        let epe_ind = mu_ind * d_ind.norm_cdf() + sig_v * d_ind.norm_pdf();
        let epe_w = mu_wwm * d_wwm.norm_cdf() + sig_v * d_wwm.norm_pdf();
        let df = (-spec.discount_rate * tk).exp();
        let dg = (prev_g - g_k).max_val(T::zero());
        ind_cva += spec.lgd * df * epe_ind * dg;
        wwm_cva += spec.lgd * df * epe_w * dg;

        surv_probs.push(g_k);
        drift_shifts.push(shift);
        wwm_epe.push(epe_w);
        prev_g = g_k;
    }

    BrigoVrinsWwmCvaReport {
        wwm_cva,
        independent_cva: ind_cva,
        wwr_cva_surcharge: wwm_cva - ind_cva,
        survival_probabilities: surv_probs,
        wwm_drift_shifts: drift_shifts,
        wwm_epe_profile: wwm_epe,
    }
}

#[derive(Clone, Debug)]
pub struct AsymmetricRateXvaBandSpec<T: Real = f64> {
    pub maturity: f64,
    pub num_steps: usize,
    pub clean_claim_value: T,
    pub clean_delta_notional: T,
    pub asset_vol: T,
    pub collateral_alpha: T,
    pub dividend_discount_rate: T,
    pub repo_lend_rate: T,
    pub repo_borrow_rate: T,
    pub funding_lend_rate: T,
    pub funding_borrow_rate: T,
    pub collateral_lend_rate: T,
    pub collateral_borrow_rate: T,
    pub own_bond_yield_mu_i: T,
    pub cpty_bond_yield_mu_c: T,
    pub own_lgd: T,
    pub cpty_lgd: T,
}

#[derive(Copy, Clone, Debug)]
pub struct AsymmetricRateXvaBandReport<T: Real = f64> {
    pub buyer_sub_hedging_xva: T,
    pub seller_super_hedging_xva: T,
    pub symmetric_benchmark_xva: T,
    pub buyer_reservation_price: T,
    pub seller_reservation_price: T,
    pub no_arbitrage_band_width: T,
    pub replicating_stock_delta_adj: T,
    pub replicating_own_bond_units: T,
    pub replicating_cpty_bond_units: T,
    pub no_arbitrage_satisfied: bool,
}

fn smooth_pos<T: Real>(x: T) -> T {
    let eps = T::from_f64(1e-6);
    (x + (x * x + eps * eps).sqrt()) * 0.5
}

/// Bichuch, Capponi & Sturm (2018) Arbitrage-Free Bilateral XVA Band $[V^-, V^+]$ & Risky-Bond
/// Replication (`papers/56_Bichuch_Capponi_Sturm_2018_Arbitrage_Free_XVA_Funding.pdf`, Theorems 4.7, 4.10 & Prop. 5.3).
pub fn compute_bichuch_capponi_sturm_xva_band<T: Real>(
    s: &AsymmetricRateXvaBandSpec<T>,
) -> AsymmetricRateXvaBandReport<T> {
    let t_mat = s.maturity.max(1e-4);
    let n_steps = s.num_steps.max(4);
    let dt = t_mat / n_steps as f64;
    let v_hat = s.clean_claim_value;
    let z_hat = s.clean_delta_notional * s.asset_vol;
    let alpha = s.collateral_alpha.max_val(T::zero()).min_val(T::one());
    let one_m_alpha = T::one() - alpha;

    let rf_plus = s.funding_lend_rate;
    let rf_minus = s.funding_borrow_rate;
    let rr_plus = s.repo_lend_rate;
    let rr_minus = s.repo_borrow_rate;
    let rc_plus = s.collateral_lend_rate;
    let rc_minus = s.collateral_borrow_rate;
    let mu_i = s.own_bond_yield_mu_i;
    let mu_c = s.cpty_bond_yield_mu_c;

    let no_arbitrage_satisfied = rr_plus.val() <= rr_minus.val() + 1e-12
        && rf_plus.val() <= rf_minus.val() + 1e-12
        && rr_plus.val() <= rf_plus.val() + 1e-12
        && rf_plus.val() <= rr_minus.val() + 1e-12
        && rf_minus.val() <= mu_i.val().min(mu_c.val()) + 1e-12;

    // Proposition 5.3 Symmetric-Rate Closed-Form Benchmark & Replicating Portfolio
    let rf_mid = (rf_plus + rf_minus) * 0.5;
    let rr_mid = (rr_plus + rr_minus) * 0.5;
    let rc_mid = (rc_plus + rc_minus) * 0.5;
    let eta = mu_i + mu_c - rf_mid;
    let rate_diff = (eta - rr_mid).max_val(T::from_f64(1e-8));
    let a_factor = (T::one() - (-rate_diff * t_mat).exp()) / rate_diff;
    let uncoll_pos = smooth_pos(one_m_alpha * v_hat);
    let uncoll_neg = smooth_pos(-one_m_alpha * v_hat);
    let sym_flow = ((rr_mid - rf_mid) + alpha * (rf_mid - rc_mid)) * v_hat
        + (mu_c - rf_mid) * s.cpty_lgd * uncoll_neg
        - (mu_i - rf_mid) * s.own_lgd * uncoll_pos;
    let sym_xva = a_factor * sym_flow;

    let ind_pos = smooth_pos(v_hat) / (v_hat.abs() + T::from_f64(1e-6));
    let ind_neg = smooth_pos(-v_hat) / (v_hat.abs() + T::from_f64(1e-6));
    let xi_stock = a_factor
        * ((rr_mid - rf_mid)
            + alpha * (rf_mid - rc_mid)
            + (mu_c - rf_mid) * s.cpty_lgd * one_m_alpha * ind_neg
            - (mu_i - rf_mid) * s.own_lgd * one_m_alpha * ind_pos)
        * s.clean_delta_notional;
    let xi_own_bond = (sym_xva + s.own_lgd * uncoll_pos) / (-mu_i * t_mat).exp();
    let xi_cpty_bond = (sym_xva - s.cpty_lgd * uncoll_neg) / (-mu_c * t_mat).exp();

    // Theorem 4.7 Reduced BSDE Drivers f^+(u) & f^-(u) integrated via exact exponential-Euler Picard recursion
    let rf_spread_half = (rf_minus - rf_plus).max_val(T::zero()) * 0.5;
    let rr_spread_half = (rr_minus - rr_plus).max_val(T::zero()) * 0.5;
    let rc_spread_half = (rc_minus - rc_plus).max_val(T::zero()) * 0.5;
    let theta_net = s.cpty_lgd * uncoll_neg - s.own_lgd * uncoll_pos;
    let abs_delta = (z_hat / s.asset_vol.max_val(T::from_f64(1e-6))).abs();
    let abs_coll = (alpha * v_hat).abs();
    let step_decay = (-rate_diff * dt).exp();
    let step_int = (T::one() - step_decay) / rate_diff;

    let mut u_plus = T::zero();
    let mut u_minus = T::zero();
    for _ in 0..n_steps {
        let mut next_p = u_plus;
        let mut next_m = u_minus;
        for _ in 0..6 {
            let gamma_p = theta_net - next_p + one_m_alpha * v_hat;
            let abs_gamma_p = smooth_pos(gamma_p) + smooth_pos(-gamma_p);
            let asym_p = rf_spread_half * abs_gamma_p
                + rr_spread_half * abs_delta
                + rc_spread_half * abs_coll;
            next_p = u_plus * step_decay + (sym_flow + asym_p) * step_int;

            let gamma_m = theta_net - next_m + one_m_alpha * v_hat;
            let abs_gamma_m = smooth_pos(gamma_m) + smooth_pos(-gamma_m);
            let asym_m = rf_spread_half * abs_gamma_m
                + rr_spread_half * abs_delta
                + rc_spread_half * abs_coll;
            next_m = u_minus * step_decay + (sym_flow - asym_m) * step_int;
        }
        u_plus = next_p;
        u_minus = next_m;
    }

    let seller_super = u_plus.max_val(u_minus);
    let buyer_sub = u_minus.min_val(u_plus);
    let width = (seller_super - buyer_sub).max_val(T::zero());

    AsymmetricRateXvaBandReport {
        buyer_sub_hedging_xva: buyer_sub,
        seller_super_hedging_xva: seller_super,
        symmetric_benchmark_xva: sym_xva,
        buyer_reservation_price: v_hat + buyer_sub,
        seller_reservation_price: v_hat + seller_super,
        no_arbitrage_band_width: width,
        replicating_stock_delta_adj: xi_stock,
        replicating_own_bond_units: xi_own_bond,
        replicating_cpty_bond_units: xi_cpty_bond,
        no_arbitrage_satisfied,
    }
}

#[derive(Clone, Debug)]
pub struct CrepeyBieleckiMarshallOlkinBilateralXvaSpec<T: Real = f64> {
    pub times: Vec<T>,
    pub expected_positive_exposure: Vec<T>,
    pub expected_negative_exposure: Vec<T>,
    pub risk_free_rate: T,
    pub lambda_institution_idio: T,
    pub lambda_counterparty_idio: T,
    pub lambda_common_shock_12: T,
    pub recovery_institution: T,
    pub recovery_counterparty: T,
    pub recovery_counterparty_simultaneous: T,
    pub gap_exposure_multiplier: T,
    pub funding_borrowing_spread: T,
    pub funding_lending_spread: T,
}
#[derive(Clone, Debug)]
pub struct CrepeyBieleckiMarshallOlkinBilateralXvaReport<T: Real = f64> {
    pub total_first_to_default_intensity: T,
    pub kendall_tau_default_correlation: T,
    pub simultaneous_default_prob_to_maturity: T,
    pub cva_idiosyncratic: T,
    pub cva_common_shock_contagion: T,
    pub bilateral_cva_total: T,
    pub bilateral_dva: T,
    pub bilateral_fva: T,
    pub total_bilateral_xva: T,
    pub independent_copula_cva_bias: T,
}
pub fn compute_crepey_bielecki_marshall_olkin_bilateral_xva<T: Real>(
    s: &CrepeyBieleckiMarshallOlkinBilateralXvaSpec<T>,
) -> CrepeyBieleckiMarshallOlkinBilateralXvaReport<T> {
    let lb = s.lambda_institution_idio.max_val(T::zero());
    let lc = s.lambda_counterparty_idio.max_val(T::zero());
    let l12 = s.lambda_common_shock_12.max_val(T::zero());
    let lt = (lb + lc + l12).max_val(T::from_f64(1e-12));
    let tau = l12 / lt;
    let mut ci = T::zero();
    let mut cc = T::zero();
    let mut dva = T::zero();
    let mut fva = T::zero();
    let mut indep = T::zero();
    let mut prev = T::zero();
    let n = s
        .times
        .len()
        .min(s.expected_positive_exposure.len())
        .min(s.expected_negative_exposure.len());
    for k in 0..n {
        let t = s.times[k];
        let dt = (t - prev).max_val(T::zero());
        let sk = (-lt * t).exp();
        let sp = (-lt * prev).exp();
        let dp = sp - sk;
        let d = (-s.risk_free_rate * t).exp();
        let e = s.expected_positive_exposure[k];
        let en = s.expected_negative_exposure[k];
        ci += (T::one() - s.recovery_counterparty) * d * dp * (lc / lt) * e;
        cc += (T::one() - s.recovery_counterparty_simultaneous)
            * d
            * dp
            * (l12 / lt)
            * (T::one() + s.gap_exposure_multiplier)
            * e;
        dva += (T::one() - s.recovery_institution) * d * dp * ((lb + l12) / lt) * en;
        fva += d * sk * dt * (s.funding_borrowing_spread * e - s.funding_lending_spread * en);
        let li = (lb + l12 + lc + l12).max_val(T::from_f64(1e-12));
        let si = (-li * t).exp();
        let spi = (-li * prev).exp();
        indep += (T::one() - s.recovery_counterparty) * d * (spi - si) * ((lc + l12) / li) * e;
        prev = t;
    }
    let ctot = ci + cc;
    let indep_bias = ctot - indep;
    let maturity = s.times.last().copied().unwrap_or(T::zero());
    CrepeyBieleckiMarshallOlkinBilateralXvaReport {
        total_first_to_default_intensity: lt,
        kendall_tau_default_correlation: tau,
        simultaneous_default_prob_to_maturity: (l12 / lt) * (T::one() - (-lt * maturity).exp()),
        cva_idiosyncratic: ci,
        cva_common_shock_contagion: cc,
        bilateral_cva_total: ctot,
        bilateral_dva: dva,
        bilateral_fva: fva,
        total_bilateral_xva: ctot - dva + fva,
        independent_copula_cva_bias: indep_bias,
    }
}
pub fn evaluate_crepey_bielecki_brigo_marshall_olkin_xva<T: Real>(
    s: &CrepeyBieleckiMarshallOlkinBilateralXvaSpec<T>,
) -> CrepeyBieleckiMarshallOlkinBilateralXvaReport<T> {
    compute_crepey_bielecki_marshall_olkin_bilateral_xva(s)
}

#[cfg(test)]
mod phase23_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn ctd_bounds() {
        let s = BrigoPallaviciniCtdCollateralXvaSpec {
            notional_mtm: 1.0,
            maturity: 1.0,
            num_steps: 8,
            base_ois_rate: 0.02,
            csa_spread_ccy1: 0.001,
            csa_spread_ccy2: 0.002,
            csa_spread_ccy3: 0.003,
            csa_mean_reversion: 1.0,
            csa_vol_ccy1: 0.01,
            csa_vol_ccy2: 0.01,
            csa_vol_ccy3: 0.01,
            csa_correlation: 0.2,
            cpty_hazard_rate: 0.02,
            cpty_lgd: 0.6,
            own_hazard_rate: 0.01,
            own_lgd: 0.6,
            exposure_volatility: 0.1,
        };
        let r = compute_brigo_pallavicini_ctd_collateral_and_closeout_xva(&s);
        assert!(r.ctd_stochastic_csa_discount_factor <= r.intrinsic_max_csa_discount_factor);
        let same = BrigoPallaviciniCtdCollateralXvaSpec {
            csa_spread_ccy1: 0.002,
            csa_spread_ccy2: 0.002,
            csa_spread_ccy3: 0.002,
            csa_correlation: 0.999999,
            ..s
        };
        let same_r = compute_brigo_pallavicini_ctd_collateral_and_closeout_xva(&same);
        assert!(same_r.ctd_collateral_switch_option_value < 1e-5);
    }

    #[test]
    fn dual_fd_parity() {
        fn mk<T: Real>(x: T) -> BrigoPallaviciniCtdCollateralXvaSpec<T> {
            BrigoPallaviciniCtdCollateralXvaSpec {
                notional_mtm: T::one(),
                maturity: T::one(),
                num_steps: 8,
                base_ois_rate: T::from_f64(0.02),
                csa_spread_ccy1: T::from_f64(0.001),
                csa_spread_ccy2: T::from_f64(0.002),
                csa_spread_ccy3: T::from_f64(0.003),
                csa_mean_reversion: T::one(),
                csa_vol_ccy1: T::from_f64(0.01),
                csa_vol_ccy2: T::from_f64(0.01),
                csa_vol_ccy3: T::from_f64(0.01),
                csa_correlation: T::from_f64(0.2),
                cpty_hazard_rate: T::from_f64(0.02),
                cpty_lgd: T::from_f64(0.6),
                own_hazard_rate: T::from_f64(0.01),
                own_lgd: T::from_f64(0.6),
                exposure_volatility: x,
            }
        }
        let h = 1e-5;
        let fd = (compute_brigo_pallavicini_ctd_collateral_and_closeout_xva(&mk(0.1 + h))
            .replacement_cost_closeout_bcva
            .val()
            - compute_brigo_pallavicini_ctd_collateral_and_closeout_xva(&mk(0.1 - h))
                .replacement_cost_closeout_bcva
                .val())
            / (2.0 * h);
        let ad = compute_brigo_pallavicini_ctd_collateral_and_closeout_xva(&mk(
            Dual::<1>::variable(0.1, 0),
        ))
        .replacement_cost_closeout_bcva;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase14_tests {
    use super::*;
    use crate::NpvCube;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;
    #[test]
    fn phase14_correlated_ftd_is_finite() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let curve = FlatForwardCurve::new(d, 0.02);
        let mut cube = NpvCube::new(1, vec![0.5, 1.0], 2);
        cube.set(0, 0, 0, 10.0);
        cube.set(0, 0, 1, 12.0);
        cube.set(0, 1, 0, 11.0);
        cube.set(0, 1, 1, 13.0);
        let r = compute_correlated_first_to_default_bilateral_xva(
            &cube, &curve, None, 0.02, 0.6, 0.03, 0.6, 0.01, 0.0,
        );
        assert!(r.cva.is_finite() && r.dva.is_finite() && r.fva.is_finite());
    }
}

#[cfg(test)]
mod phase13_tests {
    use super::*;
    #[test]
    fn contingent_empty_safe() {
        let c = NpvCube::<f64>::new(0, vec![], 0);
        let d = crate::price_contingent_cds(
            &c,
            0,
            None,
            &ql_curves::FlatForwardCurve::new(ql_time::Date(0), 0.02),
            0.01,
            1.,
            0.01,
        );
        assert_eq!(d.0, 0.);
    }
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn ctd_colva_clark_recursion_and_turlakov_wwr_ad() {
        let c1 = compute_multicurrency_ctd_colva(&MultiCurrencyCsaSpec {
            eligible_currency_spreads: vec![0.01, 0.01, 0.01],
            spread_vols: vec![0.008, 0.008, 0.008],
            cross_currency_corr: 0.2,
            collateral_balance_profile: vec![100.0, 100.0],
            payment_times: vec![0.5, 1.0],
        });
        let c2 = compute_multicurrency_ctd_colva(&MultiCurrencyCsaSpec {
            eligible_currency_spreads: vec![0.01, 0.01],
            spread_vols: vec![0.008, 0.008],
            cross_currency_corr: 0.2,
            collateral_balance_profile: vec![100.0, 100.0],
            payment_times: vec![0.5, 1.0],
        });
        assert!(c2.switch_option_colva > 0.0);
        assert!(c1.switch_option_colva > c2.switch_option_colva);

        let wwr_spec = TurlakovWwrSpec {
            times: vec![0.5, 1.0, 2.0],
            base_epe: vec![10.0, 15.0, 12.0],
            stressed_epe: vec![25.0, 35.0, 30.0],
            base_mtm: vec![8.0, 10.0, 9.0],
            stressed_mtm: vec![20.0, 25.0, 22.0],
            sovereign_hazard_rate: 0.03,
            counterparty_hazard_rate: 0.02,
            crisis_prob_threshold: 0.08,
            coupling_lambda: 0.5,
            lgd: 0.6,
            base_funding_spread: 0.01,
            stressed_funding_spread: 0.025,
            discount_rate: 0.03,
        };
        let rep = compute_turlakov_wwr_cva_fva(&wwr_spec);
        assert!(rep.wwr_cva > rep.base_cva);
        assert!(rep.wwr_fva > rep.base_fva);

        let h = 1e-5;
        let eval_lam = |lam: f64| {
            let mut sp = wwr_spec.clone();
            sp.coupling_lambda = lam;
            compute_turlakov_wwr_cva_fva(&sp).wwr_cva
        };
        let fd = (eval_lam(0.5 + h) - eval_lam(0.5 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sd = TurlakovWwrSpec {
            times: vec![0.5, 1.0, 2.0],
            base_epe: vec![
                D::from_value(10.0),
                D::from_value(15.0),
                D::from_value(12.0),
            ],
            stressed_epe: vec![
                D::from_value(25.0),
                D::from_value(35.0),
                D::from_value(30.0),
            ],
            base_mtm: vec![D::from_value(8.0), D::from_value(10.0), D::from_value(9.0)],
            stressed_mtm: vec![
                D::from_value(20.0),
                D::from_value(25.0),
                D::from_value(22.0),
            ],
            sovereign_hazard_rate: D::from_value(0.03),
            counterparty_hazard_rate: D::from_value(0.02),
            crisis_prob_threshold: D::from_value(0.08),
            coupling_lambda: D::variable(0.5, 0),
            lgd: D::from_value(0.6),
            base_funding_spread: D::from_value(0.01),
            stressed_funding_spread: D::from_value(0.025),
            discount_rate: D::from_value(0.03),
        };
        let ad = compute_turlakov_wwr_cva_fva(&sd).wwr_cva;
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase20_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn brigo_vrins_wwm_zero_corr_collapse_and_bichuch_capponi_sturm_band_ad() {
        let base_wwm = BrigoVrinsWwmCvaSpec {
            times: vec![0.5, 1.0, 2.0, 3.0],
            discount_rate: 0.03,
            base_mtm_profile: vec![10.0, 12.0, 14.0, 11.0],
            exposure_vol_nu: 8.0,
            is_brownian_bridge_swap: false,
            y0: 0.02,
            cir_kappa: 0.35,
            cir_theta: 0.025,
            cir_vol: 0.08,
            intensity_shift_psi: 0.005,
            wwr_correlation: 0.0,
            lgd: 0.60,
        };
        let r0 = compute_brigo_vrins_wwm_cva(&base_wwm);
        assert!((r0.wwm_cva - r0.independent_cva).abs() < 1e-12);
        assert!(r0.wwr_cva_surcharge.abs() < 1e-12);

        let pos_wwm = BrigoVrinsWwmCvaSpec {
            wwr_correlation: 0.65,
            ..base_wwm.clone()
        };
        let r_pos = compute_brigo_vrins_wwm_cva(&pos_wwm);
        assert!(r_pos.wwm_cva > r_pos.independent_cva);
        assert!(r_pos.wwr_cva_surcharge > 0.0);

        // Bichuch-Capponi-Sturm (2018): symmetric rates => zero band width; asymmetric borrowing > lending => positive band width
        let sym_bcs = AsymmetricRateXvaBandSpec {
            maturity: 1.0,
            num_steps: 10,
            clean_claim_value: 25.0,
            clean_delta_notional: 40.0,
            asset_vol: 0.20,
            collateral_alpha: 0.50,
            dividend_discount_rate: 0.02,
            repo_lend_rate: 0.02,
            repo_borrow_rate: 0.02,
            funding_lend_rate: 0.020,
            funding_borrow_rate: 0.020,
            collateral_lend_rate: 0.02,
            collateral_borrow_rate: 0.02,
            own_bond_yield_mu_i: 0.04,
            cpty_bond_yield_mu_c: 0.045,
            own_lgd: 0.60,
            cpty_lgd: 0.60,
        };
        let rep_sym = compute_bichuch_capponi_sturm_xva_band(&sym_bcs);
        assert!(rep_sym.no_arbitrage_satisfied);
        assert!(rep_sym.no_arbitrage_band_width.abs() < 1e-10);

        let asym_bcs = AsymmetricRateXvaBandSpec {
            funding_lend_rate: 0.022,
            funding_borrow_rate: 0.038,
            repo_lend_rate: 0.020,
            repo_borrow_rate: 0.025,
            ..sym_bcs
        };
        let rep_asym = compute_bichuch_capponi_sturm_xva_band(&asym_bcs);
        assert!(rep_asym.no_arbitrage_satisfied);
        assert!(rep_asym.no_arbitrage_band_width > 0.0);
        assert!(rep_asym.seller_reservation_price > rep_asym.buyer_reservation_price);

        let h = 1e-5;
        let eval_rho = |rho: f64| {
            let mut s = pos_wwm.clone();
            s.wwr_correlation = rho;
            compute_brigo_vrins_wwm_cva(&s).wwm_cva
        };
        let fd = (eval_rho(0.65 + h) - eval_rho(0.65 - h)) / (2.0 * h);
        type D = Dual<1>;
        let spec_d = BrigoVrinsWwmCvaSpec {
            times: vec![0.5, 1.0, 2.0, 3.0],
            discount_rate: D::from_value(0.03),
            base_mtm_profile: vec![
                D::from_value(10.0),
                D::from_value(12.0),
                D::from_value(14.0),
                D::from_value(11.0),
            ],
            exposure_vol_nu: D::from_value(8.0),
            is_brownian_bridge_swap: false,
            y0: D::from_value(0.02),
            cir_kappa: D::from_value(0.35),
            cir_theta: D::from_value(0.025),
            cir_vol: D::from_value(0.08),
            intensity_shift_psi: D::from_value(0.005),
            wwr_correlation: D::variable(0.65, 0),
            lgd: D::from_value(0.60),
        };
        let ad = compute_brigo_vrins_wwm_cva(&spec_d).wwm_cva;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase21_tests {
    use super::*;

    #[test]
    fn gregory_euler_allocations_add_up() {
        let s = GregoryNettingSetAllocationSpec {
            trade_weights: vec![1.0, 0.8],
            trade_expected_mtms: vec![1.0, -0.3],
            trade_exposure_vols: vec![0.2, 0.25],
            correlation_matrix: vec![vec![1.0, 0.2], vec![0.2, 1.0]],
            collateral_fraction: 0.1,
            horizon: 1.0,
            num_steps: 8,
            risk_free_rate: 0.02,
            cpty_hazard_rate: 0.02,
            cpty_lgd: 0.6,
            own_hazard_rate: 0.01,
            own_lgd: 0.6,
            funding_spread: 0.01,
            colva_spread: 0.005,
        };
        let r = allocate_gregory_netting_set_xva(&s);
        let sum = r.euler_allocated_total_xva.iter().copied().sum::<f64>();
        assert!((sum - r.netted_total_xva).abs() < 1e-13);
        assert!(r.netted_total_xva <= r.standalone_total_xva + 1e-12);
    }
}

#[cfg(test)]
mod phase24_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(x: T) -> CrepeyBieleckiMarshallOlkinBilateralXvaSpec<T> {
        CrepeyBieleckiMarshallOlkinBilateralXvaSpec {
            times: vec![T::from_f64(0.5), T::from_f64(1.)],
            expected_positive_exposure: vec![T::from_f64(1.), T::from_f64(1.)],
            expected_negative_exposure: vec![T::from_f64(0.5), T::from_f64(0.5)],
            risk_free_rate: T::from_f64(0.02),
            lambda_institution_idio: T::from_f64(0.01),
            lambda_counterparty_idio: T::from_f64(0.02),
            lambda_common_shock_12: x,
            recovery_institution: T::from_f64(0.4),
            recovery_counterparty: T::from_f64(0.4),
            recovery_counterparty_simultaneous: T::from_f64(0.3),
            gap_exposure_multiplier: T::from_f64(0.2),
            funding_borrowing_spread: T::from_f64(0.01),
            funding_lending_spread: T::from_f64(0.005),
        }
    }
    #[test]
    fn phase24() {
        let r = compute_crepey_bielecki_marshall_olkin_bilateral_xva(&mk(0.0));
        assert_eq!(r.kendall_tau_default_correlation, 0.0);
        assert_eq!(r.cva_common_shock_contagion, 0.0);
        assert!(
            compute_crepey_bielecki_marshall_olkin_bilateral_xva(&mk(0.01))
                .cva_common_shock_contagion
                > 0.0
        );
        let h = 1e-5;
        let fd = (compute_crepey_bielecki_marshall_olkin_bilateral_xva(&mk(0.01 + h))
            .total_bilateral_xva
            - compute_crepey_bielecki_marshall_olkin_bilateral_xva(&mk(0.01 - h))
                .total_bilateral_xva)
            / (2.0 * h);
        let ad =
            compute_crepey_bielecki_marshall_olkin_bilateral_xva(&mk(Dual::<1>::variable(0.01, 0)))
                .total_bilateral_xva;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase27_tests {
    use super::*;
    use ql_ad::Dual;
    fn make<T: Real>(x: T) -> AlbaneseCrepeyMporBalanceSheetSpec<T> {
        AlbaneseCrepeyMporBalanceSheetSpec {
            exposure_times: vec![T::one(), T::from_f64(2.0)],
            clean_mtm_profile: vec![x, x],
            portfolio_vol_profile: vec![T::from_f64(0.2), T::from_f64(0.2)],
            scheduled_cashflows_in_mpor: vec![T::from_f64(0.01), T::from_f64(0.02)],
            mpor_years: T::from_f64(0.05),
            csa_threshold: T::from_f64(0.1),
            risk_free_rate: T::from_f64(0.03),
            cpty_hazard_rate: T::from_f64(0.02),
            cpty_lgd: T::from_f64(0.6),
            bank_hazard_rate: T::from_f64(0.01),
            bank_lgd: T::from_f64(0.6),
            economic_capital_profile: vec![T::from_f64(0.05), T::from_f64(0.05)],
            hurdle_rate: T::from_f64(0.1),
            initial_margin_profile: vec![T::from_f64(0.1), T::from_f64(0.1)],
            im_funding_spread: T::from_f64(0.01),
        }
    }
    #[test]
    fn dual_fd_parity() {
        let h = 1e-5;
        let fd = (compute_albanese_crepey_mpor_and_balance_sheet_xva(&make(1.0 + h))
            .sticky_mpor_cva
            - compute_albanese_crepey_mpor_and_balance_sheet_xva(&make(1.0 - h)).sticky_mpor_cva)
            / (2.0 * h);
        let d =
            compute_albanese_crepey_mpor_and_balance_sheet_xva(&make(Dual::<1>::variable(1.0, 0)))
                .sticky_mpor_cva;
        assert!(d.d[0].abs() > 1e-6 && (d.d[0] - fd).abs() < 1e-5);
    }

    #[test]
    fn phase27_invariants() {
        let r = compute_albanese_crepey_mpor_and_balance_sheet_xva(&make(1.0));
        assert!(r.sticky_mpor_cva > r.classical_mpor_cva);
        assert!(
            (r.albanese_crepey_shareholder_fca + r.creditor_fda_wealth_transfer
                - r.unadjusted_fca_without_capital_offset)
                .abs()
                < 1e-12
        );
    }
}
