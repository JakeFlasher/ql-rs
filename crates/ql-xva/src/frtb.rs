use ql_ad::Real;
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum FrtbRiskClass {
    Girr,
    CsrNonSec,
    Equity,
    Commodity,
    Fx,
}
#[derive(Clone, Debug)]
pub struct FrtbSensitivityBucket<T: Real> {
    pub risk_class: FrtbRiskClass,
    pub bucket_id: usize,
    pub delta_sensitivities: Vec<T>,
    pub delta_risk_weights: Vec<f64>,
    pub vega_sensitivities: Vec<T>,
    pub vega_risk_weights: Vec<f64>,
    pub curvature_up_pnl: Vec<T>,
    pub curvature_down_pnl: Vec<T>,
    pub curvature_shock_weights: Vec<f64>,
    pub intra_bucket_corr: Vec<Vec<f64>>,
}
#[derive(Clone, Debug)]
pub struct FrtbDrcObligor<T: Real> {
    pub bucket_id: usize,
    pub notional: T,
    pub mark_to_market_pnl: T,
    pub lgd: f64,
    pub maturity: f64,
    pub default_risk_weight: f64,
}
#[derive(Clone, Debug)]
pub struct FrtbSaSpec<T: Real> {
    pub buckets: Vec<FrtbSensitivityBucket<T>>,
    pub inter_bucket_corr: Vec<(FrtbRiskClass, usize, usize, f64)>,
    pub drc_obligors: Vec<FrtbDrcObligor<T>>,
    pub exotic_residual_notionals: Vec<T>,
    pub other_residual_notionals: Vec<T>,
}
#[derive(Clone, Debug)]
pub struct FrtbSaReport<T: Real> {
    pub sbm_low: T,
    pub sbm_med: T,
    pub sbm_high: T,
    pub sbm_charge: T,
    pub drc_charge: T,
    pub rrao_charge: T,
    pub total_capital: T,
}
fn corr_scen(rho: f64, scale: f64) -> f64 {
    if scale < 1.0 {
        (2.0 * rho - 1.0).max(0.75 * rho).clamp(-1.0, 1.0)
    } else if scale > 1.0 {
        (1.25 * rho).clamp(-1.0, 1.0)
    } else {
        rho.clamp(-1.0, 1.0)
    }
}

fn bucket_quad<T: Real>(sens: &[T], weights: &[f64], corr: &[Vec<f64>], scale: f64) -> (T, T) {
    let n = sens.len().min(weights.len());
    let mut ws = Vec::with_capacity(n);
    let mut s_sum = T::zero();
    for i in 0..n {
        let w = sens[i] * T::from_f64(weights[i]);
        s_sum += w;
        ws.push(w);
    }
    let mut q = T::zero();
    for i in 0..n {
        for j in 0..n {
            let c = if i == j {
                1.0
            } else {
                corr_scen(
                    corr.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0),
                    scale,
                )
            };
            q += T::from_f64(c) * ws[i] * ws[j];
        }
    }
    (q.max_val(T::zero()).sqrt(), s_sum)
}

fn bucket_curv<T: Real>(b: &FrtbSensitivityBucket<T>, scale: f64) -> (T, T) {
    let n = b
        .curvature_up_pnl
        .len()
        .min(b.curvature_down_pnl.len())
        .min(b.curvature_shock_weights.len());
    if n == 0 {
        return (T::zero(), T::zero());
    }
    let mut c_up = Vec::with_capacity(n);
    let mut c_dn = Vec::with_capacity(n);
    for k in 0..n {
        let rw = T::from_f64(b.curvature_shock_weights[k]);
        let d = b.delta_sensitivities.get(k).copied().unwrap_or(T::zero());
        c_up.push(-b.curvature_up_pnl[k] + rw * d);
        c_dn.push(-b.curvature_down_pnl[k] - rw * d);
    }
    let eval_dir = |cvr: &[T]| -> (T, T) {
        let s_sum = cvr.iter().copied().fold(T::zero(), |a, z| a + z);
        let mut q = T::zero();
        for i in 0..n {
            let pos = cvr[i].max_val(T::zero());
            q += pos * pos;
            for j in 0..n {
                if i != j && !(cvr[i].val() < 0.0 && cvr[j].val() < 0.0) {
                    let rho = corr_scen(
                        b.intra_bucket_corr
                            .get(i)
                            .and_then(|r| r.get(j))
                            .copied()
                            .unwrap_or(0.0),
                        scale,
                    );
                    q += T::from_f64(rho * rho) * cvr[i] * cvr[j];
                }
            }
        }
        (q.max_val(T::zero()).sqrt(), s_sum)
    };
    let (k_up, s_up) = eval_dir(&c_up);
    let (k_dn, s_dn) = eval_dir(&c_dn);
    if k_up.val() >= k_dn.val() {
        (k_up, s_up)
    } else {
        (k_dn, s_dn)
    }
}
fn sbm_scenario_charge<T: Real>(s: &FrtbSaSpec<T>, scale: f64) -> T {
    let classes = [
        FrtbRiskClass::Girr,
        FrtbRiskClass::CsrNonSec,
        FrtbRiskClass::Equity,
        FrtbRiskClass::Commodity,
        FrtbRiskClass::Fx,
    ];
    let mut total = T::zero();
    for rc in classes {
        let idx: Vec<usize> = s
            .buckets
            .iter()
            .enumerate()
            .filter_map(|(i, b)| (b.risk_class == rc).then_some(i))
            .collect();
        if idx.is_empty() {
            continue;
        }
        let mut kd = Vec::with_capacity(idx.len());
        let mut kv = Vec::with_capacity(idx.len());
        let mut kc = Vec::with_capacity(idx.len());
        for &bi in &idx {
            let b = &s.buckets[bi];
            kd.push(bucket_quad(
                &b.delta_sensitivities,
                &b.delta_risk_weights,
                &b.intra_bucket_corr,
                scale,
            ));
            kv.push(bucket_quad(
                &b.vega_sensitivities,
                &b.vega_risk_weights,
                &b.intra_bucket_corr,
                scale,
            ));
            kc.push(bucket_curv(b, scale));
        }
        for (pairs, is_curv) in [(&kd, false), (&kv, false), (&kc, true)] {
            let eval_agg = |use_cap: bool| -> T {
                let mut q = T::zero();
                for (u, &(k_u, s_u0)) in pairs.iter().enumerate() {
                    q += k_u * k_u;
                    let s_u = if use_cap {
                        s_u0.max_val(-k_u).min_val(k_u)
                    } else {
                        s_u0
                    };
                    for (v, &(k_v, s_v0)) in pairs.iter().enumerate() {
                        if u != v {
                            let s_v = if use_cap {
                                s_v0.max_val(-k_v).min_val(k_v)
                            } else {
                                s_v0
                            };
                            if is_curv && s_u.val() < 0.0 && s_v.val() < 0.0 {
                                continue;
                            }
                            let b_u = s.buckets[idx[u]].bucket_id;
                            let b_v = s.buckets[idx[v]].bucket_id;
                            let g0 = s
                                .inter_bucket_corr
                                .iter()
                                .find_map(|&(c_rc, id1, id2, g)| {
                                    (c_rc == rc
                                        && ((id1 == b_u && id2 == b_v)
                                            || (id1 == b_v && id2 == b_u)))
                                        .then_some(g)
                                })
                                .unwrap_or(0.0);
                            let g = corr_scen(g0, scale);
                            let g_eff = if is_curv { g * g } else { g };
                            q += T::from_f64(g_eff) * s_u * s_v;
                        }
                    }
                }
                q
            };
            let mut q = eval_agg(false);
            if !is_curv && q.val() < 0.0 {
                q = eval_agg(true);
            }
            total += q.max_val(T::zero()).sqrt();
        }
    }
    total
}

pub fn compute_frtb_sa<T: Real>(s: &FrtbSaSpec<T>) -> FrtbSaReport<T> {
    let sc = [
        sbm_scenario_charge(s, 0.75),
        sbm_scenario_charge(s, 1.0),
        sbm_scenario_charge(s, 1.25),
    ];
    let mut drc = T::zero();
    if !s.drc_obligors.is_empty() {
        let max_b = s
            .drc_obligors
            .iter()
            .map(|o| o.bucket_id)
            .max()
            .unwrap_or(0);
        for bucket in 0..=max_b {
            let mut gl = T::zero();
            let mut gs = T::zero();
            let mut wl = T::zero();
            let mut ws = T::zero();
            for o in s.drc_obligors.iter().filter(|o| o.bucket_id == bucket) {
                let wm = T::from_f64(o.maturity.clamp(0.25, 1.0));
                let rw = T::from_f64(o.default_risk_weight);
                if o.notional.val() >= 0.0 {
                    let jtd = (o.notional * o.lgd + o.mark_to_market_pnl).max_val(T::zero()) * wm;
                    gl += jtd;
                    wl += jtd * rw;
                } else {
                    let jtd =
                        ((-o.notional) * o.lgd - o.mark_to_market_pnl).max_val(T::zero()) * wm;
                    gs += jtd;
                    ws += jtd * rw;
                }
            }
            if (gl + gs).val() > 0.0 {
                let hbr = gl / (gl + gs).max_val(T::from_f64(1e-14));
                drc += (wl - hbr * ws).max_val(T::zero());
            }
        }
    }
    let rrao = s
        .exotic_residual_notionals
        .iter()
        .copied()
        .fold(T::zero(), |a, x| a + x.abs())
        * 0.01
        + s.other_residual_notionals
            .iter()
            .copied()
            .fold(T::zero(), |a, x| a + x.abs())
            * 0.001;
    let ch = sc[0].max_val(sc[1]).max_val(sc[2]);
    FrtbSaReport {
        sbm_low: sc[0],
        sbm_med: sc[1],
        sbm_high: sc[2],
        sbm_charge: ch,
        drc_charge: drc,
        rrao_charge: rrao,
        total_capital: ch + drc + rrao,
    }
}

pub fn compute_frtb_ima<T: Real>(s: &FrtbImaSpec<T>) -> FrtbImaReport<T> {
    let k_len = s.position_weights.len();
    let mults = [1.0, 1.0, 2.0, 2.0, 6.0];
    let mut h = [T::zero(); 5];
    let mut h_grad = vec![[0.0_f64; 5]; k_len];
    for i in 0..5 {
        let scen = &s.horizon_scenario_pnls[i];
        if scen.is_empty() {
            continue;
        }
        let mut p: Vec<(usize, T)> = scen
            .iter()
            .enumerate()
            .map(|(idx, r)| {
                let pnl = r
                    .iter()
                    .zip(&s.position_weights)
                    .map(|(x, w)| *x * *w)
                    .fold(T::zero(), |a, z| a + z);
                (idx, pnl)
            })
            .collect();
        p.sort_by(|a, b| a.1.val().total_cmp(&b.1.val()));
        let alpha_m = ((1.0 - s.confidence) * p.len() as f64).max(1e-12);
        let whole = (alpha_m.floor() as usize).min(p.len());
        let rem = alpha_m - whole as f64;
        let mut es_i = T::zero();
        for &(idx, pnl) in p.iter().take(whole) {
            es_i -= pnl / alpha_m;
            for (k, x) in scen[idx].iter().enumerate().take(k_len) {
                h_grad[k][i] -= x.val() / alpha_m;
            }
        }
        if whole < p.len() && rem > 1e-12 {
            let (idx, pnl) = p[whole];
            es_i -= pnl * (rem / alpha_m);
            for (k, x) in scen[idx].iter().enumerate().take(k_len) {
                h_grad[k][i] -= x.val() * (rem / alpha_m);
            }
        }
        h[i] = es_i;
    }
    let mut q = T::zero();
    for i in 0..5 {
        q += h[i] * h[i] * T::from_f64(mults[i]);
    }
    let liq = q.max_val(T::zero()).sqrt();
    let denom = liq.val().max(1e-14);
    let mut marginal_es = vec![0.0; k_len];
    let mut euler_component_es = vec![0.0; k_len];
    for k in 0..k_len {
        let mut g = 0.0;
        for i in 0..5 {
            g += mults[i] * h[i].val() * h_grad[k][i] / denom;
        }
        marginal_es[k] = g;
        euler_component_es[k] = s.position_weights[k].val() * g;
    }
    FrtbImaReport {
        horizon_es: h,
        liquidity_adjusted_es: liq,
        total_ima_capital: liq + s.ses_addon,
        marginal_es,
        euler_component_es,
    }
}

#[derive(Clone, Debug)]
pub struct FrtbImaSpec<T: Real> {
    pub position_weights: Vec<T>,
    pub horizon_scenario_pnls: [Vec<Vec<T>>; 5],
    pub confidence: f64,
    pub ses_addon: T,
}
#[derive(Clone, Debug)]
pub struct FrtbImaReport<T: Real> {
    pub horizon_es: [T; 5],
    pub liquidity_adjusted_es: T,
    pub total_ima_capital: T,
    pub marginal_es: Vec<f64>,
    pub euler_component_es: Vec<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct BcbsD424OutputFloorAndSmaOpRiskSpec<T: Real = f64> {
    pub irb_credit_rwa: T,
    pub ima_market_rwa: T,
    pub standardized_credit_rwa: T,
    pub standardized_market_rwa: T,
    pub standardized_cva_rwa: T,
    pub output_floor_multiplier: T,
    pub interest_lease_dividend_component: T,
    pub services_component: T,
    pub financial_component: T,
    pub average_annual_operational_loss: T,
    pub cet1_capital: T,
    pub at1_and_tier2_capital: T,
    pub eligible_tlac_senior_debt: T,
}
#[derive(Clone, Copy, Debug)]
pub struct BcbsD424OutputFloorAndSmaOpRiskReport<T: Real = f64> {
    pub business_indicator_bi: T,
    pub business_indicator_component_bic: T,
    pub internal_loss_multiplier_ilm: T,
    pub sma_operational_risk_capital: T,
    pub sma_operational_risk_rwa: T,
    pub unfloored_total_rwa: T,
    pub standardized_total_rwa: T,
    pub floored_total_rwa: T,
    pub output_floor_rwa_add_on: T,
    pub floored_cet1_ratio: T,
    pub tlac_to_floored_rwa_ratio: T,
}
pub fn evaluate_bcbs_d424_output_floor_and_sma_op_risk<T: Real>(
    s: &BcbsD424OutputFloorAndSmaOpRiskSpec<T>,
) -> BcbsD424OutputFloorAndSmaOpRiskReport<T> {
    let bi = s.interest_lease_dividend_component + s.services_component + s.financial_component;
    let smooth_min = |a: T, b: T| {
        if (a - b).val().abs() > 1e-4 {
            a.min_val(b)
        } else {
            (a + b - ((a - b) * (a - b) + T::from_f64(1e-8)).sqrt()) / T::from_f64(2.0)
        }
    };
    let smooth_max = |a: T, b: T| {
        if (a - b).val().abs() > 1e-4 {
            a.max_val(b)
        } else {
            (a + b + ((a - b) * (a - b) + T::from_f64(1e-8)).sqrt()) / T::from_f64(2.0)
        }
    };
    let b1 = smooth_min(bi, T::one());
    let b2 = smooth_min(smooth_max(bi - T::one(), T::zero()), T::from_f64(29.0));
    let b3 = smooth_max(bi - T::from_f64(30.0), T::zero());
    let bic = b1 * T::from_f64(0.12) + b2 * T::from_f64(0.15) + b3 * T::from_f64(0.18);
    let loss_cap = s.average_annual_operational_loss * T::from_f64(15.0);
    let ratio = (loss_cap / bic.max_val(T::from_f64(1e-10))).max_val(T::from_f64(1e-10));
    let ilm = (T::from_f64(std::f64::consts::E - 1.0) + ratio.powf(0.8)).ln();
    let capital = bic * ilm;
    let op_rwa = capital * 12.5;
    let unfloored = s.irb_credit_rwa + s.ima_market_rwa + op_rwa;
    let standardized =
        s.standardized_credit_rwa + s.standardized_market_rwa + s.standardized_cva_rwa + op_rwa;
    let floor = standardized * s.output_floor_multiplier.max_val(T::zero());
    let floored = unfloored.max_val(floor);
    BcbsD424OutputFloorAndSmaOpRiskReport {
        business_indicator_bi: bi,
        business_indicator_component_bic: bic,
        internal_loss_multiplier_ilm: ilm,
        sma_operational_risk_capital: capital,
        sma_operational_risk_rwa: op_rwa,
        unfloored_total_rwa: unfloored,
        standardized_total_rwa: standardized,
        floored_total_rwa: floored,
        output_floor_rwa_add_on: (floored - unfloored).max_val(T::zero()),
        floored_cet1_ratio: s.cet1_capital / floored.max_val(T::from_f64(1e-10)),
        tlac_to_floored_rwa_ratio: (s.cet1_capital
            + s.at1_and_tier2_capital
            + s.eligible_tlac_senior_debt)
            / floored.max_val(T::from_f64(1e-10)),
    }
}

#[cfg(test)]
mod phase29_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn bcbs_invariants_and_dual() {
        fn mk<T: Real>(x: T) -> BcbsD424OutputFloorAndSmaOpRiskSpec<T> {
            BcbsD424OutputFloorAndSmaOpRiskSpec {
                irb_credit_rwa: T::from_f64(100.0) + x,
                ima_market_rwa: T::from_f64(50.0),
                standardized_credit_rwa: T::from_f64(120.0),
                standardized_market_rwa: T::from_f64(60.0),
                standardized_cva_rwa: T::from_f64(20.0),
                output_floor_multiplier: T::from_f64(0.725),
                interest_lease_dividend_component: T::from_f64(20.0),
                services_component: T::from_f64(15.0),
                financial_component: T::from_f64(30.0),
                average_annual_operational_loss: T::from_f64(2.0),
                cet1_capital: T::from_f64(20.0),
                at1_and_tier2_capital: T::from_f64(10.0),
                eligible_tlac_senior_debt: T::from_f64(40.0),
            }
        }
        let r = evaluate_bcbs_d424_output_floor_and_sma_op_risk(&mk(0.0));
        assert!(r.floored_total_rwa >= r.unfloored_total_rwa);
        let h = 1e-5;
        let fd = (evaluate_bcbs_d424_output_floor_and_sma_op_risk(&mk(h)).unfloored_total_rwa
            - evaluate_bcbs_d424_output_floor_and_sma_op_risk(&mk(-h)).unfloored_total_rwa)
            / (2.0 * h);
        let ad = evaluate_bcbs_d424_output_floor_and_sma_op_risk(&mk(Dual::<1>::variable(0.0, 0)))
            .unfloored_total_rwa;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        let mut ilms = mk(0.0);
        ilms.average_annual_operational_loss = 10.77 / 15.0;
        let r0 = evaluate_bcbs_d424_output_floor_and_sma_op_risk(&ilms);
        assert!((r0.internal_loss_multiplier_ilm - 1.0).abs() < 1e-12);
    }
}
