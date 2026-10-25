use ql_ad::Real;
use ql_math::inv_norm_cdf;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct GreenSaCcrLifetimeKvaMvaSpec<T: Real = f64> {
    pub time_grid: Vec<f64>,
    pub trade_notional: T,
    pub trade_maturity: f64,
    pub initial_mtm: T,
    pub mtm_volatility: T,
    pub collateral_threshold: T,
    pub initial_margin_ratio: T,
    pub counterparty_pd: T,
    pub counterparty_lgd: T,
    pub cva_risk_weight: T,
    pub capital_ratio_requirement: T,
    pub hurdle_rate_gamma_k: T,
    pub funding_spread_sf: T,
    pub mva_spread_sm: T,
    pub risk_free_rate: T,
}
#[derive(Clone, Debug)]
pub struct GreenSaCcrLifetimeKvaMvaReport<T: Real = f64> {
    pub ead_profile: Vec<T>,
    pub k_ccr_profile: Vec<T>,
    pub k_cva_profile: Vec<T>,
    pub total_capital_profile: Vec<T>,
    pub lifetime_kva: T,
    pub capital_fva_offset: T,
    pub net_kva_after_funding_benefit: T,
    pub lifetime_mva: T,
}
pub fn compute_green_saccr_lifetime_kva_mva<T: Real>(
    s: &GreenSaCcrLifetimeKvaMvaSpec<T>,
) -> GreenSaCcrLifetimeKvaMvaReport<T> {
    let mut ead = Vec::new();
    let mut kccr = Vec::new();
    let mut kcva = Vec::new();
    let mut kt = Vec::new();
    let mut kva = T::zero();
    let mut off = T::zero();
    let mut mva = T::zero();
    let pd = s.counterparty_pd.max_val(T::from_f64(1e-10));
    let lgd = s.counterparty_lgd.max_val(T::from_f64(1e-8));
    let lambda = -(T::one() - pd.min_val(T::from_f64(0.999))).ln();
    let mut zpd = T::from_f64(inv_norm_cdf(pd.val()));
    for _ in 0..3 {
        zpd -= (zpd.norm_cdf() - pd) / zpd.norm_pdf().max_val(T::from_f64(1e-12));
    }
    for (idx, &t) in s.time_grid.iter().enumerate() {
        if t >= s.trade_maturity {
            ead.push(T::zero());
            kccr.push(T::zero());
            kcva.push(T::zero());
            kt.push(T::zero());
            continue;
        }
        let tau = (s.trade_maturity - t).max(1e-4);
        let sd = T::from_f64((1.0 - (-0.05 * tau).exp()) / 0.05);
        let addon = s.trade_notional * sd * 0.005;
        let sig = s.mtm_volatility
            * T::from_f64((t * tau / s.trade_maturity.max(1e-8)).max(0.0).sqrt())
            + T::from_f64(1e-6);
        let mu = s.initial_mtm * T::from_f64(tau / s.trade_maturity);
        let z = (mu - s.collateral_threshold) / sig;
        let rc = (mu - s.collateral_threshold) * z.norm_cdf() + sig * z.norm_pdf();
        let mult = (T::from_f64(0.05)
            + T::from_f64(0.95)
                * ((mu - s.collateral_threshold).min_val(T::zero())
                    / (addon.abs().max_val(T::from_f64(1e-8)) * 1.9))
                    .exp())
        .min_val(T::one());
        let e = (rc + mult * addon) * 1.4;
        let _rho = T::from_f64(0.12) * (T::one() - (-pd * 50.0).exp())
            / T::from_f64(1.0 - (-50.0f64).exp())
            + T::from_f64(0.24)
                * (T::one()
                    - (T::one() - (-pd * 50.0).exp()) / T::from_f64(1.0 - (-50.0f64).exp()));
        let mm = T::from_f64(tau.clamp(1.0, 5.0));
        let b = (T::from_f64(0.11852) - pd.ln() * 0.05478).powf(2.0);
        let rho = _rho;
        let q = (zpd + rho.sqrt() * T::from_f64(3.090232306))
            / ((T::one() - rho).max_val(T::from_f64(1e-8)).sqrt());
        let w = lgd * (q.norm_cdf() - pd) * (T::one() + (mm - 2.5) * b)
            / (T::one() - b * 1.5).max_val(T::from_f64(1e-8))
            * 1.06;
        let kc = s.capital_ratio_requirement * w * e * 12.5;
        let kv = s.capital_ratio_requirement * s.cva_risk_weight * sd * e * 0.65 * 12.5 / 1.4;
        let total = kc + kv;
        ead.push(e);
        kccr.push(kc);
        kcva.push(kv);
        kt.push(total);
        let next_t = s
            .time_grid
            .get(idx + 1)
            .copied()
            .unwrap_or(s.trade_maturity)
            .min(s.trade_maturity);
        let dt = (next_t - t).max(0.0);
        let disc = (-(s.risk_free_rate + lambda) * t).exp() * dt;
        kva += disc * s.hurdle_rate_gamma_k * total;
        off -= disc * s.funding_spread_sf * total;
        mva += disc * s.mva_spread_sm * s.initial_margin_ratio * addon;
    }
    GreenSaCcrLifetimeKvaMvaReport {
        ead_profile: ead,
        k_ccr_profile: kccr,
        k_cva_profile: kcva,
        total_capital_profile: kt,
        lifetime_kva: kva,
        capital_fva_offset: off,
        net_kva_after_funding_benefit: kva + off,
        lifetime_mva: mva,
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum SaCcrAssetClass {
    InterestRate {
        ccy: String,
        start_time: f64,
        end_time: f64,
    },
    Fx {
        pair: String,
    },
    Credit {
        entity: String,
        is_index: bool,
        supervisory_factor: f64,
    },
    Equity {
        entity: String,
        is_index: bool,
    },
    Commodity {
        hedging_set: String,
        sub_class: String,
        supervisory_factor: f64,
    },
}
#[derive(Copy, Clone, Debug)]
pub struct SaCcrOptionSpec<T: Real> {
    pub underlying_price: T,
    pub strike: T,
    pub expiry: f64,
    pub supervisory_vol: f64,
    pub shift: f64,
    pub is_call: bool,
    pub is_bought: bool,
}
#[derive(Clone, Debug)]
pub struct SaCcrTrade<T: Real> {
    pub mtm: T,
    pub notional: T,
    pub maturity: f64,
    pub asset_class: SaCcrAssetClass,
    pub is_long_primary: bool,
    pub option_spec: Option<SaCcrOptionSpec<T>>,
}
#[derive(Clone, Debug)]
pub struct SaCcrNettingSet<T: Real> {
    pub trades: Vec<SaCcrTrade<T>>,
    pub collateral: T,
    pub is_margined: bool,
    pub threshold: T,
    pub mta: T,
    pub nica: T,
    pub mpor_years: f64,
}
#[derive(Clone, Debug)]
pub struct SaCcrReport<T: Real> {
    pub rc: T,
    pub add_on_ir: T,
    pub add_on_fx: T,
    pub add_on_credit: T,
    pub add_on_equity: T,
    pub add_on_commodity: T,
    pub add_on_aggregate: T,
    pub multiplier: T,
    pub pfe: T,
    pub ead: T,
}

fn sd(start: f64, end: f64) -> f64 {
    let s = start.max(0.0);
    let e = end.max(s + 10.0 / 250.0);
    ((-0.05 * s).exp() - (-0.05 * e).exp()) / 0.05
}
fn signed_delta<T: Real>(trade: &SaCcrTrade<T>) -> T {
    if let Some(o) = trade.option_spec {
        let phi = if o.is_call { 1.0 } else { -1.0 };
        let omega = if o.is_bought { 1.0 } else { -1.0 };
        let t = o.expiry.max(10.0 / 250.0);
        let st = o.supervisory_vol * t.sqrt();
        if st <= 1e-14 {
            return T::from_f64(omega * phi);
        }
        let f = o.underlying_price + o.shift;
        let k = o.strike + o.shift;
        let d = ((f / k).ln() + 0.5 * st * st) / st;
        T::from_f64(omega * phi) * (T::from_f64(phi) * d).norm_cdf()
    } else if trade.is_long_primary {
        T::one()
    } else {
        -T::one()
    }
}
fn mf<T: Real>(trade: &SaCcrTrade<T>, ns: &SaCcrNettingSet<T>) -> T {
    if ns.is_margined {
        T::from_f64(1.5 * ns.mpor_years.max(10.0 / 250.0).sqrt())
    } else {
        T::from_f64(trade.maturity.clamp(10.0 / 250.0, 1.0).sqrt())
    }
}
fn one_factor<T: Real>(items: Vec<(T, f64, f64)>) -> T {
    let mut common = T::zero();
    let mut idio = T::zero();
    for (x, sf, rho) in items {
        let a = x * sf;
        common += a * rho;
        idio += a * a * (1.0 - rho * rho);
    }
    let q = common * common + idio;
    if q.val() <= 0.0 { T::zero() } else { q.sqrt() }
}

pub fn compute_sa_ccr_ead<T: Real>(ns: &SaCcrNettingSet<T>) -> SaCcrReport<T> {
    let v = ns.trades.iter().fold(T::zero(), |a, x| a + x.mtm);
    let collateral = ns.collateral;
    let rc = if ns.is_margined {
        (v - collateral)
            .max_val(ns.threshold + ns.mta - ns.nica)
            .max_val(T::zero())
    } else {
        (v - collateral).max_val(T::zero())
    };
    let mut ir: HashMap<String, [T; 3]> = HashMap::new();
    let mut fx: HashMap<String, T> = HashMap::new();
    let mut cr: HashMap<String, (T, f64, f64)> = HashMap::new();
    let mut eq: HashMap<String, (T, f64, f64)> = HashMap::new();
    let mut com: HashMap<(String, String), T> = HashMap::new();
    for tr in &ns.trades {
        let d = signed_delta(tr) * tr.notional * mf(tr, ns);
        match &tr.asset_class {
            SaCcrAssetClass::InterestRate {
                ccy,
                start_time,
                end_time,
            } => {
                let bucket = if *end_time < 1.0 {
                    0
                } else if *end_time <= 5.0 {
                    1
                } else {
                    2
                };
                ir.entry(ccy.clone()).or_insert([T::zero(); 3])[bucket] +=
                    d * sd(*start_time, *end_time);
            }
            SaCcrAssetClass::Fx { pair } => *fx.entry(pair.clone()).or_insert(T::zero()) += d,
            SaCcrAssetClass::Credit {
                entity,
                is_index,
                supervisory_factor,
            } => {
                let rho = if *is_index { 0.8 } else { 0.5 };
                cr.entry(entity.clone())
                    .and_modify(|x| x.0 += d * sd(0.0, tr.maturity))
                    .or_insert((d * sd(0.0, tr.maturity), *supervisory_factor, rho));
            }
            SaCcrAssetClass::Equity { entity, is_index } => {
                let (sf, rho) = if *is_index { (0.20, 0.8) } else { (0.32, 0.5) };
                eq.entry(entity.clone())
                    .and_modify(|x| x.0 += d)
                    .or_insert((d, sf, rho));
            }
            SaCcrAssetClass::Commodity {
                hedging_set,
                sub_class,
                supervisory_factor,
            } => {
                *com.entry((hedging_set.clone(), sub_class.clone()))
                    .or_insert(T::zero()) += d * *supervisory_factor;
            }
        }
    }
    let mut add_ir = T::zero();
    for b in ir.values() {
        let q = b[0] * b[0]
            + b[1] * b[1]
            + b[2] * b[2]
            + T::from_f64(1.4) * b[0] * b[1]
            + T::from_f64(1.4) * b[1] * b[2]
            + T::from_f64(0.6) * b[0] * b[2];
        if q.val() > 0.0 {
            add_ir += q.sqrt() * 0.005;
        }
    }
    let add_fx = fx.values().fold(T::zero(), |a, x| a + x.abs() * 0.04);
    let add_cr = one_factor(cr.values().map(|(x, s, r)| (*x, *s, *r)).collect());
    let add_eq = one_factor(eq.values().map(|(x, s, r)| (*x, *s, *r)).collect());
    let mut by_set: HashMap<String, Vec<T>> = HashMap::new();
    for ((h, _), x) in com {
        by_set.entry(h).or_default().push(x);
    }
    let add_com = by_set
        .values()
        .map(|xs| {
            let sum = xs.iter().copied().fold(T::zero(), |a, x| a + x);
            let q = sum * 0.4 * sum * 0.4
                + xs.iter()
                    .map(|x| *x * *x * 0.84)
                    .fold(T::zero(), |a, x| a + x);
            if q.val() > 0.0 { q.sqrt() } else { T::zero() }
        })
        .fold(T::zero(), |a, x| a + x);
    let agg = add_ir + add_fx + add_cr + add_eq + add_com;
    let multiplier = if (v - collateral).val() >= 0.0 || agg.val() <= 1e-15 {
        T::one()
    } else {
        T::one()
            .min_val(T::from_f64(0.05) + T::from_f64(0.95) * ((v - collateral) / (agg * 1.9)).exp())
    };
    let pfe = multiplier * agg;
    let ead = T::from_f64(1.4) * (rc + pfe);
    SaCcrReport {
        rc,
        add_on_ir: add_ir,
        add_on_fx: add_fx,
        add_on_credit: add_cr,
        add_on_equity: add_eq,
        add_on_commodity: add_com,
        add_on_aggregate: agg,
        multiplier,
        pfe,
        ead,
    }
}

#[derive(Clone, Debug)]
pub struct BaCvaSingleNameHedge<T: Real = f64> {
    pub notional: T,
    pub maturity: T,
    pub hedge_correlation: f64,
}

#[derive(Clone, Debug)]
pub struct BaCvaCounterpartySpec<T: Real = f64> {
    pub risk_weight: T,
    pub effective_maturity: T,
    pub ead: T,
    pub single_name_hedges: Vec<BaCvaSingleNameHedge<T>>,
}

#[derive(Clone, Debug)]
pub struct BaCvaIndexHedge<T: Real = f64> {
    pub risk_weight: T,
    pub maturity: T,
    pub notional: T,
}

#[derive(Clone, Debug)]
pub struct BaCvaPortfolioSpec<T: Real = f64> {
    pub counterparties: Vec<BaCvaCounterpartySpec<T>>,
    pub index_hedges: Vec<BaCvaIndexHedge<T>>,
    pub hurdle_rate: T,
    pub risk_free_rate: T,
}

#[derive(Clone, Debug)]
pub struct BaCvaReport<T: Real> {
    pub scva_by_counterparty: Vec<T>,
    pub k_reduced: T,
    pub k_hedged: T,
    pub k_full: T,
    pub capital_charge: T,
    pub kva_charge: T,
}

fn ba_cva_supervisory_df<T: Real>(maturity: T) -> T {
    let m = maturity.max_val(T::zero());
    let x = m * 0.05;
    if x.val().abs() < 1e-2 {
        T::one() - x * 0.5 + x * x / 6.0 - x * x * x / 24.0
    } else {
        (T::one() - (-x).exp()) / x
    }
}

pub fn compute_ba_cva_capital<T: Real>(spec: &BaCvaPortfolioSpec<T>) -> BaCvaReport<T> {
    if spec.counterparties.is_empty() {
        return BaCvaReport {
            scva_by_counterparty: Vec::new(),
            k_reduced: T::zero(),
            k_hedged: T::zero(),
            k_full: T::zero(),
            capital_charge: T::zero(),
            kva_charge: T::zero(),
        };
    }
    let rho = 0.50;
    let one_minus_rho2 = 1.0 - rho * rho;
    let inv_alpha = 1.0 / 1.4;
    let mut scvas = Vec::with_capacity(spec.counterparties.len());
    let mut sum_scva = T::zero();
    let mut sum_scva2 = T::zero();
    let mut sum_net = T::zero();
    let mut sum_net2 = T::zero();
    let mut sum_hma = T::zero();
    let mut avg_m = T::zero();

    for cpty in &spec.counterparties {
        let df_c = ba_cva_supervisory_df(cpty.effective_maturity);
        let scva = cpty.risk_weight * cpty.effective_maturity * cpty.ead * df_c * inv_alpha;
        scvas.push(scva);
        sum_scva += scva;
        sum_scva2 += scva * scva;
        avg_m += cpty.effective_maturity;

        let mut snh = T::zero();
        let mut hma = T::zero();
        for h in &cpty.single_name_hedges {
            let df_h = ba_cva_supervisory_df(h.maturity);
            let h_term = cpty.risk_weight * h.maturity * h.notional * df_h;
            let r_hc = h.hedge_correlation.clamp(0.0, 1.0);
            snh += h_term * r_hc;
            hma += h_term * h_term * (1.0 - r_hc * r_hc);
        }
        let net_c = scva - snh;
        sum_net += net_c;
        sum_net2 += net_c * net_c;
        sum_hma += hma;
    }
    avg_m /= T::from_f64(spec.counterparties.len() as f64);

    let mut ih = T::zero();
    for ind in &spec.index_hedges {
        let df_i = ba_cva_supervisory_df(ind.maturity);
        ih += ind.risk_weight * ind.maturity * ind.notional * df_i;
    }

    let sys_red = sum_scva * rho;
    let k_reduced = (sys_red * sys_red + sum_scva2 * one_minus_rho2)
        .max_val(T::from_f64(1e-16))
        .sqrt();

    let sys_hedged = sum_net * rho - ih;
    let k_hedged = (sys_hedged * sys_hedged + sum_net2 * one_minus_rho2 + sum_hma)
        .max_val(T::from_f64(1e-16))
        .sqrt();

    let k_full = k_reduced * 0.25 + k_hedged * 0.75;
    let capital_charge = k_full * 0.65;
    let spread_cost = (spec.hurdle_rate - spec.risk_free_rate).max_val(T::zero());
    let kva_df = ba_cva_supervisory_df(avg_m);
    let kva_charge = capital_charge * spread_cost * avg_m * kva_df;

    BaCvaReport {
        scva_by_counterparty: scvas,
        k_reduced,
        k_hedged,
        k_full,
        capital_charge,
        kva_charge,
    }
}

#[derive(Clone, Debug)]
pub struct GreenFrtbCvaHedgeCapitalSpec<T: Real = f64> {
    pub counterparty_sacr_eads: Vec<T>,
    pub counterparty_effective_maturities: Vec<T>,
    pub counterparty_risk_weights: Vec<T>,
    pub single_name_hedge_notionals: Vec<T>,
    pub single_name_hedge_maturities: Vec<T>,
    pub hedge_reference_correlations_rhc: Vec<T>,
    pub index_hedge_notional: T,
    pub index_hedge_maturity: T,
    pub index_risk_weight: T,
    pub sa_cva_spread_deltas: Vec<T>,
    pub sa_cva_vega_sensitivities: Vec<T>,
    pub hurdle_rate_gamma_k: T,
    pub hedge_margin_funding_spread: T,
    pub discount_rate: T,
}
#[derive(Clone, Debug)]
pub struct GreenFrtbCvaHedgeCapitalReport<T: Real = f64> {
    pub ba_cva_k_reduced_unhedged: T,
    pub ba_cva_k_hedged: T,
    pub ba_cva_misalignment_basis_charge: T,
    pub ba_cva_total_capital_charge: T,
    pub sa_cva_delta_capital: T,
    pub sa_cva_vega_capital: T,
    pub sa_cva_total_capital_charge: T,
    pub unhedged_lifetime_kva: T,
    pub hedged_lifetime_kva: T,
    pub net_kva_savings_after_hedge_mva: T,
}
pub fn evaluate_green_frtb_cva_hedge_relief_and_kva<T: Real>(
    s: &GreenFrtbCvaHedgeCapitalSpec<T>,
) -> GreenFrtbCvaHedgeCapitalReport<T> {
    let n = s
        .counterparty_sacr_eads
        .len()
        .min(s.counterparty_effective_maturities.len())
        .min(s.counterparty_risk_weights.len());
    let rho = T::from_f64(0.5);
    let mut sum_scva = T::zero();
    let mut sum_scva_sq = T::zero();
    for i in 0..n {
        let m = s.counterparty_effective_maturities[i].max_val(T::from_f64(1e-8));
        let df_m = (T::one() - (-T::from_f64(0.05) * m).exp()) / (T::from_f64(0.05) * m);
        let scva = s.counterparty_risk_weights[i] * m * s.counterparty_sacr_eads[i] * df_m;
        sum_scva += scva;
        sum_scva_sq += scva * scva;
    }
    let kr = (rho * rho * sum_scva * sum_scva + (T::one() - rho * rho) * sum_scva_sq).sqrt();
    let mut sum_snh = T::zero();
    let mut sum_snh_sq = T::zero();
    let mut hma = T::zero();
    for i in 0..s.single_name_hedge_notionals.len() {
        let rhc = s
            .hedge_reference_correlations_rhc
            .get(i)
            .copied()
            .unwrap_or_else(T::zero)
            .abs()
            .min_val(T::one());
        let mh = s
            .single_name_hedge_maturities
            .get(i)
            .copied()
            .unwrap_or_else(T::one)
            .max_val(T::from_f64(1e-8));
        let rw = s
            .counterparty_risk_weights
            .get(i)
            .copied()
            .unwrap_or_else(|| T::from_f64(0.02));
        let df_h = (T::one() - (-T::from_f64(0.05) * mh).exp()) / (T::from_f64(0.05) * mh);
        let snh = rw * mh * s.single_name_hedge_notionals[i].abs() * df_h;
        sum_snh += rhc * snh;
        sum_snh_sq += rhc * rhc * snh * snh;
        hma += (T::one() - rhc * rhc) * snh * snh;
    }
    let mi = s.index_hedge_maturity.max_val(T::from_f64(1e-8));
    let df_i = (T::one() - (-T::from_f64(0.05) * mi).exp()) / (T::from_f64(0.05) * mi);
    let ih = s.index_risk_weight.abs() * mi * s.index_hedge_notional.abs() * df_i;
    let sys_diff = rho * (sum_scva - sum_snh) - ih;
    let idio_sq = (T::one() - rho * rho) * (sum_scva_sq - sum_snh_sq).max_val(T::zero());
    let kh = (sys_diff * sys_diff + idio_sq + hma).sqrt();
    let mis = hma.sqrt();
    let beta = T::from_f64(0.25);
    let total = T::from_f64(0.65) * (beta * kr + (T::one() - beta) * kh);
    let mut sd = T::zero();
    for x in &s.sa_cva_spread_deltas {
        sd += *x * *x;
    }
    let mut sv = T::zero();
    for x in &s.sa_cva_vega_sensitivities {
        sv += *x * *x;
    }
    let sd = sd.sqrt();
    let sv = sv.sqrt();
    let sa = sd + sv;
    let un =
        kr * s.hurdle_rate_gamma_k.max_val(T::zero()) / s.discount_rate.max_val(T::from_f64(1e-8));
    let hed =
        kh * s.hurdle_rate_gamma_k.max_val(T::zero()) / s.discount_rate.max_val(T::from_f64(1e-8));
    let mva = (mis + ih) * s.hedge_margin_funding_spread.max_val(T::zero());
    GreenFrtbCvaHedgeCapitalReport {
        ba_cva_k_reduced_unhedged: kr,
        ba_cva_k_hedged: kh,
        ba_cva_misalignment_basis_charge: mis,
        ba_cva_total_capital_charge: total,
        sa_cva_delta_capital: sd,
        sa_cva_vega_capital: sv,
        sa_cva_total_capital_charge: sa,
        unhedged_lifetime_kva: un,
        hedged_lifetime_kva: hed,
        net_kva_savings_after_hedge_mva: (un - hed - mva).max_val(T::zero()),
    }
}

#[cfg(test)]
mod phase17_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn ba_cva_single_name_identity_hedge_reduction_and_ad() {
        let unhedged = BaCvaPortfolioSpec {
            counterparties: vec![BaCvaCounterpartySpec {
                risk_weight: 0.05,
                effective_maturity: 3.0,
                ead: 100.0,
                single_name_hedges: vec![],
            }],
            index_hedges: vec![],
            hurdle_rate: 0.10,
            risk_free_rate: 0.03,
        };
        let rep0 = compute_ba_cva_capital(&unhedged);
        assert!((rep0.k_reduced - rep0.scva_by_counterparty[0]).abs() < 1e-12);
        assert!((rep0.k_full - rep0.k_reduced).abs() < 1e-12);

        let mut hedged = unhedged.clone();
        hedged.counterparties[0]
            .single_name_hedges
            .push(BaCvaSingleNameHedge {
                notional: 40.0,
                maturity: 3.0,
                hedge_correlation: 1.0,
            });
        let rep_h = compute_ba_cva_capital(&hedged);
        assert!(rep_h.capital_charge < rep0.capital_charge);

        let h = 1e-5;
        let eval_ead = |e: f64| {
            let mut sp = hedged.clone();
            sp.counterparties[0].ead = e;
            compute_ba_cva_capital(&sp).capital_charge
        };
        let fd = (eval_ead(100.0 + h) - eval_ead(100.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sp_d = BaCvaPortfolioSpec {
            counterparties: vec![BaCvaCounterpartySpec {
                risk_weight: D::from_value(0.05),
                effective_maturity: D::from_value(3.0),
                ead: D::variable(100.0, 0),
                single_name_hedges: vec![BaCvaSingleNameHedge {
                    notional: D::from_value(40.0),
                    maturity: D::from_value(3.0),
                    hedge_correlation: 1.0,
                }],
            }],
            index_hedges: vec![],
            hurdle_rate: D::from_value(0.10),
            risk_free_rate: D::from_value(0.03),
        };
        let ad = compute_ba_cva_capital(&sp_d).capital_charge;
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase21_tests {
    use super::*;

    #[test]
    fn green_lifetime_profiles_are_finite_and_nonnegative() {
        let s = GreenSaCcrLifetimeKvaMvaSpec {
            time_grid: vec![0.0, 1.0, 2.0],
            trade_notional: 1e6,
            trade_maturity: 2.5,
            initial_mtm: 1e4,
            mtm_volatility: 0.1,
            collateral_threshold: 0.0,
            initial_margin_ratio: 0.1,
            counterparty_pd: 0.02,
            counterparty_lgd: 0.6,
            cva_risk_weight: 0.5,
            capital_ratio_requirement: 0.08,
            hurdle_rate_gamma_k: 0.1,
            funding_spread_sf: 0.01,
            mva_spread_sm: 0.01,
            risk_free_rate: 0.03,
        };
        let r = compute_green_saccr_lifetime_kva_mva(&s);
        assert!(r.lifetime_kva.val().is_finite());
        assert!(r.ead_profile.iter().all(|x| x.val() >= 0.0));
    }
}

#[cfg(test)]
mod phase26_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(ead: T) -> GreenFrtbCvaHedgeCapitalSpec<T> {
        GreenFrtbCvaHedgeCapitalSpec {
            counterparty_sacr_eads: vec![ead, T::from_f64(80.0)],
            counterparty_effective_maturities: vec![T::from_f64(3.0), T::from_f64(5.0)],
            counterparty_risk_weights: vec![T::from_f64(0.02), T::from_f64(0.03)],
            single_name_hedge_notionals: vec![],
            single_name_hedge_maturities: vec![],
            hedge_reference_correlations_rhc: vec![],
            index_hedge_notional: T::zero(),
            index_hedge_maturity: T::from_f64(5.0),
            index_risk_weight: T::from_f64(0.02),
            sa_cva_spread_deltas: vec![T::one()],
            sa_cva_vega_sensitivities: vec![T::from_f64(0.5)],
            hurdle_rate_gamma_k: T::from_f64(0.1),
            hedge_margin_funding_spread: T::from_f64(0.01),
            discount_rate: T::from_f64(0.03),
        }
    }
    #[test]
    fn frtb_cva_invariants_and_dual_fd() {
        let r = evaluate_green_frtb_cva_hedge_relief_and_kva(&spec(100.0));
        assert!((r.ba_cva_total_capital_charge - 0.65 * r.ba_cva_k_reduced_unhedged).abs() < 1e-12);
        let h = 1e-5;
        let fd = (evaluate_green_frtb_cva_hedge_relief_and_kva(&spec(100.0 + h))
            .ba_cva_k_reduced_unhedged
            - evaluate_green_frtb_cva_hedge_relief_and_kva(&spec(100.0 - h))
                .ba_cva_k_reduced_unhedged)
            / (2.0 * h);
        let ad = evaluate_green_frtb_cva_hedge_relief_and_kva(&spec(Dual::<1>::variable(100.0, 0)))
            .ba_cva_k_reduced_unhedged;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
