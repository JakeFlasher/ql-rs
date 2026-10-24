use ql_ad::Real;
use std::collections::HashMap;

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
    (common * common + idio).max_val(T::zero()).sqrt()
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
            SaCcrAssetClass::InterestRate { ccy, end_time, .. } => {
                let bucket = if *end_time < 1.0 {
                    0
                } else if *end_time <= 5.0 {
                    1
                } else {
                    2
                };
                ir.entry(ccy.clone()).or_insert([T::zero(); 3])[bucket] += d * sd(0.0, *end_time);
            }
            SaCcrAssetClass::Fx { pair } => *fx.entry(pair.clone()).or_insert(T::zero()) += d,
            SaCcrAssetClass::Credit {
                entity,
                is_index,
                supervisory_factor,
            } => {
                let rho = if *is_index { 0.8 } else { 0.5 };
                cr.entry(entity.clone())
                    .and_modify(|x| x.0 += d)
                    .or_insert((d, *supervisory_factor, rho));
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
        add_ir += (b[0] * b[0]
            + b[1] * b[1]
            + b[2] * b[2]
            + T::from_f64(1.4) * b[0] * b[1]
            + T::from_f64(1.4) * b[1] * b[2]
            + T::from_f64(0.6) * b[0] * b[2])
            .max_val(T::zero())
            .sqrt()
            * 0.005;
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
            (sum * 0.4 * sum * 0.4
                + xs.iter()
                    .map(|x| *x * *x * 0.84)
                    .fold(T::zero(), |a, x| a + x))
            .sqrt()
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
