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
fn charge<T: Real>(b: &FrtbSensitivityBucket<T>, scale: f64) -> T {
    let mut q = T::zero();
    for i in 0..b.delta_sensitivities.len() {
        for j in 0..b.delta_sensitivities.len() {
            let c = if i == j {
                1.
            } else {
                (b.intra_bucket_corr
                    .get(i)
                    .and_then(|r| r.get(j))
                    .copied()
                    .unwrap_or(0.)
                    * scale)
                    .clamp(-1., 1.)
            };
            q += T::from_f64(c)
                * b.delta_sensitivities[i]
                * T::from_f64(b.delta_risk_weights[i])
                * b.delta_sensitivities[j]
                * T::from_f64(b.delta_risk_weights[j]);
        }
    }
    q.max_val(T::zero()).sqrt()
}
pub fn compute_frtb_sa<T: Real>(s: &FrtbSaSpec<T>) -> FrtbSaReport<T> {
    let mut sc = [T::zero(); 3];
    for b in &s.buckets {
        sc[0] += charge(b, 0.75);
        sc[1] += charge(b, 1.);
        sc[2] += charge(b, 1.25);
    }
    let mut drc = T::zero();
    for bucket in 0..=s
        .drc_obligors
        .iter()
        .map(|o| o.bucket_id)
        .max()
        .unwrap_or(0)
    {
        let mut p = T::zero();
        let mut n = T::zero();
        for o in s.drc_obligors.iter().filter(|o| o.bucket_id == bucket) {
            let wm = T::from_f64(o.maturity.clamp(0.25, 1.));
            let j = if o.notional.val() >= 0. {
                (o.notional * o.lgd + o.mark_to_market_pnl).max_val(T::zero()) * wm
            } else {
                ((-o.notional) * o.lgd - o.mark_to_market_pnl).max_val(T::zero()) * wm
            } * T::from_f64(o.default_risk_weight);
            if o.notional.val() >= 0. {
                p += j
            } else {
                n += j
            }
        }
        let h = p / (p + n).max_val(T::from_f64(1e-14));
        drc += (p - h * n).max_val(T::zero());
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
fn es<T: Real>(x: &[T], conf: f64) -> T {
    let mut v = x.to_vec();
    v.sort_by(|a, b| a.val().total_cmp(&b.val()));
    let n = ((1. - conf) * v.len() as f64).ceil().max(1.) as usize;
    v.iter().take(n).copied().fold(T::zero(), |a, z| a - z) / n as f64
}
pub fn compute_frtb_ima<T: Real>(s: &FrtbImaSpec<T>) -> FrtbImaReport<T> {
    let mut h = [T::zero(); 5];
    for i in 0..5 {
        let p = s.horizon_scenario_pnls[i]
            .iter()
            .map(|r| {
                r.iter()
                    .zip(&s.position_weights)
                    .map(|(x, w)| *x * *w)
                    .fold(T::zero(), |a, z| a + z)
            })
            .collect::<Vec<_>>();
        h[i] = es(&p, s.confidence);
    }
    let mut q = h[0] * h[0];
    for i in 1..5 {
        q += h[i] * h[i] * T::from_f64([1., 1., 2., 2., 6.][i]);
    }
    let liq = q.max_val(T::zero()).sqrt();
    FrtbImaReport {
        horizon_es: h,
        liquidity_adjusted_es: liq,
        total_ima_capital: liq + s.ses_addon,
        marginal_es: vec![0.; s.position_weights.len()],
        euler_component_es: vec![
            liq.val() / s.position_weights.len().max(1) as f64;
            s.position_weights.len()
        ],
    }
}
