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
    pub fn collateralized_exposure<T: Real>(&self, x: T) -> T {
        if x.val() > self.threshold_cpty + self.mta {
            x - T::from_f64(self.threshold_cpty)
        } else if x.val() < -(self.threshold_self + self.mta) {
            x + T::from_f64(self.threshold_self)
        } else {
            x
        }
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
                    value += cube.get(ns, k, path).val();
                }
                sum += value.max(0.0);
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
                    value += cube.get(ns, k, path).val();
                }
                values.push(value.max(0.0));
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
        let mut vals = Vec::with_capacity(cube.n_paths);
        for p in 0..cube.n_paths {
            let mut x = T::zero();
            for ns in 0..cube.n_netting_sets {
                x += cube.get(ns, k, p);
            }
            vals.push(csa.map(|c| c.collateralized_exposure(x)).unwrap_or(x));
        }
        let e = vals
            .iter()
            .map(|x| x.max_val(T::zero()))
            .fold(T::zero(), |a, b| a + b)
            / T::from_f64(cube.n_paths.max(1) as f64);
        let n = vals
            .iter()
            .map(|x| x.min_val(T::zero()))
            .fold(T::zero(), |a, b| a + b)
            / T::from_f64(cube.n_paths.max(1) as f64);
        let mut raw = vals.iter().map(|x| x.val().max(0.)).collect::<Vec<_>>();
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
    let epe =
        ee.iter().copied().fold(T::zero(), |a, b| a + b) / T::from_f64(ee.len().max(1) as f64);
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
