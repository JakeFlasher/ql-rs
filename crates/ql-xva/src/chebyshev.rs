use crate::NettingSetSwap;
use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct AlbaneseCrepeyZeronBalanceSheetXvaSpec<T: Real = f64> {
    pub maturity: T,
    pub num_time_steps: usize,
    pub risk_free_rate: T,
    pub cpty_hazard_rate: T,
    pub cpty_recovery: T,
    pub bank_hazard_rate: T,
    pub bank_recovery: T,
    pub bank_liquidity_basis: T,
    pub im_funding_spread: T,
    pub shareholder_hurdle_rate: T,
    pub base_unmargined_epe: T,
    pub base_ene: T,
    pub simm_delta_weight: T,
    pub simm_vega_weight: T,
    pub simm_delta_vega_corr: T,
    pub regulatory_capital_ratio: T,
}

#[derive(Clone, Debug)]
pub struct AlbaneseCrepeyZeronBalanceSheetXvaReport<T: Real = f64> {
    pub chebyshev_expected_initial_margin_t0: T,
    pub chebyshev_peak_initial_margin: T,
    pub ftd_cva: T,
    pub ftd_dva: T,
    pub naive_fva: T,
    pub shareholder_fva_capital_deducted: T,
    pub mva_chebyshev_simm: T,
    pub kva_retained_earnings: T,
    pub double_counting_overlap_eliminated: T,
    pub shareholder_funds_transfer_price_ava: T,
    pub balance_sheet_wealth_transfer_to_creditors: T,
}

pub fn evaluate_albanese_crepey_zeron_balance_sheet_xva<T: Real>(
    s: &AlbaneseCrepeyZeronBalanceSheetXvaSpec<T>,
) -> AlbaneseCrepeyZeronBalanceSheetXvaReport<T> {
    let n = s.num_time_steps.max(1);
    let mut im0 = T::zero();
    let mut peak = T::zero();
    for i in 0..=n {
        let u = T::from_f64(i as f64 / n as f64);
        let basis = (T::one() - u) * s.base_unmargined_epe + u * s.base_ene.abs();
        let im = (s.simm_delta_weight * basis
            + s.simm_vega_weight * basis
            + T::from_f64(2.0)
                * s.simm_delta_vega_corr
                * s.simm_delta_weight
                * s.simm_vega_weight
                * basis
                * basis)
            .max_val(T::zero())
            .sqrt();
        if i == 0 {
            im0 = im;
        }
        peak = peak.max_val(im);
    }
    let df = (-s.risk_free_rate * s.maturity).exp();
    let cva =
        s.base_unmargined_epe * (T::one() - s.cpty_recovery) * s.cpty_hazard_rate * s.maturity * df;
    let dva =
        s.base_ene.abs() * (T::one() - s.bank_recovery) * s.bank_hazard_rate * s.maturity * df;
    let fva = s.base_unmargined_epe * s.im_funding_spread * s.maturity * df;
    let capital = s.regulatory_capital_ratio * (s.base_unmargined_epe + peak);
    let shareholder = fva - capital * s.shareholder_hurdle_rate * s.maturity;
    let mva = im0 * s.im_funding_spread * s.maturity;
    let kva = capital * s.shareholder_hurdle_rate * s.maturity;
    let overlap = (mva + kva).min_val(fva.abs()).max_val(T::zero());
    AlbaneseCrepeyZeronBalanceSheetXvaReport {
        chebyshev_expected_initial_margin_t0: im0,
        chebyshev_peak_initial_margin: peak,
        ftd_cva: cva,
        ftd_dva: dva,
        naive_fva: fva,
        shareholder_fva_capital_deducted: shareholder,
        mva_chebyshev_simm: mva,
        kva_retained_earnings: kva,
        double_counting_overlap_eliminated: overlap,
        shareholder_funds_transfer_price_ava: shareholder * T::from_f64(0.5),
        balance_sheet_wealth_transfer_to_creditors: cva + dva + overlap,
    }
}
use ql_curves::FlatForwardCurve;
use ql_curves::YieldTermStructure;
use ql_models::HullWhite1F;
#[derive(Clone, Debug)]
pub struct ChebyshevInterpolant1D {
    pub a: f64,
    pub b: f64,
    pub nodes: Vec<f64>,
    pub values: Vec<f64>,
    pub weights: Vec<f64>,
    pub coeffs: Vec<f64>,
}
impl ChebyshevInterpolant1D {
    #[allow(clippy::needless_range_loop)]
    pub fn from_fn(a: f64, b: f64, degree: usize, f: impl Fn(f64) -> f64) -> Self {
        let n = degree;
        let mut nodes = Vec::with_capacity(n + 1);
        let mut values = Vec::with_capacity(n + 1);
        if n == 0 {
            let x = 0.5 * (a + b);
            nodes.push(x);
            values.push(f(x));
            return Self {
                a,
                b,
                nodes,
                values: values.clone(),
                weights: vec![1.0],
                coeffs: values,
            };
        }
        for j in 0..=n {
            let x =
                (a + b) * 0.5 + (b - a) * 0.5 * (std::f64::consts::PI * j as f64 / n as f64).cos();
            nodes.push(x);
            values.push(f(x));
        }
        let weights = (0..=n)
            .map(|j| {
                let s = if j % 2 == 0 { 1. } else { -1. };
                s * if j == 0 || j == n { 0.5 } else { 1. }
            })
            .collect::<Vec<_>>();
        let mut coeffs = vec![0.; n + 1];
        for k in 0..=n {
            let ak = if k == 0 || k == n { 0.5 } else { 1. };
            coeffs[k] = 2. * ak / n as f64
                * (0..=n)
                    .map(|j| {
                        let aj = if j == 0 || j == n { 0.5 } else { 1. };
                        aj * values[j]
                            * (std::f64::consts::PI * j as f64 * k as f64 / n as f64).cos()
                    })
                    .sum::<f64>();
        }
        Self {
            a,
            b,
            nodes,
            values,
            weights,
            coeffs,
        }
    }
    pub fn eval_barycentric(&self, x: f64) -> f64 {
        for (j, &u) in self.nodes.iter().enumerate() {
            if (x - u).abs() < 1e-14 {
                return self.values[j];
            }
        }
        let mut num = 0.;
        let mut den = 0.;
        for j in 0..self.nodes.len() {
            let q = self.weights[j] / (x - self.nodes[j]);
            num += q * self.values[j];
            den += q;
        }
        num / den
    }
    pub fn eval_clenshaw<T: Real>(&self, x: T) -> T {
        let u = (x - T::from_f64((self.a + self.b) * 0.5)) / T::from_f64((self.b - self.a) * 0.5);
        let n = self.coeffs.len() - 1;
        let mut b1 = T::zero();
        let mut b2 = T::zero();
        for k in (1..=n).rev() {
            let b0 = T::from_f64(2.) * u * b1 - b2 + T::from_f64(self.coeffs[k]);
            b2 = b1;
            b1 = b0;
        }
        u * b1 - b2 + T::from_f64(self.coeffs[0])
    }
}
#[derive(Clone, Debug)]
pub struct ChebyshevTensor2D {
    pub x: ChebyshevInterpolant1D,
    pub y: ChebyshevInterpolant1D,
    pub values: Vec<Vec<f64>>,
    pub coeffs: Vec<Vec<f64>>,
}
impl ChebyshevTensor2D {
    #[allow(clippy::needless_range_loop)]
    pub fn from_fn(
        ax: f64,
        bx: f64,
        ay: f64,
        by: f64,
        degree: usize,
        f: impl Fn(f64, f64) -> f64,
    ) -> Self {
        let x = ChebyshevInterpolant1D::from_fn(ax, bx, degree, |_| 0.);
        let y = ChebyshevInterpolant1D::from_fn(ay, by, degree, |_| 0.);
        let values: Vec<Vec<f64>> = x
            .nodes
            .iter()
            .map(|&xx| y.nodes.iter().map(|&yy| f(xx, yy)).collect())
            .collect();
        let nx = x.nodes.len() - 1;
        let ny = y.nodes.len() - 1;
        let mut coeffs = vec![vec![0.; ny + 1]; nx + 1];
        if nx == 0 || ny == 0 {
            coeffs[0][0] = values[0][0];
        } else {
            for k in 0..=nx {
                for l in 0..=ny {
                    let ak = if k == 0 || k == nx { 0.5 } else { 1. };
                    let al = if l == 0 || l == ny { 0.5 } else { 1. };
                    coeffs[k][l] = 4. * ak * al / (nx as f64 * ny as f64)
                        * (0..=nx)
                            .map(|i| {
                                let ai = if i == 0 || i == nx { 0.5 } else { 1. };
                                (0..=ny)
                                    .map(|j| {
                                        let aj = if j == 0 || j == ny { 0.5 } else { 1. };
                                        ai * aj
                                            * values[i][j]
                                            * (std::f64::consts::PI * i as f64 * k as f64
                                                / nx as f64)
                                                .cos()
                                            * (std::f64::consts::PI * j as f64 * l as f64
                                                / ny as f64)
                                                .cos()
                                    })
                                    .sum::<f64>()
                            })
                            .sum::<f64>();
                }
            }
        }
        Self {
            x,
            y,
            values,
            coeffs,
        }
    }
    pub fn eval_barycentric(&self, x: f64, y: f64) -> f64 {
        let mut vals = Vec::with_capacity(self.x.nodes.len());
        for row in &self.values {
            let it = ChebyshevInterpolant1D {
                a: self.y.a,
                b: self.y.b,
                nodes: self.y.nodes.clone(),
                values: row.clone(),
                weights: self.y.weights.clone(),
                coeffs: self.y.coeffs.clone(),
            };
            vals.push(it.eval_barycentric(y));
        }
        let it = ChebyshevInterpolant1D {
            a: self.x.a,
            b: self.x.b,
            nodes: self.x.nodes.clone(),
            values: vals,
            weights: self.x.weights.clone(),
            coeffs: self.x.coeffs.clone(),
        };
        it.eval_barycentric(x)
    }
    pub fn eval_clenshaw<T: Real>(&self, x: T, y: T) -> T {
        let uy = (y - T::from_f64((self.y.a + self.y.b) * 0.5))
            / T::from_f64((self.y.b - self.y.a) * 0.5);
        let ux = (x - T::from_f64((self.x.a + self.x.b) * 0.5))
            / T::from_f64((self.x.b - self.x.a) * 0.5);
        let mut vals = Vec::with_capacity(self.coeffs.len());
        for row in &self.coeffs {
            let mut b1 = T::zero();
            let mut b2 = T::zero();
            for l in (1..row.len()).rev() {
                let b0 = T::from_f64(2.) * uy * b1 - b2 + T::from_f64(row[l]);
                b2 = b1;
                b1 = b0;
            }
            vals.push(uy * b1 - b2 + T::from_f64(row[0]));
        }
        let mut b1 = T::zero();
        let mut b2 = T::zero();
        for k in (1..vals.len()).rev() {
            let b0 = T::from_f64(2.) * ux * b1 - b2 + vals[k];
            b2 = b1;
            b1 = b0;
        }
        ux * b1 - b2 + vals[0]
    }
}
#[derive(Clone, Debug)]
pub struct DimMvaReport {
    pub times: Vec<f64>,
    pub eim_profile: Vec<f64>,
    pub pfe_im_95: Vec<f64>,
    pub mva: f64,
}
#[allow(clippy::needless_range_loop)]
fn simulate_hw1f_short_rate_paths(
    model: &HullWhite1F<f64>,
    curve: &FlatForwardCurve<f64>,
    times: &[f64],
    n_paths: usize,
    seed: u64,
) -> Vec<Vec<f64>> {
    let mut states = vec![vec![0.; n_paths]; times.len()];
    let rng = ql_mc::PhiloxRng::new(seed);
    for p in 0..n_paths {
        let mut x = 0.;
        let mut prev = 0.;
        let sign = if p % 2 == 0 { 1. } else { -1. };
        let id = (p / 2) as u64;
        for (k, &tm) in times.iter().enumerate() {
            let dt = tm - prev;
            let z = rng.normal4(id, k as u32)[0] * sign;
            let e = (-model.a * dt).exp();
            let var = if model.a.abs() < 1e-10 {
                dt
            } else {
                (1. - (-2. * model.a * dt).exp()) / (2. * model.a)
            };
            x = e * x + model.sigma * var.sqrt() * z;
            let b = if model.a.abs() < 1e-12 {
                tm
            } else {
                (1. - (-model.a * tm).exp()) / model.a
            };
            states[k][p] =
                curve.instantaneous_forward(tm) + 0.5 * model.sigma * model.sigma * b * b + x;
            prev = tm;
        }
    }
    states
}
fn swap_rate_sensitivity_dual(
    model: &HullWhite1F<f64>,
    curve: &FlatForwardCurve<f64>,
    swap: &NettingSetSwap,
    t: f64,
    r: f64,
) -> f64 {
    let md = HullWhite1F {
        a: ql_ad::Dual::<1>::from_value(model.a),
        sigma: ql_ad::Dual::<1>::from_value(model.sigma),
    };
    let cd = FlatForwardCurve::new(
        curve.reference_date,
        ql_ad::Dual::<1>::from_value(curve.rate),
    );
    let rd = ql_ad::Dual::<1>::variable(r, 0);
    swap.npv_at(&md, &cd, t, rd).d[0] * 1e-4
}

fn swap_rate_sensitivity_fast(
    model: &HullWhite1F<f64>,
    curve: &FlatForwardCurve<f64>,
    swap: &NettingSetSwap,
    t: f64,
    r: f64,
) -> f64 {
    // Differentiate the affine Hull--White bond representation analytically.
    let mut fixed_dv = 0.0;
    let mut floating_dv = 0.0;
    let mut last = None;
    for (i, &tm) in swap.payment_times.iter().enumerate() {
        if tm > t {
            let p = model.zcb_price(curve, t, tm, r);
            let b = model.b_factor(t, tm);
            let accrual = swap.accruals.get(i).copied().unwrap_or(1.0);
            fixed_dv += swap.fixed_rate * accrual * b * p;
            last = Some((b, p));
        }
    }
    if let Some((b, p)) = last {
        floating_dv = b * p;
    }
    let sign = if swap.is_payer { 1.0 } else { -1.0 };
    // A payer swap has NPV = (1 - P_last) - fixed_leg.  Since
    // dP/dr = -B P, both -P_last and -fixed_leg contribute positively to
    // the short-rate derivative; the receiver simply reverses the sign.
    sign * swap.notional * (floating_dv + fixed_dv) * 1e-4
}
fn simm_from_swap_sensitivities(sens: &[f64]) -> f64 {
    simm_from_iter(sens.iter().copied())
}

/// Closed form of the one-bucket aggregation used by the DIM routines.
///
/// The interpolation profiles supplied here represent sensitivities to one
/// common pillar/factor.  Net those sensitivities *before* applying the SIMM
/// concentration multiplier and squaring; treating each trade as an
/// independent factor would incorrectly remove portfolio netting.  Keeping
/// this explicit one-factor aggregation also makes the convention clear to
/// callers that build profiles from several trades.
fn simm_from_iter<I>(sensitivities: I) -> f64
where
    I: IntoIterator<Item = f64>,
{
    let net: f64 = sensitivities.into_iter().sum();
    let concentration = (net.abs() / 1e6).sqrt().max(1.0);
    (net * 50.0 * concentration).abs()
}
fn report_from_sens(
    times: &[f64],
    curve: &FlatForwardCurve<f64>,
    profiles: &[Vec<f64>],
    spread: f64,
) -> DimMvaReport {
    let mut eim = Vec::new();
    let mut pfe = Vec::new();
    for v in profiles {
        if v.is_empty() {
            eim.push(0.0);
            pfe.push(0.0);
            continue;
        }
        let mut a = v.clone();
        a.sort_by(f64::total_cmp);
        eim.push(a.iter().sum::<f64>() / a.len().max(1) as f64);
        let i = (0.95 * a.len().saturating_sub(1) as f64).round() as usize;
        pfe.push(a[i.min(a.len().saturating_sub(1))]);
    }
    let mut mva = 0.;
    for k in 0..times.len() {
        let prev = if k == 0 { 0. } else { times[k - 1] };
        mva += spread * curve.discount_factor(times[k]) * eim[k] * (times[k] - prev);
    }
    DimMvaReport {
        times: times.to_vec(),
        eim_profile: eim,
        pfe_im_95: pfe,
        mva,
    }
}
#[allow(clippy::too_many_arguments)]
pub fn compute_bruteforce_dim_and_mva(
    model: &HullWhite1F<f64>,
    curve: &FlatForwardCurve<f64>,
    swaps: &[NettingSetSwap],
    times: &[f64],
    n_paths: usize,
    _degree: usize,
    spread: f64,
    seed: u64,
) -> DimMvaReport {
    let states = simulate_hw1f_short_rate_paths(model, curve, times, n_paths, seed);
    let mut prof = vec![vec![0.; n_paths]; times.len()];
    for k in 0..times.len() {
        for p in 0..n_paths {
            let sens = swaps
                .iter()
                .map(|s| swap_rate_sensitivity_dual(model, curve, s, times[k], states[k][p]))
                .collect::<Vec<_>>();
            prof[k][p] = simm_from_swap_sensitivities(&sens);
        }
    }
    report_from_sens(times, curve, &prof, spread)
}
#[allow(clippy::too_many_arguments)]
pub fn compute_chebyshev_dim_and_mva(
    model: &HullWhite1F<f64>,
    curve: &FlatForwardCurve<f64>,
    swaps: &[NettingSetSwap],
    times: &[f64],
    n_paths: usize,
    degree: usize,
    spread: f64,
    seed: u64,
) -> DimMvaReport {
    let states = simulate_hw1f_short_rate_paths(model, curve, times, n_paths, seed);
    let mut prof = vec![vec![0.; n_paths]; times.len()];
    for k in 0..times.len() {
        let rmin = states[k].iter().copied().fold(f64::INFINITY, f64::min);
        let mut rmax = states[k].iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mut rmin = rmin;
        if (rmax - rmin).abs() < 1e-12 {
            rmin -= 1e-6;
            rmax += 1e-6;
        }
        let interps = swaps
            .iter()
            .map(|s| {
                ChebyshevInterpolant1D::from_fn(rmin, rmax, degree, |r| {
                    swap_rate_sensitivity_fast(model, curve, s, times[k], r)
                })
            })
            .collect::<Vec<_>>();
        for p in 0..n_paths {
            prof[k][p] = simm_from_iter(interps.iter().map(|it| it.eval_barycentric(states[k][p])));
        }
    }
    report_from_sens(times, curve, &prof, spread)
}

#[cfg(test)]
mod phase28_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(x: T) -> AlbaneseCrepeyZeronBalanceSheetXvaSpec<T> {
        AlbaneseCrepeyZeronBalanceSheetXvaSpec {
            maturity: T::from_f64(2.0),
            num_time_steps: 8,
            risk_free_rate: T::from_f64(0.03),
            cpty_hazard_rate: T::from_f64(0.02),
            cpty_recovery: T::from_f64(0.4),
            bank_hazard_rate: T::from_f64(0.01),
            bank_recovery: T::from_f64(0.4),
            bank_liquidity_basis: T::from_f64(0.005),
            im_funding_spread: T::from_f64(0.01),
            shareholder_hurdle_rate: T::from_f64(0.1),
            base_unmargined_epe: T::from_f64(10.0) + x,
            base_ene: T::from_f64(4.0),
            simm_delta_weight: T::from_f64(0.6),
            simm_vega_weight: T::from_f64(0.4),
            simm_delta_vega_corr: T::from_f64(0.25),
            regulatory_capital_ratio: T::from_f64(0.08),
        }
    }
    #[test]
    fn xva_invariants_and_dual() {
        let r = evaluate_albanese_crepey_zeron_balance_sheet_xva(&mk(0.0));
        assert!(r.chebyshev_peak_initial_margin >= r.chebyshev_expected_initial_margin_t0);
        assert!(r.shareholder_fva_capital_deducted < r.naive_fva);
        assert!(r.double_counting_overlap_eliminated > 0.0);
        let h = 1e-5;
        let fd = (evaluate_albanese_crepey_zeron_balance_sheet_xva(&mk(h)).ftd_cva
            - evaluate_albanese_crepey_zeron_balance_sheet_xva(&mk(-h)).ftd_cva)
            / (2.0 * h);
        let ad = evaluate_albanese_crepey_zeron_balance_sheet_xva(&mk(Dual::<1>::variable(0.0, 0)))
            .ftd_cva;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
