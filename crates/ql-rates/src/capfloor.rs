use ql_ad::{Dual, Real, TapeScope, TapeVar};
use ql_curves::{DiscountCurve, YieldTermStructure};
use ql_math::brent_solve;

#[derive(Clone, Debug)]
pub struct CapFloorLet {
    pub notional: f64,
    pub strike: f64,
    pub fixing_time: f64,
    pub start_time: f64,
    pub end_time: f64,
    pub pay_time: f64,
    pub accrual: f64,
    pub is_cap: bool,
}
#[derive(Clone, Debug)]
pub struct IborCapFloor {
    pub caplets: Vec<CapFloorLet>,
}
impl IborCapFloor {
    pub fn from_schedule(notional: f64, strike: f64, period_times: &[f64], is_cap: bool) -> Self {
        let mut caplets = Vec::new();
        for (j, w) in period_times.windows(2).enumerate() {
            let tau = (w[1] - w[0]).max(1e-12);
            caplets.push(CapFloorLet {
                notional,
                strike,
                fixing_time: if j == 0 && w[0] > 0.0 {
                    w[0].max(1e-4)
                } else {
                    w[0]
                },
                start_time: w[0],
                end_time: w[1],
                pay_time: w[1],
                accrual: tau,
                is_cap,
            });
        }
        Self { caplets }
    }
    pub fn price_bachelier<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        ois: &D,
        fwd: &F,
        vols: &[T],
    ) -> T {
        self.caplets
            .iter()
            .enumerate()
            .fold(T::zero(), |sum, (j, c)| {
                let tau = c.accrual;
                let f = (fwd.discount_factor(c.start_time) / fwd.discount_factor(c.end_time)
                    - T::one())
                    / tau;
                let vol = vols.get(j).copied().unwrap_or_else(T::zero);
                let st = vol * T::from_f64(c.fixing_time.sqrt());
                let omega = if c.is_cap { 1.0 } else { -1.0 };
                if st.val().abs() <= 1e-14 {
                    return sum
                        + T::from_f64(c.notional * c.accrual)
                            * ois.discount_factor(c.pay_time)
                            * (T::from_f64(omega) * (f - c.strike)).max_val(T::zero());
                }
                let d = if st.val().abs() > 1e-14 {
                    (f - c.strike) / st
                } else {
                    T::from_f64(if (f.val() - c.strike) * omega > 0.0 {
                        1e3
                    } else {
                        -1e3
                    })
                };
                let unit =
                    T::from_f64(omega) * (f - c.strike) * (T::from_f64(omega) * d).norm_cdf()
                        + st * d.norm_pdf();
                sum + T::from_f64(c.notional * c.accrual) * ois.discount_factor(c.pay_time) * unit
            })
    }
    pub fn price_shifted_black76<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        ois: &D,
        fwd: &F,
        vols: &[T],
        shift: f64,
    ) -> T {
        self.caplets
            .iter()
            .enumerate()
            .fold(T::zero(), |sum, (j, c)| {
                let tau = c.accrual;
                let f = (fwd.discount_factor(c.start_time) / fwd.discount_factor(c.end_time)
                    - T::one())
                    / tau;
                let v = vols.get(j).copied().unwrap_or_else(T::zero);
                let ft = f + shift;
                let kt = T::from_f64(c.strike + shift);
                let st = v * T::from_f64(c.fixing_time.sqrt());
                let omega = if c.is_cap { 1.0 } else { -1.0 };
                if st.val().abs() <= 1e-14 || ft.val() <= 0.0 || kt.val() <= 0.0 {
                    return sum
                        + T::from_f64(c.notional * c.accrual)
                            * ois.discount_factor(c.pay_time)
                            * (T::from_f64(omega) * (f - c.strike)).max_val(T::zero());
                }
                let d1 = ((ft / kt).ln() + T::from_f64(0.5) * st * st) / st;
                let d2 = d1 - st;
                let unit = T::from_f64(omega)
                    * (ft * (T::from_f64(omega) * d1).norm_cdf()
                        - kt * (T::from_f64(omega) * d2).norm_cdf());
                sum + T::from_f64(c.notional * c.accrual) * ois.discount_factor(c.pay_time) * unit
            })
    }
}
#[derive(Clone, Debug)]
pub struct CalibratedCapletVolCurve {
    pub caplet_vols: Vec<f64>,
    pub ift_jacobian: Vec<Vec<f64>>,
    pub quotes: Vec<(IborCapFloor, f64)>,
}
impl CalibratedCapletVolCurve {
    #[allow(clippy::needless_range_loop)]
    pub fn bootstrap_bachelier(
        ois: &DiscountCurve<f64>,
        fwd: &DiscountCurve<f64>,
        quotes: &[(IborCapFloor, f64)],
    ) -> Result<Self, String> {
        if quotes.is_empty() {
            return Err("empty quote strip".into());
        }
        let m = quotes.len();
        let mut vols = vec![0.01; m];
        for i in 0..m {
            let prior = vols.clone();
            let (cap, flat) = (&quotes[i].0, quotes[i].1);
            let f = |x: f64| {
                let mut v = prior.clone();
                v[i] = x;
                let pv = cap.price_bachelier(ois, fwd, &v);
                let fv = cap.price_bachelier(ois, fwd, &vec![flat; cap.caplets.len()]);
                pv - fv
            };
            vols[i] = brent_solve(f, 1e-6, 0.10, 1e-12, 200).map_err(|e| format!("{e}"))?;
        }
        let mut j = vec![vec![0.; m]; m];
        let mut s = vec![0.; m];
        for col in 0..m {
            let advols: Vec<Dual<1>> = vols
                .iter()
                .enumerate()
                .map(|(i, &x)| {
                    if i == col {
                        Dual::variable(x, 0)
                    } else {
                        Dual::from_value(x)
                    }
                })
                .collect();
            for row in 0..m {
                let adois = ois.to_ad_curve(
                    &ois.discounts
                        .iter()
                        .map(|&x| Dual::from_value(x))
                        .collect::<Vec<_>>(),
                );
                let adf = fwd.to_ad_curve(
                    &fwd.discounts
                        .iter()
                        .map(|&x| Dual::from_value(x))
                        .collect::<Vec<_>>(),
                );
                j[row][col] = quotes[row].0.price_bachelier(&adois, &adf, &advols).d[0];
            }
        }
        for row in 0..m {
            let flat: Dual<1> = Dual::variable(quotes[row].1, 0);
            let adois = ois.to_ad_curve(
                &ois.discounts
                    .iter()
                    .map(|&x| Dual::from_value(x))
                    .collect::<Vec<_>>(),
            );
            let adf = fwd.to_ad_curve(
                &fwd.discounts
                    .iter()
                    .map(|&x| Dual::from_value(x))
                    .collect::<Vec<_>>(),
            );
            s[row] = quotes[row]
                .0
                .price_bachelier(&adois, &adf, &vec![flat; quotes[row].0.caplets.len()])
                .d[0];
        }
        let mut ift = vec![vec![0.; m]; m];
        for col in 0..m {
            for row in col..m {
                let diag = j[row][row];
                if diag.abs() > 1e-14 {
                    let mut rhs = if row == col { s[row] } else { 0.0 };
                    for k in col..row {
                        rhs -= j[row][k] * ift[k][col];
                    }
                    // F_i = V_i(sigma) - V_i^flat(Sigma_i), hence
                    // ∂F/∂Sigma_i = -s_i.  The IFT map is
                    // -J^{-1}F_Sigma = +J^{-1}s, so the triangular solve
                    // uses the positive right-hand side here.
                    ift[row][col] = rhs / diag;
                }
            }
        }
        Ok(Self {
            caplet_vols: vols,
            ift_jacobian: ift,
            quotes: quotes.to_vec(),
        })
    }
    pub fn value_and_par_cap_vega<F: FnOnce(&[TapeVar<f64>]) -> TapeVar<f64>>(
        &self,
        pricer: F,
    ) -> (f64, Vec<f64>) {
        let guard = TapeScope::guard_with_capacity(self.caplet_vols.len() * 64);
        let vars: Vec<_> = self
            .caplet_vols
            .iter()
            .map(|&x| TapeScope::with_active(|t| t.input(x)))
            .collect();
        let out = pricer(&vars);
        let adj = TapeScope::gradient_from(out);
        let g: Vec<f64> = vars.iter().map(|v| adj.gradient(*v)).collect();
        let _ = guard.stop();
        let mut q = vec![0.; self.ift_jacobian.len()];
        for j in 0..g.len() {
            for k in 0..q.len() {
                q[k] += g[j] * self.ift_jacobian[j][k];
            }
        }
        (out.val, q)
    }
}
pub type CapletVolStripper = CalibratedCapletVolCurve;

#[derive(Clone, Debug)]
pub struct RangeAccrualSpec<T: Real = f64> {
    pub notional: T,
    pub coupon_rate: T,
    pub accrual_fraction: f64,
    pub pay_time: f64,
    pub fixing_times: Vec<f64>,
    pub rate_tenor: f64,
    pub lower_bound: T,
    pub upper_bound: T,
    pub atm_normal_vol: T,
    pub smile_skew: T,
    pub timing_correlation: T,
}

pub fn price_range_accrual_coupon<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
    spec: &RangeAccrualSpec<T>,
    ois: &D,
    fwd: &F,
) -> T {
    if spec.fixing_times.is_empty() {
        return T::zero();
    }
    let delta = spec.rate_tenor.max(1e-4);
    let mut prob = T::zero();
    for &tm in &spec.fixing_times {
        let f = (fwd.discount_factor(tm) / fwd.discount_factor(tm + delta) - 1.0) / delta;
        let adj = f + T::from_f64(delta) / (T::one() + f * delta).max_val(T::from_f64(1e-8))
            * spec.timing_correlation
            * spec.atm_normal_vol
            * spec.atm_normal_vol
            * tm
            * ((tm + delta - spec.pay_time) / delta);
        let digital = |k: T| {
            let sig =
                (spec.atm_normal_vol + spec.smile_skew * (k - adj)).max_val(T::from_f64(1e-6));
            let st = sig * tm.max(1e-8).sqrt();
            let d = (adj - k) / st;
            (d.norm_cdf() - d.norm_pdf() * tm.max(1e-8).sqrt() * spec.smile_skew)
                .max_val(T::zero())
                .min_val(T::one())
        };
        prob += (digital(spec.lower_bound) - digital(spec.upper_bound)).max_val(T::zero());
    }
    spec.notional
        * ois.discount_factor(spec.pay_time)
        * spec.accrual_fraction
        * spec.coupon_rate
        * prob
        / spec.fixing_times.len() as f64
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum PathDependentCapFloorKind {
    Ratchet,
    StickyCap,
    StickyFloor,
}

#[derive(Clone, Debug)]
pub struct StickyRatchetCapFloorSpec<T: Real = f64> {
    pub notional: T,
    pub period_times: Vec<f64>,
    pub initial_coupon: T,
    pub spread: T,
    pub max_step_up: T,
    pub max_step_down: T,
    pub smoothing_eps: f64,
    pub kind: PathDependentCapFloorKind,
}

pub fn price_sticky_ratchet_capfloor_hw1f<
    T: Real,
    D: YieldTermStructure<T>,
    F: YieldTermStructure<T>,
>(
    spec: &StickyRatchetCapFloorSpec<T>,
    model: &ql_models::HullWhite1F<T>,
    ois: &D,
    fwd: &F,
    n_paths: usize,
    seed: u64,
) -> (T, Vec<T>) {
    let eps = spec.smoothing_eps.max(1e-6);
    let smooth_pos = |u: T| ((u + (u * u + eps * eps).sqrt() - eps).max_val(T::zero())) * 0.5;
    let smooth_min = |a: T, b: T| a - smooth_pos(a - b);
    let smooth_max = |a: T, b: T| a + smooth_pos(b - a);
    let m = spec.period_times.len().saturating_sub(1);
    if m == 0 {
        return (T::zero(), Vec::new());
    }
    let rng = ql_mc::PhiloxRng::new(seed);
    let mut total = T::zero();
    let mut path_vals = Vec::with_capacity(n_paths);
    for pth in 0..n_paths.max(1) {
        let anti = pth % 2 == 1;
        let id = (pth / 2) as u64;
        let mut x = T::zero();
        let mut pv = T::zero();
        let mut coupon = spec.initial_coupon;
        let mut prev_t = 0.0;
        let mut log_df = T::zero();
        for k in 1..=m {
            let t0 = spec.period_times[k - 1];
            let t1 = spec.period_times[k];
            let dt = (t1 - prev_t).max(0.0);
            let z = rng.normal4(id, (k - 1) as u32)[0] * if anti { -1.0 } else { 1.0 };
            let var = if model.a.val().abs() < 1e-8 {
                model.sigma * model.sigma * dt
            } else {
                model.sigma * model.sigma * (T::one() - (-model.a * dt * 2.0).exp())
                    / (model.a * 2.0)
            };
            x = x * (-model.a * dt).exp() + var.max_val(T::zero()).sqrt() * z;
            let p0 = ois.discount_factor(t0);
            let p1 = model.zcb_price(ois, t0, t1, x);
            let ratio = fwd.discount_factor(t0)
                / fwd.discount_factor(t1)
                / (ois.discount_factor(t0) / ois.discount_factor(t1));
            let tau = (t1 - t0).max(1e-12);
            let l = (ratio / p1 - 1.0) / tau;
            let cf = match spec.kind {
                PathDependentCapFloorKind::Ratchet => {
                    coupon = smooth_min(
                        smooth_max(l + spec.spread, coupon - spec.max_step_down),
                        coupon + spec.max_step_up,
                    );
                    spec.notional * tau * coupon
                }
                PathDependentCapFloorKind::StickyCap => {
                    let kcap = coupon + spec.spread;
                    let v = spec.notional * tau * smooth_pos(l - kcap);
                    coupon = smooth_min(l, kcap);
                    v
                }
                PathDependentCapFloorKind::StickyFloor => {
                    let kfloor = coupon - spec.spread;
                    let v = spec.notional * tau * smooth_pos(kfloor - l);
                    coupon = smooth_max(l, kfloor);
                    v
                }
            };
            // Normalize the simulated discount factor to the market bond ratio.
            log_df += -((ois.discount_factor(t1) / ois.discount_factor(prev_t.max(0.0))).ln());
            let _ = p0;
            pv += cf * log_df.exp();
            prev_t = t1;
        }
        total += pv;
        path_vals.push(pv);
    }
    (total / n_paths.max(1) as f64, path_vals)
}

#[cfg(test)]
mod phase15_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;
    #[test]
    fn range_accrual_and_sticky_ad() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let c = FlatForwardCurve::new(d, 0.03);
        let mk = |lo, hi| RangeAccrualSpec {
            notional: 1.0,
            coupon_rate: 0.05,
            accrual_fraction: 0.5,
            pay_time: 1.0,
            fixing_times: vec![0.25, 0.5, 0.75],
            rate_tenor: 0.25,
            lower_bound: lo,
            upper_bound: hi,
            atm_normal_vol: 0.01,
            smile_skew: 0.1,
            timing_correlation: 0.2,
        };
        let a = price_range_accrual_coupon(&mk(0.0, 0.03), &c, &c);
        let b = price_range_accrual_coupon(&mk(0.03, 0.06), &c, &c);
        let ab = price_range_accrual_coupon(&mk(0.0, 0.06), &c, &c);
        assert!((a + b - ab).abs() < 1e-12);
        let s = StickyRatchetCapFloorSpec {
            notional: 1.0,
            period_times: vec![0.0, 0.5, 1.0],
            initial_coupon: 0.03,
            spread: 0.0,
            max_step_up: 0.01,
            max_step_down: 0.01,
            smoothing_eps: 1e-4,
            kind: PathDependentCapFloorKind::Ratchet,
        };
        let m = ql_models::HullWhite1F {
            a: 0.1,
            sigma: 0.01,
        };
        let v = price_sticky_ratchet_capfloor_hw1f(&s, &m, &c, &c, 64, 7).0;
        let sd = ql_models::HullWhite1F::<Dual<1>> {
            a: Dual::from_value(0.1),
            sigma: Dual::variable(0.01, 0),
        };
        let cd = FlatForwardCurve::new(d, Dual::<1>::from_value(0.03));
        let ad = price_sticky_ratchet_capfloor_hw1f(
            &StickyRatchetCapFloorSpec {
                notional: Dual::from_value(1.0),
                period_times: s.period_times.clone(),
                initial_coupon: Dual::from_value(0.03),
                spread: Dual::from_value(0.0),
                max_step_up: Dual::variable(0.01, 0),
                max_step_down: Dual::from_value(0.01),
                smoothing_eps: 1e-4,
                kind: s.kind,
            },
            &sd,
            &cd,
            &cd,
            64,
            7,
        )
        .0;
        assert!(v.is_finite() && ad.d[0].is_finite());
    }
}
