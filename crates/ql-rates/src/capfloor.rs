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
