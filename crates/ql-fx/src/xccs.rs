//! Cross-currency basis, FX forwards, and CSA-collateralized XCCS curves.

use crate::FxError;
use nalgebra::{DMatrix, DVector};
use ql_ad::{Dual, Real, TapeScope, TapeVar};
use ql_curves::{DiscountCurve, YieldTermStructure};
use ql_math::brent_solve;
use ql_time::Date;

pub fn fx_forward_rate<T: Real>(spot_fx: T, dom_csa_df: T, for_csa_df: T) -> T {
    spot_fx * for_csa_df / dom_csa_df
}

pub fn fx_swap_points<T: Real>(spot_fx: T, dom_csa_df: T, for_csa_df: T) -> T {
    fx_forward_rate(spot_fx, dom_csa_df, for_csa_df) - spot_fx
}

pub fn implied_foreign_csa_df_from_fx_forward<T: Real>(
    spot_fx: T,
    forward_fx: T,
    dom_csa_df: T,
) -> T {
    dom_csa_df * forward_fx / spot_fx
}

#[derive(Clone, Debug)]
pub struct XccsQuote {
    pub maturity_years: f64,
    pub basis_spread: f64,
}

#[derive(Clone, Debug)]
pub struct ConstNotionalXccs {
    pub dom_notional: f64,
    pub for_notional: f64,
    pub spot_fx: f64,
    pub payment_times: Vec<f64>,
    pub accrual_fractions: Vec<f64>,
    pub for_basis_spread: f64,
    pub dom_spread: f64,
    pub receive_foreign: bool,
}

impl ConstNotionalXccs {
    pub fn npv<T: Real>(
        &self,
        dom_ois: &impl YieldTermStructure<T>,
        for_ois: &impl YieldTermStructure<T>,
        for_csa_discount: &impl YieldTermStructure<T>,
    ) -> T {
        let mut dom = T::from_f64(-self.dom_notional);
        let mut foreign = T::from_f64(-self.spot_fx * self.for_notional);
        let mut prev = 0.0;
        for (j, &t) in self.payment_times.iter().enumerate() {
            let dt = self.accrual_fractions.get(j).copied().unwrap_or(t - prev);
            let df_d0 = dom_ois.discount_factor(prev);
            let df_d1 = dom_ois.discount_factor(t);
            let df_f0 = for_ois.discount_factor(prev);
            let df_f1 = for_ois.discount_factor(t);
            let f_d = (df_d0 / df_d1 - T::one()) / dt;
            let f_f = (df_f0 / df_f1 - T::one()) / dt;
            dom += T::from_f64(self.dom_notional)
                * (f_d + self.dom_spread)
                * dt
                * dom_ois.discount_factor(t);
            foreign += T::from_f64(self.for_notional * self.spot_fx)
                * (f_f + self.for_basis_spread)
                * dt
                * for_csa_discount.discount_factor(t);
            prev = t;
        }
        if let Some(&last) = self.payment_times.last() {
            dom += T::from_f64(self.dom_notional) * dom_ois.discount_factor(last);
            foreign += T::from_f64(self.for_notional * self.spot_fx)
                * for_csa_discount.discount_factor(last);
        }
        if self.receive_foreign {
            foreign - dom
        } else {
            dom - foreign
        }
    }

    pub fn fair_foreign_basis_spread<T: Real>(
        &self,
        dom_ois: &impl YieldTermStructure<T>,
        for_ois: &impl YieldTermStructure<T>,
        for_csa_discount: &impl YieldTermStructure<T>,
    ) -> T {
        let mut zero = self.clone();
        zero.for_basis_spread = 0.0;
        let mut one = zero.clone();
        one.for_basis_spread = 1.0;
        let p0 = zero.npv(dom_ois, for_ois, for_csa_discount);
        let p1 = one.npv(dom_ois, for_ois, for_csa_discount);
        -p0 / (p1 - p0)
    }
}

#[derive(Clone, Debug)]
pub struct CalibratedXccsCurve {
    pub ref_date: Date,
    pub spot_fx: f64,
    pub dom_ois: DiscountCurve<f64>,
    pub for_ois: DiscountCurve<f64>,
    pub for_csa_curve: DiscountCurve<f64>,
    pub quotes: Vec<XccsQuote>,
    pub basis_jacobian: DMatrix<f64>,
}

impl CalibratedXccsCurve {
    pub fn bootstrap(
        ref_date: Date,
        spot_fx: f64,
        dom_ois: DiscountCurve<f64>,
        for_ois: DiscountCurve<f64>,
        quotes: &[XccsQuote],
    ) -> Result<Self, FxError> {
        if quotes.is_empty()
            || quotes
                .windows(2)
                .any(|w| w[1].maturity_years <= w[0].maturity_years)
        {
            return Err(FxError::InvalidQuote);
        }
        let times: Vec<f64> = quotes.iter().map(|q| q.maturity_years).collect();
        let mut discounts = Vec::with_capacity(times.len());
        for (k, q) in quotes.iter().enumerate() {
            let prior = discounts.clone();
            let times_k = times.clone();
            let f = |x: f64| {
                let mut ds = prior.clone();
                ds.push(x);
                let tail = ds.len();
                let curve = DiscountCurve::new(ref_date, times_k[..tail].to_vec(), ds);
                let inst = ConstNotionalXccs {
                    dom_notional: spot_fx,
                    for_notional: 1.0,
                    spot_fx,
                    payment_times: (1..=q.maturity_years.round() as usize)
                        .map(|i| i as f64)
                        .collect(),
                    accrual_fractions: vec![1.0; q.maturity_years.round() as usize],
                    for_basis_spread: q.basis_spread,
                    dom_spread: 0.0,
                    receive_foreign: true,
                };
                inst.npv(&dom_ois, &for_ois, &curve)
            };
            let guess = for_ois.discount_factor(q.maturity_years);
            let root = brent_solve(
                f,
                (guess * 0.25).max(1e-8),
                (guess * 1.75).min(1.25),
                1e-14,
                200,
            )
            .map_err(|_| FxError::InvalidQuote)?;
            discounts.push(root);
            let _ = k;
        }
        let curve = DiscountCurve::new(ref_date, times.clone(), discounts.clone());
        let m = quotes.len();
        let mut jt = DMatrix::zeros(m, m);
        let mut jq = DMatrix::zeros(m, m);
        for i in 0..m {
            for j in 0..m {
                let nodes: Vec<Dual<1>> = discounts
                    .iter()
                    .enumerate()
                    .map(|(k, &x)| {
                        if k == j {
                            Dual::variable(x, 0)
                        } else {
                            Dual::from_value(x)
                        }
                    })
                    .collect();
                let dc = DiscountCurve::new(ref_date, times.clone(), nodes);
                let q = &quotes[i];
                let inst = ConstNotionalXccs {
                    dom_notional: spot_fx,
                    for_notional: 1.0,
                    spot_fx,
                    payment_times: (1..=q.maturity_years.round() as usize)
                        .map(|x| x as f64)
                        .collect(),
                    accrual_fractions: vec![1.0; q.maturity_years.round() as usize],
                    for_basis_spread: q.basis_spread,
                    dom_spread: 0.0,
                    receive_foreign: true,
                };
                jt[(i, j)] = inst
                    .npv(
                        &dom_ois.to_ad_curve(
                            &dom_ois
                                .discounts
                                .iter()
                                .map(|&x| Dual::from_value(x))
                                .collect::<Vec<_>>(),
                        ),
                        &for_ois.to_ad_curve(
                            &for_ois
                                .discounts
                                .iter()
                                .map(|&x| Dual::from_value(x))
                                .collect::<Vec<_>>(),
                        ),
                        &dc,
                    )
                    .d[0];
            }
            let q = &quotes[i];
            // The quote derivative is explicit and linear in the foreign
            // basis coupons (the receive-foreign orientation is positive).
            let n = q.maturity_years.round() as usize;
            jq[(i, i)] = spot_fx
                * (1..=n)
                    .map(|k| curve.discount_factor(k as f64))
                    .sum::<f64>();
        }
        let basis_jacobian = jt.lu().solve(&(-jq)).ok_or(FxError::InvalidQuote)?;
        Ok(Self {
            ref_date,
            spot_fx,
            dom_ois,
            for_ois,
            for_csa_curve: curve,
            quotes: quotes.to_vec(),
            basis_jacobian,
        })
    }

    pub fn value_and_par_basis_dv01<F>(&self, pricer: F) -> (f64, Vec<f64>)
    where
        F: Fn(&DiscountCurve<TapeVar<f64>>) -> TapeVar<f64>,
    {
        TapeScope::start_recording_with_capacity(8192);
        let vars: Vec<_> = self
            .for_csa_curve
            .discounts
            .iter()
            .map(|&x| TapeScope::with_active(|t| t.input(x)))
            .collect();
        let curve = DiscountCurve::new(
            self.ref_date,
            self.for_csa_curve.times.clone(),
            vars.clone(),
        );
        let out = pricer(&curve);
        let adj = TapeScope::gradient_from(out);
        let grad = DVector::from_vec(vars.iter().map(|v| adj.gradient(*v)).collect());
        TapeScope::stop_recording();
        let q = self.basis_jacobian.transpose() * grad;
        (out.val, q.iter().copied().collect())
    }
}
