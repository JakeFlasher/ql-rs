//! Cross-currency basis, FX forwards, and CSA-collateralized XCCS curves.

use crate::FxError;
use nalgebra::{DMatrix, DVector};
use ql_ad::{Dual, Real, TapeScope, TapeVar};
use ql_curves::{DiscountCurve, YieldTermStructure};
use ql_math::brent_solve;
use ql_time::Date;

#[derive(Clone, Debug)]
pub struct BaranBianchettiMtmXccyBasisSpec<T: Real = f64> {
    pub spot_fx_dom_per_for: T,
    pub foreign_notional: T,
    pub payment_times: Vec<T>,
    pub dom_ois_discount_factors: Vec<T>,
    pub for_xccy_discount_factors: Vec<T>,
    pub dom_forward_rates: Vec<T>,
    pub for_forward_rates: Vec<T>,
    pub quoted_dom_basis_spread: T,
    pub dom_rate_vol: T,
    pub fx_spot_vol: T,
    pub corr_dom_rate_fx: T,
}

#[derive(Clone, Debug)]
pub struct BaranBianchettiMtmXccyBasisReport<T: Real = f64> {
    pub mtm_domestic_leg_pv: T,
    pub foreign_leg_pv_in_dom: T,
    pub interim_notional_reset_pv: T,
    pub fx_reset_convexity_adjustment_pv: T,
    pub mtm_swap_npv_in_dom: T,
    pub par_mtm_basis_spread_bps: T,
    pub par_constant_notional_basis_spread_bps: T,
    pub mtm_vs_constant_notional_wedge_bps: T,
}

pub fn evaluate_baran_bianchetti_mtm_xccy_basis_swap<T: Real>(
    s: &BaranBianchettiMtmXccyBasisSpec<T>,
) -> BaranBianchettiMtmXccyBasisReport<T> {
    let n = s
        .payment_times
        .len()
        .min(s.dom_ois_discount_factors.len())
        .min(s.for_xccy_discount_factors.len());
    assert!(n > 0);
    let mut dom_cpn_zero_basis = T::zero();
    let mut const_dom_zero_basis = T::zero();
    let mut for_cpn_pv = T::zero();
    let mut reset = T::zero();
    let mut convex = T::zero();
    let mut ann_mtm = T::zero();
    let mut const_ann = T::zero();
    let mut prev_t = T::zero();
    let mut dd_prev = T::one();
    let mut df_prev = T::one();
    for i in 0..n {
        let t = s.payment_times[i];
        let dt = (t - prev_t).max_val(T::from_f64(1e-8));
        let dd = s.dom_ois_discount_factors[i];
        let df = s.for_xccy_discount_factors[i];
        let rd = s.dom_forward_rates.get(i).copied().unwrap_or(T::zero());
        let rf = s.for_forward_rates.get(i).copied().unwrap_or(T::zero());
        let f_fx = s.spot_fx_dom_per_for * (df_prev / dd_prev.max_val(T::from_f64(1e-12)));
        let mid_t = (prev_t + t) * T::from_f64(0.5);
        let c = (-s.corr_dom_rate_fx * s.dom_rate_vol * s.fx_spot_vol * mid_t).exp();
        let nd = s.foreign_notional * f_fx;
        reset += nd * (dd_prev - dd);
        dom_cpn_zero_basis += nd * dt * rd * dd;
        ann_mtm += nd * c * dt * dd;
        convex += nd * (c - T::one()) * ((dd_prev - dd) + dt * rd * dd);
        for_cpn_pv += s.foreign_notional * s.spot_fx_dom_per_for * dt * rf * df;
        const_dom_zero_basis += s.foreign_notional * s.spot_fx_dom_per_for * dt * rd * dd;
        const_ann += s.foreign_notional * s.spot_fx_dom_per_for * dt * dd;
        prev_t = t;
        dd_prev = dd;
        df_prev = df;
    }
    let foreign_leg =
        for_cpn_pv + s.foreign_notional * s.spot_fx_dom_per_for * (T::one() - df_prev);
    let mtm_dom = dom_cpn_zero_basis + reset + convex + s.quoted_dom_basis_spread * ann_mtm;
    let npv = mtm_dom - foreign_leg;
    let par_mtm =
        (foreign_leg - (dom_cpn_zero_basis + reset + convex)) / ann_mtm.max_val(T::from_f64(1e-12));
    let const_principal = s.foreign_notional * s.spot_fx_dom_per_for * (T::one() - dd_prev);
    let par_const = (foreign_leg - (const_dom_zero_basis + const_principal))
        / const_ann.max_val(T::from_f64(1e-12));
    BaranBianchettiMtmXccyBasisReport {
        mtm_domestic_leg_pv: mtm_dom,
        foreign_leg_pv_in_dom: foreign_leg,
        interim_notional_reset_pv: reset,
        fx_reset_convexity_adjustment_pv: convex,
        mtm_swap_npv_in_dom: npv,
        par_mtm_basis_spread_bps: par_mtm * T::from_f64(1e4),
        par_constant_notional_basis_spread_bps: par_const * T::from_f64(1e4),
        mtm_vs_constant_notional_wedge_bps: (par_mtm - par_const) * T::from_f64(1e4),
    }
}

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
        let guard = TapeScope::guard_with_capacity(8192);
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
        let _ = guard.stop();
        let q = self.basis_jacobian.transpose() * grad;
        (out.val, q.iter().copied().collect())
    }
}

#[derive(Clone, Debug)]
pub struct MtmResettingXccs<T: Real = f64> {
    pub for_notional: T,
    pub spot_fx: T,
    pub payment_times: Vec<f64>,
    pub accrual_fractions: Vec<f64>,
    pub for_basis_spread: T,
    pub dom_spread: T,
    pub receive_foreign: bool,
}

#[derive(Clone, Debug)]
pub struct MtmXccsValuation<T: Real> {
    pub npv: T,
    pub fair_foreign_basis_spread: T,
    pub foreign_leg_pv: T,
    pub resetting_domestic_leg_pv: T,
    pub reset_notional_profile: Vec<T>,
}

pub fn price_mtm_resetting_xccs<T: Real>(
    swap: &MtmResettingXccs<T>,
    dom_ois: &impl YieldTermStructure<T>,
    for_ois: &impl YieldTermStructure<T>,
    for_csa_discount: &impl YieldTermStructure<T>,
) -> MtmXccsValuation<T> {
    let n = swap.payment_times.len().min(swap.accrual_fractions.len());
    if n == 0 {
        return MtmXccsValuation {
            npv: T::zero(),
            fair_foreign_basis_spread: T::zero(),
            foreign_leg_pv: T::zero(),
            resetting_domestic_leg_pv: T::zero(),
            reset_notional_profile: Vec::new(),
        };
    }
    let mut prof = Vec::with_capacity(n);
    let mut dom = T::zero();
    let mut f0 = -swap.for_notional * swap.spot_fx;
    let mut bpv = T::zero();
    let mut prev = 0.0;
    for i in 0..n {
        let t = swap.payment_times[i].max(prev + 1e-8);
        let tau = swap.accrual_fractions[i].max(0.0);
        let ddp = if prev <= 1e-12 {
            T::one()
        } else {
            dom_ois.discount_factor(prev)
        };
        let ddc = dom_ois.discount_factor(t).max_val(T::from_f64(1e-12));
        let fcp = if prev <= 1e-12 {
            T::one()
        } else {
            for_csa_discount.discount_factor(prev)
        };
        let fcc = for_csa_discount.discount_factor(t);
        let nf = swap.for_notional * swap.spot_fx * fcp / ddp.max_val(T::from_f64(1e-12));
        prof.push(nf);
        let dom_flt_tau = ddp / ddc - T::one();
        // Explicit period-by-period MtM principal draw (-nf at prev) + floating coupon + spread + principal repayment at t
        dom += -nf * ddp + nf * (T::one() + dom_flt_tau + swap.dom_spread * tau) * ddc;

        let pp = if prev <= 1e-12 {
            T::one()
        } else {
            for_ois.discount_factor(prev)
        };
        let pc = for_ois.discount_factor(t).max_val(T::from_f64(1e-12));
        f0 += swap.for_notional * swap.spot_fx * (pp / pc - T::one()) * fcc;
        bpv += swap.for_notional * swap.spot_fx * fcc * tau;
        if i + 1 == n {
            f0 += swap.for_notional * swap.spot_fx * fcc;
        }
        prev = t;
    }
    let fair = (dom - f0) / bpv.max_val(T::from_f64(1e-12));
    let fp = f0 + swap.for_basis_spread * bpv;
    let npv = if swap.receive_foreign {
        fp - dom
    } else {
        dom - fp
    };
    MtmXccsValuation {
        npv,
        fair_foreign_basis_spread: fair,
        foreign_leg_pv: fp,
        resetting_domestic_leg_pv: dom,
        reset_notional_profile: prof,
    }
}

#[cfg(test)]
mod domain_suite_30_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;

    #[test]
    fn mtm_xccs_telescoping_fair_spread_and_ad() {
        let d = Date(0);
        let c_dom = FlatForwardCurve::new(d, 0.04);
        let c_for = FlatForwardCurve::new(d, 0.02);
        let c_for_csa = FlatForwardCurve::new(d, 0.025);
        let mut s = MtmResettingXccs {
            for_notional: 1_000_000.0,
            spot_fx: 1.10,
            payment_times: vec![0.5, 1.0, 1.5, 2.0],
            accrual_fractions: vec![0.5; 4],
            for_basis_spread: 0.0,
            dom_spread: 0.0,
            receive_foreign: true,
        };
        let v_same = price_mtm_resetting_xccs(&s, &c_dom, &c_for, &c_for);
        assert!(v_same.resetting_domestic_leg_pv.abs() < 1e-9);
        assert!(v_same.fair_foreign_basis_spread.abs() < 1e-12);

        let v_csa = price_mtm_resetting_xccs(&s, &c_dom, &c_for, &c_for_csa);
        s.for_basis_spread = v_csa.fair_foreign_basis_spread;
        let v_par = price_mtm_resetting_xccs(&s, &c_dom, &c_for, &c_for_csa);
        assert!(v_par.npv.abs() < 1e-9);

        let h = 1e-5;
        let eval_fx = |fx: f64| {
            let mut sp = s.clone();
            sp.for_basis_spread = -0.0015;
            sp.spot_fx = fx;
            price_mtm_resetting_xccs(&sp, &c_dom, &c_for, &c_for_csa).npv
        };
        let fd = (eval_fx(1.10 + h) - eval_fx(1.10 - h)) / (2.0 * h);
        type D = Dual<1>;
        let cd_d = FlatForwardCurve::new(d, D::from_value(0.04));
        let cf_d = FlatForwardCurve::new(d, D::from_value(0.02));
        let cfc_d = FlatForwardCurve::new(d, D::from_value(0.025));
        let sd = MtmResettingXccs {
            for_notional: D::from_value(1_000_000.0),
            spot_fx: D::variable(1.10, 0),
            payment_times: vec![0.5, 1.0, 1.5, 2.0],
            accrual_fractions: vec![0.5; 4],
            for_basis_spread: D::from_value(-0.0015),
            dom_spread: D::zero(),
            receive_foreign: true,
        };
        let ad = price_mtm_resetting_xccs(&sd, &cd_d, &cf_d, &cfc_d).npv;
        assert!((ad.d[0] - fd).abs() < 2e-5);
    }
}

#[cfg(test)]
mod domain_suite_40_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(x: T) -> BaranBianchettiMtmXccyBasisSpec<T> {
        BaranBianchettiMtmXccyBasisSpec {
            spot_fx_dom_per_for: T::from_f64(1.1) + x,
            foreign_notional: T::from_f64(100.0),
            payment_times: vec![T::from_f64(0.5), T::one(), T::from_f64(1.5)],
            dom_ois_discount_factors: vec![T::from_f64(0.99), T::from_f64(0.97), T::from_f64(0.95)],
            for_xccy_discount_factors: vec![
                T::from_f64(0.985),
                T::from_f64(0.96),
                T::from_f64(0.93),
            ],
            dom_forward_rates: vec![T::from_f64(0.03); 3],
            for_forward_rates: vec![T::from_f64(0.02); 3],
            quoted_dom_basis_spread: T::from_f64(0.001),
            dom_rate_vol: T::from_f64(0.2),
            fx_spot_vol: T::from_f64(0.1),
            corr_dom_rate_fx: T::from_f64(-0.3),
        }
    }
    #[test]
    fn mtm_xccy_invariants_and_dual() {
        let mut s = mk(0.0);
        s.corr_dom_rate_fx = 0.0;
        let r = evaluate_baran_bianchetti_mtm_xccy_basis_swap(&s);
        assert!(r.fx_reset_convexity_adjustment_pv.abs() < 1e-12);
        let mut p = mk(0.0);
        let q = evaluate_baran_bianchetti_mtm_xccy_basis_swap(&p).par_mtm_basis_spread_bps * 1e-4;
        p.quoted_dom_basis_spread = q;
        assert!(
            evaluate_baran_bianchetti_mtm_xccy_basis_swap(&p)
                .mtm_swap_npv_in_dom
                .abs()
                < 1e-8
        );
        let h = 1e-5;
        let fd = (evaluate_baran_bianchetti_mtm_xccy_basis_swap(&mk(h)).mtm_domestic_leg_pv
            - evaluate_baran_bianchetti_mtm_xccy_basis_swap(&mk(-h)).mtm_domestic_leg_pv)
            / (2.0 * h);
        let ad = evaluate_baran_bianchetti_mtm_xccy_basis_swap(&mk(Dual::<1>::variable(0.0, 0)))
            .mtm_domestic_leg_pv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

mod evaluate_baba_packer_cip_cross_currency_basis_wedge_module {
    // Baba, Packer & Nagano (BIS, 2008) & Du–Tepper–Verdelhan (2018, JF) Covered Interest
    // Parity (`CIP`) Cross-Currency Basis Wedge, Synthetic USD Rate & G-SIB Balance-Sheet Bound Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BabaPackerCipCrossCurrencyBasisWedgeSpec<T: Real = f64> {
        pub spot_fx_s0: T,
        pub market_outright_forward_fx_f: T,
        pub domestic_usd_ois_rate_rd: T,
        pub foreign_ois_rate_rf: T,
        pub gsib_leverage_ratio_capital_charge_bps: T,
        pub fx_swap_bid_ask_and_funding_haircut_bps: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BabaPackerCipCrossCurrencyBasisWedgeReport<T: Real = f64> {
        pub theoretical_cip_forward_fx: T,
        pub fx_swap_implied_synthetic_usd_rate: T,
        pub covered_interest_parity_basis_bps: T,
        pub dealer_no_arbitrage_band_width_bps: T,
        pub net_executable_cip_arbitrage_bps: T,
        pub forward_points_dislocation_pips: T,
    }

    pub fn evaluate_baba_packer_cip_cross_currency_basis_wedge<T: Real>(
        spec: &BabaPackerCipCrossCurrencyBasisWedgeSpec<T>,
    ) -> BabaPackerCipCrossCurrencyBasisWedgeReport<T> {
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(0.05));
        let f_mkt = spec.market_outright_forward_fx_f.max_val(T::from_f64(0.05));
        let rd = spec.domestic_usd_ois_rate_rd;
        let rf = spec.foreign_ois_rate_rf;
        let cap_bps = spec
            .gsib_leverage_ratio_capital_charge_bps
            .max_val(T::zero());
        let hair_bps = spec
            .fx_swap_bid_ask_and_funding_haircut_bps
            .max_val(T::zero());
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        // Theoretical CIP frictionless forward: F_CIP = S0 * exp((rd - rf) * T)
        let f_cip = s0 * ((rd - rf) * t).exp();

        // Synthetic USD rate implied by FX swap market: r_syn = rf + (1/T) * ln(F_mkt / S0)
        let r_usd_syn = rf + (f_mkt / s0).ln() / t;

        // Du-Tepper-Verdelhan (2018) cross-currency basis x_{t, T} = r_d - r_d^{syn} (in bps):
        let cip_basis_bps = (rd - r_usd_syn) * T::from_f64(10_000.0);
        let band_bps = cap_bps + hair_bps;

        let abs_basis = (cip_basis_bps * cip_basis_bps + T::from_f64(1e-4)).sqrt();
        let raw_arb = abs_basis - band_bps;
        let eps = T::from_f64(1e-2);
        let net_arb_bps = T::from_f64(0.5) * (raw_arb + (raw_arb * raw_arb + eps * eps).sqrt());
        let pips_disloc = (f_mkt - f_cip) * T::from_f64(10_000.0);

        BabaPackerCipCrossCurrencyBasisWedgeReport {
            theoretical_cip_forward_fx: f_cip,
            fx_swap_implied_synthetic_usd_rate: r_usd_syn,
            covered_interest_parity_basis_bps: cip_basis_bps,
            dealer_no_arbitrage_band_width_bps: band_bps,
            net_executable_cip_arbitrage_bps: net_arb_bps,
            forward_points_dislocation_pips: pips_disloc,
        }
    }

    #[cfg(test)]
    mod evaluate_baba_packer_cip_cross_currency_basis_wedge_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn baba_packer_cip_cross_currency_basis_wedge_invariants_and_ad_parity() {
            let spec = BabaPackerCipCrossCurrencyBasisWedgeSpec {
                spot_fx_s0: 1.0800,
                market_outright_forward_fx_f: 1.1015,
                domestic_usd_ois_rate_rd: 0.0480,
                foreign_ois_rate_rf: 0.0310,
                gsib_leverage_ratio_capital_charge_bps: 14.0,
                fx_swap_bid_ask_and_funding_haircut_bps: 4.5,
                maturity_years_t: 1.0,
            };
            let rep = evaluate_baba_packer_cip_cross_currency_basis_wedge(&spec);
            assert!(rep.covered_interest_parity_basis_bps < 0.0);
            assert!(rep.net_executable_cip_arbitrage_bps > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_baba_packer_cip_cross_currency_basis_wedge(
                &BabaPackerCipCrossCurrencyBasisWedgeSpec {
                    spot_fx_s0: Dual::<7>::variable(spec.spot_fx_s0, 0),
                    market_outright_forward_fx_f: Dual::<7>::variable(
                        spec.market_outright_forward_fx_f,
                        1,
                    ),
                    domestic_usd_ois_rate_rd: Dual::<7>::variable(spec.domestic_usd_ois_rate_rd, 2),
                    foreign_ois_rate_rf: Dual::<7>::variable(spec.foreign_ois_rate_rf, 3),
                    gsib_leverage_ratio_capital_charge_bps: Dual::<7>::variable(
                        spec.gsib_leverage_ratio_capital_charge_bps,
                        4,
                    ),
                    fx_swap_bid_ask_and_funding_haircut_bps: Dual::<7>::variable(
                        spec.fx_swap_bid_ask_and_funding_haircut_bps,
                        5,
                    ),
                    maturity_years_t: Dual::<7>::variable(spec.maturity_years_t, 6),
                },
            );
            for (idx, &ad) in r_ad.net_executable_cip_arbitrage_bps.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_fx_s0 += h;
                        sm.spot_fx_s0 -= h;
                    }
                    1 => {
                        sp.market_outright_forward_fx_f += h;
                        sm.market_outright_forward_fx_f -= h;
                    }
                    2 => {
                        sp.domestic_usd_ois_rate_rd += h;
                        sm.domestic_usd_ois_rate_rd -= h;
                    }
                    3 => {
                        sp.foreign_ois_rate_rf += h;
                        sm.foreign_ois_rate_rf -= h;
                    }
                    4 => {
                        sp.gsib_leverage_ratio_capital_charge_bps += h;
                        sm.gsib_leverage_ratio_capital_charge_bps -= h;
                    }
                    5 => {
                        sp.fx_swap_bid_ask_and_funding_haircut_bps += h;
                        sm.fx_swap_bid_ask_and_funding_haircut_bps -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_baba_packer_cip_cross_currency_basis_wedge(&sp)
                    .net_executable_cip_arbitrage_bps
                    - evaluate_baba_packer_cip_cross_currency_basis_wedge(&sm)
                        .net_executable_cip_arbitrage_bps)
                    / (2.0 * h);
                assert!((ad - fd).abs() < 1e-2);
            }
        }
    }
}
pub use evaluate_baba_packer_cip_cross_currency_basis_wedge_module::*;

mod evaluate_wystup_dual_currency_deposit_yield_enhancement_module {
    // Wystup (2017, FX Options and Structured Products, 2nd ed., Ch. 2.7)
    // FX Dual-Currency Deposit (DCD) Yield Enhancement & Embedded Short-Put Strike-Conversion Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct WystupDualCurrencyDepositYieldSpec<T: Real = f64> {
        pub spot_fx_s0: T,
        pub conversion_strike_fx_k: T,
        pub domestic_deposit_rate_rd: T,
        pub foreign_deposit_rate_rf: T,
        pub fx_volatility_sigma: T,
        pub deposit_tenor_years_t: T,
        pub bank_structuring_margin_pct: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct WystupDualCurrencyDepositYieldReport<T: Real = f64> {
        pub embedded_foreign_put_pv_pct: T,
        pub net_investor_option_premium_pct: T,
        pub annualized_yield_enhancement_bps: T,
        pub enhanced_dcd_coupon_rate: T,
        pub risk_neutral_conversion_prob: T,
        pub expected_redemption_value_pct: T,
    }

    pub fn evaluate_wystup_dual_currency_deposit_yield_enhancement<T: Real>(
        spec: &WystupDualCurrencyDepositYieldSpec<T>,
    ) -> WystupDualCurrencyDepositYieldReport<T> {
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-4));
        let k = spec.conversion_strike_fx_k.max_val(T::from_f64(1e-4));
        let rd = spec.domestic_deposit_rate_rd;
        let rf = spec.foreign_deposit_rate_rf;
        let sig = spec.fx_volatility_sigma.max_val(T::from_f64(0.01));
        let t = spec.deposit_tenor_years_t.max_val(T::from_f64(0.02));
        let margin = spec.bank_structuring_margin_pct.max_val(T::zero());

        let std = sig * t.sqrt();
        let d1 = ((s0 / k).ln() + (rd - rf + T::from_f64(0.5) * sig * sig) * t) / std;
        let d2 = d1 - std;

        // In a DCD with domestic principal N_d = 1, the investor sells (1 / K) units of a Garman-Kohlhagen put:
        let gk_put_per_unit =
            k * (-rd * t).exp() * norm_cdf(-d2) - s0 * (-rf * t).exp() * norm_cdf(-d1);
        let put_pv_pct = (gk_put_per_unit / k) * T::from_f64(100.0);
        let net_prem_pct = (put_pv_pct - margin).max_val(T::from_f64(0.01));

        // Future-valued annualized yield enhancement over deposit tenor T:
        let df_d = (-rd * t).exp();
        let yield_enh = (net_prem_pct / T::from_f64(100.0)) / (df_d * t);
        let yield_enh_bps = yield_enh * T::from_f64(10_000.0);
        let enhanced_rate = rd + yield_enh;
        let conv_prob = norm_cdf(-d2);

        let exp_redemption_pct = (T::one() + enhanced_rate * t
            - (put_pv_pct / T::from_f64(100.0)) / df_d)
            * T::from_f64(100.0);

        WystupDualCurrencyDepositYieldReport {
            embedded_foreign_put_pv_pct: put_pv_pct,
            net_investor_option_premium_pct: net_prem_pct,
            annualized_yield_enhancement_bps: yield_enh_bps,
            enhanced_dcd_coupon_rate: enhanced_rate,
            risk_neutral_conversion_prob: conv_prob,
            expected_redemption_value_pct: exp_redemption_pct,
        }
    }

    #[cfg(test)]
    mod evaluate_wystup_dual_currency_deposit_yield_enhancement_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn wystup_dual_currency_deposit_yield_enhancement_invariants_and_ad_parity() {
            let spec = WystupDualCurrencyDepositYieldSpec {
                spot_fx_s0: 1.0850,
                conversion_strike_fx_k: 1.0700,
                domestic_deposit_rate_rd: 0.045,
                foreign_deposit_rate_rf: 0.030,
                fx_volatility_sigma: 0.105,
                deposit_tenor_years_t: 0.25,
                bank_structuring_margin_pct: 0.15,
            };
            let rep = evaluate_wystup_dual_currency_deposit_yield_enhancement(&spec);
            assert!(rep.annualized_yield_enhancement_bps > 50.0);
            assert!(rep.enhanced_dcd_coupon_rate > spec.domestic_deposit_rate_rd);
            assert!(
                rep.risk_neutral_conversion_prob > 0.0 && rep.risk_neutral_conversion_prob < 1.0
            );

            let h = 1e-5;
            let r_ad = evaluate_wystup_dual_currency_deposit_yield_enhancement(
                &WystupDualCurrencyDepositYieldSpec {
                    spot_fx_s0: Dual::<7>::variable(spec.spot_fx_s0, 0),
                    conversion_strike_fx_k: Dual::<7>::variable(spec.conversion_strike_fx_k, 1),
                    domestic_deposit_rate_rd: Dual::<7>::variable(spec.domestic_deposit_rate_rd, 2),
                    foreign_deposit_rate_rf: Dual::<7>::variable(spec.foreign_deposit_rate_rf, 3),
                    fx_volatility_sigma: Dual::<7>::variable(spec.fx_volatility_sigma, 4),
                    deposit_tenor_years_t: Dual::<7>::variable(spec.deposit_tenor_years_t, 5),
                    bank_structuring_margin_pct: Dual::<7>::variable(
                        spec.bank_structuring_margin_pct,
                        6,
                    ),
                },
            );
            for (idx, &ad) in r_ad.enhanced_dcd_coupon_rate.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_fx_s0 += h;
                        sm.spot_fx_s0 -= h;
                    }
                    1 => {
                        sp.conversion_strike_fx_k += h;
                        sm.conversion_strike_fx_k -= h;
                    }
                    2 => {
                        sp.domestic_deposit_rate_rd += h;
                        sm.domestic_deposit_rate_rd -= h;
                    }
                    3 => {
                        sp.foreign_deposit_rate_rf += h;
                        sm.foreign_deposit_rate_rf -= h;
                    }
                    4 => {
                        sp.fx_volatility_sigma += h;
                        sm.fx_volatility_sigma -= h;
                    }
                    5 => {
                        sp.deposit_tenor_years_t += h;
                        sm.deposit_tenor_years_t -= h;
                    }
                    _ => {
                        sp.bank_structuring_margin_pct += h;
                        sm.bank_structuring_margin_pct -= h;
                    }
                }
                let fd = (evaluate_wystup_dual_currency_deposit_yield_enhancement(&sp)
                    .enhanced_dcd_coupon_rate
                    - evaluate_wystup_dual_currency_deposit_yield_enhancement(&sm)
                        .enhanced_dcd_coupon_rate)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_wystup_dual_currency_deposit_yield_enhancement_module::*;
