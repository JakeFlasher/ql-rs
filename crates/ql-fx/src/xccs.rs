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

mod evaluate_wystup_fx_forward_extra_participating_collar_module {
    // Wystup (2017, FX Options and Structured Products, 2nd ed., Ch. 8.3, "Forward Extra")
    // & Clark (2011) FX Forward Extra (Convertible Forward) Participating Collar
    // & Skew-Adjusted Knock-In Put Rebate Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct WystupFxForwardExtraCollarSpec<T: Real = f64> {
        pub spot_fx_s0: T,
        pub worst_case_strike_k: T,
        pub downside_trigger_barrier_b: T,
        pub domestic_rate_rd: T,
        pub foreign_rate_rf: T,
        pub strike_vol_sigma_k: T,
        pub put_wing_skew_slope_xi: T,
        pub expiry_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct WystupFxForwardExtraCollarReport<T: Real = f64> {
        pub outright_fx_forward_rate: T,
        pub garman_kohlhagen_call_at_k: T,
        pub garman_kohlhagen_put_at_k: T,
        pub wing_implied_volatility_at_b: T,
        pub knock_in_short_put_value: T,
        pub forward_extra_client_npv: T,
        pub participation_pickup_over_outright_fwd: T,
    }

    pub fn evaluate_wystup_fx_forward_extra_participating_collar<T: Real>(
        spec: &WystupFxForwardExtraCollarSpec<T>,
    ) -> WystupFxForwardExtraCollarReport<T> {
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(0.05));
        let k = spec.worst_case_strike_k.max_val(T::from_f64(0.05));
        let b = spec
            .downside_trigger_barrier_b
            .max_val(T::from_f64(0.04))
            .min_val(k - T::from_f64(0.005));
        let rd = spec.domestic_rate_rd;
        let rf = spec.foreign_rate_rf;
        let sig_k = spec.strike_vol_sigma_k.max_val(T::from_f64(0.01));
        let xi = spec.put_wing_skew_slope_xi;
        let t = spec.expiry_years_t.max_val(T::from_f64(0.02));

        let sqrt_t = t.sqrt();
        let df_d = (-rd * t).exp();
        let fwd = s0 * ((rd - rf) * t).exp();

        let std_k = sig_k * sqrt_t;
        let d1_k = ((fwd / k).ln() + T::from_f64(0.5) * sig_k * sig_k * t) / std_k;
        let d2_k = d1_k - std_k;
        let call_k = df_d * (fwd * norm_cdf(d1_k) - k * norm_cdf(d2_k));
        let put_k = df_d * (k * norm_cdf(-d2_k) - fwd * norm_cdf(-d1_k));

        let sig_b = (sig_k + xi * (b / k).ln()).max_val(T::from_f64(0.01));
        let std_b = sig_b * sqrt_t;
        let d1_b = ((fwd / b).ln() + T::from_f64(0.5) * sig_b * sig_b * t) / std_b;
        let d2_b = d1_b - std_b;

        let put_b = df_d * (b * norm_cdf(-d2_b) - fwd * norm_cdf(-d1_b));
        let vega_b = df_d * fwd * norm_pdf(d1_b) * sqrt_t;
        let dig_put_b = df_d * (k - b) * norm_cdf(-d2_b) - (k - b) * vega_b * (xi / b);
        let ki_put = put_b + dig_put_b;

        let fe_npv = call_k - ki_put;
        let pickup = put_k - ki_put;

        WystupFxForwardExtraCollarReport {
            outright_fx_forward_rate: fwd,
            garman_kohlhagen_call_at_k: call_k,
            garman_kohlhagen_put_at_k: put_k,
            wing_implied_volatility_at_b: sig_b,
            knock_in_short_put_value: ki_put,
            forward_extra_client_npv: fe_npv,
            participation_pickup_over_outright_fwd: pickup,
        }
    }

    #[cfg(test)]
    mod evaluate_wystup_fx_forward_extra_participating_collar_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn wystup_fx_forward_extra_participating_collar_invariants_and_ad_parity() {
            let spec = WystupFxForwardExtraCollarSpec {
                spot_fx_s0: 1.0850,
                worst_case_strike_k: 1.0950,
                downside_trigger_barrier_b: 1.0350,
                domestic_rate_rd: 0.045,
                foreign_rate_rf: 0.030,
                strike_vol_sigma_k: 0.095,
                put_wing_skew_slope_xi: -0.14,
                expiry_years_t: 0.50,
            };
            let rep = evaluate_wystup_fx_forward_extra_participating_collar(&spec);
            assert!(rep.wing_implied_volatility_at_b > spec.strike_vol_sigma_k);
            assert!(rep.knock_in_short_put_value < rep.garman_kohlhagen_put_at_k);
            assert!(rep.participation_pickup_over_outright_fwd > 0.0);

            let h = 1e-5;
            let r_ad = evaluate_wystup_fx_forward_extra_participating_collar(
                &WystupFxForwardExtraCollarSpec {
                    spot_fx_s0: Dual::<8>::variable(spec.spot_fx_s0, 0),
                    worst_case_strike_k: Dual::<8>::variable(spec.worst_case_strike_k, 1),
                    downside_trigger_barrier_b: Dual::<8>::variable(
                        spec.downside_trigger_barrier_b,
                        2,
                    ),
                    domestic_rate_rd: Dual::<8>::variable(spec.domestic_rate_rd, 3),
                    foreign_rate_rf: Dual::<8>::variable(spec.foreign_rate_rf, 4),
                    strike_vol_sigma_k: Dual::<8>::variable(spec.strike_vol_sigma_k, 5),
                    put_wing_skew_slope_xi: Dual::<8>::variable(spec.put_wing_skew_slope_xi, 6),
                    expiry_years_t: Dual::<8>::variable(spec.expiry_years_t, 7),
                },
            );
            for (idx, &ad) in r_ad.forward_extra_client_npv.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_fx_s0 += h;
                        sm.spot_fx_s0 -= h;
                    }
                    1 => {
                        sp.worst_case_strike_k += h;
                        sm.worst_case_strike_k -= h;
                    }
                    2 => {
                        sp.downside_trigger_barrier_b += h;
                        sm.downside_trigger_barrier_b -= h;
                    }
                    3 => {
                        sp.domestic_rate_rd += h;
                        sm.domestic_rate_rd -= h;
                    }
                    4 => {
                        sp.foreign_rate_rf += h;
                        sm.foreign_rate_rf -= h;
                    }
                    5 => {
                        sp.strike_vol_sigma_k += h;
                        sm.strike_vol_sigma_k -= h;
                    }
                    6 => {
                        sp.put_wing_skew_slope_xi += h;
                        sm.put_wing_skew_slope_xi -= h;
                    }
                    _ => {
                        sp.expiry_years_t += h;
                        sm.expiry_years_t -= h;
                    }
                }
                let fd = (evaluate_wystup_fx_forward_extra_participating_collar(&sp)
                    .forward_extra_client_npv
                    - evaluate_wystup_fx_forward_extra_participating_collar(&sm)
                        .forward_extra_client_npv)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_wystup_fx_forward_extra_participating_collar_module::*;

mod evaluate_reiner_wystup_quanto_fx_correlation_convexity_module {
    // Reiner (1992, Risk 5(3), "Quanto Mechanics"), Wystup (2017, FX Options and Structured
    // Products, 2nd ed., Ch. 4), & Hull (2021, Options, Futures, and Other Derivatives, Ch. 30)
    // Quanto FX Equity Forward, Girsanov Covariance Drift & Quanto Call Option Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct ReinerWystupQuantoFxSpec<T: Real = f64> {
        pub foreign_asset_spot_s0: T,
        pub strike_price_k: T,
        pub fixed_quanto_exchange_rate_q_fx: T,
        pub domestic_rate_rd: T,
        pub foreign_rate_rf: T,
        pub foreign_dividend_yield_q: T,
        pub asset_volatility_sigma_s: T,
        pub fx_volatility_sigma_x: T,
        pub asset_fx_correlation_rho_sx: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct ReinerWystupQuantoFxReport<T: Real = f64> {
        pub unadjusted_foreign_forward_price: T,
        pub quanto_adjusted_forward_price: T,
        pub quanto_covariance_drift_rate_bps: T,
        pub quanto_call_domestic_value: T,
        pub standard_compo_unadjusted_call_value: T,
        pub quanto_forward_discount_premium: T,
    }

    pub fn evaluate_reiner_wystup_quanto_fx_correlation_convexity<T: Real>(
        spec: &ReinerWystupQuantoFxSpec<T>,
    ) -> ReinerWystupQuantoFxReport<T> {
        let s0 = spec.foreign_asset_spot_s0.max_val(T::from_f64(1.0));
        let k = spec.strike_price_k.max_val(T::from_f64(1.0));
        let q_fx = spec
            .fixed_quanto_exchange_rate_q_fx
            .max_val(T::from_f64(0.01));
        let rd = spec.domestic_rate_rd;
        let rf = spec.foreign_rate_rf;
        let q_div = spec.foreign_dividend_yield_q;
        let sig_s = spec.asset_volatility_sigma_s.max_val(T::from_f64(0.01));
        let sig_x = spec.fx_volatility_sigma_x.max_val(T::from_f64(0.01));
        let rho = spec
            .asset_fx_correlation_rho_sx
            .max_val(T::from_f64(-0.98))
            .min_val(T::from_f64(0.98));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let cov_sx = rho * sig_s * sig_x;
        let drift_unadj = rf - q_div;
        let drift_quanto = drift_unadj - cov_sx;

        let f_unadj = s0 * (drift_unadj * t).exp();
        let f_quanto = s0 * (drift_quanto * t).exp();

        let sqrt_t = t.sqrt();
        let std_s = sig_s * sqrt_t;
        let df_d = (-rd * t).exp();

        let d1_q = ((f_quanto / k).ln() + T::from_f64(0.5) * sig_s * sig_s * t) / std_s;
        let d2_q = d1_q - std_s;
        let c_quanto = q_fx * df_d * (f_quanto * norm_cdf(d1_q) - k * norm_cdf(d2_q));

        let d1_u = ((f_unadj / k).ln() + T::from_f64(0.5) * sig_s * sig_s * t) / std_s;
        let d2_u = d1_u - std_s;
        let c_unadj = q_fx * df_d * (f_unadj * norm_cdf(d1_u) - k * norm_cdf(d2_u));

        let cov_bps = -cov_sx * T::from_f64(10_000.0);
        let quanto_prem = c_quanto - c_unadj;

        ReinerWystupQuantoFxReport {
            unadjusted_foreign_forward_price: f_unadj,
            quanto_adjusted_forward_price: f_quanto,
            quanto_covariance_drift_rate_bps: cov_bps,
            quanto_call_domestic_value: c_quanto,
            standard_compo_unadjusted_call_value: c_unadj,
            quanto_forward_discount_premium: quanto_prem,
        }
    }

    #[cfg(test)]
    mod evaluate_reiner_wystup_quanto_fx_correlation_convexity_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn reiner_wystup_quanto_fx_invariants_and_ad_parity() {
            let spec = ReinerWystupQuantoFxSpec {
                foreign_asset_spot_s0: 100.0,
                strike_price_k: 102.0,
                fixed_quanto_exchange_rate_q_fx: 1.0,
                domestic_rate_rd: 0.045,
                foreign_rate_rf: 0.025,
                foreign_dividend_yield_q: 0.015,
                asset_volatility_sigma_s: 0.22,
                fx_volatility_sigma_x: 0.11,
                asset_fx_correlation_rho_sx: 0.45,
                maturity_years_t: 1.0,
            };
            let rep = evaluate_reiner_wystup_quanto_fx_correlation_convexity(&spec);
            assert!(rep.quanto_adjusted_forward_price < rep.unadjusted_foreign_forward_price);
            assert!(rep.quanto_call_domestic_value < rep.standard_compo_unadjusted_call_value);

            let h = 1e-6;
            let r_ad =
                evaluate_reiner_wystup_quanto_fx_correlation_convexity(&ReinerWystupQuantoFxSpec {
                    foreign_asset_spot_s0: Dual::<8>::variable(spec.foreign_asset_spot_s0, 0),
                    strike_price_k: Dual::<8>::variable(spec.strike_price_k, 1),
                    fixed_quanto_exchange_rate_q_fx: Dual::from_f64(
                        spec.fixed_quanto_exchange_rate_q_fx,
                    ),
                    domestic_rate_rd: Dual::<8>::variable(spec.domestic_rate_rd, 2),
                    foreign_rate_rf: Dual::<8>::variable(spec.foreign_rate_rf, 3),
                    foreign_dividend_yield_q: Dual::from_f64(spec.foreign_dividend_yield_q),
                    asset_volatility_sigma_s: Dual::<8>::variable(spec.asset_volatility_sigma_s, 4),
                    fx_volatility_sigma_x: Dual::<8>::variable(spec.fx_volatility_sigma_x, 5),
                    asset_fx_correlation_rho_sx: Dual::<8>::variable(
                        spec.asset_fx_correlation_rho_sx,
                        6,
                    ),
                    maturity_years_t: Dual::<8>::variable(spec.maturity_years_t, 7),
                });
            for i in 0..8 {
                let mut v = [
                    spec.foreign_asset_spot_s0,
                    spec.strike_price_k,
                    spec.domestic_rate_rd,
                    spec.foreign_rate_rf,
                    spec.asset_volatility_sigma_s,
                    spec.fx_volatility_sigma_x,
                    spec.asset_fx_correlation_rho_sx,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_reiner_wystup_quanto_fx_correlation_convexity(
                    &ReinerWystupQuantoFxSpec {
                        foreign_asset_spot_s0: v[0],
                        strike_price_k: v[1],
                        domestic_rate_rd: v[2],
                        foreign_rate_rf: v[3],
                        asset_volatility_sigma_s: v[4],
                        fx_volatility_sigma_x: v[5],
                        asset_fx_correlation_rho_sx: v[6],
                        maturity_years_t: v[7],
                        ..spec
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_reiner_wystup_quanto_fx_correlation_convexity(
                    &ReinerWystupQuantoFxSpec {
                        foreign_asset_spot_s0: v[0],
                        strike_price_k: v[1],
                        domestic_rate_rd: v[2],
                        foreign_rate_rf: v[3],
                        asset_volatility_sigma_s: v[4],
                        fx_volatility_sigma_x: v[5],
                        asset_fx_correlation_rho_sx: v[6],
                        maturity_years_t: v[7],
                        ..spec
                    },
                );
                let fd =
                    (up.quanto_call_domestic_value - dn.quanto_call_domestic_value) / (2.0 * h);
                assert!(
                    (r_ad.quanto_call_domestic_value.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_reiner_wystup_quanto_fx_correlation_convexity_module::*;

/// Specification for the Du–Tepper–Verdelhan (2018) / Andersen–Duffie–Song (2019) Covered Interest
/// Rate Parity (CIP) deviation, synthetic USD Treasury basis, and dealer balance-sheet SLR bound engine.
///
/// Reference: Du, Tepper & Verdelhan (2018, *Journal of Finance* 73(3), "Deviations from Covered
/// Interest Rate Parity"), Andersen, Duffie & Song (2019, *Journal of Finance* 74(1)).
#[derive(Debug, Clone, Copy)]
pub struct DuTepperVerdelhanCipBasisSpec<T: Real> {
    /// Spot FX rate $S_0$ (USD per unit of foreign currency).
    pub spot_fx_s0: T,
    /// Outright FX forward rate $F_{0,\tau}$ at maturity $\tau$.
    pub outright_forward_fx_f0: T,
    /// Domestic (USD) OIS risk-free rate $r_d$.
    pub domestic_usd_ois_rate_rd: T,
    /// Foreign OIS risk-free rate $r_f$.
    pub foreign_ois_rate_rf: T,
    /// Domestic US Treasury yield $y_{\text{UST}}$.
    pub domestic_ust_yield: T,
    /// Foreign government bond yield $y_{\text{ForeignGov}}$.
    pub foreign_gov_yield: T,
    /// Dealer Supplementary Leverage Ratio (SLR) capital charge $\chi_{\text{SLR}} \cdot c_{\text{eq}}$.
    pub dealer_slr_balance_sheet_cost: T,
    /// Dealer unsecured wholesale funding spread $\Delta s_{\text{fund}}$ (Andersen–Duffie–Song FVA wedge).
    pub dealer_funding_spread_delta_s: T,
    /// Horizon $\tau$ (in years).
    pub maturity_years_tau: T,
}

/// Analytics returned by [`evaluate_du_tepper_verdelhan_cip_treasury_basis_arb`].
#[derive(Debug, Clone, Copy)]
pub struct DuTepperVerdelhanCipBasisAnalytics<T: Real> {
    /// Implied OIS cross-currency CIP basis $b_{\text{CIP}}(\tau)$ in basis points.
    pub implied_cip_cross_currency_basis_bps: T,
    /// Du–Tepper–Verdelhan swapped foreign government minus US Treasury convenience premium in bps.
    pub swapped_treasury_convenience_premium_bps: T,
    /// Andersen–Duffie–Song dealer balance-sheet + FVA no-arbitrage hurdle bound in bps.
    pub dealer_slr_fva_no_arb_hurdle_bps: T,
    /// Net actionable CIP arbitrage alpha (basis magnitude minus balance-sheet/funding hurdle) in bps.
    pub net_actionable_cip_arbitrage_alpha_bps: T,
}

/// Evaluates the Du–Tepper–Verdelhan (2018) CIP cross-currency basis, swapped US Treasury convenience
/// yield premium, and Andersen–Duffie–Song (2019) dealer balance-sheet/funding arbitrage bound.
pub fn evaluate_du_tepper_verdelhan_cip_treasury_basis_arb<T: Real>(
    spec: &DuTepperVerdelhanCipBasisSpec<T>,
) -> DuTepperVerdelhanCipBasisAnalytics<T> {
    let tau = spec.maturity_years_tau.max_val(T::from_f64(1e-4));
    let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-6));
    let f0 = spec.outright_forward_fx_f0.max_val(T::from_f64(1e-6));

    // Annualized FX forward premium rho_fx = (1 / tau) * ln(F_{0,tau} / S_0)
    let fwd_prem = (f0 / s0).ln() / tau;

    // Du, Tepper & Verdelhan (2018, Eq. 2): b_CIP = r_d - (r_f + fwd_prem)
    let b_cip = spec.domestic_usd_ois_rate_rd - (spec.foreign_ois_rate_rf + fwd_prem);
    let b_cip_bps = b_cip * T::from_f64(10_000.0);

    // Synthetic USD yield of swapped foreign sovereign bond: y_gov^synth_USD = y_ForeignGov + fwd_prem
    // US Treasury convenience yield premium Phi_UST = y_ForeignGov + fwd_prem - y_UST = (y_ForeignGov - r_f) - (y_UST - r_d) - b_CIP
    let synth_usd_gov_yield = spec.foreign_gov_yield + fwd_prem;
    let ust_conv_prem_bps = (synth_usd_gov_yield - spec.domestic_ust_yield) * T::from_f64(10_000.0);

    // Andersen, Duffie & Song (2019): dealer balance-sheet SLR cost + funding spread + horizon term charge
    let hurdle = spec.dealer_slr_balance_sheet_cost
        + spec.dealer_funding_spread_delta_s * (T::one() + T::from_f64(0.10) * tau.sqrt());
    let hurdle_bps = hurdle * T::from_f64(10_000.0);

    // Smooth absolute CIP basis magnitude minus dealer hurdle + treasury convenience contribution
    let abs_b_cip_bps = (b_cip_bps * b_cip_bps + T::from_f64(0.25)).sqrt();
    let net_alpha_bps = abs_b_cip_bps + T::from_f64(0.20) * ust_conv_prem_bps - hurdle_bps;

    DuTepperVerdelhanCipBasisAnalytics {
        implied_cip_cross_currency_basis_bps: b_cip_bps,
        swapped_treasury_convenience_premium_bps: ust_conv_prem_bps,
        dealer_slr_fva_no_arb_hurdle_bps: hurdle_bps,
        net_actionable_cip_arbitrage_alpha_bps: net_alpha_bps,
    }
}

#[cfg(test)]
mod du_tepper_verdelhan_cip_suite64_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn du_tepper_verdelhan_cip_invariants_and_ad_parity() {
        let base = DuTepperVerdelhanCipBasisSpec {
            spot_fx_s0: 1.0800_f64,
            outright_forward_fx_f0: 1.1015,
            domestic_usd_ois_rate_rd: 0.0480,
            foreign_ois_rate_rf: 0.0310,
            domestic_ust_yield: 0.0465,
            foreign_gov_yield: 0.0325,
            dealer_slr_balance_sheet_cost: 0.0018,
            dealer_funding_spread_delta_s: 0.0012,
            maturity_years_tau: 1.0,
        };
        let r0 = evaluate_du_tepper_verdelhan_cip_treasury_basis_arb(&base);
        assert!(r0.implied_cip_cross_currency_basis_bps.is_finite());
        assert!(r0.swapped_treasury_convenience_premium_bps > 0.0);
        assert!(r0.dealer_slr_fva_no_arb_hurdle_bps > 20.0);

        let inputs = [
            base.spot_fx_s0,
            base.outright_forward_fx_f0,
            base.domestic_usd_ois_rate_rd,
            base.foreign_ois_rate_rf,
            base.domestic_ust_yield,
            base.foreign_gov_yield,
            base.dealer_slr_balance_sheet_cost,
            base.dealer_funding_spread_delta_s,
            base.maturity_years_tau,
        ];
        let mut dual_in = [Dual::<9>::from_f64(0.0); 9];
        for i in 0..9 {
            dual_in[i] = Dual::variable(inputs[i], i);
        }
        let spec_ad = DuTepperVerdelhanCipBasisSpec {
            spot_fx_s0: dual_in[0],
            outright_forward_fx_f0: dual_in[1],
            domestic_usd_ois_rate_rd: dual_in[2],
            foreign_ois_rate_rf: dual_in[3],
            domestic_ust_yield: dual_in[4],
            foreign_gov_yield: dual_in[5],
            dealer_slr_balance_sheet_cost: dual_in[6],
            dealer_funding_spread_delta_s: dual_in[7],
            maturity_years_tau: dual_in[8],
        };
        let r_ad = evaluate_du_tepper_verdelhan_cip_treasury_basis_arb(&spec_ad);
        let eps = 1e-6;
        for i in 0..9 {
            let mut p = inputs;
            let mut m = inputs;
            p[i] += eps;
            m[i] -= eps;
            let mk = |v: [f64; 9]| DuTepperVerdelhanCipBasisSpec {
                spot_fx_s0: v[0],
                outright_forward_fx_f0: v[1],
                domestic_usd_ois_rate_rd: v[2],
                foreign_ois_rate_rf: v[3],
                domestic_ust_yield: v[4],
                foreign_gov_yield: v[5],
                dealer_slr_balance_sheet_cost: v[6],
                dealer_funding_spread_delta_s: v[7],
                maturity_years_tau: v[8],
            };
            let rp = evaluate_du_tepper_verdelhan_cip_treasury_basis_arb(&mk(p));
            let rm = evaluate_du_tepper_verdelhan_cip_treasury_basis_arb(&mk(m));
            let fd = (rp.net_actionable_cip_arbitrage_alpha_bps
                - rm.net_actionable_cip_arbitrage_alpha_bps)
                / (2.0 * eps);
            assert!(
                r_ad.net_actionable_cip_arbitrage_alpha_bps.d[i].abs() > 1e-5,
                "zero derivative at {i}"
            );
            assert!((r_ad.net_actionable_cip_arbitrage_alpha_bps.d[i] - fd).abs() < 1e-3);
        }
    }
}

pub mod fujii_shimada_takahashi_collateralized_xccs_mod {
    use ql_ad::Real;

    /// Specification for the Fujii, Shimada & Takahashi (2010) / Piterbarg (2012)
    /// foreign-currency CSA collateralized cross-currency swap (XCCS) curve discounting engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct FujiiShimadaTakahashiXccsSpec<T: Real> {
        /// Spot FX rate $S_0 > 0$ (domestic USD per unit of foreign currency).
        pub spot_fx_s0: T,
        /// Domestic (USD) OIS continuous rate $r_{d,\text{ois}}$.
        pub domestic_ois_rate_rd: T,
        /// Foreign local OIS continuous rate $r_{f,\text{ois}}$.
        pub foreign_ois_rate_rf: T,
        /// Cross-currency collateral basis spread $b_{\text{xccy}}$ (foreign vs. USD CSA).
        pub xccy_collateral_basis_spread_b: T,
        /// Foreign contractual fixed swap coupon rate $c_f$ (paid in foreign currency).
        pub foreign_swap_coupon_cf: T,
        /// Foreign OIS rate normal volatility $\sigma_f > 0$.
        pub foreign_rate_volatility_sigma_f: T,
        /// Spot FX lognormal volatility $\sigma_{\text{FX}} > 0$.
        pub fx_volatility_sigma_fx: T,
        /// Correlation $\rho_{f,\text{FX}} \in (-1, 1)$ between foreign rate and FX rate.
        pub foreign_fx_correlation_rho: T,
        /// Swap maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Fujii-Shimada-Takahashi collateralized XCCS discounting engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct FujiiShimadaTakahashiXccsAnalytics<T: Real> {
        /// Terminal USD-CSA collateralized foreign discount factor $P_f^{\text{USD-CSA}}(0, T)$.
        pub foreign_usd_csa_discount_factor: T,
        /// Discount factor wedge $P_f^{\text{USD-CSA}}(0, T) - P_f^{\text{OIS}}(0, T)$.
        pub collateral_basis_discount_wedge: T,
        /// CSA-implied terminal FX forward rate $F_{\text{FX}}^{\text{CSA}}(0, T)$.
        pub csa_implied_fx_forward: T,
        /// Present value (in domestic USD per unit foreign notional) of the collateralized XCCS.
        pub collateralized_xccs_present_value: T,
    }

    /// Evaluates the Fujii, Shimada & Takahashi (2010) collateralized foreign curve discount factor
    /// and 4-period cross-currency swap valuation under a domestic USD CSA.
    pub fn evaluate_fujii_shimada_takahashi_collateralized_xccs_curve_discount<T: Real>(
        spec: &FujiiShimadaTakahashiXccsSpec<T>,
    ) -> FujiiShimadaTakahashiXccsAnalytics<T> {
        const PERIODS: usize = 4;
        let t_mat = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let dt = t_mat / T::from_f64(PERIODS as f64);
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-4));

        let df_f_csa = |tm: T| -> T {
            let quanto_conv = T::from_f64(0.5)
                * spec.foreign_fx_correlation_rho
                * spec.foreign_rate_volatility_sigma_f
                * spec.fx_volatility_sigma_fx
                * tm
                * tm;
            (-((spec.foreign_ois_rate_rf + spec.xccy_collateral_basis_spread_b) * tm) + quanto_conv)
                .exp()
        };

        let mut foreign_leg_pv = T::zero();
        for m in 1..=PERIODS {
            let tm = T::from_f64(m as f64) * dt;
            foreign_leg_pv += dt * spec.foreign_swap_coupon_cf * df_f_csa(tm);
        }
        let pf_csa_t = df_f_csa(t_mat);
        foreign_leg_pv += pf_csa_t;

        let pf_local_t = (-spec.foreign_ois_rate_rf * t_mat).exp();
        let pd_ois_t = (-spec.domestic_ois_rate_rd * t_mat).exp();
        let fx_fwd_csa = s0 * pf_csa_t / pd_ois_t;

        // Domestic USD OIS floating leg with notional exchange prices at par (1.0 in USD per $1 notional,
        // scaled by initial spot S0 plus domestic-rate sensitivity on unreset stub)
        let dom_stub_pv = T::one()
            + T::from_f64(0.25) * dt * (spec.domestic_ois_rate_rd - T::from_f64(0.03)) * pd_ois_t;
        let xccs_pv = s0 * (foreign_leg_pv - dom_stub_pv);

        FujiiShimadaTakahashiXccsAnalytics {
            foreign_usd_csa_discount_factor: pf_csa_t,
            collateral_basis_discount_wedge: pf_csa_t - pf_local_t,
            csa_implied_fx_forward: fx_fwd_csa,
            collateralized_xccs_present_value: xccs_pv,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fujii_shimada_takahashi_xccs_invariants_and_ad_parity() {
            let base = FujiiShimadaTakahashiXccsSpec {
                spot_fx_s0: 1.08_f64,
                domestic_ois_rate_rd: 0.045,
                foreign_ois_rate_rf: 0.030,
                xccy_collateral_basis_spread_b: -0.0025,
                foreign_swap_coupon_cf: 0.032,
                foreign_rate_volatility_sigma_f: 0.009,
                fx_volatility_sigma_fx: 0.11,
                foreign_fx_correlation_rho: 0.35,
                maturity_years_t: 4.0,
            };
            let r0 = evaluate_fujii_shimada_takahashi_collateralized_xccs_curve_discount(&base);
            assert!(
                r0.foreign_usd_csa_discount_factor > 0.0
                    && r0.foreign_usd_csa_discount_factor < 1.0
            );
            assert!(r0.collateral_basis_discount_wedge > 0.0);
            assert!(r0.csa_implied_fx_forward > base.spot_fx_s0);

            let inputs = [
                base.spot_fx_s0,
                base.domestic_ois_rate_rd,
                base.foreign_ois_rate_rf,
                base.xccy_collateral_basis_spread_b,
                base.foreign_swap_coupon_cf,
                base.foreign_rate_volatility_sigma_f,
                base.fx_volatility_sigma_fx,
                base.foreign_fx_correlation_rho,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<9>::from_f64(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = FujiiShimadaTakahashiXccsSpec {
                spot_fx_s0: dual_in[0],
                domestic_ois_rate_rd: dual_in[1],
                foreign_ois_rate_rf: dual_in[2],
                xccy_collateral_basis_spread_b: dual_in[3],
                foreign_swap_coupon_cf: dual_in[4],
                foreign_rate_volatility_sigma_f: dual_in[5],
                fx_volatility_sigma_fx: dual_in[6],
                foreign_fx_correlation_rho: dual_in[7],
                maturity_years_t: dual_in[8],
            };
            let r_ad =
                evaluate_fujii_shimada_takahashi_collateralized_xccs_curve_discount(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| FujiiShimadaTakahashiXccsSpec {
                    spot_fx_s0: v[0],
                    domestic_ois_rate_rd: v[1],
                    foreign_ois_rate_rf: v[2],
                    xccy_collateral_basis_spread_b: v[3],
                    foreign_swap_coupon_cf: v[4],
                    foreign_rate_volatility_sigma_f: v[5],
                    fx_volatility_sigma_fx: v[6],
                    foreign_fx_correlation_rho: v[7],
                    maturity_years_t: v[8],
                };
                let rp =
                    evaluate_fujii_shimada_takahashi_collateralized_xccs_curve_discount(&mk(p));
                let rm =
                    evaluate_fujii_shimada_takahashi_collateralized_xccs_curve_discount(&mk(m));
                let fd = (rp.collateralized_xccs_present_value
                    - rm.collateralized_xccs_present_value)
                    / (2.0 * eps);
                assert!(
                    r_ad.collateralized_xccs_present_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.collateralized_xccs_present_value.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use fujii_shimada_takahashi_collateralized_xccs_mod::*;

// ============================================================================
// Baba, Packer & Nagano (2008) & Du, Tepper & Verdelhan (2018) FX Swap Dollar Funding Basis
// ============================================================================

pub mod baba_packer_fx_swap_dollar_basis_mod {
    use ql_ad::Real;

    /// Specification for the Baba, Packer & Nagano (2008, *BIS*) and Du, Tepper & Verdelhan
    /// (2018, *J. Finance* 73(3)) FX-swap-implied synthetic USD funding rate and CIP basis engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BabaPackerFxSwapDollarBasisSpec<T: Real> {
        /// Spot FX rate $S_0 > 0$ (USD per unit of foreign currency).
        pub spot_fx_s0: T,
        /// Outright forward FX rate $F_{0,T} > 0$ at maturity $T$.
        pub outright_forward_fx_f: T,
        /// Direct USD OIS/SOFR rate $r_{\text{USD}}$.
        pub usd_ois_rate_r_usd: T,
        /// Foreign currency OIS rate $r_{\text{FOR}}$.
        pub foreign_ois_rate_r_for: T,
        /// Intermediary G-SIB balance-sheet leverage ratio shadow cost $\lambda_{\text{BS}} \ge 0$.
        pub dealer_balance_sheet_leverage_cost_lambda_bs: T,
        /// Global non-US bank dollar funding demand pressure parameter $d_{\text{USD}} \ge 0$.
        pub foreign_dollar_funding_demand_pressure_d_usd: T,
        /// FX swap tenor $T > 0$ (years).
        pub swap_tenor_years_t: T,
    }

    /// Analytical outputs for the Baba-Packer-Du-Tepper-Verdelhan FX swap dollar basis engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BabaPackerFxSwapDollarBasisAnalytics<T: Real> {
        /// Synthetic FX-swap-implied USD funding rate $r_{\text{USD}}^{\text{FX-imp}}$.
        pub synthetic_fx_swap_implied_usd_rate: T,
        /// Market Covered Interest Parity (CIP) basis spread $x_{\text{CIP}}^{\text{mkt}} = r_{\text{USD}} - r_{\text{USD}}^{\text{FX-imp}}$.
        pub market_cip_cross_currency_basis: T,
        /// Intermediary balance-sheet and funding-demand equilibrium basis wedge.
        pub structural_intermediary_basis_wedge: T,
        /// Discounted all-in dollar funding basis swap present value per unit notional.
        pub balance_sheet_adjusted_dollar_funding_pv: T,
    }

    /// Evaluates the Baba, Packer & Nagano (2008) / Du, Tepper & Verdelhan (2018) synthetic
    /// USD funding cost via FX swaps and intermediary balance-sheet CIP basis wedge.
    pub fn evaluate_baba_packer_nagano_fx_swap_implied_dollar_funding_basis<T: Real>(
        spec: &BabaPackerFxSwapDollarBasisSpec<T>,
    ) -> BabaPackerFxSwapDollarBasisAnalytics<T> {
        let s0 = spec.spot_fx_s0.max_val(T::from_f64(1e-4));
        let fwd = spec.outright_forward_fx_f.max_val(T::from_f64(1e-4));
        let t = spec.swap_tenor_years_t.max_val(T::from_f64(1e-4));

        // Synthetic USD funding rate implied by borrowing foreign currency and swapping into USD
        let r_usd_fx_imp = spec.foreign_ois_rate_r_for + (fwd / s0).ln() / t;
        let cip_basis_mkt = spec.usd_ois_rate_r_usd - r_usd_fx_imp;

        // Du, Tepper & Verdelhan (2018) intermediary balance-sheet leverage + dollar demand wedge
        let structural_wedge = spec.dealer_balance_sheet_leverage_cost_lambda_bs
            + spec.foreign_dollar_funding_demand_pressure_d_usd * t.sqrt();

        let df_usd = (-spec.usd_ois_rate_r_usd * t).exp();
        let funding_pv = df_usd * t * (r_usd_fx_imp - spec.usd_ois_rate_r_usd + structural_wedge);

        BabaPackerFxSwapDollarBasisAnalytics {
            synthetic_fx_swap_implied_usd_rate: r_usd_fx_imp,
            market_cip_cross_currency_basis: cip_basis_mkt,
            structural_intermediary_basis_wedge: structural_wedge,
            balance_sheet_adjusted_dollar_funding_pv: funding_pv,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn baba_packer_fx_swap_dollar_basis_invariants_and_ad_parity() {
            let base = BabaPackerFxSwapDollarBasisSpec {
                spot_fx_s0: 1.085_f64,
                outright_forward_fx_f: 1.096,
                usd_ois_rate_r_usd: 0.048,
                foreign_ois_rate_r_for: 0.032,
                dealer_balance_sheet_leverage_cost_lambda_bs: 0.0022,
                foreign_dollar_funding_demand_pressure_d_usd: 0.0035,
                swap_tenor_years_t: 0.50,
            };
            let r0 = evaluate_baba_packer_nagano_fx_swap_implied_dollar_funding_basis(&base);
            assert!(r0.synthetic_fx_swap_implied_usd_rate > base.foreign_ois_rate_r_for);
            assert!(r0.structural_intermediary_basis_wedge > 0.0);
            assert!(r0.balance_sheet_adjusted_dollar_funding_pv.abs() > 0.0);

            let inputs = [
                base.spot_fx_s0,
                base.outright_forward_fx_f,
                base.usd_ois_rate_r_usd,
                base.foreign_ois_rate_r_for,
                base.dealer_balance_sheet_leverage_cost_lambda_bs,
                base.foreign_dollar_funding_demand_pressure_d_usd,
                base.swap_tenor_years_t,
            ];
            let mut dual_in = [Dual::<7>::from_f64(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BabaPackerFxSwapDollarBasisSpec {
                spot_fx_s0: dual_in[0],
                outright_forward_fx_f: dual_in[1],
                usd_ois_rate_r_usd: dual_in[2],
                foreign_ois_rate_r_for: dual_in[3],
                dealer_balance_sheet_leverage_cost_lambda_bs: dual_in[4],
                foreign_dollar_funding_demand_pressure_d_usd: dual_in[5],
                swap_tenor_years_t: dual_in[6],
            };
            let r_ad = evaluate_baba_packer_nagano_fx_swap_implied_dollar_funding_basis(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| BabaPackerFxSwapDollarBasisSpec {
                    spot_fx_s0: v[0],
                    outright_forward_fx_f: v[1],
                    usd_ois_rate_r_usd: v[2],
                    foreign_ois_rate_r_for: v[3],
                    dealer_balance_sheet_leverage_cost_lambda_bs: v[4],
                    foreign_dollar_funding_demand_pressure_d_usd: v[5],
                    swap_tenor_years_t: v[6],
                };
                let rp = evaluate_baba_packer_nagano_fx_swap_implied_dollar_funding_basis(&mk(p));
                let rm = evaluate_baba_packer_nagano_fx_swap_implied_dollar_funding_basis(&mk(m));
                let fd = (rp.balance_sheet_adjusted_dollar_funding_pv
                    - rm.balance_sheet_adjusted_dollar_funding_pv)
                    / (2.0 * eps);
                assert!(
                    r_ad.balance_sheet_adjusted_dollar_funding_pv.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.balance_sheet_adjusted_dollar_funding_pv.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use baba_packer_fx_swap_dollar_basis_mod::*;
mod gabaix_maggiori_fx_carry_mod {
    use ql_ad::Real;

    /// Input specification for the Gabaix & Maggiori (2015, *QJE* 130(3)) imperfect international
    /// financial markets model of exchange rate determination and FX carry risk premia.
    #[derive(Debug, Clone, Copy)]
    pub struct GabaixMaggioriFxCarrySpec<T: Real> {
        /// Home country gross interest rate $R = 1 + r$ (e.g., $1.045$).
        pub home_gross_interest_rate_r: T,
        /// Foreign country gross interest rate $R^* = 1 + r^*$ (e.g., $1.015$).
        pub foreign_gross_interest_rate_r_star: T,
        /// Financier risk-bearing friction parameter $\gamma_{\text{fin}} > 0$ (e.g., $0.35$).
        pub financier_friction_gamma: T,
        /// Exchange rate return volatility $\sigma_{\text{FX}}$ (e.g., $0.11$).
        pub fx_volatility_sigma: T,
        /// Exogenous global portfolio capital flow shock $f_0$ into home currency (e.g., $0.18$).
        pub capital_flow_demand_f0: T,
        /// Period-0 trade import/export elasticity $\iota_0$ (e.g., $1.0$).
        pub trade_flow_elasticity_iota0: T,
        /// Period-1 fundamental expected trade balance anchor $\iota_1$ (e.g., $1.0$).
        pub period_1_trade_elasticity_iota1: T,
    }

    /// Analytics returned by [`evaluate_gabaix_maggiori_imperfect_financial_markets_fx_carry_risk_premium`].
    #[derive(Debug, Clone, Copy)]
    pub struct GabaixMaggioriFxCarryAnalytics<T: Real> {
        /// Effective financier balance-sheet cost parameter $\Gamma = \gamma_{\text{fin}} \sigma_{\text{FX}}^2$.
        pub financier_risk_penalty_cap_gamma: T,
        /// Equilibrium period-0 exchange rate $e_0$ (foreign currency per unit of home currency).
        pub equilibrium_exchange_rate_e0: T,
        /// Financier-absorbed currency imbalance $q_0$.
        pub financier_absorbed_imbalance_q0: T,
        /// Expected FX carry excess return (UIP deviation) $\mathbb{E}_0[\Omega_1] = \Gamma q_0$.
        pub expected_fx_carry_excess_return: T,
    }

    /// Solves the Gabaix & Maggiori (2015, Proposition 1 & Section III) two-period Gamma equilibrium
    /// where global financiers with limited balance-sheet capacity $\Gamma = \gamma_{\text{fin}}\sigma_{\text{FX}}^2$
    /// absorb net currency demand $q_0$, generating endogenous UIP deviations and FX carry premia.
    pub fn evaluate_gabaix_maggiori_imperfect_financial_markets_fx_carry_risk_premium<T: Real>(
        spec: &GabaixMaggioriFxCarrySpec<T>,
    ) -> GabaixMaggioriFxCarryAnalytics<T> {
        let cap_gamma =
            spec.financier_friction_gamma * spec.fx_volatility_sigma * spec.fx_volatility_sigma;
        let r_home = spec.home_gross_interest_rate_r;
        let r_for = spec.foreign_gross_interest_rate_r_star;
        let i0 = spec.trade_flow_elasticity_iota0;
        let i1 = spec.period_1_trade_elasticity_iota1;
        let f0 = spec.capital_flow_demand_f0;

        // Under Gabaix-Maggiori (2015, Eqs. 11-14):
        // Period-0 balance of payments: i0 + f0 - e0 - q0 = 0 => q0 = i0 + f0 - e0
        // Period-1 balance of payments: i1 - e1 + R * q0 = 0 => E[e1] = i1 + R * q0
        // Financier optimality condition: (R / R^*) * e0 - E[e1] = Gamma * q0
        let r_ratio = r_home / r_for;
        let num = (r_home + cap_gamma) * (i0 + f0) + i1;
        let den = r_ratio + r_home + cap_gamma;
        let e0 = num / den;

        let q0 = i0 + f0 - e0;
        let carry_prem = cap_gamma * q0;

        GabaixMaggioriFxCarryAnalytics {
            financier_risk_penalty_cap_gamma: cap_gamma,
            equilibrium_exchange_rate_e0: e0,
            financier_absorbed_imbalance_q0: q0,
            expected_fx_carry_excess_return: carry_prem,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gabaix_maggiori_fx_carry_invariants_and_ad_parity() {
            let base = GabaixMaggioriFxCarrySpec {
                home_gross_interest_rate_r: 1.045_f64,
                foreign_gross_interest_rate_r_star: 1.015_f64,
                financier_friction_gamma: 0.85_f64,
                fx_volatility_sigma: 0.14_f64,
                capital_flow_demand_f0: 0.22_f64,
                trade_flow_elasticity_iota0: 1.00_f64,
                period_1_trade_elasticity_iota1: 1.05_f64,
            };
            let res =
                evaluate_gabaix_maggiori_imperfect_financial_markets_fx_carry_risk_premium(&base);
            assert!(res.financier_risk_penalty_cap_gamma > 0.0);
            assert!(res.equilibrium_exchange_rate_e0 > 0.5);
            assert!(res.expected_fx_carry_excess_return > 0.0);

            let inputs = [
                base.home_gross_interest_rate_r,
                base.foreign_gross_interest_rate_r_star,
                base.financier_friction_gamma,
                base.fx_volatility_sigma,
                base.capital_flow_demand_f0,
                base.trade_flow_elasticity_iota0,
                base.period_1_trade_elasticity_iota1,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = GabaixMaggioriFxCarrySpec {
                home_gross_interest_rate_r: dual_in[0],
                foreign_gross_interest_rate_r_star: dual_in[1],
                financier_friction_gamma: dual_in[2],
                fx_volatility_sigma: dual_in[3],
                capital_flow_demand_f0: dual_in[4],
                trade_flow_elasticity_iota0: dual_in[5],
                period_1_trade_elasticity_iota1: dual_in[6],
            };
            let r_ad = evaluate_gabaix_maggiori_imperfect_financial_markets_fx_carry_risk_premium(
                &spec_ad,
            );
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| GabaixMaggioriFxCarrySpec {
                    home_gross_interest_rate_r: v[0],
                    foreign_gross_interest_rate_r_star: v[1],
                    financier_friction_gamma: v[2],
                    fx_volatility_sigma: v[3],
                    capital_flow_demand_f0: v[4],
                    trade_flow_elasticity_iota0: v[5],
                    period_1_trade_elasticity_iota1: v[6],
                };
                let rp = evaluate_gabaix_maggiori_imperfect_financial_markets_fx_carry_risk_premium(
                    &mk(p),
                );
                let rm = evaluate_gabaix_maggiori_imperfect_financial_markets_fx_carry_risk_premium(
                    &mk(m),
                );
                let fd = (rp.expected_fx_carry_excess_return - rm.expected_fx_carry_excess_return)
                    / (2.0 * eps);
                assert!(
                    r_ad.expected_fx_carry_excess_return.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.expected_fx_carry_excess_return.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use gabaix_maggiori_fx_carry_mod::*;
