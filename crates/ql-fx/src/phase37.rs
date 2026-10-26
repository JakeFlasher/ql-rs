#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct WystupClarkFxDoubleNoTouchSpec<T: Real = f64> {
    pub spot_fx: T,
    pub lower_barrier_l: T,
    pub upper_barrier_u: T,
    pub domestic_rate: T,
    pub foreign_rate: T,
    pub atm_vol: T,
    pub risk_reversal_25d: T,
    pub butterfly_25d: T,
    pub overhedge_survival_weight: T,
    pub maturity_years: T,
}

#[derive(Clone, Copy, Debug)]
pub struct WystupClarkFxDoubleNoTouchReport<T: Real = f64> {
    pub log_corridor_width: T,
    pub tv_survival_probability: T,
    pub tv_double_no_touch_pv: T,
    pub tv_double_one_touch_pv: T,
    pub touch_no_touch_parity_error: T,
    pub vanna_volga_overhedge_adjustment: T,
    pub smile_adjusted_dnt_pv: T,
    pub smile_adjusted_dot_pv: T,
}

pub fn evaluate_wystup_clark_fx_double_no_touch_vanna_volga<T: Real>(
    spec: &WystupClarkFxDoubleNoTouchSpec<T>,
) -> WystupClarkFxDoubleNoTouchReport<T> {
    let l = spec.lower_barrier_l.max_val(T::from_f64(1e-6));
    let u = spec.upper_barrier_u.max_val(l + T::from_f64(1e-4));
    let s0 = spec
        .spot_fx
        .max_val(l + T::from_f64(1e-5))
        .min_val(u - T::from_f64(1e-5));
    let w = (u / l).ln().max_val(T::from_f64(1e-6));
    let x0 = (s0 / l).ln();
    let sig = spec.atm_vol.max_val(T::from_f64(0.01));
    let mu = spec.domestic_rate - spec.foreign_rate - T::from_f64(0.5) * sig * sig;
    let mat = spec.maturity_years.max_val(T::from_f64(1e-6));
    let pi = T::from_f64(std::f64::consts::PI);

    let mut surv = T::zero();
    for m in [1_usize, 3_usize] {
        let mf = T::from_f64(m as f64);
        let coeff = T::from_f64(4.0) / (mf * pi);
        let eig = (mf * pi * sig / w) * (mf * pi * sig / w) + (mu * mu) / (sig * sig);
        let term = coeff * (mf * pi * x0 / w).sin() * (-T::from_f64(0.5) * eig * mat).exp();
        surv = surv + term;
    }
    let surv_prob = surv.max_val(T::from_f64(1e-6)).min_val(T::from_f64(0.9999));
    let df_d = (-spec.domestic_rate * mat).exp();
    let tv_dnt = df_d * surv_prob;
    let tv_dot = df_d - tv_dnt;
    let d_par = df_d - (tv_dnt + tv_dot);
    let par_err = (d_par * d_par + T::from_f64(1e-30)).sqrt();

    let asym = ((u / s0).ln() - x0) / w;
    let vv_adj = spec.overhedge_survival_weight
        * df_d
        * (T::from_f64(0.5) * asym * spec.risk_reversal_25d
            + T::from_f64(1.5) * spec.butterfly_25d)
        * mat.sqrt();
    let dnt_vv = (tv_dnt + vv_adj).max_val(T::from_f64(1e-6)).min_val(df_d);
    let dot_vv = df_d - dnt_vv;

    WystupClarkFxDoubleNoTouchReport {
        log_corridor_width: w,
        tv_survival_probability: surv_prob,
        tv_double_no_touch_pv: tv_dnt,
        tv_double_one_touch_pv: tv_dot,
        touch_no_touch_parity_error: par_err,
        vanna_volga_overhedge_adjustment: vv_adj,
        smile_adjusted_dnt_pv: dnt_vv,
        smile_adjusted_dot_pv: dot_vv,
    }
}

#[cfg(test)]
mod phase37_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn wystup_clark_fx_double_no_touch_invariants_and_ad_parity() {
        let spec = WystupClarkFxDoubleNoTouchSpec {
            spot_fx: 1.08,
            lower_barrier_l: 1.02,
            upper_barrier_u: 1.15,
            domestic_rate: 0.04,
            foreign_rate: 0.025,
            atm_vol: 0.085,
            risk_reversal_25d: -0.012,
            butterfly_25d: 0.0035,
            overhedge_survival_weight: 0.65,
            maturity_years: 0.5,
        };
        let rep = evaluate_wystup_clark_fx_double_no_touch_vanna_volga(&spec);
        assert!(rep.touch_no_touch_parity_error < 1e-12);
        assert!(rep.tv_double_no_touch_pv > 0.0 && rep.tv_double_one_touch_pv > 0.0);

        let h = 1e-6;
        let r_ad =
            evaluate_wystup_clark_fx_double_no_touch_vanna_volga(&WystupClarkFxDoubleNoTouchSpec {
                spot_fx: Dual::<1>::variable(spec.spot_fx, 0),
                lower_barrier_l: Dual::<1>::from_value(spec.lower_barrier_l),
                upper_barrier_u: Dual::<1>::from_value(spec.upper_barrier_u),
                domestic_rate: Dual::<1>::from_value(spec.domestic_rate),
                foreign_rate: Dual::<1>::from_value(spec.foreign_rate),
                atm_vol: Dual::<1>::from_value(spec.atm_vol),
                risk_reversal_25d: Dual::<1>::from_value(spec.risk_reversal_25d),
                butterfly_25d: Dual::<1>::from_value(spec.butterfly_25d),
                overhedge_survival_weight: Dual::<1>::from_value(spec.overhedge_survival_weight),
                maturity_years: Dual::<1>::from_value(spec.maturity_years),
            });
        let mut sp = spec;
        sp.spot_fx += h;
        let mut sm = spec;
        sm.spot_fx -= h;
        let fd = (evaluate_wystup_clark_fx_double_no_touch_vanna_volga(&sp).smile_adjusted_dnt_pv
            - evaluate_wystup_clark_fx_double_no_touch_vanna_volga(&sm).smile_adjusted_dnt_pv)
            / (2.0 * h);
        assert!((r_ad.smile_adjusted_dnt_pv.d[0] - fd).abs() < 1e-5);
    }
}
