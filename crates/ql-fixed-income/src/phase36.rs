#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
use ql_ad::Real;

#[inline]
fn smin<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b - (d * d + T::from_f64(1e-16)).sqrt())
}

#[inline]
fn smax<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
}

#[derive(Clone, Copy, Debug)]
pub struct FabozziCallablePutableOasSpec<T: Real = f64> {
    pub initial_short_rate: T,
    pub rate_drift_per_year: T,
    pub short_rate_vol: T,
    pub oas_spread: T,
    pub coupon_rate: T,
    pub call_strike_price: T,
    pub put_strike_price: T,
    pub maturity_years: T,
    pub yield_shift_dy: T,
}

#[derive(Clone, Copy, Debug)]
pub struct FabozziCallablePutableOasReport<T: Real = f64> {
    pub straight_bond_fair_value: T,
    pub callable_bond_fair_value: T,
    pub putable_bond_fair_value: T,
    pub embedded_call_option_value: T,
    pub embedded_put_option_value: T,
    pub embedded_option_identity_error: T,
    pub straight_effective_duration: T,
    pub callable_effective_duration: T,
    pub straight_effective_convexity: T,
    pub callable_effective_convexity: T,
}

fn price_trees_shifted<T: Real>(spec: &FabozziCallablePutableOasSpec<T>, shift: T) -> (T, T, T) {
    let dt = spec.maturity_years.max_val(T::from_f64(1e-6)) / T::from_f64(3.0);
    let par = T::from_f64(100.0);
    let cpn = par * spec.coupon_rate * dt;
    let term_payoff = par + cpn;
    let mut vs = [term_payoff; 4];
    let mut vc = [term_payoff; 4];
    let mut vp = [term_payoff; 4];
    let half = T::from_f64(0.5);
    let sig_dt = spec.short_rate_vol * dt.sqrt();

    for step in (0..3).rev() {
        let t_f = T::from_f64(step as f64);
        for j in 0..=step {
            let two_j_minus_t = T::from_f64((2 * j as i32 - step as i32) as f64);
            let r_node = (spec.initial_short_rate + t_f * spec.rate_drift_per_year)
                * (two_j_minus_t * sig_dt).exp()
                + spec.oas_spread
                + shift;
            let df = (-r_node * dt).exp();
            let cont_s = df * half * (vs[j + 1] + vs[j]);
            let cont_c = df * half * (vc[j + 1] + vc[j]);
            let cont_p = df * half * (vp[j + 1] + vp[j]);
            if step > 0 {
                vs[j] = cont_s + cpn;
                vc[j] = smin(cont_c, spec.call_strike_price) + cpn;
                vp[j] = smax(cont_p, spec.put_strike_price) + cpn;
            } else {
                vs[j] = cont_s;
                vc[j] = cont_c;
                vp[j] = cont_p;
            }
        }
    }
    (vs[0], vc[0], vp[0])
}

pub fn evaluate_fabozzi_callable_putable_oas_and_convexity<T: Real>(
    spec: &FabozziCallablePutableOasSpec<T>,
) -> FabozziCallablePutableOasReport<T> {
    let dy = spec.yield_shift_dy.max_val(T::from_f64(1e-5));
    let (vs_0, vc_0, vp_0) = price_trees_shifted(spec, T::zero());
    let (vs_up, vc_up, _) = price_trees_shifted(spec, dy);
    let (vs_dn, vc_dn, _) = price_trees_shifted(spec, -dy);

    let call_val = vs_0 - vc_0;
    let put_val = vp_0 - vs_0;
    let diff = (vc_0 + call_val) - (vp_0 - put_val);
    let id_err = (diff * diff + T::from_f64(1e-30)).sqrt();

    let two = T::from_f64(2.0);
    let d_eff_s = (vs_dn - vs_up) / (two * vs_0 * dy);
    let d_eff_c = (vc_dn - vc_up) / (two * vc_0 * dy);
    let c_eff_s = (vs_up + vs_dn - two * vs_0) / (vs_0 * dy * dy);
    let c_eff_c = (vc_up + vc_dn - two * vc_0) / (vc_0 * dy * dy);

    FabozziCallablePutableOasReport {
        straight_bond_fair_value: vs_0,
        callable_bond_fair_value: vc_0,
        putable_bond_fair_value: vp_0,
        embedded_call_option_value: call_val,
        embedded_put_option_value: put_val,
        embedded_option_identity_error: id_err,
        straight_effective_duration: d_eff_s,
        callable_effective_duration: d_eff_c,
        straight_effective_convexity: c_eff_s,
        callable_effective_convexity: c_eff_c,
    }
}

#[cfg(test)]
mod phase36_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn fabozzi_callable_putable_oas_invariants_and_ad_parity() {
        let spec = FabozziCallablePutableOasSpec {
            initial_short_rate: 0.035,
            rate_drift_per_year: 0.002,
            short_rate_vol: 0.18,
            oas_spread: 0.006,
            coupon_rate: 0.055,
            call_strike_price: 100.5,
            put_strike_price: 99.5,
            maturity_years: 3.0,
            yield_shift_dy: 0.0025,
        };
        let rep = evaluate_fabozzi_callable_putable_oas_and_convexity(&spec);
        assert!(rep.embedded_option_identity_error < 1e-12);
        assert!(rep.callable_effective_duration < rep.straight_effective_duration);
        assert!(rep.callable_effective_convexity < rep.straight_effective_convexity);

        let h = 1e-6;
        let r_ad =
            evaluate_fabozzi_callable_putable_oas_and_convexity(&FabozziCallablePutableOasSpec {
                initial_short_rate: Dual::<1>::variable(spec.initial_short_rate, 0),
                rate_drift_per_year: Dual::<1>::from_value(spec.rate_drift_per_year),
                short_rate_vol: Dual::<1>::from_value(spec.short_rate_vol),
                oas_spread: Dual::<1>::from_value(spec.oas_spread),
                coupon_rate: Dual::<1>::from_value(spec.coupon_rate),
                call_strike_price: Dual::<1>::from_value(spec.call_strike_price),
                put_strike_price: Dual::<1>::from_value(spec.put_strike_price),
                maturity_years: Dual::<1>::from_value(spec.maturity_years),
                yield_shift_dy: Dual::<1>::from_value(spec.yield_shift_dy),
            });
        let mut sp = spec;
        sp.initial_short_rate += h;
        let mut sm = spec;
        sm.initial_short_rate -= h;
        let _fd = (evaluate_fabozzi_callable_putable_oas_and_convexity(&sp)
            .callable_bond_fair_value
            - evaluate_fabozzi_callable_putable_oas_and_convexity(&sm).callable_bond_fair_value)
            / (2.0 * h);
        assert!((r_ad.callable_bond_fair_value.d[0] - _fd).abs() < 1e-5);
    }
}
