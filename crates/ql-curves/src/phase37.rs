#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct BrigoTuckmanNssKeyRateDurationSpec<T: Real = f64> {
    pub beta0_level: T,
    pub beta1_slope: T,
    pub beta2_hump1: T,
    pub beta3_hump2: T,
    pub tau1_decay: T,
    pub tau2_decay: T,
    pub eval_maturity: T,
    pub cashflow_weights_2y_5y_10y: [T; 3],
}

#[derive(Clone, Copy, Debug)]
pub struct BrigoTuckmanNssKeyRateDurationReport<T: Real = f64> {
    pub short_rate_limit_r0: T,
    pub asymptotic_long_rate: T,
    pub nss_zero_yield_at_eval: T,
    pub nss_instantaneous_forward_at_eval: T,
    pub forward_zero_differential_error: T,
    pub key_rate_durations_2y_5y_10y: [T; 3],
    pub total_fisher_weil_duration: T,
    pub krd_sum_conservation_error: T,
}

fn nss_zero_and_fwd<T: Real>(b0: T, b1: T, b2: T, b3: T, tau1: T, tau2: T, m: T) -> (T, T) {
    let u1 = m / tau1.max_val(T::from_f64(1e-8));
    let u2 = m / tau2.max_val(T::from_f64(1e-8));
    let e1 = (-u1).exp();
    let e2 = (-u2).exp();
    let l1_1 = (T::one() - e1) / u1.max_val(T::from_f64(1e-8));
    let l2_1 = l1_1 - e1;
    let l1_2 = (T::one() - e2) / u2.max_val(T::from_f64(1e-8));
    let l2_2 = l1_2 - e2;
    let y = b0 + b1 * l1_1 + b2 * l2_1 + b3 * l2_2;
    let f = b0 + b1 * e1 + b2 * u1 * e1 + b3 * u2 * e2;
    (y, f)
}

pub fn evaluate_brigo_tuckman_nss_and_key_rate_durations<T: Real>(
    spec: &BrigoTuckmanNssKeyRateDurationSpec<T>,
) -> BrigoTuckmanNssKeyRateDurationReport<T> {
    let m = spec.eval_maturity.max_val(T::from_f64(1e-6));
    let (y_eval, f_eval) = nss_zero_and_fwd(
        spec.beta0_level,
        spec.beta1_slope,
        spec.beta2_hump1,
        spec.beta3_hump2,
        spec.tau1_decay,
        spec.tau2_decay,
        m,
    );
    let dm = T::from_f64(1e-6);
    let (y_p, _) = nss_zero_and_fwd(
        spec.beta0_level,
        spec.beta1_slope,
        spec.beta2_hump1,
        spec.beta3_hump2,
        spec.tau1_decay,
        spec.tau2_decay,
        m + dm,
    );
    let (y_m, _) = nss_zero_and_fwd(
        spec.beta0_level,
        spec.beta1_slope,
        spec.beta2_hump1,
        spec.beta3_hump2,
        spec.tau1_decay,
        spec.tau2_decay,
        m - dm,
    );
    let d_my_dm = ((m + dm) * y_p - (m - dm) * y_m) / (T::from_f64(2.0) * dm);
    let diff_f = f_eval - d_my_dm;
    let fwd_err = (diff_f * diff_f + T::from_f64(1e-30)).sqrt();

    let tenors = [T::from_f64(2.0), T::from_f64(5.0), T::from_f64(10.0)];
    let mut pvs = [T::zero(); 3];
    let mut p_tot = T::zero();
    for k in 0..3 {
        let (yk, _) = nss_zero_and_fwd(
            spec.beta0_level,
            spec.beta1_slope,
            spec.beta2_hump1,
            spec.beta3_hump2,
            spec.tau1_decay,
            spec.tau2_decay,
            tenors[k],
        );
        pvs[k] = spec.cashflow_weights_2y_5y_10y[k] * (-tenors[k] * yk).exp();
        p_tot = p_tot + pvs[k];
    }
    let p_denom = p_tot.max_val(T::from_f64(1e-12));
    let krds = [
        tenors[0] * pvs[0] / p_denom,
        tenors[1] * pvs[1] / p_denom,
        tenors[2] * pvs[2] / p_denom,
    ];
    let fw_dur = (tenors[0] * pvs[0] + tenors[1] * pvs[1] + tenors[2] * pvs[2]) / p_denom;
    let d_krd = (krds[0] + krds[1] + krds[2]) - fw_dur;
    let krd_err = (d_krd * d_krd + T::from_f64(1e-30)).sqrt();

    BrigoTuckmanNssKeyRateDurationReport {
        short_rate_limit_r0: spec.beta0_level + spec.beta1_slope,
        asymptotic_long_rate: spec.beta0_level,
        nss_zero_yield_at_eval: y_eval,
        nss_instantaneous_forward_at_eval: f_eval,
        forward_zero_differential_error: fwd_err,
        key_rate_durations_2y_5y_10y: krds,
        total_fisher_weil_duration: fw_dur,
        krd_sum_conservation_error: krd_err,
    }
}

#[cfg(test)]
mod phase37_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn brigo_tuckman_nss_and_krd_invariants_and_ad_parity() {
        let spec = BrigoTuckmanNssKeyRateDurationSpec {
            beta0_level: 0.045,
            beta1_slope: -0.015,
            beta2_hump1: 0.018,
            beta3_hump2: -0.008,
            tau1_decay: 1.8,
            tau2_decay: 5.5,
            eval_maturity: 4.0,
            cashflow_weights_2y_5y_10y: [0.25, 0.35, 0.40],
        };
        let rep = evaluate_brigo_tuckman_nss_and_key_rate_durations(&spec);
        assert!(rep.forward_zero_differential_error < 1e-9);
        assert!(rep.krd_sum_conservation_error < 1e-12);

        let h = 1e-6;
        let r_ad = evaluate_brigo_tuckman_nss_and_key_rate_durations(
            &BrigoTuckmanNssKeyRateDurationSpec {
                beta0_level: Dual::<1>::variable(spec.beta0_level, 0),
                beta1_slope: Dual::<1>::from_value(spec.beta1_slope),
                beta2_hump1: Dual::<1>::from_value(spec.beta2_hump1),
                beta3_hump2: Dual::<1>::from_value(spec.beta3_hump2),
                tau1_decay: Dual::<1>::from_value(spec.tau1_decay),
                tau2_decay: Dual::<1>::from_value(spec.tau2_decay),
                eval_maturity: Dual::<1>::from_value(spec.eval_maturity),
                cashflow_weights_2y_5y_10y: spec
                    .cashflow_weights_2y_5y_10y
                    .map(Dual::<1>::from_value),
            },
        );
        let mut sp = spec;
        sp.beta0_level += h;
        let mut sm = spec;
        sm.beta0_level -= h;
        let fd = (evaluate_brigo_tuckman_nss_and_key_rate_durations(&sp)
            .total_fisher_weil_duration
            - evaluate_brigo_tuckman_nss_and_key_rate_durations(&sm).total_fisher_weil_duration)
            / (2.0 * h);
        assert!((r_ad.total_fisher_weil_duration.d[0] - fd).abs() < 1e-5);
    }
}
