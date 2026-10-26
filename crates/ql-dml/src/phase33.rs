use ql_ad::Real;
#[derive(Clone, Debug)]
pub struct WooldridgeGreene2slsDidSpec<T: Real = f64> {
    pub endogenous_x: Vec<T>,
    pub instrument_z: Vec<T>,
    pub outcome_y: Vec<T>,
    pub treated_pre_mean: T,
    pub treated_post_mean: T,
    pub control_pre_mean: T,
    pub control_post_mean: T,
}
#[derive(Clone, Copy, Debug)]
pub struct WooldridgeGreene2slsDidReport<T: Real = f64> {
    pub ols_slope_beta: T,
    pub first_stage_instrument_pi1: T,
    pub first_stage_f_statistic: T,
    pub two_stage_ls_iv_beta: T,
    pub hausman_endogeneity_wedge: T,
    pub control_function_residual_coef: T,
    pub did_parallel_trend_counterfactual: T,
    pub did_treatment_effect_att: T,
}
pub fn evaluate_wooldridge_greene_2sls_iv_and_did<T: Real>(
    s: &WooldridgeGreene2slsDidSpec<T>,
) -> WooldridgeGreene2slsDidReport<T> {
    let n = s.endogenous_x.len();
    let mx = s.endogenous_x.iter().copied().fold(T::zero(), |a, x| a + x) / (n as f64);
    let mz = s.instrument_z.iter().copied().fold(T::zero(), |a, x| a + x) / (n as f64);
    let my = s.outcome_y.iter().copied().fold(T::zero(), |a, x| a + x) / (n as f64);
    let mut xx = T::zero();
    let mut zz = T::zero();
    let mut xz = T::zero();
    let mut xy = T::zero();
    let mut zy = T::zero();
    for i in 0..n {
        let x = s.endogenous_x[i] - mx;
        let z = s.instrument_z[i] - mz;
        let y = s.outcome_y[i] - my;
        xx += x * x;
        zz += z * z;
        xz += x * z;
        xy += x * y;
        zy += z * y;
    }
    let ols = xy / xx;
    let pi = xz / zz;
    let mut rss = T::zero();
    let mut vu = T::zero();
    let iv = zy / xz;
    for i in 0..n {
        let fit = mx + pi * (s.instrument_z[i] - mz);
        let vhat = s.endogenous_x[i] - fit;
        let uhat = (s.outcome_y[i] - my) - iv * (s.endogenous_x[i] - mx);
        rss += vhat * vhat;
        vu += vhat * uhat;
    }
    let rss_safe = rss.max_val(T::from_f64(1e-12));
    let f = T::from_f64((n.saturating_sub(2).max(1)) as f64) * pi * pi * zz / rss_safe;
    let cf = vu / rss_safe;
    let counter = s.treated_pre_mean + s.control_post_mean - s.control_pre_mean;
    let att = s.treated_post_mean - counter;
    WooldridgeGreene2slsDidReport {
        ols_slope_beta: ols,
        first_stage_instrument_pi1: pi,
        first_stage_f_statistic: f,
        two_stage_ls_iv_beta: iv,
        hausman_endogeneity_wedge: ols - iv,
        control_function_residual_coef: cf,
        did_parallel_trend_counterfactual: counter,
        did_treatment_effect_att: att,
    }
}
#[cfg(test)]
mod phase33_tests {
    use super::*;
    #[test]
    fn phase33_dml_smoke() {
        let s = WooldridgeGreene2slsDidSpec {
            endogenous_x: vec![1.0, 2.0, 3.0, 4.0],
            instrument_z: vec![1.0, 2.1, 2.9, 4.1],
            outcome_y: vec![1.2, 2.1, 3.2, 4.2],
            treated_pre_mean: 1.0,
            treated_post_mean: 1.2,
            control_pre_mean: 1.0,
            control_post_mean: 1.1,
        };
        let r = evaluate_wooldridge_greene_2sls_iv_and_did(&s);
        assert!(r.first_stage_f_statistic > 10.0);
        assert!(
            (r.did_treatment_effect_att
                - (s.treated_post_mean - r.did_parallel_trend_counterfactual))
                .abs()
                < 1e-12
        );
    }
}
#[cfg(test)]
mod phase33_tests_ad {
    use super::*;
    use ql_ad::Dual;
    fn s<T: Real>(x: T) -> WooldridgeGreene2slsDidSpec<T> {
        WooldridgeGreene2slsDidSpec {
            endogenous_x: vec![
                T::one() + x,
                T::from_f64(2.),
                T::from_f64(3.),
                T::from_f64(4.),
            ],
            instrument_z: vec![
                T::one(),
                T::from_f64(2.1),
                T::from_f64(2.9),
                T::from_f64(4.1),
            ],
            outcome_y: vec![
                T::from_f64(1.2),
                T::from_f64(2.1),
                T::from_f64(3.2),
                T::from_f64(4.2),
            ],
            treated_pre_mean: T::one(),
            treated_post_mean: T::from_f64(1.2),
            control_pre_mean: T::one(),
            control_post_mean: T::from_f64(1.1),
        }
    }
    #[test]
    fn dual_fd() {
        let h = 1e-6;
        let fd = (evaluate_wooldridge_greene_2sls_iv_and_did(&s(h)).two_stage_ls_iv_beta
            - evaluate_wooldridge_greene_2sls_iv_and_did(&s(-h)).two_stage_ls_iv_beta)
            / (2. * h);
        let ad = evaluate_wooldridge_greene_2sls_iv_and_did(&s(Dual::<1>::variable(0., 0)))
            .two_stage_ls_iv_beta;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
