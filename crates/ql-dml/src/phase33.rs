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
    for i in 0..n {
        let fit = mx + pi * (s.instrument_z[i] - mz);
        let e = s.endogenous_x[i] - fit;
        rss += e * e;
    }
    let f = T::from_f64((n - 2) as f64) * pi * pi * zz / rss;
    let iv = zy / xz;
    let cf = ols - iv;
    let counter = s.treated_pre_mean + s.control_post_mean - s.control_pre_mean;
    let att = s.treated_post_mean - counter;
    WooldridgeGreene2slsDidReport {
        ols_slope_beta: ols,
        first_stage_instrument_pi1: pi,
        first_stage_f_statistic: f,
        two_stage_ls_iv_beta: iv,
        hausman_endogeneity_wedge: cf,
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
            endogenous_x: vec![1.0, 2.0, 3.0],
            instrument_z: vec![1.0, 2.1, 2.9],
            outcome_y: vec![1.2, 2.1, 3.2],
            treated_pre_mean: 1.0,
            treated_post_mean: 1.2,
            control_pre_mean: 1.0,
            control_post_mean: 1.1,
        };
        assert!(evaluate_wooldridge_greene_2sls_iv_and_did(&s).first_stage_f_statistic > 0.0);
    }
}
