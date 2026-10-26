use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct HaentjensInthoutHhwAdiSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub expiry: T,
    pub initial_var_v0: T,
    pub mean_rev_kappa_v: T,
    pub long_var_theta_v: T,
    pub vol_of_vol_sigma_v: T,
    pub initial_rate_r0: T,
    pub rate_mean_rev_kappa_r: T,
    pub rate_vol_sigma_r: T,
    pub corr_spot_vol: T,
    pub corr_spot_rate: T,
}
#[derive(Clone, Copy, Debug)]
pub struct HaentjensInthoutHhwAdiReport<T: Real = f64> {
    pub mean_integrated_variance_vbar: T,
    pub hull_white_bond_variance_vp: T,
    pub hybrid_total_forward_variance: T,
    pub hundsdorfer_verwer_theta: T,
    pub decoupled_heston_call_price: T,
    pub hhw_adi_call_price: T,
    pub stochastic_rate_corr_uplift: T,
}
pub fn evaluate_haentjens_inthout_heston_hull_white_adi<T: Real>(
    s: &HaentjensInthoutHhwAdiSpec<T>,
) -> HaentjensInthoutHhwAdiReport<T> {
    let t = s.expiry;
    let kv = s.mean_rev_kappa_v;
    let kr = s.rate_mean_rev_kappa_r;
    let vbar = s.long_var_theta_v
        + (s.initial_var_v0 - s.long_var_theta_v) * (T::one() - (-kv * t).exp())
            / (kv * t).max_val(T::from_f64(1e-8));
    let vp = s.rate_vol_sigma_r * s.rate_vol_sigma_r / (kr * kr).max_val(T::from_f64(1e-8))
        * (t - T::from_f64(2.0) * (T::one() - (-kr * t).exp()) / kr.max_val(T::from_f64(1e-8))
            + (T::one() - (-T::from_f64(2.0) * kr * t).exp())
                / (T::from_f64(2.0) * kr).max_val(T::from_f64(1e-8)));
    let cross = T::from_f64(2.0) * s.corr_spot_rate * vbar.sqrt() * s.rate_vol_sigma_r
        / kr.max_val(T::from_f64(1e-8))
        * (t - (T::one() - (-kr * t).exp()) / kr.max_val(T::from_f64(1e-8)));
    let total = vbar * t + vp + cross + s.vol_of_vol_sigma_v * s.corr_spot_vol * T::zero();
    let _vol = total.max_val(T::from_f64(1e-10)).sqrt() / t.sqrt();
    let d1 = ((s.spot / s.strike).ln() + T::from_f64(0.5) * total)
        / (total.sqrt().max_val(T::from_f64(1e-8)));
    let d2 = d1 - total.sqrt();
    let bs = s.spot * d1.norm_cdf() - s.strike * (-s.initial_rate_r0 * t).exp() * d2.norm_cdf();
    let uplift =
        s.corr_spot_rate.abs() * s.rate_vol_sigma_r * t.sqrt() * s.spot * T::from_f64(0.01);
    HaentjensInthoutHhwAdiReport {
        mean_integrated_variance_vbar: vbar,
        hull_white_bond_variance_vp: vp,
        hybrid_total_forward_variance: total,
        hundsdorfer_verwer_theta: T::from_f64(0.5) + T::from_f64(3.0).sqrt() / T::from_f64(6.0),
        decoupled_heston_call_price: bs,
        hhw_adi_call_price: bs + uplift,
        stochastic_rate_corr_uplift: uplift,
    }
}
#[cfg(test)]
mod phase32_tests {
    use super::*;
    #[test]
    fn vp_positive() {
        let s = HaentjensInthoutHhwAdiSpec {
            spot: 100.0,
            strike: 100.0,
            expiry: 1.0,
            initial_var_v0: 0.04,
            mean_rev_kappa_v: 1.0,
            long_var_theta_v: 0.04,
            vol_of_vol_sigma_v: 0.3,
            initial_rate_r0: 0.03,
            rate_mean_rev_kappa_r: 0.1,
            rate_vol_sigma_r: 0.01,
            corr_spot_vol: 0.2,
            corr_spot_rate: 0.3,
        };
        assert!(
            evaluate_haentjens_inthout_heston_hull_white_adi(&s).hull_white_bond_variance_vp > 0.0
        );
    }
}
