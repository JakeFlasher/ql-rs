#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Clone, Copy, Debug)]
pub struct ContTankovMertonJumpDiffusionSpec<T: Real = f64> {
    pub spot_price: T,
    pub strike_price: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub diffusion_vol_sigma: T,
    pub jump_intensity_lambda: T,
    pub log_jump_mean_mu: T,
    pub log_jump_std_delta: T,
    pub maturity_years: T,
}

#[derive(Clone, Copy, Debug)]
pub struct ContTankovMertonJumpDiffusionReport<T: Real = f64> {
    pub martingale_jump_compensator_k: T,
    pub risk_neutral_jump_intensity: T,
    pub merton_jump_diffusion_call: T,
    pub pure_diffusion_bs_call: T,
    pub jump_tail_option_uplift: T,
    pub total_log_return_variance: T,
    pub jump_induced_excess_kurtosis: T,
}

fn bs_call_mjd<T: Real>(s0: T, k: T, sig: T, r: T, q: T, mat: T) -> T {
    let st = (sig * mat.sqrt()).max_val(T::from_f64(1e-6));
    let d1 = ((s0 / k).ln() + (r - q + T::from_f64(0.5) * sig * sig) * mat) / st;
    let d2 = d1 - st;
    s0 * (-q * mat).exp() * norm_cdf(d1) - k * (-r * mat).exp() * norm_cdf(d2)
}

pub fn evaluate_cont_tankov_merton_jump_diffusion<T: Real>(
    spec: &ContTankovMertonJumpDiffusionSpec<T>,
) -> ContTankovMertonJumpDiffusionReport<T> {
    let s0 = spec.spot_price.max_val(T::from_f64(1e-4));
    let k = spec.strike_price.max_val(T::from_f64(1e-4));
    let r = spec.risk_free_rate;
    let q = spec.dividend_yield;
    let sig = spec.diffusion_vol_sigma.max_val(T::from_f64(0.01));
    let lam = spec.jump_intensity_lambda.max_val(T::from_f64(1e-4));
    let mu_j = spec.log_jump_mean_mu;
    let del_j = spec.log_jump_std_delta.max_val(T::from_f64(1e-4));
    let mat = spec.maturity_years.max_val(T::from_f64(1e-4));

    let k_j = (mu_j + T::from_f64(0.5) * del_j * del_j).exp() - T::one();
    let lam_prime = lam * (T::one() + k_j);
    let ln_1pk = mu_j + T::from_f64(0.5) * del_j * del_j;

    let mut p_n = (-lam_prime * mat).exp();
    let mut p_sum = T::zero();
    let mut call_sum = T::zero();
    let bs_diff = bs_call_mjd(s0, k, sig, r, q, mat);

    for n in 0..7 {
        let nf = T::from_f64(n as f64);
        if n > 0 {
            p_n = p_n * (lam_prime * mat) / nf;
        }
        let r_n = r - lam * k_j + (nf * ln_1pk) / mat;
        let sig_n = (sig * sig + (nf * del_j * del_j) / mat).sqrt();
        let c_n = bs_call_mjd(s0, k, sig_n, r_n, q, mat);
        p_sum = p_sum + p_n;
        call_sum = call_sum + p_n * c_n;
    }
    let mjd_call = call_sum / p_sum.max_val(T::from_f64(1e-12));

    let mu2 = mu_j * mu_j;
    let del2 = del_j * del_j;
    let var_mjd = (sig * sig + lam * (mu2 + del2)) * mat;
    let fourth_jump =
        lam * (mu2 * mu2 + T::from_f64(6.0) * mu2 * del2 + T::from_f64(3.0) * del2 * del2) * mat;
    let kurt_ex = fourth_jump / (var_mjd * var_mjd).max_val(T::from_f64(1e-12));

    ContTankovMertonJumpDiffusionReport {
        martingale_jump_compensator_k: k_j,
        risk_neutral_jump_intensity: lam_prime,
        merton_jump_diffusion_call: mjd_call,
        pure_diffusion_bs_call: bs_diff,
        jump_tail_option_uplift: mjd_call - bs_diff,
        total_log_return_variance: var_mjd,
        jump_induced_excess_kurtosis: kurt_ex,
    }
}

#[cfg(test)]
mod phase38_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn cont_tankov_merton_jump_diffusion_invariants_and_ad_parity() {
        let spec = ContTankovMertonJumpDiffusionSpec {
            spot_price: 100.0,
            strike_price: 105.0,
            risk_free_rate: 0.035,
            dividend_yield: 0.01,
            diffusion_vol_sigma: 0.16,
            jump_intensity_lambda: 0.75,
            log_jump_mean_mu: -0.06,
            log_jump_std_delta: 0.14,
            maturity_years: 1.0,
        };
        let rep = evaluate_cont_tankov_merton_jump_diffusion(&spec);
        assert!(rep.jump_induced_excess_kurtosis > 0.0);
        assert!(rep.merton_jump_diffusion_call > rep.pure_diffusion_bs_call);

        let h = 1e-6;
        let r_ad = evaluate_cont_tankov_merton_jump_diffusion(&ContTankovMertonJumpDiffusionSpec {
            spot_price: Dual::<1>::variable(spec.spot_price, 0),
            strike_price: Dual::<1>::from_value(spec.strike_price),
            risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
            dividend_yield: Dual::<1>::from_value(spec.dividend_yield),
            diffusion_vol_sigma: Dual::<1>::from_value(spec.diffusion_vol_sigma),
            jump_intensity_lambda: Dual::<1>::from_value(spec.jump_intensity_lambda),
            log_jump_mean_mu: Dual::<1>::from_value(spec.log_jump_mean_mu),
            log_jump_std_delta: Dual::<1>::from_value(spec.log_jump_std_delta),
            maturity_years: Dual::<1>::from_value(spec.maturity_years),
        });
        let mut sp = spec;
        sp.spot_price += h;
        let mut sm = spec;
        sm.spot_price -= h;
        let fd = (evaluate_cont_tankov_merton_jump_diffusion(&sp).merton_jump_diffusion_call
            - evaluate_cont_tankov_merton_jump_diffusion(&sm).merton_jump_diffusion_call)
            / (2.0 * h);
        assert!((r_ad.merton_jump_diffusion_call.d[0] - fd).abs() < 1e-5);
    }
}
