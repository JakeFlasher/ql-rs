//! Gregory (2020) & Brigo–Morini–Pallavicini (2013) Hull–White Wrong-Way / Right-Way Risk
//! (`WWR`/`RWR`) Hazard-Exposure Co-Movement & CVA Multiplier Engine.

use ql_ad::Real;
use ql_math::{norm_cdf, norm_pdf};

#[derive(Debug, Clone, Copy)]
pub struct GregoryHullWhiteWwrCvaSpec<T: Real = f64> {
    pub baseline_exposure_mu: T,
    pub exposure_volatility_sigma_x: T,
    pub baseline_hazard_rate_lambda0: T,
    pub log_intensity_vol_sigma_lambda: T,
    pub exposure_credit_corr_rho: T,
    pub recovery_rate: T,
    pub risk_free_rate: T,
    pub maturity_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct GregoryHullWhiteWwrCvaReport<T: Real = f64> {
    pub independent_epe: T,
    pub wwr_conditional_epe: T,
    pub independent_cva: T,
    pub wwr_adjusted_cva: T,
    pub wwr_cva_multiplier: T,
    pub effective_wwr_hazard_rate: T,
}

fn bachelier_epe<T: Real>(mean: T, std: T) -> T {
    let s = std.max_val(T::from_f64(1e-8));
    let d = mean / s;
    mean * norm_cdf(d) + s * norm_pdf(d)
}

pub fn evaluate_gregory_hull_white_wrong_way_risk_cva<T: Real>(
    spec: &GregoryHullWhiteWwrCvaSpec<T>,
) -> GregoryHullWhiteWwrCvaReport<T> {
    let t = spec.maturity_years_t.max_val(T::from_f64(0.25));
    let half_t = T::from_f64(0.5) * t;
    let std_x = spec.exposure_volatility_sigma_x * half_t.sqrt();

    // Independent Expected Positive Exposure at representative midpoint t/2:
    let epe_ind = bachelier_epe(spec.baseline_exposure_mu, std_x);

    // Under Hull-White / Brigo lognormal hazard linkage lambda(t) = lambda_0 exp(sigma_lambda W^lambda_t - 0.5 sigma_lambda^2 t^2),
    // conditioning on default at time u shifts the Gaussian exposure drift by
    // Delta_mu(u) = rho_{X,lambda} * sigma_X * sigma_lambda * u:
    let drift_shift = spec.exposure_credit_corr_rho
        * spec.exposure_volatility_sigma_x
        * spec.log_intensity_vol_sigma_lambda
        * half_t;
    let epe_wwr = bachelier_epe(spec.baseline_exposure_mu + drift_shift, std_x);

    let lgd = T::one() - spec.recovery_rate;
    let df = (-spec.risk_free_rate * half_t).exp();
    let pd = T::one() - (-spec.baseline_hazard_rate_lambda0 * t).exp();

    let cva_ind = df * lgd * pd * epe_ind;
    let cva_wwr = df * lgd * pd * epe_wwr;
    let mult = cva_wwr / cva_ind.max_val(T::from_f64(1e-12));
    let eff_lambda = spec.baseline_hazard_rate_lambda0 * mult;

    GregoryHullWhiteWwrCvaReport {
        independent_epe: epe_ind,
        wwr_conditional_epe: epe_wwr,
        independent_cva: cva_ind,
        wwr_adjusted_cva: cva_wwr,
        wwr_cva_multiplier: mult,
        effective_wwr_hazard_rate: eff_lambda,
    }
}

#[cfg(test)]
mod phase39_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn gregory_hull_white_wwr_cva_invariants_and_ad_parity() {
        let spec = GregoryHullWhiteWwrCvaSpec {
            baseline_exposure_mu: 12.0,
            exposure_volatility_sigma_x: 18.0,
            baseline_hazard_rate_lambda0: 0.025,
            log_intensity_vol_sigma_lambda: 0.65,
            exposure_credit_corr_rho: 0.55,
            recovery_rate: 0.40,
            risk_free_rate: 0.035,
            maturity_years_t: 3.0,
        };
        let rep = evaluate_gregory_hull_white_wrong_way_risk_cva(&spec);
        // Positive exposure-credit correlation (Wrong-Way Risk) strictly amplifies EPE and CVA:
        assert!(rep.wwr_conditional_epe > rep.independent_epe);
        assert!(rep.wwr_cva_multiplier > 1.0);

        let h = 1e-6;
        let r_ad = evaluate_gregory_hull_white_wrong_way_risk_cva(&GregoryHullWhiteWwrCvaSpec {
            baseline_exposure_mu: Dual::<1>::from_value(spec.baseline_exposure_mu),
            exposure_volatility_sigma_x: Dual::<1>::from_value(spec.exposure_volatility_sigma_x),
            baseline_hazard_rate_lambda0: Dual::<1>::from_value(spec.baseline_hazard_rate_lambda0),
            log_intensity_vol_sigma_lambda: Dual::<1>::from_value(
                spec.log_intensity_vol_sigma_lambda,
            ),
            exposure_credit_corr_rho: Dual::<1>::variable(spec.exposure_credit_corr_rho, 0),
            recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
            risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
            maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
        });
        let mut sp = spec;
        sp.exposure_credit_corr_rho += h;
        let mut sm = spec;
        sm.exposure_credit_corr_rho -= h;
        let fd = (evaluate_gregory_hull_white_wrong_way_risk_cva(&sp).wwr_adjusted_cva
            - evaluate_gregory_hull_white_wrong_way_risk_cva(&sm).wwr_adjusted_cva)
            / (2.0 * h);
        assert!((r_ad.wwr_adjusted_cva.d[0] - fd).abs() < 1e-6);
    }
}
