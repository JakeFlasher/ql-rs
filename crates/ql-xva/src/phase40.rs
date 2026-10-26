//! Gregory (2020) & Basel III/IV `SA-CCR` Replacement Cost (`RC`),
//! Exponential PFE Multiplier & Supervisory Delta-Adjusted `EAD` Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct GregorySaCcrRegulatoryEadSpec<T: Real = f64> {
    pub trade_notional_d: T,
    pub current_mtm_value_v: T,
    pub posted_collateral_c: T,
    pub swap_start_years_s: T,
    pub swap_end_years_e: T,
    pub forward_rate_p: T,
    pub option_strike_k: T,
    pub supervisory_option_vol: T,
    pub supervisory_factor_sf: T,
}

#[derive(Debug, Clone, Copy)]
pub struct GregorySaCcrRegulatoryEadReport<T: Real = f64> {
    pub supervisory_duration_sd: T,
    pub supervisory_option_delta: T,
    pub aggregate_addon: T,
    pub replacement_cost_rc: T,
    pub pfe_overcollateralization_multiplier: T,
    pub potential_future_exposure_pfe: T,
    pub sa_ccr_exposure_at_default_ead: T,
}

pub fn evaluate_gregory_sa_ccr_regulatory_ead<T: Real>(
    spec: &GregorySaCcrRegulatoryEadSpec<T>,
) -> GregorySaCcrRegulatoryEadReport<T> {
    let s = spec.swap_start_years_s.max_val(T::from_f64(0.05));
    let e = spec.swap_end_years_e.max_val(s + T::from_f64(0.25));

    // BCBS 279 supervisory duration SD = (exp(-0.05 * S) - exp(-0.05 * E)) / 0.05
    let r05 = T::from_f64(0.05);
    let sd = ((-r05 * s).exp() - (-r05 * e).exp()) / r05;

    // Supervisory Black-76 option delta for payer swaption:
    let std = spec.supervisory_option_vol.max_val(T::from_f64(1e-4)) * s.sqrt();
    let d1 = ((spec.forward_rate_p.max_val(T::from_f64(1e-4))
        / spec.option_strike_k.max_val(T::from_f64(1e-4)))
    .ln()
        + T::from_f64(0.5) * std * std)
        / std;
    let delta_sup = norm_cdf(d1);

    let addon = spec.trade_notional_d * sd * delta_sup * spec.supervisory_factor_sf;

    // Smooth Replacement Cost RC = max(V - C, 0):
    let net_vm = spec.current_mtm_value_v - spec.posted_collateral_c;
    let rc = T::from_f64(0.5) * (net_vm + (net_vm * net_vm + T::from_f64(1e-8)).sqrt());

    // BCBS 279 exponential PFE multiplier with 5% floor:
    // mult = min(1, 0.05 + 0.95 * exp((V - C) / (2 * 0.95 * AddOn)))
    let floor = T::from_f64(0.05);
    let one_m_floor = T::from_f64(0.95);
    let exp_arg = (net_vm / (T::from_f64(2.0) * one_m_floor * addon.max_val(T::from_f64(1e-8))))
        .min_val(T::zero());
    let mult = floor + one_m_floor * exp_arg.exp();
    let pfe = mult * addon;

    // SA-CCR Regulatory EAD = alpha * (RC + PFE), alpha = 1.4:
    let ead = T::from_f64(1.4) * (rc + pfe);

    GregorySaCcrRegulatoryEadReport {
        supervisory_duration_sd: sd,
        supervisory_option_delta: delta_sup,
        aggregate_addon: addon,
        replacement_cost_rc: rc,
        pfe_overcollateralization_multiplier: mult,
        potential_future_exposure_pfe: pfe,
        sa_ccr_exposure_at_default_ead: ead,
    }
}

#[cfg(test)]
mod phase40_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn gregory_sa_ccr_regulatory_ead_invariants_and_ad_parity() {
        let spec = GregorySaCcrRegulatoryEadSpec {
            trade_notional_d: 1000.0,
            current_mtm_value_v: -8.0,
            posted_collateral_c: 4.0,
            swap_start_years_s: 1.0,
            swap_end_years_e: 6.0,
            forward_rate_p: 0.038,
            option_strike_k: 0.035,
            supervisory_option_vol: 0.50,
            supervisory_factor_sf: 0.005,
        };
        let rep = evaluate_gregory_sa_ccr_regulatory_ead(&spec);
        // Over-collateralized / negative net MtM (V - C < 0) attenuates the PFE multiplier into (0.05, 1.0):
        assert!(
            rep.pfe_overcollateralization_multiplier > 0.05
                && rep.pfe_overcollateralization_multiplier < 1.0
        );
        assert!(rep.sa_ccr_exposure_at_default_ead > 0.0);

        let h = 1e-6;
        let r_ad = evaluate_gregory_sa_ccr_regulatory_ead(&GregorySaCcrRegulatoryEadSpec {
            trade_notional_d: Dual::<1>::from_value(spec.trade_notional_d),
            current_mtm_value_v: Dual::<1>::variable(spec.current_mtm_value_v, 0),
            posted_collateral_c: Dual::<1>::from_value(spec.posted_collateral_c),
            swap_start_years_s: Dual::<1>::from_value(spec.swap_start_years_s),
            swap_end_years_e: Dual::<1>::from_value(spec.swap_end_years_e),
            forward_rate_p: Dual::<1>::from_value(spec.forward_rate_p),
            option_strike_k: Dual::<1>::from_value(spec.option_strike_k),
            supervisory_option_vol: Dual::<1>::from_value(spec.supervisory_option_vol),
            supervisory_factor_sf: Dual::<1>::from_value(spec.supervisory_factor_sf),
        });
        let mut sp = spec;
        sp.current_mtm_value_v += h;
        let mut sm = spec;
        sm.current_mtm_value_v -= h;
        let fd = (evaluate_gregory_sa_ccr_regulatory_ead(&sp).sa_ccr_exposure_at_default_ead
            - evaluate_gregory_sa_ccr_regulatory_ead(&sm).sa_ccr_exposure_at_default_ead)
            / (2.0 * h);
        assert!((r_ad.sa_ccr_exposure_at_default_ead.d[0] - fd).abs() < 1e-5);
    }
}
