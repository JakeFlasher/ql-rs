//! Gregory (2020) & Green (2015) Capital Valuation Adjustment (`KVA`) under SA-CCR
//! & Hurdle-Rate Regulatory Capital Profile Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct GregorySaccrKvaRegulatoryCapitalSpec<T: Real = f64> {
    pub initial_sa_ccr_ead0: T,
    pub counterparty_risk_weight_rw: T,
    pub pillar_1_capital_ratio_k: T,
    pub cva_capital_add_on_multiplier: T,
    pub shareholder_hurdle_rate_hk: T,
    pub risk_free_discount_rate_r: T,
    pub exposure_amortization_decay_mu: T,
    pub portfolio_maturity_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct GregorySaccrKvaRegulatoryCapitalReport<T: Real = f64> {
    pub initial_risk_weighted_assets_rwa0: T,
    pub initial_required_regulatory_capital_k0: T,
    pub default_risk_capital_charge_k_ccr: T,
    pub cva_risk_capital_charge_k_cva: T,
    pub effective_capital_annuity_factor: T,
    pub lifetime_capital_valuation_adjustment_kva: T,
    pub kva_running_spread_equivalent_bps: T,
}

pub fn evaluate_gregory_saccr_kva_regulatory_capital<T: Real>(
    spec: &GregorySaccrKvaRegulatoryCapitalSpec<T>,
) -> GregorySaccrKvaRegulatoryCapitalReport<T> {
    let ead0 = spec.initial_sa_ccr_ead0.max_val(T::from_f64(1e-4));
    let rw = spec.counterparty_risk_weight_rw.max_val(T::from_f64(0.01));
    let cap_ratio = spec.pillar_1_capital_ratio_k.max_val(T::from_f64(0.08));
    let cva_mult = spec.cva_capital_add_on_multiplier.max_val(T::zero());
    let hk = spec.shareholder_hurdle_rate_hk.max_val(T::from_f64(0.01));
    let r = spec.risk_free_discount_rate_r.max_val(T::zero());
    let mu = spec
        .exposure_amortization_decay_mu
        .max_val(T::from_f64(1e-3));
    let t = spec.portfolio_maturity_years_t.max_val(T::from_f64(0.25));

    // Basel RWA = 12.5 * K_CCR = EAD * RW * (1 + CVA_multiplier):
    let rwa0 = ead0 * rw * (T::one() + cva_mult);
    let k_ccr = ead0 * rw * cap_ratio;
    let k_cva = k_ccr * cva_mult;
    let k0 = k_ccr + k_cva;

    // Integrated discounted capital profile over [0, T] under decay rate (r + mu):
    let tot_decay = r + mu;
    let eff_annuity = (T::one() - (-tot_decay * t).exp()) / tot_decay;
    let kva = hk * k0 * eff_annuity;
    let kva_bps = (kva / (ead0 * eff_annuity)) * T::from_f64(1e4);

    GregorySaccrKvaRegulatoryCapitalReport {
        initial_risk_weighted_assets_rwa0: rwa0,
        initial_required_regulatory_capital_k0: k0,
        default_risk_capital_charge_k_ccr: k_ccr,
        cva_risk_capital_charge_k_cva: k_cva,
        effective_capital_annuity_factor: eff_annuity,
        lifetime_capital_valuation_adjustment_kva: kva,
        kva_running_spread_equivalent_bps: kva_bps,
    }
}

#[cfg(test)]
mod phase41_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn gregory_saccr_kva_regulatory_capital_invariants_and_ad_parity() {
        let spec = GregorySaccrKvaRegulatoryCapitalSpec {
            initial_sa_ccr_ead0: 150.0,
            counterparty_risk_weight_rw: 0.75,
            pillar_1_capital_ratio_k: 0.105,
            cva_capital_add_on_multiplier: 0.35,
            shareholder_hurdle_rate_hk: 0.11,
            risk_free_discount_rate_r: 0.038,
            exposure_amortization_decay_mu: 0.18,
            portfolio_maturity_years_t: 5.0,
        };
        let rep = evaluate_gregory_saccr_kva_regulatory_capital(&spec);
        assert!(rep.initial_required_regulatory_capital_k0 > rep.default_risk_capital_charge_k_ccr);
        assert!(rep.lifetime_capital_valuation_adjustment_kva > 0.0);
        assert!(rep.kva_running_spread_equivalent_bps > 0.0);

        let h = 1e-6;
        let r_ad =
            evaluate_gregory_saccr_kva_regulatory_capital(&GregorySaccrKvaRegulatoryCapitalSpec {
                initial_sa_ccr_ead0: Dual::<1>::variable(spec.initial_sa_ccr_ead0, 0),
                counterparty_risk_weight_rw: Dual::<1>::from_value(
                    spec.counterparty_risk_weight_rw,
                ),
                pillar_1_capital_ratio_k: Dual::<1>::from_value(spec.pillar_1_capital_ratio_k),
                cva_capital_add_on_multiplier: Dual::<1>::from_value(
                    spec.cva_capital_add_on_multiplier,
                ),
                shareholder_hurdle_rate_hk: Dual::<1>::from_value(spec.shareholder_hurdle_rate_hk),
                risk_free_discount_rate_r: Dual::<1>::from_value(spec.risk_free_discount_rate_r),
                exposure_amortization_decay_mu: Dual::<1>::from_value(
                    spec.exposure_amortization_decay_mu,
                ),
                portfolio_maturity_years_t: Dual::<1>::from_value(spec.portfolio_maturity_years_t),
            });
        let mut sp = spec;
        sp.initial_sa_ccr_ead0 += h;
        let mut sm = spec;
        sm.initial_sa_ccr_ead0 -= h;
        let fd = (evaluate_gregory_saccr_kva_regulatory_capital(&sp)
            .lifetime_capital_valuation_adjustment_kva
            - evaluate_gregory_saccr_kva_regulatory_capital(&sm)
                .lifetime_capital_valuation_adjustment_kva)
            / (2.0 * h);
        assert!((r_ad.lifetime_capital_valuation_adjustment_kva.d[0] - fd).abs() < 1e-6);
    }
}
