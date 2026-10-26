use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct KollerMckinseyKvdEvaValuationSpec<T: Real = f64> {
    pub invested_capital_0: T,
    pub noplat_1: T,
    pub ronic: T,
    pub wacc: T,
    pub perpetual_growth_g: T,
    pub net_debt: T,
    pub diluted_shares: T,
}

#[derive(Clone, Copy, Debug)]
pub struct KollerMckinseyKvdEvaValuationReport<T: Real = f64> {
    pub existing_roic_1: T,
    pub reinvestment_rate_ir: T,
    pub free_cash_flow_1: T,
    pub economic_profit_eva_1: T,
    pub kvd_dcf_enterprise_value: T,
    pub economic_profit_eva_enterprise_value: T,
    pub dcf_eva_equivalence_error: T,
    pub ev_to_invested_capital_multiple: T,
    pub intrinsic_share_price: T,
}

pub fn evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation<T: Real>(
    s: &KollerMckinseyKvdEvaValuationSpec<T>,
) -> KollerMckinseyKvdEvaValuationReport<T> {
    let roic1 = s.noplat_1 / s.invested_capital_0;
    let ir = s.perpetual_growth_g / s.ronic;
    let fcf1 = s.noplat_1 * (T::one() - ir);
    let spread = (s.wacc - s.perpetual_growth_g).max_val(T::from_f64(1e-8));
    let v_dcf = fcf1 / spread;
    let eva1 = s.invested_capital_0 * (roic1 - s.wacc);
    let v_eva = s.invested_capital_0
        + eva1 / s.wacc
        + (s.noplat_1 * ir * (s.ronic - s.wacc)) / (s.wacc * spread);
    let diff = v_dcf - v_eva;
    KollerMckinseyKvdEvaValuationReport {
        existing_roic_1: roic1,
        reinvestment_rate_ir: ir,
        free_cash_flow_1: fcf1,
        economic_profit_eva_1: eva1,
        kvd_dcf_enterprise_value: v_dcf,
        economic_profit_eva_enterprise_value: v_eva,
        dcf_eva_equivalence_error: (diff * diff + T::from_f64(1e-30)).sqrt(),
        ev_to_invested_capital_multiple: v_dcf / s.invested_capital_0,
        intrinsic_share_price: (v_dcf - s.net_debt) / s.diluted_shares,
    }
}

#[cfg(test)]
mod phase35_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> KollerMckinseyKvdEvaValuationSpec<T> {
        KollerMckinseyKvdEvaValuationSpec {
            invested_capital_0: T::from_f64(500.0),
            noplat_1: T::from_f64(75.0) + x,
            ronic: T::from_f64(0.16),
            wacc: T::from_f64(0.09),
            perpetual_growth_g: T::from_f64(0.04),
            net_debt: T::from_f64(150.0),
            diluted_shares: T::from_f64(50.0),
        }
    }
    #[test]
    fn phase35_models_invariants_and_dual_fd() {
        let r = evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation(&spec(0.0));
        assert!(r.dcf_eva_equivalence_error < 1e-10);
        assert!(r.ev_to_invested_capital_multiple > 1.0);
        let h = 1e-6;
        let fd = (evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation(&spec(h))
            .intrinsic_share_price
            - evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation(&spec(-h)).intrinsic_share_price)
            / (2.0 * h);
        let ad =
            evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation(&spec(Dual::<1>::variable(0.0, 0)))
                .intrinsic_share_price;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
