use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct PenmanReoiOhlsonJuettnerSpec<T: Real = f64> {
    pub net_operating_assets_0: T,
    pub net_financial_obligations_0: T,
    pub operating_income_after_tax_1: T,
    pub operating_income_after_tax_2: T,
    pub net_financial_expense_after_tax_1: T,
    pub wacc: T,
    pub terminal_reoi_growth: T,
    pub current_share_price: T,
    pub forward_eps_1: T,
    pub forward_eps_2: T,
    pub forward_dps_1: T,
    pub perpetual_aeg_growth: T,
}

#[derive(Clone, Copy, Debug)]
pub struct PenmanReoiOhlsonJuettnerReport<T: Real = f64> {
    pub return_on_net_operating_assets_rnoa: T,
    pub net_borrowing_cost_nbc: T,
    pub financial_leverage_flev: T,
    pub operating_spread: T,
    pub penman_roe: T,
    pub dupont_identity_conservation_error: T,
    pub residual_operating_income_1: T,
    pub reoi_enterprise_value: T,
    pub reoi_intrinsic_equity_value: T,
    pub ohlson_juettner_implied_cost_of_equity: T,
}

pub fn evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation<T: Real>(
    s: &PenmanReoiOhlsonJuettnerSpec<T>,
) -> PenmanReoiOhlsonJuettnerReport<T> {
    let noa = s.net_operating_assets_0.max_val(T::from_f64(1e-12));
    let nfo = s.net_financial_obligations_0;
    let cse = (noa - nfo).max_val(T::from_f64(1e-12));
    let rnoa = s.operating_income_after_tax_1 / noa;
    let nbc = s.net_financial_expense_after_tax_1 / nfo.abs().max_val(T::from_f64(1e-12));
    let flev = nfo / cse;
    let spread = rnoa - nbc;
    let roe = rnoa + flev * spread;
    let roe_direct = (s.operating_income_after_tax_1 - s.net_financial_expense_after_tax_1) / cse;
    let reoi1 = s.operating_income_after_tax_1 - s.wacc * noa;
    let noa1 = noa * (T::one() + s.terminal_reoi_growth);
    let reoi2 = s.operating_income_after_tax_2 - s.wacc * noa1;
    let rw = s.wacc.max_val(T::from_f64(1e-8));
    let ev = noa
        + reoi1 / (T::one() + rw)
        + reoi2 / ((T::one() + rw) * (rw - s.terminal_reoi_growth).max_val(T::from_f64(1e-8)));
    let g = s.perpetual_aeg_growth;
    let a = (g + s.forward_dps_1 / s.current_share_price.max_val(T::from_f64(1e-12)))
        / T::from_f64(2.0);
    let growth_component =
        (s.forward_eps_2 - s.forward_eps_1) / s.forward_eps_1.max_val(T::from_f64(1e-12)) - g;
    let disc = (a * a
        + s.forward_eps_1 / s.current_share_price.max_val(T::from_f64(1e-12)) * growth_component)
        .max_val(T::zero())
        .sqrt();
    PenmanReoiOhlsonJuettnerReport {
        return_on_net_operating_assets_rnoa: rnoa,
        net_borrowing_cost_nbc: nbc,
        financial_leverage_flev: flev,
        operating_spread: spread,
        penman_roe: roe,
        dupont_identity_conservation_error: (roe_direct - roe).abs(),
        residual_operating_income_1: reoi1,
        reoi_enterprise_value: ev,
        reoi_intrinsic_equity_value: ev - nfo,
        ohlson_juettner_implied_cost_of_equity: a + disc,
    }
}

#[cfg(test)]
mod phase31_tests {
    use super::*;
    use ql_ad::Dual;

    fn spec<T: Real>(x: T) -> PenmanReoiOhlsonJuettnerSpec<T> {
        PenmanReoiOhlsonJuettnerSpec {
            net_operating_assets_0: T::from_f64(100.0) + x,
            net_financial_obligations_0: T::from_f64(30.0),
            operating_income_after_tax_1: T::from_f64(12.0),
            operating_income_after_tax_2: T::from_f64(13.0),
            net_financial_expense_after_tax_1: T::from_f64(2.0),
            wacc: T::from_f64(0.06),
            terminal_reoi_growth: T::from_f64(0.02),
            current_share_price: T::from_f64(100.0),
            forward_eps_1: T::from_f64(5.0),
            forward_eps_2: T::from_f64(5.5),
            forward_dps_1: T::from_f64(2.0),
            perpetual_aeg_growth: T::from_f64(0.03),
        }
    }

    #[test]
    fn phase31_invariants_and_dual_fd() {
        let r = evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation(&spec(0.0));
        assert!(r.dupont_identity_conservation_error < 1e-12);
        assert!(r.ohlson_juettner_implied_cost_of_equity > 0.03);
        let h = 1e-6;
        let fd = (evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation(&spec(h))
            .reoi_enterprise_value
            - evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation(&spec(-h))
                .reoi_enterprise_value)
            / (2.0 * h);
        let ad = evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation(&spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .reoi_enterprise_value;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
