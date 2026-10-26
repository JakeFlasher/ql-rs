#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Clone, Copy, Debug)]
pub struct TuckmanFabozziTipsDeflationFloorSpec<T: Real = f64> {
    pub base_cpi_index: T,
    pub current_cpi_index: T,
    pub nominal_treasury_yield: T,
    pub real_tips_yield: T,
    pub real_coupon_rate: T,
    pub cpi_volatility: T,
    pub inflation_risk_premium: T,
    pub tips_liquidity_premium: T,
    pub maturity_years: T,
}

#[derive(Clone, Copy, Debug)]
pub struct TuckmanFabozziTipsDeflationFloorReport<T: Real = f64> {
    pub canadian_index_ratio: T,
    pub inflation_accreted_principal: T,
    pub exact_fisher_breakeven_inflation: T,
    pub risk_adjusted_expected_inflation: T,
    pub bei_decomposition_error: T,
    pub unfloored_tips_dirty_price: T,
    pub embedded_deflation_floor_put_value: T,
    pub floored_tips_fair_value: T,
    pub floor_adjusted_real_yield: T,
}

pub fn evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei<T: Real>(
    spec: &TuckmanFabozziTipsDeflationFloorSpec<T>,
) -> TuckmanFabozziTipsDeflationFloorReport<T> {
    let ir_t = spec.current_cpi_index / spec.base_cpi_index.max_val(T::from_f64(1e-6));
    let par = T::from_f64(100.0);
    let acc_prin = par * ir_t;

    let bei = (T::one() + spec.nominal_treasury_yield)
        / (T::one() + spec.real_tips_yield).max_val(T::from_f64(1e-6))
        - T::one();
    let exp_inf = bei - spec.inflation_risk_premium + spec.tips_liquidity_premium;
    let d_bei = bei - (exp_inf + spec.inflation_risk_premium - spec.tips_liquidity_premium);
    let bei_err = (d_bei * d_bei + T::from_f64(1e-30)).sqrt();

    let mat = spec.maturity_years.max_val(T::from_f64(1e-6));
    let y_r = spec.real_tips_yield.max_val(T::from_f64(1e-6));
    let df_real = (-y_r * mat).exp();
    let ann_real = (T::one() - df_real) / y_r;
    let unfloored_pv = acc_prin * (spec.real_coupon_rate * ann_real + df_real);

    let f_ir = ir_t * (bei * mat).exp();
    let k_fl = T::one();
    let st = (spec.cpi_volatility * mat.sqrt()).max_val(T::from_f64(1e-6));
    let d1 = ((f_ir / k_fl).ln() + T::from_f64(0.5) * st * st) / st;
    let d2 = d1 - st;
    let df_nom = (-spec.nominal_treasury_yield * mat).exp();
    let floor_put = par * df_nom * (k_fl * norm_cdf(-d2) - f_ir * norm_cdf(-d1));
    let floored_pv = unfloored_pv + floor_put;
    let floor_adj_yield = spec.real_tips_yield - floor_put / (acc_prin * mat);

    TuckmanFabozziTipsDeflationFloorReport {
        canadian_index_ratio: ir_t,
        inflation_accreted_principal: acc_prin,
        exact_fisher_breakeven_inflation: bei,
        risk_adjusted_expected_inflation: exp_inf,
        bei_decomposition_error: bei_err,
        unfloored_tips_dirty_price: unfloored_pv,
        embedded_deflation_floor_put_value: floor_put,
        floored_tips_fair_value: floored_pv,
        floor_adjusted_real_yield: floor_adj_yield,
    }
}

#[cfg(test)]
mod phase37_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn tuckman_fabozzi_tips_deflation_floor_invariants_and_ad_parity() {
        let spec = TuckmanFabozziTipsDeflationFloorSpec {
            base_cpi_index: 295.0,
            current_cpi_index: 300.0,
            nominal_treasury_yield: 0.042,
            real_tips_yield: 0.019,
            real_coupon_rate: 0.0175,
            cpi_volatility: 0.025,
            inflation_risk_premium: 0.003,
            tips_liquidity_premium: 0.0015,
            maturity_years: 5.0,
        };
        let rep = evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei(&spec);
        assert!(rep.bei_decomposition_error < 1e-12);
        assert!(rep.embedded_deflation_floor_put_value > 0.0);
        assert!(rep.floored_tips_fair_value > rep.unfloored_tips_dirty_price);

        let h = 1e-6;
        let r_ad = evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei(
            &TuckmanFabozziTipsDeflationFloorSpec {
                base_cpi_index: Dual::<1>::from_value(spec.base_cpi_index),
                current_cpi_index: Dual::<1>::variable(spec.current_cpi_index, 0),
                nominal_treasury_yield: Dual::<1>::from_value(spec.nominal_treasury_yield),
                real_tips_yield: Dual::<1>::from_value(spec.real_tips_yield),
                real_coupon_rate: Dual::<1>::from_value(spec.real_coupon_rate),
                cpi_volatility: Dual::<1>::from_value(spec.cpi_volatility),
                inflation_risk_premium: Dual::<1>::from_value(spec.inflation_risk_premium),
                tips_liquidity_premium: Dual::<1>::from_value(spec.tips_liquidity_premium),
                maturity_years: Dual::<1>::from_value(spec.maturity_years),
            },
        );
        let mut sp = spec;
        sp.current_cpi_index += h;
        let mut sm = spec;
        sm.current_cpi_index -= h;
        let fd = (evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei(&sp)
            .floored_tips_fair_value
            - evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei(&sm).floored_tips_fair_value)
            / (2.0 * h);
        assert!((r_ad.floored_tips_fair_value.d[0] - fd).abs() < 1e-5);
    }
}
