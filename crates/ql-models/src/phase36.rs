#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
use ql_ad::Real;

#[inline]
fn smax<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
}

#[inline]
fn smooth_ind<T: Real>(z: T) -> T {
    T::one() / (T::one() + (-z).exp())
}

#[derive(Clone, Copy, Debug)]
pub struct PenmanNissimDuPontReoiValuationSpec<T: Real = f64> {
    pub net_operating_assets_0: T,
    pub net_financial_obligations_0: T,
    pub sales_revenue_1: T,
    pub after_tax_operating_income_1: T,
    pub after_tax_net_financial_expense_1: T,
    pub required_return_operations_rho: T,
    pub perpetual_growth_g: T,
}

#[derive(Clone, Copy, Debug)]
pub struct PenmanNissimDuPontReoiValuationReport<T: Real = f64> {
    pub operating_profit_margin_pm: T,
    pub asset_turnover_ato: T,
    pub return_on_net_operating_assets_rnoa: T,
    pub net_borrowing_cost_nbc: T,
    pub financial_leverage_flev: T,
    pub operating_spread: T,
    pub return_on_common_equity_roe: T,
    pub dupont_roe_conservation_error: T,
    pub residual_operating_income_1: T,
    pub reoi_enterprise_value: T,
    pub fcf_enterprise_value: T,
    pub reoi_fcf_equivalence_error: T,
    pub intrinsic_common_equity_value: T,
}

pub fn evaluate_penman_nissim_dupont_and_reoi_valuation<T: Real>(
    spec: &PenmanNissimDuPontReoiValuationSpec<T>,
) -> PenmanNissimDuPontReoiValuationReport<T> {
    let noa0 = spec.net_operating_assets_0.max_val(T::from_f64(1e-6));
    let nfo0 = spec.net_financial_obligations_0.max_val(T::from_f64(1e-6));
    let cse0 = (noa0 - nfo0).max_val(T::from_f64(1e-6));
    let sales1 = spec.sales_revenue_1.max_val(T::from_f64(1e-6));
    let oi1 = spec.after_tax_operating_income_1;
    let nfe1 = spec.after_tax_net_financial_expense_1;

    let pm = oi1 / sales1;
    let ato = sales1 / noa0;
    let rnoa = pm * ato;
    let nbc = nfe1 / nfo0;
    let flev = nfo0 / cse0;
    let spread = rnoa - nbc;
    let roe_direct = (oi1 - nfe1) / cse0;
    let roe_np = rnoa + flev * spread;
    let d_roe = roe_direct - roe_np;
    let roe_err = (d_roe * d_roe + T::from_f64(1e-30)).sqrt();

    let denom =
        (spec.required_return_operations_rho - spec.perpetual_growth_g).max_val(T::from_f64(1e-6));
    let reoi1 = oi1 - spec.required_return_operations_rho * noa0;
    let v_reoi = noa0 + reoi1 / denom;
    let fcf1 = oi1 - spec.perpetual_growth_g * noa0;
    let v_fcf = fcf1 / denom;
    let d_val = v_reoi - v_fcf;
    let val_err = (d_val * d_val + T::from_f64(1e-30)).sqrt();

    PenmanNissimDuPontReoiValuationReport {
        operating_profit_margin_pm: pm,
        asset_turnover_ato: ato,
        return_on_net_operating_assets_rnoa: rnoa,
        net_borrowing_cost_nbc: nbc,
        financial_leverage_flev: flev,
        operating_spread: spread,
        return_on_common_equity_roe: roe_direct,
        dupont_roe_conservation_error: roe_err,
        residual_operating_income_1: reoi1,
        reoi_enterprise_value: v_reoi,
        fcf_enterprise_value: v_fcf,
        reoi_fcf_equivalence_error: val_err,
        intrinsic_common_equity_value: v_reoi - nfo0,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct LjungqvistSargentMccallSearchSpec<T: Real = f64> {
    pub offer_states: [T; 4],
    pub offer_probabilities: [T; 4],
    pub waiting_compensation_c: T,
    pub discount_factor_beta: T,
}

#[derive(Clone, Copy, Debug)]
pub struct LjungqvistSargentMccallSearchReport<T: Real = f64> {
    pub reservation_threshold_w_bar: T,
    pub bellman_lifetime_value: T,
    pub option_value_of_waiting: T,
    pub acceptance_probability_per_period: T,
    pub expected_stopping_duration_periods: T,
    pub bellman_reservation_residual: T,
}

pub fn evaluate_ljungqvist_sargent_mccall_bellman_search<T: Real>(
    spec: &LjungqvistSargentMccallSearchSpec<T>,
) -> LjungqvistSargentMccallSearchReport<T> {
    let p_sum = (spec.offer_probabilities[0]
        + spec.offer_probabilities[1]
        + spec.offer_probabilities[2]
        + spec.offer_probabilities[3])
        .max_val(T::from_f64(1e-8));
    let p = [
        spec.offer_probabilities[0] / p_sum,
        spec.offer_probabilities[1] / p_sum,
        spec.offer_probabilities[2] / p_sum,
        spec.offer_probabilities[3] / p_sum,
    ];
    let beta = spec.discount_factor_beta.min_val(T::from_f64(0.999));
    let one_m_b = (T::one() - beta).max_val(T::from_f64(1e-6));
    let mut w_bar = p[0] * spec.offer_states[0]
        + p[1] * spec.offer_states[1]
        + p[2] * spec.offer_states[2]
        + p[3] * spec.offer_states[3];

    for _ in 0..80 {
        let mut exp_max = T::zero();
        for i in 0..4 {
            exp_max += p[i] * smax(spec.offer_states[i], w_bar);
        }
        w_bar = one_m_b * spec.waiting_compensation_c + beta * exp_max;
    }

    let mut exp_max_final = T::zero();
    let mut exp_excess = T::zero();
    let mut acc_prob = T::zero();
    for i in 0..4 {
        exp_max_final += p[i] * smax(spec.offer_states[i], w_bar);
        exp_excess += p[i] * smax(spec.offer_states[i] - w_bar, T::zero());
        acc_prob += p[i] * smooth_ind(T::from_f64(30.0) * (spec.offer_states[i] - w_bar));
    }
    let rhs = one_m_b * spec.waiting_compensation_c + beta * exp_max_final;
    let d_res = w_bar - rhs;
    let res_err = (d_res * d_res + T::from_f64(1e-30)).sqrt();

    LjungqvistSargentMccallSearchReport {
        reservation_threshold_w_bar: w_bar,
        bellman_lifetime_value: w_bar / one_m_b,
        option_value_of_waiting: (beta / one_m_b) * exp_excess,
        acceptance_probability_per_period: acc_prob,
        expected_stopping_duration_periods: T::one() / acc_prob.max_val(T::from_f64(1e-6)),
        bellman_reservation_residual: res_err,
    }
}

#[cfg(test)]
mod phase36_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn penman_and_mccall_models_invariants_and_ad_parity() {
        let h = 1e-6;

        let pen_spec = PenmanNissimDuPontReoiValuationSpec {
            net_operating_assets_0: 500.0,
            net_financial_obligations_0: 180.0,
            sales_revenue_1: 650.0,
            after_tax_operating_income_1: 72.0,
            after_tax_net_financial_expense_1: 9.0,
            required_return_operations_rho: 0.09,
            perpetual_growth_g: 0.03,
        };
        let pen_rep = evaluate_penman_nissim_dupont_and_reoi_valuation(&pen_spec);
        assert!(pen_rep.dupont_roe_conservation_error < 1e-12);
        assert!(pen_rep.reoi_fcf_equivalence_error < 1e-10);
        assert!(pen_rep.return_on_common_equity_roe > pen_rep.return_on_net_operating_assets_rnoa);

        let pen_ad = evaluate_penman_nissim_dupont_and_reoi_valuation(
            &PenmanNissimDuPontReoiValuationSpec {
                net_operating_assets_0: Dual::<1>::from_value(pen_spec.net_operating_assets_0),
                net_financial_obligations_0: Dual::<1>::from_value(
                    pen_spec.net_financial_obligations_0,
                ),
                sales_revenue_1: Dual::<1>::from_value(pen_spec.sales_revenue_1),
                after_tax_operating_income_1: Dual::<1>::variable(
                    pen_spec.after_tax_operating_income_1,
                    0,
                ),
                after_tax_net_financial_expense_1: Dual::<1>::from_value(
                    pen_spec.after_tax_net_financial_expense_1,
                ),
                required_return_operations_rho: Dual::<1>::from_value(
                    pen_spec.required_return_operations_rho,
                ),
                perpetual_growth_g: Dual::<1>::from_value(pen_spec.perpetual_growth_g),
            },
        );
        let mut pp = pen_spec;
        pp.after_tax_operating_income_1 += h;
        let mut pm = pen_spec;
        pm.after_tax_operating_income_1 -= h;
        let _fd_pen = (evaluate_penman_nissim_dupont_and_reoi_valuation(&pp)
            .intrinsic_common_equity_value
            - evaluate_penman_nissim_dupont_and_reoi_valuation(&pm).intrinsic_common_equity_value)
            / (2.0 * h);
        assert!((pen_ad.intrinsic_common_equity_value.d[0] - _fd_pen).abs() < 1e-5);

        let mc_spec = LjungqvistSargentMccallSearchSpec {
            offer_states: [40.0, 55.0, 72.0, 95.0],
            offer_probabilities: [0.25, 0.35, 0.25, 0.15],
            waiting_compensation_c: 32.0,
            discount_factor_beta: 0.90,
        };
        let mc_rep = evaluate_ljungqvist_sargent_mccall_bellman_search(&mc_spec);
        assert!(mc_rep.bellman_reservation_residual < 1e-9);
        assert!(mc_rep.reservation_threshold_w_bar > mc_spec.waiting_compensation_c);
        assert!(mc_rep.expected_stopping_duration_periods >= 1.0);

        let mc_ad =
            evaluate_ljungqvist_sargent_mccall_bellman_search(&LjungqvistSargentMccallSearchSpec {
                offer_states: mc_spec.offer_states.map(Dual::<1>::from_value),
                offer_probabilities: mc_spec.offer_probabilities.map(Dual::<1>::from_value),
                waiting_compensation_c: Dual::<1>::variable(mc_spec.waiting_compensation_c, 0),
                discount_factor_beta: Dual::<1>::from_value(mc_spec.discount_factor_beta),
            });
        let mut mp = mc_spec;
        mp.waiting_compensation_c += h;
        let mut mm = mc_spec;
        mm.waiting_compensation_c -= h;
        let _fd_mc = (evaluate_ljungqvist_sargent_mccall_bellman_search(&mp)
            .reservation_threshold_w_bar
            - evaluate_ljungqvist_sargent_mccall_bellman_search(&mm).reservation_threshold_w_bar)
            / (2.0 * h);
        assert!((mc_ad.reservation_threshold_w_bar.d[0] - _fd_mc).abs() < 1e-5);
    }
}
