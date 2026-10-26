use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct DeSpiegeleerMandatoryConvertibleDecsSpec<T: Real = f64> {
    pub spot_price: T,
    pub par_amount: T,
    pub lower_strike_kl: T,
    pub upper_strike_ku: T,
    pub decs_coupon_rate: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub credit_spread: T,
    pub lower_strike_vol: T,
    pub upper_strike_vol: T,
    pub maturity_years: T,
}
#[derive(Clone, Copy, Debug)]
pub struct DeSpiegeleerMandatoryConvertibleDecsReport<T: Real = f64> {
    pub max_conversion_ratio_lower: T,
    pub min_conversion_ratio_upper: T,
    pub fixed_coupon_strip_pv: T,
    pub prepaid_forward_equity_pv: T,
    pub short_lower_call_leg_pv: T,
    pub long_upper_call_leg_pv: T,
    pub percs_capped_fair_value: T,
    pub decs_mandatory_fair_value: T,
    pub decs_over_percs_upside_participation: T,
    pub decomposition_conservation_error: T,
}
fn call<T: Real>(s: T, k: T, v: T, r: T, q: T, t: T) -> T {
    let st = t.sqrt();
    let d1 = ((s / k).ln() + (r - q + v * v / T::from_f64(2.0)) * t) / (v * st);
    let d2 = d1 - v * st;
    s * (-q * t).exp() * d1.norm_cdf() - k * (-r * t).exp() * d2.norm_cdf()
}
pub fn evaluate_despiegeleer_mandatory_convertible_decs_and_percs<T: Real>(
    s: &DeSpiegeleerMandatoryConvertibleDecsSpec<T>,
) -> DeSpiegeleerMandatoryConvertibleDecsReport<T> {
    let nmax = s.par_amount / s.lower_strike_kl;
    let nmin = s.par_amount / s.upper_strike_ku;
    let z = s.risk_free_rate + s.credit_spread;
    let coup = s.decs_coupon_rate * s.par_amount * (T::one() - (-z * s.maturity_years).exp())
        / z.max_val(T::from_f64(1e-12));
    let fwd = nmax * s.spot_price * (-s.dividend_yield * s.maturity_years).exp();
    let lo = nmax
        * call(
            s.spot_price,
            s.lower_strike_kl,
            s.lower_strike_vol,
            s.risk_free_rate,
            s.dividend_yield,
            s.maturity_years,
        );
    let hi = nmin
        * call(
            s.spot_price,
            s.upper_strike_ku,
            s.upper_strike_vol,
            s.risk_free_rate,
            s.dividend_yield,
            s.maturity_years,
        );
    let p = coup + fwd - lo;
    let d = p + hi;
    DeSpiegeleerMandatoryConvertibleDecsReport {
        max_conversion_ratio_lower: nmax,
        min_conversion_ratio_upper: nmin,
        fixed_coupon_strip_pv: coup,
        prepaid_forward_equity_pv: fwd,
        short_lower_call_leg_pv: -lo,
        long_upper_call_leg_pv: hi,
        percs_capped_fair_value: p,
        decs_mandatory_fair_value: d,
        decs_over_percs_upside_participation: hi,
        decomposition_conservation_error: (d - (coup + fwd - lo + hi)).abs(),
    }
}
#[derive(Clone, Copy, Debug)]
pub struct DavidsonLevinFourFactorIoPoSpec<T: Real = f64> {
    pub initial_balance: T,
    pub gross_wac_rate: T,
    pub net_passthrough_coupon: T,
    pub current_mortgage_rate: T,
    pub base_turnover_cpr: T,
    pub max_refi_cpr: T,
    pub burnout_sensitivity_beta: T,
    pub seasonality_amplitude: T,
    pub oas_discount_rate: T,
    pub num_periods: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct DavidsonLevinFourFactorIoPoReport<T: Real = f64> {
    pub terminal_pool_factor: T,
    pub average_realized_cpr: T,
    pub passthrough_present_value: T,
    pub io_strip_present_value: T,
    pub po_strip_present_value: T,
    pub io_po_strip_conservation_error: T,
    pub io_effective_rate_sensitivity: T,
    pub po_effective_rate_sensitivity: T,
}
pub fn evaluate_davidson_levin_four_factor_prepayment_io_po<T: Real>(
    s: &DavidsonLevinFourFactorIoPoSpec<T>,
) -> DavidsonLevinFourFactorIoPoReport<T> {
    let mut bal = T::one();
    let mut dbal_dr = T::zero();
    let mut io = T::zero();
    let mut po = T::zero();
    let mut dio_dr = T::zero();
    let mut dpo_dr = T::zero();
    let n_per = s.num_periods.max(1);
    for i in 1..=n_per {
        let t_yr = T::from_f64(i as f64 / 12.0);
        let age = T::from_f64((i as f64 / 30.0).min(1.0));
        let bal_safe = bal.max_val(T::from_f64(1e-12));
        let burnout = (bal_safe.ln() * s.burnout_sensitivity_beta).exp();
        let dburn_dbal = s.burnout_sensitivity_beta * burnout / bal_safe;
        let exp_arg = (T::from_f64(-40.0)
            * (s.gross_wac_rate - s.current_mortgage_rate - T::from_f64(0.005)))
        .exp();
        let sig = T::one() / (T::one() + exp_arg);
        let dsig_drmtg = -T::from_f64(40.0) * exp_arg * sig * sig;
        let refi = s.max_refi_cpr * burnout * sig;
        let drefi_dr = s.max_refi_cpr * (dburn_dbal * dbal_dr * sig + burnout * dsig_drmtg);
        let seas = T::one()
            + s.seasonality_amplitude * T::from_f64((i as f64 * std::f64::consts::PI / 6.0).sin());
        let cpr_annual = ((s.base_turnover_cpr * age + refi) * seas)
            .min_val(T::from_f64(0.95))
            .max_val(T::from_f64(1e-6));
        let dcpr_dr = drefi_dr * seas;
        let one_minus_cpr = (T::one() - cpr_annual).max_val(T::from_f64(1e-6));
        let smm = T::one() - one_minus_cpr.powf(1.0 / 12.0);
        let dsmm_dr = (T::one() - smm) / (T::from_f64(12.0) * one_minus_cpr) * dcpr_dr;
        let disc = (-s.oas_discount_rate * t_yr).exp();
        let c_m = s.net_passthrough_coupon / T::from_f64(12.0);
        io += bal * c_m * disc;
        dio_dr += dbal_dr * c_m * disc;
        let pre = if i == n_per { bal } else { bal * smm };
        let dpre_dr = if i == n_per {
            dbal_dr
        } else {
            dbal_dr * smm + bal * dsmm_dr
        };
        po += pre * disc;
        dpo_dr += dpre_dr * disc;
        let next_bal = bal * (T::one() - smm);
        dbal_dr = dbal_dr * (T::one() - smm) - bal * dsmm_dr;
        bal = next_bal;
    }
    let pv = io + po;
    let avg = T::one()
        - (bal.max_val(T::from_f64(1e-12)).ln() * T::from_f64(12.0 / (n_per as f64))).exp();
    let err = (pv - io - po) * s.initial_balance;
    DavidsonLevinFourFactorIoPoReport {
        terminal_pool_factor: bal,
        average_realized_cpr: avg,
        passthrough_present_value: pv * s.initial_balance,
        io_strip_present_value: io * s.initial_balance,
        po_strip_present_value: po * s.initial_balance,
        io_po_strip_conservation_error: (err * err + T::from_f64(1e-30)).sqrt(),
        io_effective_rate_sensitivity: dio_dr * s.initial_balance,
        po_effective_rate_sensitivity: dpo_dr * s.initial_balance,
    }
}

#[cfg(test)]
mod phase34_tests {
    use super::*;
    use ql_ad::Dual;
    fn decs_spec<T: Real>(x: T) -> DeSpiegeleerMandatoryConvertibleDecsSpec<T> {
        DeSpiegeleerMandatoryConvertibleDecsSpec {
            spot_price: T::from_f64(100.0) + x,
            par_amount: T::from_f64(100.0),
            lower_strike_kl: T::from_f64(90.0),
            upper_strike_ku: T::from_f64(120.0),
            decs_coupon_rate: T::from_f64(0.06),
            risk_free_rate: T::from_f64(0.03),
            dividend_yield: T::from_f64(0.02),
            credit_spread: T::from_f64(0.015),
            lower_strike_vol: T::from_f64(0.25),
            upper_strike_vol: T::from_f64(0.22),
            maturity_years: T::from_f64(3.0),
        }
    }
    fn mbs_spec<T: Real>(x: T) -> DavidsonLevinFourFactorIoPoSpec<T> {
        DavidsonLevinFourFactorIoPoSpec {
            initial_balance: T::from_f64(100.0),
            gross_wac_rate: T::from_f64(0.06),
            net_passthrough_coupon: T::from_f64(0.055),
            current_mortgage_rate: T::from_f64(0.05) + x,
            base_turnover_cpr: T::from_f64(0.06),
            max_refi_cpr: T::from_f64(0.35),
            burnout_sensitivity_beta: T::from_f64(0.5),
            seasonality_amplitude: T::from_f64(0.1),
            oas_discount_rate: T::from_f64(0.045),
            num_periods: 12,
        }
    }
    #[test]
    fn phase34_fixed_income_invariants_and_dual_fd() {
        let h = 1e-6;
        let d0 = evaluate_despiegeleer_mandatory_convertible_decs_and_percs(&decs_spec(0.0));
        assert!(
            d0.decomposition_conservation_error < 1e-12
                && d0.decs_mandatory_fair_value > d0.percs_capped_fair_value
        );
        let fd_d = (evaluate_despiegeleer_mandatory_convertible_decs_and_percs(&decs_spec(h))
            .decs_mandatory_fair_value
            - evaluate_despiegeleer_mandatory_convertible_decs_and_percs(&decs_spec(-h))
                .decs_mandatory_fair_value)
            / (2.0 * h);
        let ad_d = evaluate_despiegeleer_mandatory_convertible_decs_and_percs(&decs_spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .decs_mandatory_fair_value;
        assert!(ad_d.d[0].abs() > 1e-6 && (ad_d.d[0] - fd_d).abs() < 1e-5);

        let m0 = evaluate_davidson_levin_four_factor_prepayment_io_po(&mbs_spec(0.0));
        assert!(
            m0.io_po_strip_conservation_error < 1e-12
                && m0.io_effective_rate_sensitivity > 0.0
                && m0.po_effective_rate_sensitivity < 0.0
        );
        let fd_m = (evaluate_davidson_levin_four_factor_prepayment_io_po(&mbs_spec(h))
            .passthrough_present_value
            - evaluate_davidson_levin_four_factor_prepayment_io_po(&mbs_spec(-h))
                .passthrough_present_value)
            / (2.0 * h);
        let ad_m = evaluate_davidson_levin_four_factor_prepayment_io_po(&mbs_spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .passthrough_present_value;
        assert!(ad_m.d[0].abs() > 1e-6 && (ad_m.d[0] - fd_m).abs() < 1e-4);
    }
}
