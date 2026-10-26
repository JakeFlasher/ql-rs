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
    let mut io = T::zero();
    let mut po = T::zero();
    for i in 1..=s.num_periods {
        let age = T::from_f64(i as f64 / 30.).min_val(T::one());
        let burnout = (bal.max_val(T::from_f64(1e-12)).ln() * s.burnout_sensitivity_beta).exp();
        let refi = s.max_refi_cpr
            * burnout
            * (T::one()
                + (T::from_f64(-40.0)
                    * (s.gross_wac_rate - s.current_mortgage_rate - T::from_f64(0.005)))
                .exp())
            .powf(-1.0);
        let cpr_annual = (s.base_turnover_cpr * age + refi)
            * (T::one()
                + s.seasonality_amplitude
                    * (T::from_f64(i as f64) * T::from_f64(std::f64::consts::PI / 6.).sin()))
            .min_val(T::from_f64(0.95));
        let smm = T::one() - (T::one() - cpr_annual).powf(1.0 / 12.0);
        let pre = bal * smm;
        let disc = (-s.oas_discount_rate * T::from_f64(i as f64 / 12.)).exp();
        io += bal * s.net_passthrough_coupon / T::from_f64(12.) * disc;
        po += pre * disc;
        bal *= T::one() - smm;
    }
    let pv = io + po;
    let avg = T::one()
        - (bal.max_val(T::from_f64(1e-12)).ln() * (T::one() / T::from_f64(s.num_periods as f64)))
            .exp();
    DavidsonLevinFourFactorIoPoReport {
        terminal_pool_factor: bal,
        average_realized_cpr: avg,
        passthrough_present_value: pv * s.initial_balance,
        io_strip_present_value: io * s.initial_balance,
        po_strip_present_value: po * s.initial_balance,
        io_po_strip_conservation_error: (pv - io - po).abs(),
        io_effective_rate_sensitivity: io * T::from_f64(1e-4),
        po_effective_rate_sensitivity: -po * T::from_f64(1e-4),
    }
}
