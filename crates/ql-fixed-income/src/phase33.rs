use ql_ad::Real;
#[derive(Clone, Copy, Debug)]
pub struct CalamosDeSpiegeleerAscotGammaSpec<T: Real = f64> {
    pub spot_price: T,
    pub reference_spot: T,
    pub bond_par: T,
    pub coupon_rate: T,
    pub conversion_ratio: T,
    pub implied_vol: T,
    pub realized_vol: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub base_hazard_rate: T,
    pub equity_credit_power_p: T,
    pub recovery_rate: T,
    pub maturity_years: T,
    pub holding_period_years: T,
}
#[derive(Clone, Copy, Debug)]
pub struct CalamosDeSpiegeleerAscotGammaReport<T: Real = f64> {
    pub stock_dependent_hazard_rate: T,
    pub risky_credit_bond_floor: T,
    pub conversion_parity_value: T,
    pub convertible_bond_fair_value: T,
    pub ascot_option_strip_value: T,
    pub cb_equity_delta: T,
    pub cb_dollar_gamma: T,
    pub expected_gamma_scalp_pnl: T,
}
pub fn evaluate_calamos_despiegeleer_ascot_and_gamma_scalp<T: Real>(
    s: &CalamosDeSpiegeleerAscotGammaSpec<T>,
) -> CalamosDeSpiegeleerAscotGammaReport<T> {
    let lambda =
        s.base_hazard_rate * (s.reference_spot / s.spot_price).pow(s.equity_credit_power_p);
    let rr = s.risk_free_rate + (T::one() - s.recovery_rate) * lambda;
    let disc = (-rr * s.maturity_years).exp();
    let floor = s.bond_par * disc + s.coupon_rate * s.bond_par * (T::one() - disc) / rr;
    let strike = s.bond_par / s.conversion_ratio;
    let st = s.spot_price;
    let v = s.implied_vol;
    let t = s.maturity_years;
    let d1 = ((st / strike).ln()
        + (s.risk_free_rate - s.dividend_yield + T::from_f64(0.5) * v * v) * t)
        / (v * t.sqrt());
    let d2 = d1 - v * t.sqrt();
    let call = st * (-s.dividend_yield * t).exp() * d1.norm_cdf()
        - strike * (-s.risk_free_rate * t).exp() * d2.norm_cdf();
    let cb = floor + s.conversion_ratio * call;
    let delta = s.conversion_ratio * (-s.dividend_yield * t).exp() * d1.norm_cdf();
    let gamma =
        s.conversion_ratio * (-s.dividend_yield * t).exp() * d1.norm_pdf() / (st * v * t.sqrt());
    CalamosDeSpiegeleerAscotGammaReport {
        stock_dependent_hazard_rate: lambda,
        risky_credit_bond_floor: floor,
        conversion_parity_value: s.conversion_ratio * st,
        convertible_bond_fair_value: cb,
        ascot_option_strip_value: (cb - floor).max_val(T::zero()),
        cb_equity_delta: delta,
        cb_dollar_gamma: gamma * st * st,
        expected_gamma_scalp_pnl: T::from_f64(0.5)
            * gamma
            * st
            * st
            * (s.realized_vol * s.realized_vol - s.implied_vol * s.implied_vol)
            * s.holding_period_years,
    }
}
#[derive(Clone, Copy, Debug)]
pub struct TuckmanSpecifiedPoolDollarRollSpec<T: Real = f64> {
    pub net_coupon_rate: T,
    pub discount_yield: T,
    pub effective_maturity_years: T,
    pub ctd_tba_cpr: T,
    pub llb_pool_cpr_multiplier: T,
    pub tba_front_month_price: T,
    pub tba_back_month_price: T,
    pub gc_repo_rate: T,
}
#[derive(Clone, Copy, Debug)]
pub struct TuckmanSpecifiedPoolDollarRollReport<T: Real = f64> {
    pub ctd_tba_present_value: T,
    pub specified_pool_present_value: T,
    pub specified_pool_payup_points: T,
    pub specified_pool_payup_32nds: T,
    pub dollar_roll_price_drop_32nds: T,
    pub monthly_foregone_carry: T,
    pub dollar_roll_implied_repo_rate: T,
    pub dollar_roll_specialness_bps: T,
}
fn mbs_pv<T: Real>(c: T, y: T, t: T, cpr: T) -> T {
    let lam = -(T::one() - cpr).ln();
    let k = y + lam;
    T::from_f64(100.0) * (c + lam) / k * (T::one() - (-k * t).exp())
        + T::from_f64(100.0) * (-k * t).exp()
}
pub fn evaluate_tuckman_specified_pool_payup_and_dollar_roll<T: Real>(
    s: &TuckmanSpecifiedPoolDollarRollSpec<T>,
) -> TuckmanSpecifiedPoolDollarRollReport<T> {
    let tba = mbs_pv(
        s.net_coupon_rate,
        s.discount_yield,
        s.effective_maturity_years,
        s.ctd_tba_cpr,
    );
    let sp = mbs_pv(
        s.net_coupon_rate,
        s.discount_yield,
        s.effective_maturity_years,
        s.ctd_tba_cpr * s.llb_pool_cpr_multiplier,
    );
    let pay = sp - tba;
    let drop = s.tba_front_month_price - s.tba_back_month_price;
    let carry = s.net_coupon_rate * T::from_f64(100.0) / T::from_f64(12.0);
    let imp = T::from_f64(12.0) / s.tba_front_month_price * (carry - drop);
    TuckmanSpecifiedPoolDollarRollReport {
        ctd_tba_present_value: tba,
        specified_pool_present_value: sp,
        specified_pool_payup_points: pay,
        specified_pool_payup_32nds: pay * T::from_f64(32.0),
        dollar_roll_price_drop_32nds: drop * T::from_f64(32.0),
        monthly_foregone_carry: carry,
        dollar_roll_implied_repo_rate: imp,
        dollar_roll_specialness_bps: (s.gc_repo_rate - imp) * T::from_f64(10000.0),
    }
}
#[cfg(test)]
mod phase33_tests {
    use super::*;
    #[test]
    fn phase33_fixed_income_smoke() {
        let s = TuckmanSpecifiedPoolDollarRollSpec {
            net_coupon_rate: 0.03,
            discount_yield: 0.02,
            effective_maturity_years: 5.0,
            ctd_tba_cpr: 0.08,
            llb_pool_cpr_multiplier: 0.7,
            tba_front_month_price: 100.0,
            tba_back_month_price: 99.8,
            gc_repo_rate: 0.03,
        };
        assert!(
            evaluate_tuckman_specified_pool_payup_and_dollar_roll(&s).specified_pool_payup_32nds
                > 0.0
        );
    }
}
