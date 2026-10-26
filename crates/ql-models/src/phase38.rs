#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Clone, Copy, Debug)]
pub struct BroadieGlassermanKouDiscreteBarrierSpec<T: Real = f64> {
    pub spot_price: T,
    pub strike_price: T,
    pub down_barrier_h: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub volatility: T,
    pub maturity_years: T,
    pub monitoring_steps_m: T,
}

#[derive(Clone, Copy, Debug)]
pub struct BroadieGlassermanKouDiscreteBarrierReport<T: Real = f64> {
    pub shifted_effective_barrier: T,
    pub continuous_down_and_out_call: T,
    pub discrete_bgk_down_and_out_call: T,
    pub continuity_correction_price_uplift: T,
    pub continuous_lookback_call: T,
    pub discrete_bgk_lookback_call: T,
}

fn bs_call_generic<T: Real>(s: T, k: T, r: T, q: T, sig: T, mat: T) -> T {
    let st = (sig * mat.sqrt()).max_val(T::from_f64(1e-6));
    let d1 = ((s / k).ln() + (r - q + T::from_f64(0.5) * sig * sig) * mat) / st;
    let d2 = d1 - st;
    s * (-q * mat).exp() * norm_cdf(d1) - k * (-r * mat).exp() * norm_cdf(d2)
}

fn down_and_out_call<T: Real>(s: T, k: T, b: T, r: T, q: T, sig: T, mat: T) -> T {
    let c_bs = bs_call_generic(s, k, r, q, sig, mat);
    let sig2 = (sig * sig).max_val(T::from_f64(1e-6));
    let lam = (r - q - T::from_f64(0.5) * sig2) / sig2;
    let refl_s = (b * b) / s;
    let factor = ((b / s).ln() * T::from_f64(2.0) * lam).exp();
    let c_refl = bs_call_generic(refl_s, k, r, q, sig, mat);
    (c_bs - factor * c_refl).max_val(T::from_f64(1e-6))
}

pub fn evaluate_broadie_glasserman_kou_discrete_barrier_lookback<T: Real>(
    spec: &BroadieGlassermanKouDiscreteBarrierSpec<T>,
) -> BroadieGlassermanKouDiscreteBarrierReport<T> {
    let beta1 = T::from_f64(0.582_597_157_939_010_7);
    let s = spec.spot_price.max_val(T::from_f64(1e-4));
    let k = spec.strike_price.max_val(T::from_f64(1e-4));
    let h_bar = spec.down_barrier_h.max_val(T::from_f64(1e-4));
    let sig = spec.volatility.max_val(T::from_f64(0.01));
    let mat = spec.maturity_years.max_val(T::from_f64(1e-4));
    let m_steps = spec.monitoring_steps_m.max_val(T::from_f64(1.0));
    let dt = mat / m_steps;
    let shift = (-beta1 * sig * dt.sqrt()).exp();
    let h_eff = h_bar * shift;

    let c_cont = down_and_out_call(
        s,
        k,
        h_bar,
        spec.risk_free_rate,
        spec.dividend_yield,
        sig,
        mat,
    );
    let c_disc = down_and_out_call(
        s,
        k,
        h_eff,
        spec.risk_free_rate,
        spec.dividend_yield,
        sig,
        mat,
    );

    let atm_call = bs_call_generic(s, s, spec.risk_free_rate, spec.dividend_yield, sig, mat);
    let lb_cont = atm_call + T::from_f64(0.5) * s * sig * mat.sqrt();
    let lb_disc = shift * lb_cont + (T::one() - shift) * atm_call;

    BroadieGlassermanKouDiscreteBarrierReport {
        shifted_effective_barrier: h_eff,
        continuous_down_and_out_call: c_cont,
        discrete_bgk_down_and_out_call: c_disc,
        continuity_correction_price_uplift: c_disc - c_cont,
        continuous_lookback_call: lb_cont,
        discrete_bgk_lookback_call: lb_disc,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct OhlsonPenmanLinearInfoDynamicsSpec<T: Real = f64> {
    pub current_book_value_bt: T,
    pub current_earnings_xt: T,
    pub current_dividends_dt: T,
    pub other_information_shock_nu: T,
    pub cost_of_equity_re: T,
    pub abnormal_earnings_persistence_omega: T,
    pub other_info_persistence_gamma: T,
}

#[derive(Clone, Copy, Debug)]
pub struct OhlsonPenmanLinearInfoDynamicsReport<T: Real = f64> {
    pub prior_clean_surplus_book_value: T,
    pub current_abnormal_earnings: T,
    pub abnormal_earnings_multiplier_alpha1: T,
    pub other_information_multiplier_alpha2: T,
    pub earnings_valuation_weight_k: T,
    pub lid_intrinsic_equity_value: T,
    pub weighted_earnings_book_value: T,
    pub ohlson_weighted_book_earnings_error: T,
    pub implied_price_to_book_ratio: T,
}

pub fn evaluate_ohlson_penman_linear_information_dynamics<T: Real>(
    spec: &OhlsonPenmanLinearInfoDynamicsSpec<T>,
) -> OhlsonPenmanLinearInfoDynamicsReport<T> {
    let bt = spec.current_book_value_bt.max_val(T::from_f64(1e-4));
    let xt = spec.current_earnings_xt;
    let dt = spec.current_dividends_dt;
    let nu = spec.other_information_shock_nu;
    let re = spec.cost_of_equity_re.max_val(T::from_f64(1e-4));
    let big_r = T::one() + re;
    let omega = spec.abnormal_earnings_persistence_omega;
    let gamma = spec.other_info_persistence_gamma;

    let b_tm1 = bt + dt - xt;
    let x_a = xt - re * b_tm1;

    let r_m_om = (big_r - omega).max_val(T::from_f64(1e-6));
    let r_m_gm = (big_r - gamma).max_val(T::from_f64(1e-6));
    let alpha1 = omega / r_m_om;
    let alpha2 = big_r / (r_m_om * r_m_gm);
    let k_wt = re * alpha1;
    let phi = big_r / re;

    let v_lid = bt + alpha1 * x_a + alpha2 * nu;
    let v_weighted = k_wt * (phi * xt - dt) + (T::one() - k_wt) * bt + alpha2 * nu;
    let d_eq = v_lid - v_weighted;
    let eq_err = (d_eq * d_eq + T::from_f64(1e-30)).sqrt();

    OhlsonPenmanLinearInfoDynamicsReport {
        prior_clean_surplus_book_value: b_tm1,
        current_abnormal_earnings: x_a,
        abnormal_earnings_multiplier_alpha1: alpha1,
        other_information_multiplier_alpha2: alpha2,
        earnings_valuation_weight_k: k_wt,
        lid_intrinsic_equity_value: v_lid,
        weighted_earnings_book_value: v_weighted,
        ohlson_weighted_book_earnings_error: eq_err,
        implied_price_to_book_ratio: v_lid / bt,
    }
}

#[cfg(test)]
mod phase38_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn bgk_barrier_and_ohlson_lid_invariants_and_ad_parity() {
        let h = 1e-6;

        let bgk_spec = BroadieGlassermanKouDiscreteBarrierSpec {
            spot_price: 100.0,
            strike_price: 100.0,
            down_barrier_h: 90.0,
            risk_free_rate: 0.04,
            dividend_yield: 0.01,
            volatility: 0.22,
            maturity_years: 1.0,
            monitoring_steps_m: 50.0,
        };
        let bgk_rep = evaluate_broadie_glasserman_kou_discrete_barrier_lookback(&bgk_spec);
        assert!(bgk_rep.shifted_effective_barrier < bgk_spec.down_barrier_h);
        assert!(bgk_rep.discrete_bgk_down_and_out_call > bgk_rep.continuous_down_and_out_call);
        let bgk_ad = evaluate_broadie_glasserman_kou_discrete_barrier_lookback(
            &BroadieGlassermanKouDiscreteBarrierSpec {
                spot_price: Dual::<1>::variable(bgk_spec.spot_price, 0),
                strike_price: Dual::<1>::from_value(bgk_spec.strike_price),
                down_barrier_h: Dual::<1>::from_value(bgk_spec.down_barrier_h),
                risk_free_rate: Dual::<1>::from_value(bgk_spec.risk_free_rate),
                dividend_yield: Dual::<1>::from_value(bgk_spec.dividend_yield),
                volatility: Dual::<1>::from_value(bgk_spec.volatility),
                maturity_years: Dual::<1>::from_value(bgk_spec.maturity_years),
                monitoring_steps_m: Dual::<1>::from_value(bgk_spec.monitoring_steps_m),
            },
        );
        let mut bp = bgk_spec;
        bp.spot_price += h;
        let mut bm = bgk_spec;
        bm.spot_price -= h;
        let fd_bgk = (evaluate_broadie_glasserman_kou_discrete_barrier_lookback(&bp)
            .discrete_bgk_down_and_out_call
            - evaluate_broadie_glasserman_kou_discrete_barrier_lookback(&bm)
                .discrete_bgk_down_and_out_call)
            / (2.0 * h);
        assert!((bgk_ad.discrete_bgk_down_and_out_call.d[0] - fd_bgk).abs() < 1e-5);

        let ohl_spec = OhlsonPenmanLinearInfoDynamicsSpec {
            current_book_value_bt: 100.0,
            current_earnings_xt: 14.5,
            current_dividends_dt: 4.5,
            other_information_shock_nu: 1.8,
            cost_of_equity_re: 0.09,
            abnormal_earnings_persistence_omega: 0.62,
            other_info_persistence_gamma: 0.45,
        };
        let ohl_rep = evaluate_ohlson_penman_linear_information_dynamics(&ohl_spec);
        assert!(ohl_rep.ohlson_weighted_book_earnings_error < 1e-10);
        assert!(ohl_rep.abnormal_earnings_multiplier_alpha1 > 0.0);
        let ohl_ad = evaluate_ohlson_penman_linear_information_dynamics(
            &OhlsonPenmanLinearInfoDynamicsSpec {
                current_book_value_bt: Dual::<1>::from_value(ohl_spec.current_book_value_bt),
                current_earnings_xt: Dual::<1>::variable(ohl_spec.current_earnings_xt, 0),
                current_dividends_dt: Dual::<1>::from_value(ohl_spec.current_dividends_dt),
                other_information_shock_nu: Dual::<1>::from_value(
                    ohl_spec.other_information_shock_nu,
                ),
                cost_of_equity_re: Dual::<1>::from_value(ohl_spec.cost_of_equity_re),
                abnormal_earnings_persistence_omega: Dual::<1>::from_value(
                    ohl_spec.abnormal_earnings_persistence_omega,
                ),
                other_info_persistence_gamma: Dual::<1>::from_value(
                    ohl_spec.other_info_persistence_gamma,
                ),
            },
        );
        let mut op = ohl_spec;
        op.current_earnings_xt += h;
        let mut om = ohl_spec;
        om.current_earnings_xt -= h;
        let fd_ohl = (evaluate_ohlson_penman_linear_information_dynamics(&op)
            .lid_intrinsic_equity_value
            - evaluate_ohlson_penman_linear_information_dynamics(&om).lid_intrinsic_equity_value)
            / (2.0 * h);
        assert!((ohl_ad.lid_intrinsic_equity_value.d[0] - fd_ohl).abs() < 1e-5);
    }
}
