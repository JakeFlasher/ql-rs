use ql_ad::Real;
use ql_math::{norm_cdf, norm_pdf};

#[inline]
fn smin<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b - (d * d + T::from_f64(1e-16)).sqrt())
}

#[inline]
fn smax<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
}

#[derive(Clone, Copy, Debug)]
pub struct FabozziTreasuryFuturesCtdSwitchSpec<T: Real = f64> {
    pub bond_coupons: [T; 2],
    pub bond_maturities: [T; 2],
    pub bond_clean_prices: [T; 2],
    pub bond_modified_durations: [T; 2],
    pub futures_price: T,
    pub delivery_horizon_years: T,
    pub repo_rate: T,
    pub yield_volatility: T,
}

#[derive(Clone, Copy, Debug)]
pub struct FabozziTreasuryFuturesCtdSwitchReport<T: Real = f64> {
    pub conversion_factors: [T; 2],
    pub gross_basis: [T; 2],
    pub net_basis_bnoc: [T; 2],
    pub implied_repo_rates: [T; 2],
    pub ctd_implied_repo_rate: T,
    pub ctd_net_basis_bnoc: T,
    pub quality_switch_option_value: T,
    pub switch_adjusted_theoretical_futures: T,
}

pub fn evaluate_fabozzi_treasury_futures_ctd_and_switch_option<T: Real>(
    spec: &FabozziTreasuryFuturesCtdSwitchSpec<T>,
) -> FabozziTreasuryFuturesCtdSwitchReport<T> {
    let mut cfs = [T::zero(); 2];
    let mut gbasis = [T::zero(); 2];
    let mut nbasis = [T::zero(); 2];
    let mut irrs = [T::zero(); 2];
    let mut f_conv = [T::zero(); 2];
    let tau = spec.delivery_horizon_years.max_val(T::from_f64(1e-6));
    for k in 0..2 {
        let disc = T::from_f64(1.03).pow(-T::from_f64(2.0) * spec.bond_maturities[k]);
        let cf = (spec.bond_coupons[k] / T::from_f64(0.06)) * (T::one() - disc) + disc;
        let p = spec.bond_clean_prices[k];
        let gb = p - spec.futures_price * cf;
        let carry = (T::from_f64(100.0) * spec.bond_coupons[k] - p * spec.repo_rate) * tau;
        let nb = gb - carry;
        let irr = (spec.futures_price * cf + T::from_f64(100.0) * spec.bond_coupons[k] * tau - p)
            / (p * tau);
        cfs[k] = cf;
        gbasis[k] = gb;
        nbasis[k] = nb;
        irrs[k] = irr;
        f_conv[k] = (p - carry) / cf;
    }
    let ctd_irr = smax(irrs[0], irrs[1]);
    let ctd_bnoc = smin(nbasis[0], nbasis[1]);
    let dv0 = spec.bond_modified_durations[0] * spec.bond_clean_prices[0] / cfs[0];
    let dv1 = spec.bond_modified_durations[1] * spec.bond_clean_prices[1] / cfs[1];
    let d_dv = dv1 - dv0;
    let sigma_s =
        ((d_dv * d_dv + T::from_f64(1e-6)) * spec.yield_volatility * spec.yield_volatility * tau)
            .sqrt();
    let diff_f = f_conv[1] - f_conv[0];
    let d = -((diff_f * diff_f + T::from_f64(1e-12)).sqrt()) / sigma_s;
    let switch_opt = (-spec.repo_rate * tau).exp()
        * (sigma_s * norm_pdf(d) + d * sigma_s * norm_cdf(d)).max_val(T::from_f64(1e-8));
    let min_conv = smin(f_conv[0], f_conv[1]);
    FabozziTreasuryFuturesCtdSwitchReport {
        conversion_factors: cfs,
        gross_basis: gbasis,
        net_basis_bnoc: nbasis,
        implied_repo_rates: irrs,
        ctd_implied_repo_rate: ctd_irr,
        ctd_net_basis_bnoc: ctd_bnoc,
        quality_switch_option_value: switch_opt,
        switch_adjusted_theoretical_futures: min_conv - switch_opt,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct AyacheForsythVetzalConvertibleSpec<T: Real = f64> {
    pub spot_price: T,
    pub bond_par: T,
    pub coupon_rate: T,
    pub conversion_ratio: T,
    pub maturity_years: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub equity_volatility: T,
    pub hazard_rate: T,
    pub bond_recovery_rate: T,
    pub equity_default_jump_loss_eta: T,
}

#[derive(Clone, Copy, Debug)]
pub struct AyacheForsythVetzalConvertibleReport<T: Real = f64> {
    pub tf_cash_only_bond_floor: T,
    pub tf_convertible_fair_value: T,
    pub afv_jump_compensated_equity_drift: T,
    pub afv_expected_default_recovery_pv: T,
    pub afv_survival_continuation_pv: T,
    pub afv_convertible_fair_value: T,
    pub tf_vs_afv_splitting_wedge: T,
}

fn bs_call<T: Real>(s: T, k: T, sigma: T, r: T, q: T, t: T) -> T {
    let st = sigma * t.sqrt();
    let d1 = ((s / k).ln() + (r - q + T::from_f64(0.5) * sigma * sigma) * t) / st;
    let d2 = d1 - st;
    s * (-q * t).exp() * norm_cdf(d1) - k * (-r * t).exp() * norm_cdf(d2)
}

pub fn evaluate_ayache_forsyth_vetzal_default_coupled_cb<T: Real>(
    spec: &AyacheForsythVetzalConvertibleSpec<T>,
) -> AyacheForsythVetzalConvertibleReport<T> {
    let sc = (T::one() - spec.bond_recovery_rate) * spec.hazard_rate;
    let rc = (spec.risk_free_rate + sc).max_val(T::from_f64(1e-8));
    let t = spec.maturity_years.max_val(T::from_f64(1e-6));
    let b_tf = spec.bond_par * (-rc * t).exp()
        + spec.coupon_rate * spec.bond_par * (T::one() - (-rc * t).exp()) / rc;
    let k_tf = b_tf / spec.conversion_ratio;
    let v_tf = b_tf
        + spec.conversion_ratio
            * bs_call(
                spec.spot_price,
                k_tf,
                spec.equity_volatility,
                spec.risk_free_rate,
                spec.dividend_yield,
                t,
            );
    let drift_afv = spec.risk_free_rate - spec.dividend_yield
        + spec.hazard_rate * spec.equity_default_jump_loss_eta;
    let q_afv = spec.dividend_yield - spec.hazard_rate * spec.equity_default_jump_loss_eta;
    let r_rf = spec.risk_free_rate.max_val(T::from_f64(1e-8));
    let b_rf = spec.bond_par * (-r_rf * t).exp()
        + spec.coupon_rate * spec.bond_par * (T::one() - (-r_rf * t).exp()) / r_rf;
    let surv = (-spec.hazard_rate * t).exp()
        * (b_rf
            + spec.conversion_ratio
                * bs_call(
                    spec.spot_price,
                    b_rf / spec.conversion_ratio,
                    spec.equity_volatility,
                    spec.risk_free_rate,
                    q_afv,
                    t,
                ));
    let rl = (spec.risk_free_rate + spec.hazard_rate).max_val(T::from_f64(1e-8));
    let def_claim = smax(
        spec.bond_recovery_rate * spec.bond_par,
        spec.conversion_ratio * (T::one() - spec.equity_default_jump_loss_eta) * spec.spot_price,
    );
    let rec_pv = (spec.hazard_rate / rl) * (T::one() - (-rl * t).exp()) * def_claim;
    let v_afv = surv + rec_pv;
    let wedge_diff = v_afv - v_tf;
    AyacheForsythVetzalConvertibleReport {
        tf_cash_only_bond_floor: b_tf,
        tf_convertible_fair_value: v_tf,
        afv_jump_compensated_equity_drift: drift_afv,
        afv_expected_default_recovery_pv: rec_pv,
        afv_survival_continuation_pv: surv,
        afv_convertible_fair_value: v_afv,
        tf_vs_afv_splitting_wedge: (wedge_diff * wedge_diff + T::from_f64(1e-24)).sqrt(),
    }
}

#[cfg(test)]
mod phase35_tests {
    use super::*;
    use ql_ad::Dual;
    fn ctd_spec<T: Real>(x: T) -> FabozziTreasuryFuturesCtdSwitchSpec<T> {
        FabozziTreasuryFuturesCtdSwitchSpec {
            bond_coupons: [T::from_f64(0.06), T::from_f64(0.045)],
            bond_maturities: [T::from_f64(7.0), T::from_f64(9.5)],
            bond_clean_prices: [T::from_f64(101.0) + x, T::from_f64(94.5)],
            bond_modified_durations: [T::from_f64(5.8), T::from_f64(7.6)],
            futures_price: T::from_f64(100.2),
            delivery_horizon_years: T::from_f64(0.25),
            repo_rate: T::from_f64(0.04),
            yield_volatility: T::from_f64(0.01),
        }
    }
    fn afv_spec<T: Real>(x: T) -> AyacheForsythVetzalConvertibleSpec<T> {
        AyacheForsythVetzalConvertibleSpec {
            spot_price: T::from_f64(100.0) + x,
            bond_par: T::from_f64(100.0),
            coupon_rate: T::from_f64(0.04),
            conversion_ratio: T::one(),
            maturity_years: T::from_f64(3.0),
            risk_free_rate: T::from_f64(0.03),
            dividend_yield: T::from_f64(0.015),
            equity_volatility: T::from_f64(0.25),
            hazard_rate: T::from_f64(0.025),
            bond_recovery_rate: T::from_f64(0.4),
            equity_default_jump_loss_eta: T::from_f64(0.7),
        }
    }
    #[test]
    fn phase35_fixed_income_invariants_and_dual_fd() {
        let c0 = evaluate_fabozzi_treasury_futures_ctd_and_switch_option(&ctd_spec(0.0));
        assert!((c0.conversion_factors[0] - 1.0).abs() < 1e-12);
        assert!(c0.quality_switch_option_value > 0.0);
        let h = 1e-6;
        let fd_c = (evaluate_fabozzi_treasury_futures_ctd_and_switch_option(&ctd_spec(h))
            .switch_adjusted_theoretical_futures
            - evaluate_fabozzi_treasury_futures_ctd_and_switch_option(&ctd_spec(-h))
                .switch_adjusted_theoretical_futures)
            / (2.0 * h);
        let ad_c = evaluate_fabozzi_treasury_futures_ctd_and_switch_option(&ctd_spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .switch_adjusted_theoretical_futures;
        assert!(ad_c.d[0].abs() > 1e-6 && (ad_c.d[0] - fd_c).abs() < 1e-5);

        let a0 = evaluate_ayache_forsyth_vetzal_default_coupled_cb(&afv_spec(0.0));
        assert!(a0.afv_jump_compensated_equity_drift > 0.03 - 0.015);
        assert!(a0.afv_convertible_fair_value >= a0.tf_cash_only_bond_floor);
        let fd_a = (evaluate_ayache_forsyth_vetzal_default_coupled_cb(&afv_spec(h))
            .afv_convertible_fair_value
            - evaluate_ayache_forsyth_vetzal_default_coupled_cb(&afv_spec(-h))
                .afv_convertible_fair_value)
            / (2.0 * h);
        let ad_a = evaluate_ayache_forsyth_vetzal_default_coupled_cb(&afv_spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .afv_convertible_fair_value;
        assert!(ad_a.d[0].abs() > 1e-6 && (ad_a.d[0] - fd_a).abs() < 1e-5);
    }
}
