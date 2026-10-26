use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct DavidsonLevinCmoPacCompanionSpec<T: Real = f64> {
    pub initial_pool_balance: T,
    pub gross_wac_rate: T,
    pub net_coupon_rate: T,
    pub num_periods: usize,
    pub pac_lower_psa: T,
    pub pac_upper_psa: T,
    pub realized_psa: T,
    pub discount_yield: T,
}

#[derive(Clone, Copy, Debug)]
pub struct DavidsonLevinCmoPacCompanionReport<T: Real = f64> {
    pub pac_initial_par_share: T,
    pub support_initial_par_share: T,
    pub pac_present_value: T,
    pub support_present_value: T,
    pub pac_wal_years: T,
    pub support_wal_years: T,
    pub pac_schedule_deviation_l1: T,
    pub support_lockout_broken_collar_amount: T,
    pub principal_waterfall_conservation_error: T,
}

fn mortgage_principal_schedule<T: Real>(s: &DavidsonLevinCmoPacCompanionSpec<T>, psa: T) -> Vec<T> {
    let n = s.num_periods.max(1);
    let r = s.gross_wac_rate / T::from_f64(12.0);
    let one = T::one();
    let payment = if r.abs() < T::from_f64(1e-12) {
        s.initial_pool_balance / T::from_f64(n as f64)
    } else {
        s.initial_pool_balance * r / (one - ((one + r).ln() * T::from_f64(-(n as f64))).exp())
    };
    let mut b = s.initial_pool_balance;
    let mut out = Vec::with_capacity(n);
    for t in 1..=n {
        let interest = b * r;
        let sched = (payment - interest).max_val(T::zero()).min_val(b);
        let cpr = (psa / T::from_f64(100.0))
            * T::from_f64(0.06)
            * T::from_f64(t as f64).min_val(T::from_f64(30.0))
            / T::from_f64(30.0);
        let smm = one - (one - cpr).powf(1.0 / 12.0);
        let prepay = smm * (b - sched);
        let principal = if t == n {
            b
        } else {
            (sched + prepay).min_val(b)
        };
        out.push(principal);
        b -= principal;
    }
    out
}

pub fn evaluate_davidson_levin_cmo_pac_companion_waterfall<T: Real>(
    s: &DavidsonLevinCmoPacCompanionSpec<T>,
) -> DavidsonLevinCmoPacCompanionReport<T> {
    let low = mortgage_principal_schedule(s, s.pac_lower_psa);
    let high = mortgage_principal_schedule(s, s.pac_upper_psa);
    let realized = mortgage_principal_schedule(s, s.realized_psa);
    let eps = T::from_f64(1e-8);
    let pac_sink: Vec<T> = low
        .iter()
        .zip(high.iter())
        .map(|(&a, &b)| (a + b - ((a - b) * (a - b) + eps * eps).sqrt()) / T::from_f64(2.0))
        .collect();
    let pac_par = pac_sink
        .iter()
        .fold(T::zero(), |a, &b| a + b)
        .min_val(s.initial_pool_balance);
    let support_par = (s.initial_pool_balance - pac_par).max_val(T::zero());
    let mut pac_bal = pac_par;
    let mut supp_bal = support_par;
    let mut pac_prin = Vec::with_capacity(realized.len());
    let mut supp_prin = Vec::with_capacity(realized.len());
    let mut broken = T::zero();
    for (t, &coll) in realized.iter().enumerate() {
        let target = pac_sink[t];
        let p = coll.min_val(target).min_val(pac_bal);
        let excess = (coll - p).max_val(T::zero());
        let sup = excess.min_val(supp_bal);
        let overflow = (excess - sup).max_val(T::zero());
        broken += overflow;
        pac_bal -= p;
        supp_bal -= sup;
        pac_prin.push(p + overflow);
        supp_prin.push(sup);
    }
    let mut pac_wal_num = T::zero();
    let mut sup_wal_num = T::zero();
    let mut pac_pv = T::zero();
    let mut sup_pv = T::zero();
    let mut dev = T::zero();
    let mut total = T::zero();
    for t in 0..realized.len() {
        let disc = (-(s.discount_yield / T::from_f64(12.0)) * T::from_f64((t + 1) as f64)).exp();
        pac_wal_num += T::from_f64((t + 1) as f64) / T::from_f64(12.0) * pac_prin[t];
        sup_wal_num += T::from_f64((t + 1) as f64) / T::from_f64(12.0) * supp_prin[t];
        pac_pv += pac_prin[t] * disc;
        sup_pv += supp_prin[t] * disc;
        dev += (pac_prin[t] - pac_sink[t]).abs();
        total += pac_prin[t] + supp_prin[t];
    }
    DavidsonLevinCmoPacCompanionReport {
        pac_initial_par_share: pac_par,
        support_initial_par_share: support_par,
        pac_present_value: pac_pv,
        support_present_value: sup_pv,
        pac_wal_years: pac_wal_num / pac_par.max_val(T::from_f64(1e-12)),
        support_wal_years: sup_wal_num / support_par.max_val(T::from_f64(1e-12)),
        pac_schedule_deviation_l1: dev,
        support_lockout_broken_collar_amount: broken,
        principal_waterfall_conservation_error: (s.initial_pool_balance - total).abs(),
    }
}

#[derive(Clone, Copy, Debug)]
pub struct AsquithSirbuShreveCbGameSpec<T: Real = f64> {
    pub firm_value: T,
    pub stock_price: T,
    pub conversion_ratio: T,
    pub annual_coupon_amount: T,
    pub annual_dividend_per_share: T,
    pub corporate_tax_rate: T,
    pub effective_call_price: T,
    pub call_notice_period_years: T,
    pub firm_volatility: T,
    pub risk_free_rate: T,
    pub firm_payout_rate: T,
    pub dilution_gamma: T,
}

#[derive(Clone, Copy, Debug)]
pub struct AsquithSirbuShreveCbGameReport<T: Real = f64> {
    pub asquith_daily_after_tax_coupon_minus_dividend: T,
    pub asquith_safety_call_trigger_stock_price: T,
    pub conversion_value_to_call_ratio: T,
    pub sirbu_shreve_positive_root_beta: T,
    pub sirbu_shreve_negative_root_beta: T,
    pub endogenous_default_firm_boundary: T,
    pub optimal_conversion_firm_boundary: T,
    pub perpetual_dynkin_cb_fair_value: T,
    pub issuer_call_option_reduction_value: T,
}

pub fn evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game<T: Real>(
    s: &AsquithSirbuShreveCbGameSpec<T>,
) -> AsquithSirbuShreveCbGameReport<T> {
    let daily = ((T::one() - s.corporate_tax_rate) * s.annual_coupon_amount
        - s.conversion_ratio * s.annual_dividend_per_share)
        / T::from_f64(365.0);
    let z = T::from_f64(1.645);
    let trigger = (s.effective_call_price / s.conversion_ratio.max_val(T::from_f64(1e-12)))
        * (s.firm_volatility * s.call_notice_period_years.sqrt() * z).exp();
    let ratio =
        s.conversion_ratio * s.stock_price / s.effective_call_price.max_val(T::from_f64(1e-12));
    let r = s.risk_free_rate.max_val(T::from_f64(1e-8));
    let sig2 = s.firm_volatility * s.firm_volatility;
    let b = s.risk_free_rate - s.firm_payout_rate - T::from_f64(0.5) * sig2;
    let disc = (b * b + T::from_f64(2.0) * sig2 * r).sqrt();
    let beta_p = (-b + disc) / sig2.max_val(T::from_f64(1e-12));
    let beta_m = (-b - disc) / sig2.max_val(T::from_f64(1e-12));
    let vd = beta_m / (beta_m - T::one()) * s.annual_coupon_amount / r;
    let vconv = beta_p / (beta_p - T::one()) * s.annual_coupon_amount
        / (s.dilution_gamma * s.firm_payout_rate).max_val(T::from_f64(1e-12));
    let straight = s.annual_coupon_amount / r
        + (-(s.annual_coupon_amount / r))
            * (s.firm_value / vd.max_val(T::from_f64(1e-12))).pow(beta_m);
    let conversion = s.dilution_gamma * s.firm_value;
    let fair = straight
        .max_val(conversion)
        .min_val(s.effective_call_price.max_val(conversion));
    AsquithSirbuShreveCbGameReport {
        asquith_daily_after_tax_coupon_minus_dividend: daily,
        asquith_safety_call_trigger_stock_price: trigger,
        conversion_value_to_call_ratio: ratio,
        sirbu_shreve_positive_root_beta: beta_p,
        sirbu_shreve_negative_root_beta: beta_m,
        endogenous_default_firm_boundary: vd,
        optimal_conversion_firm_boundary: vconv,
        perpetual_dynkin_cb_fair_value: fair,
        issuer_call_option_reduction_value: (conversion - fair).max_val(T::zero()),
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FasbAsu2024InducedConversionCpdiSpec<T: Real = f64> {
    pub net_income: T,
    pub basic_shares: T,
    pub bond_par: T,
    pub stated_coupon_rate: T,
    pub comparable_noncontingent_yield: T,
    pub tax_rate: T,
    pub original_conversion_shares: T,
    pub induced_sweetener_extra_shares: T,
    pub induced_sweetener_cash: T,
    pub stock_price_at_acceptance: T,
    pub maturity_years: T,
}

#[derive(Clone, Copy, Debug)]
pub struct FasbAsu2024InducedConversionCpdiReport<T: Real = f64> {
    pub basic_eps: T,
    pub asu2020_06_if_converted_diluted_eps: T,
    pub asu2020_06_instrument_c_net_share_eps: T,
    pub asu2024_04_induced_conversion_expense: T,
    pub cpdi_projected_contingent_terminal_payment: T,
    pub cpdi_year1_oid_accrual: T,
    pub cpdi_terminal_positive_adjustment: T,
}

pub fn evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid<T: Real>(
    s: &FasbAsu2024InducedConversionCpdiSpec<T>,
) -> FasbAsu2024InducedConversionCpdiReport<T> {
    let basic_eps = s.net_income / s.basic_shares.max_val(T::from_f64(1e-12));
    let after_tax_interest = s.bond_par * s.stated_coupon_rate * (T::one() - s.tax_rate);
    let ifc = (s.net_income + after_tax_interest)
        / (s.basic_shares + s.original_conversion_shares).max_val(T::from_f64(1e-12));
    let spread = ((s.stock_price_at_acceptance
        - s.bond_par / s.original_conversion_shares.max_val(T::from_f64(1e-12)))
        / s.stock_price_at_acceptance.max_val(T::from_f64(1e-12)))
    .max_val(T::zero())
        * s.original_conversion_shares;
    let inst = (s.net_income / (s.basic_shares + spread).max_val(T::from_f64(1e-12))).max_val(ifc);
    let expense =
        s.induced_sweetener_extra_shares * s.stock_price_at_acceptance + s.induced_sweetener_cash;
    let one_plus = T::one() + s.comparable_noncontingent_yield;
    let growth = one_plus.ln() * s.maturity_years;
    let terminal = s.bond_par * growth.exp()
        - s.stated_coupon_rate * s.bond_par * (growth.exp() - T::one())
            / s.comparable_noncontingent_yield.max_val(T::from_f64(1e-12));
    FasbAsu2024InducedConversionCpdiReport {
        basic_eps,
        asu2020_06_if_converted_diluted_eps: ifc,
        asu2020_06_instrument_c_net_share_eps: inst,
        asu2024_04_induced_conversion_expense: expense,
        cpdi_projected_contingent_terminal_payment: terminal,
        cpdi_year1_oid_accrual: s.bond_par
            * (s.comparable_noncontingent_yield - s.stated_coupon_rate),
        cpdi_terminal_positive_adjustment: terminal.max_val(T::zero()),
    }
}

#[cfg(test)]
mod phase31_tests {
    use super::*;
    use ql_ad::Dual;

    fn cmo<T: Real>(x: T) -> DavidsonLevinCmoPacCompanionSpec<T> {
        DavidsonLevinCmoPacCompanionSpec {
            initial_pool_balance: T::from_f64(100.0) + x,
            gross_wac_rate: T::from_f64(0.05),
            net_coupon_rate: T::from_f64(0.045),
            num_periods: 24,
            pac_lower_psa: T::from_f64(100.0),
            pac_upper_psa: T::from_f64(300.0),
            realized_psa: T::from_f64(180.0),
            discount_yield: T::from_f64(0.04),
        }
    }
    fn cb<T: Real>(x: T) -> AsquithSirbuShreveCbGameSpec<T> {
        AsquithSirbuShreveCbGameSpec {
            firm_value: T::from_f64(130.0),
            stock_price: T::from_f64(100.0) + x,
            conversion_ratio: T::one(),
            annual_coupon_amount: T::from_f64(3.0),
            annual_dividend_per_share: T::one(),
            corporate_tax_rate: T::from_f64(0.25),
            effective_call_price: T::from_f64(105.0),
            call_notice_period_years: T::from_f64(30.0 / 365.0),
            firm_volatility: T::from_f64(0.2),
            risk_free_rate: T::from_f64(0.03),
            firm_payout_rate: T::from_f64(0.02),
            dilution_gamma: T::from_f64(0.8),
        }
    }
    fn asu<T: Real>(x: T) -> FasbAsu2024InducedConversionCpdiSpec<T> {
        FasbAsu2024InducedConversionCpdiSpec {
            net_income: T::from_f64(12.0) + x,
            basic_shares: T::from_f64(100.0),
            bond_par: T::from_f64(100.0),
            stated_coupon_rate: T::from_f64(0.03),
            comparable_noncontingent_yield: T::from_f64(0.04),
            tax_rate: T::from_f64(0.25),
            original_conversion_shares: T::one(),
            induced_sweetener_extra_shares: T::from_f64(0.05),
            induced_sweetener_cash: T::one(),
            stock_price_at_acceptance: T::from_f64(110.0),
            maturity_years: T::from_f64(5.0),
        }
    }

    #[test]
    fn phase31_invariants_and_dual_fd() {
        let c = evaluate_davidson_levin_cmo_pac_companion_waterfall(&cmo(0.0));
        assert!(c.principal_waterfall_conservation_error < 1e-10);
        assert!(c.pac_schedule_deviation_l1 < 1e-3);
        let h = 1e-6;
        let fd = (evaluate_davidson_levin_cmo_pac_companion_waterfall(&cmo(h)).pac_present_value
            - evaluate_davidson_levin_cmo_pac_companion_waterfall(&cmo(-h)).pac_present_value)
            / (2.0 * h);
        let ad =
            evaluate_davidson_levin_cmo_pac_companion_waterfall(&cmo(Dual::<1>::variable(0.0, 0)))
                .pac_present_value;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        let g = evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game(&cb(0.0));
        assert!(g.sirbu_shreve_positive_root_beta > 1.0 && g.sirbu_shreve_negative_root_beta < 0.0);
        assert!(g.optimal_conversion_firm_boundary > g.endogenous_default_firm_boundary);
        let fd_cb = (evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game(&cb(h))
            .conversion_value_to_call_ratio
            - evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game(&cb(-h))
                .conversion_value_to_call_ratio)
            / (2.0 * h);
        let ad_cb = evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game(&cb(
            Dual::<1>::variable(0.0, 0),
        ))
        .conversion_value_to_call_ratio;
        assert!(ad_cb.d[0].abs() > 1e-6 && (ad_cb.d[0] - fd_cb).abs() < 1e-5);
        let a = evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid(&asu(0.0));
        assert!(a.asu2020_06_instrument_c_net_share_eps >= a.asu2020_06_if_converted_diluted_eps);
        assert!(a.asu2024_04_induced_conversion_expense > 0.0);
        let fd_asu = (evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid(&asu(h)).basic_eps
            - evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid(&asu(-h)).basic_eps)
            / (2.0 * h);
        let ad_asu = evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid(&asu(
            Dual::<1>::variable(0.0, 0),
        ))
        .basic_eps;
        assert!(ad_asu.d[0].abs() > 1e-6 && (ad_asu.d[0] - fd_asu).abs() < 1e-5);
    }
}
