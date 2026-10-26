#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
use ql_ad::Real;
use ql_math::{norm_cdf, norm_pdf};

#[inline]
fn smax<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
}

#[derive(Clone, Copy, Debug)]
pub struct BrigoGregoryMultiCurrencyCsaCtdSpec<T: Real = f64> {
    pub domestic_ois_rate: T,
    pub foreign_ois_rate: T,
    pub xccy_basis_spread: T,
    pub csa_spread_volatility: T,
    pub spot_fx: T,
    pub maturity_years: T,
    pub promised_cashflow: T,
}

#[derive(Clone, Copy, Debug)]
pub struct BrigoGregoryMultiCurrencyCsaCtdReport<T: Real = f64> {
    pub domestic_ois_df: T,
    pub foreign_csa_equiv_rate: T,
    pub foreign_csa_equiv_df: T,
    pub basis_adjusted_fx_forward: T,
    pub intrinsic_ctd_collateral_rate: T,
    pub expected_ctd_collateral_rate: T,
    pub csa_switch_option_rate_pickup_bps: T,
    pub ois_discounted_pv: T,
    pub ctd_csa_discounted_pv: T,
    pub ctd_valuation_discount: T,
}

pub fn evaluate_brigo_gregory_multicurrency_csa_ctd<T: Real>(
    spec: &BrigoGregoryMultiCurrencyCsaCtdSpec<T>,
) -> BrigoGregoryMultiCurrencyCsaCtdReport<T> {
    let mat = spec.maturity_years.max_val(T::from_f64(1e-8));
    let rd = spec.domestic_ois_rate;
    let rf = spec.foreign_ois_rate;
    let bx = spec.xccy_basis_spread;
    let r_fx_csa = rd - bx;
    let df_ois = (-rd * mat).exp();
    let df_fx_csa = (-r_fx_csa * mat).exp();
    let fx_fwd = spec.spot_fx * ((rd - rf - bx) * mat).exp();
    let intrinsic_ctd = smax(rd, r_fx_csa);
    let sig_t = (spec.csa_spread_volatility * mat.sqrt()).max_val(T::from_f64(1e-8));
    let d = (r_fx_csa - rd) / sig_t;
    let expected_ctd = rd + (r_fx_csa - rd) * norm_cdf(d) + sig_t * norm_pdf(d);
    let pickup_bps = smax(expected_ctd - intrinsic_ctd, T::from_f64(1e-8)) * T::from_f64(1e4);
    let ois_pv = spec.promised_cashflow * df_ois;
    let ctd_pv = spec.promised_cashflow * (-expected_ctd * mat).exp();
    let discount = ois_pv - ctd_pv;
    BrigoGregoryMultiCurrencyCsaCtdReport {
        domestic_ois_df: df_ois,
        foreign_csa_equiv_rate: r_fx_csa,
        foreign_csa_equiv_df: df_fx_csa,
        basis_adjusted_fx_forward: fx_fwd,
        intrinsic_ctd_collateral_rate: intrinsic_ctd,
        expected_ctd_collateral_rate: expected_ctd,
        csa_switch_option_rate_pickup_bps: pickup_bps,
        ois_discounted_pv: ois_pv,
        ctd_csa_discounted_pv: ctd_pv,
        ctd_valuation_discount: discount,
    }
}

#[cfg(test)]
mod phase36_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn brigo_gregory_multicurrency_csa_ctd_invariants_and_ad_parity() {
        let spec = BrigoGregoryMultiCurrencyCsaCtdSpec {
            domestic_ois_rate: 0.042,
            foreign_ois_rate: 0.028,
            xccy_basis_spread: -0.0025,
            csa_spread_volatility: 0.0065,
            spot_fx: 1.085,
            maturity_years: 5.0,
            promised_cashflow: 100.0,
        };
        let rep = evaluate_brigo_gregory_multicurrency_csa_ctd(&spec);
        assert!(rep.csa_switch_option_rate_pickup_bps > 0.0);
        assert!(rep.ctd_csa_discounted_pv < rep.ois_discounted_pv);

        let h = 1e-6;
        let r_ad =
            evaluate_brigo_gregory_multicurrency_csa_ctd(&BrigoGregoryMultiCurrencyCsaCtdSpec {
                domestic_ois_rate: Dual::<1>::variable(spec.domestic_ois_rate, 0),
                foreign_ois_rate: Dual::<1>::from_value(spec.foreign_ois_rate),
                xccy_basis_spread: Dual::<1>::from_value(spec.xccy_basis_spread),
                csa_spread_volatility: Dual::<1>::from_value(spec.csa_spread_volatility),
                spot_fx: Dual::<1>::from_value(spec.spot_fx),
                maturity_years: Dual::<1>::from_value(spec.maturity_years),
                promised_cashflow: Dual::<1>::from_value(spec.promised_cashflow),
            });
        let mut sp = spec;
        sp.domestic_ois_rate += h;
        let mut sm = spec;
        sm.domestic_ois_rate -= h;
        let _fd = (evaluate_brigo_gregory_multicurrency_csa_ctd(&sp).ctd_csa_discounted_pv
            - evaluate_brigo_gregory_multicurrency_csa_ctd(&sm).ctd_csa_discounted_pv)
            / (2.0 * h);
        assert!((r_ad.ctd_csa_discounted_pv.d[0] - _fd).abs() < 1e-5);
    }
}
