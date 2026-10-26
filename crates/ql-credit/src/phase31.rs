use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct RsmAsc470DebtModificationSpec<T: Real = f64> {
    pub carrying_amount: T,
    pub original_par: T,
    pub original_coupon_rate: T,
    pub original_effective_rate: T,
    pub modified_par: T,
    pub modified_coupon_rate: T,
    pub remaining_years: T,
    pub lender_fee_paid: T,
    pub conversion_option_fv_change: T,
    pub market_exit_yield: T,
}

#[derive(Clone, Copy, Debug)]
pub struct RsmAsc470DebtModificationReport<T: Real = f64> {
    pub pv_original_cashflows_at_orig_eir: T,
    pub pv_modified_cashflows_at_orig_eir: T,
    pub asc470_50_pv_difference_pct: T,
    pub smooth_extinguishment_weight: T,
    pub extinguishment_gain_or_loss: T,
    pub prospective_modified_eir: T,
    pub tdr_creditor_concession_loss: T,
}

fn annuity<T: Real>(r: T, n: T) -> T {
    let rr = r.max_val(T::from_f64(1e-10));
    (T::one() - (-((T::one() + rr).ln() * n)).exp()) / rr
}

fn pv_mod<T: Real>(s: &RsmAsc470DebtModificationSpec<T>, r: T) -> T {
    let rr = r.max_val(T::from_f64(1e-10));
    s.modified_coupon_rate * s.modified_par * annuity(rr, s.remaining_years)
        + s.modified_par * (-((T::one() + rr).ln() * s.remaining_years)).exp()
        + s.lender_fee_paid
        + s.conversion_option_fv_change
}

pub fn evaluate_rsm_asc470_debt_modification_and_restructuring<T: Real>(
    s: &RsmAsc470DebtModificationSpec<T>,
) -> RsmAsc470DebtModificationReport<T> {
    let pv_orig = s.original_coupon_rate
        * s.original_par
        * annuity(s.original_effective_rate, s.remaining_years)
        + s.original_par
            * (-((T::one() + s.original_effective_rate).ln() * s.remaining_years)).exp();
    let pv_modified = pv_mod(s, s.original_effective_rate);
    let pct = (pv_modified - pv_orig).abs() / pv_orig.abs().max_val(T::from_f64(1e-12));
    let weight = T::one() / (T::one() + (-(T::from_f64(50.0) * (pct - T::from_f64(0.10)))).exp());
    let mut r = s.market_exit_yield.max_val(T::from_f64(1e-6));
    for _ in 0..8 {
        let pv = pv_mod(s, r);
        let rr = r.max_val(T::from_f64(1e-8));
        let n = s.remaining_years;
        let disc = (-((T::one() + rr).ln() * n)).exp();
        let d_ann = (n * disc * rr / (T::one() + rr) - (T::one() - disc)) / (rr * rr);
        let deriv = s.modified_coupon_rate * s.modified_par * d_ann
            - s.modified_par * n * disc / (T::one() + rr);
        let safe_deriv = if deriv.abs() < T::from_f64(1e-10) {
            T::from_f64(-1e-10)
        } else {
            deriv
        };
        r -= (pv - s.carrying_amount) / safe_deriv;
    }
    let concession = (s.carrying_amount - pv_modified).max_val(T::zero());
    RsmAsc470DebtModificationReport {
        pv_original_cashflows_at_orig_eir: pv_orig,
        pv_modified_cashflows_at_orig_eir: pv_modified,
        asc470_50_pv_difference_pct: pct,
        smooth_extinguishment_weight: weight,
        extinguishment_gain_or_loss: s.carrying_amount - pv_modified,
        prospective_modified_eir: r,
        tdr_creditor_concession_loss: concession,
    }
}

#[cfg(test)]
mod phase31_tests {
    use super::*;
    use ql_ad::Dual;

    fn spec<T: Real>(x: T) -> RsmAsc470DebtModificationSpec<T> {
        RsmAsc470DebtModificationSpec {
            carrying_amount: T::from_f64(100.0) + x,
            original_par: T::from_f64(100.0),
            original_coupon_rate: T::from_f64(0.03),
            original_effective_rate: T::from_f64(0.04),
            modified_par: T::from_f64(100.0),
            modified_coupon_rate: T::from_f64(0.025),
            remaining_years: T::from_f64(5.0),
            lender_fee_paid: T::one(),
            conversion_option_fv_change: T::from_f64(0.5),
            market_exit_yield: T::from_f64(0.05),
        }
    }

    #[test]
    fn phase31_invariants_and_dual_fd() {
        let r = evaluate_rsm_asc470_debt_modification_and_restructuring(&spec(0.0));
        assert!(r.asc470_50_pv_difference_pct >= 0.0);
        assert!(r.prospective_modified_eir.is_finite());
        let h = 1e-6;
        let fd = (evaluate_rsm_asc470_debt_modification_and_restructuring(&spec(h))
            .prospective_modified_eir
            - evaluate_rsm_asc470_debt_modification_and_restructuring(&spec(-h))
                .prospective_modified_eir)
            / (2.0 * h);
        let ad = evaluate_rsm_asc470_debt_modification_and_restructuring(&spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .prospective_modified_eir;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
