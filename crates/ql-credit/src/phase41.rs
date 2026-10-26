//! O'Kane (2008) & Brigo–Mercurio (2006) ISDA Standard CDS Upfront,
//! Dirty RPV01 & CS01/IR01/JTD Triangle Hedge Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct OKaneIsdaStandardCdsUpfrontSpec<T: Real = f64> {
    pub notional_amount: T,
    pub market_par_spread_s: T,
    pub standard_running_coupon_c: T,
    pub recovery_rate_r: T,
    pub ois_zero_rate_r: T,
    pub maturity_years_t: T,
    pub accrued_fraction_alpha: T,
}

#[derive(Debug, Clone, Copy)]
pub struct OKaneIsdaStandardCdsUpfrontReport<T: Real = f64> {
    pub credit_triangle_hazard_lambda: T,
    pub clean_risky_pv01: T,
    pub dirty_risky_pv01: T,
    pub clean_upfront_fraction: T,
    pub cash_settlement_upfront_amount: T,
    pub credit_cs01_per_bp: T,
    pub interest_rate_ir01_per_bp: T,
    pub jump_to_default_jtd_exposure: T,
}

fn clean_rpv01<T: Real>(r: T, lam: T, t: T) -> T {
    let denom = (r + lam).max_val(T::from_f64(1e-6));
    (T::one() - (-(r + lam) * t).exp()) / denom
}

pub fn evaluate_okane_isda_standard_cds_upfront_triangle<T: Real>(
    spec: &OKaneIsdaStandardCdsUpfrontSpec<T>,
) -> OKaneIsdaStandardCdsUpfrontReport<T> {
    let lgd = (T::one() - spec.recovery_rate_r).max_val(T::from_f64(0.05));
    let s = spec.market_par_spread_s.max_val(T::from_f64(1e-5));
    let c = spec.standard_running_coupon_c;
    let r = spec.ois_zero_rate_r;
    let t = spec.maturity_years_t.max_val(T::from_f64(0.25));
    let alpha = spec
        .accrued_fraction_alpha
        .max_val(T::zero())
        .min_val(T::one());

    // Credit triangle hazard rate lambda = S / (1 - R):
    let lam = s / lgd;
    let rpv01_clean = clean_rpv01(r, lam, t);
    let accrual_on_default = T::from_f64(0.5) * T::from_f64(0.25) * lam * rpv01_clean;
    let rpv01_dirty = rpv01_clean + accrual_on_default;

    let clean_upfront_frac = (s - c) * rpv01_clean;
    let accrued_coupon = c * T::from_f64(0.25) * alpha;
    let cash_upfront = spec.notional_amount * (clean_upfront_frac - accrued_coupon);

    // Analytic finite-step CS01 and IR01 sensitivities per 1bp (1e-4):
    let bp = T::from_f64(1e-4);
    let s_up = s + bp;
    let lam_up = s_up / lgd;
    let uf_s_up = (s_up - c) * clean_rpv01(r, lam_up, t);
    let cs01 = spec.notional_amount * (uf_s_up - clean_upfront_frac);

    let uf_r_up = (s - c) * clean_rpv01(r + bp, lam, t);
    let ir01 = spec.notional_amount * (uf_r_up - clean_upfront_frac);

    // Jump-to-default (JTD) net protection-buyer exposure = N * (LGD - CleanUpfrontFrac):
    let jtd = spec.notional_amount * (lgd - clean_upfront_frac);

    OKaneIsdaStandardCdsUpfrontReport {
        credit_triangle_hazard_lambda: lam,
        clean_risky_pv01: rpv01_clean,
        dirty_risky_pv01: rpv01_dirty,
        clean_upfront_fraction: clean_upfront_frac,
        cash_settlement_upfront_amount: cash_upfront,
        credit_cs01_per_bp: cs01,
        interest_rate_ir01_per_bp: ir01,
        jump_to_default_jtd_exposure: jtd,
    }
}

#[cfg(test)]
mod phase41_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn okane_isda_standard_cds_upfront_invariants_and_ad_parity() {
        let spec = OKaneIsdaStandardCdsUpfrontSpec {
            notional_amount: 10_000_000.0,
            market_par_spread_s: 0.0185,
            standard_running_coupon_c: 0.0100,
            recovery_rate_r: 0.40,
            ois_zero_rate_r: 0.038,
            maturity_years_t: 5.0,
            accrued_fraction_alpha: 0.35,
        };
        let rep = evaluate_okane_isda_standard_cds_upfront_triangle(&spec);
        assert!(rep.clean_upfront_fraction > 0.0);
        assert!(rep.credit_cs01_per_bp > 0.0);
        assert!(rep.jump_to_default_jtd_exposure > 0.0);

        let h = 1e-6;
        let r_ad =
            evaluate_okane_isda_standard_cds_upfront_triangle(&OKaneIsdaStandardCdsUpfrontSpec {
                notional_amount: Dual::<1>::from_value(spec.notional_amount),
                market_par_spread_s: Dual::<1>::variable(spec.market_par_spread_s, 0),
                standard_running_coupon_c: Dual::<1>::from_value(spec.standard_running_coupon_c),
                recovery_rate_r: Dual::<1>::from_value(spec.recovery_rate_r),
                ois_zero_rate_r: Dual::<1>::from_value(spec.ois_zero_rate_r),
                maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                accrued_fraction_alpha: Dual::<1>::from_value(spec.accrued_fraction_alpha),
            });
        let mut sp = spec;
        sp.market_par_spread_s += h;
        let mut sm = spec;
        sm.market_par_spread_s -= h;
        let fd = (evaluate_okane_isda_standard_cds_upfront_triangle(&sp)
            .cash_settlement_upfront_amount
            - evaluate_okane_isda_standard_cds_upfront_triangle(&sm)
                .cash_settlement_upfront_amount)
            / (2.0 * h);
        assert!((r_ad.cash_settlement_upfront_amount.d[0] - fd).abs() < 1e-3);
    }
}
