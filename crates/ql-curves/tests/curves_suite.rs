use ql_ad::Real;
use ql_curves::*;
use ql_time::*;

#[test]
fn discount_short_end_log_linear_anchor() {
    let c = DiscountCurve::new(
        Date::from_ymd(2024, 1, 1).unwrap(),
        vec![1.0, 2.0],
        vec![0.95, 0.90],
    );
    assert!((c.discount_factor(0.5) - 0.95f64.sqrt()).abs() < 1e-12);
}

#[test]
fn multicurve_bootstrap_and_ift_par_dv01_gamma_parity() {
    let ref_date = Date::from_ymd(2024, 1, 1).unwrap();
    let ois_quotes = [0.04, 0.042, 0.045];
    let fwd_quotes = [0.045, 0.047, 0.050];

    let build_system = |oq: &[f64; 3], fq: &[f64; 3]| {
        let d1 = DepositInstrument::new(
            ref_date,
            ref_date.add_days(365),
            oq[0],
            DayCounter::Act365Fixed,
        );
        let o2 = OisSwapInstrument {
            start: ref_date,
            end: ref_date.add_days(365 * 2),
            fixed_rate: oq[1],
            frequency: Frequency::Annual,
            day_counter: DayCounter::Act365Fixed,
        };
        let o3 = OisSwapInstrument {
            start: ref_date,
            end: ref_date.add_days(365 * 3),
            fixed_rate: oq[2],
            frequency: Frequency::Annual,
            day_counter: DayCounter::Act365Fixed,
        };
        let ois_helpers: [&dyn SingleCurveHelper; 3] = [&d1, &o2, &o3];

        let irs = vec![
            IrsSwapInstrument {
                start: ref_date,
                end: ref_date.add_days(365),
                fixed_rate: fq[0],
                frequency: Frequency::Annual,
                day_counter: DayCounter::Act365Fixed,
            },
            IrsSwapInstrument {
                start: ref_date,
                end: ref_date.add_days(365 * 2),
                fixed_rate: fq[1],
                frequency: Frequency::Annual,
                day_counter: DayCounter::Act365Fixed,
            },
            IrsSwapInstrument {
                start: ref_date,
                end: ref_date.add_days(365 * 3),
                fixed_rate: fq[2],
                frequency: Frequency::Annual,
                day_counter: DayCounter::Act365Fixed,
            },
        ];
        let sys =
            CalibratedCurveSystem::from_multicurve_helpers(ref_date, &ois_helpers, &irs).unwrap();
        for h in ois_helpers {
            assert!(h.residual_curve(&sys.curves[0]).abs() < 1e-12);
        }
        for h in &irs {
            assert!(
                h.residual_with_discount(&sys.curves[1], &sys.curves[0])
                    .abs()
                    < 1e-12
            );
        }
        sys
    };

    let sys = build_system(&ois_quotes, &fwd_quotes);
    let times = vec![1.0, 2.0, 3.0];

    let price_swap = |params: &[ql_ad::TapeVar<f64>]| {
        let ois_c = DiscountCurve::new(ref_date, times.clone(), params[0..3].to_vec());
        let fwd_c = DiscountCurve::new(ref_date, times.clone(), params[3..6].to_vec());
        let mut npv = ql_ad::TapeVar::from_f64(0.0);
        for i in 1..=3 {
            let t = i as f64;
            let t_prev = (i - 1) as f64;
            let df = ois_c.discount_factor(t);
            let fwd = fwd_c.discount_factor(t_prev) / fwd_c.discount_factor(t) - 1.0;
            npv += (fwd - 0.048) * df * 1_000_000.0;
        }
        npv
    };

    let (_base_npv, dv01, gamma) = sys.value_par_dv01_and_gamma(price_swap);
    assert_eq!(dv01.len(), 6);
    assert_eq!(gamma.nrows(), 6);
    assert!(gamma.norm() > 1.0);

    let bump = 1e-5;
    for idx in 0..6 {
        let mut oq_p = ois_quotes;
        let mut fq_p = fwd_quotes;
        let mut oq_m = ois_quotes;
        let mut fq_m = fwd_quotes;
        if idx < 3 {
            oq_p[idx] += bump;
            oq_m[idx] -= bump;
        } else {
            fq_p[idx - 3] += bump;
            fq_m[idx - 3] -= bump;
        }
        let sys_p = build_system(&oq_p, &fq_p);
        let sys_m = build_system(&oq_m, &fq_m);
        let (vp, _) = sys_p.value_and_par_dv01(price_swap);
        let (vm, _) = sys_m.value_and_par_dv01(price_swap);
        let fd_dv01 = (vp - vm) / (2.0 * bump);
        assert!(
            (dv01[idx] - fd_dv01).abs() < 1e-2,
            "idx {idx}: ift={} fd={}",
            dv01[idx],
            fd_dv01
        );
    }
}
