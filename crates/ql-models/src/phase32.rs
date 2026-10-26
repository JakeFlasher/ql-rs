use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct GairatPolalaCheyetteCapletSpec<T: Real = f64> {
    pub forward_rate: T,
    pub strike_rate: T,
    pub fixing_time_t1: T,
    pub payment_time_t2: T,
    pub mean_reversion_kappa: T,
    pub initial_sigma0: T,
    pub initial_skew_b: T,
    pub target_atm_normal_vol: T,
    pub target_otm_normal_vol: T,
    pub discount_factor_t2: T,
}
#[derive(Clone, Copy, Debug)]
pub struct GairatPolalaCheyetteCapletReport<T: Real = f64> {
    pub deterministic_auxiliary_y_bar: T,
    pub bond_vol_factor_b: T,
    pub pre_calibration_normal_vol: T,
    pub calibrated_sigma0: T,
    pub calibrated_skew_b: T,
    pub calibrated_normal_vol: T,
    pub bachelier_caplet_pv: T,
    pub calibration_rmse: T,
}
pub fn evaluate_gairat_polala_cheyette_local_vol_caplet<T: Real>(
    s: &GairatPolalaCheyetteCapletSpec<T>,
) -> GairatPolalaCheyetteCapletReport<T> {
    let k = s.mean_reversion_kappa;
    let t = s.fixing_time_t1.max_val(T::from_f64(1e-8));
    let tau = (s.payment_time_t2 - s.fixing_time_t1).max_val(T::from_f64(1e-8));
    let e = (-T::from_f64(2.0) * k * t).exp();
    let q = (T::one() - e) / (T::from_f64(2.0) * k).max_val(T::from_f64(1e-8));
    let y = s.initial_sigma0
        * s.initial_sigma0
        * q
        * (T::one()
            + T::from_f64(1.5)
                * s.initial_skew_b
                * s.initial_skew_b
                * s.initial_sigma0
                * s.initial_sigma0
                * q);
    let b = (T::one() - (-k * tau).exp()) / k.max_val(T::from_f64(1e-8));
    let base = (T::one() + tau * s.forward_rate) / tau * b * (y / t).sqrt();
    let skew = T::one()
        + T::from_f64(0.5) * s.initial_skew_b * b * (s.strike_rate - s.forward_rate)
            / s.forward_rate.max_val(T::from_f64(1e-8));
    let pre = base * skew;
    let f_safe = s.forward_rate.max_val(T::from_f64(1e-8));
    let dk_rel = (s.strike_rate - s.forward_rate) / f_safe;
    let dk_safe = if dk_rel.abs() < T::from_f64(1e-6) {
        T::from_f64(1e-6)
    } else {
        dk_rel
    };
    let atm_target = s.target_atm_normal_vol.max_val(T::from_f64(1e-8));
    let cal_skew = T::from_f64(2.0) * (s.target_otm_normal_vol / atm_target - T::one())
        / (b * dk_safe).max_val(T::from_f64(1e-12));
    let scale_factor = (T::one() + tau * s.forward_rate) / tau * b / t.sqrt();
    let y_target = (atm_target / scale_factor.max_val(T::from_f64(1e-8))).powf(2.0);
    let disc_y = (T::one() + T::from_f64(6.0) * cal_skew * cal_skew * y_target).sqrt();
    let u_cal = T::from_f64(2.0) * y_target / (T::one() + disc_y);
    let cal_sigma = (u_cal / q.max_val(T::from_f64(1e-12)))
        .max_val(T::from_f64(1e-12))
        .sqrt();
    let y_cal = u_cal * (T::one() + T::from_f64(1.5) * cal_skew * cal_skew * u_cal);
    let atm_cal = scale_factor * y_cal.max_val(T::from_f64(1e-16)).sqrt();
    let cal_vol = atm_cal * (T::one() + T::from_f64(0.5) * cal_skew * b * dk_rel);
    let d = (s.forward_rate - s.strike_rate) / (cal_vol * tau.sqrt()).max_val(T::from_f64(1e-8));
    let pv = s.discount_factor_t2
        * tau
        * ((s.forward_rate - s.strike_rate) * d.norm_cdf() + cal_vol * tau.sqrt() * d.norm_pdf());
    let rmse = (((cal_vol - s.target_otm_normal_vol) * (cal_vol - s.target_otm_normal_vol)
        + (atm_cal - s.target_atm_normal_vol) * (atm_cal - s.target_atm_normal_vol))
        / T::from_f64(2.0)
        + T::from_f64(1e-24))
    .sqrt();
    GairatPolalaCheyetteCapletReport {
        deterministic_auxiliary_y_bar: y,
        bond_vol_factor_b: b,
        pre_calibration_normal_vol: pre,
        calibrated_sigma0: cal_sigma,
        calibrated_skew_b: cal_skew,
        calibrated_normal_vol: cal_vol,
        bachelier_caplet_pv: pv,
        calibration_rmse: rmse,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct GatheralJusselinQuadraticRoughHestonSpec<T: Real = f64> {
    pub forward: T,
    pub strike: T,
    pub expiry: T,
    pub hurst_h: T,
    pub weight_a: T,
    pub asymmetry_center_b: T,
    pub floor_variance_c: T,
    pub initial_state_z0: T,
    pub kernel_eta: T,
    pub mean_reversion_lambda: T,
}
#[derive(Clone, Copy, Debug)]
pub struct GatheralJusselinQuadraticRoughHestonReport<T: Real = f64> {
    pub instantaneous_variance_v0: T,
    pub variance_feedback_slope: T,
    pub zumbach_time_reversal_asymmetry: T,
    pub qrh_atm_implied_vol: T,
    pub qrh_atm_skew: T,
    pub qrh_atm_curvature: T,
    pub qrh_strike_implied_vol: T,
    pub qrh_call_price: T,
}
pub fn evaluate_gatheral_jusselin_quadratic_rough_heston<T: Real>(
    s: &GatheralJusselinQuadraticRoughHestonSpec<T>,
) -> GatheralJusselinQuadraticRoughHestonReport<T> {
    let v = s.weight_a
        * (s.initial_state_z0 - s.asymmetry_center_b)
        * (s.initial_state_z0 - s.asymmetry_center_b)
        + s.floor_variance_c;
    let slope = T::from_f64(2.0) * s.weight_a * (s.initial_state_z0 - s.asymmetry_center_b);
    let h_plus = s.hurst_h + T::from_f64(0.5);
    let t_h = s.expiry.max_val(T::from_f64(1e-8)).pow(h_plus);
    let mr_damp = T::one()
        / (T::one()
            + s.mean_reversion_lambda.max_val(T::zero()) * t_h
            + T::from_f64(0.5) * s.mean_reversion_lambda * s.mean_reversion_lambda * t_h * t_h);
    let gamma_h32 = T::from_f64(0.886226925452758)
        * (T::one() + T::from_f64(0.45) * s.hurst_h + T::from_f64(0.15) * s.hurst_h * s.hurst_h);
    let asym = s.weight_a * s.asymmetry_center_b * s.kernel_eta * t_h * mr_damp
        / (gamma_h32 * v.max_val(T::from_f64(1e-8)));
    let v_eff = s.floor_variance_c + (v - s.floor_variance_c) * mr_damp;
    let atm = v_eff.max_val(T::from_f64(1e-12)).sqrt();
    let skew = slope * s.kernel_eta * mr_damp / (T::from_f64(2.0) * atm.max_val(T::from_f64(1e-8)));
    let curv = T::from_f64(2.0) * s.weight_a * s.kernel_eta * s.kernel_eta * mr_damp;
    let k = (s.strike / s.forward.max_val(T::from_f64(1e-8))).ln();
    let iv =
        (atm * (T::one() + skew * k + T::from_f64(0.5) * curv * k * k)).max_val(T::from_f64(1e-8));
    let d1 = ((s.forward / s.strike.max_val(T::from_f64(1e-8))).ln()
        + T::from_f64(0.5) * iv * iv * s.expiry)
        / (iv * s.expiry.sqrt()).max_val(T::from_f64(1e-8));
    let d2 = d1 - iv * s.expiry.sqrt();
    let call = s.forward * d1.norm_cdf() - s.strike * d2.norm_cdf();
    GatheralJusselinQuadraticRoughHestonReport {
        instantaneous_variance_v0: v,
        variance_feedback_slope: slope,
        zumbach_time_reversal_asymmetry: asym,
        qrh_atm_implied_vol: atm,
        qrh_atm_skew: skew,
        qrh_atm_curvature: curv,
        qrh_strike_implied_vol: iv,
        qrh_call_price: call,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SchilitSloanDechowDichevQualitySpec<T: Real = f64> {
    pub revenue_0: T,
    pub revenue_1: T,
    pub cogs_0: T,
    pub cogs_1: T,
    pub accounts_receivable_0: T,
    pub accounts_receivable_1: T,
    pub inventory_0: T,
    pub inventory_1: T,
    pub accounts_payable_0: T,
    pub accounts_payable_1: T,
    pub net_operating_assets_0: T,
    pub net_operating_assets_1: T,
    pub net_income_1: T,
    pub cfo_prev: T,
    pub cfo_curr: T,
    pub cfo_next: T,
}
#[derive(Clone, Copy, Debug)]
pub struct SchilitSloanDechowDichevQualityReport<T: Real = f64> {
    pub days_sales_outstanding_1: T,
    pub days_inventory_on_hand_1: T,
    pub days_payable_outstanding_1: T,
    pub cash_conversion_cycle_0: T,
    pub cash_conversion_cycle_1: T,
    pub ccc_deterioration_days: T,
    pub sloan_balance_sheet_accrual_ratio: T,
    pub cash_flow_accrual_ratio: T,
    pub dechow_dichev_accrual_residual: T,
    pub schilit_shenanigan_composite_z: T,
}
pub fn evaluate_schilit_sloan_dechow_dichev_earnings_quality<T: Real>(
    s: &SchilitSloanDechowDichevQualitySpec<T>,
) -> SchilitSloanDechowDichevQualityReport<T> {
    let dso1 =
        T::from_f64(365.0) * s.accounts_receivable_1 / s.revenue_1.max_val(T::from_f64(1e-8));
    let dio1 = T::from_f64(365.0) * s.inventory_1 / s.cogs_1.max_val(T::from_f64(1e-8));
    let dpo1 = T::from_f64(365.0) * s.accounts_payable_1 / s.cogs_1.max_val(T::from_f64(1e-8));
    let ccc1 = dso1 + dio1 - dpo1;
    let dso0 =
        T::from_f64(365.0) * s.accounts_receivable_0 / s.revenue_0.max_val(T::from_f64(1e-8));
    let dio0 = T::from_f64(365.0) * s.inventory_0 / s.cogs_0.max_val(T::from_f64(1e-8));
    let dpo0 = T::from_f64(365.0) * s.accounts_payable_0 / s.cogs_0.max_val(T::from_f64(1e-8));
    let ccc0 = dso0 + dio0 - dpo0;
    let noa_avg = (s.net_operating_assets_1 + s.net_operating_assets_0) / T::from_f64(2.0);
    let sloan = (s.net_operating_assets_1 - s.net_operating_assets_0)
        / noa_avg.abs().max_val(T::from_f64(1e-8));
    let cfo_ratio = (s.net_income_1 - s.cfo_curr)
        / (s.net_income_1.abs() + noa_avg.abs()).max_val(T::from_f64(1e-8));
    let delta_wc = (s.accounts_receivable_1 - s.accounts_receivable_0)
        + (s.inventory_1 - s.inventory_0)
        - (s.accounts_payable_1 - s.accounts_payable_0);
    let wc_accrual_scaled = delta_wc / noa_avg.abs().max_val(T::from_f64(1e-8));
    let dd_pred = (T::from_f64(0.20) * s.cfo_prev - T::from_f64(0.50) * s.cfo_curr
        + T::from_f64(0.20) * s.cfo_next)
        / noa_avg.abs().max_val(T::from_f64(1e-8));
    let resid =
        ((wc_accrual_scaled - dd_pred) * (wc_accrual_scaled - dd_pred) + T::from_f64(1e-24)).sqrt();
    let ccc_diff = ((ccc1 - ccc0) * (ccc1 - ccc0) + T::from_f64(1e-24)).sqrt();
    let sloan_abs = (sloan * sloan + T::from_f64(1e-24)).sqrt();
    let cfo_abs = (cfo_ratio * cfo_ratio + T::from_f64(1e-24)).sqrt();
    let z = ccc_diff / T::from_f64(365.0) + sloan_abs + cfo_abs + resid;
    SchilitSloanDechowDichevQualityReport {
        days_sales_outstanding_1: dso1,
        days_inventory_on_hand_1: dio1,
        days_payable_outstanding_1: dpo1,
        cash_conversion_cycle_0: ccc0,
        cash_conversion_cycle_1: ccc1,
        ccc_deterioration_days: ccc1 - ccc0,
        sloan_balance_sheet_accrual_ratio: sloan,
        cash_flow_accrual_ratio: cfo_ratio,
        dechow_dichev_accrual_residual: resid,
        schilit_shenanigan_composite_z: z,
    }
}
#[cfg(test)]
mod phase32_tests {
    use super::*;
    use ql_ad::Dual;

    fn caplet_spec<T: Real>(x: T) -> GairatPolalaCheyetteCapletSpec<T> {
        GairatPolalaCheyetteCapletSpec {
            forward_rate: T::from_f64(0.03) + x,
            strike_rate: T::from_f64(0.032),
            fixing_time_t1: T::one(),
            payment_time_t2: T::from_f64(1.5),
            mean_reversion_kappa: T::from_f64(0.1),
            initial_sigma0: T::from_f64(0.015),
            initial_skew_b: T::from_f64(-0.2),
            target_atm_normal_vol: T::from_f64(0.01),
            target_otm_normal_vol: T::from_f64(0.0105),
            discount_factor_t2: T::from_f64(0.96),
        }
    }

    fn qrh_spec<T: Real>(x: T) -> GatheralJusselinQuadraticRoughHestonSpec<T> {
        GatheralJusselinQuadraticRoughHestonSpec {
            forward: T::from_f64(100.0) + x,
            strike: T::from_f64(102.0),
            expiry: T::one(),
            hurst_h: T::from_f64(0.1),
            weight_a: T::from_f64(0.2),
            asymmetry_center_b: T::from_f64(0.3),
            floor_variance_c: T::from_f64(0.02),
            initial_state_z0: T::from_f64(0.1),
            kernel_eta: T::from_f64(0.5),
            mean_reversion_lambda: T::one(),
        }
    }

    fn qual_spec<T: Real>(x: T) -> SchilitSloanDechowDichevQualitySpec<T> {
        SchilitSloanDechowDichevQualitySpec {
            revenue_0: T::from_f64(100.0),
            revenue_1: T::from_f64(105.0) + x,
            cogs_0: T::from_f64(60.0),
            cogs_1: T::from_f64(62.0),
            accounts_receivable_0: T::from_f64(10.0),
            accounts_receivable_1: T::from_f64(14.0),
            inventory_0: T::from_f64(8.0),
            inventory_1: T::from_f64(11.0),
            accounts_payable_0: T::from_f64(7.0),
            accounts_payable_1: T::from_f64(7.0),
            net_operating_assets_0: T::from_f64(50.0),
            net_operating_assets_1: T::from_f64(58.0),
            net_income_1: T::from_f64(10.0),
            cfo_prev: T::from_f64(8.0),
            cfo_curr: T::from_f64(6.0),
            cfo_next: T::from_f64(8.5),
        }
    }

    #[test]
    fn phase32_invariants_and_dual_fd() {
        let h = 1e-6;
        let c0 = evaluate_gairat_polala_cheyette_local_vol_caplet(&caplet_spec(0.0));
        assert!(c0.deterministic_auxiliary_y_bar > 0.0 && c0.calibration_rmse < 1e-4);
        let fd_c = (evaluate_gairat_polala_cheyette_local_vol_caplet(&caplet_spec(h))
            .bachelier_caplet_pv
            - evaluate_gairat_polala_cheyette_local_vol_caplet(&caplet_spec(-h))
                .bachelier_caplet_pv)
            / (2.0 * h);
        let ad_c = evaluate_gairat_polala_cheyette_local_vol_caplet(&caplet_spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .bachelier_caplet_pv;
        assert!(ad_c.d[0].abs() > 1e-6 && (ad_c.d[0] - fd_c).abs() < 1e-5);

        let q0 = evaluate_gatheral_jusselin_quadratic_rough_heston(&qrh_spec(0.0));
        assert!(q0.instantaneous_variance_v0 >= 0.02 && q0.qrh_atm_skew < 0.0);
        let fd_q = (evaluate_gatheral_jusselin_quadratic_rough_heston(&qrh_spec(h)).qrh_call_price
            - evaluate_gatheral_jusselin_quadratic_rough_heston(&qrh_spec(-h)).qrh_call_price)
            / (2.0 * h);
        let ad_q = evaluate_gatheral_jusselin_quadratic_rough_heston(&qrh_spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .qrh_call_price;
        assert!(ad_q.d[0].abs() > 1e-6 && (ad_q.d[0] - fd_q).abs() < 1e-5);

        let s0 = evaluate_schilit_sloan_dechow_dichev_earnings_quality(&qual_spec(0.0));
        assert!(
            (s0.cash_conversion_cycle_1
                - (s0.days_sales_outstanding_1 + s0.days_inventory_on_hand_1
                    - s0.days_payable_outstanding_1))
                .abs()
                < 1e-12
        );
        assert!(s0.ccc_deterioration_days > 0.0);
        let fd_s = (evaluate_schilit_sloan_dechow_dichev_earnings_quality(&qual_spec(h))
            .schilit_shenanigan_composite_z
            - evaluate_schilit_sloan_dechow_dichev_earnings_quality(&qual_spec(-h))
                .schilit_shenanigan_composite_z)
            / (2.0 * h);
        let ad_s = evaluate_schilit_sloan_dechow_dichev_earnings_quality(&qual_spec(
            Dual::<1>::variable(0.0, 0),
        ))
        .schilit_shenanigan_composite_z;
        assert!(ad_s.d[0].abs() > 1e-6 && (ad_s.d[0] - fd_s).abs() < 1e-5);
    }
}
