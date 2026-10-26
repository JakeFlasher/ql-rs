use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct GregorySaCvaEulerNettingSpec<T: Real = f64> {
    pub cva_delta_sensitivities: [T; 2],
    pub cva_hedge_sensitivities: [T; 2],
    pub delta_risk_weights: [T; 2],
    pub intra_bucket_corr_rho: T,
    pub cva_vega_sensitivity: T,
    pub vega_risk_weight: T,
    pub subportfolio_exposures: [T; 3],
    pub exposure_corr_3x3: [[T; 3]; 3],
    pub cpty_hazard_rate: T,
    pub cpty_lgd: T,
    pub horizon_years: T,
}

#[derive(Clone, Copy, Debug)]
pub struct GregorySaCvaEulerNettingReport<T: Real = f64> {
    pub sa_cva_delta_capital: T,
    pub sa_cva_vega_capital: T,
    pub total_sa_cva_capital: T,
    pub standalone_cva_sum: T,
    pub netted_portfolio_cva: T,
    pub netting_cva_reduction_benefit: T,
    pub euler_allocated_cva: [T; 3],
    pub netting_set_euler_allocation_error: T,
}

pub fn evaluate_gregory_sa_cva_and_euler_netting_allocation<T: Real>(
    s: &GregorySaCvaEulerNettingSpec<T>,
) -> GregorySaCvaEulerNettingReport<T> {
    let ws0 =
        s.delta_risk_weights[0] * (s.cva_delta_sensitivities[0] - s.cva_hedge_sensitivities[0]);
    let ws1 =
        s.delta_risk_weights[1] * (s.cva_delta_sensitivities[1] - s.cva_hedge_sensitivities[1]);
    let h0 = s.delta_risk_weights[0] * s.cva_hedge_sensitivities[0];
    let h1 = s.delta_risk_weights[1] * s.cva_hedge_sensitivities[1];
    let kd = T::from_f64(1.25)
        * (ws0 * ws0
            + ws1 * ws1
            + T::from_f64(2.0) * s.intra_bucket_corr_rho * ws0 * ws1
            + T::from_f64(0.01) * (h0 * h0 + h1 * h1))
            .sqrt();
    let v_prod = s.vega_risk_weight * s.cva_vega_sensitivity;
    let kv = T::from_f64(1.25) * (v_prod * v_prod + T::from_f64(1e-24)).sqrt();
    let mut quad = T::zero();
    for i in 0..3 {
        for j in 0..3 {
            quad += s.exposure_corr_3x3[i][j]
                * s.subportfolio_exposures[i]
                * s.subportfolio_exposures[j];
        }
    }
    let e_net = quad.max_val(T::from_f64(1e-12)).sqrt();
    let alpha_cva = s.cpty_lgd * (T::one() - (-s.cpty_hazard_rate * s.horizon_years).exp());
    let mut euler = [T::zero(); 3];
    for i in 0..3 {
        let mut row = T::zero();
        for j in 0..3 {
            row += s.exposure_corr_3x3[i][j] * s.subportfolio_exposures[j];
        }
        euler[i] = alpha_cva * s.subportfolio_exposures[i] * row / e_net;
    }
    let cva_net = alpha_cva * e_net;
    let cva_stand = alpha_cva
        * (s.subportfolio_exposures[0] + s.subportfolio_exposures[1] + s.subportfolio_exposures[2]);
    let diff = cva_net - (euler[0] + euler[1] + euler[2]);
    GregorySaCvaEulerNettingReport {
        sa_cva_delta_capital: kd,
        sa_cva_vega_capital: kv,
        total_sa_cva_capital: kd + kv,
        standalone_cva_sum: cva_stand,
        netted_portfolio_cva: cva_net,
        netting_cva_reduction_benefit: cva_stand - cva_net,
        euler_allocated_cva: euler,
        netting_set_euler_allocation_error: (diff * diff + T::from_f64(1e-30)).sqrt(),
    }
}

#[cfg(test)]
mod phase35_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> GregorySaCvaEulerNettingSpec<T> {
        GregorySaCvaEulerNettingSpec {
            cva_delta_sensitivities: [T::from_f64(10.0) + x, T::from_f64(8.0)],
            cva_hedge_sensitivities: [T::from_f64(6.0), T::from_f64(5.0)],
            delta_risk_weights: [T::from_f64(0.05), T::from_f64(0.06)],
            intra_bucket_corr_rho: T::from_f64(0.5),
            cva_vega_sensitivity: T::from_f64(4.0),
            vega_risk_weight: T::from_f64(0.1),
            subportfolio_exposures: [T::from_f64(20.0) + x, T::from_f64(15.0), T::from_f64(10.0)],
            exposure_corr_3x3: [
                [T::one(), T::from_f64(0.3), T::from_f64(0.1)],
                [T::from_f64(0.3), T::one(), T::from_f64(0.2)],
                [T::from_f64(0.1), T::from_f64(0.2), T::one()],
            ],
            cpty_hazard_rate: T::from_f64(0.02),
            cpty_lgd: T::from_f64(0.6),
            horizon_years: T::from_f64(3.0),
        }
    }
    #[test]
    fn phase35_xva_invariants_and_dual_fd() {
        let r = evaluate_gregory_sa_cva_and_euler_netting_allocation(&spec(0.0));
        assert!(r.netting_set_euler_allocation_error < 1e-12);
        assert!(r.netting_cva_reduction_benefit > 0.0);
        let h = 1e-6;
        let fd = (evaluate_gregory_sa_cva_and_euler_netting_allocation(&spec(h))
            .total_sa_cva_capital
            - evaluate_gregory_sa_cva_and_euler_netting_allocation(&spec(-h)).total_sa_cva_capital)
            / (2.0 * h);
        let ad = evaluate_gregory_sa_cva_and_euler_netting_allocation(&spec(Dual::<1>::variable(
            0.0, 0,
        )))
        .total_sa_cva_capital;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
