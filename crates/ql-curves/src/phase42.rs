//! Darbyshire (2016) & Henrard (2014) Central-Bank Meeting-Dated Step-OIS Curve
//! & Turn-of-Year (`ToY`) Repo Spike Decomposition Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct DarbyshireCentralBankStepOisTurnSpec<T: Real = f64> {
    pub current_overnight_ois_r0: T,
    pub meeting_1_time_years_t1: T,
    pub meeting_1_expected_step_bps: T,
    pub meeting_2_time_years_t2: T,
    pub meeting_2_expected_step_bps: T,
    pub turn_of_year_spike_bps: T,
    pub turn_effective_day_weight_tau: T,
    pub horizon_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct DarbyshireCentralBankStepOisTurnReport<T: Real = f64> {
    pub terminal_policy_ois_rate: T,
    pub integrated_step_ois_zero_rate: T,
    pub turn_of_year_zero_rate_addon_bps: T,
    pub all_in_discount_factor_p0t: T,
    pub strip_without_turn_discount_factor: T,
    pub meeting_1_mp01_per_bp: T,
    pub meeting_2_mp01_per_bp: T,
}

pub fn evaluate_darbyshire_central_bank_step_ois_and_turn<T: Real>(
    spec: &DarbyshireCentralBankStepOisTurnSpec<T>,
) -> DarbyshireCentralBankStepOisTurnReport<T> {
    let r0 = spec.current_overnight_ois_r0;
    let t = spec.horizon_years_t.max_val(T::from_f64(0.10));
    let t1 = spec
        .meeting_1_time_years_t1
        .max_val(T::from_f64(0.02))
        .min_val(t);
    let t2 = spec.meeting_2_time_years_t2.max_val(t1).min_val(t);

    let step1 = spec.meeting_1_expected_step_bps * T::from_f64(1e-4);
    let step2 = spec.meeting_2_expected_step_bps * T::from_f64(1e-4);
    let r1 = r0 + step1;
    let r2 = r1 + step2;

    // Integrated step-function OIS rate over [0, T]:
    // int_0^T r(u) du = r0 * t + step1 * (t - t1) + step2 * (t - t2)
    let int_step = r0 * t + step1 * (t - t1) + step2 * (t - t2);
    let zero_step = int_step / t;

    // Turn-of-year (ToY) spike contribution over effective year fraction tau_turn:
    let tau_turn = spec
        .turn_effective_day_weight_tau
        .max_val(T::from_f64(1.0 / 365.0));
    let int_turn = (spec.turn_of_year_spike_bps * T::from_f64(1e-4)) * tau_turn;
    let turn_zero_bps = (int_turn / t) * T::from_f64(1e4);

    let df_no_turn = (-int_step).exp();
    let df_all_in = (-(int_step + int_turn)).exp();

    // Meeting-date Jacobian sensitivities (MP01 per 1bp step on a 10,000 notional OIS):
    let mp01_1 = (t - t1) * df_all_in;
    let mp01_2 = (t - t2) * df_all_in;

    DarbyshireCentralBankStepOisTurnReport {
        terminal_policy_ois_rate: r2,
        integrated_step_ois_zero_rate: zero_step,
        turn_of_year_zero_rate_addon_bps: turn_zero_bps,
        all_in_discount_factor_p0t: df_all_in,
        strip_without_turn_discount_factor: df_no_turn,
        meeting_1_mp01_per_bp: mp01_1,
        meeting_2_mp01_per_bp: mp01_2,
    }
}

#[cfg(test)]
mod phase42_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn darbyshire_central_bank_step_ois_and_turn_invariants_and_ad_parity() {
        let spec = DarbyshireCentralBankStepOisTurnSpec {
            current_overnight_ois_r0: 0.045,
            meeting_1_time_years_t1: 0.15,
            meeting_1_expected_step_bps: 22.0,
            meeting_2_time_years_t2: 0.35,
            meeting_2_expected_step_bps: 18.0,
            turn_of_year_spike_bps: 35.0,
            turn_effective_day_weight_tau: 3.0 / 360.0,
            horizon_years_t: 0.50,
        };
        let rep = evaluate_darbyshire_central_bank_step_ois_and_turn(&spec);
        assert!(rep.terminal_policy_ois_rate > spec.current_overnight_ois_r0);
        assert!(rep.all_in_discount_factor_p0t < rep.strip_without_turn_discount_factor);
        assert!(rep.meeting_1_mp01_per_bp > rep.meeting_2_mp01_per_bp);

        let h = 1e-6;
        let r_ad = evaluate_darbyshire_central_bank_step_ois_and_turn(
            &DarbyshireCentralBankStepOisTurnSpec {
                current_overnight_ois_r0: Dual::<1>::variable(spec.current_overnight_ois_r0, 0),
                meeting_1_time_years_t1: Dual::<1>::from_value(spec.meeting_1_time_years_t1),
                meeting_1_expected_step_bps: Dual::<1>::from_value(
                    spec.meeting_1_expected_step_bps,
                ),
                meeting_2_time_years_t2: Dual::<1>::from_value(spec.meeting_2_time_years_t2),
                meeting_2_expected_step_bps: Dual::<1>::from_value(
                    spec.meeting_2_expected_step_bps,
                ),
                turn_of_year_spike_bps: Dual::<1>::from_value(spec.turn_of_year_spike_bps),
                turn_effective_day_weight_tau: Dual::<1>::from_value(
                    spec.turn_effective_day_weight_tau,
                ),
                horizon_years_t: Dual::<1>::from_value(spec.horizon_years_t),
            },
        );
        let mut sp = spec;
        sp.current_overnight_ois_r0 += h;
        let mut sm = spec;
        sm.current_overnight_ois_r0 -= h;
        let fd = (evaluate_darbyshire_central_bank_step_ois_and_turn(&sp)
            .all_in_discount_factor_p0t
            - evaluate_darbyshire_central_bank_step_ois_and_turn(&sm).all_in_discount_factor_p0t)
            / (2.0 * h);
        assert!((r_ad.all_in_discount_factor_p0t.d[0] - fd).abs() < 1e-7);
    }
}
