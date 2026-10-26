//! Piazzesi & Swanson (2008) & Gürkaynak–Sack–Wright (2007) FOMC Fed-Funds Futures
//! Rate-Hike Probability Tree & Risk-Adjusted Policy Path Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct PiazzesiSwansonFomcFuturesHikeSpec<T: Real = f64> {
    pub current_target_fed_funds_r0: T,
    pub futures_implied_month_avg_rate_f0: T,
    pub meeting_day_of_month_d: T,
    pub total_days_in_month_n: T,
    pub standard_hike_step_delta: T,
    pub piazzesi_swanson_term_premium_bps: T,
    pub next_meeting_futures_rate_f1: T,
}

#[derive(Debug, Clone, Copy)]
pub struct PiazzesiSwansonFomcFuturesHikeReport<T: Real = f64> {
    pub implied_post_meeting_rate_q: T,
    pub risk_neutral_hike_probability_q: T,
    pub physical_post_meeting_rate_p: T,
    pub physical_hike_probability_p: T,
    pub risk_premium_probability_wedge: T,
    pub policy_rate_bernoulli_std_bps: T,
    pub cumulative_two_meeting_expected_hikes: T,
}

pub fn evaluate_piazzesi_swanson_fomc_futures_hike_tree<T: Real>(
    spec: &PiazzesiSwansonFomcFuturesHikeSpec<T>,
) -> PiazzesiSwansonFomcFuturesHikeReport<T> {
    let r0 = spec.current_target_fed_funds_r0;
    let f0 = spec.futures_implied_month_avg_rate_f0;
    let n_days = spec.total_days_in_month_n.max_val(T::from_f64(28.0));
    let d_meet = spec
        .meeting_day_of_month_d
        .max_val(T::one())
        .min_val(n_days - T::one());
    let step = spec.standard_hike_step_delta.max_val(T::from_f64(0.0005));

    // Gürkaynak-Sack-Wright (2007) intra-month day-weight decomposition:
    // F_0 = (d / N) * r_0 + ((N - d) / N) * E^Q[r_1]
    let w_pre = d_meet / n_days;
    let w_post = (n_days - d_meet) / n_days;
    let r1_q = (f0 - w_pre * r0) / w_post;
    let prob_q = ((r1_q - r0) / step).max_val(T::zero()).min_val(T::one());

    // Piazzesi-Swanson (2008) risk-premium adjustment: E^P[r_1] = E^Q[r_1] - xi_t:
    let xi = spec.piazzesi_swanson_term_premium_bps * T::from_f64(1e-4);
    let r1_p = r1_q - xi;
    let prob_p = ((r1_p - r0) / step).max_val(T::zero()).min_val(T::one());
    let wedge = prob_q - prob_p;

    let var_bernoulli = (prob_p * (T::one() - prob_p)).max_val(T::from_f64(1e-8));
    let std_bps = step * var_bernoulli.sqrt() * T::from_f64(1e4);
    let two_meeting_hikes = (spec.next_meeting_futures_rate_f1 - xi - r0) / step;

    PiazzesiSwansonFomcFuturesHikeReport {
        implied_post_meeting_rate_q: r1_q,
        risk_neutral_hike_probability_q: prob_q,
        physical_post_meeting_rate_p: r1_p,
        physical_hike_probability_p: prob_p,
        risk_premium_probability_wedge: wedge,
        policy_rate_bernoulli_std_bps: std_bps,
        cumulative_two_meeting_expected_hikes: two_meeting_hikes,
    }
}

#[cfg(test)]
mod phase41_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn piazzesi_swanson_fomc_futures_hike_tree_invariants_and_ad_parity() {
        let spec = PiazzesiSwansonFomcFuturesHikeSpec {
            current_target_fed_funds_r0: 0.0475,
            futures_implied_month_avg_rate_f0: 0.0486,
            meeting_day_of_month_d: 15.0,
            total_days_in_month_n: 30.0,
            standard_hike_step_delta: 0.0025,
            piazzesi_swanson_term_premium_bps: 4.5,
            next_meeting_futures_rate_f1: 0.0512,
        };
        let rep = evaluate_piazzesi_swanson_fomc_futures_hike_tree(&spec);
        assert!(rep.risk_neutral_hike_probability_q > rep.physical_hike_probability_p);
        assert!(rep.physical_hike_probability_p > 0.0 && rep.physical_hike_probability_p < 1.0);
        assert!(rep.cumulative_two_meeting_expected_hikes > 1.0);

        let h = 1e-6;
        let r_ad =
            evaluate_piazzesi_swanson_fomc_futures_hike_tree(&PiazzesiSwansonFomcFuturesHikeSpec {
                current_target_fed_funds_r0: Dual::<1>::from_value(
                    spec.current_target_fed_funds_r0,
                ),
                futures_implied_month_avg_rate_f0: Dual::<1>::variable(
                    spec.futures_implied_month_avg_rate_f0,
                    0,
                ),
                meeting_day_of_month_d: Dual::<1>::from_value(spec.meeting_day_of_month_d),
                total_days_in_month_n: Dual::<1>::from_value(spec.total_days_in_month_n),
                standard_hike_step_delta: Dual::<1>::from_value(spec.standard_hike_step_delta),
                piazzesi_swanson_term_premium_bps: Dual::<1>::from_value(
                    spec.piazzesi_swanson_term_premium_bps,
                ),
                next_meeting_futures_rate_f1: Dual::<1>::from_value(
                    spec.next_meeting_futures_rate_f1,
                ),
            });
        let mut sp = spec;
        sp.futures_implied_month_avg_rate_f0 += h;
        let mut sm = spec;
        sm.futures_implied_month_avg_rate_f0 -= h;
        let fd = (evaluate_piazzesi_swanson_fomc_futures_hike_tree(&sp)
            .physical_hike_probability_p
            - evaluate_piazzesi_swanson_fomc_futures_hike_tree(&sm).physical_hike_probability_p)
            / (2.0 * h);
        assert!((r_ad.physical_hike_probability_p.d[0] - fd).abs() < 1e-6);
    }
}
