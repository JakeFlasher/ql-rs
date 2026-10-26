#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct LandoJltMarkovMigrationSpec<T: Real = f64> {
    pub ig_to_hy_rate: T,
    pub ig_to_default_rate: T,
    pub hy_to_ig_rate: T,
    pub hy_to_ccc_rate: T,
    pub hy_to_default_rate: T,
    pub ccc_to_hy_rate: T,
    pub ccc_to_default_rate: T,
    pub risk_premium_multiplier_pi: T,
    pub recovery_rate: T,
    pub risk_free_rate: T,
    pub horizon_years: T,
}

#[derive(Clone, Copy, Debug)]
pub struct LandoJltMarkovMigrationReport<T: Real = f64> {
    pub ig_transition_probs_q: [T; 4],
    pub hy_transition_probs_q: [T; 4],
    pub ccc_transition_probs_q: [T; 4],
    pub row_sum_stochastic_error: T,
    pub ig_risky_zero_bond_price: T,
    pub hy_risky_zero_bond_price: T,
    pub ccc_risky_zero_bond_price: T,
    pub ig_credit_spread_bps: T,
    pub hy_credit_spread_bps: T,
    pub ig_to_hy_downgrade_price_loss: T,
}

fn matmul4<T: Real>(a: &[[T; 4]; 4], b: &[[T; 4]; 4]) -> [[T; 4]; 4] {
    let mut c = [[T::zero(); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            let mut s = T::zero();
            for k in 0..4 {
                s += a[i][k] * b[k][j];
            }
            c[i][j] = s;
        }
    }
    c
}

pub fn evaluate_lando_jlt_markov_rating_migration<T: Real>(
    spec: &LandoJltMarkovMigrationSpec<T>,
) -> LandoJltMarkovMigrationReport<T> {
    let pi = spec.risk_premium_multiplier_pi;
    let l01 = pi * spec.ig_to_hy_rate;
    let l03 = pi * spec.ig_to_default_rate;
    let l10 = pi * spec.hy_to_ig_rate;
    let l12 = pi * spec.hy_to_ccc_rate;
    let l13 = pi * spec.hy_to_default_rate;
    let l21 = pi * spec.ccc_to_hy_rate;
    let l23 = pi * spec.ccc_to_default_rate;

    let z = T::zero();
    let generator = [
        [-(l01 + l03), l01, z, l03],
        [l10, -(l10 + l12 + l13), l12, l13],
        [z, l21, -(l21 + l23), l23],
        [z, z, z, z],
    ];
    let scale = spec.horizon_years.max_val(T::from_f64(1e-6)) * T::from_f64(0.25);
    let mut a = [[z; 4]; 4];
    let mut e0 = [[z; 4]; 4];
    let mut term = [[z; 4]; 4];
    for i in 0..4 {
        e0[i][i] = T::one();
        term[i][i] = T::one();
        for j in 0..4 {
            a[i][j] = generator[i][j] * scale;
        }
    }
    for m in 1..=8 {
        let next = matmul4(&term, &a);
        let inv_m = T::one() / T::from_f64(m as f64);
        for i in 0..4 {
            for j in 0..4 {
                term[i][j] = next[i][j] * inv_m;
                e0[i][j] += term[i][j];
            }
        }
    }
    let e1 = matmul4(&e0, &e0);
    let p_q = matmul4(&e1, &e1);

    let mut max_err = T::zero();
    for i in 0..4 {
        let row_sum = p_q[i][0] + p_q[i][1] + p_q[i][2] + p_q[i][3];
        let d = T::one() - row_sum;
        let err_i = (d * d + T::from_f64(1e-30)).sqrt();
        max_err = max_err.max_val(err_i);
    }

    let mat = spec.horizon_years.max_val(T::from_f64(1e-6));
    let df_rf = (-spec.risk_free_rate * mat).exp();
    let lgd = T::one() - spec.recovery_rate;
    let surv_ig = (T::one() - lgd * p_q[0][3]).max_val(T::from_f64(1e-8));
    let surv_hy = (T::one() - lgd * p_q[1][3]).max_val(T::from_f64(1e-8));
    let surv_ccc = (T::one() - lgd * p_q[2][3]).max_val(T::from_f64(1e-8));

    let b_ig = df_rf * surv_ig;
    let b_hy = df_rf * surv_hy;
    let b_ccc = df_rf * surv_ccc;
    let s_ig_bps = -surv_ig.ln() / mat * T::from_f64(1e4);
    let s_hy_bps = -surv_hy.ln() / mat * T::from_f64(1e4);

    LandoJltMarkovMigrationReport {
        ig_transition_probs_q: p_q[0],
        hy_transition_probs_q: p_q[1],
        ccc_transition_probs_q: p_q[2],
        row_sum_stochastic_error: max_err,
        ig_risky_zero_bond_price: b_ig,
        hy_risky_zero_bond_price: b_hy,
        ccc_risky_zero_bond_price: b_ccc,
        ig_credit_spread_bps: s_ig_bps,
        hy_credit_spread_bps: s_hy_bps,
        ig_to_hy_downgrade_price_loss: b_ig - b_hy,
    }
}

#[cfg(test)]
mod phase36_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn lando_jlt_markov_migration_invariants_and_ad_parity() {
        let spec = LandoJltMarkovMigrationSpec {
            ig_to_hy_rate: 0.045,
            ig_to_default_rate: 0.002,
            hy_to_ig_rate: 0.030,
            hy_to_ccc_rate: 0.080,
            hy_to_default_rate: 0.025,
            ccc_to_hy_rate: 0.040,
            ccc_to_default_rate: 0.140,
            risk_premium_multiplier_pi: 1.45,
            recovery_rate: 0.40,
            risk_free_rate: 0.035,
            horizon_years: 3.0,
        };
        let rep = evaluate_lando_jlt_markov_rating_migration(&spec);
        assert!(rep.row_sum_stochastic_error < 1e-12);
        assert!(rep.ig_credit_spread_bps < rep.hy_credit_spread_bps);
        assert!(rep.ig_to_hy_downgrade_price_loss > 0.0);

        let h = 1e-6;
        let r_ad = evaluate_lando_jlt_markov_rating_migration(&LandoJltMarkovMigrationSpec {
            ig_to_hy_rate: Dual::<1>::variable(spec.ig_to_hy_rate, 0),
            ig_to_default_rate: Dual::<1>::from_value(spec.ig_to_default_rate),
            hy_to_ig_rate: Dual::<1>::from_value(spec.hy_to_ig_rate),
            hy_to_ccc_rate: Dual::<1>::from_value(spec.hy_to_ccc_rate),
            hy_to_default_rate: Dual::<1>::from_value(spec.hy_to_default_rate),
            ccc_to_hy_rate: Dual::<1>::from_value(spec.ccc_to_hy_rate),
            ccc_to_default_rate: Dual::<1>::from_value(spec.ccc_to_default_rate),
            risk_premium_multiplier_pi: Dual::<1>::from_value(spec.risk_premium_multiplier_pi),
            recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
            risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
            horizon_years: Dual::<1>::from_value(spec.horizon_years),
        });
        let mut sp = spec;
        sp.ig_to_hy_rate += h;
        let mut sm = spec;
        sm.ig_to_hy_rate -= h;
        let _fd = (evaluate_lando_jlt_markov_rating_migration(&sp).ig_credit_spread_bps
            - evaluate_lando_jlt_markov_rating_migration(&sm).ig_credit_spread_bps)
            / (2.0 * h);
        assert!((r_ad.ig_credit_spread_bps.d[0] - _fd).abs() < 1e-5);
    }
}
