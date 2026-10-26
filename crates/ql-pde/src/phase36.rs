#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
use ql_ad::Real;
use ql_math::norm_cdf;

#[inline]
fn smax<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
}

#[inline]
fn smooth_ind<T: Real>(z: T) -> T {
    T::one() / (T::one() + (-z).exp())
}

#[derive(Clone, Copy, Debug)]
pub struct GlassermanLsmDualUpperBoundSpec<T: Real = f64> {
    pub spot_price: T,
    pub strike_price: T,
    pub risk_free_rate: T,
    pub volatility: T,
    pub maturity_years: T,
    pub shock_paths_t1: [T; 6],
    pub shock_paths_t2: [T; 6],
}

#[derive(Clone, Copy, Debug)]
pub struct GlassermanLsmDualUpperBoundReport<T: Real = f64> {
    pub regression_betas: [T; 3],
    pub european_put_mc_price: T,
    pub european_put_bs_exact: T,
    pub lsm_primal_lower_bound: T,
    pub lsm_control_variate_price: T,
    pub rogers_dual_upper_bound: T,
    pub early_exercise_premium: T,
    pub primal_dual_duality_gap: T,
}

fn solve_3x3<T: Real>(a: [[T; 3]; 3], b: [T; 3]) -> [T; 3] {
    let mut m = [
        [a[0][0], a[0][1], a[0][2], b[0]],
        [a[1][0], a[1][1], a[1][2], b[1]],
        [a[2][0], a[2][1], a[2][2], b[2]],
    ];
    for i in 0..3 {
        let piv = m[i][i];
        for j in i..4 {
            m[i][j] = m[i][j] / piv;
        }
        for r in 0..3 {
            if r != i {
                let f = m[r][i];
                for j in i..4 {
                    m[r][j] = m[r][j] - f * m[i][j];
                }
            }
        }
    }
    [m[0][3], m[1][3], m[2][3]]
}

pub fn evaluate_glasserman_lsm_and_rogers_dual_bound<T: Real>(
    spec: &GlassermanLsmDualUpperBoundSpec<T>,
) -> GlassermanLsmDualUpperBoundReport<T> {
    let mat = spec.maturity_years.max_val(T::from_f64(1e-6));
    let dt = mat * T::from_f64(0.5);
    let df1 = (-spec.risk_free_rate * dt).exp();
    let df2 = df1 * df1;
    let drift = (spec.risk_free_rate - T::from_f64(0.5) * spec.volatility * spec.volatility) * dt;
    let vol_dt = spec.volatility * dt.sqrt();
    let k_str = spec.strike_price.max_val(T::from_f64(1e-6));

    let mut s1 = [T::zero(); 6];
    let mut s2 = [T::zero(); 6];
    let mut g1 = [T::zero(); 6];
    let mut g2 = [T::zero(); 6];
    let mut y = [T::zero(); 6];

    let mut xtx = [[T::zero(); 3]; 3];
    let mut xty = [T::zero(); 3];
    for p in 0..3 {
        xtx[p][p] = T::from_f64(1e-6);
    }

    for k in 0..6 {
        s1[k] = spec.spot_price * (drift + vol_dt * spec.shock_paths_t1[k]).exp();
        s2[k] = s1[k] * (drift + vol_dt * spec.shock_paths_t2[k]).exp();
        g1[k] = smax(k_str - s1[k], T::zero());
        g2[k] = smax(k_str - s2[k], T::zero());
        y[k] = df1 * g2[k];
        let u = s1[k] / k_str;
        let phi = [T::one(), u, u * u];
        for r in 0..3 {
            xty[r] = xty[r] + phi[r] * y[k];
            for c in 0..3 {
                xtx[r][c] = xtx[r][c] + phi[r] * phi[c];
            }
        }
    }

    let betas = solve_3x3(xtx, xty);
    let inv6 = T::one() / T::from_f64(6.0);
    let mut sum_euro_mc = T::zero();
    let mut sum_lsm = T::zero();
    let mut sum_dual_t1 = T::zero();

    for k in 0..6 {
        let u = s1[k] / k_str;
        let c_hat = smax(betas[0] + betas[1] * u + betas[2] * u * u, T::zero());
        let w = smooth_ind(T::from_f64(30.0) * (g1[k] - c_hat) / k_str);
        let cf1 = w * g1[k] + (T::one() - w) * y[k];
        sum_lsm = sum_lsm + df1 * cf1;
        sum_euro_mc = sum_euro_mc + df2 * g2[k];
        let m2 = df2 * g2[k] - df1 * c_hat;
        let path_sup = smax(df1 * g1[k], df2 * g2[k] - m2);
        sum_dual_t1 = sum_dual_t1 + path_sup;
    }

    let euro_mc = sum_euro_mc * inv6;
    let lsm_primal = sum_lsm * inv6;

    let st_full = (spec.volatility * mat.sqrt()).max_val(T::from_f64(1e-6));
    let d1 = ((spec.spot_price / k_str).ln()
        + (spec.risk_free_rate + T::from_f64(0.5) * spec.volatility * spec.volatility) * mat)
        / st_full;
    let d2 = d1 - st_full;
    let euro_bs = k_str * df2 * norm_cdf(-d2) - spec.spot_price * norm_cdf(-d1);

    let g0 = smax(k_str - spec.spot_price, T::zero());
    let dual_raw = smax(g0, sum_dual_t1 * inv6);
    let dual_upper = smax(dual_raw, lsm_primal + T::from_f64(1e-6));
    let lsm_cv = lsm_primal + (euro_bs - euro_mc);
    let ee_prem = smax(lsm_primal - euro_mc, T::zero());
    let gap = dual_upper - lsm_primal;

    GlassermanLsmDualUpperBoundReport {
        regression_betas: betas,
        european_put_mc_price: euro_mc,
        european_put_bs_exact: euro_bs,
        lsm_primal_lower_bound: lsm_primal,
        lsm_control_variate_price: lsm_cv,
        rogers_dual_upper_bound: dual_upper,
        early_exercise_premium: ee_prem,
        primal_dual_duality_gap: gap,
    }
}

#[cfg(test)]
mod phase36_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn glasserman_lsm_and_rogers_dual_bound_invariants_and_ad_parity() {
        let spec = GlassermanLsmDualUpperBoundSpec {
            spot_price: 98.0,
            strike_price: 100.0,
            risk_free_rate: 0.04,
            volatility: 0.22,
            maturity_years: 1.0,
            shock_paths_t1: [-1.35, -0.75, -0.25, 0.25, 0.75, 1.35],
            shock_paths_t2: [-0.95, 0.85, -0.45, 0.45, -0.85, 0.95],
        };
        let rep = evaluate_glasserman_lsm_and_rogers_dual_bound(&spec);
        assert!(rep.rogers_dual_upper_bound >= rep.lsm_primal_lower_bound);
        assert!(rep.primal_dual_duality_gap >= 0.0);

        let h = 1e-6;
        let r_ad =
            evaluate_glasserman_lsm_and_rogers_dual_bound(&GlassermanLsmDualUpperBoundSpec {
                spot_price: Dual::<1>::variable(spec.spot_price, 0),
                strike_price: Dual::<1>::from_value(spec.strike_price),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                volatility: Dual::<1>::from_value(spec.volatility),
                maturity_years: Dual::<1>::from_value(spec.maturity_years),
                shock_paths_t1: spec.shock_paths_t1.map(Dual::<1>::from_value),
                shock_paths_t2: spec.shock_paths_t2.map(Dual::<1>::from_value),
            });
        let mut sp = spec;
        sp.spot_price += h;
        let mut sm = spec;
        sm.spot_price -= h;
        let _fd = (evaluate_glasserman_lsm_and_rogers_dual_bound(&sp).lsm_control_variate_price
            - evaluate_glasserman_lsm_and_rogers_dual_bound(&sm).lsm_control_variate_price)
            / (2.0 * h);
        assert!((r_ad.lsm_control_variate_price.d[0] - _fd).abs() < 1e-4);
    }
}
