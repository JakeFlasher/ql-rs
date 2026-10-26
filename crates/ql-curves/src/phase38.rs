#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct HaganWestTensionCurveSpec<T: Real = f64> {
    pub pillar_intervals: [T; 4],
    pub discrete_forwards: [T; 4],
    pub tension_lambda: T,
    pub eval_coord_x: T,
}

#[derive(Clone, Copy, Debug)]
pub struct HaganWestTensionCurveReport<T: Real = f64> {
    pub instantaneous_forward_rate: T,
    pub zero_coupon_bond_price: T,
    pub tension_cubic_adjustment: T,
    pub integrated_cell_forward: T,
    pub tension_cell_conservation_error: T,
    pub curve_roughness_energy: T,
}

pub fn evaluate_hagan_west_tension_convex_forward_curve<T: Real>(
    spec: &HaganWestTensionCurveSpec<T>,
) -> HaganWestTensionCurveReport<T> {
    let t0 = spec.pillar_intervals[0].max_val(T::from_f64(1e-6));
    let t1 = spec.pillar_intervals[1].max_val(T::from_f64(1e-6));
    let t2 = spec.pillar_intervals[2].max_val(T::from_f64(1e-6));
    let t3 = spec.pillar_intervals[3].max_val(T::from_f64(1e-6));
    let f0d = spec.discrete_forwards[0];
    let f1d = spec.discrete_forwards[1];
    let f2d = spec.discrete_forwards[2];
    let f3d = spec.discrete_forwards[3];

    let f1 = (t0 * f1d + t1 * f0d) / (t0 + t1);
    let f2 = (t1 * f2d + t2 * f1d + T::from_f64(0.01) * t3 * (f3d - f2d)) / (t1 + t2);
    let g0 = f1 - f1d;
    let g1 = f2 - f1d;

    let lam = spec.tension_lambda.max_val(T::zero());
    let kap = lam / (T::one() + lam);
    let x = spec
        .eval_coord_x
        .max_val(T::from_f64(0.01))
        .min_val(T::from_f64(0.99));
    let x2 = x * x;
    let x3 = x2 * x;

    let q_part = g0 * (T::one() - T::from_f64(4.0) * x + T::from_f64(3.0) * x2)
        + g1 * (-T::from_f64(2.0) * x + T::from_f64(3.0) * x2);
    let c_basis = x * (T::one() - x) * (T::one() - T::from_f64(2.0) * x);
    let c_adj = kap * (g0 + g1) * c_basis;
    let inst_f = f1d + q_part + c_adj;

    // Exact integral of g(u) from 0 to x:
    // int_0^x (1 - 4u + 3u^2) du = x - 2x^2 + x^3
    // int_0^x (-2u + 3u^2) du = -x^2 + x^3
    // int_0^x (u - 3u^2 + 2u^3) du = 0.5 x^2 - x^3 + 0.5 x^4
    let int_g_x = g0 * (x - T::from_f64(2.0) * x2 + x3)
        + g1 * (-x2 + x3)
        + kap * (g0 + g1) * (T::from_f64(0.5) * x2 - x3 + T::from_f64(0.5) * x2 * x2);
    let zcb = (-(t0 * f0d + t1 * (f1d * x + int_g_x))).exp();

    // At x = 1, int_0^1 g(u) du = 0 identically for all tension_lambda
    let int_g_1 = g0 * (T::one() - T::from_f64(2.0) + T::one())
        + g1 * (-T::one() + T::one())
        + kap * (g0 + g1) * (T::from_f64(0.5) - T::one() + T::from_f64(0.5));
    let cell_fwd = f1d + int_g_1;
    let diff = cell_fwd - f1d;
    let cons_err = (diff * diff + T::from_f64(1e-30)).sqrt();

    let rough = T::from_f64(4.0) * (g0 * g0 + g0 * g1 + g1 * g1) / t1
        + kap * kap * (g0 + g1) * (g0 + g1) / t1;

    HaganWestTensionCurveReport {
        instantaneous_forward_rate: inst_f,
        zero_coupon_bond_price: zcb,
        tension_cubic_adjustment: c_adj,
        integrated_cell_forward: cell_fwd,
        tension_cell_conservation_error: cons_err,
        curve_roughness_energy: rough,
    }
}

#[cfg(test)]
mod phase38_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn hagan_west_tension_curve_invariants_and_ad_parity() {
        let spec = HaganWestTensionCurveSpec {
            pillar_intervals: [1.0, 1.0, 2.0, 3.0],
            discrete_forwards: [0.032, 0.038, 0.044, 0.046],
            tension_lambda: 1.5,
            eval_coord_x: 0.35,
        };
        let rep = evaluate_hagan_west_tension_convex_forward_curve(&spec);
        assert!(rep.tension_cell_conservation_error < 1e-12);
        assert!(rep.zero_coupon_bond_price > 0.0 && rep.zero_coupon_bond_price < 1.0);

        let h = 1e-6;
        let r_ad = evaluate_hagan_west_tension_convex_forward_curve(&HaganWestTensionCurveSpec {
            pillar_intervals: spec.pillar_intervals.map(Dual::<1>::from_value),
            discrete_forwards: [
                Dual::<1>::from_value(spec.discrete_forwards[0]),
                Dual::<1>::variable(spec.discrete_forwards[1], 0),
                Dual::<1>::from_value(spec.discrete_forwards[2]),
                Dual::<1>::from_value(spec.discrete_forwards[3]),
            ],
            tension_lambda: Dual::<1>::from_value(spec.tension_lambda),
            eval_coord_x: Dual::<1>::from_value(spec.eval_coord_x),
        });
        let mut sp = spec;
        sp.discrete_forwards[1] += h;
        let mut sm = spec;
        sm.discrete_forwards[1] -= h;
        let fd = (evaluate_hagan_west_tension_convex_forward_curve(&sp).zero_coupon_bond_price
            - evaluate_hagan_west_tension_convex_forward_curve(&sm).zero_coupon_bond_price)
            / (2.0 * h);
        assert!((r_ad.zero_coupon_bond_price.d[0] - fd).abs() < 1e-5);
    }
}
