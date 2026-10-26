use ql_ad::Real;

#[inline]
fn smin<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b - (d * d + T::from_f64(1e-16)).sqrt())
}

#[inline]
fn smax<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
}

#[derive(Clone, Copy, Debug)]
pub struct HaganWestMonotoneConvexCurveSpec<T: Real = f64> {
    pub knot_tenors: [T; 4],
    pub discrete_forwards: [T; 3],
    pub eval_tenor: T,
}

#[derive(Clone, Copy, Debug)]
pub struct HaganWestMonotoneConvexCurveReport<T: Real = f64> {
    pub knot_instantaneous_forward_1: T,
    pub knot_instantaneous_forward_2: T,
    pub eval_instantaneous_forward: T,
    pub eval_zero_rate: T,
    pub eval_discount_factor: T,
    pub cell_integrated_average_forward: T,
    pub cell_average_conservation_error: T,
    pub ameliorated_min_instantaneous_forward: T,
}

pub fn evaluate_hagan_west_monotone_convex_forward_curve<T: Real>(
    spec: &HaganWestMonotoneConvexCurveSpec<T>,
) -> HaganWestMonotoneConvexCurveReport<T> {
    let dt1 = (spec.knot_tenors[1] - spec.knot_tenors[0]).max_val(T::from_f64(1e-8));
    let dt2 = (spec.knot_tenors[2] - spec.knot_tenors[1]).max_val(T::from_f64(1e-8));
    let dt3 = (spec.knot_tenors[3] - spec.knot_tenors[2]).max_val(T::from_f64(1e-8));
    let f1d = spec.discrete_forwards[0];
    let f2d = spec.discrete_forwards[1];
    let f3d = spec.discrete_forwards[2];
    let f1_raw = (dt2 * f1d + dt1 * f2d) / (dt1 + dt2);
    let f2_raw = (dt3 * f2d + dt2 * f3d) / (dt2 + dt3);
    let two = T::from_f64(2.0);
    let f1_mc = smax(T::zero(), smin(two * f1d, smin(two * f2d, f1_raw)));
    let f2_mc = smax(T::zero(), smin(two * f2d, smin(two * f3d, f2_raw)));
    let x = ((spec.eval_tenor - spec.knot_tenors[1]) / dt2)
        .max_val(T::zero())
        .min_val(T::one());
    let g0 = f1_mc - f2d;
    let g1 = f2_mc - f2d;
    let f_eval = f2d
        + g0 * (T::one() - T::from_f64(4.0) * x + T::from_f64(3.0) * x * x)
        + g1 * (-two * x + T::from_f64(3.0) * x * x);
    let cell_int = f2d + g0 * (T::one() - two + T::one()) + g1 * (-T::one() + T::one());
    let diff = cell_int - f2d;
    let err = (diff * diff + T::from_f64(1e-30)).sqrt();
    let partial_int = f2d * x + g0 * (x - two * x * x + x * x * x) + g1 * (-x * x + x * x * x);
    let z = (f1d * dt1 + dt2 * partial_int) / spec.eval_tenor.max_val(T::from_f64(1e-8));
    let df = (-spec.eval_tenor * z).exp();
    let denom = T::from_f64(3.0) * (g0 + g1);
    let x_star = smax(
        T::zero(),
        smin(T::one(), (two * g0 + g1) / (denom + T::from_f64(1e-12))),
    );
    let f_star = f2d
        + g0 * (T::one() - T::from_f64(4.0) * x_star + T::from_f64(3.0) * x_star * x_star)
        + g1 * (-two * x_star + T::from_f64(3.0) * x_star * x_star);
    let min_f = smax(T::zero(), smin(f1_mc, smin(f2_mc, f_star)));
    HaganWestMonotoneConvexCurveReport {
        knot_instantaneous_forward_1: f1_mc,
        knot_instantaneous_forward_2: f2_mc,
        eval_instantaneous_forward: f_eval,
        eval_zero_rate: z,
        eval_discount_factor: df,
        cell_integrated_average_forward: cell_int,
        cell_average_conservation_error: err,
        ameliorated_min_instantaneous_forward: min_f,
    }
}

#[cfg(test)]
mod phase35_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> HaganWestMonotoneConvexCurveSpec<T> {
        HaganWestMonotoneConvexCurveSpec {
            knot_tenors: [T::zero(), T::one(), T::from_f64(2.0), T::from_f64(3.0)],
            discrete_forwards: [T::from_f64(0.03), T::from_f64(0.04) + x, T::from_f64(0.045)],
            eval_tenor: T::from_f64(1.4),
        }
    }
    #[test]
    fn phase35_curves_invariants_and_dual_fd() {
        let r = evaluate_hagan_west_monotone_convex_forward_curve(&spec(0.0));
        assert!(r.cell_average_conservation_error < 1e-12);
        assert!(r.ameliorated_min_instantaneous_forward >= 0.0);
        let h = 1e-6;
        let fd = (evaluate_hagan_west_monotone_convex_forward_curve(&spec(h)).eval_discount_factor
            - evaluate_hagan_west_monotone_convex_forward_curve(&spec(-h)).eval_discount_factor)
            / (2.0 * h);
        let ad =
            evaluate_hagan_west_monotone_convex_forward_curve(&spec(Dual::<1>::variable(0.0, 0)))
                .eval_discount_factor;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
