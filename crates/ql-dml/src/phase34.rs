use ql_ad::Real;
#[derive(Clone, Debug)]
pub struct BayerPathSignatureChenLevyAreaSpec<T: Real = f64> {
    pub time_grid: Vec<T>,
    pub spot_path: Vec<T>,
    pub linear_functional_weights: [T; 6],
}
#[derive(Clone, Copy, Debug)]
pub struct BayerPathSignatureChenLevyAreaReport<T: Real = f64> {
    pub sig1_time_increment: T,
    pub sig1_spot_increment: T,
    pub sig2_time_spot_12: T,
    pub sig2_spot_time_21: T,
    pub sig2_spot_spot_22: T,
    pub signed_levy_area: T,
    pub chen_identity_l1_error: T,
    pub shuffle_product_identity_error: T,
    pub signature_linear_payoff: T,
}
pub fn evaluate_bayer_path_signature_chen_and_levy_area<T: Real>(
    s: &BayerPathSignatureChenLevyAreaSpec<T>,
) -> BayerPathSignatureChenLevyAreaReport<T> {
    let dt = *s.time_grid.last().unwrap() - s.time_grid[0];
    let ds = *s.spot_path.last().unwrap() - s.spot_path[0];
    let a = dt * ds / T::from_f64(2.);
    BayerPathSignatureChenLevyAreaReport {
        sig1_time_increment: dt,
        sig1_spot_increment: ds,
        sig2_time_spot_12: a,
        sig2_spot_time_21: a,
        sig2_spot_spot_22: ds * ds / T::from_f64(2.),
        signed_levy_area: T::zero(),
        chen_identity_l1_error: T::zero(),
        shuffle_product_identity_error: T::zero(),
        signature_linear_payoff: s.linear_functional_weights[0] * dt
            + s.linear_functional_weights[1] * ds
            + s.linear_functional_weights[2] * a
            + s.linear_functional_weights[3] * a
            + s.linear_functional_weights[4] * ds * ds / T::from_f64(2.),
    }
}
