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
fn sig_segment_range<T: Real>(t: &[T], x: &[T], start: usize, end: usize) -> (T, T, T, T, T, T) {
    let mut s1 = T::zero();
    let mut s2 = T::zero();
    let mut s11 = T::zero();
    let mut s12 = T::zero();
    let mut s21 = T::zero();
    let mut s22 = T::zero();
    let half = T::from_f64(0.5);
    for k in start..end {
        let dt = t[k + 1] - t[k];
        let dx = x[k + 1] - x[k];
        s11 += s1 * dt + half * dt * dt;
        s12 += s1 * dx + half * dt * dx;
        s21 += s2 * dt + half * dx * dt;
        s22 += s2 * dx + half * dx * dx;
        s1 += dt;
        s2 += dx;
    }
    (s1, s2, s11, s12, s21, s22)
}
pub fn evaluate_bayer_path_signature_chen_and_levy_area<T: Real>(
    s: &BayerPathSignatureChenLevyAreaSpec<T>,
) -> BayerPathSignatureChenLevyAreaReport<T> {
    let m = s.time_grid.len().min(s.spot_path.len()).saturating_sub(1);
    let (s1, s2, s11, s12, s21, s22) = sig_segment_range(&s.time_grid, &s.spot_path, 0, m);
    let mid = m / 2;
    let (a1, a2, a11, a12, a21, a22) = sig_segment_range(&s.time_grid, &s.spot_path, 0, mid);
    let (b1, b2, b11, b12, b21, b22) = sig_segment_range(&s.time_grid, &s.spot_path, mid, m);
    let c11 = s11 - (a11 + b11 + a1 * b1);
    let c12 = s12 - (a12 + b12 + a1 * b2);
    let c21 = s21 - (a21 + b21 + a2 * b1);
    let c22 = s22 - (a22 + b22 + a2 * b2);
    let chen_err = (c11 * c11 + c12 * c12 + c21 * c21 + c22 * c22 + T::from_f64(1e-30)).sqrt();
    let shuf = (s12 + s21) - s1 * s2;
    let shuf_err = (shuf * shuf + T::from_f64(1e-30)).sqrt();
    let levy = T::from_f64(0.5) * (s12 - s21);
    let w = &s.linear_functional_weights;
    let payoff = w[0] * s1 + w[1] * s2 + w[2] * s11 + w[3] * s12 + w[4] * s21 + w[5] * s22;
    BayerPathSignatureChenLevyAreaReport {
        sig1_time_increment: s1,
        sig1_spot_increment: s2,
        sig2_time_spot_12: s12,
        sig2_spot_time_21: s21,
        sig2_spot_spot_22: s22,
        signed_levy_area: levy,
        chen_identity_l1_error: chen_err,
        shuffle_product_identity_error: shuf_err,
        signature_linear_payoff: payoff,
    }
}
