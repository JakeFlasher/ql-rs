use crate::Real;
use nalgebra::{DMatrix, DVector, LU};
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdError {
    SingularJacobian,
    NonFiniteJacobian,
}
impl fmt::Display for AdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SingularJacobian => f.write_str("singular IFT Jacobian"),
            Self::NonFiniteJacobian => f.write_str("non-finite IFT Jacobian"),
        }
    }
}
impl std::error::Error for AdError {}

/// Lift a converged scalar root into an AD scalar. The residual is evaluated once
/// at the converged root; `f_x` is the scalar derivative with respect to the root.
pub fn ift_lift_1d<T: Real, F: Fn(T) -> T>(x_star: f64, residual: F, f_x: T, f_xx: f64) -> T {
    // A singular/non-finite Jacobian cannot be lifted reliably.  Returning a
    // primal constant is preferable to injecting infinities/NaNs into every
    // sensitivity lane.
    let fx0 = f_x.val();
    if !fx0.is_finite() || fx0.abs() <= 1e-14 {
        return T::from_f64(x_star);
    }

    // Polish the scalar root once, but only using primal values.  This removes
    // residual noise from the AD lift while retaining a useful correction when
    // a numerical root solver stopped just short of convergence.
    let r0 = residual(T::from_f64(x_star)).val();
    let x_polished = if r0.is_finite() && r0.abs() > 1e-15 {
        let step = (-r0 / fx0).clamp(-0.5 * x_star.abs().max(1.0), 0.5 * x_star.abs().max(1.0));
        let candidate = x_star + step;
        let r_candidate = residual(T::from_f64(candidate)).val();
        if r_candidate.is_finite() && r_candidate.abs() < r0.abs() {
            candidate
        } else {
            x_star
        }
    } else {
        x_star
    };

    let x = T::from_f64(x_polished);
    let r_eval = residual(x);
    // Keep only the tangent/adjoint part of the residual.  Its primal value is
    // identically zero, so d has zero primal and d*d cannot create a spurious
    // first-order term from 2*d0*d1 in Dual/Dual2 arithmetic.
    let r_diff = r_eval - T::from_f64(r_eval.val());
    let d = -r_diff / f_x;
    x + d - T::from_f64(0.5 * f_xx) * d * d / f_x
}

/// Fallible counterpart to [`ift_lift_1d`].  Existing callers retain the
/// backwards-compatible infallible API above; new solvers can surface a typed
/// domain error instead of silently returning a constant when the Jacobian is
/// unusable.
pub fn try_ift_lift_1d<T: Real, F: Fn(T) -> T>(
    x_star: f64,
    residual: F,
    f_x: T,
    f_xx: f64,
) -> Result<T, AdError> {
    let fx0 = f_x.val();
    if !fx0.is_finite() {
        return Err(AdError::NonFiniteJacobian);
    }
    if fx0.abs() <= 1e-14 {
        return Err(AdError::SingularJacobian);
    }
    Ok(ift_lift_1d(x_star, residual, f_x, f_xx))
}

#[derive(Debug, Clone)]
pub struct IftSystemSensitivity {
    pub j_theta: DMatrix<f64>,
    pub j_q: DMatrix<f64>,
    pub lu: LU<f64, nalgebra::Dyn, nalgebra::Dyn>,
    pub lu_t: LU<f64, nalgebra::Dyn, nalgebra::Dyn>,
}
impl IftSystemSensitivity {
    pub fn new(j_theta: DMatrix<f64>, j_q: DMatrix<f64>) -> Self {
        assert_eq!(j_theta.nrows(), j_theta.ncols(), "J_theta must be square");
        assert_eq!(j_theta.nrows(), j_q.nrows(), "IFT Jacobian row mismatch");
        let lu = j_theta.clone().lu();
        let lu_t = j_theta.transpose().lu();
        Self {
            j_theta,
            j_q,
            lu,
            lu_t,
        }
    }
    pub fn from_matrices(
        j_theta: &[f64],
        rows: usize,
        cols: usize,
        j_q: &[f64],
        q_cols: usize,
    ) -> Self {
        Self::new(
            DMatrix::from_row_slice(rows, cols, j_theta),
            DMatrix::from_row_slice(rows, q_cols, j_q),
        )
    }
    pub fn vjp(&self, pillar_adjoints: &[f64]) -> Vec<f64> {
        assert_eq!(pillar_adjoints.len(), self.j_theta.nrows());
        let b = DVector::from_column_slice(pillar_adjoints);
        // Reuse the factorization cached at construction.  Solving the
        // transposed system through a freshly-built LU both wastes work and
        // can produce different pivoting near singularity.
        let w = self.lu_t.solve(&b).expect("singular IFT Jacobian");
        (-self.j_q.transpose() * w).iter().copied().collect()
    }
    pub fn jvp(&self, quote_tangents: &[f64]) -> Vec<f64> {
        assert_eq!(quote_tangents.len(), self.j_q.ncols());
        let rhs = -&self.j_q * DVector::from_column_slice(quote_tangents);
        self.lu
            .solve(&rhs)
            .expect("singular IFT Jacobian")
            .iter()
            .copied()
            .collect()
    }
    pub fn jacobian_dtheta_dq(&self) -> DMatrix<f64> {
        let rhs = -self.j_q.clone();
        self.lu.solve(&rhs).expect("singular IFT Jacobian")
    }
}
