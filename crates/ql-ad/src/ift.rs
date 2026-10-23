use crate::Real;
use nalgebra::{DMatrix, DVector, LU};

/// Lift a converged scalar root into an AD scalar. The residual is evaluated once
/// at the converged root; `f_x` is the scalar derivative with respect to the root.
pub fn ift_lift_1d<T: Real, F: Fn(T) -> T>(x_star: f64, residual: F, f_x: T, f_xx: f64) -> T {
    let x = T::from_f64(x_star);
    let d = -residual(x) / f_x;
    x + d - T::from_f64(0.5 * f_xx) * d * d / f_x
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
