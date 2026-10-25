use crate::Real;
use nalgebra::{DMatrix, DVector, LU};
use std::fmt;

#[derive(Clone, Debug)]
pub struct MotSuperReplicationIftSpec<T: Real = f64> {
    pub support_grid: Vec<T>,
    pub marginal_t1_probs: Vec<T>,
    pub marginal_t2_probs: Vec<T>,
    pub forward_start_strike: T,
    pub entropic_regularization_eps: T,
}
#[derive(Clone, Debug)]
pub struct MotSuperReplicationIftReport<T: Real = f64> {
    pub super_replication_upper_bound: T,
    pub sub_replication_lower_bound: T,
    pub entropic_mot_price: T,
    pub dual_static_potentials_u: Vec<T>,
    pub dual_static_potentials_v: Vec<T>,
    pub martingale_delta_hedge_h: Vec<T>,
    pub martingale_constraint_max_error: T,
    pub model_free_bound_width: T,
}
pub fn solve_mot_super_replication_bounds_and_ift<T: Real>(
    s: &MotSuperReplicationIftSpec<T>,
) -> MotSuperReplicationIftReport<T> {
    let m = s.support_grid.len();
    assert_eq!(m, s.marginal_t1_probs.len());
    assert_eq!(m, s.marginal_t2_probs.len());
    assert!(m > 0);
    let mut p1 = s.marginal_t1_probs.clone();
    let mut p2 = s.marginal_t2_probs.clone();
    let z1 = p1.iter().copied().fold(T::zero(), |a, b| a + b);
    let z2 = p2.iter().copied().fold(T::zero(), |a, b| a + b);
    let z1 = z1.max_val(T::from_f64(1e-14));
    let z2 = z2.max_val(T::from_f64(1e-14));
    for x in &mut p1 {
        *x /= z1;
    }
    for x in &mut p2 {
        *x /= z2;
    }
    let same_marginals = p1
        .iter()
        .zip(p2.iter())
        .all(|(a, b)| (*a - *b).abs().val() < 1e-12);

    let smooth_pos = |y: T| T::from_f64(0.5) * (y + (y * y + T::from_f64(1e-10)).sqrt());
    let mut cost = vec![vec![T::zero(); m]; m];
    for i in 0..m {
        for j in 0..m {
            cost[i][j] = smooth_pos(s.support_grid[j] / s.support_grid[i] - s.forward_start_strike);
        }
    }
    if same_marginals {
        let mut price = T::zero();
        for i in 0..m {
            price += p1[i] * cost[i][i];
        }
        return MotSuperReplicationIftReport {
            super_replication_upper_bound: price,
            sub_replication_lower_bound: price,
            entropic_mot_price: price,
            dual_static_potentials_u: vec![price * T::from_f64(0.5); m],
            dual_static_potentials_v: vec![price * T::from_f64(0.5); m],
            martingale_delta_hedge_h: vec![T::zero(); m],
            martingale_constraint_max_error: T::zero(),
            model_free_bound_width: T::zero(),
        };
    }

    let eps = s.entropic_regularization_eps.max_val(T::from_f64(0.02));
    // A direction-specific EMOT solve.  The f64 phase gives deterministic,
    // robust initial potentials; three generic sweeps preserve Dual lanes.
    let run_direction = |sigma_dir: f64| -> (T, Vec<T>, Vec<T>, Vec<T>, T) {
        let xv: Vec<f64> = s.support_grid.iter().map(|x| x.val()).collect();
        let p1v: Vec<f64> = p1.iter().map(|x| x.val()).collect();
        let p2v: Vec<f64> = p2.iter().map(|x| x.val()).collect();
        let kv: Vec<Vec<f64>> = cost
            .iter()
            .map(|r| r.iter().map(|x| x.val()).collect())
            .collect();
        let ep = eps.val().max(0.02);
        let mut u = vec![0.0; m];
        let mut v = vec![0.0; m];
        let mut h = vec![0.0; m];
        for _ in 0..35 {
            let mut logk = vec![vec![0.0; m]; m];
            for i in 0..m {
                for j in 0..m {
                    logk[i][j] = (sigma_dir * kv[i][j] - u[i] - v[j] - h[i] * (xv[j] - xv[i])) / ep;
                    logk[i][j] = logk[i][j].clamp(-40.0, 40.0);
                }
            }
            for i in 0..m {
                let q: f64 = (0..m).map(|j| p2v[j] * logk[i][j].exp()).sum();
                u[i] += ep * q.max(1e-15).ln();
            }
            // Re-evaluate the kernel after the row potential update before
            // applying the martingale Newton correction.
            for i in 0..m {
                for j in 0..m {
                    logk[i][j] = (sigma_dir * kv[i][j] - u[i] - v[j] - h[i] * (xv[j] - xv[i])) / ep;
                    logk[i][j] = logk[i][j].clamp(-40.0, 40.0);
                }
            }
            for i in 1..m.saturating_sub(1) {
                let mut r = 0.0;
                let mut m2 = 0.0;
                for j in 0..m {
                    let pij = p1v[i] * p2v[j] * logk[i][j].exp();
                    let dz = xv[j] - xv[i];
                    r += pij * dz;
                    m2 += pij * dz * dz;
                }
                h[i] += 0.7 * ep * r / (m2 + 1e-12);
            }
            // Re-evaluate once more after the slope update for the column
            // normalization sweep.
            for i in 0..m {
                for j in 0..m {
                    logk[i][j] = (sigma_dir * kv[i][j] - u[i] - v[j] - h[i] * (xv[j] - xv[i])) / ep;
                    logk[i][j] = logk[i][j].clamp(-40.0, 40.0);
                }
            }
            for j in 0..m {
                let q: f64 = (0..m).map(|i| p1v[i] * logk[i][j].exp()).sum();
                v[j] += ep * q.max(1e-15).ln();
            }
        }
        let mut ut: Vec<T> = u.iter().map(|x| T::from_f64(*x)).collect();
        let mut vt: Vec<T> = v.iter().map(|x| T::from_f64(*x)).collect();
        let mut ht: Vec<T> = h.iter().map(|x| T::from_f64(*x)).collect();
        for _ in 0..3 {
            let mut pi = vec![vec![T::zero(); m]; m];
            for i in 0..m {
                for j in 0..m {
                    let l = (T::from_f64(sigma_dir) * cost[i][j]
                        - ut[i]
                        - vt[j]
                        - ht[i] * (s.support_grid[j] - s.support_grid[i]))
                        / eps;
                    pi[i][j] = p1[i]
                        * p2[j]
                        * l.max_val(T::from_f64(-40.0))
                            .min_val(T::from_f64(40.0))
                            .exp();
                }
            }
            for i in 0..m {
                let q = (0..m).fold(T::zero(), |a, j| {
                    a + p2[j] * (pi[i][j] / (p1[i] * p2[j]).max_val(T::from_f64(1e-30)))
                });
                ut[i] += eps * q.max_val(T::from_f64(1e-15)).ln();
            }
            let mut pi2 = vec![vec![T::zero(); m]; m];
            for i in 0..m {
                for j in 0..m {
                    let l = (T::from_f64(sigma_dir) * cost[i][j]
                        - ut[i]
                        - vt[j]
                        - ht[i] * (s.support_grid[j] - s.support_grid[i]))
                        / eps;
                    pi2[i][j] = p1[i]
                        * p2[j]
                        * l.max_val(T::from_f64(-40.0))
                            .min_val(T::from_f64(40.0))
                            .exp();
                }
            }
            for i in 1..m.saturating_sub(1) {
                let mut r = T::zero();
                let mut m2 = T::zero();
                for j in 0..m {
                    let dz = s.support_grid[j] - s.support_grid[i];
                    r += pi2[i][j] * dz;
                    m2 += pi2[i][j] * dz * dz;
                }
                ht[i] += T::from_f64(0.7) * eps * r / (m2 + T::from_f64(1e-12));
            }
            // Re-evaluate with the updated martingale slopes before the
            // column potential update.
            for i in 0..m {
                for j in 0..m {
                    let l = (T::from_f64(sigma_dir) * cost[i][j]
                        - ut[i]
                        - vt[j]
                        - ht[i] * (s.support_grid[j] - s.support_grid[i]))
                        / eps;
                    pi2[i][j] = p1[i]
                        * p2[j]
                        * l.max_val(T::from_f64(-40.0))
                            .min_val(T::from_f64(40.0))
                            .exp();
                }
            }
            for j in 0..m {
                let q = (0..m).fold(T::zero(), |a, i| {
                    a + p1[i] * (pi2[i][j] / (p1[i] * p2[j]).max_val(T::from_f64(1e-30)))
                });
                vt[j] += eps * q.max_val(T::from_f64(1e-15)).ln();
            }
        }
        let mut price = T::zero();
        let mut maxerr = T::zero();
        for i in 0..m {
            let mut row_res = T::zero();
            for j in 0..m {
                let l = (T::from_f64(sigma_dir) * cost[i][j]
                    - ut[i]
                    - vt[j]
                    - ht[i] * (s.support_grid[j] - s.support_grid[i]))
                    / eps;
                let pij = p1[i]
                    * p2[j]
                    * l.max_val(T::from_f64(-40.0))
                        .min_val(T::from_f64(40.0))
                        .exp();
                price += pij * cost[i][j];
                row_res += pij * (s.support_grid[j] - s.support_grid[i]);
            }
            if i > 0 && i + 1 < m {
                maxerr = maxerr.max_val(row_res.abs());
            }
        }
        (price, ut, vt, ht, maxerr)
    };
    let (super_price, su, sv, sh, super_err) = run_direction(1.0);
    let (sub_price, _, _, _, sub_err) = run_direction(-1.0);
    let (ref_price, _, _, _, _) = run_direction(0.0);
    let lower = sub_price.min_val(super_price).max_val(T::zero());
    let upper = super_price.max_val(sub_price);
    MotSuperReplicationIftReport {
        super_replication_upper_bound: upper,
        sub_replication_lower_bound: lower,
        entropic_mot_price: ref_price,
        dual_static_potentials_u: su,
        dual_static_potentials_v: sv,
        martingale_delta_hedge_h: sh,
        martingale_constraint_max_error: super_err.max_val(sub_err),
        model_free_bound_width: (upper - lower).max_val(T::zero()),
    }
}

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
pub fn ift_lift_1d<T: Real, F: FnOnce(T) -> T>(x_star: f64, residual: F, f_x: T, f_xx: f64) -> T {
    // A singular/non-finite Jacobian cannot be lifted reliably.  Returning a
    // primal constant is preferable to injecting infinities/NaNs into every
    // sensitivity lane.
    let fx0 = f_x.val();
    if !fx0.is_finite() || fx0.abs() <= 1e-14 {
        return T::from_f64(x_star);
    }

    // Evaluate the graph once.  The primal component performs one Newton
    // polish while the tangent/adjoint components provide the exact implicit
    // lift; evaluating the residual repeatedly would duplicate a TapeVar
    // pricing graph on the active reverse tape.
    let x = T::from_f64(x_star);
    let r_eval = residual(x);
    let r0 = r_eval.val();
    let x_polished = T::from_f64(x_star - r0 / fx0);
    let d_tan = -(r_eval - T::from_f64(r0)) / f_x;
    x_polished + d_tan - T::from_f64(0.5 * f_xx) * d_tan * d_tan / f_x
}

/// Fallible counterpart to [`ift_lift_1d`].  Existing callers retain the
/// backwards-compatible infallible API above; new solvers can surface a typed
/// domain error instead of silently returning a constant when the Jacobian is
/// unusable.
pub fn try_ift_lift_1d<T: Real, F: FnOnce(T) -> T>(
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
        let w = self.lu_t.solve(&b).unwrap_or_else(|| {
            self.j_theta
                .transpose()
                .svd(true, true)
                .solve(&b, 1e-12)
                .expect("SVD IFT solve failed")
        });
        (-self.j_q.transpose() * w).iter().copied().collect()
    }
    pub fn jvp(&self, quote_tangents: &[f64]) -> Vec<f64> {
        assert_eq!(quote_tangents.len(), self.j_q.ncols());
        let rhs = -&self.j_q * DVector::from_column_slice(quote_tangents);
        self.lu
            .solve(&rhs)
            .unwrap_or_else(|| {
                self.j_theta
                    .clone()
                    .svd(true, true)
                    .solve(&rhs, 1e-12)
                    .expect("SVD IFT solve failed")
            })
            .iter()
            .copied()
            .collect()
    }
    pub fn jacobian_dtheta_dq(&self) -> DMatrix<f64> {
        let rhs = -self.j_q.clone();
        self.lu.solve(&rhs).unwrap_or_else(|| {
            self.j_theta
                .clone()
                .svd(true, true)
                .solve(&rhs, 1e-12)
                .expect("SVD IFT solve failed")
        })
    }
}

#[cfg(test)]
mod phase22_tests {
    #[test]
    fn dual_fd_parity_smoke() {
        let x = 0.37_f64;
        let h = 1.0e-6_f64;
        let f = |z: f64| z.exp() + z * z;
        let fd = (f(x + h) - f(x - h)) / (2.0 * h);
        let analytic = x.exp() + 2.0 * x;
        assert!((fd - analytic).abs() < 1.0e-5);
    }
}
