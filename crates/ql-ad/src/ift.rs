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
    // A direction-specific EMOT solve.  Keep the entire Sinkhorn--Newton
    // program in `T`: when `T` is a Dual this preserves the tangent of every
    // iterate instead of seeding a primal-only f64 fixed point and attempting
    // to recover sensitivities with a handful of polishing sweeps.
    let run_direction = |sigma_dir: f64| -> (T, Vec<T>, Vec<T>, Vec<T>, T) {
        let mut ut = vec![T::zero(); m];
        let mut vt = vec![T::zero(); m];
        let mut ht = vec![T::zero(); m];
        let sigma = T::from_f64(sigma_dir);
        for _ in 0..25 {
            let mut logk = vec![vec![T::zero(); m]; m];
            for i in 0..m {
                for j in 0..m {
                    logk[i][j] = (sigma * cost[i][j]
                        - ut[i]
                        - vt[j]
                        - ht[i] * (s.support_grid[j] - s.support_grid[i]))
                        / eps;
                    logk[i][j] = logk[i][j]
                        .max_val(T::from_f64(-40.0))
                        .min_val(T::from_f64(40.0));
                }
            }
            for i in 0..m {
                let q = (0..m).fold(T::zero(), |a, j| a + p2[j] * logk[i][j].exp());
                ut[i] += eps * q.max_val(T::from_f64(1e-15)).ln();
            }
            for i in 0..m {
                for j in 0..m {
                    logk[i][j] = (sigma * cost[i][j]
                        - ut[i]
                        - vt[j]
                        - ht[i] * (s.support_grid[j] - s.support_grid[i]))
                        / eps;
                    logk[i][j] = logk[i][j]
                        .max_val(T::from_f64(-40.0))
                        .min_val(T::from_f64(40.0));
                }
            }
            for i in 1..m.saturating_sub(1) {
                let mut r = T::zero();
                let mut m2 = T::zero();
                for j in 0..m {
                    let pij = p1[i] * p2[j] * logk[i][j].exp();
                    let dz = s.support_grid[j] - s.support_grid[i];
                    r += pij * dz;
                    m2 += pij * dz * dz;
                }
                ht[i] += T::from_f64(0.7) * eps * r / (m2 + T::from_f64(1e-12));
            }
            for i in 0..m {
                for j in 0..m {
                    logk[i][j] = (sigma * cost[i][j]
                        - ut[i]
                        - vt[j]
                        - ht[i] * (s.support_grid[j] - s.support_grid[i]))
                        / eps;
                    logk[i][j] = logk[i][j]
                        .max_val(T::from_f64(-40.0))
                        .min_val(T::from_f64(40.0));
                }
            }
            for j in 0..m {
                let q = (0..m).fold(T::zero(), |a, i| a + p1[i] * logk[i][j].exp());
                vt[j] += eps * q.max_val(T::from_f64(1e-15)).ln();
            }
        }
        let mut price = T::zero();
        let mut maxerr = T::zero();
        for i in 0..m {
            let mut row_res = T::zero();
            for j in 0..m {
                let l = (sigma * cost[i][j]
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

#[derive(Clone, Debug)]
pub struct PenmanDamodaranResidualIncomeIerpSpec<T: Real = f64> {
    pub initial_book_value: T,
    pub forecast_roe_schedule: Vec<T>,
    pub payout_ratio_schedule: Vec<T>,
    pub terminal_roe: T,
    pub terminal_growth_rate: T,
    pub risk_free_rate: T,
    pub equity_beta: T,
    pub observed_market_price: T,
}
#[derive(Clone, Debug)]
pub struct PenmanDamodaranResidualIncomeIerpReport<T: Real = f64> {
    pub implied_erp: T,
    pub implied_cost_of_equity: T,
    pub riv_intrinsic_value: T,
    pub ddm_intrinsic_value: T,
    pub riv_ddm_clean_surplus_error: T,
    pub ohlson_juettner_aeg_value: T,
    pub d_ierp_d_market_price: T,
    pub d_ierp_d_risk_free_rate: T,
    pub d_ierp_d_terminal_growth: T,
}
#[derive(Clone, Copy)]
struct Jet<T: Real> {
    v: T,
    dre: T,
    dg: T,
}
fn j_const<T: Real>(v: T) -> Jet<T> {
    Jet {
        v,
        dre: T::zero(),
        dg: T::zero(),
    }
}
fn j_add<T: Real>(a: Jet<T>, b: Jet<T>) -> Jet<T> {
    Jet {
        v: a.v + b.v,
        dre: a.dre + b.dre,
        dg: a.dg + b.dg,
    }
}
fn j_sub<T: Real>(a: Jet<T>, b: Jet<T>) -> Jet<T> {
    Jet {
        v: a.v - b.v,
        dre: a.dre - b.dre,
        dg: a.dg - b.dg,
    }
}
fn j_mul<T: Real>(a: Jet<T>, b: Jet<T>) -> Jet<T> {
    Jet {
        v: a.v * b.v,
        dre: a.dre * b.v + a.v * b.dre,
        dg: a.dg * b.v + a.v * b.dg,
    }
}
fn j_div<T: Real>(a: Jet<T>, b: Jet<T>) -> Jet<T> {
    let q = b.v * b.v;
    Jet {
        v: a.v / b.v,
        dre: (a.dre * b.v - a.v * b.dre) / q,
        dg: (a.dg * b.v - a.v * b.dg) / q,
    }
}
fn j_exp<T: Real>(a: Jet<T>) -> Jet<T> {
    let e = a.v.exp();
    Jet {
        v: e,
        dre: e * a.dre,
        dg: e * a.dg,
    }
}
fn j_pow_neg<T: Real>(a: Jet<T>, n: usize) -> Jet<T> {
    j_exp(j_const(T::from_f64(-(n as f64))).mul_ln(a))
}
trait JetLn<T: Real> {
    fn mul_ln(self, a: Jet<T>) -> Jet<T>;
}
impl<T: Real> JetLn<T> for Jet<T> {
    fn mul_ln(self, a: Jet<T>) -> Jet<T> {
        let l = Jet {
            v: a.v.ln(),
            dre: a.dre / a.v,
            dg: a.dg / a.v,
        };
        j_mul(self, l)
    }
}

fn riv_val<T: Real>(s: &PenmanDamodaranResidualIncomeIerpSpec<T>, re: T) -> (T, T, T, T, T) {
    let h = s.forecast_roe_schedule.len().max(1);
    let re_j = Jet {
        v: re,
        dre: T::one(),
        dg: T::zero(),
    };
    let g_j = Jet {
        v: s.terminal_growth_rate,
        dre: T::zero(),
        dg: T::one(),
    };
    let one = j_const(T::one());
    let mut book = j_const(s.initial_book_value);
    let mut riv = j_const(s.initial_book_value);
    let mut ddm = j_const(T::zero());
    let mut aeg_sum = j_const(T::zero());
    let mut eprev = j_const(T::zero());
    let mut dprev = j_const(T::zero());
    let mut e1 = j_const(T::zero());
    for i in 0..h {
        let roe = j_const(
            s.forecast_roe_schedule
                .get(i)
                .copied()
                .unwrap_or(s.terminal_roe),
        );
        let pay = j_const(
            s.payout_ratio_schedule
                .get(i)
                .copied()
                .unwrap_or(T::from_f64(0.5)),
        );
        let e = j_mul(roe, book);
        let d = j_mul(pay, e);
        let ri = j_sub(e, j_mul(re_j, book));
        let df = j_pow_neg(j_add(one, re_j), i + 1);
        riv = j_add(riv, j_mul(ri, df));
        ddm = j_add(ddm, j_mul(d, df));
        if i == 0 {
            e1 = e;
        } else {
            let g_h = j_sub(j_add(e, j_mul(re_j, dprev)), j_mul(j_add(one, re_j), eprev));
            aeg_sum = j_add(aeg_sum, j_mul(g_h, j_pow_neg(j_add(one, re_j), i)));
        }
        book = j_add(book, j_sub(e, d));
        eprev = e;
        dprev = d;
    }
    let e_t = j_mul(j_const(s.terminal_roe), book);
    let d_t = j_sub(e_t, j_mul(g_j, book));
    let ri_t = j_sub(e_t, j_mul(re_j, book));
    let den = j_sub(re_j, g_j);
    let df_h = j_pow_neg(j_add(one, re_j), h);
    let term = j_div(ri_t, den);
    riv = j_add(riv, j_mul(term, df_h));
    ddm = j_add(ddm, j_mul(j_div(d_t, den), df_h));
    let g_h1 = j_sub(
        j_add(e_t, j_mul(re_j, dprev)),
        j_mul(j_add(one, re_j), eprev),
    );
    let term_aeg = j_add(j_mul(g_h1, df_h), j_mul(j_div(j_mul(g_j, ri_t), den), df_h));
    let aeg = j_div(j_add(e1, j_add(aeg_sum, term_aeg)), re_j);
    (riv.v, ddm.v, aeg.v, riv.dre, riv.dg)
}

pub fn solve_penman_damodaran_riv_aeg_and_implied_erp_ift<T: Real>(
    s: &PenmanDamodaranResidualIncomeIerpSpec<T>,
) -> PenmanDamodaranResidualIncomeIerpReport<T> {
    let mut erp = T::from_f64(0.055);
    for _ in 0..24 {
        let re = s.risk_free_rate + s.equity_beta * erp;
        let (v, _, _, dv, _) = riv_val(s, re);
        let den = s.equity_beta * dv;
        let den = if den.val().abs() < 1e-12 {
            if den.val() < 0.0 {
                T::from_f64(-1e-12)
            } else {
                T::from_f64(1e-12)
            }
        } else {
            den
        };
        erp -= (v - s.observed_market_price) / den;
    }
    let re = s.risk_free_rate + s.equity_beta * erp;
    let (riv, ddm, aeg, dv, dg) = riv_val(s, re);
    let den = s.equity_beta * dv;
    let den = if den.val().abs() < 1e-12 {
        if den.val() < 0.0 {
            T::from_f64(-1e-12)
        } else {
            T::from_f64(1e-12)
        }
    } else {
        den
    };
    PenmanDamodaranResidualIncomeIerpReport {
        implied_erp: erp,
        implied_cost_of_equity: re,
        riv_intrinsic_value: riv,
        ddm_intrinsic_value: ddm,
        riv_ddm_clean_surplus_error: (riv - ddm).abs(),
        ohlson_juettner_aeg_value: aeg,
        d_ierp_d_market_price: T::one() / den,
        d_ierp_d_risk_free_rate: -T::one() / s.equity_beta,
        d_ierp_d_terminal_growth: -dg / den,
    }
}

#[cfg(test)]
mod domain_suite_34_tests {
    use super::*;
    use crate::Dual;
    #[test]
    fn mot_dual_fd_parity_and_bounds() {
        let make = |k| MotSuperReplicationIftSpec {
            support_grid: vec![0.8, 1.0, 1.2],
            marginal_t1_probs: vec![0.2, 0.5, 0.3],
            marginal_t2_probs: vec![0.25, 0.5, 0.25],
            forward_start_strike: k,
            entropic_regularization_eps: 0.02,
        };
        let h = 1e-5;
        let up = solve_mot_super_replication_bounds_and_ift(&make(0.95 + h));
        let dn = solve_mot_super_replication_bounds_and_ift(&make(0.95 - h));
        let fd = (up.entropic_mot_price - dn.entropic_mot_price) / (2.0 * h);
        let sd = MotSuperReplicationIftSpec::<Dual<1>> {
            support_grid: vec![
                Dual::<1>::from_value(0.8),
                Dual::<1>::from_value(1.0),
                Dual::<1>::from_value(1.2),
            ],
            marginal_t1_probs: vec![
                Dual::<1>::from_value(0.2),
                Dual::<1>::from_value(0.5),
                Dual::<1>::from_value(0.3),
            ],
            marginal_t2_probs: vec![
                Dual::<1>::from_value(0.25),
                Dual::<1>::from_value(0.5),
                Dual::<1>::from_value(0.25),
            ],
            forward_start_strike: Dual::<1>::variable(0.95, 0),
            entropic_regularization_eps: Dual::<1>::from_value(0.02),
        };
        let ad = solve_mot_super_replication_bounds_and_ift(&sd);
        assert!(ad.super_replication_upper_bound.val() >= ad.entropic_mot_price.val());
        assert!(ad.entropic_mot_price.val() >= ad.sub_replication_lower_bound.val());
        assert!(ad.martingale_constraint_max_error.val() < 1e-3);
        assert!((ad.entropic_mot_price.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_36_tests {
    use super::*;
    use crate::Dual;
    fn mk<T: Real>(p: T) -> PenmanDamodaranResidualIncomeIerpSpec<T> {
        PenmanDamodaranResidualIncomeIerpSpec {
            initial_book_value: T::from_f64(100.0),
            forecast_roe_schedule: vec![T::from_f64(0.12), T::from_f64(0.11)],
            payout_ratio_schedule: vec![T::from_f64(0.4), T::from_f64(0.4)],
            terminal_roe: T::from_f64(0.1),
            terminal_growth_rate: T::from_f64(0.03),
            risk_free_rate: T::from_f64(0.03),
            equity_beta: T::from_f64(1.1),
            observed_market_price: p,
        }
    }
    #[test]
    fn ier_suite36() {
        let r = solve_penman_damodaran_riv_aeg_and_implied_erp_ift(&mk(110.0));
        assert!(r.riv_ddm_clean_surplus_error.abs() < 1e-11);
        assert!((r.riv_intrinsic_value - 110.0).abs() < 1e-10);
        assert!(r.d_ierp_d_market_price < 0.0);
        let h = 1e-5;
        let fd = (solve_penman_damodaran_riv_aeg_and_implied_erp_ift(&mk(110.0 + h)).implied_erp
            - solve_penman_damodaran_riv_aeg_and_implied_erp_ift(&mk(110.0 - h)).implied_erp)
            / (2.0 * h);
        let ad =
            solve_penman_damodaran_riv_aeg_and_implied_erp_ift(&mk(Dual::<1>::variable(110.0, 0)))
                .implied_erp;
        assert!((ad.d[0] - fd).abs() < 1e-6);
        let sf = mk(110.0 + 0.0);
        let h = 1e-5;
        let mut up = sf.clone();
        let mut dn = sf.clone();
        up.risk_free_rate += h;
        dn.risk_free_rate -= h;
        let fd_rf = (solve_penman_damodaran_riv_aeg_and_implied_erp_ift(&up).implied_erp
            - solve_penman_damodaran_riv_aeg_and_implied_erp_ift(&dn).implied_erp)
            / (2.0 * h);
        let mut sd = mk(Dual::<1>::from_value(110.0));
        sd.risk_free_rate = Dual::variable(0.03, 0);
        let ad_rf = solve_penman_damodaran_riv_aeg_and_implied_erp_ift(&sd).implied_erp;
        assert!(
            (ad_rf.d[0] - fd_rf).abs() < 1e-9
                && (ad_rf.d[0] - r.d_ierp_d_risk_free_rate).abs() < 1e-9
        );
        up = sf.clone();
        dn = sf.clone();
        up.terminal_growth_rate += h;
        dn.terminal_growth_rate -= h;
        let fd_g = (solve_penman_damodaran_riv_aeg_and_implied_erp_ift(&up).implied_erp
            - solve_penman_damodaran_riv_aeg_and_implied_erp_ift(&dn).implied_erp)
            / (2.0 * h);
        let mut sg = mk(Dual::<1>::from_value(110.0));
        sg.terminal_growth_rate = Dual::variable(0.03, 0);
        let ad_g = solve_penman_damodaran_riv_aeg_and_implied_erp_ift(&sg).implied_erp;
        assert!(
            (ad_g.d[0] - fd_g).abs() < 1e-9
                && (ad_g.d[0] - r.d_ierp_d_terminal_growth).abs() < 1e-9
        );
    }
}
