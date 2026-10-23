use nalgebra::{DMatrix, DVector};
use ql_ad::{Real, ift_lift_1d};
use thiserror::Error;
#[derive(Debug, Error, Clone, Copy)]
pub enum MathError {
    #[error("root is not bracketed")]
    NotBracketed,
    #[error("solver failed to converge")]
    MaxIterations,
    #[error("singular Jacobian")]
    Singular,
}
pub fn brent_solve<F: Fn(f64) -> f64>(
    f: F,
    mut a: f64,
    mut b: f64,
    tol: f64,
    max_iter: usize,
) -> Result<f64, MathError> {
    let mut fa = f(a);
    let mut fb = f(b);
    if fa == 0.0 {
        return Ok(a);
    }
    if fb == 0.0 {
        return Ok(b);
    }
    if fa * fb > 0.0 {
        return Err(MathError::NotBracketed);
    };
    let mut c = b;
    let mut fc = fb;
    let mut d = b - a;
    let mut e = d;
    for _ in 0..max_iter {
        if fb.signum() == fc.signum() {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol1 = 2.0 * f64::EPSILON * b.abs() + 0.5 * tol;
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol1 || fb == 0.0 {
            return Ok(b);
        }
        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (p0, q0) = if a == c {
                (2.0 * xm * s, 1.0 - s)
            } else {
                let q = fa / fc;
                let r = fb / fc;
                (
                    s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0)),
                    (q - 1.0) * (r - 1.0) * (s - 1.0),
                )
            };
            let (p, q) = if p0 > 0.0 { (p0, -q0) } else { (-p0, q0) };
            if 2.0 * p < (3.0 * xm * q - (tol1 * q).abs()).min((e * q).abs()) {
                e = d;
                d = p / q
            } else {
                d = xm;
                e = d
            }
        } else {
            d = xm;
            e = d
        }
        a = b;
        fa = fb;
        if d.abs() > tol1 {
            b += d
        } else {
            b += tol1.copysign(xm)
        }
        fb = f(b);
    }
    Err(MathError::MaxIterations)
}
pub fn newton_solve<F: Fn(f64) -> (f64, f64)>(
    f: F,
    mut x: f64,
    tol: f64,
    max_iter: usize,
) -> Result<f64, MathError> {
    for _ in 0..max_iter {
        let (v, d) = f(x);
        if v.abs() < tol {
            return Ok(x);
        };
        if d.abs() < f64::EPSILON {
            return Err(MathError::Singular);
        };
        let n = x - v / d;
        if (n - x).abs() < tol {
            return Ok(n);
        };
        x = n;
    }
    Err(MathError::MaxIterations)
}
pub fn halley_solve<F: Fn(f64) -> (f64, f64, f64)>(
    f: F,
    mut x: f64,
    tol: f64,
    max_iter: usize,
) -> Result<f64, MathError> {
    for _ in 0..max_iter {
        let (v, d, dd) = f(x);
        if v.abs() < tol {
            return Ok(x);
        }
        let den = 2.0 * d * d - v * dd;
        if den.abs() < f64::EPSILON {
            return Err(MathError::Singular);
        };
        let n = x - 2.0 * v * d / den;
        if (n - x).abs() < tol {
            return Ok(n);
        }
        x = n;
    }
    Err(MathError::MaxIterations)
}
pub fn brent_solve_ad<T: Real, F: Fn(T) -> T>(
    f: F,
    a: f64,
    b: f64,
    tol: f64,
    max_iter: usize,
) -> Result<T, MathError> {
    let root = brent_solve(|x| f(T::from_f64(x)).val(), a, b, tol, max_iter)?;
    let h = (1e-5_f64).max(root.abs() * 1e-5);
    let fp_ad = f(T::from_f64(root + h));
    let fm_ad = f(T::from_f64(root - h));
    let f0_ad = f(T::from_f64(root));
    let fx = (fp_ad - fm_ad) / (2.0 * h);
    let fxx = (fp_ad.val() - 2.0 * f0_ad.val() + fm_ad.val()) / (h * h);
    Ok(ift_lift_1d(root, f, fx, fxx))
}
pub fn levenberg_marquardt<F>(
    mut theta: DVector<f64>,
    residual: F,
    tol: f64,
    max_iter: usize,
) -> Result<DVector<f64>, MathError>
where
    F: Fn(&DVector<f64>) -> (DVector<f64>, DMatrix<f64>),
{
    let mut lambda: f64 = 1e-3;
    for _ in 0..max_iter {
        let (r, j) = residual(&theta);
        let cost = r.norm_squared();
        if cost < tol * tol {
            return Ok(theta.clone());
        }
        // A rejected step must be retried with stronger damping.  The old
        // implementation increased lambda based on the *previous* iterate,
        // then silently discarded a rejected step, which is not a
        // Levenberg--Marquardt trust-region update and can stall calibration.
        let jt = j.transpose();
        let rhs = -&jt * &r;
        let mut accepted = false;
        for _ in 0..32 {
            let a = &jt * &j + DMatrix::identity(theta.len(), theta.len()) * lambda;
            // QR is more robust than an LU solve of the normal equations near
            // rank-deficient calibration Jacobians.  Fall back to SVD when QR
            // detects a singular factorization.
            let d = a
                .clone()
                .qr()
                .solve(&rhs)
                .or_else(|| a.svd(true, true).solve(&rhs, 1e-12).ok())
                .ok_or(MathError::Singular)?;
            if d.norm() < tol {
                return Ok(theta);
            }
            let trial = &theta + &d;
            let (rt, _) = residual(&trial);
            let trial_cost = rt.norm_squared();
            if trial_cost < cost {
                theta = trial;
                lambda = (lambda * 0.3).max(1e-12);
                accepted = true;
                break;
            }
            lambda = (lambda * 10.0).min(1e16);
        }
        if !accepted && lambda >= 1e16 {
            return Err(MathError::MaxIterations);
        }
    }
    Err(MathError::MaxIterations)
}
