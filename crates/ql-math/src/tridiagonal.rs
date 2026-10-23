use core::simd::Simd;
use ql_ad::Real;
use std::marker::PhantomData;
pub fn solve_tridiagonal<T: Real>(
    lower: &[T],
    diag: &[T],
    upper: &[T],
    rhs: &[T],
    out: &mut [T],
    scratch: &mut [T],
) {
    let n = diag.len();
    assert!(
        n > 0
            && lower.len() >= n.saturating_sub(1)
            && upper.len() >= n.saturating_sub(1)
            && rhs.len() >= n
            && out.len() >= n
            && scratch.len() >= n
    );
    scratch[0] = if n > 1 { upper[0] / diag[0] } else { T::zero() };
    out[0] = rhs[0] / diag[0];
    for i in 1..n {
        let den = diag[i] - lower[i - 1] * scratch[i - 1];
        scratch[i] = if i + 1 < n { upper[i] / den } else { T::zero() };
        out[i] = (rhs[i] - lower[i - 1] * out[i - 1]) / den;
    }
    for i in (0..n.saturating_sub(1)).rev() {
        out[i] -= scratch[i] * out[i + 1];
    }
}
#[derive(Debug, Clone, Copy)]
pub struct SimdThomasSolver<const LANES: usize> {
    pub n: usize,
    _p: PhantomData<[f64; LANES]>,
}
impl<const LANES: usize> SimdThomasSolver<LANES> {
    pub const fn new(n: usize) -> Self {
        Self { n, _p: PhantomData }
    }
}
macro_rules! impl_simd {
    ($lanes:literal) => {
        impl SimdThomasSolver<$lanes> {
            pub fn solve_in_place(
                &self,
                lower: &[Simd<f64, $lanes>],
                diag: &[Simd<f64, $lanes>],
                upper: &[Simd<f64, $lanes>],
                rhs: &[Simd<f64, $lanes>],
                out: &mut [Simd<f64, $lanes>],
                scratch: &mut [Simd<f64, $lanes>],
            ) {
                let n = self.n;
                assert!(
                    n > 0
                        && lower.len() >= n.saturating_sub(1)
                        && upper.len() >= n.saturating_sub(1)
                        && diag.len() >= n
                        && rhs.len() >= n
                        && out.len() >= n
                        && scratch.len() >= n
                );
                scratch[0] = if n > 1 {
                    upper[0] / diag[0]
                } else {
                    Simd::splat(0.0)
                };
                out[0] = rhs[0] / diag[0];
                for i in 1..n {
                    let den = diag[i] - lower[i - 1] * scratch[i - 1];
                    scratch[i] = if i + 1 < n {
                        upper[i] / den
                    } else {
                        Simd::splat(0.0)
                    };
                    out[i] = (rhs[i] - lower[i - 1] * out[i - 1]) / den;
                }
                for i in (0..n.saturating_sub(1)).rev() {
                    out[i] -= scratch[i] * out[i + 1];
                }
            }
        }
    };
}
impl_simd!(4);
impl_simd!(8);
