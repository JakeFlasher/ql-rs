use core::fmt::{Debug, Display};
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use core::simd::Select;
use core::simd::Simd;
use core::simd::prelude::{SimdFloat, SimdInt, SimdPartialOrd};
use std::simd::StdFloat;

pub trait Real:
    Copy
    + Clone
    + Default
    + Send
    + Sync
    + PartialEq
    + PartialOrd
    + Debug
    + Display
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Div<Self, Output = Self>
    + Add<f64, Output = Self>
    + Sub<f64, Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
{
    fn from_f64(v: f64) -> Self;
    fn val(self) -> f64;
    fn to_f64(self) -> f64 {
        self.val()
    }
    fn zero() -> Self {
        Self::from_f64(0.0)
    }
    fn one() -> Self {
        Self::from_f64(1.0)
    }
    fn abs(self) -> Self;
    fn sqrt(self) -> Self;
    fn exp(self) -> Self;
    fn ln(self) -> Self;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn tan(self) -> Self;
    fn atan(self) -> Self;
    fn sinh(self) -> Self;
    fn cosh(self) -> Self;
    fn tanh(self) -> Self;
    fn powf(self, p: f64) -> Self;
    fn pow(self, p: Self) -> Self;
    fn mul_add(self, a: Self, b: Self) -> Self;
    fn erf(self) -> Self;
    fn erfc(self) -> Self {
        Self::one() - self.erf()
    }
    fn norm_cdf(self) -> Self;
    fn norm_pdf(self) -> Self;
    fn min_val(self, other: Self) -> Self;
    fn max_val(self, other: Self) -> Self;
    fn select_ge_zero(self, if_true: Self, if_false: Self) -> Self {
        if self.val() >= 0.0 { if_true } else { if_false }
    }
    fn smooth_heaviside(self, threshold: Self, half_width: f64) -> Self {
        let diff = self - threshold;
        if half_width <= 0.0 {
            return if diff.val() >= 0.0 {
                Self::one()
            } else {
                Self::zero()
            };
        }
        let z = diff.val();
        if z <= -half_width {
            Self::zero()
        } else if z >= half_width {
            Self::one()
        } else {
            // Keep the polynomial in the AD domain; converting through `val`
            // would silently discard all derivatives in the transition band.
            let u = diff / half_width;
            Self::from_f64(0.5) + u * 0.75 - u * u * u * 0.25
        }
    }
}

pub trait SimdReal: Copy + Clone + Send + Sync + Debug {
    fn exp(self) -> Self;
    fn ln(self) -> Self;
    fn sqrt(self) -> Self;
    fn norm_cdf(self) -> Self;
    fn norm_pdf(self) -> Self;
    fn select_ge_zero(self, t: Self, f: Self) -> Self;
}
pub type F64x4 = Simd<f64, 4>;
pub type F64x8 = Simd<f64, 8>;

fn norm_cdf_f64(x: f64) -> f64 {
    0.5 * (1.0 + libm::erf(x * std::f64::consts::FRAC_1_SQRT_2))
}
fn norm_pdf_f64(x: f64) -> f64 {
    std::f64::consts::TAU.sqrt().recip() * (-0.5 * x * x).exp()
}

macro_rules! impl_real_f64 {
    () => {
        impl Real for f64 {
            fn from_f64(v: f64) -> Self {
                v
            }
            fn val(self) -> f64 {
                self
            }
            fn abs(self) -> Self {
                self.abs()
            }
            fn sqrt(self) -> Self {
                self.sqrt()
            }
            fn exp(self) -> Self {
                self.exp()
            }
            fn ln(self) -> Self {
                self.ln()
            }
            fn sin(self) -> Self {
                self.sin()
            }
            fn cos(self) -> Self {
                self.cos()
            }
            fn tan(self) -> Self {
                self.tan()
            }
            fn atan(self) -> Self {
                self.atan()
            }
            fn sinh(self) -> Self {
                self.sinh()
            }
            fn cosh(self) -> Self {
                self.cosh()
            }
            fn tanh(self) -> Self {
                self.tanh()
            }
            fn powf(self, p: f64) -> Self {
                self.powf(p)
            }
            fn pow(self, p: Self) -> Self {
                self.powf(p)
            }
            fn mul_add(self, a: Self, b: Self) -> Self {
                self.mul_add(a, b)
            }
            fn erf(self) -> Self {
                libm::erf(self)
            }
            fn erfc(self) -> Self {
                libm::erfc(self)
            }
            fn norm_cdf(self) -> Self {
                norm_cdf_f64(self)
            }
            fn norm_pdf(self) -> Self {
                norm_pdf_f64(self)
            }
            fn min_val(self, o: Self) -> Self {
                self.min(o)
            }
            fn max_val(self, o: Self) -> Self {
                self.max(o)
            }
        }
    };
}
impl_real_f64!();

#[inline(always)]
fn simd_exp<const N: usize>(x: Simd<f64, N>) -> Simd<f64, N> {
    // Cody--Waite range reduction followed by Estrin evaluation of the degree-11
    // Taylor polynomial on [-ln(2)/2, ln(2)/2] (max relative error < 6.2e-15).
    // Uses IEEE-754 magic-number bit manipulation (`3 * 2^51`) to construct 2^k
    // without `f64 -> i64` casts (which scalarize on baseline SSE2/AVX2 targets)
    // and avoids `StdFloat::mul_add` (which lowers to scalar `libm::fma` calls
    // when compiled without `-C target-feature=+fma`).
    const LOG2_E: f64 = std::f64::consts::LOG2_E;
    const MAGIC: f64 = 6_755_399_441_055_744.0; // 3 * 2^51
    let x = x.simd_clamp(Simd::splat(-708.0), Simd::splat(709.0));
    let u = x * Simd::splat(LOG2_E) + Simd::splat(MAGIC);
    let k = u - Simd::splat(MAGIC);
    let r = (x - k * Simd::splat(0.6931471803691238)) - k * Simd::splat(1.9082149292705877e-10);

    let r2 = r * r;
    let r4 = r2 * r2;
    let p01 = Simd::splat(1.0) + r;
    let p23 = Simd::splat(0.5) + r * Simd::splat(1.0 / 6.0);
    let p45 = Simd::splat(1.0 / 24.0) + r * Simd::splat(1.0 / 120.0);
    let p67 = Simd::splat(1.0 / 720.0) + r * Simd::splat(1.0 / 5_040.0);
    let p89 = Simd::splat(1.0 / 40_320.0) + r * Simd::splat(1.0 / 362_880.0);
    let p10_11 = Simd::splat(1.0 / 3_628_800.0) + r * Simd::splat(1.0 / 39_916_800.0);
    let p03 = p01 + r2 * p23;
    let p47 = p45 + r2 * p67;
    let p8_11 = p89 + r2 * p10_11;
    let p = p03 + r4 * (p47 + r4 * p8_11);

    let k_i: Simd<i64, N> = k.cast();
    let biased: Simd<u64, N> = (k_i + Simd::splat(1023_i64)).cast();
    let exp_bits = biased << Simd::<u64, N>::splat(52_u64);
    let pow2k: Simd<f64, N> = Simd::from_bits(exp_bits);
    p * pow2k
}

#[inline(always)]
fn simd_exp4(x: Simd<f64, 4>) -> Simd<f64, 4> {
    simd_exp(x)
}

#[inline(always)]
fn simd_exp8(x: Simd<f64, 8>) -> Simd<f64, 8> {
    simd_exp(x)
}

fn simd_erf4(x: Simd<f64, 4>) -> Simd<f64, 4> {
    let sign = x.simd_lt(Simd::splat(0.0));
    let a = x.abs();
    let t = Simd::splat(1.0) / (Simd::splat(1.0) + Simd::splat(0.3275911) * a);
    let p = ((((Simd::splat(1.061405429) * t - Simd::splat(1.453152027)) * t
        + Simd::splat(1.421413741))
        * t
        - Simd::splat(0.284496736))
        * t
        + Simd::splat(0.254829592))
        * t;
    let y = Simd::splat(1.0) - p * simd_exp4(-a * a);
    sign.select(-y, y)
}
fn simd_erf8(x: Simd<f64, 8>) -> Simd<f64, 8> {
    let sign = x.simd_lt(Simd::splat(0.0));
    let a = x.abs();
    let t = Simd::splat(1.0) / (Simd::splat(1.0) + Simd::splat(0.3275911) * a);
    let p = ((((Simd::splat(1.061405429) * t - Simd::splat(1.453152027)) * t
        + Simd::splat(1.421413741))
        * t
        - Simd::splat(0.284496736))
        * t
        + Simd::splat(0.254829592))
        * t;
    let y = Simd::splat(1.0) - p * simd_exp8(-a * a);
    sign.select(-y, y)
}
impl SimdReal for Simd<f64, 4> {
    #[inline(always)]
    fn exp(self) -> Self {
        simd_exp4(self)
    }
    fn ln(self) -> Self {
        StdFloat::ln(self)
    }
    fn sqrt(self) -> Self {
        StdFloat::sqrt(self)
    }
    fn norm_cdf(self) -> Self {
        Simd::splat(0.5)
            * (Simd::splat(1.0) + simd_erf4(self * Simd::splat(std::f64::consts::FRAC_1_SQRT_2)))
    }
    fn norm_pdf(self) -> Self {
        Simd::splat((2.0 * std::f64::consts::PI).sqrt().recip())
            * simd_exp4(Simd::splat(-0.5) * self * self)
    }
    fn select_ge_zero(self, t: Self, f: Self) -> Self {
        self.simd_ge(Simd::splat(0.0)).select(t, f)
    }
}
impl SimdReal for Simd<f64, 8> {
    #[inline(always)]
    fn exp(self) -> Self {
        simd_exp8(self)
    }
    fn ln(self) -> Self {
        StdFloat::ln(self)
    }
    fn sqrt(self) -> Self {
        StdFloat::sqrt(self)
    }
    fn norm_cdf(self) -> Self {
        Simd::splat(0.5)
            * (Simd::splat(1.0) + simd_erf8(self * Simd::splat(std::f64::consts::FRAC_1_SQRT_2)))
    }
    fn norm_pdf(self) -> Self {
        Simd::splat((2.0 * std::f64::consts::PI).sqrt().recip())
            * simd_exp8(Simd::splat(-0.5) * self * self)
    }
    fn select_ge_zero(self, t: Self, f: Self) -> Self {
        self.simd_ge(Simd::splat(0.0)).select(t, f)
    }
}
