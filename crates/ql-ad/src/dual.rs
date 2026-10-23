use crate::scalar::Real;
use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::cmp::Ordering;

fn erf0(x: f64) -> f64 {
    libm::erf(x)
}
fn cdf0(x: f64) -> f64 {
    0.5 * libm::erfc(-x * std::f64::consts::FRAC_1_SQRT_2)
}
fn pdf0(x: f64) -> f64 {
    (-0.5 * x * x).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

#[repr(C, align(32))]
#[derive(Copy, Clone, Debug)]
pub struct Dual<const N: usize> {
    pub val: f64,
    pub d: [f64; N],
}
impl<const N: usize> Default for Dual<N> {
    fn default() -> Self {
        Self {
            val: 0.0,
            d: [0.0; N],
        }
    }
}
impl<const N: usize> PartialEq for Dual<N> {
    /// Equality of an AD scalar is equality of its primal value.  This keeps
    /// `PartialEq` consistent with `PartialOrd`: derivative lanes are
    /// metadata, not part of the mathematical scalar being compared.
    fn eq(&self, other: &Self) -> bool {
        self.val == other.val
    }
}
impl<const N: usize> PartialOrd for Dual<N> {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        self.val.partial_cmp(&o.val)
    }
}
impl<const N: usize> fmt::Display for Dual<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.val)
    }
}
impl<const N: usize> Dual<N> {
    pub fn new(val: f64, index: usize) -> Self {
        Self::variable(val, index)
    }
    pub fn seed(val: f64, index: usize) -> Self {
        Self::variable(val, index)
    }
    pub fn variable(val: f64, index: usize) -> Self {
        let mut d = [0.0; N];
        if index < N {
            d[index] = 1.0;
        }
        Self { val, d }
    }
    pub fn from_value(val: f64) -> Self {
        Self { val, d: [0.0; N] }
    }
    #[inline]
    fn scale(&self, v: f64, s: f64) -> Self {
        let mut d = [0.0; N];
        for (out, &x) in d.iter_mut().zip(self.d.iter()) {
            // An inactive lane must remain inactive even when the scalar
            // derivative is infinite (e.g. 0^p at p < 1).  IEEE 754's
            // 0*inf is NaN, so explicitly preserve zero lanes.
            *out = if x == 0.0 { 0.0 } else { x * s };
        }
        Self { val: v, d }
    }
    #[inline]
    fn unary(self, g: f64, g1: f64) -> Self {
        self.scale(g, g1)
    }
}
impl<const N: usize> Add for Dual<N> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = self.d[i] + o.d[i];
        }
        Self {
            val: self.val + o.val,
            d,
        }
    }
}
impl<const N: usize> Sub for Dual<N> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = self.d[i] - o.d[i];
        }
        Self {
            val: self.val - o.val,
            d,
        }
    }
}
impl<const N: usize> Mul for Dual<N> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let mut d = [0.0; N];
        for i in 0..N {
            let lhs = if self.d[i] == 0.0 {
                0.0
            } else {
                self.d[i] * o.val
            };
            let rhs = if o.d[i] == 0.0 {
                0.0
            } else {
                self.val * o.d[i]
            };
            d[i] = lhs + rhs;
        }
        Self {
            val: self.val * o.val,
            d,
        }
    }
}
impl<const N: usize> Div for Dual<N> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let inv = 1.0 / o.val;
        let mut d = [0.0; N];
        for i in 0..N {
            let lhs = if self.d[i] == 0.0 {
                0.0
            } else {
                self.d[i] * o.val
            };
            let rhs = if o.d[i] == 0.0 {
                0.0
            } else {
                self.val * o.d[i]
            };
            d[i] = (lhs - rhs) * inv * inv;
        }
        Self {
            val: self.val * inv,
            d,
        }
    }
}
impl<const N: usize> Neg for Dual<N> {
    type Output = Self;
    fn neg(self) -> Self {
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = -self.d[i];
        }
        Self { val: -self.val, d }
    }
}
impl<const N: usize> Add<f64> for Dual<N> {
    type Output = Self;
    fn add(self, o: f64) -> Self {
        Self {
            val: self.val + o,
            d: self.d,
        }
    }
}
impl<const N: usize> Sub<f64> for Dual<N> {
    type Output = Self;
    fn sub(self, o: f64) -> Self {
        Self {
            val: self.val - o,
            d: self.d,
        }
    }
}
impl<const N: usize> Mul<f64> for Dual<N> {
    type Output = Self;
    fn mul(self, o: f64) -> Self {
        let mut d = self.d;
        for x in &mut d {
            if *x != 0.0 {
                *x *= o;
            }
        }
        Self {
            val: self.val * o,
            d,
        }
    }
}
impl<const N: usize> Div<f64> for Dual<N> {
    type Output = Self;
    fn div(self, o: f64) -> Self {
        self * (1.0 / o)
    }
}
impl<const N: usize> AddAssign for Dual<N> {
    fn add_assign(&mut self, o: Self) {
        *self = *self + o
    }
}
impl<const N: usize> SubAssign for Dual<N> {
    fn sub_assign(&mut self, o: Self) {
        *self = *self - o
    }
}
impl<const N: usize> MulAssign for Dual<N> {
    fn mul_assign(&mut self, o: Self) {
        *self = *self * o
    }
}
impl<const N: usize> DivAssign for Dual<N> {
    fn div_assign(&mut self, o: Self) {
        *self = *self / o
    }
}
impl<const N: usize> Real for Dual<N> {
    fn from_f64(v: f64) -> Self {
        Self::from_value(v)
    }
    fn val(self) -> f64 {
        self.val
    }
    fn abs(self) -> Self {
        let s = if self.val >= 0.0 { 1.0 } else { -1.0 };
        self.unary(self.val.abs(), s)
    }
    fn sqrt(self) -> Self {
        let g = self.val.sqrt();
        if g == 0.0 {
            // The classical derivative is singular at the origin.  Returning
            // the finite subgradient keeps AD expressions well-defined for
            // non-negative quantities that can legitimately hit zero.
            self.unary(g, 0.0)
        } else {
            self.unary(g, 0.5 / g)
        }
    }
    fn exp(self) -> Self {
        let g = self.val.exp();
        self.unary(g, g)
    }
    fn ln(self) -> Self {
        self.unary(self.val.ln(), 1.0 / self.val)
    }
    fn sin(self) -> Self {
        self.unary(self.val.sin(), self.val.cos())
    }
    fn cos(self) -> Self {
        self.unary(self.val.cos(), -self.val.sin())
    }
    fn tan(self) -> Self {
        let c = self.val.cos();
        self.unary(self.val.tan(), 1.0 / (c * c))
    }
    fn atan(self) -> Self {
        self.unary(self.val.atan(), 1.0 / (1.0 + self.val * self.val))
    }
    fn sinh(self) -> Self {
        self.unary(self.val.sinh(), self.val.cosh())
    }
    fn cosh(self) -> Self {
        self.unary(self.val.cosh(), self.val.sinh())
    }
    fn tanh(self) -> Self {
        let g = self.val.tanh();
        self.unary(g, 1.0 - g * g)
    }
    fn powf(self, p: f64) -> Self {
        if p == 0.0 {
            return Self::one();
        }
        if p == 1.0 {
            return self;
        }
        if p == 2.0 {
            return self * self;
        }
        let g = self.val.powf(p);
        let g1 = if self.val == 0.0 && p > 1.0 {
            0.0
        } else {
            p * self.val.powf(p - 1.0)
        };
        self.unary(g, g1)
    }
    fn pow(self, p: Self) -> Self {
        // A constant exponent can use the guarded scalar-power path.  In
        // particular this avoids evaluating ln(0) and 0.powf(-1) for the
        // common x^0/x^1/x^2 cases.
        if p.d.iter().all(|&x| x == 0.0) {
            return self.powf(p.val);
        }
        let g = self.val.powf(p.val);
        let ln_x = self.val.ln();
        let mut d = [0.0; N];
        for i in 0..N {
            let a = if p.d[i] == 0.0 { 0.0 } else { p.d[i] * ln_x };
            let b = if self.d[i] == 0.0 {
                0.0
            } else {
                p.val * self.d[i] / self.val
            };
            let q = a + b;
            d[i] = if q == 0.0 { 0.0 } else { g * q };
        }
        Self { val: g, d }
    }
    fn mul_add(self, a: Self, b: Self) -> Self {
        self * a + b
    }
    fn erf(self) -> Self {
        self.unary(
            erf0(self.val),
            2.0 / std::f64::consts::PI.sqrt() * (-self.val * self.val).exp(),
        )
    }
    fn erfc(self) -> Self {
        self.unary(
            libm::erfc(self.val),
            -2.0 / std::f64::consts::PI.sqrt() * (-self.val * self.val).exp(),
        )
    }
    fn norm_cdf(self) -> Self {
        self.unary(cdf0(self.val), pdf0(self.val))
    }
    fn norm_pdf(self) -> Self {
        self.unary(pdf0(self.val), -self.val * pdf0(self.val))
    }
    fn min_val(self, o: Self) -> Self {
        if self.val <= o.val { self } else { o }
    }
    fn max_val(self, o: Self) -> Self {
        if self.val >= o.val { self } else { o }
    }
}

#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct Dual2<const N: usize> {
    pub val: f64,
    pub d1: [f64; N],
    pub d2: [[f64; N]; N],
}
impl<const N: usize> Default for Dual2<N> {
    fn default() -> Self {
        Self {
            val: 0.0,
            d1: [0.0; N],
            d2: [[0.0; N]; N],
        }
    }
}
impl<const N: usize> PartialEq for Dual2<N> {
    fn eq(&self, other: &Self) -> bool {
        self.val == other.val
    }
}
impl<const N: usize> PartialOrd for Dual2<N> {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        self.val.partial_cmp(&o.val)
    }
}
impl<const N: usize> fmt::Display for Dual2<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.val)
    }
}
impl<const N: usize> Dual2<N> {
    pub fn new(val: f64, index: usize) -> Self {
        Self::variable(val, index)
    }
    pub fn seed(val: f64, index: usize) -> Self {
        Self::variable(val, index)
    }
    pub fn variable(val: f64, index: usize) -> Self {
        let mut x = Self {
            val,
            d1: [0.0; N],
            d2: [[0.0; N]; N],
        };
        if index < N {
            x.d1[index] = 1.0;
        }
        x
    }
    pub fn from_value(val: f64) -> Self {
        Self {
            val,
            d1: [0.0; N],
            d2: [[0.0; N]; N],
        }
    }
    fn unary(self, g: f64, g1: f64, g2: f64) -> Self {
        let mut d1 = [0.0; N];
        let mut d2 = [[0.0; N]; N];
        for i in 0..N {
            d1[i] = if self.d1[i] == 0.0 {
                0.0
            } else {
                g1 * self.d1[i]
            };
            for j in 0..N {
                let first = if self.d2[i][j] == 0.0 {
                    0.0
                } else {
                    g1 * self.d2[i][j]
                };
                let second = if self.d1[i] == 0.0 || self.d1[j] == 0.0 {
                    0.0
                } else {
                    g2 * self.d1[i] * self.d1[j]
                };
                d2[i][j] = first + second;
            }
        }
        Self { val: g, d1, d2 }
    }
}
impl<const N: usize> Add for Dual2<N> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        let mut d1 = [0.0; N];
        let mut d2 = [[0.0; N]; N];
        for i in 0..N {
            d1[i] = self.d1[i] + o.d1[i];
            for j in 0..N {
                d2[i][j] = self.d2[i][j] + o.d2[i][j];
            }
        }
        Self {
            val: self.val + o.val,
            d1,
            d2,
        }
    }
}
impl<const N: usize> Sub for Dual2<N> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self + (-o)
    }
}
impl<const N: usize> Mul for Dual2<N> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let mut d1 = [0.0; N];
        let mut d2 = [[0.0; N]; N];
        for i in 0..N {
            let lhs = if self.d1[i] == 0.0 {
                0.0
            } else {
                self.d1[i] * o.val
            };
            let rhs = if o.d1[i] == 0.0 {
                0.0
            } else {
                self.val * o.d1[i]
            };
            d1[i] = lhs + rhs;
            for j in 0..N {
                let a = if self.d2[i][j] == 0.0 {
                    0.0
                } else {
                    self.d2[i][j] * o.val
                };
                let b = if self.d1[i] == 0.0 || o.d1[j] == 0.0 {
                    0.0
                } else {
                    self.d1[i] * o.d1[j]
                };
                let c = if self.d1[j] == 0.0 || o.d1[i] == 0.0 {
                    0.0
                } else {
                    self.d1[j] * o.d1[i]
                };
                let d = if o.d2[i][j] == 0.0 {
                    0.0
                } else {
                    self.val * o.d2[i][j]
                };
                d2[i][j] = a + b + c + d;
            }
        }
        Self {
            val: self.val * o.val,
            d1,
            d2,
        }
    }
}
impl<const N: usize> Div for Dual2<N> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let inv = 1.0 / o.val;
        let inv2 = inv * inv;
        self * o.unary(inv, -inv2, 2.0 * inv2 * inv)
    }
}
impl<const N: usize> Neg for Dual2<N> {
    type Output = Self;
    fn neg(self) -> Self {
        let mut d1 = [0.0; N];
        let mut d2 = [[0.0; N]; N];
        for i in 0..N {
            d1[i] = -self.d1[i];
            for j in 0..N {
                d2[i][j] = -self.d2[i][j];
            }
        }
        Self {
            val: -self.val,
            d1,
            d2,
        }
    }
}
impl<const N: usize> Add<f64> for Dual2<N> {
    type Output = Self;
    fn add(self, o: f64) -> Self {
        Self {
            val: self.val + o,
            d1: self.d1,
            d2: self.d2,
        }
    }
}
impl<const N: usize> Sub<f64> for Dual2<N> {
    type Output = Self;
    fn sub(self, o: f64) -> Self {
        self + (-o)
    }
}
impl<const N: usize> Mul<f64> for Dual2<N> {
    type Output = Self;
    fn mul(self, o: f64) -> Self {
        let mut d1 = self.d1;
        let mut d2 = self.d2;
        for i in 0..N {
            if d1[i] != 0.0 {
                d1[i] *= o;
            }
            for x in &mut d2[i] {
                if *x != 0.0 {
                    *x *= o;
                }
            }
        }
        Self {
            val: self.val * o,
            d1,
            d2,
        }
    }
}
impl<const N: usize> Div<f64> for Dual2<N> {
    type Output = Self;
    fn div(self, o: f64) -> Self {
        self * (1.0 / o)
    }
}
impl<const N: usize> AddAssign for Dual2<N> {
    fn add_assign(&mut self, o: Self) {
        *self = *self + o
    }
}
impl<const N: usize> SubAssign for Dual2<N> {
    fn sub_assign(&mut self, o: Self) {
        *self = *self - o
    }
}
impl<const N: usize> MulAssign for Dual2<N> {
    fn mul_assign(&mut self, o: Self) {
        *self = *self * o
    }
}
impl<const N: usize> DivAssign for Dual2<N> {
    fn div_assign(&mut self, o: Self) {
        *self = *self / o
    }
}
impl<const N: usize> Real for Dual2<N> {
    fn from_f64(v: f64) -> Self {
        Self::from_value(v)
    }
    fn val(self) -> f64 {
        self.val
    }
    fn abs(self) -> Self {
        self.unary(
            self.val.abs(),
            if self.val >= 0.0 { 1.0 } else { -1.0 },
            0.0,
        )
    }
    fn sqrt(self) -> Self {
        let g = self.val.sqrt();
        if g == 0.0 {
            self.unary(g, 0.0, 0.0)
        } else {
            self.unary(g, 0.5 / g, -0.25 / (g * g * g))
        }
    }
    fn exp(self) -> Self {
        let g = self.val.exp();
        self.unary(g, g, g)
    }
    fn ln(self) -> Self {
        self.unary(self.val.ln(), 1.0 / self.val, -1.0 / (self.val * self.val))
    }
    fn sin(self) -> Self {
        self.unary(self.val.sin(), self.val.cos(), -self.val.sin())
    }
    fn cos(self) -> Self {
        self.unary(self.val.cos(), -self.val.sin(), -self.val.cos())
    }
    fn tan(self) -> Self {
        let c = self.val.cos();
        let g = 1.0 / (c * c);
        self.unary(self.val.tan(), g, 2.0 * self.val.tan() * g)
    }
    fn atan(self) -> Self {
        self.unary(
            self.val.atan(),
            1.0 / (1.0 + self.val * self.val),
            -2.0 * self.val / (1.0 + self.val * self.val).powi(2),
        )
    }
    fn sinh(self) -> Self {
        self.unary(self.val.sinh(), self.val.cosh(), self.val.sinh())
    }
    fn cosh(self) -> Self {
        self.unary(self.val.cosh(), self.val.sinh(), self.val.cosh())
    }
    fn tanh(self) -> Self {
        let g = self.val.tanh();
        self.unary(g, 1.0 - g * g, -2.0 * g * (1.0 - g * g))
    }
    fn powf(self, p: f64) -> Self {
        if p == 0.0 {
            return Self::one();
        }
        if p == 1.0 {
            return self;
        }
        if p == 2.0 {
            return self * self;
        }
        let g = self.val.powf(p);
        let g1 = if self.val == 0.0 && p > 1.0 {
            0.0
        } else {
            p * self.val.powf(p - 1.0)
        };
        let g2 = if self.val == 0.0 && p > 2.0 {
            0.0
        } else {
            p * (p - 1.0) * self.val.powf(p - 2.0)
        };
        self.unary(g, g1, g2)
    }
    fn pow(self, p: Self) -> Self {
        if p.d1.iter().all(|&x| x == 0.0) && p.d2.iter().flatten().all(|&x| x == 0.0) {
            return self.powf(p.val);
        }
        (p * self.ln()).exp()
    }
    fn mul_add(self, a: Self, b: Self) -> Self {
        self * a + b
    }
    fn erf(self) -> Self {
        let e = (-self.val * self.val).exp();
        self.unary(
            erf0(self.val),
            2.0 / std::f64::consts::PI.sqrt() * e,
            -4.0 * self.val / std::f64::consts::PI.sqrt() * e,
        )
    }
    fn erfc(self) -> Self {
        let e = (-self.val * self.val).exp();
        self.unary(
            libm::erfc(self.val),
            -2.0 / std::f64::consts::PI.sqrt() * e,
            4.0 * self.val / std::f64::consts::PI.sqrt() * e,
        )
    }
    fn norm_cdf(self) -> Self {
        self.unary(cdf0(self.val), pdf0(self.val), -self.val * pdf0(self.val))
    }
    fn norm_pdf(self) -> Self {
        let p = pdf0(self.val);
        self.unary(p, -self.val * p, (self.val * self.val - 1.0) * p)
    }
    fn min_val(self, o: Self) -> Self {
        if self.val <= o.val { self } else { o }
    }
    fn max_val(self, o: Self) -> Self {
        if self.val >= o.val { self } else { o }
    }
}

// Scalar-left convenience operators used in generic formulas.
impl<const N: usize> Add<Dual<N>> for f64 {
    type Output = Dual<N>;
    fn add(self, r: Dual<N>) -> Dual<N> {
        r + self
    }
}
impl<const N: usize> Sub<Dual<N>> for f64 {
    type Output = Dual<N>;
    fn sub(self, r: Dual<N>) -> Dual<N> {
        Dual::from_value(self) - r
    }
}
impl<const N: usize> Mul<Dual<N>> for f64 {
    type Output = Dual<N>;
    fn mul(self, r: Dual<N>) -> Dual<N> {
        r * self
    }
}
impl<const N: usize> Div<Dual<N>> for f64 {
    type Output = Dual<N>;
    fn div(self, r: Dual<N>) -> Dual<N> {
        Dual::from_value(self) / r
    }
}
impl<const N: usize> Add<Dual2<N>> for f64 {
    type Output = Dual2<N>;
    fn add(self, r: Dual2<N>) -> Dual2<N> {
        r + self
    }
}
impl<const N: usize> Sub<Dual2<N>> for f64 {
    type Output = Dual2<N>;
    fn sub(self, r: Dual2<N>) -> Dual2<N> {
        Dual2::from_value(self) - r
    }
}
impl<const N: usize> Mul<Dual2<N>> for f64 {
    type Output = Dual2<N>;
    fn mul(self, r: Dual2<N>) -> Dual2<N> {
        r * self
    }
}
impl<const N: usize> Div<Dual2<N>> for f64 {
    type Output = Dual2<N>;
    fn div(self, r: Dual2<N>) -> Dual2<N> {
        Dual2::from_value(self) / r
    }
}
