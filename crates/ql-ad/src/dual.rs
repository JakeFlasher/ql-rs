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
#[derive(Copy, Clone, Debug, PartialEq)]
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
impl<const N: usize> Eq for Dual<N> {}
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
    fn unary(self, g: f64, g1: f64) -> Self {
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = g1 * self.d[i];
        }
        Self { val: g, d }
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
            d[i] = self.d[i] * o.val + self.val * o.d[i];
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
            d[i] = (self.d[i] * o.val - self.val * o.d[i]) * inv * inv;
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
            *x *= o;
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
        let g = self.val.powf(p);
        self.unary(g, p * self.val.powf(p - 1.0))
    }
    fn pow(self, p: Self) -> Self {
        let g = self.val.powf(p.val);
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = g * (p.d[i] * self.val.ln() + p.val * self.d[i] / self.val);
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
#[derive(Copy, Clone, Debug, PartialEq)]
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
impl<const N: usize> Eq for Dual2<N> {}
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
            d1[i] = g1 * self.d1[i];
            for j in 0..N {
                d2[i][j] = g1 * self.d2[i][j] + g2 * self.d1[i] * self.d1[j];
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
            d1[i] = self.d1[i] * o.val + self.val * o.d1[i];
            for j in 0..N {
                d2[i][j] = self.d2[i][j] * o.val
                    + self.d1[i] * o.d1[j]
                    + self.d1[j] * o.d1[i]
                    + self.val * o.d2[i][j];
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
        self * o.powf(-1.0)
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
            d1[i] *= o;
            for x in &mut d2[i] {
                *x *= o;
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
        let g = self.val.powf(p);
        self.unary(
            g,
            p * self.val.powf(p - 1.0),
            p * (p - 1.0) * self.val.powf(p - 2.0),
        )
    }
    fn pow(self, p: Self) -> Self {
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
