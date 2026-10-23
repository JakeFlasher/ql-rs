use ql_ad::Real;
use std::ops::{Add, Div, Mul, Neg, Sub};
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Complex<T: Real> {
    pub re: T,
    pub im: T,
}
impl<T: Real> Complex<T> {
    pub fn new(re: T, im: T) -> Self {
        Self { re, im }
    }
    pub fn from_real(re: T) -> Self {
        Self { re, im: T::zero() }
    }
    pub fn from_f64(re: f64, im: f64) -> Self {
        Self {
            re: T::from_f64(re),
            im: T::from_f64(im),
        }
    }
    pub fn i() -> Self {
        Self::from_f64(0., 1.)
    }
    pub fn zero() -> Self {
        Self::from_real(T::zero())
    }
    pub fn one() -> Self {
        Self::from_real(T::one())
    }
    pub fn norm_sqr(self) -> T {
        self.re * self.re + self.im * self.im
    }
    pub fn abs(self) -> T {
        self.norm_sqr().sqrt()
    }
    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }
    pub fn exp(self) -> Self {
        let r = self.re.exp();
        Self::new(r * self.im.cos(), r * self.im.sin())
    }
    pub fn sqrt(self) -> Self {
        let r = self.abs();
        if self.re.val() == 0. && self.im.val() == 0. {
            return Self::zero();
        }
        if self.re.val() >= 0. {
            let u = ((r + self.re) * 0.5).sqrt();
            Self::new(u, self.im / (u * 2.))
        } else {
            let sign = if self.im.val() >= 0. { 1. } else { -1. };
            let v = ((r - self.re) * 0.5).sqrt() * sign;
            Self::new(self.im / (v * 2.), v)
        }
    }
    pub fn arg(self) -> T {
        if self.re.val() > 0. {
            (self.im / self.re).atan()
        } else if self.re.val() < 0. && self.im.val() >= 0. {
            (self.im / self.re).atan() + std::f64::consts::PI
        } else if self.re.val() < 0. {
            (self.im / self.re).atan() - std::f64::consts::PI
        } else if self.im.val() > 0. {
            // Smooth atan2 limit on the imaginary axis.  `pi/2` by itself
            // drops the derivative with respect to a perturbation in Re.
            T::from_f64(std::f64::consts::FRAC_PI_2) - self.re / self.im
        } else {
            T::from_f64(-std::f64::consts::FRAC_PI_2) - self.re / self.im
        }
    }
    pub fn ln(self) -> Self {
        Self::new(self.abs().ln(), self.arg())
    }
}
impl<T: Real> Add for Complex<T> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }
}
impl<T: Real> Sub for Complex<T> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }
}
impl<T: Real> Mul for Complex<T> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}
impl<T: Real> Div for Complex<T> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let d = o.norm_sqr();
        Self::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
}
impl<T: Real> Neg for Complex<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}
impl<T: Real> Add<T> for Complex<T> {
    type Output = Self;
    fn add(self, x: T) -> Self {
        Self::new(self.re + x, self.im)
    }
}
impl<T: Real> Sub<T> for Complex<T> {
    type Output = Self;
    fn sub(self, x: T) -> Self {
        Self::new(self.re - x, self.im)
    }
}
impl<T: Real> Mul<T> for Complex<T> {
    type Output = Self;
    fn mul(self, x: T) -> Self {
        Self::new(self.re * x, self.im * x)
    }
}
impl<T: Real> Div<T> for Complex<T> {
    type Output = Self;
    fn div(self, x: T) -> Self {
        Self::new(self.re / x, self.im / x)
    }
}
