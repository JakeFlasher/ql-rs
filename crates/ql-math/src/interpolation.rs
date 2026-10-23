use crate::tridiagonal::solve_tridiagonal;
use ql_ad::Real;

pub trait Interpolator1D<T: Real>: Send + Sync {
    fn eval(&self, x: T) -> T;
    fn derivative(&self, x: T) -> T;
}
fn locate(xs: &[f64], z: f64) -> usize {
    if xs.len() < 2 || z <= xs[0] {
        0
    } else if z >= xs[xs.len() - 1] {
        xs.len() - 2
    } else {
        xs.partition_point(|v| *v <= z)
            .saturating_sub(1)
            .min(xs.len() - 2)
    }
}
#[derive(Debug, Clone)]
pub struct LinearInterpolator<T: Real> {
    pub x: Vec<f64>,
    pub y: Vec<T>,
    pub flat_extrapolation: bool,
}
impl<T: Real> LinearInterpolator<T> {
    pub fn new(x: Vec<f64>, y: Vec<T>) -> Self {
        assert!(x.len() >= 2 && x.len() == y.len());
        Self {
            x,
            y,
            flat_extrapolation: false,
        }
    }
    pub fn with_flat_extrapolation(mut self, v: bool) -> Self {
        self.flat_extrapolation = v;
        self
    }
}
impl<T: Real> Interpolator1D<T> for LinearInterpolator<T> {
    fn eval(&self, z: T) -> T {
        let i = locate(&self.x, z.val());
        if self.flat_extrapolation && (z.val() < self.x[0] || z.val() > self.x[self.x.len() - 1]) {
            return if z.val() < self.x[0] {
                self.y[0]
            } else {
                *self.y.last().unwrap()
            };
        };
        let h = self.x[i + 1] - self.x[i];
        self.y[i] + (z - T::from_f64(self.x[i])) * (self.y[i + 1] - self.y[i]) / h
    }
    fn derivative(&self, z: T) -> T {
        let i = locate(&self.x, z.val());
        if self.flat_extrapolation && (z.val() < self.x[0] || z.val() > self.x[self.x.len() - 1]) {
            T::zero()
        } else {
            (self.y[i + 1] - self.y[i]) / (self.x[i + 1] - self.x[i])
        }
    }
}
#[derive(Debug, Clone)]
pub struct LogLinearInterpolator<T: Real> {
    pub inner: LinearInterpolator<T>,
}
impl<T: Real> LogLinearInterpolator<T> {
    pub fn new(x: Vec<f64>, y: Vec<T>) -> Self {
        Self {
            inner: LinearInterpolator::new(x, y.into_iter().map(|v| v.ln()).collect()),
        }
    }
}
impl<T: Real> Interpolator1D<T> for LogLinearInterpolator<T> {
    fn eval(&self, z: T) -> T {
        self.inner.eval(z).exp()
    }
    fn derivative(&self, z: T) -> T {
        let e = self.inner.eval(z).exp();
        e * self.inner.derivative(z)
    }
}
#[derive(Debug, Clone)]
pub struct FlatForwardInterpolator<T: Real> {
    pub x: Vec<f64>,
    pub log_y: Vec<T>,
}
impl<T: Real> FlatForwardInterpolator<T> {
    pub fn new(x: Vec<f64>, y: Vec<T>) -> Self {
        Self {
            x,
            log_y: y.into_iter().map(|v| v.ln()).collect(),
        }
    }
}
impl<T: Real> Interpolator1D<T> for FlatForwardInterpolator<T> {
    fn eval(&self, z: T) -> T {
        let i = locate(&self.x, z.val());
        let slope = (self.log_y[i + 1] - self.log_y[i]) / (self.x[i + 1] - self.x[i]);
        (self.log_y[i] + (z - T::from_f64(self.x[i])) * slope).exp()
    }
    fn derivative(&self, z: T) -> T {
        let i = locate(&self.x, z.val());
        let slope = (self.log_y[i + 1] - self.log_y[i]) / (self.x[i + 1] - self.x[i]);
        slope * (self.log_y[i] + (z - T::from_f64(self.x[i])) * slope).exp()
    }
}
#[derive(Debug, Clone)]
pub struct MonotoneCubicInterpolator<T: Real> {
    pub x: Vec<f64>,
    pub y: Vec<T>,
    pub slopes: Vec<T>,
}
impl<T: Real> MonotoneCubicInterpolator<T> {
    pub fn new(x: Vec<f64>, y: Vec<T>) -> Self {
        let n = x.len();
        assert!(n >= 2 && n == y.len());
        let mut sec = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            sec.push((y[i + 1] - y[i]) / (x[i + 1] - x[i]));
        }
        let mut m = vec![T::zero(); n];
        m[0] = sec[0];
        m[n - 1] = sec[n - 2];
        for i in 1..n - 1 {
            m[i] = if sec[i - 1].val() * sec[i].val() <= 0.0 {
                T::zero()
            } else {
                (sec[i - 1] + sec[i]) * 0.5
            };
        }
        for i in 0..n - 1 {
            if sec[i].val() == 0.0 {
                m[i] = T::zero();
                m[i + 1] = T::zero();
            } else {
                let a = m[i] / sec[i];
                let b = m[i + 1] / sec[i];
                let s = a * a + b * b;
                if s.val() > 9.0 {
                    let t = T::from_f64(3.0) / s.sqrt();
                    m[i] = t * a * sec[i];
                    m[i + 1] = t * b * sec[i];
                }
            }
        }
        Self { x, y, slopes: m }
    }
}
impl<T: Real> Interpolator1D<T> for MonotoneCubicInterpolator<T> {
    fn eval(&self, z: T) -> T {
        let i = locate(&self.x, z.val());
        let h = self.x[i + 1] - self.x[i];
        let t = (z - T::from_f64(self.x[i])) / h;
        let t2 = t * t;
        let t3 = t2 * t;
        let h00 = T::from_f64(2.0) * t3 - T::from_f64(3.0) * t2 + T::one();
        let h10 = t3 - T::from_f64(2.0) * t2 + t;
        let h01 = -T::from_f64(2.0) * t3 + T::from_f64(3.0) * t2;
        let h11 = t3 - t2;
        h00 * self.y[i]
            + h10 * (self.slopes[i] * h)
            + h01 * self.y[i + 1]
            + h11 * (self.slopes[i + 1] * h)
    }
    fn derivative(&self, z: T) -> T {
        let i = locate(&self.x, z.val());
        let h = self.x[i + 1] - self.x[i];
        let t = (z - T::from_f64(self.x[i])) / h;
        let t2 = t * t;
        let a = (T::from_f64(6.0) * t2 - T::from_f64(6.0) * t) * self.y[i]
            + (T::from_f64(3.0) * t2 - T::from_f64(4.0) * t + T::one()) * (self.slopes[i] * h)
            + (-T::from_f64(6.0) * t2 + T::from_f64(6.0) * t) * self.y[i + 1]
            + (T::from_f64(3.0) * t2 - T::from_f64(2.0) * t) * (self.slopes[i + 1] * h);
        a / h
    }
}
#[derive(Debug, Clone)]
pub struct NaturalCubicSpline<T: Real> {
    pub x: Vec<f64>,
    pub y: Vec<T>,
    pub second: Vec<T>,
}
impl<T: Real> NaturalCubicSpline<T> {
    pub fn new(x: Vec<f64>, y: Vec<T>) -> Self {
        let n = x.len();
        assert!(n >= 2 && n == y.len());
        let mut second = vec![T::zero(); n];
        if n > 2 {
            let mut lo = vec![T::zero(); n - 2];
            let mut di = vec![T::zero(); n - 2];
            let mut up = vec![T::zero(); n - 2];
            let mut rhs = vec![T::zero(); n - 2];
            for i in 1..n - 1 {
                let h0 = x[i] - x[i - 1];
                let h1 = x[i + 1] - x[i];
                let j = i - 1;
                if j > 0 {
                    lo[j] = T::from_f64(h0)
                }
                di[j] = T::from_f64(2.0 * (h0 + h1));
                if j + 1 < n - 2 {
                    up[j] = T::from_f64(h1)
                }
                rhs[j] = T::from_f64(6.0) * ((y[i + 1] - y[i]) / h1 - (y[i] - y[i - 1]) / h0);
            }
            let mut out = vec![T::zero(); n - 2];
            let mut scratch = vec![T::zero(); n - 2];
            solve_tridiagonal(&lo, &di, &up, &rhs, &mut out, &mut scratch);
            for (j, v) in out.into_iter().enumerate() {
                second[j + 1] = v;
            }
        }
        Self { x, y, second }
    }
}
impl<T: Real> Interpolator1D<T> for NaturalCubicSpline<T> {
    fn eval(&self, z: T) -> T {
        let i = locate(&self.x, z.val());
        let h = self.x[i + 1] - self.x[i];
        let a = (T::from_f64(self.x[i + 1]) - z) / h;
        let b = (z - T::from_f64(self.x[i])) / h;
        a * self.y[i]
            + b * self.y[i + 1]
            + ((a * a * a - a) * self.second[i] + (b * b * b - b) * self.second[i + 1])
                * (h * h / 6.0)
    }
    fn derivative(&self, z: T) -> T {
        let i = locate(&self.x, z.val());
        let h = self.x[i + 1] - self.x[i];
        let a = (T::from_f64(self.x[i + 1]) - z) / h;
        let b = (z - T::from_f64(self.x[i])) / h;
        (self.y[i + 1] - self.y[i]) / h
            + ((-T::from_f64(3.0) * a * a + T::one()) * self.second[i]
                + (T::from_f64(3.0) * b * b - T::one()) * self.second[i + 1])
                * h
                / 6.0
    }
}
