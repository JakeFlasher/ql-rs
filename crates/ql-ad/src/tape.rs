use crate::scalar::Real;
use std::cell::RefCell;
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::atomic::{AtomicU16, Ordering};

pub const CONST_IDX: u32 = u32::MAX;
static NEXT_TAPE_ID: AtomicU16 = AtomicU16::new(1);
#[derive(Copy, Clone, Debug, Default)]
pub struct TapeVar<T: Copy = f64> {
    pub val: T,
    pub idx: u32,
    pub tape_id: u16,
}
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct OpRecord<T: Copy> {
    pub lhs: u32,
    pub rhs: u32,
    pub w_lhs: T,
    pub w_rhs: T,
}
#[derive(Debug, Clone)]
pub struct AdjointResult {
    pub adjoints: Vec<f64>,
    pub output: u32,
}
impl AdjointResult {
    pub fn gradient(&self, var: TapeVar<f64>) -> f64 {
        if var.idx == CONST_IDX {
            0.0
        } else {
            self.adjoints.get(var.idx as usize).copied().unwrap_or(0.0)
        }
    }
    pub fn get(&self, i: usize) -> f64 {
        self.adjoints.get(i).copied().unwrap_or(0.0)
    }
}

#[derive(Debug, Clone)]
pub struct ArenaTape<T: Copy = f64> {
    pub ops: Vec<OpRecord<T>>,
    pub mark: usize,
    pub active: bool,
    pub tape_id: u16,
}
impl<T> Default for ArenaTape<T>
where
    T: Copy + Default,
{
    fn default() -> Self {
        Self::new()
    }
}
impl<T> ArenaTape<T>
where
    T: Copy + Default,
{
    pub fn new() -> Self {
        let tape_id = loop {
            let id = NEXT_TAPE_ID.fetch_add(1, Ordering::Relaxed);
            if id != 0 {
                break id;
            }
        };
        Self {
            ops: Vec::new(),
            mark: 0,
            active: true,
            tape_id,
        }
    }
    pub fn with_capacity(n: usize) -> Self {
        let mut x = Self::new();
        x.ops.reserve(n);
        x
    }
    pub fn len(&self) -> usize {
        self.ops.len()
    }
    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }
    pub fn mark(&self) -> usize {
        self.mark
    }
    pub fn set_mark(&mut self) {
        self.mark = self.ops.len();
    }
    pub fn rewind_to_mark(&mut self) {
        self.ops.truncate(self.mark);
    }
    pub fn variable(&mut self, v: T) -> TapeVar<T> {
        self.input(v)
    }
    pub fn input(&mut self, v: T) -> TapeVar<T> {
        let idx = self.ops.len() as u32;
        self.ops.push(OpRecord {
            lhs: CONST_IDX,
            rhs: CONST_IDX,
            w_lhs: T::default(),
            w_rhs: T::default(),
        });
        TapeVar {
            val: v,
            idx,
            tape_id: self.tape_id,
        }
    }
    pub fn constant(&self, v: T) -> TapeVar<T> {
        TapeVar {
            val: v,
            idx: CONST_IDX,
            tape_id: self.tape_id,
        }
    }
    pub fn record_unary(&mut self, a: TapeVar<T>, v: T, w: T) -> TapeVar<T> {
        assert!(
            a.idx == CONST_IDX || a.tape_id == self.tape_id,
            "TapeVar belongs to a different tape"
        );
        let idx = self.ops.len() as u32;
        self.ops.push(OpRecord {
            lhs: a.idx,
            rhs: CONST_IDX,
            w_lhs: w,
            w_rhs: T::default(),
        });
        TapeVar {
            val: v,
            idx,
            tape_id: self.tape_id,
        }
    }
    pub fn record_binary(
        &mut self,
        a: TapeVar<T>,
        b: TapeVar<T>,
        v: T,
        wa: T,
        wb: T,
    ) -> TapeVar<T> {
        assert!(
            (a.idx == CONST_IDX || a.tape_id == self.tape_id)
                && (b.idx == CONST_IDX || b.tape_id == self.tape_id),
            "TapeVar belongs to a different tape"
        );
        let idx = self.ops.len() as u32;
        self.ops.push(OpRecord {
            lhs: a.idx,
            rhs: b.idx,
            w_lhs: wa,
            w_rhs: wb,
        });
        TapeVar {
            val: v,
            idx,
            tape_id: self.tape_id,
        }
    }
    pub fn propagate_range(&self, adjoints: &mut [T], from: usize, to: usize)
    where
        T: AddAssign + Mul<Output = T>,
    {
        let end = from.min(self.ops.len());
        let start = to.min(end);
        for (k, op) in self.ops[start..end].iter().enumerate().rev() {
            let i = start + k;
            let a = adjoints[i];
            if op.lhs != CONST_IDX {
                adjoints[op.lhs as usize] += a * op.w_lhs;
            }
            if op.rhs != CONST_IDX {
                adjoints[op.rhs as usize] += a * op.w_rhs;
            }
        }
    }
    pub fn propagate_to_mark(&self, adjoints: &mut [T])
    where
        T: AddAssign + Mul<Output = T>,
    {
        let start = self.mark.min(self.ops.len());
        for (i, op) in self.ops[start..].iter().enumerate().rev() {
            let k = start + i;
            if k >= adjoints.len() {
                continue;
            }
            let a = adjoints[k];
            if op.lhs != CONST_IDX {
                adjoints[op.lhs as usize] += a * op.w_lhs;
            }
            if op.rhs != CONST_IDX {
                adjoints[op.rhs as usize] += a * op.w_rhs;
            }
        }
    }
    pub fn propagate_mark_to_start(&self, adjoints: &mut [T])
    where
        T: AddAssign + Mul<Output = T>,
    {
        self.propagate_range(adjoints, self.mark, 0)
    }
}
impl ArenaTape<f64> {
    pub fn add(&mut self, a: TapeVar<f64>, b: TapeVar<f64>) -> TapeVar<f64> {
        self.record_binary(a, b, a.val + b.val, 1.0, 1.0)
    }
    pub fn sub(&mut self, a: TapeVar<f64>, b: TapeVar<f64>) -> TapeVar<f64> {
        self.record_binary(a, b, a.val - b.val, 1.0, -1.0)
    }
    pub fn mul(&mut self, a: TapeVar<f64>, b: TapeVar<f64>) -> TapeVar<f64> {
        self.record_binary(a, b, a.val * b.val, b.val, a.val)
    }
    pub fn div(&mut self, a: TapeVar<f64>, b: TapeVar<f64>) -> TapeVar<f64> {
        self.record_binary(a, b, a.val / b.val, 1.0 / b.val, -a.val / (b.val * b.val))
    }
    pub fn exp(&mut self, a: TapeVar<f64>) -> TapeVar<f64> {
        let v = a.val.exp();
        self.record_unary(a, v, v)
    }
    pub fn ln(&mut self, a: TapeVar<f64>) -> TapeVar<f64> {
        self.record_unary(a, a.val.ln(), 1.0 / a.val)
    }
    pub fn sin(&mut self, a: TapeVar<f64>) -> TapeVar<f64> {
        self.record_unary(a, a.val.sin(), a.val.cos())
    }
}

impl ArenaTape<f64> {
    pub fn gradient_from(&self, output: TapeVar<f64>) -> AdjointResult {
        assert!(
            output.idx == CONST_IDX || output.tape_id == self.tape_id,
            "TapeVar output belongs to a different tape"
        );
        let mut a = vec![0.0; self.ops.len()];
        if output.idx != CONST_IDX && (output.idx as usize) < a.len() {
            a[output.idx as usize] = 1.0;
        }
        for (i, op) in self.ops.iter().enumerate().rev() {
            let g = a[i];
            if op.lhs != CONST_IDX {
                a[op.lhs as usize] += g * op.w_lhs;
            }
            if op.rhs != CONST_IDX {
                a[op.rhs as usize] += g * op.w_rhs;
            }
        }
        AdjointResult {
            adjoints: a,
            output: output.idx,
        }
    }
}

thread_local! {static ACTIVE:RefCell<Option<ArenaTape<f64>>>=const{RefCell::new(None)};}
pub struct TapeScope;
impl TapeScope {
    pub fn start_recording() {
        ACTIVE.with(|a| *a.borrow_mut() = Some(ArenaTape::new()));
    }
    pub fn start_recording_with_capacity(capacity: usize) {
        ACTIVE.with(|a| *a.borrow_mut() = Some(ArenaTape::with_capacity(capacity)));
    }
    pub fn stop_recording() -> Option<ArenaTape<f64>> {
        ACTIVE.with(|a| a.borrow_mut().take())
    }
    pub fn with_active<R, F: FnOnce(&mut ArenaTape<f64>) -> R>(f: F) -> R {
        let had = ACTIVE.with(|a| a.borrow().is_some());
        if !had {
            Self::start_recording();
        }
        ACTIVE.with(|a| f(a.borrow_mut().as_mut().expect("active tape")))
    }
    pub fn set_mark() {
        ACTIVE.with(|a| {
            if let Some(t) = a.borrow_mut().as_mut() {
                t.set_mark()
            }
        })
    }
    pub fn rewind_to_mark() {
        ACTIVE.with(|a| {
            if let Some(t) = a.borrow_mut().as_mut() {
                t.rewind_to_mark()
            }
        })
    }
    pub fn gradient_from(output: TapeVar<f64>) -> AdjointResult {
        ACTIVE.with(|a| {
            a.borrow()
                .as_ref()
                .expect("active tape")
                .gradient_from(output)
        })
    }
}
pub fn start_recording() {
    TapeScope::start_recording()
}
pub fn stop_recording() -> Option<ArenaTape<f64>> {
    TapeScope::stop_recording()
}
pub fn set_mark() {
    TapeScope::set_mark()
}
pub fn rewind_to_mark() {
    TapeScope::rewind_to_mark()
}
pub fn gradient_from(output: TapeVar<f64>) -> AdjointResult {
    TapeScope::gradient_from(output)
}

fn active_binary(a: TapeVar<f64>, b: TapeVar<f64>, v: f64, wa: f64, wb: f64) -> TapeVar<f64> {
    if a.idx == CONST_IDX && b.idx == CONST_IDX {
        return TapeVar {
            val: v,
            idx: CONST_IDX,
            tape_id: 0,
        };
    }
    ACTIVE.with(|cell| {
        if let Some(t) = cell.borrow_mut().as_mut() {
            t.record_binary(a, b, v, wa, wb)
        } else if a.idx != CONST_IDX || b.idx != CONST_IDX {
            panic!("attempted to operate on TapeVar without an active tape")
        } else {
            TapeVar {
                val: v,
                idx: CONST_IDX,
                tape_id: 0,
            }
        }
    })
}
fn active_unary(a: TapeVar<f64>, v: f64, w: f64) -> TapeVar<f64> {
    if a.idx == CONST_IDX {
        return TapeVar {
            val: v,
            idx: CONST_IDX,
            tape_id: 0,
        };
    }
    ACTIVE.with(|cell| {
        if let Some(t) = cell.borrow_mut().as_mut() {
            t.record_unary(a, v, w)
        } else {
            panic!("attempted to operate on TapeVar without an active tape")
        }
    })
}
impl fmt::Display for TapeVar<f64> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.val)
    }
}
impl PartialEq for TapeVar<f64> {
    fn eq(&self, o: &Self) -> bool {
        self.val == o.val
    }
}
impl PartialOrd for TapeVar<f64> {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        self.val.partial_cmp(&o.val)
    }
}
impl Add for TapeVar<f64> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        active_binary(self, o, self.val + o.val, 1.0, 1.0)
    }
}
impl Sub for TapeVar<f64> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        active_binary(self, o, self.val - o.val, 1.0, -1.0)
    }
}
impl Mul for TapeVar<f64> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        active_binary(self, o, self.val * o.val, o.val, self.val)
    }
}
impl Div for TapeVar<f64> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        active_binary(
            self,
            o,
            self.val / o.val,
            1.0 / o.val,
            -self.val / (o.val * o.val),
        )
    }
}
impl Neg for TapeVar<f64> {
    type Output = Self;
    fn neg(self) -> Self {
        active_unary(self, -self.val, -1.0)
    }
}
impl Add<f64> for TapeVar<f64> {
    type Output = Self;
    fn add(self, o: f64) -> Self {
        active_binary(
            self,
            TapeVar {
                val: o,
                idx: CONST_IDX,
                tape_id: self.tape_id,
            },
            self.val + o,
            1.0,
            0.0,
        )
    }
}
impl Sub<f64> for TapeVar<f64> {
    type Output = Self;
    fn sub(self, o: f64) -> Self {
        active_binary(
            self,
            TapeVar {
                val: o,
                idx: CONST_IDX,
                tape_id: self.tape_id,
            },
            self.val - o,
            1.0,
            0.0,
        )
    }
}
impl Mul<f64> for TapeVar<f64> {
    type Output = Self;
    fn mul(self, o: f64) -> Self {
        active_unary(self, self.val * o, o)
    }
}
impl Div<f64> for TapeVar<f64> {
    type Output = Self;
    fn div(self, o: f64) -> Self {
        active_unary(self, self.val / o, 1.0 / o)
    }
}
impl AddAssign for TapeVar<f64> {
    fn add_assign(&mut self, o: Self) {
        *self = *self + o
    }
}
impl SubAssign for TapeVar<f64> {
    fn sub_assign(&mut self, o: Self) {
        *self = *self - o
    }
}
impl MulAssign for TapeVar<f64> {
    fn mul_assign(&mut self, o: Self) {
        *self = *self * o
    }
}
impl DivAssign for TapeVar<f64> {
    fn div_assign(&mut self, o: Self) {
        *self = *self / o
    }
}
impl Add<TapeVar<f64>> for f64 {
    type Output = TapeVar<f64>;
    fn add(self, o: TapeVar<f64>) -> TapeVar<f64> {
        o + self
    }
}
impl Sub<TapeVar<f64>> for f64 {
    type Output = TapeVar<f64>;
    fn sub(self, o: TapeVar<f64>) -> TapeVar<f64> {
        active_binary(
            TapeVar {
                val: self,
                idx: CONST_IDX,
                tape_id: o.tape_id,
            },
            o,
            self - o.val,
            0.0,
            -1.0,
        )
    }
}
impl Mul<TapeVar<f64>> for f64 {
    type Output = TapeVar<f64>;
    fn mul(self, o: TapeVar<f64>) -> TapeVar<f64> {
        o * self
    }
}
impl Div<TapeVar<f64>> for f64 {
    type Output = TapeVar<f64>;
    fn div(self, o: TapeVar<f64>) -> TapeVar<f64> {
        active_binary(
            TapeVar {
                val: self,
                idx: CONST_IDX,
                tape_id: o.tape_id,
            },
            o,
            self / o.val,
            0.0,
            -self / (o.val * o.val),
        )
    }
}
impl Real for TapeVar<f64> {
    fn from_f64(v: f64) -> Self {
        TapeVar {
            val: v,
            idx: CONST_IDX,
            tape_id: 0,
        }
    }
    fn val(self) -> f64 {
        self.val
    }
    fn abs(self) -> Self {
        active_unary(
            self,
            self.val.abs(),
            if self.val >= 0.0 { 1.0 } else { -1.0 },
        )
    }
    fn sqrt(self) -> Self {
        let g = self.val.sqrt();
        active_unary(self, g, if g == 0.0 { 0.0 } else { 0.5 / g })
    }
    fn exp(self) -> Self {
        let g = self.val.exp();
        active_unary(self, g, g)
    }
    fn ln(self) -> Self {
        active_unary(self, self.val.ln(), 1.0 / self.val)
    }
    fn sin(self) -> Self {
        active_unary(self, self.val.sin(), self.val.cos())
    }
    fn cos(self) -> Self {
        active_unary(self, self.val.cos(), -self.val.sin())
    }
    fn tan(self) -> Self {
        let c = self.val.cos();
        active_unary(self, self.val.tan(), 1.0 / (c * c))
    }
    fn atan(self) -> Self {
        active_unary(self, self.val.atan(), 1.0 / (1.0 + self.val * self.val))
    }
    fn sinh(self) -> Self {
        active_unary(self, self.val.sinh(), self.val.cosh())
    }
    fn cosh(self) -> Self {
        active_unary(self, self.val.cosh(), self.val.sinh())
    }
    fn tanh(self) -> Self {
        let g = self.val.tanh();
        active_unary(self, g, 1.0 - g * g)
    }
    fn powf(self, p: f64) -> Self {
        let g = self.val.powf(p);
        active_unary(self, g, p * self.val.powf(p - 1.0))
    }
    fn pow(self, p: Self) -> Self {
        (p * self.ln()).exp()
    }
    fn mul_add(self, a: Self, b: Self) -> Self {
        self * a + b
    }
    fn erf(self) -> Self {
        let g = (2.0 / std::f64::consts::PI.sqrt()) * (-self.val * self.val).exp();
        active_unary(
            self,
            crate::dual::Dual::<1>::from_value(self.val).erf().val(),
            g,
        )
    }
    fn erfc(self) -> Self {
        Self::one() - self.erf()
    }
    fn norm_cdf(self) -> Self {
        active_unary(
            self,
            0.5 * (1.0
                + crate::dual::Dual::<1>::from_value(self.val * std::f64::consts::FRAC_1_SQRT_2)
                    .erf()
                    .val()),
            (-0.5 * self.val * self.val).exp() / (2.0 * std::f64::consts::PI).sqrt(),
        )
    }
    fn norm_pdf(self) -> Self {
        let p = (-0.5 * self.val * self.val).exp() / (2.0 * std::f64::consts::PI).sqrt();
        active_unary(self, p, -self.val * p)
    }
    fn min_val(self, o: Self) -> Self {
        if self.val <= o.val { self } else { o }
    }
    fn max_val(self, o: Self) -> Self {
        if self.val >= o.val { self } else { o }
    }
}
