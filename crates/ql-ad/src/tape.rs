use crate::scalar::Real;
use std::cell::RefCell;
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::atomic::{AtomicU32, Ordering};

pub const CONST_IDX: u32 = u32::MAX;
static NEXT_TAPE_ID: AtomicU32 = AtomicU32::new(1);
#[derive(Copy, Clone, Debug, Default)]
pub struct TapeVar<T: Copy = f64> {
    pub val: T,
    pub idx: u32,
    pub tape_id: u32,
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
    pub tape_id: u32,
}
impl AdjointResult {
    pub fn gradient(&self, var: TapeVar<f64>) -> f64 {
        assert!(
            var.idx == CONST_IDX || var.tape_id == self.tape_id,
            "cross-tape variable"
        );
        if var.idx == CONST_IDX {
            0.0
        } else {
            self.adjoints.get(var.idx as usize).copied().unwrap_or(0.0)
        }
    }
    pub fn try_gradient(&self, var: TapeVar<f64>) -> Result<f64, &'static str> {
        if var.idx != CONST_IDX && var.tape_id != self.tape_id {
            return Err("cross-tape variable");
        }
        Ok(if var.idx == CONST_IDX {
            0.0
        } else {
            self.adjoints.get(var.idx as usize).copied().unwrap_or(0.0)
        })
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
    pub tape_id: u32,
}
impl<T> Default for ArenaTape<T>
where
    T: Copy + Default + PartialEq,
{
    fn default() -> Self {
        Self::new()
    }
}
impl<T> ArenaTape<T>
where
    T: Copy + Default + PartialEq,
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
            if a == T::default() {
                continue;
            }
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
            if a == T::default() {
                continue;
            }
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
            if g == 0.0 {
                continue;
            }
            if op.lhs != CONST_IDX && op.w_lhs != 0.0 {
                a[op.lhs as usize] += g * op.w_lhs;
            }
            if op.rhs != CONST_IDX && op.w_rhs != 0.0 {
                a[op.rhs as usize] += g * op.w_rhs;
            }
        }
        AdjointResult {
            adjoints: a,
            output: output.idx,
            tape_id: self.tape_id,
        }
    }
}

thread_local! {static ACTIVE:RefCell<Option<ArenaTape<f64>>>=const{RefCell::new(None)};}
pub struct TapeScope;
/// A scoped active tape.  Creating a guard replaces the current active tape;
/// dropping it restores the previous tape (or leaves recording inactive).
/// Restoration is performed during unwinding as well, so early returns and
/// panics cannot leak a tape into unrelated pricing work.
pub struct TapeGuard {
    previous: Option<ArenaTape<f64>>,
    active: bool,
}
impl TapeGuard {
    /// Finish this scope and return the tape recorded by the guard.  The
    /// previous tape, if any, is restored before returning.
    pub fn stop(mut self) -> Option<ArenaTape<f64>> {
        self.active = false;
        ACTIVE.with(|a| {
            let mut slot = a.borrow_mut();
            let current = slot.take();
            *slot = self.previous.take();
            current
        })
    }
}
impl Drop for TapeGuard {
    fn drop(&mut self) {
        if !self.active {
            return;
        }
        ACTIVE.with(|a| {
            let mut slot = a.borrow_mut();
            let _current = slot.take();
            *slot = self.previous.take();
        });
    }
}
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
        ACTIVE.with(|a| {
            f(a.borrow_mut().as_mut().expect(
                "no active tape; call TapeScope::start_recording or create a TapeScope::guard",
            ))
        })
    }
    pub fn try_with_active<R, F: FnOnce(&mut ArenaTape<f64>) -> R>(f: F) -> Option<R> {
        ACTIVE.with(|a| a.borrow_mut().as_mut().map(f))
    }
    pub fn guard() -> TapeGuard {
        Self::guard_with_capacity(0)
    }
    pub fn guard_with_capacity(capacity: usize) -> TapeGuard {
        let previous = ACTIVE.with(|a| a.borrow_mut().take());
        let tape = if capacity == 0 {
            ArenaTape::new()
        } else {
            ArenaTape::with_capacity(capacity)
        };
        ACTIVE.with(|a| *a.borrow_mut() = Some(tape));
        TapeGuard {
            previous,
            active: true,
        }
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
        active_unary(self, self.val + o, 1.0)
    }
}
impl Sub<f64> for TapeVar<f64> {
    type Output = Self;
    fn sub(self, o: f64) -> Self {
        active_unary(self, self.val - o, 1.0)
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
        -o + self
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
        active_unary(o, self / o.val, -self / (o.val * o.val))
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
        if p == 0.0 {
            return TapeVar::from_f64(1.0);
        }
        if p == 1.0 {
            return self;
        }
        if p == 2.0 {
            return self * self;
        }
        if self.val == 0.0 {
            let g = self.val.powf(p);
            let w = if p > 1.0 {
                0.0
            } else {
                p * self.val.powf(p - 1.0)
            };
            return active_unary(self, g, w);
        }
        let g = self.val.powf(p);
        active_unary(self, g, p * self.val.powf(p - 1.0))
    }
    fn pow(self, p: Self) -> Self {
        if p.idx == CONST_IDX {
            return self.powf(p.val);
        }
        if self.val == 0.0 {
            return self.powf(p.val);
        }
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
        active_unary(
            self,
            libm::erfc(self.val),
            -2.0 / std::f64::consts::PI.sqrt() * (-self.val * self.val).exp(),
        )
    }
    fn norm_cdf(self) -> Self {
        active_unary(
            self,
            0.5 * libm::erfc(-self.val * std::f64::consts::FRAC_1_SQRT_2),
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

#[derive(Clone, Debug)]
pub struct SavineRevolveCheckpointSpec<T: Real = f64> {
    pub num_time_steps: usize,
    pub num_checkpoints: usize,
    pub initial_spot: T,
    pub drift_rate: T,
    pub volatility: T,
    pub mean_reversion_kappa: T,
    pub dt: T,
    pub strike: T,
}
#[derive(Clone, Debug)]
pub struct SavineRevolveCheckpointReport<T: Real = f64> {
    pub terminal_option_payoff_value: T,
    pub checkpointed_delta: T,
    pub checkpointed_vega: T,
    pub checkpointed_rho: T,
    pub full_tape_reference_delta: T,
    pub full_tape_reference_vega: T,
    pub max_adjoint_discrepancy: T,
    pub revolve_repetition_tier_t: usize,
    pub first_revolve_checkpoint_step: usize,
    pub memory_compression_ratio: T,
}
pub fn execute_savine_revolve_checkpointed_adjoint<T: Real>(
    s: &SavineRevolveCheckpointSpec<T>,
) -> SavineRevolveCheckpointReport<T> {
    let n = s.num_time_steps.max(2);
    let c = s.num_checkpoints.clamp(1, n);
    // Griewank--Walther binomial Revolve schedule.  The repetition tier t is
    // the smallest integer with beta(c,t) >= n, beta(c,t)=C(c+t,t).
    fn beta(c: usize, t: usize) -> usize {
        let mut b = 1usize;
        for k in 1..=t {
            b = b.saturating_mul(c.saturating_add(k)) / k;
        }
        b
    }
    let mut tier = 0usize;
    while beta(c, tier) < n {
        tier += 1;
    }
    // First split is beta(c-1,t)=C(c-1+t,t), clamped to the interior.
    let first = beta(c.saturating_sub(1), tier).clamp(1, n - 1);
    let dt = s.dt;
    let step = |sn: T, i: usize| -> (T, T, T, T, T) {
        let z = (0.73 * (i + 1) as f64).sin() + 0.5 * (1.19 * ((i + 1) * (i + 1)) as f64).cos();
        let e = ((s.drift_rate
            - s.mean_reversion_kappa * (sn - s.initial_spot) / s.initial_spot
            - T::from_f64(0.5) * s.volatility * s.volatility)
            * dt
            + s.volatility * dt.sqrt() * z)
            .exp();
        let out = sn * e;
        let ds = e * (T::one() - s.mean_reversion_kappa * sn / s.initial_spot * dt);
        let dr = out * dt;
        let dv = out * (-s.volatility * dt + dt.sqrt() * z);
        let d0 = out * s.mean_reversion_kappa * sn / (s.initial_spot * s.initial_spot) * dt;
        (out, ds, dr, dv, d0)
    };
    let mut states = Vec::with_capacity(n + 1);
    let mut x = s.initial_spot;
    states.push(x);
    for i in 0..n {
        x = step(x, i).0;
        states.push(x);
    }
    // Griewank--Walther binomial Revolve forward checkpoint recurrence.  The
    // cursor advances through the complete [0,n) interval; in particular the
    // final segment is no longer silently omitted when `n` is not itself a
    // checkpoint index.
    let mut checkpoint_indices = vec![0usize];
    let mut cur = 0usize;
    for rem_c in (1..c).rev() {
        if cur + 1 >= n {
            break;
        }
        let rem_steps = n - cur;
        let mut t_rem = 0usize;
        while beta(rem_c + 1, t_rem) < rem_steps {
            t_rem += 1;
        }
        let adv = if cur == 0 {
            first
        } else {
            beta(rem_c, t_rem).clamp(1, rem_steps - 1)
        };
        cur += adv;
        checkpoint_indices.push(cur);
    }
    let checkpoints: Vec<T> = checkpoint_indices.iter().map(|&i| states[i]).collect();
    let y = x - s.strike;
    let root = (y * y + T::from_f64(1e-8)).sqrt();
    let df = (-s.drift_rate * T::from_f64(n as f64) * dt).exp();
    let payoff = T::from_f64(0.5) * (y + root) * df;
    let terminal_adj = T::from_f64(0.5) * (T::one() + y / root) * df;
    let mut a_s = terminal_adj;
    let mut a0 = T::zero();
    let mut av = T::zero();
    let mut ar = -T::from_f64(n as f64) * dt * payoff;
    for i in (0..n).rev() {
        let sn = states[i];
        let (_, ds, dr, dv, d0) = step(sn, i);
        ar += a_s * dr;
        av += a_s * dv;
        a0 += a_s * d0;
        a_s *= ds;
    }
    let full_delta = a_s + a0;
    let mut cs = terminal_adj;
    let mut c0 = T::zero();
    let mut cv = T::zero();
    let mut cr = -T::from_f64(n as f64) * dt * payoff;
    // Re-simulate each exact Revolve segment and sweep its adjoint in reverse.
    for seg in (0..checkpoint_indices.len()).rev() {
        let begin = checkpoint_indices[seg];
        let end = if seg + 1 < checkpoint_indices.len() {
            checkpoint_indices[seg + 1]
        } else {
            n
        };
        let mut local = Vec::with_capacity(end - begin + 1);
        let mut z = checkpoints[seg];
        local.push(z);
        for i in begin..end {
            z = step(z, i).0;
            local.push(z);
        }
        for off in (0..end - begin).rev() {
            let i = begin + off;
            let (_, ds, dr, dv, d0) = step(local[off], i);
            cr += cs * dr;
            cv += cs * dv;
            c0 += cs * d0;
            cs *= ds;
        }
    }
    let cp_delta = cs + c0;
    let discrepancy = (cp_delta - full_delta)
        .abs()
        .max_val((cv - av).abs())
        .max_val((cr - ar).abs());
    SavineRevolveCheckpointReport {
        terminal_option_payoff_value: payoff,
        checkpointed_delta: cp_delta,
        checkpointed_vega: cv,
        checkpointed_rho: cr,
        full_tape_reference_delta: full_delta,
        full_tape_reference_vega: av,
        max_adjoint_discrepancy: discrepancy,
        revolve_repetition_tier_t: tier,
        first_revolve_checkpoint_step: first,
        memory_compression_ratio: T::from_f64(n as f64 / c as f64),
    }
}

#[cfg(test)]
mod domain_suite_34_tests {
    use super::*;
    #[test]
    fn revolve_binomial_schedule_and_compression() {
        let s = SavineRevolveCheckpointSpec {
            num_time_steps: 32,
            num_checkpoints: 4,
            initial_spot: 100.0,
            drift_rate: 0.02,
            volatility: 0.2,
            mean_reversion_kappa: 0.1,
            dt: 1.0 / 32.0,
            strike: 100.0,
        };
        let r = execute_savine_revolve_checkpointed_adjoint(&s);
        assert!(r.first_revolve_checkpoint_step >= 1 && r.first_revolve_checkpoint_step < 32);
        assert!((r.memory_compression_ratio.val() - 8.0).abs() < 1e-12);
        assert!(r.max_adjoint_discrepancy.val() < 1e-8);
    }
}
