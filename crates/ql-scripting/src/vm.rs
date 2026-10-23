use crate::ast::*;
use ql_ad::Real;
#[derive(Copy, Clone, Debug)]
pub enum OpCode {
    BeginStep {
        step_idx: usize,
    },
    LoadConst(f64),
    LoadSpot,
    LoadVar(usize),
    StoreVar(usize),
    Add,
    Sub,
    Mul,
    Div,
    Max,
    Min,
    SmoothStep(f64),
    Pay {
        discount_time: f64,
    },
    Autocall {
        barrier: f64,
        half_width: f64,
        discount_time: f64,
    },
}
#[derive(Clone, Debug)]
pub struct BytecodeProgram {
    pub ops: Vec<OpCode>,
    pub times: Vec<f64>,
    pub num_vars: usize,
}
impl BytecodeProgram {
    pub fn compile(p: &ScriptedProduct) -> Self {
        let mut ops = Vec::new();
        for (i, e) in p.event_steps.iter().enumerate() {
            ops.push(OpCode::BeginStep { step_idx: i });
            for st in &e.statements {
                compile_stmt(st, &mut ops);
            }
        }
        Self {
            ops,
            times: p.event_steps.iter().map(|x| x.time).collect(),
            num_vars: p.num_vars,
        }
    }
    pub fn eval_path<T: Real>(
        &self,
        spots: &[T],
        rate: T,
        vars: &mut [T],
        stack: &mut Vec<T>,
    ) -> T {
        vars.fill(T::zero());
        stack.clear();
        let mut pv = T::zero();
        let mut alive = T::one();
        let mut cur = T::zero();
        for op in &self.ops {
            match *op {
                OpCode::BeginStep { step_idx } => {
                    cur = *spots
                        .get(step_idx)
                        .expect("path is shorter than program times")
                }
                OpCode::LoadConst(x) => stack.push(T::from_f64(x)),
                OpCode::LoadSpot => stack.push(cur),
                OpCode::LoadVar(i) => stack.push(vars[i]),
                OpCode::StoreVar(i) => vars[i] = stack.pop().expect("invalid bytecode stack"),
                OpCode::Add => bin(stack, |a, b| a + b),
                OpCode::Sub => bin(stack, |a, b| a - b),
                OpCode::Mul => bin(stack, |a, b| a * b),
                OpCode::Div => bin(stack, |a, b| a / b),
                OpCode::Max => bin(stack, |a, b| a.max_val(b)),
                OpCode::Min => bin(stack, |a, b| a.min_val(b)),
                OpCode::SmoothStep(w) => {
                    let th = stack.pop().expect("invalid bytecode stack");
                    let x = stack.pop().expect("invalid bytecode stack");
                    stack.push(x.smooth_heaviside(th, w));
                }
                OpCode::Pay { discount_time } => {
                    let a = stack.pop().expect("invalid bytecode stack");
                    pv += alive * a * (-rate * T::from_f64(discount_time)).exp();
                }
                OpCode::Autocall {
                    barrier,
                    half_width,
                    discount_time,
                } => {
                    let red = stack.pop().expect("invalid bytecode stack");
                    let x = stack.pop().expect("invalid bytecode stack");
                    let w = x.smooth_heaviside(T::from_f64(barrier), half_width);
                    pv += alive * w * red * (-rate * T::from_f64(discount_time)).exp();
                    alive *= T::one() - w;
                }
            }
        }
        pv
    }
}
fn bin<T: Real, F: Fn(T, T) -> T>(s: &mut Vec<T>, f: F) {
    let b = s.pop().expect("invalid bytecode stack");
    let a = s.pop().expect("invalid bytecode stack");
    s.push(f(a, b));
}
fn compile_expr(e: &Expr, o: &mut Vec<OpCode>) {
    match e {
        Expr::Const(x) => o.push(OpCode::LoadConst(*x)),
        Expr::Spot => o.push(OpCode::LoadSpot),
        Expr::Var(i) => o.push(OpCode::LoadVar(*i)),
        Expr::Add(a, b) => {
            compile_expr(a, o);
            compile_expr(b, o);
            o.push(OpCode::Add)
        }
        Expr::Sub(a, b) => {
            compile_expr(a, o);
            compile_expr(b, o);
            o.push(OpCode::Sub)
        }
        Expr::Mul(a, b) => {
            compile_expr(a, o);
            compile_expr(b, o);
            o.push(OpCode::Mul)
        }
        Expr::Div(a, b) => {
            compile_expr(a, o);
            compile_expr(b, o);
            o.push(OpCode::Div)
        }
        Expr::Max(a, b) => {
            compile_expr(a, o);
            compile_expr(b, o);
            o.push(OpCode::Max)
        }
        Expr::Min(a, b) => {
            compile_expr(a, o);
            compile_expr(b, o);
            o.push(OpCode::Min)
        }
        Expr::SmoothStep {
            x,
            threshold,
            half_width,
        } => {
            compile_expr(x, o);
            compile_expr(threshold, o);
            o.push(OpCode::SmoothStep(*half_width))
        }
    }
}
fn compile_stmt(s: &Stmt, o: &mut Vec<OpCode>) {
    match s {
        Stmt::Assign { var, expr } => {
            compile_expr(expr, o);
            o.push(OpCode::StoreVar(*var))
        }
        Stmt::Pay {
            amount,
            discount_time,
        } => {
            compile_expr(amount, o);
            o.push(OpCode::Pay {
                discount_time: *discount_time,
            })
        }
        Stmt::Autocall {
            spot_expr,
            barrier,
            half_width,
            redemption,
            discount_time,
        } => {
            compile_expr(spot_expr, o);
            compile_expr(redemption, o);
            o.push(OpCode::Autocall {
                barrier: *barrier,
                half_width: *half_width,
                discount_time: *discount_time,
            })
        }
    }
}
#[derive(Copy, Clone, Debug)]
pub struct ScriptMcEngine {
    pub n_paths: usize,
    pub seed: u64,
}
impl ScriptMcEngine {
    pub fn price_gbm<T: Real>(&self, p: &BytecodeProgram, spot: T, vol: T, rate: T, div: T) -> T {
        let rng = ql_mc::PhiloxRng::new(self.seed);
        let n = self.n_paths.max(1);
        let mut sum = T::zero();
        let mut vars = vec![T::zero(); p.num_vars];
        let mut stack = Vec::with_capacity(16);
        let mut pos = vec![T::zero(); p.times.len()];
        let mut neg = vec![T::zero(); p.times.len()];
        for path in 0..n {
            let mut sp = spot;
            let mut sn = spot;
            let mut prev = 0.;
            for (j, &tm) in p.times.iter().enumerate() {
                let dt = tm - prev;
                // Consume every component of each Philox block exactly once;
                // this keeps `n_paths` equal to the number of antithetic pairs
                // while avoiding the fourfold waste of repeatedly selecting
                // component zero from the same block.
                let block = (path / 4) as u64;
                let lane = path % 4;
                let z = T::from_f64(rng.normal4(block, j as u32)[lane]);
                let drift = (rate - div - T::from_f64(0.5) * vol * vol) * dt;
                let vol_sqrt_dt = vol * dt.sqrt();
                let diff = vol_sqrt_dt * z;
                sp *= (drift + diff).exp();
                sn *= (drift - diff).exp();
                pos[j] = sp;
                neg[j] = sn;
                prev = tm;
            }
            let a = p.eval_path(&pos, rate, &mut vars, &mut stack);
            let b = p.eval_path(&neg, rate, &mut vars, &mut stack);
            sum += (a + b) * T::from_f64(0.5);
        }
        sum / T::from_f64(n as f64)
    }
}
