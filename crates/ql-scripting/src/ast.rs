#[derive(Clone, Debug)]
pub enum Expr {
    Const(f64),
    Spot,
    Var(usize),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
    Max(Box<Expr>, Box<Expr>),
    Min(Box<Expr>, Box<Expr>),
    SmoothStep {
        x: Box<Expr>,
        threshold: Box<Expr>,
        half_width: f64,
    },
}
#[derive(Clone, Debug)]
pub enum Stmt {
    Assign {
        var: usize,
        expr: Expr,
    },
    Pay {
        amount: Expr,
        discount_time: f64,
    },
    Autocall {
        spot_expr: Expr,
        barrier: f64,
        half_width: f64,
        redemption: Expr,
        discount_time: f64,
    },
}
#[derive(Clone, Debug)]
pub struct EventStep {
    pub time: f64,
    pub statements: Vec<Stmt>,
}
#[derive(Clone, Debug)]
pub struct ScriptedProduct {
    pub event_steps: Vec<EventStep>,
    pub num_vars: usize,
}
impl ScriptedProduct {
    pub fn asian_call(times: Vec<f64>, strike: f64, notional: f64) -> Self {
        let n = times.len();
        let mut st = Vec::new();
        for (i, &t) in times.iter().enumerate() {
            let e = if i == 0 {
                Expr::Spot
            } else {
                Expr::Add(Box::new(Expr::Var(0)), Box::new(Expr::Spot))
            };
            let mut ss = vec![Stmt::Assign { var: 0, expr: e }];
            if i + 1 == n {
                ss.push(Stmt::Pay {
                    amount: Expr::Mul(
                        Box::new(Expr::Const(notional)),
                        Box::new(Expr::Max(
                            Box::new(Expr::Sub(
                                Box::new(Expr::Div(
                                    Box::new(Expr::Var(0)),
                                    Box::new(Expr::Const(n as f64)),
                                )),
                                Box::new(Expr::Const(strike)),
                            )),
                            Box::new(Expr::Const(0.)),
                        )),
                    ),
                    discount_time: t,
                });
            }
            st.push(EventStep {
                time: t,
                statements: ss,
            });
        }
        Self {
            event_steps: st,
            num_vars: 1,
        }
    }
    pub fn autocallable(
        times: Vec<f64>,
        ab: f64,
        cb: f64,
        cr: f64,
        n: f64,
        put: f64,
        w: f64,
    ) -> Self {
        let mut st = Vec::new();
        for &t in &times {
            st.push(EventStep {
                time: t,
                statements: vec![
                    Stmt::Pay {
                        amount: Expr::Mul(
                            Box::new(Expr::Const(n * cr)),
                            Box::new(Expr::SmoothStep {
                                x: Box::new(Expr::Spot),
                                threshold: Box::new(Expr::Const(cb)),
                                half_width: w,
                            }),
                        ),
                        discount_time: t,
                    },
                    Stmt::Autocall {
                        spot_expr: Expr::Spot,
                        barrier: ab,
                        half_width: w,
                        redemption: Expr::Const(n),
                        discount_time: t,
                    },
                ],
            });
        }
        if let Some(last) = st.last_mut() {
            let t = last.time;
            let put_loss = Expr::Sub(Box::new(Expr::Const(put)), Box::new(Expr::Spot));
            // Outside the fuzzy band this is exactly max(put - spot, 0).
            // Inside the band it keeps the terminal payoff differentiable,
            // which is essential for pathwise fuzzy AAD.
            let smooth_put_loss = Expr::Mul(
                Box::new(put_loss.clone()),
                Box::new(Expr::SmoothStep {
                    x: Box::new(put_loss),
                    threshold: Box::new(Expr::Const(0.0)),
                    half_width: w,
                }),
            );
            last.statements.push(Stmt::Pay {
                amount: Expr::Sub(
                    Box::new(Expr::Const(n)),
                    Box::new(Expr::Mul(
                        Box::new(Expr::Const(n / put)),
                        Box::new(smooth_put_loss),
                    )),
                ),
                discount_time: t,
            });
        }
        Self {
            event_steps: st,
            num_vars: 0,
        }
    }
    pub fn tarn(times: Vec<f64>, strike: f64, target: f64, n: f64, w: f64) -> Self {
        let mut st = Vec::new();
        for &t in &times {
            let c = Expr::Max(
                Box::new(Expr::Sub(
                    Box::new(Expr::Spot),
                    Box::new(Expr::Const(strike)),
                )),
                Box::new(Expr::Const(0.)),
            );
            let rem = Expr::Max(
                Box::new(Expr::Sub(
                    Box::new(Expr::Const(target)),
                    Box::new(Expr::Var(0)),
                )),
                Box::new(Expr::Const(0.)),
            );
            let a = Expr::Min(Box::new(c), Box::new(rem));
            st.push(EventStep {
                time: t,
                statements: vec![
                    Stmt::Assign { var: 1, expr: a },
                    Stmt::Pay {
                        amount: Expr::Mul(Box::new(Expr::Const(n)), Box::new(Expr::Var(1))),
                        discount_time: t,
                    },
                    Stmt::Assign {
                        var: 0,
                        expr: Expr::Add(Box::new(Expr::Var(0)), Box::new(Expr::Var(1))),
                    },
                    Stmt::Autocall {
                        // A capped coupon reaches the target exactly.  Shift the
                        // fuzzy trigger by one half-width so that reaching the
                        // cap is an unambiguous full autocall (rather than the
                        // 50% weight assigned to the centre of a fuzzy step).
                        spot_expr: Expr::Add(Box::new(Expr::Var(0)), Box::new(Expr::Const(w))),
                        barrier: target,
                        half_width: w,
                        redemption: Expr::Const(n),
                        discount_time: t,
                    },
                ],
            });
        }
        Self {
            event_steps: st,
            num_vars: 2,
        }
    }
}
