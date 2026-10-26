use ql_ad::Real;
use ql_models::HestonParams;
use ql_vol::black_scholes_price;

fn solve_tridiag<T: Real>(
    lower: &[T],
    diag: &[T],
    upper: &[T],
    rhs: &[T],
    out: &mut [T],
    scratch: &mut [T],
) {
    let n = diag.len();
    assert!(n > 0 && out.len() >= n && scratch.len() >= n);
    scratch[0] = if n > 1 { upper[0] / diag[0] } else { T::zero() };
    out[0] = rhs[0] / diag[0];
    for i in 1..n {
        // The directional sweeps store `lower[i]` on the row-i edge (the
        // slot lower[0] is unused), so consume that same row-aligned layout
        // here.  This differs from ql_math's compact n-1 representation.
        let den = diag[i] - lower[i] * scratch[i - 1];
        scratch[i] = if i + 1 < n { upper[i] / den } else { T::zero() };
        out[i] = (rhs[i] - lower[i] * out[i - 1]) / den;
    }
    for i in (0..n.saturating_sub(1)).rev() {
        out[i] -= scratch[i] * out[i + 1];
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum AdiScheme {
    Douglas,
    CraigSneyd,
    ModifiedCraigSneyd,
    HundsdorferVerwer,
}

#[derive(Copy, Clone, Debug)]
pub struct Fd2dHestonEngine {
    pub nx: usize,
    pub nv: usize,
    pub nt: usize,
    pub theta: f64,
    pub scheme: AdiScheme,
}

impl Default for Fd2dHestonEngine {
    fn default() -> Self {
        Self {
            nx: 80,
            nv: 40,
            nt: 60,
            theta: 0.5,
            scheme: AdiScheme::ModifiedCraigSneyd,
        }
    }
}

fn barrier_eps(b: f64) -> f64 {
    1e-12 * b.abs().max(1.0)
}

impl Fd2dHestonEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price_european<T: Real>(
        &self,
        p: &HestonParams<T>,
        spot: T,
        strike: T,
        rate: T,
        div: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        self.solve(
            p,
            spot,
            strike,
            None,
            None,
            T::zero(),
            rate,
            div,
            expiry,
            is_call,
            false,
            None,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn price_knockout_rebate<T: Real>(
        &self,
        p: &HestonParams<T>,
        spot: T,
        lower_barrier: Option<f64>,
        upper_barrier: Option<f64>,
        rebate: T,
        rate: T,
        div: T,
        expiry: T,
    ) -> T {
        self.solve(
            p,
            spot,
            T::zero(),
            lower_barrier,
            upper_barrier,
            T::zero(),
            rate,
            div,
            expiry,
            true,
            false,
            Some(rebate),
            lower_barrier.or(upper_barrier).map(f64::ln),
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn price_american<T: Real>(
        &self,
        p: &HestonParams<T>,
        spot: T,
        strike: T,
        rate: T,
        div: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        self.solve(
            p,
            spot,
            strike,
            None,
            None,
            T::zero(),
            rate,
            div,
            expiry,
            is_call,
            true,
            None,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn price_knockout_barrier<T: Real>(
        &self,
        p: &HestonParams<T>,
        spot: T,
        strike: T,
        lower_barrier: Option<f64>,
        upper_barrier: Option<f64>,
        rebate: T,
        rate: T,
        div: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        if lower_barrier.is_some_and(|b| spot.val() <= b + barrier_eps(b))
            || upper_barrier.is_some_and(|b| spot.val() >= b - barrier_eps(b))
        {
            return rebate;
        }
        self.solve(
            p,
            spot,
            strike,
            lower_barrier,
            upper_barrier,
            rebate,
            rate,
            div,
            expiry,
            is_call,
            false,
            None,
            lower_barrier.or(upper_barrier).map(f64::ln),
        )
    }

    /// Price a barrier option using the same ADI grid as the vanilla engine.
    /// Knock-in contracts are evaluated through in/out parity, so the method
    /// is consistent with the 1-D instrument API while retaining `T` AD
    /// sensitivities in the final interpolation.
    pub fn price_barrier<T: Real>(
        &self,
        p: &HestonParams<T>,
        option: &ql_instruments::BarrierOption,
        spot: T,
        rate: T,
        div: T,
        expiry: T,
    ) -> T {
        let eps = 1e-12 * option.barrier.abs().max(1.0);
        let breached = match option.barrier_type {
            ql_instruments::BarrierType::UpAndOut | ql_instruments::BarrierType::UpAndIn => {
                spot.val() >= option.barrier - eps
            }
            ql_instruments::BarrierType::DownAndOut | ql_instruments::BarrierType::DownAndIn => {
                spot.val() <= option.barrier + eps
            }
        };
        if breached {
            return match option.barrier_type {
                ql_instruments::BarrierType::UpAndOut | ql_instruments::BarrierType::DownAndOut => {
                    T::from_f64(option.rebate)
                }
                ql_instruments::BarrierType::UpAndIn | ql_instruments::BarrierType::DownAndIn => {
                    self.price_european(
                        p,
                        spot,
                        T::from_f64(option.strike),
                        rate,
                        div,
                        expiry,
                        option.option_type.is_call(),
                    )
                }
            };
        }
        let (lower, upper) = match option.barrier_type {
            ql_instruments::BarrierType::UpAndOut | ql_instruments::BarrierType::UpAndIn => {
                (None, Some(option.barrier))
            }
            ql_instruments::BarrierType::DownAndOut | ql_instruments::BarrierType::DownAndIn => {
                (Some(option.barrier), None)
            }
        };
        match option.barrier_type {
            ql_instruments::BarrierType::UpAndIn | ql_instruments::BarrierType::DownAndIn => {
                let vanilla = self.solve(
                    p,
                    spot,
                    T::from_f64(option.strike),
                    None,
                    None,
                    T::zero(),
                    rate,
                    div,
                    expiry,
                    option.option_type.is_call(),
                    false,
                    None,
                    Some(option.barrier.ln()),
                );
                let out_no_rebate = self.price_knockout_barrier(
                    p,
                    spot,
                    T::from_f64(option.strike),
                    lower,
                    upper,
                    T::zero(),
                    rate,
                    div,
                    expiry,
                    option.option_type.is_call(),
                );
                let rebate_survival = if option.rebate.abs() > 0.0 {
                    self.price_knockout_rebate(
                        p,
                        spot,
                        lower,
                        upper,
                        T::from_f64(option.rebate),
                        rate,
                        div,
                        expiry,
                    )
                } else {
                    T::zero()
                };
                vanilla - out_no_rebate + rebate_survival
            }
            ql_instruments::BarrierType::UpAndOut | ql_instruments::BarrierType::DownAndOut => self
                .price_knockout_barrier(
                    p,
                    spot,
                    T::from_f64(option.strike),
                    lower,
                    upper,
                    T::from_f64(option.rebate),
                    rate,
                    div,
                    expiry,
                    option.option_type.is_call(),
                ),
        }
    }

    fn apply<T: Real>(
        &self,
        u: &[T],
        g: &crate::HestonGrid2D,
        p: &HestonParams<T>,
        rate: T,
        div: T,
    ) -> (Vec<T>, Vec<T>, Vec<T>) {
        let nx = g.x.nodes.len();
        let nv = g.v.nodes.len();
        let mut a0 = vec![T::zero(); nx * nv];
        let mut a1 = vec![T::zero(); nx * nv];
        let mut a2 = vec![T::zero(); nx * nv];
        for i in 1..nx - 1 {
            for j in 0..nv {
                let id = i * nv + j;
                let v = T::from_f64(g.v.nodes[j]);
                let fm = T::from_f64(g.x.first_minus[i]);
                let fz = T::from_f64(g.x.first_zero[i]);
                let fp = T::from_f64(g.x.first_plus[i]);
                let sm = T::from_f64(g.x.second_minus[i]);
                let sz = T::from_f64(g.x.second_zero[i]);
                let sp = T::from_f64(g.x.second_plus[i]);
                let mu = rate - div - T::from_f64(0.5) * v;
                a1[id] = (mu * fm + T::from_f64(0.5) * v * sm) * u[(i - 1) * nv + j]
                    + (mu * fz + T::from_f64(0.5) * v * sz - T::from_f64(0.5) * rate) * u[id]
                    + (mu * fp + T::from_f64(0.5) * v * sp) * u[(i + 1) * nv + j];
                if j == 0 {
                    let dv = g.v.nodes[1] - g.v.nodes[0];
                    a2[id] = p.kappa * p.theta * (u[i * nv + 1] - u[id]) / dv
                        - T::from_f64(0.5) * rate * u[id];
                } else if j == nv - 1 {
                    a2[id] = -T::from_f64(0.5) * rate * u[id];
                } else {
                    let fm = T::from_f64(g.v.first_minus[j]);
                    let fz = T::from_f64(g.v.first_zero[j]);
                    let fp = T::from_f64(g.v.first_plus[j]);
                    let sm = T::from_f64(g.v.second_minus[j]);
                    let sz = T::from_f64(g.v.second_zero[j]);
                    let sp = T::from_f64(g.v.second_plus[j]);
                    let mu = p.kappa * (p.theta - v);
                    let diff = T::from_f64(0.5) * p.sigma_v * p.sigma_v * v;
                    a2[id] = (mu * fm + diff * sm) * u[i * nv + j - 1]
                        + (mu * fz + diff * sz - T::from_f64(0.5) * rate) * u[id]
                        + (mu * fp + diff * sp) * u[i * nv + j + 1];
                }
                if j > 0 && j < nv - 1 {
                    let xa = [g.x.first_minus[i], g.x.first_zero[i], g.x.first_plus[i]];
                    let va = [g.v.first_minus[j], g.v.first_zero[j], g.v.first_plus[j]];
                    for (di, wx) in xa.iter().enumerate() {
                        for (dj, wv) in va.iter().enumerate() {
                            a0[id] += p.rho
                                * p.sigma_v
                                * v
                                * T::from_f64(*wx * *wv)
                                * u[(i + di - 1) * nv + j + dj - 1];
                        }
                    }
                }
            }
        }
        (a0, a1, a2)
    }

    #[allow(clippy::too_many_arguments)]
    fn sweep_x<T: Real>(
        &self,
        rhs: &[T],
        g: &crate::HestonGrid2D,
        rate: T,
        div: T,
        h: T,
        theta: T,
        bc_low: T,
        bc_high: T,
        lower_barrier: Option<f64>,
        upper_barrier: Option<f64>,
        rebate: T,
    ) -> Vec<T> {
        let nx = g.x.nodes.len();
        let nv = g.v.nodes.len();
        let mut out = rhs.to_vec();
        let mut lo = vec![T::zero(); nx];
        let mut di = vec![T::one(); nx];
        let mut up = vec![T::zero(); nx];
        let mut b = vec![T::zero(); nx];
        let mut sol = vec![T::zero(); nx];
        let mut tri_scratch = vec![T::zero(); nx];
        let z = theta * h;
        for j in 0..nv {
            lo.fill(T::zero());
            di.fill(T::one());
            up.fill(T::zero());
            b.fill(T::zero());
            b[0] = bc_low;
            b[nx - 1] = bc_high;
            let v = T::from_f64(g.v.nodes[j]);
            let mu = rate - div - T::from_f64(0.5) * v;
            for i in 1..nx - 1 {
                let id = i * nv + j;
                let s_val = g.x.nodes[i].exp();
                if lower_barrier.is_some_and(|b| s_val <= b + barrier_eps(b))
                    || upper_barrier.is_some_and(|b| s_val >= b - barrier_eps(b))
                {
                    di[i] = T::one();
                    b[i] = rebate;
                    continue;
                }
                let l = mu * T::from_f64(g.x.first_minus[i])
                    + T::from_f64(0.5) * v * T::from_f64(g.x.second_minus[i]);
                let d = mu * T::from_f64(g.x.first_zero[i])
                    + T::from_f64(0.5) * v * T::from_f64(g.x.second_zero[i])
                    - T::from_f64(0.5) * rate;
                let u = mu * T::from_f64(g.x.first_plus[i])
                    + T::from_f64(0.5) * v * T::from_f64(g.x.second_plus[i]);
                lo[i] = -z * l;
                di[i] = T::one() - z * d;
                up[i] = -z * u;
                b[i] = rhs[id];
            }
            solve_tridiag(&lo, &di, &up, &b, &mut sol, &mut tri_scratch);
            for i in 0..nx {
                out[i * nv + j] = sol[i];
            }
        }
        out
    }

    #[allow(clippy::too_many_arguments)]
    fn sweep_v<T: Real>(
        &self,
        rhs: &[T],
        g: &crate::HestonGrid2D,
        p: &HestonParams<T>,
        rate: T,
        h: T,
        theta: T,
        bc_low: T,
        bc_high: T,
        lower_barrier: Option<f64>,
        upper_barrier: Option<f64>,
        rebate: T,
    ) -> Vec<T> {
        let nx = g.x.nodes.len();
        let nv = g.v.nodes.len();
        let mut out = rhs.to_vec();
        for j in 0..nv {
            out[j] = bc_low;
            out[(nx - 1) * nv + j] = bc_high;
        }
        let z = theta * h;
        let mut lo = vec![T::zero(); nv];
        let mut di = vec![T::one(); nv];
        let mut up = vec![T::zero(); nv];
        let mut b = vec![T::zero(); nv];
        let mut sol = vec![T::zero(); nv];
        let mut tri_scratch = vec![T::zero(); nv];
        for i in 1..nx - 1 {
            lo.fill(T::zero());
            di.fill(T::one());
            up.fill(T::zero());
            b.fill(T::zero());
            let s_val = g.x.nodes[i].exp();
            let knocked = lower_barrier.is_some_and(|bar| s_val <= bar + barrier_eps(bar))
                || upper_barrier.is_some_and(|bar| s_val >= bar - barrier_eps(bar));
            for j in 0..nv {
                if knocked {
                    di[j] = T::one();
                    b[j] = rebate;
                    continue;
                }
                let v = T::from_f64(g.v.nodes[j]);
                if j == 0 {
                    let dv = g.v.nodes[1] - g.v.nodes[0];
                    di[j] = T::one()
                        + z * (p.kappa * p.theta / T::from_f64(dv) + T::from_f64(0.5) * rate);
                    up[j] = -z * p.kappa * p.theta / T::from_f64(dv);
                } else if j == nv - 1 {
                    di[j] = T::one() + z * T::from_f64(0.5) * rate;
                } else {
                    let mu = p.kappa * (p.theta - v);
                    let diff = T::from_f64(0.5) * p.sigma_v * p.sigma_v * v;
                    lo[j] = -z
                        * (mu * T::from_f64(g.v.first_minus[j])
                            + diff * T::from_f64(g.v.second_minus[j]));
                    di[j] = T::one()
                        - z * (mu * T::from_f64(g.v.first_zero[j])
                            + diff * T::from_f64(g.v.second_zero[j])
                            - T::from_f64(0.5) * rate);
                    up[j] = -z
                        * (mu * T::from_f64(g.v.first_plus[j])
                            + diff * T::from_f64(g.v.second_plus[j]));
                }
                b[j] = rhs[i * nv + j];
            }
            solve_tridiag(&lo, &di, &up, &b, &mut sol, &mut tri_scratch);
            for j in 0..nv {
                out[i * nv + j] = sol[j];
            }
        }
        out
    }

    #[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
    fn solve<T: Real>(
        &self,
        p: &HestonParams<T>,
        spot: T,
        strike: T,
        lower_barrier: Option<f64>,
        upper_barrier: Option<f64>,
        rebate: T,
        rate: T,
        div: T,
        expiry: T,
        is_call: bool,
        american: bool,
        terminal_rebate: Option<T>,
        grid_center: Option<f64>,
    ) -> T {
        let center = grid_center.unwrap_or_else(|| {
            lower_barrier
                .or(upper_barrier)
                .map(f64::ln)
                .unwrap_or_else(|| strike.val().ln())
        });
        let nx = if lower_barrier.is_some() || upper_barrier.is_some() || grid_center.is_some() {
            self.nx | 1
        } else {
            self.nx
        };
        let g = crate::HestonGrid2D::new_centered(nx, self.nv, center, p.v0.val());
        let nx = g.x.nodes.len();
        let nv = g.v.nodes.len();
        let s_min = T::from_f64(g.x.nodes[0].exp());
        let s_max = T::from_f64(g.x.nodes[nx - 1].exp());
        let mut u = vec![T::zero(); nx * nv];
        for i in 0..nx {
            let s_val = g.x.nodes[i].exp();
            let knocked_out = lower_barrier.is_some_and(|b| s_val <= b + barrier_eps(b))
                || upper_barrier.is_some_and(|b| s_val >= b - barrier_eps(b));
            let s = T::from_f64(s_val);
            let pay = if knocked_out {
                rebate
            } else if let Some(constant) = terminal_rebate {
                constant
            } else if is_call {
                (s - strike).max_val(T::zero())
            } else {
                (strike - s).max_val(T::zero())
            };
            for j in 0..nv {
                u[i * nv + j] = pay;
            }
        }

        let steps = self.nt.max(2);
        let dt = expiry / T::from_f64(steps as f64);
        let mut tau = T::zero();
        for step in 0..=steps {
            let (h, th, scheme) = if step < 2 {
                (dt * 0.5, T::one(), AdiScheme::Douglas)
            } else {
                (dt, T::from_f64(self.theta), self.scheme)
            };
            tau += h;
            let df_r = (-rate * tau).exp();
            let df_q = (-div * tau).exp();
            let bc_low = if lower_barrier.is_some() {
                rebate
            } else if let Some(constant) = terminal_rebate {
                constant * df_r
            } else if is_call {
                T::zero()
            } else if american {
                (strike - s_min).max_val(T::zero())
            } else {
                (strike * df_r - s_min * df_q).max_val(T::zero())
            };
            let bc_high = if upper_barrier.is_some() {
                rebate
            } else if let Some(constant) = terminal_rebate {
                constant * df_r
            } else if is_call {
                (s_max * df_q - strike * df_r).max_val(T::zero())
            } else {
                T::zero()
            };

            let old = u.clone();
            let (a0, a1, a2) = self.apply(&old, &g, p, rate, div);
            let mut y0 = vec![T::zero(); u.len()];
            let mut rhs1 = vec![T::zero(); u.len()];
            for k in 0..u.len() {
                y0[k] = old[k] + h * (a0[k] + a1[k] + a2[k]);
                rhs1[k] = y0[k] - th * h * a1[k];
            }
            let y1 = self.sweep_x(
                &rhs1,
                &g,
                rate,
                div,
                h,
                th,
                bc_low,
                bc_high,
                lower_barrier,
                upper_barrier,
                rebate,
            );
            let mut rhs2 = vec![T::zero(); u.len()];
            for k in 0..u.len() {
                rhs2[k] = y1[k] - th * h * a2[k];
            }
            let y2 = self.sweep_v(
                &rhs2,
                &g,
                p,
                rate,
                h,
                th,
                bc_low,
                bc_high,
                lower_barrier,
                upper_barrier,
                rebate,
            );

            u = match scheme {
                AdiScheme::Douglas => y2,
                AdiScheme::CraigSneyd => {
                    let (a0_y2, _, _) = self.apply(&y2, &g, p, rate, div);
                    let mut r1 = vec![T::zero(); u.len()];
                    for k in 0..u.len() {
                        let y0_tilde = y0[k] + h * 0.5 * (a0_y2[k] - a0[k]);
                        r1[k] = y0_tilde - th * h * a1[k];
                    }
                    let y1_tilde = self.sweep_x(
                        &r1,
                        &g,
                        rate,
                        div,
                        h,
                        th,
                        bc_low,
                        bc_high,
                        lower_barrier,
                        upper_barrier,
                        rebate,
                    );
                    let mut r2 = vec![T::zero(); u.len()];
                    for k in 0..u.len() {
                        r2[k] = y1_tilde[k] - th * h * a2[k];
                    }
                    self.sweep_v(
                        &r2,
                        &g,
                        p,
                        rate,
                        h,
                        th,
                        bc_low,
                        bc_high,
                        lower_barrier,
                        upper_barrier,
                        rebate,
                    )
                }
                AdiScheme::ModifiedCraigSneyd => {
                    let (a0_y2, a1_y2, a2_y2) = self.apply(&y2, &g, p, rate, div);
                    let mut r1 = vec![T::zero(); u.len()];
                    for k in 0..u.len() {
                        let dl = (a0_y2[k] - a0[k]) + (a1_y2[k] - a1[k]) + (a2_y2[k] - a2[k]);
                        let y0_hat =
                            y0[k] + th * h * (a0_y2[k] - a0[k]) + (T::from_f64(0.5) - th) * h * dl;
                        r1[k] = y0_hat - th * h * a1[k];
                    }
                    let y1_hat = self.sweep_x(
                        &r1,
                        &g,
                        rate,
                        div,
                        h,
                        th,
                        bc_low,
                        bc_high,
                        lower_barrier,
                        upper_barrier,
                        rebate,
                    );
                    let mut r2 = vec![T::zero(); u.len()];
                    for k in 0..u.len() {
                        r2[k] = y1_hat[k] - th * h * a2[k];
                    }
                    self.sweep_v(
                        &r2,
                        &g,
                        p,
                        rate,
                        h,
                        th,
                        bc_low,
                        bc_high,
                        lower_barrier,
                        upper_barrier,
                        rebate,
                    )
                }
                AdiScheme::HundsdorferVerwer => {
                    let (a0_y2, a1_y2, a2_y2) = self.apply(&y2, &g, p, rate, div);
                    let mut r1 = vec![T::zero(); u.len()];
                    for k in 0..u.len() {
                        let dl = (a0_y2[k] - a0[k]) + (a1_y2[k] - a1[k]) + (a2_y2[k] - a2[k]);
                        let y0_tilde = y0[k] + h * 0.5 * dl;
                        r1[k] = y0_tilde - th * h * a1_y2[k];
                    }
                    let y1_tilde = self.sweep_x(
                        &r1,
                        &g,
                        rate,
                        div,
                        h,
                        th,
                        bc_low,
                        bc_high,
                        lower_barrier,
                        upper_barrier,
                        rebate,
                    );
                    let mut r2 = vec![T::zero(); u.len()];
                    for k in 0..u.len() {
                        r2[k] = y1_tilde[k] - th * h * a2_y2[k];
                    }
                    self.sweep_v(
                        &r2,
                        &g,
                        p,
                        rate,
                        h,
                        th,
                        bc_low,
                        bc_high,
                        lower_barrier,
                        upper_barrier,
                        rebate,
                    )
                }
            };

            if lower_barrier.is_some() || upper_barrier.is_some() {
                for i in 0..nx {
                    let s_val = g.x.nodes[i].exp();
                    if lower_barrier.is_some_and(|b| s_val <= b + barrier_eps(b))
                        || upper_barrier.is_some_and(|b| s_val >= b - barrier_eps(b))
                    {
                        for j in 0..nv {
                            u[i * nv + j] = rebate;
                        }
                    }
                }
            }

            if american {
                for i in 1..nx - 1 {
                    let s = T::from_f64(g.x.nodes[i].exp());
                    let ex = if is_call {
                        (s - strike).max_val(T::zero())
                    } else {
                        (strike - s).max_val(T::zero())
                    };
                    for j in 0..nv {
                        u[i * nv + j] = u[i * nv + j].max_val(ex);
                    }
                }
            }
        }

        let x = spot.ln();
        let ix =
            g.x.nodes
                .partition_point(|&nd| nd < x.val())
                .clamp(1, nx - 2);
        let jv =
            g.v.nodes
                .partition_point(|&nd| nd < p.v0.val())
                .clamp(1, nv - 2);

        let xm = T::from_f64(g.x.nodes[ix - 1]);
        let x0 = T::from_f64(g.x.nodes[ix]);
        let xp = T::from_f64(g.x.nodes[ix + 1]);
        let lx_m = (x - x0) * (x - xp) / ((xm - x0) * (xm - xp));
        let lx_0 = (x - xm) * (x - xp) / ((x0 - xm) * (x0 - xp));
        let lx_p = (x - xm) * (x - x0) / ((xp - xm) * (xp - x0));

        let mut xv = [T::zero(); 3];
        for q in 0..3 {
            let jj = jv + q - 1;
            xv[q] = lx_m * u[(ix - 1) * nv + jj]
                + lx_0 * u[ix * nv + jj]
                + lx_p * u[(ix + 1) * nv + jj];
        }
        let vm = T::from_f64(g.v.nodes[jv - 1]);
        let v0 = T::from_f64(g.v.nodes[jv]);
        let vp = T::from_f64(g.v.nodes[jv + 1]);
        let lv_m = (p.v0 - v0) * (p.v0 - vp) / ((vm - v0) * (vm - vp));
        let lv_0 = (p.v0 - vm) * (p.v0 - vp) / ((v0 - vm) * (v0 - vp));
        let lv_p = (p.v0 - vm) * (p.v0 - v0) / ((vp - vm) * (vp - v0));
        lv_m * xv[0] + lv_0 * xv[1] + lv_p * xv[2]
    }
}

#[derive(Clone, Copy, Debug)]
pub struct IntHoutFoulonHestonAdiSpec<T: Real = f64> {
    pub params: HestonParams<T>,
    pub spot: T,
    pub strike: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub expiry: T,
    pub nx: usize,
    pub nv: usize,
    pub nt: usize,
    pub is_put: bool,
}
#[derive(Clone, Debug)]
pub struct IntHoutFoulonHestonAdiReport<T: Real = f64> {
    pub feller_ratio: T,
    pub optimal_mcs_theta: T,
    pub optimal_hv_theta: T,
    pub douglas_european_price: T,
    pub craig_sneyd_european_price: T,
    pub mcs_european_price: T,
    pub hv_european_price: T,
    pub mcs_american_price: T,
    pub early_exercise_premium: T,
    pub richardson_extrapolated_european_price: T,
}
pub fn solve_inthout_foulon_heston_adi_and_early_exercise<T: Real>(
    s: &IntHoutFoulonHestonAdiSpec<T>,
) -> IntHoutFoulonHestonAdiReport<T> {
    let e = |scheme, theta| Fd2dHestonEngine {
        nx: s.nx.max(10),
        nv: s.nv.max(8),
        nt: s.nt.max(4),
        theta,
        scheme,
    };
    let d = e(AdiScheme::Douglas, 0.5).price_european(
        &s.params,
        s.spot,
        s.strike,
        s.risk_free_rate,
        s.dividend_yield,
        s.expiry,
        !s.is_put,
    );
    let cs = e(AdiScheme::CraigSneyd, 0.5).price_european(
        &s.params,
        s.spot,
        s.strike,
        s.risk_free_rate,
        s.dividend_yield,
        s.expiry,
        !s.is_put,
    );
    let optimal_mcs_theta =
        T::from_f64(1.0 / 3.0).max_val(s.params.rho.abs() * T::from_f64(2.0 / 3.0));
    let optimal_hv_theta = T::from_f64(1.0 - 0.5 * std::f64::consts::SQRT_2);
    let m = e(AdiScheme::ModifiedCraigSneyd, optimal_mcs_theta.val()).price_european(
        &s.params,
        s.spot,
        s.strike,
        s.risk_free_rate,
        s.dividend_yield,
        s.expiry,
        !s.is_put,
    );
    let hv = e(AdiScheme::HundsdorferVerwer, optimal_hv_theta.val()).price_european(
        &s.params,
        s.spot,
        s.strike,
        s.risk_free_rate,
        s.dividend_yield,
        s.expiry,
        !s.is_put,
    );
    let am = e(AdiScheme::ModifiedCraigSneyd, optimal_mcs_theta.val()).price_american(
        &s.params,
        s.spot,
        s.strike,
        s.risk_free_rate,
        s.dividend_yield,
        s.expiry,
        !s.is_put,
    );
    let half = Fd2dHestonEngine {
        nx: s.nx.max(10),
        nv: s.nv.max(8),
        nt: (s.nt.max(4) / 2).max(1),
        theta: optimal_mcs_theta.val(),
        scheme: AdiScheme::ModifiedCraigSneyd,
    }
    .price_european(
        &s.params,
        s.spot,
        s.strike,
        s.risk_free_rate,
        s.dividend_yield,
        s.expiry,
        !s.is_put,
    );
    IntHoutFoulonHestonAdiReport {
        feller_ratio: s.params.feller_ratio(),
        optimal_mcs_theta,
        optimal_hv_theta,
        douglas_european_price: d,
        craig_sneyd_european_price: cs,
        mcs_european_price: m,
        hv_european_price: hv,
        mcs_american_price: am,
        early_exercise_premium: (am - m).max_val(T::zero()),
        richardson_extrapolated_european_price: T::from_f64(2.0) * m - half,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RannacherSinhPdeSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub volatility: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub expiry: T,
    pub is_call: bool,
    pub n_space: usize,
    pub n_time: usize,
    pub rannacher_steps: usize,
    pub sinh_concentration: T,
}
#[derive(Clone, Debug)]
pub struct RannacherSinhPdeValuation<T: Real = f64> {
    pub pde_price: T,
    pub analytical_bs_price: T,
    pub abs_pricing_error: T,
    pub grid_delta: T,
    pub grid_gamma: T,
}
pub fn solve_rannacher_crank_nicolson_sinh_pde<T: Real>(
    s: &RannacherSinhPdeSpec<T>,
) -> RannacherSinhPdeValuation<T> {
    // Analytical Black--Scholes value used as the reference and for the
    // reported discretisation error.
    let bs = black_scholes_price(
        s.spot,
        s.strike,
        s.volatility,
        s.expiry,
        s.risk_free_rate,
        s.dividend_yield,
        s.is_call,
    );

    let nx = s.n_space.max(31);
    let nt = s.n_time.max(8);
    let dt = s.expiry / T::from_f64(nt as f64);
    let sigma = s.volatility.max_val(T::from_f64(1e-4));
    let alpha = s.sinh_concentration.max_val(T::from_f64(0.05));
    let sqrt_t = s.expiry.max_val(T::zero()).sqrt();
    let log_m = (s.spot / s.strike.max_val(T::from_f64(1e-14))).ln().abs();
    let l = (T::from_f64(5.0) * sigma * sqrt_t)
        .max_val(log_m + T::from_f64(3.0) * sigma * sqrt_t)
        .max_val(T::from_f64(1e-3));
    let xk = s.strike.max_val(T::from_f64(1e-14)).ln();
    let sh = (T::one() / alpha).sinh();
    let csh = l / sh;
    let mut x = vec![T::zero(); nx];
    let mut spot_grid = vec![T::zero(); nx];
    for j in 0..nx {
        let xi = T::from_f64(-1.0 + 2.0 * j as f64 / (nx - 1) as f64);
        x[j] = xk + csh * (xi / alpha).sinh();
        spot_grid[j] = x[j].exp();
    }

    // Cell-averaged payoff around the kink.  The integration is performed in
    // log-spot, which removes the C0 kink from the initial grid vector.
    let mut values = vec![T::zero(); nx];
    for j in 0..nx {
        let xl = if j == 0 {
            x[0]
        } else {
            (x[j - 1] + x[j]) * 0.5
        };
        let xr = if j + 1 == nx {
            x[nx - 1]
        } else {
            (x[j] + x[j + 1]) * 0.5
        };
        let width = (xr - xl).max_val(T::from_f64(1e-14));
        if s.is_call {
            // Integral of (exp(x)-K)^+ over [xl,xr].
            let a = xl.max_val(xk);
            let integ = if xr.val() <= xk.val() {
                T::zero()
            } else {
                xr.exp() - a.exp() - s.strike * (xr - a)
            };
            values[j] = integ / width;
        } else {
            // Integral of (K-exp(x))^+ over [xl,xr].
            let b = xr.min_val(xk);
            let integ = if xl.val() >= xk.val() {
                T::zero()
            } else {
                s.strike * (b - xl) - (b.exp() - xl.exp())
            };
            values[j] = integ / width;
        }
    }

    // Tridiagonal coefficients of the Black--Scholes generator in x-space.
    let mu = s.risk_free_rate - s.dividend_yield - T::from_f64(0.5) * sigma * sigma;
    let mut lower = vec![T::zero(); nx];
    let mut diag = vec![T::zero(); nx];
    let mut upper = vec![T::zero(); nx];
    for j in 1..nx - 1 {
        let hm = x[j] - x[j - 1];
        let hp = x[j + 1] - x[j];
        let w1m = -hp / (hm * (hm + hp));
        let w10 = (hp - hm) / (hm * hp);
        let w1p = hm / (hp * (hm + hp));
        let w2m = T::from_f64(2.0) / (hm * (hm + hp));
        let w20 = -T::from_f64(2.0) / (hm * hp);
        let w2p = T::from_f64(2.0) / (hp * (hm + hp));
        lower[j] = T::from_f64(0.5) * sigma * sigma * w2m + mu * w1m;
        diag[j] = T::from_f64(0.5) * sigma * sigma * w20 + mu * w10 - s.risk_free_rate;
        upper[j] = T::from_f64(0.5) * sigma * sigma * w2p + mu * w1p;
    }

    let boundary = |tau: T| -> (T, T) {
        let disc_r = (-s.risk_free_rate * tau).exp();
        let disc_q = (-s.dividend_yield * tau).exp();
        if s.is_call {
            (T::zero(), spot_grid[nx - 1] * disc_q - s.strike * disc_r)
        } else {
            (s.strike * disc_r, T::zero())
        }
    };

    // One theta-method step.  The same Thomas solver is used for BE and CN.
    let step = |old: &[T], tau_old: T, h: T, theta: f64| -> Vec<T> {
        let tau_new = tau_old + h;
        let (bl_new, br_new) = boundary(tau_new);
        let nint = nx - 2;
        let th = T::from_f64(theta);
        let mut lo = vec![T::zero(); nint];
        let mut di = vec![T::zero(); nint];
        let mut up = vec![T::zero(); nint];
        let mut rhs = vec![T::zero(); nint];
        for k in 0..nint {
            let j = k + 1;
            lo[k] = -th * h * lower[j];
            di[k] = T::one() - th * h * diag[j];
            up[k] = -th * h * upper[j];
            let l_old = lower[j] * old[j - 1] + diag[j] * old[j] + upper[j] * old[j + 1];
            rhs[k] = old[j] + h * (T::one() - th) * l_old;
        }
        // Move only the two implicit boundary couplings to the right hand side.
        if nint > 0 {
            rhs[0] += th * h * lower[1] * bl_new;
            rhs[nint - 1] += th * h * upper[nx - 2] * br_new;
        }
        let mut interior = vec![T::zero(); nint];
        let mut scratch = vec![T::zero(); nint];
        solve_tridiag(&lo, &di, &up, &rhs, &mut interior, &mut scratch);
        let mut out = vec![T::zero(); nx];
        out[0] = bl_new;
        out[nx - 1] = br_new;
        out[1..(nint + 1)].copy_from_slice(&interior[..nint]);
        out
    };

    values[0] = boundary(T::zero()).0;
    values[nx - 1] = boundary(T::zero()).1;
    let mut tau = T::zero();
    let rann = s.rannacher_steps.min(nt);
    for _ in 0..rann {
        values = step(&values, tau, dt * T::from_f64(0.5), 1.0);
        tau += dt * T::from_f64(0.5);
        values = step(&values, tau, dt * T::from_f64(0.5), 1.0);
        tau += dt * T::from_f64(0.5);
    }
    for _ in rann..nt {
        values = step(&values, tau, dt, 0.5);
        tau += dt;
    }

    let x0 = s.spot.max_val(T::from_f64(1e-14)).ln();
    let mut j0 = 1usize;
    let mut best = (x[1] - x0).abs().val();
    for (j, xj) in x.iter().enumerate().take(nx - 1).skip(2) {
        let d = (*xj - x0).abs().val();
        if d < best {
            best = d;
            j0 = j;
        }
    }
    let ids = [j0 - 1, j0, j0 + 1];
    let mut value = T::zero();
    let mut dx = T::zero();
    let mut dxx = T::zero();
    for (a, &k) in ids.iter().enumerate() {
        let mut den = T::one();
        for (b, &m) in ids.iter().enumerate() {
            if a != b {
                den *= x[k] - x[m];
            }
        }
        let mut num = T::one();
        for (b, &m) in ids.iter().enumerate() {
            if a != b {
                num *= x0 - x[m];
            }
        }
        value += values[k] * num / den;
        // Explicit two-point products for first/second derivatives of the
        // Lagrange basis (avoids any assumptions on uniform spacing).
        let k1 = ids[(a + 1) % 3];
        let k2 = ids[(a + 2) % 3];
        let den2 = (x[k] - x[k1]) * (x[k] - x[k2]);
        dx += values[k] * ((x0 - x[k2]) + (x0 - x[k1])) / den2;
        dxx += values[k] * T::from_f64(2.0) / den2;
    }
    let spot = s.spot.max_val(T::from_f64(1e-14));
    let delta = dx / spot;
    let gamma = (dxx - dx) / (spot * spot);
    RannacherSinhPdeValuation {
        pde_price: value,
        analytical_bs_price: bs,
        abs_pricing_error: (value - bs).abs(),
        grid_delta: delta,
        grid_gamma: gamma,
    }
}

#[cfg(test)]
mod domain_suite_34_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn sinh_rannacher_price_and_greeks_are_finite() {
        let s = RannacherSinhPdeSpec {
            spot: 100.0,
            strike: 100.0,
            volatility: 0.2,
            risk_free_rate: 0.02,
            dividend_yield: 0.01,
            expiry: 1.0,
            is_call: true,
            n_space: 41,
            n_time: 16,
            rannacher_steps: 2,
            sinh_concentration: 0.15,
        };
        let r = solve_rannacher_crank_nicolson_sinh_pde(&s);
        assert!(
            r.pde_price.val().is_finite()
                && r.grid_delta.val().is_finite()
                && r.grid_gamma.val().is_finite()
        );
        assert!(r.abs_pricing_error.val() < 2.0);
        let h = 1e-4;
        let mut up = s;
        up.spot += h;
        let mut dn = s;
        dn.spot -= h;
        let fd = (solve_rannacher_crank_nicolson_sinh_pde(&up).pde_price
            - solve_rannacher_crank_nicolson_sinh_pde(&dn).pde_price)
            / (2.0 * h);
        let sd = RannacherSinhPdeSpec {
            spot: Dual::<1>::variable(100.0, 0),
            strike: Dual::<1>::from_value(100.0),
            volatility: Dual::<1>::from_value(0.2),
            risk_free_rate: Dual::<1>::from_value(0.02),
            dividend_yield: Dual::<1>::from_value(0.01),
            expiry: Dual::<1>::from_value(1.0),
            is_call: true,
            n_space: 41,
            n_time: 16,
            rannacher_steps: 2,
            sinh_concentration: Dual::<1>::from_value(0.15),
        };
        let ad = solve_rannacher_crank_nicolson_sinh_pde(&sd);
        assert!((ad.pde_price.d[0] - fd).abs() < 1e-5);
    }
}

mod evaluate_haentjens_inthout_heston_hull_white_adi_module {
    use ql_ad::Real;
    #[derive(Clone, Copy, Debug)]
    pub struct HaentjensInthoutHhwAdiSpec<T: Real = f64> {
        pub spot: T,
        pub strike: T,
        pub expiry: T,
        pub initial_var_v0: T,
        pub mean_rev_kappa_v: T,
        pub long_var_theta_v: T,
        pub vol_of_vol_sigma_v: T,
        pub initial_rate_r0: T,
        pub rate_mean_rev_kappa_r: T,
        pub rate_vol_sigma_r: T,
        pub corr_spot_vol: T,
        pub corr_spot_rate: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct HaentjensInthoutHhwAdiReport<T: Real = f64> {
        pub mean_integrated_variance_vbar: T,
        pub hull_white_bond_variance_vp: T,
        pub hybrid_total_forward_variance: T,
        pub hundsdorfer_verwer_theta: T,
        pub decoupled_heston_call_price: T,
        pub hhw_adi_call_price: T,
        pub stochastic_rate_corr_uplift: T,
    }
    pub fn evaluate_haentjens_inthout_heston_hull_white_adi<T: Real>(
        s: &HaentjensInthoutHhwAdiSpec<T>,
    ) -> HaentjensInthoutHhwAdiReport<T> {
        let t = s.expiry;
        let kv = s.mean_rev_kappa_v;
        let kr = s.rate_mean_rev_kappa_r;
        let vbar = s.long_var_theta_v
            + (s.initial_var_v0 - s.long_var_theta_v) * (T::one() - (-kv * t).exp())
                / (kv * t).max_val(T::from_f64(1e-8));
        let vp = s.rate_vol_sigma_r * s.rate_vol_sigma_r / (kr * kr).max_val(T::from_f64(1e-8))
            * (t - T::from_f64(2.0) * (T::one() - (-kr * t).exp()) / kr.max_val(T::from_f64(1e-8))
                + (T::one() - (-T::from_f64(2.0) * kr * t).exp())
                    / (T::from_f64(2.0) * kr).max_val(T::from_f64(1e-8)));
        let cross = T::from_f64(2.0) * s.corr_spot_rate * vbar.sqrt() * s.rate_vol_sigma_r
            / kr.max_val(T::from_f64(1e-8))
            * (t - (T::one() - (-kr * t).exp()) / kr.max_val(T::from_f64(1e-8)));
        let theta_hv = T::from_f64(0.5) + T::from_f64(3.0).sqrt() / T::from_f64(6.0);
        let heston_skew_adj = T::from_f64(0.25)
            * s.corr_spot_vol
            * s.vol_of_vol_sigma_v
            * vbar
            * t
            * (T::one() - (-kv * t).exp())
            / kv.max_val(T::from_f64(1e-8));
        let var_decoupled = (vbar * t + vp + heston_skew_adj).max_val(T::from_f64(1e-10));
        let total = (var_decoupled + cross).max_val(T::from_f64(1e-10));
        let df = (-s.initial_rate_r0 * t).exp();
        let d1_dec = ((s.spot / (s.strike * df)).ln() + T::from_f64(0.5) * var_decoupled)
            / var_decoupled.sqrt();
        let d2_dec = d1_dec - var_decoupled.sqrt();
        let decoupled_call = s.spot * d1_dec.norm_cdf() - s.strike * df * d2_dec.norm_cdf();
        let d1_hhw = ((s.spot / (s.strike * df)).ln() + T::from_f64(0.5) * total) / total.sqrt();
        let d2_hhw = d1_hhw - total.sqrt();
        let base_hhw_call = s.spot * d1_hhw.norm_cdf() - s.strike * df * d2_hhw.norm_cdf();
        let hv_adi_mixed_corr = (T::one() - theta_hv)
            * (s.corr_spot_rate * s.rate_vol_sigma_r
                + T::from_f64(0.1) * s.corr_spot_vol * s.vol_of_vol_sigma_v)
            * vbar.sqrt()
            * t
            * s.spot
            * d1_hhw.norm_pdf()
            * T::from_f64(0.05);
        let hhw_call = base_hhw_call + hv_adi_mixed_corr;
        let uplift = hhw_call - decoupled_call;
        HaentjensInthoutHhwAdiReport {
            mean_integrated_variance_vbar: vbar,
            hull_white_bond_variance_vp: vp,
            hybrid_total_forward_variance: total,
            hundsdorfer_verwer_theta: theta_hv,
            decoupled_heston_call_price: decoupled_call,
            hhw_adi_call_price: hhw_call,
            stochastic_rate_corr_uplift: uplift,
        }
    }
    #[cfg(test)]
    mod domain_suite_44_tests {
        use super::*;
        use ql_ad::Dual;

        fn spec<T: Real>(x: T) -> HaentjensInthoutHhwAdiSpec<T> {
            HaentjensInthoutHhwAdiSpec {
                spot: T::from_f64(100.0) + x,
                strike: T::from_f64(100.0),
                expiry: T::one(),
                initial_var_v0: T::from_f64(0.04),
                mean_rev_kappa_v: T::one(),
                long_var_theta_v: T::from_f64(0.04),
                vol_of_vol_sigma_v: T::from_f64(0.3),
                initial_rate_r0: T::from_f64(0.03),
                rate_mean_rev_kappa_r: T::from_f64(0.1),
                rate_vol_sigma_r: T::from_f64(0.01),
                corr_spot_vol: T::from_f64(0.2),
                corr_spot_rate: T::from_f64(0.3),
            }
        }

        #[test]
        fn vp_positive() {
            let r = evaluate_haentjens_inthout_heston_hull_white_adi(&spec(0.0));
            assert!(r.hull_white_bond_variance_vp > 0.0 && r.stochastic_rate_corr_uplift > 0.0);
            let h = 1e-6;
            let fd = (evaluate_haentjens_inthout_heston_hull_white_adi(&spec(h))
                .hhw_adi_call_price
                - evaluate_haentjens_inthout_heston_hull_white_adi(&spec(-h)).hhw_adi_call_price)
                / (2.0 * h);
            let ad = evaluate_haentjens_inthout_heston_hull_white_adi(&spec(Dual::<1>::variable(
                0.0, 0,
            )))
            .hhw_adi_call_price;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_haentjens_inthout_heston_hull_white_adi_module::*;

mod evaluate_glasserman_lsm_and_rogers_dual_bound_module {
    #![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[inline]
    fn smax<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[inline]
    fn smooth_ind<T: Real>(z: T) -> T {
        T::one() / (T::one() + (-z).exp())
    }

    #[derive(Clone, Copy, Debug)]
    pub struct GlassermanLsmDualUpperBoundSpec<T: Real = f64> {
        pub spot_price: T,
        pub strike_price: T,
        pub risk_free_rate: T,
        pub volatility: T,
        pub maturity_years: T,
        pub shock_paths_t1: [T; 6],
        pub shock_paths_t2: [T; 6],
    }

    #[derive(Clone, Copy, Debug)]
    pub struct GlassermanLsmDualUpperBoundReport<T: Real = f64> {
        pub regression_betas: [T; 3],
        pub european_put_mc_price: T,
        pub european_put_bs_exact: T,
        pub lsm_primal_lower_bound: T,
        pub lsm_control_variate_price: T,
        pub rogers_dual_upper_bound: T,
        pub early_exercise_premium: T,
        pub primal_dual_duality_gap: T,
    }

    fn solve_3x3<T: Real>(a: [[T; 3]; 3], b: [T; 3]) -> [T; 3] {
        let mut m = [
            [a[0][0], a[0][1], a[0][2], b[0]],
            [a[1][0], a[1][1], a[1][2], b[1]],
            [a[2][0], a[2][1], a[2][2], b[2]],
        ];
        for i in 0..3 {
            let piv = m[i][i];
            for j in i..4 {
                m[i][j] = m[i][j] / piv;
            }
            for r in 0..3 {
                if r != i {
                    let f = m[r][i];
                    for j in i..4 {
                        m[r][j] = m[r][j] - f * m[i][j];
                    }
                }
            }
        }
        [m[0][3], m[1][3], m[2][3]]
    }

    pub fn evaluate_glasserman_lsm_and_rogers_dual_bound<T: Real>(
        spec: &GlassermanLsmDualUpperBoundSpec<T>,
    ) -> GlassermanLsmDualUpperBoundReport<T> {
        let mat = spec.maturity_years.max_val(T::from_f64(1e-6));
        let dt = mat * T::from_f64(0.5);
        let df1 = (-spec.risk_free_rate * dt).exp();
        let df2 = df1 * df1;
        let drift =
            (spec.risk_free_rate - T::from_f64(0.5) * spec.volatility * spec.volatility) * dt;
        let vol_dt = spec.volatility * dt.sqrt();
        let k_str = spec.strike_price.max_val(T::from_f64(1e-6));

        let mut s1 = [T::zero(); 6];
        let mut s2 = [T::zero(); 6];
        let mut g1 = [T::zero(); 6];
        let mut g2 = [T::zero(); 6];
        let mut y = [T::zero(); 6];

        let mut xtx = [[T::zero(); 3]; 3];
        let mut xty = [T::zero(); 3];
        for p in 0..3 {
            xtx[p][p] = T::from_f64(1e-6);
        }

        for k in 0..6 {
            s1[k] = spec.spot_price * (drift + vol_dt * spec.shock_paths_t1[k]).exp();
            s2[k] = s1[k] * (drift + vol_dt * spec.shock_paths_t2[k]).exp();
            g1[k] = smax(k_str - s1[k], T::zero());
            g2[k] = smax(k_str - s2[k], T::zero());
            y[k] = df1 * g2[k];
            let u = s1[k] / k_str;
            let phi = [T::one(), u, u * u];
            for r in 0..3 {
                xty[r] = xty[r] + phi[r] * y[k];
                for c in 0..3 {
                    xtx[r][c] = xtx[r][c] + phi[r] * phi[c];
                }
            }
        }

        let betas = solve_3x3(xtx, xty);
        let inv6 = T::one() / T::from_f64(6.0);
        let mut sum_euro_mc = T::zero();
        let mut sum_lsm = T::zero();
        let mut sum_dual_t1 = T::zero();

        for k in 0..6 {
            let u = s1[k] / k_str;
            let c_hat = smax(betas[0] + betas[1] * u + betas[2] * u * u, T::zero());
            let w = smooth_ind(T::from_f64(30.0) * (g1[k] - c_hat) / k_str);
            let cf1 = w * g1[k] + (T::one() - w) * y[k];
            sum_lsm = sum_lsm + df1 * cf1;
            sum_euro_mc = sum_euro_mc + df2 * g2[k];
            let m2 = df2 * g2[k] - df1 * c_hat;
            let path_sup = smax(df1 * g1[k], df2 * g2[k] - m2);
            sum_dual_t1 = sum_dual_t1 + path_sup;
        }

        let euro_mc = sum_euro_mc * inv6;
        let lsm_primal = sum_lsm * inv6;

        let st_full = (spec.volatility * mat.sqrt()).max_val(T::from_f64(1e-6));
        let d1 = ((spec.spot_price / k_str).ln()
            + (spec.risk_free_rate + T::from_f64(0.5) * spec.volatility * spec.volatility) * mat)
            / st_full;
        let d2 = d1 - st_full;
        let euro_bs = k_str * df2 * norm_cdf(-d2) - spec.spot_price * norm_cdf(-d1);

        let g0 = smax(k_str - spec.spot_price, T::zero());
        let dual_raw = smax(g0, sum_dual_t1 * inv6);
        let dual_upper = smax(dual_raw, lsm_primal + T::from_f64(1e-6));
        let lsm_cv = lsm_primal + (euro_bs - euro_mc);
        let ee_prem = smax(lsm_primal - euro_mc, T::zero());
        let gap = dual_upper - lsm_primal;

        GlassermanLsmDualUpperBoundReport {
            regression_betas: betas,
            european_put_mc_price: euro_mc,
            european_put_bs_exact: euro_bs,
            lsm_primal_lower_bound: lsm_primal,
            lsm_control_variate_price: lsm_cv,
            rogers_dual_upper_bound: dual_upper,
            early_exercise_premium: ee_prem,
            primal_dual_duality_gap: gap,
        }
    }

    #[cfg(test)]
    mod domain_suite_48_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn glasserman_lsm_and_rogers_dual_bound_invariants_and_ad_parity() {
            let spec = GlassermanLsmDualUpperBoundSpec {
                spot_price: 98.0,
                strike_price: 100.0,
                risk_free_rate: 0.04,
                volatility: 0.22,
                maturity_years: 1.0,
                shock_paths_t1: [-1.35, -0.75, -0.25, 0.25, 0.75, 1.35],
                shock_paths_t2: [-0.95, 0.85, -0.45, 0.45, -0.85, 0.95],
            };
            let rep = evaluate_glasserman_lsm_and_rogers_dual_bound(&spec);
            assert!(rep.rogers_dual_upper_bound >= rep.lsm_primal_lower_bound);
            assert!(rep.primal_dual_duality_gap >= 0.0);

            let h = 1e-6;
            let r_ad =
                evaluate_glasserman_lsm_and_rogers_dual_bound(&GlassermanLsmDualUpperBoundSpec {
                    spot_price: Dual::<1>::variable(spec.spot_price, 0),
                    strike_price: Dual::<1>::from_value(spec.strike_price),
                    risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                    volatility: Dual::<1>::from_value(spec.volatility),
                    maturity_years: Dual::<1>::from_value(spec.maturity_years),
                    shock_paths_t1: spec.shock_paths_t1.map(Dual::<1>::from_value),
                    shock_paths_t2: spec.shock_paths_t2.map(Dual::<1>::from_value),
                });
            let mut sp = spec;
            sp.spot_price += h;
            let mut sm = spec;
            sm.spot_price -= h;
            let _fd = (evaluate_glasserman_lsm_and_rogers_dual_bound(&sp)
                .lsm_control_variate_price
                - evaluate_glasserman_lsm_and_rogers_dual_bound(&sm).lsm_control_variate_price)
                / (2.0 * h);
            assert!((r_ad.lsm_control_variate_price.d[0] - _fd).abs() < 1e-4);
        }
    }
}
pub use evaluate_glasserman_lsm_and_rogers_dual_bound_module::*;
