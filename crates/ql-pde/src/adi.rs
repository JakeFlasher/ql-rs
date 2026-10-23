use ql_ad::Real;
use ql_models::HestonParams;

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
        let (lower, upper) = match option.barrier_type {
            ql_instruments::BarrierType::UpAndOut | ql_instruments::BarrierType::UpAndIn => {
                (None, Some(option.barrier))
            }
            ql_instruments::BarrierType::DownAndOut | ql_instruments::BarrierType::DownAndIn => {
                (Some(option.barrier), None)
            }
        };
        let out = self.price_knockout_barrier(
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
        );
        match option.barrier_type {
            ql_instruments::BarrierType::UpAndIn | ql_instruments::BarrierType::DownAndIn => {
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
                self.price_european(
                    p,
                    spot,
                    T::from_f64(option.strike),
                    rate,
                    div,
                    expiry,
                    option.option_type.is_call(),
                ) - out_no_rebate
                    + (out - out_no_rebate)
            }
            ql_instruments::BarrierType::UpAndOut | ql_instruments::BarrierType::DownAndOut => out,
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
                if lower_barrier.is_some_and(|b| s_val <= b)
                    || upper_barrier.is_some_and(|b| s_val >= b)
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
            let knocked = lower_barrier.is_some_and(|bar| s_val <= bar)
                || upper_barrier.is_some_and(|bar| s_val >= bar);
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
    ) -> T {
        let center = lower_barrier
            .or(upper_barrier)
            .map(f64::ln)
            .unwrap_or_else(|| strike.val().ln());
        let g = crate::HestonGrid2D::new_centered(self.nx, self.nv, center, p.v0.val());
        let nx = g.x.nodes.len();
        let nv = g.v.nodes.len();
        let s_min = T::from_f64(g.x.nodes[0].exp());
        let s_max = T::from_f64(g.x.nodes[nx - 1].exp());
        let mut u = vec![T::zero(); nx * nv];
        for i in 0..nx {
            let s_val = g.x.nodes[i].exp();
            let knocked_out = lower_barrier.is_some_and(|b| s_val <= b)
                || upper_barrier.is_some_and(|b| s_val >= b);
            let s = T::from_f64(s_val);
            let pay = if knocked_out {
                rebate
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
            } else if is_call {
                T::zero()
            } else if american {
                (strike - s_min).max_val(T::zero())
            } else {
                (strike * df_r - s_min * df_q).max_val(T::zero())
            };
            let bc_high = if upper_barrier.is_some() {
                rebate
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
                    if lower_barrier.is_some_and(|b| s_val <= b)
                        || upper_barrier.is_some_and(|b| s_val >= b)
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
