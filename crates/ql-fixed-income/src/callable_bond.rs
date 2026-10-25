#![allow(clippy::excessive_precision, clippy::needless_range_loop)]
use ql_ad::Real;
use ql_vol::black76_price;
#[derive(Clone, Copy, Debug)]
pub struct BondCallability<T: Real> {
    pub exercise_time: T,
    pub clean_strike: T,
    pub is_call: bool,
}
#[derive(Clone, Debug)]
pub struct CallableFixedRateBond<T: Real> {
    pub face_amount: T,
    pub coupon_rate: T,
    pub coupon_times: Vec<T>,
    pub maturity: T,
    pub schedule: Vec<BondCallability<T>>,
}
impl<T: Real> CallableFixedRateBond<T> {
    /// Coupon amount at a schedule index.  `coupon_rate` is an annualized
    /// fixed rate, so irregular/semiannual schedules accrue over the interval
    /// since the preceding coupon date.
    fn coupon_amount(&self, i: usize) -> T {
        let end = self.coupon_times[i];
        let start = if i == 0 {
            T::zero()
        } else {
            self.coupon_times[i - 1]
        };
        self.face_amount * self.coupon_rate * (end - start)
    }

    pub fn straight_bond_npv(&self, rate: T) -> T {
        self.coupon_times
            .iter()
            .enumerate()
            .map(|(i, t)| {
                let cf = self.coupon_amount(i)
                    + if i + 1 == self.coupon_times.len() {
                        self.face_amount
                    } else {
                        T::zero()
                    };
                cf * (-rate * *t).exp()
            })
            .fold(T::zero(), |a, x| a + x)
    }

    fn cashflow_at(&self, i: usize) -> T {
        self.coupon_amount(i)
            + if i + 1 == self.coupon_times.len() {
                self.face_amount
            } else {
                T::zero()
            }
    }

    /// Forward clean price, yield and modified duration of the remaining
    /// bond cashflows at an exercise date under a flat yield `rate`.
    fn forward_bond_data(&self, exercise: T, rate: T) -> (T, T, T) {
        let mut price0 = T::zero();
        for (i, &t) in self.coupon_times.iter().enumerate() {
            if t.val() > exercise.val() {
                price0 += self.cashflow_at(i) * (-rate * (t - exercise)).exp();
            }
        }
        let mut y = rate;
        for _ in 0..24 {
            let mut p = T::zero();
            let mut dp = T::zero();
            for (i, &t) in self.coupon_times.iter().enumerate() {
                if t.val() > exercise.val() {
                    let tau = t - exercise;
                    let e = (-y * tau).exp();
                    p += self.cashflow_at(i) * e;
                    dp -= self.cashflow_at(i) * tau * e;
                }
            }
            if dp.val().abs() < 1e-14 {
                break;
            }
            let step = (p - price0) / dp;
            y -= step;
            if step.val().abs() < 1e-13 {
                break;
            }
        }
        let mut weighted = T::zero();
        let mut p = T::zero();
        for (i, &t) in self.coupon_times.iter().enumerate() {
            if t.val() > exercise.val() {
                let tau = t - exercise;
                let e = (-y * tau).exp();
                let cf = self.cashflow_at(i);
                p += cf * e;
                weighted += tau * cf * e;
            }
        }
        let duration = if p.val().abs() > 1e-14 {
            weighted / p
        } else {
            T::zero()
        };
        (price0, y, duration)
    }
    pub fn black_callable_bond_price(&self, rate: T, yield_vol: T) -> T {
        let p0 = self.straight_bond_npv(rate);
        let mut call_value = T::zero();
        let mut put_value = T::zero();
        for c in &self.schedule {
            let (forward, y_fwd, duration) = self.forward_bond_data(c.exercise_time, rate);
            let sigma_p = yield_vol * duration * y_fwd.abs();
            let option = black76_price(
                forward,
                c.clean_strike,
                sigma_p,
                c.exercise_time,
                (-rate * c.exercise_time).exp(),
                c.is_call,
            );
            if c.is_call {
                call_value = call_value.max_val(option);
            } else {
                put_value = put_value.max_val(option);
            }
        }
        p0 - call_value + put_value
    }
    pub fn tree_callable_bond_price(&self, short_rate: T, rate_vol: T, oas: T, steps: usize) -> T {
        let n = steps.max(1);
        let dt = self.maturity / n as f64;
        let vol_step = rate_vol * dt.sqrt();
        // Calibrate the deterministic lattice shift at every time level to
        // the initial flat discount curve.  For a flat curve this is zero,
        // but retaining it explicitly makes the node short rate a single
        // well-defined calibrated state rather than two ad-hoc edge rates.
        let theta_at = |j: usize| {
            let t0 = dt * j as f64;
            let t1 = dt * (j + 1) as f64;
            let p0 = (-(short_rate + oas) * t0).exp();
            let p1 = (-(short_rate + oas) * t1).exp();
            -(p1 / p0).ln() / dt - short_rate - oas
        };
        let mut values = vec![T::zero(); n + 1];
        // Principal is paid at maturity; the final coupon is added at the
        // maturity node just like every other coupon below.
        values.fill(self.face_amount);
        for (idx, &t) in self.coupon_times.iter().enumerate() {
            let step_idx = ((t / dt).val().round() as usize).min(n);
            if step_idx == n {
                let cf = self.cashflow_at(idx) - self.face_amount;
                values.fill(self.face_amount + cf);
            }
        }
        let exercise_at = |j: usize| -> Vec<&BondCallability<T>> {
            self.schedule
                .iter()
                .filter(|c| ((c.exercise_time / dt).val().round() as usize).min(n) == j)
                .collect()
        };
        for j in (0..n).rev() {
            let mut next = vec![T::zero(); j + 1];
            let theta_j = theta_at(j);
            for i in 0..=j {
                let short_node =
                    short_rate + oas + theta_j + vol_step * (2.0 * i as f64 - j as f64);
                let disc = (-short_node * dt).exp();
                next[i] = (values[i] + values[i + 1]) * disc * 0.5;
            }
            // Callability strikes are clean prices: a coupon paid on the
            // exercise date is received in addition to the call/put amount.
            // Apply the exercise condition to the ex-coupon continuation and
            // add any date cashflow afterwards.  (Applying the cap after the
            // coupon would incorrectly cap the coupon as well.)
            for c in exercise_at(j) {
                for value in &mut next {
                    *value = if c.is_call {
                        value.min_val(c.clean_strike)
                    } else {
                        value.max_val(c.clean_strike)
                    };
                }
            }
            for (idx, &ct) in self.coupon_times.iter().enumerate() {
                let step_idx = ((ct / dt).val().round() as usize).min(n);
                if step_idx == j {
                    let coupon = self.cashflow_at(idx)
                        - if idx + 1 == self.coupon_times.len() {
                            self.face_amount
                        } else {
                            T::zero()
                        };
                    next.iter_mut().for_each(|value| *value += coupon);
                }
            }
            values = next;
        }
        values[0]
    }
    pub fn solve_callable_bond_oas(
        &self,
        market: f64,
        short: f64,
        vol: f64,
        steps: usize,
    ) -> (f64, f64, f64) {
        let mut x = 0.;
        for _ in 0..50 {
            let p = self
                .tree_callable_bond_price(
                    T::from_f64(short),
                    T::from_f64(vol),
                    T::from_f64(x),
                    steps,
                )
                .val()
                - market;
            let h = 1e-5;
            let d = (self
                .tree_callable_bond_price(
                    T::from_f64(short),
                    T::from_f64(vol),
                    T::from_f64(x + h),
                    steps,
                )
                .val()
                - self
                    .tree_callable_bond_price(
                        T::from_f64(short),
                        T::from_f64(vol),
                        T::from_f64(x - h),
                        steps,
                    )
                    .val())
                / (2. * h);
            if d.abs() < 1e-12 {
                break;
            }
            x -= p / d;
        }
        let h = 1e-4;
        let pp = self
            .tree_callable_bond_price(
                T::from_f64(short),
                T::from_f64(vol),
                T::from_f64(x + h),
                steps,
            )
            .val();
        let pm = self
            .tree_callable_bond_price(
                T::from_f64(short),
                T::from_f64(vol),
                T::from_f64(x - h),
                steps,
            )
            .val();
        let p0 = self
            .tree_callable_bond_price(T::from_f64(short), T::from_f64(vol), T::from_f64(x), steps)
            .val();
        let dur = -(pp - pm) / (2. * h * p0);
        let conv = (pp - 2. * p0 + pm) / (h * h * p0);
        (x, dur, conv)
    }
}

/// Free-function spelling for callers that keep the bond as a value object.
pub fn black_callable_bond_price<T: Real>(
    bond: &CallableFixedRateBond<T>,
    rate: T,
    yield_vol: T,
) -> T {
    bond.black_callable_bond_price(rate, yield_vol)
}

pub fn tree_callable_bond_price<T: Real>(
    bond: &CallableFixedRateBond<T>,
    short_rate: T,
    rate_vol: T,
    oas: T,
    steps: usize,
) -> T {
    bond.tree_callable_bond_price(short_rate, rate_vol, oas, steps)
}

pub fn solve_callable_bond_oas(
    bond: &CallableFixedRateBond<f64>,
    market_price: f64,
    short_rate: f64,
    rate_vol: f64,
    steps: usize,
) -> (f64, f64, f64) {
    bond.solve_callable_bond_oas(market_price, short_rate, rate_vol, steps)
}

#[derive(Clone, Debug)]
pub struct CallableBondOasSpec<T: Real = f64> {
    pub market_dirty_price: T,
    pub face_amount: T,
    pub coupon_rate: T,
    pub payment_times: Vec<f64>,
    pub call_price: Option<T>,
    pub put_price: Option<T>,
    pub short_rate: T,
    pub short_rate_vol: T,
    pub tree_steps: usize,
    pub smoothing_eps: f64,
}

#[derive(Copy, Clone, Debug)]
pub struct CallableBondOasReport<T: Real> {
    pub oas: T,
    pub straight_bond_price: T,
    pub embedded_option_value: T,
    pub effective_duration: T,
    pub effective_convexity: T,
}

fn smooth_pos<T: Real>(x: T, eps: f64) -> T {
    let e = eps.max(1e-8);
    (x + (x * x + T::from_f64(e * e)).sqrt()) * 0.5
}

pub fn price_callable_tree_generic<T: Real>(s: &CallableBondOasSpec<T>, oas: T) -> (T, T) {
    if s.payment_times.is_empty() {
        return (T::zero(), T::zero());
    }
    let maturity = s
        .payment_times
        .iter()
        .copied()
        .fold(0.0_f64, f64::max)
        .max(1e-6);
    let n = s.tree_steps.max(4);
    let dt = maturity / (n as f64);
    let vol_step = s.short_rate_vol.max_val(T::zero()) * dt.sqrt();
    let eps = s.smoothing_eps.max(1e-8);

    let mut step_coupons = vec![T::zero(); n + 1];
    for (idx, &t) in s.payment_times.iter().enumerate() {
        let prev_t = if idx == 0 {
            0.0
        } else {
            s.payment_times[idx - 1]
        };
        let tau = (t - prev_t).max(0.0);
        let step_idx = ((t / dt).round() as usize).min(n);
        step_coupons[step_idx] += s.face_amount * s.coupon_rate * tau;
    }

    let term_val = s.face_amount + step_coupons[n];
    let mut callable_vals = vec![term_val; n + 1];
    let mut straight_vals = vec![term_val; n + 1];

    for j in (0..n).rev() {
        let mut next_call = vec![T::zero(); j + 1];
        let mut next_str = vec![T::zero(); j + 1];
        for i in 0..=j {
            let state_shift = 2.0 * (i as f64) - (j as f64);
            let r_node = s.short_rate + oas + vol_step * state_shift;
            let df = (-r_node * dt).exp();
            let cont_str = (straight_vals[i] + straight_vals[i + 1]) * df * 0.5;
            let mut cont_call = (callable_vals[i] + callable_vals[i + 1]) * df * 0.5;
            if j > 0 {
                if let Some(k_call) = s.call_price {
                    cont_call -= smooth_pos(cont_call - k_call, eps);
                }
                if let Some(k_put) = s.put_price {
                    cont_call += smooth_pos(k_put - cont_call, eps);
                }
            }
            next_str[i] = cont_str + step_coupons[j];
            next_call[i] = cont_call + step_coupons[j];
        }
        straight_vals = next_str;
        callable_vals = next_call;
    }
    (callable_vals[0], straight_vals[0])
}

pub fn analyze_callable_bond_oas_ad<T: Real>(
    s: &CallableBondOasSpec<T>,
) -> CallableBondOasReport<T> {
    let s_f64 = CallableBondOasSpec::<f64> {
        market_dirty_price: s.market_dirty_price.val(),
        face_amount: s.face_amount.val(),
        coupon_rate: s.coupon_rate.val(),
        payment_times: s.payment_times.clone(),
        call_price: s.call_price.map(|v| v.val()),
        put_price: s.put_price.map(|v| v.val()),
        short_rate: s.short_rate.val(),
        short_rate_vol: s.short_rate_vol.val(),
        tree_steps: s.tree_steps,
        smoothing_eps: s.smoothing_eps,
    };
    let eval_f64 =
        |z: f64| -> f64 { price_callable_tree_generic(&s_f64, z).0 - s_f64.market_dirty_price };
    let mut lo = -0.15_f64;
    let mut hi = 0.50_f64;
    let f_lo = eval_f64(lo);
    for _ in 0..52 {
        let mid = 0.5 * (lo + hi);
        let f_mid = eval_f64(mid);
        if (f_mid >= 0.0) == (f_lo >= 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mut x_root = 0.5 * (lo + hi);
    let h = 1e-5;
    for _ in 0..4 {
        let r0 = eval_f64(x_root);
        let dr = (eval_f64(x_root + h) - eval_f64(x_root - h)) / (2.0 * h);
        if dr.abs() > 1e-12 {
            x_root = (x_root - r0 / dr).clamp(-0.25, 1.0);
        }
    }
    let r_up = eval_f64(x_root + h);
    let r_0 = eval_f64(x_root);
    let r_dn = eval_f64(x_root - h);
    let f_x = (r_up - r_dn) / (2.0 * h);
    let f_xx = (r_up - 2.0 * r_0 + r_dn) / (h * h);

    let residual_at =
        |z_t: T| -> T { price_callable_tree_generic(s, z_t).0 - s.market_dirty_price };
    let oas = ql_ad::ift_lift_1d(x_root, residual_at, T::from_f64(f_x), f_xx);
    let (cp, sp) = price_callable_tree_generic(s, oas);
    let dr = 1e-4;
    let mut s_up = s.clone();
    s_up.short_rate += T::from_f64(dr);
    let mut s_dn = s.clone();
    s_dn.short_rate -= T::from_f64(dr);
    let p_plus = price_callable_tree_generic(&s_up, oas).0;
    let p_minus = price_callable_tree_generic(&s_dn, oas).0;
    let cp_safe = cp.max_val(T::from_f64(1e-12));
    let dur = -(p_plus - p_minus) / (cp_safe * (2.0 * dr));
    let conv = (p_plus + p_minus - cp * 2.0) / (cp_safe * (dr * dr));
    CallableBondOasReport {
        oas,
        straight_bond_price: sp,
        embedded_option_value: sp - cp,
        effective_duration: dur,
        effective_convexity: conv,
    }
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn callable_bond_tree_oas_vega_and_ad() {
        let s_high_call = CallableBondOasSpec {
            market_dirty_price: 100.0,
            face_amount: 100.0,
            coupon_rate: 0.04,
            payment_times: vec![1.0, 2.0, 3.0],
            call_price: Some(1000.0),
            put_price: None,
            short_rate: 0.03,
            short_rate_vol: 0.015,
            tree_steps: 12,
            smoothing_eps: 1e-4,
        };
        let r_high = analyze_callable_bond_oas_ad(&s_high_call);
        assert!(r_high.embedded_option_value.abs() < 1e-5);

        let mut s_call = s_high_call.clone();
        s_call.call_price = Some(100.5);
        let r_call = analyze_callable_bond_oas_ad(&s_call);
        assert!(r_call.embedded_option_value > 0.0);
        assert!(r_call.oas < r_high.oas);

        let h = 1e-5;
        let eval_v = |vol: f64| {
            let mut sp = s_call.clone();
            sp.short_rate_vol = vol;
            analyze_callable_bond_oas_ad(&sp).oas
        };
        let fd = (eval_v(0.015 + h) - eval_v(0.015 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sd = CallableBondOasSpec {
            market_dirty_price: D::from_value(100.0),
            face_amount: D::from_value(100.0),
            coupon_rate: D::from_value(0.04),
            payment_times: vec![1.0, 2.0, 3.0],
            call_price: Some(D::from_value(100.5)),
            put_price: None,
            short_rate: D::from_value(0.03),
            short_rate_vol: D::variable(0.015, 0),
            tree_steps: 12,
            smoothing_eps: 1e-4,
        };
        let ad = analyze_callable_bond_oas_ad(&sd).oas;
        assert!(ad.d[0] < 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-4);
    }
}

#[derive(Clone, Debug)]
pub struct BlackKarasinskiCallableFrnSpec<T: Real = f64> {
    pub face_amount: T,
    pub quoted_margin: T,
    pub cap_rate: T,
    pub floor_rate: T,
    pub maturity: f64,
    pub coupon_steps: usize,
    pub mean_reversion: T,
    pub log_rate_vol: T,
    pub flat_zero_rate: T,
    pub call_price: T,
    pub first_call_time: f64,
    pub target_market_price: T,
    pub smoothing_eps: f64,
}

#[derive(Copy, Clone, Debug)]
pub struct BlackKarasinskiCallableFrnReport<T: Real = f64> {
    pub oas_spread: T,
    pub model_price_at_oas: T,
    pub straight_frn_price: T,
    pub capped_floored_no_call_price: T,
    pub embedded_collar_value: T,
    pub embedded_call_option_value: T,
    pub effective_duration: T,
    pub effective_convexity: T,
}

fn bk_tree_frn_prices<T: Real>(s: &BlackKarasinskiCallableFrnSpec<T>, oas: T) -> (T, T, T) {
    let n = s.coupon_steps.max(2);
    let dt = s.maturity.max(0.25) / n as f64;
    let a = s.mean_reversion.max_val(T::from_f64(1e-6));
    let sigma = s.log_rate_vol.max_val(T::from_f64(1e-6));
    let r0 = s.flat_zero_rate.max_val(T::from_f64(1e-4));
    let dx = sigma * (3.0 * dt).sqrt();
    let j_max = ((0.184 / (a.val() * dt).max(1e-4)).ceil() as isize)
        .max(2)
        .min(n as isize);
    let eps = T::from_f64(s.smoothing_eps.max(1e-4));

    let idx_of = |j: isize, m: isize| -> usize { (j + m) as usize };
    let branch_probs = |j: isize, m: isize| -> (isize, T, T, T) {
        let aj = a * (j as f64 * dt);
        if j >= m && m > 0 {
            let pu = T::from_f64(7.0 / 6.0) + (aj * aj - aj * 3.0) * 0.5;
            let pm = T::from_f64(-1.0 / 3.0) - aj * aj + aj * 2.0;
            let pd = T::from_f64(1.0 / 6.0) + (aj * aj - aj) * 0.5;
            (j - 1, pu, pm, pd)
        } else if j <= -m && m > 0 {
            let pu = T::from_f64(1.0 / 6.0) + (aj * aj + aj) * 0.5;
            let pm = T::from_f64(-1.0 / 3.0) - aj * aj - aj * 2.0;
            let pd = T::from_f64(7.0 / 6.0) + (aj * aj + aj * 3.0) * 0.5;
            (j + 1, pu, pm, pd)
        } else {
            let pu = T::from_f64(1.0 / 6.0) + (aj * aj - aj) * 0.5;
            let pm = T::from_f64(2.0 / 3.0) - aj * aj;
            let pd = T::from_f64(1.0 / 6.0) + (aj * aj + aj) * 0.5;
            (j, pu, pm, pd)
        }
    };

    // Forward Arrow-Debreu calibration of Black-Karasinski median rate r*(t_i) = exp(alpha_i)
    let mut r_star = vec![r0; n];
    let mut q_state = vec![T::one()];
    for i in 0..n {
        let m = (i as isize).min(j_max);
        let p_target = (-r0 * ((i + 1) as f64 * dt)).exp();
        let mut alpha = r0.ln();
        for _ in 0..6 {
            let mut f = -p_target;
            let mut df = T::zero();
            for j in -m..=m {
                let qj = q_state[idx_of(j, m)];
                let r_ij = (alpha + dx * (j as f64)).exp();
                let disc = (-r_ij * dt).exp();
                f += qj * disc;
                df -= qj * disc * r_ij * dt;
            }
            alpha -= f / df.min_val(T::from_f64(-1e-12));
        }
        r_star[i] = alpha.exp();
        if i + 1 < n {
            let m_next = ((i + 1) as isize).min(j_max);
            let mut q_next = vec![T::zero(); (2 * m_next + 1) as usize];
            for j in -m..=m {
                let qj = q_state[idx_of(j, m)];
                let r_ij = r_star[i] * (dx * (j as f64)).exp();
                let disc = (-r_ij * dt).exp();
                let (k_mid, pu, pm, pd) = branch_probs(j, m);
                q_next[idx_of(k_mid + 1, m_next)] += qj * disc * pu;
                q_next[idx_of(k_mid, m_next)] += qj * disc * pm;
                q_next[idx_of(k_mid - 1, m_next)] += qj * disc * pd;
            }
            q_state = q_next;
        }
    }

    let m_last = ((n - 1) as isize).min(j_max);
    let mut v_call = vec![s.face_amount; (2 * m_last + 1) as usize];
    let mut v_cap = vec![s.face_amount; (2 * m_last + 1) as usize];
    let mut v_str = vec![s.face_amount; (2 * m_last + 1) as usize];

    for i in (0..n).rev() {
        let t_i = i as f64 * dt;
        let m = (i as isize).min(j_max);
        let m_next = if i + 1 < n {
            ((i + 1) as isize).min(j_max)
        } else {
            m
        };
        let mut next_call = vec![T::zero(); (2 * m + 1) as usize];
        let mut next_cap = vec![T::zero(); (2 * m + 1) as usize];
        let mut next_str = vec![T::zero(); (2 * m + 1) as usize];
        for j in -m..=m {
            let idx = idx_of(j, m);
            let r_ij = r_star[i] * (dx * (j as f64)).exp();
            let disc_oas = (-(r_ij + oas) * dt).exp();
            let raw_rate = r_ij + s.quoted_margin;
            let floor_diff = raw_rate - s.floor_rate;
            let floored =
                (raw_rate + s.floor_rate + (floor_diff * floor_diff + eps * eps).sqrt()) * 0.5;
            let cap_diff = floored - s.cap_rate;
            let capped_rate =
                (floored + s.cap_rate - (cap_diff * cap_diff + eps * eps).sqrt()) * 0.5;

            let cpn_str = s.face_amount * raw_rate * dt;
            let cpn_cap = s.face_amount * capped_rate * dt;

            let (ev_call, ev_cap, ev_str) = if i + 1 == n {
                (s.face_amount, s.face_amount, s.face_amount)
            } else {
                let (k_mid, pu, pm, pd) = branch_probs(j, m);
                let iu = idx_of(k_mid + 1, m_next);
                let im = idx_of(k_mid, m_next);
                let id = idx_of(k_mid - 1, m_next);
                (
                    pu * v_call[iu] + pm * v_call[im] + pd * v_call[id],
                    pu * v_cap[iu] + pm * v_cap[im] + pd * v_cap[id],
                    pu * v_str[iu] + pm * v_str[im] + pd * v_str[id],
                )
            };
            let cont_str = (ev_str + cpn_str) * disc_oas;
            let cont_cap = (ev_cap + cpn_cap) * disc_oas;
            let cont_ex_cpn = ev_call * disc_oas;
            let called_ex_cpn = if i > 0 && t_i + 1e-9 >= s.first_call_time {
                let d_call = cont_ex_cpn - s.call_price;
                (cont_ex_cpn + s.call_price - (d_call * d_call + eps * eps).sqrt()) * 0.5
            } else {
                cont_ex_cpn
            };
            next_str[idx] = cont_str;
            next_cap[idx] = cont_cap;
            next_call[idx] = called_ex_cpn + cpn_cap * disc_oas;
        }
        v_str = next_str;
        v_cap = next_cap;
        v_call = next_call;
    }
    (v_call[0], v_cap[0], v_str[0])
}

/// Black-Karasinski (1991) lognormal short-rate trinomial tree OAS, collar decomposition,
/// and effective duration/convexity engine with single-graph `ift_lift_1d` AD.
pub fn analyze_black_karasinski_callable_frn_oas_ad<T: Real>(
    s: &BlackKarasinskiCallableFrnSpec<T>,
) -> BlackKarasinskiCallableFrnReport<T> {
    let s_f64 = BlackKarasinskiCallableFrnSpec::<f64> {
        face_amount: s.face_amount.val(),
        quoted_margin: s.quoted_margin.val(),
        cap_rate: s.cap_rate.val(),
        floor_rate: s.floor_rate.val(),
        maturity: s.maturity,
        coupon_steps: s.coupon_steps,
        mean_reversion: s.mean_reversion.val(),
        log_rate_vol: s.log_rate_vol.val(),
        flat_zero_rate: s.flat_zero_rate.val(),
        call_price: s.call_price.val(),
        first_call_time: s.first_call_time,
        target_market_price: s.target_market_price.val(),
        smoothing_eps: s.smoothing_eps,
    };
    let eval_f64 =
        |oas: f64| -> f64 { bk_tree_frn_prices(&s_f64, oas).0 - s_f64.target_market_price };
    let mut lo = -0.08_f64;
    let mut hi = 0.25_f64;
    let f_lo = eval_f64(lo);
    for _ in 0..32 {
        let mid = 0.5 * (lo + hi);
        if (eval_f64(mid) >= 0.0) == (f_lo >= 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mut oas_root = 0.5 * (lo + hi);
    let h = 1e-5;
    for _ in 0..5 {
        let r0 = eval_f64(oas_root);
        let dr = (eval_f64(oas_root + h) - eval_f64(oas_root - h)) / (2.0 * h);
        if dr.abs() > 1e-12 {
            oas_root = (oas_root - r0 / dr).clamp(-0.08, 0.25);
        }
    }
    let r_up = eval_f64(oas_root + h);
    let r_0 = eval_f64(oas_root);
    let r_dn = eval_f64(oas_root - h);
    let f_xx = (r_up - 2.0 * r_0 + r_dn) / (h * h);
    let oas_t0 = T::from_f64(oas_root);
    let h_t = T::from_f64(h);
    let f_x =
        (bk_tree_frn_prices(s, oas_t0 + h_t).0 - bk_tree_frn_prices(s, oas_t0 - h_t).0) / (2.0 * h);
    let oas_spread = ql_ad::ift_lift_1d(
        oas_root,
        |z| bk_tree_frn_prices(s, z).0 - s.target_market_price,
        f_x,
        f_xx,
    );
    let (model_price_at_oas, capped_floored_no_call_price, straight_frn_price) =
        bk_tree_frn_prices(s, oas_spread);
    let embedded_collar_value = capped_floored_no_call_price - straight_frn_price;
    let embedded_call_option_value = capped_floored_no_call_price - model_price_at_oas;

    let dr = 1e-4;
    let mut s_up = s.clone();
    s_up.flat_zero_rate += T::from_f64(dr);
    let mut s_dn = s.clone();
    s_dn.flat_zero_rate -= T::from_f64(dr);
    let p_up = bk_tree_frn_prices(&s_up, oas_spread).0;
    let p_dn = bk_tree_frn_prices(&s_dn, oas_spread).0;
    let p_safe = model_price_at_oas.max_val(T::from_f64(1e-12));
    let effective_duration = -(p_up - p_dn) / (p_safe * (2.0 * dr));
    let effective_convexity = (p_up + p_dn - model_price_at_oas * 2.0) / (p_safe * (dr * dr));

    BlackKarasinskiCallableFrnReport {
        oas_spread,
        model_price_at_oas,
        straight_frn_price,
        capped_floored_no_call_price,
        embedded_collar_value,
        embedded_call_option_value,
        effective_duration,
        effective_convexity,
    }
}

#[cfg(test)]
mod phase20_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn black_karasinski_callable_frn_call_option_oas_and_ad() {
        let base = BlackKarasinskiCallableFrnSpec {
            face_amount: 100.0,
            quoted_margin: 0.008,
            cap_rate: 0.055,
            floor_rate: 0.020,
            maturity: 3.0,
            coupon_steps: 9,
            mean_reversion: 0.15,
            log_rate_vol: 0.22,
            flat_zero_rate: 0.035,
            call_price: 100.2,
            first_call_time: 1.0,
            target_market_price: 99.8,
            smoothing_eps: 1e-3,
        };
        let rep = analyze_black_karasinski_callable_frn_oas_ad(&base);
        assert!((rep.model_price_at_oas - 99.8).abs() < 1e-6);
        assert!(rep.embedded_call_option_value > 0.0);

        let h = 1e-5;
        let eval_v = |v: f64| {
            let mut s = base.clone();
            s.log_rate_vol = v;
            analyze_black_karasinski_callable_frn_oas_ad(&s).oas_spread
        };
        let fd = (eval_v(0.22 + h) - eval_v(0.22 - h)) / (2.0 * h);
        type D = Dual<1>;
        let spec_d = BlackKarasinskiCallableFrnSpec {
            face_amount: D::from_value(100.0),
            quoted_margin: D::from_value(0.008),
            cap_rate: D::from_value(0.055),
            floor_rate: D::from_value(0.020),
            maturity: 3.0,
            coupon_steps: 9,
            mean_reversion: D::from_value(0.15),
            log_rate_vol: D::variable(0.22, 0),
            flat_zero_rate: D::from_value(0.035),
            call_price: D::from_value(100.2),
            first_call_time: 1.0,
            target_market_price: D::from_value(99.8),
            smoothing_eps: 1e-3,
        };
        let ad = analyze_black_karasinski_callable_frn_oas_ad(&spec_d).oas_spread;
        assert!((ad.d[0] - fd).abs() < 1e-4);
    }
}
